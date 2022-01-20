/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_implics.c
 * @brief  implics presolver
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_implics.h"


#define PRESOL_NAME            "implics"
#define PRESOL_DESC            "implication graph aggregator"
#define PRESOL_PRIORITY          -10000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_MEDIUM /* timing of the presolver (fast, medium, or exhaustive) */


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyImplics)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolImplics(scip) );

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecImplics)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR** bdchgvars;
   SCIP_BOUNDTYPE* bdchgtypes;
   SCIP_Real* bdchgvals;
   SCIP_VAR** aggrvars;
   SCIP_VAR** aggraggvars;
   SCIP_Real* aggrcoefs;
   SCIP_Real* aggrconsts;
   int nbdchgs;
   int naggregations;
   int nbinvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* initialize fixing and aggregation storages */
   bdchgvars = NULL;
   bdchgtypes = NULL;
   bdchgvals = NULL;
   nbdchgs = 0;
   aggrvars = NULL;
   aggraggvars = NULL;
   aggrcoefs = NULL;
   aggrconsts = NULL;
   naggregations = 0;

   /* get active binary problem variables */
   vars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);

   /* look for variable implications in x == 0 and x == 1 with the same implied variable:
    *  x = 0 -> y = lb, and x = 1 -> y = lb: fix y to lb
    *  x = 0 -> y = lb, and x = 1 -> y = ub: aggregate y == lb + (ub-lb)x
    *  x = 0 -> y = ub, and x = 1 -> y = lb: aggregate y == ub - (ub-lb)x
    *  x = 0 -> y = ub, and x = 1 -> y = ub: fix y to ub
    * the fixings and aggregations are stored in a buffer and applied afterwards, because fixing and aggregation
    * would modify the vars array and the implication arrays
    */
   for( v = 0; v < nbinvars; ++v )
   {
      SCIP_VAR** implvars[2];
      SCIP_BOUNDTYPE* impltypes[2];
      SCIP_Real* implbounds[2];
      int nimpls[2];
      int varfixing;
      int i0;
      int i1;

      /* don't perform presolving operations on deleted variables */
      if( SCIPvarIsDeleted(vars[v]) )
         continue;

      /* get implications for given variable */
      for( varfixing = 0; varfixing < 2; ++varfixing )
      {
         implvars[varfixing] = SCIPvarGetImplVars(vars[v], (SCIP_Bool)varfixing);
         impltypes[varfixing] = SCIPvarGetImplTypes(vars[v], (SCIP_Bool)varfixing);
         implbounds[varfixing] = SCIPvarGetImplBounds(vars[v], (SCIP_Bool)varfixing);
         nimpls[varfixing] = SCIPvarGetNImpls(vars[v], (SCIP_Bool)varfixing);
      }

      /* scan implication arrays for equal variables */
      i0 = 0;
      i1 = 0;
      while( i0 < nimpls[0] && i1 < nimpls[1] )
      {
         int index0;
         int index1;

         /* scan the binary or non-binary part of the implication arrays */
         index0 = SCIPvarGetIndex(implvars[0][i0]);
         index1 = SCIPvarGetIndex(implvars[1][i1]);
         while( index0 < index1 )
         {
            i0++;
            if( i0 == nimpls[0] )
            {
               index0 = -1;
               break;
            }
            index0 = SCIPvarGetIndex(implvars[0][i0]);
         }
         while( index1 < index0 )
         {
            i1++;
            if( i1 == nimpls[1] )
            {
               index1 = -1;
               break;
            }
            index1 = SCIPvarGetIndex(implvars[1][i1]);
         }
         /**@todo for all implied binary variables y, check the cliques of x == !varfixing if y is contained */

         if( index0 == index1 )
         {
            assert(index0 >= 0);
            assert(i0 < nimpls[0]);
            assert(i1 < nimpls[1]);
            assert(implvars[0][i0] == implvars[1][i1]);

            if( impltypes[0][i0] == impltypes[1][i1] )
            {
               /* found implication x = 0 -> y >= b / y <= b  and  x = 1 -> y >= c / y <= c
                *   =>  change bound y >= min(b,c) / y <= max(b,c)
                */
               SCIP_CALL( SCIPreallocBufferArray(scip, &bdchgvars, nbdchgs+1) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &bdchgtypes, nbdchgs+1) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &bdchgvals, nbdchgs+1) );
               bdchgvars[nbdchgs] = implvars[0][i0];
               bdchgtypes[nbdchgs] = impltypes[0][i0];
               if( impltypes[0][i0] == SCIP_BOUNDTYPE_LOWER )
                  bdchgvals[nbdchgs] = MIN(implbounds[0][i0], implbounds[1][i1]);
               else
                  bdchgvals[nbdchgs] = MAX(implbounds[0][i0], implbounds[1][i1]);

               SCIPdebugMsg(scip, " -> <%s> = 0 -> <%s> %s %g, and <%s> = 1 -> <%s> %s %g:  tighten <%s> %s %g\n",
                  SCIPvarGetName(vars[v]), SCIPvarGetName(implvars[0][i0]),
                  impltypes[0][i0] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", implbounds[0][i0],
                  SCIPvarGetName(vars[v]), SCIPvarGetName(implvars[1][i1]),
                  impltypes[1][i1] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", implbounds[1][i1],
                  SCIPvarGetName(bdchgvars[nbdchgs]), bdchgtypes[nbdchgs] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
                  bdchgvals[nbdchgs]);

               nbdchgs++;
            }
            else
            {
               SCIP_Real implvarlb;
               SCIP_Real implvarub;

               implvarlb = SCIPvarGetLbGlobal(implvars[0][i0]);
               implvarub = SCIPvarGetUbGlobal(implvars[0][i0]);

               if( impltypes[0][i0] == SCIP_BOUNDTYPE_UPPER
                  && SCIPisEQ(scip, implbounds[0][i0], implvarlb)
                  && SCIPisEQ(scip, implbounds[1][i1], implvarub) )
               {
                  /* found implication x = 0 -> y = lb and x = 1 -> y = ub  =>  aggregate y = lb + (ub-lb) * x */
                  SCIP_CALL( SCIPreallocBufferArray(scip, &aggrvars, naggregations+1) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &aggraggvars, naggregations+1) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &aggrcoefs, naggregations+1) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &aggrconsts, naggregations+1) );
                  aggrvars[naggregations] = implvars[0][i0];
                  aggraggvars[naggregations] = vars[v];
                  aggrcoefs[naggregations] = implvarub - implvarlb;
                  aggrconsts[naggregations] = implvarlb;

                  SCIPdebugMsg(scip, " -> <%s> = 0 -> <%s> = %g, and <%s> = 1 -> <%s> = %g:  aggregate <%s> = %g %+g<%s>\n",
                     SCIPvarGetName(vars[v]), SCIPvarGetName(implvars[0][i0]), implbounds[0][i0],
                     SCIPvarGetName(vars[v]), SCIPvarGetName(implvars[1][i1]), implbounds[1][i1],
                     SCIPvarGetName(aggrvars[naggregations]), aggrconsts[naggregations], aggrcoefs[naggregations],
                     SCIPvarGetName(aggraggvars[naggregations]));

                  naggregations++;
               }
               else if( impltypes[0][i0] == SCIP_BOUNDTYPE_LOWER
                  && SCIPisEQ(scip, implbounds[0][i0], implvarub)
                  && SCIPisEQ(scip, implbounds[1][i1], implvarlb) )
               {
                  /* found implication x = 0 -> y = ub and x = 1 -> y = lb  =>  aggregate y = ub - (ub-lb) * x */
                  SCIP_CALL( SCIPreallocBufferArray(scip, &aggrvars, naggregations+1) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &aggraggvars, naggregations+1) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &aggrcoefs, naggregations+1) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &aggrconsts, naggregations+1) );
                  aggrvars[naggregations] = implvars[0][i0];
                  aggraggvars[naggregations] = vars[v];
                  aggrcoefs[naggregations] = implvarlb - implvarub;
                  aggrconsts[naggregations] = implvarub;

                  SCIPdebugMsg(scip, " -> <%s> = 0 -> <%s> = %g, and <%s> = 1 -> <%s> = %g:  aggregate <%s> = %g %+g<%s>\n",
                     SCIPvarGetName(vars[v]), SCIPvarGetName(implvars[0][i0]), implbounds[0][i0],
                     SCIPvarGetName(vars[v]), SCIPvarGetName(implvars[1][i1]), implbounds[1][i1],
                     SCIPvarGetName(aggrvars[naggregations]), aggrconsts[naggregations], aggrcoefs[naggregations],
                     SCIPvarGetName(aggraggvars[naggregations]));

                  naggregations++;
               }
            }

            /* process the next implications */
            i0++;
            i1++;
         }
      }
   }

   /**@todo check cliques of x == 0 and x == 1 for equal entries y == b -> fix y == !b */

   /* perform the bound changes
    *
    * note, that we cannot assume y to be active (see var.c: varRemoveImplicsVbs()), but it should not cause any 
    * troubles as this case seems to be handled correctly in SCIPtightenVarLb/Ub().
    */
   for( v = 0; v < nbdchgs && *result != SCIP_CUTOFF; ++v )
   {
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      assert(bdchgtypes != NULL);
      assert(bdchgvars != NULL);
      assert(bdchgvals != NULL);

      if( bdchgtypes[v] == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPtightenVarLb(scip, bdchgvars[v], bdchgvals[v], FALSE, &infeasible, &tightened) );
      }
      else
      {
         SCIP_CALL( SCIPtightenVarUb(scip, bdchgvars[v], bdchgvals[v], FALSE, &infeasible, &tightened) );
      }

      if( infeasible )
      {
         SCIPdebugMsg(scip, " -> infeasible bound change <%s> %s %g\n", SCIPvarGetName(bdchgvars[v]),
            bdchgtypes[v] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", bdchgvals[v]);
         *result = SCIP_CUTOFF;
      }
      else if( tightened )
      {
         (*nchgbds)++;
         *result = SCIP_SUCCESS;
      }
   }

   /* perform the aggregations
    * 
    * note, that we cannot assume y to be active (see var.c: varRemoveImplicsVbs()), but it should not cause any 
    * troubles as this case seems to be handled correctly in SCIPaggregateVars().
    */
   for( v = 0; v < naggregations && *result != SCIP_CUTOFF; ++v )
   {
      SCIP_Bool infeasible;
      SCIP_Bool redundant;
      SCIP_Bool aggregated;

      assert(aggrvars != NULL);
      assert(aggraggvars != NULL);
      assert(aggrcoefs != NULL);
      assert(aggrconsts != NULL);

      /* aggregation y = const + coef * x  =>  y - coef * x = const */
      SCIP_CALL( SCIPaggregateVars(scip, aggrvars[v], aggraggvars[v], 1.0, -aggrcoefs[v], aggrconsts[v],
            &infeasible, &redundant, &aggregated) );
      if( infeasible )
      {
         SCIPdebugMsg(scip, " -> infeasible aggregation <%s> = %g %+g<%s>\n",
            SCIPvarGetName(aggrvars[v]), aggrconsts[v], aggrcoefs[v], SCIPvarGetName(aggraggvars[v]));
         *result = SCIP_CUTOFF;
      }
      else if( aggregated )
      {
         (*naggrvars)++;
         *result = SCIP_SUCCESS;
      }
   }

   /* free the storage buffers */
   SCIPfreeBufferArrayNull(scip, &aggrconsts);
   SCIPfreeBufferArrayNull(scip, &aggrcoefs);
   SCIPfreeBufferArrayNull(scip, &aggraggvars);
   SCIPfreeBufferArrayNull(scip, &aggrvars);
   SCIPfreeBufferArrayNull(scip, &bdchgvals);
   SCIPfreeBufferArrayNull(scip, &bdchgtypes);
   SCIPfreeBufferArrayNull(scip, &bdchgvars);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the implics presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolImplics(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presolptr;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presolptr, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING, presolExecImplics, NULL) );

   assert(presolptr != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presolptr, presolCopyImplics) );

   return SCIP_OKAY;
}
