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

/**@file   presol_inttobinary.c
 * @brief  presolver that converts integer variables with domain [a,a+1] to binaries
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_inttobinary.h"


#define PRESOL_NAME            "inttobinary"
#define PRESOL_DESC            "converts integer variables with domain [a,a+1] to binaries"
#define PRESOL_PRIORITY        +7000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolver (fast, medium, or exhaustive) */

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyInttobinary)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolInttobinary(scip) );

   return SCIP_OKAY;
}


/** presolving execution method */
static
SCIP_DECL_PRESOLEXEC(presolExecInttobinary)
{  /*lint --e{715}*/
   SCIP_VAR** scipvars;
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPdoNotAggr(scip) )
      return SCIP_OKAY;

   /* get the problem variables */
   scipvars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   nintvars = SCIPgetNIntVars(scip);
   if( nintvars == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* copy the integer variables into an own array, since adding binary variables affects the left-most slots in the
    * array and thereby interferes with our search loop
    */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, &scipvars[nbinvars], nintvars) );

   /* scan the integer variables for possible conversion into binaries;
    * we have to collect the variables first in an own 
    */
   for( v = 0; v < nintvars; ++v )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER);

      /* get variable's bounds */
      lb = SCIPvarGetLbGlobal(vars[v]);
      ub = SCIPvarGetUbGlobal(vars[v]);

      /* check if bounds are exactly one apart */
      if( SCIPisEQ(scip, lb, ub - 1.0) )
      {
         SCIP_VAR* binvar;
         char binvarname[SCIP_MAXSTRLEN];
         SCIP_Bool infeasible;
         SCIP_Bool redundant;
         SCIP_Bool aggregated;

         SCIPdebugMsg(scip, "converting <%s>[%g,%g] into binary variable\n", SCIPvarGetName(vars[v]), lb, ub);

         /* create binary variable */
         (void) SCIPsnprintf(binvarname, SCIP_MAXSTRLEN, "%s_bin", SCIPvarGetName(vars[v]));
         SCIP_CALL( SCIPcreateVar(scip, &binvar, binvarname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
               SCIPvarIsInitial(vars[v]), SCIPvarIsRemovable(vars[v]), NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, binvar) );

         /* aggregate integer and binary variable */
         SCIP_CALL( SCIPaggregateVars(scip, vars[v], binvar, 1.0, -1.0, lb, &infeasible, &redundant, &aggregated) );

         /* release binary variable */
         SCIP_CALL( SCIPreleaseVar(scip, &binvar) );

         /* it can be the case that this aggregation detects an infeasibility; for example, during the copy of the
          * variable bounds from the integer variable to the binary variable, infeasibility can be detected; this can
          * happen because an upper bound or a lower bound of such a variable bound variable was "just" changed and the
          * varbound constraint handler, who would detect that infeasibility (since it was creating it from a varbound
          * constraint), was called before that bound change was detected due to the presolving priorities;
          */
         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            break;
         }

         assert(redundant);
         assert(aggregated);
         (*nchgvartypes)++;
         ++(*naggrvars);
         *result = SCIP_SUCCESS;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the inttobinary presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolInttobinary(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presolptr;

   /* create inttobinary presolver data */
   presoldata = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presolptr, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecInttobinary,
         presoldata) );

   assert(presolptr != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presolptr, presolCopyInttobinary) );

   return SCIP_OKAY;
}
