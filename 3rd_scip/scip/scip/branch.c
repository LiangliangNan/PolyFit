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

/**@file   branch.c
 * @brief  methods for branching rules and branching candidate storage
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Stefan Heinz
 * @author Michael Winkler
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/sepastore.h"
#include "scip/scip.h"
#include "scip/branch.h"
#include "scip/solve.h"

#include "scip/struct_branch.h"

/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that lpcands array can store at least num entries */
static
SCIP_RETCODE ensureLpcandsSize(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(branchcand->nlpcands <= branchcand->lpcandssize);

   if( num > branchcand->lpcandssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->lpcands, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->lpcandssol, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->lpcandsfrac, newsize) );
      branchcand->lpcandssize = newsize;
   }
   assert(num <= branchcand->lpcandssize);

   return SCIP_OKAY;
}

/** ensures, that pseudocands array can store at least num entries */
static
SCIP_RETCODE ensurePseudocandsSize(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(branchcand->npseudocands <= branchcand->pseudocandssize);

   if( num > branchcand->pseudocandssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->pseudocands, newsize) );
      branchcand->pseudocandssize = newsize;
   }
   assert(num <= branchcand->pseudocandssize);

   return SCIP_OKAY;
}

/** ensures, that externcands array can store at least num entries */
static
SCIP_RETCODE ensureExterncandsSize(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(branchcand->nexterncands <= branchcand->externcandssize);

   if( num > branchcand->externcandssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->externcands, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->externcandsscore, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&branchcand->externcandssol, newsize) );
      branchcand->externcandssize = newsize;
   }
   assert(num <= branchcand->externcandssize);

   return SCIP_OKAY;
}



/*
 * branching candidate storage methods
 */

/** creates a branching candidate storage */
SCIP_RETCODE SCIPbranchcandCreate(
   SCIP_BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   )
{
   assert(branchcand != NULL);

   SCIP_ALLOC( BMSallocMemory(branchcand) );
   (*branchcand)->lpcands = NULL;
   (*branchcand)->lpcandssol = NULL;
   (*branchcand)->lpcandsfrac = NULL;
   (*branchcand)->externcands = NULL;
   (*branchcand)->externcandssol = NULL;
   (*branchcand)->externcandsscore = NULL;
   (*branchcand)->pseudocands = NULL;
   (*branchcand)->lpcandssize = 0;
   (*branchcand)->nlpcands = 0;
   (*branchcand)->nimpllpfracs = 0;
   (*branchcand)->npriolpcands = 0;
   (*branchcand)->npriolpbins = 0;
   (*branchcand)->lpmaxpriority = INT_MIN;
   (*branchcand)->externcandssize = 0;
   (*branchcand)->nexterncands = 0;
   (*branchcand)->nprioexterncands = 0;
   (*branchcand)->nprioexternbins = 0;
   (*branchcand)->nprioexternints = 0;
   (*branchcand)->nprioexternimpls = 0;
   (*branchcand)->externmaxpriority = INT_MIN;
   (*branchcand)->pseudocandssize = 0;
   (*branchcand)->npseudocands = 0;
   (*branchcand)->npriopseudocands = 0;
   (*branchcand)->npriopseudobins = 0;
   (*branchcand)->npriopseudoints = 0;
   (*branchcand)->pseudomaxpriority = INT_MIN;

   SCIPbranchcandInvalidate(*branchcand);

   return SCIP_OKAY;
}

/** frees branching candidate storage */
SCIP_RETCODE SCIPbranchcandFree(
   SCIP_BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   )
{
   assert(branchcand != NULL);

   BMSfreeMemoryArrayNull(&(*branchcand)->lpcands);
   BMSfreeMemoryArrayNull(&(*branchcand)->lpcandssol);
   BMSfreeMemoryArrayNull(&(*branchcand)->lpcandsfrac);
   BMSfreeMemoryArrayNull(&(*branchcand)->pseudocands);
   BMSfreeMemoryArrayNull(&(*branchcand)->externcands);
   BMSfreeMemoryArrayNull(&(*branchcand)->externcandsscore);
   BMSfreeMemoryArrayNull(&(*branchcand)->externcandssol);
   BMSfreeMemory(branchcand);

   return SCIP_OKAY;
}

/** resets branching candidates storage */
void SCIPbranchcandInvalidate(
   SCIP_BRANCHCAND*      branchcand          /**< pointer to store branching candidate storage */
   )
{
   assert(branchcand != NULL);

   branchcand->validlpcandslp = -1;
}

/** calculates branching candidates for LP solution branching (fractional variables) */
static
SCIP_RETCODE branchcandCalcLPCands(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(branchcand != NULL);
   assert(stat != NULL);
   assert(branchcand->validlpcandslp <= stat->lpcount);
   assert(lp != NULL);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

   SCIPsetDebugMsg(set, "calculating LP branching candidates: validlp=%" SCIP_LONGINT_FORMAT ", lpcount=%" SCIP_LONGINT_FORMAT "\n",
      branchcand->validlpcandslp, stat->lpcount);

   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      branchcand->lpmaxpriority = INT_MIN / 2;
      branchcand->nlpcands = 0;
      branchcand->npriolpcands = 0;
      branchcand->npriolpbins = 0;
      branchcand->nimpllpfracs = 0;
      branchcand->validlpcandslp = stat->lpcount;

      SCIPsetDebugMsg(set, " LP is unbounded -> no branching candidates\n");
      return SCIP_OKAY;
   }

   /* check, if the current LP branching candidate array is invalid */
   if( branchcand->validlpcandslp < stat->lpcount )
   {
      SCIP_COL** cols;
      SCIP_VAR* var;
      SCIP_COL* col;
      SCIP_Real primsol;
      SCIP_Real frac;
      SCIP_VARTYPE vartype;
      int branchpriority;
      int ncols;
      int c;
      int insertpos;

      SCIPsetDebugMsg(set, " -> recalculating LP branching candidates\n");

      cols = SCIPlpGetCols(lp);
      ncols = SCIPlpGetNCols(lp);

      /* construct the LP branching candidate set, moving the candidates with maximal priority to the front */
      SCIP_CALL( ensureLpcandsSize(branchcand, set, ncols) );

      branchcand->lpmaxpriority = INT_MIN / 2;
      branchcand->nlpcands = 0;
      branchcand->nimpllpfracs = 0;
      branchcand->npriolpcands = 0;
      branchcand->npriolpbins = 0;
      for( c = 0; c < ncols; ++c )
      {
         col = cols[c];
         assert(col != NULL);
         assert(col->lppos == c);
         assert(col->lpipos >= 0);

         primsol = SCIPcolGetPrimsol(col);
         assert(primsol < SCIP_INVALID);
         assert(SCIPsetIsInfinity(set, -col->lb) || SCIPsetIsFeasGE(set, primsol, col->lb));
         assert(SCIPsetIsInfinity(set, col->ub) || SCIPsetIsFeasLE(set, primsol, col->ub));

         var = col->var;
         assert(var != NULL);
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(var) == col);

         /* LP branching candidates are fractional binary and integer variables; implicit variables are kept at the end
          * of the candidates array for some rounding heuristics
          */
         vartype = SCIPvarGetType(var);
         if( vartype == SCIP_VARTYPE_CONTINUOUS )
            continue;

         /* ignore fixed variables (due to numerics, it is possible, that the LP solution of a fixed integer variable
          * (with large fixed value) is fractional in terms of absolute feasibility measure)
          */
         if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
            continue;

         /* check, if the LP solution value is fractional */
         frac = SCIPsetFeasFrac(set, primsol);

         /* The fractionality should not be smaller than -feastol, however, if the primsol is large enough
          * and close to an integer, fixed precision floating point arithmetic might give us values slightly
          * smaller than -feastol. Originally, the "frac >= -feastol"-check was within SCIPsetIsFeasFracIntegral(),
          * however, we relaxed it to "frac >= -2*feastol" and have the stricter check here for small-enough primsols.
          */
         assert(SCIPsetIsGE(set, frac, -SCIPsetFeastol(set)) || (primsol > 1e14 * SCIPsetFeastol(set)));

         if( SCIPsetIsFeasFracIntegral(set, frac) )
            continue;

         /* insert candidate in candidate list */
         branchpriority = SCIPvarGetBranchPriority(var);
         insertpos = branchcand->nlpcands + branchcand->nimpllpfracs;
         assert(insertpos < branchcand->lpcandssize);

         if( vartype == SCIP_VARTYPE_IMPLINT )
            branchpriority = INT_MIN;

         assert(vartype == SCIP_VARTYPE_IMPLINT || branchpriority >= INT_MIN/2);
         /* ensure that implicit variables are stored at the end of the array */
         if( vartype != SCIP_VARTYPE_IMPLINT && branchcand->nimpllpfracs > 0 )
         {
            assert(branchcand->lpcands[branchcand->nlpcands] != NULL
                  && SCIPvarGetType(branchcand->lpcands[branchcand->nlpcands]) == SCIP_VARTYPE_IMPLINT );

            branchcand->lpcands[insertpos] = branchcand->lpcands[branchcand->nlpcands];
            branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[branchcand->nlpcands];
            branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[branchcand->nlpcands];

            insertpos = branchcand->nlpcands;
         }

         if( branchpriority > branchcand->lpmaxpriority )
         {
            /* candidate has higher priority than the current maximum:
             * move it to the front and declare it to be the single best candidate
             */
            if( insertpos != 0 )
            {
               branchcand->lpcands[insertpos] = branchcand->lpcands[0];
               branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[0];
               branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[0];
               insertpos = 0;
            }
            branchcand->npriolpcands = 1;
            branchcand->npriolpbins = (vartype == SCIP_VARTYPE_BINARY ? 1 : 0);
            branchcand->lpmaxpriority = branchpriority;
         }
         else if( branchpriority == branchcand->lpmaxpriority )
         {
            /* candidate has equal priority as the current maximum:
             * move away the first non-maximal priority candidate, move the current candidate to the correct
             * slot (binaries first) and increase the number of maximal priority candidates
             */
            if( insertpos != branchcand->npriolpcands )
            {
               branchcand->lpcands[insertpos] = branchcand->lpcands[branchcand->npriolpcands];
               branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[branchcand->npriolpcands];
               branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[branchcand->npriolpcands];
               insertpos = branchcand->npriolpcands;
            }
            branchcand->npriolpcands++;
            if( vartype == SCIP_VARTYPE_BINARY )
            {
               if( insertpos != branchcand->npriolpbins )
               {
                  branchcand->lpcands[insertpos] = branchcand->lpcands[branchcand->npriolpbins];
                  branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[branchcand->npriolpbins];
                  branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[branchcand->npriolpbins];
                  insertpos = branchcand->npriolpbins;
               }
               branchcand->npriolpbins++;
            }
         }
         /* insert variable at the correct position of the candidates storage */
         branchcand->lpcands[insertpos] = var;
         branchcand->lpcandssol[insertpos] = primsol;
         branchcand->lpcandsfrac[insertpos] = frac;

         /* increase the counter depending on the variable type */
         if( vartype != SCIP_VARTYPE_IMPLINT )
            branchcand->nlpcands++;
         else
            branchcand->nimpllpfracs++;

         SCIPsetDebugMsg(set, " -> candidate %d: var=<%s>, sol=%g, frac=%g, prio=%d (max: %d) -> pos %d\n",
            branchcand->nlpcands, SCIPvarGetName(var), primsol, frac, branchpriority, branchcand->lpmaxpriority,
            insertpos);
      }

#ifndef NDEBUG
      /* in debug mode we assert that the variables are positioned correctly (binaries and integers first,
       * implicit integers last)
       */
      for( c = 0; c < branchcand->nlpcands + branchcand->nimpllpfracs; ++c )
      {
         assert(c >= branchcand->nlpcands || SCIPvarGetType(branchcand->lpcands[c]) != SCIP_VARTYPE_IMPLINT);
         assert(c < branchcand->nlpcands || SCIPvarGetType(branchcand->lpcands[c]) == SCIP_VARTYPE_IMPLINT);
      }
#endif

      branchcand->validlpcandslp = stat->lpcount;
   }
   assert(0 <= branchcand->npriolpcands && branchcand->npriolpcands <= branchcand->nlpcands);

   SCIPsetDebugMsg(set, " -> %d fractional variables (%d of maximal priority)\n", branchcand->nlpcands, branchcand->npriolpcands);

   return SCIP_OKAY;
}

/** gets branching candidates for LP solution branching (fractional variables) */
SCIP_RETCODE SCIPbranchcandGetLPCands(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   SCIP_Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   SCIP_Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*                  nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*                  npriolpcands,       /**< pointer to store the number of candidates with maximal priority, or NULL */
   int*                  nfracimplvars       /**< pointer to store the number of implicit fractional variables, or NULL */
   )
{
   /* calculate branching candidates */
   SCIP_CALL( branchcandCalcLPCands(branchcand, set, stat, lp) );

   /* assign return values */
   if( lpcands != NULL )
      *lpcands = branchcand->lpcands;
   if( lpcandssol != NULL )
      *lpcandssol = branchcand->lpcandssol;
   if( lpcandsfrac != NULL )
      *lpcandsfrac = branchcand->lpcandsfrac;
   if( nlpcands != NULL )
      *nlpcands = branchcand->nlpcands;
   if( npriolpcands != NULL )
      *npriolpcands = (set->branch_preferbinary && branchcand->npriolpbins > 0 ? branchcand->npriolpbins
         : branchcand->npriolpcands);
   if( nfracimplvars != NULL )
      *nfracimplvars = branchcand->nimpllpfracs;

   return SCIP_OKAY;
}

/** gets external branching candidates */
SCIP_RETCODE SCIPbranchcandGetExternCands(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR***           externcands,        /**< pointer to store the array of external branching candidates, or NULL */
   SCIP_Real**           externcandssol,     /**< pointer to store the array of external candidate solution values, or NULL */
   SCIP_Real**           externcandsscore,   /**< pointer to store the array of external candidate scores, or NULL */
   int*                  nexterncands,       /**< pointer to store the number of external branching candidates, or NULL */
   int*                  nprioexterncands,   /**< pointer to store the number of candidates with maximal priority, or NULL */
   int*                  nprioexternbins,    /**< pointer to store the number of binary candidates with maximal priority, or NULL */
   int*                  nprioexternints,    /**< pointer to store the number of integer candidates with maximal priority, or NULL */
   int*                  nprioexternimpls    /**< pointer to store the number of implicit integer candidates with maximal priority, 
                                              *   or NULL */
   )
{
   assert(branchcand != NULL);

   /* assign return values */
   if( externcands != NULL )
      *externcands = branchcand->externcands;
   if( externcandssol != NULL )
      *externcandssol = branchcand->externcandssol;
   if( externcandsscore != NULL )
      *externcandsscore = branchcand->externcandsscore;
   if( nexterncands != NULL )
      *nexterncands = branchcand->nexterncands;
   if( nprioexterncands != NULL )
      *nprioexterncands = branchcand->nprioexterncands;
   if( nprioexternbins != NULL )
      *nprioexternbins = branchcand->nprioexternbins;
   if( nprioexternints != NULL )
      *nprioexternints = branchcand->nprioexternints;
   if( nprioexternimpls != NULL )
      *nprioexternimpls = branchcand->nprioexternimpls;

   return SCIP_OKAY;
}

/** gets maximal branching priority of LP branching candidates */
int SCIPbranchcandGetLPMaxPrio(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->lpmaxpriority;
}

/** gets number of LP branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioLPCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriolpcands;
}

/** gets maximal branching priority of external branching candidates */
int SCIPbranchcandGetExternMaxPrio(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->externmaxpriority;
}

/** gets number of external branching candidates */
int SCIPbranchcandGetNExternCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->nexterncands;
}

/** gets number of external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->nprioexterncands;
}

/** gets number of binary external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternBins(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->nprioexternbins;
}

/** gets number of integer external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternInts(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->nprioexternints;
}

/** gets number of implicit integer external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternImpls(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->nprioexternimpls;
}

/** gets number of continuous external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternConts(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->nprioexterncands - branchcand->nprioexternbins - branchcand->nprioexternints - branchcand->nprioexternimpls;
}

/** insert variable, its score and its solution value into the external branching candidate storage
 * the absolute difference of the current lower and upper bounds of the variable must be at least epsilon
 */
SCIP_RETCODE SCIPbranchcandAddExternCand(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to insert */
   SCIP_Real             score,              /**< score of external candidate, e.g. infeasibility */
   SCIP_Real             solval              /**< value of the variable in the current solution */
   )
{
   SCIP_VARTYPE vartype;
   int branchpriority;
   int insertpos;

   assert(branchcand != NULL);
   assert(var != NULL);
   assert(!SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))); /* the variable should not be fixed yet */
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || !SCIPsetIsEQ(set, SCIPsetCeil(set, SCIPvarGetLbLocal(var)), SCIPsetFloor(set, SCIPvarGetUbLocal(var)))); /* a discrete variable should also not be fixed, even after rounding bounds to integral values */
   assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR || !SCIPsetIsEQ(set, SCIPvarGetMultaggrLbLocal(var, set), SCIPvarGetMultaggrUbLocal(var, set))); /* also the current bounds of a multi-aggregated variable should not be fixed yet */
   assert(branchcand->nprioexterncands <= branchcand->nexterncands);
   assert(branchcand->nexterncands <= branchcand->externcandssize);

   vartype = SCIPvarGetType(var);
   branchpriority = SCIPvarGetBranchPriority(var);
   insertpos = branchcand->nexterncands;

   SCIP_CALL( ensureExterncandsSize(branchcand, set, branchcand->nexterncands+1) );

   SCIPsetDebugMsg(set, "inserting external candidate <%s> of type %d and priority %d into candidate set (maxprio: %d), score = %g, solval = %g\n",
      SCIPvarGetName(var), vartype, branchpriority, branchcand->externmaxpriority, score, solval);

   /* insert the variable into externcands, making sure, that the highest priority candidates are at the front
    * and ordered binaries, integers, implicit integers, continuous
    */
   if( branchpriority > branchcand->externmaxpriority )
   {
      /* candidate has higher priority than the current maximum:
       * move it to the front and declare it to be the single best candidate
       */
      branchcand->externcands[insertpos] = branchcand->externcands[0];
      branchcand->externcandsscore[insertpos] = branchcand->externcandsscore[0];
      branchcand->externcandssol[insertpos] = branchcand->externcandssol[0];

      insertpos = 0;

      branchcand->nprioexterncands = 1;
      branchcand->nprioexternbins = (vartype == SCIP_VARTYPE_BINARY ? 1 : 0);
      branchcand->nprioexternints = (vartype == SCIP_VARTYPE_INTEGER ? 1 : 0);
      branchcand->nprioexternimpls = (vartype == SCIP_VARTYPE_IMPLINT ? 1 : 0);
      branchcand->externmaxpriority = branchpriority;
   }
   else if( branchpriority == branchcand->externmaxpriority )
   {
      /* candidate has equal priority as the current maximum:
       * move away the first non-maximal priority candidate, move the current candidate to the correct
       * slot (binaries first, integers next, implicit integers next, continuous last) and increase the number 
       * of maximal priority candidates
       */
      if( insertpos != branchcand->nprioexterncands )
      {
         branchcand->externcands[insertpos] = branchcand->externcands[branchcand->nprioexterncands];
         branchcand->externcandsscore[insertpos] = branchcand->externcandsscore[branchcand->nprioexterncands];
         branchcand->externcandssol[insertpos] = branchcand->externcandssol[branchcand->nprioexterncands];

         insertpos = branchcand->nprioexterncands;
      }
      branchcand->nprioexterncands++;
      if( vartype == SCIP_VARTYPE_BINARY || vartype == SCIP_VARTYPE_INTEGER || vartype == SCIP_VARTYPE_IMPLINT )
      {
         if( insertpos != branchcand->nprioexternbins + branchcand->nprioexternints + branchcand->nprioexternimpls )
         {
            branchcand->externcands[insertpos] =
               branchcand->externcands[branchcand->nprioexternbins + branchcand->nprioexternints + branchcand->nprioexternimpls];
            branchcand->externcandsscore[insertpos] =
               branchcand->externcandsscore[branchcand->nprioexternbins + branchcand->nprioexternints + branchcand->nprioexternimpls];
            branchcand->externcandssol[insertpos] =
               branchcand->externcandssol[branchcand->nprioexternbins + branchcand->nprioexternints + branchcand->nprioexternimpls];

            insertpos = branchcand->nprioexternbins + branchcand->nprioexternints + branchcand->nprioexternimpls;
         }
         branchcand->nprioexternimpls++;


         if( vartype == SCIP_VARTYPE_BINARY || vartype == SCIP_VARTYPE_INTEGER )
         {
            if( insertpos != branchcand->nprioexternbins + branchcand->nprioexternints )
            {
               branchcand->externcands[insertpos] = 
                  branchcand->externcands[branchcand->nprioexternbins + branchcand->nprioexternints];
               branchcand->externcandsscore[insertpos] = 
                  branchcand->externcandsscore[branchcand->nprioexternbins + branchcand->nprioexternints];
               branchcand->externcandssol[insertpos] = 
                  branchcand->externcandssol[branchcand->nprioexternbins + branchcand->nprioexternints];

               insertpos = branchcand->nprioexternbins + branchcand->nprioexternints;
            }
            branchcand->nprioexternints++;
            branchcand->nprioexternimpls--;


            if( vartype == SCIP_VARTYPE_BINARY )
            {
               if( insertpos != branchcand->nprioexternbins )
               {
                  branchcand->externcands[insertpos] = branchcand->externcands[branchcand->nprioexternbins];
                  branchcand->externcandsscore[insertpos] = branchcand->externcandsscore[branchcand->nprioexternbins];
                  branchcand->externcandssol[insertpos] = branchcand->externcandssol[branchcand->nprioexternbins];

                  insertpos = branchcand->nprioexternbins;
               }
               branchcand->nprioexternbins++;
               branchcand->nprioexternints--;
            }
         }
      }
   }
   branchcand->externcands[insertpos] = var;
   branchcand->externcandsscore[insertpos] = score;
   branchcand->externcandssol[insertpos] = solval;
   branchcand->nexterncands++;

   SCIPsetDebugMsg(set, " -> inserted at position %d (nprioexterncands=%d)\n", insertpos, branchcand->nprioexterncands);

   assert(0 <= branchcand->nprioexterncands && branchcand->nprioexterncands <= branchcand->nexterncands);
   assert(0 <= branchcand->nprioexternbins && branchcand->nprioexternbins <= branchcand->nprioexterncands);
   assert(0 <= branchcand->nprioexternints && branchcand->nprioexternints <= branchcand->nprioexterncands);
   assert(0 <= branchcand->nprioexternimpls && branchcand->nprioexternimpls <= branchcand->nprioexterncands);

   return SCIP_OKAY;
}

/** removes all external candidates from the storage for external branching */
void SCIPbranchcandClearExternCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   branchcand->nexterncands = 0;
   branchcand->nprioexterncands = 0;
   branchcand->nprioexternbins = 0;
   branchcand->nprioexternints = 0;
   branchcand->nprioexternimpls = 0;
   branchcand->externmaxpriority = INT_MIN;
}

/** checks whether the given variable is contained in the candidate storage for external branching */
SCIP_Bool SCIPbranchcandContainsExternCand(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR*             var                 /**< variable to look for */
   )
{
   int branchpriority;
   int i;

   assert(branchcand != NULL);
   assert(var != NULL);
   assert(branchcand->nprioexterncands <= branchcand->nexterncands);
   assert(branchcand->nexterncands <= branchcand->externcandssize);

   branchpriority = SCIPvarGetBranchPriority(var);

   /* look for the variable in the externcands, using the fact, that the highest priority candidates are at the front
    * and ordered binaries, integers, implicit integers, continuous
    */
   if( branchpriority > branchcand->externmaxpriority )
   {
      /* the branching priority of the variable is higher than the maximal priority contained in the array;
       * the variable can thus not be contained */
      return FALSE;
   }
   if( branchpriority == branchcand->externmaxpriority )
   {
      SCIP_VARTYPE vartype;

      vartype = SCIPvarGetType(var);

      /* variable has equal priority as the current maximum:
       * look for it in the correct slot (binaries first, integers next, implicit integers next, continuous last)
       */
      if( vartype == SCIP_VARTYPE_BINARY )
      {
         /* the variable is binary, look at the first branchcand->nprioexternbins slots */
         for( i = 0; i < branchcand->nprioexternbins; i++ )
            if( branchcand->externcands[i] == var )
               return TRUE;
         return FALSE;
      }
      if( vartype == SCIP_VARTYPE_INTEGER )
      {
         /* the variable is integer, look at the slots containing integers */
         for( i = 0; i < branchcand->nprioexternints; i++ )
            if( branchcand->externcands[branchcand->nprioexternbins + i] == var )
               return TRUE;
         return FALSE;
      }
      if( vartype == SCIP_VARTYPE_IMPLINT )
      {
         /* the variable is implicit integer, look at the slots containing implicit integers */
         for( i = 0; i < branchcand->nprioexternimpls; i++ )
            if( branchcand->externcands[branchcand->nprioexternbins + branchcand->nprioexternints + i] == var )
               return TRUE;
         return FALSE;
      }
      /* the variable is continuous, look at the slots containing continuous variables */
      assert(vartype == SCIP_VARTYPE_CONTINUOUS);
      for( i = branchcand->nprioexternbins + branchcand->nprioexternints + branchcand->nprioexternimpls; 
           i < branchcand->nprioexterncands; i++ )
         if( branchcand->externcands[i] == var )
            return TRUE;
      return FALSE;
   }
   /* the branching priority of the variable is lower than the maximal priority contained in the array;
    * look for the variable in the candidates which do not have maximal priority */
   assert(branchpriority < branchcand->externmaxpriority);

   for( i = branchcand->nprioexterncands; i < branchcand->nexterncands; i++ )
      if( branchcand->externcands[i] == var )
         return TRUE;
   return FALSE;
}

/** gets branching candidates for pseudo solution branching (non-fixed variables) */
SCIP_RETCODE SCIPbranchcandGetPseudoCands(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*                  npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*                  npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   assert(branchcand != NULL);

#ifndef NDEBUG
   /* check, if the current pseudo branching candidate array is correct */
   {
      SCIP_VAR* var;
      int npcs;
      int v;

      assert(prob != NULL);

      /* pseudo branching candidates are non-fixed binary, integer, and implicit integer variables */
      npcs = 0;
      for( v = 0; v < prob->nbinvars + prob->nintvars + prob->nimplvars; ++v )
      {
         var = prob->vars[v];
         assert(var != NULL);
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY
            || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER
            || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
         assert(SCIPsetIsFeasIntegral(set, SCIPvarGetLbLocal(var)));
         assert(SCIPsetIsFeasIntegral(set, SCIPvarGetUbLocal(var)));
         assert(SCIPsetIsLE(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

         if( SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            assert(0 <= var->pseudocandindex && var->pseudocandindex < branchcand->npseudocands);
            assert(branchcand->pseudocands[var->pseudocandindex] == var);
            npcs++;
         }
         else
         {
            assert(var->pseudocandindex == -1);
         }
      }
      assert(branchcand->npseudocands == npcs);
      for (v = 0; v < branchcand->npriopseudocands; ++v)
         assert( branchcand->pseudocands[v]->branchpriority == branchcand->pseudomaxpriority );
   }
#endif

   /* assign return values */
   if( pseudocands != NULL )
      *pseudocands = branchcand->pseudocands;
   if( npseudocands != NULL )
      *npseudocands = branchcand->npseudocands;
   if( npriopseudocands != NULL )
      *npriopseudocands = (set->branch_preferbinary && branchcand->npriopseudobins > 0 ? branchcand->npriopseudobins
         : branchcand->npriopseudocands);

   return SCIP_OKAY;
}

/** gets number of branching candidates for pseudo solution branching (non-fixed variables) */
int SCIPbranchcandGetNPseudoCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npseudocands;
}

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriopseudocands;
}

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoBins(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriopseudobins;
}

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoInts(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriopseudoints;
}

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoImpls(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriopseudocands - branchcand->npriopseudobins - branchcand->npriopseudoints;
}

/** insert pseudocand at given position, or to the first positions of the maximal priority candidates, using the
 *  given position as free slot for the other candidates
 */
static
void branchcandInsertPseudoCand(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR*             var,                /**< variable to insert */
   int                   insertpos           /**< free position to insert the variable */
   )
{
   SCIP_VARTYPE vartype;
   int branchpriority;

   assert(branchcand != NULL);
   assert(var != NULL);
   assert(branchcand->npriopseudocands <= insertpos && insertpos < branchcand->npseudocands);
   assert(branchcand->npseudocands <= branchcand->pseudocandssize);

   vartype = SCIPvarGetType(var);
   branchpriority = SCIPvarGetBranchPriority(var);

   SCIPdebugMessage("inserting pseudo candidate <%s> of type %d and priority %d into candidate set at position %d (maxprio: %d)\n",
      SCIPvarGetName(var), vartype, branchpriority, insertpos, branchcand->pseudomaxpriority);

   /* insert the variable into pseudocands, making sure, that the highest priority candidates are at the front
    * and ordered binaries, integers, implicit integers
    */
   if( branchpriority > branchcand->pseudomaxpriority )
   {
      /* candidate has higher priority than the current maximum:
       * move it to the front and declare it to be the single best candidate
       */
      if( insertpos != 0 )
      {
         branchcand->pseudocands[insertpos] = branchcand->pseudocands[0];
         branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
         insertpos = 0;
      }
      branchcand->npriopseudocands = 1;
      branchcand->npriopseudobins = (vartype == SCIP_VARTYPE_BINARY ? 1 : 0);
      branchcand->npriopseudoints = (vartype == SCIP_VARTYPE_INTEGER ? 1 : 0);
      branchcand->pseudomaxpriority = branchpriority;
   }
   else if( branchpriority == branchcand->pseudomaxpriority )
   {
      /* candidate has equal priority as the current maximum:
       * move away the first non-maximal priority candidate, move the current candidate to the correct
       * slot (binaries first, integers next, implicit integers last) and increase the number of maximal priority candidates
       */
      if( insertpos != branchcand->npriopseudocands )
      {
         branchcand->pseudocands[insertpos] = branchcand->pseudocands[branchcand->npriopseudocands];
         branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
         insertpos = branchcand->npriopseudocands;
      }
      branchcand->npriopseudocands++;
      if( vartype == SCIP_VARTYPE_BINARY || vartype == SCIP_VARTYPE_INTEGER )
      {
         if( insertpos != branchcand->npriopseudobins + branchcand->npriopseudoints )
         {
            branchcand->pseudocands[insertpos] =
               branchcand->pseudocands[branchcand->npriopseudobins + branchcand->npriopseudoints];
            branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
            insertpos = branchcand->npriopseudobins + branchcand->npriopseudoints;
         }
         branchcand->npriopseudoints++;

         if( vartype == SCIP_VARTYPE_BINARY )
         {
            if( insertpos != branchcand->npriopseudobins )
            {
               branchcand->pseudocands[insertpos] = branchcand->pseudocands[branchcand->npriopseudobins];
               branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
               insertpos = branchcand->npriopseudobins;
            }
            branchcand->npriopseudobins++;
            branchcand->npriopseudoints--;
         }
      }
   }
   branchcand->pseudocands[insertpos] = var;
   var->pseudocandindex = insertpos;

   SCIPdebugMessage(" -> inserted at position %d (npriopseudocands=%d)\n", insertpos, branchcand->npriopseudocands);

   assert(0 <= branchcand->npriopseudocands && branchcand->npriopseudocands <= branchcand->npseudocands);
   assert(0 <= branchcand->npriopseudobins && branchcand->npriopseudobins <= branchcand->npriopseudocands);
   assert(0 <= branchcand->npriopseudoints && branchcand->npriopseudoints <= branchcand->npriopseudocands);
}

/** sorts the pseudo branching candidates, such that the candidates of maximal priority are at the front,
 *  ordered by binaries, integers, implicit integers
 */
static
void branchcandSortPseudoCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   SCIP_VAR* var;
   int i;

   assert(branchcand != NULL);
   assert(branchcand->npriopseudocands == 0); /* is only be called after removal of last maximal candidate */
   assert(branchcand->npriopseudobins == 0);
   assert(branchcand->npriopseudoints == 0);

   SCIPdebugMessage("resorting pseudo candidates\n");

   branchcand->pseudomaxpriority = INT_MIN;

   for( i = 0; i < branchcand->npseudocands; ++i )
   {
      var = branchcand->pseudocands[i];
      assert(var->pseudocandindex == i);

      if( SCIPvarGetBranchPriority(var) >= branchcand->pseudomaxpriority )
         branchcandInsertPseudoCand(branchcand, var, i);
   }

   assert(0 <= branchcand->npriopseudocands && branchcand->npriopseudocands <= branchcand->npseudocands);
   assert(0 <= branchcand->npriopseudobins && branchcand->npriopseudobins <= branchcand->npriopseudocands);
   assert(0 <= branchcand->npriopseudoints && branchcand->npriopseudoints <= branchcand->npriopseudocands);
#ifndef NDEBUG
   {
      for (i = 0; i < branchcand->npriopseudocands; ++i)
         assert( branchcand->pseudocands[i]->branchpriority == branchcand->pseudomaxpriority );
   }
#endif
}

/** removes pseudo candidate from pseudocands array
 */
static
void branchcandRemovePseudoCand(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR*             var                 /**< variable to remove */
   )
{
   int freepos;

   assert(branchcand != NULL);
   assert(var != NULL);
   assert(var->pseudocandindex < branchcand->npseudocands);
   assert(branchcand->pseudocands[var->pseudocandindex] == var);
   assert(branchcand->pseudocands[branchcand->npseudocands-1] != NULL);

   /* Note that the branching priority of the variable to be removed is not necessarily equal to pseudomaxpriority, since
    * the status of the variable might have changed, leading to a change in the branching priority. Moreover, if the
    * variable was part of an aggregation, even other variables might at this point have different priorities. */
   SCIPdebugMessage("removing pseudo candidate <%s> of type %d and priority %d at %d from candidate set (maxprio: %d)\n",
      SCIPvarGetName(var), SCIPvarGetType(var), SCIPvarGetBranchPriority(var), var->pseudocandindex,
      branchcand->pseudomaxpriority);

   /* delete the variable from pseudocands, making sure, that the highest priority candidates are at the front
    * and ordered binaries, integers, implicit integers
    */
   freepos = var->pseudocandindex;
   var->pseudocandindex = -1;
   assert(0 <= freepos && freepos < branchcand->npseudocands);

   if( freepos < branchcand->npriopseudobins )
   {
      /* a binary candidate of maximal priority was removed */
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      if( freepos != branchcand->npriopseudobins - 1 )
      {
         branchcand->pseudocands[freepos] = branchcand->pseudocands[branchcand->npriopseudobins - 1];
         branchcand->pseudocands[freepos]->pseudocandindex = freepos;
         freepos = branchcand->npriopseudobins - 1;
      }
      branchcand->npriopseudobins--;
      branchcand->npriopseudoints++;
   }

   if( freepos < branchcand->npriopseudobins + branchcand->npriopseudoints )
   {
      /* a binary or integer candidate of maximal priority was removed */
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);
      if( freepos != branchcand->npriopseudobins + branchcand->npriopseudoints - 1 )
      {
         branchcand->pseudocands[freepos] =
            branchcand->pseudocands[branchcand->npriopseudobins + branchcand->npriopseudoints - 1];
         branchcand->pseudocands[freepos]->pseudocandindex = freepos;
         freepos = branchcand->npriopseudobins + branchcand->npriopseudoints - 1;
      }
      branchcand->npriopseudoints--;
   }

   if( freepos < branchcand->npriopseudocands )
   {
      /* a candidate of maximal priority was removed */
      if( freepos != branchcand->npriopseudocands - 1 )
      {
         branchcand->pseudocands[freepos] = branchcand->pseudocands[branchcand->npriopseudocands - 1];
         branchcand->pseudocands[freepos]->pseudocandindex = freepos;
         freepos = branchcand->npriopseudocands - 1;
      }
      branchcand->npriopseudocands--;
   }
   if( freepos != branchcand->npseudocands - 1 )
   {
      branchcand->pseudocands[freepos] = branchcand->pseudocands[branchcand->npseudocands - 1];
      branchcand->pseudocands[freepos]->pseudocandindex = freepos;
   }
   branchcand->npseudocands--;

   assert(0 <= branchcand->npriopseudocands && branchcand->npriopseudocands <= branchcand->npseudocands);
   assert(0 <= branchcand->npriopseudobins && branchcand->npriopseudobins <= branchcand->npriopseudocands);
   assert(0 <= branchcand->npriopseudoints && branchcand->npriopseudoints <= branchcand->npriopseudocands);

   /* if all maximal priority candidates were removed, resort the array s.t. the new maximal priority candidates
    * are at the front
    */
   if( branchcand->npriopseudocands == 0 )
      branchcandSortPseudoCands(branchcand);
}

/** removes variable from branching candidate list */
SCIP_RETCODE SCIPbranchcandRemoveVar(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR*             var                 /**< variable that changed its bounds */
   )
{
   assert(var != NULL);

   /* make sure, variable is not member of the pseudo branching candidate list */
   if( var->pseudocandindex >= 0 )
   {
      branchcandRemovePseudoCand(branchcand, var);
   }

   return SCIP_OKAY;
}

/** updates branching candidate list for a given variable */
SCIP_RETCODE SCIPbranchcandUpdateVar(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that changed its bounds */
   )
{
   assert(branchcand != NULL);
   assert(var != NULL);

   if( (SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN)
      && SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS
      && SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      /* variable is neither continuous nor fixed and has non-empty domain: make sure it is member of the pseudo branching candidate list */
      if( var->pseudocandindex == -1 )
      {
         SCIP_CALL( ensurePseudocandsSize(branchcand, set, branchcand->npseudocands+1) );

         branchcand->npseudocands++;
         branchcandInsertPseudoCand(branchcand, var, branchcand->npseudocands-1);
      }
   }
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED
         || SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
         || SCIPsetIsGE(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

      /* variable is continuous or fixed or has empty domain: make sure it is not member of the pseudo branching candidate list */
      SCIP_CALL( SCIPbranchcandRemoveVar(branchcand, var) );
   }

   return SCIP_OKAY;
}

/** updates branching priority of the given variable and update the pseudo candidate array if needed */
SCIP_RETCODE SCIPbranchcandUpdateVarBranchPriority(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable that changed its bounds */
   int                   branchpriority      /**< branch priority of the variable */
   )
{
   int oldbranchpriority;
   int pseudomaxpriority;

   assert(branchcand != NULL);

   oldbranchpriority = SCIPvarGetBranchPriority(var);

   if( oldbranchpriority == branchpriority )
      return SCIP_OKAY;

   pseudomaxpriority = branchcand->pseudomaxpriority;

   /* if the variable currently belongs to the priority set or the new branching priority is larger than the current one,
    * remove it from the pseudo branch candidate array */
   if( oldbranchpriority == pseudomaxpriority || branchpriority > pseudomaxpriority )
   {
      SCIP_CALL( SCIPbranchcandRemoveVar(branchcand, var) );
      assert(var->pseudocandindex == -1);
   }

   /* change the branching priority of the variable */
   SCIP_CALL( SCIPvarChgBranchPriority(var, branchpriority) );

   /* if the variable is not part of the pseudo branching candidate array, check if it is a pseudo branching candidate
    * and add it if so */
   SCIP_CALL( SCIPbranchcandUpdateVar(branchcand, set, var) );

   return SCIP_OKAY;
}



/*
 * branching rule methods
 */

/** compares two branching rules w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPbranchruleComp)
{  /*lint --e{715}*/
   return ((SCIP_BRANCHRULE*)elem2)->priority - ((SCIP_BRANCHRULE*)elem1)->priority;
}

/** comparison method for sorting branching rules w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPbranchruleCompName)
{
   return strcmp(SCIPbranchruleGetName((SCIP_BRANCHRULE*)elem1), SCIPbranchruleGetName((SCIP_BRANCHRULE*)elem2));
}

/** method to call, when the priority of a branching rule was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdBranchrulePriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetBranchrulePriority() to mark the branchrules unsorted */
   SCIP_CALL( SCIPsetBranchrulePriority(scip, (SCIP_BRANCHRULE*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given branchrule to a new scip */
SCIP_RETCODE SCIPbranchruleCopyInclude(
   SCIP_BRANCHRULE*      branchrule,         /**< branchrule */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( branchrule->branchcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including branching rule %s in subscip %p\n", SCIPbranchruleGetName(branchrule), (void*)set->scip);
      SCIP_CALL( branchrule->branchcopy(set->scip, branchrule) );
   }

   return SCIP_OKAY;
}

/** creates a branching rule */
SCIP_RETCODE SCIPbranchruleCreate(
   SCIP_BRANCHRULE**     branchrule,         /**< pointer to store branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of branching rule */
   const char*           desc,               /**< description of branching rule */
   int                   priority,           /**< priority of the branching rule */
   int                   maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_DECL_BRANCHCOPY  ((*branchcopy)),    /**< copy method of branching rule */
   SCIP_DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol)),/**< solving process initialization method of branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol)),/**< solving process deinitialization method of branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   SCIP_DECL_BRANCHEXECEXT((*branchexecext)),/**< branching execution method for external solutions */
   SCIP_DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(branchrule != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_ALLOC( BMSallocMemory(branchrule) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*branchrule)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*branchrule)->desc, desc, strlen(desc)+1) );
   (*branchrule)->priority = priority;
   (*branchrule)->maxdepth = maxdepth;
   (*branchrule)->maxbounddist = maxbounddist;
   (*branchrule)->branchcopy = branchcopy;
   (*branchrule)->branchfree = branchfree;
   (*branchrule)->branchinit = branchinit;
   (*branchrule)->branchexit = branchexit;
   (*branchrule)->branchinitsol = branchinitsol;
   (*branchrule)->branchexitsol = branchexitsol;
   (*branchrule)->branchexeclp = branchexeclp;
   (*branchrule)->branchexecext = branchexecext;
   (*branchrule)->branchexecps = branchexecps;
   (*branchrule)->branchruledata = branchruledata;
   SCIP_CALL( SCIPclockCreate(&(*branchrule)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*branchrule)->branchclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*branchrule)->nlpcalls = 0;
   (*branchrule)->nexterncalls = 0;
   (*branchrule)->npseudocalls = 0;
   (*branchrule)->ncutoffs = 0;
   (*branchrule)->ncutsfound = 0;
   (*branchrule)->nconssfound = 0;
   (*branchrule)->ndomredsfound = 0;
   (*branchrule)->nchildren = 0;
   (*branchrule)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "branching/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of branching rule <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*branchrule)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
         paramChgdBranchrulePriority, (SCIP_PARAMDATA*)(*branchrule)) ); /*lint !e740*/
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "branching/%s/maxdepth", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "maximal depth level, up to which branching rule <%s> should be used (-1 for no limit)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*branchrule)->maxdepth, FALSE, maxdepth, -1, SCIP_MAXTREEDEPTH,
         NULL, NULL) ); /*lint !e740*/
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "branching/%s/maxbounddist", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying branching rule (0.0: only on current best node, 1.0: on all nodes)");
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*branchrule)->maxbounddist, FALSE, maxbounddist, 0.0, 1.0,
         NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** frees memory of branching rule */
SCIP_RETCODE SCIPbranchruleFree(
   SCIP_BRANCHRULE**     branchrule,         /**< pointer to branching rule data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(branchrule != NULL);
   assert(*branchrule != NULL);
   assert(!(*branchrule)->initialized);
   assert(set != NULL);

   /* call destructor of branching rule */
   if( (*branchrule)->branchfree != NULL )
   {
      SCIP_CALL( (*branchrule)->branchfree(set->scip, *branchrule) );
   }

   SCIPclockFree(&(*branchrule)->branchclock);
   SCIPclockFree(&(*branchrule)->setuptime);
   BMSfreeMemoryArray(&(*branchrule)->name);
   BMSfreeMemoryArray(&(*branchrule)->desc);
   BMSfreeMemory(branchrule);

   return SCIP_OKAY;
}

/** initializes branching rule */
SCIP_RETCODE SCIPbranchruleInit(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);

   if( branchrule->initialized )
   {
      SCIPerrorMessage("branching rule <%s> already initialized\n", branchrule->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(branchrule->setuptime);
      SCIPclockReset(branchrule->branchclock);
      branchrule->nlpcalls = 0;
      branchrule->nexterncalls = 0;
      branchrule->npseudocalls = 0;
      branchrule->ncutoffs = 0;
      branchrule->ncutsfound = 0;
      branchrule->nconssfound = 0;
      branchrule->ndomredsfound = 0;
      branchrule->nchildren = 0;
   }

   if( branchrule->branchinit != NULL )
   {
      /* start timing */
      SCIPclockStart(branchrule->setuptime, set);

      SCIP_CALL( branchrule->branchinit(set->scip, branchrule) );

      /* stop timing */
      SCIPclockStop(branchrule->setuptime, set);
   }
   branchrule->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes branching rule */
SCIP_RETCODE SCIPbranchruleExit(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);

   if( !branchrule->initialized )
   {
      SCIPerrorMessage("branching rule <%s> not initialized\n", branchrule->name);
      return SCIP_INVALIDCALL;
   }

   if( branchrule->branchexit != NULL )
   {
      /* start timing */
      SCIPclockStart(branchrule->setuptime, set);

      SCIP_CALL( branchrule->branchexit(set->scip, branchrule) );

      /* stop timing */
      SCIPclockStop(branchrule->setuptime, set);
   }
   branchrule->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs branching rule that the branch and bound process is being started */
SCIP_RETCODE SCIPbranchruleInitsol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);

   /* call solving process initialization method of branching rule */
   if( branchrule->branchinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(branchrule->setuptime, set);

      SCIP_CALL( branchrule->branchinitsol(set->scip, branchrule) );

      /* stop timing */
      SCIPclockStop(branchrule->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs branching rule that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbranchruleExitsol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of branching rule */
   if( branchrule->branchexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(branchrule->setuptime, set);

      SCIP_CALL( branchrule->branchexitsol(set->scip, branchrule) );

      /* stop timing */
      SCIPclockStop(branchrule->setuptime, set);
   }

   return SCIP_OKAY;
}

/** executes branching rule for fractional LP solution */
SCIP_RETCODE SCIPbranchruleExecLPSol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->focusnode != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexeclp != NULL
      && (branchrule->maxdepth == -1 || branchrule->maxdepth >= SCIPtreeGetCurrentDepth(tree)) )
   {
      SCIP_Real loclowerbound;
      SCIP_Real glblowerbound;
      SCIP_Bool runbranchrule;

      loclowerbound = SCIPnodeGetLowerbound(tree->focusnode);
      glblowerbound = SCIPtreeGetLowerbound(tree, set);

      /* we distinguish between finite and infinite global lower bounds to avoid comparisons between different values > SCIPinfinity() */
      if( SCIPsetIsInfinity(set, -glblowerbound) )
         runbranchrule = SCIPsetIsInfinity(set, -loclowerbound) || SCIPsetIsGE(set, branchrule->maxbounddist, 1.0);
      else
      {
         assert(!SCIPsetIsInfinity(set, -loclowerbound));
         runbranchrule = SCIPsetIsLE(set, loclowerbound - glblowerbound, branchrule->maxbounddist * (cutoffbound - glblowerbound));
      }

      if( runbranchrule )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         SCIP_Longint oldnactiveconss;
         int oldncuts;

         SCIPsetDebugMsg(set, "executing LP branching rule <%s>\n", branchrule->name);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
         oldncuts = SCIPsepastoreGetNCuts(sepastore);
         oldnactiveconss = stat->nactiveconssadded;

         /* start timing */
         SCIPclockStart(branchrule->branchclock, set);

         /* call external method */
         SCIP_CALL( branchrule->branchexeclp(set->scip, branchrule, allowaddcons, result) );

         /* stop timing */
         SCIPclockStop(branchrule->branchclock, set);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_BRANCHED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN )
         {
            SCIPerrorMessage("branching rule <%s> returned invalid result code <%d> from LP solution branching\n",
               branchrule->name, *result);
            return SCIP_INVALIDRESULT;
         }
         if( *result == SCIP_CONSADDED && !allowaddcons )
         {
            SCIPerrorMessage("branching rule <%s> added a constraint in LP solution branching without permission\n",
               branchrule->name);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            branchrule->nlpcalls++;
         if( *result == SCIP_CUTOFF )
            branchrule->ncutoffs++;
         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);

            /* update domain reductions; therefore remove the domain
             * reduction counts which were generated in probing mode */
            branchrule->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            branchrule->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

            branchrule->ncutsfound += SCIPsepastoreGetNCuts(sepastore) - oldncuts; /*lint !e776*/
            branchrule->nconssfound += stat->nactiveconssadded - oldnactiveconss; /*lint !e776*/
         }
         else
            branchrule->nchildren += tree->nchildren;
      }
   }

   return SCIP_OKAY;
}

/** executes branching rule for external branching candidates */
SCIP_RETCODE SCIPbranchruleExecExternSol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->focusnode != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexecext != NULL
      && (branchrule->maxdepth == -1 || branchrule->maxdepth >= SCIPtreeGetCurrentDepth(tree)) )
   {
      SCIP_Real loclowerbound;
      SCIP_Real glblowerbound;
      SCIP_Bool runbranchrule;

      loclowerbound = SCIPnodeGetLowerbound(tree->focusnode);
      glblowerbound = SCIPtreeGetLowerbound(tree, set);
      assert(!SCIPsetIsInfinity(set, loclowerbound));

      /* we distinguish between finite and infinite global lower bounds to avoid comparisons between different values > SCIPinfinity() */
      if( SCIPsetIsInfinity(set, -glblowerbound) )
         runbranchrule = SCIPsetIsInfinity(set, -loclowerbound) || SCIPsetIsGE(set, branchrule->maxbounddist, 1.0);
      else
      {
         assert(!SCIPsetIsInfinity(set, -loclowerbound));
         runbranchrule = SCIPsetIsLE(set, loclowerbound - glblowerbound, branchrule->maxbounddist * (cutoffbound - glblowerbound));
      }

      if( runbranchrule )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         int oldncuts;
         int oldnactiveconss;

         SCIPsetDebugMsg(set, "executing external solution branching rule <%s>\n", branchrule->name);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
         oldncuts = SCIPsepastoreGetNCuts(sepastore);
         oldnactiveconss = stat->nactiveconss;

         /* start timing */
         SCIPclockStart(branchrule->branchclock, set);

         /* call external method */
         SCIP_CALL( branchrule->branchexecext(set->scip, branchrule, allowaddcons, result) );

         /* stop timing */
         SCIPclockStop(branchrule->branchclock, set);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_BRANCHED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN )
         {
            SCIPerrorMessage("branching rule <%s> returned invalid result code <%d> from external solution branching\n",
               branchrule->name, *result);
            return SCIP_INVALIDRESULT;
         }
         if( *result == SCIP_CONSADDED && !allowaddcons )
         {
            SCIPerrorMessage("branching rule <%s> added a constraint in external solution branching without permission\n",
               branchrule->name);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            branchrule->nexterncalls++;
         if( *result == SCIP_CUTOFF )
            branchrule->ncutoffs++;
         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);

            /* update domain reductions; therefore remove the domain
             * reduction counts which were generated in probing mode */
            branchrule->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            branchrule->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

            branchrule->ncutsfound += SCIPsepastoreGetNCuts(sepastore) - oldncuts; /*lint !e776*/
            branchrule->nconssfound += stat->nactiveconss - oldnactiveconss; /*lint !e776*/
         }
         else
            branchrule->nchildren += tree->nchildren;
      }
   }
   return SCIP_OKAY;
}

/** executes branching rule for not completely fixed pseudo solution */
SCIP_RETCODE SCIPbranchruleExecPseudoSol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexecps != NULL
      && (branchrule->maxdepth == -1 || branchrule->maxdepth >= SCIPtreeGetCurrentDepth(tree)) )
   {
      SCIP_Real loclowerbound;
      SCIP_Real glblowerbound;
      SCIP_Bool runbranchrule;

      loclowerbound = SCIPnodeGetLowerbound(tree->focusnode);
      glblowerbound = SCIPtreeGetLowerbound(tree, set);

      /* we distinguish between finite and infinite global lower bounds to avoid comparisons between different values > SCIPinfinity() */
      if( SCIPsetIsInfinity(set, -glblowerbound) )
         runbranchrule = SCIPsetIsInfinity(set, -loclowerbound) || SCIPsetIsGE(set, branchrule->maxbounddist, 1.0);
      else
      {
         assert(!SCIPsetIsInfinity(set, -loclowerbound));
         runbranchrule = SCIPsetIsLE(set, loclowerbound - glblowerbound, branchrule->maxbounddist * (cutoffbound - glblowerbound));
      }

      if( runbranchrule )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         SCIP_Longint oldnactiveconss;

         SCIPsetDebugMsg(set, "executing pseudo branching rule <%s>\n", branchrule->name);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
         oldnactiveconss = stat->nactiveconss;

         /* start timing */
         SCIPclockStart(branchrule->branchclock, set);

         /* call external method */
         SCIP_CALL( branchrule->branchexecps(set->scip, branchrule, allowaddcons, result) );

         /* stop timing */
         SCIPclockStop(branchrule->branchclock, set);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_BRANCHED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN )
         {
            SCIPerrorMessage("branching rule <%s> returned invalid result code <%d> from pseudo solution branching\n",
               branchrule->name, *result);
            return SCIP_INVALIDRESULT;
         }
         if( *result == SCIP_CONSADDED && !allowaddcons )
         {
            SCIPerrorMessage("branching rule <%s> added a constraint in pseudo solution branching without permission\n",
               branchrule->name);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            branchrule->npseudocalls++;
         if( *result == SCIP_CUTOFF )
            branchrule->ncutoffs++;
         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);

            /* update domain reductions; therefore remove the domain
             * reduction counts which were generated in probing mode */
            branchrule->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            branchrule->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

            branchrule->nconssfound += stat->nactiveconss - oldnactiveconss;
         }
         else
            branchrule->nchildren += tree->nchildren;
      }
   }

   return SCIP_OKAY;
}

/** gets user data of branching rule */
SCIP_BRANCHRULEDATA* SCIPbranchruleGetData(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->branchruledata;
}

/** sets user data of branching rule; user has to free old data in advance! */
void SCIPbranchruleSetData(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< new branching rule user data */
   )
{
   assert(branchrule != NULL);

   branchrule->branchruledata = branchruledata;
}

/** sets copy method of branching rule */
void SCIPbranchruleSetCopy(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHCOPY  ((*branchcopy))     /**< copy method of branching rule or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(branchrule != NULL);

   branchrule->branchcopy = branchcopy;
}

/** sets destructor method of branching rule */
void SCIPbranchruleSetFree(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHFREE  ((*branchfree))     /**< destructor of branching rule */
   )
{
   assert(branchrule != NULL);

   branchrule->branchfree = branchfree;
}

/** sets initialization method of branching rule */
void SCIPbranchruleSetInit(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit))     /**< initialize branching rule */
   )
{
   assert(branchrule != NULL);

   branchrule->branchinit = branchinit;
}

/** sets deinitialization method of branching rule */
void SCIPbranchruleSetExit(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit))     /**< deinitialize branching rule */
   )
{
   assert(branchrule != NULL);

   branchrule->branchexit = branchexit;
}

/** sets solving process initialization method of branching rule */
void SCIPbranchruleSetInitsol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol)) /**< solving process initialization method of branching rule */
   )
{
   assert(branchrule != NULL);

   branchrule->branchinitsol = branchinitsol;
}

/** sets solving process deinitialization method of branching rule */
void SCIPbranchruleSetExitsol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol)) /**< solving process deinitialization method of branching rule */
   )
{
   assert(branchrule != NULL);

   branchrule->branchexitsol = branchexitsol;
}



/** sets branching execution method for fractional LP solutions */
void SCIPbranchruleSetExecLp(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp))   /**< branching execution method for fractional LP solutions */
   )
{
   assert(branchrule != NULL);

   branchrule->branchexeclp = branchexeclp;
}

/** sets branching execution method for external candidates  */
void SCIPbranchruleSetExecExt(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECEXT((*branchexecext)) /**< branching execution method for external candidates */
   )
{
   assert(branchrule != NULL);

   branchrule->branchexecext = branchexecext;
}

/** sets branching execution method for not completely fixed pseudo solutions */
void SCIPbranchruleSetExecPs(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECPS((*branchexecps))   /**< branching execution method for not completely fixed pseudo solutions */
   )
{
   assert(branchrule != NULL);

   branchrule->branchexecps = branchexecps;
}

/** gets name of branching rule */
const char* SCIPbranchruleGetName(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->name;
}

/** gets description of branching rule */
const char* SCIPbranchruleGetDesc(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->desc;
}

/** gets priority of branching rule */
int SCIPbranchruleGetPriority(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->priority;
}

/** sets priority of branching rule */
void SCIPbranchruleSetPriority(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the branching rule */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);

   branchrule->priority = priority;
   set->branchrulessorted = FALSE;
}

/** gets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
int SCIPbranchruleGetMaxdepth(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->maxdepth;
}

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
void SCIPbranchruleSetMaxdepth(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   int                   maxdepth            /**< new maxdepth of the branching rule */
   )
{
   assert(branchrule != NULL);
   assert(maxdepth >= -1);

   branchrule->maxdepth = maxdepth;
}

/** gets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
SCIP_Real SCIPbranchruleGetMaxbounddist(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->maxbounddist;
}

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
void SCIPbranchruleSetMaxbounddist(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_Real             maxbounddist        /**< new maxbounddist of the branching rule */
   )
{
   assert(branchrule != NULL);
   assert(maxbounddist >= -1);

   branchrule->maxbounddist = maxbounddist;
}

/** enables or disables all clocks of \p branchrule, depending on the value of the flag */
void SCIPbranchruleEnableOrDisableClocks(
   SCIP_BRANCHRULE*      branchrule,         /**< the branching rule for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the branching rule be enabled? */
   )
{
   assert(branchrule != NULL);

   SCIPclockEnableOrDisable(branchrule->setuptime, enable);
   SCIPclockEnableOrDisable(branchrule->branchclock, enable);
}

/** gets time in seconds used in this branching rule for setting up for next stages */
SCIP_Real SCIPbranchruleGetSetupTime(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return SCIPclockGetTime(branchrule->setuptime);
}

/** gets time in seconds used in this branching rule */
SCIP_Real SCIPbranchruleGetTime(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return SCIPclockGetTime(branchrule->branchclock);
}

/** gets the total number of times, the branching rule was called on an LP solution */
SCIP_Longint SCIPbranchruleGetNLPCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->nlpcalls;
}

/** gets the total number of times, the branching rule was called on an external solution */
SCIP_Longint SCIPbranchruleGetNExternCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->nexterncalls;
}

/** gets the total number of times, the branching rule was called on a pseudo solution */
SCIP_Longint SCIPbranchruleGetNPseudoCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->npseudocalls;
}

/** gets the total number of times, the branching rule detected a cutoff */
SCIP_Longint SCIPbranchruleGetNCutoffs(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->ncutoffs;
}

/** gets the total number of cuts, the branching rule separated */
SCIP_Longint SCIPbranchruleGetNCutsFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->ncutsfound;
}

/** gets the total number of constraints, the branching rule added to the respective local nodes (not counting constraints
 *  that were added to the child nodes as branching decisions)
 */
SCIP_Longint SCIPbranchruleGetNConssFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->nconssfound;
}

/** gets the total number of domain reductions, the branching rule found */
SCIP_Longint SCIPbranchruleGetNDomredsFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->ndomredsfound;
}

/** gets the total number of children, the branching rule created */
SCIP_Longint SCIPbranchruleGetNChildren(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->nchildren;
}

/** is branching rule initialized? */
SCIP_Bool SCIPbranchruleIsInitialized(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->initialized;
}




/*
 * branching methods
 */

/** calculates the branching score out of the gain predictions for a binary branching */
SCIP_Real SCIPbranchGetScore(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_Real             downgain,           /**< prediction of objective gain for rounding downwards */
   SCIP_Real             upgain              /**< prediction of objective gain for rounding upwards */
   )
{
   SCIP_Real score;
   SCIP_Real eps;

   assert(set != NULL);

   /* adjust scores near zero to always yield product score greater than 0 */
   eps = SCIPsetSumepsilon(set);
   if( set->branch_sumadjustscore )
   {
      /* adjust scores by adding eps to keep near zero score differences between variables */
      downgain = downgain + eps;
      upgain = upgain + eps;
   }
   else
   {
      /* disregard near zero score differences between variables and consider the branching factor for them */
      downgain = MAX(downgain, eps);
      upgain = MAX(upgain, eps);
   }

   switch( set->branch_scorefunc )
   {
   case 's':  /* linear sum score function */
      /* weigh the two child nodes with branchscorefac and 1-branchscorefac */
      if( downgain > upgain )
         score = set->branch_scorefac * downgain + (1.0-set->branch_scorefac) * upgain;
      else
         score = set->branch_scorefac * upgain + (1.0-set->branch_scorefac) * downgain;
      break;

   case 'p':  /* product score function */
      score = downgain * upgain;
      break;
   case 'q': /* quotient score function */
      if( downgain > upgain )
         score = upgain * upgain / downgain;
      else
         score = downgain * downgain / upgain;
      break;
   default:
      SCIPerrorMessage("invalid branching score function <%c>\n", set->branch_scorefunc);
      SCIPABORT();
      score = 0.0;
   }

   /* apply the branch factor of the variable */
   if( var != NULL )
      score *= SCIPvarGetBranchFactor(var);

   return score;
}

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children */
SCIP_Real SCIPbranchGetScoreMultiple(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int                   nchildren,          /**< number of children that the branching will create */
   SCIP_Real*            gains               /**< prediction of objective gain for each child */
   )
{
   SCIP_Real min1;
   SCIP_Real min2;
   int c;

   assert(nchildren == 0 || gains != NULL);

   /* search for the two minimal gains in the child list and use these to calculate the branching score */
   min1 = SCIPsetInfinity(set);
   min2 = SCIPsetInfinity(set);
   for( c = 0; c < nchildren; ++c )
   {
      if( gains[c] < min1 )
      {
         min2 = min1;
         min1 = gains[c];
      }
      else if( gains[c] < min2 )
         min2 = gains[c];
   }

   return SCIPbranchGetScore(set, var, min1, min2);
}

/** computes a branching point for a (not necessarily discrete) variable
 * a suggested branching point is first projected onto the box
 * if no point is suggested, then the value in the current LP or pseudo solution is used
 * if this value is at infinity, then 0.0 projected onto the bounds and then moved inside the interval is used 
 * for a discrete variable, it is ensured that the returned value is fractional
 * for a continuous variable, the parameter branching/clamp defines how far a branching point need to be from the bounds of a variable
 * the latter is only applied if no point has been suggested, or the suggested point is not inside the variable's interval
 */
SCIP_Real SCIPbranchGetBranchingPoint(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable, of which the branching point should be computed */
   SCIP_Real             suggestion          /**< suggestion for branching point, or SCIP_INVALID if no suggestion */
   )
{
   SCIP_Real branchpoint;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(set != NULL);
   assert(var != NULL);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   /* for a fixed variable, we cannot branch further */
   assert(!SCIPsetIsEQ(set, lb, ub));

   if( !SCIPsetIsInfinity(set, REALABS(suggestion)) )
   {
      /* use user suggested branching point */

      /* first, project it onto the current domain */
      branchpoint = MAX(lb, MIN(suggestion, ub));

      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {
         /* if it is a discrete variable, then round it down and up and accept this choice */
         if( SCIPsetIsEQ(set, branchpoint, ub) )
         {
            /* in the right branch, variable is fixed to its current upper bound */
            return SCIPsetFloor(set, branchpoint) - 0.5;
         }
         else
         {
            /* in the  left branch, variable is fixed to its current lower bound */
            return SCIPsetFloor(set, branchpoint) + 0.5;
         }
      }
      else if( (SCIPsetIsInfinity(set, -lb) || SCIPsetIsRelGT(set, branchpoint, lb)) &&
                (SCIPsetIsInfinity(set,  ub) || SCIPsetIsRelLT(set, branchpoint, ub)) )
      {
         /* if it is continuous and inside the box, then accept it */ 
         return branchpoint;
      }
      /* if it is continuous and suggestion is on of the bounds, continue below */
   }
   else
   {
      /* if no point is suggested and the value in LP solution is not too big, try the LP or pseudo LP solution
       * otherwise, if the value in the LP or pseudosolution is large (here 1e+12), use 0.0
       * in both cases, project onto bounds of course
       */
      branchpoint = SCIPvarGetSol(var, SCIPtreeHasCurrentNodeLP(tree));
      if( REALABS(branchpoint) > 1e+12 )
         branchpoint = 0.0;
      branchpoint = MAX(lb, MIN(branchpoint, ub));
   }

   /* if value is at +/- infty, then choose some value a bit off from bounds or 0.0 */
   if( SCIPsetIsInfinity(set, branchpoint) )
   {
      /* if value is at +infty, then the upper bound should be at infinity */
      assert(SCIPsetIsInfinity(set, ub));

      /* choose 0.0 or something above lower bound if lower bound > 0 */
      if( SCIPsetIsPositive(set, lb) )
         branchpoint = lb + 1000.0;
      else
         branchpoint = 0.0;
   }
   else if( SCIPsetIsInfinity(set, -branchpoint) )
   { 
      /* if value is at -infty, then the lower bound should be at -infinity */
      assert(SCIPsetIsInfinity(set, -lb));

      /* choose 0.0 or something below upper bound if upper bound < 0 */
      if( SCIPsetIsNegative(set, ub) )
         branchpoint = ub - 1000.0;
      else
         branchpoint = 0.0;
   }

   assert(SCIPsetIsInfinity(set,  ub) || SCIPsetIsLE(set, branchpoint, ub));
   assert(SCIPsetIsInfinity(set, -lb) || SCIPsetIsGE(set, branchpoint, lb));

   if( SCIPvarGetType(var) >= SCIP_VARTYPE_IMPLINT )
   {
      if( !SCIPsetIsInfinity(set, -lb) || !SCIPsetIsInfinity(set, ub) )
      {
         /* if one bound is missing, we are temporarily guessing the other one, so we can apply the clamp below */
         if( SCIPsetIsInfinity(set, ub) )
         {
            ub = lb + MIN(MAX(0.5 * REALABS(lb), 1000), 0.9 * (SCIPsetInfinity(set) - lb)); /*lint !e666*/
         }
         else if( SCIPsetIsInfinity(set, -lb) )
         {
            lb = ub - MIN(MAX(0.5 * REALABS(ub), 1000), 0.9 * (SCIPsetInfinity(set) + ub)); /*lint !e666*/
         }

         /* if branching point is too close to the bounds, move more into the middle of the interval */
         if( SCIPrelDiff(ub, lb) <= 2.02 * SCIPsetEpsilon(set) )
         {
            /* for very tiny intervals we set it exactly into the middle */
            branchpoint = (lb+ub)/2.0;
         }
         else
         {
            /* otherwise we project it away from the bounds */
            SCIP_Real minbrpoint;
            SCIP_Real maxbrpoint;
            SCIP_Real scale;
            SCIP_Real lbabs;
            SCIP_Real ubabs;

            lbabs = REALABS(lb);
            ubabs = REALABS(ub);

            scale = MAX3(lbabs, ubabs, 1.0);

            /* the minimal branching point should be
             * - set->clamp away from the lower bound - relative to the local domain size
             * - SCIPsetEpsilon(set)*scale above the lower bound - in absolute value
             */
            minbrpoint = (1.0 - set->branch_clamp) * lb + set->branch_clamp * ub;
            minbrpoint = MAX(lb + 1.01*SCIPsetEpsilon(set)*scale, minbrpoint);  /*lint !e666*/

            /* the maximal branching point should be
             * - set->clamp away from the upper bound - relative to the local domain size
             * - SCIPsetEpsilon(set)*scale below the upper bound - in absolute value
             */
            maxbrpoint = set->branch_clamp * lb + (1.0 - set->branch_clamp) * ub;
            maxbrpoint = MIN(ub - 1.01*SCIPsetEpsilon(set)*scale, maxbrpoint);  /*lint !e666*/

            /* project branchpoint into [minbrpoint, maxbrpoint] */
            branchpoint = MAX(minbrpoint, MIN(branchpoint, maxbrpoint));

            /* if selected branching point is close to 0.0 and bounds are away from 0.0, it often makes sense to branch exactly on 0.0 */
            if( SCIPsetIsFeasZero(set, branchpoint) && SCIPsetIsFeasNegative(set, lb) && SCIPsetIsFeasPositive(set, ub) )
               branchpoint = 0.0;

            assert(SCIPsetIsRelLT(set, SCIPvarGetLbLocal(var), branchpoint));
            assert(SCIPsetIsRelLT(set, branchpoint, SCIPvarGetUbLocal(var)));
         }
      }

      /* ensure fractional branching point for implicit integer variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT && SCIPsetIsIntegral(set, branchpoint) )
      {
         /* if branchpoint is integral but not on bounds, then it should be one of the value {lb+1, ..., ub-1} */
         assert(SCIPsetIsGE(set, SCIPsetRound(set, branchpoint), lb + 1.0));
         assert(SCIPsetIsLE(set, SCIPsetRound(set, branchpoint), ub - 1.0));
         /* if branchpoint is integral, create one branch with x <= x'-1 and one with x >= x'
          * @todo could in the same way be x <= x' and x >= x'+1; is there some easy way to know which is better?
          */
         return branchpoint - 0.5;
      }

      return branchpoint;
   }
   else
   {
      /* discrete variables */
      if( branchpoint <= lb + 0.5 )
      {
         /* if branchpoint is on lower bound, create one branch with x = lb and one with x >= lb+1 */
         return lb + 0.5;
      }
      else if( branchpoint >= ub - 0.5 )
      {
         /* if branchpoint is on upper bound, create one branch with x = ub and one with x <= ub-1 */
         return ub - 0.5;
      }
      else if( SCIPsetIsIntegral(set, branchpoint) )
      {
         /* if branchpoint is integral but not on bounds, then it should be one of the value {lb+1, ..., ub-1} */
         assert(SCIPsetIsGE(set, SCIPsetRound(set, branchpoint), lb + 1.0));
         assert(SCIPsetIsLE(set, SCIPsetRound(set, branchpoint), ub - 1.0));
         /* if branchpoint is integral, create one branch with x <= x'-1 and one with x >= x'
          * @todo could in the same way be x <= x' and x >= x'+1; is there some easy way to know which is better? */
         return branchpoint - 0.5;
      }
      else
      {
         /* branchpoint is somewhere between bounds and fractional, so just round down and up */
         return branchpoint;
      }
   }
}

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
SCIP_RETCODE SCIPbranchExecLP(
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   int i;
   int nalllpcands;  /* sum of binary, integer, and implicit branching candidates */

   assert(branchcand != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* calculate branching candidates */
   SCIP_CALL( branchcandCalcLPCands(branchcand, set, stat, lp) );
   assert(0 <= branchcand->npriolpcands && branchcand->npriolpcands <= branchcand->nlpcands);
   assert((branchcand->npriolpcands == 0) == (branchcand->nlpcands == 0));

   SCIPsetDebugMsg(set, "branching on LP solution with %d (+%d) fractional (+implicit fractional) variables (%d of maximal priority)\n",
      branchcand->nlpcands, branchcand->nimpllpfracs, branchcand->npriolpcands);

   nalllpcands = branchcand->nlpcands + branchcand->nimpllpfracs;
   /* do nothing, if no fractional variables exist */
   if( nalllpcands == 0 )
      return SCIP_OKAY;

   /* if there is a non-fixed variable with higher priority than the maximal priority of the fractional candidates,
    * use pseudo solution branching instead
    */
   if( branchcand->pseudomaxpriority > branchcand->lpmaxpriority )
   {
      SCIP_CALL( SCIPbranchExecPseudo(blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cutoffbound,
            allowaddcons, result) );
      assert(*result != SCIP_DIDNOTRUN && *result != SCIP_DIDNOTFIND);
      return SCIP_OKAY;
   }

   /* sort the branching rules by priority */
   SCIPsetSortBranchrules(set);

   /* try all branching rules until one succeeded to branch */
   for( i = 0; i < set->nbranchrules && (*result == SCIP_DIDNOTRUN || *result == SCIP_DIDNOTFIND) && !SCIPsolveIsStopped(set, stat, FALSE); ++i )
   {
      SCIP_CALL( SCIPbranchruleExecLPSol(set->branchrules[i], set, stat, tree, sepastore, cutoffbound, allowaddcons, result) );
   }

   if( *result == SCIP_DIDNOTRUN || *result == SCIP_DIDNOTFIND )
   {
      SCIP_VAR* var;
      SCIP_Real factor;
      SCIP_Real bestfactor;
      int priority;
      int bestpriority;
      int bestcand;

      /* no branching method succeeded in choosing a branching: just branch on the first fractional variable with maximal
       * priority, and out of these on the one with maximal branch factor
       */
      bestcand = -1;
      bestpriority = INT_MIN;
      bestfactor = SCIP_REAL_MIN;
      for( i = 0; i < nalllpcands; ++i )
      {
         priority = SCIPvarGetBranchPriority(branchcand->lpcands[i]);
         factor = SCIPvarGetBranchFactor(branchcand->lpcands[i]);
         if( priority > bestpriority || (priority == bestpriority && factor > bestfactor) )
         {
            bestcand = i;
            bestpriority = priority;
            bestfactor = factor;
         }
      }
      assert(0 <= bestcand && bestcand < nalllpcands);

      var = branchcand->lpcands[bestcand];
      assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
      assert(branchcand->nlpcands == 0 || SCIPvarGetType(var) != SCIP_VARTYPE_IMPLINT);

      assert(!SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

      SCIP_CALL( SCIPtreeBranchVar(tree, reopt, blkmem, set, stat, transprob, origprob, lp, branchcand, eventqueue, var, SCIP_INVALID,
            NULL, NULL, NULL) );

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** calls branching rules to branch on an external solution; if no external branching candidates exist, the result is SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPbranchExecExtern(
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   int i;

   assert(branchcand != NULL);
   assert(result != NULL);
   assert(0 <= branchcand->nprioexterncands && branchcand->nprioexterncands <= branchcand->nexterncands);
   assert((branchcand->nprioexterncands == 0) == (branchcand->nexterncands == 0));

   *result = SCIP_DIDNOTRUN;

   SCIPsetDebugMsg(set, "branching on external solution with %d branching candidates (%d of maximal priority)\n",
      branchcand->nexterncands, branchcand->nprioexterncands);

   /* do nothing, if no external candidates exist */
   if( branchcand->nexterncands == 0 )
      return SCIP_OKAY;

   /* if there is a non-fixed variable with higher priority than the maximal priority of the external candidates,
    * use pseudo solution branching instead
    */
   if( branchcand->pseudomaxpriority > branchcand->externmaxpriority )
   {
      /* @todo: adjust this, that also LP branching might be called, if lpmaxpriority != externmaxpriority.
       * Therefor, it has to be clear which of both has the higher priority
       */
      SCIP_CALL( SCIPbranchExecPseudo(blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cutoffbound,
            allowaddcons, result) );
      assert(*result != SCIP_DIDNOTRUN && *result != SCIP_DIDNOTFIND);
      return SCIP_OKAY;
   }

   /* sort the branching rules by priority */
   SCIPsetSortBranchrules(set);

   /* try all branching rules until one succeeded to branch */
   for( i = 0; i < set->nbranchrules && (*result == SCIP_DIDNOTRUN || *result == SCIP_DIDNOTFIND); ++i )
   {
      SCIP_CALL( SCIPbranchruleExecExternSol(set->branchrules[i], set, stat, tree, sepastore, cutoffbound, allowaddcons, result) );
   }

   if( *result == SCIP_DIDNOTRUN || *result == SCIP_DIDNOTFIND )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Real bestfactor;
      SCIP_Real bestdomain;
      int bestpriority;
      int bestcand;

      /* if all branching rules did nothing, then they should also not have cleared all branching candidates */
      assert(branchcand->nexterncands > 0);

      /* no branching method succeeded in choosing a branching: just branch on the first branching candidates with maximal
       * priority, and out of these on the one with maximal branch factor, and out of these on the one with largest domain
       */
      bestcand = -1;
      bestpriority = INT_MIN;
      bestfactor = SCIP_REAL_MIN;
      bestdomain = 0.0;
      for( i = 0; i < branchcand->nexterncands; ++i )
      {
         SCIP_VAR* cand;
         SCIP_Real domain;
         SCIP_Real factor;
         int priority;

         cand = branchcand->externcands[i];
         priority = SCIPvarGetBranchPriority(cand);
         factor = SCIPvarGetBranchFactor(cand);

         /* the domain size is infinite, iff one of the bounds is infinite */
         if( SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(cand)) || SCIPsetIsInfinity(set, SCIPvarGetUbLocal(cand)) )
            domain = SCIPsetInfinity(set);
         else
            domain = SCIPvarGetUbLocal(cand) - SCIPvarGetLbLocal(cand);

         /* choose variable with higher priority, higher factor, larger domain (in that order) */
         if( priority > bestpriority || (priority == bestpriority && factor > bestfactor) || (priority == bestpriority && factor == bestfactor && domain > bestdomain) ) /*lint !e777*/
         {
            bestcand = i;
            bestpriority = priority;
            bestfactor = factor;
            bestdomain = domain;
         }
      }
      assert(0 <= bestcand && bestcand < branchcand->nexterncands);

      var = branchcand->externcands[bestcand];
      val = SCIPbranchGetBranchingPoint(set, tree, var, branchcand->externcandssol[bestcand]);
      assert(!SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
      assert(SCIPrelDiff(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var)) <= 2.02 * SCIPsetEpsilon(set)
         || SCIPsetIsLT(set, SCIPvarGetLbLocal(var), val));
      assert(SCIPrelDiff(SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var)) <= 2.02 * SCIPsetEpsilon(set)
         || SCIPsetIsLT(set, val, SCIPvarGetUbLocal(var)));

      SCIPsetDebugMsg(set, "no branching method succeeded; fallback selected to branch on variable <%s> with bounds [%g, %g] on value %g\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), val);

      SCIP_CALL( SCIPtreeBranchVar(tree, reopt, blkmem, set, stat, transprob, origprob, lp, branchcand, eventqueue, var, val,
            NULL, NULL, NULL) );

      if( tree->nchildren >= 1 )
         *result = SCIP_BRANCHED;
      /* if the bounds are too close, it may happen that we cannot branch but rather fix the variable */
      else
      {
         assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
         *result = SCIP_REDUCEDDOM;
      }
   }

   return SCIP_OKAY;
}

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPbranchExecPseudo(
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   int i;

   assert(branchcand != NULL);
   assert(result != NULL);

   SCIPsetDebugMsg(set, "branching on pseudo solution with %d unfixed variables\n", branchcand->npseudocands);

   *result = SCIP_DIDNOTRUN;

   /* do nothing, if no unfixed variables exist */
   if( branchcand->npseudocands == 0 )
      return SCIP_OKAY;

   /* sort the branching rules by priority */
   SCIPsetSortBranchrules(set);

   /* try all branching rules until one succeeded to branch */
   for( i = 0; i < set->nbranchrules && (*result == SCIP_DIDNOTRUN || *result == SCIP_DIDNOTFIND); ++i )
   {
      SCIP_CALL( SCIPbranchruleExecPseudoSol(set->branchrules[i], set, stat, tree, cutoffbound, allowaddcons, result) );
   }

   if( *result == SCIP_DIDNOTRUN || *result == SCIP_DIDNOTFIND )
   {
      SCIP_VAR* var;
      SCIP_Real factor;
      SCIP_Real bestfactor;
      int priority;
      int bestpriority;
      int bestcand;

      /* no branching method succeeded in choosing a branching: just branch on the first unfixed variable with maximal
       * priority, and out of these on the one with maximal branch factor
       */
      bestcand = -1;
      bestpriority = INT_MIN;
      bestfactor = SCIP_REAL_MIN;
      for( i = 0; i < branchcand->npseudocands; ++i )
      {
         priority = SCIPvarGetBranchPriority(branchcand->pseudocands[i]);
         factor = SCIPvarGetBranchFactor(branchcand->pseudocands[i]);
         if( priority > bestpriority || (priority == bestpriority && factor > bestfactor) )
         {
            bestcand = i;
            bestpriority = priority;
            bestfactor = factor;
         }
      }
      assert(0 <= bestcand && bestcand < branchcand->npseudocands);

      var = branchcand->pseudocands[bestcand];
      assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
      assert(!SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

      SCIP_CALL( SCIPtreeBranchVar(tree, reopt, blkmem, set, stat, transprob, origprob, lp, branchcand, eventqueue, var, SCIP_INVALID,
            NULL, NULL, NULL) );

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

