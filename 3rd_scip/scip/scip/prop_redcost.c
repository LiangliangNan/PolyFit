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

/**@file   prop_redcost.c
 * @brief  propagator using the LP reduced cost and the cutoff bound
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Matthias Miltenberger
 * @author Michael Winkler
 *
 * This propagator uses the reduced cost of an optimal solved LP relaxation to propagate the variables against the
 * cutoff bound.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_redcost.h"


/**@name Propagator properties
 *
 * @{
 */

#define PROP_NAME              "redcost"
#define PROP_DESC              "reduced cost strengthening propagator"
#define PROP_TIMING             SCIP_PROPTIMING_DURINGLPLOOP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY          +1000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

/**@} */


/**@name Default parameter values
 *
 * @{
 */

#define DEFAULT_CONTINUOUS        FALSE /**< should reduced cost fixing be also applied to continuous variables? */
#define DEFAULT_USEIMPLICS         TRUE /**< should implications be used to strength the reduced cost for binary variables? */
#define DEFAULT_FORCE             FALSE /**< should the propagator be forced even if active pricer are present? Note that
                                         *   the reductions are always valid, but installing an upper bound on priced
                                         *   variables may lead to problems in pricing (existing variables at their upper
                                         *   bound may be priced again since they may have negative reduced costs) */

/**@} */


/*
 * Data structures
 */


/** propagator data */
struct SCIP_PropData
{
   SCIP_Bool             continuous;         /**< should reduced cost fixing be also applied to continuous variables? */
   SCIP_Real             maxredcost;         /**< maximum reduced cost of a single binary variable */
   SCIP_Bool             usefullimplics;     /**< are the implied reduced cost useful */
   SCIP_Bool             useimplics;         /**< should implications be used to strength the reduced cost for binary variables? */
   SCIP_Bool             force;              /**< should the propagator be forced even if active pricer are present? */
};


/**@name Local methods
 *
 * @{
 */

/** propagate the given binary variable/column using the root reduced cost stored in the SCIP internal data structures
 *  and check if the implications can be useful. Depending on that implications are used or not used during the search to
 *  strength the reduced costs.
 */
static
SCIP_RETCODE propagateRootRedcostBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data structure */
   SCIP_VAR*             var,                /**< variable to use for propagation */
   SCIP_COL*             col,                /**< LP column of the variable */
   SCIP_Real             cutoffbound,        /**< the current cutoff bound */
   int*                  nchgbds             /**< pointer to count the number of bound changes */
   )
{
   SCIP_Real rootredcost;
   SCIP_Real rootsol;
   SCIP_Real rootlpobjval;

   assert(SCIPgetDepth(scip) == 0);

   /* skip binary variable if it is locally fixed */
   if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5 )
      return SCIP_OKAY;

   rootredcost = SCIPvarGetBestRootRedcost(var);
   rootsol = SCIPvarGetBestRootSol(var);
   rootlpobjval = SCIPvarGetBestRootLPObjval(var);

   if( SCIPisDualfeasZero(scip, rootredcost) )
      return SCIP_OKAY;

   assert(rootlpobjval != SCIP_INVALID); /*lint !e777*/

   if( rootsol > 0.5 )
   {
      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasPositive(scip, rootredcost));

      /* update maximum reduced cost of a single binary variable */
      propdata->maxredcost = MAX(propdata->maxredcost, -rootredcost);

      if( rootlpobjval - rootredcost > cutoffbound )
      {
         SCIPdebugMsg(scip, "globally fix binary variable <%s> to 1.0\n", SCIPvarGetName(var));

         SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
         (*nchgbds)++;
         return SCIP_OKAY;
      }
   }
   else
   {
      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasNegative(scip, rootredcost));

      /* update maximum reduced cost of a single binary variable */
      propdata->maxredcost = MAX(propdata->maxredcost, rootredcost);

      if( rootlpobjval + rootredcost > cutoffbound )
      {
         SCIPdebugMsg(scip, "globally fix binary variable <%s> to 0.0\n", SCIPvarGetName(var));

         SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
         (*nchgbds)++;
         return SCIP_OKAY;
      }
   }

   /* evaluate if the implications are useful; the implications are seen to be useful if they provide an increase for
    * the root reduced costs
    */
   if( !propdata->usefullimplics )
   {
      SCIP_Real lbredcost;
      SCIP_Real ubredcost;

      lbredcost = SCIPgetVarImplRedcost(scip, var, FALSE);
      assert(!SCIPisDualfeasPositive(scip, lbredcost));

      ubredcost = SCIPgetVarImplRedcost(scip, var, TRUE);
      assert(!SCIPisDualfeasNegative(scip, ubredcost));

      switch( SCIPcolGetBasisStatus(col) )
      {
      case SCIP_BASESTAT_LOWER:
         ubredcost -= SCIPgetVarRedcost(scip, var);

         /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
          * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
          * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
          * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
          */
         assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasNegative(scip, ubredcost));
         break;

      case SCIP_BASESTAT_UPPER:
         lbredcost -= SCIPgetVarRedcost(scip, var);

         /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
          * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
          * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
          * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
          */
         assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasPositive(scip, lbredcost));
         break;

      case SCIP_BASESTAT_BASIC:
      case SCIP_BASESTAT_ZERO:
      default:
         break;
      }

      propdata->usefullimplics = (lbredcost < 0.0) || (ubredcost > 0.0);
   }

   return SCIP_OKAY;
}

/** propagate the given binary variable/column using the reduced cost */
static
SCIP_RETCODE propagateRedcostBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data structure */
   SCIP_VAR*             var,                /**< variable to use for propagation */
   SCIP_COL*             col,                /**< LP column of the variable */
   SCIP_Real             requiredredcost,    /**< required reduset cost to be able to fix a binary variable */
   int*                  nchgbds,            /**< pointer to count the number of bound changes */
   SCIP_Bool*            cutoff              /**< pointer to store if an cutoff was detected */
   )
{
   SCIP_Real lbredcost;
   SCIP_Real ubredcost;
   SCIP_Real redcost;

   /* skip binary variable if it is locally fixed */
   if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5 )
      return SCIP_OKAY;

   /* first use the redcost cost to fix the binary variable */
   switch( SCIPcolGetBasisStatus(col) )
   {
   case SCIP_BASESTAT_LOWER:
      redcost = SCIPgetVarRedcost(scip, var);

      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasNegative(scip, redcost));

      if( redcost > requiredredcost )
      {
         SCIPdebugMsg(scip, "variable <%s>: fixed 0.0 (requiredredcost <%g>, redcost <%g>)\n",
            SCIPvarGetName(var), requiredredcost, redcost);

         SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
         (*nchgbds)++;
         return SCIP_OKAY;
      }
      break;

   case SCIP_BASESTAT_UPPER:
      redcost = SCIPgetVarRedcost(scip, var);

      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasPositive(scip, redcost));

      if( -redcost > requiredredcost )
      {
         SCIPdebugMsg(scip, "variable <%s>: fixed 1.0 (requiredredcost <%g>, redcost <%g>)\n",
            SCIPvarGetName(var), requiredredcost, redcost);

         SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
         (*nchgbds)++;
         return SCIP_OKAY;
      }
      break;

   case SCIP_BASESTAT_BASIC:
      return SCIP_OKAY;

   case SCIP_BASESTAT_ZERO:
      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || SCIPisDualfeasZero(scip, SCIPgetColRedcost(scip, col)));
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid basis state\n");
      return SCIP_INVALIDDATA;
   }

   /* second, if the implications should be used and if the implications are seen to be promising use the implied
    * reduced costs to fix the binary variable
    */
   if( propdata->useimplics && propdata->usefullimplics )
   {
      /* collect implied reduced costs if the variable would be fixed to its lower bound */
      lbredcost = SCIPgetVarImplRedcost(scip, var, FALSE);

      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasPositive(scip, lbredcost)
            || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );

      /* collect implied reduced costs if the variable would be fixed to its upper bound */
      ubredcost = SCIPgetVarImplRedcost(scip, var, TRUE);
      assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasNegative(scip, ubredcost)
            || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );

      if( -lbredcost > requiredredcost && ubredcost > requiredredcost )
      {
         SCIPdebugMsg(scip, "variable <%s>: cutoff (requiredredcost <%g>, lbredcost <%g>, ubredcost <%g>)\n",
            SCIPvarGetName(var), requiredredcost, lbredcost, ubredcost);

         (*cutoff) = TRUE;
      }
      else if( -lbredcost > requiredredcost )
      {
         SCIPdebugMsg(scip, "variable <%s>: fixed 1.0 (requiredredcost <%g>, redcost <%g>, lbredcost <%g>)\n",
            SCIPvarGetName(var), requiredredcost, redcost, lbredcost);

         SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
         (*nchgbds)++;
      }
      else if( ubredcost > requiredredcost )
      {
         SCIPdebugMsg(scip, "variable <%s>: fixed 0.0 (requiredredcost <%g>, redcost <%g>, ubredcost <%g>)\n",
            SCIPvarGetName(var), requiredredcost, redcost, ubredcost);

         SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
         (*nchgbds)++;
      }

      /* update maximum reduced cost of a single binary variable */
      propdata->maxredcost = MAX3(propdata->maxredcost, -lbredcost, ubredcost);
   }

   return SCIP_OKAY;
}

/** propagate the given none binary variable/column using the reduced cost */
static
SCIP_RETCODE propagateRedcostVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to use for propagation */
   SCIP_COL*             col,                /**< LP column of the variable */
   SCIP_Real             lpobjval,           /**< objective value of the current LP */
   SCIP_Real             cutoffbound,        /**< the current cutoff bound */
   int*                  nchgbds             /**< pointer to count the number of bound changes */
   )
{
   SCIP_Real redcost;

   switch( SCIPcolGetBasisStatus(col) )
   {
   case SCIP_BASESTAT_LOWER:
      redcost = SCIPgetColRedcost(scip, col);

      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasNegative(scip, redcost)
            || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );

      if( SCIPisDualfeasPositive(scip, redcost) )
      {
         SCIP_Real oldlb;
         SCIP_Real oldub;

         oldlb = SCIPvarGetLbLocal(var);
         oldub = SCIPvarGetUbLocal(var);
         assert(SCIPisEQ(scip, oldlb, SCIPcolGetLb(col)));
         assert(SCIPisEQ(scip, oldub, SCIPcolGetUb(col)));

         if( SCIPisFeasLT(scip, oldlb, oldub) )
         {
            SCIP_Real newub;
            SCIP_Bool strengthen;

            /* calculate reduced cost based bound */
            newub = (cutoffbound - lpobjval) / redcost + oldlb;

            /* check, if new bound is good enough:
             *  - integer variables: take all possible strengthenings
             *  - continuous variables: strengthening must cut part of the variable's dynamic range, and
             *                          at least 20% of the current domain
             */
            if( SCIPvarIsIntegral(var) )
            {
               newub = SCIPadjustedVarUb(scip, var, newub);
               strengthen = (newub < oldub - 0.5);
            }
            else
               strengthen = (newub < SCIPcolGetMaxPrimsol(col) && newub <= 0.2 * oldlb + 0.8 * oldub);

            if( strengthen )
            {
               /* strengthen upper bound */
               SCIPdebugMsg(scip, "redcost strengthening upper bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                  SCIPvarGetName(var), oldlb, oldub, oldlb, newub, cutoffbound, lpobjval, redcost);
               SCIP_CALL( SCIPchgVarUb(scip, var, newub) );
               (*nchgbds)++;
            }
         }
      }
      break;

   case SCIP_BASESTAT_BASIC:
      break;

   case SCIP_BASESTAT_UPPER:
      redcost = SCIPgetColRedcost(scip, col);

      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasPositive(scip, redcost)
            || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );

      if( SCIPisDualfeasNegative(scip, redcost) )
      {
         SCIP_Real oldlb;
         SCIP_Real oldub;

         oldlb = SCIPvarGetLbLocal(var);
         oldub = SCIPvarGetUbLocal(var);
         assert(SCIPisEQ(scip, oldlb, SCIPcolGetLb(col)));
         assert(SCIPisEQ(scip, oldub, SCIPcolGetUb(col)));

         if( SCIPisFeasLT(scip, oldlb, oldub) )
         {
            SCIP_Real newlb;
            SCIP_Bool strengthen;

            /* calculate reduced cost based bound */
            newlb = (cutoffbound - lpobjval) / redcost + oldub;

            /* check, if new bound is good enough:
             *  - integer variables: take all possible strengthenings
             *  - continuous variables: strengthening must cut part of the variable's dynamic range, and
             *                          at least 20% of the current domain
             */
            if( SCIPvarIsIntegral(var) )
            {
               newlb = SCIPadjustedVarLb(scip, var, newlb);
               strengthen = (newlb > oldlb + 0.5);
            }
            else
               strengthen = (newlb > SCIPcolGetMinPrimsol(col) && newlb >= 0.8 * oldlb + 0.2 * oldub);

            /* check, if new bound is good enough: at least 20% strengthening for continuous variables */
            if( strengthen )
            {
               /* strengthen lower bound */
               SCIPdebugMsg(scip, "redcost strengthening lower bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                  SCIPvarGetName(var), oldlb, oldub, newlb, oldub, cutoffbound, lpobjval, redcost);
               SCIP_CALL( SCIPchgVarLb(scip, var, newlb) );
               (*nchgbds)++;
            }
         }
      }
      break;

   case SCIP_BASESTAT_ZERO:
      /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
       * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
       * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
       * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
       */
      assert(!SCIPisLPDualReliable(scip) || SCIPisDualfeasZero(scip, SCIPgetColRedcost(scip, col)));
      break;

   default:
      SCIPerrorMessage("invalid basis state\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of propagator
 *
 * @{
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyRedcost)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludePropRedcost(scip) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
/**! [SnippetPropFreeRedcost] */
static
SCIP_DECL_PROPFREE(propFreeRedcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeBlockMemory(scip, &propdata);

   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}
/**! [SnippetPropFreeRedcost] */

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolRedcost)
{
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->usefullimplics = FALSE;
   propdata->maxredcost = 0.0;

   return SCIP_OKAY;
}

/** reduced cost propagation method for an LP solution */
static
SCIP_DECL_PROPEXEC(propExecRedcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_COL** cols;
   SCIP_Real requiredredcost;
   SCIP_Real cutoffbound;
   SCIP_Real lpobjval;
   SCIP_Bool propbinvars;
   SCIP_Bool cutoff;
   int nchgbds;
   int ncols;
   int c;

   *result = SCIP_DIDNOTRUN;

   /* in case we have a zero objective function, we skip the reduced cost propagator */
   if( SCIPgetNObjVars(scip) == 0 )
      return SCIP_OKAY;

   /* propagator can only be applied during solving stage */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* we cannot apply reduced cost fixing, if we want to solve exactly */
   /**@todo implement reduced cost fixing with interval arithmetics */
   if( SCIPisExactSolve(scip) )
      return SCIP_OKAY;

   /* only call propagator, if the current node has an LP */
   if( !SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   /* only call propagator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call propagator, if the current LP is a valid relaxation */
   if( !SCIPisLPRelax(scip) )
      return SCIP_OKAY;

   /* we cannot apply reduced cost strengthening, if no simplex basis is available */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* do not run if propagation w.r.t. objective is not allowed */
   if( !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   /* get current cutoff bound */
   cutoffbound = SCIPgetCutoffbound(scip);

   /* reduced cost strengthening can only be applied, if we have a finite cutoff */
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   /* get LP columns */
   cols = SCIPgetLPCols(scip);
   ncols = SCIPgetNLPCols(scip);

   /* do nothing if the LP has no columns (is empty) */
   if( ncols == 0 )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* check if all integral variables are fixed and the continuous variables should not be propagated */
   if( !propdata->continuous && SCIPgetNPseudoBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* get LP objective value */
   lpobjval = SCIPgetLPObjval(scip);

   /* check if binary variables should be propagated */
   propbinvars = (SCIPgetDepth(scip) == 0) || (cutoffbound - lpobjval < 5 * propdata->maxredcost);

   /* skip the propagator if the problem has only binary variables and those should not be propagated */
   if( !propbinvars && SCIPgetNVars(scip) == SCIPgetNBinVars(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;
   cutoff = FALSE;
   nchgbds = 0;

   /* compute the required reduced cost which are needed for a binary variable to be fixed */
   requiredredcost = cutoffbound - lpobjval;

   SCIPdebugMsg(scip, "lpobjval <%g>, cutoffbound <%g>, max reduced <%g>, propgate binary %u, use implics %u\n",
      lpobjval, cutoffbound, propdata->maxredcost, propbinvars, propdata->usefullimplics);

   /* check reduced costs for non-basic columns */
   for( c = 0; c < ncols && !cutoff; ++c )
   {
      SCIP_VAR* var;

      var = SCIPcolGetVar(cols[c]);

      /* skip continuous variables in case the corresponding parameter is set */
      if( !propdata->continuous && !SCIPvarIsIntegral(var) )
         continue;

      if( SCIPvarIsBinary(var) )
      {
         if( propbinvars )
         {
            if( SCIPgetDepth(scip) == 0 )
            {
               SCIP_CALL( propagateRootRedcostBinvar(scip, propdata, var, cols[c], cutoffbound, &nchgbds) );
            }
            else
            {
               SCIP_CALL( propagateRedcostBinvar(scip, propdata, var, cols[c], requiredredcost, &nchgbds, &cutoff) );
            }
         }
      }
      else
      {
         SCIP_CALL( propagateRedcostVar(scip, var, cols[c], lpobjval, cutoffbound, &nchgbds) );
      }
   }

   if( cutoff )
   {
      *result = SCIP_CUTOFF;

      SCIPdebugMsg(scip, "node %" SCIP_LONGINT_FORMAT ": detected cutoff\n",
         SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
   }
   else if( nchgbds > 0 )
   {
      *result = SCIP_REDUCEDDOM;

      SCIPdebugMsg(scip, "node %" SCIP_LONGINT_FORMAT ": %d bound changes (max redcost <%g>)\n",
         SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) , nchgbds, propdata->maxredcost);
   }

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the redcost propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropRedcost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create redcost propagator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecRedcost, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyRedcost) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolRedcost) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeRedcost) );

   /* add redcost propagator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/continuous",
         "should reduced cost fixing be also applied to continuous variables?",
         &propdata->continuous, FALSE, DEFAULT_CONTINUOUS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/useimplics",
         "should implications be used to strength the reduced cost for binary variables?",
         &propdata->useimplics, FALSE, DEFAULT_USEIMPLICS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/force",
         "should the propagator be forced even if active pricer are present?",
         &propdata->force, TRUE, DEFAULT_FORCE, NULL, NULL) );

   return SCIP_OKAY;
}

/**@} */
