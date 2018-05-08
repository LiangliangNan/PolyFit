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

/**@file   prop_rootredcost.c
 * @brief  reduced cost strengthening using root node reduced costs and the cutoff bound
 * @author Tobias Achterberg
 * @author Stefan Heinz
 *
 * This propagator uses the root reduced cost to (globally) propagate against the cutoff bound. The propagator checks if
 * the variables with non-zero root reduced cost can exceed the cutoff bound. If this is the case the corresponding
 * bound can be tightened.
 *
 * The propagate is performed during the search any time a new cutoff bound (primal solution) is found.
 *
 * @todo do not sort the variables; just store the cutoff bound which leads to a fixing. If that appears loop over all
 *       variables and fix and store the next cutoff bound which leads to an fixing
 * @todo resolve the root LP in case of repropagation and update root reduced costs use root LP counter to check if new
 *       best root combinations might be available
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_rootredcost.h"

/**@name Propagator properties
 *
 * @{
 */

#define PROP_NAME              "rootredcost"
#define PROP_DESC              "reduced cost strengthening using root node reduced costs and the cutoff bound"
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY         +10000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

/**@} */

/**@name Default parameter values
 *
 * @{
 */
#define DEFAULT_ONLYBINARY        FALSE /**< should only binary variables be propagated? */
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
   SCIP_VAR**            redcostvars;        /**< variables with non-zero root reduced cost */
   SCIP_Real             lastcutoffbound;    /**< cutoff bound for which the root reduced costs were already processed */
   int                   nredcostvars;       /**< number of variables with non-zero root reduced cost */
   int                   nredcostbinvars;    /**< number of binary variables with non-zero root reduced cost */
   int                   glbfirstnonfixed;   /**< index of first globally non-fixed binary variable */
   SCIP_Bool             initialized;        /**< is the propagator data initialized */
   SCIP_Bool             onlybinary;         /**< should only binary variables be propagated? */
   SCIP_Bool             force;              /**< should the propagator be forced even if active pricer are present? */
};


/**@name Local methods
 *
 * @{
 */

/** reset structure memember of propagator data structure */
static
void propdataReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data to reset */
   )
{
   propdata->redcostvars = NULL;
   propdata->lastcutoffbound = SCIP_INVALID;
   propdata->nredcostvars = 0;
   propdata->nredcostbinvars = 0;
   propdata->glbfirstnonfixed = 0;
   propdata->initialized = FALSE;
}

/** compare variables with non-zero reduced cost w.r.t.
 *  (i)  the cutoff bound which would lead to a fixing
 *  (ii) variable problem index;
 */
static
SCIP_DECL_SORTPTRCOMP(varCompRedcost)
{
   SCIP_VAR* var1 = (SCIP_VAR*)elem1;
   SCIP_VAR* var2 = (SCIP_VAR*)elem2;
   SCIP_Real key1;
   SCIP_Real key2;

   assert(SCIPvarIsBinary(var1));
   assert(SCIPvarGetBestRootRedcost(var1) != 0.0);

   assert(SCIPvarIsBinary(var2));
   assert(SCIPvarGetBestRootRedcost(var2) != 0.0);

   /* collect sorting key for both variables */
   key1 = REALABS(SCIPvarGetBestRootRedcost(var1)) + SCIPvarGetBestRootLPObjval(var1);
   key2 = REALABS(SCIPvarGetBestRootRedcost(var2)) + SCIPvarGetBestRootLPObjval(var2);

   if( key1 < key2 )
      return -1;
   else if( key1 > key2 )
      return +1;

   /* second criteria use the problem index
    *
    * @note The problem index is unique. That means the resulting sorting is unique.
    */
   return SCIPvarCompare(var1, var2);
}

/** create propagator data structure */
static
SCIP_RETCODE propdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< pointer to store the created propagator data */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, propdata) );

   propdataReset(scip, *propdata);

   return SCIP_OKAY;
}

/** counts the number of variables with non-zero root reduced cost */
static
int countNonZeroRootRedcostVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< number of variables */
   )
{
   int count;
   int v;

   count = 0;

   /* count number of variables with non-zero root reduced cost */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real redcost;

      assert(vars[v] != NULL);

      redcost = SCIPvarGetBestRootRedcost(vars[v]);
      if( !SCIPisDualfeasZero(scip, redcost) )
         count++;
   }

   return count;
}

/** free propagator data */
static
SCIP_RETCODE propdataExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int v;

   /* release all variables */
   for( v = 0; v < propdata->nredcostvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &propdata->redcostvars[v]) );
   }

   /* free memory for non-zero reduced cost variables */
   SCIPfreeBlockMemoryArrayNull(scip, &propdata->redcostvars, propdata->nredcostvars);

   propdataReset(scip, propdata);

   return SCIP_OKAY;
}

/** initializate the propagator */
static
SCIP_RETCODE propdataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nredcostvars;
   int nredcostbinvars;
   int v;

   assert(scip != NULL);
   assert(propdata != NULL);

   /* check if the propagator data structure is already initialized */
   if( propdata->initialized )
      return SCIP_OKAY;

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nbinvars = SCIPgetNBinVars(scip);

   /* count binary variables with non-zero root reduced cost */
   nredcostbinvars = countNonZeroRootRedcostVars(scip, vars, nbinvars);
   SCIPdebugMsg(scip, "There are %d (poor) binary variables with non-zero root reduced cost <%s>.\n", nredcostbinvars, SCIPgetProbName(scip));

   /* count non-binary variables with non-zero root reduced cost */
   nredcostvars = countNonZeroRootRedcostVars(scip, &vars[nbinvars], nvars - nbinvars);

   nredcostvars += nredcostbinvars;

   /* collect the variables with non-zero reduced costs */
   if( nredcostvars > 0 )
   {
      int k;

      k = 0;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->redcostvars, nredcostvars) );

      SCIPdebugMsg(scip, "Store non-zero root reduced cost variables at address <%p>.\n", (void*)propdata->redcostvars);

      for( v = 0; v < nvars; ++v )
      {
         SCIP_Real redcost;
         SCIP_VAR* var;

         var = vars[v];
         redcost = SCIPvarGetBestRootRedcost(var);

         if( SCIPisDualfeasZero(scip, redcost) )
            continue;

         assert(k < nredcostvars);

         /* check if one of the non-binary variables is implicit binary */
         if( k >= nredcostbinvars && SCIPvarIsBinary(var) )
         {
            /* move the first non-binary variable to end of the array */
            propdata->redcostvars[k] = propdata->redcostvars[nredcostbinvars];

            /* place the binary variable at the end of the binary section */
            propdata->redcostvars[nredcostbinvars] = var;
            nredcostbinvars++;
         }
         else
            propdata->redcostvars[k] = var;

         /* captures the variable */
         SCIP_CALL( SCIPcaptureVar(scip, var) ) ;

         k++;

         /* check if already visited all variable with non-zero redcostective coefficient */
         if( k == nredcostvars )
            break;
      }

      /* sort binary variables with respect to their cutoff bound which would lead to an fixing; this order can be used
       * to efficiently propagate the binary variables
       */
      SCIPsortDownPtr((void**)propdata->redcostvars, varCompRedcost, nredcostbinvars);

      assert(k == nredcostvars);

      SCIPdebugMsg(scip, "variables with non-zero redcostective coefficient: %d binaries, %d non-binaries\n", nredcostbinvars, nredcostvars - nredcostbinvars);
   }

   propdata->nredcostvars = nredcostvars;
   propdata->nredcostbinvars = nredcostbinvars;
   propdata->glbfirstnonfixed = 0;
   propdata->lastcutoffbound = SCIPinfinity(scip);
   propdata->initialized = TRUE;

   return SCIP_OKAY;
}

/** propagates the root reduced cost against the cutoff bound for the given variable */
static
SCIP_RETCODE propagateRootRedcostVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to propagate */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store if the bound was tightened */
   )
{
   SCIP_Real rootsol;
   SCIP_Real rootredcost;
   SCIP_Real rootlpobjval;
   SCIP_Real newbd;

   rootredcost = SCIPvarGetBestRootRedcost(var);
   assert(rootredcost != SCIP_INVALID); /*lint !e777*/

   /* SCIPisLPDualReliable should always return TRUE if the dual feasibility check is enabled and the LP claims to
    * have a dual feasible solution. if the check is disabled the dual solution might be incorrect and the assert
    * might fail. however, if the user decides to disable the dual feasibility check (which also can lead to wrong
    * cutoffs) we don't want to skip propagating with reduced costs as an unexpected side-effect.
    */
   assert(!SCIPisLPDualReliable(scip) || !SCIPisDualfeasZero(scip, rootredcost));

   rootsol = SCIPvarGetBestRootSol(var);
   rootlpobjval = SCIPvarGetBestRootLPObjval(var);

   /* calculate reduced cost based bound */
   newbd = rootsol + (cutoffbound - rootlpobjval) / rootredcost;

   if( SCIPisDualfeasPositive(scip, rootredcost) )
   {
      assert(SCIPisFeasLE(scip, rootsol, SCIPvarGetLbGlobal(var))); /* lb might have been increased in the meantime */

      /* strengthen upper bound */
      SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newbd, FALSE, infeasible, tightened) );
   }
   else
   {
      assert(!SCIPisLPDualReliable(scip) || SCIPisDualfeasNegative(scip, rootredcost));
      assert(SCIPisFeasGE(scip, rootsol, SCIPvarGetUbGlobal(var))); /* ub might have been decreased in the meantime */

      /* strengthen lower bound */
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newbd, FALSE, infeasible, tightened) );
   }

   return SCIP_OKAY;
}

/** propagate binary variables with non-zero root reduced cost */
static
SCIP_RETCODE propagateBinaryBestRootRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data structure */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   int*                  nchgbds,            /**< pointer to store the number of bound changes */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_VAR** redcostvars;
   int v;

   assert(!(*cutoff));

   /* the binary variables are stored in the beginning of the variable array; these variables are sorted w.r.t. cutoff
    * bound which would lead to a fixing; that give us an abort criteria (see below)
    */
   redcostvars = propdata->redcostvars;
   assert(redcostvars != NULL);

#ifndef NDEBUG
   /* check that the binary variables are correctly sorted
    *
    * @note In case the assertion fails it indicates that a new LP root solving arose after we initialized the
    *       propagator; The new LP solution destroyed the sorting of the binary variables since the reduced cost of the
    *       variables changed. This could lead to potentially missing a domain reductions. Currently, it is not possible to
    *       check if a new root LP was solved, changing the root reduced costs. This case, however, should not happen in the
    *       current SCIP version.
    */
   for( v = 1; v < propdata->nredcostbinvars; ++v )
      assert(varCompRedcost(redcostvars[v-1], redcostvars[v]) == 1);

   /* check that the variables before glbfirstnonfixed are globally fixed */
   for( v = 0; v < propdata->glbfirstnonfixed; ++v )
   {
      SCIP_VAR* var;

      var =  redcostvars[v];
      assert(var != NULL);

      assert(SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5);
   }
#endif

   /* propagate binary variables */
   for( v = propdata->glbfirstnonfixed; v < propdata->nredcostbinvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool tightened;

      var =  redcostvars[v];
      assert(var != NULL);

      /* check if the variables is already globally fixed; if so continue with the next potential candidate */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
         continue;

      /* try to tighten the variable bound */
      SCIP_CALL( propagateRootRedcostVar(scip, var, cutoffbound, cutoff, &tightened) );

      if( tightened )
      {
         /* @note The variable might not be globally fixed right away since this would destroy the local internal data
          *       structure of a search node; the bound change is in that case pending; hence we cannot assert that the
          *       variable is globally fixed
          */
         assert(!(*cutoff));

         SCIPdebugMsg(scip, "globally fixed binary variable <%s> [%g,%g] bestroot sol <%g>, redcost <%g>, lpobj <%g>\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
            SCIPvarGetBestRootSol(var), SCIPvarGetBestRootRedcost(var), SCIPvarGetBestRootLPObjval(var) );

         (*nchgbds)++;
      }
      else
      {
         assert(!tightened);

         /* The binary variables are sorted in non-increasing manner w.r.t. their cutoff bound which would lead to a
          * global fixing; That is, abs(rootredcost) + rootlpobjval. Depending on the sign of the reduced cost the
          * following two cases can arise for binary variables which are not fixed globally yet:
          *
          * - redcost > 0 --> newub = 0.0 + (cutoffbound - lpobjval) / redcost --> newub = 0 <=> cutoffbound < redcost + lpobjval = sorting key
          * - redcost < 0 --> newlb = 1.0 + (cutoffbound - lpobjval) / redcost --> newlb = 1 <=> cutoffbound < -redcost + lpobjval = sorting key
          *
          * Due to the order of the binary variables it follows if one binary variable cannot be fixed anymore all the
          * remaining once can also not be fixed since these have only an smaller or equal cutoff bound which would lead
          * to global fixing. Hence, we can break that loop.
          *
          * Note that variables with non-zero reduced cost are sitting at one of their bound. That is the lower one if
          * the reduced cost are positive and the upper bound if the reduced cost are negative. For binary variables
          * that is 0 for the lower bound and 1 for the upper bound.
          */
         SCIPdebugMsg(scip, "interrupt propagation for binary variables after %d from %d binary variables\n",
            v, propdata->nredcostbinvars);

         if( *cutoff )
         {
            SCIPdebugMsg(scip, "detected cutoff: binary variable <%s> [%g,%g], redcost <%g>, rootsol <%g>, rootlpobjval <%g>\n",
               SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
               SCIPvarGetBestRootRedcost(var), SCIPvarGetBestRootSol(var), SCIPvarGetBestRootLPObjval(var));
         }

         break;
      }
   }
   /* store the index of the variable which is not globally fixed */
   propdata->glbfirstnonfixed = v;

#if 0 /* due to numerics it might be that the abort criteria did not work correctly, because the sorting mechanism may
       * have evaluated variables with a really small difference in their reduced cost values but with really huge
       * lpobjval as the same
       */
#ifndef NDEBUG
   /* check that the abort criteria works; that means none of the remaining binary variables can be fixed */
   for( ; v < propdata->nredcostbinvars && !(*cutoff); ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool tightened;

      var =  redcostvars[v];
      assert(var != NULL);

      /* check if the variables is already globally fixed; if so continue with the potential candidate */
      if( SCIPvarGetLbGlobal(var) > 0.5 || SCIPvarGetUbGlobal(var) < 0.5)
         continue;

      /* try to tighten the variable bound */
      SCIP_CALL( propagateRootRedcostVar(scip, var, cutoffbound, cutoff, &tightened) );
      assert(!tightened);
      assert(!(*cutoff));
   }
#endif
#endif

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods of propagator
 *
 * @{
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyRootredcost)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropRootredcost(scip) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeRootredcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->redcostvars == NULL);

   SCIPfreeBlockMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolRootredcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* reset propagator data structure */
   SCIP_CALL( propdataExit(scip, propdata) );

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecRootredcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** redcostvars;
   SCIP_Real cutoffbound;
   SCIP_Real lpobjval;
   SCIP_Bool cutoff;
   int nredcostvars;
   int nchgbds;
   int v;

   *result = SCIP_DIDNOTRUN;

   /* in case we have a zero objective fucntion, we skip the root reduced cost propagator */
   if( SCIPgetNObjVars(scip) == 0 )
      return SCIP_OKAY;

   /* propagator can only be applied during solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* the propagator should run in all nodes except the root node; for the root node the redcost propagator does
    * the job already
    */
   if( SCIPgetDepth(scip) < 1 )
      return SCIP_OKAY;

   /* first check root LP objective value if it exists */
   lpobjval = SCIPgetLPRootObjval(scip);
   if( lpobjval == SCIP_INVALID ) /*lint !e777*/
      return SCIP_OKAY;

   /* do not run in probing mode since this propagator changes global variable bounds */
   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* do not run if propagation w.r.t. objective is not allowed */
   if( !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* do nothing if active pricer are present and force flag is not TRUE */
   if( !propdata->force && SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* get current cutoff bound */
   cutoffbound = SCIPgetCutoffbound(scip);

   /* reduced cost strengthening can only be applied, if we have a finite upper bound on the LP value */
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   /* initialize propagator data structure */
   SCIP_CALL( propdataInit(scip, propdata) );
   assert(cutoffbound <= propdata->lastcutoffbound);

   if( cutoffbound == propdata->lastcutoffbound ) /*lint !e777*/
      return SCIP_OKAY;

   /* get variables */
   redcostvars = propdata->redcostvars;
   nredcostvars = propdata->nredcostvars;

   /* since no variables has non-zero reduced cost do nothing */
   if( nredcostvars == 0 )
      return SCIP_OKAY;

   /* store cutoff bound to remember later that for that particular cutoff bound the propagation was already
    * preformed
    */
   propdata->lastcutoffbound = cutoffbound;

   SCIPdebugMsg(scip, "searching for root reduced cost fixings\n");
   SCIPdebugMsg(scip, "-> cutoffbound <%g>\n", cutoffbound);
   SCIPdebugMsg(scip, "-> LP objective value <%g>\n", lpobjval);

   *result = SCIP_DIDNOTFIND;
   nchgbds = 0;
   cutoff = FALSE;

   /* propagate the binary variables with non-zero root reduced cost */
   SCIP_CALL( propagateBinaryBestRootRedcost(scip, propdata, cutoffbound, &nchgbds, &cutoff) );

   if( !propdata->onlybinary )
   {
      /* check reduced costs for non-binary variables that were columns of the root LP */
      for( v = propdata->nredcostbinvars; v < nredcostvars && !cutoff; ++v )
      {
         SCIP_VAR* var;
         SCIP_Bool tightened;

         var = redcostvars[v];
         assert(var != NULL);

         /* try to tighten the variable bound */
         SCIP_CALL( propagateRootRedcostVar(scip, var, cutoffbound, &cutoff, &tightened) );

         if( tightened )
            nchgbds++;
      }
   }

   /* evaluate propagation results */
   if( cutoff )
   {
      /* we are done with solving since the cutoff is w.r.t. a global bound change; cutoff root node */
      SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
      (*result) = SCIP_CUTOFF;
   }
   else if( nchgbds > 0 )
      (*result) = SCIP_REDUCEDDOM;

   SCIPdebugMsg(scip, "tightened %d variable domains (%u cutoff)\n", nchgbds, cutoff);

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the root node reduced cost strengthening propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropRootredcost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create rootredcost propagator data */
   SCIP_CALL( propdataCreate(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecRootredcost, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyRootredcost) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeRootredcost) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolRootredcost) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/onlybinary",
         "should only binary variables be propagated?",
         &propdata->onlybinary, TRUE, DEFAULT_ONLYBINARY, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/force",
         "should the propagator be forced even if active pricer are present?",
         &propdata->force, TRUE, DEFAULT_FORCE, NULL, NULL) );

   return SCIP_OKAY;
}

/**@} */
