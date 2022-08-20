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

/**@file   prop_dualfix.c
 * @brief  fixing roundable variables to best bound
 * @author Tobias Achterberg
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_dualfix.h"


#define PROP_NAME                  "dualfix"
#define PROP_DESC                  "roundable variables dual fixing"
#define PROP_TIMING                 SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY               +8000000 /**< propagation priority */
#define PROP_FREQ                          0 /**< propagation frequency */
#define PROP_DELAY                     FALSE /**< should propagation method be delayed, if other propagators found
                                             *   reductions? */
#define PROP_PRESOL_PRIORITY        +8000000 /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
#define PROP_PRESOL_MAXROUNDS             -1 /**< maximal number of propving rounds the propver participates in (-1: no limit) */
#define PROP_PRESOLTIMING           SCIP_PRESOLTIMING_FAST /* timing of the presolving method (fast, medium, or exhaustive) */


/**@name Local methods
 *
 * @{
 */

/** perform dual presolving */
static
SCIP_RETCODE performDualfix(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nfixedvars,         /**< pointer to store number of fixed variables */
   SCIP_Bool*            unbounded,          /**< pointer to store if an unboundness was detected */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   /* get active problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* look for fixable variables
    * loop backwards, since a variable fixing can change the current and the subsequent slots in the vars array
    */
   for( v = nvars - 1; v >= 0; --v )
   {
      SCIP_VAR* var;
      SCIP_Real bound;
      SCIP_Real obj;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      var = vars[v];
      assert(var != NULL);

      /* don't perform dual presolving operations on deleted variables */
      if( SCIPvarIsDeleted(var) )
         continue;

      /* ignore already fixed variables (use feasibility tolerance since this is used in SCIPfixVar() */
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         continue;

      obj = SCIPvarGetObj(var);

      /* if the objective coefficient of the variable is 0 and it may be rounded both
       * up and down, then fix it to the closest feasible value to 0 */
      if( SCIPisZero(scip, obj) && SCIPvarMayRoundDown(var) && SCIPvarMayRoundUp(var) )
      {
         SCIP_Real roundbound;

         bound = SCIPvarGetLbGlobal(var);
         if( SCIPisLT(scip, bound, 0.0) )
         {
            if( SCIPisLE(scip, 0.0, SCIPvarGetUbGlobal(var)) )
               bound = 0.0;
            else
            {
               /* try to take an integer value, only for polishing */
               roundbound = SCIPfloor(scip, SCIPvarGetUbGlobal(var));

               if( roundbound < bound )
                  bound = SCIPvarGetUbGlobal(var);
               else
                  bound = roundbound;
            }
         }
         else
         {
            /* try to take an integer value, only for polishing */
            roundbound = SCIPceil(scip, bound);

            if( roundbound < SCIPvarGetUbGlobal(var) )
               bound = roundbound;
         }
         SCIPdebugMsg(scip, "fixing variable <%s> with objective 0 to %g\n", SCIPvarGetName(var), bound);
      }
      else
      {
         /* if it is always possible to round variable in direction of objective value, fix it to its proper bound */
         if( SCIPvarMayRoundDown(var) && !SCIPisNegative(scip, obj) )
         {
            bound = SCIPvarGetLbGlobal(var);
            if ( SCIPisInfinity(scip, -bound) )
            {
               /* variable can be fixed to -infinity */
               if ( SCIPgetStage(scip) > SCIP_STAGE_PRESOLVING )
               {
                  /* Fixing variables to infinity is not allowed after presolving, since LP-solvers cannot handle this
                   * consistently. We thus have to ignore this (should better be handled in presolving). */
                  continue;
               }
               if ( SCIPisZero(scip, obj) && SCIPvarGetNLocksUp(var) == 1 )
               {
                  /* Variable is only contained in one constraint: we hope that the corresponding constraint handler is
                   * clever enough to set/aggregate the variable to something more useful than -infinity and do nothing
                   * here. */
                  continue;
               }
            }
            SCIPdebugMsg(scip, "fixing variable <%s> with objective %g and %d uplocks to lower bound %g\n",
               SCIPvarGetName(var), SCIPvarGetObj(var), SCIPvarGetNLocksUp(var), bound);
         }
         else if( SCIPvarMayRoundUp(var) && !SCIPisPositive(scip, obj) )
         {
            bound = SCIPvarGetUbGlobal(var);
            if ( SCIPisInfinity(scip, bound) )
            {
               /* variable can be fixed to infinity */
               if ( SCIPgetStage(scip) > SCIP_STAGE_PRESOLVING )
               {
                  /* Fixing variables to infinity is not allowed after presolving, since LP-solvers cannot handle this
                   * consistently. We thus have to ignore this (should better be handled in presolving). */
                  continue;
               }
               if ( SCIPisZero(scip, obj) && SCIPvarGetNLocksDown(var) == 1 )
               {
                  /* Variable is only contained in one constraint: we hope that the corresponding constraint handler is
                   * clever enough to set/aggregate the variable to something more useful than +infinity and do nothing
                   * here */
                  continue;
               }
            }
            SCIPdebugMsg(scip, "fixing variable <%s> with objective %g and %d downlocks to upper bound %g\n",
               SCIPvarGetName(var), SCIPvarGetObj(var), SCIPvarGetNLocksDown(var), bound);
         }
         else
            continue;
      }

      if( SCIPisInfinity(scip, REALABS(bound)) && !SCIPisZero(scip, obj) )
      {
         SCIPdebugMsg(scip, " -> unbounded fixing\n");
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
            "problem infeasible or unbounded: variable <%s> with objective %.15g can be made infinitely %s\n",
            SCIPvarGetName(var), SCIPvarGetObj(var), bound < 0.0 ? "small" : "large");
         *unbounded = TRUE;
         return SCIP_OKAY;
      }

      /* apply the fixing */
      SCIPdebugMsg(scip, "apply fixing of variable %s to %g\n", SCIPvarGetName(var), bound);
      SCIP_CALL( SCIPfixVar(scip, var, bound, &infeasible, &fixed) );

      if( infeasible )
      {
         SCIPdebugMsg(scip, " -> infeasible fixing\n");
         *cutoff = TRUE;
         return SCIP_OKAY;
      }

      assert(fixed || (SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPisFeasEQ(scip, bound, SCIPvarGetLbLocal(var))
            && SCIPisFeasEQ(scip, bound, SCIPvarGetUbLocal(var))));
      (*nfixedvars)++;
   }

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods
 *
 * @{
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyDualfix)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropDualfix(scip) );

   return SCIP_OKAY;
}

/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolDualfix)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   SCIP_Bool unbounded;
   int oldnfixedvars;

   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   cutoff = FALSE;
   unbounded = FALSE;
   oldnfixedvars = *nfixedvars;

   SCIP_CALL( performDualfix(scip, nfixedvars, &unbounded, &cutoff) );

   /* evaluate propagation result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( unbounded )
      *result = SCIP_UNBOUNDED;
   else if( *nfixedvars > oldnfixedvars )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecDualfix)
{  /*lint --e{715}*/
   int nfixedvars;
   SCIP_Bool cutoff;
   SCIP_Bool unbounded;

   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /** @warning Don't run in probing or in repropagation since this can lead to wrong conclusion
    *
    *  do not run if propagation w.r.t. current objective is not allowed
    */
   if( SCIPinProbing(scip) || SCIPinRepropagation(scip) || !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   cutoff = FALSE;
   unbounded = FALSE;
   nfixedvars = 0;

   SCIP_CALL( performDualfix(scip, &nfixedvars, &unbounded, &cutoff) );

   /* evaluate propagation result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( unbounded )
      *result = SCIP_UNBOUNDED;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the dual fixing propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropDualfix(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, propExecDualfix, NULL) );
   assert(prop != NULL);

   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyDualfix) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolDualfix, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );

   return SCIP_OKAY;
}

/**@} */
