/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_probing.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for the probing mode
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/debug.h"
#include "scip/heur.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_relax.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/relax.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/struct_lp.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_tree.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"

/** returns whether we are in probing mode; probing mode is activated via SCIPstartProbing() and stopped
 *  via SCIPendProbing()
 *
 *  @return TRUE, if SCIP is currently in probing mode, otherwise FALSE
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Bool SCIPinProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPinProbing", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeProbing(scip->tree);
}

/** initiates probing, making methods SCIPnewProbingNode(), SCIPbacktrackProbing(), SCIPchgVarLbProbing(),
 *  SCIPchgVarUbProbing(), SCIPfixVarProbing(), SCIPpropagateProbing(), and SCIPsolveProbingLP() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note The collection of variable statistics is turned off during probing. If these statistics should be collected
 *        during probing use the method SCIPenableVarHistory() to turn the collection explicitly on.
 */
SCIP_RETCODE SCIPstartProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPstartProbing", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("already in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   if( scip->lp != NULL && SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("cannot start probing while in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* use a different separation storage for probing mode; otherwise SCIP will remove the cuts that are currently in the
    * separation storage after solving an LP in probing mode
    */
   if( scip->sepastore != NULL )
   {
      assert(scip->sepastoreprobing != NULL);
      SCIPswapPointers((void**)&scip->sepastore, (void**)&scip->sepastoreprobing);
   }

   SCIP_CALL( SCIPtreeStartProbing(scip->tree, scip->mem->probmem, scip->set, scip->lp, scip->relaxation, scip->transprob, FALSE) );

   /* disables the collection of any statistic for a variable */
   SCIPstatDisableVarHistory(scip->stat);

   return SCIP_OKAY;
}

/** creates a new probing sub node, whose changes can be undone by backtracking to a higher node in the probing path
 *  with a call to SCIPbacktrackProbing();
 *  using a sub node for each set of probing bound changes can improve conflict analysis
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPnewProbingNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPnewProbingNode", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   retcode = SCIPtreeCreateProbingNode(scip->tree, scip->mem->probmem, scip->set, scip->lp);

   if( retcode == SCIP_MAXDEPTHLEVEL )
   {
      SCIPwarningMessage(scip, "probing reached maximal depth; it should be stopped\n");
   }
   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/** returns the current probing depth
 *
 *  @return the probing depth, i.e. the number of probing sub nodes existing in the probing path
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetProbingDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetProbingDepth", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      SCIPABORT();
      return -1; /*lint !e527*/
   }

   return SCIPtreeGetProbingDepth(scip->tree);
}

/** undoes all changes to the problem applied in probing up to the given probing depth;
 *  the changes of the probing node of the given probing depth are the last ones that remain active;
 *  changes that were applied before calling SCIPnewProbingNode() cannot be undone
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPbacktrackProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbacktrackProbing", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   if( probingdepth < 0 || probingdepth > SCIPtreeGetProbingDepth(scip->tree) )
   {
      SCIPerrorMessage("backtracking probing depth %d out of current probing range [0,%d]\n",
         probingdepth, SCIPtreeGetProbingDepth(scip->tree));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtreeBacktrackProbing(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
         scip->origprob, scip->lp, scip->primal, scip->branchcand, scip->eventqueue, scip->eventfilter,
         scip->cliquetable, probingdepth) );

   return SCIP_OKAY;
}

/** quits probing and resets bounds and constraints to the focus node's environment
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPendProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPendProbing", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   /* switch back from probing to normal operation mode and restore variables and constraints to focus node */
   SCIP_CALL( SCIPtreeEndProbing(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
         scip->transprob, scip->origprob, scip->lp, scip->relaxation, scip->primal,
         scip->branchcand, scip->eventqueue, scip->eventfilter, scip->cliquetable) );

   /* enables the collection of statistics for a variable */
   SCIPstatEnableVarHistory(scip->stat);

   /* switch to the original separation storage */
   if( scip->sepastore != NULL )
   {
      assert(scip->sepastoreprobing != NULL);
      SCIPswapPointers((void**)&scip->sepastore, (void**)&scip->sepastoreprobing);
      assert(SCIPsepastoreGetNCuts(scip->sepastoreprobing) == 0);
   }

   return SCIP_OKAY;
}

/** injects a change of variable's lower bound into current probing node; the same can also be achieved with a call to
 *  SCIPchgVarLb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarLbProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarLbProbing", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPnodeGetType(SCIPtreeGetCurrentNode(scip->tree)) == SCIP_NODETYPE_PROBINGNODE);

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* ignore tightenings of lower bounds to +infinity during solving process */
   if( SCIPisInfinity(scip, newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore lower bound tightening for %s from %e to +infinity\n", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var));
#endif
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
         scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable,
         var, newbound, SCIP_BOUNDTYPE_LOWER, TRUE) );

   return SCIP_OKAY;
}

/** injects a change of variable's upper bound into current probing node; the same can also be achieved with a call to
 *  SCIPchgVarUb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarUbProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarUbProbing", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPnodeGetType(SCIPtreeGetCurrentNode(scip->tree)) == SCIP_NODETYPE_PROBINGNODE);

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* ignore tightenings of upper bounds to -infinity during solving process */
   if( SCIPisInfinity(scip, -newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore upper bound tightening for %s from %e to -infinity\n", SCIPvarGetName(var),
         SCIPvarGetUbLocal(var));
#endif
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
         scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable,
         var, newbound, SCIP_BOUNDTYPE_UPPER, TRUE) );

   return SCIP_OKAY;
}

/** gets variable's objective value in current probing
 *
 *  @return the variable's objective value in current probing.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetVarObjProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarObjProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALID;
   }

   return SCIPvarGetObjLP(var);
}

/** injects a change of variable's bounds into current probing node to fix the variable to the specified value;
 *  the same can also be achieved with a call to SCIPfixVar(), but in this case, the bound changes would be treated
 *  like deductions instead of branching decisions
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPfixVarProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             fixedval            /**< value to fix variable to */
   )
{
   SCIP_Real fixlb;
   SCIP_Real fixub;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPfixVarProbing", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPnodeGetType(SCIPtreeGetCurrentNode(scip->tree)) == SCIP_NODETYPE_PROBINGNODE);

   /* we adjust the fixing value here and compare the old bound with the adjusted values because otherwise,
    * it might happen that the unadjusted value is better and we add the boundchange,
    * but within SCIPnodeAddBoundchg() the bounds are adjusted - using the feasibility epsilon for integer variables -
    * and it is asserted, that the bound is still better than the old one which might then be incorrect.
    */
   fixlb = fixedval;
   fixub = fixedval;
   SCIPvarAdjustLb(var, scip->set, &fixlb);
   SCIPvarAdjustUb(var, scip->set, &fixub);
   assert(SCIPsetIsEQ(scip->set, fixlb, fixub));

   if( SCIPsetIsGT(scip->set, fixlb, SCIPvarGetLbLocal(var)) )
   {
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
            scip->cliquetable, var, fixlb, SCIP_BOUNDTYPE_LOWER, TRUE) );
   }
   if( SCIPsetIsLT(scip->set, fixub, SCIPvarGetUbLocal(var)) )
   {
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable,
            var, fixub, SCIP_BOUNDTYPE_UPPER, TRUE) );
   }

   return SCIP_OKAY;
}

/** changes (column) variable's objective value during probing mode
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @pre The variable needs to be a column variable.
 */
SCIP_RETCODE SCIPchgVarObjProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective for */
   SCIP_Real             newobj              /**< new objective function value */
   )
{
   SCIP_NODE* node;
   SCIP_Real oldobj;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarObjProbing", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   /* get current probing node */
   node = SCIPtreeGetCurrentNode(scip->tree);
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

   /* get old objective function value */
   oldobj = SCIPvarGetObj(var);

   if( SCIPisEQ(scip, oldobj, newobj) )
      return SCIP_OKAY;

   if( node->data.probingnode->nchgdobjs == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &node->data.probingnode->origobjvars, 1) ); /*lint !e506*/
      SCIP_CALL( SCIPallocMemoryArray(scip, &node->data.probingnode->origobjvals, 1) ); /*lint !e506*/
   }
   else
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &node->data.probingnode->origobjvars, node->data.probingnode->nchgdobjs + 1) ); /*lint !e776*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &node->data.probingnode->origobjvals, node->data.probingnode->nchgdobjs + 1) ); /*lint !e776*/
   }

   node->data.probingnode->origobjvars[node->data.probingnode->nchgdobjs] = var;
   node->data.probingnode->origobjvals[node->data.probingnode->nchgdobjs] = oldobj;
   ++node->data.probingnode->nchgdobjs;
   ++scip->tree->probingsumchgdobjs;

   assert(SCIPtreeProbingObjChanged(scip->tree) == SCIPlpDivingObjChanged(scip->lp));

   /* inform tree and LP that the objective was changed and invalidate the LP's cutoff bound, since this has nothing to
    * do with the current objective value anymore; the cutoff bound is reset in SCIPendProbing()
    */
   if( !SCIPtreeProbingObjChanged(scip->tree) )
   {
      SCIP_CALL( SCIPlpSetCutoffbound(scip->lp, scip->set, scip->transprob, SCIPsetInfinity(scip->set)) );

      SCIPtreeMarkProbingObjChanged(scip->tree);
      SCIPlpMarkDivingObjChanged(scip->lp);
   }
   assert(SCIPisInfinity(scip, scip->lp->cutoffbound));

   /* perform the objective change */
   SCIP_CALL( SCIPvarChgObj(var, scip->mem->probmem, scip->set,  scip->transprob, scip->primal, scip->lp, scip->eventqueue, newobj) );

   return SCIP_OKAY;
}

/** returns whether the objective function has changed during probing mode
 *
 *  @return \ref TRUE if objective has changed, \ref FALSE otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Bool SCIPisObjChangedProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisObjChangedProbing", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->tree != NULL && SCIPinProbing(scip) && SCIPtreeProbingObjChanged(scip->tree);
}

/** applies domain propagation on the probing sub problem, that was changed after SCIPstartProbing() was called;
 *  the propagated domains of the variables can be accessed with the usual bound accessing calls SCIPvarGetLbLocal()
 *  and SCIPvarGetUbLocal(); the propagation is only valid locally, i.e. the local bounds as well as the changed
 *  bounds due to SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), and SCIPfixVarProbing() are used for propagation
 *
 *  @note Conflict analysis can run if the propagation finds infeasibilities. SCIPpropagateProbing can even find
 *  globally valid bound changes. For this reason, the function restores the original objective (i.e. undoes the changes
 *  done by SCIPchgVarObjProbing before performing the propagation, as the propagators don't know that the objective
 *  might have changed. Thus, SCIPpropagateProbing can have an effect on the problem after probing ends.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPpropagateProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the probing node can be cut off */
   SCIP_Longint*         ndomredsfound       /**< pointer to store the number of domain reductions found, or NULL */
   )
{
   SCIP_VAR** objchgvars;
   SCIP_Real* objchgvals;
   SCIP_Bool changedobj;
   int nobjchg;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPpropagateProbing", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   objchgvars = NULL;
   objchgvals = NULL;
   changedobj = FALSE;
   nobjchg = 0;

   /* undo objective changes if we want to propagate during probing */
   if( scip->tree->probingobjchanged )
   {
      SCIP_VAR** vars;
      int nvars;
      int i;

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      SCIP_CALL( SCIPallocBufferArray(scip, &objchgvals, MIN(nvars, scip->tree->probingsumchgdobjs)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &objchgvars, MIN(nvars, scip->tree->probingsumchgdobjs)) );
      nobjchg = 0;

      for( i = 0; i < nvars; ++i )
      {
         if( !SCIPisEQ(scip, vars[i]->unchangedobj, SCIPgetVarObjProbing(scip, vars[i])) )
         {
            objchgvars[nobjchg] = vars[i];
            objchgvals[nobjchg] = SCIPgetVarObjProbing(scip, vars[i]);
            ++nobjchg;

            SCIP_CALL( SCIPvarChgObj(vars[i], scip->mem->probmem, scip->set, scip->transprob, scip->primal, scip->lp,
                  scip->eventqueue, vars[i]->unchangedobj) );
         }
      }
      assert(nobjchg <= scip->tree->probingsumchgdobjs);

      SCIPlpUnmarkDivingObjChanged(scip->lp);
      scip->tree->probingobjchanged = FALSE;
      changedobj = TRUE;
   }

   if( ndomredsfound != NULL )
      *ndomredsfound = -(scip->stat->nprobboundchgs + scip->stat->nprobholechgs);

   SCIP_CALL( SCIPpropagateDomains(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->conflict, scip->cliquetable,
         SCIPgetDepth(scip), maxproprounds, SCIP_PROPTIMING_ALWAYS, cutoff) );

   if( ndomredsfound != NULL )
      *ndomredsfound += scip->stat->nprobboundchgs + scip->stat->nprobholechgs;

   /* restore old objective function */
   if( changedobj )
   {
      int i;

      assert(objchgvars != NULL);
      assert(objchgvals != NULL);

      SCIPlpMarkDivingObjChanged(scip->lp);
      scip->tree->probingobjchanged = TRUE;

      for( i = 0; i < nobjchg; ++i )
      {
         SCIP_CALL( SCIPvarChgObj(objchgvars[i], scip->mem->probmem, scip->set,  scip->transprob, scip->primal,
               scip->lp, scip->eventqueue, objchgvals[i]) );
      }

      SCIPfreeBufferArray(scip, &objchgvars);
      SCIPfreeBufferArray(scip, &objchgvals);
   }

   return SCIP_OKAY;
}

/** applies domain propagation on the probing sub problem, that was changed after SCIPstartProbing() was called;
 *  only propagations of the binary variables fixed at the current probing node that are triggered by the implication
 *  graph and the clique table are applied;
 *  the propagated domains of the variables can be accessed with the usual bound accessing calls SCIPvarGetLbLocal()
 *  and SCIPvarGetUbLocal(); the propagation is only valid locally, i.e. the local bounds as well as the changed
 *  bounds due to SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), and SCIPfixVarProbing() are used for propagation
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPpropagateProbingImplications(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing node can be cut off */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPpropagateProbingImplications", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnodePropagateImplics(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
         scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, cutoff) );

   return SCIP_OKAY;
}

/** solves the LP at the current probing node (cannot be applied at preprocessing stage) with or without pricing */
static
SCIP_RETCODE solveProbingLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool             pricing,            /**< should pricing be applied? */
   SCIP_Bool             pretendroot,        /**< should the pricers be called as if we are at the root node? */
   SCIP_Bool             displayinfo,        /**< should info lines be displayed after each pricing round? */
   int                   maxpricerounds,     /**< maximal number of pricing rounds (-1: no limit) */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing LP was infeasible or the objective
                                              *   limit was reached (or NULL, if not needed) */
   )
{
   SCIP_Bool initcutoff;

   assert(lperror != NULL);
   assert(SCIPtreeIsFocusNodeLPConstructed(scip->tree));

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPtreeGetCurrentDepth(scip->tree) > 0);

   SCIP_CALL( SCIPinitConssLP(scip->mem->probmem, scip->set, scip->sepastore, scip->cutpool, scip->stat, scip->transprob,
         scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->eventfilter,
         scip->cliquetable, FALSE, FALSE, &initcutoff) );

   if( initcutoff )
   {
      if( cutoff != NULL )
         *cutoff = TRUE;

      return SCIP_OKAY;
   }
   else if( cutoff != NULL )
      *cutoff = FALSE;

   /* load the LP state (if necessary) */
   SCIP_CALL( SCIPtreeLoadProbingLPState(scip->tree, scip->mem->probmem, scip->set, scip->transprob, scip->eventqueue, scip->lp) );

   SCIPlpSetIsRelax(scip->lp, TRUE);

   /* solve probing LP */
   SCIP_CALL( SCIPlpSolveAndEval(scip->lp, scip->set, scip->messagehdlr, scip->mem->probmem, scip->stat,
         scip->eventqueue, scip->eventfilter, scip->transprob, (SCIP_Longint)itlim, FALSE, FALSE, FALSE, lperror) );

   assert((*lperror) || SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_NOTSOLVED);

   /* mark the probing node to have a solved LP */
   if( !(*lperror) )
   {
      SCIP_CALL( SCIPtreeMarkProbingNodeHasLP(scip->tree, scip->mem->probmem, scip->lp) );

      /* call pricing */
      if( pricing )
      {
         SCIP_Bool mustsepa;
         int npricedcolvars;
         SCIP_Bool result;

            mustsepa = FALSE;
            SCIP_CALL( SCIPpriceLoop(scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->transprob,
                  scip->origprob, scip->primal, scip->tree, scip->reopt, scip->lp, scip->pricestore, scip->sepastore, scip->cutpool,
                  scip->branchcand, scip->eventqueue, scip->eventfilter, scip->cliquetable, pretendroot, displayinfo,
                  maxpricerounds, &npricedcolvars, &mustsepa, lperror, &result) );

         /* mark the probing node again to update the LP size in the node and the tree path */
         if( !(*lperror) )
         {
            SCIP_CALL( SCIPtreeMarkProbingNodeHasLP(scip->tree, scip->mem->probmem, scip->lp) );
         }
      }
   }

   /* remember that probing might have changed the LPi state; this holds even if solving returned with an LP error */
   scip->tree->probingsolvedlp = TRUE;

   /* the LP is infeasible or the objective limit was reached */
   if( !(*lperror) && (SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_INFEASIBLE
         || SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_OBJLIMIT ||
         (SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_OPTIMAL
            && SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))) )
   {
      /* analyze the infeasible LP (only if all columns are in the LP and no external pricers exist) */
      if( !scip->set->misc_exactsolve && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->tree->probingobjchanged )
      {
         SCIP_CALL( SCIPconflictAnalyzeLP(scip->conflict, scip->conflictstore, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, NULL) );
      }

      if( cutoff != NULL )
         *cutoff = TRUE;
   }

   return SCIP_OKAY;
}

/** solves the LP at the current probing node (cannot be applied at preprocessing stage);
 *  no separation or pricing is applied
 *
 *  The LP has to be constructed before (you can use SCIPisLPConstructed() or SCIPconstructLP()).
 *
 *  @note if the LP is infeasible or the objective limit is reached, and if all columns are in the LP and no external
 *  pricers exist then conflict analysis will be run. This can have an effect on the problem after probing ends.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsolveProbingLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing LP was infeasible or the objective
                                              *   limit was reached (or NULL, if not needed) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveProbingLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( solveProbingLP(scip, itlim, FALSE, FALSE, FALSE, -1, lperror, cutoff) );

   return SCIP_OKAY;
}

/** solves the LP at the current probing node (cannot be applied at preprocessing stage) and applies pricing
 *  until the LP is solved to optimality; no separation is applied
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsolveProbingLPWithPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             pretendroot,        /**< should the pricers be called as if we were at the root node? */
   SCIP_Bool             displayinfo,        /**< should info lines be displayed after each pricing round? */
   int                   maxpricerounds,     /**< maximal number of pricing rounds (-1: no limit);
                                              *   a finite limit means that the LP might not be solved to optimality! */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing LP was infeasible or the objective
                                              *   limit was reached (or NULL, if not needed) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveProbingLPWithPricing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( solveProbingLP(scip, -1, TRUE, pretendroot, displayinfo, maxpricerounds, lperror, cutoff) );

   return SCIP_OKAY;
}

/** sets the LP state for the current probing node
 *
 *  @note state and norms are stored at the node and later released by SCIP; therefore, the pointers are set
 *        to NULL by the method
 *
 *  @note the pointers to state and norms must not be NULL; however, they may point to a NULL pointer if the
 *        respective information should not be set
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetProbingLPState(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LPISTATE**       lpistate,           /**< pointer to LP state information (like basis information) */
   SCIP_LPINORMS**       lpinorms,           /**< pointer to LP pricing norms information */
   SCIP_Bool             primalfeas,         /**< primal feasibility when LP state information was stored */
   SCIP_Bool             dualfeas            /**< dual feasibility when LP state information was stored */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetProbingLPState", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPtreeSetProbingLPState(scip->tree, scip->mem->probmem, scip->lp, lpistate, lpinorms, primalfeas, dualfeas) );

   return SCIP_OKAY;
}

/** adds a row to the LP in the current probing node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPaddRowProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row to be added */
   )
{
   SCIP_NODE* node;
   int depth;

   assert(scip != NULL);
   assert(row != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddRowProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   /* get depth of current node */
   node = SCIPtreeGetCurrentNode(scip->tree);
   assert(node != NULL);
   depth = SCIPnodeGetDepth(node);

   SCIP_CALL( SCIPlpAddRow(scip->lp, scip->mem->probmem, scip->set, scip->eventqueue, scip->eventfilter, row, depth) );

   return SCIP_OKAY;
}


/** applies the cuts in the separation storage to the LP and clears the storage afterwards;
 *  this method can only be applied during probing; the user should resolve the probing LP afterwards
 *  in order to get a new solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPapplyCutsProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPapplyCutsProbing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsepastoreApplyCuts(scip->sepastore, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
         scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->eventfilter,
         scip->cliquetable, FALSE, SCIP_EFFICIACYCHOICE_LP, cutoff) );

   return SCIP_OKAY;
}

/** solves relaxation(s) at the current probing node (cannot be applied at preprocessing stage);
 *  no separation or pricing is applied
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsolveProbingRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether a relaxation was infeasible or the objective
                                              *   limit was reached (or NULL, if not needed) */
   )
{
   SCIP_SET* set;
   int r;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveProbingRelax", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if ( ! SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   assert( SCIPtreeGetCurrentDepth(scip->tree) > 0 );

   assert( cutoff != NULL );
   *cutoff = FALSE;

   set = scip->set;

   /* sort relaxators by priority */
   SCIPsetSortRelaxs(set);

   /* solve relaxations */
   for (r = 0; r < set->nrelaxs && !(*cutoff); ++r)
   {
      SCIP_RELAX* relax;
      SCIP_Real lowerbound;
      SCIP_RESULT result;

      lowerbound = -SCIPinfinity(scip);

      relax = set->relaxs[r];
      assert( relax != NULL );

      SCIP_CALL( SCIPrelaxExec(relax, set, scip->tree, scip->stat, SCIPtreeGetCurrentDepth(scip->tree), &lowerbound, &result) );

      switch( result )
      {
      case SCIP_CUTOFF:
         *cutoff = TRUE;
         SCIPdebugMsg(scip, " -> relaxator <%s> detected cutoff\n", SCIPrelaxGetName(relax));
         break;

      case SCIP_CONSADDED:
      case SCIP_REDUCEDDOM:
      case SCIP_SEPARATED:
      case SCIP_SUSPENDED:
         SCIPerrorMessage("The relaxator should not return <%d> within probing mode.\n", result);
         break;

      case SCIP_SUCCESS:
      case SCIP_DIDNOTRUN:
         break;

      default:
         SCIPerrorMessage("Invalid result code <%d> of relaxator <%s>\n", result, SCIPrelaxGetName(relax));
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** print statistics of probing */
char* SCIPsnprintfProbingStats(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 strbuf,             /**< string buffer */
   int                   len                 /**< length of string buffer */
   )
{
   char* ptr = strbuf;
   const int nvartypes = 4;

   assert(scip != NULL);
   assert(strbuf != NULL);

   if( SCIPinProbing(scip) )
   {
      SCIP_VAR** vars;
      int nbinvars = SCIPgetNBinVars(scip);
      int nintvars = SCIPgetNIntVars(scip);
      int nimplvars = SCIPgetNImplVars(scip);
      int nvars = SCIPgetNVars(scip);
      int vartypeend[] = {
            nbinvars,
            nbinvars + nintvars,
            nbinvars + nintvars + nimplvars,
            nvars
      };
      const char* vartypenames[] = {
            "binary",
            "integer",
            "implicit integer",
            "continuous"
      };
      int nvartypefixed[4];
      int nvarsfixed = 0;
      int depth;
      int probingdepth;
      int vartypestart = 0;
      int v;
      int p;

      vars = SCIPgetVars(scip);
      BMSclearMemoryArray(nvartypefixed, nvartypes);

      /* loop over vartypes and count fixings */
      for( p = 0; p < nvartypes; ++p )
      {
         for( v = vartypestart; v < vartypeend[p]; ++v )
         {
            if( SCIPisEQ(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v])) )
               ++nvartypefixed[p];
         }
         nvarsfixed += nvartypefixed[p];
         vartypestart = vartypeend[p];
      }

      depth = SCIPgetDepth(scip);
      probingdepth = SCIPgetProbingDepth(scip);

      ptr += SCIPsnprintf(ptr, len, "Depth: (%d total, %d probing) ", depth, probingdepth);
      ptr += SCIPsnprintf(ptr, len, "Fixed/Variables: %d / %d (", nvarsfixed, vartypeend[nvartypes - 1]);

      for( p = 0; p < nvartypes; ++p )
      {
         int ntypevars = vartypeend[p] - (p == 0 ? 0 : vartypeend[p - 1]);
         ptr += SCIPsnprintf(ptr, len, "%d / %d %s%s", nvartypefixed[p], ntypevars, vartypenames[p], p < (nvartypes - 1) ? ", " : ")");
      }
   }
   else
   {
      (void) SCIPsnprintf(strbuf, len, "Not in probing");
   }

   return strbuf;
}

/** gets the candidate score and preferred rounding direction for a candidate variable */
SCIP_RETCODE SCIPgetDivesetScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< general diving settings */
   SCIP_DIVETYPE         divetype,           /**< represents different methods for a dive set to explore the next children */
   SCIP_VAR*             divecand,           /**< the candidate for which the branching direction is requested */
   SCIP_Real             divecandsol,        /**< LP solution value of the candidate */
   SCIP_Real             divecandfrac,       /**< fractionality of the candidate */
   SCIP_Real*            candscore,          /**< pointer to store the candidate score */
   SCIP_Bool*            roundup             /**< pointer to store whether preferred direction for diving is upwards */
   )
{
   assert(scip != NULL);
   assert(candscore != NULL);
   assert(roundup != NULL);
   assert(divecand != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetDivesetScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPdivesetGetScore(diveset, scip->set, divetype, divecand, divecandsol, divecandfrac, candscore,
         roundup) );

   return SCIP_OKAY;
}

/** update diveset LP statistics, should be called after every LP solved by this diving heuristic */
void SCIPupdateDivesetLPStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_Longint          niterstoadd,        /**< additional number of LP iterations to be added */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   )
{
   assert(scip != NULL);
   assert(diveset != NULL);

   SCIPdivesetUpdateLPStats(diveset, scip->stat, niterstoadd, divecontext);
}

/** update diveset statistics and global diveset statistics */
void SCIPupdateDivesetStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< diveset to be reset */
   int                   nprobingnodes,      /**< the number of probing nodes explored this time */
   int                   nbacktracks,        /**< the number of backtracks during probing this time */
   SCIP_Longint          nsolsfound,         /**< the number of solutions found */
   SCIP_Longint          nbestsolsfound,     /**< the number of best solutions found */
   SCIP_Longint          nconflictsfound,    /**< number of new conflicts found this time */
   SCIP_Bool             leavewassol,        /**< was a solution found at the leaf? */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   )
{
   assert(scip != NULL);
   assert(diveset != NULL);
   assert(SCIPinProbing(scip));

   SCIPdivesetUpdateStats(diveset, scip->stat, SCIPgetDepth(scip), nprobingnodes, nbacktracks, nsolsfound,
         nbestsolsfound, nconflictsfound, leavewassol, divecontext);
}

/** enforces a probing/diving solution by suggesting bound changes that maximize the score w.r.t. the current diving settings
 *
 *  the process is guided by the enforcement priorities of the constraint handlers and the scoring mechanism provided by
 *  the dive set.
 *  Constraint handlers may suggest diving bound changes in decreasing order of their enforcement priority, based on the
 *  solution values in the solution @p sol and the current local bounds of the variables. A diving bound change
 *  is a triple (variable,branching direction,value) and is used inside SCIPperformGenericDivingAlgorithm().
 *
 *  After a successful call, SCIP holds two arrays of suggested dive bound changes, one for the preferred child
 *  and one for the alternative.
 *
 *  @see SCIPgetDiveBoundChangeData() for retrieving the dive bound change suggestions.
 *
 *  The method stops after the first constraint handler was successful
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetDiveBoundChanges(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< diving settings to control scoring */
   SCIP_SOL*             sol,                /**< current solution of diving mode */
   SCIP_Bool*            success,            /**< pointer to store whether constraint handler successfully found a variable */
   SCIP_Bool*            infeasible          /**< pointer to store whether the current node was detected to be infeasible */
   )
{
   int i;

   assert(scip != NULL);
   assert(diveset != NULL);
   assert(SCIPinProbing(scip));
   assert(infeasible != NULL);
   assert(success != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetDiveBoundChanges", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *success = FALSE;
   *infeasible = FALSE;

   /* we invalidate the previously stored bound changes */
   SCIPclearDiveBoundChanges(scip);

   /* loop over constraint handlers until a constraint handler successfully found a variable/value assignment for proceeding
    * or a constraint handler detected the infeasibility of the local node
    */
   for( i = 0; i < scip->set->nconshdlrs && !(*success || *infeasible); ++i )
   {
      SCIP_CALL( SCIPconshdlrGetDiveBoundChanges(scip->set->conshdlrs_enfo[i], scip->set, diveset, sol,
            success, infeasible) );
   }

#ifndef NDEBUG
   /* check if the constraint handler correctly assigned values to the dive set */
   if( *success )
   {
      SCIP_VAR** bdchgvars;
      SCIP_BRANCHDIR* bdchgdirs;
      SCIP_Real* values;
      int nbdchanges;
      SCIPtreeGetDiveBoundChangeData(scip->tree, &bdchgvars, &bdchgdirs, &values, &nbdchanges, TRUE);
      assert(nbdchanges > 0);
      SCIPtreeGetDiveBoundChangeData(scip->tree, &bdchgvars, &bdchgdirs, &values, &nbdchanges, FALSE);
      assert(nbdchanges > 0);
   }
#endif

   return SCIP_OKAY;
}

/** adds a diving bound change to the diving bound change storage of SCIP together with the information if this is a
 *  bound change for the preferred direction or not
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPaddDiveBoundChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to apply the bound change to */
   SCIP_BRANCHDIR        dir,                /**< direction of the bound change */
   SCIP_Real             value,              /**< value to adjust this variable bound to */
   SCIP_Bool             preferred           /**< is this a bound change for the preferred child? */
   )
{
   assert(scip->tree != NULL);
   assert(scip->mem->probmem != NULL);
   assert(SCIPinProbing(scip));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddDiveBoundChange", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPtreeAddDiveBoundChange(scip->tree, scip->mem->probmem, var, dir, value, preferred) );

   return SCIP_OKAY;
}

/** get the dive bound change data for the preferred or the alternative direction
 *
 *   @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void SCIPgetDiveBoundChangeData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           variables,          /**< pointer to store variables for the specified direction */
   SCIP_BRANCHDIR**      directions,         /**< pointer to store the branching directions */
   SCIP_Real**           values,             /**< pointer to store bound change values */
   int*                  ndivebdchgs,        /**< pointer to store the number of dive bound changes */
   SCIP_Bool             preferred           /**< should the dive bound changes for the preferred child be output? */
   )
{
   assert(variables != NULL);
   assert(directions != NULL);
   assert(values != NULL);
   assert(ndivebdchgs != NULL);
   assert(SCIPinProbing(scip));

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetDiveBoundChangeData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPtreeGetDiveBoundChangeData(scip->tree, variables, directions, values, ndivebdchgs, preferred);
}

/** clear the dive bound change data structures
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void SCIPclearDiveBoundChanges(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip->tree != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPclearDiveBoundChanges", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPtreeClearDiveBoundChanges(scip->tree);
}
