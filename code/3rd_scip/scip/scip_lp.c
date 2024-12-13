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

/**@file   scip_lp.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for the LP relaxation, rows and columns
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
#include "lpi/lpi.h"
#include "scip/conflict.h"
#include "scip/debug.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_tree.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/set.h"
#include "scip/solve.h"
#include "scip/struct_lp.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_tree.h"
#include "scip/tree.h"
#include "scip/var.h"

/** returns, whether the LP was or is to be solved in the current node
 *
 *  @return whether the LP was or is to be solved in the current node.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPhasCurrentNodeLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPhasCurrentNodeLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeHasCurrentNodeLP(scip->tree);
}

/** returns, whether the LP of the current node is already constructed
 *
 *  @return whether the LP of the current node is already constructed.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPisLPConstructed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisLPConstructed", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeIsFocusNodeLPConstructed(scip->tree);
}

/** makes sure that the LP of the current node is loaded and may be accessed through the LP information methods
 *
 *  @warning Contructing the LP might change the amount of variables known in the transformed problem and therefore also
 *           the variables array of SCIP (returned by SCIPgetVars() and SCIPgetVarsData()), so it might be necessary to
 *           call one of the later method after this one
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPconstructLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPconstructLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconstructCurrentLP(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->tree, scip->reopt, scip->lp, scip->pricestore, scip->sepastore, scip->cutpool, scip->branchcand,
         scip->eventqueue, scip->eventfilter, scip->cliquetable, FALSE, cutoff) );

   return SCIP_OKAY;
}

/** makes sure that the LP of the current node is flushed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPflushLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPflushLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpFlush(scip->lp, scip->mem->probmem, scip->set, scip->transprob, scip->eventqueue) );

   return SCIP_OKAY;
}

/** gets solution status of current LP
 *
 *  @return the solution status of current LP.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_LPSOLSTAT SCIPgetLPSolstat(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPSolstat", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetSolstat(scip->lp);
   else
      return SCIP_LPSOLSTAT_NOTSOLVED;
}

/** returns whether the current LP solution passed the primal feasibility check
 *
 *  @return whether the current LP solution passed the primal feasibility check.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPisLPPrimalReliable(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisLPPrimalReliable", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpIsPrimalReliable(scip->lp);
}

/** returns whether the current LP solution passed the dual feasibility check
 *
 *  @returns whether the current LP solution passed the dual feasibility check.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPisLPDualReliable(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisLPDualReliable", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpIsDualReliable(scip->lp);
}

/** returns whether the current lp is a relaxation of the current problem and its optimal objective value is a local lower bound
 *
 *  @return whether the current lp is a relaxation of the current problem and its optimal objective value is a local lower bound.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPisLPRelax(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisLPRelax", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpIsRelax(scip->lp);
}

/** gets objective value of current LP (which is the sum of column and loose objective value)
 *
 *  @return the objective value of current LP (which is the sum of column and loose objective value).
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This method returns the objective value of the current LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the LP solution status returned by
 *        SCIPgetLPSolstat() is SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetLPObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetObjval(scip->lp, scip->set, scip->transprob);
}

/** gets part of objective value of current LP that results from COLUMN variables only
 *
 *  @return the part of objective value of current LP that results from COLUMN variables only.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetLPColumnObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPColumnObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetColumnObjval(scip->lp);
}

/** gets part of objective value of current LP that results from LOOSE variables only
 *
 *  @return part of objective value of current LP that results from LOOSE variables only.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetLPLooseObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPLooseObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetLooseObjval(scip->lp, scip->set, scip->transprob);
}

/** gets the global pseudo objective value; that is all variables set to their best  (w.r.t. the objective
 *  function) global bound
 *
 *  @return the global pseudo objective value; that is all variables set to their best  (w.r.t. the objective
 *  function) global bound.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetGlobalPseudoObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetGlobalPseudoObjval", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetGlobalPseudoObjval(scip->lp, scip->set, scip->transprob);
}

/** gets the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound
 *
 *  @return the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetPseudoObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPseudoObjval", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetPseudoObjval(scip->lp, scip->set, scip->transprob);
}

/** returns whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound
 *
 *  @return whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPisRootLPRelax(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisRootLPRelax", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpIsRootLPRelax(scip->lp);
}

/** gets the objective value of the root node LP or SCIP_INVALID if the root node LP was not (yet) solved
 *
 *  @return the objective value of the root node LP or SCIP_INVALID if the root node LP was not (yet) solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetLPRootObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPRootObjval", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetRootObjval(scip->lp);
}

/** gets part of the objective value of the root node LP that results from COLUMN variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 *
 *  @return the part of the objective value of the root node LP that results from COLUMN variables only;
 *  or SCIP_INVALID if the root node LP was not (yet) solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetLPRootColumnObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPRootColumnObjval", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetRootColumnObjval(scip->lp);
}

/** gets part of the objective value of the root node LP that results from LOOSE variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 *
 *  @return the part of the objective value of the root node LP that results from LOOSE variables only;
 *  or SCIP_INVALID if the root node LP was not (yet) solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetLPRootLooseObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPRootLooseObjval", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetRootLooseObjval(scip->lp);
}

/** gets current primal feasibility tolerance of LP */
SCIP_Real SCIPgetLPFeastol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPFeastol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpGetFeastol(scip->lp);
}

/** sets primal feasibility tolerance of LP */
void SCIPsetLPFeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newfeastol          /**< new primal feasibility tolerance for LP */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPsetLPFeastol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPlpSetFeastol(scip->lp, scip->set, newfeastol);
}

/** resets primal feasibility tolerance of LP
 *
 * Sets primal feasibility tolerance to min of numerics/lpfeastolfactor * numerics/feastol and relaxfeastol.
 */
void SCIPresetLPFeastol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPresetLPFeastol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPlpResetFeastol(scip->lp, scip->set);
}

/** gets current LP columns along with the current number of LP columns
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPColsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*                  ncols               /**< pointer to store the number of LP columns, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPColsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      if( cols != NULL )
         *cols = SCIPlpGetCols(scip->lp);
      if( ncols != NULL )
         *ncols = SCIPlpGetNCols(scip->lp);
   }
   else
   {
      if( cols != NULL )
         *cols = NULL;
      if( ncols != NULL )
         *ncols = 0;
   }

   return SCIP_OKAY;
}

/** gets current LP columns
 *
 *  @return the current LP columns.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_COL** SCIPgetLPCols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetCols(scip->lp);
   else
      return NULL;
}

/** gets current number of LP columns
 *
 *  @return the current number of LP columns.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNLPCols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetNCols(scip->lp);
   else
      return 0;
}

/** gets current number of unfixed LP columns
 *
 *  @return the current number of unfixed LP columns.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNUnfixedLPCols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNUnfixedLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetNUnfixedCols(scip->lp, scip->set->num_epsilon);
   else
      return 0;
}

/** gets current LP rows along with the current number of LP rows
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPRowsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW***           rows,               /**< pointer to store the array of LP rows, or NULL */
   int*                  nrows               /**< pointer to store the number of LP rows, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPRowsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      if( rows != NULL )
         *rows = SCIPlpGetRows(scip->lp);
      if( nrows != NULL )
         *nrows = SCIPlpGetNRows(scip->lp);
   }
   else
   {
      if( rows != NULL )
         *rows = NULL;
      if( nrows != NULL )
         *nrows = 0;
   }

   return SCIP_OKAY;
}

/** gets current LP rows
 *
 *  @return the current LP rows.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_ROW** SCIPgetLPRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetRows(scip->lp);
   else
      return NULL;
}

/** gets current number of LP rows
 *
 *  @return the current number of LP rows.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetNLPRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetNRows(scip->lp);
   else
      return 0;
}

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 *
 *  @return TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPallColsInLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPallColsInLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp);
}

/** returns whether the current LP solution is basic, i.e. is defined by a valid simplex basis
 *
 *  @return whether the current LP solution is basic, i.e. is defined by a valid simplex basis.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPisLPSolBasic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisLPSolBasic", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPlpIsSolBasic(scip->lp);
}

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPBasisInd(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  basisind            /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPBasisInd", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBasisInd(scip->lp, basisind) );

   return SCIP_OKAY;
}

/** gets a row from the inverse basis matrix B^-1
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPBInvRow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coefs,              /**< array to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPBInvRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBInvRow(scip->lp, r, coefs, inds, ninds) );

   /* debug check if the coef is the r-th line of the inverse matrix B^-1 */
   SCIP_CALL( SCIPdebugCheckBInvRow(scip, r, coefs) ); /*lint !e506 !e774*/

   return SCIP_OKAY;
}

/** gets a column from the inverse basis matrix B^-1
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPBInvCol(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP
                                              *   returned by SCIPcolGetLPPos(); you have to call SCIPgetBasisInd()
                                              *   to get the array which links the B^-1 column numbers to the row and
                                              *   column numbers of the LP! c must be between 0 and nrows-1, since the
                                              *   basis has the size nrows * nrows */
   SCIP_Real*            coefs,              /**< array to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPBInvCol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBInvCol(scip->lp, c, coefs, inds, ninds) );

   return SCIP_OKAY;
}

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPBInvARow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            binvrow,            /**< row in B^-1 from prior call to SCIPgetLPBInvRow(), or NULL */
   SCIP_Real*            coefs,              /**< array to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPBInvARow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBInvARow(scip->lp, r, binvrow, coefs, inds, ninds) );

   return SCIP_OKAY;
}

/** gets a column from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A),
 *  i.e., it computes B^-1 * A_c with A_c being the c'th column of A
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPBInvACol(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   c,                  /**< column number which can be accessed by SCIPcolGetLPPos() */
   SCIP_Real*            coefs,              /**< array to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPBInvACol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBInvACol(scip->lp, c, coefs, inds, ninds) );

   return SCIP_OKAY;
}

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPsumLPRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   SCIP_Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   SCIP_Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsumLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpSumRows(scip->lp, scip->set, scip->transprob, weights, sumcoef, sumlhs, sumrhs) );

   return SCIP_OKAY;
}

/** interrupts or disables the interrupt of the currently ongoing lp solve; if the lp is not currently constructed just returns with no effect
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called in any SCIP stage
 */
SCIP_RETCODE SCIPinterruptLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             interrupt           /**< TRUE if interrupt should be set, FALSE if it should be disabled */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPinterruptLP", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( scip->lp == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPlpInterrupt(scip->lp, interrupt) );
   if( interrupt )
      scip->stat->userinterrupt = TRUE;

   return SCIP_OKAY;
}

/** writes current LP to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPwriteLP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
   SCIP_Bool cutoff;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      SCIP_CALL( SCIPconstructCurrentLP(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->tree, scip->reopt, scip->lp, scip->pricestore, scip->sepastore, scip->cutpool, scip->branchcand,
            scip->eventqueue, scip->eventfilter, scip->cliquetable, FALSE, &cutoff) );
   }

   /* we need a flushed lp to write the current lp */
   SCIP_CALL( SCIPlpFlush(scip->lp, scip->mem->probmem, scip->set, scip->transprob, scip->eventqueue) );

   SCIP_CALL( SCIPlpWrite(scip->lp, filename) );

   return SCIP_OKAY;
}

/** writes MIP relaxation of the current branch-and-bound node to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPwriteMIP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< file name */
   SCIP_Bool             genericnames,       /**< should generic names like x_i and row_j be used in order to avoid
                                              *   troubles with reserved symbols? */
   SCIP_Bool             origobj,            /**< should the original objective function be used? */
   SCIP_Bool             lazyconss           /**< output removable rows as lazy constraints? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteMIP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* we need a flushed lp to write the current mip */
   SCIP_CALL( SCIPlpFlush(scip->lp, scip->mem->probmem, scip->set, scip->transprob, scip->eventqueue) );

   SCIP_CALL( SCIPlpWriteMip(scip->lp, scip->set, scip->messagehdlr, filename, genericnames,
         origobj, scip->origprob->objsense, scip->transprob->objscale, scip->transprob->objoffset, lazyconss) );

   return SCIP_OKAY;
}

/** gets the LP interface of SCIP;
 *  with the LPI you can use all of the methods defined in lpi/lpi.h;
 *
 *  @warning You have to make sure, that the full internal state of the LPI does not change or is recovered completely
 *           after the end of the method that uses the LPI. In particular, if you manipulate the LP or its solution
 *           (e.g. by calling one of the SCIPlpiAdd...() or one of the SCIPlpiSolve...() methods), you have to check in
 *           advance with SCIPlpiWasSolved() whether the LP is currently solved. If this is the case, you have to make
 *           sure, the internal solution status is recovered completely at the end of your method. This can be achieved
 *           by getting the LPI state before applying any LPI manipulations with SCIPlpiGetState() and restoring it
 *           afterwards with SCIPlpiSetState() and SCIPlpiFreeState(). Additionally you have to resolve the LP with the
 *           appropriate SCIPlpiSolve...() call in order to reinstall the internal solution status.
 *
 *  @warning Make also sure, that all parameter values that you have changed are set back to their original values.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetLPI(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LPI**            lpi                 /**< pointer to store the LP interface */
   )
{
   assert(lpi != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPI", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   *lpi = SCIPlpGetLPI(scip->lp);

   return SCIP_OKAY;
}

/** displays quality information about the current LP solution. An LP solution need to be available; information printed
 *  is subject to what the LP solver supports
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_FREE
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note The printing process is done via the message handler system.
 */
SCIP_RETCODE SCIPprintLPSolutionQuality(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_LPI* lpi;
   SCIP_Real quality;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintLPSolutionQuality", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE) );

   switch( scip->set->stage )
   {
      case SCIP_STAGE_INIT:
      case SCIP_STAGE_PROBLEM:
      case SCIP_STAGE_TRANSFORMED:
      case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
      case SCIP_STAGE_PRESOLVED:
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "Problem not solving yet, no LP available.\n");
         return SCIP_OKAY;

      case SCIP_STAGE_SOLVING:
      case SCIP_STAGE_SOLVED:
         break;

      default:
         SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
         return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* note that after diving mode, the LPI may only have the basis information, but SCIPlpiWasSolved() can be false; in
    * this case, we will (depending on the LP solver) probably not obtain the quality measure; one solution would be to
    * store the results of SCIPlpiGetRealSolQuality() within the SCIP_LP after each LP solve; this would have the added
    * advantage, that we reduce direct access to the LPI, but it sounds potentially expensive
    */
   lpi = SCIPlpGetLPI(scip->lp);
   assert(lpi != NULL);

   SCIP_CALL( SCIPlpiGetRealSolQuality(lpi, SCIP_LPSOLQUALITY_ESTIMCONDITION, &quality) );
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Basis matrix condition (estimated): ");
   if( quality != SCIP_INVALID ) /*lint !e777*/
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "%.6e\n", quality);
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "not available\n");

   SCIP_CALL( SCIPlpiGetRealSolQuality(lpi, SCIP_LPSOLQUALITY_EXACTCONDITION, &quality) );
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Basis matrix condition (exact):     ");
   if( quality != SCIP_INVALID ) /*lint !e777*/
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "%.6e\n", quality);
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "not available\n");

   return SCIP_OKAY;
}

/** compute relative interior point to current LP
 *  @see SCIPlpComputeRelIntPoint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcomputeLPRelIntPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             relaxrows,          /**< should the rows be relaxed */
   SCIP_Bool             inclobjcutoff,      /**< should a row for the objective cutoff be included */
   SCIP_Real             timelimit,          /**< time limit for LP solver */
   int                   iterlimit,          /**< iteration limit for LP solver */
   SCIP_SOL**            point               /**< relative interior point on exit */
   )
{
   SCIP_Real* pointvals;
   SCIP_Bool success;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcomputeLPRelIntPoint", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert(scip != NULL);
   assert(scip->lp != NULL);
   assert(point != NULL);

   *point = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &pointvals, SCIPlpGetNCols(scip->lp)) );

   SCIP_CALL( SCIPlpComputeRelIntPoint(scip->set, scip->messagehdlr, scip->lp, scip->transprob,
         relaxrows, inclobjcutoff, timelimit, iterlimit, pointvals, &success) );

   /* if successful, create new solution with point values */
   if( success )
   {
      int i;

      SCIP_CALL( SCIPcreateSol(scip, point, NULL) );

      for( i = 0; i < SCIPlpGetNCols(scip->lp); ++i )
      {
         SCIP_CALL( SCIPsetSolVal(scip, *point, SCIPcolGetVar(SCIPlpGetCols(scip->lp)[i]), pointvals[i]) );
      }
   }

   SCIPfreeBufferArray(scip, &pointvals);

   return SCIP_OKAY;
}

/*
 * LP column methods
 */

/** returns the reduced costs of a column in the last (feasible) LP
 *
 *  @return the reduced costs of a column in the last (feasible) LP
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @note calling this method in SCIP_STAGE_SOLVED is only recommended to experienced users and should only be called
 *        for pure LP instances (without presolving)
 *
 *  @note The return value of this method should be used carefully if the dual feasibility check was explictely disabled.
 */
SCIP_Real SCIPgetColRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetColRedcost", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("cannot get reduced costs, because node LP is not processed\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }

   return SCIPcolGetRedcost(col, scip->stat, scip->lp);
}


/** returns the Farkas coefficient of a column in the last (infeasible) LP
 *
 *  @return the Farkas coefficient of a column in the last (infeasible) LP
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetColFarkasCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetColFarkasCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("cannot get Farkas coeff, because node LP is not processed\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }

   return SCIPcolGetFarkasCoef(col, scip->stat, scip->lp);
}

/** marks a column to be not removable from the LP in the current node
 *
 *  @pre this method can be called in the following stage of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void SCIPmarkColNotRemovableLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPmarkColNotRemovableLocal", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPcolMarkNotRemovableLocal(col, scip->stat);
}

/*
 * LP row methods
 */

/** creates and captures an LP row from a constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateRowConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler that creates the row */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateRowConshdlr", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, len, cols, vals, lhs, rhs, SCIP_ROWORIGINTYPE_CONSHDLR, (void*) conshdlr, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row from a constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateRowCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   SCIP_CONS*            cons,               /**< constraint that creates the row */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateRowCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, len, cols, vals, lhs, rhs, SCIP_ROWORIGINTYPE_CONS, (void*) cons, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row from a separator
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateRowSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   SCIP_SEPA*            sepa,               /**< separator that creates the row */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(sepa != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateRowSepa", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, len, cols, vals, lhs, rhs, SCIP_ROWORIGINTYPE_SEPA, (void*) sepa, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row from an unspecified source
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateRowUnspec(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateRowUnspec", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, len, cols, vals, lhs, rhs, SCIP_ROWORIGINTYPE_UNSPEC, NULL, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @deprecated Please use SCIPcreateRowConshdlr() or SCIPcreateRowSepa() when calling from a constraint handler or separator in order
 *              to facilitate correct statistics. If the call is from neither a constraint handler or separator, use SCIPcreateRowUnspec().
 */
SCIP_RETCODE SCIPcreateRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcreateRowUnspec(scip, row, name, len, cols, vals, lhs, rhs, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row without any coefficients from a constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateEmptyRowConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler that creates the row */
   const char*           name,               /**< name of row */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRowConshdlr", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, 0, NULL, NULL, lhs, rhs, SCIP_ROWORIGINTYPE_CONSHDLR, (void*) conshdlr, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row without any coefficients from a constraint
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateEmptyRowCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   SCIP_CONS*            cons,               /**< constraint that creates the row */
   const char*           name,               /**< name of row */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(cons != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRowCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, 0, NULL, NULL, lhs, rhs, SCIP_ROWORIGINTYPE_CONS, (void*) cons, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row without any coefficients from a separator
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateEmptyRowSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   SCIP_SEPA*            sepa,               /**< separator that creates the row */
   const char*           name,               /**< name of row */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRowSepa", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, 0, NULL, NULL, lhs, rhs, SCIP_ROWORIGINTYPE_SEPA, (void*) sepa, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row without any coefficients from an unspecified source
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateEmptyRowUnspec(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   const char*           name,               /**< name of row */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRowUnspec", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, 0, NULL, NULL, lhs, rhs, SCIP_ROWORIGINTYPE_UNSPEC, NULL, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row without any coefficients
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @deprecated Please use SCIPcreateEmptyRowConshdlr() or SCIPcreateEmptyRowSepa() when calling from a constraint handler or separator in order
 *              to facilitate correct statistics. If the call is from neither a constraint handler or separator, use SCIPcreateEmptyRowUnspec().
 */
SCIP_RETCODE SCIPcreateEmptyRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   const char*           name,               /**< name of row */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, row, name, lhs, rhs, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** increases usage counter of LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcaptureRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row to capture */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcaptureRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIProwCapture(row);

   return SCIP_OKAY;
}

/** decreases usage counter of LP row, and frees memory if necessary
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPreleaseRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row                 /**< pointer to LP row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPreleaseRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIProwRelease(row, scip->mem->probmem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** changes left hand side of LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgRowLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgRowLhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(!SCIPlpDiving(scip->lp) || (row->lppos == -1));

   SCIP_CALL( SCIProwChgLhs(row, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgRowRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgRowRhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(!SCIPlpDiving(scip->lp) || (row->lppos == -1));

   SCIP_CALL( SCIProwChgRhs(row, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, rhs) );

   return SCIP_OKAY;
}

/** informs row, that all subsequent additions of variables to the row should be cached and not directly applied;
 *  after all additions were applied, SCIPflushRowExtensions() must be called;
 *  while the caching of row extensions is activated, information methods of the row give invalid results;
 *  caching should be used, if a row is build with SCIPaddVarToRow() calls variable by variable to increase
 *  the performance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcacheRowExtensions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcacheRowExtensions", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   return SCIP_OKAY;
}

/** flushes all cached row extensions after a call of SCIPcacheRowExtensions() and merges coefficients with
 *  equal columns into a single coefficient
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPflushRowExtensions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPflushRowExtensions", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* force the row sorting, and merge equal column entries */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** resolves variable to columns and adds them with the coefficient to the row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @attention If the absolute value of val is below the SCIP epsilon tolerance, the variable will not added.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note In case calling this method in the enforcement process of an lp solution, it might be that some variables,
 *        that were not yet in the LP (e.g. dynamic columns) will change their lp solution value returned by SCIP.
 *        For example, a variable, which has a negative objective value, that has no column in the lp yet, is in the lp solution
 *        on its upper bound (variables with status SCIP_VARSTATUS_LOOSE are in an lp solution on it's best bound), but
 *        creating the column, changes the solution value (variable than has status SCIP_VARSTATUS_COLUMN, and the
 *        initialization sets the lp solution value) to 0.0. (This leads to the conclusion that, if a constraint was
 *        violated, the linear relaxation might not be violated anymore.)
 *
 *  @note If the variable being added is FIXED (as given by the status SCIP_VARSTATUS_FIXED), then the variable is not
 *        added to the row, but the corresponding constant is added. Similarly, if the input variable is aggregated (as
 *        given by the status SCIP_VARSTATUS_AGGREGATED), then the input variable is substituted with its aggregation.
 *        For other cases, and to better understand the function behavior, please check the code of SCIPvarAddToRow.
 */
SCIP_RETCODE SCIPaddVarToRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarToRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarAddToRow(var, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lp, row, val) );

   return SCIP_OKAY;
}

/** resolves variables to columns and adds them with the coefficients to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @attention If a coefficients absolute value is below the SCIP epsilon tolerance, the variable with its value is not added.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarsToRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real*            vals                /**< values of coefficients */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarsToRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* resize the row to be able to store all variables (at least, if they are COLUMN variables) */
   SCIP_CALL( SCIProwEnsureSize(row, scip->mem->probmem, scip->set, SCIProwGetNNonz(row) + nvars) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarAddToRow(vars[v], scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lp,
            row, vals[v]) );
   }

   /* force the row sorting */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** resolves variables to columns and adds them with the same single coefficient to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @attention If the absolute value of val is below the SCIP epsilon tolerance, the variables will not added.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarsToRowSameCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real             val                 /**< unique value of all coefficients */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarsToRowSameCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* resize the row to be able to store all variables (at least, if they are COLUMN variables) */
   SCIP_CALL( SCIProwEnsureSize(row, scip->mem->probmem, scip->set, SCIProwGetNNonz(row) + nvars) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarAddToRow(vars[v], scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lp,
            row, val) );
   }

   /* force the row sorting */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** tries to find a value, such that all row coefficients, if scaled with this value become integral
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcalcRowIntegralScalar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcalcRowIntegralScalar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCalcIntegralScalar(row, scip->set, mindelta, maxdelta, maxdnom, maxscale,
         usecontvars, intscalar, success) );

   return SCIP_OKAY;
}

/** tries to scale row, s.t. all coefficients (of integer variables) become integral
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPmakeRowIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal value to scale row with */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Bool*            success             /**< stores whether row could be made rational */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPmakeRowIntegral", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwMakeIntegral(row, scip->mem->probmem, scip->set, scip->eventqueue, scip->stat, scip->lp, mindelta, maxdelta, maxdnom, maxscale,
         usecontvars, success) );

   return SCIP_OKAY;
}

/** marks a row to be not removable from the LP in the current node
 *
 *  @pre this method can be called in the following stage of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void SCIPmarkRowNotRemovableLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPmarkRowNotRemovableLocal", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIProwMarkNotRemovableLocal(row, scip->stat);
}

/** returns number of integral columns in the row
 *
 *  @return number of integral columns in the row
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetRowNumIntCols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowNumIntCols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetNumIntCols(row, scip->set);
}

/** returns minimal absolute value of row vector's non-zero coefficients
 *
 *  @return minimal absolute value of row vector's non-zero coefficients
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowMinCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowMinCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetMinval(row, scip->set);
}

/** returns maximal absolute value of row vector's non-zero coefficients
 *
 *  @return maximal absolute value of row vector's non-zero coefficients
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowMaxCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowMaxCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetMaxval(row, scip->set);
}

/** returns the minimal activity of a row w.r.t. the column's bounds
 *
 *  @return the minimal activity of a row w.r.t. the column's bounds
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowMinActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowMinActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetMinActivity(row, scip->set, scip->stat);
}

/** returns the maximal activity of a row w.r.t. the column's bounds
 *
 *  @return the maximal activity of a row w.r.t. the column's bounds
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowMaxActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowMaxActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetMaxActivity(row, scip->set, scip->stat);
}

/** recalculates the activity of a row in the last LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPrecalcRowLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPrecalcRowLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIProwRecalcLPActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row in the last LP solution
 *
 *  @return activity of a row in the last LP solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetLPActivity(row, scip->set, scip->stat, scip->lp);
}

/** returns the feasibility of a row in the last LP solution
 *
 *  @return the feasibility of a row in the last LP solution: negative value means infeasibility
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowLPFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowLPFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetLPFeasibility(row, scip->set, scip->stat, scip->lp);
}

/** recalculates the activity of a row for the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPrecalcRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPrecalcRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIProwRecalcPseudoActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row for the current pseudo solution
 *
 *  @return the activity of a row for the current pseudo solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetPseudoActivity(row, scip->set, scip->stat);
}

/** returns the feasibility of a row for the current pseudo solution: negative value means infeasibility
 *
 *  @return the feasibility of a row for the current pseudo solution: negative value means infeasibility
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowPseudoFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowPseudoFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetPseudoFeasibility(row, scip->set, scip->stat);
}

/** recalculates the activity of a row in the last LP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPrecalcRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPrecalcRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      SCIProwRecalcLPActivity(row, scip->stat);
   else
      SCIProwRecalcPseudoActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row in the last LP or pseudo solution
 *
 *  @return the activity of a row in the last LP or pseudo solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPActivity(row, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoActivity(row, scip->set, scip->stat);
}

/** returns the feasibility of a row in the last LP or pseudo solution
 *
 *  @return the feasibility of a row in the last LP or pseudo solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPFeasibility(row, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoFeasibility(row, scip->set, scip->stat);
}

/** returns the activity of a row for the given primal solution
 *
 *  @return the activitiy of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowSolActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowSolActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
      return SCIProwGetSolActivity(row, scip->set, scip->stat, sol);
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPActivity(row, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoActivity(row, scip->set, scip->stat);
}

/** returns the feasibility of a row for the given primal solution
 *
 *  @return the feasibility of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowSolFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowSolFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
      return SCIProwGetSolFeasibility(row, scip->set, scip->stat, sol);
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPFeasibility(row, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoFeasibility(row, scip->set, scip->stat);
}

/** returns the parallelism of row with objective function
 *
 *  @return 1 is returned if the row is parallel to the objective function and 0 if it is orthogonal
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRowObjParallelism(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowObjParallelism", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIProwGetObjParallelism(row, scip->set, scip->lp);
}

/** output row to file stream via the message handler system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPprintRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(row != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIProwPrint(row, scip->messagehdlr, file);

   return SCIP_OKAY;
}

/** initiates LP diving, making methods SCIPchgVarObjDive(), SCIPchgVarLbDive(), and SCIPchgVarUbDive() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note diving is allowed even if the current LP is not flushed, not solved, or not solved to optimality; be aware
 *  that solving the (first) diving LP may take longer than expect and that the latter two cases could stem from
 *  numerical troubles during the last LP solve; because of this, most users will want to call this method only if
 *  SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL
 */
SCIP_RETCODE SCIPstartDive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPstartDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   assert(SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_FOCUSNODE);

   if( SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("already in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("cannot start diving while being in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   if( !SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      SCIPerrorMessage("cannot start diving if LP has not been constructed\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPtreeHasCurrentNodeLP(scip->tree));

   SCIP_CALL( SCIPlpStartDive(scip->lp, scip->mem->probmem, scip->set, scip->stat) );

   /* remember the relaxation solution to reset it later */
   if( SCIPisRelaxSolValid(scip) )
   {
      SCIP_CALL( SCIPtreeStoreRelaxSol(scip->tree, scip->set, scip->relaxation, scip->transprob) );
   }

   return SCIP_OKAY;
}

/** quits LP diving and resets bounds and objective values of columns to the current node's values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPendDive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPendDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* unmark the diving flag in the LP and reset all variables' objective and bound values */
   SCIP_CALL( SCIPlpEndDive(scip->lp, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->eventqueue, scip->eventfilter,
         scip->transprob, scip->transprob->vars, scip->transprob->nvars) );

   /* the lower bound may have changed slightly due to LP resolve in SCIPlpEndDive() */
   if( !scip->lp->resolvelperror && scip->tree->focusnode != NULL && SCIPlpIsRelax(scip->lp) && SCIPlpIsSolved(scip->lp) )
   {
      assert(SCIPtreeIsFocusNodeLPConstructed(scip->tree));
      SCIP_CALL( SCIPnodeUpdateLowerboundLP(scip->tree->focusnode, scip->set, scip->stat, scip->tree, scip->transprob,
            scip->origprob, scip->lp) );
   }
   /* reset the probably changed LP's cutoff bound */
   SCIP_CALL( SCIPlpSetCutoffbound(scip->lp, scip->set, scip->transprob, scip->primal->cutoffbound) );
   assert(scip->lp->cutoffbound == scip->primal->cutoffbound); /*lint !e777*/

   /* if a new best solution was created, the cutoff of the tree was delayed due to diving;
    * the cutoff has to be done now.
    */
   if( scip->tree->cutoffdelayed )
   {
      SCIP_CALL( SCIPtreeCutoff(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
            scip->eventqueue, scip->lp, scip->primal->cutoffbound) );
   }

   /* if a relaxation was stored before diving, restore it now */
   if( scip->tree->probdiverelaxstored )
   {
      SCIP_CALL( SCIPtreeRestoreRelaxSol(scip->tree, scip->set, scip->relaxation, scip->transprob) );
   }

   return SCIP_OKAY;
}

/** changes cutoffbound in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPchgCutoffboundDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newcutoffbound      /**< new cutoffbound */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgCutoffboundDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpSetCutoffbound(scip->lp, scip->set, scip->transprob, newcutoffbound) );

   return SCIP_OKAY;
}

/** changes variable's objective value in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPchgVarObjDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarObjDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* invalidate the LP's cutoff bound, since this has nothing to do with the current objective value anymore;
    * the cutoff bound is reset in SCIPendDive()
    */
   SCIP_CALL( SCIPlpSetCutoffbound(scip->lp, scip->set, scip->transprob, SCIPsetInfinity(scip->set)) );

   /* mark the LP's objective function invalid */
   SCIPlpMarkDivingObjChanged(scip->lp);

   /* change the objective value of the variable in the diving LP */
   SCIP_CALL( SCIPvarChgObjDive(var, scip->set, scip->lp, newobj) );

   return SCIP_OKAY;
}

/** changes variable's lower bound in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPchgVarLbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPvarChgLbDive(var, scip->set, scip->lp, newbound) );

   return SCIP_OKAY;
}

/** changes variable's upper bound in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPchgVarUbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPvarChgUbDive(var, scip->set, scip->lp, newbound) );

   return SCIP_OKAY;
}

/** adds a row to the LP in current dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPaddRowDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row to be added */
   )
{
   SCIP_NODE* node;
   int depth;

   assert(scip != NULL);
   assert(row != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddRowDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* get depth of current node */
   node = SCIPtreeGetCurrentNode(scip->tree);
   assert(node != NULL);
   depth = SCIPnodeGetDepth(node);

   SCIP_CALL( SCIPlpAddRow(scip->lp, scip->mem->probmem, scip->set, scip->eventqueue, scip->eventfilter, row, depth) );

   return SCIP_OKAY;
}

/** changes row lhs in current dive, change will be undone after diving ends, for permanent changes use SCIPchgRowLhs()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPchgRowLhsDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< row to change the lhs for */
   SCIP_Real             newlhs              /**< new value for lhs */
   )
{
   assert(scip != NULL);
   assert(row != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgRowLhsDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpRecordOldRowSideDive(scip->lp, row, SCIP_SIDETYPE_LEFT) );
   SCIP_CALL( SCIProwChgLhs(row, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, newlhs) );

   return SCIP_OKAY;
}

/** changes row rhs in current dive, change will be undone after diving ends, for permanent changes use SCIPchgRowRhs()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPchgRowRhsDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< row to change the lhs for */
   SCIP_Real             newrhs              /**< new value for rhs */
   )
{
   assert(scip != NULL);
   assert(row != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgRowRhsDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpRecordOldRowSideDive(scip->lp, row, SCIP_SIDETYPE_RIGHT) );
   SCIP_CALL( SCIProwChgRhs(row, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, newrhs) );

   return SCIP_OKAY;
}

/** gets variable's objective value in current dive
 *
 *  @return the variable's objective value in current dive.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetVarObjDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarObjDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }

   return SCIPvarGetObjLP(var);
}

/** gets variable's lower bound in current dive
 *
 *  @return the variable's lower bound in current dive.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetVarLbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }

   return SCIPvarGetLbLP(var, scip->set);
}

/** gets variable's upper bound in current dive
 *
 *  @return the variable's upper bound in current dive.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Real SCIPgetVarUbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }

   return SCIPvarGetUbLP(var, scip->set);
}

/** solves the LP of the current dive; no separation or pricing is applied
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note be aware that the LP solve may take longer than expected if SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL,
 *  compare the explanation of SCIPstartDive()
 */
SCIP_RETCODE SCIPsolveDiveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the diving LP was infeasible or the objective
                                              *   limit was reached (or NULL, if not needed) */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveDiveLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( cutoff != NULL )
      *cutoff = FALSE;

   /* solve diving LP */
   SCIP_CALL( SCIPlpSolveAndEval(scip->lp, scip->set, scip->messagehdlr, scip->mem->probmem, scip->stat,
         scip->eventqueue, scip->eventfilter, scip->transprob, (SCIP_Longint)itlim, FALSE, FALSE, FALSE, lperror) );

   /* the LP is infeasible or the objective limit was reached */
   if( SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_OBJLIMIT
      || (SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_OPTIMAL &&
         SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip))) )
   {
      /* analyze the infeasible LP (only if the objective was not changed, all columns are in the LP, and no external
       * pricers exist) */
      if( !scip->set->misc_exactsolve && !(SCIPlpDivingObjChanged(scip->lp) || SCIPlpDivingRowsChanged(scip->lp))
         && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) )
      {
         SCIP_CALL( SCIPconflictAnalyzeLP(scip->conflict, scip->conflictstore, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, NULL) );
      }

      if( cutoff != NULL )
         *cutoff = TRUE;
   }

   return SCIP_OKAY;
}

/** returns the number of the node in the current branch and bound run, where the last LP was solved in diving
 *  or probing mode
 *
 *  @return the number of the node in the current branch and bound run, where the last LP was solved in diving
 *  or probing mode.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Longint SCIPgetLastDivenode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLastDivenode", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->lastdivenode;
}

/** returns whether we are in diving mode
 *
 *  @return whether we are in diving mode.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPinDive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPinDive", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPlpDiving(scip->lp);
}

/** computes two measures for dual degeneracy (dual degeneracy rate and variable-constraint ratio)
 *  based on the changes applied when reducing the problem to the optimal face
 *
 *  returns the dual degeneracy rate, i.e., the share of nonbasic variables with reduced cost 0
 *  and the variable-constraint ratio, i.e., the number of unfixed variables in relation to the basis size
 */
SCIP_RETCODE SCIPgetLPDualDegeneracy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            degeneracy,         /**< pointer to store the dual degeneracy rate */
   SCIP_Real*            varconsratio        /**< pointer to store the variable-constraint ratio */
   )
{
   assert(scip != NULL);
   assert(degeneracy != NULL);
   assert(varconsratio != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLPDualDegeneracy", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpGetDualDegeneracy(scip->lp, scip->set, scip->stat, degeneracy, varconsratio) );

   return SCIP_OKAY;
}
