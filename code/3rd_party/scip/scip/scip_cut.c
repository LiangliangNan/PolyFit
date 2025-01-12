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

/**@file   scip_cut.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for cuts and aggregation rows
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

#include "scip/cutpool.h"
#include "scip/debug.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pub_cutpool.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cut.h"
#include "scip/scip_numerics.h"
#include "scip/scip_tree.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/solve.h"
#include "scip/struct_lp.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/tree.h"

/** returns row's cutoff distance in the direction of the given primal solution
 *
 *  @return the cutoff distance of the cut with respect to the LP solution in the direction of the given primal solution
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetCutLPSolCutoffDistance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to compute direction for cutoff distance; must not be NULL */
   SCIP_ROW*             cut                 /**< separated cut */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetCutLPSolCutoffDistance", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(sol != NULL);

   return SCIProwGetLPSolCutoffDistance(cut, scip->set, scip->stat, sol, scip->lp);
}

/** returns efficacy of the cut with respect to the given primal solution or the current LP solution:
 *  e = -feasibility/norm
 *
 *  @return the efficacy of the cut with respect to the given primal solution or the current LP solution:
 *          e = -feasibility/norm
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetCutEfficacy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for current LP solution */
   SCIP_ROW*             cut                 /**< separated cut */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetCutEfficacy", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( sol == NULL )
      return SCIProwGetLPEfficacy(cut, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetSolEfficacy(cut, scip->set, scip->stat, sol);
}

/** returns whether the cut's efficacy with respect to the given primal solution or the current LP solution is greater
 *  than the minimal cut efficacy
 *
 *  @return TRUE if the cut's efficacy with respect to the given primal solution or the current LP solution is greater
 *          than the minimal cut efficacy, otherwise FALSE
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Bool SCIPisCutEfficacious(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for current LP solution */
   SCIP_ROW*             cut                 /**< separated cut */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisCutEfficacious", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( sol == NULL )
      return SCIProwIsLPEfficacious(cut, scip->set, scip->stat, scip->lp, (SCIPtreeGetCurrentDepth(scip->tree) == 0));
   else
      return SCIProwIsSolEfficacious(cut, scip->set, scip->stat, sol, (SCIPtreeGetCurrentDepth(scip->tree) == 0));
}

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy
 *
 *  @return TRUE if the given cut's efficacy is larger than the minimal cut efficacy, otherwise FALSE
 */
SCIP_Bool SCIPisEfficacious(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             efficacy            /**< efficacy of the cut */
   )
{
   assert(scip != NULL);

   return SCIPsetIsEfficacious(scip->set, (SCIPtreeGetCurrentDepth(scip->tree) == 0), efficacy);
}

/** calculates the efficacy norm of the given vector, which depends on the "separating/efficacynorm" parameter
 *
 *  @return the efficacy norm of the given vector, which depends on the "separating/efficacynorm" parameter
 */
SCIP_Real SCIPgetVectorEfficacyNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvals               /**< number of values */
   )
{
   SCIP_Real norm;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   norm = 0.0;
   switch( scip->set->sepa_efficacynorm )
   {
   case 'e':
      for( i = 0; i < nvals; ++i )
         norm += SQR(vals[i]);
      norm = SQRT(norm);
      break;
   case 'm':
      for( i = 0; i < nvals; ++i )
      {
         SCIP_Real absval;

         absval = REALABS(vals[i]);
         norm = MAX(norm, absval);
      }
      break;
   case 's':
      for( i = 0; i < nvals; ++i )
         norm += REALABS(vals[i]);
      break;
   case 'd':
      for( i = 0; i < nvals; ++i )
      {
         if( !SCIPisZero(scip, vals[i]) )
         {
            norm = 1.0;
            break;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", scip->set->sepa_efficacynorm);
      assert(FALSE); /*lint !e506*/
   }

   return norm;
}

/** indicates whether a cut is applicable, i.e., will modify the LP when applied
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @return whether the cut is modifiable, not a bound change, or a bound change that changes bounds by at least epsilon
 */
SCIP_Bool SCIPisCutApplicable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             cut                 /**< separated cut */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisCutApplicable", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPsepastoreIsCutApplicable(scip->set, cut);
}

/** adds cut to separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @deprecated Please use SCIPaddRow() instead, or, if the row is a global cut, add it only to the global cutpool.
 */
SCIP_DEPRECATED
SCIP_RETCODE SCIPaddCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that was separated, or NULL for LP solution */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   SCIP_Bool*            infeasible          /**< pointer to store whether cut has been detected to be infeasible for local bounds */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_UNUSED(sol);

   return SCIPaddRow(scip, cut, forcecut, infeasible);
}

/** adds row to separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< row */
   SCIP_Bool             forcecut,           /**< should the row be forced to enter the LP? */
   SCIP_Bool*            infeasible          /**< pointer to store whether row has been detected to be infeasible for local bounds */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(SCIPtreeGetCurrentNode(scip->tree) != NULL);

   SCIP_CALL( SCIPsepastoreAddCut(scip->sepastore, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
         scip->eventfilter, scip->lp, row, forcecut, (SCIPtreeGetCurrentDepth(scip->tree) == 0), infeasible) );

   /* possibly run conflict analysis */
   if ( *infeasible && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && SCIPisConflictAnalysisApplicable(scip) )
   {
      SCIP_Real act;
      SCIP_VAR* var;
      SCIP_Real val;
      int ncols;
      int j;

      /* initialize conflict analysis */
      SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

      if ( ! SCIPisInfinity(scip, -row->lhs) )
      {
         act = SCIProwGetMaxActivity(row, scip->set, scip->stat);
         if ( SCIPisLT(scip, act, row->lhs) )
         {
            ncols = SCIProwGetNNonz(row);
            for (j = 0; j < ncols; ++j)
            {
               val = row->vals[j];
               if ( ! SCIPisZero(scip, val) )
               {
                  var = SCIPcolGetVar(row->cols[j]);
                  assert( var != NULL );

                  if ( val > 0.0 )
                  {
                     SCIP_CALL( SCIPaddConflictUb(scip, var, NULL) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPaddConflictLb(scip, var, NULL) );
                  }
               }
            }
         }
      }
      else if ( ! SCIPisInfinity(scip, row->rhs) )
      {
         act = SCIProwGetMinActivity(row, scip->set, scip->stat);
         if ( SCIPisGT(scip, act, row->rhs) )
         {
            ncols = SCIProwGetNNonz(row);
            for (j = 0; j < ncols; ++j)
            {
               val = row->vals[j];
               if ( ! SCIPisZero(scip, val) )
               {
                  var = SCIPcolGetVar(row->cols[j]);
                  assert( var != NULL );

                  if ( val > 0.0 )
                  {
                     SCIP_CALL( SCIPaddConflictLb(scip, var, NULL) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPaddConflictUb(scip, var, NULL) );
                  }
               }
            }
         }
      }

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, SCIPgetDepth(scip), NULL) );
   }

   return SCIP_OKAY;
}

/** checks if cut is already existing in global cutpool
 *
 *  @return TRUE is returned if the cut is not already existing in the global cutpool, FALSE otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Bool SCIPisCutNew(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisCutNew", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPcutpoolIsCutNew(scip->cutpool, scip->set, row);
}

/** if not already existing, adds row to global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row to remove */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddPoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolAddRow(scip->cutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** removes the row from the global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPdelPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdelPoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolDelRow(scip->cutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** gets current cuts in the global cut pool
 *
 *  @return the current cuts in the global cut pool
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_CUT** SCIPgetPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPoolCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPcutpoolGetCuts(scip->cutpool);
}

/** gets current number of rows in the global cut pool
 *
 *  @return the current number of rows in the global cut pool
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPoolCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPcutpoolGetNCuts(scip->cutpool);
}

/** gets the global cut pool used by SCIP
 *
 *  @return the global cut pool used by SCIP
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_CUTPOOL* SCIPgetGlobalCutpool(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetGlobalCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->cutpool;
}

/** creates a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
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
 */
SCIP_RETCODE SCIPcreateCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   int                   agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateCutpool", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolCreate(cutpool, scip->mem->probmem, scip->set, agelimit, FALSE) );

   return SCIP_OKAY;
}

/** frees a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
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
 */
SCIP_RETCODE SCIPfreeCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL**        cutpool             /**< pointer to store cut pool */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeCutpool", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPcutpoolFree(cutpool, scip->mem->probmem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** if not already existing, adds row to a cut pool and captures it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolAddRow(cutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** adds row to a cut pool and captures it; doesn't check for multiple cuts
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddNewRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddNewRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolAddNewRow(cutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** removes the LP row from a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPdelRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< row to remove */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdelRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolDelRow(cutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** separates cuts from a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPseparateCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPseparateCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(SCIPtreeGetCurrentNode(scip->tree) != NULL);

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("cannot add cuts, because node LP is not processed\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPcutpoolSeparate(cutpool, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter,
         scip->lp, scip->sepastore, NULL, FALSE, (SCIPtreeGetCurrentDepth(scip->tree) == 0), result) );

   return SCIP_OKAY;
}

/** separates cuts w.r.t. given solution from a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPseparateSolCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_Bool             pretendroot,        /**< should the cut separators be called as if we are at the root node? */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPseparateSolCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(SCIPtreeGetCurrentNode(scip->tree) != NULL);

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("cannot add cuts, because node LP is not processed\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPcutpoolSeparate(cutpool, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter,
         scip->lp, scip->sepastore, sol, FALSE, pretendroot, result) );

   return SCIP_OKAY;
}

/** if not already existing, adds row to delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddDelayedPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddDelayedPoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolAddRow(scip->delayedcutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** removes the row from the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPdelDelayedPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPdelDelayedPoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolDelRow(scip->delayedcutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** gets current cuts in the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_CUT** SCIPgetDelayedPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetDelayedPoolCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPcutpoolGetCuts(scip->delayedcutpool);
}

/** gets current number of rows in the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
int SCIPgetNDelayedPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNDelayedPoolCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPcutpoolGetNCuts(scip->delayedcutpool);
}

/** gets the delayed global cut pool used by SCIP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_CUTPOOL* SCIPgetDelayedGlobalCutpool(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetDelayedGlobalCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->delayedcutpool;
}

/** separates the given primal solution or the current LP solution by calling the separators and constraint handlers'
 *  separation methods;
 *  the generated cuts are stored in the separation storage and can be accessed with the methods SCIPgetCuts() and
 *  SCIPgetNCuts();
 *  after evaluating the cuts, you have to call SCIPclearCuts() in order to remove the cuts from the
 *  separation storage;
 *  it is possible to call SCIPseparateSol() multiple times with different solutions and evaluate the found cuts
 *  afterwards
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPseparateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_Bool             pretendroot,        /**< should the cut separators be called as if we are at the root node? */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             onlydelayed,        /**< should only separators be called that were delayed in the previous round? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   int actdepth;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPseparateSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* get current depth */
   actdepth = (pretendroot ? 0 : SCIPtreeGetCurrentDepth(scip->tree));

   /* apply separation round */
   SCIP_CALL( SCIPseparationRound(scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->eventqueue,
         scip->eventfilter, scip->transprob, scip->primal, scip->tree, scip->lp, scip->sepastore,
         sol, actdepth, allowlocal, onlydelayed, delayed, cutoff) );

   return SCIP_OKAY;
}

/** gets the array of cuts currently stored in the separation storage
 *
 *  @return the array of cuts currently stored in the separation storage
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_ROW** SCIPgetCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPsepastoreGetCuts(scip->sepastore);
}

/** get current number of cuts in the separation storage
 *
 *  @return the current number of cuts in the separation storage
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPsepastoreGetNCuts(scip->sepastore);
}

/** clears the separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPclearCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPclearCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsepastoreClearCuts(scip->sepastore, scip->mem->probmem, scip->set, scip->eventqueue, scip->eventfilter, scip->lp) );

   return SCIP_OKAY;
}

/** removes cuts that are inefficacious w.r.t. the current LP solution from separation storage without adding the cuts to the LP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPremoveInefficaciousCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool isroot = FALSE;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPremoveInefficaciousCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeGetCurrentDepth(scip->tree) == 0 )
      isroot = TRUE;

   SCIP_CALL( SCIPsepastoreRemoveInefficaciousCuts(scip->sepastore, scip->mem->probmem, scip->set, scip->stat,
         scip->eventqueue, scip->eventfilter, scip->lp, isroot, SCIP_EFFICIACYCHOICE_LP) );

   return SCIP_OKAY;
}
