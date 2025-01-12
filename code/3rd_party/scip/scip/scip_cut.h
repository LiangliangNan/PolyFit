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

/**@file   scip_cut.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for cuts and aggregation rows
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_CUT_H__
#define __SCIP_SCIP_CUT_H__


#include "scip/def.h"
#include "scip/type_cutpool.h"
#include "scip/type_lp.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCutMethods
 *
 * @{
 */

/** returns row's cutoff distance in the direction of the given primal solution
 *
 *  @return the cutoff distance of the cut with respect to the LP solution in the direction of the given primal solution
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_Real SCIPgetCutLPSolCutoffDistance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to compute direction for cutoff distance; must not be NULL */
   SCIP_ROW*             cut                 /**< separated cut */
   );

/** returns efficacy of the cut with respect to the given primal solution or the current LP solution:
 *  e = -feasibility/norm
 *
 *  @return the efficacy of the cut with respect to the given primal solution or the current LP solution:
 *          e = -feasibility/norm
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_Real SCIPgetCutEfficacy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for current LP solution */
   SCIP_ROW*             cut                 /**< separated cut */
   );

/** returns whether the cut's efficacy with respect to the given primal solution or the current LP solution is greater
 *  than the minimal cut efficacy
 *
 *  @return TRUE if the cut's efficacy with respect to the given primal solution or the current LP solution is greater
 *          than the minimal cut efficacy, otherwise FALSE
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_Bool SCIPisCutEfficacious(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for current LP solution */
   SCIP_ROW*             cut                 /**< separated cut */
   );

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy
 *
 *  @return TRUE if the given cut's efficacy is larger than the minimal cut efficacy, otherwise FALSE
 */
SCIP_EXPORT
SCIP_Bool SCIPisEfficacious(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             efficacy            /**< efficacy of the cut */
   );

/** calculates the efficacy norm of the given vector, which depends on the "separating/efficacynorm" parameter
 *
 *  @return the efficacy norm of the given vector, which depends on the "separating/efficacynorm" parameter
 */
SCIP_EXPORT
SCIP_Real SCIPgetVectorEfficacyNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvals               /**< number of values */
   );

/** indicates whether a cut is applicable
 *
 *  If the cut has only one variable and this method returns FALSE, it may
 *  still be possible that the cut can be added to the LP (as a row instead
 *  of a boundchange), but it will be a very weak cut. The user is asked
 *  to avoid such cuts.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @return whether the cut is modifiable, not a bound change, or a bound change that changes bounds by at least epsilon
 */
SCIP_EXPORT
SCIP_Bool SCIPisCutApplicable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             cut                 /**< separated cut */
   );

/** adds cut to separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @deprecated Please use SCIPaddRow() instead, or, if the row is a global cut and it might be useful to keep it for future use,
 *  consider adding it to the global cutpool with SCIPaddPoolCut().
 */
SCIP_DEPRECATED SCIP_EXPORT
SCIP_RETCODE SCIPaddCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that was separated, or NULL for LP solution */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   SCIP_Bool*            infeasible          /**< pointer to store whether cut has been detected to be infeasible for local bounds */
   );

/** adds row to separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< row */
   SCIP_Bool             forcecut,           /**< should the row be forced to enter the LP? */
   SCIP_Bool*            infeasible          /**< pointer to store whether row has been detected to be infeasible for local bounds */
   );

/** checks if cut is already existing in global cutpool
 *
 *  @return TRUE is returned if the cut is not already existing in the global cutpool, FALSE otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_Bool SCIPisCutNew(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** if not already existing, adds row to global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** removes the row from the global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdelPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row to remove */
   );

/** gets current cuts in the global cut pool
 *
 *  @return the current cuts in the global cut pool
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_EXPORT
SCIP_CUT** SCIPgetPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current number of rows in the global cut pool
 *
 *  @return the current number of rows in the global cut pool
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_EXPORT
int SCIPgetNPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the global cut pool used by SCIP
 *
 *  @return the global cut pool used by SCIP
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_EXPORT
SCIP_CUTPOOL* SCIPgetGlobalCutpool(
   SCIP*                 scip                /**< SCIP data structure */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPcreateCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   int                   agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPfreeCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL**        cutpool             /**< pointer to store cut pool */
   );

/** if not already existing, adds row to a cut pool and captures it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** adds row to a cut pool and captures it; doesn't check for multiple cuts
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddNewRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPdelRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< row to remove */
   );

/** separates cuts from a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPseparateCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   );

/** separates cuts w.r.t. given solution from a cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPseparateSolCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_Bool             pretendroot,        /**< should the cut separators be called as if we are at the root node? */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   );

/** if not already existing, adds row to the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddDelayedPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** removes the row from the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdelDelayedPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** gets current cuts in the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_CUT** SCIPgetDelayedPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current number of rows in the delayed global cut pool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
int SCIPgetNDelayedPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the delayed global cut pool used by SCIP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is the stages \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_CUTPOOL* SCIPgetDelayedGlobalCutpool(
   SCIP*                 scip                /**< SCIP data structure */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPseparateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_Bool             pretendroot,        /**< should the cut separators be called as if we are at the root node? */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             onlydelayed,        /**< should only separators be called that were delayed in the previous round? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   );

/** gets the array of cuts currently stored in the separation storage
 *
 *  @return the array of cuts currently stored in the separation storage
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
SCIP_ROW** SCIPgetCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** get current number of cuts in the separation storage
 *
 *  @return the current number of cuts in the separation storage
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
int SCIPgetNCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** clears the separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPclearCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** removes cuts that are inefficacious w.r.t. the current LP solution from separation storage without adding the cuts to the LP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPremoveInefficaciousCuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
