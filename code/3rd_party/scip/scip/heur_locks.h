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

/**@file   heur_locks.h
 * @ingroup PRIMALHEURISTICS
 * @brief  locks primal heuristic
 * @author Michael Winkler
 * @author Gerald Gamrath
 *
 * The locks heuristic is a start heuristic that first tries to fix all binary variables, then solves the resulting LP
 * and tries to round the solution and finally solves a sub-MIP on the remaining problem if the LP solution could not be
 * rounded. The fixing works as follows: First, all variables are sorted by their total number of rounding locks (up-
 * and down-locks summed up). Then, looking at the variable with the highest number of locks first, the variable is
 * fixed to the bound where there are fewer locks (in case of ties, the bound which is better w.r.t. the objective
 * function). This fix is propagated and the activities of all LP rows are updated. If any LP row becomes redundant
 * w.r.t. the updated bounds, we adjust the rounding locks.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LOCKS_H__
#define __SCIP_HEUR_LOCKS_H__

#include "scip/def.h"
#include "scip/type_heur.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the locks primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurLocks(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** apply fix-and-propagate scheme based on variable locks
 *
 *  @note probing mode of SCIP needs to be enabled before
 */
SCIP_EXPORT
SCIP_RETCODE SCIPapplyLockFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            cutoff,             /**< pointer to store if a cutoff was detected */
   SCIP_Bool*            allrowsfulfilled    /**< pointer to store if all rows became redundant */
   );

#ifdef __cplusplus
}
#endif

#endif
