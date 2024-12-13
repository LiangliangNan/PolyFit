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

/**@file   heur_sync.h
 * @ingroup PRIMALHEURISTICS
 * @brief  primal heuristic that adds given solutions
 * @author Leona Gottwald
 *
 * This heuristic takes solutions from somewhere else via the function SCIPheurSyncPassSol(). It
 * then tries to commit this solution. It is used by the concurrent solvers, when solutions are
 * communicated between solvers, but cannot directly submitted because SCIP might be in a stage where
 * this is not allowed.
 * If multiple solutions are passed it will keep the best N solutions depending on the parameter setting
 * "concsolvers/sync/maxnsols"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_SYNC_H__
#define __SCIP_HEUR_SYNC_H__

#include "scip/def.h"
#include "scip/type_sol.h"
#include "scip/type_scip.h"
#include "scip/type_heur.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the sync primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurSync(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
  *
  * @{
  */

/** pass solution to sync heuristic */
SCIP_EXPORT
SCIP_RETCODE SCIPheurSyncPassSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< sync heuristic */
   SCIP_SOL*             sol                 /**< solution to be passed */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
