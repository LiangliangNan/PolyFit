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

/**@file   concurrent.h
 * @ingroup PARALLEL
 * @brief  helper functions for concurrent scip solvers
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/type_concurrent.h"
#include "scip/type_scip.h"
#include "scip/type_concsolver.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"
#include "scip/type_syncstore.h"
#include "scip/def.h"

#ifndef __SCIP_CONCURRENT_H__
#define __SCIP_CONCURRENT_H__

#ifdef __cplusplus
extern "C" {
#endif

/** create concurrent data */
SCIP_RETCODE SCIPcreateConcurrent(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver of given SCIP instance */
   int*                  varperm             /**< permutation of variables for communication */
   );

/** get number of initialized concurrent solvers */
int SCIPgetNConcurrentSolvers(
   SCIP*                 scip                /**< SCIP datastructure */
   );

/** gets the concurrent solvers */
SCIP_CONCSOLVER** SCIPgetConcurrentSolvers(
   SCIP*                 scip                /**< SCIP datastructure */
   );

/** adds a concurrent solver */
SCIP_RETCODE SCIPaddConcurrentSolver(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver of given SCIP instance */
   );

/** frees concurrent data */
SCIP_RETCODE SCIPfreeConcurrent(
   SCIP*                 scip                /**< SCIP datastructure */
   );

/** increments the time counter for synchronization */
SCIP_RETCODE SCIPincrementConcurrentTime(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real             val                 /**< value by which the time counter for synchronization is incremented */
   );

/** synchronize with other concurrent solvers */
SCIP_RETCODE SCIPsynchronize(
   SCIP*                 scip                /**< SCIP datastructure */
   );

/** pass a solution to the given SCIP instance using that was received via synchronization by using
 * the sync heuristic */
SCIP_RETCODE SCIPaddConcurrentSol(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_SOL*             sol                 /**< solution */
   );

/** adds a global boundchange to the given SCIP, by passing it to the sync propagator */
SCIP_RETCODE SCIPaddConcurrentBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for bound */
   SCIP_Real             val,                /**< value of bound */
   SCIP_BOUNDTYPE        bndtype             /**< type of bound */
   );

/** copy the nodenumber, depth, time, and runnumber of one solution to another one */
SCIP_RETCODE SCIPcopySolStats(
   SCIP_SOL*             source,             /**< source for solution statistics */
   SCIP_SOL*             target              /**< target for solution statistics */
   );

/** copy solving statistics */
SCIP_RETCODE SCIPcopyConcurrentSolvingStats(
   SCIP*                 source,             /**< SCIP data structure */
   SCIP*                 target              /**< target SCIP data structure */
   );

/** get variable index of original variable that is the same between concurrent solvers */
int SCIPgetConcurrentVaridx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable */
   );

/** has the solution been created after the last synchronization point */
SCIP_Bool SCIPIsConcurrentSolNew(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< the solution */
   );

/** gets the global bound changes since the last synchronization point */
SCIP_BOUNDSTORE* SCIPgetConcurrentGlobalBoundChanges(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** start solving in parallel using the given set of concurrent solvers */
SCIP_RETCODE SCIPconcurrentSolve(
   SCIP*                 scip                /**< pointer to scip datastructure */
   );

/** disables storing global bound changes */
void SCIPdisableConcurrentBoundStorage(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** enables storing global bound changes */
void SCIPenableConcurrentBoundStorage(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets total memory usage of all concurrent solvers together */
SCIP_Longint SCIPgetConcurrentMemTotal(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the dualbound in the last synchronization */
SCIP_Real SCIPgetConcurrentDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the primalbound in the last synchronization */
SCIP_Real SCIPgetConcurrentPrimalbound(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the gap in the last synchronization */
SCIP_Real SCIPgetConcurrentGap(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gives the total number of tightened bounds received from other concurrent solvers */
SCIP_Longint SCIPgetConcurrentNTightenedBnds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gives the total number of tightened bounds for integer variables received from
 *  other concurrent solvers */
SCIP_Longint SCIPgetConcurrentNTightenedIntBnds(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
