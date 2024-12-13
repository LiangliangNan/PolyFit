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

/**@file   conflictstore.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing conflicts
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONFLICTSTORE_H__
#define __SCIP_CONFLICTSTORE_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_conflictstore.h"
#include "scip/type_retcode.h"
#include "scip/type_cons.h"
#include "scip/type_event.h"
#include "scip/type_conflict.h"
#include "scip/type_prob.h"
#include "scip/type_reopt.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates separation storage */
SCIP_RETCODE SCIPconflictstoreCreate(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict store */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** frees separation storage */
SCIP_RETCODE SCIPconflictstoreFree(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** clears conflict store */
SCIP_RETCODE SCIPconflictstoreClear(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** cleans up conflict store */
SCIP_RETCODE SCIPconflictstoreClean(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** adds a constraint to the pool of proof constraints based on dual rays
 *
 *  @note this methods captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddDualraycons(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_CONS*            dualproof,          /**< constraint based on a dual ray */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_Bool             hasrelaxvar         /**< does the dual proof contain at least one variable that exists in
                                               *  the current relaxation only? */
   );

/** adds a constraint to the pool of proof constraints based on dual solutions
 *
 *  @note this methods captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddDualsolcons(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_CONS*            dualproof,          /**< constraint based on a dual solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_Real             scale,              /**< scaling factor that needs to be considered when updating the side */
   SCIP_Bool             updateside,         /**< should the side be updated if a new incumbent is found */
   SCIP_Bool             hasrelaxvar         /**< does the dual proof contain at least one variable that exists in
                                               *  the current relaxation only? */
   );

/** adds a conflict to the conflict store
 *
 *  @note this method captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddConflict(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree (or NULL for an original constraint) */
   SCIP_PROB*            transprob,          /**< transformed problem (or NULL for an original constraint) */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_CONS*            cons,               /**< constraint representing the conflict */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             cutoffinvolved,     /**< is a cutoff bound involved in this conflict */
   SCIP_Real             primalbound         /**< primal bound the conflict depend on (or -SCIPinfinity) */
   );

/** deletes all conflicts depending on a cutoff bound larger than the given bound */
SCIP_RETCODE SCIPconflictstoreCleanNewIncumbent(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem*/
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_Real             cutoffbound         /**< current cutoff bound */
   );

/** returns the maximal size of the conflict pool */
int SCIPconflictstoreGetMaxPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   );

/** returns the initial size of the conflict pool */
int SCIPconflictstoreGetInitPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   );

/** returns the number of stored conflicts on the conflict pool
 *
 *  @note the number of active conflicts can be less
 */
int SCIPconflictstoreGetNConflictsInStore(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   );

/** returns all active conflicts stored in the conflict store */
SCIP_RETCODE SCIPconflictstoreGetConflicts(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_CONS**           conflicts,          /**< array to store conflicts */
   int                   conflictsize,       /**< size of the conflict array */
   int*                  nconflicts          /**< pointer to store the number of conflicts */
   );

/** transforms all original conflicts into transformed conflicts */
SCIP_RETCODE SCIPconflictstoreTransform(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** returns the average number of non-zeros over all stored dual ray constraints */
SCIP_Real SCIPconflictstoreGetAvgNnzDualInfProofs(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   );

/** return the number of stored dualray constraints */
int SCIPconflictstoreGetNDualInfProofs(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   );

/** returns the average number of non-zeros over all stored boundexceeding proofs */
SCIP_Real SCIPconflictstoreGetAvgNnzDualBndProofs(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   );

/** returns the number of all stored boundexceeding proofs */
int SCIPconflictstoreGetNDualBndProofs(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   );

#ifdef __cplusplus
}
#endif

#endif
