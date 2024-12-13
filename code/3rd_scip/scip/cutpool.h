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

/**@file   cutpool.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTPOOL_H__
#define __SCIP_CUTPOOL_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_event.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_sol.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_sepastore.h"
#include "scip/type_cutpool.h"
#include "scip/pub_cutpool.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates cut pool */
SCIP_RETCODE SCIPcutpoolCreate(
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   agelimit,           /**< maximum age a cut can reach before it is deleted from the pool */
   SCIP_Bool             globalcutpool       /**< is this the global cut pool of SCIP? */
   );

/** frees cut pool */
SCIP_RETCODE SCIPcutpoolFree(
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** removes all rows from the cut pool */
SCIP_RETCODE SCIPcutpoolClear(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** checks if cut is already existing */
SCIP_Bool SCIPcutpoolIsCutNew(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** if not already existing, adds row to cut pool and captures it */
SCIP_RETCODE SCIPcutpoolAddRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** adds row to cut pool and captures it; doesn't check for multiple cuts */
SCIP_RETCODE SCIPcutpoolAddNewRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** removes the LP row from the cut pool */
SCIP_RETCODE SCIPcutpoolDelRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< row to remove */
   );

/** separates cuts of the cut pool */
SCIP_RETCODE SCIPcutpoolSeparate(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL for LP-solution) */
   SCIP_Bool             cutpoolisdelayed,   /**< is the cutpool delayed (count cuts found)? */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   );

/** adds the maximum number of cuts that were stored in the pool;
 *  this is primarily used to keep statistics when SCIP performs a restart */
void SCIPcutpoolAddMaxNCuts(
   SCIP_CUTPOOL*         cutpool,             /**< cut pool */
   SCIP_Longint          ncuts             /**< number of cuts to add */
   );

/** sets time in seconds used for separating cuts from the pool;
 *  this is primarily used to keep statistics when SCIP performs a restart */
void SCIPcutpoolSetTime(
   SCIP_CUTPOOL*         cutpool,             /**< cut pool */
   SCIP_Real             time                 /**< poolclock time */
   );

/** adds the number of times the cut pool was separated;
 *  this is primarily used to keep statistics when SCIP performs a restart */
void SCIPcutpoolAddNCalls(
   SCIP_CUTPOOL*         cutpool,             /**< cut pool */
   SCIP_Longint          ncalls               /**< ncalls */
   );

/** adds the number of times the cut pool was separated at the root;
 *  this is primarily used to keep statistics when SCIP performs a restart */
void SCIPcutpoolAddNRootCalls(
   SCIP_CUTPOOL*         cutpool,             /**< cut pool */
   SCIP_Longint          nrootcalls           /**< nrootcalls */
);

/** adds the total number of cuts that were added to the pool;
 *  this is primarily used to keep statistics when SCIP performs a restart */
void SCIPcutpoolAddNCutsFound(
   SCIP_CUTPOOL*         cutpool,             /**< cut pool */
   SCIP_Longint          ncutsfound           /**< total number of cuts added to cut pool */
   );

/** adds the total number of cuts that were separated from the pool;
 *  this is primarily used to keep statistics when SCIP performs a restart */
void SCIPcutpoolAddNCutsAdded(
   SCIP_CUTPOOL*         cutpool,             /**< cut pool */
   SCIP_Longint          ncutsadded           /**< total number of cuts added from cut pool to sepastore */
   );

#ifdef __cplusplus
}
#endif

#endif
