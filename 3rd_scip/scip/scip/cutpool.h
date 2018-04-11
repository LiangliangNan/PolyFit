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
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_sepastore.h"
#include "scip/type_cutpool.h"
#include "scip/pub_cutpool.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates cut pool */
extern
SCIP_RETCODE SCIPcutpoolCreate(
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   agelimit,           /**< maximum age a cut can reach before it is deleted from the pool */
   SCIP_Bool             globalcutpool       /**< is this the global cut pool of SCIP? */
   );

/** frees cut pool */
extern
SCIP_RETCODE SCIPcutpoolFree(
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** removes all rows from the cut pool */
extern
SCIP_RETCODE SCIPcutpoolClear(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** checks if cut is already existing */
extern
SCIP_Bool SCIPcutpoolIsCutNew(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** if not already existing, adds row to cut pool and captures it */
extern
SCIP_RETCODE SCIPcutpoolAddRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** adds row to cut pool and captures it; doesn't check for multiple cuts */
extern
SCIP_RETCODE SCIPcutpoolAddNewRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< cutting plane to add */
   );

/** removes the LP row from the cut pool */
extern
SCIP_RETCODE SCIPcutpoolDelRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< row to remove */
   );

/** separates cuts of the cut pool */
extern
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

#ifdef __cplusplus
}
#endif

#endif
