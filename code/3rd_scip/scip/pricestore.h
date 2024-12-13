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

/**@file   pricestore.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing priced variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRICESTORE_H__
#define __SCIP_PRICESTORE_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_pricestore.h"
#include "scip/type_branch.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates pricing storage */
SCIP_RETCODE SCIPpricestoreCreate(
   SCIP_PRICESTORE**     pricestore          /**< pointer to store pricing storage */
   );

/** frees pricing storage */
SCIP_RETCODE SCIPpricestoreFree(
   SCIP_PRICESTORE**     pricestore          /**< pointer to store pricing storage */
   );

/** informs pricing storage, that the setup of the initial LP starts now */
void SCIPpricestoreStartInitialLP(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

/** informs pricing storage, that the setup of the initial LP is now finished */
void SCIPpricestoreEndInitialLP(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

/** adds variable to pricing storage and capture it */
SCIP_RETCODE SCIPpricestoreAddVar(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_VAR*             var,                /**< priced variable */
   SCIP_Real             score,              /**< pricing score of variable (the larger, the better the variable) */
   SCIP_Bool             root                /**< are we at the root node? */
   );

/** adds variable where zero violates the bounds to pricing storage, capture it */
SCIP_RETCODE SCIPpricestoreAddBdviolvar(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var                 /**< variable, where zero violates the bounds */
   );

/** adds problem variables with negative reduced costs to pricing storage */
SCIP_RETCODE SCIPpricestoreAddProbVars(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** adds priced variables to the LP */
SCIP_RETCODE SCIPpricestoreApplyVars(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< LP data */
   );

/** reset variables' bounds violated by zero to its original value */
SCIP_RETCODE SCIPpricestoreResetBounds(
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** gets number of variables in pricing storage */
int SCIPpricestoreGetNVars(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

/** gets number of variables in pricing storage whose bounds must be reset */
int SCIPpricestoreGetNBoundResets(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

/** gets time needed to price existing problem variables */
SCIP_Real SCIPpricestoreGetProbPricingTime(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

/** gets total number of calls to problem variable pricing */
int SCIPpricestoreGetNProbPricings(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

/** gets total number of times, a problem variable was priced in */
int SCIPpricestoreGetNProbvarsFound(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

/** get total number of variables found so far in pricing */
int SCIPpricestoreGetNVarsFound(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

/** get total number of variables priced into the LP so far */
int SCIPpricestoreGetNVarsApplied(
   SCIP_PRICESTORE*      pricestore          /**< pricing storage */
   );

#ifdef __cplusplus
}
#endif

#endif
