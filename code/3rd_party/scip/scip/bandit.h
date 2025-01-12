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

/**@file   bandit.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for bandit algorithms
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BANDIT_H__
#define __SCIP_BANDIT_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_primal.h"
#include "scip/type_bandit.h"
#include "scip/stat.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates and resets bandit algorithm */
SCIP_RETCODE SCIPbanditCreate(
   SCIP_BANDIT**         bandit,             /**< pointer to bandit algorithm data structure */
   SCIP_BANDITVTABLE*    banditvtable,       /**< virtual table for this bandit algorithm */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   int                   nactions,           /**< the positive number of actions for this bandit */
   unsigned int          initseed,           /**< initial seed for random number generation */
   SCIP_BANDITDATA*      banditdata          /**< algorithm specific bandit data */
   );

/** calls destructor and frees memory of bandit algorithm */
SCIP_RETCODE SCIPbanditFree(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_BANDIT**         bandit              /**< pointer to bandit algorithm data structure */
   );

/** reset the bandit algorithm */
SCIP_RETCODE SCIPbanditReset(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDIT*          bandit,             /**< pointer to bandit algorithm data structure */
   SCIP_Real*            priorities,         /**< priorities for every action, or NULL if not needed */
   unsigned int          seed                /**< initial random seed for bandit selection */
   );

/** get data of this bandit algorithm */
SCIP_BANDITDATA* SCIPbanditGetData(
   SCIP_BANDIT*          bandit              /**< pointer to bandit algorithm data structure */
   );

/** set the data of this bandit algorithm */
void SCIPbanditSetData(
   SCIP_BANDIT*          bandit,             /**< bandit algorithm data structure */
   SCIP_BANDITDATA*      banditdata          /**< bandit algorihm specific data */
   );

/** create a bandit VTable for bandit algorithm callback functions */
SCIP_RETCODE SCIPbanditvtableCreate(
   SCIP_BANDITVTABLE**   banditvtable,       /**< pointer to virtual table for bandit algorithm */
   const char*           name,               /**< a name for the algorithm represented by this vtable */
   SCIP_DECL_BANDITFREE  ((*banditfree)),    /**< callback to free bandit specific data structures */
   SCIP_DECL_BANDITSELECT((*banditselect)),  /**< selection callback for bandit selector */
   SCIP_DECL_BANDITUPDATE((*banditupdate)),  /**< update callback for bandit algorithms */
   SCIP_DECL_BANDITRESET ((*banditreset))    /**< update callback for bandit algorithms */
   );

/** free a bandit vTable for bandit algorithm callback functions */
void SCIPbanditvtableFree(
   SCIP_BANDITVTABLE**   banditvtable        /**< pointer to virtual table for bandit algorithm */
   );

#ifdef __cplusplus
}
#endif

#endif
