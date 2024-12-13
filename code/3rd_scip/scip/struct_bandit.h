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

/**@file   struct_bandit.h
 * @ingroup INTERNALAPI
 * @brief  data structures for bandit selection algorithms
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BANDIT_H__
#define __SCIP_STRUCT_BANDIT_H__


#include "scip/def.h"
#include "misc.h"
#include "scip/type_clock.h"
#include "scip/type_bandit.h"

#ifdef __cplusplus
extern "C" {
#endif

/** virtual function table for bandit selection algorithms */
struct SCIP_BanditVTable
{
   const char*           name;               /**< name of the represented bandit algorithm */
   SCIP_DECL_BANDITFREE  ((*banditfree));    /**< callback to free bandit specific data structures */
   SCIP_DECL_BANDITSELECT((*banditselect));  /**< selection callback for bandit selector */
   SCIP_DECL_BANDITUPDATE((*banditupdate));  /**< update callback for bandit algorithms */
   SCIP_DECL_BANDITRESET ((*banditreset));   /**< update callback for bandit algorithms */
};

/** data structure for bandit algorithms */
struct SCIP_Bandit
{
   SCIP_BANDITVTABLE*    vtable;             /**< virtual function table for callbacks */
   SCIP_RANDNUMGEN*      rng;                /**< random number generator for randomized selection */
   int                   nactions;           /**< the number of actions to select from */
   SCIP_BANDITDATA*      data;               /**< specific data for bandit algorithm implementations */
};
#ifdef __cplusplus
}
#endif

#endif
