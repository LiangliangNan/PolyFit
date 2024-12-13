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

/**@file   type_bandit.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for bandit selection algorithms
 * @author Gregor Hendel
 *
 *  This file defines the interface for bandit selection algorithms implemented in C.
 *  see \ref PublicBanditMethods for all publicly available bandit methods.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_BANDIT_H__
#define __SCIP_TYPE_BANDIT_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "scip/type_timing.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure for bandit algorithms */
typedef struct SCIP_Bandit SCIP_BANDIT;

/** virtual function table for bandit callbacks */
typedef struct SCIP_BanditVTable SCIP_BANDITVTABLE;

/** data structure for specific bandit algorithm implementation */
typedef struct SCIP_BanditData SCIP_BANDITDATA;

/*
 * callbacks for bandit virtual function table
 */

/** callback to free bandit specific data structures */
#define SCIP_DECL_BANDITFREE(x) SCIP_RETCODE x (  \
   BMS_BLKMEM*           blkmem,                  \
   SCIP_BANDIT*          bandit                   \
)

/** selection callback for bandit selector */
#define SCIP_DECL_BANDITSELECT(x) SCIP_RETCODE x ( \
   SCIP_BANDIT*          bandit,                   \
   int*                  selection                 \
)

/** update callback for bandit algorithms */
#define SCIP_DECL_BANDITUPDATE(x) SCIP_RETCODE x ( \
   SCIP_BANDIT*          bandit,                   \
   int                   selection,                \
   SCIP_Real             score                     \
)

/** reset callback for bandit algorithms */
#define SCIP_DECL_BANDITRESET(x) SCIP_RETCODE x (  \
   BMS_BUFMEM*           bufmem,                   \
   SCIP_BANDIT*          bandit,                   \
   SCIP_Real*            priorities                \
)

#ifdef __cplusplus
}
#endif

#endif
