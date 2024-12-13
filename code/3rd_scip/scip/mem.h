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

/**@file   mem.h
 * @ingroup INTERNALAPI
 * @brief  methods for block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MEM_H__
#define __SCIP_MEM_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_mem.h"
#include "scip/struct_mem.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates block and buffer memory structures */
SCIP_RETCODE SCIPmemCreate(
   SCIP_MEM**            mem                 /**< pointer to block and buffer memory structure */
   );

/** frees block and buffer memory structures */
SCIP_RETCODE SCIPmemFree(
   SCIP_MEM**            mem                 /**< pointer to block and buffer memory structure */
   );

/** returns the total number of bytes used in block and buffer memory */
SCIP_Longint SCIPmemGetUsed(
   SCIP_MEM*             mem                 /**< pointer to block and buffer memory structure */
   );

/** returns the total number of bytes in block and buffer memory */
SCIP_Longint SCIPmemGetTotal(
   SCIP_MEM*             mem                 /**< pointer to block and buffer memory structure */
   );

/** returns the maximal number of used bytes in block memory */
SCIP_Longint SCIPmemGetUsedBlockmemoryMax(
   SCIP_MEM*             mem                 /**< pointer to block and buffer memory structure */
   );

/** returns the maximal number of allocated but not used bytes in block memory */
SCIP_Longint SCIPmemGetUnusedBlockmemoryMax(
   SCIP_MEM*             mem                 /**< pointer to block and buffer memory structure */
   );

/** returns the maximal number of allocated bytes in block memory */
SCIP_Longint SCIPmemGetAllocatedBlockmemoryMax(
   SCIP_MEM*             mem                 /**< pointer to block and buffer memory structure */
   );

#ifdef __cplusplus
}
#endif

#endif
