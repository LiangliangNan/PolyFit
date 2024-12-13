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

/**@file   struct_mem.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_MEM_H__
#define __SCIP_STRUCT_MEM_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_mem.h"

#ifdef __cplusplus
extern "C" {
#endif

/** various block memory buffers */
struct SCIP_Mem
{
   BMS_BLKMEM*           setmem;             /**< memory blocks for parameter settings */
   BMS_BLKMEM*           probmem;            /**< memory blocks for original problem and solution process: preprocessing, bab-tree, ... */
   BMS_BUFMEM*           buffer;             /**< memory buffers for short living temporary objects */
   BMS_BUFMEM*           cleanbuffer;        /**< memory buffers for short living temporary objects, initialized to all zero */
};

#ifdef __cplusplus
}
#endif

#endif
