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

/**@file   scip_mem.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for memory management
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_MEM_H__
#define __SCIP_SCIP_MEM_H__


#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicMemoryMethods
 *
 * @{
 */

/* Standard Memory Management Macros */

#define SCIPallocMemory(scip,ptr)               ( (BMSallocMemory((ptr)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocClearMemory(scip,ptr)          ( (BMSallocClearMemory((ptr)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocMemoryArray(scip,ptr,num)      ( (BMSallocMemoryArray((ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocClearMemoryArray(scip,ptr,num) ( (BMSallocClearMemoryArray((ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocMemorySize(scip,ptr,size)      ( (BMSallocMemorySize((ptr), (size)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocMemoryArray(scip,ptr,newnum) ( (BMSreallocMemoryArray((ptr), (newnum)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocMemorySize(scip,ptr,newsize) ( (BMSreallocMemorySize((ptr), (newsize)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateMemory(scip, ptr, source)  ( (BMSduplicateMemory((ptr), (source)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateMemoryArray(scip, ptr, source, num) ( (BMSduplicateMemoryArray((ptr), (source), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPfreeMemory(scip,ptr)                BMSfreeMemory(ptr)
#define SCIPfreeMemoryNull(scip,ptr)            BMSfreeMemoryNull(ptr)
#define SCIPfreeMemoryArray(scip,ptr)           BMSfreeMemoryArray(ptr)
#define SCIPfreeMemoryArrayNull(scip,ptr)       BMSfreeMemoryArrayNull(ptr)
#define SCIPfreeMemorySize(scip,ptr)            BMSfreeMemorySize(ptr)
#define SCIPfreeMemorySizeNull(scip,ptr)        BMSfreeMemorySizeNull(ptr)

/* Block Memory Management Macros
 *
 */

#define SCIPallocBlockMemory(scip,ptr)          ( (BMSallocBlockMemory(SCIPblkmem(scip), (ptr)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocClearBlockMemory(scip,ptr)     ( (BMSallocClearBlockMemory(SCIPblkmem(scip), (ptr)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocBlockMemoryArray(scip,ptr,num) ( (BMSallocBlockMemoryArray(SCIPblkmem(scip), (ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocBlockMemorySize(scip,ptr,size) ( (BMSallocBlockMemorySize(SCIPblkmem(scip), (ptr), (size)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocClearBlockMemoryArray(scip,ptr,num) ( (BMSallocClearBlockMemoryArray(SCIPblkmem(scip), (ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocBlockMemoryArray(scip,ptr,oldnum,newnum) ( (BMSreallocBlockMemoryArray(SCIPblkmem(scip), (ptr), (oldnum), (newnum)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocBlockMemorySize(scip,ptr,oldsize,newsize) ( (BMSreallocBlockMemorySize(SCIPblkmem(scip), (ptr), (oldsize), (newsize)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateBlockMemory(scip, ptr, source) ( (BMSduplicateBlockMemory(SCIPblkmem(scip), (ptr), (source)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateBlockMemoryArray(scip, ptr, source, num) ( (BMSduplicateBlockMemoryArray(SCIPblkmem(scip), (ptr), (source), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPensureBlockMemoryArray(scip,ptr,arraysizeptr,minsize) ( (SCIPensureBlockMemoryArray_call((scip), (void**)(ptr), sizeof(**(ptr)), (arraysizeptr), (minsize))) )
#define SCIPfreeBlockMemory(scip,ptr)           BMSfreeBlockMemory(SCIPblkmem(scip), (ptr))
#define SCIPfreeBlockMemoryNull(scip,ptr)       BMSfreeBlockMemoryNull(SCIPblkmem(scip), (ptr))
#define SCIPfreeBlockMemoryArray(scip,ptr,num)  BMSfreeBlockMemoryArray(SCIPblkmem(scip), (ptr), (num))
#define SCIPfreeBlockMemoryArrayNull(scip,ptr,num) BMSfreeBlockMemoryArrayNull(SCIPblkmem(scip), (ptr), (num))
#define SCIPfreeBlockMemorySize(scip,ptr,size)  BMSfreeBlockMemorySize(SCIPblkmem(scip), (ptr), (size))
#define SCIPfreeBlockMemorySizeNull(scip,ptr,size) BMSfreeBlockMemorySizeNull(SCIPblkmem(scip), (ptr), (size))


/* Buffer Memory Management Macros
 *
 *
 */


#define SCIPallocBuffer(scip,ptr)               ( (BMSallocBufferMemory(SCIPbuffer(scip), (ptr)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocBufferArray(scip,ptr,num)      ( (BMSallocBufferMemoryArray(SCIPbuffer(scip), (ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocClearBufferArray(scip,ptr,num) ( (BMSallocClearBufferMemoryArray(SCIPbuffer(scip), (ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPreallocBufferArray(scip,ptr,num)    ( (BMSreallocBufferMemoryArray(SCIPbuffer(scip), (ptr), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateBuffer(scip,ptr,source)    ( (BMSduplicateBufferMemory(SCIPbuffer(scip), (ptr), (source), (size_t)sizeof(**(ptr))) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPduplicateBufferArray(scip,ptr,source,num) ( (BMSduplicateBufferMemoryArray(SCIPbuffer(scip), (ptr), (source), (num)) == NULL) \
                                                       ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPfreeBuffer(scip,ptr)                BMSfreeBufferMemorySize(SCIPbuffer(scip), (ptr))
#define SCIPfreeBufferNull(scip,ptr)            BMSfreeBufferMemoryNull(SCIPbuffer(scip), (ptr))
#define SCIPfreeBufferArray(scip,ptr)           BMSfreeBufferMemoryArray(SCIPbuffer(scip), (ptr))
#define SCIPfreeBufferArrayNull(scip,ptr)       BMSfreeBufferMemoryArrayNull(SCIPbuffer(scip), (ptr))


#define SCIPallocCleanBuffer(scip,ptr)          ( (BMSallocBufferMemory(SCIPcleanbuffer(scip), (ptr)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPallocCleanBufferArray(scip,ptr,num) ( (BMSallocBufferMemoryArray(SCIPcleanbuffer(scip), (ptr), (num)) == NULL) \
                                                  ? SCIP_NOMEMORY : SCIP_OKAY )
#define SCIPfreeCleanBuffer(scip,ptr)           BMSfreeBufferMemorySize(SCIPcleanbuffer(scip), (ptr))
#define SCIPfreeCleanBufferNull(scip,ptr)       BMSfreeBufferMemoryNull(SCIPcleanbuffer(scip), (ptr))
#define SCIPfreeCleanBufferArray(scip,ptr)      BMSfreeBufferMemoryArray(SCIPcleanbuffer(scip), (ptr))
#define SCIPfreeCleanBufferArrayNull(scip,ptr)  BMSfreeBufferMemoryArrayNull(SCIPcleanbuffer(scip), (ptr))


/* Memory Management Functions
 *
 *
 */

/** returns block memory to use at the current time
 *
 *  @return the block memory to use at the current time.
 */
SCIP_EXPORT
BMS_BLKMEM* SCIPblkmem(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns buffer memory for short living temporary objects
 *
 *  @return the buffer memory for short living temporary objects
 */
SCIP_EXPORT
BMS_BUFMEM* SCIPbuffer(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns clean buffer memory for short living temporary objects initialized to all zero
 *
 *  @return the buffer memory for short living temporary objects initialized to all zero
 */
SCIP_EXPORT
BMS_BUFMEM* SCIPcleanbuffer(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the total number of bytes used in block and buffer memory
 *
 *  @return the total number of bytes used in block and buffer memory.
 */
SCIP_EXPORT
SCIP_Longint SCIPgetMemUsed(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the total number of bytes in block and buffer memory
 *
 *  @return the total number of bytes in block and buffer memory.
 */
SCIP_EXPORT
SCIP_Longint SCIPgetMemTotal(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the estimated number of bytes used by external software, e.g., the LP solver
 *
 *  @return the estimated number of bytes used by external software, e.g., the LP solver.
 */
SCIP_EXPORT
SCIP_Longint SCIPgetMemExternEstim(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** calculate memory size for dynamically allocated arrays
 *
 *  @return the memory size for dynamically allocated arrays.
 */
SCIP_EXPORT
int SCIPcalcMemGrowSize(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   num                 /**< minimum number of entries to store */
   );

/** extends a dynamically allocated block memory array to be able to store at least the given number of elements;
 *  use SCIPensureBlockMemoryArray() define to call this method!
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPensureBlockMemoryArray_call(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                arrayptr,           /**< pointer to dynamically sized array */
   size_t                elemsize,           /**< size in bytes of each element in array */
   int*                  arraysize,          /**< pointer to current array size */
   int                   minsize             /**< required minimal array size */
   );

/** prints output about used memory */
SCIP_EXPORT
void SCIPprintMemoryDiagnostic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
