/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the library                         */
/*          BMS --- Block Memory Shell                                       */
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

/**@file   memory.h
 * @brief  memory allocation routines
 * @author Tobias Achterberg
 * @author Gerald Gamrath
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BMS_MEMORY_H__
#define __BMS_MEMORY_H__

#include <limits.h>
#include <stdlib.h>
#include <stddef.h>

/*
 * include build configuration flags
 */


#ifdef __cplusplus


/* special thanks to Daniel Junglas for following template and macros */

template<typename T> T* docast(T*, void *v);
template<typename T> T* docast(T*, void *v) { return reinterpret_cast<T*>(v); }

/* For C++11, we can easily check whether the types for memory functions like BMSduplicateXYZArray() are equal. */
#if __cplusplus > 199711L
#include <type_traits>

/* the following adds a type check for the parameters, used in ASSIGNCHECK below */
template<typename T1, typename T2> T1* docastcheck(T1* v1, void* v, T2* v2)
{
   typedef typename std::remove_const<T1>::type t1;
   typedef typename std::remove_const<T2>::type t2;
   static_assert(std::is_same<t1, t2>::value, "need equal types");
   return reinterpret_cast<T1*>(v);
}
#else
/* for older compilers do nothing */
template<typename T1, typename T2> T1* docastcheck(T1* v1, void* v, T2* v2) { return reinterpret_cast<T1*>(v); }
#endif


extern "C" {

#define ASSIGN(pointerstarstar, voidstarfunction) (*(pointerstarstar) = docast(*(pointerstarstar), (voidstarfunction)))
#define ASSIGNCHECK(pointerstarstar, voidstarfunction, origpointer) (*(pointerstarstar) = docastcheck(*(pointerstarstar), (voidstarfunction), (origpointer)))

#else

#define ASSIGN(pointerstarstar, voidstarfunction) (*(pointerstarstar) = (voidstarfunction))
#define ASSIGNCHECK(pointerstarstar, voidstarfunction, origpointer) (*(pointerstarstar) = (voidstarfunction))

#endif

/*
 * Define the macro SCIP_EXPORT depending if the OS is Windows or not
 */
#ifndef SCIP_EXPORT

#if defined(_WIN32) || defined(_WIN64)
#define SCIP_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
#define SCIP_EXPORT __attribute__((__visibility__("default")))
#else
#define SCIP_EXPORT
#endif

#endif

/* define if not already existing to make file independent from def.h */
#ifndef SCIP_UNUSED
#define SCIP_UNUSED(x) ((void) (x))
#endif


/*************************************************************************************
 * Standard Memory Management
 *
 * In debug mode, these methods extend malloc() and free() by logging all currently
 * allocated memory elements in an allocation list. This can be used as a simple leak
 * detection.
 *************************************************************************************/

/* Note: values that are passed as a size_t parameter are first converted to ptrdiff_t to be sure that negative numbers
 * are extended to the larger size. Then they are converted to size_t. Thus, negative numbers are converted to very
 * large size_t values. This is then checked within the functions. */

#define BMSallocMemory(ptr)                   ASSIGN((ptr), BMSallocMemory_call( sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSallocClearMemory(ptr)              ASSIGN((ptr), BMSallocClearMemory_call((size_t)(1), sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSallocMemorySize(ptr,size)          ASSIGN((ptr), BMSallocMemory_call( (size_t)(ptrdiff_t)(size), __FILE__, __LINE__ ))
#define BMSallocMemoryCPP(size)               BMSallocMemory_call( (size_t)(ptrdiff_t)(size), __FILE__, __LINE__ )
#define BMSallocClearMemorySize(ptr,size)     ASSIGN((ptr), BMSallocClearMemory_call((size_t)(1), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__ ))
#define BMSallocMemoryArray(ptr,num)          ASSIGN((ptr), BMSallocMemoryArray_call((size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSallocMemoryArrayCPP(num,size)      BMSallocMemoryArray_call( (size_t)(ptrdiff_t)(num), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__ )
#define BMSallocClearMemoryArray(ptr,num)     ASSIGN((ptr), BMSallocClearMemory_call((size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSreallocMemorySize(ptr,size)        ASSIGN((ptr), BMSreallocMemory_call((void*)(*(ptr)), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__ ))
#define BMSreallocMemoryArray(ptr,num)        ASSIGN((ptr), BMSreallocMemoryArray_call( *(ptr), (size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__ ))

#define BMSclearMemory(ptr)                   BMSclearMemory_call( (void*)(ptr), sizeof(*(ptr)) )
#define BMSclearMemoryArray(ptr, num)         BMSclearMemory_call( (void*)(ptr), (size_t)(ptrdiff_t)(num)*sizeof(*(ptr)) )
#define BMSclearMemorySize(ptr, size)         BMSclearMemory_call( (void*)(ptr), (size_t)(ptrdiff_t)(size) )

#define BMScopyMemory(ptr, source)            BMScopyMemory_call( (void*)(ptr), (const void*)(source), sizeof(*(ptr)) )
#define BMScopyMemoryArray(ptr, source, num)  BMScopyMemory_call( (void*)(ptr), (const void*)(source), (size_t)(ptrdiff_t)(num)*sizeof(*(ptr)) )
#define BMScopyMemorySize(ptr, source, size)  BMScopyMemory_call( (void*)(ptr), (const void*)(source), (size_t)(ptrdiff_t)(size) )

#define BMSmoveMemory(ptr, source)            BMSmoveMemory_call( (void*)(ptr), (const void*)(source), sizeof(*(ptr)) )
#define BMSmoveMemoryArray(ptr, source, num)  BMSmoveMemory_call( (void*)(ptr), (const void*)(source), (size_t)(ptrdiff_t)(num) * sizeof(*(ptr)) )
#define BMSmoveMemorySize(ptr, source, size)  BMSmoveMemory_call( (void*)(ptr), (const void*)(source), (size_t)(ptrdiff_t)(size) )

#define BMSduplicateMemory(ptr, source)       ASSIGN((ptr), BMSduplicateMemory_call( (const void*)(source), sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSduplicateMemorySize(ptr, source, size) ASSIGN((ptr), BMSduplicateMemory_call( (const void*)(source), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__ ))
#define BMSduplicateMemoryArray(ptr, source, num) ASSIGNCHECK((ptr), BMSduplicateMemoryArray_call( (const void*)(source), (size_t)(ptrdiff_t)(num), \
                                                  sizeof(**(ptr)), __FILE__, __LINE__ ), source)
#define BMSfreeMemory(ptr)                    BMSfreeMemory_call( (void**)(ptr), __FILE__, __LINE__ )
#define BMSfreeMemoryNull(ptr)                BMSfreeMemoryNull_call( (void**)(ptr), __FILE__, __LINE__ )
#define BMSfreeMemoryArray(ptr)               BMSfreeMemory_call( (void**)(ptr), __FILE__, __LINE__ )
#define BMSfreeMemoryArrayNull(ptr)           BMSfreeMemoryNull_call( (void**)(ptr), __FILE__, __LINE__ )
#define BMSfreeMemorySize(ptr)                BMSfreeMemory_call( (void**)(ptr), __FILE__, __LINE__ )
#define BMSfreeMemorySizeNull(ptr)            BMSfreeMemoryNull_call( (void**)(ptr), __FILE__, __LINE__ )

#ifndef NDEBUG
#define BMSgetPointerSize(ptr)                BMSgetPointerSize_call(ptr)
#define BMSdisplayMemory()                    BMSdisplayMemory_call()
#define BMScheckEmptyMemory()                 BMScheckEmptyMemory_call()
#define BMSgetMemoryUsed()                    BMSgetMemoryUsed_call()
#else
#define BMSgetPointerSize(ptr)                0
#define BMSdisplayMemory()                    /**/
#define BMScheckEmptyMemory()                 /**/
#define BMSgetMemoryUsed()                    0LL
#endif

/** allocates array and initializes it with 0; returns NULL if memory allocation failed */
SCIP_EXPORT
void* BMSallocClearMemory_call(
   size_t                num,                /**< number of memory element to allocate */
   size_t                typesize,           /**< size of memory element to allocate */
   const char*           filename,           /**< source file where the allocation is performed */
   int                   line                /**< line number in source file where the allocation is performed */
   );

/** allocates memory; returns NULL if memory allocation failed */
SCIP_EXPORT
void* BMSallocMemory_call(
   size_t                size,               /**< size of memory element to allocate */
   const char*           filename,           /**< source file where the allocation is performed */
   int                   line                /**< line number in source file where the allocation is performed */
   );

/** allocates array; returns NULL if memory allocation failed */
SCIP_EXPORT
void* BMSallocMemoryArray_call(
   size_t                num,                /**< number of components of array to allocate */
   size_t                typesize,           /**< size of each component */
   const char*           filename,           /**< source file where the allocation is performed */
   int                   line                /**< line number in source file where the allocation is performed */
   );

/** allocates memory; returns NULL if memory allocation failed */
SCIP_EXPORT
void* BMSreallocMemory_call(
   void*                 ptr,                /**< pointer to memory to reallocate */
   size_t                size,               /**< new size of memory element */
   const char*           filename,           /**< source file where the reallocation is performed */
   int                   line                /**< line number in source file where the reallocation is performed */
   );

/** reallocates array; returns NULL if memory allocation failed */
SCIP_EXPORT
void* BMSreallocMemoryArray_call(
   void*                 ptr,                /**< pointer to memory to reallocate */
   size_t                num,                /**< number of components of array to allocate */
   size_t                typesize,           /**< size of each component */
   const char*           filename,           /**< source file where the reallocation is performed */
   int                   line                /**< line number in source file where the reallocation is performed */
   );

/** clears a memory element (i.e. fills it with zeros) */
SCIP_EXPORT
void BMSclearMemory_call(
   void*                 ptr,                /**< pointer to memory element */
   size_t                size                /**< size of memory element */
   );

/** copies the contents of one memory element into another memory element */
SCIP_EXPORT
void BMScopyMemory_call(
   void*                 ptr,                /**< pointer to target memory element */
   const void*           source,             /**< pointer to source memory element */
   size_t                size                /**< size of memory element to copy */
   );

/** moves the contents of one memory element into another memory element, should be used if both elements overlap,
 *  otherwise BMScopyMemory is faster
 */
SCIP_EXPORT
void BMSmoveMemory_call(
   void*                 ptr,                /**< pointer to target memory element */
   const void*           source,             /**< pointer to source memory element */
   size_t                size                /**< size of memory element to copy */
   );

/** allocates memory and copies the contents of the given memory element into the new memory element */
SCIP_EXPORT
void* BMSduplicateMemory_call(
   const void*           source,             /**< pointer to source memory element */
   size_t                size,               /**< size of memory element to copy */
   const char*           filename,           /**< source file where the duplication is performed */
   int                   line                /**< line number in source file where the duplication is performed */
   );

/** allocates array and copies the contents of the given source array into the new array */
SCIP_EXPORT
void* BMSduplicateMemoryArray_call(
   const void*           source,             /**< pointer to source memory element */
   size_t                num,                /**< number of components of array to allocate */
   size_t                typesize,           /**< size of each component */
   const char*           filename,           /**< source file where the duplication is performed */
   int                   line                /**< line number in source file where the duplication is performed */
   );

/** frees an allocated memory element and sets pointer to NULL */
SCIP_EXPORT
void BMSfreeMemory_call(
   void**                ptr,                /**< pointer to pointer to memory element */
   const char*           filename,           /**< source file where the deallocation is performed */
   int                   line                /**< line number in source file where the deallocation is performed */
   );

/** frees an allocated memory element if pointer is not NULL and sets pointer to NULL */
SCIP_EXPORT
void BMSfreeMemoryNull_call(
   void**                ptr,                /**< pointer to pointer to memory element */
   const char*           filename,           /**< source file where the deallocation is performed */
   int                   line                /**< line number in source file where the deallocation is performed */
   );

/** returns the size of an allocated memory element */
SCIP_EXPORT
size_t BMSgetPointerSize_call(
   const void*           ptr                 /**< pointer to allocated memory */
   );

/** outputs information about currently allocated memory to the screen */
SCIP_EXPORT
void BMSdisplayMemory_call(
   void
   );

/** displays a warning message on the screen, if allocated memory exists */
SCIP_EXPORT
void BMScheckEmptyMemory_call(
   void
   );

/** returns total number of allocated bytes */
SCIP_EXPORT
long long BMSgetMemoryUsed_call(
   void
   );




/********************************************************************
 * Chunk Memory Management
 *
 * Efficient memory management for multiple objects of the same size
 ********************************************************************/

typedef struct BMS_ChkMem BMS_CHKMEM;           /**< collection of memory chunks of the same element size */


#ifndef BMS_NOBLOCKMEM

#define BMScreateChunkMemory(sz,isz,gbf)      BMScreateChunkMemory_call( (sz), (isz), (gbf), __FILE__, __LINE__ )
#define BMSclearChunkMemory(mem)              BMSclearChunkMemory_call( (mem), __FILE__, __LINE__ )
#define BMSdestroyChunkMemory(mem)            BMSdestroyChunkMemory_call( (mem), __FILE__, __LINE__ )

#define BMSallocChunkMemory(mem,ptr)          ASSIGN((ptr), BMSallocChunkMemory_call((mem), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSduplicateChunkMemory(mem, ptr, source) ASSIGN((ptr), BMSduplicateChunkMemory_call((mem), (const void*)(source), \
                                                sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSfreeChunkMemory(mem,ptr)           BMSfreeChunkMemory_call( (mem), (void**)(ptr), sizeof(**(ptr)), __FILE__, __LINE__ )
#define BMSfreeChunkMemoryNull(mem,ptr)       BMSfreeChunkMemoryNull_call( (mem), (void**)(ptr), sizeof(**(ptr)), __FILE__, __LINE__ )
#define BMSgarbagecollectChunkMemory(mem)     BMSgarbagecollectChunkMemory_call(mem)
#define BMSgetChunkMemoryUsed(mem)            BMSgetChunkMemoryUsed_call(mem)

#else

/* block memory management mapped to standard memory management */

#define BMScreateChunkMemory(sz,isz,gbf)           (void*)(0x01) /* dummy to not return a NULL pointer */
#define BMSclearChunkMemory(mem)                   /**/
#define BMSclearChunkMemoryNull(mem)               /**/
#define BMSdestroyChunkMemory(mem)                 /**/
#define BMSdestroyChunkMemoryNull(mem)             /**/
#define BMSallocChunkMemory(mem,ptr)               BMSallocMemory(ptr)
#define BMSduplicateChunkMemory(mem, ptr, source)  BMSduplicateMemory(ptr,source)
#define BMSfreeChunkMemory(mem,ptr)                BMSfreeMemory(ptr)
#define BMSfreeChunkMemoryNull(mem,ptr)            BMSfreeMemoryNull(ptr)
#define BMSgarbagecollectChunkMemory(mem)          /**/
#define BMSgetChunkMemoryUsed(mem)                 0LL

#endif


/** aligns the given byte size corresponding to the minimal alignment for chunk and block memory */
SCIP_EXPORT
void BMSalignMemsize(
   size_t*               size                /**< pointer to the size to align */
   );

/** checks whether the given size meets the alignment conditions for chunk and block memory  */
SCIP_EXPORT
int BMSisAligned(
   size_t                size                /**< size to check for alignment */
   );

/** creates a new chunk block data structure */
SCIP_EXPORT
BMS_CHKMEM* BMScreateChunkMemory_call(
   size_t                size,               /**< element size of the chunk block */
   int                   initchunksize,      /**< number of elements in the first chunk of the chunk block */
   int                   garbagefactor,      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** clears a chunk block data structure */
SCIP_EXPORT
void BMSclearChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** destroys and frees a chunk block data structure */
SCIP_EXPORT
void BMSdestroyChunkMemory_call(
   BMS_CHKMEM**          chkmem,             /**< pointer to chunk block */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates a memory element of the given chunk block */
SCIP_EXPORT
void* BMSallocChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** duplicates a given memory element by allocating a new element of the same chunk block and copying the data */
SCIP_EXPORT
void* BMSduplicateChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   const void*           source,             /**< source memory element */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees a memory element of the given chunk block and sets pointer to NULL */
SCIP_EXPORT
void BMSfreeChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   void**                ptr,                /**< pointer to pointer to memory element to free */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees a memory element of the given chunk block if pointer is not NULL and sets pointer to NULL */
SCIP_EXPORT
void BMSfreeChunkMemoryNull_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   void**                ptr,                /**< pointer to pointer to memory element to free */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** calls garbage collection of chunk block and frees chunks without allocated memory elements */
SCIP_EXPORT
void BMSgarbagecollectChunkMemory_call(
   BMS_CHKMEM*           chkmem              /**< chunk block */
   );

/** returns the number of allocated bytes in the chunk block */
SCIP_EXPORT
long long BMSgetChunkMemoryUsed_call(
   const BMS_CHKMEM*     chkmem              /**< chunk block */
   );




/***********************************************************
 * Block Memory Management
 *
 * Efficient memory management for objects of varying sizes
 ***********************************************************/

typedef struct BMS_BlkMem BMS_BLKMEM;           /**< block memory: collection of chunk blocks */

#ifndef BMS_NOBLOCKMEM

/* block memory methods for faster memory access */

/* Note: values that are passed as a size_t parameter are first converted to ptrdiff_t to be sure that negative numbers
 * are extended to the larger size. Then they are converted to size_t. Thus, negative numbers are converted to very
 * large size_t values. This is then checked within the functions. */

#define BMScreateBlockMemory(csz,gbf)         BMScreateBlockMemory_call( (csz), (gbf), __FILE__, __LINE__ )
#define BMSclearBlockMemory(mem)              BMSclearBlockMemory_call( (mem), __FILE__, __LINE__ )
#define BMSdestroyBlockMemory(mem)            BMSdestroyBlockMemory_call( (mem), __FILE__, __LINE__ )

#define BMSallocBlockMemory(mem,ptr)          ASSIGN((ptr), BMSallocBlockMemory_call((mem), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSallocClearBlockMemory(mem,ptr)     ASSIGN((ptr), BMSallocClearBlockMemory_call((mem), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSallocBlockMemorySize(mem,ptr,size) ASSIGN((ptr), BMSallocBlockMemory_call((mem), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__))
#define BMSallocBlockMemoryArray(mem,ptr,num) ASSIGN((ptr), BMSallocBlockMemoryArray_call((mem), (size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSallocClearBlockMemoryArray(mem,ptr,num) ASSIGN((ptr), BMSallocClearBlockMemoryArray_call((mem), (size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSreallocBlockMemorySize(mem,ptr,oldsize,newsize) ASSIGN((ptr), BMSreallocBlockMemory_call((mem), (void*)(*(ptr)), \
                                                (size_t)(ptrdiff_t)(oldsize), (size_t)(ptrdiff_t)(newsize), __FILE__, __LINE__))
#define BMSreallocBlockMemoryArray(mem,ptr,oldnum,newnum) ASSIGN((ptr), BMSreallocBlockMemoryArray_call((mem), (void*)(*(ptr)), \
                                                (size_t)(ptrdiff_t)(oldnum), (size_t)(ptrdiff_t)(newnum), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSduplicateBlockMemory(mem, ptr, source) ASSIGN((ptr), BMSduplicateBlockMemory_call((mem), (const void*)(source), \
                                                sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSduplicateBlockMemoryArray(mem, ptr, source, num) ASSIGNCHECK((ptr), BMSduplicateBlockMemoryArray_call( (mem), (const void*)(source), \
                                                (size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__ ), source)

#define BMSfreeBlockMemory(mem,ptr)           BMSfreeBlockMemory_call( (mem), (void**)(ptr), sizeof(**(ptr)), __FILE__, __LINE__ )
#define BMSfreeBlockMemoryNull(mem,ptr)       BMSfreeBlockMemoryNull_call( (mem), (void**)(ptr), sizeof(**(ptr)), __FILE__, __LINE__ )
#define BMSfreeBlockMemoryArray(mem,ptr,num)  BMSfreeBlockMemory_call( (mem), (void**)(ptr), (num)*sizeof(**(ptr)), __FILE__, __LINE__ )
#define BMSfreeBlockMemoryArrayNull(mem,ptr,num) BMSfreeBlockMemoryNull_call( (mem), (void**)(ptr), (num)*sizeof(**(ptr)), __FILE__, __LINE__ )
#define BMSfreeBlockMemorySize(mem,ptr,size)  BMSfreeBlockMemory_call( (mem), (void**)(ptr), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__ )
#define BMSfreeBlockMemorySizeNull(mem,ptr,size) BMSfreeBlockMemory_call( (mem), (void**)(ptr), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__ )

#define BMSgarbagecollectBlockMemory(mem)     BMSgarbagecollectBlockMemory_call(mem)
#define BMSgetBlockMemoryAllocated(mem)       BMSgetBlockMemoryAllocated_call(mem)
#define BMSgetBlockMemoryUsed(mem)            BMSgetBlockMemoryUsed_call(mem)
#define BMSgetBlockMemoryUnused(mem)          BMSgetBlockMemoryUnused_call(mem)
#define BMSgetBlockMemoryUsedMax(mem)         BMSgetBlockMemoryUsedMax_call(mem)
#define BMSgetBlockMemoryUnusedMax(mem)       BMSgetBlockMemoryUnusedMax_call(mem)
#define BMSgetBlockMemoryAllocatedMax(mem)    BMSgetBlockMemoryAllocatedMax_call(mem)
#define BMSgetBlockPointerSize(mem,ptr)       BMSgetBlockPointerSize_call((mem), (ptr))
#define BMSdisplayBlockMemory(mem)            BMSdisplayBlockMemory_call(mem)
#define BMSblockMemoryCheckEmpty(mem)         BMScheckEmptyBlockMemory_call(mem)

#else

/* block memory management mapped to standard memory management */

#define BMScreateBlockMemory(csz,gbf)                        (SCIP_UNUSED(csz), SCIP_UNUSED(gbf), (void*)(0x01)) /* dummy to not return a NULL pointer */
#define BMSclearBlockMemory(mem)                             SCIP_UNUSED(mem)
#define BMSclearBlockMemoryNull(mem)                         SCIP_UNUSED(mem)
#define BMSdestroyBlockMemory(mem)                           SCIP_UNUSED(mem)
#define BMSdestroyBlockMemoryNull(mem)                       SCIP_UNUSED(mem)
#define BMSallocBlockMemory(mem,ptr)                         (SCIP_UNUSED(mem), BMSallocMemory(ptr))
#define BMSallocClearBlockMemory(mem,ptr)                    (SCIP_UNUSED(mem), BMSallocClearMemory(ptr))
#define BMSallocBlockMemoryArray(mem,ptr,num)                (SCIP_UNUSED(mem), BMSallocMemoryArray(ptr,num))
#define BMSallocClearBlockMemoryArray(mem,ptr,num)           (SCIP_UNUSED(mem), BMSallocClearMemoryArray(ptr,num))
#define BMSallocBlockMemorySize(mem,ptr,size)                (SCIP_UNUSED(mem), BMSallocMemorySize(ptr,size))
#define BMSreallocBlockMemoryArray(mem,ptr,oldnum,newnum)    (SCIP_UNUSED(mem), SCIP_UNUSED(oldnum), BMSreallocMemoryArray(ptr,newnum))
#define BMSreallocBlockMemorySize(mem,ptr,oldsize,newsize)   (SCIP_UNUSED(mem), SCIP_UNUSED(oldsize), BMSreallocMemorySize(ptr,newsize))
#define BMSduplicateBlockMemory(mem, ptr, source)            (SCIP_UNUSED(mem), BMSduplicateMemory(ptr,source))
#define BMSduplicateBlockMemoryArray(mem, ptr, source, num)  (SCIP_UNUSED(mem), BMSduplicateMemoryArray(ptr,source,num))
#define BMSfreeBlockMemory(mem,ptr)                          (SCIP_UNUSED(mem), BMSfreeMemory(ptr))
#define BMSfreeBlockMemoryNull(mem,ptr)                      (SCIP_UNUSED(mem), BMSfreeMemoryNull(ptr))
#define BMSfreeBlockMemoryArray(mem,ptr,num)                 (SCIP_UNUSED(mem), SCIP_UNUSED(num), BMSfreeMemoryArray(ptr))
#define BMSfreeBlockMemoryArrayNull(mem,ptr,num)             (SCIP_UNUSED(mem), SCIP_UNUSED(num), BMSfreeMemoryArrayNull(ptr))
#define BMSfreeBlockMemorySize(mem,ptr,size)                 (SCIP_UNUSED(mem), SCIP_UNUSED(size), BMSfreeMemory(ptr))
#define BMSfreeBlockMemorySizeNull(mem,ptr,size)             (SCIP_UNUSED(mem), SCIP_UNUSED(size), BMSfreeMemoryNull(ptr))
#define BMSgarbagecollectBlockMemory(mem)                    SCIP_UNUSED(mem)
#define BMSgetBlockMemoryAllocated(mem)                      (SCIP_UNUSED(mem), 0LL)
#define BMSgetBlockMemoryUsed(mem)                           (SCIP_UNUSED(mem), 0LL)
#define BMSgetBlockMemoryUnused(mem)                         (SCIP_UNUSED(mem), 0LL)
#define BMSgetBlockMemoryUsedMax(mem)                        (SCIP_UNUSED(mem), 0LL)
#define BMSgetBlockMemoryUnusedMax(mem)                      (SCIP_UNUSED(mem), 0LL)
#define BMSgetBlockMemoryAllocatedMax(mem)                   (SCIP_UNUSED(mem), 0LL)
#define BMSgetBlockPointerSize(mem,ptr)                      (SCIP_UNUSED(mem), SCIP_UNUSED(ptr), 0)
#define BMSdisplayBlockMemory(mem)                           SCIP_UNUSED(mem)
#define BMSblockMemoryCheckEmpty(mem)                        (SCIP_UNUSED(mem), 0LL)

#endif


/** creates a block memory allocation data structure */
SCIP_EXPORT
BMS_BLKMEM* BMScreateBlockMemory_call(
   int                   initchunksize,      /**< number of elements in the first chunk of each chunk block */
   int                   garbagefactor,      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees all chunk blocks in the block memory */
SCIP_EXPORT
void BMSclearBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** clears and deletes block memory */
SCIP_EXPORT
void BMSdestroyBlockMemory_call(
   BMS_BLKMEM**          blkmem,             /**< pointer to block memory */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates memory in the block memory pool */
SCIP_EXPORT
void* BMSallocBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   size_t                size,               /**< size of memory element to allocate */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates memory in the block memory pool and clears it */
SCIP_EXPORT
void* BMSallocClearBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   size_t                size,               /**< size of memory element to allocate */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates array in the block memory pool */
SCIP_EXPORT
void* BMSallocBlockMemoryArray_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   size_t                num,                /**< size of array to be allocated */
   size_t                typesize,           /**< size of each component */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates array in the block memory pool and clears it */
SCIP_EXPORT
void* BMSallocClearBlockMemoryArray_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   size_t                num,                /**< size of array to be allocated */
   size_t                typesize,           /**< size of each component */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** resizes memory element in the block memory pool and copies the data */
SCIP_EXPORT
void* BMSreallocBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 ptr,                /**< memory element to reallocated */
   size_t                oldsize,            /**< old size of memory element */
   size_t                newsize,            /**< new size of memory element */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** resizes array in the block memory pool and copies the data */
SCIP_EXPORT
void* BMSreallocBlockMemoryArray_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 ptr,                /**< memory element to reallocated */
   size_t                oldnum,             /**< old size of array */
   size_t                newnum,             /**< new size of array */
   size_t                typesize,           /**< size of each component */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** duplicates memory element in the block memory pool and copies the data */
SCIP_EXPORT
void* BMSduplicateBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const void*           source,             /**< memory element to duplicate */
   size_t                size,               /**< size of memory elements */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** duplicates array in the block memory pool and copies the data */
SCIP_EXPORT
void* BMSduplicateBlockMemoryArray_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const void*           source,             /**< memory element to duplicate */
   size_t                num,                /**< size of array to be duplicated */
   size_t                typesize,           /**< size of each component */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees memory element in the block memory pool and sets pointer to NULL */
SCIP_EXPORT
void BMSfreeBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void**                ptr,                /**< pointer to pointer to memory element to free */
   size_t                size,               /**< size of memory element */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees memory element in the block memory pool if pointer is not NULL and sets pointer to NULL */
SCIP_EXPORT
void BMSfreeBlockMemoryNull_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void**                ptr,                /**< pointer to pointer to memory element to free */
   size_t                size,               /**< size of memory element */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** calls garbage collection of block memory, frees chunks without allocated memory elements, and frees
 *  chunk blocks without any chunks
 */
SCIP_EXPORT
void BMSgarbagecollectBlockMemory_call(
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns the number of allocated bytes in the block memory */
SCIP_EXPORT
long long BMSgetBlockMemoryAllocated_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** returns the number of used bytes in the block memory */
SCIP_EXPORT
long long BMSgetBlockMemoryUsed_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** returns the number of allocated but not used bytes in the block memory */
SCIP_EXPORT
long long BMSgetBlockMemoryUnused_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** returns the maximal number of used bytes in the block memory */
SCIP_EXPORT
long long BMSgetBlockMemoryUsedMax_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** returns the maximal number of allocated but not used bytes in the block memory */
SCIP_EXPORT
long long BMSgetBlockMemoryUnusedMax_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** returns the maximal number of allocated bytes in the block memory */
long long BMSgetBlockMemoryAllocatedMax_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** returns the size of the given memory element; returns 0, if the element is not member of the block memory */
SCIP_EXPORT
size_t BMSgetBlockPointerSize_call(
   const BMS_BLKMEM*     blkmem,             /**< block memory */
   const void*           ptr                 /**< memory element */
   );

/** outputs allocation diagnostics of block memory */
SCIP_EXPORT
void BMSdisplayBlockMemory_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** outputs error messages, if there are allocated elements in the block memory and returns number of unfreed bytes */
SCIP_EXPORT
long long BMScheckEmptyBlockMemory_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );





/***********************************************************
 * Buffer Memory Management
 *
 * Efficient memory management for temporary objects
 ***********************************************************/

typedef struct BMS_BufMem BMS_BUFMEM;        /**< buffer memory for temporary objects */

/* Note: values that are passed as a size_t parameter are first converted to ptrdiff_t to be sure that negative numbers
 * are extended to the larger size. Then they are converted to size_t. Thus, negative numbers are converted to very
 * large size_t values. This is then checked within the functions. */

#define BMSallocBufferMemory(mem,ptr)        ASSIGN((ptr), BMSallocBufferMemory_call((mem), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSallocBufferMemorySize(mem,ptr,size) ASSIGN((ptr), BMSallocBufferMemory_call((mem), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__))
#define BMSreallocBufferMemorySize(mem,ptr,size) \
                                             ASSIGN((ptr), BMSreallocBufferMemory_call((mem), (void*)(*(ptr)), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__))
#define BMSallocBufferMemoryArray(mem,ptr,num) ASSIGN((ptr), BMSallocBufferMemoryArray_call((mem), (size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSallocClearBufferMemoryArray(mem,ptr,num) ASSIGN((ptr), BMSallocClearBufferMemoryArray_call((mem), (size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSreallocBufferMemoryArray(mem,ptr,num) ASSIGN((ptr), BMSreallocBufferMemoryArray_call((mem), (void*)(*(ptr)), (size_t)(ptrdiff_t)(num), \
                                                 sizeof(**(ptr)), __FILE__, __LINE__))
#define BMSduplicateBufferMemory(mem,ptr,source,size) \
                                             ASSIGN((ptr), BMSduplicateBufferMemory_call((mem), (const void*)(source), (size_t)(ptrdiff_t)(size), __FILE__, __LINE__))
#define BMSduplicateBufferMemoryArray(mem,ptr,source,num) ASSIGNCHECK((ptr), BMSduplicateBufferMemoryArray_call((mem), \
                                                 (const void*)(source), (size_t)(ptrdiff_t)(num), sizeof(**(ptr)), __FILE__, __LINE__), source)

#define BMSfreeBufferMemory(mem,ptr)         BMSfreeBufferMemory_call((mem), (void**)(ptr), __FILE__, __LINE__)
#define BMSfreeBufferMemoryNull(mem,ptr)     BMSfreeBufferMemoryNull_call((mem), (void**)(ptr), __FILE__, __LINE__)
#define BMSfreeBufferMemoryArray(mem,ptr)    BMSfreeBufferMemory_call((mem), (void**)(ptr), __FILE__, __LINE__)
#define BMSfreeBufferMemoryArrayNull(mem,ptr) BMSfreeBufferMemoryNull_call((mem), (void**)(ptr), __FILE__, __LINE__)
#define BMSfreeBufferMemorySize(mem,ptr)     BMSfreeBufferMemory_call((mem), (void**)(ptr), __FILE__, __LINE__);
#define BMSfreeBufferMemorySizeNull(mem,ptr) BMSfreeBufferMemoryNull_call((mem), (void**)(ptr), __FILE__, __LINE__)

#define BMScreateBufferMemory(fac,init,clean) BMScreateBufferMemory_call((fac), (init), (clean), __FILE__, __LINE__)
#define BMSdestroyBufferMemory(mem)          BMSdestroyBufferMemory_call((mem), __FILE__, __LINE__)


/** creates memory buffer storage */
SCIP_EXPORT
BMS_BUFMEM* BMScreateBufferMemory_call(
   double                arraygrowfac,       /**< memory growing factor for dynamically allocated arrays */
   int                   arraygrowinit,      /**< initial size of dynamically allocated arrays */
   unsigned int          clean,              /**< should the memory blocks in the buffer be initialized to zero? */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** destroys buffer memory */
SCIP_EXPORT
void BMSdestroyBufferMemory_call(
   BMS_BUFMEM**          buffer,             /**< pointer to memory buffer storage */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** set arraygrowfac */
SCIP_EXPORT
void BMSsetBufferMemoryArraygrowfac(
   BMS_BUFMEM*           buffer,             /**< pointer to memory buffer storage */
   double                arraygrowfac        /**< memory growing factor for dynamically allocated arrays */
   );

/** set arraygrowinit */
SCIP_EXPORT
void BMSsetBufferMemoryArraygrowinit(
   BMS_BUFMEM*           buffer,             /**< pointer to memory buffer storage */
   int                   arraygrowinit       /**< initial size of dynamically allocated arrays */
   );

/** allocates the next unused buffer */
SCIP_EXPORT
void* BMSallocBufferMemory_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   size_t                size,               /**< minimal required size of the buffer */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates the next unused buffer array */
SCIP_EXPORT
void* BMSallocBufferMemoryArray_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   size_t                num,                /**< size of array to be allocated */
   size_t                typesize,           /**< size of components */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates the next unused buffer and clears it */
SCIP_EXPORT
void* BMSallocClearBufferMemoryArray_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   size_t                num,                /**< size of array to be allocated */
   size_t                typesize,           /**< size of components */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** reallocates the buffer to at least the given size */
SCIP_EXPORT
void* BMSreallocBufferMemory_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   void*                 ptr,                /**< pointer to the allocated memory buffer */
   size_t                size,               /**< minimal required size of the buffer */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** reallocates an array in the buffer to at least the given size */
SCIP_EXPORT
void* BMSreallocBufferMemoryArray_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   void*                 ptr,                /**< pointer to the allocated memory buffer */
   size_t                num,                /**< size of array to be allocated */
   size_t                typesize,           /**< size of components */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates the next unused buffer and copies the given memory into the buffer */
SCIP_EXPORT
void* BMSduplicateBufferMemory_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   const void*           source,             /**< memory block to copy into the buffer */
   size_t                size,               /**< minimal required size of the buffer */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates an array in the next unused buffer and copies the given memory into the buffer */
SCIP_EXPORT
void* BMSduplicateBufferMemoryArray_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   const void*           source,             /**< memory block to copy into the buffer */
   size_t                num,                /**< size of array to be allocated */
   size_t                typesize,           /**< size of components */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees a buffer and sets pointer to NULL */
SCIP_EXPORT
void BMSfreeBufferMemory_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   void**                ptr,                /**< pointer to pointer to the allocated memory buffer */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees a buffer if pointer is not NULL and sets pointer to NULL */
SCIP_EXPORT
void BMSfreeBufferMemoryNull_call(
   BMS_BUFMEM*           buffer,             /**< memory buffer storage */
   void**                ptr,                /**< pointer to pointer to the allocated memory buffer */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** gets number of used buffers */
SCIP_EXPORT
size_t BMSgetNUsedBufferMemory(
   BMS_BUFMEM*           buffer              /**< memory buffer storage */
   );

/** returns the number of allocated bytes in the buffer memory */
SCIP_EXPORT
long long BMSgetBufferMemoryUsed(
   const BMS_BUFMEM*     bufmem              /**< buffer memory */
   );

/** outputs statistics about currently allocated buffers to the screen */
SCIP_EXPORT
void BMSprintBufferMemory(
   BMS_BUFMEM*           buffer              /**< memory buffer storage */
   );


#ifdef __cplusplus
}
#endif

#endif
