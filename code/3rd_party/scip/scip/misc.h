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

/**@file   misc.h
 * @ingroup INTERNALAPI
 * @brief  internal miscellaneous methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MISC_H__
#define __SCIP_MISC_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_misc.h"
#include "scip/pub_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Dynamic Arrays
 */

/** creates a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayCreate(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayCopy(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REALARRAY*       sourcerealarray     /**< dynamic real array to copy */
   );

/** frees a dynamic array of real values */
SCIP_RETCODE SCIPrealarrayFree(
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPrealarrayExtend(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic real array */
SCIP_RETCODE SCIPrealarrayClear(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** gets value of entry in dynamic array */
SCIP_Real SCIPrealarrayGetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPrealarraySetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPrealarrayIncVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
int SCIPrealarrayGetMinIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** returns the maximal index of all stored non-zero elements */
int SCIPrealarrayGetMaxIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** creates a dynamic array of int values */
SCIP_RETCODE SCIPintarrayCreate(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the int array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of int values */
SCIP_RETCODE SCIPintarrayCopy(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the copied int array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_INTARRAY*        sourceintarray      /**< dynamic int array to copy */
   );

/** frees a dynamic array of int values */
SCIP_RETCODE SCIPintarrayFree(
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPintarrayExtend(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic int array */
SCIP_RETCODE SCIPintarrayClear(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** gets value of entry in dynamic array */
int SCIPintarrayGetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPintarraySetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPintarrayIncVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
int SCIPintarrayGetMinIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** returns the maximal index of all stored non-zero elements */
int SCIPintarrayGetMaxIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** creates a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayCreate(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the bool array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayCopy(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the copied bool array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_BOOLARRAY*       sourceboolarray     /**< dynamic bool array to copy */
   );

/** frees a dynamic array of bool values */
SCIP_RETCODE SCIPboolarrayFree(
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPboolarrayExtend(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic bool array */
SCIP_RETCODE SCIPboolarrayClear(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** gets value of entry in dynamic array */
SCIP_Bool SCIPboolarrayGetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPboolarraySetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   );

/** returns the minimal index of all stored non-zero elements */
int SCIPboolarrayGetMinIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** returns the maximal index of all stored non-zero elements */
int SCIPboolarrayGetMaxIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** creates a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayCreate(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the ptr array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayCopy(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the copied ptr array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PTRARRAY*        sourceptrarray      /**< dynamic ptr array to copy */
   );

/** frees a dynamic array of pointer values */
SCIP_RETCODE SCIPptrarrayFree(
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the ptr array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPptrarrayExtend(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic pointer array */
SCIP_RETCODE SCIPptrarrayClear(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

/** gets value of entry in dynamic array */
void* SCIPptrarrayGetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPptrarraySetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   );

/** returns the minimal index of all stored non-zero elements */
int SCIPptrarrayGetMinIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

/** returns the maximal index of all stored non-zero elements */
int SCIPptrarrayGetMaxIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );


/* SCIP disjoint set data structure
 *
 * internal disjoint set functions (see \ref DisjointSet for public methods)
 */

/** creates a disjoint set (union find) structure \p djset for \p ncomponents many components (of size one) */
SCIP_RETCODE SCIPdisjointsetCreate(
   SCIP_DISJOINTSET**    djset,              /**< disjoint set (union find) data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncomponents         /**< number of components */
   );

/** frees the disjoint set (union find) data structure */
void SCIPdisjointsetFree(
   SCIP_DISJOINTSET**    djset,              /**< pointer to disjoint set (union find) data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );


/** SCIP digraph functions
 *
 * internal digraph functions (see \ref DirectedGraph for public digraph methods)
 */

/** creates directed graph structure */
SCIP_RETCODE SCIPdigraphCreate(
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   BMS_BLKMEM*           blkmem,             /**< block memory to store the data */
   int                   nnodes              /**< number of nodes */
   );

/** copies directed graph structure
 *
 *  @note The data in nodedata is copied verbatim. This possibly has to be adapted by the user.
 */
SCIP_RETCODE SCIPdigraphCopy(
   SCIP_DIGRAPH**        targetdigraph,      /**< pointer to store the copied directed graph */
   SCIP_DIGRAPH*         sourcedigraph,      /**< source directed graph */
   BMS_BLKMEM*           targetblkmem        /**< block memory to store the target block memory, or NULL to use the same
                                              *  the same block memory as used for the \p sourcedigraph */
   );

/*
 * Additional math functions
 */

/** negates a number
 *
 * negation of a number that can be used to avoid that a negation is optimized away by a compiler
 */
SCIP_Real SCIPnegateReal(
   SCIP_Real             x                   /**< value to negate */
   );

/** internal random number generator methods
 *
 * see \ref RandomNumbers for public random number generator methods
 */

/** creates and initializes a random number generator */
SCIP_EXPORT
SCIP_RETCODE SCIPrandomCreate(
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          initialseed         /**< initial random seed */
   );

/** frees a random number generator */
SCIP_EXPORT
void SCIPrandomFree(
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** initializes a random number generator with a given start seed */
SCIP_EXPORT
void SCIPrandomSetSeed(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   unsigned int          initseed            /**< initial random seed */
   );

#ifdef __cplusplus
}
#endif

#endif
