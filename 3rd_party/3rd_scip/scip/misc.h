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
extern
SCIP_RETCODE SCIPrealarrayCreate(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the real array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of real values */
extern
SCIP_RETCODE SCIPrealarrayCopy(
   SCIP_REALARRAY**      realarray,          /**< pointer to store the copied real array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REALARRAY*       sourcerealarray     /**< dynamic real array to copy */
   );

/** frees a dynamic array of real values */
extern
SCIP_RETCODE SCIPrealarrayFree(
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPrealarrayExtend(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic real array */
extern
SCIP_RETCODE SCIPrealarrayClear(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** gets value of entry in dynamic array */
extern
SCIP_Real SCIPrealarrayGetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPrealarraySetVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
SCIP_RETCODE SCIPrealarrayIncVal(
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPrealarrayGetMinIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPrealarrayGetMaxIdx(
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   );

/** creates a dynamic array of int values */
extern
SCIP_RETCODE SCIPintarrayCreate(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the int array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of int values */
extern
SCIP_RETCODE SCIPintarrayCopy(
   SCIP_INTARRAY**       intarray,           /**< pointer to store the copied int array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_INTARRAY*        sourceintarray      /**< dynamic int array to copy */
   );

/** frees a dynamic array of int values */
extern
SCIP_RETCODE SCIPintarrayFree(
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPintarrayExtend(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic int array */
extern
SCIP_RETCODE SCIPintarrayClear(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** gets value of entry in dynamic array */
extern
int SCIPintarrayGetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPintarraySetVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
SCIP_RETCODE SCIPintarrayIncVal(
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPintarrayGetMinIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPintarrayGetMaxIdx(
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   );

/** creates a dynamic array of bool values */
extern
SCIP_RETCODE SCIPboolarrayCreate(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the bool array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of bool values */
extern
SCIP_RETCODE SCIPboolarrayCopy(
   SCIP_BOOLARRAY**      boolarray,          /**< pointer to store the copied bool array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_BOOLARRAY*       sourceboolarray     /**< dynamic bool array to copy */
   );

/** frees a dynamic array of bool values */
extern
SCIP_RETCODE SCIPboolarrayFree(
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPboolarrayExtend(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic bool array */
extern
SCIP_RETCODE SCIPboolarrayClear(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** gets value of entry in dynamic array */
extern
SCIP_Bool SCIPboolarrayGetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPboolarraySetVal(
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPboolarrayGetMinIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPboolarrayGetMaxIdx(
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** creates a dynamic array of pointer values */
extern
SCIP_RETCODE SCIPptrarrayCreate(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the ptr array */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a copy of a dynamic array of pointer values */
extern
SCIP_RETCODE SCIPptrarrayCopy(
   SCIP_PTRARRAY**       ptrarray,           /**< pointer to store the copied ptr array */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PTRARRAY*        sourceptrarray      /**< dynamic ptr array to copy */
   );

/** frees a dynamic array of pointer values */
extern
SCIP_RETCODE SCIPptrarrayFree(
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the ptr array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
SCIP_RETCODE SCIPptrarrayExtend(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic pointer array */
extern
SCIP_RETCODE SCIPptrarrayClear(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

/** gets value of entry in dynamic array */
extern
void* SCIPptrarrayGetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
SCIP_RETCODE SCIPptrarraySetVal(
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic ptr array */
   int                   arraygrowinit,      /**< initial size of array */
   SCIP_Real             arraygrowfac,       /**< growing factor of array */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   );

/** returns the minimal index of all stored non-zero elements */
extern
int SCIPptrarrayGetMinIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

/** returns the maximal index of all stored non-zero elements */
extern
int SCIPptrarrayGetMaxIdx(
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   );

/** creates a disjoint set (union find) structure \p uf for \p ncomponents many components (of size one) */
extern
SCIP_RETCODE SCIPdisjointsetCreate(
   SCIP_DISJOINTSET**    djset,              /**< disjoint set (union find) data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncomponents         /**< number of components */
   );

/** frees the disjoint set (union find) data structure */
extern
void SCIPdisjointsetFree(
   SCIP_DISJOINTSET**    djset,              /**< pointer to disjoint set (union find) data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** SCIP digraph functions
 *
 * internal digraph functions (see \ref DirectedGraph for public digraph methods)
 */

/** creates directed graph structure */
extern
SCIP_RETCODE SCIPdigraphCreate(
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   BMS_BLKMEM*           blkmem,             /**< block memory to store the data */
   int                   nnodes              /**< number of nodes */
   );

/** copies directed graph structure
 *
 *  @note The data in nodedata is copied verbatim. This possibly has to be adapted by the user.
 */
extern
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
extern
SCIP_Real SCIPnegateReal(
   SCIP_Real             x                   /**< value to negate */
   );

/** internal random number generator methods
 *
 * see \ref RandomNumbers for public random number generator methods
 */

/** creates and initializes a random number generator */
extern
SCIP_RETCODE SCIPrandomCreate(
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          initialseed         /**< initial random seed */
   );

/** frees a random number generator */
extern
void SCIPrandomFree(
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** initializes a random number generator with a given start seed */
extern
void SCIPrandomSetSeed(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   unsigned int          initseed            /**< initial random seed */
   );

#ifdef __cplusplus
}
#endif

#endif
