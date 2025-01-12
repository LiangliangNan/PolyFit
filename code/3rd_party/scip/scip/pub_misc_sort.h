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

/**@file   pub_misc_sort.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods for sorting joint arrays of various types
 * @author Gregor Hendel
 *
 * This file contains methods for sorting joint arrays of various types.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MISC_SORT_H__
#define __SCIP_PUB_MISC_SORT_H__

#include "scip/def.h"
#include "type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Sorting algorithms
 */

/**@defgroup SortingAlgorithms Sorting Algorithms
 * @ingroup MiscellaneousMethods
 * @brief public methods for in place sorting of arrays
 *
 * Below are the public methods for in place sorting of up to six arrays of joint data.
 *
 * @{
 */

/** default comparer for integers */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsortCompInt);

/* first all upwards-sorting methods */

/** sort an indexed element set in non-decreasing order, resulting in a permutation index array */
SCIP_EXPORT
void SCIPsort(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-decreasing order */
SCIP_EXPORT
void SCIPsortInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );


/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrRealRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/Reals/SCIP_Bools/SCIP_Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrRealRealBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/Reals/ints/SCIP_Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrRealRealIntBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrRealReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Reals in non-decreasing order */
SCIP_EXPORT
void SCIPsortReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort array of ints in non-decreasing order */
SCIP_EXPORT
void SCIPsortInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntRealLong(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Longints in non-decreasing order */
SCIP_EXPORT
void SCIPsortLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/* now all downwards-sorting methods */

/** sort an indexed element set in non-increasing order, resulting in a permutation index array */
SCIP_EXPORT
void SCIPsortDown(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Reals in non-increasing order */
SCIP_EXPORT
void SCIPsortDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< integer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );


/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort array of ints in non-increasing order */
SCIP_EXPORT
void SCIPsortDownInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Longints in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/*
 * Sorted vectors
 */

/* upwards insertion */

/** insert a new element into an index array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtr(
   void**                ptrarray,           /**< pointer to the pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrRealRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/Reals/SCIP_Bools/SCIP_Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrRealRealBoolBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/Reals/ints/SCIP_Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrRealRealIntBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrRealBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Bool             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Bools/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an arrays of Reals, sorted in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted  */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   intval,             /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< pointer array where an element is to be inserted */
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   intval1,            /**< additional value of new element */
   int                   intval2,            /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Longint          field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealRealIntInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   void*                 field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   void*                 field5val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of ints in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int*                  intarray3,          /**< third int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/SCIP_Real/SCIP_Longint, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntRealLong(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntIntPtr(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntIntReal(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/pointers/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntPtrReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntIntIntPtr(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int*                  intarray3,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntIntIntReal(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int*                  intarray3,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntPtrIntReal(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Longints, sorted in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be inserted  */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into six joint arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be inserted  */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   SCIP_Bool             field5val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecInsertIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   SCIP_Bool             field5val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );


/* downwards insertion */

/** insert a new element into an index array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrRealBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Bool             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );


/** insert a new element into four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Reals, sorted in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );


/** insert a new element into three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Reals/Pointer/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   intval,             /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< pointer array where an element is to be inserted */
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   intval1,            /**< additional value of new element */
   int                   intval2,            /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array  where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Longint          field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   void*                 field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   void*                 field5val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of ints in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int*                  intarray3,          /**< third int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   int*                  intarray3,          /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/int/ints/reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   int*                  intarray3,          /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Longints, sorted in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be inserted  */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into six joint arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be inserted  */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   SCIP_Bool             field5val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecInsertDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increased order */
SCIP_EXPORT
void SCIPsortedvecInsertDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   SCIP_Bool             field4val,          /**< additional value of new element */
   SCIP_Bool             field5val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/* upwards position deletion */

/** delete the element at the given position from an index array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/RealsReals//ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrRealRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/RealsReals/SCIP_Bools/SCIP_Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrRealRealBoolBool(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array where an element is to be deleted */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/RealsReals/ints/SCIP_Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrRealRealIntBool(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrRealBool(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Bools/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an arrays of Reals, sorted in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array  where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array where an element is to be deleted */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of ints in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int*                  intarray3,          /**< third int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/SCIP_Real/SCIP_Longint, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntRealLong(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/pointers/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntPtrReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted by in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from six joint arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/* downwards position deletion */

/** delete the element at the given position from an index array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrRealBool(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Reals, sorted in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< integer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array where an element is to be deleted */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of ints in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int*                  intarray3,          /**< third int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two arrays of Long/pointer, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from six joint arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
SCIP_EXPORT
void SCIPsortedvecDelPosDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/* upwards binary search */

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindInd(
   int*                  indarray,           /**< index array to be searched */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindPtr(
   void**                ptrarray,           /**< pointer array to be searched */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be searched */
   SCIP_Real             val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindInt(
   int*                  intarray,           /**< int array to be searched */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be searched */
   SCIP_Longint          val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );


/* downwards binary search */

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindDownInd(
   int*                  indarray,           /**< index array to be searched */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindDownPtr(
   void**                ptrarray,           /**< pointer array to be searched */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be searched */
   SCIP_Real             val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindDownInt(
   int*                  intarray,           /**< int array to be searched */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_EXPORT
SCIP_Bool SCIPsortedvecFindDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be searched */
   SCIP_Longint          val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
