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
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPsortCompInt);

/* first all upwards-sorting methods */

/** sort an indexed element set in non-decreasing order, resulting in a permutation index array */
EXTERN
void SCIPsort(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-decreasing order */
EXTERN
void SCIPsortInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-decreasing order */
EXTERN
void SCIPsortPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );


/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrRealRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrRealReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/Bools, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
void SCIPsortReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
void SCIPsortInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntRealLong(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Longints in non-decreasing order */
EXTERN
void SCIPsortLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order */
EXTERN
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
EXTERN
void SCIPsortLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortDown(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-increasing order */
EXTERN
void SCIPsortDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-increasing order */
EXTERN
void SCIPsortDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
void SCIPsortDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< integer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );


/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
void SCIPsortDownInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Longints in non-increasing order */
EXTERN
void SCIPsortDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order */
EXTERN
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
EXTERN
void SCIPsortDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-decreasing order */
EXTERN
void SCIPsortedvecInsertPtr(
   void**                ptrarray,           /**< pointer to the pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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

/** insert a new element into three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an arrays of Reals, sorted in non-decreasing order */
EXTERN
void SCIPsortedvecInsertReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortedvecInsertLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
void SCIPsortedvecInsertDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/pointers, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
void SCIPsortedvecInsertDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertDownInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecInsertDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortedvecInsertDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecDelPosInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/RealsReals//ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrRealRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrRealBool(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecDelPosRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an arrays of Reals, sorted in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
void SCIPsortedvecDelPosRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array  where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecDelPosInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int*                  intarray3,          /**< third int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/SCIP_Real/SCIP_Longint, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntRealLong(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/pointers/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntPtrReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted by in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecDelPosLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecDelPosDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrRealBool(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecDelPosDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< integer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
void SCIPsortedvecDelPosDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecDelPosDownInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/ints, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int*                  intarray3,          /**< third int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosDownIntIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order */
EXTERN
void SCIPsortedvecDelPosDownIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three two arrays of Long/pointer, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order */
EXTERN
void SCIPsortedvecDelPosDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
void SCIPsortedvecDelPosDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
