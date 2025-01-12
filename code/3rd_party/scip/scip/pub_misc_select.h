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

/**@file   pub_misc_select.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods for selecting (weighted) k-medians
 * @author Gregor Hendel
 *
 * This file contains headers for selecting (weighted) k-medians
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MISC_SELECT_H__
#define __SCIP_PUB_MISC_SELECT_H__

#include "scip/def.h"
#include "type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Selection and weighted selection algorithms
 */

/**@defgroup SelectionAlgorithms Algorithms for (Weighted) Median Selection
 * @ingroup MiscellaneousMethods
 * @brief public methods for the selection of (weighted) k-median.
 *
 * The methods in this group perform a selection of the (weighted)  \f$ k \f$-median from an unsorted array of elements.
 * The necessary element swaps are performed in-place on the array of keys.
 * The necessary permutations are also performed on up to six associated arrays.
 *
 * For methods that perform complete in place sorting, see \ref SortingAlgorithms.
 *
 * For an array a containing n elements \f$ a[0], ..., a[n-1] \f$ and an integer \f$ 0 \leq k \leq n - 1 \f$ , we call an element
 *  \f$ a[i] \f$   \f$ k \f$-median if
 * there exists a permutation  \f$ \pi \f$  of the array indices such that  \f$ \pi(i) = k  \f$
 * and  \f$ a[\pi^{-1}(j)] \leq a[i]  \f$
 * for  \f$ j = 0, \dots, k-1  \f$ and  \f$ a[\pi^{-1}(j)] > a[i]  \f$ for  \f$ j = k + 1,\dots,n - 1 \f$.
 * The  \f$ k \f$-median is hence an element that would appear at position  \f$ k  \f$ after sorting the input array.
 * Note that there may exist several  \f$ k \f$-medians if the array elements are not unique, only its key value  \f$ a[i] \f$.
 *
 * In order to determine the  \f$ k \f$-median, the algorithm selects a pivot element and determines the array position for
 * this pivot like quicksort. In contrast to quicksort, however, one recursion can be saved during the selection process.
 * After a single iteration that placed the pivot at position  \f$ p \f$ , the algorithm either terminates if  \f$ p = k \f$,
 * or it continues in the left half of the array if  \f$ p > k \f$, or in the right half of the array if  \f$ p < k \f$.
 *
 * After the algorithm terminates, the  \f$ k \f$-median can be accessed by accessing the array element at position  \f$ k \f$.
 *
 * A weighted median denotes the generalization of the  \f$ k \f$-median to arbitrary, nonnegative associated
 * weights \f$ w[0], \dots, w[n-1] \in \mathbb{R}\f$ and a capacity  \f$ 0 \leq C \in \mathbb{R} \f$. An element  \f$ a[i] \f$
 * is called weighted median if there exists a permutation that satisfies the same weak sorting as above and in addition
 * \f$ W:= \sum_{j = 0}^{k - 1}w[\pi^{-1}(j)] < C\f$, but  \f$ W + w[i] \geq C\f$. In other words, the weighted median
 * is the first element in the weak sorting such that its weight together with the sum of all preceding item weights
 * reach or exceed the given capacity \f$ C \f$. If all weights are equal to \f$ 1 \f$ and the capacity is  \f$ C = k + 0.5\f$,
 * the weighted median becomes the  \f$ k \f$-median.
 *
 * @{
 */

/** partial sort an index array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an index array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of an array of pointers in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of an array of pointers in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/Reals/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrRealRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/Reals/SCIP_Bools/SCIP_Bools, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrRealRealBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/Reals/ints/SCIP_Bools, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrRealRealIntBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrRealRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/Reals/Reals/SCIP_Bools/SCIP_Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrRealRealBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/Reals/Reals/ints/SCIP_Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrRealRealIntBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/Reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrRealReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrRealReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/Bools, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an array of Reals in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Reals in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort array of ints in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort array of ints in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedInt(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntRealLong(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntRealLong(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an array of Longints in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Longints in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an index array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an index array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of an array of pointers in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of an array of pointers in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an array of Reals in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Reals in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< integer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< integer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );

/** partial sort of three joint arrays of Reals/Reals/Pointer/Pointer, sorted by first array in non-increasing order around the \p k-th element */
SCIP_EXPORT
void SCIPselectDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );

/** partial sort of three joint arrays of Reals/Reals/Pointer/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
SCIP_EXPORT
void SCIPselectWeightedDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );

/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort array of ints in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort array of ints in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownInt(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an array of Longints in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Longints in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT
void SCIPselectWeightedDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
