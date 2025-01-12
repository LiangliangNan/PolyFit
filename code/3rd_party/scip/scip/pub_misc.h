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

/**@file   pub_misc.h
 * @ingroup PUBLICCOREAPI
 * @brief  public data structures and miscellaneous methods
 * @author Tobias Achterberg
 * @author Gerald Gamrath
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * This file contains a bunch of data structures and miscellaneous methods:
 *
 * - \ref DataStructures "Data structures"
 * - \ref MiscellaneousMethods "Miscellaneous Methods"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MISC_H__
#define __SCIP_PUB_MISC_H__

/* on SunOS, the function finite(a) (for the SCIPisFinite macro below) is declared in ieeefp.h */
#ifdef __sun
#include <ieeefp.h>
#endif
#include <math.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_misc.h"
#include "scip/type_message.h"
#include "scip/type_var.h"
#include "scip/pub_misc_select.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_misc_linear.h"
#include "scip/pub_misc_rowprep.h"

/* in optimized mode some of the function are handled via defines, for that the structs are needed */
#ifdef NDEBUG
#include "scip/struct_misc.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * methods for statistical tests
 */

/**@defgroup STATISTICALTESTS Statistical tests
 * @ingroup MiscellaneousMethods
 * @brief public methods for statistical tests
 *
 * Below are the public methods for statistical tests inside of \SCIP
 *
 * @{
 */

/** get critical value of a Student-T distribution for a given number of degrees of freedom at a confidence level */
SCIP_EXPORT
SCIP_Real SCIPstudentTGetCriticalValue(
   SCIP_CONFIDENCELEVEL  clevel,             /**< (one-sided) confidence level */
   int                   df                  /**< degrees of freedom */
   );

/** compute a t-value for the hypothesis that x and y are from the same population; Assuming that
 *  x and y represent normally distributed random samples with equal variance, the returned value
 *  comes from a Student-T distribution with countx + county - 2 degrees of freedom; this
 *  value can be compared with a critical value (see also SCIPstudentTGetCriticalValue()) at
 *  a predefined confidence level for checking if x and y significantly differ in location
 */
SCIP_EXPORT
SCIP_Real SCIPcomputeTwoSampleTTestValue(
   SCIP_Real             meanx,              /**< the mean of the first distribution */
   SCIP_Real             meany,              /**< the mean of the second distribution */
   SCIP_Real             variancex,          /**< the variance of the x-distribution */
   SCIP_Real             variancey,          /**< the variance of the y-distribution */
   SCIP_Real             countx,             /**< number of samples of x */
   SCIP_Real             county              /**< number of samples of y */
   );

/** returns the value of the Gauss error function evaluated at a given point */
SCIP_EXPORT
SCIP_Real SCIPerf(
   SCIP_Real             x                   /**< value to evaluate */
   );

/** get critical value of a standard normal distribution  at a given confidence level */
SCIP_EXPORT
SCIP_Real SCIPnormalGetCriticalValue(
   SCIP_CONFIDENCELEVEL  clevel              /**< (one-sided) confidence level */
   );

/** calculates the cumulative distribution P(-infinity <= x <= value) that a normally distributed
 *  random variable x takes a value between -infinity and parameter \p value.
 *
 *  The distribution is given by the respective mean and deviation. This implementation
 *  uses the error function erf().
 */
SCIP_EXPORT
SCIP_Real SCIPnormalCDF(
   SCIP_Real             mean,               /**< the mean value of the distribution */
   SCIP_Real             variance,           /**< the square of the deviation of the distribution */
   SCIP_Real             value               /**< the upper limit of the calculated distribution integral */
   );

/**@} */

/**@defgroup Regression Linear Regression
 * @ingroup MiscellaneousMethods
 * @brief methods for linear regression
 *
 * Below are the public methods for incremental linear regression of observations pairs \f$(X_i,Y_i), i=1\dots,n\f$
 *
 * @{
 */

/** returns the number of observations of this regression */
SCIP_EXPORT
int SCIPregressionGetNObservations(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   );

/** return the current slope of the regression */
SCIP_EXPORT
SCIP_Real SCIPregressionGetSlope(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   );

/** get the current y-intercept of the regression */
SCIP_EXPORT
SCIP_Real SCIPregressionGetIntercept(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   );

/** removes an observation (x,y) from the regression */
SCIP_EXPORT
void SCIPregressionRemoveObservation(
   SCIP_REGRESSION*      regression,         /**< regression data structure */
   SCIP_Real             x,                  /**< X of observation */
   SCIP_Real             y                   /**< Y of the observation */
   );

/** update regression by a new observation (x,y) */
SCIP_EXPORT
void SCIPregressionAddObservation(
   SCIP_REGRESSION*      regression,         /**< regression data structure */
   SCIP_Real             x,                  /**< X of observation */
   SCIP_Real             y                   /**< Y of the observation */
   );

/** reset regression data structure */
SCIP_EXPORT
void SCIPregressionReset(
   SCIP_REGRESSION*      regression          /**< regression data structure */
   );

/** creates and resets a regression */
SCIP_EXPORT
SCIP_RETCODE SCIPregressionCreate(
   SCIP_REGRESSION**     regression          /**< regression data structure */
   );

/** frees a regression */
SCIP_EXPORT
void SCIPregressionFree(
   SCIP_REGRESSION**     regression          /**< regression data structure */
   );

/**@} */

/*
 */

/**@defgroup GMLgraph GML Graphical Printing
 * @ingroup MiscellaneousMethods
 * @brief GML graph printing methods
 *
 * For a detailed format decription see http://docs.yworks.com/yfiles/doc/developers-guide/gml.html
 *
 * @{
 */


/** writes a node section to the given graph file */
SCIP_EXPORT
void SCIPgmlWriteNode(
   FILE*                 file,               /**< file to write to */
   unsigned int          id,                 /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node, or NULL */
   const char*           fillcolor,          /**< color of the node's interior, or NULL */
   const char*           bordercolor         /**< color of the node's border, or NULL */
   );

/** writes a node section including weight to the given graph file */
SCIP_EXPORT
void SCIPgmlWriteNodeWeight(
   FILE*                 file,               /**< file to write to */
   unsigned int          id,                 /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node, or NULL */
   const char*           fillcolor,          /**< color of the node's interior, or NULL */
   const char*           bordercolor,        /**< color of the node's border, or NULL */
   SCIP_Real             weight              /**< weight of node */
   );

/** writes an edge section to the given graph file */
SCIP_EXPORT
void SCIPgmlWriteEdge(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   );

/** writes an arc section to the given graph file */
SCIP_EXPORT
void SCIPgmlWriteArc(
   FILE*                 file,               /**< file to write to */
   unsigned int          source,             /**< source node id of the node */
   unsigned int          target,             /**< target node id of the edge */
   const char*           label,              /**< label of the edge, or NULL */
   const char*           color               /**< color of the edge, or NULL */
   );

/** writes the starting line to a GML graph file, does not open a file */
SCIP_EXPORT
void SCIPgmlWriteOpening(
   FILE*                 file,               /**< file to write to */
   SCIP_Bool             directed            /**< is the graph directed */
   );

/** writes the ending lines to a GML graph file, does not close a file */
SCIP_EXPORT
void SCIPgmlWriteClosing(
   FILE*                 file                /**< file to close */
   );

/**@} */

/*
 * Sparse solution
 */

/**@defgroup SparseSol Sparse Solution
 * @ingroup DataStructures
 * @brief sparse storage for multiple integer solutions
 *
 * @{
 */

/** creates a sparse solution */
SCIP_EXPORT
SCIP_RETCODE SCIPsparseSolCreate(
   SCIP_SPARSESOL**      sparsesol,          /**< pointer to store the created sparse solution */
   SCIP_VAR**            vars,               /**< variables in the sparse solution, must not contain continuous variables */
   int                   nvars,              /**< number of variables to store, size of the lower and upper bound arrays */
   SCIP_Bool             cleared             /**< should the lower and upper bound arrays be cleared (entries set to 0) */
   );

/** frees sparse solution */
SCIP_EXPORT
void SCIPsparseSolFree(
   SCIP_SPARSESOL**      sparsesol           /**< pointer to a sparse solution */
   );

/** returns the variables in the given sparse solution */
SCIP_EXPORT
SCIP_VAR** SCIPsparseSolGetVars(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   );

/** returns the number of variables in the given sparse solution */
SCIP_EXPORT
int SCIPsparseSolGetNVars(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   );

/** returns the the lower bound array for all variables for a given sparse solution */
SCIP_EXPORT
SCIP_Longint* SCIPsparseSolGetLbs(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   );

/** returns the the upper bound array for all variables for a given sparse solution */
SCIP_EXPORT
SCIP_Longint* SCIPsparseSolGetUbs(
   SCIP_SPARSESOL*       sparsesol           /**< a sparse solution */
   );

/** constructs the first solution of sparse solution (all variables are set to their lower bound value */
SCIP_EXPORT
void SCIPsparseSolGetFirstSol(
   SCIP_SPARSESOL*       sparsesol,          /**< sparse solutions */
   SCIP_Longint*         sol,                /**< array to store the first solution */
   int                   nvars               /**< number of variables */
   );

/** constructs the next solution of the sparse solution and return whether there was one more or not */
SCIP_EXPORT
SCIP_Bool SCIPsparseSolGetNextSol(
   SCIP_SPARSESOL*       sparsesol,          /**< sparse solutions */
   SCIP_Longint*         sol,                /**< current solution array which get changed to the next solution */
   int                   nvars               /**< number of variables */
   );

/**@} */


/*
 * Queue
 */

/**@defgroup Queue Queue
 * @ingroup DataStructures
 * @brief circular FIFO queue
 *
 * @{
 */


/** creates a (circular) queue, best used if the size will be fixed or will not be increased that much */
SCIP_EXPORT
SCIP_RETCODE SCIPqueueCreate(
   SCIP_QUEUE**          queue,              /**< pointer to the new queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac             /**< memory growing factor applied, if more element slots are needed */
   );


/** frees queue, but not the data elements themselves */
SCIP_EXPORT
void SCIPqueueFree(
   SCIP_QUEUE**          queue               /**< pointer to a queue */
   );

/** clears the queue, but doesn't free the data elements themselves */
SCIP_EXPORT
void SCIPqueueClear(
   SCIP_QUEUE*           queue               /**< queue */
   );

/** inserts pointer element at the end of the queue */
SCIP_EXPORT
SCIP_RETCODE SCIPqueueInsert(
   SCIP_QUEUE*           queue,              /**< queue */
   void*                 elem                /**< element to be inserted */
   );

/** inserts unsigned integer element at the end of the queue */
SCIP_EXPORT
SCIP_RETCODE SCIPqueueInsertUInt(
   SCIP_QUEUE*           queue,              /**< queue */
   unsigned int          elem                /**< element to be inserted */
   );

/** removes and returns the first element of the queue, or NULL if no element exists */
SCIP_EXPORT
void* SCIPqueueRemove(
   SCIP_QUEUE*           queue               /**< queue */
   );

/** removes and returns the first unsigned integer element of the queue, or UNIT_MAX if no element exists */
SCIP_EXPORT
unsigned int SCIPqueueRemoveUInt(
   SCIP_QUEUE*           queue               /**< queue */
   );

/** returns the first element of the queue without removing it, or NULL if no element exists */
SCIP_EXPORT
void* SCIPqueueFirst(
   SCIP_QUEUE*           queue               /**< queue */
   );

/** returns the first unsigned integer element of the queue without removing it, or UINT_MAX if no element exists */
SCIP_EXPORT
unsigned int SCIPqueueFirstUInt(
   SCIP_QUEUE*           queue               /**< queue */
   );

/** returns whether the queue is empty */
SCIP_EXPORT
SCIP_Bool SCIPqueueIsEmpty(
   SCIP_QUEUE*           queue               /**< queue */
   );

/** returns the number of elements in the queue */
SCIP_EXPORT
int SCIPqueueNElems(
   SCIP_QUEUE*           queue               /**< queue */
   );

/**@} */

/*
 * Priority Queue
 */

/**@defgroup PriorityQueue Priority Queue
 * @ingroup DataStructures
 * @brief priority queue with O(1) access to the minimum element
 *
 * @{
 */

/** creates priority queue */
SCIP_EXPORT
SCIP_RETCODE SCIPpqueueCreate(
   SCIP_PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_DECL_PQUEUEELEMCHGPOS((*elemchgpos)) /**< callback to act on position change of elem in priority queue, or NULL */
   );

/** frees priority queue, but not the data elements themselves */
SCIP_EXPORT
void SCIPpqueueFree(
   SCIP_PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

/** clears the priority queue, but doesn't free the data elements themselves */
SCIP_EXPORT
void SCIPpqueueClear(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** inserts element into priority queue */
SCIP_EXPORT
SCIP_RETCODE SCIPpqueueInsert(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   );

/** delete element at specified position, maintaining the heap property */
SCIP_EXPORT
void SCIPpqueueDelPos(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   int                   pos                 /**< position of element that should be deleted */
   );

/** removes and returns best element from the priority queue */
SCIP_EXPORT
void* SCIPpqueueRemove(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the best element of the queue without removing it */
SCIP_EXPORT
void* SCIPpqueueFirst(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the number of elements in the queue */
SCIP_EXPORT
int SCIPpqueueNElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
SCIP_EXPORT
void** SCIPpqueueElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** return the position of @p elem in the priority queue, or -1 if element is not found */
SCIP_EXPORT
int SCIPpqueueFind(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   );

/**@} */


/*
 * Hash Table
 */

/**@defgroup HashTable Hash Table
 * @ingroup DataStructures
 * @brief hash table that resolves conflicts by probing
 *
 *@{
 */

/* fast 2-universal hash functions for two to seven 32bit elements with 32bit output */

#define SCIPhashSignature64(a)              (UINT64_C(0x8000000000000000)>>((UINT32_C(0x9e3779b9) * ((uint32_t)(a)))>>26))

#define SCIPhashTwo(a, b)                   ((uint32_t)((((uint32_t)(a) + 0xd37e9a1ce2148403ULL) * ((uint32_t)(b) + 0xe5fcc163aef32782ULL) )>>32))

#define SCIPhashThree(a, b, c)            ((uint32_t)((((uint32_t)(a) + 0xbd5c89185f082658ULL) * ((uint32_t)(b) + 0xe5fcc163aef32782ULL) + \
                                                        (uint32_t)(c) * 0xd37e9a1ce2148403ULL)>>32 ))

#define SCIPhashFour(a, b, c, d)            ((uint32_t)((((uint32_t)(a) + 0xbd5c89185f082658ULL) * ((uint32_t)(b) + 0xe5fcc163aef32782ULL) + \
                                                         ((uint32_t)(c) + 0xd37e9a1ce2148403ULL) * ((uint32_t)(d) + 0x926f2d4dc4a67218ULL))>>32 ))

#define SCIPhashFive(a, b, c, d, e)            ((uint32_t)((((uint32_t)(a) + 0xbd5c89185f082658ULL) * ((uint32_t)(b) + 0xe5fcc163aef32782ULL) + \
                                                         ((uint32_t)(c) + 0xd37e9a1ce2148403ULL) * ((uint32_t)(d) + 0x926f2d4dc4a67218ULL) + \
                                                           (uint32_t)(e) * 0xf48d4cd331e14327ULL)>>32 ))

#define SCIPhashSix(a, b, c, d, e, f)            ((uint32_t)((((uint32_t)(a) + 0xbd5c89185f082658ULL) * ((uint32_t)(b) + 0xe5fcc163aef32782ULL) + \
                                                         ((uint32_t)(c) + 0xd37e9a1ce2148403ULL) * ((uint32_t)(d) + 0x926f2d4dc4a67218ULL) + \
                                                         ((uint32_t)(e) + 0xf48d4cd331e14327ULL) * ((uint32_t)(f) + 0x80791a4edfc44c75ULL))>>32 ))

#define SCIPhashSeven(a, b, c, d, e, f, g)            ((uint32_t)((((uint32_t)(a) + 0xbd5c89185f082658ULL) * ((uint32_t)(b) + 0xe5fcc163aef32782ULL) + \
                                                         ((uint32_t)(c) + 0xd37e9a1ce2148403ULL) * ((uint32_t)(d) + 0x926f2d4dc4a67218ULL) + \
                                                         ((uint32_t)(e) + 0xf48d4cd331e14327ULL) * ((uint32_t)(f) + 0x80791a4edfc44c75ULL) + \
                                                         (uint32_t)(g) * 0x7f497d9ba3bd83c0ULL)>>32 ))

/** computes a hashcode for double precision floating point values containing
 *  15 significant bits, the sign and the exponent
 */
INLINE static
uint32_t SCIPrealHashCode(double x)
{
   int theexp;
   return (((uint32_t)(uint16_t)(int16_t)ldexp(frexp(x, &theexp), 15))<<16) | (uint32_t)(uint16_t)theexp;
}

/** creates a hash table */
SCIP_EXPORT
SCIP_RETCODE SCIPhashtableCreate(
   SCIP_HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr             /**< user pointer */
   );

/** frees the hash table */
SCIP_EXPORT
void SCIPhashtableFree(
   SCIP_HASHTABLE**      hashtable           /**< pointer to the hash table */
   );

/** removes all elements of the hash table
 *
 *  @note From a performance point of view you should not fill and clear a hash table too often since the clearing can
 *        be expensive. Clearing is done by looping over all buckets and removing the hash table lists one-by-one.
 *
 *  @deprecated Please use SCIPhashtableRemoveAll()
 */
SCIP_EXPORT
SCIP_DEPRECATED
void SCIPhashtableClear(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** inserts element in hash table (multiple inserts of same element override the previous entry) */
SCIP_EXPORT
SCIP_RETCODE SCIPhashtableInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   );

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
SCIP_EXPORT
SCIP_RETCODE SCIPhashtableSafeInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   );

/** retrieve element with key from hash table, returns NULL if not existing */
SCIP_EXPORT
void* SCIPhashtableRetrieve(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 key                 /**< key to retrieve */
   );

/** returns whether the given element exists in the table */
SCIP_EXPORT
SCIP_Bool SCIPhashtableExists(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to search in the table */
   );

/** removes element from the hash table, if it exists */
SCIP_EXPORT
SCIP_RETCODE SCIPhashtableRemove(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to remove from the table */
   );

/** removes all elements of the hash table */
SCIP_EXPORT
void SCIPhashtableRemoveAll(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** returns number of hash table elements */
SCIP_EXPORT
SCIP_Longint SCIPhashtableGetNElements(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** gives the number of entries in the internal arrays of a hash table */
SCIP_EXPORT
int SCIPhashtableGetNEntries(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** gives the element at the given index or NULL if entry at that index has no element */
SCIP_EXPORT
void* SCIPhashtableGetEntry(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   int                   entryidx            /**< index of hash table entry */
   );

/** returns the load of the given hash table in percentage */
SCIP_EXPORT
SCIP_Real SCIPhashtableGetLoad(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** prints statistics about hash table usage */
SCIP_EXPORT
void SCIPhashtablePrintStatistics(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/**@} */

/*
 * MultiHash Table
 */

/**@defgroup MultiHash Multi Hash table
 * @ingroup DataStructures
 * @brief hash table that resolves conflicts by queueing, thereby allowing for duplicate entries
 *
 *@{
 */

/** returns a reasonable hash table size (a prime number) that is at least as large as the specified value */
SCIP_EXPORT
int SCIPcalcMultihashSize(
   int                   minsize             /**< minimal size of the hash table */
   );

/** creates a multihash table */
SCIP_EXPORT
SCIP_RETCODE SCIPmultihashCreate(
   SCIP_MULTIHASH**      multihash,          /**< pointer to store the created multihash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store multihash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr             /**< user pointer */
   );

/** frees the multihash table */
SCIP_EXPORT
void SCIPmultihashFree(
   SCIP_MULTIHASH**      multihash           /**< pointer to the multihash table */
   );

/** inserts element in multihash table (multiple inserts of same element possible)
 *
 *  @note A pointer to a multihashlist returned by SCIPmultihashRetrieveNext() might get invalid when adding an element
 *        to the hash table, due to dynamic resizing.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPmultihashInsert(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 element             /**< element to insert into the table */
   );

/** inserts element in multihash table (multiple insertion of same element is checked and results in an error)
 *
 *  @note A pointer to a multihashlist returned by SCIPmultihashRetrieveNext() might get invalid when adding a new
 *        element to the multihash table, due to dynamic resizing.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPmultihashSafeInsert(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 element             /**< element to insert into the table */
   );

/** retrieve element with key from multihash table, returns NULL if not existing */
SCIP_EXPORT
void* SCIPmultihashRetrieve(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 key                 /**< key to retrieve */
   );

/** retrieve element with key from multihash table, returns NULL if not existing
 *  can be used to retrieve all entries with the same key (one-by-one)
 *
 *  @note The returned multimultihashlist pointer might get invalid when adding a new element to the multihash table.
 */
SCIP_EXPORT
void* SCIPmultihashRetrieveNext(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   SCIP_MULTIHASHLIST**  multihashlist,      /**< input: entry in hash table list from which to start searching, or NULL
                                              *   output: entry in hash table list corresponding to element after
                                              *           retrieved one, or NULL */
   void*                 key                 /**< key to retrieve */
   );

/** returns whether the given element exists in the multihash table */
SCIP_EXPORT
SCIP_Bool SCIPmultihashExists(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 element             /**< element to search in the table */
   );

/** removes element from the multihash table, if it exists */
SCIP_EXPORT
SCIP_RETCODE SCIPmultihashRemove(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   void*                 element             /**< element to remove from the table */
   );

/** removes all elements of the multihash table
 *
 *  @note From a performance point of view you should not fill and clear a hash table too often since the clearing can
 *        be expensive. Clearing is done by looping over all buckets and removing the hash table lists one-by-one.
 */
SCIP_EXPORT
void SCIPmultihashRemoveAll(
   SCIP_MULTIHASH*       multihash           /**< multihash table */
   );

/** returns number of multihash table elements */
SCIP_EXPORT
SCIP_Longint SCIPmultihashGetNElements(
   SCIP_MULTIHASH*       multihash           /**< multihash table */
   );

/** returns the load of the given multihash table in percentage */
SCIP_EXPORT
SCIP_Real SCIPmultihashGetLoad(
   SCIP_MULTIHASH*       multihash           /**< multihash table */
   );

/** prints statistics about multihash table usage */
SCIP_EXPORT
void SCIPmultihashPrintStatistics(
   SCIP_MULTIHASH*       multihash,          /**< multihash table */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** standard hash key comparator for string keys */
SCIP_EXPORT
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqString);

/** standard hashing function for string keys */
SCIP_EXPORT
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValString);

/** gets the element as the key */
SCIP_EXPORT
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyStandard);

/** returns TRUE iff both keys(pointer) are equal */
SCIP_EXPORT
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqPtr);

/** returns the hash value of the key */
SCIP_EXPORT
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValPtr);

/**@} */


/*
 * Hash Map
 */

/**@defgroup HashMap Hash Map
 * @ingroup DataStructures
 * @brief hash map to store key-value pairs (called \p origin and \p image)
 *
 * @{
 */

/** creates a hash map mapping pointers to pointers */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapCreate(
   SCIP_HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash map entries */
   int                   mapsize             /**< size of the hash map */
   );

/** frees the hash map */
SCIP_EXPORT
void SCIPhashmapFree(
   SCIP_HASHMAP**        hashmap             /**< pointer to the hash map */
   );

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapInsert(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   );

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapInsertInt(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   int                   image               /**< new image for origin */
   );

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapInsertReal(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   SCIP_Real             image               /**< new image for origin */
   );

/** retrieves image of given origin from the hash map, or NULL if no image exists */
SCIP_EXPORT
void* SCIPhashmapGetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   );

/** retrieves image of given origin from the hash map, or INT_MAX if no image exists */
SCIP_EXPORT
int SCIPhashmapGetImageInt(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   );

/** retrieves image of given origin from the hash map, or SCIP_INVALID if no image exists */
SCIP_EXPORT
SCIP_Real SCIPhashmapGetImageReal(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   );

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapSetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   );

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapSetImageInt(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   int                   image               /**< new image for origin */
   );

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapSetImageReal(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   SCIP_Real             image               /**< new image for origin */
   );

/** checks whether an image to the given origin exists in the hash map */
SCIP_EXPORT
SCIP_Bool SCIPhashmapExists(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to search for */
   );

/** removes origin->image pair from the hash map, if it exists */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapRemove(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to remove from the list */
   );

/** prints statistics about hash map usage */
SCIP_EXPORT
void SCIPhashmapPrintStatistics(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** indicates whether a hash map has no entries */
SCIP_EXPORT
SCIP_Bool SCIPhashmapIsEmpty(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   );

/** gives the number of elements in a hash map */
SCIP_EXPORT
int SCIPhashmapGetNElements(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   );

/** gives the number of entries in the internal arrays of a hash map */
SCIP_EXPORT
int SCIPhashmapGetNEntries(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   );

/** gives the hashmap entry at the given index or NULL if entry has no element */
SCIP_EXPORT
SCIP_HASHMAPENTRY* SCIPhashmapGetEntry(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   int                   entryidx            /**< index of hash map entry */
   );

/** gives the origin of the hashmap entry */
SCIP_EXPORT
void* SCIPhashmapEntryGetOrigin(
   SCIP_HASHMAPENTRY*    entry               /**< hash map entry */
   );

/** gives the image of the hashmap entry */
SCIP_EXPORT
void* SCIPhashmapEntryGetImage(
   SCIP_HASHMAPENTRY*    entry               /**< hash map entry */
   );

/** gives the image of the hashmap entry */
SCIP_EXPORT
int SCIPhashmapEntryGetImageInt(
   SCIP_HASHMAPENTRY*    entry               /**< hash map entry */
   );

/** gives the image of the hashmap entry */
SCIP_EXPORT
SCIP_Real SCIPhashmapEntryGetImageReal(
   SCIP_HASHMAPENTRY*    entry               /**< hash map entry */
   );

/** sets pointer image of a hashmap entry */
SCIP_EXPORT
void SCIPhashmapEntrySetImage(
   SCIP_HASHMAPENTRY*    entry,              /**< hash map entry */
   void*                 image               /**< new image */
   );

/** sets integer image of a hashmap entry */
SCIP_EXPORT
void SCIPhashmapEntrySetImageInt(
   SCIP_HASHMAPENTRY*    entry,              /**< hash map entry */
   int                   image               /**< new image */
   );

/** sets real image of a hashmap entry */
SCIP_EXPORT
void SCIPhashmapEntrySetImageReal(
   SCIP_HASHMAPENTRY*    entry,              /**< hash map entry */
   SCIP_Real             image               /**< new image */
   );

/** removes all entries in a hash map. */
SCIP_EXPORT
SCIP_RETCODE SCIPhashmapRemoveAll(
   SCIP_HASHMAP*         hashmap             /**< hash map */
   );

/**@} */


/*
 * Hash Set
 */

/**@defgroup HashSet Hash Set
 * @ingroup DataStructures
 * @brief very lightweight hash set of pointers
 *
 * @{
 */

/** creates a hash set of pointers */
SCIP_EXPORT
SCIP_RETCODE SCIPhashsetCreate(
   SCIP_HASHSET**        hashset,            /**< pointer to store the created hash set */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash set entries */
   int                   size                /**< initial size of the hash set; it is guaranteed that the set is not
                                              *   resized if at most that many elements are inserted */
   );

/** frees the hash set */
SCIP_EXPORT
void SCIPhashsetFree(
   SCIP_HASHSET**        hashset,            /**< pointer to the hash set */
   BMS_BLKMEM*           blkmem              /**< block memory used to store hash set entries */
   );

/** inserts new element into the hash set */
SCIP_EXPORT
SCIP_RETCODE SCIPhashsetInsert(
   SCIP_HASHSET*         hashset,            /**< hash set */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash set entries */
   void*                 element             /**< element to insert */
   );

/** checks whether an element exists in the hash set */
SCIP_EXPORT
SCIP_Bool SCIPhashsetExists(
   SCIP_HASHSET*         hashset,            /**< hash set */
   void*                 element             /**< element to search for */
   );

/** removes an element from the hash set, if it exists */
SCIP_EXPORT
SCIP_RETCODE SCIPhashsetRemove(
   SCIP_HASHSET*         hashset,            /**< hash set */
   void*                 element             /**< origin to remove from the list */
   );

/** prints statistics about hash set usage */
SCIP_EXPORT
void SCIPhashsetPrintStatistics(
   SCIP_HASHSET*         hashset,            /**< hash set */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** indicates whether a hash set has no entries */
SCIP_EXPORT
SCIP_Bool SCIPhashsetIsEmpty(
   SCIP_HASHSET*         hashset             /**< hash set */
   );

/** gives the number of elements in a hash set */
SCIP_EXPORT
int SCIPhashsetGetNElements(
   SCIP_HASHSET*         hashset             /**< hash set */
   );

/** gives the number of slots of a hash set */
SCIP_EXPORT
int SCIPhashsetGetNSlots(
   SCIP_HASHSET*         hashset             /**< hash set */
   );

/** gives the array of hash set slots; contains all elements in indetermined order and may contain NULL values */
SCIP_EXPORT
void** SCIPhashsetGetSlots(
   SCIP_HASHSET*         hashset             /**< hash set */
   );

/** removes all entries in a hash set. */
SCIP_EXPORT
void SCIPhashsetRemoveAll(
   SCIP_HASHSET*         hashset             /**< hash set */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPhashsetIsEmpty(hashset)        ((hashset)->nelements == 0)
#define SCIPhashsetGetNElements(hashset)   ((hashset)->nelements)
#define SCIPhashsetGetNSlots(hashset)      (1u << (64 - (hashset)->shift))
#define SCIPhashsetGetSlots(hashset)       ((hashset)->slots)

#endif

/**@} */


/*
 * Activity
 */

/**@defgroup ResourceActivity Resource Activity
 * @ingroup DataStructures
 * @brief ressource activity data structure
 *
 * @{
 */

/** create a resource activity */
SCIP_EXPORT
SCIP_RETCODE SCIPactivityCreate(
   SCIP_RESOURCEACTIVITY** activity,         /**< pointer to store the resource activity */
   SCIP_VAR*             var,                /**< start time variable of the activity */
   int                   duration,           /**< duration of the activity */
   int                   demand              /**< demand of the activity */
   );

/** frees a resource activity */
SCIP_EXPORT
void SCIPactivityFree(
   SCIP_RESOURCEACTIVITY** activity          /**< pointer to the resource activity */
   );

#ifndef NDEBUG

/** returns the start time variable of the resource activity */
SCIP_EXPORT
SCIP_VAR* SCIPactivityGetVar(
   SCIP_RESOURCEACTIVITY* activity           /**< resource activity */
   );

/** returns the duration of the resource activity */
SCIP_EXPORT
int SCIPactivityGetDuration(
   SCIP_RESOURCEACTIVITY* activity           /**< resource activity */
   );

/** returns the demand of the resource activity */
SCIP_EXPORT
int SCIPactivityGetDemand(
   SCIP_RESOURCEACTIVITY* activity           /**< resource activity */
   );

/** returns the energy of the resource activity */
SCIP_EXPORT
int SCIPactivityGetEnergy(
   SCIP_RESOURCEACTIVITY* activity           /**< resource activity */
   );

#else

#define SCIPactivityGetVar(activity)         ((activity)->var)
#define SCIPactivityGetDuration(activity)    ((activity)->duration)
#define SCIPactivityGetDemand(activity)      ((activity)->demand)
#define SCIPactivityGetEnergy(activity)      ((activity)->duration * (activity)->demand)

#endif

/**@} */


/*
 * Resource Profile
 */

/**@defgroup ResourceProfile Resource Profile
 * @ingroup DataStructures
 * @brief ressource profile data structure
 *
 * @{
 */

/** creates resource profile */
SCIP_EXPORT
SCIP_RETCODE SCIPprofileCreate(
   SCIP_PROFILE**        profile,            /**< pointer to store the resource profile */
   int                   capacity            /**< resource capacity */
   );

/** frees given resource profile */
SCIP_EXPORT
void SCIPprofileFree(
   SCIP_PROFILE**        profile             /**< pointer to the resource profile */
   );

/** output of the given resource profile */
SCIP_EXPORT
void SCIPprofilePrint(
   SCIP_PROFILE*         profile,            /**< resource profile to output */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** returns the capacity of the resource profile */
SCIP_EXPORT
int SCIPprofileGetCapacity(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the number time points of the resource profile */
SCIP_EXPORT
int SCIPprofileGetNTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the time points of the resource profile */
SCIP_EXPORT
int* SCIPprofileGetTimepoints(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the loads of the resource profile */
SCIP_EXPORT
int* SCIPprofileGetLoads(
   SCIP_PROFILE*         profile             /**< resource profile to use */
   );

/** returns the time point for given position of the resource profile */
SCIP_EXPORT
int SCIPprofileGetTime(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   pos                 /**< position */
   );

/** returns the loads of the resource profile at the given position */
SCIP_EXPORT
int SCIPprofileGetLoad(
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   pos                 /**< position */
   );

/** returns if the given time point exists in the resource profile and stores the position of the given time point if it
 *  exists; otherwise the position of the next smaller existing time point is stored
 */
SCIP_EXPORT
SCIP_Bool SCIPprofileFindLeft(
   SCIP_PROFILE*         profile,            /**< resource profile to search */
   int                   timepoint,          /**< time point to search for */
   int*                  pos                 /**< pointer to store the position */
   );

/** insert a core into resource profile; if the core is non-empty the resource profile will be updated otherwise nothing
 *  happens
 */
SCIP_EXPORT
SCIP_RETCODE SCIPprofileInsertCore(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height,             /**< height of the core */
   int*                  pos,                /**< pointer to store the first position were it gets infeasible */
   SCIP_Bool*            infeasible          /**< pointer to store if the core does not fit due to capacity */
   );

/** subtracts the height from the resource profile during core time */
SCIP_EXPORT
SCIP_RETCODE SCIPprofileDeleteCore(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height              /**< height of the core */
   );

/** return the earliest possible starting point within the time interval [lb,ub] for a given core (given by its height
 *  and duration)
 */
SCIP_EXPORT
int SCIPprofileGetEarliestFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   est,                /**< earliest starting time of the given core */
   int                   lst,                /**< latest starting time of the given core */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the corer cannot be inserted */
   );

/** return the latest possible starting point within the time interval [lb,ub] for a given core (given by its height and
 *  duration)
 */
SCIP_EXPORT
int SCIPprofileGetLatestFeasibleStart(
   SCIP_PROFILE*         profile,            /**< resource profile to use */
   int                   lb,                 /**< earliest possible start point */
   int                   ub,                 /**< latest possible start point */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the core cannot be inserted */
   );

/**@} */

/*
 * Directed graph
 */

/**@addtogroup DirectedGraph
 *
 * @{
 */

/** resize directed graph structure */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphResize(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   nnodes              /**< new number of nodes */
   );

/** sets the sizes of the successor lists for the nodes in a directed graph and allocates memory for the lists */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphSetSizes(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int*                  sizes               /**< sizes of the successor lists */
   );

/** frees given directed graph structure */
SCIP_EXPORT
void SCIPdigraphFree(
   SCIP_DIGRAPH**        digraph             /**< pointer to the directed graph */
   );

/** add (directed) arc and a related data to the directed graph structure
 *
 *  @note if the arc is already contained, it is added a second time
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphAddArc(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   );

/** add (directed) arc to the directed graph structure, if it is not contained, yet
 *
 * @note if there already exists an arc from startnode to endnode, the new arc is not added,
 *       even if its data is different
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphAddArcSafe(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the arc */
   int                   endnode,            /**< start node of the arc */
   void*                 data                /**< data that should be stored for the arc; or NULL */
   );

/** sets the number of successors to a given value */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphSetNSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node,               /**< node for which the number of successors has to be changed */
   int                   nsuccessors         /**< new number of successors */
   );

/** returns the number of nodes of the given digraph */
SCIP_EXPORT
int SCIPdigraphGetNNodes(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the node data, or NULL if no data exist */
SCIP_EXPORT
void* SCIPdigraphGetNodeData(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the node data is returned */
   );

/** sets the node data */
SCIP_EXPORT
void SCIPdigraphSetNodeData(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   void*                 dataptr,            /**< user node data pointer, or NULL */
   int                   node                /**< node for which the node data is returned */
   );

/** returns the total number of arcs in the given digraph */
SCIP_EXPORT
int SCIPdigraphGetNArcs(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the number of successor nodes of the given node */
SCIP_EXPORT
int SCIPdigraphGetNSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the number of outgoing arcs is returned */
   );

/** returns the array of indices of the successor nodes; this array must not be changed from outside */
SCIP_EXPORT
int* SCIPdigraphGetSuccessors(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the array of outgoing arcs is returned */
   );

/** returns the array of data corresponding to the arcs originating at the given node, or NULL if no data exist; this
 *  array must not be changed from outside
 */
SCIP_EXPORT
void** SCIPdigraphGetSuccessorsData(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the data corresponding to the outgoing arcs is returned */
   );

/** identifies the articulation points in a given directed graph
 *  uses the helper recursive function findArticulationPointsUtil
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphGetArticulationPoints(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int**                 articulations,      /**< array to store the sorted node indices of the computed articulation points, or NULL */
   int*                  narticulations      /**< number of the computed articulation points, or NULL */
   );

/** Compute undirected connected components on the given graph.
 *
 *  @note For each arc, its reverse is added, so the graph does not need to be the directed representation of an
 *        undirected graph.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphComputeUndirectedComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   minsize,            /**< all components with less nodes are ignored */
   int*                  components,         /**< array with as many slots as there are nodes in the directed graph
                                              *   to store for each node the component to which it belongs
                                              *   (components are numbered 0 to ncomponents - 1); or NULL, if components
                                              *   are accessed one-by-one using SCIPdigraphGetComponent() */
   int*                  ncomponents         /**< pointer to store the number of components; or NULL, if the
                                              *   number of components is accessed by SCIPdigraphGetNComponents() */
   );

/** Computes all strongly connected components of an undirected connected component with Tarjan's Algorithm.
 *  The resulting strongly connected components are sorted topologically (starting from the end of the
 *  strongcomponents array).
 *
 *  @note In general a topological sort of the strongly connected components is not unique.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphComputeDirectedComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   compidx,            /**< number of the undirected connected component */
   int*                  strongcomponents,   /**< array to store the strongly connected components
                                              *   (length >= size of the component) */
   int*                  strongcompstartidx, /**< array to store the start indices of the strongly connected
                                              *   components (length >= size of the component) */
   int*                  nstrongcomponents   /**< pointer to store the number of strongly connected
                                              *   components */
   );

/** Performes an (almost) topological sort on the undirected components of the given directed graph. The undirected
 *  components should be computed before using SCIPdigraphComputeUndirectedComponents().
 *
 *  @note In general a topological sort is not unique.  Note, that there might be directed cycles, that are randomly
 *        broken, which is the reason for having only almost topologically sorted arrays.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdigraphTopoSortComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the number of previously computed undirected components for the given directed graph */
SCIP_EXPORT
int SCIPdigraphGetNComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** Returns the previously computed undirected component of the given number for the given directed graph.
 *  If the components were sorted using SCIPdigraphTopoSortComponents(), the component is (almost) topologically sorted.
 */
SCIP_EXPORT
void SCIPdigraphGetComponent(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   compidx,            /**< number of the component to return */
   int**                 nodes,              /**< pointer to store the nodes in the component; or NULL, if not needed */
   int*                  nnodes              /**< pointer to store the number of nodes in the component;
                                              *   or NULL, if not needed */
   );

/** frees the component information for the given directed graph */
SCIP_EXPORT
void SCIPdigraphFreeComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** output of the given directed graph via the given message handler */
SCIP_EXPORT
void SCIPdigraphPrint(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** prints the given directed graph structure in GML format into the given file */
SCIP_EXPORT
void SCIPdigraphPrintGml(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   FILE*                 file                /**< file to write to */
   );


/** output of the given directed graph via the given message handler */
SCIP_EXPORT
void SCIPdigraphPrintComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/**@} */

/*
 * Binary search tree
 */

/**@defgroup BinaryTree Binary Search Tree
 * @ingroup DataStructures
 * @brief binary search tree data structure
 *@{
 */

/** creates a binary tree node with sorting value and user data */
SCIP_EXPORT
SCIP_RETCODE SCIPbtnodeCreate(
   SCIP_BT*              tree,               /**< binary search tree */
   SCIP_BTNODE**         node,               /**< pointer to store the created search node */
   void*                 dataptr             /**< user node data pointer, or NULL */
   );

/** frees the binary node including the rooted subtree
 *
 *  @note The user pointer (object) is not freed. If needed, it has to be done by the user.
 */
SCIP_EXPORT
void SCIPbtnodeFree(
   SCIP_BT*              tree,               /**< binary tree */
   SCIP_BTNODE**         node                /**< node to be freed */
   );

/** returns the user data pointer stored in that node */
SCIP_EXPORT
void* SCIPbtnodeGetData(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns the parent which can be NULL if the given node is the root */
SCIP_EXPORT
SCIP_BTNODE* SCIPbtnodeGetParent(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns left child which can be NULL if the given node is a leaf */
SCIP_EXPORT
SCIP_BTNODE* SCIPbtnodeGetLeftchild(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns right child which can be NULL if the given node is a leaf */
SCIP_EXPORT
SCIP_BTNODE* SCIPbtnodeGetRightchild(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns the sibling of the node or NULL if does not exist */
SCIP_EXPORT
SCIP_BTNODE* SCIPbtnodeGetSibling(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns whether the node is a root node */
SCIP_EXPORT
SCIP_Bool SCIPbtnodeIsRoot(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns whether the node is a leaf */
SCIP_EXPORT
SCIP_Bool SCIPbtnodeIsLeaf(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns TRUE if the given node is left child */
SCIP_EXPORT
SCIP_Bool SCIPbtnodeIsLeftchild(
   SCIP_BTNODE*          node                /**< node */
   );

/** returns TRUE if the given node is right child */
SCIP_EXPORT
SCIP_Bool SCIPbtnodeIsRightchild(
   SCIP_BTNODE*          node                /**< node */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPbtnodeGetData(node)               ((node)->dataptr)
#define SCIPbtnodeGetParent(node)             ((node)->parent)
#define SCIPbtnodeGetLeftchild(node)          ((node)->left)
#define SCIPbtnodeGetRightchild(node)         ((node)->right)
#define SCIPbtnodeGetSibling(node)            ((node)->parent == NULL ? NULL : \
                                               (node)->parent->left == (node) ? (node)->parent->right : (node)->parent->left)
#define SCIPbtnodeIsRoot(node)                ((node)->parent == NULL)
#define SCIPbtnodeIsLeaf(node)                ((node)->left == NULL && (node)->right == NULL)
#define SCIPbtnodeIsLeftchild(node)           ((node)->parent == NULL ? FALSE : (node)->parent->left == (node) ? TRUE : FALSE)
#define SCIPbtnodeIsRightchild(node)          ((node)->parent == NULL ? FALSE : (node)->parent->right == (node) ? TRUE : FALSE)

#endif

/** sets the give node data
 *
 *  @note The old user pointer is not freed.
 */
SCIP_EXPORT
void SCIPbtnodeSetData(
   SCIP_BTNODE*          node,               /**< node */
   void*                 dataptr             /**< node user data pointer */
   );

/** sets parent node
 *
 *  @note The old parent including the rooted subtree is not delete.
 */
SCIP_EXPORT
void SCIPbtnodeSetParent(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          parent              /**< new parent node, or NULL */
   );

/** sets left child
 *
 *  @note The old left child including the rooted subtree is not delete.
 */
SCIP_EXPORT
void SCIPbtnodeSetLeftchild(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          left                /**< new left child, or NULL */
   );

/** sets right child
 *
 *  @note The old right child including the rooted subtree is not delete.
 */
SCIP_EXPORT
void SCIPbtnodeSetRightchild(
   SCIP_BTNODE*          node,               /**< node */
   SCIP_BTNODE*          right               /**< new right child, or NULL */
   );

/** creates an binary tree */
SCIP_EXPORT
SCIP_RETCODE SCIPbtCreate(
   SCIP_BT**             tree,               /**< pointer to store the created binary tree */
   BMS_BLKMEM*           blkmem              /**< block memory used to create nodes */
   );

/** frees binary tree
 *
 *  @note The user pointers (object) of the search nodes are not freed. If needed, it has to be done by the user.
 */
SCIP_EXPORT
void SCIPbtFree(
   SCIP_BT**             tree                /**< pointer to binary tree */
   );

/** prints the binary tree in GML format into the given file */
SCIP_EXPORT
void SCIPbtPrintGml(
   SCIP_BT*              tree,               /**< binary tree */
   FILE*                 file                /**< file to write to */
   );

/** returns whether the binary tree is empty (has no nodes) */
SCIP_EXPORT
SCIP_Bool SCIPbtIsEmpty(
   SCIP_BT *             tree                /**< binary tree */
   );

/** returns the root node of the binary tree or NULL if the binary tree is empty */
SCIP_EXPORT
SCIP_BTNODE* SCIPbtGetRoot(
   SCIP_BT*              tree                /**< tree to be evaluated */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPbtIsEmpty(tree) (tree->root == NULL)
#define SCIPbtGetRoot(tree) (tree->root)

#endif

/** sets root node
 *
 *  @note The old root including the rooted subtree is not delete.
 */
SCIP_EXPORT
void SCIPbtSetRoot(
   SCIP_BT*              tree,               /**< tree to be evaluated */
   SCIP_BTNODE*          root                /**< new root, or NULL */
   );

/**@} */

/**@addtogroup DisjointSet
 *
 * @{
 */

/*
 * disjoint set data structure
 */

/** clears the disjoint set (union find) structure \p djset */
SCIP_EXPORT
void SCIPdisjointsetClear(
   SCIP_DISJOINTSET*     djset               /**< disjoint set (union find) data structure */
   );

/** finds and returns the component identifier of this \p element */
SCIP_EXPORT
int SCIPdisjointsetFind(
   SCIP_DISJOINTSET*     djset,              /**< disjoint set (union find) data structure */
   int                   element             /**< element to be found */
   );

/** merges the components containing the elements \p p and \p q */
SCIP_EXPORT
void SCIPdisjointsetUnion(
   SCIP_DISJOINTSET*     djset,              /**< disjoint set (union find) data structure */
   int                   p,                  /**< first element */
   int                   q,                  /**< second element */
   SCIP_Bool             forcerepofp         /**< force representative of p to be new representative */
   );

/** returns the number of independent components in this disjoint set (union find) data structure */
SCIP_EXPORT
int SCIPdisjointsetGetComponentCount(
   SCIP_DISJOINTSET*     djset               /**< disjoint set (union find) data structure */
   );

/** returns the size (number of nodes) of this disjoint set (union find) data structure */
SCIP_EXPORT
int SCIPdisjointsetGetSize(
   SCIP_DISJOINTSET*     djset               /**< disjoint set (union find) data structure */
   );

/** @} */

/*
 * Numerical methods
 */

/**@defgroup NumericalMethods Numerical Methods
 * @ingroup MiscellaneousMethods
 * @brief commonly used numerical methods
 *
 * @{
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
SCIP_EXPORT
SCIP_Real SCIPcalcMachineEpsilon(
   void
   );

/** returns the next representable value of from in the direction of to */
SCIP_EXPORT
SCIP_Real SCIPnextafter(
   SCIP_Real             from,               /**< value from which the next representable value should be returned */
   SCIP_Real             to                  /**< direction in which the next representable value should be returned */
   );

/** calculates the greatest common divisor of the two given values */
SCIP_EXPORT
SCIP_Longint SCIPcalcGreComDiv(
   SCIP_Longint          val1,               /**< first value of greatest common devisor calculation */
   SCIP_Longint          val2                /**< second value of greatest common devisor calculation */
   );

/** calculates the smallest common multiple of the two given values */
SCIP_EXPORT
SCIP_Longint SCIPcalcSmaComMul(
   SCIP_Longint          val1,               /**< first value of smallest common multiple calculation */
   SCIP_Longint          val2                /**< second value of smallest common multiple calculation */
   );

/** calculates a binomial coefficient n over m, choose m elements out of n, maximal value will be 33 over 16 (because
 *  the n=33 is the last line in the Pascal's triangle where each entry fits in a 4 byte value), an error occurs due to
 *  big numbers or an negative value m (and m < n) and -1 will be returned
 */
SCIP_EXPORT
SCIP_Longint SCIPcalcBinomCoef(
   int                   n,                  /**< number of different elements */
   int                   m                   /**< number to choose out of the above */
   );

/** calculates hash for floating-point number by using Fibonacci hashing */
SCIP_EXPORT
unsigned int SCIPcalcFibHash(
   SCIP_Real             v                   /**< number to hash */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcalcFibHash(v)   ((v) >= 0 ? ((unsigned long long)((v) * 2654435769)) % UINT_MAX : ((unsigned long long)(-(v) * 683565275)) % UINT_MAX )

#endif

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
SCIP_EXPORT
SCIP_Bool SCIPrealToRational(
   SCIP_Real             val,                /**< real value r to convert into rational number */
   SCIP_Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   SCIP_Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   );

/** tries to find a value, such that all given values, if scaled with this value become integral in relative allowed
 *  difference in between mindelta and maxdelta
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcalcIntegralScalar(
   SCIP_Real*            vals,               /**< values to scale */
   int                   nvals,              /**< number of values to scale */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   );

/** given a (usually very small) interval, tries to find a rational number with simple denominator (i.e. a small
 *  number, probably multiplied with powers of 10) out of this interval; returns TRUE iff a valid rational
 *  number inside the interval was found
 */
SCIP_EXPORT
SCIP_Bool SCIPfindSimpleRational(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed for resulting rational number */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   );

/** given a (usually very small) interval, selects a value inside this interval; it is tried to select a rational number
 *  with simple denominator (i.e. a small number, probably multiplied with powers of 10);
 *  if no valid rational number inside the interval was found, selects the central value of the interval
 */
SCIP_EXPORT
SCIP_Real SCIPselectSimpleValue(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom             /**< maximal denominator allowed for resulting rational number */
   );

/** Performs the Newton Procedure from a given starting point to compute a root of the given function with
 *  specified precision and maximum number of iterations. If the procedure fails, SCIP_INVALID is returned.
 */
SCIP_EXPORT
SCIP_Real SCIPcalcRootNewton(
   SCIP_DECL_NEWTONEVAL((*function)),       /**< pointer to function for which roots are computed */
   SCIP_DECL_NEWTONEVAL((*derivative)),      /**< pointer to derivative of above function */
   SCIP_Real*            params,            /**< parameters needed for function (can be NULL) */
   int                   nparams,           /**< number of parameters (can be 0) */
   SCIP_Real             x,                 /**< starting point */
   SCIP_Real             eps,               /**< tolerance */
   int                   k                  /**< iteration limit */
   );

/* The C99 standard defines the function (or macro) isfinite.
 * On MacOS X, isfinite is also available.
 * From the BSD world, there comes a function finite.
 * On SunOS, finite is also available.
 * In the MS compiler world, there is a function _finite.
 * As last resort, we check whether x == x does not hold, but this works only for NaN's, not for infinities!
 */
#if _XOPEN_SOURCE >= 600 || defined(_ISOC99_SOURCE) || _POSIX_C_SOURCE >= 200112L || defined(__APPLE__)
#define SCIPisFinite isfinite
#elif defined(_BSD_SOURCE) || defined(__sun)
#define SCIPisFinite finite
#elif defined(_MSC_VER)
#define SCIPisFinite _finite
#else
#define SCIPisFinite(x) ((x) == (x))
#endif

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
SCIP_EXPORT
SCIP_Real SCIPrelDiff(
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPrelDiff(val1, val2)         ( ((val1)-(val2))/(MAX3(1.0,REALABS(val1),REALABS(val2))) )

#endif

/** computes the gap from the primal and the dual bound */
SCIP_EXPORT
SCIP_Real SCIPcomputeGap(
   SCIP_Real             eps,                /**< the value treated as zero */
   SCIP_Real             inf,                /**< the value treated as infinity */
   SCIP_Real             primalbound,        /**< the primal bound */
   SCIP_Real             dualbound           /**< the dual bound */
   );

/**@} */


/*
 * Random Numbers
 */

/**@defgroup RandomNumbers Random Numbers
 * @ingroup MiscellaneousMethods
 * @brief structures and methods for pseudo random number generation
 *
 *@{
 */

/** returns a random integer between minrandval and maxrandval
 *
 *  @deprecated Please use SCIPrandomGetInt() to request a random integer.
 */
SCIP_EXPORT
SCIP_DEPRECATED
int SCIPgetRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   );


/** returns a random integer between minrandval and maxrandval */
SCIP_EXPORT
int SCIPrandomGetInt(
   SCIP_RANDNUMGEN*      randgen,            /**< random number generator data */
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval          /**< maximal value to return */
   );

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
SCIP_EXPORT
SCIP_RETCODE SCIPrandomGetSubset(
   SCIP_RANDNUMGEN*      randgen,            /**< random number generator */
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems           /**< number of elements that should be drawn and stored */
   );

/** returns a random real between minrandval and maxrandval */
SCIP_EXPORT
SCIP_Real SCIPrandomGetReal(
   SCIP_RANDNUMGEN*      randgen,            /**< random number generator data */
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval          /**< maximal value to return */
   );

/** returns a random real between minrandval and maxrandval
 *
 *  @deprecated Please use SCIPrandomGetReal() to request a random real.
 */
SCIP_EXPORT
SCIP_DEPRECATED
SCIP_Real SCIPgetRandomReal(
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   );

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 *
 *  @deprecated Please use SCIPrandomGetSubset()
 */
SCIP_EXPORT
SCIP_DEPRECATED
SCIP_RETCODE SCIPgetRandomSubset(
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems,          /**< number of elements that should be drawn and stored */
   unsigned int          randseed            /**< seed value for random generator */
   );

/**@} */

/*
 * Permutations / Shuffling
 */

/**@defgroup PermutationsShuffling Permutations Shuffling
 * @ingroup MiscellaneousMethods
 * @brief methods for shuffling arrays
 *
 * @{
 */

/** swaps two ints */
SCIP_EXPORT
void SCIPswapInts(
   int*                  value1,             /**< pointer to first integer */
   int*                  value2              /**< pointer to second integer */
   );

/** swaps two real values */
SCIP_EXPORT
void SCIPswapReals(
   SCIP_Real*            value1,             /**< pointer to first real value */
   SCIP_Real*            value2              /**< pointer to second real value */
);

/** swaps the addresses of two pointers */
SCIP_EXPORT
void SCIPswapPointers(
   void**                pointer1,           /**< first pointer */
   void**                pointer2            /**< second pointer */
   );

/** randomly shuffles parts of an integer array using the Fisher-Yates algorithm
 *
 *  @deprecated Please use SCIPrandomPermuteIntArray()
 */
SCIP_EXPORT
SCIP_DEPRECATED
void SCIPpermuteIntArray(
   int*                  array,              /**< array to be shuffled */
   int                   begin,              /**< first included index that should be subject to shuffling
                                              *   (0 for first array entry)
                                              */
   int                   end,                /**< first excluded index that should not be subject to shuffling
                                              *   (array size for last array entry)
                                              */
   unsigned int*         randseed            /**< seed value for the random generator */
   );

/** randomly shuffles parts of an integer array using the Fisher-Yates algorithm */
SCIP_EXPORT
void SCIPrandomPermuteIntArray(
   SCIP_RANDNUMGEN*      randgen,            /**< random number generator */
   int*                  array,              /**< array to be shuffled */
   int                   begin,              /**< first included index that should be subject to shuffling
                                              *   (0 for first array entry)
                                              */
   int                   end                 /**< first excluded index that should not be subject to shuffling
                                              *   (array size for last array entry)
                                              */
   );

/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
SCIP_EXPORT
void SCIPrandomPermuteArray(
   SCIP_RANDNUMGEN*      randgen,            /**< random number generator */
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first included index that should be subject to shuffling
                                              *   (0 for first array entry)
                                              */
   int                   end                 /**< first excluded index that should not be subject to shuffling
                                              *   (array size for last array entry)
                                              */
   );

/** randomly shuffles parts of an array using the Fisher-Yates algorithm
 *
 *  @deprecated Please use SCIPrandomPermuteArray()
 */
SCIP_EXPORT
SCIP_DEPRECATED
void SCIPpermuteArray(
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first included index that should be subject to shuffling
                                              *   (0 for first array entry)
                                              */
   int                   end,                /**< first excluded index that should not be subject to shuffling
                                              *   (array size for last array entry)
                                              */
   unsigned int*         randseed            /**< pointer to seed value for the random generator */
   );

/**@} */


/*
 * Arrays
 */

/**@defgroup Arrays Arrays
 * @ingroup MiscellaneousMethods
 * @brief miscellaneous methods for arrays
 *
 * @{
 */


/** computes set intersection (duplicates removed) of two integer arrays that are ordered ascendingly
 *
 * @deprecated Switch to SCIPcomputeArraysIntersectionInt().
 */
SCIP_DEPRECATED
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeArraysIntersection(
   int*                  array1,             /**< first array (in ascending order) */
   int                   narray1,            /**< number of entries of first array */
   int*                  array2,             /**< second array (in ascending order) */
   int                   narray2,            /**< number of entries of second array */
   int*                  intersectarray,     /**< intersection of array1 and array2
                                              *   (note: it is possible to use array1 for this input argument) */
   int*                  nintersectarray     /**< pointer to store number of entries of intersection array
                                              *   (note: it is possible to use narray1 for this input argument) */
   );

/** computes set intersection (duplicates removed) of two integer arrays that are ordered ascendingly */
SCIP_EXPORT
void SCIPcomputeArraysIntersectionInt(
   int*                  array1,             /**< first array (in ascending order) */
   int                   narray1,            /**< number of entries of first array */
   int*                  array2,             /**< second array (in ascending order) */
   int                   narray2,            /**< number of entries of second array */
   int*                  intersectarray,     /**< intersection of array1 and array2
                                              *   (note: it is possible to use array1 for this input argument) */
   int*                  nintersectarray     /**< pointer to store number of entries of intersection array
                                              *   (note: it is possible to use narray1 for this input argument) */
   );

/** computes set intersection (duplicates removed) of two void-pointer arrays that are ordered ascendingly */
SCIP_EXPORT
void SCIPcomputeArraysIntersectionPtr(
   void**                array1,             /**< first array (in ascending order) */
   int                   narray1,            /**< number of entries of first array */
   void**                array2,             /**< second array (in ascending order) */
   int                   narray2,            /**< number of entries of second array */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void**                intersectarray,     /**< intersection of array1 and array2
                                              *   (note: it is possible to use array1 for this input argument) */
   int*                  nintersectarray     /**< pointer to store number of entries of intersection array
                                              *   (note: it is possible to use narray1 for this input argument) */
);

/** computes set difference (duplicates removed) of two integer arrays that are ordered ascendingly
 *
 * @deprecated Switch to SCIPcomputeArraysSetminusInt().
 */
SCIP_DEPRECATED
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeArraysSetminus(
   int*                  array1,             /**< first array (in ascending order) */
   int                   narray1,            /**< number of entries of first array */
   int*                  array2,             /**< second array (in ascending order) */
   int                   narray2,            /**< number of entries of second array */
   int*                  setminusarray,      /**< array to store entries of array1 that are not an entry of array2
                                              *   (note: it is possible to use array1 for this input argument) */
   int*                  nsetminusarray      /**< pointer to store number of entries of setminus array
                                              *   (note: it is possible to use narray1 for this input argument) */
   );

/** computes set difference (duplicates removed) of two integer arrays that are ordered ascendingly */
SCIP_EXPORT
void SCIPcomputeArraysSetminusInt(
   int*                  array1,             /**< first array (in ascending order) */
   int                   narray1,            /**< number of entries of first array */
   int*                  array2,             /**< second array (in ascending order) */
   int                   narray2,            /**< number of entries of second array */
   int*                  setminusarray,      /**< array to store entries of array1 that are not an entry of array2
                                              *   (note: it is possible to use array1 for this input argument) */
   int*                  nsetminusarray      /**< pointer to store number of entries of setminus array
                                              *   (note: it is possible to use narray1 for this input argument) */
   );

/**@} */


/*
 * Strings
 */

/**@defgroup StringMethods String Methods
 * @ingroup MiscellaneousMethods
 * @brief commonly used methods for strings
 *
 *@{
 */

/** copies characters from 'src' to 'dest', copying is stopped when either the 'stop' character is reached or after
 *  'cnt' characters have been copied, whichever comes first.
 *
 *  @note undefined behaviuor on overlapping arrays
 */
SCIP_EXPORT
int SCIPmemccpy(
   char*                 dest,               /**< destination pointer to copy to */
   const char*           src,                /**< source pointer to copy to */
   char                  stop,               /**< character when found stop copying */
   unsigned int          cnt                 /**< maximal number of characters to copy too */
   );

/** prints an error message containing of the given string followed by a string describing the current system error;
 *  prefers to use the strerror_r method, which is threadsafe; on systems where this method does not exist,
 *  NO_STRERROR_R should be defined (see INSTALL), in this case, srerror is used which is not guaranteed to be
 *  threadsafe (on SUN-systems, it actually is)
 */
SCIP_EXPORT
void SCIPprintSysError(
   const char*           message             /**< first part of the error message, e.g. the filename */
   );

/** extracts tokens from strings - wrapper method for strtok_r() */
SCIP_EXPORT
char* SCIPstrtok(
   char*                 s,                  /**< string to parse */
   const char*           delim,              /**< delimiters for parsing */
   char**                ptrptr              /**< pointer to working char pointer - must stay the same while parsing */
   );

/** translates the given string into a string where symbols ", ', and spaces are escaped with a \ prefix */
SCIP_EXPORT
void SCIPescapeString(
   char*                 t,                  /**< target buffer to store escaped string */
   int                   bufsize,            /**< size of buffer t */
   const char*           s                   /**< string to transform into escaped string */
   );

/** safe version of snprintf */
SCIP_EXPORT
int SCIPsnprintf(
   char*                 t,                  /**< target string */
   int                   len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
   );

/** safe version of strncpy
 *
 *  Copies string in s to t using at most @a size-1 nonzero characters (strncpy copies size characters). It always adds
 *  a terminating zero char. Does not pad the remaining string with zero characters (unlike strncpy). Returns the number
 *  of copied nonzero characters, if the length of s is at most size - 1, and returns size otherwise. Thus, the original
 *  string was truncated if the return value is size.
 */
SCIP_EXPORT
int SCIPstrncpy(
   char*                 t,                  /**< target string */
   const char*           s,                  /**< source string */
   int                   size                /**< maximal size of t */
   );

/** extract the next token as a integer value if it is one; in case no value is parsed the endptr is set to @p str
 *
 *  @return Returns TRUE if a value could be extracted, otherwise FALSE
 */
SCIP_EXPORT
SCIP_Bool SCIPstrToIntValue(
   const char*           str,                /**< string to search */
   int*                  value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   );

/** extract the next token as a double value if it is one; in case a value is parsed the endptr is set to @p str
 *
 *  @return Returns TRUE if a value could be extracted, otherwise FALSE
 */
SCIP_EXPORT
SCIP_Bool SCIPstrToRealValue(
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   );

/** copies the first size characters between a start and end character of str into token, if no error occurred endptr
 *  will point to the position after the read part, otherwise it will point to @p str
 */
SCIP_EXPORT
void SCIPstrCopySection(
   const char*           str,                /**< string to search */
   char                  startchar,          /**< character which defines the beginning */
   char                  endchar,            /**< character which defines the ending */
   char*                 token,              /**< string to store the copy */
   int                   size,               /**< size of the token char array */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   );

/** checks whether a given string t appears at the beginning of the string s (up to spaces at beginning) */
SCIP_EXPORT
SCIP_Bool SCIPstrAtStart(
        const char*           s,                  /**< string to search in */
        const char*           t,                  /**< string to search for */
        size_t                tlen                /**< length of t */
);

/**@} */

/*
 * File methods
 */

/**@defgroup FileMethods File Methods
 * @ingroup MiscellaneousMethods
 * @brief commonly used file methods
 *
 * @{
 */

/** returns, whether the given file exists */
SCIP_EXPORT
SCIP_Bool SCIPfileExists(
   const char*           filename            /**< file name */
   );

/** splits filename into path, name, and extension */
SCIP_EXPORT
void SCIPsplitFilename(
   char*                 filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**                path,               /**< pointer to store path, or NULL if not needed */
   char**                name,               /**< pointer to store name, or NULL if not needed */
   char**                extension,          /**< pointer to store extension, or NULL if not needed */
   char**                compression         /**< pointer to store compression extension, or NULL if not needed */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
