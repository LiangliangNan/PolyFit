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

/**@file   type_misc.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for miscellaneous datastructures
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_MISC_H__
#define __SCIP_TYPE_MISC_H__

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** represents different confidence levels for (one-sided) hypothesis testing; in order to obtain two-sided confidence
 *  levels, calculate 2 * c - 1, i.e., if the one-sided confidence level is 90 %, the two-sided level is 80 %
 */
enum SCIP_Confidencelevel
{
   SCIP_CONFIDENCELEVEL_MIN = 0,    /**< one-sided confidence level 75 %, two-sided 50 % */
   SCIP_CONFIDENCELEVEL_LOW = 1,    /**< (one-sided) confidence level 87.5 %, two-sided 75 % */
   SCIP_CONFIDENCELEVEL_MEDIUM = 2, /**< (one-sided) confidence level 90 %, two-sided 80 % */
   SCIP_CONFIDENCELEVEL_HIGH = 3,   /**< (one-sided) confidence level 95 %, two-sided 90 % */
   SCIP_CONFIDENCELEVEL_MAX = 4     /**< (one-sided) confidence level 97.5 %, two-sided 95 % */
};
typedef enum SCIP_Confidencelevel SCIP_CONFIDENCELEVEL;

/** type of hashmap: are pointers, reals or ints stored, or unknown */
enum SCIP_Hashmaptype
{
   SCIP_HASHMAPTYPE_UNKNOWN = 0,    /**< the hashmap did not store a single element yet, type unknown */
   SCIP_HASHMAPTYPE_POINTER = 1,    /**< the hashmap stores pointers % */
   SCIP_HASHMAPTYPE_REAL    = 2,    /**< the hashmap stores reals */
   SCIP_HASHMAPTYPE_INT     = 3     /**< the hashmap stores ints */
};
typedef enum SCIP_Hashmaptype SCIP_HASHMAPTYPE;

/** Sparse solution data structure
 *
 *  - \ref SparseSol "List of all available methods"
 */
typedef struct SCIP_SparseSol SCIP_SPARSESOL;

/** (circular) Queue data structure
 *
 *  - \ref Queue "List of all available methods"
 */
typedef struct SCIP_Queue SCIP_QUEUE;

/** Priority queue data structure
 *
 *  - \ref PriorityQueue "List of all available methods"
 */
typedef struct SCIP_PQueue SCIP_PQUEUE;

/** Hash table data structure
 *
 *  - \ref HashTable "List of all available methods"
 */
typedef struct SCIP_HashTable SCIP_HASHTABLE;

/** Hash table data structure which allows multiple occurences of an element
 *
 *  - \ref MultiHash "List of all available methods"
 */
typedef struct SCIP_MultiHash SCIP_MULTIHASH;

/** Hash table element list to store single elements of a multi hash table */
typedef struct SCIP_MultiHashList SCIP_MULTIHASHLIST;

/** Hash map entry */
typedef struct SCIP_HashMapEntry SCIP_HASHMAPENTRY;

/** Hash map data structure
 *
 *  - \ref HashMap "List of all available methods"
 */
typedef struct SCIP_HashMap SCIP_HASHMAP;

/** Hash set data structure
 *
 *  - \ref HashMap "List of all available methods"
 */
typedef struct SCIP_HashSet SCIP_HASHSET;

/** dynamic array for storing SCIP_Real values */
typedef struct SCIP_RealArray SCIP_REALARRAY;

/** dynamic array for storing int values */
typedef struct SCIP_IntArray SCIP_INTARRAY;

/** dynamic array for storing SCIP_Bool values */
typedef struct SCIP_BoolArray SCIP_BOOLARRAY;

/** dynamic array for storing pointers */
typedef struct SCIP_PtrArray SCIP_PTRARRAY;

/** random number generator */
typedef struct SCIP_RandNumGen SCIP_RANDNUMGEN;

/** Resource activity data structure
 *
 *  - \ref ResourceActivity "List of all available methods"
 */
typedef struct SCIP_ResourceActivity SCIP_RESOURCEACTIVITY;

/** Resource profile data structure
 *
 *  - \ref ResourceProfile "List of all available methods"
 */
typedef struct SCIP_Profile SCIP_PROFILE;

/** Directed graph data structure (stored as adjacency list)
 *
 *  - \ref DirectedGraph "List of all available methods"
 */
typedef struct SCIP_Digraph SCIP_DIGRAPH;

/** Binary tree data structure
 *
 *  - \ref BinaryTree "List of all available methods"
 */
typedef struct SCIP_Bt SCIP_BT;

/** search node of \ref SCIP_BT "binary tree" */
typedef struct SCIP_BtNode SCIP_BTNODE;

/** regression data structure to compute an incremental linear regression of paired observations
 *
 *  - \ref Regression "List of all available methods"
 */
typedef struct SCIP_Regression SCIP_REGRESSION;

/** disjoint set (disjoint set (union find)) data structure for querying and updating connectedness of a graph with integer vertices 0,...,n - 1
 *
 *  - \ref DisjointSet "List of available methods"
 */
typedef struct SCIP_DisjointSet SCIP_DISJOINTSET;

/** a linear inequality row in preparation to become a SCIP_ROW
 *
 * Used to assemble data that could eventually make a SCIP_ROW.
 * @note Only one-sided rows are allowed here.
 */
typedef struct SCIP_RowPrep SCIP_ROWPREP;

/** compares two element indices
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
#define SCIP_DECL_SORTINDCOMP(x) int x (void* dataptr, int ind1, int ind2)

/** compares two data element pointers
 *  result:
 *    < 0: elem1 comes before (is better than) elem2
 *    = 0: both elements have the same value
 *    > 0: elem2 comes after (is worse than) elem2
 */
#define SCIP_DECL_SORTPTRCOMP(x) int x (void* elem1, void* elem2)

/** gets the key of the given element */
#define SCIP_DECL_HASHGETKEY(x) void* x (void* userptr, void* elem)

/** returns TRUE iff both keys are equal */
#define SCIP_DECL_HASHKEYEQ(x) SCIP_Bool x (void* userptr, void* key1, void* key2)

/** returns the hash value of the key */
#define SCIP_DECL_HASHKEYVAL(x) uint64_t x (void* userptr, void* key)

/** evaluates the real function at the given point, used for Newton Procedure
 *  input:
 *      point: the point where the function is evaluated
 *      params: an array of parameters that the function depends on (e.g. f(x) = a*x + b)
 *      nparams: the number of parameters stored in params
 */
#define SCIP_DECL_NEWTONEVAL(x) SCIP_Real x (SCIP_Real point, SCIP_Real* params, int nparams)

/** callback to act on position change of @p elem in priority queue */
#define SCIP_DECL_PQUEUEELEMCHGPOS(x) void x (void* elem, int oldpos, int newpos)

#ifdef __cplusplus
}
#endif

#endif
