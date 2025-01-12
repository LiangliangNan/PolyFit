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

/**@file   struct_misc.h
 * @ingroup INTERNALAPI
 * @brief  miscellaneous datastructures
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_MISC_H__
#define __SCIP_STRUCT_MISC_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure for sparse solutions */
struct SCIP_SparseSol
{
   SCIP_VAR**            vars;               /**< variables */
   SCIP_Longint*         lbvalues;           /**< array of lower bounds */
   SCIP_Longint*         ubvalues;           /**< array of upper bounds */
   int                   nvars;              /**< number of variables */
};

typedef union {
   void*                 ptr;                /**< pointer element */
   unsigned int          uinteger;           /**< unsigned integer element */
} SCIP_QUEUEELEMENT;

/** (circular) Queue data structure */
struct SCIP_Queue
{
   SCIP_Real             sizefac;            /**< memory growing factor */
   SCIP_QUEUEELEMENT*    slots;              /**< array of element slots */
   int                   firstfree;          /**< first free slot */
   int                   firstused;          /**< first used slot */
   int                   size;               /**< total number of available element slots */
};

/** priority queue data structure
 *  Elements are stored in an array, which grows dynamically in size as new elements are added to the queue.
 *  The ordering is done through a pointer comparison function.
 *  The array is organized as follows. The root element (that is the "best" element \f$ r \f$ with \f$ r \leq x \f$ for all \f$ x \f$ )
 *  is stored in position 0. The children of an element at position \f$ p \f$ are stored at positions \f$ q_1 = 2*p+1 \f$ and
 *  \f$ q_2 = 2*p+2 \f$ . That means, the parent of the element at position \f$ q \f$ is at position \f$ p = (q-1)/2 \f$ .
 *  At any time, the condition holds that \f$ p \leq q \f$ for each parent \f$ p \f$ and its children \f$ q \f$ .
 *  Insertion and removal of single elements needs time \f$ O(log n) \f$ .
 */
struct SCIP_PQueue
{
   SCIP_Real             sizefac;            /**< memory growing factor */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp));        /**< compares two data elements */
   SCIP_DECL_PQUEUEELEMCHGPOS((*elemchgpos));/**< callback to act on position change of elem in priority queue, or NULL */
   void**                slots;              /**< array of element slots */
   int                   len;                /**< number of used element slots */
   int                   size;               /**< total number of available element slots */
};

/** hash table data structure */
struct SCIP_HashTable
{
   SCIP_DECL_HASHGETKEY((*hashgetkey));      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq));       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval));      /**< returns the hash value of the key */
   BMS_BLKMEM*           blkmem;             /**< block memory used to store hash map entries */
   void*                 userptr;            /**< user pointer */
   void**                slots;              /**< slots of the hash table */
   uint32_t*             hashes;             /**< hash values of elements stored in slots */
   uint32_t              shift;              /**< power such that 2^(32-shift) == nslots */
   uint32_t              mask;               /**< mask used for fast modulo, i.e. nslots - 1 */
   uint32_t              nelements;          /**< number of elements in the hashtable */
};

/** element list to store single elements of a hash table */
struct SCIP_MultiHashList
{
   void*                 element;            /**< this element */
   SCIP_MULTIHASHLIST*   next;               /**< rest of the hash table list */
};

/** multihash table data structure */
struct SCIP_MultiHash
{
   SCIP_DECL_HASHGETKEY((*hashgetkey));      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq));       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval));      /**< returns the hash value of the key */
   BMS_BLKMEM*           blkmem;             /**< block memory used to store hash map entries */
   SCIP_MULTIHASHLIST**  lists;              /**< multihash table lists of the hash table */
   int                   nlists;             /**< number of lists stored in the hash table */
   void*                 userptr;            /**< user pointer */
   SCIP_Longint          nelements;          /**< number of elements in the hashtable */
};

typedef union {
   void*                 ptr;                /**< pointer image */
   int                   integer;            /**< integer image */
   SCIP_Real             real;               /**< real image */
} SCIP_HASHMAPIMAGE;

/** hash map entry */
struct SCIP_HashMapEntry
{
   void*                 origin;             /**< origin of element */
   SCIP_HASHMAPIMAGE     image;              /**< image of element */
};

/** hash map data structure to map pointers on pointers */
struct SCIP_HashMap
{
   BMS_BLKMEM*           blkmem;             /**< block memory used to store hash map entries */
   SCIP_HASHMAPENTRY*    slots;              /**< buffer for hashmap entries */
   uint32_t*             hashes;             /**< hashes of elements */
   uint32_t              shift;              /**< power such that 2^(32-shift) == nslots */
   uint32_t              mask;               /**< mask used for fast modulo, i.e. nslots - 1 */
   uint32_t              nelements;          /**< number of elements in the hashtable */
   SCIP_HASHMAPTYPE      hashmaptype;        /**< type of entries */
};

/** lightweight hash set data structure to map pointers on pointers */
struct SCIP_HashSet
{
   void**                slots;              /**< buffer for hashmap entries */
   uint32_t              shift;              /**< power such that 2^(64-shift) == nslots */
   uint32_t              nelements;          /**< number of elements in the hashtable */
};

/** dynamic array for storing real values */
struct SCIP_RealArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   SCIP_Real*            vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

/** dynamic array for storing int values */
struct SCIP_IntArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   int*                  vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

/** dynamic array for storing bool values */
struct SCIP_BoolArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   SCIP_Bool*            vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

/** dynamic array for storing pointers */
struct SCIP_PtrArray
{
   BMS_BLKMEM*           blkmem;             /**< block memory that stores the vals array */
   void**                vals;               /**< array values */
   int                   valssize;           /**< size of vals array */
   int                   firstidx;           /**< index of first element in vals array */
   int                   minusedidx;         /**< index of first non zero element in vals array */
   int                   maxusedidx;         /**< index of last non zero element in vals array */
};

/** resource activity */
struct SCIP_ResourceActivity
{
   SCIP_VAR*             var;                /**< start time variable of the activity */
   int                   duration;           /**< duration of the activity */
   int                   demand;             /**< demand of the activity */
};

/** resource profile */
struct SCIP_Profile
{
   int*                  timepoints;         /**< time point array */
   int*                  loads;              /**< array holding the load for each time point */
   int                   capacity;           /**< capacity of the resource profile */
   int                   ntimepoints;        /**< current number of entries */
   int                   arraysize;          /**< current array size */
};

/** digraph structure to store and handle graphs */
struct SCIP_Digraph
{
   BMS_BLKMEM*           blkmem;             /**< block memory pointer to store the data */
   int**                 successors;         /**< adjacency list: for each node (first dimension) list of all successors */
   void***               arcdata;            /**< arc data corresponding to the arcs to successors given by the successors array  */
   void**                nodedata;           /**< data for each node of graph */
   int*                  successorssize;     /**< sizes of the successor lists for the nodes */
   int*                  nsuccessors;        /**< number of successors stored in the adjacency lists of the nodes */
   int*                  components;         /**< array to store the node indices of the components, one component after the other */
   int*                  componentstarts;    /**< array to store the start indices of the components in the components array */
   int*                  articulations;      /**< array  of size narticulations to store the node indices of the articulation points */
   int                   ncomponents;        /**< number of undirected components stored */
   int                   componentstartsize; /**< size of array componentstarts */
   int                   nnodes;             /**< number of nodes, nodes should be numbered from 0 to nnodes-1 */
   int                   narticulations;     /**< number of articulation points among the graph nodes */
   SCIP_Bool             articulationscheck; /**< TRUE if the (computed) articulation nodes are up-to-date and FALSE otherwise */
};

/** binary node data structure for binary tree */
struct SCIP_BtNode
{
   SCIP_BTNODE*          parent;             /**< pointer to the parent node */
   SCIP_BTNODE*          left;               /**< pointer to the left child node */
   SCIP_BTNODE*          right;              /**< pointer to the right child node */
   void*                 dataptr;            /**< user pointer */
};

/** binary search tree data structure */
struct SCIP_Bt
{
   SCIP_BTNODE*          root;               /**< pointer to the dummy root node; root is left child */
   BMS_BLKMEM*           blkmem;             /**< block memory used to store tree nodes */
};

/** data structure for incremental linear regression of data points (X_i, Y_i)  */
struct SCIP_Regression
{
   SCIP_Real             intercept;          /**< the current axis intercept of the regression */
   SCIP_Real             slope;              /**< the current slope of the regression */
   SCIP_Real             meanx;              /**< mean of all X observations */
   SCIP_Real             meany;              /**< mean of all Y observations */
   SCIP_Real             sumxy;              /**< accumulated sum of all products X * Y */
   SCIP_Real             variancesumx;       /**< incremental variance term for X observations  */
   SCIP_Real             variancesumy;       /**< incremental variance term for Y observations */
   SCIP_Real             corrcoef;           /**< correlation coefficient of X and Y */
   int                   nobservations;      /**< number of observations so far */
};

/** random number generator data */
struct SCIP_RandNumGen
{
   uint32_t              seed;               /**< start seed */
   uint32_t              xor_seed;           /**< Xorshift seed */
   uint32_t              mwc_seed;           /**< Multiply-with-carry seed */
   uint32_t              cst_seed;           /**< constant seed */
};

/** disjoint set (disjoint set (union find)) data structure for querying and updating connectedness in a graph with integer vertices 0,...,n - 1 */
struct SCIP_DisjointSet
{
   int*                  parents;            /**< array to store the parent node index for every vertex */
   int*                  sizes;              /**< array to store the size of the subtree rooted at each vertex */
   int                   size;               /**< the number of vertices in the graph */
   int                   componentcount;     /**< counter for the number of connected components of the graph */
};

/** a linear inequality row in preparation to become a SCIP_ROW */
struct SCIP_RowPrep
{
   SCIP_VAR**            vars;               /**< variables */
   SCIP_Real*            coefs;              /**< coefficients of variables */
   int                   nvars;              /**< number of variables (= number of coefficients) */
   int                   varssize;           /**< length of variables array (= lengths of coefficients array) */
   SCIP_Real             side;               /**< side */
   SCIP_SIDETYPE         sidetype;           /**< type of side */
   SCIP_Bool             local;              /**< whether the row is only locally valid (i.e., for the current node) */
   char                  name[SCIP_MAXSTRLEN]; /**< row name */

   SCIP_Bool             recordmodifications;/**< whether to remember variables which coefficients were modified during cleanup */
   SCIP_VAR**            modifiedvars;       /**< variables which coefficient were modified by cleanup */
   int                   nmodifiedvars;      /**< number of variables which coefficient was modified */
   int                   modifiedvarssize;   /**< length of `modifiedvars` array */
   SCIP_Bool             modifiedside;       /**< whether the side was modified (relaxed) by cleanup */
};

#ifdef __cplusplus
}
#endif

#endif
