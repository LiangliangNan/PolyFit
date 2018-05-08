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

/**@file   rbtree.h
 * @brief  intrusive red black tree datastructure
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RB_TREE_H__
#define __SCIP_RB_TREE_H__

#include "scip/def.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_RBTreeNode SCIP_RBTREENODE;

struct SCIP_RBTreeNode
{
   uintptr_t             parent;
   SCIP_RBTREENODE*      child[2];
};

/* macro to make any structure a node in a rbtree.
 * It is very important that this is the first member of the structure.
 *
 * Example usage:
 * struct SomeStruct
 * {
 *    SCIP_RBTREE_HOOKS;
 *    OTHERDATA*   mydata;
 *    int          othermember;
 * };
 *
 */
#define SCIP_RBTREE_HOOKS SCIP_RBTREENODE _rbtreenode

/* convenience macros that automtically cast the given arguments to SCIP_RBTREENODE */
#define SCIPrbtreeFirst(root)  SCIPrbtreeFirst_call((SCIP_RBTREENODE*)(root))
#define SCIPrbtreeLast(root)  SCIPrbtreeLast_call((SCIP_RBTREENODE*)(root))
#define SCIPrbtreeSuccessor(x)  SCIPrbtreeSuccessor_call((SCIP_RBTREENODE*)(x))
#define SCIPrbtreePredecessor(x)  SCIPrbtreePredecessor_call((SCIP_RBTREENODE*)(x))
#define SCIPrbtreeDelete(root, node)     SCIPrbtreeDelete_call((SCIP_RBTREENODE**)(root), (SCIP_RBTREENODE*)(node))
#define SCIPrbtreeInsert(r,p,c,n)  SCIPrbtreeInsert_call((SCIP_RBTREENODE**)(r), (SCIP_RBTREENODE*)(p), (c), (SCIP_RBTREENODE*)(n) )
#define SCIPrbtreeFindInt(r,k,n)  SCIPrbtreeFindInt_call((SCIP_RBTREENODE*)(r),(k),(SCIP_RBTREENODE**)(n))
#define SCIPrbtreeFindReal(r,k,n)  SCIPrbtreeFindReal_call((SCIP_RBTREENODE*)(r),(k),(SCIP_RBTREENODE**)(n))
#define SCIPrbtreeFindPtr(c,r,k,n)  SCIPrbtreeFindPtr_call((c),(SCIP_RBTREENODE*)(r),(void*)(k),(SCIP_RBTREENODE**)(n))
#define SCIPrbtreeFindElem(c,r,k,n)  SCIPrbtreeFindElem_call((c),(SCIP_RBTREENODE*)(r),(SCIP_RBTREENODE*)(k),(SCIP_RBTREENODE**)(n))


#define FOR_EACH_NODE(type, n, r, body) \
   { \
     type n; \
     type __next;  \
     n = (type) SCIPrbtreeFirst(r); \
     while( n != NULL ) \
     { \
        __next = (type) SCIPrbtreeSuccessor(n); \
        body \
        n = __next; \
     } \
   }

#define SCIP_DEF_RBTREE_FIND(NAME, KEYTYPE, NODETYPE, LT, GT) \
   /** Searches for an element in the tree given by it's root. If a node equal to the given element exists in the tree, \
    *  (*node) will point to that node upon termination of this function and 0 will be returned. \
    *  If the tree is empty (*node) will be NULL. Otherwise (*node) will point to the predecessor or \
    *  successor of the given element and -1 or 1 will be returned respectively. The return value and the \
    *  predecessor or successor can then be passed to the insert function to insert the element but only if \
    *  it is not in the tree already, i.e. this function did not return 0. \
    * \
    *  @returns 0 if the key was found, then *node is the node with the key and \
    *           -1 or 1 if the node was node found, then *node is the node with the \
    *           next smaller, or next larger key repectively. \
    */ \
   int NAME( \
      NODETYPE*             root,            /**< root of the tree */ \
      KEYTYPE               key,             /**< key to search for */ \
      NODETYPE**            node             /**< pointer to return node */ \
      ) \
   { \
      SCIP_RBTREENODE* x; \
      *node = NULL; \
      x = (SCIP_RBTREENODE*) root;  \
      while( x != NULL ) \
      { \
         *node = (NODETYPE*) x; \
         if( LT(key, ((NODETYPE*)x)) ) \
            x = x->child[0]; \
         else if( GT(key, ((NODETYPE*)x)) ) \
            x = x->child[1]; \
         else \
            return 0; \
      } \
      if( *node != NULL && LT(key, ((NODETYPE*)(*node)) ) ) \
         return 1; \
      return -1; \
   }


/** get first element in tree with respect to sorting key */
EXTERN
SCIP_RBTREENODE* SCIPrbtreeFirst_call(
   SCIP_RBTREENODE*      root                /**< root of the tree */
   );

/** get last element in tree with respect to sorting key */
EXTERN
SCIP_RBTREENODE* SCIPrbtreeLast_call(
   SCIP_RBTREENODE*      root                /**< root of the tree */
   );

/** get successor of given element in the tree */
EXTERN
SCIP_RBTREENODE* SCIPrbtreeSuccessor_call(
   SCIP_RBTREENODE*      x                   /**< element to get successor for */
   );

/** get predecessor of given element in the tree */
EXTERN
SCIP_RBTREENODE* SCIPrbtreePredecessor_call(
   SCIP_RBTREENODE*      x                   /**< element to get predecessor for */
   );

/** delete the given node from the tree given by it's root node.
 *  The node must be contained in the tree rooted at root.
 */
EXTERN
void SCIPrbtreeDelete_call(
   SCIP_RBTREENODE**     root,               /**< root of the tree */
   SCIP_RBTREENODE*      node                /**< node to delete from tree */
   );

/** insert node into the tree given by it's root. Requires the
 *  future parent and the position of the parent as returned by
 *  the tree's find function defined using the SCIP_DEF_RBTREE_FIND
 *  macro.
 */
EXTERN
void SCIPrbtreeInsert_call(
   SCIP_RBTREENODE**     root,               /**< root of the tree */
   SCIP_RBTREENODE*      parent,             /**< future parent of node to be inserted */
   int                   pos,                /**< position of parent with respect to node, i.e.
                                              *   -1 if the parent's key comes before node and 1
                                              *   if the parent's key comes after the node
                                              */
   SCIP_RBTREENODE*      node                /**< node to insert into the tree */
   );

#ifdef __cplusplus
}
#endif

#endif
