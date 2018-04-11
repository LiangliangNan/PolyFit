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

/**@file   pub_tree.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for branch and bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_TREE_H__
#define __SCIP_PUB_TREE_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"

#ifdef NDEBUG
#include "scip/struct_tree.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Node methods
 */

/**@addtogroup PublicNodeMethods
 *
 * @{
 */

/** node comparator for best lower bound */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPnodeCompLowerbound);

/** returns the set of variable branchings that were performed in the parent node to create this node */
EXTERN
void SCIPnodeGetParentBranchings(
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR**            branchvars,         /**< array of variables on which the branching has been performed in the parent node */
   SCIP_Real*            branchbounds,       /**< array of bounds which the branching in the parent node set */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which the branching in the parent node set */
   int*                  nbranchvars,        /**< number of variables on which branching has been performed in the parent node
                                              *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int                   branchvarssize      /**< available slots in arrays */
   );

/** returns the set of variable branchings that were performed in all ancestor nodes (nodes on the path to the root) to create this node */
EXTERN
void SCIPnodeGetAncestorBranchings(
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR**            branchvars,         /**< array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds,       /**< array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which the branchings in all ancestors set */
   int*                  nbranchvars,        /**< number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int                   branchvarssize      /**< available slots in arrays */
   );

/** returns the set of variable branchings that were performed between the given @p node and the given @p parent node. */
EXTERN
void SCIPnodeGetAncestorBranchingsPart(
   SCIP_NODE*            node,               /**< node data */
   SCIP_NODE*            parent,             /**< node data */
   SCIP_VAR**            branchvars,         /**< array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds,       /**< array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which the branchings in all ancestors set */
   int*                  nbranchvars,        /**< number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int                   branchvarssize      /**< available slots in arrays */
   );

/** outputs the path into given file stream in GML format */
EXTERN
SCIP_RETCODE SCIPnodePrintAncestorBranchings(
   SCIP_NODE*            node,               /**< node data */
   FILE*                 file                /**< file to output the path */
   );

/** returns the set of variable branchings that were performed in all ancestor nodes (nodes on the path to the root) to create this node
 *  sorted by the nodes, starting from the current node going up to the root
 */
EXTERN
void SCIPnodeGetAncestorBranchingPath(
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR**            branchvars,         /**< array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds,       /**< array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which the branchings in all ancestors set */
   int*                  nbranchvars,        /**< number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   branchvarssize,     /**< available slots in arrays */
   int*                  nodeswitches,       /**< marks, where in the arrays the branching decisions of the next node on the path
                                              *   start; branchings performed at the parent of node always start at position 0.
                                              *   For single variable branching, nodeswitches[i] = i holds */
   int*                  nnodes,             /**< number of nodes in the nodeswitch array */
   int                   nodeswitchsize      /**< available slots in node switch array */
   );


/** checks for two nodes whether they share the same root path, i.e., whether one is an ancestor of the other */
EXTERN
SCIP_Bool SCIPnodesSharePath(
   SCIP_NODE*            node1,              /**< node data */
   SCIP_NODE*            node2               /**< node data */
   );

/** finds the common ancestor node of two given nodes */
EXTERN
SCIP_NODE* SCIPnodesGetCommonAncestor(
   SCIP_NODE*            node1,              /**< node data */
   SCIP_NODE*            node2               /**< node data */
   );


/** gets the type of the node */
EXTERN
SCIP_NODETYPE SCIPnodeGetType(
   SCIP_NODE*            node                /**< node */
   );

/** gets successively assigned number of the node */
EXTERN
SCIP_Longint SCIPnodeGetNumber(
   SCIP_NODE*            node                /**< node */
   );

/** gets the depth of the node */
EXTERN
int SCIPnodeGetDepth(
   SCIP_NODE*            node                /**< node */
   );

/** gets the lower bound of the node */
EXTERN
SCIP_Real SCIPnodeGetLowerbound(
   SCIP_NODE*            node                /**< node */
   );

/** gets the estimated value of the best feasible solution in subtree of the node */
EXTERN
SCIP_Real SCIPnodeGetEstimate(
   SCIP_NODE*            node                /**< node */
   );


/** gets the reoptimization type of a node */
EXTERN
SCIP_REOPTTYPE SCIPnodeGetReopttype(
   SCIP_NODE*            node                /**< node */
   );

/** gets the unique id to identify the node during reoptimization; id is 0 if the node is the root or not part of the
 * reoptimization tree
 */
EXTERN
unsigned int SCIPnodeGetReoptID(
   SCIP_NODE*            node                /**< node */
   );

/** sets the reoptimization type of the node */
EXTERN
void SCIPnodeSetReopttype(
   SCIP_NODE*            node,               /**< node */
   SCIP_REOPTTYPE        reopttype           /**< reoptimization type */
   );

/** sets a unique id to identify the node during reoptimization */
EXTERN
void SCIPnodeSetReoptID(
   SCIP_NODE*            node,               /**< node */
   unsigned int          id                  /**< unique id */
   );

/** counts the number of bound changes due to branching, constraint propagation, and propagation */
EXTERN
void SCIPnodeGetNDomchg(
   SCIP_NODE*            node,               /**< node */
   int*                  nbranchings,        /**< pointer to store number of branchings (or NULL if not needed) */
   int*                  nconsprop,          /**< pointer to store number of constraint propagations (or NULL if not needed) */
   int*                  nprop               /**< pointer to store number of propagations (or NULL if not needed) */
   );

/** gets the domain change information of the node, i.e., the information about the differences in the
 *  variables domains to the parent node
 */
EXTERN
SCIP_DOMCHG* SCIPnodeGetDomchg(
   SCIP_NODE*            node                /**< node */
   );

/** gets the parent node of a node in the branch-and-bound tree, if any */
EXTERN
SCIP_NODE* SCIPnodeGetParent(
   SCIP_NODE*            node                /**< node */
   );

/** returns whether node is in the path to the current node */
EXTERN
SCIP_Bool SCIPnodeIsActive(
   SCIP_NODE*            node                /**< node */
   );

/** returns whether the node is marked to be propagated again */
EXTERN
SCIP_Bool SCIPnodeIsPropagatedAgain(
   SCIP_NODE*            node                /**< node data */
   );

/* returns the set of changed constraints for a particular node */
EXTERN
SCIP_CONSSETCHG* SCIPnodeGetConssetchg(
   SCIP_NODE*            node                /**< node data */
   );


#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPnodeGetType(node)           ((SCIP_NODETYPE)(node)->nodetype)
#define SCIPnodeGetNumber(node)         ((node)->number)
#define SCIPnodeGetDepth(node)          ((int) (node)->depth)
#define SCIPnodeGetLowerbound(node)     ((node)->lowerbound)
#define SCIPnodeGetEstimate(node)       ((node)->estimate)
#define SCIPnodeGetDomchg(node)         ((node)->domchg)
#define SCIPnodeGetParent(node)         ((node)->parent)
#define SCIPnodeIsActive(node)          ((node)->active)
#define SCIPnodeIsPropagatedAgain(node) ((node)->reprop)
#define SCIPnodeGetConssetchg(node)    ((node)->conssetchg)

#endif

/* @} */

#ifdef __cplusplus
}
#endif

#endif
