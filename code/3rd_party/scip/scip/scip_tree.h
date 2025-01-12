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

/**@file   scip_tree.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for the branch-and-bound tree
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_TREE_H__
#define __SCIP_SCIP_TREE_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicTreeMethods
 *
 * @{
 */

/** gets focus node in the tree
 *
 *  if we are in probing/diving mode this method returns the node in the tree where the probing/diving mode was started.
 *
 *  @return the current node of the search tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetFocusNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current node in the tree
 *
 *  @return the current node of the search tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetCurrentNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets depth of current node, or -1 if no current node exists; in probing, the current node is the last probing node,
 *  such that the depth includes the probing path
 *
 *  @return the depth of current node, or -1 if no current node exists; in probing, the current node is the last probing node,
 *  such that the depth includes the probing path
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_EXPORT
int SCIPgetDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets depth of the focus node, or -1 if no focus node exists; the focus node is the currently processed node in the
 *  branching tree, excluding the nodes of the probing path
 *
 *  @return the depth of the focus node, or -1 if no focus node exists; the focus node is the currently processed node in the
 *  branching tree, excluding the nodes of the probing path
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_EXPORT
int SCIPgetFocusDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets current plunging depth (successive times, a child was selected as next node)
 *
 *  @return the current plunging depth (successive times, a child was selected as next node)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
int SCIPgetPlungeDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the root node of the tree
 *
 *  @return the root node of the search tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the effective root depth, i.e., the depth of the deepest node which is part of all paths from the root node
 *  to the unprocessed nodes.
 *
 *  @return effective root depth
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
int SCIPgetEffectiveRootDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether the current node is already solved and only propagated again
 *
 *  @return TRUE is returned if \SCIP performance repropagation, otherwise FALSE.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_Bool SCIPinRepropagation(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets children of focus node along with the number of children
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          children,           /**< pointer to store children array, or NULL if not needed */
   int*                  nchildren           /**< pointer to store number of children, or NULL if not needed */
   );

/** gets number of children of focus node
 *
 *  @return number of children of the focus node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
int SCIPgetNChildren(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets siblings of focus node along with the number of siblings
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetSiblings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          siblings,           /**< pointer to store siblings array, or NULL if not needed */
   int*                  nsiblings           /**< pointer to store number of siblings, or NULL if not needed */
   );

/** gets number of siblings of focus node
 *
 *  @return the number of siblings of focus node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
int SCIPgetNSiblings(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets leaves of the tree along with the number of leaves
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetLeaves(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          leaves,             /**< pointer to store leaves array, or NULL if not needed */
   int*                  nleaves             /**< pointer to store number of leaves, or NULL if not needed */
   );

/** gets number of leaves in the tree
 *
 *  @return the number of leaves in the tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
int SCIPgetNLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets number of nodes left in the tree (children + siblings + leaves)
 *
 *  @return the number of nodes left in the tree (children + siblings + leaves)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
int SCIPgetNNodesLeft(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @return the best child of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetPrioChild(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @return the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetPrioSibling(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best child of the focus node w.r.t. the node selection strategy
 *
 *  @return the best child of the focus node w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetBestChild(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best sibling of the focus node w.r.t. the node selection strategy
 *
 *  @return the best sibling of the focus node w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetBestSibling(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best leaf from the node queue w.r.t. the node selection strategy
 *
 *  @return the best leaf from the node queue w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetBestLeaf(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy
 *
 *  @return the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetBestNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf)
 *
 *  @return the node with smallest lower bound from the tree (child, sibling, or leaf)
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_NODE* SCIPgetBestboundNode(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** access to all data of open nodes (leaves, children, and siblings)
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetOpenNodesData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          leaves,             /**< pointer to store the leaves, or NULL if not needed */
   SCIP_NODE***          children,           /**< pointer to store the children, or NULL if not needed */
   SCIP_NODE***          siblings,           /**< pointer to store the siblings, or NULL if not needed */
   int*                  nleaves,            /**< pointer to store the number of leaves, or NULL */
   int*                  nchildren,          /**< pointer to store the number of children, or NULL */
   int*                  nsiblings           /**< pointer to store the number of siblings, or NULL */
   );

/** marks node and whole sub tree to be cut off from branch and bound tree
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcutoffNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be cut off */
   );

/** removes all nodes from branch and bound tree that were marked to be cut off via SCIPcutoffNode()
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note In diving mode, the removal of nodes is delayed until diving ends.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPpruneTree(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** marks the given node to be propagated again the next time a node of its subtree is processed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPrepropagateNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be propagated again */
   );

/** returns depth of first node in active path that is marked being cutoff
 *
 *  @return depth of first node in active path that is marked being cutoff
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
int SCIPgetCutoffdepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns depth of first node in active path that has to be propagated again
 *
 *  @return depth of first node in active path that has to be propagated again
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
int SCIPgetRepropdepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints all branching decisions on variables from the root to the given node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPprintNodeRootPath(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** sets whether the LP should be solved at the focus node
 *
 *  @note In order to have an effect, this method needs to be called after a node is focused but before the LP is
 *        solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
void SCIPsetFocusnodeLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             solvelp             /**< should the LP be solved? */
   );


/** query if node was the last parent of a branching of the tree
 *
 *  @return TRUE if node was the last parent of a branching of the tree
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
SCIP_Bool SCIPwasNodeLastBranchParent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< tree node, or NULL to check focus node */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
