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

/**@file   scip_tree.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for the branch-and-bound tree
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/debug.h"
#include "scip/nodesel.h"
#include "scip/pub_message.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_mem.h"
#include "scip/scip_tree.h"
#include "scip/scip_numerics.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/struct_tree.h"
#include "scip/tree.h"

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
SCIP_NODE* SCIPgetFocusNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetFocusNode", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetFocusNode(scip->tree);
}

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
SCIP_NODE* SCIPgetCurrentNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetCurrentNode", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetCurrentNode(scip->tree);
}

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
SCIP_NODE* SCIPgetRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRootNode", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetRootNode(scip->tree);
}

/** gets the effective root depth, i.e., the depth of the deepest node which is part of all paths from the root node
 *  to the unprocessed nodes.
 *
 *  @return effective root depth
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetEffectiveRootDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetEffectiveRootDepth", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetEffectiveRootDepth(scip->tree);
}

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
SCIP_Bool SCIPinRepropagation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPinRepropagation", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeInRepropagation(scip->tree);
}

/** gets children of focus node along with the number of children
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          children,           /**< pointer to store children array, or NULL if not needed */
   int*                  nchildren           /**< pointer to store number of children, or NULL if not needed */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetChildren", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( children != NULL )
      *children = scip->tree->children;
   if( nchildren != NULL )
      *nchildren = scip->tree->nchildren;

   return SCIP_OKAY;
}

/** gets number of children of focus node
 *
 *  @return number of children of the focus node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNChildren(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNChildren", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->tree->nchildren;
}

/** gets siblings of focus node along with the number of siblings
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetSiblings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          siblings,           /**< pointer to store siblings array, or NULL if not needed */
   int*                  nsiblings           /**< pointer to store number of siblings, or NULL if not needed */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetSiblings", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( siblings != NULL )
      *siblings = scip->tree->siblings;
   if( nsiblings != NULL )
      *nsiblings = scip->tree->nsiblings;

   return SCIP_OKAY;
}

/** gets number of siblings of focus node
 *
 *  @return the number of siblings of focus node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNSiblings(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNSiblings", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->tree->nsiblings;
}

/** gets leaves of the tree along with the number of leaves
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetLeaves(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          leaves,             /**< pointer to store leaves array, or NULL if not needed */
   int*                  nleaves             /**< pointer to store number of leaves, or NULL if not needed */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetLeaves", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( leaves != NULL )
      *leaves = SCIPnodepqNodes(scip->tree->leaves);
   if( nleaves != NULL )
      *nleaves = SCIPnodepqLen(scip->tree->leaves);

   return SCIP_OKAY;
}

/** gets number of leaves in the tree
 *
 *  @return the number of leaves in the tree
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLeaves", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPnodepqLen(scip->tree->leaves);
}

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @return the best child of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NODE* SCIPgetPrioChild(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPrioChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetPrioChild(scip->tree);
}

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @return the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NODE* SCIPgetPrioSibling(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPrioSibling", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetPrioSibling(scip->tree);
}

/** gets the best child of the focus node w.r.t. the node selection strategy
 *
 *  @return the best child of the focus node w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NODE* SCIPgetBestChild(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBestChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestChild(scip->tree, scip->set);
}

/** gets the best sibling of the focus node w.r.t. the node selection strategy
 *
 *  @return the best sibling of the focus node w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NODE* SCIPgetBestSibling(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBestSibling", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestSibling(scip->tree, scip->set);
}

/** gets the best leaf from the node queue w.r.t. the node selection strategy
 *
 *  @return the best leaf from the node queue w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NODE* SCIPgetBestLeaf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBestLeaf", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestLeaf(scip->tree);
}

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy
 *
 *  @return the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NODE* SCIPgetBestNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBestNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestNode(scip->tree, scip->set);
}

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf)
 *
 *  @return the node with smallest lower bound from the tree (child, sibling, or leaf)
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NODE* SCIPgetBestboundNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBestboundNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetLowerboundNode(scip->tree, scip->set);
}

/** access to all data of open nodes (leaves, children, and siblings)
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetOpenNodesData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          leaves,             /**< pointer to store the leaves, or NULL if not needed */
   SCIP_NODE***          children,           /**< pointer to store the children, or NULL if not needed */
   SCIP_NODE***          siblings,           /**< pointer to store the siblings, or NULL if not needed */
   int*                  nleaves,            /**< pointer to store the number of leaves, or NULL */
   int*                  nchildren,          /**< pointer to store the number of children, or NULL */
   int*                  nsiblings           /**< pointer to store the number of siblings, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetOpenNodesData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( leaves != NULL )
      *leaves = SCIPnodepqNodes(scip->tree->leaves);
   if( children != NULL )
      *children = scip->tree->children;
   if( siblings != NULL )
      *siblings = scip->tree->siblings;
   if( nleaves != NULL )
      *nleaves = SCIPnodepqLen(scip->tree->leaves);
   if( nchildren != NULL )
      *nchildren = SCIPtreeGetNChildren(scip->tree);
   if( nsiblings != NULL )
      *nsiblings = SCIPtreeGetNSiblings(scip->tree);

   return SCIP_OKAY;
}

/** cuts off node and whole sub tree from branch and bound tree
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcutoffNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be cut off */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcutoffNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnodeCutoff(node, scip->set, scip->stat, scip->tree, scip->transprob, scip->origprob, scip->reopt,
         scip->lp, scip->mem->probmem) );

   return SCIP_OKAY;
}

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
SCIP_RETCODE SCIPpruneTree(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPpruneTree", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPtreeCutoff(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
         scip->eventqueue, scip->lp, SCIPinfinity(scip)) );

   return SCIP_OKAY;
}

/** marks the given node to be propagated again the next time a node of its subtree is processed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPrepropagateNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be propagated again */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPrepropagateNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPnodePropagateAgain(node, scip->set, scip->stat, scip->tree);

   return SCIP_OKAY;
}

/** returns depth of first node in active path that is marked being cutoff
 *
 *  @return depth of first node in active path that is marked being cutoff
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetCutoffdepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetCutoffdepth", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return scip->tree->cutoffdepth;
}

/** returns depth of first node in active path that has to be propagated again
 *
 *  @return depth of first node in active path that has to be propagated again
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetRepropdepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRepropdepth", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return scip->tree->repropdepth;
}

/** prints all branching decisions on variables from the root to the given node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPprintNodeRootPath(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR**            branchvars;         /* array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real*            branchbounds;       /* array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE*       boundtypes;         /* array of boundtypes which the branchings in all ancestors set */
   int*                  nodeswitches;       /* marks, where in the arrays the branching decisions of the next node on the path start
                                              * branchings performed at the parent of node always start at position 0. For single variable branching,
                                              * nodeswitches[i] = i holds */
   int                   nbranchvars;        /* number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int                   branchvarssize;     /* available slots in arrays */
   int                   nnodes;             /* number of nodes in the nodeswitch array */
   int                   nodeswitchsize;     /* available slots in node switch array */

   branchvarssize = SCIPnodeGetDepth(node);
   nodeswitchsize = branchvarssize;

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &branchbounds, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeswitches, nodeswitchsize) );

   SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize );

   /* if the arrays were to small, we have to reallocate them and recall SCIPnodeGetAncestorBranchingPath */
   if( nbranchvars > branchvarssize || nnodes > nodeswitchsize )
   {
      branchvarssize = nbranchvars;
      nodeswitchsize = nnodes;

      /* memory reallocation */
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchvars, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchbounds, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &boundtypes, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &nodeswitches, nodeswitchsize) );

      SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize);
      assert(nbranchvars == branchvarssize);
   }

   /* we only want to create output, if branchings were performed */
   if( nbranchvars >= 1 )
   {
      int i;
      int j;

      /* print all nodes, starting from the root, which is last in the arrays */
      for( j = nnodes-1; j >= 0; --j)
      {
         int end;
         if(j == nnodes-1)
            end =  nbranchvars;
         else
            end =  nodeswitches[j+1];

         for( i = nodeswitches[j]; i < end; ++i )
         {
            if( i > nodeswitches[j] )
               SCIPmessageFPrintInfo(scip->messagehdlr, file, " AND ");
            SCIPmessageFPrintInfo(scip->messagehdlr, file, "<%s> %s %.1f",SCIPvarGetName(branchvars[i]), boundtypes[i] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", branchbounds[i]);
         }
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");
         if( j > 0 )
         {
            if(  nodeswitches[j]-nodeswitches[j-1] != 1 )
               SCIPmessageFPrintInfo(scip->messagehdlr, file, " |\n |\n");
            else if( boundtypes[i-1] == SCIP_BOUNDTYPE_LOWER )
               SCIPmessageFPrintInfo(scip->messagehdlr, file, "\\ \n \\\n");
            else
               SCIPmessageFPrintInfo(scip->messagehdlr, file, " /\n/ \n");
         }
      }
   }

   /* free all local memory */
   SCIPfreeBufferArray(scip, &nodeswitches);
   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &branchbounds);
   SCIPfreeBufferArray(scip, &branchvars);

   return SCIP_OKAY;
}

/** sets whether the LP should be solved at the focus node
 *
 *  @note In order to have an effect, this method needs to be called after a node is focused but before the LP is
 *        solved.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
void SCIPsetFocusnodeLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             solvelp             /**< should the LP be solved? */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPsetFocusnodeLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPtreeSetFocusNodeLP(scip->tree, solvelp);
}

/** gets number of nodes left in the tree (children + siblings + leaves)
 *
 *  @return the number of nodes left in the tree (children + siblings + leaves)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNNodesLeft(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNodesLeft", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetNNodes(scip->tree);
}

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
int SCIPgetDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetDepth", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetCurrentDepth(scip->tree);
}

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
int SCIPgetFocusDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetFocusDepth", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetFocusDepth(scip->tree);
}

/** gets current plunging depth (successive times, a child was selected as next node)
 *
 *  @return the current plunging depth (successive times, a child was selected as next node)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetPlungeDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPlungeDepth", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return scip->stat->plungedepth;
}

/** query if node was the last parent of a branching of the tree
 *
 *  @return TRUE if node was the last parent of a branching of the tree
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Bool SCIPwasNodeLastBranchParent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< tree node, or NULL to check focus node */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPwasNodeLastBranchParent", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeWasNodeLastBranchParent(scip->tree, node);
}
