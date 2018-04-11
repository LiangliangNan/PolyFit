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

/**@file   nodesel_uct.c
 * @brief  uct node selector which balances exploration and exploitation by considering node visits
 * @author Gregor Hendel
 *
 * the UCT node selection rule selects the next leaf according to a mixed score of the node's actual lower bound
 * and the number of times it has been visited so far compared to its parent node.
 *
 * The idea of UCT node selection for MIP appeared in:
 * Ashish Sabharwal and Horst Samulowitz
 * Guiding Combinatorial Optimization with UCT (2011)
 *
 * The authors adapted a game-tree exploration scheme called UCB to MIP trees. Starting from the root node as current node,
 * the algorithm selects the current node's child \f$N_i\f$ which maximizes the UCT score
 *
 * \f$ \mbox{score}(N_i) := -\mbox{estimate}_{N_i} + \mbox{weight} \cdot \frac{\mbox{visits}(\mbox{parent}(N_i))}{\mbox{visits}(N_i)}
 * \f$
 *
 * where \f$\mbox{estimate}\f$ is the node's lower bound normalized by the root lower bound, and \f$\mbox{visits}\f$
 * denotes the number of times a leaf in the subtree rooted at this node has been explored so far.
 *
 * The selected node in the sense of the SCIP node selection is the leaf reached by the above criterion.
 *
 * The authors suggest that this node selection rule is particularly useful at the beginning of the solving process, but
 * to switch to a different node selection after a number of nodes has been explored to reduce computational overhead.
 * Our implementation uses only information available from the original SCIP tree which does not support the
 * forward path mechanism needed for the most efficient node selection. Instead, the algorithm selects the next leaf
 * by looping over all leaves and comparing the best leaf found so far with the next one. Two leaves l_1, l_2 are compared
 * by following their paths back upwards until their deepest common ancestor \f$a\f$ is reached, together with the two
 * children of \f$a\f$ representing the two paths to l_1, l_2. The leaf represented by the child of \f$a\f$
 * with higher UCT score is a candidate for the next selected leaf.
 *
 * The node selector features several parameters:
 *
 * the nodelimit delimits the number of explored nodes before UCT selection is turned off
 * the weight parameter changes the relevance of the visits quotient in the UCT score (see above score formula)
 * useestimate determines whether the node's estimate or lower bound is taken as estimate
 *
 * @note It should be avoided to switch to uct node selection after the branch and bound process has begun because
 *       the central UCT score information how often a path was taken is not collected if UCT is inactive. A safe use of
 *       UCT is to switch it on before SCIP starts optimization.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/nodesel_uct.h"

#define NODESEL_NAME            "uct"
#define NODESEL_DESC            "node selector which balances exploration and exploitation "
#define NODESEL_STDPRIORITY     10
#define NODESEL_MEMSAVEPRIORITY 0

/** default values for user parameters */
#define DEFAULT_WEIGHT          0.1     /**< weight of node visits in UCT score */
#define DEFAULT_NODELIMIT       31      /**< limit of node selections after which UCT node selection is turned off */
#define DEFAULT_USEESTIMATE     FALSE   /**< should the estimate (TRUE) or the lower bound of a node be used for UCT score? */
#define INITIALSIZE             1024    /**< initial size of node visits array (increased dynamically if required) */
#define MAXNODELIMIT            1000000 /**< the maximum value for user parameter nodelimit */
/*
 * Data structures
 */

/** node selector data */
struct SCIP_NodeselData
{
   int*                  nodevisits;         /**< array to store the number of node visits so far for every node */
   SCIP_Real             weight;             /**< weight of node visits in UCT score */
   int                   nodelimit;          /**< limit of node selections after which UCT node selection is turned off */
   int                   sizenodevisits;     /**< the size of the visits array */
   int                   nselections;        /**< counter for the number of node selections */
   int                   origstdpriority;    /**< priority of node selector when starting branch and bound */
   SCIP_Bool             useestimate;        /**< should the estimate (TRUE) or the lower bound of a node be used for UCT score? */
};

/*
 * Local methods
 */

/** get the number times @p node has been visited so far */
static
int nodeGetVisits(
   SCIP_NODESELDATA*     nodeseldata,        /**< node selector data */
   SCIP_NODE*            node                /**< the node in question */
   )
{
   int nodenumber;

   assert(nodeseldata != NULL);
   assert(node != NULL);

   /* nodes numbers start with 1 for the root node */
   nodenumber = (int)(SCIPnodeGetNumber(node) - 1);
   assert(nodenumber >= 0);

   if( nodenumber >= nodeseldata->sizenodevisits )
      return 0;
   else
      return nodeseldata->nodevisits[nodenumber];
}

/** increases the visits counter along the path from @p node to the root node */
static
void updateVisits(
   SCIP_NODESELDATA*     nodeseldata,        /**< node selector data */
   SCIP_NODE*            node                /**< leaf node of path along which the visits are backpropagated */
   )
{
   int* visits;

   assert(nodeseldata != NULL);
   assert(node != NULL);

   visits = nodeseldata->nodevisits;
   assert(visits != NULL);

   /* increase visits counter of all nodes along the path until root node is reached (which has NULL as parent) */
   do
   {
      int nodenumber;

      nodenumber = (int)(SCIPnodeGetNumber(node) - 1);
      if( nodenumber < nodeseldata->sizenodevisits )
         ++(visits[nodenumber]);

      assert(SCIPnodeGetParent(node) == NULL || SCIPnodeGetDepth(node) >= 1);
      node = SCIPnodeGetParent(node);
   }
   while( node != NULL );
}

/** switches to a different node selection rule by assigning the lowest priority of all node selectors to uct */
static
SCIP_RETCODE turnoffNodeSelector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel             /**< the node selector to be turned off */
   )
{
   SCIP_NODESEL** nodesels;
   int nnodesels;
   int newpriority;
   int prio;
   int n;

   nodesels = SCIPgetNodesels(scip);
   nnodesels = SCIPgetNNodesels(scip);
   newpriority = SCIPnodeselGetStdPriority(nodesel);

   /* loop over node selectors to find minimum priority */
   for( n = 0; n < nnodesels; ++n )
   {
      prio = SCIPnodeselGetStdPriority(nodesels[n]);
      newpriority = MIN(newpriority, prio);
   }
   newpriority = MAX(newpriority, INT_MIN + 1);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Reached node limit of UCT node selection rule -> switching to default\n");
   SCIP_CALL( SCIPsetNodeselStdPriority(scip, nodesel, newpriority - 1) );

   return SCIP_OKAY;
}

/** returns UCT score of @p node; the UCT score is a mixture of the node's lower bound or estimate and the number of times
 *  it has been visited so far in relation with the number of times its parent has been visited so far
 */
static
SCIP_Real nodeGetUctScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< the node for which UCT score is requested */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   SCIP_NODE* parent;
   SCIP_Real rootlowerbound;
   SCIP_Real score;
   int parentvisits;

   rootlowerbound = SCIPgetLowerboundRoot(scip);

   /* the objective part of the UCT score uses the (negative) gap between node estimate and root lower bound */
   score = nodeseldata->useestimate ? SCIPnodeGetEstimate(node) : SCIPnodeGetLowerbound(node);

   /* if the root lower bound is infinite due to LP errors, we ignore the gap part of the UCT score */
   if( !SCIPisInfinity(scip, REALABS(rootlowerbound)) && !SCIPisEQ(scip, score, rootlowerbound) )
   {
      SCIP_Real absscore;
      SCIP_Real absrootlowerbound;
      SCIP_Real minabs;

      assert(SCIPisGE(scip, score, rootlowerbound));
      absscore = REALABS(score);
      absrootlowerbound = REALABS(rootlowerbound);
      minabs = MIN(absscore, absrootlowerbound);
      minabs = MAX(minabs, 1.0);

      score = (rootlowerbound - score) / minabs;
   }
   else
      score = 0.0;

   /* the visits part of the UCT score function */
   parent = SCIPnodeGetParent(node);
   assert(parent != NULL);
   parentvisits = nodeGetVisits(nodeseldata, parent);

   if( parentvisits > 0 )
   {
      int visits;

      visits = nodeGetVisits(nodeseldata, node);
      score += nodeseldata->weight * parentvisits / (SCIP_Real)(1 + visits);
   }

   return score;
}

/** compares two leaf nodes by comparing the UCT scores of the two children of their deepest common ancestor */
static
int compareNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESELDATA*     nodeseldata,        /**< node selector data */
   SCIP_NODE*            node1,              /**< first node for comparison */
   SCIP_NODE*            node2               /**< second node for comparisons */
   )
{
   SCIP_Real score1;
   SCIP_Real score2;

   assert(node1 != node2);

   /* go back in the tree to find the two shallowest ancestors of node1 and node2 which share the same parent */
   while( SCIPnodeGetParent(node1) != SCIPnodeGetParent(node2) )
   {
      /* if the nodes have the same depth but not the same parent both pointers can be updated, otherwise only the deeper
       * node pointer is moved
       */
      if( SCIPnodeGetDepth(node1) == SCIPnodeGetDepth(node2) )
      {
         node1 = SCIPnodeGetParent(node1);
         node2 = SCIPnodeGetParent(node2);
      }
      else if( SCIPnodeGetDepth(node1) > SCIPnodeGetDepth(node2) )
         node1 = SCIPnodeGetParent(node1);
      else if( SCIPnodeGetDepth(node1) < SCIPnodeGetDepth(node2) )
         node2 = SCIPnodeGetParent(node2);

      assert(node1 != NULL);
      assert(node2 != NULL);
   }

   /* get UCT scores for both nodes */
   score1 = nodeGetUctScore(scip, node1, nodeseldata);
   score2 = nodeGetUctScore(scip, node2, nodeseldata);

   if( (SCIPisInfinity(scip, score1) && SCIPisInfinity(scip, score2)) ||
      (SCIPisInfinity(scip, -score1) && SCIPisInfinity(scip, -score2)) ||
      SCIPisEQ(scip, score1, score2) )
   {
      return 0;
   }
   else if( SCIPisLT(scip, score1, score2) )
      return -1;
   else
   {
      assert(SCIPisGT(scip, score1, score2));
      return 1;
   }
}

/** selects the best node among @p nodes with respect to UCT score */
static
void selectBestNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE**           selnode,            /**< pointer to store the selected node, needs not be empty */
   SCIP_NODESELDATA*     nodeseldata,        /**< node selector data */
   SCIP_NODE**           nodes,              /**< array of nodes to select from */
   int                   nnodes              /**< size of the nodes array */
   )
{
   int n;

   assert(nnodes == 0 || nodes != NULL);
   assert(nnodes >= 0);
   assert(selnode != NULL);

   if( nnodes == 0 )
      return;

   /* loop over nodes, always keeping reference to the best found node so far */
   for( n = 0; n < nnodes; ++n )
   {
      assert(nodes[n] != NULL);
      /* update the selected node if the current node has a higher score */
      if( *selnode == NULL || compareNodes(scip, nodeseldata, *selnode, nodes[n]) < 0 )
         *selnode = nodes[n];
   }
}

/** keeps visits array large enough to save visits for all nodes in the branch and bound tree */
static
SCIP_RETCODE ensureMemorySize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   assert(nodeseldata != NULL);

   /* if array has not been allocated yet, do this now with default initial capacity */
   if( nodeseldata->nodevisits == NULL )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &nodeseldata->nodevisits, INITIALSIZE) ); /*lint !e506*/
      nodeseldata->sizenodevisits = INITIALSIZE;
   }

   assert(nodeseldata->nodelimit >= SCIPgetNNodes(scip));

   /* if user node limit has not been reached yet, resize the visits array if necessary */
   if( nodeseldata->sizenodevisits < 2 * nodeseldata->nodelimit && nodeseldata->sizenodevisits < (int)(2 * SCIPgetNNodes(scip)))
   {
      int newcapacity;
      newcapacity = MIN(2 * nodeseldata->sizenodevisits, 2 * nodeseldata->nodelimit);

      SCIPdebugMsg(scip, "Resizing node visits array, old capacity: %d new capacity : %d\n", nodeseldata->sizenodevisits, newcapacity);
      assert(newcapacity > nodeseldata->sizenodevisits);

      SCIP_CALL( SCIPreallocMemoryArray(scip, &nodeseldata->nodevisits, newcapacity) );
      BMSclearMemoryArray(&(nodeseldata->nodevisits[nodeseldata->sizenodevisits]), newcapacity - nodeseldata->sizenodevisits); /*lint !e866*/

      nodeseldata->sizenodevisits = newcapacity;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of node selector
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyUct)
{  /*lint --e{715}*/
   assert(scip != NULL);
   SCIP_CALL( SCIPincludeNodeselUct(scip) );

   return SCIP_OKAY;
}

/** solving process initialization method of node selector (called when branch and bound process is about to begin) */
static
SCIP_DECL_NODESELINITSOL(nodeselInitsolUct)
{
   SCIP_NODESELDATA* nodeseldata;
   assert(scip != NULL);
   assert(nodesel != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);

   assert(nodeseldata != NULL);
   nodeseldata->nselections = 0;
   nodeseldata->sizenodevisits = 0;
   nodeseldata->origstdpriority = SCIPnodeselGetStdPriority(nodesel);

   return SCIP_OKAY;
}

/** solving process deinitialization method of node selector (called when branch and bound process data gets freed) */
static
SCIP_DECL_NODESELEXITSOL(nodeselExitsolUct)
{
   SCIP_NODESELDATA* nodeseldata;
   assert(scip != NULL);
   assert(nodesel != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);

   assert(nodeseldata != NULL);

   if( nodeseldata->sizenodevisits > 0 )
   {
      assert(nodeseldata->nodevisits != NULL);
      SCIPfreeMemoryArray(scip, &nodeseldata->nodevisits);
   }
   nodeseldata->sizenodevisits = 0;
   nodeseldata->nselections = 0;

   /* reset node selector priority to its original value (before turning it off) */
   SCIP_CALL( SCIPsetNodeselStdPriority(scip, nodesel, nodeseldata->origstdpriority) );

   return SCIP_OKAY;
}

/** destructor of node selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_NODESELFREE(nodeselFreeUct)
{
   SCIP_NODESELDATA* nodeseldata;
   assert(scip != NULL);
   assert(nodesel != NULL);

   nodeseldata = SCIPnodeselGetData(nodesel);
   if( nodeseldata->sizenodevisits > 0 )
   {
      assert(nodeseldata->nodevisits != NULL);
      SCIPfreeMemoryArray(scip, &nodeseldata->nodevisits);
   }
   SCIPfreeBlockMemory(scip, &nodeseldata);

   SCIPnodeselSetData(nodesel, NULL);

   return SCIP_OKAY;
}

/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectUct)
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODE** leaves;
   SCIP_NODE** children;
   SCIP_NODE** siblings;
   int nleaves;
   int nsiblings;
   int nchildren;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   *selnode = NULL;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   if( nodeseldata->nodelimit < SCIPgetNNodes(scip) )
   {
      SCIPerrorMessage("UCT node limit exceeded\n");
      return SCIP_INVALIDCALL;
   }

   /* collect leaves, children and siblings data */
   SCIP_CALL( SCIPgetOpenNodesData(scip, &leaves, &children, &siblings, &nleaves, &nchildren, &nsiblings) );
   assert(nleaves + nchildren + nsiblings == SCIPgetNNodesLeft(scip));

   if( SCIPgetNNodesLeft(scip) == 0 )
      return SCIP_OKAY;

   /* make sure that UCT node selection data is large enough to store node visits */
   SCIP_CALL( ensureMemorySize(scip, nodeseldata) );

   /* select next node as best node with respect to UCT-based comparison method */
   selectBestNode(scip, selnode, nodeseldata, children, nchildren);
   selectBestNode(scip, selnode, nodeseldata, siblings, nsiblings);
   selectBestNode(scip, selnode, nodeseldata, leaves, nleaves);

   if( *selnode == NULL )
   {
      SCIPerrorMessage("Node selection rule UCT could not select a node.\n");
      return SCIP_INVALIDCALL;
   }

   /* increase the number of selections */
   ++nodeseldata->nselections;

   /* switch back to default node selection rule if the node limit is exceeded */
   if( nodeseldata->nselections == nodeseldata->nodelimit )
   {
      SCIP_CALL( turnoffNodeSelector(scip, nodesel) );
   }
   else
   {
      /* trigger update of visits along the path from the selected node to the root node */
      SCIPdebugMsg(scip, "updating node visits from node number %" SCIP_LONGINT_FORMAT "\n", SCIPnodeGetNumber(*selnode));
      updateVisits(nodeseldata, *selnode);
   }

   return SCIP_OKAY;
}

/** node comparison method of UCT node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompUct)
{  /*lint --e{715}*/
   SCIP_Real lowerbound1;
   SCIP_Real lowerbound2;

   lowerbound1 = SCIPnodeGetLowerbound(node1);
   lowerbound2 = SCIPnodeGetLowerbound(node2);

   if( SCIPisLT(scip, lowerbound1, lowerbound2) )
      return -1;
   else if( SCIPisGT(scip, lowerbound1, lowerbound2) )
      return 1;
   else
      return 0;
}

/*
 * node selector specific interface methods
 */

/** creates the uct node selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselUct(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODESEL* nodesel;

   /* create uct node selector data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nodeseldata) );

   nodesel = NULL;
   nodeseldata->nodevisits = NULL;
   nodeseldata->nselections = 0;
   nodeseldata->sizenodevisits = 0;
   nodeseldata->origstdpriority = NODESEL_STDPRIORITY;

   /* use SCIPincludeNodeselBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY,
          NODESEL_MEMSAVEPRIORITY, nodeselSelectUct, nodeselCompUct, nodeseldata) );

   assert(nodesel != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyUct) );
   SCIP_CALL( SCIPsetNodeselInitsol(scip, nodesel, nodeselInitsolUct) );
   SCIP_CALL( SCIPsetNodeselFree(scip, nodesel, nodeselFreeUct) );
   SCIP_CALL( SCIPsetNodeselExitsol(scip, nodesel, nodeselExitsolUct) );

   /* add uct node selector parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "nodeselection/" NODESEL_NAME "/nodelimit",
         "maximum number of nodes before switching to default rule",
         &nodeseldata->nodelimit, TRUE, DEFAULT_NODELIMIT, 0, MAXNODELIMIT, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "nodeselection/" NODESEL_NAME "/weight",
         "weight for visit quotient of node selection rule",
         &nodeseldata->weight, TRUE, DEFAULT_WEIGHT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "nodeselection/" NODESEL_NAME "/useestimate",
         "should the estimate (TRUE) or lower bound of a node be used for UCT score?",
         &nodeseldata->useestimate, TRUE, DEFAULT_USEESTIMATE, NULL, NULL) );

   return SCIP_OKAY;
}
