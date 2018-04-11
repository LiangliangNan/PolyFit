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

/**@file   nodesel.c
 * @brief  methods for node selectors
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/stat.h"
#include "scip/visual.h"
#include "scip/paramset.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/lp.h"
#include "scip/scip.h"
#include "scip/nodesel.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_nodesel.h"
#include "scip/struct_scip.h"

/* 
 * node priority queue methods
 */

#define PQ_PARENT(q) (((q)+1)/2-1)
#define PQ_LEFTCHILD(p) (2*(p)+1)
#define PQ_RIGHTCHILD(p) (2*(p)+2)


/** node comparator for node numbers */
static
SCIP_DECL_SORTPTRCOMP(nodeCompNumber)
{  /*lint --e{715}*/
   assert(elem1 != NULL);
   assert(elem2 != NULL);

   if( ((SCIP_NODE*)elem1)->number < ((SCIP_NODE*)elem2)->number )
      return -1;
   else if( ((SCIP_NODE*)elem1)->number > ((SCIP_NODE*)elem2)->number )
      return +1;
   else
   {
      /* no such two nodes should have the same node number */
      assert(elem1 == elem2);
      return 0;
   }
}


/** resizes node memory to hold at least the given number of nodes */
static
SCIP_RETCODE nodepqResize(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   minsize             /**< minimal number of storable nodes */
   )
{
   assert(nodepq != NULL);

   if( minsize <= nodepq->size )
      return SCIP_OKAY;

   nodepq->size = SCIPsetCalcTreeGrowSize(set, minsize);
   SCIP_ALLOC( BMSreallocMemoryArray(&nodepq->slots, nodepq->size) );
   SCIP_ALLOC( BMSreallocMemoryArray(&nodepq->bfsposs, nodepq->size) );
   SCIP_ALLOC( BMSreallocMemoryArray(&nodepq->bfsqueue, nodepq->size) );

   return SCIP_OKAY;
}

/** creates node priority queue */
SCIP_RETCODE SCIPnodepqCreate(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   )
{  /*lint --e{715}*/
   assert(nodepq != NULL);

   SCIP_ALLOC( BMSallocMemory(nodepq) );
   (*nodepq)->nodesel = nodesel;
   (*nodepq)->slots = NULL;
   (*nodepq)->bfsposs = NULL;
   (*nodepq)->bfsqueue = NULL;
   (*nodepq)->len = 0;
   (*nodepq)->size = 0;
   (*nodepq)->lowerboundsum = 0.0;

   return SCIP_OKAY;
}

/** frees node priority queue, but not the data nodes themselves */
void SCIPnodepqDestroy(
   SCIP_NODEPQ**         nodepq              /**< pointer to a node priority queue */
   )
{
   assert(nodepq != NULL);
   assert(*nodepq != NULL);

   BMSfreeMemoryArrayNull(&(*nodepq)->slots);
   BMSfreeMemoryArrayNull(&(*nodepq)->bfsposs);
   BMSfreeMemoryArrayNull(&(*nodepq)->bfsqueue);
   BMSfreeMemory(nodepq);
}

/** frees node priority queue and all nodes in the queue */
SCIP_RETCODE SCIPnodepqFree(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(nodepq != NULL);
   assert(*nodepq != NULL);

   /* free the nodes of the queue */
   SCIP_CALL( SCIPnodepqClear(*nodepq, blkmem, set, stat, eventqueue, tree, lp) );

   /* free the queue data structure */
   SCIPnodepqDestroy(nodepq);

   return SCIP_OKAY;
}

/** deletes all nodes in the node priority queue */
SCIP_RETCODE SCIPnodepqClear(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(nodepq != NULL);

   if( nodepq->len > 0 )
   {
      /* sort the sorts downwards after their number to increase speed when freeing in debug mode */
      /* @todo: if a node is freed, the parent will also be freed, if no children are left; maybe we want to free all
       *        nodes in the order of decreasing node numbers
       */
      SCIPsortDownPtr((void**)nodepq->slots, nodeCompNumber, nodepq->len);

      /* free the nodes of the queue */
      for( i = 0; i < nodepq->len; ++i )
      {
         assert(nodepq->slots[i] != NULL);
         assert(SCIPnodeGetType(nodepq->slots[i]) == SCIP_NODETYPE_LEAF);

         SCIP_CALL( SCIPnodeFree(&nodepq->slots[i], blkmem, set, stat, eventqueue, tree, lp) );
      }
   }

   /* reset data */
   nodepq->len = 0;
   nodepq->lowerboundsum = 0.0;

   return SCIP_OKAY;
}

/** returns the node selector associated with the given node priority queue */
SCIP_NODESEL* SCIPnodepqGetNodesel(
   SCIP_NODEPQ*          nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);

   return nodepq->nodesel;
}

/** sets the node selector used for sorting the nodes in the queue, and resorts the queue if necessary */
SCIP_RETCODE SCIPnodepqSetNodesel(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   )
{
   SCIP_NODEPQ* newnodepq;
   SCIP_RETCODE retcode;
   int i;

   assert(nodepq != NULL);
   assert(*nodepq != NULL);
   assert((*nodepq)->len >= 0);
   assert(nodesel != NULL);
   assert(nodesel->nodeselcomp != NULL);

   if( (*nodepq)->nodesel == nodesel )
      return SCIP_OKAY;

   /* create new node priority queue */
   SCIP_CALL( SCIPnodepqCreate(&newnodepq, set, nodesel) );

   /* resize the new node priority queue to be able to store all nodes */
   retcode = nodepqResize(newnodepq, set, (*nodepq)->len);

   /* insert all nodes in the new node priority queue */
   for( i = 0; i < (*nodepq)->len && retcode == SCIP_OKAY; ++i )
   {
      retcode = SCIPnodepqInsert(newnodepq, set, (*nodepq)->slots[i]);
   }

   if( retcode != SCIP_OKAY )
   {
      SCIPnodepqDestroy(&newnodepq);

      return retcode;
   }

   /* destroy the old node priority queue without freeing the nodes */
   SCIPnodepqDestroy(nodepq);

   /* use the new node priority queue */
   *nodepq = newnodepq;

   return SCIP_OKAY;
}

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
int SCIPnodepqCompare(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node1,              /**< first node to compare */
   SCIP_NODE*            node2               /**< second node to compare */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->nodesel != NULL);
   assert(nodepq->nodesel->nodeselcomp != NULL);
   assert(set != NULL);

   return SCIPnodeselCompare(nodepq->nodesel, set, node1, node2);
}

/** inserts node into node priority queue */
SCIP_RETCODE SCIPnodepqInsert(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to be inserted */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE** slots;
   int* bfsposs;
   int* bfsqueue;
   SCIP_Real lowerbound;
   int pos;
   int bfspos;

   assert(nodepq != NULL);
   assert(nodepq->len >= 0);
   assert(set != NULL);
   assert(node != NULL);

   nodesel = nodepq->nodesel;
   assert(nodesel != NULL);
   assert(nodesel->nodeselcomp != NULL);

   SCIP_CALL( nodepqResize(nodepq, set, nodepq->len+1) );
   slots = nodepq->slots;
   bfsposs = nodepq->bfsposs;
   bfsqueue = nodepq->bfsqueue;

   /* insert node as leaf in the tree, move it towards the root as long it is better than its parent */
   nodepq->len++;
   nodepq->lowerboundsum += SCIPnodeGetLowerbound(node);
   pos = nodepq->len-1;
   while( pos > 0 && nodesel->nodeselcomp(set->scip, nodesel, node, slots[PQ_PARENT(pos)]) < 0 )
   {
      slots[pos] = slots[PQ_PARENT(pos)];
      bfsposs[pos] = bfsposs[PQ_PARENT(pos)];
      bfsqueue[bfsposs[pos]] = pos;
      pos = PQ_PARENT(pos);
   }
   slots[pos] = node;

   /* insert the final position into the bfs index queue */
   lowerbound = SCIPnodeGetLowerbound(node);
   bfspos = nodepq->len-1;
   while( bfspos > 0 && lowerbound < SCIPnodeGetLowerbound(slots[bfsqueue[PQ_PARENT(bfspos)]]) )
   {
      bfsqueue[bfspos] = bfsqueue[PQ_PARENT(bfspos)];
      bfsposs[bfsqueue[bfspos]] = bfspos;
      bfspos = PQ_PARENT(bfspos);
   }
   bfsqueue[bfspos] = pos;
   bfsposs[pos] = bfspos;

   SCIPsetDebugMsg(set, "inserted node %p[%g] at pos %d and bfspos %d of node queue\n", (void*)node, lowerbound, pos, bfspos);

   return SCIP_OKAY;
}

/** deletes node at given position from the node priority queue; returns TRUE, if the parent fell down to the
 *  free position
 */
static
SCIP_Bool nodepqDelPos(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   rempos              /**< queue position of node to remove */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE** slots;
   int* bfsposs;
   int* bfsqueue;
   SCIP_NODE* lastnode;
   int lastbfspos;
   int lastbfsqueueidx;
   int freepos;
   int freebfspos;
   SCIP_Bool parentfelldown;
   SCIP_Bool bfsparentfelldown;

   assert(nodepq != NULL);
   assert(nodepq->len > 0);
   assert(set != NULL);
   assert(0 <= rempos && rempos < nodepq->len);

   nodesel = nodepq->nodesel;
   assert(nodesel != NULL);
   assert(nodesel->nodeselcomp != NULL);

   slots = nodepq->slots;
   bfsposs = nodepq->bfsposs;
   bfsqueue = nodepq->bfsqueue;

   nodepq->lowerboundsum -= SCIPnodeGetLowerbound(slots[rempos]);
   freepos = rempos;
   freebfspos = bfsposs[rempos];
   assert(0 <= freebfspos && freebfspos < nodepq->len);

   SCIPsetDebugMsg(set, "delete node %p[%g] at pos %d and bfspos %d of node queue\n",
      (void*)slots[freepos], SCIPnodeGetLowerbound(slots[freepos]), freepos, freebfspos);

   /* remove node of the tree and get a free slot,
    * if the removed node was the last node of the queue
    *  - do nothing
    * if the last node of the queue is better than the parent of the removed node:
    *  - move the parent to the free slot, until the last node can be placed in the free slot
    * if the last node of the queue is not better than the parent of the free slot:
    *  - move the better child to the free slot until the last node can be placed in the free slot
    */
   nodepq->len--;

   /* process the slots queue ordered by the node selection comparator */
   lastnode = slots[nodepq->len];
   lastbfspos = bfsposs[nodepq->len];
   parentfelldown = FALSE;
   if( freepos < nodepq->len )
   {
      int parentpos;

      /* try to move parents downwards to insert last node */
      parentpos = PQ_PARENT(freepos);
      while( freepos > 0 && nodesel->nodeselcomp(set->scip, nodesel, lastnode, slots[parentpos]) < 0 )
      {
         slots[freepos] = slots[parentpos];
         bfsposs[freepos] = bfsposs[parentpos];
         bfsqueue[bfsposs[freepos]] = freepos;
         freepos = parentpos;
         parentpos = PQ_PARENT(freepos);
         parentfelldown = TRUE;
      }
      if( !parentfelldown )
      {
         /* downward moving of parents was not successful -> move children upwards */
         while( freepos <= PQ_PARENT(nodepq->len-1) ) /* as long as free slot has children... */
         {
            int childpos;
            int brotherpos;

            /* select the better child of free slot */
            childpos = PQ_LEFTCHILD(freepos);
            assert(childpos < nodepq->len);
            brotherpos = PQ_RIGHTCHILD(freepos);
            if( brotherpos < nodepq->len
               && nodesel->nodeselcomp(set->scip, nodesel, slots[brotherpos], slots[childpos]) < 0 )
               childpos = brotherpos;

            /* exit search loop if better child is not better than last node */
            if( nodesel->nodeselcomp(set->scip, nodesel, lastnode, slots[childpos]) <= 0 )
               break;

            /* move better child upwards, free slot is now the better child's slot */
            slots[freepos] = slots[childpos];
            bfsposs[freepos] = bfsposs[childpos];
            bfsqueue[bfsposs[freepos]] = freepos;
            freepos = childpos;
         }
      }
      assert(0 <= freepos && freepos < nodepq->len);
      assert(!parentfelldown || PQ_LEFTCHILD(freepos) < nodepq->len);
      slots[freepos] = lastnode;
      bfsposs[freepos] = lastbfspos;
      bfsqueue[lastbfspos] = freepos;
   }

   /* process the bfs queue ordered by the lower bound */
   lastbfsqueueidx = bfsqueue[nodepq->len];
   bfsparentfelldown = FALSE;
   if( freebfspos < nodepq->len )
   {
      SCIP_Real lastlowerbound;
      int parentpos;

      /* try to move parents downwards to insert last queue index */
      lastlowerbound = SCIPnodeGetLowerbound(slots[lastbfsqueueidx]);
      parentpos = PQ_PARENT(freebfspos);
      while( freebfspos > 0 && lastlowerbound < SCIPnodeGetLowerbound(slots[bfsqueue[parentpos]]) )
      {
         bfsqueue[freebfspos] = bfsqueue[parentpos];
         bfsposs[bfsqueue[freebfspos]] = freebfspos;
         freebfspos = parentpos;
         parentpos = PQ_PARENT(freebfspos);
         bfsparentfelldown = TRUE;
      }
      if( !bfsparentfelldown )
      {
         /* downward moving of parents was not successful -> move children upwards */
         while( freebfspos <= PQ_PARENT(nodepq->len-1) ) /* as long as free slot has children... */
         {
            int childpos;
            int brotherpos;

            /* select the better child of free slot */
            childpos = PQ_LEFTCHILD(freebfspos);
            assert(childpos < nodepq->len);
            brotherpos = PQ_RIGHTCHILD(freebfspos);
            if( brotherpos < nodepq->len
               && SCIPnodeGetLowerbound(slots[bfsqueue[brotherpos]]) < SCIPnodeGetLowerbound(slots[bfsqueue[childpos]]) )
               childpos = brotherpos;

            /* exit search loop if better child is not better than last node */
            if( lastlowerbound <= SCIPnodeGetLowerbound(slots[bfsqueue[childpos]]) )
               break;

            /* move better child upwards, free slot is now the better child's slot */
            bfsqueue[freebfspos] = bfsqueue[childpos];
            bfsposs[bfsqueue[freebfspos]] = freebfspos;
            freebfspos = childpos;
         }
      }
      assert(0 <= freebfspos && freebfspos < nodepq->len);
      assert(!bfsparentfelldown || PQ_LEFTCHILD(freebfspos) < nodepq->len);
      bfsqueue[freebfspos] = lastbfsqueueidx;
      bfsposs[lastbfsqueueidx] = freebfspos;
   }

   return parentfelldown;
}

/** returns the position of given node in the priority queue, or -1 if not existing */
static
int nodepqFindNode(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to find */
   )
{
   int pos;

   assert(nodepq != NULL);
   assert(nodepq->len >= 0);
   assert(set != NULL);
   assert(node != NULL);

   /* search the node in the queue */
   for( pos = 0; pos < nodepq->len && node != nodepq->slots[pos]; ++pos )
   {}

   if( pos == nodepq->len )
      pos = -1;

   return pos;
}

/** removes node from the node priority queue */
SCIP_RETCODE SCIPnodepqRemove(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to remove */
   )
{
   int pos;

   pos = nodepqFindNode(nodepq, set, node);
   if( pos == -1 )
   {
      SCIPerrorMessage("node doesn't exist in node priority queue\n");
      return SCIP_INVALIDDATA;
   }

   (void)nodepqDelPos(nodepq, set, pos);

   return SCIP_OKAY;
}

/** returns the best node of the queue without removing it */
SCIP_NODE* SCIPnodepqFirst(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->len >= 0);

   if( nodepq->len == 0 )
      return NULL;

   assert(nodepq->slots[0] != NULL);

   return nodepq->slots[0];
}

/** returns the nodes array of the queue */
SCIP_NODE** SCIPnodepqNodes(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);

   return nodepq->slots;
}

/** returns the number of nodes stored in the node priority queue */
int SCIPnodepqLen(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->len >= 0);

   return nodepq->len;
}

/** gets the minimal lower bound of all nodes in the queue */
SCIP_Real SCIPnodepqGetLowerbound(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->nodesel != NULL);
   assert(set != NULL);

   if( nodepq->len > 0 )
   {
      int bfspos;

      bfspos = nodepq->bfsqueue[0];
      assert(0 <= bfspos && bfspos < nodepq->len);
      assert(nodepq->slots[bfspos] != NULL);
      return SCIPnodeGetLowerbound(nodepq->slots[bfspos]);
   }
   else
      return SCIPsetInfinity(set);
}

/** gets the node with minimal lower bound of all nodes in the queue */
SCIP_NODE* SCIPnodepqGetLowerboundNode(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodepq != NULL);
   assert(nodepq->nodesel != NULL);
   assert(set != NULL);

   /* the node selector's compare method sorts the minimal lower bound to the front */
   if( nodepq->len > 0 )
   {
      int bfspos;

      bfspos = nodepq->bfsqueue[0];
      assert(0 <= bfspos && bfspos < nodepq->len);
      assert(nodepq->slots[bfspos] != NULL);
      return nodepq->slots[bfspos];
   }
   else
      return NULL;
}

/** gets the sum of lower bounds of all nodes in the queue */
SCIP_Real SCIPnodepqGetLowerboundSum(
   SCIP_NODEPQ*          nodepq              /**< node priority queue */
   )
{
   assert(nodepq != NULL);

   return nodepq->lowerboundsum;
}

/** free all nodes from the queue that are cut off by the given upper bound */
SCIP_RETCODE SCIPnodepqBound(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   )
{
   SCIP_NODE* node;
   int pos;
   SCIP_Bool parentfelldown;

   assert(nodepq != NULL);

   SCIPsetDebugMsg(set, "bounding node queue of length %d with cutoffbound=%g\n", nodepq->len, cutoffbound);
   pos = nodepq->len-1;
   while( pos >= 0 )
   {
      assert(pos < nodepq->len);
      node = nodepq->slots[pos];
      assert(node != NULL);
      assert(SCIPnodeGetType(node) == SCIP_NODETYPE_LEAF);
      if( SCIPsetIsGE(set, SCIPnodeGetLowerbound(node), cutoffbound) )
      {
         SCIPsetDebugMsg(set, "free node in slot %d (len=%d) at depth %d with lowerbound=%g\n",
            pos, nodepq->len, SCIPnodeGetDepth(node), SCIPnodeGetLowerbound(node));

         /* cut off node; because we looped from back to front, the existing children of the node must have a smaller
          * lower bound than the cut off value
          */
         assert(PQ_LEFTCHILD(pos) >= nodepq->len
            || SCIPsetIsLT(set, SCIPnodeGetLowerbound(nodepq->slots[PQ_LEFTCHILD(pos)]), cutoffbound));
         assert(PQ_RIGHTCHILD(pos) >= nodepq->len
            || SCIPsetIsLT(set, SCIPnodeGetLowerbound(nodepq->slots[PQ_RIGHTCHILD(pos)]), cutoffbound));

         /* free the slot in the node PQ */
         parentfelldown = nodepqDelPos(nodepq, set, pos);

         /* - if the slot was occupied by the parent, we have to check this slot (the parent) again; unfortunately,
          *   we will check the node which occupied the parent's slot again, even though it cannot be cut off;
          * - otherwise, the slot was the last slot or it was occupied by a node with a position greater than
          *   the current position; this node was already checked and we can decrease the position
          */
         if( !parentfelldown )
            pos--;

         SCIPvisualCutoffNode(stat->visual, set, stat, node, FALSE);

         if( set->reopt_enable )
         {
            assert(reopt != NULL);
            SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, node, SCIP_EVENTTYPE_NODEINFEASIBLE, lp,
                  SCIPlpGetSolstat(lp), SCIPnodeGetDepth(node) == 0, SCIPtreeGetFocusNode(tree) == node,
                  SCIPnodeGetLowerbound(node), SCIPtreeGetEffectiveRootDepth(tree)));
         }

         /* free memory of the node */
         SCIP_CALL( SCIPnodeFree(&node, blkmem, set, stat, eventqueue, tree, lp) );
      }
      else
         pos--;
   }
   SCIPsetDebugMsg(set, " -> bounded node queue has length %d\n", nodepq->len);

   return SCIP_OKAY;
}




/*
 * node selector methods 
 */

/** method to call, when the standard mode priority of a node selector was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNodeselStdPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetNodeselStdPriority() to mark the nodesels unsorted */
   SCIP_CALL( SCIPsetNodeselStdPriority(scip, (SCIP_NODESEL*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** method to call, when the memory saving mode priority of a node selector was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNodeselMemsavePriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetNodeselMemsavePriority() to mark the nodesels unsorted */
   SCIP_CALL( SCIPsetNodeselMemsavePriority(scip, (SCIP_NODESEL*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given node selector to a new scip */
SCIP_RETCODE SCIPnodeselCopyInclude(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( nodesel->nodeselcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including node selector %s in subscip %p\n", SCIPnodeselGetName(nodesel), (void*)set->scip);
      SCIP_CALL( nodesel->nodeselcopy(set->scip, nodesel) );
   }
   return SCIP_OKAY;
}

/** creates a node selector */
SCIP_RETCODE SCIPnodeselCreate(
   SCIP_NODESEL**        nodesel,            /**< pointer to store node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy)),   /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   SCIP_DECL_NODESELINITSOL((*nodeselinitsol)),/**< solving process initialization method of node selector */
   SCIP_DECL_NODESELEXITSOL((*nodeselexitsol)),/**< solving process deinitialization method of node selector */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(nodesel != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(nodeselselect != NULL);
   assert(nodeselcomp != NULL);

   SCIP_ALLOC( BMSallocMemory(nodesel) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*nodesel)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*nodesel)->desc, desc, strlen(desc)+1) );
   (*nodesel)->stdpriority = stdpriority;
   (*nodesel)->memsavepriority = memsavepriority;
   (*nodesel)->nodeselcopy = nodeselcopy;
   (*nodesel)->nodeselfree = nodeselfree;
   (*nodesel)->nodeselinit = nodeselinit;
   (*nodesel)->nodeselexit = nodeselexit;
   (*nodesel)->nodeselinitsol = nodeselinitsol;
   (*nodesel)->nodeselexitsol = nodeselexitsol;
   (*nodesel)->nodeselselect = nodeselselect;
   (*nodesel)->nodeselcomp = nodeselcomp;
   (*nodesel)->nodeseldata = nodeseldata;
   (*nodesel)->initialized = FALSE;
   /* create clocks */
   SCIP_CALL( SCIPclockCreate(&(*nodesel)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*nodesel)->nodeseltime, SCIP_CLOCKTYPE_DEFAULT) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "nodeselection/%s/stdpriority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of node selection rule <%s> in standard mode", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*nodesel)->stdpriority, FALSE, stdpriority, INT_MIN/4, INT_MAX/2,
                  paramChgdNodeselStdPriority, (SCIP_PARAMDATA*)(*nodesel)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "nodeselection/%s/memsavepriority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of node selection rule <%s> in memory saving mode", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*nodesel)->memsavepriority, TRUE, memsavepriority, INT_MIN/4, INT_MAX/4,
                  paramChgdNodeselMemsavePriority, (SCIP_PARAMDATA*)(*nodesel)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** frees memory of node selector */
SCIP_RETCODE SCIPnodeselFree(
   SCIP_NODESEL**        nodesel,            /**< pointer to node selector data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(*nodesel != NULL);
   assert(!(*nodesel)->initialized);
   assert(set != NULL);

   /* call destructor of node selector */
   if( (*nodesel)->nodeselfree != NULL )
   {
      SCIP_CALL( (*nodesel)->nodeselfree(set->scip, *nodesel) );
   }

   /* free clocks */
   SCIPclockFree(&(*nodesel)->nodeseltime);
   SCIPclockFree(&(*nodesel)->setuptime);

   BMSfreeMemoryArray(&(*nodesel)->name);
   BMSfreeMemoryArray(&(*nodesel)->desc);
   BMSfreeMemory(nodesel);

   return SCIP_OKAY;
}

/** initializes node selector */
SCIP_RETCODE SCIPnodeselInit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   if( nodesel->initialized )
   {
      SCIPerrorMessage("node selector <%s> already initialized", nodesel->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(nodesel->setuptime);
      SCIPclockReset(nodesel->nodeseltime);
   }

   if( nodesel->nodeselinit != NULL )
   {
      /* start timing */
      SCIPclockStart(nodesel->setuptime, set);

      SCIP_CALL( nodesel->nodeselinit(set->scip, nodesel) );

      /* stop timing */
      SCIPclockStop(nodesel->setuptime, set);
   }
   nodesel->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes node selector */
SCIP_RETCODE SCIPnodeselExit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   if( !nodesel->initialized )
   {
      SCIPerrorMessage("node selector <%s> not initialized", nodesel->name);
      return SCIP_INVALIDCALL;
   }

   if( nodesel->nodeselexit != NULL )
   {
      /* start timing */
      SCIPclockStart(nodesel->setuptime, set);

      SCIP_CALL( nodesel->nodeselexit(set->scip, nodesel) );

      /* stop timing */
      SCIPclockStop(nodesel->setuptime, set);
   }
   nodesel->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs node selector that the branch and bound process is being started */
SCIP_RETCODE SCIPnodeselInitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   /* call solving process initialization method of node selector */
   if( nodesel->nodeselinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(nodesel->setuptime, set);

      SCIP_CALL( nodesel->nodeselinitsol(set->scip, nodesel) );

      /* stop timing */
      SCIPclockStop(nodesel->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs node selector that the branch and bound process data is being freed */
SCIP_RETCODE SCIPnodeselExitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of node selector */
   if( nodesel->nodeselexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(nodesel->setuptime, set);

      SCIP_CALL( nodesel->nodeselexitsol(set->scip, nodesel) );

      /* stop timing */
      SCIPclockStop(nodesel->setuptime, set);
   }

   return SCIP_OKAY;
}

/** select next node to be processed */
SCIP_RETCODE SCIPnodeselSelect(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE**           selnode             /**< pointer to store node to be processed next */
   )
{

   assert(nodesel != NULL);
   assert(nodesel->nodeselselect != NULL);
   assert(set != NULL);
   assert(selnode != NULL);

   /* start timing */
   SCIPclockStart(nodesel->nodeseltime, set);

   SCIP_CALL( nodesel->nodeselselect(set->scip, nodesel, selnode) );

   /* stop timing */
   SCIPclockStop(nodesel->nodeseltime, set);

   return SCIP_OKAY;
}

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
int SCIPnodeselCompare(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node1,              /**< first node to compare */
   SCIP_NODE*            node2               /**< second node to compare */
   )
{
   assert(nodesel != NULL);
   assert(nodesel->nodeselcomp != NULL);
   assert(set != NULL);
   assert(node1 != NULL);
   assert(node2 != NULL);

   return nodesel->nodeselcomp(set->scip, nodesel, node1, node2);
}

/** gets name of node selector */
const char* SCIPnodeselGetName(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->name;
}

/** gets description of node selector */
const char* SCIPnodeselGetDesc(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->desc;
}

/** gets priority of node selector in standard mode */
int SCIPnodeselGetStdPriority(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->stdpriority;
}

/** gets priority of node selector in standard mode */
void SCIPnodeselSetStdPriority(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the node selector */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   nodesel->stdpriority = priority;
   set->nodesel = NULL;
}

/** gets priority of node selector in memory saving mode */
int SCIPnodeselGetMemsavePriority(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->memsavepriority;
}

/** sets priority of node selector in memory saving mode */
void SCIPnodeselSetMemsavePriority(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the node selector */
   )
{
   assert(nodesel != NULL);
   assert(set != NULL);

   nodesel->memsavepriority = priority;
   set->nodesel = NULL;
}

/** gets user data of node selector */
SCIP_NODESELDATA* SCIPnodeselGetData(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->nodeseldata;
}

/** sets user data of node selector; user has to free old data in advance! */
void SCIPnodeselSetData(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_NODESELDATA*     nodeseldata         /**< new node selector user data */
   )
{
   assert(nodesel != NULL);

   nodesel->nodeseldata = nodeseldata;
}

/* new callback/method setter methods */

/** sets copy method of node selector */
void SCIPnodeselSetCopy(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy))    /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(nodesel != NULL);

   nodesel->nodeselcopy = nodeselcopy;
}

/** sets destructor method of node selector */
void SCIPnodeselSetFree(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELFREE ((*nodeselfree))    /**< destructor of node selector */
   )
{
   assert(nodesel != NULL);

   nodesel->nodeselfree = nodeselfree;
}

/** sets initialization method of node selector */
void SCIPnodeselSetInit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit))    /**< initialize node selector */
   )
{
   assert(nodesel != NULL);

   nodesel->nodeselinit = nodeselinit;
}

/** sets deinitialization method of node selector */
void SCIPnodeselSetExit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit))    /**< deinitialize node selector */
   )
{
   assert(nodesel != NULL);

   nodesel->nodeselexit = nodeselexit;
}

/** sets solving process initialization method of node selector */
void SCIPnodeselSetInitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINITSOL ((*nodeselinitsol))/**< solving process initialization method of node selector */
   )
{
   assert(nodesel != NULL);

   nodesel->nodeselinitsol = nodeselinitsol;
}

/** sets solving process deinitialization method of node selector */
void SCIPnodeselSetExitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXITSOL ((*nodeselexitsol))/**< solving process deinitialization method of node selector */
   )
{
   assert(nodesel != NULL);

   nodesel->nodeselexitsol = nodeselexitsol;
}

/** is node selector initialized? */
SCIP_Bool SCIPnodeselIsInitialized(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return nodesel->initialized;
}

/** enables or disables all clocks of \p nodesel, depending on the value of the flag */
void SCIPnodeselEnableOrDisableClocks(
   SCIP_NODESEL*         nodesel,            /**< the node selector for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the node selector be enabled? */
   )
{
   assert(nodesel != NULL);

   SCIPclockEnableOrDisable(nodesel->setuptime, enable);
   SCIPclockEnableOrDisable(nodesel->nodeseltime, enable);
}

/** gets time in seconds used in this node selector for setting up for next stages */
SCIP_Real SCIPnodeselGetSetupTime(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return SCIPclockGetTime(nodesel->setuptime);
}

/** gets time in seconds used in this node selector */
SCIP_Real SCIPnodeselGetTime(
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   assert(nodesel != NULL);

   return SCIPclockGetTime(nodesel->nodeseltime);
}
