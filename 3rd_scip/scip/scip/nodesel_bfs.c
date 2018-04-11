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

/**@file   nodesel_bfs.c
 * @brief  node selector for best first search
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/nodesel_bfs.h"


#define NODESEL_NAME             "bfs"
#define NODESEL_DESC             "best first search"
#define NODESEL_STDPRIORITY      100000
#define NODESEL_MEMSAVEPRIORITY       0


/*
 * Default parameter settings
 */

#define MINPLUNGEDEPTH               -1 /**< minimal plunging depth, before new best node may be selected (-1 for dynamic setting) */
#define MAXPLUNGEDEPTH               -1 /**< maximal plunging depth, before new best node is forced to be selected (-1 for dynamic setting) */
#define MAXPLUNGEQUOT              0.25 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where plunging is performed */


/** node selector data for best first search node selection */
struct SCIP_NodeselData
{
   SCIP_Real             maxplungequot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where plunging is performed */
   int                   minplungedepth;     /**< minimal plunging depth, before new best node may be selected
                                              *   (-1 for dynamic setting) */
   int                   maxplungedepth;     /**< maximal plunging depth, before new best node is forced to be selected
                                              *   (-1 for dynamic setting) */
};


/*
 * Callback methods
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyBfs)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeNodeselBfs(scip) );

   return SCIP_OKAY;
}

/** destructor of node selector to free user data (called when SCIP is exiting) */
/**! [SnippetNodeselFreeBfs] */
static
SCIP_DECL_NODESELFREE(nodeselFreeBfs)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   /* free user data of node selector */
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   SCIPfreeBlockMemory(scip, &nodeseldata);
   SCIPnodeselSetData(nodesel, nodeseldata);

   return SCIP_OKAY;
}
/**! [SnippetNodeselFreeBfs] */


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectBfs)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;
   int minplungedepth;
   int maxplungedepth;
   int plungedepth;
   SCIP_Real maxplungequot;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   *selnode = NULL;

   /* get node selector user data */
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   /* calculate minimal and maximal plunging depth */
   minplungedepth = nodeseldata->minplungedepth;
   maxplungedepth = nodeseldata->maxplungedepth;
   maxplungequot = nodeseldata->maxplungequot;
   if( minplungedepth == -1 )
   {
      minplungedepth = SCIPgetMaxDepth(scip)/10;
      if( SCIPgetNStrongbranchLPIterations(scip) > 2*SCIPgetNNodeLPIterations(scip) )
        minplungedepth += 10;
      if( maxplungedepth >= 0 )
         minplungedepth = MIN(minplungedepth, maxplungedepth);
   }
   if( maxplungedepth == -1 )
      maxplungedepth = SCIPgetMaxDepth(scip)/2;
   maxplungedepth = MAX(maxplungedepth, minplungedepth);

   /* check, if we exceeded the maximal plunging depth */
   plungedepth = SCIPgetPlungeDepth(scip);
   if( plungedepth >= maxplungedepth )
   {
      /* we don't want to plunge again: select best node from the tree */
      SCIPdebugMsg(scip, "plungedepth: [%d,%d], cur: %d -> abort plunging\n", minplungedepth, maxplungedepth, plungedepth);
      *selnode = SCIPgetBestNode(scip);
      SCIPdebugMsg(scip, "  -> best node   : lower=%g\n",
         *selnode != NULL ? SCIPnodeGetLowerbound(*selnode) : SCIPinfinity(scip));
   }
   else
   {
      SCIP_NODE* node;
      SCIP_Real maxbound;

      /* check, if plunging is forced at the current depth */
      if( plungedepth < minplungedepth )
      {
         maxbound = SCIPinfinity(scip);
         SCIPdebugMsg(scip, "plungedepth: [%d,%d], cur: %d => maxbound: infinity\n",
            minplungedepth, maxplungedepth, plungedepth);
      }
      else
      {
         SCIP_Real lowerbound;
         SCIP_Real cutoffbound;
         /* get global lower and cutoff bound */
         lowerbound = SCIPgetLowerbound(scip);
         cutoffbound = SCIPgetCutoffbound(scip);

         /* if we didn't find a solution yet, the cutoff bound is usually very bad:
          * use only 20% of the gap as cutoff bound
          */
         if( SCIPgetNSolsFound(scip) == 0 )
            cutoffbound = lowerbound + 0.2 * (cutoffbound - lowerbound);
         /* calculate maximal plunging bound */
         maxbound = lowerbound + maxplungequot * (cutoffbound - lowerbound);

         SCIPdebugMsg(scip, "plungedepth: [%d,%d], cur: %d, bounds: [%g,%g], maxbound: %g\n",
            minplungedepth, maxplungedepth, plungedepth, lowerbound, cutoffbound, maxbound);         
      }

      /* we want to plunge again: prefer children over siblings, and siblings over leaves,
       * but only select a child or sibling, if its dual bound is small enough;
       * prefer using nodes with higher node selection priority assigned by the branching rule
       */
      node = SCIPgetPrioChild(scip);
      if( node != NULL && SCIPnodeGetLowerbound(node) < maxbound )
      {
         *selnode = node;
         SCIPdebugMsg(scip, "  -> selected prio child: lower=%g\n", SCIPnodeGetLowerbound(*selnode));
      }
      else
      {
         node = SCIPgetBestChild(scip);
         if( node != NULL && SCIPnodeGetLowerbound(node) < maxbound )
         {
            *selnode = node;
            SCIPdebugMsg(scip, "  -> selected best child: lower=%g\n", SCIPnodeGetLowerbound(*selnode));
         }
         else
         {
            node = SCIPgetPrioSibling(scip);
            if( node != NULL && SCIPnodeGetLowerbound(node) < maxbound )
            {
               *selnode = node;
               SCIPdebugMsg(scip, "  -> selected prio sibling: lower=%g\n", SCIPnodeGetLowerbound(*selnode));
            }
            else
            {
               node = SCIPgetBestSibling(scip);
               if( node != NULL && SCIPnodeGetLowerbound(node) < maxbound )
               {
                  *selnode = node;
                  SCIPdebugMsg(scip, "  -> selected best sibling: lower=%g\n", SCIPnodeGetLowerbound(*selnode));
               }
               else
               {
                  *selnode = SCIPgetBestNode(scip);
                  SCIPdebugMsg(scip, "  -> selected best leaf: lower=%g\n",
                     *selnode != NULL ? SCIPnodeGetLowerbound(*selnode) : SCIPinfinity(scip));
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompBfs)
{  /*lint --e{715}*/
   SCIP_Real lowerbound1;
   SCIP_Real lowerbound2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   lowerbound1 = SCIPnodeGetLowerbound(node1);
   lowerbound2 = SCIPnodeGetLowerbound(node2);
   if( SCIPisLT(scip, lowerbound1, lowerbound2) )
      return -1;
   else if( SCIPisGT(scip, lowerbound1, lowerbound2) )
      return +1;
   else
   {
      SCIP_Real estimate1;
      SCIP_Real estimate2;

      estimate1 = SCIPnodeGetEstimate(node1);
      estimate2 = SCIPnodeGetEstimate(node2);
      if( (SCIPisInfinity(scip,  estimate1) && SCIPisInfinity(scip,  estimate2)) ||
          (SCIPisInfinity(scip, -estimate1) && SCIPisInfinity(scip, -estimate2)) ||
          SCIPisEQ(scip, estimate1, estimate2) )
      {
         SCIP_NODETYPE nodetype1;
         SCIP_NODETYPE nodetype2;

         nodetype1 = SCIPnodeGetType(node1);
         nodetype2 = SCIPnodeGetType(node2);
         if( nodetype1 == SCIP_NODETYPE_CHILD && nodetype2 != SCIP_NODETYPE_CHILD )
            return -1;
         else if( nodetype1 != SCIP_NODETYPE_CHILD && nodetype2 == SCIP_NODETYPE_CHILD )
            return +1;
         else if( nodetype1 == SCIP_NODETYPE_SIBLING && nodetype2 != SCIP_NODETYPE_SIBLING )
            return -1;
         else if( nodetype1 != SCIP_NODETYPE_SIBLING && nodetype2 == SCIP_NODETYPE_SIBLING )
            return +1;
         else
         {
            int depth1;
            int depth2;

            depth1 = SCIPnodeGetDepth(node1);
            depth2 = SCIPnodeGetDepth(node2);
            if( depth1 < depth2 )
               return -1;
            else if( depth1 > depth2 )
               return +1;
            else
               return 0;
         }
      }

      if( SCIPisLT(scip, estimate1, estimate2) )
         return -1;

      assert(SCIPisGT(scip, estimate1, estimate2));
      return +1;
   }
}


/*
 * bfs specific interface methods
 */

/** creates the node selector for best first search and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselBfs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODESEL* nodesel;

   /* allocate and initialize node selector data; this has to be freed in the destructor */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nodeseldata) );

   /* include node selector */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
          nodeselSelectBfs, nodeselCompBfs, nodeseldata) );

   assert(nodesel != NULL);

   SCIP_CALL( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyBfs) );
   SCIP_CALL( SCIPsetNodeselFree(scip, nodesel, nodeselFreeBfs) );

   /* add node selector parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "nodeselection/bfs/minplungedepth",
         "minimal plunging depth, before new best node may be selected (-1 for dynamic setting)",
         &nodeseldata->minplungedepth, TRUE, MINPLUNGEDEPTH, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "nodeselection/bfs/maxplungedepth",
         "maximal plunging depth, before new best node is forced to be selected (-1 for dynamic setting)",
         &nodeseldata->maxplungedepth, TRUE, MAXPLUNGEDEPTH, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "nodeselection/bfs/maxplungequot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where plunging is performed",
         &nodeseldata->maxplungequot, TRUE, MAXPLUNGEQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

