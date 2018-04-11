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

/**@file   nodesel_restartdfs.c
 * @brief  node selector for depth first search with periodical selection of the best node
 * @author Tobias Achterberg
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/nodesel_restartdfs.h"


#define NODESEL_NAME             "restartdfs"
#define NODESEL_DESC             "depth first search with periodical selection of the best node"
#define NODESEL_STDPRIORITY       10000
#define NODESEL_MEMSAVEPRIORITY   50000


/*
 * Default parameter settings
 */

#define SELECTBESTFREQ              100 /**< frequency for selecting the best node instead of the deepest one */
#define COUNTONLYLEAVES            TRUE /**< only count leaf nodes or all nodes */


/** node selector data for best first search node selection */
struct SCIP_NodeselData
{
   SCIP_Longint          lastrestart;        /**< node number where the last best node was selected */
   SCIP_Longint          nprocessedleaves;   /**< number of processed leafs since the last restart */
   int                   selectbestfreq;     /**< frequency for selecting the best node instead of the deepest one */
   SCIP_Bool             countonlyleaves;    /**< only count leaf nodes or all nodes */
};


/*
 * Callback methods
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyRestartdfs)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeNodeselRestartdfs(scip) );

   return SCIP_OKAY;
}

/** destructor of node selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_NODESELFREE(nodeselFreeRestartdfs)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);

   /* free user data of node selector */
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   SCIPfreeBlockMemory(scip, &nodeseldata);
   SCIPnodeselSetData(nodesel, nodeseldata);

   return SCIP_OKAY;
}


/** solving process initialization method of node selector (called when branch and bound process is about to begin) */
static
SCIP_DECL_NODESELINITSOL(nodeselInitsolRestartdfs)
{
   SCIP_NODESELDATA* nodeseldata;

   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   /* reset counters */
   nodeseldata->lastrestart = 0;
   nodeseldata->nprocessedleaves = 0;

   return SCIP_OKAY;
}


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectRestartdfs)
{  /*lint --e{715}*/
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(selnode != NULL);

   /* decide if we want to select the node with lowest bound or the deepest node; finish the current dive in any case */
   *selnode = SCIPgetPrioChild(scip);
   if( *selnode == NULL )
   {
      SCIP_NODESELDATA* nodeseldata;
      SCIP_Longint nnodes;

      /* get node selector user data */
      nodeseldata = SCIPnodeselGetData(nodesel);
      assert(nodeseldata != NULL);

      /* increase the number of processed leafs since we are in a leaf */
      nodeseldata->nprocessedleaves++;

      nnodes = SCIPgetNNodes(scip);

      /* check if in case of "only leaves" the number processed leaves exceeds the frequency or in the other case the
       * number of processed node does it 
       */
      if( (nodeseldata->countonlyleaves && nodeseldata->nprocessedleaves >= nodeseldata->selectbestfreq) 
         || (!nodeseldata->countonlyleaves && nnodes - nodeseldata->lastrestart >= nodeseldata->selectbestfreq ) )
      {
         nodeseldata->lastrestart = nnodes;
         nodeseldata->nprocessedleaves = 0;
         *selnode = SCIPgetBestboundNode(scip);
      }
      else
      {
         *selnode = SCIPgetPrioSibling(scip);
         if( *selnode == NULL )
            *selnode = SCIPgetBestLeaf(scip);
      }
   }

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompRestartdfs)
{  /*lint --e{715}*/
   return (int)(SCIPnodeGetNumber(node2) - SCIPnodeGetNumber(node1));
}


/*
 * restartdfs specific interface methods
 */

/** creates the node selector for restarting depth first search and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselRestartdfs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODESEL* nodesel;

   /* allocate and initialize node selector data; this has to be freed in the destructor */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nodeseldata) );
   nodeseldata->lastrestart = 0;
   nodeseldata->nprocessedleaves = 0;
   nodeseldata->selectbestfreq = SELECTBESTFREQ;
   nodeseldata->countonlyleaves = COUNTONLYLEAVES;

   /* include node selector */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
          nodeselSelectRestartdfs, nodeselCompRestartdfs, nodeseldata) );

   assert(nodesel != NULL);

   SCIP_CALL( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyRestartdfs) );
   SCIP_CALL( SCIPsetNodeselFree(scip, nodesel, nodeselFreeRestartdfs) );
   SCIP_CALL( SCIPsetNodeselInitsol(scip, nodesel, nodeselInitsolRestartdfs) );

   /* add node selector parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "nodeselection/restartdfs/selectbestfreq",
         "frequency for selecting the best node instead of the deepest one",
         &nodeseldata->selectbestfreq, FALSE, SELECTBESTFREQ, 0, INT_MAX, NULL, NULL) );

   /* add node selector parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "nodeselection/restartdfs/countonlyleaves",
         "count only leaf nodes (otherwise all nodes)?",
         &nodeseldata->countonlyleaves, FALSE, COUNTONLYLEAVES, NULL, NULL) );

   return SCIP_OKAY;
}

