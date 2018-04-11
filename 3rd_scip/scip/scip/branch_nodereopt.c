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

/**@file   branch_nodereopt.c
 * @brief  branching rule to reconstruct the search tree
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>

#include "scip/branch_nodereopt.h"
#include "scip/branch_relpscost.h"
#include "scip/cons_logicor.h"
#include "scip/scip.h"
#include "scip/tree.h"
#include "scip/pub_reopt.h"

#define BRANCHRULE_NAME            "nodereopt"
#define BRANCHRULE_DESC            "branching rule for node reoptimization"
#define BRANCHRULE_PRIORITY        -9000000
#define BRANCHRULE_MAXDEPTH            -1
#define BRANCHRULE_MAXBOUNDDIST         1.0

/*
 * Data structures
 */


/** execute the branching of nodes with additional constraints */
static
SCIP_RETCODE Exec(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_REOPTNODE* reoptnode;
   SCIP_NODE* curnode;
   SCIP_REOPTTYPE reopttype;
   SCIP_Bool localrestart;
   unsigned int* childids;
   unsigned int curid;
   int naddedconss;
   int nchilds;
   int childnodessize;
   int ncreatednodes;
   int c;

   assert(scip != NULL );
   assert(SCIPisReoptEnabled(scip));

   curnode = SCIPgetCurrentNode(scip);
   assert(curnode != NULL);

   curid = SCIPnodeGetReoptID(curnode);
   assert(curid >= 1 || SCIPgetRootNode(scip) == curnode);

   /* calculate local similarity and delete the induced subtree if the similarity is to low */
   localrestart = FALSE;
   SCIP_CALL( SCIPcheckReoptRestart(scip, curnode, &localrestart) );

   ncreatednodes = 0;

   if( localrestart )
   {
      *result = SCIP_DIDNOTRUN;
      goto TERMINATE;
   }

   SCIPdebugMsg(scip, "current node is %lld, ID %u:\n", SCIPnodeGetNumber(curnode), curid);

   /* get the corresponding node of the reoptimization tree */
   reoptnode = SCIPgetReoptnode(scip, curid);
   assert(reoptnode != NULL);
   reopttype = (SCIP_REOPTTYPE)SCIPreoptnodeGetType(reoptnode);


   /* The current node is equal to the root and dual reductions were performed. Since the root has a special role
    * within the reoptimiziation we have to split the root node into several nodes and move all stored child nodes to
    * the one representing the root node including all dual reductions as before.
    *
    * @note If the type is infsubtree, there cannot exist a child node and the method SCIPapplyReopt adds a global valid
    * constraint only.
    */
   if( curid == 0 )
   {
      if( reopttype == SCIP_REOPTTYPE_STRBRANCHED || reopttype == SCIP_REOPTTYPE_INFSUBTREE )
      {
         int ncreatedchilds;

         /* apply the reoptimization at the root node */
         SCIP_CALL( SCIPsplitReoptRoot(scip, &ncreatedchilds, &naddedconss) );

         if( reopttype == SCIP_REOPTTYPE_INFSUBTREE )
         {
            assert(ncreatedchilds == 0);
            assert(naddedconss == 1);

            /* there is nothing to do */
            *result = SCIP_DIDNOTRUN;

            goto TERMINATE;
         }

         assert(reopttype == SCIP_REOPTTYPE_STRBRANCHED);
         assert(ncreatedchilds >= 2);

         ncreatednodes += ncreatedchilds;

         /* We decrease the counter by one because after splitting the root node and moving all children to the node
          * representing the original root with all fixings (caused by dual reductions), we continue reactivating the
          * original children nodes of the root. Thus, the node containing all the fixings can be replaced by the children
          * nodes
          */
         --ncreatednodes;
      }

      goto REVIVE;
   }

   /* if we reach this part of the code the current has to be different to the root node */
   assert(curid >= 1);

  REVIVE:

   /* get the IDs of all child nodes */
   childnodessize = SCIPreoptnodeGetNChildren(reoptnode);
   SCIP_CALL( SCIPallocBufferArray(scip, &childids, childnodessize) );
   SCIP_CALL( SCIPgetReoptChildIDs(scip, curnode, childids, childnodessize, &nchilds) );

   if( childnodessize < nchilds )
   {
      childnodessize = SCIPreoptnodeGetNChildren(reoptnode);
      SCIP_CALL( SCIPreallocBufferArray(scip, &childids, childnodessize) );
      SCIP_CALL( SCIPgetReoptChildIDs(scip, curnode, childids, childnodessize, &nchilds) );
   }
   assert(nchilds <= childnodessize);

   naddedconss = 0;

   for(c = 0; c < nchilds; c++)
   {
      SCIP_NODE** childnodes;
      SCIP_Bool success;
      unsigned int childid;
      int ncreatedchilds;

      childid = childids[c];
      assert(childid >= 1);

      SCIPdebugMsg(scip, "process child at ID %u\n", childid);

      reoptnode = SCIPgetReoptnode(scip, childid);
      assert(reoptnode != NULL);

      reopttype = (SCIP_REOPTTYPE)SCIPreoptnodeGetType(reoptnode);
      ncreatedchilds = 0;

      /* check whether node need to be split */
      if( reopttype == SCIP_REOPTTYPE_STRBRANCHED || reopttype == SCIP_REOPTTYPE_INFSUBTREE )
      {
         /* by default we assume the node get split into two node (because using a constraint to split the node is
          * the default case */
         childnodessize = 2;
      }
      else
      {
         /* we only need to reconstruct the node */
         childnodessize = 1;
      }

      /* allocate buffer */
      SCIP_CALL( SCIPallocBufferArray(scip, &childnodes, childnodessize) );

      /* apply the reoptimization */
      SCIP_CALL( SCIPapplyReopt(scip, reoptnode, childid, SCIPnodeGetEstimate(curnode), childnodes, &ncreatedchilds,
            &naddedconss, childnodessize, &success) );

      if( !success )
      {
         assert(ncreatedchilds > childnodessize);

         /* reallocate buffer memory */
         childnodessize = ncreatedchilds+1;
         SCIP_CALL( SCIPreallocBufferArray(scip, &childnodes, childnodessize) );

         /* apply the reoptimization */
         SCIP_CALL( SCIPapplyReopt(scip, reoptnode, childid, SCIPnodeGetEstimate(curnode), childnodes, &ncreatedchilds,
               &naddedconss, childnodessize, &success) );
      }

      assert(success);

      /* free buffer memory */
      SCIPfreeBufferArray(scip, &childnodes);

      ncreatednodes += ncreatedchilds;
   }

   if( ncreatednodes == 0 )
      *result = SCIP_DIDNOTRUN;
   else
      *result = SCIP_BRANCHED;

   /* free the buffer memory */
   SCIPfreeBufferArray(scip, &childids);

  TERMINATE:

   SCIPdebugMsg(scip, "**** finish reoptimizing %d child nodes of node %lld ****\n", ncreatednodes, SCIPnodeGetNumber(curnode));

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyNodereopt)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleNodereopt(scip) );

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpNodereopt)
{/*lint --e{715}*/
   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   *result = SCIP_DIDNOTRUN;

   if( SCIPisReoptEnabled(scip) && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {
      SCIP_VAR** branchcands;
      SCIP_Real* branchcandssol;
      SCIP_Real* branchcandsfrac;
      SCIP_Real objsimrootlp;
      SCIP_Bool sbinit;
      int nbranchcands;

      assert(SCIPgetNReoptRuns(scip) > 1);

      SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/strongbranchinginit", &sbinit) );
      SCIP_CALL( SCIPgetRealParam(scip, "reoptimization/objsimrootLP", &objsimrootlp) );

      if( sbinit && SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip)
       && SCIPgetReoptSimilarity(scip, SCIPgetNReoptRuns(scip)-1, SCIPgetNReoptRuns(scip)) <= objsimrootlp ) /* check objsimrootlp */
      {
         /* get branching candidates */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &branchcands, &branchcandssol, &branchcandsfrac, NULL, &nbranchcands, NULL) );

         /* run strong branching initialization */
         if( nbranchcands > 0 )
         {
            SCIP_CALL( SCIPexecRelpscostBranching(scip, branchcands, branchcandssol, branchcandsfrac, nbranchcands, FALSE, result) );
            assert(*result == SCIP_DIDNOTRUN || *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM || *result == SCIP_CONSADDED);
         }
      }

      if( *result != SCIP_CUTOFF && *result != SCIP_REDUCEDDOM && *result != SCIP_CONSADDED )
      {
         assert((SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) == 0 && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
              || 1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

         SCIP_CALL( Exec(scip, result) );
      }
   }

   return SCIP_OKAY;
}

/** branching execution method for external candidates */
static SCIP_DECL_BRANCHEXECEXT(branchExecextNodereopt)
{/*lint --e{715}*/
   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   *result = SCIP_DIDNOTRUN;

   if( SCIPisReoptEnabled(scip) && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {
      assert((SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) == 0 && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
           || 1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

      SCIP_CALL( Exec(scip, result) );
   }

   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static SCIP_DECL_BRANCHEXECPS(branchExecpsNodereopt)
{/*lint --e{715}*/
   assert(branchrule != NULL );
   assert(*result != SCIP_BRANCHED);

   *result = SCIP_DIDNOTRUN;

   if( SCIPisReoptEnabled(scip) && SCIPreoptimizeNode(scip, SCIPgetCurrentNode(scip)) )
   {
      assert((SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)) == 0 && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
           || 1 <= SCIPnodeGetReoptID(SCIPgetCurrentNode(scip)));

      SCIP_CALL( Exec(scip, result) );
   }

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the nodereopt branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleNodereopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULE* branchrule;

   assert(scip != NULL );

   /* include nodereopt branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC,
         BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, NULL));

   assert(branchrule != NULL );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetBranchruleCopy(scip, branchrule, branchCopyNodereopt));
   SCIP_CALL(SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpNodereopt));
   SCIP_CALL(SCIPsetBranchruleExecExt(scip, branchrule, branchExecextNodereopt));
   SCIP_CALL(SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsNodereopt));

   return SCIP_OKAY;
}
