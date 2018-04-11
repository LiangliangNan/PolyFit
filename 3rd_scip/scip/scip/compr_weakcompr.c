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

/**@file   compr_weakcompr.c
 * @brief  weakcompr tree compression
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/mem.h"
#include "scip/misc.h"
#include "scip/compr_weakcompr.h"
#include "scip/compr.h"
#include "scip/pub_reopt.h"

#define COMPR_NAME             "weakcompr"
#define COMPR_DESC             "reduce the search frontier to k+1 or max{2, |C|+1} nodes."
#define COMPR_PRIORITY         1000
#define COMPR_MINNNODES        50

#define DEFAULT_MEM_REPR       2 /* since we cannot convert the added constraints into node currently, we choose 2 as default value */
/*
 * Data structures
 */

/** tree compression data */
struct SCIP_ComprData
{
   /* representative data */
   SCIP_REOPTNODE**      representatives;    /**< list of representatives */
   int                   nrepresentatives;   /**< number of representatives */
   int                   representativessize;/**< size of array representatives */
   SCIP_Bool             initialized;        /**< was compressor data initialized? */

   /* parameter */
   SCIP_Bool             convertconss;       /**< convert added logic-or constraints of size k into k nodes */
};


/*
 * Local methods
 */

/** sort the ids of child nodes by their dual bound of the last iteration */
static
SCIP_RETCODE sortIDs(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int*         childids,           /**< array of child ids */
   int                   nchildids           /**< first position */
   )
{
   SCIP_Real* lowerbounds;
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &lowerbounds, nchildids) );

   for( i = 0; i < nchildids; i++ )
   {
      lowerbounds[i] = SCIPreoptnodeGetLowerbound(SCIPgetReoptnode(scip, childids[i]));
   }

   /* sort the ids in decreasing order */
   SCIPsortDownRealInt(lowerbounds, (signed int*)childids, nchildids);

   /* free buffer memory */
   SCIPfreeBufferArray(scip, &lowerbounds);

   return SCIP_OKAY;
}

/** check of enough memory is allocated and reallocate of necessary. */
static
SCIP_RETCODE checkMemSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPRDATA*       comprdata,          /**< compression data */
   int                   nrepresentatives    /**< number of representatives */
   )
{
   assert(scip != NULL);
   assert(comprdata != NULL);

   if( comprdata->representativessize < nrepresentatives )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &comprdata->representatives, comprdata->representativessize, nrepresentatives) );
      comprdata->representativessize = nrepresentatives;
   }

   return SCIP_OKAY;
}

/** try to find a good representation
 *
 *  @todo implement:
 *      1) using k nodes without added constraint;
 *      2) resolve the added nods via some kind of interdiction branching
 */
static
SCIP_RETCODE constructCompression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< compression method */
   SCIP_COMPRDATA*       comprdata,          /**< compression data */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_NODE* currentnode;
   SCIP_VAR**** conss_var;
   SCIP_VAR*** vars;
   SCIP_Real*** conss_val;
   SCIP_Real** vals;
   SCIP_BOUNDTYPE** boundtypes;
   SCIP_BOUNDTYPE*** conss_boundtypes;
   int** conss_nvars;
   unsigned int* leaveids;
   int* nconss;
   int* nvars;
   int mem_vars;
   int nids;
   int nleaveids;
   int pos_repr_fix;
   int size;
   int k;
   int r;

   assert(scip != NULL);
   assert(comprdata != NULL);

   *result = SCIP_DIDNOTRUN;

   size = 1;
   currentnode = SCIPgetStage(scip) <= SCIP_STAGE_PRESOLVED ? NULL : SCIPgetCurrentNode(scip);

   if( SCIPgetStage(scip) <= SCIP_STAGE_PRESOLVED )
      nleaveids = SCIPgetNReoptLeaves(scip, currentnode);
   else
   {
      assert(currentnode != NULL);
      nleaveids = SCIPgetNReoptLeaves(scip, currentnode);
   }

   SCIPdebugMsg(scip, ">> start <%s> at node %llu (nleaves: %d, depth: %d)\n", COMPR_NAME,
      SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED ? 0 : SCIPnodeGetNumber(SCIPgetCurrentNode(scip)),
      nleaveids, SCIPgetStage(scip) <= SCIP_STAGE_PRESOLVED ? 0 : SCIPnodeGetDepth(currentnode));

   if( SCIPcomprGetMinNodes(compr) > nleaveids )
   {
      SCIPdebugMsg(scip, "-> skip compression (min. leaves = %d)\n", SCIPcomprGetMinNodes(compr));
      return SCIP_OKAY;
   }

   if( nleaveids == 0 )
   {
      SCIPdebugMsg(scip, "-> skip compression (k = %d, nleaves = %d)\n", 1, nleaveids);
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "-> try compression with %d node(s)\n", 1);

   *result = SCIP_DIDNOTFIND;

   /* collect the nodes to compress */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &leaveids, nleaveids) );

   SCIP_CALL( SCIPgetReoptLeaveIDs(scip, currentnode, leaveids, nleaveids, &nids) );
   assert(nids == nleaveids);

   /* sort the ids */
   SCIP_CALL( sortIDs(scip, leaveids, nleaveids) );

   /* allocate memory to store the old tree */

   mem_vars = 2*SCIPgetNVars(scip);

   /* we use normal instead of buffer memory because we may need to reallocate the 2-dimensional arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, size) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vals, size) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &boundtypes, size) );

   /* allocate buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &conss_var, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conss_val, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conss_boundtypes, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conss_nvars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nvars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nconss, size) );

   /* get data of nodes */
   for( k = size-1; k < 1; k++ )
   {
      SCIP_REOPTNODE* reoptnode;
      int mem_conss;
      int nvars2;
      int nafterdualvars;
      SCIPdebug(int c);

      /* we use normal instead of buffer memory because we may need to reallocate the 2-dimensional arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars[k], mem_vars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vals[k], mem_vars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &boundtypes[k], mem_vars) );

      /* get the branching path */
      reoptnode = SCIPgetReoptnode(scip, leaveids[k]);
      assert(reoptnode != NULL);

      SCIPgetReoptnodePath(scip, reoptnode, vars[k], vals[k], boundtypes[k], mem_vars, &nvars2, &nafterdualvars);
      assert(mem_vars >= nvars2 + nafterdualvars);

      nvars[k] = nvars2 + nafterdualvars;

      /* get the constraints */
      mem_conss = SCIPreoptnodeGetNConss(reoptnode);

      SCIP_CALL( SCIPallocBufferArray(scip, &conss_var[k], mem_conss) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &conss_val[k], mem_conss) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &conss_boundtypes[k], mem_conss) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &conss_nvars[k], mem_conss) ); /*lint !e866*/

      SCIPreoptnodeGetConss(reoptnode, conss_var[k], conss_val[k], conss_boundtypes[k], mem_conss, &nconss[k],
         conss_nvars[k]);
      assert(mem_conss == nconss[k]);

#ifdef SCIP_DEBUG
      for( c = 0; c < mem_conss; c++ )
         assert(conss_nvars[k][c] <= SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip));
#endif

#ifndef NDEBUG
      reoptnode = SCIPgetReoptnode(scip, leaveids[k]);
      assert(reoptnode != NULL);

      SCIPdebugMsg(scip, "-> use node at id %u, %d vars, %d conss, lowerbound = %.g\n", leaveids[k], nvars[k],
            SCIPreoptnodeGetNConss(reoptnode), SCIPreoptnodeGetLowerbound(reoptnode));
#endif
   }

   /* perform the compression */
   assert(comprdata->nrepresentatives == 0);

   pos_repr_fix = 1;

   /* calculate the number of representatives */
   comprdata->nrepresentatives = (nvars[0] > 0 ? 2 : 1);
   comprdata->nrepresentatives += nconss[0];

   /* check memory size */
   SCIP_CALL( checkMemSize(scip, comprdata, comprdata->nrepresentatives) );
   assert(comprdata->nrepresentatives <= comprdata->representativessize);

   /* initialize the representatives */
   SCIP_CALL( SCIPinitRepresentation(scip, comprdata->representatives, comprdata->nrepresentatives) );

   /* create 2 candidates for the fixed variables */
   if( nvars[0] >= 1 )
   {
      SCIP_Bool linear;
      int v;

      assert(pos_repr_fix < comprdata->nrepresentatives);

      linear = TRUE; /* todo: we have to adapt the compression to handle integer variables */

      /* create a representative at position 1 with fixed branching path */
      assert(SCIPreoptnodeGetNVars(comprdata->representatives[pos_repr_fix]) == 0);
      for( r = pos_repr_fix; r < comprdata->nrepresentatives; r++ )
      {
         /* copy the branching path to all representatives */
         assert(comprdata->representatives[r] != NULL);

         for( v = 0; v < nvars[0]; v++ )
         {
            SCIP_CALL( SCIPaddReoptnodeBndchg(scip, comprdata->representatives[r], vars[0][v],
                  vals[0][v], SCIPisFeasEQ(scip, vals[0][v], 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER) );
         }
      }

      /* create a representative at position 0 with an added constraint corresponding
       * to the branching path of the node
       */
      assert(comprdata->representatives[pos_repr_fix-1] != NULL);
      SCIP_CALL( SCIPaddReoptnodeCons(scip, comprdata->representatives[pos_repr_fix-1], vars[0], vals[0], boundtypes[k],
            1.0, SCIPinfinity(scip), nvars[0], REOPT_CONSTYPE_DUALREDS, linear) );

   }

   assert(0 <= pos_repr_fix && pos_repr_fix < comprdata->nrepresentatives);

   /* create nconss[0] nodes for the added constraints */
   for( k = 0; k < nconss[0]; k++ )
   {
      SCIP_Bool linear;
      int v;

      assert(pos_repr_fix < comprdata->nrepresentatives);

      linear = TRUE; /* todo: we have to adapt the compression to handle integer variables */

      /* create a node with fixed bounds corresponding to constraint at position k */

      /* fix the branching path */
      for( v = 0; v < conss_nvars[0][k]; v++ )
      {
         SCIP_CALL( SCIPaddReoptnodeBndchg(scip, comprdata->representatives[pos_repr_fix], conss_var[0][k][v],
               conss_val[0][k][v], SCIPisFeasEQ(scip, conss_val[0][k][v], 1.0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER) );
      }

      /* add this constraint to all further representatives */
      for( r = pos_repr_fix + 1; r < comprdata->nrepresentatives; r++ )
      {
         SCIP_CALL( SCIPaddReoptnodeCons(scip, comprdata->representatives[r], conss_var[0][k], conss_val[0][k],
               conss_boundtypes[0][k], 1.0, SCIPinfinity(scip), conss_nvars[0][k], REOPT_CONSTYPE_DUALREDS, linear) );
      }

      pos_repr_fix++;
   }

   /* @todo use more than one node */

   *result = SCIP_SUCCESS;

   SCIPdebugMsg(scip, "-> found representation of size %d.\n", comprdata->nrepresentatives);

   /* free memory */
   for( k = size-1; k >= 0; k-- )
   {
      SCIPfreeBufferArray(scip, &conss_nvars[k]);
      SCIPfreeBufferArray(scip, &conss_val[k]);
      SCIPfreeBufferArray(scip, &conss_var[k]);
      SCIPfreeBlockMemoryArray(scip, &boundtypes[k], mem_vars);
      SCIPfreeBlockMemoryArray(scip, &vals[k], mem_vars);
      SCIPfreeBlockMemoryArray(scip, &vars[k], mem_vars);
   }

   SCIPfreeBufferArray(scip, &nconss);
   SCIPfreeBufferArray(scip, &nvars);
   SCIPfreeBufferArray(scip, &conss_nvars);
   SCIPfreeBufferArray(scip, &conss_val);
   SCIPfreeBufferArray(scip, &conss_var);
   SCIPfreeBlockMemoryArray(scip, &boundtypes, size);
   SCIPfreeBlockMemoryArray(scip, &vals, size);
   SCIPfreeBlockMemoryArray(scip, &vars, size);

   SCIPfreeBlockMemoryArray(scip, &leaveids, nleaveids);

   return SCIP_OKAY;
}

/** apply the stored representation to the reopttree */
static
SCIP_RETCODE applyCompression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< compression method */
   SCIP_COMPRDATA*       comprdata,          /**< compression data */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_Bool success;
   int r;

   assert(scip != NULL);
   assert(compr != NULL);
   assert(comprdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( comprdata->nrepresentatives == 0 )
      return SCIP_OKAY;

   /* set references to the root node */
   for( r = 0; r < comprdata->nrepresentatives; r++ )
      SCIPreoptnodeSetParentID(comprdata->representatives[r], 0);

   success = FALSE;
   SCIP_CALL( SCIPsetReoptCompression(scip, comprdata->representatives, comprdata->nrepresentatives, &success) );

   if( success )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * Callback methods of tree compression
 */

/** copy method for tree compression plugins (called when SCIP copies plugins) */
static
SCIP_DECL_COMPRCOPY(comprCopyWeakcompr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(compr != NULL);
   assert(strcmp(SCIPcomprGetName(compr), COMPR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeComprWeakcompr(scip) );

   return SCIP_OKAY;
}

/** destructor of tree compression to free user data (called when SCIP is exiting) */
static
SCIP_DECL_COMPRFREE(comprFreeWeakcompr)
{
   SCIP_COMPRDATA* comprdata;

   assert(scip != NULL);
   assert(compr != NULL);

   comprdata = SCIPcomprGetData(compr);
   assert(comprdata != NULL);

   SCIPfreeBlockMemory(scip, &comprdata);
   SCIPcomprSetData(compr, NULL);

   return SCIP_OKAY;
}

/** deinitialization method of tree compression (called before transformed problem is freed) */
static
SCIP_DECL_COMPREXIT(comprExitWeakcompr)
{
   SCIP_COMPRDATA* comprdata;

   assert(scip != NULL);
   assert(compr != NULL);

   comprdata = SCIPcomprGetData(compr);
   assert(comprdata != NULL);

   if( comprdata->initialized )
   {
      int r;

      for( r = 0; r < comprdata->nrepresentatives; r++ )
      {
         SCIP_CALL( SCIPdeleteReoptnode(scip, &comprdata->representatives[r]) );
      }

      if( comprdata->representativessize > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &comprdata->representatives, comprdata->representativessize);
      }

      comprdata->representatives = NULL;
      comprdata->representativessize = 0;
      comprdata->nrepresentatives = 0;
      comprdata->initialized = FALSE;
   }

   return SCIP_OKAY;
}

/** execution method of tree compression */
static
SCIP_DECL_COMPREXEC(comprExecWeakcompr)
{
   SCIP_COMPRDATA* comprdata;

   assert(SCIPcomprIsInitialized(compr));

   comprdata = SCIPcomprGetData(compr);
   assert(comprdata != NULL);

   if( !comprdata->initialized )
   {
      SCIPdebugMsg(scip, ">> initializing <%s>\n", COMPR_NAME);

      comprdata->representativessize = DEFAULT_MEM_REPR;
      comprdata->nrepresentatives = 0;
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &comprdata->representatives, comprdata->representativessize) );
      comprdata->initialized = TRUE;
   }

   /* try to find a representation */
   SCIP_CALL( constructCompression(scip, compr, comprdata, result) );

   assert(*result == SCIP_DIDNOTRUN || *result == SCIP_DIDNOTFIND || *result == SCIP_SUCCESS);

   /* apply the representation, if some was found */
   if( *result == SCIP_SUCCESS )
   {
      SCIP_CALL( applyCompression(scip, compr, comprdata, result) );
      assert(*result == SCIP_DIDNOTRUN || *result == SCIP_SUCCESS);

      SCIPdebugMsg(scip, "->%s apply compression.\n", *result == SCIP_DIDNOTRUN ? " did not" : "");
   }

   return SCIP_OKAY;
}


/*
 * tree compression specific interface methods
 */

/** creates the weakcompr tree compression and includes it in SCIP */
SCIP_RETCODE SCIPincludeComprWeakcompr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_COMPRDATA* comprdata;
   SCIP_COMPR* compr;

   /* create weakcompr tree compression data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &comprdata) );
   assert(comprdata != NULL);
   comprdata->initialized = FALSE;

   /* include tree compression */
   SCIP_CALL( SCIPincludeComprBasic(scip, &compr,COMPR_NAME, COMPR_DESC, COMPR_PRIORITY, COMPR_MINNNODES,
         comprExecWeakcompr, comprdata) );

   assert(compr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetComprCopy(scip, compr, comprCopyWeakcompr) );
   SCIP_CALL( SCIPsetComprExit(scip, compr, comprExitWeakcompr) );
   SCIP_CALL( SCIPsetComprFree(scip, compr, comprFreeWeakcompr) );

   /* add weakcompr tree compression parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "compression/" COMPR_NAME "/convertconss", "convert constraints into nodes", &comprdata->convertconss, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
