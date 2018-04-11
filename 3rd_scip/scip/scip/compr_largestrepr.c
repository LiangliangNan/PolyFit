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

/**@file   compr_largestrepr.c
 * @brief  largestrepr tree compression
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/compr_largestrepr.h"
#include "scip/compr.h"
#include "scip/pub_reopt.h"


#define COMPR_NAME             "largestrepr"
#define COMPR_DESC             "heuristic searching for large common representatives"
#define COMPR_PRIORITY         2000
#define COMPR_MINNNODES        20

#define DEFAUL_MEM_REPR        10
#define DEFAULT_ITERS           5
#define DEFAULT_MINCOMMONVARS   3

/*
 * Data structures
 */

/** tree compression data */
struct SCIP_ComprData
{
   /* representative data */
   SCIP_REOPTNODE**      representatives;    /**< list of representatives */
   int                   nrepresentatives;   /**< number of representatives */
   int                   representativessize;/**< allocated memory for representatives */
   SCIP_Bool             initialized;        /**< was compressor data initialized? */

   /* statictics */
   SCIP_Real             rate;               /**< rate of compression */
   SCIP_Real             score;              /**< score of the best representation found */
   int                   nnodes;             /**< number of nodes after compressing */

   /* parameters */
   int                   mincomvars;         /**< minimal number of common variables */
   int                   niters;             /**< number of runs in the constrained part */
};


/*
 * Local methods
 */

/** calculate a signature of variables fixed to 0 and 1 by using binary shift
 *  and or operations. we calculate the signature on the basis of SCIPvarGetProbindex() % 64
 */
static
void calcSignature(
   SCIP_VAR**            vars,               /**< variable array */
   SCIP_Real*            vals,               /**< value array */
   int                   nvars,              /**< number of variables */
   SCIP_Longint*         signature0,         /**< pointer to store the signatures of variables fixed to 0 */
   SCIP_Longint*         signature1          /**< pointer to store the signatures of variables fixed to 1 */
   )
{
   int v;

   (*signature0) = 0;
   (*signature1) = 0;

   for( v = nvars - 1; v >= 0; --v )
   {
      if( vals[v] == 0 )
         (*signature0) |= ((unsigned int)1 << ((unsigned int)SCIPvarGetProbindex(vars[v]) % 64)); /*lint !e647*/
      else
         (*signature1) |= ((unsigned int)1 << ((unsigned int)SCIPvarGetProbindex(vars[v]) % 64)); /*lint !e647*/
   }

   return;
}

/** try to find a representation of the current search frontier.
 *
 *  We use the signatures of variables fixed to 0 and 1 to decide if there is
 *  definitely no intersection or if the intersection is potentially non-empty.
 *
 *  To find a good representation we start the procedure with a node and choose the best one.
 *  the heuristic tries to find a representation of size 2 in each iteration, i.e., runs in the
 *  constrained part.
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
   SCIP_VAR*** repvars;
   SCIP_Real** repvals;
   int* nrepvars;
   SCIP_VAR*** vars;
   SCIP_Real** vals;
   SCIP_BOUNDTYPE** bounds;
   SCIP_Real* lowerbounds;
   SCIP_Bool* covered;
   const char** varnames;
   SCIP_Real score;
   int nreps;
   SCIP_Longint* signature0;
   SCIP_Longint* signature1;
   int* common_vars;
   unsigned int* leaveids;
   int* nvars;
   int nids;
   int nleaveids;
   int k;
   int current_id;
   int start_id;

   assert(scip != NULL);
   assert(comprdata != NULL);

   *result = SCIP_DIDNOTRUN;

   assert(SCIPgetStage(scip) <= SCIP_STAGE_PRESOLVED);

   currentnode = NULL;
   nleaveids = SCIPgetNReoptLeaves(scip, currentnode);

   SCIPdebugMsg(scip, ">> start <%s> (nleaves: %d)\n", COMPR_NAME, nleaveids);

   if( SCIPcomprGetMinNodes(compr) > nleaveids )
   {
      SCIPdebugMsg(scip, "-> skip compression (min. leaves = %d)\n", SCIPcomprGetMinNodes(compr));
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "-> try compression with %d node(s)\n", nleaveids);

   /* collect the nodes to compress */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &leaveids, nleaveids) );

   SCIP_CALL( SCIPgetReoptLeaveIDs(scip, currentnode, leaveids, nleaveids, &nids) );
   assert(nids == nleaveids);

   /* allocate memory to store the old tree */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nleaveids) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nleaveids) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nleaveids) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nvars, nleaveids) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lowerbounds, nleaveids) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varnames, SCIPgetNOrigVars(scip)) );

   /* allocate memory to store the signatures */
   SCIP_CALL( SCIPallocBufferArray(scip, &signature0, nleaveids) );
   SCIP_CALL( SCIPallocBufferArray(scip, &signature1, nleaveids) );

   /* get data of nodes */
   for( k = 0; k < nleaveids; k++ )
   {
      SCIP_REOPTNODE* reoptnode;
      int mem_vars;
      int nvars2;
      int nafterdualvars;

      mem_vars = SCIPgetNBinVars(scip);

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars[k], mem_vars) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &vals[k], mem_vars) ); /*lint !e866*/
      SCIP_CALL( SCIPallocBufferArray(scip, &bounds[k], mem_vars) ); /*lint !e866*/

      /* get the branching path */
      reoptnode = SCIPgetReoptnode(scip, leaveids[k]);
      SCIPgetReoptnodePath(scip, reoptnode, vars[k], vals[k], bounds[k], mem_vars, &nvars2, &nafterdualvars);
      lowerbounds[k] = SCIPreoptnodeGetLowerbound(reoptnode);
      nvars[k] = nvars2 + nafterdualvars;

      /* calculate the signature*/
      calcSignature(vars[k], vals[k], nvars[k], &signature0[k], &signature1[k]);
   }

   for( start_id = 0; start_id < nleaveids; start_id++ )
   {
      nreps = -1;
      score = 0.0;

      /* we want to find an intersection that merges at least 3 nodes */
      if( nvars[start_id] <= comprdata->mincomvars + 1 )
         continue;

      /* initialize the covered-array with FALSE */
      SCIP_CALL( SCIPallocClearBufferArray(scip, &covered, nleaveids) );

      current_id = start_id;

      /* initialize storage for representatives */
      SCIP_CALL( SCIPallocBufferArray(scip, &repvars, 2+comprdata->niters) );
      SCIP_CALL( SCIPallocBufferArray(scip, &repvals, 2+comprdata->niters) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nrepvars, 2+comprdata->niters) );

      SCIPdebugMsg(scip, "+---+ start round %d +---+\n", start_id + 1);

      /* try to find common representatives */
      while( nreps-1 <= comprdata->niters && (nreps == -1 || (current_id % nleaveids) != start_id) )
      {
         int* idx_common_vars;
         int* idx_non_zero;
         int* covered_ids;
         int ncovered;
         int ncommon_vars;
         int nnon_zero_vars;
         int next_id;
         int nnemptyinters;
         int v;

         /* find the first id which is not covered */
         while( covered[current_id % nleaveids] && (nreps == -1 || (current_id % nleaveids) != start_id) )
            current_id++;

         current_id %= nleaveids;

         if( nreps > 0 && current_id == start_id )
            goto TERMINATE;

         /* if the this is the fist round with a new start_id we set the number of representatives to 0 */
         nreps = MAX(0, nreps);

         nnemptyinters = 0;

         /* mark the node as covered */
         covered[current_id] = TRUE;

         /* find the next not covered id */
         next_id = (current_id + 1) % nleaveids ;
         while( covered[next_id] && next_id != current_id )
            next_id = (next_id + 1) % nleaveids;

         if( next_id == current_id )
            goto TERMINATE;

         /* we want to find an intersection that merges at least 3 nodes */
         if( nvars[next_id] <= comprdata->mincomvars + 1 )
            continue;

         /* get a clear array of size SCIPgetNOrigVars */
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &common_vars, SCIPgetNOrigVars(scip)) );

         /* allocate buffer */
         nnon_zero_vars = 0;
         SCIP_CALL( SCIPallocBufferArray(scip, &idx_common_vars, nvars[current_id]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &idx_non_zero, nvars[current_id]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &covered_ids, nleaveids) );
         ncovered = 0;

         /* initialize common vars:
          *   vars[i] = 0 -> variable with idx i is not fixed
          *   vars[i] = 1 -> variable with idx i is fixed to zero
          *   vars[i] = 2 -> variable with idx i is fixed to one */
         for( v = 0; v < nvars[current_id]; v++ )
         {
            if( SCIPisFeasEQ(scip, vals[current_id][v], 0.0) )
               common_vars[SCIPvarGetProbindex(vars[current_id][v])] = 1;
            else
            {
               assert(SCIPisFeasEQ(scip, vals[current_id][v], 1.0));
               common_vars[SCIPvarGetProbindex(vars[current_id][v])] = 2;
            }

            varnames[SCIPvarGetProbindex(vars[current_id][v])] = SCIPvarGetName(vars[current_id][v]);

            idx_non_zero[nnon_zero_vars] = SCIPvarGetProbindex(vars[current_id][v]);
            nnon_zero_vars++;
         }

         SCIPdebugMsg(scip, "start with ID %u, %d fixed variables\n", leaveids[current_id], nnon_zero_vars);

         covered_ids[ncovered] = current_id;
         ncovered++;

         while( nnon_zero_vars >= comprdata->mincomvars )
         {
            /* find the next id which is not covered */
            while( covered[next_id % nleaveids] && (next_id % nleaveids) != current_id )
            {
               /* go to the next node if the intersection is empty */
               if( (signature0[current_id] & signature0[next_id % nleaveids]) == 0
                  && (signature1[current_id] & signature1[next_id % nleaveids]) == 0 )
                  next_id++;
               else
                  break;
            }

            if( (next_id % nleaveids) == current_id )
               break;

            next_id %= nleaveids;

            if( covered[next_id] )
               goto EMPTY;

            ncommon_vars = 0;

            /* calculate the intersection */
            for( v = 0; v < nvars[next_id]; v++ )
            {
               if( common_vars[SCIPvarGetProbindex(vars[next_id][v])] == (vals[next_id][v] == 0 ? 1 : 2) )
               {
                  idx_common_vars[ncommon_vars] = SCIPvarGetProbindex(vars[next_id][v]);
                  ncommon_vars++;
               }
            }

            /* the number of common variables should be at least mincomvars */
            if( ncommon_vars < comprdata->mincomvars )
               goto EMPTY;

            /* clear all non-zero entries which are not part of the intersection */
            for( v = 0; v < nnon_zero_vars; )
            {
               int v2;
               for( v2 = 0; v2 < ncommon_vars; v2++ )
               {
                  if( idx_non_zero[v] == idx_common_vars[v2] )
                     break;
               }

               /* the variable is not part of the intersection */
               if( v2 == ncommon_vars )
               {
                  common_vars[idx_non_zero[v]] = 0;

                  /* replace the idx with the last */
                  idx_non_zero[v] = idx_non_zero[nnon_zero_vars-1];
                  nnon_zero_vars--;
               }
               else
                  v++;
            }

            /* mark the node as covered */
            if( nnon_zero_vars > 0 )
            {
               covered[next_id] = TRUE;
               nnemptyinters++;

               SCIPdebugMessage("-> found intersection with ID %u, %d/%d common variables\n", leaveids[next_id],
                     nnon_zero_vars, MAX(nvars[current_id], nvars[next_id]));

               covered_ids[ncovered] = next_id;
               ncovered++;
            }

           EMPTY:
            next_id++;

            if( next_id % nleaveids == (current_id-1) % nleaveids )
               break;
         }

         if( nnemptyinters > 0 )
         {
            /* increase the number of representatives */
            if( nreps == 0 )
               nreps += 2;
            else
               nreps++;

            SCIP_CALL( SCIPallocBufferArray(scip, &repvars[nreps-2], nnon_zero_vars) ); /*lint !e866*/
            SCIP_CALL( SCIPallocBufferArray(scip, &repvals[nreps-2], nnon_zero_vars) ); /*lint !e866*/
            nrepvars[nreps-2] = nnon_zero_vars;

            /* set the common variables and bounds (use non-zero idx)*/
            for( k = 0; k < nnon_zero_vars; k++ )
            {
               repvars[nreps-2][k] = SCIPfindVar(scip, varnames[idx_non_zero[k]]);
               repvals[nreps-2][k] = common_vars[idx_non_zero[k]] - 1;
               assert(repvals[nreps-2][k] == 0 || repvals[nreps-2][k] == 1);
            }
         }
         else
         {
            SCIPdebugMsg(scip, "-> did not found a intersection larger than %d\n", comprdata->mincomvars);
            covered[current_id] = FALSE;
         }

         /* calculate the score */
         score += (SCIP_Real) ncovered * nnon_zero_vars;

         SCIPdebugMessage("-> current representation is of size %d with score = %.1f\n", nreps, score);

         current_id = (current_id + 1) % nleaveids;

         /* free memory */
         SCIPfreeBlockMemoryArray(scip, &common_vars, SCIPgetNOrigVars(scip));

         SCIPfreeBufferArray(scip, &covered_ids);
         SCIPfreeBufferArray(scip, &idx_non_zero);
         SCIPfreeBufferArray(scip, &idx_common_vars);
      }

     TERMINATE:

      /* add the number of variables of uncovered nodes to the loss of information */
      SCIPdebugMessage("-> final representation is of size %d with score = %.1f\n", nreps, score);

      /* We found a better representation, i.e., with less loss of information.
       * 1. reset the previous represenation
       * 2. check if we need to reallocate the memory
       * 3. set the new representation
       */
      if( SCIPisSumGT(scip, score, comprdata->score) )
      {
         /* reset the current representation */
         SCIP_CALL( SCIPresetRepresentation(scip, comprdata->representatives, comprdata->nrepresentatives) );

         /* ensure that enough memory is allocated */
         if( comprdata->representativessize < nreps )
         {
            /* free the representatoin */
            SCIP_CALL( SCIPfreeRepresentation(scip, comprdata->representatives, comprdata->nrepresentatives) );

            /* reallocate memory */
            SCIP_CALL( SCIPreallocMemoryArray(scip, &comprdata->representatives, nreps) );
            comprdata->representativessize = nreps;

            /* initialize the representation */
            SCIP_CALL( SCIPinitRepresentation(scip, comprdata->representatives, comprdata->representativessize) );
         }

         /* set the new representation
          *
          * copy the new representation. we skip the last representative because it is implicitly given by the union of
          * the logic-or constraints of all previous representatives.
          */
         comprdata->score = score;
         comprdata->nrepresentatives = nreps;

         for( k = 0; k < nreps-1; k++ )
         {
            int v;

            for( v = 0; v < nrepvars[k]; v++ )
            {
               SCIP_CALL( SCIPaddReoptnodeBndchg(scip, comprdata->representatives[k], repvars[k][v],
                     repvals[k][v], repvals[k][v] == 0 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER) );
            }
         }

         *result = SCIP_SUCCESS;
      }

      /* free representatives storage */
      for( k = 0; k <= nreps-2; k++ )
      {
         SCIPfreeBufferArray(scip, &repvals[k]);
         SCIPfreeBufferArray(scip, &repvars[k]);
      }

      SCIPfreeBufferArray(scip, &nrepvars);
      SCIPfreeBufferArray(scip, &repvals);
      SCIPfreeBufferArray(scip, &repvars);

      /* free covered array */
      SCIPfreeBufferArray(scip, &covered);
   }

   /* check if we have found a representation and construct the missing constraints */
   if( comprdata->nrepresentatives > 0 )
   {
      SCIPdebugMessage("best representation found has %d leaf nodes and score = %g\n", comprdata->nrepresentatives, comprdata->score);

      /* iterate over all representatives */
      for( k = 0; k < comprdata->nrepresentatives-1; k++ )
      {
         int r;

         /* add a constraint (corresponding to the branching path of k) to all representatives
          * in the subtree induced by the sibling of k
          */
         for( r = k + 1; r < comprdata->nrepresentatives; r++ )
         {
            SCIP_VAR** pathvars;
            SCIP_Real* pathvals;
            SCIP_BOUNDTYPE* pathboundtypes;
            SCIP_Real lhs;
            SCIP_Bool linear;
            int pathvarssize;
            int npathvars;
            int npathafterdualvars;
            int i;

            pathvarssize = SCIPreoptnodeGetNVars(comprdata->representatives[k]);

            /* allocate buffer */
            SCIP_CALL( SCIPallocBufferArray(scip, &pathvars, pathvarssize) );
            SCIP_CALL( SCIPallocBufferArray(scip, &pathvals, pathvarssize) );
            SCIP_CALL( SCIPallocBufferArray(scip, &pathboundtypes, pathvarssize) );

            /* get the stored path */
            SCIPgetReoptnodePath(scip, comprdata->representatives[k], pathvars, pathvals, pathboundtypes, pathvarssize,
                  &npathvars, &npathafterdualvars);

            lhs = 1.0;
            linear = TRUE; /* todo: we have to adapt the compression to handle integer variables */

            /* negate the branching path */
            for( i = 0; i < npathvars; i++ )
            {
               assert(SCIPvarIsOriginal(pathvars[i]));

               /* we have to construct a linear constraint that can be upgraded to a logic-or constraint
                *
                * each variable i with pathvals[i] == 0 and pathboundtypes[i] == SCIP_BOUNDTYPE_UPPER needs a coefficient
                * of 1.0, all remaining variables (i.e., pathvals[i] == 1 and pathboundtypes[i] == SCIP_BOUNDTYPE_LOWER)
                * need a -1.0 and we have to reduce the lhs by -1.0.
                *
                *        sum_{i : pathvals[i] == 0.0} x_i + sum_{j : pathvals[j] == 1.0} (1.0-x_{j}) >= 1.0
                *  <==>  sum_{i : pathvals[i] == 0.0} x_i + sum_{j : pathvals[j] == 1.0}     -x_{j}  >= 1.0 - sum_{j : pathvals[j] == 1.0} 1.0
                */
               if( SCIPisEQ(scip, pathvals[i], 0.0) )
               {
                  assert(pathboundtypes[i] == SCIP_BOUNDTYPE_UPPER);

                  pathvals[i] = 1.0;
               }
               else
               {
                  assert(SCIPisEQ(scip, pathvals[i], 1.0));
                  assert(pathboundtypes[i] == SCIP_BOUNDTYPE_LOWER);

                  pathvals[i] = -1.0;
                  lhs -= 1.0;
               }
            }

            SCIP_CALL( SCIPaddReoptnodeCons(scip, comprdata->representatives[r], pathvars, pathvals, NULL, lhs,
                  SCIPinfinity(scip), npathvars, REOPT_CONSTYPE_DUALREDS, linear) );

            /* free buffer */
            SCIPfreeBufferArray(scip, &pathboundtypes);
            SCIPfreeBufferArray(scip, &pathvals);
            SCIPfreeBufferArray(scip, &pathvars);
         }
      }
   }

   /* free memory */
   for( k = nleaveids-1; k >= 0; k-- )
   {
      SCIPfreeBufferArray(scip, &bounds[k]);
      SCIPfreeBufferArray(scip, &vals[k]);
      SCIPfreeBufferArray(scip, &vars[k]);
   }

   SCIPfreeBufferArray(scip, &signature1);
   SCIPfreeBufferArray(scip, &signature0);

   SCIPfreeBufferArray(scip, &varnames);
   SCIPfreeBufferArray(scip, &lowerbounds);
   SCIPfreeBufferArray(scip, &nvars);
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   SCIPfreeBlockMemoryArray(scip, &leaveids, nleaveids);

   return SCIP_OKAY;
}

/** apply the found representation to the reopttree. */
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

   SCIP_CALL( SCIPfreeRepresentation(scip, comprdata->representatives, comprdata->representativessize) );

   if( success )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * Callback methods of tree compression
 */

/** copy method for tree compression plugins (called when SCIP copies plugins) */
static
SCIP_DECL_COMPRCOPY(comprCopyLargestrepr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(compr != NULL);
   assert(strcmp(SCIPcomprGetName(compr), COMPR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeComprLargestrepr(scip) );

   return SCIP_OKAY;
}

/** destructor of tree compression to free user data (called when SCIP is exiting) */
static
SCIP_DECL_COMPRFREE(comprFreeLargestrepr)
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
SCIP_DECL_COMPREXIT(comprExitLargestrepr)
{
   SCIP_COMPRDATA* comprdata;

   assert(scip != NULL);
   assert(compr != NULL);

   comprdata = SCIPcomprGetData(compr);
   assert(comprdata != NULL);

   if( comprdata->initialized )
   {
      if( comprdata->representativessize > 0 )
      {
         SCIPfreeMemoryArray(scip, &comprdata->representatives);
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
SCIP_DECL_COMPREXEC(comprExecLargestrepr)
{
   SCIP_COMPRDATA* comprdata;

   comprdata = SCIPcomprGetData(compr);
   assert(comprdata != NULL);

   if( !comprdata->initialized )
   {
      SCIPdebugMsg(scip, ">> initializing <%s>\n", COMPR_NAME);

      comprdata->representativessize = DEFAUL_MEM_REPR;
      comprdata->nrepresentatives = 0;
      comprdata->rate = 0.0;
      comprdata->score = 0.0;
      comprdata->nnodes = 0;
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &comprdata->representatives, comprdata->representativessize) );

      /* initialize the representation */
      SCIP_CALL( SCIPinitRepresentation(scip, comprdata->representatives, comprdata->representativessize) );

      comprdata->initialized = TRUE;
   }

   *result = SCIP_DIDNOTRUN;

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
   else
   {
      SCIP_CALL( SCIPfreeRepresentation(scip, comprdata->representatives, comprdata->representativessize) );
   }

   return SCIP_OKAY;
}

/*
 * tree compression specific interface methods
 */

/** creates the largestrepr tree compression and includes it in SCIP */
SCIP_RETCODE SCIPincludeComprLargestrepr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_COMPRDATA* comprdata;
   SCIP_COMPR* compr;

   /* create largestrepr tree compression data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &comprdata) );
   assert(comprdata != NULL);
   comprdata->initialized = FALSE;

   /* include tree compression */
   SCIP_CALL( SCIPincludeComprBasic(scip, &compr, COMPR_NAME, COMPR_DESC, COMPR_PRIORITY, COMPR_MINNNODES,
         comprExecLargestrepr, comprdata) );

   assert(compr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetComprCopy(scip, compr, comprCopyLargestrepr) );
   SCIP_CALL( SCIPsetComprExit(scip, compr, comprExitLargestrepr) );
   SCIP_CALL( SCIPsetComprFree(scip, compr, comprFreeLargestrepr) );

   /* add largestrepr tree compression parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "compression/" COMPR_NAME "/iterations", "number of runs in the constrained part.",
         &comprdata->niters, FALSE, DEFAULT_ITERS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "compression/" COMPR_NAME "/mincommonvars", "minimal number of common variables.",
         &comprdata->mincomvars, FALSE, DEFAULT_MINCOMMONVARS, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
