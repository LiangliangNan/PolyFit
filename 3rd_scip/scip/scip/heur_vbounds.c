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

/**@file   heur_vbounds.c
 * @brief  LNS heuristic uses the variable lower and upper bounds to determine the search neighborhood
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Jens Schulz
 * @author Gerald Gamrath
 *
 * @todo allow smaller fixing rate for probing LP?
 * @todo allow smaller fixing rate after presolve if total number of variables is small (<= 1000)?
 *
 * More details about the heuristic can be found in@n
 * Structure-Based Primal Heuristics for Mixed Integer Programming@n
 * Gerald Gamrath, Timo Berthold, Stefan Heinz, and Michael Winkler@n
 * Optimization in the Real World, Volume 13 of the series Mathematics for Industry, pp 37-53@n
 * Preliminary version available as <a href="https://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/5551">ZIB-Report 15-26</a>.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_vbounds.h"
#include "scip/heur_locks.h"

#ifdef SCIP_STATISTIC
#include "scip/clock.h"
#endif

#define VBOUNDVARIANT_NOOBJ      0x001u
#define VBOUNDVARIANT_BESTBOUND  0x002u
#define VBOUNDVARIANT_WORSTBOUND 0x004u

#define HEUR_NAME             "vbounds"
#define HEUR_DESC             "LNS heuristic uses the variable lower and upper bounds to determine the search neighborhood"
#define HEUR_DISPCHAR         'V'
#define HEUR_PRIORITY         2500
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL    /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MININTFIXINGRATE 0.65    /**< minimum percentage of integer variables that have to be fixed */
#define DEFAULT_MINMIPFIXINGRATE 0.65    /**< minimuskipobjm percentage of variables that have to be fixed within sub-SCIP
                                         *   (integer and continuous) */
#define DEFAULT_MINIMPROVE    0.01      /**< factor by which vbounds heuristic should at least improve the
                                         *   incumbent */
#define DEFAULT_MINNODES      500LL     /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS      500LL     /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT     0.1       /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_MAXPROPROUNDS 2         /**< maximum number of propagation rounds during probing */
#define DEFAULT_MAXBACKTRACKS 10        /**< maximum number of backtracks during the fixing process */
#define DEFAULT_COPYCUTS      TRUE      /**< should all active cuts from the cutpool of the original scip be copied to
                                         *   constraints of the subscip? */
#define DEFAULT_USELOCKFIXINGS FALSE    /**< should more variables be fixed based on variable locks if
                                         *   the fixing rate was not reached?
                                         */

/** which variants of the vbounds heuristic that try to stay feasible should be called? */
#define DEFAULT_FEASVARIANT   (VBOUNDVARIANT_BESTBOUND | VBOUNDVARIANT_WORSTBOUND)

/** which tightening variants of the vbounds heuristic should be called? */
#define DEFAULT_TIGHTENVARIANT   (VBOUNDVARIANT_NOOBJ | VBOUNDVARIANT_BESTBOUND | VBOUNDVARIANT_WORSTBOUND)


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_VAR**            vbvars;             /**< topological sorted variables with respect to the variable bounds */
   SCIP_BOUNDTYPE*       vbbounds;           /**< topological sorted variables with respect to the variable bounds */
   int                   nvbvars;            /**< number of variables in variable lower bound array */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Longint          usednodes;          /**< nodes already used by vbounds heuristic in earlier calls */
   SCIP_Real             minintfixingrate;   /**< minimum percentage of integer variables that have to be fixed */
   SCIP_Real             minmipfixingrate;   /**< minimum percentage of variables that have to be fixed within sub-SCIP
                                              *   (integer and continuous) */
   SCIP_Real             minimprove;         /**< factor by which vbounds heuristic should at least improve the incumbent */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             cutoffbound;
   int                   maxproprounds;      /**< maximum number of propagation rounds during probing */
   int                   maxbacktracks;      /**< maximum number of backtracks during the fixing process */
   int                   feasvariant;        /**< which variants of the vbounds heuristic that try to stay feasible
                                              *   should be called? */
   int                   tightenvariant;     /**< which tightening variants of the vbounds heuristic should be called? */
   SCIP_Bool             initialized;        /**< is the candidate list initialized? */
   SCIP_Bool             applicable;         /**< is the heuristic applicable? */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   subproblem? */
   SCIP_Bool             uselockfixings;     /**< should more variables be fixed based on variable locks if
                                              *   the fixing rate was not reached?
                                              */


};

/**@name Heuristic defines
 *
 * @{
 *
 * The heuristic works on indices representing a bound of a variable. This index will be called bound index in the
 * following. For a given active variable with problem index i (note that active variables have problem indices
 * between 0 and nactivevariable - 1), the bound index of its lower bound is 2*i, the bound index of its upper
 * bound is 2*i + 1. The other way around, a given bound index i corresponds to the variable with problem index
 * i/2 (rounded down), and to the lower bound, if i is even, to the upper bound if i is odd.
 * The following macros can be used to convert bound index into variable problem index and boundtype and vice versa.
 */
#define getLbIndex(idx) (2*(idx))
#define getUbIndex(idx) (2*(idx)+1)
#define getVarIndex(idx) ((idx)/2)
#define getBoundtype(idx) (((idx) % 2 == 0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER)
#define isIndexLowerbound(idx) ((idx) % 2 == 0)
#define getOtherBoundIndex(idx) (((idx) % 2 == 0) ? (idx) + 1 : (idx) - 1)


/*
 * Local methods
 */

/** reset heuristic data structure */
static
void heurdataReset(
   SCIP_HEURDATA*        heurdata            /**< structure containing heurdata */
   )
{
   heurdata->vbvars = NULL;
   heurdata->vbbounds = NULL;
   heurdata->nvbvars = 0;
   heurdata->initialized = FALSE;
   heurdata->applicable = FALSE;
}


/** performs depth-first-search in the implicitly given directed graph from the given start index */
static
SCIP_RETCODE dfs(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   startnode,          /**< node to start the depth-first-search */
   SCIP_Shortbool*       visited,            /**< array to store for each node, whether it was already visited */
   int*                  dfsstack,           /**< array of size number of nodes to store the stack;
                                              *   only needed for performance reasons */
   int*                  stacknextedge,      /**< array of size number of nodes to store the number of adjacent nodes
                                              *   already visited for each node on the stack; only needed for
                                              *   performance reasons */
   int*                  stacknextcliquevar, /**< array of size number of nodes to store the number of variables
                                              *   already evaluated for the clique currently being evaluated */
   int*                  cliqueexit,         /**< exit node when entering a clique */
   int*                  dfsnodes,           /**< array of nodes that can be reached starting at startnode, in reverse
                                              *   dfs order */
   int*                  ndfsnodes           /**< pointer to store number of nodes that can be reached starting at
                                              *   startnode */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_VAR** vbvars;
   SCIP_Real* coefs;
   SCIP_Bool lower;
   SCIP_Bool found;
   int maxstacksize;
   int stacksize;
   int curridx;
   int idx;
   int nvbvars;
   int i;

   assert(startnode >= 0);
   assert(startnode < 2 * SCIPgetNVars(scip));
   assert(visited != NULL);
   assert(visited[startnode] == FALSE);
   assert(dfsstack != NULL);
   assert(dfsnodes != NULL);
   assert(ndfsnodes != NULL);

   vars = SCIPgetVars(scip);

   /* put start node on the stack */
   dfsstack[0] = startnode;
   stacknextcliquevar[0] = 0;
   stacknextedge[0] = 0;
   maxstacksize = 1;
   stacksize = 1;
   idx = -1;

   /* we run until no more bounds indices are on the stack */
   while( stacksize > 0 )
   {
      /* get next node from stack */
      curridx = dfsstack[stacksize - 1];

      /* mark current node as visited */
      assert(visited[curridx] == (stacknextedge[stacksize - 1] != 0));
      visited[curridx] = TRUE;
      found = FALSE;

      startvar = vars[getVarIndex(curridx)];
      lower = isIndexLowerbound(curridx);

      if( stacknextedge[stacksize - 1] >= 0 )
      {
         /* go over edges corresponding to varbounds */
         if( lower )
         {
            vbvars = SCIPvarGetVlbVars(startvar);
            coefs = SCIPvarGetVlbCoefs(startvar);
            nvbvars = SCIPvarGetNVlbs(startvar);
         }
         else
         {
            vbvars = SCIPvarGetVubVars(startvar);
            coefs = SCIPvarGetVubCoefs(startvar);
            nvbvars = SCIPvarGetNVubs(startvar);
         }

         /* iterate over all vbounds for the given bound */
         for( i = stacknextedge[stacksize - 1]; i < nvbvars; ++i )
         {
            if( !SCIPvarIsActive(vbvars[i]) )
               continue;

            idx = (SCIPisPositive(scip, coefs[i]) == lower) ? getLbIndex(SCIPvarGetProbindex(vbvars[i])) : getUbIndex(SCIPvarGetProbindex(vbvars[i]));
            assert(idx >= 0);

            /* break when the first unvisited node is reached */
            if( !visited[idx] )
               break;
         }

         /* we stopped because we found an unhandled node and not because we reached the end of the list */
         if( i < nvbvars )
         {
            assert(!visited[idx]);

            /* put the adjacent node onto the stack */
            dfsstack[stacksize] = idx;
            stacknextedge[stacksize] = 0;
            stacknextcliquevar[stacksize] = 0;
            stacknextedge[stacksize - 1] = i + 1;
            stacksize++;
            assert(stacksize <= 2* SCIPgetNVars(scip));

            /* restart while loop, get next index from stack */
            continue;
         }
      }

      stacknextedge[stacksize - 1] = -1;

      /* treat cliques */
      if( SCIPvarIsBinary(startvar) )
      {
         SCIP_CLIQUE** cliques = SCIPvarGetCliques(startvar, !lower);
         int ncliques = SCIPvarGetNCliques(startvar, !lower);
         int j;

         /* iterate over all not yet handled cliques and search for an unvisited node */
         for( j = -stacknextedge[stacksize - 1] - 1; j < ncliques; ++j )
         {
            SCIP_VAR** cliquevars;
            SCIP_Bool* cliquevals;
            int ncliquevars;

            /* the first time we evaluate this clique for the current node */
            if( stacknextcliquevar[stacksize - 1] == 0 )
            {
               if( cliqueexit[SCIPcliqueGetIndex(cliques[j])] > 0 )
               {
                  if( !visited[cliqueexit[SCIPcliqueGetIndex(cliques[j])] - 1] &&
                     cliqueexit[SCIPcliqueGetIndex(cliques[j])] - 1 != curridx )
                  {
                     stacknextedge[stacksize - 1] = -j - 2;
                     stacknextcliquevar[stacksize - 1] = 0;
                     idx = cliqueexit[SCIPcliqueGetIndex(cliques[j])] - 1;
                     cliqueexit[SCIPcliqueGetIndex(cliques[j])] = -1;
                     found = TRUE;
                  }
                  else
                     continue;
               }
               else if( cliqueexit[SCIPcliqueGetIndex(cliques[j])] == 0 )
               {
                  cliqueexit[SCIPcliqueGetIndex(cliques[j])] = getOtherBoundIndex(curridx) + 1;
               }
               else
                  continue;
            }
            if( !found )
            {
               cliquevars = SCIPcliqueGetVars(cliques[j]);
               cliquevals = SCIPcliqueGetValues(cliques[j]);
               ncliquevars = SCIPcliqueGetNVars(cliques[j]);

               for( i = 0; i < ncliquevars; ++i )
               {
                  if( cliquevars[i] == startvar )
                     continue;

                  if( SCIPvarGetIndex(cliquevars[i]) < 0 )
                     continue;

                  if( cliquevals[i] )
                     idx = getLbIndex(SCIPvarGetProbindex(cliquevars[i]));
                  else
                     idx = getUbIndex(SCIPvarGetProbindex(cliquevars[i]));

                  assert(idx >= 0 && idx < 2 * SCIPgetNVars(scip));

                  /* break when the first unvisited node is reached */
                  if( idx >= 0 && !visited[idx] )
                  {
                     if( i < ncliquevars - 1 )
                     {
                        stacknextedge[stacksize - 1] = -j - 1;
                        stacknextcliquevar[stacksize - 1] = i + 1;
                     }
                     else
                     {
                        stacknextedge[stacksize - 1] = -j - 2;
                        stacknextcliquevar[stacksize - 1] = 0;
                     }
                     found = TRUE;
                     break;
                  }
               }
            }
            if( found )
            {
               assert(!visited[idx]);

               /* put the adjacent node onto the stack */
               dfsstack[stacksize] = idx;
               stacknextedge[stacksize] = 0;
               stacknextcliquevar[stacksize] = 0;
               stacksize++;
               assert(stacksize <= 2* SCIPgetNVars(scip));

               break;
            }
         }
         /* restart while loop, get next index from stack */
         if( found )
            continue;
      }

      maxstacksize = MAX(maxstacksize, stacksize);

      /* the current node was completely handled, remove it from the stack */
      stacksize--;

      if( (maxstacksize > 1) && SCIPvarGetType(startvar) != SCIP_VARTYPE_CONTINUOUS )
      {
         /* store node in the sorted nodes array */
         dfsnodes[(*ndfsnodes)] = curridx;
         (*ndfsnodes)++;
      }
      else
         visited[curridx] = FALSE;
   }

   return SCIP_OKAY;
}


/** sort the bounds of variables topologically */
static
SCIP_RETCODE topologicalSort(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  vbvars,             /**< array to store variable bounds in topological order */
   int*                  nvbvars             /**< pointer to store number of variable bounds in the graph */
   )
{
   int* dfsstack;
   int* stacknextedge;
   int* stacknextcliquevar;
   int* cliqueexit;
   SCIP_Shortbool* visited;
   int nbounds;
   int i;

   assert(scip != NULL);

   nbounds = 2 * SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &dfsstack, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stacknextedge, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stacknextcliquevar, nbounds) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &cliqueexit, SCIPgetNCliques(scip)) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &visited, nbounds) );


   /* while there are unvisited nodes, run dfs on the inverse graph starting from one of these nodes; the dfs orders are
    * stored in the topoorder array, later dfs calls are just appended after the stacks of previous dfs calls, which
    * gives us a topological order
    */
   for( i = 0; i < nbounds; ++i )
   {
      if( !visited[i] )
      {
         SCIP_CALL( dfs(scip, i, visited, dfsstack, stacknextedge, stacknextcliquevar, cliqueexit, vbvars, nvbvars) );
      }
   }
   assert(*nvbvars <= nbounds);

   SCIPfreeBufferArray(scip, &visited);
   SCIPfreeBufferArray(scip, &cliqueexit);
   SCIPfreeBufferArray(scip, &stacknextcliquevar);
   SCIPfreeBufferArray(scip, &stacknextedge);
   SCIPfreeBufferArray(scip, &dfsstack);

   return SCIP_OKAY;
}

/** initialize candidate lists */
static
SCIP_RETCODE initializeCandsLists(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< structure containing heurdata */
   )
{
   SCIP_VAR** vars;
   int* vbs;
   int nvars;
   int nvbs;
   int v;

   SCIPdebugMsg(scip, "initialize variable bound heuristic (%s)\n", SCIPgetProbName(scip));

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNIntVars(scip) + SCIPgetNBinVars(scip) + SCIPgetNImplVars(scip);
   nvbs = 0;

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->initialized = TRUE;

   if( nvars == 0 )
      return SCIP_OKAY;

   /* allocate memory for the arrays of the heurdata */
   SCIP_CALL( SCIPallocBufferArray(scip, &vbs, 2 * nvars) );

   /* create the topological sorted variable array with respect to the variable bounds */
   SCIP_CALL( topologicalSort(scip, vbs, &nvbs) );

   /* check if the candidate list contains enough candidates */
   if( nvbs > 0 && nvbs >= 0.1 * heurdata->minintfixingrate * nvars )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->vbvars, nvbs) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->vbbounds, nvbs) );

      /* capture variable candidate list */
      for( v = 0; v < nvbs; ++v )
      {
         heurdata->vbvars[v] = vars[getVarIndex(vbs[v])];
         heurdata->vbbounds[v] = getBoundtype(vbs[v]);
         assert(SCIPvarIsIntegral(heurdata->vbvars[v]));

         SCIP_CALL( SCIPcaptureVar(scip, heurdata->vbvars[v]) );
      }

      heurdata->nvbvars = nvbs;
      heurdata->applicable = TRUE;
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &vbs);

   SCIPstatisticMessage("vbvars %.3g (%s)\n",
      (nvbs * 100.0) / nvars, SCIPgetProbName(scip));

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      SCIP_Real minimprove;
      SCIP_Real cutoffbound;

      minimprove = heurdata->minimprove;
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
      {
         cutoffbound = (1-minimprove) * SCIPgetUpperbound(scip) + minimprove * SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound ( scip ) >= 0 )
            cutoffbound = (1 - minimprove) * SCIPgetUpperbound(scip);
         else
            cutoffbound = (1 + minimprove) * SCIPgetUpperbound(scip);
      }
      heurdata->cutoffbound = MIN(upperbound, cutoffbound);
   }
   else
      heurdata->cutoffbound = SCIPinfinity(scip);
   return SCIP_OKAY;
}

/** apply variable bound fixing during probing */
static
SCIP_RETCODE applyVboundsFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< structure containing heurdata */
   SCIP_VAR**            vars,               /**< variables to fix during probing */
   int                   nvbvars,            /**< number of variables in the variable bound graph */
   SCIP_Bool             tighten,            /**< should variables be fixed to cause other fixings? */
   int                   obj,                /**< should the objective be taken into account? */
   SCIP_Bool*            allobj1,            /**< pointer to store whether all variables were fixed according to obj=1 scheme */
   SCIP_Bool*            allobj2,            /**< pointer to store whether all variables were fixed according to obj=2 scheme */
   SCIP_Bool*            backtracked,        /**< was backtracking performed at least once? */
   SCIP_Bool*            infeasible          /**< pointer to store whether propagation detected infeasibility */
   )
{
   SCIP_VAR* lastvar;
   SCIP_VAR* var;
   SCIP_Real lastfixval;
   SCIP_Bool lastfixedlb;
   SCIP_Bool fixtolower;
   SCIP_BOUNDTYPE bound;
   int nbacktracks = 0;
   int v;

   /* loop over variables in topological order */
   for( v = 0; v < nvbvars && !(*infeasible); ++v )
   {
      var = vars[v];
      bound = heurdata->vbbounds[v];

      /*SCIPdebugMsg(scip, "topoorder[%d]: %s(%s) (%s) [%g,%g] (obj=%g)\n", v,
         bound == SCIP_BOUNDTYPE_UPPER ? "ub" : "lb", SCIPvarGetName(var),
         SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS ? "c" : "d",
         SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetObj(var));*/

      /* only check integer or binary variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         continue;

      /* skip variables which are already fixed */
      if( SCIPvarGetLbLocal(var) + 0.5 > SCIPvarGetUbLocal(var) )
         continue;


      /* there are two cases for tighten:
       * 1) tighten == TRUE:  we go through the list of variables and fix variables to force propagation;
       *                      this is be obtained by fixing the variable to the other bound (which means
       *                      that the current bound is changed and so, much propagation is triggered
       *                      since we are starting with the bounds which are most influential).
       * 2) tighten == FALSE: we fix variables to avoid too much propagation in order to avoid reaching
       *                      infeasibility. Therefore, we fix the variable to the current bound, so that
       *                      this bound is not changed and does not propagate. The other bound is changed
       *                      and propagates, but is later in the order, so less influential.
       */
      fixtolower = (tighten == (bound == SCIP_BOUNDTYPE_UPPER));

      /* if we want to take into account the objective function coefficients, we only perform a fixing if the variable
       *  would be fixed to its best bound; otherwise, we just continue
       */
      if( ((SCIPvarGetObj(var) >= 0) != fixtolower) )
      {
         if( obj == 1 )
            continue;
         else
            *allobj1 = FALSE;
      }
      /* if we want to take into account the objective function coefficients but reverted, we only perform a fixing if the variable
       *  would be fixed to its worst bound; otherwise, we just continue
       */
      if( ((SCIPvarGetObj(var) >= 0) == fixtolower) )
      {
         if( obj == 2 )
            continue;
         else
            *allobj2 = FALSE;
      }
      lastvar = var;

      /* fix the variable to its bound */
      if( fixtolower )
      {
         /* we cannot fix to infinite bounds */
         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
            continue;

         /* only open a new probing node if we will not exceed the maximal tree depth */
         if( SCIP_MAXTREEDEPTH > SCIPgetDepth(scip) )
         {
            SCIP_CALL( SCIPnewProbingNode(scip) );
         }

         /* fix variable to lower bound */
         SCIP_CALL( SCIPfixVarProbing(scip, var, SCIPvarGetLbLocal(var)) );
         SCIPdebugMsg(scip, "fixing %d: variable <%s> (obj=%g) to lower bound <%g> (%d pseudo cands)\n",
            v, SCIPvarGetName(var), SCIPvarGetObj(var), SCIPvarGetLbLocal(var), SCIPgetNPseudoBranchCands(scip));
         lastfixedlb = TRUE;
         lastfixval = SCIPvarGetLbLocal(var);
      }
      else
      {
         /* we cannot fix to infinite bounds */
         if( SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
            continue;

         /* only open a new probing node if we will not exceed the maximal tree depth */
         if( SCIP_MAXTREEDEPTH > SCIPgetDepth(scip) )
         {
            SCIP_CALL( SCIPnewProbingNode(scip) );
         }

         /* fix variable to upper bound */
         SCIP_CALL( SCIPfixVarProbing(scip, var, SCIPvarGetUbLocal(var)) );
         SCIPdebugMsg(scip, "fixing %d: variable <%s> (obj=%g) to upper bound <%g> (%d pseudo cands)\n",
            v, SCIPvarGetName(var), SCIPvarGetObj(var), SCIPvarGetUbLocal(var), SCIPgetNPseudoBranchCands(scip));
         lastfixedlb = FALSE;
         lastfixval = SCIPvarGetUbLocal(var);
      }

      /* check if problem is already infeasible */
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, infeasible, NULL) );

      /* probing detected infeasibility: backtrack */
      if( *infeasible )
      {
         assert(lastvar != NULL);

         SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );
         ++nbacktracks;
         *infeasible = FALSE;

         /* increase the lower bound of the variable which caused the infeasibility */
         if( lastfixedlb && lastfixval + 0.5 < SCIPvarGetUbLocal(lastvar) )
         {
            if( lastfixval + 0.5 > SCIPvarGetLbLocal(lastvar) )
            {
               SCIP_CALL( SCIPchgVarLbProbing(scip, lastvar, lastfixval + 1.0) );
            }
         }
         else if( !lastfixedlb && lastfixval - 0.5 > SCIPvarGetLbLocal(lastvar) )
         {
            if( lastfixval - 0.5 < SCIPvarGetUbLocal(lastvar) )
            {
               SCIP_CALL( SCIPchgVarUbProbing(scip, lastvar, lastfixval - 1.0) );
            }
         }
         /* because of the limited number of propagation rounds, it may happen that conflict analysis finds a valid
          * global bound for the last fixed variable that conflicts with applying the reverse bound change after backtracking;
          * in that case, we ran into a deadend and stop
          */
         else
         {
            *infeasible = TRUE;
         }
         lastvar = NULL;

         if( !(*infeasible) )
         {
            /* propagate fixings */
            SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, infeasible, NULL) );

            SCIPdebugMessage("backtrack %d was %sfeasible\n", nbacktracks, (*infeasible ? "in" : ""));
         }

         if( *infeasible )
         {
            SCIPdebugMsg(scip, "probing was infeasible after %d backtracks\n", nbacktracks);

            break;
         }
         else if( nbacktracks > heurdata->maxbacktracks )
         {
            SCIPdebugMsg(scip, "interrupt probing after %d backtracks\n", nbacktracks);
            break;
         }
      }
   }

   *backtracked = (nbacktracks > 0);

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_SOL*             newsol,             /**< working solution */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( subvars != NULL );
   assert( subsol != NULL );

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert( nvars <= SCIPgetNOrigVars(subscip) );

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** copy problem to sub-SCIP, solve it, and add solutions */
static
SCIP_RETCODE setupAndSolveSubscip(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VAR**            vars,               /**< variables of the main SCIP */
   int                   nvars,              /**< number of variables of the main SCIP */
   SCIP_SOL*             sol,                /**< working solution */
   SCIP_Longint          nstallnodes,        /**< stalling node limit for the sub-SCIP */
   SCIP_Real             lowerbound,         /**< lower bound of the main SCIP / current subproblem */
   int*                  nprevars,           /**< pointer to store the number of presolved variables */
   SCIP_Bool*            wasfeas,            /**< pointer to store if a feasible solution was found */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_VAR** subvars;
   SCIP_HASHMAP* varmap;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heurdata != NULL);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), nvars) );

   SCIP_CALL( SCIPcopyConsCompression(scip, subscip, varmap, NULL, "_vbounds", NULL, NULL, 0, FALSE, FALSE, TRUE, NULL) );

   if( heurdata->copycuts )
   {
      /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
      SCIP_CALL( SCIPcopyCuts(scip, subscip, varmap, NULL, FALSE, NULL) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmap);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* set limits for the subproblem */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* forbid call of heuristics and separators solving sub-CIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handlers; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no decutions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
   }

   /* set a cutoff bound */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      SCIP_Real minimprove;
      SCIP_Real cutoffbound;

      minimprove = heurdata->minimprove;
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip, -1.0 * lowerbound) )
      {
         cutoffbound = (1-minimprove) * SCIPgetUpperbound(scip) + minimprove * lowerbound;
      }
      else
      {
         if( SCIPgetUpperbound ( scip ) >= 0 )
            cutoffbound = (1 - minimprove) * SCIPgetUpperbound(scip);
         else
            cutoffbound = (1 + minimprove) * SCIPgetUpperbound(scip);
      }
      heurdata->cutoffbound = MIN(upperbound, cutoffbound);
   }

   if( !SCIPisInfinity(scip, heurdata->cutoffbound) )
   {
      SCIP_CALL( SCIPsetObjlimit(subscip, heurdata->cutoffbound) );
      SCIPdebugMsg(scip, "setting objlimit for subscip to %g\n", heurdata->cutoffbound);
   }

   SCIPdebugMsg(scip, "starting solving vbound-submip at time %g\n", SCIPgetSolvingTime(scip));

   /* solve the subproblem */
   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   SCIP_CALL_ABORT( SCIPpresolve(subscip) );

   SCIPdebugMsg(scip, "vbounds heuristic presolved subproblem at time %g : %d vars, %d cons; fixing value = %g\n",
      SCIPgetSolvingTime(scip), SCIPgetNVars(subscip), SCIPgetNConss(subscip),
      ((nvars - SCIPgetNVars(subscip)) / (SCIP_Real)nvars));

   *nprevars = SCIPgetNVars(subscip);

   /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
    * to ensure that not only the MIP but also the LP relaxation is easy enough
    */
   if( ((nvars - SCIPgetNVars(subscip)) / (SCIP_Real)nvars) >= heurdata->minmipfixingrate )
   {
      SCIP_SOL** subsols;
      SCIP_Bool success;
      int nsubsols;

      SCIPdebugMsg(scip, "solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->maxnodes);

      SCIP_CALL_ABORT( SCIPsolve(subscip) );

      SCIPdebugMsg(scip, "ending solving vbounds-submip at time %g, status = %d\n", SCIPgetSolvingTime(scip), SCIPgetStatus(subscip));

      /* check, whether a solution was found; due to numerics, it might happen that not all solutions are feasible ->
       * try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      *wasfeas = FALSE;

      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, sol, subsols[i], &success) );
         if( !(*wasfeas) )
         {
            SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, wasfeas) );
            if( (*wasfeas) )
               SCIPdebugMsg(scip, "found feasible solution in sub-MIP: %16.9g\n", SCIPgetSolOrigObj(scip, sol));
         }
      }
      if( success )
      {
         *result = SCIP_FOUNDSOL;
      }
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}

/** main procedure of the vbounds heuristic */
static
SCIP_RETCODE applyVbounds(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            vbvars,             /**< variables to fix during probing */
   int                   nvbvars,            /**< number of variables to fix */
   SCIP_Bool             tighten,            /**< should variables be fixed to cause other fixings? */
   int                   obj,                /**< should the objective be taken into account? */
   SCIP_Bool*            skipobj1,           /**< pointer to store whether the run with obj=1 can be skipped, or NULL */
   SCIP_Bool*            skipobj2,           /**< pointer to store whether the run with obj=2 can be skipped, or NULL */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIPstatistic( SCIP_CLOCK* clock; )
   SCIP_VAR** vars;
   SCIP_SOL* sol = NULL;
   SCIP_Longint nstallnodes;
   SCIP_LPSOLSTAT lpstatus;
   SCIP_Real lowerbound;
   SCIP_Bool wasfeas = FALSE;
   SCIP_Bool cutoff;
   SCIP_Bool lperror;
   SCIP_Bool solvelp;
   SCIP_Bool allobj1 = TRUE;
   SCIP_Bool allobj2 = TRUE;
   SCIP_Bool backtracked = TRUE;
   int oldnpscands;
   int npscands;
   int nvars;
   int nprevars;

   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(nvbvars > 0);

   /* initialize default values */
   cutoff = FALSE;

   if( skipobj1 != NULL )
      *skipobj1 = FALSE;
   if( skipobj2 != NULL )
      *skipobj2 = FALSE;

   /* get variable data of original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIPstatistic( nprevars = nvars; )

   if( nvbvars < nvars * heurdata->minintfixingrate )
      return SCIP_OKAY;

   if( *result == SCIP_DIDNOTRUN )
      *result = SCIP_DIDNOTFIND;

   lowerbound = SCIPgetLowerbound(scip);

   oldnpscands = SCIPgetNPseudoBranchCands(scip);

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward variable bounds heuristic if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   SCIPdebugMsg(scip, "apply variable bounds heuristic at node %lld on %d variable bounds, tighten: %d obj: %d\n",
      SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), nvbvars, tighten, obj);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping " HEUR_NAME ": nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIPstatistic( SCIP_CALL( SCIPcreateClock(scip, &clock) ) );
   SCIPstatistic( SCIP_CALL( SCIPstartClock(scip, clock) ) );

   /* check whether the LP should be solved at the current node in the tree to determine whether the heuristic
    * is allowed to solve an LP
    */
   solvelp = SCIPhasCurrentNodeLP(scip);

   if( !SCIPisLPConstructed(scip) && solvelp )
   {
      SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

      /* manually cut off the node if the LP construction detected infeasibility (heuristics cannot return such a result) */
      if( cutoff )
      {
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetCurrentNode(scip)) );
         goto TERMINATE;
      }

      SCIP_CALL( SCIPflushLP(scip) );
   }

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

#ifdef COLLECTSTATISTICS
   SCIPenableVarHistory(scip);
#endif

   /* create temporary solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   /* apply the variable fixings */
   SCIP_CALL( applyVboundsFixings(scip, heurdata, vbvars, nvbvars, tighten, obj, &allobj1, &allobj2, &backtracked, &cutoff) );

   if( skipobj1 != NULL )
      *skipobj1 = allobj1;

   if( skipobj2 != NULL )
      *skipobj2 = allobj2;

   if( cutoff || SCIPisStopped(scip) )
      goto TERMINATE;

   /* check that we had enough fixings */
   npscands = SCIPgetNPseudoBranchCands(scip);

   SCIPdebugMsg(scip, "npscands=%d, oldnpscands=%d, heurdata->minintfixingrate=%g\n", npscands, oldnpscands, heurdata->minintfixingrate);

   /* check fixing rate */
   if( npscands > oldnpscands * (1.0 - heurdata->minintfixingrate) )
   {
      if( heurdata->uselockfixings && npscands <= 2.0 * oldnpscands * (1.0 - heurdata->minintfixingrate) )
      {
         SCIP_Bool allrowsfulfilled = FALSE;

         SCIP_CALL( SCIPapplyLockFixings(scip, NULL, &cutoff, &allrowsfulfilled) );

         if( cutoff || SCIPisStopped(scip) )
         {
            SCIPdebugMsg(scip, "cutoff or timeout in locks fixing\n");
            goto TERMINATE;
         }

         npscands = SCIPgetNPseudoBranchCands(scip);

         SCIPdebugMsg(scip, "after lockfixings: npscands=%d, oldnpscands=%d, allrowsfulfilled=%u, heurdata->minintfixingrate=%g\n",
            npscands, oldnpscands, allrowsfulfilled, heurdata->minintfixingrate);

         if( !allrowsfulfilled && npscands > oldnpscands * (1 - heurdata->minintfixingrate) )
         {
            SCIPdebugMsg(scip, "--> too few fixings\n");

            goto TERMINATE;
         }
      }
      else
      {
         SCIPdebugMsg(scip, "--> too few fixings\n");

         goto TERMINATE;
      }
   }

   assert(!cutoff);

   /*************************** Probing LP Solving ***************************/
   lpstatus = SCIP_LPSOLSTAT_ERROR;
   lperror = FALSE;
   /* solve lp only if the problem is still feasible */
   if( solvelp )
   {
      SCIPdebugMsg(scip, "starting solving vbound-lp at time %g\n", SCIPgetSolvingTime(scip));

      /* solve LP; errors in the LP solver should not kill the overall solving process, if the LP is just needed for a
       * heuristic.  hence in optimized mode, the return code is caught and a warning is printed, only in debug mode,
       * SCIP will stop.
       */
#ifdef NDEBUG
      {
         SCIP_Bool retstat;
         retstat = SCIPsolveProbingLP(scip, -1, &lperror, NULL);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while solving LP in vbound heuristic; LP solve terminated with code <%d>\n",
               retstat);
         }
      }
#else
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
#endif
      SCIPdebugMsg(scip, "ending solving vbound-lp at time %g\n", SCIPgetSolvingTime(scip));

      lpstatus = SCIPgetLPSolstat(scip);

      SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
      SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, lpstatus);
   }

   /* check if this is a feasible solution */
   if( lpstatus == SCIP_LPSOLSTAT_OPTIMAL && !lperror )
   {
      SCIP_Bool stored;
      SCIP_Bool success;

      lowerbound = SCIPgetLPObjval(scip);

      /* copy the current LP solution to the working solution */
      SCIP_CALL( SCIPlinkLPSol(scip, sol) );

      SCIP_CALL( SCIProundSol(scip, sol, &success) );

      if( success )
      {
         SCIPdebugMsg(scip, "vbound heuristic found roundable primal solution: obj=%g\n",
            SCIPgetSolOrigObj(scip, sol));

         /* check solution for feasibility, and add it to solution store if possible.
          * Neither integrality nor feasibility of LP rows have to be checked, because they
          * are guaranteed by the heuristic at this stage.
          */
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPtrySol(scip, sol, TRUE, TRUE, TRUE, TRUE, TRUE, &stored) );
#else
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, FALSE, FALSE, &stored) );
#endif

#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, &wasfeas) );
         assert(wasfeas);
         SCIPdebugMsg(scip, "found feasible solution by LP rounding: %16.9g\n", SCIPgetSolOrigObj(scip, sol));
#endif

         if( stored )
         {
            *result = SCIP_FOUNDSOL;
         }

         /* we found a solution, so we are done */
         goto TERMINATE;
      }
   }
   /*************************** END Probing LP Solving ***************************/

   /*************************** Start Subscip Solving ***************************/
   /* if no solution has been found --> fix all other variables by subscip if necessary */
   if( !lperror && lpstatus != SCIP_LPSOLSTAT_INFEASIBLE && lpstatus != SCIP_LPSOLSTAT_OBJLIMIT )
   {
      SCIP* subscip;
      SCIP_RETCODE retcode;
      SCIP_Bool valid;

      /* check whether there is enough time and memory left */
      SCIP_CALL( SCIPcheckCopyLimits(scip, &valid) );

      if( !valid )
         goto TERMINATE;

      /* create subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      retcode = setupAndSolveSubscip(scip, subscip, heurdata, vars, nvars, sol, nstallnodes, lowerbound,
         &nprevars, &wasfeas, result);

      SCIP_CALL( SCIPfree(&subscip) );

      SCIP_CALL( retcode );
   }

   /*************************** End Subscip Solving ***************************/

 TERMINATE:
#ifdef SCIP_STATISTIC
   SCIP_CALL( SCIPstopClock(scip, clock) );
   SCIPstatisticMessage("vbound: tighten=%u obj=%d nvars=%d presolnvars=%d ratio=%.2f infeas=%u found=%d time=%.4f\n",
      tighten, obj, nvars, nprevars, (nvars - nprevars) / (SCIP_Real)nvars, cutoff,
      wasfeas ? 1 : 0, SCIPclockGetTime(clock) );
#endif

   SCIPstatistic( SCIP_CALL( SCIPfreeClock(scip, &clock) ) );

   /* free solution */
   if( sol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &sol) );
   }

   /* exit probing mode */
   if( SCIPinProbing(scip) )
   {
      SCIP_CALL( SCIPendProbing(scip) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyVbounds)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of heuristic */
   SCIP_CALL( SCIPincludeHeurVbounds(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeVbounds)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolVbounds)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int v;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* release all variables */
   for( v = 0; v < heurdata->nvbvars; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &heurdata->vbvars[v]) );
   }

   /* free varbounds array */
   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->vbbounds, heurdata->nvbvars);
   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->vbvars, heurdata->nvbvars);

   /* reset heuristic data structure */
   heurdataReset(heurdata);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecVbounds)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool skipobj1;
   SCIP_Bool skipobj2;
#ifdef NOCONFLICT
   SCIP_Bool enabledconflicts;
#endif

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNPseudoBranchCands(scip) == 0 )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( !heurdata->initialized )
   {
      SCIP_CALL( initializeCandsLists(scip, heurdata) );
   }

   if( !heurdata->applicable )
      return SCIP_OKAY;

#ifdef NOCONFLICT
   /* disable conflict analysis */
   SCIP_CALL( SCIPgetBoolParam(scip, "conflict/enable", &enabledconflicts) );

   if( !SCIPisParamFixed(scip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", FALSE) );
   }
#endif

   /* try variable bounds */
   skipobj1 = FALSE;
   skipobj2 = FALSE;
   if( (heurdata->feasvariant & VBOUNDVARIANT_NOOBJ) != 0 )
   {
      SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, FALSE, 0,
            &skipobj1, &skipobj2, result) );
   }
   if( !skipobj1 && (heurdata->feasvariant & VBOUNDVARIANT_BESTBOUND) != 0)
   {
      SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, FALSE, 1, NULL, NULL, result) );
   }
   if( !skipobj2 && (heurdata->feasvariant & VBOUNDVARIANT_WORSTBOUND) != 0)
   {
      SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, FALSE, 2, NULL, NULL, result) );
   }

   skipobj1 = FALSE;
   skipobj2 = FALSE;
   if( (heurdata->tightenvariant & VBOUNDVARIANT_NOOBJ) != 0 )
   {
      SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, TRUE, 0,
            &skipobj1, &skipobj2, result) );
   }
   if( !skipobj1 && (heurdata->tightenvariant & VBOUNDVARIANT_BESTBOUND) != 0)
   {
      SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, TRUE, 1, NULL, NULL, result) );
   }
   if( !skipobj2 && (heurdata->tightenvariant & VBOUNDVARIANT_WORSTBOUND) != 0)
   {
      SCIP_CALL( applyVbounds(scip, heur, heurdata, heurdata->vbvars, heurdata->nvbvars, TRUE, 2, NULL, NULL, result) );
   }

#ifdef NOCONFLICT
   /* reset the conflict analysis */
   if( !SCIPisParamFixed(scip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", enabledconflicts) );
   }
#endif

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the vbounds primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create vbounds primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   heurdataReset(heurdata);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecVbounds, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyVbounds) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeVbounds) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolVbounds) );

   /* add variable bounds primal heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minintfixingrate",
         "minimum percentage of integer variables that have to be fixed",
         &heurdata->minintfixingrate, FALSE, DEFAULT_MININTFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minmipfixingrate",
         "minimum percentage of variables that have to be fixed within sub-SCIP (integer and continuous)",
         &heurdata->minmipfixingrate, FALSE, DEFAULT_MINMIPFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which " HEUR_NAME " heuristic should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxproprounds",
         "maximum number of propagation rounds during probing (-1 infinity)",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX/4, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselockfixings",
         "should more variables be fixed based on variable locks if the fixing rate was not reached?",
         &heurdata->uselockfixings, TRUE, DEFAULT_USELOCKFIXINGS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxbacktracks",
         "maximum number of backtracks during the fixing process",
         &heurdata->maxbacktracks, TRUE, DEFAULT_MAXBACKTRACKS, -1, INT_MAX/4, NULL, NULL) );

      SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/feasvariant",
         "which variants of the vbounds heuristic that try to stay feasible should be called? (0: off, 1: w/o looking at obj, 2: only fix to best bound, 4: only fix to worst bound",
         &heurdata->feasvariant, TRUE, DEFAULT_FEASVARIANT, 0, 7, NULL, NULL) );

      SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/tightenvariant",
         "which tightening variants of the vbounds heuristic should be called? (0: off, 1: w/o looking at obj, 2: only fix to best bound, 4: only fix to worst bound",
         &heurdata->tightenvariant, TRUE, DEFAULT_TIGHTENVARIANT, 0, 7, NULL, NULL) );

   return SCIP_OKAY;
}
