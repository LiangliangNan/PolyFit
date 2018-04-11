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

/**@file   heur_clique.c
 * @brief  LNS heuristic using a clique partition to restrict the search neighborhood
 * @brief  clique primal heuristic
 * @author Stefan Heinz
 * @author Michael Winkler
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
#include "scip/pub_implics.h"
#include "scip/heur_clique.h"
#include "scip/heur_locks.h"
#include "scip/cons_logicor.h"
#include "scip/pub_misc.h"


#define HEUR_NAME             "clique"
#define HEUR_DESC             "LNS heuristic using a clique partition to restrict the search neighborhood"
#define HEUR_DISPCHAR         'Q'
#define HEUR_PRIORITY         5000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE                       /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL                     /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MININTFIXINGRATE 0.65                    /**< minimum percentage of integer variables that have to be fixed */
#define DEFAULT_MINMIPFIXINGRATE 0.65                    /**< minimum percentage of variables that have to be fixed within sub-SCIP
                                                          *   (integer and continuous) */
#define DEFAULT_MINIMPROVE    0.01                       /**< factor by which clique heuristic should at least improve the
                                                          *   incumbent
                                                          */
#define DEFAULT_MINNODES      500LL                      /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS      500LL                      /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT     0.1                        /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_MAXPROPROUNDS 2                          /**< maximum number of propagation rounds during probing */
#define DEFAULT_MAXBACKTRACKS 10                         /**< maximum number of backtracks during the fixing process */
#define DEFAULT_COPYCUTS      TRUE                       /**< should all active cuts from the cutpool of the
                                                          *   original scip be copied to constraints of the subscip
                                                          */
#define DEFAULT_USELOCKFIXINGS FALSE                     /**< should more variables be fixed based on variable locks if
                                                          *   the fixing rate was not reached?
                                                          */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Longint          usednodes;          /**< nodes already used by clique heuristic in earlier calls */
   SCIP_Real             minintfixingrate;   /**< minimum percentage of integer variables that have to be fixed */
   SCIP_Real             minmipfixingrate;   /**< minimum percentage of variables that have to be fixed within sub-SCIP
                                              *   (integer and continuous) */
   SCIP_Real             minimprove;         /**< factor by which clique heuristic should at least improve the incumbent */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   int                   maxproprounds;      /**< maximum number of propagation rounds during probing */
   int                   maxbacktracks;      /**< maximum number of backtracks during the fixing process */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   subproblem?
                                              */
   SCIP_Bool             uselockfixings;     /**< should more variables be fixed based on variable locks if
                                              *   the fixing rate was not reached?
                                              */
};

/*
 * Local methods
 */

/** comparison method for sorting cliques by their size */
static
SCIP_DECL_SORTINDCOMP(compCliquesSize)
{
   int* cliquesizes = (int*)dataptr;

   return cliquesizes[ind2] - cliquesizes[ind1];
}

static
int getCliqueUnfixedVars(
   SCIP_CLIQUE*          clique
   )
{
   SCIP_VAR** cliquevars;
   SCIP_VAR* var;
   int ncliquevars;
   int nunfixed = 0;
   int v;

   ncliquevars = SCIPcliqueGetNVars(clique);
   cliquevars = SCIPcliqueGetVars(clique);

   for( v = 0; v < ncliquevars; ++v )
   {
      var = cliquevars[v];

      /* is variable unfixed? */
      if( SCIPvarGetUbLocal(var) > SCIPvarGetLbLocal(var) + 0.5 )
         ++nunfixed;
   }

   return nunfixed;
}

/** apply clique fixing using probing */
static
SCIP_RETCODE applyCliqueFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< structure containing heurdata */
   SCIP_VAR**            onefixvars,         /**< array to store all variables which are fixed to one in the cliques */
   SCIP_Shortbool*       onefixvals,         /**< array to store the values of all variables fixed to one in the cliques */
   int*                  nonefixvars,        /**< pointer to store the number of variables fixed to one */
   SCIP_Bool*            cutoff              /**< pointer to store whether the propagation stopped with infeasibility */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_CLIQUE* clique;
   SCIP_VAR** cliquevars;
   SCIP_VAR* var;
   SCIP_Bool* cliquevals;
   SCIP_Bool* propagated;
   int* cliquesizes;
   int* permutation;
   SCIP_Real bestobj;
   SCIP_Real obj;
   SCIP_Bool alreadyone;
   SCIP_Bool newnode;
   int probingdepthofonefix;
   int ncliquevars;
   int ncliques;
   int bestpos;
   int firstclique;
   int bestclique;
   int cliquesize;
   int bestcliquesize;
   int nbacktracks = 0;
   int v = 0;
   int c;
   int i;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(onefixvars != NULL);
   assert(nonefixvars != NULL);
   assert(cutoff != NULL);

   cliques = SCIPgetCliques(scip);
   ncliques = SCIPgetNCliques(scip);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquesizes, ncliques) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permutation, ncliques) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &propagated, ncliques) );

   for( c = ncliques - 1; c >= 0; --c )
   {
      cliquesizes[c] = SCIPcliqueGetNVars(cliques[c]);
   }

   SCIPsort(permutation, compCliquesSize, (void*)cliquesizes, ncliques);

#ifndef NDEBUG
   for( c = ncliques - 1; c >= 1; --c )
   {
      assert(cliquesizes[permutation[c]] <= cliquesizes[permutation[c-1]]);
   }
#endif

   *cutoff = FALSE;
   probingdepthofonefix = 0;
   firstclique = 0;

   SCIP_CALL( SCIPnewProbingNode(scip) );

   /* @todo maybe try to fix more than one variable to one in each probing node, to gain faster results */
   for( c = 0; c < ncliques; ++c )
   {
      bestpos = -1;
      bestobj = SCIPinfinity(scip);
      alreadyone = FALSE;
      newnode = FALSE;

      bestclique = firstclique;

      if( bestclique >= ncliques )
         break;

      bestcliquesize = getCliqueUnfixedVars(cliques[permutation[bestclique]]);
      assert(!propagated[permutation[bestclique]]);

      for( i = firstclique + 1; i < ncliques; ++i)
      {
         if( cliquesizes[permutation[i]] < bestcliquesize )
            break;

         if( propagated[permutation[i]] )
            continue;

         cliquesize = getCliqueUnfixedVars(cliques[permutation[i]]);

         if( cliquesize > bestcliquesize )
         {
            bestclique = i;
            bestcliquesize = cliquesize;
         }
         else if( cliquesize == 0 )
         {
            propagated[permutation[i]] = TRUE;
         }
      }
      clique = cliques[permutation[bestclique]];
      propagated[permutation[bestclique]] = TRUE;

      while( firstclique < ncliques && propagated[permutation[firstclique]] )
         ++firstclique;

      ncliquevars = SCIPcliqueGetNVars(clique);
      cliquevars = SCIPcliqueGetVars(clique);
      cliquevals = SCIPcliqueGetValues(clique);

      for( v = 0; v < ncliquevars; ++v )
      {
         var = cliquevars[v];

         /* variable is already fixed */
         if( SCIPvarGetUbLocal(var) < SCIPvarGetLbLocal(var) + 0.5 )
         {
            SCIPdebugMessage("<%s> is already fixed to %g\n", SCIPvarGetName(var), SCIPvarGetUbLocal(var));

            /* clique variable is fixed to 1 */
            if( cliquevals[v] == (SCIPvarGetLbLocal(var) > 0.5) )
            {
               assert(!alreadyone);
               alreadyone = TRUE;
               break;
            }
            continue;
         }

         obj = cliquevals[v] ? SCIPvarGetObj(var) : -SCIPvarGetObj(var);

         /* @todo use a tiebreaker (locks?) */
         if( obj < bestobj )
         {
            /* variable is not the best one in the clique anymore, fix it to 0 */
            if( bestpos >= 0 )
            {
               if( cliquevals[bestpos] )
               {
                  SCIP_CALL( SCIPfixVarProbing(scip, cliquevars[bestpos], 0.0) );
               }
               else
               {
                  SCIP_CALL( SCIPfixVarProbing(scip, cliquevars[bestpos], 1.0) );
               }
               SCIPdebugMessage("fixed <%s> to %g\n", SCIPvarGetName(cliquevars[bestpos]), SCIPvarGetUbLocal(cliquevars[bestpos]));
               newnode = TRUE;
            }

            bestobj = obj;
            bestpos = v;
         }
         /* variable is not the best one in the clique, fix it to 0 */
         else
         {
            assert(bestpos >= 0);

            if( cliquevals[v] )
            {
               SCIP_CALL( SCIPfixVarProbing(scip, var, 0.0) );
            }
            else
            {
               SCIP_CALL( SCIPfixVarProbing(scip, var, 1.0) );
            }
            SCIPdebugMessage("fixed <%s> to %g\n", SCIPvarGetName(var), SCIPvarGetUbLocal(var));
            newnode = TRUE;
         }
      }
      /* we found a variable in the clique which is already fixed to 1 */
      if( alreadyone )
      {
         /* fix (so far) best candidate to 0 */
         if( bestpos >= 0 )
         {
            if( cliquevals[bestpos] )
            {
               SCIP_CALL( SCIPfixVarProbing(scip, cliquevars[bestpos], 0.0) );
            }
            else
            {
               SCIP_CALL( SCIPfixVarProbing(scip, cliquevars[bestpos], 1.0) );
            }
            SCIPdebugMessage("fixed <%s> to %g\n", SCIPvarGetName(cliquevars[bestpos]), SCIPvarGetUbLocal(cliquevars[bestpos]));
            newnode = TRUE;
         }

         /* fix all variables not yet processed to 0 */
         for( ; v < ncliquevars; ++v )
         {
            var = cliquevars[v];

            if( SCIPvarGetUbLocal(var) < SCIPvarGetLbLocal(var) + 0.5 )
               continue;

            if( cliquevals[v] )
            {
               SCIP_CALL( SCIPfixVarProbing(scip, var, 0.0) );
            }
            else
            {
               SCIP_CALL( SCIPfixVarProbing(scip, var, 1.0) );
            }
            SCIPdebugMessage("fixed <%s> to %g\n", SCIPvarGetName(var), SCIPvarGetUbLocal(var));
            newnode = TRUE;
         }
      }
      /* fix the best variable to 1 */
      else if( bestpos >= 0 )
      {
         assert(bestpos <= ncliquevars);

         probingdepthofonefix = SCIPgetProbingDepth(scip);
         onefixvars[(*nonefixvars)] = cliquevars[bestpos];

         /* @todo should we even fix the best candidate to 1? */
         if( cliquevals[bestpos] )
         {
            SCIP_CALL( SCIPfixVarProbing(scip, cliquevars[bestpos], 1.0) );
            onefixvals[(*nonefixvars)] = 1;
         }
         else
         {
            SCIP_CALL( SCIPfixVarProbing(scip, cliquevars[bestpos], 0.0) );
            onefixvals[(*nonefixvars)] = 0;
         }
         SCIPdebugMessage("fixed <%s> to %g*\n", SCIPvarGetName(cliquevars[bestpos]), SCIPvarGetUbLocal(cliquevars[bestpos]));
         ++(*nonefixvars);
         newnode = TRUE;
      }

      if( newnode )
      {
         /* propagate fixings */
         SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, cutoff, NULL) );

         SCIPdebugMessage("propagate fixings of clique %d: cutoff=%u\n", c, *cutoff);

         if( SCIPisStopped(scip) )
            break;

         /* stop if we reached the depth limit */
         if( SCIP_MAXTREEDEPTH <= SCIPgetDepth(scip) )
            break;

         /* probing detected infeasibility: backtrack */
         if( *cutoff )
         {
            if( *nonefixvars > 0 )
            {
               if( probingdepthofonefix > 0 )
               {

                  SCIP_CALL( SCIPbacktrackProbing(scip, probingdepthofonefix - 1) );
                  probingdepthofonefix = 0;
                  ++nbacktracks;

                  /* because of the limited number of propagation rounds, it may happen that conflict analysis finds a
                   * valid global fixing for the last fixed variable that conflicts with applying the reverse fixing
                   * after backtracking; in that case, we ran into a deadend and stop
                   */
                  if( SCIPvarGetLbLocal(onefixvars[*nonefixvars - 1]) < 1.5 - onefixvals[*nonefixvars - 1]
                     && SCIPvarGetUbLocal(onefixvars[*nonefixvars - 1]) > 0.5 - onefixvals[*nonefixvars - 1] )
                  {
                     /* fix the last variable, which was fixed to 1 and led to the cutoff, to 0 */
                     SCIP_CALL( SCIPfixVarProbing(scip, onefixvars[*nonefixvars - 1], 1.0 - onefixvals[*nonefixvars - 1]) );
                     --(*nonefixvars);

                     /* propagate fixings */
                     SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, cutoff, NULL) );

                     SCIPdebugMessage("backtrack %d was %sfeasible\n", nbacktracks, (*cutoff ? "in" : ""));
                  }
#ifndef NDEBUG
                  else
                     assert(*cutoff == TRUE);
#endif
               }
               if( *cutoff )
               {
#ifndef NOCONFLICT
                  SCIP_CONS* conflictcons;
                  char consname[SCIP_MAXSTRLEN];

                  /* create own conflict */
                  (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "conf%" SCIP_LONGINT_FORMAT "", SCIPgetNNodes(scip));

                  /* get variables for the conflict */
                  for( i = 0; i < *nonefixvars; ++i )
                  {
                     /* if the variable was fixed to 1 by the heuristic, get its negated variable */
                     if( onefixvals[i] )
                     {
                        SCIP_CALL( SCIPgetNegatedVar(scip, onefixvars[i], &onefixvars[i]) );
                     }
                  }

                  SCIPdebugMsg(scip, "probing was infeasible after %d backtracks\n", nbacktracks);

                  /* create conflict constraint */
                  SCIP_CALL( SCIPcreateConsLogicor(scip, &conflictcons, consname, *nonefixvars, onefixvars,
                        FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );
                  SCIP_CALL( SCIPaddConsNode(scip, SCIPgetFocusNode(scip), conflictcons, NULL) );
                  SCIPdebugPrintCons(scip, conflictcons, NULL);
                  SCIP_CALL( SCIPreleaseCons(scip, &conflictcons) );
#endif
                  break;
               }
               else if( nbacktracks > heurdata->maxbacktracks )
               {
                  SCIPdebugMsg(scip, "interrupt probing after %d backtracks\n", nbacktracks);
                  break;
               }
            }
            /* we had a cutoff without a single one-fixing, so the current problem seems to be infeasible already */
            else
               break;
         }

         SCIP_CALL( SCIPnewProbingNode(scip) );
      }
   }
   assert((*nonefixvars > 0) || probingdepthofonefix == 0 );

   SCIPfreeBufferArray(scip, &propagated);
   SCIPfreeBufferArray(scip, &permutation);
   SCIPfreeBufferArray(scip, &cliquesizes);

   SCIPdebugMsg(scip, "fixed %d of %d variables in probing\n", v, SCIPgetNBinVars(scip));
   SCIPdebugMsg(scip, "applied %d of %d cliques in probing\n", c, ncliques);
   SCIPdebugMsg(scip, "probing was %sfeasible\n", (*cutoff) ? "in" : "");

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_SOL*             newsol,             /**< working solution */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   int nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);
   assert(success != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyClique)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurClique(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeClique)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitClique)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* reset heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->usednodes = 0;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecClique)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   SCIP_Real lowerbound;
   int nvars;
   int nbinvars;
   int oldnpscands;
   int npscands;
   int i;
   SCIP_Bool cutoff;
   SCIP_Bool lperror;

   SCIP_VAR** onefixvars;
   SCIP_Shortbool* onefixvals;
   int nonefixvars;
   SCIP_Bool enabledconflicts;
   SCIP_LPSOLSTAT lpstatus;
   SCIP_CONS* conflictcons;
   SCIP_Bool solvelp;
   char consname[SCIP_MAXSTRLEN];

   SCIP_Longint nstallnodes;

   SCIP_SOL* sol;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   nbinvars = SCIPgetNBinVars(scip);

   if( nbinvars < 2 )
      return SCIP_OKAY;

   /* check for necessary information to apply this heuristic */
   if( SCIPgetNCliques(scip) == 0 )
      return SCIP_OKAY;

   lowerbound = SCIPgetLowerbound(scip);

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward clique heuristic if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping " HEUR_NAME ": nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   oldnpscands = SCIPgetNPseudoBranchCands(scip);
   onefixvars = NULL;
   onefixvals = NULL;
   sol = NULL;

   /* disable conflict analysis, because we can it better than SCIP itself, cause we have more information */
   SCIP_CALL( SCIPgetBoolParam(scip, "conflict/enable", &enabledconflicts) );

   if( !SCIPisParamFixed(scip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", FALSE) );
   }

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

   *result = SCIP_DIDNOTFIND;

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

#ifdef COLLECTSTATISTICS
   SCIPenableVarHistory(scip);
#endif

   /* create a solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   /* allocate memory for all variables which will be fixed to one during probing */
   SCIP_CALL(SCIPallocBufferArray(scip, &onefixvars, nbinvars) );
   SCIP_CALL(SCIPallocBufferArray(scip, &onefixvals, nbinvars) );
   nonefixvars = 0;

   /* apply fixings due to clique information */
   SCIP_CALL( applyCliqueFixings(scip, heurdata, onefixvars, onefixvals, &nonefixvars, &cutoff) );

   if( cutoff || SCIPisStopped(scip) )
      goto TERMINATE;

   /* check that we had enough fixings */
   npscands = SCIPgetNPseudoBranchCands(scip);

   SCIPdebugMsg(scip, "npscands=%d, oldnpscands=%d, heurdata->minintfixingrate=%g\n", npscands, oldnpscands, heurdata->minintfixingrate);

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

   /*************************** Probing LP Solving ***************************/

   lpstatus = SCIP_LPSOLSTAT_ERROR;
   lperror = FALSE;

   /* solve lp only if the problem is still feasible */
   if( solvelp )
   {
      SCIPdebugMsg(scip, "starting solving clique-lp at time %g\n", SCIPgetSolvingTime(scip));

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
            SCIPwarningMessage(scip, "Error while solving LP in clique heuristic; LP solve terminated with code <%d>\n",
               retstat);
         }
      }
#else
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
#endif
      SCIPdebugMsg(scip, "ending solving clique-lp at time %g\n", SCIPgetSolvingTime(scip));

      lpstatus = SCIPgetLPSolstat(scip);

      SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
      SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, lpstatus);
   }

   /* check if this is a feasible solution */
   if( lpstatus == SCIP_LPSOLSTAT_OPTIMAL && !lperror )
   {
      SCIP_Bool stored;
      SCIP_Bool success;

      assert(!cutoff);

      lowerbound = SCIPgetLPObjval(scip);

      /* copy the current LP solution to the working solution */
      SCIP_CALL( SCIPlinkLPSol(scip, sol) );

      SCIP_CALL( SCIProundSol(scip, sol, &success) );

      if( success )
      {
         SCIPdebugMsg(scip, "clique heuristic found roundable primal solution: obj=%g\n",
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

         if( stored )
         {
            SCIPdebugMsg(scip, "found feasible solution:\n");
            SCIPdebug( SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) ) );
            *result = SCIP_FOUNDSOL;
         }

         /* we found a solution, so we are done */
         goto TERMINATE;
      }
   }
   /*************************** END Probing LP Solving ***************************/


   /*************************** Create Conflict ***************************/
   if( SCIPallColsInLP(scip) && (lpstatus == SCIP_LPSOLSTAT_INFEASIBLE || lpstatus == SCIP_LPSOLSTAT_OBJLIMIT) )
   {
#ifndef NOCONFLICT
      /* create own conflict */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "conf%" SCIP_LONGINT_FORMAT "", SCIPgetNNodes(scip));

      /* get variables for the conflict */
      for( i = 0; i < nonefixvars; ++i )
      {
         /* if the variable was fixed to 1 by the heuristic, get its negated variable */
         if( onefixvals[i] )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, onefixvars[i], &onefixvars[i]) );
         }
      }

      /* create conflict constraint */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &conflictcons, consname, nonefixvars, onefixvars,
            FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddConsNode(scip, SCIPgetFocusNode(scip), conflictcons, NULL) );
      SCIPdebugPrintCons(scip, conflictcons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &conflictcons) );
#endif
      goto TERMINATE;
   }
   /*************************** End Conflict ***************************/

   /*************************** Start Subscip Solving ***************************/
   /* no solution has been found yet and the subproblem is still feasible --> fix all other variables by subscip if
    * necessary
    */
   if( !lperror )
   {
      SCIP* subscip;
      SCIP_VAR** subvars;
      SCIP_HASHMAP* varmap;
      SCIP_Bool valid;

      /* check whether there is enough time and memory left */
      SCIP_CALL( SCIPcheckCopyLimits(scip, &valid) );

      if( !valid )
         goto TERMINATE;

      /* get all variables */
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

      /* create subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      /* allocate temporary memory for subscip variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), nvars) );

      SCIP_CALL( SCIPcopyConsCompression(scip, subscip, varmap, NULL, "_clique", NULL, NULL, 0, FALSE, FALSE, TRUE, &valid) );

      if( heurdata->copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmap, NULL, FALSE, NULL) );
      }

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

      /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
       * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
       * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
       * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
       * made for the original SCIP
       */
      if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
      }

      /* if there is already a solution, add an objective cutoff */
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
         cutoffbound = MIN(upperbound, cutoffbound);
         SCIP_CALL( SCIPsetObjlimit(subscip, cutoffbound) );
         SCIPdebugMsg(scip, "setting objlimit for subscip to %g\n", cutoffbound);
      }

      SCIPdebugMsg(scip, "starting solving clique-submip at time %g\n", SCIPgetSolvingTime(scip));

      /* solve the subproblem */
      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
      SCIP_CALL_ABORT( SCIPpresolve(subscip) );

      SCIPdebugMsg(scip, "clique heuristic presolved subproblem at time %g : %d vars, %d cons; fixing value = %g\n", SCIPgetSolvingTime(scip), SCIPgetNVars(subscip), SCIPgetNConss(subscip), ((nvars - SCIPgetNVars(subscip)) / (SCIP_Real)nvars));

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

         SCIPdebugMsg(scip, "ending solving clique-submip at time %g, status = %d\n", SCIPgetSolvingTime(scip), SCIPgetStatus(subscip));

         /* check, whether a solution was found; due to numerics, it might happen that not all solutions are feasible ->
          * try all solutions until one was accepted
          */
         nsubsols = SCIPgetNSols(subscip);
         subsols = SCIPgetSols(subscip);
         success = FALSE;

         for( i = 0; i < nsubsols && !success; ++i )
         {
            SCIP_CALL( createNewSol(scip, subscip, subvars, sol, subsols[i], &success) );
         }
         if( success )
            *result = SCIP_FOUNDSOL;
#ifndef NOCONFLICT
         /* if subscip was infeasible, add a conflict */
         if( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE )
         {
            /* create own conflict */
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "conf%" SCIP_LONGINT_FORMAT "", SCIPgetNNodes(scip));

            /* get variables for the conflict */
            for( i = 0; i < nonefixvars; ++i )
            {
               /* if the variable was fixed to 1 by the heuristic, get its negated variable */
               if( onefixvals[i] )
               {
                  SCIP_CALL( SCIPgetNegatedVar(scip, onefixvars[i], &onefixvars[i]) );
               }
            }

            /* create conflict constraint */
            SCIP_CALL( SCIPcreateConsLogicor(scip, &conflictcons, consname, nonefixvars, onefixvars,
                  FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddConsNode(scip, SCIPgetFocusNode(scip), conflictcons, NULL) );
            SCIPdebugPrintCons(scip, conflictcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &conflictcons) );
         }
#endif
      }

#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      /* free subproblem */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
   }

   /*************************** End Subscip Solving ***************************/

 TERMINATE:

   /* reset the conflict analysis */
   if( !SCIPisParamFixed(scip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", enabledconflicts) );
   }

   /* free conflict variables */
   SCIPfreeBufferArrayNull(scip, &onefixvars);
   SCIPfreeBufferArrayNull(scip, &onefixvals);

   /* freeing solution */
   if( sol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &sol) );
   }

   /* end probing */
   if( SCIPinProbing(scip) )
   {
      SCIP_CALL( SCIPendProbing(scip) );
   }

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the clique primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurClique(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create clique primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecClique, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyClique) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeClique) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitClique) );

   /* add clique primal heuristic parameters */

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minintfixingrate",
         "minimum percentage of integer variables that have to be fixable",
         &heurdata->minintfixingrate, FALSE, DEFAULT_MININTFIXINGRATE, 0.0, 1.0, NULL, NULL) );

      SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minmipfixingrate",
         "minimum percentage of fixed variables in the sub-MIP",
         &heurdata->minmipfixingrate, FALSE, DEFAULT_MINMIPFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

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

   return SCIP_OKAY;
}
