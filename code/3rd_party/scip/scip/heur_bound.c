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

/**@file   heur_bound.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  heuristic which fixes all integer variables to a bound (lower/upper) and solves the remaining LP
 * @author Gerald Gamrath
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/heur_bound.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include <string.h>

#ifdef SCIP_STATISTIC
#include "scip/clock.h"
#endif

#define HEUR_NAME             "bound"
#define HEUR_DESC             "heuristic which fixes all integer variables to a bound and solves the remaining LP"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_PROP
#define HEUR_PRIORITY         -1107000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE     /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_ONLYWITHOUTSOL   TRUE   /**< Should heuristic only be executed if no primal solution was found, yet? */
#define DEFAULT_MAXPROPROUNDS    0      /* maximum number of propagation rounds during probing */
#define DEFAULT_BOUND           'l'     /**< to which bound should integer variables be fixed? */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Bool             onlywithoutsol;     /**< Should heuristic only be executed if no primal solution was found, yet? */
   int                   maxproprounds;      /**< maximum number of propagation rounds during probing */
   char                  bound;              /**< to which bound should integer variables be fixed? */
};

/*
 * Local methods
 */

/** main procedure of the bound heuristic */
static
SCIP_RETCODE applyBoundHeur(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_Bool             lower,              /**< should integer variables be fixed to their lower bound? */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Bool infeasible = FALSE;
   int maxproprounds;
   int nbinvars;
   int nintvars;
   int nvars;
   int v;

   /* get variable data of original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   maxproprounds = heurdata->maxproprounds;
   if( maxproprounds == -2 )
      maxproprounds = 0;

   /* only look at binary and integer variables */
   nvars = nbinvars + nintvars;

   /* stop if we would have infinite fixings */
   if( lower )
   {
      for( v = 0; v < nvars; ++v )
      {
         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(vars[v])) )
            return SCIP_OKAY;
      }
   }
   else
   {
      for( v = 0; v < nvars; ++v )
      {
         if( SCIPisInfinity(scip, SCIPvarGetUbLocal(vars[v])) )
            return SCIP_OKAY;
      }
   }

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];

      assert(SCIPvarGetType(var) < SCIP_VARTYPE_IMPLINT);

      /* skip variables which are already fixed */
      if( SCIPvarGetLbLocal(var) + 0.5 > SCIPvarGetUbLocal(var) )
         continue;

      /* fix variable to lower bound */
      if( lower )
      {
         SCIP_CALL( SCIPfixVarProbing(scip, var, SCIPvarGetLbLocal(var)) );
         SCIPdebugMsg(scip, "fixing %d: variable <%s> to lower bound <%g> (%d pseudo cands)\n",
            v, SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPgetNPseudoBranchCands(scip));
      }
      /* fix variable to upper bound */
      else
      {
         SCIP_CALL( SCIPfixVarProbing(scip, var, SCIPvarGetUbLocal(var)) );
         SCIPdebugMsg(scip, "fixing %d: variable <%s> to upper bound <%g> (%d pseudo cands)\n",
            v, SCIPvarGetName(var), SCIPvarGetUbLocal(var), SCIPgetNPseudoBranchCands(scip));
      }

      /* propagate fixings */
      if( heurdata->maxproprounds != 0 )
      {
         SCIP_CALL( SCIPpropagateProbing(scip, maxproprounds, &infeasible, NULL) );
      }

      /* try to repair probing */
      if( infeasible )
      {
#if 0
         SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );

         /* fix the last variable, which was fixed the reverse bound */
         SCIP_CALL( SCIPfixVarProbing(scip, var, SCIPvarGetUbLocal(var)) );

         /* propagate fixings */
         SCIP_CALL( SCIPpropagateProbing(scip, maxproprounds, &infeasible, NULL) );

         SCIPdebugMsg(scip, "backtracking ended with %sfeasible problem\n", (infeasible ? "in" : ""));

         if( infeasible )
#endif
            break;
      }
   }

   SCIPdebugMsg(scip, "probing ended with %sfeasible problem\n", infeasible ? "in" : "");

   /*************************** Probing LP Solving ***************************/

   /* solve lp only if the problem is still feasible */
   if( !infeasible )
   {
      char strbuf[SCIP_MAXSTRLEN];
      SCIP_LPSOLSTAT lpstatus;
      SCIP_Bool lperror;
      int ncols;

      /* print message if relatively large LP is solved from scratch, since this could lead to a longer period during
       * which the user sees no output; more detailed probing stats only in debug mode */
      ncols = SCIPgetNLPCols(scip);
      if( !SCIPisLPSolBasic(scip) && ncols > 1000 )
      {
         int nunfixedcols = SCIPgetNUnfixedLPCols(scip);

         if( nunfixedcols > 0.5 * ncols )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
               "Heuristic " HEUR_NAME " solving LP from scratch with %.1f %% unfixed columns (%d of %d) ...\n",
               100.0 * (nunfixedcols / (SCIP_Real)ncols), nunfixedcols, ncols);
         }
      }
      SCIPdebugMsg(scip, "Heuristic " HEUR_NAME " probing LP: %s\n",
         SCIPsnprintfProbingStats(scip, strbuf, SCIP_MAXSTRLEN));

      /* solve LP; errors in the LP solver should not kill the overall solving process, if the LP is just needed for a
       * heuristic.  hence in optimized mode, the return code is caught and a warning is printed, only in debug mode,
       * SCIP will stop.
       */
      SCIPdebugMsg(scip, "starting solving bound-heur LP at time %g, LP iterations: %" SCIP_LONGINT_FORMAT "\n",
         SCIPgetSolvingTime(scip), SCIPgetNLPIterations(scip));
#ifdef NDEBUG
      {
         SCIP_Bool retstat;
         retstat = SCIPsolveProbingLP(scip, -1, &lperror, NULL);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while solving LP in bound heuristic; LP solve terminated with code <%d>\n",
               retstat);
         }
      }
#else
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
#endif
      SCIPdebugMsg(scip, "ending solving bound-heur LP at time %g\n", SCIPgetSolvingTime(scip));

      lpstatus = SCIPgetLPSolstat(scip);

      SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
      SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, lpstatus);

      /* check if this is a feasible solution */
      if( lpstatus == SCIP_LPSOLSTAT_OPTIMAL && !lperror )
      {
         SCIP_SOL* newsol;
         SCIP_Bool stored;
         SCIP_Bool success;

         /* create temporary solution */
         SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

         /* copy the current LP solution to the working solution */
         SCIP_CALL( SCIPlinkLPSol(scip, newsol) );

         SCIP_CALL( SCIProundSol(scip, newsol, &success) );

         if( success )
         {
            SCIPdebugMsg(scip, "bound heuristic found roundable primal solution: obj=%g\n",
               SCIPgetSolOrigObj(scip, newsol));

            /* check solution for feasibility, and add it to solution store if possible.
             * Neither integrality nor feasibility of LP rows have to be checked, because they
             * are guaranteed by the heuristic at this stage.
             */
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPtrySol(scip, newsol, TRUE, TRUE, TRUE, TRUE, TRUE, &stored) );
#else
            SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, TRUE, FALSE, FALSE, &stored) );
#endif

            if( stored )
            {
               SCIPdebugMsg(scip, "found feasible solution:\n");
               *result = SCIP_FOUNDSOL;
            }
         }

         /* free solution */
         SCIP_CALL( SCIPfreeSol(scip, &newsol) );
      }
   }

   /* exit probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyBound)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of heuristic */
   SCIP_CALL( SCIPincludeHeurBound(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeBound)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecBound)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNPseudoBranchCands(scip) == 0 )
      return SCIP_OKAY;

   if( !SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DIDNOTFIND;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* stop execution method if there is already a primal feasible solution at hand */
   if( SCIPgetBestSol(scip) != NULL && heurdata->onlywithoutsol )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "apply bound heuristic at node %lld\n",
      SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_Bool cutoff;

      SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

      /* manually cut off the node if the LP construction detected infeasibility (heuristics cannot return such a result) */
      if( cutoff )
      {
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetCurrentNode(scip)) );
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPflushLP(scip) );
   }

   if( heurdata->bound == 'l' || heurdata->bound == 'b' )
   {
      SCIP_CALL(applyBoundHeur(scip, heur, heurdata, TRUE, result) );
   }
   if( heurdata->bound == 'u' || heurdata->bound == 'b' )
   {
      SCIP_CALL(applyBoundHeur(scip, heur, heurdata, FALSE, result) );
   }

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the bound primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurBound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create bound primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecBound, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyBound) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeBound) );

   /* add bound heuristic parameters */

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/onlywithoutsol",
         "Should heuristic only be executed if no primal solution was found, yet?",
         &heurdata->onlywithoutsol, TRUE, DEFAULT_ONLYWITHOUTSOL, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxproprounds",
         "maximum number of propagation rounds during probing (-1 infinity, -2 parameter settings)",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX/4, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/bound",
         "to which bound should integer variables be fixed? ('l'ower, 'u'pper, or 'b'oth)",
         &heurdata->bound, FALSE, DEFAULT_BOUND, "lub", NULL, NULL) );

   return SCIP_OKAY;
}
