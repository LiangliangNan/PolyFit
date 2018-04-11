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

/**@file   heur_bound.c
 * @brief  heuristic which fixes all integer variables to a bound (lower/upper) and solves the remaining LP
 * @author Gerald Gamrath
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/heur_bound.h"

#ifdef SCIP_STATISTIC
#include "scip/clock.h"
#endif

#define HEUR_NAME             "bound"
#define HEUR_DESC             "heuristic which fixes all integer variables to a bound and solves the remaining LP"
#define HEUR_DISPCHAR         'H'
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
      SCIP_LPSOLSTAT lpstatus;
      SCIP_Bool lperror;

      SCIPdebugMsg(scip, "starting solving bound-heur LP at time %g, LP iterations: %" SCIP_LONGINT_FORMAT "\n",
         SCIPgetSolvingTime(scip), SCIPgetNLPIterations(scip));

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
