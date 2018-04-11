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

/**@file   heur_fixandinfer.c
 * @brief  fix-and-infer primal heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_fixandinfer.h"


#define HEUR_NAME             "fixandinfer"
#define HEUR_DESC             "iteratively fixes variables and propagates inferences"
#define HEUR_DISPCHAR         'I'
#define HEUR_PRIORITY         -500000
#define HEUR_FREQ             -1        /* at the moment, the heuristic seems to be useless */
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define MAXDIVEDEPTH          100


/*
 * Default parameter settings
 */

#define DEFAULT_PROPROUNDS           0  /**< maximal number of propagation rounds in probing subproblems */
#define DEFAULT_MINFIXINGS         100  /**< minimal number of fixings to apply before dive may be aborted */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   proprounds;         /**< maximal number of propagation rounds in probing subproblems */
   int                   minfixings;         /**< minimal number of fixings to apply before dive may be aborted */
};


/*
 * Local methods
 */

/** selects a variable and fixes it to its current pseudo solution value */
static
SCIP_RETCODE fixVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            pseudocands,        /**< array of unfixed variables */
   int                   npseudocands,       /**< number of unfixed variables */
   SCIP_Real             large               /**< large value to be used instead of infinity */
   )
{
   SCIP_VAR* var;
   SCIP_Real bestscore;
   SCIP_Real score;
   SCIP_Real solval;
   int bestcand;
   int ncands;
   int c;

   assert(pseudocands != NULL);
   assert(npseudocands > 0);

   /* if existing, choose one of the highest priority binary variables; if no high priority binary variables
    * exist, choose a variable among all unfixed integral variables
    */
   ncands = SCIPgetNPrioPseudoBranchBins(scip);
   if( ncands == 0 )
      ncands = npseudocands;

   /* select variable to tighten the domain for */
   bestscore = -SCIPinfinity(scip);
   bestcand = -1;
   for( c = 0; c < ncands; ++c )
   {
      score = SCIPgetVarAvgInferenceScore(scip, pseudocands[c]);
      if( score > bestscore )
      {
         bestscore = score;
         bestcand = c;
      }
   }
   assert(bestcand != -1);

   /* fix variable to its current pseudo solution value */
   var = pseudocands[bestcand];
   solval = SCIPgetVarSol(scip, var);

   /* adapt solution value if it is infinite */
   if( SCIPisInfinity(scip, solval) )
   {
      SCIP_Real lb;
      assert(SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
      lb = SCIPvarGetLbLocal(var);

      /* adapt fixing value by changing it to a large value */
      if( SCIPisInfinity(scip, -lb) )
         solval = SCIPceil(scip, large);
      else if( !SCIPisInfinity(scip, SCIPceil(scip, lb+large)) )
         solval = SCIPceil(scip, lb+large);
   }
   else if( SCIPisInfinity(scip, -solval) )
   {
      SCIP_Real ub;
      assert(SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
      ub = SCIPvarGetUbLocal(var);

      /* adapt fixing value by changing it to a large negative value */
      if( SCIPisInfinity(scip, ub) )
         solval = SCIPfloor(scip, -large);
      else if( !SCIPisInfinity(scip, -SCIPfloor(scip, ub-large)) )
         solval = SCIPfloor(scip, ub-large);
   }

   assert(SCIPisFeasIntegral(scip, solval)); /* in probing, we always have the pseudo solution */
   SCIPdebugMsg(scip, " -> fixed variable <%s>[%g,%g] = %g (%d candidates left)\n",
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), solval, npseudocands - 1);
   SCIP_CALL( SCIPfixVarProbing(scip, var, solval) );

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyFixandinfer)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurFixandinfer(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeFixandinfer) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecFixandinfer)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** cands;
   int ncands;
   int startncands;
   int divedepth;
   SCIP_Bool cutoff;
   SCIP_Real large;

   *result = SCIP_DIDNOTRUN;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* we cannot run on problems with continuous variables */
   if( SCIPgetNContVars(scip) > 0 )
      return SCIP_OKAY;

   /* get unfixed variables */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, &ncands, NULL) );
   if( ncands == 0 )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* fix variables and propagate inferences as long as the problem is still feasible and there are
    * unfixed integral variables
    */
   cutoff = FALSE;
   divedepth = 0;
   startncands = ncands;

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   if( SCIP_MAXTREEDEPTH <= SCIPgetDepth(scip) )
   {
      SCIP_CALL( SCIPendProbing(scip) );
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "starting fix-and-infer heuristic with %d unfixed integral variables\n", ncands);

   *result = SCIP_DIDNOTFIND;

   /* create next probing node */
   SCIP_CALL( SCIPnewProbingNode(scip) );

   /* determine large value to set variables to */
   large = SCIPinfinity(scip);
   if( !SCIPisInfinity(scip, 0.1 / SCIPfeastol(scip)) )
      large = 0.1 / SCIPfeastol(scip);

   while( !cutoff && ncands > 0
      && (divedepth < heurdata->minfixings || (startncands - ncands) * 2 * MAXDIVEDEPTH >= startncands * divedepth)
      && !SCIPisStopped(scip) )
   {
      divedepth++;

      /* fix next variable */
      SCIP_CALL( fixVariable(scip, cands, ncands, large) );

      /* propagate the fixing */
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->proprounds, &cutoff, NULL) );

      /* get remaining unfixed variables */
      if( !cutoff )
      {
         SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, &ncands, NULL) );
      }
   }

   /* check, if we are still feasible */
   if( cutoff )
   {
      SCIPdebugMsg(scip, "propagation detected a cutoff\n");
   }
   else if( ncands == 0 )
   {
      SCIP_Bool success;

      success = FALSE;

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtryCurrentSol(scip, heur, FALSE, FALSE, FALSE, TRUE, &success) );

      if( success )
      {
         SCIPdebugMsg(scip, "found primal feasible solution\n");
         *result = SCIP_FOUNDSOL;
      }
      else
      {
         SCIPdebugMsg(scip, "primal solution was rejected\n");
      }
   }
   else
   {
      SCIPdebugMsg(scip, "probing was aborted (probing depth: %d, fixed: %d/%d)", divedepth, startncands - ncands, startncands);
   }

   /* end probing */
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the fix-and-infer primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurFixandinfer(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

    /* create Fixandinfer primal heuristic data */
    SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

    /* include primal heuristic */
    SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
          HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
          HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecFixandinfer, heurdata) );

    assert(heur != NULL);

    /* set non-NULL pointers to callback methods */
    SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyFixandinfer) );
    SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeFixandinfer) );

   /* fixandinfer heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/fixandinfer/proprounds",
         "maximal number of propagation rounds in probing subproblems (-1: no limit, 0: auto)",
         &heurdata->proprounds, TRUE, DEFAULT_PROPROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/fixandinfer/minfixings",
         "minimal number of fixings to apply before dive may be aborted",
         &heurdata->minfixings, TRUE, DEFAULT_MINFIXINGS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
