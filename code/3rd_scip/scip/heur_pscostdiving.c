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

/**@file   heur_pscostdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the pseudo cost values
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_pscostdiving.h"

#define HEUR_NAME             "pscostdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the pseudo cost values"
#define HEUR_DISPCHAR         'p'
#define HEUR_PRIORITY         -1002000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          2
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */
#define DIVESET_DIVETYPES     SCIP_DIVETYPE_INTEGRALITY /**< bit mask that represents all supported dive types */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.15 /**< percentage of immediate domain changes during probing to trigger LP resolve */
#define DEFAULT_LPSOLVEFREQ           0 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS TRUE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_RANDSEED            103 /**< initial random seed */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
};

/*
 * local methods
 */

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyPscostdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurPscostdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreePscostdiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);


   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   *result = SCIP_DIDNOTRUN;

   /* terminate if there are no integer variables (note that, e.g., SOS1 variables may be present) */
   if( SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );

   return SCIP_OKAY;
}


/** returns a score for the given candidate -- the best candidate maximizes the diving score */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScorePscostdiving)
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;
   SCIP_Real pscostquot;

   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;

   mayrounddown = SCIPvarMayRoundDown(cand);
   mayroundup = SCIPvarMayRoundUp(cand);


   /* bound fractions to not prefer variables that are nearly integral */
   candsfrac = MAX(candsfrac, 0.1);
   candsfrac = MIN(candsfrac, 0.9);

   pscostdown = SCIPgetVarPseudocostVal(scip, cand, 0.0 - candsfrac);
   pscostup = SCIPgetVarPseudocostVal(scip, cand, 1.0 - candsfrac);

   /*  determine the candidate direction. if the variable may be trivially rounded in one direction, take the other direction;
    *  otherwise, consider first the direction from the root solution, second the candidate fractionality, and
    *  last the direction of smaller pseudo costs
    *
    *  to avoid performance variability caused by numerics we use random numbers to decide whether we want to roundup or
    *  round down if the values to compare are equal within tolerances.
    */
   assert(pscostdown >= 0.0 && pscostup >= 0.0);
   if( mayrounddown != mayroundup )
      *roundup = mayrounddown;
   else if( SCIPisLT(scip, candsol, SCIPvarGetRootSol(cand) - 0.4)
         || (SCIPisEQ(scip, candsol, SCIPvarGetRootSol(cand) - 0.4) && SCIPrandomGetInt(SCIPdivesetGetRandnumgen(diveset), 0, 1) == 0) )
      *roundup = FALSE;
   else if( SCIPisGT(scip, candsol, SCIPvarGetRootSol(cand) + 0.4)
         || (SCIPisEQ(scip, candsol, SCIPvarGetRootSol(cand) + 0.4) && SCIPrandomGetInt(SCIPdivesetGetRandnumgen(diveset), 0, 1) == 0) )
      *roundup = TRUE;
   else if( SCIPisLT(scip, candsfrac, 0.3)
         || (SCIPisEQ(scip, candsfrac, 0.3) && SCIPrandomGetInt(SCIPdivesetGetRandnumgen(diveset), 0, 1) == 0) )
      *roundup = FALSE;
   else if( SCIPisGT(scip, candsfrac, 0.7)
         || (SCIPisEQ(scip, candsfrac, 0.7) && SCIPrandomGetInt(SCIPdivesetGetRandnumgen(diveset), 0, 1) == 0) )
      *roundup = TRUE;
   else if( SCIPisEQ(scip, pscostdown, pscostup) )
      *roundup = (SCIPrandomGetInt(SCIPdivesetGetRandnumgen(diveset), 0, 1) == 0);
   else if( pscostdown > pscostup )
      *roundup = TRUE;
   else
      *roundup = FALSE;

   if( *roundup )
      pscostquot = sqrt(candsfrac) * (1.0 + pscostdown) / (1.0 + pscostup);
   else
      pscostquot = sqrt(1.0 - candsfrac) * (1.0 + pscostup) / (1.0 + pscostdown);

   /* prefer decisions on binary variables */
   if( SCIPvarIsBinary(cand) && !(SCIPvarMayRoundDown(cand) || SCIPvarMayRoundUp(cand)))
      pscostquot *= 1000.0;

   assert(pscostquot >= 0);
   *score = pscostquot;

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

/** creates the pscostdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurPscostdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Pscostdiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecPscostdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyPscostdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreePscostdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitPscostdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitPscostdiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_LPRESOLVEDOMCHGQUOT,
         DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED, DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS, DIVESET_DIVETYPES, divesetGetScorePscostdiving) );

   return SCIP_OKAY;
}

