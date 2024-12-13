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

/**@file   heur_conflictdiving.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  LP diving heuristic that chooses fixings w.r.t. conflict locks
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/heur_conflictdiving.h"
#include "scip/heuristics.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include <string.h>

#define HEUR_NAME                    "conflictdiving"
#define HEUR_DESC                    "LP diving heuristic that chooses fixings w.r.t. conflict locks"
#define HEUR_DISPCHAR                SCIP_HEURDISPCHAR_DIVING
#define HEUR_PRIORITY                -1000100
#define HEUR_FREQ                    10
#define HEUR_FREQOFS                 0
#define HEUR_MAXDEPTH                -1
#define HEUR_TIMING                  SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP             FALSE  /**< does the heuristic use a secondary SCIP instance? */
#define DIVESET_DIVETYPES            SCIP_DIVETYPE_INTEGRALITY | SCIP_DIVETYPE_SOS1VARIABLE /**< bit mask that represents all supported dive types */
#define DIVESET_ISPUBLIC             FALSE  /**< is this dive set publicly available (ie., can be used by other primal heuristics?) */
#define DEFAULT_RANDSEED             151 /**< default random seed */

/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.15 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.15/**< percentage of immediate domain changes during probing to trigger LP resolve */
#define DEFAULT_LPSOLVEFREQ           0 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS FALSE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_LOCKWEIGHT         0.75 /**< weight used in a convex combination of conflict and variable locks */
#define DEFAULT_LIKECOEF          FALSE /**< perform rounding like coefficient diving */
#define DEFAULT_MAXVIOL            TRUE /**< prefer rounding direction with most violation */
#define DEFAULT_MINCONFLICTLOCKS      5 /**< threshold for penalizing the score */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             lockweight;         /**< weight factor to combine conflict and variable locks */
   SCIP_Bool             likecoefdiving;     /**< use the same rounding strategy like coefficent diving */
   SCIP_Bool             maxviol;            /**< rounding into potentially infeasible direction */
   int                   minconflictlocks;   /**< threshold for penalizing the score */
};

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyConflictdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeHeurConflictdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeConflictdiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitConflictdiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXIT(heurExitConflictdiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXEC(heurExecConflictdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   *result = SCIP_DELAYED;

   /* don't run if no conflict constraints where found */
   if( SCIPgetNConflictConssFound(scip) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible, -1L, SCIP_DIVECONTEXT_SINGLE) );

   return SCIP_OKAY;
}

#define MIN_RAND 1e-06
#define MAX_RAND 1e-05

/** calculate score variant 1: use rounding strategy like coefficent diving */
static
SCIP_RETCODE getScoreLikeCoefdiving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_RANDNUMGEN*      rng,                /**< random number generator of the diveset */
   SCIP_DIVETYPE         divetype,           /**< divetype of the heuristic */
   SCIP_VAR*             cand,               /**< diving candidate */
   SCIP_Real             candsol,            /**< diving candidate solution */
   SCIP_Real             candsfrac,          /**< fractionality of the candidate solution */
   SCIP_Real*            score,              /**< pointer to store the score */
   SCIP_Bool*            roundup             /**< pointer to store whether the candidate should be rounded upwards */
   )
{
   SCIP_Real upweight;
   SCIP_Real downweight;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   int nconflictlocksdown;
   int nconflictlocksup;
   int nlocksdown;
   int nlocksup;

   /* get conflict locks */
   nconflictlocksup = SCIPvarGetNLocksUpType(cand, SCIP_LOCKTYPE_CONFLICT);
   nconflictlocksdown = SCIPvarGetNLocksDownType(cand, SCIP_LOCKTYPE_CONFLICT);

   /* get variable locks */
   nlocksup = SCIPvarGetNLocksUpType(cand, SCIP_LOCKTYPE_MODEL);
   nlocksdown = SCIPvarGetNLocksDownType(cand, SCIP_LOCKTYPE_MODEL);

   /* combine conflict and variable locks */
   upweight = heurdata->lockweight * nconflictlocksup +  (1.0 - heurdata->lockweight) * nlocksup;
   downweight = heurdata->lockweight * nconflictlocksdown + (1.0 - heurdata->lockweight) * nlocksdown;

   /* check whether there exists a direction w/o any locks */
   mayrounddown = SCIPisZero(scip, upweight);
   mayroundup = SCIPisZero(scip, downweight);

   if( mayrounddown || mayroundup )
   {
      /* choose rounding direction:
       * - if variable may be rounded in both directions, round corresponding to the fractionality
       * - otherwise, round in the infeasible direction
       */
      if( mayrounddown && mayroundup )
      {
         assert(divetype != SCIP_DIVETYPE_SOS1VARIABLE || heurdata->lockweight > 0);

         /* try to avoid variability; decide randomly if the LP solution can contain some noise */
         if( SCIPisEQ(scip, candsfrac, 0.5) )
            *roundup = (SCIPrandomGetInt(rng, 0, 1) == 0);
         else
            *roundup = (candsfrac > 0.5);
      }
      else
         *roundup = mayrounddown;
   }
   else
   {
      /* the candidate may not be rounded */
      *roundup = (SCIPisGT(scip, downweight, upweight) || (SCIPisEQ(scip, downweight, upweight) && candsfrac > 0.5));
   }

   if( *roundup )
   {
      switch( divetype )
      {
         case SCIP_DIVETYPE_INTEGRALITY:
            candsfrac = 1.0 - candsfrac;
            break;
         case SCIP_DIVETYPE_SOS1VARIABLE:
            if( SCIPisFeasPositive(scip, candsol) )
               candsfrac = 1.0 - candsfrac;
            break;
         default:
            SCIPerrorMessage("Error: Unsupported diving type\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
      } /*lint !e788*/

      /* add some noise to avoid ties */
      *score = upweight + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND);
   }
   else
   {
      if( divetype == SCIP_DIVETYPE_SOS1VARIABLE && SCIPisFeasNegative(scip, candsol) )
         candsfrac = 1.0 - candsfrac;

      /* add some noise to avoid ties */
      *score = downweight + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND);
   }

   /* penalize too small fractions */
   if( SCIPisEQ(scip, candsfrac, 0.01) )
   {
      /* try to avoid variability; decide randomly if the LP solution can contain some noise.
       * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for scaling the score
       */
      if( SCIPrandomGetInt(rng, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
         (*score) *= 0.01;
   }
   else if( candsfrac < 0.01 )
      (*score) *= 0.1;

   /* prefer decisions on binary variables */
   if( !SCIPvarIsBinary(cand) )
      *score = -1.0 / *score;

   return SCIP_OKAY;
}

/** calculate score variant 2: use a rounding strategy that tends towards infeasibility */
static
SCIP_RETCODE getScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_RANDNUMGEN*      rng,                /**< random number generator of the diveset */
   SCIP_DIVETYPE         divetype,           /**< divetype of the heuristic */
   SCIP_VAR*             cand,               /**< diving candidate */
   SCIP_Real             candsol,            /**< diving candidate solution */
   SCIP_Real             candsfrac,          /**< fractionality of the candidate solution */
   SCIP_Real*            score,              /**< pointer to store the score */
   SCIP_Bool*            roundup             /**< pointer to store whether the candidate should be rounded upwards */
   )
{
   SCIP_Real conflictlocksum;
   SCIP_Real upweight;
   SCIP_Real downweight;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   int nlocksup;
   int nlocksdown;
   int nconflictlocksup;
   int nconflictlocksdown;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(rng != NULL);

   /* get conflict locks */
   nconflictlocksup = SCIPvarGetNLocksUpType(cand, SCIP_LOCKTYPE_CONFLICT);
   nconflictlocksdown = SCIPvarGetNLocksDownType(cand, SCIP_LOCKTYPE_CONFLICT);
   conflictlocksum = nconflictlocksup + nconflictlocksdown;

   /* get variable locks */
   nlocksup = SCIPvarGetNLocksUpType(cand, SCIP_LOCKTYPE_MODEL);
   nlocksdown = SCIPvarGetNLocksDownType(cand, SCIP_LOCKTYPE_MODEL);

   /* combine conflict and variable locks */
   upweight = heurdata->lockweight * nconflictlocksup +  (1.0 - heurdata->lockweight) * nlocksup;
   downweight = heurdata->lockweight * nconflictlocksdown + (1.0 - heurdata->lockweight) * nlocksdown;

   /* check whether there exists a rounding direction w/o any locks */
   mayrounddown = SCIPisZero(scip, upweight);
   mayroundup = SCIPisZero(scip, downweight);

   /* variable can be rounded in exactly one direction and we try to go into the feasible direction */
   if( mayrounddown || mayroundup )
   {
      /* choose rounding direction:
       * - if variable may be rounded in both directions, round corresponding to the fractionality
       * - otherwise, round in the feasible direction
       */
      if( mayrounddown && mayroundup )
      {
         assert(divetype != SCIP_DIVETYPE_SOS1VARIABLE || heurdata->lockweight > 0);

         /* try to avoid variability; decide randomly if the LP solution can contain some noise */
         if( SCIPisEQ(scip, candsfrac, 0.5) )
            *roundup = (SCIPrandomGetInt(rng, 0, 1) == 0);
         else
            *roundup = (candsfrac > 0.5);
      }
      else
         *roundup = mayroundup;
   }
   else
   {
      assert(!mayrounddown);

      /* both rounding directions have a different amount of locks */
      if( !SCIPisEQ(scip, upweight, downweight) )
      {
         *roundup = (heurdata->maxviol ? SCIPisGT(scip, upweight, downweight) : SCIPisLT(scip, upweight, downweight));
      }
      /* break ties with lp fractionality != 0.5 */
      else if( !SCIPisEQ(scip, candsfrac, 0.5) )
      {
         *roundup = (candsfrac > 0.5);
      }
      /* break tie randomly */
      else
      {
         *roundup = (SCIPrandomGetInt(rng, 0, 1) == 1);
      }
   }

   if( *roundup )
   {
      switch( divetype )
      {
         case SCIP_DIVETYPE_INTEGRALITY:
            candsfrac = 1.0 - candsfrac;
            break;
         case SCIP_DIVETYPE_SOS1VARIABLE:
            if( SCIPisFeasPositive(scip, candsol) )
               candsfrac = 1.0 - candsfrac;
            break;
         default:
            SCIPerrorMessage("Error: Unsupported diving type\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
      } /*lint !e788*/

      /* add some noise to avoid ties */
      *score = upweight + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND);
   }
   else
   {
      if( divetype == SCIP_DIVETYPE_SOS1VARIABLE && SCIPisFeasNegative(scip, candsol) )
         candsfrac = 1.0 - candsfrac;

      /* add some noise to avoid ties */
      *score = downweight + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND);
   }

   /* penalize too few conflict locks */
   if( conflictlocksum > 0 && conflictlocksum < heurdata->minconflictlocks )
      (*score) *= 0.1;

   /* penalize if no conflict locks exist at all */
   if( conflictlocksum == 0 )
      (*score) *= 0.01;

   /* penalize too small fractions */
   if( SCIPisEQ(scip, candsfrac, 0.01) )
   {
      /* try to avoid variability; decide randomly if the LP solution can contain some noise.
       * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for scaling the score
       */
      if( SCIPrandomGetInt(rng, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
         (*score) *= 0.01;
   }
   else if( candsfrac < 0.01 )
      (*score) *= 0.01;

   /* prefer decisions on binary variables */
   if( !SCIPvarIsBinary(cand) )
      *score = -1.0 / *score;

   return SCIP_OKAY;
}


/** returns a score for the given candidate -- the best candidate maximizes the diving score */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreConflictdiving)
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   SCIP_RANDNUMGEN* rng;

   rng = SCIPdivesetGetRandnumgen(diveset);
   assert(rng != NULL);

   heur = SCIPdivesetGetHeur(diveset);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( heurdata->likecoefdiving )
   {
      SCIP_CALL( getScoreLikeCoefdiving(scip, heurdata, rng, divetype, cand, candsol, candsfrac, score, roundup) );
   }
   else
   {
      SCIP_CALL( getScore(scip, heurdata, rng, divetype, cand, candsol, candsfrac, score, roundup) );
   }

   /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
   assert( (0.0 < candsfrac && candsfrac < 1.0) || SCIPvarIsBinary(cand) || divetype == SCIP_DIVETYPE_SOS1VARIABLE );

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

#define divesetAvailableConflictdiving NULL

/** creates the conflictdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurConflictdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create conflictdiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ,
         HEUR_FREQOFS, HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecConflictdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyConflictdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeConflictdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitConflictdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitConflictdiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_LPRESOLVEDOMCHGQUOT,
         DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED, DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS,
         DIVESET_ISPUBLIC, DIVESET_DIVETYPES, divesetGetScoreConflictdiving, divesetAvailableConflictdiving) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/maxviol", "try to maximize the violation",
         &heurdata->maxviol, TRUE, DEFAULT_MAXVIOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/likecoef",
         "perform rounding like coefficient diving",
         &heurdata->likecoefdiving, TRUE, DEFAULT_LIKECOEF, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minconflictlocks",
         "minimal number of conflict locks per variable",
         &heurdata->minconflictlocks, TRUE, DEFAULT_MINCONFLICTLOCKS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lockweight",
         "weight used in a convex combination of conflict and variable locks",
         &heurdata->lockweight, TRUE, DEFAULT_LOCKWEIGHT, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
