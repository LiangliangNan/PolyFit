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

/**@file   heur_randrounding.c
 * @brief  randomized LP rounding heuristic which also generates conflicts via an auxiliary probing tree
 * @author Gregor Hendel
 *
 * Randomized LP rounding uses a random variable from a uniform distribution
 * over [0,1] to determine whether the fractional LP value x should be rounded
 * up with probability x - floor(x) or down with probability ceil(x) - x.
 *
 * This implementation uses domain propagation techniques to tighten the variable domains after every
 * rounding step.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_randrounding.h"
#include "scip/pub_misc.h"


#define HEUR_NAME             "randrounding"
#define HEUR_DESC             "fast LP rounding heuristic"
#define HEUR_DISPCHAR         'G'
#define HEUR_PRIORITY         -200
#define HEUR_FREQ             20
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_ONCEPERNODE   FALSE          /**< should the heuristic only be called once per node? */
#define DEFAULT_RANDSEED         23          /**< default random seed */
#define DEFAULT_USESIMPLEROUNDING    FALSE   /**< should the heuristic apply the variable lock strategy of simple rounding,
                                               *  if possible? */
#define DEFAULT_MAXPROPROUNDS 1              /**< limit of rounds for each propagation call */
#define DEFAULT_PROPAGATEONLYROOT TRUE       /**< should the probing part of the heuristic be applied exclusively at the root node */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generation */
   SCIP_Longint          lastlp;             /**< last LP number where the heuristic was applied */
   int                   maxproprounds;      /**< limit of rounds for each propagation call */
   SCIP_Bool             oncepernode;        /**< should the heuristic only be called once per node? */
   SCIP_Bool             usesimplerounding;  /**< should the heuristic apply the variable lock strategy of simple rounding,
                                               *  if possible? */
   SCIP_Bool             propagateonlyroot;  /**< should the probing part of the heuristic be applied exclusively at the root node? */
};

/*
 * Local methods
 */

/** perform randomized rounding of the given solution. Domain propagation is optionally applied after every rounding
 *  step
 */
static
SCIP_RETCODE performRandRounding(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_SOL*             sol,                /**< solution to round */
   SCIP_VAR**            cands,              /**< candidate variables */
   int                   ncands,             /**< number of candidates */
   SCIP_Bool             propagate,          /**< should the rounding be propagated? */
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   int c;
   SCIP_Bool stored;
   SCIP_VAR** permutedcands;
   SCIP_Bool cutoff;

   assert(heurdata != NULL);

   /* start probing tree before rounding begins */
   if( propagate )
   {
      SCIP_CALL( SCIPstartProbing(scip) );
      SCIPenableVarHistory(scip);
   }

   /* copy and permute the candidate array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &permutedcands, cands, ncands) );

   assert(permutedcands != NULL);

   SCIPrandomPermuteArray(heurdata->randnumgen, (void **)permutedcands, 0, ncands);
   cutoff = FALSE;

   /* loop over candidates and perform randomized rounding and optionally probing. */
   for (c = 0; c < ncands && !cutoff; ++c)
   {
      SCIP_VAR* var;
      SCIP_Real oldsolval;
      SCIP_Real newsolval;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Longint ndomreds;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real ceilval;
      SCIP_Real floorval;

      /* get next variable from permuted candidate array */
      var = permutedcands[c];
      oldsolval = SCIPgetSolVal(scip, sol, var);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      assert( ! SCIPisFeasIntegral(scip, oldsolval) );
      assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );

      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      ceilval = SCIPfeasCeil(scip, oldsolval);
      floorval = SCIPfeasFloor(scip, oldsolval);

      SCIPdebugMsg(scip, "rand rounding heuristic: var <%s>, val=%g, rounddown=%u, roundup=%u\n",
         SCIPvarGetName(var), oldsolval, mayrounddown, mayroundup);

      /* abort if rounded ceil and floor value lie outside the variable domain. Otherwise, check if
       * bounds allow only one rounding direction, anyway */
      if( lb > ceilval + 0.5 || ub < floorval - 0.5 )
      {
         cutoff = TRUE;
         break;
      }
      else if( SCIPisFeasEQ(scip, lb, ceilval) )
      {
         /* only rounding up possible */
         assert(SCIPisFeasGE(scip, ub, ceilval));
         newsolval = ceilval;
      }
      else if( SCIPisFeasEQ(scip, ub, floorval) )
      {
         /* only rounding down possible */
         assert(SCIPisFeasLE(scip,lb, floorval));
         newsolval = floorval;
      }
      else if( !heurdata->usesimplerounding || !(mayroundup || mayrounddown) )
      {
         /* the standard randomized rounding */
         SCIP_Real randnumber;

         randnumber = SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0);
         if( randnumber <= oldsolval - floorval )
            newsolval = ceilval;
         else
            newsolval = floorval;
      }
      /* choose rounding direction, if possible, or use the only direction guaranteed to be feasible */
      else if( mayrounddown && mayroundup )
      {
         /* we can round in both directions: round in objective function direction */
         if ( SCIPvarGetObj(var) >= 0.0 )
            newsolval = floorval;
         else
            newsolval = ceilval;
      }
      else if( mayrounddown )
         newsolval = floorval;
      else
      {
         assert(mayroundup);
         newsolval = ceilval;
      }

      assert(SCIPisFeasLE(scip, lb, newsolval));
      assert(SCIPisFeasGE(scip, ub, newsolval));

      /* if propagation is enabled, fix the candidate variable to its rounded value and propagate the solution */
      if( propagate )
      {
         SCIP_Bool lbadjust;
         SCIP_Bool ubadjust;

         lbadjust = SCIPisGT(scip, newsolval, lb);
         ubadjust = SCIPisLT(scip, newsolval, ub);

         assert( lbadjust || ubadjust || SCIPisFeasEQ(scip, lb, ub));

         /* enter a new probing node if the variable was not already fixed before */
         if( lbadjust || ubadjust )
         {
            if( SCIPisStopped(scip) )
               break;

            /* We only want to create a new probing node if we do not exceeed the maximal tree depth,
             * otherwise we finish at this point.
             * @todo: Maybe we want to continue with the same node because we do not backtrack.
             */
            if( SCIPgetDepth(scip) < SCIP_MAXTREEDEPTH )
            {
               SCIP_CALL( SCIPnewProbingNode(scip) );
            }
            else
               break;

            /* tighten the bounds to fix the variable for the probing node */
            if( lbadjust )
            {
               SCIP_CALL( SCIPchgVarLbProbing(scip, var, newsolval) );
            }
            if( ubadjust )
            {
               SCIP_CALL( SCIPchgVarUbProbing(scip, var, newsolval) );
            }

            /* call propagation routines for the reduced problem */
            SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, &ndomreds) );
         }
      }
      /* store new solution value */
      SCIP_CALL( SCIPsetSolVal(scip, sol, var, newsolval) );
   }

   /* if no cutoff was detected, the solution is a candidate to be checked for feasibility */
   if( !cutoff && ! SCIPisStopped(scip) )
   {
      if( SCIPallColsInLP(scip) )
      {
         /* check solution for feasibility, and add it to solution store if possible
          * neither integrality nor feasibility of LP rows has to be checked, because all fractional
          * variables were already moved in feasible direction to the next integer
          */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, TRUE, &stored) );
      }
      else
      {
         /* if there are variables which are not present in the LP, e.g., for
          * column generation, we need to check their bounds
          */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, FALSE, TRUE, &stored) );
      }

      if( stored )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "found feasible rounded solution:\n");
         SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
         *result = SCIP_FOUNDSOL;
      }
   }

   assert( !propagate || SCIPinProbing(scip) );

   /* exit probing mode and free locally allocated memory */
   if( propagate )
   {
      SCIP_CALL( SCIPendProbing(scip) );
   }

   SCIPfreeBufferArray(scip, &permutedcands);

   return SCIP_OKAY;
}

/** perform randomized LP-rounding */
static
SCIP_RETCODE performLPRandRounding(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing mask */
   SCIP_Bool             propagate,          /**< should the heuristic apply SCIP's propagation? */
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** lpcands;
   SCIP_Longint nlps;
   int nlpcands;

   /* only call heuristic, if an optimal LP solution is at hand */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get fractional variables, that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, &nlpcands, NULL, NULL) );

   /* only call heuristic, if LP solution is fractional; except we are called during pricing, in this case we
    * want to detect a (mixed) integer (LP) solution which is primal feasible */
   if ( nlpcands == 0  && heurtiming != SCIP_HEURTIMING_DURINGPRICINGLOOP )
      return SCIP_OKAY;

   /* get the working solution from heuristic's local data */
   sol = heurdata->sol;
   assert( sol != NULL );

   /* copy the current LP solution to the working solution */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   /* don't call heuristic, if we have already processed the current LP solution */
   nlps = SCIPgetNLPs(scip);
   if( nlps == heurdata->lastlp )
      return SCIP_OKAY;
   heurdata->lastlp = nlps;

   *result = SCIP_DIDNOTFIND;

   /* perform random rounding */
   SCIPdebugMsg(scip, "executing rand LP-rounding heuristic: %d fractionals\n", nlpcands);
   SCIP_CALL( performRandRounding(scip, heurdata, sol, lpcands, nlpcands, propagate, result) );

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRandrounding)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRandrounding(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRandrounding) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitRandrounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create heuristic data */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   heurdata->lastlp = -1;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRandrounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolRandrounding)
{
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   heurdata->lastlp = -1;

   /* change the heuristic's timingmask, if it should be called only once per node */
   if( heurdata->oncepernode )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_AFTERLPNODE);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolRandrounding)
{
   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRandrounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool propagate;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand or if relaxation solution is available */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* don't call heuristic, if we have already processed the current LP solution but no relaxation solution is available */
   if ( SCIPgetNLPs(scip) == heurdata->lastlp && ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;

   propagate = (!heurdata->propagateonlyroot || SCIPgetDepth(scip) == 0);

   /* try to round LP solution */
   SCIP_CALL( performLPRandRounding(scip, heurdata, heurtiming, propagate, result) );

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

/** creates the rand rounding heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRandrounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRandrounding, heurdata) );
   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRandrounding) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRandrounding) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRandrounding) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolRandrounding) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolRandrounding) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRandrounding) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/oncepernode",
         "should the heuristic only be called once per node?",
         &heurdata->oncepernode, TRUE, DEFAULT_ONCEPERNODE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usesimplerounding", "should the heuristic apply the variable lock strategy of simple rounding, if possible?",
         &heurdata->usesimplerounding, TRUE, DEFAULT_USESIMPLEROUNDING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/propagateonlyroot",
         "should the probing part of the heuristic be applied exclusively at the root node?",
         &heurdata->propagateonlyroot, TRUE, DEFAULT_PROPAGATEONLYROOT, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxproprounds",
         "limit of rounds for each propagation call",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS,
         -1, INT_MAX, NULL, NULL) );
   return SCIP_OKAY;
}
