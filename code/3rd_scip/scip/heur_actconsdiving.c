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

/**@file   heur_actconsdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the active constraints the variable appear in
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_actconsdiving.h"

#define HEUR_NAME             "actconsdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the active constraints"
#define HEUR_DISPCHAR         'a'
#define HEUR_PRIORITY         -1003700
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          5
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
#define DEFAULT_MAXDIVEUBQUOTNOSOL  1.0 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 1.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.15 /**< percentage of immediate domain changes during probing to trigger LP resolve */
#define DEFAULT_LPSOLVEFREQ           0 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS  TRUE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_RANDSEED            149 /**< default random seed */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
};


/*
 * local methods
 */

/** returns a score value for the given variable based on the active constraints that the variable appears in */
static
SCIP_Real getNActiveConsScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< working solution */
   SCIP_VAR*             var,                /**< variable to get the score value for */
   SCIP_Real*            downscore,          /**< pointer to store the score for branching downwards */
   SCIP_Real*            upscore             /**< pointer to store the score for branching upwards */
   )
{
   SCIP_COL* col;
   SCIP_ROW** rows;
   SCIP_Real* vals;
   int nrows;
   int r;
   int nactrows;
   SCIP_Real nlprows;
   SCIP_Real downcoefsum;
   SCIP_Real upcoefsum;
   SCIP_Real score;

   assert(downscore != NULL);
   assert(upscore != NULL);

   *downscore = 0.0;
   *upscore = 0.0;
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return 0.0;

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   rows = SCIPcolGetRows(col);
   vals = SCIPcolGetVals(col);
   nrows = SCIPcolGetNLPNonz(col);
   nactrows = 0;
   downcoefsum = 0.0;
   upcoefsum = 0.0;
   for( r = 0; r < nrows; ++r )
   {
      SCIP_ROW* row;
      SCIP_Real activity;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real dualsol;

      row = rows[r];
      /* calculate number of active constraint sides, i.e., count equations as two */
      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);

      /* @todo this is suboptimal because activity is calculated by looping over all nonzeros of this row, need to
       * store LP activities instead (which cannot be retrieved if no LP was solved at this node)
       */
      activity = SCIPgetRowSolActivity(scip, row, sol);

      dualsol = SCIProwGetDualsol(row);
      if( SCIPisFeasEQ(scip, activity, lhs) )
      {
         SCIP_Real coef;

         nactrows++;
         coef = vals[r] / SCIProwGetNorm(row);
         if( SCIPisFeasPositive(scip, dualsol) )
         {
            if( coef > 0.0 )
               downcoefsum += coef;
            else
               upcoefsum -= coef;
         }
      }
      else if( SCIPisFeasEQ(scip, activity, rhs) )
      {
         SCIP_Real coef;

         nactrows++;
         coef = vals[r] / SCIProwGetNorm(row);
         if( SCIPisFeasNegative(scip, dualsol) )
         {
            if( coef > 0.0 )
               upcoefsum += coef;
            else
               downcoefsum -= coef;
         }
      }
   }

   /* use the number of LP rows for normalization */
   nlprows = (SCIP_Real)SCIPgetNLPRows(scip);
   upcoefsum /= nlprows;
   downcoefsum /= nlprows;

   /* calculate the score using SCIP's branch score. Pass NULL as variable to not have the var branch factor influence
    * the result
    */
   score = nactrows / nlprows + SCIPgetBranchScore(scip, NULL, downcoefsum, upcoefsum);

   assert(score <= 3.0);
   assert(score >= 0.0);

   *downscore = downcoefsum;
   *upscore = upcoefsum;

   return score;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyActconsdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurActconsdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeActconsdiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitActconsdiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXIT(heurExitActconsdiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXEC(heurExecActconsdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   heurdata = SCIPheurGetData(heur);

   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   *result = SCIP_DIDNOTRUN;

   /* if there are no integer variables, stop execution */
   if( SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );

   return SCIP_OKAY;
}

/** calculate score and preferred rounding direction for the candidate variable; the best candidate maximizes the
 *  score
 */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreActconsdiving)
{
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Real downscore;
   SCIP_Real upscore;

   mayrounddown = SCIPvarMayRoundDown(cand);
   mayroundup = SCIPvarMayRoundUp(cand);

   /* first, calculate the variable score */
   assert(SCIPdivesetGetWorkSolution(diveset) != NULL);
   *score = getNActiveConsScore(scip, SCIPdivesetGetWorkSolution(diveset), cand, &downscore, &upscore);

   /* get the rounding direction: prefer an unroundable direction */
   if( mayrounddown && mayroundup )
   {
      /* try to avoid variability; decide randomly if the LP solution can contain some noise */
      if( SCIPisEQ(scip, candsfrac, 0.5) )
         *roundup = (SCIPrandomGetInt(SCIPdivesetGetRandnumgen(diveset), 0, 1) == 0);
      else
         *roundup = (candsfrac > 0.5);
   }
   else if( mayrounddown || mayroundup )
      *roundup = mayrounddown;
   else
      *roundup = (downscore > upscore);

   if( *roundup )
      candsfrac = 1.0 - candsfrac;

   /* penalize too small fractions */
   if( SCIPisEQ(scip, candsfrac, 0.01) )
   {
      /* try to avoid variability; decide randomly if the LP solution can contain some noise.
       * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for scaling the score
       */
      if( SCIPrandomGetInt(SCIPdivesetGetRandnumgen(diveset), 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
         (*score) *= 0.01;
   }
   else if( candsfrac < 0.01 )
      (*score) *= 0.01;

   /* prefer decisions on binary variables */
   if( !SCIPvarIsBinary(cand) )
      (*score) *= 0.01;

   /* penalize variable if it may be rounded */
   if( mayrounddown || mayroundup )
      *score -= 3.0;

   assert(!(mayrounddown || mayroundup) || *score <= 0.0);

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

/** creates the actconsdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurActconsdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create actconsdiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecActconsdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyActconsdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeActconsdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitActconsdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitActconsdiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_LPRESOLVEDOMCHGQUOT,
         DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED,
         DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS, DIVESET_DIVETYPES, divesetGetScoreActconsdiving) );

   return SCIP_OKAY;
}

