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

/**@file   heur_nlpdiving.c
 * @brief  NLP diving heuristic that chooses fixings w.r.t. the fractionalities
 * @author Timo Berthold
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_nlpdiving.h"
#include "scip/heur_subnlp.h" /* for NLP initialization */
#include "scip/heur_undercover.h" /* for cover computation */
#include "nlpi/nlpi.h" /* for NLP statistics, currently */


#define HEUR_NAME             "nlpdiving"
#define HEUR_DESC             "NLP diving heuristic that chooses fixings w.r.t. the fractionalities"
#define HEUR_DISPCHAR         'd'
#define HEUR_PRIORITY         -1003000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          3
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

/* event handler properties */
#define EVENTHDLR_NAME         "Nlpdiving"
#define EVENTHDLR_DESC         "bound change event handler for " HEUR_NAME " heuristic"


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXNLPITERABS       200 /**< minimial absolute number of allowed NLP iterations */
#define DEFAULT_MAXNLPITERREL        10 /**< additional allowed number of NLP iterations relative to successfully found solutions */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MINSUCCQUOT         0.1 /**< heuristic will not run if less then this percentage of calls succeeded (0.0: no limit) */
#define DEFAULT_MAXFEASNLPS          10 /**< maximal number of NLPs with feasible solution to solve during one dive */
#define DEFAULT_FIXQUOT             0.2 /**< percentage of fractional variables that should be fixed before the next NLP solve */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_LP                FALSE /**< should the LP relaxation be solved before the NLP relaxation? */
#define DEFAULT_PREFERLPFRACS     FALSE /**< prefer variables that are also fractional in LP solution? */
#define DEFAULT_PREFERCOVER        TRUE /**< should variables in a minimal cover be preferred? */
#define DEFAULT_SOLVESUBMIP       FALSE /**< should a sub-MIP be solved if all cover variables are fixed? */
#define DEFAULT_NLPSTART            's' /**< which point should be used as starting point for the NLP solver? */
#define DEFAULT_VARSELRULE          'd' /**< which variable selection should be used? ('f'ractionality, 'c'oefficient,
                                         *   'p'seudocost, 'g'uided, 'd'ouble)
                                         */
#define DEFAULT_NLPFASTFAIL        TRUE /**< should the NLP solver stop early if it converges slow? */
#define DEFAULT_RANDSEED             97 /**< initial random seed */

#define MINNLPITER                   10 /**< minimal number of NLP iterations allowed in each NLP solving call */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   int                   maxnlpiterabs;      /**< minimial absolute number of allowed NLP iterations */
   int                   maxnlpiterrel;      /**< additional allowed number of NLP iterations relative to successfully found solutions */
   SCIP_Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   int                   maxfeasnlps;        /**< maximal number of NLPs with feasible solution to solve during one dive */
   SCIP_Real             minsuccquot;        /**< heuristic will not run if less then this percentage of calls succeeded (0.0: no limit) */
   SCIP_Real             fixquot;            /**< percentage of fractional variables that should be fixed before the next NLP solve */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Bool             lp;                 /**< should the LP relaxation be solved before the NLP relaxation? */
   SCIP_Bool             preferlpfracs;      /**< prefer variables that are also fractional in LP solution? */
   SCIP_Bool             prefercover;        /**< should variables in a minimal cover be preferred? */
   SCIP_Bool             solvesubmip;        /**< should a sub-MIP be solved if all cover variables are fixed? */
   SCIP_Bool             nlpfastfail;        /**< should the NLP solver stop early if it converges slow? */
   char                  nlpstart;           /**< which point should be used as starting point for the NLP solver? */
   char                  varselrule;         /**< which variable selection should be used? ('f'ractionality, 'c'oefficient,
                                              *   'p'seudocost, 'g'uided, 'd'ouble)
                                              */

   int                   nnlpiterations;     /**< NLP iterations used in this heuristic */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
   int                   nfixedcovervars;    /**< number of variables in the cover that are already fixed */
#ifdef SCIP_STATISTIC
   int                   nnlpsolves;         /**< number of NLP solves */
   int                   nfailcutoff;        /**< number of fails due to cutoff */
   int                   nfaildepth;         /**< number of fails due to too deep */
   int                   nfailnlperror;      /**< number of fails due to NLP error */
#endif
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
};


/*
 * local methods
 */

/** gets fractional variables of last NLP solution along with solution values and fractionalities
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
static
SCIP_RETCODE getNLPFracVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR***           nlpcands,           /**< pointer to store the array of NLP fractional variables, or NULL */
   SCIP_Real**           nlpcandssol,        /**< pointer to store the array of NLP fractional variables solution values, or NULL */
   SCIP_Real**           nlpcandsfrac,       /**< pointer to store the array of NLP fractional variables fractionalities, or NULL */
   int*                  nnlpcands           /**< pointer to store the number of NLP fractional variables , or NULL */
   )
{
   int c;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nlpcands != NULL);
   assert(nlpcandssol != NULL);
   assert(nlpcandsfrac != NULL);
   assert(nnlpcands != NULL);

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetNLPFracVars(scip, nlpcands, nlpcandssol, nlpcandsfrac, nnlpcands, NULL) );

   /* values may be outside the domain in exact arithmetic, but inside the domain within relative tolerance, and still
    * slightly fractional, because SCIPisFeasIntegral() uses absolute tolerance; project value onto domain to avoid this
    * (example: primsol=29.99999853455704, lower bound = 30)
    */
   for( c = 0; c < *nnlpcands; ++c )
   {
      assert(!SCIPisFeasIntegral(scip, (*nlpcandssol)[c]));

      if( (*nlpcandssol)[c] < SCIPvarGetLbLocal((*nlpcands)[c]) || (*nlpcandssol)[c] > SCIPvarGetUbLocal((*nlpcands)[c]) )
      {
         SCIP_Real newval;

         newval = ((*nlpcandssol)[c] < SCIPvarGetLbLocal((*nlpcands)[c]))
            ? SCIPvarGetLbLocal((*nlpcands)[c]) - 0.5*SCIPfeastol(scip)
            : SCIPvarGetUbLocal((*nlpcands)[c]) + 0.5*SCIPfeastol(scip);

         assert(SCIPisFeasIntegral(scip, newval));

         SCIP_CALL( SCIPsetSolVal(scip, heurdata->sol, (*nlpcands)[c], newval) );

         (*nnlpcands)--;

         if( c < *nnlpcands )
         {
            (*nlpcands)[c] = (*nlpcands)[*nnlpcands];
            (*nlpcandssol)[c] = (*nlpcandssol)[*nnlpcands];
            (*nlpcandsfrac)[c] = (*nlpcandsfrac)[*nnlpcands];
         }
      }
   }

   /* prefer decisions on variables which are also fractional in LP solution */
   if( heurdata->preferlpfracs && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      for( c = 0; c < *nnlpcands; ++c )
      {
         if( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, (*nlpcands)[c])) )
            (*nlpcandsfrac)[c] *= 100.0;
      }
   }

   return SCIP_OKAY;
}

/** finds best candidate variable w.r.t. fractionality:
 * - prefer variables that may not be rounded without destroying NLP feasibility:
 *   - of these variables, round least fractional variable in corresponding direction
 * - if all remaining fractional variables may be rounded without destroying NLP feasibility:
 *   - round variable with least increasing objective value
 * - binary variables are prefered
 * - variables in a minimal cover or variables that are also fractional in an optimal LP solution might
 *   also be prefered if a correpsonding parameter is set
 */
static
SCIP_RETCODE chooseFracVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            nlpcands,           /**< array of NLP fractional variables */
   SCIP_Real*            nlpcandssol,        /**< array of NLP fractional variables solution values */
   SCIP_Real*            nlpcandsfrac,       /**< array of NLP fractional variables fractionalities */
   int                   nnlpcands,          /**< number of NLP fractional variables */
   SCIP_HASHMAP*         varincover,         /**< hash map for variables */
   SCIP_Bool             covercomputed,      /**< has a minimal cover been computed? */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround,   /**< pointer to store whether best candidate is trivially roundable */
   SCIP_Bool*            bestcandroundup     /**< pointer to store whether best candidate should be rounded up */
   )
{
   SCIP_Real bestobjgain;
   SCIP_Real bestfrac;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nlpcands != NULL);
   assert(nlpcandssol != NULL);
   assert(nlpcandsfrac != NULL);
   assert(covercomputed == (varincover != NULL));
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);
   assert(bestcandroundup != NULL);

   bestcandmayrounddown = TRUE;
   bestcandmayroundup = TRUE;
   bestobjgain = SCIPinfinity(scip);
   bestfrac = SCIP_INVALID;

   for( c = 0; c < nnlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Bool roundup;
      SCIP_Real frac;
      SCIP_Real obj;
      SCIP_Real objgain;

      var = nlpcands[c];

      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      frac = nlpcandsfrac[c];
      obj = SCIPvarGetObj(var);

      /* since we are not solving the NLP after each fixing, the old NLP solution might be outside the propagated bounds */
      if( SCIPisLT(scip, nlpcandssol[c], SCIPvarGetLbLocal(var)) || SCIPisGT(scip, nlpcandssol[c], SCIPvarGetUbLocal(var)) )
         continue;

      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( bestcandmayrounddown || bestcandmayroundup )
         {
            /* choose rounding direction:
             * - if variable may be rounded in both directions, round corresponding to the fractionality
             * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
             *   the current fractional solution
             */
            if( mayrounddown && mayroundup )
            {
               if( SCIPisEQ(scip, frac, 0.5) )
                  roundup = (SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0);
               else
                  roundup = (frac > 0.5);
            }
            else
               roundup = mayrounddown;

            if( roundup )
            {
               frac = 1.0 - frac;
               objgain = frac*obj;
            }
            else
               objgain = -frac*obj;

            /* penalize too small fractions */
            if( SCIPisEQ(scip, frac, 0.01) )
            {
               /* try to avoid variability; decide randomly if the LP solution can contain some noise.
                * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for increasing the fractionality, i.e., the score.
                */
               if( SCIPrandomGetInt(heurdata->randnumgen, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
                  objgain *= 1000.0;
            }
            else if( frac < 0.01 )
               objgain *= 1000.0;

            /* prefer decisions on binary variables */
            if( !SCIPvarIsBinary(var) )
               objgain *= 1000.0;

            /* prefer decisions on cover variables */
            if( covercomputed && heurdata->prefercover && !SCIPhashmapExists(varincover, var) )
               objgain *= 1000.0;

            /* check, if candidate is new best candidate */
            if( SCIPisLT(scip, objgain, bestobjgain) || (SCIPisEQ(scip, objgain, bestobjgain) && frac < bestfrac) )
            {
               *bestcand = c;
               bestobjgain = objgain;
               bestfrac = frac;
               bestcandmayrounddown = mayrounddown;
               bestcandmayroundup = mayroundup;
               *bestcandroundup = roundup;
            }
         }
      }
      else
      {
         /* the candidate may not be rounded */
         if( SCIPisEQ(scip, frac, 0.5) )
            roundup = (SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0);
         else if( frac < 0.5 )
            roundup = FALSE;
         else
            roundup = TRUE;

         /* adjust fractional part */
         if( roundup )
            frac = 1.0 - frac;

         /* penalize too small fractions */
         if( SCIPisEQ(scip, frac, 0.01) )
         {
            /* try to avoid variability; decide randomly if the LP solution can contain some noise.
             * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for increasing the fractionality, i.e., the score.
             */
            if( SCIPrandomGetInt(heurdata->randnumgen, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
               frac += 10.0;
         }
         else if( frac < 0.01 )
            frac += 10.0;

         /* prefer decisions on binary variables */
         if( !SCIPvarIsBinary(var) )
            frac *= 1000.0;

         /* prefer decisions on cover variables */
         if( covercomputed && heurdata->prefercover && !SCIPhashmapExists(varincover, var) )
            frac *= 1000.0;

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         if( bestcandmayrounddown || bestcandmayroundup || frac < bestfrac )
         {
            *bestcand = c;
            bestfrac = frac;
            bestcandmayrounddown = FALSE;
            bestcandmayroundup = FALSE;
            *bestcandroundup = roundup;
         }
         assert(bestfrac < SCIP_INVALID);
      }
   }

   *bestcandmayround = bestcandmayroundup || bestcandmayrounddown;

   return SCIP_OKAY;
}

/** finds best candidate variable w.r.t. vector length:
 * - round variable with a small ratio between the increase in the objective and the locking numbers
 * - binary variables are prefered
 * - variables in a minimal cover or variables that are also fractional in an optimal LP solution might
 *   also be prefered if a corresponding parameter is set
 */
static
SCIP_RETCODE chooseVeclenVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            nlpcands,           /**< array of NLP fractional variables */
   SCIP_Real*            nlpcandssol,        /**< array of NLP fractional variables solution values */
   SCIP_Real*            nlpcandsfrac,       /**< array of NLP fractional variables fractionalities */
   int                   nnlpcands,          /**< number of NLP fractional variables */
   SCIP_HASHMAP*         varincover,         /**< hash map for variables */
   SCIP_Bool             covercomputed,      /**< has a minimal cover been computed? */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround,   /**< pointer to store whether best candidate is trivially roundable */
   SCIP_Bool*            bestcandroundup     /**< pointer to store whether best candidate should be rounded up */
   )
{
   SCIP_Real bestscore;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nlpcands != NULL);
   assert(nlpcandsfrac != NULL);
   assert(nlpcandssol != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);
   assert(bestcandroundup != NULL);

   *bestcandmayround = TRUE;
   bestscore = SCIP_REAL_MAX;

   /* get best candidate */
   for( c = 0; c < nnlpcands; ++c )
   {
      SCIP_VAR* var;

      SCIP_Real obj;
      SCIP_Real objdelta;
      SCIP_Real frac;
      SCIP_Real score;

      int nlocks;

      SCIP_Bool roundup;

      var = nlpcands[c];

      /* since we are not solving the NLP after each fixing, the old NLP solution might be outside the propagated bounds */
      if( SCIPisLT(scip, nlpcandssol[c], SCIPvarGetLbLocal(var)) || SCIPisGT(scip, nlpcandssol[c], SCIPvarGetUbLocal(var)) )
         continue;

      frac = nlpcandsfrac[c];
      obj = SCIPvarGetObj(var);
      roundup = (obj >= 0.0);
      objdelta = (roundup ? (1.0-frac)*obj : -frac * obj);
      assert(objdelta >= 0.0);

      /* check whether the variable is roundable */
      *bestcandmayround = *bestcandmayround && (SCIPvarMayRoundDown(var) || SCIPvarMayRoundUp(var));
      nlocks = SCIPvarGetNLocksDown(var) + SCIPvarGetNLocksUp(var);

      /* smaller score is better */
      score = (objdelta + SCIPsumepsilon(scip))/((SCIP_Real)nlocks+1.0);

      /* prefer decisions on binary variables */
      if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
         score *= 1000.0;

      /* prefer decisions on cover variables */
      if( covercomputed && heurdata->prefercover && !SCIPhashmapExists(varincover, var) )
         score *= 1000;

      /* check, if candidate is new best candidate */
      if( score < bestscore )
      {
         *bestcand = c;
         bestscore = score;
         *bestcandroundup = roundup;
      }
   }

   return SCIP_OKAY;
}


/** finds best candidate variable w.r.t. locking numbers:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round variable with least number of locks in corresponding direction
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable with least number of locks in opposite of its feasible rounding direction
 * - binary variables are prefered
 * - variables in a minimal cover or variables that are also fractional in an optimal LP solution might
 *   also be prefered if a correpsonding parameter is set
 */
static
SCIP_RETCODE chooseCoefVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            nlpcands,           /**< array of NLP fractional variables */
   SCIP_Real*            nlpcandssol,        /**< array of NLP fractional variables solution values */
   SCIP_Real*            nlpcandsfrac,       /**< array of NLP fractional variables fractionalities */
   int                   nnlpcands,          /**< number of NLP fractional variables */
   SCIP_HASHMAP*         varincover,         /**< hash map for variables */
   SCIP_Bool             covercomputed,      /**< has a minimal cover been computed? */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround,   /**< pointer to store whether best candidate is trivially roundable */
   SCIP_Bool*            bestcandroundup     /**< pointer to store whether best candidate should be rounded up */
   )
{
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   int bestnviolrows;             /* number of violated rows for best candidate */
   SCIP_Real bestcandfrac;        /* fractionality of best candidate */
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nlpcands != NULL);
   assert(nlpcandsfrac != NULL);
   assert(nlpcandssol != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);
   assert(bestcandroundup != NULL);

   bestcandmayrounddown = TRUE;
   bestcandmayroundup = TRUE;
   bestnviolrows = INT_MAX;
   bestcandfrac = SCIP_INVALID;

   /* get best candidate */
   for( c = 0; c < nnlpcands; ++c )
   {
      SCIP_VAR* var;

      int nlocksdown;
      int nlocksup;
      int nviolrows;

      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Bool roundup;
      SCIP_Real frac;

      var = nlpcands[c];
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      frac = nlpcandsfrac[c];

      /* since we are not solving the NLP after each fixing, the old NLP solution might be outside the propagated bounds */
      if( SCIPisLT(scip, nlpcandssol[c], SCIPvarGetLbLocal(var)) || SCIPisGT(scip, nlpcandssol[c], SCIPvarGetUbLocal(var)) )
         continue;

      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( bestcandmayrounddown || bestcandmayroundup )
         {
            /* choose rounding direction:
             * - if variable may be rounded in both directions, round corresponding to the fractionality
             * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
             *   the current fractional solution
             */
            if( mayrounddown && mayroundup )
            {
               if( SCIPisEQ(scip, frac, 0.5) )
                  roundup = (SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0);
               else
                  roundup = (frac > 0.5);
            }
            else
               roundup = mayrounddown;

            if( roundup )
            {
               frac = 1.0 - frac;
               nviolrows = SCIPvarGetNLocksUp(var);
            }
            else
               nviolrows = SCIPvarGetNLocksDown(var);

            /* penalize too small fractions */
            if( SCIPisEQ(scip, frac, 0.01) )
            {
               /* try to avoid variability; decide randomly if the LP solution can contain some noise.
                * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for increasing the fractionality, i.e., the score.
                */
               if( SCIPrandomGetInt(heurdata->randnumgen, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
                  nviolrows *= 100;
            }
            else if( frac < 0.01 )
               nviolrows *= 100;

            /* prefer decisions on binary variables */
            if( !SCIPvarIsBinary(var) )
               nviolrows *= 1000;

            /* prefer decisions on cover variables */
            if( covercomputed && heurdata->prefercover && !SCIPhashmapExists(varincover, var) )
               nviolrows *= 1000;

            /* check, if candidate is new best candidate */
            assert( (0.0 < frac && frac < 1.0) || SCIPvarIsBinary(var) );
            if( nviolrows + frac < bestnviolrows + bestcandfrac )
            {
               *bestcand = c;
               bestnviolrows = nviolrows;
               bestcandfrac = frac;
               bestcandmayrounddown = mayrounddown;
               bestcandmayroundup = mayroundup;
               *bestcandroundup = roundup;
            }
         }
      }
      else
      {
         /* the candidate may not be rounded */
         nlocksdown = SCIPvarGetNLocksDown(var);
         nlocksup = SCIPvarGetNLocksUp(var);

         roundup = (nlocksdown > nlocksup);
         if( !roundup )
         {
            roundup = (nlocksdown == nlocksup);
            if( SCIPisEQ(scip, frac, 0.5) )
               roundup = (roundup && (SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0));
            else
               roundup = (roundup && frac > 0.5);
         }

         if( roundup )
         {
            nviolrows = nlocksup;
            frac = 1.0 - frac;
         }
         else
            nviolrows = nlocksdown;

         /* penalize too small fractions */
         if( SCIPisEQ(scip, frac, 0.01) )
         {
            /* try to avoid variability; decide randomly if the LP solution can contain some noise.
             * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for increasing the fractionality, i.e., the score.
             */
            if( SCIPrandomGetInt(heurdata->randnumgen, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
               nviolrows *= 100;
         }
         else if( frac < 0.01 )
            nviolrows *= 100;

         /* prefer decisions on binary variables */
         if( !SCIPvarIsBinary(var) )
            nviolrows *= 100;

         /* prefer decisions on cover variables */
         if( covercomputed && heurdata->prefercover && !SCIPhashmapExists(varincover, var) )
            nviolrows *= 1000;

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         assert((0.0 < frac && frac < 1.0) || SCIPvarIsBinary(var));
         if( bestcandmayrounddown || bestcandmayroundup || nviolrows + frac < bestnviolrows + bestcandfrac )
         {
            *bestcand = c;
            bestnviolrows = nviolrows;
            bestcandfrac = frac;
            bestcandmayrounddown = FALSE;
            bestcandmayroundup = FALSE;
            *bestcandroundup = roundup;
         }
         assert(bestcandfrac < SCIP_INVALID);
      }
   }

   *bestcandmayround = bestcandmayroundup || bestcandmayrounddown;

   return SCIP_OKAY;
}

/** calculates the pseudocost score for a given variable w.r.t. a given solution value and a given rounding direction */
static
void calcPscostQuot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             primsol,            /**< primal solution of variable */
   SCIP_Real             frac,               /**< fractionality of variable */
   int                   rounddir,           /**< -1: round down, +1: round up, 0: select due to pseudo cost values */
   SCIP_Real*            pscostquot,         /**< pointer to store pseudo cost quotient */
   SCIP_Bool*            roundup,            /**< pointer to store whether the variable should be rounded up */
   SCIP_Bool             prefvar             /**< should this variable be preferred because it is in a minimal cover? */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   assert(heurdata != NULL);
   assert(pscostquot != NULL);
   assert(roundup != NULL);
   assert(SCIPisEQ(scip, frac, primsol - SCIPfeasFloor(scip, primsol)));

   /* bound fractions to not prefer variables that are nearly integral */
   frac = MAX(frac, 0.1);
   frac = MIN(frac, 0.9);

   /* get pseudo cost quotient */
   pscostdown = SCIPgetVarPseudocostVal(scip, var, 0.0-frac);
   pscostup = SCIPgetVarPseudocostVal(scip, var, 1.0-frac);
   assert(pscostdown >= 0.0 && pscostup >= 0.0);

   /* choose rounding direction
    *
    * to avoid performance variability caused by numerics we use random numbers to decide whether we want to roundup or
    * round down if the values to compare are equal within tolerances.
    */
   if( rounddir == -1 )
      *roundup = FALSE;
   else if( rounddir == +1 )
      *roundup = TRUE;
   else if( SCIPisLT(scip, primsol, SCIPvarGetRootSol(var) - 0.4)
         || (SCIPisEQ(scip, primsol, SCIPvarGetRootSol(var) - 0.4) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
      *roundup = FALSE;
   else if( SCIPisGT(scip, primsol, SCIPvarGetRootSol(var) + 0.4)
         || (SCIPisEQ(scip, primsol, SCIPvarGetRootSol(var) + 0.4) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
      *roundup = TRUE;
   else if( SCIPisLT(scip, frac, 0.3) || (SCIPisEQ(scip, frac, 0.3) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
      *roundup = FALSE;
   else if( SCIPisGT(scip, frac, 0.7) || (SCIPisEQ(scip, frac, 0.7) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
      *roundup = TRUE;
   else if( SCIPisLT(scip, pscostdown, pscostup)
         || (SCIPisEQ(scip, pscostdown, pscostup) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0))
      *roundup = FALSE;
   else
      *roundup = TRUE;

   /* calculate pseudo cost quotient */
   if( *roundup )
      *pscostquot = sqrt(frac) * (1.0+pscostdown) / (1.0+pscostup);
   else
      *pscostquot = sqrt(1.0-frac) * (1.0+pscostup) / (1.0+pscostdown);

   /* prefer decisions on binary variables */
   if( SCIPvarIsBinary(var) )
      (*pscostquot) *= 1000.0;

   /* prefer decisions on cover variables */
   if( prefvar )
      (*pscostquot) *= 1000.0;
}

/** finds best candidate variable w.r.t. pseudo costs:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round variable with largest rel. difference of pseudo cost values in corresponding
 *     direction
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable in the objective value direction
 * - binary variables are prefered
 * - variables in a minimal cover or variables that are also fractional in an optimal LP solution might
 *   also be prefered if a correpsonding parameter is set
 */
static
SCIP_RETCODE choosePscostVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            nlpcands,           /**< array of NLP fractional variables */
   SCIP_Real*            nlpcandssol,        /**< array of NLP fractional variables solution values */
   SCIP_Real*            nlpcandsfrac,       /**< array of NLP fractional variables fractionalities */
   int                   nnlpcands,          /**< number of NLP fractional variables */
   SCIP_HASHMAP*         varincover,         /**< hash map for variables */
   SCIP_Bool             covercomputed,      /**< has a minimal cover been computed? */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround,   /**< pointer to store whether best candidate is trivially roundable */
   SCIP_Bool*            bestcandroundup     /**< pointer to store whether best candidate should be rounded up */
   )
{
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   SCIP_Real bestpscostquot;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nlpcands != NULL);
   assert(nlpcandsfrac != NULL);
   assert(nlpcandssol != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);
   assert(bestcandroundup != NULL);

   bestcandmayrounddown = TRUE;
   bestcandmayroundup = TRUE;
   bestpscostquot = -1.0;

   for( c = 0; c < nnlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real primsol;

      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Bool roundup;
      SCIP_Bool prefvar;
      SCIP_Real frac;
      SCIP_Real pscostquot;

      var = nlpcands[c];
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      primsol = nlpcandssol[c];
      frac = nlpcandsfrac[c];
      prefvar = covercomputed && heurdata->prefercover && SCIPhashmapExists(varincover, var);
      pscostquot = SCIP_INVALID;

      if( SCIPisLT(scip, nlpcandssol[c], SCIPvarGetLbLocal(var)) || SCIPisGT(scip, nlpcandssol[c], SCIPvarGetUbLocal(var)) )
         continue;

      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( bestcandmayrounddown || bestcandmayroundup )
         {
            /* choose rounding direction:
             * - if variable may be rounded in both directions, round corresponding to the pseudo cost values
             * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
             *   the current fractional solution
             */
            roundup = FALSE;
            if( mayrounddown && mayroundup )
               calcPscostQuot(scip, heurdata, var, primsol, frac, 0, &pscostquot, &roundup, prefvar);
            else if( mayrounddown )
               calcPscostQuot(scip, heurdata, var, primsol, frac, +1, &pscostquot, &roundup, prefvar);
            else
               calcPscostQuot(scip, heurdata, var, primsol, frac, -1, &pscostquot, &roundup, prefvar);

            assert(!SCIPisInfinity(scip,ABS(pscostquot)));

            /* check, if candidate is new best candidate */
            if( pscostquot > bestpscostquot )
            {
               *bestcand = c;
               bestpscostquot = pscostquot;
               bestcandmayrounddown = mayrounddown;
               bestcandmayroundup = mayroundup;
               *bestcandroundup = roundup;
            }
         }
      }
      else
      {
         /* the candidate may not be rounded: calculate pseudo cost quotient and preferred direction */
         calcPscostQuot(scip, heurdata, var, primsol, frac, 0, &pscostquot, &roundup, prefvar);
         assert(!SCIPisInfinity(scip,ABS(pscostquot)));

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         if( bestcandmayrounddown || bestcandmayroundup || pscostquot > bestpscostquot )
         {
            *bestcand = c;
            bestpscostquot = pscostquot;
            bestcandmayrounddown = FALSE;
            bestcandmayroundup = FALSE;
            *bestcandroundup = roundup;
         }
      }
   }

   *bestcandmayround = bestcandmayroundup || bestcandmayrounddown;

   return SCIP_OKAY;
}

/** finds best candidate variable w.r.t. the incumbent solution:
 * - prefer variables that may not be rounded without destroying LP feasibility:
 *   - of these variables, round a variable to its value in direction of incumbent solution, and choose the
 *     variable that is closest to its rounded value
 * - if all remaining fractional variables may be rounded without destroying LP feasibility:
 *   - round variable in direction that destroys LP feasibility (other direction is checked by SCIProundSol())
 *   - round variable with least increasing objective value
 * - binary variables are prefered
 * - variables in a minimal cover or variables that are also fractional in an optimal LP solution might
 *   also be prefered if a correpsonding parameter is set
 */
static
SCIP_RETCODE chooseGuidedVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            nlpcands,           /**< array of NLP fractional variables */
   SCIP_Real*            nlpcandssol,        /**< array of NLP fractional variables solution values */
   SCIP_Real*            nlpcandsfrac,       /**< array of NLP fractional variables fractionalities */
   int                   nnlpcands,          /**< number of NLP fractional variables */
   SCIP_SOL*             bestsol,            /**< incumbent solution */
   SCIP_HASHMAP*         varincover,         /**< hash map for variables */
   SCIP_Bool             covercomputed,      /**< has a minimal cover been computed? */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Bool*            bestcandmayround,   /**< pointer to store whether best candidate is trivially roundable */
   SCIP_Bool*            bestcandroundup     /**< pointer to store whether best candidate should be rounded up */
   )
{
   SCIP_Real bestobjgain;
   SCIP_Real bestfrac;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nlpcands != NULL);
   assert(nlpcandsfrac != NULL);
   assert(nlpcandssol != NULL);
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);
   assert(bestcandroundup != NULL);

   bestcandmayrounddown = TRUE;
   bestcandmayroundup = TRUE;
   bestobjgain = SCIPinfinity(scip);
   bestfrac = SCIP_INVALID;

   for( c = 0; c < nnlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real bestsolval;
      SCIP_Real solval;
      SCIP_Real obj;
      SCIP_Real frac;
      SCIP_Real objgain;

      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Bool roundup;

      var = nlpcands[c];
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      solval = nlpcandssol[c];
      frac = nlpcandsfrac[c];
      obj = SCIPvarGetObj(var);
      bestsolval = SCIPgetSolVal(scip, bestsol, var);

      /* since we are not solving the NLP after each fixing, the old NLP solution might be outside the propagated bounds */
      if( SCIPisLT(scip, solval, SCIPvarGetLbLocal(var)) || SCIPisGT(scip, solval, SCIPvarGetUbLocal(var)) )
         continue;

      /* select default rounding direction
       * try to avoid variability; decide randomly if the LP solution can contain some noise
       */
      if( SCIPisEQ(scip, solval, bestsolval) )
         roundup = (SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0);
      else
         roundup = (solval < bestsolval);

      if( mayrounddown || mayroundup )
      {
         /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
         if( bestcandmayrounddown || bestcandmayroundup )
         {
            /* choose rounding direction:
             * - if variable may be rounded in both directions, round corresponding to its value in incumbent solution
             * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
             *   the current fractional solution with SCIProundSol()
             */
            if( !mayrounddown || !mayroundup )
               roundup = mayrounddown;

            if( roundup )
            {
               frac = 1.0 - frac;
               objgain = frac*obj;
            }
            else
               objgain = -frac*obj;

            /* penalize too small fractions */
            if( SCIPisEQ(scip, frac, 0.01) )
            {
               /* try to avoid variability; decide randomly if the LP solution can contain some noise.
                * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for increasing the fractionality, i.e., the score.
                */
               if( SCIPrandomGetInt(heurdata->randnumgen, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
                  objgain *= 1000.0;
            }
            else if( frac < 0.01 )
               objgain *= 1000.0;

            /* prefer decisions on binary variables */
            if( !SCIPvarIsBinary(var) )
               objgain *= 1000.0;

            /* prefer decisions on cover variables */
            if( covercomputed && heurdata->prefercover && !SCIPhashmapExists(varincover, var) )
               objgain *= 1000.0;

            /* check, if candidate is new best candidate */
            if( SCIPisLT(scip, objgain, bestobjgain) || (SCIPisEQ(scip, objgain, bestobjgain) && frac < bestfrac) )
            {
               *bestcand = c;
               bestobjgain = objgain;
               bestfrac = frac;
               bestcandmayrounddown = mayrounddown;
               bestcandmayroundup = mayroundup;
               *bestcandroundup = roundup;
            }
         }
      }
      else
      {
         /* the candidate may not be rounded */
         if( roundup )
            frac = 1.0 - frac;

         /* penalize too small fractions */
         if( SCIPisEQ(scip, frac, 0.01) )
         {
            /* try to avoid variability; decide randomly if the LP solution can contain some noise.
             * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for increasing the fractionality, i.e., the score.
             */
            if( SCIPrandomGetInt(heurdata->randnumgen, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
               frac += 10.0;
         }
         else if( frac < 0.01 )
            frac += 10.0;

         /* prefer decisions on binary variables */
         if( !SCIPvarIsBinary(var) )
            frac *= 1000.0;

         /* prefer decisions on cover variables */
         if( covercomputed && heurdata->prefercover && !SCIPhashmapExists(varincover, var) )
            frac *= 1000.0;

         /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
         if( bestcandmayrounddown || bestcandmayroundup || frac < bestfrac )
         {
            *bestcand = c;
            bestfrac = frac;
            bestcandmayrounddown = FALSE;
            bestcandmayroundup = FALSE;
            *bestcandroundup = roundup;
         }
      }
   }

   *bestcandmayround = bestcandmayroundup || bestcandmayrounddown;

   return SCIP_OKAY;
}

/** finds best candidate variable w.r.t. both, the LP and the NLP solution:
 * - choose a variable for which the sum of the distances from the relaxations' solutions to a common
 *   integer value is minimal
 * - binary variables are prefered
 * - variables in a minimal cover might be prefered if a corresponding parameter is set
 */
static
SCIP_RETCODE chooseDoubleVar(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR**            pseudocands,        /**< array of non-fixed variables */
   SCIP_Real*            pseudocandsnlpsol,  /**< array of NLP solution values */
   SCIP_Real*            pseudocandslpsol,   /**< array of LP solution values */
   int                   npseudocands,       /**< number of NLP fractional variables */
   SCIP_HASHMAP*         varincover,         /**< hash map for variables */
   SCIP_Bool             covercomputed,      /**< has a minimal cover been computed? */
   int*                  bestcand,           /**< pointer to store the index of the best candidate variable */
   SCIP_Real*            bestboundval,       /**< pointer to store the bound, the best candidate should be rounded to */
   SCIP_Bool*            bestcandmayround,   /**< pointer to store whether best candidate is trivially roundable */
   SCIP_Bool*            bestcandroundup     /**< pointer to store whether best candidate should be rounded up */
   )
{
   SCIP_Real bestfrac;
   int c;

   /* check preconditions */
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(pseudocands != NULL);
   assert(pseudocandsnlpsol != NULL);
   assert(pseudocandslpsol != NULL);
   assert(covercomputed == (varincover != NULL));
   assert(bestcand != NULL);
   assert(bestcandmayround != NULL);
   assert(bestcandroundup != NULL);

   bestfrac = SCIP_INVALID;

   for( c = 0; c < npseudocands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Bool mayround;
      SCIP_Bool roundup;

      SCIP_Real frac;
      SCIP_Real lpsol;
      SCIP_Real nlpsol;
      SCIP_Real lpsolfloor;
      SCIP_Real nlpsolfloor;
      SCIP_Real lpsolceil;
      SCIP_Real nlpsolceil;
      SCIP_Real boundval;
      SCIP_Real floorval;
      SCIP_Real ceilval;

      var = pseudocands[c];
      lpsol = pseudocandslpsol[c];
      nlpsol = pseudocandsnlpsol[c];

      assert(SCIPvarGetUbLocal(var)-SCIPvarGetLbLocal(var) > 0.5);
      assert(SCIPisLE(scip, SCIPvarGetLbLocal(var), lpsol) && SCIPisLE(scip, lpsol, SCIPvarGetUbLocal(var)));

      /* since we are not solving the NLP after each fixing, the old NLP solution might be outside the propagated bounds */
      if( SCIPisLT(scip, nlpsol, SCIPvarGetLbLocal(var)) || SCIPisGT(scip, nlpsol, SCIPvarGetUbLocal(var)) )
         continue;

      mayround = SCIPvarMayRoundDown(var) || SCIPvarMayRoundUp(var);

      /* if this candidate is trivially roundable, and we already know a candidate that is not, continue */
      if( mayround && !(*bestcandmayround) )
         continue;

      if( SCIPisFeasEQ(scip, lpsol, nlpsol) && SCIPisFeasIntegral(scip, lpsol))
         continue;

      lpsolfloor = SCIPfeasFloor(scip, lpsol);
      nlpsolfloor =  SCIPfeasFloor(scip, nlpsol);
      lpsolceil = SCIPfeasCeil(scip, lpsol);
      nlpsolceil =  SCIPfeasCeil(scip, nlpsol);
      floorval = MIN(lpsolfloor,nlpsolfloor);
      ceilval =  MAX(lpsolceil,nlpsolceil);
      /* if both values are in the same interval, find out which integer is (in sum) the closer one, this will be the
       * new bound. The minima and maxima are necessary since one or both values with be integer
       */
      if( SCIPvarIsBinary(var) || ceilval-floorval < 1.5 )
      {

         frac = 0.33*(lpsol-floorval) + 0.67*(nlpsol-floorval);
         if( frac < 0.5 )
         {
            roundup = FALSE;
            boundval = MIN(lpsolfloor,nlpsolfloor);
         }
         else
         {
            roundup = TRUE;
            frac = 1.0-frac;
            boundval = MAX(nlpsolceil,lpsolceil);
         }
      }
      else
      {
         /* determine new bound in the middle of both relaxations, such that the NLP stays feasible */
         SCIP_Real midval;
         midval = (nlpsol+lpsol)/2.0;
         roundup = nlpsol > lpsol;
         frac = ABS(nlpsol-lpsol);

         if( roundup )
            boundval = SCIPfeasCeil(scip, midval);
         else
            boundval = SCIPfeasFloor(scip, midval);

         assert(roundup == SCIPisGT(scip, nlpsol, boundval));
      }

      /* penalize too small fractions */
      if( SCIPisEQ(scip, frac, 0.01) )
      {
         /* try to avoid variability; decide randomly if the LP solution can contain some noise.
          * use a 1:SCIP_PROBINGSCORE_PENALTYRATIO chance for increasing the fractionality, i.e., the score.
          */
         if( SCIPrandomGetInt(heurdata->randnumgen, 0, SCIP_PROBINGSCORE_PENALTYRATIO) == 0 )
            frac += 10.0;
      }
      else if( frac < 0.01 )
         frac += 10.0;

      /* prefer decisions on binary variables */
      if( !SCIPvarIsBinary(var) )
         frac *= 1000.0;

      /* prefer decisions on cover variables */
      if( covercomputed && heurdata->prefercover && !SCIPhashmapExists(varincover, var) )
         frac *= 1000.0;

      /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
      if( frac < bestfrac || (*bestcandmayround && !mayround) )
      {
         *bestcand = c;
         bestfrac = frac;
         *bestcandmayround = FALSE;
         *bestcandroundup = roundup;
         *bestboundval = boundval;
      }
      assert(bestfrac < SCIP_INVALID);
   }

   if( *bestcandroundup )
      *bestboundval -= 0.5;
   else
      *bestboundval += 0.5;

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_HEUR*            heur,               /**< heuristic structure                                 */
   SCIP_HASHMAP*         varmap,             /**< hash map for variables */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   SCIP_VAR** subvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   int        nvars;                         /* the original problem's number of variables      */
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subvars);
   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** todo setup and solve the subMIP */
static
SCIP_RETCODE doSolveSubMIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< NLP diving subscip */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_VAR**            covervars,          /**< variables in the cover, should be fixed locally */
   int                   ncovervars,         /**< number of variables in the cover */
   SCIP_Bool*            success             /**< pointer to store whether a solution was found */
   )
{
   SCIP_HASHMAP* varmap;
   SCIP_SOL** subsols;
   int c;
   int nsubsols;

   assert(subscip != NULL);
   assert(scip != NULL);
   assert(heur != NULL);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), SCIPgetNVars(scip)) );

   *success = FALSE;

   /* copy original problem to subproblem; do not copy pricers */
   SCIP_CALL( SCIPcopyConsCompression(scip, subscip, varmap, NULL, "undercoversub", NULL, NULL, 0, FALSE, FALSE, TRUE, NULL) );

   /* assert that cover variables are fixed in source and target SCIP */
   for( c = 0; c < ncovervars; c++)
   {
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(covervars[c]), SCIPvarGetUbLocal(covervars[c])));
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal((SCIP_VAR*) SCIPhashmapGetImage(varmap, covervars[c])),
            SCIPvarGetUbGlobal((SCIP_VAR*) SCIPhashmapGetImage(varmap, covervars[c]))));
   }

   /* set parameters for sub-SCIP */

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
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", (SCIP_Longint)100) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", (SCIP_Longint)500) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* enable conflict analysis, disable analysis of boundexceeding LPs, and restrict conflict pool */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetCharParam(subscip, "conflict/useboundlp", 'o') );
   }
   if( !SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }

   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      SCIP_Real cutoffbound;
      SCIP_Real minimprove;

      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      minimprove = 0.01;

      if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
      {
         cutoffbound = (1-minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound(scip) >= 0 )
            cutoffbound = (1 - minimprove)*SCIPgetUpperbound(scip);
         else
            cutoffbound = (1 + minimprove)*SCIPgetUpperbound(scip);
      }
      cutoffbound = MIN(upperbound, cutoffbound);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoffbound) );
   }

   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   for( c = 0; c < nsubsols && !(*success); ++c )
   {
      SCIP_CALL( createNewSol(scip, subscip, heur, varmap, subsols[c], success) );
   }

   SCIPhashmapFree(&varmap);

   return SCIP_OKAY;
}


/** solves subproblem and passes best feasible solution to original SCIP instance */
static
SCIP_RETCODE solveSubMIP(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_VAR**            covervars,          /**< variables in the cover, should be fixed locally */
   int                   ncovervars,         /**< number of variables in the cover */
   SCIP_Bool*            success             /**< pointer to store whether a solution was found */
   )
{
   SCIP* subscip;

   SCIP_RETCODE retcode;


   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, success) );

   if( !(*success) )
      return SCIP_OKAY;

   /* create subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   retcode = doSolveSubMIP(scip, subscip, heur, covervars, ncovervars, success);

   /* free sub-SCIP even if an error occurred during the subscip solve */
   SCIP_CALL( SCIPfree(&subscip) );

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * We update the number of variables fixed in the cover
 */
static
SCIP_DECL_EVENTEXEC(eventExecNlpdiving)
{
   SCIP_EVENTTYPE eventtype;
   SCIP_HEURDATA* heurdata;
   SCIP_VAR* var;

   SCIP_Real oldbound;
   SCIP_Real newbound;
   SCIP_Real otherbound;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);
   assert(0 <= heurdata->nfixedcovervars && heurdata->nfixedcovervars <= SCIPgetNVars(scip));

   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);
   var = SCIPeventGetVar(event);

   eventtype = SCIPeventGetType(event);
   otherbound = (eventtype & SCIP_EVENTTYPE_LBCHANGED) ? SCIPvarGetUbLocal(var) : SCIPvarGetLbLocal(var);

   switch( eventtype )
   {
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
      /* if cover variable is now fixed */
      if( SCIPisFeasEQ(scip, newbound, otherbound) && !SCIPisFeasEQ(scip, oldbound, otherbound) )
      {
         assert(!SCIPisEQ(scip, oldbound, otherbound));
         ++(heurdata->nfixedcovervars);
      }
      break;
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBRELAXED:
      /* if cover variable is now unfixed */
      if( SCIPisFeasEQ(scip, oldbound, otherbound) && !SCIPisFeasEQ(scip, newbound, otherbound) )
      {
         assert(!SCIPisEQ(scip, newbound, otherbound));
         --(heurdata->nfixedcovervars);
      }
      break;
   default:
      SCIPerrorMessage("invalid event type.\n");
      return SCIP_INVALIDDATA;
   }
   assert(0 <= heurdata->nfixedcovervars && heurdata->nfixedcovervars <= SCIPgetNVars(scip));

   /* SCIPdebugMsg(scip, "changed bound of cover variable <%s> from %f to %f (nfixedcovervars: %d).\n", SCIPvarGetName(var),
      oldbound, newbound, heurdata->nfixedcovervars); */

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyNlpdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurNlpdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeNlpdiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitNlpdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, DEFAULT_RANDSEED) );

   /* initialize data */
   heurdata->nnlpiterations = 0;
   heurdata->nsuccess = 0;
   heurdata->nfixedcovervars = 0;
   SCIPstatistic(
      heurdata->nnlpsolves = 0;
      heurdata->nfailcutoff = 0;
      heurdata->nfaildepth = 0;
      heurdata->nfailnlperror = 0;
      );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitNlpdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   SCIPstatistic(
      if( strstr(SCIPgetProbName(scip), "_covering") == NULL && SCIPheurGetNCalls(heur) > 0 )
      {
         SCIPstatisticMessage("%-30s %5" SCIP_LONGINT_FORMAT " sols in %5" SCIP_LONGINT_FORMAT " runs, %6.1fs, %7d NLP iters in %5d NLP solves, %5.1f avg., %3d%% success %3d%% cutoff %3d%% depth %3d%% nlperror\n",
            SCIPgetProbName(scip), SCIPheurGetNSolsFound(heur), SCIPheurGetNCalls(heur), SCIPheurGetTime(heur),
            heurdata->nnlpiterations, heurdata->nnlpsolves, heurdata->nnlpiterations/MAX(1.0,(SCIP_Real)heurdata->nnlpsolves),
            (100*heurdata->nsuccess) / (int)SCIPheurGetNCalls(heur), (100*heurdata->nfailcutoff) / (int)SCIPheurGetNCalls(heur), (100*heurdata->nfaildepth) / (int)SCIPheurGetNCalls(heur), (100*heurdata->nfailnlperror) / (int)SCIPheurGetNCalls(heur)
            );
      }
      );

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolNlpdiving)
{  /*lint --e{715}*/
   SCIP_HEUR* nlpheur;

   if( !SCIPisNLPConstructed(scip) )
      return SCIP_OKAY;

   /* find NLP local search heuristic */
   nlpheur = SCIPfindHeur(scip, "subnlp");

   /* add global linear constraints to NLP relaxation */
   if( nlpheur != NULL )
   {
      SCIP_CALL( SCIPaddLinearConsToNlpHeurSubNlp(scip, nlpheur, TRUE, TRUE) );
   }

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecNlpdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_NLPSOLSTAT nlpsolstat;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_SOL* nlpstartsol;
   SCIP_SOL* bestsol;
   SCIP_VAR** nlpcands;
   SCIP_VAR** covervars;
   SCIP_Real* nlpcandssol;
   SCIP_Real* nlpcandsfrac;
   SCIP_Real* pseudocandslpsol;
   SCIP_Real* pseudocandsnlpsol;
   SCIP_HASHMAP* varincover;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real objval;
   SCIP_Real oldobjval;
   SCIP_Real fixquot;
   SCIP_Real bestboundval;
   SCIP_Real timelim;
   SCIP_Bool bestcandmayround;
   SCIP_Bool bestcandroundup;
   SCIP_Bool nlperror;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool backtracked;
   SCIP_Bool solvenlp;
   SCIP_Bool covercomputed;
   SCIP_Bool solvesubmip;
   SCIP_Longint ncalls;
   SCIP_Longint nsolsfound;
   int avgnnlpiterations;
   int maxnnlpiterations;
   int npseudocands;
   int nlpbranchcands;
   int ncovervars;
   int nnlpcands;
   int startnnlpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int lastnlpsolvedepth;
   int nfeasnlps;
   int bestcand;
   int origiterlim;
   int origfastfail;
   int c;
   int       backtrackdepth;   /* depth where to go when backtracking */
   SCIP_VAR* backtrackvar;     /* (first) variable to fix differently in backtracking */
   SCIP_Real backtrackvarval;  /* (fractional) value of backtrack variable */
   SCIP_Bool backtrackroundup; /* whether variable should be rounded up in backtracking */

   backtrackdepth = -1;
   backtrackvar = NULL;
   backtrackvarval = 0.0;
   backtrackroundup = FALSE;
   bestsol = NULL;
   pseudocandsnlpsol = NULL;
   pseudocandslpsol = NULL;
   covervars = NULL;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   /* assert(SCIPhasCurrentNodeLP(scip)); */

   *result = SCIP_DIDNOTRUN;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* only call heuristic, if an NLP relaxation has been constructed */
   if( !SCIPisNLPConstructed(scip) || SCIPgetNNlpis(scip) == 0 )
      return SCIP_OKAY;

   /* only call heuristic, if the current node will not be cutoff, e.g., due to a (integer and NLP-)feasible LP solution */
   if( SCIPisFeasGE(scip, SCIPgetLocalLowerbound(scip), SCIPgetUpperbound(scip)) )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* do not call heuristic, if it barely succeded */
   if( (SCIPheurGetNSolsFound(heur) + 1.0) / (SCIP_Real)(SCIPheurGetNCalls(heur) + 1.0) < heurdata->minsuccquot )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < heurdata->minreldepth*maxdepth || depth > heurdata->maxreldepth*maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of NLP iterations until heuristic is aborted
    * maximal number is maxnlpiterabs plus a success-depending multiplier of maxnlpiterrel
    */
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + heurdata->nsuccess;
   maxnnlpiterations = heurdata->maxnlpiterabs;
   maxnnlpiterations += (int)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxnlpiterrel);

   /* don't try to dive, if we took too many NLP iterations during diving */
   if( heurdata->nnlpiterations >= maxnnlpiterations )
      return SCIP_OKAY;

   /* allow at least a bit more than the so far average number of NLP iterations per dive */
   avgnnlpiterations = (int)(heurdata->nnlpiterations / MAX(ncalls, 1.0));
   maxnnlpiterations = (int)MAX((SCIP_Real) maxnnlpiterations, (SCIP_Real) heurdata->nnlpiterations + 1.2*avgnnlpiterations);

   /* don't try to dive, if there are no unfixed discrete variables */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, NULL, &npseudocands, NULL) );
   if( npseudocands == 0 )
      return SCIP_OKAY;

   /* set time limit for NLP solver */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelim) );
   if( !SCIPisInfinity(scip, timelim) )
      timelim -= SCIPgetSolvingTime(scip);
   /* possibly exit if time is up (need to check here, since the paramter has to be >= 0) */
   if ( timelim <= 0.0 )
      return SCIP_OKAY;
   SCIP_CALL( SCIPsetNLPRealPar(scip, SCIP_NLPPAR_TILIM, timelim) );

   *result = SCIP_DIDNOTFIND;

#ifdef SCIP_DEBUG
   /* SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, 1) ); */
#endif

   /* set iteration limit */
   SCIP_CALL( SCIPgetNLPIntPar(scip, SCIP_NLPPAR_ITLIM, &origiterlim) );
   SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_ITLIM, maxnnlpiterations - heurdata->nnlpiterations) );

   /* set whether NLP solver should fail fast */
   SCIP_CALL( SCIPgetNLPIntPar(scip, SCIP_NLPPAR_FASTFAIL, &origfastfail) );
   SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_FASTFAIL, (int)heurdata->nlpfastfail) );

   /* set starting point to lp solution */
   SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );

   /* solve NLP relaxation, if not solved already */
   nlpsolstat = SCIPgetNLPSolstat(scip);
   if( nlpsolstat > SCIP_NLPSOLSTAT_FEASIBLE )
   {
      SCIP_NLPSTATISTICS* nlpstatistics;

      SCIP_CALL( SCIPsolveNLP(scip) );
      SCIPstatistic( ++heurdata->nnlpsolves );

      /* update iteration count */
      if( SCIPgetNLPTermstat(scip) < SCIP_NLPTERMSTAT_NUMERR )
      {
         SCIP_CALL( SCIPnlpStatisticsCreate(SCIPblkmem(scip), &nlpstatistics) );
         SCIP_CALL( SCIPgetNLPStatistics(scip, nlpstatistics) );
         heurdata->nnlpiterations += SCIPnlpStatisticsGetNIterations(nlpstatistics);
         SCIPnlpStatisticsFree(SCIPblkmem(scip), &nlpstatistics);
      }

      nlpsolstat = SCIPgetNLPSolstat(scip);

      /* give up, if no feasible solution found */
      if( nlpsolstat >= SCIP_NLPSOLSTAT_LOCINFEASIBLE )
      {
         SCIPdebugMsg(scip, "initial NLP infeasible or not solvable --> stop\n");

         SCIPstatistic(
            if( SCIPgetNLPTermstat(scip) < SCIP_NLPTERMSTAT_NUMERR )
               heurdata->nfailcutoff++;
            else
               heurdata->nfailnlperror++;
         )

         /* reset changed NLP parameters */
         SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_ITLIM, origiterlim) );
         SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_FASTFAIL, origfastfail) );

         return SCIP_OKAY;
      }
   }

   /* get NLP solution */
   SCIP_CALL( SCIPlinkNLPSol(scip, heurdata->sol) );

   /* get fractional variables that should be integral */
   SCIP_CALL( getNLPFracVars(scip, heurdata, &nlpcands, &nlpcandssol, &nlpcandsfrac, &nnlpcands) );
   assert(nnlpcands <= npseudocands);

   /* get LP candidates if LP solution is optimal */
   lpsolstat = SCIPgetLPSolstat(scip);
   if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      nlpbranchcands = SCIPgetNLPBranchCands(scip);
   else
      nlpbranchcands = 0;

   /* don't try to dive, if there are no fractional variables */
   if( nnlpcands == 0 )
   {
      SCIP_Bool success;

      /* check, if solution was feasible and good enough
       *
       * Note that even if the NLP solver found a feasible solution it does not mean that is satisfy the integrality
       * conditions for fixed variables. This happens because the NLP solver uses relative tolerances for the bound
       * constraints but SCIP uses absolute tolerances for checking the integrality conditions.
       */
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, TRUE, TRUE, FALSE, TRUE, TRUE, &success) );
#else
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, TRUE, TRUE, &success) );
#endif
      if( success )
      {
         SCIPdebugMsg(scip, " -> solution of first NLP was integral, feasible, and good enough\n");
         *result = SCIP_FOUNDSOL;
      }

      /* reset changed NLP parameters */
      SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_ITLIM, origiterlim) );
      SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_FASTFAIL, origfastfail) );

      return SCIP_OKAY;
   }

   /* for guided diving: don't dive, if no feasible solutions exist */
   if( heurdata->varselrule == 'g' && SCIPgetNSols(scip) == 0 )
      return SCIP_OKAY;

   /* for guided diving: get best solution that should guide the search; if this solution lives in the original variable space,
    * we cannot use it since it might violate the global bounds of the current problem
    */
   if( heurdata->varselrule == 'g' && SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
      return SCIP_OKAY;

   nlpstartsol = NULL;
   assert(nlpcandsfrac != NULL);
   assert(nlpcands != NULL);
   assert(nlpcandssol != NULL);

   /* save solution of first NLP, if we may use it later */
   if( heurdata->nlpstart != 'n' )
   {
      SCIP_CALL( SCIPcreateNLPSol(scip, &nlpstartsol, heur) );
      SCIP_CALL( SCIPunlinkSol(scip, nlpstartsol) );
   }

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      if( heurdata->maxdiveubquotnosol > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquotnosol * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquotnosol > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   else
   {
      if( heurdata->maxdiveubquot > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquot > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   searchbound = MIN(searchubbound, searchavgbound);
   if( SCIPisObjIntegral(scip) )
      searchbound = SCIPceil(scip, searchbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;

   covercomputed = FALSE;
   varincover = NULL;

   /* compute cover, if required */
   if( heurdata->prefercover || heurdata->solvesubmip )
   {
      SCIP_Real timelimit;
      SCIP_Real memorylimit;

      /* get limits */
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
      if( !SCIPisInfinity(scip, timelimit) )
         timelimit -= SCIPgetSolvingTime(scip);

      /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
      if( !SCIPisInfinity(scip, memorylimit) )
      {
         memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
         memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
      }

      /* compute cover */
      ncovervars = -1;
      SCIP_CALL( SCIPallocBufferArray(scip, &covervars, SCIPgetNVars(scip)) );
      if( memorylimit > 2.0*SCIPgetMemExternEstim(scip)/1048576.0 && timelimit > 0.05 )
      {
         SCIP_CALL( SCIPcomputeCoverUndercover(scip, &ncovervars, covervars, timelimit, memorylimit, SCIPinfinity(scip), FALSE, FALSE, FALSE, 'u', &covercomputed) );
      }

      if( covercomputed )
      {
         /* a cover can be empty, if the cover computation reveals that all nonlinear constraints are linear w.r.t. current variable fixations */
         assert(ncovervars >= 0);

         /* create hash map */
         SCIP_CALL( SCIPhashmapCreate(&varincover, SCIPblkmem(scip), ncovervars) );

         /* process variables in the cover */
         for( c = 0; c < ncovervars; c++ )
         {
            /* insert variable into hash map */
            if( SCIPvarGetType(covervars[c]) < SCIP_VARTYPE_IMPLINT )
            {
               assert(!SCIPhashmapExists(varincover, covervars[c]));
               SCIP_CALL( SCIPhashmapInsert(varincover, covervars[c], (void*) (size_t) (c+1)) );
            }

            /* catch bound change events of cover variables */
            assert(heurdata->eventhdlr != NULL);
            SCIP_CALL( SCIPcatchVarEvent(scip, covervars[c], SCIP_EVENTTYPE_BOUNDCHANGED, heurdata->eventhdlr,
                  (SCIP_EVENTDATA*) heurdata, NULL) );
            assert(!SCIPisFeasEQ(scip, SCIPvarGetLbLocal(covervars[c]), SCIPvarGetUbLocal(covervars[c])));
         }
      }
   }
   else
   {
      covervars = NULL;
      ncovervars = 0;
   }

   /* start diving */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   /* get NLP objective value*/
   objval = SCIPgetNLPObjval(scip);

   SCIPdebugMsg(scip, "(node %" SCIP_LONGINT_FORMAT ") executing nlpdiving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nnlpcands, SCIPgetDualbound(scip), SCIPretransformObj(scip, searchbound));

   /* store a copy of the best solution, if guided diving should be used */
   if( heurdata->varselrule == 'g' )
   {
      assert(SCIPgetNSols(scip) > 0);
      assert(!SCIPsolIsOriginal(SCIPgetBestSol(scip)));

      SCIP_CALL( SCIPcreateSolCopy(scip, &bestsol, SCIPgetBestSol(scip)) );
   }

   /* if double diving should be used, create arrays to hold to entire LP and NLP solution */
   if( heurdata->varselrule == 'd' )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pseudocandslpsol, npseudocands) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pseudocandsnlpsol, npseudocands) );
   }

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   nlperror = FALSE;
   lperror = FALSE;
   cutoff = FALSE;
   divedepth = 0;
   lastnlpsolvedepth = 0;
   backtracked = FALSE;    /* whether we are in backtracking */
   fixquot = heurdata->fixquot;
   nfeasnlps = 1;
   startnnlpcands = nnlpcands;
   solvesubmip = heurdata->solvesubmip;

   while( !nlperror && !cutoff && (nlpsolstat <= SCIP_NLPSOLSTAT_FEASIBLE || nlpsolstat == SCIP_NLPSOLSTAT_UNKNOWN) && nnlpcands > 0
      && (nfeasnlps < heurdata->maxfeasnlps
         || nnlpcands <= startnnlpcands - divedepth/2
         || (nfeasnlps < maxdivedepth && heurdata->nnlpiterations < maxnnlpiterations && objval < searchbound))
      && !SCIPisStopped(scip) )
   {
      SCIP_VAR* var;
      SCIP_Bool updatepscost;

      /* open a new probing node if this will not exceed the maximal tree depth, otherwise stop here */
      if( SCIPgetDepth(scip) < SCIP_MAXTREEDEPTH )
      {
         SCIP_CALL( SCIPnewProbingNode(scip) );
         divedepth++;
      }
      else
         break;

      bestcand = -1;
      bestcandmayround = TRUE;
      bestcandroundup = FALSE;
      bestboundval = SCIP_INVALID;
      updatepscost = TRUE;
      var = NULL;

      /* find best candidate variable */
      switch( heurdata->varselrule )
      {
      case 'c':
         SCIP_CALL( chooseCoefVar(scip, heurdata, nlpcands, nlpcandssol, nlpcandsfrac, nnlpcands, varincover, covercomputed,
               &bestcand, &bestcandmayround, &bestcandroundup) );
         if( bestcand >= 0 )
         {
            var = nlpcands[bestcand];
            bestboundval = nlpcandssol[bestcand];
         }
         break;
      case 'v':
         SCIP_CALL( chooseVeclenVar(scip, heurdata, nlpcands, nlpcandssol, nlpcandsfrac, nnlpcands, varincover, covercomputed,
               &bestcand, &bestcandmayround, &bestcandroundup) );
         if( bestcand >= 0 )
         {
            var = nlpcands[bestcand];
            bestboundval = nlpcandssol[bestcand];
         }
         break;
      case 'p':
         SCIP_CALL( choosePscostVar(scip, heurdata, nlpcands, nlpcandssol, nlpcandsfrac, nnlpcands, varincover, covercomputed,
               &bestcand, &bestcandmayround, &bestcandroundup) );
         if( bestcand >= 0 )
         {
            var = nlpcands[bestcand];
            bestboundval = nlpcandssol[bestcand];
         }
         break;
      case 'g':
         SCIP_CALL( chooseGuidedVar(scip, heurdata, nlpcands, nlpcandssol, nlpcandsfrac, nnlpcands, bestsol, varincover, covercomputed,
               &bestcand, &bestcandmayround, &bestcandroundup) );
         if( bestcand >= 0 )
         {
            var = nlpcands[bestcand];
            bestboundval = nlpcandssol[bestcand];
         }
         break;
      case 'd':
         /* double diving only works if we have both relaxations at hand, otherwise we fall back to fractional diving */
         if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIP_VAR** pseudocands;

            SCIP_CALL( SCIPgetPseudoBranchCands(scip, &pseudocands, &npseudocands, NULL) );
            assert(backtrackdepth > 0 || nnlpcands <= npseudocands);
            assert(SCIPgetNLPBranchCands(scip) <= npseudocands);
            SCIP_CALL( SCIPgetSolVals(scip, NULL, npseudocands, pseudocands, pseudocandslpsol) );
            SCIP_CALL( SCIPgetSolVals(scip, heurdata->sol, npseudocands, pseudocands, pseudocandsnlpsol) );
            SCIP_CALL( chooseDoubleVar(scip, heurdata, pseudocands, pseudocandsnlpsol, pseudocandslpsol, npseudocands,
                  varincover, covercomputed, &bestcand, &bestboundval, &bestcandmayround, &bestcandroundup) );
            if( bestcand >= 0 )
               var = pseudocands[bestcand];
            break;
         }
         else
            updatepscost = FALSE;
         /*lint -fallthrough*/
      case 'f':
         SCIP_CALL( chooseFracVar(scip, heurdata, nlpcands, nlpcandssol, nlpcandsfrac, nnlpcands, varincover, covercomputed,
               &bestcand, &bestcandmayround, &bestcandroundup) );
         if( bestcand >= 0 )
         {
            var = nlpcands[bestcand];
            bestboundval = nlpcandssol[bestcand];
         }
         break;
      default:
         SCIPerrorMessage("invalid variable selection rule\n");
         return SCIP_INVALIDDATA;
      }

      /* if all candidates are roundable, try to round the solution
       * if var == NULL (i.e., bestcand == -1), then all solution candidates are outside bounds
       *   this should only happen if they are slightly outside bounds (i.e., still within feastol, relative tolerance),
       *   but far enough out to be considered as fractional (within feastol, but using absolute tolerance)
       *   in this case, we also try our luck with rounding
       */
      if( (var == NULL || bestcandmayround) && backtrackdepth == -1 )
      {
         SCIP_Bool success;

         /* create solution from diving NLP and try to round it */
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if( success )
         {
            SCIPdebugMsg(scip, "nlpdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

            /* try to add solution to SCIP */
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPtrySol(scip, heurdata->sol, TRUE, TRUE, FALSE, FALSE, TRUE, &success) );
#else
            SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );
#endif

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMsg(scip, " -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      /* if all variables have been found to be essentially integral (even though there is some numerical doubt, see comment above), then stop */
      if( var == NULL )
         break;

      do
      {
         SCIP_Real frac;
         frac = SCIP_INVALID;

         if( backtracked && backtrackdepth > 0 )
         {
            assert(backtrackvar != NULL);

            /* if the variable is already fixed or if the solution value is outside the domain, numerical troubles may have
             * occured or variable was fixed by propagation while backtracking => Abort diving!
             */
            if( SCIPvarGetLbLocal(backtrackvar) >= SCIPvarGetUbLocal(backtrackvar) - 0.5 )
            {
               SCIPdebugMsg(scip, "Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
                  SCIPvarGetName(backtrackvar), SCIPvarGetLbLocal(backtrackvar), SCIPvarGetUbLocal(backtrackvar), backtrackvarval);
               cutoff = TRUE;
               break;
            }
            if( SCIPisFeasLT(scip, backtrackvarval, SCIPvarGetLbLocal(backtrackvar)) || SCIPisFeasGT(scip, backtrackvarval, SCIPvarGetUbLocal(backtrackvar)) )
            {
               SCIPdebugMsg(scip, "selected variable's <%s> solution value is outside the domain [%g,%g] (solval: %.9f), diving aborted\n",
                  SCIPvarGetName(backtrackvar), SCIPvarGetLbLocal(backtrackvar), SCIPvarGetUbLocal(backtrackvar), backtrackvarval);
               assert(backtracked);
               break;
            }

            /* round backtrack variable up or down */
            if( backtrackroundup )
            {
               SCIPdebugMsg(scip, "  dive %d/%d, NLP iter %d/%d: var <%s>, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
                  divedepth, maxdivedepth, heurdata->nnlpiterations, maxnnlpiterations,
                  SCIPvarGetName(backtrackvar), backtrackvarval, SCIPvarGetLbLocal(backtrackvar), SCIPvarGetUbLocal(backtrackvar),
                  SCIPfeasCeil(scip, backtrackvarval), SCIPvarGetUbLocal(backtrackvar));
               SCIP_CALL( SCIPchgVarLbProbing(scip, backtrackvar, SCIPfeasCeil(scip, backtrackvarval)) );
            }
            else
            {
               SCIPdebugMsg(scip, "  dive %d/%d, NLP iter %d/%d: var <%s>, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
                  divedepth, maxdivedepth, heurdata->nnlpiterations, maxnnlpiterations,
                  SCIPvarGetName(backtrackvar), backtrackvarval, SCIPvarGetLbLocal(backtrackvar), SCIPvarGetUbLocal(backtrackvar),
                  SCIPvarGetLbLocal(backtrackvar), SCIPfeasFloor(scip, backtrackvarval));
               SCIP_CALL( SCIPchgVarUbProbing(scip, backtrackvar, SCIPfeasFloor(scip, backtrackvarval)) );
            }

            /* forget about backtrack variable */
            backtrackdepth = -1;

            /* for pseudo cost computation */
            bestcandroundup = backtrackroundup;
            frac = SCIPfrac(scip, backtrackvarval);
            var = backtrackvar;
         }
         else
         {
            assert(var != NULL);

            /* if the variable is already fixed or if the solution value is outside the domain, numerical troubles may have
             * occured or variable was fixed by propagation while backtracking => Abort diving!
             */
            if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
            {
               SCIPdebugMsg(scip, "Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
                  SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), bestboundval);
               cutoff = TRUE;
               break;
            }
            if( SCIPisFeasLT(scip, bestboundval, SCIPvarGetLbLocal(var)) || SCIPisFeasGT(scip, bestboundval, SCIPvarGetUbLocal(var)) )
            {
               SCIPdebugMsg(scip, "selected variable's <%s> solution value is outside the domain [%g,%g] (solval: %.9f), diving aborted\n",
                  SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), bestboundval);
               assert(backtracked);
               break;
            }

            /* apply rounding of best candidate */
            if( bestcandroundup == !backtracked )
            {
               /* round variable up */
               SCIPdebugMsg(scip, "  dive %d/%d, NLP iter %d/%d: var <%s>, round=%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
                  divedepth, maxdivedepth, heurdata->nnlpiterations, maxnnlpiterations,
                  SCIPvarGetName(var), bestcandmayround,
                  bestboundval, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                  SCIPfeasCeil(scip, bestboundval), SCIPvarGetUbLocal(var));
               SCIP_CALL( SCIPchgVarLbProbing(scip, var, SCIPfeasCeil(scip, bestboundval)) );

               /* remember variable for backtracking, if we have none yet (e.g., we are just after NLP solve) or we are half way to the next NLP solve */
               if( backtrackdepth == -1 || (divedepth - lastnlpsolvedepth == (int)(MIN(fixquot * nnlpcands, nlpbranchcands)/2.0)) )
               {
                  backtrackdepth   = divedepth;
                  backtrackvar     = var;
                  backtrackvarval  = bestboundval;
                  backtrackroundup = FALSE;
               }
            }
            else
            {
               /* round variable down */
               SCIPdebugMsg(scip, "  dive %d/%d, NLP iter %d/%d: var <%s>, round=%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
                  divedepth, maxdivedepth, heurdata->nnlpiterations, maxnnlpiterations,
                  SCIPvarGetName(var), bestcandmayround,
                  bestboundval, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
                  SCIPvarGetLbLocal(var), SCIPfeasFloor(scip, bestboundval));
               SCIP_CALL( SCIPchgVarUbProbing(scip, var, SCIPfeasFloor(scip, bestboundval)) );

               /* remember variable for backtracking, if we have none yet (e.g., we are just after NLP solve) or we are half way to the next NLP solve */
               if( backtrackdepth == -1 || (divedepth - lastnlpsolvedepth == (int)(MIN(fixquot * nnlpcands, nlpbranchcands)/2.0)) )
               {
                  backtrackdepth   = divedepth;
                  backtrackvar     = var;
                  backtrackvarval  = bestboundval;
                  backtrackroundup = TRUE;
               }
            }

            /* for pseudo-cost computation */
            if( updatepscost && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               if( heurdata->varselrule == 'd' )
               {
                  assert(pseudocandsnlpsol != NULL);
                  assert(0 <= bestcand && bestcand < npseudocands);
                  frac = SCIPfrac(scip, pseudocandsnlpsol[bestcand]);
               }
               else
                  frac = nlpcandsfrac[bestcand];
            }
         }

         /* apply domain propagation */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         if( cutoff )
         {
            SCIPdebugMsg(scip, "  *** cutoff detected in propagation at level %d\n", SCIPgetProbingDepth(scip));
         }

         /* if all variables in the cover are fixed or there is no fractional variable in the cover,
          * then solve a sub-MIP
          */
         if( !cutoff && solvesubmip && covercomputed &&
            (heurdata->nfixedcovervars == ncovervars ||
               (heurdata->nfixedcovervars >= (ncovervars+1)/2 && !SCIPhashmapExists(varincover, var))) )
         {
            int probingdepth;

            solvesubmip = FALSE;
            probingdepth = SCIPgetProbingDepth(scip);
            assert(probingdepth >= 1);
            assert(covervars != NULL);

            if( heurdata->nfixedcovervars != ncovervars )
            {
               /* fix all remaining cover variables */
               for( c = 0; c < ncovervars && !cutoff ; c++ )
               {
                  SCIP_Real lb;
                  SCIP_Real ub;
                  lb = SCIPvarGetLbLocal(covervars[c]);
                  ub = SCIPvarGetUbLocal(covervars[c]);
                  if( !SCIPisFeasEQ(scip, lb, ub) )
                  {
                     SCIP_Real nlpsolval;

                     /* adopt lpsolval w.r.t. intermediate bound changes by propagation */
                     nlpsolval = SCIPvarGetNLPSol(covervars[c]);
                     nlpsolval = MIN(nlpsolval,ub);
                     nlpsolval = MAX(nlpsolval,lb);
                     assert(SCIPvarGetType(covervars[c]) == SCIP_VARTYPE_CONTINUOUS || SCIPisFeasIntegral(scip, nlpsolval));

                     /* open a new probing node if this will not exceed the maximal tree depth,
                      * otherwise fix all the remaining variables at the same probing node
                      * @todo do we need a new probing node for each fixing? if one of these fixings leads to a cutoff
                      *       we backtrack to the last probing node before we started to fix the covervars (and we do
                      *       not solve the probing LP). thus, it would be less work load in SCIPendProbing
                      *       and SCIPbacktrackProbing.
                      */
                     if( SCIP_MAXTREEDEPTH > SCIPgetDepth(scip) )
                     {
                        SCIP_CALL( SCIPnewProbingNode(scip) );
                     }

                     /* fix and propagate */
                     assert(SCIPisLbBetter(scip, nlpsolval, lb, ub) || SCIPisUbBetter(scip, nlpsolval, lb, ub));

                     if( SCIPisLbBetter(scip, nlpsolval, lb, ub) )
                     {
                        SCIP_CALL( SCIPchgVarLbProbing(scip, covervars[c], nlpsolval) );
                     }
                     if( SCIPisUbBetter(scip, nlpsolval, lb, ub) )
                     {
                        SCIP_CALL( SCIPchgVarUbProbing(scip, covervars[c], nlpsolval) );
                     }

                     SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
                  }
               }
            }

            /* solve sub-MIP or return to standard diving */
            if( cutoff )
            {
               SCIP_CALL( SCIPbacktrackProbing(scip, probingdepth) );
            }
            else
            {
               SCIP_Bool success;
               success = FALSE;

               SCIP_CALL( solveSubMIP(scip, heur, covervars, ncovervars, &success));
               if( success )
                  *result = SCIP_FOUNDSOL;
               backtracked = TRUE; /* to avoid backtracking */
               nnlpcands = 0; /* to force termination */
               cutoff = TRUE;
            }
         }

         /* resolve the diving LP */
         if( !cutoff && !lperror && (heurdata->lp || heurdata->varselrule == 'd')
            && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL && SCIPisLPSolBasic(scip) )
         {
            SCIP_CALL( SCIPsolveProbingLP(scip, 100, &lperror, &cutoff) );

            /* get LP solution status, objective value, and fractional variables, that should be integral */
            lpsolstat = SCIPgetLPSolstat(scip);
            assert(cutoff || (lpsolstat != SCIP_LPSOLSTAT_OBJLIMIT && lpsolstat != SCIP_LPSOLSTAT_INFEASIBLE &&
                  (lpsolstat != SCIP_LPSOLSTAT_OPTIMAL || SCIPisLT(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))));

            if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
            {
               nlpbranchcands = SCIPgetNLPBranchCands(scip);

               /* get new objective value */
               oldobjval = objval;
               objval = SCIPgetLPObjval(scip);

               /* update pseudo cost values */
               if( updatepscost && SCIPisGT(scip, objval, oldobjval) )
               {
                  assert(frac != SCIP_INVALID);  /*lint !e777*/
                  if( bestcandroundup )
                  {
                     SCIP_CALL( SCIPupdateVarPseudocost(scip, var, 1.0-frac, objval - oldobjval, 1.0) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPupdateVarPseudocost(scip, var, 0.0-frac, objval - oldobjval, 1.0) );
                  }
               }
            }
            else
            {
               nlpbranchcands = 0;
            }

            if( cutoff )
            {
               SCIPdebugMsg(scip, "  *** cutoff detected in LP solving at level %d, lpsolstat = %d\n", SCIPgetProbingDepth(scip), lpsolstat);
            }
         }
         else
            lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

         /* check whether we want to solve the NLP, which is the case if
          * - we are in backtracking, or
          * - we have (actively) fixed/rounded fixquot*nnlpcands variables
          * - all fractional variables were rounded/fixed (due to fixing and domain propagation)
          */
         solvenlp = backtracked;
         if( !solvenlp && !cutoff )
         {
            solvenlp = (lastnlpsolvedepth < divedepth - fixquot * nnlpcands);
            if( !solvenlp )
            {
               /* check if fractional NLP variables are left (some may have been fixed by propagation) */
               for( c = 0; c < nnlpcands; ++c )
               {
                  var = nlpcands[c];
                  if( SCIPisLT(scip, nlpcandssol[c], SCIPvarGetLbLocal(var)) || SCIPisGT(scip, nlpcandssol[c], SCIPvarGetUbLocal(var)) )
                     continue;
                  else
                     break;
               }
               if( c == nnlpcands )
                  solvenlp = TRUE;
            }
         }

         nlpsolstat = SCIP_NLPSOLSTAT_UNKNOWN;

         /* resolve the diving NLP */
         if( !cutoff && solvenlp )
         {
            SCIP_NLPTERMSTAT termstat;
            SCIP_NLPSTATISTICS* nlpstatistics;

            /* set iteration limit, allow at least MINNLPITER many iterations */
            SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_ITLIM, MAX(maxnnlpiterations - heurdata->nnlpiterations, MINNLPITER)) );

            /* set time limit for NLP solver */
            SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelim) );
            if( !SCIPisInfinity(scip, timelim) )
               timelim = MAX(0.0, timelim-SCIPgetSolvingTime(scip));/*lint !e666*/
            SCIP_CALL( SCIPsetNLPRealPar(scip, SCIP_NLPPAR_TILIM, timelim) );

            /* set start solution, if we are in backtracking (previous NLP solve was infeasible) */
            if( heurdata->nlpstart != 'n' && backtracked )
            {
               assert(nlpstartsol != NULL);

               SCIPdebugMsg(scip, "setting NLP initial guess\n");

               SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, nlpstartsol) );
            }

            SCIP_CALL( SCIPsolveNLP(scip) );
            SCIPstatistic( ++heurdata->nnlpsolves );

            termstat = SCIPgetNLPTermstat(scip);
            if( termstat >= SCIP_NLPTERMSTAT_NUMERR )
            {
               if( termstat >= SCIP_NLPTERMSTAT_LICERR )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL,
                     "Error while solving NLP in nlpdiving heuristic; NLP solve terminated with code <%d>\n", termstat);
               }
               nlperror = TRUE;
               break;
            }

            /* update iteration count */
            SCIP_CALL( SCIPnlpStatisticsCreate(SCIPblkmem(scip), &nlpstatistics) );
            SCIP_CALL( SCIPgetNLPStatistics(scip, nlpstatistics) );
            heurdata->nnlpiterations += SCIPnlpStatisticsGetNIterations(nlpstatistics);
            SCIPnlpStatisticsFree(SCIPblkmem(scip), &nlpstatistics);

            /* get NLP solution status, objective value, and fractional variables, that should be integral */
            nlpsolstat = SCIPgetNLPSolstat(scip);
            cutoff = (nlpsolstat > SCIP_NLPSOLSTAT_FEASIBLE);

            if( cutoff )
            {
               SCIPdebugMsg(scip, "  *** cutoff detected in NLP solving at level %d, nlpsolstat: %d\n", SCIPgetProbingDepth(scip), nlpsolstat);
            }
            else
            {
               SCIP_CALL( SCIPlinkNLPSol(scip, heurdata->sol) );

               /* remember that we have solve NLP on this depth successfully */
               lastnlpsolvedepth = divedepth;
               /* forget previous backtrack variable, we will never go back to a depth before the current one */
               backtrackdepth = -1;
               /* store NLP solution for warmstarting, if nlpstart is 'f' */
               if( heurdata->nlpstart == 'f' )
               {
                  assert(nlpstartsol != NULL);

                  /* copy NLP solution values into nlpstartsol, is there a better way to do this???? */
                  SCIP_CALL( SCIPlinkNLPSol(scip, nlpstartsol) );
                  SCIP_CALL( SCIPunlinkSol(scip, nlpstartsol) );
               }
               /* increase counter on number of NLP solves with feasible solution */
               ++nfeasnlps;
            }
         }

         /* perform backtracking if a cutoff was detected */
         if( cutoff && !backtracked && heurdata->backtrack )
         {
            if( backtrackdepth == -1 )
            {
               /* backtrack one step */
               SCIPdebugMsg(scip, "  *** cutoff detected at level %d - backtracking one step\n", SCIPgetProbingDepth(scip));
               SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );

               /* after backtracking there has to be at least one open node without exceeding the maximal tree depth */
               assert(SCIP_MAXTREEDEPTH > SCIPgetDepth(scip));

               SCIP_CALL( SCIPnewProbingNode(scip) );
            }
            else
            {
               /* if we have a stored a depth for backtracking, go there */
               SCIPdebugMsg(scip, "  *** cutoff detected at level %d - backtracking to depth %d\n", SCIPgetProbingDepth(scip), backtrackdepth);
               SCIP_CALL( SCIPbacktrackProbing(scip, backtrackdepth-1) );

               /* after backtracking there has to be at least one open node without exceeding the maximal tree depth */
               assert(SCIP_MAXTREEDEPTH > SCIPgetDepth(scip));

               SCIP_CALL( SCIPnewProbingNode(scip) );
               divedepth = backtrackdepth;

               /* do not update pseudocosts if backtracking by more than one level */
               updatepscost = FALSE;

               /* in case, we are feasible after backtracking, fix less variables at once in continuing diving
                * @todo should we remember the fixquot in heurdata for the next run?
                */
               fixquot *= 0.5;
            }
            /* remember that we are backtracking now */
            backtracked = TRUE;
         }
         else
            backtracked = FALSE;
      }
      while( backtracked );

      if( !nlperror && !cutoff && nlpsolstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         /* get new fractional variables */
         SCIP_CALL( getNLPFracVars(scip, heurdata, &nlpcands, &nlpcandssol, &nlpcandsfrac, &nnlpcands) );
      }
      SCIPdebugMsg(scip, "   -> nlpsolstat=%d, objval=%g/%g, nfrac nlp=%d lp=%d\n", nlpsolstat, objval, searchbound, nnlpcands, nlpbranchcands);
   }

   /*lint --e{774}*/
   SCIPdebugMsg(scip, "NLP nlpdiving ABORT due to ");
   if( nlperror || (nlpsolstat > SCIP_NLPSOLSTAT_LOCINFEASIBLE && nlpsolstat != SCIP_NLPSOLSTAT_UNKNOWN) )
   {
      SCIPdebugMsgPrint(scip, "NLP bad status - nlperror: %ud nlpsolstat: %d \n", nlperror, nlpsolstat);
      SCIPstatistic( heurdata->nfailnlperror++ );
   }
   else if( SCIPisStopped(scip) || cutoff )
   {
      SCIPdebugMsgPrint(scip, "LIMIT hit - stop: %ud cutoff: %ud \n", SCIPisStopped(scip), cutoff);
      SCIPstatistic( heurdata->nfailcutoff++ );
   }
   else if(! (divedepth < 10
         || nnlpcands <= startnnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nnlpiterations < maxnnlpiterations && objval < searchbound) ) )
   {
      SCIPdebugMsgPrint(scip, "TOO DEEP - divedepth: %4d cands halfed: %d ltmaxdepth: %d ltmaxiter: %d bound: %d\n", divedepth,
         (nnlpcands > startnnlpcands - divedepth/2), (divedepth >= maxdivedepth), (heurdata->nnlpiterations >= maxnnlpiterations),
         (objval >= searchbound));
      SCIPstatistic( heurdata->nfaildepth++ );
   }
   else if( nnlpcands == 0 && !nlperror && !cutoff && nlpsolstat <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      SCIPdebugMsgPrint(scip, "SUCCESS\n");
   }
   else
   {
      SCIPdebugMsgPrint(scip, "UNKNOWN, very mysterical reason\n");  /* see also special case var == NULL (bestcand == -1) after choose*Var above */
   }

   /* check if a solution has been found */
   if( nnlpcands == 0 && !nlperror && !cutoff && nlpsolstat <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      SCIP_Bool success;

      /* create solution from diving NLP */
      SCIPdebugMsg(scip, "nlpdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP
       *
       * Note that even if the NLP solver found a feasible solution it does not mean that is satisfy the integrality
       * conditions for fixed variables. This happens because the NLP solver uses relative tolerances for the bound
       * constraints but SCIP uses absolute tolerances for checking the integrality conditions.
       */
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, TRUE, TRUE, FALSE, TRUE, TRUE, &success) );
#else
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, TRUE, TRUE, &success) );
#endif

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMsg(scip, " -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
      else
      {
         SCIPdebugMsg(scip, " -> solution was not accepted\n");
      }
   }

   /* end diving */
   SCIP_CALL( SCIPendProbing(scip) );

   /* free hash map and drop variable bound change events */
   if( covercomputed )
   {
      assert(heurdata->eventhdlr != NULL);
      assert(heurdata->nfixedcovervars >= 0); /* variables might have been globally fixed in propagation */
      assert(varincover != NULL);
      assert(covervars != NULL);

      SCIPhashmapFree(&varincover);

      /* drop bound change events of cover variables */
      for( c = 0; c < ncovervars; c++ )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, covervars[c], SCIP_EVENTTYPE_BOUNDCHANGED, heurdata->eventhdlr, (SCIP_EVENTDATA*)heurdata, -1) );
      }
   }
   else
      assert(varincover == NULL);

   /* free array of cover variables */
   if( heurdata->prefercover || heurdata->solvesubmip )
   {
      assert(covervars != NULL || !covercomputed);
      if( covervars != NULL )
         SCIPfreeBufferArray(scip, &covervars);
   }
   else
      assert(covervars == NULL);

   /* free NLP start solution */
   if( nlpstartsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &nlpstartsol) );
   }

   /* reset changed NLP parameters */
   SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_ITLIM, origiterlim) );
   SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_FASTFAIL, origfastfail) );

   /* free copied best solution */
   if( heurdata->varselrule == 'g' )
   {
      assert(bestsol != NULL);
      SCIP_CALL( SCIPfreeSol(scip, &bestsol) );
   }
   else
      assert(bestsol == NULL);

   /* free arrays of LP and NLP solution */
   if( heurdata->varselrule == 'd' )
   {
      assert(pseudocandsnlpsol != NULL);
      assert(pseudocandsnlpsol != NULL);
      SCIPfreeBufferArray(scip, &pseudocandslpsol);
      SCIPfreeBufferArray(scip, &pseudocandsnlpsol);
   }
   else
   {
      assert(pseudocandsnlpsol == NULL);
      assert(pseudocandsnlpsol == NULL);
   }

   if( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

   SCIPdebugMsg(scip, "nlpdiving heuristic finished\n");

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the nlpdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurNlpdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur = NULL;

   /* create heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecNlpdiving, heurdata) );

   assert(heur != NULL);
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyNlpdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeNlpdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitNlpdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitNlpdiving) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolNlpdiving) );

   /* get event handler for bound change events */
   heurdata->eventhdlr = NULL;
   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &heurdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecNlpdiving, NULL) );
   if ( heurdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* nlpdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxnlpiterabs",
         "minimial absolute number of allowed NLP iterations",
         &heurdata->maxnlpiterabs, FALSE, DEFAULT_MAXNLPITERABS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxnlpiterrel",
         "additional allowed number of NLP iterations relative to successfully found solutions",
         &heurdata->maxnlpiterrel, FALSE, DEFAULT_MAXNLPITERREL, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveubquot, TRUE, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/maxdiveavgquot",
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveavgquot, TRUE, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/maxdiveubquotnosol",
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveubquotnosol, TRUE, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/maxdiveavgquotnosol",
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveavgquotnosol, TRUE, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/" HEUR_NAME "/maxfeasnlps",
         "maximal number of NLPs with feasible solution to solve during one dive",
         &heurdata->maxfeasnlps, FALSE, DEFAULT_MAXFEASNLPS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/backtrack",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/lp",
         "should the LP relaxation be solved before the NLP relaxation?",
         &heurdata->lp, TRUE, DEFAULT_LP, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/preferlpfracs",
         "prefer variables that are also fractional in LP solution?",
         &heurdata->preferlpfracs, TRUE, DEFAULT_PREFERLPFRACS, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/minsuccquot",
         "heuristic will not run if less then this percentage of calls succeeded (0.0: no limit)",
         &heurdata->minsuccquot, FALSE, DEFAULT_MINSUCCQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/fixquot",
         "percentage of fractional variables that should be fixed before the next NLP solve",
         &heurdata->fixquot, FALSE, DEFAULT_FIXQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/prefercover",
         "should variables in a minimal cover be preferred?",
         &heurdata->prefercover, FALSE, DEFAULT_PREFERCOVER, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/solvesubmip",
         "should a sub-MIP be solved if all cover variables are fixed?",
         &heurdata->solvesubmip, FALSE, DEFAULT_SOLVESUBMIP, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/nlpfastfail",
         "should the NLP solver stop early if it converges slow?",
         &heurdata->nlpfastfail, FALSE, DEFAULT_NLPFASTFAIL, NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip,
         "heuristics/" HEUR_NAME "/nlpstart",
         "which point should be used as starting point for the NLP solver? ('n'one, last 'f'easible, from dive's'tart)",
         &heurdata->nlpstart, TRUE, DEFAULT_NLPSTART, "fns", NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip,
         "heuristics/" HEUR_NAME "/varselrule",
         "which variable selection should be used? ('f'ractionality, 'c'oefficient, 'p'seudocost, 'g'uided, 'd'ouble, 'v'eclen)",
         &heurdata->varselrule, FALSE, DEFAULT_VARSELRULE, "fcpgdv", NULL, NULL) );

   return SCIP_OKAY;
}
