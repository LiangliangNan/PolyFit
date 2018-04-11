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

/**@file   heur_objpscostdiving.c
 * @brief  LP diving heuristic that changes variable's objective value instead of bounds, using pseudo cost values as guide
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_objpscostdiving.h"


#define HEUR_NAME             "objpscostdiving"
#define HEUR_DESC             "LP diving heuristic that changes variable's objective values instead of bounds, using pseudo costs as guide"
#define HEUR_DISPCHAR         'o'
#define HEUR_PRIORITY         -1004000
#define HEUR_FREQ             20
#define HEUR_FREQOFS          4
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.01 /**< maximal fraction of diving LP iterations compared to total iteration number */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXSOLS              -1 /**< total number of feasible solutions found up to which heuristic is called
                                         *   (-1: no limit) */
#define DEFAULT_DEPTHFAC            0.5 /**< maximal diving depth: number of binary/integer variables times depthfac */
#define DEFAULT_DEPTHFACNOSOL       2.0 /**< maximal diving depth factor if no feasible solution was found yet */
#define DEFAULT_RANDSEED            139 /**< initial random seed */

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to total iteration number */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   int                   maxsols;            /**< total number of feasible solutions found up to which heuristic is called
                                              *   (-1: no limit) */
   SCIP_Real             depthfac;           /**< maximal diving depth: number of binary/integer variables times depthfac */
   SCIP_Real             depthfacnosol;      /**< maximal diving depth factor if no feasible solution was found yet */
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this heuristic */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
};


/*
 * local methods
 */

static
void calcPscostQuot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             primsol,            /**< primal solution of variable */
   SCIP_Real             frac,               /**< fractionality of variable */
   int                   rounddir,           /**< -1: round down, +1: round up, 0: select due to pseudo cost values */
   SCIP_Real*            pscostquot,         /**< pointer to store pseudo cost quotient */
   SCIP_Bool*            roundup             /**< pointer to store whether the variable should be rounded up */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   assert(heurdata != NULL);
   assert(pscostquot != NULL);
   assert(roundup != NULL);

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
   else if( SCIPisLT(scip, frac, 0.3) || (SCIPisEQ(scip, frac, 0.3) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
      *roundup = FALSE;
   else if( SCIPisGT(scip, frac, 0.7) || (SCIPisEQ(scip, frac, 0.7) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
      *roundup = TRUE;
   else if( SCIPisLT(scip, primsol, SCIPvarGetRootSol(var) - 0.4)
         || (SCIPisEQ(scip, primsol, SCIPvarGetRootSol(var) - 0.4) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
      *roundup = FALSE;
   else if( SCIPisGT(scip, primsol, SCIPvarGetRootSol(var) + 0.4)
         || (SCIPisEQ(scip, primsol, SCIPvarGetRootSol(var) + 0.4) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
      *roundup = TRUE;
   else if( SCIPisLT(scip, pscostdown, pscostup)
         || (SCIPisEQ(scip, pscostdown, pscostup) && SCIPrandomGetInt(heurdata->randnumgen, 0, 1) == 0) )
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
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyObjpscostdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurObjpscostdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeObjpscostdiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitObjpscostdiving) /*lint --e{715}*/
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
   heurdata->nlpiterations = 0;
   heurdata->nsuccess = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitObjpscostdiving) /*lint --e{715}*/
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

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecObjpscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_VAR* var;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   SCIP_Real primsol;
   SCIP_Real frac;
   SCIP_Real pscostquot;
   SCIP_Real bestpscostquot;
   SCIP_Real oldobj;
   SCIP_Real newobj;
   SCIP_Real objscale;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   SCIP_Bool bestcandroundup;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Bool roundup;
   SCIP_Bool lperror;
   SCIP_Longint ncalls;
   SCIP_Longint nsolsfound;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   int* roundings;
   int nvars;
   int varidx;
   int nlpcands;
   int startnlpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestcand;
   int c;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only apply heuristic, if only a few solutions have been found */
   if( heurdata->maxsols >= 0 && SCIPgetNSolsFound(scip) >= heurdata->maxsols )
      return SCIP_OKAY;

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < heurdata->minreldepth*maxdepth || depth > heurdata->maxreldepth*maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + heurdata->nsuccess;
   maxnlpiterations = (SCIP_Longint)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations);
   maxnlpiterations += heurdata->maxlpiterofs;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + MINLPITER);

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );

   /* don't try to dive, if there are no fractional variables */
   if( nlpcands == 0 )
      return SCIP_OKAY;

   /* calculate the maximal diving depth */
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   if( SCIPgetNSolsFound(scip) == 0 )
      maxdivedepth = (int)(heurdata->depthfacnosol * nvars);
   else
      maxdivedepth = (int)(heurdata->depthfac * nvars);
   maxdivedepth = MIN(maxdivedepth, 10*maxdepth);


   *result = SCIP_DIDNOTFIND;

   /* get temporary memory for remembering the current soft roundings */
   SCIP_CALL( SCIPallocBufferArray(scip, &roundings, nvars) );
   BMSclearMemoryArray(roundings, nvars);

   /* start diving */
   SCIP_CALL( SCIPstartDive(scip) );

   SCIPdebugMsg(scip, "(node %" SCIP_LONGINT_FORMAT ") executing objpscostdiving heuristic: depth=%d, %d fractionals, dualbound=%g, maxnlpiterations=%" SCIP_LONGINT_FORMAT ", maxdivedepth=%d\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), maxnlpiterations, maxdivedepth);

   /* dive as long we are in the given diving depth and iteration limits and fractional variables exist, but
    * - if the last objective change was in a direction, that corresponds to a feasible rounding, we continue in any case
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   lperror = FALSE;
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   divedepth = 0;
   startnlpcands = nlpcands;
   while( !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && nlpcands <= startnlpcands - divedepth/10
            && heurdata->nlpiterations < maxnlpiterations)) && !SCIPisStopped(scip) )
   {
      SCIP_RETCODE retcode;

      divedepth++;

      /* choose variable for objective change:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, change objective value of variable with largest rel. difference of pseudo cost values
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - change objective value of variable with largest rel. difference of pseudo cost values
       */
      bestcand = -1;
      bestpscostquot = -1.0;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestcandroundup = FALSE;
      for( c = 0; c < nlpcands; ++c )
      {
         var = lpcands[c];
         mayrounddown = SCIPvarMayRoundDown(var);
         mayroundup = SCIPvarMayRoundUp(var);
         primsol = lpcandssol[c];
         frac = lpcandsfrac[c];
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
                  calcPscostQuot(scip, heurdata, var, primsol, frac, 0, &pscostquot, &roundup);
               else if( mayrounddown )
                  calcPscostQuot(scip, heurdata, var, primsol, frac, +1, &pscostquot, &roundup);
               else
                  calcPscostQuot(scip, heurdata, var, primsol, frac, -1, &pscostquot, &roundup);

               /* prefer variables, that have already been soft rounded but failed to get integral */
               varidx = SCIPvarGetProbindex(var);
               assert(0 <= varidx && varidx < nvars);
               if( roundings[varidx] != 0 )
                  pscostquot *= 1000.0;

               /* check, if candidate is new best candidate */
               if( pscostquot > bestpscostquot )
               {
                  bestcand = c;
                  bestpscostquot = pscostquot;
                  bestcandmayrounddown = mayrounddown;
                  bestcandmayroundup = mayroundup;
                  bestcandroundup = roundup;
               }
            }
         }
         else
         {
            /* the candidate may not be rounded: calculate pseudo cost quotient and preferred direction */
            calcPscostQuot(scip, heurdata, var, primsol, frac, 0, &pscostquot, &roundup);

            /* prefer variables, that have already been soft rounded but failed to get integral */
            varidx = SCIPvarGetProbindex(var);
            assert(0 <= varidx && varidx < nvars);
            if( roundings[varidx] != 0 )
               pscostquot *= 1000.0;

            /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
            if( bestcandmayrounddown || bestcandmayroundup || pscostquot > bestpscostquot )
            {
               bestcand = c;
               bestpscostquot = pscostquot;
               bestcandmayrounddown = FALSE;
               bestcandmayroundup = FALSE;
               bestcandroundup = roundup;
            }
         }
      }
      assert(bestcand != -1);

      /* if all candidates are roundable, try to round the solution */
      if( bestcandmayrounddown || bestcandmayroundup )
      {
         SCIP_Bool success;

         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if( success )
         {
            SCIPdebugMsg(scip, "objpscostdiving found roundable primal solution: obj=%g\n",
               SCIPgetSolOrigObj(scip, heurdata->sol));

            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMsg(scip, " -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      var = lpcands[bestcand];

      /* check, if the best candidate was already subject to soft rounding */
      varidx = SCIPvarGetProbindex(var);
      assert(0 <= varidx && varidx < nvars);
      if( roundings[varidx] == +1 )
      {
         /* variable was already soft rounded upwards: hard round it downwards */
         SCIP_CALL( SCIPchgVarUbDive(scip, var, SCIPfeasFloor(scip, lpcandssol[bestcand])) );
         SCIPdebugMsg(scip, "  dive %d/%d: var <%s>, round=%u/%u, sol=%g, was already soft rounded upwards -> bounds=[%g,%g]\n",
            divedepth, maxdivedepth, SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var));
      }
      else if( roundings[varidx] == -1 )
      {
         /* variable was already soft rounded downwards: hard round it upwards */
         SCIP_CALL( SCIPchgVarLbDive(scip, var, SCIPfeasCeil(scip, lpcandssol[bestcand])) );
         SCIPdebugMsg(scip, "  dive %d/%d: var <%s>, round=%u/%u, sol=%g, was already soft rounded downwards -> bounds=[%g,%g]\n",
            divedepth, maxdivedepth, SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var));
      }
      else
      {
         assert(roundings[varidx] == 0);

         /* apply soft rounding of best candidate via a change in the objective value */
         objscale = divedepth * 1000.0;
         oldobj = SCIPgetVarObjDive(scip, var);
         if( bestcandroundup )
         {
            /* soft round variable up: make objective value (more) negative */
            if( oldobj < 0.0 )
               newobj = objscale * oldobj;
            else
               newobj = -objscale * oldobj;
            newobj = MIN(newobj, -objscale);

            /* remember, that this variable was soft rounded upwards */
            roundings[varidx] = +1;
         }
         else
         {
            /* soft round variable down: make objective value (more) positive */
            if( oldobj > 0.0 )
               newobj = objscale * oldobj;
            else
               newobj = -objscale * oldobj;
            newobj = MAX(newobj, objscale);

            /* remember, that this variable was soft rounded downwards */
            roundings[varidx] = -1;
         }
         SCIP_CALL( SCIPchgVarObjDive(scip, var, newobj) );
         SCIPdebugMsg(scip, "  dive %d/%d, LP iter %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ": var <%s>, round=%u/%u, sol=%g, bounds=[%g,%g], obj=%g, newobj=%g\n",
            divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
            SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var), oldobj, newobj);
      }

      /* resolve the diving LP */
      nlpiterations = SCIPgetNLPIterations(scip);
      retcode =  SCIPsolveDiveLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror, NULL);
      lpsolstat = SCIPgetLPSolstat(scip);

      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
      if( retcode != SCIP_OKAY )
      {
#ifndef NDEBUG
         if( lpsolstat != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
         {
            SCIP_CALL( retcode );
         }
#endif
         SCIPwarningMessage(scip, "Error while solving LP in Objpscostdiving heuristic; LP solve terminated with code <%d>\n", retcode);
         SCIPwarningMessage(scip, "This does not affect the remaining solution procedure --> continue\n");
      }

      if( lperror )
         break;

      /* update iteration count */
      heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;

      /* get LP solution status  and fractional variables, that should be integral */
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new fractional variables */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );
      }
      SCIPdebugMsg(scip, "   -> lpsolstat=%d, nfrac=%d\n", lpsolstat, nlpcands);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIPdebugMsg(scip, "objpscostdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMsg(scip, " -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   SCIP_CALL( SCIPendDive(scip) );

   if( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

   /* free temporary memory for remembering the current soft roundings */
   SCIPfreeBufferArray(scip, &roundings);

   SCIPdebugMsg(scip, "objpscostdiving heuristic finished\n");

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the objpscostdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurObjpscostdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Objpscostdiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecObjpscostdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyObjpscostdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeObjpscostdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitObjpscostdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitObjpscostdiving) );

   /* objpscostdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/objpscostdiving/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/objpscostdiving/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/objpscostdiving/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to total iteration number",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/objpscostdiving/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/objpscostdiving/maxsols",
         "total number of feasible solutions found up to which heuristic is called (-1: no limit)",
         &heurdata->maxsols, TRUE, DEFAULT_MAXSOLS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/objpscostdiving/depthfac",
         "maximal diving depth: number of binary/integer variables times depthfac",
         &heurdata->depthfac, TRUE, DEFAULT_DEPTHFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/objpscostdiving/depthfacnosol",
         "maximal diving depth factor if no feasible solution was found yet",
         &heurdata->depthfacnosol, TRUE, DEFAULT_DEPTHFACNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

