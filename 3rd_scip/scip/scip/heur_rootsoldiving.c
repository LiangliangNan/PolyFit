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

/**@file   heur_rootsoldiving.c
 * @brief  LP diving heuristic that changes variable's objective values using root LP solution as guide
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_rootsoldiving.h"

#define HEUR_NAME         "rootsoldiving"
#define HEUR_DESC         "LP diving heuristic that changes variable's objective values using root LP solution as guide"
#define HEUR_DISPCHAR     'S'
#define HEUR_PRIORITY     -1005000
#define HEUR_FREQ         20
#define HEUR_FREQOFS       5
#define HEUR_MAXDEPTH     -1
#define HEUR_TIMING       SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.01 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXSOLS              -1 /**< total number of feasible solutions found up to which heuristic is called
                                              *   (-1: no limit) */
#define DEFAULT_DEPTHFAC            0.5 /**< maximal diving depth: number of binary/integer variables times depthfac */
#define DEFAULT_DEPTHFACNOSOL       2.0 /**< maximal diving depth factor if no feasible solution was found yet */

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */
#define DEFAULT_ALPHA               0.9 /**< soft rounding factor to fade out objective coefficients */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   int                   maxsols;            /**< total number of feasible solutions found up to which heuristic is called
                                              *   (-1: no limit) */
   SCIP_Real             depthfac;           /**< maximal diving depth: number of binary/integer variables times depthfac */
   SCIP_Real             depthfacnosol;      /**< maximal diving depth factor if no feasible solution was found yet */
   SCIP_Real             alpha;              /**< soft rounding factor to fade out objective coefficients */
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this heuristic */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
};


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRootsoldiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRootsoldiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRootsoldiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitRootsoldiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data */
   heurdata->nlpiterations = 0;
   heurdata->nsuccess = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRootsoldiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXEC(heurExecRootsoldiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   SCIP_Real* rootsol;
   SCIP_Real* objchgvals;
   int* softroundings;
   int* intvalrounds;
   int nvars;
   int nbinvars;
   int nintvars;
   int nlpcands;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_Real absstartobjval;
   SCIP_Real objstep;
   SCIP_Real alpha;
   SCIP_Real oldobj;
   SCIP_Real newobj;
   SCIP_Bool lperror;
   SCIP_Bool lpsolchanged;
   SCIP_Longint nsolsfound;
   SCIP_Longint ncalls;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int startnlpcands;
   int ncycles;
   int i;

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

   /* get number of fractional variables, that should be integral */
   nlpcands = SCIPgetNLPBranchCands(scip);

   /* don't try to dive, if there are no fractional variables */
   if( nlpcands == 0 )
      return SCIP_OKAY;

   /* calculate the maximal diving depth */
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   if( SCIPgetNSolsFound(scip) == 0 )
      maxdivedepth = (int)(heurdata->depthfacnosol * nvars);
   else
      maxdivedepth = (int)(heurdata->depthfac * nvars);
   maxdivedepth = MAX(maxdivedepth, 10);

   *result = SCIP_DIDNOTFIND;

   /* get all variables of LP */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* get root solution value of all binary and integer variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &rootsol, nbinvars + nintvars) );
   for( i = 0; i < nbinvars + nintvars; i++ )
      rootsol[i] = SCIPvarGetRootSol(vars[i]);

   /* get current LP objective value, and calculate length of a single step in an objective coefficient */
   absstartobjval = SCIPgetLPObjval(scip);
   absstartobjval = ABS(absstartobjval);
   absstartobjval = MAX(absstartobjval, 1.0);
   objstep = absstartobjval / 10.0;

   /* initialize array storing the preferred soft rounding directions and counting the integral value rounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &softroundings, nbinvars + nintvars) );
   BMSclearMemoryArray(softroundings, nbinvars + nintvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &intvalrounds, nbinvars + nintvars) );
   BMSclearMemoryArray(intvalrounds, nbinvars + nintvars);

   /* allocate temporary memory for buffering objective changes */
   SCIP_CALL( SCIPallocBufferArray(scip, &objchgvals, nbinvars + nintvars) );

   /* start diving */
   SCIP_CALL( SCIPstartDive(scip) );

   SCIPdebugMsg(scip, "(node %" SCIP_LONGINT_FORMAT ") executing rootsoldiving heuristic: depth=%d, %d fractionals, dualbound=%g, maxnlpiterations=%" SCIP_LONGINT_FORMAT ", maxdivedepth=%d, LPobj=%g, objstep=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), maxnlpiterations, maxdivedepth,
      SCIPgetLPObjval(scip), objstep);

   lperror = FALSE;
   divedepth = 0;
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   alpha = heurdata->alpha;
   ncycles = 0;
   lpsolchanged = TRUE;
   startnlpcands = nlpcands;
   while( !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0 && ncycles < 10
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations))
      && !SCIPisStopped(scip) )
   {
      SCIP_Bool success;
      int hardroundingidx;
      int hardroundingdir;
      SCIP_Real hardroundingoldbd;
      SCIP_Real hardroundingnewbd;
      SCIP_Bool boundschanged;

      SCIP_RETCODE retcode;

      /* create solution from diving LP and try to round it */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

      if( success )
      {
         SCIPdebugMsg(scip, "rootsoldiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

         /* try to add solution to SCIP */
         SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

         /* check, if solution was feasible and good enough */
         if( success )
         {
            SCIPdebugMsg(scip, " -> solution was feasible and good enough\n");
            *result = SCIP_FOUNDSOL;
         }
      }

      divedepth++;
      hardroundingidx = -1;
      hardroundingdir = 0;
      hardroundingoldbd = 0.0;
      hardroundingnewbd = 0.0;
      boundschanged = FALSE;

      SCIPdebugMsg(scip, "dive %d/%d, LP iter %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ":\n", divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations);

      /* round solution x* from diving LP:
       *   - x~_j = down(x*_j)    if x*_j is integer or binary variable and x*_j <= root solution_j
       *   - x~_j = up(x*_j)      if x*_j is integer or binary variable and x*_j  > root solution_j
       *   - x~_j = x*_j          if x*_j is continuous variable
       * change objective function in diving LP:
       *   - if x*_j is integral, or j is a continuous variable, set obj'_j = alpha * obj_j
       *   - otherwise, set obj'_j = alpha * obj_j + sign(x*_j - x~_j)
       */
      for( i = 0; i < nbinvars + nintvars; i++ )
      {
         SCIP_VAR* var;
         SCIP_Real solval;

         var = vars[i];
         oldobj = SCIPgetVarObjDive(scip, var);
         newobj = oldobj;

         solval =  SCIPvarGetLPSol(var);
         if( SCIPisFeasIntegral(scip, solval) )
         {
            /* if the variable became integral after a soft rounding, count the rounds; after a while, fix it to its
             * current integral value;
             * otherwise, fade out the objective value
             */
            if( softroundings[i] != 0 && lpsolchanged )
            {
               intvalrounds[i]++;
               if( intvalrounds[i] == 5 && SCIPgetVarLbDive(scip, var) < SCIPgetVarUbDive(scip, var) - 0.5 )
               {
                  /* use exact integral value, if the variable is only integral within numerical tolerances */
                  solval = SCIPfloor(scip, solval+0.5);
                  SCIPdebugMsg(scip, " -> fixing <%s> = %g\n", SCIPvarGetName(var), solval);
                  SCIP_CALL( SCIPchgVarLbDive(scip, var, solval) );
                  SCIP_CALL( SCIPchgVarUbDive(scip, var, solval) );
                  boundschanged = TRUE;
               }
            }
            else
               newobj = alpha * oldobj;
         }
         else if( solval <= rootsol[i] )
         {
            /* if the variable was soft rounded most of the time downwards, round it downwards by changing the bounds;
             * otherwise, apply soft rounding by changing the objective value
             */
            softroundings[i]--;
            if( softroundings[i] <= -10 && hardroundingidx == -1 )
            {
               SCIPdebugMsg(scip, " -> hard rounding <%s>[%g] <= %g\n",
                  SCIPvarGetName(var), solval, SCIPfeasFloor(scip, solval));
               hardroundingidx = i;
               hardroundingdir = -1;
               hardroundingoldbd = SCIPgetVarUbDive(scip, var);
               hardroundingnewbd = SCIPfeasFloor(scip, solval);
               SCIP_CALL( SCIPchgVarUbDive(scip, var, hardroundingnewbd) );
               boundschanged = TRUE;
            }
            else
               newobj = alpha * oldobj + objstep;
         }
         else
         {
            /* if the variable was soft rounded most of the time upwards, round it upwards by changing the bounds;
             * otherwise, apply soft rounding by changing the objective value
             */
            softroundings[i]++;
            if( softroundings[i] >= +10 && hardroundingidx == -1 )
            {
               SCIPdebugMsg(scip, " -> hard rounding <%s>[%g] >= %g\n",
                  SCIPvarGetName(var), solval, SCIPfeasCeil(scip, solval));
               hardroundingidx = i;
               hardroundingdir = +1;
               hardroundingoldbd = SCIPgetVarLbDive(scip, var);
               hardroundingnewbd = SCIPfeasCeil(scip, solval);
               SCIP_CALL( SCIPchgVarLbDive(scip, var, hardroundingnewbd) );
               boundschanged = TRUE;
            }
            else
               newobj = alpha * oldobj - objstep;
         }

         /* remember the objective change */
         objchgvals[i] = newobj;
      }

      /* apply objective changes if there was no bound change */
      if( !boundschanged )
      {
         /* apply cached changes on integer variables */
         for( i = 0; i < nbinvars + nintvars; ++i )
         {
            SCIP_VAR* var;

            var = vars[i];
            SCIPdebugMsg(scip, " -> i=%d  var <%s>, solval=%g, rootsol=%g, oldobj=%g, newobj=%g\n",
               i, SCIPvarGetName(var), SCIPvarGetLPSol(var), rootsol[i], SCIPgetVarObjDive(scip, var), objchgvals[i]);

            SCIP_CALL( SCIPchgVarObjDive(scip, var, objchgvals[i]) );
         }

         /* fade out the objective values of the continuous variables */
         for( i = nbinvars + nintvars; i < nvars; i++ )
         {
            SCIP_VAR* var;

            var = vars[i];
            oldobj = SCIPgetVarObjDive(scip, var);
            newobj = alpha * oldobj;

            SCIPdebugMsg(scip, " -> i=%d  var <%s>, solval=%g, oldobj=%g, newobj=%g\n",
               i, SCIPvarGetName(var), SCIPvarGetLPSol(var), oldobj, newobj);

            SCIP_CALL( SCIPchgVarObjDive(scip, var, newobj) );
         }
      }

   SOLVEAGAIN:
      /* resolve the diving LP */
      nlpiterations = SCIPgetNLPIterations(scip);

      retcode = SCIPsolveDiveLP(scip,  MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror, NULL);
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
         SCIPwarningMessage(scip, "Error while solving LP in Rootsoldiving heuristic; LP solve terminated with code <%d>\n", retcode);
         SCIPwarningMessage(scip, "This does not affect the remaining solution procedure --> continue\n");
      }

      if( lperror )
         break;

      /* update iteration count */
      heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;

      /* if no LP iterations were performed, we stayed at the same solution -> count this cycling */
      lpsolchanged = (SCIPgetNLPIterations(scip) != nlpiterations);
      if( lpsolchanged )
         ncycles = 0;
      else if( !boundschanged ) /* do not count if integral variables have been fixed */
         ncycles++;

      /* get LP solution status and number of fractional variables, that should be integral */
      if( lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE && hardroundingidx != -1 )
      {
         SCIP_VAR* var;

         var = vars[hardroundingidx];

         /* round the hard rounded variable to the opposite direction and resolve the LP */
         if( hardroundingdir == -1 )
         {
            SCIPdebugMsg(scip, " -> opposite hard rounding <%s> >= %g\n", SCIPvarGetName(var), hardroundingnewbd + 1.0);
            SCIP_CALL( SCIPchgVarUbDive(scip, var, hardroundingoldbd) );
            SCIP_CALL( SCIPchgVarLbDive(scip, var, hardroundingnewbd + 1.0) );
         }
         else
         {
            SCIPdebugMsg(scip, " -> opposite hard rounding <%s> <= %g\n", SCIPvarGetName(var), hardroundingnewbd - 1.0);
            SCIP_CALL( SCIPchgVarLbDive(scip, var, hardroundingoldbd) );
            SCIP_CALL( SCIPchgVarUbDive(scip, var, hardroundingnewbd - 1.0) );
         }
         hardroundingidx = -1;
         goto SOLVEAGAIN;
      }
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
         nlpcands = SCIPgetNLPBranchCands(scip);
      SCIPdebugMsg(scip, "   -> lpsolstat=%d, nfrac=%d\n", lpsolstat, nlpcands);
   }

   SCIPdebugMsg(scip, "---> diving finished: lpsolstat = %d, depth %d/%d, LP iter %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT "\n",
      lpsolstat, divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations);

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIPdebugMsg(scip, "rootsoldiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

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

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &objchgvals);
   SCIPfreeBufferArray(scip, &intvalrounds);
   SCIPfreeBufferArray(scip, &softroundings);
   SCIPfreeBufferArray(scip, &rootsol);

   SCIPdebugMsg(scip, "rootsoldiving heuristic finished\n");

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the rootsoldiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRootsoldiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Rootsoldiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRootsoldiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRootsoldiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRootsoldiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRootsoldiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRootsoldiving) );

   /* rootsoldiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/rootsoldiving/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/rootsoldiving/maxsols",
         "total number of feasible solutions found up to which heuristic is called (-1: no limit)",
         &heurdata->maxsols, TRUE, DEFAULT_MAXSOLS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/depthfac",
         "maximal diving depth: number of binary/integer variables times depthfac",
         &heurdata->depthfac, TRUE, DEFAULT_DEPTHFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/depthfacnosol",
         "maximal diving depth factor if no feasible solution was found yet",
         &heurdata->depthfacnosol, TRUE, DEFAULT_DEPTHFACNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/rootsoldiving/alpha",
         "soft rounding factor to fade out objective coefficients",
         &heurdata->alpha, TRUE, DEFAULT_ALPHA, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}

