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

/**@file   heur_oneopt.c
 * @brief  improvement heuristic that alters single variable values
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_oneopt.h"

/* @note If the heuristic runs in the root node, the timing is changed to (SCIP_HEURTIMING_DURINGLPLOOP |
 *       SCIP_HEURTIMING_BEFORENODE), see SCIP_DECL_HEURINITSOL callback.
 */

#define HEUR_NAME             "oneopt"
#define HEUR_DESC             "1-opt heuristic which tries to improve setting of single integer variables"
#define HEUR_DISPCHAR         'b'
#define HEUR_PRIORITY         -20000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_WEIGHTEDOBJ   TRUE           /**< should the objective be weighted with the potential shifting value when sorting the shifting candidates? */
#define DEFAULT_DURINGROOT    TRUE           /**< should the heuristic be called before and during the root node? */
#define DEFAULT_BEFOREPRESOL  FALSE          /**< should the heuristic be called before presolving */
#define DEFAULT_FORCELPCONSTRUCTION FALSE    /**< should the construction of the LP be forced even if LP solving is deactivated? */
#define DEFAULT_USELOOP       TRUE           /**< should the heuristic continue to run as long as improvements are found? */
/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   lastsolindex;       /**< index of the last solution for which oneopt was performed */
   SCIP_Bool             weightedobj;        /**< should the objective be weighted with the potential shifting value when sorting the shifting candidates? */
   SCIP_Bool             duringroot;         /**< should the heuristic be called before and during the root node? */
   SCIP_Bool             forcelpconstruction;/**< should the construction of the LP be forced even if LP solving is deactivated? */
   SCIP_Bool             beforepresol;       /**< should the heuristic be called before presolving */
   SCIP_Bool             useloop;            /**< should the heuristic continue to run as long as improvements are found? */
};


/*
 * Local methods
 */

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< zeroobj heuristic structure                         */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;                         /* the original problem's number of variables      */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** compute value by which the solution of variable @p var can be shifted */
static
SCIP_Real calcShiftVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable that should be shifted */
   SCIP_Real             solval,             /**< current solution value */
   SCIP_Real*            activities          /**< LP row activities */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_Real shiftval;

   SCIP_COL* col;
   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   SCIP_Bool shiftdown;

   int ncolrows;
   int i;


   /* get variable's solution value, global bounds and objective coefficient */
   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);
   obj = SCIPvarGetObj(var);
   shiftdown = TRUE;

   /* determine shifting direction and maximal possible shifting w.r.t. corresponding bound */
   if( obj > 0.0 && SCIPisFeasGE(scip, solval - 1.0, lb) )
      shiftval = SCIPfeasFloor(scip, solval - lb);
   else if( obj < 0.0 && SCIPisFeasLE(scip, solval + 1.0, ub) )
   {
      shiftval = SCIPfeasFloor(scip, ub - solval);
      shiftdown = FALSE;
   }
   else
      return 0.0;


   SCIPdebugMsg(scip, "Try to shift %s variable <%s> with\n", shiftdown ? "down" : "up", SCIPvarGetName(var) );
   SCIPdebugMsg(scip, "    lb:<%g> <= val:<%g> <= ub:<%g> and obj:<%g> by at most: <%g>\n", lb, solval, ub, obj, shiftval);

   /* get data of LP column */
   col = SCIPvarGetCol(var);
   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   ncolrows = SCIPcolGetNLPNonz(col);

   assert(ncolrows == 0 || (colrows != NULL && colvals != NULL));

   /* find minimal shift value, st. all rows stay valid */
   for( i = 0; i < ncolrows && shiftval > 0.0; ++i )
   {
      SCIP_ROW* row;
      int rowpos;

      row = colrows[i];
      rowpos = SCIProwGetLPPos(row);
      assert(-1 <= rowpos && rowpos < SCIPgetNLPRows(scip) );

      /* only global rows need to be valid */
      if( rowpos >= 0 && !SCIProwIsLocal(row) )
      {
         SCIP_Real shiftvalrow;

         assert(SCIProwIsInLP(row));

         if( shiftdown == (colvals[i] > 0) )
            shiftvalrow = SCIPfeasFloor(scip, (activities[rowpos] - SCIProwGetLhs(row)) / ABS(colvals[i]));
         else
            shiftvalrow = SCIPfeasFloor(scip, (SCIProwGetRhs(row) -  activities[rowpos]) / ABS(colvals[i]));
#ifdef SCIP_DEBUG
         if( shiftvalrow < shiftval )
         {
            SCIPdebugMsg(scip, " -> The shift value had to be reduced to <%g>, because of row <%s>.\n",
               shiftvalrow, SCIProwGetName(row));
            SCIPdebugMsg(scip, "    lhs:<%g> <= act:<%g> <= rhs:<%g>, colval:<%g>\n",
               SCIProwGetLhs(row), activities[rowpos], SCIProwGetRhs(row), colvals[i]);
         }
#endif
         shiftval = MIN(shiftval, shiftvalrow);
         /* shiftvalrow might be negative, if we detected infeasibility -> make sure that shiftval is >= 0 */
         shiftval = MAX(shiftval, 0.0);
      }
   }
   if( shiftdown )
      shiftval *= -1.0;

   /* we must not shift variables to infinity */
   if( SCIPisInfinity(scip, solval + shiftval) )
      shiftval = 0.0;

   return shiftval;
}


/** update row activities after a variable's solution value changed */
static
SCIP_RETCODE updateRowActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            activities,         /**< LP row activities */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             shiftval            /**< value that is added to variable */
   )
{
   SCIP_Real* colvals;
   SCIP_ROW** colrows;
   SCIP_COL* col;

   int ncolrows;
   int i;

   assert(activities != NULL);

   /* get data of column associated to variable */
   col = SCIPvarGetCol(var);
   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   ncolrows = SCIPcolGetNLPNonz(col);
   assert(ncolrows == 0 || (colrows != NULL && colvals != NULL));

   /* enumerate all rows with nonzero entry in this column */
   for( i = 0; i < ncolrows; ++i )
   {
      SCIP_ROW* row;
      int rowpos;

      row = colrows[i];
      rowpos = SCIProwGetLPPos(row);
      assert(-1 <= rowpos && rowpos < SCIPgetNLPRows(scip) );

      /* update row activity, only regard global rows in the LP */
      if( rowpos >= 0 && !SCIProwIsLocal(row) )
      {
         activities[rowpos] +=  shiftval * colvals[i];

         if( SCIPisInfinity(scip, activities[rowpos]) )
            activities[rowpos] = SCIPinfinity(scip);
         else if( SCIPisInfinity(scip, -activities[rowpos]) )
            activities[rowpos] = -SCIPinfinity(scip);
      }
   }

   return SCIP_OKAY;
}

/** setup and solve oneopt sub-SCIP */
static
SCIP_RETCODE setupAndSolveSubscipOneopt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_HEUR*            heur,               /**< mutation heuristic */
   SCIP_VAR**            vars,               /**< SCIP variables */
   SCIP_VAR**            subvars,            /**< subproblem's variables */
   SCIP_SOL*             bestsol,            /**< incumbent solution */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   SCIP_Bool*            valid               /**< pointer to store the valid value */
   )
{
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_SOL** subsols;
   SCIP_SOL* startsol;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   int nsubsols;
   int nvars;                                /* number of original problem's variables          */
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);

   nvars = SCIPgetNVars(scip);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

   /* copy complete SCIP instance */
   *valid = FALSE;
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "oneopt", TRUE, FALSE, TRUE, valid) );
   SCIP_CALL( SCIPtransformProb(subscip) );

   /* get variable image */
   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* copy the solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, bestsol, nvars, vars, subsolvals) );

   /* create start solution for the subproblem */
   SCIP_CALL( SCIPcreateOrigSol(subscip, &startsol, NULL) );
   SCIP_CALL( SCIPsetSolVals(subscip, startsol, nvars, subvars, subsolvals) );

   /* try to add new solution to sub-SCIP and free it immediately */
   *valid = FALSE;
   SCIP_CALL( SCIPtrySolFree(subscip, &startsol, FALSE, FALSE, FALSE, FALSE, FALSE, valid) );
   SCIPfreeBufferArray(scip, &subsolvals);
   SCIPhashmapFree(&varmapfw);

   /* deactivate basically everything except oneopt in the sub-SCIP */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(subscip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* set limits for the subproblem */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", 1LL) );

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

   /* if necessary, some of the parameters have to be unfixed first */
   if( SCIPisParamFixed(subscip, "lp/solvefreq") )
   {
      SCIPwarningMessage(scip, "unfixing parameter lp/solvefreq in subscip of oneopt heuristic\n");
      SCIP_CALL( SCIPunfixParam(subscip, "lp/solvefreq") );
   }
   SCIP_CALL( SCIPsetIntParam(subscip, "lp/solvefreq", -1) );

   if( SCIPisParamFixed(subscip, "heuristics/oneopt/freq") )
   {
      SCIPwarningMessage(scip, "unfixing parameter heuristics/oneopt/freq in subscip of oneopt heuristic\n");
      SCIP_CALL( SCIPunfixParam(subscip, "heuristics/oneopt/freq") );
   }
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/oneopt/freq", 1) );

   if( SCIPisParamFixed(subscip, "heuristics/oneopt/forcelpconstruction") )
   {
      SCIPwarningMessage(scip, "unfixing parameter heuristics/oneopt/forcelpconstruction in subscip of oneopt heuristic\n");
      SCIP_CALL( SCIPunfixParam(subscip, "heuristics/oneopt/forcelpconstruction") );
   }
   SCIP_CALL( SCIPsetBoolParam(subscip, "heuristics/oneopt/forcelpconstruction", TRUE) );

   /* avoid recursive call, which would lead to an endless loop */
   if( SCIPisParamFixed(subscip, "heuristics/oneopt/beforepresol") )
   {
      SCIPwarningMessage(scip, "unfixing parameter heuristics/oneopt/beforepresol in subscip of oneopt heuristic\n");
      SCIP_CALL( SCIPunfixParam(subscip, "heuristics/oneopt/beforepresol") );
   }
   SCIP_CALL( SCIPsetBoolParam(subscip, "heuristics/oneopt/beforepresol", FALSE) );

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   if( *valid )
   {
      /* errors in solving the subproblem should not kill the overall solving process;
       * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
      SCIP_CALL_ABORT( SCIPsolve(subscip) );

#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      *valid = FALSE;
      for( i = 0; i < nsubsols && !(*valid); ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], valid) );
         if( *valid )
            *result = SCIP_FOUNDSOL;
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyOneopt)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurOneopt(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeOneopt)
{   /*lint --e{715}*/
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


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolOneopt)
{
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* create heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* if the heuristic is called at the root node, we may want to be called during the cut-and-price loop and even before the first LP solve */
   if( heurdata->duringroot && SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_BEFORENODE);

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolOneopt)
{
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitOneopt)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize last solution index */
   heurdata->lastsolindex = -1;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecOneopt)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;
   SCIP_SOL* bestsol;                        /* incumbent solution                   */
   SCIP_SOL* worksol;                        /* heuristic's working solution         */
   SCIP_VAR** vars;                          /* SCIP variables                       */
   SCIP_VAR** shiftcands;                    /* shiftable variables                  */
   SCIP_ROW** lprows;                        /* SCIP LP rows                         */
   SCIP_Real* activities;                    /* row activities for working solution  */
   SCIP_Real* shiftvals;
   SCIP_Bool shifted;

   SCIP_RETCODE retcode;

   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool localrows;
   SCIP_Bool valid;
   int nchgbound;
   int nbinvars;
   int nintvars;
   int nvars;
   int nlprows;
   int i;
   int nshiftcands;
   int shiftcandssize;
   int nsuccessfulshifts;
   int niterations;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DELAYED;

   /* we only want to process each solution once */
   bestsol = SCIPgetBestSol(scip);
   if( bestsol == NULL || heurdata->lastsolindex == SCIPsolGetIndex(bestsol) )
      return SCIP_OKAY;

   /* reset the timing mask to its default value (at the root node it could be different) */
   if( SCIPgetNNodes(scip) > 1 )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);

   /* get problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nintvars += nbinvars;

   /* do not run if there are no discrete variables */
   if( nintvars == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL )
   {
      SCIP*                 subscip;            /* the subproblem created by oneopt                */
      SCIP_VAR**            subvars;            /* subproblem's variables                          */

      SCIP_Bool success;

      if( !heurdata->beforepresol )
         return SCIP_OKAY;

      /* check whether there is enough time and memory left */
      SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

      if( !success )
         return SCIP_OKAY;

      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

      /* initialize the subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      /* setup and solve the subproblem and catch the return code */
      retcode = setupAndSolveSubscipOneopt(scip, subscip, heur, vars, subvars, bestsol, result, &valid);

      /* free the subscip in any case */
      SCIP_CALL( SCIPfree(&subscip) );
      SCIP_CALL( retcode );

      SCIPfreeBufferArray(scip, &subvars);

      return SCIP_OKAY;
   }

   /* we can only work on solutions valid in the transformed space */
   if( SCIPsolIsOriginal(bestsol) )
      return SCIP_OKAY;

   if( heurtiming == SCIP_HEURTIMING_BEFORENODE && (SCIPhasCurrentNodeLP(scip) || heurdata->forcelpconstruction) )
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

      /* get problem variables again, SCIPconstructLP() might have added new variables */
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
      nintvars += nbinvars;
   }

   /* we need an LP */
   if( SCIPgetNLPRows(scip) == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   heurdata->lastsolindex = SCIPsolGetIndex(bestsol);
   SCIP_CALL( SCIPcreateSolCopy(scip, &worksol, bestsol) );
   SCIPsolSetHeur(worksol,heur);

   SCIPdebugMsg(scip, "Starting bound adjustment in 1-opt heuristic\n");

   nchgbound = 0;
   /* change solution values due to possible global bound changes first */
   for( i = nvars - 1; i >= 0; --i )
   {
      SCIP_VAR* var;
      SCIP_Real solval;

      var = vars[i];
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      solval = SCIPgetSolVal(scip, worksol, var);
      /* old solution value is smaller than the actual lower bound */
      if( SCIPisFeasLT(scip, solval, lb) )
      {
         /* set the solution value to the global lower bound */
         SCIP_CALL( SCIPsetSolVal(scip, worksol, var, lb) );
         ++nchgbound;
         SCIPdebugMsg(scip, "var <%s> type %d, old solval %g now fixed to lb %g\n", SCIPvarGetName(var), SCIPvarGetType(var), solval, lb);
      }
      /* old solution value is greater than the actual upper bound */
      else if( SCIPisFeasGT(scip, solval, SCIPvarGetUbGlobal(var)) )
      {
         /* set the solution value to the global upper bound */
         SCIP_CALL( SCIPsetSolVal(scip, worksol, var, ub) );
         ++nchgbound;
         SCIPdebugMsg(scip, "var <%s> type %d, old solval %g now fixed to ub %g\n", SCIPvarGetName(var), SCIPvarGetType(var), solval, ub);
      }
   }

   SCIPdebugMsg(scip, "number of bound changes (due to global bounds) = %d\n", nchgbound);

   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activities, nlprows) );

   localrows = FALSE;
   valid = TRUE;

   /* initialize LP row activities */
   for( i = 0; i < nlprows; ++i )
   {
      SCIP_ROW* row;

      row = lprows[i];
      assert(SCIProwGetLPPos(row) == i);

      if( !SCIProwIsLocal(row) )
      {
         activities[i] = SCIPgetRowSolActivity(scip, row, worksol);
         SCIPdebugMsg(scip, "Row <%s> has activity %g\n", SCIProwGetName(row), activities[i]);
         if( SCIPisFeasLT(scip, activities[i], SCIProwGetLhs(row)) || SCIPisFeasGT(scip, activities[i], SCIProwGetRhs(row)) )
         {
            valid = FALSE;
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            SCIPdebugMsg(scip, "row <%s> activity %g violates bounds, lhs = %g, rhs = %g\n", SCIProwGetName(row), activities[i], SCIProwGetLhs(row), SCIProwGetRhs(row));
            break;
         }
      }
      else
         localrows = TRUE;
   }

   if( !valid )
   {
      /** @todo try to correct lp rows */
      SCIPdebugMsg(scip, "Some global bound changes were not valid in lp rows.\n");

      SCIPfreeBufferArray(scip, &activities);
      SCIP_CALL( SCIPfreeSol(scip, &worksol) );

      return SCIP_OKAY;
   }

   /* allocate buffer storage for possible shift candidates */
   shiftcandssize = 8;
   SCIP_CALL( SCIPallocBufferArray(scip, &shiftcands, shiftcandssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &shiftvals, shiftcandssize) );
   nsuccessfulshifts = 0;
   niterations = 0;
   do
   {
      /* initialize data */
      shifted = FALSE;
      nshiftcands = 0;
      ++niterations;
      SCIPdebugMsg(scip, "Starting 1-opt heuristic iteration #%d\n", niterations);

      /* enumerate all integer variables and find out which of them are shiftable */
      /* @todo if useloop=TRUE store for each variable which constraint blocked it and only iterate over those variables
       *       in the following rounds for which the constraint slack was increased by previous shifts
       */
      for( i = 0; i < nintvars; i++ )
      {
         if( SCIPvarGetStatus(vars[i]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_Real shiftval;
            SCIP_Real solval;

            /* find out whether the variable can be shifted */
            solval = SCIPgetSolVal(scip, worksol, vars[i]);
            shiftval = calcShiftVal(scip, vars[i], solval, activities);

            /* insert the variable into the list of shifting candidates */
            if( !SCIPisFeasZero(scip, shiftval) )
            {
               SCIPdebugMsg(scip, " -> Variable <%s> can be shifted by <%1.1f> \n", SCIPvarGetName(vars[i]), shiftval);

               if( nshiftcands == shiftcandssize)
               {
                  shiftcandssize *= 8;
                  SCIP_CALL( SCIPreallocBufferArray(scip, &shiftcands, shiftcandssize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &shiftvals, shiftcandssize) );
               }
               shiftcands[nshiftcands] = vars[i];
               shiftvals[nshiftcands] = shiftval;
               nshiftcands++;
            }
         }
      }

      /* if at least one variable can be shifted, shift variables sorted by their objective */
      if( nshiftcands > 0 )
      {
         SCIP_Real shiftval;
         SCIP_Real solval;
         SCIP_VAR* var;

         /* the case that exactly one variable can be shifted is slightly easier */
         if( nshiftcands == 1 )
         {
            var = shiftcands[0];
            assert(var != NULL);
            solval = SCIPgetSolVal(scip, worksol, var);
            shiftval = shiftvals[0];
            assert(!SCIPisFeasZero(scip,shiftval));
            SCIPdebugMsg(scip, " Only one shiftcand found, var <%s>, which is now shifted by<%1.1f> \n",
               SCIPvarGetName(var), shiftval);
            SCIP_CALL( SCIPsetSolVal(scip, worksol, var, solval+shiftval) );
            SCIP_CALL( updateRowActivities(scip, activities, var, shiftval) );
            ++nsuccessfulshifts;
         }
         else
         {
            SCIP_Real* objcoeffs;

            SCIP_CALL( SCIPallocBufferArray(scip, &objcoeffs, nshiftcands) );

            SCIPdebugMsg(scip, " %d shiftcands found \n", nshiftcands);

            /* sort the variables by their objective, optionally weighted with the shiftval */
            if( heurdata->weightedobj )
            {
               for( i = 0; i < nshiftcands; ++i )
                  objcoeffs[i] = SCIPvarGetObj(shiftcands[i])*shiftvals[i];
            }
            else
            {
               for( i = 0; i < nshiftcands; ++i )
                  objcoeffs[i] = SCIPvarGetObj(shiftcands[i]);
            }

            /* sort arrays with respect to the first one */
            SCIPsortRealPtr(objcoeffs, (void**)shiftcands, nshiftcands);

            /* try to shift each variable -> Activities have to be updated */
            for( i = 0; i < nshiftcands; ++i )
            {
               var = shiftcands[i];
               assert(var != NULL);
               solval = SCIPgetSolVal(scip, worksol, var);
               shiftval = calcShiftVal(scip, var, solval, activities);
               assert(i > 0 || !SCIPisFeasZero(scip, shiftval));
               assert(SCIPisFeasGE(scip, solval+shiftval, SCIPvarGetLbGlobal(var)) && SCIPisFeasLE(scip, solval+shiftval, SCIPvarGetUbGlobal(var)));

               /* update data structures for nonzero shift value */
               if( ! SCIPisFeasZero(scip, shiftval) )
               {
                  SCIPdebugMsg(scip, " -> Variable <%s> is now shifted by <%1.1f> \n", SCIPvarGetName(vars[i]), shiftval);
                  SCIP_CALL( SCIPsetSolVal(scip, worksol, var, solval+shiftval) );
                  SCIP_CALL( updateRowActivities(scip, activities, var, shiftval) );
                  ++nsuccessfulshifts;
               }
            }

            SCIPfreeBufferArray(scip, &objcoeffs);
         }
         shifted = TRUE;
      }

   } while( heurdata->useloop && shifted );

   if( nsuccessfulshifts > 0 )
   {
      /* if the problem is a pure IP, try to install the solution, if it is a MIP, solve LP again to set the continuous
       * variables to the best possible value
       */
      if( nvars == nintvars || !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIP_Bool success;

         /* since we ignore local rows, we cannot guarantee their feasibility and have to set the checklprows flag to
          * TRUE if local rows are present
          */
         SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, localrows, &success) );

         if( success )
         {
            SCIPdebugMsg(scip, "found feasible shifted solution:\n");
            SCIPdebug( SCIP_CALL( SCIPprintSol(scip, worksol, NULL, FALSE) ) );
            *result = SCIP_FOUNDSOL;
         }
      }
      else
      {
         SCIP_Bool lperror;
#ifdef NDEBUG
         SCIP_RETCODE retstat;
#endif

         SCIPdebugMsg(scip, "shifted solution should be feasible -> solve LP to fix continuous variables to best values\n");

         /* start diving to calculate the LP relaxation */
         SCIP_CALL( SCIPstartDive(scip) );

         /* set the bounds of the variables: fixed for integers, global bounds for continuous */
         for( i = 0; i < nvars; ++i )
         {
            if( SCIPvarGetStatus(vars[i]) == SCIP_VARSTATUS_COLUMN )
            {
               SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], SCIPvarGetLbGlobal(vars[i])) );
               SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], SCIPvarGetUbGlobal(vars[i])) );
            }
         }
         /* apply this after global bounds to not cause an error with intermediate empty domains */
         for( i = 0; i < nintvars; ++i )
         {
            if( SCIPvarGetStatus(vars[i]) == SCIP_VARSTATUS_COLUMN )
            {
               SCIP_Real solval;
               solval = SCIPgetSolVal(scip, worksol, vars[i]);
               SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], solval) );
               SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], solval) );
            }
         }

         /* solve LP */
         SCIPdebugMsg(scip, " -> old LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));

         /**@todo in case of an MINLP, if SCIPisNLPConstructed() is TRUE, say, rather solve the NLP instead of the LP */
         /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
          * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
          */
#ifdef NDEBUG
         retstat = SCIPsolveDiveLP(scip, -1, &lperror, NULL);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while solving LP in 1-opt heuristic; LP solve terminated with code <%d>\n",retstat);
         }
#else
         SCIP_CALL( SCIPsolveDiveLP(scip, -1, &lperror, NULL) );
#endif

         SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
         SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, SCIPgetLPSolstat(scip));

         /* check if this is a feasible solution */
         if( !lperror && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIP_Bool success;

            /* copy the current LP solution to the working solution */
            SCIP_CALL( SCIPlinkLPSol(scip, worksol) );
            SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check solution for feasibility */
            if( success )
            {
               SCIPdebugMsg(scip, "found feasible shifted solution:\n");
               SCIPdebug( SCIP_CALL( SCIPprintSol(scip, worksol, NULL, FALSE) ) );
               *result = SCIP_FOUNDSOL;
            }
         }

         /* terminate the diving */
         SCIP_CALL( SCIPendDive(scip) );
      }
   }

   /* heuristic should not rerun on this incumbent because the heuristic loop finishes only after no further
    * improvements of the incumbent solution are possible
    */
   if( heurdata->useloop )
      heurdata->lastsolindex = SCIPsolGetIndex(SCIPgetBestSol(scip));

   SCIPfreeBufferArray(scip, &shiftvals);
   SCIPfreeBufferArray(scip, &shiftcands);
   SCIPfreeBufferArray(scip, &activities);

   SCIP_CALL( SCIPfreeSol(scip, &worksol) );

   SCIPdebugMsg(scip, "Finished 1-opt heuristic\n");

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the oneopt primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurOneopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Oneopt primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecOneopt, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyOneopt) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeOneopt) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolOneopt) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolOneopt) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitOneopt) );

   /* add oneopt primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/oneopt/weightedobj",
         "should the objective be weighted with the potential shifting value when sorting the shifting candidates?",
         &heurdata->weightedobj, TRUE, DEFAULT_WEIGHTEDOBJ, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/oneopt/duringroot",
         "should the heuristic be called before and during the root node?",
         &heurdata->duringroot, TRUE, DEFAULT_DURINGROOT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/oneopt/forcelpconstruction",
         "should the construction of the LP be forced even if LP solving is deactivated?",
         &heurdata->forcelpconstruction, TRUE, DEFAULT_FORCELPCONSTRUCTION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/oneopt/beforepresol",
         "should the heuristic be called before presolving?",
         &heurdata->beforepresol, TRUE, DEFAULT_BEFOREPRESOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/oneopt/useloop",
         "should the heuristic continue to run as long as improvements are found?",
         &heurdata->useloop, TRUE, DEFAULT_USELOOP, NULL, NULL) );

   return SCIP_OKAY;
}
