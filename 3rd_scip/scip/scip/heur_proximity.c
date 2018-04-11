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

/**@file   heur_proximity.c
 * @brief  improvement heuristic which uses an auxiliary objective instead of the original objective function which
 *         is itself added as a constraint to a sub-SCIP instance. The heuristic was presented by Matteo Fischetti
 *         and Michele Monaci.
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_proximity.h"
#include "scip/cons_linear.h"

#define HEUR_NAME             "proximity"
#define HEUR_DESC             "heuristic trying to improve the incumbent by an auxiliary proximity objective function"
#define HEUR_DISPCHAR         'P'
#define HEUR_PRIORITY         -2000000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

/* event handler properties */
#define EVENTHDLR_NAME        "Proximity"
#define EVENTHDLR_DESC        "LP event handler for " HEUR_NAME " heuristic"

/* default values for proximity-specific parameters */
/* todo refine these values */
#define DEFAULT_MAXNODES      10000LL    /**< maximum number of nodes to regard in the subproblem                        */
#define DEFAULT_MINIMPROVE    0.02       /**< factor by which proximity should at least improve the incumbent            */
#define DEFAULT_MINGAP        0.01       /**< minimum primal-dual gap for which the heuristic is executed                */
#define DEFAULT_MINNODES      1LL        /**< minimum number of nodes to regard in the subproblem                        */
#define DEFAULT_MINLPITERS    200LL      /**< minimum number of LP iterations to perform in one sub-mip                  */
#define DEFAULT_MAXLPITERS    100000LL   /**< maximum number of LP iterations to be performed in the subproblem          */
#define DEFAULT_NODESOFS      50LL       /**< number of nodes added to the contingent of the total nodes                 */
#define DEFAULT_WAITINGNODES  100LL      /**< default waiting nodes since last incumbent before heuristic is executed    */
#define DEFAULT_NODESQUOT     0.1        /**< default quotient of sub-MIP nodes with respect to number of processed nodes*/
#define DEFAULT_USELPROWS     FALSE      /**< should subproblem be constructed based on LP row information? */
#define DEFAULT_BINVARQUOT    0.1        /**< default threshold for percentage of binary variables required to start     */
#define DEFAULT_RESTART       TRUE       /**< should the heuristic immediately run again on its newly found solution? */
#define DEFAULT_USEFINALLP    FALSE      /**< should the heuristic solve a final LP in case of continuous objective variables? */
#define DEFAULT_LPITERSQUOT   0.2        /**< default quotient of sub-MIP LP iterations with respect to LP iterations so far */
#define DEFAULT_USEUCT        FALSE      /**< should uct node selection be used at the beginning of the search?     */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          maxlpiters;         /**< maximum number of LP iterations to be performed in the subproblem   */
   SCIP_Longint          nusedlpiters;       /**< number of actually performed LP iterations                          */
   SCIP_Longint          minlpiters;         /**< minimum number of LP iterations to perform in one sub-mip           */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;          /**< nodes already used by proximity in earlier calls                    */
   SCIP_Longint          waitingnodes;       /**< waiting nodes since last incumbent before heuristic is executed     */
   SCIP_Real             lpitersquot;        /**< quotient of sub-MIP LP iterations with respect to LP iterations so far */
   SCIP_Real             minimprove;         /**< factor by which proximity should at least improve the incumbent     */
   SCIP_Real             mingap;             /**< minimum primal-dual gap for which the heuristic is executed         */
   SCIP_Real             nodesquot;          /**< quotient of sub-MIP nodes with respect to number of processed nodes */
   SCIP_Real             binvarquot;         /**<  threshold for percantage of binary variables required to start     */

   SCIP*                 subscip;            /**< the subscip used by the heuristic                                   */
   SCIP_HASHMAP*         varmapfw;           /**< map between scip variables and subscip variables                    */
   SCIP_VAR**            subvars;            /**< variables in subscip                                                */
   SCIP_CONS*            objcons;            /**< the objective cutoff constraint of the subproblem                   */

   int                   nsubvars;           /**< the number of subvars                                               */
   int                   lastsolidx;         /**< index of last solution on which the heuristic was processed         */
   int                   subprobidx;         /**< counter for the subproblem index to be solved by proximity */

   SCIP_Bool             uselprows;          /**< should subproblem be constructed based on LP row information? */
   SCIP_Bool             restart;            /**< should the heuristic immediately run again on its newly found solution? */
   SCIP_Bool             usefinallp;         /**< should the heuristic solve a final LP in case of continuous objective variables? */
   SCIP_Bool             useuct;             /**< should uct node selection be used at the beginning of the search?  */
};


/*
 * Local methods
 */

/** optimizes the continuous variables in an LP diving by fixing all integer variables to the given solution values */
static
SCIP_RETCODE solveLp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< candidate solution for which continuous variables should be optimized */
   SCIP_Bool*            success             /**< was the dive successful? */
   )
{
   SCIP_VAR** vars;
   SCIP_RETCODE retstat;

   int v;
   int nvars;
   int ncontvars;
   int nintvars;

   SCIP_Bool lperror;
   SCIP_Bool requiresnlp;

   assert(success != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );

   nintvars = nvars - ncontvars;

   /**@todo in case of an MINLP, if SCIPisNLPConstructed() is TRUE rather solve the NLP instead of the LP */
   requiresnlp = SCIPisNLPConstructed(scip);
   if( requiresnlp || ncontvars == 0 )
      return SCIP_OKAY;

   /* start diving to calculate the LP relaxation */
   SCIP_CALL( SCIPstartDive(scip) );

   /* set the bounds of the variables: fixed for integers, global bounds for continuous */
   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPchgVarLbDive(scip, vars[v], SCIPvarGetLbGlobal(vars[v])) );
         SCIP_CALL( SCIPchgVarUbDive(scip, vars[v], SCIPvarGetUbGlobal(vars[v])) );
      }
   }

   /* apply this after global bounds to not cause an error with intermediate empty domains */
   for( v = 0; v < nintvars; ++v )
   {
      if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_Real solval;

         solval = SCIPgetSolVal(scip, sol, vars[v]);
         SCIP_CALL( SCIPchgVarLbDive(scip, vars[v], solval) );
         SCIP_CALL( SCIPchgVarUbDive(scip, vars[v], solval) );
      }
   }

   /* solve LP */
   SCIPdebugMsg(scip, " -> old LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));

   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   retstat = SCIPsolveDiveLP(scip, -1, &lperror, NULL);
   if( retstat != SCIP_OKAY )
   {
#ifdef NDEBUG
      SCIPwarningMessage(scip, "Error while solving LP in Proximity heuristic; LP solve terminated with code <%d>\n",retstat);
#else
      SCIP_CALL( retstat );
#endif
   }

   SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
   SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, SCIPgetLPSolstat(scip));
   if( !lperror && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_CALL( SCIPlinkLPSol(scip, sol) );
      SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );
   }

   /* terminate diving mode */
   SCIP_CALL( SCIPendDive(scip) );

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< proximity heuristic structure                       */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool             usefinallp,         /**< should continuous variables be optimized by a final LP */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;                         /* the original problem's number of variables      */
   int        ncontvars;                     /* the original problem's number of continuous variables */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);
   assert(success != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );

   /* The sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP. */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   *success = FALSE;

   /* solve an LP with all integer variables fixed to improve solution quality */
   if( ncontvars > 0 && usefinallp && SCIPisLPConstructed(scip) )
   {
      int v;
      int ncontobjvars = 0;    /* does the problem instance have continuous variables with nonzero objective coefficients? */
      SCIP_Real sumofobjsquares = 0.0;

      /* check if continuous variables with nonzero objective coefficient are present */
      for( v = nvars - 1; v >= nvars - ncontvars; --v )
      {
         SCIP_VAR* var;

         var = vars[v];
         assert(vars[v] != NULL);
         assert(!SCIPvarIsIntegral(var));

         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && !SCIPisZero(scip, SCIPvarGetObj(var)) )
         {
            ++ncontobjvars;
            sumofobjsquares += SCIPvarGetObj(var) * SCIPvarGetObj(var);
         }
      }

      SCIPstatisticMessage(" Continuous Objective variables: %d, Euclidean OBJ: %g total, %g continuous\n", ncontobjvars, SCIPgetObjNorm(scip), sumofobjsquares);

      /* solve a final LP to optimize solution values of continuous problem variables */
      SCIPstatisticMessage("Solution Value before LP resolve: %g\n", SCIPgetSolOrigObj(scip, newsol));
      SCIP_CALL( solveLp(scip, newsol, success) );

      /* if the LP solve was not successful, reset the solution */
      if( !*success )
      {
         for( v = nvars - 1; v >= nvars - ncontvars; --v )
         {
            SCIP_CALL( SCIPsetSolVal(scip, newsol, vars[v], subsolvals[v]) );
         }

      }
   }

   /* try to add new solution to SCIP and free it immediately */
   if( !*success )
   {
      SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );
   }
   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** sets solving parameters for the subproblem created by the heuristic */
static
SCIP_RETCODE setupSubproblem(
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP*                 subscip             /**< copied SCIP data structure */
   )
{
   assert(subscip != NULL);

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

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* use restart dfs node selection */
   if( SCIPfindNodesel(subscip, "restartdfs") != NULL && !SCIPisParamFixed(subscip, "nodeselection/restartdfs/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/restartdfs/stdpriority", INT_MAX/4) );
   }

   /* activate uct node selection at the top of the tree */
   if( heurdata->useuct && SCIPfindNodesel(subscip, "uct") != NULL && !SCIPisParamFixed(subscip, "nodeselection/uct/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/uct/stdpriority", INT_MAX/2) );
   }

   /* disable expensive presolving
    * todo maybe presolving can be entirely turned off here - parameter???
    */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) ); */
   if( !SCIPisParamFixed(subscip, "presolving/maxrounds") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "presolving/maxrounds", 50) );
   }

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* todo: check branching rule in sub-SCIP */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable feasibility pump and fractional diving */
   if( !SCIPisParamFixed(subscip, "heuristics/feaspump/freq") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/feaspump/freq", -1) );
   }
   if( !SCIPisParamFixed(subscip, "heuristics/fracdiving/freq") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/fracdiving/freq", -1) );
   }

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 500) );
   }

   /* todo check if
    * SCIP_CALL( SCIPsetEmphasis(subscip, SCIP_PARAMEMPHASIS_FEASIBILITY, TRUE) );
    * improves performance */

   return SCIP_OKAY;
}

/** frees the subproblem */
static
SCIP_RETCODE deleteSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   /* free remaining memory from heuristic execution */
   if( heurdata->subscip != NULL )
   {
      assert(heurdata->varmapfw != NULL);
      assert(heurdata->subvars != NULL);
      assert(heurdata->objcons != NULL);

      SCIPdebugMsg(scip, "Freeing subproblem of proximity heuristic\n");
      SCIPfreeBlockMemoryArray(scip, &heurdata->subvars, heurdata->nsubvars);
      SCIPhashmapFree(&heurdata->varmapfw);
      SCIP_CALL( SCIPreleaseCons(heurdata->subscip, &heurdata->objcons) );
      SCIP_CALL( SCIPfree(&heurdata->subscip) );

      heurdata->subscip = NULL;
      heurdata->varmapfw = NULL;
      heurdata->subvars = NULL;
      heurdata->objcons = NULL;
   }
   return SCIP_OKAY;
}

/* ---------------- Callback methods of event handler ---------------- */

/** exec the event handler
 *
 *  We interrupt the solution process.
 */
static
SCIP_DECL_EVENTEXEC(eventExecProximity)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_NODESOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP
    * todo adjust interruption limit */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_ITERLIMIT || SCIPgetNLPIterations(scip) >= heurdata->maxlpiters )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}


/* ---------------- Callback methods of primal heuristic ---------------- */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyProximity)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurProximity(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeProximity)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitProximity)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->usednodes = 0LL;
   heurdata->lastsolidx = -1;
   heurdata->nusedlpiters = 0LL;
   heurdata->subprobidx = 0;

   heurdata->subscip = NULL;
   heurdata->varmapfw = NULL;
   heurdata->subvars = NULL;
   heurdata->objcons = NULL;

   heurdata->nsubvars = 0;

   return SCIP_OKAY;
}

/** solution process exiting method of proximity heuristic */
static
SCIP_DECL_HEUREXITSOL(heurExitsolProximity)
{
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   SCIP_CALL( deleteSubproblem(scip, heurdata) );

   assert(heurdata->subscip == NULL && heurdata->varmapfw == NULL && heurdata->subvars == NULL && heurdata->objcons == NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecProximity)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata; /* heuristic's data                            */
   SCIP_Longint nnodes;     /* number of stalling nodes for the subproblem */
   SCIP_Longint nlpiters;   /* lp iteration limit for the subproblem       */
   SCIP_Bool foundsol = FALSE;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* do not run heuristic when there are only few binary varables */
   if( SCIPgetNBinVars(scip) < heurdata->binvarquot * SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* calculate branching node limit for sub problem */
   /* todo maybe treat root node differently */
   nnodes = (SCIP_Longint) (heurdata->nodesquot * SCIPgetNNodes(scip));
   nnodes += heurdata->nodesofs;

   /* determine the node and LP iteration limit for the solve of the sub-SCIP */
   nnodes -= heurdata->usednodes;
   nnodes = MIN(nnodes, heurdata->maxnodes);

   nlpiters = (SCIP_Longint) (heurdata->lpitersquot * SCIPgetNRootFirstLPIterations(scip));
   nlpiters = MIN(nlpiters, heurdata->maxlpiters);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping proximity: nnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* do not run proximity, if the problem does not have an objective function anyway */
   if( SCIPgetNObjVars(scip) == 0 )
   {
      SCIPdebugMsg(scip, "skipping proximity: pure feasibility problem anyway\n");
      return SCIP_OKAY;
   }

   do
   {
      /* main loop of proximity: in every iteration, a new subproblem is set up and solved until no improved solution
       * is found or one of the heuristic limits on nodes or LP iterations is hit
       * heuristic performs only one iteration if restart parameter is set to FALSE
       */
      SCIP_Longint nusednodes = 0LL;
      SCIP_Longint nusedlpiters = 0LL;

      nlpiters = MAX(nlpiters, heurdata->minlpiters);

      /* define and solve the proximity subproblem */
      SCIP_CALL( SCIPapplyProximity(scip, heur, result, heurdata->minimprove, nnodes, nlpiters, &nusednodes, &nusedlpiters, FALSE) );

      /* adjust node limit and LP iteration limit for future iterations */
      assert(nusednodes <= nnodes);
      heurdata->usednodes += nusednodes;
      nnodes -= nusednodes;

      nlpiters -= nusedlpiters;
      heurdata->nusedlpiters += nusedlpiters;

      /* memorize if a new solution has been found in at least one iteration */
      if( *result == SCIP_FOUNDSOL )
         foundsol = TRUE;
   }
   while( *result == SCIP_FOUNDSOL && heurdata->restart && !SCIPisStopped(scip) && nnodes > 0 );

   /* reset result pointer if solution has been found in previous iteration */
   if( foundsol )
      *result = SCIP_FOUNDSOL;

   /* free the occupied memory */
   if( heurdata->subscip != NULL )
   {
      /* just for testing the library method, in debug mode, we call the wrapper method for the actual delete method */
#ifndef NDEBUG
      SCIP_CALL( SCIPdeleteSubproblemProximity(scip) );
#else
      SCIP_CALL( deleteSubproblem(scip, heurdata) );
#endif
   }
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** frees the sub-MIP created by proximity */
SCIP_RETCODE SCIPdeleteSubproblemProximity(
   SCIP*                 scip                /** SCIP data structure */
   )
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);

   heur = SCIPfindHeur(scip, HEUR_NAME);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   if( heurdata != NULL )
   {
      SCIP_CALL( deleteSubproblem(scip, heurdata) );
   }

   return SCIP_OKAY;
}

/** main procedure of the proximity heuristic, creates and solves a sub-SCIP
 *
 *  @note The method can be applied in an iterative way, keeping the same subscip in between. If the @p freesubscip
 *        parameter is set to FALSE, the heuristic will keep the subscip data structures. Always set this parameter
 *        to TRUE, or call SCIPdeleteSubproblemProximity() afterwards.
 */
SCIP_RETCODE SCIPapplyProximity(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                            */
   SCIP_RESULT*          result,             /**< result data structure                                               */
   SCIP_Real             minimprove,         /**< factor by which proximity should at least improve the incumbent     */
   SCIP_Longint          nnodes,             /**< node limit for the subproblem                                       */
   SCIP_Longint          nlpiters,           /**< LP iteration limit for the subproblem                               */
   SCIP_Longint*         nusednodes,         /**< pointer to store number of used nodes in subscip                    */
   SCIP_Longint*         nusedlpiters,       /**< pointer to store number of used LP iterations in subscip            */
   SCIP_Bool             freesubscip         /**< should the created sub-MIP be freed at the end of the method?       */
   )
{
   SCIP*                 subscip;            /* the subproblem created by proximity              */
   SCIP_HASHMAP*         varmapfw;           /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR**            vars;               /* original problem's variables                    */
   SCIP_VAR**            subvars;            /* subproblem's variables                          */
   SCIP_HEURDATA*        heurdata;           /* heuristic's private data structure              */
   SCIP_EVENTHDLR*       eventhdlr;          /* event handler for LP events                     */

   SCIP_SOL* incumbent;
   SCIP_CONS* objcons;
   SCIP_Longint iterlim;

   SCIP_Real large;
   SCIP_Real inf;

   SCIP_Real bestobj;
   SCIP_Real objcutoff;
   SCIP_Real lowerbound;

   int nvars;                                /* number of original problem's variables          */
   int nfixedvars;
   int nsubsols;
   int solidx;
   int i;

   SCIP_Bool valid;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   assert(nnodes >= 0);
   assert(0.0 <= minimprove && minimprove <= 1.0);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only call the heuristic if we have an incumbent  */
   if( SCIPgetNSolsFound(scip) == 0 )
      return SCIP_OKAY;

   /* do not use heuristic on problems without binary variables */
   if( SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;

   incumbent = SCIPgetBestSol(scip);
   assert(incumbent != NULL);

   /* make sure that the incumbent is valid for the transformed space, otherwise terminate */
   if( SCIPsolIsOriginal(incumbent) )
      return SCIP_OKAY;

   solidx = SCIPsolGetIndex(incumbent);

   if( heurdata->lastsolidx == solidx )
      return SCIP_OKAY;

   /* only call heuristic, if the best solution does not come from trivial heuristic */
   if( SCIPsolGetHeur(incumbent) != NULL && strcmp(SCIPheurGetName(SCIPsolGetHeur(incumbent)), "trivial") == 0 )
      return SCIP_OKAY;

   /* waitingnodes parameter defines the minimum number of nodes to wait before a new incumbent is processed */
   if( SCIPgetNNodes(scip) > 1 && SCIPgetNNodes(scip) - SCIPsolGetNodenum(incumbent) < heurdata->waitingnodes )
      return SCIP_OKAY;

   bestobj = SCIPgetSolTransObj(scip, incumbent);
   lowerbound = SCIPgetLowerbound(scip);

   /* use knowledge about integrality of objective to round up lower bound */
   if( SCIPisObjIntegral(scip) )
   {
      SCIPdebugMsg(scip, " Rounding up lower bound: %f --> %f \n", lowerbound, SCIPfeasCeil(scip, lowerbound));
      lowerbound = SCIPfeasCeil(scip, lowerbound);
   }

   /* do not trigger heuristic if primal and dual bound are already close together */
   if( SCIPisFeasLE(scip, bestobj, lowerbound) || SCIPgetGap(scip) <= heurdata->mingap )
      return SCIP_OKAY;

   /* calculate the minimum improvement for a heuristic solution in terms of the distance between incumbent objective
    * and the lower bound */
   if( SCIPisInfinity(scip, REALABS(lowerbound)) )
   {
      if( SCIPisZero(scip, bestobj) )
         objcutoff = bestobj - 1;
      else
         objcutoff = (1 - minimprove) * bestobj;
   }
   else
      objcutoff = minimprove * lowerbound + (1 - minimprove) * (bestobj);

   /* use integrality of the objective function to round down (and thus strengthen) the objective cutoff */
   if( SCIPisObjIntegral(scip) )
      objcutoff = SCIPfeasFloor(scip, objcutoff);

   if( SCIPisFeasLT(scip, objcutoff, lowerbound) )
      objcutoff = lowerbound;

   /* exit execution if the right hand side of the objective constraint does not change (suggests that the heuristic
    * was not successful in a previous iteration) */
   if( heurdata->objcons != NULL && SCIPisFeasEQ(scip, SCIPgetRhsLinear(heurdata->subscip, heurdata->objcons), objcutoff) )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &valid) );

   if( ! valid )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   heurdata->lastsolidx = solidx;

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create a subscip and copy the original scip instance into it */
   if( heurdata->subscip == NULL )
   {
      assert(heurdata->varmapfw == NULL);
      assert(heurdata->objcons == NULL);

      /* initialize the subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &subvars, nvars) );

      /* copy complete SCIP instance */
      valid = FALSE;

      /* create a problem copy as sub SCIP */
      SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "proximity", NULL, NULL, 0, heurdata->uselprows, TRUE,
            &success, &valid) );

      SCIPdebugMsg(scip, "Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");

      /* create event handler for LP events */
      eventhdlr = NULL;
      SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecProximity, NULL) );
      if( eventhdlr == NULL )
      {
         SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
         return SCIP_PLUGINNOTFOUND;
      }

      /* set up parameters for the copied instance */
      SCIP_CALL( setupSubproblem(heurdata, subscip) );

      /* create the objective constraint in the sub scip, first without variables and values which will be added later */
      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &objcons, "objbound_of_origscip", 0, NULL, NULL, -SCIPinfinity(subscip), SCIPinfinity(subscip)) );

      /* determine large value to set variable bounds to, safe-guard to avoid fixings to infinite values */
      large = SCIPinfinity(scip);
      if( !SCIPisInfinity(scip, 0.1 / SCIPfeastol(scip)) )
         large = 0.1 / SCIPfeastol(scip);
      inf = SCIPinfinity(subscip);

      /* get variable image and change objective to proximity function (Manhattan distance) in sub-SCIP */
      for( i = 0; i < nvars; i++ )
      {
         SCIP_Real adjustedbound;
         SCIP_Real lb;
         SCIP_Real ub;

         subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

         SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], 0.0) );

         lb = SCIPvarGetLbGlobal(subvars[i]);
         ub = SCIPvarGetUbGlobal(subvars[i]);

         /* adjust infinite bounds in order to avoid that variables with non-zero objective
          * get fixed to infinite value in proximity subproblem
          */
         if( SCIPisInfinity(subscip, ub) )
         {
            adjustedbound = MAX(large, lb + large);
            adjustedbound = MIN(adjustedbound, inf);
            SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], adjustedbound) );
         }
         if( SCIPisInfinity(subscip, -lb) )
         {
            adjustedbound = MIN(-large, ub - large);
            adjustedbound = MAX(adjustedbound, -inf);
            SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], adjustedbound) );
         }

         /* add all nonzero objective coefficients to the objective constraint */
         if( !SCIPisFeasZero(subscip, SCIPvarGetObj(vars[i])) )
         {
            SCIP_CALL( SCIPaddCoefLinear(subscip, objcons, subvars[i], SCIPvarGetObj(vars[i])) );
         }
      }

      /* add objective constraint to the subscip */
      SCIP_CALL( SCIPaddCons(subscip, objcons) );
   }
   else
   {
      /* the instance, event handler, hash map and variable array were already copied in a previous iteration
       * and stored in heuristic data
       */
      assert(heurdata->varmapfw != NULL);
      assert(heurdata->subvars != NULL);
      assert(heurdata->objcons != NULL);

      subscip = heurdata->subscip;
      varmapfw = heurdata->varmapfw;
      subvars = heurdata->subvars;
      objcons = heurdata->objcons;

      eventhdlr = SCIPfindEventhdlr(subscip, EVENTHDLR_NAME);
      assert(eventhdlr != NULL);
   }

   SCIP_CALL( SCIPchgRhsLinear(subscip, objcons, objcutoff) );

   for( i = 0; i < SCIPgetNBinVars(scip); ++i )
   {
      SCIP_Real solval;

      /* objective coefficients are only set for binary variables of the problem */
      assert(SCIPvarIsBinary(subvars[i]));

      solval = SCIPgetSolVal(scip, incumbent, vars[i]);
      assert(SCIPisFeasGE(scip, solval, 0.0));
      assert(SCIPisFeasLE(scip, solval, 1.0));
      assert(SCIPisFeasIntegral(scip, solval));

      if( solval < 0.5 )
      {
         SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], -1.0) );
      }
   }

   /* set limits for the subproblem */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nnodes) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/solutions", 1) );

   /* restrict LP iterations */
   /* todo set iterations limit depending on the number of iterations of the original problem root */
   iterlim = nlpiters;
   SCIP_CALL( SCIPsetLongintParam(subscip, "lp/iterlim", MAX(1, iterlim / MIN(10, nnodes))) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "lp/rootiterlim", iterlim) );

   /* catch LP events of sub-SCIP */
   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   SCIPstatisticMessage("solving subproblem at Node: %" SCIP_LONGINT_FORMAT " "
         "nnodes: %" SCIP_LONGINT_FORMAT " "
         "iterlim: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNNodes(scip), nnodes, iterlim);

   /* solve the subproblem with all previously adjusted parameters */
   nfixedvars = SCIPgetNFixedVars(subscip);

   SCIP_CALL( SCIPpresolve(subscip) );

   nfixedvars = SCIPgetNFixedVars(subscip) - nfixedvars;
   assert(nfixedvars >= 0);
   SCIPstatisticMessage("presolve fixings %d: %d\n", ++(heurdata->subprobidx), nfixedvars);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );
   SCIPstatisticMessage("solve of subscip %d:"
         "usednodes: %" SCIP_LONGINT_FORMAT " "
         "lp iters: %" SCIP_LONGINT_FORMAT " "
         "root iters: %" SCIP_LONGINT_FORMAT " "
         "Presolving Time: %.2f\n", heurdata->subprobidx,
         SCIPgetNNodes(subscip), SCIPgetNLPIterations(subscip), SCIPgetNRootLPIterations(subscip), SCIPgetPresolvingTime(subscip));

   SCIPstatisticMessage("Solving Time %d: %.2f\n", heurdata->subprobidx, SCIPgetSolvingTime(subscip) );

   /* drop LP events of sub-SCIP */
   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* keep track of relevant information for future runs of heuristic */
   if( nusednodes != NULL )
      *nusednodes = SCIPgetNNodes(subscip);
   if( nusedlpiters != NULL )
      *nusedlpiters = SCIPgetNLPIterations(subscip);

   /* check whether a solution was found */
   nsubsols = SCIPgetNSols(subscip);
   incumbent = SCIPgetBestSol(subscip);
   assert(nsubsols == 0 || incumbent != NULL);

   SCIPstatisticMessage("primal bound before subproblem %d: %g\n", heurdata->subprobidx, SCIPgetPrimalbound(scip));
   if( nsubsols > 0 )
   {
      /* try to translate the sub problem solution to the original scip instance */
      success = FALSE;
      SCIP_CALL( createNewSol(scip, subscip, subvars, heur, incumbent, heurdata->usefinallp, &success) );

      if( success )
         *result = SCIP_FOUNDSOL;
   }
   SCIPstatisticMessage("primal bound after subproblem %d: %g\n", heurdata->subprobidx, SCIPgetPrimalbound(scip));

   /* free the transformed subproblem data */
   SCIP_CALL( SCIPfreeTransform(subscip) );

   /* save subproblem in heuristic data for subsequent runs if it has been successful, otherwise free subproblem */
   heurdata->subscip = subscip;
   heurdata->varmapfw = varmapfw;
   heurdata->subvars = subvars;
   heurdata->objcons = objcons;
   heurdata->nsubvars = nvars;

   /* delete the sub problem */
   if( freesubscip )
   {
      SCIP_CALL( deleteSubproblem(scip, heurdata) );
   }

   return SCIP_OKAY;
}


/** creates the proximity primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurProximity(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur = NULL;

   /* create heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecProximity, heurdata) );
   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyProximity) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeProximity) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitProximity) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolProximity) );

   /* add proximity primal heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be constructed based on LP row information?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/restart",
         "should the heuristic immediately run again on its newly found solution?",
         &heurdata->restart, TRUE, DEFAULT_RESTART, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usefinallp",
         "should the heuristic solve a final LP in case of continuous objective variables?",
         &heurdata->usefinallp, TRUE, DEFAULT_USEFINALLP, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, TRUE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxlpiters",
         "maximum number of LP iterations to be performed in the subproblem",
         &heurdata->maxlpiters, TRUE, DEFAULT_MAXLPITERS, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minlpiters",
         "minimum number of LP iterations performed in subproblem",
         &heurdata->minlpiters, TRUE, DEFAULT_MINLPITERS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/waitingnodes",
          "waiting nodes since last incumbent before heuristic is executed",
         &heurdata->waitingnodes, TRUE, DEFAULT_WAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which proximity should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "sub-MIP node limit w.r.t number of original nodes",
         &heurdata->nodesquot, TRUE, DEFAULT_NODESQUOT, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/binvarquot",
         "threshold for percentage of binary variables required to start",
         &heurdata->binvarquot, TRUE, DEFAULT_BINVARQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lpitersquot",
            "quotient of sub-MIP LP iterations with respect to LP iterations so far",
            &heurdata->lpitersquot, TRUE, DEFAULT_LPITERSQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/mingap",
         "minimum primal-dual gap for which the heuristic is executed",
         &heurdata->mingap, TRUE, DEFAULT_MINGAP, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useuct",
         "should uct node selection be used at the beginning of the search?",
         &heurdata->useuct, TRUE, DEFAULT_USEUCT, NULL, NULL) );

   return SCIP_OKAY;
}
