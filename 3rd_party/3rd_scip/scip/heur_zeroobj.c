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

/**@file   heur_zeroobj.c
 * @brief  heuristic that tries to solve the problem without objective. In Gurobi, this heuristic is known as "Hail Mary"
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_zeroobj.h"
#include "scip/cons_linear.h"

#define HEUR_NAME             "zeroobj"
#define HEUR_DESC             "heuristic trying to solve the problem without objective"
#define HEUR_DISPCHAR         'Z'
#define HEUR_PRIORITY         100
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_BEFOREPRESOL
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

/* event handler properties */
#define EVENTHDLR_NAME         "Zeroobj"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"

/* default values for zeroobj-specific plugins */
#define DEFAULT_MAXNODES      1000LL    /* maximum number of nodes to regard in the subproblem                       */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which zeroobj should at least improve the incumbent             */
#define DEFAULT_MINNODES      100LL     /* minimum number of nodes to regard in the subproblem                       */
#define DEFAULT_MAXLPITERS    5000LL    /* maximum number of LP iterations to be performed in the subproblem         */
#define DEFAULT_NODESOFS      100LL     /* number of nodes added to the contingent of the total nodes                */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem             */
#define DEFAULT_ADDALLSOLS    FALSE     /* should all subproblem solutions be added to the original SCIP?            */
#define DEFAULT_ONLYWITHOUTSOL   TRUE   /**< should heuristic only be executed if no primal solution was found, yet? */
#define DEFAULT_USEUCT        FALSE     /* should uct node selection be used at the beginning of the search?     */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          maxlpiters;         /**< maximum number of LP iterations to be performed in the subproblem   */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;          /**< nodes already used by zeroobj in earlier calls                      */
   SCIP_Real             minimprove;         /**< factor by which zeroobj should at least improve the incumbent       */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_Bool             addallsols;         /**< should all subproblem solutions be added to the original SCIP?      */
   SCIP_Bool             onlywithoutsol;     /**< should heuristic only be executed if no primal solution was found, yet? */
   SCIP_Bool             useuct;             /**< should uct node selection be used at the beginning of the search?  */
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
   SCIP_HEUR*            heur,               /**< zeroobj heuristic structure                            */
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

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecZeroobj)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_NODESOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_ITERLIMIT || SCIPgetNLPIterations(scip) >= heurdata->maxlpiters )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}
/* ---------------- Callback methods of primal heuristic ---------------- */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyZeroobj)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurZeroobj(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeZeroobj)
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
SCIP_DECL_HEURINIT(heurInitZeroobj)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->usednodes = 0;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecZeroobj)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_Longint nnodes;                 /* number of stalling nodes for the subproblem */

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward zeroobj if it succeeded often */
   nnodes = (SCIP_Longint)(nnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
   nnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nnodes -= heurdata->usednodes;
   nnodes = MIN(nnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping zeroobj: nnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* do not run zeroobj, if the problem does not have an objective function anyway */
   if( SCIPgetNObjVars(scip) == 0 )
   {
      SCIPdebugMsg(scip, "skipping zeroobj: pure feasibility problem anyway\n");
      return SCIP_OKAY;
   }

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPapplyZeroobj(scip, heur, result, heurdata->minimprove, nnodes) );

   return SCIP_OKAY;
}

/** setup and solve subscip */
static
SCIP_RETCODE setupAndSolveSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Real             minimprove,         /**< factor by which zeroobj should at least improve the incumbent */
   SCIP_Longint          nnodes              /**< node limit for the subproblem */
   )
{
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem             */
   SCIP_Real large;
   SCIP_HASHMAP*         varmapfw;           /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR**            vars;               /* original problem's variables                    */
   SCIP_VAR**            subvars;            /* subproblem's variables                          */
   SCIP_SOL** subsols;
   SCIP_HEURDATA*        heurdata;           /* heuristic's private data structure              */
   SCIP_EVENTHDLR*       eventhdlr;          /* event handler for LP events                     */

   int nsubsols;
   int nvars;                                /* number of original problem's variables          */
   int i;
   SCIP_Bool success;
   SCIP_Bool valid;


   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* different methods to create sub-problem: either copy LP relaxation or the CIP with all constraints */
   valid = FALSE;

   /* copy complete SCIP instance */
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "zeroobj",  TRUE, FALSE, TRUE, &valid) );
   SCIPdebugMsg(scip, "Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");

   /* create event handler for LP events */
   eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecZeroobj, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* determine large value to set variables to */
   large = SCIPinfinity(scip);
   if( !SCIPisInfinity(scip, 0.1 / SCIPfeastol(scip)) )
      large = 0.1 / SCIPfeastol(scip);

   /* get variable image and change to 0.0 in sub-SCIP */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real adjustedbound;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real inf;

      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);
      SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], 0.0) );

      lb = SCIPvarGetLbGlobal(subvars[i]);
      ub = SCIPvarGetUbGlobal(subvars[i]);
      inf = SCIPinfinity(subscip);

      /* adjust infinite bounds in order to avoid that variables with non-zero objective 
       * get fixed to infinite value in zeroobj subproblem
       */
      if( SCIPisInfinity(subscip, ub ) )
      {
         adjustedbound = MAX(large, lb+large);
         adjustedbound = MIN(adjustedbound, inf);
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], adjustedbound) );
      }
      if( SCIPisInfinity(subscip, -lb ) )
      {
         adjustedbound = MIN(-large, ub-large);
         adjustedbound = MAX(adjustedbound, -inf);
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], adjustedbound) );
      }
   }

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

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
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nnodes) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/solutions", 1) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable expensive techniques that merely work on the dual bound */

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );
   if( !SCIPisParamFixed(subscip, "presolving/maxrounds") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "presolving/maxrounds", 50) );
   }

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
   /* use least infeasible branching */
   if( SCIPfindBranchrule(subscip, "leastinf") != NULL && !SCIPisParamFixed(subscip, "branching/leastinf/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/leastinf/priority", INT_MAX/4) );
   }

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
   }

   /* disable feaspump and fracdiving */
   if( !SCIPisParamFixed(subscip, "heuristics/feaspump/freq") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/feaspump/freq", -1) );
   }
   if( !SCIPisParamFixed(subscip, "heuristics/fracdiving/freq") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/fracdiving/freq", -1) );
   }

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* restrict LP iterations */
   SCIP_CALL( SCIPsetLongintParam(subscip, "lp/iterlim", 2*heurdata->maxlpiters / MAX(1,nnodes)) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "lp/rootiterlim", heurdata->maxlpiters) );

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      SCIP_CONS* origobjcons;
#ifndef NDEBUG
      int nobjvars;
      nobjvars = 0;
#endif

      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
      {
         cutoff = (1-minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound(scip) >= 0 )
            cutoff = ( 1 - minimprove ) * SCIPgetUpperbound ( scip );
         else
            cutoff = ( 1 + minimprove ) * SCIPgetUpperbound ( scip );
      }
      cutoff = MIN(upperbound, cutoff);

      SCIP_CALL( SCIPcreateConsLinear(subscip, &origobjcons, "objbound_of_origscip", 0, NULL, NULL, -SCIPinfinity(subscip), cutoff,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      for( i = 0; i < nvars; ++i)
      {
         if( !SCIPisFeasZero(subscip, SCIPvarGetObj(vars[i])) )
         {
            SCIP_CALL( SCIPaddCoefLinear(subscip, origobjcons, subvars[i], SCIPvarGetObj(vars[i])) );
#ifndef NDEBUG
            nobjvars++;
#endif
         }
      }
      SCIP_CALL( SCIPaddCons(subscip, origobjcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &origobjcons) );
      assert(nobjvars == SCIPgetNObjVars(scip));
   }

   /* catch LP events of sub-SCIP */
   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   SCIPdebugMsg(scip, "solving subproblem: nnodes=%" SCIP_LONGINT_FORMAT "\n", nnodes);


   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* drop LP events of sub-SCIP */
   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   success = FALSE;
   for( i = 0; i < nsubsols && (!success || heurdata->addallsols); ++i )
   {
      SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      if( success )
         *result = SCIP_FOUNDSOL;
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */


/** main procedure of the zeroobj heuristic, creates and solves a sub-SCIP */
SCIP_RETCODE SCIPapplyZeroobj(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                            */
   SCIP_RESULT*          result,             /**< result data structure                                               */
   SCIP_Real             minimprove,         /**< factor by which zeroobj should at least improve the incumbent      */
   SCIP_Longint          nnodes              /**< node limit for the subproblem                                       */
   )
{
   SCIP*                 subscip;            /* the subproblem created by zeroobj              */
   SCIP_HEURDATA*        heurdata;           /* heuristic's private data structure              */
   SCIP_Bool             success;
   SCIP_RETCODE          retcode;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   assert(nnodes >= 0);
   assert(0.0 <= minimprove && minimprove <= 1.0);

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic once at the root */
   if( SCIPgetDepth(scip) <= 0 && SCIPheurGetNCalls(heur) > 0 )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only call the heuristic if we do not have an incumbent  */
   if( SCIPgetNSolsFound(scip) > 0 && heurdata->onlywithoutsol )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if( !success )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;


   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   retcode = setupAndSolveSubscip(scip, subscip, heur, result, minimprove, nnodes);

   SCIP_CALL( SCIPfree(&subscip) );

   return retcode;
}


/** creates the zeroobj primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurZeroobj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   heur = NULL;
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecZeroobj, heurdata) );
   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyZeroobj) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeZeroobj) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitZeroobj) );

   /* add zeroobj primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxlpiters",
         "maximum number of LP iterations to be performed in the subproblem",
         &heurdata->maxlpiters, TRUE, DEFAULT_MAXLPITERS, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which zeroobj should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/addallsols",
         "should all subproblem solutions be added to the original SCIP?",
         &heurdata->addallsols, TRUE, DEFAULT_ADDALLSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/onlywithoutsol",
         "should heuristic only be executed if no primal solution was found, yet?",
         &heurdata->onlywithoutsol, TRUE, DEFAULT_ONLYWITHOUTSOL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useuct",
         "should uct node selection be used at the beginning of the search?",
         &heurdata->useuct, TRUE, DEFAULT_USEUCT, NULL, NULL) );

   return SCIP_OKAY;
}
