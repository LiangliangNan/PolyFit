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

/**@file   heur_ofins.c
 * @brief  OFINS - Objective Function Induced Neighborhood Search - a primal heuristic for reoptimization
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "scip/heur_ofins.h"
#include "scip/scipdefplugins.h"       /* needed for the secondary SCIP instance */
#include "scip/pub_misc.h"

#define HEUR_NAME             "ofins"
#define HEUR_DESC             "primal heuristic for reoptimization, objective function induced neighborhood search"
#define HEUR_DISPCHAR         'A'
#define HEUR_PRIORITY         60000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

/* default values for OFINS-specific plugins */
#define DEFAULT_MAXNODES      5000LL    /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MAXCHGRATE    0.50      /**< maximum percentage of changed objective coefficients */
#define DEFAULT_COPYCUTS      TRUE      /**< if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         *   of the original scip be copied to constraints of the subscip */
#define DEFAULT_MAXCHANGE     0.04      /**< maximal rate of change per coefficient to get fixed */
#define DEFAULT_MINIMPROVE    0.01      /**< factor by which OFINS should at least improve the incumbent */
#define DEFAULT_ADDALLSOLS   FALSE      /**< should all subproblem solutions be added to the original SCIP? */
#define DEFAULT_MINNODES      50LL      /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS      500LL     /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT     0.1       /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_LPLIMFAC      2.0       /**< factor by which the limit on the number of LP depends on the node limit */

/* event handler properties */
#define EVENTHDLR_NAME         "Ofins"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Real             maxchangerate;      /**< maximal rate of changed coefficients in the objective function */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in subproblem? */
   SCIP_Bool             addallsols;         /**< should all subproblem solutions be added to the original SCIP? */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Real             maxchange;          /**< maximal rate of change per coefficient to get fixed */
   SCIP_Real             minimprove;         /**< factor by which OFINS should at least improve the incumbent */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
};

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecOfins)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetNLPs(scip) > heurdata->lplimfac * heurdata->nodelimit )
   {
      SCIPdebugMsg(scip, "interrupt after %" SCIP_LONGINT_FORMAT " LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< RENS heuristic structure                            */
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
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, TRUE, TRUE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/* setup and solve the sub-SCIP */
static
SCIP_RETCODE setupAndSolve(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_HEURDATA*        heurdata,           /**< euristic's private data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_Bool*            chgcoeffs           /**< array of changed coefficients */

   )
{
   SCIP_HASHMAP* varmapfw;
   SCIP_VAR** vars;
   SCIP_VAR** subvars;
   SCIP_EVENTHDLR* eventhdlr;

   SCIP_SOL* sol;
   SCIP_VAR** fixedvars;
   SCIP_Real* fixedvals;
   int nfixedvars;

   int nvars;
   int nintvars;
   int i;

   SCIP_SOL** subsols;
   int nsubsols = 0;

   SCIP_Bool success;
   SCIP_RETCODE retcode;
   SCIP_STATUS status;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);
   assert(chgcoeffs != NULL);

   SCIPdebugMsg(scip, "+---+ Start OFINS heuristic +---+\n");

   /* get variable data */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

   /* get optimal solution of the last iteration */
   sol = SCIPgetReoptLastOptSol(scip);

   /* if the solution is NULL the last problem was infeasible */
   if( sol == NULL )
      return SCIP_OKAY;

   nintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, nvars) );

   /* determine variables to fix in the sub-SCIP */
   nfixedvars = 0;
   for( i = 0; i < nintvars; i++ )
   {
      if( !chgcoeffs[i] )
      {
         fixedvars[nfixedvars] = vars[i];
         fixedvals[nfixedvars] = SCIPgetSolVal(scip, sol, vars[i]);
         ++nfixedvars;
      }
   }

   /* create a problem copy as sub SCIP */
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "ofins", fixedvars, fixedvals, nfixedvars, FALSE,
         FALSE, &success, NULL) );
   assert(success);

   SCIPfreeBufferArrayNull(scip, &fixedvals);
   SCIPfreeBufferArrayNull(scip, &fixedvars);

   /* create event handler for LP events */
   eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecOfins, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   for( i = 0; i < nvars; i++ )
   {
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);
     assert(subvars[i] != NULL);
   }

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* set an objective limit */
   SCIPdebugMsg(scip, "set objective limit of %g to sub-SCIP\n", SCIPgetUpperbound(scip));
   SCIP_CALL( SCIPsetObjlimit(subscip, SCIPgetUpperbound(scip)) );

   SCIPdebugMsg(scip, "OFINS subproblem: %d vars, %d cons\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip));

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
   heurdata->nodelimit = heurdata->maxnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );

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

   /* disable conflict analysis */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", FALSE) );
   }

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* presolve the subproblem */
   retcode = SCIPpresolve(subscip);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while presolving subproblem in %s heuristic; sub-SCIP terminated with code <%d>\n", HEUR_NAME, retcode);

      SCIPABORT(); /*lint --e{527}*/

      /* free */
      SCIPfreeBufferArray(scip, &subvars);
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "%s presolved subproblem: %d vars, %d cons\n", HEUR_NAME, SCIPgetNVars(subscip), SCIPgetNConss(subscip));

   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   /* solve the subproblem */
   SCIPdebugMsg(scip, "solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->maxnodes);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   status = SCIPgetStatus(subscip);

   switch (status) {
   case SCIP_STATUS_INFEASIBLE:
      break;
   case SCIP_STATUS_INFORUNBD:
   case SCIP_STATUS_UNBOUNDED:
   {
      int nsubvars;

      nsubvars = SCIPgetNOrigVars(subscip);

      /* transfer the primal ray from the sub-SCIP to the main SCIP */
      if( SCIPhasPrimalRay(subscip) )
      {
         SCIP_SOL* primalray;

         SCIP_CALL( SCIPcreateSol(scip, &primalray, heur) );

         /* transform the ray into the space of the source scip */
         for( i = 0; i < nsubvars; i++ )
         {
            SCIP_CALL( SCIPsetSolVal(scip, primalray, vars[SCIPvarGetProbindex(subvars[i])],
                  SCIPgetPrimalRayVal(subscip, subvars[i])) );
         }

         SCIPdebug( SCIP_CALL( SCIPprintRay(scip, primalray, 0, FALSE) ); );

         /* update the primal ray of the source scip */
         SCIP_CALL( SCIPupdatePrimalRay(scip, primalray) );
         SCIP_CALL( SCIPfreeSol(scip, &primalray) );

         *result = SCIP_UNBOUNDED;
      }

      break;
   }
   default:
      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      for( i = 0; i < nsubsols && (!success || heurdata->addallsols); i++ )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
         if( success )
            *result = SCIP_FOUNDSOL;
      }
      break;
   } /*lint !e788*/

   SCIPstatisticPrintf("%s statistic: fixed %6.3f integer variables, needed %6.1f seconds, %" SCIP_LONGINT_FORMAT " nodes, solution %10.4f found at node %" SCIP_LONGINT_FORMAT "\n",
      HEUR_NAME, 0.0, SCIPgetSolvingTime(subscip), SCIPgetNNodes(subscip), success ? SCIPgetPrimalbound(scip) : SCIPinfinity(scip),
      nsubsols > 0 ? SCIPsolGetNodenum(SCIPgetBestSol(subscip)) : -1 );

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}

/** main procedure of the OFINS heuristic, creates and solves a sub-SCIP */
static
SCIP_RETCODE applyOfins(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_HEURDATA*        heurdata,           /**< euristic's private data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_Bool*            chgcoeffs           /**< array of changed coefficients */
   )
{
   SCIP* subscip;
   SCIP_RETCODE retcode;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);
   assert(chgcoeffs != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if( !success )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* do not run, if no solution was found */
   if ( SCIPgetReoptLastOptSol(scip) == NULL )
      return SCIP_OKAY;

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   retcode = setupAndSolve(scip, subscip, heur, heurdata, result, nstallnodes, chgcoeffs);

   SCIP_CALL( SCIPfree(&subscip) );

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}




/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyOfins)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurOfins(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeOfins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecOfins)
{/*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   SCIP_Bool* chgcoeffs;
   SCIP_Longint nstallnodes;
   int nchgcoefs;
   int nvars;
   int v;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* only call heuristic, if reoptimization is enabled */
   if( !SCIPisReoptEnabled(scip) )
      return SCIP_OKAY;

   /* only call the heuristic, if we are in run >= 2 */
   if( SCIPgetNReoptRuns(scip) <= 1 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward OFINS if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping OFINS: nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* get variable data and check which coefficient has changed  */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);
   nchgcoefs = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &chgcoeffs, nvars) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real newcoef;
      SCIP_Real oldcoef;
      SCIP_Real newcoefabs;
      SCIP_Real oldcoefabs;
      SCIP_Real frac;

      /* we only want to count variables that are unfixed after the presolving */
      assert(SCIPvarGetStatus(vars[v]) != SCIP_VARSTATUS_ORIGINAL);
      assert(SCIPvarIsActive(vars[v]));

      SCIP_CALL( SCIPgetReoptOldObjCoef(scip, vars[v], SCIPgetNReoptRuns(scip), &newcoef) );
      SCIP_CALL( SCIPgetReoptOldObjCoef(scip, vars[v], SCIPgetNReoptRuns(scip)-1, &oldcoef) );
      newcoefabs = REALABS(newcoef);
      oldcoefabs = REALABS(oldcoef);

      /* if both coefficients are zero nothing has changed */
      if( SCIPisZero(scip, newcoef) && SCIPisZero(scip, oldcoef) )
      {
         frac = 0;
      }
      /* if exactly one coefficient is zero, the other need to be close to zero */
      else if( SCIPisZero(scip, newcoef) || SCIPisZero(scip, oldcoef) )
      {
         assert(SCIPisZero(scip, newcoef) != SCIPisZero(scip, oldcoef));
         if( !SCIPisZero(scip, newcoef) )
            frac = MIN(1, newcoefabs);
         else
            frac = MIN(1, oldcoefabs);
      }
      /* if both coefficients have the same sign we calculate the quotient
       * MIN(newcoefabs, oldcoefabs)/MAX(newcoefabs, oldcoefabs)
       */
      else if( SCIPisPositive(scip, newcoef) == SCIPisPositive(scip, oldcoef) )
      {
         frac = 1.0 - MIN(newcoefabs, oldcoefabs)/MAX(newcoefabs, oldcoefabs);
      }
      /* if both coefficients have a different sign, we set frac = 1 */
      else
      {
         assert((SCIPisPositive(scip, newcoef) && SCIPisNegative(scip, oldcoef))
             || (SCIPisNegative(scip, newcoef) && SCIPisPositive(scip, oldcoef)));

         frac = 1;
      }

      if( frac > heurdata->maxchange )
      {
         chgcoeffs[v] = TRUE;
         nchgcoefs++;
      }
      else
         chgcoeffs[v] = FALSE;
   }

   SCIPdebugMsg(scip, "%d (rate %.4f) changed coefficients\n", nchgcoefs, nchgcoefs/((SCIP_Real)nvars));

   /* we only want to run the heuristic, if there at least 3 changed coefficients.
    * if the number of changed coefficients is 2 the trivialnegation heuristic will construct an
    * optimal solution without solving a MIP.
    */
   if( nchgcoefs < 3 )
      goto TERMINATE;

   /* run the heuristic, if not too many coefficients have changed */
   if( nchgcoefs/((SCIP_Real)nvars) > heurdata->maxchangerate )
      goto TERMINATE;

   SCIP_CALL( applyOfins(scip, heur, heurdata, result, nstallnodes, chgcoeffs) );

  TERMINATE:
   SCIPfreeBufferArray(scip, &chgcoeffs);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the ofins primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurOfins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create ofins primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   assert(heurdata != NULL);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecOfins, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyOfins) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeOfins) );

   /* add ofins primal heuristic parameters */

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxchangerate",
         "maximal rate of changed coefficients",
         &heurdata->maxchangerate, FALSE, DEFAULT_MAXCHGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxchange",
         "maximal rate of change per coefficient to get fixed",
         &heurdata->maxchange, FALSE, DEFAULT_MAXCHANGE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/addallsols",
         "should all subproblem solutions be added to the original SCIP?",
         &heurdata->addallsols, TRUE, DEFAULT_ADDALLSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which RENS should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
