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

/**@file   heur_localbranching.c
 * @brief  Local branching heuristic according to Fischetti and Lodi
 * @author Timo Berthold
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_localbranching.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "localbranching"
#define HEUR_DESC             "local branching heuristic by Fischetti and Lodi"
#define HEUR_DISPCHAR         'L'
#define HEUR_PRIORITY         -1102000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NEIGHBORHOODSIZE  18    /* radius of the incumbents neighborhood to be searched                     */
#define DEFAULT_NODESOFS      1000      /* number of nodes added to the contingent of the total nodes               */
#define DEFAULT_MAXNODES      10000     /* maximum number of nodes to regard in the subproblem                      */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which localbranching should at least improve the incumbent     */
#define DEFAULT_MINNODES      1000      /* minimum number of nodes required to start the subproblem                 */
#define DEFAULT_NODESQUOT     0.05      /* contingent of sub problem nodes in relation to original nodes            */
#define DEFAULT_LPLIMFAC      1.5       /* factor by which the limit on the number of LP depends on the node limit  */
#define DEFAULT_NWAITINGNODES 200       /* number of nodes without incumbent change that heuristic should wait      */
#define DEFAULT_USELPROWS     FALSE     /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used    */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip
                                         */
#define DEFAULT_BESTSOLLIMIT   3         /* limit on number of improving incumbent solutions in sub-CIP            */
#define DEFAULT_USEUCT        FALSE     /* should uct node selection be used at the beginning of the search?     */

/* event handler properties */
#define EVENTHDLR_NAME         "Localbranching"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"


#define EXECUTE               0
#define WAITFORNEWSOL         1


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait  */
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes           */
   int                   minnodes;           /**< minimum number of nodes required to start the subproblem             */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem                  */
   SCIP_Longint          usednodes;          /**< amount of nodes local branching used during all calls                */
   SCIP_Real             nodesquot;          /**< contingent of sub problem nodes in relation to original nodes        */
   SCIP_Real             minimprove;         /**< factor by which localbranching should at least improve the incumbent */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   int                   neighborhoodsize;   /**< radius of the incumbent's neighborhood to be searched                */
   int                   callstatus;         /**< current status of localbranching heuristic                           */
   SCIP_SOL*             lastsol;            /**< the last incumbent localbranching used as reference point            */
   int                   curneighborhoodsize;/**< current neighborhoodsize                                             */
   int                   curminnodes;        /**< current minimal number of nodes required to start the subproblem     */
   int                   emptyneighborhoodsize;/**< size of neighborhood that was proven to be empty                   */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?         */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?
                                              */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP            */
   SCIP_Bool             useuct;             /**< should uct node selection be used at the beginning of the search?  */
};


/*
 * Local methods
 */

/** create the extra constraint of local branching and add it to subscip */
static
SCIP_RETCODE addLocalBranchingConstraint(
   SCIP*                 scip,               /**< SCIP data structure of the original problem   */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem         */
   SCIP_VAR**            subvars,            /**< variables of the subproblem                   */
   SCIP_HEURDATA*        heurdata            /**< heuristic's data structure                    */
   )
{
   SCIP_CONS* cons;                        /* local branching constraint to create */
   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;

   int nbinvars;
   int i;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real* consvals;
   char consname[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_localbranchcons", SCIPgetProbName(scip));

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, NULL, NULL, NULL) );
   bestsol = SCIPgetBestSol(scip);
   assert( bestsol != NULL );

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nbinvars) );

   /* set initial left and right hand sides of local branching constraint */
   lhs = (SCIP_Real)heurdata->emptyneighborhoodsize + 1.0;
   rhs = (SCIP_Real)heurdata->curneighborhoodsize;

   /* create the distance (to incumbent) function of the binary variables */
   for( i = 0; i < nbinvars; i++ )
   {
      SCIP_Real solval;

      solval = SCIPgetSolVal(scip, bestsol, vars[i]);
      assert( SCIPisFeasIntegral(scip,solval) );

      /* is variable i  part of the binary support of bestsol? */
      if( SCIPisFeasEQ(scip,solval,1.0) )
      {
         consvals[i] = -1.0;
         rhs -= 1.0;
         lhs -= 1.0;
      }
      else
         consvals[i] = 1.0;
      consvars[i] = subvars[i];
      assert( SCIPvarGetType(consvars[i]) == SCIP_VARTYPE_BINARY );
   }

   /* creates localbranching constraint and adds it to subscip */
   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, consname, nbinvars, consvars, consvals,
         lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   /* free local memory */
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< SCIP data structure  of the original problem      */
   SCIP*                 subscip,            /**< SCIP data structure  of the subproblem            */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< the Localbranching heuristic                      */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< pointer to store, whether new solution was found  */
   )
{
   SCIP_VAR** vars;
   int nvars;
   SCIP_SOL* newsol;
   SCIP_Real* subsolvals;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( subvars != NULL );
   assert( subsol != NULL );

   /* copy the solution */
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
SCIP_DECL_EVENTEXEC(eventExecLocalbranching)
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
      SCIPdebugMsg(scip, "interrupt after  %" SCIP_LONGINT_FORMAT " LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLocalbranching)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurLocalbranching(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLocalbranching)
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
SCIP_DECL_HEURINIT(heurInitLocalbranching)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* with a little abuse we initialize the heurdata as if localbranching would have finished its last step regularly */
   heurdata->callstatus = WAITFORNEWSOL;
   heurdata->lastsol = NULL;
   heurdata->usednodes = 0;
   heurdata->curneighborhoodsize = heurdata->neighborhoodsize;
   heurdata->curminnodes = heurdata->minnodes;
   heurdata->emptyneighborhoodsize = 0;

   return SCIP_OKAY;
}

/** todo setup And Solve Subscip */
static
SCIP_RETCODE setupAndSolveSubscipLocalbranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< the subproblem created by localbranching */
   SCIP_HEUR*            heur,               /**< localbranching heuristic */
   SCIP_Longint          nsubnodes,          /**< nodelimit for subscip */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_VAR** subvars;                       /* subproblem's variables                                */
   SCIP_EVENTHDLR* eventhdlr;                /* event handler for LP events                     */
   SCIP_HEURDATA* heurdata;
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR** vars;

   SCIP_Real cutoff;                         /* objective cutoff for the subproblem                   */
   SCIP_Real upperbound;

   int nvars;
   int i;

   SCIP_Bool success;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
   success = FALSE;

   /* create a problem copy as sub SCIP */
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "localbranching", NULL, NULL, 0, heurdata->uselprows,
         heurdata->copycuts, &success, NULL) );

   SCIPdebugMsg(scip, "Copying SCIP was %ssuccessful.\n", success ? "" : "not ");

   /* if the subproblem could not be created, free memory and return */
   if( !success )
   {
      *result = SCIP_DIDNOTRUN;
      goto TERMINATE;
   }

   /* create event handler for LP events */
   eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecLocalbranching, NULL) );
   if( eventhdlr == NULL )
   {
      /* free hash map */
      SCIPhashmapFree(&varmapfw);

      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   for (i = 0; i < nvars; ++i)
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

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
   heurdata->nodelimit = nsubnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nsubnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", MAX(10, nsubnodes/10)) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsollimit) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
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

   /* activate uct node selection at the top of the tree */
   if( heurdata->useuct && SCIPfindNodesel(subscip, "uct") != NULL && !SCIPisParamFixed(subscip, "nodeselection/uct/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/uct/stdpriority", INT_MAX/2) );
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

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

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

   SCIP_CALL( addLocalBranchingConstraint(scip, subscip, subvars, heurdata) );

   /* add an objective cutoff */
   assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

   upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
   if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
   {
      cutoff = (1-heurdata->minimprove)*SCIPgetUpperbound(scip) + heurdata->minimprove*SCIPgetLowerbound(scip);
   }
   else
   {
      if( SCIPgetUpperbound ( scip ) >= 0 )
         cutoff = ( 1 - heurdata->minimprove ) * SCIPgetUpperbound ( scip );
      else
         cutoff = ( 1 + heurdata->minimprove ) * SCIPgetUpperbound ( scip );
   }
   cutoff = MIN(upperbound, cutoff );
   SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );

   /* catch LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );
   }

   /* solve the subproblem */
   SCIPdebugMsg(scip, "solving local branching subproblem with neighborhoodsize %d and maxnodes %" SCIP_LONGINT_FORMAT "\n",
      heurdata->curneighborhoodsize, nsubnodes);


   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* drop LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   heurdata->usednodes += SCIPgetNNodes(subscip);
   SCIPdebugMsg(scip, "local branching used %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT " nodes\n",
      SCIPgetNNodes(subscip), nsubnodes);

   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      }
      if( success )
      {
         SCIPdebugMsg(scip, "-> accepted solution of value %g\n", SCIPgetSolOrigObj(subscip, subsols[i]));
         *result = SCIP_FOUNDSOL;
      }
   }

   /* check the status of the sub-MIP */
   switch( SCIPgetStatus(subscip) )
   {
   case SCIP_STATUS_OPTIMAL:
   case SCIP_STATUS_BESTSOLLIMIT:
      heurdata->callstatus = WAITFORNEWSOL; /* new solution will immediately be installed at next call */
      SCIPdebugMsg(scip, " -> found new solution\n");
      break;

   case SCIP_STATUS_NODELIMIT:
   case SCIP_STATUS_STALLNODELIMIT:
   case SCIP_STATUS_TOTALNODELIMIT:
      heurdata->callstatus = EXECUTE;
      heurdata->curneighborhoodsize = (heurdata->emptyneighborhoodsize + heurdata->curneighborhoodsize)/2;
      heurdata->curminnodes *= 2;
      SCIPdebugMsg(scip, " -> node limit reached: reduced neighborhood to %d, increased minnodes to %d\n",
         heurdata->curneighborhoodsize, heurdata->curminnodes);
      if( heurdata->curneighborhoodsize <= heurdata->emptyneighborhoodsize )
      {
         heurdata->callstatus = WAITFORNEWSOL;
         SCIPdebugMsg(scip, " -> new neighborhood was already proven to be empty: wait for new solution\n");
      }
      break;

   case SCIP_STATUS_INFEASIBLE:
   case SCIP_STATUS_INFORUNBD:
      heurdata->emptyneighborhoodsize = heurdata->curneighborhoodsize;
      heurdata->curneighborhoodsize += heurdata->curneighborhoodsize/2;
      heurdata->curneighborhoodsize = MAX(heurdata->curneighborhoodsize, heurdata->emptyneighborhoodsize + 2);
      heurdata->callstatus = EXECUTE;
      SCIPdebugMsg(scip, " -> neighborhood is empty: increased neighborhood to %d\n", heurdata->curneighborhoodsize);
      break;

   case SCIP_STATUS_UNKNOWN:
   case SCIP_STATUS_USERINTERRUPT:
   case SCIP_STATUS_TIMELIMIT:
   case SCIP_STATUS_MEMLIMIT:
   case SCIP_STATUS_GAPLIMIT:
   case SCIP_STATUS_SOLLIMIT:
   case SCIP_STATUS_RESTARTLIMIT:
   case SCIP_STATUS_UNBOUNDED:
   default:
      heurdata->callstatus = WAITFORNEWSOL;
      SCIPdebugMsg(scip, " -> unexpected sub-MIP status <%d>: waiting for new solution\n", SCIPgetStatus(subscip));
      break;
   }

 TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLocalbranching)
{  /*lint --e{715}*/
   SCIP_Longint maxnnodes;                   /* maximum number of subnodes                            */
   SCIP_Longint nsubnodes;                   /* nodelimit for subscip                                 */

   SCIP_HEURDATA* heurdata;
   SCIP* subscip;                            /* the subproblem created by localbranching              */

   SCIP_SOL* bestsol;                        /* best solution so far                                  */

   SCIP_Bool success;
   SCIP_RETCODE retcode;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* there should be enough binary variables that a local branching constraint makes sense */
   if( SCIPgetNBinVars(scip) < 2*heurdata->neighborhoodsize )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   /* only call heuristic, if an IP solution is at hand */
   if( SCIPgetNSols(scip) <= 0  )
      return SCIP_OKAY;

   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   /* only call heuristic, if the best solution comes from transformed problem */
   if( SCIPsolIsOriginal(bestsol) )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip, bestsol)  < heurdata->nwaitingnodes)
      return SCIP_OKAY;

   /* only call heuristic, if the best solution does not come from trivial heuristic */
   if( SCIPsolGetHeur(bestsol) != NULL && strcmp(SCIPheurGetName(SCIPsolGetHeur(bestsol)), "trivial") == 0 )
      return SCIP_OKAY;

   /* reset neighborhood and minnodes, if new solution was found */
   if( heurdata->lastsol != bestsol )
   {
      heurdata->curneighborhoodsize = heurdata->neighborhoodsize;
      heurdata->curminnodes = heurdata->minnodes;
      heurdata->emptyneighborhoodsize = 0;
      heurdata->callstatus = EXECUTE;
      heurdata->lastsol = bestsol;
   }

   /* if no new solution was found and local branching also seems to fail, just keep on waiting */
   if( heurdata->callstatus == WAITFORNEWSOL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward local branching if it succeeded often */
   maxnnodes = (SCIP_Longint)(maxnnodes * (1.0 + 2.0*(SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur)+1.0)));
   maxnnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nsubnodes = maxnnodes - heurdata->usednodes;
   nsubnodes = MIN(nsubnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call sub problem solving */
   if( nsubnodes < heurdata->curminnodes )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   /* abort if no time is left or not enough memory to create a copy of SCIP */
   if( !success )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "running localbranching heuristic ...\n");

   SCIP_CALL( SCIPcreate(&subscip) );

   retcode = setupAndSolveSubscipLocalbranching(scip, subscip, heur, nsubnodes, result);

   SCIP_CALL( SCIPfree(&subscip) );

   return retcode;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the localbranching primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLocalbranching(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Localbranching primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLocalbranching, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLocalbranching) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLocalbranching) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLocalbranching) );

   /* add localbranching primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/neighborhoodsize",
         "radius (using Manhattan metric) of the incumbent's neighborhood to be searched",
         &heurdata->neighborhoodsize, FALSE, DEFAULT_NEIGHBORHOODSIZE, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which localbranching should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/bestsollimit",
         "limit on number of improving incumbent solutions in sub-CIP",
         &heurdata->bestsollimit, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useuct",
         "should uct node selection be used at the beginning of the search?",
         &heurdata->useuct, TRUE, DEFAULT_USEUCT, NULL, NULL) );

   return SCIP_OKAY;
}
