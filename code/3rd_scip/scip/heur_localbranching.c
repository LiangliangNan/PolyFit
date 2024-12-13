/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_localbranching.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  Local branching heuristic according to Fischetti and Lodi
 * @author Timo Berthold
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/heuristics.h"
#include "scip/heur_localbranching.h"
#include "scip/pub_event.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include <string.h>

#define HEUR_NAME             "localbranching"
#define HEUR_DESC             "local branching heuristic by Fischetti and Lodi"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
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
};


/*
 * Local methods
 */

/** create the extra constraint of local branching and add it to subscip */
static
SCIP_RETCODE addLocalbranchingConstraintAndObjcutoff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< the subproblem created by localbranching */
   SCIP_HEUR*            heur,               /**< the local branching heuristic */
   SCIP_VAR**            subvars             /**< the subproblem variables */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_CONS* cons;                        /* local branching constraint to create */
   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;

   int nbinvars;
   int nconsvars;
   int i;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real* consvals;
   char consname[SCIP_MAXSTRLEN];

   SCIP_Real cutoff;
   SCIP_Real upperbound;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_localbranchcons", SCIPgetProbName(scip));

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, NULL, NULL, NULL) );
   bestsol = SCIPgetBestSol(scip);
   assert( bestsol != NULL );

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nbinvars) );
   nconsvars = 0;

   /* set initial left and right hand sides of local branching constraint */
   lhs = (SCIP_Real)heurdata->emptyneighborhoodsize + 1.0;
   rhs = (SCIP_Real)heurdata->curneighborhoodsize;

   /* create the distance (to incumbent) function of the binary variables */
   for( i = 0; i < nbinvars; i++ )
   {
      SCIP_Real solval;

      if( subvars[i] == NULL )
         continue;

      solval = SCIPgetSolVal(scip, bestsol, vars[i]);
      assert( SCIPisFeasIntegral(scip, solval) );

      /* is variable i  part of the binary support of bestsol? */
      if( SCIPisFeasEQ(scip, solval, 1.0) )
      {
         consvals[nconsvars] = -1.0;
         rhs -= 1.0;
         lhs -= 1.0;
      }
      else
         consvals[nconsvars] = 1.0;

      consvars[nconsvars] = subvars[i];
      assert( SCIPvarGetType(consvars[nconsvars]) == SCIP_VARTYPE_BINARY );

      ++nconsvars;
   }

   /* creates localbranching constraint and adds it to subscip */
   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, consname, nconsvars, consvars, consvals,
         lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

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

   /* free local memory */
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/** event handler execution callback to interrupt the solution process */
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

/** setup And solve local branching subscip */
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

   heurdata->nodelimit = nsubnodes;
   SCIP_CALL( SCIPsetCommonSubscipParams(scip, subscip, nsubnodes, MAX(10, nsubnodes/10), heurdata->bestsollimit) );

   /* adds the local branching constraint and the objective cutoff to the auxiliary problem */
   SCIP_CALL( addLocalbranchingConstraintAndObjcutoff(scip, subscip, heur, subvars) );

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

   /* checks the solutions of the sub SCIP and adds them to the main SCIP if feasible */
   SCIP_CALL( SCIPtranslateSubSols(scip, subscip, heur, subvars, &success, NULL) );

   if( success )
      *result = SCIP_FOUNDSOL;

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
   case SCIP_STATUS_TERMINATE:
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

   return SCIP_OKAY;
}
