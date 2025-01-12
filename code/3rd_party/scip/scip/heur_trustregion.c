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

/**@file   heur_trustregion.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  Large neighborhood search heuristic for Benders' decomposition based on trust region methods
 * @author Stephen J. Maher
 *
 * The Trust Region heuristic draws upon trust region methods for solving optimization problems, especially in the
 * context of Benders' decomposition. This heuristic has been developed to improve the heuristic performance of the
 * Benders' decomposition algorithm within SCIP.
 *
 * The Trust Region heuristic copies the original SCIP instance and adds a constraint to penalize changes from the
 * incumbent solution. Consider a problem that includes a set of binary variables \f$\mathcal{B}\f$. Given a feasible
 * solution \f$\hat{x}\f$ to the original problem, we define the set \f$\mathcal{B}^{+}\f$ as the index set for the
 * binary variables that are 1 in the input solution and \f$\mathcal{B}^{-}\f$ as the index set for binary variables
 * that are 0. The trust region constraint, which is added to the sub-SCIP, is given by
 *
 * \f[
 *    \sum_{i \in \mathcal{B}^{+}}(1 - x_{i}) + \sum_{i \in \mathcal{B}^{-}}x_{i} \le \theta
 * \f]
 *
 * The variable \f$\theta\f$ measure the distance, in terms of the binary variables, of candidate solutions to the input
 * solution.
 *
 * In addition, an upper bounding constraint is explicitly added to enforce a minimum improvement from the heuristic,
 * given by \f$f(x) \le f(\hat{x}) - \epsilon\f$. The parameter \f$\epsilon \ge 0\f$ denotes the minimum improvement
 * that must be achieved by the heuristic.
 *
 * The objective function is then modified to \f$f(x) + M\theta\f$, where \f$M\f$ is a parameter for penalizing the
 * distance of solutions from the input solution \f$\hat{x}\f$.
 *
 * If a new incumbent solution is found by this heuristic, then the Trust Region heuristic is immediately
 * re-executed with this new incumbent solution.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/heuristics.h"
#include "scip/heur_trustregion.h"
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
#include "scip/scip_var.h"
#include <string.h>

#define HEUR_NAME             "trustregion"
#define HEUR_DESC             "LNS heuristic for Benders' decomposition based on trust region methods"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY         -1102010
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MINBINVARS    10        /**< the minimum number of binary variables necessary to run the heuristic */
#define DEFAULT_NODESOFS      1000      /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_MAXNODES      10000     /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINNODES      100       /**< minimum number of nodes required to start the subproblem */
#define DEFAULT_NODESQUOT     0.05      /**< contingent of sub problem nodes in relation to original nodes */
#define DEFAULT_LPLIMFAC      1.5       /**< factor by which the limit on the number of LP depends on the node limit */
#define DEFAULT_NWAITINGNODES 1         /**< number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_USELPROWS     FALSE     /**< should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE      /**< if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip */
#define DEFAULT_BESTSOLLIMIT   3         /**< limit on number of improving incumbent solutions in sub-CIP */

#define DEFAULT_VIOLPENALTY   100.0     /**< the penalty for violating the trust region */
#define DEFAULT_OBJMINIMPROVE 1e-2      /**< the minimum absolute improvement in the objective function value */

/* event handler properties */
#define EVENTHDLR_NAME         "Trustregion"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"


#define EXECUTE               0
#define WAITFORNEWSOL         1


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             lastsol;            /**< the last incumbent trustregion used as reference point */
   SCIP_Longint          usednodes;          /**< amount of nodes trust region used during all calls */
   SCIP_Real             nodesquot;          /**< contingent of sub problem nodes in relation to original nodes */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   SCIP_Real             violpenalty;        /**< the penalty for violating the trust region */
   SCIP_Real             objminimprove;      /**< the minimum absolute improvement in the objective function value */
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   int                   minnodes;           /**< minimum number of nodes required to start the subproblem */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   int                   minbinvars;         /**< minimum number of binary variables necessary to run the heuristic */
   int                   callstatus;         /**< current status of trustregion heuristic */
   int                   curminnodes;        /**< current minimal number of nodes required to start the subproblem */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows? */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem? */
};


/*
 * Local methods
 */

/** create the extra constraint of trust region and add it to \p subscip */
static
SCIP_RETCODE addTrustRegionConstraints(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem */
   SCIP_VAR**            subvars,            /**< variables of the subproblem */
   SCIP_HEURDATA*        heurdata            /**< heuristic's data structure */
   )
{
   SCIP_CONS* cons;                        /* trust region constraint to create */
   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;

   int nvars;
   int nbinvars;
   int nconsvars;
   int i;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real* consvals;
   char name[SCIP_MAXSTRLEN];

   /* adding the neighborhood constraint for the trust region heuristic */
   SCIP_CALL( SCIPaddTrustregionNeighborhoodConstraint(scip, subscip, subvars, heurdata->violpenalty) );

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   bestsol = SCIPgetBestSol(scip);
   assert( bestsol != NULL );

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars + 1) );
   nconsvars = 0;

   /* create the upper bounding constraint. An absolute minimum improvement is used for this heuristic. This is
    * different to other LNS heuristics, where a relative improvement is used. The absolute improvement tries to take
    * into account problem specific information that is available to the user, such as a minimum step in the objective
    * limit if the objective function is integer
    */
   lhs = -SCIPinfinity(subscip);
   rhs = SCIPgetSolTransObj(scip, bestsol) - heurdata->objminimprove;

   /* if the objective function is integer, then the floor of the RHS is taken */
   if( SCIPisObjIntegral(scip) )
      rhs = SCIPfeasFloor(scip, rhs);

   /* adding the coefficients to the upper bounding constraint */
   for( i = 0; i < nvars; i++ )
   {
      if( subvars[i] == NULL )
         continue;
      consvals[nconsvars] = SCIPvarGetObj(subvars[i]);
      consvars[nconsvars] = subvars[i];
      ++nconsvars;
   }

   /* creates trustregion constraint and adds it to subscip */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_upperboundcons", SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, name, nconsvars, consvars, consvals,
         lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   /* free local memory */
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/** event handler execution callback to interrupt the solution process */
static
SCIP_DECL_EVENTEXEC(eventExecTrustregion)
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
SCIP_DECL_HEURCOPY(heurCopyTrustregion)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTrustregion(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTrustregion)
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
SCIP_DECL_HEURINIT(heurInitTrustregion)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* with a little abuse we initialize the heurdata as if trustregion would have finished its last step regularly */
   heurdata->callstatus = WAITFORNEWSOL;
   heurdata->lastsol = NULL;
   heurdata->usednodes = 0;
   heurdata->curminnodes = heurdata->minnodes;

   return SCIP_OKAY;
}

/** sets up and solves the sub SCIP for the Trust Region heuristic */
static
SCIP_RETCODE setupAndSolveSubscipTrustregion(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< the subproblem created by trustregion */
   SCIP_HEUR*            heur,               /**< trustregion heuristic */
   SCIP_Longint          nsubnodes,          /**< nodelimit for subscip */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_VAR** subvars;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_HEURDATA* heurdata;
   SCIP_HASHMAP* varmapfw;
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
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "trustregion", NULL, NULL, 0, heurdata->uselprows,
         heurdata->copycuts, &success, NULL) );

   SCIPdebugMsg(scip, "Copying SCIP was %s successful.\n", success ? "" : "not ");

   /* if the subproblem could not be created, free memory and return */
   if( !success )
   {
      *result = SCIP_DIDNOTRUN;
      goto TERMINATE;
   }

   /* create event handler for LP events */
   eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecTrustregion, NULL) );
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

   SCIP_CALL( addTrustRegionConstraints(scip, subscip, subvars, heurdata) );

   /* catch LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );
   }

   /* solve the subproblem */
   SCIPdebugMsg(scip, "solving trust region subproblem with maxnodes %" SCIP_LONGINT_FORMAT "\n", nsubnodes);

   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/trysol/priority", 100000) );

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
   SCIPdebugMsg(scip, "trust region used %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT " nodes\n",
      SCIPgetNNodes(subscip), nsubnodes);

   /* checks the solutions of the sub SCIP and adds them to the main SCIP if feasible */
   SCIP_CALL( SCIPtranslateSubSols(scip, subscip, heur, subvars, &success, NULL) );

   if( success )
      *result = SCIP_FOUNDSOL;

   /* checking the status of the subscip */
   heurdata->callstatus = WAITFORNEWSOL;
   if( SCIPgetStatus(subscip) == SCIP_STATUS_NODELIMIT || SCIPgetStatus(subscip) == SCIP_STATUS_STALLNODELIMIT
      || SCIPgetStatus(subscip) == SCIP_STATUS_TOTALNODELIMIT )
   {
      heurdata->callstatus = EXECUTE;
      heurdata->curminnodes *= 2;
   }

 TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTrustregion)
{  /*lint --e{715}*/
   SCIP_Longint maxnnodes;
   SCIP_Longint nsubnodes;

   SCIP_HEURDATA* heurdata;
   SCIP* subscip;

   SCIP_SOL* bestsol;

   SCIP_Bool success;
   SCIP_RETCODE retcode;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* there should be enough binary variables that a trust region constraint makes sense */
   if( SCIPgetNBinVars(scip) < heurdata->minbinvars )
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

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward trust region if it found solutions often.
    * In this case, the trust region heuristic is designed for Benders' decomposition and solutions found may not be
    * added by this heuristic but by trysol. So we don't reward finding best solutions, but finding any solution. */
   maxnnodes = (SCIP_Longint)(maxnnodes * (1.0 + 2.0*(SCIPheurGetNSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur)+1.0)));
   maxnnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
   maxnnodes += heurdata->nodesofs;

   *result = SCIP_DIDNOTRUN;

   /* we continue to execute the trust region heuristic until no new best solution is found */
   do
   {
      SCIP_RESULT heurresult;

      /* storing the best solution again since it is needed for the execution loop */
      bestsol = SCIPgetBestSol(scip);

      /* reset minnodes if new solution was found */
      if( heurdata->lastsol != bestsol )
      {
         heurdata->curminnodes = heurdata->minnodes;
         heurdata->callstatus = EXECUTE;
         heurdata->lastsol = bestsol;
      }

      /* if no new solution was found and trust region also seems to fail, just keep on waiting */
      if( heurdata->callstatus == WAITFORNEWSOL )
         return SCIP_OKAY;

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

      /* abort if no time is left or there is not enough memory to create a copy of SCIP */
      if( !success )
         return SCIP_OKAY;

      heurresult = SCIP_DIDNOTFIND;

      SCIPdebugMsg(scip, "running trust region heuristic ...\n");

      SCIP_CALL( SCIPcreate(&subscip) );

      retcode = setupAndSolveSubscipTrustregion(scip, subscip, heur, nsubnodes, &heurresult);

      SCIP_CALL( SCIPfree(&subscip) );

      /* if the result is FOUNDSOL, this means that a solution was found during a previous execution of the heuristic.
       * So the heuristic result should only be updated if the result is not FOUNDSOL.
       */
      if( *result != SCIP_FOUNDSOL )
         *result = heurresult;
   }
   while( bestsol != SCIPgetBestSol(scip) && retcode == SCIP_OKAY );

   return retcode;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the trustregion primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTrustregion(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Trustregion primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTrustregion, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTrustregion) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTrustregion) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitTrustregion) );

   /* add trustregion primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minbinvars",
         "the number of binary variables necessary to run the heuristic",
         &heurdata->minbinvars, FALSE, DEFAULT_MINBINVARS, 1, INT_MAX, NULL, NULL) );

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

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/bestsollimit",
         "limit on number of improving incumbent solutions in sub-CIP",
         &heurdata->bestsollimit, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/violpenalty",
         "the penalty for each change in the binary variables from the candidate solution",
         &heurdata->violpenalty, FALSE, DEFAULT_VIOLPENALTY, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/objminimprove",
         "the minimum absolute improvement in the objective function value",
         &heurdata->objminimprove, FALSE, DEFAULT_OBJMINIMPROVE, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
