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

/**@file   heur_rens.c
 * @brief  LNS heuristic that finds the optimal rounding to a given point
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/heur_rens.h"
#include "scip/scipdefplugins.h"       /* needed for the secondary SCIP instance */
#include "scip/cons_linear.h"          /* needed if the LP relaxation gets copied into linear constraints */
#include "scip/pub_misc.h"

/* default values for standard parameters that every primal heuristic has in SCIP */
#define HEUR_NAME             "rens"
#define HEUR_DESC             "LNS exploring fractional neighborhood of relaxation's optimum"
#define HEUR_DISPCHAR         'E'
#define HEUR_PRIORITY         -1100000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

/* default values for RENS-specific plugins */
#define DEFAULT_BINARYBOUNDS  TRUE      /* should general integers get binary bounds [floor(.),ceil(.)] ?      */
#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.5       /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which RENS should at least improve the incumbent          */
#define DEFAULT_MINNODES      50LL      /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_NODESOFS      500LL     /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_LPLIMFAC      2.0       /* factor by which the limit on the number of LP depends on the node limit  */
#define DEFAULT_STARTSOL      'l'       /* solution that is used for fixing values                             */
#define STARTSOL_CHOICES      "nl"      /* possible values for startsol ('l'p relaxation, 'n'lp relaxation)    */
#define DEFAULT_USELPROWS     FALSE     /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip
                                         */
#define DEFAULT_EXTRATIME    FALSE      /* should the RENS sub-CIP get its own full time limit? This is only
                                         * implemented for testing and not recommended to be used!
                                         */
#define DEFAULT_ADDALLSOLS   FALSE      /* should all subproblem solutions be added to the original SCIP?       */

#define DEFAULT_FULLSCALE    FALSE      /* should the RENS sub-CIP be solved with full-scale SCIP settings, including
                                         * techniques that merely work on the dual bound, e.g., cuts?  This is only
                                         * implemented for testing and not recommended to be used!
                                         */
#define DEFAULT_BESTSOLLIMIT   -1       /* limit on number of improving incumbent solutions in sub-CIP            */
#define DEFAULT_USEUCT         FALSE    /* should uct node selection be used at the beginning of the search?     */

/* event handler properties */
#define EVENTHDLR_NAME         "Rens"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;          /**< nodes already used by RENS in earlier calls                         */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove;         /**< factor by which RENS should at least improve the incumbent          */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   char                  startsol;           /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   SCIP_Bool             binarybounds;       /**< should general integers get binary bounds [floor(.),ceil(.)] ?      */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?        */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem? */
   SCIP_Bool             extratime;          /**< should the RENS sub-CIP get its own full time limit? This is only
                                              *   implemented for testing and not recommended to be used! */
   SCIP_Bool             addallsols;         /**< should all subproblem solutions be added to the original SCIP? */
   SCIP_Bool             fullscale;          /**< should the RENS sub-CIP be solved with full-scale SCIP settings,
                                              * including techniques that merely work on the dual bound, e.g., cuts?
                                              * This is only implemented for testing and not recommended to be used! */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP            */
   SCIP_Bool             useuct;             /**< should uct node selection be used at the beginning of the search?  */
};


/*
 * Local methods
 */

/** compute the number of initial fixings and check whether the fixing rate exceeds the minimum fixing rate */
static
SCIP_RETCODE computeFixingrate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            fixedvars,          /**< array to store source SCIP variables whose copies should be fixed in the sub-SCIP */
   SCIP_Real*            fixedvals,          /**< array to store solution values for variable fixing */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int                   fixedvarssize,      /**< size of the arrays to store fixing variables */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed */
   char*                 startsol,           /**< pointer to solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   SCIP_Real*            fixingrate,         /**< percentage of integers that get actually fixed */
   SCIP_Bool*            success             /**< pointer to store whether minimum fixingrate is exceeded */
   )
{
   SCIP_VAR** vars;
   int nintvars;
   int nbinvars;
   int i;

   assert(fixedvars != NULL);
   assert(fixedvals != NULL);
   assert(nfixedvars != NULL);

   *fixingrate = 1.0;
   *success = FALSE;

   /* if there is no NLP relaxation available (e.g., because the presolved problem is linear), use LP relaxation */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "no NLP present, use LP relaxation instead\n");
      (*startsol) = 'l';
   }

   /* get required variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );
   assert(fixedvarssize >= nbinvars + nintvars);
   (*nfixedvars) = 0;

   /* try to solve NLP relaxation */
   if( (*startsol) == 'n' )
   {
      SCIP_NLPSOLSTAT stat;
      SCIPdebug( int nlpverblevel; )

      /* only call this function if NLP relaxation is available */
      assert(SCIPisNLPConstructed(scip));

      /* activate NLP solver output if we are in SCIP's debug mode */
      SCIPdebug( SCIP_CALL( SCIPgetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, &nlpverblevel) ) );
      SCIPdebug( SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, MAX(1,nlpverblevel)) ) );

      SCIPdebugMsg(scip, "try to solve NLP relaxation to obtain fixing values\n");

      /* set starting point to LP solution */
      SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );

      /* solve NLP relaxation */
      SCIP_CALL( SCIPsolveNLP(scip) );

      /* get solution status of NLP solver */
      stat = SCIPgetNLPSolstat(scip);
      *success = (stat == SCIP_NLPSOLSTAT_GLOBOPT) || (stat == SCIP_NLPSOLSTAT_LOCOPT) || stat == (SCIP_NLPSOLSTAT_FEASIBLE);
      SCIPdebugMsg(scip, "solving NLP relaxation was %s successful (stat=%d)\n", *success ? "" : "not", stat);

      /* reset NLP verblevel to the value it had before */
      SCIPdebug( SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, nlpverblevel) ) );

      /* it the NLP was not successfully solved we stop the heuristic right away */
      if( !(*success) )
         return SCIP_OKAY;
   }
   else
   {
      assert(*startsol == 'l');
   }

   /* count the number of variables with integral solution values in the current NLP or LP solution */
   for( i = 0; i < nbinvars + nintvars; ++i )
   {
      SCIP_Real solval;

      /* get solution value in the relaxation in question */
      solval = (*startsol == 'l') ? SCIPvarGetLPSol(vars[i]) : SCIPvarGetNLPSol(vars[i]);

      /* append variable to the buffer storage for integer variables with integer solution values */
      if( SCIPisFeasIntegral(scip, solval) )
      {
         /* fix variables to current LP/NLP solution if it is integral,
          * use exact integral value, if the variable is only integral within numerical tolerances
          */
         solval = SCIPfloor(scip, solval+0.5);
         fixedvars[(*nfixedvars)] = vars[i];
         fixedvals[(*nfixedvars)] = solval;
         (*nfixedvars)++;
      }
   }

   /* abort, if all integer variables were fixed (which should not happen for MIP),
    * but frequently happens for MINLPs using an LP relaxation
    */
   if( (*nfixedvars) == nbinvars + nintvars )
      return SCIP_OKAY;

   *fixingrate = (*nfixedvars) / (SCIP_Real)(MAX(nbinvars + nintvars, 1));

   /* abort, if the amount of fixed variables is insufficient */
   if( *fixingrate < minfixingrate )
      return SCIP_OKAY;

   *success = TRUE;
   return SCIP_OKAY;
}

/** fixes bounds of unfixed integer variables to binary bounds */
static
SCIP_RETCODE restrictToBinaryBounds(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                              */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                                     */
   char                  startsol            /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   )
{
   SCIP_VAR** vars;                          /* original SCIP variables */

   int nbinvars;
   int nintvars;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);

   assert(startsol == 'l' || startsol == 'n');

   /* get required variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* change bounds of integer variables of the subproblem */
   for( i = nbinvars; i < nbinvars + nintvars; i++ )
   {
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;

      /* get the current LP/NLP solution for each variable */
      if( startsol == 'l')
         solval = SCIPvarGetLPSol(vars[i]);
      else
         solval = SCIPvarGetNLPSol(vars[i]);

      /* restrict bounds to nearest integers if the solution value is not already integer */
      if( !SCIPisFeasIntegral(scip, solval) )
      {
         lb = SCIPfeasFloor(scip,solval);
         ub = SCIPfeasCeil(scip,solval);

         /* perform the bound change */
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], lb) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], ub) );
      }
      else
      {
         /* the variable bounds should be already fixed to this solution value */
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(subvars[i]), SCIPfloor(scip, solval+0.5)));
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(subvars[i]), SCIPfloor(scip, solval+0.5)));
      }
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
SCIP_DECL_EVENTEXEC(eventExecRens)
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

/** setup and solve the RENS sub-SCIP */
static
SCIP_RETCODE setupAndSolveSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub SCIP data structure */
   SCIP_RESULT*          result,             /**< result pointer */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_VAR**            fixedvars,          /**< array of variables that should be fixed */
   SCIP_Real*            fixedvals,          /**< array of fixing values */
   int                   nfixedvars,         /**< number of variables that should be fixed */
   SCIP_Real             intfixingrate,      /**< percentage of integer variables fixed           */
   SCIP_Real             minfixingrate,      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove,         /**< factor by which RENS should at least improve the incumbent          */
   SCIP_Longint          maxnodes,           /**< maximum number of  nodes for the subproblem                         */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                         */
   char                  startsol,           /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)]?       */
   SCIP_Bool             uselprows           /**< should subproblem be created out of the rows in the LP rows? */
   )
{
   SCIP_VAR** vars;                          /* original problem's variables                    */
   SCIP_VAR** subvars;                       /* subproblem's variables                          */
   SCIP_HEURDATA* heurdata;                  /* heuristic data */
   SCIP_EVENTHDLR* eventhdlr;                /* event handler for LP events                     */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem             */
   SCIP_Real allfixingrate;                  /* percentage of all variables fixed               */
   SCIP_Bool success;
   int i;
   int nvars;              /**< number of original problem's variables */
   SCIP_RETCODE retcode;

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

   /* create a problem copy as sub SCIP */
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "rens", fixedvars, fixedvals, nfixedvars, uselprows,
         heurdata->copycuts, &success, NULL) );

   eventhdlr = NULL;
   /* create event handler for LP events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecRens, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* copy subproblem variables into the same order as the source SCIP variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* restrict the integer variables to binary bounds  */
   if( binarybounds )
   {
      SCIP_CALL( restrictToBinaryBounds(scip, subscip, subvars, startsol) );
   }

   SCIPdebugMsg(scip, "RENS subproblem: %d vars, %d cons\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip));

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
   heurdata->nodelimit = maxnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", maxnodes) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsollimit) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable expensive techniques that merely work on the dual bound */
   if( !heurdata->fullscale )
   {
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
   }

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
      {
         cutoff = (1 - minimprove) * SCIPgetUpperbound(scip)
                       + minimprove * SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound(scip) >= 0 )
            cutoff = (1 - minimprove) * SCIPgetUpperbound(scip);
         else
            cutoff = (1 + minimprove) * SCIPgetUpperbound(scip);
      }
      cutoff = MIN(upperbound, cutoff);
      SCIP_CALL(SCIPsetObjlimit(subscip, cutoff));
   }

   /* presolve the subproblem */
   retcode = SCIPpresolve(subscip);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while presolving subproblem in RENS heuristic; sub-SCIP terminated with code <%d>\n", retcode);
      SCIPABORT();  /*lint --e{527}*/

      /* free sub problem data */
      SCIPfreeBufferArray(scip, &subvars);

      return retcode;
   }

   SCIPdebugMsg(scip, "RENS presolved subproblem: %d vars, %d cons, success=%u\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip), success);

   allfixingrate = (SCIPgetNOrigVars(subscip) - SCIPgetNVars(subscip)) / (SCIP_Real)SCIPgetNOrigVars(subscip);

   /* additional variables added in presolving may lead to the subSCIP having more variables than the original */
   allfixingrate = MAX(allfixingrate, 0.0);

   /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
    * to ensure that not only the MIP but also the LP relaxation is easy enough
    */
   if( allfixingrate >= minfixingrate / 2.0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      /* catch LP events of sub-SCIP */
      assert(eventhdlr != NULL);
      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

      /* solve the subproblem */
      SCIPdebugMsg(scip, "solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, maxnodes);
      retcode = SCIPsolve(subscip);

      /* drop LP events of sub-SCIP */
      SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

      /* errors in solving the subproblem should not kill the overall solving process;
       * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
      if( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving subproblem in RENS heuristic; sub-SCIP terminated with code <%d>\n", retcode);
         SCIPABORT();

         /* free sub problem data */
         SCIPfreeBufferArray(scip, &subvars);

         return retcode;
      }
      else
      {
         /* transfer variable statistics from sub-SCIP */
         SCIP_CALL( SCIPmergeVariableStatistics(subscip, scip, subvars, vars, nvars) );
      }

      /* print solving statistics of subproblem if we are in SCIP's debug mode */
      SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

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

      SCIPstatisticPrintf("RENS statistic: fixed %6.3f integer variables, %6.3f all variables, needed %6.1f seconds, %" SCIP_LONGINT_FORMAT " nodes, solution %10.4f found at node %" SCIP_LONGINT_FORMAT "\n",
            intfixingrate, allfixingrate, SCIPgetSolvingTime(subscip), SCIPgetNNodes(subscip), success ? SCIPgetPrimalbound(scip) : SCIPinfinity(scip),
                  nsubsols > 0 ? SCIPsolGetNodenum(SCIPgetBestSol(subscip)) : -1 );
   }
   else
   {
      SCIPstatisticPrintf("RENS statistic: fixed only %6.3f integer variables, %6.3f all variables --> abort \n", intfixingrate, allfixingrate);
   }

   /* free sub problem data */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}

/* ---------------- external methods of RENS heuristic ---------------- */

/** main procedure of the RENS heuristic, creates and solves a sub-SCIP */
SCIP_RETCODE SCIPapplyRens(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                            */
   SCIP_RESULT*          result,             /**< result data structure                                               */
   SCIP_Real             minfixingrate,      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove,         /**< factor by which RENS should at least improve the incumbent          */
   SCIP_Longint          maxnodes,           /**< maximum number of  nodes for the subproblem                         */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                         */
   char                  startsol,           /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)]?       */
   SCIP_Bool             uselprows           /**< should subproblem be created out of the rows in the LP rows?        */
   )
{
   SCIP* subscip;                            /* the subproblem created by RENS                  */

   SCIP_Real intfixingrate;                  /* percentage of integer variables fixed           */

   SCIP_VAR** fixedvars;
   SCIP_Real* fixedvals;
   int nfixedvars;
   int fixedvarssize;
   int nbinvars;
   int nintvars;

   SCIP_Bool success;
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   assert(maxnodes >= 0);
   assert(nstallnodes >= 0);

   assert(0.0 <= minfixingrate && minfixingrate <= 1.0);
   assert(0.0 <= minimprove && minimprove <= 1.0);
   assert(startsol == 'l' || startsol == 'n');

   *result = SCIP_DIDNOTRUN;

   nbinvars = SCIPgetNBinVars(scip);
   nintvars = SCIPgetNIntVars(scip);

   /* allocate buffer storage to keep fixings for the variables in the sub SCIP */
   fixedvarssize = nbinvars + nintvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, fixedvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, fixedvarssize) );
   nfixedvars = 0;

   /* compute the number of initial fixings and check if the fixing rate exceeds the minimum fixing rate */
   SCIP_CALL( computeFixingrate(scip, fixedvars, fixedvals, &nfixedvars, fixedvarssize, minfixingrate, &startsol, &intfixingrate, &success) );

   if( !success )
   {
      SCIPstatisticPrintf("RENS statistic: fixed only %5.2f integer variables --> abort \n", intfixingrate);
      goto TERMINATE;
   }

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if( !success )
      goto TERMINATE;

   *result = SCIP_DIDNOTFIND;

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   retcode = setupAndSolveSubscip(scip, subscip, result, heur, fixedvars, fixedvals, nfixedvars, intfixingrate, minfixingrate, minimprove, maxnodes, nstallnodes, startsol, binarybounds, uselprows);

   SCIP_CALL( SCIPfree(&subscip) );

   SCIP_CALL( retcode );

TERMINATE:
   /* free buffer storage for variable fixings */
   SCIPfreeBufferArray(scip, &fixedvals);
   SCIPfreeBufferArray(scip, &fixedvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRens)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRens(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRens)
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
SCIP_DECL_HEURINIT(heurInitRens)
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
SCIP_DECL_HEUREXEC(heurExecRens)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );
   assert( SCIPhasCurrentNodeLP(scip) );

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* only call heuristic, if an optimal LP solution is at hand */
   if( heurdata->startsol == 'l' && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( heurdata->startsol == 'l' && SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* only continue with some fractional variables */
   if( heurdata->startsol == 'l' && SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* do not proceed, when we should use the NLP relaxation, but there is no NLP solver included in SCIP */
   if( heurdata->startsol == 'n' && SCIPgetNNlpis(scip) == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward RENS if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping RENS: nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   if( SCIPisStopped(scip) && !heurdata->extratime )
      return SCIP_OKAY;

   SCIP_CALL( SCIPapplyRens(scip, heur, result, heurdata->minfixingrate, heurdata->minimprove,
         heurdata->maxnodes, nstallnodes, heurdata->startsol, heurdata->binarybounds, heurdata->uselprows) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the rens primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRens(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Rens primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRens, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRens) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRens) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRens) );

   /* add rens primal heuristic parameters */

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "minimum percentage of integer variables that have to be fixable",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which RENS should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/startsol",
         "solution that is used for fixing values ('l'p relaxation, 'n'lp relaxation)",
         &heurdata->startsol, FALSE, DEFAULT_STARTSOL, STARTSOL_CHOICES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/binarybounds",
         "should general integers get binary bounds [floor(.),ceil(.)] ?",
         &heurdata->binarybounds, TRUE, DEFAULT_BINARYBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/extratime",
         "should the RENS sub-CIP get its own full time limit? This is only for tesing and not recommended!",
         &heurdata->extratime, TRUE, DEFAULT_EXTRATIME, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/addallsols",
         "should all subproblem solutions be added to the original SCIP?",
         &heurdata->addallsols, TRUE, DEFAULT_ADDALLSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/fullscale",
         "should the RENS sub-CIP be solved with cuts, conflicts, strong branching,... This is only for tesing and not recommended!",
         &heurdata->fullscale, TRUE, DEFAULT_FULLSCALE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/bestsollimit",
         "limit on number of improving incumbent solutions in sub-CIP",
         &heurdata->bestsollimit, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useuct",
         "should uct node selection be used at the beginning of the search?",
         &heurdata->useuct, TRUE, DEFAULT_USEUCT, NULL, NULL) );

   return SCIP_OKAY;
}
