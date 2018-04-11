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

/**@file   heur_mutation.c
 * @brief  LNS heuristic that tries to randomly mutate the incumbent solution
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/heur_mutation.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "mutation"
#define HEUR_DESC             "mutation heuristic randomly fixing variables"
#define HEUR_DISPCHAR         'M'
#define HEUR_PRIORITY         -1103000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          8
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NODESOFS      500           /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_MAXNODES      5000          /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01          /* factor by which Mutation should at least improve the incumbent      */
#define DEFAULT_MINNODES      500           /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.8           /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_NODESQUOT     0.1           /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_NWAITINGNODES 200           /* number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_USELPROWS     FALSE         /* should subproblem be created out of the rows in the LP rows,
                                             * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE          /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the
                                             * cutpool of the original scip be copied to constraints of the subscip */
#define DEFAULT_BESTSOLLIMIT   -1           /* limit on number of improving incumbent solutions in sub-CIP            */
#define DEFAULT_USEUCT         FALSE        /* should uct node selection be used at the beginning of the search?     */
#define DEFAULT_RANDSEED       19           /* initial random seed */
/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   int                   minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed       */
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   SCIP_Real             minimprove;         /**< factor by which Mutation should at least improve the incumbent      */
   SCIP_Longint          usednodes;          /**< nodes already used by Mutation in earlier calls                     */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator                              */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?        */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?
                                              */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP            */
   SCIP_Bool             useuct;             /**< should uct node selection be used at the beginning of the search?  */
};


/*
 * Local methods
 */

/** determine variables and values which should be fixed in the mutation subproblem */
static
SCIP_RETCODE determineVariableFixings(
   SCIP*                 scip,               /**< original SCIP data structure                                  */
   SCIP_VAR**            fixedvars,          /**< array to store the variables that should be fixed in the subproblem */
   SCIP_Real*            fixedvals,          /**< array to store the fixing values to fix variables in the subproblem */
   int*                  nfixedvars,         /**< pointer to store the number of variables that should be fixed */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed         */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator                                       */
   SCIP_Bool*            success             /**< used to store whether the creation of the subproblem worked   */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                    */
   SCIP_SOL* sol;                            /* pool of solutions                          */

   int nvars;
   int nbinvars;
   int nintvars;
   int ndiscretevars;
   int i;

   assert(fixedvars != NULL);
   assert(fixedvals != NULL);

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   sol = SCIPgetBestSol(scip);
   assert(sol != NULL);

   /* compute the number of variables that should be fixed in the subproblem */
   *nfixedvars = (int)(minfixingrate * (nbinvars + nintvars));

   /* avoid the two corner cases that no or all discrete variables should be fixed */
   if( *nfixedvars == 0 || *nfixedvars == nbinvars + nintvars )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   assert(*nfixedvars < nbinvars + nintvars);

   ndiscretevars = nbinvars + nintvars;
   /* copy the binary and integer variables into fixedvars */
   BMScopyMemoryArray(fixedvars, vars, ndiscretevars);

   /* shuffle the array randomly */
   SCIPrandomPermuteArray(randnumgen, (void **)fixedvars, 0, nbinvars + nintvars);

   *success = TRUE;
   /* store the fixing values for the subset of variables that should be fixed */
   for( i = 0; i < *nfixedvars; ++i )
   {
      /* fix all randomly marked variables */
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;

      solval = SCIPgetSolVal(scip, sol, fixedvars[i]);
      lb = SCIPvarGetLbGlobal(fixedvars[i]);
      ub = SCIPvarGetUbGlobal(fixedvars[i]);
      assert(SCIPisLE(scip, lb, ub));

      /* due to dual reductions, it may happen that the solution value is not in
            the variable's domain anymore */
      if( SCIPisLT(scip, solval, lb) )
         solval = lb;
      else if( SCIPisGT(scip, solval, ub) )
         solval = ub;

      /* we cannot fix to infinite solution values, better break in this case */
      if( SCIPisInfinity(scip, REALABS(solval)) )
      {
         *success = FALSE;
         break;
      }

      /* store the possibly adjusted solution value as fixing value */
      fixedvals[i] = solval;
   }

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< mutation heuristic structure                        */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
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

/** setup and solve mutation sub-SCIP */
static
SCIP_RETCODE setupAndSolveSubscipMutation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_HEUR*            heur,               /**< mutation heuristic */
   SCIP_VAR**            fixedvars,          /**< array to store the variables that should be fixed in the subproblem */
   SCIP_Real*            fixedvals,          /**< array to store the fixing values to fix variables in the subproblem */
   int                   nfixedvars,         /**< the number of variables that should be fixed */
   SCIP_Longint          nsubnodes,          /**< node limit for the subproblem */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_VAR** subvars;                       /* subproblem's variables                              */
   SCIP_VAR** vars;                          /* original problem's variables                        */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_HEURDATA* heurdata;
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem                 */
   SCIP_Real upperbound;
   int nvars;                                /* number of original problem's variables              */
   int i;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);
   assert(fixedvars != NULL);
   assert(fixedvals != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

   /* create a problem copy as sub SCIP */
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "mutation", fixedvars, fixedvals, nfixedvars,
         heurdata->uselprows, heurdata->copycuts, &success, NULL) );

   for( i = 0; i < nvars; i++ )
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
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nsubnodes) );
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

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handlers; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no decutions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
   }

   /* add an objective cutoff */
   assert( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) );

   upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
   if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
   {
      cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip)
                     + heurdata->minimprove * SCIPgetLowerbound(scip);
   }
   else
   {
      if( SCIPgetUpperbound(scip) >= 0 )
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
      else
         cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
   }
   cutoff = MIN(upperbound, cutoff);
   SCIP_CALL(SCIPsetObjlimit(subscip, cutoff));

   /* solve the subproblem
    *
    * Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught but only in debug mode, SCIP will stop.
    */
   SCIPdebugMsg(scip, "Solve Mutation subMIP\n");
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* transfer variable statistics from sub-SCIP */
   SCIP_CALL( SCIPmergeVariableStatistics(subscip, scip, subvars, vars, nvars) );

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   heurdata->usednodes += SCIPgetNNodes(subscip);

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
         *result = SCIP_FOUNDSOL;
   }

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyMutation)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurMutation(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMutation)
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

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitMutation)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic */
static
SCIP_DECL_HEUREXIT(heurExitMutation)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMutation)
{  /*lint --e{715}*/
   SCIP_Longint maxnnodes;
   SCIP_Longint nsubnodes;                   /* node limit for the subproblem                       */

   SCIP_HEURDATA* heurdata;                  /* heuristic's data                                    */
   SCIP* subscip;                            /* the subproblem created by mutation                  */
   SCIP_VAR** fixedvars;                     /* array to store variables that should be fixed in the subproblem */
   SCIP_Real* fixedvals;                     /* array to store fixing values for the variables */

   SCIP_Real maxnnodesr;

   int nfixedvars;
   int nbinvars;
   int nintvars;

   SCIP_Bool success;

   SCIP_RETCODE retcode;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DELAYED;

   /* only call heuristic, if feasible solution is available */
   if( SCIPgetNSols(scip) <= 0 )
      return SCIP_OKAY;

   /* only call heuristic, if the best solution comes from transformed problem */
   assert(SCIPgetBestSol(scip) != NULL);
   if( SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip,SCIPgetBestSol(scip))  < heurdata->nwaitingnodes)
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( SCIPgetVarsData(scip, NULL, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* only call heuristic, if discrete variables are present */
   if( nbinvars + nintvars == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodesr = heurdata->nodesquot * SCIPgetNNodes(scip);

   /* reward mutation if it succeeded often, count the setup costs for the sub-MIP as 100 nodes */
   maxnnodesr *= 1.0 + 2.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0);
   maxnnodes = (SCIP_Longint) maxnnodesr - 100 * SCIPheurGetNCalls(heur);
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nsubnodes = maxnnodes - heurdata->usednodes;
   nsubnodes = MIN(nsubnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nsubnodes < heurdata->minnodes )
       return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if( !success )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, nbinvars + nintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, nbinvars + nintvars) );

   /* determine variables that should be fixed in the mutation subproblem */
   SCIP_CALL( determineVariableFixings(scip, fixedvars, fixedvals, &nfixedvars, heurdata->minfixingrate, heurdata->randnumgen, &success) );

   /* terminate if it is not possible to create the subproblem */
   if( !success )
   {
      SCIPdebugMsg(scip, "Could not create the subproblem -> skip call\n");
      goto TERMINATE;
   }

   *result = SCIP_DIDNOTFIND;

   /* initializing the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* setup and solve the subproblem and catch the return code */
   retcode = setupAndSolveSubscipMutation(scip, subscip, heur, fixedvars, fixedvals, nfixedvars, nsubnodes, result);

   /* free the subscip in any case */
   SCIP_CALL( SCIPfree(&subscip) );
   SCIP_CALL( retcode );

   /* free storage for subproblem fixings */
 TERMINATE:
   SCIPfreeBufferArray(scip, &fixedvals);
   SCIPfreeBufferArray(scip, &fixedvars);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the mutation primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMutation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Mutation primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMutation, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyMutation) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMutation) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitMutation) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitMutation) );

   /* add mutation primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, SCIPsumepsilon(scip), 1.0-SCIPsumepsilon(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which " HEUR_NAME " should at least improve the incumbent",
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
