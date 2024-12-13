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

/**@file   sepa_rapidlearning.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  rapidlearning separator
 * @author Timo Berthold
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#ifndef NDEBUG
#include <string.h>
#endif

#include "scip/sepa_rapidlearning.h"
#include "scip/scipdefplugins.h"
#include "scip/heuristics.h"
#include "scip/pub_var.h"

#define SEPA_NAME              "rapidlearning"
#define SEPA_DESC               "rapid learning heuristic and separator"
#define SEPA_PRIORITY          -1200000
#define SEPA_FREQ                     5
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP           TRUE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_APPLYCONFLICTS     TRUE /**< should the found conflicts be applied in the original SCIP? */
#define DEFAULT_APPLYBDCHGS        TRUE /**< should the found global bound deductions be applied in the original SCIP?
                                         *   apply only if conflicts and incumbent solution will be copied too
                                         */
#define DEFAULT_APPLYINFERVALS     TRUE /**< should the inference values be used as initialization in the original SCIP? */
#define DEFAULT_REDUCEDINFER      FALSE /**< should the inference values only be used when rapid learning found other reductions? */
#define DEFAULT_APPLYPRIMALSOL     TRUE /**< should the incumbent solution be copied to the original SCIP? */
#define DEFAULT_APPLYSOLVED        TRUE /**< should a solved status be copied to the original SCIP? */

#define DEFAULT_CHECKEXEC          TRUE /**< check whether rapid learning should be executed */
#define DEFAULT_CHECKDEGANERACY    TRUE /**< should local LP degeneracy be checked? */
#define DEFAULT_CHECKDUALBOUND    FALSE /**< should the progress on the dual bound be checked? */
#define DEFAULT_CHECKLEAVES       FALSE /**< should the ratio of leaves proven to be infeasible and exceeding the
                                         *   cutoff bound be checked? */
#define DEFAULT_CHECKOBJ          FALSE /**< should the local objection function be checked? */
#define DEFAULT_CHECKNSOLS         TRUE /**< should the number of solutions found so far be checked? */
#define DEFAULT_MINDEGENERACY       0.7 /**< minimal degeneracy threshold to allow local rapid learning */
#define DEFAULT_MININFLPRATIO      10.0 /**< minimal threshold of inf/obj leaves to allow local rapid learning */
#define DEFAULT_MINVARCONSRATIO     2.0 /**< minimal ratio of unfixed variables in relation to basis size to
                                         *   allow local rapid learning */
#define DEFAULT_NWAITINGNODES      100L /**< number of nodes that should be processed before rapid learning is
                                         *   executed locally based on the progress of the dualbound */

#define DEFAULT_MAXNVARS          10000 /**< maximum problem size (variables) for which rapid learning will be called */
#define DEFAULT_MAXNCONSS         10000 /**< maximum problem size (constraints) for which rapid learning will be called */
#define DEFAULT_MAXCALLS            100 /**< maximum number of overall calls */

#define DEFAULT_MINNODES            500 /**< minimum number of nodes considered in rapid learning run */
#define DEFAULT_MAXNODES           5000 /**< maximum number of nodes considered in rapid learning run */

#define DEFAULT_CONTVARS          FALSE /**< should rapid learning be applied when there are continuous variables? */
#define DEFAULT_CONTVARSQUOT        0.3 /**< maximal portion of continuous variables to apply rapid learning */
#define DEFAULT_LPITERQUOT          0.2 /**< maximal fraction of LP iterations compared to node LP iterations */
#define DEFAULT_COPYCUTS           TRUE /**< should all active cuts from the cutpool of the
                                         *   original scip be copied to constraints of the subscip */


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_Real             lpiterquot;         /**< maximal fraction of LP iterations compared to node LP iterations */
   SCIP_Real             mindegeneracy;      /**< minimal degeneracy threshold to allow local rapid learning */
   SCIP_Real             mininflpratio;      /**< minimal threshold of inf/obj leaves to allow local rapid learning */
   SCIP_Real             minvarconsratio;    /**< minimal ratio of unfixed variables in relation to basis size to
                                              *   allow local rapid learning */
   int                   maxnvars;           /**< maximum problem size (variables) for which rapid learning will be called */
   int                   maxnconss;          /**< maximum problem size (constraints) for which rapid learning will be called */
   int                   maxcalls;           /**< maximum number of overall calls */
   int                   minnodes;           /**< minimum number of nodes considered in rapid learning run */
   int                   maxnodes;           /**< maximum number of nodes considered in rapid learning run */
   SCIP_Longint          nwaitingnodes;      /**< number of nodes that should be processed before rapid learning is executed locally
                                              *   based on the progress of the dualbound */
   SCIP_Bool             applybdchgs;        /**< should the found global bound deductions be applied in the original SCIP? */
   SCIP_Bool             applyconflicts;     /**< should the found conflicts be applied in the original SCIP? */
   SCIP_Bool             applyinfervals;     /**< should the inference values be used as initialization in the original SCIP? */
   SCIP_Bool             applyprimalsol;     /**< should the incumbent solution be copied to the original SCIP? */
   SCIP_Bool             applysolved;        /**< should a solved status ba copied to the original SCIP? */
   SCIP_Bool             checkdegeneracy;    /**< should local LP degeneracy be checked? */
   SCIP_Bool             checkdualbound;     /**< should the progress on the dual bound be checked? */
   SCIP_Bool             checkleaves;        /**< should the ratio of leaves proven to be infeasible and exceeding the
                                              *   cutoff bound be checked? */
   SCIP_Bool             checkexec;          /**< check whether rapid learning should be executed */
   SCIP_Bool             checkobj;           /**< should the (local) objective function be checked? */
   SCIP_Bool             checknsols;         /**< should number if solutions found so far be checked? */
   SCIP_Bool             contvars;           /**< should rapid learning be applied when there are continuous variables? */
   SCIP_Real             contvarsquot;       /**< maximal portion of continuous variables to apply rapid learning */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   subproblem? */
   SCIP_Bool             reducedinfer;       /**< should the inference values only be used when rapid learning found other reductions? */
};

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyRapidlearning)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaRapidlearning(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeRapidlearning)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   SCIPfreeBlockMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** setup and solve sub-SCIP */
static
SCIP_RETCODE setupAndSolveSubscipRapidlearning(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< subSCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   int                   randseed,           /**< global seed shift used in the sub-SCIP */
   SCIP_Bool             global,             /**< should rapid learning run on the global problem? */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_VAR** vars;                          /* original problem's variables */
   SCIP_VAR** subvars;                       /* subproblem's variables */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_HASHMAP* varmapbw = NULL;            /* mapping of sub-SCIP variables to SCIP variables */

   SCIP_CONSHDLR** conshdlrs = NULL;         /* array of constraint handler's that might that might obtain conflicts */
   int* oldnconss = NULL;                    /* number of constraints without rapid learning conflicts */

   SCIP_Longint nodelimit;                   /* node limit for the subproblem */

   int nconshdlrs;                           /* size of conshdlr and oldnconss array */
   int nvars;                                /* number of variables */
   int nbinvars;
   int nintvars;
   int nimplvars;
   int implstart;
   int implend;
   int restartnum;                           /* maximal number of conflicts that should be created */
   int i;                                    /* counter */

   SCIP_Bool success;                        /* was problem creation / copying constraint successful? */

   SCIP_Bool cutoff;                         /* detected infeasibility */
   int nconflicts;                           /* statistic: number of conflicts applied */
   int nbdchgs;                              /* statistic: number of bound changes applied */

   SCIP_Bool soladded = FALSE;               /* statistic: was a new incumbent found? */
   SCIP_Bool dualboundchg;                   /* statistic: was a new dual bound found? */
   SCIP_Bool disabledualreductions;          /* TRUE, if dual reductions in sub-SCIP are not valid for original SCIP,
                                              * e.g., because a constraint could not be copied or a primal solution
                                              * could not be copied back */
   int initseed;
   int seedshift;
   SCIP_Bool valid;

#ifdef SCIP_DEBUG
   int n1startinfers = 0;                    /* statistic: number of one side infer values */
   int n2startinfers = 0;                    /* statistic: number of both side infer values */
#endif

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, NULL) );

   /* initializing the subproblem */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
   valid = FALSE;

   /* copy the subproblem */
   SCIP_CALL( SCIPcopyConsCompression(scip, subscip, varmapfw, NULL, "rapid", NULL, NULL, 0, global, FALSE, FALSE, TRUE, &valid) );

   if( sepadata->copycuts )
   {
      /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
      SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, global, NULL) );
   }

   /* fill subvars array in the order of the variables of the main SCIP */
   for( i = 0; i < nvars; i++ )
   {
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);
   }
   SCIPhashmapFree(&varmapfw);

   /* change implicit integer variables to integer type */
   implstart = nbinvars + nintvars;
   implend = nbinvars + nintvars + nimplvars;
   for( i = implstart; i < implend; i++ )
   {
      SCIP_Bool infeasible;

      if( subvars[i] == NULL )
         continue;

      assert(SCIPvarGetType(subvars[i]) == SCIP_VARTYPE_IMPLINT);
      SCIP_CALL( SCIPchgVarType(subscip, subvars[i], SCIP_VARTYPE_INTEGER, &infeasible) );
      assert(!infeasible);
   }

   /* This avoids dual presolving.
    *
    * If the copy is not valid, it should be a relaxation of the problem (constraints might have failed to be copied,
    * but no variables should be missing because we stop earlier anyway if pricers are present).
    * By disabling dual presolving, conflicts and bound changes found in a relaxation are still valid for the original problem.
    */
   if( ! valid )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/allowweakdualreds", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/allowstrongdualreds", FALSE) );
   }

   SCIPdebugMsg(scip, "Copying SCIP was%s valid.\n", valid ? "" : " not");

   /* mimic an FD solver: DFS, no LP solving, 1-FUIP instead of all-FUIP, ... */
   if( SCIPisParamFixed(subscip, "lp/solvefreq") )
   {
      SCIPwarningMessage(scip, "unfixing parameter lp/solvefreq in subscip of rapidlearning\n");
      SCIP_CALL( SCIPunfixParam(subscip, "lp/solvefreq") );
   }
   if( SCIPisParamFixed(subscip, "nodeselection/dfs/stdpriority") )
   {
      SCIPwarningMessage(scip, "unfixing parameter nodeselection/dfs/stdpriority in subscip of rapidlearning\n");
      SCIP_CALL( SCIPunfixParam(subscip, "nodeselection/dfs/stdpriority") );
   }
   SCIP_CALL( SCIPsetEmphasis(subscip, SCIP_PARAMEMPHASIS_CPSOLVER, TRUE) );

   /* turn off pseudo objective propagation */
   if( !SCIPisParamFixed(subscip, "propagating/pseudoobj/freq") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "propagating/pseudoobj/freq", -1) );
   }

   /* use classic inference branching */
   if( !SCIPisParamFixed(subscip, "branching/inference/useweightedsum") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "branching/inference/useweightedsum", FALSE) );
   }

   /* only create short conflicts */
   if( !SCIPisParamFixed(subscip, "conflict/maxvarsfac") )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "conflict/maxvarsfac", 0.05) );
   }

   /* set node limit for the subproblem based on the number of LP iterations per node,
    * which are a determistic measure for the node processing time.
    *
    * Note: We scale by number of LPs + 1 because the counter is increased after solving the LP.
    */
   nodelimit = SCIPgetNLPIterations(scip) / (SCIPgetNLPs(scip) + 1);
   nodelimit = MAX(sepadata->minnodes, nodelimit);
   nodelimit = MIN(sepadata->maxnodes, nodelimit);

   /* change global random seed */
   assert(randseed >= 0);
   SCIP_CALL( SCIPgetIntParam(scip, "randomization/randomseedshift", &seedshift) );

   initseed = ((randseed + seedshift) % INT_MAX);
   SCIP_CALL( SCIPsetIntParam(subscip, "randomization/randomseedshift", initseed) );

   restartnum = 1000;

 #ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", -1) );
 #else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
 #endif

   /* set limits for the subproblem */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit/5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/restarts", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "conflict/restartnum", restartnum) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* add an objective cutoff */
   SCIP_CALL( SCIPsetObjlimit(subscip, SCIPgetUpperbound(scip)) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapbw, SCIPblkmem(scip), nvars) );

   /* store reversing mapping of variables */
   SCIP_CALL( SCIPtransformProb(subscip) );
   for( i = 0; i < nvars; ++i)
   {
      if( subvars[i] != NULL )
      {
         SCIP_CALL( SCIPhashmapInsert(varmapbw, SCIPvarGetTransVar(subvars[i]), vars[i]) );
      }
   }

   /* allocate memory for constraints storage. Each constraint that will be created from now on will be a conflict.
    * Therefore, we need to remember oldnconss to get the conflicts from the FD search.
    */
   nconshdlrs = 4;
   SCIP_CALL( SCIPallocBufferArray(scip, &conshdlrs, nconshdlrs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oldnconss, nconshdlrs) );

   /* store number of constraints before rapid learning search */
   conshdlrs[0] = SCIPfindConshdlr(subscip, "setppc");
   conshdlrs[1] = SCIPfindConshdlr(subscip, "logicor");
   conshdlrs[2] = SCIPfindConshdlr(subscip, "linear");
   conshdlrs[3] = SCIPfindConshdlr(subscip, "bounddisjunction");

   /* redundant constraints might be eliminated in presolving */
   SCIP_CALL( SCIPpresolve(subscip) );

   for( i = 0; i < nconshdlrs; ++i)
   {
      if( conshdlrs[i] != NULL )
         oldnconss[i] = SCIPconshdlrGetNConss(conshdlrs[i]);
   }

   /* solve the subproblem, abort after errors in debug mode */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* if problem was already solved do not increase limits to run again */
   if( SCIPgetStage(subscip) == SCIP_STAGE_SOLVED )
   {
      SCIPdebugMsg(scip, "Subscip was completely solved, status %d.\n", SCIPgetStatus(subscip));
   }
   /* abort solving, if limit of applied conflicts is reached */
   else if( SCIPgetNConflictConssApplied(subscip) >= restartnum )
   {
      SCIPdebugMsg(scip, "finish after %" SCIP_LONGINT_FORMAT " successful conflict calls.\n", SCIPgetNConflictConssApplied(subscip));
   }
   /* if the first 20% of the solution process were successful, proceed */
   else if( (sepadata->applyprimalsol && SCIPgetNSols(subscip) > 0 && SCIPisFeasLT(scip, SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip) ) )
      || (sepadata->applybdchgs && SCIPgetNRootboundChgs(subscip) > 0 )
      || (sepadata->applyconflicts && SCIPgetNConflictConssApplied(subscip) > 0) )
   {
      SCIPdebugMsg(scip, "proceed solving after the first 20%% of the solution process, since:\n");

      if( SCIPgetNSols(subscip) > 0 && SCIPisFeasLE(scip, SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip) ) )
      {
         SCIPdebugMsg(scip, "   - there was a better solution (%f < %f)\n",SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip));
      }
      if( SCIPgetNRootboundChgs(subscip) > 0 )
      {
         SCIPdebugMsg(scip, "   - there were %d changed variables bounds\n", SCIPgetNRootboundChgs(subscip) );
      }
      if( SCIPgetNConflictConssFound(subscip) > 0 )
      {
         SCIPdebugMsg(scip, "   - there were %" SCIP_LONGINT_FORMAT " conflict constraints created\n", SCIPgetNConflictConssApplied(subscip));
      }

      /* set node limit to 100% */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );

      /* solve the subproblem, abort after errors in debug mode */
      SCIP_CALL_ABORT( SCIPsolve(subscip) );
   }
   else
   {
      SCIPdebugMsg(scip, "do not proceed solving after the first 20%% of the solution process.\n");
   }

 #ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
 #endif

   if( SCIPallowStrongDualReds(scip) )
      disabledualreductions = FALSE;
   else
      disabledualreductions = TRUE;

   /* check, whether a solution was found */
   if( sepadata->applyprimalsol && SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until was declared to be feasible
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      soladded = FALSE;

      /* try adding solution from subSCIP to SCIP, until finding one that is accepted */
      for( i = 0; i < nsubsols && !soladded; ++i )
      {
         SCIP_SOL* newsol;

         SCIP_CALL( SCIPtranslateSubSol(scip, subscip, subsols[i], NULL, subvars, &newsol) );
         SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &soladded) );
      }
      if( !soladded || !SCIPisEQ(scip, SCIPgetSolOrigObj(subscip, subsols[i-1]), SCIPgetSolOrigObj(subscip, subsols[0])) )
         disabledualreductions = TRUE;
   }

   /* if the sub problem was solved completely, we update the dual bound */
   dualboundchg = FALSE;
   if( sepadata->applysolved && !disabledualreductions
      && (SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL || SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE) )
   {
      /* we need to multiply the dualbound with the scaling factor and add the offset,
       * because this information has been disregarded in the sub-SCIP
       */
      SCIPdebugMsg(scip, "Update old dualbound %g to new dualbound %g.\n",
         SCIPgetDualbound(scip), SCIPretransformObj(scip, SCIPgetDualbound(subscip)));

      SCIP_CALL( SCIPupdateLocalDualbound(scip, SCIPretransformObj(scip, SCIPgetDualbound(subscip))) );
      dualboundchg = TRUE;
   }

   /* check, whether conflicts were created */
   nconflicts = 0;
   if( sepadata->applyconflicts && !disabledualreductions && SCIPgetNConflictConssApplied(subscip) > 0 )
   {
      SCIP_HASHMAP* consmap;
      int hashtablesize;
      int nmaxconfs;

      assert(SCIPgetNConflictConssApplied(subscip) < (SCIP_Longint) INT_MAX);
      hashtablesize = (int) SCIPgetNConflictConssApplied(subscip);
      assert(hashtablesize < INT_MAX/5);

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), hashtablesize) );

      SCIP_CALL( SCIPgetIntParam(scip, "conflict/maxconss", &nmaxconfs) );
      if( global )
         nmaxconfs *= 20;

      /* loop over all constraint handlers that might contain conflict constraints
       * @todo select promising constraints and not greedy
       */
      for( i = 0; i < nconshdlrs && nconflicts < nmaxconfs; ++i)
      {
         /* copy constraints that have been created in FD run */
         if( conshdlrs[i] != NULL && SCIPconshdlrGetNConss(conshdlrs[i]) > oldnconss[i] )
         {
            SCIP_CONS** conss;
            int c;
            int nconss;

            nconss = SCIPconshdlrGetNConss(conshdlrs[i]);
            conss = SCIPconshdlrGetConss(conshdlrs[i]);

            /* loop over all constraints that have been added in sub-SCIP run, these are the conflicts */
            for( c = oldnconss[i]; c < nconss && nconflicts < nmaxconfs; ++c)
            {
               SCIP_CONS* cons;
               SCIP_CONS* conscopy;

               cons = conss[c];
               assert(cons != NULL);

               success = FALSE;

               /* @todo assert that flags are as they should be for conflicts */
               SCIP_CALL( SCIPgetConsCopy(subscip, scip, cons, &conscopy, conshdlrs[i], varmapbw, consmap, NULL,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), !global, FALSE, SCIPconsIsDynamic(cons),
                     SCIPconsIsRemovable(cons), FALSE, TRUE, &success) );

               if( success )
               {
                  nconflicts++;

                  SCIP_CALL( SCIPaddConflict(scip, global ? NULL : SCIPgetCurrentNode(scip), conscopy, NULL,
                     SCIP_CONFTYPE_UNKNOWN, FALSE) );
               }
               else
               {
                  SCIPdebugMsg(scip, "failed to copy conflict constraint %s back to original SCIP\n", SCIPconsGetName(cons));
               }
            }
         }
      }
      SCIPhashmapFree(&consmap);
   }

   /* check, whether tighter (global) bounds were detected */
   cutoff = FALSE;
   nbdchgs = 0;
   if( sepadata->applybdchgs && !disabledualreductions )
   {
      for( i = 0; i < nvars; ++i )
      {
         SCIP_Bool tightened;

         if( subvars[i] == NULL )
            continue;

         assert(SCIPisLE(scip, SCIPvarGetLbGlobal(vars[i]), SCIPvarGetLbGlobal(subvars[i])));
         assert(SCIPisLE(scip, SCIPvarGetLbGlobal(subvars[i]), SCIPvarGetUbGlobal(subvars[i])));
         assert(SCIPisLE(scip, SCIPvarGetUbGlobal(subvars[i]), SCIPvarGetUbGlobal(vars[i])));

         /* update the bounds of the original SCIP, if a better bound was proven in the sub-SCIP */
         if( global )
         {
#ifndef NDEBUG
            assert(SCIPgetEffectiveRootDepth(scip) == SCIPgetDepth(scip));
#else
            if( SCIPgetEffectiveRootDepth(scip) < SCIPgetDepth(scip) )
               return SCIP_INVALIDCALL;
#endif
            tightened = FALSE;

            SCIP_CALL( SCIPtightenVarUbGlobal(scip, vars[i], SCIPvarGetUbGlobal(subvars[i]), FALSE, &cutoff, &tightened) );

            if( cutoff )
               break;

            if( tightened )
               nbdchgs++;

            tightened = FALSE;

            SCIP_CALL( SCIPtightenVarLbGlobal(scip, vars[i], SCIPvarGetLbGlobal(subvars[i]), FALSE, &cutoff, &tightened) );

            if( cutoff )
               break;

            if( tightened )
               nbdchgs++;
         }
         else
         {
            tightened = FALSE;

            SCIP_CALL( SCIPtightenVarUb(scip, vars[i], SCIPvarGetUbGlobal(subvars[i]), FALSE, &cutoff, &tightened) );

            if( cutoff )
               break;

            if( tightened )
               nbdchgs++;

            tightened = FALSE;

            SCIP_CALL( SCIPtightenVarLb(scip, vars[i], SCIPvarGetLbGlobal(subvars[i]), FALSE, &cutoff, &tightened) );

            if( cutoff )
               break;

            if( tightened )
               nbdchgs++;
         }
      }
   }

   /* install start values for inference branching */
   /* @todo use different nbranching counters for pseudo cost and inference values and update inference values in the tree */
   if( sepadata->applyinfervals && global && (!sepadata->reducedinfer || soladded || nbdchgs + nconflicts > 0) )
   {
      for( i = 0; i < nvars; ++i )
      {
         SCIP_Real downinfer;
         SCIP_Real upinfer;
         SCIP_Real downvsids;
         SCIP_Real upvsids;
         SCIP_Real downconflen;
         SCIP_Real upconflen;

         if( subvars[i] == NULL )
            continue;

         /* copy downwards branching statistics */
         downvsids = SCIPgetVarVSIDS(subscip, subvars[i], SCIP_BRANCHDIR_DOWNWARDS);
         downconflen = SCIPgetVarAvgConflictlength(subscip, subvars[i], SCIP_BRANCHDIR_DOWNWARDS);
         downinfer = SCIPgetVarAvgInferences(subscip, subvars[i], SCIP_BRANCHDIR_DOWNWARDS);

         /* copy upwards branching statistics */
         upvsids = SCIPgetVarVSIDS(subscip, subvars[i], SCIP_BRANCHDIR_UPWARDS);
         upconflen = SCIPgetVarAvgConflictlength(subscip, subvars[i], SCIP_BRANCHDIR_UPWARDS);
         upinfer = SCIPgetVarAvgInferences(subscip, subvars[i], SCIP_BRANCHDIR_UPWARDS);

#ifdef SCIP_DEBUG
         /* memorize statistics */
         if( downinfer+downconflen+downvsids > 0.0 || upinfer+upconflen+upvsids != 0 )
            n1startinfers++;

         if( downinfer+downconflen+downvsids > 0.0 && upinfer+upconflen+upvsids != 0 )
            n2startinfers++;
#endif

         SCIP_CALL( SCIPinitVarBranchStats(scip, vars[i], 0.0, 0.0, downvsids, upvsids, downconflen, upconflen, downinfer, upinfer, 0.0, 0.0) );
      }
   }

#ifdef SCIP_DEBUG
   if( cutoff )
   {
      SCIPdebugMsg(scip, "Rapidlearning detected %s infeasibility.\n", global ? "global" : "local");
   }

   SCIPdebugMsg(scip, "Rapidlearning added %d %s conflicts, changed %d bounds, %s primal solution, %s dual bound improvement.\n",
      nconflicts, global ? "global" : "local", nbdchgs, soladded ? "found" : "no",  dualboundchg ? "found" : "no");

   SCIPdebugMsg(scip, "YYY Infervalues initialized on one side: %5.2f %% of variables, %5.2f %% on both sides\n",
      100.0 * n1startinfers/(SCIP_Real)nvars, 100.0 * n2startinfers/(SCIP_Real)nvars);
#endif

   /* change result pointer */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nconflicts > 0 || dualboundchg )
      *result = SCIP_CONSADDED;
   else if( nbdchgs > 0 )
      *result = SCIP_REDUCEDDOM;

   /* free local data */
   assert(oldnconss != NULL);
   assert(conshdlrs != NULL);
   assert(varmapbw != NULL);
   SCIPfreeBufferArray(scip, &oldnconss);
   SCIPfreeBufferArray(scip, &conshdlrs);
   SCIPhashmapFree(&varmapbw);

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}

/** returns whether rapid learning is allowed to run locally */
static
SCIP_RETCODE checkExec(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator's private data */
   SCIP_Bool*            run                 /**< pointer to store whether rapid learning is allowed to run */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);

   *run = FALSE;

   /* return TRUE if local exec should not be checked */
   if( !sepadata->checkexec )
   {
      *run = TRUE;
   }

   /* problem has zero objective function, i.e., it is a pure feasibility problem */
   if( !(*run) && sepadata->checkobj && SCIPgetNObjVars(scip) == 0 )
   {
         SCIPdebugMsg(scip, "-> allow local rapid learning due to global zero objective\n");

         *run = TRUE;
   }

   /* check whether a solution was found */
   if( !(*run) && sepadata->checknsols && SCIPgetNSolsFound(scip) == 0 )
   {
      SCIPdebugMsg(scip, "-> allow local rapid learning due to no solution found so far\n");

      *run = TRUE;
   }

   /* check whether the dual bound has not changed since the root node */
   if( !(*run) && sepadata->checkdualbound && sepadata->nwaitingnodes < SCIPgetNNodes(scip) )
   {
      SCIP_Real rootdualbound;
      SCIP_Real locdualbound;

      rootdualbound = SCIPgetLowerboundRoot(scip);
      locdualbound = SCIPgetLocalLowerbound(scip);

      if( SCIPisEQ(scip, rootdualbound, locdualbound) )
      {
         SCIPdebugMsg(scip, "-> allow local rapid learning due to equal dualbound\n");

         *run = TRUE;
      }
   }

   /* check leaf nodes */
   if( !(*run) && sepadata->checkleaves )
   {
      SCIP_Real ratio = (SCIPgetNInfeasibleLeaves(scip) + 1.0) / (SCIPgetNObjlimLeaves(scip) + 1.0);

      if( SCIPisLE(scip, sepadata->mininflpratio, ratio) )
      {
         SCIPdebugMsg(scip, "-> allow local rapid learning due to inf/obj leaves ratio\n");

         *run = TRUE;
      }
   }

   /* check whether all undecided integer variables have zero objective coefficient */
   if( !(*run) && sepadata->checkobj )
   {
      SCIP_Bool allzero;
      SCIP_VAR** vars;
      int ndiscvars;
      int i;

      allzero = TRUE;
      vars = SCIPgetVars(scip);
      ndiscvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);

      for( i = 0; i < ndiscvars; i++ )
      {
         assert(SCIPvarIsIntegral(vars[i]));

         /* skip locally fixed variables */
         if( SCIPisEQ(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) )
            continue;

         if( !SCIPisZero(scip, SCIPvarGetObj(vars[i])) )
         {
            allzero = FALSE;
            break;
         }
      }

      if( allzero )
      {
         SCIPdebugMsg(scip, "-> allow local rapid learning due to local zero objective\n");

         *run = TRUE;
      }
   }

   /* check degeneracy */
   if( !(*run) && sepadata->checkdegeneracy )
   {
      SCIP_Real degeneracy;
      SCIP_Real varconsratio;

      SCIP_CALL( SCIPgetLPDualDegeneracy(scip, &degeneracy, &varconsratio) );

      SCIPdebugMsg(scip, "degeneracy: %.2f ratio: %.2f\n", degeneracy, varconsratio);

      if( degeneracy >= sepadata->mindegeneracy || varconsratio >= sepadata->minvarconsratio )
      {
         SCIPdebugMsg(scip, "-> allow local rapid learning due to degeneracy\n");

         *run = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpRapidlearning)
{/*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP* subscip;
   SCIP_SEPADATA* sepadata;
   SCIP_Bool global;
   SCIP_Bool run;
   SCIP_Bool success;
   SCIP_RETCODE retcode;
   int ndiscvars;
   int i;

   assert(sepa != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   ndiscvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);

   /* only run when still not fixed binary variables exists */
   if( ndiscvars == 0 )
      return SCIP_OKAY;

   /* get separator's data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* call separator at most maxcalls times */
   if( SCIPsepaGetNCalls(sepa) >= sepadata->maxcalls )
      return SCIP_OKAY;

   /* only run for integer programs */
   if( !sepadata->contvars && ndiscvars != SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* only run if there are few enough continuous variables */
   if( sepadata->contvars && SCIPgetNContVars(scip) > sepadata->contvarsquot * SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* do not run if pricers are present */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* if the separator should be exclusive to the root node, this prevents multiple calls due to restarts */
   if( SCIPsepaGetFreq(sepa) == 0 && SCIPsepaGetNCalls(sepa) > 0 )
      return SCIP_OKAY;

   /* call separator at most once per node */
   if( SCIPsepaGetNCallsAtNode(sepa) > 0 )
      return SCIP_OKAY;

   /* the information deduced from rapid learning is globally valid only if we are at the root node; thus we can't use
    * the depth argument of the callback
    */
   global = (SCIPgetDepth(scip) <= SCIPgetEffectiveRootDepth(scip));

   /* check if rapid learning should be applied locally */
   SCIP_CALL( checkExec(scip, sepadata, &run) );

   /* @todo check whether we want to run at the root node again, e.g., inf/obj ratio is large enough */
   if( !run )
      return SCIP_OKAY;

   /* do not call rapid learning, if the problem is too big */
   if( SCIPgetNVars(scip) > sepadata->maxnvars || SCIPgetNConss(scip) > sepadata->maxnconss )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if( !success)
      return SCIP_OKAY;

   /* skip rapid learning when the sub-SCIP would contain an integer variable with an infinite bound in direction of the
    * objective function; this might lead to very bad branching decisions when enforcing a pseudo solution (#1439)
    */
   vars = SCIPgetVars(scip);
   for( i = SCIPgetNBinVars(scip); i < ndiscvars; i++ )
   {
      SCIP_Real lb = SCIPvarGetLbLocal(vars[i]);
      SCIP_Real ub = SCIPvarGetUbLocal(vars[i]);
      SCIP_Real obj = SCIPvarGetObj(vars[i]);

      if( (SCIPisNegative(scip, obj) && SCIPisInfinity(scip, ub))
         || (SCIPisPositive(scip, obj) && SCIPisInfinity(scip, -lb)) )
      {
         SCIPdebugMsg(scip, "unbounded integer variable %s (in [%g,%g]) with objective %g -> skip rapid learning\n",
            SCIPvarGetName(vars[i]), lb, ub, obj);
         return SCIP_OKAY;
      }
   }

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPcreate(&subscip) );

   retcode = setupAndSolveSubscipRapidlearning(scip, subscip, sepadata, (int)SCIPsepaGetNCalls(sepa)+1, global, result);

   SCIP_CALL( SCIPfree(&subscip) );

   return retcode;
}


/*
 * separator specific interface methods
 */

/** creates the rapidlearning separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaRapidlearning(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create rapidlearning separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpRapidlearning, NULL,
         sepadata) );

   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyRapidlearning) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeRapidlearning) );

   /* add rapidlearning separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/applyconflicts",
         "should the found conflicts be applied in the original SCIP?",
         &sepadata->applyconflicts, TRUE, DEFAULT_APPLYCONFLICTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/applybdchgs",
         "should the found global bound deductions be applied in the original SCIP?",
         &sepadata->applybdchgs, TRUE, DEFAULT_APPLYBDCHGS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/applyinfervals",
         "should the inference values be used as initialization in the original SCIP?",
         &sepadata->applyinfervals, TRUE, DEFAULT_APPLYINFERVALS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/reducedinfer",
         "should the inference values only be used when " SEPA_NAME " found other reductions?",
         &sepadata->reducedinfer, TRUE, DEFAULT_REDUCEDINFER, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/applyprimalsol",
         "should the incumbent solution be copied to the original SCIP?",
         &sepadata->applyprimalsol, TRUE, DEFAULT_APPLYPRIMALSOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/applysolved",
         "should a solved status be copied to the original SCIP?",
         &sepadata->applysolved, TRUE, DEFAULT_APPLYSOLVED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/checkdegeneracy",
         "should local LP degeneracy be checked?",
         &sepadata->checkdegeneracy, TRUE, DEFAULT_CHECKDEGANERACY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/checkdualbound",
         "should the progress on the dual bound be checked?",
         &sepadata->checkdualbound, TRUE, DEFAULT_CHECKDUALBOUND, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/checkleaves",
         "should the ratio of leaves proven to be infeasible and exceeding the cutoff bound be checked?",
         &sepadata->checkleaves, TRUE, DEFAULT_CHECKLEAVES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/checkexec",
         "check whether rapid learning should be executed",
         &sepadata->checkexec, TRUE, DEFAULT_CHECKEXEC, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/checkobj",
         "should the (local) objective function be checked?",
         &sepadata->checkobj, TRUE, DEFAULT_CHECKOBJ, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/checknsols",
         "should the number of solutions found so far be checked?",
         &sepadata->checknsols, TRUE, DEFAULT_CHECKNSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/contvars",
         "should rapid learning be applied when there are continuous variables?",
         &sepadata->contvars, TRUE, DEFAULT_CONTVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/contvarsquot",
         "maximal portion of continuous variables to apply rapid learning",
         &sepadata->contvarsquot, TRUE, DEFAULT_CONTVARSQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/lpiterquot",
         "maximal fraction of LP iterations compared to node LP iterations",
         &sepadata->lpiterquot, TRUE, DEFAULT_LPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/mindegeneracy",
         "minimal degeneracy threshold to allow local rapid learning",
         &sepadata->mindegeneracy, TRUE, DEFAULT_MINDEGENERACY, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/mininflpratio",
         "minimal threshold of inf/obj leaves to allow local rapid learning",
         &sepadata->mininflpratio, TRUE, DEFAULT_MININFLPRATIO, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/minvarconsratio",
         "minimal ratio of unfixed variables in relation to basis size to allow local rapid learning",
         &sepadata->minvarconsratio, TRUE, DEFAULT_MINVARCONSRATIO, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxnvars",
         "maximum problem size (variables) for which rapid learning will be called",
         &sepadata->maxnvars, TRUE, DEFAULT_MAXNVARS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxnconss",
         "maximum problem size (constraints) for which rapid learning will be called",
         &sepadata->maxnconss, TRUE, DEFAULT_MAXNCONSS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxcalls",
         "maximum number of overall calls",
         &sepadata->maxcalls, TRUE, DEFAULT_MAXCALLS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxnodes",
         "maximum number of nodes considered in rapid learning run",
         &sepadata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/minnodes",
         "minimum number of nodes considered in rapid learning run",
         &sepadata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "separating/" SEPA_NAME "/nwaitingnodes",
         "number of nodes that should be processed before rapid learning is executed locally based on the progress of the dualbound",
         &sepadata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &sepadata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
