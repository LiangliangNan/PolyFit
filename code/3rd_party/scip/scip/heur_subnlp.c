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

/**@file    heur_subnlp.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief   NLP local search primal heuristic using sub-SCIPs
 * @author  Stefan Vigerske
 * 
 * @todo reconstruct sub-SCIP if problem has changed
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/nlpi_ipopt.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_setppc.h"
#include "scip/heur_subnlp.h"
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
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_nlpi.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_presol.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"
#include <string.h>

#define HEUR_NAME        "subnlp"
#define HEUR_DESC        "primal heuristic that performs a local search in an NLP after fixing integer variables and presolving"
#define HEUR_DISPCHAR    SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY    -2000010
#define HEUR_FREQ        1
#define HEUR_FREQOFS     0
#define HEUR_MAXDEPTH    -1
#define HEUR_TIMING      SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP FALSE               /**< does the heuristic use a secondary SCIP instance? we set this to FALSE because we want this heuristic to also run within other heuristics */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP*                 subscip;            /**< copy of CIP where presolving and NLP solving is done */
   SCIP_Bool             triedsetupsubscip;  /**< whether we have tried to setup a sub-SCIP */
   SCIP_Bool             subscipisvalid;     /**< whether all constraints have been copied */
   SCIP_Bool             continuous;         /**< whether problem was continuous when sub-SCIP was created */
   int                   nseriousnlpierror;  /**< number of consecutive serious NLP solver failures (memout, ...) */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for global bound change events */

   int                   nvars;              /**< number of active transformed variables in SCIP */
   int                   nsubvars;           /**< number of original variables in sub-SCIP */
   SCIP_VAR**            var_subscip2scip;   /**< mapping variables in sub-SCIP to SCIP variables */
   SCIP_VAR**            var_scip2subscip;   /**< mapping variables in SCIP to sub-SCIP variables */

   SCIP_SOL*             startcand;          /**< candidate for start point for heuristic */
   SCIP_Real             startcandviol;      /**< violation of start point candidate w.r.t. constraint that reported this candidate */
   SCIP_SOL*             lastsol;            /**< pointer to last found solution (or NULL if none), not captured, thus may be dangling */

   int                   nlpverblevel;       /**< verbosity level of NLP solver */
   SCIP_Real             opttol;             /**< optimality tolerance to use for NLP solves */
   SCIP_Real             feastolfactor;      /**< factor on SCIP feasibility tolerance for NLP solves if resolving when NLP solution not feasible in CIP */
   SCIP_Real             feastol;            /**< feasibility tolerance for NLP solves */
   SCIP_Bool             tighterfeastolfailed;/**< whether we tried to use a tighter feasibility tolerance but the NLP solution was still not accepted */
   int                   maxpresolverounds;  /**< limit on number of presolve rounds in sub-SCIP */
   int                   presolveemphasis;   /**< presolve emphasis in sub-SCIP */
   SCIP_Bool             setcutoff;          /**< whether to set cutoff in sub-SCIP to current primal bound */
   SCIP_Bool             forbidfixings;      /**< whether to add constraints that forbid specific fixations that turned out to be infeasible */
   SCIP_Bool             keepcopy;           /**< whether to keep SCIP copy or to create new copy each time heuristic is applied */
   SCIP_Real             expectinfeas;       /**< when to tell NLP solver that an infeasible NLP is not unexpected */

   SCIP_Longint          iterused;           /**< number of iterations used so far (+ number of heuristic runs + number of presolve runs in subscip) */
   SCIP_Longint          iterusedokay;       /**< number of iterations used so far when NLP stopped with status okay */
   SCIP_Longint          iterusediterlim;    /**< maximal number of iterations used when NLP stopped due to iteration limit */
   int                   nnlpsolves;         /**< number of NLP solves */
   int                   nnlpsolvesokay;     /**< number of NLP solves with status okay */
   int                   nnlpsolvesiterlim;  /**< number of NLP solves that hit an iteration limit */
   int                   nnlpsolvesinfeas;   /**< number of NLP solves with status okay and infeasible */
   int                   nodesoffset;        /**< number of nodes added to the actual number of nodes when computing itercontingent */
   SCIP_Real             nodesfactor;        /**< factor to apply to number of nodes in SCIP to compute initial itercontingent */
   SCIP_Real             successrateexp;     /**< exponent for power of success rate to be multiplied with itercontingent */
   int                   iterinit;           /**< number of iterations used for initial NLP solves */
   int                   ninitsolves;        /**< number of successful NLP solves until switching to iterlimit guess and using success rate */
   int                   itermin;            /**< minimal number of iterations for NLP solves */
};


/*
 * Local methods
 */

/** indicates whether the heuristic should be running, i.e., whether we expect something nonlinear after fixing all discrete variables */
static
SCIP_RETCODE runHeuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            runheur             /**< buffer to store whether to run heuristic */
   )
{
   assert(scip != NULL);
   assert(runheur != NULL);

   /* do not run heuristic if no NLP solver is available */
   if( SCIPgetNNlpis(scip) <= 0 )
   {
      *runheur = FALSE;
      return SCIP_OKAY;
   }

   /* do not run heuristic if no NLP */
   if( !SCIPisNLPConstructed(scip) )
   {
      *runheur = FALSE;
      return SCIP_OKAY;
   }

   /* do not run heuristic if no continuous nonlinear variables in NLP */
   SCIP_CALL( SCIPhasNLPContinuousNonlinearity(scip, runheur) );

   return SCIP_OKAY;
}

/** free sub-SCIP data structure */
static
SCIP_RETCODE freeSubSCIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_VAR** subvars;
   int        nsubvars;
   int        i;
   SCIP_VAR*  var;
   SCIP_VAR*  subvar;

   assert(scip != NULL);
   assert(heurdata != NULL);

   assert(heurdata->subscip != NULL);

   SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, NULL, NULL, NULL, NULL) );
   assert(nsubvars == heurdata->nsubvars);

   /* drop global bound change events
    * release variables in SCIP and sub-SCIP
    */
   for( i = 0; i < heurdata->nsubvars; ++i )
   {
      subvar = subvars[i];
      assert(subvar != NULL);
      assert(SCIPvarGetProbindex(subvar) == i);

      var = heurdata->var_subscip2scip[SCIPvarGetProbindex(subvar)];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) <= heurdata->nvars);
      assert(!SCIPvarIsActive(var) || heurdata->var_scip2subscip[SCIPvarGetProbindex(var)] == subvar);

      SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GBDCHANGED, heurdata->eventhdlr, (SCIP_EVENTDATA*)heurdata, -1) );

      SCIP_CALL( SCIPreleaseVar(heurdata->subscip, &subvar) );
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* free variable mappings subscip -> scip and scip -> subscip */
   SCIPfreeBlockMemoryArray(scip, &heurdata->var_subscip2scip, heurdata->nsubvars);
   SCIPfreeBlockMemoryArray(scip, &heurdata->var_scip2subscip, heurdata->nvars);
   heurdata->nsubvars = 0;
   heurdata->nvars = 0;

   /* free sub-SCIP */
   SCIP_CALL( SCIPfree(&heurdata->subscip) );

   return SCIP_OKAY;
}

/** creates copy of CIP from problem in SCIP */
static
SCIP_RETCODE createSubSCIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   int nvars;
   SCIP_VAR** vars;
   SCIP_VAR** subvars;
   SCIP_VAR*  var;
   SCIP_VAR*  subvar;
   SCIP_Bool success;
   char probname[SCIP_MAXSTRLEN];
   int i;
   SCIP_HASHMAP* varsmap;
   SCIP_HASHMAP* conssmap;

   assert(heurdata != NULL);
   assert(heurdata->subscip == NULL);

   heurdata->triedsetupsubscip = TRUE;

   /* initializing the subproblem */
   SCIP_CALL( SCIPcreate(&heurdata->subscip) );

   /* create sub-SCIP copy of CIP */

   /* copy interesting plugins */
   success = TRUE;
   SCIP_CALL( SCIPcopyPlugins(scip, heurdata->subscip,
         FALSE, /* readers */
         FALSE, /* pricers */
         TRUE,  /* conshdlrs */
         FALSE, /* conflicthdlrs */
         TRUE,  /* presolvers */
         FALSE, /* relaxators */
         FALSE, /* separators */
         FALSE, /* cutselectors */
         TRUE,  /* propagators */
         FALSE, /* heuristics */
         TRUE,  /* eventhandler */
         TRUE,  /* nodeselectors (SCIP gives an error if there is none) */
         FALSE, /* branchrules */
         TRUE,  /* displays */
         FALSE, /* tables */
         FALSE, /* dialogs */
         TRUE,  /* expression handlers */
         TRUE,  /* nlpis */
         TRUE,  /* message handler */
         &success) );
   if( !success )
   {
      SCIPdebugMsg(scip, "failed to copy some plugins to sub-SCIP, continue anyway\n");
   }

   /* check if we still have NLPI's in subscip */
   if( SCIPgetNNlpis(heurdata->subscip) <= 0 )
   {
      SCIPdebugMsg(scip, "none of the NLPIs from main SCIP copied into sub-SCIP, give up heuristic.\n");
      SCIP_CALL( SCIPfree(&heurdata->subscip) );

      return SCIP_OKAY;
   }

   /* copy parameter settings */
   SCIP_CALL( SCIPcopyParamSettings(scip, heurdata->subscip) );

   /* create problem in sub-SCIP */
   /* get name of the original problem and add "subnlp" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_subnlp", SCIPgetProbName(scip));
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPhashmapCreate(&varsmap, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&conssmap, SCIPblkmem(scip), SCIPgetNConss(scip)) );
   SCIP_CALL( SCIPcopyProb(scip, heurdata->subscip, varsmap, conssmap, TRUE, probname) );

   /* copy all variables */
   SCIP_CALL( SCIPcopyVars(scip, heurdata->subscip, varsmap, conssmap, NULL, NULL, 0, TRUE) );

   /* copy as many constraints as possible */
   SCIP_CALL( SCIPcopyConss(scip, heurdata->subscip, varsmap, conssmap, TRUE, FALSE, &heurdata->subscipisvalid) );
   SCIPhashmapFree(&conssmap);
   if( !heurdata->subscipisvalid )
   {
      SCIPdebugMsg(scip, "failed to copy some constraints to sub-SCIP, continue anyway\n");
   }

   /* create arrays translating scip transformed vars to subscip original vars, and vice versa
    * capture variables in SCIP and sub-SCIP
    * catch global bound change events
    */

   SCIP_CALL( SCIPgetVarsData(heurdata->subscip, &subvars, &heurdata->nsubvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &heurdata->var_subscip2scip, heurdata->nsubvars) );

   heurdata->nvars = nvars;
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &heurdata->var_scip2subscip, heurdata->nvars) );

   /* we need to get all subscip variables, also those which are copies of fixed variables from the main scip
    * therefore we iterate over the hashmap
    */
   for( i = 0; i < SCIPhashmapGetNEntries(varsmap); ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      entry = SCIPhashmapGetEntry(varsmap, i);
      if( entry != NULL )
      {
         var    = (SCIP_VAR*) SCIPhashmapEntryGetOrigin(entry);
         subvar = (SCIP_VAR*) SCIPhashmapEntryGetImage(entry);
         assert(subvar != NULL);
         assert(SCIPvarGetProbindex(subvar) >= 0);
         assert(SCIPvarGetProbindex(subvar) <= heurdata->nsubvars);

         if( SCIPvarIsActive(var) )
         {
            assert(SCIPvarGetProbindex(var) <= heurdata->nvars);
            assert(heurdata->var_scip2subscip[SCIPvarGetProbindex(var)] == NULL);  /* assert that we have no mapping for this var yet */
            heurdata->var_scip2subscip[SCIPvarGetProbindex(var)] = subvar;
         }

         assert(heurdata->var_subscip2scip[SCIPvarGetProbindex(subvar)] == NULL);  /* assert that we have no mapping for this subvar yet */
         heurdata->var_subscip2scip[SCIPvarGetProbindex(subvar)] = var;
      }
   }

   for( i = 0; i < heurdata->nsubvars; ++i )
   {
      subvar = SCIPgetVars(heurdata->subscip)[i];
      assert(SCIPvarGetProbindex(subvar) == i);
      var = heurdata->var_subscip2scip[i];

      SCIP_CALL( SCIPcaptureVar(scip, var) );
      SCIP_CALL( SCIPcaptureVar(heurdata->subscip, subvar) );

      assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetLbGlobal(subvar)));
      assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), SCIPvarGetUbGlobal(subvar)));

      SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GBDCHANGED, heurdata->eventhdlr, (SCIP_EVENTDATA*)heurdata, NULL) );
   }

#ifndef NDEBUG
   for( i = 0; i < heurdata->nvars; ++i )
   {
      assert(heurdata->var_scip2subscip[i] == NULL || (SCIP_VAR*)SCIPhashmapGetImage(varsmap, (void*)vars[i]) == heurdata->var_scip2subscip[i]);
   }
   for( i = 0; i < heurdata->nsubvars; ++i )
   {
      assert(heurdata->var_subscip2scip[i] != NULL);
      assert((SCIP_VAR*)SCIPhashmapGetImage(varsmap, (void*)heurdata->var_subscip2scip[i]) == subvars[i]);
   }
#endif

   /* do not need hashmap anymore */
   SCIPhashmapFree(&varsmap);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(heurdata->subscip, "misc/catchctrlc", FALSE) );

   /* disable keeping solutions from one subscip solve for next solve (with usually different fixings) */
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "limits/maxorigsol", 0) );

#ifdef SCIP_DEBUG
   /* for debugging, enable SCIP output */
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "display/verblevel", 5) );
#else
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "display/verblevel", 0) );
#endif

   /* reset some limits to default values, in case users changed them in main scip (SCIPcopy copies parameter values :-() */
   SCIP_CALL( SCIPresetParam(heurdata->subscip, "limits/absgap") );
   SCIP_CALL( SCIPresetParam(heurdata->subscip, "limits/bestsol") );
   SCIP_CALL( SCIPresetParam(heurdata->subscip, "limits/gap") );
   SCIP_CALL( SCIPresetParam(heurdata->subscip, "limits/restarts") );
   SCIP_CALL( SCIPresetParam(heurdata->subscip, "limits/solutions") );
   SCIP_CALL( SCIPresetParam(heurdata->subscip, "limits/time") );
   SCIP_CALL( SCIPresetParam(heurdata->subscip, "limits/totalnodes") );

   /* we remember here which way (continuous or not) we went, in case all binary and integer vars get fixed in root */
   heurdata->continuous = SCIPgetNBinVars(heurdata->subscip) == 0 && SCIPgetNIntVars(heurdata->subscip) == 0;
   if( !heurdata->continuous )
   {
      /* set presolve maxrounds and emphasis; always disable components presolver
       * heuristics and separators were not copied into subscip, so should not need to switch off
       */
      if( !SCIPisParamFixed(heurdata->subscip, "presolving/maxrounds") )
      {
         SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "presolving/maxrounds", heurdata->maxpresolverounds) );
      }
      SCIP_CALL( SCIPsetPresolving(heurdata->subscip, (SCIP_PARAMSETTING)heurdata->presolveemphasis, TRUE) );
      if( !SCIPisParamFixed(heurdata->subscip, "constraints/components/maxprerounds") )
      {
         SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "constraints/components/maxprerounds", 0) );
      }
   }
   else
   {
      /* for continuous problems, disable presolve and move subscip into a stage where it has a NLP
       * the only reason why we don't solve the NLP in the main SCIP is that we want global variable bounds for the NLP
       */
      SCIP_RETCODE retcode;

      SCIP_CALL( SCIPtransformProb(heurdata->subscip) );

      SCIP_CALL( SCIPsetPresolving(heurdata->subscip, SCIP_PARAMSETTING_OFF, TRUE) );
      SCIP_CALL( SCIPpresolve(heurdata->subscip) );

      if( SCIPgetStage(heurdata->subscip) != SCIP_STAGE_PRESOLVED || SCIPgetNVars(heurdata->subscip) == 0 )
      {
         /* presolve found problem infeasible, solved it, or stopped due to some limit
          * all a bit strange, since problem should be the same as original, presolve was disabled, and we didn't set any limits
          * we will give up and not run the heuristic
          */
         SCIP_CALL( freeSubSCIP(scip, heurdata) );
         return SCIP_OKAY;
      }

      /* do initial solve, i.e., "solve" root node with node limit 0 (should do scip.c::initSolve and then stop immediately in solve.c::SCIPsolveCIP) */
      SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 0LL) );
      retcode = SCIPsolve(heurdata->subscip);

      /* errors in solving the subproblem should not kill the overall solving process
       * hence, the return code is caught and a warning is printed
       */
      if( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while initializing subproblem in subnlp heuristic; sub-SCIP terminated with code <%d>\n", retcode);
         SCIP_CALL( freeSubSCIP(scip, heurdata) );
         return SCIP_OKAY;
      }

      /* If we are in stage "solved" (strange) or have no NLP (also strange), then do not run heuristic, too */
      if( SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED || !SCIPisNLPConstructed(heurdata->subscip) )
      {
         SCIP_CALL( freeSubSCIP(scip, heurdata) );
         return SCIP_OKAY;
      }

      assert(SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVING);
      assert(SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_NODELIMIT);
      assert(SCIPisNLPConstructed(heurdata->subscip));
   }

   return SCIP_OKAY;
}

/** process variable global bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   int            idx;

   assert(scip      != NULL);
   assert(event     != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata  != NULL);

   var = SCIPeventGetVar(event);
   assert(var != NULL);

   idx = SCIPvarGetProbindex(var);
   /* if event corresponds to an active variable, we can easily look up the corresponding subvar
    * if it is an inactive variable that has been copied to the subproblem,
    * then we need to check the subscip2scip mapping
    * @todo we could do this faster if we keep the variables mapping from SCIPcopy around
    */
   if( idx >= 0 )
   {
      assert(idx < heurdata->nvars);

      subvar = heurdata->var_scip2subscip[idx];
   }
   else
   {
      for( idx = 0; idx < heurdata->nsubvars; ++idx )
      {
         if( heurdata->var_subscip2scip[idx] == var )
            break;
      }
      assert(idx < heurdata->nsubvars);
      subvar = SCIPgetVars(heurdata->subscip)[idx];
   }
   assert(subvar != NULL);

   if( SCIPeventGetType(event) & SCIP_EVENTTYPE_GLBCHANGED )
   {
      SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, SCIPvarGetLbGlobal(var)) );
   }

   if( SCIPeventGetType(event) & SCIP_EVENTTYPE_GUBCHANGED )
   {
      SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, SCIPvarGetUbGlobal(var)) );
   }

   return SCIP_OKAY;
}

/* creates a SCIP_SOL in our SCIP space out of the solution from NLP solver in sub-SCIP */
static
SCIP_RETCODE createSolFromNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_SOL**            sol,                /**< buffer to store solution value; if pointing to NULL, then a new solution is created, otherwise values in the given one are overwritten */
   SCIP_HEUR*            authorheur          /**< the heuristic which should be registered as author of the solution */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR**     vars;
   int            nvars;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   SCIP_Real      solval;
   int            i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol  != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(SCIPhasNLPSolution(heurdata->subscip));

   if( *sol == NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, sol, authorheur) );
   }
   else
   {
      SCIPsolSetHeur(*sol, authorheur);
   }

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   assert(nvars >= heurdata->nvars);
   for( i = 0; i < heurdata->nvars; ++i )
   {
      var = vars[i];
      assert(var != NULL);
      assert(SCIPvarIsActive(var));  /* SCIPgetVarsData should have given us only active vars */

      subvar = heurdata->var_scip2subscip[i];
      if( subvar == NULL )
         solval = MIN(MAX(0.0, SCIPvarGetLbLocal(var)), SCIPvarGetUbLocal(var));  /*lint !e666*/
      else
         solval = SCIPvarGetNLPSol(subvar);

      assert(solval != SCIP_INVALID);  /*lint !e777*/
      SCIP_CALL( SCIPsetSolVal(scip, *sol, var, solval) );
   }

   for( ; i < nvars; ++i )
   {
      var = vars[i];
      assert(var != NULL);
      assert(SCIPvarIsActive(var));  /* SCIPgetVarsData should have given us only active vars */

      solval = MIN(MAX(0.0, SCIPvarGetLbLocal(var)), SCIPvarGetUbLocal(var));  /*lint !e666*/
      SCIP_CALL( SCIPsetSolVal(scip, *sol, var, solval) );
   }

   return SCIP_OKAY;
}

/** creates SCIP solution from NLP and tries adding to SCIP or only checks feasibility */
static
SCIP_RETCODE processNLPSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_HEUR*            authorheur,         /**< the heuristic that should be the author of solution, if any */
   SCIP_RESULT*          result,             /**< buffer to store result FOUNDSOL if a solution has been found and accepted */
   SCIP_SOL*             resultsol           /**< a solution where to store found solution values, if any, or NULL if to try adding to SCIP */
)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assert(SCIPhasNLPSolution(heurdata->subscip));

   if( resultsol == NULL )
   {
      /* resultsol NULL means we should try adding the sol to SCIP */
      if( SCIPisLE(scip, SCIPgetNLPObjval(heurdata->subscip), SCIPgetUpperbound(scip)) )
      {
         /* solution is feasible and should improve upper bound, so try adding it to SCIP */
         SCIP_SOL*  sol;
         SCIP_Bool  stored;

         sol = NULL;
         SCIP_CALL( createSolFromNLP(scip, heur, &sol, authorheur) );

         heurdata->lastsol = sol; /* remember just the pointer so we might recognize if this solution comes back as startingpoint */
#ifdef SCIP_DEBUG
         /* print the infeasibilities to stdout */
         SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, FALSE, TRUE, &stored) );
#else
         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, FALSE, TRUE, &stored) );
#endif

         if( stored )
         {
            /* SCIP stored solution (yippi!), so we are done */
            if( heurdata->nlpverblevel >= 1 )
            {
               SCIPinfoMessage(scip, NULL, "SCIP stored solution from NLP solve\n");
            }
            else
            {
               SCIPdebugMsg(scip, "SCIP stored solution from NLP solve\n");
            }

            *result = SCIP_FOUNDSOL;
         }
         else
         {
            if( heurdata->nlpverblevel >= 1 )
            {
               SCIPinfoMessage(scip, NULL, "solution reported by NLP solver not stored by SCIP\n");
            }
            else
            {
               SCIPdebugMsg(scip, "solution reported by NLP solver not stored by SCIP\n");
            }
         }
      }
      else if( heurdata->nlpverblevel >= 1 )
      {
         SCIPinfoMessage(scip, NULL, "subnlp solution objval %e is above the primal bound %e\n",
            SCIPgetNLPObjval(heurdata->subscip), SCIPgetUpperbound(scip));
      }
   }
   else
   {
      /* only create a solution and pass it back in resultsol, do not add to SCIP */
      SCIP_Bool feasible;

      SCIP_CALL( createSolFromNLP(scip, heur, &resultsol, authorheur) );

      heurdata->lastsol = resultsol;
#ifdef SCIP_DEBUG
      /* print the infeasibilities to stdout */
      SCIP_CALL( SCIPcheckSol(scip, resultsol, TRUE, TRUE, TRUE, FALSE, TRUE, &feasible) );
#else
      SCIP_CALL( SCIPcheckSol(scip, resultsol, FALSE, FALSE, TRUE, FALSE, TRUE, &feasible) );
#endif
      if( feasible )
      {
         /* SCIP find solution feasible, so we are done */
         if( heurdata->nlpverblevel >= 1 )
         {
            SCIPinfoMessage(scip, NULL, "solution reported by NLP solver feasible for SCIP\n");
         }
         else
         {
            SCIPdebugMsg(scip, "solution reported by NLP solver feasible for SCIP\n");
         }
         *result = SCIP_FOUNDSOL;
      }
      else
      {
         if( heurdata->nlpverblevel >= 1 )
         {
            SCIPinfoMessage(scip, NULL, "solution reported by NLP solver not feasible for SCIP\n");
         }
         else
         {
            SCIPdebugMsg(scip, "solution reported by NLP solver not feasible for SCIP\n");
         }
      }
   }

   return SCIP_OKAY;
}

/* creates a SCIP_SOL in our SCIP space out of the SCIP_SOL from a sub-SCIP */
static
SCIP_RETCODE createSolFromSubScipSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_SOL**            sol,                /**< buffer to store solution value; if pointing to NULL, then a new solution is created, otherwise values in the given one are overwritten */
   SCIP_SOL*             subsol,             /**< solution of sub-SCIP */
   SCIP_HEUR*            authorheur          /**< the heuristic which should be registered as author of the solution */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR**     vars;
   int            nvars;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   SCIP_Real      solval;
   int            i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol  != NULL);
   assert(subsol != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( *sol == NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, sol, authorheur) );
   }
   else
   {
      SCIPsolSetHeur(*sol, authorheur);
   }

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   assert(nvars >= heurdata->nvars);
   for( i = 0; i < heurdata->nvars; ++i )
   {
      var = vars[i];
      assert(var != NULL);
      assert(SCIPvarIsActive(var));

      subvar = heurdata->var_scip2subscip[i];
      if( subvar == NULL )
         solval = MIN(MAX(0.0, SCIPvarGetLbLocal(var)), SCIPvarGetUbLocal(var));  /*lint !e666*/
      else
         solval = SCIPgetSolVal(heurdata->subscip, subsol, subvar);

      assert(solval != SCIP_INVALID);  /*lint !e777*/
      SCIP_CALL( SCIPsetSolVal(scip, *sol, var, solval) );
   }

   for( ; i < nvars; ++i )
   {
      var = vars[i];
      assert(var != NULL);
      assert(SCIPvarIsActive(var));

      solval = MIN(MAX(0.0, SCIPvarGetLbLocal(var)), SCIPvarGetUbLocal(var));  /*lint !e666*/
      SCIP_CALL( SCIPsetSolVal(scip, *sol, var, solval) );
   }

   return SCIP_OKAY;
}

/** finds an iteration limit */  /*lint --e{715}*/
static
int calcIterLimit(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   /* if we hit more often an iterlimit than we were successful (termstatus=okay), then allow for more iterations:
    * take twice the maximal iterusage on solves that hit the iterlimit
    */
   if( heurdata->nnlpsolvesiterlim > heurdata->nnlpsolvesokay )
      return MAX(heurdata->itermin, 2 * heurdata->iterusediterlim);  /*lint !e712*/

   /* if we had sufficiently many successful solves, then take twice the average of previous iterusages on successful solves */
   if( heurdata->nnlpsolvesokay >= heurdata->ninitsolves )
      return MAX(heurdata->itermin, 2 * heurdata->iterusedokay / heurdata->nnlpsolvesokay);  /*lint !e712*/

   /* if we had too few successful solves, then still ensure that we allow for at least iterinit iterations */
   if( heurdata->nnlpsolvesokay > 0 )
      return MAX3(heurdata->itermin, heurdata->iterinit, 2 * heurdata->iterusedokay / heurdata->nnlpsolvesokay);  /*lint !e712*/

   /* if we had no successful solve so far and none that hit an iterlimit, e.g., we are at the first NLP solve, then use iterinit */
   return MAX(heurdata->itermin, heurdata->iterinit);
}

/** solves the subNLP specified in subscip */
static
SCIP_RETCODE solveSubNLP(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< buffer to store result, DIDNOTFIND, FOUNDSOL, or CUTOFF        */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_SOL*             resultsol           /**< a solution where to store found solution values, if any, or NULL if to try adding to SCIP */
   )
{
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   SCIP_RETCODE   retcode;
   SCIP_Real*     startpoint;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   int            i;
   SCIP_HEUR*     authorheur;   /* the heuristic which will be the author of a solution, if found */
   SCIP_Real      timelimit;
   SCIP_Bool      expectinfeas;
   SCIP_NLPSTATISTICS nlpstatistics;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);
   assert(SCIPisTransformed(heurdata->subscip));

   /* get remaining SCIP solve time; if no time left, then stop */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit <= 0.0 )
         return SCIP_OKAY;
   }
   /* set timelimit for NLP solve and in case presolve is unexpectedly expensive */
   SCIP_CALL( SCIPsetRealParam(heurdata->subscip, "limits/time", timelimit) );

   /* if the refpoint comes from a heuristic, then make it the author of a found solution,
    * otherwise let the subNLP heuristic claim authorship
    * TODO: I doubt that this has much effect; for the statistics, the number of solutions found by a heuristic
    *       seems to be computed as the increase in number of solutions before and after a heuristic is run
    *       check this and maybe change
    */
   if( refpoint == NULL || SCIPsolGetHeur(refpoint) == NULL )
      authorheur = heur;
   else
      authorheur = SCIPsolGetHeur(refpoint);

   if( !heurdata->continuous )
   {
      /* presolve sub-SCIP
       *  set node limit to 1 so that presolve can go
       */
      SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 1LL) );
      SCIP_CALL( SCIPpresolve(heurdata->subscip) );

      /* count one presolve round as on NLP iteration for now
       * plus one extra for all the setup cost
       * this is mainly to avoid that the primal heuristics runs all the time on instances that are solved in the subscip-presolve
       */
      heurdata->iterused += 1 + SCIPgetNPresolRounds(scip);  /*lint !e776*/

      if( SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED )
      {
         /* presolve probably found the subproblem infeasible */
         SCIPdebugMsg(scip, "SCIP returned from presolve in stage solved with status %d and %d sols\n", SCIPgetStatus(heurdata->subscip), SCIPgetNSols(heurdata->subscip));
         /* if presolve found subproblem infeasible, report this to caller by setting *result to cutoff */
         if( SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_INFEASIBLE )
            *result = SCIP_CUTOFF;
      }
      else if( SCIPgetStage(heurdata->subscip) == SCIP_STAGE_PRESOLVING )
      {
         /* presolve was stopped because some still existing limit was hit (e.g., memory) */
         SCIPdebugMsg(scip, "SCIP returned from presolve in stage presolving with status %d and %d sols\n", SCIPgetStatus(heurdata->subscip), SCIPgetNSols(heurdata->subscip));
         /* if presolve found subproblem infeasible, report this to caller by setting *result to cutoff */
         if( SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_INFEASIBLE )
            *result = SCIP_CUTOFF;
      }
      else
      {
         assert(SCIPgetStage(heurdata->subscip) == SCIP_STAGE_PRESOLVED);

         if( SCIPgetNVars(heurdata->subscip) > 0 )
         {
            /* do initial solve, i.e., "solve" root node with node limit 0 (should do scip.c::initSolve and then stop immediately in solve.c::SCIPsolveCIP) */
            SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 0LL) );
            retcode = SCIPsolve(heurdata->subscip);

            /* If no NLP was constructed, then there were no nonlinearities after presolve.
             * So we increase the nodelimit to 1 and hope that SCIP will find some solution to this probably linear subproblem.
             */
            if( retcode == SCIP_OKAY && SCIPgetStage(heurdata->subscip) != SCIP_STAGE_SOLVED && !SCIPisNLPConstructed(heurdata->subscip) )
            {
               SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 1LL) );
               retcode = SCIPsolve(heurdata->subscip);
            }
         }
         else
         {
            /* If all variables were removed by presolve, but presolve did not end with status SOLVED,
             * then we run solve, still with nodelimit=1, and hope to find some (maybe trivial) solution.
             */
            retcode = SCIPsolve(heurdata->subscip);
         }

         /* errors in solving the subproblem should not kill the overall solving process
          * hence, the return code is caught and a warning is printed
          */
         if( retcode != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while solving subproblem in subnlp heuristic; sub-SCIP terminated with code <%d>\n", retcode);
            return SCIP_OKAY;
         }
      }

      /* we should either have variables, or the problem was trivial, in which case it should have been presolved or solved */
      assert(SCIPgetNVars(heurdata->subscip) > 0 || SCIPgetStage(heurdata->subscip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED);

      SCIPdebug( SCIP_CALL( SCIPprintStatistics(heurdata->subscip, NULL) ); )

      /* if sub-SCIP found solutions already, then pass them to main scip */
      for( i = 0; i < SCIPgetNSols(heurdata->subscip); ++i )
      {
         if( resultsol == NULL )
         {
            SCIP_Bool stored;
            SCIP_SOL* sol;

            sol = NULL;
            SCIP_CALL( createSolFromSubScipSol(scip, heur, &sol, SCIPgetSols(heurdata->subscip)[i], authorheur) );

            heurdata->lastsol = sol; /* remember just the pointer so we might recognize if this solution comes back as startingpoint */
            SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, FALSE, TRUE, &stored) );
            if( stored )
            {
               if( heurdata->nlpverblevel >= 1 )
               {
                  SCIPinfoMessage(scip, NULL, "SCIP stored solution from sub-SCIP root node\n");
               }
               else
               {
                  SCIPdebugMsg(scip, "SCIP stored solution from sub-SCIP root node\n");
               }
               *result = SCIP_FOUNDSOL;
               break;
            }
            else
            {
               if( heurdata->nlpverblevel >= 1 )
               {
                  SCIPinfoMessage(scip, NULL, "SCIP did not store sub-SCIP optimal solution\n");
               }
               else
               {
                  SCIPdebugMsg(scip, "SCIP did not store sub-SCIP optimal solution\n");
               }
            }
         }
         else
         {
            SCIP_Bool feasible;

            SCIP_CALL( createSolFromSubScipSol(scip, heur, &resultsol, SCIPgetSols(heurdata->subscip)[i], authorheur) );

            heurdata->lastsol = resultsol;
            SCIP_CALL( SCIPcheckSol(scip, resultsol, FALSE, FALSE, TRUE, FALSE, TRUE, &feasible) );
            if( feasible )
            {
               if( heurdata->nlpverblevel >= 1 )
               {
                  SCIPinfoMessage(scip, NULL, "SCIP solution from sub-SCIP root node is feasible\n");
               }
               else
               {
                  SCIPdebugMsg(scip, "SCIP solution from sub-SCIP root node is feasible\n");
               }
               *result = SCIP_FOUNDSOL;
               break;
            }
            else
            {
               if( heurdata->nlpverblevel >= 1 )
               {
                  SCIPinfoMessage(scip, NULL, "SCIP solution from sub-SCIP root node is not feasible\n");
               }
               else
               {
                  SCIPdebugMsg(scip, "SCIP solution from sub-SCIP root node is not feasible\n");
               }
            }
         }
      }

      /* if subscip is infeasible here, we signal this to the caller */
      if( SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_INFEASIBLE )
      {
         if( heurdata->nlpverblevel >= 1 )
         {
            SCIPinfoMessage(scip, NULL, "sub-SCIP detected infeasibility\n");
         }
         else
         {
            SCIPdebugMsg(scip, "sub-SCIP detected infeasibility\n");
         }

         assert(SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED);
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* if we stopped for some other reason, or there is no NLP, we also stop */
      if( SCIPgetStage(heurdata->subscip) <= SCIP_STAGE_PRESOLVED || SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED || !SCIPisNLPConstructed(heurdata->subscip) )
         return SCIP_OKAY;

      /* in most cases, the status should be nodelimit
       * in some cases, if the sub-SCIP is very easy, it may report optimal, so we do not need invoke an NLP solver
       * if the presolve found the problem infeasible, then there is no use in solving an NLP
       * if the user interrupted or a timelimit was reached, then we should also stop here
       * unbounded is very unlikely to happen, in most cases, it should have been concluded in the main scip already
       */
      switch( SCIPgetStatus(heurdata->subscip) )
      {
         case SCIP_STATUS_NODELIMIT:
            break; /* this is the status that is most likely happening */
         case SCIP_STATUS_TOTALNODELIMIT:
         case SCIP_STATUS_STALLNODELIMIT:
         case SCIP_STATUS_GAPLIMIT:
         case SCIP_STATUS_SOLLIMIT:
         case SCIP_STATUS_BESTSOLLIMIT:
            /* these should not happen, but if one does, it's safe to return */
            SCIPABORT();    /*lint -fallthrough*/
         case SCIP_STATUS_OPTIMAL:
         case SCIP_STATUS_INFEASIBLE:
         case SCIP_STATUS_USERINTERRUPT:
         case SCIP_STATUS_TIMELIMIT:
         case SCIP_STATUS_MEMLIMIT:
         case SCIP_STATUS_UNBOUNDED:
         case SCIP_STATUS_INFORUNBD:
            return SCIP_OKAY;
         default:
            SCIPerrorMessage("unexpected status of sub-SCIP: <%d>\n", SCIPgetStatus(heurdata->subscip));
            return SCIP_ERROR;
      } /*lint !e788*/
   }
   else
   {
      /* for continuous problem, createSubSCIP() should have put us into a state where we can invoke the NLP solver */
      assert(SCIPisNLPConstructed(heurdata->subscip));
      assert(SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVING);
      assert(SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_NODELIMIT);
   }

   /* set starting values (=refpoint, if not NULL; otherwise LP solution (or pseudo solution)) */
   SCIP_CALL( SCIPallocBufferArray(scip, &startpoint, SCIPgetNNLPVars(heurdata->subscip)) );

   if( heurdata->nlpverblevel >= 3 )
   {
      SCIPinfoMessage(scip, NULL, "set NLP starting point\n");
   }

   for( i = 0; i < SCIPgetNNLPVars(heurdata->subscip); ++i )
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      subvar = SCIPgetNLPVars(heurdata->subscip)[i];

      /* gets corresponding original variable */
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&subvar, &scalar, &constant) );
      if( subvar == NULL )
      {
         startpoint[i] = constant;

         if( heurdata->nlpverblevel >= 3 && !SCIPisZero(heurdata->subscip, startpoint[i]) )
         {
            SCIPinfoMessage(scip, NULL, "%s = %e\n", SCIPvarGetName(SCIPgetNLPVars(heurdata->subscip)[i]), startpoint[i]);
         }

         continue;
      }

      assert(SCIPvarGetProbindex(subvar) >= 0);
      assert(SCIPvarGetProbindex(subvar) <  heurdata->nsubvars);
      var = heurdata->var_subscip2scip[SCIPvarGetProbindex(subvar)];
      if( var == NULL || REALABS(SCIPgetSolVal(scip, refpoint, var)) > 1.0e+12 )
         startpoint[i] = MIN(MAX(0.0, SCIPvarGetLbGlobal(subvar)), SCIPvarGetUbGlobal(subvar));  /*lint !e666*/
      else
         /* scalar*subvar+constant corresponds to nlpvar[i], so nlpvar[i] gets value scalar*varval+constant */
         startpoint[i] = scalar * SCIPgetSolVal(scip, refpoint, var) + constant;

      if( heurdata->nlpverblevel >= 3 && !SCIPisZero(heurdata->subscip, startpoint[i]) )
      {
         SCIPinfoMessage(scip, NULL, "%s = %e\n", SCIPvarGetName(SCIPgetNLPVars(heurdata->subscip)[i]), startpoint[i]);
      }
   }
   SCIP_CALL( SCIPsetNLPInitialGuess(heurdata->subscip, startpoint) );

   SCIPfreeBufferArray(scip, &startpoint);

   *result = SCIP_DIDNOTFIND;

   /* if we had many (fraction > expectinfeas) infeasible NLPs, then tell NLP solver to expect an infeasible problem */
   expectinfeas = FALSE;
   if( heurdata->expectinfeas == 0.0 )  /* to keep original behavior on default settings */
      expectinfeas = TRUE;
   else if( heurdata->nnlpsolvesokay > heurdata->ninitsolves && heurdata->nnlpsolvesinfeas > heurdata->expectinfeas * heurdata->nnlpsolvesokay )
      expectinfeas = TRUE;

   /* let the NLP solver do its magic */
   SCIPdebugMsg(scip, "start NLP solve with iteration limit %d\n", calcIterLimit(scip, heurdata));
   SCIP_CALL( SCIPsolveNLP(heurdata->subscip,
      .iterlimit = calcIterLimit(scip, heurdata),
      .opttol = heurdata->opttol,
      .feastol = heurdata->feastol,
      .verblevel = (unsigned short)heurdata->nlpverblevel,
      .expectinfeas = expectinfeas
   ) );  /*lint !e666*/

   SCIPdebugMsg(scip, "NLP solver returned with termination status %d and solution status %d, objective value is %g\n",
      SCIPgetNLPTermstat(heurdata->subscip), SCIPgetNLPSolstat(heurdata->subscip), SCIPgetNLPObjval(heurdata->subscip));

   /* add NLP solve statistics from subscip to main SCIP, so they show up in final statistics
    * for continuous problems, we also ask to reset statistics, since we do not retransform subSCIP in the next run (which would reset all stats)
    * (merging statistics once in exitsol is too late, since they may be printed before)
    */
   SCIPmergeNLPIStatistics(heurdata->subscip, scip, heurdata->continuous);

   if( SCIPgetNLPTermstat(heurdata->subscip) >= SCIP_NLPTERMSTAT_OUTOFMEMORY )
   {
      /* oops, something did not go well at all */
     if( heurdata->nlpverblevel >= 1 )
     {
        SCIPinfoMessage(scip, NULL, "NLP solver in subNLP heuristic for problem <%s> returned with bad termination status %d.\n",
           SCIPgetProbName(scip), SCIPgetNLPTermstat(heurdata->subscip));
     }

     ++(heurdata->nseriousnlpierror);
     SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL,
         "NLP solver in subNLP heuristic for problem <%s> returned with bad termination status %d. This was the %d%s successive time.\n",
         SCIPgetProbName(scip), SCIPgetNLPTermstat(heurdata->subscip), heurdata->nseriousnlpierror,
         heurdata->nseriousnlpierror == 1 ? "st" : heurdata->nseriousnlpierror == 2 ? "nd" : heurdata->nseriousnlpierror == 3 ? "rd" : "th");
      if( heurdata->nseriousnlpierror >= 5 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Will not run subNLP heuristic again for this run.\n");
         SCIP_CALL( freeSubSCIP(scip, heurdata) );
      }
      return SCIP_OKAY;
   }
   heurdata->nseriousnlpierror = 0;

   SCIP_CALL( SCIPgetNLPStatistics(heurdata->subscip, &nlpstatistics) );

   SCIPdebugMsg(scip, "NLP solver used %d iterations and %g seconds; violation cons %g, bounds %g\n",
      nlpstatistics.niterations, nlpstatistics.totaltime, nlpstatistics.consviol, nlpstatistics.boundviol);

   heurdata->iterused += nlpstatistics.niterations;
   ++heurdata->nnlpsolves;
   if( SCIPgetNLPTermstat(heurdata->subscip) == SCIP_NLPTERMSTAT_OKAY )
   {
      ++heurdata->nnlpsolvesokay;
      heurdata->iterusedokay += nlpstatistics.niterations;

      if( (SCIPgetNLPSolstat(heurdata->subscip) == SCIP_NLPSOLSTAT_GLOBINFEASIBLE) || (SCIPgetNLPSolstat(heurdata->subscip) == SCIP_NLPSOLSTAT_LOCINFEASIBLE) )
         ++heurdata->nnlpsolvesinfeas;
   }
   else if( SCIPgetNLPTermstat(heurdata->subscip) == SCIP_NLPTERMSTAT_ITERLIMIT )
   {
      ++heurdata->nnlpsolvesiterlim;
      heurdata->iterusediterlim = MAX(heurdata->iterusediterlim, nlpstatistics.niterations);
   }

   if( SCIPgetNLPSolstat(heurdata->subscip) > SCIP_NLPSOLSTAT_FEASIBLE )
      return SCIP_OKAY;

   /* create SCIP solution, check whether feasible, and try adding to SCIP (if resultsol==NULL) */
   SCIP_CALL( processNLPSol(scip, heur, authorheur, result, resultsol) );

   if( *result == SCIP_FOUNDSOL || !SCIPisLE(scip, SCIPgetNLPObjval(heurdata->subscip), SCIPgetUpperbound(scip)) )
      return SCIP_OKAY;

   /* if solution was not added to SCIP, then either
    * - the objective function value was not good enough,
    * - the NLP was missing some constraints of the original CIP, or
    * - the solution is feasible for the presolved CIP, but slightly infeasible for the unpresolved problem
    *
    * The first case we can check easily (see if() above).
    * For the last case, we try whether tightening the feasibility tolerance for the NLP solve may help.
    * If that doesn't help, we guess that we are in the second case and will not try a tighter feastol anymore.
    */

   /* if we tried with a tighter feastol before, but solution was still not accepted, then don't try again */
   if( heurdata->tighterfeastolfailed )
      return SCIP_OKAY;

   /* if resolve with tighter feastol is disabled, then don't do anything */
   if( heurdata->feastolfactor == 1.0 )
      return SCIP_OKAY;

   /* if we have already used a tighter feastol, then give up */
   if( heurdata->feastol < SCIPfeastol(scip) )
      return SCIP_OKAY;

   /* if original CIP is continuous, then we have not done any presolve, so it shouldn't have caused problems */
   if( heurdata->continuous )
      return SCIP_OKAY;

   /* if solution is NLP-feasible for a tightened tolerance already, then there is no use in resolving with that tighter feastol */
   if( MAX(nlpstatistics.consviol, nlpstatistics.boundviol) <= heurdata->feastolfactor * heurdata->feastol )
      return SCIP_OKAY;

   /* let the NLP solver redo its magic
    * as iterlimit, we use the number of iterations it took for the first solve, or itermin
    */
   SCIPdebugMsg(scip, "start NLP solve with iteration limit %d\n", calcIterLimit(scip, heurdata));
   SCIP_CALL( SCIPsolveNLP(heurdata->subscip,
      .iterlimit = MAX(heurdata->itermin, nlpstatistics.niterations),
      .opttol = heurdata->opttol,
      .feastol = heurdata->feastolfactor * heurdata->feastol,
      .verblevel = (unsigned short)heurdata->nlpverblevel,
      .warmstart = TRUE
   ) );  /*lint !e666*/

   SCIPdebugMsg(scip, "NLP solver returned with termination status %d and solution status %d, objective value is %g\n",
      SCIPgetNLPTermstat(heurdata->subscip), SCIPgetNLPSolstat(heurdata->subscip), SCIPgetNLPObjval(heurdata->subscip));

   /* add NLP solve statistics from subscip to main SCIP again, so they show up in final statistics */
   SCIPmergeNLPIStatistics(heurdata->subscip, scip, heurdata->continuous);

   /* some serious problem: just pretend it didn't happen */
   if( SCIPgetNLPTermstat(heurdata->subscip) >= SCIP_NLPTERMSTAT_OUTOFMEMORY )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetNLPStatistics(heurdata->subscip, &nlpstatistics) );
   SCIPdebugMsg(scip, "NLP solver used %d iterations and %g seconds; violation cons %g, bounds %g\n",
      nlpstatistics.niterations, nlpstatistics.totaltime, nlpstatistics.consviol, nlpstatistics.boundviol);

   /* we account only the extra iterations for this unusual NLP solve, but don't add anything else to our statistics (nnlpsolved, etc) */
   heurdata->iterused += nlpstatistics.niterations;

   /* if failed to get a feasible NLP solution now, then nothing to do */
   if( SCIPgetNLPSolstat(heurdata->subscip) > SCIP_NLPSOLSTAT_FEASIBLE )
      return SCIP_OKAY;

   SCIP_CALL( processNLPSol(scip, heur, authorheur, result, resultsol) );

   /* if successful, then use tighter feastol for all NLP solves from now on
    * if still not accepted, then don't try this again
    * (maybe the NLP is incomplete; we could give up on running this heur completely, but for now let the successrate factor in heurExec take care of running it less often)
    */
   if( *result == SCIP_FOUNDSOL )
      heurdata->feastol *= heurdata->feastolfactor;
   else
      heurdata->tighterfeastolfailed = TRUE;

   return SCIP_OKAY;
}


/** adds a set covering or bound disjunction constraint to the original problem */
static
SCIP_RETCODE forbidFixation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   SCIP_VAR** subvars;
   int nsubvars;
   int nsubbinvars;
   int nsubintvars;
   SCIP_VAR* var;
   SCIP_VAR* subvar;
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   int nconsvars;
   char name[SCIP_MAXSTRLEN];
   int i;
   SCIP_Real fixval;

   assert(scip != NULL);

   SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );
   assert(nsubvars == heurdata->nsubvars);

   if( nsubbinvars == 0 && nsubintvars == 0 )
   {
      /* If we did not fix any discrete variables but found the "sub"CIP infeasible, then also the CIP is infeasible. */
      SCIPdebugMsg(scip, "heur_subnlp found subCIP infeasible after fixing no variables, something is strange here...\n");
      return SCIP_OKAY;
   }

   /* initialize */
   cons = NULL;
   consvars = NULL;

   /* create constraint name */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "subnlp_cutoff");

   /* if all discrete variables in the CIP are binary, then we create a set covering constraint
    *  sum_{x_i fixed at 0} x_i + sum_{x_i fixed at 1} ~x_i >= 1
    */
   if( nsubintvars == 0 )
   {
      /* allocate memory for constraint variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nsubbinvars) );

      /* get fixations of discrete variables 
       * to be sure, we take the values that were put into the subCIP before
       */
      nconsvars = 0;
      for( i = nsubbinvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL || SCIPisEQ(scip, SCIPvarGetLbGlobal(subvar), SCIPvarGetUbGlobal(subvar))); /* otherwise we should have exited in the variable fixation loop */
         if( var == NULL )
            continue;

         fixval = SCIPvarGetLbGlobal(subvar);
         assert(fixval == SCIPvarGetUbGlobal(subvar)); /* variable should be fixed in sub-SCIP */  /*lint !e777*/
         assert(fixval == 0.0 || fixval == 1.0); /* we have rounded values before fixing */

         if( fixval == 0.0 )
         {
            /* variable fixed at lower bound */
            consvars[nconsvars] = var;
         }
         else
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, var, &consvars[nconsvars]) );
         }

         ++nconsvars;
      }

      /* create conflict constraint
       * In undercover, ConsLogicor is used, since then the inequality is not added to the LP.
       * However, I may want to use Setcover to avoid that the same fixing is computed by some LP based heuristic again.
       */
      SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, name, nconsvars, consvars,
            FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   }
   else
   {
      /* if there are also integer variable, then create a bound disjunction constraint
       * x_1 >= fixval_1 + 1 || x_1 <= fixval_1 - 1 || x_2 >= fixval_2 + 1 || x_2 <= fixval_2 - 1 || ...
       */
      SCIP_BOUNDTYPE* boundtypes;
      SCIP_Real* bounds;

      /* allocate memory for constraint variables, boundtypes, and bounds
       * (there should be at most two literals for each integer variable)
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars,   nsubbinvars + 2*nsubintvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, nsubbinvars + 2*nsubintvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &bounds,     nsubbinvars + 2*nsubintvars) );

      /* get fixations of discrete variables 
       * to be sure, we take the values that were put into the subCIP before
       */
      nconsvars = 0;
      for( i = nsubbinvars + nsubintvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL || SCIPisEQ(scip, SCIPvarGetLbGlobal(subvar), SCIPvarGetUbGlobal(subvar))); /* otherwise we should have exited in the variable fixation loop */

         if( var == NULL )
            continue;

         fixval = SCIPvarGetLbGlobal(subvar);
         assert(fixval == SCIPvarGetUbGlobal(subvar)); /* variable should be fixed in sub-SCIP */   /*lint !e777*/
         assert(SCIPceil(scip, fixval - 0.5) == fixval); /* we have rounded values before fixing */ /*lint !e777*/
         assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY || SCIPvarGetLbGlobal(var) == fixval || SCIPvarGetUbGlobal(var) == fixval); /* for binaries, the fixval should be either 0.0 or 1.0 */  /*lint !e777*/

         if( SCIPvarGetLbGlobal(var) < fixval )
         {
            assert(nconsvars < nsubbinvars + 2*nsubintvars);

            /* literal x_i <= fixval-1 */
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_UPPER;
            bounds[nconsvars]     = fixval - 1.0;
            consvars[nconsvars]   = var;
            ++nconsvars;
         }

         if( SCIPvarGetUbGlobal(var) > fixval )
         {
            assert(nconsvars < nsubbinvars + 2*nsubintvars);

            /* literal x_i >= fixval+1 */
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_LOWER;
            bounds[nconsvars]     = fixval + 1.0;
            consvars[nconsvars]   = var;
            ++nconsvars;
         }
      }

      /* create conflict constraint */
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, nconsvars, consvars, boundtypes, bounds,
            FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

      SCIPfreeBufferArray(scip, &bounds);
      SCIPfreeBufferArray(scip, &boundtypes);
      SCIPfreeBufferArray(scip, &consvars);
   }

   /* add and release constraint if created successfully */
   if( cons != NULL )
   {
      SCIPdebugMsg(scip, "adding constraint to forbid fixation in main problem\n");
      /* SCIPdebugPrintCons(scip, cons, NULL); */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &consvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySubNlp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSubNlp(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSubNlp)
{
   SCIP_HEURDATA* heurdata;
   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->subscip == NULL);
   assert(heurdata->var_subscip2scip == NULL);
   assert(heurdata->var_scip2subscip == NULL);
   assert(heurdata->startcand == NULL);

   SCIPfreeBlockMemory(scip, &heurdata);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSubNlp)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->subscip == NULL);

   /* reset or initialize some flags and counters */
   heurdata->feastol = SCIPfeastol(scip);
   heurdata->tighterfeastolfailed = FALSE;
   heurdata->triedsetupsubscip = FALSE;
   heurdata->nseriousnlpierror = 0;
   heurdata->iterused = 0;
   heurdata->iterusedokay = 0;
   heurdata->iterusediterlim = 0;
   heurdata->nnlpsolves = 0;
   heurdata->nnlpsolvesokay = 0;
   heurdata->nnlpsolvesiterlim = 0;
   heurdata->nnlpsolvesinfeas = 0;

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolSubNlp)
{
   assert(scip != NULL);
   assert(heur != NULL);

   /* if the heuristic is called at the root node, we want to be called directly after the initial root LP solve */
   if( SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP | HEUR_TIMING);

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolSubNlp)
{
   SCIP_HEURDATA* heurdata;
   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic's data */  
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( heurdata->subscip != NULL )
   {
      SCIP_CALL( freeSubSCIP(scip, heurdata) );
      heurdata->triedsetupsubscip = FALSE;
   }

   /* free start candidate */
   if( heurdata->startcand != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
   }

   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSubNlp)
{  /*lint --e{666,715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool      runheur;
   SCIP_Real      itercontingent;

   assert(scip != NULL);
   assert(heur != NULL);

   /* obviously, we did not do anything yet */
   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* if triedsetupsubscip and keepcopy and subscip == NULL, then we tried to setup a subSCIP before, but failed due to some serious error
    * thus, we should do not need to try again
    *
    * otherwise, we continue and let SCIPapplyHeurSubNlp try to create subscip
    */
   if( heurdata->subscip == NULL && heurdata->keepcopy && heurdata->triedsetupsubscip )
      return SCIP_OKAY;

   /* before we run the heuristic for the first time, check whether we want to run the heuristic at all */
   if( SCIPheurGetNCalls(heur) == 0 )
   {
      SCIP_CALL( runHeuristic(scip, &runheur) );
      if( !runheur )
         return SCIP_OKAY;
   }

   if( heurdata->startcand == NULL )
   {
      /* if no start candidate is given, we consider the LP solution of the current node */

      /* however, if the node was already detected to be infeasible, then there is no point to look at its LP solution */
      if( nodeinfeasible )
         return SCIP_OKAY;

      /* at least if we are not called the first time, we call the heuristic only if an optimal LP solution is available 
       * if we are called the first time and the LP is unbounded, then we are quite desperate and still give the NLP a try
       */
      if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      {
         if( SCIPgetNNodes(scip) > 1 || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
         {
            *result = SCIP_DELAYED;
            SCIPdebugMsg(scip, "NLP heuristic delayed because no start candidate given and no LP solution available; LP status = %d\n", SCIPgetLPSolstat(scip));
            return SCIP_OKAY;
         }
         else
         {
            SCIPdebugMsg(scip, "LP is unbounded in root node, so we are quite desperate; run NLP heuristic and pray\n");
         }
      }
      else if( SCIPgetNLPBranchCands(scip) > 0 )
      {
         /* only call heuristic, if there are no fractional variables */
         *result = SCIP_DELAYED;
         SCIPdebugMsg(scip, "NLP heuristic delayed because no start candidate given and current LP solution is fractional\n");
         return SCIP_OKAY;
      }
      else if( !SCIPisInfinity(scip, SCIPgetPrimalbound(scip)) && SCIPisEQ(scip, SCIPgetLocalDualbound(scip), SCIPgetPrimalbound(scip)) )
      {
         /* only call heuristic, if there is still room for improvement in the current node */
         SCIPdebugMsg(scip, "NLP heuristic delayed because lower and upper bound coincide in current node\n");
         return SCIP_OKAY;
      }
      SCIPdebugMsg(scip, "using current LP solution as startcand\n");
   }
   else
   {
      SCIPdebugMsg(scip, "have startcand from heur %s\n", SCIPsolGetHeur(heurdata->startcand) ? SCIPheurGetName(SCIPsolGetHeur(heurdata->startcand)) : "NULL");
   }

   /* check if enough nodes have been processed so that we want to run the heuristic again */

   /* compute the contingent on number of iterations that the NLP solver is allowed to use
    * we make it depending on the current number of processed nodes
    */
   itercontingent = heurdata->nodesfactor * (SCIPgetNNodes(scip) + heurdata->nodesoffset);
   /* weight by previous success of heuristic if we have been running already
    * require at least ninitsolves many runs that didn't run into the NLP iterlimit
    * (so if we are still in the phase of finding a good iterlimit, do not consider success rate so far)
    */
   if( heurdata->successrateexp > 0.0 && SCIPheurGetNCalls(heur) - heurdata->nnlpsolvesiterlim >= heurdata->ninitsolves )
      itercontingent *= pow((SCIPheurGetNSolsFound(heur) + 1.0) / (SCIPheurGetNCalls(heur) + 1.0), heurdata->successrateexp);
   /* subtract the number of iterations used for all NLP solves so far */
   itercontingent -= heurdata->iterused;

   /* check whether the itercontingent is sufficient for the iteration limit we would use */
   if( itercontingent < calcIterLimit(scip, heurdata) )
   {
      /* not enough iterations left to start NLP solver */
      SCIPdebugMsg(scip, "skip NLP heuristic; contingent=%f; iterlimit=%d; success ratio=%g\n",
         itercontingent, calcIterLimit(scip, heurdata), pow((SCIPheurGetNSolsFound(heur) + 1.0) / (SCIPheurGetNCalls(heur) + 1.0), heurdata->successrateexp));
      return SCIP_OKAY;
   }

   /* so far we have not found any solution, but now we are willing to search for one */
   *result = SCIP_DIDNOTFIND;

   if( heurdata->nlpverblevel >= 1 )
   {
      SCIPinfoMessage(scip, NULL, "calling subnlp heuristic\n");
   }

   SCIP_CALL( SCIPapplyHeurSubNlp(scip, heur, result, heurdata->startcand, NULL) );

   /* SCIP does not like cutoff as return, so we say didnotfind, since we did not find a solution */
   if( *result == SCIP_CUTOFF )
      *result = SCIP_DIDNOTFIND;

   /* forget startcand */
   if( heurdata->startcand != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
   }

   /* reset timing, if it was changed temporary (at the root node) */
   if( heurtiming != HEUR_TIMING )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the NLP local search primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSubNlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Nlp primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* include variable event handler */
   heurdata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &heurdata->eventhdlr, HEUR_NAME, "propagates a global bound change to the sub-SCIP",
         processVarEvent, NULL) );
   assert(heurdata->eventhdlr != NULL);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSubNlp, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySubNlp) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSubNlp) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSubNlp) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolSubNlp) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolSubNlp) );

   /* add Nlp primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/" HEUR_NAME "/nlpverblevel",
         "verbosity level of NLP solver",
         &heurdata->nlpverblevel, FALSE, 0, 0, USHRT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/" HEUR_NAME "/nodesoffset",
         "number of nodes added to the current number of nodes when computing itercontingent (higher value runs heuristic more often in early search)",
         &heurdata->nodesoffset, FALSE, 1600, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesfactor",
         "factor on number of nodes in SCIP (plus nodesoffset) to compute itercontingent (higher value runs heuristics more frequently)",
         &heurdata->nodesfactor, FALSE, 0.3, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/successrateexp",
         "exponent for power of success rate to be multiplied with itercontingent (lower value decreases impact of success rate)",
         &heurdata->successrateexp, FALSE, 1.0, 0.0, DBL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/" HEUR_NAME "/iterinit",
         "number of iterations used for initial NLP solves",
         &heurdata->iterinit, FALSE, 300, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/" HEUR_NAME "/ninitsolves",
         "number of successful NLP solves until switching to iterlimit guess and using success rate",
         &heurdata->ninitsolves, FALSE, 2, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/" HEUR_NAME "/itermin",
         "minimal number of iterations for NLP solves",
         &heurdata->itermin, FALSE, 20, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/opttol",
         "absolute optimality tolerance to use for NLP solves",
         &heurdata->opttol, TRUE, SCIPdualfeastol(scip), 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/feastolfactor",
         "factor on SCIP feasibility tolerance for NLP solves if resolving when NLP solution not feasible in CIP",
         &heurdata->feastolfactor, FALSE, 0.1, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxpresolverounds",
         "limit on number of presolve rounds in sub-SCIP (-1 for unlimited, 0 for no presolve)",
         &heurdata->maxpresolverounds, FALSE, -1, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/presolveemphasis",
         "presolve emphasis in sub-SCIP (0: default, 1: aggressive, 2: fast, 3: off)",
         &heurdata->presolveemphasis, FALSE, (int)SCIP_PARAMSETTING_FAST, (int)SCIP_PARAMSETTING_DEFAULT, (int)SCIP_PARAMSETTING_OFF, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam (scip, "heuristics/" HEUR_NAME "/setcutoff",
         "whether to set cutoff in sub-SCIP to current primal bound",
         &heurdata->setcutoff, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam (scip, "heuristics/" HEUR_NAME "/forbidfixings",
         "whether to add constraints that forbid specific fixings that turned out to be infeasible",
         &heurdata->forbidfixings, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam (scip, "heuristics/" HEUR_NAME "/keepcopy",
         "whether to keep SCIP copy or to create new copy each time heuristic is applied",
         &heurdata->keepcopy, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/expectinfeas",
         "percentage of NLP solves with infeasible status required to tell NLP solver to expect an infeasible NLP",
         &heurdata->expectinfeas, FALSE, 0.0, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}

/** main procedure of the subNLP heuristic */
SCIP_RETCODE SCIPapplyHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< pointer to store result of: did not run, solution found, no solution found, or fixing is infeasible (cutoff) */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_SOL*             resultsol           /**< a solution where to store found solution values, if any, or NULL if to try adding to SCIP */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   int            i;
   SCIP_Real      cutoff = SCIPinfinity(scip);

   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* try to setup NLP if not tried before */
   if( heurdata->subscip == NULL && !heurdata->triedsetupsubscip )
   {
      SCIP_CALL( createSubSCIP(scip, heurdata) );
   }

   *result = SCIP_DIDNOTRUN;

   /* if subSCIP could not be created, then do not run */
   if( heurdata->subscip == NULL )
      return SCIP_OKAY;

   assert(heurdata->nsubvars > 0);
   assert(heurdata->var_subscip2scip != NULL);

   /* fix discrete variables in sub-SCIP */
   if( !heurdata->continuous )
   {
      SCIP_Real  fixval;
      SCIP_VAR** subvars;
      int        nsubvars;
      int        nsubbinvars;
      int        nsubintvars;
      SCIP_Bool  infeas;
      SCIP_Bool  tightened;

      /* transform sub-SCIP, so variable fixing are easily undone by free-transform */
      assert(!SCIPisTransformed(heurdata->subscip));
      SCIP_CALL( SCIPtransformProb(heurdata->subscip) );

      SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );
      assert(nsubvars == heurdata->nsubvars);

      /* fix discrete variables to values in startpoint */
      for( i = nsubbinvars + nsubintvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL);

         /* at this point, variables in subscip and in our scip should have same bounds */
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(subvar), SCIPvarGetLbGlobal(var)));
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(subvar), SCIPvarGetUbGlobal(var)));

         fixval = SCIPgetSolVal(scip, refpoint, var);

         /* only run heuristic on integer feasible points (unless we are on an unbounded LP) */
         if( !SCIPisFeasIntegral(scip, fixval) )
         {
            if( refpoint != NULL || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIPdebugMsg(scip, "skip NLP heuristic because start candidate not integer feasible: var <%s> has value %g\n", SCIPvarGetName(var), fixval);
               goto CLEANUP;
            }
         }
         /* if we do not really have a startpoint, then we should take care that we do not fix variables to very large values
          *  thus, we set to 0.0 here and project on bounds below
          */
         if( REALABS(fixval) > 1E+10 && refpoint == NULL && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
            fixval = 0.0;

         /* fixing variables to infinity causes problems, we should not have been passed such a solution as refpoint */
         assert(!SCIPisInfinity(scip, REALABS(fixval)));

         /* round fractional variables to the nearest integer */
         fixval = SCIPround(scip, fixval);

         /* adjust value to the global bounds of the corresponding SCIP variable */
         fixval = MAX(fixval, SCIPvarGetLbGlobal(var));  /*lint !e666*/
         fixval = MIN(fixval, SCIPvarGetUbGlobal(var));  /*lint !e666*/

         /* SCIPdebugMsg(scip, "fix variable <%s> to %g\n", SCIPvarGetName(var), fixval); */
         SCIP_CALL( SCIPtightenVarLb(heurdata->subscip, subvar, fixval, TRUE, &infeas, &tightened) );
         if( !infeas )
         {
            SCIP_CALL( SCIPtightenVarUb(heurdata->subscip, subvar, fixval, TRUE, &infeas, &tightened) );
         }
         if( infeas )
         {
            SCIPdebugMsg(scip, "skip NLP heuristic because start candidate not feasible: fixing var <%s> to value %g is infeasible\n", SCIPvarGetName(var), fixval);
            goto CLEANUP;
         }
      }

      /* if there is already a solution, possibly add an objective cutoff in sub-SCIP
       * we do this here only for problems with discrete variables, since the cutoff may be useful when presolving the subscip
       * for the NLP solver, a cutoff is useless at best
       */
      if( SCIPgetNSols(scip) > 0 && heurdata->setcutoff )
      {
         cutoff = SCIPgetUpperbound(scip);
         assert( !SCIPisInfinity(scip, cutoff) );

         SCIP_CALL( SCIPsetObjlimit(heurdata->subscip, cutoff) );
         SCIPdebugMsg(scip, "set objective limit %g\n", cutoff);
      }
   }
   else
   {
      /* for continuous problems, we should already be in the transformed stage */
      assert(SCIPisTransformed(heurdata->subscip));
   }

   /* solve the subNLP and try to add solution to SCIP */
   SCIP_CALL( solveSubNLP(scip, heur, result, refpoint, resultsol) );

   if( heurdata->subscip == NULL )
   {
      /* something horrible must have happened that we decided to give up completely on this heuristic */
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   if( *result == SCIP_CUTOFF )
   {
      if( heurdata->subscipisvalid && SCIPgetNActivePricers(scip) == 0 )
      {
         /* if the subNLP is valid and turned out to be globally infeasible (i.e., proven by SCIP), then we forbid this fixation in the main problem */
         if( SCIPisInfinity(scip, cutoff) && heurdata->forbidfixings )
         {
            SCIP_CALL( forbidFixation(scip, heurdata) );
         }
      }
      else
      {
         /* if the subNLP turned out to be globally infeasible but we are not sure that we have a valid copy, we change to DIDNOTFIND */
         *result = SCIP_DIDNOTFIND;
      }
   }

 CLEANUP:
   if( !heurdata->continuous )
   {
      SCIP_CALL( SCIPfreeTransform(heurdata->subscip) );
   }

   /* if the heuristic was applied before solving has started, then destroy subSCIP, since EXITSOL may not be called
    * also if keepcopy is disabled, then destroy subSCIP
    */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING || !heurdata->keepcopy )
   {
      SCIP_CALL( freeSubSCIP(scip, heurdata) );
      heurdata->triedsetupsubscip = FALSE;
   }

   return SCIP_OKAY;
}

/** updates the starting point for the NLP heuristic
 * 
 * Is called by a constraint handler that handles nonlinear constraints when a check on feasibility of a solution fails.
 */
SCIP_RETCODE SCIPupdateStartpointHeurSubNlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< NLP heuristic */
   SCIP_SOL*             solcand,            /**< solution candidate */
   SCIP_Real             violation           /**< constraint violation of solution candidate */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(solcand != NULL);
   assert(SCIPisPositive(scip, violation));

   /* too early or the game is over already: no more interest in starting points */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( heurdata->subscip == NULL )
   {
      /* if we do not have a sub-SCIP, but tried to set one up before or will never create a subSCIP, then do not need a starting point */
      SCIP_Bool runheur;
      if( heurdata->triedsetupsubscip )
         return SCIP_OKAY;
      if( SCIPheurGetFreq(heur) < 0 )
         return SCIP_OKAY;
      SCIP_CALL( runHeuristic(scip, &runheur) );
      if( !runheur )
         return SCIP_OKAY;
   }

   /* if the solution is the one we created (last), then it is useless to use it as starting point again
    * (we cannot check SCIPsolGetHeur()==heur, as subnlp may not be registered as author of the solution)
    */
   if( heurdata->lastsol == solcand )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "consider solution candidate with violation %g and objective %g from %s\n",
      violation, SCIPgetSolTransObj(scip, solcand), SCIPsolGetHeur(solcand) ? SCIPheurGetName(SCIPsolGetHeur(solcand)) : "tree");

   /* if we have no point yet, or the new point has a lower constraint violation, or it has a better objective function value, then take the new point */
   if( heurdata->startcand == NULL || violation < heurdata->startcandviol ||
      SCIPisRelGT(scip, SCIPgetSolTransObj(scip, heurdata->startcand), SCIPgetSolTransObj(scip, solcand)) )
   {
      if( heurdata->startcand != NULL )
      {
         SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
      }
      SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->startcand, solcand) );
      SCIP_CALL( SCIPunlinkSol(scip, heurdata->startcand) );
      heurdata->startcandviol = violation;

      /* remember which heuristic proposed the candidate */
      SCIPsolSetHeur(heurdata->startcand, SCIPgetSolHeur(scip, solcand));
   }

   return SCIP_OKAY;
}

/** gets startpoint candidate to be used in next call to NLP heuristic, or NULL if none */
SCIP_SOL* SCIPgetStartCandidateHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->startcand;
}
