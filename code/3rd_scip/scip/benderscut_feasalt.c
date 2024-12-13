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

/**@file   benderscut_feasalt.c
 * @brief  Alternative feasibility cuts for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/pub_expr.h"
#include "scip/benderscut_feasalt.h"
#include "scip/benderscut_opt.h"
#include "scip/cons_linear.h"
#include "scip/pub_benderscut.h"
#include "scip/pub_benders.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_linear.h"
#include "scip/pub_nlp.h"
#include "scip/pub_var.h"
#include "scip/scip_benders.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_nlpi.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"

#define BENDERSCUT_NAME             "feasalt"
#define BENDERSCUT_DESC             "Alternative feasibility cuts for Benders' decomposition"
#define BENDERSCUT_PRIORITY     10001
#define BENDERSCUT_LPCUT         TRUE

#define SCIP_DEFAULT_DISPLAYFREQ 20
#define SLACKVAR_NAME    "##bendersslackvar" /** the name for the Benders' slack variables added to each
                                              *  constraints in the subproblems */

struct SCIP_BenderscutData
{
   SCIP_NLPI*            nlpi;               /**< nlpi used to create the nlpi problem */
   SCIP_NLPIPROBLEM*     nlpiprob;           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_HASHMAP*         var2idx;            /**< mapping the variable to the index in the NLPI problem */
   SCIP_HASHMAP*         row2idx;            /**< mapping the rows to the index in the NLPI problem */
   SCIP_VAR**            nlpivars;           /**< the variables in the NLPI problem */
   SCIP_NLROW**          nlpirows;           /**< the rows in the NLPI problem */
   int                   nlpinvars;          /**< the number of variables in the NPLI problem */
   int                   nlpinrows;          /**< the number of rows in the NLPI problem */
   int                   nlpinslackvars;     /**< the number of slack variables in the NLPI problem */
   int                   nlpiprobsubprob;    /**< the index of the subproblem that the nonlinear problem belongs to */

   SCIP_Real*            slackvarlbs;        /**< an array of zeros for the slack variable lower bounds*/
   SCIP_Real*            slackvarubs;        /**< an array of infinity for the slack variable upper bounds*/
   int*                  slackvarinds;       /**< array of indices for the slack variables */
};

/*
 * Local methods
 */

/** frees the non linear problem */
static
SCIP_RETCODE freeNonlinearProblem(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERSCUT*      benderscut          /**< the Benders' decomposition structure */
   )
{
   SCIP_BENDERSCUTDATA* benderscutdata;

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benderscut != NULL);

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);

   if( benderscutdata->nlpiprob != NULL )
   {
      assert(benderscutdata->nlpi != NULL);

      SCIPfreeBlockMemoryArray(masterprob, &benderscutdata->slackvarinds, benderscutdata->nlpinvars);
      SCIPfreeBlockMemoryArray(masterprob, &benderscutdata->slackvarubs, benderscutdata->nlpinvars);
      SCIPfreeBlockMemoryArray(masterprob, &benderscutdata->slackvarlbs, benderscutdata->nlpinvars);
      SCIPfreeBlockMemoryArray(masterprob, &benderscutdata->nlpirows, benderscutdata->nlpinrows);
      SCIPfreeBlockMemoryArray(masterprob, &benderscutdata->nlpivars, benderscutdata->nlpinvars);
      SCIPhashmapFree(&benderscutdata->row2idx);
      SCIPhashmapFree(&benderscutdata->var2idx);

      SCIP_CALL( SCIPfreeNlpiProblem(subproblem, benderscutdata->nlpi, &benderscutdata->nlpiprob) );

      benderscutdata->nlpinslackvars = 0;
      benderscutdata->nlpinrows = 0;
      benderscutdata->nlpinvars = 0;

      benderscutdata->nlpi = NULL;
   }

   return SCIP_OKAY;
}

/** solves the auxiliary feasibility subproblem.
 *
 *  @note: the variable fixings need to be setup before calling this function
 */
static
SCIP_RETCODE solveFeasibilityNonlinearSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUTDATA*  benderscutdata,     /**< Benders' cut data */
   SCIP_Bool*            success             /**< returns whether solving the feasibility problem was successful */
   )
{
   SCIP_NLPSOLSTAT nlpsolstat;

   assert(scip != NULL);
   assert(benderscutdata != NULL);

   (*success) = TRUE;

   SCIP_CALL( SCIPsolveNlpi(scip, benderscutdata->nlpi, benderscutdata->nlpiprob, .iterlimit = 3000) );  /*lint !e666*/
   SCIPdebugMsg(scip, "NLP solstat = %d\n", SCIPgetNlpiSolstat(scip, benderscutdata->nlpi, benderscutdata->nlpiprob));

   nlpsolstat = SCIPgetNlpiSolstat(scip, benderscutdata->nlpi, benderscutdata->nlpiprob);

   /* if the feasibility NLP is not feasible, then it is not possible to generate a Benders' cut. This is also an error,
    * since the NLP should always be feasible. In debug mode, an ABORT will be thrown.
    */
   if( nlpsolstat > SCIP_NLPSOLSTAT_FEASIBLE )
      (*success) = FALSE;

   return SCIP_OKAY;
}

/** builds the non-linear problem to resolve to generate a cut for the infeasible subproblem */
static
SCIP_RETCODE createAuxiliaryNonlinearSubproblem(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERSCUT*      benderscut          /**< the benders' decomposition cut method */
   )
{
   SCIP_BENDERSCUTDATA* benderscutdata;
   SCIP_Real* obj;
   int i;

   assert(masterprob != NULL);

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);

   /* first freeing the non-linear problem if it exists */
   SCIP_CALL( freeNonlinearProblem(masterprob, subproblem, benderscut) );

   assert(benderscutdata->nlpi == NULL);
   assert(benderscutdata->nlpiprob == NULL);

   benderscutdata->nlpinvars = SCIPgetNVars(subproblem);
   benderscutdata->nlpinrows = SCIPgetNNLPNlRows(subproblem);
   benderscutdata->nlpi = SCIPgetNlpis(subproblem)[0];
   assert(benderscutdata->nlpi != NULL);

   SCIP_CALL( SCIPhashmapCreate(&benderscutdata->var2idx, SCIPblkmem(masterprob), benderscutdata->nlpinvars) );
   SCIP_CALL( SCIPhashmapCreate(&benderscutdata->row2idx, SCIPblkmem(masterprob), benderscutdata->nlpinrows) );

   SCIP_CALL( SCIPduplicateBlockMemoryArray(masterprob, &benderscutdata->nlpivars, SCIPgetVars(subproblem),
         benderscutdata->nlpinvars) ); /*lint !e666*/
   SCIP_CALL( SCIPduplicateBlockMemoryArray(masterprob, &benderscutdata->nlpirows, SCIPgetNLPNlRows(subproblem),
         benderscutdata->nlpinrows) ); /*lint !e666*/

   SCIP_CALL( SCIPcreateNlpiProblemFromNlRows(subproblem, benderscutdata->nlpi, &benderscutdata->nlpiprob, "benders-feascutalt-nlp",
         SCIPgetNLPNlRows(subproblem), benderscutdata->nlpinrows, benderscutdata->var2idx, benderscutdata->row2idx, NULL, SCIPinfinity(subproblem), FALSE,
         FALSE) );

   /* storing the slack variable bounds and indices */
   SCIP_CALL( SCIPallocBufferArray(masterprob, &obj, benderscutdata->nlpinvars) );

   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &benderscutdata->slackvarlbs, benderscutdata->nlpinvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &benderscutdata->slackvarubs, benderscutdata->nlpinvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(masterprob, &benderscutdata->slackvarinds, benderscutdata->nlpinvars) );
   benderscutdata->nlpinslackvars = 0;
   for( i = 0; i < benderscutdata->nlpinvars; i++ )
   {
      if( strstr(SCIPvarGetName(benderscutdata->nlpivars[i]), SLACKVAR_NAME) )
      {
         benderscutdata->slackvarlbs[benderscutdata->nlpinslackvars] = 0.0;
         benderscutdata->slackvarubs[benderscutdata->nlpinslackvars] = SCIPinfinity(subproblem);
         benderscutdata->slackvarinds[benderscutdata->nlpinslackvars] = SCIPhashmapGetImageInt(benderscutdata->var2idx,
            (void*)benderscutdata->nlpivars[i]);

         obj[benderscutdata->nlpinslackvars] = 1.0;

         benderscutdata->nlpinslackvars++;
      }
   }

   /* setting the objective function */
   SCIP_CALL( SCIPsetNlpiObjective(subproblem, benderscutdata->nlpi, benderscutdata->nlpiprob, benderscutdata->nlpinslackvars,
         benderscutdata->slackvarinds, obj, NULL, 0.0) );

   /* unfixing the slack variables */
   SCIP_CALL( SCIPchgNlpiVarBounds(subproblem, benderscutdata->nlpi, benderscutdata->nlpiprob, benderscutdata->nlpinslackvars,
         benderscutdata->slackvarinds, benderscutdata->slackvarlbs, benderscutdata->slackvarubs) );

   SCIPfreeBufferArray(masterprob, &obj);

   return SCIP_OKAY;
}

/** updates the non-linear problem that is resolved to generate a cut for the infeasible subproblem */
static
SCIP_RETCODE updateAuxiliaryNonlinearSubproblem(
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERSCUT*      benderscut          /**< the benders' decomposition cut method */
   )
{
   SCIP_BENDERSCUTDATA* benderscutdata;

   assert(subproblem != NULL);
   assert(benderscut != NULL);

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);
   assert(benderscutdata->nlpi != NULL);
   assert(benderscutdata->nlpiprob != NULL);
   assert(benderscutdata->var2idx != NULL);
   assert(benderscutdata->row2idx != NULL);

   /* setting the variable bounds to that from the current subproblem */
   SCIP_CALL( SCIPupdateNlpiProblem(subproblem, benderscutdata->nlpi, benderscutdata->nlpiprob, benderscutdata->var2idx,
         benderscutdata->nlpivars, benderscutdata->nlpinvars, SCIPinfinity(subproblem)) );

   /* unfixing the slack variables */
   SCIP_CALL( SCIPchgNlpiVarBounds(subproblem, benderscutdata->nlpi, benderscutdata->nlpiprob, benderscutdata->nlpinslackvars,
         benderscutdata->slackvarinds, benderscutdata->slackvarlbs, benderscutdata->slackvarubs) );

   return SCIP_OKAY;
}

/** generates and applies Benders' cuts */
static
SCIP_RETCODE generateAndApplyBendersCuts(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< the benders' decomposition cut method */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the pricing problem */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_RESULT*          result              /**< the result from solving the subproblems */
   )
{
   SCIP_BENDERSCUTDATA* benderscutdata;
   SCIP_Real* primalvals;
   SCIP_Real* consdualvals;
   SCIP_Real* varlbdualvals;
   SCIP_Real* varubdualvals;
   SCIP_Real obj;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Bool success;
#ifdef SCIP_EVENMOREDEBUG
   int i;
#endif

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(result != NULL);

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);

   /* creating or updating the NLPI problem */
   if( benderscutdata->nlpiprob == NULL || benderscutdata->nlpiprobsubprob != probnumber )
   {
      SCIP_CALL( createAuxiliaryNonlinearSubproblem(masterprob, subproblem, benderscut) );
      benderscutdata->nlpiprobsubprob = probnumber;
   }
   else
   {
      SCIP_CALL( updateAuxiliaryNonlinearSubproblem(subproblem, benderscut) );
   }

   /* solving the NLPI problem to get the minimum infeasible solution */
   SCIP_CALL( solveFeasibilityNonlinearSubproblem(subproblem, benderscutdata, &success) );

   if( !success )
   {
      (*result) = SCIP_DIDNOTFIND;
      SCIPdebugMsg(masterprob, "Error in generating Benders' feasibility cut for problem %d. "
         "The feasibility subproblem failed to solve with a feasible solution.\n", probnumber);
      return SCIP_OKAY;
   }

   /* getting the solution from the NLPI problem */
   SCIP_CALL( SCIPgetNlpiSolution(subproblem, benderscutdata->nlpi, benderscutdata->nlpiprob, &primalvals, &consdualvals,
         &varlbdualvals, &varubdualvals, &obj) );

#ifdef SCIP_EVENMOREDEBUG
   SCIPdebugMsg(masterprob, "NLP Feasibility problem solution.\n");
   SCIPdebugMsg(masterprob, "Objective: %g.\n", obj);
   for( i = 0; i < benderscutdata->nlpinvars; i++ )
   {
      int varindex;
      SCIP_Real solval;
      if( SCIPhashmapExists(benderscutdata->var2idx, benderscutdata->nlpivars[i]) )
      {
         varindex = SCIPhashmapGetImageInt(benderscutdata->var2idx, benderscutdata->nlpivars[i]);
         solval = primalvals[varindex];

         if( !SCIPisZero(masterprob, solval) )
         {
            SCIPdebugMsg(masterprob, "%s (obj: %g): %20g\n", SCIPvarGetName(benderscutdata->nlpivars[i]),
               SCIPvarGetObj(benderscutdata->nlpivars[i]), solval);
         }
      }
   }
#endif

   /* setting the name of the generated cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "altfeasibilitycut_%d_%" SCIP_LONGINT_FORMAT, probnumber,
      SCIPbenderscutGetNFound(benderscut) );

   /* generating a Benders' decomposition cut using the classical optimality cut methods */
   SCIP_CALL( SCIPgenerateAndApplyBendersOptCut(masterprob, subproblem, benders, benderscut,
         sol, probnumber, cutname, obj, primalvals, consdualvals, varlbdualvals, varubdualvals, benderscutdata->row2idx,
         benderscutdata->var2idx, type, FALSE, TRUE, result) );

   if( (*result) == SCIP_CONSADDED )
   {
      if( SCIPisInfinity(masterprob, -SCIPgetDualbound(masterprob))
         && SCIPbenderscutGetNFound(benderscut) % SCIP_DEFAULT_DISPLAYFREQ == 0 )
      {
         if( SCIPgetStage(masterprob) == SCIP_STAGE_SOLVING )
         {
            SCIP_CALL( SCIPprintDisplayLine(masterprob, NULL, SCIP_VERBLEVEL_NORMAL, TRUE) );
            SCIPverbMessage(masterprob, SCIP_VERBLEVEL_NORMAL, NULL,
               "Benders' Decomposition: Master problem LP is infeasible. Added %" SCIP_LONGINT_FORMAT " feasibility cuts.\n",
               SCIPbenderscutGetNFound(benderscut));
         }
      }
      SCIPdebugMsg(masterprob, "Constraint <%s> has been added to the master problem.\n", cutname);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of Benders' decomposition cuts
 */

/** deinitialization method of Benders' decomposition cuts (called before transformed problem is freed) */
static
SCIP_DECL_BENDERSCUTEXIT(benderscutExitFeasalt)
{  /*lint --e{715}*/
   assert( benderscut != NULL );
   assert( strcmp(SCIPbenderscutGetName(benderscut), BENDERSCUT_NAME) == 0 );

   return SCIP_OKAY;
}

/** destructor of the Benders' decomposition cut to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BENDERSCUTFREE(benderscutFreeFeasalt)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   assert(scip != NULL);
   assert(benderscut != NULL);

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);

   SCIPfreeBlockMemory(scip, &benderscutdata);

   return SCIP_OKAY;
}

/** execution method of Benders' decomposition cuts */
static
SCIP_DECL_BENDERSCUTEXEC(benderscutExecFeasalt)
{  /*lint --e{715}*/
   SCIP* subproblem;
   SCIP_Bool nlprelaxation;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(benderscut != NULL);
   assert(result != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* setting a flag to indicate whether the NLP relaxation should be used to generate cuts */
   nlprelaxation = SCIPisNLPConstructed(subproblem) && SCIPgetNNlpis(subproblem)
      && SCIPbendersGetSubproblemType(benders, probnumber) <= SCIP_BENDERSSUBTYPE_CONVEXDIS;

   /* only generate feasibility cuts if the subproblem LP or NLP is infeasible,
    * since we use the farkas proof from the LP or the dual solution of the NLP to construct the feasibility cut
    */
   if( SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING &&
       (nlprelaxation && (SCIPgetNLPSolstat(subproblem) == SCIP_NLPSOLSTAT_LOCINFEASIBLE || SCIPgetNLPSolstat(subproblem) == SCIP_NLPSOLSTAT_GLOBINFEASIBLE)) )
   {
      /* generating a cut for a given subproblem */
      SCIP_CALL( generateAndApplyBendersCuts(scip, subproblem, benders, benderscut, sol, probnumber, type, result) );

      /* TODO this was in benderscutExitFeasalt, but freeNonlinearProblem now needs subproblem, which didn't seem to be easily available there */
      /* freeing the non-linear problem information */
      SCIP_CALL( freeNonlinearProblem(scip, subproblem, benderscut) );
   }

   return SCIP_OKAY;
}


/*
 * Benders' decomposition cuts specific interface methods
 */

/** creates the Alternative Feasibility Benders' decomposition cuts and includes it in SCIP */
SCIP_RETCODE SCIPincludeBenderscutFeasalt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_BENDERSCUT* benderscut;
   SCIP_BENDERSCUTDATA* benderscutdata;

   assert(benders != NULL);

   benderscut = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &benderscutdata) );
   BMSclearMemory(benderscutdata);
   benderscutdata->nlpiprobsubprob = -1;

   /* include Benders' decomposition cuts */
   SCIP_CALL( SCIPincludeBenderscutBasic(scip, benders, &benderscut, BENDERSCUT_NAME, BENDERSCUT_DESC,
         BENDERSCUT_PRIORITY, BENDERSCUT_LPCUT, benderscutExecFeasalt, benderscutdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBenderscutFree(scip, benderscut, benderscutFreeFeasalt) );
   SCIP_CALL( SCIPsetBenderscutExit(scip, benderscut, benderscutExitFeasalt) );

   assert(benderscut != NULL);

   return SCIP_OKAY;
}
