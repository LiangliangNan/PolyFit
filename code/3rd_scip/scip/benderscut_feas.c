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

/**@file   benderscut_feas.c
 * @ingroup OTHER_CFILES
 * @brief  Standard feasibility cuts for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_expr.h"
#include "scip/benderscut_feas.h"
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
#include "scip/scip_prob.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"

#define BENDERSCUT_NAME             "feas"
#define BENDERSCUT_DESC             "Standard feasibility cuts for Benders' decomposition"
#define BENDERSCUT_PRIORITY     10000
#define BENDERSCUT_LPCUT         TRUE

/*
 * Local methods
 */

/** adds a variable and value to the constraint/row arrays */
static
SCIP_RETCODE addVariableToArray(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_VAR***           vars,               /**< pointer to array of variables in the generated cut with non-zero coefficient */
   SCIP_Real**           vals,               /**< pointer to array of coefficients of the variables in the generated cut */
   SCIP_VAR*             addvar,             /**< the variable that will be added to the array */
   SCIP_Real             addval,             /**< the value that will be added to the array */
   int*                  nvars,              /**< the number of variables in the variable array */
   int*                  varssize            /**< the length of the variable size */
   )
{
   assert(masterprob != NULL);
   assert(vars != NULL);
   assert(*vars != NULL);
   assert(vals != NULL);
   assert(*vals != NULL);
   assert(addvar != NULL);
   assert(nvars != NULL);
   assert(varssize != NULL);

   if( *nvars >= *varssize )
   {
      *varssize = SCIPcalcMemGrowSize(masterprob, *varssize + 1);
      SCIP_CALL( SCIPreallocBufferArray(masterprob, vars, *varssize) );
      SCIP_CALL( SCIPreallocBufferArray(masterprob, vals, *varssize) );
   }
   assert(*nvars < *varssize);

   (*vars)[*nvars] = addvar;
   (*vals)[*nvars] = addval;
   (*nvars)++;

   return SCIP_OKAY;
}

/** computing as standard Benders' feasibility cut from the dual solutions of the LP */
static
SCIP_RETCODE computeStandardLPFeasibilityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_VAR***           vars,               /**< pointer to array of variables in the generated cut with non-zero coefficient */
   SCIP_Real**           vals,               /**< pointer to array of coefficients of the variables in the generated cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   int*                  nvars,              /**< the number of variables in the cut */
   int*                  varssize,           /**< the number of variables in the array */
   SCIP_Bool*            success             /**< was the cut generation successful? */
   )
{
   SCIP_VAR** subvars;
   int nsubvars;
   int nrows;
   SCIP_Real dualsol;
   SCIP_Real addval;    /* the value that must be added to the lhs */
   int i;

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE);

   (*success) = FALSE;

   /* looping over all LP rows and setting the coefficients of the cut */
   nrows = SCIPgetNLPRows(subproblem);
   for( i = 0; i < nrows; i++ )
   {
      SCIP_ROW* lprow;

      lprow = SCIPgetLPRows(subproblem)[i];
      assert(lprow != NULL);

      dualsol = SCIProwGetDualfarkas(lprow);
      assert( !SCIPisInfinity(subproblem, dualsol) && !SCIPisInfinity(subproblem, -dualsol) );

      if( SCIPisDualfeasZero(subproblem, dualsol) )
         continue;

      if( dualsol > 0.0 )
         addval = dualsol*SCIProwGetLhs(lprow);
      else
         addval = dualsol*SCIProwGetRhs(lprow);

      *lhs += addval;

      /* if the bound becomes infinite, then the cut generation terminates. */
      if( SCIPisInfinity(masterprob, *lhs) || SCIPisInfinity(masterprob, -*lhs)
         || SCIPisInfinity(masterprob, addval) || SCIPisInfinity(masterprob, -addval))
      {
         (*success) = FALSE;
         SCIPdebugMsg(masterprob, "Infinite bound when generating feasibility cut.\n");
         return SCIP_OKAY;
      }
   }

   nsubvars = SCIPgetNVars(subproblem);
   subvars = SCIPgetVars(subproblem);

   /* looping over all variables to update the coefficients in the computed cut. */
   for( i = 0; i < nsubvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_VAR* mastervar;

      var = subvars[i];

      /* retrieving the master problem variable for the given subproblem variable. */
      SCIP_CALL( SCIPgetBendersMasterVar(masterprob, benders, var, &mastervar) );

      dualsol = SCIPgetVarFarkasCoef(subproblem, var);

      if( SCIPisZero(subproblem, dualsol) )
         continue;

      /* checking whether the original variable is a linking variable.
       * If this is the case, then the corresponding master variable is added to the generated cut.
       * If the pricing variable is not a linking variable, then the farkas dual value is added to the lhs
       */
      if( mastervar != NULL )
      {
         SCIPdebugMsg(masterprob ,"Adding coeffs to feasibility cut: <%s> dualsol %g\n", SCIPvarGetName(mastervar), dualsol);

         /* adding the variable to the storage */
         SCIP_CALL( addVariableToArray(masterprob, vars, vals, mastervar, dualsol, nvars, varssize) );
      }
      else
      {
         addval = 0;

         if( SCIPisPositive(subproblem, dualsol) )
            addval = dualsol*SCIPvarGetUbGlobal(var);
         else if( SCIPisNegative(subproblem, dualsol) )
            addval = dualsol*SCIPvarGetLbGlobal(var);

         *lhs -= addval;

         /* if the bound becomes infinite, then the cut generation terminates. */
         if( SCIPisInfinity(masterprob, *lhs) || SCIPisInfinity(masterprob, -*lhs)
            || SCIPisInfinity(masterprob, addval) || SCIPisInfinity(masterprob, -addval))
         {
            (*success) = FALSE;
            SCIPdebugMsg(masterprob, "Infinite bound when generating feasibility cut.\n");
            return SCIP_OKAY;
         }
      }
   }

   (*success) = TRUE;

   return SCIP_OKAY;
}


/** computing as standard Benders' feasibility cut from the dual solutions of the NLP
 *
 *  NOTE: The cut must be created before being passed to this function
 */
static
SCIP_RETCODE computeStandardNLPFeasibilityCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_VAR***           vars,               /**< pointer to array of variables in the generated cut with non-zero coefficient */
   SCIP_Real**           vals,               /**< pointer to array of coefficients of the variables in the generated cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   int*                  nvars,              /**< the number of variables in the cut */
   int*                  varssize,           /**< the number of variables in the array */
   SCIP_Bool*            success             /**< was the cut generation successful? */
   )
{
   int nrows;
   SCIP_Real activity;
   SCIP_Real dirderiv;
   SCIP_Real dualsol;
   int i;

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(SCIPisNLPConstructed(subproblem));
   assert(SCIPgetNLPSolstat(subproblem) == SCIP_NLPSOLSTAT_LOCINFEASIBLE || SCIPgetNLPSolstat(subproblem) == SCIP_NLPSOLSTAT_GLOBINFEASIBLE);

   (*success) = FALSE;

   *lhs = 0.0;
   dirderiv = 0.0;

   /* looping over all NLP rows and setting the corresponding coefficients of the cut */
   nrows = SCIPgetNNLPNlRows(subproblem);
   for( i = 0; i < nrows; i++ )
   {
      SCIP_NLROW* nlrow;

      nlrow = SCIPgetNLPNlRows(subproblem)[i];
      assert(nlrow != NULL);

      dualsol = SCIPnlrowGetDualsol(nlrow);
      assert( !SCIPisInfinity(subproblem, dualsol) && !SCIPisInfinity(subproblem, -dualsol) );

      if( SCIPisZero(subproblem, dualsol) )
         continue;

      SCIP_CALL( SCIPaddNlRowGradientBenderscutOpt(masterprob, subproblem, benders, nlrow, -dualsol,
            NULL, NULL, &dirderiv, vars, vals, nvars, varssize) );

      SCIP_CALL( SCIPgetNlRowActivity(subproblem, nlrow, &activity) );

      if( dualsol > 0.0 )
      {
         assert(!SCIPisInfinity(subproblem, SCIPnlrowGetRhs(nlrow)));
         *lhs += dualsol * (activity - SCIPnlrowGetRhs(nlrow));
      }
      else
      {
         assert(!SCIPisInfinity(subproblem, -SCIPnlrowGetLhs(nlrow)));
         *lhs += dualsol * (activity - SCIPnlrowGetLhs(nlrow));
      }
   }

   *lhs += dirderiv;

   /* if the side became infinite or dirderiv was infinite, then the cut generation terminates. */
   if( SCIPisInfinity(masterprob, *lhs) || SCIPisInfinity(masterprob, -*lhs)
      || SCIPisInfinity(masterprob, dirderiv) || SCIPisInfinity(masterprob, -dirderiv))
   {
      (*success) = FALSE;
      SCIPdebugMsg(masterprob, "Infinite bound when generating feasibility cut. lhs = %g dirderiv = %g.\n", *lhs, dirderiv);
      return SCIP_OKAY;
   }

   (*success) = TRUE;

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
   SCIP_RESULT*          result              /**< the result from solving the subproblems */
   )
{
   SCIP_CONS* cut;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real activity;
   int nvars;
   int varssize;
   int nmastervars;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Bool success;

   assert(masterprob != NULL);
   assert(subproblem != NULL);
   assert(benders != NULL);
   assert(result != NULL);

   /* allocating memory for the variable and values arrays */
   nmastervars = SCIPgetNVars(masterprob) + SCIPgetNFixedVars(masterprob);
   SCIP_CALL( SCIPallocClearBufferArray(masterprob, &vars, nmastervars) );
   SCIP_CALL( SCIPallocClearBufferArray(masterprob, &vals, nmastervars) );
   lhs = 0.0;
   nvars = 0;
   varssize = nmastervars;

   /* setting the name of the generated cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "feasibilitycut_%d_%" SCIP_LONGINT_FORMAT, probnumber,
      SCIPbenderscutGetNFound(benderscut) );

   if( SCIPisNLPConstructed(subproblem) && SCIPgetNNlpis(subproblem) )
   {
      /* computing the coefficients of the feasibility cut from the NLP */
      SCIP_CALL( computeStandardNLPFeasibilityCut(masterprob, subproblem, benders, &vars, &vals, &lhs, &nvars, &varssize,
            &success) );
   }
   else
   {
      if( SCIPgetNLPIterations(subproblem) == 0 )
      {
         SCIPverbMessage(masterprob, SCIP_VERBLEVEL_FULL, NULL, "There were no iterations in pricing problem %d. "
           "A Benders' decomposition feasibility cut will be generated from the presolved LP data.\n", probnumber);
      }

      /* computing the coefficients of the feasibility cut from the LP */
      SCIP_CALL( computeStandardLPFeasibilityCut(masterprob, subproblem, benders, &vars, &vals, &lhs, &nvars, &varssize,
            &success) );
   }

   /* if success is FALSE, then there was an error in generating the feasibility cut. No cut will be added to the master
    * problem. Otherwise, the constraint is added to the master problem.
    */
   if( !success )
   {
      (*result) = SCIP_DIDNOTFIND;
      SCIPdebugMsg(masterprob, "Error in generating Benders' feasibility cut for problem %d.\n", probnumber);
   }
   else
   {
      /* creating a constraint with the variables and coefficients previously stored */
      SCIP_CALL( SCIPcreateConsBasicLinear(masterprob, &cut, cutname, nvars, vars, vals, lhs, SCIPinfinity(masterprob)) );
      SCIP_CALL( SCIPsetConsDynamic(masterprob, cut, TRUE) );
      SCIP_CALL( SCIPsetConsRemovable(masterprob, cut, TRUE) );

      assert(SCIPisInfinity(masterprob, SCIPgetRhsLinear(masterprob, cut)));

      /* the activity of the cut should be less than the lhs. This will ensure that the evaluated solution will be cut off.
       * It is possible that the activity is greater than the lhs. This could be caused by numerical difficulties. In this
       * case, no cut will be generated.
       */
      lhs = SCIPgetLhsLinear(masterprob, cut);
      activity = SCIPgetActivityLinear(masterprob, cut, sol);
      if( SCIPisGE(masterprob, activity, lhs) )
      {
         success = FALSE;
         SCIPdebugMsg(masterprob ,"Invalid feasibility cut - activity is greater than lhs %g >= %g.\n", activity, lhs);
#ifdef SCIP_DEBUG
         SCIPABORT();
#endif
      }

      assert(cut != NULL);

      if( success )
      {
         /* adding the constraint to the master problem */
         SCIP_CALL( SCIPaddCons(masterprob, cut) );

         SCIPdebugPrintCons(masterprob, cut, NULL);

         (*result) = SCIP_CONSADDED;
      }

      SCIP_CALL( SCIPreleaseCons(masterprob, &cut) );
   }

   SCIPfreeBufferArray(masterprob, &vals);
   SCIPfreeBufferArray(masterprob, &vars);

   return SCIP_OKAY;
}

/*
 * Callback methods of Benders' decomposition cuts
 */

/** execution method of Benders' decomposition cuts */
static
SCIP_DECL_BENDERSCUTEXEC(benderscutExecFeas)
{  /*lint --e{715}*/
   SCIP* subproblem;
   SCIP_Bool nlprelaxation;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(benderscut != NULL);
   assert(result != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   if( subproblem == NULL )
   {
      SCIPdebugMsg(scip, "The subproblem %d is set to NULL. The <%s> Benders' decomposition cut can not be executed.\n",
         probnumber, BENDERSCUT_NAME);

      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* setting a flag to indicate whether the NLP relaxation should be used to generate cuts */
   nlprelaxation = SCIPisNLPConstructed(subproblem) && SCIPgetNNlpis(subproblem);

   /* only generate feasibility cuts if the subproblem LP or NLP is infeasible,
    * since we use the farkas proof from the LP or the dual solution of the NLP to construct the feasibility cut
    */
   if( SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING &&
      ((!nlprelaxation && SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE) ||
       (nlprelaxation && (SCIPgetNLPSolstat(subproblem) == SCIP_NLPSOLSTAT_LOCINFEASIBLE || SCIPgetNLPSolstat(subproblem) == SCIP_NLPSOLSTAT_GLOBINFEASIBLE))) )
   {
      /* generating a cut for a given subproblem */
      SCIP_CALL( generateAndApplyBendersCuts(scip, subproblem, benders, benderscut,
            sol, probnumber, result) );
   }

   return SCIP_OKAY;
}


/*
 * Benders' decomposition cuts specific interface methods
 */

/** creates the Standard Feasibility Benders' decomposition cuts and includes it in SCIP */
SCIP_RETCODE SCIPincludeBenderscutFeas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_BENDERSCUT* benderscut;

   assert(benders != NULL);

   benderscut = NULL;

   /* include Benders' decomposition cuts */
   SCIP_CALL( SCIPincludeBenderscutBasic(scip, benders, &benderscut, BENDERSCUT_NAME, BENDERSCUT_DESC,
         BENDERSCUT_PRIORITY, BENDERSCUT_LPCUT, benderscutExecFeas, NULL) );

   assert(benderscut != NULL);

   return SCIP_OKAY;
}
