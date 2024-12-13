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

/**@file   scip_sol.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

#include "blockmemshell/memory.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/debug.h"
#include "scip/lp.h"
#include "scip/nlp.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/pub_cons.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/relax.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/struct_lp.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "xml/xml.h"

/** checks solution for feasibility in original problem without adding it to the solution store; to improve the
 *  performance we use the following order when checking for violations:
 *
 *  1. variable bounds
 *  2. constraint handlers with positive or zero priority that don't need constraints (e.g. integral constraint handler)
 *  3. original constraints
 *  4. constraint handlers with negative priority that don't need constraints (e.g. Benders' decomposition constraint handler)
 */
static
SCIP_RETCODE checkSolOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked if printreason is true? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             checkmodifiable     /**< have modifiable constraint to be checked? */
   )
{
   SCIP_RESULT result;
   int v;
   int c;
   int h;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(feasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "checkSolOrig", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   *feasible = TRUE;

   SCIPsolResetViolations(sol);

   if( !printreason )
      completely = FALSE;

   /* check bounds */
   if( checkbounds )
   {
      for( v = 0; v < scip->origprob->nvars; ++v )
      {
         SCIP_VAR* var;
         SCIP_Real solval;
         SCIP_Real lb;
         SCIP_Real ub;

         var = scip->origprob->vars[v];
         solval = SCIPsolGetVal(sol, scip->set, scip->stat, var);

         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);

         SCIPupdateSolBoundViolation(scip, sol, lb - solval, SCIPrelDiff(lb, solval));
         SCIPupdateSolBoundViolation(scip, sol, solval - ub, SCIPrelDiff(solval, ub));

         if( SCIPsetIsFeasLT(scip->set, solval, lb) || SCIPsetIsFeasGT(scip->set, solval, ub) )
         {
            *feasible = FALSE;

            if( printreason )
            {
               SCIPmessagePrintInfo(scip->messagehdlr, "solution violates original bounds of variable <%s> [%g,%g] solution value <%g>\n",
                  SCIPvarGetName(var), lb, ub, solval);
            }

            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   /* call constraint handlers with positive or zero check priority that don't need constraints */
   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      if( SCIPconshdlrGetCheckPriority(scip->set->conshdlrs[h]) >= 0 )
      {
         if( !SCIPconshdlrNeedsCons(scip->set->conshdlrs[h]) )
         {
            SCIP_CALL( SCIPconshdlrCheck(scip->set->conshdlrs[h], scip->mem->probmem, scip->set, scip->stat, sol,
                  checkintegrality, checklprows, printreason, completely, &result) );

            if( result != SCIP_FEASIBLE )
            {
               *feasible = FALSE;

               if( !completely )
                  return SCIP_OKAY;
            }
         }
      }
      /* constraint handlers are sorted by priority, so we can break when reaching the first one with negative priority */
      else
         break;
   }

   /* check original constraints
    *
    * in general modifiable constraints can not be checked, because the variables to fulfill them might be missing in
    * the original problem; however, if the solution comes from a heuristic during presolving modifiable constraints
    * have to be checked;
    */
   for( c = 0; c < scip->origprob->nconss; ++c )
   {
      if( SCIPconsIsChecked(scip->origprob->conss[c]) && (checkmodifiable || !SCIPconsIsModifiable(scip->origprob->conss[c])) )
      {
         /* check solution */
         SCIP_CALL( SCIPconsCheck(scip->origprob->conss[c], scip->set, sol,
               checkintegrality, checklprows, printreason, &result) );

         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;

            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   /* call constraint handlers with negative check priority that don't need constraints;
    * continue with the first constraint handler with negative priority which caused us to break in the above loop */
   for( ; h < scip->set->nconshdlrs; ++h )
   {
      assert(SCIPconshdlrGetCheckPriority(scip->set->conshdlrs[h]) < 0);
      if( !SCIPconshdlrNeedsCons(scip->set->conshdlrs[h]) )
      {
         SCIP_CALL( SCIPconshdlrCheck(scip->set->conshdlrs[h], scip->mem->probmem, scip->set, scip->stat, sol,
               checkintegrality, checklprows, printreason, completely, &result) );

         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;

            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** update integrality violation of a solution */
void SCIPupdateSolIntegralityViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol             /**< absolute violation */
   )
{
   if( SCIPprimalUpdateViolations(scip->origprimal) )
      SCIPsolUpdateIntegralityViolation(sol, absviol);
}

/** update bound violation of a solution */
void SCIPupdateSolBoundViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation */
   SCIP_Real             relviol             /**< relative violation */
   )
{
   if( SCIPprimalUpdateViolations(scip->origprimal) )
      SCIPsolUpdateBoundViolation(sol, absviol, relviol);
}

/** update LP row violation of a solution */
void SCIPupdateSolLPRowViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation */
   SCIP_Real             relviol             /**< relative violation */
   )
{
   if( SCIPprimalUpdateViolations(scip->origprimal) )
      SCIPsolUpdateLPRowViolation(sol, absviol, relviol);
}

/** update constraint violation of a solution */
void SCIPupdateSolConsViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation */
   SCIP_Real             relviol             /**< relative violation */
   )
{
   if( SCIPprimalUpdateViolations(scip->origprimal) )
      SCIPsolUpdateConsViolation(sol, absviol, relviol);
}

/** update LP row and constraint violations of a solution */
void SCIPupdateSolLPConsViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation */
   SCIP_Real             relviol             /**< relative violation */
   )
{
   if( SCIPprimalUpdateViolations(scip->origprimal) )
      SCIPsolUpdateLPConsViolation(sol, absviol, relviol);
}

/** allow violation updates */
void SCIPactivateSolViolationUpdates(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIPprimalSetUpdateViolations(scip->origprimal, TRUE);
}

/** disallow violation updates */
void SCIPdeactivateSolViolationUpdates(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIPprimalSetUpdateViolations(scip->origprimal, FALSE);
}

/** creates a primal solution, initialized to zero
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPsolCreateOriginal(sol, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->origprimal, NULL, heur) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPsolCreate(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, heur) );
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDDATA;
   }  /*lint !e788*/
}

/** creates a primal solution, initialized to the current LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("LP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolCreateLPSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
         scip->tree, scip->lp, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateNLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateNLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPerrorMessage("NLP does not exist\n");
      return SCIP_INVALIDCALL;
   }
   assert(scip->nlp != NULL);

   if( !SCIPnlpHasSolution(scip->nlp) )
   {
      SCIPerrorMessage("NLP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolCreateNLPSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, scip->nlp,
         heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current relaxation solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateRelaxSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateRelaxSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPrelaxationIsSolValid(scip->relaxation) )
   {
      SCIPerrorMessage("relaxation solution is not valid\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolCreateRelaxSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, scip->relaxation, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreatePseudoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreatePseudoSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolCreatePseudoSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
         scip->tree, scip->lp, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current LP or pseudo solution, depending on whether the LP was solved
 *  at the current node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolCreateCurrentSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
         scip->tree, scip->lp, heur) );

   return SCIP_OKAY;
}

/** creates a partial primal solution, initialized to unknown values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPcreatePartialSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreatePartialSol", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolCreatePartial(sol, scip->mem->probmem, scip->set, scip->stat, scip->origprimal, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to unknown values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateUnknownSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateUnknownSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolCreateUnknown(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution living in the original problem space, initialized to zero;
 *  a solution in original space allows to set original variables to values that would be invalid in the
 *  transformed problem due to preprocessing fixings or aggregations
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcreateOrigSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateOrigSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPsolCreateOriginal(sol, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->origprimal, NULL, heur) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      SCIP_CALL( SCIPsolCreateOriginal(sol, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->primal, scip->tree, heur) );
      return SCIP_OKAY;

   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** creates a copy of a primal solution; note that a copy of a linked solution is also linked and needs to be unlinked
 *  if it should stay unaffected from changes in the LP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_FREETRANS
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcreateSolCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateSolCopy", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* check if we want to copy the current solution, which is the same as creating a current solution */
   if( sourcesol == NULL )
   {
      SCIP_CALL( SCIPcreateCurrentSol(scip, sol, NULL) );
   }
   else
   {
      SCIP_CALL( SCIPsolCopy(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, sourcesol) );
   }

   return SCIP_OKAY;
}

/** creates a copy of a solution in the original primal solution space
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPcreateSolCopyOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateSolCopyOrig", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   /* check if we want to copy the current solution, which is the same as creating a current solution */
   if( sourcesol == NULL )
   {
      SCIP_CALL( SCIPcreateCurrentSol(scip, sol, NULL) );
   }
   else
   {
      switch( scip->set->stage )
      {
      case SCIP_STAGE_PROBLEM:
      case SCIP_STAGE_FREETRANS:
      case SCIP_STAGE_SOLVED:
      case SCIP_STAGE_TRANSFORMING:
      case SCIP_STAGE_TRANSFORMED:
      case SCIP_STAGE_INITPRESOLVE:
      case SCIP_STAGE_PRESOLVING:
      case SCIP_STAGE_EXITPRESOLVE:
      case SCIP_STAGE_PRESOLVED:
      case SCIP_STAGE_INITSOLVE:
      case SCIP_STAGE_SOLVING:
         SCIP_CALL( SCIPsolCopy(sol, scip->mem->probmem, scip->set, scip->stat, scip->origprimal, sourcesol) );
         break;
      default:
         assert(FALSE); /*lint !e506*/
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** helper method that sets up and solves the sub-SCIP for removing infinite values from solutions */
static
SCIP_RETCODE setupAndSolveFiniteSolSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure of sub-SCIP*/
   SCIP_VAR**            origvars,           /**< original problem variables of main SCIP */
   int                   norigvars,          /**< number of original problem variables of main SCIP */
   SCIP_Real*            solvals,            /**< array with solution values of variables; infinite ones are replaced */
   SCIP_Bool*            success             /**< pointer to store if removing infinite values was successful */
   )
{
   SCIP_HASHMAP* varmap;
   SCIP_VAR* varcopy;
   SCIP_Real fixval;
   SCIP_Bool valid;
   SCIP_SOL* bestsol;
   int v;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(origvars != NULL);
   assert(solvals != NULL);
   assert(success != NULL);

   /* copy the original problem to the sub-SCIP */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), norigvars) );
   SCIP_CALL( SCIPcopyOrig(scip, subscip, varmap, NULL, "removeinffixings", TRUE, FALSE, TRUE, &valid) );

   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );

   /* in the sub-SCIP, we try to minimize the absolute values of all variables with infinite values in the solution
    * and fix all other variables to the value they have in the solution
    */
   for( v = 0; v < norigvars; ++v )
   {
      varcopy = (SCIP_VAR*) SCIPhashmapGetImage(varmap, (void*)origvars[v]);
      assert(varcopy != NULL);

      fixval = solvals[v];

      if( SCIPisInfinity(scip, fixval) || SCIPisInfinity(scip, -fixval) )
      {
         /* If a variable with a finite finite lower bound was set to +infinity, we just change its objective to 1.0
          * to minimize its value; if a variable with a finite finite upper bound was set to -infinity, we just
          * change its objective to -1.0 to maximize its value; if a variable is free, we split the variable into
          * positive and negative part by creating two new non-negative variables and one constraint linking those
          * variables.
          */
         if( SCIPisInfinity(scip, fixval) && !SCIPisInfinity(scip, -SCIPvarGetLbLocal(varcopy)) )
         {
            SCIP_CALL( SCIPchgVarObj(subscip, varcopy, 1.0) );
         }
         else if( SCIPisInfinity(scip, -fixval) && !SCIPisInfinity(scip, SCIPvarGetUbLocal(varcopy)) )
         {
            SCIP_CALL( SCIPchgVarObj(subscip, varcopy, -1.0) );
         }
         else
         {
            char name[SCIP_MAXSTRLEN];
            SCIP_VAR* posvar;
            SCIP_VAR* negvar;
            SCIP_CONS* linkcons;

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s", SCIPvarGetName(varcopy), "run");
            SCIP_CALL( SCIPcreateVar(subscip, &posvar, name, 0.0, SCIPinfinity(scip), 1.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(subscip, posvar) );

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s", SCIPvarGetName(varcopy), "neg");
            SCIP_CALL( SCIPcreateVar(subscip, &negvar, name, 0.0, SCIPinfinity(scip), 1.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(subscip, negvar) );

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s", SCIPvarGetName(varcopy), "linkcons");
            SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &linkcons, name, 0, NULL, NULL, 0.0, 0.0 ) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, linkcons, varcopy, 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, linkcons, posvar, -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, linkcons, negvar, 1.0) );
            SCIP_CALL( SCIPaddCons(subscip, linkcons) );

            SCIP_CALL( SCIPreleaseCons(subscip, &linkcons) );
            SCIP_CALL( SCIPreleaseVar(subscip, &posvar) );
            SCIP_CALL( SCIPreleaseVar(subscip, &negvar) );

            SCIP_CALL( SCIPchgVarObj(subscip, varcopy, 0.0) );
         }
      }
      else
      {
         SCIP_Bool infeasible;
         SCIP_Bool fixed;

         if( SCIPisFeasLT(scip, solvals[v], SCIPvarGetLbLocal(varcopy)) || SCIPisFeasGT(scip, solvals[v], SCIPvarGetUbLocal(varcopy)) )
         {
            SCIP_CALL( SCIPchgVarType(subscip, varcopy, SCIP_VARTYPE_CONTINUOUS, &infeasible) );
            assert(!infeasible);
         }

         /* fix variable to its value in the solution */
         SCIP_CALL( SCIPfixVar(subscip, varcopy, fixval, &infeasible, &fixed) );
         assert(!infeasible);
      }
   }

   SCIP_CALL( SCIPsolve(subscip) );

   bestsol = SCIPgetBestSol(subscip);

   if( bestsol != NULL )
   {
      /* change the stored solution values for variables fixed to infinite values */
      for( v = 0; v < norigvars; ++v )
      {
         varcopy = (SCIP_VAR*) SCIPhashmapGetImage(varmap, (void*)origvars[v]);
         assert(varcopy != NULL);

         if( (SCIPisInfinity(scip, solvals[v]) || SCIPisInfinity(scip, -solvals[v])) )
         {
            solvals[v] = SCIPgetSolVal(subscip, bestsol, varcopy);
         }
      }
   }
   else
   {
      *success = FALSE;
   }

   SCIPhashmapFree(&varmap);

   return SCIP_OKAY;
}


/** creates a copy of a primal solution, thereby replacing infinite fixings of variables by finite values;
 *  the copy is always defined in the original variable space;
 *  success indicates whether the objective value of the solution was changed by removing infinite values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPcreateFiniteSolCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_SOL*             sourcesol,          /**< primal CIP solution to copy */
   SCIP_Bool*            success             /**< does the finite solution have the same objective value? */
   )
{
   SCIP_VAR** fixedvars;
   SCIP_VAR** origvars;
   SCIP_Real* solvals;
   SCIP_VAR* var;
   int nfixedvars;
   int norigvars;
   int v;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateFiniteSolCopy", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   assert(scip != NULL);
   assert(sol != NULL);
   assert(sourcesol != NULL);
   assert(success != NULL);

   *success = TRUE;
   *sol = NULL;

   fixedvars = SCIPgetFixedVars(scip);
   nfixedvars = SCIPgetNFixedVars(scip);
   assert(fixedvars != NULL || nfixedvars == 0);

   /* get original variables and their values in the optimal solution */
   SCIP_CALL( SCIPgetOrigVarsData(scip, &origvars, &norigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, norigvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sourcesol, norigvars, origvars, solvals) );

   /* check whether there are variables fixed to an infinite value */
   for( v = 0; v < nfixedvars; ++v )
   {
      var = fixedvars[v]; /*lint !e613*/

      /* skip (multi-)aggregated variables */
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_FIXED )
         continue;

      assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)));

      if( (SCIPisInfinity(scip, SCIPvarGetLbGlobal(var)) || SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var))) )
      {
         SCIPdebugMsg(scip, "var <%s> is fixed to infinite value %g\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
         break;
      }
   }

   /* there were variables fixed to infinite values */
   if( v < nfixedvars )
   {
      SCIP* subscip;
      SCIP_RETCODE retcode;

      /* if one of the variables was fixed to infinity in the original problem, we stop here */
      for( v = 0; v < norigvars; ++v )
      {
         var = origvars[v];

         if( SCIPisInfinity(scip, SCIPvarGetLbOriginal(var)) || SCIPisInfinity(scip, -SCIPvarGetUbOriginal(var)) )
         {
            assert(SCIPisEQ(scip, SCIPvarGetLbOriginal(var), SCIPvarGetUbOriginal(var)));

            SCIPdebugMsg(scip, "--> var <%s> is fixed to infinite value %g in the original problem, stop making solution finite\n",
               SCIPvarGetName(var), SCIPvarGetLbOriginal(var));

            *success = FALSE;

            goto TERMINATE;
         }
      }

      /* create sub-SCIP */
      SCIP_CALL( SCIPcreate(&subscip) );

      retcode = setupAndSolveFiniteSolSubscip(scip, subscip, origvars, norigvars, solvals, success);

      /* free sub-SCIP */
      SCIP_CALL( SCIPfree(&subscip) );

      SCIP_CALL( retcode );
   }

   /* create original solution and set the solution values */
   if( *success )
   {
      SCIP_CALL( SCIPcreateOrigSol(scip, sol, NULL) );
      for( v = 0; v < norigvars; ++v )
      {
         SCIP_CALL( SCIPsetSolVal(scip, *sol, origvars[v], solvals[v]) );
      }
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "created finites solution copy:\n");
   SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif

   /* the solution of the sub-SCIP should have the same objective value */
   if( *success && !SCIPisEQ(scip, SCIPgetSolOrigObj(scip, *sol), SCIPgetSolOrigObj(scip, sourcesol)) )
   {
      /* @todo how should we avoid numerical trobles here for large objective values? */
      if( (SCIPgetSolOrigObj(scip, *sol) / SCIPepsilon(scip)) < 1e+15 ||
         REALABS(SCIPgetSolOrigObj(scip, *sol) - SCIPgetSolOrigObj(scip, sourcesol)) > 1e-12 * SCIPgetSolOrigObj(scip, *sol) )
         *success = FALSE;
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}

/** frees primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPfreeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol                 /**< pointer to the solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPsolFree(sol, scip->mem->probmem, scip->origprimal) );
      break;
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      SCIP_CALL( SCIPsolFree(sol, scip->mem->probmem, scip->primal) );
      break;
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** links a primal solution to the current LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPlinkLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPlinkLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolved(scip->lp) )
   {
      SCIPerrorMessage("LP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolLinkLPSol(sol, scip->set, scip->stat, scip->transprob, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** links a primal solution to the current NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPlinkNLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPlinkNLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP does not exist\n");
      return SCIP_INVALIDCALL;
   }

   if( SCIPnlpGetSolstat(scip->nlp) > SCIP_NLPSOLSTAT_FEASIBLE )
   {
      SCIPerrorMessage("NLP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolLinkNLPSol(sol, scip->stat, scip->tree, scip->nlp) );

   return SCIP_OKAY;
}

/** links a primal solution to the current relaxation solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPlinkRelaxSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPlinkRelaxSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPrelaxationIsSolValid(scip->relaxation) )
   {
      SCIPerrorMessage("relaxation solution is not valid\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolLinkRelaxSol(sol, scip->set, scip->stat, scip->tree, scip->relaxation) );

   return SCIP_OKAY;
}

/** links a primal solution to the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPlinkPseudoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPlinkPseudoSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolLinkPseudoSol(sol, scip->set, scip->stat, scip->transprob, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** links a primal solution to the current LP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPlinkCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPlinkCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolLinkCurrentSol(sol, scip->set, scip->stat, scip->transprob, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** clears a primal solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPclearSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPclearSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPsolClear(sol, scip->stat, scip->tree) );

   return SCIP_OKAY;
}

/** stores solution values of variables in solution's own array
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPunlinkSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPunlinkSol", FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPsolUnlink(sol, scip->set, scip->transprob) );

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPsetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Real             val                 /**< solution value of variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSolVal", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert( var->scip == scip );

   if( SCIPsolIsOriginal(sol) && SCIPvarIsTransformed(var) )
   {
      SCIPerrorMessage("cannot set value of transformed variable <%s> in original space solution\n",
         SCIPvarGetName(var));
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolSetVal(sol, scip->set, scip->stat, scip->tree, var, val) );

   return SCIP_OKAY;
}

/** sets values of multiple variables in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPsetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   int                   nvars,              /**< number of variables to set solution value for */
   SCIP_VAR**            vars,               /**< array with variables to add to solution */
   SCIP_Real*            vals                /**< array with solution values of variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSolVals", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPsolIsOriginal(sol) )
   {
      for( v = 0; v < nvars; ++v )
      {
         if( SCIPvarIsTransformed(vars[v]) )
         {
            SCIPerrorMessage("cannot set value of transformed variable <%s> in original space solution\n",
               SCIPvarGetName(vars[v]));
            return SCIP_INVALIDCALL;
         }
      }
   }

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPsolSetVal(sol, scip->set, scip->stat, scip->tree, vars[v], vals[v]) );
   }

   return SCIP_OKAY;
}

/** increases value of variable in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPincSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Real             incval              /**< increment for solution value of variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPincSolVal", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert( var->scip == scip );

   if( SCIPsolIsOriginal(sol) && SCIPvarIsTransformed(var) )
   {
      SCIPerrorMessage("cannot increase value of transformed variable <%s> in original space solution\n",
         SCIPvarGetName(var));
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolIncVal(sol, scip->set, scip->stat, scip->tree, var, incval) );

   return SCIP_OKAY;
}

/** returns value of variable in primal CIP solution, or in current LP/pseudo solution
 *
 *  @return value of variable in primal CIP solution, or in current LP/pseudo solution
 *
 *  @pre In case the solution pointer @p sol is @b NULL, that means it is asked for the LP or pseudo solution, this method
 *       can only be called if @p scip is in the solving stage \ref SCIP_STAGE_SOLVING. In any other case, this method
 *       can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Real SCIPgetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolVal", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert( var->scip == scip );

   if( sol != NULL )
      return SCIPsolGetVal(sol, scip->set, scip->stat, var);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolVal(sol==NULL)", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPvarGetSol(var, SCIPtreeHasCurrentNodeLP(scip->tree));
}

/** gets values of multiple variables in primal CIP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPgetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   )
{
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetSolVals", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( sol != NULL )
   {
      int v;

      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPsolGetVal(sol, scip->set, scip->stat, vars[v]);
   }
   else
   {
      SCIP_CALL( SCIPgetVarSols(scip, nvars, vars, vals) );
   }

   return SCIP_OKAY;
}

/** returns objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value
 *
 *  @return objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Real SCIPgetSolOrigObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   )
{
   /* for original solutions, an original objective value is already available in SCIP_STAGE_PROBLEM
    * for all other solutions, we should be at least in SCIP_STAGE_TRANSFORMING
    */
   if( sol != NULL && SCIPsolIsOriginal(sol) )
   {
      SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolOrigObj", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

      return SCIPsolGetOrigObj(sol);
   }

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolOrigObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( sol != NULL )
      return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPsolGetObj(sol, scip->set, scip->transprob, scip->origprob));
   else
   {
      SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolOrigObj(sol==NULL)", \
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      if( SCIPtreeHasCurrentNodeLP(scip->tree) )
         return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPlpGetObjval(scip->lp, scip->set, scip->transprob));
      else
         return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPlpGetPseudoObjval(scip->lp, scip->set, scip->transprob));
   }
}

/** returns transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value
 *
 *  @return transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Real SCIPgetSolTransObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolTransObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( sol != NULL )
      return SCIPsolGetObj(sol, scip->set, scip->transprob, scip->origprob);
   else
   {
      SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolTransObj(sol==NULL)", \
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      if( SCIPtreeHasCurrentNodeLP(scip->tree) )
         return SCIPlpGetObjval(scip->lp, scip->set, scip->transprob);
      else
         return SCIPlpGetPseudoObjval(scip->lp, scip->set, scip->transprob);
   }
}

/** recomputes the objective value of an original solution, e.g., when transferring solutions
 *  from the solution pool (objective coefficients might have changed in the meantime)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 */
SCIP_RETCODE SCIPrecomputeSolObj(
   SCIP*                 scip,
   SCIP_SOL*             sol
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPrecomputeSolObj", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPsolRecomputeObj(sol, scip->set, scip->stat, scip->origprob);

   return SCIP_OKAY;
}

/** maps original space objective value into transformed objective value
 *
 *  @return transformed objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPtransformObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             obj                 /**< original space objective value to transform */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPtransformObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobInternObjval(scip->transprob, scip->origprob, scip->set, obj);
}

/** maps transformed objective value into original space
 *
 *  @return objective value into original space
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPretransformObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             obj                 /**< transformed objective value to retransform in original space */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPretransformObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, obj);
}

/** gets clock time, when this solution was found
 *
 *  @return clock time, when this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Real SCIPgetSolTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolTime", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsolGetTime(sol);
}

/** gets branch and bound run number, where this solution was found
 *
 *  @return branch and bound run number, where this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
int SCIPgetSolRunnum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolRunnum", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsolGetRunnum(sol);
}

/** gets node number of the specific branch and bound run, where this solution was found
 *
 *  @return node number of the specific branch and bound run, where this solution was found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Longint SCIPgetSolNodenum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolNodenum", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsolGetNodenum(sol);
}

/** gets heuristic, that found this solution (or NULL if it's from the tree)
 *
 *  @return heuristic, that found this solution (or NULL if it's from the tree)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_HEUR* SCIPgetSolHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolHeur", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsolGetHeur(sol);
}

/** returns whether two given solutions are exactly equal
 *
 *  @return returns whether two given solutions are exactly equal
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Bool SCIPareSolsEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2                /**< second primal CIP solution */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPareSolsEqual", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsolsAreEqual(sol1, sol2, scip->set, scip->stat, scip->origprob, scip->transprob);
}

/** adjusts solution values of implicit integer variables in handed solution. Solution objective value is not
 *  deteriorated by this method.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPadjustImplicitSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             uselprows           /**< should LP row information be considered for none-objective variables */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPadjustImplicitSolVals", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(sol != NULL);
   SCIP_CALL( SCIPsolAdjustImplicitSolVals(sol, scip->set, scip->stat, scip->transprob, scip->tree, uselprows) );

   return SCIP_OKAY;
}

/** outputs non-zero variables of solution in original problem space to the given file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre In case the solution pointer @p sol is NULL (asking for the current LP/pseudo solution), this method can be
 *       called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @pre In case the solution pointer @p sol is @b not NULL, this method can be called if @p scip is in one of the
 *       following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPprintSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Real objvalue;
   SCIP_Bool currentsol;
   SCIP_Bool oldquiet = FALSE;

   assert(SCIPisTransformed(scip) || sol != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   currentsol = (sol == NULL);
   if( currentsol )
   {
      SCIP_CALL( SCIPcheckStage(scip, "SCIPprintSol(sol==NULL)", \
            FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

      /* create a temporary solution that is linked to the current solution */
      SCIP_CALL( SCIPsolCreateCurrentSol(&sol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
            scip->tree, scip->lp, NULL) );
   }

   if( file != NULL && scip->messagehdlr != NULL )
   {
      oldquiet = SCIPmessagehdlrIsQuiet(scip->messagehdlr);
      SCIPmessagehdlrSetQuiet(scip->messagehdlr, FALSE);
   }

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "objective value:                 ");

   if( SCIPsolIsPartial(sol) )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "unknown\n");
   }
   else
   {
      if( SCIPsolIsOriginal(sol) )
         objvalue = SCIPsolGetOrigObj(sol);
      else
         objvalue = SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPsolGetObj(sol, scip->set, scip->transprob, scip->origprob));

      SCIPprintReal(scip, file, objvalue, 20, 15);
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");
   }

   SCIP_CALL( SCIPsolPrint(sol, scip->set, scip->messagehdlr, scip->stat, scip->origprob, scip->transprob, file, FALSE,
         printzeros) );

   if( file != NULL && scip->messagehdlr != NULL )
   {
      SCIPmessagehdlrSetQuiet(scip->messagehdlr, oldquiet);
   }

   if( currentsol )
   {
      /* free temporary solution */
      SCIP_CALL( SCIPsolFree(&sol, scip->mem->probmem, scip->primal) );
   }

   return SCIP_OKAY;
}

/** outputs non-zero variables of solution in transformed problem space to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPprintTransSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Bool currentsol;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintTransSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   currentsol = (sol == NULL);
   if( currentsol )
   {
      /* create a temporary solution that is linked to the current solution */
      SCIP_CALL( SCIPsolCreateCurrentSol(&sol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
            scip->tree, scip->lp, NULL) );
   }

   if( SCIPsolIsOriginal(sol) )
   {
      SCIPerrorMessage("cannot print original space solution as transformed solution\n");
      return SCIP_INVALIDCALL;
   }

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "objective value:                 ");
   SCIPprintReal(scip, file, SCIPsolGetObj(sol, scip->set, scip->transprob, scip->origprob), 20, 9);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");

   SCIP_CALL( SCIPsolPrint(sol, scip->set, scip->messagehdlr, scip->stat, scip->transprob, NULL, file, FALSE, printzeros) );

   if( currentsol )
   {
      /* free temporary solution */
      SCIP_CALL( SCIPsolFree(&sol, scip->mem->probmem, scip->primal) );
   }

   return SCIP_OKAY;
}

/** outputs discrete variables of solution in original problem space to the given file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPprintMIPStart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_Real objvalue;
   SCIP_Bool oldquiet = FALSE;

   assert(sol != NULL);
   assert(!SCIPsolIsPartial(sol));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintMIPStart", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( file != NULL && scip->messagehdlr != NULL )
   {
      oldquiet = SCIPmessagehdlrIsQuiet(scip->messagehdlr);
      SCIPmessagehdlrSetQuiet(scip->messagehdlr, FALSE);
   }

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "objective value:                 ");

   if( SCIPsolIsOriginal(sol) )
      objvalue = SCIPsolGetOrigObj(sol);
   else
      objvalue = SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPsolGetObj(sol, scip->set, scip->transprob, scip->origprob));

   SCIPprintReal(scip, file, objvalue, 20, 15);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");

   SCIP_CALL( SCIPsolPrint(sol, scip->set, scip->messagehdlr, scip->stat, scip->origprob, scip->transprob, file, TRUE,
         TRUE) );

   if( file != NULL && scip->messagehdlr != NULL )
   {
      SCIPmessagehdlrSetQuiet(scip->messagehdlr, oldquiet);
   }

   return SCIP_OKAY;
}

/** returns dual solution value of a constraint */
SCIP_RETCODE SCIPgetDualSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which the dual solution should be returned */
   SCIP_Real*            dualsolval,         /**< pointer to store the dual solution value */
   SCIP_Bool*            boundconstraint     /**< pointer to store whether the constraint is a bound constraint (or NULL) */
   )
{
   SCIP_CONS* transcons;
   int nvars;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(dualsolval != NULL);

   assert(SCIPconsGetHdlr(cons) != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear" ) == 0);

   SCIP_CALL( SCIPconsGetNVars(cons, scip->set, &nvars, &success) );
   assert(success);  /* is always successful, since we only have linear constraints */

   if( boundconstraint != NULL )
      *boundconstraint = (nvars == 1);

   if( SCIPconsIsTransformed(cons) )
      transcons = cons;
   else
      transcons = SCIPconsGetTransformed(cons);

   /* it can happen that a transformed constraints gets deleted due to redundancy. by complementary slackness the
    * corresponding dual solution value would be zero. however, if the constraint contains exactly one variable we need
    * to check the reduced costs of the variable.
    */
   if( nvars == 0 || (nvars > 1 && transcons == NULL) )
      (*dualsolval) = 0.0;
   else
   {
      if( nvars > 1 )
         (*dualsolval) = SCIPgetDualsolLinear(scip, transcons);
      else
      {
         /* the constraint is a bound constraint */
         SCIP_VAR** vars;
         SCIP_Real* vals;
         SCIP_Real activity;

         vars = SCIPgetVarsLinear(scip, cons);
         vals = SCIPgetValsLinear(scip, cons);

         activity = SCIPvarGetLPSol(vars[0]) * vals[0];

         /* return the reduced cost of the variable if the constraint would be tight */
         if( SCIPsetIsEQ(scip->set, activity, SCIPgetRhsLinear(scip, cons))
          || SCIPsetIsEQ(scip->set, activity, SCIPgetLhsLinear(scip, cons)) )
            (*dualsolval) = SCIPgetVarRedcost(scip, vars[0]);
         else
            (*dualsolval) = 0.0;
      }
   }
   assert(*dualsolval != SCIP_INVALID); /*lint !e777*/

   /* dual values are coming from the LP solver that is always solving a minimization problem */
   if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE )
      (*dualsolval) *= -1.0;

   return SCIP_OKAY;
}

/** outputs dual solution from LP solver to file stream */
static
SCIP_RETCODE printDualSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Bool boundconstraint;
   int c;

   assert(scip->lp != NULL);
   assert(scip->lp->solved);
   assert(scip->lp->dualfeasible);

   /* print dual solution values of all constraints */
   for( c = 0; c < scip->origprob->nconss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_Real solval;

      cons = scip->origprob->conss[c];
      assert(cons != NULL);

      SCIP_CALL( SCIPgetDualSolVal(scip, cons, &solval, &boundconstraint) );

      if( printzeros || !SCIPisZero(scip, solval) )
      {
         SCIP_MESSAGEHDLR* messagehdlr = scip->messagehdlr;

         SCIPmessageFPrintInfo(messagehdlr, file, "%-32s", SCIPconsGetName(cons));

         if( SCIPisInfinity(scip, solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            +infinity\n");
         else if( SCIPisInfinity(scip, -solval) )
            SCIPmessageFPrintInfo(messagehdlr, file, "            -infinity\n");
         else
         {
            if( boundconstraint )
               SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g*\n", solval);
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " %20.15g\n", solval);
         }
      }
   }

   return SCIP_OKAY;
}

/** check whether the dual solution is available
 *
 * @note This is used when calling \ref SCIPprintDualSol()
 *
 * @return is dual solution available?
 *
 * @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Bool SCIPisDualSolAvailable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             printreason         /**< print warning message if dualsol is not available? */
   )
{
   int c;

   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisDualSolAvailable", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVED )
   {
      if( printreason )
         SCIPmessageFPrintInfo(scip->messagehdlr, NULL, "No dual solution available.\n");
      return FALSE;
   }

   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);

   /* dual solution only useful when no presolving was performed */
   if( scip->stat->performpresol )
   {
      if( printreason )
         SCIPwarningMessage(scip, "No dual information available when presolving was performed.\n");
      return FALSE;
   }

   /* dual solution is created by LP solver and therefore only available for pure LPs */
   if( scip->transprob->nvars != scip->transprob->ncontvars )
   {
      if( printreason )
         SCIPwarningMessage(scip, "Dual information only available for pure LPs (only continuous variables).\n");
      return FALSE;
   }

   /* dual solution is created by LP solver and therefore only available for linear constraints */
   for( c = scip->transprob->nconss - 1; c >= 0; --c )
   {
      SCIP_CONSHDLR* conshdlr;

      conshdlr = SCIPconsGetHdlr(scip->transprob->conss[c]);
      assert(conshdlr != NULL);

      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear" ) != 0 )
      {
         if( printreason )
            SCIPwarningMessage(scip, "Dual information only available for pure LPs (only linear constraints).\n");
         return FALSE;
      }
   }

   return TRUE;
}

/** outputs dual solution from LP solver to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called in all stages but only prints dual information when called in \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPprintDualSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   if( SCIPisDualSolAvailable(scip, TRUE) )
   {
      /* print dual solution */
      SCIP_CALL( printDualSol(scip, file, printzeros) );
   }

   return SCIP_OKAY;
}


/** outputs non-zero variables of solution representing a ray in original problem space to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPprintRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution representing ray */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   assert(scip != NULL);
   assert(sol != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintRay", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolPrintRay(sol, scip->set, scip->messagehdlr, scip->stat, scip->origprob, scip->transprob, file, printzeros) );

   return SCIP_OKAY;
}

/** gets number of feasible primal solutions stored in the solution storage in case the problem is transformed;
 *  in case the problem stage is SCIP_STAGE_PROBLEM, the number of solution in the original solution candidate
 *  storage is returned
 *
 *  @return number of feasible primal solutions stored in the solution storage in case the problem is transformed; or
 *          number of solution in the original solution candidate storage if the problem stage is SCIP_STAGE_PROBLEM
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNSols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNSols", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprimal->nsols;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      return scip->primal->nsols;

   case SCIP_STAGE_INIT:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return -1; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets array of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *  if the problem stage is in SCIP_STAGE_PROBLEM, it returns the number array of solution candidate stored
 *
 *  @return array of feasible primal solutions
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_SOL** SCIPgetSols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSols", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprimal->sols;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      return scip->primal->sols;

   case SCIP_STAGE_INIT:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return NULL;
   }  /*lint !e788*/
}

/** gets best feasible primal solution found so far if the problem is transformed; in case the problem is in
 *  SCIP_STAGE_PROBLEM it returns the best solution candidate, or NULL if no solution has been found or the candidate
 *  store is empty;
 *
 *  @return best feasible primal solution so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_SOL* SCIPgetBestSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetBestSol", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
      return NULL;
   case SCIP_STAGE_PROBLEM:
      assert(scip->origprimal != NULL);
      if(  scip->origprimal->nsols > 0 )
      {
         assert(scip->origprimal->sols != NULL);
         assert(scip->origprimal->sols[0] != NULL);
         return scip->origprimal->sols[0];
      }
      break;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
      assert(scip->primal != NULL);
      if(  scip->primal->nsols > 0 )
      {
         assert(scip->primal->sols != NULL);
         assert(scip->primal->sols[0] != NULL);
         return scip->primal->sols[0];
      }
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return NULL;
   }

   return NULL;
}

/** outputs best feasible primal solution found so far to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPprintBestSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_SOL* sol;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintBestSol", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   sol = SCIPgetBestSol(scip);

   if( sol == NULL )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "no solution available\n");
   else
   {
      SCIP_CALL( SCIPprintSol(scip, sol, file, printzeros) );
   }

   return SCIP_OKAY;
}

/** outputs best feasible primal solution found so far in transformed variables to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPprintBestTransSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_SOL* sol;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintBestTransSol", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   sol = SCIPgetBestSol(scip);

   if( sol != NULL && SCIPsolIsOriginal(sol) )
   {
      SCIPerrorMessage("best solution is defined in original space - cannot print it as transformed solution\n");
      return SCIP_INVALIDCALL;
   }

   if( sol == NULL )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "no solution available\n");
   else
   {
      SCIP_CALL( SCIPprintTransSol(scip, sol, file, printzeros) );
   }

   return SCIP_OKAY;
}

/** try to round given solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIProundSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Bool*            success             /**< pointer to store whether rounding was successful */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIProundSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPsolIsOriginal(sol) )
   {
      SCIPerrorMessage("cannot round original space solution\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolRound(sol, scip->set, scip->stat, scip->transprob, scip->tree, success) );

   return SCIP_OKAY;
}

/** retransforms solution to original problem space
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPretransformSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPretransformSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch ( SCIPsolGetOrigin(sol) )
   {
   case SCIP_SOLORIGIN_ORIGINAL:
      /* nothing to do */
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_LPSOL:
   case SCIP_SOLORIGIN_NLPSOL:
   case SCIP_SOLORIGIN_RELAXSOL:
   case SCIP_SOLORIGIN_PSEUDOSOL:

      /* first unlink solution */
      SCIP_CALL( SCIPunlinkSol(scip, sol) );

      /*lint -fallthrough*/
   case SCIP_SOLORIGIN_ZERO:
   {
      SCIP_Bool hasinfval;

      SCIP_CALL( SCIPsolRetransform(sol, scip->set, scip->stat, scip->origprob, scip->transprob, &hasinfval) );
      break;
   }
   case SCIP_SOLORIGIN_PARTIAL:
   case SCIP_SOLORIGIN_UNKNOWN:
      SCIPerrorMessage("unknown solution origin.\n");
      return SCIP_INVALIDCALL;

   default:
      /* note that this is in an internal SCIP error since all solution origins are covert in the switch above */
      SCIPerrorMessage("invalid solution origin <%d>\n", SCIPsolGetOrigin(sol));
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** reads a given solution file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPreadSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPreadSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* we pass the reading of the solution file on to reader_sol via the following call */
   SCIP_CALL( SCIPreadProb(scip, filename, "sol") );

   return SCIP_OKAY;
}

/** reads a given solution file and store the solution values in the given solution pointer */
static
SCIP_RETCODE readSolFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of the input file */
   SCIP_SOL*             sol,                /**< solution pointer */
   SCIP_Bool*            partial,            /**< pointer to store if the solution is partial (or NULL, if not needed) */
   SCIP_Bool*            error               /**< pointer store if an error occured */
   )
{
   SCIP_FILE* file;
   SCIP_Bool unknownvariablemessage;
   SCIP_Bool localpartial;
   int lineno;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(error != NULL);

   /* open input file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   *error = FALSE;
   localpartial = SCIPsolIsPartial(sol);

   unknownvariablemessage = FALSE;
   lineno = 0;

   /* read the file */
   while( !SCIPfeof(file) && !(*error) )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char valuestring[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      char format[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real value;
      int nread;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* there are some lines which may preceed the solution information */
      if( strncasecmp(buffer, "solution status:", 16) == 0 || strncasecmp(buffer, "objective value:", 16) == 0 ||
         strncasecmp(buffer, "Log started", 11) == 0 || strncasecmp(buffer, "Variable Name", 13) == 0 ||
         strncasecmp(buffer, "All other variables", 19) == 0 || strncasecmp(buffer, "\n", 1) == 0 ||
         strncasecmp(buffer, "NAME", 4) == 0 || strncasecmp(buffer, "ENDATA", 6) == 0 )    /* allow parsing of SOL-format on the MIPLIB 2003 pages */
         continue;

      /* parse the line */
      (void) SCIPsnprintf(format, SCIP_MAXSTRLEN, "%%%ds %%%ds %%%ds\n", SCIP_MAXSTRLEN, SCIP_MAXSTRLEN, SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, varname, valuestring, objstring);
      if( nread < 2 )
      {
         SCIPerrorMessage("Invalid input line %d in solution file <%s>: <%s>.\n", lineno, filename, buffer);
         *error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "unknown variable <%s> in line %d of solution file <%s>\n",
               varname, lineno, filename);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* cast the value */
      if( strncasecmp(valuestring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(valuestring, "+inf", 4) == 0 || strncasecmp(valuestring, "inf", 3) == 0 )
         value = SCIPinfinity(scip);
      else if( strncasecmp(valuestring, "-inf", 4) == 0 )
         value = -SCIPinfinity(scip);
      else if( strncasecmp(valuestring, "unknown", 7) == 0 )
      {
         value = SCIP_UNKNOWN;
         localpartial = TRUE;
      }
      else
      {
         /* coverity[secure_coding] */
         nread = sscanf(valuestring, "%lf", &value);
         if( nread != 1 )
         {
            SCIPerrorMessage("Invalid solution value <%s> for variable <%s> in line %d of solution file <%s>.\n",
               valuestring, varname, lineno, filename);
            *error = TRUE;
            break;
         }
      }

      /* set the solution value of the variable, if not multiaggregated */
      if( SCIPisTransformed(scip) && SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_MULTAGGR )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored solution value for multiaggregated variable <%s>\n", SCIPvarGetName(var));
      }
      else
      {
         SCIP_RETCODE retcode;

         retcode = SCIPsetSolVal(scip, sol, var, value);

         if( retcode == SCIP_INVALIDDATA )
         {
            if( SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_FIXED )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored conflicting solution value for fixed variable <%s>\n",
                  SCIPvarGetName(var));
            }
            else
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored solution value for multiaggregated variable <%s>\n",
                  SCIPvarGetName(var));
            }
         }
         else
         {
            SCIP_CALL_FINALLY( retcode, SCIPfclose(file) );
         }
      }
   }

   /* close input file */
   SCIPfclose(file);

   if( localpartial && !SCIPsolIsPartial(sol) )
   {
      if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
      {
         SCIP_CALL( SCIPsolMarkPartial(sol, scip->set, scip->stat, scip->origprob->vars, scip->origprob->nvars) );
      }
      else
         *error = TRUE;
   }

   if( partial != NULL )
      *partial = localpartial;

   return SCIP_OKAY;
}

/** reads a given xml solution file and store the solution values in the given solution pointer */
static
SCIP_RETCODE readXmlSolFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of the input file */
   SCIP_SOL*             sol,                /**< solution pointer */
   SCIP_Bool*            partial,            /**< pointer to store if the solution is partial (or NULL if not needed) */
   SCIP_Bool*            error               /**< pointer store if an error occured */
   )
{
   SCIP_Bool unknownvariablemessage;
   SCIP_Bool localpartial;
   XML_NODE* start;
   const XML_NODE* varsnode;
   const XML_NODE* varnode;
   const char* tag;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(error != NULL);

   /* read xml file */
   start = xmlProcess(filename);

   if( start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing the XML solution file.\n");
      return SCIP_READERROR;
   }

   *error = FALSE;
   localpartial = SCIPsolIsPartial(sol);

   /* find variable sections */
   tag = "variables";
   varsnode = xmlFindNodeMaxdepth(start, tag, 0, 3);
   if( varsnode == NULL )
   {
      /* free xml data */
      xmlFreeNode(start);

      SCIPerrorMessage("Variable section not found.\n");
      return SCIP_READERROR;
   }

   /* loop through all variables */
   unknownvariablemessage = FALSE;
   for( varnode = xmlFirstChild(varsnode); varnode != NULL; varnode = xmlNextSibl(varnode) )
   {
      SCIP_VAR* var;
      const char* varname;
      const char* valuestring;
      SCIP_Real value;
      int nread;

      /* find variable name */
      varname = xmlGetAttrval(varnode, "name");
      if( varname == NULL )
      {
         SCIPerrorMessage("Attribute \"name\" of variable not found.\n");
         *error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "unknown variable <%s> of solution file <%s>\n",
               varname, filename);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* find value of variable */
      valuestring = xmlGetAttrval(varnode, "value");
      if( valuestring == NULL )
      {
         SCIPerrorMessage("Attribute \"value\" of variable not found.\n");
         *error = TRUE;
         break;
      }

      /* cast the value */
      if( strncasecmp(valuestring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(valuestring, "+inf", 4) == 0 || strncasecmp(valuestring, "inf", 3) == 0 )
         value = SCIPinfinity(scip);
      else if( strncasecmp(valuestring, "-inf", 4) == 0 )
         value = -SCIPinfinity(scip);
      else if( strncasecmp(valuestring, "unknown", 7) == 0 )
      {
         value = SCIP_UNKNOWN;
         localpartial = TRUE;
      }
      else
      {
         /* coverity[secure_coding] */
         nread = sscanf(valuestring, "%lf", &value);
         if( nread != 1 )
         {
            SCIPwarningMessage(scip, "invalid solution value <%s> for variable <%s> in XML solution file <%s>\n", valuestring, varname, filename);
            *error = TRUE;
            break;
         }
      }

      /* set the solution value of the variable, if not multiaggregated */
      if( SCIPisTransformed(scip) && SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_MULTAGGR )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored solution value for multiaggregated variable <%s>\n", SCIPvarGetName(var));
      }
      else
      {
         SCIP_RETCODE retcode;
         retcode = SCIPsetSolVal(scip, sol, var, value);

         if( retcode == SCIP_INVALIDDATA )
         {
            if( SCIPvarGetStatus(SCIPvarGetProbvar(var)) == SCIP_VARSTATUS_FIXED )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored conflicting solution value for fixed variable <%s>\n",
                  SCIPvarGetName(var));
            }
            else
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ignored solution value for multiaggregated variable <%s>\n",
                  SCIPvarGetName(var));
            }
         }
         else
         {
            SCIP_CALL( retcode );
         }
      }
   }

   /* free xml data */
   xmlFreeNode(start);

   if( localpartial && !SCIPsolIsPartial(sol)  )
   {
      if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
      {
         SCIP_CALL( SCIPsolMarkPartial(sol, scip->set, scip->stat, scip->origprob->vars, scip->origprob->nvars) );
      }
      else
         *error = TRUE;
   }

   if( partial != NULL )
      *partial = localpartial;

   return SCIP_OKAY;
}

/** reads a given solution file and store the solution values in the given solution pointer
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPreadSolFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of the input file */
   SCIP_SOL*             sol,                /**< solution pointer */
   SCIP_Bool             xml,                /**< true, iff the given solution in written in XML */
   SCIP_Bool*            partial,            /**< pointer to store if the solution is partial */
   SCIP_Bool*            error               /**< pointer store if an error occured */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPreadSolFile", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( xml )
   {
      SCIP_CALL( readXmlSolFile(scip, filename, sol, partial, error) );
   }
   else
   {
      SCIP_CALL( readSolFile(scip, filename, sol, partial, error) );
   }

   return SCIP_OKAY;
}

/** adds feasible primal solution to solution storage by copying it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note Do not call during propagation, use heur_trysol instead.
 */
SCIP_RETCODE SCIPaddSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_FREETRANS:
      assert(SCIPsolIsOriginal(sol));
      SCIP_CALL( SCIPprimalAddOrigSol(scip->origprimal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, sol, stored) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
      /* if the solution is added during presolving and it is not defined on original variables,
       * presolving operations will destroy its validity, so we retransform it to the original space
       */
      if( !SCIPsolIsOriginal(sol) )
      {
         SCIP_SOL* bestsol = SCIPgetBestSol(scip);
         SCIP_SOL* tmpsol = sol;
         SCIP_Bool hasinfval;

         SCIP_CALL( SCIPcreateSolCopy(scip, &tmpsol, sol) );

         SCIP_CALL( SCIPsolUnlink(tmpsol, scip->set, scip->transprob) );
         SCIP_CALL( SCIPsolRetransform(tmpsol, scip->set, scip->stat, scip->origprob, scip->transprob, &hasinfval) );

         SCIP_CALL( SCIPprimalAddSolFree(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
               scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter,
               &tmpsol, stored) );

         if( *stored && (bestsol != SCIPgetBestSol(scip)) )
         {
            SCIPstoreSolutionGap(scip);
         }

         return SCIP_OKAY;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   {
      SCIP_SOL* bestsol = SCIPgetBestSol(scip);

      SCIP_CALL( SCIPprimalAddSol(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
            scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter, sol,
            stored) );

      /* @todo use solution index rather than pointer */
      if( *stored && (bestsol != SCIPgetBestSol(scip)) )
      {
         SCIPstoreSolutionGap(scip);
      }

      return SCIP_OKAY;
   }
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** adds primal solution to solution storage, frees the solution afterwards
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note Do not call during propagation, use heur_trysol instead.
 */
SCIP_RETCODE SCIPaddSolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddSolFree", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_FREETRANS:
      assert(SCIPsolIsOriginal(*sol));
      SCIP_CALL( SCIPprimalAddOrigSolFree(scip->origprimal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, sol, stored) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
      /* if the solution is added during presolving and it is not defined on original variables,
       * presolving operations will destroy its validity, so we retransform it to the original space
       */
      if( !SCIPsolIsOriginal(*sol) )
      {
         SCIP_Bool hasinfval;

         SCIP_CALL( SCIPsolUnlink(*sol, scip->set, scip->transprob) );
         SCIP_CALL( SCIPsolRetransform(*sol, scip->set, scip->stat, scip->origprob, scip->transprob, &hasinfval) );
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   {
      SCIP_SOL* bestsol = SCIPgetBestSol(scip);

      SCIP_CALL( SCIPprimalAddSolFree(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
            scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter,
            sol, stored) );

      if( *stored )
      {
         if( bestsol != SCIPgetBestSol(scip) )
         {
            assert(SCIPgetBestSol(scip) != NULL);
            SCIPstoreSolutionGap(scip);
         }
      }

      return SCIP_OKAY;
   }
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** adds current LP/pseudo solution to solution storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_SOL* bestsol;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   bestsol = SCIPgetBestSol(scip);

   SCIP_CALL( SCIPprimalAddCurrentSol(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
         scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter, heur,
         stored) );

   if( *stored )
   {
      if( bestsol != SCIPgetBestSol(scip) )
         SCIPstoreSolutionGap(scip);
   }

   return SCIP_OKAY;
}

/** checks solution for feasibility; if possible, adds it to storage by copying
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note Do not call during propagation, use heur_trysol instead.
 */
SCIP_RETCODE SCIPtrySol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violation be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked if printreason is true? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   SCIP_SOL* bestsol;

   assert(sol != NULL);
   assert(stored != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtrySol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   bestsol = SCIPgetBestSol(scip);

   if( !printreason )
      completely = FALSE;

   /* we cannot check partial solutions */
   if( SCIPsolIsPartial(sol) )
   {
      SCIPerrorMessage("Cannot check feasibility of partial solutions.\n");
      return SCIP_INVALIDDATA;
   }

   /* if the solution is added during presolving and it is not defined on original variables,
    * presolving operations will destroy its validity, so we retransform it to the original space
    */
   if( scip->set->stage == SCIP_STAGE_PRESOLVING && !SCIPsolIsOriginal(sol) )
   {
      SCIP_Bool hasinfval;

      SCIP_CALL( SCIPsolUnlink(sol, scip->set, scip->transprob) );
      SCIP_CALL( SCIPsolRetransform(sol, scip->set, scip->stat, scip->origprob, scip->transprob, &hasinfval) );
   }

   if( SCIPsolIsOriginal(sol) )
   {
      SCIP_Bool feasible;

      /* SCIPprimalTrySol() can only be called on transformed solutions; therefore check solutions in original problem
       * including modifiable constraints */
      SCIP_CALL( checkSolOrig(scip, sol, &feasible, printreason, completely, checkbounds, checkintegrality, checklprows, TRUE) );
      if( feasible )
      {
         SCIP_CALL( SCIPprimalAddSol(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
               scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter,
               sol, stored) );

         if( *stored )
         {
            if( bestsol != SCIPgetBestSol(scip) )
               SCIPstoreSolutionGap(scip);
         }
      }
      else
         *stored = FALSE;
   }
   else
   {
      SCIP_CALL( SCIPprimalTrySol(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->origprob,
            scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter, sol, printreason,
            completely, checkbounds, checkintegrality, checklprows, stored) );

      if( *stored )
      {
         if( bestsol != SCIPgetBestSol(scip) )
         {
#ifdef SCIP_DEBUG_ABORTATORIGINFEAS
            SCIP_Bool feasible;
            SCIP_CALL( checkSolOrig(scip, sol, &feasible, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

            if( ! feasible )
            {
               SCIPerrorMessage("Accepted solution not feasible for original problem\n");
               SCIPABORT();
            }
#endif
            SCIPstoreSolutionGap(scip);
         }
      }
   }

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note Do not call during propagation, use heur_trysol instead.
 */
SCIP_RETCODE SCIPtrySolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed */
   SCIP_Bool             completely,         /**< Should all violations be checked if printreason is true? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   )
{
   SCIP_SOL* bestsol;

   assert(stored != NULL);
   assert(sol != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtrySolFree", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   bestsol = SCIPgetBestSol(scip);

   if( !printreason )
      completely = FALSE;

   /* we cannot check partial solutions */
   if( SCIPsolIsPartial(*sol) )
   {
      SCIPerrorMessage("Cannot check feasibility of partial solutions.\n");
      return SCIP_INVALIDDATA;
   }

   /* if the solution is added during presolving and it is not defined on original variables,
    * presolving operations will destroy its validity, so we retransform it to the original space
    */
   if( scip->set->stage == SCIP_STAGE_PRESOLVING && !SCIPsolIsOriginal(*sol) )
   {
      SCIP_Bool hasinfval;

      SCIP_CALL( SCIPsolUnlink(*sol, scip->set, scip->transprob) );
      SCIP_CALL( SCIPsolRetransform(*sol, scip->set, scip->stat, scip->origprob, scip->transprob, &hasinfval) );
   }

   if( SCIPsolIsOriginal(*sol) )
   {
      SCIP_Bool feasible;

      /* SCIPprimalTrySol() can only be called on transformed solutions; therefore check solutions in original problem
       * including modifiable constraints
       */
      SCIP_CALL( checkSolOrig(scip, *sol, &feasible, printreason, completely, checkbounds, checkintegrality, checklprows, TRUE) );

      if( feasible )
      {
         SCIP_CALL( SCIPprimalAddSolFree(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
               scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter,
               sol, stored) );

         if( *stored )
         {
            if( bestsol != SCIPgetBestSol(scip) )
               SCIPstoreSolutionGap(scip);
         }
      }
      else
      {
         SCIP_CALL( SCIPsolFree(sol, scip->mem->probmem, scip->primal) );
         *stored = FALSE;
      }
   }
   else
   {
      SCIP_CALL( SCIPprimalTrySolFree(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
            scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter,
            sol, printreason, completely, checkbounds, checkintegrality, checklprows, stored) );

      if( *stored )
      {
         if( bestsol != SCIPgetBestSol(scip) )
         {
#ifdef SCIP_DEBUG_ABORTATORIGINFEAS
            SCIP_Bool feasible;
            SCIP_CALL( checkSolOrig(scip, SCIPgetBestSol(scip), &feasible, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

            if( ! feasible )
            {
               SCIPerrorMessage("Accepted incumbent not feasible for original problem\n");
               SCIPABORT();
            }
#endif
            SCIPstoreSolutionGap(scip);
         }
      }
   }

   return SCIP_OKAY;
}

/** checks current LP/pseudo solution for feasibility; if possible, adds it to storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPtryCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked if printreason is true? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   SCIP_SOL* bestsol;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtryCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   bestsol = SCIPgetBestSol(scip);

   if( !printreason )
      completely = FALSE;

   SCIP_CALL( SCIPprimalTryCurrentSol(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
         scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter, heur,
         printreason, completely, checkintegrality, checklprows, stored) );

   if( *stored )
   {
      if( bestsol != SCIPgetBestSol(scip) )
      {
#ifdef SCIP_DEBUG_ABORTATORIGINFEAS
         SCIP_Bool feasible;
         SCIP_CALL( checkSolOrig(scip, SCIPgetBestSol(scip), &feasible, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

         if( ! feasible )
         {
            SCIPerrorMessage("Accepted incumbent not feasible for original problem\n");
            SCIPABORT();
         }
#endif
         SCIPstoreSolutionGap(scip);
      }
   }

   return SCIP_OKAY;
}

/** returns all partial solutions
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_SOL** SCIPgetPartialSols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPartialSols", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->origprimal->partialsols;
}

/** returns number of partial solutions
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNPartialSols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPartialSols", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->origprimal->npartialsols;
}

/** checks solution for feasibility without adding it to the solution store
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcheckSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked if printreason is true? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            feasible            /**< stores whether given solution is feasible */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcheckSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* return immediately if the solution is of type partial */
   if( SCIPsolIsPartial(sol) )
   {
      SCIPerrorMessage("Cannot check feasibility of partial solutions.");
      return SCIP_INVALIDDATA;
   }

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || scip->set->misc_exactsolve;

   if( !printreason )
      completely = FALSE;

   if( SCIPsolIsOriginal(sol) )
   {
      /* SCIPsolCheck() can only be called on transformed solutions */
      SCIP_CALL( checkSolOrig(scip, sol, feasible, printreason, completely, checkbounds, checkintegrality, checklprows, FALSE) );
   }
   else
   {
      SCIP_CALL( SCIPsolCheck(sol, scip->set, scip->messagehdlr, scip->mem->probmem, scip->stat, scip->transprob,
            printreason, completely, checkbounds, checkintegrality, checklprows, feasible) );
   }

   return SCIP_OKAY;
}

/** checks solution for feasibility in original problem without adding it to the solution store;
 *  this method is used to double check a solution in order to validate the presolving process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcheckSolOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool             completely          /**< Should all violations be checked if printreason is true? */
   )
{
   assert(scip != NULL);
   assert(sol != NULL);
   assert(feasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcheckSolOrig", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* return immediately if the solution is of type partial */
   if( SCIPsolIsPartial(sol) )
   {
      SCIPerrorMessage("Cannot check feasibility of partial solutions.");
      return SCIP_INVALIDDATA;
   }

   if( !printreason )
      completely = FALSE;

   /* check solution in original problem; that includes bounds, integrality, and non modifiable constraints */
   SCIP_CALL( checkSolOrig(scip, sol, feasible, printreason, completely, TRUE, TRUE, TRUE, FALSE) );

   return SCIP_OKAY;
}

/** return whether a primal ray is stored that proves unboundedness of the LP relaxation
 *
 *  @return return whether a primal ray is stored that proves unboundedness of the LP relaxation
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Bool SCIPhasPrimalRay(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPhasPrimalRay", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->primal->primalray != NULL;
}

/** gets value of given variable in primal ray causing unboundedness of the LP relaxation;
 *  should only be called if such a ray is stored (check with SCIPhasPrimalRay())
 *
 *  @return value of given variable in primal ray causing unboundedness of the LP relaxation
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetPrimalRayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPrimalRayVal", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert(var != NULL);
   assert(scip->primal->primalray != NULL);
   assert(var->scip == scip);

   return SCIPsolGetRayVal(scip->primal->primalray, scip->set, scip->stat, var);
}

/** updates the primal ray thats proves unboundedness
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPupdatePrimalRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             primalray           /**< the new primal ray */
   )
{
   assert(scip != NULL);
   assert(primalray != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdatePrimalRay", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPprimalUpdateRay(scip->primal, scip->set, scip->stat, primalray, scip->mem->probmem) );

   return SCIP_OKAY;
}
