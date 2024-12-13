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

/**@file   scip_solve.c
 * @ingroup OTHER_CFILES
 * @brief  public solving methods
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

#include "blockmemshell/memory.h"
#include "scip/branch.h"
#include "scip/clock.h"
#include "scip/compr.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cutpool.h"
#include "scip/dcmp.h"
#include "scip/debug.h"
#include "scip/event.h"
#include "scip/implics.h"
#include "scip/interrupt.h"
#include "scip/lp.h"
#include "scip/nlp.h"
#include "scip/presol.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/pub_branch.h"
#include "scip/pub_compr.h"
#include "scip/pub_cons.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_select.h"
#include "scip/pub_presol.h"
#include "scip/pub_prop.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/relax.h"
#include "scip/reopt.h"
#include "scip/scip_benders.h"
#include "scip/scip_branch.h"
#include "scip/scip_concurrent.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/struct_event.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_tree.h"
#include "scip/syncstore.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/visual.h"

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
   SCIP_Bool             completely,         /**< Should all violations be checked? */
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

/** calculates number of nonzeros in problem */
static
SCIP_RETCODE calcNonZeros(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint*         nchecknonzeros,     /**< pointer to store number of non-zeros in all check constraints */
   SCIP_Longint*         nactivenonzeros,    /**< pointer to store number of non-zeros in all active constraints */
   SCIP_Bool*            approxchecknonzeros,/**< pointer to store if the number of non-zeros in all check constraints
                                              *   is only a lowerbound
                                              */
   SCIP_Bool*            approxactivenonzeros/**< pointer to store if the number of non-zeros in all active constraints
                                              *   is only a lowerbound
                                              */
   )
{
   SCIP_CONS** conss;
   SCIP_Bool success;
   SCIP_Bool ischeck;
   int nconss;
   int nvars;
   int c;
   int h;

   *nchecknonzeros = 0LL;
   *nactivenonzeros = 0LL;
   *approxchecknonzeros = FALSE;
   *approxactivenonzeros = FALSE;

   /* computes number of non-zeros over all active constraints */
   for( h = scip->set->nconshdlrs - 1; h >= 0; --h )
   {
      nconss = SCIPconshdlrGetNActiveConss(scip->set->conshdlrs[h]);

      if( nconss > 0 )
      {
         conss = SCIPconshdlrGetConss(scip->set->conshdlrs[h]);

         /* calculate all active constraints */
         for( c = nconss - 1; c >= 0; --c )
         {
            SCIP_CALL( SCIPconsGetNVars(conss[c], scip->set, &nvars, &success) );
            ischeck = SCIPconsIsChecked(conss[c]);

            if( !success )
            {
               *approxactivenonzeros = TRUE;
               if( ischeck )
                  *approxchecknonzeros = TRUE;
            }
            else
            {
               *nactivenonzeros += nvars;
               if( ischeck )
                  *nchecknonzeros += nvars;
            }
         }
      }

      /* add nonzeros on inactive check constraints */
      nconss = SCIPconshdlrGetNCheckConss(scip->set->conshdlrs[h]);
      if( nconss > 0 )
      {
         conss = SCIPconshdlrGetCheckConss(scip->set->conshdlrs[h]);

         for( c = nconss - 1; c >= 0; --c )
         {
            if( !SCIPconsIsActive(conss[c]) )
            {
               SCIP_CALL( SCIPconsGetNVars(conss[c], scip->set, &nvars, &success) );

               if( !success )
                  *approxchecknonzeros = TRUE;
               else
                  *nchecknonzeros += nvars;
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** initializes solving data structures and transforms problem
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
 *       - \ref SCIP_STAGE_FREETRANS
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post When calling this method in the \ref SCIP_STAGE_PROBLEM stage, the \SCIP stage is changed to \ref
 *        SCIP_STAGE_TRANSFORMED; otherwise, the stage is not changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPtransformProb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Longint oldnsolsfound;
   int nfeassols;
   int ncandsols;
   int h;
   int s;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtransformProb", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* check, if the problem was already transformed */
   if( scip->set->stage >= SCIP_STAGE_TRANSFORMED )
      return SCIP_OKAY;

   assert(scip->stat->status == SCIP_STATUS_UNKNOWN);

   /* check, if a node selector exists */
   if( SCIPsetGetNodesel(scip->set, scip->stat) == NULL )
   {
      SCIPerrorMessage("no node selector available\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* call garbage collector on original problem and parameter settings memory spaces */
   BMSgarbagecollectBlockMemory(scip->mem->setmem);
   BMSgarbagecollectBlockMemory(scip->mem->probmem);

   /* remember number of constraints */
   SCIPprobMarkNConss(scip->origprob);

   /* switch stage to TRANSFORMING */
   scip->set->stage = SCIP_STAGE_TRANSFORMING;

   /* mark statistics before solving */
   SCIPstatMark(scip->stat);

   /* init solve data structures */
   SCIP_CALL( SCIPeventfilterCreate(&scip->eventfilter, scip->mem->probmem) );
   SCIP_CALL( SCIPeventqueueCreate(&scip->eventqueue) );
   SCIP_CALL( SCIPbranchcandCreate(&scip->branchcand) );
   SCIP_CALL( SCIPlpCreate(&scip->lp, scip->set, scip->messagehdlr, scip->stat, SCIPprobGetName(scip->origprob)) );
   SCIP_CALL( SCIPprimalCreate(&scip->primal) );
   SCIP_CALL( SCIPtreeCreate(&scip->tree, scip->mem->probmem, scip->set, SCIPsetGetNodesel(scip->set, scip->stat)) );
   SCIP_CALL( SCIPrelaxationCreate(&scip->relaxation, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree) );
   SCIP_CALL( SCIPconflictCreate(&scip->conflict, scip->mem->probmem, scip->set) );
   SCIP_CALL( SCIPcliquetableCreate(&scip->cliquetable, scip->set, scip->mem->probmem) );

   /* copy problem in solve memory */
   SCIP_CALL( SCIPprobTransform(scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree,
         scip->reopt, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue, scip->conflictstore,
         &scip->transprob) );

   /* switch stage to TRANSFORMED */
   scip->set->stage = SCIP_STAGE_TRANSFORMED;

   /* check, whether objective value is always integral by inspecting the problem, if it is the case adjust the
    * cutoff bound if primal solution is already known
    */
   SCIP_CALL( SCIPprobCheckObjIntegral(scip->transprob, scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->primal,
	 scip->tree, scip->reopt, scip->lp, scip->eventfilter, scip->eventqueue) );

   /* if possible, scale objective function such that it becomes integral with gcd 1 */
   SCIP_CALL( SCIPprobScaleObj(scip->transprob, scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->primal,
	 scip->tree, scip->reopt, scip->lp, scip->eventfilter, scip->eventqueue) );

   /* check solution of solution candidate storage */
   nfeassols = 0;
   ncandsols = scip->origprimal->nsols;
   oldnsolsfound = 0;

   /* update upper bound and cutoff bound due to objective limit in primal data */
   SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
         scip->eventqueue, scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp) );

   if( !scip->set->reopt_enable && scip->set->nactivebenders == 0 )
   {
      oldnsolsfound = scip->primal->nsolsfound;
      for( s = scip->origprimal->nsols - 1; s >= 0; --s )
      {
         SCIP_Bool feasible;
         SCIP_SOL* sol;

         sol =  scip->origprimal->sols[s];

         /* recompute objective function, since the objective might have changed in the meantime */
         SCIPsolRecomputeObj(sol, scip->set, scip->stat, scip->origprob);

         /* SCIPprimalTrySol() can only be called on transformed solutions; therefore check solutions in original problem
          * including modifiable constraints
          */
         SCIP_CALL( checkSolOrig(scip, sol, &feasible,
               (scip->set->disp_verblevel >= SCIP_VERBLEVEL_HIGH ? scip->set->misc_printreason : FALSE),
               FALSE, TRUE, TRUE, TRUE, TRUE) );

         if( feasible )
         {
            SCIP_Real abssolobj;

            abssolobj = REALABS(SCIPsolGetObj(sol, scip->set, scip->transprob, scip->origprob));

            /* we do not want to add solutions with objective value +infinity */
            if( !SCIPisInfinity(scip, abssolobj) )
            {
               SCIP_SOL* bestsol = SCIPgetBestSol(scip);
               SCIP_Bool stored;

               /* add primal solution to solution storage by copying it */
               SCIP_CALL( SCIPprimalAddSol(scip->primal, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->origprob, scip->transprob,
                     scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter, sol, &stored) );

               if( stored )
               {
                  nfeassols++;

                  if( bestsol != SCIPgetBestSol(scip) )
                     SCIPstoreSolutionGap(scip);
               }
            }
         }

         SCIP_CALL( SCIPsolFree(&sol, scip->mem->probmem, scip->origprimal) );
         scip->origprimal->nsols--;
      }
   }

   assert(scip->origprimal->nsols == 0);

   scip->stat->nexternalsolsfound += scip->primal->nsolsfound - oldnsolsfound;

   if( nfeassols > 0 )
   {
      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "%d/%d feasible solution%s given by solution candidate storage, new primal bound %.6e\n\n",
         nfeassols, ncandsols, (nfeassols > 1 ? "s" : ""), SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip)));
   }
   else if( ncandsols > 0 && !scip->set->reopt_enable )
   {
      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "all %d solutions given by solution candidate storage are infeasible\n\n", ncandsols);
   }

   /* print transformed problem statistics */
   SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
      "transformed problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
      scip->transprob->nvars, scip->transprob->nbinvars, scip->transprob->nintvars, scip->transprob->nimplvars,
      scip->transprob->ncontvars, scip->transprob->nconss);

   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      int nactiveconss;

      nactiveconss = SCIPconshdlrGetNActiveConss(scip->set->conshdlrs[h]);
      if( nactiveconss > 0 )
      {
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "%7d constraints of type <%s>\n", nactiveconss, SCIPconshdlrGetName(scip->set->conshdlrs[h]));
      }
   }
   SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL, "\n");

   {
      SCIP_Real maxnonzeros;
      SCIP_Longint nchecknonzeros;
      SCIP_Longint nactivenonzeros;
      SCIP_Bool approxchecknonzeros;
      SCIP_Bool approxactivenonzeros;

      /* determine number of non-zeros */
      maxnonzeros = (SCIP_Real)SCIPgetNConss(scip) * SCIPgetNVars(scip);
      maxnonzeros = MAX(maxnonzeros, 1.0);
      SCIP_CALL( calcNonZeros(scip, &nchecknonzeros, &nactivenonzeros, &approxchecknonzeros, &approxactivenonzeros) );
      scip->stat->nnz = nactivenonzeros;
      scip->stat->avgnnz = (SCIPgetNConss(scip) == 0 ? 0.0 : (SCIP_Real) nactivenonzeros / ((SCIP_Real) SCIPgetNConss(scip)));

      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "original problem has %s%" SCIP_LONGINT_FORMAT " active (%g%%) nonzeros and %s%" SCIP_LONGINT_FORMAT " (%g%%) check nonzeros\n",
         approxactivenonzeros ? "more than " : "", nactivenonzeros, nactivenonzeros/maxnonzeros * 100,
         approxchecknonzeros ? "more than " : "", nchecknonzeros, nchecknonzeros/maxnonzeros * 100);
      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL, "\n");
   }

   /* call initialization methods of plugins */
   SCIP_CALL( SCIPsetInitPlugins(scip->set, scip->mem->probmem, scip->stat) );

   /* in case the permutation seed is different to 0, permute the transformed problem */
   if( scip->set->random_permutationseed > 0 )
   {
      SCIP_Bool permuteconss;
      SCIP_Bool permutevars;
      int permutationseed;

      permuteconss = scip->set->random_permuteconss;
      permutevars = scip->set->random_permutevars;
      permutationseed = scip->set->random_permutationseed;

      SCIP_CALL( SCIPpermuteProb(scip, (unsigned int)permutationseed, permuteconss, permutevars, permutevars, permutevars, permutevars) );
   }

   if( scip->set->misc_estimexternmem )
   {
      /* the following formula was estimated empirically using linear regression */
      scip->stat->externmemestim = (SCIP_Longint) (MAX(1, 8.5e-04 * SCIPgetNConss(scip) + 7.6e-04 * SCIPgetNVars(scip) + 3.5e-05 * scip->stat->nnz) * 1048576.0); /*lint !e666*/
      SCIPdebugMsg(scip, "external memory usage estimated to %" SCIP_LONGINT_FORMAT " byte\n", scip->stat->externmemestim);
   }

   return SCIP_OKAY;
}

/** initializes presolving */
static
SCIP_RETCODE initPresolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#ifndef NDEBUG
   size_t nusedbuffers;
   size_t nusedcleanbuffers;
#endif

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);
   assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);

   /* retransform all existing solutions to original problem space, because the transformed problem space may
    * get modified in presolving and the solutions may become invalid for the transformed problem
    */
   SCIP_CALL( SCIPprimalRetransformSolutions(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
         scip->eventqueue, scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp) );

   /* reset statistics for presolving and current branch and bound run */
   SCIPstatResetPresolving(scip->stat, scip->set, scip->transprob, scip->origprob);

   /* increase number of branch and bound runs */
   scip->stat->nruns++;

   /* remember problem size of previous run */
   scip->stat->prevrunnvars = scip->transprob->nvars;

   /* switch stage to INITPRESOLVE */
   scip->set->stage = SCIP_STAGE_INITPRESOLVE;

   /* create temporary presolving root node */
   SCIP_CALL( SCIPtreeCreatePresolvingRoot(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->messagehdlr,
         scip->stat, scip->transprob, scip->origprob, scip->primal, scip->lp, scip->branchcand, scip->conflict,
         scip->conflictstore, scip->eventfilter, scip->eventqueue, scip->cliquetable) );

   /* GCG wants to perform presolving during the reading process of a file reader;
    * hence the number of used buffers does not need to be zero, however, it should not
    * change by calling SCIPsetInitprePlugins()
    */
#ifndef NDEBUG
   nusedbuffers = BMSgetNUsedBufferMemory(SCIPbuffer(scip));
   nusedcleanbuffers = BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip));
#endif

   /* inform plugins that the presolving is abound to begin */
   SCIP_CALL( SCIPsetInitprePlugins(scip->set, scip->mem->probmem, scip->stat) );
   assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
   assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

   /* delete the variables from the problems that were marked to be deleted */
   SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->cliquetable, scip->lp, scip->branchcand) );

   /* switch stage to PRESOLVING */
   scip->set->stage = SCIP_STAGE_PRESOLVING;

   return SCIP_OKAY;
}

/** deinitializes presolving */
static
SCIP_RETCODE exitPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             solved,             /**< is problem already solved? */
   SCIP_Bool*            infeasible          /**< pointer to store if the clique clean up detects an infeasibility */
   )
{
#ifndef NDEBUG
   size_t nusedbuffers;
   size_t nusedcleanbuffers;
#endif

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);
   assert(scip->set->stage == SCIP_STAGE_PRESOLVING);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* switch stage to EXITPRESOLVE */
   scip->set->stage = SCIP_STAGE_EXITPRESOLVE;

   if( !solved )
   {
      SCIP_VAR** vars;
      int nvars;
      int v;

      /* flatten all variables */
      vars = SCIPgetFixedVars(scip);
      nvars = SCIPgetNFixedVars(scip);
      assert(nvars == 0 || vars != NULL);

      for( v = nvars - 1; v >= 0; --v )
      {
	 SCIP_VAR* var;
#ifndef NDEBUG
	 SCIP_VAR** multvars;
	 int i;
#endif
	 var = vars[v]; /*lint !e613*/
	 assert(var != NULL);

	 if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
	 {
	    /* flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later-on */
	    SCIP_CALL( SCIPvarFlattenAggregationGraph(var, scip->mem->probmem, scip->set, scip->eventqueue) );

#ifndef NDEBUG
	    multvars = SCIPvarGetMultaggrVars(var);
	    for( i = SCIPvarGetMultaggrNVars(var) - 1; i >= 0; --i)
	       assert(SCIPvarGetStatus(multvars[i]) != SCIP_VARSTATUS_MULTAGGR);
#endif
	 }
      }
   }

   /* exitPresolve() might be called during the reading process of a file reader;
    * hence the number of used buffers does not need to be zero, however, it should not
    * change by calling SCIPsetExitprePlugins() or SCIPprobExitPresolve()
    */
#ifndef NDEBUG
   nusedbuffers = BMSgetNUsedBufferMemory(SCIPbuffer(scip));
   nusedcleanbuffers = BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip));
#endif

   /* inform plugins that the presolving is finished, and perform final modifications */
   SCIP_CALL( SCIPsetExitprePlugins(scip->set, scip->mem->probmem, scip->stat) );
   assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
   assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

   /* remove empty and single variable cliques from the clique table, and convert all two variable cliques
    * into implications
    * delete the variables from the problems that were marked to be deleted
    */
   if( !solved )
   {
      int nlocalbdchgs = 0;

      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
            scip->cliquetable, scip->lp, scip->branchcand) );

      SCIP_CALL( SCIPcliquetableCleanup(scip->cliquetable, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, &nlocalbdchgs,
            infeasible) );

      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "clique table cleanup detected %d bound changes%s\n", nlocalbdchgs, *infeasible ? " and infeasibility" : "");
   }

   /* exit presolving */
   SCIP_CALL( SCIPprobExitPresolve(scip->transprob,  scip->set) );
   assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
   assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

   if( !solved )
   {
      /* check, whether objective value is always integral by inspecting the problem, if it is the case adjust the
       * cutoff bound if primal solution is already known
       */
      SCIP_CALL( SCIPprobCheckObjIntegral(scip->transprob, scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->primal,
	    scip->tree, scip->reopt, scip->lp, scip->eventfilter, scip->eventqueue) );

      /* if possible, scale objective function such that it becomes integral with gcd 1 */
      SCIP_CALL( SCIPprobScaleObj(scip->transprob, scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->primal,
            scip->tree, scip->reopt, scip->lp, scip->eventfilter, scip->eventqueue) );

      scip->stat->lastlowerbound = SCIPprobInternObjval(scip->transprob, scip->origprob, scip->set, scip->transprob->dualbound);

      /* we need to update the primal dual integral here to update the last{upper/dual}bound values after a restart */
      if( scip->set->misc_calcintegral )
      {
         SCIPstatUpdatePrimalDualIntegrals(scip->stat, scip->set, scip->transprob, scip->origprob, SCIPgetUpperbound(scip), SCIPgetLowerbound(scip) );
      }
   }

   /* free temporary presolving root node */
   SCIP_CALL( SCIPtreeFreePresolvingRoot(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->messagehdlr,
         scip->stat, scip->transprob, scip->origprob, scip->primal, scip->lp, scip->branchcand, scip->conflict,
         scip->conflictstore, scip->eventfilter, scip->eventqueue, scip->cliquetable) );

   /* switch stage to PRESOLVED */
   scip->set->stage = SCIP_STAGE_PRESOLVED;

   return SCIP_OKAY;
}

/** applies one round of presolving with the given presolving timing
 *
 *  This method will always be called with presoltiming fast first. It iterates over all presolvers, propagators, and
 *  constraint handlers and calls their presolving callbacks with timing fast.  If enough reductions are found, it
 *  returns and the next presolving round will be started (again with timing fast).  If the fast presolving does not
 *  find enough reductions, this methods calls itself recursively with presoltiming medium.  Again, it calls the
 *  presolving callbacks of all presolvers, propagators, and constraint handlers with timing medium.  If enough
 *  reductions are found, it returns and the next presolving round will be started (with timing fast).  Otherwise, it is
 *  called recursively with presoltiming exhaustive. In exhaustive presolving, presolvers, propagators, and constraint
 *  handlers are called w.r.t. their priority, but this time, we stop as soon as enough reductions were found and do not
 *  necessarily call all presolving methods. If we stop, we return and another presolving round is started with timing
 *  fast.
 *
 *  @todo check if we want to do the following (currently disabled):
 *  In order to avoid calling the same expensive presolving methods again and again (which is possibly ineffective
 *  for the current instance), we continue the loop for exhaustive presolving where we stopped it the last time.  The
 *  {presol/prop/cons}start pointers are used to this end: they provide the plugins to start the loop with in the
 *  current presolving round (if we reach exhaustive presolving), and are updated in this case to the next ones to be
 *  called in the next round. In case we reach the end of the loop in exhaustive presolving, we call the method again
 *  with exhaustive timing, now starting with the first presolving steps in the loop until we reach the ones we started
 *  the last call with.  This way, we won't stop until all exhaustive presolvers were called without finding enough
 *  reductions (in sum).
 */
static
SCIP_RETCODE presolveRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLTIMING*    timing,             /**< pointer to current presolving timing */
   SCIP_Bool*            unbounded,          /**< pointer to store whether presolving detected unboundedness */
   SCIP_Bool*            infeasible,         /**< pointer to store whether presolving detected infeasibility */
   SCIP_Bool             lastround,          /**< is this the last presolving round due to a presolving round limit? */
   int*                  presolstart,        /**< pointer to get the presolver to start exhaustive presolving with in
                                              *   the current round and store the one to start with in the next round */
   int                   presolend,          /**< last presolver to treat in exhaustive presolving */
   int*                  propstart,          /**< pointer to get the propagator to start exhaustive presolving with in
                                              *   the current round and store the one to start with in the next round */
   int                   propend,            /**< last propagator to treat in exhaustive presolving */
   int*                  consstart,          /**< pointer to get the constraint handler to start exhaustive presolving with in
                                              *   the current round and store the one to start with in the next round */
   int                   consend             /**< last constraint handler to treat in exhaustive presolving */
   )
{
   SCIP_RESULT result;
   SCIP_EVENT event;
   SCIP_Bool aborted;
   SCIP_Bool lastranpresol;
#if 0
   int oldpresolstart = 0;
   int oldpropstart = 0;
   int oldconsstart = 0;
#endif
   int priopresol;
   int prioprop;
   int i;
   int j;
   int k;
#ifndef NDEBUG
   size_t nusedbuffers;
   size_t nusedcleanbuffers;
#endif

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(unbounded != NULL);
   assert(infeasible != NULL);
   assert(presolstart != NULL);
   assert(propstart != NULL);
   assert(consstart != NULL);

   assert((presolend == scip->set->npresols && propend == scip->set->nprops && consend == scip->set->nconshdlrs)
      || (*presolstart == 0 && *propstart == 0 && *consstart == 0));

   *unbounded = FALSE;
   *infeasible = FALSE;
   aborted = FALSE;

   assert( scip->set->propspresolsorted );

   /* GCG wants to perform presolving during the reading process of a file reader;
    * hence the number of used buffers does not need to be zero, however, it should not
    * change by calling the presolving callbacks
    */
#ifndef NDEBUG
   nusedbuffers = BMSgetNUsedBufferMemory(SCIPbuffer(scip));
   nusedcleanbuffers = BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip));
#endif

   if( *timing == SCIP_PRESOLTIMING_EXHAUSTIVE )
   {
      /* In exhaustive presolving, we continue the loop where we stopped last time to avoid calling the same
       * (possibly ineffective) presolving step again and again. If we reach the end of the arrays of presolvers,
       * propagators, and constraint handlers without having made enough reductions, we start again from the beginning
       */
      i = *presolstart;
      j = *propstart;
      k = *consstart;
#if 0
      oldpresolstart = i;
      oldpropstart = j;
      oldconsstart = k;
#endif
      if( i >= presolend && j >= propend && k >= consend )
         return SCIP_OKAY;

      if( i == 0 && j == 0 && k == 0 )
         ++(scip->stat->npresolroundsext);
   }
   else
   {
      /* in fast and medium presolving, we always iterate over all presolvers, propagators, and constraint handlers */
      assert(presolend == scip->set->npresols);
      assert(propend == scip->set->nprops);
      assert(consend == scip->set->nconshdlrs);

      i = 0;
      j = 0;
      k = 0;

      if( *timing == SCIP_PRESOLTIMING_FAST )
         ++(scip->stat->npresolroundsfast);
      if( *timing == SCIP_PRESOLTIMING_MEDIUM )
         ++(scip->stat->npresolroundsmed);
   }

   SCIPdebugMsg(scip, "starting presolving round %d (%d/%d/%d), timing = %u\n",
      scip->stat->npresolrounds, scip->stat->npresolroundsfast, scip->stat->npresolroundsmed,
      scip->stat->npresolroundsext, *timing);

   /* call included presolvers with nonnegative priority */
   while( !(*unbounded) && !(*infeasible) && !aborted && (i < presolend || j < propend) )
   {
      if( i < presolend )
         priopresol = SCIPpresolGetPriority(scip->set->presols[i]);
      else
         priopresol = -1;

      if( j < propend )
         prioprop = SCIPpropGetPresolPriority(scip->set->props_presol[j]);
      else
         prioprop = -1;

      /* call next propagator */
      if( prioprop >= priopresol )
      {
         /* only presolving methods which have non-negative priority will be called before constraint handlers */
         if( prioprop < 0 )
            break;

         SCIPdebugMsg(scip, "executing presolving of propagator <%s>\n", SCIPpropGetName(scip->set->props_presol[j]));
         SCIP_CALL( SCIPpropPresol(scip->set->props_presol[j], scip->set, *timing, scip->stat->npresolrounds,
               &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
               &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
               &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs,
               &scip->stat->npresolchgsides, &result) );
         assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
         assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

         lastranpresol = FALSE;
         ++j;
      }
      /* call next presolver */
      else
      {
         /* only presolving methods which have non-negative priority will be called before constraint handlers */
         if( priopresol < 0 )
            break;

         SCIPdebugMsg(scip, "executing presolver <%s>\n", SCIPpresolGetName(scip->set->presols[i]));
         SCIP_CALL( SCIPpresolExec(scip->set->presols[i], scip->set, *timing, scip->stat->npresolrounds,
               &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
               &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
               &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs,
               &scip->stat->npresolchgsides, &result) );
         assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
         assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

         lastranpresol = TRUE;
         ++i;
      }

      if( result == SCIP_CUTOFF )
      {
         *infeasible = TRUE;

         if( lastranpresol )
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "presolver <%s> detected infeasibility\n", SCIPpresolGetName(scip->set->presols[i-1]));
         else
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "propagator <%s> detected infeasibility\n", SCIPpropGetName(scip->set->props_presol[j-1]));
      }
      else if( result == SCIP_UNBOUNDED )
      {
         *unbounded = TRUE;

         if( lastranpresol )
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "presolver <%s> detected unboundedness (or infeasibility)\n", SCIPpresolGetName(scip->set->presols[i-1]));
         else
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "propagator <%s> detected  unboundedness (or infeasibility)\n", SCIPpropGetName(scip->set->props_presol[j-1]));
      }

      /* delete the variables from the problems that were marked to be deleted */
      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->cliquetable, scip->lp,
            scip->branchcand) );

      SCIPdebugMsg(scip, "presolving callback returned result <%d>\n", result);

      /* if we work off the exhaustive presolvers, we stop immediately if a reduction was found */
      if( (*timing == SCIP_PRESOLTIMING_EXHAUSTIVE) && !lastround && !SCIPisPresolveFinished(scip) )
      {
         assert(*consstart == 0);

         if( lastranpresol )
         {
            *presolstart = i + 1;
            *propstart = j;
         }
         else
         {
            *presolstart = i;
            *propstart = j + 1;
         }
         aborted = TRUE;

         break;
      }
   }

   /* call presolve methods of constraint handlers */
   while( k < consend && !(*unbounded) && !(*infeasible) && !aborted )
   {
      SCIPdebugMsg(scip, "executing presolve method of constraint handler <%s>\n",
         SCIPconshdlrGetName(scip->set->conshdlrs[k]));
      SCIP_CALL( SCIPconshdlrPresolve(scip->set->conshdlrs[k], scip->mem->probmem, scip->set, scip->stat,
            *timing, scip->stat->npresolrounds,
            &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
            &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
            &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs,
            &scip->stat->npresolchgsides, &result) );
      assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
      assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

      ++k;

      if( result == SCIP_CUTOFF )
      {
         *infeasible = TRUE;
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "constraint handler <%s> detected infeasibility\n", SCIPconshdlrGetName(scip->set->conshdlrs[k-1]));
      }
      else if( result == SCIP_UNBOUNDED )
      {
         *unbounded = TRUE;
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "constraint handler <%s> detected unboundedness (or infeasibility)\n",
            SCIPconshdlrGetName(scip->set->conshdlrs[k-1]));
      }

      /* delete the variables from the problems that were marked to be deleted */
      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->cliquetable, scip->lp,
            scip->branchcand) );

      SCIPdebugMsg(scip, "presolving callback returned with result <%d>\n", result);

      /* if we work off the exhaustive presolvers, we stop immediately if a reduction was found */
      if( (*timing == SCIP_PRESOLTIMING_EXHAUSTIVE) && !lastround && !SCIPisPresolveFinished(scip) )
      {
         *presolstart = i;
         *propstart = j;
         *consstart = k + 1;
         aborted = TRUE;

         break;
      }
   }

   assert( scip->set->propspresolsorted );

   /* call included presolvers with negative priority */
   while( !(*unbounded) && !(*infeasible) && !aborted && (i < presolend || j < propend) )
   {
      if( i < scip->set->npresols )
         priopresol = SCIPpresolGetPriority(scip->set->presols[i]);
      else
         priopresol = -INT_MAX;

      if( j < scip->set->nprops )
         prioprop = SCIPpropGetPresolPriority(scip->set->props_presol[j]);
      else
         prioprop = -INT_MAX;

      /* choose presolving */
      if( prioprop >= priopresol )
      {
         assert(prioprop <= 0);

         SCIPdebugMsg(scip, "executing presolving of propagator <%s>\n", SCIPpropGetName(scip->set->props_presol[j]));
         SCIP_CALL( SCIPpropPresol(scip->set->props_presol[j], scip->set, *timing, scip->stat->npresolrounds,
               &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
               &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
               &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs,
               &scip->stat->npresolchgsides, &result) );
         assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
         assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

         lastranpresol = FALSE;
         ++j;
      }
      else
      {
         assert(priopresol < 0);

         SCIPdebugMsg(scip, "executing presolver <%s>\n", SCIPpresolGetName(scip->set->presols[i]));
         SCIP_CALL( SCIPpresolExec(scip->set->presols[i], scip->set, *timing, scip->stat->npresolrounds,
               &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
               &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
               &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs,
               &scip->stat->npresolchgsides, &result) );
         assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
         assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

         lastranpresol = TRUE;
         ++i;
      }

      if( result == SCIP_CUTOFF )
      {
         *infeasible = TRUE;

         if( lastranpresol )
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "presolver <%s> detected infeasibility\n", SCIPpresolGetName(scip->set->presols[i-1]));
         else
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "propagator <%s> detected infeasibility\n", SCIPpropGetName(scip->set->props_presol[j-1]));
      }
      else if( result == SCIP_UNBOUNDED )
      {
         *unbounded = TRUE;

         if( lastranpresol )
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "presolver <%s> detected unboundedness (or infeasibility)\n", SCIPpresolGetName(scip->set->presols[i-1]));
         else
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "propagator <%s> detected  unboundedness (or infeasibility)\n", SCIPpropGetName(scip->set->props_presol[j-1]));
      }

      /* delete the variables from the problems that were marked to be deleted */
      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->cliquetable, scip->lp,
            scip->branchcand) );

      SCIPdebugMsg(scip, "presolving callback return with result <%d>\n", result);

      /* if we work off the exhaustive presolvers, we stop immediately if a reduction was found */
      if( (*timing == SCIP_PRESOLTIMING_EXHAUSTIVE) && !lastround && !SCIPisPresolveFinished(scip) )
      {
         assert(k == consend);

         if( lastranpresol )
         {
            *presolstart = i + 1;
            *propstart = j;
         }
         else
         {
            *presolstart = i;
            *propstart = j + 1;
         }
         *consstart = k;

         break;
      }
   }

   /* remove empty and single variable cliques from the clique table */
   if( !(*unbounded) && !(*infeasible) )
   {
      int nlocalbdchgs = 0;

      SCIP_CALL( SCIPcliquetableCleanup(scip->cliquetable, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, &nlocalbdchgs,
            infeasible) );

      if( nlocalbdchgs > 0 || *infeasible )
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "clique table cleanup detected %d bound changes%s\n", nlocalbdchgs, *infeasible ? " and infeasibility" : "");

      scip->stat->npresolfixedvars += nlocalbdchgs;

      if( !*infeasible && scip->set->nheurs > 0 )
      {
         /* call primal heuristics that are applicable during presolving */
         SCIP_Bool foundsol;

         SCIPdebugMsg(scip, "calling primal heuristics during presolving\n");

         /* call primal heuristics */
         SCIP_CALL( SCIPprimalHeuristics(scip->set, scip->stat, scip->transprob, scip->primal, NULL, NULL, NULL,
               SCIP_HEURTIMING_DURINGPRESOLLOOP, FALSE, &foundsol, unbounded) );

         /* output a message, if a solution was found */
         if( foundsol )
         {
            SCIP_SOL* sol;

            assert(SCIPgetNSols(scip) > 0);
            sol = SCIPgetBestSol(scip);
            assert(sol != NULL);
            assert(SCIPgetSolOrigObj(scip,sol) != SCIP_INVALID); /*lint !e777*/

            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
               "feasible solution found by %s heuristic after %.1f seconds, objective value %.6e\n",
               SCIPheurGetName(SCIPsolGetHeur(sol)), SCIPgetSolvingTime(scip), SCIPgetSolOrigObj(scip, sol));
         }
      }
   }

   if( !(*unbounded) && !(*infeasible) )
   {
      /* call more expensive presolvers */
      if( (SCIPisPresolveFinished(scip) || lastround) )
      {
         if( *timing != SCIP_PRESOLTIMING_FINAL )
         {
            assert((*timing == SCIP_PRESOLTIMING_FAST) || (*timing == SCIP_PRESOLTIMING_MEDIUM) || (*timing == SCIP_PRESOLTIMING_EXHAUSTIVE));

            SCIPdebugMsg(scip, "not enough reductions in %s presolving, running %s presolving now...\n",
               *timing == SCIP_PRESOLTIMING_FAST ? "fast" : *timing == SCIP_PRESOLTIMING_MEDIUM ? "medium" : "exhaustive",
               *timing == SCIP_PRESOLTIMING_FAST ? "medium" : *timing == SCIP_PRESOLTIMING_MEDIUM ? "exhaustive" : "final");

            /* increase timing */
            *timing = ((*timing == SCIP_PRESOLTIMING_FAST) ? SCIP_PRESOLTIMING_MEDIUM : (*timing == SCIP_PRESOLTIMING_MEDIUM) ? SCIP_PRESOLTIMING_EXHAUSTIVE : SCIP_PRESOLTIMING_FINAL);

            /* computational experiments showed that always starting the loop of exhaustive presolvers from the beginning
             * performs better than continuing from the last processed presolver. Therefore, we start from 0, but keep
             * the mechanisms to possibly change this back later.
             * @todo try starting from the last processed exhaustive presolver
             */
            *presolstart = 0;
            *propstart = 0;
            *consstart = 0;

            SCIP_CALL( presolveRound(scip, timing, unbounded, infeasible, lastround, presolstart, presolend,
                  propstart, propend, consstart, consend) );
         }
#if 0
         /* run remaining exhaustive presolvers (if we did not start from the beginning anyway) */
         else if( (oldpresolstart > 0 || oldpropstart > 0 || oldconsstart > 0) && presolend == scip->set->npresols
            && propend == scip->set->nprops && consend == scip->set->nconshdlrs )
         {
            int newpresolstart = 0;
            int newpropstart = 0;
            int newconsstart = 0;

            SCIPdebugMsg(scip, "reached end of exhaustive presolving loop, starting from the beginning...\n");

            SCIP_CALL( presolveRound(scip, timing, unbounded, infeasible, lastround, &newpresolstart,
                  oldpresolstart, &newpropstart, oldpropstart, &newconsstart, oldconsstart) );

            *presolstart = newpresolstart;
            *propstart = newpropstart;
            *consstart = newconsstart;
         }
#endif
      }
   }

   /* issue PRESOLVEROUND event */
   SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_PRESOLVEROUND) );
   SCIP_CALL( SCIPeventProcess(&event, scip->set, NULL, NULL, NULL, scip->eventfilter) );

   return SCIP_OKAY;
}


/** loops through the included presolvers and constraint's presolve methods, until changes are too few */
static
SCIP_RETCODE presolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            unbounded,          /**< pointer to store whether presolving detected unboundedness */
   SCIP_Bool*            infeasible,         /**< pointer to store whether presolving detected infeasibility */
   SCIP_Bool*            vanished            /**< pointer to store whether the problem vanished in presolving */
   )
{
   SCIP_PRESOLTIMING presoltiming;
   SCIP_Bool finished;
   SCIP_Bool stopped;
   SCIP_Bool lastround;
   int presolstart = 0;
   int propstart = 0;
   int consstart = 0;
#ifndef NDEBUG
   size_t nusedbuffers;
   size_t nusedcleanbuffers;
#endif

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->primal != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);
   assert(scip->set->stage == SCIP_STAGE_TRANSFORMED || scip->set->stage == SCIP_STAGE_PRESOLVING);
   assert(unbounded != NULL);
   assert(infeasible != NULL);

   *unbounded = FALSE;
   *vanished = FALSE;

   /* GCG wants to perform presolving during the reading process of a file reader;
    * hence the number of used buffers does not need to be zero, however, it should
    * be the same again after presolve is finished
    */
#ifndef NDEBUG
   nusedbuffers = BMSgetNUsedBufferMemory(SCIPbuffer(scip));
   nusedcleanbuffers = BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip));
#endif

   /* switch status to unknown */
   scip->stat->status = SCIP_STATUS_UNKNOWN;

   /* update upper bound and cutoff bound due to objective limit in primal data */
   SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
         scip->eventqueue, scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp) );

   /* start presolving timer */
   SCIPclockStart(scip->stat->presolvingtime, scip->set);
   SCIPclockStart(scip->stat->presolvingtimeoverall, scip->set);

   /* initialize presolving */
   if( scip->set->stage == SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( initPresolve(scip) );
   }
   assert(scip->set->stage == SCIP_STAGE_PRESOLVING);

   /* call primal heuristics that are applicable before presolving */
   if( scip->set->nheurs > 0 )
   {
      SCIP_Bool foundsol;

      SCIPdebugMsg(scip, "calling primal heuristics before presolving\n");

      /* call primal heuristics */
      SCIP_CALL( SCIPprimalHeuristics(scip->set, scip->stat, scip->transprob, scip->primal, NULL, NULL, NULL,
            SCIP_HEURTIMING_BEFOREPRESOL, FALSE, &foundsol, unbounded) );

      /* output a message, if a solution was found */
      if( foundsol )
      {
         SCIP_SOL* sol;

         assert(SCIPgetNSols(scip) > 0);
         sol = SCIPgetBestSol(scip);
         assert(sol != NULL);
         assert(SCIPgetSolOrigObj(scip,sol) != SCIP_INVALID);  /*lint !e777*/

         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "feasible solution found by %s heuristic after %.1f seconds, objective value %.6e\n",
            SCIPheurGetName(SCIPsolGetHeur(sol)), SCIPgetSolvingTime(scip), SCIPgetSolOrigObj(scip, sol));
      }
   }

   SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "presolving:\n");

   *infeasible = FALSE;
   *unbounded = (*unbounded) || (SCIPgetNSols(scip) > 0 && SCIPisInfinity(scip, -SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip))));

   finished = (scip->set->presol_maxrounds != -1 && scip->stat->npresolrounds >= scip->set->presol_maxrounds)
         || (*unbounded) || (scip->set->reopt_enable && scip->stat->nreoptruns >= 1);
   stopped = SCIPsolveIsStopped(scip->set, scip->stat, TRUE);

   /* perform presolving rounds */
   while( !finished && !stopped )
   {
      /* store current number of reductions */
      scip->stat->lastnpresolfixedvars = scip->stat->npresolfixedvars;
      scip->stat->lastnpresolaggrvars = scip->stat->npresolaggrvars;
      scip->stat->lastnpresolchgvartypes = scip->stat->npresolchgvartypes;
      scip->stat->lastnpresolchgbds = scip->stat->npresolchgbds;
      scip->stat->lastnpresoladdholes = scip->stat->npresoladdholes;
      scip->stat->lastnpresoldelconss = scip->stat->npresoldelconss;
      scip->stat->lastnpresoladdconss = scip->stat->npresoladdconss;
      scip->stat->lastnpresolupgdconss = scip->stat->npresolupgdconss;
      scip->stat->lastnpresolchgcoefs = scip->stat->npresolchgcoefs;
      scip->stat->lastnpresolchgsides = scip->stat->npresolchgsides;
#ifdef SCIP_DISABLED_CODE
      scip->stat->lastnpresolimplications = scip->stat->nimplications;
      scip->stat->lastnpresolcliques = SCIPcliquetableGetNCliques(scip->cliquetable);
#endif

      /* set presolving flag */
      scip->stat->performpresol = TRUE;

      /* sort propagators */
      SCIPsetSortPropsPresol(scip->set);

      /* sort presolvers by priority */
      SCIPsetSortPresols(scip->set);

      /* check if this will be the last presolving round (in that case, we want to run all presolvers) */
      lastround = (scip->set->presol_maxrounds == -1 ? FALSE : (scip->stat->npresolrounds + 1 >= scip->set->presol_maxrounds));

      presoltiming = SCIP_PRESOLTIMING_FAST;

      /* perform the presolving round by calling the presolvers, propagators, and constraint handlers */
      assert(!(*unbounded));
      assert(!(*infeasible));
      SCIP_CALL( presolveRound(scip, &presoltiming, unbounded, infeasible, lastround,
            &presolstart, scip->set->npresols, &propstart, scip->set->nprops, &consstart, scip->set->nconshdlrs) );

      /* check, if we should abort presolving due to not enough changes in the last round */
      finished = SCIPisPresolveFinished(scip) || presoltiming == SCIP_PRESOLTIMING_FINAL;

      SCIPdebugMsg(scip, "presolving round %d returned with unbounded = %u, infeasible = %u, finished = %u\n", scip->stat->npresolrounds, *unbounded, *infeasible, finished);

      /* check whether problem is infeasible or unbounded */
      finished = finished || *unbounded || *infeasible;

      /* increase round number */
      scip->stat->npresolrounds++;

      if( !finished )
      {
         /* print presolving statistics */
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "(round %d, %-11s %d del vars, %d del conss, %d add conss, %d chg bounds, %d chg sides, %d chg coeffs, %d upgd conss, %d impls, %d clqs\n",
            scip->stat->npresolrounds, ( presoltiming == SCIP_PRESOLTIMING_FAST ? "fast)" :
               (presoltiming == SCIP_PRESOLTIMING_MEDIUM ? "medium)" :
                  (presoltiming == SCIP_PRESOLTIMING_EXHAUSTIVE ?"exhaustive)" :
                     "final)")) ),
            scip->stat->npresolfixedvars + scip->stat->npresolaggrvars,
            scip->stat->npresoldelconss, scip->stat->npresoladdconss,
            scip->stat->npresolchgbds, scip->stat->npresolchgsides,
            scip->stat->npresolchgcoefs, scip->stat->npresolupgdconss,
            scip->stat->nimplications, SCIPcliquetableGetNCliques(scip->cliquetable));
      }

      /* abort if time limit was reached or user interrupted */
      stopped = SCIPsolveIsStopped(scip->set, scip->stat, TRUE);
   }

   /* first change status of scip, so that all plugins in their exitpre callbacks can ask SCIP for the correct status */
   if( *infeasible )
   {
      /* switch status to OPTIMAL */
      if( scip->primal->nlimsolsfound > 0 )
      {
         scip->stat->status = SCIP_STATUS_OPTIMAL;
      }
      else /* switch status to INFEASIBLE */
         scip->stat->status = SCIP_STATUS_INFEASIBLE;
   }
   else if( *unbounded )
   {
      if( scip->primal->nsols >= 1 ) /* switch status to UNBOUNDED */
         scip->stat->status = SCIP_STATUS_UNBOUNDED;
      else /* switch status to INFORUNBD */
         scip->stat->status = SCIP_STATUS_INFORUNBD;
   }
   /* if no variables and constraints are present, we try to add the empty solution (constraint handlers with needscons
    * flag FALSE could theoretically reject it); if no active pricers could create variables later, we conclude
    * optimality or infeasibility */
   else if( scip->transprob->nvars == 0 && scip->transprob->nconss == 0 )
   {
      SCIP_SOL* sol;
      SCIP_Bool stored;

      SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, FALSE, FALSE, FALSE, &stored) );

      if( scip->set->nactivepricers == 0 )
      {
         if( scip->primal->nlimsolsfound > 0 )
            scip->stat->status = SCIP_STATUS_OPTIMAL;
         else
            scip->stat->status = SCIP_STATUS_INFEASIBLE;

         *vanished = TRUE;
      }
   }

   /* deinitialize presolving */
   if( finished && (!stopped || *unbounded || *infeasible || *vanished) )
   {
      SCIP_Real maxnonzeros;
      SCIP_Longint nchecknonzeros;
      SCIP_Longint nactivenonzeros;
      SCIP_Bool approxchecknonzeros;
      SCIP_Bool approxactivenonzeros;
      SCIP_Bool infeas;

      SCIP_CALL( exitPresolve(scip, *unbounded || *infeasible || *vanished, &infeas) );
      *infeasible = *infeasible || infeas;

      assert(scip->set->stage == SCIP_STAGE_PRESOLVED);

      /* resort variables if we are not already done (unless variable permutation was explicitly activated) */
      if( !scip->set->random_permutevars && !(*infeasible) && !(*unbounded) && !(*vanished) )
      {
         /* (Re)Sort the variables, which appear in the four categories (binary, integer, implicit, continuous) after
          * presolve with respect to their original index (within their categories). Adjust the problem index afterwards
          * which is supposed to reflect the position in the variable array. This additional (re)sorting is supposed to
          * get more robust against the order presolving fixed variables. (We also reobtain a possible block structure
          * induced by the user model)
          */
         SCIPprobResortVars(scip->transprob);
      }

      /* determine number of non-zeros */
      maxnonzeros = (SCIP_Real)SCIPgetNConss(scip) * SCIPgetNVars(scip);
      maxnonzeros = MAX(maxnonzeros, 1.0);
      SCIP_CALL( calcNonZeros(scip, &nchecknonzeros, &nactivenonzeros, &approxchecknonzeros, &approxactivenonzeros) );
      scip->stat->nnz = nactivenonzeros;

      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL, "\n");
      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "presolved problem has %s%" SCIP_LONGINT_FORMAT " active (%g%%) nonzeros and %s%" SCIP_LONGINT_FORMAT " (%g%%) check nonzeros\n",
         approxactivenonzeros ? "more than " : "", nactivenonzeros, nactivenonzeros/maxnonzeros * 100,
         approxchecknonzeros ? "more than " : "", nchecknonzeros, nchecknonzeros/maxnonzeros * 100);
      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL, "\n");
   }
   assert(BMSgetNUsedBufferMemory(SCIPbuffer(scip)) == nusedbuffers);
   assert(BMSgetNUsedBufferMemory(SCIPcleanbuffer(scip)) == nusedcleanbuffers);

   /* stop presolving time */
   SCIPclockStop(scip->stat->presolvingtime, scip->set);
   SCIPclockStop(scip->stat->presolvingtimeoverall, scip->set);

   /* print presolving statistics */
   SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      "presolving (%d rounds: %d fast, %d medium, %d exhaustive):\n", scip->stat->npresolrounds,
      scip->stat->npresolroundsfast, scip->stat->npresolroundsmed, scip->stat->npresolroundsext);
   SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      " %d deleted vars, %d deleted constraints, %d added constraints, %d tightened bounds, %d added holes, %d changed sides, %d changed coefficients\n",
      scip->stat->npresolfixedvars + scip->stat->npresolaggrvars, scip->stat->npresoldelconss, scip->stat->npresoladdconss,
      scip->stat->npresolchgbds, scip->stat->npresoladdholes, scip->stat->npresolchgsides, scip->stat->npresolchgcoefs);
   SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      " %d implications, %d cliques\n", scip->stat->nimplications, SCIPcliquetableGetNCliques(scip->cliquetable));

   /* remember number of constraints */
   SCIPprobMarkNConss(scip->transprob);

   return SCIP_OKAY;
}

/** tries to transform original solutions to the transformed problem space */
static
SCIP_RETCODE transformSols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SOL** sols;
   SCIP_SOL** scipsols;
   SCIP_SOL* sol;
   SCIP_Real* solvals;
   SCIP_Bool* solvalset;
   SCIP_Bool added;
   SCIP_Longint oldnsolsfound;
   int nsols;
   int ntransvars;
   int naddedsols;
   int s;

   nsols = SCIPgetNSols(scip);
   oldnsolsfound = scip->primal->nsolsfound;

   /* no solution to transform */
   if( nsols == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "try to transfer %d original solutions into the transformed problem space\n", nsols);

   ntransvars = scip->transprob->nvars;
   naddedsols = 0;

   /* It might happen, that the added transferred solution does not equal the corresponding original one, which might
    * result in the array of solutions being changed.  Thus we temporarily copy the array and traverse it in reverse
    * order to ensure that the regarded solution in the copied array was not already freed when new solutions were added
    * and the worst solutions were freed.
    */
   scipsols = SCIPgetSols(scip);
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sols, scipsols, nsols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvalset, ntransvars) );

   for( s = nsols-1; s >= 0; --s )
   {
      sol = sols[s];

      /* it might happen that a transferred original solution has a better objective than its original counterpart
       * (e.g., because multi-aggregated variables get another value, but the solution is still feasible);
       * in this case, it might happen that the solution is not an original one and we just skip this solution
       */
      if( !SCIPsolIsOriginal(sol) )
         continue;

      SCIP_CALL( SCIPprimalTransformSol(scip->primal, sol, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
            scip->origprob, scip->transprob, scip->tree, scip->reopt, scip->lp, scip->eventqueue, scip->eventfilter, solvals,
            solvalset, ntransvars, &added) );

      if( added )
         ++naddedsols;
   }

   if( naddedsols > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "transformed %d/%d original solutions to the transformed problem space\n",
         naddedsols, nsols);

      scip->stat->nexternalsolsfound += scip->primal->nsolsfound - oldnsolsfound;
   }

   SCIPfreeBufferArray(scip, &solvalset);
   SCIPfreeBufferArray(scip, &solvals);
   SCIPfreeBufferArray(scip, &sols);

   return SCIP_OKAY;
}

/** initializes solution process data structures */
static
SCIP_RETCODE initSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             solved              /**< is problem already solved? */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->nlp == NULL);
   assert(scip->set->stage == SCIP_STAGE_PRESOLVED);

   /**@todo check whether other methodscan be skipped if problem has been solved */
   /* if problem has been solved, several time consuming tasks must not be performed */
   if( !solved )
   {
      /* reset statistics for current branch and bound run */
      SCIPstatResetCurrentRun(scip->stat, scip->set, scip->transprob, scip->origprob, solved);
      SCIPstatEnforceLPUpdates(scip->stat);

      /* LP is empty anyway; mark empty LP to be solved and update validsollp counter */
      SCIP_CALL( SCIPlpReset(scip->lp, scip->mem->probmem, scip->set, scip->transprob, scip->stat, scip->eventqueue, scip->eventfilter) );

      /* update upper bound and cutoff bound due to objective limit in primal data */
      SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
       scip->eventqueue, scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp) );
   }

   /* switch stage to INITSOLVE */
   scip->set->stage = SCIP_STAGE_INITSOLVE;

   /* initialize NLP if there are nonlinearities */
   if( scip->transprob->nlpenabled && !scip->set->nlp_disable )
   {
      SCIPdebugMsg(scip, "constructing empty NLP\n");

      SCIP_CALL( SCIPnlpCreate(&scip->nlp, scip->mem->probmem, scip->set, scip->stat, SCIPprobGetName(scip->transprob), scip->transprob->nvars) );
      assert(scip->nlp != NULL);

      SCIP_CALL( SCIPnlpAddVars(scip->nlp, scip->mem->probmem, scip->set, scip->transprob->nvars, scip->transprob->vars) );

      /* Adjust estimation of external memory: SCIPtransformProb() estimated the memory used for the LP-solver. As a
       * very crude approximation just double this number. Only do this once in the first run. */
      if( scip->set->misc_estimexternmem && scip->stat->nruns <= 1 )
      {
         scip->stat->externmemestim *= 2;
         SCIPdebugMsg(scip, "external memory usage estimated to %" SCIP_LONGINT_FORMAT " byte\n", scip->stat->externmemestim);
      }
   }

   /* possibly create visualization output file */
   SCIP_CALL( SCIPvisualInit(scip->stat->visual, scip->mem->probmem, scip->set, scip->messagehdlr) );

   /* initialize solution process data structures */
   SCIP_CALL( SCIPpricestoreCreate(&scip->pricestore) );
   SCIP_CALL( SCIPsepastoreCreate(&scip->sepastore, scip->mem->probmem, scip->set) );
   SCIP_CALL( SCIPsepastoreCreate(&scip->sepastoreprobing, scip->mem->probmem, scip->set) );
   SCIP_CALL( SCIPcutpoolCreate(&scip->cutpool, scip->mem->probmem, scip->set, scip->set->sepa_cutagelimit, TRUE) );
   SCIP_CALL( SCIPcutpoolCreate(&scip->delayedcutpool, scip->mem->probmem, scip->set, scip->set->sepa_cutagelimit, FALSE) );
   SCIP_CALL( SCIPtreeCreateRoot(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter, scip->eventqueue,
         scip->lp) );

   /* update dual bound of the root node if a valid dual bound is at hand */
   if( scip->transprob->dualbound < SCIP_INVALID )
   {
      SCIP_Real internobjval = SCIPprobInternObjval(scip->transprob, scip->origprob, scip->set, scip->transprob->dualbound);

      scip->stat->lastlowerbound = internobjval;

      SCIPnodeUpdateLowerbound(SCIPtreeGetRootNode(scip->tree), scip->stat, scip->set, scip->tree, scip->transprob,
         scip->origprob, internobjval);
   }

   /* try to transform original solutions to the transformed problem space */
   if( scip->set->misc_transorigsols )
   {
      SCIP_CALL( transformSols(scip) );
   }

   /* inform the transformed problem that the branch and bound process starts now */
   SCIP_CALL( SCIPprobInitSolve(scip->transprob, scip->set) );

   /* transform the decomposition storage */
   SCIP_CALL( SCIPtransformDecompstore(scip) );

   /* inform plugins that the branch and bound process starts now */
   SCIP_CALL( SCIPsetInitsolPlugins(scip->set, scip->mem->probmem, scip->stat) );

   /* remember number of constraints */
   SCIPprobMarkNConss(scip->transprob);

   /* if all variables are known, calculate a trivial primal bound by setting all variables to their worst bound */
   if( scip->set->nactivepricers == 0 )
   {
      SCIP_VAR* var;
      SCIP_Real obj;
      SCIP_Real objbound;
      SCIP_Real bd;
      int v;

      objbound = 0.0;
      for( v = 0; v < scip->transprob->nvars && !SCIPsetIsInfinity(scip->set, objbound); ++v )
      {
         var = scip->transprob->vars[v];
         obj = SCIPvarGetObj(var);
         if( !SCIPsetIsZero(scip->set, obj) )
         {
            bd = SCIPvarGetWorstBoundGlobal(var);
            if( SCIPsetIsInfinity(scip->set, REALABS(bd)) )
               objbound = SCIPsetInfinity(scip->set);
            else
               objbound += obj * bd;
         }
      }

      /* adjust primal bound, such that solution with worst bound may be found */
      if( objbound + SCIPsetCutoffbounddelta(scip->set) != objbound ) /*lint !e777*/
         objbound += SCIPsetCutoffbounddelta(scip->set);
      /* if objbound is very large, adding the cutoffbounddelta may not change the number; in this case, we are using
       * SCIPnextafter to ensure that the cutoffbound is really larger than the best possible solution value
       */
      else
         objbound = SCIPnextafter(objbound, SCIP_REAL_MAX);

      /* update cutoff bound */
      if( !SCIPsetIsInfinity(scip->set, objbound) && SCIPsetIsLT(scip->set, objbound, scip->primal->cutoffbound) )
      {
         /* adjust cutoff bound */
         SCIP_CALL( SCIPprimalSetCutoffbound(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
               scip->eventqueue, scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, objbound, FALSE) );
      }
   }

   /* switch stage to SOLVING */
   scip->set->stage = SCIP_STAGE_SOLVING;

   return SCIP_OKAY;
}

/** frees solution process data structures */
static
SCIP_RETCODE freeSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             restart             /**< was this free solve call triggered by a restart? */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->set->stage == SCIP_STAGE_SOLVING || scip->set->stage == SCIP_STAGE_SOLVED);

   /* mark that we are currently restarting */
   if( restart )
   {
      scip->stat->inrestart = TRUE;

      /* copy the current dual bound into the problem data structure such that it can be used initialize the new search
       * tree
       */
      SCIPprobUpdateDualbound(scip->transprob, SCIPgetDualbound(scip));
   }

   /* remove focus from the current focus node */
   if( SCIPtreeGetFocusNode(scip->tree) != NULL )
   {
      SCIP_NODE* node = NULL;
      SCIP_Bool cutoff;

      SCIP_CALL( SCIPnodeFocus(&node, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->transprob,
            scip->origprob, scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->conflict,
            scip->conflictstore, scip->eventfilter, scip->eventqueue, scip->cliquetable, &cutoff, FALSE, TRUE) );
      assert(!cutoff);
   }

   /* switch stage to EXITSOLVE */
   scip->set->stage = SCIP_STAGE_EXITSOLVE;

   /* cleanup the conflict storage */
   SCIP_CALL( SCIPconflictstoreClean(scip->conflictstore, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->reopt) );

   /* inform plugins that the branch and bound process is finished */
   SCIP_CALL( SCIPsetExitsolPlugins(scip->set, scip->mem->probmem, scip->stat, restart) );

   /* free the NLP, if there is one, and reset the flags indicating nonlinearity */
   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpFree(&scip->nlp, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
   }
   scip->transprob->nlpenabled = FALSE;

   /* clear the LP, and flush the changes to clear the LP of the solver */
   SCIP_CALL( SCIPlpReset(scip->lp, scip->mem->probmem, scip->set, scip->transprob, scip->stat, scip->eventqueue, scip->eventfilter) );
   SCIPlpInvalidateRootObjval(scip->lp);

   /* resets the debug environment */
   SCIP_CALL( SCIPdebugReset(scip->set) ); /*lint !e506 !e774*/

   /* clear all row references in internal data structures */
   SCIP_CALL( SCIPcutpoolClear(scip->cutpool, scip->mem->probmem, scip->set, scip->lp) );
   SCIP_CALL( SCIPcutpoolClear(scip->delayedcutpool, scip->mem->probmem, scip->set, scip->lp) );

   /* we have to clear the tree prior to the problem deinitialization, because the rows stored in the forks and
    * subroots have to be released
    */
   SCIP_CALL( SCIPtreeClear(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter, scip->eventqueue, scip->lp) );

   SCIPexitSolveDecompstore(scip);

   /* deinitialize transformed problem */
   SCIP_CALL( SCIPprobExitSolve(scip->transprob, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, restart) );

   /* free solution process data structures */
   SCIP_CALL( SCIPcutpoolFree(&scip->cutpool, scip->mem->probmem, scip->set, scip->lp) );
   SCIP_CALL( SCIPcutpoolFree(&scip->delayedcutpool, scip->mem->probmem, scip->set, scip->lp) );
   SCIP_CALL( SCIPsepastoreFree(&scip->sepastoreprobing, scip->mem->probmem) );
   SCIP_CALL( SCIPsepastoreFree(&scip->sepastore, scip->mem->probmem) );
   SCIP_CALL( SCIPpricestoreFree(&scip->pricestore) );

   /* possibly close visualization output file */
   SCIPvisualExit(scip->stat->visual, scip->set, scip->messagehdlr);

   /* reset statistics for current branch and bound run */
   if( scip->stat->status == SCIP_STATUS_INFEASIBLE || scip->stat->status == SCIP_STATUS_OPTIMAL || scip->stat->status == SCIP_STATUS_UNBOUNDED || scip->stat->status == SCIP_STATUS_INFORUNBD )
      SCIPstatResetCurrentRun(scip->stat, scip->set, scip->transprob, scip->origprob, TRUE);
   else
      SCIPstatResetCurrentRun(scip->stat, scip->set, scip->transprob, scip->origprob, FALSE);

   /* switch stage to TRANSFORMED */
   scip->set->stage = SCIP_STAGE_TRANSFORMED;

   /* restart finished */
   assert( ! restart || scip->stat->inrestart );
   scip->stat->inrestart = FALSE;

   return SCIP_OKAY;
}

/** frees solution process data structures when reoptimization is used
 *
 *  in contrast to a freeSolve() this method will preserve the transformed problem such that another presolving round
 *  after changing the problem (modifying the objective function) is not necessary.
 */
static
SCIP_RETCODE freeReoptSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->set->stage == SCIP_STAGE_SOLVING || scip->set->stage == SCIP_STAGE_SOLVED);

   /* remove focus from the current focus node */
   if( SCIPtreeGetFocusNode(scip->tree) != NULL )
   {
      SCIP_NODE* node = NULL;
      SCIP_Bool cutoff;

      SCIP_CALL( SCIPnodeFocus(&node, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->transprob,
            scip->origprob, scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->conflict,
            scip->conflictstore, scip->eventfilter, scip->eventqueue, scip->cliquetable, &cutoff, FALSE, TRUE) );
      assert(!cutoff);
   }

   /* mark current stats, such that new solve begins with the var/col/row indices from the previous run */
   SCIPstatMark(scip->stat);

   /* switch stage to EXITSOLVE */
   scip->set->stage = SCIP_STAGE_EXITSOLVE;

   /* deinitialize conflict store */
   SCIP_CALL( SCIPconflictstoreClear(scip->conflictstore, scip->mem->probmem, scip->set, scip->stat, scip->reopt) );

   /* invalidate the dual bound */
   SCIPprobInvalidateDualbound(scip->transprob);

   /* inform plugins that the branch and bound process is finished */
   SCIP_CALL( SCIPsetExitsolPlugins(scip->set, scip->mem->probmem, scip->stat, FALSE) );

   /* call exit methods of plugins */
   SCIP_CALL( SCIPsetExitPlugins(scip->set, scip->mem->probmem, scip->stat) );

   /* free the NLP, if there is one, and reset the flags indicating nonlinearity */
   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpFree(&scip->nlp, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
   }
   scip->transprob->nlpenabled = FALSE;

   /* clear the LP, and flush the changes to clear the LP of the solver */
   SCIP_CALL( SCIPlpReset(scip->lp, scip->mem->probmem, scip->set, scip->transprob, scip->stat, scip->eventqueue, scip->eventfilter) );
   SCIPlpInvalidateRootObjval(scip->lp);

   /* resets the debug environment */
   SCIP_CALL( SCIPdebugReset(scip->set) ); /*lint !e506 !e774*/

   /* clear all row references in internal data structures */
   SCIP_CALL( SCIPcutpoolClear(scip->cutpool, scip->mem->probmem, scip->set, scip->lp) );
   SCIP_CALL( SCIPcutpoolClear(scip->delayedcutpool, scip->mem->probmem, scip->set, scip->lp) );

   /* we have to clear the tree prior to the problem deinitialization, because the rows stored in the forks and
    * subroots have to be released
    */
   SCIP_CALL( SCIPtreeClear(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter, scip->eventqueue, scip->lp) );

   /* deinitialize transformed problem */
   SCIP_CALL( SCIPprobExitSolve(scip->transprob, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, FALSE) );

   /* free solution process data structures */
   SCIP_CALL( SCIPrelaxationFree(&scip->relaxation) );

   SCIP_CALL( SCIPcutpoolFree(&scip->cutpool, scip->mem->probmem, scip->set, scip->lp) );
   SCIP_CALL( SCIPcutpoolFree(&scip->delayedcutpool, scip->mem->probmem, scip->set, scip->lp) );
   SCIP_CALL( SCIPsepastoreFree(&scip->sepastoreprobing, scip->mem->probmem) );
   SCIP_CALL( SCIPsepastoreFree(&scip->sepastore, scip->mem->probmem) );
   SCIP_CALL( SCIPpricestoreFree(&scip->pricestore) );

   /* possibly close visualization output file */
   SCIPvisualExit(scip->stat->visual, scip->set, scip->messagehdlr);

   /* reset statistics for current branch and bound run */
   SCIPstatResetCurrentRun(scip->stat, scip->set, scip->transprob, scip->origprob, FALSE);

   /* switch stage to PRESOLVED */
   scip->set->stage = SCIP_STAGE_PRESOLVED;

   /* restart finished */
   scip->stat->inrestart = FALSE;

   /* reset solving specific paramters */
   if( scip->set->reopt_enable )
   {
      assert(scip->reopt != NULL);
      SCIP_CALL( SCIPreoptReset(scip->reopt, scip->set, scip->mem->probmem) );
   }

   /* free the debug solution which might live in transformed primal data structure */
   SCIP_CALL( SCIPprimalClear(&scip->primal, scip->mem->probmem) );

   if( scip->set->misc_resetstat )
   {
      /* reset statistics to the point before the problem was transformed */
      SCIPstatReset(scip->stat, scip->set, scip->transprob, scip->origprob);
   }
   else
   {
      /* even if statistics are not completely reset, a partial reset of the primal-dual integral is necessary */
      SCIPstatResetPrimalDualIntegrals(scip->stat, scip->set, TRUE);
   }

   /* reset objective limit */
   SCIP_CALL( SCIPsetObjlimit(scip, SCIP_INVALID) );

   return SCIP_OKAY;
}

/** free transformed problem */
static
SCIP_RETCODE freeTransform(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool reducedfree;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->stat != NULL);
   assert(scip->set->stage == SCIP_STAGE_TRANSFORMED || scip->set->stage == SCIP_STAGE_PRESOLVING ||
         (scip->set->stage == SCIP_STAGE_PRESOLVED && scip->set->reopt_enable));

   /* If the following evaluates to true, SCIPfreeReoptSolve() has already called the exit-callbacks of the plugins.
    * We can skip calling some of the following methods. This can happen if a new objective function was
    * installed but the solve was not started.
    */
   reducedfree = (scip->set->stage == SCIP_STAGE_PRESOLVED && scip->set->reopt_enable);

   if( !reducedfree )
   {
      /* call exit methods of plugins */
      SCIP_CALL( SCIPsetExitPlugins(scip->set, scip->mem->probmem, scip->stat) );
   }

   /* copy best primal solutions to original solution candidate list */
   if( !scip->set->reopt_enable && scip->set->limit_maxorigsol > 0 && scip->set->misc_transsolsorig && scip->set->nactivebenders == 0 )
   {
      SCIP_Bool stored;
      SCIP_Bool hasinfval;
      int maxsols;
      int nsols;
      int s;

      assert(scip->origprimal->nsols == 0);

      nsols = scip->primal->nsols;
      maxsols = scip->set->limit_maxorigsol;
      stored = TRUE;
      s = 0;

      /* iterate over all solutions as long as the original solution candidate store size limit is not reached */
      while( s < nsols && scip->origprimal->nsols < maxsols )
      {
         SCIP_SOL* sol;

         sol = scip->primal->sols[s];
         assert(sol != NULL);

         if( !SCIPsolIsOriginal(sol) )
         {
            /* retransform solution into the original problem space */
            SCIP_CALL( SCIPsolRetransform(sol, scip->set, scip->stat, scip->origprob, scip->transprob, &hasinfval) );
         }
         else
            hasinfval = FALSE;

         /* removing infinite fixings is turned off by the corresponding parameter */
         if( !scip->set->misc_finitesolstore )
            hasinfval = FALSE;

         if( !hasinfval )
         {
            /* add solution to original candidate solution storage */
            SCIP_CALL( SCIPprimalAddOrigSol(scip->origprimal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, sol, &stored) );
         }
         else
         {
            SCIP_SOL* newsol;
            SCIP_Bool success;

            SCIP_CALL( SCIPcreateFiniteSolCopy(scip, &newsol, sol, &success) );

            /* infinite fixing could be removed */
            if( newsol != NULL )
            {
               /* add solution to original candidate solution storage; we must not use SCIPprimalAddOrigSolFree()
                * because we want to create a copy of the solution in the origprimal solution store, but newsol was
                * created in the (transformed) primal
                */
               SCIP_CALL( SCIPprimalAddOrigSol(scip->origprimal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, newsol, &stored) );

               /* free solution in (transformed) primal where it was created */
               SCIP_CALL( SCIPsolFree(&newsol, scip->mem->probmem, scip->primal) );
            }
         }
         ++s;
      }

      if( scip->origprimal->nsols > 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "stored the %d best primal solutions in the original solution candidate list\n", scip->origprimal->nsols);
      }
      else if( scip->origprimal->nsols == 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "stored the best primal solution in the original solution candidate list\n");
      }
   }

   /* switch stage to FREETRANS */
   scip->set->stage = SCIP_STAGE_FREETRANS;

   /* reset solving specific paramters */
   assert(!scip->set->reopt_enable || scip->reopt != NULL);
   if( scip->set->reopt_enable && scip->reopt != NULL )
   {
      SCIP_CALL( SCIPreoptReset(scip->reopt, scip->set, scip->mem->probmem) );
   }

   if( !reducedfree )
   {
      /* clear the conflict store
       *
       * since the conflict store can contain transformed constraints we need to remove them. the store will be finally
       * freed in SCIPfreeProb().
       */
      SCIP_CALL( SCIPconflictstoreClear(scip->conflictstore, scip->mem->probmem, scip->set, scip->stat, scip->reopt) );
   }

   /* free transformed problem data structures */
   SCIP_CALL( SCIPprobFree(&scip->transprob, scip->messagehdlr, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
   SCIP_CALL( SCIPcliquetableFree(&scip->cliquetable, scip->mem->probmem) );
   SCIP_CALL( SCIPconflictFree(&scip->conflict, scip->mem->probmem) );

   if( !reducedfree )
   {
      SCIP_CALL( SCIPrelaxationFree(&scip->relaxation) );
   }
   SCIP_CALL( SCIPtreeFree(&scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter, scip->eventqueue, scip->lp) );

   /* free the debug solution which might live in transformed primal data structure */
   SCIP_CALL( SCIPdebugFreeSol(scip->set) ); /*lint !e506 !e774*/
   SCIP_CALL( SCIPprimalFree(&scip->primal, scip->mem->probmem) );

   SCIP_CALL( SCIPlpFree(&scip->lp, scip->mem->probmem, scip->set, scip->eventqueue, scip->eventfilter) );
   SCIP_CALL( SCIPbranchcandFree(&scip->branchcand) );
   SCIP_CALL( SCIPeventfilterFree(&scip->eventfilter, scip->mem->probmem, scip->set) );
   SCIP_CALL( SCIPeventqueueFree(&scip->eventqueue) );

   if( scip->set->misc_resetstat && !reducedfree )
   {
      /* reset statistics to the point before the problem was transformed */
      SCIPstatReset(scip->stat, scip->set, scip->transprob, scip->origprob);
   }
   else
   {
      /* even if statistics are not completely reset, a partial reset of the primal-dual integral is necessary */
      SCIPstatResetPrimalDualIntegrals(scip->stat, scip->set, TRUE);
   }

   /* switch stage to PROBLEM */
   scip->set->stage = SCIP_STAGE_PROBLEM;

   /* reset objective limit */
   SCIP_CALL( SCIPsetObjlimit(scip, SCIP_INVALID) );

   /* reset original variable's local and global bounds to their original values */
   SCIP_CALL( SCIPprobResetBounds(scip->origprob, scip->mem->probmem, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** free transformed problem in case an error occurs during transformation and return to SCIP_STAGE_PROBLEM */
static
SCIP_RETCODE freeTransforming(
   SCIP*                 scip                /**< SCIP data structure */
   )
{

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->stat != NULL);
   assert(scip->set->stage == SCIP_STAGE_TRANSFORMING);

   /* switch stage to FREETRANS */
   scip->set->stage = SCIP_STAGE_FREETRANS;

   /* free transformed problem data structures */
   SCIP_CALL( SCIPprobFree(&scip->transprob, scip->messagehdlr, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
   SCIP_CALL( SCIPcliquetableFree(&scip->cliquetable, scip->mem->probmem) );
   SCIP_CALL( SCIPconflictFree(&scip->conflict, scip->mem->probmem) );
   SCIP_CALL( SCIPrelaxationFree(&scip->relaxation) );
   SCIP_CALL( SCIPtreeFree(&scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter, scip->eventqueue, scip->lp) );

   /* free the debug solution which might live in transformed primal data structure */
   SCIP_CALL( SCIPdebugFreeSol(scip->set) ); /*lint !e506 !e774*/
   SCIP_CALL( SCIPprimalFree(&scip->primal, scip->mem->probmem) );

   SCIP_CALL( SCIPlpFree(&scip->lp, scip->mem->probmem, scip->set, scip->eventqueue, scip->eventfilter) );
   SCIP_CALL( SCIPbranchcandFree(&scip->branchcand) );
   SCIP_CALL( SCIPeventfilterFree(&scip->eventfilter, scip->mem->probmem, scip->set) );
   SCIP_CALL( SCIPeventqueueFree(&scip->eventqueue) );

   if( scip->set->misc_resetstat )
   {
      /* reset statistics to the point before the problem was transformed */
      SCIPstatReset(scip->stat, scip->set, scip->transprob, scip->origprob);
   }
   else
   {
      /* even if statistics are not completely reset, a partial reset of the primal-dual integral is necessary */
      SCIPstatResetPrimalDualIntegrals(scip->stat, scip->set, TRUE);
   }

   /* switch stage to PROBLEM */
   scip->set->stage = SCIP_STAGE_PROBLEM;

   return SCIP_OKAY;
}

/** displays most relevant statistics after problem was solved */
static
SCIP_RETCODE displayRelevantStats(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   /* display most relevant statistics */
   if( scip->set->disp_verblevel >= SCIP_VERBLEVEL_NORMAL && scip->set->disp_relevantstats )
   {
      SCIP_Bool objlimitreached = FALSE;

      /* We output that the objective limit has been reached if the problem has been solved, no solution respecting the
       * objective limit has been found (nlimsolsfound == 0) and the primal bound is finite. Note that it still might be
       * that the original problem is infeasible, even without the objective limit, i.e., we cannot be sure that we
       * actually reached the objective limit. */
      if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED && scip->primal->nlimsolsfound == 0 && ! SCIPisInfinity(scip, SCIPgetPrimalbound(scip)) )
         objlimitreached = TRUE;

      SCIPmessagePrintInfo(scip->messagehdlr, "\n");
      SCIPmessagePrintInfo(scip->messagehdlr, "SCIP Status        : ");
      SCIP_CALL( SCIPprintStage(scip, NULL) );
      SCIPmessagePrintInfo(scip->messagehdlr, "\n");
      if( scip->set->reopt_enable )
         SCIPmessagePrintInfo(scip->messagehdlr, "Solving Time (sec) : %.2f (over %d runs: %.2f)\n", SCIPclockGetTime(scip->stat->solvingtime), scip->stat->nreoptruns, SCIPclockGetTime(scip->stat->solvingtimeoverall));
      else
         SCIPmessagePrintInfo(scip->messagehdlr, "Solving Time (sec) : %.2f\n", SCIPclockGetTime(scip->stat->solvingtime));
      if( scip->stat->nruns > 1 )
         SCIPmessagePrintInfo(scip->messagehdlr, "Solving Nodes      : %" SCIP_LONGINT_FORMAT " (total of %" SCIP_LONGINT_FORMAT " nodes in %d runs)\n",
            scip->stat->nnodes, scip->stat->ntotalnodes, scip->stat->nruns);
      else if( scip->set->reopt_enable )
      {
         SCIP_BRANCHRULE* branchrule;

         branchrule = SCIPfindBranchrule(scip, "nodereopt");
         assert(branchrule != NULL);

         SCIPmessagePrintInfo(scip->messagehdlr, "Solving Nodes      : %" SCIP_LONGINT_FORMAT " (%" SCIP_LONGINT_FORMAT " reactivated)\n", scip->stat->nnodes, SCIPbranchruleGetNChildren(branchrule));
      }
      else
         SCIPmessagePrintInfo(scip->messagehdlr, "Solving Nodes      : %" SCIP_LONGINT_FORMAT "\n", scip->stat->nnodes);
      if( scip->set->stage >= SCIP_STAGE_TRANSFORMED && scip->set->stage <= SCIP_STAGE_EXITSOLVE )
      {
         if( objlimitreached )
         {
            SCIPmessagePrintInfo(scip->messagehdlr, "Primal Bound       : %+.14e (objective limit, %" SCIP_LONGINT_FORMAT " solutions",
               SCIPgetPrimalbound(scip), scip->primal->nsolsfound);
            if( scip->primal->nsolsfound > 0 )
            {
               SCIPmessagePrintInfo(scip->messagehdlr, ", best solution %+.14e", SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip)));
            }
            SCIPmessagePrintInfo(scip->messagehdlr, ")\n");
         }
         else
         {
            char limsolstring[SCIP_MAXSTRLEN];
            if( scip->primal->nsolsfound != scip->primal->nlimsolsfound )
               (void) SCIPsnprintf(limsolstring, SCIP_MAXSTRLEN, ", %" SCIP_LONGINT_FORMAT " respecting the objective limit", scip->primal->nlimsolsfound);
            else
               (void) SCIPsnprintf(limsolstring, SCIP_MAXSTRLEN,"");

            SCIPmessagePrintInfo(scip->messagehdlr, "Primal Bound       : %+.14e (%" SCIP_LONGINT_FORMAT " solutions%s)\n",
               SCIPgetPrimalbound(scip), scip->primal->nsolsfound, limsolstring);
         }
      }
      if( scip->set->stage >= SCIP_STAGE_SOLVING && scip->set->stage <= SCIP_STAGE_SOLVED )
      {
         SCIPmessagePrintInfo(scip->messagehdlr, "Dual Bound         : %+.14e\n", SCIPgetDualbound(scip));

         SCIPmessagePrintInfo(scip->messagehdlr, "Gap                : ");
         if( SCIPsetIsInfinity(scip->set, SCIPgetGap(scip)) )
            SCIPmessagePrintInfo(scip->messagehdlr, "infinite\n");
         else
            SCIPmessagePrintInfo(scip->messagehdlr, "%.2f %%\n", 100.0*SCIPgetGap(scip));
      }

      /* check solution for feasibility in original problem */
      if( scip->set->stage >= SCIP_STAGE_TRANSFORMED )
      {
         SCIP_SOL* sol;

         sol = SCIPgetBestSol(scip);
         if( sol != NULL )
         {
            SCIP_Real checkfeastolfac;
            SCIP_Real oldfeastol;
            SCIP_Bool dispallviols;
            SCIP_Bool feasible;

            oldfeastol = SCIPfeastol(scip);
            SCIP_CALL( SCIPgetRealParam(scip, "numerics/checkfeastolfac", &checkfeastolfac) );
            SCIP_CALL( SCIPgetBoolParam(scip, "display/allviols", &dispallviols) );

            /* scale feasibility tolerance by set->num_checkfeastolfac */
            if( !SCIPisEQ(scip, checkfeastolfac, 1.0) )
            {
               SCIP_CALL( SCIPchgFeastol(scip, oldfeastol * checkfeastolfac) );
            }

            SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, dispallviols) );

            /* restore old feasibilty tolerance */
            if( !SCIPisEQ(scip, checkfeastolfac, 1.0) )
            {
               SCIP_CALL( SCIPchgFeastol(scip, oldfeastol) );
            }

            if( !feasible )
            {
               SCIPmessagePrintInfo(scip->messagehdlr, "best solution is not feasible in original problem\n");
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** calls compression based on the reoptimization structure after the presolving */
static
SCIP_RETCODE compressReoptTree(
   SCIP*                 scip                /**< global SCIP settings */
   )
{
   SCIP_RESULT result;
   int c;
   int noldnodes;
   int nnewnodes;

   result = SCIP_DIDNOTFIND;

   noldnodes = SCIPreoptGetNNodes(scip->reopt, scip->tree->root);

   /* do not run if there exists only the root node */
   if( noldnodes <= 1 )
      return SCIP_OKAY;

   /* do not run a tree compression if the problem contains (implicit) integer variables */
   if( scip->transprob->nintvars > 0 || scip->transprob->nimplvars > 0 )
      return SCIP_OKAY;

   SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "tree compression:\n");
   SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "  given tree has %d nodes.\n", noldnodes);

   /* sort compressions by priority */
   SCIPsetSortComprs(scip->set);

   for(c = 0; c < scip->set->ncomprs; c++)
   {
      assert(result == SCIP_DIDNOTFIND || result == SCIP_DIDNOTRUN);

      /* call tree compression technique */
      SCIP_CALL( SCIPcomprExec(scip->set->comprs[c], scip->set, scip->reopt, &result) );

      if( result == SCIP_SUCCESS )
      {
         nnewnodes = SCIPreoptGetNNodes(scip->reopt, scip->tree->root);
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
               "  <%s> compressed the search tree to %d nodes (rate %g).\n", SCIPcomprGetName(scip->set->comprs[c]),
               nnewnodes, ((SCIP_Real)nnewnodes)/noldnodes);

         break;
      }
   }

   if( result != SCIP_SUCCESS )
   {
      assert(result == SCIP_DIDNOTFIND || result == SCIP_DIDNOTRUN);
      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "  search tree could not be compressed.\n");
   }

   return SCIP_OKAY;
}

/* prepare all plugins and data structures for a reoptimization run */
static
SCIP_RETCODE prepareReoptimization(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool reoptrestart;

   assert(scip != NULL);
   assert(scip->set->reopt_enable);

   /* @ todo: we could check if the problem is feasible, eg, by backtracking */

   /* increase number of reopt_runs */
   ++scip->stat->nreoptruns;

   /* inform the reoptimization plugin that a new iteration starts */
   SCIP_CALL( SCIPreoptAddRun(scip->reopt, scip->set, scip->mem->probmem, scip->origprob->vars,
         scip->origprob->nvars, scip->set->limit_maxsol) );

   /* check whether we need to add globally valid constraints */
   if( scip->set->reopt_sepaglbinfsubtrees || scip->set->reopt_sepabestsol )
   {
      SCIP_CALL( SCIPreoptApplyGlbConss(scip, scip->reopt, scip->set, scip->stat, scip->mem->probmem) );
   }

   /* after presolving the problem the first time we remember all global bounds and active constraints. bounds and
    * constraints will be restored within SCIPreoptInstallBounds() and SCIPreoptResetActiveConss().
    */
   if( scip->stat->nreoptruns == 1 )
   {
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED || scip->set->stage == SCIP_STAGE_SOLVED);

      SCIP_CALL( SCIPreoptSaveGlobalBounds(scip->reopt, scip->transprob, scip->mem->probmem) );

      SCIP_CALL( SCIPreoptSaveActiveConss(scip->reopt, scip->set, scip->transprob, scip->mem->probmem) );
   }
   /* we are at least in the second run */
   else
   {
      assert(scip->transprob != NULL);

      SCIP_CALL( SCIPreoptMergeVarHistory(scip->reopt, scip->set, scip->stat, scip->origprob->vars, scip->origprob->nvars) );

      SCIP_CALL( SCIPrelaxationCreate(&scip->relaxation, scip->mem->probmem, scip->set, scip->stat, scip->primal,
            scip->tree) );

      /* mark statistics before solving */
      SCIPstatMark(scip->stat);

      SCIPbranchcandInvalidate(scip->branchcand);

      SCIP_CALL( SCIPreoptResetActiveConss(scip->reopt, scip->set, scip->stat) );

      /* check whether we want to restart the tree search */
      SCIP_CALL( SCIPreoptCheckRestart(scip->reopt, scip->set, scip->mem->probmem, NULL, scip->transprob->vars,
            scip->transprob->nvars, &reoptrestart) );

      /* call initialization methods of plugins */
      SCIP_CALL( SCIPsetInitPlugins(scip->set, scip->mem->probmem, scip->stat) );

      /* install globally valid lower and upper bounds */
      SCIP_CALL( SCIPreoptInstallBounds(scip->reopt, scip->set, scip->stat, scip->transprob, scip->lp, scip->branchcand,
            scip->eventqueue, scip->cliquetable, scip->mem->probmem) );

      /* check, whether objective value is always integral by inspecting the problem, if it is the case adjust the
       * cutoff bound if primal solution is already known
       */
      SCIP_CALL( SCIPprobCheckObjIntegral(scip->transprob, scip->origprob, scip->mem->probmem, scip->set, scip->stat,
            scip->primal, scip->tree, scip->reopt, scip->lp, scip->eventfilter, scip->eventqueue) );

      /* if possible, scale objective function such that it becomes integral with gcd 1 */
      SCIP_CALL( SCIPprobScaleObj(scip->transprob, scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->primal,
            scip->tree, scip->reopt, scip->lp, scip->eventfilter, scip->eventqueue) );

      SCIPlpRecomputeLocalAndGlobalPseudoObjval(scip->lp, scip->set, scip->transprob);
   }

   /* try to compress the search tree */
   if( scip->set->compr_enable )
   {
      SCIP_CALL( compressReoptTree(scip) );
   }

   return SCIP_OKAY;
}

/** transforms and presolves problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages:
 *        - \ref SCIP_STAGE_PRESOLVING if the presolving process was interrupted
 *        - \ref SCIP_STAGE_PRESOLVED if the presolving process was finished and did not solve the problem
 *        - \ref SCIP_STAGE_SOLVED if the problem was solved during presolving
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPpresolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool unbounded;
   SCIP_Bool infeasible;
   SCIP_Bool vanished;
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPpresolve", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* start solving timer */
   SCIPclockStart(scip->stat->solvingtime, scip->set);
   SCIPclockStart(scip->stat->solvingtimeoverall, scip->set);

   /* capture the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptCapture(scip->interrupt);

   /* reset the user interrupt flag */
   scip->stat->userinterrupt = FALSE;
   SCIP_CALL( SCIPinterruptLP(scip, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* initialize solving data structures and transform problem */
      retcode =  SCIPtransformProb(scip);
      if( retcode != SCIP_OKAY )
      {
         SCIP_CALL( SCIPfreeTransform(scip) );
         return retcode;
      }

      assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);

      /*lint -fallthrough*/

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      /* presolve problem */
      SCIP_CALL( presolve(scip, &unbounded, &infeasible, &vanished) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED || scip->set->stage == SCIP_STAGE_PRESOLVING);

      if( infeasible || unbounded || vanished )
      {
         assert(scip->set->stage == SCIP_STAGE_PRESOLVED);

         /* initialize solving process data structures to be able to switch to SOLVED stage */
         SCIP_CALL( initSolve(scip, TRUE) );

         /* switch stage to SOLVED */
         scip->set->stage = SCIP_STAGE_SOLVED;

         /* print solution message */
         switch( scip->stat->status )/*lint --e{788}*/
         {
         case SCIP_STATUS_OPTIMAL:
            /* remove the root node from the tree, s.t. the lower bound is set to +infinity ???????????? (see initSolve())*/
            SCIP_CALL( SCIPtreeClear(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter, scip->eventqueue, scip->lp) );
            break;

         case SCIP_STATUS_INFEASIBLE:
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
               "presolving detected infeasibility\n");
            break;

         case SCIP_STATUS_UNBOUNDED:
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
               "presolving detected unboundedness\n");
            break;

         case SCIP_STATUS_INFORUNBD:
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
               "presolving detected unboundedness (or infeasibility)\n");
            break;

         default:
            /* note that this is in an internal SCIP error since the status is corrupted */
            SCIPerrorMessage("invalid SCIP status <%d>\n", scip->stat->status);
            SCIPABORT();
            return SCIP_ERROR; /*lint !e527*/
         }
      }
      else if( scip->set->stage == SCIP_STAGE_PRESOLVED )
      {
         int h;

         /* print presolved problem statistics */
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
            "presolved problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
            scip->transprob->nvars, scip->transprob->nbinvars, scip->transprob->nintvars, scip->transprob->nimplvars,
            scip->transprob->ncontvars, scip->transprob->nconss);

         for( h = 0; h < scip->set->nconshdlrs; ++h )
         {
            int nactiveconss;

            nactiveconss = SCIPconshdlrGetNActiveConss(scip->set->conshdlrs[h]);
            if( nactiveconss > 0 )
            {
               SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
                  "%7d constraints of type <%s>\n", nactiveconss, SCIPconshdlrGetName(scip->set->conshdlrs[h]));
            }
         }

         if( SCIPprobIsObjIntegral(scip->transprob) )
         {
            SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
               "transformed objective value is always integral (scale: %.15g)\n", scip->transprob->objscale);
         }
      }
      else
      {
         assert(scip->set->stage == SCIP_STAGE_PRESOLVING);
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "presolving was interrupted.\n");
      }

      /* display timing statistics */
      SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "Presolving Time: %.2f\n", SCIPclockGetTime(scip->stat->presolvingtime));
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVED:
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* release the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptRelease(scip->interrupt);

   /* stop solving timer */
   SCIPclockStop(scip->stat->solvingtime, scip->set);
   SCIPclockStop(scip->stat->solvingtimeoverall, scip->set);

   if( scip->set->stage == SCIP_STAGE_SOLVED )
   {
      /* display most relevant statistics */
      SCIP_CALL( displayRelevantStats(scip) );
   }

   return SCIP_OKAY;
}

/** transforms, presolves, and solves problem
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *        - \ref SCIP_STAGE_PRESOLVING if the solution process was interrupted during presolving
 *        - \ref SCIP_STAGE_SOLVING if the solution process was interrupted during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the solving process was not interrupted
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPsolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Longint cutpoolncutsfoundbeforerestart = 0;
   SCIP_Longint cutpoolncutsaddedbeforerestart = 0;
   SCIP_Longint cutpoolncallsbeforerestart = 0;
   SCIP_Longint cutpoolnrootcallsbeforerestart = 0;
   SCIP_Longint cutpoolmaxncutsbeforerestart = 0;
   SCIP_Real cutpooltimebeforerestart = 0;
   SCIP_Bool statsprinted = FALSE;
   SCIP_Bool restart;
   SCIP_Bool transferstatistics = FALSE;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolve", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* if the stage is already SCIP_STAGE_SOLVED do nothing */
   if( scip->set->stage == SCIP_STAGE_SOLVED )
      return SCIP_OKAY;

   if( scip->stat->status == SCIP_STATUS_INFEASIBLE || scip->stat->status == SCIP_STATUS_OPTIMAL || scip->stat->status == SCIP_STATUS_UNBOUNDED || scip->stat->status == SCIP_STATUS_INFORUNBD )
   {
      SCIPwarningMessage(scip, "SCIPsolve() was called, but problem is already solved\n");
      return SCIP_OKAY;
   }

   /* check, if a node selector exists */
   if( SCIPsetGetNodesel(scip->set, scip->stat) == NULL )
   {
      SCIPerrorMessage("no node selector available\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* check, if an integrality constraint handler exists if there are integral variables */
   if( (SCIPgetNBinVars(scip) >= 0 || SCIPgetNIntVars(scip) >= 0) && SCIPfindConshdlr(scip, "integral") == NULL )
   {
      SCIPwarningMessage(scip, "integrality constraint handler not available\n");
   }

   /* initialize presolving flag (may be modified in SCIPpresolve()) */
   scip->stat->performpresol = FALSE;

   /* if a decomposition exists and Benders' decomposition has been enabled, then a decomposition is performed */
   if( scip->set->stage == SCIP_STAGE_PROBLEM && SCIPdecompstoreGetNOrigDecomps(scip->decompstore) > 0
      && scip->set->decomp_applybenders && SCIPgetNActiveBenders(scip) == 0 )
   {
      int decompindex = 0;

      /* applying the Benders' decomposition */
      SCIP_CALL( SCIPapplyBendersDecomposition(scip, decompindex) );
   }

   /* start solving timer */
   SCIPclockStart(scip->stat->solvingtime, scip->set);
   SCIPclockStart(scip->stat->solvingtimeoverall, scip->set);

   /* capture the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptCapture(scip->interrupt);

   /* reset the user interrupt flag */
   scip->stat->userinterrupt = FALSE;
   SCIP_CALL( SCIPinterruptLP(scip, FALSE) );

   /* automatic restarting loop */
   restart = scip->stat->userrestart;

   do
   {
      if( restart )
      {
         transferstatistics = TRUE;
         cutpoolncutsfoundbeforerestart = SCIPcutpoolGetNCutsFound(scip->cutpool);
         cutpoolncutsaddedbeforerestart = SCIPcutpoolGetNCutsAdded(scip->cutpool);
         cutpooltimebeforerestart = SCIPcutpoolGetTime(scip->cutpool);
         cutpoolncallsbeforerestart = SCIPcutpoolGetNCalls(scip->cutpool);
         cutpoolnrootcallsbeforerestart = SCIPcutpoolGetNRootCalls(scip->cutpool);
         cutpoolmaxncutsbeforerestart = SCIPcutpoolGetMaxNCuts(scip->cutpool);

         /* free the solving process data in order to restart */
         assert(scip->set->stage == SCIP_STAGE_SOLVING);
         if( scip->stat->userrestart )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
               "(run %d, node %" SCIP_LONGINT_FORMAT ") performing user restart\n",
               scip->stat->nruns, scip->stat->nnodes);
         else
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
               "(run %d, node %" SCIP_LONGINT_FORMAT ") restarting after %d global fixings of integer variables\n",
               scip->stat->nruns, scip->stat->nnodes, scip->stat->nrootintfixingsrun);
         /* an extra blank line should be printed separately since the buffer message handler only handles up to one line
          * correctly */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
         /* reset relaxation solution, so that the objective value is recomputed from scratch next time, using the new
          * fixings which may be produced during the presolving after the restart */
         SCIP_CALL( SCIPclearRelaxSolVals(scip, NULL) );

         SCIP_CALL( freeSolve(scip, TRUE) );
         assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);
      }
      restart = FALSE;
      scip->stat->userrestart = FALSE;

      switch( scip->set->stage )
      {
      case SCIP_STAGE_PROBLEM:
      case SCIP_STAGE_TRANSFORMED:
      case SCIP_STAGE_PRESOLVING:
         /* initialize solving data structures, transform and problem */

         SCIP_CALL( SCIPpresolve(scip) );
         /* remember that we already printed the relevant statistics */
         if( scip->set->stage == SCIP_STAGE_SOLVED )
            statsprinted = TRUE;

         if( scip->set->stage == SCIP_STAGE_SOLVED || scip->set->stage == SCIP_STAGE_PRESOLVING )
         {
            if ( scip->set->reopt_enable )
            {
               SCIP_CALL( prepareReoptimization(scip) );
            }
            break;
         }
         assert(scip->set->stage == SCIP_STAGE_PRESOLVED);

         if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
            break;
         /*lint -fallthrough*/

      case SCIP_STAGE_PRESOLVED:
         /* check if reoptimization is enabled and global constraints are saved */
         if( scip->set->reopt_enable )
         {
            SCIP_CALL( prepareReoptimization(scip) );
         }

         /* initialize solving process data structures */
         SCIP_CALL( initSolve(scip, FALSE) );
         assert(scip->set->stage == SCIP_STAGE_SOLVING);
         SCIPmessagePrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL, "\n");

         /*lint -fallthrough*/

      case SCIP_STAGE_SOLVING:
         /* reset display */
         SCIPstatResetDisplay(scip->stat);

         /* remember cutpool statistics after restart */
         if( transferstatistics )
         {
            SCIPcutpoolAddNCutsFound(scip->cutpool, cutpoolncutsfoundbeforerestart);
            SCIPcutpoolAddNCutsAdded(scip->cutpool, cutpoolncutsaddedbeforerestart);
            SCIPcutpoolSetTime(scip->cutpool, cutpooltimebeforerestart);
            SCIPcutpoolAddNCalls(scip->cutpool, cutpoolncallsbeforerestart);
            SCIPcutpoolAddNRootCalls(scip->cutpool, cutpoolnrootcallsbeforerestart);
            SCIPcutpoolAddMaxNCuts(scip->cutpool, cutpoolmaxncutsbeforerestart);
         }

         /* continue solution process */
         SCIP_CALL( SCIPsolveCIP(scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->mem, scip->origprob, scip->transprob,
               scip->primal, scip->tree, scip->reopt, scip->lp, scip->relaxation, scip->pricestore, scip->sepastore,
               scip->cutpool, scip->delayedcutpool, scip->branchcand, scip->conflict, scip->conflictstore,
               scip->eventfilter, scip->eventqueue, scip->cliquetable, &restart) );

         /* detect, whether problem is solved */
         if( SCIPtreeGetNNodes(scip->tree) == 0 && SCIPtreeGetCurrentNode(scip->tree) == NULL )
         {
            assert(scip->stat->status == SCIP_STATUS_OPTIMAL
               || scip->stat->status == SCIP_STATUS_INFEASIBLE
               || scip->stat->status == SCIP_STATUS_UNBOUNDED
               || scip->stat->status == SCIP_STATUS_INFORUNBD);
            assert(!restart);

            /* tree is empty, and no current node exists -> problem is solved */
            scip->set->stage = SCIP_STAGE_SOLVED;
         }
         break;

      case SCIP_STAGE_SOLVED:
         assert(scip->stat->status == SCIP_STATUS_OPTIMAL
            || scip->stat->status == SCIP_STATUS_INFEASIBLE
            || scip->stat->status == SCIP_STATUS_UNBOUNDED
            || scip->stat->status == SCIP_STATUS_INFORUNBD);

         break;

      default:
         SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
         return SCIP_INVALIDCALL;
      }  /*lint !e788*/
   }
   while( restart && !SCIPsolveIsStopped(scip->set, scip->stat, TRUE) );

   /* we have to store all unprocessed nodes if reoptimization is enabled */
   if( scip->set->reopt_enable && scip->set->stage != SCIP_STAGE_PRESOLVING
    && SCIPsolveIsStopped(scip->set, scip->stat, TRUE) )
   {
      /* save unprocessed nodes */
      if( SCIPgetNNodesLeft(scip) > 0 )
      {
         SCIP_NODE** leaves;
         SCIP_NODE** children;
         SCIP_NODE** siblings;
         int nleaves;
         int nchildren;
         int nsiblings;

         /* get all open leave nodes */
         SCIP_CALL( SCIPgetLeaves(scip, &leaves, &nleaves) );

         /* get all open children nodes */
         SCIP_CALL( SCIPgetChildren(scip, &children, &nchildren) );

         /* get all open sibling nodes */
         SCIP_CALL( SCIPgetSiblings(scip, &siblings, &nsiblings) );

         /* add all open node to the reoptimization tree */
         SCIP_CALL( SCIPreoptSaveOpenNodes(scip->reopt, scip->set, scip->lp, scip->mem->probmem, leaves, nleaves,
               children, nchildren, siblings, nsiblings) );
      }
   }

   /* release the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptRelease(scip->interrupt);

   if( scip->set->reopt_enable )
   {
      /* save found solutions */
      int nsols;
      int s;

      nsols = scip->set->reopt_savesols == -1 ? INT_MAX : MAX(scip->set->reopt_savesols, 1);
      nsols = MIN(scip->primal->nsols, nsols);

      for( s = 0; s < nsols; s++ )
      {
         SCIP_SOL* sol;
         SCIP_Bool added;

         sol = scip->primal->sols[s];
         assert(sol != NULL);

         if( !SCIPsolIsOriginal(sol) )
         {
            SCIP_Bool hasinfval;

            /* retransform solution into the original problem space */
            SCIP_CALL( SCIPsolRetransform(sol, scip->set, scip->stat, scip->origprob, scip->transprob, &hasinfval) );
         }

         if( SCIPsolGetNodenum(sol) > 0 || SCIPsolGetHeur(sol) != NULL || (s == 0 && scip->set->reopt_sepabestsol) )
         {
            /* if the best solution should be separated, we must not store it in the solution tree */
            if( s == 0 && scip->set->reopt_sepabestsol )
            {
               SCIP_CALL( SCIPreoptAddOptSol(scip->reopt, sol, scip->mem->probmem, scip->set, scip->stat, scip->origprimal,
                     scip->origprob->vars, scip->origprob->nvars) );
            }
            /* add solution to solution tree */
            else
            {
               SCIPdebugMsg(scip, "try to add solution to the solution tree:\n");
               SCIPdebug( SCIP_CALL( SCIPsolPrint(sol, scip->set, scip->messagehdlr, scip->stat, scip->origprob, \
                     scip->transprob, NULL, FALSE, FALSE) ); );

               SCIP_CALL( SCIPreoptAddSol(scip->reopt, scip->set, scip->stat, scip->origprimal, scip->mem->probmem,
                     sol, s == 0, &added, scip->origprob->vars, scip->origprob->nvars, scip->stat->nreoptruns) );
            }
         }
      }

      SCIPdebugMsg(scip, "-> saved %d solution.\n", nsols);

      /* store variable history */
      if( scip->set->reopt_storevarhistory )
      {
         SCIP_CALL( SCIPreoptUpdateVarHistory(scip->reopt, scip->set, scip->stat, scip->mem->probmem,
               scip->origprob->vars, scip->origprob->nvars) );
      }
   }

   /* stop solving timer */
   SCIPclockStop(scip->stat->solvingtime, scip->set);
   SCIPclockStop(scip->stat->solvingtimeoverall, scip->set);

   /* decrease time limit during reoptimization */
   if( scip->set->reopt_enable && scip->set->reopt_commontimelimit )
   {
      SCIP_Real timelimit;
      SCIP_Real usedtime;

      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      usedtime = SCIPgetSolvingTime(scip);
      timelimit = timelimit - usedtime;
      timelimit = MAX(0, timelimit);

      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit) );
   }

   if( !statsprinted )
   {
      /* display most relevant statistics */
      SCIP_CALL( displayRelevantStats(scip) );
   }

   return SCIP_OKAY;
}

/** transforms, presolves, and solves problem using the configured concurrent solvers
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *        - \ref SCIP_STAGE_PRESOLVING if the solution process was interrupted during presolving
 *        - \ref SCIP_STAGE_SOLVING if the solution process was interrupted during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the solving process was not interrupted
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @deprecated Please use SCIPsolveConcurrent() instead.
 */
SCIP_RETCODE SCIPsolveParallel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveParallel", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPsolveConcurrent(scip);
}

/** transforms, presolves, and solves problem using the configured concurrent solvers
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *        - \ref SCIP_STAGE_PRESOLVING if the solution process was interrupted during presolving
 *        - \ref SCIP_STAGE_SOLVING if the solution process was interrupted during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the solving process was not interrupted
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPsolveConcurrent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
#ifdef TPI_NONE
   SCIPinfoMessage(scip, NULL, "SCIP was compiled without task processing interface. Parallel solve not possible\n");
   return SCIP_OKAY;
#else
   SCIP_RETCODE     retcode;
   int              i;
   SCIP_RANDNUMGEN* rndgen;
   int              minnthreads;
   int              maxnthreads;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveConcurrent", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsetIntParam(scip, "timing/clocktype", SCIP_CLOCKTYPE_WALL) );

   minnthreads = scip->set->parallel_minnthreads;
   maxnthreads = scip->set->parallel_maxnthreads;

   if( minnthreads > maxnthreads )
   {
      SCIPerrorMessage("minimum number of threads greater than maximum number of threads\n");
      return SCIP_INVALIDDATA;
   }
   if( scip->concurrent == NULL )
   {
      int                   nconcsolvertypes;
      SCIP_CONCSOLVERTYPE** concsolvertypes;
      SCIP_Longint          nthreads;
      SCIP_Real             memorylimit;
      int*                  solvertypes;
      SCIP_Longint*         weights;
      SCIP_Real*            prios;
      int                   ncandsolvertypes;
      SCIP_Real             prefpriosum;

      /* check if concurrent solve is configured to presolve the problem
       * before setting up the concurrent solvers
       */
      if( scip->set->concurrent_presolvebefore )
      {
         /* if yes, then presolve the problem */
         SCIP_CALL( SCIPpresolve(scip) );
         if( SCIPgetStatus(scip) >= SCIP_STATUS_OPTIMAL )
            return SCIP_OKAY;
      }
      else
      {
         SCIP_Bool infeas;

         /* if not, transform the problem and switch stage to presolved */
         SCIP_CALL( SCIPtransformProb(scip) );
         SCIP_CALL( initPresolve(scip) );
         SCIP_CALL( exitPresolve(scip, TRUE, &infeas) );
         assert(!infeas);
      }

      /* the presolving must have run into a limit, so we stop here */
      if( scip->set->stage < SCIP_STAGE_PRESOLVED )
      {
         SCIP_CALL( displayRelevantStats(scip) );
         return SCIP_OKAY;
      }

      nthreads = INT_MAX;
      /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
      memorylimit = scip->set->limit_memory;
      if( memorylimit < SCIP_MEM_NOLIMIT )
      {
         memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
         memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
         /* estimate maximum number of copies that be created based on memory limit */
         if( !scip->set->misc_avoidmemout )
         {
            nthreads = MAX(1, memorylimit / (4.0*SCIPgetMemExternEstim(scip)/1048576.0));
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "estimated a maximum of %lli threads based on memory limit\n", nthreads);
         }
         else
         {
            nthreads = minnthreads;
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "ignoring memory limit; all threads can be created\n");
         }
      }
      nconcsolvertypes = SCIPgetNConcsolverTypes(scip);
      concsolvertypes = SCIPgetConcsolverTypes(scip);

      if( minnthreads > nthreads )
      {
         SCIP_CALL( initSolve(scip, TRUE) );
         scip->stat->status = SCIP_STATUS_MEMLIMIT;
         SCIPsyncstoreSetSolveIsStopped(SCIPgetSyncstore(scip), TRUE);
         SCIPwarningMessage(scip, "requested minimum number of threads could not be satisfied with given memory limit\n");
         SCIP_CALL( displayRelevantStats(scip) );
         return SCIP_OKAY;
      }

      if( nthreads == 1 )
      {
         SCIPwarningMessage(scip, "can only use 1 thread, doing sequential solve instead\n");
         SCIP_CALL( SCIPfreeConcurrent(scip) );
         return SCIPsolve(scip);
      }
      nthreads = MIN(nthreads, maxnthreads);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "using %lli threads for concurrent solve\n", nthreads);

      /* now set up nthreads many concurrent solvers that will be used for the concurrent solve
       * using the preferred priorities of each concurrent solver
       */
      prefpriosum = 0.0;
      for( i = 0; i < nconcsolvertypes; ++i )
         prefpriosum += SCIPconcsolverTypeGetPrefPrio(concsolvertypes[i]);

      ncandsolvertypes = 0;
      SCIP_CALL( SCIPallocBufferArray(scip, &solvertypes, nthreads + nconcsolvertypes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &weights, nthreads + nconcsolvertypes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &prios, nthreads + nconcsolvertypes) );
      for( i = 0; i < nconcsolvertypes; ++i )
      {
         int j;
         SCIP_Real prio;
         prio = nthreads * SCIPconcsolverTypeGetPrefPrio(concsolvertypes[i]) / prefpriosum;
         while( prio > 0.0 )
         {
            j = ncandsolvertypes++;
            assert(j < 2*nthreads);
            weights[j] = 1;
            solvertypes[j] = i;
            prios[j] = MIN(1.0, prio);
            prio = prio - 1.0;
         }
      }
      /* select nthreads many concurrent solver types to create instances
       * according to the preferred prioriteis the user has set
       * This basically corresponds to a knapsack problem
       * with unit weights and capacity nthreads, where the profits are
       * the unrounded fraction of the total number of threads to be used.
       */
      SCIPselectDownRealInt(prios, solvertypes, nthreads, ncandsolvertypes);

      SCIP_CALL( SCIPcreateRandom(scip, &rndgen, (unsigned) scip->set->concurrent_initseed, TRUE) );
      for( i = 0; i < nthreads; ++i )
      {
         SCIP_CONCSOLVER* concsolver;

         SCIP_CALL( SCIPconcsolverCreateInstance(scip->set, concsolvertypes[solvertypes[i]], &concsolver) );
         if( scip->set->concurrent_changeseeds && SCIPgetNConcurrentSolvers(scip) > 1 )
            SCIP_CALL( SCIPconcsolverInitSeeds(concsolver, SCIPrandomGetInt(rndgen, 0, INT_MAX)) );
      }
      SCIPfreeRandom(scip, &rndgen);
      SCIPfreeBufferArray(scip, &prios);
      SCIPfreeBufferArray(scip, &weights);
      SCIPfreeBufferArray(scip, &solvertypes);

      assert(SCIPgetNConcurrentSolvers(scip) == nthreads);

      SCIP_CALL( SCIPsyncstoreInit(scip) );
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVED )
   {
      /* switch stage to solving */
      SCIP_CALL( initSolve(scip, TRUE) );
   }

   SCIPclockStart(scip->stat->solvingtime, scip->set);
   retcode = SCIPconcurrentSolve(scip);
   SCIPclockStop(scip->stat->solvingtime, scip->set);
   SCIP_CALL( displayRelevantStats(scip) );

   return retcode;
#endif
}

/** include specific heuristics and branching rules for reoptimization
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPenableReoptimization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             enable              /**< enable reoptimization (TRUE) or disable it (FALSE) */
   )
{
   assert(scip != NULL);

   /* we want to skip if nothing has changed */
   if( (enable && scip->set->reopt_enable && scip->reopt != NULL)
    || (!enable && !scip->set->reopt_enable && scip->reopt == NULL) )
      return SCIP_OKAY;

   /* check stage and throw an error if we try to disable reoptimization during the solving process.
    *
    * @note the case that we will disable the reoptimization and have already performed presolving can only happen if
    *       we are try to solve a general MIP
    *
    * @note this fix is only for the bugfix release 3.2.1, in the next major release reoptimization can be used for
    *       general MIPs, too.
    */
   if( scip->set->stage > SCIP_STAGE_PROBLEM && !(!enable && scip->set->stage == SCIP_STAGE_PRESOLVED) )
   {
      SCIPerrorMessage("Reoptimization cannot be %s after starting the (pre)solving process.\n", enable ? "enabled" : "disabled");
      return SCIP_INVALIDCALL;
   }

   /* if the current stage is SCIP_STAGE_PROBLEM we have to include the heuristics and branching rule */
   if( scip->set->stage == SCIP_STAGE_PROBLEM || (!enable && scip->set->stage == SCIP_STAGE_PRESOLVED) )
   {
      /* initialize all reoptimization data structures */
      if( enable && scip->reopt == NULL )
      {
         /* set enable flag */
         scip->set->reopt_enable = enable;

         SCIP_CALL( SCIPreoptCreate(&scip->reopt, scip->set, scip->mem->probmem) );
         SCIP_CALL( SCIPsetSetReoptimizationParams(scip->set, scip->messagehdlr) );
      }
      /* disable all reoptimization plugins and free the structure if necessary */
      else if( (!enable && scip->reopt != NULL) || (!enable && scip->set->reopt_enable && scip->reopt == NULL) )
      {
         /* set enable flag */
         scip->set->reopt_enable = enable;

         if( scip->reopt != NULL )
         {
            SCIP_CALL( SCIPreoptFree(&(scip->reopt), scip->set, scip->origprimal, scip->mem->probmem) );
            assert(scip->reopt == NULL);
         }
         SCIP_CALL( SCIPsetSetReoptimizationParams(scip->set, scip->messagehdlr) );
      }
   }
   else
   {
      /* set enable flag */
      scip->set->reopt_enable = enable;
   }

   return SCIP_OKAY;
}

/** save bound change based on dual information in the reoptimization tree
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPaddReoptDualBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             newbound,           /**< new bound of the variable */
   SCIP_Real             oldbound            /**< old bound of the variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddReoptDualBndchg", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert(SCIPsetIsFeasLT(scip->set, newbound, oldbound) || SCIPsetIsFeasGT(scip->set, newbound, oldbound));

   SCIP_CALL( SCIPreoptAddDualBndchg(scip->reopt, scip->set, scip->mem->probmem, node, var, newbound, oldbound) );

   return SCIP_OKAY;
}

/** returns the optimal solution of the last iteration or NULL of none exists */
SCIP_SOL* SCIPgetReoptLastOptSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SOL* sol;

   assert(scip != NULL);

   sol = NULL;

   if( scip->set->reopt_enable && scip->stat->nreoptruns > 1 )
   {
      sol = SCIPreoptGetLastBestSol(scip->reopt);
   }

   return sol;
}

/** returns the objective coefficent of a given variable in a previous iteration
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetReoptOldObjCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable */
   int                   run,                /**< number of the run */
   SCIP_Real*            objcoef             /**< pointer to store the objective coefficient */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(0 < run && run <= scip->stat->nreoptruns);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetReoptOldObjCoef", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPvarIsOriginal(var) )
      *objcoef = SCIPreoptGetOldObjCoef(scip->reopt, run, SCIPvarGetIndex(var));
   else
   {
      SCIP_VAR* origvar;
      SCIP_Real constant;
      SCIP_Real scalar;

      assert(SCIPvarIsActive(var));

      origvar = var;
      constant = 0.0;
      scalar = 1.0;

      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );
      assert(origvar != NULL);
      assert(SCIPvarIsOriginal(origvar));

      *objcoef = SCIPreoptGetOldObjCoef(scip->reopt, run, SCIPvarGetIndex(origvar));
   }
   return SCIP_OKAY;
}

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post If this method is called in \SCIP stage \ref SCIP_STAGE_INIT or \ref SCIP_STAGE_PROBLEM, the stage of
 *        \SCIP is not changed; otherwise, the \SCIP stage is changed to \ref SCIP_STAGE_TRANSFORMED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPfreeSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             restart             /**< should certain data be preserved for improved restarting? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeSolve", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PROBLEM:
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   {
      SCIP_Bool infeasible;

      assert(scip->stat->status != SCIP_STATUS_INFEASIBLE);
      assert(scip->stat->status != SCIP_STATUS_INFORUNBD);
      assert(scip->stat->status != SCIP_STATUS_UNBOUNDED);
      assert(scip->stat->status != SCIP_STATUS_OPTIMAL);

      /* exit presolving */
      SCIP_CALL( exitPresolve(scip, FALSE, &infeasible) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED);
   }

   /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
      /* switch stage to TRANSFORMED */
      scip->set->stage = SCIP_STAGE_TRANSFORMED;
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* free solution process data structures */
      SCIP_CALL( freeSolve(scip, restart) );
      assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post If this method is called in \SCIP stage \ref SCIP_STAGE_INIT, \ref SCIP_STAGE_TRANSFORMED or \ref SCIP_STAGE_PROBLEM,
 *        the stage of \SCIP is not changed; otherwise, the \SCIP stage is changed to \ref SCIP_STAGE_PRESOLVED.
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPfreeReoptSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeReoptSolve", TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_PROBLEM:
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   {
      SCIP_Bool infeasible;

      assert(scip->stat->status != SCIP_STATUS_INFEASIBLE);
      assert(scip->stat->status != SCIP_STATUS_INFORUNBD);
      assert(scip->stat->status != SCIP_STATUS_UNBOUNDED);
      assert(scip->stat->status != SCIP_STATUS_OPTIMAL);

      /* exit presolving */
      SCIP_CALL( exitPresolve(scip, FALSE, &infeasible) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED);

      return SCIP_OKAY;
   }

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* free solution process data structures */
      SCIP_CALL( freeReoptSolve(scip) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED);
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** frees all solution process data including presolving and transformed problem, only original problem is kept
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @post After calling this method \SCIP reaches one of the following stages:
 *        - \ref SCIP_STAGE_INIT if the method was called from \SCIP stage \ref SCIP_STAGE_INIT
 *        - \ref SCIP_STAGE_PROBLEM if the method was called from any other of the allowed stages
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPfreeTransform(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPfreeTransform", TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* release variables and constraints captured by reoptimization */
   if( scip->reopt != NULL )
   {
      SCIP_CALL( SCIPreoptReleaseData(scip->reopt, scip->set, scip->mem->probmem) );
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   {
      SCIP_Bool infeasible;

      assert(scip->stat->status != SCIP_STATUS_INFEASIBLE);
      assert(scip->stat->status != SCIP_STATUS_INFORUNBD);
      assert(scip->stat->status != SCIP_STATUS_UNBOUNDED);
      assert(scip->stat->status != SCIP_STATUS_OPTIMAL);

      /* exit presolving */
      SCIP_CALL( exitPresolve(scip, FALSE, &infeasible) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED);
   }

   /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* the solve was already freed, we directly go to freeTransform() */
      if( !scip->set->reopt_enable || scip->set->stage != SCIP_STAGE_PRESOLVED )
      {
         /* free solution process data */
         SCIP_CALL( SCIPfreeSolve(scip, FALSE) );
         assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);
      }
      /*lint -fallthrough*/

   case SCIP_STAGE_TRANSFORMED:
      /* free transformed problem data structures */
      SCIP_CALL( freeTransform(scip) );
      assert(scip->set->stage == SCIP_STAGE_PROBLEM);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
      assert(scip->set->stage == SCIP_STAGE_TRANSFORMING);
      SCIP_CALL( freeTransforming(scip) );
      assert(scip->set->stage == SCIP_STAGE_PROBLEM);
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** informs \SCIP that the solving process should be interrupted as soon as possible (e.g., after the current node has
 *   been solved)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note the \SCIP stage does not get changed
 */
SCIP_RETCODE SCIPinterruptSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPinterruptSolve", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   /* set the userinterrupt flag */
   scip->stat->userinterrupt = TRUE;

   return SCIP_OKAY;
}

/** indicates whether \SCIP has been informed that the solving process should be interrupted as soon as possible
 *
 *  This function returns whether SCIPinterruptSolve() has been called, which is different from SCIPinterrupted(),
 *  which returns whether a SIGINT signal has been received by the SCIP signal handler.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note the \SCIP stage does not get changed
 */
SCIP_Bool SCIPisSolveInterrupted(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisSolveInterrupted", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->userinterrupt;
}

/** informs SCIP that the solving process should be restarted as soon as possible (e.g., after the current node has
 *  been solved)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the \SCIP stage does not get changed
 */
SCIP_RETCODE SCIPrestartSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPrestartSolve", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* set the userrestart flag */
   scip->stat->userrestart = TRUE;

   return SCIP_OKAY;
}

/** returns whether reoptimization is enabled or not */
SCIP_Bool SCIPisReoptEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->set->reopt_enable;
}

/** returns the stored solutions corresponding to a given run */
SCIP_RETCODE SCIPgetReoptSolsRun(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   run,                /**< number of the run */
   SCIP_SOL**            sols,               /**< array to store solutions */
   int                   solssize,           /**< size of the array */
   int*                  nsols               /**< pointer to store number of solutions */
   )
{
   assert(scip != NULL);
   assert(sols != NULL);
   assert(solssize > 0);

   if( scip->set->reopt_enable )
   {
      assert(run > 0 && run <= scip->stat->nreoptruns);
      SCIP_CALL( SCIPreoptGetSolsRun(scip->reopt, run, sols, solssize, nsols) );
   }
   else
   {
      *nsols = 0;
   }

   return SCIP_OKAY;
}

/** mark all stored solutions as not updated */
void SCIPresetReoptSolMarks(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set->reopt_enable);
   assert(scip->reopt != NULL);

   if( scip->set->reopt_enable )
   {
      assert(scip->reopt != NULL);
      SCIPreoptResetSolMarks(scip->reopt);
   }
}

/** check if the reoptimization process should be restarted
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcheckReoptRestart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< current node of the branch and bound tree (or NULL) */
   SCIP_Bool*            restart             /**< pointer to store of the reoptimitation process should be restarted */
   )
{
   assert(scip != NULL);
   assert(scip->set->reopt_enable);
   assert(scip->reopt != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcheckReoptRestart", FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPreoptCheckRestart(scip->reopt, scip->set, scip->mem->probmem, node, scip->transprob->vars,
         scip->transprob->nvars, restart) );

   return SCIP_OKAY;
}

/** returns whether we are in the restarting phase
 *
 *  @return TRUE, if we are in the restarting phase; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
SCIP_Bool SCIPisInRestart(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisInRestart", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   /* return the restart status */
   return scip->stat->inrestart;
}
