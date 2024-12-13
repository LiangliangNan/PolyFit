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

/**@file   scip_copy.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for problem copies
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
#include "scip/benders.h"
#include "scip/clock.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/dcmp.h"
#include "scip/debug.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/pub_cons.h"
#include "scip/pub_cutpool.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_nlpi.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/syncstore.h"
#include "scip/var.h"

/** returns true if the @p cut matches the selection criterium for copying */
static
SCIP_Bool takeCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUT*             cut,                /**< a cut */
   char                  cutsel              /**< cut selection for sub SCIPs  ('a'ge, activity 'q'uotient) */
   )
{
   SCIP_Bool takecut;

   assert(cut != NULL);

   if( !SCIProwIsInLP(SCIPcutGetRow(cut)) )
      return FALSE;

   switch( cutsel )
   {
   case 'a':
      takecut = (SCIPcutGetAge(cut) == 0);
      break;
   case 'q':
      takecut = (SCIPcutGetLPActivityQuot(cut) >= scip->set->sepa_minactivityquot);
      break;
   default:
      SCIPerrorMessage("unknown cut selection strategy %c, must be either 'a' or 'q'\n", cutsel);
      SCIPABORT();
      takecut = FALSE;  /*lint !e527*/
      break;
   }

   return takecut;
}

/** copy active and tight cuts from one SCIP instance to linear constraints of another SCIP instance */
static
SCIP_RETCODE copyCuts(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_CUT**            cuts,               /**< cuts to copy */
   int                   ncuts,              /**< number of cuts to copy */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   int*                  ncutsadded          /**< pointer to store number of copied cuts */
   )
{
   int c;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(cuts != NULL || ncuts == 0);
   assert(ncutsadded != NULL);

   for( c = 0; c < ncuts; ++c )
   {
      SCIP_ROW* row;
      SCIP_Bool takecut;

      assert( cuts[c] != NULL ); /*lint !e613*/
      row = SCIPcutGetRow(cuts[c]); /*lint !e613*/
      assert(!SCIProwIsLocal(row));
      assert(!SCIProwIsModifiable(row));

      /* in case of a restart, convert the cuts with a good LP activity quotient; in other cases, e.g., when heuristics
       * copy cuts into subscips, take only currently active ones
       */
      if( sourcescip == targetscip )
      {
         assert( SCIPisInRestart(sourcescip) );
         takecut = takeCut(sourcescip, cuts[c], sourcescip->set->sepa_cutselrestart); /*lint !e613*/
      }
      else
         takecut = takeCut(sourcescip, cuts[c], sourcescip->set->sepa_cutselsubscip); /*lint !e613*/

      /* create a linear constraint out of the cut */
      if( takecut )
      {
         char name[SCIP_MAXSTRLEN];
         SCIP_CONS* cons;
         SCIP_COL** cols;
         SCIP_VAR** vars;
         int ncols;
         int i;

         cols = SCIProwGetCols(row);
         ncols = SCIProwGetNNonz(row);

         /* get all variables of the row */
         SCIP_CALL( SCIPallocBufferArray(targetscip, &vars, ncols) );
         for( i = 0; i < ncols && takecut; ++i )
         {
            vars[i] = SCIPcolGetVar(cols[i]);
            takecut = !SCIPvarIsRelaxationOnly(vars[i]);
         }

         /* discard cut if it contains a variable which is invalid after a restart */
         if( !takecut )
         {
            /* free temporary memory */
            SCIPfreeBufferArray(targetscip, &vars);
            continue;
         }

         /* get corresponding variables in targetscip if necessary */
         if( sourcescip != targetscip )
         {
            SCIP_Bool success;

            for( i = 0; i < ncols; ++i )
            {
               SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, vars[i], &vars[i], varmap, consmap, global, &success) );

               if( !success )
               {
                  SCIPdebugMsg(sourcescip, "Converting cuts to constraints failed.\n");

                  /* free temporary memory */
                  SCIPfreeBufferArray(targetscip, &vars);
                  return SCIP_OKAY;
               }
            }
         }

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d", SCIProwGetName(row), SCIPgetNRuns(sourcescip));
         SCIP_CALL( SCIPcreateConsLinear(targetscip, &cons, name, ncols, vars, SCIProwGetVals(row),
               SCIProwGetLhs(row) - SCIProwGetConstant(row), SCIProwGetRhs(row) - SCIProwGetConstant(row),
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(targetscip, cons) );

         SCIPdebugMsg(sourcescip, "Converted cut <%s> to constraint <%s>.\n", SCIProwGetName(row), SCIPconsGetName(cons));
         SCIPdebugPrintCons(targetscip, cons, NULL);
         SCIP_CALL( SCIPreleaseCons(targetscip, &cons) );

         /* free temporary memory */
         SCIPfreeBufferArray(targetscip, &vars);

         ++(*ncutsadded);
      }
   }

   return SCIP_OKAY;
}

/** copies plugins from sourcescip to targetscip; in case that a constraint handler which does not need constraints
 *  cannot be copied, valid will return FALSE. All plugins can declare that, if their copy process failed, the
 *  copied SCIP instance might not represent the same problem semantics as the original.
 *  Note that in this case dual reductions might be invalid.
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *        Also, 'passmessagehdlr' should be set to FALSE.
 *
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @note This method does not copy Benders' plugins. To this end, the method SCIPcopyBenders() must be called
 *        separately.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method targetscip reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyPlugins(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_Bool             copyreaders,        /**< should the file readers be copied */
   SCIP_Bool             copypricers,        /**< should the variable pricers be copied */
   SCIP_Bool             copyconshdlrs,      /**< should the constraint handlers be copied */
   SCIP_Bool             copyconflicthdlrs,  /**< should the conflict handlers be copied */
   SCIP_Bool             copypresolvers,     /**< should the presolvers be copied */
   SCIP_Bool             copyrelaxators,     /**< should the relaxation handlers be copied */
   SCIP_Bool             copyseparators,     /**< should the separators be copied */
   SCIP_Bool             copycutselectors,   /**< should the cut selectors be copied */
   SCIP_Bool             copypropagators,    /**< should the propagators be copied */
   SCIP_Bool             copyheuristics,     /**< should the heuristics be copied */
   SCIP_Bool             copyeventhdlrs,     /**< should the event handlers be copied */
   SCIP_Bool             copynodeselectors,  /**< should the node selectors be copied */
   SCIP_Bool             copybranchrules,    /**< should the branchrules be copied */
   SCIP_Bool             copydisplays,       /**< should the display columns be copied */
   SCIP_Bool             copydialogs,        /**< should the dialogs be copied */
   SCIP_Bool             copytables,         /**< should the statistics tables be copied */
   SCIP_Bool             copyexprhdlrs,      /**< should the expression handlers be copied */
   SCIP_Bool             copynlpis,          /**< should the NLPIs be copied */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether plugins, in particular all constraint
                                              *   handlers which do not need constraints were validly copied */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourcescip->set != NULL);
   assert(targetscip->set != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyPlugins", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyPlugins", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   /* passes the message handler of the source SCIP to the target SCIP, also if NULL */
   if( passmessagehdlr )
   {
      SCIP_CALL( SCIPsetMessagehdlr(targetscip, SCIPgetMessagehdlr(sourcescip)) );
   }

   SCIP_CALL( SCIPsetCopyPlugins(sourcescip->set, targetscip->set,
         copyreaders, copypricers, copyconshdlrs, copyconflicthdlrs, copypresolvers, copyrelaxators, copyseparators, copycutselectors, copypropagators,
         copyheuristics, copyeventhdlrs, copynodeselectors, copybranchrules, copydisplays, copydialogs, copytables, copyexprhdlrs, copynlpis, valid) );

   return SCIP_OKAY;
}

/** copies all Benders' decomposition plugins
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note the 'threadsafe' parameter should only be set to TRUE if you are absolutely certain that the source and target
 *        SCIP instances will be solved in parallel. The usual case is to set this to FALSE, since thread safety
 *        typically incurs a performance cost.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @post After calling this method targetscip reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyBenders(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables; if NULL the transfer of cuts is not possible */
   SCIP_Bool             threadsafe,         /**< FALSE, if data can be safely shared between the source and target
                                              *   SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool*            valid               /**< pointer to store whether all plugins were validly copied */
   )
{
   /* TODO: If the Benders' decomposition is not copied, then cons_benders needs to be deactivated. */
   SCIP_Bool copybendersvalid;
   int p;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourcescip != targetscip);
   assert(sourcescip->set != NULL);
   assert(targetscip->set != NULL);
   assert(valid != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyBenders", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyBenders", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   *valid = TRUE;

   if( sourcescip->set->benders != NULL )
   {
      for( p = sourcescip->set->nbenders - 1; p >= 0; --p )
      {
         copybendersvalid = FALSE;
         SCIP_CALL( SCIPbendersCopyInclude(sourcescip->set->benders[p], sourcescip->set, targetscip->set, varmap,
               threadsafe, &copybendersvalid) );
         *valid = *valid && copybendersvalid;
      }
   }

   return SCIP_OKAY;
}

/** create a problem by copying the problem data of the source SCIP */
static
SCIP_RETCODE copyProb(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL  */
   SCIP_Bool             original,           /**< should the original problem be copied? */
   SCIP_Bool             global,             /**< create a global or a local copy? (set to TRUE for original copy) */
   const char*           name                /**< problem name of target */
   )
{
   SCIP_PROB* sourceprob;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(!original || global);
   assert(original || SCIPisTransformed(sourcescip));

   /* free old problem */
   SCIP_CALL( SCIPfreeProb(targetscip) );
   assert(targetscip->set->stage == SCIP_STAGE_INIT);

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPgetNVars(sourcescip)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPgetNConss(sourcescip)) );
   }
   else
      localconsmap = consmap;

   /* switch stage to PROBLEM */
   targetscip->set->stage = SCIP_STAGE_PROBLEM;

   if( original )
      sourceprob = sourcescip->origprob;
   else
      sourceprob = sourcescip->transprob;

   /* create the statistics data structure */
   SCIP_CALL( SCIPstatCreate(&targetscip->stat, targetscip->mem->probmem, targetscip->set, targetscip->transprob, targetscip->origprob, targetscip->messagehdlr) );
   targetscip->stat->subscipdepth = sourcescip->stat->subscipdepth + 1;

   /* create the problem by copying the source problem */
   SCIP_CALL( SCIPprobCopy(&targetscip->origprob, targetscip->mem->probmem, targetscip->set, name, sourcescip, sourceprob, localvarmap, localconsmap, global) );

   /* creating the solution candidates storage */
   /**@todo copy solution of source SCIP as candidates for the target SCIP */
   SCIP_CALL( SCIPprimalCreate(&targetscip->origprimal) );

   /* create conflict store to store conflict constraints */
   SCIP_CALL( SCIPconflictstoreCreate(&targetscip->conflictstore, targetscip->set) );

   SCIP_CALL( SCIPdecompstoreCreate(&targetscip->decompstore, SCIPblkmem(targetscip), SCIP_DECOMPSTORE_CAPA) );

   SCIP_CALL( SCIPdebugSolDataCreate(&targetscip->set->debugsoldata) );

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}


/** create a problem by copying the problem data of the source SCIP
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
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
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method targetscip reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyProb(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   const char*           name                /**< problem name of target */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyProb", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyProb", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );

   SCIP_CALL( copyProb(sourcescip, targetscip, varmap, consmap, FALSE, global, name) );

   return SCIP_OKAY;
}

/** create a problem by copying the original problem data of the source SCIP
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method targetscip reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyOrigProb(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL  */
   const char*           name                /**< problem name of target */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyOrigProb", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyOrigProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   SCIP_CALL( copyProb(sourcescip, targetscip, varmap, consmap, TRUE, TRUE, name) );

   /* set the correct objective sense; necessary if we maximize in the original problem */
   SCIP_CALL( SCIPsetObjsense(targetscip, SCIPgetObjsense(sourcescip)) );

   /* set the objective offset */
   SCIP_CALL( SCIPaddOrigObjoffset(targetscip, SCIPgetOrigObjoffset(sourcescip)) );

   return SCIP_OKAY;
}

/** enables constraint compression.
 *
 *  If constraint compression is enabled, fixed variables will be treated as constants
 *  by all constraints that are copied after calling this method.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPenableConsCompression(
   SCIP*                 scip                /**< source SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->origprob != NULL);

   /* check stage */
   SCIP_CALL( SCIPcheckStage(scip, "SCIPenableConsCompression", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* enable problem compression */
   SCIPprobEnableConsCompression(scip->origprob);

   return SCIP_OKAY;
}

/** is constraint compression enabled?
 *
 *  If constraint compression is enabled, fixed variables can be treated as constants
 *  by all constraints that are copied after calling this method.
 *
 *  @return TRUE if problem constraint compression is enabled, otherwise FALSE
 *
 *  @pre This method can be called if scip is in one of the following stages:
  *      - \ref SCIP_STAGE_PROBLEM
  *      - \ref SCIP_STAGE_TRANSFORMING
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPisConsCompressionEnabled(
   SCIP*                 scip                /**< source SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->origprob != NULL);

   /* check stage */
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisConsCompressionEnabled", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   /* is problem compression enabled */
   return SCIPprobIsConsCompressionEnabled(scip->origprob) && SCIPgetStage(scip) == SCIP_STAGE_PROBLEM;
}

/** returns copy of the source variable; if there already is a copy of the source variable in the variable hash map,
 *  it is just returned as target variable; otherwise, if the variables it not marked as relaxation-only, a new variable
 *  will be created and added to the target SCIP; this created variable is added to the variable hash map and returned as target variable;
 *  relaxation-only variables are not copied and FALSE is returned in *success
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *  @note if a new variable was created, this variable will be added to the target-SCIP, but it is not captured
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note targetscip stage does not get changed
 *
 *  @note sourcescip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetVarCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_VAR*             sourcevar,          /**< source variable */
   SCIP_VAR**            targetvar,          /**< pointer to store the target variable */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< should global or local bounds be used? */
   SCIP_Bool*            success             /**< pointer to store whether the copying was successful or not */
   )
{
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_VAR* var;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourcevar != NULL);
   assert(targetvar != NULL);
   assert(sourcevar->scip == sourcescip);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPgetVarCopy", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPgetVarCopy", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);
   *success = TRUE;

   /* try to retrieve copied variable from hashmap */
   if( !uselocalvarmap )
   {
      *targetvar = (SCIP_VAR*) SCIPhashmapGetImage(varmap, sourcevar);
      if( *targetvar != NULL )
         return SCIP_OKAY;
   }

   /* reject copying of relaxation-only variables */
   if( SCIPvarIsRelaxationOnly(sourcevar) )
   {
      *success = FALSE;
      *targetvar = NULL;
   }

   /* if the target SCIP is already in solving stage we currently are not copying the variable!
    * this has to be done because we cannot simply add variables to SCIP during solving and thereby enlarge the search
    * space.
    * unlike column generation we cannot assume here that the variable could be implicitly set to zero in all prior
    * computations
    */
   if( SCIPgetStage(targetscip) > SCIP_STAGE_PROBLEM )
   {
      *success = FALSE;
      *targetvar = NULL;

      return SCIP_OKAY;
   }

   /* create the variable mapping hash map */
   if( uselocalvarmap )
   {
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPgetNVars(sourcescip)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPgetNConss(sourcescip)) );
   }
   else
      localconsmap = consmap;

   /* if variable does not exist yet in target SCIP, create it */
   switch( SCIPvarGetStatus(sourcevar) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      SCIP_CALL( SCIPvarCopy(&var, targetscip->mem->probmem, targetscip->set, targetscip->stat,
            sourcescip, sourcevar, localvarmap, localconsmap, global) );
      break;

   case SCIP_VARSTATUS_AGGREGATED:
   {
      SCIP_CONS* cons;
      char name[SCIP_MAXSTRLEN];

      SCIP_VAR* sourceaggrvar;
      SCIP_VAR* targetaggrvar;
      SCIP_Real aggrcoef;
      SCIP_Real constant;

      /* get aggregation data */
      sourceaggrvar = SCIPvarGetAggrVar(sourcevar);
      aggrcoef = SCIPvarGetAggrScalar(sourcevar);
      constant = SCIPvarGetAggrConstant(sourcevar);

      /* get copy of the aggregation variable */
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourceaggrvar, &targetaggrvar, localvarmap, localconsmap, global, success) );
      assert(*success);

      /* create copy of the aggregated variable */
      SCIP_CALL( SCIPvarCopy(&var, targetscip->mem->probmem, targetscip->set, targetscip->stat,
            sourcescip, sourcevar, localvarmap, localconsmap, global) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_aggr", SCIPvarGetName(sourcevar));

      /* add aggregation x = a*y + c as linear constraint x - a*y = c */
      SCIP_CALL( SCIPcreateConsLinear(targetscip, &cons, name, 0, NULL, NULL, constant,
            constant, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCoefLinear(targetscip, cons, var, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(targetscip, cons, targetaggrvar, -aggrcoef) );

      SCIP_CALL( SCIPaddCons(targetscip, cons) );
      SCIP_CALL( SCIPreleaseCons(targetscip, &cons) );

      break;
   }
   case SCIP_VARSTATUS_MULTAGGR:
   {
      SCIP_CONS* cons;
      char name[SCIP_MAXSTRLEN];

      SCIP_VAR** sourceaggrvars;
      SCIP_VAR** targetaggrvars;
      SCIP_Real* aggrcoefs;
      SCIP_Real constant;

      int naggrvars;
      int i;

      /* get the active representation */
      SCIP_CALL( SCIPflattenVarAggregationGraph(sourcescip, sourcevar) );

      /* get multi-aggregation data */
      naggrvars = SCIPvarGetMultaggrNVars(sourcevar);
      sourceaggrvars = SCIPvarGetMultaggrVars(sourcevar);
      aggrcoefs = SCIPvarGetMultaggrScalars(sourcevar);
      constant = SCIPvarGetMultaggrConstant(sourcevar);

      SCIP_CALL( SCIPallocBufferArray(targetscip, &targetaggrvars, naggrvars) );

      /* get copies of the active variables of the multi-aggregation */
      for( i = 0; i < naggrvars; ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourceaggrvars[i], &targetaggrvars[i], localvarmap, localconsmap, global, success) );
         assert(*success);
      }

      /* create copy of the multi-aggregated variable */
      SCIP_CALL( SCIPvarCopy(&var, targetscip->mem->probmem, targetscip->set, targetscip->stat,
            sourcescip, sourcevar, localvarmap, localconsmap, global) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_multaggr", SCIPvarGetName(sourcevar));

      /* add multi-aggregation x = a^T y + c as linear constraint a^T y - x = -c */
      SCIP_CALL( SCIPcreateConsLinear(targetscip, &cons, name, naggrvars, targetaggrvars, aggrcoefs, -constant,
            -constant, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCoefLinear(targetscip, cons, var, -1.0) );
      SCIP_CALL( SCIPaddCons(targetscip, cons) );
      SCIP_CALL( SCIPreleaseCons(targetscip, &cons) );

      SCIPfreeBufferArray(targetscip, &targetaggrvars);

      break;
   }
   case SCIP_VARSTATUS_NEGATED:
   {
      SCIP_VAR* sourcenegatedvar;
      SCIP_VAR* targetnegatedvar;

      /* get negated source variable */
      sourcenegatedvar = SCIPvarGetNegationVar(sourcevar);
      assert(sourcenegatedvar != NULL);
      assert(SCIPvarGetStatus(sourcenegatedvar) != SCIP_VARSTATUS_NEGATED);

      /* get copy of negated source variable */
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourcenegatedvar, &targetnegatedvar, localvarmap, localconsmap, global, success) );
      assert(*success);
      assert(SCIPvarGetStatus(targetnegatedvar) != SCIP_VARSTATUS_NEGATED);

      /* get negation of copied negated source variable, this is the target variable */
      SCIP_CALL( SCIPgetNegatedVar(targetscip, targetnegatedvar, targetvar) );
      assert(SCIPvarGetStatus(*targetvar) == SCIP_VARSTATUS_NEGATED);

      /* free local hash maps if necessary */
      if( uselocalvarmap )
         SCIPhashmapFree(&localvarmap);

      if( uselocalconsmap )
         SCIPhashmapFree(&localconsmap);

      /* we have to return right away, to avoid adding the negated variable to the problem since the "not negated"
       * variable was already added */
      return SCIP_OKAY;
   }
   default:
      /* note that this is in an internal SCIP error since the variable status is only handled by the core */
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_ERROR; /*lint !e527*/
   }

   /* add the (new) target variable to the target problem */
   SCIP_CALL( SCIPaddVar(targetscip, var) );

   *targetvar = var;

   /* remove the variable capture which was done due to the creation of the variable */
   SCIP_CALL( SCIPreleaseVar(targetscip, &var) );

   /* free local hash maps if necessary */
   if( uselocalvarmap )
      SCIPhashmapFree(&localvarmap);

   if( uselocalconsmap )
      SCIPhashmapFree(&localconsmap);

   return SCIP_OKAY;
}

/** copies all original or active variables from source-SCIP except those that are marked as relaxation-only, fixed, or aggregated
 *  and adds these variable to the target-SCIP
 *
 *  the mapping between these variables are stored in the variable hashmap
 *  target-SCIP has to be in problem creation stage
 */
static
SCIP_RETCODE copyVars(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             original,           /**< should original variables be copied? */
   SCIP_Bool             global              /**< should global or local bounds be used? (for original=FALSE) */
   )
{
   SCIP_VAR** sourcevars;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   int nsourcevars;
#ifndef NDEBUG
   int nrelaxonlybinvars = 0;
   int nrelaxonlyintvars = 0;
   int nrelaxonlyimplvars = 0;
   int nrelaxonlycontvars = 0;
#endif
   int i;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(nfixedvars == 0 || fixedvars != NULL);
   assert(nfixedvars == 0 || fixedvals != NULL);

   if( original )
   {
      /* get original variables of the source SCIP */
      SCIP_CALL( SCIPgetOrigVarsData(sourcescip, &sourcevars, &nsourcevars, NULL, NULL, NULL, NULL) );
   }
   else
   {
      /* get active variables of the source SCIP */
      SCIP_CALL( SCIPgetVarsData(sourcescip, &sourcevars, &nsourcevars, NULL, NULL, NULL, NULL) );
   }

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPgetNVars(sourcescip)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPgetNConss(sourcescip)) );
   }
   else
      localconsmap = consmap;

   /* create the variables of the target SCIP */
   for( i = 0; i < nsourcevars; ++i )
   {
      SCIP_Bool success;
      SCIP_VAR* targetvar;

      if( SCIPvarIsRelaxationOnly(sourcevars[i]) )
      {
#ifndef NDEBUG
         switch( SCIPvarGetType(sourcevars[i]) )
         {
         case SCIP_VARTYPE_BINARY:
            nrelaxonlybinvars++;
            break;
         case SCIP_VARTYPE_INTEGER:
            nrelaxonlyintvars++;
            break;
         case SCIP_VARTYPE_IMPLINT:
            nrelaxonlyimplvars++;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            nrelaxonlycontvars++;
            break;
         default:
            SCIPerrorMessage("unknown variable type\n");
            return SCIP_INVALIDDATA;
         }
#endif
         continue;
      }

      /* copy variable and add this copy to the target SCIP if the copying was valid */
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourcevars[i], &targetvar, localvarmap, localconsmap, global, &success) );
      assert(success);
      assert(targetvar != NULL);
   }

   /* fix the variables that should be fixed right away */
   for( i = 0; i < nfixedvars; ++i )
   {
      SCIP_VAR* targetvar;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      if( SCIPvarIsRelaxationOnly(sourcevars[i]) )
         continue;

      /* retrieve target variable as image of the source variable */
      targetvar = (SCIP_VAR*) SCIPhashmapGetImage(localvarmap, (void *)fixedvars[i]);
      assert(targetvar != NULL);

      /* fix the variable to the specified value */
      infeasible = fixed = FALSE;
      SCIP_CALL( SCIPfixVar(targetscip, targetvar, fixedvals[i], &infeasible, &fixed) );

      assert(!infeasible);
      assert(fixed);
   }

   /* integer variables that are fixed to zero or one or have bounds [0,1] will be converted to binaries */
#ifndef NDEBUG
   if( original )
   {
      /* TODO : account for integers converted to binaries
      assert(SCIPgetNOrigBinVars(sourcescip) == SCIPgetNOrigBinVars(targetscip));
      assert(SCIPgetNOrigIntVars(sourcescip) == SCIPgetNOrigIntVars(targetscip));
      assert(SCIPgetNOrigImplVars(sourcescip) == SCIPgetNOrigImplVars(targetscip));
      assert(SCIPgetNOrigContVars(sourcescip) == SCIPgetNOrigContVars(targetscip));
      */
   }
   else
   {
      SCIP_VAR** sourcefixedvars;
      int nsourcefixedvars;
      int nfixedbinvars;
      int nfixedintvars;
      int nfixedimplvars;
      int nfixedcontvars;

      sourcefixedvars = SCIPgetFixedVars(sourcescip);
      nsourcefixedvars = SCIPgetNFixedVars(sourcescip);
      nfixedbinvars = 0;
      nfixedintvars = 0;
      nfixedimplvars = 0;
      nfixedcontvars = 0;

      /* count number of fixed variables for all variable types */
      for( i = 0; i < nsourcefixedvars; ++i )
      {
         switch( SCIPvarGetType(sourcefixedvars[i]) )
         {
         case SCIP_VARTYPE_BINARY:
            nfixedbinvars++;
            break;
         case SCIP_VARTYPE_INTEGER:
            nfixedintvars++;
            break;
         case SCIP_VARTYPE_IMPLINT:
            nfixedimplvars++;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            nfixedcontvars++;
            break;
         default:
            SCIPerrorMessage("unknown variable type\n");
            return SCIP_INVALIDDATA;
         }
      }
      assert(nsourcefixedvars == nfixedbinvars + nfixedintvars + nfixedimplvars + nfixedcontvars);
      assert(SCIPgetNBinVars(sourcescip) <= SCIPgetNBinVars(targetscip) + nrelaxonlybinvars);
      assert(SCIPgetNIntVars(sourcescip) + SCIPgetNBinVars(sourcescip) <= SCIPgetNIntVars(targetscip) + nrelaxonlyintvars + SCIPgetNBinVars(targetscip) + nrelaxonlybinvars);
      assert(SCIPgetNIntVars(targetscip) + nrelaxonlyintvars + SCIPgetNBinVars(targetscip) + nrelaxonlybinvars <= SCIPgetNIntVars(sourcescip) + SCIPgetNBinVars(sourcescip) + nfixedbinvars + nfixedintvars );
      assert(SCIPgetNImplVars(sourcescip) <= SCIPgetNImplVars(targetscip) + nrelaxonlyimplvars);
      assert(SCIPgetNImplVars(targetscip) + nrelaxonlyimplvars <= SCIPgetNImplVars(sourcescip) + nfixedimplvars);
      assert(SCIPgetNContVars(sourcescip) <= SCIPgetNContVars(targetscip) + nrelaxonlycontvars);
      assert(SCIPgetNContVars(targetscip) + nrelaxonlycontvars <= SCIPgetNContVars(sourcescip) + nfixedcontvars);
   }
#endif

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}

/** Copies all active (thus unfixed) variables from source-SCIP, except those that are marked as relaxation only,
 *  and adds these variable to the target-SCIP.
 *
 *  The mapping between these variables are stored in the variable hashmap.
 *
 *  The target-SCIP has to be in problem creation stage.
 *
 *  @note the variables are added to the target-SCIP but not captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyVars(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             global              /**< should global or local bounds be used? */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyVars", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( copyVars(sourcescip, targetscip, varmap, consmap, fixedvars, fixedvals, nfixedvars, FALSE, global) );

   return SCIP_OKAY;
}

/** copies all original variables from source-SCIP and adds these variable to the target-SCIP; the mapping between these
 *  variables are stored in the variable hashmap, target-SCIP has to be in problem creation stage, fixed and aggregated
 *  variables do not get copied
 *
 *  @note the variables are added to the target-SCIP but not captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyOrigVars(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars          /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyOrigVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyOrigVars", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( copyVars(sourcescip, targetscip, varmap, consmap, fixedvars, fixedvals, nfixedvars, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** merges the histories of variables from a source SCIP into a target SCIP. The two data structures should point to
 *  different SCIP instances.
 *
 *  @note the notion of source and target is inverted here; \p sourcescip usually denotes a copied SCIP instance, whereas
 *        \p targetscip denotes the original instance
 */
SCIP_RETCODE SCIPmergeVariableStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_VAR**            sourcevars,         /**< source variables for history merge, NULL entries are ignored */
   SCIP_VAR**            targetvars,         /**< target variables for history merge, NULL entries are ignored */
   int                   nvars               /**< number of variables in both variable arrays */
   )
{
   int i;

   /* check if target scip has been set to allow merging variable statistics */
   if( !targetscip->set->history_allowmerge )
      return SCIP_OKAY;

   assert(nvars == 0 || (sourcevars != NULL && targetvars != NULL));
   assert(sourcescip != targetscip);

   /* we do not want to copy statistics from a scip that has not really started solving */
   if( SCIPgetStage(sourcescip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* if the transformation of the source was subject to scaling, the history information cannot be just copied */
   if( !SCIPsetIsEQ(targetscip->set, 1.0, SCIPgetOrigObjscale(sourcescip))
      || !SCIPsetIsEQ(targetscip->set, 0.0, SCIPgetOrigObjoffset(sourcescip)) )
      return SCIP_OKAY;

   /* merge histories of the targetSCIP-variables to the SCIP variables. */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VARSTATUS sourcevarstatus;

      if( sourcevars[i] == NULL || targetvars[i] == NULL )
         continue;

      assert(sourcevars[i]->scip == sourcescip);
      assert(targetvars[i]->scip == targetscip);

      sourcevarstatus = SCIPvarGetStatus(sourcevars[i]);

      /* depending on the variable status, we use either the transformed variable history or the history of the col itself */
      switch( sourcevarstatus )
      {
      case SCIP_VARSTATUS_ORIGINAL:
         assert(NULL != SCIPvarGetTransVar(sourcevars[i]));
         SCIPvarMergeHistories(targetvars[i], SCIPvarGetTransVar(sourcevars[i]), targetscip->stat);
         break;
      case SCIP_VARSTATUS_COLUMN:
         SCIPvarMergeHistories(targetvars[i], sourcevars[i], targetscip->stat);
         break;
      default:
         /* other variable status are currently not supported for the merging */
         break;
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** merges the statistics of NLPIs from a source SCIP into a target SCIP
 *
 * The two SCIP instances should point to different SCIP instances.
 *
 *  @note the notion of source and target is inverted here; \p sourcescip usually denotes a copied SCIP instance, whereas
 *        \p targetscip denotes the original instance
 */
void SCIPmergeNLPIStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_Bool             reset               /**< whether to reset statistics in sourcescip */
   )
{
   int i;

   assert(sourcescip != targetscip);

   for( i = 0; i < sourcescip->set->nnlpis; ++i )
   {
      SCIP_NLPI* sourcenlpi;
      SCIP_NLPI* targetnlpi;

      sourcenlpi = sourcescip->set->nlpis[i];
      /* probably NLPI is on same position in target and source, otherwise do search */
      if( strcmp(SCIPnlpiGetName(targetscip->set->nlpis[i]), SCIPnlpiGetName(sourcenlpi)) == 0 )
         targetnlpi = targetscip->set->nlpis[i];
      else
         targetnlpi = SCIPsetFindNlpi(targetscip->set, SCIPnlpiGetName(sourcenlpi));

      if( targetnlpi != NULL )
         SCIPnlpiMergeStatistics(targetnlpi, sourcenlpi, reset);
      else
      {
         SCIPdebugMsg(targetscip, "NLPI <%s> from source SCIP not available in target SCIP\n", SCIPnlpiGetName(sourcenlpi));
      }
   }
}

/** provides values of a solution from a subscip according to the variable in the main scip
 *
 * Given a subscip solution, fills an array with solution values, matching the variables given by SCIPgetVars().
 * Variables that are relaxation-only in the master SCIP are set to 0 or the bound closest to 0. Such variables
 * are represented as NULL entry in the \p subvars array.
 */
static
SCIP_RETCODE translateSubSol(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables from the subproblem in the same order as the main \p scip */
   SCIP_Real*            solvals             /**< array where to set values taken from subsol, must have length at least SCIPgetNVars(scip) */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subsol != NULL);
   assert(subvars != NULL);
   assert(solvals != NULL);

   /* copy the solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* copy the solution */
   for( i = 0; i < nvars; ++i )
   {
      if( subvars[i] != NULL )
         solvals[i] = SCIPgetSolVal(subscip, subsol, subvars[i]);
      else
         solvals[i] = MIN(MAX(0.0, SCIPvarGetLbLocal(vars[i])), SCIPvarGetUbLocal(vars[i]));  /*lint !e666*/
   }

   return SCIP_OKAY;
}

/** translates a solution from a subscip to the main scip
 *
 * Variables that are relaxation-only in the master SCIP are set to 0 or the bound closest to 0. Such variables
 * are represented as NULL entry in the \p subvars array.
 *
 * @note This method allocates a new solution of the main \p scip that needs to be freed by the user.
 */
SCIP_RETCODE SCIPtranslateSubSol(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_VAR**            subvars,            /**< the variables from the subproblem in the same order as the main \p scip */
   SCIP_SOL**            newsol              /**< buffer to store pointer to created solution in main SCIP */
   )
{
   SCIP_VAR** vars;
   int nvars;
   SCIP_Real* subsolvals;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subsol != NULL);
   assert(subvars != NULL);
   assert(newsol != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* get the solution values */
   SCIP_CALL( translateSubSol(scip, subscip, subsol, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, *newsol, nvars, vars, subsolvals) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** checks the solutions from the subscip and adds the first one that is found feasible to the master SCIP
 *
 * Variables that are relaxation-only in the master SCIP are set to 0 or the bound closest to 0. Such variables
 * are represented as NULL entry in the \p subvars array.
 */
SCIP_RETCODE SCIPtranslateSubSols(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_VAR**            subvars,            /**< the variables from the subproblem in the same order as the main \p scip */
   SCIP_Bool*            success,            /**< pointer to store, whether new solution was found */
   int*                  solindex            /**< pointer to store solution index of stored solution, or NULL if not of interest */
   )
{
   SCIP_SOL* newsol = NULL;
   SCIP_SOL** subsols;
   int nsubsols;
   int i;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Real* solvals;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);
   assert(subvars != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check, whether a solution was found */
   if( SCIPgetNSols(subscip) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   for( i = 0; i < nsubsols; ++i )
   {
      /* better do not copy unbounded solutions as this will mess up the SCIP solution status */
      if( SCIPisInfinity(scip, -SCIPgetSolOrigObj(subscip, subsols[i])) )
         continue;

      if( newsol == NULL )
      {
         SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
         if( solindex != NULL )
            *solindex = SCIPsolGetIndex(newsol);
      }

      /* put values from subsol into newsol */
      SCIP_CALL( translateSubSol(scip, subscip, subsols[i], subvars, solvals) );
      SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, solvals) );

      /* check whether feasible */
      SCIP_CALL( SCIPcheckSol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );
      if( *success )
      {
         /* if feasible, then there is a good chance that we can add it
          * we use SCIPaddSolFree to make sure that newsol is indeed added and not some copy, so *solindex stays valid
          */
         SCIP_CALL( SCIPaddSolFree(scip, &newsol, success) );
         if( *success )
         {
            SCIPdebugMsg(scip, "-> accepted solution of value %g\n", SCIPgetSolOrigObj(subscip, subsols[i]));
            break;
         }
         else
         {
            /* continue with next subsol
             * as we have used addSolFree, newsol should be NULL now
             */
            assert(newsol == NULL);
         }
      }
   }

   SCIPfreeBufferArray(scip, &solvals);

   if( newsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
   }

   return SCIP_OKAY;
}

/** returns copy of the source constraint; if there already is a copy of the source constraint in the constraint hash
 *  map, it is just returned as target constraint; elsewise a new constraint will be created; this created constraint is
 *  added to the constraint hash map and returned as target constraint; the variable map is used to map the variables of
 *  the source SCIP to the variables of the target SCIP
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 *
 *  @note The constraint is not added to the target SCIP. You can check whether a constraint is added by calling
 *        SCIPconsIsAdded(). (If you mix SCIPgetConsCopy() with SCIPcopyConss() you should pay attention to what you add
 *        explicitly and what is already added.)
 *
 *  @note The constraint is always captured, either during the creation of the copy or after finding the copy of the
 *        constraint in the constraint hash map
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetConsCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_CONS*            sourcecons,         /**< source constraint of the source SCIP */
   SCIP_CONS**           targetcons,         /**< pointer to store the created target constraint */
   SCIP_CONSHDLR*        sourceconshdlr,     /**< source constraint handler for this constraint */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           name,               /**< name of constraint, or NULL if the name of the source constraint should be used */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid or not */
   )
{
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;

   assert(targetcons != NULL);
   assert(sourceconshdlr != NULL);

   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPgetConsCopy", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPgetConsCopy", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   /* a variables map and a constraint map is needed to avoid infinite recursion */
   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPgetNVars(sourcescip)) );
   }
   else
      localvarmap = varmap;

   *targetcons = NULL;
   if( uselocalconsmap )
   {
      /* create local constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPgetNConss(sourcescip)) );
   }
   else
   {
      /* use global map and try to retrieve copied constraint */
      localconsmap = consmap;
      *targetcons = (SCIP_CONS*) SCIPhashmapGetImage(localconsmap, sourcecons);
   }

   if( *targetcons != NULL )
   {
      /* if found capture existing copy of the constraint */
      SCIP_CALL( SCIPcaptureCons(targetscip, *targetcons) );
      *valid = TRUE;
   }
   else
   {
      /* otherwise create a copy of the constraint */
      SCIP_CALL( SCIPconsCopy(targetcons, targetscip->set, name, sourcescip, sourceconshdlr, sourcecons, localvarmap, localconsmap,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );

      /* it is possible for the constraint handler to declare the copy valid although no target constraint was created */
      assert(*targetcons == NULL || *valid);

      /* if a target constraint was created */
      if( *targetcons != NULL && !uselocalconsmap )
      {
         /* insert constraint into mapping between source SCIP and the target SCIP */
         SCIP_CALL( SCIPhashmapInsert(consmap, sourcecons, *targetcons) );
      }
   }

   /* free locally allocated hash maps */
   if( uselocalvarmap )
   {
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}

/** copies constraints from the source-SCIP and adds these to the target-SCIP; for mapping the
 *  variables between the source and the target SCIP a hash map can be given; if the variable hash
 *  map is NULL or necessary variable mapping is missing, the required variables are created in the
 *  target-SCIP and added to the hash map, if not NULL; all variables which are created are added to
 *  the target-SCIP but not (user) captured; if the constraint hash map is not NULL the mapping
 *  between the constraints of the source and target-SCIP is stored
 *
 *  *valid is set to TRUE iff all constraints that are marked as checked or enforced were copied successfully.
 *  If other constraints could not be copied, *valid can still be set to TRUE.
 *
 *  @note the constraints are added to the target-SCIP but are not (user) captured in the target SCIP. (If you mix
 *        SCIPgetConsCopy() with SCIPcopyConss() you should pay attention to what you add explicitly and what is already
 *        added.) You can check whether a constraint is added by calling SCIPconsIsAdded().
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyConss(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance?
                                              *   If TRUE, the modifiable flag of constraints will be copied. */
   SCIP_Bool*            valid               /**< pointer to store whether all checked or enforced constraints were validly copied */
   )
{
   SCIP_CONSHDLR** sourceconshdlrs;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   int nsourceconshdlrs;
   int i;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(valid      != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyConss", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check if we locally need to create a variable or constraint hash map */
   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPgetNVars(sourcescip)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPgetNConss(sourcescip)) );
   }
   else
      localconsmap = consmap;

   nsourceconshdlrs = SCIPgetNConshdlrs(sourcescip);
   sourceconshdlrs = SCIPgetConshdlrs(sourcescip);
   assert(nsourceconshdlrs == 0 || sourceconshdlrs != NULL);

   *valid = TRUE;

   /* copy constraints: loop through all (source) constraint handlers */
   for( i = 0; i < nsourceconshdlrs; ++i )
   {
      SCIP_CONS** sourceconss;
      SCIP_CONS* targetcons;
      int nsourceconss;
      int c;

      assert(sourceconshdlrs[i] != NULL);

      /* constraint handlers have to explicitly set the valid pointer to TRUE for every single constraint */

      /* Get all active constraints for copying; this array contains all active constraints;
       * constraints are active if they are globally valid and not deleted after presolving OR they
       * were locally added during the search and we are currently in a node which belongs to the
       * corresponding subtree.
       */
      nsourceconss = SCIPconshdlrGetNActiveConss(sourceconshdlrs[i]);
      sourceconss = SCIPconshdlrGetConss(sourceconshdlrs[i]);

#ifdef SCIP_DISABLED_CODE
      /* @todo using the following might reduce the number of copied constraints - check whether this is better */
      /* Get all checked constraints for copying; this included local constraints */
      if( !global )
      {
         nsourceconss = SCIPconshdlrGetNCheckConss(sourceconshdlrs[i]);
         sourceconss = SCIPconshdlrGetCheckConss(sourceconshdlrs[i]);
      }
#endif

      assert(nsourceconss == 0 || sourceconss != NULL);

      if( nsourceconss > 0 )
      {
         SCIPdebugMsg(sourcescip, "Attempting to copy %d %s constraints\n", nsourceconss, SCIPconshdlrGetName(sourceconshdlrs[i]));
      }

      /* copy all constraints of one constraint handler */
      for( c = 0; c < nsourceconss; ++c )
      {
         SCIP_Bool singlevalid = FALSE;
         /* all constraints have to be active */
         assert(sourceconss[c] != NULL);
         assert(SCIPconsIsActive(sourceconss[c]));
         assert(!SCIPconsIsDeleted(sourceconss[c]));

         /* in case of copying the global problem we have to ignore the local constraints which are active */
         if( global && SCIPconsIsLocal(sourceconss[c]) )
         {
            SCIPdebugMsg(sourcescip, "did not copy local constraint <%s> when creating global copy\n", SCIPconsGetName(sourceconss[c]));
            continue;
         }

         /* use the copy constructor of the constraint handler and creates and captures the constraint if possible */
         targetcons = NULL;
         SCIP_CALL( SCIPgetConsCopy(sourcescip, targetscip, sourceconss[c], &targetcons, sourceconshdlrs[i], localvarmap, localconsmap, NULL,
               SCIPconsIsInitial(sourceconss[c]), SCIPconsIsSeparated(sourceconss[c]),
               SCIPconsIsEnforced(sourceconss[c]), SCIPconsIsChecked(sourceconss[c]),
               SCIPconsIsPropagated(sourceconss[c]), FALSE, SCIPconsIsModifiable(sourceconss[c]),
               SCIPconsIsDynamic(sourceconss[c]), SCIPconsIsRemovable(sourceconss[c]), FALSE, global, &singlevalid) );

         /* it is possible for a constraint handler to declare the copy valid, although no target constraint was created */
         assert(targetcons == NULL || singlevalid);

         /* add the copied constraint to target SCIP if the copying process created a constraint */
         if( targetcons != NULL )
         {
            if( !enablepricing )
               SCIPconsSetModifiable(targetcons, FALSE);

            /* add constraint to target SCIP */
            SCIP_CALL( SCIPaddCons(targetscip, targetcons) );

            /* add the conflict constraint to the store of targetscip */
            if( SCIPconsIsConflict(sourceconss[c]) )
            {
               /* add the constraint as a conflict to the conflict pool of targetscip */
               SCIP_CALL( SCIPconflictstoreAddConflict(targetscip->conflictstore, targetscip->mem->probmem, targetscip->set,
                     targetscip->stat, NULL, NULL, targetscip->reopt, targetcons, SCIP_CONFTYPE_UNKNOWN, FALSE, -SCIPinfinity(targetscip)) );
            }

            /* release constraint once for the creation capture */
            SCIP_CALL( SCIPreleaseCons(targetscip, &targetcons) );
         }
         else
         {
            /* if an enforced or checked constraint could not be copied, then the copy is not valid, i.e.,
             * the feasible set may be larger; for other constraints, it should be safe if they are omitted
             * from the copy
             */
            if( SCIPconsIsEnforced(sourceconss[c]) || SCIPconsIsChecked(sourceconss[c]) )
               *valid = FALSE;
            SCIPdebugMsg(sourcescip, "Constraint %s not copied, copy is %svalid\n",
                  SCIPconsGetName(sourceconss[c]), *valid ? "" : "not ");
         }
      }
   }

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}

/** copies all original constraints from the source-SCIP and adds these to the target-SCIP; for mapping the
 *  variables between the source and the target SCIP a hash map can be given; if the variable hash
 *  map is NULL or necessary variable mapping is missing, the required variables are created in the
 *  target-SCIP and added to the hash map, if not NULL; all variables which are created are added to
 *  the target-SCIP but not (user) captured; if the constraint hash map is not NULL the mapping
 *  between the constraints of the source and target-SCIP is stored
 *
 *  @note the constraints are added to the target-SCIP but are not (user) captured in the target SCIP. (If you mix
 *        SCIPgetConsCopy() with SCIPcopyConss() you should pay attention to what you add explicitly and what is already
 *        added.) You can check whether a constraint is added by calling SCIPconsIsAdded().
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyOrigConss(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance?
                                              *   If TRUE, the modifiable flag of constraints will be copied. */
   SCIP_Bool*            valid               /**< pointer to store whether all constraints were validly copied */
   )
{
   SCIP_CONS** sourceconss;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   int nsourceconss;
   int c;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(valid      != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyOrigConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyOrigConss", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check if we locally need to create a variable or constraint hash map */
   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPgetNVars(sourcescip)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPgetNConss(sourcescip)) );
   }
   else
      localconsmap = consmap;

   sourceconss = SCIPgetOrigConss(sourcescip);
   nsourceconss = SCIPgetNOrigConss(sourcescip);

   *valid = TRUE;

   SCIPdebugMsg(sourcescip, "Attempting to copy %d original constraints\n", nsourceconss);

   /* copy constraints: loop through all (source) constraint handlers */
   for( c = 0; c < nsourceconss; ++c )
   {
      SCIP_CONS* targetcons;
      SCIP_Bool success;

      /* constraint handlers have to explicitly set the success pointer to TRUE */
      success = FALSE;

      /* all constraints have to be active */
      assert(sourceconss[c] != NULL);
      assert(SCIPconsIsOriginal(sourceconss[c]));

      /* use the copy constructor of the constraint handler and creates and captures the constraint if possible */
      targetcons = NULL;
      SCIP_CALL( SCIPgetConsCopy(sourcescip, targetscip, sourceconss[c], &targetcons, SCIPconsGetHdlr(sourceconss[c]), localvarmap, localconsmap, NULL,
            SCIPconsIsInitial(sourceconss[c]), SCIPconsIsSeparated(sourceconss[c]),
            SCIPconsIsEnforced(sourceconss[c]), SCIPconsIsChecked(sourceconss[c]),
            SCIPconsIsPropagated(sourceconss[c]), FALSE, SCIPconsIsModifiable(sourceconss[c]),
            SCIPconsIsDynamic(sourceconss[c]), SCIPconsIsRemovable(sourceconss[c]), FALSE, TRUE, &success) );

      /* add the copied constraint to target SCIP if the copying process was valid */
      if( success )
      {
         assert(targetcons != NULL);

         if( !enablepricing )
            SCIPconsSetModifiable(targetcons, FALSE);

         /* add constraint to target SCIP */
         SCIP_CALL( SCIPaddCons(targetscip, targetcons) );

         /* release constraint once for the creation capture */
         SCIP_CALL( SCIPreleaseCons(targetscip, &targetcons) );
      }
      else
      {
         *valid = FALSE;
         SCIPdebugMsg(sourcescip, "failed to copy constraint %s\n", SCIPconsGetName(sourceconss[c]));
      }
   }

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}


/** convert all active cuts from cutpool to linear constraints
 *
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note SCIP stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPconvertCutsToConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   int*                  ncutsadded          /**< pointer to store number of added cuts, or NULL */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* check stages for the SCIP data structure */
   SCIP_CALL( SCIPcheckStage(scip, "SCIPconvertCutsToConss", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   /* if we do not have any cuts, nothing can be converted */
   if( scip->set->stage < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* create out of all active cuts in cutpool linear constraints in targetscip */
   SCIP_CALL( SCIPcopyCuts(scip, scip, varmap, consmap, global, ncutsadded) );

   return SCIP_OKAY;
}

/** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip
 *
 *  Cuts that contain variables that are marked as relaxation-only are skipped.
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyCuts(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   int*                  ncutsadded          /**< pointer to store number of copied cuts, or NULL */
   )
{
   SCIP_CUT** cuts;
   int ncuts;
   int nlocalcutsadded;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyCuts", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyCuts", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   if ( ncutsadded != NULL )
      *ncutsadded = 0;
   nlocalcutsadded = 0;

   /* if we do not have any cuts, nothing can be converted */
   if( sourcescip->set->stage < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   if( SCIPfindConshdlr(targetscip, "linear") == NULL )
   {
      SCIPdebugMsg(sourcescip, "No linear constraint handler available. Cannot convert cuts.\n");
      return SCIP_OKAY;
   }

   /* convert cut from global cut pool */
   cuts = SCIPgetPoolCuts(sourcescip);
   ncuts = SCIPgetNPoolCuts(sourcescip);

   SCIP_CALL( copyCuts(sourcescip, targetscip, cuts, ncuts, varmap, consmap, global, &nlocalcutsadded) );

   SCIPdebugMsg(sourcescip, "Converted %d active cuts to constraints.\n", nlocalcutsadded);

   /* convert delayed cuts from global delayed cut pool */
   cuts = SCIPgetDelayedPoolCuts(sourcescip);
   ncuts = SCIPgetNDelayedPoolCuts(sourcescip);

   SCIP_CALL( copyCuts(sourcescip, targetscip, cuts, ncuts, varmap, consmap, global, &nlocalcutsadded) );

   if( ncutsadded != NULL )
      *ncutsadded = nlocalcutsadded;

   SCIPdebugMsg(sourcescip, "Converted %d active cuts to constraints.\n", nlocalcutsadded);

   return SCIP_OKAY;
}

/** copies all active conflicts from the conflict pool of sourcescip and adds them as linear constraints to targetscip
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note sourcescip stage does not change
 *
 *  @note targetscip stage does not change
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyConflicts(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance?
                                              *   If TRUE, the modifiable flag of constraints will be copied. */
   SCIP_Bool*            valid               /**< pointer to store whether all constraints were validly copied */
   )
{
   SCIP_CONS** sourceconfs;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   SCIP_Bool success;
   int sourceconfssize;
   int nsourceconfs;
   int c;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyConss", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check if we locally need to create a variable or constraint hash map */
   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPgetNVars(sourcescip)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPgetNConss(sourcescip)) );
   }
   else
      localconsmap = consmap;

   /* get number of conflicts stored in the conflict pool */
   sourceconfssize = SCIPconflictstoreGetNConflictsInStore(sourcescip->conflictstore);

   /* allocate buffer */
   SCIP_CALL( SCIPallocBufferArray(sourcescip, &sourceconfs, sourceconfssize) );

   /* get all conflicts stored in the conflict pool */
   SCIP_CALL( SCIPconflictstoreGetConflicts(sourcescip->conflictstore, sourceconfs, sourceconfssize, &nsourceconfs) );
   assert(nsourceconfs <= sourceconfssize);

   /* copy conflicts */
   for( c = 0; c < nsourceconfs; ++c )
   {
      SCIP_CONS* targetcons;

      /* all constraints have to be active */
      assert(sourceconfs[c] != NULL);
      assert(SCIPconsIsActive(sourceconfs[c]));
      assert(!SCIPconsIsDeleted(sourceconfs[c]));
      assert(SCIPconsIsConflict(sourceconfs[c]));

      /* in case of copying the global problem we have to ignore the local constraints which are active */
      if( global && SCIPconsIsLocal(sourceconfs[c]) )
      {
         SCIPdebugMsg(sourcescip, "did not copy local constraint <%s> when creating global copy\n", SCIPconsGetName(sourceconfs[c]));
         continue;
      }

      /* use the copy constructor of the constraint handler and creates and captures the constraint if possible */
      targetcons = NULL;
      SCIP_CALL( SCIPgetConsCopy(sourcescip, targetscip, sourceconfs[c], &targetcons, SCIPconsGetHdlr(sourceconfs[c]),
            localvarmap, localconsmap, NULL, SCIPconsIsInitial(sourceconfs[c]), SCIPconsIsSeparated(sourceconfs[c]),
            SCIPconsIsEnforced(sourceconfs[c]), SCIPconsIsChecked(sourceconfs[c]),
            SCIPconsIsPropagated(sourceconfs[c]), FALSE, SCIPconsIsModifiable(sourceconfs[c]),
            SCIPconsIsDynamic(sourceconfs[c]), SCIPconsIsRemovable(sourceconfs[c]), FALSE, global, &success) );

      /* add the copied constraint to target SCIP if the copying process was valid */
      if( success )
      {
         assert(targetcons != NULL);

         if( !enablepricing )
            SCIPconsSetModifiable(targetcons, FALSE);

         /* add constraint to target SCIP */
         SCIP_CALL( SCIPaddCons(targetscip, targetcons) );

         /* release constraint once for the creation capture */
         SCIP_CALL( SCIPreleaseCons(targetscip, &targetcons) );
      }
      else
      {
         *valid = FALSE;
         SCIPdebugMsg(sourcescip, "failed to copy constraint %s\n", SCIPconsGetName(sourceconfs[c]));
      }
   }

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}

/** copies implications and cliques of sourcescip to targetscip
 *
 *  This function should be called for a targetscip in transformed stage. It can save time in presolving of the
 *  targetscip, since implications and cliques are copied.
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyImplicationsCliques(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs,            /**< pointer to store the number of performed bound changes, or NULL */
   int*                  ncopied             /**< pointer to store number of copied implications and cliques, or NULL */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_VAR** sourcevars;
   SCIP_Bool success;
   int nvars;
   int nbinvars;
   int ncliques;
   int j;
   int c;

   assert( sourcescip != NULL );
   assert( targetscip != NULL );
   assert( sourcescip != targetscip );
   assert( infeasible != NULL );

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyImplicationsCliques", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyImplicationsCliques", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if ( ncopied != NULL )
      *ncopied = 0;
   if ( nbdchgs != NULL )
      *nbdchgs = 0;

   /* get all active variables */
   SCIP_CALL( SCIPgetVarsData(sourcescip, &sourcevars, &nvars, &nbinvars, NULL, NULL, NULL) );

   /* stop if no possible variables for cliques exist */
   if ( nbinvars == 0 )
      return SCIP_OKAY;

   /* get cliques */
   ncliques = SCIPgetNCliques(sourcescip);
   if ( ncliques > 0 )
   {
      SCIP_VAR** targetclique;

      /* get space for target cliques */
      SCIP_CALL( SCIPallocBufferArray(targetscip, &targetclique, nvars) );
      cliques = SCIPgetCliques(sourcescip);

      /* loop through all cliques */
      for (c = 0; c < ncliques; ++c)
      {
         SCIP_VAR** cliquevars;
         SCIP_Bool* cliquevals;
         int cliquesize;
         int nboundchg = 0;

         assert( cliques[c] != NULL );
         cliquevals = SCIPcliqueGetValues(cliques[c]);
         cliquevars = SCIPcliqueGetVars(cliques[c]);
         cliquesize = SCIPcliqueGetNVars(cliques[c]);

         /* get target variables of clique */
         for (j = 0; j < cliquesize; ++j)
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, cliquevars[j], &targetclique[j], varmap, consmap, global, &success) );
            if ( ! success )
            {
               SCIPdebugMsg(sourcescip, "Getting copy for variable <%s> failed.\n", SCIPvarGetName(cliquevars[j]));
               SCIPfreeBufferArray(targetscip, &targetclique);
               return SCIP_OKAY;
            }
         }

         /* create clique */
         SCIP_CALL( SCIPaddClique(targetscip, targetclique, cliquevals, cliquesize, SCIPcliqueIsEquation(cliques[c]),
               infeasible, &nboundchg) );

         if ( *infeasible )
         {
            SCIPfreeBufferArray(targetscip, &targetclique);
            return SCIP_OKAY;
         }
         if ( ncopied != NULL )
            ++(*ncopied);
         if ( nbdchgs != NULL )
            *nbdchgs += nboundchg;
      }
      SCIPfreeBufferArray(targetscip, &targetclique);
   }

   /* create binary implications */
   for (j = 0; j < nbinvars; ++j)
   {
      SCIP_VAR* sourcevar;
      SCIP_VAR* targetvar;
      SCIP_Bool d;

      sourcevar = sourcevars[j];
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourcevar, &targetvar, varmap, consmap, global, &success) );
      if ( ! success )
      {
         SCIPdebugMsg(sourcescip, "Getting copy for variable <%s> failed.\n", SCIPvarGetName(sourcevar));
         return SCIP_OKAY;
      }

      /* consider both possible implications */
      for (d = 0; d <= 1; ++d)
      {
         SCIP_BOUNDTYPE* impltypes;
         SCIP_VAR** implvars;
         SCIP_Real* implbounds;
         int nimpls;
         int l;

         nimpls = SCIPvarGetNImpls(sourcevar, d);
         if ( nimpls == 0 )
            continue;

         impltypes = SCIPvarGetImplTypes(sourcevar, d);
         implvars = SCIPvarGetImplVars(sourcevar, d);
         implbounds = SCIPvarGetImplBounds(sourcevar, d);

         /* create implications */
         for (l = 0; l < nimpls; ++l)
         {
            SCIP_VAR* implvar;
            int nboundchg = 0;

            SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, implvars[l], &implvar, varmap, consmap, global, &success) );
            if ( ! success )
            {
               SCIPdebugMsg(sourcescip, "Getting copy for variable <%s> failed.\n", SCIPvarGetName(implvars[l]));
               return SCIP_OKAY;
            }

            SCIP_CALL( SCIPaddVarImplication(targetscip, targetvar, d, implvar, impltypes[l], implbounds[l], infeasible, &nboundchg) );
            if ( *infeasible )
               return SCIP_OKAY;
            if ( ncopied != NULL )
               ++(*ncopied);
            if ( nbdchgs != NULL )
               *nbdchgs += nboundchg;
	 }
      }
   }

   return SCIP_OKAY;
}

/** copies parameter settings from sourcescip to targetscip
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_FREE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyParamSettings(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourcescip->set != NULL);
   assert(targetscip->set != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyParamSettings", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyParamSettings", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPsetCopyParams(sourcescip->set, targetscip->set, targetscip->messagehdlr) );

   return SCIP_OKAY;
}

/** gets depth of current scip instance (increased by each copy call)
 *
 *  @return Depth of subscip of SCIP is returned.
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
 *
 *  @note SCIP stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
int SCIPgetSubscipDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );
   assert( scip->stat != NULL );

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSubscipDepth", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->subscipdepth;
}

/** sets depth of scip instance
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note SCIP stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void SCIPsetSubscipDepth(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   newdepth            /**< new subscip depth */
   )
{
   assert( scip != NULL );
   assert( newdepth > 0 );

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPsetSubscipDepth", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert( scip->stat != NULL );
   scip->stat->subscipdepth = newdepth;
}

/** copies source SCIP data into target SCIP data structure
 *
 * distinguishes between
 * - local and global copies
 * - copies of the original or transformed problem
 *
 * Allows for constraint compression by specifying a number of source variables
 * and values that should be fixed in the copy.
 */
static
SCIP_RETCODE doCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           suffix,             /**< optional suffix for problem name inside the target SCIP */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             useconscompression, /**< should constraint compression be used when constraints are created? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             original,           /**< copy original or transformed problem? if TRUE, a copy using local bounds is not possible */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? If TRUE, pricer
                                              *   plugins will be copied and activated, and the modifiable flag of
                                              *   constraints will be respected. If FALSE, valid will be set to FALSE, when
                                              *   there are pricers present */
   SCIP_Bool             threadsafe,         /**< FALSE, if data can be safely shared between the source and target
                                              *   SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid or not, or NULL */
   )
{
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Real startcopytime;
   SCIP_Real copytime;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   SCIP_Bool consscopyvalid;
   SCIP_Bool benderscopyvalid;
   SCIP_Bool localvalid;
   SCIP_Bool msghdlrquiet;
   char name[SCIP_MAXSTRLEN];

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(suffix != NULL);

   /* copy the original problem if we are in SCIP_STAGE_PROBLEM stage */
   if( SCIPgetStage(sourcescip) == SCIP_STAGE_PROBLEM )
      original = TRUE;

   /* global must be TRUE for the original problem */
   assert(global || !original);

   /* get time before start of copy procedure */
   startcopytime = SCIPclockGetTime(sourcescip->stat->copyclock);

   /* start time measuring */
   SCIPclockStart(sourcescip->stat->copyclock, sourcescip->set);

   /* copy all plugins */
   SCIP_CALL( SCIPcopyPlugins(sourcescip, targetscip, TRUE, enablepricing, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, passmessagehdlr, &localvalid) );

   /* in case there are active pricers and pricing is disabled, targetscip will not be a valid copy of sourcescip */
   if( ! enablepricing && SCIPgetNActivePricers(sourcescip) > 0 )
      localvalid = FALSE;

   SCIPdebugMsg(sourcescip, "Copying plugins was%s valid.\n", localvalid ? "" : " not");

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPgetNVars(sourcescip)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPgetNConss(sourcescip)) );
   }
   else
      localconsmap = consmap;

   /* construct name for the target SCIP using the source problem name and the given suffix string */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s", SCIPgetProbName(sourcescip), suffix);

   /* store the quiet state of the message handler, if existent */
   msghdlrquiet = SCIPmessagehdlrIsQuiet(targetscip->messagehdlr);

   /* explicitly suppress output when copying parameters */
   SCIPsetMessagehdlrQuiet(targetscip, TRUE);

   /* copy all settings */
   SCIP_CALL( SCIPcopyParamSettings(sourcescip, targetscip) );

   /* restore original quiet state */
   SCIPsetMessagehdlrQuiet(targetscip, msghdlrquiet);

   /* create problem in the target SCIP copying the source original or transformed problem data */
   if( original )
   {
      SCIP_CALL( SCIPcopyOrigProb(sourcescip, targetscip, localvarmap, localconsmap, name) );
   }
   else
   {
      SCIP_CALL( SCIPcopyProb(sourcescip, targetscip, localvarmap, localconsmap, global, name) );
   }

   /* copy original or transformed variables and perform fixings if needed */
   SCIP_CALL( copyVars(sourcescip, targetscip, localvarmap, localconsmap, fixedvars, fixedvals, nfixedvars, original, global) );

   /* if fixed variables are directly specified or inferred from local bounds, enable constraint compression */
   if( useconscompression && (nfixedvars > 0 || !global) )
   {
      SCIP_CALL( SCIPenableConsCompression(targetscip) );

      /* domain reductions yield a copy that is no longer guaranteed to be valid */
      localvalid = FALSE;
   }

   /* copy all (original) constraints */
   if( original )
   {
      SCIP_CALL( SCIPcopyOrigConss(sourcescip, targetscip, localvarmap, localconsmap, enablepricing, &consscopyvalid) );
   }
   else
   {
      SCIP_CALL( SCIPcopyConss(sourcescip, targetscip, localvarmap, localconsmap, global, enablepricing, &consscopyvalid) );
   }

   SCIPdebugMsg(sourcescip, "Copying constraints was%s valid.\n", consscopyvalid ? "" : " not");

   localvalid = localvalid && consscopyvalid;

   /* copy the Benders' decomposition plugins explicitly, because it requires the variable mapping hash map */
   SCIP_CALL( SCIPcopyBenders(sourcescip, targetscip, localvarmap, threadsafe, &benderscopyvalid) );

   SCIPdebugMsg(sourcescip, "Copying Benders' decomposition plugins was%s valid.\n", benderscopyvalid ? "" : " not");

   localvalid = localvalid && benderscopyvalid;

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   /* stop time measuring */
   SCIPclockStop(sourcescip->stat->copyclock, sourcescip->set);

   /* get time after copying procedure */
   copytime = SCIPclockGetTime(sourcescip->stat->copyclock) - startcopytime;

   if( copytime > sourcescip->stat->maxcopytime )
      sourcescip->stat->maxcopytime = copytime;
   if( copytime < sourcescip->stat->mincopytime )
      sourcescip->stat->mincopytime = copytime;

   /* increase copy counter */
   ++(sourcescip->stat->ncopies);

   targetscip->concurrent = sourcescip->concurrent;
   SCIP_CALL( SCIPsyncstoreRelease(&targetscip->syncstore) );
   targetscip->syncstore = sourcescip->syncstore;
   SCIP_CALL( SCIPsyncstoreCapture(targetscip->syncstore) );

   /* return the information about a valid copy to the user */
   if( valid != NULL )
      *valid = localvalid;

   return SCIP_OKAY;
}

/** copies source SCIP to target SCIP; the copying process is done in the following order:
 *  1) copy the plugins
 *  2) copy the settings
 *  3) create problem data in target-SCIP and copy the problem data of the source-SCIP
 *  4) copy all active variables except those that are marked as relaxation-only
 *  5) copy all constraints
 *
 *  The source problem depends on the stage of the \p sourcescip - In SCIP_STAGE_PROBLEM, the original problem is copied,
 *  otherwise, the transformed problem is copied. For an explicit copy of the original problem, use SCIPcopyOrig().
 *
 *  @note all variables and constraints which are created in the target-SCIP are not (user) captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *        Also, 'passmessagehdlr' should be set to FALSE.
 *  @note the 'threadsafe' parameter should only be set to TRUE if you are absolutely certain that the source and target
 *        SCIP instances will be solved in parallel. The usual case is to set this to FALSE, since thread safety
 *        typically incurs a performance cost.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           suffix,             /**< optional suffix for problem name inside the target SCIP */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? If TRUE, pricer
                                              *   plugins will be copied and activated, and the modifiable flag of
                                              *   constraints will be respected. If FALSE, valid will be set to FALSE, when
                                              *   there are pricers present */
   SCIP_Bool             threadsafe,         /**< FALSE, if data can be safely shared between the source and target
                                              *   SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   )
{
   SCIP_VAR** fixedvars = NULL;
   SCIP_Real* fixedvals = NULL;
   int nfixedvars = 0;
   SCIP_Bool original = FALSE;
   SCIP_Bool useconscompression = FALSE;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(suffix != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopy", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopy", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   /* copy source SCIP data into target SCIP data structure */
   SCIP_CALL( doCopy(sourcescip, targetscip, varmap, consmap, suffix, fixedvars, fixedvals, nfixedvars,
         useconscompression, global, original, enablepricing, threadsafe, passmessagehdlr, valid) );

   return SCIP_OKAY;
}

/** copies source SCIP to target SCIP but compresses constraints
 *
 *  constraint compression is performed by removing fixed variables immediately
 *  during constraint creation if the involved constraint handlers support
 *  compression
 *
 *  the copying process is done in the following order:
 *  1) copy the plugins
 *  2) copy the settings
 *  3) create problem data in target-SCIP and copy the problem data of the source-SCIP
 *  4) copy all active variables except those that are marked as relaxation-only
 *     a) fix all variable copies specified by \p fixedvars, \p fixedvals, and \p nfixedvars
 *     b) enable constraint compression
 *  5) copy all constraints
 *
 * The source problem depends on the stage of the \p sourcescip - In SCIP_STAGE_PROBLEM, the original problem is copied,
 * otherwise, the transformed problem is copied. For an explicit copy of the original problem, use SCIPcopyOrigConsCompression().
 *
 *  @note: in case that a combination of local bounds and explicit fixing values should be used,
 *         the fixing value of a variable is preferred if local bounds and fixing value disagree.
 *
 *  @note all variables and constraints which are created in the target-SCIP are not (user) captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *        Also, 'passmessagehdlr' should be set to FALSE.
 *  @note the 'threadsafe' parameter should only be set to TRUE if you are absolutely certain that the source and target
 *        SCIP instances will be solved in parallel. The usual case is to set this to FALSE, since thread safety
 *        typically incurs a performance cost.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyConsCompression(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           suffix,             /**< optional suffix for problem name inside the target SCIP */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? If TRUE, pricer
                                              *   plugins will be copied and activated, and the modifiable flag of
                                              *   constraints will be respected. If FALSE, valid will be set to FALSE, when
                                              *   there are pricers present */
   SCIP_Bool             threadsafe,         /**< FALSE, if data can be safely shared between the source and target
                                              *   SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   )
{
   SCIP_Bool original = FALSE;
   SCIP_Bool useconscompression = TRUE;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(suffix != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyConsCompression", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyConsCompression", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   /* copy the source problem data */
   SCIP_CALL( doCopy(sourcescip, targetscip, varmap, consmap, suffix, fixedvars, fixedvals, nfixedvars,
         useconscompression, global, original, enablepricing, threadsafe, passmessagehdlr, valid) );

   return SCIP_OKAY;
}


/** copies source SCIP original problem to target SCIP; the copying process is done in the following order:
 *  1) copy the plugins
 *  2) copy the settings
 *  3) create problem data in target-SCIP and copy the original problem data of the source-SCIP
 *  4) copy all original variables
 *  5) copy all original constraints
 *
 *  @note all variables and constraints which are created in the target-SCIP are not (user) captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *        Also, 'passmessagehdlr' should be set to FALSE.
 *  @note the 'threadsafe' parameter should only be set to TRUE if you are absolutely certain that the source and target
 *        SCIP instances will be solved in parallel. The usual case is to set this to FALSE, since thread safety
 *        typically incurs a performance cost.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyOrig(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           suffix,             /**< suffix which will be added to the names of the target SCIP, might be empty */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? If TRUE, pricer
                                              *   plugins will be copied and activated, and the modifiable flag of
                                              *   constraints will be respected. If FALSE, valid will be set to FALSE, when
                                              *   there are pricers present */
   SCIP_Bool             threadsafe,         /**< FALSE, if data can be safely shared between the source and target
                                              *   SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   )
{
   SCIP_VAR** fixedvars = NULL;
   SCIP_Real* fixedvals = NULL;
   int nfixedvars = 0;
   SCIP_Bool global = TRUE;
   SCIP_Bool original = TRUE;
   SCIP_Bool useconscompression = FALSE;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(suffix != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyOrig", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyOrig", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   SCIP_CALL( doCopy(sourcescip, targetscip, varmap, consmap, suffix, fixedvars, fixedvals, nfixedvars,
         useconscompression, global, original, enablepricing, threadsafe, passmessagehdlr, valid) );

   return SCIP_OKAY;
}

/** copies source SCIP original problem to target SCIP but compresses constraints
 *
 *  constraint compression is performed by removing fixed variables immediately
 *  during constraint creation if the involved constraint handlers support
 *  compression
 *
 *  the copying process is done in the following order:
 *  1) copy the plugins
 *  2) copy the settings
 *  3) create problem data in target-SCIP and copy the problem data of the source-SCIP
 *  4) copy all original variables
 *     a) fix all variable copies specified by \p fixedvars, \p fixedvals, and \p nfixedvars
 *     b) enable constraint compression
 *  5) copy all constraints
 *
 *  @note all variables and constraints which are created in the target-SCIP are not (user) captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *        Also, 'passmessagehdlr' should be set to FALSE.
 *  @note the 'threadsafe' parameter should only be set to TRUE if you are absolutely certain that the source and target
 *        SCIP instances will be solved in parallel. The usual case is to set this to FALSE, since thread safety
 *        typically incurs a performance cost.
 *  @note Do not change the source SCIP environment during the copying process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  @pre This method can be called if targetscip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @note sourcescip stage does not get changed
 *
 *  @note targetscip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyOrigConsCompression(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           suffix,             /**< optional suffix for problem name inside the target SCIP */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? If TRUE, pricer
                                              *   plugins will be copied and activated, and the modifiable flag of
                                              *   constraints will be respected. If FALSE, valid will be set to FALSE, when
                                              *   there are pricers present */
   SCIP_Bool             threadsafe,         /**< FALSE, if data can be safely shared between the source and target
                                              *   SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   )
{
   SCIP_Bool original = TRUE;
   SCIP_Bool global = TRUE;
   SCIP_Bool useconscompression = TRUE;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(suffix != NULL);

   /* check stages for both, the source and the target SCIP data structure */
    SCIP_CALL( SCIPcheckStage(sourcescip, "SCIPcopyOrigConsCompression", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
    SCIP_CALL( SCIPcheckStage(targetscip, "SCIPcopyOrigConsCompression", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   /* copy the source problem data */
   SCIP_CALL( doCopy(sourcescip, targetscip, varmap, consmap, suffix, fixedvars, fixedvals, nfixedvars,
         useconscompression, global, original, enablepricing, threadsafe, passmessagehdlr, valid) );

   SCIP_CALL( SCIPsyncstoreRelease(&targetscip->syncstore) );
   targetscip->syncstore = sourcescip->syncstore;
   SCIP_CALL( SCIPsyncstoreCapture(targetscip->syncstore) );

   return SCIP_OKAY;
}

/** return updated time limit for a sub-SCIP */
static
SCIP_RETCODE getCopyTimelimit(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_Real*            timelimit           /**< pointer to store sub-SCIP time limit */
   )
{
   SCIP_CALL( SCIPgetRealParam(sourcescip, "limits/time", timelimit) );
   if( !SCIPisInfinity(sourcescip, *timelimit) )
      (*timelimit) -= SCIPgetSolvingTime(sourcescip);

   return SCIP_OKAY;
}

/** set updated time limit for a sub-SCIP */
static
SCIP_RETCODE copySofttimelimit(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   if( SCIPgetParam(targetscip, "limits/softtime") == NULL )
      return SCIP_OKAY;
   else
   {
      SCIP_Real timelimit = -1.0;

      SCIP_CALL( SCIPgetRealParam(sourcescip, "limits/softtime", &timelimit) );
      if( !SCIPisNegative(sourcescip, timelimit) )
      {
         timelimit -= SCIPgetSolvingTime(sourcescip);
         timelimit = MAX(0.0, timelimit);
      }

      SCIP_CALL( SCIPsetRealParam(targetscip, "limits/softtime", timelimit) );
   }
   return SCIP_OKAY;
}

/** return updated memory limit for a sub-SCIP */
static
SCIP_RETCODE getCopyMemlimit(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_Real*            memorylimit         /**< pointer to store sub-SCIP memory limit */
   )
{
   SCIP_CALL( SCIPgetRealParam(sourcescip, "limits/memory", memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(sourcescip, *memorylimit) )
      (*memorylimit) -= (SCIPgetMemUsed(sourcescip) + SCIPgetMemExternEstim(sourcescip))/1048576.0;

   return SCIP_OKAY;
}

/** checks if there is enough time and memory left for copying the sourcescip into a sub-SCIP and solve the sub-SCIP
 *
 *  This is the case if the time and memory limit that would be passed to the sub-SCIP are larger than 0.0
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcheckCopyLimits(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_Bool*            success             /**< pointer to store whether there is time and memory left to copy the
                                              *   problem and run the sub-SCIP */
   )
{
   SCIP_Real timelimit;

   SCIP_CALL( getCopyTimelimit(sourcescip, &timelimit) );

   if( sourcescip->set->misc_avoidmemout )
   {
      SCIP_Real memorylimit;

      /* try to avoid running into memory limit */
      SCIP_CALL( getCopyMemlimit(sourcescip, &memorylimit) );
      *success = timelimit > 0.0 && memorylimit > 2.0 * SCIPgetMemExternEstim(sourcescip) / 1048576.0;
   }
   else
      *success = timelimit > 0.0;

   return SCIP_OKAY;
}

/** copies limits from source SCIP to target SCIP
 *
 *  @note time and memory limit are reduced by the amount already spent in the source SCIP before installing the limit
 *        in the target SCIP
 *  @note all other limits are disabled and need to be enabled afterwards, if needed
 *
 *  @pre This method can be called if sourcescip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcopyLimits(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   SCIP_Real timelimit;
   SCIP_Real memorylimit;

   SCIP_CALL( getCopyTimelimit(sourcescip, &timelimit) );
   SCIP_CALL( getCopyMemlimit(sourcescip, &memorylimit) );

   /* avoid negative limits */
   if( timelimit < 0.0 )
      timelimit = 0.0;
   if( memorylimit < 0.0 )
      memorylimit = 0.0;

   /* set time and memory limit to the adjusted values */
   SCIP_CALL( SCIPsetRealParam(targetscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(targetscip, "limits/memory", memorylimit) );

   /* copy and adjust soft time limit (or disable it) */
   SCIP_CALL( copySofttimelimit(sourcescip, targetscip) );

   /* disable all other limits */
   SCIP_CALL( SCIPsetRealParam(targetscip, "limits/absgap", 0.0) );
   SCIP_CALL( SCIPsetIntParam(targetscip, "limits/bestsol", -1) );
   SCIP_CALL( SCIPsetRealParam(targetscip, "limits/gap", 0.0) );
   SCIP_CALL( SCIPsetLongintParam(targetscip, "limits/nodes", -1LL) );
   SCIP_CALL( SCIPsetIntParam(targetscip, "limits/restarts", -1) );
   SCIP_CALL( SCIPsetIntParam(targetscip, "limits/solutions", -1) );
   SCIP_CALL( SCIPsetLongintParam(targetscip, "limits/stallnodes", -1LL) );
   SCIP_CALL( SCIPsetLongintParam(targetscip, "limits/totalnodes", -1LL) );

   return SCIP_OKAY;
}

/** sets the working limits as well as common search parameters for the auxiliary problem
 *
 *  @note memory and time limits are not affected, and must be set using SCIPcopyLimits() instead
 */
SCIP_RETCODE SCIPsetCommonSubscipParams(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 subscip,            /**< target SCIP data structure, often a copy of \p sourcescip */
   SCIP_Longint          nsubnodes,          /**< nodelimit for subscip, or -1 for no limit */
   SCIP_Longint          nstallnodes,        /**< stall node limit for subscip, or -1 for no limit */
   int                   bestsollimit        /**< the limit on the number of best solutions found, or -1 for no limit */
   )
{
   SCIP_Bool useuct;

   assert(sourcescip != NULL);
   assert(subscip != NULL);

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
   SCIP_CALL( SCIPcopyLimits(sourcescip, subscip) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nsubnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", bestsollimit) );

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
   SCIP_CALL( SCIPgetBoolParam(sourcescip, "heuristics/useuctsubscip", &useuct) );
   if( useuct && SCIPfindNodesel(subscip, "uct") != NULL && !SCIPisParamFixed(subscip, "nodeselection/uct/stdpriority") )
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

   return SCIP_OKAY;
}
