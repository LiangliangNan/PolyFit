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

/**@file   scip_copy.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for problem copies
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_COPY_H__
#define __SCIP_SCIP_COPY_H__


#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup CopyMethods
 *
 * @{
 */

/** copies plugins from sourcescip to targetscip; in case that a constraint handler which does not need constraints
 *  cannot be copied, valid will return FALSE. All plugins can declare that, if their copy process failed, the
 *  copied SCIP instance might not represent the same problem semantics as the original.
 *  Note that in this case dual reductions might be invalid.
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *        Also, 'passmessagehdlr' should be set to FALSE.
 *
 *  @note Do not change the source SCIP environment during the copying process.
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
SCIP_EXPORT
SCIP_RETCODE SCIPcopyPlugins(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_Bool             copyreaders,        /**< should the file readers be copied */
   SCIP_Bool             copypricers,        /**< should the variable pricers be copied */
   SCIP_Bool             copyconshdlrs,      /**< should the constraint handlers be copied */
   SCIP_Bool             copyconflicthdlrs,  /**< should the conflict handlers be copied */
   SCIP_Bool             copypresolvers,     /**< should the presolvers be copied */
   SCIP_Bool             copyrelaxators,     /**< should the relaxation handler be copied */
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
   );

/** copies all Benders' decomposition plugins
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note the 'threadsafe' parameter must be set to TRUE if you are absolutely certain that the source and target
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
 *
 *  @post After calling this method targetscip reaches one of the following stages depending on if and when the solution
 *        process was interrupted:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note sourcescip stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcopyBenders(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables; if NULL the transfer of cuts is not possible */
   SCIP_Bool             threadsafe,         /**< FALSE, if data can be safely shared between the source and target
                                                  SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool*            valid               /**< pointer to store whether all plugins were validly copied */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPcopyProb(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   const char*           name                /**< problem name */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPcopyOrigProb(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL  */
   const char*           name                /**< problem name of target */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPenableConsCompression(
   SCIP*                 scip                /**< source SCIP data structure */
   );

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
SCIP_EXPORT
SCIP_Bool SCIPisConsCompressionEnabled(
   SCIP*                 scip                /**< source SCIP data structure */
   );

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
SCIP_EXPORT
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
   );

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
SCIP_EXPORT
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
   );

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
SCIP_EXPORT
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
   );

/** merges the histories of variables from a source SCIP into a target SCIP. The two data structures should point to
 *  different SCIP instances.
 *
 *  @note the notion of source and target is inverted here; \p sourcescip usually denotes a copied SCIP instance, whereas
 *        \p targetscip denotes the original instance
 */

SCIP_EXPORT
SCIP_RETCODE SCIPmergeVariableStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_VAR**            sourcevars,         /**< source variables for history merge, NULL entries are ignored */
   SCIP_VAR**            targetvars,         /**< target variables for history merge, NULL entries are ignored */
   int                   nvars               /**< number of variables in both variable arrays */
   );

/** merges the statistics of NLPIs from a source SCIP into a target SCIP
 *
 * The two SCIP instances should point to different SCIP instances.
 *
 *  @note the notion of source and target is inverted here; \p sourcescip usually denotes a copied SCIP instance, whereas
 *        \p targetscip denotes the original instance
 */
SCIP_EXPORT
void SCIPmergeNLPIStatistics(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_Bool             reset               /**< whether to reset statistics in sourcescip */
   );

/** translates a solution from a subscip to the main scip
 *
 * Variables that are relaxation-only in the master SCIP are set to 0 or the bound closest to 0. Such variables
 * are represented as NULL entry in the \p subvars array.
 *
 * @note This method allocates a new solution of the main \p scip that needs to be freed by the user.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPtranslateSubSol(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_VAR**            subvars,            /**< the variables from the subproblem in the same order as the main \p scip */
   SCIP_SOL**            newsol              /**< buffer to store pointer to created solution in main SCIP */
   );

/** checks the solutions from the subscip and adds the first one that is found feasible to the master SCIP
 *
 * Variables that are relaxation-only in the master SCIP are set to 0 or the bound closest to 0. Such variables
 * are represented as NULL entry in the \p subvars array.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPtranslateSubSols(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_VAR**            subvars,            /**< the variables from the subproblem in the same order as the main \p scip */
   SCIP_Bool*            success,            /**< pointer to store, whether new solution was found */
   int*                  solindex            /**< pointer to store solution index of stored solution, or NULL if not of interest */
   );

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
SCIP_EXPORT
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
   );

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
SCIP_EXPORT
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
   );

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
SCIP_EXPORT
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
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPconvertCutsToConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   int*                  ncutsadded          /**< pointer to store number of added cuts, or NULL */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPcopyCuts(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   int*                  ncutsadded          /**< pointer to store number of copied cuts, or NULL */
   );

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
SCIP_EXPORT
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
   );

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
SCIP_EXPORT
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
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPcopyParamSettings(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   );

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
SCIP_EXPORT
int SCIPgetSubscipDepth(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets depth of scip instance
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note SCIP stage does not get changed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
void SCIPsetSubscipDepth(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   newdepth            /**< new subscip depth */
   );

/** copies source SCIP to target SCIP; the copying process is done in the following order:
 *  1) copy the plugins
 *  2) copy the settings
 *  3) create problem data in target-SCIP and copy the problem data of the source-SCIP
 *  4) copy all active variables except those are marked as relaxation-only
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
SCIP_EXPORT
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
                                                  SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   );

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
 *  4) copy all active variables except those are marked as relaxation-only
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
SCIP_EXPORT
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
                                                  SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   );

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
SCIP_EXPORT
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
                                                  SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   );

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
SCIP_EXPORT
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
                                                  SCIP, otherwise TRUE. This is usually set to FALSE */
   SCIP_Bool             passmessagehdlr,    /**< should the message handler be passed */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   );

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
SCIP_EXPORT
SCIP_RETCODE SCIPcheckCopyLimits(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_Bool*            success             /**< pointer to store whether there is time and memory left to copy the
                                              *   problem and run the sub-SCIP */
   );

/** copies limits from source SCIP to target SCIP
 *
 *  @note time and memory limit are reduced by the amount already spent in the source SCIP before installing the limit
 *        in the target SCIP
 *  @note all other limits are disabled and need to be enabled afterwards, if needed
 *
 *  @see SCIPsetCommonSubscipParams() to set further working limits and other parameters commonly used for auxiliary problems
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
SCIP_EXPORT
SCIP_RETCODE SCIPcopyLimits(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   );


/** sets the working limits as well as common search parameters for the auxiliary problem
 *
 *  @note memory and time limits are not affected, and must be set using SCIPcopyLimits() instead
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsetCommonSubscipParams(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 subscip,            /**< target SCIP data structure, often a copy of \p sourcescip */
   SCIP_Longint          nsubnodes,          /**< nodelimit for subscip, or -1 for no limit */
   SCIP_Longint          nstallnodes,        /**< stall node limit for subscip, or -1 for no limit */
   int                   bestsollimit        /**< the limit on the number of best solutions found, or -1 for no limit */
   );


/**@} */

#ifdef __cplusplus
}
#endif

#endif
