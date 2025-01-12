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

/**@file   cons_logicor.h
 * @ingroup CONSHDLRS
 * @brief  Constraint handler for logicor constraints \f$1^T x \ge 1\f$ (equivalent to set covering, but algorithms are suited for depth first search).
 * @author Tobias Achterberg
 * @author Michael Winkler
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_LOGICOR_H__
#define __SCIP_CONS_LOGICOR_H__


#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for logic or constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrLogicor(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Logicor Constraints
 *
 * @{
 *
 * This constraint handler handles a special type of linear constraints, namely
 * logic or constraints. These are equivalent to set covering constraints, but
 * are handled by special algorithms which are better suited for depth first search.
 * For a set of binary variables \f$x_i, i=1,\dots,n\f$, a logic or constraint has the form
 * \f[
 *   \sum_{i=1}^n x_i \ge 1.
 * \f]
 */

/** creates and captures a logic or constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures a logicor constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLogicor(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLogicor() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars                /**< array with variables of constraint entries */
   );

/** adds coefficient in logic or constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPaddCoefLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< logicor constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   );

/** gets number of variables in logic or constraint */
SCIP_EXPORT
int SCIPgetNVarsLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array of variables in logic or constraint */
SCIP_EXPORT
SCIP_VAR** SCIPgetVarsLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the dual solution of the logic or constraint in the current LP */
SCIP_EXPORT
SCIP_Real SCIPgetDualsolLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the dual Farkas value of the logic or constraint in the current infeasible LP */
SCIP_EXPORT
SCIP_Real SCIPgetDualfarkasLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the linear relaxation of the given logic or constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_EXPORT
SCIP_ROW* SCIPgetRowLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** cleans up (multi-)aggregations and fixings from logicor constraints */
SCIP_EXPORT
SCIP_RETCODE SCIPcleanupConssLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             onlychecked,        /**< should only checked constraints be cleaned up? */
   int*                  naddconss,          /**< pointer to count number of added (linear) constraints */
   int*                  ndelconss,          /**< pointer to count number of deleted (logicor) constraints */
   int*                  nchgcoefs           /**< pointer to count number of changed coefficients */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
