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

/**@file   cons_indicator.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for indicator constraints
 * @author Marc Pfetsch
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_INDICATOR_H__
#define __SCIP_CONS_INDICATOR_H__


#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for indicator constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrIndicator(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Indicator Constraints
 *
 * @{
 *
 * An indicator constraint is given by a binary variable \f$z\f$ and an inequality \f$ax \leq
 * b\f$. It states that if \f$z = 1\f$ then \f$ax \leq b\f$ holds.
 *
 * This constraint is handled by adding a slack variable \f$s:\; ax - s \leq b\f$ with \f$s \geq
 * 0\f$. The constraint is enforced by fixing \f$s\f$ to 0 if \f$z = 1\f$.
 *
 * @note The constraint only implements an implication not an equivalence, i.e., it does not ensure
 * that \f$z = 1\f$ if \f$ax \leq b\f$ or equivalently if \f$s = 0\f$ holds.
 *
 * This constraint is equivalent to a linear constraint \f$ax - s \leq b\f$ and an SOS1 constraint on
 * \f$z\f$ and \f$s\f$ (at most one should be nonzero). In the indicator context we can, however,
 * separate more inequalities.
 */

/** creates and captures an indicator constraint
 *
 *  @note @a binvar is checked to be binary only later. This enables a change of the type in
 *  procedures reading an instance.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint (indicator or quadratic) */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable (or NULL) */
   int                   nvars,              /**< number of variables in the inequality */
   SCIP_VAR**            vars,               /**< array with variables of inequality (or NULL) */
   SCIP_Real*            vals,               /**< values of variables in inequality (or NULL) */
   SCIP_Real             rhs,                /**< rhs of the inequality */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? Usually set to TRUE. */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures an indicator constraint in its most basic version, i. e., all constraint flags are set to their
 *  basic value as explained for the method SCIPcreateConsIndicator(); all flags can be set via
 *  SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsIndicator() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint (indicator or quadratic) */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable (or NULL) */
   int                   nvars,              /**< number of variables in the inequality */
   SCIP_VAR**            vars,               /**< array with variables of inequality (or NULL) */
   SCIP_Real*            vals,               /**< values of variables in inequality (or NULL) */
   SCIP_Real             rhs                 /**< rhs of the inequality */
   );

/** creates and captures a indicator constraint in a more generic version.
 *
 *  The key difference from SCIPcreateConsIndicator() is the activeone and lessthanineq Booleans.
 *  If \f$z = o\f$, with \f$o\f$ the activeone flag, then:
 *  if lessthanineq then \f$a^T x \leq b\f$ holds, else the passed vectors are assumed to be of the form \f$a^T x \geq b\f$.
 *  The underlying linear constraint is always created as a less-than inequality.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsIndicatorGeneric(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint (indicator or quadratic) */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable (or NULL) */
   int                   nvars,              /**< number of variables in the inequality */
   SCIP_VAR**            vars,               /**< array with variables of inequality (or NULL) */
   SCIP_Real*            vals,               /**< values of variables in inequality (or NULL) */
   SCIP_Real             rhs,                /**< rhs of the inequality */
   SCIP_Bool             activeone,          /**< is the constraint active when the binary is 1? */
   SCIP_Bool             lessthanineq,       /**< is the linear constraint a less than RHS (TRUE) or greater than RHS (FALSE)? */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? Usually set to TRUE. */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures an indicator constraint with given linear constraint and slack variable
 *
 *  @note @a binvar is checked to be binary only later. This enables a change of the type in
 *  procedures reading an instance.
 *
 *  @note we assume that @a slackvar actually appears in @a lincons and we also assume that it takes
 *  the role of a slack variable!
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsIndicatorLinCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable (or NULL) */
   SCIP_CONS*            lincons,            /**< linear constraint */
   SCIP_VAR*             slackvar,           /**< slack variable */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? Usually set to TRUE. */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures an indicator constraint with given linear constraint and slack variable
 *  in a generic version, i. e., with a flag activeone indicating whether the constraint is active on
 *  value 1 or 0 of the binary variable.

 *  @note @a binvar is checked to be binary only later. This enables a change of the type in
 *  procedures reading an instance.
 *
 *  @note we assume that @a slackvar actually appears in @a lincons and we also assume that it takes
 *  the role of a slack variable!
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 *
 *  @see SCIPcreateConsIndicatorLinCons() for information about the basic constraint flag configuration
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsIndicatorGenericLinCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable (or NULL) */
   SCIP_CONS*            lincons,            /**< linear constraint */
   SCIP_VAR*             slackvar,           /**< slack variable */
   SCIP_Bool             activeone,          /**< is the constraint active when the binary is 1? */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? Usually set to TRUE. */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures an indicator constraint with given linear constraint and slack variable
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsIndicator(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h

 *  @note @a binvar is checked to be binary only later. This enables a change of the type in
 *  procedures reading an instance.
 *
 *  @note we assume that @a slackvar actually appears in @a lincons and we also assume that it takes
 *  the role of a slack variable!
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 *
 *  @see SCIPcreateConsIndicatorLinCons() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicIndicatorLinCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable (or NULL) */
   SCIP_CONS*            lincons,            /**< linear constraint */
   SCIP_VAR*             slackvar            /**< slack variable */
   );

/** adds variable to the inequality of the indicator constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPaddVarIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_VAR*             var,                /**< variable to add to the inequality */
   SCIP_Real             val                 /**< value of variable */
   );

/** gets the linear constraint corresponding to the indicator constraint (may be NULL) */
SCIP_EXPORT
SCIP_CONS* SCIPgetLinearConsIndicator(
   SCIP_CONS*            cons                /**< indicator constraint */
   );

/** sets the linear constraint corresponding to the indicator constraint (may be NULL) */
SCIP_EXPORT
SCIP_RETCODE SCIPsetLinearConsIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_CONS*            lincons             /**< linear constraint */
   );

/** sets binary indicator variable for indicator constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPsetBinaryVarIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_VAR*             binvar              /**< binary variable to add to the inequality */
   );

/** gets activation value of an indicator constraint, TRUE for active on 1, FALSE for active on 0 */
SCIP_EXPORT
SCIP_Bool SCIPgetActiveOnIndicator(
   SCIP_CONS*            cons                /**< indicator constraint */
   );

/** gets binary variable corresponding to indicator constraint. Returns the negative of the original binary variable if activeone was set to false */
SCIP_EXPORT
SCIP_VAR* SCIPgetBinaryVarIndicator(
   SCIP_CONS*            cons                /**< indicator constraint */
   );

/** similar to SCIPgetBinaryVarIndicator but returns the original binary variable passed by the user. */
SCIP_EXPORT
SCIP_VAR* SCIPgetBinaryVarIndicatorGeneric(
   SCIP_CONS*            cons                /**< indicator constraint */
   );

/** gets slack variable corresponding to indicator constraint */
SCIP_EXPORT
SCIP_VAR* SCIPgetSlackVarIndicator(
   SCIP_CONS*            cons                /**< indicator constraint */
   );

/** sets upper bound for slack variable corresponding to indicator constraint
 *
 *  Use with care if you know that the maximal violation of the corresponding constraint is at most @p ub. This bound
 *  might be improved automatically during the solution process.
 *
 *  @pre This method should only be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsetSlackVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_Real             ub                  /**< upper bound for slack variable */
   );

/** checks whether indicator constraint is violated w.r.t. sol */
SCIP_EXPORT
SCIP_Bool SCIPisViolatedIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   );

/** based on values of other variables, computes slack and binary variable to turn constraint feasible */
SCIP_EXPORT
SCIP_RETCODE SCIPmakeIndicatorFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            changed             /**< pointer to store whether the solution has been changed */
   );

/** based on values of other variables, computes slack and binary variable to turn all constraints feasible */
SCIP_EXPORT
SCIP_RETCODE SCIPmakeIndicatorsFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< indicator constraint handler */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            changed             /**< pointer to store whether the solution has been changed */
   );

/** adds additional linear constraint that is not connected with an indicator constraint, but can be used for separation */
SCIP_EXPORT
SCIP_RETCODE SCIPaddLinearConsIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< indicator constraint handler */
   SCIP_CONS*            lincons             /**< linear constraint */
   );

/** adds additional globally valid row that is not connected with an indicator constraint, but can be used for separation */
SCIP_EXPORT
SCIP_RETCODE SCIPaddRowIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< indicator constraint handler */
   SCIP_ROW*             row                 /**< row to add */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
