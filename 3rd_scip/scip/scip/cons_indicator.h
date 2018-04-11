/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for indicator constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
SCIP_RETCODE SCIPcreateConsBasicIndicatorLinCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< binary indicator variable (or NULL) */
   SCIP_CONS*            lincons,            /**< linear constraint */
   SCIP_VAR*             slackvar            /**< slack variable */
   );

/** adds variable to the inequality of the indicator constraint */
EXTERN
SCIP_RETCODE SCIPaddVarIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_VAR*             var,                /**< variable to add to the inequality */
   SCIP_Real             val                 /**< value of variable */
   );

/** gets the linear constraint corresponding to the indicator constraint (may be NULL) */
EXTERN
SCIP_CONS* SCIPgetLinearConsIndicator(
   SCIP_CONS*            cons                /**< indicator constraint */
   );

/** sets the linear constraint corresponding to the indicator constraint (may be NULL) */
EXTERN
SCIP_RETCODE SCIPsetLinearConsIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_CONS*            lincons             /**< linear constraint */
   );

/** sets binary indicator variable for indicator constraint */
EXTERN
SCIP_RETCODE SCIPsetBinaryVarIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_VAR*             binvar              /**< binary variable to add to the inequality */
   );

/** gets binary variable corresponding to indicator constraint */
EXTERN
SCIP_VAR* SCIPgetBinaryVarIndicator(
   SCIP_CONS*            cons                /**< indicator constraint */
   );

/** gets slack variable corresponding to indicator constraint */
EXTERN
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
EXTERN
SCIP_RETCODE SCIPsetSlackVarUb(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_CONS*            cons,                /**< indicator constraint */
   SCIP_Real             ub                   /**< upper bound for slack variable */
   );

/** checks whether indicator constraint is violated w.r.t. sol */
EXTERN
SCIP_Bool SCIPisViolatedIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   );

/** based on values of other variables, computes slack and binary variable to turn constraint feasible */
EXTERN
SCIP_RETCODE SCIPmakeIndicatorFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< indicator constraint */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            changed             /**< pointer to store whether the solution has been changed */
   );

/** based on values of other variables, computes slack and binary variable to turn all constraints feasible */
EXTERN
SCIP_RETCODE SCIPmakeIndicatorsFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< indicator constraint handler */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            changed             /**< pointer to store whether the solution has been changed */
   );

/** adds additional linear constraint that is not connected with an indicator constraint, but can be used for separation */
EXTERN
SCIP_RETCODE SCIPaddLinearConsIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< indicator constraint handler */
   SCIP_CONS*            lincons             /**< linear constraint */
   );

/** adds additional globally valid row that is not connected with an indicator constraint, but can be used for separation */
EXTERN
SCIP_RETCODE SCIPaddRowIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< indicator constraint handler */
   SCIP_ROW*             row                 /**< row to add */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
