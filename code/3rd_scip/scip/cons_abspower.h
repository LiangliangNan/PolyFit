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

/**@file   cons_abspower.h
 * @ingroup CONSHDLRS
 * @brief  Constraint handler for absolute power constraints \f$\textrm{lhs} \leq \textrm{sign}(x+a) |x+a|^n + c z \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_ABSPOWER_H__
#define __SCIP_CONS_ABSPOWER_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for absolute power constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrAbspower(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Absolute Power Constraints
 *
 * @{
 *
 * This constraint handler handles constraints of the form
 * \f[
 *   \textrm{lhs} \leq \textrm{sign}(x+a) |x+a|^n + c z \leq \textrm{rhs}
 * \f]
 * for \f$n > 1.0\f$ a rational number, \f$a\f$ and \f$c\f$ arbitrary, and \f$x\f$ and \f$z\f$ variables.
 * Note that \f$x\f$ can have \f$-a\f$ in the interior of its domain.
 *
 * Constraints are enforced by separation, domain propagation, and spatial branching.
 *
 * Cuts that separate on the convex hull of the graph of \f$\textrm{sign}(x+a) |x+a|^n\f$ are generated as long as they separate the relaxation solution.
 * Otherwise, spatial branching on \f$x\f$ is applied.
 *
 * Further, domain propagation is implemented to propagate bound changes on \f$x\f$ onto \f$z\f$, and vice versa, and
 * repropagation is implemented to allow for conflict analysis.
 * During presolve, a pairwise comparison of absolute power constraints may allow to fix or aggregate some variables.
 * See also
 *
 * @par
 * Stefan Vigerske@n
 * Decomposition of Multistage Stochastic Programs and a Constraint Integer Programming Approach to Mixed-Integer Nonlinear Programming@n
 * PhD Thesis, Humboldt-University Berlin, 2012, submitted.
 *
 */

/** creates and captures a absolute power constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             x,                  /**< nonlinear variable x in constraint */
   SCIP_VAR*             z,                  /**< linear variable z in constraint */
   SCIP_Real             exponent,           /**< exponent n of |x+offset|^n term in constraint */
   SCIP_Real             xoffset,            /**< offset in |x+offset|^n term in constraint */
   SCIP_Real             zcoef,              /**< coefficient of z in constraint */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
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
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures an absolute power constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsAbspower(); all flags can be set via SCIPconsSetFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsAbspower() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             x,                  /**< nonlinear variable x in constraint */
   SCIP_VAR*             z,                  /**< linear variable z in constraint */
   SCIP_Real             exponent,           /**< exponent n of |x+offset|^n term in constraint */
   SCIP_Real             xoffset,            /**< offset in |x+offset|^n term in constraint */
   SCIP_Real             zcoef,              /**< coefficient of z in constraint */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** gets the absolute power constraint as a nonlinear row representation */
EXTERN
SCIP_RETCODE SCIPgetNlRowAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< a buffer where to store pointer to nonlinear row */
   );

/** gets nonlinear variable x in absolute power constraint */
EXTERN
SCIP_VAR* SCIPgetNonlinearVarAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   );

/** gets linear variable z in absolute power constraint */
EXTERN
SCIP_VAR* SCIPgetLinearVarAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   );

/** gets exponent in power term in absolute power constraint */
EXTERN
SCIP_Real SCIPgetExponentAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   );

/** gets offset in power term in absolute power constraint */
EXTERN
SCIP_Real SCIPgetOffsetAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   );

/** gets coefficient of linear variable in absolute power constraint */
EXTERN
SCIP_Real SCIPgetCoefLinearAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   );

/** gets left hand side in absolute power constraint */
EXTERN
SCIP_Real SCIPgetLhsAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   );

/** gets right hand side in absolute power constraint */
EXTERN
SCIP_Real SCIPgetRhsAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   );

/** gets the absolute violation of a absolute power constraint by a solution */
EXTERN
SCIP_Real SCIPgetViolationAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< absolute power constraint */
   SCIP_SOL*             sol                 /**< LP solution */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
