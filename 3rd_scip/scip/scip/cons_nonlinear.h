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

/**@file   cons_nonlinear.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for nonlinear constraints \f$\textrm{lhs} \leq \sum_{i=1}^n a_ix_i + \sum_{j=1}^m c_jf_j(x) \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_NONLINEAR_H__
#define __SCIP_CONS_NONLINEAR_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** upgrading method for nonlinear constraints into more specific constraints
 *
 * the method might upgrade a nonlinear constraint into a set of upgrade constraints
 * the caller provided an array upgdconss to store upgrade constraints
 * the length of upgdconss is given by upgdconsssize
 * if an upgrade is not possible, set *nupgdconss to zero
 * if more than upgdconsssize many constraints shall replace cons, the function
 * should return the required number as negated value in *nupgdconss
 * i.e., if cons should be replaced by 3 constraints, the function should set
 * *nupgdconss to -3 and return with SCIP_OKAY
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cons            : the nonlinear constraint to upgrade
 *  - nupgdconss      : pointer to store number of constraints that replace this constraint
 *  - upgdconss       : array to store constraints that replace this constraint
 *  - upgdconsssize   : length of the provided upgdconss array
 */
#define SCIP_DECL_NONLINCONSUPGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS* cons, \
      int* nupgdconss, SCIP_CONS** upgdconss, int upgdconsssize)

/** reformulation method for expression graph nodes
 *
 * The method might reformulate a node in an expression graph by adding
 * auxiliary constraints and/or variables.
 * The caller provided an expression graph node which is to be reformulated.
 * If the method takes action, it has to return the node that should replace
 * the given node in *reformnode. The caller will then ensure that all parents of
 * node will use *reformnode, so node may be freed.
 * If the method does not do any reformulation, it shall return NULL in *reformnode.
 * The counter naddcons can be used to setup the names of added variables/constraints.
 * The method should increase this counter by the number of added constraints.
 * The method has to ensure that the reformulated node, if still valid,
 * has valid bound and curvature information.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - exprgraph       : the expression graph which node to reformulate
 *  - node            : the expression graph node to reformulate
 *  - naddcons        : counter on number of added constraints so far
 *
 *  output:
 *  - naddcons        : to be increased by number of additionally added constraints
 *  - reformnode      : reformulated node to replace node with, or NULL if no reformulation
 */
#define SCIP_DECL_EXPRGRAPHNODEREFORM(x) SCIP_RETCODE x (SCIP* scip,    \
      SCIP_EXPRGRAPH* exprgraph, SCIP_EXPRGRAPHNODE* node,              \
      int* naddcons, SCIP_EXPRGRAPHNODE** reformnode)

/** creates the handler for nonlinear constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Nonlinear Constraints
 *
 * @{
 *
 * This constraint handler handles constraints of the form
 * \f[
 *   \textrm{lhs} \leq \sum_{i=1}^n a_ix_i + \sum_{j=1}^m c_jf_j(x) \leq \textrm{rhs},
 * \f]
 * where \f$a_i\f$ and \f$c_j\f$ are coefficients and
 * \f$f_j(x)\f$ are nonlinear functions (given as expression tree).
 *
 * Constraints are enforced by separation, domain propagation, and spatial branching.
 *
 * For convex or concave \f$f_j(x)\f$, cuts that separate on the convex hull of the function graph are implemented.
 * For \f$f_j(x)\f$ that are not known to be convex or concave, a simple variant of linear estimation based on interval gradients is implemented.
 *
 * Branching is performed for variables in nonconvex terms, if the relaxation solution cannot be separated.
 *
 * This header offers the upgrade functionality to upgrade a general nonlinear constraint into a more specific constraint
 * via SCIP_DECL_NONLINCONSUPGD().
 *
 * Furthermore, the definition of callbacks used to reformulate an expression graph is offered by
 * SCIP_DECL_EXPRGRAPHNODEREFORM().
 *
 * Further, the function representation is stored in an expression graph, which allows to propagate variable domains
 * and constraint sides and offers a simple convexity check.
 * During presolve, the expression graph is reformulated, whereby new variables and constraints are created
 * such that for the remaining nonlinear constraints the functions \f$f_j(x)\f$ are known to be convex or concave.
 * See also
 *
 * @par
 * Stefan Vigerske@n
 * Decomposition of Multistage Stochastic Programs and a Constraint Integer Programming Approach to Mixed-Integer Nonlinear Programming@n
 * PhD Thesis, Humboldt-University Berlin, 2012, submitted.
 */

/** includes a nonlinear constraint upgrade method into the nonlinear constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeNonlinconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_NONLINCONSUPGD((*nonlinconsupgd)),/**< method to call for upgrading nonlinear constraint, or NULL */
   SCIP_DECL_EXPRGRAPHNODEREFORM((*nodereform)),/**< method to call for reformulating expression graph node, or NULL */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   );

/** creates and captures a nonlinear constraint
 * this variant takes expression trees as input
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   int                   nexprtrees,         /**< number of expression trees for nonlinear part of constraint */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees for nonlinear part of constraint */
   SCIP_Real*            nonlincoefs,        /**< coefficients for expression trees for nonlinear part, or NULL if all 1.0 */
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

/** creates and captures a nonlinear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsNonlinear(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  this variant takes expression trees as input
 *
 *  @see SCIPcreateConsNonlinear() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   int                   nexprtrees,         /**< number of expression trees for nonlinear part of constraint */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees for nonlinear part of constraint */
   SCIP_Real*            nonlincoefs,        /**< coefficients for expression trees for nonlinear part, or NULL if all 1.0 */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** creates and captures a nonlinear constraint
 * this variant takes a node of the expression graph as input and can only be used during presolving
 * it is assumed that the nonlinear constraint will be added to the transformed problem short after creation
 * the given exprgraphnode is captured in this method
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsNonlinear2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   SCIP_EXPRGRAPHNODE*   exprgraphnode,      /**< expression graph node associated to nonlinear expression */
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

/** creates and captures a nonlinear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsNonlinear2(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  this variant takes a node of the expression graph as input and can only be used during presolving
 *  it is assumed that the nonlinear constraint will be added to the transformed problem short after creation
 *  the given exprgraphnode is captured in this method
 *
 *  @see SCIPcreateConsNonlinear2() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicNonlinear2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   SCIP_EXPRGRAPHNODE*   exprgraphnode,      /**< expression graph node associated to nonlinear expression */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** adds a linear variable with coefficient to a nonlinear constraint */
EXTERN
SCIP_RETCODE SCIPaddLinearVarNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient of variable */
   );

/** sets the expression trees in a nonlinear constraint
 * constraint must not be active yet
 */
EXTERN
SCIP_RETCODE SCIPsetExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< new expression trees, or NULL if nexprtrees is 0 */
   SCIP_Real*            coefs               /**< coefficients of expression trees, or NULL if all 1.0 */
   );

/** adds expression trees to a nonlinear constraint
 * constraint must not be active yet
 */
EXTERN
SCIP_RETCODE SCIPaddExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< new expression trees, or NULL if nexprtrees is 0 */
   SCIP_Real*            coefs               /**< coefficients of expression trees, or NULL if all 1.0 */
   );

/** gets the nonlinear constraint as a nonlinear row representation */
EXTERN
SCIP_RETCODE SCIPgetNlRowNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   );

/** gets the number of variables in the linear term of a nonlinear constraint */
EXTERN
int SCIPgetNLinearVarsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the variables in the linear part of a nonlinear constraint */
EXTERN
SCIP_VAR** SCIPgetLinearVarsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the coefficients in the linear part of a nonlinear constraint */
EXTERN
SCIP_Real* SCIPgetLinearCoefsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the number of expression trees of a nonlinear constraint */
EXTERN
int SCIPgetNExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the expression trees of a nonlinear constraint */
EXTERN
SCIP_EXPRTREE** SCIPgetExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the coefficients of the expression trees of a nonlinear constraint */
EXTERN
SCIP_Real* SCIPgetExprtreeCoefsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the expression graph node of a nonlinear constraint */
EXTERN
SCIP_EXPRGRAPHNODE* SCIPgetExprgraphNodeNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the left hand side of a nonlinear constraint */
EXTERN
SCIP_Real SCIPgetLhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the right hand side of a nonlinear constraint */
EXTERN
SCIP_Real SCIPgetRhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** check the function of a nonlinear constraint for convexity/concavity, if not done yet */
EXTERN
SCIP_RETCODE SCIPcheckCurvatureNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the curvature of the nonlinear function of a nonlinear constraint
 *
 * The curvature is computed by summing up the curvature for each nonlinear summand.
 * To get the curvature for single summands, use SCIPgetExprtreeCurvaturesNonlinear().
 */
EXTERN
SCIP_RETCODE SCIPgetCurvatureNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             checkcurv,          /**< whether to check constraint curvature, if not checked before */
   SCIP_EXPRCURV*        curvature           /**< pointer to store curvature of constraint */
   );

/** gets the curvature of the expression trees (multiplied by their coefficient) of a nonlinear constraint */
EXTERN
SCIP_RETCODE SCIPgetExprtreeCurvaturesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             checkcurv,          /**< whether to check constraint curvature, if not checked before */
   SCIP_EXPRCURV**       curvatures          /**< buffer to store curvatures of exprtrees */
   );

/** computes the violation of a nonlinear constraint by a solution */
EXTERN
SCIP_RETCODE SCIPgetViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution which violation to calculate, or NULL for LP solution */
   SCIP_Real*            violation           /**< pointer to store violation of constraint */
   );

/** get index of a linear variable of a nonlinear constraint that may be decreased without making any other constraint infeasible, or -1 if none */
EXTERN
int SCIPgetLinvarMayDecreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** get index of a linear variable of a nonlinear constraint that may be increased without making any other constraint infeasible, or -1 if none */
EXTERN
int SCIPgetLinvarMayIncreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets expression graph of nonlinear constraint handler */
EXTERN
SCIP_EXPRGRAPH* SCIPgetExprgraphNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   );

/** given three points, constructs coefficient of equation for hyperplane generated by these three points
 * Three points a, b, and c are given.
 * Computes coefficients alpha, beta, gamma, and delta, such that a, b, and c, satisfy
 * alpha * x1 + beta * x2 + gamma * x3 = delta and gamma >= 0.0.
 */
EXTERN
SCIP_RETCODE SCIPcomputeHyperplaneThreePoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             a1,                 /**< first coordinate of a */
   SCIP_Real             a2,                 /**< second coordinate of a */
   SCIP_Real             a3,                 /**< third coordinate of a */
   SCIP_Real             b1,                 /**< first coordinate of b */
   SCIP_Real             b2,                 /**< second coordinate of b */
   SCIP_Real             b3,                 /**< third coordinate of b */
   SCIP_Real             c1,                 /**< first coordinate of c */
   SCIP_Real             c2,                 /**< second coordinate of c */
   SCIP_Real             c3,                 /**< third coordinate of c */
   SCIP_Real*            alpha,              /**< coefficient of first coordinate */
   SCIP_Real*            beta,               /**< coefficient of second coordinate */
   SCIP_Real*            gamma_,             /**< coefficient of third coordinate */
   SCIP_Real*            delta               /**< constant right-hand side */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
