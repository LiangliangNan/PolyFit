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

/**@file   cons_linear.h
 * @ingroup CONSHDLRS
 * @brief  Constraint handler for linear constraints in their most general form, \f$lhs <= a^T x <= rhs\f$.
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Kati Wolter
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_LINEAR_H__
#define __SCIP_CONS_LINEAR_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif



/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrLinear(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Linear Constraints
 *
 * This constraint handler handles linear constraints in their most general form. That is,
 * \f[
 *   lhs \leq \sum_{i=1}^n a_i x_i \leq rhs
 * \f]
 * with \f$a_i \in Q, i = 1,\dots,n\f$, \f$lhs\in Q \cup \{-\infty\}\f$, \f$rhs\in Q \cup \{\infty\}\f$,
 * and decision variables \f$x_i, i = 1,\dots,n\f$ which can be binary, integer, or continuous.
 *
 * Furthermore, this header offers the upgrade functionality of a general linear constraint into a more specific
 * constraint, such as a knapsack constraint, via SCIP_DECL_LINCONSUPGD() and SCIPincludeLinconsUpgrade()
 *
 * @{
 */

typedef struct SCIP_LinConsUpgrade SCIP_LINCONSUPGRADE; /**< linear constraint update method */

/** upgrading method for linear constraints into more specific constraints
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cons            : the linear constraint to upgrade
 *  - nvars           : number of variables in the constraint
 *  - vars            : array with constraint variables
 *  - vals            : array with constraint coefficients
 *  - lhs             : left hand side of linear constraint
 *  - rhs             : right hand side of linear constraint
 *  - nposbin         : number of binary variables with positive coefficient
 *  - nnegbin         : number of binary variables with negative coefficient
 *  - nposint         : number of integer variables with positive coefficient
 *  - nnegint         : number of integer variables with negative coefficient
 *  - nposimpl        : number of implicit integer variables with positive coefficient (including implicit binary variables)
 *  - nnegimpl        : number of implicit integer variables with negative coefficient (including implicit binary variables)
 *  - nposimplbin     : number of implicit binary variables with positive coefficient
 *  - nnegimplbin     : number of implicit binary variables with negative coefficient
 *  - nposcont        : number of continuous variables with positive coefficient
 *  - nnegcont        : number of continuous variables with negative coefficient
 *  - ncoeffspone     : number of +1 coefficients
 *  - ncoeffsnone     : number of -1 coefficients
 *  - ncoeffspint     : number of positive integral coefficients other than +1
 *  - ncoeffsnint     : number of negative integral coefficients other than -1
 *  - ncoeffspfrac    : number of positive fractional coefficients
 *  - ncoeffsnfrac    : number of negative fractional coefficients
 *  - poscoeffsum     : sum of all positive coefficients
 *  - negcoeffsum     : sum of all negative coefficients
 *  - integral        : TRUE iff constraints activity value is always integral
 *  - upgdcons        : pointer to store the upgraded constraint
 */
#define SCIP_DECL_LINCONSUPGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS* cons, int nvars, SCIP_VAR** vars, SCIP_Real* vals, SCIP_Real lhs, SCIP_Real rhs, \
      int nposbin, int nnegbin, int nposint, int nnegint, int nposimpl, int nnegimpl, int nposimplbin, int nnegimplbin, int nposcont, int nnegcont, \
      int ncoeffspone, int ncoeffsnone, int ncoeffspint, int ncoeffsnint, int ncoeffspfrac, int ncoeffsnfrac, \
      SCIP_Real poscoeffsum, SCIP_Real negcoeffsum, SCIP_Bool integral, SCIP_CONS** upgdcons)

/** includes a linear constraint update method into the linear constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeLinconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int                   priority,           /**< priority of upgrading method */
   const char*           conshdlrname        /**< name of the constraint handler */
   );

/** creates and captures a linear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures a linear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLinear(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLinear() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** creates by copying and captures a linear constraint */
EXTERN
SCIP_RETCODE SCIPcopyConsLinear(
   SCIP*                 scip,               /**< target SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to store the created target constraint */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in source variable array */
   SCIP_VAR**            sourcevars,         /**< source variables of the linear constraints */
   SCIP_Real*            sourcecoefs,        /**< coefficient array of the linear constraint, or NULL if all coefficients are one */
   SCIP_Real             lhs,                /**< left hand side of the linear constraint */
   SCIP_Real             rhs,                /**< right hand side of the linear constraint */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
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
   SCIP_Bool*            valid               /**< pointer to store if the copying was valid */
   );

/** adds coefficient to linear constraint (if it is not zero) */
EXTERN
SCIP_RETCODE SCIPaddCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   );

/** changes coefficient of variable in linear constraint; deletes the variable if coefficient is zero; adds variable if
 *  not yet contained in the constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint and variable.
 *
 *  @note This method requires linear time to search for occurences of the variable in the constraint data.
 */
EXTERN
SCIP_RETCODE SCIPchgCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< new coefficient of constraint entry */
   );

/** deletes variable from linear constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint and variable.
 *
 *  @note This method requires linear time to search for occurences of the variable in the constraint data.
 */
EXTERN
SCIP_RETCODE SCIPdelCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   );

/** gets left hand side of linear constraint */
EXTERN
SCIP_Real SCIPgetLhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets right hand side of linear constraint */
EXTERN
SCIP_Real SCIPgetRhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** changes left hand side of linear constraint */
EXTERN
SCIP_RETCODE SCIPchgLhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of linear constraint */
EXTERN
SCIP_RETCODE SCIPchgRhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   );

/** gets the number of variables in the linear constraint */
EXTERN
int SCIPgetNVarsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the array of variables in the linear constraint; the user must not modify this array! */
EXTERN
SCIP_VAR** SCIPgetVarsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
EXTERN
SCIP_Real* SCIPgetValsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the activity of the linear constraint in the given solution
 *
 *  @note if the solution contains values at infinity, this method will return SCIP_INVALID in case the activity
 *        comprises positive and negative infinity contributions
 */
EXTERN
SCIP_Real SCIPgetActivityLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   );

/** gets the feasibility of the linear constraint in the given solution */
EXTERN
SCIP_Real SCIPgetFeasibilityLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   );

/** gets the dual solution of the linear constraint in the current LP */
EXTERN
SCIP_Real SCIPgetDualsolLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the dual Farkas value of the linear constraint in the current infeasible LP */
EXTERN
SCIP_Real SCIPgetDualfarkasLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
EXTERN
SCIP_ROW* SCIPgetRowLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** tries to automatically convert a linear constraint into a more specific and more specialized constraint */
EXTERN
SCIP_RETCODE SCIPupgradeConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_CONS**           upgdcons            /**< pointer to store upgraded constraint, or NULL if not successful */
   );

/** performs linear constraint type classification as used for MIPLIB
 *
 *  iterates through all linear constraints and stores relevant statistics in the linear constraint statistics \p linconsstats.
 *
 *  @note only constraints are iterated that belong to the linear constraint handler. If the problem has been presolved already,
 *  constraints that were upgraded to more special types such as, e.g., varbound constraints, will not be shown correctly anymore.
 *  Similarly, if specialized constraints were created through the API, these are currently not present.
 */
EXTERN
SCIP_RETCODE SCIPclassifyConstraintTypesLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LINCONSSTATS*    linconsstats        /**< linear constraint type classification */
   );


/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
