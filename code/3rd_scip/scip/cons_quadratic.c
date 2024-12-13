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

/**@file   cons_quadratic.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  some API functions of removed constraint handler for quadratic constraints \f$\textrm{lhs} \leq \sum_{i,j} a_{i,j} x_i x_j + \sum_i b_i x_i \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_nonlinear.h"
#include "scip/cons_quadratic.h"
#include "scip/expr_var.h"
#include "scip/expr_pow.h"
#include "scip/expr_product.h"

/** Creates and captures a quadratic nonlinear constraint.
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_j z_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 *
 *  @deprecated Use SCIPcreateConsQuadraticNonlinear() instead.
 */
SCIP_RETCODE SCIPcreateConsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< variables in linear part (x_i) or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< coefficients of variables in linear part (b_i) or NULL if nlinvars == 0 */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) or NULL if nquadterms == 0 */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) or NULL if nquadterms == 0 */
   SCIP_Real*            quadcoeffs,         /**< array with coefficients of quadratic terms (a_j) or NULL if nquadterms == 0 */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (l) */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation (u) */
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
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CALL( SCIPcreateConsQuadraticNonlinear(scip, cons, name, nlinvars, linvars, lincoefs, nquadterms, quadvars1, quadvars2, quadcoeffs, lhs, rhs,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   return SCIP_OKAY;
}

/** creates and captures a quadratic nonlinear constraint
 *  in its most basic variant, i.e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME()
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_jz_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @see SCIPcreateConsQuadratic() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 *
 *  @deprecated Use SCIPcreateConsBasicQuadraticNonlinear instead.
 */
SCIP_RETCODE SCIPcreateConsBasicQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (ell) */
   SCIP_Real             rhs                 /**< right hand side of quadratic equation (u) */
   )
{
   SCIP_CALL( SCIPcreateConsBasicQuadraticNonlinear(scip, cons, name, nlinvars, linvars, lincoefs, nquadterms, quadvars1, quadvars2, quadcoefs, lhs, rhs) );

   return SCIP_OKAY;
}

/** Adds a constant to the constraint function, that is, subtracts a constant from both sides
 *
 * @deprecated Use SCIPchgLhsNonlinear() and SCIPchgRhsNonlinear() instead.
 */
void SCIPaddConstantQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             constant            /**< constant to subtract from both sides */
   )
{
   SCIP_Real side;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   side = SCIPgetLhsNonlinear(cons);
   if( !SCIPisInfinity(scip, -side) )
   {
      SCIP_CALL_ABORT( SCIPchgLhsNonlinear(scip, cons, side-constant) );
   }

   side = SCIPgetRhsNonlinear(cons);
   if( !SCIPisInfinity(scip, side) )
   {
      SCIP_CALL_ABORT( SCIPchgRhsNonlinear(scip, cons, side-constant) );
   }
}

/** Adds a linear variable with coefficient to a quadratic constraint.
 *
 * @deprecated Use SCIPaddLinearVarNonlinear() instead.
 */
SCIP_RETCODE SCIPaddLinearVarQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient of variable */
   )
{
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   SCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, var, coef) );

   return SCIP_OKAY;
}

/** Adds a quadratic variable with linear and square coefficient to a quadratic constraint.
 *
 * @deprecated Use SCIPaddLinearVarNonlinear() and SCIPaddExprNonlinear() instead.
 */
SCIP_RETCODE SCIPaddQuadVarQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             lincoef,            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef             /**< square coefficient of variable */
   )
{
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   if( lincoef != 0.0 )
   {
      SCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, var, lincoef) );
   }

   if( sqrcoef != 0.0 )
   {
      SCIP_EXPR* varexpr;
      SCIP_EXPR* sqrexpr;

      SCIP_CALL( SCIPcreateExprVar(scip, &varexpr, var, NULL, NULL) );
      SCIP_CALL( SCIPcreateExprPow(scip, &sqrexpr, varexpr, 2.0, NULL, NULL) );

      SCIP_CALL( SCIPaddExprNonlinear(scip, cons, sqrexpr, sqrcoef) );

      SCIP_CALL( SCIPreleaseExpr(scip, &sqrexpr) );
      SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );
   }

   return SCIP_OKAY;
}

/** Adds a linear coefficient for a quadratic variable.
 *
 * Variable will be added with square coefficient 0.0 if not existing yet.
 *
 * @deprecated Use SCIPaddLinearVarNonlinear() instead.
 */
SCIP_RETCODE SCIPaddQuadVarLinearCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value to add to linear coefficient of variable */
   )
{
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   SCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, var, coef) );

   return SCIP_OKAY;
}

/** Adds a square coefficient for a quadratic variable.
 *
 * Variable will be added with linear coefficient 0.0 if not existing yet.
 *
 * @deprecated Use SCIPaddExprNonlinear() instead.
 */
SCIP_RETCODE SCIPaddSquareCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value to add to square coefficient of variable */
   )
{
   SCIP_EXPR* varexpr;
   SCIP_EXPR* sqrexpr;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   SCIP_CALL( SCIPcreateExprVar(scip, &varexpr, var, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &sqrexpr, varexpr, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPaddExprNonlinear(scip, cons, sqrexpr, coef) );

   SCIP_CALL( SCIPreleaseExpr(scip, &sqrexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );

   return SCIP_OKAY;
}

/** Adds a bilinear term to a quadratic constraint.
 *
 * Variables will be added with linear and square coefficient 0.0 if not existing yet.
 * If variables are equal, only the square coefficient of the variable is updated.
 *
 * @deprecated Use SCIPaddExprNonlinear() instead.
 */
SCIP_RETCODE SCIPaddBilinTermQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Real             coef                /**< coefficient of bilinear term */
   )
{
   SCIP_EXPR* varexprs[2];
   SCIP_EXPR* prodexpr;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[0], var1, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[1], var2, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &prodexpr, 2, varexprs, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddExprNonlinear(scip, cons, prodexpr, coef) );

   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[0]) );

   return SCIP_OKAY;
}

/** Gets the quadratic constraint as a nonlinear row representation.
 *
 * @deprecated Use SCIPgetNlRowNonlinear() instead.
 */
SCIP_RETCODE SCIPgetNlRowQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   )
{
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   SCIP_CALL( SCIPgetNlRowNonlinear(scip, cons, nlrow) );

   return SCIP_OKAY;
}

/** sets the left hand side of a quadratic constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint.
 *
 *  @deprecated Use SCIPchgLhsNonlinear() instead.
 */
SCIP_RETCODE SCIPchgLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   SCIP_CALL( SCIPchgLhsNonlinear(scip, cons, lhs) );

   return SCIP_OKAY;
}

/** sets the right hand side of a quadratic constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint.
 *
 *  @deprecated Use SCIPchgRhsNonlinear() instead.
 */
SCIP_RETCODE SCIPchgRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "nonlinear") == 0);

   SCIP_CALL( SCIPchgRhsNonlinear(scip, cons, rhs) );

   return SCIP_OKAY;
}
