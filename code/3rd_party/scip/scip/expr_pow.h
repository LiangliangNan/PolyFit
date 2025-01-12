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

/**@file   expr_pow.h
 * @ingroup EXPRHDLRS
 * @brief  power and signed power expression handlers
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_POW_H__
#define __SCIP_EXPR_POW_H__

#include "scip/scip.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup EXPRHDLRS
 *
 * @{
 *
 * @name Power and signed power expression
 *
 * These expression handler provide the power function, that is,
 * \f[
 *   x \mapsto \begin{cases}
 *     x^e & \textrm{if}\; x \geq 0\; \textrm{or}\; e\in\mathbb{Z}, \\
 *     \textrm{undefined}, & \textrm{otherwise}.
 *     \end{cases}
 * \f]
 * and the signed power function, that is,
 * \f[
 *   x \mapsto \textrm{sign}(e) |x|^e
 * \f]
 * for some exponent \f$e\f$.
 *
 * @{
 */

/** creates a power expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_Real             exponent,           /**< exponent of the power expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** creates a signpower expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprSignpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_Real             exponent,           /**< exponent of the power expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** indicates whether expression is of signpower-type */
SCIP_EXPORT
SCIP_Bool SCIPisExprSignpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** @}
  * @}
  */

/** creates the handler for power expression and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrPow(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates the handler for signed power expression and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrSignpower(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** computes coefficients of linearization of a square term in a reference point */
SCIP_EXPORT
void SCIPaddSquareLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             refpoint,           /**< point where to linearize */
   SCIP_Bool             isint,              /**< whether corresponding variable is a discrete variable, and thus linearization could be moved */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   );

/** computes coefficients of secant of a square term */
SCIP_EXPORT
void SCIPaddSquareSecant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             lb,                 /**< lower bound on variable */
   SCIP_Real             ub,                 /**< upper bound on variable */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of secant */
   SCIP_Real*            linconstant,        /**< buffer to add constant of secant */
   SCIP_Bool*            success             /**< buffer to set to FALSE if secant has failed due to large numbers or unboundedness */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_POW_H__ */
