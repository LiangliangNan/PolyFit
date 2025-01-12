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

/**@file   nlhdlr_bilinear.h
 * @ingroup NLHDLRS
 * @brief  bilinear nonlinear handler
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_BILINEAR_H__
#define __SCIP_NLHDLR_BILINEAR_H__

#include "scip/scip.h"
#include "scip/pub_nlhdlr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes bilinear nonlinear handler in nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrBilinear(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLHDLRS
 * @{
 *
 * @name Bilinear nonlinear handler
 *
 * This nonlinear handler detects and collects bilinear terms and provides specialized propagation and estimation functionality.
 *
 * @{
 */

/** returns an array of expressions that have been detected by the bilinear nonlinear handler */
SCIP_EXPORT
SCIP_EXPR** SCIPgetExprsBilinear(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** returns an array of nonlinear handler expressions data of expressions that have been detected by the bilinear nonlinear handler */
SCIP_EXPORT
SCIP_NLHDLREXPRDATA** SCIPgetExprsdataBilinear(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** returns the total number of expressions that have been detected by the bilinear nonlinear handler */
SCIP_EXPORT
int SCIPgetNExprsBilinear(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** adds a globally valid inequality of the form \f$\text{xcoef}\cdot x \leq \text{ycoef} \cdot y + \text{constant}\f$ to a product expression of the form \f$x\cdot y\f$ */
SCIP_EXPORT
SCIP_RETCODE SCIPaddIneqBilinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_EXPR*            expr,               /**< product expression */
   SCIP_Real             xcoef,              /**< x coefficient */
   SCIP_Real             ycoef,              /**< y coefficient */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool*            success             /**< buffer to store whether inequality has been accepted */
   );

/** computes coefficients of linearization of a bilinear term in a reference point */
SCIP_EXPORT
void SCIPaddBilinLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             refpointx,          /**< point where to linearize first  variable */
   SCIP_Real             refpointy,          /**< point where to linearize second variable */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   );

/** computes coefficients of McCormick under- or overestimation of a bilinear term */
SCIP_EXPORT
void SCIPaddBilinMcCormick(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   );

/** computes coefficients of linearization of a bilinear term in a reference point when given a linear inequality
 *  involving only the variables of the bilinear term
 *
 *  @note the formulas are extracted from "Convex envelopes of bivariate functions through the solution of KKT systems"
 *        by Marco Locatelli
 */
SCIP_EXPORT
void SCIPcomputeBilinEnvelope1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real             xcoef,              /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             ycoef,              /**< y coefficient of linear inequality */
   SCIP_Real             constant,           /**< constant of linear inequality */
   SCIP_Real* RESTRICT   lincoefx,           /**< buffer to store coefficient of first  variable in linearization */
   SCIP_Real* RESTRICT   lincoefy,           /**< buffer to store coefficient of second variable in linearization */
   SCIP_Real* RESTRICT   linconstant,        /**< buffer to store constant of linearization */
   SCIP_Bool* RESTRICT   success             /**< buffer to store whether linearization was successful */
   );

/** computes coefficients of linearization of a bilinear term in a reference point when given two linear inequalities
 *  involving only the variables of the bilinear term
 *
 *  @note the formulas are extracted from "Convex envelopes of bivariate functions through the solution of KKT systems"
 *        by Marco Locatelli
 */
SCIP_EXPORT
void SCIPcomputeBilinEnvelope2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real             alpha1,             /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             beta1,              /**< y coefficient of linear inequality */
   SCIP_Real             gamma1,             /**< constant of linear inequality */
   SCIP_Real             alpha2,             /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             beta2,              /**< y coefficient of linear inequality */
   SCIP_Real             gamma2,             /**< constant of linear inequality */
   SCIP_Real* RESTRICT   lincoefx,           /**< buffer to store coefficient of first  variable in linearization */
   SCIP_Real* RESTRICT   lincoefy,           /**< buffer to store coefficient of second variable in linearization */
   SCIP_Real* RESTRICT   linconstant,        /**< buffer to store constant of linearization */
   SCIP_Bool* RESTRICT   success             /**< buffer to store whether linearization was successful */
   );

/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_BILINEAR_H__ */
