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

/**@file   nlhdlr_soc.h
 * @ingroup NLHDLRS
 * @brief  soc nonlinear handler
 *
 * @author Benjamin Mueller
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_SOC_H__
#define __SCIP_NLHDLR_SOC_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes SOC nonlinear handler in nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrSoc(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLHDLRS
 *
 * @{
 *
 * @name SOC nonlinear handler
 *
 * This nonlinear handler detects second-order cone constraints in the extended formulation and provides specialized separation functionality.
 *
 * @{
 */

/** checks whether constraint is SOC representable in original variables and returns the SOC representation
 *
 * The SOC representation has the form:
 * \f$\sqrt{\sum_{i=1}^{n} (v_i^T x + \beta_i)^2} - v_{n+1}^T x - \beta_{n+1} \lessgtr 0\f$,
 * where \f$n+1 = \text{nterms}\f$ and the inequality type is given by sidetype (`SCIP_SIDETYPE_RIGHT` if inequality
 * is \f$\leq\f$, `SCIP_SIDETYPE_LEFT` if \f$\geq\f$).
 *
 * For each term (i.e. for each \f$i\f$ in the above notation as well as \f$n+1\f$), the constant \f$\beta_i\f$ is given by the
 * corresponding element `offsets[i-1]` and `termbegins[i-1]` is the starting position of the term in arrays
 * `transcoefs` and `transcoefsidx`. The overall number of nonzeros is `termbegins[nterms]`.
 *
 * Arrays `transcoefs` and `transcoefsidx` have size `termbegins[nterms]` and define the linear expressions \f$v_i^T x\f$
 * for each term. For a term \f$i\f$ in the above notation, the nonzeroes are given by elements
 * `termbegins[i-1]...termbegins[i]` of `transcoefs` and `transcoefsidx`. There may be no nonzeroes for some term (i.e.,
 * constant terms are possible). `transcoefs` contains the coefficients \f$v_i\f$ and `transcoefsidx` contains positions of
 * variables in the `vars` array.
 *
 * The `vars` array has size `nvars` and contains \f$x\f$ variables; each variable is included at most once.
 *
 * The arrays should be freed by calling SCIPfreeSOCArraysNonlinear().
 *
 * This function uses the methods that are used in the detection algorithm of the SOC nonlinear handler.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPisSOCNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_Bool             compeigenvalues,    /**< whether eigenvalues should be computed to detect complex cases */
   SCIP_Bool*            success,            /**< pointer to store whether SOC structure has been detected */
   SCIP_SIDETYPE*        sidetype,           /**< pointer to store which side of cons is SOC representable; only
                                               valid when success is TRUE */
   SCIP_VAR***           vars,               /**< variables (x) that appear on both sides; no duplicates are allowed */
   SCIP_Real**           offsets,            /**< offsets of both sides (beta_i) */
   SCIP_Real**           transcoefs,         /**< non-zeros of linear transformation vectors (v_i) */
   int**                 transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int**                 termbegins,         /**< starting indices of transcoefs for each term */
   int*                  nvars,              /**< total number of variables appearing (i.e. size of vars) */
   int*                  nterms              /**< number of summands in the SQRT +1 for RHS (n+1) */
   );

/** frees arrays created by SCIPisSOCNonlinear() */
SCIP_EXPORT
void SCIPfreeSOCArraysNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variables that appear on both sides (x) */
   SCIP_Real**           offsets,            /**< offsets of both sides (beta_i) */
   SCIP_Real**           transcoefs,         /**< non-zeros of linear transformation vectors (v_i) */
   int**                 transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int**                 termbegins,         /**< starting indices of transcoefs for each term */
   int                   nvars,              /**< total number of variables appearing */
   int                   nterms              /**< number of summands in the SQRT +1 for RHS (n+1) */
   );

/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_SOC_H__ */
