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

/**@file   benderscut_int.h
 * @ingroup BENDERSCUTS
 * @brief  Generates a Laporte and Louveaux Benders' decomposition integer cut
 * @author Stephen J. Maher
 *
 * The classical Benders' decomposition algorithm is only applicable to problems with continuous second stage variables.
 * Laporte and Louveaux (1993) developed a method for generating cuts when Benders' decomposition is applied to problem
 * with discrete second stage variables. However, these cuts are only applicable when the master problem is a pure
 * binary problem.
 *
 * The integer optimality cuts are a point-wise underestimator of the optimal subproblem objective function value.
 * Similar to benderscuts_opt.c, an auxiliary variable, \f$\varphi\f$. is required in the master problem as a lower
 * bound on the optimal objective function value for the Benders' decomposition subproblem.
 *
 * Consider the Benders' decomposition subproblem that takes the master problem solution \f$\bar{x}\f$ as input:
 * \f[
 * z(\bar{x}) = \min\{d^{T}y : Ty \geq h - H\bar{x}, y \mbox{ integer}\}
 * \f]
 * If the subproblem is feasible, and \f$z(\bar{x}) > \varphi\f$ (indicating that the current underestimators are not
 * optimal) then the Benders' decomposition integer optimality cut can be generated from the optimal solution of the
 * subproblem. Let \f$S_{r}\f$ be the set of indicies for master problem variables that are 1 in \f$\bar{x}\f$ and
 * \f$L\f$ a known lowerbound on the subproblem objective function value.
 *
 * The resulting cut is:
 * \f[
 * \varphi \geq (z(\bar{x}) - L)(\sum_{i \in S_{r}}(x_{i} - 1) + \sum_{i \notin S_{r}}x_{i} + 1)
 * \f]
 *
 * Laporte, G. & Louveaux, F. V. The integer L-shaped method for stochastic integer programs with complete recourse
 * Operations Research Letters, 1993, 13, 133-142
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_INT_H__
#define __SCIP_BENDERSCUT_INT_H__


#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the integer optimality cut for Benders' decomposition cut and includes it in SCIP
 *
 *  @ingroup BenderscutIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBenderscutInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

#ifdef __cplusplus
}
#endif

#endif
