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

/**@file   benderscut_feas.h
 * @ingroup BENDERSCUTS
 * @brief  Standard feasibility cuts for Benders' decomposition
 * @author Stephen J. Maher
 *
 * The classical Benders' decomposition feasibility cuts arise from an infeasible instance of the Benders' decomposition
 * subproblem.
 * Consider the linear Benders' decomposition subproblem that takes the master problem solution \f$\bar{x}\f$ as input:
 * \f[
 * z(\bar{x}) = \min\{d^{T}y : Ty \geq h - H\bar{x}, y \geq 0\}
 * \f]
 * If the subproblem is infeasible as a result of the solution \f$\bar{x}\f$, then the Benders' decomposition
 * feasibility cut can be generated from the dual ray. Let \f$w\f$ be the vector corresponding to the dual ray of the
 * Benders' decomposition subproblem. The resulting cut is:
 * \f[
 * 0 \geq w^{T}(h - Hx)
 * \f]
 *
 * Next, consider the nonlinear Benders' decomposition subproblem that takes the master problem solution \f$\bar{x}\f$ as input:
 * \f[
 * z(\bar{x}) = \min\{d^{T}y : g(\bar{x}, y) \leq 0, y \geq 0\}
 * \f]
 * If the subproblem is infeasible as a result of the solution \f$\bar{x}\f$, then the Benders' decomposition
 * feasibility cut can be generated from a minimal infeasible solution, i.e., a solution of the NLP
 * \f[
 *   \min\left\{\sum_i u_i : g(\bar{x}, y) \leq u, y \geq 0, u \geq 0\right\}
 * \f]
 * Let \f$\bar{y}\f$, \f$w\f$ be the vectors corresponding to the primal and dual solution of this auxiliary NLP.
 * The resulting cut is:
 * \f[
 * 0 \geq w^{T}\left(g(\bar{x},\bar{y}) + \nabla_x g(\bar{x},\bar{y}) (x - \bar{x})\right)
 * \f]
 * Note, that usually NLP solvers already provide a minimal infeasible solution when declaring the Benders'
 * decomposition subproblem as infeasible.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_FEAS_H__
#define __SCIP_BENDERSCUT_FEAS_H__


#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the Standard Feasibility Benders' decomposition cuts and includes it in SCIP
 *
 *  @ingroup BenderscutIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBenderscutFeas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

#ifdef __cplusplus
}
#endif

#endif
