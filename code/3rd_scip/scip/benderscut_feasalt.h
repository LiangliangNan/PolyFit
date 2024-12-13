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

/**@file   benderscut_feasalt.h
 * @ingroup BENDERSCUTS
 * @brief  Alternative feasibility cuts for Benders' decomposition
 * @author Stephen J. Maher
 *
 * The alternative feasibility cut for Benders' decomposition uses the optimality cut generation code to generate a cut
 * that minimises the violation of the constraints.
 * Consider the linear Benders' decomposition subproblem that takes the master problem solution \f$\bar{x}\f$ as input:
 * \f[
 * z(\bar{x}) = \min\{d^{T}y : Ty \geq h - H\bar{x}, y \geq 0\}
 * \f]
 * If the subproblem is infeasible as a result of the solution \f$\bar{x}\f$, then some of the constraints are violated.
 * In this case, we define an alternative/auxiliary subproblem to find a solution that minimises the constraint
 * violations. Such a problem is given by
 * \f[
 * \min\{\mathbb{1}{T}v : Ty + v \geq h - H\bar{x}, y \geq 0, v \geq 0\}
 * \f]
 *
 * This auxiliary problem is guaranteed to always be feasible. Given a solution to this problem, it is possible to
 * generate a classical Benders' optimality cut. For such a cut, the reader is referred to \ref benderscut_opt.h.
 *
 * If the Benders' decomposition subproblem contains non-linear constraints, an equivalent auxiliary subproblem can be
 * formed to generate an alternative feasibility cut.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_FEASALT_H__
#define __SCIP_BENDERSCUT_FEASALT_H__


#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the Alternative Feasibility Benders' decomposition cuts and includes it in SCIP
 *
 *  @ingroup BenderscutIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBenderscutFeasalt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

#ifdef __cplusplus
}
#endif

#endif
