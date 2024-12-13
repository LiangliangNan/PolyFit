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

/**@file   prop_nlobbt.h
 * @ingroup PROPAGATORS
 * @brief  nonlinear OBBT propagator
 * @author Benjamin Mueller
 *
 * In Nonlinear Optimization-Based Bound Tightening (NLOBBT), we solve auxiliary NLPs of the form
 * \f[
 *      \min / \max \, x_i \\
 * \f]
 * \f[
 *      s.t. \; g_j(x) \le 0 \, \forall j=1,\ldots,m \\
 * \f]
 * \f[
 *      c'x \le \mathcal{U}
 * \f]
 * \f[
 *      x \in [\ell,u]
 * \f]
 *
 * where each \f$ g_j \f$ is a convex function and \f$ \mathcal{U} \f$ the solution value of the current
 * incumbent. Clearly, the optimal objective value of this nonlinear program provides a valid lower/upper bound on
 * variable \f$ x_i \f$.
 *
 * The propagator sorts all variables w.r.t. their occurrences in convex nonlinear constraints and solves sequentially
 * all convex NLPs. Variables which could be successfully tightened by the propagator will be prioritized in the next
 * call of a new node in the branch-and-bound tree. By default, the propagator requires at least one nonconvex
 * constraints to be executed. For purely convex problems, the benefit of having tighter bounds is negligible.
 *
 * By default, NLOBBT is only applied for non-binary variables. A reason for this can be found <a
 * href="http://dx.doi.org/10.1007/s10898-016-0450-4">here </a>. Variables which do not appear non-linearly in the
 * nonlinear constraints will not be considered even though they might lead to additional tightenings.
 *
 * After solving the NLP to optimize \f$ x_i \f$ we try to exploit the dual information to generate a globally valid
 * inequality, called Generalized Variable Bound (see @ref prop_genvbounds.h). Let \f$ \lambda_j \f$, \f$ \mu \f$, \f$
 * \alpha \f$, and \f$ \beta \f$ be the dual multipliers for the constraints of the NLP where \f$ \alpha \f$ and \f$
 * \beta \f$ correspond to the variable bound constraints. Because of the convexity of \f$ g_j \f$ we know that
 *
 * \f[
 *      g_j(x) \ge g_j(x^*) + \nabla g_j(x^*)(x-x^*)
 * \f]
 *
 * holds for every \f$ x^* \in [\ell,u] \f$. Let \f$ x^* \f$ be the optimal solution after solving the NLP for the case
 * of minimizing \f$ x_i \f$ (similiar for the case of maximizing \f$ x_i \f$). Since the NLP is convex we know that the
 * KKT conditions
 *
 * \f[
 *      e_i + \lambda' \nabla g(x^*) + \mu' c + \alpha - \beta = 0
 * \f]
 * \f[
 *      \lambda_j g_j(x^*) = 0
 * \f]
 *
 * hold. Aggregating the inequalities \f$ x_i \ge x_i \f$ and \f$ g_j(x) \le 0 \f$ leads to the inequality
 *
 * \f[
 *      x_i \ge x_i + \sum_{j} g_j(x_i)
 * \f]
 *
 * Instead of calling the (expensive) propagator during the tree search we can use this inequality to derive further
 * reductions on \f$ x_i \f$. Multiplying the first KKT condition by \f$ (x - x^*) \f$ and using the fact that each
 * \f$ g_j \f$ is convex we can rewrite the previous inequality to
 *
 * \f[
 *      x_i \ge (\beta - \alpha)'x + (e_i + \alpha - \beta) x^* + \mu \mathcal{U}.
 * \f]
 *
 * which is passed to the genvbounds propagator. Note that if \f$ \alpha_i \neq \beta_i \f$ we know that the bound of
 * \f$ x_i \f$ is the proof for optimality and thus no useful genvbound can be found.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_NLOBBT_H__
#define __SCIP_PROP_NLOBBT_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the nlobbt propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropNlobbt(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
