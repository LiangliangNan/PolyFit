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

/**@file   sepa_eccuts.h
 * @ingroup SEPARATORS
 * @brief  edge concave cut separator
 * @author Benjamin Mueller
 *
 * We call \f$ f \f$ an edge-concave function on a polyhedron \f$P\f$ iff it is concave in all edge directions of
 * \f$P\f$. For the special case \f$ P = [\ell,u]\f$ this is equivalent to \f$f\f$ being concave componentwise.
 *
 * Since the convex envelope of an edge-concave function is a polytope, the value of the convex envelope for a
 * \f$ x \in [\ell,u] \f$ can be obtained by solving the following LP:
 *
 * \f{align}{
 *     \min \, & \sum_i \lambda_i f(v_i)  \\
 *     s.t. \, & \sum_i \lambda_i v_i = x \\
 *             & \sum_i \lambda_i = 1
 * \f}
 * where \f$ \{ v_i \} \f$ are the vertices of the domain \f$ [\ell,u] \f$. Let \f$ (\alpha, \alpha_0) \f$ be the dual
 * solution of this LP. It can be shown that \f$ \alpha' x + \alpha_0 \f$ is a facet of the convex envelope of \f$ f \f$
 * if \f$ x \f$ is in the interior of \f$ [\ell,u] \f$.
 *
 * We use this as follows:  We transform the problem to the unit box \f$ [0,1]^n \f$ by using a linear affine
 * transformation \f$ T(x) = Ax + b \f$ and perturb \f$ T(x) \f$ if it is not an interior point.
 * This has the advantage that we do not have to update the matrix of the LP for different edge-concave functions.
 *
 * For a given quadratic constraint \f$ g(x) := x'Qx + b'x + c \le 0 \f$ we decompose \f$ g \f$ into several
 * edge-concave aggregations and a remaining part, e.g.,
 *
 * \f[
 *             g(x) = \sum_{i = 1}^k f_i(x) + r(x)
 * \f]
 *
 * where each \f$ f_i \f$ is edge-concave. To separate a given solution \f$ x \f$, we compute a facet of the convex
 * envelope \f$ \tilde f \f$ for each edge-concave function \f$ f_i \f$ and an underestimator \f$ \tilde r \f$
 * for \f$ r \f$. The resulting cut looks like:
 *
 * \f[
 *             \tilde f_i(x) + \tilde r(x) \le 0
 * \f]
 *
 * We solve auxiliary MIP problems to identify good edge-concave aggregations. From the literature it is known that the
 * convex envelope of a bilinear edge-concave function \f$ f_i \f$ differs from McCormick iff in the graph
 * representation of \f$ f_i \f$ there exist a cycle with an odd number of positive weighted edges. We look for a
 * subgraph of the graph representation of the quadratic function \f$ g(x) \f$ with the previous property using a model
 * based on binary flow arc variables.
 *
 * This separator is currently disabled by default. It requires additional
 * tuning to be enabled by default. However, it may be useful to enable
 * it on instances with nonconvex quadratic constraints, in particular boxQPs.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_ECCUTS_H__
#define __SCIP_SEPA_ECCUTS_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the edge-concave separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaEccuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
