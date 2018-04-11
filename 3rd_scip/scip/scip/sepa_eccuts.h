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

/**@file   sepa_eccuts.h
 * @ingroup SEPARATORS
 * @brief  edge concave cut separator
 * @author Benjamin Müller
 *
 * We call \f$ f \f$ an edge-concave function on a polyhedron \f$P\f$ iff it is concave in all edge directions of
 * \f$P\f$. For the special case \f$ P = [\ell,u]\f$ this is equivalent to \f$f\f$ being concave componentwise.
 *
 * Since the convex envelope of an edge-concave function is a polytope, the value of the convex envelope for a
 * \f$ x \in [\ell,u] \f$ can be obtained by solving the following LP:
 *
 * \f[
 *              \min \, \sum_i \lambda_i f(v_i)
 * \f]
 * \f[
       s.t. \; \sum_i \lambda_i v_i = x
 * \f]
 * \f[
 *             \sum_i \lambda_i = 1
 * \f]
 *
 * where \f$ \{ v_i \} \f$ are the vertices of the domain \f$ [\ell,u] \f$. Let \f$ (\alpha, \alpha_0) \f$ be the dual
 * solution of this LP. It can be shown that \f$ \alpha' x + \alpha_0 \f$ is a facet of the convex envelope of \f$ f \f$
 * if \f$ x \f$ is in the interior of \f$ [\ell,u] \f$.
 *
 * We use this as follows:  We transform the problem to the unit box \f$ [0,1]^n \f$ by using an linear affine
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
 * convex envelope of an bilinear edge-concave function \f$ f_i \f$ differs from McCormick iff in the graph
 * representation of \f$ f_i \f$ there exist a cycle with an odd number of positive weighted edges. We look for a
 * subgraph of the graph representation of the quadratic function \f$ g(x) \f$ with the previous property using a model
 * based on binary flow arc variables.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_ECCUTS_H__
#define __SCIP_SEPA_ECCUTS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the edge-concave separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeSepaEccuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
