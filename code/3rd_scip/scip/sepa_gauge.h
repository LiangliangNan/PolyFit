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

/**@file   sepa_gauge.h
 * @ingroup SEPARATORS
 * @brief  gauge separator
 * @author Felipe Serrano
 *
 * This separator receives a point \f$ x_0 \f$ to separate and, given an interior point \f$ \bar x \f$, finds the
 * intersection between the boundary of a convex relaxation of the current problem and the segment joining \f$ x_0 \f$
 * and \f$ \bar x \f$. Then it generates gradient cuts at the intersection.
 *
 * The interior point \f$ \bar x \f$ is computed only once, by solving
 * \f{align}{
 *      \min \; & t \\
 *      s.t. \; & g_j(x) \le t & \forall j=1,\ldots,m \\
 *      & l_k(x) \le 0 & \forall k=1,\ldots,p
 * \f}
 * where each \f$ g_j \f$ is a convex function and \f$ l_k \f$ is a linear function and
 * \f[
 *      C = \{ x \colon g_j(x) \le 0 \, \forall j=1,\ldots,m, l_k(x) \le 0 \, \forall k=1,\ldots,p \}
 * \f]
 * is a convex relaxation of the current problem.
 * If we can not find an interior solution, the separator will not be executed again.
 *
 * Note that we do not try to push the linear constraints into the interior, i.e. we use \f$ l_k(x) \le 0 \f$ instead
 * of \f$ l_k(x) \le t \f$, since some of the inequalities might actually be equalities, forcing \f$ t \f$ to zero.
 * We also use an arbitrary lower bound on \f$ t \f$ to handle the case when \f$ C \f$ is unbounded.
 *
 * By default, the separator, if enabled, runs only if the convex relaxation has at least two nonlinear convex constraints.
 *
 * In order to compute the boundary point, we consider only nonlinear convex constraints that are violated by the point
 * we want to separate. These constraints define a convex region for which \f$ \bar x \f$ is an interior point. Then,
 * a binary search is perform on the segment \f$[\bar x, x_0]\f$ in order to find the boundary point. Gradient cuts are
 * computed for each of these nonlinear convex constraints which are active at the boundary point.
 *
 * Technical details:
 *    - We consider a constraint for the binary search only when its violation is larger than \f$ 10^{-4} \f$, see
 *    MIN_VIOLATION in sepa_gauge.c. The reason is that if the violation is too small, chances are that the point in the
 *    boundary is in the interior for this constraint and we wouldn't generate a cut for it anyway. On the other hand,
 *    even if we generate a cut for this constraint, it is likely that the boundary point is very close to the point to
 *    separate. Hence the cut generated would be very similar to the gradient cut at the point to separate.
 *    - Before separating, if a slight perturbation of the interior point in the direction of the point to separate
 *    gives a point outside the region, we do not separate. The reason is that the interior point we computed could be
 *    almost at the boundary and the segment \f$[\bar x, x_0]\f$ could be tangent to the region. In that case, the cuts
 *    we generate will not separate \f$ x_0 \f$ from the feasible region.
 *
 * This separator is currently disabled by default. It requires additional
 * tuning to be enabled by default. However, it may be useful to enable
 * it on instances with convex nonlinear constraints if SCIP spends
 * many iterations in the separation loop without doing sufficient progress.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_GAUGE_H__
#define __SCIP_SEPA_GAUGE_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the gauge separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaGauge(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
