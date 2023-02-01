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

/**@file   sepa_convexproj.h
 * @ingroup SEPARATORS
 * @brief  convexproj separator
 * @author Felipe Serrano
 *
 * This separator receives a point \f$ x_0 \f$ to separate, projects it onto a convex relaxation
 * of the current problem and then generates gradient cuts at the projection.
 *
 * In more detail, the separators builds and stores a convex relaxation of the problem
 * \f[
 *      C = \{ x \colon g_j(x) \le 0 \, \forall j=1,\ldots,m \}
 * \f]
 * where each \f$ g_j \f$ is a convex function and computes the projection by solving
 * \f[
 *      \min || x - x_0 ||^2 \\
 * \f]
 * \f[
 *      s.t. \; g_j(x) \le 0 \, \forall j=1,\ldots,m \\
 * \f]
 *
 * By default, the separator runs only if the convex relaxation has at least one nonlinear convex function
 *
 * The separator generates cuts for constraints which were violated by the solution we want to separate and active
 * at the projection. If the projection problem is not solved to optimality, it still tries to add a cut at the
 * best solution found. In case that the projection problem is solved to optimality, it is guaranteed that a cut
 * separates the point. To see this, remember that \f$ z \f$ is the projection if and only if
 * \f[
 *      \langle x - z, z - x_0 \rangle \ge 0 \, \forall x \in C \\
 * \f]
 * This inequality is violated by for \f$ x = x_0 \f$. On the other hand, one of the optimality conditions of the
 * projection problem at the optimum looks like
 * \f[
 *      2 (z - x_0) + \sum_j \lambda_j \nabla g_j(z) = 0
 * \f]
 * Now suppose that the no gradient cut at \f$ z \f$ separates \f$ x_0 \f$, i.e.,
 * \f[
 *      g_j(z) + \langle \nabla g_j(z), x_0 - z \rangle \le 0
 * \f]
 * Multiplying each inequality with \f$ \lambda_j \ge 0 \f$ and summing up, we get the following contradiction:
 * \f[
 *      \langle -2(z - x_0), x_0 - z \rangle \le 0
 * \f]
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_CONVEXPROJ_H__
#define __SCIP_SEPA_CONVEXPROJ_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the convexproj separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeSepaConvexproj(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
