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

/**@file   heur_multistart.h
 * @ingroup PRIMALHEURISTICS
 * @brief  multistart heuristic for convex and nonconvex MINLPs
 * @author Benjamin Mueller
 *
 * The heuristic applies multiple NLP local searches to a mixed-integer nonlinear program with, probably nonconvex,
 * constraints of the form \f$g_j(x) \le 0\f$. The algorithm tries to identify clusters which approximate the boundary
 * of the feasible set of the continuous relaxation by sampling and improving randomly generated points. For each
 * cluster we use a local search heuristic to find feasible solutions. The algorithm consists of the following four
 * steps:
 *
 * 1. sample points
 *
 *    Sample random points \f$ x^1, \ldots, x^n \f$ in the box \f$ [\ell,u] \f$. For an unbounded variable \f$ x_i \f$
 *    we consider \f$ [\ell_i,\ell_i + \alpha], [u_i - \alpha,u_i], \f$ or \f$ [-\alpha / 2, \alpha / 2]\f$ for an \f$
 *    \alpha > 0 \f$ depending on which bound is infinite.
 *
 * 2. reduce infeasibility
 *
 *   For each point \f$ x^i \f$ we use a gradient descent method to reduce the maximum infeasibility. We first compute
 *
 *    \f[
 *        d_j = -\frac{g_j(x^i)}{||\nabla g_j(x^i)||^2} \nabla g_j(x^i)
 *    \f]
 *
 *    and update the current point \f$ x^i \f$ with
 *
 *    \f[
 *        x^i := x^i + \frac{1}{n_j} \sum_{j} d_j
 *    \f]
 *
 *    where \f$ n_j \f$ is the number of strictly positive \f$ d_j \f$. The algorithm is called Constraint Consensus
 *    Method and has been introduced by <a
 *    href="http://www.sce.carleton.ca/faculty/chinneck/docs/ConstraintConsensusJoC.pdf">here </a>.
 *
 * 3. cluster points
 *
 *    We use a greedy algorithm to all of the resulting points of step 3. to find clusters which (hopefully) approximate
 *    the boundary of the feasible set locally. Points with a too large violations will be filtered.
 *
 * 4. solve sub-problems
 *
 *    Depending on the current setting, we solve a sub-problem for each identified cluster. The default strategy is to
 *    compute a starting point for the sub-NLP heuristic (see @ref heur_subnlp.h) by using a linear combination of the
 *    points in a cluster \f$ C \f$, i.e.,
 *
 *    \f[
 *        s := \sum_{x \in C} \lambda_x x
 *    \f]
 *
 *    Since the sub-NLP heuristic requires a starting point which is integer feasible we round each fractional
 *    value \f$ s_i \f$ to its closest integer.
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_MULTISTART_H__
#define __SCIP_HEUR_MULTISTART_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the multistart primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurMultistart(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
