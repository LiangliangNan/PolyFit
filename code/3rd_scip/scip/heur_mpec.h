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

/**@file   heur_mpec.h
 * @ingroup PRIMALHEURISTICS
 * @brief  mpec primal heuristic
 * @author Felipe Serrano
 * @author Benjamin Mueller
 *
 * This heuristic is based on the paper:
 * @par
 * Lars Schewe and Martin Schmidt@n
 * [Computing Feasible Points for MINLPs with MPECs](http://www.optimization-online.org/DB_HTML/2016/12/5778.html)
 *
 * An MPEC is a mathematical program with complementarity constraint.
 * For example, the constraint \f$x \in \{0, 1\}\f$ as \f$x (1-x) = 0\f$
 * can be modeled as complementarity constraint \f$x = 0\f$ or \f$x = 1\f$.
 *
 * This heuristic applies only to mixed binary nonlinear problems.
 * The idea is to rewrite the MBNLP as MPEC and solve the MPEC instead (to a
 * a local optimum) by replacing each integrality constraint by the
 * complementarity constraint \f$x = 0\f$ or \f$x = 1\f$.
 * In principle, this MPEC can be reformulated to a NLP by rewriting this
 * constraint as equation \f$x (1-x) = 0\f$.
 * However, solving this NLP reformulation with a generic NLP solver will
 * often fail. One issue is that the reformulated complementarity constraints
 * will not, in general, satisfy constraint qualifications (for instance,
 * Slater's conditions, which requires the existence of a relative interior
 * point, will not be satisfied).
 * In order to increase the chances of solving the NLP reformulation of
 * the MPEC by a NLP solver, the heuristic applies a regularization
 * (proposed by Scholtes): it relaxes \f$x(1-x) = 0\f$ to
 * \f$x(1-x) \leq \theta\f$.
 *
 * So the heuristic proceeds as follows.
 * - Build the regularized NLP (rNLP) with a starting \f$\theta \in (0, \tfrac{1}{4}\f$.
 * - Give the current LP solution as starting point to the NLP solver.
 * - Solves rNLP and let \f$x^*\f$ be the best point found (if there is no point, abort).
 *   - If feasible, then reduce \f$\theta\f$ by a factor \f$\sigma\f$ and use
 *     its solution as the starting point for the next iteration.
 *
 *   - If the rNLP is found infeasible, but the regularization constraints are feasible, abort.
 *
 *   - If some variable violates the regularization constraint, i.e.,
 *   \f$x^*_i(1-x^*_i) > \tau\f$ then solve the rNLP again using its starting solution
 *   modified by \f$x_i = 0\f$ if \f$x^*_i > 0.5\f$ and \f$x_i = 1\f$ if \f$x^*_i < 0.5\f$.
 *   One possible explanation for this choice is that, assuming \f$x^*_i > 0.5\f$,
 *   if really \f$x_i = 1\f$ were a solution, then the NLP solver should not have had troubles
 *   pushing \f$x_i\f$ towards 1, making at least the regularization constraint feasible.
 *   Instead, it might be that there is a solution with \f$x_i = 0\f$, but since \f$x^*_i > 0.5\f$
 *   the NLP solver is having more problems pushing it to 0.
 *
 *   - If the modification of the starting point did not help finding a feasible solution,
 *   solve the problem again, but now fixing the problematic variables using the same criteria.
 *
 *   - If still we do not get a feasible solution, abort (note that the paper suggests to backtrack,
 *   but this might be just too expensive).
 *
 * - If the maximum binary infeasibility is small enough, call sub-NLP heuristic
 *   with binary variables fixed to the value suggested by \f$x^*\f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_MPEC_H__
#define __SCIP_HEUR_MPEC_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the mpec primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurMpec(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
