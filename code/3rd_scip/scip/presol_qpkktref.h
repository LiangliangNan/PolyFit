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

/**@file   presol_qpkktref.h
 * @ingroup PRESOLVERS
 * @brief  qpkktref presolver
 * @author Tobias Fischer
 *
 * This presolver tries to add the KKT conditions as additional (redundant) constraints to the (mixed-binary) quadratic
 * program
 *  \f[
 *  \begin{array}{ll}
 *  \min         & x^T Q x + c^T x + d \\
 *               & A x \leq b, \\
 *               & x \in \{0, 1\}^{p} \times R^{n-p}.
 * \end{array}
 * \f]
 *
 * We first check if the structure of the program is like (QP), see the documentation of the function
 * checkConsQuadraticProblem().
 *
 * If the problem is known to be bounded (all variables have finite lower and upper bounds), then we add the KKT
 * conditions. For a continuous QPs the KKT conditions have the form
 * \f[
 *  \begin{array}{ll}
 *   Q x + c + A^T \mu = 0,\\
 *   Ax \leq b,\\
 *   \mu_i \cdot (Ax - b)_i = 0,    & i \in \{1, \dots, m\},\\
 *   \mu \geq 0.
 * \end{array}
 * \f]
 * where \f$\mu\f$ are the Lagrangian variables. Each of the complementarity constraints \f$\mu_i \cdot (Ax - b)_i = 0\f$
 * is enforced via an SOS1 constraint for \f$\mu_i\f$ and an additional slack variable \f$s_i = (Ax - b)_i\f$.
 *
 * For mixed-binary QPs, the KKT-like conditions are
 * \f[
 *  \begin{array}{ll}
 *   Q x + c + A^T \mu + I_J \lambda = 0,\\
 *   Ax \leq b,\\
 *   x_j \in \{0,1\}                    & j \in J,\\
 *   (1 - x_j) \cdot z_j = 0            & j \in J,\\
 *   x_j \cdot (z_j - \lambda_j) = 0    & j \in J,\\
 *   \mu_i \cdot (Ax - b)_i = 0         & i \in \{1, \dots, m\},\\
 *   \mu \geq 0,
 * \end{array}
 * \f]
 * where \f$J = \{1,\dots, p\}\f$, \f$\mu\f$ and \f$\lambda\f$ are the Lagrangian variables, and \f$I_J\f$ is the
 * submatrix of the \f$n\times n\f$ identity matrix with columns indexed by \f$J\f$. For the derivation of the KKT-like
 * conditions, see
 *
 *  Branch-And-Cut for Complementarity and Cardinality Constrained Linear Programs,@n
 *  Tobias Fischer, PhD Thesis (2016)
 *
 * Algorithmically:
 *
 * - we handle the quadratic term variables of the quadratic constraint like in the method
 *   presolveAddKKTQuadQuadraticTerms()
 * - we handle the bilinear term variables of the quadratic constraint like in the method presolveAddKKTQuadBilinearTerms()
 * - we handle the linear term variables of the quadratic constraint like in the method presolveAddKKTQuadLinearTerms()
 * - we handle linear constraints in the method presolveAddKKTLinearConss()
 * - we handle aggregated variables in the method presolveAddKKTAggregatedVars()
 *
 * we have a hashmap from each variable to the index of the dual constraint in the KKT conditions.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_QPKKTREF_H__
#define __SCIP_PRESOL_QPKKTREF_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the QP KKT reformulation presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresolQPKKTref(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
