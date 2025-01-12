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

/**@file   presol_tworowbnd.h
 * @ingroup DEFPLUGINS_PRESOL
 * @brief  do bound tightening by using two rows
 * @author Dieter Weninger
 * @author Patrick Gemander
 *
 * Perform bound tightening on two inequalities with some common variables.
 * Two possible methods are being used.
 *
 * 1. LP-bound
 * Let two constraints be given:
 * \f{eqnarray*}{
 *   A_{iR} x_R + A_{iS} x_S              \geq b_i\\
 *   A_{kR} x_R              + A_{kT} x_T \geq b_k
 * \f}
 * with \f$N\f$ the set of variable indexes, \f$R \subseteq N\f$, \f$S \subseteq N\f$, \f$T \subseteq N\f$,
 * \f$R \cap S = \emptyset\f$, \f$R \cap T = \emptyset\f$, \f$S \cap T = \emptyset\f$ and row indices \f$i \not= k\f$.
 *
 * Let \f$\ell\f$ and \f$u\f$ be bound vectors for x and solve the following two LPs
 * \f{eqnarray*}{
 *   L = \min \{ A_{kR} x_R : A_{iR} x_R + A_{iS} x_S \geq b_i, \ell \leq x \leq u \}\\
 *   U = \max \{ A_{kR} x_R : A_{iR} x_R + A_{iS} x_S \geq b_i, \ell \leq x \leq u \}
 * \f}
 * and use \f$L\f$ and \f$U\f$ for getting bounds on \f$x_T\f$.
 *
 * If \f$L + \mbox{infimum}(A_{kT}x_T) \geq b_k\f$, then the second constraint above is redundant.
 *
 * 2. ConvComb with clique-extension
 * Given two constraints
 * \f{eqnarray*}{
 *     A_{i\cdot} x \geq b_i \\
 *     A_{k\cdot} x \geq b_k \\
 *     \ell \leq x \leq u \\
 * \f}
 * this method determines promising values for \f$\lambda \in (0,1)\f$ and
 * applies feasibility-based bound tightening to the convex combinations
 *
 * \f$(\lambda A_{i\cdot} + (1 - \lambda) A_{k\cdot}) x \geq \lambda b_i + (1 - \lambda) b_k\f$.
 *
 * Additionally, cliques drawn from the SCIPcliqueTable are used
 * to further strengthen the above bound tightening. Full details can be found in
 *    - Belotti P. "Bound reduction using pairs of linear inequalities"
 *    - Chen W. et. al "Two-row and two-column mixed-integer presolve using hashing-based pairing methods"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_TWOROWBND_H__
#define __SCIP_PRESOL_TWOROWBND_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the tworowbnd presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresolTworowbnd(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
