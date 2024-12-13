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

/**@file   sepa_minor.h
 * @ingroup SEPARATORS
 * @brief  principal minor separator
 * @author Benjamin Mueller
 *
 * This separator detects all principal minors of the matrix \f$ xx' \f$ for which all auxiliary variables \f$ X \f$
 * exist, i.e., two indices \f$ i \neq j \f$ such that \f$ X_{ii} \f$, \f$ X_{jj} \f$, and \f$ X_{ij} \f$ exist. Because
 * \f$ X - xx' \f$ is required to be positive semi-definite, it follows that the matrix
 *
 * \f[
 *    A(x,X) = \begin{bmatrix} 1 & x_i & x_j \\ x_i & X_{ii} & X_{ij} \\ x_j & X_{ij} & X_{jj} \end{bmatrix}
 * \f]
 *
 * is also required to be positive semi-definite. Let \f$ v \f$ be a negative eigenvector for \f$ A(x^*,X^*) \f$ in a
 * point \f$ (x^*,X^*) \f$, which implies that \f$ v' A(x^*,X^*) v < 0 \f$. To cut off \f$ (x^*,X^*) \f$, the separator
 * computes the globally valid linear inequality \f$ v' A(x,X) v \ge 0 \f$.
 *
 *
 * To identify which entries of the matrix X exist, we (the separator) iterate over the available nonlinear constraints.
 * For each constraint, we explore its expression and collect all nodes (expressions) of the form
 * - \f$x^2\f$
 * - \f$y \cdot z\f$
 *
 * Then, we go through the found bilinear terms \f$(yz)\f$ and if the corresponding \f$y^2\f$ and \f$z^2\f$ exist, then we have found
 * a minor.
 *
 * For circle packing instances, the minor cuts are not really helpful (see [Packing circles in a square: a theoretical
 * comparison of various convexification techniques](http://www.optimization-online.org/DB_HTML/2017/03/5911.html)).
 * Furthermore, the performance was negatively affected, thus circle packing constraint are identified and ignored in
 * the above algorithm. This behavior is controlled with the parameter "separating/minor/ignorepackingconss".
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_MINOR_H__
#define __SCIP_SEPA_MINOR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the minor separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaMinor(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup SEPARATORS
 *
 * @{
 */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
