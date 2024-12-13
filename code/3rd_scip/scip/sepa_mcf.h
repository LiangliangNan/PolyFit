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

/**@file   sepa_mcf.h
 * @ingroup SEPARATORS
 * @brief  multi-commodity-flow network cut separator
 * @author Tobias Achterberg
 * @author Christian Raack
 *
 * We try to identify a multi-commodity flow structure in the LP relaxation of the
 * following type:
 *
 *  (1)  sum_{a in delta^+(v)} f_a^k  - sum_{a in delta^-(v)} f_a^k  <=  -d_v^k   for all v in V and k in K
 *  (2)  sum_{k in K} f_a^k - c_a x_a                                <=  0        for all a in A
 *
 * Constraints (1) are flow conservation constraints, which say that for each commodity k and node v the
 * outflow (delta^+(v)) minus the inflow (delta^-(v)) of a node v must not exceed the negative of the demand of
 * node v in commodity k. To say it the other way around, inflow minus outflow must be at least equal to the demand.
 * Constraints (2) are the arc capacity constraints, which say that the sum of all flow over an arc a must not
 * exceed its capacity c_a x_a, with x being a binary or integer variable.
 * c_a x_a does not need to be a single product of a capacity and an integer variable; we also accept general scalar
 * products.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_MCF_H__
#define __SCIP_SEPA_MCF_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#ifdef __cplusplus
extern "C" {
#endif

/** creates the mcf separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaMcf(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
