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

/**@file   sepa_aggregation.h
 * @ingroup SEPARATORS
 * @brief  flow cover and complemented mixed integer rounding cuts separator (Marchand's version)
 * @author Leona Gottwald
 * @author Kati Wolter
 * @author Tobias Achterberg
 *
 * For an overview see:
 *
 * Marchand, H., & Wolsey, L. A. (2001).@n
 * Aggregation and mixed integer rounding to solve MIPs.@n
 * Operations research, 49(3), 363-371.
 *
 * Some remarks:
 * - In general, continuous variables are less prefered than integer variables, since their cut
 *   coefficient is worse.
 * - We seek for aggregations that project out continuous variables that are far away from their bound,
 *   since if it is at its bound then it doesn't contribute to the violation
 * - These aggregations are also useful for the flowcover separation, so after building an aggregation
 *   we try to generate a MIR cut and a flowcover cut.
 * - We only keep the best cut.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_CMIR_H__
#define __SCIP_SEPA_CMIR_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the aggregation separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaAggregation(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
