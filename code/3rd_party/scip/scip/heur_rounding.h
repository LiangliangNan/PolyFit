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

/**@file   heur_rounding.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LP rounding heuristic that tries to recover from intermediate infeasibilities
 * @author Tobias Achterberg
 *
 * Rounding heuristic that starts from an LP-feasible point and reduces the number of fractional variables by one in
 * each step. As long as no LP row is violated, the algorithm iterates over the fractional variables and applies a
 * rounding into the direction of fewer locks, updating the activities of the LP rows after each step.  If there is a
 * violated LP row, the heuristic will try to find a fractional variable that can be rounded in a direction such that
 * the violation of the constraint is decreased, using the number of up- and down-locks as a tie breaker.  If no
 * rounding can decrease the violation of the constraint, the procedure is aborted.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_ROUNDING_H__
#define __SCIP_HEUR_ROUNDING_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the rounding heuristic with infeasibility recovering and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurRounding(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
