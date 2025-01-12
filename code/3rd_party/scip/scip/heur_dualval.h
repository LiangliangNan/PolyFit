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

/**@file   heur_dualval.h
 * @ingroup PRIMALHEURISTICS
 * @brief  primal heuristic that uses dualvalues for successive switching variable values
 * @author Tobias Buchwald
 *
 * This heuristic tries to find solutions by taking the LP or NLP, rounding solution values, fixing the variables to the
 * rounded values and then changing some of the values.To determine which variable is changed we give each variable a
 * ranking dependent on its dualvalue.  We work with a transformed problem that is always feasible and has objective = 0
 * iff the original problem is also feasible. Thus we cannot expect to find really good solutions.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_DUALVAL_H__
#define __SCIP_HEUR_DUALVAL_H__

#include "scip/def.h"
#include "scip/type_heur.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the dualVal primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurDualval(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
 *
 * @{
 */

/** main procedure of the dualval heuristic */
SCIP_EXPORT
SCIP_RETCODE SCIPapplyHeurDualval(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_RESULT*          result,             /**< pointer to store result of: did not run, solution found, no solution
                                              *   found, or fixing is infeasible (cutoff) */
   SCIP_SOL*             refpoint            /**< point to take fixation of discrete variables from; if NULL, then LP
                                              *   solution is used */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
