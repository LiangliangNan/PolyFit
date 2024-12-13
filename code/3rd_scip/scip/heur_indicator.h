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

/**@file   heur_indicator.h
 * @ingroup PRIMALHEURISTICS
 * @brief  handle partial solutions for linear problems with indicators and otherwise continuous variables
 * @author Marc Pfetsch
 *
 * For linear problems with indicators and otherwise continuous variables, the indicator constraint handler can produce
 * partial solutions, i.e., values for the indicator variables. This partial solution can be passed to this heuristic,
 * which then fixes these values and solves an LP. Additionally a local search for a better solution is added.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_INDICATOR_H__
#define __SCIP_HEUR_INDICATOR_H__

#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_heur.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the indicator primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurIndicator(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
 *
 * @{
 */

/** pass partial solution for indicator variables to heuristic */
SCIP_EXPORT
SCIP_RETCODE SCIPheurPassIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< indicator heuristic */
   int                   nindconss,          /**< number of indicator constraints */
   SCIP_CONS**           indconss,           /**< indicator constraints */
   SCIP_Bool*            solcand,            /**< values for indicator variables in partial solution */
   SCIP_Real             obj                 /**< objective of solution */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
