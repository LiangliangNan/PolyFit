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

/**@file   heur_rens.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LNS heuristic that finds the optimal rounding to a given point
 * @author Timo Berthold
 *
 * RENS is a large neighborhood search start heuristic, i.e., unlike other LNS heuristics, it does not need a known
 * feasible solution. It solves a sub-SCIP that is created by fixing variables which take an integral value in a given
 * LP or NLP solution. For the remaining integer variables, the bounds get tightened to the two nearest integral values.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_RENS_H__
#define __SCIP_HEUR_RENS_H__

#include "scip/def.h"
#include "scip/type_heur.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates RENS primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurRens(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
 *
 * @{
 */

/** main procedure of the RNS heuristic, creates and solves a sub-SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPapplyRens(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                            */
   SCIP_RESULT*          result,             /**< result data structure                                               */
   SCIP_Real             minfixingrate,      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove,         /**< factor by which RENS should at least improve the incumbent          */
   SCIP_Longint          maxnodes,           /**< maximum number of  nodes for the subproblem                         */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                         */
   char                  startsol,           /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)]?       */
   SCIP_Bool             uselprows           /**< should subproblem be created out of the rows in the LP rows?        */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
