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

/**@file   event_estim.h
 * @ingroup EVENTS
 * @brief  event handler for tree size estimation and restarts
 *
 * This event handler plugin provides different methods for approximating the current fraction of the search
 * that has already been completed and for estimating the total tree size at completion.
 * It can trigger restarts of the current run if the current run seems hopeless.
 *
 * For details about the available approximations of search completion, please see
 *
 * Anderson, Hendel, Le Bodic, Pfetsch
 * Estimating The Size of Branch-and-Bound Trees
 * under preparation
 *
 * This code is a largely enriched version of a code that was used for clairvoyant restarts, see
 *
 * Anderson, Hendel, Le Bodic, Viernickel
 * Clairvoyant Restarts in Branch-and-Bound Search Using Online Tree-Size Estimation
 * AAAI-19: Proceedings of the Thirty-Third AAAI Conference on Artificial Intelligence, 2018
 *
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_ESTIM_H__
#define __SCIP_EVENT_ESTIM_H__


#include "scip/type_scip.h"
#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for tree size estimation */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeEventHdlrEstim(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* return an estimation of the final tree size */
SCIP_EXPORT
SCIP_Real SCIPgetTreesizeEstimation(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
