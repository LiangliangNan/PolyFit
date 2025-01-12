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

/**@file   branch_vanillafullstrong.h
 * @ingroup BRANCHINGRULES
 * @brief  vanilla full strong LP branching rule
 * @author Tobias Achterberg
 * @author Maxime Gasse
 *
 * The vanilla full strong branching rule is a purged implementation of full strong branching, for academic purposes.
 * It implements full strong branching with the following specific features:
 * - no cutoff or domain reduction: only branching.
 * - idempotent (optional): leave SCIP, as much as possible, in the same state before / after the strong branching
 *   calls. Basically, do not update any statistic.
 * - donotbranch (optional): do no perform branching. So that the brancher can be called as an oracle only (on which
 *   variable would you branch ? But do not branch please).
 * - scoreall (optional): continue scoring variables, even if infeasibility is detected along the way.
 * - collectscores (optional): store the candidate scores from the last call, which can then be retrieved by calling
 *   SCIPgetVanillafullstrongData().
 * - integralcands (optional): get candidates from SCIPgetPseudoBranchCands() instead of SCIPgetLPBranchCands(), i.e.,
 *   consider all non-fixed variables as branching candidates, not only fractional ones.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_VANILLAFULLSTRONG_H__
#define __SCIP_BRANCH_VANILLAFULLSTRONG_H__


#include "scip/def.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the vanilla full strong branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBranchruleVanillafullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** recovers candidate variables and their scores from last vanilla full strong branching call */
SCIP_EXPORT
SCIP_RETCODE SCIPgetVanillafullstrongData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           cands,              /**< pointer to store candidate variables; or NULL */
   SCIP_Real**           candscores,         /**< pointer to store candidate scores; or NULL */
   int*                  ncands,             /**< pointer to store number of candidates; or NULL */
   int*                  npriocands,         /**< pointer to store number of priority candidates; or NULL */
   int*                  bestcand            /**< pointer to store best branching candidate; or NULL */
   );


#ifdef __cplusplus
}
#endif

#endif
