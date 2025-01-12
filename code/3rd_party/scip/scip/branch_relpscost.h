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

/**@file   branch_relpscost.h
 * @ingroup BRANCHINGRULES
 * @brief  reliable pseudo costs branching rule
 * @author Tobias Achterberg
 *
 * The reliable pseudo costs branching rule uses the notion of pseudo costs to measure the expected
 * gain in the dual bound when branching on a particular variable.
 * The pseudo cost information is collected during the branch-and-bound search in the same manner as for
 * the pseudo costs branching rule.
 *
 * The reliable pseudo costs branching rule, however, uses a limited number of look-ahead LP-iterations
 * at the beginning of the search in order to obtain better pseudo cost estimates and make branching decisions in a
 * sense more "reliable" at an early stage of the search,
 * at the price of a higher computational cost at the beginning of the search.
 *
 * For a more mathematical description and a comparison between the reliable pseudo costs rule and other branching rules
 * in SCIP, we refer to
 *
 * @par
 * Tobias Achterberg@n
 * Constraint Integer Programming@n
 * PhD Thesis, Technische Universität Berlin, 2007@n
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_RELPSCOST_H__
#define __SCIP_BRANCH_RELPSCOST_H__


#include "scip/def.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the reliable pseudo cost branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBranchruleRelpscost(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup BRANCHINGRULES
 *
 * @{
 */

/** execution reliability pseudo cost branching with the given branching candidates */
SCIP_EXPORT
SCIP_RETCODE SCIPexecRelpscostBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            branchcands,        /**< branching candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   SCIP_Real*            branchcandsfrac,    /**< fractional part of the branching candidates */
   int                   nbranchcands,       /**< number of branching candidates */
   SCIP_Bool             executebranching,   /**< perform a branching step after probing */
   SCIP_RESULT*          result              /**< pointer to the result of the execution */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
