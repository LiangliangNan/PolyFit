/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the reliable pseudo cost branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleRelpscost(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup BRANCHINGRULES
 *
 * @{
 */

/** execution reliability pseudo cost branching with the given branching candidates */
EXTERN
SCIP_RETCODE SCIPexecRelpscostBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            branchcands,        /**< branching candidates */
   SCIP_Real*            branchcandssol,     /**< solution value for the branching candidates */
   SCIP_Real*            branchcandsfrac,    /**< fractional part of the branching candidates */
   int                   nbranchcands,       /**< number of branching candidates */
   SCIP_Bool             executebranching,   /**< perform a branching step after probing */
   SCIP_RESULT*          result              /**< pointer to the result of the execution */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
