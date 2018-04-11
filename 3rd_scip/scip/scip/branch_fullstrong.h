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

/**@file   branch_fullstrong.h
 * @ingroup BRANCHINGRULES
 * @brief  full strong LP branching rule
 * @author Tobias Achterberg
 *
 * The full strong branching rule applies strong branching to every fractional variable of the LP solution
 * at the current node of the branch-and-bound search. The rule selects the candidate
 * which will cause the highest gain of the dual bound in the created sub-tree among all branching variables.
 *
 * For calculating the gain, a look-ahead is performed by solving the child node LPs which will result
 * from branching on a variable.
 *
 * For a more mathematical description and a comparison between the strong branching rule and other branching rules
 * in SCIP, we refer to
 *
 * @par
 * Tobias Achterberg@n
 * Constraint Integer Programming@n
 * PhD Thesis, Technische Universität Berlin, 2007@n
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_FULLSTRONG_H__
#define __SCIP_BRANCH_FULLSTRONG_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the full strong LP branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleFullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup BRANCHINGRULES
 *
 * @{
 */

/**
 * Selects a variable from a set of candidates by strong branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 * @note The variables in the lpcands array must have a fractional value in the current LP solution
 */
EXTERN
SCIP_RETCODE SCIPselectVarStrongBranching(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP_VAR**            lpcands,            /**< branching candidates                                */
   SCIP_Real*            lpcandssol,         /**< solution values of the branching candidates         */
   SCIP_Real*            lpcandsfrac,        /**< fractional values of the branching candidates       */
   SCIP_Bool*            skipdown,           /**< should down branchings be skipped? */
   SCIP_Bool*            skipup,             /**< should up branchings be skipped? */
   int                   nlpcands,           /**< number of branching candidates                      */
   int                   npriolpcands,       /**< number of priority branching candidates             */
   int                   ncomplete,          /**< number of branching candidates without skip         */
   int*                  start,              /**< starting index in lpcands                           */
   int                   maxproprounds,      /**< maximum number of propagation rounds to be performed during strong
                                              *   branching before solving the LP (-1: no limit, -2: parameter settings) */
   SCIP_Bool             probingbounds,      /**< should valid bounds be identified in a probing-like fashion during
                                              *   strong branching (only with propagation)? */
   SCIP_Bool             forcestrongbranch,  /**< should strong branching be applied even if there is just a single candidate? */
   int*                  bestcand,           /**< best candidate for branching                        */
   SCIP_Real*            bestdown,           /**< objective value of the down branch for bestcand     */
   SCIP_Real*            bestup,             /**< objective value of the up branch for bestcand       */
   SCIP_Real*            bestscore,          /**< score for bestcand                                  */
   SCIP_Bool*            bestdownvalid,      /**< is bestdown a valid dual bound for the down branch? */
   SCIP_Bool*            bestupvalid,        /**< is bestup a valid dual bound for the up branch?     */
   SCIP_Real*            provedbound,        /**< proved dual bound for current subtree               */
   SCIP_RESULT*          result              /**< result pointer                                      */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
