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

/**@file   heur_locks.h
 * @ingroup PRIMALHEURISTICS
 * @brief  locks primal heuristic
 * @author Michael Winkler
 * @author Gerald Gamrath
 *
 * The locks heuristic is a start heuristic that first tries to fix all binary variables, then solves the resulting LP
 * and tries to round the solution and finally solves a sub-MIP on the remaining problem if the LP solution could not be
 * rounded. The fixing works as follows: First, all variables are sorted by their total number of rounding locks (up-
 * and down-locks summed up). Then, looking at the variable with the highest number of locks first, the variable is
 * fixed to the bound where there are fewer locks (in case of ties, the bound which is better w.r.t. the objective
 * function). This fix is propagated and the activities of all LP rows are updated. If any LP row becomes redundant
 * w.r.t. the updated bounds, we adjust the rounding locks.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LOCKS_H__
#define __SCIP_HEUR_LOCKS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the locks primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurLocks(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** apply fix-and-propagate scheme based on variable locks
 *
 *  @note probing mode of SCIP needs to be enabled before
 */
EXTERN
SCIP_RETCODE SCIPapplyLockFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            cutoff,             /**< pointer to store if a cutoff was detected */
   SCIP_Bool*            allrowsfulfilled    /**< pointer to store if all rows became redundant */
   );

#ifdef __cplusplus
}
#endif

#endif
