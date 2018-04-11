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


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the rounding heuristic with infeasibility recovering and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurRounding(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
