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

/**@file   heur_coefdiving.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LP diving heuristic that chooses fixings w.r.t. the matrix coefficients
 * @author Tobias Achterberg
 *
 * Diving heuristic: Iteratively fixes some fractional variable and resolves the LP-relaxation, thereby simulating a
 * depth-first-search in the tree. Coefficient Diving chooses the variable with the fewest locking number in any
 * direction and rounds it into this direction. One-level backtracking is applied: If the LP gets infeasible, the last
 * fixing is undone, and the opposite fixing is tried. If this is infeasible, too, the procedure aborts.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_COEFDIVING_H__
#define __SCIP_HEUR_COEFDIVING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the coefdiving heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurCoefdiving(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
