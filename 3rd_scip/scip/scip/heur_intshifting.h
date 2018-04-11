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

/**@file   heur_intshifting.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LP rounding heuristic that tries to recover from intermediate infeasibilities, shifts integer variables, and
 *         solves a final LP to calculate feasible values for continuous variables
 * @author Tobias Achterberg
 *
 * This heuristic is similar to the Shifting heuristic (see @ref heur_shifting.h), but it ignores continuous variables
 * during the shifting phase and solves a final LP to find feasible (and optimal w.r.t. the integer fixings) values for
 * the continuous variables.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_INTSHIFTING_H__
#define __SCIP_HEUR_INTSHIFTING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the intshifting heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurIntshifting(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
