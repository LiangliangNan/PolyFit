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

/**@file   heur_shifting.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LP rounding heuristic that tries to recover from intermediate infeasibilities and shifts continuous variables
 * @author Tobias Achterberg
 *
 * This heuristic is similar to the Rounding heuristic (see @ref heur_rounding.h), but it tries to continue in the case
 * that no rounding can decrease the violation of a linear constraint.  In this case, the value of a continuous variable
 * or an integer variable with integral value will be shifted in order to decrease the violation of the constraint.  To
 * avoid cycling, the procedure terminates after a certain number of non-improving shifts.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_SHIFTING_H__
#define __SCIP_HEUR_SHIFTING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the shifting heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurShifting(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
