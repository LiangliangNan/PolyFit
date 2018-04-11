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

/**@file   heur_simplerounding.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Simple and fast LP rounding heuristic
 * @author Tobias Achterberg
 *
 * Simple rounding is a very cheap heuristic that iterates over the set of fractional variables of an LP-feasible
 * point. It only performs roundings for variable that have zero up- or downlocks. Hence, they are guaranteed to keep
 * all constraints satisfied.  If all fractional variables can be rounded that way, the resulting solution will bve
 * integral and LP-feasible.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_SIMPLEROUNDING_H__
#define __SCIP_HEUR_SIMPLEROUNDING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the simple rounding heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurSimplerounding(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
