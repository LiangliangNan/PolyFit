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

/**@file   heur_zirounding.h
 * @ingroup PRIMALHEURISTICS
 * @brief  ZI Round primal heuristic
 * @author Gregor Hendel
 *
 * ZI Round (C. Wallace, Journal of Heuristics 2009) reduces the integer infeasibility of an LP solution step-by-step by
 * shifting fractional values towards integrality, but not necessarily rounding them.  For each integer variable with
 * fractional solution value, the heuristic calculates bounds for both possible rounding directions such that the
 * obtained solution stays LP-feasible. The solution value is then shifted by the corresponding bound into the direction
 * which reduces the fractionality most.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_ZIROUNDING_H__
#define __SCIP_HEUR_ZIROUNDING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the zirounding primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurZirounding(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
