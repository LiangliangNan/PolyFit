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

/**@file   heur_oneopt.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Improvement heuristic that alters single variable values
 * @author Timo Berthold
 *
 * Oneopt is a straightforward improvement heuristic: given a feasible MIP solution, the value of
 * an integer variable x<sub>j</sub> can be decreased for c<sub>j</sub> > 0 or increased for c<sub>j</sub> < 0
 * if the resulting solution is still feasible.  If more than one variable can be shifted, they are sorted by
 * non-decreasing impact on the objective and sequentially shifted until no more improvements can be
 * obtained. 
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_ONEOPT_H__
#define __SCIP_HEUR_ONEOPT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the oneopt primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurOneopt(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
