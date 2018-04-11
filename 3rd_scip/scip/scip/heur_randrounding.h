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

/**@file   heur_randrounding.h
 * @ingroup PRIMALHEURISTICS
 * @brief  randomized LP rounding heuristic which also generates conflicts via an auxiliary probing tree
 * @author Gregor Hendel
 *
 * Randomized LP rounding uses a random variable from a uniform distribution
 * over [0,1] to determine whether the fractional LP value x should be rounded
 * up with probability x - floor(x) or down with probability ceil(x) - x.
 *
 * This implementation uses domain propagation techniques to tighten the variable domains after every
 * rounding step.
 *
 * @see: The most relevant publication is Raghavan & Thompson,
 *       "Randomized rounding: A technique for provably good algorithms and algorithmic proofs",
 *        Combinatorica 7 (4): 365–374
 *        1987
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_RANDROUNDING_H__
#define __SCIP_HEUR_RANDROUNDING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the rand rounding heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurRandrounding(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
