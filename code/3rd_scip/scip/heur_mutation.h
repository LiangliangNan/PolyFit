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

/**@file   heur_mutation.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LNS heuristic that tries to randomly mutate the incumbent solution
 * @author Timo Berthold
 *
 * Mutation is a large neighborhood search improvement heuristic that is inspired by genetic algorithms and requires a
 * known feasible solution. It randomly fixes variables to their value in the incumbent solution and solves a sub-SCIP,
 * consisting of the remaining variables. See also @ref heur_crossover.h
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_MUTATION_H__
#define __SCIP_HEUR_MUTATION_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the mutation primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurMutation(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
