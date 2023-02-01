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

/**@file   heur_twoopt.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Primal heuristic to improve incumbent solution by flipping pairs of variables
 * @author Timo Berthold
 * @author Gregor Hendel
 *
 * The Twoopt heuristic attempts to improve a feasible MIP solution by altering the solution values of pairs of
 * variables. Only variables which share a pre-defined ratio of LP rows are considered as pairs. Each step of the
 * heuristic consists of improving the objective value by shifting one variable, and then compensating the resulting
 * infeasibilities by shifting a second variable, without completely losing the objective improvement. Similarly to
 * Oneopt (see @ref heur_oneopt.h), pairs are processed in non-decreasing order of their impact on the objective.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_TWOOPT_H__
#define __SCIP_HEUR_TWOOPT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the twoopt primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurTwoopt(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
