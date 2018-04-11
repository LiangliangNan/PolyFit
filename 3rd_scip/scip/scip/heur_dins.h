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

/**@file   heur_dins.h
 * @ingroup PRIMALHEURISTICS
 * @brief  DINS primal heuristic
 * @author Robert Waniek
 *
 * DINS combines the ideas of RINS (see @ref heur_rins.h) and Local Branching (see @ref heur_localbranching.h). It
 * defines the neighborhood by introducing a distance function between the incumbent solution and the optimum of the LP
 * relaxation. When applied during a branch-and-bound search, it further takes into account how variables change their
 * values at different nodes of the tree.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_DINS_H__
#define __SCIP_HEUR_DINS_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the DINS primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurDins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
