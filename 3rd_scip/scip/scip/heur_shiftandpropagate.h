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

/**@file   heur_shiftandpropagate.h
 * @ingroup PRIMALHEURISTICS
 * @brief  preroot heuristic that alternatingly fixes variables and propagates domains
 * @author Timo Berthold
 * @author Gregor Hendel
 *
 * Preroot primal heuristic that fixes variables and propagates these fixings. In each step, the heuristic fixes a
 * variable such that the number of violated LP rows gets maximally reduced.  The fixing is then propagated to reduce
 * further variable domains.  In case that the domain propagation detects the infeasibility of the current partial
 * solution, the domain is reset to its previous state and the variable is postponed.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_SHIFTANDPROPAGATE_H__
#define __SCIP_HEUR_SHIFTANDPROPAGATE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the shiftandpropagate primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurShiftandpropagate(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
