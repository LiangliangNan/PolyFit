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

/**@file   heur_nlpdiving.h
 * @ingroup PRIMALHEURISTICS
 * @brief  NLP diving heuristic that chooses fixings w.r.t. the fractionalities
 * @author Timo Berthold
 * @author Stefan Vigerske
 *
 * Diving heuristic: Iteratively fixes some fractional variable and resolves the NLP-relaxation, thereby simulating a
 * depth-first-search in the tree. Fractional Diving chooses the variable with the highest fractionality and rounds it to the
 * nearest integer. One-level backtracking is applied: If the NLP gets infeasible, the last fixing is undone, and the
 * opposite fixing is tried. If this is infeasible, too, the procedure aborts.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_NLPDIVING_H__
#define __SCIP_HEUR_NLPDIVING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the fracdiving heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurNlpdiving(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
