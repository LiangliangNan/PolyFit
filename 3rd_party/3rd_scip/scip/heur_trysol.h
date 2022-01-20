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

/**@file   heur_trysol.h
 * @ingroup PRIMALHEURISTICS
 * @brief  primal heuristic that tries a given solution
 * @author Marc Pfetsch
 *
 * This heuristic takes a solution from somewhere else via the function SCIPheurPassSolTrySol(). It
 * then tries to commit this solution. It is mainly used by cons_indicator, which tries to correct a
 * given solution, but cannot directly submit this solution, because it is a constraint handler and
 * not a heuristic.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_TRYSOL_H__
#define __SCIP_HEUR_TRYSOL_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the trysol primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurTrySol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
 *
 * @{
 */

/** pass solution to trysol heuristic */
EXTERN
SCIP_RETCODE SCIPheurPassSolTrySol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< trysol heuristic */
   SCIP_SOL*             sol                 /**< solution to be passed */
   );

/** pass solution to trysol heuristic which just gets added (without checking feasibility */
EXTERN
SCIP_RETCODE SCIPheurPassSolAddSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< trysol heuristic */
   SCIP_SOL*             sol                 /**< solution to be passed */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
