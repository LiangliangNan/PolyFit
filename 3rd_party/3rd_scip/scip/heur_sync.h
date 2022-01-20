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

/**@file   heur_sync.h
 * @ingroup PRIMALHEURISTICS
 * @brief  primal heuristic that adds given solutions
 * @author Robert Lion Gottwald
 *
 * This heuristic takes solutions from somewhere else via the function SCIPheurSyncPassSol(). It
 * then tries to commit this solution. It is used by the concurrent solvers, when solutions are
 * communicated between solvers, but cannot directly submitted because SCIP might be in a stage where
 * this is not allowed.
 * If multiple solutions are passed it will keep the best N solutions depending on the parameter setting
 * "concsolvers/sync/maxnsols"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_SYNC_H__
#define __SCIP_HEUR_SYNC_H__

#include "scip/def.h"
#include "scip/type_sol.h"
#include "scip/type_scip.h"
#include "scip/type_heur.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the sync primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurSync(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
  *
  * @{
  */

/** pass solution to sync heuristic */
EXTERN
SCIP_RETCODE SCIPheurSyncPassSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< sync heuristic */
   SCIP_SOL*             sol                 /**< solution to be passed */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
