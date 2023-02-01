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

/**@file   heur_rins.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LNS heuristic that combines the incumbent with the LP optimum
 * @author Timo Berthold
 *
 * RINS is a large neighborhood search improvement heuristic, i.e., it requires a known feasible solution. It solves a
 * sub-SCIP that is created by fixing variables which take the same value in the incumbent and the current node's LP
 * relaxation.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_RINS_H__
#define __SCIP_HEUR_RINS_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates RINS primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurRins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
