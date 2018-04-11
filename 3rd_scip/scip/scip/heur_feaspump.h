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

/**@file   heur_feaspump.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Objective Feasibility Pump 2.0
 * @author Timo Berthold
 * @author Domenico Salvagnin
 *
 * The fundamental idea of the Feasibility Pump is to construct two sequences of points which hopefully converge to a
 * feasible solution. One sequence consists of LP-feasiblepoints, the other one of integer feasible points.  They are
 * produced by alternately rounding an LP-feasible point and solvng an LP that finds a point on the LP polyhedron which
 * is closest to the rounded, integral point (w.r.t. Manhattan distance).
 *
 * The version implemented in SCIP supports using an Objective Feasibility Pump that uses a convex combination of the
 * Manhattan distance and the original LP objective for reoptimization. It further features Feasibility Pump 2.0
 * capabilities, hence propagating the fixings for a faster convergence.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_FEASPUMP_H__
#define __SCIP_HEUR_FEASPUMP_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the feaspump primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurFeaspump(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
