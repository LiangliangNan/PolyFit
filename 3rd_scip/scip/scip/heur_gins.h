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

/**@file   heur_gins.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LNS heuristic that tries to delimit the search region to a neighborhood in the constraint graph
 * @author Gregor Hendel
 *
 *
 * Graph Induced Neighborhood Search (GINS) is a Large Neighborhood Search Heuristic that attempts to improve
 * an incumbent solution by fixing a suitable percentage of integer variables to the incumbent and
 * solving the resulting, smaller and presumably easier sub-MIP.
 *
 * Its search neighborhoods are based on distances in a bipartite graph \f$G\f$ with the variables and constraints as nodes and
 * an edge between a variable and a constraint, if the variable is part of the constraint.
 * Given an integer \f$k\f$, the \f$k\f$-neighborhood of a variable \f$v\f$ in \f$G\f$ is the set of variables, whose nodes
 * are connected to \f$v\f$ by a path not longer than \f$2 \cdot k\f$. Intuitively, a judiciously chosen neighborhood size
 * allows to consider a local portion of the overall problem.
 *
 * An initial variable selection is made by randomly sampling different neighborhoods across the whole main problem.
 * The neighborhood that offers the largest potential for improvement is selected to become the local search neighborhood,
 * while all variables outside the neighborhood are fixed to their incumbent solution values.
 *
 * GINS also supports a rolling horizon approach, during which several local neighborhoods are considered
 * with increasing distance to the variable selected for the initial sub-problem. The rolling horizon approach ends
 * if no improvement could be found or a sufficient part of the problem component variables has been part of
 * at least one neighborhood.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_GINS_H__
#define __SCIP_HEUR_GINS_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the gins primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurGins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
