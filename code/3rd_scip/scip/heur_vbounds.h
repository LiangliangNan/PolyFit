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

/**@file   heur_vbounds.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LNS heuristic uses the variable lower and upper bounds to determine the search neighborhood
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Jens Schulz
 * @author Gerald Gamrath
 *
 * More details about the heuristic can be found in "Structure-Based Primal Heuristics for Mixed Integer Programming"
 * by Gamrath, Berthold, Heinz, Winkler: http://link.springer.com/chapter/10.1007%2F978-4-431-55420-2_3
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_VBOUNDS_H__
#define __SCIP_HEUR_VBOUNDS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the variable bounds primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
