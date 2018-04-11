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

/**@file   heur_proximity.h
 * @ingroup PRIMALHEURISTICS
 * @brief  improvement heuristic which uses an auxiliary objective instead of the original objective function which
 *         is itself added as a constraint to a sub-SCIP instance. The heuristic was presented by Matteo Fischetti
 *         and Michele Monaci
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_PROXIMITY_H__
#define __SCIP_HEUR_PROXIMITY_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the proximity primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurProximity(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
 *
 * @{
 */

/** main procedure of the proximity heuristic, creates and solves a sub-SCIP
 *
 *  @note the method can be applied in an iterative way, keeping the same subscip in between. If the @p freesubscip
 *        parameter is set to FALSE, the heuristic will keep the subscip data structures. Always set this parameter
 *        to TRUE, or call SCIPdeleteSubproblemProximity() afterwards
 */
EXTERN
SCIP_RETCODE SCIPapplyProximity(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                            */
   SCIP_RESULT*          result,             /**< result data structure                                               */
   SCIP_Real             minimprove,         /**< factor by which proximity should at least improve the incumbent     */
   SCIP_Longint          nnodes,             /**< node limit for the subproblem                                       */
   SCIP_Longint          nlpiters,           /**< LP iteration limit for the subproblem                               */
   SCIP_Longint*         nusednodes,         /**< pointer to store number of used nodes in subscip                    */
   SCIP_Longint*         nusedlpiters,       /**< pointer to store number of used LP iterations in subscip            */
   SCIP_Bool             freesubscip         /**< should the created sub-MIP be freed at the end of the method?       */
   );

/** frees the sub-MIP created by proximity */
EXTERN
SCIP_RETCODE SCIPdeleteSubproblemProximity(
   SCIP*                 scip                /** SCIP data structure */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
