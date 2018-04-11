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

/**@file   heur_undercover.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Undercover primal heuristic for MINLPs
 * @author Timo Berthold
 * @author Ambros Gleixner
 *
 * The undercover heuristic is designed for mixed-integer nonlinear programs and tries to fix a subset of variables such
 * as to make each constraint linear or convex. For this purpose it solves a binary program to automatically determine
 * the minimum number of variable fixings necessary. As fixing values, we use values from the LP relaxation, the NLP
 * relaxation, or the incumbent solution.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_UNDERCOVER_H__
#define __SCIP_HEUR_UNDERCOVER_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the undercover primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurUndercover(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
 *
 * @{
 */

/** computes a minimal set of covering variables */
EXTERN
SCIP_RETCODE SCIPcomputeCoverUndercover(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  coversize,          /**< size of the computed cover */
   SCIP_VAR**            cover,              /**< pointer to store the variables (of the original SCIP) in the computed cover
                                              *   (should be ready to hold SCIPgetNVars(scip) entries) */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Real             objlimit,           /**< objective limit: upper bound on coversize */
   SCIP_Bool             globalbounds,       /**< should global bounds on variables be used instead of local bounds at focus node? */
   SCIP_Bool             onlyconvexify,      /**< should we only fix/dom.red. variables creating nonconvexity? */
   SCIP_Bool             coverbd,            /**< should bounddisjunction constraints be covered (or just copied)? */
   char                  coveringobj,        /**< objective function of the covering problem ('b'ranching status,
                                              *   influenced nonlinear 'c'onstraints/'t'erms, 'd'omain size, 'l'ocks,
                                              *   'm'in of up/down locks, 'u'nit penalties, constraint 'v'iolation) */
   SCIP_Bool*            success             /**< feasible cover found? */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
