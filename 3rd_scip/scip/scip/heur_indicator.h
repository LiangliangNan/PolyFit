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

/**@file   heur_indicator.h
 * @ingroup PRIMALHEURISTICS
 * @brief  handle partial solutions for linear problems with indicators and otherwise continuous variables
 * @author Marc Pfetsch
 *
 * For linear problems with indicators and otherwise continuous variables, the indicator constraint handler can produce
 * partial solutions, i.e., values for the indicator variables. This partial solution can be passed to this heuristic,
 * which then fixes these values and solves an LP. Additionally a local search for a better solution is added.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_INDICATOR_H__
#define __SCIP_HEUR_INDICATOR_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the indicator primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurIndicator(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
 *
 * @{
 */

/** pass partial solution for indicator variables to heuristic */
EXTERN
SCIP_RETCODE SCIPheurPassIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< indicator heuristic */
   int                   nindconss,          /**< number of indicator constraints */
   SCIP_CONS**           indconss,           /**< indicator constraints */
   SCIP_Bool*            solcand,            /**< values for indicator variables in partial solution */
   SCIP_Real             obj                 /**< objective of solution */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
