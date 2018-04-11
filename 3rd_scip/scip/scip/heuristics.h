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

/**@file   heuristics.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods commonly used by primal heuristics
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEURISTICS_H__
#define __SCIP_HEURISTICS_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_heur.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@defgroup PublicSpecialHeuristicMethods Special Methods
 * @ingroup PublicHeuristicMethods
 * @brief  methods commonly used by primal heuristics
 *
 * @{
 */

/** performs a diving within the limits of the diveset parameters
 *
 *  This method performs a diving according to the settings defined by the diving settings @p diveset; Contrary to the
 *  name, SCIP enters probing mode (not diving mode) and dives along a path into the tree. Domain propagation
 *  is applied at every node in the tree, whereas probing LPs might be solved less frequently.
 *
 *  Starting from the current LP solution, the algorithm selects candidates which maximize the
 *  score defined by the @p diveset and whose solution value has not yet been rendered infeasible by propagation,
 *  and propagates the bound change on this candidate.
 *
 *  The algorithm iteratively selects the the next (unfixed) candidate in the list, until either enough domain changes
 *  or the resolve frequency of the LP trigger an LP resolve (and hence, the set of potential candidates changes),
 *  or the last node is proven to be infeasible. It optionally backtracks and tries the
 *  other branching direction.
 *
 *  After the set of remaining candidates is empty or the targeted depth is reached, the node LP is
 *  solved, and the old candidates are replaced by the new LP candidates.
 *
 *  @see heur_guideddiving.c for an example implementation of a dive set controlling the diving algorithm.
 *
 *  @note the node from where the algorithm is called is checked for a basic LP solution. If the solution
 *        is non-basic, e.g., when barrier without crossover is used, the method returns without performing a dive.
 *
 *  @note currently, when multiple diving heuristics call this method and solve an LP at the same node, only the first
 *        call will be executed, @see SCIPgetLastDiveNode().
 */
EXTERN
SCIP_RETCODE SCIPperformGenericDivingAlgorithm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< settings for diving */
   SCIP_SOL*             worksol,            /**< non-NULL working solution */
   SCIP_HEUR*            heur,               /**< the calling primal heuristic */
   SCIP_RESULT*          result,             /**< SCIP result pointer */
   SCIP_Bool             nodeinfeasible      /**< is the current node known to be infeasible? */
   );

/** get a sub-SCIP copy of the transformed problem */
EXTERN
SCIP_RETCODE SCIPcopyLargeNeighborhoodSearch(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP used by the heuristic */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables */
   const char*           suffix,             /**< suffix for the problem name */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             uselprows,          /**< should the linear relaxation of the problem defined by LP rows be copied? */
   SCIP_Bool             copycuts,           /**< should cuts be copied (only if uselprows == FALSE) */
   SCIP_Bool*            success,            /**< was the copying successful? */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
