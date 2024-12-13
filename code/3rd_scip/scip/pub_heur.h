/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_heur.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_HEUR_H__
#define __SCIP_PUB_HEUR_H__

#include "scip/def.h"
#include "scip/type_heur.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_timing.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicHeuristicMethods
 *
 * @{
 */



/** compares two heuristics w.r.t. to their delay positions and priorities */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPheurComp);

/** compares two heuristics w.r.t. to their priority values */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPheurCompPriority);

/** comparison method for sorting heuristics w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPheurCompName);

/** gets user data of primal heuristic */
SCIP_EXPORT
SCIP_HEURDATA* SCIPheurGetData(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** sets user data of primal heuristic; user has to free old data in advance! */
SCIP_EXPORT
void SCIPheurSetData(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< new primal heuristic user data */
   );

/** gets name of primal heuristic */
SCIP_EXPORT
const char* SCIPheurGetName(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets description of primal heuristic */
SCIP_EXPORT
const char* SCIPheurGetDesc(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets display character of primal heuristic */
SCIP_EXPORT
char SCIPheurGetDispchar(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** returns the timing mask of the heuristic */
SCIP_EXPORT
SCIP_HEURTIMING SCIPheurGetTimingmask(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** sets new timing mask for heuristic */
SCIP_EXPORT
void SCIPheurSetTimingmask(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_HEURTIMING       timingmask          /**< new timing mask of heuristic */
   );

/** does the heuristic use a secondary SCIP instance? */
SCIP_EXPORT
SCIP_Bool SCIPheurUsesSubscip(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets priority of primal heuristic */
SCIP_EXPORT
int SCIPheurGetPriority(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets frequency of primal heuristic */
SCIP_EXPORT
int SCIPheurGetFreq(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** sets frequency of primal heuristic */
SCIP_EXPORT
void SCIPheurSetFreq(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   freq                /**< new frequency of heuristic */
   );

/** gets frequency offset of primal heuristic */
SCIP_EXPORT
int SCIPheurGetFreqofs(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets maximal depth level for calling primal heuristic (returns -1, if no depth limit exists) */
SCIP_EXPORT
int SCIPheurGetMaxdepth(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets the number of times, the heuristic was called and tried to find a solution */
SCIP_EXPORT
SCIP_Longint SCIPheurGetNCalls(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets the number of primal feasible solutions found by this heuristic */
SCIP_EXPORT
SCIP_Longint SCIPheurGetNSolsFound(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets the number of new best primal feasible solutions found by this heuristic */
SCIP_EXPORT
SCIP_Longint SCIPheurGetNBestSolsFound(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** is primal heuristic initialized? */
SCIP_EXPORT
SCIP_Bool SCIPheurIsInitialized(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets time in seconds used in this heuristic for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPheurGetSetupTime(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** gets time in seconds used in this heuristic */
SCIP_EXPORT
SCIP_Real SCIPheurGetTime(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** returns array of divesets of this primal heuristic, or NULL if it has no divesets */
SCIP_EXPORT
SCIP_DIVESET** SCIPheurGetDivesets(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** returns the number of divesets of this primal heuristic */
SCIP_EXPORT
int SCIPheurGetNDivesets(
   SCIP_HEUR*            heur                /**< primal heuristic */
   );

/** @} */

/** get the heuristic to which this diving setting belongs */
SCIP_EXPORT
SCIP_HEUR* SCIPdivesetGetHeur(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the working solution of this dive set */
SCIP_EXPORT
SCIP_SOL* SCIPdivesetGetWorkSolution(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** set the working solution for this dive set */
SCIP_EXPORT
void SCIPdivesetSetWorkSolution(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_SOL*             sol                 /**< new working solution for this dive set, or NULL */
   );

/**@addtogroup PublicDivesetMethods
 *
 * @{
 */

/** get the name of the dive set */
SCIP_EXPORT
const char* SCIPdivesetGetName(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the minimum relative depth of the diving settings */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetMinRelDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the maximum relative depth of the diving settings */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetMaxRelDepth(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the number of successful runs of the diving settings */
SCIP_EXPORT
SCIP_Longint SCIPdivesetGetSolSuccess(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the number of calls to this dive set */
SCIP_EXPORT
int SCIPdivesetGetNCalls(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the number of calls successfully terminated at a feasible leaf node */
SCIP_EXPORT
int SCIPdivesetGetNSolutionCalls(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the minimum depth reached by this dive set */
SCIP_EXPORT
int SCIPdivesetGetMinDepth(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the maximum depth reached by this dive set */
SCIP_EXPORT
int SCIPdivesetGetMaxDepth(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the average depth this dive set reached during execution */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetAvgDepth(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the minimum depth at which this dive set found a solution */
SCIP_EXPORT
int SCIPdivesetGetMinSolutionDepth(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the maximum depth at which this dive set found a solution */
SCIP_EXPORT
int SCIPdivesetGetMaxSolutionDepth(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the average depth at which this dive set found a solution */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetAvgSolutionDepth(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the total number of LP iterations used by this dive set */
SCIP_EXPORT
SCIP_Longint SCIPdivesetGetNLPIterations(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the total number of probing nodes used by this dive set */
SCIP_EXPORT
SCIP_Longint SCIPdivesetGetNProbingNodes(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the total number of backtracks performed by this dive set */
SCIP_EXPORT
SCIP_Longint SCIPdivesetGetNBacktracks(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the total number of conflicts found by this dive set */
SCIP_EXPORT
SCIP_Longint SCIPdivesetGetNConflicts(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the total number of solutions (leaf and rounded solutions) found by the dive set */
SCIP_EXPORT
SCIP_Longint SCIPdivesetGetNSols(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the maximum LP iterations quotient of the diving settings */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetMaxLPIterQuot(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the maximum LP iterations offset of the diving settings */
SCIP_EXPORT
int SCIPdivesetGetMaxLPIterOffset(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the maximum upper bound quotient parameter of the diving settings if no solution is available */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetUbQuotNoSol(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the average quotient parameter of the diving settings if no solution is available */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetAvgQuotNoSol(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the maximum upper bound quotient parameter of the diving settings if an incumbent solution exists */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetUbQuot(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** get the average upper bound quotient parameter of the diving settings if an incumbent solution exists */
SCIP_EXPORT
SCIP_Real SCIPdivesetGetAvgQuot(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** should backtracking be applied? */
SCIP_EXPORT
SCIP_Bool SCIPdivesetUseBacktrack(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** returns the LP solve frequency for diving LPs (0: dynamically based on number of intermediate domain reductions) */
SCIP_EXPORT
int SCIPdivesetGetLPSolveFreq(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** returns the domain reduction quotient for triggering an immediate resolve of the diving LP (0.0: always resolve)*/
SCIP_EXPORT
SCIP_Real SCIPdivesetGetLPResolveDomChgQuot(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** should only LP branching candidates be considered instead of the slower but
 *  more general constraint handler diving variable selection?
 */
SCIP_EXPORT
SCIP_Bool SCIPdivesetUseOnlyLPBranchcands(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** returns TRUE if dive set supports diving of the specified type */
SCIP_EXPORT
SCIP_Bool SCIPdivesetSupportsType(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_DIVETYPE         divetype            /**< bit mask that represents the supported dive types by this dive set */
   );

/** returns the random number generator of this \p diveset for tie-breaking */
SCIP_EXPORT
SCIP_RANDNUMGEN* SCIPdivesetGetRandnumgen(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** is this dive set publicly available (ie., can be used by other primal heuristics?) */
SCIP_EXPORT
SCIP_Bool SCIPdivesetIsPublic(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** @} */

/**@defgroup PublicVariableGraphMethods Public Variable Graph Methods
 * @ingroup MiscellaneousMethods
 *
 * @brief  methods to create a variable graph and perform breadth-first search
 *
 * @{
 */

/** Perform breadth-first (BFS) search on the variable constraint graph.
 *
 *  The result of the algorithm is that the \p distances array contains the correct distances for
 *  every variable from the start variables. The distance of a variable can then be accessed through its
 *  problem index (calling SCIPvarGetProbindex()).
 *  Hence, The method assumes that the length of \p distances is at least
 *  SCIPgetNVars().
 *  Variables that are not connected through constraints to the start variables have a distance of -1.
 *
 *  Limits can be provided to further restrict the breadth-first search. If a distance limit is given,
 *  the search will be performed until the first variable at this distance is popped from the queue, i.e.,
 *  all variables with a distance < maxdistance have been labeled by the search.
 *  If a variable limit is given, the search stops after it completes the distance level at which
 *  the limit was reached. Hence, more variables may be actually labeled.
 *  The start variables are accounted for those variable limits.
 *
 *  If no variable variable constraint graph is provided, the method will create one and free it at the end
 *  This is useful for a single use of the variable constraint graph. For several consecutive uses,
 *  it is advised to create a variable constraint graph via SCIPvariableGraphCreate().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPvariablegraphBreadthFirst(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VGRAPH*          vargraph,           /**< pointer to the variable graph, or NULL to let the function create a local graph */
   SCIP_VAR**            startvars,          /**< array of start variables to calculate distance from */
   int                   nstartvars,         /**< number of starting variables, at least 1 */
   int*                  distances,          /**< array to keep distance in vargraph from start variables for every variable */
   int                   maxdistance,        /**< maximum distance >= 0 from start variable (INT_MAX for complete BFS) */
   int                   maxvars,            /**< maximum number of variables to compute distance for */
   int                   maxbinintvars       /**< maximum number of binary or integer variables to compute distance for */
   );

/** initialization method of variable graph data structure */
SCIP_EXPORT
SCIP_RETCODE SCIPvariableGraphCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VGRAPH**         vargraph,           /**< pointer to the variable graph */
   SCIP_Bool             relaxdenseconss,    /**< should dense constraints (at least as dense as \p density) be
                                              *   ignored by connectivity graph? */
   SCIP_Real             relaxdensity,       /**< density (with respect to number of variables) to relax constraint from graph */
   int*                  nrelaxedconstraints  /**< pointer to store the number of constraints that were relaxed, or NULL if not needed */
   );

/** deinitialization method of variable graph data structure */
SCIP_EXPORT
void SCIPvariableGraphFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VGRAPH**         vargraph            /**< pointer to the variable graph */
   );

/** @} */


#ifdef __cplusplus
}
#endif

#endif
