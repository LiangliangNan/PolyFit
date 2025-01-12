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

/**@file   cutsel_hybrid.h
 * @ingroup CUTSELECTORS
 * @brief  hybrid cut selector
 * @author Leona Gottwald
 * @author Felipe Serrano
 * @author Mark Turner
 *
 * The hybrid cut selector scores cuts by using a weighted sum of the efficacy, directed cutoff distance, objective
 * parallelism, and integer support of the cuts. Afterwards, it selects the cuts using the score and filtering for
 * parallelism after selecting each cut.
 *
 * If a cut is given by \f$ a^T x \leq b \f$, then
 *  - the efficacy is defined as the distance between the LP solution and the hyperplane \f$ a^T x = b \f$;
 *  - the directed cutoff distance is defined as the distance between the LP solution and the hyperplane \f$ a^T x = b \f$
 *    restricted to the line segment joining the LP solution to the currently best primal solution; therefore, it is only
 *    defined when a primal solution is available;
 *  - the objective parallelism is how parallel the vector \f$ a \f$ is w.r.t. the objective function \f$ c \f$. That
 *    is, the objective parallelism is given by \f$ \frac{a^T c}{\|a\| \|c\|} \f$. Notice that the vectors are parallel
 *    when this formula returns 1;
 *  - the integer support of a cut is the ratio between the number of nonzero integer columns and the number of nonzero
 *    columns.
 *
 * These features of a cut can be recovered and/or computed with the functions @ref SCIPgetCutEfficacy(), @ref
 * SCIPgetCutLPSolCutoffDistance(), @ref SCIPgetRowObjParallelism(), and @ref SCIPgetRowNumIntCols(), @ref
 * SCIProwGetNNonz().
 *
 * The filtering step works as follows.
 * After computing the scores, these are divided in two groups: good scores and bad scores.  Any score larger or equal
 * to 90% of the largest score is considered a good score.
 *
 * First, the forced cuts --- cuts that are going to enter the LP no matter what --- are used to filter the non-forced
 * cuts. This means that for each forced cut, @p fcut, the parallelism between @p fcut and
 * every non-forced cut, @p cut, is computed (the parallelism between two cuts \f$ a^T x \leq b \f$ and \f$ d^T x \leq e\f$
 * is \f$ \frac{a^T d}{\|a\| \|d\|} \f$).
 * If the score of cut is good, then cut is dropped if its parallelism with @p fcut is larger or equal than the maximum
 * between \f$ \frac{1}{2} \f$ and 1 - minimum orthogonality.
 * If the score of cut is not good, then cut is dropped if its parallelism with @p fcut is larger or equal than 1 - minimum
 * orthogonality.
 *
 * @note The minimum orthogonality is a parameter that can be set, as well as the weights for the score.
 *
 * @note In the case of no primal solution, the weight assigned to the directed cutoff distance is transfered to the
 * efficacy.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTSEL_HYBRID_H__
#define __SCIP_CUTSEL_HYBRID_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the hybrid separator and includes it in SCIP
 *
 * @ingroup CutSelectorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeCutselHybrid(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CUTSELECTORS
 *
 * @{
 */

/** perform a cut selection algorithm for the given array of cuts
 *
 *  This is the selection method of the hybrid cut selector which uses a weighted sum of the
 *  efficacy, parallelism, directed cutoff distance, and the integral support.
 *  The input cuts array gets resorted s.t the selected cuts come first and the remaining
 *  ones are the end.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPselectCutsHybrid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_ROW**            forcedcuts,         /**< array with forced cuts */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator for tie-breaking, or NULL */
   SCIP_Real             goodscorefac,       /**< factor of best score among the given cuts to consider a cut good
                                              *   and filter with less strict settings of the maximum parallelism */
   SCIP_Real             badscorefac,        /**< factor of best score among the given cuts to consider a cut bad
                                              *   and discard it regardless of its parallelism to other cuts */
   SCIP_Real             goodmaxparall,      /**< maximum parallelism for good cuts */
   SCIP_Real             maxparall,          /**< maximum parallelism for non-good cuts */
   SCIP_Real             dircutoffdistweight,/**< weight of directed cutoff distance in cut score calculation */
   SCIP_Real             efficacyweight,     /**< weight of efficacy in cut score calculation */
   SCIP_Real             objparalweight,     /**< weight of objective parallelism in cut score calculation */
   SCIP_Real             intsupportweight,   /**< weight of integral support in cut score calculation */
   int                   ncuts,              /**< number of cuts in cuts array */
   int                   nforcedcuts,        /**< number of forced cuts */
   int                   maxselectedcuts,    /**< maximal number of cuts from cuts array to select */
   int*                  nselectedcuts       /**< pointer to return number of selected cuts from cuts array */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
