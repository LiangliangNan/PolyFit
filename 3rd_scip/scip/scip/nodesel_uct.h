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

/**@file   nodesel_uct.h
 * @ingroup NODESELECTORS
 * @brief  uct node selector which balances exploration and exploitation by considering node visits
 * @author Gregor Hendel
 *
 * the UCT node selection rule selects the next leaf according to a mixed score of the node's actual lower bound
 * and the number of times it has been visited so far compared to its parent node.
 *
 * The idea of UCT node selection for MIP appeared in:
 * Ashish Sabharwal and Horst Samulowitz
 * Guiding Combinatorial Optimization with UCT (2011)
 *
 * The authors adapted a game-tree exploration scheme called UCB to MIP trees. Starting from the root node as current node,
 * the algorithm selects the current node's child \f$N_i\f$ which maximizes the UCT score
 *
 * \f$ \mbox{score}(N_i) := -\mbox{estimate}_{N_i} + \mbox{weight} \cdot \frac{\mbox{visits}(\mbox{parent}(N_i))}{\mbox{visits}(N_i)}
 * \f$
 *
 * where \f$\mbox{estimate}\f$ is the node's lower bound normalized by the root lower bound, and \f$\mbox{visits}\f$
 * denotes the number of times a leaf in the subtree rooted at this node has been explored so far.
 *
 * The selected node in the sense of the SCIP node selection is the leaf reached by the above criterion.
 *
 * The authors suggest that this node selection rule is particularly useful at the beginning of the solving process, but
 * to switch to a different node selection after a number of nodes has been explored to reduce computational overhead.
 * Our implementation uses only information available from the original SCIP tree which does not support the
 * forward path mechanism needed for the most efficient node selection. Instead, the algorithm selects the next leaf
 * by looping over all leaves and comparing the best leaf found so far with the next one. Two leaves l_1, l_2 are compared
 * by following their paths back upwards until their deepest common ancestor \f$a\f$ is reached, together with the two
 * children of \f$a\f$ representing the two paths to l_1, l_2. The leaf represented by the child of \f$a\f$
 * with higher UCT score is a candidate for the next selected leaf.
 *
 * The node selector features several parameters:
 *
 * the nodelimit delimits the number of explored nodes before UCT selection is turned off
 * the weight parameter changes the relevance of the visits quotient in the UCT score (see above score formula)
 * useestimate determines whether the node's estimate or lower bound is taken as estimate
 *
 * @note It should be avoided to switch to uct node selection after the branch and bound process has begun because
 *       the central UCT score information how often a path was taken is not collected if UCT is inactive. A safe use of
 *       UCT is to switch it on before SCIP starts optimization.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NODESEL_UCT_H__
#define __SCIP_NODESEL_UCT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the uct node selector and includes it in SCIP
 *
 *  @ingroup NodeSelectorIncludes
 */
extern
SCIP_RETCODE SCIPincludeNodeselUct(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
