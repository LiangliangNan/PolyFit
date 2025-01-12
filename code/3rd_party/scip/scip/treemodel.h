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

/**@file   treemodel.h
 * @ingroup PUBLICCOREAPI
 * @brief  Branching rules based on the Single-Variable-Branching (SVB) model
 * @author Daniel Anderson
 * @author Pierre Le Bodic
 *
 * The Single-Variable-Branching (SVB) model is a simplified model of
 * Branch & Bound trees, from which several nontrivial variable selection
 * rules arise. The Treemodel branching rule complements SCIP's hybrid
 * branching by suggesting improved branching variables given the current
 * pseudocosts and the current dual gap.
 *
 * Given a variable with dual bound changes (l, r) (both positive)
 * and an absolute gap G, the SVB model describes the tree that needs to be
 * built by branching on that same variable at every node until the value G
 * is reached at every leaf, starting from 0 at the root node.
 * If we do so for every variable, we can select the variable that produces
 * the smallest tree.
 * In the case where the gap is not known, then we can compute the growth rate
 * of the tree, which we call the ratio.
 * The ratio of a variable (l, r) is the factor by which the size of the tree
 * built using (l, r) that closes a gap G must be multiplied by to close a gap
 * G+1. This ratio is not constant for all gaps, but when G tends to infinity,
 * it converges to a fixed value we can compute numerically using a root finding
 * algorithm (e.g. Laguerre).
 * The ratio is used when the gap is too large (e.g. no primal bound known) or
 * to help approximate the size of the SVB tree for that variable.
 *
 * See the following publication for more detail:
 *
 * @par
 * Pierre Le Bodic and George Nemhauser@n
 * An abstract model for branching and its application to mixed integer programming@n
 * Mathematical Programming, 2017@n
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TREEMODEL_H__
#define __SCIP_TREEMODEL_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** initialises the Treemodel parameter data structure */
SCIP_EXPORT
SCIP_RETCODE SCIPtreemodelInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL**      treemodel           /**< Treemodel parameter data structure */
);

/** frees the Treemodel parameter data structure */
SCIP_EXPORT
SCIP_RETCODE SCIPtreemodelFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL**      treemodel           /**< Treemodel parameter data structure */
);

/** returns TRUE if the Treemodel branching rules are enabled */
SCIP_EXPORT
SCIP_Bool SCIPtreemodelIsEnabled(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel           /**< Treemodel parameter data structure */
);

/** apply the Treemodel branching rules to attempt to select a better
 *  branching candidate than the one selected by pseudocost branching */
SCIP_EXPORT
SCIP_RETCODE SCIPtreemodelSelectCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR**            branchcands,        /**< branching candidate storage */
   SCIP_Real*            mingains,           /**< minimum gain of rounding downwards or upwards */
   SCIP_Real*            maxgains,           /**< maximum gain of rounding downwards or upwards */
   SCIP_Real*            tiebreakerscore,    /**< scores to use for tie breaking */
   int                   nbranchcands,       /**< the number of branching candidates */
   int*                  bestcand            /**< the best branching candidate found before the call,
					          and the best candidate after the call (possibly the same) */
);

#ifdef __cplusplus
}
#endif

#endif
