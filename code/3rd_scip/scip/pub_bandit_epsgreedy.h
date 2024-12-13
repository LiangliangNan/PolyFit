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

/**@file   pub_bandit_epsgreedy.h
 * @ingroup PublicBanditMethods
 * @brief  public methods for the epsilon greedy bandit selector
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SRC_SCIP_PUB_BANDIT_EPSGREEDY_H_
#define SRC_SCIP_PUB_BANDIT_EPSGREEDY_H_


#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_bandit.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBanditMethods
 *
 * ## Epsilon greedy
 *
 * Epsilon greedy is a randomized algorithm for the multi-armed bandit problem.
 *
 * In every iteration, it either
 * selects an action uniformly at random with
 * probability \f$ \varepsilon_t\f$
 * or it greedily exploits the best action seen so far with
 * probability \f$ 1 - \varepsilon_t \f$.
 * In this implementation, \f$ \varepsilon_t \f$ decreases over time
 * (number of selections performed), controlled by the epsilon parameter.
 *
 * @{
 */

/** create and resets an epsilon greedy bandit algorithm */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateBanditEpsgreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         epsgreedy,          /**< pointer to store the epsilon greedy bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             eps,                /**< parameter to increase probability for exploration between all actions */
   SCIP_Bool             preferrecent,       /**< should the weights be updated in an exponentially decaying way? */
   SCIP_Real             decayfactor,        /**< the factor to reduce the weight of older observations if exponential decay is enabled */
   int                   avglim,             /**< nonnegative limit on observation number before the exponential decay starts,
                                               *  only relevant if exponential decay is enabled
                                               */
   int                   nactions,           /**< the number of possible actions */
   unsigned int          initseed            /**< initial seed for random number generation */
   );

/** get weights array of epsilon greedy bandit algorithm */
SCIP_EXPORT
SCIP_Real* SCIPgetWeightsEpsgreedy(
   SCIP_BANDIT*          epsgreedy           /**< epsilon greedy bandit algorithm */
   );

/** set epsilon parameter of epsilon greedy bandit algorithm */
SCIP_EXPORT
void SCIPsetEpsilonEpsgreedy(
   SCIP_BANDIT*          epsgreedy,          /**< epsilon greedy bandit algorithm */
   SCIP_Real             eps                 /**< parameter to increase probability for exploration between all actions */
   );

/** @} */



#ifdef __cplusplus
}
#endif

#endif
