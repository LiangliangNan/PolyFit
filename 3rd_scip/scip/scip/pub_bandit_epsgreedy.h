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

/**@file   pub_bandit_epsgreedy.h
 * @ingroup PublicBanditMethods
 * @brief  public methods for the epsilon greedy bandit selector
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SRC_SCIP_PUB_BANDIT_EPSGREEDY_H_
#define SRC_SCIP_PUB_BANDIT_EPSGREEDY_H_


#include "scip/scip.h"

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
EXTERN
SCIP_RETCODE SCIPcreateBanditEpsgreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         epsgreedy,          /**< pointer to store the epsilon greedy bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             eps,                /**< parameter to increase probability for exploration between all actions */
   int                   nactions,           /**< the number of possible actions */
   unsigned int          initseed            /**< initial seed for random number generation */
   );

/** get weights array of epsilon greedy bandit algorithm */
EXTERN
SCIP_Real* SCIPgetWeightsEpsgreedy(
   SCIP_BANDIT*          epsgreedy           /**< epsilon greedy bandit algorithm */
   );

/** set epsilon parameter of epsilon greedy bandit algorithm */
EXTERN
void SCIPsetEpsilonEpsgreedy(
   SCIP_BANDIT*          epsgreedy,          /**< epsilon greedy bandit algorithm */
   SCIP_Real             eps                 /**< parameter to increase probability for exploration between all actions */
   );

/* @} */



#ifdef __cplusplus
}
#endif

#endif
