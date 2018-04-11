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

/**@file   pub_bandit_ucb.h
 * @ingroup PublicBanditMethods
 * @brief  public methods for UCB bandit selection
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SRC_SCIP_PUB_BANDIT_UCB_H_
#define SRC_SCIP_PUB_BANDIT_UCB_H_


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup PublicBanditMethods
 *
 * ## Upper Confidence Bounds (UCB)
 *
 * UCB (Upper confidence bounds) is a deterministic
 * selection algorithm for the multi-armed bandit problem.
 * In every iteration, UCB selects the action that maximizes
 * a tradeoff between its performance in the past
 * and a variance term.
 * The influence of the variance (confidence width) can be
 * controlled by the parameter \f$ \alpha \f$.
 *
 * @{
 */


/** create and reset UCB bandit algorithm */
EXTERN
SCIP_RETCODE SCIPcreateBanditUcb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         ucb,                /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             alpha,              /**< parameter to increase confidence width */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial random number seed */
   );

/** returns the upper confidence bound of a selected action */
EXTERN
SCIP_Real SCIPgetConfidenceBoundUcb(
   SCIP_BANDIT*          ucb,                /**< UCB bandit algorithm */
   int                   action              /**< index of the queried action */
   );

/** return start permutation of the UCB bandit algorithm */
EXTERN
int* SCIPgetStartPermutationUcb(
   SCIP_BANDIT*          ucb                 /**< UCB bandit algorithm */
   );

/** @}*/


#ifdef __cplusplus
}
#endif

#endif
