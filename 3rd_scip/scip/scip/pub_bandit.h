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

/**@file   pub_bandit.h
 * @ingroup PublicBanditMethods
 * @brief  public methods for bandit algorithms
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_BANDIT_H__
#define __SCIP_PUB_BANDIT_H__

#include "scip/def.h"
#include "scip/pub_bandit_epsgreedy.h"
#include "scip/pub_bandit_exp3.h"
#include "scip/pub_bandit_ucb.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBanditMethods
 *
 * @{
 */

/** select the next action */
EXTERN
SCIP_RETCODE SCIPbanditSelect(
   SCIP_BANDIT*          bandit,             /**< bandit algorithm data structure */
   int*                  action              /**< pointer to store the selected action */
   );

/** update the score of the selected action */
EXTERN
SCIP_RETCODE SCIPbanditUpdate(
   SCIP_BANDIT*          bandit,             /**< bandit algorithm data structure */
   int                   action,             /**< index of action for which the score should be updated */
   SCIP_Real             score               /**< observed gain of the i'th action */
   );

/** return the name of this bandit virtual function table */
EXTERN
const char* SCIPbanditvtableGetName(
   SCIP_BANDITVTABLE*    banditvtable        /**< virtual table for bandit algorithm */
   );

/** return the random number generator of a bandit algorithm */
EXTERN
SCIP_RANDNUMGEN* SCIPbanditGetRandnumgen(
   SCIP_BANDIT*          bandit              /**< bandit algorithm data structure */
   );

/** return number of actions of this bandit algorithm */
EXTERN
int SCIPbanditGetNActions(
   SCIP_BANDIT*          bandit              /**< bandit algorithm data structure */
   );

/* @} */


#ifdef __cplusplus
}
#endif

#endif
