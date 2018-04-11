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

/**@file   bandit_exp3.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for Exp.3 bandit algorithm
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BANDIT_EXP3_H__
#define __SCIP_BANDIT_EXP3_H__


#include "scip/scip.h"
#include "scip/bandit.h"

#ifdef __cplusplus
extern "C" {
#endif

/** include virtual function table for Exp.3 bandit algorithms */
extern
SCIP_RETCODE SCIPincludeBanditvtableExp3(
   SCIP*                 scip                /**< SCIP data structure */
   );

/*
 * callback methods for Exp.3 bandit algorithm to create a virtual function table
 */

/** callback to free bandit specific data structures */
extern
SCIP_DECL_BANDITFREE(SCIPbanditFreeExp3);

/** selection callback for bandit selector */
extern
SCIP_DECL_BANDITSELECT(SCIPbanditSelectExp3);

/** update callback for bandit algorithm */
extern
SCIP_DECL_BANDITUPDATE(SCIPbanditUpdateExp3);

/** reset callback for bandit algorithm */
extern
SCIP_DECL_BANDITRESET(SCIPbanditResetExp3);

/** direct bandit creation method for the core where no SCIP pointer is available */
extern
SCIP_RETCODE SCIPbanditCreateExp3(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_BANDITVTABLE*    vtable,             /**< virtual function table for callback functions of Exp.3 */
   SCIP_BANDIT**         exp3,               /**< pointer to store bandit algorithm */
   SCIP_Real*            priorities,         /**< nonnegative priorities for each action, or NULL if not needed */
   SCIP_Real             gammaparam,         /**< weight between uniform (gamma ~ 1) and weight driven (gamma ~ 0) probability distribution */
   SCIP_Real             beta,               /**< gain offset between 0 and 1 at every observation */
   int                   nactions,           /**< the positive number of actions for this bandit algorithm */
   unsigned int          initseed            /**< initial random seed */
   );

#ifdef __cplusplus
}
#endif

#endif
