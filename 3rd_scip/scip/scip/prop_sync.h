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

/**@file   prop_sync.h
 * @ingroup PROPAGATORS
 * @brief  propagator for applying global bound changes that were communicated by other
 *         concurrent solvers
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_SYNC_H__
#define __SCIP_PROP_SYNC_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the sync propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePropSync(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PROPAGATORS
  *
  * @{
  */

/** adds a boundchange to the sync propagator */
EXTERN
SCIP_RETCODE SCIPpropSyncAddBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< sync propagator */
   SCIP_VAR*             var,                /**< variable for bound */
   SCIP_Real             val,                /**< value of bound */
   SCIP_BOUNDTYPE        bndtype             /**< type of bound */
   );

/** gives the total number of tightened bounds found by the sync propagator */
EXTERN
SCIP_Longint SCIPpropSyncGetNTightenedBnds(
   SCIP_PROP*            prop                /**< sync propagator */
   );

/** gives the total number of tightened bounds for integer variables found by the sync propagator */
EXTERN
SCIP_Longint SCIPpropSyncGetNTightenedIntBnds(
   SCIP_PROP*            prop                /**< sync propagator */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
