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

/**@file   prop_sync.h
 * @ingroup PROPAGATORS
 * @brief  propagator for applying global bound changes that were communicated by other
 *         concurrent solvers
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_SYNC_H__
#define __SCIP_PROP_SYNC_H__

#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_prop.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the sync propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropSync(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PROPAGATORS
  *
  * @{
  */

/** adds a boundchange to the sync propagator */
SCIP_EXPORT
SCIP_RETCODE SCIPpropSyncAddBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< sync propagator */
   SCIP_VAR*             var,                /**< variable for bound */
   SCIP_Real             val,                /**< value of bound */
   SCIP_BOUNDTYPE        bndtype             /**< type of bound */
   );

/** gives the total number of tightened bounds found by the sync propagator */
SCIP_EXPORT
SCIP_Longint SCIPpropSyncGetNTightenedBnds(
   SCIP_PROP*            prop                /**< sync propagator */
   );

/** gives the total number of tightened bounds for integer variables found by the sync propagator */
SCIP_EXPORT
SCIP_Longint SCIPpropSyncGetNTightenedIntBnds(
   SCIP_PROP*            prop                /**< sync propagator */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
