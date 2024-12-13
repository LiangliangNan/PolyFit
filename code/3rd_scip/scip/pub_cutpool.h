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

/**@file   pub_cutpool.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CUTPOOL_H__
#define __SCIP_PUB_CUTPOOL_H__


#include "scip/def.h"
#include "scip/type_cutpool.h"
#include "scip/type_lp.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCutMethods
 *
 * @{
 */

/** gets the row of the cut */
SCIP_EXPORT
SCIP_ROW* SCIPcutGetRow(
   SCIP_CUT*             cut                 /**< cut */
   );

/** gets the age of the cut: the number of consecutive cut pool separation rounds where the cut was neither in the LP nor violated */
SCIP_EXPORT
int SCIPcutGetAge(
   SCIP_CUT*             cut                 /**< cut */
   );

/** returns the ratio of LPs where the row belonging to this cut was active in an LP solution, i.e.
 *  where the age of its row has not been increased
 *
 *  @see SCIPcutGetAge() to get the age of a cut
 */
SCIP_EXPORT
SCIP_Real SCIPcutGetLPActivityQuot(
   SCIP_CUT*             cut                 /**< cut */
   );

/** gets array of cuts in the cut pool */
SCIP_EXPORT
SCIP_CUT** SCIPcutpoolGetCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get number of cuts in the cut pool */
SCIP_EXPORT
int SCIPcutpoolGetNCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get maximum number of cuts that were stored in the cut pool at the same time */
SCIP_EXPORT
SCIP_Longint SCIPcutpoolGetMaxNCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** gets time in seconds used for separating cuts from the pool */
SCIP_EXPORT
SCIP_Real SCIPcutpoolGetTime(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get number of times the cut pool was separated */
SCIP_EXPORT
SCIP_Longint SCIPcutpoolGetNCalls(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get number of times the cut pool was separated at the root */
SCIP_EXPORT
SCIP_Longint SCIPcutpoolGetNRootCalls(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get total number of cuts that were separated from the cut pool */
SCIP_EXPORT
SCIP_Longint SCIPcutpoolGetNCutsFound(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get total number of cuts that were added from the cut pool to sepastore */
SCIP_EXPORT
SCIP_Longint SCIPcutpoolGetNCutsAdded(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );
/** @} */

#ifdef __cplusplus
}
#endif

#endif
