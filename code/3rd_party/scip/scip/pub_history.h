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

/**@file   pub_history.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for branching and inference history structure
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_HISTORY_H__
#define __SCIP_PUB_HISTORY_H__

#include "scip/def.h"
#include "scip/type_history.h"

#ifdef NDEBUG
#include "scip/struct_history.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** gets the conflict score of the history entry */
SCIP_EXPORT
SCIP_Real SCIPhistoryGetVSIDS(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   );

/** gets the average conflict length of the history entry */
SCIP_EXPORT
SCIP_Real SCIPhistoryGetAvgConflictlength(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   );

/** get number of cutoffs counter */
SCIP_EXPORT
SCIP_Real SCIPhistoryGetCutoffSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** get number of inferences counter */
SCIP_EXPORT
SCIP_Real SCIPhistoryGetInferenceSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** return the number of (domain) values for which a history exists */
SCIP_EXPORT
int SCIPvaluehistoryGetNValues(
   SCIP_VALUEHISTORY*    valuehistory        /**< value based history */
   );

/** return the array containing the histories for the individual (domain) values */
SCIP_EXPORT
SCIP_HISTORY** SCIPvaluehistoryGetHistories(
   SCIP_VALUEHISTORY*    valuehistory        /**< value based history */
   );

/** return the array containing the (domain) values for which a history exists */
SCIP_EXPORT
SCIP_Real* SCIPvaluehistoryGetValues(
   SCIP_VALUEHISTORY*    valuehistory        /**< value based history */
   );

#ifdef NDEBUG

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPhistoryGetVSIDS(history,dir)   ((history)->vsids[dir])
#define SCIPhistoryGetAvgConflictlength(history,dir) ((history)->conflengthsum[dir] > 0.0 \
      ? (SCIP_Real)(history)->nactiveconflicts[dir]/(SCIP_Real)(history)->conflengthsum[dir] : 0.0)
#define SCIPhistoryGetCutoffSum(history,dir)        ((history)->cutoffsum[dir])
#define SCIPhistoryGetInferenceSum(history,dir)     ((history)->inferencesum[dir])
#define SCIPvaluehistoryGetNValues(valuehistory)     (valuehistory)->nvalues
#define SCIPvaluehistoryGetHistories(valuehistory)      (valuehistory)->histories
#define SCIPvaluehistoryGetValues(valuehistory)      (valuehistory)->values

#endif


#ifdef __cplusplus
}
#endif

#endif
