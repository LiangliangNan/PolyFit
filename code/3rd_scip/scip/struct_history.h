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

/**@file   struct_history.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_HISTORY_H__
#define __SCIP_STRUCT_HISTORY_H__


#include "scip/def.h"
#include "scip/type_history.h"

#ifdef __cplusplus
extern "C" {
#endif

/** branching and inference history information for single variable independent of the domain value */
struct SCIP_History
{
   SCIP_Real             pscostcount[2];     /**< nr of (partial) summands in down/upwards pseudo costs (may be fractional) */
   SCIP_Real             pscostweightedmean[2]; /**< weighted mean of (partial) pseudo cost values for down/upwards branching */
   SCIP_Real             pscostvariance[2];  /**< weighted variance of (partial) pseudo cost history for down/upwards branching */
   SCIP_Real             vsids[2];           /**< degree of how often the variable was reason for a conflict */
   SCIP_Real             conflengthsum[2];   /**< overall length of all active conflicts for which the variable gave reason */
   SCIP_Real             inferencesum[2];    /**< degree of how often branching on the variable lead to inference of another bound */
   SCIP_Real             cutoffsum[2];       /**< degree of how often branching on the variable lead to an infeasible sub problem */
   SCIP_Real             ratio;              /**< Most recently computed value of the unpowered ratio (phi^l) in the Treemodel rules */
   SCIP_Real             balance;            /**< Most recently value of r/l used to compute a ratio in the Treemodel rules */
   SCIP_Bool             ratiovalid;         /**< Whether the ratio value is valid for the Treemodel rules */
   SCIP_Longint          nactiveconflicts[2];/**< number of active conflicts for which the variable gave reason */
   SCIP_Longint          nbranchings[2];     /**< nr of times, the variable changed its bounds due to branching */
   SCIP_Longint          branchdepthsum[2];  /**< sum of depth levels, at which the branching bound changes took place */
};

/** branching and inference history information for single variable dependent on the domain value */
struct SCIP_ValueHistory
{
   SCIP_HISTORY**        histories;          /**< branching and inference history information for individual domain values */
   SCIP_Real*            values;             /**< (sorted) array containing the domain values */
   int                   nvalues;            /**< number of different domain values */
   int                   sizevalues;         /**< size of the both arrays */
};

#ifdef __cplusplus
}
#endif

#endif
