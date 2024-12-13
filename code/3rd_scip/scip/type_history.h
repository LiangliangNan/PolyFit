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

/**@file   type_history.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_HISTORY_H__
#define __SCIP_TYPE_HISTORY_H__

#ifdef __cplusplus
extern "C" {
#endif

/** branching direction for branching on variables */
enum SCIP_BranchDir
{
   SCIP_BRANCHDIR_DOWNWARDS = 0,        /**< downwards branching: decreasing upper bound  */
   SCIP_BRANCHDIR_UPWARDS   = 1,        /**< upwards branching: increasing lower bound    */
   SCIP_BRANCHDIR_FIXED     = 2,        /**< fixed branching: both bounds changed         */
   SCIP_BRANCHDIR_AUTO      = 3         /**< automatic setting for choosing bound changes */
};
typedef enum SCIP_BranchDir SCIP_BRANCHDIR;       /**< branching direction for branching on variables */

typedef struct SCIP_History SCIP_HISTORY;         /**< branching and inference history information for single variable */

/** Value history data structure
 *
 *  branching and inference history informations for single variable dependent on the domain value
 *
 *  - \ref SCIP_VALUEHISTORY "List of all available methods"
 */
typedef struct SCIP_ValueHistory SCIP_VALUEHISTORY;


#ifdef __cplusplus
}
#endif

#endif
