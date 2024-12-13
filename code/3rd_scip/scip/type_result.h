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

/**@file   type_result.h
 * @brief  result codes for SCIP callback methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_RESULT_H__
#define __SCIP_TYPE_RESULT_H__

#ifdef __cplusplus
extern "C" {
#endif

/** result codes for SCIP callback methods */
enum SCIP_Result
{
   SCIP_DIDNOTRUN   =   1,            /**< the method was not executed */
   SCIP_DELAYED     =   2,            /**< the method was not executed, but should be called again later */
   SCIP_DIDNOTFIND  =   3,            /**< the method was executed, but failed finding anything */
   SCIP_FEASIBLE    =   4,            /**< no infeasibility could be found */
   SCIP_INFEASIBLE  =   5,            /**< an infeasibility was detected */
   SCIP_UNBOUNDED   =   6,            /**< an unboundedness was detected */
   SCIP_CUTOFF      =   7,            /**< the current node is infeasible and can be cut off */
   SCIP_SEPARATED   =   8,            /**< the method added a cutting plane */
   SCIP_NEWROUND    =   9,            /**< the method added a cutting plane and a new separation round should immediately start */
   SCIP_REDUCEDDOM  =  10,            /**< the method reduced the domain of a variable */
   SCIP_CONSADDED   =  11,            /**< the method added a constraint */
   SCIP_CONSCHANGED =  12,            /**< the method changed a constraint */
   SCIP_BRANCHED    =  13,            /**< the method created a branching */
   SCIP_SOLVELP     =  14,            /**< the current node's LP must be solved */
   SCIP_FOUNDSOL    =  15,            /**< the method found a feasible primal solution */
   SCIP_SUSPENDED   =  16,            /**< the method interrupted its execution, but can continue if needed */
   SCIP_SUCCESS     =  17,            /**< the method was successfully executed */
   SCIP_DELAYNODE   =  18             /**< the processing of the branch-and-bound node should stopped and continued later */
};
typedef enum SCIP_Result SCIP_RESULT;           /**< result codes for SCIP callback methods */

#ifdef __cplusplus
}
#endif

#endif
