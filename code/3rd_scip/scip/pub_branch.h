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

/**@file   pub_branch.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for branching rules
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_BRANCH_H__
#define __SCIP_PUB_BRANCH_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_branch.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBranchRuleMethods
 *
 * @{
 */

/** compares two branching rules w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPbranchruleComp);

/** comparison method for sorting branching rules w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPbranchruleCompName);

/** gets user data of branching rule */
SCIP_EXPORT
SCIP_BRANCHRULEDATA* SCIPbranchruleGetData(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** sets user data of branching rule; user has to free old data in advance! */
SCIP_EXPORT
void SCIPbranchruleSetData(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< new branching rule user data */
   );

/** gets name of branching rule */
SCIP_EXPORT
const char* SCIPbranchruleGetName(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets description of branching rule */
SCIP_EXPORT
const char* SCIPbranchruleGetDesc(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets priority of branching rule */
SCIP_EXPORT
int SCIPbranchruleGetPriority(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
SCIP_EXPORT
int SCIPbranchruleGetMaxdepth(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
SCIP_EXPORT
SCIP_Real SCIPbranchruleGetMaxbounddist(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets time in seconds used in this branching rule for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPbranchruleGetSetupTime(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets time in seconds used in this branching rule */
SCIP_EXPORT
SCIP_Real SCIPbranchruleGetTime(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of times, the branching rule was called on an LP solution */
SCIP_EXPORT
SCIP_Longint SCIPbranchruleGetNLPCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of times, the branching rule was called on external candidates */
SCIP_EXPORT
SCIP_Longint SCIPbranchruleGetNExternCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of times, the branching rule was called on a pseudo solution */
SCIP_EXPORT
SCIP_Longint SCIPbranchruleGetNPseudoCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of times, the branching rule detected a cutoff */
SCIP_EXPORT
SCIP_Longint SCIPbranchruleGetNCutoffs(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of cuts, the branching rule separated */
SCIP_EXPORT
SCIP_Longint SCIPbranchruleGetNCutsFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of constraints, the branching rule added to the respective local nodes (not counting constraints
 *  that were added to the child nodes as branching decisions)
 */
SCIP_EXPORT
SCIP_Longint SCIPbranchruleGetNConssFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of domain reductions, the branching rule found */
SCIP_EXPORT
SCIP_Longint SCIPbranchruleGetNDomredsFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of children, the branching rule created */
SCIP_EXPORT
SCIP_Longint SCIPbranchruleGetNChildren(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** is branching rule initialized? */
SCIP_EXPORT
SCIP_Bool SCIPbranchruleIsInitialized(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
