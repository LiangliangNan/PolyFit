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

/**@file   branch_multaggr.h
 * @ingroup BRANCHINGRULES
 * @brief  fullstrong branching on fractional and multi-aggregated variables
 * @author Anna Melchiori
 * @author Gerald Gamrath
 *
 * This branching rule uses all fractional binary and integer variables as candidates,
 * as well as fractional multiaggregated binary and integer variables. Although not
 * directly contained in the presolved problem anymore, the multi-aggregation provides
 * an affine linear sum of integer variables, on which branching can be performed.
 *
 * For more details, see
 * G.Gamrath, A.Melchiori, T.Berthold, A.M.Gleixner, D.Salvagnin: Branching on Multi-aggregated Variables
 * (http://dx.doi.org/10.1007/978-3-319-18008-3_10)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_MULTAGGR_H__
#define __SCIP_BRANCH_MULTAGGR_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the multi-aggregated branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBranchruleMultAggr(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
