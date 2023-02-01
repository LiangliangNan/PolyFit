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
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPbranchruleComp);

/** comparison method for sorting branching rules w.r.t. to their name */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPbranchruleCompName);

/** gets user data of branching rule */
EXTERN
SCIP_BRANCHRULEDATA* SCIPbranchruleGetData(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** sets user data of branching rule; user has to free old data in advance! */
EXTERN
void SCIPbranchruleSetData(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< new branching rule user data */
   );

/** gets name of branching rule */
EXTERN
const char* SCIPbranchruleGetName(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets description of branching rule */
EXTERN
const char* SCIPbranchruleGetDesc(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets priority of branching rule */
EXTERN
int SCIPbranchruleGetPriority(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
EXTERN
int SCIPbranchruleGetMaxdepth(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
EXTERN
SCIP_Real SCIPbranchruleGetMaxbounddist(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets time in seconds used in this branching rule for setting up for next stages */
EXTERN
SCIP_Real SCIPbranchruleGetSetupTime(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets time in seconds used in this branching rule */
EXTERN
SCIP_Real SCIPbranchruleGetTime(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of times, the branching rule was called on an LP solution */
EXTERN
SCIP_Longint SCIPbranchruleGetNLPCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of times, the branching rule was called on external candidates */
EXTERN
SCIP_Longint SCIPbranchruleGetNExternCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of times, the branching rule was called on a pseudo solution */
EXTERN
SCIP_Longint SCIPbranchruleGetNPseudoCalls(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of times, the branching rule detected a cutoff */
EXTERN
SCIP_Longint SCIPbranchruleGetNCutoffs(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of cuts, the branching rule separated */
EXTERN
SCIP_Longint SCIPbranchruleGetNCutsFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of constraints, the branching rule added to the respective local nodes (not counting constraints
 *  that were added to the child nodes as branching decisions)
 */
EXTERN
SCIP_Longint SCIPbranchruleGetNConssFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of domain reductions, the branching rule found */
EXTERN
SCIP_Longint SCIPbranchruleGetNDomredsFound(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets the total number of children, the branching rule created */
EXTERN
SCIP_Longint SCIPbranchruleGetNChildren(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/** is branching rule initialized? */
EXTERN
SCIP_Bool SCIPbranchruleIsInitialized(
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
