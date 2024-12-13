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

/**@file   pub_reopt.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for reoptimization
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_REOPT_H__
#define __SCIP_PUB_REOPT_H__

#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_reopt.h"
#include "scip/type_var.h"

#ifdef NDEBUG
#include "scip/struct_reopt.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * ReoptNode methods
 */

/** returns the number of bound changes stored in the reoptnode */
SCIP_EXPORT
int SCIPreoptnodeGetNVars(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   );

/** returns the number of bound changes at the node stored at ID id */
SCIP_EXPORT
int SCIPreoptnodeGetNConss(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   );

/** returns the number of stored bound changes based on dual information in the reopttree at ID id */
SCIP_EXPORT
int SCIPreoptnodeGetNDualBoundChgs(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   );

/** returns the number of child nodes of @p reoptnode */
SCIP_EXPORT
int SCIPreoptnodeGetNChildren(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimizzation tree */
   );

/* return the lower bound stored at @p ID id */
SCIP_EXPORT
SCIP_Real SCIPreoptnodeGetLowerbound(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   );

/** returns the type of the @p reoptnode */
SCIP_EXPORT
SCIP_REOPTTYPE SCIPreoptnodeGetType(
   SCIP_REOPTNODE*       reoptnode           /**< node of the reoptimization tree */
   );

/** create the constraint which splits the node stored at ID id on the basis of the stored dual information. */
SCIP_EXPORT
void SCIPreoptnodeGetSplitCons(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array to store the variables of the constraint */
   SCIP_Real*            vals,               /**< array to store the coefficients of the variables */
   REOPT_CONSTYPE*       constype,           /**< type of the constraint */
   int                   conssize,           /**< size of the arrays */
   int*                  nvars               /**< pointer to store the size of the constraints */
   );

/** returns all added constraints at ID id */
SCIP_EXPORT
void SCIPreoptnodeGetConss(
   SCIP_REOPTNODE*       reoptnode,          /**< reoptimization data structure */
   SCIP_VAR***           vars,               /**< 2-dim array of variables */
   SCIP_Real**           bounds,             /**< 2-dim array of bounds */
   SCIP_BOUNDTYPE**      boundtypes,         /**< 2-dim array of boundtypes */
   int                   mem,                /**< allocated memory for constraints */
   int*                  nconss,             /**< pointer to store the number of constraints */
   int*                  nvars               /**< pointer to store the number of variables */
   );

/** set the parent id */
SCIP_EXPORT
void SCIPreoptnodeSetParentID(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   unsigned int          parentid            /**< id of the parent node */
   );

/*
 * Reopt methods
 */

/** returns the number of global restarts */
SCIP_EXPORT
int SCIPreoptGetNRestartsGlobal(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** returns the number of local restarts in the current run */
int SCIPreoptGetNRestartsLocal(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   );

/** returns the number of local restarts over all runs */
int SCIPreoptGetNTotalRestartsLocal(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   );

/** returns the number of iteration with the first global restarts */
SCIP_EXPORT
int SCIPreoptGetFirstRestarts(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   );

/** returns the number of iteration with the last global restarts */
SCIP_EXPORT
int SCIPreoptGetLastRestarts(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   );

/** returns the number of nodes providing an improving feasible LP solution in the current run */
SCIP_EXPORT
int SCIPreoptGetNFeasNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   );

/** returns the number of nodes providing an improving feasible LP solution over all runs */
SCIP_EXPORT
int SCIPreoptGetNTotalFeasNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   );

/** returns the number of nodes that exceeded the cutoff bound in the current run */
SCIP_EXPORT
int SCIPreoptGetNPrunedNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   );

/** returns the number of nodes that exceeded the cutoff bound over all runs */
SCIP_EXPORT
int SCIPreoptGetNTotalPrunedNodes(
   SCIP_REOPT*           reopt               /**< reoptimization data structure */
   );

/** returns the number of reoptimized nodes that were cut off in the current run */
SCIP_EXPORT
int SCIPreoptGetNCutoffReoptnodes(
   SCIP_REOPT*           reopt               /*< reoptimization data structure */
   );

/** returns the number of reoptimized nodes that were cut off over all runs */
SCIP_EXPORT
int SCIPreoptGetNTotalCutoffReoptnodes(
   SCIP_REOPT*           reopt               /*< reoptimization data structure */
   );

/** returns the number of stored nodes with an infeasible LP in the current run */
SCIP_EXPORT
int SCIPreoptGetNInfNodes(
   SCIP_REOPT*           reopt               /*< reoptimization data structure */
   );

/** returns the number of stored nodes with an infeasible LP over all runs */
SCIP_EXPORT
int SCIPreoptGetNTotalInfNodes(
   SCIP_REOPT*           reopt               /*< reoptimization data structure */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPreoptnodeGetNVars(reoptnode)          (reoptnode->nvars)
#define SCIPreoptnodeGetNConss(reoptnode)         (reoptnode->nconss)
#define SCIPreoptnodeGetNDualBoundChgs(reoptnode) (reoptnode->dualconscur->nvars)
#define SCIPreoptnodeGetNChildren(reoptnode)      (reoptnode->nchilds)
#define SCIPreoptnodeGetLowerbound(reoptnode)     (reoptnode->lowerbound)
#define SCIPreoptnodeGetType(reoptnode)           (reoptnode->reopttype)

#define SCIPreoptGetNRestartsGlobal(reopt)        (reopt->nglbrestarts)
#define SCIPreoptGetNRestartsLocal(reopt)         (reopt->nlocrestarts)
#define SCIPreoptGetNTotalRestartsLocal(reopt)    (reopt->ntotallocrestarts)
#define SCIPreoptGetFirstRestarts(reopt)          (reopt->firstrestart)
#define SCIPreoptGetLastRestarts(reopt)           (reopt->lastrestart)
#define SCIPreoptGetNFeasNodes(reopt)             (reopt->reopttree->nfeasnodes)
#define SCIPreoptGetNTotalFeasNodes(reopt)        (reopt->reopttree->ntotalfeasnodes)
#define SCIPreoptGetNPrunedNodes(reopt)           (reopt->reopttree->nprunednodes)
#define SCIPreoptGetNTotalPrunedNodes(reopt)      (reopt->reopttree->ntotalprunednodes)
#define SCIPreoptGetNCutoffReoptnodes(reopt)      (reopt->reopttree->ncutoffreoptnodes)
#define SCIPreoptGetNTotalCutoffReoptnodes(reopt) (reopt->reopttree->ntotalcutoffreoptnodes)
#define SCIPreoptGetNInfNodes(reopt)              (reopt->reopttree->ninfsubtrees)
#define SCIPreoptGetNTotalInfNodes(reopt)         (reopt->reopttree->ntotalinfnodes)

#endif

#ifdef __cplusplus
}
#endif

#endif
