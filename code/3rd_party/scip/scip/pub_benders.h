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

/**@file   pub_benders.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_BENDERS_H__
#define __SCIP_PUB_BENDERS_H__

#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_benderscut.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "scip/type_stat.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBendersMethods
 *
 * @{
 */

/** compares two benderss w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPbendersComp);

/** comparison method for sorting benderss w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPbendersCompName);

/** gets user data of Benders' decomposition */
SCIP_EXPORT
SCIP_BENDERSDATA* SCIPbendersGetData(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** sets user data of Benders' decomposition; user has to free old data in advance! */
SCIP_EXPORT
void SCIPbendersSetData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSDATA*     bendersdata         /**< new Benders' decomposition user data */
   );

/** gets name of Benders' decomposition */
SCIP_EXPORT
const char* SCIPbendersGetName(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets description of Benders' decomposition */
SCIP_EXPORT
const char* SCIPbendersGetDesc(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets priority of Benders' decomposition */
SCIP_EXPORT
int SCIPbendersGetPriority(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of subproblems for the Benders' decomposition */
SCIP_EXPORT
int SCIPbendersGetNSubproblems(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   );

/** returns the SCIP instance for a given subproblem */
SCIP_EXPORT
SCIP* SCIPbendersSubproblem(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber          /**< the subproblem number */
   );

/** gets the number of times, the Bender' decomposition was called and tried to find a violated second stage constraint */
SCIP_EXPORT
int SCIPbendersGetNCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of optimality cuts found by the collection of Benders' decomposition subproblems */
SCIP_EXPORT
int SCIPbendersGetNCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of cuts found from the strengthening round */
SCIP_EXPORT
int SCIPbendersGetNStrengthenCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of calls to the strengthening round */
SCIP_EXPORT
int SCIPbendersGetNStrengthenCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets the number of calls to the strengthening round that fail */
SCIP_EXPORT
int SCIPbendersGetNStrengthenFails(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets time in seconds used in this Benders' decomposition for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPbendersGetSetupTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** gets execution time in seconds used in this Benders' decomposition */
SCIP_EXPORT
SCIP_Real SCIPbendersGetTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Is Benders' decomposition initialized? */
SCIP_EXPORT
SCIP_Bool SCIPbendersIsInitialized(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** returns whether the given Benders' decomposition is in use in the current problem */
SCIP_EXPORT
SCIP_Bool SCIPbendersIsActive(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   );

/** Returns whether only the convex relaxations will be checked in this solve loop
 *  when Benders' is used in the LNS heuristics, only the convex relaxations of the master/subproblems are checked,
 *  i.e. no integer cuts are generated. In this case, then Benders' decomposition is performed under the assumption
 *  that all subproblems are convex relaxations.
 */
SCIP_EXPORT
SCIP_Bool SCIPbendersOnlyCheckConvexRelax(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_Bool             subscipsoff         /**< flag indicating whether plugins using sub-SCIPs are deactivated */
   );

/** Are Benders' cuts generated from the LP solutions? */
SCIP_EXPORT
SCIP_Bool SCIPbendersCutLP(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Are Benders' cuts generated from the pseudo solutions? */
SCIP_EXPORT
SCIP_Bool SCIPbendersCutPseudo(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Are Benders' cuts generated from the relaxation solutions? */
SCIP_EXPORT
SCIP_Bool SCIPbendersCutRelaxation(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Should this Benders' use the auxiliary variables from the highest priority Benders'? */
SCIP_EXPORT
SCIP_Bool SCIPbendersShareAuxVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** sets the subproblem setup flag */
SCIP_EXPORT
void SCIPbendersSetSubproblemIsSetup(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             issetup             /**< flag to indicate whether the subproblem has been setup */
   );

/** returns the subproblem setup flag */
SCIP_EXPORT
SCIP_Bool SCIPbendersSubproblemIsSetup(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the auxiliary variable for the given subproblem */
SCIP_EXPORT
SCIP_VAR* SCIPbendersGetAuxiliaryVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns all auxiliary variables */
SCIP_EXPORT
SCIP_VAR** SCIPbendersGetAuxiliaryVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** stores the objective function value of the subproblem for use in cut generation */
SCIP_EXPORT
void SCIPbendersSetSubproblemObjval(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             objval              /**< the objective function value for the subproblem */
   );

/** returns the objective function value of the subproblem for use in cut generation */
SCIP_EXPORT
SCIP_Real SCIPbendersGetSubproblemObjval(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the number of cuts that have been added for storage */
SCIP_EXPORT
int SCIPbendersGetNStoredCuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition cut */
   );

/** returns the data for the cuts that have been added by the Benders' cut plugin */
SCIP_EXPORT
SCIP_RETCODE SCIPbendersGetStoredCutData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition cut */
   int                   cutidx,             /**< the index for the cut data that is requested */
   SCIP_VAR***           vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real**           vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   SCIP_Real*            rhs,                /**< the right hand side of the cut */
   int*                  nvars               /**< the number of variables with non-zero coefficients in the cut */
   );

/** returns the original problem data for the cuts that have been added by the Benders' cut plugin. The stored
 *  variables and values will populate the input vars and vals arrays. Thus, memory must be allocated for the vars and
 *  vals arrays
 */
SCIP_EXPORT
SCIP_RETCODE SCIPbendersGetStoredCutOrigData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition cut */
   int                   cutidx,             /**< the index for the cut data that is requested */
   SCIP_VAR***           vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real**           vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   SCIP_Real*            rhs,                /**< the right hand side of the cut */
   int*                  nvars,              /**< the number of variables with non-zero coefficients in the cut */
   int                   varssize            /**< the available slots in the array */
   );

/*
 * Public functions associated with Benders' cuts
 */

/** returns the Benders' cut of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_BENDERSCUT* SCIPfindBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   const char*           name                /**< name of Benderscut' decomposition */
   );


/** returns the array of currently available Benders' cuts; active Benders' decomposition are in the first slots of
 * the array
 */
SCIP_EXPORT
SCIP_BENDERSCUT** SCIPbendersGetBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );


/** returns the number of currently available Benders' cuts */
SCIP_EXPORT
int SCIPbendersGetNBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** sets the priority of a Benders' decomposition */
SCIP_EXPORT
SCIP_RETCODE SCIPbendersSetBenderscutPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' cut */
   int                   priority            /**< new priority of the Benders' decomposition */
   );

/** returns whether the solution has non-zero slack variables */
SCIP_EXPORT
SCIP_RETCODE SCIPbendersSolSlackVarsActive(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_Bool*            activeslack         /**< flag to indicate whether a slack variable is active */
   );

/** sets the subproblem type
 *
 * The subproblem types are:
 *    - Convex constraints with continuous variables
 *    - Convex constraints with discrete variables
 *    - Non-convex constraints with continuous variables
 *    - Non-convex constraints with discrete variables
 */
SCIP_EXPORT
void SCIPbendersSetSubproblemType(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSSUBTYPE   subprobtype         /**< the subproblem type */
   );

/** returns the type of the subproblem
 *
 *  This type is used to determine whether the duals of the problem can be used to generate cuts
 */
SCIP_EXPORT
SCIP_BENDERSSUBTYPE SCIPbendersGetSubproblemType(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** sets the flag indicating whether a subproblem is convex
 *
 *  It is possible that this can change during the solving process. One example is when the three-phase method is
 *  employed, where the first phase solves the convex relaxation of both the master and subproblems, the second phase
 *  reintroduces the integrality constraints to the master problem and the third phase then reintroduces integrality
 *  constraints to the subproblems.
 */
SCIP_EXPORT
void SCIPbendersSetSubproblemIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isconvex            /**< flag to indicate whether the subproblem is convex */
   );

/** returns whether the subproblem is convex
 *
 *  This means that the dual solution can be used to generate cuts.
 */
SCIP_EXPORT
SCIP_Bool SCIPbendersSubproblemIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the number of subproblems that are convex */
SCIP_EXPORT
int SCIPbendersGetNConvexSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** sets the flag indicating whether a subproblem contains non-linear constraints */
SCIP_EXPORT
void SCIPbendersSetSubproblemIsNonlinear(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isnonlinear         /**< flag to indicate whether the subproblem contains non-linear constraints */
   );

/** returns whether the subproblem contains non-linear constraints. */
SCIP_EXPORT
SCIP_Bool SCIPbendersSubproblemIsNonlinear(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the number of subproblems that contain non-linear constraints  */
SCIP_EXPORT
int SCIPbendersGetNNonlinearSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** sets the flag indicating whether the master problem contains non-linear constraints */
SCIP_EXPORT
void SCIPbendersSetMasterIsNonlinear(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_Bool             isnonlinear         /**< flag to indicate whether the subproblem contains non-linear constraints */
   );

/** returns whether the master problem contains non-linear constraints. */
SCIP_EXPORT
SCIP_Bool SCIPbendersMasterIsNonlinear(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** returns the flag indicating that Benders' decomposition is in a cut strengthening round */
SCIP_EXPORT
SCIP_Bool SCIPbendersInStrengthenRound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** solves the LP of the Benders' decomposition subproblem
 *
 *  This requires that the subproblem is in probing mode.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPbendersSolveSubproblemLP(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_STATUS*          solvestatus,        /**< status of subproblem solve */
   SCIP_Real*            objective           /**< optimal value of subproblem, if solved to optimality */
   );

/** solves the Benders' decomposition subproblem */
SCIP_EXPORT
SCIP_RETCODE SCIPbendersSolveSubproblemCIP(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_STATUS*          solvestatus,        /**< status of subproblem solve */
   SCIP_Bool             solvecip            /**< directly solve the CIP subproblem */
   );

/** returns the number of cuts that have been transferred from sub SCIPs to the master SCIP */
SCIP_EXPORT
int SCIPbendersGetNTransferredCuts(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   );

/** updates the lower bound for the subproblem. If the lower bound is not greater than the previously stored lowerbound,
 * then no update occurs.
 */
SCIP_EXPORT
void SCIPbendersUpdateSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             lowerbound          /**< the lower bound */
   );

/** returns the stored lower bound for the given subproblem */
SCIP_EXPORT
SCIP_Real SCIPbendersGetSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** sets the independent subproblem flag */
SCIP_EXPORT
void SCIPbendersSetSubproblemIsIndependent(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isindep             /**< flag to indicate whether the subproblem is independent */
   );

/** returns whether the subproblem is independent */
SCIP_EXPORT
SCIP_Bool SCIPbendersSubproblemIsIndependent(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** returns whether the subproblem is enabled, i.e. the subproblem is still solved in the solving loop. */
SCIP_EXPORT
SCIP_Bool SCIPbendersSubproblemIsEnabled(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
