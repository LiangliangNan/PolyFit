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

/**@file   benders.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERS_H__
#define __SCIP_BENDERS_H__

#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_benderscut.h"
#include "scip/type_dcmp.h"
#include "scip/type_message.h"
#include "scip/type_misc.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_set.h"
#include "scip/type_sol.h"
#include "scip/type_stat.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif


/** copies the given Benders' decomposition to a new scip */
SCIP_RETCODE SCIPbendersCopyInclude(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             sourceset,          /**< SCIP_SET of SCIP to copy from */
   SCIP_SET*             targetset,          /**< SCIP_SET of SCIP to copy to */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables; must not be NULL */
   SCIP_Bool             copysubproblems,    /**< must the subproblems be copied with the Benders' decomposition copy */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   );

/** creates a Benders' decomposition */
SCIP_RETCODE SCIPbendersCreate(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_Bool             shareauxvars,       /**< should this Benders' use the highest priority Benders aux vars */
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre)),/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre)),/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve)),/**< called prior to the subproblem solving loop */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex)),/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve)),/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub)),/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   );

/** calls destructor and frees memory of Benders' decomposition */
SCIP_RETCODE SCIPbendersFree(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes Benders' decomposition */
SCIP_RETCODE SCIPbendersInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of Benders' decomposition */
SCIP_RETCODE SCIPbendersExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs the Benders' decomposition that the presolving process is being started */
SCIP_RETCODE SCIPbendersInitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** informs the Benders' decomposition that the presolving process has completed */
SCIP_RETCODE SCIPbendersExitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** informs Benders' decomposition that the branch and bound process is being started */
SCIP_RETCODE SCIPbendersInitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs Benders' decomposition that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbendersExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** activates Benders' decomposition such that it is called in LP solving loop */
SCIP_RETCODE SCIPbendersActivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nsubproblems        /**< the number subproblems used in this decomposition */
   );

/** deactivates Benders' decomposition such that it is no longer called in LP solving loop */
SCIP_RETCODE SCIPbendersDeactivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** enables or disables all clocks of Benders' decomposition depending on the value of the flag */
void SCIPbendersEnableOrDisableClocks(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the Benders' decomposition be enabled? */
   );

/** solves the subproblem using the current master problem solution.
 *
 *  The checkint flag indicates whether integer feasibility can be assumed. If it is not assumed, i.e. checkint ==
 *  FALSE, then only the convex relaxations of the subproblems are solved. If integer feasibility is assumed, i.e.
 *  checkint == TRUE, then the convex relaxations and the full CIP are solved to generate Benders' cuts and check
 *  solution feasibility.
 */
SCIP_RETCODE SCIPbendersExec(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_Bool*            auxviol,            /**< set to TRUE only if the solution is feasible but the aux vars are violated */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should the integer solution be checked by the subproblems */
   );

/** Executes the subproblem solving process. */
SCIP_RETCODE SCIPbendersExecSubproblemSolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnum,            /**< the subproblem number */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the solve loop iteration. The first iter is for LP, the second for IP */
   SCIP_Bool             enhancement,        /**< is the solve performed as part of an enhancement? */
   SCIP_Bool*            solved,             /**< flag to indicate whether the subproblem was solved */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   );

/** sets up the subproblem using the solution to the master problem  */
SCIP_RETCODE SCIPbendersSetupSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   );

/** Solve a Benders' decomposition subproblems. This will either call the user defined method or the generic solving
 *  methods. If the generic method is called, then the subproblem must be set up before calling this method. */
SCIP_RETCODE SCIPbendersSolveSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_Bool             solvecip,           /**< directly solve the CIP subproblem */
   SCIP_Real*            objective           /**< the objective function value of the subproblem, can be NULL */
   );

/** frees the subproblems */
SCIP_RETCODE SCIPbendersFreeSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   );

/** compares the subproblem objective value with the auxiliary variable value for optimality */
SCIP_Bool SCIPbendersSubproblemIsOptimal(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber          /**< the subproblem number */
   );

/** returns the value of the auxiliary variable value in a master problem solution */
SCIP_Real SCIPbendersGetAuxiliaryVarVal(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber          /**< the subproblem number */
   );

/** Solves an independent subproblem to identify its lower bound. The lower bound is then used to update the bound on
 *  the auxiliary variable.
 */
SCIP_RETCODE SCIPbendersComputeSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber,         /**< the subproblem to be evaluated */
   SCIP_Real*            lowerbound,         /**< the lowerbound for the subproblem */
   SCIP_Bool*            infeasible          /**< was the subproblem found to be infeasible? */
   );

/** merges a subproblem into the master problem. This process just adds a copy of the subproblem variables and
 *  constraints to the master problem, but keeps the subproblem stored in the Benders' decomposition data structure.
 *  The reason for keeping the subproblem available is for when it is queried for solutions after the problem is solved.
 *
 *  Once the subproblem is merged into the master problem, then the subproblem is flagged as disabled. This means that
 *  it will not be solved in the subsequent subproblem solving loops.
 *
 *  The associated auxiliary variables are kept in the master problem. The objective function of the merged subproblem
 *  is added as an underestimator constraint.
 */
SCIP_RETCODE SCIPbendersMergeSubproblemIntoMaster(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of subproblem variables corresponding
                                              *   to the newly created master variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of subproblem constraints to the
                                                  corresponding newly created constraints, or NULL */
   int                   probnumber          /**< the number of the subproblem that will be merged into the master problem*/
   );

/** Applies a Benders' decomposition to the problem based upon the decomposition selected from the storage */
extern
SCIP_RETCODE SCIPbendersApplyDecomposition(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DECOMP*          decomp              /**< the decomposition to apply to the problem */
   );

/** sets priority of Benders' decomposition */
void SCIPbendersSetPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the Benders' decomposition */
   );

/** sets copy callback of Benders' decomposition */
void SCIPbendersSetCopy(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSCOPY ((*benderscopy))    /**< copy callback of Benders' decomposition */
   );

/** sets destructor callback of Benders' decomposition */
void SCIPbendersSetFree(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREE ((*bendersfree))    /**< destructor of Benders' decomposition */
   );

/** sets initialization callback of Benders' decomposition */
void SCIPbendersSetInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINIT((*bendersinit))     /**< initialize Benders' decomposition */
   );

/** sets deinitialization callback of Benders' decomposition */
void SCIPbendersSetExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXIT((*bendersexit))     /**< deinitialize Benders' decomposition */
   );

/** sets presolving initialization callback of Benders' decomposition */
void SCIPbendersSetInitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre))/**< initialize presolving for Benders' decomposition */
   );

/** sets presolving deinitialization callback of Benders' decomposition */
void SCIPbendersSetExitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre))/**< deinitialize presolving for Benders' decomposition */
   );

/** sets solving process initialization callback of Benders' decomposition */
void SCIPbendersSetInitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol))/**< solving process initialization callback of Benders' decomposition */
   );

/** sets solving process deinitialization callback of Benders' decomposition */
void SCIPbendersSetExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol))/**< solving process deinitialization callback of Benders' decomposition */
   );

/** sets the pre subproblem solve callback of Benders' decomposition */
void SCIPbendersSetPresubsolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve))/**< called prior to the subproblem solving loop */
   );

/** sets convex solve callback of Benders' decomposition */
void SCIPbendersSetSolvesubconvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex))/**< solving method for the convex Benders' decomposition subproblem */
   );

/** sets solve callback of Benders' decomposition */
void SCIPbendersSetSolvesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub))/**< solving method for a Benders' decomposition subproblem */
   );

/** sets post-solve callback of Benders' decomposition */
void SCIPbendersSetPostsolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve))/**< solving process deinitialization callback of Benders' decomposition */
   );

/** sets post-solve callback of Benders' decomposition */
void SCIPbendersSetSubproblemComp(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_SORTPTRCOMP((*benderssubcomp))  /**< a comparator for defining the solving order of the subproblems */
   );

/** sets free subproblem callback of Benders' decomposition */
void SCIPbendersSetFreesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub))/**< the freeing callback for the subproblem */
   );

/** Returns the corresponding master or subproblem variable for the given variable.
 *  This provides a call back for the variable mapping between the master and subproblems. */
SCIP_RETCODE SCIPbendersGetVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< the variable for which the corresponding variable is desired */
   SCIP_VAR**            mappedvar,          /**< the variable that is mapped to var */
   int                   probnumber          /**< the problem number for the desired variable, -1 for the master problem */
   );

/** adds a subproblem to the Benders' decomposition data */
SCIP_RETCODE SCIPbendersAddSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP*                 subproblem          /**< subproblem to be added to the data storage */
   );

/** removes the subproblems from the Benders' decomposition data */
void SCIPbendersRemoveSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** Sets whether the subproblem is enabled or disabled. A subproblem is disabled if it has been merged into the master
 *  problem.
 */
void SCIPbendersSetSubproblemEnabled(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             enabled             /**< flag to indicate whether the subproblem is enabled */
   );

/** changes all of the master problem variables in the given subproblem to continuous */
SCIP_RETCODE SCIPbendersChgMastervarsToCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   );

/** sets a flag to indicate whether the master variables are all set to continuous */
SCIP_RETCODE SCIPbendersSetMastervarsCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             arecont             /**< flag to indicate whether the master problem variables are continuous */
   );

/** returns whether the master variables are all set to continuous */
SCIP_Bool SCIPbendersGetMastervarsCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   );

/** adds the data for the generated cuts to the Benders' cut storage */
SCIP_RETCODE SCIPbendersStoreCut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition cut */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real*            vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real             lhs,                /**< the left hand side of the cut */
   SCIP_Real             rhs,                /**< the right hand side of the cut */
   int                   nvars               /**< the number of variables with non-zero coefficients in the cut */
   );

/** inserts a Benders' cut algorithm plugin into the Benders' cuts plugin list */
SCIP_RETCODE SCIPbendersIncludeBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERSCUT*      benderscut          /**< Benders' cut */
   );

/** sets the Benders' cuts sorted flags in the Benders' decomposition */
void SCIPbendersSetBenderscutsSorted(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_Bool             sorted              /**< the value to set the sorted flag to */
   );

/** sorts Benders' decomposition cuts by priorities */
void SCIPbendersSortBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** sorts Benders' decomposition cuts by name */
void SCIPbendersSortBenderscutsName(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

#ifdef __cplusplus
}
#endif

#endif
