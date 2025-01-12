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

/**@file   conflict.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONFLICT_H__
#define __SCIP_CONFLICT_H__


#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_conflictstore.h"
#include "scip/type_event.h"
#include "scip/type_implics.h"
#include "scip/type_lp.h"
#include "scip/type_message.h"
#include "scip/type_prob.h"
#include "scip/type_reopt.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Conflict Handler
 */

/** copies the given conflict handler to a new scip */
SCIP_RETCODE SCIPconflicthdlrCopyInclude(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a conflict handler */
SCIP_RETCODE SCIPconflicthdlrCreate(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of conflict handler */
   const char*           desc,               /**< description of conflict handler */
   int                   priority,           /**< priority of the conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy)),  /**< copy method of conflict handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol)),/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol)),/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   );

/** calls destructor and frees memory of conflict handler */
SCIP_RETCODE SCIPconflicthdlrFree(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls init method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrInit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs conflict handler that the branch and bound process is being started */
SCIP_RETCODE SCIPconflicthdlrInitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs conflict handler that the branch and bound process data is being freed */
SCIP_RETCODE SCIPconflicthdlrExitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExec(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node to add conflict constraint to */
   SCIP_NODE*            validnode,          /**< node at which the constraint is valid */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change resembling the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict */
   int                   nbdchginfos,        /**< number of bound changes in the conflict set */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             usescutoffbound,    /**< depends the conflict on the cutoff bound? */
   SCIP_Bool             resolved,           /**< was the conflict set already used to create a constraint? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of conflict handler */
void SCIPconflicthdlrSetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the conflict handler */
   );

/** set copy method of conflict handler */
void SCIPconflicthdlrSetCopy(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy))   /**< copy method of the conflict handler */
   );

/** set destructor of conflict handler */
void SCIPconflicthdlrSetFree(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTFREE((*conflictfree))   /**< destructor of conflict handler */
   );

/** set initialization method of conflict handler */
void SCIPconflicthdlrSetInit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit))   /**< initialization method conflict handler */
   );

/** set deinitialization method of conflict handler */
void SCIPconflicthdlrSetExit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit))   /**< deinitialization method conflict handler */
   );

/** set solving process initialization method of conflict handler */
void SCIPconflicthdlrSetInitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol))/**< solving process initialization method of conflict handler */
   );

/** set solving process deinitialization method of conflict handler */
void SCIPconflicthdlrSetExitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol))/**< solving process deinitialization method of conflict handler */
   );

/** enables or disables all clocks of \p conflicthdlr, depending on the value of the flag */
void SCIPconflicthdlrEnableOrDisableClocks(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< the conflict handler for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the conflict handler be enabled? */
   );

/*
 * Conflict Analysis
 */

/** return TRUE if conflict analysis is applicable; In case the function return FALSE there is no need to initialize the
 *  conflict analysis since it will not be applied
 */
SCIP_Bool SCIPconflictApplicable(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** creates conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictCreate(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** frees conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictFree(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   );

/** initializes the propagation conflict analysis by clearing the conflict candidate queue */
SCIP_RETCODE SCIPconflictInit(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             usescutoffbound     /**< depends the conflict on a cutoff bound? */
   );

/** adds variable's bound to conflict candidate queue */
SCIP_RETCODE SCIPconflictAddBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   );

/** adds variable's bound to conflict candidate queue with the additional information of a relaxed bound */
SCIP_RETCODE SCIPconflictAddRelaxedBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Real             relaxedbd           /**< the relaxed bound */
   );

/** checks if the given variable is already part of the current conflict set or queued for resolving with the same or
 *  even stronger bound
 */
SCIP_RETCODE SCIPconflictIsVarUsed(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for which the score should be increased */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Bool*            used                /**< pointer to store if the variable is already used */
   );

/** returns the conflict lower bound if the variable is present in the current conflict set; otherwise the global lower
 *  bound
 */
SCIP_Real SCIPconflictGetVarLb(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the conflict upper bound if the variable is present in the current conflict set; otherwise the global upper
 *  bound
 */
SCIP_Real SCIPconflictGetVarUb(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound() and
 *  SCIPconflictAddRelaxedBound(), and on success, calls the conflict handlers to create a conflict constraint out of
 *  the resulting conflict set; updates statistics for propagation conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyze(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/** adds the collected conflict constraints to the corresponding nodes; the best set->conf_maxconss conflict constraints
 *  are added to the node of their validdepth; additionally (if not yet added, and if repropagation is activated), the
 *  conflict constraint that triggers the earliest repropagation is added to the node of its validdepth
 */
SCIP_RETCODE SCIPconflictFlushConss(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** returns the current number of conflict sets in the conflict set storage */
int SCIPconflictGetNConflicts(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** returns the total number of conflict constraints that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** returns the total number of literals in conflict constraints that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** returns the total number of global bound changes applied by the conflict analysis */
SCIP_Longint SCIPconflictGetNGlobalChgBds(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** returns the total number of conflict constraints that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** returns the total number of literals in conflict constraints that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** returns the total number of local bound changes applied by the conflict analysis */
SCIP_Longint SCIPconflictGetNLocalChgBds(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** returns the total number of conflict constraints that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** returns the total number of literals in conflict constraints that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets time in seconds used for preprocessing global conflict constraint before appliance */
SCIP_Real SCIPconflictGetGlobalApplTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets time in seconds used for analyzing propagation conflicts */
SCIP_Real SCIPconflictGetPropTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to propagation conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNPropSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict constraints detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict constraints created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence constraints detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence constraints created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );




/*
 * Infeasible LP Conflict Analysis
 */

/** analyzes an infeasible or bound exceeding LP to find out the bound changes on variables that were responsible for the
 *  infeasibility or for exceeding the primal bound;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible or bound exceeding LP conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/** gets time in seconds used for analyzing infeasible LP conflicts */
SCIP_Real SCIPconflictGetInfeasibleLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to infeasible LP conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNInfeasibleLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict constraints detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict constraints created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence constraints detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence constraints created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of LP iterations in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets time in seconds used for analyzing bound exceeding LP conflicts */
SCIP_Real SCIPconflictGetBoundexceedingLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to bound exceeding LP conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNBoundexceedingLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict constraints detected in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict constraints created in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence constraints detected in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence constraints created in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of LP iterations in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );




/*
 * infeasible strong branching conflict analysis
 */

/** analyses infeasible strong branching sub problems for conflicts */
SCIP_RETCODE SCIPconflictAnalyzeStrongbranch(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_COL*             col,                /**< LP column with at least one infeasible strong branching subproblem */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict          /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   );

/** gets time in seconds used for analyzing infeasible strong branching conflicts */
SCIP_Real SCIPconflictGetStrongbranchTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of successful calls to dual proof analysis derived from infeasible LPs */
SCIP_Longint SCIPconflictGetNDualproofsInfSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of globally valid dual proof constraints derived from infeasible LPs */
SCIP_Longint SCIPconflictGetNDualproofsInfGlobal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of locally valid dual proof constraints derived from infeasible LPs */
SCIP_Longint SCIPconflictGetNDualproofsInfLocal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets average length of dual proof constraints derived from infeasible LPs */
SCIP_Longint SCIPconflictGetNDualproofsInfNonzeros(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of successfully analyzed dual proofs derived from bound exceeding LPs */
SCIP_Longint SCIPconflictGetNDualproofsBndSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of globally applied dual proofs derived from bound exceeding LPs */
SCIP_Longint SCIPconflictGetNDualproofsBndGlobal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of locally applied dual proofs derived from bound exceeding LPs */
SCIP_Longint SCIPconflictGetNDualproofsBndLocal(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets average length of dual proofs derived from bound exceeding LPs */
SCIP_Longint SCIPconflictGetNDualproofsBndNonzeros(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to infeasible strong branching conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNStrongbranchSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict constraints detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict constraints created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence constraints detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence constraints created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of LP iterations in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );




/*
 * pseudo solution conflict analysis
 */

/** analyzes a pseudo solution with objective value exceeding the current cutoff to find out the bound changes on
 *  variables that were responsible for the objective value degradation;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for pseudo solution conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyzePseudo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   );

/** gets time in seconds used for analyzing pseudo solution conflicts */
SCIP_Real SCIPconflictGetPseudoTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to pseudo solution conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNPseudoSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of conflict constraints detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in conflict constraints created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of reconvergence constraints detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets total number of literals in reconvergence constraints created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** enables or disables all clocks of \p conflict, depending on the value of the flag */
void SCIPconflictEnableOrDisableClocks(
   SCIP_CONFLICT*        conflict,           /**< the conflict analysis data for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the conflict analysis data be enabled? */
   );

#ifdef __cplusplus
}
#endif

#endif
