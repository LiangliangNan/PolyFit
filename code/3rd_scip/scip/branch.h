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

/**@file   branch.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for branching rules and branching candidate storage
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_H__
#define __SCIP_BRANCH_H__


#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/type_branch.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_message.h"
#include "scip/type_prob.h"
#include "scip/type_reopt.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sepastore.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * branching candidate storage methods
 */

/** creates a branching candidate storage */
SCIP_RETCODE SCIPbranchcandCreate(
   SCIP_BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   );

/** frees branching candidate storage */
SCIP_RETCODE SCIPbranchcandFree(
   SCIP_BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   );

/** invalidates branching candidates storage */
void SCIPbranchcandInvalidate(
   SCIP_BRANCHCAND*      branchcand          /**< pointer to store branching candidate storage */
   );

/** gets branching candidates for LP solution branching (fractional variables) */
SCIP_RETCODE SCIPbranchcandGetLPCands(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   SCIP_Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   SCIP_Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*                  nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*                  npriolpcands,       /**< pointer to store the number of candidates with maximal priority, or NULL */
   int*                  nfracimplvars       /**< pointer to store the number of implicit fractional variables, or NULL */
   );


/** gets external branching candidates */
SCIP_RETCODE SCIPbranchcandGetExternCands(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR***           externcands,        /**< pointer to store the array of external branching candidates, or NULL */
   SCIP_Real**           externcandssol,     /**< pointer to store the array of external candidate solution values, or NULL */
   SCIP_Real**           externcandsscore,   /**< pointer to store the array of external candidate scores, or NULL */
   int*                  nexterncands,       /**< pointer to store the number of external branching candidates, or NULL */
   int*                  nprioexterncands,   /**< pointer to store the number of candidates with maximal priority, or NULL */
   int*                  nprioexternbins,    /**< pointer to store the number of binary candidates with maximal priority, or NULL */
   int*                  nprioexternints,    /**< pointer to store the number of integer candidates with maximal priority, or NULL */
   int*                  nprioexternimpls    /**< pointer to store the number of implicit integer candidates with maximal priority, 
                                              *   or NULL */
   );

/** gets maximal branching priority of LP branching candidates */
int SCIPbranchcandGetLPMaxPrio(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of LP branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioLPCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets maximal branching priority of external branching candidates */
int SCIPbranchcandGetExternMaxPrio(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of external branching candidates */
int SCIPbranchcandGetNExternCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of binary external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternBins(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of integer external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternInts(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of implicit integer external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternImpls(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of continuous external branching candidates with maximal branch priority */
int SCIPbranchcandGetNPrioExternConts(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** insert variable, its score and its solution value into the external branching candidate storage
 * the absolute difference of the current lower and upper bounds of the variable must be at least epsilon
 */
SCIP_RETCODE SCIPbranchcandAddExternCand(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to insert */
   SCIP_Real             score,              /**< score of external candidate, e.g. infeasibility */
   SCIP_Real             solval              /**< value of the variable in the current solution */
   );

/** removes all external candidates from the storage for external branching */
void SCIPbranchcandClearExternCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** checks whether the given variable is contained in the candidate storage for external branching */
SCIP_Bool SCIPbranchcandContainsExternCand(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR*             var                 /**< variable to look for */
   );

/** gets branching candidates for pseudo solution branching (non-fixed variables) */
SCIP_RETCODE SCIPbranchcandGetPseudoCands(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*                  npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*                  npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   );

/** gets number of branching candidates for pseudo solution branching (non-fixed variables) */
int SCIPbranchcandGetNPseudoCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoCands(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoBins(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoInts(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoImpls(
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   );

/** removes variable from branching candidate list */
SCIP_RETCODE SCIPbranchcandRemoveVar(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR*             var                 /**< variable that changed its bounds */
   );

/** updates branching candidate list for a given variable */
SCIP_RETCODE SCIPbranchcandUpdateVar(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that changed its bounds */
   );

/** updates branching priority of the given variable and update the pseudo candidate array if needed */
SCIP_RETCODE SCIPbranchcandUpdateVarBranchPriority(
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable that changed its bounds */
   int                   branchpriority      /**< branch priority of the variable */
   );




/*
 * branching rules
 */

/** copies the given branchrule to a new scip */
SCIP_RETCODE SCIPbranchruleCopyInclude(
   SCIP_BRANCHRULE*      branchrule,         /**< branchrule */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a branching rule */
SCIP_RETCODE SCIPbranchruleCreate(
   SCIP_BRANCHRULE**     branchrule,         /**< pointer to store branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of branching rule */
   const char*           desc,               /**< description of branching rule */
   int                   priority,           /**< priority of the branching rule */
   int                   maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_DECL_BRANCHCOPY  ((*branchcopy)),    /**< copy method of branching rule */
   SCIP_DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol)),/**< solving process initialization method of branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol)),/**< solving process deinitialization method of branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   SCIP_DECL_BRANCHEXECEXT((*branchexecext)),/**< branching execution method for external solutions */
   SCIP_DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   );

/** frees memory of branching rule */
SCIP_RETCODE SCIPbranchruleFree(
   SCIP_BRANCHRULE**     branchrule,         /**< pointer to branching rule data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes branching rule */
SCIP_RETCODE SCIPbranchruleInit(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** deinitializes branching rule */
SCIP_RETCODE SCIPbranchruleExit(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs branching rule that the branch and bound process is being started */
SCIP_RETCODE SCIPbranchruleInitsol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs branching rule that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbranchruleExitsol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** executes branching rule for fractional LP solution */
SCIP_RETCODE SCIPbranchruleExecLPSol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** executes branching rule for external branching candidates */
SCIP_RETCODE SCIPbranchruleExecExternSol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** executes branching rule for not completely fixed pseudo solution */
SCIP_RETCODE SCIPbranchruleExecPseudoSol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of branching rule */
void SCIPbranchruleSetPriority(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the branching rule */
   );

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
void SCIPbranchruleSetMaxdepth(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   int                   maxdepth            /**< new maxdepth of the branching rule */
   );

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
void SCIPbranchruleSetMaxbounddist(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_Real             maxbounddist        /**< new maxbounddist of the branching rule */
   );

/** sets copy method of branching rule */
void SCIPbranchruleSetCopy(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHCOPY  ((*branchcopy))     /**< copy method of branching rule or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of branching rule */
void SCIPbranchruleSetFree(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHFREE  ((*branchfree))     /**< destructor of branching rule */
   );

/** sets initialization method of branching rule */
void SCIPbranchruleSetInit(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit))     /**< initialize branching rule */
   );

/** sets deinitialization method of branching rule */
void SCIPbranchruleSetExit(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit))     /**< deinitialize branching rule */
   );

/** sets solving process initialization method of branching rule */
void SCIPbranchruleSetInitsol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol)) /**< solving process initialization method of branching rule */
   );

/** sets solving process deinitialization method of branching rule */
void SCIPbranchruleSetExitsol(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol)) /**< solving process deinitialization method of branching rule */
   );

/** sets branching execution method for fractional LP solutions */
void SCIPbranchruleSetExecLp(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp))   /**< branching execution method for fractional LP solutions */
   );

/** sets branching execution method for external candidates  */
void SCIPbranchruleSetExecExt(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECEXT((*branchexecext)) /**< branching execution method for external candidates */
   );

/** sets branching execution method for not completely fixed pseudo solutions */
void SCIPbranchruleSetExecPs(
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_DECL_BRANCHEXECPS((*branchexecps))   /**< branching execution method for not completely fixed pseudo solutions */
   );

/** enables or disables all clocks of \p branchrule, depending on the value of the flag */
void SCIPbranchruleEnableOrDisableClocks(
   SCIP_BRANCHRULE*      branchrule,         /**< the branching rule for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the branching rule be enabled? */
   );

/*
 * branching methods
 */

/** calculates the branching score out of the gain predictions for a binary branching */
SCIP_Real SCIPbranchGetScore(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_Real             downgain,           /**< prediction of objective gain for rounding downwards */
   SCIP_Real             upgain              /**< prediction of objective gain for rounding upwards */
   );

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children */
SCIP_Real SCIPbranchGetScoreMultiple(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int                   nchildren,          /**< number of children that the branching will create */
   SCIP_Real*            gains               /**< prediction of objective gain for each child */
   );

/** computes a branching point for a (not necessarily discrete) variable
 * a suggested branching point is first projected onto the box
 * if no point is suggested, then the value in the current LP or pseudo solution is used
 * if this value is at infinity, then 0.0 projected onto the bounds and then moved inside the interval is used 
 * for a discrete variable, it is ensured that the returned value is fractional
 * for a continuous variable, the parameter branching/clamp defines how far a branching point need to be from the bounds of a variable
 * the latter is only applied if no point has been suggested, or the suggested point is not inside the variable's interval
 */
SCIP_Real SCIPbranchGetBranchingPoint(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable, of which the branching point should be computed */
   SCIP_Real             suggestion          /**< suggestion for branching point, or SCIP_INVALID if no suggestion */
   );

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
SCIP_RETCODE SCIPbranchExecLP(
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching */
   );

/** calls branching rules to branch on an external solution; if no external branching candidates exist, the result is SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPbranchExecExtern(
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching */
   );

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPbranchExecPseudo(
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             cutoffbound,        /**< global upper cutoff bound */
   SCIP_Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching */
   );

#ifdef __cplusplus
}
#endif

#endif
