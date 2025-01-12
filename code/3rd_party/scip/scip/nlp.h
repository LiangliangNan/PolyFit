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

/**@file   nlp.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for NLP management
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLP_H__
#define __SCIP_NLP_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_event.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_misc.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"
#include "scip/pub_nlp.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@name Nonlinear row methods */
/**@{ */

/** create a new nonlinear row
 *
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreate(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             constant,           /**< constant */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< linear variables, or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< linear coefficients, or NULL if nlinvars == 0 */
   SCIP_EXPR*            expr,               /**< expression, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_EXPRCURV         curvature           /**< curvature of the nonlinear row */
   );

/** create a nonlinear row that is a copy of a given row
 *
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreateCopy(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           sourcenlrow         /**< nonlinear row to copy */
   );

/** create a new nonlinear row from a linear row
 *
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreateFromRow(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_ROW*             row                 /**< the linear row to copy */
   );

/** output nonlinear row to file stream */
SCIP_RETCODE SCIPnlrowPrint(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** increases usage counter of nonlinear row */
void SCIPnlrowCapture(
   SCIP_NLROW*           nlrow               /**< nonlinear row to capture */
   );

/** decreases usage counter of nonlinear row */
SCIP_RETCODE SCIPnlrowRelease(
   SCIP_NLROW**          nlrow,              /**< nonlinear row to free */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** ensures, that linear coefficient array of nonlinear row can store at least num entries */
SCIP_RETCODE SCIPnlrowEnsureLinearSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** adds a previously non existing linear coefficient to a nonlinear row */
SCIP_RETCODE SCIPnlrowAddLinearCoef(
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             val                 /**< value of coefficient */
   );

/** deletes linear coefficient from nonlinear row */
SCIP_RETCODE SCIPnlrowDelLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var                 /**< coefficient to be deleted */
   );

/** changes or adds a linear coefficient to a nonlinear row */
SCIP_RETCODE SCIPnlrowChgLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< new value of coefficient */
   );

/** replaces or deletes an expression in a nonlinear row */
SCIP_RETCODE SCIPnlrowChgExpr(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_EXPR*            expr                /**< new expression, or NULL to delete current one */
   );

/** changes constant of nonlinear row */
SCIP_RETCODE SCIPnlrowChgConstant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             constant            /**< new constant */
   );

/** changes left hand side of nonlinear row */
SCIP_RETCODE SCIPnlrowChgLhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of nonlinear row */
SCIP_RETCODE SCIPnlrowChgRhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             rhs                 /**< new right hand side */
   );

/** removes (or substitutes) all fixed, negated, aggregated, multi-aggregated variables from the linear and nonlinear part of a nonlinear row and simplifies its expression */
SCIP_RETCODE SCIPnlrowSimplify(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** recalculates the current activity of a nonlinear row in the current NLP solution */
SCIP_RETCODE SCIPnlrowRecalcNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives the activity of a nonlinear row in the current NLP solution */
SCIP_RETCODE SCIPnlrowGetNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            activity            /**< buffer to store activity value */
   );

/** gives the feasibility of a nonlinear row in the current NLP solution: negative value means infeasibility */
SCIP_RETCODE SCIPnlrowGetNLPFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   );

/** calculates the current pseudo activity of a nonlinear row */
SCIP_RETCODE SCIPnlrowRecalcPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< SCIP problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< SCIP LP */
   );

/** returns the pseudo activity of a nonlinear row in the current pseudo solution */
SCIP_RETCODE SCIPnlrowGetPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< SCIP problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< SCIP LP */
   SCIP_Real*            pseudoactivity      /**< buffer to store pseudo activity value */
   );

/** returns the pseudo feasibility of a nonlinear row in the current pseudo solution: negative value means infeasibility */
SCIP_RETCODE SCIPnlrowGetPseudoFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< SCIP problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< SCIP LP */
   SCIP_Real*            pseudofeasibility   /**< buffer to store pseudo feasibility value */
   );

/** returns the activity of a nonlinear row for a given solution */
SCIP_RETCODE SCIPnlrowGetSolActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            activity            /**< buffer to store activity value */
   );

/** returns the feasibility of a nonlinear row for the given solution */
SCIP_RETCODE SCIPnlrowGetSolFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   );

/** returns the minimal activity of a nonlinear row w.r.t. the variables' bounds */
SCIP_RETCODE SCIPnlrowGetActivityBounds(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Real*            minactivity,        /**< buffer to store minimal activity, or NULL */
   SCIP_Real*            maxactivity         /**< buffer to store maximal activity, or NULL */
   );

/** returns whether the nonlinear row is redundant w.r.t. the variables' bounds */
SCIP_RETCODE SCIPnlrowIsRedundant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Bool*            isredundant         /**< buffer to store whether row is redundant */
   );

/**@} */

/**@name NLP methods */
/**@{ */

/** includes event handler that is used by NLP */
SCIP_RETCODE SCIPnlpInclude(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** construct a new empty NLP */
SCIP_RETCODE SCIPnlpCreate(
   SCIP_NLP**            nlp,                /**< NLP handler, call by reference */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< problem name */
   int                   nvars_estimate      /**< an estimate on the number of variables that may be added to the NLP later */
   );

/** frees NLP data object */
SCIP_RETCODE SCIPnlpFree(
   SCIP_NLP**            nlp,                /**< pointer to NLP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   );

/** resets the NLP to the empty NLP by removing all variables and rows from NLP,
 *  releasing all rows, and flushing the changes to the NLP solver
 */
SCIP_RETCODE SCIPnlpReset(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   );

/** currently a dummy function that always returns TRUE */
SCIP_Bool SCIPnlpHasCurrentNodeNLP(
   SCIP_NLP*             nlp                 /**< NLP data */
   );

/** ensures, that variables array of NLP can store at least num entries */
SCIP_RETCODE SCIPnlpEnsureVarsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** adds a variable to the NLP and captures the variable */
SCIP_RETCODE SCIPnlpAddVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable */
   );

/** adds a set of variables to the NLP and captures the variables */
SCIP_RETCODE SCIPnlpAddVars(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables to add */
   SCIP_VAR**            vars                /**< variables to add */
   );

/** deletes a variable from the NLP and releases the variable */
SCIP_RETCODE SCIPnlpDelVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< SCIP LP, needed to release variable */
   SCIP_VAR*             var                 /**< variable */
   );

/** ensures, that nonlinear rows array of NLP can store at least num entries */
SCIP_RETCODE SCIPnlpEnsureNlRowsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** adds a nonlinear row to the NLP and captures it
 *
 * all variables of the row need to be present in the NLP
 */
SCIP_RETCODE SCIPnlpAddNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   );

/** adds nonlinear rows to the NLP and captures them
 *
 * all variables of the row need to be present in the NLP
 */
SCIP_RETCODE SCIPnlpAddNlRows(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   nnlrows,            /**< number of rows to add */
   SCIP_NLROW**          nlrows              /**< rows to add */
   );

/** deletes a nonlinear row from the NLP
 *
 * does nothing if nonlinear row is not in NLP
 */
SCIP_RETCODE SCIPnlpDelNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   );

/** applies all cached changes to the NLP solver */
SCIP_RETCODE SCIPnlpFlush(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** solves the NLP or diving NLP */
SCIP_RETCODE SCIPnlpSolve(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLPPARAM*        nlpparam            /**< NLP solve parameters */
   );

/** gets objective value of current NLP */
SCIP_Real SCIPnlpGetObjval(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives current pseudo objective value */
SCIP_RETCODE SCIPnlpGetPseudoObjval(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< SCIP problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< SCIP LP */
   SCIP_Real*            pseudoobjval        /**< buffer to store pseudo objective value */
   );

/** gets fractional variables of last NLP solution along with solution values and fractionalities
 */
SCIP_RETCODE SCIPnlpGetFracVars(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR***           fracvars,           /**< pointer to store the array of NLP fractional variables, or NULL */
   SCIP_Real**           fracvarssol,        /**< pointer to store the array of NLP fractional variables solution values, or NULL */
   SCIP_Real**           fracvarsfrac,       /**< pointer to store the array of NLP fractional variables fractionalities, or NULL */
   int*                  nfracvars,          /**< pointer to store the number of NLP fractional variables , or NULL */
   int*                  npriofracvars       /**< pointer to store the number of NLP fractional variables with maximal branching priority, or NULL */
   );

/** removes all redundant nonlinear rows */
SCIP_RETCODE SCIPnlpRemoveRedundantNlRows(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** set initial guess (approximate primal solution) for next solve
 *
 *  array initguess must be NULL or have length at least SCIPnlpGetNVars()
 */
SCIP_RETCODE SCIPnlpSetInitialGuess(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_Real*            initguess           /**< new initial guess, or NULL to clear previous one */
   );

/** writes NLP to a file */
SCIP_RETCODE SCIPnlpWrite(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           fname               /**< file name */
   );

/*
 * NLP diving methods
 */

/** signals start of diving */
SCIP_RETCODE SCIPnlpStartDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** resets the bound and objective changes made during diving and disables diving mode */
SCIP_RETCODE SCIPnlpEndDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** changes coefficient of variable in diving NLP */
SCIP_RETCODE SCIPnlpChgVarObjDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             coef                /**< new linear coefficient of variable in objective */
   );

/** changes bounds of variable in diving NLP */
SCIP_RETCODE SCIPnlpChgVarBoundsDive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable which bounds to change */
   SCIP_Real             lb,                 /**< new lower bound of variable */
   SCIP_Real             ub                  /**< new upper bound of variable */
   );

/** changes bounds of a set of variables in diving NLP */
SCIP_RETCODE SCIPnlpChgVarsBoundsDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables which bounds to change */
   SCIP_VAR**            vars,               /**< variables which bounds to change */
   SCIP_Real*            lbs,                /**< new lower bounds of variables */
   SCIP_Real*            ubs                 /**< new upper bounds of variables */
   );

/** returns whether the objective function has been changed during diving */
SCIP_Bool SCIPnlpIsDivingObjChanged(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets array with variables of the NLP */
SCIP_VAR** SCIPnlpGetVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets current number of variables in NLP */
int SCIPnlpGetNVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** computes for each variables the number of NLP rows in which the variable appears in a nonlinear var */
SCIP_RETCODE SCIPnlpGetVarsNonlinearity(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   int*                  nlcount             /**< an array of length at least SCIPnlpGetNVars() to store nonlinearity counts of variables */
   );

/** indicates whether there exists a row that contains a continuous variable in a nonlinear term
 *
 * @note The method may have to touch every row and nonlinear term to compute its result.
 */
SCIP_RETCODE SCIPnlpHasContinuousNonlinearity(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            result              /**< buffer to store whether continuous variable present in an expression of any row */
   );

/** gives dual solution values associated with lower bounds of NLP variables */
SCIP_Real* SCIPnlpGetVarsLbDualsol(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives dual solution values associated with upper bounds of NLP variables */
SCIP_Real* SCIPnlpGetVarsUbDualsol(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets array with nonlinear rows of the NLP */
SCIP_NLROW** SCIPnlpGetNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets current number of nonlinear rows in NLP */
int SCIPnlpGetNNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets the NLP solver interface */
SCIP_NLPI* SCIPnlpGetNLPI(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets the NLP problem in the solver interface */
SCIP_NLPIPROBLEM* SCIPnlpGetNLPIProblem(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** indicates whether NLP is currently in diving mode */
SCIP_Bool SCIPnlpIsDiving(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets solution status of current NLP */
SCIP_NLPSOLSTAT SCIPnlpGetSolstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets termination status of last NLP solve */
SCIP_NLPTERMSTAT SCIPnlpGetTermstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives statistics (number of iterations, solving time, ...) of last NLP solve */
SCIP_RETCODE SCIPnlpGetStatistics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   );

/** indicates whether a solution for the current NLP is available
 *
 * The solution may be optimal, feasible, or infeasible.
 * Thus, returns whether the NLP solution status is at most \ref SCIP_NLPSOLSTAT_LOCINFEASIBLE.
 */
SCIP_Bool SCIPnlpHasSolution(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLP_H__ */
