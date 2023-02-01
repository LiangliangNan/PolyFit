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

/**@name Expressions and Expression tree methods */
/**@{ */

/** removes fixed variables from an expression tree, so that at exit all variables are active */
extern
SCIP_RETCODE SCIPexprtreeRemoveFixedVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool*            changed,            /**< buffer to store whether the tree was changed, i.e., whether there was a fixed variable */
   int*                  varpos,             /**< array of length at least tree->nvars to store new indices of previously existing variables in expression tree, or -1 if variable was removed; set to NULL if not of interest */
   int*                  newvarsstart        /**< buffer to store index in tree->vars array where new variables begin, or NULL if not of interest */
   );

/**@} */

/**@name Nonlinear row methods */
/**@{ */

/** create a new nonlinear row
 * the new row is already captured
 */
extern
SCIP_RETCODE SCIPnlrowCreate(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             constant,           /**< constant */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< linear variables, or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< linear coefficients, or NULL if nlinvars == 0 */
   int                   nquadvars,          /**< number variables in quadratic terms */
   SCIP_VAR**            quadvars,           /**< variables in quadratic terms, or NULL if nquadvars == 0 */
   int                   nquadelems,         /**< number of entries in quadratic term matrix */
   SCIP_QUADELEM*        quadelems,          /**< elements of quadratic term matrix, or NULL if nquadelems == 0 */
   SCIP_EXPRTREE*        exprtree,           /**< expression tree, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_EXPRCURV         curvature           /**< curvature of the nonlinear row */
   );

/** create a nonlinear row that is a copy of a given row
 * the new row is already captured
 */
extern
SCIP_RETCODE SCIPnlrowCreateCopy(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           sourcenlrow         /**< nonlinear row to copy */
   );

/** create a new nonlinear row from a linear row
 * the new row is already captured
 */
extern
SCIP_RETCODE SCIPnlrowCreateFromRow(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< the linear row to copy */
   );

/** frees a nonlinear row */
extern
SCIP_RETCODE SCIPnlrowFree(
   SCIP_NLROW**          nlrow,              /**< pointer to NLP row */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** output nonlinear row to file stream */
extern
SCIP_RETCODE SCIPnlrowPrint(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** increases usage counter of NLP nonlinear row */
extern
void SCIPnlrowCapture(
   SCIP_NLROW*           nlrow               /**< nonlinear row to capture */
   );

/** decreases usage counter of NLP nonlinear row */
extern
SCIP_RETCODE SCIPnlrowRelease(
   SCIP_NLROW**          nlrow,              /**< nonlinear row to free */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** ensures, that linear coefficient array of nonlinear row can store at least num entries */
extern
SCIP_RETCODE SCIPnlrowEnsureLinearSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** adds a previously non existing linear coefficient to an NLP nonlinear row */
extern
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
extern
SCIP_RETCODE SCIPnlrowDelLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var                 /**< coefficient to be deleted */
   );

/** changes or adds a linear coefficient to a nonlinear row */
extern
SCIP_RETCODE SCIPnlrowChgLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< new value of coefficient */
   );

/** ensures, that quadratic variables array of nonlinear row can store at least num entries */
extern
SCIP_RETCODE SCIPnlrowEnsureQuadVarsSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** adds variable to quadvars array of row */
extern
SCIP_RETCODE SCIPnlrowAddQuadVar(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable to search for */
   );

/** ensures, that quadratic elements array of nonlinear row can store at least num entries */
extern
SCIP_RETCODE SCIPnlrowEnsureQuadElementsSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** adds a previously non existing quadratic element to an NLP nonlinear row */
extern
SCIP_RETCODE SCIPnlrowAddQuadElement(
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_QUADELEM         elem                /**< quadratic element to add */
   );

/** deletes quadratic element from nonlinear row */
extern
SCIP_RETCODE SCIPnlrowDelQuadElement(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   idx1,               /**< index of first variable in element */
   int                   idx2                /**< index of second variable in element */
   );

/** changes or adds a quadratic element to a nonlinear row */
extern
SCIP_RETCODE SCIPnlrowChgQuadElem(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_QUADELEM         elem                /**< new quadratic element */
   );

/** replaces or deletes an expression tree in nonlinear row */
extern
SCIP_RETCODE SCIPnlrowChgExprtree(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_EXPRTREE*        exprtree            /**< new expression tree, or NULL to delete current one */
   );

/** changes a parameter in an expression of a nonlinear row */
extern
SCIP_RETCODE SCIPnlrowChgExprtreeParam(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   paramidx,           /**< index of parameter in expression tree's parameter array */
   SCIP_Real             paramval            /**< new value of parameter */
   );

/** changes all parameters in an expression of a nonlinear row */
extern
SCIP_RETCODE SCIPnlrowChgExprtreeParams(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            paramvals           /**< new values of parameters */
   );

/** changes constant of nonlinear row */
extern
SCIP_RETCODE SCIPnlrowChgConstant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             constant            /**< new constant */
   );

/** changes left hand side of nonlinear row */
extern
SCIP_RETCODE SCIPnlrowChgLhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of nonlinear row */
extern
SCIP_RETCODE SCIPnlrowChgRhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             rhs                 /**< new right hand side */
   );

/** removes (or substitutes) all fixed, negated, aggregated, multi-aggregated variables from the linear, quadratic, and nonquadratic terms of a nonlinear row */
extern
SCIP_RETCODE SCIPnlrowRemoveFixedVars(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** recalculates the current activity of a nonlinear row in the current NLP solution */
extern
SCIP_RETCODE SCIPnlrowRecalcNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives the activity of a nonlinear row in the current NLP solution */
extern
SCIP_RETCODE SCIPnlrowGetNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            activity            /**< buffer to store activity value */
   );

/** gives the feasibility of a nonlinear row in the current NLP solution: negative value means infeasibility */
extern
SCIP_RETCODE SCIPnlrowGetNLPFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   );

/** calculates the current pseudo activity of a nonlinear row */
extern
SCIP_RETCODE SCIPnlrowRecalcPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** returns the pseudo activity of a nonlinear row in the current pseudo solution */
extern
SCIP_RETCODE SCIPnlrowGetPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudoactivity      /**< buffer to store pseudo activity value */
   );

/** returns the pseudo feasibility of a nonlinear row in the current pseudo solution: negative value means infeasibility */
extern
SCIP_RETCODE SCIPnlrowGetPseudoFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudofeasibility   /**< buffer to store pseudo feasibility value */
   );

/** returns the activity of a nonlinear row for a given solution */
extern
SCIP_RETCODE SCIPnlrowGetSolActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            activity            /**< buffer to store activity value */
   );

/** returns the feasibility of a nonlinear row for the given solution */
extern
SCIP_RETCODE SCIPnlrowGetSolFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   );

/** returns the minimal activity of a nonlinear row w.r.t. the variables' bounds */
extern
SCIP_RETCODE SCIPnlrowGetActivityBounds(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Real*            minactivity,        /**< buffer to store minimal activity, or NULL */
   SCIP_Real*            maxactivity         /**< buffer to store maximal activity, or NULL */
   );

/** returns whether the nonlinear row is redundant w.r.t. the variables' bounds */
extern
SCIP_RETCODE SCIPnlrowIsRedundant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
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
extern
SCIP_RETCODE SCIPnlpFree(
   SCIP_NLP**            nlp,                /**< pointer to NLP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   );

/** resets the NLP to the empty NLP by removing all variables and rows from NLP,
 *  releasing all rows, and flushing the changes to the NLP solver
 */
extern
SCIP_RETCODE SCIPnlpReset(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   );

/** currently a dummy function that always returns TRUE */
extern
SCIP_Bool SCIPnlpHasCurrentNodeNLP(
   SCIP_NLP*             nlp                 /**< NLP data */
   );

/** ensures, that variables array of NLP can store at least num entries */
extern
SCIP_RETCODE SCIPnlpEnsureVarsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** adds a variable to the NLP and captures the variable */
extern
SCIP_RETCODE SCIPnlpAddVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable */
   );

/** adds a set of variables to the NLP and captures the variables */
extern
SCIP_RETCODE SCIPnlpAddVars(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables to add */
   SCIP_VAR**            vars                /**< variables to add */
   );

/** deletes a variable from the NLP and releases the variable */
extern
SCIP_RETCODE SCIPnlpDelVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< SCIP LP, needed to release variable */
   SCIP_VAR*             var                 /**< variable */
   );

/** ensures, that nonlinear rows array of NLP can store at least num entries */
extern
SCIP_RETCODE SCIPnlpEnsureNlRowsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** adds a nonlinear row to the NLP and captures it
 * all variables of the row need to be present in the NLP */
extern
SCIP_RETCODE SCIPnlpAddNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   );

/** adds nonlinear rows to the NLP and captures them
 * all variables of the row need to be present in the NLP */
extern
SCIP_RETCODE SCIPnlpAddNlRows(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   nnlrows,            /**< number of rows to add */
   SCIP_NLROW**          nlrows              /**< rows to add */
   );

/** deletes a nonlinear row from the NLP
 * does nothing if nonlinear row is not in NLP */
extern
SCIP_RETCODE SCIPnlpDelNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   );

/** applies all cached changes to the NLP solver */
extern
SCIP_RETCODE SCIPnlpFlush(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** solves the NLP */
extern
SCIP_RETCODE SCIPnlpSolve(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** gets objective value of current NLP */
extern
SCIP_Real SCIPnlpGetObjval(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives current pseudo objective value */
extern
SCIP_RETCODE SCIPnlpGetPseudoObjval(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudoobjval        /**< buffer to store pseudo objective value */
   );

/** gets fractional variables of last NLP solution along with solution values and fractionalities
 */
extern
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
extern
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
extern
SCIP_RETCODE SCIPnlpSetInitialGuess(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_Real*            initguess           /**< new initial guess, or NULL to clear previous one */
   );

/** writes NLP to a file */
extern
SCIP_RETCODE SCIPnlpWrite(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           fname               /**< file name */
   );

/*
 * NLP diving methods
 */

/** signals start of diving */
extern
SCIP_RETCODE SCIPnlpStartDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** resets the bound and objective changes made during diving and disables diving mode */
extern
SCIP_RETCODE SCIPnlpEndDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** changes coefficient of variable in diving NLP */
extern
SCIP_RETCODE SCIPnlpChgVarObjDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             coef                /**< new linear coefficient of variable in objective */
   );

/** changes bounds of variable in diving NLP */
extern
SCIP_RETCODE SCIPnlpChgVarBoundsDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable which bounds to change */
   SCIP_Real             lb,                 /**< new lower bound of variable */
   SCIP_Real             ub                  /**< new upper bound of variable */
   );

/** changes bounds of a set of variables in diving NLP */
extern
SCIP_RETCODE SCIPnlpChgVarsBoundsDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables which bounds to change */
   SCIP_VAR**            vars,               /**< variables which bounds to change */
   SCIP_Real*            lbs,                /**< new lower bounds of variables */
   SCIP_Real*            ubs                 /**< new upper bounds of variables */
   );

/** returns whether the objective function has been changed during diving */
extern
SCIP_Bool SCIPnlpIsDivingObjChanged(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** solves diving NLP */
extern
SCIP_RETCODE SCIPnlpSolveDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** gets array with variables of the NLP */
extern
SCIP_VAR** SCIPnlpGetVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets current number of variables in NLP */
extern
int SCIPnlpGetNVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** computes for each variables the number of NLP rows in which the variable appears in a nonlinear var */
extern
SCIP_RETCODE SCIPnlpGetVarsNonlinearity(
   SCIP_NLP*             nlp,                /**< current NLP data */
   int*                  nlcount             /**< an array of length at least SCIPnlpGetNVars() to store nonlinearity counts of variables */
   );

/** indicates whether there exists a row that contains a continuous variable in a nonlinear term
 *
 * @note The method may have to touch every row and nonlinear term to compute its result.
 */
extern
SCIP_Bool SCIPnlpHasContinuousNonlinearity(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives dual solution values associated with lower bounds of NLP variables */
extern
SCIP_Real* SCIPnlpGetVarsLbDualsol(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives dual solution values associated with upper bounds of NLP variables */
extern
SCIP_Real* SCIPnlpGetVarsUbDualsol(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets array with nonlinear rows of the NLP */
extern
SCIP_NLROW** SCIPnlpGetNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets current number of nonlinear rows in NLP */
extern
int SCIPnlpGetNNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets the NLP solver interface */
extern
SCIP_NLPI* SCIPnlpGetNLPI(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets the NLP problem in the solver interface */
extern
SCIP_NLPIPROBLEM* SCIPnlpGetNLPIProblem(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** indicates whether NLP is currently in diving mode */
extern
SCIP_Bool SCIPnlpIsDiving(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets solution status of current NLP */
extern
SCIP_NLPSOLSTAT SCIPnlpGetSolstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets termination status of last NLP solve */
extern
SCIP_NLPTERMSTAT SCIPnlpGetTermstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gives statistics (number of iterations, solving time, ...) of last NLP solve */
extern
SCIP_RETCODE SCIPnlpGetStatistics(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   );

/** indicates whether a feasible solution for the current NLP is available
 * thus, returns whether the solution status <= feasible  */
extern
SCIP_Bool SCIPnlpHasSolution(
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** gets integer parameter of NLP */
extern
SCIP_RETCODE SCIPnlpGetIntPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   );

/** sets integer parameter of NLP */
extern
SCIP_RETCODE SCIPnlpSetIntPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** gets floating point parameter of NLP */
extern
SCIP_RETCODE SCIPnlpGetRealPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   );

/** sets floating point parameter of NLP */
extern
SCIP_RETCODE SCIPnlpSetRealPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   );

/** gets string parameter of NLP */
extern
SCIP_RETCODE SCIPnlpGetStringPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value */
   );

/** sets string parameter of NLP */
extern
SCIP_RETCODE SCIPnlpSetStringPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLP_H__ */
