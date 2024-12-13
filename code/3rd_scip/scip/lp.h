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

/**@file   lp.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for LP management
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LP_H__
#define __SCIP_LP_H__


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
#include "scip/type_branch.h"
#include "scip/pub_lp.h"

#include "scip/struct_lp.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Column methods
 */

/** creates an LP column */
SCIP_RETCODE SCIPcolCreate(
   SCIP_COL**            col,                /**< pointer to column data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable, this column represents */
   int                   len,                /**< number of nonzeros in the column */
   SCIP_ROW**            rows,               /**< array with rows of column entries */
   SCIP_Real*            vals,               /**< array with coefficients of column entries */
   SCIP_Bool             removable           /**< should the column be removed from the LP due to aging or cleanup? */
   );

/** frees an LP column */
SCIP_RETCODE SCIPcolFree(
   SCIP_COL**            col,                /**< pointer to LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** output column to file stream */
void SCIPcolPrint(
   SCIP_COL*             col,                /**< LP column */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** adds a previously non existing coefficient to an LP column */
SCIP_RETCODE SCIPcolAddCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val                 /**< value of coefficient */
   );

/** deletes coefficient from column */
SCIP_RETCODE SCIPcolDelCoef(
   SCIP_COL*             col,                /**< column to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< coefficient to be deleted */
   );

/** changes or adds a coefficient to an LP column */
SCIP_RETCODE SCIPcolChgCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val                 /**< value of coefficient */
   );

/** increases value of an existing or nonexisting coefficient in an LP column */
SCIP_RETCODE SCIPcolIncCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             incval              /**< value to add to the coefficient */
   );

/** changes objective value of column */
SCIP_RETCODE SCIPcolChgObj(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newobj              /**< new objective value */
   );

/** changes lower bound of column */
SCIP_RETCODE SCIPcolChgLb(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newlb               /**< new lower bound value */
   );

/** changes upper bound of column */
SCIP_RETCODE SCIPcolChgUb(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newub               /**< new upper bound value */
   );

/** calculates the reduced costs of a column using the given dual solution vector */
SCIP_Real SCIPcolCalcRedcost(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            dualsol             /**< dual solution vector for current LP rows */
   );

/** gets the reduced costs of a column in last LP or after recalculation */
SCIP_Real SCIPcolGetRedcost(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets the feasibility of (the dual row of) a column in last LP or after recalculation */
SCIP_Real SCIPcolGetFeasibility(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** calculates the Farkas coefficient y^T A_i of a column i using the given dual Farkas vector y */
SCIP_Real SCIPcolCalcFarkasCoef(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            dualfarkas          /**< dense dual Farkas vector for current LP rows */
   );

/** gets the Farkas coefficient y^T A_i of a column i in last LP (which must be infeasible) */
SCIP_Real SCIPcolGetFarkasCoef(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets the Farkas value of a column in last LP (which must be infeasible), i.e. the Farkas coefficient y^T A_i times
 *  the best bound for this coefficient, i.e. max{y^T A_i x_i | lb <= x_i <= ub}
 */
SCIP_Real SCIPcolGetFarkasValue(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpStartStrongbranch(
   SCIP_LP*              lp                  /**< LP data */
   );

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpEndStrongbranch(
   SCIP_LP*              lp                  /**< LP data */
   );

/** sets strong branching information for a column variable */
void SCIPcolSetStrongbranchData(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real             lpobjval,           /**< objective value of the current LP */
   SCIP_Real             primsol,            /**< primal solution value of the column in the current LP */
   SCIP_Real             sbdown,             /**< dual bound after branching column down */
   SCIP_Real             sbup,               /**< dual bound after branching column up */
   SCIP_Bool             sbdownvalid,        /**< is the returned down value a valid dual bound? */
   SCIP_Bool             sbupvalid,          /**< is the returned up value a valid dual bound? */
   SCIP_Longint          iter,               /**< total number of strong branching iterations */
   int                   itlim               /**< iteration limit applied to the strong branching call */
   );

/** invalidates strong branching information for a column variable */
void SCIPcolInvalidateStrongbranchData(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp                  /**< LP data */
   );

/** gets strong branching information on a column variable */
SCIP_RETCODE SCIPcolGetStrongbranch(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Bool             integral,           /**< should integral strong branching be performed? */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Bool             updatecol,          /**< should col be updated, or should it stay in its current state ? */
   SCIP_Bool             updatestat,         /**< should stat be updated, or should it stay in its current state ? */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   );

/** gets strong branching information on column variables */
SCIP_RETCODE SCIPcolGetStrongbranches(
   SCIP_COL**            cols,               /**< LP columns */
   int                   ncols,              /**< number of columns */
   SCIP_Bool             integral,           /**< should integral strong branching be performed? */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   );

/** gets last strong branching information available for a column variable;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given column;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
void SCIPcolGetStrongbranchLast(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            down,               /**< stores dual bound after branching column down, or NULL */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up, or NULL */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Real*            solval,             /**< stores LP solution value of column at last strong branching call, or NULL */
   SCIP_Real*            lpobjval            /**< stores LP objective value at last strong branching call, or NULL */
   );

/** if strong branching was already applied on the column at the current node, returns the number of LPs solved after
 *  the LP where the strong branching on this column was applied;
 *  if strong branching was not yet applied on the column at the current node, returns INT_MAX
 */
SCIP_Longint SCIPcolGetStrongbranchLPAge(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** marks a column to be not removable from the LP in the current node because it became obsolete */
void SCIPcolMarkNotRemovableLocal(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat                /**< problem statistics */
   );


/*
 * Row methods
 */

/** creates and captures an LP row */
SCIP_RETCODE SCIProwCreate(
   SCIP_ROW**            row,                /**< pointer to LP row data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_ROWORIGINTYPE    origintype,         /**< type of origin of row */
   void*                 origin,             /**< pointer to constraint handler or separator who created the row (NULL if unkown) */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   );

/** frees an LP row */
SCIP_RETCODE SCIProwFree(
   SCIP_ROW**            row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** output row to file stream */
void SCIProwPrint(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** ensures, that column array of row can store at least num entries */
SCIP_RETCODE SCIProwEnsureSize(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   );

/** increases usage counter of LP row */
void SCIProwCapture(
   SCIP_ROW*             row                 /**< LP row */
   );

/** decreases usage counter of LP row, and frees memory if necessary */
SCIP_RETCODE SCIProwRelease(
   SCIP_ROW**            row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** enables delaying of row sorting */
void SCIProwDelaySort(
   SCIP_ROW*             row                 /**< LP row */
   );

/** disables delaying of row sorting, sorts row and merges coefficients with equal columns */
void SCIProwForceSort(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** adds a previously non existing coefficient to an LP row */
SCIP_RETCODE SCIProwAddCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val                 /**< value of coefficient */
   );

/** deletes coefficient from row */
SCIP_RETCODE SCIProwDelCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col                 /**< coefficient to be deleted */
   );

/** changes or adds a coefficient to an LP row */
SCIP_RETCODE SCIProwChgCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val                 /**< value of coefficient */
   );

/** increases value of an existing or nonexisting coefficient in an LP column */
SCIP_RETCODE SCIProwIncCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             incval              /**< value to add to the coefficient */
   );

/** changes constant value of a row */
SCIP_RETCODE SCIProwChgConstant(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             constant            /**< new constant value */
   );

/** add constant value to a row */
SCIP_RETCODE SCIProwAddConstant(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             addval              /**< constant value to add to the row */
   );

/** changes left hand side of LP row */
SCIP_RETCODE SCIProwChgLhs(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of LP row */
SCIP_RETCODE SCIProwChgRhs(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             rhs                 /**< new right hand side */
   );

/** changes the local flag of LP row */
SCIP_RETCODE SCIProwChgLocal(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Bool             local               /**< new value for local flag */
   );

/** tries to find a value, such that all row coefficients, if scaled with this value become integral */
SCIP_RETCODE SCIProwCalcIntegralScalar(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   );

/** tries to scale row, s.t. all coefficients become integral */
SCIP_RETCODE SCIProwMakeIntegral(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal value to scale row with */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Bool*            success             /**< stores whether row could be made rational */
   );

/** recalculates the current activity of a row */
void SCIProwRecalcLPActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** returns the activity of a row in the current LP solution */
SCIP_Real SCIProwGetLPActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** returns the feasibility of a row in the current LP solution: negative value means infeasibility */
SCIP_Real SCIProwGetLPFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** returns the feasibility of a row in the current relaxed solution: negative value means infeasibility */
SCIP_Real SCIProwGetRelaxFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** returns the feasibility of a row in the current NLP solution: negative value means infeasibility */
SCIP_Real SCIProwGetNLPFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** calculates the current pseudo activity of a row */
void SCIProwRecalcPseudoActivity(
   SCIP_ROW*             row,                /**< row data */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** returns the pseudo activity of a row in the current pseudo solution */
SCIP_Real SCIProwGetPseudoActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** returns the pseudo feasibility of a row in the current pseudo solution: negative value means infeasibility */
SCIP_Real SCIProwGetPseudoFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** returns the activity of a row for a given solution */
SCIP_Real SCIProwGetSolActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns the feasibility of a row for the given solution */
SCIP_Real SCIProwGetSolFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns the minimal activity of a row w.r.t. the columns' bounds */
SCIP_Real SCIProwGetMinActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** returns the maximal activity of a row w.r.t. the columns' bounds */
SCIP_Real SCIProwGetMaxActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** returns whether the row is unmodifiable and redundant w.r.t. the columns' bounds */
SCIP_Bool SCIProwIsRedundant(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** gets maximal absolute value of row vector coefficients */
SCIP_Real SCIProwGetMaxval(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets minimal absolute value of row vector's non-zero coefficients */
SCIP_Real SCIProwGetMinval(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets maximal column index of row entries */
int SCIProwGetMaxidx(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets minimal column index of row entries */
int SCIProwGetMinidx(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets number of integral columns in row */
int SCIProwGetNumIntCols(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns row's cutoff distance in the direction of the given primal solution */
SCIP_Real SCIProwGetLPSolCutoffDistance(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< solution to compute direction for cutoff distance; must not be NULL */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** returns row's efficacy with respect to the current LP solution: e = -feasibility/norm */
SCIP_Real SCIProwGetLPEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** returns whether the row's efficacy with respect to the current LP solution is greater than the minimal cut efficacy */
SCIP_Bool SCIProwIsLPEfficacious(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             root                /**< should the root's minimal cut efficacy be used? */
   );

/** returns row's efficacy with respect to the given primal solution: e = -feasibility/norm */
SCIP_Real SCIProwGetSolEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns whether the row's efficacy with respect to the given primal solution is greater than the minimal cut
 *  efficacy
 */
SCIP_Bool SCIProwIsSolEfficacious(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             root                /**< should the root's minimal cut efficacy be used? */
   );

/** returns row's efficacy with respect to the relaxed solution: e = -feasibility/norm */
SCIP_Real SCIProwGetRelaxEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** returns row's efficacy with respect to the NLP solution: e = -feasibility/norm */
SCIP_Real SCIProwGetNLPEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   );

/** gets parallelism of row with objective function: if the returned value is 1, the row is parallel to the objective
 *  function, if the value is 0, it is orthogonal to the objective function
 */
SCIP_Real SCIProwGetObjParallelism(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** includes event handler with given data in row's event filter */
SCIP_RETCODE SCIProwCatchEvent(
   SCIP_ROW*             row,                /**< row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   );

/** deletes event handler with given data from row's event filter */
SCIP_RETCODE SCIProwDropEvent(
   SCIP_ROW*             row,                /**< row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int                   filterpos           /**< position of event filter entry returned by SCIPvarCatchEvent(), or -1 */
   );

/** marks a row to be not removable from the LP in the current node */
void SCIProwMarkNotRemovableLocal(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   );


/*
 * LP methods
 */

/** creates empty LP data object */
SCIP_RETCODE SCIPlpCreate(
   SCIP_LP**             lp,                 /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name                /**< problem name */
   );

/** frees LP data object */
SCIP_RETCODE SCIPlpFree(
   SCIP_LP**             lp,                 /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** resets the LP to the empty LP by removing all columns and rows from LP, releasing all rows, and flushing the
 *  changes to the LP solver
 */
SCIP_RETCODE SCIPlpReset(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** adds a column to the LP and captures the variable */
SCIP_RETCODE SCIPlpAddCol(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col,                /**< LP column */
   int                   depth               /**< depth in the tree where the column addition is performed */
   );

/** adds a row to the LP and captures it */
SCIP_RETCODE SCIPlpAddRow(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_ROW*             row,                /**< LP row */
   int                   depth               /**< depth in the tree where the row addition is performed */
   );

/** removes all columns after the given number of columns from the LP */
SCIP_RETCODE SCIPlpShrinkCols(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newncols            /**< new number of columns in the LP */
   );

/** removes and releases all rows after the given number of rows from the LP */
SCIP_RETCODE SCIPlpShrinkRows(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   int                   newnrows            /**< new number of rows in the LP */
   );

/** removes all columns and rows from LP, releases all rows */
SCIP_RETCODE SCIPlpClear(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** remembers number of columns and rows to track the newly added ones */
void SCIPlpMarkSize(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** sets the remembered number of columns and rows to the given values */
void SCIPlpSetSizeMark(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   nrows,              /**< number of rows to set the size marker to */
   int                   ncols               /**< number of columns to set the size marker to */
   );

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
SCIP_RETCODE SCIPlpGetBasisInd(
   SCIP_LP*              lp,                 /**< LP data */
   int*                  basisind            /**< pointer to store basis indices ready to keep number of rows entries */
   );

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpGetBase(
   SCIP_LP*              lp,                 /**< LP data */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   );

/** gets a row from the inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpGetBInvRow(
   SCIP_LP*              lp,                 /**< LP data */
   int                   r,                  /**< row number */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   );

/** gets a column from the inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpGetBInvCol(
   SCIP_LP*              lp,                 /**< LP data */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP
                                              *   returned by SCIPcolGetLPPos(); you have to call SCIPgetBasisInd()
                                              *   to get the array which links the B^-1 column numbers to the row and
                                              *   column numbers of the LP! c must be between 0 and nrows-1, since the
                                              *   basis has the size nrows * nrows */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   );

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
SCIP_RETCODE SCIPlpGetBInvARow(
   SCIP_LP*              lp,                 /**< LP data */
   int                   r,                  /**< row number */
   SCIP_Real*            binvrow,            /**< row in B^-1 from prior call to SCIPlpGetBInvRow(), or NULL */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   );

/** gets a column from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A),
 *  i.e., it computes B^-1 * A_c with A_c being the c'th column of A
 */
SCIP_RETCODE SCIPlpGetBInvACol(
   SCIP_LP*              lp,                 /**< LP data */
   int                   c,                  /**< column number which can be accessed by SCIPcolGetLPPos() */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   );

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
SCIP_RETCODE SCIPlpSumRows(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   SCIP_Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   SCIP_Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   );

/** stores LP state (like basis information) into LP state object */
SCIP_RETCODE SCIPlpGetState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** loads LP state (like basis information) into solver */
SCIP_RETCODE SCIPlpSetState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPISTATE*        lpistate,           /**< LP state information (like basis information) */
   SCIP_Bool             wasprimfeas,        /**< primal feasibility when LP state information was stored */
   SCIP_Bool             wasprimchecked,     /**< true if the LP solution has passed the primal feasibility check */
   SCIP_Bool             wasdualfeas,        /**< dual feasibility when LP state information was stored */
   SCIP_Bool             wasdualchecked      /**< true if the LP solution has passed the dual feasibility check */
   );

/** frees LP state information */
SCIP_RETCODE SCIPlpFreeState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** interrupts the currently ongoing lp solve or disables the interrupt */
SCIP_RETCODE SCIPlpInterrupt(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             interrupt           /**< TRUE if interrupt should be set, FALSE if it should be disabled */
   );

/** stores pricing norms into LP norms object */
SCIP_RETCODE SCIPlpGetNorms(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LP pricing norms information */
   );

/** loads pricing norms from LP norms object into solver */
SCIP_RETCODE SCIPlpSetNorms(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS*        lpinorms            /**< LP pricing norms information */
   );

/** frees pricing norms information */
SCIP_RETCODE SCIPlpFreeNorms(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LP pricing norms information */
   );

/** return the current cutoff bound of the lp */
SCIP_Real SCIPlpGetCutoffbound(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** sets the upper objective limit of the LP solver */
SCIP_RETCODE SCIPlpSetCutoffbound(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             cutoffbound         /**< new upper objective limit */
   );

/** gets current primal feasibility tolerance of LP solver */
SCIP_Real SCIPlpGetFeastol(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** sets primal feasibility tolerance of LP solver */
void SCIPlpSetFeastol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newfeastol          /**< new primal feasibility tolerance for LP */
   );

/** resets primal feasibility tolerance of LP solver
 *
 * Sets primal feasibility tolerance to min of numerics/lpfeastolfactor * numerics/feastol and relaxfeastol.
 */
void SCIPlpResetFeastol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** applies all cached changes to the LP solver */
SCIP_RETCODE SCIPlpFlush(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** marks the LP to be flushed, even if the LP thinks it is not flushed */
SCIP_RETCODE SCIPlpMarkFlushed(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** solves the LP with simplex algorithm, and copy the solution into the column's data */
SCIP_RETCODE SCIPlpSolveAndEval(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Longint          itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool             limitresolveiters,  /**< should LP iterations for resolving calls be limited?
                                              *   (limit is computed within the method w.r.t. the average LP iterations) */
   SCIP_Bool             aging,              /**< should aging and removal of obsolete cols/rows be applied? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   );

/** gets solution status of current LP */
SCIP_LPSOLSTAT SCIPlpGetSolstat(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** sets whether the root LP is a relaxation of the problem and its optimal objective value is a global lower bound */
void SCIPlpSetRootLPIsRelax(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             isrelax             /**< is the root lp a relaxation of the problem? */
   );

/** returns whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound */
SCIP_Bool SCIPlpIsRootLPRelax(
   SCIP_LP*              lp                  /**< LP data */
   );

/** gets objective value of current LP
 *
 *  @note This method returns the objective value of the current LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the LP solution status is
 *        SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 */
SCIP_Real SCIPlpGetObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   );

/** gets part of objective value of current LP that results from COLUMN variables only */
SCIP_Real SCIPlpGetColumnObjval(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets part of objective value of current LP that results from LOOSE variables only */
SCIP_Real SCIPlpGetLooseObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   );

/** remembers the current LP objective value as root solution value */
void SCIPlpStoreRootObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   );

/** invalidates the root LP solution value */
void SCIPlpInvalidateRootObjval(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets the global pseudo objective value; that is all variables set to their best (w.r.t. the objective function)
 *  global bound
 */
SCIP_Real SCIPlpGetGlobalPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   );

/** recomputes local and global pseudo objective values */
void SCIPlpRecomputeLocalAndGlobalPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   );

/** gets the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound
 */
SCIP_Real SCIPlpGetPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   );

/** gets pseudo objective value, if a bound of the given variable would be modified in the given way */
SCIP_Real SCIPlpGetModifiedPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** gets pseudo objective value, if a bound of the given variable would be modified in the given way;
 *  perform calculations with interval arithmetic to get an exact lower bound
 */
SCIP_Real SCIPlpGetModifiedProvedPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** updates current pseudo and loose objective value for a change in a variable's objective coefficient */
SCIP_RETCODE SCIPlpUpdateVarObj(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldobj,             /**< old objective coefficient of variable */
   SCIP_Real             newobj              /**< new objective coefficient of variable */
   );

/** updates current root pseudo objective value for a global change in a variable's lower bound */
SCIP_RETCODE SCIPlpUpdateVarLbGlobal(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb               /**< new lower bound of variable */
   );

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
SCIP_RETCODE SCIPlpUpdateVarLb(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb               /**< new lower bound of variable */
   );

/** updates current root pseudo objective value for a global change in a variable's upper bound */
SCIP_RETCODE SCIPlpUpdateVarUbGlobal(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub               /**< new upper bound of variable */
   );

/** updates current pseudo objective value for a change in a variable's upper bound */
SCIP_RETCODE SCIPlpUpdateVarUb(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub               /**< new upper bound of variable */
   );

/** informs LP, that given variable was added to the problem */
SCIP_RETCODE SCIPlpUpdateAddVar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that is now a LOOSE problem variable */
   );

/** informs LP, that given variable is to be deleted from the problem */
SCIP_RETCODE SCIPlpUpdateDelVar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that will be deleted from the problem */
   );

/** informs LP, that given formerly loose problem variable is now a column variable */
SCIP_RETCODE SCIPlpUpdateVarColumn(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   );

/** informs LP, that given formerly column problem variable is now again a loose variable */
SCIP_RETCODE SCIPlpUpdateVarLoose(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   );

/** decrease the number of loose variables by one */
void SCIPlpDecNLoosevars(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** stores the LP solution in the columns and rows */
SCIP_RETCODE SCIPlpGetSol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            dualfeasible        /**< pointer to store whether the solution is dual feasible, or NULL */
   );

/** stores LP solution with infinite objective value in the columns and rows */
SCIP_RETCODE SCIPlpGetUnboundedSol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            rayfeasible         /**< pointer to store whether the primal ray is a feasible unboundedness proof, or NULL */
   );

/** returns primal ray proving the unboundedness of the current LP */
SCIP_RETCODE SCIPlpGetPrimalRay(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            ray                 /**< array for storing primal ray values, they are stored w.r.t. the problem index of the variables,
                                              *   so the size of this array should be at least number of active variables
                                              *   (all entries have to be initialized to 0 before) */
   );

/** stores the dual Farkas multipliers for infeasibility proof in rows. besides, the proof is checked for validity if
 *  lp/checkfarkas = TRUE.
 *
 *  @note the check will not be performed if @p valid is NULL.
 */
SCIP_RETCODE SCIPlpGetDualfarkas(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            valid               /**< pointer to store whether the Farkas proof is valid  or NULL */
   );

/** get number of iterations used in last LP solve */
SCIP_RETCODE SCIPlpGetIterations(
   SCIP_LP*              lp,                 /**< current LP data */
   int*                  iterations          /**< pointer to store the iteration count */
   );

/** increases age of columns with solution value 0.0 and rows with activity not at its bounds,
 *  resets age of non-zero columns and sharp rows
 */
SCIP_RETCODE SCIPlpUpdateAges(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** removes all non-basic columns and basic rows in the part of the LP created at the current node, that are too old */
SCIP_RETCODE SCIPlpRemoveNewObsoletes(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** removes all non-basic columns and basic rows in whole LP, that are too old */
SCIP_RETCODE SCIPlpRemoveAllObsoletes(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** removes all non-basic columns at 0.0 and basic rows in the part of the LP created at the current node */
SCIP_RETCODE SCIPlpCleanupNew(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_Bool             root                /**< are we at the root node? */
   );

/** removes all non-basic columns at 0.0 and basic rows in the whole LP */
SCIP_RETCODE SCIPlpCleanupAll(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_Bool             root                /**< are we at the root node? */
   );

/** removes all redundant rows that were added at the current node */
SCIP_RETCODE SCIPlpRemoveRedundantRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** initiates LP diving */
SCIP_RETCODE SCIPlpStartDive(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
SCIP_RETCODE SCIPlpEndDive(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_VAR**            vars,               /**< array with all active variables */
   int                   nvars               /**< number of active variables */
   );

/** records a current row side such that any change will be undone after diving */
SCIP_RETCODE SCIPlpRecordOldRowSideDive(
   SCIP_LP*              lp,                 /**< LP data object */
   SCIP_ROW*             row,                /**< row affected by the change */
   SCIP_SIDETYPE         sidetype            /**< side type */
   );

/** informs the LP that probing mode was initiated */
SCIP_RETCODE SCIPlpStartProbing(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** informs the LP that probing mode was finished */
SCIP_RETCODE SCIPlpEndProbing(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** informs the LP that the probing mode is now used for strongbranching */
void SCIPlpStartStrongbranchProbing(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** informs the LP that the probing mode is not used for strongbranching anymore */
void SCIPlpEndStrongbranchProbing(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets proven lower (dual) bound of last LP solution */
SCIP_RETCODE SCIPlpGetProvedLowerbound(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            bound               /**< pointer to store proven dual bound */
   );

/** gets proven dual bound of last LP solution */
SCIP_RETCODE SCIPlpIsInfeasibilityProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool*            proved              /**< pointer to store whether infeasibility is proven */
   );

/** writes LP to a file */
SCIP_RETCODE SCIPlpWrite(
   SCIP_LP*              lp,                 /**< current LP data */
   const char*           fname               /**< file name */
   );

/** writes MIP to a file */
SCIP_RETCODE SCIPlpWriteMip(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           fname,              /**< file name */
   SCIP_Bool             genericnames,       /**< should generic names like x_i and row_j be used in order to avoid
                                              *   troubles with reserved symbols? */
   SCIP_Bool             origobj,            /**< should the original objective function be used? */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< objective scaling factor */
   SCIP_Real             objoffset,          /**< objective offset, e.g., caused by variable fixings in presolving */
   SCIP_Bool             lazyconss           /**< output removable rows as lazy constraints? */
   );

/** recalculates Euclidean norm of objective function vector of column variables if it have gotten unreliable during calculation */
void SCIPlpRecalculateObjSqrNorm(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< LP data */
   );

/** compute relative interior point */
SCIP_RETCODE SCIPlpComputeRelIntPoint(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             relaxrows,          /**< should the rows be relaxed */
   SCIP_Bool             inclobjcutoff,      /**< should a row for the objective cutoff be included */
   SCIP_Real             timelimit,          /**< time limit for LP solver */
   int                   iterlimit,          /**< iteration limit for LP solver */
   SCIP_Real*            point,              /**< array to store relative interior point on exit */
   SCIP_Bool*            success             /**< buffer to indicate whether interior point was successfully computed */
   );

/** computes two measures for dual degeneracy (dual degeneracy rate and variable-constraint ratio)
 *  based on the changes applied when reducing the problem to the optimal face
 *
 *  returns the dual degeneracy rate, i.e., the share of nonbasic variables with reduced cost 0
 *  and the variable-constraint ratio, i.e., the number of unfixed variables in relation to the basis size
 */
SCIP_RETCODE SCIPlpGetDualDegeneracy(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            degeneracy,         /**< pointer to store the dual degeneracy rate */
   SCIP_Real*            varconsratio        /**< pointer to store the variable-constraint ratio */
   );

/** gets array with columns of the LP */
SCIP_COL** SCIPlpGetCols(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets current number of columns in LP */
int SCIPlpGetNCols(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets current number of unfixed columns in LP */
int SCIPlpGetNUnfixedCols(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             eps                 /**< numerical tolerance */
   );

/** gets array with rows of the LP */
SCIP_ROW** SCIPlpGetRows(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets current number of rows in LP */
int SCIPlpGetNRows(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets array with newly added columns after the last mark */
SCIP_COL** SCIPlpGetNewcols(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets number of newly added columns after the last mark */
int SCIPlpGetNNewcols(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets array with newly added rows after the last mark */
SCIP_ROW** SCIPlpGetNewrows(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets number of newly added rows after the last mark */
int SCIPlpGetNNewrows(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** gets Euclidean norm of objective function vector of column variables, only use this method if
 *  lp->objsqrnormunreliable == FALSE, so probably you have to call SCIPlpRecalculateObjSqrNorm before */
SCIP_Real SCIPlpGetObjNorm(
   SCIP_LP*              lp                  /**< LP data */
   );

/** gets the objective value of the root node LP; returns SCIP_INVALID if the root node LP was not (yet) solved */
SCIP_Real SCIPlpGetRootObjval(
   SCIP_LP*              lp                  /**< LP data */
   );

/** gets part of the objective value of the root node LP that results from COLUMN variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 */
SCIP_Real SCIPlpGetRootColumnObjval(
   SCIP_LP*              lp                  /**< LP data */
   );

/** gets part of the objective value of the root node LP that results from LOOSE variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 */
SCIP_Real SCIPlpGetRootLooseObjval(
   SCIP_LP*              lp                  /**< LP data */
   );

/** gets the LP solver interface */
SCIP_LPI* SCIPlpGetLPI(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** sets whether the current lp is a relaxation of the current problem and its optimal objective value is a local lower bound */
void SCIPlpSetIsRelax(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             relax               /**< is the current lp a relaxation? */
   );

/** returns whether the current LP is a relaxation of the problem for which it has been solved and its 
 *  solution value a valid local lower bound? 
 */
SCIP_Bool SCIPlpIsRelax(
   SCIP_LP*              lp                  /**< LP data */
   );

/** returns whether the current LP is flushed and solved */
SCIP_Bool SCIPlpIsSolved(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** return whether the current LP solution passed the primal feasibility check */
SCIP_Bool SCIPlpIsPrimalReliable(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** return whether the current LP solution passed the dual feasibility check */
SCIP_Bool SCIPlpIsDualReliable(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** returns whether the current LP solution is a basic solution */
SCIP_Bool SCIPlpIsSolBasic(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** returns whether the LP is in diving mode */
SCIP_Bool SCIPlpDiving(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** returns whether the LP is in diving mode and the objective value of at least one column was changed */
SCIP_Bool SCIPlpDivingObjChanged(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** marks the diving LP to have a changed objective function */
void SCIPlpMarkDivingObjChanged(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** marks the diving LP to not have a changed objective function anymore */
void SCIPlpUnmarkDivingObjChanged(
   SCIP_LP*              lp                  /**< current LP data */
   );

/* returns TRUE if at least one left/right hand side of an LP row was changed during diving mode */
SCIP_Bool SCIPlpDivingRowsChanged(
   SCIP_LP*              lp                  /**< current LP data */
   );

/** checks, if absolute difference of values is in range of LP primal feastol */
SCIP_Bool SCIPlpIsFeasEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if absolute difference of val1 and val2 is lower than LP primal feastol */
SCIP_Bool SCIPlpIsFeasLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if absolute difference of val1 and val2 is not greater than LP primal feastol */
SCIP_Bool SCIPlpIsFeasLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if absolute difference of val1 and val2 is greater than LP primal feastol */
SCIP_Bool SCIPlpIsFeasGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if absolute difference of val1 and val2 is not lower than -LP primal feastol */
SCIP_Bool SCIPlpIsFeasGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range LP primal feasibility tolerance of 0.0 */
SCIP_Bool SCIPlpIsFeasZero(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than LP primal feasibility tolerance */
SCIP_Bool SCIPlpIsFeasPositive(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -LP primal feasibility tolerance */
SCIP_Bool SCIPlpIsFeasNegative(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             val                 /**< value to be compared against zero */
   );


#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPlpGetCols(lp)               ((lp)->cols)
#define SCIPlpGetNCols(lp)              ((lp)->ncols)
#define SCIPlpGetRows(lp)               ((lp)->rows)
#define SCIPlpGetNRows(lp)              ((lp)->nrows)
#define SCIPlpGetNewcols(lp)            (&((lp)->cols[(lp)->firstnewcol]))
#define SCIPlpGetNNewcols(lp)           ((lp)->ncols - (lp)->firstnewcol)
#define SCIPlpGetNewrows(lp)            (&((lp)->rows[(lp)->firstnewrow]))
#define SCIPlpGetNNewrows(lp)           ((lp)->nrows - (lp)->firstnewrow)
#define SCIPlpGetObjNorm(lp)            (SQRT((lp)->objsqrnorm))
#define SCIPlpGetRootObjval(lp)         (MIN((lp)->rootlpobjval + (lp)->rootlooseobjval, SCIP_INVALID))
#define SCIPlpGetRootColumnObjval(lp)   ((lp)->rootlpobjval)
#define SCIPlpGetRootLooseObjval(lp)    ((lp)->rootlooseobjval)
#define SCIPlpGetLPI(lp)                (lp)->lpi
#define SCIPlpSetIsRelax(lp,relax)      ((lp)->isrelax = relax)
#define SCIPlpIsRelax(lp)               (lp)->isrelax
#define SCIPlpIsSolved(lp)              ((lp)->flushed && (lp)->solved)
#define SCIPlpIsSolBasic(lp)            ((lp)->solisbasic)
#define SCIPlpDiving(lp)                (lp)->diving
#define SCIPlpDivingObjChanged(lp)      (lp)->divingobjchg
#define SCIPlpMarkDivingObjChanged(lp)  ((lp)->divingobjchg = TRUE)
#define SCIPlpUnmarkDivingObjChanged(lp) ((lp)->divingobjchg = FALSE)
#define SCIPlpDivingRowsChanged(lp)     ((lp)->ndivechgsides > 0)

#define SCIPlpIsFeasEQ(set, lp, val1, val2) ( EPSEQ(val1, val2, (lp)->feastol) )
#define SCIPlpIsFeasLT(set, lp, val1, val2) ( EPSLT(val1, val2, (lp)->feastol) )
#define SCIPlpIsFeasLE(set, lp, val1, val2) ( EPSLE(val1, val2, (lp)->feastol) )
#define SCIPlpIsFeasGT(set, lp, val1, val2) ( EPSGT(val1, val2, (lp)->feastol) )
#define SCIPlpIsFeasGE(set, lp, val1, val2) ( EPSGE(val1, val2, (lp)->feastol) )
#define SCIPlpIsFeasZero(lp, val)        ( EPSZ(val, (lp)->feastol) )
#define SCIPlpIsFeasPositive(lp, val)    ( EPSP(val, (lp)->feastol) )
#define SCIPlpIsFeasNegative(lp, val)    ( EPSN(val, (lp)->feastol) )

#endif

#ifdef __cplusplus
}
#endif

#endif
