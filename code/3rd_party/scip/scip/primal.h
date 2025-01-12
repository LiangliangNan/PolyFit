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

/**@file   primal.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRIMAL_H__
#define __SCIP_PRIMAL_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/type_heur.h"

#include "scip/struct_primal.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates primal data */
SCIP_RETCODE SCIPprimalCreate(
   SCIP_PRIMAL**         primal              /**< pointer to primal data */
   );

/** frees primal data */
SCIP_RETCODE SCIPprimalFree(
   SCIP_PRIMAL**         primal,             /**< pointer to primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** clears primal data */
SCIP_RETCODE SCIPprimalClear(
   SCIP_PRIMAL**         primal,             /**< pointer to primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** sets the cutoff bound in primal data and in LP solver */
SCIP_RETCODE SCIPprimalSetCutoffbound(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound,        /**< new cutoff bound */
   SCIP_Bool             useforobjlimit      /**< should the cutoff bound be used to update the objective limit, if
                                              *   better? */
   );

/** sets upper bound in primal data and in LP solver */
SCIP_RETCODE SCIPprimalSetUpperbound(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             upperbound          /**< new upper bound */
   );

/** updates upper bound and cutoff bound in primal data after a tightening of the problem's objective limit */
SCIP_RETCODE SCIPprimalUpdateObjlimit(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** recalculates upper bound and cutoff bound in primal data after a change of the problem's objective offset */
SCIP_RETCODE SCIPprimalUpdateObjoffset(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** adds additional objective offset in origanal space to all existing solution (in original space) */
void SCIPprimalAddOrigObjoffset(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             addval              /**< additional objective offset in original space */
   );

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 */
SCIP_Bool SCIPprimalUpperboundIsSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob            /**< original problem data */
   );

/** returns the primal ray thats proves unboundedness */
SCIP_SOL* SCIPprimalGetRay(
   SCIP_PRIMAL*          primal              /**< primal data */
   );

/** update the primal ray thats proves unboundedness */
SCIP_RETCODE SCIPprimalUpdateRay(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_SOL*             primalray,          /**< the new primal ray */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** adds primal solution to solution storage by copying it */
SCIP_RETCODE SCIPprimalAddSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds primal solution to solution storage, frees the solution afterwards */
SCIP_RETCODE SCIPprimalAddSolFree(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds primal solution to solution candidate storage of original problem space */
SCIP_RETCODE SCIPprimalAddOrigSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< original problem data */
   SCIP_SOL*             sol,                /**< primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds primal solution to solution candidate storage of original problem space, frees the solution afterwards */
SCIP_RETCODE SCIPprimalAddOrigSolFree(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< original problem data */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** adds current LP/pseudo solution to solution storage */
SCIP_RETCODE SCIPprimalAddCurrentSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution (or NULL if it's from the tree) */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage by copying it */
SCIP_RETCODE SCIPprimalTrySol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   );

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
SCIP_RETCODE SCIPprimalTrySolFree(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   );

/** checks current LP/pseudo solution; if feasible, adds it to storage */
SCIP_RETCODE SCIPprimalTryCurrentSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution (or NULL if it's from the tree) */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   );

/** inserts solution into the global array of all existing primal solutions */
SCIP_RETCODE SCIPprimalSolCreated(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** removes solution from the global array of all existing primal solutions */
void SCIPprimalSolFreed(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** updates all existing primal solutions after a change in a variable's objective value */
void SCIPprimalUpdateVarObj(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldobj,             /**< old objective value */
   SCIP_Real             newobj              /**< new objective value */
   );

/** retransforms all existing solutions to original problem space
 *
 * @note as a side effect, the objective value of the solutions can change (numerical errors)
 * so we update the objective cutoff value and upper bound accordingly
 */
SCIP_RETCODE SCIPprimalRetransformSolutions(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** tries to transform original solution to the transformed problem space */
SCIP_RETCODE SCIPprimalTransformSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sol,                /**< primal solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_Real*            solvals,            /**< array for internal use to store solution values, or NULL;
                                              *   if the method is called multiple times in a row, an array with size >=
                                              *   number of active variables should be given for performance reasons */
   SCIP_Bool*            solvalset,          /**< array for internal use to store which solution values were set, or NULL;
                                              *   if the method is called multiple times in a row, an array with size >=
                                              *   number of active variables should be given for performance reasons */
   int                   solvalssize,        /**< size of solvals and solvalset arrays, should be >= number of active
                                              *   variables */
   SCIP_Bool*            added               /**< pointer to store whether the solution was added */
   );

/** is the updating of violations enabled for this problem? */
SCIP_Bool SCIPprimalUpdateViolations(
   SCIP_PRIMAL*          primal              /**< problem data */
   );

/** set whether the updating of violations is turned on */
void SCIPprimalSetUpdateViolations(
   SCIP_PRIMAL*          primal,             /**< problem data */
   SCIP_Bool             updateviolations    /**< TRUE to enable violation updates, FALSE otherwise */
   );

#ifdef __cplusplus
}
#endif

#endif
