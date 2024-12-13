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

/**@file   sol.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SOL_H__
#define __SCIP_SOL_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_heur.h"
#include "scip/pub_sol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates primal CIP solution, initialized to zero */
SCIP_RETCODE SCIPsolCreate(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution in original problem space, initialized to the offset in the original problem */
SCIP_RETCODE SCIPsolCreateOriginal(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates a copy of a primal CIP solution */
SCIP_RETCODE SCIPsolCopy(
   SCIP_SOL**            sol,                /**< pointer to store the copy of the primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   );

/** transformes given original solution to the transformed space; a corresponding transformed solution has to be given
 *  which is copied into the existing solution and freed afterwards
 */
SCIP_RETCODE SCIPsolTransform(
   SCIP_SOL*             sol,                /**< primal CIP solution to change, living in original space */
   SCIP_SOL**            transsol,           /**< pointer to corresponding transformed primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal              /**< primal data */
   );

/** adjusts solution values of implicit integer variables in handed solution. Solution objective value is not
 *  deteriorated by this method.
 */
SCIP_RETCODE SCIPsolAdjustImplicitSolVals(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< either original or transformed problem, depending on sol origin */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             uselprows           /**< should LP row information be considered for none-objective variables */
   );

/** creates primal CIP solution, initialized to the current LP solution */
SCIP_RETCODE SCIPsolCreateLPSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current NLP solution */
SCIP_RETCODE SCIPsolCreateNLPSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current relaxation solution */
SCIP_RETCODE SCIPsolCreateRelaxSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current pseudo solution */
SCIP_RETCODE SCIPsolCreatePseudoSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to the current solution */
SCIP_RETCODE SCIPsolCreateCurrentSol(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates partial primal CIP solution, initialized to unknown values */
SCIP_RETCODE SCIPsolCreatePartial(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** creates primal CIP solution, initialized to unknown values */
SCIP_RETCODE SCIPsolCreateUnknown(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** frees primal CIP solution */
SCIP_RETCODE SCIPsolFree(
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PRIMAL*          primal              /**< primal data */
   );

/** copies current LP solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkLPSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** copies current NLP solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkNLPSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp                 /**< current NLP data */
   );

/** copies current relaxation solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkRelaxSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   );

/** copies current pseudo solution into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkPseudoSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** copies current solution (LP or pseudo solution) into CIP solution by linking */
SCIP_RETCODE SCIPsolLinkCurrentSol(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** clears primal CIP solution */
SCIP_RETCODE SCIPsolClear(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** declares all entries in the primal CIP solution to be unknown */
SCIP_RETCODE SCIPsolSetUnknown(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** stores solution values of variables in solution's own array */
SCIP_RETCODE SCIPsolUnlink(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< transformed problem data */
   );

/** sets value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolSetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Real             val                 /**< solution value of variable */
   );

/** increases value of variable in primal CIP solution */
SCIP_RETCODE SCIPsolIncVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Real             incval              /**< increment for solution value of variable */
   );

/** returns value of variable in primal CIP solution */
SCIP_Real SCIPsolGetVal(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   );

/** returns value of variable in primal ray represented by primal CIP solution */
SCIP_Real SCIPsolGetRayVal(
   SCIP_SOL*             sol,                /**< primal CIP solution, representing a primal ray */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var                 /**< variable to get value for */
   );


/** gets objective value of primal CIP solution in transformed problem */
SCIP_Real SCIPsolGetObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob            /**< original problem data */
   );

/** updates primal solutions after a change in a variable's objective value */
void SCIPsolUpdateVarObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldobj,             /**< old objective value */
   SCIP_Real             newobj              /**< new objective value */
   );

/* mark the given solution as partial solution */
SCIP_RETCODE SCIPsolMarkPartial(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars               /**< number of problem variables */
   );

/** checks primal CIP solution for feasibility
 *
 *  @note The difference between SCIPsolCheck() and SCIPcheckSolOrig() is that modifiable constraints are handled
 *        differently. There might be some variables which do not have an original counter part (e.g. in
 *        branch-and-price). Therefore, modifiable constraints can not be double-checked in the original space.
 */
SCIP_RETCODE SCIPsolCheck(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            feasible            /**< stores whether solution is feasible */
   );

/** try to round given solution */
SCIP_RETCODE SCIPsolRound(
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool*            success             /**< pointer to store whether rounding was successful */
   );

/** updates the solution value sums in variables by adding the value in the given solution */
void SCIPsolUpdateVarsum(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_Real             weight              /**< weight of solution in weighted average */
   );

/** retransforms solution to original problem space */
SCIP_RETCODE SCIPsolRetransform(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_Bool*            hasinfval           /**< pointer to store whether the solution has infinite values */
   );

/** recomputes the objective value of an original solution, e.g., when transferring solutions
 *  from the solution pool (objective coefficients might have changed in the meantime)
 */
void SCIPsolRecomputeObj(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob            /**< original problem */
   );


/** returns whether the given solutions in transformed space are equal */
SCIP_Bool SCIPsolsAreEqual(
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2,               /**< second primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob           /**< transformed problem after presolve */
   );

/** outputs non-zero elements of solution to file stream */
SCIP_RETCODE SCIPsolPrint(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             mipstart,           /**< should only discrete variables be printed? */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );

/** outputs non-zero elements of solution representing a ray to file stream */
SCIP_RETCODE SCIPsolPrintRay(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data (original or transformed) */
   SCIP_PROB*            transprob,          /**< transformed problem data or NULL (to display priced variables) */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   );


/** reset violations of a solution */
void SCIPsolResetViolations(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** update integrality violation of a solution */
void SCIPsolUpdateIntegralityViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolintegrality  /**< absolute violation of integrality */
   );

/** update bound violation of a solution */
void SCIPsolUpdateBoundViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolbounds,      /**< absolute violation of bounds */
   SCIP_Real             relviolbounds       /**< relative violation of bounds */
   );

/** update LP row violation of a solution */
void SCIPsolUpdateLPRowViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviollprows,      /**< absolute violation of LP rows */
   SCIP_Real             relviollprows       /**< relative violation of LP rows */
   );

/** update constraint violation of a solution */
void SCIPsolUpdateConsViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviolcons,        /**< absolute violation of constraint */
   SCIP_Real             relviolcons         /**< relative violation of constraint */
   );

/** update violation of a constraint that is represented in the LP */
void SCIPsolUpdateLPConsViolation(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             absviol,            /**< absolute violation of constraint */
   SCIP_Real             relviol             /**< relative violation of constraint */
   );



/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** adds value to the objective value of a given original primal CIP solution */
void SCIPsolOrigAddObjval(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real             addval              /**< offset value to add */
   );

/** gets current position of solution in array of existing solutions of primal data */
int SCIPsolGetPrimalIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** sets current position of solution in array of existing solutions of primal data */
void SCIPsolSetPrimalIndex(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   primalindex         /**< new primal index of solution */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsolOrigAddObjval(sol, addval) ((sol)->obj += (addval))
#define SCIPsolGetPrimalIndex(sol)      ((sol)->primalindex)
#define SCIPsolSetPrimalIndex(sol,idx)  { (sol)->primalindex = idx; }

#endif

#ifdef __cplusplus
}
#endif

#endif
