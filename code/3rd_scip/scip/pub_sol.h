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

/**@file   pub_sol.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for primal CIP solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_SOL_H__
#define __SCIP_PUB_SOL_H__


#include "scip/def.h"
#include "scip/type_sol.h"
#include "scip/type_heur.h"
#include "scip/type_relax.h"

#ifdef NDEBUG
#include "scip/struct_sol.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicSolutionMethods
 *
 * @{
 */


/** gets origin of solution */
SCIP_EXPORT
SCIP_SOLORIGIN SCIPsolGetOrigin(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns whether the given solution is defined on original variables */
SCIP_EXPORT
SCIP_Bool SCIPsolIsOriginal(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns whether the given solution is partial */
SCIP_EXPORT
SCIP_Bool SCIPsolIsPartial(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets objective value of primal CIP solution which lives in the original problem space */
SCIP_EXPORT
SCIP_Real SCIPsolGetOrigObj(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets clock time, when this solution was found */
SCIP_EXPORT
SCIP_Real SCIPsolGetTime(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets branch and bound run number, where this solution was found */
SCIP_EXPORT
int SCIPsolGetRunnum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets node number of the specific branch and bound run, where this solution was found */
SCIP_EXPORT
SCIP_Longint SCIPsolGetNodenum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets node's depth, where this solution was found */
SCIP_EXPORT
int SCIPsolGetDepth(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets information if solution was found by the LP, a primal heuristic, or a custom relaxator */
SCIP_EXPORT
SCIP_SOLTYPE SCIPsolGetType(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets heuristic that found this solution, or NULL if solution has type different than SCIP_SOLTYPE_HEUR */
SCIP_EXPORT
SCIP_HEUR* SCIPsolGetHeur(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets relaxation handler that found this solution, or NULL if solution has different type than SCIP_SOLTYPE_RELAX */
SCIP_EXPORT
SCIP_RELAX* SCIPsolGetRelax(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** informs the solution that it now belongs to the given primal heuristic. For convenience and backwards compatibility,
 *  the method accepts NULL as input for \p heur, in which case the solution type is set to SCIP_SOLTYPE_LPRELAX.
 *
 *  @note Relaxation handlers should use SCIPsolSetRelax() instead.
 */
SCIP_EXPORT
void SCIPsolSetHeur(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_HEUR*            heur                /**< primal heuristic that found the solution, or NULL for LP solutions */
   );

/** informs the solution that it now belongs to the given relaxation handler */
SCIP_EXPORT
void SCIPsolSetRelax(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RELAX*           relax               /**< relaxator that found the solution */
   );

/** informs the solution that it is an LP relaxation solution */
SCIP_EXPORT
void SCIPsolSetLPRelaxation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** informs the solution that it is a solution found during strong branching */
SCIP_EXPORT
void SCIPsolSetStrongbranching(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** informs the solution that it originates from a pseudo solution */
SCIP_EXPORT
void SCIPsolSetPseudo(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns unique index of given solution */
SCIP_EXPORT
int SCIPsolGetIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum absolute bound violation of solution */
SCIP_EXPORT
SCIP_Real SCIPsolGetAbsBoundViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum relative bound violation of solution */
SCIP_EXPORT
SCIP_Real SCIPsolGetRelBoundViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum absolute integrality violation of solution */
SCIP_EXPORT
SCIP_Real SCIPsolGetAbsIntegralityViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum absolute LP row violation of solution */
SCIP_EXPORT
SCIP_Real SCIPsolGetAbsLPRowViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum relative LP row violation of solution */
SCIP_EXPORT
SCIP_Real SCIPsolGetRelLPRowViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum absolute constraint violation of solution */
SCIP_EXPORT
SCIP_Real SCIPsolGetAbsConsViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum relative constraint violation of solution */
SCIP_EXPORT
SCIP_Real SCIPsolGetRelConsViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsolGetOrigin(sol)           ((sol)->solorigin)
#define SCIPsolIsOriginal(sol)          ((sol)->solorigin == SCIP_SOLORIGIN_ORIGINAL || (sol)->solorigin == SCIP_SOLORIGIN_PARTIAL)
#define SCIPsolGetOrigObj(sol)          (sol)->obj
#define SCIPsolGetTime(sol)             (sol)->time
#define SCIPsolGetNodenum(sol)          (sol)->nodenum
#define SCIPsolGetRunnum(sol)           (sol)->runnum
#define SCIPsolGetDepth(sol)            (sol)->depth
#define SCIPsolGetHeur(sol)             ((sol)->type == SCIP_SOLTYPE_HEUR ? (sol)->creator.heur : NULL)
#define SCIPsolGetRelax(sol)            ((sol)->type == SCIP_SOLTYPE_RELAX ? (sol)->creator.relax : NULL)
#define SCIPsolGetIndex(sol)            (sol)->index
#define SCIPsolGetType(sol)             (sol)->type
#define SCIPsolSetLPRelaxation(sol)     ((sol)->type = SCIP_SOLTYPE_LPRELAX)
#define SCIPsolSetStrongbranching(sol)  ((sol)->type = SCIP_SOLTYPE_STRONGBRANCH)
#define SCIPsolSetPseudo(sol)           ((sol)->type = SCIP_SOLTYPE_PSEUDO)
#endif

/** @} */

#ifdef __cplusplus
}
#endif

#endif
