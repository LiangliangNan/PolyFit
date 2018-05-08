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
EXTERN
SCIP_SOLORIGIN SCIPsolGetOrigin(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns whether the given solution is defined on original variables */
EXTERN
SCIP_Bool SCIPsolIsOriginal(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** returns whether the given solution is partial */
EXTERN
SCIP_Bool SCIPsolIsPartial(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets objective value of primal CIP solution which lives in the original problem space */
EXTERN
SCIP_Real SCIPsolGetOrigObj(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets clock time, when this solution was found */
EXTERN
SCIP_Real SCIPsolGetTime(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets branch and bound run number, where this solution was found */
EXTERN
int SCIPsolGetRunnum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets node number of the specific branch and bound run, where this solution was found */
EXTERN
SCIP_Longint SCIPsolGetNodenum(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets node's depth, where this solution was found */
EXTERN
int SCIPsolGetDepth(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
EXTERN
SCIP_HEUR* SCIPsolGetHeur(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** informs the solution that it now belongs to the given primal heuristic */
EXTERN
void SCIPsolSetHeur(
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   );

/** returns unique index of given solution */
EXTERN
int SCIPsolGetIndex(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum absolute bound violation of solution */
EXTERN
SCIP_Real SCIPsolGetAbsBoundViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum relative bound violation of solution */
EXTERN
SCIP_Real SCIPsolGetRelBoundViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum absolute integrality violation of solution */
EXTERN
SCIP_Real SCIPsolGetAbsIntegralityViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum absolute LP row violation of solution */
EXTERN
SCIP_Real SCIPsolGetAbsLPRowViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum relative LP row violation of solution */
EXTERN
SCIP_Real SCIPsolGetRelLPRowViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum absolute constraint violation of solution */
EXTERN
SCIP_Real SCIPsolGetAbsConsViolation(
   SCIP_SOL*             sol                 /**< primal CIP solution */
   );

/** get maximum relative constraint violation of solution */
EXTERN
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
#define SCIPsolGetHeur(sol)             (sol)->heur
#define SCIPsolGetIndex(sol)            (sol)->index
#define SCIPsolSetHeur(sol,newheur)     ((sol)->heur = (newheur))
#endif

/* @} */

#ifdef __cplusplus
}
#endif

#endif
