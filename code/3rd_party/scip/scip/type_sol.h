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

/**@file   type_sol.h
 * @brief  type definitions for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SOL_H__
#define __SCIP_TYPE_SOL_H__

#ifdef __cplusplus
extern "C" {
#endif

/** origin of solution: where to retrieve uncached elements */
enum SCIP_SolOrigin
{
   SCIP_SOLORIGIN_ORIGINAL  = 0,        /**< solution describes original variables; non-cached elements are zero */
   SCIP_SOLORIGIN_ZERO      = 1,        /**< all non-cached elements in solution are equal to zero */
   SCIP_SOLORIGIN_LPSOL     = 2,        /**< all non-cached elements in solution are equal to current LP solution */
   SCIP_SOLORIGIN_NLPSOL    = 3,        /**< all non-cached elements in solution are equal to current NLP solution */
   SCIP_SOLORIGIN_RELAXSOL  = 4,        /**< all non-cached elements in solution are equal to current relaxation solution */
   SCIP_SOLORIGIN_PSEUDOSOL = 5,        /**< all non-cached elements in solution are equal to current pseudo solution */
   SCIP_SOLORIGIN_PARTIAL   = 6,        /**< solution describes original solution; all non-cached elements in solution
                                         *   are treated as being an arbitrary value in the variable's bounds
                                         */
   SCIP_SOLORIGIN_UNKNOWN   = 7         /**< all non-cached elements in solution are unknown; they have to be treated
                                         *   as being an arbitrary value in the variable's bounds
                                         */
};
typedef enum SCIP_SolOrigin SCIP_SOLORIGIN;

typedef struct SCIP_Sol SCIP_SOL;                 /**< primal CIP solution */

typedef struct SCIP_Viol SCIP_VIOL;               /**< maximum violations of problem constraints */

/** type of solution: heuristic or (LP) relaxation solution, or unspecified origin */
enum SCIP_SolType
{
   SCIP_SOLTYPE_UNKNOWN   = 0,          /**< type of solution unspecified (the default) */
   SCIP_SOLTYPE_HEUR      = 1,          /**< solution was found by a heuristic */
   SCIP_SOLTYPE_RELAX     = 2,          /**< solution was found by a relaxation */
   SCIP_SOLTYPE_LPRELAX   = 3,          /**< solution was found by the LP relaxation */
   SCIP_SOLTYPE_STRONGBRANCH = 4,       /**< solution was found during strong branching */
   SCIP_SOLTYPE_PSEUDO    = 5           /**< solution originates from a pseudo solution */
};
typedef enum SCIP_SolType SCIP_SOLTYPE;

#ifdef __cplusplus
}
#endif

#endif
