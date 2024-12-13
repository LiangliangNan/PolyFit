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

/**@file   struct_primal.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PRIMAL_H__
#define __SCIP_STRUCT_PRIMAL_H__


#include "scip/def.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"

#ifdef __cplusplus
extern "C" {
#endif

/** primal data and solution storage */
struct SCIP_Primal
{
   SCIP_Longint          nsolsfound;         /**< number of primal CIP solutions found up to now */
   SCIP_Longint          nlimsolsfound;      /**< number of primal CIP solutions respecting the objective limit found
                                              *   up to now */
   SCIP_Longint          nbestsolsfound;     /**< number of new best primal CIP solutions found up to now */
   SCIP_Longint          nlimbestsolsfound;  /**< number of new best primal CIP solutions respecting the objective limit
                                              *   found up to now */
   SCIP_Real             upperbound;         /**< upper (primal) bound of CIP: objective value of best solution or user bound */
   SCIP_Real             cutoffbound;        /**< upper bound for better primal solutions (if objective value is always
                                              *   integral, cutoffbound is equal to ceil(upperbound) - 1.0 (+eps) */
   SCIP_SOL**            sols;               /**< primal CIP solutions */
   SCIP_SOL**            partialsols;        /**< partial solutions */
   SCIP_SOL**            existingsols;       /**< all existing primal solutions (feasible, partial, and infeasible) */
   SCIP_SOL*             currentsol;         /**< internal solution for temporarily storing the current solution */
   SCIP_SOL*             primalray;          /**< solution representing the primal ray for (infeasible or) unbounded problems;
                                              *   warning: this does not have to be a feasible solution */
   int                   solssize;           /**< size of sols array */
   int                   partialsolssize;    /**< size of partialsols array */
   int                   nsols;              /**< number of primal CIP solutions stored in sols array */
   int                   npartialsols;       /**< number of partial solutions stored in partialsol array */
   int                   existingsolssize;   /**< size of existingsols array */
   int                   nexistingsols;      /**< number of primal CIP solutions stored in existingsols array */

   SCIP_Bool             updateviolations;   /**< marks whether the updating of violations is turned on */
};

#ifdef __cplusplus
}
#endif

#endif
