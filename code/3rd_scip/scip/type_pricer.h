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

/**@file   type_pricer.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_PRICER_H__
#define __SCIP_TYPE_PRICER_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Pricer SCIP_PRICER;           /**< variable pricer data */
typedef struct SCIP_PricerData SCIP_PRICERDATA;   /**< locally defined variable pricer data */


/** copy method for pricer plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 *  - valid           : was the copying process valid? 
 */
#define SCIP_DECL_PRICERCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRICER* pricer, SCIP_Bool* valid)

/** destructor of variable pricer to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define SCIP_DECL_PRICERFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRICER* pricer)

/** initialization method of variable pricer (called after problem was transformed and pricer is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define SCIP_DECL_PRICERINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRICER* pricer)

/** deinitialization method of variable pricer (called before transformed problem is freed and pricer is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define SCIP_DECL_PRICEREXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRICER* pricer)

/** solving process initialization method of variable pricer (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The variable pricer may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define SCIP_DECL_PRICERINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRICER* pricer)

/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The variable pricer should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define SCIP_DECL_PRICEREXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRICER* pricer)

/** reduced cost pricing method of variable pricer for feasible LPs
 *
 *  Searches for variables that can contribute to improve the current LP's solution value.
 *  In standard branch-and-price, these are variables with negative dual feasibility, that is negative
 *  reduced costs for non-negative variables, positive reduced costs for non-positive variables,
 *  and non-zero reduced costs for variables that can be negative and positive.
 *
 *  The method is called in the LP solving loop after an LP was proven to be feasible.
 *
 *  Whenever the pricer finds a variable with negative dual feasibility, it should call SCIPcreateVar()
 *  and SCIPaddPricedVar() to add the variable to the problem. Furthermore, it should call the appropriate
 *  methods of the constraint handlers to add the necessary variable entries to the constraints.
 *
 *  In the usual case that the pricer either adds a new variable or ensures that there are no further variables with
 *  negative dual feasibility, the result pointer should be set to SCIP_SUCCESS. Only if the pricer aborts pricing
 *  without creating a new variable, but there might exist additional variables with negative dual feasibility, the
 *  result pointer should be set to SCIP_DIDNOTRUN.  In this case, which sometimes is referred to as "early branching",
 *  the lp solution will not be used as a lower bound. The pricer can, however, store a valid lower bound in the
 *  lowerbound pointer.  If you use your own branching rule (e.g., to branch on constraints), make sure that it is able
 *  to branch on pseudo solutions. Otherwise, SCIP will use its default branching rules (which all branch on
 *  variables). This could disturb the pricing problem or branching might not even be possible, e.g., if all yet created
 *  variables have already been fixed.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 *  - lowerbound      : pointer to store a lower bound found by the pricer
 *  - stopearly       : should pricing be stopped, although new variables were added? (doing early branching)
 *  - result          : pointer to store the result of the pricer call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
 *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP solution is
 *                      optimal
 *
 */
#define SCIP_DECL_PRICERREDCOST(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRICER* pricer, SCIP_Real* lowerbound, SCIP_Bool* stopearly, SCIP_RESULT* result)

/** Farkas pricing method of variable pricer for infeasible LPs
 *
 *  Searches for variables that can contribute to the feasibility of the current LP.
 *  In standard branch-and-price, these are variables with positive Farkas values:
 *
 *  The LP was proven infeasible, so we have an infeasibility proof by the dual Farkas multipliers y.
 *  With the values of y, an implicit inequality  y^T A x >= y^T b  is associated, with b given
 *  by the sides of the LP rows and the sign of y:
 *   - if y_i is positive, b_i is the left hand side of the row,
 *   - if y_i is negative, b_i is the right hand side of the row.
 *
 *  y is chosen in a way, such that the valid inequality  y^T A x >= y^T b  is violated by all x,
 *  especially by the (for this inequality least infeasible solution) x' defined by 
 *     x'_i := ub_i, if y^T A_i >= 0
 *     x'_i := lb_i, if y^T A_i < 0.
 *  Pricing in this case means to add variables i with positive Farkas value, i.e. y^T A_i x'_i > 0.
 *
 *  The method is called in the LP solving loop after an LP was proven to be infeasible.
 *
 *  Whenever the pricer finds a variable with positive Farkas value, it should call SCIPcreateVar()
 *  and SCIPaddPricedVar() to add the variable to the problem. Furthermore, it should call the appropriate
 *  methods of the constraint handlers to add the necessary variable entries to the constraints.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 *  - result          : pointer to store the result of the pricer call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : at least one variable was found, which can contribute to the feasibility of the current LP,
 *                      or it is ensured that no such variable exists
 *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP is indeed infeasible
 */
#define SCIP_DECL_PRICERFARKAS(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRICER* pricer, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
