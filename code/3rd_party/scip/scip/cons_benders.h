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

/**@file   cons_benders.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_BENDERS_H__
#define __SCIP_CONS_BENDERS_H__


#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_cons.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for Benders' decomposition and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrBenders(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Benders Constraints
 *
 * Two constraint handlers are implemented for the generation of Benders' decomposition cuts. When included in a
 * problem, the Benders' decomposition constraint handlers generate cuts during the enforcement of LP and relaxation
 * solutions. Additionally, Benders' decomposition cuts can be generated when checking the feasibility of solutions with
 * respect to the subproblem constraints.
 *
 * This constraint handler has an enforcement priority that is less than the integer constraint handler. This means that
 * only integer feasible solutions from the LP solver are enforced by this constraint handler. This is the traditional
 * behaviour of the branch-and-check approach to Benders' decomposition. Additionally, the check priority is set low,
 * such that this expensive constraint handler is only called as a final check on primal feasible solutions.
 *
 * This constraint handler in the standard constraint handler that should be added when using Benders' decomposition.
 * Additionally, there is a flag in SCIPincludeConshdlrBenders that permits the addition of the LP constraint handler,
 * cons_benderslp. The use of both cons_benders and cons_benderslp allows the user to perform a multiphase Benders'
 * decomposition algorithm.
 *
 * @{
 */

/** enforces Benders' constraints for given solution
 *
 *  This method is called from cons_benderslp and cons_benders. If the method is called from cons_benderslp, then the
 *  solutions are not guaranteed to be integer feasible. This is because the default priority is set greater than the
 *  integer constraint handler. If this method is called from cons_benders, then, because the default enforcement
 *  priority is set less than that of the integer constraint handler, then it can be assumed that the solutions are
 *  integer feasible.
 *
 *  The checkint flag indicates whether integer feasibility can be assumed. If it is not assumed, i.e. checkint ==
 *  FALSE, then only the convex relaxations of the subproblems are solved. If integer feasibility is assumed, i.e.
 *  checkint == TRUE, then the convex relaxations and the full CIP are solved to generate Benders' cuts and check
 *  solution feasibility.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPconsBendersEnforceSolution(
   SCIP*                 scip,               /**< the SCIP instance */
   SCIP_SOL*             sol,                /**< the primal solution to enforce, or NULL for the current LP/pseudo sol */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_RESULT*          result,             /**< the result of the enforcement */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should integrality be considered when checking the subproblems */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
