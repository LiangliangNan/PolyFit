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

/**@file   scip_validation.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for validation
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_VALIDATION_H__
#define __SCIP_SCIP_VALIDATION_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicValidationMethods
 *
 * @{
 */

/** validate the result of the solve
 *
 *  the validation includes
 *
 *  - checking the feasibility of the incumbent solution in the original problem (using SCIPcheckSolOrig())
 *
 *  - checking if the objective bounds computed by SCIP agree with external primal and dual reference bounds.
 *
 *  All external reference bounds the original problem space and the original objective sense.
 *
 *  For infeasible problems, +/-SCIPinfinity() should be passed as reference bounds depending on the objective sense
 *  of the original problem.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPvalidateSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             primalreference,    /**< external primal reference value for the problem, or SCIP_UNKNOWN */
   SCIP_Real             dualreference,      /**< external dual reference value for the problem, or SCIP_UNKNOWN */
   SCIP_Real             reftol,             /**< relative tolerance for acceptable violation of reference values */
   SCIP_Bool             quiet,              /**< TRUE if no status line should be printed */
   SCIP_Bool*            feasible,           /**< pointer to store if the best solution is feasible in the original problem,
                                               *  or NULL */
   SCIP_Bool*            primalboundcheck,   /**< pointer to store if the primal bound respects the given dual reference
                                               *  value, or NULL */
   SCIP_Bool*            dualboundcheck      /**< pointer to store if the dual bound respects the given primal reference
                                               *  value, or NULL */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
