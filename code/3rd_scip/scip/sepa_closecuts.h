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

/**@file   sepa_closecuts.h
 * @ingroup SEPARATORS
 * @brief  closecuts meta separator
 * @author Marc Pfetsch
 *
 * This separator generates a convex combination of the current LP solution and either the best
 * primal feasible solution or an interior point of the LP relaxation. If the convex combination is
 * proper, the new point is closer to the convex hull of the feasible points. The separator then
 * calls all other separators to separate this point. The idea is that in this way possibly "deeper"
 * cuts are generated. Note, however, that the new point is not a basic solution, i.e., separators
 * relying basis information, e.g., Gomory cuts, will not work.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_CLOSECUTS_H__
#define __SCIP_SEPA_CLOSECUTS_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the closecuts separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaClosecuts(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup SEPARATORS
 *
 * @{
 */

/** sets point to be used as base point for computing the point to be separated
 *
 *  The point is only stored if separation of relative interior points is used. The solution is copied.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsetBasePointClosecuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< base point solution */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
