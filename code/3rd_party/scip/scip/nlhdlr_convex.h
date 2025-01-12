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

/**@file   nlhdlr_convex.h
 * @ingroup NLHDLRS
 * @brief  nonlinear handlers for convex and concave expressions, respectively
 * @author Benjamin Mueller
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_CONVEX_H__
#define __SCIP_NLHDLR_CONVEX_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes convex nonlinear handler in nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrConvex(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes concave nonlinear handler in nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrConcave(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLHDLRS
 * @{
 *
 * @name Convex and concave nonlinear handlers
 *
 * These nonlinear handler detect convex and concave subexpressions and provide specialized estimation functionality.
 *
 * @{
 */

/** checks whether a given expression is convex or concave w.r.t. the original variables
 *
 * This function uses the methods that are used in the detection algorithm of the convex nonlinear handler.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPhasExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRCURV         curv,               /**< curvature to check for */
   SCIP_Bool*            success,            /**< buffer to store whether expression has curvature curv (w.r.t. original variables) */
   SCIP_HASHMAP*         assumevarfixed      /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   );

/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_CONVEX_H__ */
