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

/**@file   expr_erf.h
 * @brief  handler for Gaussian error function expressions
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_ERF_H__
#define __SCIP_EXPR_ERF_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @{
 *
 * @name Gaussian error function expression
 *
 * This expression handler provides the Gaussian error function, that is
 *
 * \f[
 *   x \mapsto \frac{2}{\sqrt{\pi}}\int_0^x \exp(-t^2) dt.
 * \f]
 *
 * @attention The implementation of this expression handler is incomplete.
 * It is not usable for most use cases so far.
 * @{
 */

/** creates an erf expression
 *
 * @attention The implementation of `erf` expressions is incomplete.
 * They are not usable for most use cases so far.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprErf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< child expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** indicates whether expression is of erf-type */
SCIP_EXPORT
SCIP_Bool SCIPisExprErf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** @}
  * @}
  */

/** creates the handler for erf expressions and includes it into SCIP
 *
 * @attention The implementation of this expression handler is incomplete.
 * It is not usable for most use cases so far.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrErf(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_ERF_H__ */
