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

/**@file   prop_pseudoobj.h
 * @ingroup PROPAGATORS
 * @brief  Pseudo objective propagator
 * @author Tobias Achterberg
 * @author Stefan Heinz
 *
 * This propagator propagates the objective function using the cutoff bound and the pseudo objective value. The pseudo
 * objective value can be seen as minimum activity of the linear objective function. Using this, this propagator checks
 * if variables with non-zero objective coefficients can exceed the cutoff bound. If this is the case the corresponding
 * bound can be tightened.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_PSEUDOOBJ_H__
#define __SCIP_PROP_PSEUDOOBJ_H__

#include "scip/def.h"
#include "scip/type_prop.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the pseudo objective function propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropPseudoobj(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PROPAGATORS
  *
  * @{
  */

/** propagates the cutoff bound for the given variables */
SCIP_EXPORT
SCIP_RETCODE SCIPpropagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variables to propagate */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   SCIP_Bool*            tightened           /**< pointer to if the domain was tightened */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
