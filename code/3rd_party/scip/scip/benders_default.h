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

/**@file   benders_default.h
 * @ingroup BENDERSDECOMPOSITION
 * @brief  default Benders' decomposition plugin
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERS_DEFAULT_H__
#define __SCIP_BENDERS_DEFAULT_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the default Benders' decomposition and includes it in SCIP
 *
 *  @ingroup BendersIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBendersDefault(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**@addtogroup BENDERS
 *
 * @{
 */

/** Creates a default Benders' decomposition algorithm and activates it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateBendersDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
