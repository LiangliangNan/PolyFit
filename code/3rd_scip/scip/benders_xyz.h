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

/**@file   benders_xyz.h
 * @ingroup BENDERSDECOMPOSITION
 * @brief  xyz Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERS_XYZ_H__
#define __SCIP_BENDERS_XYZ_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the xyz Benders' decomposition and includes it in SCIP
 *
 *  @ingroup BendersIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBendersXyz(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup BENDERSS
 *
 * @{
 */

/** TODO: add public methods to this group for documentation purposes */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
