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

/**@file   presol_boundshift.h
 * @ingroup PRESOLVERS
 * @brief  presolver that converts integer variables with domain [a,b] to integer variables with domain [0,b-a]
 * @author Tobias Achterberg
 * @author Michael Winkler
 *
 * This presolver converts all integer variables with domain \f$[a,b]\f$ to integer variables with domain
 * \f$[0,b-a]\f$. This is done by creating a new integer variable \f$y\f$ which will be aggregated to the old variable
 * \f$x\f$ such that
 * \f[
 * x = y + a
 * \f]
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_BOUNDSHIFT_H__
#define __SCIP_PRESOL_BOUNDSHIFT_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the boundshift presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresolBoundshift(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
