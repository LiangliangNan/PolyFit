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

/**@file   reader_sol.h
 * @ingroup FILEREADERS
 * @brief  file reader for primal solutions
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 *
 * This reader handles solutions in two formats:
 *
 * - <b>SCIP raw format</b>@n
 *   The format is as follows:@n@n
 *   line 1: "solution status: <status>"@n
 *   line 2: "objective value: <value>"@n
 *   line 3+i: \<variable name\> \<value\> (obj: \<objective coefficient of variable\>)
 *   @n@n
 *   Only nonzero values need to be listed.
 *   @par
 *   Example:
 *   @code
 *     solution status: optimal
 *     objective value: 1
 *     x1  1 (obj:1)
 *     x2  1 (obj:0)
 *   @endcode
 * - <b>XML format</b>@n
 *   This format is used by CPLEX, for example. For reading we require a section of @p
 *   \<variables\>. Each entry in this section consists of@n
 *   \<variable name="<name>" index="<number>" value="<value>"/>
 *   @par
 *   Example:
 *   @code
 *   <?xml version = "1.0" standalone="yes"?>
 *   <variables>
 *      <variable name="x1" index="1" value="1"/>
 *      <variable name="x2" index="2" value="1"/>
 *   </variables>
 *   </xml>
 *   @endcode
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_SOL_H__
#define __SCIP_READER_SOL_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the sol file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeReaderSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
