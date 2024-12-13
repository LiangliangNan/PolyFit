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

/**@file   presol_gateextraction.h
 * @ingroup PRESOLVERS
 * @brief  gateextraction presolver
 * @author Michael Winkler
 */

/* This presolver tries to extract gate-constraints meaning and-constraints and set-partitioning constraints (and could
 * be expanded to find xor-constraints too). This is done by detecting linearizations or systems of inequalities which
 * form an and-constraint or a set-partitioning constraint. An example:
 *
 * we have a logicor constraint of the form:                x + y + z >= 1
 *
 * and we also have the following set-packing constraints: (x + y <= 1 and x + z <= 1) <=> (~x + ~y >= 1 and ~x + ~z >= 1)
 *
 * - these three constraints form an and-constraint:        x = ~y * ~z (x = AND(~y,~z))
 *
 * if an additional set-packing constraint exists:          y + z <= 1
 *
 * - these four constraints form a set-partitioning cons.:  x + y + z = 1
 *
 * some information can be found:
 *
 *  http://www.cs.ubc.ca/~hutter/earg/papers07/cnf-structure.pdf
 *  http://www.cadence.com/cn/cadence/cadence_labs/Documents/niklas_SAT_2005_Effective.pdf
 *
 * We also do some check for logicor and set-packing/-partitioning constraint with the same variables to upgrade these
 * both constraints into one. For example:
 *
 *  x + y + z >= 1 and x + y + z <= 1 form x + y + z = 1
 *
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_GATEEXTRACTION_H__
#define __SCIP_PRESOL_GATEEXTRACTION_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the gateextraction presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresolGateextraction(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
