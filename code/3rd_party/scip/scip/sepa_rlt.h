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

/**@file   sepa_rlt.h
 * @ingroup SEPARATORS
 * @brief  reformulation-linearization technique separator
 * @author Fabian Wegscheider
 * @author Ksenia Bestuzheva
 *
 *
 * This separator generates a collection of cuts constructed by the reformulation-linearization technique (RLT).
 * For an LP row L and a variable x in [lb,ub], L is multiplied either with (ub-x) or with (x-lb). All known terms that
 * appear in the product are replaced by their respective auxiliary variable and all unknown terms are replaced by a
 * suitable linear relaxation, e.g., McCormick. In general, the separator computes four different cuts for a row with
 * finite sides and a variable with finite bounds.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_RLT_H__
#define __SCIP_SEPA_RLT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the RLT separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaRlt(
   SCIP*                 scip                /**< SCIP data structure */
);

/**@addtogroup SEPARATORS
 *
 * @{
 */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
