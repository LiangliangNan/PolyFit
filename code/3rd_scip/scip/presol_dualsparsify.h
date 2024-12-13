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

/**@file   presol_dualsparsify.h
 * @brief  cancel nonzeros of the constraint matrix based on the columns
 * @author Dieter Weninger
 * @author Leona Gottwald
 * @author Ambros Gleixner
 * @author Weikun Chen
 *
 * This presolver attempts to cancel non-zero entries of the constraint
 * matrix by adding scaled columns to other columns.
 *
 * In more detail, for two columns A_{j.} and A_{k.}, suppose for a given value s, we have
 *
 *                  | A_{j.} | - | A_{j.} - s*A_{k.} | > eta,
 *
 * where eta is an nonnegative integer. Then we introduce a new variable y := s*x_j + x_k
 * and aggregate the variable x_k = y - s*x_j. After aggregation, the column of the variable
 * x_j is A_{j.} + s*A_{j.} which is sparser than A_{j.}. In the case that x_k is no implied
 * free variable, we need to add a new constraint l_k <= y - weight*x_j <= u_k into the problem
 * to keep the bounds constraints of variable x_k.
 *
 * Further information can be found in
 * Chen et al. "Two-row and two-column mixed-integer presolve using hasing-based pairing methods".
 *
 * @todo add infrastructure to SCIP for handling aggregated binary variables
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_DUALSPARSIFY_H__
#define __SCIP_PRESOL_DUALSPARSIFY_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the dual sparsify presolver and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresolDualsparsify(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
