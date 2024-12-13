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

/**@file   reader_dec.h
 * @ingroup FILEREADERS
 * @brief  file reader for decompositions in the constraint based dec-file format.
 * @author Gregor Hendel
 *
 *
 * This reader allows to read a file containing decompositions for constraints of the current original problem. The
 * standard line ending for this format is '.dec'. The content of the file should obey the following format
 *
 *     \\ place for comments and statistics
 *     NBLOCKS
 *     2
 *     BLOCK 0
 *     consA
 *     consB
 *     [...]
 *     BLOCK 1
 *     consC
 *     consD
 *     [...]
 *     MASTERCONSS
 *     linkingcons
 *
 * A block in a problem decomposition is a set of constraints that are independent from all other blocks after removing
 * the special blocks of linking constraints denoted as MASTERCONSS.
 *
 * Imagine the following example, which involves seven variables
 * and the five constraints from the file above. The asterisks (*) indicate that the variable affects the feasibility
 * of the constraint. In the special case of a linear optimization problem, the asterisks correspond to the
 * nonzero entries of the constraint matrix.
 *
 *                     x1  x2  x3  x4  x5  x6  x7
 *            consA     *       *                 \ BLOCK 0
 *            consB     *   *                     /
 *            consC                 *   *         \ BLOCK 1
 *            consD                     *   *     /
 *     linkingconss     *   *   *   *   *   *   * > MASTERCONSS
 *
 * The nonzero pattern has been chosen in a way that after the removal of the last constraint 'linkingcons', the remaining problem
 * consists of two independent parts, namely the blocks '0' and '1'.
 *
 * The corresponding variable labels are inferred from the constraint labels. A variable is assigned the label
 *
 * - of its unique block, if it only occurs in exactly 1 named block, and probably in MASTERCONSS.
 * - the special label of a linking variable if it occurs only in the master constraints or in 2 or even more named blocks.
 *
 * @note A trivial decomposition is to assign all constraints of a problem to the MASTERCONSS.
 *
 * @note The number of blocks is the number of named blocks: a trivial decomposition should have 0 blocks
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_DEC_H__
#define __SCIP_READER_DEC_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the decomposition file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeReaderDec(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
