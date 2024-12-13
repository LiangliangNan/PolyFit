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

/**@file   pub_dcmp.h
 * @ingroup DecompMethods
 * @brief  public methods for decompositions
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SCIP_PUB_DECOMP_H_
#define SCIP_PUB_DECOMP_H_


#include "blockmemshell/memory.h"
#include "scip/type_cons.h"
#include "scip/type_dcmp.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup DecompMethods
 *
 * @{
 */

/** creates a decomposition */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompCreate(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   nblocks,            /**< the number of blocks (without the linking block) */
   SCIP_Bool             original,           /**< is this a decomposition in the original (TRUE) or transformed space? */
   SCIP_Bool             benderslabels       /**< should the variables be labeled for the application of Benders' decomposition */
   );

/** frees a decomposition */
SCIP_EXPORT
void SCIPdecompFree(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns TRUE if decomposition is in the original space */
SCIP_EXPORT
SCIP_Bool SCIPdecompIsOriginal(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** sets the parameter that indicates whether the variables must be labeled for the application of Benders'
 * decomposition
 */
SCIP_EXPORT
void SCIPdecompSetUseBendersLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             benderslabels       /**< whether Benders' variable labels should be used */
   );

/** returns TRUE if the variables must be labeled for the application of Benders' decomposition */
SCIP_EXPORT
SCIP_Bool SCIPdecompUseBendersLabels(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of blocks of this decomposition */
SCIP_EXPORT
int SCIPdecompGetNBlocks(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets area score of this decomposition */
SCIP_EXPORT
SCIP_Real SCIPdecompGetAreaScore(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets modularity of this decomposition */
SCIP_EXPORT
SCIP_Real SCIPdecompGetModularity(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets variable size for each block, sorted by increasing block label
 *
 * To get all variable sizes, set nlabels to SCIPdecompGetNBlocks() + 1.
 * The first entry corresponds to the number of border variables.
 *
 * @note Ensure that SCIPcomputeDecompStats() has been called before.
 *       If the decomposition was read from a file, this was done automatically.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompGetVarsSize(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   int*                  varssize,           /**< array to store variable sizes of blocks*/
   int                   nblocks             /**< length of variable sizes array */
   );

/** gets constraint size for each block, sorted by increasing block label
 *
 * To get all constraint sizes, set nlabels to SCIPdecompGetNBlocks() + 1.
 * The first entry corresponds to the number of border constraints.
 *
 * @note Ensure that SCIPcomputeDecompStats() has been called before.
 *       If the decomposition was read from a file, this was done automatically.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompGetConssSize(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   int*                  consssize,          /**< array to store constraint sizes of blocks*/
   int                   nblocks             /**< length of constraint sizes array */
   );

/** gets number of border variables of this decomposition
 *
 * @note Ensure that SCIPcomputeDecompStats() has been called before.
 *       If the decomposition was read from a file, this was done automatically.
 */
SCIP_EXPORT
int SCIPdecompGetNBorderVars(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of border constraints of this decomposition
 *
 * @note Ensure that SCIPcomputeDecompStats() has been called before.
 *       If the decomposition was read from a file, this was done automatically.
 */
SCIP_EXPORT
int SCIPdecompGetNBorderConss(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of edges in the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetNBlockGraphEdges(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of connected components in the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetNBlockGraphComponents(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets number of articulation points in the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetNBlockGraphArticulations(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets the maximum degree of the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetBlockGraphMaxDegree(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** gets the minimum degree of the block-decomposition graph of this decomposition */
SCIP_EXPORT
int SCIPdecompGetBlockGraphMinDegree(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   );

/** sets labels for an array of variables */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompSetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< array of labels, one per variable */
   int                   nvars               /**< length of variables array */
   );

/** queries labels for an array of variables */
SCIP_EXPORT
void SCIPdecompGetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< buffer to store labels, one per variable */
   int                   nvars               /**< length of variables array */
   );

/** sets labels for an array of constraints */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompSetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   );

/** queries labels for an array of constraints */
SCIP_EXPORT
void SCIPdecompGetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   );

/** clears the corresponding labeling (constraints, variables, or both) of this decomposition */
SCIP_EXPORT
SCIP_RETCODE SCIPdecompClear(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             clearvarlabels,     /**< should the variable labels be cleared? */
   SCIP_Bool             clearconslabels     /**< should the constraint labels be cleared? */
   );

/** prints decomposition statistics into string buffer */
SCIP_EXPORT
char* SCIPdecompPrintStats(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   char*                 strbuf              /**< string buffer storage */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
