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

/**@file   type_implics.h
 * @brief  type definitions for implications, variable bounds, and cliques
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_IMPLICS_H__
#define __SCIP_TYPE_IMPLICS_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_VBounds SCIP_VBOUNDS;         /**< variable bounds of a variable x in the form x <= c*y or x >= c*y */
typedef struct SCIP_Implics SCIP_IMPLICS;         /**< implications in the form x <= 0 or x >= 1 ==> y <= b or y >= b for x binary, NULL if x nonbinary */
typedef struct SCIP_Clique SCIP_CLIQUE;           /**< single clique, stating that at most one of the binary variables can be fixed
                                              *   to the corresponding value */
typedef struct SCIP_CliqueTable SCIP_CLIQUETABLE; /**< collection of cliques */
typedef struct SCIP_CliqueList SCIP_CLIQUELIST;   /**< list of cliques for a single variable */

#ifdef __cplusplus
}
#endif

#endif
