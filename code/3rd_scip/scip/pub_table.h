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

/**@file   pub_table.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for displaying statistic tables
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_TABLE_H__
#define __SCIP_PUB_TABLE_H__


#include "scip/def.h"
#include "scip/type_set.h"
#include "scip/type_table.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicTableMethods
 *
 * @{
 */

/** gets user data of statistics table */
SCIP_EXPORT
SCIP_TABLEDATA* SCIPtableGetData(
   SCIP_TABLE*           table               /**< statistics table */
   );

/** sets user data of statistics table; user has to free old data in advance! */
SCIP_EXPORT
void SCIPtableSetData(
   SCIP_TABLE*           table,              /**< statistics table */
   SCIP_TABLEDATA*       tabledata           /**< new statistics table user data */
   );

/** gets name of statistics table */
SCIP_EXPORT
const char* SCIPtableGetName(
   SCIP_TABLE*           table               /**< statistics table */
   );

/** gets description of statistics table */
SCIP_EXPORT
const char* SCIPtableGetDesc(
   SCIP_TABLE*           table               /**< statistics table */
   );

/** gets position of statistics table */
SCIP_EXPORT
int SCIPtableGetPosition(
   SCIP_TABLE*           table               /**< statistics table */
   );

/** gets earliest stage of statistics table */
SCIP_EXPORT
SCIP_STAGE SCIPtableGetEarliestStage(
   SCIP_TABLE*           table               /**< statistics table */
   );

/** is statistics table currently active? */
SCIP_EXPORT
SCIP_Bool SCIPtableIsActive(
   SCIP_TABLE*           table               /**< statistics table */
   );

/** is statistics table initialized? */
SCIP_EXPORT
SCIP_Bool SCIPtableIsInitialized(
   SCIP_TABLE*           table               /**< statistics table */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
