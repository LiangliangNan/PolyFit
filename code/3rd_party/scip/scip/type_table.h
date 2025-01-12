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

/**@file   type_table.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for displaying statistics tables
 * @author Tristan Gally
 *
 *  This file defines the interface for statistics tables implemented in C.
 *
 * - \ref TABLE "Instructions for implementing a statistics table"
 * - \ref TABLES "List of available statistics tables"
 * - \ref scip::ObjTable "C++ wrapper class
 */

/** @defgroup DEFPLUGINS_TABLE Default Tables
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default statistics tables of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_TABLE_H__
#define __SCIP_TYPE_TABLE_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Table SCIP_TABLE;        /**< statistics table data structure */
typedef struct SCIP_TableData SCIP_TABLEDATA; /**< statistics table specific data */


/**  copy method for statistics table plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - table           : the statistics table itself
 */
#define SCIP_DECL_TABLECOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_TABLE* table)

/** destructor of statistics table to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - table           : the statistics table itself
 */
#define SCIP_DECL_TABLEFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_TABLE* table)

/** initialization method of statistics table (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - table           : the statistics table itself
 */
#define SCIP_DECL_TABLEINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_TABLE* table)

/** deinitialization method of statistics table (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - table           : the statistics table itself
 */
#define SCIP_DECL_TABLEEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_TABLE* table)

/** solving process initialization method of statistics table (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The statistics table may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - table           : the statistics table itself
 */
#define SCIP_DECL_TABLEINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_TABLE* table)

/** solving process deinitialization method of statistics table (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The statistics table should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - table           : the display column itself
 */
#define SCIP_DECL_TABLEEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_TABLE* table)

/** output method of statistics table to output file stream 'file'
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - table           : the statistics table itself
 *  - file            : file stream for output
 */
#define SCIP_DECL_TABLEOUTPUT(x) SCIP_RETCODE x (SCIP* scip, SCIP_TABLE* table, FILE* file)

#ifdef __cplusplus
}
#endif

#endif
