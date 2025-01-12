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

/**@file   reader_sto.h
 * @ingroup FILEREADERS
 * @brief  STO file reader - the stochastic information of an instance in SMPS format
 * @author Stephen J. Maher
 *
 * This is a reader for the stochastic information of a stochastic programming instance in SMPS format.
 * The three files that must be read are:
 * - .cor
 * - .tim
 * - .sto
 *
 * Alternatively, it is possible to create a .smps file with the relative path to the .cor, .tim and .sto files.
 * A file reader is available for the .smps file.
 *
 * Details regarding the SMPS file format can be found at:
 * Birge, J. R.; Dempster, M. A.; Gassmann, H. I.; Gunn, E.; King, A. J. & Wallace, S. W.
 * A standard input format for multiperiod stochastic linear programs
 * IIASA, Laxenburg, Austria, WP-87-118, 1987
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_STO_H__
#define __SCIP_READER_STO_H__

#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_prob.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the sto file reader into SCIP
 *
 *  @ingroup FileReaderIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeReaderSto(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup FILEREADERS
 *
 * @{
 */

/** reads the stochastic information for a stochastic program that is in SMPS format */
SCIP_EXPORT
SCIP_RETCODE SCIPreadSto(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   );

/** writes problem to file */
SCIP_EXPORT
SCIP_RETCODE SCIPwriteSto(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< scalar applied to objective function; external objective value is
                                              * extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real             objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**            vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                   nvars,              /**< number of active variables in the problem */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars,           /**< number of general integer variables */
   int                   nimplvars,          /**< number of implicit integer variables */
   int                   ncontvars,          /**< number of continuous variables */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   );

/** returns the total number of scenarios added to the problem */
SCIP_EXPORT
int SCIPstoGetNScenarios(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
