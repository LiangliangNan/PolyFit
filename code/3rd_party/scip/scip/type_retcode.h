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

/**@file   type_retcode.h
 * @brief  type definitions for return codes for SCIP methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_RETCODE_H__
#define __SCIP_TYPE_RETCODE_H__

#ifdef __cplusplus
extern "C" {
#endif

/** return codes for SCIP methods: non-positive return codes are errors */
enum SCIP_Retcode
{
   SCIP_OKAY               =  +1,       /**< normal termination */
   SCIP_ERROR              =   0,       /**< unspecified error */
   SCIP_NOMEMORY           =  -1,       /**< insufficient memory error */
   SCIP_READERROR          =  -2,       /**< read error */
   SCIP_WRITEERROR         =  -3,       /**< write error */
   SCIP_NOFILE             =  -4,       /**< file not found error */
   SCIP_FILECREATEERROR    =  -5,       /**< cannot create file */
   SCIP_LPERROR            =  -6,       /**< error in LP solver */
   SCIP_NOPROBLEM          =  -7,       /**< no problem exists */
   SCIP_INVALIDCALL        =  -8,       /**< method cannot be called at this time in solution process */
   SCIP_INVALIDDATA        =  -9,       /**< error in input data */
   SCIP_INVALIDRESULT      = -10,       /**< method returned an invalid result code */
   SCIP_PLUGINNOTFOUND     = -11,       /**< a required plugin was not found */
   SCIP_PARAMETERUNKNOWN   = -12,       /**< the parameter with the given name was not found */
   SCIP_PARAMETERWRONGTYPE = -13,       /**< the parameter is not of the expected type */
   SCIP_PARAMETERWRONGVAL  = -14,       /**< the value is invalid for the given parameter */
   SCIP_KEYALREADYEXISTING = -15,       /**< the given key is already existing in table */
   SCIP_MAXDEPTHLEVEL      = -16,       /**< maximal branching depth level exceeded */
   SCIP_BRANCHERROR        = -17,       /**< no branching could be created */
   SCIP_NOTIMPLEMENTED     = -18        /**< function not implemented */
};
typedef enum SCIP_Retcode SCIP_RETCODE;           /**< return code for SCIP method */

#ifdef __cplusplus
}
#endif

#endif
