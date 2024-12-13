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

/**@file   type_paramset.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for handling parameter settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_PARAMSET_H__
#define __SCIP_TYPE_PARAMSET_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** possible parameter types */
enum SCIP_ParamType
{
   SCIP_PARAMTYPE_BOOL    = 0,               /**< bool values: TRUE or FALSE */
   SCIP_PARAMTYPE_INT     = 1,               /**< integer values */
   SCIP_PARAMTYPE_LONGINT = 2,               /**< long integer values */
   SCIP_PARAMTYPE_REAL    = 3,               /**< real values */
   SCIP_PARAMTYPE_CHAR    = 4,               /**< characters */
   SCIP_PARAMTYPE_STRING  = 5                /**< strings: arrays of characters */
};
typedef enum SCIP_ParamType SCIP_PARAMTYPE;

/** possible parameter settings - used to determine the behavior of different SCIP components, e.g., heuristics, separators, ... */
enum SCIP_ParamSetting
{
   SCIP_PARAMSETTING_DEFAULT     = 0,        /**< use default values */

   SCIP_PARAMSETTING_AGGRESSIVE  = 1,        /**< set to aggressive settings */
   SCIP_PARAMSETTING_FAST        = 2,        /**< set to fast settings */
   SCIP_PARAMSETTING_OFF         = 3         /**< turn off */
};
typedef enum SCIP_ParamSetting SCIP_PARAMSETTING;

/** possible parameter emphases - used to determine the general SCIP behavior */
enum SCIP_ParamEmphasis
{
   SCIP_PARAMEMPHASIS_DEFAULT     = 0,        /**< use default values */

   SCIP_PARAMEMPHASIS_CPSOLVER    = 1,        /**< get CP like search (e.g. no LP relaxation) */
   SCIP_PARAMEMPHASIS_EASYCIP     = 2,        /**< solve easy problems fast */
   SCIP_PARAMEMPHASIS_FEASIBILITY = 3,        /**< detect feasibility fast */
   SCIP_PARAMEMPHASIS_HARDLP      = 4,        /**< be capable to handle hard LPs */
   SCIP_PARAMEMPHASIS_OPTIMALITY  = 5,        /**< prove optimality fast */
   SCIP_PARAMEMPHASIS_COUNTER     = 6,        /**< get a feasible and "fast" counting process */
   SCIP_PARAMEMPHASIS_PHASEFEAS   = 7,        /**< feasibility phase settings during 3-phase solving approach */
   SCIP_PARAMEMPHASIS_PHASEIMPROVE= 8,        /**< improvement phase settings during 3-phase solving approach */
   SCIP_PARAMEMPHASIS_PHASEPROOF  = 9,        /**< proof phase settings during 3-phase solving approach */
   SCIP_PARAMEMPHASIS_NUMERICS    = 10,       /**< emphasis parameters for increased numerical safety */
   SCIP_PARAMEMPHASIS_BENCHMARK   = 11        /**< do not try to avoid running into memory limit */
};
typedef enum SCIP_ParamEmphasis SCIP_PARAMEMPHASIS;

typedef struct SCIP_Param SCIP_PARAM;             /**< single parameter */
typedef struct SCIP_ParamData SCIP_PARAMDATA;     /**< locally defined parameter specific data */
typedef struct SCIP_ParamSet SCIP_PARAMSET;       /**< set of parameters */


/** information method for changes in the parameter
 *
 *  Method is called if the parameter was changed through a SCIPparamsetSetXyz() call
 *  (which is called by SCIPsetXyzParam()).
 *  It will not be called, if the parameter was changed directly by changing the value
 *  in the memory location.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    param           : the changed parameter (already set to its new value)
 */
#define SCIP_DECL_PARAMCHGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_PARAM* param)

#ifdef __cplusplus
}
#endif

#endif
