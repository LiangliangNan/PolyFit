/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_paramset.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for handling parameter settings
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PARAMSET_H__
#define __SCIP_STRUCT_PARAMSET_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_paramset.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data for SCIP_Bool parameters */
struct SCIP_BoolParam
{
   SCIP_Bool*            valueptr;           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   SCIP_Bool             defaultvalue;       /**< default value of the parameter */
};
typedef struct SCIP_BoolParam SCIP_BOOLPARAM;

/** data for int parameters */
struct SCIP_IntParam
{
   int*                  valueptr;           /**< pointer to store the current parameter value, or NULL */
   int                   curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   int                   defaultvalue;       /**< default value of the parameter */
   int                   minvalue;           /**< minimum value for parameter */
   int                   maxvalue;           /**< maximum value for parameter */
};
typedef struct SCIP_IntParam SCIP_INTPARAM;

/** data for SCIP_Longint parameters */
struct SCIP_LongintParam
{
   SCIP_Longint          curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   SCIP_Longint          defaultvalue;       /**< default value of the parameter */
   SCIP_Longint          minvalue;           /**< minimum value for parameter */
   SCIP_Longint          maxvalue;           /**< maximum value for parameter */
   SCIP_Longint*         valueptr;           /**< pointer to store the current parameter value, or NULL */
};
typedef struct SCIP_LongintParam SCIP_LONGINTPARAM;

/** data for SCIP_Real parameters */
struct SCIP_RealParam
{
   SCIP_Real             curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   SCIP_Real             defaultvalue;       /**< default value of the parameter */
   SCIP_Real             minvalue;           /**< minimum value for parameter */
   SCIP_Real             maxvalue;           /**< maximum value for parameter */
   SCIP_Real*            valueptr;           /**< pointer to store the current parameter value, or NULL */
};
typedef struct SCIP_RealParam SCIP_REALPARAM;

/** data for char parameters */
struct SCIP_CharParam
{
   char*                 valueptr;           /**< pointer to store the current parameter value, or NULL */
   char*                 allowedvalues;      /**< array with possible parameter values, or NULL if not restricted */
   char                  curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   char                  defaultvalue;       /**< default value of the parameter */
};
typedef struct SCIP_CharParam SCIP_CHARPARAM;

/** data for char* parameters */
struct SCIP_StringParam
{
   char**                valueptr;           /**< pointer to store the current parameter value, or NULL */
   char*                 curvalue;           /**< stores the current parameter value if it is not stored in *valueptr */
   char*                 defaultvalue;       /**< default value of the parameter */
};
typedef struct SCIP_StringParam SCIP_STRINGPARAM;

/** single parameter */
struct SCIP_Param
{
   union
   {
      SCIP_BOOLPARAM     boolparam;          /**< data for SCIP_Bool parameters */
      SCIP_INTPARAM      intparam;           /**< data for int parameters */
      SCIP_LONGINTPARAM  longintparam;       /**< data for SCIP_Longint parameters */
      SCIP_REALPARAM     realparam;          /**< data for SCIP_Real parameters */
      SCIP_CHARPARAM     charparam;          /**< data for char parameters */
      SCIP_STRINGPARAM   stringparam;        /**< data for char* parameters */
   } data;
   char*                 name;               /**< name of the parameter */
   char*                 desc;               /**< description of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd));     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata;          /**< locally defined parameter specific data */
   unsigned int          isadvanced:1;       /**< is this parameter an advanced parameter? */
   unsigned int          isfixed:1;          /**< is this parameter fixed? */
   SCIP_PARAMTYPE        paramtype;          /**< type of this parameter */
};

/** set of parameters */
struct SCIP_ParamSet
{
   SCIP_HASHTABLE*       hashtable;          /**< hash table to store the parameters */
   SCIP_PARAM**          params;             /**< array with parameters */
   int                   nparams;            /**< number of parameters */
   int                   paramssize;         /**< size of params array */
};

#ifdef __cplusplus
}
#endif

#endif
