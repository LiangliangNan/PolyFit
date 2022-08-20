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

/**@file   pub_paramset.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for handling parameter settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PARAMSET_H__
#define __SCIP_PUB_PARAMSET_H__


#include "scip/def.h"
#include "scip/type_paramset.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns type of parameter */
EXTERN
SCIP_PARAMTYPE SCIPparamGetType(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns name of parameter */
EXTERN
const char* SCIPparamGetName(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns description of parameter */
EXTERN
const char* SCIPparamGetDesc(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns locally defined parameter specific data */
EXTERN
SCIP_PARAMDATA* SCIPparamGetData(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns whether parameter is advanced */
EXTERN
SCIP_Bool SCIPparamIsAdvanced(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns whether parameter is fixed */
EXTERN
SCIP_Bool SCIPparamIsFixed(
   SCIP_PARAM*           param               /**< parameter */
   );

/** sets fixing status of given parameter */
EXTERN
void SCIPparamSetFixed(
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Bool             fixed               /**< new fixing status of the parameter */
   );

/** returns value of SCIP_Bool parameter */
EXTERN
SCIP_Bool SCIPparamGetBool(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of SCIP_Bool parameter */
EXTERN
SCIP_Bool SCIPparamGetBoolDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of int parameter */
EXTERN
int SCIPparamGetInt(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns minimal value of int parameter */
EXTERN
int SCIPparamGetIntMin(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns maximal value of int parameter */
EXTERN
int SCIPparamGetIntMax(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of int parameter */
EXTERN
int SCIPparamGetIntDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of SCIP_Longint parameter */
EXTERN
SCIP_Longint SCIPparamGetLongint(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns minimal value of longint parameter */
EXTERN
SCIP_Longint SCIPparamGetLongintMin(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns maximal value of longint parameter */
EXTERN
SCIP_Longint SCIPparamGetLongintMax(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of SCIP_Longint parameter */
EXTERN
SCIP_Longint SCIPparamGetLongintDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of SCIP_Real parameter */
EXTERN
SCIP_Real SCIPparamGetReal(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns minimal value of real parameter */
EXTERN
SCIP_Real SCIPparamGetRealMin(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns maximal value of real parameter */
EXTERN
SCIP_Real SCIPparamGetRealMax(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of SCIP_Real parameter */
EXTERN
SCIP_Real SCIPparamGetRealDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of char parameter */
EXTERN
char SCIPparamGetChar(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns allowed values of char parameter, or NULL if everything is allowed */
EXTERN
char* SCIPparamGetCharAllowedValues(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of char parameter */
EXTERN
char SCIPparamGetCharDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns value of string parameter */
EXTERN
char* SCIPparamGetString(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns default value of String parameter */
EXTERN
char* SCIPparamGetStringDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

/** returns whether the parameter is on its default setting */
EXTERN
SCIP_Bool SCIPparamIsDefault(
   SCIP_PARAM*           param               /**< parameter */
   );

#ifdef __cplusplus
}
#endif

#endif
