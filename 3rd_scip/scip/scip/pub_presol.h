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

/**@file   pub_presol.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PRESOL_H__
#define __SCIP_PUB_PRESOL_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_presol.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicPresolverMethods
 *
 * @{
 */

/** compares two presolvers w. r. to their priority */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPpresolComp);

/** comparison method for sorting presolvers w.r.t. to their name */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPpresolCompName);

/** gets user data of presolver */
EXTERN
SCIP_PRESOLDATA* SCIPpresolGetData(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** sets user data of presolver; user has to free old data in advance! */
EXTERN
void SCIPpresolSetData(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_PRESOLDATA*      presoldata          /**< new presolver user data */
   );

/** gets name of presolver */
EXTERN
const char* SCIPpresolGetName(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets description of presolver */
EXTERN
const char* SCIPpresolGetDesc(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets priority of presolver */
EXTERN
int SCIPpresolGetPriority(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets round limit of presolver */
EXTERN
int SCIPpresolGetMaxrounds(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets the timing mask of the presolver */
EXTERN
SCIP_PRESOLTIMING SCIPpresolGetTiming(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** sets the timing mask of the presolver */
EXTERN
void SCIPpresolSetTiming(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_PRESOLTIMING     timing              /**< timing mask of the presolver */
   );

/** is presolver initialized? */
EXTERN
SCIP_Bool SCIPpresolIsInitialized(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets time in seconds used in this presolver for setting up for next stages */
EXTERN
SCIP_Real SCIPpresolGetSetupTime(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets time in seconds used in this presolver */
EXTERN
SCIP_Real SCIPpresolGetTime(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of variables fixed in presolver */
EXTERN
int SCIPpresolGetNFixedVars(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of variables aggregated in presolver */
EXTERN
int SCIPpresolGetNAggrVars(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of variable types changed in presolver */
EXTERN
int SCIPpresolGetNChgVarTypes(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of bounds changed in presolver */
EXTERN
int SCIPpresolGetNChgBds(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of holes added to domains of variables in presolver */
EXTERN
int SCIPpresolGetNAddHoles(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of constraints deleted in presolver */
EXTERN
int SCIPpresolGetNDelConss(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of constraints added in presolver */
EXTERN
int SCIPpresolGetNAddConss(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of constraints upgraded in presolver */
EXTERN
int SCIPpresolGetNUpgdConss(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of coefficients changed in presolver */
EXTERN
int SCIPpresolGetNChgCoefs(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of constraint sides changed in presolver */
EXTERN
int SCIPpresolGetNChgSides(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/** gets number of times the presolver was called and tried to find reductions */
EXTERN
int SCIPpresolGetNCalls(
   SCIP_PRESOL*          presol              /**< presolver */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
