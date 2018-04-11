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

/**@file   pub_prop.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PROP_H__
#define __SCIP_PUB_PROP_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_prop.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicPropagatorMethods
 *
 * @{
 */

/** compares two propagators w. r. to their priority */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPpropComp);

/** compares two propagators w. r. to their presolving priority */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPpropCompPresol);

/** comparison method for sorting propagators w.r.t. to their name */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPpropCompName);

/** gets user data of propagator */
EXTERN
SCIP_PROPDATA* SCIPpropGetData(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets user data of propagator; user has to free old data in advance! */
EXTERN
void SCIPpropSetData(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROPDATA*        propdata            /**< new propagator user data */
   );

/** gets name of propagator */
EXTERN
const char* SCIPpropGetName(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets description of propagator */
EXTERN
const char* SCIPpropGetDesc(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets priority of propagator */
EXTERN
int SCIPpropGetPriority(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets presolving priority of propagator */
EXTERN
int SCIPpropGetPresolPriority(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets frequency of propagator */
EXTERN
int SCIPpropGetFreq(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used for setting up this propagator for new stages */
EXTERN
SCIP_Real SCIPpropGetSetupTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets frequency of propagator */
EXTERN
void SCIPpropSetFreq(
   SCIP_PROP*            prop,               /**< propagator */
   int                   freq                /**< new frequency of propagator */
   );

/** gets time in seconds used in this propagator */
EXTERN
SCIP_Real SCIPpropGetTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator during strong branching */
EXTERN
SCIP_Real SCIPpropGetStrongBranchPropTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator for resolve propagation */
EXTERN
SCIP_Real SCIPpropGetRespropTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator for presolving */
EXTERN
SCIP_Real SCIPpropGetPresolTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets the total number of times, the propagator was called */
EXTERN
SCIP_Longint SCIPpropGetNCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets the total number of times, the propagator was called for resolving a propagation */
EXTERN
SCIP_Longint SCIPpropGetNRespropCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets total number of times, this propagator detected a cutoff */
EXTERN
SCIP_Longint SCIPpropGetNCutoffs(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets total number of domain reductions found by this propagator */
EXTERN
SCIP_Longint SCIPpropGetNDomredsFound(
   SCIP_PROP*            prop                /**< propagator */
   );

/** should propagator be delayed, if other propagators found reductions? */
EXTERN
SCIP_Bool SCIPpropIsDelayed(
   SCIP_PROP*            prop                /**< propagator */
   );

/** was propagator delayed at the last call? */
EXTERN
SCIP_Bool SCIPpropWasDelayed(
   SCIP_PROP*            prop                /**< propagator */
   );

/** is propagator initialized? */
EXTERN
SCIP_Bool SCIPpropIsInitialized(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variables fixed during presolving of propagator */
EXTERN
int SCIPpropGetNFixedVars(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variables aggregated during presolving of propagator  */
EXTERN
int SCIPpropGetNAggrVars(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variable types changed during presolving of propagator  */
EXTERN
int SCIPpropGetNChgVarTypes(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of bounds changed during presolving of propagator  */
EXTERN
int SCIPpropGetNChgBds(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of holes added to domains of variables during presolving of propagator  */
EXTERN
int SCIPpropGetNAddHoles(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints deleted during presolving of propagator */
EXTERN
int SCIPpropGetNDelConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints added during presolving of propagator */
EXTERN
int SCIPpropGetNAddConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints upgraded during presolving of propagator  */
EXTERN
int SCIPpropGetNUpgdConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of coefficients changed during presolving of propagator */
EXTERN
int SCIPpropGetNChgCoefs(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraint sides changed during presolving of propagator */
EXTERN
int SCIPpropGetNChgSides(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of times the propagator was called in presolving and tried to find reductions */
EXTERN
int SCIPpropGetNPresolCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** returns the timing mask of the propagator */
EXTERN
SCIP_PROPTIMING SCIPpropGetTimingmask(
   SCIP_PROP*            prop                /**< propagator */
   );

/** does the propagator perform presolving? */
EXTERN
SCIP_Bool SCIPpropDoesPresolve(
   SCIP_PROP*            prop                /**< propagator */
   );

/** returns the timing mask of the presolving method of the propagator */
EXTERN
SCIP_PRESOLTIMING SCIPpropGetPresolTiming(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets the timing mask of the presolving method of the propagator */
EXTERN
void SCIPpropSetPresolTiming(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PRESOLTIMING     presoltiming        /** timing mask to be set */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
