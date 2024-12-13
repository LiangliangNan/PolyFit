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
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPpropComp);

/** compares two propagators w. r. to their presolving priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPpropCompPresol);

/** comparison method for sorting propagators w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPpropCompName);

/** gets user data of propagator */
SCIP_EXPORT
SCIP_PROPDATA* SCIPpropGetData(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets user data of propagator; user has to free old data in advance! */
SCIP_EXPORT
void SCIPpropSetData(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROPDATA*        propdata            /**< new propagator user data */
   );

/** gets name of propagator */
SCIP_EXPORT
const char* SCIPpropGetName(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets description of propagator */
SCIP_EXPORT
const char* SCIPpropGetDesc(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets priority of propagator */
SCIP_EXPORT
int SCIPpropGetPriority(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets presolving priority of propagator */
SCIP_EXPORT
int SCIPpropGetPresolPriority(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets frequency of propagator */
SCIP_EXPORT
int SCIPpropGetFreq(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used for setting up this propagator for new stages */
SCIP_EXPORT
SCIP_Real SCIPpropGetSetupTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets frequency of propagator */
SCIP_EXPORT
void SCIPpropSetFreq(
   SCIP_PROP*            prop,               /**< propagator */
   int                   freq                /**< new frequency of propagator */
   );

/** gets time in seconds used in this propagator */
SCIP_EXPORT
SCIP_Real SCIPpropGetTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator during strong branching */
SCIP_EXPORT
SCIP_Real SCIPpropGetStrongBranchPropTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator for resolve propagation */
SCIP_EXPORT
SCIP_Real SCIPpropGetRespropTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets time in seconds used in this propagator for presolving */
SCIP_EXPORT
SCIP_Real SCIPpropGetPresolTime(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets the total number of times, the propagator was called */
SCIP_EXPORT
SCIP_Longint SCIPpropGetNCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets the total number of times, the propagator was called for resolving a propagation */
SCIP_EXPORT
SCIP_Longint SCIPpropGetNRespropCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets total number of times, this propagator detected a cutoff */
SCIP_EXPORT
SCIP_Longint SCIPpropGetNCutoffs(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets total number of domain reductions found by this propagator */
SCIP_EXPORT
SCIP_Longint SCIPpropGetNDomredsFound(
   SCIP_PROP*            prop                /**< propagator */
   );

/** should propagator be delayed, if other propagators found reductions? */
SCIP_EXPORT
SCIP_Bool SCIPpropIsDelayed(
   SCIP_PROP*            prop                /**< propagator */
   );

/** was propagator delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPpropWasDelayed(
   SCIP_PROP*            prop                /**< propagator */
   );

/** is propagator initialized? */
SCIP_EXPORT
SCIP_Bool SCIPpropIsInitialized(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variables fixed during presolving of propagator */
SCIP_EXPORT
int SCIPpropGetNFixedVars(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variables aggregated during presolving of propagator  */
SCIP_EXPORT
int SCIPpropGetNAggrVars(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of variable types changed during presolving of propagator  */
SCIP_EXPORT
int SCIPpropGetNChgVarTypes(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of bounds changed during presolving of propagator  */
SCIP_EXPORT
int SCIPpropGetNChgBds(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of holes added to domains of variables during presolving of propagator  */
SCIP_EXPORT
int SCIPpropGetNAddHoles(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints deleted during presolving of propagator */
SCIP_EXPORT
int SCIPpropGetNDelConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints added during presolving of propagator */
SCIP_EXPORT
int SCIPpropGetNAddConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraints upgraded during presolving of propagator  */
SCIP_EXPORT
int SCIPpropGetNUpgdConss(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of coefficients changed during presolving of propagator */
SCIP_EXPORT
int SCIPpropGetNChgCoefs(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of constraint sides changed during presolving of propagator */
SCIP_EXPORT
int SCIPpropGetNChgSides(
   SCIP_PROP*            prop                /**< propagator */
   );

/** gets number of times the propagator was called in presolving and tried to find reductions */
SCIP_EXPORT
int SCIPpropGetNPresolCalls(
   SCIP_PROP*            prop                /**< propagator */
   );

/** returns the timing mask of the propagator */
SCIP_EXPORT
SCIP_PROPTIMING SCIPpropGetTimingmask(
   SCIP_PROP*            prop                /**< propagator */
   );

/** does the propagator perform presolving? */
SCIP_EXPORT
SCIP_Bool SCIPpropDoesPresolve(
   SCIP_PROP*            prop                /**< propagator */
   );

/** returns the timing mask of the presolving method of the propagator */
SCIP_EXPORT
SCIP_PRESOLTIMING SCIPpropGetPresolTiming(
   SCIP_PROP*            prop                /**< propagator */
   );

/** sets the timing mask of the presolving method of the propagator */
SCIP_EXPORT
void SCIPpropSetPresolTiming(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PRESOLTIMING     presoltiming        /** timing mask to be set */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
