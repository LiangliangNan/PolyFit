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

/**@file   pub_cutsel.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for cut selectors
 * @author Mark Turner
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CUTSEL_H__
#define __SCIP_PUB_CUTSEL_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_cutsel.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCutSelectorMethods
 *
 * @{
 */

/** gets name of cut selector */
SCIP_EXPORT
const char* SCIPcutselGetName(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** gets user data of cut selector */
SCIP_EXPORT
SCIP_CUTSELDATA* SCIPcutselGetData(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** gets description of cut selector */
SCIP_EXPORT
const char* SCIPcutselGetDesc(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** gets priority of cut selector */
SCIP_EXPORT
int SCIPcutselGetPriority(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** sets user data of cut selector; user has to free old data in advance! */
SCIP_EXPORT
void SCIPcutselSetData(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_CUTSELDATA*      cutseldata          /**< new cut selector user data */
   );

/** is cut selector initialized? */
SCIP_EXPORT
SCIP_Bool SCIPcutselIsInitialized(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** gets time in seconds used in this cut selector for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPcutselGetSetupTime(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** gets time in seconds used in this cut selector */
SCIP_EXPORT
SCIP_Real SCIPcutselGetTime(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** get number of times the cutselector was called */
SCIP_Longint SCIPcutselGetNCalls(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** get number of times the cutselector was called at the root */
SCIP_Longint SCIPcutselGetNRootCalls(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** get total number of cuts that were selected at the root */
SCIP_Longint SCIPcutselGetNRootCuts(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** get total number of forced cuts that were selected at the root */
SCIP_Longint SCIPcutselGetNRootForcedCuts(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** get total number of root cuts that were filtered */
SCIP_Longint SCIPcutselGetNRootCutsFiltered(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** get total number of local cuts that were selected */
SCIP_Longint SCIPcutselGetNLocalCuts(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** get total number of forced local cuts that were selected */
SCIP_Longint SCIPcutselGetNLocalForcedCuts(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** get total number of local cuts that were filtered */
SCIP_Longint SCIPcutselGetNLocalCutsFiltered(
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   );

/** compares two cut selectors w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPcutselComp);

/** @} */

#ifdef __cplusplus
}
#endif

#endif
