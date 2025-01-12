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

/**@file   pub_sepa.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_SEPA_H__
#define __SCIP_PUB_SEPA_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_sepa.h"

#ifdef __cplusplus
extern "C" {
#endif

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@addtogroup PublicSeparatorMethods
 *
 * @{
 */


/** compares two separators w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsepaComp);

/** comparison method for sorting separators w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsepaCompName);

/** gets user data of separator */
SCIP_EXPORT
SCIP_SEPADATA* SCIPsepaGetData(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** sets user data of separator; user has to free old data in advance! */
SCIP_EXPORT
void SCIPsepaSetData(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata            /**< new separator user data */
   );

/** gets name of separator */
SCIP_EXPORT
const char* SCIPsepaGetName(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets description of separator */
SCIP_EXPORT
const char* SCIPsepaGetDesc(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets priority of separator */
SCIP_EXPORT
int SCIPsepaGetPriority(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets frequency of separator */
SCIP_EXPORT
int SCIPsepaGetFreq(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** sets frequency of separator */
SCIP_EXPORT
void SCIPsepaSetFreq(
   SCIP_SEPA*            sepa,               /**< separator */
   int                   freq                /**< new frequency of separator */
   );

/** get maximal bound distance at which the separator is called */
SCIP_EXPORT
SCIP_Real SCIPsepaGetMaxbounddist(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** does the separator use a secondary SCIP instance? */
SCIP_EXPORT
SCIP_Bool SCIPsepaUsesSubscip(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets time in seconds used in this separator for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPsepaGetSetupTime(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets time in seconds used in this separator */
SCIP_EXPORT
SCIP_Real SCIPsepaGetTime(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of times the separator was called */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCalls(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of times the separator was called at the root */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNRootCalls(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the number of times, the separator was called at the current node */
SCIP_EXPORT
int SCIPsepaGetNCallsAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of times, the separator detected a cutoff */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutoffs(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes added from the separator to the cut pool
 *  and to the sepastore directly */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes added from the separator to the sepastore;
 *  equal to the sum of added cuts directly and via the pool. */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsAdded(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the number of cutting planes from the separator added from the cut pool */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsAddedViaPool(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the number of cutting planes from the separator added directly to the sepastore */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsAddedDirect(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes from the separator applied to the LP */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsApplied(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes from the separator applied to the LP from the cutpool */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsAppliedViaPool(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes from the separator applied directly to the LP */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsAppliedDirect(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the number of cutting planes found by this separator at the current node */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNCutsFoundAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of additional constraints added by this separator */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNConssFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of domain reductions found by this separator */
SCIP_EXPORT
SCIP_Longint SCIPsepaGetNDomredsFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** should separator be delayed, if other separators found cuts? */
SCIP_EXPORT
SCIP_Bool SCIPsepaIsDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** was separation of the LP solution delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPsepaWasLPDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** was separation of the primal solution delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPsepaWasSolDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** is separator initialized? */
SCIP_EXPORT
SCIP_Bool SCIPsepaIsInitialized(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets whether separator is a parent separator */
SCIP_EXPORT
SCIP_Bool SCIPsepaIsParentsepa(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets parent separator (or NULL) */
SCIP_EXPORT
SCIP_SEPA* SCIPsepaGetParentsepa(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
