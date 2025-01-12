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

/**@file   pub_conflict.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for conflict analysis handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CONFLICT_H__
#define __SCIP_PUB_CONFLICT_H__



#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_conflict.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicConflictMethods
 *
 * @{
 */

/** compares two conflict handlers w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrComp);

/** comparison method for sorting conflict handler w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrCompName);

/** gets user data of conflict handler */
SCIP_EXPORT
SCIP_CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** sets user data of conflict handler; user has to free old data in advance! */
SCIP_EXPORT
void SCIPconflicthdlrSetData(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   );

/** gets name of conflict handler */
SCIP_EXPORT
const char* SCIPconflicthdlrGetName(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets description of conflict handler */
SCIP_EXPORT
const char* SCIPconflicthdlrGetDesc(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets priority of conflict handler */
SCIP_EXPORT
int SCIPconflicthdlrGetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** is conflict handler initialized? */
SCIP_EXPORT
SCIP_Bool SCIPconflicthdlrIsInitialized(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets time in seconds used in this conflict handler for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPconflicthdlrGetSetupTime(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets time in seconds used in this conflict handler */
SCIP_EXPORT
SCIP_Real SCIPconflicthdlrGetTime(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
