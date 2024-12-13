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

/**@file   event_globalbnd.h
 * @ingroup EVENTS
 * @brief  eventhdlr for storing all global bound changes
 * @author Leona Gottwald
 *
 * the bound changes are stored so that they can be shared with other threads
 * in a concurrent solve.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_GLOBALBND_H__
#define __SCIP_EVENT_GLOBALBND_H__

#include "scip/def.h"
#include "scip/type_event.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_syncstore.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates event handler for global bound event */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeEventHdlrGlobalbnd(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the global bound changes stored in the eventhandler */
SCIP_EXPORT
SCIP_BOUNDSTORE* SCIPeventGlobalbndGetBoundChanges(
   SCIP_EVENTHDLR*       eventhdlr           /**< the globalbound eventhandler */
   );

/** enables storing of bound changes */
SCIP_EXPORT
void SCIPeventGlobalbndEnableBoundStorage(
   SCIP_EVENTHDLR*       eventhdlr           /**< the globalbound eventhandler */
   );

/** disables storing of bound changes */
SCIP_EXPORT
void SCIPeventGlobalbndDisableBoundStorage(
   SCIP_EVENTHDLR*       eventhdlr           /**< the globalbound eventhandler */
   );

/** clears all bound changes stored in the eventhandler */
SCIP_EXPORT
void SCIPeventGlobalbndClearBoundChanges(
   SCIP_EVENTHDLR*       eventhdlr           /**< the globalbound eventhandler */
   );


#ifdef __cplusplus
}
#endif

#endif
