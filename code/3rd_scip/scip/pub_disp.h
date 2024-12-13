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

/**@file   pub_disp.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for displaying runtime statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_DISP_H__
#define __SCIP_PUB_DISP_H__


#include <stdio.h>

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_disp.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicDisplayMethods
 *
 * @{
 */

/** gets user data of display column */
SCIP_EXPORT
SCIP_DISPDATA* SCIPdispGetData(
   SCIP_DISP*            disp                /**< display column */
   );

/** sets user data of display column; user has to free old data in advance! */
SCIP_EXPORT
void SCIPdispSetData(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_DISPDATA*        dispdata            /**< new display column user data */
   );

/** gets name of display column */
SCIP_EXPORT
const char* SCIPdispGetName(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets description of display column */
SCIP_EXPORT
const char* SCIPdispGetDesc(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets head line of display column */
SCIP_EXPORT
const char* SCIPdispGetHeader(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets width of display column */
SCIP_EXPORT
int SCIPdispGetWidth(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets priority of display column */
SCIP_EXPORT
int SCIPdispGetPriority(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets position of display column */
SCIP_EXPORT
int SCIPdispGetPosition(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets status of display column */
SCIP_EXPORT
SCIP_DISPSTATUS SCIPdispGetStatus(
   SCIP_DISP*            disp                /**< display column */
   );

/** is display column initialized? */
SCIP_EXPORT
SCIP_Bool SCIPdispIsInitialized(
   SCIP_DISP*            disp                /**< display column */
   );

/** displays a long integer in decimal form fitting in a given width */
SCIP_EXPORT
void SCIPdispLongint(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   SCIP_Longint          val,                /**< value to display */
   int                   width               /**< width to fit into */
   );

/** displays an integer in decimal form fitting in a given width */
SCIP_EXPORT
void SCIPdispInt(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   int                   val,                /**< value to display */
   int                   width               /**< width to fit into */
   );

/** displays a time value fitting in a given width */
SCIP_EXPORT
void SCIPdispTime(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   SCIP_Real             val,                /**< value in seconds to display */
   int                   width               /**< width to fit into */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
