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
EXTERN
SCIP_DISPDATA* SCIPdispGetData(
   SCIP_DISP*            disp                /**< display column */
   );

/** sets user data of display column; user has to free old data in advance! */
EXTERN
void SCIPdispSetData(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_DISPDATA*        dispdata            /**< new display column user data */
   );

/** gets name of display column */
EXTERN
const char* SCIPdispGetName(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets description of display column */
EXTERN
const char* SCIPdispGetDesc(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets head line of display column */
EXTERN
const char* SCIPdispGetHeader(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets width of display column */
EXTERN
int SCIPdispGetWidth(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets priority of display column */
EXTERN
int SCIPdispGetPriority(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets position of display column */
EXTERN
int SCIPdispGetPosition(
   SCIP_DISP*            disp                /**< display column */
   );

/** gets status of display column */
EXTERN
SCIP_DISPSTATUS SCIPdispGetStatus(
   SCIP_DISP*            disp                /**< display column */
   );

/** is display column initialized? */
EXTERN
SCIP_Bool SCIPdispIsInitialized(
   SCIP_DISP*            disp                /**< display column */
   );

/** displays a long integer in decimal form fitting in a given width */
EXTERN
void SCIPdispLongint(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   SCIP_Longint          val,                /**< value to display */
   int                   width               /**< width to fit into */
   );

/** displays an integer in decimal form fitting in a given width */
EXTERN
void SCIPdispInt(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   int                   val,                /**< value to display */
   int                   width               /**< width to fit into */
   );

/** displays a time value fitting in a given width */
EXTERN
void SCIPdispTime(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   SCIP_Real             val,                /**< value in seconds to display */
   int                   width               /**< width to fit into */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
