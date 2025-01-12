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

/**@file   scip_message.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for message handling
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_MESSAGE_H__
#define __SCIP_SCIP_MESSAGE_H__


#include "scip/def.h"
#include "scip/type_message.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup MessageOutputMethods
 *
 * @{
 */

/* if we have a C99 compiler */
#ifdef SCIP_HAVE_VARIADIC_MACROS

/** prints a debugging message if SCIP_DEBUG flag is set */
#ifdef SCIP_DEBUG
#define SCIPdebugMsg(scip, ...)         SCIPprintDebugMessage(scip, __FILE__, __LINE__, __VA_ARGS__)
#define SCIPdebugMsgPrint(scip, ...)    SCIPdebugMessagePrint(scip, __VA_ARGS__)
#else
#define SCIPdebugMsg(scip, ...)         while ( FALSE ) SCIPprintDebugMessage(scip, __FILE__, __LINE__, __VA_ARGS__)
#define SCIPdebugMsgPrint(scip, ...)    while ( FALSE ) SCIPdebugMessagePrint(scip, __VA_ARGS__)
#endif

#else
/* if we do not have a C99 compiler, use a workaround that prints a message, but not the file and linenumber */

/** prints a debugging message if SCIP_DEBUG flag is set */
#ifdef SCIP_DEBUG
#define SCIPdebugMsg                    printf("debug: "), SCIPdebugMessagePrint
#define SCIPdebugMsgPrint               SCIPdebugMessagePrint
#else
#define SCIPdebugMsg                    while ( FALSE ) SCIPdebugMessagePrint
#define SCIPdebugMsgPrint               while ( FALSE ) SCIPdebugMessagePrint
#endif

#endif


/** installs the given message handler, such that all messages are passed to this handler. A messages handler can be
 *  created via SCIPmessagehdlrCreate().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note The currently installed messages handler gets freed if this SCIP instance is its last user (w.r.t. capture/release).
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsetMessagehdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   );

/** returns the currently installed message handler
 *
 *  @return the currently installed message handler, or NULL if messages are currently suppressed
 */
SCIP_EXPORT
SCIP_MESSAGEHDLR* SCIPgetMessagehdlr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the log file name for the currently installed message handler */
SCIP_EXPORT
void SCIPsetMessagehdlrLogfile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of log file, or NULL (no log) */
   );

/** sets the currently installed message handler to be quiet (or not) */
SCIP_EXPORT
void SCIPsetMessagehdlrQuiet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   );

/** prints a warning message via the message handler */
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
SCIP_EXPORT
void SCIPwarningMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a debug message */
#ifdef __GNUC__
__attribute__((format(printf, 4, 5)))
#endif
SCIP_EXPORT
void SCIPprintDebugMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a debug message without precode */
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
SCIP_EXPORT
void SCIPdebugMessagePrint(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction or is a direct response to a user interactive command */
#ifdef __GNUC__
__attribute__((format(printf, 3, 4)))
#endif
SCIP_EXPORT
void SCIPdialogMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message */
#ifdef __GNUC__
__attribute__((format(printf, 3, 4)))
#endif
SCIP_EXPORT
void SCIPinfoMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message depending on the verbosity level */
#ifdef __GNUC__
__attribute__((format(printf, 4, 5)))
#endif
SCIP_EXPORT
void SCIPverbMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** returns the current message verbosity level
 *
 *  @return message verbosity level of SCIP
 *
 *  @see \ref SCIP_VerbLevel "SCIP_VERBLEVEL" for a list of all verbosity levels
 */
SCIP_EXPORT
SCIP_VERBLEVEL SCIPgetVerbLevel(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**@} */

#ifdef __cplusplus
}
#endif

#endif
