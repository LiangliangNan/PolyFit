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

/**@file   pub_message.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for message output
 * @author Tobias Achterberg
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MESSAGE_H__
#define __SCIP_PUB_MESSAGE_H__

#include <stdarg.h>
#include <string.h>

#include "scip/def.h"
#include "scip/type_message.h"

#ifdef NDEBUG
#include "scip/struct_message.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** define to identify SCIP version with thread-safe version of message handlers */
#define SCIP_THREADSAFE_MESSAGEHDLRS

/** define to get the filename of __FILE__ */
#if defined(_WIN32) || defined(_WIN64)
/*lint -e613*/
#define __FILENAME__ (strrchr("\\" __FILE__, '\\') + 1)
#else
/*lint -e613*/
#define __FILENAME__ (strrchr("/" __FILE__, '/') + 1)
#endif

/** prints an error message */
#define SCIPerrorMessage                SCIPmessagePrintErrorHeader(__FILENAME__, __LINE__); \
                                        SCIPmessagePrintError

/** define used in blockmemshell/memory.c */
#define printErrorHeader                SCIPmessagePrintErrorHeader
#define printError                      SCIPmessagePrintError

#ifdef SCIP_DEBUG

/** executes command only if SCIP_DEBUG flag is set */
#define SCIPdebug(x)                        x

/** prints a debugging message if SCIP_DEBUG flag is set - also consider using SCIPdebugMsg/SCIPsetDebugMsg */
#define SCIPdebugMessage                printf("[%s:%d] debug: ", __FILENAME__, __LINE__), printf

/** executes printf command only if SCIP_DEBUG flag is set */
#define SCIPdebugPrintf                 printf

/** executes SCIPprintCons() and prints termination symbol ";\n" only if SCIP_DEBUG flag is set */
#define SCIPdebugPrintCons(scip,cons,file) do                                                                   \
                                           {                                                                    \
                                              SCIP_CALL_ABORT( SCIPprintCons((scip), (cons), (file)) );         \
                                              SCIPinfoMessage((scip), (file), ";\n");                           \
                                           }                                                                    \
                                           while( FALSE )

#else

/** executes command only if SCIP_DEBUG flag is set */
#define SCIPdebug(x)                        /**/

/** prints a debugging message if SCIP_DEBUG flag is set - also consider using SCIPdebugMsg/SCIPsetDebugMsg */
#define SCIPdebugMessage                while( FALSE ) /*lint -e{530}*/ printf

/** executes printf command only if SCIP_DEBUG flag is set */
#define SCIPdebugPrintf                 while( FALSE ) /*lint -e{530}*/ printf

/** executes SCIPprintCons() and prints termination symbol ";\n" only if SCIP_DEBUG flag is set */
#define SCIPdebugPrintCons(x,y,z)           /**/

#endif

#ifdef SCIP_STATISTIC

/** executes command only if SCIP_STATISTIC flag is set */
#define SCIPstatistic(x)                        x

/** prints a statistic message if SCIP_STATISTIC flag is set */
#define SCIPstatisticMessage                printf("[%s:%d] statistic: ", __FILENAME__, __LINE__), printf

/** executes printf command only if SCIP_STATISTIC flag is set */
#define SCIPstatisticPrintf                 printf

#else

/** executes command only if SCIP_STATISTIC flag is set */
#define SCIPstatistic(x)                        /**/

/** prints a statistic message if SCIP_STATISTIC flag is set */
#define SCIPstatisticMessage                while( FALSE ) /*lint -e{530}*/ printf

/** executes printf command only if SCIP_STATISTIC flag is set */
#define SCIPstatisticPrintf                 while( FALSE ) /*lint -e{530}*/ printf

#endif


/** Creates and captures a message handler which deals with warning, information, and dialog (interactive shell) methods.
 *
 *  Use SCIPsetMessagehdlr() to make SCIP aware of the created message handler.
 *  @note The message handler does not handle error messages. For that see SCIPmessageSetErrorPrinting()
 *  @note Creating a message handler automatically captures it.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPmessagehdlrCreate(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   SCIP_Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   const char*           filename,           /**< name of log file, or NULL for no log */
   SCIP_Bool             quiet,              /**< should screen messages be suppressed? */
   SCIP_DECL_MESSAGEWARNING((*messagewarning)),/**< warning message print method of message handler */
   SCIP_DECL_MESSAGEDIALOG((*messagedialog)),/**< dialog message print method of message handler */
   SCIP_DECL_MESSAGEINFO ((*messageinfo)),   /**< info message print method of message handler */
   SCIP_DECL_MESSAGEHDLRFREE((*messagehdlrfree)), /**< destructor of message handler to free message handler data */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< message handler data */
   );

/** captures message handler */
SCIP_EXPORT
void SCIPmessagehdlrCapture(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler, or NULL */
   );

/** releases message handler */
SCIP_EXPORT
SCIP_RETCODE SCIPmessagehdlrRelease(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   );

/** sets the user data of the message handler */
SCIP_EXPORT
SCIP_RETCODE SCIPmessagehdlrSetData(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler; must not be NULL */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< new message handler data to attach to the handler */
   );

/** sets the log file name for the message handler */
SCIP_EXPORT
void SCIPmessagehdlrSetLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename            /**< log file name where to copy messages into, or NULL */
   );

/** sets the messages handler to be quiet */
SCIP_EXPORT
void SCIPmessagehdlrSetQuiet(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   );

/** prints a message, acting like the printf() command */
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
SCIP_EXPORT
void SCIPmessagePrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message, acting like the vprintf() command */
SCIP_EXPORT
void SCIPmessageVPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message into a file, acting like the fprintf() command */
#ifdef __GNUC__
__attribute__((format(printf, 3, 4)))
#endif
SCIP_EXPORT
void SCIPmessageFPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message into a file, acting like the vfprintf() command */
SCIP_EXPORT
void SCIPmessageVFPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a warning message, acting like the printf() command */
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
SCIP_EXPORT
void SCIPmessagePrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a warning message, acting like the vprintf() command */
SCIP_EXPORT
void SCIPmessageVPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a warning message into a file, acting like the fprintf() command */
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
SCIP_EXPORT
void SCIPmessageFPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a warning message into a file, acting like the vfprintf() command */
SCIP_EXPORT
void SCIPmessageVFPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a dialog message that requests user interaction, acting like the printf() command */
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
SCIP_EXPORT
void SCIPmessagePrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction, acting like the vprintf() command */
SCIP_EXPORT
void SCIPmessageVPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a dialog message that requests user interaction into a file, acting like the fprintf() command */
#ifdef __GNUC__
__attribute__((format(printf, 3, 4)))
#endif
SCIP_EXPORT
void SCIPmessageFPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a dialog message that requests user interaction into a file, acting like the vfprintf() command */
SCIP_EXPORT
void SCIPmessageVFPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message depending on the verbosity level, acting like the printf() command */
#ifdef __GNUC__
__attribute__((format(printf, 4, 5)))
#endif
SCIP_EXPORT
void SCIPmessagePrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message depending on the verbosity level, acting like the vprintf() command */
SCIP_EXPORT
void SCIPmessageVPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints a message into a file depending on the verbosity level, acting like the fprintf() command */
#ifdef __GNUC__
__attribute__((format(printf, 5, 6)))
#endif
SCIP_EXPORT
void SCIPmessageFPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints a message into a file depending on the verbosity level, acting like the vfprintf() command */
SCIP_EXPORT
void SCIPmessageVFPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** prints the header with source file location for an error message using the static message handler */
SCIP_EXPORT
void SCIPmessagePrintErrorHeader(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   );

/** prints an error message, acting like the printf() command using the static message handler */
#ifdef __GNUC__
__attribute__((format(printf, 1, 2)))
#endif
SCIP_EXPORT
void SCIPmessagePrintError(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   );

/** prints an error message, acting like the vprintf() command using the static message handler */
SCIP_EXPORT
void SCIPmessageVPrintError(
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   );

/** Method to set the error printing method. Setting the error printing method to NULL will suspend all error methods.
 *
 *  @note The error printing method is a static variable. This means that all occurring errors are handled via this method.
 */
SCIP_EXPORT
void SCIPmessageSetErrorPrinting(
   SCIP_DECL_ERRORPRINTING((*errorPrinting)),/**< error message print method of message handler, or NULL */
   void*                 data                /**< data pointer which will be passed to the error printing method, or NULL */
   );

/** Method to set the error printing method to default version prints everything the stderr.
 *
 *  @note The error printing method is a static variable. This means that all occurring errors are handled via this method.
 */
SCIP_EXPORT
void SCIPmessageSetErrorPrintingDefault(
   void
   );


/** returns the user data of the message handler */
SCIP_EXPORT
SCIP_MESSAGEHDLRDATA* SCIPmessagehdlrGetData(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** returns the log file or NULL for stdout */
SCIP_EXPORT
FILE* SCIPmessagehdlrGetLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** returns TRUE if the message handler is set to be quiet */
SCIP_EXPORT
SCIP_Bool SCIPmessagehdlrIsQuiet(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPmessagehdlrGetData(messagehdlr)     ((messagehdlr) != NULL) ? messagehdlr->messagehdlrdata : NULL
#define SCIPmessagehdlrGetLogfile(messagehdlr)  ((messagehdlr) == NULL ? NULL : (messagehdlr)->logfile)
#define SCIPmessagehdlrIsQuiet(messagehdlr)     ((messagehdlr) == NULL || (messagehdlr)->quiet)

#endif

#ifdef __cplusplus
}
#endif

#endif
