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

/**@file   scip_message.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for message handling
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/scip_message.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"

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
SCIP_RETCODE SCIPsetMessagehdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetMessagehdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE) );

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPmessagehdlrCapture(messagehdlr);

   SCIP_CALL( SCIPmessagehdlrRelease(&scip->messagehdlr) );
   assert(scip->messagehdlr == NULL);

   scip->messagehdlr = messagehdlr;

   return SCIP_OKAY;
}

/** returns the currently installed message handler
 *
 *  @return the currently installed message handler, or NULL if messages are currently suppressed
 */
SCIP_MESSAGEHDLR* SCIPgetMessagehdlr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return scip->messagehdlr;
}

/** sets the log file name for the currently installed message handler */
void SCIPsetMessagehdlrLogfile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of log file, or NULL (no log) */
   )
{
   if( scip->messagehdlr != NULL )
   {
      SCIPmessagehdlrSetLogfile(scip->messagehdlr, filename);
   }
}

/** sets the currently installed message handler to be quiet (or not) */
void SCIPsetMessagehdlrQuiet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   )
{
   if( scip->messagehdlr != NULL )
   {
      SCIPmessagehdlrSetQuiet(scip->messagehdlr, quiet);
   }
}

/** prints a warning message via the message handler */
void SCIPwarningMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(scip != NULL);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintWarning(scip->messagehdlr, formatstr, ap);
   va_end(ap);
}

/** prints a debug message */
void SCIPprintDebugMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   const char* filename;
   int subscipdepth = 0;
   va_list ap;

   assert( sourcefile != NULL );
   assert( scip != NULL );

   /* strip directory from filename */
#if defined(_WIN32) || defined(_WIN64)
   filename = strrchr(sourcefile, '\\');
#else
   filename = strrchr(sourcefile, '/');
#endif
   if ( filename == NULL )
      filename = sourcefile;
   else
      ++filename;

   if ( scip->stat != NULL )
      subscipdepth = scip->stat->subscipdepth;
   if ( subscipdepth > 0 )
      SCIPmessageFPrintInfo(scip->messagehdlr, NULL, "%d: [%s:%d] debug: ", subscipdepth, filename, sourceline);
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, NULL, "[%s:%d] debug: ", filename, sourceline);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(scip->messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a debug message without precode */
void SCIPdebugMessagePrint(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert( scip != NULL );

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(scip->messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction or is a direct response to a user interactive command */
void SCIPdialogMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(scip != NULL);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintDialog(scip->messagehdlr, file, formatstr, ap);
   va_end(ap);
}

/** prints a message */
void SCIPinfoMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(scip != NULL);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(scip->messagehdlr, file, formatstr, ap);
   va_end(ap);
}

/** prints a message depending on the verbosity level */
void SCIPverbMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert(scip != NULL);
   assert(scip->set != NULL);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintVerbInfo(scip->messagehdlr, scip->set->disp_verblevel, msgverblevel, file, formatstr, ap);
   va_end(ap);
}

/** returns the current message verbosity level
 *
 *  @return message verbosity level of SCIP
 *
 *  @see \ref SCIP_VerbLevel "SCIP_VERBLEVEL" for a list of all verbosity levels
 */
SCIP_VERBLEVEL SCIPgetVerbLevel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->disp_verblevel;
}
