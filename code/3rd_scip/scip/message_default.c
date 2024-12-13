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

/**@file   message_default.c
 * @ingroup PUBLICMETHODS
 * @brief  default message handler
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/message_default.h"
#include "scip/struct_message.h"

/*
 * Local methods
 */

/** prints a message to the given file stream and writes the same messate to the log file */
static
void logMessage(
   FILE*                 file,               /**< file stream to print message into */
   const char*           msg                 /**< message to print (or NULL to flush) */
   )
{
   if ( msg != NULL )
      fputs(msg, file);
   fflush(file);
}

/*
 * Callback methods of message handler
 */

/** warning message print method of message handler */
static
SCIP_DECL_MESSAGEWARNING(messageWarningDefault)
{  /*lint --e{715}*/
   if ( msg != NULL && msg[0] != '\0' && msg[0] != '\n' )
      fputs("WARNING: ", file);

   logMessage(file, msg);
}

/** dialog message print method of message handler */
static
SCIP_DECL_MESSAGEDIALOG(messageDialogDefault)
{  /*lint --e{715}*/
   logMessage(file, msg);
}

/** info message print method of message handler */
static
SCIP_DECL_MESSAGEINFO(messageInfoDefault)
{  /*lint --e{715}*/
   logMessage(file, msg);
}

/** Create default message handler. To free the message handler use SCIPmessagehdlrRelease(). */
SCIP_RETCODE SCIPcreateMessagehdlrDefault(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store message handler */
   SCIP_Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   const char*           filename,           /**< name of log file, or NULL (stdout) */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   )
{
   /* create message handler */
   SCIP_CALL( SCIPmessagehdlrCreate(messagehdlr, bufferedoutput, filename, quiet,
         messageWarningDefault, messageDialogDefault, messageInfoDefault,
         NULL, NULL) );

   return SCIP_OKAY;
}
