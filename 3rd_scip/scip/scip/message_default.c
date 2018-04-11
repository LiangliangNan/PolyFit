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
