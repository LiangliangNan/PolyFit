/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This1 file is part of the program and library             */
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

/**@file   message.c
 * @brief  message output methods
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdarg.h>
#include <stdio.h>
#include <assert.h>

#include "scip/type_message.h"
#include "scip/struct_message.h"
#include "scip/def.h"
#include "scip/pub_misc.h"
#include "blockmemshell/memory.h"


#ifndef va_copy
#define va_copy(dest, src) do { BMScopyMemory(&dest, &src); } while( 0 )
#endif

/* do defines for windows directly her to make the lpi more independent*/
#if defined(_WIN32) || defined(_WIN64)
#define snprintf _snprintf
#define vsnprintf _vsnprintf
#endif

/** handles the output of the given message */
static
void handleMessage(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_DECL_MESSAGEOUTPUTFUNC(outputfunc),  /**< message handler function used for output */
   FILE*                 file1,              /**< file stream to print into, or NULL for stdout */
   SCIP_Bool             usefile1,           /**< Should file1 be used? */
   FILE*                 file2,              /**< file stream to print into */
   SCIP_Bool             usefile2,           /**< Should file2 be used? */
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   char*                 buffer,             /**< message buffer */
   int*                  bufferlen           /**< pointer to the currently used entries in the message buffer */
   )
{
   const char* s;

   assert( messagehdlr != NULL );
   assert( outputfunc != NULL );
   assert( !usefile2 || file2 != NULL );
   assert( buffer == NULL || bufferlen != NULL );

   /* if we do not have a buffer directly output the message */
   if ( buffer == NULL )
   {
      /* we do not have a buffer, so it makes no sense to flush it if msg == NULL */
      if ( msg != NULL )
      {
         if ( usefile1 )
            outputfunc(messagehdlr, file1, msg);
         if ( usefile2 )
            outputfunc(messagehdlr, file2, msg);
      }
      return;
   }
   assert(bufferlen != NULL);

   /* should the buffer be flushed? */
   if ( msg == NULL )
   {
      assert( *bufferlen < SCIP_MAXSTRLEN );
      assert( buffer[*bufferlen] == '\0' );
      if ( usefile1 )
         outputfunc(messagehdlr, file1, buffer);
      if ( usefile2 )
         outputfunc(messagehdlr, file2, buffer);
      *bufferlen = 0;
      buffer[0] = '\0';
      return;
   }
   assert( msg != NULL && buffer != NULL );

   /* if no output is activated, to not copy message into buffer */
   if ( ! usefile1 && ! usefile2 )
      return;

   /* determine message length and last newline (if any) */
   s = msg;
   while ( *s != '\0' )
   {
      /* if we reached a newline or the size limit, empty buffer and reset (need possibly space for newline and '\0') */
      if ( *s == '\n' || *bufferlen >= SCIP_MAXSTRLEN-2 )
      {
         if ( *s == '\n' )
            buffer[(*bufferlen)++] = *(s++);
         buffer[*bufferlen] = '\0';

         if ( usefile1 )
            outputfunc(messagehdlr, file1, buffer);
         if ( usefile2 )
            outputfunc(messagehdlr, file2, buffer);
         *bufferlen = 0;
         buffer[0] = '\0';
      }
      else
         buffer[(*bufferlen)++] = *(s++);
   }
   buffer[*bufferlen] = '\0';

   return;
}

/** default error printing method which is used to print all occurring errors */
static
SCIP_DECL_ERRORPRINTING(errorPrintingDefault)
{  /*lint --e{715}*/
   if ( msg != NULL )
   {
      if ( file != NULL )
         fputs(msg, file);
      else
         fputs(msg, stderr);
   }
   fflush(stderr);
}

/** static variable which holds the error printing method */
static SCIP_DECL_ERRORPRINTING((*staticErrorPrinting)) = errorPrintingDefault;

/** static variable which holds a data pointer for the error prinint callback */
static void* staticErrorPrintingData = NULL;

/** prints error message with the current static message handler */
static
void messagePrintError(
   FILE*                 file,               /**< file stream to print error, or NULL for stderr */
   const char*           msg                 /**< message to print; NULL to flush the output buffer */
   )
{
   if( staticErrorPrinting != NULL )
      staticErrorPrinting(staticErrorPrintingData, file, msg);
}

/** prints warning message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           msg                 /**< message to print; NULL to flush the output buffer */
   )
{  /*lint --e{715}*/
   if ( messagehdlr != NULL && messagehdlr->messagewarning != NULL && (! messagehdlr->quiet || messagehdlr->logfile != NULL) )
   {
      handleMessage(messagehdlr, messagehdlr->messagewarning, stderr, ! messagehdlr->quiet, messagehdlr->logfile, (messagehdlr->logfile != NULL),
         msg, messagehdlr->warningbuffer, &messagehdlr->warningbufferlen);
   }
}

/** prints dialog message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           msg                 /**< message to print; NULL to flush the output buffer */
   )
{  /*lint --e{715}*/
   if ( messagehdlr != NULL && messagehdlr->messagedialog != NULL )
   {
      if ( (file == NULL || file == stdout) && ! messagehdlr->quiet )
      {
         handleMessage(messagehdlr, messagehdlr->messagedialog, (file == NULL) ? stdout : file, TRUE, messagehdlr->logfile, (messagehdlr->logfile != NULL),
            msg, messagehdlr->dialogbuffer, &messagehdlr->dialogbufferlen);
      }
      else if ( msg != NULL )
      {
         /* file output cannot be buffered because the output file may change */
         if ( *msg != '\0' )
         {
            handleMessage(messagehdlr, messagehdlr->messagedialog, file, !messagehdlr->quiet || (file != NULL && file != stdout), messagehdlr->logfile, (messagehdlr->logfile != NULL), msg, NULL, NULL);
         }
      }
   }
}

/** prints info message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           msg                 /**< message to print; NULL to flush the output buffer */
   )
{  /*lint --e{715}*/
   if ( messagehdlr != NULL && messagehdlr->messageinfo != NULL )
   {
      if ( (file == NULL || file == stdout) && ! messagehdlr->quiet )
      {
         handleMessage(messagehdlr, messagehdlr->messageinfo, (file == NULL) ? stdout : file, TRUE, messagehdlr->logfile, (messagehdlr->logfile != NULL),
            msg, messagehdlr->infobuffer, &messagehdlr->infobufferlen);
      }
      else if ( msg != NULL )
      {
         /* file output cannot be buffered because the output file may change or the message is to long */
         if ( *msg != '\0' )
         {
            handleMessage(messagehdlr, messagehdlr->messagedialog, file, !messagehdlr->quiet || (file != NULL && file != stdout), messagehdlr->logfile, (messagehdlr->logfile != NULL), msg, NULL, NULL);
         }
      }
   }
}

/** if the given file is not NULL a log file is opened */
static
void messagehdlrOpenLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename            /**< name of log file, or NULL (stdout) */
   )
{
   if( filename != NULL )
   {
      messagehdlr->logfile = fopen(filename, "a"); /* append to log file */

      if( messagehdlr->logfile == NULL )
      {
         SCIPerrorMessage("cannot open log file <%s> for writing\n", filename);
      }
   }
   else
      messagehdlr->logfile = NULL;
}

/** frees message handler */
static
SCIP_RETCODE messagehdlrFree(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   )
{
   assert(messagehdlr != NULL);

   if( *messagehdlr != NULL )
   {
      /* flush message buffers */
      messagePrintWarning(*messagehdlr, NULL);
      messagePrintDialog(*messagehdlr, NULL, NULL);
      messagePrintInfo(*messagehdlr, NULL, NULL);

      if( (*messagehdlr)->messagehdlrfree != NULL )
      {
         /* call destructor method of message handler to free the message handler data */
         SCIP_CALL( (*messagehdlr)->messagehdlrfree(*messagehdlr) );
      }

      /* close the log file if one exists */
      if( (*messagehdlr)->logfile != NULL )
      {
         fclose((*messagehdlr)->logfile);
      }

      /* free buffer arrays */
      BMSfreeMemoryArrayNull(&(*messagehdlr)->warningbuffer);
      BMSfreeMemoryArrayNull(&(*messagehdlr)->dialogbuffer);
      BMSfreeMemoryArrayNull(&(*messagehdlr)->infobuffer);
      BMSfreeMemory(messagehdlr);
   }

   return SCIP_OKAY;
}

/** Creates and captures a message handler which deals with warning, information, and dialog (interactive shell) methods.
 *
 *  @note The message handler does not handle error messages; see SCIPmessageSetErrorPrinting()
 */
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
   )
{
   SCIP_ALLOC( BMSallocMemory(messagehdlr) );
   (*messagehdlr)->messagewarning = messagewarning;
   (*messagehdlr)->messagedialog = messagedialog;
   (*messagehdlr)->messageinfo = messageinfo;
   (*messagehdlr)->messagehdlrfree = messagehdlrfree;
   (*messagehdlr)->messagehdlrdata = messagehdlrdata;
   (*messagehdlr)->warningbuffer = NULL;
   (*messagehdlr)->dialogbuffer = NULL;
   (*messagehdlr)->infobuffer = NULL;
   (*messagehdlr)->warningbufferlen = 0;
   (*messagehdlr)->dialogbufferlen = 0;
   (*messagehdlr)->infobufferlen = 0;
   (*messagehdlr)->nuses = 1;

   (*messagehdlr)->quiet = quiet;
   messagehdlrOpenLogfile(*messagehdlr, filename);

   /* allocate buffer for buffered output */
   if( bufferedoutput )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->warningbuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->dialogbuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->infobuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      (*messagehdlr)->warningbuffer[0] = '\0';
      (*messagehdlr)->dialogbuffer[0] = '\0';
      (*messagehdlr)->infobuffer[0] = '\0';
   }

   return SCIP_OKAY;
}

/** captures message handler */
void SCIPmessagehdlrCapture(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler, or NULL */
   )
{
   if( messagehdlr != NULL )
      ++messagehdlr->nuses;
}

/** releases message handler */
SCIP_RETCODE SCIPmessagehdlrRelease(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   )
{
   assert(messagehdlr != NULL);

   if( *messagehdlr == NULL )
      return SCIP_OKAY;

   assert((*messagehdlr)->nuses >= 1);

   /* decrement usage counter */
   --(*messagehdlr)->nuses;

   /* the last one turns the light off */
   if( (*messagehdlr)->nuses == 0 )
   {
      SCIP_CALL( messagehdlrFree(messagehdlr) );
      assert(*messagehdlr == NULL);
   }
   else
   {
      *messagehdlr = NULL;
   }

   return SCIP_OKAY;
}

/** sets the user data of the message handler */
SCIP_RETCODE SCIPmessagehdlrSetData(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler; must not be NULL */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< new message handler data to attach to the handler */
   )
{
   assert(messagehdlr != NULL);

   if( messagehdlr == NULL ) /*lint !e774*/
      return SCIP_INVALIDDATA;

   messagehdlr->messagehdlrdata = messagehdlrdata;

   return SCIP_OKAY;
}

/** sets the log file name for the message handler */
void SCIPmessagehdlrSetLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename            /**< log file name where to copy messages into, or NULL */
   )
{
   assert(messagehdlr != NULL);

   /* close the old log file if one exists */
   if( messagehdlr->logfile != NULL )
   {
      fclose(messagehdlr->logfile);
   }

   /* opens the log file */
   messagehdlrOpenLogfile(messagehdlr, filename);
}

/** sets the messages handler to be quiet */
void SCIPmessagehdlrSetQuiet(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   )
{
   assert(messagehdlr != NULL);

   /* flush message buffers in order to not loose information */
   messagePrintWarning(messagehdlr, NULL);
   messagePrintDialog(messagehdlr, NULL, NULL);
   messagePrintInfo(messagehdlr, NULL, NULL);

   messagehdlr->quiet = quiet;
}

/** prints a warning message, acting like the printf() command */
void SCIPmessagePrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintWarning(messagehdlr, formatstr, ap);
   va_end(ap);
}

/** prints a warning message, acting like the vprintf() command */
void SCIPmessageVPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintWarning(messagehdlr, formatstr, ap);
}

/** prints a warning message, acting like the fprintf() command */
void SCIPmessageFPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintWarning(messagehdlr, formatstr, ap);
   va_end(ap);
}

/** prints a warning message, acting like the vfprintf() command */
void SCIPmessageVFPrintWarning(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;

   va_copy(aq, ap); /*lint !e838*/

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#endif
      assert(m == n);
      va_end(aq);
      messagePrintWarning(messagehdlr, bigmsg);
      BMSfreeMemory(&bigmsg);
      return;
   }

   messagePrintWarning(messagehdlr, msg);
   va_end(aq);
}

/** prints a dialog message that requests user interaction, acting like the printf() command */
void SCIPmessagePrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintDialog(messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction, acting like the vprintf() command */
void SCIPmessageVPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintDialog(messagehdlr, NULL, formatstr, ap);
}

/** prints a dialog message that requests user interaction into a file, acting like the fprintf() command */
void SCIPmessageFPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintDialog(messagehdlr, file, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction into a file, acting like the vfprintf() command */
void SCIPmessageVFPrintDialog(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;

   va_copy(aq, ap); /*lint !e838*/

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#endif
      assert(m == n);
      va_end(aq);
      messagePrintDialog(messagehdlr, file, bigmsg);
      BMSfreeMemory(&bigmsg);
      return;
   }
   messagePrintDialog(messagehdlr, file, msg);
   va_end(aq);
}

/** prints a message, acting like the printf() command */
void SCIPmessagePrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a message, acting like the vprintf() command */
void SCIPmessageVPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintInfo(messagehdlr, NULL, formatstr, ap);
}

/** prints a message into a file, acting like the fprintf() command */
void SCIPmessageFPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(messagehdlr, file, formatstr, ap);
   va_end(ap);
}

/** prints a message into a file, acting like the vfprintf() command */
void SCIPmessageVFPrintInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;

   va_copy(aq, ap); /*lint !e838*/

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#endif
      assert(m == n);
      va_end(aq);
      messagePrintInfo(messagehdlr, file, bigmsg);
      BMSfreeMemory(&bigmsg);
      return;
   }
   messagePrintInfo(messagehdlr, file, msg);
   va_end(aq);
}

/** prints a message depending on the verbosity level, acting like the printf() command */
void SCIPmessagePrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintVerbInfo(messagehdlr, verblevel, msgverblevel, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a message depending on the verbosity level, acting like the vprintf() command */
void SCIPmessageVPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintVerbInfo(messagehdlr, verblevel, msgverblevel, NULL, formatstr, ap);
}

/** prints a message into a file depending on the verbosity level, acting like the fprintf() command */
void SCIPmessageFPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintVerbInfo(messagehdlr, verblevel, msgverblevel, file, formatstr, ap);
   va_end(ap);
}

/** prints a message into a file depending on the verbosity level, acting like the vfprintf() command */
void SCIPmessageVFPrintVerbInfo(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   assert(msgverblevel > SCIP_VERBLEVEL_NONE);
   assert(msgverblevel <= SCIP_VERBLEVEL_FULL);
   assert(verblevel <= SCIP_VERBLEVEL_FULL);

   if( msgverblevel <= verblevel )
   {
      char msg[SCIP_MAXSTRLEN];
      int n;
      va_list aq;

      va_copy(aq, ap); /*lint !e838*/

      n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
      if( n < 0 )
         msg[SCIP_MAXSTRLEN-1] = '\0';
      else if( n >= SCIP_MAXSTRLEN )
      {
         char* bigmsg;
#ifndef NDEBUG
         int m;
#endif

         if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
         {
            va_end(aq);
            return;
         }

#ifndef NDEBUG
         m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#else
         vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#endif
         assert(m == n);
         va_end(aq);
         messagePrintInfo(messagehdlr, file, bigmsg);
         BMSfreeMemory(&bigmsg);
         return;
      }
      messagePrintInfo(messagehdlr, file, msg);
      va_end(aq);
   }
}

/** prints the header with source file location for an error message using the static message handler */
void SCIPmessagePrintErrorHeader(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   )
{
   char msg[SCIP_MAXSTRLEN];

   /* safe string printing - do not use SCIPsnprintf() since message.c should be independent */
   (void) snprintf(msg, SCIP_MAXSTRLEN, "[%s:%d] ERROR: ", sourcefile, sourceline);
   msg[SCIP_MAXSTRLEN-1] = '\0';
   messagePrintError(NULL, msg);
}

/** prints a error message, acting like the printf() command */
void SCIPmessagePrintError(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVPrintError(formatstr, ap);
   va_end(ap);
}

/** prints an error message, acting like the vprintf() command using the static message handler */
void SCIPmessageVPrintError(
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;

   va_copy(aq, ap); /*lint !e838*/

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e571*/
#endif
      assert(m == n);
      va_end(aq);
      messagePrintError(NULL, bigmsg);
      BMSfreeMemory(&bigmsg);
      return;
   }

   messagePrintError(NULL, msg);
   va_end(aq);
}

/** Method to set the error printing method. Setting the error printing method to NULL will suspend all error methods.
 *
 *  @note The error printing method is static variable. That means all occurring errors are handled via that methods
 */
void SCIPmessageSetErrorPrinting(
   SCIP_DECL_ERRORPRINTING((*errorPrinting)),/**< error message print method of message handler, or NULL */
   void*                 data                /**< data pointer which will be passed to the error printing method, or NULL */
   )
{
   staticErrorPrinting = errorPrinting;
   staticErrorPrintingData = data;
}

/** Method to set the error printing method to default version prints everything the stderr.
 *
 *  @note The error printing method is a static variable. This means that all occurring errors are handled via this method.
 */
void SCIPmessageSetErrorPrintingDefault(
   void
   )
{
   staticErrorPrinting = errorPrintingDefault;
   staticErrorPrintingData = NULL;
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPmessagehdlrGetData
#undef SCIPmessagehdlrGetLogfile
#undef SCIPmessagehdlrIsQuiet

/** returns the user data of the message handler */
SCIP_MESSAGEHDLRDATA* SCIPmessagehdlrGetData(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   if( messagehdlr != NULL )
      return messagehdlr->messagehdlrdata;
   else
      return NULL;
}


/** returns the log file or NULL for stdout */
FILE* SCIPmessagehdlrGetLogfile(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   if( messagehdlr == NULL )
      return NULL;

   return messagehdlr->logfile;
}

/** returns TRUE if the message handler is set to be quiet */
SCIP_Bool SCIPmessagehdlrIsQuiet(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   return (messagehdlr == NULL || messagehdlr->quiet);
}
