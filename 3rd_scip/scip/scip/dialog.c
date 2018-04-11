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

/**@file   dialog.c
 * @brief  methods for user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#ifdef WITH_READLINE
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include "scip/scip.h"
#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/pub_misc.h"
#include "scip/dialog.h"

#include "scip/struct_dialog.h"




/*
 * read line methods
 */

#ifdef WITH_READLINE

/** reads a line of input from stdin */
static
SCIP_RETCODE readLine(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   const char*           prompt,             /**< prompt to display */
   SCIP_Bool*            endoffile           /**< pointer to store whether the end of the input file was reached */
   )
{
   char* s;

   assert(endoffile != NULL);

   s = readline(prompt);
   if( s != NULL )
   {
      (void)strncpy(&dialoghdlr->buffer[dialoghdlr->bufferpos], s,
         (unsigned int)(dialoghdlr->buffersize - dialoghdlr->bufferpos));
      free(s);
      *endoffile = FALSE;
   }
   else
      *endoffile = TRUE;

   return SCIP_OKAY;
}

/** puts the given string on the command history */
static
SCIP_RETCODE addHistory(
   const char*           s                   /**< string to add to the command history */
   )
{
   add_history(s);

   return SCIP_OKAY;
}

/** returns the current length of the history list */
static
int getHistoryLength(
   void
   )
{
#ifndef NO_REMOVE_HISTORY
   return history_length;
#else
   return 0;
#endif
}

/** removes a single element from the history list */
static
SCIP_RETCODE removeHistory(
   int                   pos                 /**< list position of history entry to remove */
   )
{
#ifndef NO_REMOVE_HISTORY
   HIST_ENTRY* entry;

   entry = remove_history(pos);

   /* Free readline/history storage: there seem to be differences in the versions (and the amount of
    * data to be freed). The following should be a good approximation; if it doesn't work define
    * NO_REMOVE_HISTORY - see the INSTALL file. This will produce minor memory leaks.
    */
#if RL_VERSION_MAJOR >= 5
   (void)free_history_entry(entry);
#else
   if( entry != NULL )
   {
      free((void*)entry->line);
      free(entry);
   }
#endif
#endif

   return SCIP_OKAY;
}

/** writes command history into file of the specified name */
static
SCIP_RETCODE writeHistory(
   const char*           filename            /**< name of file to (over)write history to */
   )
{
   int retval = write_history(filename);

   if( retval == 0 )
      return SCIP_OKAY;
   else
      return SCIP_FILECREATEERROR;
}

#else

/** reads a line of input from stdin */
static
SCIP_RETCODE readLine(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   const char*           prompt,             /**< prompt to display */
   SCIP_Bool*            endoffile           /**< pointer to store whether the end of the input file was reached */
   )
{
   char* s;

   assert(dialoghdlr != NULL);
   assert(dialoghdlr->buffer != NULL);
   assert(dialoghdlr->bufferpos < dialoghdlr->buffersize);
   assert(dialoghdlr->buffer[dialoghdlr->bufferpos] == '\0');
   assert(endoffile != NULL);

   /* check for EOF (due to CTRL-D or unexpected end of piped-in file) */
   if( feof(stdin) )
      *endoffile = TRUE;
   else
   {
      char* result;

      /* display prompt */
      printf("%s", prompt);

      /* read line from stdin */
      result = fgets(&dialoghdlr->buffer[dialoghdlr->bufferpos], dialoghdlr->buffersize - dialoghdlr->bufferpos, stdin);
      assert(result != NULL);
      (void) result; /* disable compiler warning [-Wunused-result] */

      /* replace newline with \0 */
      s = strchr(&dialoghdlr->buffer[dialoghdlr->bufferpos], '\n');
      if( s != NULL )
         *s = '\0';
      *endoffile = FALSE;
   }

   return SCIP_OKAY;
}

/** puts the given string on the command history */
static
SCIP_RETCODE addHistory(
   const char*           s                   /**< string to add to the command history */
   )
{  /*lint --e{715}*/
   /* nothing to do here */
   return SCIP_OKAY;
}

/** returns the current length of the history list */
static
int getHistoryLength(
   void
   )
{
   return 0;
}

/** removes a single element from the history list */
static
SCIP_RETCODE removeHistory(
   int                   pos                 /**< list position of history entry to remove */
   )
{  /*lint --e{715}*/
   /* nothing to do here */
   return SCIP_OKAY;
}


/** writes command history into file of the specified name */
static
SCIP_RETCODE writeHistory(
   const char*           filename            /**< name of file to (over)write history to */
   )
{  /*lint --e{715}*/
   /* nothing to do here */
   return SCIP_OKAY;
}

#endif

/** frees a single linelist entry, but not its successors */
static
void linelistFree(
   SCIP_LINELIST**       linelist            /**< pointer to line list */
   )
{
   assert(linelist != NULL);

   BMSfreeMemoryArray(&(*linelist)->inputline);
   BMSfreeMemory(linelist);
}

/** frees a linelist entry and all of its successors */
static
void linelistFreeAll(
   SCIP_LINELIST**       linelist            /**< pointer to line list */
   )
{
   assert(linelist != NULL);

   while( *linelist != NULL )
   {
      SCIP_LINELIST* nextline;

      nextline = (*linelist)->nextline;
      linelistFree(linelist);
      *linelist = nextline;
   }
}

/** reads a line of input from stdin or from the stored input lines in the input list */
static
SCIP_RETCODE readInputLine(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   const char*           prompt,             /**< prompt to display */
   SCIP_Bool*            endoffile           /**< pointer to store whether the end of the input file was reached */
   )
{
   assert(dialoghdlr != NULL);
   assert(dialoghdlr->buffer != NULL);
   assert(dialoghdlr->bufferpos < dialoghdlr->buffersize);
   assert(dialoghdlr->buffer[dialoghdlr->bufferpos] == '\0');
   assert(endoffile != NULL);

   *endoffile = FALSE;

   if( dialoghdlr->inputlist == NULL )
   {
      /* read a line from stdin */
      SCIP_CALL( readLine(dialoghdlr, prompt, endoffile) );
   }
   else
   {
      SCIP_LINELIST* nextline;

      /* copy the next input line into the input buffer */
      (void)strncpy(&dialoghdlr->buffer[dialoghdlr->bufferpos], dialoghdlr->inputlist->inputline,
         (size_t)(dialoghdlr->buffersize - dialoghdlr->bufferpos)); /*lint !e571 !e776*/
      dialoghdlr->buffer[dialoghdlr->buffersize-1] = '\0';

      /* free the input line */
      nextline = dialoghdlr->inputlist->nextline;
      if( dialoghdlr->inputlistptr == &(dialoghdlr->inputlist->nextline) )
         dialoghdlr->inputlistptr = &dialoghdlr->inputlist;
      linelistFree(&dialoghdlr->inputlist);
      dialoghdlr->inputlist = nextline;
      assert(dialoghdlr->inputlistptr != NULL);
      assert(*dialoghdlr->inputlistptr == NULL);
   }

   return SCIP_OKAY;
}




/*
 * dialog handler
 */

/** copies the given dialog to a new scip */
SCIP_RETCODE SCIPdialogCopyInclude(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(dialog != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( dialog->dialogcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including dialog %s in subscip %p\n", SCIPdialogGetName(dialog), (void*)set->scip);
      SCIP_CALL( dialog->dialogcopy(set->scip, dialog) );
   }
   return SCIP_OKAY;
}

/** creates a dialog handler */
SCIP_RETCODE SCIPdialoghdlrCreate(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DIALOGHDLR**     dialoghdlr          /**< pointer to store dialog handler */
   )
{  /*lint --e{715}*/
#ifdef WITH_READLINE
   char readlineversion[20];
#endif

   assert(dialoghdlr != NULL);

   SCIP_ALLOC( BMSallocMemory(dialoghdlr) );
   (*dialoghdlr)->rootdialog = NULL;
   (*dialoghdlr)->inputlist = NULL;
   (*dialoghdlr)->inputlistptr = &(*dialoghdlr)->inputlist;
   (*dialoghdlr)->buffersize = SCIP_MAXSTRLEN;
   (*dialoghdlr)->nprotectedhistelems = -1;
   SCIP_ALLOC( BMSallocMemoryArray(&(*dialoghdlr)->buffer, (*dialoghdlr)->buffersize) );

   SCIPdialoghdlrClearBuffer(*dialoghdlr);

#ifdef WITH_READLINE
   (void) SCIPsnprintf(readlineversion, sizeof(readlineversion), "Readline %s", rl_library_version);
   SCIP_CALL( SCIPsetIncludeExternalCode(set, readlineversion, "GNU library for command line editing (gnu.org/s/readline)") );
#endif

   return SCIP_OKAY;
}

/** frees a dialog handler and it's dialog tree */
SCIP_RETCODE SCIPdialoghdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOGHDLR**     dialoghdlr          /**< pointer to dialog handler */
   )
{
   assert(dialoghdlr != NULL);
   if( *dialoghdlr == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPdialoghdlrSetRoot(scip, *dialoghdlr, NULL) );
   linelistFreeAll(&(*dialoghdlr)->inputlist);
   BMSfreeMemoryArray(&(*dialoghdlr)->buffer);
   BMSfreeMemory(dialoghdlr);

   return SCIP_OKAY;
}

/** executes the root dialog of the dialog handler */
SCIP_RETCODE SCIPdialoghdlrExec(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_DIALOG* dialog;

   assert(dialoghdlr != NULL);
   assert(dialoghdlr->buffer != NULL);

   /* clear the buffer, start with the root dialog */
   SCIPdialoghdlrClearBuffer(dialoghdlr);
   dialog = dialoghdlr->rootdialog;

   /* execute dialogs until a NULL is returned as next dialog */
   while( dialog != NULL )
   {
      SCIP_CALL( SCIPdialogExec(dialog, set, dialoghdlr, &dialog) );

      /* reset buffer, it is was consumed completely */
      if( dialoghdlr->buffer[dialoghdlr->bufferpos] == '\0' )
         SCIPdialoghdlrClearBuffer(dialoghdlr);
   }

   return SCIP_OKAY;
}

/** makes given dialog the root dialog of dialog handler; captures dialog and releases former root dialog */
SCIP_RETCODE SCIPdialoghdlrSetRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG*          dialog              /**< dialog to be the root */
   )
{
   assert(dialoghdlr != NULL);

   if( dialoghdlr->rootdialog != NULL )
   {
      SCIP_CALL( SCIPdialogRelease(scip, &dialoghdlr->rootdialog) );
   }
   assert(dialoghdlr->rootdialog == NULL);

   dialoghdlr->rootdialog = dialog;

   if( dialog != NULL )
      SCIPdialogCapture(dialog);

   return SCIP_OKAY;
}

/** returns the root dialog of the dialog handler */
SCIP_DIALOG* SCIPdialoghdlrGetRoot(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   )
{
   assert(dialoghdlr != NULL);

   return dialoghdlr->rootdialog;
}

/** clears the input command buffer of the dialog handler */
void SCIPdialoghdlrClearBuffer(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   )
{
   assert(dialoghdlr != NULL);

   dialoghdlr->buffer[0] = '\0';
   dialoghdlr->bufferpos = 0;
}

/** returns TRUE iff input command buffer is empty */
SCIP_Bool SCIPdialoghdlrIsBufferEmpty(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   )
{
   assert(dialoghdlr != NULL);
   assert(dialoghdlr->bufferpos < dialoghdlr->buffersize);

   return (dialoghdlr->buffer[dialoghdlr->bufferpos] == '\0');
}

/** returns the next line in the handler's command buffer; if the buffer is empty, displays the given prompt or the
 *  current dialog's path and asks the user for further input; the user must not free or modify the returned string
 */
SCIP_RETCODE SCIPdialoghdlrGetLine(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG*          dialog,             /**< current dialog */
   const char*           prompt,             /**< prompt to display, or NULL to display the current dialog's path */
   char**                inputline,          /**< pointer to store the complete line in the handler's command buffer */
   SCIP_Bool*            endoffile           /**< pointer to store whether the end of the input file was reached */
   )
{
   char path[SCIP_MAXSTRLEN];
   char p[SCIP_MAXSTRLEN];

   assert(dialoghdlr != NULL);
   assert(dialoghdlr->buffer != NULL);
   assert(dialoghdlr->bufferpos < dialoghdlr->buffersize);
   assert(inputline != NULL);

   /* get input from the user, if the buffer is empty */
   if( SCIPdialoghdlrIsBufferEmpty(dialoghdlr) )
   {
      int len;

      /* clear the buffer */
      SCIPdialoghdlrClearBuffer(dialoghdlr);

      if( prompt == NULL )
      {
         /* use current dialog's path as prompt */
         SCIPdialogGetPath(dialog, '/', path);
         (void) SCIPsnprintf(p, SCIP_MAXSTRLEN, "%s> ", path);
         prompt = p;
      }

      /* read command line from stdin or from the input line list */
      SCIP_CALL( readInputLine(dialoghdlr, prompt, endoffile) );

      /* strip trailing spaces */
      len = (int)strlen(&dialoghdlr->buffer[dialoghdlr->bufferpos]);
      if( len > 0 )
      {
         while( isspace((unsigned char)dialoghdlr->buffer[dialoghdlr->bufferpos + len - 1]) )
         {
            dialoghdlr->buffer[dialoghdlr->bufferpos + len - 1] = '\0';
            len--;
         }
      }

      /* insert command in command history */
      if( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' )
      {
         SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, NULL, &dialoghdlr->buffer[dialoghdlr->bufferpos], FALSE) );
      }
   }

   /* the last character in the buffer must be a '\0' */
   dialoghdlr->buffer[dialoghdlr->buffersize-1] = '\0';


   /* skip leading spaces: find start of first word */
   while( isspace((unsigned char)dialoghdlr->buffer[dialoghdlr->bufferpos]) )
      dialoghdlr->bufferpos++;

   /* copy the complete line */
   *inputline = &dialoghdlr->buffer[dialoghdlr->bufferpos];

   /* go to the end of the line */
   dialoghdlr->bufferpos += (int)strlen(&dialoghdlr->buffer[dialoghdlr->bufferpos]);

   if( dialoghdlr->buffer[dialoghdlr->buffersize-1] == '\0' )
      *endoffile = TRUE;

   return SCIP_OKAY;
}


/** returns the next word in the handler's command buffer; if the buffer is empty, displays the given prompt or the
 *  current dialog's path and asks the user for further input; the user must not free or modify the returned string
 */
SCIP_RETCODE SCIPdialoghdlrGetWord(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG*          dialog,             /**< current dialog */
   const char*           prompt,             /**< prompt to display, or NULL to display the current dialog's path */
   char**                inputword,          /**< pointer to store the next word in the handler's command buffer */
   SCIP_Bool*            endoffile           /**< pointer to store whether the end of the input file was reached */
   )
{
   char path[SCIP_MAXSTRLEN];
   char p[SCIP_MAXSTRLEN];
   char* firstword;
   int pos;

   assert(dialoghdlr != NULL);
   assert(dialoghdlr->buffer != NULL);
   assert(dialoghdlr->bufferpos < dialoghdlr->buffersize);
   assert(inputword != NULL);
   assert(endoffile != NULL);

   *endoffile = FALSE;

   /* get input from the user, if the buffer is empty */
   if( SCIPdialoghdlrIsBufferEmpty(dialoghdlr) )
   {
      int len;

      /* clear the buffer */
      SCIPdialoghdlrClearBuffer(dialoghdlr);

      if( prompt == NULL )
      {
         /* use current dialog's path as prompt */
         SCIPdialogGetPath(dialog, '/', path);
         (void) SCIPsnprintf(p, SCIP_MAXSTRLEN, "%s> ", path);
         prompt = p;
      }

      /* read command line from stdin or from the input line list */
      SCIP_CALL( readInputLine(dialoghdlr, prompt, endoffile) );

      /* strip trailing spaces */
      len = (int)strlen(&dialoghdlr->buffer[dialoghdlr->bufferpos]);
      if( len > 0 )
      {
         while( isspace((unsigned char)dialoghdlr->buffer[dialoghdlr->bufferpos + len - 1]) )
         {
            dialoghdlr->buffer[dialoghdlr->bufferpos + len - 1] = '\0';
            len--;
         }
      }

      /* insert command in command history */
      if( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' )
      {
         SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, NULL, &dialoghdlr->buffer[dialoghdlr->bufferpos], FALSE) );
      }
   }

   /* the last character in the buffer must be a '\0' */
   dialoghdlr->buffer[dialoghdlr->buffersize-1] = '\0';

   /* skip leading spaces: find start of first word */
   while( isspace((unsigned char)dialoghdlr->buffer[dialoghdlr->bufferpos]) )
      dialoghdlr->bufferpos++;
   firstword = &dialoghdlr->buffer[dialoghdlr->bufferpos];

   pos = dialoghdlr->bufferpos;
   while( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' && !isspace((unsigned char)dialoghdlr->buffer[dialoghdlr->bufferpos]) )
   {
      assert(pos <= dialoghdlr->bufferpos);

      switch( dialoghdlr->buffer[dialoghdlr->bufferpos] )
      {
      case '"':
         dialoghdlr->bufferpos++;
         /* read characters as they are until the next " */
         while( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' && dialoghdlr->buffer[dialoghdlr->bufferpos] != '"' )
         {
            /* watch out for \" and \\ which should be treated as " and \, respectively */
            if( dialoghdlr->buffer[dialoghdlr->bufferpos] == '\\'
               && (dialoghdlr->buffer[dialoghdlr->bufferpos+1] == '"'
                  || dialoghdlr->buffer[dialoghdlr->bufferpos+1] == '\\') )
            {
               dialoghdlr->bufferpos++;
            }
            dialoghdlr->buffer[pos] = dialoghdlr->buffer[dialoghdlr->bufferpos];
            pos++;
            dialoghdlr->bufferpos++;
         }
         if( dialoghdlr->buffer[dialoghdlr->bufferpos] == '"' )
            dialoghdlr->bufferpos++; /* skip final " */
         break;
      case '\'':
         dialoghdlr->bufferpos++;
         /* read characters as they are until the next ' */
         while( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' && dialoghdlr->buffer[dialoghdlr->bufferpos] != '\'' )
         {
            /* watch out for \' and \\ which should be treated as ' and \, respectively */
            if( dialoghdlr->buffer[dialoghdlr->bufferpos] == '\\'
               && (dialoghdlr->buffer[dialoghdlr->bufferpos+1] == '\''
                  || dialoghdlr->buffer[dialoghdlr->bufferpos+1] == '\\') )
            {
               dialoghdlr->bufferpos++;
            }
            dialoghdlr->buffer[pos] = dialoghdlr->buffer[dialoghdlr->bufferpos];
            pos++;
            dialoghdlr->bufferpos++;
         }
         if( dialoghdlr->buffer[dialoghdlr->bufferpos] == '\'' )
            dialoghdlr->bufferpos++; /* skip final ' */
         break;
      case '\\':
         /* if the next character is a space, a ", or a ', read next character as it is;
          * otherwise, treat the \ as normal character
          */
         if( dialoghdlr->buffer[dialoghdlr->bufferpos+1] == ' '
            || dialoghdlr->buffer[dialoghdlr->bufferpos+1] == '"'
            || dialoghdlr->buffer[dialoghdlr->bufferpos+1] == '\'' )
         {
            dialoghdlr->bufferpos++;
         }
         /*lint -fallthrough*/
      default:
         dialoghdlr->buffer[pos] = dialoghdlr->buffer[dialoghdlr->bufferpos];
         pos++;
         dialoghdlr->bufferpos++;
         break;
      }
   }
   assert(pos <= dialoghdlr->bufferpos);

   /* move buffer to the next position */
   if( dialoghdlr->buffer[dialoghdlr->bufferpos] != '\0' )
      dialoghdlr->bufferpos++;

   /* truncate the command word in the buffer */
   if( dialoghdlr->buffer[pos] != '\0' )
      dialoghdlr->buffer[pos] = '\0';

   /* remove additional spaces */
   while( isspace((unsigned char)dialoghdlr->buffer[dialoghdlr->bufferpos]) )
      dialoghdlr->bufferpos++;

   *inputword = firstword;

   SCIPdebugMessage("next word: <%s>\n", *inputword);

   return SCIP_OKAY;
}

/** adds a single line of input to the dialog handler which is treated as if the user entered the command line */
SCIP_RETCODE SCIPdialoghdlrAddInputLine(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   const char*           inputline           /**< input line to add */
   )
{
   SCIP_LINELIST* linelist;
   SCIP_RETCODE retcode = SCIP_OKAY;

   assert(dialoghdlr != NULL);
   assert(dialoghdlr->inputlistptr != NULL);
   assert(*dialoghdlr->inputlistptr == NULL);
   assert(inputline != NULL);

   SCIP_ALLOC( BMSallocMemory(&linelist) );
   SCIP_ALLOC_TERMINATE( retcode, BMSduplicateMemoryArray(&linelist->inputline, inputline, strlen(inputline)+1), TERMINATE );
   linelist->nextline = NULL;
   *dialoghdlr->inputlistptr = linelist;
   dialoghdlr->inputlistptr = &linelist->nextline;

 TERMINATE:
   if( retcode != SCIP_OKAY )
      BMSfreeMemory(&linelist);

   return retcode;
}

/** adds a command to the command history of the dialog handler; if a dialog is given, the command is preceeded
 *  by the dialog's command path; if no command is given, only the path to the dialog is added to the command history
 */
SCIP_RETCODE SCIPdialoghdlrAddHistory(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG*          dialog,             /**< current dialog, or NULL */
   const char*           command,            /**< command string to add to the command history, or NULL */
   SCIP_Bool             escapecommand       /**< should special characters in command be prefixed by an escape char? */
   )
{
   char s[SCIP_MAXSTRLEN];
   char h[SCIP_MAXSTRLEN];
   SCIP_Bool cleanuphistory;

   assert(dialoghdlr != NULL);

   /* the current history list should be cleaned up if a dialog is given (i.e. the command is not partial) */
   cleanuphistory = (dialog != NULL);

   /* generate the string to add to the history */
   s[SCIP_MAXSTRLEN-1] = '\0';
   h[SCIP_MAXSTRLEN-1] = '\0';

   if( command != NULL )
   {
      if( escapecommand )
         SCIPescapeString(h, SCIP_MAXSTRLEN, command);
      else
         (void)strncpy(h, command, SCIP_MAXSTRLEN-1);
   }
   else
      h[0] = '\0';

   while( dialog != NULL && dialog != dialoghdlr->rootdialog )
   {
      if( h[0] == '\0' )
         (void)strncpy(h, dialog->name, SCIP_MAXSTRLEN-1);
      else
      {
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%s %s", dialog->name, h);
         (void)strncpy(h, s, SCIP_MAXSTRLEN-1);
      }
      dialog = dialog->parent;
   }

   /* clean up the unmarked history entries */
   if( cleanuphistory )
   {
      int i;

      for( i = getHistoryLength()-1; i >= dialoghdlr->nprotectedhistelems; --i )
      {
         SCIP_CALL( removeHistory(i) );
      }
   }

   /* add command to history */
   if( h[0] != '\0' )
   {
      SCIP_CALL( addHistory(h) );
   }

   /* if the history string was a full command line, protect the history entry from future cleanups */
   if( cleanuphistory )
   {
      dialoghdlr->nprotectedhistelems = getHistoryLength();
   }

   return SCIP_OKAY;
}




/*
 * dialog
 */

/** ensures, that sub-dialogs array can store at least the given number of sub-dialogs */
static
SCIP_RETCODE ensureSubdialogMem(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal storage size for sub-dialogs */
   )
{
   assert(dialog != NULL);

   if( num > dialog->subdialogssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&(dialog->subdialogs), newsize) );
      dialog->subdialogssize = newsize;
   }
   assert(num <= dialog->subdialogssize);

   return SCIP_OKAY;
}

/** creates and captures a user interface dialog */
SCIP_RETCODE SCIPdialogCreate(
   SCIP_DIALOG**         dialog,             /**< pointer to store the dialog */
   SCIP_DECL_DIALOGCOPY  ((*dialogcopy)),    /**< copy method of dialog or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   SCIP_DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   SCIP_DECL_DIALOGFREE  ((*dialogfree)),    /**< destructor of dialog to free user data, or NULL */
   const char*           name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*           desc,               /**< description of dialog used if description output method is NULL */
   SCIP_Bool             issubmenu,          /**< is the dialog a sub-menu? */
   SCIP_DIALOGDATA*      dialogdata          /**< user defined dialog data */
   )
{
   SCIP_RETCODE retcode;

   assert(dialog != NULL);
   assert(name != NULL);

   retcode = SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemory(dialog) );
   (*dialog)->dialogcopy = dialogcopy;
   (*dialog)->dialogexec = dialogexec;
   (*dialog)->dialogdesc = dialogdesc;
   (*dialog)->dialogfree = dialogfree;

   SCIP_ALLOC_TERMINATE( retcode, BMSduplicateMemoryArray(&(*dialog)->name, name, strlen(name)+1), TERMINATE );
   if( desc != NULL )
   {
      SCIP_ALLOC_TERMINATE( retcode, BMSduplicateMemoryArray(&(*dialog)->desc, desc, strlen(desc)+1), TERMINATE );
   }
   else
      (*dialog)->desc = NULL;

   (*dialog)->issubmenu = issubmenu;
   (*dialog)->parent = NULL;
   (*dialog)->subdialogs = NULL;
   (*dialog)->nsubdialogs = 0;
   (*dialog)->subdialogssize = 0;
   (*dialog)->nuses = 0;
   (*dialog)->dialogdata = dialogdata;

   /* capture dialog */
   SCIPdialogCapture(*dialog);

 TERMINATE:
   if( retcode != SCIP_OKAY )
   {
      BMSfreeMemoryArrayNull(&(*dialog)->name);
      BMSfreeMemory(dialog);
   }

   return retcode;
}

/** frees dialog and all of its sub-dialogs */
static
SCIP_RETCODE dialogFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog              /**< pointer to dialog */
   )
{
   int i;

   assert(dialog != NULL);
   assert(*dialog != NULL);
   assert((*dialog)->nuses == 0);

   /* call destructor of dialog */
   if( (*dialog)->dialogfree != NULL )
   {
      SCIP_CALL( (*dialog)->dialogfree(scip, *dialog) );
   }

   /* release sub-dialogs */
   for( i = 0; i < (*dialog)->nsubdialogs; ++i )
   {
      SCIP_CALL( SCIPdialogRelease(scip, &(*dialog)->subdialogs[i]) );
   }
   BMSfreeMemoryArrayNull(&(*dialog)->subdialogs);

   BMSfreeMemoryArrayNull(&(*dialog)->name);
   BMSfreeMemoryArrayNull(&(*dialog)->desc);
   BMSfreeMemory(dialog);

   return SCIP_OKAY;
}

/** captures a dialog */
void SCIPdialogCapture(
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   dialog->nuses++;
}

/** releases a dialog */
SCIP_RETCODE SCIPdialogRelease(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog              /**< pointer to dialog */
   )
{
   assert(dialog != NULL);

   (*dialog)->nuses--;
   if( (*dialog)->nuses == 0 )
   {
      SCIP_CALL( dialogFree(scip, dialog) );
   }

   return SCIP_OKAY;
}

/** executes dialog */
SCIP_RETCODE SCIPdialogExec(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG**         nextdialog          /**< pointer to store the next dialog to process */
   )
{
   assert(dialog != NULL);
   assert(dialog->dialogexec != NULL);
   assert(set != NULL);
   assert(nextdialog != NULL);

   SCIP_CALL( dialog->dialogexec(set->scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}

/** comparison method for sorting dialogs w.r.t. to their name */
static
SCIP_DECL_SORTPTRCOMP(dialogComp)
{
   return strcmp( SCIPdialogGetName((SCIP_DIALOG*)elem1), SCIPdialogGetName((SCIP_DIALOG*)elem2) );
}

/** adds a sub-dialog to the given dialog as menu entry and captures the sub-dialog */
SCIP_RETCODE SCIPdialogAddEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DIALOG*          subdialog           /**< sub-dialog to add as menu entry in dialog */
   )
{
   assert(dialog != NULL);
   assert(subdialog != NULL);

   /* check, if sub-dialog already exists */
   if( SCIPdialogHasEntry(dialog, SCIPdialogGetName(subdialog)) )
   {
      SCIPerrorMessage("dialog entry with name <%s> already exists in dialog <%s>\n",
         SCIPdialogGetName(subdialog), SCIPdialogGetName(dialog));
      return SCIP_INVALIDDATA;
   }

   /* resize the sub-dialogs array */
   SCIP_CALL( ensureSubdialogMem(dialog, set, dialog->nsubdialogs+1) );

   /* link the dialogs as parent-child pair; the sub-dialogs are sorted non-decreasing w.r.t. their name */
   SCIPsortedvecInsertPtr((void**)dialog->subdialogs, dialogComp, (void*)subdialog, &dialog->nsubdialogs, NULL);
   subdialog->parent = dialog;

   /* capture sub-dialog */
   SCIPdialogCapture(subdialog);

   return SCIP_OKAY;
}

/** returns TRUE iff a dialog entry matching exactly the given name is existing in the given dialog */
SCIP_Bool SCIPdialogHasEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   const char*           entryname           /**< name of the dialog entry to find */
   )
{
   SCIP_DIALOG** subdialogs;
   int nsubdialogs;
   int i;

   assert(dialog != NULL);
   assert(entryname != NULL);

   /* check entryname w.r.t. available dialog options */
   subdialogs = SCIPdialogGetSubdialogs(dialog);
   nsubdialogs = SCIPdialogGetNSubdialogs(dialog);
   for( i = 0; i < nsubdialogs; ++i )
   {
      /* check, if the sub-dialog's name matches entryname */
      if( strcmp(entryname, SCIPdialogGetName(subdialogs[i])) == 0 )
         return TRUE;
   }

   return FALSE;
}

/** searches the dialog for entries corresponding to the given name;
 *  If a complete match is found, the entry is returned as "subdialog" and
 *  the return value is 1.
 *  If no dialog entry completely matches the given "entryname", the number
 *  of entries with names beginning with "entryname" is returned. If this
 *  number is 1, the single match is returned as "subdialog". Otherwise,
 *  "subdialog" is set to NULL.
 */
int SCIPdialogFindEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   const char*           entryname,          /**< name of the dialog entry to find */
   SCIP_DIALOG**         subdialog           /**< pointer to store the found dialog entry */
   )
{
   SCIP_DIALOG** subdialogs;
   unsigned int namelen;
   int nsubdialogs;
   int nfound;
   int i;

   assert(dialog != NULL);
   assert(entryname != NULL);
   assert(subdialog != NULL);

   *subdialog = NULL;

   /* check entryname w.r.t. available dialog options */
   subdialogs = SCIPdialogGetSubdialogs(dialog);
   nsubdialogs = SCIPdialogGetNSubdialogs(dialog);
   namelen = (unsigned int) strlen(entryname);
   nfound = 0;
   for( i = 0; i < nsubdialogs; ++i )
   {
      /* check, if the beginning of the sub-dialog's name matches entryname */
      if( strncmp(entryname, SCIPdialogGetName(subdialogs[i]), namelen) == 0 )
      {
         *subdialog = subdialogs[i];
         nfound++;

         /* if entryname exactly matches the sub-dialog's name, use this sub-dialog */
         if( namelen == (unsigned int) strlen(SCIPdialogGetName(subdialogs[i])) )
            return 1;
      }
   }

   if( nfound != 1 )
      *subdialog = NULL;

   return nfound;
}

/** displays the dialog's menu */
SCIP_RETCODE SCIPdialogDisplayMenu(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip                /**< SCIP data structure */   
   )
{
   int i;

   assert(dialog != NULL);

   /* display the dialog's sub menus */
   for( i = 0; i < dialog->nsubdialogs; ++i )
   {
      if( SCIPdialogIsSubmenu(dialog->subdialogs[i]) )
      {
         SCIP_CALL( SCIPdialogDisplayMenuEntry(dialog->subdialogs[i], scip) );
      }
   }

   /* display the dialog's menu options */
   for( i = 0; i < dialog->nsubdialogs; ++i )
   {
      if( !SCIPdialogIsSubmenu(dialog->subdialogs[i]) )
      {
         SCIP_CALL( SCIPdialogDisplayMenuEntry(dialog->subdialogs[i], scip) );
      }
   }

   if( dialog->nsubdialogs == 0 )
      SCIPdialogMessage(scip, NULL, "<no options available>\n");

   return SCIP_OKAY;
}

/** displays the entry for the dialog in it's parent's menu */
SCIP_RETCODE SCIPdialogDisplayMenuEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip                /**< SCIP data structure */   
   )
{
   char name[SCIP_MAXSTRLEN];

   assert(dialog != NULL);

   /* display the dialog's name */
   if( dialog->issubmenu )
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "<%s>", dialog->name);
   else
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", dialog->name);
   SCIPdialogMessage(scip, NULL, "  %-21s ", name);
   if( strlen(name) > 21 )
   {
      /* break the line, and start the description in the next line */
      SCIPdialogMessage(scip, NULL, "\n                   -->  ");
   }

   /* display the dialog's description */
   if( dialog->dialogdesc != NULL )
   {
      SCIP_CALL( dialog->dialogdesc(scip, dialog) );
   }
   else
      SCIPdialogMessage(scip, NULL, "%s",dialog->desc);
   SCIPdialogMessage(scip, NULL, "\n");

   return SCIP_OKAY;
}

/** displays all dialog entries with names starting with the given "entryname" */
SCIP_RETCODE SCIPdialogDisplayCompletions(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip,               /**< SCIP data structure */   
   const char*           entryname           /**< name of the dialog entry to find */
   )
{
   SCIP_DIALOG** subdialogs;
   unsigned int namelen;
   int nsubdialogs;
   int i;

   assert(dialog != NULL);
   assert(entryname != NULL);

   /* check entryname w.r.t. available dialog options */
   subdialogs = SCIPdialogGetSubdialogs(dialog);
   nsubdialogs = SCIPdialogGetNSubdialogs(dialog);
   namelen = (unsigned int) strlen(entryname);
   for( i = 0; i < nsubdialogs; ++i )
   {
      /* check, if the beginning of the sub-dialog's name matches entryname */
      if( strncmp(entryname, SCIPdialogGetName(subdialogs[i]), namelen) == 0 )
      {
         SCIP_CALL( SCIPdialogDisplayMenuEntry(subdialogs[i], scip) );
      }
   }

   return SCIP_OKAY;
}

/** gets the name of the current path in the dialog tree, separated by the given character */
void SCIPdialogGetPath(
   SCIP_DIALOG*          dialog,             /**< dialog */
   const char            sepchar,            /**< separation character to insert in path */
   char*                 path                /**< string buffer to store the path */
   )
{
   char s[SCIP_MAXSTRLEN];

   assert(dialog != NULL);

   (void)strncpy(path, dialog->name, SCIP_MAXSTRLEN);
   path[SCIP_MAXSTRLEN - 1] = '\0';

   dialog = dialog->parent;
   while( dialog != NULL )
   {
      (void)SCIPsnprintf(s, SCIP_MAXSTRLEN, "%s%c%s", dialog->name, sepchar, path);
      (void)strncpy(path, s, SCIP_MAXSTRLEN);
      path[SCIP_MAXSTRLEN - 1] = '\0';
      dialog = dialog->parent;
   }
}

/** gets the command name of the dialog */
const char* SCIPdialogGetName(
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->name;
}

/** gets the description of the dialog */
const char* SCIPdialogGetDesc(
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->desc;
}

/** returns whether the dialog is a sub menu */
SCIP_Bool SCIPdialogIsSubmenu(
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->issubmenu;
}

/** gets the parent dialog of the given dialog */
SCIP_DIALOG* SCIPdialogGetParent(
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->parent;
}

/** gets the array of sub-dialogs associated with the given dialog */
SCIP_DIALOG** SCIPdialogGetSubdialogs(
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->subdialogs;
}

/** gets the number of sub-dialogs associated with the given dialog */
int SCIPdialogGetNSubdialogs(
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->nsubdialogs;
}

/** gets the user defined data associated with the given dialog */
SCIP_DIALOGDATA* SCIPdialogGetData(
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(dialog != NULL);

   return dialog->dialogdata;
}

/** sets user data of dialog; user has to free old data in advance! */
void SCIPdialogSetData(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP_DIALOGDATA*      dialogdata          /**< new dialog user data */
   )
{
   assert(dialog != NULL);

   dialog->dialogdata = dialogdata;
}

/** writes command history to specified filename */
SCIP_RETCODE SCIPdialogWriteHistory(
   const char*           filename            /**< file name for (over)writing history */
   )
{
   return writeHistory(filename);
}
