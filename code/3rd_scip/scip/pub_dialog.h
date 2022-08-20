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

/**@file   pub_dialog.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_DIALOG_H__
#define __SCIP_PUB_DIALOG_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_dialog.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * dialog handler
 */

/**@addtogroup PublicDialogMethods
 *
 * @{
 */
/** returns the root dialog of the dialog handler */
EXTERN
SCIP_DIALOG* SCIPdialoghdlrGetRoot(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** clears the input command buffer of the dialog handler */
EXTERN
void SCIPdialoghdlrClearBuffer(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** returns TRUE iff input command buffer is empty */
EXTERN
SCIP_Bool SCIPdialoghdlrIsBufferEmpty(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** returns the next line in the handler's command buffer; if the buffer is empty, displays the given prompt or the
 *  current dialog's path and asks the user for further input; the user must not free or modify the returned string
 */
EXTERN
SCIP_RETCODE SCIPdialoghdlrGetLine(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG*          dialog,             /**< current dialog */
   const char*           prompt,             /**< prompt to display, or NULL to display the current dialog's path */
   char**                inputline,          /**< pointer to store the complete line in the handler's command buffer */
   SCIP_Bool*            endoffile           /**< pointer to store whether the end of the input file was reached */
   );

/** returns the next word in the handler's command buffer; if the buffer is empty, displays the given prompt or the 
 *  current dialog's path and asks the user for further input; the user must not free or modify the returned string
 */
EXTERN
SCIP_RETCODE SCIPdialoghdlrGetWord(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG*          dialog,             /**< current dialog */
   const char*           prompt,             /**< prompt to display, or NULL to display the current dialog's path */
   char**                inputword,          /**< pointer to store the next word in the handler's command buffer */
   SCIP_Bool*            endoffile           /**< pointer to store whether the end of the input file was reached */
   );

/** adds a single line of input to the dialog handler which is treated as if the user entered the command line */
EXTERN
SCIP_RETCODE SCIPdialoghdlrAddInputLine(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   const char*           inputline           /**< input line to add */
   );

/** adds a command to the command history of the dialog handler; if a dialog is given, the command is preceeded
 *  by the dialog's command path; if no command is given, only the path to the dialog is added to the command history
 */
EXTERN
SCIP_RETCODE SCIPdialoghdlrAddHistory(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG*          dialog,             /**< current dialog, or NULL */
   const char*           command,            /**< command string to add to the command history, or NULL */
   SCIP_Bool             escapecommand       /**< should special characters in command be prefixed by an escape char? */
   );




/*
 * dialog
 */

/** returns TRUE iff a dialog entry matching exactly the given name is existing in the given dialog */
EXTERN
SCIP_Bool SCIPdialogHasEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   const char*           entryname           /**< name of the dialog entry to find */
   );

/** searches the dialog for entries corresponding to the given name;
 *  If a complete match is found, the entry is returned as "subdialog" and
 *  the return value is 1.
 *  If no dialog entry completely matches the given "entryname", the number
 *  of entries with names beginning with "entryname" is returned. If this
 *  number is 1, the single match is returned as "subdialog". Otherwise,
 *  "subdialog" is set to NULL.
 */
EXTERN
int SCIPdialogFindEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   const char*           entryname,          /**< name of the dialog entry to find */
   SCIP_DIALOG**         subdialog           /**< pointer to store the found dialog entry */
   );

/** displays the dialog's menu */
EXTERN
SCIP_RETCODE SCIPdialogDisplayMenu(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip                /**< SCIP data structure */   
   );

/** displays the entry for the dialog in it's parent's menu */
EXTERN
SCIP_RETCODE SCIPdialogDisplayMenuEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip                /**< SCIP data structure */   
   );

/** displays all dialog entries with names starting with the given "entryname" */
EXTERN
SCIP_RETCODE SCIPdialogDisplayCompletions(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip,               /**< SCIP data structure */   
   const char*           entryname           /**< name of the dialog entry to find */
   );

/** gets the name of the current path in the dialog tree, separated by the given character */
EXTERN
void SCIPdialogGetPath(
   SCIP_DIALOG*          dialog,             /**< dialog */
   const char            sepchar,            /**< separation character to insert in path */
   char*                 path                /**< string buffer to store the path */
   );

/** gets the command name of the dialog */
EXTERN
const char* SCIPdialogGetName(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the description of the dialog */
EXTERN
const char* SCIPdialogGetDesc(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** returns whether the dialog is a sub menu */
EXTERN
SCIP_Bool SCIPdialogIsSubmenu(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the parent dialog of the given dialog */
EXTERN
SCIP_DIALOG* SCIPdialogGetParent(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the array of sub-dialogs associated with the given dialog */
EXTERN
SCIP_DIALOG** SCIPdialogGetSubdialogs(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the number of sub-dialogs associated with the given dialog */
EXTERN
int SCIPdialogGetNSubdialogs(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the user defined data associated with the given dialog */
EXTERN
SCIP_DIALOGDATA* SCIPdialogGetData(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** sets user data of dialog; user has to free old data in advance! */
EXTERN
void SCIPdialogSetData(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP_DIALOGDATA*      dialogdata          /**< new dialog user data */
   );

/** writes command history to specified filename */
EXTERN
SCIP_RETCODE SCIPdialogWriteHistory(
   const char*           filename            /**< file name for (over)writing history */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
