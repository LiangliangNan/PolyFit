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
SCIP_EXPORT
SCIP_DIALOG* SCIPdialoghdlrGetRoot(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** clears the input command buffer of the dialog handler */
SCIP_EXPORT
void SCIPdialoghdlrClearBuffer(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** returns TRUE iff input command buffer is empty */
SCIP_EXPORT
SCIP_Bool SCIPdialoghdlrIsBufferEmpty(
   SCIP_DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** returns the next line in the handler's command buffer; if the buffer is empty, displays the given prompt or the
 *  current dialog's path and asks the user for further input; the user must not free or modify the returned string
 */
SCIP_EXPORT
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
SCIP_EXPORT
SCIP_RETCODE SCIPdialoghdlrGetWord(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG*          dialog,             /**< current dialog */
   const char*           prompt,             /**< prompt to display, or NULL to display the current dialog's path */
   char**                inputword,          /**< pointer to store the next word in the handler's command buffer */
   SCIP_Bool*            endoffile           /**< pointer to store whether the end of the input file was reached */
   );

/** adds a single line of input to the dialog handler which is treated as if the user entered the command line */
SCIP_EXPORT
SCIP_RETCODE SCIPdialoghdlrAddInputLine(
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   const char*           inputline           /**< input line to add */
   );

/** adds a command to the command history of the dialog handler; if a dialog is given, the command is preceeded
 *  by the dialog's command path; if no command is given, only the path to the dialog is added to the command history
 */
SCIP_EXPORT
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
SCIP_EXPORT
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
SCIP_EXPORT
int SCIPdialogFindEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   const char*           entryname,          /**< name of the dialog entry to find */
   SCIP_DIALOG**         subdialog           /**< pointer to store the found dialog entry */
   );

/** displays the dialog's menu */
SCIP_EXPORT
SCIP_RETCODE SCIPdialogDisplayMenu(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip                /**< SCIP data structure */   
   );

/** displays the entry for the dialog in it's parent's menu */
SCIP_EXPORT
SCIP_RETCODE SCIPdialogDisplayMenuEntry(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip                /**< SCIP data structure */   
   );

/** displays all dialog entries with names starting with the given "entryname" */
SCIP_EXPORT
SCIP_RETCODE SCIPdialogDisplayCompletions(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP*                 scip,               /**< SCIP data structure */   
   const char*           entryname           /**< name of the dialog entry to find */
   );

/** gets the name of the current path in the dialog tree, separated by the given character */
SCIP_EXPORT
void SCIPdialogGetPath(
   SCIP_DIALOG*          dialog,             /**< dialog */
   const char            sepchar,            /**< separation character to insert in path */
   char*                 path                /**< string buffer to store the path */
   );

/** gets the command name of the dialog */
SCIP_EXPORT
const char* SCIPdialogGetName(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the description of the dialog */
SCIP_EXPORT
const char* SCIPdialogGetDesc(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** returns whether the dialog is a sub menu */
SCIP_EXPORT
SCIP_Bool SCIPdialogIsSubmenu(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the parent dialog of the given dialog */
SCIP_EXPORT
SCIP_DIALOG* SCIPdialogGetParent(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the array of sub-dialogs associated with the given dialog */
SCIP_EXPORT
SCIP_DIALOG** SCIPdialogGetSubdialogs(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the number of sub-dialogs associated with the given dialog */
SCIP_EXPORT
int SCIPdialogGetNSubdialogs(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** gets the user defined data associated with the given dialog */
SCIP_EXPORT
SCIP_DIALOGDATA* SCIPdialogGetData(
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** sets user data of dialog; user has to free old data in advance! */
SCIP_EXPORT
void SCIPdialogSetData(
   SCIP_DIALOG*          dialog,             /**< dialog */
   SCIP_DIALOGDATA*      dialogdata          /**< new dialog user data */
   );

/** writes command history to specified filename */
SCIP_EXPORT
SCIP_RETCODE SCIPdialogWriteHistory(
   const char*           filename            /**< file name for (over)writing history */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
