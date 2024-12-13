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

/**@file   struct_dialog.h
 * @ingroup INTERNALAPI
 * @brief  data structures for user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_DIALOG_H__
#define __SCIP_STRUCT_DIALOG_H__


#include "scip/def.h"
#include "scip/type_dialog.h"

#ifdef __cplusplus
extern "C" {
#endif

/** user interface dialog */
struct SCIP_Dialog
{
   SCIP_DECL_DIALOGCOPY  ((*dialogcopy));    /**< copy method of dialog or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DIALOGEXEC  ((*dialogexec));    /**< execution method of dialog */
   SCIP_DECL_DIALOGDESC  ((*dialogdesc));    /**< description output method of dialog, or NULL */
   SCIP_DECL_DIALOGFREE  ((*dialogfree));    /**< destructor of dialog to free user data, or NULL */
   char*                 name;               /**< name of dialog: command name appearing in parent's dialog menu */
   char*                 desc;               /**< description of dialog used if description output method is NULL */
   SCIP_DIALOG*          parent;             /**< parent dialog of dialog */
   SCIP_DIALOG**         subdialogs;         /**< sub dialogs of dialog */
   SCIP_DIALOGDATA*      dialogdata;         /**< user defined dialog data */
   int                   nsubdialogs;        /**< number of sub dialogs */
   int                   subdialogssize;     /**< size of subdialogs array */
   int                   nuses;              /**< number of times, the dialog is used */
   SCIP_Bool             issubmenu;          /**< is the dialog a submenu? */
};

/** linked list of single input lines */
struct SCIP_Linelist
{
   char*                 inputline;          /**< single line of input */
   SCIP_LINELIST*        nextline;           /**< next input line */
};

/** dialog handler */
struct SCIP_Dialoghdlr
{
   SCIP_DIALOG*          rootdialog;         /**< main (root) dialog */
   SCIP_LINELIST*        inputlist;          /**< list of input lines that are processed before stdin inputs */
   SCIP_LINELIST**       inputlistptr;       /**< pointer to the ending nextline pointer of the list (which points to 0) */
   char*                 buffer;             /**< command buffer */
   int                   buffersize;         /**< size of command buffer */
   int                   bufferpos;          /**< position of first unprocessed character in buffer */
   int                   nprotectedhistelems;/**< number of history entries protected from cleaning up */
};

#ifdef __cplusplus
}
#endif

#endif
