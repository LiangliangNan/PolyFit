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

/**@file   type_dialog.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for user interface dialog
 * @author Tobias Achterberg
 *
 *  This file defines the interface for dialogs implemented in C.
 *
 * - \ref DIALOG "Instructions for implementing a dialog"
 * - \ref DIALOGS "List of available dialogs"
 * - \ref scip::ObjDialog "C++ wrapper class
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_DIALOG_H__
#define __SCIP_TYPE_DIALOG_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Dialog SCIP_DIALOG;           /**< user interface dialog */
typedef struct SCIP_DialogData SCIP_DIALOGDATA;   /**< user defined dialog data */
typedef struct SCIP_Dialoghdlr SCIP_DIALOGHDLR;   /**< dialog handler */
typedef struct SCIP_Linelist SCIP_LINELIST;       /**< linked list of single input lines */


/** copy method for dialog plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - dialog          : the dialog itself
 */
#define SCIP_DECL_DIALOGCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIALOG* dialog)

/** destructor of dialog to free user data (called when the dialog is not captured anymore)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - dialog          : the dialog itself
 */
#define SCIP_DECL_DIALOGFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIALOG* dialog)

/** description output method of dialog
 *
 *  This method should output (usually a single line of) information describing the meaning of the dialog.
 *  The method is called, when the help menu of the parent's dialog is displayed.
 *  If no description output method is given, the description string of the dialog is displayed instead.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - *dialog         : the dialog itself
 */
#define SCIP_DECL_DIALOGDESC(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIALOG* dialog)

/** execution method of dialog
 *
 *  This method is invoked, if the user selected the dialog's command name in the parent's dialog menu.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - dialoghdlr      : dialog handler to call for user interaction
 *  - dialog          : the dialog itself
 *
 *  output:
 *  - *nextdialog     : next dialog to process (or NULL to quit dialog processing)
 */
#define SCIP_DECL_DIALOGEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIALOG* dialog, SCIP_DIALOGHDLR* dialoghdlr, SCIP_DIALOG** nextdialog)

#ifdef __cplusplus
}
#endif

#endif
