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

/**@file   scip_dialog.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for dialog handler plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_DIALOG_H__
#define __SCIP_SCIP_DIALOG_H__


#include "scip/def.h"
#include "scip/type_dialog.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicDialogMethods
 *
 * @{
 */

/** creates and includes dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog,             /**< pointer to store the dialog */
   SCIP_DECL_DIALOGCOPY  ((*dialogcopy)),    /**< copy method of dialog or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   SCIP_DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   SCIP_DECL_DIALOGFREE  ((*dialogfree)),    /**< destructor of dialog to free user data, or NULL */
   const char*           name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*           desc,               /**< description of dialog used if description output method is NULL */
   SCIP_Bool             issubmenu,          /**< is the dialog a submenu? */
   SCIP_DIALOGDATA*      dialogdata          /**< user defined dialog data */
   );

/** returns if the dialog already exists
 *
 *  @return TRUE is returned if the dialog exists, otherwise FALSE.
 */
SCIP_EXPORT
SCIP_Bool SCIPexistsDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** captures a dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcaptureDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** releases a dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPreleaseDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog              /**< pointer to the dialog */
   );

/** makes given dialog the root dialog of SCIP's interactive user shell; captures dialog and releases former root dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsetRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog to be the root */
   );

/** returns the root dialog of SCIP's interactive user shell
 *
 *  @return the root dialog of SCIP's interactive user shell is returned.
 */
SCIP_EXPORT
SCIP_DIALOG* SCIPgetRootDialog(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds a sub dialog to the given dialog as menu entry and captures it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddDialogEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog,             /**< dialog to extend, or NULL for root dialog */
   SCIP_DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   );

/** adds a single line of input which is treated as if the user entered the command line
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddDialogInputLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           inputline           /**< input line to add */
   );

/** adds a single line of input to the command history which can be accessed with the cursor keys
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddDialogHistoryLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           inputline           /**< input line to add */
   );

/** starts interactive mode of SCIP by executing the root dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the
 *        interactive shell was closed:
 *        - \ref SCIP_STAGE_PROBLEM if the interactive shell was closed after the problem was created
 *        - \ref SCIP_STAGE_TRANSFORMED if the interactive shell was closed after the problem was transformed
 *        - \ref SCIP_STAGE_PRESOLVING if the interactive shell was closed  during presolving
 *        - \ref SCIP_STAGE_PRESOLVED if the interactive shell was closed after presolve
 *        - \ref SCIP_STAGE_SOLVING if the interactive shell was closed during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the interactive shell was closed after the problem was solved
 *        - \ref SCIP_STAGE_FREE if the interactive shell was closed after the problem was freed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPstartInteraction(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
