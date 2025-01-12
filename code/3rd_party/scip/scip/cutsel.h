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

/**@file   cutsel.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for cut selectors
 * @author Felipe Serrano
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTSEL_H__
#define __SCIP_CUTSEL_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/pub_cutsel.h"
#include "scip/lp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a cut selector */
SCIP_RETCODE SCIPcutselCreate(
   SCIP_CUTSEL**         cutsel,             /**< pointer to store cut selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of cut selector */
   const char*           desc,               /**< description of cut selector */
   int                   priority,           /**< priority of the cut selector in standard mode */
   SCIP_DECL_CUTSELCOPY ((*cutselcopy)),     /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CUTSELFREE ((*cutselfree)),     /**< destructor of cut selector */
   SCIP_DECL_CUTSELINIT ((*cutselinit)),     /**< initialize cut selector */
   SCIP_DECL_CUTSELEXIT ((*cutselexit)),     /**< deinitialize cut selector */
   SCIP_DECL_CUTSELINITSOL((*cutselinitsol)),/**< solving process initialization method of cut selector */
   SCIP_DECL_CUTSELEXITSOL((*cutselexitsol)),/**< solving process deinitialization method of cut selector */
   SCIP_DECL_CUTSELSELECT((*cutselselect)),  /**< cut selection method */
   SCIP_CUTSELDATA*      cutseldata          /**< cut selector data */
   );

/** enables or disables all clocks of @p cutsel, depending on the value of the flag */
void SCIPcutselEnableOrDisableClocks(
   SCIP_CUTSEL*          cutsel,             /**< the cut selector for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the cut selector be enabled? */
   );

/** calls cut selectors to select cuts */
SCIP_RETCODE SCIPcutselsSelect(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW**            cuts,               /**< array with cuts to select from */
   int                   ncuts,              /**< length of cuts */
   int                   nforcedcuts,        /**< number of forced cuts at start of given array */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_Bool             initiallp,          /**< is the separation storage currently being filled with the initial LP rows? */
   int                   maxnselectedcuts,   /**< maximum number of cuts to be selected */
   int*                  nselectedcuts       /**< pointer to return number of selected cuts */
   );

/** copies the given cut selector to a new scip */
SCIP_RETCODE SCIPcutselCopyInclude(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** sets copy method of cut selector */
void SCIPcutselSetCopy(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELCOPY  ((*cutselcopy))     /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** initializes cut selector */
SCIP_RETCODE SCIPcutselInit(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** deinitializes cut selector */
SCIP_RETCODE SCIPcutselExit(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** frees memory of cut selector */
SCIP_RETCODE SCIPcutselFree(
   SCIP_CUTSEL**         cutsel,             /**< pointer to cut selector data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs cut selector that the branch and bound process is being started */
SCIP_RETCODE SCIPcutselInitsol(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs cut selector that the branch and bound process is being started */
SCIP_RETCODE SCIPcutselExitsol(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sets destructor method of cut selector */
void SCIPcutselSetFree(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELFREE  ((*cutselfree))     /**< destructor of cut selector */
   );

/** sets initialization method of cut selector */
void SCIPcutselSetInit(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELINIT  ((*cutselinit))     /**< initialize cut selector */
   );

/** sets deinitialization method of cut selector */
void SCIPcutselSetExit(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELEXIT  ((*cutselexit))     /**< deinitialize cut selector */
   );

/** sets solving process initialization method of cut selector */
void SCIPcutselSetInitsol(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELINITSOL ((*cutselinitsol))/**< solving process initialization method of cut selector */
   );

/** sets solving process deinitialization method of cut selector */
void SCIPcutselSetExitsol(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELEXITSOL ((*cutselexitsol))/**< solving process deinitialization method of cut selector */
   );

/** sets priority of cut selector */
void SCIPcutselSetPriority(
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the cut selector */
   );

#ifdef __cplusplus
}
#endif

#endif
