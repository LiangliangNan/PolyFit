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

/**@file   scip_disp.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for display handler plugins
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

#ifndef __SCIP_SCIP_DISP_H__
#define __SCIP_SCIP_DISP_H__


#include "scip/def.h"
#include "scip/type_disp.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicDisplayMethods
 *
 * @{
 */

/** creates a display column and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of display column */
   const char*           desc,               /**< description of display column */
   const char*           header,             /**< head line of display column */
   SCIP_DISPSTATUS       dispstatus,         /**< display activation status of display column */
   SCIP_DECL_DISPCOPY    ((*dispcopy)),      /**< copy method of display column or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   SCIP_DECL_DISPINIT    ((*dispinit)),      /**< initialize display column */
   SCIP_DECL_DISPEXIT    ((*dispexit)),      /**< deinitialize display column */
   SCIP_DECL_DISPINITSOL ((*dispinitsol)),   /**< solving process initialization method of display column */
   SCIP_DECL_DISPEXITSOL ((*dispexitsol)),   /**< solving process deinitialization method of display column */
   SCIP_DECL_DISPOUTPUT  ((*dispoutput)),    /**< output method */
   SCIP_DISPDATA*        dispdata,           /**< display column data */
   int                   width,              /**< width of display column (no. of chars used) */
   int                   priority,           /**< priority of display column */
   int                   position,           /**< relative position of display column */
   SCIP_Bool             stripline           /**< should the column be separated with a line from its right neighbor? */
   );

/** returns the display column of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_DISP* SCIPfindDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of display column */
   );

/** returns the array of currently available display columns */
SCIP_EXPORT
SCIP_DISP** SCIPgetDisps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available display columns */
SCIP_EXPORT
int SCIPgetNDisps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** automatically selects display columns for being shown w.r.t. the display width parameter */
SCIP_EXPORT
SCIP_RETCODE SCIPautoselectDisps(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** changes the display column mode */
SCIP_EXPORT
void SCIPchgDispMode(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_DISPMODE         mode                /**< the display column mode */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
