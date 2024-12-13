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

/**@file   struct_disp.h
 * @ingroup INTERNALAPI
 * @brief  data structures for displaying runtime statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_DISP_H__
#define __SCIP_STRUCT_DISP_H__


#include "scip/def.h"
#include "scip/type_disp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** display column */
struct SCIP_Disp
{
   char*                 name;               /**< name of display column */
   char*                 desc;               /**< description of display column */
   char*                 header;             /**< head line of display column */
   SCIP_DECL_DISPCOPY    ((*dispcopy));      /**< copy method of display column or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DISPFREE    ((*dispfree));      /**< destructor of display column */
   SCIP_DECL_DISPINIT    ((*dispinit));      /**< initialize display column */
   SCIP_DECL_DISPEXIT    ((*dispexit));      /**< deinitialize display column */
   SCIP_DECL_DISPINITSOL ((*dispinitsol));   /**< solving process initialization method of display column */
   SCIP_DECL_DISPEXITSOL ((*dispexitsol));   /**< solving process deinitialization method of display column */
   SCIP_DECL_DISPOUTPUT  ((*dispoutput));    /**< output method */
   SCIP_DISPDATA*        dispdata;           /**< display column data */
   int                   width;              /**< width of display column (no. of chars used) */
   int                   priority;           /**< priority of display column */
   int                   position;           /**< relative position of display column */
   SCIP_DISPSTATUS       dispstatus;         /**< display activation status of display column */
   SCIP_Bool             stripline;          /**< should the column be separated with a line from its right neighbor? */
   SCIP_Bool             initialized;        /**< is display column initialized? */
   SCIP_Bool             active;             /**< should column be displayed to the screen? */
   SCIP_DISPMODE         mode;               /**< mode in which the display column is activated */
};

#ifdef __cplusplus
}
#endif

#endif
