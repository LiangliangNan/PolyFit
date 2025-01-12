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

/**@file   struct_cutsel.h
 * @ingroup INTERNALAPI
 * @brief  data structures for cut selectors
 * @author Felipe Serrano
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CUTSEL_H__
#define __SCIP_STRUCT_CUTSEL_H__


#include "scip/def.h"
#include "scip/type_cutsel.h"

#ifdef __cplusplus
extern "C" {
#endif

/** cut selector */
struct SCIP_Cutsel
{
   char*                 name;               /**< name of cut selector */
   char*                 desc;               /**< description of cut selector */
   SCIP_DECL_CUTSELCOPY ((*cutselcopy));     /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CUTSELFREE ((*cutselfree));     /**< destructor of cut selector */
   SCIP_DECL_CUTSELINIT ((*cutselinit));     /**< initialize cut selector */
   SCIP_DECL_CUTSELEXIT ((*cutselexit));     /**< deinitialize cut selector */
   SCIP_DECL_CUTSELINITSOL((*cutselinitsol));/**< solving process initialization method of cut selector */
   SCIP_DECL_CUTSELEXITSOL((*cutselexitsol));/**< solving process deinitialization method of cut selector */
   SCIP_DECL_CUTSELSELECT((*cutselselect));  /**< cut selection method */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this cut selector for the next stages */
   SCIP_CLOCK*           cutseltime;         /**< cut selector execution time */
   SCIP_CUTSELDATA*      cutseldata;         /**< cut selector data */
   int                   priority;           /**< priority of the cut selector */
   SCIP_Bool             initialized;        /**< is cut selector initialized? */
   SCIP_Longint          ncalls;             /**< number of times, this cutselector was called */
   SCIP_Longint          nrootcalls;         /**< number of times, this cutselector was called */
   SCIP_Longint          nrootcutsselected;  /**< number of cuts selected at the root */
   SCIP_Longint          nrootcutsforced;    /**< number of forced cuts at the root */
   SCIP_Longint          nrootcutsfiltered;  /**< number of cuts filtered at the root */
   SCIP_Longint          nlocalcutsselected; /**< number of local cuts selected */
   SCIP_Longint          nlocalcutsforced;   /**< number of forced local cuts */
   SCIP_Longint          nlocalcutsfiltered; /**< number of local cuts filtered */
};

#ifdef __cplusplus
}
#endif

#endif
