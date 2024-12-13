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

/**@file   struct_relax.h
 * @ingroup INTERNALAPI
 * @brief  data structures for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_RELAX_H__
#define __SCIP_STRUCT_RELAX_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_relax.h"

#ifdef __cplusplus
extern "C" {
#endif

/** relaxators data */
struct SCIP_Relax
{
   SCIP_Longint          ncalls;             /**< number of times, this relaxator was called */
   SCIP_Longint          ncutoffs;           /**< number of times, this relaxator cut off a node */
   SCIP_Longint          nimprbounds;        /**< number of times, this relaxator improves the node's lower bound (w.r.t. SCIPisRelGT) */
   SCIP_Longint          naddedconss;        /**< number of times, this relaxator added constraints */
   SCIP_Longint          nreduceddom;        /**< number of times, this relaxator reduced variable domains */
   SCIP_Longint          nseparated;         /**< number of times, this relaxator separated cuts */
   SCIP_Longint          lastsolvednode;     /**< last total nodes counter, where the current relaxation was solved */
   char*                 name;               /**< name of relaxator */
   char*                 desc;               /**< description of relaxator */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy));     /**< copy method of relaxator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_RELAXFREE   ((*relaxfree));     /**< destructor of relaxator */
   SCIP_DECL_RELAXINIT   ((*relaxinit));     /**< initialize relaxator */
   SCIP_DECL_RELAXEXIT   ((*relaxexit));     /**< deinitialize relaxator */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol));  /**< solving process initialization method of relaxator */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol));  /**< solving process deinitialization method of relaxator */
   SCIP_DECL_RELAXEXEC   ((*relaxexec));     /**< execution method of relaxator */
   SCIP_RELAXDATA*       relaxdata;          /**< relaxators local data */
   SCIP_CLOCK*           setuptime;          /**< time spent for setting up this relaxator for the next stages */
   SCIP_CLOCK*           relaxclock;         /**< relaxation time */
   SCIP_Real             imprtime;           /**< time spent for executing relaxator when improving bound (or cutoff) */
   int                   priority;           /**< priority of the relaxator */
   int                   freq;               /**< frequency for calling relaxator */
   SCIP_Bool             initialized;        /**< is relaxator initialized? */
};

/** relaxation information data */
struct SCIP_Relaxation
{
   SCIP_Real             relaxsolobjval;
   SCIP_Bool             relaxsolvalid;
   SCIP_Bool             relaxsolincludeslp;
   SCIP_Bool             relaxsolzero;
   SCIP_RELAX*           lastsolrelax;       /**< relaxator of the most recent relaxation solution */
};


#ifdef __cplusplus
}
#endif

#endif
