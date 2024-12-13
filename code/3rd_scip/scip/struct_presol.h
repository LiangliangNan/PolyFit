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

/**@file   struct_presol.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PRESOL_H__
#define __SCIP_STRUCT_PRESOL_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_presol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** presolver */
struct SCIP_Presol
{
   char*                 name;               /**< name of presolver */
   char*                 desc;               /**< description of presolver */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy));    /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRESOLFREE  ((*presolfree));    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   SCIP_DECL_PRESOLINIT  ((*presolinit));    /**< initialization method of presolver (called after problem was transformed) */
   SCIP_DECL_PRESOLEXIT  ((*presolexit));    /**< deinitialization method of presolver (called before transformed problem is freed) */
   SCIP_DECL_PRESOLINITPRE((*presolinitpre));/**< presolving initialization method of presolver (called when presolving is about to begin) */
   SCIP_DECL_PRESOLEXITPRE((*presolexitpre));/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   SCIP_DECL_PRESOLEXEC  ((*presolexec));    /**< execution method of presolver */
   SCIP_PRESOLDATA*      presoldata;         /**< presolver data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this presolver for the next stages */
   SCIP_CLOCK*           presolclock;        /**< presolving time */
   int                   priority;           /**< priority of the presolver */
   int                   maxrounds;          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   int                   lastnfixedvars;     /**< number of variables fixed before the last call to the presolver */
   int                   lastnaggrvars;      /**< number of variables aggregated before the last call to the presolver */
   int                   lastnchgvartypes;   /**< number of variable type changes before the last call to the presolver */
   int                   lastnchgbds;        /**< number of variable bounds tightened before the last call to the presolver */
   int                   lastnaddholes;      /**< number of domain holes added before the last call to the presolver */
   int                   lastndelconss;      /**< number of deleted constraints before the last call to the presolver */
   int                   lastnaddconss;      /**< number of added constraints before the last call to the presolver */
   int                   lastnupgdconss;     /**< number of upgraded constraints before the last call to the presolver */
   int                   lastnchgcoefs;      /**< number of changed coefficients before the last call to the presolver */
   int                   lastnchgsides;      /**< number of changed left or right hand sides before the last call */
   int                   nfixedvars;         /**< total number of variables fixed by this presolver */
   int                   naggrvars;          /**< total number of variables aggregated by this presolver */
   int                   nchgvartypes;       /**< total number of variable type changes by this presolver */
   int                   nchgbds;            /**< total number of variable bounds tightened by this presolver */
   int                   naddholes;          /**< total number of domain holes added by this presolver */
   int                   ndelconss;          /**< total number of deleted constraints by this presolver */
   int                   naddconss;          /**< total number of added constraints by this presolver */
   int                   nupgdconss;         /**< total number of upgraded constraints by this presolver */
   int                   nchgcoefs;          /**< total number of changed coefficients by this presolver */
   int                   nchgsides;          /**< total number of changed left or right hand sides by this presolver */
   int                   ncalls;             /**< number of times the presolver was called and tried to find reductions */
   SCIP_Bool             initialized;        /**< is presolver initialized? */
   SCIP_PRESOLTIMING     timing;             /**< timing of the presolver */
};

#ifdef __cplusplus
}
#endif

#endif
