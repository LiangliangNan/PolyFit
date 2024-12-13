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

/**@file   struct_sepastore.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for storing conflicts
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CONFLICTSTORE_H__
#define __SCIP_STRUCT_CONFLICTSTORE_H__


#include "scip/def.h"
#include "scip/type_conflictstore.h"

#ifdef __cplusplus
extern "C" {
#endif

/** storage for conflicts */
struct SCIP_ConflictStore
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler to catch improving solutions */
   SCIP_CONS**           conflicts;          /**< array with conflicts */
   SCIP_CONS**           dualrayconfs;       /**< array with proofs based on dual rays */
   SCIP_CONS**           dualsolconfs;       /**< array with proofs based on dual solutions */
   SCIP_CONS**           origconfs;          /**< array of original conflicts added in stage SCIP_STAGE_PROBLEM */
   SCIP_Real*            confprimalbnds;     /**< array of primal bounds valid at the time the corresponding bound exceeding
                                              *   conflict was found (-infinity if the conflict based on an infeasible LP) */
   SCIP_Real*            dualprimalbnds;     /**< array of primal bounds valid at the time the corresponding dual proof
                                              *   based on a dual solution was found */
   SCIP_Real*            scalefactors;       /**< scaling factor that needs to be considered when updating the side */
   SCIP_Bool*            updateside;         /**< array to store whether the side should be updated whenever a new incumbent is found */
   SCIP_Bool*            drayrelaxonly;      /**< array to store whether the dual proof is valid for the current relaxation only */
   SCIP_Bool*            dsolrelaxonly;      /**< array to store whether the dual proof is valid for the current relaxation only */
   SCIP_Real             avgswitchlength;    /**< average length of switched paths */
   SCIP_Real             lastcutoffbound;    /**< last cutoff bound for which the conflict store was cleaned */
   SCIP_Longint          lastnodenum;        /**< number of the last seen node */
   SCIP_Longint          ncleanups;          /**< number of storage cleanups */
   SCIP_Longint          nnzdualrays;        /**< number of non-zeros in all stored proofs based on dual rays */
   SCIP_Longint          nnzdualsols;        /**< number of non-zeros in all stored proofs based on dual solutions */
   int                   conflictsize;       /**< size of conflict array (bounded by conflict->maxpoolsize) */
   int                   origconflictsize;   /**< size of origconfs array */
   int                   nconflicts;         /**< number of stored conflicts */
   int                   ndualrayconfs;      /**< number of stored proofs based on dual rays */
   int                   ndualsolconfs;      /**< number of stored proofs based on dual solutions */
   int                   norigconfs;         /**< number of original conflicts */
   int                   ncbconflicts;       /**< number of conflicts depending on cutoff bound */
   int                   nconflictsfound;    /**< total number of conflicts found so far */
   int                   cleanupfreq;        /**< frequency to cleanup the storage if the storage is not full */
   int                   nswitches;          /**< number of path switches */
   int                   initstoresize;      /**< initial size of the storage (different to maxstoresize iff dynamic) */
   int                   storesize;          /**< current size of the storage (different to maxstoresize iff dynamic) */
   int                   maxstoresize;       /**< maximal size of the storage */
};

#ifdef __cplusplus
}
#endif

#endif
