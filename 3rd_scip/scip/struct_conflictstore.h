/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
