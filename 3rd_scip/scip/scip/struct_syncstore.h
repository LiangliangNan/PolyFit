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

/**@file   struct_syncstore.h
 * @ingroup PARALLEL
 * @brief  the struct definitions for the synchronization store
 * @author Stephen J. Maher
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_SYNCSTORE_H__
#define __STRUCT_SYNCSTORE_H__

#include "scip/type_syncstore.h"
#include "tpi/type_tpi.h"
#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct SCIP_SyncStore
{
   int                   nuses;              /**< number of uses of the synchronization store */
   SCIP_PARALLELMODE     mode;               /**< the mode for the parallel solving */
   SCIP_Bool             initialized;        /**< flag to indicate whether the syncstore has been initialized */
   int                   ninitvars;          /**< number of variables it has been initialized for */
   SCIP_SYNCDATA*        syncdata;           /**< array of size nsyncdata, containing the synchronization data
                                              *   for each active synchroization */
   SCIP_SYNCDATA*        lastsync;           /**< pointer to the last synchronization data that has been synchronized
                                              *   by all threads */

   SCIP*                 mainscip;           /**< the SCIP instance that was used for initializing the syncstore */
   SCIP_Bool             stopped;            /**< flag to indicate if the solving is stopped */
   SCIP_LOCK             lock;               /**< lock to protect the syncstore data structure from data races */

   /* SPI settings */
   int                   nsyncdata;          /**< the size of the synchronization data array */
   SCIP_Real             minsyncdelay;       /**< the minimum delay before a synchronization data may be read */
   int                   maxnsyncdelay;      /**< maximum number of synchronizations before the reading of the next
                                              *   synchronization data is enforced regardless of the minimal synchroization
                                              *   delay */
   SCIP_Real             syncfreqinit;       /**< the initial synchronization frequency which is read from the settings
                                              *   of the main SCIP when the syncstore is initialized */
   SCIP_Real             syncfreqmax;        /**< the maximum synchronization frequency */
   int                   maxnsols;           /**< maximum number of solutions that can be shared in one synchronization */
   int                   nsolvers;           /**< number of solvers synchronizing with this syncstore */
};


struct SCIP_SyncData
{
   SCIP_Real*            solobj;             /**< array with the objective value of all stored solutions */
   SCIP_Real**           sols;               /**< array with the solution values of each variable for all stored solutions */
   int*                  solsource;          /**< the solverid of the solution came from */
   int                   nsols;              /**< number of solutions currently stored in the synchronization data */
   SCIP_Real             bestlowerbound;     /**< largest lower bound on the objective value that was stored in this
                                              *   synchroization data */
   SCIP_Real             bestupperbound;     /**< smalles upper bound on the objective value that was stored in this
                                              *   synchroization data */
   SCIP_Longint          syncnum;            /**< the synchronization number of this synchronization data */
   int                   winner;             /**< the solverid of the solver with the best status */
   SCIP_STATUS           status;             /**< the best status that was stored in this synchronization data */
   SCIP_LOCK             lock;               /**< a lock to protect this synchronization data */
   int                   syncedcount;        /**< a counter of how many solvers have finished writing to this synchronization data */
   SCIP_CONDITION        allsynced;          /**< a condition variable to signal when the last solver has finished writing to this
                                              *   synchronization data */
   SCIP_BOUNDSTORE*      boundstore;         /**< a boundstore for storing all the bound changes that were added to this
                                              *   synchronization data */
   SCIP_Real             syncfreq;           /**< the synchroization frequency that was set in this synchronization data */
   SCIP_Longint          memtotal;           /**< the total amount of memory used by all solvers including the main SCIP */
};

/** struct for storing the position of avariables lower and upper bound in the boundstore */
typedef struct
{
   int                   pos[2];             /**< stores at pos[SCIP_BOUNDTYPE_LOWER] the position of the lowerbound and
                                              *   at pos[SCIP_BOUNDTYPE_UPPER] the position of the upperbound */
} BoundPos;

/** struct for storing a single boundchange in the boundstore */
typedef struct
{
   int                   varidx;             /**< the variables position in the variable array of the main scip */
   SCIP_Real             newbound;           /**< the variables new bound */
   SCIP_BOUNDTYPE        boundtype;          /**< the type of the variables new bound */
} BoundChg;

struct SCIP_BoundStore
{
   int                   nvars;              /**< the number of variables to store bounds for */
   BoundPos*             bndpos;             /**< array of size nvars to store the positions for all the bound changes
                                              *   stored in this boundstore */
   BoundChg*             bndchg;             /**< array of boundchanges */
   int                   nbndchg;            /**< the number of boundchanges stored in this bound store */
   int                   bndchgsize;         /**< the size of the bound change array */
};

#ifdef __cplusplus
}
#endif

#endif
