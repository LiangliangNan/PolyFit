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

/**@file   struct_cutpool.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CUTPOOL_H__
#define __SCIP_STRUCT_CUTPOOL_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_misc.h"
#include "scip/type_lp.h"
#include "scip/type_cutpool.h"

#ifdef __cplusplus
extern "C" {
#endif

/** datastructure for cuts in a cut pool */
struct SCIP_Cut
{
   SCIP_ROW*             row;                /**< LP row of this cut */
   SCIP_Longint          processedlp;        /**< last LP, where this cut was processed in separation of the LP solution */
   SCIP_Longint          processedlpsol;     /**< last LP, where this cut was processed in separation of other solutions */
   int                   age;                /**< age of the cut: number of successive times, the cut was not violated */
   int                   pos;                /**< position of cut in the cuts array of the cut pool */
};

/** storage for pooled cuts */
struct SCIP_Cutpool
{
   SCIP_Longint          ncalls;             /**< number of times, the cutpool was separated */
   SCIP_Longint          ncutsfound;         /**< total number of cuts that were separated from the pool */
   SCIP_CLOCK*           poolclock;          /**< separation time */
   SCIP_HASHTABLE*       hashtable;          /**< hash table to identify already stored cuts */
   SCIP_CUT**            cuts;               /**< stored cuts of the pool */
   SCIP_Longint          processedlp;        /**< last LP that has been processed for separating the LP */
   SCIP_Longint          processedlpsol;     /**< last LP that has been processed for separating other solutions */
   SCIP_Real             processedlpefficacy;/**< minimal efficacy used in last processed LP */
   SCIP_Real             processedlpsolefficacy;/**< minimal efficacy used in last processed LP for separating other solutions */
   int                   cutssize;           /**< size of cuts array */
   int                   ncuts;              /**< number of cuts stored in the pool */
   int                   nremovablecuts;     /**< number of cuts stored in the pool that are marked to be removable */
   int                   agelimit;           /**< maximum age a cut can reach before it is deleted from the pool */
   int                   firstunprocessed;   /**< first cut that has not been processed in the last LP */
   int                   firstunprocessedsol;/**< first cut that has not been processed in the last LP when separating other solutions */
   int                   maxncuts;           /**< maximal number of cuts stored in the pool at the same time */
   SCIP_Bool             globalcutpool;      /**< is this the global cut pool of SCIP? */
};

#ifdef __cplusplus
}
#endif

#endif
