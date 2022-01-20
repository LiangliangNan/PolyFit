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

/**@file   struct_concurrent.h
 * @ingroup INTERNALAPI
 * @brief  concurrent data struct
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CONCURRENT_H__
#define __SCIP_STRUCT_CONCURRENT_H__

#include "scip/def.h"
#include "scip/type_concurrent.h"
#include "scip/type_clock.h"
#include "scip/type_concsolver.h"
#include "scip/type_prop.h"
#include "scip/type_heur.h"
#include "scip/type_event.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data for concurrent solve stored in struct scip */
struct SCIP_Concurrent
{
   SCIP*                    mainscip;           /**< main scip for concurrent solver */
   SCIP_CONCSOLVER*         concsolver;         /**< the concurrent solver of the main scip */
   int*                     varperm;            /**< permutation of variables to get the position of variable in the original SCIP's
                                                 *   variable array by the index of an original variable in this concurrent's main SCIP */
   SCIP_Real                dettime;            /**< deterministic time since last sync */
   SCIP_CLOCK*              wallclock;          /**< wallclock time since last sync */
   SCIP_PROP*               propsync;          /**< sync propagator */
   SCIP_HEUR*               heursync;          /**< sync heuristic */
   SCIP_EVENTHDLR*          eventglobalbnd;     /**< global bound eventhandler */
   int                      solidx;             /**< solution index after last synchronization */
};

#ifdef __cplusplus
}
#endif

#endif
