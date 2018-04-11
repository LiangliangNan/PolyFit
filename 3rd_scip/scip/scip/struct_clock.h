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

/**@file   struct_clock.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for clocks and timing issues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CLOCK_H__
#define __SCIP_STRUCT_CLOCK_H__


#if defined(_WIN32) || defined(_WIN64)
#include <time.h>
#else
#include <sys/times.h>
#endif

#include "scip/def.h"
#include "scip/type_clock.h"

#ifdef __cplusplus
extern "C" {
#endif

/** CPU clock counter */
struct SCIP_CPUClock
{
   clock_t               user;               /**< clock ticks for user CPU time */
};

/** wall clock counter */
struct SCIP_WallClock
{
   long                  sec;                /**< seconds counter */
   long                  usec;               /**< microseconds counter */
};

/** clock timer */
struct SCIP_Clock
{
   union
   {
      SCIP_CPUCLOCK      cpuclock;           /**< CPU clock counter */
      SCIP_WALLCLOCK     wallclock;          /**< wall clock counter */
   } data;
   SCIP_Real             lasttime;           /**< last validated time of clock */
   int                   nruns;              /**< number of SCIPclockStart() calls without SCIPclockStop() calls */
   SCIP_CLOCKTYPE        clocktype;          /**< current type of clock used */
   SCIP_Bool             usedefault;         /**< should the clock's type be overruled by the default clock type? */
   SCIP_Bool             enabled;            /**< should the clock be used? */
};

#ifdef __cplusplus
}
#endif

#endif
