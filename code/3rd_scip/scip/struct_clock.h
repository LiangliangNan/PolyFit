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
