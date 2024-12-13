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

/**@file   type_clock.h
 * @brief  type definitions for clocks and timing issues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_CLOCK_H__
#define __SCIP_TYPE_CLOCK_H__

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

enum SCIP_ClockType
{
   SCIP_CLOCKTYPE_DEFAULT = 0,          /**< use default clock type */
   SCIP_CLOCKTYPE_CPU     = 1,          /**< use CPU clock */
   SCIP_CLOCKTYPE_WALL    = 2           /**< use wall clock */
};
typedef enum SCIP_ClockType SCIP_CLOCKTYPE;       /**< clock type to use */

typedef struct SCIP_Clock SCIP_CLOCK;             /**< clock timer */
typedef struct SCIP_CPUClock SCIP_CPUCLOCK;       /**< CPU clock counter */
typedef struct SCIP_WallClock SCIP_WALLCLOCK;     /**< wall clock counter */

#ifdef __cplusplus
}
#endif

#endif
