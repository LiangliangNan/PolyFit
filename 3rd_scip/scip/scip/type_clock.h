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
