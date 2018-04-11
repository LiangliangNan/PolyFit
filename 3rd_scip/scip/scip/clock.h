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

/**@file   clock.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for clocks and timing issues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CLOCK_H__
#define __SCIP_CLOCK_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_clock.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a clock and initializes it */
extern
SCIP_RETCODE SCIPclockCreate(
   SCIP_CLOCK**          clck,               /**< pointer to clock timer */
   SCIP_CLOCKTYPE        clocktype           /**< type of clock */
   );

/** frees a clock */
extern
void SCIPclockFree(
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** initializes and resets a clock */
extern
void SCIPclockInit(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_CLOCKTYPE        clocktype           /**< type of clock */
   );

/** completely stop the clock and reset the clock's counter to zero */
extern
void SCIPclockReset(
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** enables the clock */
extern
void SCIPclockEnable(
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** disables and resets the clock */
extern
void SCIPclockDisable(
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** enables or disables \p clck, depending on the value of the flag */
extern
void SCIPclockEnableOrDisable(
   SCIP_CLOCK*           clck,               /**< the clock to be disabled/enabled */
   SCIP_Bool             enable              /**< should the clock be enabled? */
   );

/** sets the type of the clock, overriding the default clock type, and resets the clock */
extern
void SCIPclockSetType(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_CLOCKTYPE        clocktype           /**< type of clock */
   );

/** starts measurement of time in the given clock, update the clock's type if it is bound to the default type */
extern
void SCIPclockStart(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** stops measurement of time in the given clock */
extern
void SCIPclockStop(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns whether the clock is currently running */
extern
SCIP_Bool SCIPclockIsRunning(
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** gets the used time of this clock in seconds */
extern
SCIP_Real SCIPclockGetTime(
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** gets the last validated time of this clock in seconds */
extern
SCIP_Real SCIPclockGetLastTime(
   SCIP_CLOCK*           clck                /**< clock timer */
   );


/** sets the used time of this clock in seconds */
extern
void SCIPclockSetTime(
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_Real             sec                 /**< time in seconds to set the clock's timer to */
   );

/** gets current time of day in seconds (standard time zone) */
extern
SCIP_Real SCIPclockGetTimeOfDay(
   void
   );

#ifdef __cplusplus
}
#endif

#endif
