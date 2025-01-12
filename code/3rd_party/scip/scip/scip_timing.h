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

/**@file   scip_timing.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for timing
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_TIMING_H__
#define __SCIP_SCIP_TIMING_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicTimingMethods
 *
 * @{
 */

/** gets current time of day in seconds (standard time zone)
 *
 *  @return the current time of day in seconds (standard time zone).
 */
SCIP_EXPORT
SCIP_Real SCIPgetTimeOfDay(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a clock using the default clock type
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** creates a clock counting the CPU user seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateCPUClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** creates a clock counting the wall clock seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateWallClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** frees a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   );

/** resets the time measurement of a clock to zero and completely stops the clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPresetClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** starts the time measurement of a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPstartClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** stops the time measurement of a clock
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPstopClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** enables or disables \p clck */
SCIP_EXPORT
void SCIPsetClockEnabled(
   SCIP_CLOCK*           clck,               /**< the clock to be disabled/enabled */
   SCIP_Bool             enable              /**< should the clock be enabled or disabled? */
   );

/** enables or disables all statistic clocks of SCIP concerning plugin statistics,
 *  LP execution time, strong branching time, etc.
 *
 *  Method reads the value of the parameter timing/statistictiming. In order to disable statistic timing,
 *  set the parameter to FALSE.
 *
 *  @note: The (pre-)solving time clocks which are relevant for the output during (pre-)solving
 *         are not affected by this method
 *
 *  @see: For completely disabling all timing of SCIP, consider setting the parameter timing/enabled to FALSE
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPenableOrDisableStatisticTiming(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** starts the current solving time
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPstartSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** stops the current solving time in seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPstopSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the measured time of a clock in seconds
 *
 *  @return the measured time of a clock in seconds.
 */
SCIP_EXPORT
SCIP_Real SCIPgetClockTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   );

/** sets the measured time of a clock to the given value in seconds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsetClockTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_Real             sec                 /**< time in seconds to set the clock's timer to */
   );

/** gets the current total SCIP time in seconds, possibly accumulated over several problems.
 *
 *  @return the current total SCIP time in seconds, ie. the total time since the SCIP instance has been created
 */
SCIP_EXPORT
SCIP_Real SCIPgetTotalTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current solving time in seconds
 *
 *  @return the current solving time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_Real SCIPgetSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current reading time in seconds
 *
 *  @return the current reading time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_Real SCIPgetReadingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the current presolving time in seconds
 *
 *  @return the current presolving time in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_Real SCIPgetPresolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the time need to solve the first LP in the root node
 *
 *  @return the solving time for the first LP in the root node in seconds.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_Real SCIPgetFirstLPTime(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
