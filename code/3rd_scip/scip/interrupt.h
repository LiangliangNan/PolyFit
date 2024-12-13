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

/**@file   interrupt.h
 * @ingroup INTERNALAPI
 * @brief  methods for catching the user CTRL-C interrupt
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_INTERRUPT_H__
#define __SCIP_INTERRUPT_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_interrupt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a CTRL-C interrupt data */
SCIP_RETCODE SCIPinterruptCreate(
   SCIP_INTERRUPT**      interrupt           /**< pointer to store the CTRL-C interrupt data */
   );

/** frees a CTRL-C interrupt data */
void SCIPinterruptFree(
   SCIP_INTERRUPT**      interrupt           /**< pointer to the CTRL-C interrupt data */
   );

/** captures the CTRL-C interrupt to call the SCIP's own interrupt handler */
void SCIPinterruptCapture(
   SCIP_INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   );

/** releases the CTRL-C interrupt and restores the old interrupt handler */
void SCIPinterruptRelease(
   SCIP_INTERRUPT*       interrupt           /**< CTRL-C interrupt data */
   );

/** returns whether the user interrupted by pressing CTRL-C */
SCIP_Bool SCIPinterrupted(
   void
   );

/** returns whether the process has received a SIGTERM */
SCIP_Bool SCIPterminated(
   void
   );

/** sends a termination signal to all SCIP processes so that they try to terminate as soon as possible
 *
 *  @note For terminating a specific SCIP process use SCIPinterruptSolve().
 */
void SCIPtryTerminate(
   void
   );

/** resets the number of interrupts to 0 */
void SCIPresetInterrupted(
   void
   );

#ifdef __cplusplus
}
#endif

#endif
