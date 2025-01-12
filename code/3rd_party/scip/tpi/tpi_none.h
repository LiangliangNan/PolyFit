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

/**@file   tpi_none.h
 * @ingroup TASKINTERFACE
 * @brief  the dummy implementation defines all functions as macros
 * @author Leona Gottwald
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifdef TPI_NONE

#include "scip/def.h"

#ifndef _TPI_NONE_H_
#define _TPI_NONE_H_

/* lock */
#define SCIPtpiInitLock(lock)                 (SCIP_UNUSED(lock), SCIP_OKAY)
#define SCIPtpiDestroyLock(lock)              SCIP_UNUSED(lock)
#define SCIPtpiAcquireLock(lock)              (SCIP_UNUSED(lock), SCIP_OKAY)
#define SCIPtpiReleaseLock(lock)              (SCIP_UNUSED(lock), SCIP_OKAY)

/* condition */
#define SCIPtpiInitCondition(condition)       (SCIP_UNUSED(condition), SCIP_OKAY)
#define SCIPtpiDestroyCondition(condition)    SCIP_UNUSED(condition)
#define SCIPtpiSignalCondition(condition)     (SCIP_UNUSED(condition), SCIP_OKAY)
#define SCIPtpiBroadcastCondition(condition)  (SCIP_UNUSED(condition), SCIP_OKAY)
#define SCIPtpiWaitCondition(condition, lock) /*lint -e505*/ (SCIP_UNUSED(condition), SCIP_UNUSED(lock), SCIP_OKAY)


#define SCIPtpiGetNumThreads()                1
#define SCIPtpiGetThreadNum()                 0

#define SCIPtpiGetNewJobID()                  0

#endif

#endif
