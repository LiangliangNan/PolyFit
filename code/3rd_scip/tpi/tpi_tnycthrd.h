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

/**@file   tpi_tnycthrd.h
 * @ingroup TASKINTERFACE
 * @brief  the tinycthreads implementation defines the lock and condition functionality as macros
 * @author Leona Gottwald
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifdef TPI_TNYC

#ifndef _TPI_TINYCTHRD_H_
#define _TPI_TINYCTHRD_H_

/* lock */
#define SCIPtpiInitLock(lock)                 ( mtx_init((lock), mtx_plain) == thrd_success ? SCIP_OKAY : SCIP_ERROR )
#define SCIPtpiDestroyLock(lock)              ( mtx_destroy(lock) )
#define SCIPtpiAcquireLock(lock)               ( mtx_lock(lock) == thrd_success ? SCIP_OKAY : SCIP_ERROR )
#define SCIPtpiReleaseLock(lock)              ( mtx_unlock(lock) == thrd_success ? SCIP_OKAY : SCIP_ERROR )

/* condition */
#define SCIPtpiInitCondition(condition)       ( cnd_init(condition) == thrd_success ? SCIP_OKAY : SCIP_ERROR )
#define SCIPtpiDestroyCondition(condition)    ( cnd_destroy(condition) )
#define SCIPtpiSignalCondition(condition)     ( cnd_signal(condition) == thrd_success ? SCIP_OKAY : SCIP_ERROR )
#define SCIPtpiBroadcastCondition(condition)  ( cnd_broadcast(condition) == thrd_success ? SCIP_OKAY : SCIP_ERROR )
#define SCIPtpiWaitCondition(condition, lock) ( cnd_wait((condition), (lock)) == thrd_success ? SCIP_OKAY: SCIP_ERROR )


extern _Thread_local int _threadnumber;

#define SCIPtpiGetThreadNum()                 /*lint -e40*/ _threadnumber

#endif

#endif
