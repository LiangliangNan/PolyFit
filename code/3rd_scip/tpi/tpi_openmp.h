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

/**@file   tpi_openmp.h
 * @ingroup TASKINTERFACE
 * @brief  the tpi_openmp redefines the lock functionality and some condition functionality as macros
 * @author Leona Gottwald
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifdef TPI_OMP

#ifndef _TPI_OPENMP_H_
#define _TPI_OPENMP_H_

/* define locks and some condition functionality as macros */

/* lock */
#define SCIPtpiInitLock(lock)     (omp_init_lock(lock), SCIP_OKAY)
#define SCIPtpiDestroyLock(lock)  (omp_destroy_lock(lock))
#define SCIPtpiAcquireLock(lock)   (omp_set_lock(lock), SCIP_OKAY)
#define SCIPtpiReleaseLock(lock)  (omp_unset_lock(lock), SCIP_OKAY)

/* condition */
#define SCIPtpiInitCondition(condition)    ( omp_init_lock(&(condition)->_lock), \
                                             (condition)->_waiters = 0, (condition)->_waitnum = 0, (condition)->_signals = 0, SCIP_OKAY )
#define SCIPtpiDestroyCondition(condition) do { assert((condition)->_waiters == 0); assert((condition)->_waitnum == 0); assert((condition)->_signals == 0); omp_destroy_lock(&(condition)->_lock); } while(0)
#endif

#endif
