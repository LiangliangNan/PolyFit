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

/**@file   tpi_openmp.h
 * @ingroup TASKINTERFACE
 * @brief  the tpi_openmp redefines the lock functionality and some condition functionality as macros
 * @author Robert Lion Gottwald
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
