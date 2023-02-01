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

/**@file   tpi_tnycthrd.h
 * @ingroup TASKINTERFACE
 * @brief  the tinycthreads implementation defines the lock and condition functionality as macros
 * @author Robert Lion Gottwald
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
