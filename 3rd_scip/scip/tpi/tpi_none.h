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

/**@file   tpi_none.h
 * @ingroup TASKINTERFACE
 * @brief  the dummy implementation defines all functions as macros
 * @author Robert Lion Gottwald
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
