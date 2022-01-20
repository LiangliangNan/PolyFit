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

/**@file   type_tpi_openmp.h
 * @ingroup TASKINTERFACE
 * @brief  the type definitions for the the locks and condition variables for the OpenMP implementation of the TPI
 * @author Robert Lion Gottwald
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifdef TPI_OMP

#ifndef __TYPE_TPI_OPENMP_H__
#define __TYPE_TPI_OPENMP_H__

#include <omp.h>

typedef omp_lock_t SCIP_LOCK;
typedef struct {
   SCIP_LOCK  _lock;
   int        _waiters;
   int        _waitnum;
   int        _signals;
} SCIP_CONDITION;

#endif

#endif
