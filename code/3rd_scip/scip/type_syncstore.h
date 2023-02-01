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

/**@file   type_syncstore.h
 * @ingroup PARALLEL
 * @brief  the type definitions for the synchronization store
 * @author Stephen J. Maher
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_SPI_H__
#define __TYPE_SPI_H__

#ifdef __cplusplus
extern "C" {
#endif

/** The parallel mode */
enum SCIP_Parallelmode
{
   SCIP_PARA_OPPORTUNISTIC   = 0,
   SCIP_PARA_DETERMINISTIC   = 1
};
typedef enum SCIP_Parallelmode SCIP_PARALLELMODE;

typedef struct SCIP_SyncStore SCIP_SYNCSTORE;   /**< structure to store information for synchronization */
typedef struct SCIP_SyncData SCIP_SYNCDATA;     /**< data for a single synchronization */
typedef struct SCIP_BoundStore SCIP_BOUNDSTORE; /**< structure to store boundchanges for synchronization */

#ifdef __cplusplus
}
#endif

#endif
