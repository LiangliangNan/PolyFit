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

/**@file   type_syncstore.h
 * @ingroup PARALLEL
 * @brief  the type definitions for the synchronization store
 * @author Stephen J. Maher
 * @author Leona Gottwald
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
