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

/**@file   concsolver.h
 * @brief  datastructures for concurrent solvers
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONCSOLVER_H__
#define __SCIP_CONCSOLVER_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_concsolver.h"
#include "scip/type_syncstore.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a concurrent solver type */
SCIP_RETCODE SCIPconcsolverTypeCreate(
   SCIP_CONCSOLVERTYPE** concsolvertype,     /**< pointer to concurrent solver data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of concurrent solver */
   SCIP_Real             prefpriodefault,    /**< the default preferred priority of this concurrent solver type */
   SCIP_DECL_CONCSOLVERCREATEINST ((*concsolvercreateinst)),/**< data copy method of concurrent solver */
   SCIP_DECL_CONCSOLVERDESTROYINST ((*concsolverdestroyinst)),/**< data copy method of concurrent solver */
   SCIP_DECL_CONCSOLVERINITSEEDS ((*concsolverinitseeds)),/**< initialize random seeds of concurrent solver */
   SCIP_DECL_CONCSOLVEREXEC ((*concsolverexec)),/**< execution method of concurrent solver */
   SCIP_DECL_CONCSOLVERCOPYSOLVINGDATA ((*concsolvercopysolvdata)),/**< method to copy solving data */
   SCIP_DECL_CONCSOLVERSTOP ((*concsolverstop)),/**< terminate solving in concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCWRITE ((*concsolversyncwrite)),/**< synchronization method of concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCREAD ((*concsolversyncread)),/**< synchronization method of concurrent solver */
   SCIP_DECL_CONCSOLVERTYPEFREEDATA ((*concsolvertypefreedata)),/**< method to free data of concurrent solver type */
   SCIP_CONCSOLVERTYPEDATA* data             /**< the concurent solver type's data */
   );

/** frees all memory of a concurrent solver type */
void SCIPconcsolverTypeFree(
   SCIP_CONCSOLVERTYPE** concsolvertype      /**< pointer to concurrent solver data structure */
   );

/** gets the data of a concurrent solver type */
SCIP_CONCSOLVERTYPEDATA* SCIPconcsolverTypeGetData(
   SCIP_CONCSOLVERTYPE*  concsolvertype      /**< concurrent solver type */
   );

/** sets the data of a concurrent solver type */
void SCIPconcsolverTypeSetData(
   SCIP_CONCSOLVERTYPE*  concsolvertype,     /**< concurrent solver type */
   SCIP_CONCSOLVERTYPEDATA* data             /**< the concurrent solver's data */
   );

/** gets the name of a concurrent solver type */
char* SCIPconcsolverTypeGetName(
   SCIP_CONCSOLVERTYPE*  concsolvertype      /**< concurrent solver type */
   );

/** gets the preferred priority from a concurrent solver type */
SCIP_Real SCIPconcsolverTypeGetPrefPrio(
   SCIP_CONCSOLVERTYPE*  concsolvertype      /**< concurrent solver type */
   );

/** creates an instance of the given concurrent solver type */
SCIP_RETCODE SCIPconcsolverCreateInstance(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONCSOLVERTYPE*  concsolvertype,     /**< concurrent solver type to create */
   SCIP_CONCSOLVER**     concsolver          /**< pointer to return concurrent solver instance */
   );

/** destroys an instance of the given concurrent solver */
SCIP_RETCODE SCIPconcsolverDestroyInstance(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONCSOLVER**     concsolver          /**< concurrent solver */
   );

/** gets the data of a concurrent solver */
SCIP_CONCSOLVERDATA* SCIPconcsolverGetData(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** sets the data of a concurrent solver */
void SCIPconcsolverSetData(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   SCIP_CONCSOLVERDATA*  data                /**< the concurrent solver's data */
   );

/** gets the name of a concurrent solver */
char* SCIPconcsolverGetName(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** initializes the random seeds of a concurrent solver */
SCIP_RETCODE SCIPconcsolverInitSeeds(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   unsigned int          seed                /**< seed for initializing the solver's internal random seeds */
   );

/** start the solving process of a concurrent solver */
SCIP_RETCODE SCIPconcsolverExec(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets solving data of concurrent solver and stores it in the given SCIP instance */
SCIP_RETCODE SCIPconcsolverGetSolvingData(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   SCIP*                 scip                /**< SCIP datastructure */
   );

/** interrupt solving in a concurrent solver */
SCIP_RETCODE SCIPconcsolverStop(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** let the given concurrent solver synchronize, i.e. pass its own solutions and bounds to
 *  the SPI.
 */
SCIP_RETCODE SCIPconcsolverSync(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets the current synchronization frequency of the concurent solver */
SCIP_Real SCIPconcsolverGetSyncFreq(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets the total memory used by the concurent solver */
SCIP_Longint SCIPconcsolverGetMemTotal(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** sets the time elapsed since the last synchronization. Must be set before the synchronization is
 *  started.
 */
void SCIPconcsolverSetTimeSinceLastSync(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   SCIP_Real             time                /**< the time passed since the last synchronization */
   );

/** gets the solving time of the concurrent solver */
SCIP_Real SCIPconcsolverGetSolvingTime(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets the time spent for synchronization for the concurrent solver */
SCIP_Real SCIPconcsolverGetSyncTime(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets the number of lp iterations the concurrent solver used */
SCIP_Longint SCIPconcsolverGetNLPIterations(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets the number of branch and bound nodes the concurrent solver used */
SCIP_Longint SCIPconcsolverGetNNodes(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets the number of solutions the concurrent solver received during synchronization */
SCIP_Longint SCIPconcsolverGetNSolsRecvd(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets the number of solutions the concurrent solver shared during synchronization */
SCIP_Longint SCIPconcsolverGetNSolsShared(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets the number of tighter global variable bounds the solver received */
SCIP_Longint SCIPconcsolverGetNTighterBnds(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets the number of tighter global variable bounds of integer variables the solver received */
SCIP_Longint SCIPconcsolverGetNTighterIntBnds(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

/** gets index of concurrent solver */
int SCIPconcsolverGetIdx(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   );

#ifdef __cplusplus
}
#endif

#endif
