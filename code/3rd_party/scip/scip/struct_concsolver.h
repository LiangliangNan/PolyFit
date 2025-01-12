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

/**@file   struct_concsolver.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for concurrent solvers
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CONCSOLVER_H__
#define __SCIP_STRUCT_CONCSOLVER_H__


#include "scip/def.h"
#include "scip/type_concsolver.h"
#include "scip/type_clock.h"

#ifdef __cplusplus
extern "C" {
#endif

/** concurrent solver data structure */
struct SCIP_ConcSolverType
{
   int                                 ninstances;                 /**< number of instances created from this concurrent solver type */
   SCIP_Real                           prefprio;                   /**< the weight of the concurrent */
   char*                               name;                       /**< name of concurrent solver */
   SCIP_CONCSOLVERTYPEDATA*            data;                       /**< user data of concurrent solver type */
   SCIP_DECL_CONCSOLVERCREATEINST      ((*concsolvercreateinst));  /**< creates an instance of the concurrent solver */
   SCIP_DECL_CONCSOLVERDESTROYINST     ((*concsolverdestroyinst)); /**< destroys an instance of the concurrent solver */
   SCIP_DECL_CONCSOLVERINITSEEDS       ((*concsolverinitseeds));   /**< initialize random seeds of concurrent solver */
   SCIP_DECL_CONCSOLVEREXEC            ((*concsolverexec));        /**< execution method of concurrent solver */
   SCIP_DECL_CONCSOLVERCOPYSOLVINGDATA ((*concsolvercopysolvdata));/**< copies the solving data */
   SCIP_DECL_CONCSOLVERSTOP            ((*concsolverstop));        /**< terminate solving in concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCWRITE       ((*concsolversyncwrite));   /**< synchronization method of concurrent solver for sharing it's data */
   SCIP_DECL_CONCSOLVERSYNCREAD        ((*concsolversyncread));    /**< synchronization method of concurrent solver for reading shared data */
   SCIP_DECL_CONCSOLVERTYPEFREEDATA    ((*concsolvertypefreedata));/**< frees user data of concurrent solver type */
};

/** concurrent solver data structure */
struct SCIP_ConcSolver
{
   SCIP_CONCSOLVERTYPE*                type;                      /**< type of this concurrent solver */
   int                                 idx;                       /**< index of initialized exernal solver */
   char*                               name;                      /**< name of concurrent solver */
   SCIP_CONCSOLVERDATA*                data;                      /**< user data of concurrent solver */
   SCIP_SYNCDATA*                      syncdata;                  /**< most recent synchronization data that has been read */
   SCIP_Longint                        nsyncs;                    /**< total number of synchronizations */
   SCIP_Real                           timesincelastsync;         /**< time since the last synchronization */
   SCIP_Real                           syncdelay;                 /**< current delay of synchronization data */
   SCIP_Real                           syncfreq;                  /**< current synchronization frequency of the concurrent solver */
   SCIP_Real                           solvingtime;               /**< solving time with wall clock */
   SCIP_Bool                           stopped;                   /**< flag to store if the concurrent solver has been stopped
                                                                   *   through the SCIPconcsolverStop function */
   SCIP_Longint                        nlpiterations;             /**< number of lp iterations the concurrent solver used */
   SCIP_Longint                        nnodes;                    /**< number of nodes the concurrent solver used */
   SCIP_Longint                        nsolsrecvd;                /**< number of solutions the concurrent solver received */
   SCIP_Longint                        nsolsshared;               /**< number of solutions the concurrent solver found */
   SCIP_Longint                        ntighterbnds;              /**< number of tighter global variable bounds the concurrent solver received */
   SCIP_Longint                        ntighterintbnds;           /**< number of tighter global variable bounds the concurrent solver received
                                                                   *   on integer variables */
   SCIP_CLOCK*                         totalsynctime;             /**< total time used for synchronization, including idle time */
};

#ifdef __cplusplus
}
#endif

#endif
