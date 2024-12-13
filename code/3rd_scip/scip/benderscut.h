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

/**@file   benderscut.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for Benders' decomposition cuts
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_H__
#define __SCIP_BENDERSCUT_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_benderscut.h"
#include "scip/type_benders.h"
#include "scip/pub_benderscut.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given Benders' decomposition cut to a new scip */
SCIP_RETCODE SCIPbenderscutCopyInclude(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition that the cuts are copied to */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutCreate(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT**     benderscut,         /**< pointer to the Benders' decomposition cut data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of the Benders' decomposition cut */
   const char*           desc,               /**< description of the Benders' decomposition cut */
   int                   priority,           /**< priority of the Benders' decomposition cut */
   SCIP_Bool             islpcut,            /**< indicates whether the cut is generated from the LP solution */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy)),/**< copy method of the Benders' decomposition cut or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree)),/**< destructor of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit)),/**< initialize Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit)),/**< deinitialize Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol)),/**< solving process initialization method of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol)),/**< solving process deinitialization method of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec)),/**< execution method of the Benders' decomposition cut */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< Benders' decomposition cut data */
   );

/** calls destructor and frees memory of the Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutFree(
   SCIP_BENDERSCUT**     benderscut,         /**< pointer to the Benders' decomposition cut data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes the Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutInit(
   SCIP_BENDERSCUT*      benderscut,         /**< the Benders' decomposition cut */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of the Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutExit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs the Benders' decomposition cut that the branch and bound process is being started */
SCIP_RETCODE SCIPbenderscutInitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs the Benders' decomposition cut that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbenderscutExitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of the Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutExec(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the subproblem for which the cut is generated */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of the Benders' decomposition cut */
void SCIPbenderscutSetPriority(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   int                   priority            /**< new priority of the Benders' decomposition cut */
   );

/** sets copy callback of the Benders' decomposition cut */
void SCIPbenderscutSetCopy(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy))/**< copy callback of the Benders' decomposition cut or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor callback of the Benders' decomposition cut */
void SCIPbenderscutSetFree(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree))/**< destructor of the Benders' decomposition cut */
   );

/** sets initialization callback of the Benders' decomposition cut */
void SCIPbenderscutSetInit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit))/**< initialize the Benders' decomposition cut */
   );

/** sets deinitialization callback of the Benders' decomposition cut */
void SCIPbenderscutSetExit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit))/**< deinitialize the Benders' decomposition cut */
   );

/** sets solving process initialization callback of the Benders' decomposition cut */
void SCIPbenderscutSetInitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol))/**< solving process initialization callback of the Benders' decomposition cut */
   );

/** sets solving process deinitialization callback of the Benders' decomposition cut */
void SCIPbenderscutSetExitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol))/**< solving process deinitialization callback of the Benders' decomposition cut */
   );

#ifdef __cplusplus
}
#endif

#endif
