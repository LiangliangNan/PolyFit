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

/**@file   presol.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_H__
#define __SCIP_PRESOL_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_presol.h"
#include "scip/pub_presol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given presolver to a new scip */
SCIP_RETCODE SCIPpresolCopyInclude(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a presolver */
SCIP_RETCODE SCIPpresolCreate(
   SCIP_PRESOL**         presol,             /**< pointer to store presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of presolver */
   const char*           desc,               /**< description of presolver */
   int                   priority,           /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
   int                   maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   SCIP_PRESOLTIMING     timing,             /**< timing mask of the presolver */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy)),    /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   SCIP_DECL_PRESOLINIT  ((*presolinit)),    /**< initialization method of presolver (called after problem was transformed) */
   SCIP_DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialization method of presolver (called before transformed problem is freed) */
   SCIP_DECL_PRESOLINITPRE((*presolinitpre)),/**< presolving initialization method of presolver (called when presolving is about to begin) */
   SCIP_DECL_PRESOLEXITPRE((*presolexitpre)),/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   SCIP_DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   );

/** frees memory of presolver */
SCIP_RETCODE SCIPpresolFree(
   SCIP_PRESOL**         presol,             /**< pointer to presolver data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes presolver */
SCIP_RETCODE SCIPpresolInit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** deinitializes presolver */
SCIP_RETCODE SCIPpresolExit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs presolver that the presolving process is being started */
SCIP_RETCODE SCIPpresolInitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs presolver that the presolving process is finished */
SCIP_RETCODE SCIPpresolExitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** executes presolver */
SCIP_RETCODE SCIPpresolExec(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRESOLTIMING     timing,             /**< current presolving timing */
   int                   nrounds,            /**< number of presolving rounds already done */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of presolver */
void SCIPpresolSetPriority(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the presolver */
   );

/** sets copy method of presolver */
void SCIPpresolSetCopy(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy))     /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of presolver */
void SCIPpresolSetFree(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLFREE  ((*presolfree))     /**< destructor of presolver */
   );

/** sets initialization method of presolver */
void SCIPpresolSetInit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLINIT  ((*presolinit))     /**< initialize presolver */
   );

/** sets deinitialization method of presolver */
void SCIPpresolSetExit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLEXIT  ((*presolexit))     /**< deinitialize presolver */
   );

/** sets solving process initialization method of presolver */
void SCIPpresolSetInitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLINITPRE ((*presolinitpre))/**< solving process initialization method of presolver */
   );

/** sets solving process deinitialization method of presolver */
void SCIPpresolSetExitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLEXITPRE ((*presolexitpre))/**< solving process deinitialization method of presolver */
   );

/** enables or disables all clocks of \p presol, depending on the value of the flag */
void SCIPpresolEnableOrDisableClocks(
   SCIP_PRESOL*          presol,             /**< the presolver for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the presolver be enabled? */
   );

#ifdef __cplusplus
}
#endif

#endif
