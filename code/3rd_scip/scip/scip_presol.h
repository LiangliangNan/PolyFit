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

/**@file   scip_presol.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for presolving plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_PRESOL_H__
#define __SCIP_SCIP_PRESOL_H__


#include "scip/def.h"
#include "scip/type_presol.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_timing.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicPresolverMethods
 *
 * @{
 */

/** creates a presolver and includes it in SCIP
 *
 *  @note method has all presolver callbacks as arguments and is thus changed every time a new
 *        callback is added
 *        in future releases; consider using SCIPincludePresolBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresol(
   SCIP*                 scip,               /**< SCIP data structure */
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

/** Creates a presolver and includes it in SCIP with its fundamental callback. All non-fundamental (or optional)
 *  callbacks as, e.g., init and exit callbacks, will be set to NULL. Optional callbacks can be set via specific setter
 *  functions. These are SCIPsetPresolCopy(), SCIPsetPresolFree(), SCIPsetPresolInit(), SCIPsetPresolExit(),
 *  SCIPsetPresolInitpre(), and SCIPsetPresolExitPre().
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludePresol() instead
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresolBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL**         presolptr,          /**< reference to presolver, or NULL */
   const char*           name,               /**< name of presolver */
   const char*           desc,               /**< description of presolver */
   int                   priority,           /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
   int                   maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   SCIP_PRESOLTIMING     timing,             /**< timing mask of the presolver */
   SCIP_DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   );

/** sets copy method of presolver */
SCIP_EXPORT
SCIP_RETCODE SCIPsetPresolCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy))     /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of presolver */
SCIP_EXPORT
SCIP_RETCODE SCIPsetPresolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLFREE ((*presolfree))      /**< destructor of presolver */
   );

/** sets initialization method of presolver */
SCIP_EXPORT
SCIP_RETCODE SCIPsetPresolInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLINIT  ((*presolinit))     /**< initialize presolver */
   );

/** sets deinitialization method of presolver */
SCIP_EXPORT
SCIP_RETCODE SCIPsetPresolExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLEXIT  ((*presolexit))     /**< deinitialize presolver */
   );

/** sets solving process initialization method of presolver */
SCIP_EXPORT
SCIP_RETCODE SCIPsetPresolInitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLINITPRE ((*presolinitpre))/**< solving process initialization method of presolver */
   );

/** sets solving process deinitialization method of presolver */
SCIP_RETCODE SCIPsetPresolExitpre(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLEXITPRE ((*presolexitpre))/**< solving process deinitialization method of presolver */
   );

/** returns the presolver of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_PRESOL* SCIPfindPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of presolver */
   );

/** returns the array of currently available presolvers */
SCIP_EXPORT
SCIP_PRESOL** SCIPgetPresols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available presolvers */
SCIP_EXPORT
int SCIPgetNPresols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a presolver */
SCIP_EXPORT
SCIP_RETCODE SCIPsetPresolPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   int                   priority            /**< new priority of the presolver */
   );

/** returns the number of presolve rounds (current or last presolve) */
SCIP_EXPORT
int SCIPgetNPresolRounds(
   SCIP*                 scip                /**< SCIP data structure */
);

/** @} */

#ifdef __cplusplus
}
#endif

#endif
