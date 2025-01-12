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

/**@file   scip_relax.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for relaxator plugins
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

#ifndef __SCIP_SCIP_RELAX_H__
#define __SCIP_SCIP_RELAX_H__


#include "scip/def.h"
#include "scip/type_relax.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicRelaxatorMethods
 *
 * @{
 */

/** creates a relaxation handler and includes it in SCIP
 *
 *  @note method has all relaxation handler callbacks as arguments and is thus changed every time a new
 *        callback is added
 *        in future releases; consider using SCIPincludeRelaxBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of relaxation handler */
   const char*           desc,               /**< description of relaxation handler */
   int                   priority,           /**< priority of the relaxation handler (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxation handler */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy)),     /**< copy method of relaxation handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxation handler */
   SCIP_DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxation handler */
   SCIP_DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxation handler */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxation handler */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxation handler */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< relaxation handler data */
   );

/** creates a relaxation handler and includes it in SCIP. All non fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetRelaxInit(), SCIPsetRelaxExit(),
 *  SCIPsetRelaxCopy(), SCIPsetRelaxFree(), SCIPsetRelaxInitsol(), and SCIPsetRelaxExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeRelax() instead
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeRelaxBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX**          relaxptr,           /**< reference to relaxation pointer, or NULL */
   const char*           name,               /**< name of relaxation handler */
   const char*           desc,               /**< description of relaxation handler */
   int                   priority,           /**< priority of the relaxation handler (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxation handler */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< relaxation handler data */
   );

/** sets copy method of relaxation handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetRelaxCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy))      /**< copy method of relaxation handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of relaxation handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetRelaxFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXFREE   ((*relaxfree))      /**< destructor of relaxation handler */
   );

/** sets initialization method of relaxation handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetRelaxInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXINIT   ((*relaxinit))      /**< initialize relaxation handler */
   );

/** sets deinitialization method of relaxation handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetRelaxExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXEXIT   ((*relaxexit))      /**< deinitialize relaxation handler */
   );

/** sets solving process initialization method of relaxation handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetRelaxInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol))   /**< solving process initialization method of relaxation handler */
   );

/** sets solving process deinitialization method of relaxation handler */
SCIP_EXPORT
SCIP_RETCODE SCIPsetRelaxExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol))   /**< solving process deinitialization method of relaxation handler */
   );

/** returns the relaxation handler of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_RELAX* SCIPfindRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of relaxation handler */
   );

/** returns the array of currently available relaxation handlers  */
SCIP_EXPORT
SCIP_RELAX** SCIPgetRelaxs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available relaxation handlers  */
SCIP_EXPORT
int SCIPgetNRelaxs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a relaxation handler*/
SCIP_EXPORT
SCIP_RETCODE SCIPsetRelaxPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   int                   priority            /**< new priority of the relaxation handler */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
