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

/**@file   sepa.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_H__
#define __SCIP_SEPA_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_sepastore.h"
#include "scip/type_sepa.h"
#include "scip/pub_sepa.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given separator to a new scip */
SCIP_RETCODE SCIPsepaCopyInclude(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a separator */
SCIP_RETCODE SCIPsepaCreate(
   SCIP_SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of separator */
   const char*           desc,               /**< description of separator */
   int                   priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling separator */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   SCIP_Bool             usessubscip,        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   SCIP_DECL_SEPACOPY    ((*sepacopy)),      /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   SCIP_DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol)),   /**< solving process initialization method of separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol)),   /**< solving process deinitialization method of separator */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp)),    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol)),   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   );

/** calls destructor and frees memory of separator */
SCIP_RETCODE SCIPsepaFree(
   SCIP_SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes separator */
SCIP_RETCODE SCIPsepaInit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of separator */
SCIP_RETCODE SCIPsepaExit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs separator that the branch and bound process is being started */
SCIP_RETCODE SCIPsepaInitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs separator that the branch and bound process data is being freed */
SCIP_RETCODE SCIPsepaExitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls LP separation method of separator */
SCIP_RETCODE SCIPsepaExecLP(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   depth,              /**< depth of current node */
   SCIP_Real             bounddist,          /**< current relative distance of local dual bound to global dual bound */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls primal solution separation method of separator */
SCIP_RETCODE SCIPsepaExecSol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of separator */
void SCIPsepaSetPriority(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the separator */
   );

/** sets copy method of separator */
void SCIPsepaSetCopy(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPACOPY    ((*sepacopy))       /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of separator */
void SCIPsepaSetFree(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAFREE    ((*sepafree))       /**< destructor of separator */
   );

/** sets initialization method of separator */
void SCIPsepaSetInit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAINIT    ((*sepainit))       /**< initialize separator */
   );

/** sets deinitialization method of separator */
void SCIPsepaSetExit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit))       /**< deinitialize separator */
   );

/** sets solving process initialization method of separator */
void SCIPsepaSetInitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol))    /**< solving process initialization method of separator */
   );

/** sets solving process deinitialization method of separator */
void SCIPsepaSetExitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol))    /**< solving process deinitialization method of separator */
   );

/** declares separator to be a parent separator */
void SCIPsepaSetIsParentsepa(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** sets the parent separator */
void SCIPsepaSetParentsepa(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPA*            parentsepa          /**< parent separator */
   );

/** enables or disables all clocks of \p sepa, depending on the value of the flag */
void SCIPsepaEnableOrDisableClocks(
   SCIP_SEPA*            sepa,               /**< the separator for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the separator be enabled? */
   );

/** increase count of applied cuts by one */
void SCIPsepaIncNCutsApplied(
   SCIP_SEPA*            sepa,                /**< separator */
   SCIP_Bool             fromcutpool          /**< whether the cuts were added from the cutpool to sepastore */
   );

/** increase count of found cuts by one */
void SCIPsepaIncNCutsAdded(
   SCIP_SEPA*            sepa,                /**< separator */
   SCIP_Bool             fromcutpool          /**< whether the cuts were added from the cutpool to sepastore */
   );

/** decrease the count of added cuts by one */
void SCIPsepaDecNCutsAdded(
   SCIP_SEPA*            sepa,                /**< separator */
   SCIP_Bool             fromcutpool          /**< whether the cuts were added from the cutpool to sepastore */
   );

/** increase count of found cuts by one */
void SCIPsepaIncNCutsFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** increase count of found cuts at current node by one */
void SCIPsepaIncNCutsFoundAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   );

#ifdef __cplusplus
}
#endif

#endif
