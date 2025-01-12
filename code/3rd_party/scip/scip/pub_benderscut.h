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

/**@file   pub_benderscut.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for Benders' decomposition cuts
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_BENDERSCUT_H__
#define __SCIP_PUB_BENDERSCUT_H__

#include "scip/def.h"
#include "scip/type_benderscut.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicBenderscutsMethods
 *
 * @{
 */

/** compares two Benders' decomposition cuts w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPbenderscutComp);

/** comparison method for sorting Benders' decomposition cuts w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPbenderscutCompName);

/** gets user data of the Benders' decomposition cut */
SCIP_EXPORT
SCIP_BENDERSCUTDATA* SCIPbenderscutGetData(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** sets user data of the Benders' decomposition cut; user has to free old data in advance! */
SCIP_EXPORT
void SCIPbenderscutSetData(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< new Benders' decomposition cut user data */
   );

/** gets name of the Benders' decomposition cut */
SCIP_EXPORT
const char* SCIPbenderscutGetName(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets description of the Benders' decomposition cut */
SCIP_EXPORT
const char* SCIPbenderscutGetDesc(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets priority of the Benders' decomposition cut */
SCIP_EXPORT
int SCIPbenderscutGetPriority(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets the number of times, the Benders' decomposition cut was called and tried to find a violated cut */
SCIP_EXPORT
SCIP_Longint SCIPbenderscutGetNCalls(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets the number of the cuts found by this Benders' decomposition cut */
SCIP_EXPORT
SCIP_Longint SCIPbenderscutGetNFound(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** is the Benders' decomposition cut initialized? */
SCIP_EXPORT
SCIP_Bool SCIPbenderscutIsInitialized(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets time in seconds used in this Benders' decomposition cut for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPbenderscutGetSetupTime(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** gets time in seconds used in this Benders' decomposition cut */
SCIP_EXPORT
SCIP_Real SCIPbenderscutGetTime(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** returns whether the Benders' cut uses the LP information */
SCIP_EXPORT
SCIP_Bool SCIPbenderscutIsLPCut(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   );

/** sets the enabled flag of the Benders' decomposition cut method */
SCIP_EXPORT
void SCIPbenderscutSetEnabled(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_Bool             enabled             /**< flag to indicate whether the Benders' decomposition cut is enabled */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
