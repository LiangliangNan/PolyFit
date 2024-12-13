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

/**@file   bendersdefcuts.c
 * @brief  default cuts for Benders' decomposition
 * @author Stephen J. Maher
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSDEFCUTS_H__
#define __SCIP_BENDERSDEFCUTS_H__

/* include header files here, such that the user only has to include bendersdefcuts.h */
#include "scip/benderscut_feas.h"
#include "scip/benderscut_feasalt.h"
#include "scip/benderscut_int.h"
#include "scip/benderscut_nogood.h"
#include "scip/benderscut_opt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes default Benders' decomposition cuts plugins into SCIP and the associated Benders' decomposition */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBendersDefaultCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition struture */
   );

#ifdef __cplusplus
}
#endif

#endif
