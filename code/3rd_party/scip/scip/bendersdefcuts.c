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
 * @ingroup OTHER_CFILES
 * @brief  default cuts for Benders' decomposition
 * @author Stephen J. Maher
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/bendersdefcuts.h"
#include "scip/pub_message.h"

/** includes default Benders' decomposition cuts plugins into SCIP and the associated Benders' decomposition */
SCIP_RETCODE SCIPincludeBendersDefaultCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition struture */
   )
{
   SCIP_CALL( SCIPincludeBenderscutFeas(scip, benders) );
   SCIP_CALL( SCIPincludeBenderscutFeasalt(scip, benders) );
   SCIP_CALL( SCIPincludeBenderscutInt(scip, benders) );
   SCIP_CALL( SCIPincludeBenderscutNogood(scip, benders) );
   SCIP_CALL( SCIPincludeBenderscutOpt(scip, benders) );

   return SCIP_OKAY;
}
