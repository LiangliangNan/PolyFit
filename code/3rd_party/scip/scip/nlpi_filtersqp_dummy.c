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

/**@file    nlpi_filtersqp_dummy.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   dummy filterSQP NLP interface for the case that FilterSQP is not available
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/nlpi_filtersqp.h"

/** create solver interface for filterSQP solver and include it into SCIP, if filterSQP is available */
SCIP_RETCODE SCIPincludeNlpSolverFilterSQP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{  /*lint --e{715}*/
   assert(scip != NULL);

   return SCIP_OKAY;
}

/** gets string that identifies filterSQP */
const char* SCIPgetSolverNameFilterSQP(void)
{
   return "";
}

/** gets string that describes filterSQP */
const char* SCIPgetSolverDescFilterSQP(void)
{
   return "";
}

/** returns whether filterSQP is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisFilterSQPAvailableFilterSQP(void)
{
   return FALSE;
}
