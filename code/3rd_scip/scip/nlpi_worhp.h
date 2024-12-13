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

/**@file    nlpi_worhp.h
 * @brief   Worhp NLP interface
 * @ingroup NLPIS
 * @author  Benjamin Mueller
 * @author  Renke Kuhlmann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_WORHP_H__
#define __SCIP_NLPI_WORHP_H__

#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create solver interface for Worhp solver and includes it into SCIP, if Worhp is available
 *
 * @ingroup NLPIIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlpSolverWorhp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             useip               /**< TRUE for using Interior Point, FALSE for SQP */
   );

/**@addtogroup NLPIS
 *
 * @{
 */

/** gets string that identifies Worhp (version number) */
SCIP_EXPORT
const char* SCIPgetSolverNameWorhp(
   void
   );

/** gets string that describes Worhp (version number) */
SCIP_EXPORT
const char* SCIPgetSolverDescWorhp(
   void
   );

/** returns whether Worhp is available, i.e., whether it has been linked in */
SCIP_EXPORT
SCIP_Bool SCIPisWorhpAvailableWorhp(
   void
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_WORHP_H__ */
