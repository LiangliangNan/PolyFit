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

/**@file   pub_nlpi.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for NLP solver interfaces
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_NLPI_H__
#define __SCIP_PUB_NLPI_H__

#include "scip/def.h"
#include "scip/type_nlpi.h"
#include "scip/type_misc.h"

#ifdef NDEBUG
#include "scip/struct_nlpi.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicNLPIInterfaceMethods
 *
 * @{
 */

/** compares two NLPIs w.r.t. their priority */
SCIP_DECL_SORTPTRCOMP(SCIPnlpiComp);

/** gets data of an NLPI */
SCIP_EXPORT
SCIP_NLPIDATA* SCIPnlpiGetData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gets NLP solver name */
SCIP_EXPORT
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gets NLP solver description */
SCIP_EXPORT
const char* SCIPnlpiGetDesc(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gets NLP solver priority */
SCIP_EXPORT
int SCIPnlpiGetPriority(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/**@name Statistics */
/**@{ */

/** gives number of problems created for NLP solver so far */
SCIP_EXPORT
int SCIPnlpiGetNProblems(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gives total time spend in problem creation/modification/freeing */
SCIP_EXPORT
SCIP_Real SCIPnlpiGetProblemTime(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** total number of NLP solves so far */
SCIP_EXPORT
int SCIPnlpiGetNSolves(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gives total time spend in NLP solves (as reported by solver) */
SCIP_EXPORT
SCIP_Real SCIPnlpiGetSolveTime(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gives total time spend in function evaluation during NLP solves
 *
 * If parameter `timing/nlpieval` is off (the default), depending on the NLP solver, this may just return 0.
 */
SCIP_EXPORT
SCIP_Real SCIPnlpiGetEvalTime(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gives total number of iterations spend by NLP solver so far */
SCIP_EXPORT
SCIP_Longint SCIPnlpiGetNIterations(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gives number of times a solve ended with a specific termination status */
SCIP_EXPORT
int SCIPnlpiGetNTermStat(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   SCIP_NLPTERMSTAT      termstatus          /**< the termination status to query for */
   );

/** gives number of times a solve ended with a specific solution status */
SCIP_EXPORT
int SCIPnlpiGetNSolStat(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   SCIP_NLPSOLSTAT       solstatus           /**< the solution status to query for */
   );

/** adds statistics from one NLPI to another */
SCIP_EXPORT
void SCIPnlpiMergeStatistics(
   SCIP_NLPI*            targetnlpi,         /**< NLP interface where to add statistics */
   SCIP_NLPI*            sourcenlpi,         /**< NLP interface from which to add statistics */
   SCIP_Bool             reset               /**< whether to reset statistics in sourcescip */
   );

#ifdef NDEBUG
/* If NDEBUG is defined, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */
#define SCIPnlpiGetData(nlpi)                  (nlpi)->nlpidata
#define SCIPnlpiGetName(nlpi)                  (nlpi)->name
#define SCIPnlpiGetDesc(nlpi)                  (nlpi)->description
#define SCIPnlpiGetPriority(nlpi)              (nlpi)->priority
#define SCIPnlpiGetNProblems(nlpi)             (nlpi)->nproblems
#define SCIPnlpiGetProblemTime(nlpi)           SCIPclockGetTime((nlpi)->problemtime)
#define SCIPnlpiGetNSolves(nlpi)               (nlpi)->nsolves
#define SCIPnlpiGetSolveTime(nlpi)             (nlpi)->solvetime
#define SCIPnlpiGetEvalTime(nlpi)              (nlpi)->evaltime
#define SCIPnlpiGetNIterations(nlpi)           (nlpi)->niter
#define SCIPnlpiGetNTermStat(nlpi, termstatus) (nlpi)->ntermstat[termstatus]
#define SCIPnlpiGetNSolStat(nlpi, solstatus)   (nlpi)->nsolstat[solstatus]
#endif

/**@} */ /* Statistics */

/**@} */ /* PublicNLPIMethods */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_PUB_NLPI_H__ */
