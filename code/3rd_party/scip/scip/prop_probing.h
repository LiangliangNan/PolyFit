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

/**@file   prop_probing.h
 * @ingroup PROPAGATORS
 * @brief  probing propagator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_PROBING_H__
#define __SCIP_PROP_PROBING_H__

#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the probing propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropProbing(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PROPAGATORS
  *
  * @{
  */

/** applies and evaluates probing of a single variable in the given direction and bound */
SCIP_EXPORT
SCIP_RETCODE SCIPapplyProbingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   probingpos,         /**< variable number to apply probing on */
   SCIP_BOUNDTYPE        boundtype,          /**< which bound should be changed */
   SCIP_Real             bound,              /**< which bound should be set */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Real*            impllbs,            /**< array to store lower bounds after applying implications and cliques */
   SCIP_Real*            implubs,            /**< array to store upper bounds after applying implications and cliques */
   SCIP_Real*            proplbs,            /**< array to store lower bounds after full propagation */
   SCIP_Real*            propubs,            /**< array to store upper bounds after full propagation */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   );

/** analyses boundchanges resulting from probing on a variable and performs deduced fixations, aggregations, and domain tightenings
 *
 *  Given a variable probingvar with domain [l,u] and bound tightening results from reducing the
 *  domain once to [l,leftub] and once to [rightlb,u], the method computes and applies resulting
 *  variable fixations, aggregations, implications, and bound changes. Variable probingvar does not
 *  need to be binary.  The whole domain of probingvar need to be covered by the left and right
 *  branches, i.e., we assume leftub >= rightlb for continuous variables or floor(leftub) >=
 *  ceil(rightlb)-1 for discrete variables.  Bounds after applying implications and cliques do not
 *  need to be provided, but if they are omitted and probingvar is a binary variable, then already
 *  existing implications may be added.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPanalyzeDeductionsProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             probingvar,         /**< the probing variable */
   SCIP_Real             leftub,             /**< upper bound of probing variable in left branch */
   SCIP_Real             rightlb,            /**< lower bound of probing variable in right branch */
   int                   nvars,              /**< number of variables which bound changes should be analyzed */
   SCIP_VAR**            vars,               /**< variables which bound changes should be analyzed */
   SCIP_Real*            leftimpllbs,        /**< lower bounds after applying implications and cliques in left branch, or NULL */
   SCIP_Real*            leftimplubs,        /**< upper bounds after applying implications and cliques in left branch, or NULL */
   SCIP_Real*            leftproplbs,        /**< lower bounds after applying domain propagation in left branch */
   SCIP_Real*            leftpropubs,        /**< upper bounds after applying domain propagation in left branch */
   SCIP_Real*            rightimpllbs,       /**< lower bounds after applying implications and cliques in right branch, or NULL */
   SCIP_Real*            rightimplubs,       /**< upper bounds after applying implications and cliques in right branch, or NULL */
   SCIP_Real*            rightproplbs,       /**< lower bounds after applying domain propagation in right branch */
   SCIP_Real*            rightpropubs,       /**< upper bounds after applying domain propagation in right branch */
   int*                  nfixedvars,         /**< pointer to counter which is increased by the number of deduced variable fixations */
   int*                  naggrvars,          /**< pointer to counter which is increased by the number of deduced variable aggregations */
   int*                  nimplications,      /**< pointer to counter which is increased by the number of deduced implications */
   int*                  nchgbds,            /**< pointer to counter which is increased by the number of deduced bound tightenings */
   SCIP_Bool*            cutoff              /**< buffer to store whether a cutoff is detected */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
