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

/**@file   benderscut_opt.h
 * @ingroup BENDERSCUTS
 * @brief  Generates a standard Benders' decomposition optimality cut
 * @author Stephen J. Maher
 *
 * The classical Benders' decomposition optimality cuts arise from a feasible instance of the Benders' decomposition
 * subproblem. The optimality cuts are an underestimator of the subproblem objective function value. Auxiliary
 * variables, \f$\varphi\f$ are added to the master problem as a lower bound on the subproblem objective function value.
 *
 * Consider a linear Benders' decomposition subproblem that takes the master problem solution \f$\bar{x}\f$ as input:
 * \f[
 * z(\bar{x}) = \min\{d^{T}y : Ty \geq h - H\bar{x}, y \geq 0\}
 * \f]
 * If the subproblem is feasible, and \f$z(\bar{x}) > \varphi\f$ (indicating that the current underestimators are not
 * optimal) then the Benders' decomposition optimality cut can be generated from the optimal dual solution of the
 * subproblem. Let \f$w\f$ be the vector corresponding to the optimal dual solution of the Benders' decomposition
 * subproblem. The resulting cut is:
 * \f[
 * \varphi \geq w^{T}(h - Hx)
 * \f]
 *
 * Next, consider a nonlinear Benders' decomposition subproblem that takes the master problem solution \f$\bar{x}\f$ as input:
 * \f[
 * z(\bar{x}) = \min\{d^{T}y : g(\bar{x},y) \leq 0, y \geq 0\}
 * \f]
 * If the subproblem is feasible, and \f$z(\bar{x}) > \varphi\f$ (indicating that the current underestimators are not
 * optimal) then the Benders' decomposition optimality cut can be generated from the optimal dual solution of the
 * subproblem. Let \f$w\f$ be the vector corresponding to the optimal dual solution of the Benders' decomposition subproblem.
 * The resulting cut is:
 * \f[
 * \varphi \geq z(\bar{x}) + w^{T} \nabla_x g(\bar{x}, y) (x-\bar{x})
 * \f]
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_OPT_H__
#define __SCIP_BENDERSCUT_OPT_H__


#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_benderscut.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_misc.h"
#include "scip/type_nlp.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_exprinterpret.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the optimality Benders' decomposition cut and includes it in SCIP
 *
 *  @ingroup BenderscutIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBenderscutOpt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/** @addtogroup BENDERSCUTS
 * @{
 */

/** Generates a classical Benders' optimality cut using the dual solutions from the subproblem or the input arrays. If
 *  the dual solutions are input as arrays, then a mapping between the array indices and the rows/variables is required.
 *  This method can also be used to generate a feasiblity, is a problem to minimise the infeasibilities has been solved
 *  to generate the dual solutions
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgenerateAndApplyBendersOptCut(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the pricing problem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< the benders' decomposition cut method */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the pricing problem */
   char*                 cutname,            /**< the name for the cut to be generated */
   SCIP_Real             objective,          /**< the objective function of the subproblem */
   SCIP_Real*            primalvals,         /**< the primal solutions for the NLP, can be NULL */
   SCIP_Real*            consdualvals,       /**< dual variables for the constraints, can be NULL */
   SCIP_Real*            varlbdualvals,      /**< the dual variables for the variable lower bounds, can be NULL */
   SCIP_Real*            varubdualvals,      /**< the dual variables for the variable upper bounds, can be NULL */
   SCIP_HASHMAP*         row2idx,            /**< mapping between the row in the subproblem to the index in the dual array, can be NULL */
   SCIP_HASHMAP*         var2idx,            /**< mapping from variable of the subproblem to the index in the dual arrays, can be NULL */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_Bool             addcut,             /**< should the Benders' cut be added as a cut or constraint */
   SCIP_Bool             feasibilitycut,     /**< is this called for the generation of a feasibility cut */
   SCIP_RESULT*          result              /**< the result from solving the subproblems */
   );

/** adds the gradient of a nonlinear row in the current NLP solution of a subproblem to a linear row or constraint in the master problem
 *
 * Only computes gradient w.r.t. master problem variables.
 * Computes also the directional derivative, that is, mult times gradient times solution.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddNlRowGradientBenderscutOpt(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP*                 subproblem,         /**< the SCIP instance of the subproblem */
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_Real             mult,               /**< multiplier */
   SCIP_Real*            primalvals,         /**< the primal solutions for the NLP, can be NULL */
   SCIP_HASHMAP*         var2idx,            /**< mapping from variable of the subproblem to the index in the dual arrays, can be NULL */
   SCIP_Real*            dirderiv,           /**< storage to add directional derivative */
   SCIP_VAR***           vars,               /**< pointer to array of variables in the generated cut with non-zero coefficient */
   SCIP_Real**           vals,               /**< pointer to array of coefficients of the variables in the generated cut */
   int*                  nvars,              /**< the number of variables in the cut */
   int*                  varssize            /**< the number of variables in the array */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
