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

/**@file   cons_pseudoboolean.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for pseudoboolean constraints
 * @author Stefan Heinz
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_PSEUDOBOOLEAN_H__
#define __SCIP_CONS_PSEUDOBOOLEAN_H__


#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ARTIFICIALVARNAMEPREFIX "andresultant_"



/** creates the handler for pseudoboolean constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrPseudoboolean(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Pseudoboolean Constraints
 *
 * @{
 *
 * The constraint handler deals with pseudo boolean constraints. These are constraints of the form
 * \f[
 * \mbox{lhs} \leq \sum_{k=0}^m c_k \cdot x_k  +  \sum_{i=0}^n c_i \cdot \prod_{j \in I_i} x_j \leq \mbox{rhs}
 * \f]
 * where all \f$x\f$ are binary.
 */

/** solution status after solving LP */
enum SCIP_LinearConsType
{
   SCIP_LINEARCONSTYPE_INVALIDCONS = -1,     /**< this is no valid linear constraint type */
   SCIP_LINEARCONSTYPE_LINEAR      =  0,     /**< this is the common linear constraint */
   SCIP_LINEARCONSTYPE_LOGICOR     =  1,     /**< this is a logicor constraint */
   SCIP_LINEARCONSTYPE_KNAPSACK    =  2,     /**< this is a knapsack constraint */
#ifndef WITHEQKNAPSACK
   SCIP_LINEARCONSTYPE_SETPPC      =  3      /**< this is a setppc constraint */
#else
   SCIP_LINEARCONSTYPE_SETPPC      =  3,     /**< this is a setppc constraint */
   SCIP_LINEARCONSTYPE_EQKNAPSACK  =  4      /**< this is a equality knapsack constraint */
#endif
};
typedef enum SCIP_LinearConsType SCIP_LINEARCONSTYPE;

/** creates and captures a pseudoboolean constraint, with given linear and and-constraints
 *
 *  @note intvar must currently be NULL
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsPseudobooleanWithConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONS*            lincons,            /**< associated linear constraint */
   SCIP_LINEARCONSTYPE   linconstype,        /**< linear constraint type of associated linear constraint */
   SCIP_CONS**           andconss,           /**< associated and-constraints */
   SCIP_Real*            andcoefs,           /**< associated coefficients of and-constraints */
   int                   nandconss,          /**< number of associated and-constraints */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< an artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant
                                              *   constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant
                                              *   constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching
                                              *   constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if
                                              *   pricing adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which are seperated as
                                              *   constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user
                                              *   cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent
                                              *   node data. */
   );

/** creates and captures a pseudoboolean constraint
 *
 *  @note linear and nonlinear terms can be added using SCIPaddCoefPseudoboolean() and SCIPaddTermPseudoboolean(),
 *        respectively
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 *
 *  @note intvar must currently be NULL
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            linvars,            /**< variables of the linear part, or NULL */
   int                   nlinvars,           /**< number of variables of the linear part */
   SCIP_Real*            linvals,            /**< coefficients of linear part, or NULL */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   SCIP_Real*            termvals,           /**< coefficients of nonlinear parts, or NULL */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< an artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures a pseudoboolean constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME() in scip.h
 *
 *  @see SCIPcreateConsPseudoboolean() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 *
 *  @note intvar must currently be NULL
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicPseudoboolean(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            linvars,            /**< variables of the linear part, or NULL */
   int                   nlinvars,           /**< number of variables of the linear part */
   SCIP_Real*            linvals,            /**< coefficients of linear part, or NULL */
   SCIP_VAR***           terms,              /**< nonlinear terms of variables, or NULL */
   int                   nterms,             /**< number of terms of variables of nonlinear term */
   int*                  ntermvars,          /**< number of variables in nonlinear terms, or NULL */
   SCIP_Real*            termvals,           /**< coefficients of nonlinear parts, or NULL */
   SCIP_VAR*             indvar,             /**< indicator variable if it's a soft constraint, or NULL */
   SCIP_Real             weight,             /**< weight of the soft constraint, if it is one */
   SCIP_Bool             issoftcons,         /**< is this a soft constraint */
   SCIP_VAR*             intvar,             /**< a artificial variable which was added only for the objective function,
                                              *   if this variable is not NULL this constraint (without this integer
                                              *   variable) describes the objective function */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** adds linear term pseudo boolean constraint (if it is not zero)
 *
 * @note you can only add a coefficient if the special type of linear constraint won't changed
 *
 * @todo if adding a coefficient would change the type of the special linear constraint, we need to erase it and
 *       create a new linear constraint
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddCoefPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_VAR* const       var,                /**< variable of constraint entry */
   SCIP_Real const       val                 /**< coefficient of constraint entry */
   );

/** adds nonlinear term to pseudo boolean constraint (if it is not zero)
 *
 * @note you can only add a coefficient if the special type of linear constraint won't changed
 *
 * @todo if adding a coefficient would change the type of the special linear constraint, we need to erase it and
 *       create a new linear constraint
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddTermPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_VAR**const       vars,               /**< variables of the nonlinear term */
   int const             nvars,              /**< number of variables of the nonlinear term */
   SCIP_Real const       val                 /**< coefficient of constraint entry */
   );

/** gets indicator variable of pseudoboolean constraint, or NULL if there is no */
SCIP_EXPORT
SCIP_VAR* SCIPgetIndVarPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   );

/** gets linear constraint of pseudoboolean constraint */
SCIP_EXPORT
SCIP_CONS* SCIPgetLinearConsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   );

/** gets type of linear constraint of pseudoboolean constraint */
SCIP_EXPORT
SCIP_LINEARCONSTYPE SCIPgetLinearConsTypePseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   );

/** gets number of linear variables without artificial terms variables of pseudoboolean constraint */
SCIP_EXPORT
int SCIPgetNLinVarsWithoutAndPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   );

/** gets linear constraint of pseudoboolean constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPgetLinDatasWithoutAndPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_VAR**const       linvars,            /**< array to store and-constraints */
   SCIP_Real*const       lincoefs,           /**< array to store and-coefficients */
   int*const             nlinvars            /**< pointer to store the required array size for and-constraints, have to
                                              *   be initialized with size of given array */
   );

/** gets and-constraints of pseudoboolean constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPgetAndDatasPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_CONS**const      andconss,           /**< array to store and-constraints */
   SCIP_Real*const       andcoefs,           /**< array to store and-coefficients */
   int*const             nandconss           /**< pointer to store the required array size for and-constraints, have to
                                              *   be initialized with size of given array */
   );

/** gets number of and constraints of pseudoboolean constraint */
SCIP_EXPORT
int SCIPgetNAndsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   );

/** changes left hand side of pseudoboolean constraint
 *
 * @note you can only change the left hand side if the special type of linear constraint won't changed
 *
 * @todo if changing the left hand side would change the type of the special linear constraint, we need to erase it
 *       and create a new linear constraint
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgLhsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_Real const       lhs                 /**< new left hand side */
   );

/** changes right hand side of pseudoboolean constraint
 *
 * @note you can only change the right hand side if the special type of linear constraint won't changed
 *
 * @todo if changing the right hand side would change the type of the special linear constraint, we need to erase it
 *       and create a new linear constraint
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgRhsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons,               /**< pseudoboolean constraint */
   SCIP_Real const       rhs                 /**< new right hand side */
   );

/** get left hand side of pseudoboolean constraint */
SCIP_EXPORT
SCIP_Real SCIPgetLhsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   );

/** get right hand side of pseudoboolean constraint */
SCIP_EXPORT
SCIP_Real SCIPgetRhsPseudoboolean(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_CONS*const       cons                /**< pseudoboolean constraint */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
