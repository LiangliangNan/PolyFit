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

/**@file   cons_superindicator.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for indicator constraints over arbitrary constraint types
 * @author Ambros Gleixner
 * @author Frederic Pythoud
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SUPERINDICATOR_H__
#define __SCIP_CONS_SUPERINDICATOR_H__


#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_dialog.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif



/*
 *  constraint-specific interface methods
 */

/** creates the handler for superindicator constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrSuperindicator(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Superindicator Constraints
 *
 * @{
 *
 * Superindicator constraints are constraints of the form
 * \f[
 *    x_i = 1 \Rightarrow C(x)
 * \f]
 * where \f$ x_i \f$ is a binary variable and \f$ C(\dot) \f$ a constraint.  The superindicator constraint is satisfied
 * if and only if x_i is zero or C is satisfied.
 */

/** creates and captures a superindicator constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< pointer to the indicator constraint  */
   SCIP_CONS*            slackcons,          /**< constraint corresponding to the handled constraint */
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
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures a superindicator constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsSuperindicator(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsSuperindicator() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< pointer to the indicator constraint  */
   SCIP_CONS*            slackcons           /**< constraint corresponding to the handled constraint */
   );

/** gets binary variable corresponding to the superindicator constraint */
SCIP_EXPORT
SCIP_VAR* SCIPgetBinaryVarSuperindicator(
   SCIP_CONS*            cons                /**< superindicator constraint */
   );

/** gets the slack constraint corresponding to the superindicator constraint */
SCIP_EXPORT
SCIP_CONS* SCIPgetSlackConsSuperindicator(
   SCIP_CONS*            cons                /**< superindicator constraint */
   );



/*
 *  constraint-dependent SCIP methods
 */

/** transforms the current problem into a MinUC problem (minimizing the number of unsatisfied constraints),
 *  a CIP generalization of the MinULR (min. unsatisfied linear relations) problem
 */
SCIP_EXPORT
SCIP_RETCODE SCIPtransformMinUC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            success             /**< pointer to store whether all constraints could be transformed */
   );



/*
 *  constraint-dependent dialog entries
 */

/** dialog execution method for the SCIPtransformMinUC() command */
SCIP_EXPORT
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeMinUC);

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
