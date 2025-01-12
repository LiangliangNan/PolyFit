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

/**@file   cons_linking.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for linking binary variables to a linking (continuous or integer) variable
 * @author Stefan Heinz
 * @author Jens Schulz
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_LINKING_H__
#define __SCIP_CONS_LINKING_H__


#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for linking constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrLinking(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Linking Constraints
 *
 * @{
 *
 * The constraints handler stores linking constraints between a linking variable (continuous or integer) and an array of binary variables. Such
 * a linking constraint has the form:
 * \f[
 * y = \sum_{i=1}^n {c_i * x_i}
 * \f]
 * with linking variable (continuous or integer) \f$ y \f$, binary variables \f$ x_1, \dots, x_n \f$ and offset \f$b \in Q\f$, and
 * with the additional side condition that exactly one binary variable has to be one (set partitioning condition).
 *
 * This constraint can be created only with the linking variable, if it is an integer variable. In this case the binary variables are only created on
 * demand. That is, whenever someone asks for the binary variables. Therefore, such constraints can be used to get a
 * "binary representation" of the domain of the linking variable which will be dynamically created.
 */

/** creates and captures a linking constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             linkvar,            /**< linking variable (continuous or integer) which should be linked */
   SCIP_VAR**            binvars,            /**< binary variables */
   SCIP_Real*            vals,               /**< coefficients of the binary variables */
   int                   nbinvars,           /**< number of binary starting variables */
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

/** creates and captures a linking constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLinking(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLinking() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             linkvar,            /**< linking variable (continuous or integer) which should be linked */
   SCIP_VAR**            binvars,            /**< binary variables, or NULL */
   SCIP_Real*            vals,               /**< coefficients of the binary variables */
   int                   nbinvars            /**< number of binary variables */
   );


/** checks if for the given linking variable (continuous or integer) a linking constraint exists */
SCIP_EXPORT
SCIP_Bool SCIPexistsConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             linkvar             /**< linking variable (continuous or integer) which should be linked */
   );

/** returns the linking constraint belonging the given linking variable (continuous or integer) or NULL if it does not exist yet */
SCIP_EXPORT
SCIP_CONS* SCIPgetConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             linkvar             /**< linking variable (continuous or integer) which should be linked */
   );

/** returns the linking variable (continuous or integer) of the linking constraint */
SCIP_EXPORT
SCIP_VAR* SCIPgetLinkvarLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

/** returns the binary variables of the linking constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPgetBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR***           binvars,            /**< pointer to store the binary variables array pointer */
   int*                  nbinvars            /**< pointer to store the number of returned binary variables */
   );

/** returns the number of binary variables of the linking constraint */
SCIP_EXPORT
int SCIPgetNBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

/** returns the coefficients of the binary variables */
SCIP_EXPORT
SCIP_Real* SCIPgetValsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

/** return all binary variable information of the linking constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPgetBinvarsDataLinking(
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR***           binvars,            /**< pointer to store binary variables, or NULL */
   SCIP_Real**           vals,               /**< pointer to store the binary coefficients, or NULL */
   int*                  nbinvars            /**< pointer to store the number of binary variables, or NULL */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
