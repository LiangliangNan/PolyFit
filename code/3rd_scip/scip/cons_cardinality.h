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

/**@file   cons_cardinality.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for cardinality constraints
 * @author Tobias Fischer
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_CARDINALITY_H__
#define __SCIP_CONS_CARDINALITY_H__


#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for cardinality constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrCardinality(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Cardinality Constraints
 *
 * @{
 *
 * This constraint handler handles cardinality constraints of the form
 * \f[
 *   |\mbox{supp}(x)| \leq b
 * \f]
 * with integer right-hand side \f$b\f$. Here, \f$|\mbox{supp}(x)|\f$ denotes the number of nonzero entries of the
 * vector \f$x\f$.
 *
 * Cardinality constraints generalize special ordered set of type one (SOS1) constraints in which \f$b = 1\f$.
 */

/** creates and captures an cardinality constraint
 *
 *  We set the constraint to not be modifable. If the weights are non
 *  NULL, the variables are ordered according to these weights (in
 *  ascending order).
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   int                   cardval,            /**< number of variables allowed to be nonzero */
   SCIP_VAR**            indvars,            /**< indicator variables to indicate which variables may be treated as nonzero
                                              *   in cardinality constraint, or NULL if indicator variables should be
                                              *   created automatically */
   SCIP_Real*            weights,            /**< weights determining the variable order, or NULL if variables should be
                                              *   ordered in the same way they were added to the constraint */
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

/** creates and captures an cardinality constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME() in scip.h
 *
 *  @see SCIPcreateConsCardinality() for the default constraint flag configuration
 *
 *  @warning Do NOT set the constraint to be modifiable manually, because this might lead
 *  to wrong results as the variable array will not be resorted
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   int                   cardval,            /**< number of variables allowed to be nonzero */
   SCIP_VAR**            indvars,            /**< indicator variables to indicate which variables may be treated as nonzero
                                              *   in cardinality constraint, or NULL if indicator variables should be
                                              *   created automatically */
   SCIP_Real*            weights             /**< weights determining the variable order, or NULL if variables should be
                                              *   ordered in the same way they were added to the constraint */
   );

/** changes cardinality value of cardinality constraint (i.e., right hand side of cardinality constraint) */
SCIP_EXPORT
SCIP_RETCODE  SCIPchgCardvalCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pointer to hold the created constraint */
   int                   cardval             /**< number of variables allowed to be nonzero */
   );

/** adds variable to cardinality constraint, the position is determined by the given weight */
SCIP_EXPORT
SCIP_RETCODE SCIPaddVarCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_VAR*             indvar,             /**< indicator variable to indicate whether variable may be treated as nonzero
                                              *   in cardinality constraint (or NULL if this variable should be created
                                              *   automatically) */
   SCIP_Real             weight              /**< weight determining position of variable */
   );

/** appends variable to cardinality constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPappendVarCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_VAR*             indvar              /**< indicator variable to indicate whether variable may be treated as nonzero
                                              *   in cardinality constraint (or NULL if this variable should be created
                                              *   automatically) */
   );

/** gets number of variables in cardinality constraint */
SCIP_EXPORT
int SCIPgetNVarsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets array of variables in cardinality constraint */
SCIP_EXPORT
SCIP_VAR** SCIPgetVarsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets cardinality value of cardinality constraint (i.e., right hand side of cardinality constraint) */
SCIP_EXPORT
int SCIPgetCardvalCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array of weights in cardinality constraint (or NULL if not existent) */
SCIP_EXPORT
SCIP_Real* SCIPgetWeightsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
