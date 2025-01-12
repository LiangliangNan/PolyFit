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

/**@file   cons_orbitope.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for (partitioning/packing/full) orbitope constraints w.r.t. the full symmetric group
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_ORBITOPE_H__
#define __SCIP_CONS_ORBITOPE_H__

#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "symmetry/type_symmetry.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for orbitope constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrOrbitope(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Orbitope Constraints
 *
 * @{
 *
 * This constraint handler can be used to handle symmetries in certain 0/1-programs. The principle
 * structure is that some variables can be ordered in matrix form, such that permuting columns does
 * not change the validity and objective function value of a solution. That is, the symmetry group
 * of the program contains the full symmetric group obtained by permuting the columns of this
 * matrix. These symmetries can be handled by so-called full orbitopes.
 *
 * Moreover, if the variables in each row are contained in set packing or partitioning
 * constraint, these symmetries can be handled by specialized packing or partitioning orbitopes.
 *
 * In more mathematical terms the structure has to be as follows: There are 0/1-variables
 * \f$x_{ij}\f$, \f$i \in \{1, \dots, p\}\f$, \f$j \in \{1, \dots, q\}\f$. The variables may be coupled
 * through set packing or partitioning constraints:
 * \f[
 *    \sum_{j = 1}^q x_{ij} \leq 1  \quad \mbox{or} \quad \sum_{j = 1}^q x_{ij} = 1 \quad \mbox{for all }i = 1, \ldots, p.
 * \f]
 * Permuting columns of \f$x\f$ does not change the validity and objective function value of any feasible solution.
 *
 * We distinguish whether an orbitope is a model constraint or not. If it is a model constraint, then
 * its information are copied to subSCIPs. Otherwise, the constraint was added just for the purpose of
 * symmetry handling and we do not copy its information to subSCIPs.
 */

/** creates and captures a orbitope constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint */
   int                   nspcons,            /**< number of set partitioning/packing constraints  <=> p */
   int                   nblocks,            /**< number of symmetric variable blocks             <=> q */
   SCIP_Bool             usedynamicprop,     /**< whether dynamic propagation should be used */
   SCIP_Bool             mayinteract,        /**< whether symmetries corresponding to orbitope might interact
                                              *   with symmetries handled by other routines */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons,        /**< whether the orbitope is a model constraint */
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

/** creates and captures an orbitope constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME() in scip.h
 *
 *  @see SCIPcreateConsOrbitope() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_ORBITOPETYPE     orbitopetype,       /**< type of orbitope constraint */
   int                   nspcons,            /**< number of set partitioning/packing constraints  <=> p */
   int                   nblocks,            /**< number of symmetric variable blocks             <=> q */
   SCIP_Bool             usedynamicprop,     /**< whether dynamic propagation should be used */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
   SCIP_Bool             ismodelcons,        /**< whether the orbitope is a model constraint */
   SCIP_Bool             mayinteract         /**< whether symmetries corresponding to orbitope might interact
                                              *   with symmetries handled by other routines */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
