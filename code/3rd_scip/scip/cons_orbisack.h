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

/**@file   cons_orbisack.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for orbisack constraints
 * @author Christopher Hojny
 *
 * The constraint works on two vectors of variables, which are interpreted as columns of a matrix such that the first
 * column is lexicographically not smaller than the second.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_ORBISACK_H__
#define __SCIP_CONS_ORBISACK_H__

#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the handler for orbisack constraints and includes it in SCIP
 *
 *  @ingroup ConshdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrOrbisack(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Orbisack Constraints
 *
 * @{
 *
 * This constraint handler can be used to handle symmetries in certain 0/1-programs. The principle
 * structure is that some variables can be ordered in matrix form with two columns, such that
 * permuting both columns does not change the validity and objective function value of a solution.
 * That is, there exists a permutation symmetry of the program that permutes the variables of the
 * first and second column row-wise.
 *
 * In more mathematical terms the structure has to be as follows: There are 0/1-variables
 * \f$x_{ij}\f$, \f$i \in \{1, \dots, n\}\f$, \f$j \in \{1, 2\}\f$. Permuting columns of
 * \f$x\f$ does not change the validity and objective function value of any feasible solution.
 */


/** checks whether a given binary solution is feasible for the orbisack */
SCIP_EXPORT
SCIP_RETCODE SCIPcheckSolutionOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to check for feasibility */
   SCIP_VAR**            vars1,              /**< variables of first column */
   SCIP_VAR**            vars2,              /**< variables of second column */
   int                   nrows,              /**< number of rows */
   SCIP_Bool             printreason,        /**< whether reason for infeasibility should be printed */
   SCIP_Bool*            feasible            /**< memory address to store whether sol is feasible */
   );

/** creates and captures a orbisack constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*const*       vars1,              /**< first column matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< second column matrix of variables on which the symmetry acts */
   int                   nrows,              /**< number of rows in variable matrix */
   SCIP_Bool             ispporbisack,       /**< whether the orbisack is a packing/partitioning orbisack */
   SCIP_Bool             isparttype,         /**< whether the orbisack is a partitioning orbisack */
   SCIP_Bool             ismodelcons,        /**< whether the orbisack is a model constraint */
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

/** creates and captures an orbisack constraint in its most basic variant
 *
 *  All constraint flags set to their default values, which can be set afterwards using SCIPsetConsFLAGNAME() in scip.h.
 *
 *  @see SCIPcreateConsOrbisack() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR**            vars2,              /**< second column of matrix of variables on which the symmetry acts */
   int                   nrows,              /**< number of rows in constraint matrix */
   SCIP_Bool             ispporbisack,       /**< whether the orbisack is a packing/partitioning orbisack */
   SCIP_Bool             isparttype,         /**< whether the orbisack is a partitioning orbisack */
   SCIP_Bool             ismodelcons         /**< whether the orbisack is a model constraint */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
