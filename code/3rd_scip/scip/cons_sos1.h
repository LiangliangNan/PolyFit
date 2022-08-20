/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_sos1.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for SOS type 1 constraints
 * @author Tobias Fischer
 * @author Marc Pfetsch
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_SOS1_H__
#define __SCIP_CONS_SOS1_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for SOS1 constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrSOS1(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Specially Ordered Set (SOS) Type 1 Constraints
 *
 * @{
 *
 * A specially ordered set of type 1 (SOS1) is a sequence of variables such that at most one
 * variable is nonzero. The special case of two variables arises, for instance, from equilibrium or
 * complementary conditions like \f$x \cdot y = 0\f$. Note that it is in principle allowed that a
 * variable appears twice, but it then can be fixed to 0.
 */

/** creates and captures an SOS1 constraint
 *
 *  We set the constraint to not be modifable. If the weights are non
 *  NULL, the variables are ordered according to these weights (in
 *  ascending order).
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
EXTERN
SCIP_RETCODE SCIPcreateConsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights,            /**< weights determining the variable order, or NULL if natural order should be used */
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

/** creates and captures an SOS1 constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME() in scip.h
 *
 *  @see SCIPcreateConsSOS1() for the default constraint flag configuration
 *
 *  @warning Do NOT set the constraint to be modifiable manually, because this might lead
 *  to wrong results as the variable array will not be resorted
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            weights             /**< weights determining the variable order, or NULL if natural order should be used */
   );

/** adds variable to SOS1 constraint, the position is determined by the given weight */
EXTERN
SCIP_RETCODE SCIPaddVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_Real             weight              /**< weight determining position of variable */
   );

/** appends variable to SOS1 constraint */
EXTERN
SCIP_RETCODE SCIPappendVarSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   );

/** gets number of variables in SOS1 constraint */
EXTERN
int SCIPgetNVarsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets array of variables in SOS1 constraint */
EXTERN
SCIP_VAR** SCIPgetVarsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array of weights in SOS1 constraint (or NULL if not existent) */
EXTERN
SCIP_Real* SCIPgetWeightsSOS1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets conflict graph of SOS1 constraints (or NULL if not existent)
 *
 *  @note The conflict graph is globally valid; local changes are not taken into account.
 */
EXTERN
SCIP_DIGRAPH* SCIPgetConflictgraphSOS1(
   SCIP_CONSHDLR*        conshdlr            /**< SOS1 constraint handler */
   );

/** gets number of problem variables that are part of the SOS1 conflict graph */
EXTERN
int SCIPgetNSOS1Vars(
   SCIP_CONSHDLR*        conshdlr            /**< SOS1 constraint handler */
   );

/** returns whether variable is part of the SOS1 conflict graph */
EXTERN
SCIP_Bool SCIPvarIsSOS1(
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_VAR*             var                 /**< variable */
   );

/** returns node of variable in the conflict graph or -1 if variable is not part of the SOS1 conflict graph */
EXTERN
int SCIPvarGetNodeSOS1(
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_VAR*             var                 /**< variable */
   );

/** returns variable that belongs to a given node from the conflict graph */
EXTERN
SCIP_VAR* SCIPnodeGetVarSOS1(
   SCIP_DIGRAPH*         conflictgraph,      /**< conflict graph */
   int                   node                /**< node from the conflict graph */
   );

/** based on solution values of the variables, fixes variables to zero to turn all SOS1 constraints feasible  */
EXTERN
SCIP_RETCODE SCIPmakeSOS1sFeasible(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< SOS1 constraint handler */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool*            changed,            /**< pointer to store whether the solution has been changed */
   SCIP_Bool*            success             /**< pointer to store whether SOS1 constraints have been turned feasible and
                                              *   solution was good enough */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
