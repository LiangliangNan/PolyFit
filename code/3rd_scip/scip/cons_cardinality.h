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

/**@file   cons_cardinality.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for cardinality constraints
 * @author Tobias Fischer
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_CARDINALITY_H__
#define __SCIP_CONS_CARDINALITY_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for cardinality constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
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
EXTERN
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
EXTERN
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
EXTERN
SCIP_RETCODE  SCIPchgCardvalCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pointer to hold the created constraint */
   int                   cardval             /**< number of variables allowed to be nonzero */
   );

/** adds variable to cardinality constraint, the position is determined by the given weight */
EXTERN
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
EXTERN
SCIP_RETCODE SCIPappendVarCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_VAR*             indvar              /**< indicator variable to indicate whether variable may be treated as nonzero
                                              *   in cardinality constraint (or NULL if this variable should be created
                                              *   automatically) */
   );

/** gets number of variables in cardinality constraint */
EXTERN
int SCIPgetNVarsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets array of variables in cardinality constraint */
EXTERN
SCIP_VAR** SCIPgetVarsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets cardinality value of cardinality constraint (i.e., right hand side of cardinality constraint) */
EXTERN
int SCIPgetCardvalCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets array of weights in cardinality constraint (or NULL if not existent) */
EXTERN
SCIP_Real* SCIPgetWeightsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
