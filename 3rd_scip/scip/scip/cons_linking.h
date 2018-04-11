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

/**@file   cons_linking.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for linking binary variables to an integer variable
 * @author Stefan Heinz
 * @author Jens Schulz
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_LINKING_H__
#define __SCIP_CONS_LINKING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for linking constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
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
 * The constraints handler stores linking constraints between an integer variable and an array of binary variables. Such
 * a linking constraint has the form:
 * \f[
 * y = \sum_{i=1}^n {c_i * x_i}
 * \f]
 * with integer variable \f$ y \f$, binary variables \f$ x_1, \dots, x_n \f$ and offset \f$b \in Q\f$, and
 * with the additional side condition that exactly one binary variable has to be one (set partitioning condition).
 *
 * This constraint can be created only with the integer variable. In this case the binary variables are only created on
 * demand. That is, whenever someone asks for the binary variables. Therefore, such constraints can be used to get a
 * "binary representation" of the domain of the integer variable which will be dynamically created.
 */

/** creates and captures a linking constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             intvar,             /**< integer variable which should be linked */
   SCIP_VAR**            binvars,            /**< binary variables */
   int*                  vals,               /**< coefficients of the binary variables */
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
EXTERN
SCIP_RETCODE SCIPcreateConsBasicLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             intvar,             /**< integer variable which should be linked */
   SCIP_VAR**            binvars,            /**< binary variables, or NULL */
   int*                  vals,               /**< coefficients of the binary variables */
   int                   nbinvars            /**< number of binary variables */
   );


/** checks if for the given integer variable a linking constraint exists */
EXTERN
SCIP_Bool SCIPexistsConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             intvar              /**< integer variable which should be linked */
   );

/** returns the linking constraint belonging the given integer variable or NULL if it does not exist yet */
EXTERN
SCIP_CONS* SCIPgetConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             intvar              /**< integer variable which should be linked */
   );

/** returns the integer variable of the linking constraint */
EXTERN
SCIP_VAR* SCIPgetIntvarLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

/** returns the binary variables of the linking constraint */
EXTERN
SCIP_RETCODE SCIPgetBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR***           binvars,            /**< pointer to store the binary variables array pointer */
   int*                  nbinvars            /**< pointer to store the number of returned binary variables */
   );

/** returns the number of binary variables of the linking constraint */
EXTERN
int SCIPgetNBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

/** returns the coefficients of the binary variables */
EXTERN
int* SCIPgetValsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
