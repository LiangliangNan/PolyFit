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

/**@file   cons_bivariate.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for bivariate nonlinear constraints \f$\textrm{lhs} \leq f(x,y) + c z \leq \textrm{rhs}\f$
 * @author Martin Ballerstein
 * @author Dennis Michaels
 * @author Stefan Vigerske
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_BIVARIATE_H__
#define __SCIP_CONS_BIVARIATE_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for bivariate constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrBivariate(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Bivariate Constraints
 *
 * This constraint handler handles constraints of the form
 * \f[
 *   \textrm{lhs} \leq f(x,y) + c z \leq \textrm{rhs}
 * \f]
 * for a bivariate nonlinear function \f$f(x,y)\f$ (given as expression tree) that has
 * a fixed convexity behaviour, that is, \f$f(x,y)\f$ has to be either jointly convex in \f$(x,y)\f$,
 * or convex in \f$x\f$ and concave in \f$y\f$,
 * or convex in \f$x\f$ and convex in \f$y\f$, but indefinite w.r.t. \f$(x,y)\f$.
 * See also
 *
 * @par
 * Martin Ballerstein, Dennis Michaels, and Stefan Vigerske@n
 * Linear Underestimators for bivariate functions with a fixed convexity behavior@n
 * ZIB Report 13-02, 2013. http://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/1764
 *
 * @{
 */

typedef enum
{
   SCIP_BIVAR_ALLCONVEX          = 0,        /* f(x,y) is convex */
   SCIP_BIVAR_1CONVEX_INDEFINITE = 1,        /* f(x,y) is 1-convex and indefinite */
   SCIP_BIVAR_CONVEX_CONCAVE     = 2,        /* f(x,y) is convex in x and concave in y */
   SCIP_BIVAR_UNKNOWN            = 3         /* unknown */
} SCIP_BIVAR_CONVEXITY;

/** creates and captures a bivariate constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPRTREE*        f,                  /**< expression tree specifying bivariate function f(x,y) */
   SCIP_BIVAR_CONVEXITY  convextype,         /**< kind of convexity of f(x,y) */
   SCIP_VAR*             z,                  /**< linear variable in constraint */
   SCIP_Real             zcoef,              /**< coefficient of linear variable */
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
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures an absolute power constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsBivariate(); all flags can be set via SCIPconsSetFLAGNAME-methods in cons.h
 *
 *  @see SCIPcreateConsBivariate() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPRTREE*        f,                  /**< expression tree specifying bivariate function f(x,y) */
   SCIP_BIVAR_CONVEXITY  convextype,         /**< kind of convexity of f(x,y) */
   SCIP_VAR*             z,                  /**< linear variable in constraint */
   SCIP_Real             zcoef,              /**< coefficient of linear variable */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** gets the linear variable of a bivariate constraint, or NULL if no such variable */
EXTERN
SCIP_VAR* SCIPgetLinearVarBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the coefficients of the linear variable of a bivariate constraint */
EXTERN
SCIP_Real SCIPgetLinearCoefBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the expression tree of a bivariate constraint */
EXTERN
SCIP_EXPRTREE* SCIPgetExprtreeBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the left hand side of a bivariate constraint */
EXTERN
SCIP_Real SCIPgetLhsBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets the right hand side of a bivariate constraint */
EXTERN
SCIP_Real SCIPgetRhsBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
