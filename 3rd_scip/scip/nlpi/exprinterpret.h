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

/**@file   exprinterpret.h
 * @brief  methods to interpret (evaluate) an expression tree "fast"
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 * Realized similar to LPI: one implementation of an interpreter is linked in.
 */

/* @todo product Gradient times vector
   @todo product Hessian times vector
   @todo product Hessian of Lagrangian times vector
   @todo sparse Hessian of expression tree
   @todo sparse Hessian of Lagrangian (sets of expression trees and quadratic parts)?
*/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPRINTERPRET_H__
#define __SCIP_EXPRINTERPRET_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "nlpi/type_expr.h"
#include "nlpi/type_exprinterpret.h"
#include "scip/intervalarith.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup EXPRINTS
 * @{
 */

/** gets name and version of expression interpreter */
EXTERN
const char* SCIPexprintGetName(void);

/** gets descriptive text of expression interpreter */
EXTERN
const char* SCIPexprintGetDesc(void);

/** gets capabilities of expression interpreter (using bitflags) */
EXTERN
SCIP_EXPRINTCAPABILITY SCIPexprintGetCapability(void);

/** creates an expression interpreter object */
EXTERN
SCIP_RETCODE SCIPexprintCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
   );

/** frees an expression interpreter object */
EXTERN
SCIP_RETCODE SCIPexprintFree(
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
   );

/** compiles an expression tree and stores compiled data in expression tree */
EXTERN
SCIP_RETCODE SCIPexprintCompile(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** gives the capability to evaluate an expression by the expression interpreter
 *
 * In cases of user-given expressions, higher order derivatives may not be available for the user-expression,
 * even if the expression interpreter could handle these. This method allows to recognize that, e.g., the
 * Hessian for an expression is not available because it contains a user expression that does not provide
 * Hessians.
 */
EXTERN
SCIP_EXPRINTCAPABILITY SCIPexprintGetExprtreeCapability(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** frees interpreter data */
EXTERN
SCIP_RETCODE SCIPexprintFreeData(
   SCIP_EXPRINTDATA**    interpreterdata     /**< interpreter data that should freed */
   );

/** notify expression interpreter that a new parameterization is used
 * this probably causes retaping by AD algorithms
 */
EXTERN
SCIP_RETCODE SCIPexprintNewParametrization(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   );

/** evaluates an expression tree */
EXTERN
SCIP_RETCODE SCIPexprintEval(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Real*            val                 /**< buffer to store value of expression */
   );

/** evaluates an expression tree on intervals */
EXTERN
SCIP_RETCODE SCIPexprintEvalInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables */
   SCIP_INTERVAL*        val                 /**< buffer to store interval value of expression */
   );

/** computes value and gradient of an expression tree */
EXTERN
SCIP_RETCODE SCIPexprintGrad(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to a point evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store expression value */
   SCIP_Real*            gradient            /**< buffer to store expression gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
   );

/** computes interval value and interval gradient of an expression tree */
EXTERN
SCIP_RETCODE SCIPexprintGradInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable interval values changed since last call to an interval evaluation routine? */
   SCIP_INTERVAL*        val,                /**< buffer to store expression interval value */
   SCIP_INTERVAL*        gradient            /**< buffer to store expression interval gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
   );

/** gives sparsity pattern of hessian
 *
 * NOTE: this function might be replaced later by something nicer.
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 */
EXTERN
SCIP_RETCODE SCIPexprintHessianSparsityDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Bool*            sparsity            /**< buffer to store sparsity pattern of Hessian, sparsity[i+n*j] indicates whether entry (i,j) is nonzero in the hessian */
   );

/** computes value and dense hessian of an expression tree
 *
 *  The full hessian is computed (lower left and upper right triangle).
 */
EXTERN
SCIP_RETCODE SCIPexprintHessianDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   SCIP_Real*            hessian             /**< buffer to store hessian values, need to have size at least n*n */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPRINTERPRET_H__ */
