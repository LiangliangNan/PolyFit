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

/**@file   expr.h
 * @brief  private functions to work with algebraic expressions
 * @ingroup INTERNALAPI
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_EXPR_H_
#define SCIP_EXPR_H_

#include "scip/pub_expr.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_clock.h"
#include "blockmemshell/memory.h"

#ifdef NDEBUG
#include "scip/struct_expr.h"
#include "scip/struct_set.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@name Expression Handler Methods */
/**@{ */

/** create expression handler */
SCIP_RETCODE SCIPexprhdlrCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRHDLR**       exprhdlr,           /**< buffer where to store created expression handler */
   const char*           name,               /**< name of expression handler (must not be NULL) */
   const char*           desc,               /**< description of expression handler (can be NULL) */
   unsigned int          precedence,         /**< precedence of expression operation (used for printing) */
   SCIP_DECL_EXPREVAL((*eval)),              /**< point evaluation callback (must not be NULL) */
   SCIP_EXPRHDLRDATA*    data                /**< data of expression handler (can be NULL) */
   );

/** frees expression handler */
SCIP_RETCODE SCIPexprhdlrFree(
   SCIP_EXPRHDLR**       exprhdlr,           /**< pointer to expression handler to be freed */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** copies the given expression handler to a new scip */
SCIP_RETCODE SCIPexprhdlrCopyInclude(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             targetset           /**< SCIP_SET of SCIP to copy to */
   );

/** initialization of expression handler (resets statistics) */
void SCIPexprhdlrInit(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls the print callback of an expression handler
 *
 * The method prints an expression.
 * It is called while iterating over the expression graph at different stages.
 *
 * @see SCIP_DECL_EXPRPRINT
 */
SCIP_RETCODE SCIPexprhdlrPrintExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRITER_STAGE   stage,              /**< stage of expression iteration */
   int                   currentchild,       /**< index of current child if in stage visitingchild or visitedchild */
   unsigned int          parentprecedence,   /**< precedence of parent */
   FILE*                 file                /**< the file to print to */
   );

/** calls the parse callback of an expression handler
 *
 * The method parses an expression.
 * It should be called when parsing an expression and an operator with the expr handler name is found.
 *
 * @see SCIP_DECL_EXPRPARSE
 */
SCIP_RETCODE SCIPexprhdlrParseExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           string,             /**< string containing expression to be parse */
   const char**          endstring,          /**< buffer to store the position of string after parsing */
   SCIP_EXPR**           expr,               /**< buffer to store the parsed expression */
   SCIP_Bool*            success,            /**< buffer to store whether the parsing was successful or not */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** calls the curvature check callback of an expression handler
 *
 * @see SCIP_DECL_EXPRCURVATURE
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcallExprCurvature() macro */
SCIP_RETCODE SCIPexprhdlrCurvatureExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to check the curvature for */
   SCIP_EXPRCURV         exprcurvature,      /**< desired curvature of this expression */
   SCIP_Bool*            success,            /**< buffer to store whether the desired curvature be obtained */
   SCIP_EXPRCURV*        childcurv           /**< array to store required curvature for each child */
   );

/** calls the monotonicity check callback of an expression handler
 *
 * @see SCIP_DECL_EXPRMONOTONICITY
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcallExprMonotonicity() macro */
SCIP_RETCODE SCIPexprhdlrMonotonicityExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to check the monotonicity for */
   int                   childidx,           /**< index of the considered child expression */
   SCIP_MONOTONE*        result              /**< buffer to store the monotonicity */
   );

/** calls the integrality check callback of an expression handler
 *
 * @see SCIP_DECL_EXPRINTEGRALITY
 */
SCIP_RETCODE SCIPexprhdlrIntegralityExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to check integrality for */
   SCIP_Bool*            isintegral          /**< buffer to store whether expression is integral */
   );

/** calls the hash callback of an expression handler
 *
 * The method hashes an expression by taking the hashes of its children into account.
 *
 * @see SCIP_DECL_EXPRHASH
 */
SCIP_RETCODE SCIPexprhdlrHashExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be hashed */
   unsigned int*         hashkey,            /**< buffer to store the hash value */
   unsigned int*         childrenhashes      /**< array with hash values of children */
   );

/** calls the compare callback of an expression handler
 *
 * The method receives two expressions, expr1 and expr2, and returns
 * - -1 if expr1 < expr2,
 * - 0  if expr1 = expr2,
 * - 1  if expr1 > expr2.
 *
 * @see SCIP_DECL_EXPRCOMPARE
 */
int SCIPexprhdlrCompareExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr1,              /**< first expression in comparison */
   SCIP_EXPR*            expr2               /**< second expression in comparison */
   );

/** calls the evaluation callback of an expression handler
 *
 * The method evaluates an expression by taking the values of its children into account.
 *
 * Further, allows to evaluate w.r.t. given expression and children values instead of those stored in children expressions.
 *
 * @see SCIP_DECL_EXPREVAL
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcallExprEval() macro */
SCIP_RETCODE SCIPexprhdlrEvalExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory, can be NULL if childrenvals is NULL */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_Real*            val,                /**< buffer to store value of expression */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*             sol                 /**< solution that is evaluated (can be NULL) */
   );

/** calls the backward derivative evaluation callback of an expression handler
 *
 * The method should compute the partial derivative of expr w.r.t its child at childidx.
 * That is, it returns
 * \f[
 *   \frac{\partial \text{expr}}{\partial \text{child}_{\text{childidx}}}
 * \f]
 *
 * Further, allows to differentiate w.r.t. given expression and children values instead of those stored in children expressions.
 *
 * @see SCIP_DECL_EXPRBWDIFF
 */
SCIP_RETCODE SCIPexprhdlrBwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory, can be NULL if childrenvals is NULL */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   int                   childidx,           /**< index of the child */
   SCIP_Real*            derivative,         /**< buffer to store the partial derivative w.r.t. the i-th children */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_Real             exprval             /**< value for expression, used only if childrenvals is not NULL */
   );

/** calls the forward differentiation callback of an expression handler
 *
 * @see SCIP_DECL_EXPRFWDIFF
 */
SCIP_RETCODE SCIPexprhdlrFwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   SCIP_Real*            dot,                /**< buffer to store derivative value */
   SCIP_SOL*             direction           /**< direction of the derivative (useful only for var expressions) */
   );

/** calls the evaluation and forward-differentiation callback of an expression handler
 *
 * The method evaluates an expression by taking the values of its children into account.
 * The method differentiates an expression by taking the values and directional derivatives of its children into account.
 *
 * Further, allows to evaluate and differentiate w.r.t. given values for children instead of those stored in children expressions.
 *
 * It probably doesn't make sense to call this function for a variable-expression if sol and/or direction are not given.
 *
 * @see SCIP_DECL_EXPREVAL
 * @see SCIP_DECL_EXPRFWDIFF
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcallExprEvalFwdiff() macro */
SCIP_RETCODE SCIPexprhdlrEvalFwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory, can be NULL if childrenvals is NULL */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_Real*            val,                /**< buffer to store value of expression */
   SCIP_Real*            dot,                /**< buffer to store derivative value */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*             sol,                /**< solution that is evaluated (can be NULL) */
   SCIP_Real*            childrendirs,       /**< directional derivatives for children, or NULL if dot-values stored in children should be used */
   SCIP_SOL*             direction           /**< direction of the derivative (useful only for var expressions, can be NULL if childrendirs is given) */
   );

/** calls the evaluation callback for Hessian directions (backward over forward) of an expression handler
 *
 * @see SCIP_DECL_EXPRBWFWDIFF
 */
SCIP_RETCODE SCIPexprhdlrBwFwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   int                   childidx,           /**< index of the child */
   SCIP_Real*            bardot,             /**< buffer to store derivative value */
   SCIP_SOL*             direction           /**< direction of the derivative (useful only for var expressions) */
   );

/** calls the interval evaluation callback of an expression handler
 *
 * @see SCIP_DECL_EXPRINTEVAL
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcallExprInteval() macro */
SCIP_RETCODE SCIPexprhdlrIntEvalExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_INTERVAL*        interval,           /**< buffer where to store interval */
   SCIP_DECL_EXPR_INTEVALVAR((*intevalvar)), /**< callback to be called when interval-evaluating a variable */
   void*                 intevalvardata      /**< data to be passed to intevalvar callback */
   );

/** calls the estimator callback of an expression handler
 *
 * @see SCIP_DECL_EXPRESTIMATE
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcallExprEstimate() macro */
SCIP_RETCODE SCIPexprhdlrEstimateExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be estimated */
   SCIP_INTERVAL*        localbounds,        /**< current bounds for children */
   SCIP_INTERVAL*        globalbounds,       /**< global bounds for children */
   SCIP_Real*            refpoint,           /**< children values for the reference point where to estimate */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_Real             targetvalue,        /**< a value that the estimator shall exceed, can be +/-infinity */
   SCIP_Real*            coefs,              /**< array to store coefficients of estimator */
   SCIP_Real*            constant,           /**< buffer to store constant part of estimator */
   SCIP_Bool*            islocal,            /**< buffer to store whether estimator is valid locally only */
   SCIP_Bool*            success,            /**< buffer to indicate whether an estimator could be computed */
   SCIP_Bool*            branchcand          /**< array to indicate which children (not) to consider for branching */
   );

/** calls the intitial estimators callback of an expression handler
 *
 * @see SCIP_DECL_EXPRINITESTIMATES
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcallExprInitestimates() macro */
SCIP_RETCODE SCIPexprhdlrInitEstimatesExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be estimated */
   SCIP_INTERVAL*        bounds,             /**< bounds for children */
   SCIP_Bool             overestimate,       /**< whether the expression shall be overestimated or underestimated */
   SCIP_Real*            coefs[SCIP_EXPR_MAXINITESTIMATES], /**< buffer to store coefficients of computed estimators */
   SCIP_Real             constant[SCIP_EXPR_MAXINITESTIMATES], /**< buffer to store constant of computed estimators */
   int*                  nreturned           /**< buffer to store number of estimators that have been computed */
   );

/** calls the simplification callback of an expression handler
 *
 * @see SCIP_DECL_EXPRSIMPLIFY
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPsimplifyExpr() and SCIPexprSimplify() macros */
SCIP_RETCODE SCIPexprhdlrSimplifyExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to simplify */
   SCIP_EXPR**           simplifiedexpr,     /**< buffer to store the simplified expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** calls the reverse propagation callback of an expression handler
 *
 * The method propagates given bounds over the children of an expression.
 *
 * @see SCIP_DECL_EXPRREVERSEPROP
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcallExprReverseprop() macro */
SCIP_RETCODE SCIPexprhdlrReversePropExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to propagate */
   SCIP_INTERVAL         bounds,             /**< the bounds on the expression that should be propagated */
   SCIP_INTERVAL*        childrenbounds,     /**< array to store computed bounds for children, initialized with current activity */
   SCIP_Bool*            infeasible          /**< buffer to store whether a children bounds were propagated to an empty interval */
   );

/**@} */


/**@name Expression Methods */
/**@{ */

/** creates and captures an expression with given expression data and children */
SCIP_RETCODE SCIPexprCreate(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data (expression assumes ownership) */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,            /**< children (can be NULL if nchildren is 0) */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** appends child to the children list of expr */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPappendExprChild() macro */
SCIP_RETCODE SCIPexprAppendChild(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPR*            child               /**< expression to be appended */
   );

/** overwrites/replaces a child of an expressions
 *
 * @note the old child is released and the newchild is captured, unless they are the same (=same pointer)
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPreplaceExprChild() macro */
SCIP_RETCODE SCIPexprReplaceChild(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression where a child is going to be replaced */
   int                   childidx,           /**< index of child being replaced */
   SCIP_EXPR*            newchild            /**< the new child */
   );

/** remove all children of expr */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPremoveExprChildren() macro */
SCIP_RETCODE SCIPexprRemoveChildren(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr                /**< expression */
   );

/** copies an expression including subexpressions
 *
 * @note If copying fails due to an expression handler not being available in the targetscip, then *targetexpr will be set to NULL.
 *
 * For all or some expressions, a mapping to an existing expression can be specified via the mapexpr callback.
 * The mapped expression (including its children) will not be copied in this case and its ownerdata will not be touched.
 * If, however, the mapexpr callback returns NULL for the targetexpr, then the expr will be copied in the usual way.
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPduplicateExpr() macro */
SCIP_RETCODE SCIPexprCopy(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             targetset,          /**< global SCIP settings data structure where target expression will live */
   SCIP_STAT*            targetstat,         /**< dynamic problem statistics in target SCIP */
   BMS_BLKMEM*           targetblkmem,       /**< block memory in target SCIP */
   SCIP_EXPR*            sourceexpr,         /**< expression to be copied */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to copy of source expression */
   SCIP_DECL_EXPR_MAPEXPR((*mapexpr)),       /**< expression mapping function, or NULL for creating new expressions */
   void*                 mapexprdata,        /**< data of expression mapping function */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** duplicates the given expression without its children */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPduplicateExprShallow() macro */
SCIP_RETCODE SCIPexprDuplicateShallow(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store (shallow) duplicate of expr */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** captures an expression (increments usage count) */
void SCIPexprCapture(
   SCIP_EXPR*            expr                /**< expression */
   );

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPreleaseExpr() macro */
SCIP_RETCODE SCIPexprRelease(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr                /**< pointer to expression */
   );

/** returns whether an expression is a variable expression */
SCIP_Bool SCIPexprIsVar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a value expression */
SCIP_Bool SCIPexprIsValue(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a sum expression */
SCIP_Bool SCIPexprIsSum(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a product expression */
SCIP_Bool SCIPexprIsProduct(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a power expression */
SCIP_Bool SCIPexprIsPower(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   );

/** print an expression as info-message */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPprintExpr() macro */
SCIP_RETCODE SCIPexprPrint(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   SCIP_EXPR*            expr                /**< expression to be printed */
   );

/** initializes printing of expressions in dot format to a give FILE* pointer */
SCIP_RETCODE SCIPexprPrintDotInit(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRPRINTDATA**  printdata,          /**< buffer to store dot printing data */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   SCIP_EXPRPRINT_WHAT   whattoprint         /**< info on what to print for each expression */
   );

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_RETCODE SCIPexprPrintDotInit2(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRPRINTDATA**  printdata,          /**< buffer to store dot printing data */
   const char*           filename,           /**< name of file to print to */
   SCIP_EXPRPRINT_WHAT   whattoprint         /**< info on what to print for each expression */
   );

/** main part of printing an expression in dot format */
SCIP_RETCODE SCIPexprPrintDot(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPRPRINTDATA*   printdata,          /**< data as initialized by \ref SCIPprintExprDotInit() */
   SCIP_EXPR*            expr                /**< expression to be printed */
   );

/** finishes printing of expressions in dot format */
SCIP_RETCODE SCIPexprPrintDotFinal(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRPRINTDATA**  printdata           /**< buffer where dot printing data has been stored */
   );

/** prints structure of an expression a la Maple's dismantle */
SCIP_RETCODE SCIPexprDismantle(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   SCIP_EXPR*            expr                /**< expression to dismantle */
   );

/** evaluate an expression in a point
 *
 * Iterates over expressions to also evaluate children, if necessary.
 * Value can be received via SCIPexprGetEvalValue().
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 *
 * If a nonzero \p soltag is passed, then only (sub)expressions are
 * reevaluated that have a different solution tag. If a soltag of 0
 * is passed, then subexpressions are always reevaluated.
 * The tag is stored together with the value and can be received via
 * SCIPexprGetEvalTag().
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPevalExpr() macro */
SCIP_RETCODE SCIPexprEval(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** evaluates gradient of an expression for a given point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPevalExprGradient() macro */
SCIP_RETCODE SCIPexprEvalGradient(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr,           /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated (NULL for the current LP solution) */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** evaluates Hessian-vector product of an expression for a given point and direction
 *
 * Evaluates children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffGradientDirNonlinear()
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPevalExprHessianDir() macro */
SCIP_RETCODE SCIPexprEvalHessianDir(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr,           /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated (NULL for the current LP solution) */
   SCIP_Longint          soltag,             /**< tag that uniquely identifies the solution (with its values), or 0. */
   SCIP_SOL*             direction           /**< direction */
   );

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is no longer uptodate.
 * If the expr owner provided a evalactivity-callback, then call this.
 * Otherwise, loop over descendants and compare activitytag with stat's domchgcount, i.e.,
 * whether some bound was changed since last evaluation, to check whether exprhdlrs INTEVAL should be called.
 *
 * @note If expression is set to be integral, then activities are tightened to integral values.
 *   Thus, ensure that the integrality information is valid (if set to TRUE; the default (FALSE) is always ok).
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPevalExprActivity() macro */
SCIP_RETCODE SCIPexprEvalActivity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr            /**< expression */
   );

/** compare expressions
 *
 * @return -1, 0 or 1 if expr1 <, =, > expr2, respectively
 * @note The given expressions are assumed to be simplified.
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcompareExpr() macro */
int SCIPexprCompare(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr1,              /**< first expression */
   SCIP_EXPR*            expr2               /**< second expression */
   );

/** simplifies an expression
 *
 * @see SCIPsimplifyExpr
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPsimplifyExpr() macro */
SCIP_RETCODE SCIPexprSimplify(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr,           /**< expression to be simplified */
   SCIP_EXPR**           simplified,         /**< buffer to store simplified expression */
   SCIP_Bool*            changed,            /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*            infeasible,         /**< buffer to store whether infeasibility has been detected */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

#ifdef NDEBUG
#define SCIPexprCapture(expr) ++(expr)->nuses
#define SCIPexprIsVar(set, expr)     ((expr)->exprhdlr == (set)->exprhdlrvar)
#define SCIPexprIsValue(set, expr)   ((expr)->exprhdlr == (set)->exprhdlrval)
#define SCIPexprIsSum(set, expr)     ((expr)->exprhdlr == (set)->exprhdlrsum)
#define SCIPexprIsProduct(set, expr) ((expr)->exprhdlr == (set)->exprhdlrproduct)
#define SCIPexprIsPower(set, expr)   ((expr)->exprhdlr == (set)->exprhdlrpow)
#endif

/**@} */

/**@name Expression Iterator Methods */
/**@{ */

/** creates an expression iterator */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcreateExpriter() macro */
SCIP_RETCODE SCIPexpriterCreate(
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRITER**       iterator            /**< buffer to store expression iterator */
   );

/** frees an expression iterator */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPfreeExpriter() macro */
void SCIPexpriterFree(
   SCIP_EXPRITER**       iterator            /**< pointer to the expression iterator */
   );

/**@} */


/**@name Quadratic expression functions */
/**@{ */

/** checks whether an expression is quadratic
 *
 * An expression is quadratic if it is either a power expression with exponent 2.0, a product of two expressions,
 * or a sum of terms where at least one is a square or a product of two.
 *
 * Use \ref SCIPexprGetQuadraticData to get data about the representation as quadratic.
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcheckExprQuadratic() macro */
SCIP_RETCODE SCIPexprCheckQuadratic(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool*            isquadratic         /**< buffer to store result */
   );

/** frees information on quadratic representation of an expression
 *
 * Reverts SCIPexprCheckQuadratic().
 * Before doing changes to an expression, it can be useful to call this function.
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPfreeExprQuadratic() macro */
void SCIPexprFreeQuadratic(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr                /**< expression */
   );

/** Checks the curvature of the quadratic function stored in quaddata
 *
 * For this, it builds the matrix Q of quadratic coefficients and computes its eigenvalues using LAPACK.
 * If Q is
 * - semidefinite positive -> curv is set to convex,
 * - semidefinite negative -> curv is set to concave,
 * - otherwise -> curv is set to unknown.
 *
 * If `assumevarfixed` is given and some expressions in quadratic terms correspond to variables present in
 * this hashmap, then the corresponding rows and columns are ignored in the matrix Q.
 */
SCIP_EXPORT  /* need SCIP_EXPORT here, because func is exposed in API via SCIPcomputeExprQuadraticCurvature() macro */
SCIP_RETCODE SCIPexprComputeQuadraticCurvature(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_EXPRCURV*        curv,               /**< pointer to store the curvature of quadratics */
   SCIP_HASHMAP*         assumevarfixed,     /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   SCIP_Bool             storeeigeninfo      /**< whether the eigenvalues and eigenvectors should be stored */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif /* SCIP_EXPR_H_ */
