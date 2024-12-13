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

/**@file   pub_expr.h
 * @ingroup PUBLICCOREAPI
 * @brief  public functions to work with algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_PUB_EXPR_H_
#define SCIP_PUB_EXPR_H_

#include "scip/def.h"
#include "scip/type_expr.h"
#include "scip/type_misc.h"

#ifdef NDEBUG
#include "scip/struct_expr.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicExprHandlerMethods
 * @{
 */

/** set the expression handler callbacks to copy and free an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetCopyFreeHdlr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOPYHDLR((*copyhdlr)),      /**< handler copy callback (can be NULL) */
   SCIP_DECL_EXPRFREEHDLR((*freehdlr))       /**< handler free callback (can be NULL) */
);

/** set the expression handler callbacks to copy and free expression data */
SCIP_EXPORT
void SCIPexprhdlrSetCopyFreeData(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOPYDATA((*copydata)),      /**< expression data copy callback (can be NULL for expressions without data) */
   SCIP_DECL_EXPRFREEDATA((*freedata))       /**< expression data free callback (can be NULL if data does not need to be freed) */
);

/** set the print callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetPrint(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRPRINT((*print))             /**< print callback (can be NULL) */
);

/** set the parse callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetParse(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRPARSE((*parse))             /**< parse callback (can be NULL) */
);

/** set the curvature detection callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetCurvature(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCURVATURE((*curvature))     /**< curvature detection callback (can be NULL) */
);

/** set the monotonicity detection callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetMonotonicity(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRMONOTONICITY((*monotonicity)) /**< monotonicity detection callback (can be NULL) */
);

/** set the integrality detection callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetIntegrality(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINTEGRALITY((*integrality)) /**< integrality detection callback (can be NULL) */
);

/** set the hash callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetHash(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRHASH((*hash))               /**< hash callback (can be NULL) */
);

/** set the compare callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetCompare(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOMPARE((*compare))         /**< compare callback (can be NULL) */
);

/** set differentiation callbacks of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetDiff(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRBWDIFF((*bwdiff)),          /**< backward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRFWDIFF((*fwdiff)),          /**< forward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRBWFWDIFF((*bwfwdiff))       /**< backward-forward derivative evaluation callback (can be NULL) */
);

/** set the interval evaluation callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetIntEval(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINTEVAL((*inteval))         /**< interval evaluation callback (can be NULL) */
);

/** set the simplify callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetSimplify(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRSIMPLIFY((*simplify))       /**< simplify callback (can be NULL) */
);

/** set the reverse propagation callback of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetReverseProp(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
);

/** set the estimation callbacks of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrSetEstimate(
   SCIP_EXPRHDLR*        exprhdlr,                /**< expression handler */
   SCIP_DECL_EXPRINITESTIMATES((*initestimates)), /**< initial estimators callback (can be NULL) */
   SCIP_DECL_EXPRESTIMATE((*estimate))            /**< estimator callback (can be NULL) */
);

/** gives the name of an expression handler */
SCIP_EXPORT
const char* SCIPexprhdlrGetName(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
);

/** gives the description of an expression handler (can be NULL) */
SCIP_EXPORT
const char* SCIPexprhdlrGetDescription(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
);

/** gives the precedence of an expression handler */
SCIP_EXPORT
unsigned int SCIPexprhdlrGetPrecedence(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
);

/** gives the data of an expression handler */
SCIP_EXPORT
SCIP_EXPRHDLRDATA* SCIPexprhdlrGetData(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
);

/** returns whether expression handler implements the print callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasPrint(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the backward differentiation callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasBwdiff(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the forward differentiation callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasFwdiff(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the interval evaluation callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasIntEval(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the estimator callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasEstimate(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the initial estimators callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasInitEstimates(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the simplification callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasSimplify(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the curvature callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasCurvature(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the monotonicity callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasMonotonicity(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** returns whether expression handler implements the reverse propagation callback */
SCIP_EXPORT
SCIP_Bool SCIPexprhdlrHasReverseProp(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** compares two expression handler w.r.t. their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPexprhdlrComp);

/**@name Expression Handler Statistics */
/**@{ */

/** gets number of times an expression has been created with given expression handler */
SCIP_EXPORT
unsigned int SCIPexprhdlrGetNCreated(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets number of times the interval evaluation callback was called */
SCIP_EXPORT
SCIP_Longint SCIPexprhdlrGetNIntevalCalls(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets time spend in interval evaluation callback */
SCIP_EXPORT
SCIP_Real SCIPexprhdlrGetIntevalTime(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets number of times the reverse propagation callback was called */
SCIP_EXPORT
SCIP_Longint SCIPexprhdlrGetNReversepropCalls(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets time spend in reverse propagation callback */
SCIP_EXPORT
SCIP_Real SCIPexprhdlrGetReversepropTime(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets number of times an empty interval was found in reverse propagation */
SCIP_EXPORT
SCIP_Longint SCIPexprhdlrGetNCutoffs(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets number of times a bound reduction was found in reverse propagation (and accepted by caller) */
SCIP_EXPORT
SCIP_Longint SCIPexprhdlrGetNDomainReductions(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** increments the domain reductions count of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrIncrementNDomainReductions(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   int                   nreductions         /**< number of reductions to add to counter */
   );

/** gets number of times the estimation callback was called */
SCIP_EXPORT
SCIP_Longint SCIPexprhdlrGetNEstimateCalls(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets time spend in estimation callback */
SCIP_EXPORT
SCIP_Real SCIPexprhdlrGetEstimateTime(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets number of times branching candidates reported by of this expression handler were used to assemble branching candidates
 *
 * that is, how often did we consider branching on a child of this expression
 */
SCIP_EXPORT
SCIP_Longint SCIPexprhdlrGetNBranchings(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** increments the branching candidates count of an expression handler */
SCIP_EXPORT
void SCIPexprhdlrIncrementNBranchings(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets number of times the simplify callback was called */
SCIP_EXPORT
SCIP_Longint SCIPexprhdlrGetNSimplifyCalls(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets time spend in simplify callback */
SCIP_EXPORT
SCIP_Real SCIPexprhdlrGetSimplifyTime(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** gets number of times the simplify callback found a simplification */
SCIP_EXPORT
SCIP_Longint SCIPexprhdlrGetNSimplifications(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   );

/** @} */ /* expression handler statistics */

#ifdef NDEBUG

/* If NDEBUG is defined, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlr_, freehdlr_) do { (exprhdlr)->copyhdlr = copyhdlr_; (exprhdlr)->freehdlr = freehdlr_; } while (FALSE)
#define SCIPexprhdlrSetCopyFreeData(exprhdlr, copydata_, freedata_) do { (exprhdlr)->copydata = copydata_; (exprhdlr)->freedata = freedata_; } while (FALSE)
#define SCIPexprhdlrSetPrint(exprhdlr, print_)               (exprhdlr)->print = print_
#define SCIPexprhdlrSetParse(exprhdlr, parse_)               (exprhdlr)->parse = parse_
#define SCIPexprhdlrSetCurvature(exprhdlr, curvature_)       (exprhdlr)->curvature = curvature_
#define SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicity_) (exprhdlr)->monotonicity = monotonicity_
#define SCIPexprhdlrSetIntegrality(exprhdlr, integrality_)   (exprhdlr)->integrality = integrality_
#define SCIPexprhdlrSetHash(exprhdlr, hash_)                 (exprhdlr)->hash = hash_
#define SCIPexprhdlrSetCompare(exprhdlr, compare_)           (exprhdlr)->compare = compare_
#define SCIPexprhdlrSetDiff(exprhdlr, bwdiff_, fwdiff_, bwfwdiff_) do { (exprhdlr)->bwdiff = bwdiff_; (exprhdlr)->fwdiff = fwdiff_; (exprhdlr)->bwfwdiff = bwfwdiff_; } while (FALSE)
#define SCIPexprhdlrSetIntEval(exprhdlr, inteval_)           (exprhdlr)->inteval = inteval_
#define SCIPexprhdlrSetSimplify(exprhdlr, simplify_)         (exprhdlr)->simplify = simplify_
#define SCIPexprhdlrSetReverseProp(exprhdlr, reverseprop_)   (exprhdlr)->reverseprop = reverseprop_
#define SCIPexprhdlrSetEstimate(exprhdlr, initestimates_, estimate_) do { (exprhdlr)->initestimates = initestimates_; (exprhdlr)->estimate = estimate_; } while (FALSE)
#define SCIPexprhdlrGetName(exprhdlr)              (exprhdlr)->name
#define SCIPexprhdlrGetDescription(exprhdlr)       (exprhdlr)->desc
#define SCIPexprhdlrGetPrecedence(exprhdlr)        (exprhdlr)->precedence
#define SCIPexprhdlrGetData(exprhdlr)              (exprhdlr)->data
#define SCIPexprhdlrHasPrint(exprhdlr)             ((exprhdlr)->print != NULL)
#define SCIPexprhdlrHasBwdiff(exprhdlr)            ((exprhdlr)->bwdiff != NULL)
#define SCIPexprhdlrHasFwdiff(exprhdlr)            ((exprhdlr)->fwdiff != NULL)
#define SCIPexprhdlrHasIntEval(exprhdlr)           ((exprhdlr)->inteval != NULL)
#define SCIPexprhdlrHasEstimate(exprhdlr)          ((exprhdlr)->estimate != NULL)
#define SCIPexprhdlrHasInitEstimates(exprhdlr)     ((exprhdlr)->initestimates != NULL)
#define SCIPexprhdlrHasSimplify(exprhdlr)          ((exprhdlr)->simplify != NULL)
#define SCIPexprhdlrHasCurvature(exprhdlr)         ((exprhdlr)->curvature != NULL)
#define SCIPexprhdlrHasMonotonicity(exprhdlr)      ((exprhdlr)->monotonicity != NULL)
#define SCIPexprhdlrHasReverseProp(exprhdlr)       ((exprhdlr)->reverseprop != NULL)
#define SCIPexprhdlrGetNCreated(exprhdlr)          (exprhdlr)->ncreated
#define SCIPexprhdlrGetNIntevalCalls(exprhdlr)     (exprhdlr)->nintevalcalls
#define SCIPexprhdlrGetIntevalTime(exprhdlr)       SCIPclockGetTime((exprhdlr)->intevaltime)
#define SCIPexprhdlrGetNReversepropCalls(exprhdlr) (exprhdlr)->npropcalls
#define SCIPexprhdlrGetReversepropTime(exprhdlr)   SCIPclockGetTime((exprhdlr)->proptime)
#define SCIPexprhdlrGetNCutoffs(exprhdlr)          (exprhdlr)->ncutoffs
#define SCIPexprhdlrGetNDomainReductions(exprhdlr) (exprhdlr)->ndomreds
#define SCIPexprhdlrIncrementNDomainReductions(exprhdlr, nreductions) (exprhdlr)->ndomreds += nreductions
#define SCIPexprhdlrGetNEstimateCalls(exprhdlr)    (exprhdlr)->nestimatecalls
#define SCIPexprhdlrGetEstimateTime(exprhdlr)      SCIPclockGetTime((exprhdlr)->estimatetime)
#define SCIPexprhdlrGetNBranchings(exprhdlr)       (exprhdlr)->nbranchscores
#define SCIPexprhdlrIncrementNBranchings(exprhdlr) ++(exprhdlr)->nbranchscores
#define SCIPexprhdlrGetNSimplifyCalls(exprhdlr)    (exprhdlr)->nsimplifycalls
#define SCIPexprhdlrGetSimplifyTime(exprhdlr)      SCIPclockGetTime((exprhdlr)->simplifytime)
#define SCIPexprhdlrGetNSimplifications(exprhdlr)  (exprhdlr)->nsimplified
#endif


/** @} */ /* expression handler methods */

/**@defgroup PublicExprMethods Expressions
 * @ingroup DataStructures
 * @brief an algebraic expression used for nonlinear constraints and NLPs
 *
 *@{
 */

/**@name Expressions */
/**@{ */

/** gets the number of times the expression is currently captured */
SCIP_EXPORT
int SCIPexprGetNUses(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the number of children of an expression */
SCIP_EXPORT
int SCIPexprGetNChildren(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the children of an expression (can be NULL if no children) */
SCIP_EXPORT
SCIP_EXPR** SCIPexprGetChildren(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gets the expression handler of an expression
 *
 * This identifies the type of the expression (sum, variable, ...).
 */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPexprGetHdlr(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gets the expression data of an expression */
SCIP_EXPORT
SCIP_EXPRDATA* SCIPexprGetData(
   SCIP_EXPR*            expr                /**< expression */
   );

/** sets the expression data of an expression
 *
 * The pointer to possible old data is overwritten and the
 * freedata-callback is not called before.
 * This function is intended to be used by expression handler only.
 */
SCIP_EXPORT
void SCIPexprSetData(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRDATA*        exprdata            /**< expression data to be set (can be NULL) */
   );

/** gets the data that the owner of an expression has stored in an expression */
SCIP_EXPORT
SCIP_EXPR_OWNERDATA* SCIPexprGetOwnerData(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error)
 *
 * @see SCIPevalExpr to evaluate the expression at a given solution.
 */
SCIP_EXPORT
SCIP_Real SCIPexprGetEvalValue(
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the evaluation tag from the last evaluation, or 0
 *
 * @see SCIPevalExpr
 */
SCIP_EXPORT
SCIP_Longint SCIPexprGetEvalTag(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the derivative stored in an expression (or SCIP_INVALID if there was an evaluation error)
 *
 * @see SCIPevalExprGradient
 */
SCIP_EXPORT
SCIP_Real SCIPexprGetDerivative(
   SCIP_EXPR*     expr              /**< expression */
   );

/** gives the value of directional derivative from the last evaluation of a directional derivative of expression
 * (or SCIP_INVALID if there was an error)
 *
 * @see SCIPevalExprHessianDir
 */
SCIP_EXPORT
SCIP_Real SCIPexprGetDot(
   SCIP_EXPR*     expr              /**< expression */
   );

/** gives the value of directional derivative from the last evaluation of a directional derivative of derivative
 * of root (or SCIP_INVALID if there was an error)
 *
 * @see SCIPevalExprHessianDir
 */
SCIP_EXPORT
SCIP_Real SCIPexprGetBardot(
   SCIP_EXPR*     expr              /**< expression */
   );

/** returns the difftag stored in an expression
 *
 * can be used to check whether partial derivative value is valid
 *
 * @see SCIPevalExprGradient
 */
SCIP_EXPORT
SCIP_Longint SCIPexprGetDiffTag(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the activity that is currently stored for an expression
 *
 * @see SCIPevalExprActivity
 */
SCIP_EXPORT
SCIP_INTERVAL SCIPexprGetActivity(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the tag associated with the activity of the expression
 *
 * It can depend on the owner of the expression how to interpret this tag.
 * SCIPevalExprActivity() compares with `stat->domchgcount`.
 *
 * @see SCIPevalExprActivity
 */
SCIP_EXPORT
SCIP_Longint SCIPexprGetActivityTag(
   SCIP_EXPR*            expr                /**< expression */
   );

/** set the activity with tag for an expression */
SCIP_EXPORT
void SCIPexprSetActivity(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_INTERVAL         activity,           /**< new activity */
   SCIP_Longint          activitytag         /**< tag associated with activity */
   );

/** returns the curvature of an expression
 *
 *  @note Call SCIPcomputeExprCurvature() before calling this function.
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPexprGetCurvature(
   SCIP_EXPR*            expr                /**< expression */
   );

/** sets the curvature of an expression */
SCIP_EXPORT
void SCIPexprSetCurvature(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRCURV         curvature           /**< curvature of the expression */
   );

/** returns whether an expression is integral */
SCIP_EXPORT
SCIP_Bool SCIPexprIsIntegral(
   SCIP_EXPR*            expr                /**< expression */
   );

/** sets the integrality flag of an expression */
SCIP_EXPORT
void SCIPexprSetIntegrality(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool             isintegral          /**< integrality of the expression */
   );

/** @} */

/**@name Quadratic Expressions */
/**@{ */

/** gives the coefficients and expressions that define a quadratic expression
 *
 * It can return the constant part, the number, arguments, and coefficients of the purely linear part
 * and the number of quadratic terms and bilinear terms.
 * Note that for arguments that appear in the quadratic part, a linear coefficient is
 * stored with the quadratic term.
 * Use SCIPexprGetQuadraticQuadTerm() and SCIPexprGetQuadraticBilinTerm()
 * to access the data for a quadratic or bilinear term.
 *
 * It can also return the eigenvalues and the eigenvectors of the matrix \f$Q\f$ when the quadratic is written
 * as \f$x^T Q x + b^T x + c^T y + d\f$, where \f$c^T y\f$ defines the purely linear part.
 * Note, however, that to have access to them one needs to call SCIPcomputeExprQuadraticCurvature()
 * with `storeeigeninfo=TRUE`. If the eigen information was not stored or it failed to be computed,
 * `eigenvalues` and `eigenvectors` will be set to NULL.
 *
 * This function returns pointers to internal data in linexprs and lincoefs.
 * The user must not change this data.
 *
 * @attention SCIPcheckExprQuadratic() needs to be called first to check whether expression is quadratic and initialize the data of the quadratic representation.
 */
SCIP_EXPORT
void SCIPexprGetQuadraticData(
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_Real*            constant,           /**< buffer to store constant term, or NULL */
   int*                  nlinexprs,          /**< buffer to store number of expressions that appear linearly, or NULL */
   SCIP_EXPR***          linexprs,           /**< buffer to store pointer to array of expressions that appear linearly, or NULL */
   SCIP_Real**           lincoefs,           /**< buffer to store pointer to array of coefficients of expressions that appear linearly, or NULL */
   int*                  nquadexprs,         /**< buffer to store number of expressions in quadratic terms, or NULL */
   int*                  nbilinexprs,        /**< buffer to store number of bilinear expressions terms, or NULL */
   SCIP_Real**           eigenvalues,        /**< buffer to store pointer to array of eigenvalues of Q, or NULL */
   SCIP_Real**           eigenvectors        /**< buffer to store pointer to array of eigenvectors of Q, or NULL */
   );

/** gives the data of a quadratic expression term
 *
 * For a term \f$a \cdot \text{expr}^2 + b \cdot \text{expr} + \sum_i (c_i \cdot \text{expr} \cdot \text{otherexpr}_i)\f$, returns
 * `expr`, \f$a\f$, \f$b\f$, the number of summands, and indices of bilinear terms in the quadratic expressions `bilinexprterms`.
 *
 * This function returns pointers to internal data in adjbilin.
 * The user must not change this data.
 */
SCIP_EXPORT
void SCIPexprGetQuadraticQuadTerm(
   SCIP_EXPR*            quadexpr,           /**< quadratic expression */
   int                   termidx,            /**< index of quadratic term */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to argument expression (the 'x') of this term, or NULL */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of variable, or NULL */
   SCIP_Real*            sqrcoef,            /**< buffer to store square coefficient of variable, or NULL */
   int*                  nadjbilin,          /**< buffer to store number of bilinear terms this variable is involved in, or NULL */
   int**                 adjbilin,           /**< buffer to store pointer to indices of associated bilinear terms, or NULL */
   SCIP_EXPR**           sqrexpr             /**< buffer to store pointer to square expression (the 'x^2') of this term or NULL if no square expression, or NULL */
   );

/** gives the data of a bilinear expression term
 *
 * For a term a*expr1*expr2, returns expr1, expr2, a, and
 * the position of the quadratic expression term that uses expr2 in the quadratic expressions `quadexprterms`.
 */
SCIP_EXPORT
void SCIPexprGetQuadraticBilinTerm(
   SCIP_EXPR*            expr,               /**< quadratic expression */
   int                   termidx,            /**< index of bilinear term */
   SCIP_EXPR**           expr1,              /**< buffer to store first factor, or NULL */
   SCIP_EXPR**           expr2,              /**< buffer to store second factor, or NULL */
   SCIP_Real*            coef,               /**< buffer to coefficient, or NULL */
   int*                  pos2,               /**< buffer to position of expr2 in quadexprterms array of quadratic expression, or NULL */
   SCIP_EXPR**           prodexpr            /**< buffer to store pointer to expression that is product if first and second factor, or NULL */
   );

/** returns whether all expressions that are used in a quadratic expression are variable expressions
 *
 * @return TRUE iff all `linexprs` and `quadexprterms[.].expr` are variable expressions
 */
SCIP_EXPORT
SCIP_Bool SCIPexprAreQuadraticExprsVariables(
   SCIP_EXPR*            expr                /**< quadratic expression */
   );

/** @} */

#ifdef NDEBUG
#define SCIPexprGetNUses(expr)                    (expr)->nuses
#define SCIPexprGetNChildren(expr)                (expr)->nchildren
#define SCIPexprGetChildren(expr)                 (expr)->children
#define SCIPexprGetHdlr(expr)                     (expr)->exprhdlr
#define SCIPexprGetData(expr)                     (expr)->exprdata
#define SCIPexprSetData(expr, exprdata_)          (expr)->exprdata = exprdata_
#define SCIPexprGetOwnerData(expr)                (expr)->ownerdata
#define SCIPexprGetEvalValue(expr)                (expr)->evalvalue
#define SCIPexprGetEvalTag(expr)                  (expr)->evaltag
#define SCIPexprGetDerivative(expr)               (expr)->derivative
#define SCIPexprGetDot(expr)                      (expr)->dot
#define SCIPexprGetBardot(expr)                   (expr)->bardot
#define SCIPexprGetDiffTag(expr)                  (expr)->difftag
#define SCIPexprGetActivity(expr)                 (expr)->activity
#define SCIPexprGetActivityTag(expr)              (expr)->activitytag
#define SCIPexprSetActivity(expr, activity_, activitytag_) do { (expr)->activity = activity_; (expr)->activitytag = activitytag_; } while (FALSE)
#define SCIPexprGetCurvature(expr)                (expr)->curvature
#define SCIPexprSetCurvature(expr, curvature_)    (expr)->curvature = curvature_
#define SCIPexprIsIntegral(expr)                  (expr)->isintegral
#define SCIPexprSetIntegrality(expr, isintegral_) expr->isintegral = isintegral_
#define SCIPexprAreQuadraticExprsVariables(expr)  (expr)->quaddata->allexprsarevars
#endif

/**@name Core Expression Handlers */
/**@{ */
/* these are here to have them accessible also in the expr core
 * so these cannot make use of SCIP pointer
 */

/** gets the variable of a variable expression */
SCIP_EXPORT
SCIP_VAR* SCIPgetVarExprVar(
   SCIP_EXPR*            expr                /**< var expression */
   );

/** gets the value of a constant value expression */
SCIP_EXPORT
SCIP_Real SCIPgetValueExprValue(
   SCIP_EXPR*            expr                /**< value expression */
   );

/** gets the coefficients of a summation expression */
SCIP_EXPORT
SCIP_Real* SCIPgetCoefsExprSum(
   SCIP_EXPR*            expr                /**< sum expression */
   );

/** gets the constant of a summation expression */
SCIP_EXPORT
SCIP_Real SCIPgetConstantExprSum(
   SCIP_EXPR*            expr                /**< sum expression */
   );

/** gets the constant coefficient of a product expression */
SCIP_EXPORT
SCIP_Real SCIPgetCoefExprProduct(
   SCIP_EXPR*            expr                /**< product expression */
   );

/** gets the exponent of a power or signed power expression */
SCIP_EXPORT
SCIP_Real SCIPgetExponentExprPow(
   SCIP_EXPR*            expr                /**< (signed) power expression */
   );

#ifdef NDEBUG
#define SCIPgetVarExprVar(expr) ((SCIP_VAR*)SCIPexprGetData(expr))
#endif

/**@} */


/**@name Expression Iterator
 *
 * @anchor SCIP_EXPRITER_DFS
 * More details on the DFS mode:
 * Many algorithms over expression trees need to traverse the tree in depth-first manner and a
 * natural way of implementing these algorithms is by using recursion.
 * In general, a function which traverses the tree in depth-first looks like
 * <pre>
 * fun( expr )
 *    enterexpr()
 *    continue skip or abort
 *       for( child in expr->children )
 *          visitingchild()
 *          continue skip or abort
 *          fun(child, data, proceed)
 *          visitedchild()
 *          continue skip or abort
 *    leaveexpr()
 * </pre>
 * Given that some expressions might be quite deep we provide this functionality in an iterative fashion.
 *
 * Consider an expression (x*y) + z + log(x-y).
 * The corresponding expression graph is
 * <pre>
 *           [+]
 *       /    |   \
 *    [*]     |    [log]
 *    / \     |      |
 *   /   \    |     [-]
 *   |   |    |     / \
 *  [x] [y]  [z]  [x] [y]
 * </pre>
 * (where [x] and [y] are actually the same expression).
 *
 * If a pointer to the [+] expression is given as root to this expression, it will iterate
 * the graph in a depth-first manner and stop at various stages.
 * - When entering an expression, it stops in the \ref SCIP_EXPRITER_ENTEREXPR stage.
 *   The SCIPexpriterGetParentDFS() function indicates from where the expression has been entered (NULL for the root expression).
 * - Before visiting a child of an expression, it stops in the \ref SCIP_EXPRITER_VISITINGCHILD stage.
 *   The SCIPexpriterGetChildIdxDFS() function returns which child will be visited (as an index in the current expr's children array).
 *   Use SCIPexpriterGetChildExprDFS() to obtain the corresponding expression.
 * - When returning from visiting a child of an expression, it stops in the \ref SCIP_EXPRITER_VISITEDCHILD stage.
 *   Again the SCIPexpriterGetChildExprDFS() function returns which child has been visited.
 * - When leaving an expression, it stops in the \ref SCIP_EXPRITER_LEAVEEXPR stage.
 *
 * Thus, for the above expression, the expression are visited in the following order and stages:
 * - `enterexpr([+])`
 * - `visitingchild([+])`, currentchild = 0
 * - `enterexpr([*])`
 * - `visitingchild([*])`, currentchild = 0
 * - `enterexpr([x])`
 * - `leaveexpr([x])`
 * - `visitedchild([*])`, currentchild = 0
 * - `visitingchild([*])`, currentchild = 1
 * - `enterexpr([y])`
 * - `leaveexpr([y])`
 * - `visitedchild([*])`, currentchild = 1
 * - `leaveexpr([*])`
 * - `visitedchild([+])`, currentchild = 0
 * - `visitingchild([+])`, currentchild = 1
 * - `enterexpr([z])`
 * - `leaveexpr([z])`
 * - `visitedchild([+])`, currentchild = 1
 * - `visitingchild([+])`, currentchild = 2
 * - `enterexpr([log])`
 * - `visitingchild([log])`, currentchild = 0
 * - `enterexpr([-])`
 * - `visitingchild([-])`,  currentchild = 0
 * - `enterexpr([x])`
 * - `leaveexpr([x])`
 * - `visitedchild([-])`, currentchild = 0
 * - `visitingchild([-])`, currentchild = 1
 * - `enterexpr([y])`
 * - `leaveexpr([y])`
 * - `visitedchild([-])`, currentchild = 1
 * - `leaveexpr([-])`
 * - `visitedchild([log])`, currentchild = 0
 * - `leaveexpr([log])`
 * - `visitedchild([+])` currentchild = 2
 * - `leaveexpr([+])`
 *
 * The caller can direct the iterator to skip parts of the tree:
 * - If calling SCIPexpriterSkipDFS() in SCIP_EXPRITER_ENTEREXPR stage, all children of that expression will be skipped. The SCIP_EXPRITER_LEAVEEXPR stage will still be next.
 * - If calling SCIPexpriterSkipDFS() in SCIP_EXPRITER_VISITINGCHILD stage, visiting the current child will be skipped.
 * - If calling SCIPexpriterSkipDFS() in SCIP_EXPRITER_VISITEDCHILD child, visiting the remaining children will be skipped.
 *
 * @{
 */

/** returns whether expression iterator is currently initialized */
SCIP_EXPORT
SCIP_Bool SCIPexpriterIsInit(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** initializes an expression iterator
 *
 * @note If `expr` is NULL, then iterator will be set into ended-state (SCIPexpriterIsEnd() is TRUE). Useful if following with SCIPexpriterRestartDFS().
 *
 * If type is DFS, then `stopstages` will be set to \ref SCIP_EXPRITER_ENTEREXPR.
 * Use `SCIPexpriterSetStagesDFS` to change this.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexpriterInit(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr,               /**< expression of the iterator, can be NULL */
   SCIP_EXPRITER_TYPE    type,               /**< type of expression iterator */
   SCIP_Bool             allowrevisit        /**< whether expression are allowed to be visited more than once */
   );

/** restarts an already initialized expression iterator in DFS mode
 *
 * The expression iterator will continue from the given expression, not revisiting expressions that
 * this iterator has already been visited (if initialized with `allowrevisit=FALSE`) and giving access
 * to the same iterator specified expression data that may have been set already.
 * Also the stop-stages are not reset.
 *
 * If revisiting is forbidden and given expr has already been visited, then the iterator will behave
 * as on the end of iteration (SCIPexpriterIsEnd() is TRUE).
 * If the enterexpr stage is not one of the stop stages, then the iterator will be moved forward
 * (SCIPexpriterGetNext() is called).
 *
 * @return The current expression.
 */
SCIP_EXPORT
SCIP_EXPR* SCIPexpriterRestartDFS(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr                /**< expression of the iterator */
   );

/** specifies in which stages to stop a DFS iterator
 *
 * Parameter `stopstages` should be a bitwise OR of different \ref SCIP_EXPRITER_STAGE values
 *
 * If the current stage is not one of the `stopstages`, then the iterator will be moved on.
 */
SCIP_EXPORT
void SCIPexpriterSetStagesDFS(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPRITER_STAGE   stopstages          /**< the stages in which to stop when iterating via DFS */
   );

/** gets the current expression that the expression iterator points to */
SCIP_EXPORT
SCIP_EXPR* SCIPexpriterGetCurrent(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** gets the current stage that the expression iterator is in when using DFS
 *
 * If the iterator has finished (SCIPexpriterIsEnd() is TRUE), then the stage is undefined.
 */
SCIP_EXPORT
SCIP_EXPRITER_STAGE SCIPexpriterGetStageDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** gets the index of the child that the expression iterator considers when in DFS mode and stage \ref SCIP_EXPRITER_VISITINGCHILD or \ref SCIP_EXPRITER_VISITEDCHILD */
SCIP_EXPORT
int SCIPexpriterGetChildIdxDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** gets the child expression that the expression iterator considers when in DFS mode and stage \ref SCIP_EXPRITER_VISITINGCHILD or \ref SCIP_EXPRITER_VISITEDCHILD */
SCIP_EXPORT
SCIP_EXPR* SCIPexpriterGetChildExprDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** gives the parent of the current expression of an expression iteration if in DFS mode
 *
 * @return the expression from which the current expression has been accessed
 */
SCIP_EXPORT
SCIP_EXPR* SCIPexpriterGetParentDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** gives the iterator specific user data of the current expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPORT
SCIP_EXPRITER_USERDATA SCIPexpriterGetCurrentUserData(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** gives the iterator specific user data of the current expressions current child
 *
 * @note The expression iterator mode must be in DFS mode and stage \ref SCIP_EXPRITER_VISITINGCHILD or \ref SCIP_EXPRITER_VISITEDCHILD
 */
SCIP_EXPORT
SCIP_EXPRITER_USERDATA SCIPexpriterGetChildUserDataDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** gives the iterator specific user data of a given expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPORT
SCIP_EXPRITER_USERDATA SCIPexpriterGetExprUserData(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr                /**< expression for which to get the userdata of this iterator */
   );

/** sets the iterator specific user data of the current expression for an expression iteration if in DFS mode
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPORT
void SCIPexpriterSetCurrentUserData(
   SCIP_EXPRITER*         iterator,          /**< expression iterator */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored */
   );

/** sets the iterator specific user data of a given expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPORT
void SCIPexpriterSetExprUserData(
   SCIP_EXPRITER*         iterator,          /**< expression iterator */
   SCIP_EXPR*             expr,              /**< expression where to set iterator data */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored in current child */
   );

/** sets the iterator specific user data of the current expressions current child
 *
 * @note The expression iterator mode must be in DFS mode and stage \ref SCIP_EXPRITER_VISITINGCHILD or \ref SCIP_EXPRITER_VISITEDCHILD
 */
SCIP_EXPORT
void SCIPexpriterSetChildUserData(
   SCIP_EXPRITER*         iterator,          /**< expression iterator */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored in current child */
   );

/** moves the iterator to the next expression according to the mode of the expression iterator
 *
 * @return the next expression, if any, and NULL otherwise
 */
SCIP_EXPORT
SCIP_EXPR* SCIPexpriterGetNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** moves a DFS iterator to one of the next expressions
 *
 * - If in \ref SCIP_EXPRITER_ENTEREXPR stage, then all children of that expression will be skipped.
 *   If \ref SCIP_EXPRITER_LEAVEEXPR is one of the `stopstages`, then it will be the next stage. Otherwise, the iterator will move further on (go to the parent, etc).
 * - If in \ref SCIP_EXPRITER_VISITINGCHILD stage, then the child that was going to be visited next will be skipped and the iterator will be moved on to the next child (if any).
 * - If in \ref SCIP_EXPRITER_VISITEDCHILD stage, then all remaining children will be skipped and we move on to the \ref SCIP_EXPRITER_LEAVEEXPR stage (if a stop stage, otherwise further on).
 * - It is not allowed to call this function when in \ref SCIP_EXPRITER_LEAVEEXPR stage.
 *
 * @return the next expression, if any, and NULL otherwise
 */
SCIP_EXPORT
SCIP_EXPR* SCIPexpriterSkipDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

/** returns whether the iterator visited all expressions already */
SCIP_EXPORT
SCIP_Bool SCIPexpriterIsEnd(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   );

#ifdef NDEBUG
#define SCIPexpriterIsInit(iterator)                           (iterator)->initialized
#define SCIPexpriterGetCurrent(iterator)                       (iterator)->curr
#define SCIPexpriterGetStageDFS(iterator)                      (iterator)->dfsstage
#define SCIPexpriterGetChildIdxDFS(iterator)                   (iterator)->curr->iterdata[(iterator)->iterindex].currentchild
#define SCIPexpriterGetChildExprDFS(iterator)                  (iterator)->curr->children[(iterator)->curr->iterdata[(iterator)->iterindex].currentchild]
#define SCIPexpriterGetParentDFS(iterator)                     (iterator)->curr->iterdata[(iterator)->iterindex].parent
#define SCIPexpriterGetCurrentUserData(iterator)               (iterator)->curr->iterdata[(iterator)->iterindex].userdata
#define SCIPexpriterGetChildUserDataDFS(iterator)              (iterator)->curr->children[(iterator)->curr->iterdata[(iterator)->iterindex].currentchild]->iterdata[(iterator)->iterindex].userdata
#define SCIPexpriterGetExprUserData(iterator, expr)            (expr)->iterdata[(iterator)->iterindex].userdata
#define SCIPexpriterSetCurrentUserData(iterator, userdata_)    (iterator)->curr->iterdata[(iterator)->iterindex].userdata = userdata_
#define SCIPexpriterSetExprUserData(iterator, expr, userdata_) (expr)->iterdata[(iterator)->iterindex].userdata = userdata_
#define SCIPexpriterSetChildUserData(iterator, userdata_)      (iterator)->curr->children[(iterator)->curr->iterdata[(iterator)->iterindex].currentchild]->iterdata[(iterator)->iterindex].userdata = userdata_
#define SCIPexpriterIsEnd(iterator)                            ((iterator)->curr == NULL)
#endif

/** @} */

/**@name Function Curvature */
/**@{ */

/** gives curvature for a sum of two functions with given curvature */
SCIP_EXPORT
SCIP_EXPRCURV SCIPexprcurvAdd(
   SCIP_EXPRCURV         curv1,              /**< curvature of first summand */
   SCIP_EXPRCURV         curv2               /**< curvature of second summand */
   );

/** gives the curvature for the negation of a function with given curvature */
SCIP_EXPORT
SCIP_EXPRCURV SCIPexprcurvNegate(
   SCIP_EXPRCURV         curvature           /**< curvature of function */
   );

/** gives curvature for a functions with given curvature multiplied by a constant factor */
SCIP_EXPORT
SCIP_EXPRCURV SCIPexprcurvMultiply(
   SCIP_Real             factor,             /**< constant factor */
   SCIP_EXPRCURV         curvature           /**< curvature of other factor */
   );

/** gives curvature for base^exponent for given bounds and curvature of base-function and constant exponent */
SCIP_EXPORT
SCIP_EXPRCURV SCIPexprcurvPower(
   SCIP_INTERVAL         basebounds,         /**< bounds on base function */
   SCIP_EXPRCURV         basecurv,           /**< curvature of base function */
   SCIP_Real             exponent            /**< exponent */
   );

/** gives required curvature for base so that base^exponent has given curvature under given bounds on base and constant exponent
 *
 * returns curvature unknown if expected curvature cannot be obtained
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPexprcurvPowerInv(
   SCIP_INTERVAL         basebounds,         /**< bounds on base function */
   SCIP_Real             exponent,           /**< exponent, must not be 0 */
   SCIP_EXPRCURV         powercurv           /**< expected curvature for power */
   );

/** gives curvature for a monomial with given curvatures and bounds for each factor
 *
 *  See Maranas and Floudas, Finding All Solutions of Nonlinearly Constrained Systems of Equations, JOGO 7, 1995
 *  for the categorization in the case that all factors are linear.
 *
 *  Exponents can also be negative or rational.
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPexprcurvMonomial(
   int                   nfactors,           /**< number of factors in monomial */
   SCIP_Real*            exponents,          /**< exponents in monomial, or NULL if all 1.0 */
   int*                  factoridxs,         /**< indices of factors, or NULL if identity mapping */
   SCIP_EXPRCURV*        factorcurv,         /**< curvature of each factor */
   SCIP_INTERVAL*        factorbounds        /**< bounds of each factor */
   );

/** for a monomial with given bounds for each factor, gives condition on the curvature of each factor, so that monomial has a requested curvature, if possible
 *
 * @return whether `monomialcurv` can be achieved
 */
SCIP_EXPORT
SCIP_Bool SCIPexprcurvMonomialInv(
   SCIP_EXPRCURV         monomialcurv,       /**< desired curvature */
   int                   nfactors,           /**< number of factors in monomial */
   SCIP_Real*            exponents,          /**< exponents in monomial, or NULL if all 1.0 */
   SCIP_INTERVAL*        factorbounds,       /**< bounds of each factor */
   SCIP_EXPRCURV*        factorcurv          /**< buffer to store required curvature of each factor */
   );

/** gives name as string for a curvature */
SCIP_EXPORT
const char* SCIPexprcurvGetName(
   SCIP_EXPRCURV         curv                /**< curvature */
   );

#ifdef NDEBUG
#define SCIPexprcurvAdd(curv1, curv2)  ((SCIP_EXPRCURV) ((curv1) & (curv2)))
#define SCIPexprcurvNegate(curvature)  (((curvature) == SCIP_EXPRCURV_CONCAVE) ? SCIP_EXPRCURV_CONVEX : ((curvature) == SCIP_EXPRCURV_CONVEX) ? SCIP_EXPRCURV_CONCAVE : (curvature))
#define SCIPexprcurvMultiply(factor, curvature) (((factor) == 0.0) ? SCIP_EXPRCURV_LINEAR : (factor) > 0.0 ? (curvature) : SCIPexprcurvNegate(curvature))
#endif

/**@} */

/**@} */

#ifdef __cplusplus
}
#endif

#endif /* SCIP_PUB_EXPR_H_ */
