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

/**@file   type_nlhdlr.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions related to nonlinear handlers of nonlinear constraints
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 *
 *  This file defines the interface for nonlinear handlers.
 *
 *  - \ref NLHDLRS "List of available nonlinear handlers"
 */

/** @defgroup DEFPLUGINS_NLHDLR Default nonlinear handlers
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default nonlinear handlers of SCIP
 */

#ifndef SCIP_TYPE_NLHDLR_H_
#define SCIP_TYPE_NLHDLR_H_

#include "scip/type_expr.h"
#include "scip/type_cons.h"
#include "scip/type_misc.h"

#define SCIP_NLHDLR_METHOD_NONE       0x0u   /**< no enforcement */
#define SCIP_NLHDLR_METHOD_SEPABELOW  0x1u   /**< separation for expr <= auxvar, thus might estimate expr from below */
#define SCIP_NLHDLR_METHOD_SEPAABOVE  0x2u   /**< separation for expr >= auxvar, thus might estimate expr from above */
#define SCIP_NLHDLR_METHOD_SEPABOTH   (SCIP_NLHDLR_METHOD_SEPABELOW | SCIP_NLHDLR_METHOD_SEPAABOVE)  /**< separation for expr == auxvar */
#define SCIP_NLHDLR_METHOD_ACTIVITY   0x4u   /**< activity computation (interval evaluation) and propagation (reverse propagation) */
#define SCIP_NLHDLR_METHOD_ALL        (SCIP_NLHDLR_METHOD_SEPABOTH | SCIP_NLHDLR_METHOD_ACTIVITY) /**< all enforcement methods */

typedef unsigned int SCIP_NLHDLR_METHOD; /**< nlhdlr methods bitflags */

/** nonlinear handler copy callback
 *
 * The method includes the nonlinear handler into a nonlinear constraint handler.
 *
 * This method is usually called when doing a copy of a nonlinear constraint handler.
 *
 * \param[in] targetscip      target SCIP main data structure
 * \param[in] targetconshdlr  target nonlinear constraint handler
 * \param[out] sourceconshdlr nonlinear constraint handler in source SCIP
 * \param[out] sourcenlhdlr   nonlinear handler in source SCIP
 */
#define SCIP_DECL_NLHDLRCOPYHDLR(x) SCIP_RETCODE x (\
   SCIP*          targetscip,     \
   SCIP_CONSHDLR* targetconshdlr, \
   SCIP_CONSHDLR* sourceconshdlr, \
   SCIP_NLHDLR*   sourcenlhdlr)

/** callback to free data of handler
 *
 * \param[in] scip       SCIP data structure
 * \param[in] nlhdlr     nonlinear handler
 * \param[in] nlhdlrdata nonlinear handler data to be freed
 */
#define SCIP_DECL_NLHDLRFREEHDLRDATA(x) SCIP_RETCODE x (\
   SCIP*             scip,   \
   SCIP_NLHDLR*      nlhdlr, \
   SCIP_NLHDLRDATA** nlhdlrdata)

/** callback to free expression specific data
 *
 * \param[in] scip           SCIP data structure
 * \param[in] nlhdlr         nonlinear handler
 * \param[in] expr           expression
 * \param[in] nlhdlrexprdata nonlinear handler expression data to be freed
 */
#define SCIP_DECL_NLHDLRFREEEXPRDATA(x) SCIP_RETCODE x (\
   SCIP*                 scip,   \
   SCIP_NLHDLR*          nlhdlr, \
   SCIP_EXPR*            expr,   \
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata)

/** callback to be called in initialization (called after problem was transformed)
 *
 * \param[in] scip   SCIP data structure
 * \param[in] nlhdlr nonlinear handler
 */
#define SCIP_DECL_NLHDLRINIT(x) SCIP_RETCODE x (\
   SCIP*        scip, \
   SCIP_NLHDLR* nlhdlr)

/** callback to be called in deinitialization (called before transformed problem is freed)
 *
 * \param[in] scip   SCIP data structure
 * \param[in] nlhdlr nonlinear handler
 */
#define SCIP_DECL_NLHDLREXIT(x) SCIP_RETCODE x (\
   SCIP*        scip, \
   SCIP_NLHDLR* nlhdlr)

/** callback to detect structure in expression
 *
 * The nonlinear handler shall analyze the current expression and decide whether it wants to contribute
 * in enforcing the relation between this expression (`expr`) and its descendants (e.g., children) via
 * linear under- or overestimation, cut generation, and/or activity computation and propagation.
 * For linear under- or overestimation and cut generation, an auxiliary variable (`auxvar`) can be assumed to
 * be associated with `expr` and auxiliary variables may be requested in descendant expressions.
 *
 * We distinguish the following enforcement methods:
 * - \ref SCIP_NLHDLR_METHOD_SEPABELOW : linear underestimation of `expr` or cut generation for the relation `expr` &le; `auxvar` (denoted as "below")
 * - \ref SCIP_NLHDLR_METHOD_SEPAABOVE : linear overestimation of `expr` or cut generation for the relation `expr` &ge; `auxvar` (denoted as "above")
 * - \ref SCIP_NLHDLR_METHOD_ACTIVITY  : domain propagation (i.e., constant under/overestimation) for the relation `expr` = `auxvar`.
 *
 * On input, parameter `enforcing` indicates for any of these methods, whether
 * - it is not necessary to have such a method, e.g., because no `auxvar` will exist for `expr`, or no one uses or sets activities of this expression,
 *   or because analysis of the expression has shown that a relation like `expr` &ge; `auxvar` is not necessary to be satisfied,
 * - or there already exists a nonlinear handler that will provide this method in an "enforcement" sense, that is,
 *   it believes that no one else could provide this method in a stronger sense. (This is mainly used by nlhdlr_default to check whether
 *   it should still reach out to the exprhdlr or whether it would be dominated by some nonlinear handler.)
 *
 * The DETECT callback shall augment the `enforcing` bitmask by setting the enforcement methods it wants to provide in an "enforcement" sense.
 *
 * Additionally, the `participating` bitmask shall be set if the nonlinear handler wants to be called on this expression at all.
 * Here, it shall set all methods that it wants to provide, which are those set in `enforcing`, but additionally those where it wants
 * to participate but leave enforcement to another nonlinear handler.
 * This can be useful for nonlinear handlers that do not implement a complete enforcement, e.g., a handler that only contributes
 * cutting planes in some situations only.
 *
 * A nonlinear handler will be called only for those callbacks that it mentioned in `participating`, which is
 * - \ref SCIP_DECL_NLHDLRENFO "ENFO" and/or \ref SCIP_DECL_NLHDLRESTIMATE "ESTIMATE" will be called with `overestimate==FALSE` if \ref SCIP_NLHDLR_METHOD_SEPABELOW has been set
 * - \ref SCIP_DECL_NLHDLRENFO "ENFO" and/or \ref SCIP_DECL_NLHDLRESTIMATE "ESTIMATE" will be called with `overestimate==TRUE`  if \ref SCIP_NLHDLR_METHOD_SEPAABOVE has been set
 * - \ref SCIP_DECL_NLHDLRINTEVAL "INTEVAL" and/or \ref SCIP_DECL_NLHDLRREVERSEPROP "REVERSEPROP" will be called if \ref SCIP_NLHDLR_METHOD_ACTIVITY has been set
 *
 * If \ref SCIP_NLHDLR_METHOD_SEPABELOW or \ref SCIP_NLHDLR_METHOD_SEPAABOVE has been set, then at least one of the
 * callbacks \ref SCIP_DECL_NLHDLRENFO "ENFO" and \ref SCIP_DECL_NLHDLRESTIMATE "ESTIMATE" needs to be implemented.
 * Also \ref SCIP_DECL_NLHDLREVALAUX "EVALAUX" will be called in this case.
 * If \ref SCIP_NLHDLR_METHOD_ACTIVITY has been set, then at least one of \ref SCIP_DECL_NLHDLRINTEVAL "INTEVAL" and
 * \ref SCIP_DECL_NLHDLRREVERSEPROP "REVERSEPROP" needs to be implemented.
 * If the nonlinear handler chooses not to participate, then it must not set `nlhdlrexprdata` and can leave `participating` at its
 * initial value (\ref SCIP_NLHDLR_METHOD_NONE).
 *
 * Additionally, a nonlinear handler that decides to participate in any of the enforcement methods must call
 * @ref SCIPregisterExprUsageNonlinear() for every subexpression that it will use and indicate whether
 * - it will use an auxiliary variable in \ref SCIP_DECL_NLHDLRENFO "ENFO" or \ref SCIP_DECL_NLHDLRESTIMATE "ESTIMATE",
 * - it will use activity for some subexpressions when computing estimators or cuts, and
 * - it will use activity for some subexpressions when in \ref SCIP_DECL_NLHDLRINTEVAL "INTEVAL" or \ref SCIP_DECL_NLHDLRREVERSEPROP "REVERSEPROP".
 *
 * @note Auxiliary variables do not exist in subexpressions during DETECT and are not created by a call to @ref SCIPregisterExprUsageNonlinear().
 *   They will be available when the \ref SCIP_DECL_NLHDLRINITSEPA "INITSEPA" callback is called.
 *
 * \param[in] scip              SCIP data structure
 * \param[in] conshdlr          nonlinear constraint handler
 * \param[in] nlhdlr            nonlinear handler
 * \param[in] expr              expression to analyze
 * \param[in] cons              the constraint that expression defines, or NULL when the expr does not define any constraint, that is, when it is not the root of an expression of a constraint
 * \param[in,out] enforcing     enforcement methods that are provided by some nonlinear handler (to be updated by detect callback)
 * \param[out] participating    enforcement methods that this nonlinear handler should be called for (to be set by detect callback), initialized to SCIP_NLHDLR_METHOD_NONE
 * \param[out] nlhdlrexprdata   nlhdlr's expr data to be stored in expr, can only be set to non-NULL if success is set to TRUE
 */
#define SCIP_DECL_NLHDLRDETECT(x) SCIP_RETCODE x (\
   SCIP*                 scip,          \
   SCIP_CONSHDLR*        conshdlr,      \
   SCIP_NLHDLR*          nlhdlr,        \
   SCIP_EXPR*            expr,          \
   SCIP_CONS*            cons,          \
   SCIP_NLHDLR_METHOD*   enforcing,     \
   SCIP_NLHDLR_METHOD*   participating, \
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata)

/** auxiliary evaluation callback of nonlinear handler
 *
 * Evaluates the expression w.r.t. the auxiliary variables that were introduced by the nonlinear handler (if any).
 * The method is used to determine the violation of the relation that the nonlinear handler attempts to enforce.
 * During enforcement, this violation value is used to decide whether estimation/separation callbacks should be called.
 *
 * It can be assumed that the expression itself has been evaluated in the given sol.
 *
 * \param[in]  scip           SCIP data structure
 * \param[in]  nlhdlr         nonlinear handler
 * \param[in]  expr           expression to evaluate
 * \param[in]  nlhdlrexprdata expression specific data of the nonlinear handler
 * \param[out] auxvalue       buffer to store value of expression w.r.t. auxiliary variables
 * \param[in]  sol            point to evaluate
 */
#define SCIP_DECL_NLHDLREVALAUX(x) SCIP_RETCODE x (\
   SCIP*                scip,           \
   SCIP_NLHDLR*         nlhdlr,         \
   SCIP_EXPR*           expr,           \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_Real*           auxvalue,       \
   SCIP_SOL*            sol)

/** nonlinear handler interval evaluation (activity computation) callback
 *
 * The method computes an interval that contains the image (range) of the expression.
 *
 * \param[in] scip           SCIP main data structure
 * \param[in] nlhdlr         nonlinear handler
 * \param[in] expr           expression
 * \param[in] nlhdlrexprdata expression specific data of the nonlinear handler
 * \param[in,out] interval   buffer where to store interval (on input: current interval for expr, on output: computed interval for expr)
 * \param[in] intevalvar     callback to be called when interval evaluating a variable
 * \param[in] intevalvardata data to be passed to intevalvar callback
 */
#define SCIP_DECL_NLHDLRINTEVAL(x) SCIP_RETCODE x (\
   SCIP*                scip,                \
   SCIP_NLHDLR*         nlhdlr,              \
   SCIP_EXPR*           expr,                \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata,      \
   SCIP_INTERVAL*       interval,            \
   SCIP_DECL_EXPR_INTEVALVAR((*intevalvar)), \
   void*                intevalvardata)

/** nonlinear handler callback for reverse propagation
 *
 * The method propagates the given bounds over the arguments of an expression.
 * The arguments of an expression are other expressions and the tighter intervals should be passed
 * to the corresponding argument (expression) via SCIPtightenExprIntervalNonlinear().
 *
 * \param[in] scip           SCIP main data structure
 * \param[in] conshdlr       nonlinear constraint handler
 * \param[in] nlhdlr         nonlinear handler
 * \param[in] expr           expression
 * \param[in] nlhdlrexprdata expression specific data of the nonlinear handler
 * \param[in] bounds         the bounds on the expression that should be propagated
 * \param[out] infeasible    buffer to store whether an expression's bounds were propagated to an empty interval
 * \param[out] nreductions   buffer to store the number of interval reductions of all children
 */
#define SCIP_DECL_NLHDLRREVERSEPROP(x) SCIP_RETCODE x (\
   SCIP*                scip,           \
   SCIP_CONSHDLR*       conshdlr,       \
   SCIP_NLHDLR*         nlhdlr,         \
   SCIP_EXPR*           expr,           \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_INTERVAL        bounds,         \
   SCIP_Bool*           infeasible,     \
   int*                 nreductions)

/** separation initialization method of a nonlinear handler (called during CONSINITLP)
 *
 * The method shall initialize the separation data of the nonlinear handler, if any, and add initial cuts to the LP relaxation.
 *
 * \param[in] scip           SCIP main data structure
 * \param[in] conshdlr       nonlinear constraint handler
 * \param[in] cons           nonlinear constraint
 * \param[in] nlhdlr         nonlinear handler
 * \param[in] nlhdlrexprdata exprdata of nonlinear handler
 * \param[in] expr           expression
 * \param[in] overestimate   whether the expression needs to be overestimated
 * \param[in] underestimate  whether the expression needs to be underestimated
 * \param[out] infeasible    buffer to store whether infeasibility was detected while building the LP
 */
#define SCIP_DECL_NLHDLRINITSEPA(x) SCIP_RETCODE x (\
   SCIP*                scip,           \
   SCIP_CONSHDLR*       conshdlr,       \
   SCIP_CONS*           cons,           \
   SCIP_NLHDLR*         nlhdlr,         \
   SCIP_EXPR*           expr,           \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata, \
   SCIP_Bool            overestimate,   \
   SCIP_Bool            underestimate,  \
   SCIP_Bool*           infeasible)

/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL)
 *
 * The method shall deinitialize the separation data of the nonlinear handler, if any.
 *
 * \param[in] scip           SCIP main data structure
 * \param[in] nlhdlr         nonlinear handler
 * \param[in] nlhdlrexprdata exprdata of nonlinear handler
 * \param[in] expr           expression
 */
#define SCIP_DECL_NLHDLREXITSEPA(x) SCIP_RETCODE x (\
   SCIP*                scip,   \
   SCIP_NLHDLR*         nlhdlr, \
   SCIP_EXPR*           expr,   \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata)

/** nonlinear handler separation and enforcement callback
 *
 * The method tries to separate the given solution from the set defined by either
 * <pre>
 *   expr - auxvar <= 0 (if !overestimate)
 * </pre>
 * or
 * <pre>
 *   expr - auxvar >= 0 (if  overestimate),
 * </pre>
 * where `auxvar = SCIPgetExprAuxVarNonlinear(expr)`.
 *
 * It can do so by
 * - separation, i.e., finding an affine hyperplane (a cut) that separates the given point,
 * - bound tightening, i.e., changing bounds on a variable so that the given point is outside the updated domain,
 * - adding branching scores to potentially split the current problem into 2 subproblems
 *
 * If parameter `inenforcement` is FALSE, then only the first option (separation) is allowed.
 *
 * If the nonlinear handler always separates by computing a linear under- or overestimator of expr,
 * then it could be advantageous to implement the \ref SCIP_DECL_NLHDLRESTIMATE "ESTIMATE" callback instead.
 *
 * Note, that the nonlinear handler may also choose to separate for a relaxation of the mentioned sets,
 * e.g., `expr` &le; upperbound(`auxvar`)  or  `expr` &ge; lowerbound(`auxvar`).
 * This is especially useful in situations where `expr` is the root expression of a constraint
 * and it is sufficient to satisfy `lhs` &le; `expr` &le; `rhs`.
 * cons_nonlinear ensures that `lhs` &le; lowerbound(`auxvar`) and upperbound(`auxvar`) &le; `rhs`.
 *
 * cons_nonlinear may call this callback first with `allowweakcuts` = FALSE and repeat later with
 * `allowweakcuts` = TRUE, if it didn't succeed to enforce a solution without using weak cuts.
 * If in enforcement and the nonlinear handler cannot enforce by separation or bound tightening, it should register
 * branching scores for those expressions where branching may help to compute tighter cuts in children.
 *
 * The nonlinear handler must set `result` to \ref SCIP_SEPARATED if it added a cut,
 * to \ref SCIP_REDUCEDDOM if it added a bound change, and
 * to \ref SCIP_BRANCHED if it added branching scores.
 * Otherwise, it may set result to \ref SCIP_DIDNOTRUN or \ref SCIP_DIDNOTFIND.
 *
 * Parameter `cons` gives the constraint that is currently enforced.
 * Note that `expr` does not need to be the root of this constraint, i.e., `SCIPgetExprNonlinear(cons)==expr` may not hold.
 * If an expression appears in several constraints, it is not well defined which constraint is given in `cons`.
 * The main purpose of `cons` is to provide a constraint source for LP rows that are added in this callback.
 *
 * \param[in] scip           SCIP main data structure
 * \param[in] conshdlr       nonlinear constraint handler
 * \param[in] cons           nonlinear constraint that is currently enforced
 * \param[in] nlhdlr         nonlinear handler
 * \param[in] expr           expression
 * \param[in] nlhdlrexprdata expression specific data of the nonlinear handler
 * \param[in] sol            solution to be separated (NULL for the LP solution)
 * \param[in] auxvalue       current value of expression w.r.t. auxiliary variables as obtained from \ref SCIP_DECL_NLHDLREVALAUX "EVALAUX"
 * \param[in] overestimate   whether the expression needs to be over- or underestimated
 * \param[in] allowweakcuts  whether we should only look for "strong" cuts, or anything that separates is fine
 * \param[in] separated      whether another nonlinear handler already added a cut for this expression
 * \param[in] inenforcement  whether we are in enforcement, or only in separation
 * \param[out] result        pointer to store the result
 */
#define SCIP_DECL_NLHDLRENFO(x) SCIP_RETCODE x (\
   SCIP*                scip,            \
   SCIP_CONSHDLR*       conshdlr,        \
   SCIP_CONS*           cons,            \
   SCIP_NLHDLR*         nlhdlr,          \
   SCIP_EXPR*           expr,            \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata,  \
   SCIP_SOL*            sol,             \
   SCIP_Real            auxvalue,        \
   SCIP_Bool            overestimate,    \
   SCIP_Bool            allowweakcuts,   \
   SCIP_Bool            separated,       \
   SCIP_Bool            addbranchscores, \
   SCIP_RESULT*         result)

/** nonlinear handler under/overestimation callback
 *
 * The method tries to compute linear under- or overestimators of `expr` that are as tight as possible at a given point.
 * If the value of the estimator in the solution is smaller (larger) than `targetvalue`
 * when underestimating (overestimating), then no estimator needs to be computed.
 * Note, that `targetvalue` can be infinite if any estimator will be accepted.
 * If successful, it shall store the estimators in the given `rowpreps` data structure and set the
 * `rowprep->local` flag accordingly (SCIProwprepSetLocal()).
 * The sidetype of a rowprep must be set to \ref SCIP_SIDETYPE_LEFT if overestimating and
 * \ref SCIP_SIDETYPE_RIGHT if underestimating.
 *
 * If the callback is required to indicate for which expression a reduction in the local bounds (usually by branching)
 * would improve the estimator, it shall do so via calls to SCIPaddExprsViolScoreNonlinear().
 *
 * \param[in] scip               SCIP main data structure
 * \param[in] conshdlr           constraint handler
 * \param[in] nlhdlr             nonlinear handler
 * \param[in] expr               expression
 * \param[in] nlhdlrexprdata     expression data of nonlinear handler
 * \param[in] sol                solution at which to estimate (NULL for the LP solution)
 * \param[in] auxvalue           current value of expression w.r.t. auxiliary variables as obtained from \ref SCIP_DECL_NLHDLREVALAUX "EVALAUX"
 * \param[in] overestimate       whether the expression needs to be over- or underestimated
 * \param[in] targetvalue        a value the estimator shall exceed, can be +/-infinity
 * \param[in] addbranchscores    indicates whether to register branching scores
 * \param[out] rowpreps          an array where to store the estimators
 * \param[out] success           buffer to indicate whether an estimator could be computed
 * \param[out] addedbranchscores buffer to store whether the branching score callback was successful
 */
#define SCIP_DECL_NLHDLRESTIMATE(x) SCIP_RETCODE x (\
   SCIP*                scip,            \
   SCIP_CONSHDLR*       conshdlr,        \
   SCIP_NLHDLR*         nlhdlr,          \
   SCIP_EXPR*           expr,            \
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata,  \
   SCIP_SOL*            sol,             \
   SCIP_Real            auxvalue,        \
   SCIP_Bool            overestimate,    \
   SCIP_Real            targetvalue,     \
   SCIP_Bool            addbranchscores, \
   SCIP_PTRARRAY*       rowpreps,        \
   SCIP_Bool*           success,         \
   SCIP_Bool*           addedbranchscores)

typedef struct SCIP_Nlhdlr         SCIP_NLHDLR;          /**< nonlinear handler */
typedef struct SCIP_NlhdlrData     SCIP_NLHDLRDATA;      /**< nonlinear handler data */
typedef struct SCIP_NlhdlrExprData SCIP_NLHDLREXPRDATA;  /**< nonlinear handler data for a specific expression */

#endif /* SCIP_TYPE_NLHDLR_H_ */
