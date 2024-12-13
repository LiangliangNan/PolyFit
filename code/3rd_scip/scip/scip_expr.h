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

/**@file   scip_expr.h
 * @ingroup PUBLICCOREAPI
 * @brief  public functions to work with algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_SCIP_EXPR_H_
#define SCIP_SCIP_EXPR_H_

#include "scip/type_scip.h"
#include "scip/type_expr.h"
#include "scip/type_misc.h"

#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_mem.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/expr.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicExprHandlerMethods
 * @{
 */

/** creates the handler for an expression handler and includes it into SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlr(
   SCIP*                 scip,         /**< SCIP data structure */
   SCIP_EXPRHDLR**       exprhdlr,     /**< buffer where to store created expression handler */
   const char*           name,         /**< name of expression handler (must not be NULL) */
   const char*           desc,         /**< description of expression handler (can be NULL) */
   unsigned int          precedence,   /**< precedence of expression operation (used for printing) */
   SCIP_DECL_EXPREVAL((*eval)),        /**< point evaluation callback (must not be NULL) */
   SCIP_EXPRHDLRDATA*    data          /**< data of expression handler (can be NULL) */
   );

/** gives expression handlers */
SCIP_EXPORT
SCIP_EXPRHDLR** SCIPgetExprhdlrs(
   SCIP*                      scip           /**< SCIP data structure */
);

/** gives number of expression handlers */
SCIP_EXPORT
int SCIPgetNExprhdlrs(
   SCIP*                      scip           /**< SCIP data structure */
);

/** returns an expression handler of a given name (or NULL if not found) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPfindExprhdlr(
   SCIP*                      scip,          /**< SCIP data structure */
   const char*                name           /**< name of expression handler */
   );

/** returns expression handler for variable expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprhdlrVar(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for constant value expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprhdlrValue(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for sum expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprhdlrSum(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for product expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprhdlrProduct(
   SCIP*                      scip           /**< SCIP data structure */
   );

/** returns expression handler for power expressions (or NULL if not included) */
SCIP_EXPORT
SCIP_EXPRHDLR* SCIPgetExprhdlrPower(
   SCIP*                      scip           /**< SCIP data structure */
   );

#ifdef NDEBUG
/* If NDEBUG is defined, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */
#define SCIPgetExprhdlrs(scip)       (scip)->set->exprhdlrs
#define SCIPgetNExprhdlrs(scip)      (scip)->set->nexprhdlrs
#define SCIPfindExprhdlr(scip, name) SCIPsetFindExprhdlr((scip)->set, name)
#define SCIPgetExprhdlrVar(scip)     (scip)->set->exprhdlrvar
#define SCIPgetExprhdlrValue(scip)   (scip)->set->exprhdlrval
#define SCIPgetExprhdlrSum(scip)     (scip)->set->exprhdlrsum
#define SCIPgetExprhdlrProduct(scip) (scip)->set->exprhdlrproduct
#define SCIPgetExprhdlrPower(scip)   (scip)->set->exprhdlrpow
#endif

/** @} */

/**@addtogroup PublicExprMethods
 * @{
 */

/**@name Expressions */
/**@{ */

/** creates and captures an expression with given expression data and children */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data (expression assumes ownership) */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children (can be NULL if nchildren is 0) */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** creates and captures an expression with given expression data and up to two children */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExpr2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data */
   SCIP_EXPR*            child1,             /**< first child (can be NULL) */
   SCIP_EXPR*            child2,             /**< second child (can be NULL) */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** creates and captures an expression representing a quadratic function */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** creates and captures an expression representing a monomial
 *
 * @note In deviation from the actual definition of monomials, we also allow for negative and rational exponents.
 * So this function actually creates an expression for a signomial that has exactly one term.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nfactors,           /**< number of factors in monomial */
   SCIP_VAR**            vars,               /**< variables in the monomial */
   SCIP_Real*            exponents,          /**< exponent in each factor, or NULL if all 1.0 */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** appends child to the children list of expr
 *
 * @attention Only use if you really know what you are doing. The expression handler of the expression needs to be able to handle an increase in the number of children.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPappendExprChild(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPR*            child               /**< expression to be appended */
   );

/** overwrites/replaces a child of an expressions
 *
 * The old child is released and the newchild is captured, unless they are the same (=same pointer).
 */
SCIP_EXPORT
SCIP_RETCODE SCIPreplaceExprChild(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPR*              expr,             /**< expression which is going to replace a child */
   int                     childidx,         /**< index of child being replaced */
   SCIP_EXPR*              newchild          /**< the new child */
   );

/** remove all children of expr
 *
 * @attention Only use if you really know what you are doing. The expression handler of the expression needs to be able to handle the removal of all children.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPremoveExprChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** duplicates the given expression and its children */
SCIP_EXPORT
SCIP_RETCODE SCIPduplicateExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_DECL_EXPR_MAPEXPR((*mapexpr)),       /**< expression mapping function, or NULL for creating new expressions */
   void*                 mapexprdata,        /**< data of expression mapping function */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** duplicates the given expression, but reuses its children */
SCIP_EXPORT
SCIP_RETCODE SCIPduplicateExprShallow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store (shallow) duplicate of expr */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** copies an expression including children to use in a (possibly different) SCIP instance */
SCIP_EXPORT
SCIP_RETCODE SCIPcopyExpr(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata,    /**< data to pass to ownercreate */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid               /**< pointer to store whether all checked or enforced constraints were validly copied */
   );

/** creates an expression from a string
 *
 * We specify the grammar that defines the syntax of an expression.
 * Loosely speaking, a `Base` will be any "block", a `Factor` is a `Base` to a power,
 * a `Term` is a product of `Factors` and an `Expression` is a sum of `Terms`.
 *
 * The actual definition:
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * Term       -> Factor { ("*" | "/" ) Factor }
 * Factor     -> Base [ "^" "number" | "^(" "number" ")" ]
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 * where `[a|b]` means `a` or `b` or none, `(a|b)` means `a` or `b`, `{a}` means 0 or more `a`.
 *
 * Note that `Op` and `OpExpression` are undefined.
 * `Op` corresponds to the name of an expression handler and `OpExpression` to whatever string the expression handler accepts (through its parse method).
 */
SCIP_EXPORT
SCIP_RETCODE SCIPparseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer to store the expr parsed */
   const char*           exprstr,            /**< string with the expr to parse */
   const char**          finalpos,           /**< buffer to store the position of exprstr where we finished reading, or NULL if not of interest */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** captures an expression (increments usage count) */
SCIP_EXPORT
void SCIPcaptureExpr(
   SCIP_EXPR*            expr                /**< expression to be captured */
   );

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_EXPORT
SCIP_RETCODE SCIPreleaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr                /**< pointer to expression to be released */
   );

/** returns whether an expression is a variable expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a value expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a sum expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a product expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns whether an expression is a power expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprPower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** print an expression as info-message */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be printed */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   );

/** initializes printing of expressions in dot format to a give FILE* pointer */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprDotInit(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata,        /**< buffer to store dot printing data */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_EXPRPRINT_WHAT     whattoprint       /**< info on what to print for each expression */
   );

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprDotInit2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata,        /**< buffer to store dot printing data */
   const char*             filename,         /**< name of file to print to */
   SCIP_EXPRPRINT_WHAT     whattoprint       /**< info on what to print for each expression */
   );

/** main part of printing an expression in dot format */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprDot(
   SCIP*                  scip,              /**< SCIP data structure */
   SCIP_EXPRPRINTDATA*    printdata,         /**< data as initialized by \ref SCIPprintExprDotInit() */
   SCIP_EXPR*             expr               /**< expression to be printed */
   );

/** finishes printing of expressions in dot format */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprDotFinal(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata         /**< buffer where dot printing data has been stored */
   );

/** shows a single expression by use of dot and gv
 *
 * This function is meant for debugging purposes.
 * It's signature is kept as simple as possible to make it
 * easily callable from gdb, for example.
 *
 * It prints the expression into a temporary file in dot format, then calls dot to create a postscript file, then calls ghostview (gv) to show the file.
 * SCIP will hold until ghostscript is closed.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPshowExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression to be printed */
   );

/** prints structure of an expression a la Maple's dismantle */
SCIP_EXPORT
SCIP_RETCODE SCIPdismantleExpr(
   SCIP*                 scip,               /**< SCIP data structure */
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
SCIP_EXPORT
SCIP_RETCODE SCIPevalExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** returns a previously unused solution tag for expression evaluation */
SCIP_EXPORT
SCIP_Longint SCIPgetExprNewSoltag(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

/** @name Differentiation
 * @anchor SCIP_EXPR_DIFF
 *
 * @par Gradients (Automatic differentiation Backward mode)
 *
 * Given a function, say, \f$f(s(x,y),t(x,y))\f$ there is a common mnemonic technique to compute its partial derivatives, using a tree diagram.
 * Suppose we want to compute the partial derivative of \f$f\f$ w.r.t. \f$x\f$.
 * Write the function as a tree:
 *
 *     f
 *     |-----|
 *     s     t
 *     |--|  |--|
 *     x  y  x  y
 *
 * The weight of an edge between two nodes represents the partial derivative of the parent w.r.t. the children, e.g.,
 *
 *     f
 *     |
 *     s
 *
 * is \f$ \partial_sf \f$.
 * The weight of a path is the product of the weight of the edges in the path.
 * The partial derivative of \f$f\f$ w.r.t. \f$x\f$ is then the sum of the weights of all paths connecting \f$f\f$ with \f$x\f$:
 * \f[ \frac{\partial f}{\partial x} = \partial_s f \cdot \partial_x s + \partial_t f \cdot \partial_x t. \f]
 *
 * We follow this method in order to compute the gradient of an expression (root) at a given point (point).
 * Note that an expression is a DAG representation of a function, but there is a 1-1 correspondence between paths
 * in the DAG and path in a tree diagram of a function.
 * Initially, we set `root->derivative` to 1.0.
 * Then, traversing the tree in Depth First (see \ref SCIPexpriterInit), for every expr that *has* children,
 * we store in its i-th child, `child[i]->derivative`, the derivative of expr w.r.t. child evaluated at point multiplied with `expr->derivative`.
 *
 * For example:
 * 1. `f->derivative` = 1.0
 * 2. `s->derivative` = \f$\partial_s f \,\cdot\f$ `f->derivative` = \f$\partial_s f\f$
 * 3. `x->derivative` = \f$\partial_x s \,\cdot\f$ `s->derivative` = \f$\partial_x s \cdot \partial_s f\f$
 *
 * However, when the child is a variable expressions, we actually need to initialize `child->derivative` to 0.0
 * and afterwards add, instead of overwrite the computed value.
 * The complete example would then be:
 *
 * 1. `f->derivative` = 1.0, `x->derivative` = 0.0, `y->derivative` = 0.0
 * 2. `s->derivative` =  \f$\partial_s f \,\cdot\f$ `f->derivative` = \f$\partial_s f\f$
 * 3. `x->derivative` += \f$\partial_x s \,\cdot\f$ `s->derivative` = \f$\partial_x s \cdot \partial_s f\f$
 * 4. `y->derivative` += \f$\partial_y s \,\cdot\f$ `s->derivative` = \f$\partial_y s \cdot \partial_s f\f$
 * 5. `t->derivative` =  \f$\partial_t f \,\cdot\f$ `f->derivative` = \f$\partial_t f\f$
 * 6. `x->derivative` += \f$\partial_x t \,\cdot\f$ `t->derivative` = \f$\partial_x t \cdot \partial_t f\f$
 * 7. `y->derivative` += \f$\partial_y t \,\cdot\f$ `t->derivative` = \f$\partial_y t \cdot \partial_t f\f$
 *
 * Note that, to compute this, we only need to know, for each expression, its partial derivatives w.r.t a given child at a point.
 * This is what the callback `SCIP_DECL_EXPRBWDIFF` should return.
 * Indeed, from "derivative of expr w.r.t. child evaluated at point multiplied with expr->derivative",
 * note that at the moment of processing a child, we already know `expr->derivative`, so the only
 * missing piece of information is "the derivative of expr w.r.t. child evaluated at point".
 *
 * An equivalent way of interpreting the procedure is that `expr->derivative` stores the derivative of the root w.r.t. expr.
 * This way, `x->derivative` and `y->derivative` will contain the partial derivatives of root w.r.t. the variable, that is, the gradient.
 * Note, however, that this analogy is only correct for leave expressions, since the derivative value of an intermediate expression gets overwritten.
 *
 *
 * \par Hessian (Automatic differentiation Backward on Forward mode)
 *
 * Computing the Hessian is more complicated since it is the derivative of the gradient, which is a function with more than one output.
 * We compute the Hessian by computing "directions" of the Hessian, that is \f$H\cdot u\f$ for different \f$u\f$.
 * This is easy in general, since it is the gradient of the *scalar* function \f$\nabla f u\f$, that is,
 * the directional derivative of \f$f\f$ in the direction \f$u\f$: \f$D_u f\f$.
 *
 * This is easily computed via the so called forward mode.
 * Just as `expr->derivative` stores the partial derivative of the root w.r.t. expr,
 * `expr->dot` stores the directional derivative of expr in the direction \f$u\f$.
 * Then, by the chain rule, `expr->dot` = \f$\sum_{c:\text{children}} \partial_c \text{expr} \,\cdot\f$ `c->dot`.
 *
 * Starting with `x[i]->dot` = \f$u_i\f$, we can compute `expr->dot` for every expression at the same time we evaluate expr.
 * Computing `expr->dot` is the purpose of the callback `SCIP_DECL_EXPRFWDIFF`.
 * Obviously, when this callback is called, the "dots" of all children are known
 * (just like evaluation, where the value of all children are known).
 *
 * Once we have this information, we compute the gradient of this function, following the same idea as before.
 * We define `expr->bardot` to be the directional derivative in direction \f$u\f$ of the partial derivative of the root w.r.t `expr`,
 * that is \f$D_u (\partial_{\text{expr}} f) = D_u\f$ (`expr->derivative`).
 *
 * This way, `x[i]->bardot` = \f$D_u (\partial_{x_i} f) = e_i^T H_f u\f$.
 * Hence `vars->bardot` contain \f$H_f u\f$.
 * By the chain rule, product rule, and definition we have
 * \f{eqnarray*}{
 * \texttt{expr->bardot} & = & D_u (\partial_{\text{expr}} f) \\
 *   & = & D_u ( \partial_{\text{parent}} f \cdot \partial_{\text{expr}} \text{parent} )  \\
 *   & = & D_u ( \texttt{parent->derivative} \cdot \partial_{\text{expr}} \text{parent} ) \\
 *   & = & \partial_{\text{expr}} \text{parent} \cdot D_u (\texttt{parent->derivative}) + \texttt{parent->derivative} \cdot D_u (\partial_{\text{expr}} \text{parent}) \\
 *   & = & \texttt{parent->bardot} \cdot \partial_{\text{expr}} \text{parent} + \texttt{parent->derivative} \cdot D_u (\partial_{\text{expr}} \text{parent})
 * \f}
 *
 * Note that we have computed `parent->bardot` and `parent->derivative` at this point,
 * while \f$\partial_{\text{expr}} \text{parent}\f$ is the return of `SCIP_DECL_EXPRBWDIFF`.
 * Hence the only information we need to compute is \f$D_u (\partial_{\text{expr}} \text{parent})\f$.
 * This is the purpose of the callback `SCIP_DECL_EXPRBWFWDIFF`.
 *
 * @{
 */

/** evaluates gradient of an expression for a given point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalExprGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   SCIP_SOL*             sol,                /**< solution to be evaluated (NULL for the current LP solution) */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   );

/** evaluates Hessian-vector product of an expression for a given point and direction
 *
 * Evaluates children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffGradientDirNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalExprHessianDir(
   SCIP*                 scip,             /**< SCIP data structure */
   SCIP_EXPR*            expr,             /**< expression to be differentiated */
   SCIP_SOL*             sol,              /**< solution to be evaluated (NULL for the current LP solution) */
   SCIP_Longint          soltag,           /**< tag that uniquely identifies the solution (with its values), or 0. */
   SCIP_SOL*             direction         /**< direction */
   );

/**@} */  /* end of differentiation methods */

/**@name Expressions
 * @{
 */

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is no longer uptodate (some bound was changed since last evaluation).
 *
 * The owner of the expression may overwrite the methods used to evaluate the activity,
 * including whether the local or global domain of variables is used.
 * By default (no owner, or owner doesn't overwrite activity evaluation),
 * the local domain of variables is used.
 *
 * @note If expression is set to be integral, then activities are tightened to integral values.
 *   Thus, ensure that the integrality information is valid (if set to TRUE; the default (FALSE) is always ok).
 */
SCIP_EXPORT
SCIP_RETCODE SCIPevalExprActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** compare expressions
 * @return -1, 0 or 1 if expr1 <, =, > expr2, respectively
 * @note The given expressions are assumed to be simplified.
 */
SCIP_EXPORT
int SCIPcompareExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr1,              /**< first expression */
   SCIP_EXPR*            expr2               /**< second expression */
   );

/** compute the hash value of an expression */
SCIP_EXPORT
SCIP_RETCODE SCIPhashExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   unsigned int*         hashval             /**< pointer to store the hash value */
   );

/** simplifies an expression
 *
 * This is largely inspired by Joel Cohen's
 * *Computer algebra and symbolic computation: Mathematical methods*,
 * in particular Chapter 3.
 * The other fountain of inspiration are the simplifying methods of expr.c in SCIP 7.
 *
 * Note: The things to keep in mind when adding simplification rules are the following.
 * I will be using the product expressions (see expr_product.c) as an example.
 * There are mainly 3 parts of the simplification process. You need to decide
 * at which stage the simplification rule makes sense.
 * 1. Simplify each factor (simplifyFactor()): At this stage we got the children of the product expression.
 *    At this point, each child is simplified when viewed as a stand-alone expression, but not necessarily when viewed as child of a product expression.
 *    Rules like SP2, SP7, etc are enforced at this point.
 * 2. Multiply the factors (mergeProductExprlist()): At this point rules like SP4, SP5 and SP14 are enforced.
 * 3. Build the actual simplified product expression (buildSimplifiedProduct()):
 *    At this point rules like SP10, SP11, etc are enforced.
 *
 * During steps 1 and 2 do not forget to set the flag `changed` to TRUE when something actually changes.
 *
 * \par Definition of simplified expressions
 *
 * An expression is simplified if it
 * - is a value expression
 * - is a var expression
 * - is a product expression such that
 *   - SP1: every child is simplified
 *   - SP2: no child is a product
 *   - SP4: no two children are the same expression (those should be multiplied)
 *   - SP5: the children are sorted [commutative rule]
 *   - SP7: no child is a value
 *   - SP8: its coefficient is 1.0 (otherwise should be written as sum)
 *   - SP10: it has at least two children
 *   - TODO?: at most one child is an `abs`
 *   - SP11: no two children are `expr*log(expr)`
 *     (TODO: we could handle more complicated stuff like \f$xy\log(x) \to - y * \mathrm{entropy}(x)\f$, but I am not sure this should happen at the simplification level;
 *            similar for \f$(xy) \log(xy)\f$, which currently simplifies to \f$xy \log(xy)\f$)
 *   - SP12: if it has two children, then neither of them is a sum (expand sums)
 *   - SP13: no child is a sum with a single term
 *   - SP14: at most one child is an `exp`
 * - is a power expression such that
 *   - POW1: exponent is not 0
 *   - POW2: exponent is not 1
 *   - POW3: its child is not a value
 *   - POW4: its child is simplified
 *   - POW5: if exponent is integer, its child is not a product
 *   - POW6: if exponent is integer, its child is not a sum with a single term (\f$(2x)^2 \to 4x^2\f$)
 *   - POW7: if exponent is 2, its child is not a sum (expand sums)
 *   - POW8: its child is not a power unless \f$(x^n)^m\f$ with \f$nm\f$ being integer and \f$n\f$ or \f$m\f$ fractional and \f$n\f$ not being even integer
 *   - POW9: its child is not a sum with a single term with a positive coefficient: \f$(25x)^{0.5} \to 5 x^{0.5}\f$
 *   - POW10: its child is not a binary variable: \f$b^e, e > 0 \to b\f$; \f$b^e, e < 0 \to b := 1\f$
 *   - POW11: its child is not an exponential: \f$\exp(\text{expr})^e \to \exp(e\cdot\text{expr})\f$
 * - is a signedpower expression such that
 *   - SPOW1: exponent is not 0
 *   - SPOW2: exponent is not 1
 *   - SPOW3: its child is not a value
 *   - SPOW4: its child is simplified
 *   - SPOW5: (TODO) do we want to distribute signpowers over products like we do for powers?
 *   - SPOW6: exponent is not an odd integer: (signpow odd expr) -> (pow odd expr)
 *   - SPOW8: if exponent is integer, its child is not a power
 *   - SPOW9: its child is not a sum with a single term: \f$\mathrm{signpow}(25x,0.5) \to 5\mathrm{signpow}(x,0.5)\f$
 *   - SPOW10: its child is not a binary variable: \f$\mathrm{signpow}(b,e), e > 0 \to b\f$; \f$\mathrm{signpow}(b,e), e < 0 \to b := 1\f$
 *   - SPOW11: its child is not an exponential: \f$\mathrm{signpow}(\exp(\text{expr}),e) \to \exp(e\cdot\text{expr})\f$
 *   - TODO: what happens when child is another signed power?
 *   - TODO: if child &ge; 0 -> transform to normal power; if child < 0 -> transform to - normal power
 *
 *   TODO: Some of these criteria are too restrictive for signed powers; for example, the exponent does not need to be
 *   an integer for signedpower to distribute over a product (SPOW5, SPOW6, SPOW8). Others can also be improved.
 * - is a sum expression such that
 *   - SS1: every child is simplified
 *   - SS2: no child is a sum
 *   - SS3: no child is a value (values should go in the constant of the sum)
 *   - SS4: no two children are the same expression (those should be summed up)
 *   - SS5: the children are sorted [commutative rule]
 *   - SS6: it has at least one child
 *   - SS7: if it consists of a single child, then either constant is != 0.0 or coef != 1
 *   - SS8: no child has coefficient 0
 *   - SS9: if a child c is a product that has an exponential expression as one of its factors, then the coefficient of c is +/-1.0
 *   - SS10: if a child c is an exponential, then the coefficient of c is +/-1.0
 * - it is a function with simplified arguments, but not all of them can be values
 * - TODO? a logarithm doesn't have a product as a child
 * - TODO? the exponent of an exponential is always 1
 *
 * \par Ordering Rules (see SCIPexprCompare())
 * \anchor EXPR_ORDER
 * These rules define a total order on *simplified* expressions.
 * There are two groups of rules, when comparing equal type expressions and different type expressions.
 *
 * Equal type expressions:
 * - OR1: u,v value expressions: u < v &hArr; val(u) < val(v)
 * - OR2: u,v var expressions: u < v &hArr; `SCIPvarGetIndex(var(u))` < `SCIPvarGetIndex(var(v))`
 * - OR3: u,v are both sum or product expression: < is a lexicographical order on the terms
 * - OR4: u,v are both pow: u < v &hArr; base(u) < base(v) or, base(u) = base(v) and expo(u) < expo(v)
 * - OR5: u,v are \f$u = f(u_1, ..., u_n), v = f(v_1, ..., v_m)\f$: u < v &hArr; For the first k such that \f$u_k \neq v_k\f$, \f$u_k < v_k\f$, or if such a \f$k\f$ doesn't exist, then \f$n < m\f$.
 *
 * Different type expressions:
 * - OR6: u value, v other: u < v always
 * - OR7: u sum, v var or func: u < v &hArr; u < 0+v;
 *        In other words, if \f$u = \sum_{i=1}^n \alpha_i u_i\f$, then u < v &hArr; \f$u_n\f$ < v or if \f$u_n\f$ = v and \f$\alpha_n\f$ < 1.
 * - OR8: u product, v pow, sum, var or func: u < v &hArr; u < 1*v;
 *        In other words, if \f$u = \prod_{i=1}^n u_i\f$, then u < v &hArr; \f$u_n\f$ < v.
 *        Note: since this applies only to simplified expressions, the form of the product is correct.
 *              Simplified products  do *not* have constant coefficients.
 * - OR9: u pow, v sum, var or func: u < v &hArr; u < v^1
 * - OR10: u var, v func: u < v always
 * - OR11: u func, v other type of func: u < v &hArr; name(type(u)) < name(type(v))
 * - OR12: none of the rules apply: u < v &hArr; ! v < u
 *
 * Examples:
 * - x < x^2 ?:  x is var and x^2 power, so none applies (OR12).
 *   Hence, we try to answer x^2 < x ?: x^2 < x &hArr; x < x or if x = x and 2 < 1 &hArr; 2 < 1 &hArr; False. So x < x^2 is True.
 * - x < x^-1 --OR12&rarr; ~(x^-1 < x) --OR9&rarr; ~(x^-1 < x^1) --OR4&rarr; ~(x < x or -1 < 1) &rarr; ~True &rarr; False
 * - x*y < x --OR8&rarr; x*y < 1*x --OR3&rarr; y < x --OR2&rarr; False
 * - x*y < y --OR8&rarr; x*y < 1*y --OR3&rarr; y < x --OR2&rarr; False
 *
 * \par Algorithm
 *
 * The recursive version of the algorithm is
 *
 *     EXPR simplify(expr)
 *        for c in 1..expr->nchildren
 *           expr->children[c] = simplify(expr->children[c])
 *        end
 *        return expr->exprhdlr->simplify(expr)
 *     end
 *
 * Important: Whatever is returned by a simplify callback **has** to be simplified.
 * Also, all children of the given expression **are** already simplified.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsimplifyExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            rootexpr,           /**< expression to be simplified */
   SCIP_EXPR**           simplified,         /**< buffer to store simplified expression */
   SCIP_Bool*            changed,            /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*            infeasible,         /**< buffer to store whether infeasibility has been detected */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** replaces common sub-expressions in a given expression graph by using a hash key for each expression
 *
 *  The algorithm consists of two steps:
 *
 *  1. traverse through all given expressions and compute for each of them a (not necessarily unique) hash
 *
 *  2. initialize an empty hash table and traverse through all expression; check for each of them if we can find a
 *     structural equivalent expression in the hash table; if yes we replace the expression by the expression inside the
 *     hash table, otherwise we add it to the hash table
 *
 *  @note the hash keys of the expressions are used for the hashing inside the hash table; to compute if two expressions
 *  (with the same hash) are structurally the same we use the function SCIPexprCompare().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPreplaceCommonSubexpressions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           exprs,              /**< expressions (possibly replaced by equivalent on output) */
   int                   nexprs,             /**< total number of expressions */
   SCIP_Bool*            replacedroot        /**< buffer to store whether any root expression (expression in exprs) was replaced */
);

/** computes the curvature of a given expression and all its subexpressions
 *
 *  @note this function also evaluates all subexpressions w.r.t. current variable bounds
 *  @note this function relies on information from the curvature callback of expression handlers only,
 *    consider using function @ref SCIPhasExprCurvature() of the convex-nlhdlr instead, as that uses more information to deduce convexity
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** computes integrality information of a given expression and all its subexpressions
 *
 * The integrality information can be accessed via SCIPexprIsIntegral().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeExprIntegrality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the total number of variable expressions in an expression
 *
 * The function counts variable expressions in common sub-expressions only once, but
 * counts variables appearing in several variable expressions multiple times.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetExprNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   int*                  nvars               /**< buffer to store the total number of variables */
   );

/** returns all variable expressions contained in a given expression
 *
 * The array to store all variable expressions needs to be at least of size
 * the number of unique variable expressions in the expression which is given by SCIPgetExprNVars().
 *
 * If every variable is represented by only one variable expression (common subexpression have been removed)
 * then SCIPgetExprNVars() can be bounded by SCIPgetNTotalVars().
 * If, in addition, non-active variables have been removed from the expression, e.g., by simplifying,
 * then SCIPgetExprNVars() can be bounded by SCIPgetNVars().
 *
 * @note function captures variable expressions
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetExprVarExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPR**           varexprs,           /**< array to store all variable expressions */
   int*                  nvarexprs           /**< buffer to store the total number of variable expressions */
   );

/** @} */

/**@name Expression Handler Callbacks
 * @{
 */

/** calls the print callback for an expression
 *
 * @see SCIP_DECL_EXPRPRINT
 */
SCIP_EXPORT
SCIP_DECL_EXPRPRINT(SCIPcallExprPrint);

/** calls the curvature callback for an expression
 *
 * @see SCIP_DECL_EXPRCURVATURE
 *
 * Returns unknown curvature if callback not implemented.
 */
SCIP_EXPORT
SCIP_DECL_EXPRCURVATURE(SCIPcallExprCurvature);

/** calls the monotonicity callback for an expression
 *
 * @see SCIP_DECL_EXPRMONOTONICITY
 *
 * Returns unknown monotonicity if callback not implemented.
 */
SCIP_EXPORT
SCIP_DECL_EXPRMONOTONICITY(SCIPcallExprMonotonicity);

/** calls the eval callback for an expression with given values for children
 *
 * Does not iterates over expressions, but requires values for children to be given.
 * Value is not stored in expression, but returned in `val`.
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to `SCIP_INVALID`.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcallExprEval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_Real*            childrenvalues,     /**< values for children */
   SCIP_Real*            val                 /**< buffer to store evaluated value */
   );

/** calls the eval and fwdiff callback of an expression with given values for children
 *
 * Does not iterates over expressions, but requires values for children and direction to be given.
 *
 * Value is not stored in expression, but returned in `val`.
 * If an evaluation error (division by zero, ...) occurs, this value will be set to `SCIP_INVALID`.
 *
 * Direction is not stored in expression, but returned in `dot`.
 * If an differentiation error (division by zero, ...) occurs, this value will be set to `SCIP_INVALID`.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcallExprEvalFwdiff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_Real*            childrenvalues,     /**< values for children */
   SCIP_Real*            direction,          /**< direction in which to differentiate */
   SCIP_Real*            val,                /**< buffer to store evaluated value */
   SCIP_Real*            dot                 /**< buffer to store derivative value */
   );

/** calls the interval evaluation callback for an expression
 *
 * @see SCIP_DECL_EXPRINTEVAL
 *
 * Returns entire interval if callback not implemented.
 */
SCIP_EXPORT
SCIP_DECL_EXPRINTEVAL(SCIPcallExprInteval);

/** calls the estimate callback for an expression
 *
 * @see SCIP_DECL_EXPRESTIMATE
 *
 * Returns without success if callback not implemented.
 */
SCIP_EXPORT
SCIP_DECL_EXPRESTIMATE(SCIPcallExprEstimate);

/** calls the initial estimators callback for an expression
 *
 * @see SCIP_DECL_EXPRINITESTIMATES
 *
 * Returns no estimators if callback not implemented.
 */
SCIP_EXPORT
SCIP_DECL_EXPRINITESTIMATES(SCIPcallExprInitestimates);

/** calls the simplify callback for an expression
 *
 * @see SCIP_DECL_EXPRSIMPLIFY
 *
 * Returns unmodified expression if simplify callback not implemented.
 *
 * Does not simplify descendants (children, etc). Use SCIPsimplifyExpr() for that.
 */
SCIP_EXPORT
SCIP_DECL_EXPRSIMPLIFY(SCIPcallExprSimplify);

/** calls the reverse propagation callback for an expression
 *
 * @see SCIP_DECL_EXPRREVERSEPROP
 *
 * Returns unmodified `childrenbounds` if reverseprop callback not implemented.
 */
SCIP_EXPORT
SCIP_DECL_EXPRREVERSEPROP(SCIPcallExprReverseprop);

#ifdef NDEBUG
#define SCIPappendExprChild(scip, expr, child)               SCIPexprAppendChild((scip)->set, (scip)->mem->probmem, expr, child)
#define SCIPreplaceExprChild(scip, expr, childidx, newchild) SCIPexprReplaceChild((scip)->set, (scip)->stat, (scip)->mem->probmem, expr, childidx, newchild)
#define SCIPremoveExprChildren(scip, expr)                   SCIPexprRemoveChildren((scip)->set, (scip)->stat, (scip)->mem->probmem, expr)
#define SCIPduplicateExpr(scip, expr, copyexpr, mapexpr, mapexprdata, ownercreate, ownercreatedata) SCIPexprCopy((scip)->set, (scip)->stat, (scip)->mem->probmem, (scip)->set, (scip)->stat, (scip)->mem->probmem, expr, copyexpr, mapexpr, mapexprdata, ownercreate, ownercreatedata)
#define SCIPduplicateExprShallow(scip, expr, copyexpr, ownercreate, ownercreatedata) SCIPexprDuplicateShallow((scip)->set, (scip)->mem->probmem, expr, copyexpr, ownercreate, ownercreatedata)
#define SCIPcaptureExpr(expr)                                SCIPexprCapture(expr)
#define SCIPreleaseExpr(scip, expr)                          SCIPexprRelease((scip)->set, (scip)->stat, (scip)->mem->probmem, expr)
#define SCIPisExprVar(scip, expr)                            SCIPexprIsVar((scip)->set, expr)
#define SCIPisExprValue(scip, expr)                          SCIPexprIsValue((scip)->set, expr)
#define SCIPisExprSum(scip, expr)                            SCIPexprIsSum((scip)->set, expr)
#define SCIPisExprProduct(scip, expr)                        SCIPexprIsProduct((scip)->set, expr)
#define SCIPisExprPower(scip, expr)                          SCIPexprIsPower((scip)->set, expr)
#define SCIPprintExpr(scip, expr, file)                      SCIPexprPrint((scip)->set, (scip)->stat, (scip)->mem->probmem, (scip)->messagehdlr, file, expr)
#define SCIPevalExpr(scip, expr, sol, soltag)                SCIPexprEval((scip)->set, (scip)->stat, (scip)->mem->probmem, expr, sol, soltag)
#define SCIPgetExprNewSoltag(scip)                           (++((scip)->stat->exprlastsoltag))
#define SCIPevalExprGradient(scip, expr, sol, soltag)        SCIPexprEvalGradient((scip)->set, (scip)->stat, (scip)->mem->probmem, expr, sol, soltag)
#define SCIPevalExprHessianDir(scip, expr, sol, soltag, direction) SCIPexprEvalHessianDir((scip)->set, (scip)->stat, (scip)->mem->probmem, expr, sol, soltag, direction)
#define SCIPevalExprActivity(scip, expr)                     SCIPexprEvalActivity((scip)->set, (scip)->stat, (scip)->mem->probmem, expr)
#define SCIPcompareExpr(scip, expr1, expr2)                  SCIPexprCompare((scip)->set, expr1, expr2)
#define SCIPsimplifyExpr(scip, rootexpr, simplified, changed, infeasible, ownercreate, ownercreatedata) SCIPexprSimplify((scip)->set, (scip)->stat, (scip)->mem->probmem, rootexpr, simplified, changed, infeasible, ownercreate, ownercreatedata)
#define SCIPcallExprCurvature(scip, expr, exprcurvature, success, childcurv) SCIPexprhdlrCurvatureExpr(SCIPexprGetHdlr(expr), (scip)->set, expr, exprcurvature, success, childcurv)
#define SCIPcallExprMonotonicity(scip, expr, childidx, result) SCIPexprhdlrMonotonicityExpr(SCIPexprGetHdlr(expr), (scip)->set, expr, childidx, result)
#define SCIPcallExprEval(scip, expr, childrenvalues, val)    SCIPexprhdlrEvalExpr(SCIPexprGetHdlr(expr), (scip)->set, (scip)->mem->buffer, expr, val, childrenvalues, NULL)
#define SCIPcallExprEvalFwdiff(scip, expr, childrenvalues, direction, val, dot) SCIPexprhdlrEvalFwDiffExpr(SCIPexprGetHdlr(expr), (scip)->set, (scip)->mem->buffer, expr, val, dot, childrenvalues, NULL, direction, NULL)
#define SCIPcallExprInteval(scip, expr, interval, intevalvar, intevalvardata) SCIPexprhdlrIntEvalExpr(SCIPexprGetHdlr(expr), (scip)->set, expr, interval, intevalvar, intevalvardata)
#define SCIPcallExprEstimate(scip, expr, localbounds, globalbounds, refpoint, overestimate, targetvalue, coefs, constant, islocal, success, branchcand) SCIPexprhdlrEstimateExpr(SCIPexprGetHdlr(expr), (scip)->set, expr, localbounds, globalbounds, refpoint, overestimate, targetvalue, coefs, constant, islocal, success, branchcand)
#define SCIPcallExprInitestimates(scip, expr, bounds, overestimate, coefs, constant, nreturned)  SCIPexprhdlrInitEstimatesExpr(SCIPexprGetHdlr(expr), (scip)->set, expr, bounds, overestimate, coefs, constant, nreturned)
#define SCIPcallExprSimplify(scip, expr, simplifiedexpr, ownercreate, ownercreatedata)  SCIPexprhdlrSimplifyExpr(SCIPexprGetHdlr(expr), (scip)->set, expr, simplifiedexpr, ownercreate, ownercreatedata)
#define SCIPcallExprReverseprop(scip, expr, bounds, childrenbounds, infeasible) SCIPexprhdlrReversePropExpr(SCIPexprGetHdlr(expr), (scip)->set, expr, bounds, childrenbounds, infeasible)
#endif

/** @} */


/**@name Expression Iterator */
/**@{ */

/** creates an expression iterator */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExpriter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRITER**       iterator            /**< buffer to store expression iterator */
   );

/** frees an expression iterator */
SCIP_EXPORT
void SCIPfreeExpriter(
   SCIP_EXPRITER**       iterator            /**< pointer to the expression iterator */
   );

#ifdef NDEBUG
#define SCIPcreateExpriter(scip, iterator)  SCIPexpriterCreate((scip)->stat, (scip)->mem->probmem, iterator)
#define SCIPfreeExpriter(iterator)          SCIPexpriterFree(iterator)
#endif

/** @} */


/**@name Quadratic Expressions */
/**@{ */

/** checks whether an expression is quadratic
 *
 * An expression is quadratic if it is either a square (of some expression), a product (of two expressions),
 * or a sum of terms where at least one is a square or a product.
 *
 * Use SCIPexprGetQuadraticData() to get data about the representation as quadratic.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcheckExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool*            isquadratic         /**< buffer to store result */
   );

/** frees information on quadratic representation of an expression
 *
 * Before doing changes to an expression, it can be useful to call this function.
 */
SCIP_EXPORT
void SCIPfreeExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** evaluates quadratic term in a solution
 *
 * \note This requires that every expressiion used in the quadratic data is a variable expression.
 */
SCIP_EXPORT
SCIP_Real SCIPevalExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL for LP solution */
   );

/** prints quadratic expression */
SCIP_EXPORT
SCIP_RETCODE SCIPprintExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< quadratic expression */
   );

/** checks the curvature of the quadratic expression
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
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeExprQuadraticCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_EXPRCURV*        curv,               /**< pointer to store the curvature of quadratics */
   SCIP_HASHMAP*         assumevarfixed,     /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   SCIP_Bool             storeeigeninfo      /**< whether the eigenvalues and eigenvectors should be stored */
   );

#ifdef NDEBUG
#define SCIPcheckExprQuadratic(scip, expr, isquadratic)  SCIPexprCheckQuadratic((scip)->set, (scip)->mem->probmem, expr, isquadratic)
#define SCIPfreeExprQuadratic(scip, expr)                SCIPexprFreeQuadratic((scip)->mem->probmem, expr)
#define SCIPcomputeExprQuadraticCurvature(scip, expr, curv, assumevarfixed, storeeigeninfo)  SCIPexprComputeQuadraticCurvature((scip)->set, (scip)->mem->probmem, (scip)->mem->buffer, (scip)->messagehdlr, expr, curv, assumevarfixed, storeeigeninfo)
#endif

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* SCIP_SCIP_EXPR_H_ */
