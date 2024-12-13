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

/**@file   cons_nonlinear.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for nonlinear constraints specified by algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 *
 * For additional documentation on this constraint handler, see also the SCIP 8 release report.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_NONLINEAR_H__
#define __SCIP_CONS_NONLINEAR_H__


#include "scip/scip.h"
#include "scip/type_nlhdlr.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup CONSHDLRS
 * @{
 *
 * @name Nonlinear Constraints
 * @{
 */

/** linear auxiliary expression of the form xy {&le;,&ge;,=} coefs[0]w + coefs[1]x + coefs[2]y + cst */
struct SCIP_ConsNonlinear_Auxexpr
{
   SCIP_Real             coefs[3];           /**< coefficients in the auxiliary expression */
   SCIP_Real             cst;                /**< constant */
   SCIP_VAR*             auxvar;             /**< auxiliary variable w in xy {&le;,&ge;,=} auxexpr(w, x, y) */
   SCIP_Bool             underestimate;      /**< whether the auxexpr underestimates the product */
   SCIP_Bool             overestimate;       /**< whether the auxexpr overestimates the product */
};
typedef struct SCIP_ConsNonlinear_Auxexpr SCIP_CONSNONLINEAR_AUXEXPR;

/** bilinear term structure
 *
 * This can represent a product which
 * - explicitly exists in the problem and is under- and/or overestimated by a single auxiliary variable
 *   stored as `var` in the union `aux` (case `nauxexprs` = 0) or
 * - is involved in bilinear relations implicitly given by linear constraints with binary variables, and
 *   is under- and/or overestimated by linear expression(s) stored as `exprs` in the union `aux` (case `nauxexprs` > 0).
 *
 * An explicitly existing product can also be involved in implicit relations, then it will be stored as in
 * the second case.
 */
struct SCIP_ConsNonlinear_BilinTerm
{
   SCIP_VAR*             x;                  /**< first variable */
   SCIP_VAR*             y;                  /**< second variable */
   union
   {
      SCIP_CONSNONLINEAR_AUXEXPR** exprs;    /**< auxiliary expressions for the implicit product of x and y */
      SCIP_VAR*          var;                /**< auxiliary variable for the explicit product of x and y */
   }                     aux;
   int                   nauxexprs;          /**< number of aux.exprs (0 for products without implicit relations) */
   int                   auxexprssize;       /**< size of the aux.exprs array */
   int                   nlockspos;          /**< number of positive expression locks */
   int                   nlocksneg;          /**< number of negative expression locks */
   SCIP_Bool             existing;           /**< does the product exist explicitly in the problem? */
};
typedef struct SCIP_ConsNonlinear_BilinTerm SCIP_CONSNONLINEAR_BILINTERM; /**< bilinear term structure */

/** evaluation callback for (vertex-polyhedral) functions used as input for facet computation of its envelopes
 *
 * \param[in] args     the point to be evaluated
 * \param[in] nargs    the number of arguments of the function (length of array `args`)
 * \param[in] funcdata user-data of function evaluation callback
 * \return value of function in point given by `args` or SCIP_INVALID if could not be evaluated
 */
#define SCIP_DECL_VERTEXPOLYFUN(f) SCIP_Real f (SCIP_Real* args, int nargs, void* funcdata)

/** maximum dimension of vertex-polyhedral function for which we can try to compute a facet of its convex or concave envelope */
#define SCIP_MAXVERTEXPOLYDIM 14

/** upgrading method for nonlinear constraints into more specific constraints
 *
 * The method might upgrade a nonlinear constraint into a set of upgrade constraints.
 * The caller provided an array `upgdconss` of size `upgdconsssize` to store upgrade constraints.
 * If an upgrade is not possible, set `*nupgdconss` to zero.
 * If more than `upgdconsssize` many constraints shall replace `cons`, the function
 * should return the required number as negated value in `*nupgdconss`,
 * e.g., if `cons` should be replaced by 3 constraints, the function should set
 * `*nupgdconss` to -3 and return with SCIP_OKAY.
 *
 * \param[in] scip           SCIP main data structure
 * \param[in] cons           the nonlinear constraint to upgrade
 * \param[in] nvarexprs      total number of variable expressions in the nonlinear constraint
 * \param[out] nupgdconss    pointer to store number of constraints that replace this constraint
 * \param[out] upgdconss     array to store constraints that replace this constraint
 * \param[in] upgdconsssize  length of the provided `upgdconss` array
 */
#define SCIP_DECL_NONLINCONSUPGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS* cons, int nvarexprs, \
      int* nupgdconss, SCIP_CONS** upgdconss, int upgdconsssize)

/** @} */
/** @} */

/** creates the handler for nonlinear constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Nonlinear Constraints
 *
 * @{
 */

/* Nonlinear Constraint Handler Methods */

/** includes a nonlinear constraint upgrade method into the nonlinear constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsUpgradeNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_NONLINCONSUPGD((*nlconsupgd)),  /**< method to call for upgrading nonlinear constraint */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   );

/** creates and captures a nonlinear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*            expr,               /**< expression of constraint (must not be NULL) */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   );

/** creates and captures a nonlinear constraint with all its constraint flags set to their default values
 *
 *  All flags can be set via SCIPconsSetFLAGNAME-methods.
 *
 *  @see SCIPcreateConsNonlinear() for information about the basic constraint flag configuration.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*            expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** creates and captures a quadratic nonlinear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsQuadraticNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation */
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
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   );

/** creates and captures a quadratic nonlinear constraint with all its constraint flags set to their default values
 *
 *  All flags can be set via SCIPconsSetFLAGNAME-methods.
 *
 *  @see SCIPcreateConsQuadraticNonlinear() for information about the basic constraint flag configuration.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicQuadraticNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs                 /**< right hand side of quadratic equation */
   );

/** creates and captures a nonlinear constraint that is a second-order cone constraint with all its constraint flags set to their default values
 *
 * \f$\sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} \leq \alpha_{n+1}\, (x_{n+1}+\beta_{n+1})\f$
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicSOCNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables on left hand side of constraint (n) */
   SCIP_VAR**            vars,               /**< array with variables on left hand side (x_i) */
   SCIP_Real*            coefs,              /**< array with coefficients of left hand side variables (alpha_i), or NULL if all 1.0 */
   SCIP_Real*            offsets,            /**< array with offsets of variables (beta_i), or NULL if all 0.0 */
   SCIP_Real             constant,           /**< constant on left hand side (gamma) */
   SCIP_VAR*             rhsvar,             /**< variable on right hand side of constraint (x_{n+1}) */
   SCIP_Real             rhscoeff,           /**< coefficient of variable on right hand side (alpha_{n+1}) */
   SCIP_Real             rhsoffset           /**< offset of variable on right hand side (beta_{n+1}) */
   );

/** creates and captures a signpower nonlinear constraint with all its constraint flags set to their default values
 *
 * \f$\textrm{lhs} \leq \textrm{sign}(x+a) |x+a|^n + c z \leq \textrm{rhs}\f$
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicSignpowerNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             x,                  /**< nonlinear variable x in constraint */
   SCIP_VAR*             z,                  /**< linear variable z in constraint */
   SCIP_Real             exponent,           /**< exponent n of |x+offset|^n term in constraint */
   SCIP_Real             xoffset,            /**< offset in |x+offset|^n term in constraint */
   SCIP_Real             zcoef,              /**< coefficient of z in constraint */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** gets tag indicating current local variable bounds */
SCIP_EXPORT
SCIP_Longint SCIPgetCurBoundsTagNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   );

/** gets the `curboundstag` from the last time where variable bounds were relaxed */
SCIP_EXPORT
SCIP_Longint SCIPgetLastBoundRelaxTagNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   );

/** increments `curboundstag` and resets `lastboundrelax` in constraint handler data
 *
 * @attention This method is not intended for normal use.
 *   These tags are maintained by the event handler for variable bound change events.
 *   This method is used by some unittests.
 */
SCIP_EXPORT
void SCIPincrementCurBoundsTagNonlinear(
   SCIP_CONSHDLR*          conshdlr,         /**< nonlinear constraint handler */
   SCIP_Bool               boundrelax        /**< indicates whether a bound was relaxed, i.e., lastboundrelax should be set too */
   );

/** returns the hashmap that is internally used to map variables to their corresponding variable expressions */
SCIP_EXPORT
SCIP_HASHMAP* SCIPgetVarExprHashmapNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   );

/** processes a rowprep for cut addition and maybe report branchscores */
SCIP_EXPORT
SCIP_RETCODE SCIPprocessRowprepNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler which provided the estimator */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_ROWPREP*         rowprep,            /**< cut to be added */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             auxvalue,           /**< current value of expression w.r.t. auxiliary variables as obtained from EVALAUX */
   SCIP_Bool             allowweakcuts,      /**< whether we should only look for "strong" cuts, or anything that separates is fine */
   SCIP_Bool             branchscoresuccess, /**< whether the estimator generation generated branching scores */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement, or only in separation */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result */
   );

/** returns whether all nonlinear constraints are assumed to be convex */
SCIP_EXPORT
SCIP_Bool SCIPassumeConvexNonlinear(
   SCIP_CONSHDLR*        conshdlr
   );

/** collects all bilinear terms for a given set of constraints
 *
 * @attention This method should only be used for unit tests that depend on SCIPgetBilinTermsNonlinear(),
 *       SCIPgetBilinTermNonlinear() or SCIPgetBilinTermIdxNonlinear().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcollectBilinTermsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_CONS**           conss,              /**< nonlinear constraints */
   int                   nconss              /**< total number of nonlinear constraints */
   );

/** returns the total number of bilinear terms that are contained in all nonlinear constraints
 *
 *  @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 */
SCIP_EXPORT
int SCIPgetNBilinTermsNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   );

/** returns all bilinear terms that are contained in all nonlinear constraints
 *
 * @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @note The value of the auxiliary variable of a bilinear term might be NULL, which indicates that the term does not have an auxiliary variable.
 */
SCIP_EXPORT
SCIP_CONSNONLINEAR_BILINTERM* SCIPgetBilinTermsNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   );

/** returns the index of the bilinear term representing the product of the two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns -1 if the variables do not appear bilinearly.
 */
SCIP_EXPORT
int SCIPgetBilinTermIdxNonlinear(
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y                   /**< second variable */
   );

/** returns the bilinear term that represents the product of two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns NULL if the variables do not appear bilinearly.
 */
SCIP_EXPORT
SCIP_CONSNONLINEAR_BILINTERM* SCIPgetBilinTermNonlinear(
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y                   /**< second variable */
   );

/** evaluates an auxiliary expression for a bilinear term */
SCIP_EXPORT
SCIP_Real SCIPevalBilinAuxExprNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable of the bilinear term */
   SCIP_VAR*             y,                  /**< second variable of the bilinear term */
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr,      /**< auxiliary expression */
   SCIP_SOL*             sol                 /**< solution at which to evaluate (can be NULL) */
   );

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertBilinearTermExistingNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   int                   nlockspos,          /**< number of positive expression locks */
   int                   nlocksneg           /**< number of negative expression locks */
   );

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertBilinearTermImplicitNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   SCIP_Real             coefx,              /**< coefficient of x in the auxiliary expression */
   SCIP_Real             coefy,              /**< coefficient of y in the auxiliary expression */
   SCIP_Real             coefaux,            /**< coefficient of auxvar in the auxiliary expression */
   SCIP_Real             cst,                /**< constant of the auxiliary expression */
   SCIP_Bool             overestimate        /**< whether the auxiliary expression overestimates the bilinear product */
   );

/** computes a facet of the convex or concave envelope of a vertex polyhedral function
 *
 * If \f$ f(x) \f$ is vertex-polyhedral, then \f$ g \f$ is a convex underestimator if and only if
 * \f$ g(v^i) \leq f(v^i), \forall i \f$, where \f$ \{ v^i \}_{i = 1}^{2^n} \subseteq \mathbb R^n \f$ are the vertices
 * of the domain of \f$ x \f$, \f$ [\ell,u] \f$. Hence, we can compute a linear underestimator by solving the following
 * LP (we don't necessarily get a facet of the convex envelope, see below):
 *
 * \f{align*}{
 *              \max \, & \alpha^T x^* + \beta \\
 *     s.t. \; & \alpha^T v^i + \beta \le f(v^i), \, \forall i = 1, \ldots, 2^n
 * \f}
 *
 * In principle, one would need to update the LP whenever the domain changes. However, \f$ [\ell,u] = T([0, 1]^n) \f$,
 * where \f$ T \f$ is an affine linear invertible transformation given by \f$ T(y)_i = (u_i - \ell_i) y_i + \ell_i \f$.
 * Working with the change of variables \f$ x = T(y) \f$ allows us to keep the constraints of the LP, even if the domain
 * changes. Indeed, after the change of variables, the problem is: find an affine underestimator \f$ g \f$ such that \f$
 * g(T(y)) \le f(T(y)) \f$, for all \f$ y \in [0, 1]^n \f$. Now \f$ f(T(y)) \f$ is componentwise affine, but still
 * satisfies that \f$ g \f$ is a valid underestimator if and only if \f$ g(T(u)) \leq f(T(u)), \forall u \in \{0, 1\}^n
 * \f$. So we now look for \f$ \bar g(y) := g(T(y)) = g(((u_i - \ell_i) y_i + \ell_i)_i) = \bar \alpha^T y + \bar \beta
 * \f$, where \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i \f$ and \f$ \bar \beta = \sum_i \alpha_i \ell_i + \beta \f$. So
 * we find \f$ \bar g \f$ by solving the LP:
 *
 * \f{align*}{
 *              \max \, & \bar \alpha^T T^{-1}(x^*) + \bar \beta \\
 *     s.t. \; & \bar \alpha^T u + \bar \beta \le f(T(u)), \, \forall u \in \{0, 1\}^n
 * \f}
 *
 * and recover \f$ g \f$ by calculating \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i, \bar \beta = \sum_i \alpha_i \ell_i +
 * \beta \f$. Notice that \f$ f(T(u^i)) = f(v^i) \f$ so the right hand side doesn't change after the change of variables.
 *
 * Furthermore, the LP has more constraints than variables, so we solve its dual:
 * \f{align*}{
 *              \min \, & \sum_i \lambda_i f(v^i) \\
 *     s.t. \; & \sum_i \lambda_i u^i = T^{-1}(x^*) \\
 *             & \sum_i \lambda_i = 1 \\
 *             & \forall i, \, \lambda_i \geq 0
 * \f}
 *
 * In case we look for an overestimate, we do exactly the same, but have to maximize in the dual LP instead
 * of minimize.
 *
 * #### Technical and implementation details
 * -# \f$ U \f$ has exponentially many variables, so we only apply this separator for \f$n\f$ &le; \ref SCIP_MAXVERTEXPOLYDIM.
 * -# If the bounds are not finite, there is no underestimator. Also, \f$ T^{-1}(x^*) \f$ must be in the domain,
 * otherwise the dual is infeasible.
 * -# After a facet is computed, we check whether it is a valid facet (i.e. we check \f$ \alpha^T v + \beta \le f(v) \f$
 *  for every vertex \f$ v \f$). If we find a violation of at most ADJUSTFACETFACTOR * SCIPlpfeastol(), then we weaken \f$
 *  \beta \f$ by this amount, otherwise, we discard the cut.
 * -# If a variable is fixed within tolerances, we replace it with its value and compute the facet of the remaining
 * expression. Note that since we are checking the cut for validity, this will never produce wrong result.
 * -# If \f$ x^* \f$ is in the boundary of the domain, then the LP has infinitely many solutions, some of which might
 * have very bad numerical properties. For this reason, we perturb \f$ x^* \f$ to be in the interior of the region.
 * Furthermore, for some interior points, there might also be infinitely many solutions (e.g. for \f$ x y \f$ in \f$
 * [0,1]^2 \f$ any point \f$ (x^*, y^*) \f$ such that \f$ y^* = 1 - x^* \f$ has infinitely many solutions). For this
 * reason, we perturb any given \f$ x^* \f$. The idea is to try to get a facet of the convex/concave envelope. This only
 * happens when the solution has \f$ n + 1 \f$ non zero \f$ \lambda \f$'s (i.e. the primal has a unique solution).
 * -# We need to compute \f$ f(v^i) \f$ for every vertex of \f$ [\ell,u] \f$. A vertex is encoded by a number between 0
 * and \f$ 2^n - 1 \f$, via its binary representation (0 bit is lower bound, 1 bit is upper bound), so we can compute
 * all these values by iterating between 0 and \f$ 2^n - 1 \f$.
 * -# To check that the computed cut is valid we do the following: we use a gray code to loop over the vertices
 * of the box domain w.r.t. unfixed variables in order to evaluate the underestimator. To ensure the validity of the
 * underestimator, we check whether \f$ \alpha v^i + \beta \le f(v^i) \f$ for every vertex \f$ v^i \f$ and adjust
 * \f$ \beta \f$ if the maximal violation is small.
 *
 * @todo the solution is a facet if all variables of the primal have positive reduced costs (i.e. the solution is
 * unique). In the dual, this means that there are \f$ n + 1 \f$ variables with positive value. Can we use this or some
 * other information to handle any of both cases (point in the boundary or point in the intersection of polytopes
 * defining different pieces of the convex envelope)? In the case where the point is in the boundary, can we use that
 * information to maybe solve another to find a facet? How do the polytopes defining the pieces where the convex
 * envelope is linear looks like, i.e, given a point in the interior of a facet of the domain, does the midpoint of the
 * segment joining \f$ x^* \f$ with the center of the domain, always belongs to the interior of one of those polytopes?
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeFacetVertexPolyhedralNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_DECL_VERTEXPOLYFUN((*function)),     /**< pointer to vertex polyhedral function */
   void*                 fundata,            /**< data for function evaluation (can be NULL) */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
   );


/* Nonlinear Constraint Methods */

/** returns the expression of the given nonlinear constraint */
SCIP_EXPORT
SCIP_EXPR* SCIPgetExprNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the left hand side of a nonlinear constraint */
SCIP_EXPORT
SCIP_Real SCIPgetLhsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the right hand side of a nonlinear constraint */
SCIP_EXPORT
SCIP_Real SCIPgetRhsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the nonlinear constraint as a nonlinear row representation. */
SCIP_EXPORT
SCIP_RETCODE SCIPgetNlRowNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   );

/** returns the curvature of the expression of a given nonlinear constraint
 *
 * @note The curvature information is computed during CONSINITSOL.
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPgetCurvatureNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   );

/** checks whether expression of constraint can be represented as quadratic form
 *
 * Only sets `*isquadratic` to TRUE if the whole expression is quadratic (in the non-extended formulation) and non-linear.
 * That is, the expression in each \ref SCIP_QUADEXPR_QUADTERM will be a variable expressions and
 * \ref SCIPgetVarExprVar() can be used to retrieve the variable.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcheckQuadraticNonlinear(
   SCIP*                    scip,               /**< SCIP data structure */
   SCIP_CONS*               cons,               /**< constraint data */
   SCIP_Bool*               isquadratic         /**< buffer to store whether constraint is quadratic */
   );

/** changes left-hand-side of a nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgLhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left-hand-side */
   );

/** changes right-hand-side of a nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgRhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right-hand-side */
   );

/** changes expression of a nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgExprNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_EXPR*            expr                /**< new expression */
   );

/** adds coef * var to nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddLinearVarNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient */
   );

/** adds coef * expr to nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddExprNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             coef                /**< coefficient */
   );

/** gets absolute violation of nonlinear constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * If this value is at most SCIPfeastol(), the constraint would be considered feasible.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetAbsViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   );

/** gets scaled violation of nonlinear constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * The scaling that is applied to the absolute violation of the constraint
 * depends on the setting of parameter constraints/nonlinear/violscale.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetRelViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   );

/** returns a variable that appears linearly that may be decreased without making any other constraint infeasible */
SCIP_EXPORT
void SCIPgetLinvarMayDecreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   );

/** returns a variable that appears linearly that may be increased without making any other constraint infeasible */
SCIP_EXPORT
void SCIPgetLinvarMayIncreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   );


/* Methods for Expressions in Nonlinear Constraints
 * All functions in this group assume that the expression is owned by a the nonlinear constraint handler.
 */

/** returns the number of positive rounding locks of an expression */
SCIP_EXPORT
int SCIPgetExprNLocksPosNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the number of negative rounding locks of an expression */
SCIP_EXPORT
int SCIPgetExprNLocksNegNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the variable used for linearizing a given expression (return value might be NULL)
 *
 * @note for variable expression it returns the corresponding variable
 */
SCIP_EXPORT
SCIP_VAR* SCIPgetExprAuxVarNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the number of enforcements for an expression */
SCIP_EXPORT
int SCIPgetExprNEnfosNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the data for one of the enforcements of an expression */
SCIP_EXPORT
void SCIPgetExprEnfoDataNonlinear(
   SCIP_EXPR*            expr,               /**< expression */
   int                   idx,                /**< position of enforcement in enfos array */
   SCIP_NLHDLR**         nlhdlr,             /**< buffer to store nlhldr */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< buffer to store nlhdlr data for expression, or NULL */
   SCIP_NLHDLR_METHOD* nlhdlrparticipation, /**< buffer to store methods where nonlinear handler participates, or NULL */
   SCIP_Bool*            sepabelowusesactivity, /**< buffer to store whether sepabelow uses activity of some expression, or NULL */
   SCIP_Bool*            sepaaboveusesactivity, /**< buffer to store whether sepaabove uses activity of some expression, or NULL */
   SCIP_Real*            auxvalue            /**< buffer to store current auxvalue, or NULL */
   );

/** sets the auxiliary value of expression for one of the enforcements of an expression */
SCIP_EXPORT
void SCIPsetExprEnfoAuxValueNonlinear(
   SCIP_EXPR*            expr,               /**< expression */
   int                   idx,                /**< position of enforcement in enfos array */
   SCIP_Real             auxvalue            /**< the new value of auxval */
   );

/** number of nonlinear handlers whose activity computation and propagation methods depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetExprNPropUsesActivityNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   );

/** number of nonlinear handlers whose separation methods (estimate or enforcement) depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetExprNSepaUsesActivityNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   );

/** number of nonlinear handlers whose separation methods (estimate or enforcement) use auxiliary variable of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
SCIP_EXPORT
unsigned int SCIPgetExprNAuxvarUsesNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   );

/** method to be called by a nlhdlr during NLHDLRDETECT to notify an expression that it will be used
 *
 * - if `useauxvar` is enabled, then ensures that an auxiliary variable will be created in INITLP
 * - if `useactivityforprop` or `useactivityforsepa{below,above}` is enabled, then ensured that activity will be updated for `expr`
 * - if `useactivityforprop` is enabled, then increments the count returned by SCIPgetExprNPropUsesActivityNonlinear()
 * - if `useactivityforsepa{below,above}` is enabled, then increments the count returned by SCIPgetExprNSepaUsesActivityNonlinear()
 *   and also increments this count for all variables in the expression.
 *
 * The distinction into `useactivityforprop` and `useactivityforsepa{below,above}` is to recognize variables which domain influences
 * under/overestimators. Domain propagation routines (like OBBT) may invest more work for these variables.
 * The distinction into `useactivityforsepabelow` and `useactivityforsepaabove` is to recognize whether a nlhdlr that called this method
 * will use activity of `expr` in enfomethod \ref SCIP_NLHDLR_METHOD_SEPABELOW or \ref SCIP_NLHDLR_METHOD_SEPAABOVE.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPregisterExprUsageNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool             useauxvar,          /**< whether an auxiliary variable will be used for estimate or cut generation */
   SCIP_Bool             useactivityforprop, /**< whether activity of expr will be used by domain propagation or activity calculation (inteval) */
   SCIP_Bool             useactivityforsepabelow, /**< whether activity of expr will be used by underestimation */
   SCIP_Bool             useactivityforsepaabove  /**< whether activity of expr will be used by overestimation */
   );

/** computes absolute violation for auxvar relation in an expression w.r.t. original variables
 *
 * Assume the expression is f(x), where x are original (i.e., not auxiliary) variables.
 * Assume that f(x) is associated with auxiliary variable z.
 *
 * If there are negative locks, then returns the violation of z &le; f(x) and sets `violover` to TRUE.
 * If there are positive locks, then returns the violation of z &ge; f(x) and sets `violunder` to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z = f(x).
 *
 * If necessary, f is evaluated in the given solution. If that fails (domain error),
 * then `viol` is set to SCIPinfinity() and both `violover` and `violunder` are set to TRUE.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetExprAbsOrigViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Longint          soltag,             /**< tag of solution */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(x) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(x) is violated, or NULL */
   );

/** computes absolute violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * If there are negative locks, then returns the violation of z &le; f(w) and sets `violover` to TRUE.
 * If there are positive locks, then returns the violation of z &ge; f(w) and sets `violunder` to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z = f(w).
 *
 * If the given value of f(w) is SCIP_INVALID, then `viol` is set to SCIPinfinity() and
 * both `violover` and `violunder` are set to TRUE.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetExprAbsAuxViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   );

/** computes relative violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * Taking the absolute violation from SCIPgetExprAbsAuxViolationNonlinear(), this function returns
 * the absolute violation divided by max(1,|f(w)|).
 *
 * If the given value of f(w) is SCIP_INVALID, then `viol` is set to SCIPinfinity() and
 * both `violover` and `violunder` are set to TRUE.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetExprRelAuxViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   );

/** returns bounds on the expression
 *
 * This gives an intersection of bounds from
 * - activity calculation (SCIPexprGetActivity()), if valid,
 * - auxiliary variable, if present,
 * - stored by SCIPtightenExprIntervalNonlinear() during domain propagation
 *
 * @note The returned interval can be empty!
 */
SCIP_EXPORT
SCIP_INTERVAL SCIPgetExprBoundsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** informs the expression about new bounds that can be used for reverse-propagation and to tighten bounds of
 * corresponding (auxiliary) variable (if any)
 *
 * @attention this function should only be called during domain propagation in cons_nonlinear
 */
SCIP_EXPORT
SCIP_RETCODE SCIPtightenExprIntervalNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be tightened */
   SCIP_INTERVAL         newbounds,          /**< new bounds for the expression */
   SCIP_Bool*            cutoff,             /**< buffer to store whether a cutoff was detected */
   int*                  ntightenings        /**< buffer to add the total number of tightenings, or NULL */
   );

/** mark constraints that include this expression to be propagated again
 *
 * This can be used by, e.g., nlhdlrs, to trigger a new propagation of constraints without
 * a change of variable bounds, e.g., because new information on the expression is available
 * that could potentially lead to tighter expression activity values.
 *
 * Note, that this call marks also constraints for propagation which only share some variable
 * with this expression.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPmarkExprPropagateNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression to propagate again */
   );

/** adds violation-branching score to an expression
 *
 * Adds a score to the expression-specific violation-branching score, thereby marking it as branching candidate.
 * The expression must either be a variable expression or have an aux-variable.
 * In the latter case, branching on auxiliary variables must have been enabled.
 * In case of doubt, use SCIPaddExprsViolScoreNonlinear(). Roughly, the difference between these functions is that the current
 * function adds `violscore` to the expression directly, while SCIPaddExprsViolScoreNonlinear() will split the
 * violation score among all the given expressions according to parameter constraints/nonlinear/branching/violsplit.
 *
 * @see SCIPaddExprsViolScoreNonlinear()
 */
SCIP_EXPORT
void SCIPaddExprViolScoreNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression where to add branching score */
   SCIP_Real             violscore           /**< violation score to add to expression */
   );

/** adds violation-branching score to a set of expressions, distributing the score among all the expressions
 *
 * Each expression must either be a variable expression or have an aux-variable.
 * If branching on aux-variables is disabled, then the violation branching score will be distributed among all
 * variables present in `exprs`.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddExprsViolScoreNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           exprs,              /**< expressions where to add branching score */
   int                   nexprs,             /**< number of expressions */
   SCIP_Real             violscore,          /**< violation score to add to expression */
   SCIP_SOL*             sol,                /**< current solution */
   SCIP_Bool*            success             /**< buffer to store whether at least one violscore was added */
   );

/** gives violation-branching score stored in expression, or 0.0 if no valid score has been stored */
SCIP_EXPORT
SCIP_Real SCIPgetExprViolScoreNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   );

/** returns the partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error)
 *
 * @see SCIPexprGetDerivative()
 */
SCIP_EXPORT
SCIP_Real SCIPgetExprPartialDiffNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< root expression of constraint used in the last SCIPevalExprGradient() call */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   );

/** returns the var's coordinate of Hu partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error)
 *
 * @see SCIPexprGetBardot()
 */
SCIP_EXPORT
SCIP_Real SCIPgetExprPartialDiffGradientDirNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< root expression of constraint used in the last SCIPevalExprHessianDir() call */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   );

/** evaluates quadratic term in a solution w.r.t. auxiliary variables
 *
 * \note This requires that for every expr used in the quadratic data, a variable or auxiliary variable is available.
 */
SCIP_EXPORT
SCIP_Real SCIPevalExprQuadraticAuxNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL for LP solution */
   );

/** @} */
/** @} */

/**@addtogroup PublicNlhdlrInterfaceMethods
 * @{
 */

/** creates a nonlinear handler and includes it into the nonlinear constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlr,             /**< buffer where to store nonlinear handler */
   const char*           name,               /**< name of nonlinear handler (must not be NULL) */
   const char*           desc,               /**< description of nonlinear handler (can be NULL) */
   int                   detectpriority,     /**< detection priority of nonlinear handler */
   int                   enfopriority,       /**< enforcement priority of nonlinear handler */
   SCIP_DECL_NLHDLRDETECT((*detect)),        /**< structure detection callback of nonlinear handler */
   SCIP_DECL_NLHDLREVALAUX((*evalaux)),      /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_NLHDLRDATA*      nlhdlrdata          /**< data of nonlinear handler (can be NULL) */
   );

/** get number of nonlinear handler */
SCIP_EXPORT
int SCIPgetNNlhdlrsNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   );

/** get nonlinear handlers */
SCIP_EXPORT
SCIP_NLHDLR** SCIPgetNlhdlrsNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   );
/** returns a nonlinear handler of a given name (or NULL if not found) */
SCIP_EXPORT
SCIP_NLHDLR* SCIPfindNlhdlrNonlinear(
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   const char*           name                /**< name of nonlinear handler */
   );

/** gives expression data that a given nonlinear handler stored in an expression
 *
 * Returns NULL if expr has not been detected by nlhdlr or nlhdlr did not store data.
 */
SCIP_EXPORT
SCIP_NLHDLREXPRDATA* SCIPgetNlhdlrExprDataNonlinear(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_EXPR*            expr                /**< expression */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
