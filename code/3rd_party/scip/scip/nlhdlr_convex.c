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

/**@file   nlhdlr_convex.c
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  nonlinear handlers for convex and concave expressions
 * @author Benjamin Mueller
 * @author Stefan Vigerske
 *
 * TODO convex: perturb reference point if separation fails due to too large numbers
 */

#include <string.h>

#include "scip/nlhdlr_convex.h"
#include "scip/pub_nlhdlr.h"
#include "scip/scip_expr.h"
#include "scip/cons_nonlinear.h"
#include "scip/expr_var.h"
#include "scip/pub_misc_rowprep.h"
#include "scip/dbldblarith.h"

/* fundamental nonlinear handler properties */
#define CONVEX_NLHDLR_NAME             "convex"
#define CONVEX_NLHDLR_DESC             "handler that identifies and estimates convex expressions"
#define CONVEX_NLHDLR_DETECTPRIORITY   50
#define CONVEX_NLHDLR_ENFOPRIORITY     50

#define CONCAVE_NLHDLR_NAME            "concave"
#define CONCAVE_NLHDLR_DESC            "handler that identifies and estimates concave expressions"
#define CONCAVE_NLHDLR_DETECTPRIORITY  40
#define CONCAVE_NLHDLR_ENFOPRIORITY    40

#define DEFAULT_DETECTSUM             FALSE
#define DEFAULT_EXTENDEDFORM          TRUE
#define DEFAULT_CVXQUADRATIC_CONVEX   TRUE
#define DEFAULT_CVXQUADRATIC_CONCAVE  FALSE
#define DEFAULT_CVXSIGNOMIAL          TRUE
#define DEFAULT_CVXPRODCOMP           TRUE
#define DEFAULT_HANDLETRIVIAL         FALSE

#define INITLPMAXVARVAL          1000.0 /**< maximal absolute value of variable for still generating a linearization cut at that point in initlp */

/*lint -e440*/
/*lint -e441*/
/*lint -e666*/
/*lint -e777*/

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_NlhdlrExprData
{
   SCIP_EXPR*            nlexpr;             /**< expression (copy) for which this nlhdlr estimates */
   SCIP_HASHMAP*         nlexpr2origexpr;    /**< mapping of our copied expression to original expression */

   int                   nleafs;             /**< number of distinct leafs of nlexpr, i.e., number of distinct (auxiliary) variables handled */
   SCIP_EXPR**           leafexprs;          /**< distinct leaf expressions (excluding value-expressions), thus variables */
};

/** nonlinear handler data */
struct SCIP_NlhdlrData
{
   SCIP_Bool             isnlhdlrconvex;     /**< whether this data is used for the convex nlhdlr (TRUE) or the concave one (FALSE) */
   SCIP_SOL*             evalsol;            /**< solution used for evaluating expression in a different point,
                                                  e.g., for facet computation of vertex-polyhedral function */

   /* parameters */
   SCIP_Bool             detectsum;          /**< whether to run detection when the root of an expression is a non-quadratic sum */
   SCIP_Bool             extendedform;       /**< whether to create extended formulations instead of looking for maximal possible subexpression */

   /* advanced parameters (maybe remove some day) */
   SCIP_Bool             cvxquadratic;       /**< whether to use convexity check on quadratics */
   SCIP_Bool             cvxsignomial;       /**< whether to use convexity check on signomials */
   SCIP_Bool             cvxprodcomp;        /**< whether to use convexity check on product composition f(h)*h */
   SCIP_Bool             handletrivial;      /**< whether to handle trivial expressions, i.e., those where all children are variables */
};

/** data struct to be be passed on to vertexpoly-evalfunction (see SCIPcomputeFacetVertexPolyhedralNonlinear) */
typedef struct
{
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata;
   SCIP_SOL*             evalsol;
   SCIP*                 scip;
} VERTEXPOLYFUN_EVALDATA;

/** stack used in constructExpr to store expressions that need to be investigated ("to do list") */
typedef struct
{
   SCIP_EXPR**           stack;              /**< stack elements */
   int                   stacksize;          /**< allocated space (in number of pointers) */
   int                   stackpos;           /**< position of top element of stack */
} EXPRSTACK;

#define DECL_CURVCHECK(x) SCIP_RETCODE x( \
   SCIP*                 scip,               /**< SCIP data structure */ \
   SCIP_EXPR*            nlexpr,             /**< nlhdlr-expr to check */ \
   SCIP_Bool             isrootexpr,         /**< whether nlexpr is the root from where detection has been started */ \
   EXPRSTACK*            stack,              /**< stack where to add generated leafs */ \
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from our expression copy to original expression */ \
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< data of nlhdlr */ \
   SCIP_HASHMAP*         assumevarfixed,     /**< hashmap containing variables that should be assumed to be fixed, or NULL */ \
   SCIP_Bool*            success             /**< whether we found something */ \
   )

/*
 * static methods
 */

/** create nlhdlr-expression
 *
 * does not create children, i.e., assumes that this will be a leaf
 */
static
SCIP_RETCODE nlhdlrExprCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from copied to original expression */
   SCIP_EXPR**           nlhdlrexpr,         /**< buffer to store created expr */
   SCIP_EXPR*            origexpr,           /**< original expression to be copied */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
   )
{
   assert(scip != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(nlhdlrexpr != NULL);
   assert(origexpr != NULL);

   if( SCIPexprGetNChildren(origexpr) == 0 )
   {
      /* for leaves, do not copy */
      *nlhdlrexpr = origexpr;
      SCIPcaptureExpr(*nlhdlrexpr);
      if( !SCIPhashmapExists(nlexpr2origexpr, (void*)*nlhdlrexpr) )
      {
         SCIP_CALL( SCIPhashmapInsert(nlexpr2origexpr, (void*)*nlhdlrexpr, (void*)origexpr) );
      }
      return SCIP_OKAY;
   }

   /* create copy of expression, but without children */
   SCIP_CALL( SCIPduplicateExprShallow(scip, origexpr, nlhdlrexpr, NULL, NULL) );
   assert(*nlhdlrexpr != NULL);  /* copies within the same SCIP must always work */

   /* store the curvature we want to get in the curvature flag of the copied expression
    * it's a bit of a misuse, but once we are done with everything, this is actually correct
    */
   SCIPexprSetCurvature(*nlhdlrexpr, curv);

   /* remember which the original expression was */
   SCIP_CALL( SCIPhashmapInsert(nlexpr2origexpr, (void*)*nlhdlrexpr, (void*)origexpr) );

   return SCIP_OKAY;
}

/** expand nlhdlr-expression by adding children according to original expression */
static
SCIP_RETCODE nlhdlrExprGrowChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from copied to original expression */
   SCIP_EXPR*            nlhdlrexpr,         /**< expression for which to create children */
   SCIP_EXPRCURV*        childrencurv        /**< curvature required for children, or NULL if to set to UNKNOWN */
   )
{
   SCIP_EXPR* origexpr;
   SCIP_EXPR* child;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(nlhdlrexpr != NULL);
   assert(SCIPexprGetNChildren(nlhdlrexpr) == 0);

   origexpr = (SCIP_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlhdlrexpr);

   nchildren = SCIPexprGetNChildren(origexpr);
   if( nchildren == 0 )
      return SCIP_OKAY;

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( nlhdlrExprCreate(scip, nlexpr2origexpr, &child, SCIPexprGetChildren(origexpr)[i],
         childrencurv != NULL ? childrencurv[i] : SCIP_EXPRCURV_UNKNOWN) );
      SCIP_CALL( SCIPappendExprChild(scip, nlhdlrexpr, child) );
      /* append captures child, so we can release the capture from nlhdlrExprCreate */
      SCIP_CALL( SCIPreleaseExpr(scip, &child) );
   }

   assert(SCIPexprGetNChildren(nlhdlrexpr) == SCIPexprGetNChildren(origexpr));

   return SCIP_OKAY;
}

/** evaluate expression at solution w.r.t. auxiliary variables */
static
SCIP_DECL_VERTEXPOLYFUN(nlhdlrExprEvalConcave)
{
   VERTEXPOLYFUN_EVALDATA* evaldata = (VERTEXPOLYFUN_EVALDATA*)funcdata;
   int i;

   assert(args != NULL);
   assert(nargs == evaldata->nlhdlrexprdata->nleafs);
   assert(evaldata != NULL);

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(evaldata->scip, "eval vertexpolyfun at\n");
#endif
   for( i = 0; i < nargs; ++i )
   {
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(evaldata->scip, "  <%s> = %g\n",
            SCIPvarGetName(SCIPgetVarExprVar(evaldata->nlhdlrexprdata->leafexprs[i])), args[i]);
#endif
      SCIP_CALL_ABORT( SCIPsetSolVal(evaldata->scip, evaldata->evalsol,
            SCIPgetVarExprVar(evaldata->nlhdlrexprdata->leafexprs[i]), args[i]) );
   }

   SCIP_CALL_ABORT( SCIPevalExpr(evaldata->scip, evaldata->nlhdlrexprdata->nlexpr, evaldata->evalsol, 0L) );

   return SCIPexprGetEvalValue(evaldata->nlhdlrexprdata->nlexpr);
}

/** initialize expression stack */
static
SCIP_RETCODE exprstackInit(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRSTACK*            exprstack,          /**< stack to initialize */
   int                   initsize            /**< initial size */
   )
{
   assert(scip != NULL);
   assert(exprstack != NULL);
   assert(initsize > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &exprstack->stack, initsize) );
   exprstack->stacksize = initsize;
   exprstack->stackpos = -1;

   return SCIP_OKAY;
}

/** free expression stack */
static
void exprstackFree(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRSTACK*            exprstack           /**< free expression stack */
   )
{
   assert(scip != NULL);
   assert(exprstack != NULL);

   SCIPfreeBufferArray(scip, &exprstack->stack);
}

/** add expressions to expression stack */
static
SCIP_RETCODE exprstackPush(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRSTACK*            exprstack,          /**< expression stack */
   int                   nexprs,             /**< number of expressions to push */
   SCIP_EXPR**           exprs               /**< expressions to push */
   )
{
   assert(scip != NULL);
   assert(exprstack != NULL);

   if( nexprs == 0 )
      return SCIP_OKAY;

   assert(exprs != NULL);

   if( exprstack->stackpos+1 + nexprs > exprstack->stacksize )  /*lint !e644*/
   {
      exprstack->stacksize = SCIPcalcMemGrowSize(scip, exprstack->stackpos+1 + nexprs);    /*lint !e644*/
      SCIP_CALL( SCIPreallocBufferArray(scip, &exprstack->stack, exprstack->stacksize) );
   }

   memcpy(exprstack->stack + (exprstack->stackpos+1), exprs, nexprs * sizeof(SCIP_EXPR*));  /*lint !e679*/ /*lint !e737*/
   exprstack->stackpos += nexprs;

   return SCIP_OKAY;
}

/** gives expression from top of expression stack and removes it from stack */
static
SCIP_EXPR* exprstackPop(
   EXPRSTACK*            exprstack           /**< expression stack */
   )
{
   assert(exprstack != NULL);
   assert(exprstack->stackpos >= 0);

   return exprstack->stack[exprstack->stackpos--];
}

/** indicate whether expression stack is empty */
static
SCIP_Bool exprstackIsEmpty(
   EXPRSTACK*            exprstack           /**< expression stack */
   )
{
   assert(exprstack != NULL);

   return exprstack->stackpos < 0;
}

/** looks whether given expression is (proper) quadratic and has a given curvature
 *
 * If having a given curvature, currently require all arguments of quadratic to be linear.
 * Hence, not using this for a simple square term, as curvCheckExprhdlr may provide a better condition on argument curvature then.
 * Also we wouldn't do anything useful for a single bilinear term.
 * Thus, run on sum's only.
 */
static
DECL_CURVCHECK(curvCheckQuadratic)
{  /*lint --e{715}*/
   SCIP_EXPR* expr;
   SCIP_EXPRCURV presentcurv;
   SCIP_EXPRCURV wantedcurv;
   SCIP_HASHSET* lonelysquares = NULL;
   SCIP_Bool isquadratic;
   int nbilinexprs;
   int nquadexprs;
   int i;

   assert(nlexpr != NULL);
   assert(stack != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( !nlhdlrdata->cvxquadratic )
      return SCIP_OKAY;

   if( !SCIPisExprSum(scip, nlexpr) )
      return SCIP_OKAY;

   wantedcurv = SCIPexprGetCurvature(nlexpr);
   if( wantedcurv == SCIP_EXPRCURV_LINEAR )
      return SCIP_OKAY;
   assert(wantedcurv == SCIP_EXPRCURV_CONVEX || wantedcurv == SCIP_EXPRCURV_CONCAVE);

   expr = (SCIP_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr);
   assert(expr != NULL);

   /* check whether quadratic */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, expr, &isquadratic) );

   /* if not quadratic, then give up here */
   if( !isquadratic )
      return SCIP_OKAY;

   SCIPexprGetQuadraticData(expr, NULL, NULL, NULL, NULL, &nquadexprs, &nbilinexprs, NULL, NULL);

   /* if only single square term (+linear), then give up here (let curvCheckExprhdlr handle this) */
   if( nquadexprs <= 1 )
      return SCIP_OKAY;

   /* if root expression is only sum of squares (+linear) and detectsum is disabled, then give up here, too */
   if( isrootexpr && !nlhdlrdata->detectsum && nbilinexprs == 0 )
      return SCIP_OKAY;

   /* get curvature of quadratic
    * TODO as we know what curvature we want, we could first do some simple checks like computing xQx for a random x
    */
   SCIP_CALL( SCIPcomputeExprQuadraticCurvature(scip, expr, &presentcurv, assumevarfixed, FALSE) );

   /* if not having desired curvature, return */
   if( presentcurv != wantedcurv )
      return SCIP_OKAY;

   *success = TRUE;

   if( !nlhdlrdata->detectsum )
   {
      /* first step towards block-decomposition of quadratic term:
       * collect all square-expressions (in original expr) which have no adjacent bilinear term
       * we will treat these x^2 as linear, i.e., add an auxvar for them, so x^2 maybe linearized
       * more efficiently (in particular if x is discrete)
       */
      SCIP_CALL( SCIPhashsetCreate(&lonelysquares, SCIPblkmem(scip), nquadexprs) );
      for( i = 0; i < nquadexprs; ++i )
      {
         int nadjbilin;
         SCIP_EXPR* sqrexpr;

         SCIPexprGetQuadraticQuadTerm(expr, i, NULL, NULL, NULL, &nadjbilin, NULL, &sqrexpr);
         if( nadjbilin == 0 )
         {
            assert(sqrexpr != NULL);
            SCIP_CALL( SCIPhashsetInsert(lonelysquares, SCIPblkmem(scip), (void*)sqrexpr) );
         }
      }
   }

   /* add immediate children to nlexpr */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, nlexpr, NULL) );
   assert(SCIPexprGetNChildren(nlexpr) == SCIPexprGetNChildren(expr));

   /* put children that are not square or product on stack
    * grow child for children that are square or product and put this child on stack
    * require all children to be linear
    */
   for( i = 0; i < SCIPexprGetNChildren(nlexpr); ++i )
   {
      SCIP_EXPR* child;
      SCIP_EXPRCURV curvlinear[2] = { SCIP_EXPRCURV_LINEAR, SCIP_EXPRCURV_LINEAR };

      child = SCIPexprGetChildren(nlexpr)[i];
      assert(child != NULL);

      assert(SCIPhashmapGetImage(nlexpr2origexpr, (void*)child) == SCIPexprGetChildren(expr)[i]);

      if( SCIPisExprPower(scip, child) && SCIPgetExponentExprPow(child) == 2.0 &&
         (lonelysquares == NULL || !SCIPhashsetExists(lonelysquares, SCIPexprGetChildren(expr)[i])) )
      {
         /* square term that isn't lonely, i.e., orig-version of child is a square-expr and nadjbilin>0 */
         SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, child, curvlinear) );
         assert(SCIPexprGetNChildren(child) == 1);
         SCIP_CALL( exprstackPush(scip, stack, 1, SCIPexprGetChildren(child)) );
      }
      else if( SCIPisExprProduct(scip, child) && SCIPexprGetNChildren(SCIPexprGetChildren(expr)[i]) == 2 )
         /* using original version of child here as NChildren(child)==0 atm */
      {
         /* bilinear term */
         SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, child, curvlinear) );
         assert(SCIPexprGetNChildren(child) == 2);
         SCIP_CALL( exprstackPush(scip, stack, 2, SCIPexprGetChildren(child)) );
      }
      else
      {
         /* linear term (or term to be considered as linear) or lonely square term
          * if we want extended formulations, then require linearity, so an auxvar will be introduced if it is nonlinear
          * if we do not want extended formulations, then the term needs to have curvature "wantedcurv"
          *   thus, if the coef is negative, then the child needs to have the curvature opposite to "wantedcurv"
          */
         if( nlhdlrdata->extendedform )
            SCIPexprSetCurvature(child, SCIP_EXPRCURV_LINEAR);
         else
            SCIPexprSetCurvature(child, SCIPexprcurvMultiply(SCIPgetCoefsExprSum(nlexpr)[i], wantedcurv));
         SCIP_CALL( exprstackPush(scip, stack, 1, &child) );
      }
   }

   if( lonelysquares != NULL )
      SCIPhashsetFree(&lonelysquares, SCIPblkmem(scip));

   return SCIP_OKAY;
}

/** looks whether top of given expression looks like a signomial that can have a given curvature
 *
 * e.g., sqrt(x)*sqrt(y) is convex if x,y >= 0 and x and y are convex
 *
 * unfortunately, doesn't work for tls, because i) it's originally sqrt(x*y), and ii) it is expanded into some sqrt(z*y+y);
 * but works for cvxnonsep_nsig
 */
static
DECL_CURVCHECK(curvCheckSignomial)
{  /*lint --e{715}*/
   SCIP_EXPR* expr;
   SCIP_EXPR* child;
   SCIP_Real* exponents;
   SCIP_INTERVAL* bounds;
   SCIP_EXPRCURV* curv;
   int nfactors;
   int i;

   assert(nlexpr != NULL);
   assert(stack != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( !nlhdlrdata->cvxsignomial )
      return SCIP_OKAY;

   if( !SCIPisExprProduct(scip, nlexpr) )
      return SCIP_OKAY;

   expr = (SCIP_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr);
   assert(expr != NULL);

   nfactors = SCIPexprGetNChildren(expr);
   if( nfactors <= 1 )  /* boooring */
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &exponents, nfactors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nfactors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &curv, nfactors) );

   for( i = 0; i < nfactors; ++i )
   {
      child = SCIPexprGetChildren(expr)[i];
      assert(child != NULL);

      if( !SCIPisExprPower(scip, child) )
      {
         exponents[i] = 1.0;
         SCIP_CALL( SCIPevalExprActivity(scip, child) );
         bounds[i] = SCIPexprGetActivity(child);
      }
      else
      {
         exponents[i] = SCIPgetExponentExprPow(child);
         SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(child)[0]) );
         bounds[i] = SCIPexprGetActivity(SCIPexprGetChildren(child)[0]);
      }
   }

   if( !SCIPexprcurvMonomialInv(SCIPexprcurvMultiply(SCIPgetCoefExprProduct(expr), SCIPexprGetCurvature(nlexpr)),
       nfactors, exponents, bounds, curv) )
      goto TERMINATE;

   /* add immediate children to nlexpr
    * some entries in curv actually apply to arguments of pow's, will correct this next
    */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, nlexpr, curv) );
   assert(SCIPexprGetNChildren(nlexpr) == nfactors);

   /* put children that are not power on stack
    * grow child for children that are power and put this child on stack
    * if extendedform, then require children to be linear
    * unless they are linear, an auxvar will be introduced for them and thus they will be handled as var here
    */
   for( i = 0; i < nfactors; ++i )
   {
      child = SCIPexprGetChildren(nlexpr)[i];
      assert(child != NULL);

      if( SCIPisExprPower(scip, child) )
      {
         SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, child, &curv[i]) );
         assert(SCIPexprGetNChildren(child) == 1);
         child = SCIPexprGetChildren(child)[0];
      }
      assert(SCIPexprGetNChildren(child) == 0);

      if( nlhdlrdata->extendedform )
      {
         SCIPexprSetCurvature(child, SCIP_EXPRCURV_LINEAR);
#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "Extendedform: Require linearity for ");
         SCIPprintExpr(scip, child, NULL);
         SCIPinfoMessage(scip, NULL, "\n");
#endif
      }

      SCIP_CALL( exprstackPush(scip, stack, 1, &child) );
   }

   *success = TRUE;

TERMINATE:
   SCIPfreeBufferArray(scip, &curv);
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &exponents);

   return SCIP_OKAY;
}

/** looks for \f$f(c h(x)+d) h(x) \cdot \text{constant}\f$ and tries to conclude conditions on curvature
 *
 * Assume \f$h\f$ is univariate:
 * - First derivative is \f$f'(c h + d) c h' h + f(c h + d) h'\f$.
 * - Second derivative is \f{align}{&f''(c h + d) c h' c h' h + f'(c h + d) (c h'' h + c h' h') + f'(c h + d) c h' h' + f(c h + d) h'' \\
 *   =& f''(c h + d) c^2 h'^2 h + f'(c h + d) c h'' h + 2 f'(c h + d) c h'^2 + f(c h + d) h''.\f}
 *   Remove always positive factors leaves \f[f''(c h + d) h,\quad f'(c h + d) c h'' h,\quad f'(c h + d) c,\quad f(c h + d) h''.\f]
 *   For convexity we want all these terms to be nonnegative. For concavity we want all of them to be nonpositive.
 *   Note, that in each term either both \f$f'(c h + d)\f$ and \f$c\f$ occur, or none of them.
 * - Thus, \f$f(c h(x) + d)h(x)\f$ is convex if \f$cf\f$ is monotonically increasing \f$(c f' \geq 0)\f$ and either
 *   - \f$f\f$ is convex \f$(f'' \geq 0)\f$ and \f$h\f$ is nonnegative \f$(h \geq 0)\f$ and \f$h\f$ is convex \f$(h'' \geq 0)\f$ and [\f$f\f$ is nonnegative \f$(f \geq 0)\f$ or \f$h\f$ is linear \f$(h''=0)\f$], or
 *   - \f$f\f$ is concave \f$(f'' \leq 0)\f$ and \f$h\f$ is nonpositive \f$(h \leq 0)\f$ and \f$h\f$ is concave \f$(h'' \leq 0)\f$ and [\f$f\f$ is nonpositive \f$(f \leq 0)\f$ or \f$h\f$ is linear \f$(h''=0)\f$].
 * - Further, \f$f(c h(x) + d)h(x)\f$ is concave if \f$cf\f$ is monotonically decreasing \f$(c f' \leq 0)\f$ and either
 *   - f is convex \f$(f'' \geq 0)\f$ and \f$h\f$ is nonpositive \f$(h \leq 0)\f$ and \f$h\f$ is concave \f$(h'' \leq 0)\f$ and [\f$f\f$ is nonnegative \f$(f \geq 0)\f$ or \f$h\f$ is linear \f$(h''=0)\f$], or
 *   - f is concave \f$(f'' \leq 0)\f$ and \f$h\f$ is nonnegative \f$(h >= 0)\f$ and \f$h\f$ is convex \f$(h'' \geq 0)\f$ and [\f$f\f$ is nonpositive \f$(f \leq 0)\f$ or \f$h\f$ is linear \f$(h''=0)\f$].
 *
 * This should hold also for multivariate and linear \f$h\f$, as things are invariant under linear transformations.
 * Similar to signomial, I'll assume that this will also hold for other multivariate \f$h\f$ (someone has a formal proof?).
 */
static
DECL_CURVCHECK(curvCheckProductComposite)
{  /*lint --e{715}*/
   SCIP_EXPR* expr;
   SCIP_EXPR* f;
   SCIP_EXPR* h = NULL;
   SCIP_Real c = 0.0;
   SCIP_EXPR* ch = NULL; /* c * h */
   SCIP_INTERVAL fbounds;
   SCIP_INTERVAL hbounds;
   SCIP_MONOTONE fmonotonicity;
   SCIP_EXPRCURV desiredcurv;
   SCIP_EXPRCURV hcurv;
   SCIP_EXPRCURV dummy;
   int fidx;

   assert(nlexpr != NULL);
   assert(stack != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( !nlhdlrdata->cvxprodcomp )
      return SCIP_OKAY;

   if( !SCIPisExprProduct(scip, nlexpr) )
      return SCIP_OKAY;

   expr = (SCIP_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr);
   assert(expr != NULL);

   if( SCIPexprGetNChildren(expr) != 2 )
      return SCIP_OKAY;

   /* check whether we have f(c * h(x)) * h(x) or h(x) * f(c * h(x)) */
   for( fidx = 0; fidx <= 1; ++fidx )
   {
      f = SCIPexprGetChildren(expr)[fidx];

      if( SCIPexprGetNChildren(f) != 1 )
         continue;

      ch = SCIPexprGetChildren(f)[0];
      c = 1.0;
      h = ch;

      /* check whether ch is of the form c*h(x), then switch h to child ch */
      if( SCIPisExprSum(scip, ch) && SCIPexprGetNChildren(ch) == 1 )
      {
         c = SCIPgetCoefsExprSum(ch)[0];
         h = SCIPexprGetChildren(ch)[0];
         assert(c != 1.0 || SCIPgetConstantExprSum(ch) != 0.0);  /* we could handle this, but it should have been simplified away */
      }

#ifndef NLHDLR_CONVEX_UNITTEST
      /* can assume that duplicate subexpressions have been identified and comparing pointer is sufficient */
      if( SCIPexprGetChildren(expr)[1-fidx] == h )
#else
      /* called from unittest -> duplicate subexpressions were not identified -> compare more expensively */
      if( SCIPcompareExpr(scip, SCIPexprGetChildren(expr)[1-fidx], h) == 0 )
#endif
         break;
   }
   if( fidx == 2 )
      return SCIP_OKAY;

#ifdef SCIP_MORE_DEBUG
   SCIPinfoMessage(scip, NULL, "f(c*h+d)*h with f = %s, c = %g, d = %g, h = ", SCIPexprhdlrGetName(SCIPexprGetHdlr(f)),
         c, h != ch ? SCIPgetConstantExprSum(ch) : 0.0);
   SCIPprintExpr(scip, h, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   assert(c != 0.0);

   SCIP_CALL( SCIPevalExprActivity(scip, f) );
   SCIP_CALL( SCIPevalExprActivity(scip, h) );
   fbounds = SCIPexprGetActivity(f);
   hbounds = SCIPexprGetActivity(h);

   /* if h has mixed sign, then cannot conclude anything */
   if( hbounds.inf < 0.0 && hbounds.sup > 0.0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcallExprMonotonicity(scip, f, 0, &fmonotonicity) );

   /* if f is not monotone, then cannot conclude anything */
   if( fmonotonicity == SCIP_MONOTONE_UNKNOWN )
      return SCIP_OKAY;

   /* curvature we want to achieve (negate if product has negative coef) */
   desiredcurv = SCIPexprcurvMultiply(SCIPgetCoefExprProduct(nlexpr), SCIPexprGetCurvature(nlexpr));

   /* now check the conditions as stated above */
   if( desiredcurv == SCIP_EXPRCURV_CONVEX )
   {
      /* f(c h(x)+d)h(x) is convex if c*f is monotonically increasing (c f' >= 0) and either
      *   - f is convex (f'' >= 0) and h is nonnegative (h >= 0) and h is convex (h'' >= 0) and [f is nonnegative (f >= 0) or h is linear (h''=0)], or
      *   - f is concave (f'' <= 0) and h is nonpositive (h <= 0) and h is concave (h'' <= 0) and [f is nonpositive (f <= 0) or h is linear (h''=0)]
      *  as the curvature requirements on f are on f only and not the composition f(h), we can ignore the requirements returned by SCIPcallExprCurvature (last arg)
      */
      if( (c > 0.0 && fmonotonicity != SCIP_MONOTONE_INC) || (c < 0.0 && fmonotonicity != SCIP_MONOTONE_DEC) )
         return SCIP_OKAY;

      /* check whether f can be convex (h>=0) or concave (h<=0), resp., and derive requirements for h */
      if( hbounds.inf >= 0 )
      {
         SCIP_CALL( SCIPcallExprCurvature(scip, f, SCIP_EXPRCURV_CONVEX, success, &dummy) );

         /* now h also needs to be convex; and if f < 0, then h actually needs to be linear */
         if( fbounds.inf < 0.0 )
            hcurv = SCIP_EXPRCURV_LINEAR;
         else
            hcurv = SCIP_EXPRCURV_CONVEX;
      }
      else
      {
         SCIP_CALL( SCIPcallExprCurvature(scip, f, SCIP_EXPRCURV_CONCAVE, success, &dummy) );

         /* now h also needs to be concave; and if f > 0, then h actually needs to be linear */
         if( fbounds.sup > 0.0 )
            hcurv = SCIP_EXPRCURV_LINEAR;
         else
            hcurv = SCIP_EXPRCURV_CONCAVE;
      }
   }
   else
   {
      /* f(c h(x)+d)*h(x) is concave if c*f is monotonically decreasing (c f' <= 0) and either
      *   - f is convex (f'' >= 0) and h is nonpositive (h <= 0) and h is concave (h'' <= 0) and [f is nonnegative (f >= 0) or h is linear (h''=0)], or
      *   - f is concave (f'' <= 0) and h is nonnegative (h >= 0) and h is convex (h'' >= 0) and [f is nonpositive (f <= 0) or h is linear (h''=0)]
      *  as the curvature requirements on f are on f only and not the composition f(h), we can ignore the requirements returned by SCIPcallExprCurvature (last arg)
      */
      if( (c > 0.0 && fmonotonicity != SCIP_MONOTONE_DEC) || (c < 0.0 && fmonotonicity != SCIP_MONOTONE_INC) )
         return SCIP_OKAY;

      /* check whether f can be convex (h<=0) or concave (h>=0), resp., and derive requirements for h */
      if( hbounds.sup <= 0 )
      {
         SCIP_CALL( SCIPcallExprCurvature(scip, f, SCIP_EXPRCURV_CONVEX, success, &dummy) );

         /* now h also needs to be concave; and if f < 0, then h actually needs to be linear */
         if( fbounds.inf < 0.0 )
            hcurv = SCIP_EXPRCURV_LINEAR;
         else
            hcurv = SCIP_EXPRCURV_CONCAVE;
      }
      else
      {
         SCIP_CALL( SCIPcallExprCurvature(scip, f, SCIP_EXPRCURV_CONCAVE, success, &dummy) );

         /* now h also needs to be convex; and if f > 0, then h actually needs to be linear */
         if( fbounds.sup > 0.0 )
            hcurv = SCIP_EXPRCURV_LINEAR;
         else
            hcurv = SCIP_EXPRCURV_CONVEX;
      }
   }

   if( !*success )
      return SCIP_OKAY;

   /* add immediate children (f and ch) to nlexpr; we set required curvature for h further below */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, nlexpr, NULL) );
   assert(SCIPexprGetNChildren(nlexpr) == 2);

   /* copy of f (and h) should have same child position in nlexpr as f (and h) has on expr (resp) */
   assert(SCIPhashmapGetImage(nlexpr2origexpr, (void*)SCIPexprGetChildren(nlexpr)[fidx]) == (void*)f);
#ifndef NLHDLR_CONVEX_UNITTEST
   assert(SCIPhashmapGetImage(nlexpr2origexpr, (void*)SCIPexprGetChildren(nlexpr)[1-fidx]) == (void*)h);
#endif
   /* push this h onto stack for further checking */
   SCIP_CALL( exprstackPush(scip, stack, 1, &(SCIPexprGetChildren(nlexpr)[1-fidx])) );

   /* if we prefer extended formulations, then we always want h() to be linear */
   if( nlhdlrdata->extendedform )
      hcurv = SCIP_EXPRCURV_LINEAR;

   /* h-child of product should have curvature hcurv */
   SCIPexprSetCurvature(SCIPexprGetChildren(nlexpr)[1-fidx], hcurv);

   if( h != ch )
   {
      /* add copy of ch as child to copy of f */
      SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, SCIPexprGetChildren(nlexpr)[fidx], NULL) );
      assert(SCIPexprGetNChildren(SCIPexprGetChildren(nlexpr)[fidx]) == 1);
      assert(SCIPhashmapGetImage(nlexpr2origexpr, (void*)SCIPexprGetChildren(SCIPexprGetChildren(nlexpr)[fidx])[0]) == (void*)ch);

      /* add copy of h (created above as child of product) as child in copy of ch */
      SCIP_CALL( SCIPappendExprChild(scip,
         SCIPexprGetChildren(SCIPexprGetChildren(nlexpr)[fidx])[0] /* copy of ch */,
         SCIPexprGetChildren(nlexpr)[1-fidx] /* copy of h */) );
   }
   else
   {
      /* add copy of h (created above as child of product) as child in copy of f */
      SCIP_CALL( SCIPappendExprChild(scip,
         SCIPexprGetChildren(nlexpr)[fidx] /* copy of f */,
         SCIPexprGetChildren(nlexpr)[1-fidx] /* copy of h */) );
   }

   return SCIP_OKAY;
}

/** use expression handlers curvature callback to check whether given curvature can be achieved */
static
DECL_CURVCHECK(curvCheckExprhdlr)
{  /*lint --e{715}*/
   SCIP_EXPR* origexpr;
   int nchildren;
   SCIP_EXPRCURV* childcurv;

   assert(nlexpr != NULL);
   assert(stack != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(success != NULL);

   origexpr = (SCIP_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, nlexpr);
   assert(origexpr != NULL);
   nchildren = SCIPexprGetNChildren(origexpr);

   if( nchildren == 0 )
   {
      /* if originally no children, then should be var or value, which should have every curvature,
       * so should always be success
       */
      SCIP_CALL( SCIPcallExprCurvature(scip, origexpr, SCIPexprGetCurvature(nlexpr), success, NULL) );
      assert(*success);

      return SCIP_OKAY;
   }

   /* ignore sums if > 1 children
    * NOTE: this means that for something like 1+f(x), even if f is a trivial convex expression, we would handle 1+f(x)
    * with this nlhdlr, instead of formulating this as 1+z and handling z=f(x) with the default nlhdlr, i.e., the exprhdlr
    * today, I prefer handling this here, as it avoids introducing an extra auxiliary variable
    */
   if( isrootexpr && !nlhdlrdata->detectsum && SCIPisExprSum(scip, nlexpr) && nchildren > 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, nchildren) );

   /* check whether and under which conditions origexpr can have desired curvature */
   SCIP_CALL( SCIPcallExprCurvature(scip, origexpr, SCIPexprGetCurvature(nlexpr), success, childcurv) );
#ifdef SCIP_MORE_DEBUG
   SCIPprintExpr(scip, origexpr, NULL);
   SCIPinfoMessage(scip, NULL, " is %s? %d\n", SCIPexprcurvGetName(SCIPexprGetCurvature(nlexpr)), *success);
#endif
   if( !*success )
      goto TERMINATE;

   /* if origexpr can have curvature curv, then don't treat it as leaf, but include its children */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, nlexpr, childcurv) );
   assert(SCIPexprGetChildren(nlexpr) != NULL);
   assert(SCIPexprGetNChildren(nlexpr) == nchildren);

   /* If we prefer extended formulations, then require all children to be linear.
    * Unless they are, auxvars will be introduced and they will be handles as variables, which can be an
    * advantage in the context of extended formulations.
    */
   if( nlhdlrdata->extendedform )
   {
      int i;
      for( i = 0; i < nchildren; ++i )
         SCIPexprSetCurvature(SCIPexprGetChildren(nlexpr)[i], SCIP_EXPRCURV_LINEAR);
#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "require linearity for children of ");
      SCIPprintExpr(scip, origexpr, NULL);
      SCIPinfoMessage(scip, NULL, "\n");
#endif
   }

   /* add children expressions to to-do list (stack) */
   SCIP_CALL( exprstackPush(scip, stack, nchildren, SCIPexprGetChildren(nlexpr)) );

TERMINATE:
   SCIPfreeBufferArray(scip, &childcurv);

   return SCIP_OKAY;
}

/** curvature check and expression-growing methods
 *
 * some day this could be plugins added by users at runtime, but for now we have a fixed list here
 * @note curvCheckExprhdlr should be last
 */
static DECL_CURVCHECK((*CURVCHECKS[])) = { curvCheckProductComposite, curvCheckSignomial, curvCheckQuadratic, curvCheckExprhdlr };
/** number of curvcheck methods */
static const int NCURVCHECKS = sizeof(CURVCHECKS) / sizeof(void*);

/** checks whether expression is a sum with more than one child and each child being a variable or going to be a variable if `expr` is a nlhdlr-specific copy
 *
 * Within constructExpr(), we can have an expression of any type which is a copy of an original expression,
 * but without children. At the end of constructExpr() (after the loop with the stack), these expressions
 * will remain as leafs and will eventually be turned into variables in collectLeafs(). Thus, we treat
 * every child that has no children as if it were a variable. Theoretically, there is still the possibility
 * that it could be a constant (value-expression), but simplify should have removed these.
 */
static
SCIP_Bool exprIsMultivarLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression to check */
   )
{
   int nchildren;
   int c;

   assert(expr != NULL);

   if( !SCIPisExprSum(scip, expr) )
      return FALSE;

   nchildren = SCIPexprGetNChildren(expr);
   if( nchildren <= 1 )
      return FALSE;

   for( c = 0; c < nchildren; ++c )
      /*if( !SCIPisExprVar(scip, SCIPexprGetChildren(expr)[c]) ) */
      if( SCIPexprGetNChildren(SCIPexprGetChildren(expr)[c]) > 0 )
         return FALSE;

   return TRUE;
}

/** constructs a subexpression (as nlhdlr-expression) of maximal size that has a given curvature
 *
 * If the curvature cannot be achieved for an expression in the original expression graph,
 * then this expression becomes a leaf in the nlhdlr-expression.
 *
 * Sets `*rootnlexpr` to NULL if failed.
 */
static
SCIP_RETCODE constructExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_EXPR**           rootnlexpr,         /**< buffer to store created expression */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from our expression copy to original expression */
   int*                  nleafs,             /**< number of leafs in constructed expression */
   SCIP_EXPR*            rootexpr,           /**< expression */
   SCIP_EXPRCURV         curv,               /**< curvature to achieve */
   SCIP_HASHMAP*         assumevarfixed,     /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   SCIP_Bool             assumecurvature,    /**< whether to assume that desired curvature is given (skips curvature checks) */
   SCIP_Bool*            curvsuccess         /**< pointer to store whether the curvature could be achieved
                                                  w.r.t. the original variables (might be NULL) */
   )
{
   SCIP_EXPR* nlexpr;
   EXPRSTACK stack; /* to do list: expressions where to check whether they can have the desired curvature when taking their children into account */
   int oldstackpos;
   SCIP_Bool isrootexpr = TRUE;

   assert(scip != NULL);
   assert(nlhdlrdata != NULL);
   assert(rootnlexpr != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(nleafs != NULL);
   assert(rootexpr != NULL);
   assert(curv == SCIP_EXPRCURV_CONVEX || curv == SCIP_EXPRCURV_CONCAVE);

   /* create root expression */
   SCIP_CALL( nlhdlrExprCreate(scip, nlexpr2origexpr, rootnlexpr, rootexpr, curv) );

   *nleafs = 0;
   if( curvsuccess != NULL )
      *curvsuccess = TRUE;

   SCIP_CALL( exprstackInit(scip, &stack, 20) );
   SCIP_CALL( exprstackPush(scip, &stack, 1, rootnlexpr) );
   while( !exprstackIsEmpty(&stack) )
   {
      /* take expression from stack */
      nlexpr = exprstackPop(&stack);
      assert(nlexpr != NULL);
      assert(SCIPexprGetNChildren(nlexpr) == 0);

      oldstackpos = stack.stackpos;
      if( nlhdlrdata->isnlhdlrconvex && !SCIPexprhdlrHasBwdiff(SCIPexprGetHdlr(nlexpr)) )
      {
         /* if bwdiff is not implemented, then we could not generate cuts in the convex nlhdlr, so "stop" (treat nlexpr as variable) */
      }
      else if( !nlhdlrdata->isnlhdlrconvex && exprIsMultivarLinear(scip, (SCIP_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr)) )
      {
         /* if we are in the concave handler, we would like to treat linear multivariate subexpressions by a new auxvar always,
          * e.g., handle log(x+y) as log(z), z=x+y, because the estimation problem will be smaller then without making the estimator worse
          * (cons_nonlinear does this, too)
          * this check takes care of this when x and y are original variables
          * however, it isn't unlikely that we will have sums that become linear after we add auxvars for some children
          * this will be handled in a postprocessing below
          * for now, the check is performed on the original expression since there is not enough information in nlexpr yet
          */
#ifdef SCIP_MORE_DEBUG
         SCIPprintExpr(scip, SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr), NULL);
         SCIPinfoMessage(scip, NULL, "... is a multivariate linear sum that we'll treat as auxvar\n");
#endif
      }
      else if( SCIPexprGetCurvature(nlexpr) != SCIP_EXPRCURV_UNKNOWN && !assumecurvature )
      {
         /* if we are here, either convexity or concavity is required; try to check for this curvature */
         SCIP_Bool success;
         int method;

         /* try through curvature check methods until one succeeds */
         for( method = 0; method < NCURVCHECKS; ++method )
         {
            SCIP_CALL( CURVCHECKS[method](scip, nlexpr, isrootexpr, &stack, nlexpr2origexpr, nlhdlrdata, assumevarfixed, &success) );
            if( success )
               break;
         }
      }
      else
      {
         /* if we don't care about curvature in this subtree anymore (very unlikely),
          * or we are told to assume that the desired curvature is present (assumecurvature==TRUE),
          * then only continue iterating this subtree to assemble leaf expressions
          */
         SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr2origexpr, nlexpr, NULL) );

         /* add children expressions, if any, to to-do list (stack) */
         SCIP_CALL( exprstackPush(scip, &stack, SCIPexprGetNChildren(nlexpr), SCIPexprGetChildren(nlexpr)) );
      }
      assert(stack.stackpos >= oldstackpos);  /* none of the methods above should have removed something from the stack */

      isrootexpr = FALSE;

      /* if nothing was added, then none of the successors of nlexpr were added to the stack
       * this is either because nlexpr was already a variable or value expressions, thus a leaf,
       * or because the desired curvature could not be achieved, so it will be handled as variables, thus a leaf
       */
      if( stack.stackpos == oldstackpos )
      {
         ++*nleafs;

         /* check whether the new leaf is not an original variable (or constant) */
         if( curvsuccess != NULL && !SCIPisExprVar(scip, nlexpr) && !SCIPisExprValue(scip, nlexpr) )
            *curvsuccess = FALSE;
      }
   }

   exprstackFree(scip, &stack);

   if( !nlhdlrdata->isnlhdlrconvex && *rootnlexpr != NULL )
   {
      /* remove multivariate linear subexpressions, that is, change some f(z1+z2) into f(z3) (z3=z1+z2 will be done by nlhdlr_default)
       * this handles the case that was not covered by the above check, which could recognize f(x+y) for x, y original variables
       */
      SCIP_EXPRITER* it;

      SCIP_CALL( SCIPcreateExpriter(scip, &it) );
      SCIP_CALL( SCIPexpriterInit(it, *rootnlexpr, SCIP_EXPRITER_DFS, FALSE) );
      SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);

      while( !SCIPexpriterIsEnd(it) )
      {
         SCIP_EXPR* child;

         child = SCIPexpriterGetChildExprDFS(it);
         assert(child != NULL);

         /* We want to change some f(x+y+z) into just f(), where f is the expression the iterator points to
          * and x+y+z is child. A child of a child, e.g., z, may not be a variable yet (these are added in collectLeafs later),
          * but an expression of some nonlinear type without children.
          */
         if( exprIsMultivarLinear(scip, child) )
         {
            /* turn child (x+y+z) into a sum without children
             * collectLeafs() should then replace this by an auxvar
             */
#ifdef SCIP_MORE_DEBUG
            SCIPprintExpr(scip, child, NULL);
            SCIPinfoMessage(scip, NULL, "... is a multivariate linear sum that we'll treat as auxvar instead (postprocess)\n");
#endif

            /* TODO remove children from nlexpr2origexpr ?
             * should also do this if they are not used somewhere else; we could check nuses for this
             * however, it shouldn't matter to have some stray entries in the hashmap either
             */
            SCIP_CALL( SCIPremoveExprChildren(scip, child) );
            assert(SCIPexprGetNChildren(child) == 0);

            (void) SCIPexpriterSkipDFS(it);
         }
         else
         {
            (void) SCIPexpriterGetNext(it);
         }
      }

      SCIPfreeExpriter(&it);
   }

   if( *rootnlexpr != NULL )
   {
      SCIP_Bool istrivial = TRUE;

      /* if handletrivial is enabled, then only require that rootnlexpr itself has required curvature (so has children; see below) and
       * that we are not a trivial sum  (because the previous implementation of this nlhdlr didn't allow this, either)
       */
      if( !nlhdlrdata->handletrivial || SCIPisExprSum(scip, *rootnlexpr) )
      {
         /* if all children do not have children, i.e., are variables, or will be replaced by auxvars, then free
          * also if rootnlexpr has no children, then free
          */
         int i;
         for( i = 0; i < SCIPexprGetNChildren(*rootnlexpr); ++i )
         {
            if( SCIPexprGetNChildren(SCIPexprGetChildren(*rootnlexpr)[i]) > 0 )
            {
               istrivial = FALSE;
               break;
            }
         }
      }
      else if( SCIPexprGetNChildren(*rootnlexpr) > 0 )  /* if handletrivial, then just require children */
            istrivial = FALSE;

      if( istrivial )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, rootnlexpr) );
      }
   }

   return SCIP_OKAY;
}

/** collects (non-value) leaf expressions and ensure that they correspond to a variable (original or auxiliary)
 *
 * For children where we could not achieve the desired curvature, get the auxvar and replace the child by a
 * var-expression that points to this auxvar.
 * Collect all leaf expressions (if not a value-expression) and index them.
 */
static
SCIP_RETCODE collectLeafs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< nlhdlr expression data */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* nlexpr;
   SCIP_HASHMAP* leaf2index;
   int i;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlexpr != NULL);
   assert(nlhdlrexprdata->nlexpr2origexpr != NULL);
   /* nleafs should be the upper bound on the number of variables given by constructExpr
    * leafexprs should be NULL, as this is what we want to setup here
    */
   assert(nlhdlrexprdata->nleafs > 0);
   assert(nlhdlrexprdata->leafexprs == NULL);

   /* collect all auxvars and collect all variables */
   SCIP_CALL( SCIPhashmapCreate(&leaf2index, SCIPblkmem(scip), nlhdlrexprdata->nleafs) );
   nlhdlrexprdata->nleafs = 0;  /* we start a new count, this time skipping value-expressions */

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, nlhdlrexprdata->nlexpr, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);

   for( nlexpr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); nlexpr = SCIPexpriterGetNext(it) )
   {
      SCIP_EXPR* child;
      SCIP_EXPR* origexpr;

      assert(nlexpr != NULL);

      child = SCIPexpriterGetChildExprDFS(it);

      /* if the to-be-visited child has children, then it doesn't need to be replaced by a new expression (representing the auxvar) */
      if( SCIPexprGetNChildren(child) > 0 )
         continue;

      origexpr = (SCIP_EXPR*)SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)child);
      assert(origexpr != NULL);

      if( SCIPexprGetNChildren(origexpr) > 0 )
      {
         SCIP_EXPR* newchild;
         int childidx;
         SCIP_VAR* var;

         /* having a child that had children in original but not in copy means that we could not achieve the desired curvature
          * thus, replace by a new child that points to the auxvar of the original expression
          * we registered in createNlhdlrExprData that we need an auxvar, so it should exist now
          */
         var = SCIPgetExprAuxVarNonlinear(origexpr);
         assert(var != NULL);

         SCIP_CALL( SCIPcreateExprVar(scip, &newchild, var, NULL, NULL) );  /* this captures newchild once */

         childidx = SCIPexpriterGetChildIdxDFS(it);
         SCIP_CALL( SCIPreplaceExprChild(scip, nlexpr, childidx, newchild) );  /* this captures newchild again */

         /* do not remove child->origexpr from hashmap, as child may appear again due to common subexprs
          * (created by curvCheckProductComposite, for example)
          * if it doesn't reappear, though, but the memory address is reused, we need to make sure it
          * points to the right origexpr
          */
         /* SCIP_CALL( SCIPhashmapRemove(nlexpr2origexpr, (void*)child) ); */
         SCIP_CALL( SCIPhashmapSetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)newchild, (void*)origexpr) );

         if( !SCIPhashmapExists(leaf2index, (void*)newchild) )
         {
            /* new leaf -> new index and remember in hashmap */
            SCIP_CALL( SCIPhashmapInsertInt(leaf2index, (void*)newchild, nlhdlrexprdata->nleafs++) );
         }

         child = newchild;
         SCIP_CALL( SCIPreleaseExpr(scip, &newchild) );  /* because it was captured by both create and replace */
      }
      else if( SCIPisExprVar(scip, child) )
      {
         /* if variable, then add to hashmap, if not already there */
         if( !SCIPhashmapExists(leaf2index, (void*)child) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(leaf2index, (void*)child, nlhdlrexprdata->nleafs++) );
         }
      }
      /* else: it's probably a value-expression, nothing to do */

      /* update integrality flag for future leaf expressions: convex nlhdlr may use this information */
      SCIP_CALL( SCIPcomputeExprIntegrality(scip, child) );
   }
   assert(nlhdlrexprdata->nleafs > 0);

   SCIPfreeExpriter(&it);

   /* assemble auxvars array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(nlhdlrexprdata->leafexprs), nlhdlrexprdata->nleafs) );
   for( i = 0; i < SCIPhashmapGetNEntries(leaf2index); ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_EXPR* leaf;
      int idx;

      entry = SCIPhashmapGetEntry(leaf2index, i);
      if( entry == NULL )
         continue;

      leaf = (SCIP_EXPR*) SCIPhashmapEntryGetOrigin(entry);
      assert(leaf != NULL);
      assert(SCIPisExprVar(scip, leaf));

      idx = SCIPhashmapEntryGetImageInt(entry);
      assert(idx >= 0);
      assert(idx < nlhdlrexprdata->nleafs);

      nlhdlrexprdata->leafexprs[idx] = leaf;

      SCIPdebugMsg(scip, "leaf %d: <%s>\n", idx, SCIPvarGetName(SCIPgetVarExprVar(leaf)));
   }

   SCIPhashmapFree(&leaf2index);

   return SCIP_OKAY;
}

/** creates nonlinear handler expression data structure and registers expr usage */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nlhdlr data */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nlhdlr expression data */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR*            nlexpr,             /**< our copy of expression */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping of expression copy to original */
   int                   nleafs,             /**< number of leafs as counted by constructExpr */
   SCIP_NLHDLR_METHOD    participating       /**< the enfo methods in which we plan to participate */
   )
{
   SCIP_EXPRITER* it;
   SCIP_Bool usingaux;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata == NULL);
   assert(nlexpr != NULL);
   assert(nlexpr2origexpr != NULL);

   assert(SCIPexprGetNChildren(nlexpr) > 0);
   assert(SCIPexprGetChildren(nlexpr) != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   (*nlhdlrexprdata)->nlexpr = nlexpr;
   (*nlhdlrexprdata)->nlexpr2origexpr = nlexpr2origexpr;
   (*nlhdlrexprdata)->nleafs = nleafs;

   usingaux = FALSE;

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, nlexpr, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);

   for( ; !SCIPexpriterIsEnd(it); (void) SCIPexpriterGetNext(it) )
   {
      SCIP_EXPR* child;
      SCIP_EXPR* origexpr;

      /* check whether to-be-visited child needs to be replaced by a new expression (representing the auxvar)
       * if child has children, then that is not the case
       * if child has no children, but also corresponding origexpr has no chilren, then this is also not the case
       */
      child = SCIPexpriterGetChildExprDFS(it);
      if( SCIPexprGetNChildren(child) > 0 )
         continue;

      origexpr = (SCIP_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)child);
      assert(origexpr != NULL);

      /* if child had children in original but not in copy means that we could not achieve the desired curvature
       * thus, we will later replace by a new child that points to the auxvar of the original expression
       * as we do not have the auxvar now, we will only register that we will need the auxvar later (if origexpr isn't a variable or constant)
       * if we are working for the concave nlhdlr, then we also indicate interest on the exprs activity for estimate (distinguish below or above)
       */
      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, origexpr,
         SCIPexprGetNChildren(origexpr) > 0, FALSE,
         !nlhdlrdata->isnlhdlrconvex && (participating & SCIP_NLHDLR_METHOD_SEPABELOW),
         !nlhdlrdata->isnlhdlrconvex && (participating & SCIP_NLHDLR_METHOD_SEPAABOVE)) );

      /* remember that we use an auxvar */
      if( SCIPexprGetNChildren(origexpr) > 0 )
         usingaux = TRUE;
   }

   SCIPfreeExpriter(&it);

#ifdef SCIP_DEBUG
   SCIPprintExpr(scip, nlexpr, NULL);
   SCIPinfoMessage(scip, NULL, " (%p) is handled as %s\n", SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr),
         SCIPexprcurvGetName(SCIPexprGetCurvature(nlexpr)));
#endif

   /* If we don't work on the extended formulation, then set curvature also in original expression
    * (in case someone wants to pick this up; this might be removed again).
    * This doesn't ensure that every convex or concave original expression is actually marked here.
    * Not only because our tests are incomprehensive, but also because we may not detect on sums,
    * prefer extended formulations (in nlhdlr_convex), or introduce auxvars for linear subexpressions
    * on purpose (in nlhdlr_concave).
    */
   if( !usingaux )
      SCIPexprSetCurvature(expr, SCIPexprGetCurvature(nlexpr));

   return SCIP_OKAY;
}

/** adds an estimator for a vertex-polyhedral (e.g., concave) function to a given rowprep
 *
 * Calls \ref SCIPcomputeFacetVertexPolyhedralNonlinear() for given function and
 * box set to local bounds of auxiliary variables.
 */
static
SCIP_RETCODE estimateVertexPolyhedral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_SOL*             sol,                /**< solution to use, unless usemidpoint is TRUE */
   SCIP_Bool             usemidpoint,        /**< whether to use the midpoint of the domain instead of sol */
   SCIP_Bool             overestimate,       /**< whether over- or underestimating */
   SCIP_Real             targetvalue,        /**< a target value to achieve; if not reachable, then can give up early */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to store estimator */
   SCIP_Bool*            success             /**< buffer to store whether successful */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   VERTEXPOLYFUN_EVALDATA evaldata;
   SCIP_Real* xstar;
   SCIP_Real* box;
   SCIP_Real facetconstant;
   SCIP_VAR* var;
   int i;
   SCIP_Bool allfixed;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* caller is responsible to have checked whether we can estimate, i.e., expression curvature and overestimate flag match */
   assert( overestimate || SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONCAVE);  /* if underestimate, then must be concave */
   assert(!overestimate || SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONVEX);   /* if overestimate, then must be convex */

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "%sestimate expression ", overestimate ? "over" : "under");
   SCIPprintExpr(scip, nlhdlrexprdata->nlexpr, NULL);
   SCIPinfoMessage(scip, NULL, " at point\n");
   for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
   {
      var = SCIPgetVarExprVar(nlhdlrexprdata->leafexprs[i]);
      assert(var != NULL);

      SCIPinfoMessage(scip, NULL, "  <%s> = %g [%g,%g]\n", SCIPvarGetName(var),
         usemidpoint ? 0.5 * (SCIPvarGetLbLocal(var) + SCIPvarGetUbLocal(var)) : SCIPgetSolVal(scip, sol, var),
        SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
   }
#endif

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( nlhdlrdata->evalsol == NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, &nlhdlrdata->evalsol, NULL) );
   }

   evaldata.nlhdlrexprdata = nlhdlrexprdata;
   evaldata.evalsol = nlhdlrdata->evalsol;
   evaldata.scip = scip;

   SCIP_CALL( SCIPallocBufferArray(scip, &xstar, nlhdlrexprdata->nleafs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &box, 2*nlhdlrexprdata->nleafs) );

   allfixed = TRUE;
   for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
   {
      var = SCIPgetVarExprVar(nlhdlrexprdata->leafexprs[i]);
      assert(var != NULL);

      box[2*i] = SCIPvarGetLbLocal(var);
      if( SCIPisInfinity(scip, -box[2*i]) )
      {
         SCIPdebugMsg(scip, "lower bound at -infinity, no estimate possible\n");
         goto TERMINATE;
      }

      box[2*i+1] = SCIPvarGetUbLocal(var);
      if( SCIPisInfinity(scip, box[2*i+1]) )
      {
         SCIPdebugMsg(scip, "upper bound at +infinity, no estimate possible\n");
         goto TERMINATE;
      }

      if( !SCIPisRelEQ(scip, box[2*i], box[2*i+1]) )
         allfixed = FALSE;

      if( usemidpoint )
         xstar[i] = 0.5 * (box[2*i] + box[2*i+1]);
      else
         xstar[i] = SCIPgetSolVal(scip, sol, var);
      assert(xstar[i] != SCIP_INVALID);
   }

   if( allfixed )
   {
      /* SCIPcomputeFacetVertexPolyhedralNonlinear prints a warning and does not succeed if all is fixed */
      SCIPdebugMsg(scip, "all variables fixed, skip estimate\n");
      goto TERMINATE;
   }

   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nlhdlrexprdata->nleafs + 1) );

   SCIP_CALL( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, conshdlr, overestimate, nlhdlrExprEvalConcave, (void*)&evaldata,
      xstar, box, nlhdlrexprdata->nleafs, targetvalue, success, SCIProwprepGetCoefs(rowprep), &facetconstant) );

   if( !*success )
   {
      SCIPdebugMsg(scip, "failed to compute facet of convex hull\n");
      goto TERMINATE;
   }

   SCIProwprepSetLocal(rowprep, TRUE);
   SCIProwprepAddConstant(rowprep, facetconstant);
   for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
   {
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetVarExprVar(nlhdlrexprdata->leafexprs[i]), SCIProwprepGetCoefs(rowprep)[i]) );
   }

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "computed estimator: ");
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

 TERMINATE:
   SCIPfreeBufferArray(scip, &box);
   SCIPfreeBufferArray(scip, &xstar);

   return SCIP_OKAY;
}

/** adds an estimator computed via a gradient to a given rowprep */
static
SCIP_RETCODE estimateGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_SOL*             sol,                /**< solution to use */
   SCIP_Real             auxvalue,           /**< value of nlexpr in sol - we may not be able to take this value
                                                  from nlexpr if it was evaluated at a different sol recently */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to store estimator */
   SCIP_Bool*            success             /**< buffer to store whether successful */
   )
{
   SCIP_EXPR* nlexpr;
   SCIP_Real QUAD(constant);
   int i;

   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "estimate expression ");
   SCIPprintExpr(scip, nlexpr, NULL);
   SCIPinfoMessage(scip, NULL, " by gradient\n");
#endif

   *success = FALSE;

   /* evaluation error -> skip */
   if( auxvalue == SCIP_INVALID )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", auxvalue, (void*)nlexpr);
      return SCIP_OKAY;
   }

   /* compute gradient (TODO: this also re-evaluates (soltag=0), which shouldn't be necessary unless we tried ConvexSecant before) */
   SCIP_CALL( SCIPevalExprGradient(scip, nlexpr, sol, 0L) );

   /* gradient evaluation error -> skip */
   if( SCIPexprGetDerivative(nlexpr) == SCIP_INVALID )
   {
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)nlexpr);
      return SCIP_OKAY;
   }

   /* add gradient underestimator to rowprep: f(sol) + (x - sol) \nabla f(sol)
    * constant will store f(sol) - sol * \nabla f(sol)
    * to avoid some cancellation errors when linear variables take huge values (like 1e20),
    * we use double-double arithemtic here
    */
   QUAD_ASSIGN(constant, SCIPexprGetEvalValue(nlexpr)); /* f(sol) */
   for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real deriv;
      SCIP_Real varval;

      assert(SCIPexprGetDiffTag(nlhdlrexprdata->leafexprs[i]) == SCIPexprGetDiffTag(nlexpr));
      deriv = SCIPexprGetDerivative(nlhdlrexprdata->leafexprs[i]);
      if( deriv == SCIP_INVALID )
      {
         SCIPdebugMsg(scip, "gradient evaluation error for component %d of %p\n", i, (void*)nlexpr);
         return SCIP_OKAY;
      }

      var = SCIPgetVarExprVar(nlhdlrexprdata->leafexprs[i]);
      assert(var != NULL);

      varval = SCIPgetSolVal(scip, sol, var);

      SCIPdebugMsg(scip, "add %g * (<%s> - %g) to rowprep\n", deriv, SCIPvarGetName(var), varval);

      /* add deriv * var to rowprep and deriv * (-varval) to constant */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, deriv) );
      SCIPquadprecSumQD(constant, constant, -deriv * varval);
   }

   SCIProwprepAddConstant(rowprep, QUAD_TO_DBL(constant));
   SCIProwprepSetLocal(rowprep, FALSE);

   *success = TRUE;

   return SCIP_OKAY;
}

/** adds an estimator generated by putting a secant through the coordinates given by the two closest integer points */
static
SCIP_RETCODE estimateConvexSecant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_SOL*             sol,                /**< solution to use, unless usemidpoint is TRUE */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to store estimator */
   SCIP_Bool*            success             /**< buffer to store whether successful */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_EXPR* nlexpr;
   SCIP_VAR* var;
   SCIP_Real x;
   SCIP_Real left, right;
   SCIP_Real fleft, fright;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nleafs == 1);
   assert(rowprep != NULL);
   assert(success != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);

   *success = FALSE;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   var = SCIPgetVarExprVar(nlhdlrexprdata->leafexprs[0]);
   assert(var != NULL);

   x = SCIPgetSolVal(scip, sol, var);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "estimate expression ");
   SCIPprintExpr(scip, nlexpr, NULL);
   SCIPinfoMessage(scip, NULL, " by secant\n");
   SCIPinfoMessage(scip, NULL, "integral variable <%s> = %g [%g,%g]\n", SCIPvarGetName(var),
         x, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
#endif

   /* find out coordinates of var left and right to sol */
   if( SCIPisIntegral(scip, x) )
   {
      x = SCIPround(scip, x);
      if( SCIPisEQ(scip, x, SCIPvarGetLbGlobal(var)) )
      {
         left = x;
         right = left + 1.0;
      }
      else
      {
         right = x;
         left = right - 1.0;
      }
   }
   else
   {
      left = SCIPfloor(scip, x);
      right = SCIPceil(scip, x);
   }
   assert(left != right);

   /* now evaluate at left and right */
   if( nlhdlrdata->evalsol == NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, &nlhdlrdata->evalsol, NULL) );
   }

   SCIP_CALL( SCIPsetSolVal(scip, nlhdlrdata->evalsol, var, left) );
   SCIP_CALL( SCIPevalExpr(scip, nlexpr, nlhdlrdata->evalsol, 0L) );

   /* evaluation error or a too large constant -> skip */
   fleft = SCIPexprGetEvalValue(nlexpr);
   if( SCIPisInfinity(scip, REALABS(fleft)) )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", SCIPexprGetEvalValue(nlexpr), (void*)nlexpr);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsetSolVal(scip, nlhdlrdata->evalsol, var, right) );
   SCIP_CALL( SCIPevalExpr(scip, nlexpr, nlhdlrdata->evalsol, 0L) );

   /* evaluation error or a too large constant -> skip */
   fright = SCIPexprGetEvalValue(nlexpr);
   if( SCIPisInfinity(scip, REALABS(fright)) )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", SCIPexprGetEvalValue(nlexpr), (void*)nlexpr);
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "f(%g)=%g, f(%g)=%g\n", left, fleft, right, fright);

   /* skip if too steep
    * for clay0204h, this resulted in a wrong cut from f(0)=1e12 f(1)=0.99998,
    * since due to limited precision, this was handled as if f(1)=1
    */
   if( (!SCIPisZero(scip, fleft)  && REALABS(fright/fleft)*SCIPepsilon(scip) > 1.0) ||
       (!SCIPisZero(scip, fright) && REALABS(fleft/fright)*SCIPepsilon(scip) > 1.0) )
   {
      SCIPdebugMsg(scip, "function is too steep, abandoning\n");
      return SCIP_OKAY;
   }

   /* now add f(left) + (f(right) - f(left)) * (x - left) as estimator to rowprep */
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, fright - fleft) );
   SCIProwprepAddConstant(rowprep, fleft - (fright - fleft) * left);
   SCIProwprepSetLocal(rowprep, FALSE);

   *success = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of convex nonlinear handler
 */

/** free handler data of convex or concave nlhdlr */
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrfreeHdlrDataConvexConcave)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nlhdlrdata != NULL);
   assert(*nlhdlrdata != NULL);
   assert((*nlhdlrdata)->evalsol == NULL);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrfreeExprDataConvexConcave)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexprdata)->leafexprs, (*nlhdlrexprdata)->nleafs);
   SCIP_CALL( SCIPreleaseExpr(scip, &(*nlhdlrexprdata)->nlexpr) );
   SCIPhashmapFree(&(*nlhdlrexprdata)->nlexpr2origexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** deinitialization of problem-specific data */
static
SCIP_DECL_NLHDLREXIT(nlhdlrExitConvex)
{
   SCIP_NLHDLRDATA* nlhdlrdata;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( nlhdlrdata->evalsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &nlhdlrdata->evalsol) );
   }

   return SCIP_OKAY;
}

/** checks whether expression (or -expression) is convex, possibly after introducing auxiliary variables */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectConvex)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_EXPR* nlexpr = NULL;
   SCIP_HASHMAP* nlexpr2origexpr;
   int nleafs = 0;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);
   assert(nlhdlrexprdata != NULL);

   /* we currently do not participate if only activity computation is required */
   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) == SCIP_NLHDLR_METHOD_SEPABOTH )
      return SCIP_OKAY;

   /* ignore pure constants and variables */
   if( SCIPexprGetNChildren(expr) == 0 )
      return SCIP_OKAY;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);
   assert(nlhdlrdata->isnlhdlrconvex);

   SCIPdebugMsg(scip, "nlhdlr_convex detect for expr %p\n", (void*)expr);

   /* initialize mapping from copied expression to original one
    * 20 is not a bad estimate for the size of convex subexpressions that we can usually discover
    * when expressions will be allowed to store "user"data, we could get rid of this hashmap (TODO)
    */
   SCIP_CALL( SCIPhashmapCreate(&nlexpr2origexpr, SCIPblkmem(scip), 20) );

   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABELOW) == 0 )  /* if no separation below yet */
   {
      SCIP_CALL( constructExpr(scip, nlhdlrdata, &nlexpr, nlexpr2origexpr, &nleafs, expr,
         SCIP_EXPRCURV_CONVEX, NULL, SCIPassumeConvexNonlinear(conshdlr), NULL) );
      if( nlexpr != NULL )
      {
         assert(SCIPexprGetNChildren(nlexpr) > 0);  /* should not be trivial */

         *participating |= SCIP_NLHDLR_METHOD_SEPABELOW;

         SCIPdebugMsg(scip, "detected expr %p to be convex -> can enforce expr <= auxvar\n", (void*)expr);
      }
      else
      {
         SCIP_CALL( SCIPhashmapRemoveAll(nlexpr2origexpr) );
      }
   }

   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPAABOVE) == 0 && nlexpr == NULL )  /* if no separation above and not convex */
   {
      SCIP_CALL( constructExpr(scip, nlhdlrdata, &nlexpr, nlexpr2origexpr, &nleafs, expr,
         SCIP_EXPRCURV_CONCAVE, NULL, SCIPassumeConvexNonlinear(conshdlr), NULL) );
      if( nlexpr != NULL )
      {
         assert(SCIPexprGetNChildren(nlexpr) > 0);  /* should not be trivial */

         *participating |= SCIP_NLHDLR_METHOD_SEPAABOVE;

         SCIPdebugMsg(scip, "detected expr %p to be concave -> can enforce expr >= auxvar\n", (void*)expr);
      }
   }

   /* everything we participate in we also enforce */
   *enforcing |= *participating;

   assert(*participating || nlexpr == NULL);
   if( !*participating )
   {
      SCIPhashmapFree(&nlexpr2origexpr);
      return SCIP_OKAY;
   }

   /* create the expression data of the nonlinear handler
    * notify conshdlr about expr for which we will require auxiliary variables
    */
   SCIP_CALL( createNlhdlrExprData(scip, nlhdlrdata, nlhdlrexprdata, expr, nlexpr, nlexpr2origexpr, nleafs, *participating) );

   return SCIP_OKAY;
}

/** auxiliary evaluation callback */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalAuxConvexConcave)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlexpr != NULL);
   assert(auxvalue != NULL);

   SCIP_CALL( SCIPevalExpr(scip, nlhdlrexprdata->nlexpr, sol, 0L) );
   *auxvalue = SCIPexprGetEvalValue(nlhdlrexprdata->nlexpr);

   return SCIP_OKAY;
}

/** init sepa callback that initializes LP */
static
SCIP_DECL_NLHDLRINITSEPA(nlhdlrInitSepaConvex)
{  /*lint --e{715}*/
   SCIP_EXPR* nlexpr;
   SCIP_EXPRCURV curvature;
   SCIP_Bool success;
   SCIP_ROWPREP* rowprep = NULL;
   SCIP_ROW* row;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real lambda;
   SCIP_SOL* sol;
   int k;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   /* setup nlhdlrexprdata->leafexprs */
   SCIP_CALL( collectLeafs(scip, nlhdlrexprdata) );

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);
   assert(SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)nlexpr) == expr);

   curvature = SCIPexprGetCurvature(nlexpr);
   assert(curvature == SCIP_EXPRCURV_CONVEX || curvature == SCIP_EXPRCURV_CONCAVE);

   /* we can only be estimating on the convex side */
   if( curvature == SCIP_EXPRCURV_CONVEX )
      overestimate = FALSE;
   else if( curvature == SCIP_EXPRCURV_CONCAVE )
      underestimate = FALSE;
   if( !overestimate && !underestimate )
      return SCIP_OKAY;

   /* linearizes at 5 different points obtained as convex combination of the lower and upper bound of the variables
    * present in the convex expression; whether more weight is given to the lower or upper bound of a variable depends
    * on whether the fixing of the variable to that value is better for the objective function
    */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   *infeasible = FALSE;

   for( k = 0; k < 5; ++k )
   {
      int i;
      lambda = 0.1 * (k+1); /* lambda = 0.1, 0.2, 0.3, 0.4, 0.5 */

      for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
      {
         SCIP_VAR* var;

         var = SCIPgetVarExprVar(nlhdlrexprdata->leafexprs[i]);

         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);

         /* make sure the absolute values of bounds are not too large */
         if( ub > -INITLPMAXVARVAL )
            lb = MAX(lb, -INITLPMAXVARVAL);
         if( lb <  INITLPMAXVARVAL )
            ub = MIN(ub,  INITLPMAXVARVAL);

         /* in the case when ub < -maxabsbnd or lb > maxabsbnd, we still want to at least make bounds finite */
         if( SCIPisInfinity(scip, -lb) )
            lb = MIN(-10.0, ub - 0.1*REALABS(ub));
         if( SCIPisInfinity(scip,  ub) )
            ub = MAX( 10.0, lb + 0.1*REALABS(lb));

         if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, lambda * ub + (1.0 - lambda) * lb) );
         else
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, lambda * lb + (1.0 - lambda) * ub) );
      }

      SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );
      SCIP_CALL( estimateGradient(scip, nlhdlrexprdata, sol, 0.0, rowprep, &success) );
      if( !success )
      {
         SCIPdebugMsg(scip, "failed to linearize for k = %d\n", k);
         SCIPfreeRowprep(scip, &rowprep);
         continue;
      }

      /* add auxiliary variable */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetExprAuxVarNonlinear(expr), -1.0) );

      /* straighten out numerics */
      SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, NULL, SCIPgetHugeValue(scip), &success) );
      if( !success )
      {
         SCIPdebugMsg(scip, "failed to cleanup rowprep numerics for k = %d\n", k);
         SCIPfreeRowprep(scip, &rowprep);
         continue;
      }

      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "%sestimate_gradient%p_initsepa_%d",
            overestimate ? "over" : "under", (void*)expr, k);
      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );
      SCIPfreeRowprep(scip, &rowprep);

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "initsepa computed row: ");
      SCIPprintRow(scip, row, NULL);
#endif

      SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );

      if( *infeasible )
         break;
   }

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}

/** estimator callback */
static
SCIP_DECL_NLHDLRESTIMATE(nlhdlrEstimateConvex)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlexpr != NULL);
   assert(rowpreps != NULL);
   assert(success != NULL);

   assert(SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)nlhdlrexprdata->nlexpr) == expr);

   /* we must be called only for the side that we indicated to participate in during DETECT */
   assert(SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONVEX
          || SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONCAVE);
   assert(!overestimate || SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONCAVE);
   assert( overestimate || SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONVEX);

   *success = FALSE;
   *addedbranchscores = FALSE;

   /* we can skip eval as nlhdlrEvalAux should have been called for same solution before */
   /* SCIP_CALL( nlhdlrExprEval(scip, nlexpr, sol) ); */
   assert(auxvalue == SCIPexprGetEvalValue(nlhdlrexprdata->nlexpr)); /* given value (originally from
         nlhdlrEvalAuxConvexConcave) should coincide with the one stored in nlexpr */

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );

   if( nlhdlrexprdata->nleafs == 1 && SCIPexprIsIntegral(nlhdlrexprdata->leafexprs[0]) )
   {
      SCIP_CALL( estimateConvexSecant(scip, nlhdlr, nlhdlrexprdata, sol, rowprep, success) );

      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "%sestimate_convexsecant%p_%s%" SCIP_LONGINT_FORMAT,
         overestimate ? "over" : "under",
         (void*)expr,
         sol != NULL ? "sol" : "lp",
         sol != NULL ? (SCIP_Longint) SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));
   }

   /* if secant method was not used or failed, then try with gradient */
   if( !*success )
   {
      SCIP_CALL( estimateGradient(scip, nlhdlrexprdata, sol, auxvalue, rowprep, success) );

      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "%sestimate_convexgradient%p_%s%" SCIP_LONGINT_FORMAT,
         overestimate ? "over" : "under",
         (void*)expr,
         sol != NULL ? "sol" : "lp",
         sol != NULL ? (SCIP_Longint) SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));
   }

   if( *success )
   {
      SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep) );
   }
   else
   {
      SCIPfreeRowprep(scip, &rowprep);
   }

   return SCIP_OKAY;
}

/** include nlhdlr in another scip instance */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrConvex)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), CONVEX_NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrConvex(targetscip) );

   return SCIP_OKAY;
}

/** includes convex nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrConvex(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLR* nlhdlr;
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   nlhdlrdata->isnlhdlrconvex = TRUE;
   nlhdlrdata->evalsol = NULL;

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, CONVEX_NLHDLR_NAME, CONVEX_NLHDLR_DESC,
      CONVEX_NLHDLR_DETECTPRIORITY, CONVEX_NLHDLR_ENFOPRIORITY, nlhdlrDetectConvex, nlhdlrEvalAuxConvexConcave, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONVEX_NLHDLR_NAME "/detectsum",
      "whether to run convexity detection when the root of an expression is a non-quadratic sum",
      &nlhdlrdata->detectsum, FALSE, DEFAULT_DETECTSUM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONVEX_NLHDLR_NAME "/extendedform",
      "whether to create extended formulations instead of looking for maximal convex expressions",
      &nlhdlrdata->extendedform, FALSE, DEFAULT_EXTENDEDFORM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONVEX_NLHDLR_NAME "/cvxquadratic",
      "whether to use convexity check on quadratics",
      &nlhdlrdata->cvxquadratic, TRUE, DEFAULT_CVXQUADRATIC_CONVEX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONVEX_NLHDLR_NAME "/cvxsignomial",
      "whether to use convexity check on signomials",
      &nlhdlrdata->cvxsignomial, TRUE, DEFAULT_CVXSIGNOMIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONVEX_NLHDLR_NAME "/cvxprodcomp",
      "whether to use convexity check on product composition f(h)*h",
      &nlhdlrdata->cvxprodcomp, TRUE, DEFAULT_CVXPRODCOMP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONVEX_NLHDLR_NAME "/handletrivial",
      "whether to also handle trivial convex expressions",
      &nlhdlrdata->handletrivial, TRUE, DEFAULT_HANDLETRIVIAL, NULL, NULL) );

   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrfreeHdlrDataConvexConcave);
   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrConvex);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrfreeExprDataConvexConcave);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaConvex, NULL, nlhdlrEstimateConvex, NULL);
   SCIPnlhdlrSetInitExit(nlhdlr, NULL, nlhdlrExitConvex);

   return SCIP_OKAY;
}

/*
 * Callback methods of concave nonlinear handler
 */

/** deinitialization of problem-specific data */
static
SCIP_DECL_NLHDLREXIT(nlhdlrExitConcave)
{
   SCIP_NLHDLRDATA* nlhdlrdata;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( nlhdlrdata->evalsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &nlhdlrdata->evalsol) );
   }

   return SCIP_OKAY;
}

/** checks whether expression (or -expression) is concave, possibly after introducing auxiliary variables */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectConcave)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_EXPR* nlexpr = NULL;
   SCIP_HASHMAP* nlexpr2origexpr;
   int nleafs = 0;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);
   assert(nlhdlrexprdata != NULL);

   /* we currently do not participate if only activity computation is required */
   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) == SCIP_NLHDLR_METHOD_SEPABOTH )
      return SCIP_OKAY;

   /* ignore pure constants and variables */
   if( SCIPexprGetNChildren(expr) == 0 )
      return SCIP_OKAY;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);
   assert(!nlhdlrdata->isnlhdlrconvex);

   SCIPdebugMsg(scip, "nlhdlr_concave detect for expr %p\n", (void*)expr);

   /* initialize mapping from copied expression to original one
    * 20 is not a bad estimate for the size of concave subexpressions that we can usually discover
    * when expressions will be allowed to store "user"data, we could get rid of this hashmap (TODO)
    */
   SCIP_CALL( SCIPhashmapCreate(&nlexpr2origexpr, SCIPblkmem(scip), 20) );

   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABELOW) == 0 )  /* if no separation below yet */
   {
      SCIP_CALL( constructExpr(scip, nlhdlrdata, &nlexpr, nlexpr2origexpr, &nleafs, expr,
         SCIP_EXPRCURV_CONCAVE, NULL, FALSE, NULL) );

      if( nlexpr != NULL && nleafs > SCIP_MAXVERTEXPOLYDIM )
      {
         SCIPdebugMsg(scip, "Too many variables (%d) in constructed expression. Will not be able to estimate. Rejecting.\n", nleafs);
         SCIP_CALL( SCIPreleaseExpr(scip, &nlexpr) );
      }

      if( nlexpr != NULL )
      {
         assert(SCIPexprGetNChildren(nlexpr) > 0);  /* should not be trivial */

         *participating |= SCIP_NLHDLR_METHOD_SEPABELOW;

         SCIPdebugMsg(scip, "detected expr %p to be concave -> can enforce expr <= auxvar\n", (void*)expr);
      }
      else
      {
         SCIP_CALL( SCIPhashmapRemoveAll(nlexpr2origexpr) );
      }
   }

   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPAABOVE) == 0 && nlexpr == NULL )  /* if no separation above and not concave */
   {
      SCIP_CALL( constructExpr(scip, nlhdlrdata, &nlexpr, nlexpr2origexpr, &nleafs, expr,
         SCIP_EXPRCURV_CONVEX, NULL, FALSE, NULL) );

      if( nlexpr != NULL && nleafs > SCIP_MAXVERTEXPOLYDIM )
      {
         SCIPdebugMsg(scip, "Too many variables (%d) in constructed expression. Will not be able to estimate. Rejecting.\n", nleafs);
         SCIP_CALL( SCIPreleaseExpr(scip, &nlexpr) );
      }

      if( nlexpr != NULL )
      {
         assert(SCIPexprGetNChildren(nlexpr) > 0);  /* should not be trivial */

         *participating |= SCIP_NLHDLR_METHOD_SEPAABOVE;

         SCIPdebugMsg(scip, "detected expr %p to be convex -> can enforce expr >= auxvar\n", (void*)expr);
      }
   }

   /* everything we participate in we also enforce (at the moment) */
   *enforcing |= *participating;

   assert(*participating || nlexpr == NULL);
   if( !*participating )
   {
      SCIPhashmapFree(&nlexpr2origexpr);
      return SCIP_OKAY;
   }

   /* create the expression data of the nonlinear handler
    * notify conshdlr about expr for which we will require auxiliary variables and use activity
    */
   SCIP_CALL( createNlhdlrExprData(scip, nlhdlrdata, nlhdlrexprdata, expr, nlexpr, nlexpr2origexpr, nleafs, *participating) );

   return SCIP_OKAY;
}

/** init sepa callback that initializes LP */
static
SCIP_DECL_NLHDLRINITSEPA(nlhdlrInitSepaConcave)
{
   SCIP_EXPR* nlexpr;
   SCIP_EXPRCURV curvature;
   SCIP_Bool success;
   SCIP_ROWPREP* rowprep = NULL;
   SCIP_ROW* row;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);
   assert(SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)nlexpr) == expr);

   /* setup nlhdlrexprdata->leafexprs */
   SCIP_CALL( collectLeafs(scip, nlhdlrexprdata) );

   curvature = SCIPexprGetCurvature(nlexpr);
   assert(curvature == SCIP_EXPRCURV_CONVEX || curvature == SCIP_EXPRCURV_CONCAVE);
   /* we can only be estimating on non-convex side */
   if( curvature == SCIP_EXPRCURV_CONCAVE )
      overestimate = FALSE;
   else if( curvature == SCIP_EXPRCURV_CONVEX )
      underestimate = FALSE;
   if( !overestimate && !underestimate )
      return SCIP_OKAY;

   /* compute estimator and store in rowprep */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );
   SCIP_CALL( estimateVertexPolyhedral(scip, conshdlr, nlhdlr, nlhdlrexprdata, NULL, TRUE, overestimate,
         overestimate ? SCIPinfinity(scip) : -SCIPinfinity(scip), rowprep, &success) );
   if( !success )
   {
      SCIPdebugMsg(scip, "failed to compute facet of convex hull\n");
      goto TERMINATE;
   }

   /* add auxiliary variable */
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetExprAuxVarNonlinear(expr), -1.0) );

   /* straighten out numerics */
   SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, NULL, SCIPgetHugeValue(scip), &success) );
   if( !success )
   {
      SCIPdebugMsg(scip, "failed to cleanup rowprep numerics\n");
      goto TERMINATE;
   }

   (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "%sestimate_concave%p_initsepa",
         overestimate ? "over" : "under", (void*)expr);
   SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "initsepa computed row: ");
   SCIPprintRow(scip, row, NULL);
#endif

   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

 TERMINATE:
   if( rowprep != NULL )
      SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** estimator callback */
static
SCIP_DECL_NLHDLRESTIMATE(nlhdlrEstimateConcave)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlexpr != NULL);
   assert(rowpreps != NULL);
   assert(success != NULL);

   assert(SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)nlhdlrexprdata->nlexpr) == expr);

   /* we must be called only for the side that we indicated to participate in during DETECT */
   assert(SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONVEX
         || SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONCAVE);
   assert(!overestimate || SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONVEX);
   assert( overestimate || SCIPexprGetCurvature(nlhdlrexprdata->nlexpr) == SCIP_EXPRCURV_CONCAVE);

   *success = FALSE;
   *addedbranchscores = FALSE;

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );

   SCIP_CALL( estimateVertexPolyhedral(scip, conshdlr, nlhdlr, nlhdlrexprdata, sol, FALSE, overestimate, targetvalue, rowprep, success) );

   if( *success )
   {
      SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep) );

      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "%sestimate_concave%p_%s%" SCIP_LONGINT_FORMAT,
         overestimate ? "over" : "under",
         (void*)expr,
         sol != NULL ? "sol" : "lp",
         sol != NULL ? (SCIP_Longint) SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));
   }
   else
   {
      SCIPfreeRowprep(scip, &rowprep);
   }

   if( addbranchscores )
   {
      SCIP_Real violation;

      /* check how much is the violation on the side that we estimate */
      if( auxvalue == SCIP_INVALID )
      {
         /* if cannot evaluate, then always branch */
         violation = SCIPinfinity(scip);
      }
      else
      {
         SCIP_Real auxval;

         /* get value of auxiliary variable of this expression */
         assert(SCIPgetExprAuxVarNonlinear(expr) != NULL);
         auxval = SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr));

         /* compute the violation
          * if we underestimate, then we enforce expr <= auxval, so violation is (positive part of) auxvalue - auxval
          * if we overestimate,  then we enforce expr >= auxval, so violation is (positive part of) auxval - auxvalue
          */
         if( !overestimate )
            violation = MAX(0.0, auxvalue - auxval);
         else
            violation = MAX(0.0, auxval - auxvalue);
      }
      assert(violation >= 0.0);

      /* add violation as branching-score to expressions; the core will take care distributing this onto variables */
      if( nlhdlrexprdata->nleafs == 1 )
      {
         SCIP_EXPR* e;
         e = (SCIP_EXPR*)SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, nlhdlrexprdata->leafexprs[0]);
         SCIP_CALL( SCIPaddExprsViolScoreNonlinear(scip, &e, 1, violation, sol, addedbranchscores) );
      }
      else
      {
         SCIP_EXPR** exprs;
         int c;

         /* map leaf expressions back to original expressions
          * TODO do this once at end of detect and store in nlhdlrexprdata
          */
         SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nlhdlrexprdata->nleafs) );
         for( c = 0; c < nlhdlrexprdata->nleafs; ++c )
               exprs[c] = (SCIP_EXPR*)SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, nlhdlrexprdata->leafexprs[c]);

         SCIP_CALL( SCIPaddExprsViolScoreNonlinear(scip, exprs, nlhdlrexprdata->nleafs, violation, sol, addedbranchscores) );

         SCIPfreeBufferArray(scip, &exprs);
      }
   }

   return SCIP_OKAY;
}

/** includes nonlinear handler in another scip instance */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrConcave)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), CONCAVE_NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrConcave(targetscip) );

   return SCIP_OKAY;
}

/** includes concave nonlinear handler in nonlinear constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrConcave(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLR* nlhdlr;
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   nlhdlrdata->isnlhdlrconvex = FALSE;
   nlhdlrdata->evalsol = NULL;

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, CONCAVE_NLHDLR_NAME, CONCAVE_NLHDLR_DESC,
      CONCAVE_NLHDLR_DETECTPRIORITY, CONCAVE_NLHDLR_ENFOPRIORITY, nlhdlrDetectConcave, nlhdlrEvalAuxConvexConcave, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONCAVE_NLHDLR_NAME "/detectsum",
      "whether to run convexity detection when the root of an expression is a sum",
      &nlhdlrdata->detectsum, FALSE, DEFAULT_DETECTSUM, NULL, NULL) );

   /* "extended" formulations of a concave expressions can give worse estimators */
   nlhdlrdata->extendedform = FALSE;

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONCAVE_NLHDLR_NAME "/cvxquadratic",
      "whether to use convexity check on quadratics",
      &nlhdlrdata->cvxquadratic, TRUE, DEFAULT_CVXQUADRATIC_CONCAVE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONCAVE_NLHDLR_NAME "/cvxsignomial",
      "whether to use convexity check on signomials",
      &nlhdlrdata->cvxsignomial, TRUE, DEFAULT_CVXSIGNOMIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONCAVE_NLHDLR_NAME "/cvxprodcomp",
      "whether to use convexity check on product composition f(h)*h",
      &nlhdlrdata->cvxprodcomp, TRUE, DEFAULT_CVXPRODCOMP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" CONCAVE_NLHDLR_NAME "/handletrivial",
      "whether to also handle trivial convex expressions",
      &nlhdlrdata->handletrivial, TRUE, DEFAULT_HANDLETRIVIAL, NULL, NULL) );

   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrfreeHdlrDataConvexConcave);
   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrConcave);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrfreeExprDataConvexConcave);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaConcave, NULL, nlhdlrEstimateConcave, NULL);
   SCIPnlhdlrSetInitExit(nlhdlr, NULL, nlhdlrExitConcave);

   return SCIP_OKAY;
}

/** checks whether a given expression is convex or concave w.r.t. the original variables
 *
 * This function uses the methods that are used in the detection algorithm of the convex nonlinear handler.
 */
SCIP_RETCODE SCIPhasExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRCURV         curv,               /**< curvature to check for */
   SCIP_Bool*            success,            /**< buffer to store whether expression has curvature curv (w.r.t. original variables) */
   SCIP_HASHMAP*         assumevarfixed      /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   )
{
   SCIP_NLHDLRDATA nlhdlrdata;
   SCIP_EXPR* rootnlexpr;
   SCIP_HASHMAP* nlexpr2origexpr;
   int nleafs;

   assert(expr != NULL);
   assert(curv != SCIP_EXPRCURV_UNKNOWN);
   assert(success != NULL);

   /* create temporary hashmap */
   SCIP_CALL( SCIPhashmapCreate(&nlexpr2origexpr, SCIPblkmem(scip), 20) );

   /* prepare nonlinear handler data */
   nlhdlrdata.isnlhdlrconvex = TRUE;
   nlhdlrdata.evalsol = NULL;
   nlhdlrdata.detectsum = TRUE;
   nlhdlrdata.extendedform = FALSE;
   nlhdlrdata.cvxquadratic = TRUE;
   nlhdlrdata.cvxsignomial = TRUE;
   nlhdlrdata.cvxprodcomp = TRUE;
   nlhdlrdata.handletrivial = TRUE;

   SCIP_CALL( constructExpr(scip, &nlhdlrdata, &rootnlexpr, nlexpr2origexpr, &nleafs, expr, curv, assumevarfixed, FALSE, success) );

   /* free created expression */
   if( rootnlexpr != NULL )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &rootnlexpr) );
   }

   /* free hashmap */
   SCIPhashmapFree(&nlexpr2origexpr);

   return SCIP_OKAY;
}
