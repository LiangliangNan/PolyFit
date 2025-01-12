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

/**@file   struct_expr.h
 * @ingroup INTERNALAPI
 * @brief  structure definitions related to algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_STRUCT_EXPR_H_
#define SCIP_STRUCT_EXPR_H_

#include "scip/type_expr.h"
#include "scip/type_clock.h"
#include "scip/type_stat.h"
#include "blockmemshell/memory.h"

/** generic data and callback methods of an expression handler */
struct SCIP_Exprhdlr
{
   char*                 name;               /**< expression handler name */
   char*                 desc;               /**< expression handler description (can be NULL) */
   SCIP_EXPRHDLRDATA*    data;               /**< data of handler */
   unsigned int          precedence;         /**< precedence of expression operation relative to other expression (used for printing) */

   /* callbacks */
   SCIP_DECL_EXPRCOPYHDLR((*copyhdlr));      /**< handler copy callback (can be NULL) */
   SCIP_DECL_EXPRFREEHDLR((*freehdlr));      /**< handler free callback (can be NULL) */
   SCIP_DECL_EXPRCOPYDATA((*copydata));      /**< data copy callback, or NULL for expressions that have no data */
   SCIP_DECL_EXPRFREEDATA((*freedata));      /**< data free callback, or NULL for expressions that have no data or which data does not need to be freed */
   SCIP_DECL_EXPRSIMPLIFY((*simplify));      /**< simplify callback (can be NULL) */
   SCIP_DECL_EXPRCOMPARE((*compare));        /**< compare callback (can be NULL) */
   SCIP_DECL_EXPRPRINT((*print));            /**< print callback (can be NULL) */
   SCIP_DECL_EXPRPARSE((*parse));            /**< parse callback (can be NULL) */
   SCIP_DECL_EXPREVAL((*eval));              /**< point evaluation callback (can never be NULL) */
   SCIP_DECL_EXPRBWDIFF((*bwdiff));          /**< backward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRFWDIFF((*fwdiff));          /**< forward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRBWFWDIFF((*bwfwdiff));      /**< backward over forward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRINTEVAL((*inteval));        /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_EXPRESTIMATE((*estimate));      /**< estimation callback (can be NULL) */
   SCIP_DECL_EXPRINITESTIMATES((*initestimates)); /**< initial estimators callback (can be NULL) */
   SCIP_DECL_EXPRREVERSEPROP((*reverseprop));/**< reverse propagation callback (can be NULL) */
   SCIP_DECL_EXPRHASH((*hash));              /**< hash callback (can be NULL) */
   SCIP_DECL_EXPRCURVATURE((*curvature));    /**< curvature detection callback (can be NULL) */
   SCIP_DECL_EXPRMONOTONICITY((*monotonicity)); /**< monotonicity detection callback (can be NULL) */
   SCIP_DECL_EXPRINTEGRALITY((*integrality));/**< integrality detection callback (can be NULL) */

   /* statistics */
   unsigned int          ncreated;           /**< number of times expression has been created */
   SCIP_Longint          nestimatecalls;     /**< number of times the estimation callback was called */
   SCIP_Longint          nintevalcalls;      /**< number of times the interval evaluation callback was called */
   SCIP_Longint          npropcalls;         /**< number of times the propagation callback was called */
   SCIP_Longint          ncutoffs;           /**< number of cutoffs found so far by this expression handler */
   SCIP_Longint          ndomreds;           /**< number of domain reductions found so far by this expression handler */
   SCIP_Longint          nsimplifycalls;     /**< number of times the simplification callback was called */
   SCIP_Longint          nsimplified;        /**< number of times the simplification callback simplified */
   SCIP_Longint          nbranchscores;      /**< number of times branching scores were added by (or for) this expression handler */

   SCIP_CLOCK*           estimatetime;       /**< time used for estimation */
   SCIP_CLOCK*           intevaltime;        /**< time used for interval evaluation */
   SCIP_CLOCK*           proptime;           /**< time used for propagation */
   SCIP_CLOCK*           simplifytime;       /**< time used for expression simplification */
};

/** expression iteration data stored in an expression */
struct SCIP_ExprIterData
{
   SCIP_EXPR*             parent;            /**< parent expression in DFS iteration */
   int                    currentchild;      /**< child that is currently visited (or will be visited next) by DFS iteration */
   SCIP_Longint           visitedtag;        /**< tag to identify whether an expression has been visited already */
   SCIP_EXPRITER_USERDATA userdata;          /**< space for iterator user to store some (temporary) data */
};

/* private types */
typedef struct SCIP_QuadExpr  SCIP_QUADEXPR;     /**< representation of expression as quadratic */
typedef struct SCIP_QuadExpr_QuadTerm  SCIP_QUADEXPR_QUADTERM;  /**< a single term associated to a quadratic variable */
typedef struct SCIP_QuadExpr_BilinTerm SCIP_QUADEXPR_BILINTERM; /**< a single bilinear term */

/** an algebraic expression */
struct SCIP_Expr
{
   SCIP_EXPRHDLR*        exprhdlr;           /**< expression type (as pointer to its handler) */
   SCIP_EXPRDATA*        exprdata;           /**< expression data */

   int                   nchildren;          /**< number of children */
   int                   childrensize;       /**< length of children array */
   SCIP_EXPR**           children;           /**< children expressions */

   int                   nuses;              /**< reference counter */
   SCIP_EXPRITERDATA     iterdata[SCIP_EXPRITER_MAXNACTIVE];  /**< data for expression iterators */

   /* owner data */
   SCIP_EXPR_OWNERDATA*  ownerdata;          /**< data stored by owner of expression */
   SCIP_DECL_EXPR_OWNERFREE((*ownerfree));   /**< callback for freeing ownerdata */
   SCIP_DECL_EXPR_OWNERPRINT((*ownerprint)); /**< callback for printing ownerdata */
   SCIP_DECL_EXPR_OWNEREVALACTIVITY((*ownerevalactivity)); /**< callback for evaluating activity */

   /* point-evaluation and differentiation*/
   SCIP_Real             evalvalue;          /**< value of expression from last evaluation (corresponding to evaltag) */
   SCIP_Real             derivative;         /**< partial derivative of a "root path" w.r.t. this expression (see \ref SCIP_EXPR_DIFF "Differentiation methods") */
   SCIP_Real             dot;                /**< directional derivative of this expr (see \ref SCIP_EXPR_DIFF "Differentiation methods") */
   SCIP_Real             bardot;             /**< directional derivative of derivative of root (strictly speaking, a path) w.r.t this expr (see \ref SCIP_EXPR_DIFF "Differentiation methods") */
   SCIP_Longint          evaltag;            /**< tag of point for which the expression has been evaluated last, or 0 */
   SCIP_Longint          difftag;            /**< when computing partial derivatives of an expression w.r.t. a variable,
                                               *  the tag allows us to decide whether the expression depends on the variable;
                                               *  the tag will be checked in SCIPgetExprPartialDiff() */

   /* interval-evaluation (activity) */
   SCIP_INTERVAL         activity;           /**< activity of expression with respect to variable bounds */
   SCIP_Longint          activitytag;        /**< tag of variable bounds for which activity is valid */

   /* curvature information */
   SCIP_EXPRCURV         curvature;           /**< curvature of the expression w.r.t. bounds that have been used in the last curvature detection */

   /* integrality information */
   SCIP_Bool             isintegral;          /**< whether expression value is integral in feasible solutions */

   /* view expression as quadratic */
   SCIP_QUADEXPR*        quaddata;            /**< representation of expression as a quadratic, if checked and being quadratic */
   SCIP_Bool             quadchecked;         /**< whether it has been checked whether the expression is quadratic */
};

/** representation of an expression as quadratic */
struct SCIP_QuadExpr
{
   SCIP_Real             constant;           /**< a constant term */

   int                   nlinexprs;          /**< number of linear terms */
   SCIP_EXPR**           linexprs;           /**< expressions of linear terms */
   SCIP_Real*            lincoefs;           /**< coefficients of linear terms */

   int                   nquadexprs;         /**< number of expressions in quadratic terms */
   SCIP_QUADEXPR_QUADTERM* quadexprterms;    /**< array with quadratic expression terms */

   int                   nbilinexprterms;    /**< number of bilinear expressions terms */
   SCIP_QUADEXPR_BILINTERM* bilinexprterms;  /**< bilinear expression terms array */

   SCIP_Bool             allexprsarevars;    /**< whether all arguments (linexprs, quadexprterms[.].expr) are variable expressions */

   SCIP_EXPRCURV         curvature;          /**< curvature of the quadratic representation of the expression */
   SCIP_Bool             curvaturechecked;   /**< whether curvature has been checked */

   /* eigen decomposition information */
   SCIP_Bool             eigeninfostored;    /**< whether the eigen information is stored */
   SCIP_Real*            eigenvalues;        /**< eigenvalues of the Q matrix: size of nquadexprs */
   SCIP_Real*            eigenvectors;       /**< eigenvalues of the Q matrix: size of nquadexprs^2 */
};

/** a quadratic term associated to an expression */
struct SCIP_QuadExpr_QuadTerm
{
   SCIP_EXPR*            expr;               /**< expression of quadratic term */
   SCIP_Real             lincoef;            /**< linear coefficient of expression */
   SCIP_Real             sqrcoef;            /**< square coefficient of expression */

   int                   nadjbilin;          /**< number of bilinear terms this expression is involved in */
   int                   adjbilinsize;       /**< size of adjacent bilinear terms array */
   int*                  adjbilin;           /**< indices of associated bilinear terms */

   SCIP_EXPR*            sqrexpr;            /**< expression that was found to be the square of expr, or NULL if no square term (sqrcoef==0) */
};

/** a single bilinear term coef * expr1 * expr2
 *
 * except for temporary reasons, we assume that the index of var1 is smaller than the index of var2
 */
struct SCIP_QuadExpr_BilinTerm
{
   SCIP_EXPR*            expr1;              /**< first factor of bilinear term */
   SCIP_EXPR*            expr2;              /**< second factor of bilinear term */
   SCIP_Real             coef;               /**< coefficient of bilinear term */
   int                   pos2;               /**< position of expr2's quadexprterm in quadexprterms */

   SCIP_EXPR*            prodexpr;           /**< expression that was found to be the product of expr1 and expr2 */
};

/** expression iterator */
struct SCIP_ExprIter
{
   BMS_BLKMEM*           blkmem;             /**< block memory */
   SCIP_STAT*            stat;               /**< dynamic problem statistics */

   SCIP_Bool             initialized;        /**< whether the iterator has been initialized, that is, is in use */
   SCIP_EXPRITER_TYPE    itertype;           /**< type of expression iterator */
   SCIP_EXPR*            curr;               /**< current expression of the iterator */
   int                   iterindex;          /**< index of iterator data in expressions, or -1 if not using iterator data in expressions */
   SCIP_Longint          visitedtag;         /**< tag to mark and recognize an expression as visited, or 0 if not avoiding multiple visits */

   /* data for rtopological mode */
   SCIP_EXPR**           dfsexprs;           /**< DFS stack */
   int*                  dfsnvisited;        /**< number of visited children for each expression in the stack */
   int                   dfsnexprs;          /**< total number of expression in stack */
   int                   dfssize;            /**< size of DFS stack */

   /* data for BFS mode */
   SCIP_QUEUE*           queue;              /**< BFS queue */

   /* data for DFS mode */
   SCIP_EXPRITER_STAGE   dfsstage;           /**< current stage */
   unsigned int          stopstages;         /**< stages in which to interrupt iterator */
};

#endif /* SCIP_STRUCT_EXPR_H_ */
