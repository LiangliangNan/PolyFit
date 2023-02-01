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

/**@file   struct_expr.h
 * @brief  data definitions for expressions and expression trees
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_EXPRESSION_H__
#define __SCIP_STRUCT_EXPRESSION_H__

#include "scip/def.h"
#include "scip/type_misc.h"
#include "nlpi/type_expr.h"
#include "nlpi/type_exprinterpret.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** operator data of an expression */
union SCIP_ExprOpData 
{
   int                   intval;             /**< index of a variable or parameter or a constant integer value */
   SCIP_Real             dbl;                /**< a constant double value */
   void*                 data;               /**< pointer to some data structure */
};

/** arithmetic expression node */
struct SCIP_Expr 
{
   SCIP_EXPROP           op;                 /**< operator of the node */
   int                   nchildren;          /**< number of children */
   SCIP_EXPR**           children;           /**< children nodes */
   SCIP_EXPROPDATA       data;               /**< operator data */
};

/** expression tree */
struct SCIP_ExprTree 
{
   BMS_BLKMEM*           blkmem;             /**< block memory data structure */
   SCIP_EXPR*            root;               /**< root node expression of expression tree */
   int                   nvars;              /**< number of variables */
   void**                vars;               /**< mapping of variable indices to user variables, may be NULL */
   int                   nparams;            /**< number of parameters (modifiable constants) in expression */
   SCIP_Real*            params;             /**< current values for parameters, or NULL if no parameters */
   SCIP_EXPRINTDATA*     interpreterdata;    /**< data of expression interpreter (evaluator) */
};

/** data of quadratic expression: sum_i coef_i x_i y_i */
struct SCIP_ExprData_Quadratic 
{
   SCIP_Real             constant;           /**< constant term */
   SCIP_Real*            lincoefs;           /**< linear coefficients of children */
   SCIP_QUADELEM*        quadelems;          /**< quadratic elements */
   int                   nquadelems;         /**< number of quadratic elements */
   SCIP_Bool             sorted;             /**< are the quadratic elements sorted? */
};

/** data of polynomial expression: constant + sum_i monom_i */
struct SCIP_ExprData_Polynomial 
{
   SCIP_Real             constant;           /**< constant term of polynomial */
   SCIP_EXPRDATA_MONOMIAL** monomials;       /**< monomials that constitute the polynomial */
   int                   monomialssize;      /**< size of monomials array */
   int                   nmonomials;         /**< number of monomials */
   SCIP_Bool             sorted;             /**< are the monomials sorted? */
};

/** data of monomial in polynomial expression: coef * prod_i child_i^exponent_i
 * we allow for real values exponents here 
 */
struct SCIP_ExprData_Monomial 
{
   SCIP_Real             coef;               /**< coefficient of monomial */
   int                   factorssize;        /**< size of factors arrays */
   int                   nfactors;           /**< number of factors */
   int*                  childidxs;          /**< children corresponding to factors */
   SCIP_Real*            exponents;          /**< value of exponent for each factor */
   SCIP_Bool             sorted;             /**< are the factors sorted (by childidx)? */
};

/** data of a user-defined expression
 */
struct SCIP_ExprData_User
{
   SCIP_USEREXPRDATA*    userdata;           /**< user data for expression */
   SCIP_EXPRINTCAPABILITY evalcapability;    /**< capabilities of evaluation functions */
   SCIP_DECL_USEREXPREVAL    ((*eval));      /**< evaluation function */
   SCIP_DECL_USEREXPRINTEVAL ((*inteval));   /**< interval evaluation function */
   SCIP_DECL_USEREXPRCURV    ((*curv));      /**< curvature check function */
   SCIP_DECL_USEREXPRPROP    ((*prop));      /**< interval propagation function */
   SCIP_DECL_USEREXPRESTIMATE ((*estimate)); /**< under-/over-estimator function */
   SCIP_DECL_USEREXPRCOPYDATA ((*copydata)); /**< expression data copy function, or NULL if nothing to copy */
   SCIP_DECL_USEREXPRFREEDATA ((*freedata)); /**< expression data free function, or NULL if nothing to free */
   SCIP_DECL_USEREXPRPRINT ((*print));        /**< expression print function, or NULL for default string "user()" */
};

/** a node in an expression graph */
struct SCIP_ExprGraphNode
{
   /* definition of expression */
   SCIP_EXPROP         op;                   /**< operand of expression */
   SCIP_EXPROPDATA     data;                 /**< data of expression */

   /* location of expression in expression graph */
   int                 depth;                /**< depth of node in expression graph */
   int                 pos;                  /**< position of node in expression graph nodes array of depth depth*/

   /* children of expression */
   int                 nchildren;            /**< number of child nodes in expression graph */
   SCIP_EXPRGRAPHNODE** children;            /**< children nodes */

   /* parents of expression */
   int                 parentssize;          /**< length of parents array */
   int                 nparents;             /**< number of parent nodes in expression graph */
   SCIP_EXPRGRAPHNODE** parents;             /**< parent nodes */
   SCIP_Bool           parentssorted;        /**< whether the parents array is sorted */

   /* domain propagation */
   SCIP_INTERVAL       bounds;               /**< bounds on expression */
   SCIP_EXPRBOUNDSTATUS boundstatus;         /**< status of bounds */

   /* value */
   SCIP_Real           value;                /**< value of expression in last evaluation call */

   /* curvature */
   SCIP_EXPRCURV       curv;                 /**< curvature of expression */

   /* miscellaneous */
   SCIP_Bool           enabled;              /**< whether the node is enabled currently */
   SCIP_Bool           simplified;           /**< whether the node has been simplified */
   int                 nuses;                /**< how often node is used */
};

/** a set of expression trees, stored in a single directed acyclic graph
 * the variables of the graph are stored at depth 0
 * for each depth, an array of nodes is stored */
struct SCIP_ExprGraph
{
   BMS_BLKMEM*          blkmem;              /**< block memory */

   int                  depth;               /**< depth of expression graph */
   int*                 nodessize;           /**< current size of nodes array for each depth */
   int*                 nnodes;              /**< number of nodes for each depth */
   SCIP_EXPRGRAPHNODE*** nodes;              /**< nodes of expression graph for each depth */

   int                  varssize;            /**< length of vars array */
   int                  nvars;               /**< number of variables in expression graph */
   void**               vars;                /**< array for variables in expression graph, having length varssize */
   SCIP_EXPRGRAPHNODE** varnodes;            /**< nodes corresponding to variables in expression graph */
   SCIP_INTERVAL*       varbounds;           /**< current bounds on variables */
   SCIP_HASHMAP*        varidxs;             /**< maps variables to variable indices */

   int                  constssize;          /**< length of consts array */
   int                  nconsts;             /**< number of constants in expression graph */
   SCIP_EXPRGRAPHNODE** constnodes;          /**< nodes corresponding to constants in expression graph */
   SCIP_Bool            constssorted;        /**< whether the constnodes array has been sorted */

   SCIP_DECL_EXPRGRAPHVARADDED((*exprgraphvaradded)); /**< callback for variable addition event */
   SCIP_DECL_EXPRGRAPHVARREMOVE((*exprgraphvarremove)); /**< callback for variable removal event */
   SCIP_DECL_EXPRGRAPHVARCHGIDX((*exprgraphvarchgidx)); /**< callback for variable index change event */
   void*                userdata;            /**< user data associated with callback methods */

   SCIP_Bool            needvarboundprop;    /**< whether variable bounds need be propagated, e.g., because new nodes have been added to the graph */
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_EXPRESSION_H__ */
