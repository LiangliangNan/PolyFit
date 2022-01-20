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

/**@file   type_expr.h
 * @brief  type definitions for expressions and expression trees
 * @ingroup TYPEDEFINITIONS
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_EXPRESSION_H__
#define __SCIP_TYPE_EXPRESSION_H__

#ifdef __cplusplus
extern "C" {
#endif

/** Operators of expressions.
 */
enum SCIP_ExprOp {
   /**@name Terminals (Leaves) */
   /**@{ */
   SCIP_EXPR_VARIDX    =  1,  /**< variable given by index (stored in data.idx) */
   SCIP_EXPR_CONST     =  2,  /**< constant (value stored in data.dbl) */
   SCIP_EXPR_PARAM     =  3,  /**< parameter = a constant that can be modified (should not be simplified away) */
   /**@} */

   /**@name Simple Operands */
   /**@{ */
   SCIP_EXPR_PLUS      =  8,  /**< addition (2 operands) */
   SCIP_EXPR_MINUS     =  9,  /**< substraction (2 operands) */
   SCIP_EXPR_MUL       = 10,  /**< multiplication (2 operands) */
   SCIP_EXPR_DIV       = 11,  /**< division (2 operands) */
   SCIP_EXPR_SQUARE    = 12,  /**< square (1 operand) */
   SCIP_EXPR_SQRT      = 13,  /**< square root (1 operand) */
   SCIP_EXPR_REALPOWER = 14,  /**< power with real exponent (1 operand!, assumed to be nonnegative, exponent is stored in expression data) */
   SCIP_EXPR_INTPOWER  = 15,  /**< power with integer exponent (1 operand!, exponent stored in expression data) */
   SCIP_EXPR_SIGNPOWER = 16,  /**< signed power (sign(x)|x|^a, 1 operand!, exponent is stored in expression data) */
   SCIP_EXPR_EXP       = 17,  /**< exponential (e^x, 1 operand) */
   SCIP_EXPR_LOG       = 18,  /**< natural logarithm (ln(x), 1 operand) */
   SCIP_EXPR_SIN       = 19,  /**< sinus (1 operand) */
   SCIP_EXPR_COS       = 20,  /**< cosinus (1 operand) */
   SCIP_EXPR_TAN       = 21,  /**< tangent (1 operand) */
   /* SCIP_EXPR_ERF       = 22, */  /**< gaussian error function (1 operand) */
   /* SCIP_EXPR_ERFI      = 23, */  /**< imaginary part of gaussian error function (1 operand) */
   SCIP_EXPR_MIN       = 24,  /**< minimum (2 operands) */
   SCIP_EXPR_MAX       = 25,  /**< maximum (2 operands) */
   SCIP_EXPR_ABS       = 26,  /**< absolute value (1 operand) */
   SCIP_EXPR_SIGN      = 27,  /**< sign of value (1 operand) */
   /**@} */

   /**@name Complex Operands
    */
   /**@{ */
   SCIP_EXPR_SUM       = 64,  /**< summation sum_{i=1}^n op_i (n operands) */
   SCIP_EXPR_PRODUCT   = 65,  /**< product prod_{i=1}^n op_i (n operands) */
   SCIP_EXPR_LINEAR    = 66,  /**< linear term sum_{i=1}^n a_i op_i (n operands) */
   SCIP_EXPR_QUADRATIC = 67,  /**< quadratic term sum_{i,j=1}^n a_{i,j} op_i op_j (n operands) */
   SCIP_EXPR_POLYNOMIAL= 68,  /**< polynomial term sum_{I} a_{I}ops^I (I a multiindex, n operands) */
   SCIP_EXPR_USER      = 69,  /**< a user defined expression */
   /**@} */

   SCIP_EXPR_LAST      = 70   /**< no expression, used for counting reasons */
};

/** Curvature types */
enum SCIP_ExprCurv
{
   SCIP_EXPRCURV_UNKNOWN    = 0,             /**< unknown curvature (or indefinite) */
   SCIP_EXPRCURV_CONVEX     = 1,             /**< convex */
   SCIP_EXPRCURV_CONCAVE    = 2,             /**< concave */
   SCIP_EXPRCURV_LINEAR     = SCIP_EXPRCURV_CONVEX | SCIP_EXPRCURV_CONCAVE/**< linear = convex and concave */
};

typedef enum   SCIP_ExprOp      SCIP_EXPROP;     /**< expression operand */
typedef union  SCIP_ExprOpData  SCIP_EXPROPDATA; /**< expression operand data */
typedef struct SCIP_Expr        SCIP_EXPR;       /**< expression */
typedef struct SCIP_ExprTree    SCIP_EXPRTREE;   /**< expression tree */
typedef enum   SCIP_ExprCurv    SCIP_EXPRCURV;   /**< curvature types */

/** An element of a quadratic term: two variable indices and a coefficient.
 * The convention is to have idx1 <= idx2.
 */
struct SCIP_QuadElement
{
   int                   idx1;             /**< index of first variable */
   int                   idx2;             /**< index of second variable */
   SCIP_Real             coef;             /**< value of coefficient at position (idx1, idx2) */
};
/* We have defined struct SCIP_QuadElement here (instead of type_expression.h) to allow fast access, allocation, and copying. (similar to SCIP_INTERVAL) */

typedef struct SCIP_QuadElement        SCIP_QUADELEM;           /**< element of a quadratic term */
typedef struct SCIP_ExprData_Quadratic SCIP_EXPRDATA_QUADRATIC; /**< the data of a quadratic expression (SCIP_EXPR_QUADRATIC) */

typedef struct SCIP_ExprData_Monomial   SCIP_EXPRDATA_MONOMIAL;   /**< a monomial as part of the data in a polynomial expression */
typedef struct SCIP_ExprData_Polynomial SCIP_EXPRDATA_POLYNOMIAL; /**< the data of a polynomial expression (SCIP_EXPR_POLYNOMIAL) */

typedef struct SCIP_ExprData_User      SCIP_EXPRDATA_USER;      /**< expression data of a user expression (not the user-data of a user expression) */

#define SCIP_EXPR_DEGREEINFINITY 65535       /**< value that stands for an infinite degree of an expression (see SCIPexprGetMaxDegree) */

/** signature of an expression (pointwise) evaluation function
 * The function should return nan, inf, or -inf in result if the function is undefined for the given arguments.
 *
 * - opdata    operand data
 * - nargs     number of arguments
 * - argvals   values of arguments
 * - varvals   values for variables
 * - paramvals values for parameters
 * - result    buffer where to store result of evaluation
 */
#define SCIP_DECL_EXPREVAL(x) SCIP_RETCODE x (SCIP_EXPROPDATA opdata, int nargs, SCIP_Real* argvals, SCIP_Real* varvals, SCIP_Real* paramvals, SCIP_Real* result)

/** signature of an expression (interval) evaluation function
 * The function should return an empty interval if the function is undefined for the given arguments.
 *
 * - infinity  value for infinity
 * - opdata    operand data
 * - nargs     number of arguments
 * - argvals   interval values of arguments
 * - varvals   interval values for variables
 * - paramvals values for parameters
 * - result    buffer where to store result of evaluation
 */
#define SCIP_DECL_EXPRINTEVAL(x) SCIP_RETCODE x (SCIP_Real infinity, SCIP_EXPROPDATA opdata, int nargs, SCIP_INTERVAL* argvals, SCIP_INTERVAL* varvals, SCIP_Real* paramvals, SCIP_INTERVAL* result)

/** signature of a simple expression curvature check function
 *
 * - infinity  value for infinity
 * - opdata    operand data
 * - nargs     number of arguments
 * - argbounds bounds on value of arguments
 * - argcurv   curvature of arguments
 * - paramvals values for parameters
 * - result    buffer where to store result of curvature check
 */
#define SCIP_DECL_EXPRCURV(x) SCIP_RETCODE x (SCIP_Real infinity, SCIP_EXPROPDATA opdata, int nargs, SCIP_INTERVAL* argbounds, SCIP_EXPRCURV* argcurv, SCIP_EXPRCURV* result)

/** signature of a expression data copy function
 *
 * - blkmem       block memory
 * - nchildren    number of children in expression
 * - opdatasource source expression data
 * - opdatatarget pointer to target expression data
 */
#define SCIP_DECL_EXPRCOPYDATA(x) SCIP_RETCODE x (BMS_BLKMEM* blkmem, int nchildren, SCIP_EXPROPDATA opdatasource, SCIP_EXPROPDATA* opdatatarget)

/** signature of a expression data free function
 *
 * - blkmem       block memory
 * - nchildren    number of children in expression
 * - opdata       expression data to free
 */
#define SCIP_DECL_EXPRFREEDATA(x) void x (BMS_BLKMEM* blkmem, int nchildren, SCIP_EXPROPDATA opdata)

typedef struct SCIP_ExprGraphNode SCIP_EXPRGRAPHNODE; /**< node in an expression graph */
typedef struct SCIP_ExprGraph     SCIP_EXPRGRAPH;     /**< an expression graph (DAG) */

/** callback method of expression graph invoked when a new variable has been added to the graph
 *
 * input:
 * - exprgraph    expression graph
 * - userdata     a pointer to user data
 * - var          variable that has been added to expression graph
 * - varnode      new expression graph node for a variable
 */
#define SCIP_DECL_EXPRGRAPHVARADDED(x) SCIP_RETCODE x (SCIP_EXPRGRAPH* exprgraph, void* userdata, void* var, SCIP_EXPRGRAPHNODE* varnode)

/** callback method of expression graph invoked when a variable is to be removed from the graph
 *
 * input:
 * - exprgraph    expression graph
 * - userdata     a pointer to user data
 * - var          variable that will be removed from the expression graph
 * - varnode      expression graph node corresponding to variable
 */
#define SCIP_DECL_EXPRGRAPHVARREMOVE(x) SCIP_RETCODE x (SCIP_EXPRGRAPH* exprgraph, void* userdata, void* var, SCIP_EXPRGRAPHNODE* varnode)

/** callback method of expression graph invoked when a variable changes its index
 *
 * input:
 * - exprgraph    expression graph
 * - userdata     a pointer to user data
 * - var          variable which will change its index
 * - varnode      expression graph node corresponding to variable
 * - oldidx       current index of variable
 * - newidx       new index the variable will have
 */
#define SCIP_DECL_EXPRGRAPHVARCHGIDX(x) SCIP_RETCODE x (SCIP_EXPRGRAPH* exprgraph, void* userdata, void* var, SCIP_EXPRGRAPHNODE* varnode, int oldidx, int newidx)

/* SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTFORCE is used to indicate that bounds in a node should be propagated down even if they are not better than the bounds given by the child nodes
 * this is useful if the expression itself has a domain that is not the whole real numbers and we want to propagate this information down to a child node
 * e.g., a x^0.3 should result in x >= 0 even if no new bounds on the expression x^0.3 have been obtained
 */
#define SCIP_EXPRBOUNDSTATUS_VALID             0x0 /**< bounds are valid, i.e., conform with bounds of children */
#define SCIP_EXPRBOUNDSTATUS_CHILDTIGHTENED    0x1 /**< a child bounds were tightened since last calculation */
#define SCIP_EXPRBOUNDSTATUS_CHILDRELAXED      0x2 /**< bounds are not valid and need to be recomputed, because the bounds in a child were relaxed */
#define SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENT 0x4 /**< bounds have been tightened by reverse propagation in a parent, they are valid as long as there has been no relaxation of bounds somewhere in the graph */
#define SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTRECENT (0x8 | SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENT) /**< bounds have recently been tightened by reverse propagation in a parent, this tightening has not been propagated further down yet */
#define SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTFORCE (0x10 | SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTRECENT) /**< bounds may have recently been tightened by reverse propagation in a parent, in any case we want to propagate bounds further down */

typedef char             SCIP_EXPRBOUNDSTATUS;     /**< bitflags that indicate the status of bounds stored in a node of an expression graph */


typedef struct SCIP_UserExprData SCIP_USEREXPRDATA; /**< the user data of a user expression */

/** signature of an user's expression under/over estimation function
 * The function should return nan, inf, or -inf in result if the function is undefined for the given arguments.
 *
 * - infinity      value for infinity
 * - data          user expression data
 * - nargs         number of arguments
 * - argvals       values of arguments
 * - argbounds     bounds on value of arguments
 * - overestimate  flag indicating whether to over- or under estimate the expression
 * - coeffs        buffer where to store resulting coeffs of arguments for the estimator
 * - constant      buffer where to store resulting constant of the estimator
 * - success       buffer to indicate whether under-/overestimation was successful
 */
#define SCIP_DECL_USEREXPRESTIMATE(x) SCIP_RETCODE x (SCIP_Real infinity, SCIP_USEREXPRDATA* data, int nargs, SCIP_Real* argvals, SCIP_INTERVAL* argbounds, SCIP_Bool overestimate, SCIP_Real* coeffs, SCIP_Real* constant, SCIP_Bool *success)

/** signature of an user's expression (pointwise) evaluation function
 * The function should return nan, inf, or -inf in result if the function is undefined for the given arguments.
 *
 * - data      user expression data
 * - nargs     number of arguments
 * - argvals   values of arguments
 * - funcvalue buffer where to store result of function evaluation
 * - gradient  buffer where to store result of gradient evaluation (NULL if not requested)
 * - hessian   buffer where to store result of Hessian evaluation (NULL if not requested, currently full dense matrix)
 */
#define SCIP_DECL_USEREXPREVAL(x) SCIP_RETCODE x (SCIP_USEREXPRDATA* data, int nargs, SCIP_Real* argvals, SCIP_Real* funcvalue, SCIP_Real* gradient, SCIP_Real* hessian)

/** signature of an user's expression (interval) evaluation function
 * The function should return an empty interval if the function is undefined for the given arguments.
 *
 * - infinity  value for infinity
 * - data      user expression data
 * - nargs     number of arguments
 * - argvals   interval values of arguments
 * - funvalue  buffer where to store result of function evaluation
 * - gradient  buffer where to store result of gradient evaluation (NULL if not requested)
 * - hessian   buffer where to store result of Hessian evaluation (NULL if not requested, currently full dense matrix)
 */
#define SCIP_DECL_USEREXPRINTEVAL(x) SCIP_RETCODE x (SCIP_Real infinity, SCIP_USEREXPRDATA* data, int nargs, SCIP_INTERVAL* argvals, SCIP_INTERVAL* funcvalue, SCIP_INTERVAL* gradient, SCIP_INTERVAL* hessian)

/** signature of a user's expression curvature check function
 *
 * - infinity  value for infinity
 * - data      user expression data
 * - nargs     number of arguments
 * - argbounds bounds on value of arguments
 * - argcurv   curvature of arguments
 * - result    buffer where to store result of curvature check
 */
#define SCIP_DECL_USEREXPRCURV(x) SCIP_RETCODE x (SCIP_Real infinity, SCIP_USEREXPRDATA* data, int nargs, SCIP_INTERVAL* argbounds, SCIP_EXPRCURV* argcurv, SCIP_EXPRCURV* result)

/** signature of an user's expression interval propagation function
 * The function should compute intervals of the arguments given an interval for the function itself and all arguments.
 *
 * - infinity    value for infinity
 * - data        user expression data
 * - nargs       number of arguments
 * - argbounds   bounds on values of arguments (on output: tightened bounds)
 * - funcbounds  bounds on function value
 * - cutoff      buffer to indicate whether an empty child interval was found
 */
#define SCIP_DECL_USEREXPRPROP(x) SCIP_RETCODE x (SCIP_Real infinity, SCIP_USEREXPRDATA* data, int nargs, SCIP_INTERVAL* argbounds, SCIP_INTERVAL funcbounds, SCIP_Bool* cutoff)

/** signature of a user's expression data copy function
 *
 * - blkmem       block memory
 * - nchildren    number of children in expression
 * - datasource   source user expression data
 * - datatarget   target user expression data
 */
#define SCIP_DECL_USEREXPRCOPYDATA(x) SCIP_RETCODE x (BMS_BLKMEM* blkmem, int nchildren, SCIP_USEREXPRDATA* datasource, SCIP_USEREXPRDATA** datatarget)

/** signature of a user's expression data free function
 *
 * - blkmem       block memory
 * - nchildren    number of children in expression
 * - data         user expression data to free
 */
#define SCIP_DECL_USEREXPRFREEDATA(x) void x (BMS_BLKMEM* blkmem, int nchildren, SCIP_USEREXPRDATA* data)

/** signature of a user's expression print function
 * The function should print the user expression's name that prepends the list of arguments "(x1,x2,...)". If not specified, only "user" is printed.
 *
 * - data         user expression data
 * - messagehdlr  SCIP message handler
 * - file         output file, or NULL if standard output should be used
 */
#define SCIP_DECL_USEREXPRPRINT(x) void x (SCIP_USEREXPRDATA* data, SCIP_MESSAGEHDLR*  messagehdlr, FILE* file)

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_EXPRESSION_H__ */
