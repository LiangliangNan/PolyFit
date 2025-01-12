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

/**@file   expr_product.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  product expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/pub_expr.h"
#include "scip/expr_product.h"
#include "scip/expr_sum.h"
#include "scip/expr_pow.h"
#include "scip/expr_value.h"
#include "scip/expr_exp.h"
#include "scip/expr_abs.h"
#include "scip/expr_entropy.h"
#include "scip/cons_nonlinear.h"
#include "scip/pub_misc.h"
#include "scip/nlhdlr_bilinear.h"

#define EXPRHDLR_NAME         "prod"
#define EXPRHDLR_DESC         "product expression"
#define EXPRHDLR_PRECEDENCE   50000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(54949.0)

/** macro to activate/deactivate debugging information of simplify method */
/*lint -emacro(681,debugSimplify) */
/*lint -emacro(506,debugSimplify) */
/*lint -emacro(774,debugSimplify) */
#ifdef SIMPLIFY_DEBUG
#define debugSimplify                   printf
#else
#define debugSimplify                   while( FALSE ) printf
#endif


/*lint -e777*/

/*
 * Data structures
 */

/** expression data */
struct SCIP_ExprData
{
   SCIP_Real             coefficient;        /**< coefficient */
};

struct SCIP_ExprhdlrData
{
   SCIP_CONSHDLR*        conshdlr;           /**< nonlinear constraint handler (to compute estimates for > 2-dim products) */
};

/** node for linked list of expressions */
struct exprnode
{
   SCIP_EXPR*            expr;               /**< expression in node */
   struct exprnode*      next;               /**< next node */
};

typedef struct exprnode EXPRNODE;

/*
 * Local methods
 */

/** evaluation callback for (vertex-polyhedral) functions used as input for facet computation of its envelopes */
static
SCIP_DECL_VERTEXPOLYFUN(prodfunction)
{
   /* funcdata is a pointer to the double holding the coefficient */
   SCIP_Real ret = *(SCIP_Real*)funcdata;
   int i;

   for( i = 0; i < nargs; ++i )
      ret *= args[i];

   return ret;
}

static
SCIP_RETCODE buildSimplifiedProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             simplifiedcoef,     /**< simplified product should be simplifiedcoef * PI simplifiedfactors */
   EXPRNODE**            simplifiedfactors,  /**< factors of simplified product */
   SCIP_Bool             changed,            /**< indicates whether some of the simplified factors was changed */
   SCIP_EXPR**           simplifiedexpr,     /**< buffer to store the simplified expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/*  methods for handling linked list of expressions */

/** inserts newnode at beginning of list */
static
void insertFirstList(
   EXPRNODE*             newnode,            /**< node to insert */
   EXPRNODE**            list                /**< list */
   )
{
   assert(list != NULL);
   assert(newnode != NULL);

   newnode->next = *list;
   *list = newnode;
}

/** removes first element of list and returns it */
static
EXPRNODE* listPopFirst(
   EXPRNODE**            list                /**< list */
   )
{
   EXPRNODE* first;

   assert(list != NULL);

   if( *list == NULL )
      return NULL;

   first = *list;
   *list = (*list)->next;
   first->next = NULL;

   return first;
}

/** returns length of list */
static
int listLength(
   EXPRNODE*             list                /**< list */
   )
{
   int length;

   if( list == NULL )
      return 0;

   length = 1;
   while( (list=list->next) != NULL )
      ++length;

   return length;
}

/** creates expression node and captures expression */
static
SCIP_RETCODE createExprNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression stored at node */
   EXPRNODE**            newnode             /**< pointer to store node */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, newnode) );

   (*newnode)->expr = expr;
   (*newnode)->next = NULL;
   SCIPcaptureExpr(expr);

   return SCIP_OKAY;
}

/** creates expression list from expressions */
static
SCIP_RETCODE createExprlistFromExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           exprs,              /**< expressions stored in list */
   int                   nexprs,             /**< number of expressions */
   EXPRNODE**            list                /**< pointer to store list */
   )
{
   int i;

   assert(*list == NULL);
   assert(nexprs > 0);

   debugSimplify("building expr list from %d expressions\n", nexprs);
   for( i = nexprs - 1; i >= 0; --i )
   {
      EXPRNODE* newnode;

      SCIP_CALL( createExprNode(scip, exprs[i], &newnode) );
      insertFirstList(newnode, list);
   }
   assert(nexprs > 1 || (*list)->next == NULL);

   return SCIP_OKAY;
}

/** frees expression node and releases expressions */
static
SCIP_RETCODE freeExprNode(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE**            node                /**< node to be freed */
   )
{
   assert(node != NULL && *node != NULL);

   SCIP_CALL( SCIPreleaseExpr(scip, &(*node)->expr) );
   SCIPfreeBlockMemory(scip, node);

   return SCIP_OKAY;
}

/** frees an expression list */
static
SCIP_RETCODE freeExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE**            exprlist            /**< list */
   )
{
   EXPRNODE* current;

   if( *exprlist == NULL )
      return SCIP_OKAY;

   current = *exprlist;
   while( current != NULL )
   {
      EXPRNODE* tofree;

      tofree = current;
      current = current->next;
      SCIP_CALL( freeExprNode(scip, &tofree) );
   }
   assert(current == NULL);
   *exprlist = NULL;

   return SCIP_OKAY;
}

/* helper functions for simplifying expressions */

/** creates a product expression with the elements of exprlist as its children */
static
SCIP_RETCODE createExprProductFromExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE*             exprlist,           /**< list containing the children of expr */
   SCIP_Real             coef,               /**< coef of expr */
   SCIP_EXPR**           expr,               /**< pointer to store the product expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   int i;
   int nchildren;
   SCIP_EXPR** children;

   /* asserts SP8 */
   assert(coef == 1.0);
   nchildren = listLength(exprlist);

   SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      children[i] = exprlist->expr;
      exprlist = exprlist->next;
   }

   assert(exprlist == NULL);

   SCIP_CALL( SCIPcreateExprProduct(scip, expr, nchildren, children, coef, ownercreate, ownercreatedata) );

   SCIPfreeBufferArray(scip, &children);

   return SCIP_OKAY;
}

/** simplifies a factor of a product expression: base, so that it is a valid children of a simplified product expr
 *
 * @note In contrast to other simplify methods, this does *not* return a simplified expression.
 * Instead, the method is intended to be called only when simplifying a product expression.
 * Since in general, base is not a simplified child of a product expression, this method returns
 * a list of expressions L, such that (prod L) = baset *and* each expression in L
 * is a valid child of a simplified product expression.
 */
static
SCIP_RETCODE simplifyFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            factor,             /**< expression to be simplified */
   SCIP_Real*            simplifiedcoef,     /**< coefficient of parent product expression */
   EXPRNODE**            simplifiedfactor,   /**< pointer to store the resulting expression node/list of nodes */
   SCIP_Bool*            changed             /**< pointer to store if some term actually got simplified */
   )
{
   assert(simplifiedfactor != NULL);
   assert(*simplifiedfactor == NULL);
   assert(factor != NULL);
   assert(changed != NULL);

   /* enforces SP7 */
   if( SCIPisExprValue(scip, factor) )
   {
      *changed = TRUE;
      *simplifiedcoef *= SCIPgetValueExprValue(factor);
      return SCIP_OKAY;
   }

   /* enforces SP2 */
   if( SCIPisExprProduct(scip, factor) )
   {
      *changed = TRUE;

      /* assert SP8 */
      assert(SCIPgetCoefExprProduct(factor) == 1.0);
      debugSimplify("[simplifyFactor] seeing a product: include its children\n");

      SCIP_CALL( createExprlistFromExprs(scip, SCIPexprGetChildren(factor), SCIPexprGetNChildren(factor), simplifiedfactor) );

      return SCIP_OKAY;
   }

   /* enforces SP13: a sum with a unique child and no constant -> take the coefficient and use its child as factor */
   if( SCIPisExprSum(scip, factor) && SCIPexprGetNChildren(factor) == 1 && SCIPgetConstantExprSum(factor) == 0.0 )
   {
      *changed = TRUE;

      /* assert SS8 and SS7 */
      assert(SCIPgetCoefsExprSum(factor)[0] != 0.0 && SCIPgetCoefsExprSum(factor)[0] != 1.0);
      debugSimplify("[simplifyFactor] seeing a sum of the form coef * child : take coef and child apart\n");

      if( SCIPisExprProduct(scip, SCIPexprGetChildren(factor)[0]) )
      {
         /* if child is a product, then add its children to exprlist */
         SCIP_CALL( createExprlistFromExprs(scip, SCIPexprGetChildren(SCIPexprGetChildren(factor)[0]), SCIPexprGetNChildren(SCIPexprGetChildren(factor)[0]), simplifiedfactor) );
         *simplifiedcoef *= SCIPgetCoefExprProduct(SCIPexprGetChildren(factor)[0]);
      }
      else
      {
         SCIP_CALL( createExprlistFromExprs(scip, SCIPexprGetChildren(factor), 1, simplifiedfactor) );
      }
      *simplifiedcoef *= SCIPgetCoefsExprSum(factor)[0];

      return SCIP_OKAY;
   }

   /* the given (simplified) expression `factor`, can be a child of a simplified product */
   assert(!SCIPisExprProduct(scip, factor));
   assert(!SCIPisExprValue(scip, factor));

   SCIP_CALL( createExprNode(scip, factor, simplifiedfactor) );

   return SCIP_OKAY;
}

/** merges tomerge into finalchildren
 *
 * Both, tomerge and finalchildren contain expressions that could be the children of a simplified product
 * (except for SP8 and SP10 which are enforced later).
 * However, the concatenation of both lists will not in general yield a simplified product expression,
 * because SP4, SP5 and SP14 could be violated.  So the purpose of this method is to enforce SP4, SP5 and SP14.
 * In the process of enforcing SP4, it could happen that SP2 is violated. Since enforcing SP2
 * could generate further violations, we remove the affected children from finalchildren
 * and include them in unsimplifiedchildren for further processing.
 * @note if tomerge has more than one element, then they are the children of a simplified product expression
 */
static
SCIP_RETCODE mergeProductExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE*             tomerge,            /**< list to merge */
   EXPRNODE**            finalchildren,      /**< pointer to store the result of merge between tomerge and *finalchildren */
   EXPRNODE**            unsimplifiedchildren,/**< the list of children that should go to the product expression;
                                               *   they are unsimplified when seen as children of a simplified product */
   SCIP_Bool*            changed,            /**< pointer to store if some term actually got simplified */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   EXPRNODE* tomergenode;
   EXPRNODE* current;
   EXPRNODE* previous;

   if( tomerge == NULL )
      return SCIP_OKAY;

   if( *finalchildren == NULL )
   {
      *finalchildren = tomerge;
      return SCIP_OKAY;
   }

   tomergenode = tomerge;
   current = *finalchildren;
   previous = NULL;

   while( tomergenode != NULL && current != NULL )
   {
      int compareres;
      EXPRNODE* aux;
      SCIP_EXPR* base1;
      SCIP_EXPR* base2;
      SCIP_Real expo1;
      SCIP_Real expo2;
      SCIP_Bool issignpower1;
      SCIP_Bool issignpower2;

      /* assert invariants */
      assert(!SCIPisExprValue(scip, tomergenode->expr));
      assert(!SCIPisExprValue(scip, current->expr));
      assert(previous == NULL || previous->next == current);

      /* we are going to multiply the two exprs: current and tomergenode
       * we first check if they are both exponentials
       * if so, we multiply them
       * otherwise, we interpret them as base^exponent
       * the base of an expr is the expr itself
       * if type(expr) != pow, otherwise it is the child of pow
       */

      /* if both are exponentials, create a new exponential with the sum of their children */
      if( SCIPisExprExp(scip, current->expr) && SCIPisExprExp(scip, tomergenode->expr) )
      {
         SCIP_EXPR* sum;
         SCIP_EXPR* simplifiedsum;
         SCIP_EXPR* expexpr;
         SCIP_EXPR* simplifiedexp;

         /* inform that expressions changed */
         *changed = TRUE;

         /* create sum */
         SCIP_CALL( SCIPcreateExprSum(scip, &sum, 1, SCIPexprGetChildren(current->expr), NULL, 0.0, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPappendExprSumExpr(scip, sum, SCIPexprGetChildren(tomergenode->expr)[0], 1.0) );

         /* simplify sum */
         SCIP_CALL( SCIPcallExprSimplify(scip, sum, &simplifiedsum, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &sum) );

         /* create exponential */
         SCIP_CALL( SCIPcreateExprExp(scip, &expexpr, simplifiedsum, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedsum) );

         /* simplify exponential */
         SCIP_CALL( SCIPcallExprSimplify(scip, expexpr, &simplifiedexp, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expexpr) );

         /* note that simplified exponential might be a product exp(x) * exp(-x + log(y*z)) -> y*z and so it is not a
          * valid child of a simplified product; therefore we add it to the unsimplifiedchildren's list
          */

         /* replace tomergenode's expression with simplifiedexp */
         /* TODO: this code repeats below; add new function to avoid duplication */
         SCIP_CALL( SCIPreleaseExpr(scip, &tomergenode->expr) );
         tomergenode->expr = simplifiedexp;

         /* move tomergenode to unsimplifiedchildren */
         aux = tomergenode;
         tomergenode = tomergenode->next;
         insertFirstList(aux, unsimplifiedchildren);

         /* remove current */
         if( current == *finalchildren )
         {
            assert(previous == NULL);
            aux = listPopFirst(finalchildren);
            assert(aux == current);
            current = *finalchildren;
         }
         else
         {
            assert(previous != NULL);
            aux = current;
            current = current->next;
            previous->next = current;
         }
         SCIP_CALL( freeExprNode(scip, &aux) );

         continue;
      }

      /* they were not exponentials, so collect bases and exponents */
      if( SCIPisExprPower(scip, current->expr) )
      {
         base1 = SCIPexprGetChildren(current->expr)[0];
         expo1 = SCIPgetExponentExprPow(current->expr);
         issignpower1 = FALSE;
      }
      else if( SCIPisExprSignpower(scip, current->expr) )
      {
         base1 = SCIPexprGetChildren(current->expr)[0];
         expo1 = SCIPgetExponentExprPow(current->expr);
         issignpower1 = TRUE;
      }
      else
      {
         base1 = current->expr;
         expo1 = 1.0;
         issignpower1 = FALSE;
      }
      if( SCIPisExprPower(scip, tomergenode->expr) )
      {
         base2 = SCIPexprGetChildren(tomergenode->expr)[0];
         expo2 = SCIPgetExponentExprPow(tomergenode->expr);
         issignpower2 = FALSE;
      }
      else if( SCIPisExprSignpower(scip, tomergenode->expr) )
      {
         base2 = SCIPexprGetChildren(tomergenode->expr)[0];
         expo2 = SCIPgetExponentExprPow(tomergenode->expr);
         issignpower2 = TRUE;
      }
      else
      {
         base2 = tomergenode->expr;
         expo2 = 1.0;
         issignpower2 = FALSE;
      }

      if( SCIPcompareExpr(scip, base1, base2) == 0 )
      {
         /* the bases are the same, so we should try to merge the multiplication of the powers */
         SCIP_EXPR* power = NULL;

         if( !issignpower1 && !issignpower2 )
         {
            /* and both are normal power, then add to unsimplifiedchildren the resulting expr of simplify(base^(expo1 + expo2)) */
#if 0  /* TODO we should not loose the implicit base >= 0 constraint, if there is one, but then we should look at bounds on base; simplify currently doesn't */
            /*
             * unless expo1 or expo2 are fractional but expo1+expo2 is not fractional, then we better keep the original
             * the reason for that is that x^fractional implies a constraint x >= 0
             */
            if( (EPSISINT(expo1, 0.0) && EPSISINT(expo2, 0.0)) || !EPSISINT(expo1+expo2, 0.0) )  /*lint !e835*/
#endif
            {
               SCIP_CALL( SCIPcreateExprPow(scip, &power, base1, expo1 + expo2, ownercreate, ownercreatedata) );
            }
         }
         else if( issignpower1 ^ issignpower2 )
         {
            /* exactly one is signpower: sign(x) |x|^expo1 x^expo2 = sign(x)^(1+expo2) |x|^(expo1+expo2),  with x = base */
            if( EPSISINT(expo2, 0.0) )  /*lint !e835*/
            {
               if( (int)expo2 % 2 == 0 )
               {
                  /* if expo2 is even, then sign(x)^(1+expo2) = sign(x), so we have signpower: sign(x) |x|^(expo1+expo2)
                   * TODO: we can remove this case distinction once the simplification of power expressions tranform
                   * |expr|^even -> expr^even, since the call to SCIPcallExprSimplify(scip, conshdlr, power,
                   * &simplifiedpower) below will take care of this.
                   */
                  SCIP_CALL( SCIPcreateExprSignpower(scip, &power, base1, expo1 + expo2, ownercreate, ownercreatedata) );
               }
               else
               {
                  /* if expo2 is odd, then sign(x)^(1+expo2) = 1, so we have |x|^(expo1+expo2) */
                  SCIP_EXPR* absbase;

                  SCIP_CALL( SCIPcreateExprAbs(scip, &absbase, base1, ownercreate, ownercreatedata) );
                  SCIP_CALL( SCIPcreateExprPow(scip, &power, absbase, expo1 + expo2, ownercreate, ownercreatedata) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &absbase) );
               }
            }
            else if( !EPSISINT(expo1+expo2, 0.0) )  /*lint !e835*/
            {
               /* if expo2 is fractional and expo1+expo2 is fractional, then we need x >= 0, so we can use x^(expo1+expo2) */
               SCIP_CALL( SCIPcreateExprPow(scip, &power, base1, expo1 + expo2, ownercreate, ownercreatedata) );
            }
            /* else: expo2 is fractional but expo1+expo2 is integral, then we better do not do anything for now
             * (leave power at NULL)
             */
         }
         else
         {
            /* if both are signpower, then we have |base|^(expo1+expo2)
             * if expo1+expo2 is even, then we can change this to base^(expo1+expo2)
             */
            if( EPSISINT(expo1+expo2, 0.0) && (int)(expo1+expo2)%2 == 0 ) /*lint !e835*/
            {
               SCIP_CALL( SCIPcreateExprPow(scip, &power, base1, expo1 + expo2, ownercreate, ownercreatedata) );
            }
            else
            {
               SCIP_EXPR* absbase;

               SCIP_CALL( SCIPcreateExprAbs(scip, &absbase, base1, ownercreate, ownercreatedata) );
               SCIP_CALL( SCIPcreateExprPow(scip, &power, absbase, expo1 + expo2, ownercreate, ownercreatedata) );
               SCIP_CALL( SCIPreleaseExpr(scip, &absbase) );
            }
         }

         if( power != NULL )
         {
            /* we have created a new power: simplify again and continue */
            SCIP_EXPR* simplifiedpower;

            /* call simplifyPow or simplifySignpower */
            SCIP_CALL( SCIPcallExprSimplify(scip, power, &simplifiedpower, ownercreate, ownercreatedata) );
            SCIP_CALL( SCIPreleaseExpr(scip, &power) );

            /* replace tomergenode's expression with simplifiedpower */
            SCIP_CALL( SCIPreleaseExpr(scip, &tomergenode->expr) );
            tomergenode->expr = simplifiedpower;

            *changed = TRUE;

            /* move tomergenode to unsimplifiedchildren */
            aux = tomergenode;
            tomergenode = tomergenode->next;
            insertFirstList(aux, unsimplifiedchildren);

            /* remove current */
            if( current == *finalchildren )
            {
               assert(previous == NULL);
               aux = listPopFirst(finalchildren);
               assert(aux == current);
               current = *finalchildren;
            }
            else
            {
               assert(previous != NULL);
               aux = current;
               current = current->next;
               previous->next = current;
            }
            SCIP_CALL( freeExprNode(scip, &aux) );

            continue;
         }
      }

      /* bases are not the same, or we do not want to merge them
       * then expressions cannot be the same
       * therefore we need to insert tomergenode in finalchildren
       * for this, we need to take care of the order
       */
      compareres = SCIPcompareExpr(scip, current->expr, tomergenode->expr);
      if( compareres == -1 )
      {
         /* current < tomergenode => move current */
         previous = current;
         current = current->next;
      }
      else
      {
         *changed = TRUE;
         assert(compareres == 1);

         /* insert: if current is the first node, then insert at beginning; otherwise, insert between previous and current */
         if( current == *finalchildren )
         {
            assert(previous == NULL);
            aux = tomergenode;
            tomergenode = tomergenode->next;
            insertFirstList(aux, finalchildren);
            previous = *finalchildren;
         }
         else
         {
            assert(previous != NULL);
            /* extract */
            aux = tomergenode;
            tomergenode = tomergenode->next;
            /* insert */
            previous->next = aux;
            aux->next = current;
            previous = aux;
         }
      }
   }

   /* if all nodes of tomerge were merged, we are done */
   if( tomergenode == NULL )
      return SCIP_OKAY;

   assert(current == NULL);

   /* if all nodes of finalchildren were cancelled by nodes of tomerge (i.e., transfered to unsimplifiedchildren),
    * then the rest of tomerge is finalchildren
    */
   if( *finalchildren == NULL )
   {
      assert(previous == NULL);
      *finalchildren = tomergenode;
      return SCIP_OKAY;
   }

   /* there are still nodes of tomerge unmerged; these nodes are larger than finalchildren, so append at end */
   assert(previous != NULL && previous->next == NULL);
   previous->next = tomergenode;

   return SCIP_OKAY;
}

/** simplifies the given (simplified) exprs so that they can be factors of a simplified product
 *
 * in particular, it will sort and multiply factors whose product leads to new expressions
 */
static
SCIP_RETCODE simplifyMultiplyChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           exprs,              /**< factors to be simplified */
   int                   nexprs,             /**< number of factors */
   SCIP_Real*            simplifiedcoef,     /**< buffer to store coefficient of PI exprs; needs to be initialized */
   EXPRNODE**            finalchildren,      /**< expr node list to store the simplified factors */
   SCIP_Bool*            changed,            /**< buffer to store whether some factor changed */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   EXPRNODE* unsimplifiedchildren;

   /* set up list of current children (when looking at each of them individually, they are simplified, but as
    * children of a product expression they might be unsimplified)
    */
   unsimplifiedchildren = NULL;
   SCIP_CALL( createExprlistFromExprs(scip, exprs, nexprs, &unsimplifiedchildren) );

   *changed = FALSE;

   /* while there are still children to process */
   *finalchildren  = NULL;
   while( unsimplifiedchildren != NULL )
   {
      EXPRNODE* tomerge;
      EXPRNODE* first;

      first = listPopFirst(&unsimplifiedchildren);
      assert(first != NULL);

#ifdef SIMPLIFY_DEBUG
      debugSimplify("simplifying factor:\n");
      SCIP_CALL( SCIPprintExpr(scip, first->expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* enforces SP2, SP7 and SP13 */
      tomerge = NULL;
      SCIP_CALL( simplifyFactor(scip, first->expr, simplifiedcoef, &tomerge, changed) );

      /* enforces SP4 and SP5 note: merge frees (or uses) the nodes of the tomerge list */
      SCIP_CALL( mergeProductExprlist(scip, tomerge, finalchildren, &unsimplifiedchildren, changed, ownercreate, ownercreatedata) );

      /* free first */
      SCIP_CALL( freeExprlist(scip, &first) );

      /* if the simplified coefficient is 0, we can return value 0 */
      if( *simplifiedcoef == 0.0 )
      {
         *changed = TRUE;
         SCIP_CALL( freeExprlist(scip, finalchildren) );
         SCIP_CALL( freeExprlist(scip, &unsimplifiedchildren) );
         assert(*finalchildren == NULL);
         break;
      }
   }
   return SCIP_OKAY;
}

/* make sure product has at least two children
 * - if it is empty; return value
 * - if it has one child and coef = 1; return child
 * - if it has one child and coef != 1; return (sum 0 coef expr)
 */
static
SCIP_RETCODE enforceSP10(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             simplifiedcoef,     /**< simplified product should be simplifiedcoef * PI simplifiedfactors */
   EXPRNODE*             finalchildren,      /**< factors of simplified product */
   SCIP_EXPR**           simplifiedexpr,     /**< buffer to store the simplified expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   /* empty list --> return value */
   if( finalchildren == NULL )
   {
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, simplifiedcoef, ownercreate, ownercreatedata) );
      return SCIP_OKAY;
   }

   /* one child and coef equal to 1 --> return child */
   if( finalchildren->next == NULL && simplifiedcoef == 1.0 )
   {
      *simplifiedexpr = finalchildren->expr;
      SCIPcaptureExpr(*simplifiedexpr);
      return SCIP_OKAY;
   }

   /* one child and coef different from 1 --> return (sum 0 coef child) */
   if( finalchildren->next == NULL )
   {
      SCIP_EXPR* sum;

      SCIP_CALL( SCIPcreateExprSum(scip, &sum, 1, &(finalchildren->expr), &simplifiedcoef, 0.0, ownercreate, ownercreatedata) );

      /* simplifying here is necessary, the product could have sums as children e.g., (prod 2 (sum 1 <x>))
       * -> (sum 0 2 (sum 1 <x>)) and that needs to be simplified to (sum 0 2 <x>)
       */
      SCIP_CALL( SCIPcallExprSimplify(scip, sum, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &sum) );
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** checks if it is entropy expression */
static
SCIP_RETCODE enforceSP11(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             simplifiedcoef,     /**< simplified product should be simplifiedcoef * PI simplifiedfactors */
   EXPRNODE*             finalchildren,      /**< factors of simplified product */
   SCIP_EXPR**           simplifiedexpr,     /**< buffer to store the simplified expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPR* entropicchild = NULL;

   if( !(finalchildren != NULL && finalchildren->next != NULL && finalchildren->next->next == NULL) )
      return SCIP_OKAY;

   /* could be log(expr) * expr, e.g., log(sin(x)) * sin(x) (OR11) */
   if( strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(finalchildren->expr)), "log") == 0 )
   {
      assert(SCIPexprGetNChildren(finalchildren->expr) == 1);
      if( 0 == SCIPcompareExpr(scip, SCIPexprGetChildren(finalchildren->expr)[0], finalchildren->next->expr) )
         entropicchild = finalchildren->next->expr;
   }
   /* could be expr * log(expr), e.g., (1 + abs(x)) log(1 + abs(x)) (OR11) */
   else if( strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(finalchildren->next->expr)), "log") == 0 )
   {
      assert(SCIPexprGetNChildren(finalchildren->next->expr) == 1);
      if( 0 == SCIPcompareExpr(scip, SCIPexprGetChildren(finalchildren->next->expr)[0], finalchildren->expr) )
         entropicchild = finalchildren->expr;
   }

   /* success --> replace finalchildren by entropy expression */
   if( entropicchild != NULL )
   {
      SCIP_EXPR* entropy;

      simplifiedcoef *= -1.0;

      SCIP_CALL( SCIPcreateExprEntropy(scip, &entropy, entropicchild, ownercreate, ownercreatedata) );

      /* enforces SP8: if simplifiedcoef != 1.0, transform it into a sum with the (simplified) entropy as child */
      if( simplifiedcoef != 1.0 )
      {
         SCIP_CALL( SCIPcreateExprSum(scip, simplifiedexpr, 1, &entropy, &simplifiedcoef, 0.0, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &entropy) );
      }
      else
         *simplifiedexpr = entropy;
   }

   return SCIP_OKAY;
}

/* expands product of two sums or one sum and another expression
 * -) two sums: (prod (sum c1 s1 ... sn) (sum c2 t1 ... tm)
 *    Builds a sum representing the expansion, where all of its children are simplified, and then simplify the sum
 *    - constant != 0 --> c1 ti or c2 * sj is simplified (ti, sj are not sums, because they are children of a simplified sum)
 *    - sj * ti may be not be simplified, so put them in a product list and simplify them from there
 * -) one sum: (prod factor (sum c s1 ... sn))
 *    - c != 0 --> c * factor is simplified (i.e. factor is not sum!)
 *    - factor * si may be not be simplified, so put them in a product list and simplify them from there
 */
static
SCIP_RETCODE enforceSP12(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             simplifiedcoef,     /**< simplified product should be simplifiedcoef * PI simplifiedfactors */
   EXPRNODE*             finalchildren,      /**< factors of simplified product */
   SCIP_EXPR**           simplifiedexpr,     /**< buffer to store the simplified expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   /* we need only two children */
   if( ! (finalchildren != NULL && finalchildren->next != NULL && finalchildren->next->next == NULL) )
      return SCIP_OKAY;

   /* handle both sums case */
   if( SCIPisExprSum(scip, finalchildren->expr) && SCIPisExprSum(scip, finalchildren->next->expr) )
   {
      SCIP_EXPR* expanded = NULL;
      SCIP_Real c1 = SCIPgetConstantExprSum(finalchildren->expr);
      SCIP_Real c2 = SCIPgetConstantExprSum(finalchildren->next->expr);
      int nchildren1 = SCIPexprGetNChildren(finalchildren->expr);
      int nchildren2 = SCIPexprGetNChildren(finalchildren->next->expr);
      int j;
      int k;

#ifdef SIMPLIFY_DEBUG
      debugSimplify("Multiplying sum1 * sum2\n");
      debugSimplify("sum1: \n");
      SCIP_CALL( SCIPprintExpr(scip, finalchildren->expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
      debugSimplify("sum2: \n");
      SCIP_CALL( SCIPprintExpr(scip, finalchildren->next->expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif
      SCIP_CALL( SCIPcreateExprSum(scip, &expanded, 0, NULL, NULL, c1 * c2 * simplifiedcoef, ownercreate, ownercreatedata) );

      /* multiply c1 * sum2 */
      if( c1 != 0.0 )
      {
         int i;

         for( i = 0; i < nchildren2; ++i )
         {
            SCIP_EXPR* term;

            term = SCIPexprGetChildren(finalchildren->next->expr)[i];
            SCIP_CALL( SCIPappendExprSumExpr(scip, expanded, term, SCIPgetCoefsExprSum(finalchildren->next->expr)[i] * c1 * simplifiedcoef) );
            /* we are just re-using a child here, so do not release term! */
#ifdef SIMPLIFY_DEBUG
            debugSimplify("Multiplying %f * summand2_i\n", c1);
            debugSimplify("summand2_i: \n");
            SCIP_CALL( SCIPprintExpr(scip, term, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
         }
      }
      /* multiply c2 * sum1 */
      if( c2 != 0.0 )
      {
         int i;

         for( i = 0; i < nchildren1; ++i )
         {
            SCIP_EXPR* term;

            term = SCIPexprGetChildren(finalchildren->expr)[i];
            SCIP_CALL( SCIPappendExprSumExpr(scip, expanded, term, SCIPgetCoefsExprSum(finalchildren->expr)[i] * c2 * simplifiedcoef) );
            /* we are just re-using a child here, so do not release term! */
#ifdef SIMPLIFY_DEBUG
            debugSimplify("Multiplying summand1_i * %f\n", c2);
            debugSimplify("summand1_i: \n");
            SCIP_CALL( SCIPprintExpr(scip, term, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
         }
      }
      /* multiply sum1 * sum2 without constants */
      for( j = 0; j < nchildren1; ++j )
      {
         SCIP_EXPR* factors[2];
         SCIP_Real coef1;

         coef1 = SCIPgetCoefsExprSum(finalchildren->expr)[j];
         factors[0] = SCIPexprGetChildren(finalchildren->expr)[j];
         for( k = 0; k < nchildren2; ++k )
         {
            EXPRNODE* finalfactors;
            SCIP_Real factorscoef;
            SCIP_Real coef2;
            SCIP_EXPR* term = NULL;
            SCIP_Bool dummy;

            coef2 = SCIPgetCoefsExprSum(finalchildren->next->expr)[k];
            factors[1] = SCIPexprGetChildren(finalchildren->next->expr)[k];

#ifdef SIMPLIFY_DEBUG
            debugSimplify("multiplying %g expr1 * %g expr2\n", coef1, coef2);
            debugSimplify("expr1:\n");
            SCIP_CALL( SCIPprintExpr(scip, factors[0], NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
            debugSimplify("expr2\n");
            SCIP_CALL( SCIPprintExpr(scip, factors[1], NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif

            factorscoef = coef1 * coef2;
            SCIP_CALL( simplifyMultiplyChildren(scip, factors, 2, &factorscoef, &finalfactors, &dummy, ownercreate, ownercreatedata) );
            assert(factorscoef != 0.0);

#ifdef SIMPLIFY_DEBUG
            {
               EXPRNODE* node;
               int i;

               debugSimplify("Building product from simplified factors\n");
               node = finalfactors;
               i = 0;
               while( node != NULL )
               {
                  debugSimplify("factor %d (nuses %d):\n", i, SCIPexprGetNUses(node->expr));
                  SCIP_CALL( SCIPprintExpr(scip, node->expr, NULL) );
                  SCIPinfoMessage(scip, NULL, "\n");
                  node = node->next;
                  i++;
               }
            }
#endif

            SCIP_CALL( buildSimplifiedProduct(scip, 1.0, &finalfactors, TRUE, &term, ownercreate, ownercreatedata) );
            assert(finalfactors == NULL);
            assert(term != NULL);

#ifdef SIMPLIFY_DEBUG
            debugSimplify("%g expr1 * %g expr2 = %g * product\n", coef1, coef2, coef1 * coef2);
            debugSimplify("product: (nused %d)\n", SCIPexprGetNUses(term));
            SCIP_CALL( SCIPprintExpr(scip, term, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif

            SCIP_CALL( SCIPappendExprSumExpr(scip, expanded, term, factorscoef * simplifiedcoef) );

            SCIP_CALL( SCIPreleaseExpr(scip, &term) );
         }
      }

      /* simplify the sum */
      SCIP_CALL( SCIPcallExprSimplify(scip, expanded, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &expanded) );

      return SCIP_OKAY;
   }

   /* handle one sum case */
   if( SCIPisExprSum(scip, finalchildren->expr) || SCIPisExprSum(scip, finalchildren->next->expr) )
   {
      SCIP_EXPR* expanded = NULL;
      SCIP_EXPR* factors[2];
      SCIP_EXPR* sum = NULL;
      SCIP_Real constant;
      int nchildren;
      int j;

      if( SCIPisExprSum(scip, finalchildren->expr) )
      {
         assert(!SCIPisExprSum(scip, finalchildren->next->expr));
         sum = finalchildren->expr;
         factors[0] = finalchildren->next->expr;
      }
      else
      {
         assert(!SCIPisExprSum(scip, finalchildren->expr));
         sum = finalchildren->next->expr;
         factors[0] = finalchildren->expr;
      }
      constant = simplifiedcoef * SCIPgetConstantExprSum(sum);
      nchildren = SCIPexprGetNChildren(sum);

      SCIP_CALL( SCIPcreateExprSum(scip, &expanded, 1, &factors[0], &constant, 0.0, ownercreate, ownercreatedata) );
      /* we are just re-using a child here, so do not release factor! */

      for( j = 0; j < nchildren; ++j )
      {
         SCIP_Real coef;
         SCIP_Real termcoef;
         SCIP_Bool dummy;
         EXPRNODE* finalfactors;
         SCIP_EXPR* term = NULL;

         coef = SCIPgetCoefsExprSum(sum)[j];
         factors[1] = SCIPexprGetChildren(sum)[j];

         termcoef = coef;
         SCIP_CALL( simplifyMultiplyChildren(scip, factors, 2, &termcoef, &finalfactors, &dummy, ownercreate, ownercreatedata) );
         assert(termcoef != 0.0);

         SCIP_CALL( buildSimplifiedProduct(scip, 1.0, &finalfactors, TRUE, &term, ownercreate, ownercreatedata) );
         assert(finalfactors == NULL);
         assert(term != NULL);

         SCIP_CALL( SCIPappendExprSumExpr(scip, expanded, term, termcoef * simplifiedcoef) );
         SCIP_CALL( SCIPreleaseExpr(scip, &term) );
      }

      /* simplify the sum */
      SCIP_CALL( SCIPcallExprSimplify(scip, expanded, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &expanded) );
   }

   return SCIP_OKAY;
}

/** builds a simplified product from simplifiedfactors
 *
 * @note this function also releases simplifiedfactors
 */
static
SCIP_RETCODE buildSimplifiedProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             simplifiedcoef,     /**< simplified product should be simplifiedcoef * PI simplifiedfactors */
   EXPRNODE**            simplifiedfactors,  /**< factors of simplified product */
   SCIP_Bool             changed,            /**< indicates whether some of the simplified factors was changed */
   SCIP_EXPR**           simplifiedexpr,     /**< buffer to store the simplified expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   EXPRNODE* finalchildren = *simplifiedfactors;

   /* build product expression from finalchildren and post-simplify */
   debugSimplify("[simplifyProduct] finalchildren has length %d\n", listLength(finalchildren));

   *simplifiedexpr = NULL;

   SCIP_CALL( enforceSP11(scip, simplifiedcoef, *simplifiedfactors, simplifiedexpr, ownercreate, ownercreatedata) );
   if( *simplifiedexpr != NULL )
      goto CLEANUP;

   SCIP_CALL( enforceSP12(scip, simplifiedcoef, *simplifiedfactors, simplifiedexpr, ownercreate, ownercreatedata) );
   if( *simplifiedexpr != NULL )
      goto CLEANUP;

   SCIP_CALL( enforceSP10(scip, simplifiedcoef, *simplifiedfactors, simplifiedexpr, ownercreate, ownercreatedata) );
   if( *simplifiedexpr != NULL )
      goto CLEANUP;

   /* enforces SP8: if simplifiedcoef != 1.0, transform it into a sum with the (simplified) product as child */
   if( simplifiedcoef != 1.0 )
   {
      SCIP_EXPR* aux;
      SCIP_EXPR* sum;

      /* create sum */
      SCIP_CALL( createExprProductFromExprlist(scip, finalchildren, 1.0, &aux, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPcreateExprSum(scip, &sum, 1, &aux, &simplifiedcoef, 0.0, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &aux) );

      /* simplify sum */
      SCIP_CALL( SCIPcallExprSimplify(scip, sum, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &sum) );

      goto CLEANUP;
   }

   /* build product expression from list */
   if( changed )
   {
      SCIP_CALL( createExprProductFromExprlist(scip, finalchildren, simplifiedcoef, simplifiedexpr, ownercreate, ownercreatedata) );
      goto CLEANUP;
   }

CLEANUP:

   SCIP_CALL( freeExprlist(scip, simplifiedfactors) );
   return SCIP_OKAY;
}

/** computes an estimator for a product as a vertex polyhedral function
 *
 * Since the product is multilinear, its convex and concave envelopes are piecewise linear.
 */
static
SCIP_RETCODE estimateVertexPolyhedralProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   int                   nfactors,           /**< number of factors */
   SCIP_INTERVAL*        bounds,             /**< bound for each factor */
   SCIP_Real             constantfactor,     /**< another constant factor */
   SCIP_Real*            refpoint,           /**< reference point where to estimate, or NULL if called from initestimates */
   SCIP_Bool             overestimate,       /**< should estimator overestimate expr (TRUE) or underestimate (FALSE) */
   SCIP_Real             targetvalue,        /**< no need to compute facet if value in xstar would be worse than target value */
   SCIP_Real*            coefs,              /**< array to store cut coefficients */
   SCIP_Real*            constant,           /**< pointer to store cut constant */
   SCIP_Bool*            success             /**< pointer to store whether estimation was successful */
   )
{
   SCIP_Real* box;
   SCIP_Real* xstar;
   int nfixed;
   int i;

   assert(conshdlr != NULL);
   assert(nfactors > 0);
   assert(bounds != NULL);
   assert(constantfactor != 0.0);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* assemble box, check for unbounded variables, assemble xstar */
   SCIP_CALL( SCIPallocBufferArray(scip, &box, 2*nfactors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xstar, nfactors) );
   for( i = 0, nfixed = 0; i < nfactors; ++i )
   {
      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bounds[i]));

      if( SCIPisInfinity(scip, -bounds[i].inf) || SCIPisInfinity(scip, bounds[i].sup) )
      {
         SCIPdebugMsg(scip, "a factor is unbounded, no cut is possible\n");
         goto CLEANUP;
      }

      box[2*i] = bounds[i].inf;
      box[2*i+1] = bounds[i].sup;

      xstar[i] = refpoint != NULL ? refpoint[i] : 0.5 * (box[2*i] + box[2*i+1]);

      if( SCIPisRelEQ(scip, box[2*i], box[2*i+1]) )
         ++nfixed;
   }

   if( nfixed < nfactors && nfactors - nfixed <= SCIP_MAXVERTEXPOLYDIM )
   {
      SCIP_CALL( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, conshdlr,
         overestimate, prodfunction, &constantfactor, xstar, box, nfactors, targetvalue, success, coefs, constant) );
   }

CLEANUP:
   SCIPfreeBufferArray(scip, &xstar);
   SCIPfreeBufferArray(scip, &box);

   return SCIP_OKAY;
}

/*
 * Callback methods of expression handler
 */

/** simplifies a product expression
 *
 * Summary: we first build a list of expressions (called finalchildren) which will be the children of the simplified product
 * and then we process this list in order to enforce SP8 and SP10.
 *
 * Description: In order to build finalchildren, we first build a list of unsimplified children (called unsimplifiedchildren)
 * with the children of the product. Each node of the list is manipulated (see simplifyFactor) in order to satisfy
 * SP2 and SP7 as follows:
 * - SP7: if the node's expression is a value, multiply the value to the products's coef
 * - SP2: if the node's expression is a product, then build a list with the child's children
 *
 * Then, we merge the built list (or the simplified node) into finalchildren. While merging, nodes from finalchildren
 * can go back to unsimplifiedchildren for further processing (see mergeProductExprlist() for more details).
 * After building finalchildren, we create the simplified product out of it, taking care that SP8 and SP10 are satisfied
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyProduct)
{  /*lint --e{715}*/
   EXPRNODE* finalchildren;
   SCIP_Real simplifiedcoef;
   SCIP_Bool changed;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);

   simplifiedcoef = SCIPgetCoefExprProduct(expr);

#ifdef SIMPLIFY_DEBUG
   debugSimplify("Simplifying expr:\n");
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   debugSimplify("First multiplying children\n");
#endif

   /* simplify and multiply factors */
   SCIP_CALL( simplifyMultiplyChildren(scip, SCIPexprGetChildren(expr), SCIPexprGetNChildren(expr), &simplifiedcoef,
         &finalchildren, &changed, ownercreate, ownercreatedata) );

#ifdef SIMPLIFY_DEBUG
   {
      EXPRNODE* node;
      int i;

      debugSimplify("Building product from simplified factors\n");
      node = finalchildren;
      i = 0;
      while( node != NULL )
      {
         debugSimplify("factor %d:\n", i);
         SCIP_CALL( SCIPprintExpr(scip, node->expr, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
         node = node->next;
         i++;
      }
   }
#endif

   /* get simplified product from simplified factors in finalchildren */
   SCIP_CALL( buildSimplifiedProduct(scip, simplifiedcoef, &finalchildren, changed, simplifiedexpr, ownercreate,
         ownercreatedata) );
   assert(finalchildren == NULL);

   if( *simplifiedexpr == NULL )
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
   }
   assert(*simplifiedexpr != NULL);

   return SCIP_OKAY;
}

/** compare two product expressions
 *
 *  The order of two product expressions, u and v, is a lexicographical order on the factors.
 *
 *  Starting from the *last*, we find the first child where they differ, say, the i-th.
 *  Then u < v <=> u_i < v_i.
 *  If there is no such children and they have different number of children, then u < v <=> nchildren(u) < nchildren(v).
 *  If all children are the same and they have the same number of children, then u < v <=> coeff(u) < coeff(v).
 *  Otherwise, they are the same.
 *
 *  Note: we are assuming expression are simplified, so within u, we have u_1 < u_2, etc.
 *
 *  Example: y * z < x * y * z
 */
static
SCIP_DECL_EXPRCOMPARE(compareProduct)
{  /*lint --e{715}*/
   int compareresult;
   int i;
   int j;
   int nchildren1;
   int nchildren2;
   SCIP_EXPR** children1;
   SCIP_EXPR** children2;

   nchildren1 = SCIPexprGetNChildren(expr1);
   nchildren2 = SCIPexprGetNChildren(expr2);
   children1 = SCIPexprGetChildren(expr1);
   children2 = SCIPexprGetChildren(expr2);

   for( i = nchildren1 - 1, j = nchildren2 - 1; i >= 0 && j >= 0; --i, --j )
   {
      compareresult = SCIPcompareExpr(scip, children1[i], children2[j]);
      if( compareresult != 0 )
         return compareresult;
      /* expressions are equal, continue */
   }

   /* all children of one expression are children of the other expression, use number of children as a tie-breaker */
   if( i < j )
   {
      assert(i == -1);
      /* expr1 has less elements, hence expr1 < expr2 */
      return -1;
   }
   if( i > j )
   {
      assert(j == -1);
      /* expr1 has more elements, hence expr1 > expr2 */
      return 1;
   }

   /* everything is equal, use coefficient as tie-breaker */
   assert(i == -1 && j == -1);
   if( SCIPgetCoefExprProduct(expr1) < SCIPgetCoefExprProduct(expr2) )
      return -1;
   if( SCIPgetCoefExprProduct(expr1) > SCIPgetCoefExprProduct(expr2) )
      return 1;

   /* they are equal */
   return 0;
}

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrProduct)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrProduct(scip) );

   return SCIP_OKAY;
}

/** expression handler free callback */
static
SCIP_DECL_EXPRFREEHDLR(freehdlrProduct)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(exprhdlr != NULL);
   assert(exprhdlrdata != NULL);
   assert(*exprhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, exprhdlrdata);
   assert(*exprhdlrdata == NULL);

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* sourceexprdata;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   sourceexprdata = SCIPexprGetData(sourceexpr);
   assert(sourceexprdata != NULL);

   SCIP_CALL( SCIPduplicateBlockMemory(targetscip, targetexprdata, sourceexprdata) );

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   SCIPfreeBlockMemory(scip, &exprdata);

   SCIPexprSetData(expr, NULL);

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_EXPRPRINT(printProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   switch( stage )
   {
      case SCIP_EXPRITER_ENTEREXPR :
      {
         /* print opening parenthesis, if necessary */
         if( EXPRHDLR_PRECEDENCE <= parentprecedence )
         {
            SCIPinfoMessage(scip, file, "(");
         }

         /* print coefficient, if not one */
         if( exprdata->coefficient != 1.0 )
         {
            if( exprdata->coefficient < 0.0 && EXPRHDLR_PRECEDENCE > parentprecedence )
            {
               SCIPinfoMessage(scip, file, "(%g)", exprdata->coefficient);
            }
            else
            {
               SCIPinfoMessage(scip, file, "%g", exprdata->coefficient);
            }
         }
         break;
      }

      case SCIP_EXPRITER_VISITINGCHILD :
      {
         /* print multiplication sign, if not first factor */
         if( exprdata->coefficient != 1.0 || currentchild > 0 )
         {
            SCIPinfoMessage(scip, file, "*");
         }
         break;
      }

      case SCIP_EXPRITER_VISITEDCHILD :
      {
         break;
      }

      case SCIP_EXPRITER_LEAVEEXPR :
      {
         /* print closing parenthesis, if necessary */
         if( EXPRHDLR_PRECEDENCE <= parentprecedence )
         {
            SCIPinfoMessage(scip, file, ")");
         }
         break;
      }

      default:
         /* all stages should have been covered above */
         SCIPABORT();
   }

   return SCIP_OKAY;
}

/** product hash callback */
static
SCIP_DECL_EXPRHASH(hashProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(exprdata->coefficient);

   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      *hashkey ^= childrenhashes[c];

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   SCIP_Real childval;
   int c;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   *val = exprdata->coefficient;
   for( c = 0; c < SCIPexprGetNChildren(expr) && (*val != 0.0); ++c )
   {
      childval = SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[c]);
      assert(childval != SCIP_INVALID);

      *val *= childval;
   }

   return SCIP_OKAY;
}

/** derivative evaluation callback computing <gradient, children.dot>
 *
 * If expr is \f$\prod_i x_i\f$, then computes \f$\sum_j \prod_{i\neq j} x_i x^{\text{dot}}_j\f$.
 */
static
SCIP_DECL_EXPRFWDIFF(fwdiffProduct)
{  /*lint --e{715}*/
   int c;

   assert(expr != NULL);
   assert(dot != NULL);

   assert(SCIPexprGetData(expr) != NULL);

   /* TODO add special handling for nchildren == 2 */

   /**! [SnippetExprFwdiffProduct] */
   *dot = 0.0;
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      SCIP_EXPR* child;

      child = SCIPexprGetChildren(expr)[c];

      assert(SCIPexprGetEvalValue(child) != SCIP_INVALID);
      assert(SCIPexprGetDot(child) != SCIP_INVALID);

      if( SCIPexprGetDot(child) == 0.0 )
         continue;

      if( SCIPexprGetEvalValue(child) != 0.0 )
         *dot += SCIPexprGetEvalValue(expr) / SCIPexprGetEvalValue(child) * SCIPexprGetDot(child);
      else
      {
         SCIP_Real partial;
         int i;

         partial = SCIPexprGetData(expr)->coefficient;
         for( i = 0; i < SCIPexprGetNChildren(expr) && (partial != 0.0); ++i )
         {
            if( i == c )
               continue;

            partial *= SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[i]);
         }
         *dot += partial * SCIPexprGetDot(child);
      }
   }
   /**! [SnippetExprFwdiffProduct] */

   return SCIP_OKAY;
}

/** expression backward forward derivative evaluation callback
 *
 * Computes \f$\frac{\partial}{\partial \text{childidx}} ( \langle \text{gradient}, \text{children.dot}\rangle )\f$.
 *
 * If expr is \f$\prod_i x_i\f$, and childidx is \f$k\f$ then computes
 *   \f$\partial_k \sum_j \prod_{i \neq j} x_i x^{\text{dot}}_j
 *   = \sum_{j \neq k} \prod_{i \neq j, k} x_i x^{\text{dot}}_j\f$
 */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffProduct)
{  /*lint --e{715}*/
   SCIP_EXPR* partialchild;
   int c;

   assert(expr != NULL);
   assert(bardot != NULL);
   assert(SCIPexprGetData(expr) != NULL);
   assert(childidx >= 0 && childidx < SCIPexprGetNChildren(expr));

   partialchild = SCIPexprGetChildren(expr)[childidx];
   assert(partialchild != NULL);
   assert(!SCIPisExprValue(scip, partialchild));
   assert(SCIPexprGetEvalValue(partialchild) != SCIP_INVALID);

   /* TODO add special handling for nchildren == 2 */

   /**! [SnippetExprBwfwdiffProduct] */
   *bardot = 0.0;
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      SCIP_EXPR* child;

      if( c == childidx )
         continue;

      child = SCIPexprGetChildren(expr)[c];

      assert(SCIPexprGetEvalValue(child) != SCIP_INVALID);
      assert(SCIPexprGetDot(child) != SCIP_INVALID);

      if( SCIPexprGetDot(child) == 0.0 )
         continue;

      if( SCIPexprGetEvalValue(child) != 0.0 && SCIPexprGetEvalValue(partialchild) != 0.0 )
         *bardot += SCIPexprGetEvalValue(expr) / (SCIPexprGetEvalValue(child) * SCIPexprGetEvalValue(partialchild)) * SCIPexprGetDot(child);
      else
      {
         SCIP_Real partial;
         int i;

         partial = SCIPexprGetData(expr)->coefficient;
         for( i = 0; i < SCIPexprGetNChildren(expr) && (partial != 0.0); ++i )
         {
            if( i == c || i == childidx )
               continue;

            partial *= SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[i]);
         }
         *bardot += partial * SCIPexprGetDot(child);
      }
   }
   /**! [SnippetExprBwfwdiffProduct] */

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffProduct)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);
   assert(childidx >= 0 && childidx < SCIPexprGetNChildren(expr));

   child = SCIPexprGetChildren(expr)[childidx];
   assert(child != NULL);
   assert(!SCIPisExprValue(scip, child));
   assert(SCIPexprGetEvalValue(child) != SCIP_INVALID);

   /* TODO add special handling for nchildren == 2 */

   /**! [SnippetExprBwdiffProduct] */
   if( !SCIPisZero(scip, SCIPexprGetEvalValue(child)) )
   {
      *val = SCIPexprGetEvalValue(expr) / SCIPexprGetEvalValue(child);
   }
   else
   {
      int i;

      *val = SCIPexprGetData(expr)->coefficient;
      for( i = 0; i < SCIPexprGetNChildren(expr) && (*val != 0.0); ++i )
      {
         if( i == childidx )
            continue;

         *val *= SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[i]);
      }
   }
   /**! [SnippetExprBwdiffProduct] */

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int c;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   /**! [SnippetExprIntevalProduct] */
   SCIPintervalSet(interval, exprdata->coefficient);

   SCIPdebugMsg(scip, "inteval %p with %d children: %.20g", (void*)expr, SCIPexprGetNChildren(expr), exprdata->coefficient);

   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      SCIP_INTERVAL childinterval;

      childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[c]);
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      {
         SCIPintervalSetEmpty(interval);
         break;
      }

      /* multiply childinterval with the so far computed interval */
      SCIPintervalMul(SCIP_INTERVAL_INFINITY, interval, *interval, childinterval);

      SCIPdebugMsgPrint(scip, " *[%.20g,%.20g]", childinterval.inf, childinterval.sup);
   }
   SCIPdebugMsgPrint(scip, " = [%.20g,%.20g]\n", interval->inf, interval->sup);
   /**! [SnippetExprIntevalProduct] */

   return SCIP_OKAY;
}

/** estimates a multilinear function of the form \f$ f(x) := a \prod_{i = 1}^n x_i  \f$
 *
 * \f$ x_i \f$ are the auxiliary variables of the children.
 * If !overestimate, then we look for an affine underestimator of \f$ f(x) \f$ which has a value above targetvalue at \f$ x^* \f$,
 * i.e., \f$ g(x) := \alpha^T x + \beta \le f(x)\f$ for all \f$ x \f$ in the domain, such that \f$ \alpha x^* + \beta > \text{targetvalue}\f$.
 *
 * Since \f$ f(x) \f$ is componentwise linear, its convex envelope is piecewise linear and its value can be computed by
 * finding the largest affine underestimator.
 * This is done either explicitly (if n=2) or by solving an LP, see SCIPcomputeFacetVertexPolyhedralNonlinear().
 */
static
SCIP_DECL_EXPRESTIMATE(estimateProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int nchildren;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(refpoint != NULL);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   *success = FALSE;
   *islocal = TRUE;

   nchildren = SCIPexprGetNChildren(expr);

   /* debug output: prints expression we are trying to estimate, bounds of variables and point */
#ifdef SCIP_DEBUG
   {
      int c;

      SCIPdebugMsg(scip, "%sestimating product with %d variables\n", overestimate ? "over": "under", SCIPexprGetNChildren(expr));
      for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      {
         SCIPdebugMsg(scip, "child %d = %g in [%g, %g]\n", c, refpoint[c], localbounds[c].inf, localbounds[c].sup);

         if( SCIPisInfinity(scip, localbounds[c].sup) || SCIPisInfinity(scip, -localbounds[c].inf) )
         {
            SCIPdebugMsg(scip, "unbounded factor related to\n");
            SCIP_CALL( SCIPdismantleExpr(scip, NULL, SCIPexprGetChildren(expr)[0]) );
         }
      }
   }
#endif

   /* bilinear term */
   if( nchildren == 2 )
   {
      SCIP_Real refpointx;
      SCIP_Real refpointy;
      SCIP_INTERVAL bndx;
      SCIP_INTERVAL bndy;

      /* collect first variable */
      refpointx = refpoint[0];
      bndx = localbounds[0];
      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bndx));

      /* collect second variable */
      refpointy = refpoint[1];
      bndy = localbounds[1];
      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bndy));

      /* adjust the reference points */
      refpointx = MIN(MAX(refpointx, bndx.inf), bndx.sup);
      refpointy = MIN(MAX(refpointy, bndy.inf), bndy.sup);

      coefs[0] = 0.0;
      coefs[1] = 0.0;
      *constant = 0.0;
      *success = TRUE;

      SCIPaddBilinMcCormick(scip, exprdata->coefficient, bndx.inf, bndx.sup, refpointx,
            bndy.inf, bndy.sup, refpointy, overestimate, &coefs[0], &coefs[1], constant,
            success);
   }
   else
   {
      SCIP_EXPRHDLRDATA* exprhdlrdata;

      exprhdlrdata = SCIPexprhdlrGetData(SCIPexprGetHdlr(expr));
      assert(exprhdlrdata != NULL);

      if( exprhdlrdata->conshdlr != NULL )
      {
         SCIP_CALL( estimateVertexPolyhedralProduct(scip, exprhdlrdata->conshdlr, nchildren, localbounds, exprdata->coefficient, refpoint, overestimate,
            targetvalue, coefs, constant, success) );
      }
      else
      {
         SCIPdebugMsg(scip, "no cons_nonlinear included in SCIP, cannot estimate vertex-polyhedral product function\n");
      }
   }

   return SCIP_OKAY;
}

/** initial estimators callback */
static
SCIP_DECL_EXPRINITESTIMATES(initestimatesProduct)
{
   SCIP_EXPRDATA* exprdata;
   SCIP_Bool success = TRUE;
   int nchildren;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(nreturned != NULL);

   nchildren = SCIPexprGetNChildren(expr);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   if( nchildren == 2 )
   {
      SCIP_INTERVAL bndx = bounds[0];
      SCIP_INTERVAL bndy = bounds[1];

      constant[0] = 0.0;
      coefs[0][0] = 0.0;
      coefs[0][1] = 0.0;

      /* build estimator */
      SCIPaddBilinMcCormick(scip, exprdata->coefficient, bndx.inf, bndx.sup, (bndx.inf + bndx.sup) / 2.0,
         bndy.inf, bndy.sup, (bndy.inf + bndy.sup ) / 2.0, overestimate, &coefs[0][0], &coefs[0][1],
         constant, &success);
   }
   else
   {
      SCIP_EXPRHDLRDATA* exprhdlrdata;

      exprhdlrdata = SCIPexprhdlrGetData(SCIPexprGetHdlr(expr));
      assert(exprhdlrdata != NULL);

      if( exprhdlrdata->conshdlr != NULL )
      {
         SCIP_CALL( estimateVertexPolyhedralProduct(scip, exprhdlrdata->conshdlr, nchildren, bounds, exprdata->coefficient, NULL, overestimate,
            overestimate ? SCIPinfinity(scip) : -SCIPinfinity(scip), coefs[0], constant, &success) );
      }
      else
      {
         SCIPdebugMsg(scip, "no cons_nonlinear included in SCIP, cannot estimate vertex-polyhedral product function\n");
      }
   }

   if( success )
      *nreturned = 1;

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   SCIP_INTERVAL childbounds;
   SCIP_INTERVAL otherfactor;
   SCIP_INTERVAL zero;
   int i;
   int j;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) > 0);
   assert(infeasible != NULL);
   assert(childrenbounds != NULL);

   *infeasible = FALSE;

   /* too expensive (runtime here is quadratic in number of children)
    * TODO implement something faster for larger numbers of factors, e.g., split product into smaller products
    */
   if( SCIPexprGetNChildren(expr) > 10 )
      return SCIP_OKAY;

   /* not possible to learn bounds on children if expression bounds are unbounded in both directions */
   if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, bounds) )
      return SCIP_OKAY;

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   /**! [SnippetExprReversepropProduct] */
   SCIPintervalSet(&zero, 0.0);

   /* f = const * prod_k c_k => c_i solves c_i * (const * prod_{j:j!=i} c_j) = f */
   for( i = 0; i < SCIPexprGetNChildren(expr) && !(*infeasible); ++i )
   {
      SCIPintervalSet(&otherfactor, exprdata->coefficient);

      /* compute prod_{j:j!=i} c_j */
      for( j = 0; j < SCIPexprGetNChildren(expr); ++j )
      {
         if( i == j )
            continue;

         /* TODO we should compute these only one time instead of repeating this for almost every i */
         childbounds = childrenbounds[j];
         if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childbounds) )
         {
            *infeasible = TRUE;
            return SCIP_OKAY;
         }

         SCIPintervalMul(SCIP_INTERVAL_INFINITY, &otherfactor, otherfactor, childbounds);
      }

      childbounds = childrenbounds[i];
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childbounds) )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }

      /* solve x*otherfactor = f for x in c_i */
      SCIPintervalSolveUnivariateQuadExpression(SCIP_INTERVAL_INFINITY, &childbounds, zero, otherfactor, bounds, childbounds);

      SCIPdebugMsg(scip, "child %d: solved [%g,%g]*x = [%g,%g] with x in [%g,%g] -> x = [%g,%g]\n", i, otherfactor.inf, otherfactor.sup,
         bounds.inf, bounds.sup,
         childrenbounds[i].inf, childrenbounds[i].sup,
         childbounds.inf, childbounds.sup);

      /* store computed bounds of the expression */
      SCIPintervalIntersect(&childrenbounds[i], childrenbounds[i], childbounds);
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childrenbounds[i]) )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
   }
   /**! [SnippetExprReversepropProduct] */

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureProduct)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) > 0);

   if( SCIPexprGetNChildren(expr) == 1 )
   {
      *childcurv = SCIPexprcurvMultiply(SCIPgetCoefExprProduct(expr), exprcurvature);
      *success = TRUE;
   }
   else
   {
      *success = FALSE;
   }

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityProduct)
{  /*lint --e{715}*/
   SCIP_Real coef;
   int i;
   int nneg;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(SCIPexprGetNChildren(expr) >= 1);
   assert(childidx >= 0);
   assert(childidx < SCIPexprGetNChildren(expr));

   coef = SCIPgetCoefExprProduct(expr);

   /* count the number of negative children (except for childidx); if some children changes sign
    * -> monotonicity unknown
    */
   nneg = 0;
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_INTERVAL interval;

      if( i == childidx )
         continue;

      assert(SCIPexprGetChildren(expr)[i] != NULL);
      SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(expr)[i]) );
      interval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[i]);

      if( SCIPintervalGetSup(interval) <= 0.0 )
         nneg++;
      else if( SCIPintervalGetInf(interval) < 0.0 )
      {
         *result = SCIP_MONOTONE_UNKNOWN;
         return SCIP_OKAY;
      }
   }

   /* note that the monotonicity depends on the sign of the coefficient */
   if( nneg % 2 == 0 )
      *result = (coef >= 0.0) ? SCIP_MONOTONE_INC : SCIP_MONOTONE_DEC;
   else
      *result = (coef >= 0.0) ? SCIP_MONOTONE_DEC : SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityProduct)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   *isintegral = EPSISINT(exprdata->coefficient, 0.0); /*lint !e835*/

   for( i = 0; i < SCIPexprGetNChildren(expr) && *isintegral; ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
      assert(child != NULL);

      *isintegral = SCIPexprIsIntegral(child);
   }

   return SCIP_OKAY;
}

/** creates the handler for product expressions and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrProduct(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLRDATA* exprhdlrdata;
   SCIP_EXPRHDLR* exprhdlr;

   /* allocate expression handler data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &exprhdlrdata) );
   exprhdlrdata->conshdlr = SCIPfindConshdlr(scip, "nonlinear");

   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalProduct,
         exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrProduct, freehdlrProduct);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataProduct, freedataProduct);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyProduct);
   SCIPexprhdlrSetCompare(exprhdlr, compareProduct);
   SCIPexprhdlrSetPrint(exprhdlr, printProduct);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalProduct);
   SCIPexprhdlrSetEstimate(exprhdlr, initestimatesProduct, estimateProduct);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropProduct);
   SCIPexprhdlrSetHash(exprhdlr, hashProduct);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffProduct, fwdiffProduct, bwfwdiffProduct);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureProduct);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityProduct);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityProduct);

   return SCIP_OKAY;
}

/** creates a product expression */
SCIP_RETCODE SCIPcreateExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children */
   SCIP_Real             coefficient,        /**< constant coefficient of product */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRDATA* exprdata;

   /**! [SnippetCreateExprProduct] */
   SCIP_CALL( SCIPallocBlockMemory(scip, &exprdata) );
   exprdata->coefficient  = coefficient;

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPgetExprhdlrProduct(scip), exprdata, nchildren, children, ownercreate, ownercreatedata) );
   /**! [SnippetCreateExprProduct] */

   return SCIP_OKAY;
}

/* from pub_expr.h */

/** gets the constant coefficient of a product expression */
SCIP_Real SCIPgetCoefExprProduct(
   SCIP_EXPR*            expr                /**< product expression */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   return exprdata->coefficient;
}
