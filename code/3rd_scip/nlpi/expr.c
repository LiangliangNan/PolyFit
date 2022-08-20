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

/**@file   nlpi/expr.c
 * @brief  methods for expressions, expression trees, expression graphs, and related
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 * @author Ingmar Vierhaus (exprparse)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "nlpi/pub_expr.h"
#include "nlpi/struct_expr.h"
#include "nlpi/exprinterpret.h"

#include "scip/intervalarith.h"
#include "scip/pub_misc.h"
#include "scip/misc.h"
#include "scip/pub_message.h"


#define SCIP_EXPRESSION_MAXCHILDEST 16       /**< estimate on maximal number of children */

/** sign of a value (-1 or +1)
 *
 *  0.0 has sign +1
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

/** ensures that a block memory array has at least a given size
 *
 *  if cursize is 0, then *array1 can be NULL
 */
#define ensureBlockMemoryArraySize(blkmem, array1, cursize, minsize)    \
   do {                                                                 \
      int __newsize;                                                    \
      assert((blkmem)  != NULL);                                        \
      if( *(cursize) >= (minsize) )                                     \
         break;                                                         \
      __newsize = calcGrowSize(minsize);                                \
      assert(__newsize >= (minsize));                                   \
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, array1, *(cursize), __newsize) ); \
      *(cursize) = __newsize;                                           \
   } while( FALSE )

#ifdef SCIP_DISABLED_CODE /* this macro is currently not used, which offends lint, so disable it */
/** ensures that two block memory arrays have at least a given size
 *
 *  if cursize is 0, then arrays can be NULL
 */
#define ensureBlockMemoryArraySize2(blkmem, array1, array2, cursize, minsize) \
   do {                                                                 \
      int __newsize;                                                    \
      assert((blkmem)  != NULL);                                        \
      if( *(cursize) >= (minsize) )                                     \
         break;                                                         \
      __newsize = calcGrowSize(minsize);                                \
      assert(__newsize >= (minsize));                                   \
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, array1, *(cursize), __newsize) ); \
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, array2, *(cursize), __newsize) ); \
      *(cursize) = __newsize;                                           \
   } while( FALSE )
#endif

/** ensures that three block memory arrays have at least a given size
 *
 *  if cursize is 0, then arrays can be NULL
 */
#define ensureBlockMemoryArraySize3(blkmem, array1, array2, array3, cursize, minsize) \
   do {                                                                 \
      int __newsize;                                                    \
      assert((blkmem)  != NULL);                                        \
      if( *(cursize) >= (minsize) )                                     \
         break;                                                         \
      __newsize = calcGrowSize(minsize);                                \
      assert(__newsize >= (minsize));                                   \
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, array1, *(cursize), __newsize) ); \
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, array2, *(cursize), __newsize) ); \
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, array3, *(cursize), __newsize) ); \
      *(cursize) = __newsize;                                           \
   } while( FALSE )

/**@name Miscellaneous private methods */
/**@{ */

/** calculate memory size for dynamically allocated arrays (copied from scip/set.c) */
static
int calcGrowSize(
   int                   num                 /**< minimum number of entries to store */
   )
{
   int size;

   /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
   size = 4;
   while( size < num )
      size = (int)(1.2 * size + 4);

   return size;
}

/** expression graph nodes comparison to use in sorting methods
 *
 * The nodes need to have been added to the expression graph (depth,pos >= 0).
 * The better node is the one with the lower depth and lower position, if depth is equal.
 */
static
SCIP_DECL_SORTPTRCOMP(exprgraphnodecomp)
{
   SCIP_EXPRGRAPHNODE* node1 = (SCIP_EXPRGRAPHNODE*)elem1;
   SCIP_EXPRGRAPHNODE* node2 = (SCIP_EXPRGRAPHNODE*)elem2;

   assert(node1 != NULL);
   assert(node2 != NULL);
   assert(node1->depth >= 0);
   assert(node1->pos >= 0);
   assert(node2->depth >= 0);
   assert(node2->pos >= 0);

   if( node1->depth != node2->depth )
      return node1->depth - node2->depth;

   /* there should be no two nodes on the same position */
   assert((node1->pos != node2->pos) || (node1 == node2));

   return node1->pos - node2->pos;
}

/** checks if a given new lower bound is tighter (w.r.t. given bound strengthening epsilon) than the old one (copied from scip/set.c) */
static
SCIP_Bool isLbBetter(
   SCIP_Real             minstrength,        /**< minimal relative improvement required to be a better bound */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   )
{
   SCIP_Real eps;

   /* nothing can be tighter than an empty interval */
   if( oldlb > oldub )
      return FALSE;

   eps = REALABS(oldlb);
   eps = MIN(oldub - oldlb, eps);
   return EPSGT(newlb, oldlb, minstrength * MAX(eps, 1e-3));
}

/** checks if a given new upper bound is tighter (w.r.t. given bound strengthening epsilon) than the old one (copied from scip/set.c) */
static
SCIP_Bool isUbBetter(
   SCIP_Real             minstrength,        /**< minimal relative improvement required to be a better bound */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   )
{
   SCIP_Real eps;

   /* nothing can be tighter than an empty interval */
   if( oldlb > oldub )
      return FALSE;

   eps = REALABS(oldub);
   eps = MIN(oldub - oldlb, eps);
   return EPSLT(newub, oldub, minstrength * MAX(eps, 1e-3));
}

/**@} */

/**@name Expression curvature methods */
/**@{ */

/** curvature names as strings */
static
const char* curvnames[4] =
   {
      "unknown",
      "convex",
      "concave",
      "linear"
   };

#undef SCIPexprcurvAdd

/** gives curvature for a sum of two functions with given curvature */
SCIP_EXPRCURV SCIPexprcurvAdd(
   SCIP_EXPRCURV         curv1,              /**< curvature of first summand */
   SCIP_EXPRCURV         curv2               /**< curvature of second summand */
   )
{
   return (SCIP_EXPRCURV) (curv1 & curv2);
}

/** gives the curvature for the negation of a function with given curvature */
SCIP_EXPRCURV SCIPexprcurvNegate(
   SCIP_EXPRCURV         curvature           /**< curvature of function */
   )
{
   switch( curvature )
   {
   case SCIP_EXPRCURV_CONCAVE:
      return SCIP_EXPRCURV_CONVEX;

   case SCIP_EXPRCURV_CONVEX:
      return SCIP_EXPRCURV_CONCAVE;

   case SCIP_EXPRCURV_LINEAR:
   case SCIP_EXPRCURV_UNKNOWN:
      /* can return curvature, do this below */
      break;

   default:
      SCIPerrorMessage("unknown curvature status.\n");
      SCIPABORT();
   }

   return curvature;
}

/** gives curvature for a functions with given curvature multiplied by a constant factor */
SCIP_EXPRCURV SCIPexprcurvMultiply(
   SCIP_Real             factor,             /**< constant factor */
   SCIP_EXPRCURV         curvature           /**< curvature of other factor */
   )
{
   if( factor == 0.0 )
      return SCIP_EXPRCURV_LINEAR;
   if( factor > 0.0 )
      return curvature;
   return SCIPexprcurvNegate(curvature);
}

/** gives curvature for base^exponent for given bounds and curvature of base-function and constant exponent */
SCIP_EXPRCURV SCIPexprcurvPower(
   SCIP_INTERVAL         basebounds,         /**< bounds on base function */
   SCIP_EXPRCURV         basecurv,           /**< curvature of base function */
   SCIP_Real             exponent            /**< exponent */
   )
{
   SCIP_Bool expisint;

   assert(basebounds.inf <= basebounds.sup);

   if( exponent == 0.0 )
      return SCIP_EXPRCURV_LINEAR;

   if( exponent == 1.0 )
      return basecurv;

   expisint = EPSISINT(exponent, 0.0); /*lint !e835*/

   /* if exponent is fractional, then power is not defined for a negative base
    * thus, consider only positive part of basebounds
    */
   if( !expisint && basebounds.inf < 0.0 )
   {
      basebounds.inf = 0.0;
      if( basebounds.sup < 0.0 )
         return SCIP_EXPRCURV_LINEAR;
   }

   /* if basebounds contains 0.0, consider negative and positive interval separately, if possible */
   if( basebounds.inf < 0.0 && basebounds.sup > 0.0 )
   {
      SCIP_INTERVAL leftbounds;
      SCIP_INTERVAL rightbounds;

      /* something like x^(-2) may look convex on each side of zero, but is not convex on the whole interval due to the singularity at 0.0 */
      if( exponent < 0.0 )
         return SCIP_EXPRCURV_UNKNOWN;

      SCIPintervalSetBounds(&leftbounds,  basebounds.inf, 0.0);
      SCIPintervalSetBounds(&rightbounds, 0.0, basebounds.sup);

      return (SCIP_EXPRCURV) (SCIPexprcurvPower(leftbounds,  basecurv, exponent) & SCIPexprcurvPower(rightbounds, basecurv, exponent));
   }
   assert(basebounds.inf >= 0.0 || basebounds.sup <= 0.0);

   /* (base^exponent)'' = exponent * ( (exponent-1) base^(exponent-2) (base')^2 + base^(exponent-1) base'' )
    *
    * if base'' is positive, i.e., base is convex, then
    * - for base > 0.0 and exponent > 1.0, the second deriv. is positive -> convex
    * - for base < 0.0 and exponent > 1.0, we can't say (first and second summand opposite signs)
    * - for base > 0.0 and 0.0 < exponent < 1.0, we can't say (first sommand negative, second summand positive)
    * - for base > 0.0 and exponent < 0.0, we can't say (first and second summand opposite signs)
    * - for base < 0.0 and exponent < 0.0 and even, the second deriv. is positive -> convex
    * - for base < 0.0 and exponent < 0.0 and odd, the second deriv. is negative -> concave
    *
    * if base'' is negative, i.e., base is concave, then
    * - for base > 0.0 and exponent > 1.0, we can't say (first summand positive, second summand negative)
    * - for base < 0.0 and exponent > 1.0 and even, the second deriv. is positive -> convex
    * - for base < 0.0 and exponent > 1.0 and odd, the second deriv. is negative -> concave
    * - for base > 0.0 and 0.0 < exponent < 1.0, the second deriv. is negative -> concave
    * - for base > 0.0 and exponent < 0.0, the second deriv. is positive -> convex
    * - for base < 0.0 and exponent < 0.0, we can't say (first and second summand opposite signs)
    *
    * if base'' is zero, i.e., base is linear, then
    *   (base^exponent)'' = exponent * (exponent-1) base^(exponent-2) (base')^2
    * - just multiply signs
    */

   if( basecurv == SCIP_EXPRCURV_LINEAR )
   {
      SCIP_Real sign;

      /* base^(exponent-2) is negative, if base < 0.0 and exponent is odd */
      sign = exponent * (exponent - 1.0);
      assert(basebounds.inf >= 0.0 || expisint);
      if( basebounds.inf < 0.0 && ((int)exponent)%2 != 0 )
         sign *= -1.0;
      assert(sign != 0.0);

      return sign > 0.0 ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_CONCAVE;
   }

   if( basecurv == SCIP_EXPRCURV_CONVEX )
   {
      if( basebounds.sup <= 0.0 && exponent < 0.0 && expisint )
         return ((int)exponent)%2 == 0 ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_CONCAVE;
      if( basebounds.inf >= 0.0 && exponent > 1.0 )
         return SCIP_EXPRCURV_CONVEX ;
      return SCIP_EXPRCURV_UNKNOWN;
   }

   if( basecurv == SCIP_EXPRCURV_CONCAVE )
   {
      if( basebounds.sup <= 0.0 && exponent > 1.0 && expisint )
         return ((int)exponent)%2 == 0 ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_CONCAVE;
      if( basebounds.inf >= 0.0 && exponent < 1.0 )
         return exponent < 0.0 ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_CONCAVE;
      return SCIP_EXPRCURV_UNKNOWN;
   }

   return SCIP_EXPRCURV_UNKNOWN;
}

/** gives curvature for a monomial with given curvatures and bounds for each factor
 *
 *  See Maranas and Floudas, Finding All Solutions of Nonlinearly Constrained Systems of Equations, JOGO 7, 1995
 *  for the categorization in the case that all factors are linear.
 */
SCIP_EXPRCURV SCIPexprcurvMonomial(
   int                   nfactors,           /**< number of factors in monomial */
   SCIP_Real*            exponents,          /**< exponents in monomial, or NULL if all 1.0 */
   int*                  factoridxs,         /**< indices of factors (but not exponents), or NULL if identity mapping */
   SCIP_EXPRCURV*        factorcurv,         /**< curvature of each factor */
   SCIP_INTERVAL*        factorbounds        /**< bounds of each factor */
   )
{
   SCIP_Real mult;
   SCIP_Real e;
   SCIP_EXPRCURV curv;
   SCIP_EXPRCURV fcurv;
   int nnegative;
   int npositive;
   SCIP_Real sum;
   SCIP_Bool expcurvpos;
   SCIP_Bool expcurvneg;
   int j;
   int f;

   assert(nfactors >= 0);
   assert(factorcurv   != NULL || nfactors == 0);
   assert(factorbounds != NULL || nfactors == 0);

   if( nfactors == 0 )
      return SCIP_EXPRCURV_LINEAR;

   if( nfactors == 1 )
   {
      f = factoridxs != NULL ? factoridxs[0] : 0;
      e = exponents != NULL ? exponents[0] : 1.0;
      /* SCIPdebugMessage("monomial [%g,%g]^%g is %s\n",
         factorbounds[f].inf, factorbounds[f].sup, e,
         SCIPexprcurvGetName(SCIPexprcurvPower(factorbounds[f], factorcurv[f], e))); */
      return SCIPexprcurvPower(factorbounds[f], factorcurv[f], e);  /*lint !e613*/
   }

   mult = 1.0;

   nnegative = 0; /* number of negative exponents */
   npositive = 0; /* number of positive exponents */
   sum = 0.0;     /* sum of exponents */
   expcurvpos = TRUE; /* whether exp_j * f_j''(x) >= 0 for all factors (assuming f_j >= 0) */
   expcurvneg = TRUE; /* whether exp_j * f_j''(x) <= 0 for all factors (assuming f_j >= 0) */

   for( j = 0; j < nfactors; ++j )
   {
      f = factoridxs != NULL ? factoridxs[j] : j;
      if( factorcurv[f] == SCIP_EXPRCURV_UNKNOWN ) /*lint !e613*/
         return SCIP_EXPRCURV_UNKNOWN;
      if( factorbounds[f].inf < 0.0 && factorbounds[f].sup > 0.0 )  /*lint !e613*/
         return SCIP_EXPRCURV_UNKNOWN;

      e = exponents != NULL ? exponents[j] : 1.0;
      if( e < 0.0 )
         ++nnegative;
      else
         ++npositive;
      sum += e;

      if( factorbounds[f].inf < 0.0 )  /*lint !e613*/
      {
         /* if argument is negative, then exponent should be integer */
         assert(EPSISINT(e, 0.0));  /*lint !e835*/

         /* flip j'th argument: (f_j)^(exp_j) = (-1)^(exp_j) (-f_j)^(exp_j) */

         /* -f_j has negated curvature of f_j */
         fcurv = SCIPexprcurvNegate(factorcurv[f]);  /*lint !e613*/

         /* negate monomial, if exponent is odd, i.e., (-1)^(exp_j) = -1 */
         if( (int)e % 2 != 0 )
            mult *= -1.0;
      }
      else
      {
         fcurv = factorcurv[f];  /*lint !e613*/
      }

      /* check if exp_j * fcurv is convex (>= 0) and/or concave */
      fcurv = SCIPexprcurvMultiply(e, fcurv);
      if( !(fcurv & SCIP_EXPRCURV_CONVEX) )
         expcurvpos = FALSE;
      if( !(fcurv & SCIP_EXPRCURV_CONCAVE) )
         expcurvneg = FALSE;
   }

   /* if all factors are linear, then a product f_j^exp_j with f_j >= 0 is convex if
    * - all exponents are negative, or
    * - all except one exponent j* are negative and exp_j* >= 1 - sum_{j!=j*}exp_j, but the latter is equivalent to sum_j exp_j >= 1
    * further, the product is concave if
    * - all exponents are positive and the sum of exponents is <= 1.0
    *
    * if factors are nonlinear, then we require additionally, that for convexity
    * - each factor is convex if exp_j >= 0, or concave if exp_j <= 0, i.e., exp_j*f_j'' >= 0
    * and for concavity, we require that
    * - all factors are concave, i.e., exp_j*f_j'' <= 0
    */

   if( nnegative == nfactors && expcurvpos )
      curv = SCIP_EXPRCURV_CONVEX;
   else if( nnegative == nfactors-1 && EPSGE(sum, 1.0, 1e-9) && expcurvpos )
      curv = SCIP_EXPRCURV_CONVEX;
   else if( npositive == nfactors && EPSLE(sum, 1.0, 1e-9) && expcurvneg )
      curv = SCIP_EXPRCURV_CONCAVE;
   else
      curv = SCIP_EXPRCURV_UNKNOWN;
   curv = SCIPexprcurvMultiply(mult, curv);

   return curv;
}

/** gives name as string for a curvature */
const char* SCIPexprcurvGetName(
   SCIP_EXPRCURV         curv                /**< curvature */
   )
{
   assert(curv <= SCIP_EXPRCURV_LINEAR);  /*lint !e685*/

   return curvnames[curv];
}

/**@} */

/**@name Quadratic expression data private methods */
/**@{ */

/** creates SCIP_EXPRDATA_QUADRATIC data structure from given quadratic elements */
static
SCIP_RETCODE quadraticdataCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRDATA_QUADRATIC** quadraticdata,  /**< buffer to store pointer to quadratic data */
   SCIP_Real             constant,           /**< constant */
   int                   nchildren,          /**< number of children */
   SCIP_Real*            lincoefs,           /**< linear coefficients of children, or NULL if all 0.0 */
   int                   nquadelems,         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems           /**< quadratic elements */
   )
{
   assert(blkmem != NULL);
   assert(quadraticdata != NULL);
   assert(quadelems != NULL || nquadelems == 0);
   assert(nchildren >= 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, quadraticdata) );

   (*quadraticdata)->constant   = constant;
   (*quadraticdata)->lincoefs   = NULL;
   (*quadraticdata)->nquadelems = nquadelems;
   (*quadraticdata)->quadelems  = NULL;
   (*quadraticdata)->sorted     = (nquadelems <= 1);

   if( lincoefs != NULL )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*quadraticdata)->lincoefs, lincoefs, nchildren) );
   }

   if( nquadelems > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*quadraticdata)->quadelems, quadelems, nquadelems) );
   }

   return SCIP_OKAY;
}

/** sorts quadratic elements in a SCIP_EXPRDATA_QUADRATIC data structure */
static
void quadraticdataSort(
   SCIP_EXPRDATA_QUADRATIC* quadraticdata    /**< quadratic data */
   )
{
   assert(quadraticdata != NULL);

   if( quadraticdata->sorted )
   {
#ifndef NDEBUG
      int i;
      for( i = 1; i < quadraticdata->nquadelems; ++i )
      {
         assert(quadraticdata->quadelems[i].idx1 <= quadraticdata->quadelems[i].idx2);
         assert(quadraticdata->quadelems[i-1].idx1 <= quadraticdata->quadelems[i].idx1);
         assert(quadraticdata->quadelems[i-1].idx1 < quadraticdata->quadelems[i].idx1 ||
            quadraticdata->quadelems[i-1].idx2 <= quadraticdata->quadelems[i].idx2);
      }
#endif
      return;
   }

   if( quadraticdata->nquadelems > 0 )
      SCIPquadelemSort(quadraticdata->quadelems, quadraticdata->nquadelems);

   quadraticdata->sorted = TRUE;
}

/**@} */

/**@name Polynomial expression data private methods */
/**@{ */

/** compares two monomials
 *
 *  gives 0 if monomials are equal */
static
SCIP_DECL_SORTPTRCOMP(monomialdataCompare)
{
   SCIP_EXPRDATA_MONOMIAL* monomial1;
   SCIP_EXPRDATA_MONOMIAL* monomial2;

   int i;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   monomial1 = (SCIP_EXPRDATA_MONOMIAL*)elem1;
   monomial2 = (SCIP_EXPRDATA_MONOMIAL*)elem2;

   /* make sure, both monomials are equal */
   SCIPexprSortMonomialFactors(monomial1);
   SCIPexprSortMonomialFactors(monomial2);

   /* for the first factor where both monomials differ,
    * we return either the difference in the child indices, if children are different
    * or the sign of the difference in the exponents
    */
   for( i = 0; i < monomial1->nfactors && i < monomial2->nfactors; ++i )
   {
      if( monomial1->childidxs[i] != monomial2->childidxs[i] )
         return monomial1->childidxs[i] - monomial2->childidxs[i];
      if( monomial1->exponents[i] > monomial2->exponents[i] )
         return 1;
      else if( monomial1->exponents[i] < monomial2->exponents[i] )
         return -1;
   }

   /* if the factors of one monomial are a proper subset of the factors of the other monomial,
    * we return the difference in the number of monomials
    */
   return monomial1->nfactors - monomial2->nfactors;
}

/** ensures that the factors arrays of a monomial have at least a given size */
static
SCIP_RETCODE monomialdataEnsureFactorsSize(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRDATA_MONOMIAL*  monomialdata,    /**< monomial data */
   int                   minsize             /**< minimal size of factors arrays */
   )
{
   assert(blkmem != NULL);
   assert(monomialdata != NULL);

   if( minsize > monomialdata->factorssize )
   {
      int newsize;

      newsize = calcGrowSize(minsize);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &monomialdata->childidxs, monomialdata->factorssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &monomialdata->exponents, monomialdata->factorssize, newsize) );
      monomialdata->factorssize = newsize;
   }
   assert(minsize <= monomialdata->factorssize);

   return SCIP_OKAY;
}

/** creates SCIP_EXPRDATA_POLYNOMIAL data structure from given monomials */
static
SCIP_RETCODE polynomialdataCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRDATA_POLYNOMIAL** polynomialdata,/**< buffer to store pointer to polynomial data */
   int                   nmonomials,         /**< number of monomials */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< monomials */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool             copymonomials       /**< whether to copy monomials, or copy only given pointers, in which case polynomialdata assumes ownership of monomial structure */
   )
{
   assert(blkmem != NULL);
   assert(polynomialdata != NULL);
   assert(monomials != NULL || nmonomials == 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, polynomialdata) );

   (*polynomialdata)->constant = constant;
   (*polynomialdata)->nmonomials  = nmonomials;
   (*polynomialdata)->monomialssize = nmonomials;
   (*polynomialdata)->monomials   = NULL;
   (*polynomialdata)->sorted   = (nmonomials <= 1);

   if( nmonomials > 0 )
   {
      int i;

      if( copymonomials )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*polynomialdata)->monomials, nmonomials) );

         for( i = 0; i < nmonomials; ++i )
         {
            assert(monomials[i] != NULL);  /*lint !e613*/
            SCIP_CALL( SCIPexprCreateMonomial(blkmem, &(*polynomialdata)->monomials[i],
                  monomials[i]->coef, monomials[i]->nfactors, monomials[i]->childidxs, monomials[i]->exponents) );  /*lint !e613*/
         }
      }
      else
      {
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*polynomialdata)->monomials, monomials, nmonomials) );
      }
   }

   return SCIP_OKAY;
}

/** creates a copy of a SCIP_EXPRDATA_POLYNOMIAL data structure */
static
SCIP_RETCODE polynomialdataCopy(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRDATA_POLYNOMIAL** polynomialdata,/**< buffer to store pointer to polynomial data */
   SCIP_EXPRDATA_POLYNOMIAL* sourcepolynomialdata /**< polynomial data to copy */
   )
{
   assert(blkmem != NULL);
   assert(polynomialdata != NULL);
   assert(sourcepolynomialdata != NULL);

   SCIP_ALLOC( BMSduplicateBlockMemory(blkmem, polynomialdata, sourcepolynomialdata) );

   (*polynomialdata)->monomialssize = sourcepolynomialdata->nmonomials;
   if( sourcepolynomialdata->nmonomials > 0 )
   {
      int i;

      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*polynomialdata)->monomials, (*polynomialdata)->monomialssize) );

      for( i = 0; i < sourcepolynomialdata->nmonomials; ++i )
      {
         assert(sourcepolynomialdata->monomials[i] != NULL);  /*lint !e613*/
         SCIP_CALL( SCIPexprCreateMonomial(blkmem, &(*polynomialdata)->monomials[i], sourcepolynomialdata->monomials[i]->coef,
               sourcepolynomialdata->monomials[i]->nfactors, sourcepolynomialdata->monomials[i]->childidxs, sourcepolynomialdata->monomials[i]->exponents) );
         (*polynomialdata)->monomials[i]->sorted = sourcepolynomialdata->monomials[i]->sorted;
      }
   }
   else
   {
      (*polynomialdata)->monomials = NULL;
   }

   return SCIP_OKAY;
}

/** frees a SCIP_EXPRDATA_POLYNOMIAL data structure */
static
void polynomialdataFree(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRDATA_POLYNOMIAL** polynomialdata /**< pointer to polynomial data to free */
   )
{
   assert(blkmem != NULL);
   assert(polynomialdata != NULL);
   assert(*polynomialdata != NULL);

   if( (*polynomialdata)->monomialssize > 0 )
   {
      int i;

      for( i = 0; i < (*polynomialdata)->nmonomials; ++i )
      {
         assert((*polynomialdata)->monomials[i] != NULL);
         SCIPexprFreeMonomial(blkmem, &(*polynomialdata)->monomials[i]);
         assert((*polynomialdata)->monomials[i] == NULL);
      }

      BMSfreeBlockMemoryArray(blkmem, &(*polynomialdata)->monomials, (*polynomialdata)->monomialssize);
   }
   assert((*polynomialdata)->monomials == NULL);

   BMSfreeBlockMemory(blkmem, polynomialdata);
}

/** ensures that the monomials array of a polynomial has at least a given size */
static
SCIP_RETCODE polynomialdataEnsureMonomialsSize(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data */
   int                   minsize             /**< minimal size of monomials array */
   )
{
   assert(blkmem != NULL);
   assert(polynomialdata != NULL);

   ensureBlockMemoryArraySize(blkmem, &polynomialdata->monomials, &polynomialdata->monomialssize, minsize);
   assert(minsize <= polynomialdata->monomialssize);

   return SCIP_OKAY;
}

/** adds an array of monomials to a polynomial */
static
SCIP_RETCODE polynomialdataAddMonomials(
   BMS_BLKMEM*           blkmem,             /**< block memory of expression */
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data */
   int                   nmonomials,         /**< number of monomials to add */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< the monomials to add */
   SCIP_Bool             copymonomials       /**< whether to copy monomials or to assume ownership */
   )
{
   int i;

   assert(blkmem != NULL);
   assert(polynomialdata != NULL);
   assert(monomials != NULL || nmonomials == 0);

   if( nmonomials == 0 )
      return SCIP_OKAY;

   SCIP_CALL( polynomialdataEnsureMonomialsSize(blkmem, polynomialdata, polynomialdata->nmonomials + nmonomials) );
   assert(polynomialdata->monomialssize >= polynomialdata->nmonomials + nmonomials);

   if( copymonomials )
   {
      for( i = 0; i < nmonomials; ++i )
      {
         assert(monomials[i] != NULL);  /*lint !e613*/
         SCIP_CALL( SCIPexprCreateMonomial(blkmem, &polynomialdata->monomials[polynomialdata->nmonomials + i],
               monomials[i]->coef, monomials[i]->nfactors, monomials[i]->childidxs, monomials[i]->exponents) );  /*lint !e613*/
      }
   }
   else
   {
      BMScopyMemoryArray(&polynomialdata->monomials[polynomialdata->nmonomials], monomials, nmonomials);  /*lint !e866*/
   }
   polynomialdata->nmonomials += nmonomials;

   polynomialdata->sorted = (polynomialdata->nmonomials <= 1);

   return SCIP_OKAY;
}

/** ensures that monomials of a polynomial are sorted */
static
void polynomialdataSortMonomials(
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata  /**< polynomial expression */
   )
{
   assert(polynomialdata != NULL);

   if( polynomialdata->sorted )
   {
#ifndef NDEBUG
      int i;

      /* a polynom with more than one monoms can only be sorted if its monoms are sorted */
      for( i = 1; i < polynomialdata->nmonomials; ++i )
      {
         assert(polynomialdata->monomials[i-1]->sorted);
         assert(polynomialdata->monomials[i]->sorted);
         assert(monomialdataCompare(polynomialdata->monomials[i-1], polynomialdata->monomials[i]) <= 0);
      }
#endif
      return;
   }

   if( polynomialdata->nmonomials > 0 )
      SCIPsortPtr((void**)polynomialdata->monomials, monomialdataCompare, polynomialdata->nmonomials);

   polynomialdata->sorted = TRUE;
}

/** merges monomials that differ only in coefficient into a single monomial
 *
 *  Eliminates monomials with coefficient between -eps and eps.
 */
static
void polynomialdataMergeMonomials(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data */
   SCIP_Real             eps,                /**< threshold under which numbers are treat as zero */
   SCIP_Bool             mergefactors        /**< whether to merge factors in monomials too */
   )
{
   int i;
   int offset;
   int oldnfactors;

   assert(polynomialdata != NULL);
   assert(eps >= 0.0);

   polynomialdataSortMonomials(polynomialdata);

   /* merge monomials by adding their coefficients, eliminate monomials with no factors or zero coefficient*/
   offset = 0;
   i = 0;
   while( i + offset < polynomialdata->nmonomials )
   {
      if( offset > 0 )
      {
         assert(polynomialdata->monomials[i] == NULL);
         assert(polynomialdata->monomials[i+offset] != NULL);
         polynomialdata->monomials[i] = polynomialdata->monomials[i+offset];
#ifndef NDEBUG
         polynomialdata->monomials[i+offset] = NULL;
#endif
      }

      if( mergefactors )
      {
         oldnfactors = polynomialdata->monomials[i]->nfactors;
         SCIPexprMergeMonomialFactors(polynomialdata->monomials[i], eps);

         /* if monomial has changed, then we cannot assume anymore that polynomial is sorted */
         if( oldnfactors != polynomialdata->monomials[i]->nfactors )
            polynomialdata->sorted = FALSE;
      }

      while( i+offset+1 < polynomialdata->nmonomials )
      {
         assert(polynomialdata->monomials[i+offset+1] != NULL);
         if( mergefactors )
         {
            oldnfactors = polynomialdata->monomials[i+offset+1]->nfactors;
            SCIPexprMergeMonomialFactors(polynomialdata->monomials[i+offset+1], eps);

            /* if monomial has changed, then we cannot assume anymore that polynomial is sorted */
            if( oldnfactors != polynomialdata->monomials[i+offset+1]->nfactors )
               polynomialdata->sorted = FALSE;
         }
         if( monomialdataCompare((void*)polynomialdata->monomials[i], (void*)polynomialdata->monomials[i+offset+1]) != 0 )
            break;
         polynomialdata->monomials[i]->coef += polynomialdata->monomials[i+offset+1]->coef;
         SCIPexprFreeMonomial(blkmem, &polynomialdata->monomials[i+offset+1]);
         ++offset;
      }

      if( polynomialdata->monomials[i]->nfactors == 0 )
      {
         /* constant monomial */
         polynomialdata->constant += polynomialdata->monomials[i]->coef;
         SCIPexprFreeMonomial(blkmem, &polynomialdata->monomials[i]);
         ++offset;
         continue;
      }

      if( EPSZ(polynomialdata->monomials[i]->coef, eps) )
      {
         SCIPexprFreeMonomial(blkmem, &polynomialdata->monomials[i]);
         ++offset;
         continue;
      }

      ++i;
   }

#ifndef NDEBUG
   for( ; i < polynomialdata->nmonomials; ++i )
      assert(polynomialdata->monomials[i] == NULL);
#endif

   polynomialdata->nmonomials -= offset;

   if( EPSZ(polynomialdata->constant, eps) )
      polynomialdata->constant = 0.0;
}

/** multiplies each summand of a polynomial by a given constant */
static
void polynomialdataMultiplyByConstant(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data */
   SCIP_Real             factor              /**< constant factor */
   )
{
   int i;

   assert(polynomialdata != NULL);

   if( factor == 1.0 )
      return;

   if( factor == 0.0 )
   {
      for( i = 0; i < polynomialdata->nmonomials; ++i )
         SCIPexprFreeMonomial(blkmem, &polynomialdata->monomials[i]);
      polynomialdata->nmonomials = 0;
   }
   else
   {
      for( i = 0; i < polynomialdata->nmonomials; ++i )
         SCIPexprChgMonomialCoef(polynomialdata->monomials[i], polynomialdata->monomials[i]->coef * factor);
   }

   polynomialdata->constant *= factor;
}

/** multiplies each summand of a polynomial by a given monomial */
static
SCIP_RETCODE polynomialdataMultiplyByMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data */
   SCIP_EXPRDATA_MONOMIAL* factor,           /**< monomial factor */
   int*                  childmap            /**< map children in factor to children in expr, or NULL for 1:1 */
   )
{
   int i;

   assert(blkmem != NULL);
   assert(factor != NULL);
   assert(polynomialdata != NULL);

   if( factor->nfactors == 0 )
   {
      polynomialdataMultiplyByConstant(blkmem, polynomialdata, factor->coef);
      return SCIP_OKAY;
   }

   /* multiply each monomial by factor */
   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      SCIP_CALL( SCIPexprMultiplyMonomialByMonomial(blkmem, polynomialdata->monomials[i], factor, childmap) );
   }

   /* add new monomial for constant multiplied by factor */
   if( polynomialdata->constant != 0.0 )
   {
      SCIP_CALL( polynomialdataEnsureMonomialsSize(blkmem, polynomialdata, polynomialdata->nmonomials+1) );
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &polynomialdata->monomials[polynomialdata->nmonomials], polynomialdata->constant, 0, NULL, NULL) );
      SCIP_CALL( SCIPexprMultiplyMonomialByMonomial(blkmem, polynomialdata->monomials[polynomialdata->nmonomials], factor, childmap) );
      ++polynomialdata->nmonomials;
      polynomialdata->sorted = FALSE;
      polynomialdata->constant = 0.0;
   }

   return SCIP_OKAY;
}

/** multiplies a polynomial by a polynomial
 *
 *  Factors need to be different.
 */
static
SCIP_RETCODE polynomialdataMultiplyByPolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data */
   SCIP_EXPRDATA_POLYNOMIAL* factordata,     /**< polynomial factor data */
   int*                  childmap            /**< map children in factor to children in polynomialdata, or NULL for 1:1 */
   )
{
   int i1;
   int i2;
   int orignmonomials;

   assert(blkmem != NULL);
   assert(polynomialdata != NULL);
   assert(factordata != NULL);
   assert(polynomialdata != factordata);

   if( factordata->nmonomials == 0 )
   {
      polynomialdataMultiplyByConstant(blkmem, polynomialdata, factordata->constant);
      return SCIP_OKAY;
   }

   if( factordata->nmonomials == 1 && factordata->constant == 0.0 )
   {
      SCIP_CALL( polynomialdataMultiplyByMonomial(blkmem, polynomialdata, factordata->monomials[0], childmap) );
      return SCIP_OKAY;
   }

   /* turn constant into a monomial, so we can assume below that constant is 0.0 */
   if( polynomialdata->constant != 0.0 )
   {
      SCIP_CALL( polynomialdataEnsureMonomialsSize(blkmem, polynomialdata, polynomialdata->nmonomials+1) );
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &polynomialdata->monomials[polynomialdata->nmonomials], polynomialdata->constant, 0, NULL, NULL) );
      ++polynomialdata->nmonomials;
      polynomialdata->sorted = FALSE;
      polynomialdata->constant = 0.0;
   }

   SCIP_CALL( polynomialdataEnsureMonomialsSize(blkmem, polynomialdata, polynomialdata->nmonomials * (factordata->nmonomials + (factordata->constant == 0.0 ? 0 : 1))) );

   /* for each monomial in factordata (except the last, if factordata->constant is 0),
    * duplicate monomials from polynomialdata and multiply them by the monomial for factordata */
   orignmonomials = polynomialdata->nmonomials;
   for( i2 = 0; i2 < factordata->nmonomials; ++i2 )
   {
      /* add a copy of original monomials to end of polynomialdata's monomials array */
      assert(polynomialdata->nmonomials + orignmonomials <= polynomialdata->monomialssize); /* reallocating in polynomialdataAddMonomials would make the polynomialdata->monomials invalid, so assert that above the monomials array was made large enough */
      SCIP_CALL( polynomialdataAddMonomials(blkmem, polynomialdata, orignmonomials, polynomialdata->monomials, TRUE) );
      assert(polynomialdata->nmonomials == (i2+2) * orignmonomials);

      /* multiply each copied monomial by current monomial from factordata */
      for( i1 = (i2+1) * orignmonomials; i1 < (i2+2) * orignmonomials; ++i1 )
      {
         SCIP_CALL( SCIPexprMultiplyMonomialByMonomial(blkmem, polynomialdata->monomials[i1], factordata->monomials[i2], childmap) );
      }

      if( factordata->constant == 0.0 && i2 == factordata->nmonomials - 2 )
      {
         ++i2;
         break;
      }
   }

   if( factordata->constant != 0.0 )
   {
      assert(i2 == factordata->nmonomials);
      /* multiply original monomials in polynomialdata by constant in factordata */
      for( i1 = 0; i1 < orignmonomials; ++i1 )
         SCIPexprChgMonomialCoef(polynomialdata->monomials[i1], polynomialdata->monomials[i1]->coef * factordata->constant);
   }
   else
   {
      assert(i2 == factordata->nmonomials - 1);
      /* multiply original monomials in polynomialdata by last monomial in factordata */
      for( i1 = 0; i1 < orignmonomials; ++i1 )
      {
         SCIP_CALL( SCIPexprMultiplyMonomialByMonomial(blkmem, polynomialdata->monomials[i1], factordata->monomials[i2], childmap) );
      }
   }

   return SCIP_OKAY;
}

/** takes a power of a polynomial
 *
 *  Exponent needs to be an integer,
 *  polynomial needs to be a monomial, if exponent is negative.
 */
static
SCIP_RETCODE polynomialdataPower(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data */
   int                   exponent            /**< exponent of power operation */
   )
{
   SCIP_EXPRDATA_POLYNOMIAL* factor;
   int i;

   assert(blkmem != NULL);
   assert(polynomialdata != NULL);

   if( exponent == 0 )
   {
      /* x^0 = 1, except if x = 0 */
      if( polynomialdata->nmonomials == 0 && polynomialdata->constant == 0.0 )
      {
         polynomialdata->constant = 0.0;
      }
      else
      {
         polynomialdata->constant = 1.0;

         for( i = 0; i < polynomialdata->nmonomials; ++i )
            SCIPexprFreeMonomial(blkmem, &polynomialdata->monomials[i]);
         polynomialdata->nmonomials = 0;
      }

      return SCIP_OKAY;
   }

   if( exponent == 1 )
      return SCIP_OKAY;

   if( polynomialdata->nmonomials == 1 && polynomialdata->constant == 0.0 )
   {
      /* polynomial is a single monomial */
      SCIPexprMonomialPower(polynomialdata->monomials[0], exponent);
      return SCIP_OKAY;
   }

   if( polynomialdata->nmonomials == 0 )
   {
      /* polynomial is a constant */
      polynomialdata->constant = pow(polynomialdata->constant, (SCIP_Real)exponent);
      return SCIP_OKAY;
   }

   assert(exponent >= 2); /* negative exponents not allowed if more than one monom */

   /* todo improve, look how SCIPintervalPowerScalar in intervalarith.c does it */

   /* get copy of our polynomial */
   SCIP_CALL( polynomialdataCopy(blkmem, &factor, polynomialdata) );

   /* do repeated multiplication */
   for( i = 2; i <= exponent; ++i )
   {
      SCIP_CALL( polynomialdataMultiplyByPolynomial(blkmem, polynomialdata, factor, NULL) );
      polynomialdataMergeMonomials(blkmem, polynomialdata, 0.0, TRUE);
   }

   /* free copy again */
   polynomialdataFree(blkmem, &factor);

   return SCIP_OKAY;
}

/** applies a mapping of child indices to the indices used in polynomial monomials */
static
void polynomialdataApplyChildmap(
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data */
   int*                  childmap            /**< mapping of child indices */
   )
{
   SCIP_EXPRDATA_MONOMIAL* monomial;
   int i;
   int j;

   assert(polynomialdata != NULL);

   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      monomial = polynomialdata->monomials[i];
      assert(monomial != NULL);

      for( j = 0; j < monomial->nfactors; ++j )
      {
         monomial->childidxs[j] = childmap[monomial->childidxs[j]];
         assert(monomial->childidxs[j] >= 0);
      }
      monomial->sorted = FALSE;
   }

   polynomialdata->sorted = FALSE;
}

/** replaces a factor in a monomial by a polynomial and expands the result */
static
SCIP_RETCODE polynomialdataExpandMonomialFactor(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata, /**< polynomial data where to expand a monomial */
   int                   monomialpos,        /**< position of monomial which factor to expand */
   int                   factorpos,          /**< position of factor in monomial to expand */
   SCIP_EXPRDATA_POLYNOMIAL* factorpolynomial,/**< polynomial that should replace factor */
   int*                  childmap,           /**< map of child indices in factorpolynomial to children of polynomial */
   int                   maxexpansionexponent,/**< maximal exponent for which polynomials (with > 1 summands) are expanded */
   SCIP_Bool*            success             /**< buffer to store whether expansion has been done */
   )
{
   SCIP_EXPRDATA_POLYNOMIAL* factorpolynomialcopy;
   SCIP_EXPRDATA_MONOMIAL* monomial;
   int i;

   assert(blkmem != NULL);
   assert(polynomialdata != NULL);
   assert(factorpolynomial != NULL);
   assert(childmap != NULL || factorpolynomial->nmonomials == 0);
   assert(success != NULL);
   assert(monomialpos >= 0);
   assert(monomialpos < polynomialdata->nmonomials);
   assert(factorpos >= 0);

   monomial = polynomialdata->monomials[monomialpos];
   assert(monomial != NULL);
   assert(factorpos < monomial->nfactors);

   *success = TRUE;

   if( factorpolynomial->nmonomials == 0 )
   {
      /* factorpolynomial is a constant */

      if( !EPSISINT(monomial->exponents[factorpos], 0.0) && factorpolynomial->constant < 0.0 )  /*lint !e835*/
      {
         /* if polynomial is a negative constant and our exponent is not integer, then cannot do expansion */
         SCIPmessagePrintWarning(messagehdlr, "got negative constant %g to the power of a noninteger exponent %g\n", factorpolynomial->constant, monomial->exponents[factorpos]);
         *success = FALSE;
         return SCIP_OKAY;
      }
      monomial->coef *= pow(factorpolynomial->constant, monomial->exponents[factorpos]);

      /* move last factor to position factorpos */
      if( factorpos < monomial->nfactors-1 )
      {
         monomial->exponents[factorpos] = monomial->exponents[monomial->nfactors-1];
         monomial->childidxs[factorpos] = monomial->childidxs[monomial->nfactors-1];
      }
      --monomial->nfactors;
      monomial->sorted = FALSE;
      polynomialdata->sorted = FALSE;

      return SCIP_OKAY;
   }

   if( factorpolynomial->constant == 0.0 && factorpolynomial->nmonomials == 1 )
   {
      /* factorpolynomial is a single monomial */
      SCIP_EXPRDATA_MONOMIAL* factormonomial;
      int childidx;
      SCIP_Real exponent;

      factormonomial = factorpolynomial->monomials[0];
      assert(factormonomial != NULL);

      if( !EPSISINT(monomial->exponents[factorpos], 0.0) )  /*lint !e835*/
      {
         if( factormonomial->coef < 0.0 )
         {
            /* if coefficient of monomial is negative and our exponent is not integer, then do not do expansion
             * @todo the only case where this could make sense is if the factors can be negative, i.e., when we have negative arguments with an odd exponent: (-x^a)^b = (-x)^(ab) for a odd
             */
            *success = FALSE;
            return SCIP_OKAY;
         }
         if( factormonomial->nfactors > 1 )
         {
            /* @todo if there is an even number of factors in factormonomial that are negative, then they always multiply to something positive
             * however, we cannot expand them as below, since we cannot compute the single powers
             * since we do not have the bounds on the factors here, we skip expansion in this case
             * MINLPLib instances tls2,4,6 are examples where we are loosing here (do not recognize convexity)
             */
            *success = FALSE;
            return SCIP_OKAY;
         }
      }

      SCIP_CALL( monomialdataEnsureFactorsSize(blkmem, monomial, monomial->nfactors + factormonomial->nfactors) );

      for( i = 0; i < factormonomial->nfactors; ++i )
      {
         childidx = childmap[factormonomial->childidxs[i]];  /*lint !e613*/
         /* can do this because monomial->exponents[factorpos] is assumed to be integer or factormonomial has positive coefficient and only one factor
          * thus, if factormonomial->exponents[i] is fractional, then we can assume that it's argument is positive
          */
         exponent = factormonomial->exponents[i] * monomial->exponents[factorpos];
         SCIP_CALL( SCIPexprAddMonomialFactors(blkmem, monomial, 1, &childidx, &exponent) );
      }

      monomial->coef *= pow(factormonomial->coef, monomial->exponents[factorpos]);

      /* move last factor to position factorpos */
      if( factorpos < monomial->nfactors-1 )
      {
         monomial->exponents[factorpos] = monomial->exponents[monomial->nfactors-1];
         monomial->childidxs[factorpos] = monomial->childidxs[monomial->nfactors-1];
      }
      --monomial->nfactors;
      monomial->sorted = FALSE;
      polynomialdata->sorted = FALSE;

      return SCIP_OKAY;
   }

   /* if exponent is negative or fractional and the polynomial is not just a monomial, then we cannot do expansion */
   if( !EPSISINT(monomial->exponents[factorpos], 0.0) || monomial->exponents[factorpos] < 0.0 )  /*lint !e835*/
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* if exponent is too large, skip expansion */
   if( monomial->exponents[factorpos] > maxexpansionexponent )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* check whether maximal degree of expansion would exceed maxexpansionexponent
    * that is, assume monomial is f1^a1 f2^a2 ... and we want to expand f1 = (g11^beta11 g12^beta12... + g21^beta21 g22^beta22 ... + ...)
    * then we do this only if all ai and all beta are > 0.0 and a1 max(beta11+beta12+..., beta21+beta22+..., ...) + a2 + ... < maxexpansionexponent
    * exception (there need to be one) is if monomial is just f1
    */
   if( maxexpansionexponent < INT_MAX && (monomial->nfactors > 1 || monomial->exponents[factorpos] != 1.0) )
   {
      SCIP_Real restdegree;
      SCIP_Real degree;
      int j;

      restdegree = -monomial->exponents[factorpos];
      for( i = 0; i < monomial->nfactors; ++i )
      {
         if( monomial->exponents[i] < 0.0 )
         {
            /* ai < 0.0 */
            SCIPdebugMessage("skip expansion because factor %d in monomial has negative exponent\n", i);
            *success = FALSE;
            return SCIP_OKAY;
         }
         restdegree += monomial->exponents[i];
      }

      for( i = 0; i < factorpolynomial->nmonomials; ++i )
      {
         degree = 0.0;
         for( j = 0; j < factorpolynomial->monomials[i]->nfactors; ++j )
         {
            if( factorpolynomial->monomials[i]->exponents[j] < 0.0 )
            {
               /* beta_ij < 0.0 */
               SCIPdebugMessage("skip expansion because %d'th factor in %d'th monomial of factorpolynomial is negative\n", i, j);
               *success = FALSE;
               return SCIP_OKAY;
            }
            degree += factorpolynomial->monomials[i]->exponents[j];
         }
         if( degree * monomial->exponents[factorpos] + restdegree > maxexpansionexponent )
         {
            /* (beta_i1+beta_i2+...)*monomial->exponents[factorpos] + rest > maxexpansion */
            SCIPdebugMessage("skip expansion because degree of %d'th monomial would yield degree %g > max = %d in expansion\n",
               i, degree * monomial->exponents[factorpos] + restdegree, maxexpansionexponent);
            *success = FALSE;
            return SCIP_OKAY;
         }
      }
   }

   /* create a copy of factor */
   SCIP_CALL( polynomialdataCopy(blkmem, &factorpolynomialcopy, factorpolynomial) );
   /* apply childmap to copy */
   polynomialdataApplyChildmap(factorpolynomialcopy, childmap);
   /* create power of factor */
   SCIP_CALL( polynomialdataPower(blkmem, factorpolynomialcopy, (int)EPSFLOOR(monomial->exponents[factorpos], 0.0)) );  /*lint !e835*/

   /* remove factor from monomial by moving last factor to position factorpos */
   if( factorpos < monomial->nfactors-1 )
   {
      monomial->exponents[factorpos] = monomial->exponents[monomial->nfactors-1];
      monomial->childidxs[factorpos] = monomial->childidxs[monomial->nfactors-1];
   }
   --monomial->nfactors;
   monomial->sorted = FALSE;

   /* multiply factor with this reduced monomial */
   SCIP_CALL( polynomialdataMultiplyByMonomial(blkmem, factorpolynomialcopy, monomial, NULL) );

   /* remove monomial from polynomial and move last monomial to monomialpos */
   SCIPexprFreeMonomial(blkmem, &polynomialdata->monomials[monomialpos]);
   if( monomialpos < polynomialdata->nmonomials-1 )
      polynomialdata->monomials[monomialpos] = polynomialdata->monomials[polynomialdata->nmonomials-1];
   --polynomialdata->nmonomials;
   polynomialdata->sorted = FALSE;

   /* add factorpolynomialcopy to polynomial */
   SCIP_CALL( polynomialdataAddMonomials(blkmem, polynomialdata, factorpolynomialcopy->nmonomials, factorpolynomialcopy->monomials, FALSE) );
   polynomialdata->constant += factorpolynomialcopy->constant;

   factorpolynomialcopy->nmonomials = 0;
   polynomialdataFree(blkmem, &factorpolynomialcopy);

   return SCIP_OKAY;
}

/**@} */

/**@name Expression operand private methods */
/**@{ */

/** a default implementation of expression interval evaluation that always gives a correct result */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntDefault )
{   /*lint --e{715}*/
   SCIPintervalSetEntire(infinity, result);

   return SCIP_OKAY;
}

/** a default implementation of expression curvature check that always gives a correct result */
static
SCIP_DECL_EXPRCURV( exprcurvDefault )
{   /*lint --e{715}*/
   *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_VAR */
static
SCIP_DECL_EXPREVAL( exprevalVar )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(varvals != NULL);

   *result = varvals[opdata.intval];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_VAR */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntVar )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(varvals != NULL);

   *result = varvals[opdata.intval];

   return SCIP_OKAY;
}

/** curvature for EXPR_VAR */
static
SCIP_DECL_EXPRCURV( exprcurvVar )
{   /*lint --e{715}*/
   assert(result  != NULL);

   *result = SCIP_EXPRCURV_LINEAR;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_CONST */
static
SCIP_DECL_EXPREVAL( exprevalConst )
{   /*lint --e{715}*/
   assert(result != NULL);

   *result = opdata.dbl;

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_CONST */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntConst )
{   /*lint --e{715}*/
   assert(result != NULL);

   SCIPintervalSet(result, opdata.dbl);

   return SCIP_OKAY;
}

/** curvature for EXPR_CONST */
static
SCIP_DECL_EXPRCURV( exprcurvConst )
{   /*lint --e{715}*/
   assert(result  != NULL);

   *result = SCIP_EXPRCURV_LINEAR;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_PARAM */
static
SCIP_DECL_EXPREVAL( exprevalParam )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(paramvals != NULL );

   *result = paramvals[opdata.intval];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_PARAM */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntParam )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(paramvals != NULL );

   SCIPintervalSet(result, paramvals[opdata.intval]);

   return SCIP_OKAY;
}

/** curvature for EXPR_PARAM */
static
SCIP_DECL_EXPRCURV( exprcurvParam )
{   /*lint --e{715}*/
   assert(result  != NULL);

   *result = SCIP_EXPRCURV_LINEAR;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_PLUS */
static
SCIP_DECL_EXPREVAL( exprevalPlus )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] + argvals[1];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_PLUS */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntPlus )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalAdd(infinity, result, argvals[0], argvals[1]);

   return SCIP_OKAY;
}

/** curvature for EXPR_PLUS */
static
SCIP_DECL_EXPRCURV( exprcurvPlus )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argcurv != NULL);

   *result = SCIPexprcurvAdd(argcurv[0], argcurv[1]);

   return SCIP_OKAY;
}

/** point evaluation for EXPR_MINUS */
static
SCIP_DECL_EXPREVAL( exprevalMinus )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] - argvals[1];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_MINUS */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntMinus )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalSub(infinity, result, argvals[0], argvals[1]);

   return SCIP_OKAY;
}

/** curvature for EXPR_MINUS */
static
SCIP_DECL_EXPRCURV( exprcurvMinus )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argcurv != NULL);

   *result = SCIPexprcurvAdd(argcurv[0], SCIPexprcurvNegate(argcurv[1]));

   return SCIP_OKAY;
}

/** point evaluation for EXPR_MUL */
static
SCIP_DECL_EXPREVAL( exprevalMult )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] * argvals[1];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_MUL */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntMult )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalMul(infinity, result, argvals[0], argvals[1]);

   return SCIP_OKAY;
}

/** curvature for EXPR_MUL */
static
SCIP_DECL_EXPRCURV( exprcurvMult )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   /* if one factor is constant, then product is
    * - linear, if constant is 0.0
    * - same curvature as other factor, if constant is positive
    * - negated curvature of other factor, if constant is negative
    *
    * if both factors are not constant, then product may not be convex nor concave
    */
   if( argbounds[1].inf == argbounds[1].sup )  /*lint !e777*/
      *result = SCIPexprcurvMultiply(argbounds[1].inf, argcurv[0]);
   else if( argbounds[0].inf == argbounds[0].sup )  /*lint !e777*/
      *result = SCIPexprcurvMultiply(argbounds[0].inf, argcurv[1]);
   else
      *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_DIV */
static
#if defined(__GNUC__) && __GNUC__ * 100 + __GNUC_MINOR__ * 10 >= 490 && !defined(__INTEL_COMPILER)
__attribute__((no_sanitize_undefined))
#endif
SCIP_DECL_EXPREVAL( exprevalDiv )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] / argvals[1];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_DIV */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntDiv )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalDiv(infinity, result, argvals[0], argvals[1]);

   return SCIP_OKAY;
}

/** curvature for EXPR_DIV */
static
SCIP_DECL_EXPRCURV( exprcurvDiv )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   /* if denominator is constant, then quotient has curvature sign(denominator) * curv(nominator)
    *
    * if nominator is a constant, then quotient is
    * - sign(nominator) * convex, if denominator is concave and positive
    * - sign(nominator) * concave, if denominator is convex and negative
    *
    * if denominator is positive but convex, then we don't know, e.g.,
    *   - 1/x^2 is convex for x>=0
    *   - 1/(1+(x-1)^2) is neither convex nor concave for x >= 0
    *
    * if both nominator and denominator are not constant, then quotient may not be convex nor concave
    */
   if( argbounds[1].inf == argbounds[1].sup )  /*lint !e777*/
   {
      /* denominator is constant */
      *result = SCIPexprcurvMultiply(argbounds[1].inf, argcurv[0]);
   }
   else if( argbounds[0].inf == argbounds[0].sup )  /*lint !e777*/
   {
      /* nominator is constant */
      if( argbounds[1].inf >= 0.0 && (argcurv[1] & SCIP_EXPRCURV_CONCAVE) )
         *result = SCIPexprcurvMultiply(argbounds[0].inf, SCIP_EXPRCURV_CONVEX);
      else if( argbounds[1].sup <= 0.0 && (argcurv[1] & SCIP_EXPRCURV_CONVEX) )
         *result = SCIPexprcurvMultiply(argbounds[0].inf, SCIP_EXPRCURV_CONCAVE);
      else
         *result = SCIP_EXPRCURV_UNKNOWN;
   }
   else
   {
      /* denominator and nominator not constant */
      *result = SCIP_EXPRCURV_UNKNOWN;
   }

   return SCIP_OKAY;
}

/** point evaluation for EXPR_SQUARE */
static
SCIP_DECL_EXPREVAL( exprevalSquare )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = argvals[0] * argvals[0];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_SQUARE */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntSquare )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalSquare(infinity, result, argvals[0]);

   return SCIP_OKAY;
}

/** curvature for EXPR_SQUARE */
static
SCIP_DECL_EXPRCURV( exprcurvSquare )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   *result = SCIPexprcurvPower(argbounds[0], argcurv[0], 2.0);

   return SCIP_OKAY;
}

/** point evaluation for EXPR_SQRT */
static
SCIP_DECL_EXPREVAL( exprevalSquareRoot )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = sqrt(argvals[0]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_SQRT */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntSquareRoot )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalSquareRoot(infinity, result, argvals[0]);

   return SCIP_OKAY;
}

/** curvature for EXPR_SQRT */
static
SCIP_DECL_EXPRCURV( exprcurvSquareRoot )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);

   /* square-root is concave, if child is concave
    * otherwise, we don't know
    */

   if( argcurv[0] & SCIP_EXPRCURV_CONCAVE )
      *result = SCIP_EXPRCURV_CONCAVE;
   else
      *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_REALPOWER */
static
SCIP_DECL_EXPREVAL( exprevalRealPower )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = pow(argvals[0], opdata.dbl);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_REALPOWER */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntRealPower )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalPowerScalar(infinity, result, argvals[0], opdata.dbl);

   return SCIP_OKAY;
}

/** curvature for EXPR_REALPOWER */
static
SCIP_DECL_EXPRCURV( exprcurvRealPower )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   *result = SCIPexprcurvPower(argbounds[0], argcurv[0], opdata.dbl);

   return SCIP_OKAY;
}

/** point evaluation for EXPR_INTPOWER */
static
#if defined(__GNUC__) && __GNUC__ * 100 + __GNUC_MINOR__ * 10 >= 490 && !defined(__INTEL_COMPILER)
__attribute__((no_sanitize_undefined))
#endif
SCIP_DECL_EXPREVAL( exprevalIntPower )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   switch( opdata.intval )
   {
   case -1:
      *result = 1.0 / argvals[0];
      return SCIP_OKAY;

   case 0:
      *result = 1.0;
      return SCIP_OKAY;

   case 1:
      *result = argvals[0];
      return SCIP_OKAY;

   case 2:
      *result = argvals[0] * argvals[0];
      return SCIP_OKAY;

   default:
      *result = pow(argvals[0], (SCIP_Real)opdata.intval);
   }

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_INTPOWER */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntIntPower )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalPowerScalar(infinity, result, argvals[0], (SCIP_Real)opdata.intval);

   return SCIP_OKAY;
}

/** curvature for EXPR_INTPOWER */
static
SCIP_DECL_EXPRCURV( exprcurvIntPower )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   *result = SCIPexprcurvPower(argbounds[0], argcurv[0], (SCIP_Real)opdata.intval);

   return SCIP_OKAY;
}

/** point evaluation for EXPR_SIGNPOWER */
static
SCIP_DECL_EXPREVAL( exprevalSignPower )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   if( argvals[0] > 0 )
      *result =  pow( argvals[0], opdata.dbl);
   else
      *result = -pow(-argvals[0], opdata.dbl);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_SIGNPOWER */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntSignPower )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalSignPowerScalar(infinity, result, argvals[0], opdata.dbl);

   return SCIP_OKAY;
}

/** curvature for EXPR_SIGNPOWER */
static
SCIP_DECL_EXPRCURV( exprcurvSignPower )
{   /*lint --e{715}*/
   SCIP_INTERVAL tmp;
   SCIP_EXPRCURV left;
   SCIP_EXPRCURV right;

   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   /* for x <= 0, signpower(x,c) = -(-x)^c
    * for x >= 0, signpower(x,c) =  ( x)^c
    *
    * thus, get curvatures for both parts and "intersect" them
    */

   if( argbounds[0].inf < 0 )
   {
      SCIPintervalSetBounds(&tmp, 0.0, -argbounds[0].inf);
      left = SCIPexprcurvNegate(SCIPexprcurvPower(tmp, SCIPexprcurvNegate(argcurv[0]), opdata.dbl));
   }
   else
   {
      left = SCIP_EXPRCURV_LINEAR;
   }

   if( argbounds[0].sup > 0 )
   {
      SCIPintervalSetBounds(&tmp, 0.0,  argbounds[0].sup);
      right = SCIPexprcurvPower(tmp, argcurv[0], opdata.dbl);
   }
   else
   {
      right = SCIP_EXPRCURV_LINEAR;
   }

   *result = (SCIP_EXPRCURV) (left & right);

   return SCIP_OKAY;
}

/** point evaluation for EXPR_EXP */
static
SCIP_DECL_EXPREVAL( exprevalExp )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = exp(argvals[0]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_EXP */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntExp )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalExp(infinity, result, argvals[0]);

   return SCIP_OKAY;
}

/** curvature for EXPR_EXP */
static
SCIP_DECL_EXPRCURV( exprcurvExp )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);

   /* expression is convex if child is convex
    * otherwise, we don't know
    */
   if( argcurv[0] & SCIP_EXPRCURV_CONVEX )
      *result = SCIP_EXPRCURV_CONVEX;
   else
      *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_LOG */
static
SCIP_DECL_EXPREVAL( exprevalLog )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = log(argvals[0]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_LOG */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntLog )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalLog(infinity, result, argvals[0]);

   return SCIP_OKAY;
}

/** curvature for EXPR_LOG */
static
SCIP_DECL_EXPRCURV( exprcurvLog )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);

   /* expression is concave if child is concave
    * otherwise, we don't know
    */
   if( argcurv[0] & SCIP_EXPRCURV_CONCAVE )
      *result = SCIP_EXPRCURV_CONCAVE;
   else
      *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_SIN */
static
SCIP_DECL_EXPREVAL( exprevalSin )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = sin(argvals[0]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_SIN */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntSin )
{   /*lint --e{715}*/
   assert(result != NULL);
   assert(argvals != NULL);
   assert(nargs == 1);

   SCIPintervalSin(infinity, result, *argvals);

   return SCIP_OKAY;
}

/* @todo implement exprcurvSin */
#define exprcurvSin exprcurvDefault

/** point evaluation for EXPR_COS */
static
SCIP_DECL_EXPREVAL( exprevalCos )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = cos(argvals[0]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_COS */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntCos )
{   /*lint --e{715}*/
   assert(result != NULL);
   assert(argvals != NULL);
   assert(nargs == 1);

   SCIPintervalCos(infinity, result, *argvals);

   return SCIP_OKAY;
}

/* @todo implement exprcurvCos */
#define exprcurvCos exprcurvDefault

/** point evaluation for EXPR_TAN */
static
SCIP_DECL_EXPREVAL( exprevalTan )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = tan(argvals[0]);

   return SCIP_OKAY;
}

/* @todo implement SCIPintervalTan */
#define exprevalIntTan exprevalIntDefault

/* @todo implement exprcurvTan */
#define exprcurvTan exprcurvDefault

/* erf and erfi do not seem to exist on every system, and we cannot really handle them anyway, so they are currently disabled */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_EXPREVAL( exprevalErf )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = erf(argvals[0]);

   return SCIP_OKAY;
}

/* @todo implement SCIPintervalErf */
#define exprevalIntErf exprevalIntDefault

/* @todo implement SCIPintervalErf */
#define exprcurvErf exprcurvDefault

static
SCIP_DECL_EXPREVAL( exprevalErfi )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   /* @TODO implement erfi evaluation */
   SCIPerrorMessage("erfi not implemented");

   return SCIP_ERROR;
}

/* @todo implement SCIPintervalErfi */
#define exprevalIntErfi NULL

#define exprcurvErfi exprcurvDefault
#endif

/** point evaluation for EXPR_MIN */
static
SCIP_DECL_EXPREVAL( exprevalMin )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = MIN(argvals[0], argvals[1]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_MIN */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntMin )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalMin(infinity, result, argvals[0], argvals[1]);

   return SCIP_OKAY;
}

/** curvature for EXPR_MIN */
static
SCIP_DECL_EXPRCURV( exprcurvMin )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argcurv != NULL);

   /* the minimum of two concave functions is concave
    * otherwise, we don't know
    */

   if( (argcurv[0] & SCIP_EXPRCURV_CONCAVE) && (argcurv[1] & SCIP_EXPRCURV_CONCAVE) )
      *result = SCIP_EXPRCURV_CONCAVE;
   else
      *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_MAX */
static
SCIP_DECL_EXPREVAL( exprevalMax )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = MAX(argvals[0], argvals[1]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_MAX */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntMax )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalMax(infinity, result, argvals[0], argvals[1]);

   return SCIP_OKAY;
}

/** curvature for EXPR_MAX */
static
SCIP_DECL_EXPRCURV( exprcurvMax )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argcurv != NULL);

   /* the maximum of two convex functions is convex
    * otherwise, we don't know
    */
   if( (argcurv[0] & SCIP_EXPRCURV_CONVEX) && (argcurv[1] & SCIP_EXPRCURV_CONVEX) )
      *result = SCIP_EXPRCURV_CONVEX;
   else
      *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_ABS */
static
SCIP_DECL_EXPREVAL( exprevalAbs )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = ABS(argvals[0]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_ABS */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntAbs )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalAbs(infinity, result, argvals[0]);

   return SCIP_OKAY;
}

/** curvature for EXPR_ABS */
static
SCIP_DECL_EXPRCURV( exprcurvAbs )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   /* if child is only negative, then abs(child) = -child
    * if child is only positive, then abs(child) = child
    * if child is both positive and negative, but also linear, then abs(child) is convex
    * otherwise, we don't know
    */
   if( argbounds[0].sup <= 0.0 )
      *result = SCIPexprcurvMultiply(-1.0, argcurv[0]);
   else if( argbounds[0].inf >= 0.0 )
      *result = argcurv[0];
   else if( argcurv[0] == SCIP_EXPRCURV_LINEAR )
      *result = SCIP_EXPRCURV_CONVEX;
   else
      *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_SIGN */
static
SCIP_DECL_EXPREVAL( exprevalSign )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   *result = SIGN(argvals[0]);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_SIGN */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntSign )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalSign(infinity, result, argvals[0]);

   return SCIP_OKAY;
}

/** curvature for EXPR_SIGN */
static
SCIP_DECL_EXPRCURV( exprcurvSign )
{   /*lint --e{715}*/
   assert(result    != NULL);
   assert(argbounds != NULL);

   /* if sign of child is clear, then sign is linear otherwise, we don't know */
   if( argbounds[0].sup <= 0.0 || argbounds[0].inf >= 0.0 )
      *result = SCIP_EXPRCURV_LINEAR;
   else
      *result = SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** point evaluation for EXPR_SUM */
static
SCIP_DECL_EXPREVAL( exprevalSum )
{   /*lint --e{715}*/
   int i;

   assert(result  != NULL);
   assert(argvals != NULL);

   *result = 0.0;
   for( i = 0; i < nargs; ++i )
      *result += argvals[i];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_SUM */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntSum )
{   /*lint --e{715}*/
   int i;

   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalSet(result, 0.0);

   for( i = 0; i < nargs; ++i )
      SCIPintervalAdd(infinity, result, *result, argvals[i]);

   return SCIP_OKAY;
}

/** curvature for EXPR_SUM */
static
SCIP_DECL_EXPRCURV( exprcurvSum )
{   /*lint --e{715}*/
   int i;

   assert(result  != NULL);
   assert(argcurv != NULL);

   *result = SCIP_EXPRCURV_LINEAR;

   for( i = 0; i < nargs; ++i )
      *result = SCIPexprcurvAdd(*result, argcurv[i]);

   return SCIP_OKAY;
}

/** point evaluation for EXPR_PRODUCT */
static
SCIP_DECL_EXPREVAL( exprevalProduct )
{   /*lint --e{715}*/
   int i;

   assert(result  != NULL);
   assert(argvals != NULL);

   *result = 1.0;
   for( i = 0; i < nargs; ++i )
      *result *= argvals[i];

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_PRODUCT */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntProduct )
{   /*lint --e{715}*/
   int i;

   assert(result  != NULL);
   assert(argvals != NULL);

   SCIPintervalSet(result, 1.0);

   for( i = 0; i < nargs; ++i )
      SCIPintervalMul(infinity, result, *result, argvals[i]);

   return SCIP_OKAY;
}

/** curvature for EXPR_PRODUCT */
static
SCIP_DECL_EXPRCURV( exprcurvProduct )
{   /*lint --e{715}*/
   SCIP_Bool hadnonconst;
   SCIP_Real constants;
   int i;

   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   /* if all factors are constant, then product is linear (even constant)
    * if only one factor is not constant, then product is curvature of this factor, multiplied by sign of product of remaining factors
    */
   *result = SCIP_EXPRCURV_LINEAR;
   hadnonconst = FALSE;
   constants = 1.0;

   for( i = 0; i < nargs; ++i )
   {
      if( argbounds[i].inf == argbounds[i].sup )  /*lint !e777*/
      {
         constants *= argbounds[i].inf;
      }
      else if( !hadnonconst )
      {
         /* first non-constant child */
         *result = argcurv[i];
         hadnonconst = TRUE;
      }
      else
      {
         /* more than one non-constant child, thus don't know curvature */
         *result = SCIP_EXPRCURV_UNKNOWN;
         break;
      }
   }

   *result = SCIPexprcurvMultiply(constants, *result);

   return SCIP_OKAY;
}

/** point evaluation for EXPR_LINEAR */
static
SCIP_DECL_EXPREVAL( exprevalLinear )
{   /*lint --e{715}*/
   SCIP_Real* coef;
   int i;

   assert(result  != NULL);
   assert(argvals != NULL || nargs == 0);
   assert(opdata.data != NULL);

   coef = &((SCIP_Real*)opdata.data)[nargs];

   *result = *coef;
   for( i = nargs-1, --coef; i >= 0; --i, --coef )
      *result += *coef * argvals[i];  /*lint !e613*/

   assert(++coef == (SCIP_Real*)opdata.data);

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_LINEAR */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntLinear )
{   /*lint --e{715}*/
   assert(result  != NULL);
   assert(argvals != NULL || nargs == 0);
   assert(opdata.data != NULL);

   SCIPintervalScalprodScalars(infinity, result, nargs, argvals, (SCIP_Real*)opdata.data);
   SCIPintervalAddScalar(infinity, result, *result, ((SCIP_Real*)opdata.data)[nargs]);

   return SCIP_OKAY;
}

/** curvature for EXPR_LINEAR */
static
SCIP_DECL_EXPRCURV( exprcurvLinear )
{   /*lint --e{715}*/
   SCIP_Real* data;
   int i;

   assert(result  != NULL);
   assert(argcurv != NULL);

   data = (SCIP_Real*)opdata.data;
   assert(data != NULL);

   *result = SCIP_EXPRCURV_LINEAR;

   for( i = 0; i < nargs; ++i )
      *result = SCIPexprcurvAdd(*result, SCIPexprcurvMultiply(data[i], argcurv[i]));

   return SCIP_OKAY;
}

/** expression data copy for EXPR_LINEAR */
static
SCIP_DECL_EXPRCOPYDATA( exprCopyDataLinear )
{  /*lint --e{715}*/
   SCIP_Real* targetdata;

   assert(blkmem != NULL);
   assert(nchildren >= 0);
   assert(opdatatarget != NULL);

   /* for a linear expression, we need to copy the array that holds the coefficients and constant term */
   assert(opdatasource.data != NULL);
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &targetdata, (SCIP_Real*)opdatasource.data, nchildren + 1) );  /*lint !e866*/
   opdatatarget->data = targetdata;

   return SCIP_OKAY;
}

/** expression data free for EXPR_LINEAR */
static
SCIP_DECL_EXPRFREEDATA( exprFreeDataLinear )
{  /*lint --e{715}*/
   SCIP_Real* freedata;

   assert(blkmem != NULL);
   assert(nchildren >= 0);

   freedata = (SCIP_Real*)opdata.data;
   assert(freedata != NULL);

   BMSfreeBlockMemoryArray(blkmem, &freedata, nchildren + 1);  /*lint !e866*/
}

/** point evaluation for EXPR_QUADRATIC */
static
SCIP_DECL_EXPREVAL( exprevalQuadratic )
{   /*lint --e{715}*/
   SCIP_EXPRDATA_QUADRATIC* quaddata;
   SCIP_Real* lincoefs;
   SCIP_QUADELEM* quadelems;
   int nquadelems;
   int i;

   assert(result  != NULL);
   assert(argvals != NULL || nargs == 0);

   quaddata = (SCIP_EXPRDATA_QUADRATIC*)opdata.data;
   assert(quaddata != NULL);

   lincoefs   = quaddata->lincoefs;
   nquadelems = quaddata->nquadelems;
   quadelems  = quaddata->quadelems;

   assert(quadelems != NULL || nquadelems == 0);
   assert(argvals != NULL || nquadelems == 0);

   *result = quaddata->constant;

   if( lincoefs != NULL )
   {
      for( i = nargs-1; i >= 0; --i )
         *result += lincoefs[i] * argvals[i];  /*lint !e613*/
   }

   for( i = 0; i < nquadelems; ++i, ++quadelems )  /*lint !e613*/
      *result += quadelems->coef * argvals[quadelems->idx1] * argvals[quadelems->idx2];  /*lint !e613*/

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_QUADRATIC */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntQuadratic )
{   /*lint --e{715}*/
   SCIP_EXPRDATA_QUADRATIC* quaddata;
   SCIP_Real* lincoefs;
   SCIP_QUADELEM* quadelems;
   int nquadelems;
   int i;
   int argidx;
   SCIP_Real sqrcoef;
   SCIP_INTERVAL lincoef;
   SCIP_INTERVAL tmp;

   assert(result  != NULL);
   assert(argvals != NULL || nargs == 0);

   quaddata = (SCIP_EXPRDATA_QUADRATIC*)opdata.data;
   assert(quaddata != NULL);

   lincoefs   = quaddata->lincoefs;
   nquadelems = quaddata->nquadelems;
   quadelems  = quaddata->quadelems;

   assert(quadelems != NULL || nquadelems == 0);
   assert(argvals   != NULL || nargs == 0);

   /* something fast for case of only one child */
   if( nargs == 1 )
   {
      SCIPintervalSet(&lincoef, lincoefs != NULL ? lincoefs[0] : 0.0);

      sqrcoef = 0.0;
      for( i = 0; i < nquadelems; ++i )
      {
         assert(quadelems[i].idx1 == 0);  /*lint !e613*/
         assert(quadelems[i].idx2 == 0);  /*lint !e613*/
         sqrcoef += quadelems[i].coef;    /*lint !e613*/
      }

      SCIPintervalQuad(infinity, result, sqrcoef, lincoef, argvals[0]);  /*lint !e613*/
      SCIPintervalAddScalar(infinity, result, *result, quaddata->constant);

      return SCIP_OKAY;
   }

   if( nargs == 2 && nquadelems > 0 )
   {
      /* if it's a bivariate quadratic expression with bilinear term, do something special */
      SCIP_Real ax;  /* square coefficient of first  child */
      SCIP_Real ay;  /* square coefficient of second child */
      SCIP_Real axy; /* bilinear coefficient */

      ax = 0.0;
      ay = 0.0;
      axy = 0.0;
      for( i = 0; i < nquadelems; ++i )
         if( quadelems[i].idx1 == 0 && quadelems[i].idx2 == 0 )       /*lint !e613*/
            ax += quadelems[i].coef;                                  /*lint !e613*/
         else if( quadelems[i].idx1 == 1 && quadelems[i].idx2 == 1 )  /*lint !e613*/
            ay += quadelems[i].coef;                                  /*lint !e613*/
         else
            axy += quadelems[i].coef;                                 /*lint !e613*/

      SCIPintervalQuadBivar(infinity, result, ax, ay, axy,
         lincoefs != NULL ? lincoefs[0] : 0.0, lincoefs != NULL ? lincoefs[1] : 0.0,
         argvals[0], argvals[1]);                                     /*lint !e613*/
      SCIPdebugMessage("%g x^2 + %g y^2 + %g x y + %g x + %g y = [%g,%g] for x = [%g,%g], y = [%g,%g]\n",
         ax, ay, axy, lincoefs != NULL ? lincoefs[0] : 0.0, lincoefs != NULL ? lincoefs[1] : 0.0,
         result->inf, result->sup, argvals[0].inf, argvals[0].sup, argvals[1].inf, argvals[1].sup);  /*lint !e613*/

      SCIPintervalAddScalar(infinity, result, *result, quaddata->constant);

      return SCIP_OKAY;
   }

   /* make sure coefficients are sorted */
   quadraticdataSort(quaddata);

   SCIPintervalSet(result, quaddata->constant);

   /* for each argument, we collect it's linear index from lincoefs, it's square coefficients and all factors from bilinear terms
    * then we compute the interval sqrcoef*x^2 + lincoef*x and add it to result
    * @todo split quadratic expression into bivariate quadratic terms and apply the above method
    */
   i = 0;
   for( argidx = 0; argidx < nargs; ++argidx )
   {
      if( i == nquadelems || quadelems[i].idx1 > argidx )  /*lint !e613*/
      {
         /* there are no quadratic terms with argidx in its first argument, that should be easy to handle */
         if( lincoefs != NULL )
         {
            SCIPintervalMulScalar(infinity, &tmp, argvals[argidx], lincoefs[argidx]);  /*lint !e613*/
            SCIPintervalAdd(infinity, result, *result, tmp);
         }
         continue;
      }

      sqrcoef = 0.0;
      SCIPintervalSet(&lincoef, lincoefs != NULL ? lincoefs[argidx] : 0.0);

      assert(i < nquadelems && quadelems[i].idx1 == argidx);  /*lint !e613*/
      do
      {
         if( quadelems[i].idx2 == argidx )  /*lint !e613*/
         {
            sqrcoef += quadelems[i].coef;   /*lint !e613*/
         }
         else
         {
            SCIPintervalMulScalar(infinity, &tmp, argvals[quadelems[i].idx2], quadelems[i].coef);  /*lint !e613*/
            SCIPintervalAdd(infinity, &lincoef, lincoef, tmp);
         }
         ++i;
      }
      while( i < nquadelems && quadelems[i].idx1 == argidx );  /*lint !e613*/
      assert(i == nquadelems || quadelems[i].idx1 > argidx);   /*lint !e613*/

      SCIPintervalQuad(infinity, &tmp, sqrcoef, lincoef, argvals[argidx]);  /*lint !e613*/
      SCIPintervalAdd(infinity, result, *result, tmp);
   }
   assert(i == nquadelems);

   return SCIP_OKAY;
}

/** curvature for EXPR_QUADRATIC */
static
SCIP_DECL_EXPRCURV( exprcurvQuadratic )
{   /*lint --e{715}*/
   SCIP_EXPRDATA_QUADRATIC* data;
   SCIP_QUADELEM* quadelems;
   int nquadelems;
   SCIP_Real* lincoefs;
   int i;

   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   data = (SCIP_EXPRDATA_QUADRATIC*)opdata.data;
   assert(data != NULL);

   lincoefs   = data->lincoefs;
   quadelems  = data->quadelems;
   nquadelems = data->nquadelems;

   *result = SCIP_EXPRCURV_LINEAR;

   if( lincoefs != NULL )
      for( i = 0; i < nargs; ++i )
         *result = SCIPexprcurvAdd(*result, SCIPexprcurvMultiply(lincoefs[i], argcurv[i]));

   /* @todo could try cholesky factorization if all children linear...
    * @todo should then cache the result
    */
   for( i = 0; i < nquadelems && *result != SCIP_EXPRCURV_UNKNOWN; ++i )
   {
      if( quadelems[i].coef == 0.0 )
         continue;

      if( argbounds[quadelems[i].idx1].inf == argbounds[quadelems[i].idx1].sup &&  /*lint !e777*/
         +argbounds[quadelems[i].idx2].inf == argbounds[quadelems[i].idx2].sup
         )  /*lint !e777*/
      {
         /* both factors are constants -> curvature does not change */
         continue;
      }

      if( argbounds[quadelems[i].idx1].inf == argbounds[quadelems[i].idx1].sup )  /*lint !e777*/
      {
         /* first factor is constant, second is not -> add curvature of second */
         *result = SCIPexprcurvAdd(*result, SCIPexprcurvMultiply(quadelems[i].coef * argbounds[quadelems[i].idx1].inf, argcurv[quadelems[i].idx2]));
      }
      else if( argbounds[quadelems[i].idx2].inf == argbounds[quadelems[i].idx2].sup )  /*lint !e777*/
      {
         /* first factor is not constant, second is -> add curvature of first */
         *result = SCIPexprcurvAdd(*result, SCIPexprcurvMultiply(quadelems[i].coef * argbounds[quadelems[i].idx2].inf, argcurv[quadelems[i].idx1]));
      }
      else if( quadelems[i].idx1 == quadelems[i].idx2 )
      {
         /* both factors not constant, but the same (square term) */
         *result = SCIPexprcurvAdd(*result, SCIPexprcurvMultiply(quadelems[i].coef, SCIPexprcurvPower(argbounds[quadelems[i].idx1], argcurv[quadelems[i].idx1], 2.0)));
      }
      else
      {
         /* two different non-constant factors -> can't tell about curvature */
         *result = SCIP_EXPRCURV_UNKNOWN;
      }
   }

   return SCIP_OKAY;
}

/** expression data copy for EXPR_QUADRATIC */
static
SCIP_DECL_EXPRCOPYDATA( exprCopyDataQuadratic )
{  /*lint --e{715}*/
   SCIP_EXPRDATA_QUADRATIC* sourcedata;

   assert(blkmem != NULL);
   assert(opdatatarget != NULL);

   sourcedata = (SCIP_EXPRDATA_QUADRATIC*)opdatasource.data;
   assert(sourcedata != NULL);

   SCIP_CALL( quadraticdataCreate(blkmem, (SCIP_EXPRDATA_QUADRATIC**)&opdatatarget->data,
         sourcedata->constant, nchildren, sourcedata->lincoefs, sourcedata->nquadelems, sourcedata->quadelems) );

   return SCIP_OKAY;
}

/** expression data free for EXPR_QUADRATIC */
static
SCIP_DECL_EXPRFREEDATA( exprFreeDataQuadratic )
{  /*lint --e{715}*/
   SCIP_EXPRDATA_QUADRATIC* quadraticdata;

   assert(blkmem != NULL);
   assert(nchildren >= 0);

   quadraticdata = (SCIP_EXPRDATA_QUADRATIC*)opdata.data;
   assert(quadraticdata != NULL);

   if( quadraticdata->lincoefs != NULL )
   {
      BMSfreeBlockMemoryArray(blkmem, &quadraticdata->lincoefs, nchildren);
   }

   if( quadraticdata->nquadelems > 0 )
   {
      assert(quadraticdata->quadelems != NULL);
      BMSfreeBlockMemoryArray(blkmem, &quadraticdata->quadelems, quadraticdata->nquadelems);
   }

   BMSfreeBlockMemory(blkmem, &quadraticdata);
}

/** point evaluation for EXPR_POLYNOMIAL */
static
SCIP_DECL_EXPREVAL( exprevalPolynomial )
{   /*lint --e{715}*/
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
   SCIP_EXPRDATA_MONOMIAL*   monomialdata;
   SCIP_Real childval;
   SCIP_Real exponent;
   SCIP_Real monomialval;
   int i;
   int j;

   assert(result != NULL);
   assert(argvals != NULL || nargs == 0);
   assert(opdata.data != NULL);

   polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)opdata.data;
   assert(polynomialdata != NULL);

   *result = polynomialdata->constant;

   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      monomialdata = polynomialdata->monomials[i];
      assert(monomialdata != NULL);

      monomialval = monomialdata->coef;
      for( j = 0; j < monomialdata->nfactors; ++j )
      {
         assert(monomialdata->childidxs[j] >= 0);
         assert(monomialdata->childidxs[j] < nargs);

         childval = argvals[monomialdata->childidxs[j]];  /*lint !e613*/
         if( childval == 1.0 )  /* 1^anything == 1 */
            continue;

         exponent = monomialdata->exponents[j];

         if( childval == 0.0 )
         {
            if( exponent > 0.0 )
            {
               /* 0^positive == 0 */
               monomialval = 0.0;
               break;
            }
            else if( exponent < 0.0 )
            {
               /* 0^negative = nan (or should it be +inf?, doesn't really matter) */
#ifdef NAN
               *result = NAN;
#else
               /* cppcheck-suppress wrongmathcall */
               *result = pow(0.0, -1.0);
#endif
               return SCIP_OKAY;
            }
            /* 0^0 == 1 */
            continue;
         }

         /* cover some special exponents separately to avoid calling expensive pow function */
         if( exponent == 0.0 )
            continue;
         if( exponent == 1.0 )
         {
            monomialval *= childval;
            continue;
         }
         if( exponent == 2.0 )
         {
            monomialval *= childval * childval;
            continue;
         }
         if( exponent == 0.5 )
         {
            monomialval *= sqrt(childval);
            continue;
         }
         if( exponent == -1.0 )
         {
            monomialval /= childval;
            continue;
         }
         if( exponent == -2.0 )
         {
            monomialval /= childval * childval;
            continue;
         }
         monomialval *= pow(childval, exponent);
      }

      *result += monomialval;
   }

   return SCIP_OKAY;
}

/** interval evaluation for EXPR_POLYNOMIAL */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntPolynomial )
{   /*lint --e{715}*/
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
   SCIP_EXPRDATA_MONOMIAL*   monomialdata;
   SCIP_INTERVAL childval;
   SCIP_INTERVAL monomialval;
   SCIP_Real exponent;
   int i;
   int j;

   assert(result != NULL);
   assert(argvals != NULL || nargs == 0);
   assert(opdata.data != NULL);

   polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)opdata.data;
   assert(polynomialdata != NULL);

   SCIPintervalSet(result, polynomialdata->constant);

   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      monomialdata = polynomialdata->monomials[i];
      assert(monomialdata != NULL);

      SCIPintervalSet(&monomialval, monomialdata->coef);
      for( j = 0; j < monomialdata->nfactors && !SCIPintervalIsEntire(infinity, monomialval); ++j )
      {
         assert(monomialdata->childidxs[j] >= 0);
         assert(monomialdata->childidxs[j] < nargs);

         childval = argvals[monomialdata->childidxs[j]];  /*lint !e613*/

         exponent = monomialdata->exponents[j];

         /* cover some special exponents separately to avoid calling expensive pow function */
         if( exponent == 0.0 )
            continue;

         if( exponent == 1.0 )
         {
            SCIPintervalMul(infinity, &monomialval, monomialval, childval);
            continue;
         }

         if( exponent == 2.0 )
         {
            SCIPintervalSquare(infinity, &childval, childval);
            SCIPintervalMul(infinity, &monomialval, monomialval, childval);
            continue;
         }

         if( exponent == 0.5 )
         {
            SCIPintervalSquareRoot(infinity, &childval, childval);
            if( SCIPintervalIsEmpty(infinity, childval) )
            {
               SCIPintervalSetEmpty(result);
               break;
            }
            SCIPintervalMul(infinity, &monomialval, monomialval, childval);
            continue;
         }
         else if( exponent == -1.0 )
         {
            SCIPintervalDiv(infinity, &monomialval, monomialval, childval);
         }
         else if( exponent == -2.0 )
         {
            SCIPintervalSquare(infinity, &childval, childval);
            SCIPintervalDiv(infinity, &monomialval, monomialval, childval);
         }
         else
         {
            SCIPintervalPowerScalar(infinity, &childval, childval, exponent);
            if( SCIPintervalIsEmpty(infinity, childval) )
            {
               SCIPintervalSetEmpty(result);
               return SCIP_OKAY;
            }
            SCIPintervalMul(infinity, &monomialval, monomialval, childval);
         }

         /* the cases in which monomialval gets empty should have been catched */
         assert(!SCIPintervalIsEmpty(infinity, monomialval));
      }

      SCIPintervalAdd(infinity, result, *result, monomialval);
   }

   return SCIP_OKAY;
}

/** curvature for EXPR_POLYNOMIAL */
static
SCIP_DECL_EXPRCURV( exprcurvPolynomial )
{   /*lint --e{715}*/
   SCIP_EXPRDATA_POLYNOMIAL* data;
   SCIP_EXPRDATA_MONOMIAL** monomials;
   SCIP_EXPRDATA_MONOMIAL* monomial;
   int nmonomials;
   int i;

   assert(result    != NULL);
   assert(argcurv   != NULL);
   assert(argbounds != NULL);

   data = (SCIP_EXPRDATA_POLYNOMIAL*)opdata.data;
   assert(data != NULL);

   monomials  = data->monomials;
   nmonomials = data->nmonomials;

   *result = SCIP_EXPRCURV_LINEAR;

   for( i = 0; i < nmonomials && *result != SCIP_EXPRCURV_UNKNOWN; ++i )
   {
      /* we assume that some simplifier was running, so that monomials do not have constants in their factors and such that all factors are different
       * (result would still be correct)
       */
      monomial = monomials[i];
      *result = SCIPexprcurvAdd(*result, SCIPexprcurvMultiply(monomial->coef, SCIPexprcurvMonomial(monomial->nfactors, monomial->exponents, monomial->childidxs, argcurv, argbounds)));
   }

   return SCIP_OKAY;
}

/** expression data copy for EXPR_POLYNOMIAL */
static
SCIP_DECL_EXPRCOPYDATA( exprCopyDataPolynomial )
{  /*lint --e{715}*/
   SCIP_EXPRDATA_POLYNOMIAL* sourcepolynomialdata;
   SCIP_EXPRDATA_POLYNOMIAL* targetpolynomialdata;

   assert(blkmem != NULL);
   assert(opdatatarget != NULL);

   sourcepolynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)opdatasource.data;
   assert(sourcepolynomialdata != NULL);

   SCIP_CALL( polynomialdataCopy(blkmem, &targetpolynomialdata, sourcepolynomialdata) );

   opdatatarget->data = (void*)targetpolynomialdata;

   return SCIP_OKAY;
}

/** expression data free for EXPR_POLYNOMIAL */
static
SCIP_DECL_EXPRFREEDATA( exprFreeDataPolynomial )
{  /*lint --e{715}*/
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;

   assert(blkmem != NULL);

   polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)opdata.data;
   assert(polynomialdata != NULL);

   polynomialdataFree(blkmem, &polynomialdata);
}

/** point evaluation for user expression */
static
SCIP_DECL_EXPREVAL( exprevalUser )
{  /*lint --e{715}*/
   SCIP_EXPRDATA_USER* exprdata;

   exprdata = (SCIP_EXPRDATA_USER*) opdata.data;

   SCIP_CALL( exprdata->eval(exprdata->userdata, nargs, argvals, result, NULL, NULL) );

   return SCIP_OKAY;
}

/** interval evaluation for user expression */
static
SCIP_DECL_EXPRINTEVAL( exprevalIntUser )
{  /*lint --e{715}*/
   SCIP_EXPRDATA_USER* exprdata;

   exprdata = (SCIP_EXPRDATA_USER*) opdata.data;

   if( exprdata->inteval != NULL )
   {
      SCIP_CALL( exprdata->inteval(infinity, exprdata->userdata, nargs, argvals, result, NULL, NULL) );
   }
   else
   {
      /* if user does not provide interval evaluation, then return a result that is always correct */
      SCIPintervalSetEntire(infinity, result);
   }

   return SCIP_OKAY;
}

/** curvature check for user expression */
static
SCIP_DECL_EXPRCURV( exprcurvUser )
{
   SCIP_EXPRDATA_USER* exprdata;

   exprdata = (SCIP_EXPRDATA_USER*) opdata.data;

   if( exprdata->curv != NULL )
   {
      SCIP_CALL( exprdata->curv(infinity, exprdata->userdata, nargs, argbounds, argcurv, result) );
   }
   else
   {
      /* if user does not provide curvature check, then return unknown (which is handled like indefinite) */
      *result = SCIP_EXPRCURV_UNKNOWN;
   }

   return SCIP_OKAY;
}

/** data copy for user expression */
static
SCIP_DECL_EXPRCOPYDATA( exprCopyDataUser )
{
   SCIP_EXPRDATA_USER* exprdatasource;
   SCIP_EXPRDATA_USER* exprdatatarget;

   assert(blkmem != NULL);
   assert(opdatatarget != NULL);

   exprdatasource = (SCIP_EXPRDATA_USER*)opdatasource.data;
   assert(exprdatasource != NULL);

   /* duplicate expression data */
   SCIP_ALLOC( BMSduplicateBlockMemory(blkmem, &exprdatatarget, exprdatasource) );

   /* duplicate user expression data, if any */
   if( exprdatasource->copydata != NULL )
   {
      SCIP_CALL( exprdatasource->copydata(blkmem, nchildren, exprdatasource->userdata, &exprdatatarget->userdata) );
   }
   else
   {
      /* if no copy function for data, then there has to be no data */
      assert(exprdatatarget->userdata == NULL);
   }

   opdatatarget->data = (void*)exprdatatarget;

   return SCIP_OKAY;
}

/** data free for user expression */
static
SCIP_DECL_EXPRFREEDATA( exprFreeDataUser )
{
   SCIP_EXPRDATA_USER* exprdata;

   assert(blkmem != NULL);

   exprdata = (SCIP_EXPRDATA_USER*)opdata.data;

   /* free user expression data, if any */
   if( exprdata->freedata != NULL )
   {
      exprdata->freedata(blkmem, nchildren, exprdata->userdata);
   }
   else
   {
      assert(exprdata->userdata == NULL);
   }

   /* free expression data */
   BMSfreeBlockMemory(blkmem, &exprdata);
}

/** element in table of expression operands */
struct exprOpTableElement
{
   const char*           name;               /**< name of operand (used for printing) */
   int                   nargs;              /**< number of arguments (negative if not fixed) */
   SCIP_DECL_EXPREVAL    ((*eval));          /**< evaluation function */
   SCIP_DECL_EXPRINTEVAL ((*inteval));       /**< interval evaluation function */
   SCIP_DECL_EXPRCURV    ((*curv));          /**< curvature check function */
   SCIP_DECL_EXPRCOPYDATA ((*copydata));     /**< expression data copy function, or NULL to only opdata union */
   SCIP_DECL_EXPRFREEDATA ((*freedata));     /**< expression data free function, or NULL if nothing to free */
};

#define EXPROPEMPTY {NULL, -1, NULL, NULL, NULL, NULL, NULL}

/** table containing for each operand the name, the number of children, and some evaluation functions */
static
struct exprOpTableElement exprOpTable[] =
   {
      EXPROPEMPTY,
      { "variable",          0, exprevalVar,        exprevalIntVar,        exprcurvVar,        NULL, NULL  },
      { "constant",          0, exprevalConst,      exprevalIntConst,      exprcurvConst,      NULL, NULL  },
      { "parameter",         0, exprevalParam,      exprevalIntParam,      exprcurvParam,      NULL, NULL  },
      EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY,
      { "plus",              2, exprevalPlus,       exprevalIntPlus,       exprcurvPlus,       NULL, NULL  },
      { "minus",             2, exprevalMinus,      exprevalIntMinus,      exprcurvMinus,      NULL, NULL  },
      { "mul",               2, exprevalMult,       exprevalIntMult,       exprcurvMult,       NULL, NULL  },
      { "div",               2, exprevalDiv,        exprevalIntDiv,        exprcurvDiv,        NULL, NULL  },
      { "sqr",               1, exprevalSquare,     exprevalIntSquare,     exprcurvSquare,     NULL, NULL  },
      { "sqrt",              1, exprevalSquareRoot, exprevalIntSquareRoot, exprcurvSquareRoot, NULL, NULL  },
      { "realpower",         1, exprevalRealPower,  exprevalIntRealPower,  exprcurvRealPower,  NULL, NULL  },
      { "intpower",          1, exprevalIntPower,   exprevalIntIntPower,   exprcurvIntPower,   NULL, NULL  },
      { "signpower",         1, exprevalSignPower,  exprevalIntSignPower,  exprcurvSignPower,  NULL, NULL  },
      { "exp",               1, exprevalExp,        exprevalIntExp,        exprcurvExp,        NULL, NULL  },
      { "log",               1, exprevalLog,        exprevalIntLog,        exprcurvLog,        NULL, NULL  },
      { "sin",               1, exprevalSin,        exprevalIntSin,        exprcurvSin,        NULL, NULL  },
      { "cos",               1, exprevalCos,        exprevalIntCos,        exprcurvCos,        NULL, NULL  },
      { "tan",               1, exprevalTan,        exprevalIntTan,        exprcurvTan,        NULL, NULL  },
      /* { "erf",               1, exprevalErf,        exprevalIntErf,        exprcurvErf,        NULL, NULL  }, */
      /* { "erfi",              1, exprevalErfi,       exprevalIntErfi        exprcurvErfi,       NULL, NULL  }, */
      EXPROPEMPTY, EXPROPEMPTY,
      { "min",               2, exprevalMin,        exprevalIntMin,        exprcurvMin,        NULL, NULL  },
      { "max",               2, exprevalMax,        exprevalIntMax,        exprcurvMax,        NULL, NULL  },
      { "abs",               1, exprevalAbs,        exprevalIntAbs,        exprcurvAbs,        NULL, NULL  },
      { "sign",              1, exprevalSign,       exprevalIntSign,       exprcurvSign,       NULL, NULL  },
      EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY,
      EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY,
      EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY,
      EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY,
      EXPROPEMPTY, EXPROPEMPTY, EXPROPEMPTY,
      { "sum",              -2, exprevalSum,        exprevalIntSum,        exprcurvSum,        NULL, NULL  },
      { "prod",             -2, exprevalProduct,    exprevalIntProduct,    exprcurvProduct,    NULL, NULL  },
      { "linear",           -2, exprevalLinear,     exprevalIntLinear,     exprcurvLinear,     exprCopyDataLinear,     exprFreeDataLinear     },
      { "quadratic",        -2, exprevalQuadratic,  exprevalIntQuadratic,  exprcurvQuadratic,  exprCopyDataQuadratic,  exprFreeDataQuadratic  },
      { "polynomial",       -2, exprevalPolynomial, exprevalIntPolynomial, exprcurvPolynomial, exprCopyDataPolynomial, exprFreeDataPolynomial },
      { "user",             -2, exprevalUser,       exprevalIntUser,       exprcurvUser,       exprCopyDataUser,       exprFreeDataUser       }
   };

/**@} */

/**@name Expression operand methods */
/**@{ */

/** gives the name of an operand as string */
const char* SCIPexpropGetName(
   SCIP_EXPROP           op                  /**< expression operand */
   )
{
   assert(op < SCIP_EXPR_LAST);

   return exprOpTable[op].name;
}

/** gives the number of children of a simple operand */
int SCIPexpropGetNChildren(
   SCIP_EXPROP           op                  /**< expression operand */
   )
{
   assert(op < SCIP_EXPR_LAST);

   return exprOpTable[op].nargs;
}

/**@} */

/**@name Expressions private methods */
/**@{ */

/** creates an expression
 *
 *  Note, that the expression is allocated but for the children only the pointer is copied.
 */
static
SCIP_RETCODE exprCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   SCIP_EXPROP           op,                 /**< operand of expression */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children */
   SCIP_EXPROPDATA       opdata              /**< operand data */
   )
{
   assert(blkmem != NULL);
   assert(expr   != NULL);
   assert(children != NULL || nchildren == 0);
   assert(children == NULL || nchildren >  0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, expr) );

   (*expr)->op        = op;
   (*expr)->nchildren = nchildren;
   (*expr)->children  = children;
   (*expr)->data      = opdata;

   return SCIP_OKAY;
}

/** tries to convert a given (operator,operatordata) pair into a polynomial operator with corresponding data
 *
 *  Does not do this for constants.
 *  If conversion is not possible or operator is already polynomial, *op and *data are
 *  left untouched.
 */
static
SCIP_RETCODE exprConvertToPolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPROP*          op,                 /**< pointer to expression operator */
   SCIP_EXPROPDATA*      data,               /**< pointer to expression data */
   int                   nchildren           /**< number of children of operator */
   )
{
   assert(blkmem != NULL);
   assert(op != NULL);
   assert(data != NULL);

   switch( *op )
   {
   case SCIP_EXPR_VARIDX:
   case SCIP_EXPR_PARAM:
   case SCIP_EXPR_CONST:
      break;

   case SCIP_EXPR_PLUS:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomials[2];
      int childidx;
      SCIP_Real exponent;

      assert(nchildren == 2);

      /* create monomial for first child */
      childidx = 0;
      exponent = 1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomials[0], 1.0, 1, &childidx, &exponent) );

      /* create monomial for second child */
      childidx = 1;
      exponent = 1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomials[1], 1.0, 1, &childidx, &exponent) );

      /* create polynomial for sum of children */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 2, monomials, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_MINUS:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomials[2];
      int childidx;
      SCIP_Real exponent;

      assert(nchildren == 2);

      /* create monomial for first child */
      childidx = 0;
      exponent = 1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomials[0],  1.0, 1, &childidx, &exponent) );

      /* create monomial for second child */
      childidx = 1;
      exponent = 1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomials[1], -1.0, 1, &childidx, &exponent) );

      /* create polynomial for difference of children */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 2, monomials, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_MUL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx[2];
      SCIP_Real exponent[2];

      assert(nchildren == 2);

      /* create monomial for product of children */
      childidx[0] = 0;
      childidx[1] = 1;
      exponent[0] = 1.0;
      exponent[1] = 1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 2, childidx, exponent) );

      /* create polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 1, &monomial, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_DIV:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx[2];
      SCIP_Real exponent[2];

      assert(nchildren == 2);

      /* create monomial for division of children */
      childidx[0] = 0;
      childidx[1] = 1;
      exponent[0] =  1.0;
      exponent[1] = -1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 2, childidx, exponent) );

      /* create polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 1, &monomial, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_SQUARE:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx;
      SCIP_Real exponent;

      assert(nchildren == 1);

      /* create monomial for square of child */
      childidx = 0;
      exponent = 2.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 1, &childidx, &exponent) );

      /* create polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 1, &monomial, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_SQRT:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx;
      SCIP_Real exponent;

      assert(nchildren == 1);

      /* create monomial for square root of child */
      childidx = 0;
      exponent = 0.5;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 1, &childidx, &exponent) );

      /* create polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 1, &monomial, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_REALPOWER:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx;

      assert(nchildren == 1);

      /* convert to child0 to the power of exponent */

      /* create monomial for power of first child */
      childidx = 0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 1, &childidx, &data->dbl) );

      /* create polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 1, &monomial, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_SIGNPOWER:
   {
      SCIP_Real exponent;

      assert(nchildren == 1);

      /* check if exponent is an odd integer */
      exponent = data->dbl;
      if( EPSISINT(exponent, 0.0) && (int)exponent % 2 != 0 )  /*lint !e835*/
      {
         /* convert to child0 to the power of exponent, since sign is kept by taking power */
         SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
         SCIP_EXPRDATA_MONOMIAL* monomial;
         int childidx;

         /* create monomial for power of first child */
         childidx = 0;
         SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 1, &childidx, &exponent) );

         /* create polynomial */
         SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 1, &monomial, 0.0, FALSE) );

         *op = SCIP_EXPR_POLYNOMIAL;
         data->data = (void*)polynomialdata;
      }
      /* if exponent is not an odd integer constant, then keep it as signpower expression */
      break;
   }

   case SCIP_EXPR_INTPOWER:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx;
      SCIP_Real exponent;

      assert(nchildren == 1);

      /* create monomial for power of child */
      childidx = 0;
      exponent = data->intval;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 1, &childidx, &exponent) );

      /* create polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 1, &monomial, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   case SCIP_EXPR_USER:
      break;

   case SCIP_EXPR_SUM:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx;
      int i;
      SCIP_Real exponent;

      /* create empty polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 0, NULL, 0.0, FALSE) );
      SCIP_CALL( polynomialdataEnsureMonomialsSize(blkmem, polynomialdata, nchildren) );
      assert(polynomialdata->monomialssize >= nchildren);

      /* add summands as monomials */
      childidx = 0;
      exponent = 1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 1, &childidx, &exponent) );
      for( i = 0; i < nchildren; ++i )
      {
         monomial->childidxs[0] = i;
         SCIP_CALL( polynomialdataAddMonomials(blkmem, polynomialdata, 1, &monomial, TRUE) );
      }
      SCIPexprFreeMonomial(blkmem, &monomial);

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_PRODUCT:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx;
      int i;
      SCIP_Real exponent;

      /* create monomial */
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 0, NULL, NULL) );
      SCIP_CALL( monomialdataEnsureFactorsSize(blkmem, monomial, nchildren) );
      exponent = 1.0;
      for( i = 0; i < nchildren; ++i )
      {
         childidx = i;
         SCIP_CALL( SCIPexprAddMonomialFactors(blkmem, monomial, 1, &childidx, &exponent) );
      }

      /* create polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 1, &monomial, 0.0, FALSE) );

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_LINEAR:
   {
      SCIP_Real* lineardata;
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      int childidx;
      int i;
      SCIP_Real exponent;

      /* get coefficients of linear term */
      lineardata = (SCIP_Real*)data->data;
      assert(lineardata != NULL);

      /* create polynomial consisting of constant from linear term */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 0, NULL, lineardata[nchildren], FALSE) );
      /* ensure space for linear coefficients */
      SCIP_CALL( polynomialdataEnsureMonomialsSize(blkmem, polynomialdata, nchildren) );
      assert(polynomialdata->monomialssize >= nchildren);

      /* add summands as monomials */
      childidx = 0;
      exponent = 1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &monomial, 1.0, 1, &childidx, &exponent) );
      for( i = 0; i < nchildren; ++i )
      {
         monomial->coef = lineardata[i];
         monomial->childidxs[0] = i;
         SCIP_CALL( polynomialdataAddMonomials(blkmem, polynomialdata, 1, &monomial, TRUE) );
      }
      SCIPexprFreeMonomial(blkmem, &monomial);

      /* free linear expression data */
      exprFreeDataLinear(blkmem, nchildren, *data);

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* quaddata;
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* squaremonomial;
      SCIP_EXPRDATA_MONOMIAL* bilinmonomial;
      SCIP_EXPRDATA_MONOMIAL* linmonomial;
      int childidx[2];
      SCIP_Real exponent[2];
      int i;

      /* get data of quadratic expression */
      quaddata = (SCIP_EXPRDATA_QUADRATIC*)data->data;
      assert(quaddata != NULL);

      /* create empty polynomial */
      SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 0, NULL, quaddata->constant, FALSE) );
      /* ensure space for linear and quadratic terms */
      SCIP_CALL( polynomialdataEnsureMonomialsSize(blkmem, polynomialdata, (quaddata->lincoefs != NULL ? nchildren : 0) + quaddata->nquadelems) );
      assert(polynomialdata->monomialssize >= quaddata->nquadelems);

      childidx[0] = 0;
      childidx[1] = 0;

      /* create monomial templates */
      exponent[0] = 2.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &squaremonomial, 1.0, 1, childidx, exponent) );
      exponent[0] = 1.0;
      exponent[1] = 1.0;
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &bilinmonomial,  1.0, 2, childidx, exponent) );
      SCIP_CALL( SCIPexprCreateMonomial(blkmem, &linmonomial,    1.0, 1, childidx, exponent) );

      /* add linear terms as monomials */
      if( quaddata->lincoefs != NULL )
         for( i = 0; i < nchildren; ++i )
            if( quaddata->lincoefs[i] != 0.0 )
            {
               linmonomial->childidxs[0] = i;
               linmonomial->coef         = quaddata->lincoefs[i];
               SCIP_CALL( polynomialdataAddMonomials(blkmem, polynomialdata, 1, &linmonomial, TRUE) );
            }

      /* add quadratic terms as monomials */
      for( i = 0; i < quaddata->nquadelems; ++i )
      {
         if( quaddata->quadelems[i].idx1 == quaddata->quadelems[i].idx2 )
         {
            squaremonomial->childidxs[0] = quaddata->quadelems[i].idx1;
            squaremonomial->coef         = quaddata->quadelems[i].coef;
            SCIP_CALL( polynomialdataAddMonomials(blkmem, polynomialdata, 1, &squaremonomial, TRUE) );
         }
         else
         {
            bilinmonomial->childidxs[0] = quaddata->quadelems[i].idx1;
            bilinmonomial->childidxs[1] = quaddata->quadelems[i].idx2;
            bilinmonomial->coef         = quaddata->quadelems[i].coef;
            SCIP_CALL( polynomialdataAddMonomials(blkmem, polynomialdata, 1, &bilinmonomial, TRUE) );
         }
      }
      SCIPexprFreeMonomial(blkmem, &squaremonomial);
      SCIPexprFreeMonomial(blkmem, &bilinmonomial);
      SCIPexprFreeMonomial(blkmem, &linmonomial);

      /* free quadratic expression data */
      exprFreeDataQuadratic(blkmem, nchildren, *data);

      *op = SCIP_EXPR_POLYNOMIAL;
      data->data = (void*)polynomialdata;

      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   case SCIP_EXPR_LAST:
      break;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** converts polynomial expression back into simpler expression, if possible */
static
SCIP_RETCODE exprUnconvertPolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPROP*          op,                 /**< pointer to expression operator */
   SCIP_EXPROPDATA*      data,               /**< pointer to expression data holding polynomial data */
   int                   nchildren,          /**< number of children of operator */
   void**                children            /**< children array */
   )
{
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
   SCIP_EXPRDATA_MONOMIAL* monomial;
   int maxdegree;
   int nlinmonomials;
   int i;
   int j;

   assert(blkmem != NULL);
   assert(op != NULL);
   assert(*op == SCIP_EXPR_POLYNOMIAL);
   assert(data != NULL);
   assert(children != NULL || nchildren == 0);

   polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)data->data;
   assert(polynomialdata != NULL);

   /* make sure monomials are sorted and merged */
   polynomialdataMergeMonomials(blkmem, polynomialdata, 0.0, TRUE);

   /* if no monomials, then leave as it is */
   if( polynomialdata->nmonomials == 0 )
      return SCIP_OKAY;

   /* check maximal degree of polynomial only - not considering children expressions
    * check number of linear monomials */
   maxdegree = 0;
   nlinmonomials = 0;
   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      int monomialdegree;

      monomial = polynomialdata->monomials[i];
      assert(monomial != NULL);

      monomialdegree = 0;
      for(j = 0; j < monomial->nfactors; ++j )
      {
         if( !EPSISINT(monomial->exponents[j], 0.0) || monomial->exponents[j] < 0.0 )  /*lint !e835*/
         {
            monomialdegree = SCIP_EXPR_DEGREEINFINITY;
            break;
         }

         monomialdegree += (int)EPSROUND(monomial->exponents[j], 0.0);  /*lint !e835*/
      }

      if( monomialdegree == SCIP_EXPR_DEGREEINFINITY )
      {
         maxdegree = SCIP_EXPR_DEGREEINFINITY;
         break;
      }

      if( monomialdegree == 1 )
         ++nlinmonomials;

      if( monomialdegree > maxdegree )
         maxdegree = monomialdegree;
   }
   assert(maxdegree > 0 );

   if( maxdegree == 1 )
   {
      /* polynomial is a linear expression in children */

      /* polynomial simplification and monomial merging should ensure that monomial i corresponds to child i and that there are not unused children */
      assert(polynomialdata->nmonomials == nchildren);
      assert(polynomialdata->nmonomials == nlinmonomials);

      if( polynomialdata->constant == 0.0 && polynomialdata->nmonomials == 2 && polynomialdata->monomials[0]->coef == 1.0 && polynomialdata->monomials[1]->coef == 1.0 )
      {
         /* polynomial is addition of two expressions, so turn into SCIP_EXPR_PLUS */
         assert(polynomialdata->monomials[0]->nfactors == 1);
         assert(polynomialdata->monomials[0]->exponents[0] == 1.0);
         assert(polynomialdata->monomials[1]->nfactors == 1);
         assert(polynomialdata->monomials[1]->exponents[0] == 1.0);

         polynomialdataFree(blkmem, &polynomialdata);
         data->data = NULL;

         /* change operator type to PLUS */
         *op = SCIP_EXPR_PLUS;

         return SCIP_OKAY;
      }

      if( polynomialdata->constant == 0.0 && polynomialdata->nmonomials == 2 && polynomialdata->monomials[0]->coef == 1.0 && polynomialdata->monomials[1]->coef == -1.0 )
      {
         /* polynomial is substraction of two expressions, so turn into SCIP_EXPR_MINUS */
         assert(polynomialdata->monomials[0]->nfactors == 1);
         assert(polynomialdata->monomials[0]->exponents[0] == 1.0);
         assert(polynomialdata->monomials[1]->nfactors == 1);
         assert(polynomialdata->monomials[1]->exponents[0] == 1.0);

         polynomialdataFree(blkmem, &polynomialdata);
         data->data = NULL;

         /* change operator type to MINUS */
         *op = SCIP_EXPR_MINUS;

         return SCIP_OKAY;
      }

      if( polynomialdata->constant == 0.0 && polynomialdata->nmonomials == 2 && polynomialdata->monomials[0]->coef == -1.0 && polynomialdata->monomials[1]->coef == 1.0 )
      {
         /* polynomial is substraction of two expressions, so turn into SCIP_EXPR_MINUS */
         void* tmp;

         assert(polynomialdata->monomials[0]->nfactors == 1);
         assert(polynomialdata->monomials[0]->exponents[0] == 1.0);
         assert(polynomialdata->monomials[1]->nfactors == 1);
         assert(polynomialdata->monomials[1]->exponents[0] == 1.0);

         polynomialdataFree(blkmem, &polynomialdata);
         data->data = NULL;

         /* swap children */
         tmp = children[1];           /*lint !e613*/
         children[1] = children[0];   /*lint !e613*/
         children[0] = tmp;           /*lint !e613*/

         /* change operator type to MINUS */
         *op = SCIP_EXPR_MINUS;

         return SCIP_OKAY;
      }

      if( polynomialdata->constant == 0.0 )
      {
         /* check if all monomials have coefficient 1.0 */
         for( i = 0; i < polynomialdata->nmonomials; ++i )
            if( polynomialdata->monomials[i]->coef != 1.0 )
               break;

         if( i == polynomialdata->nmonomials )
         {
            /* polynomial is sum of children, so turn into SCIP_EXPR_SUM */

            polynomialdataFree(blkmem, &polynomialdata);
            data->data = NULL;

            /* change operator type to SUM */
            *op = SCIP_EXPR_SUM;

            return SCIP_OKAY;
         }
      }

      /* turn polynomial into linear expression */
      {
         SCIP_Real* lindata;

         /* monomial merging should ensure that each child appears in at most one monomial,
          * that monomials are ordered according to the child index, and that constant monomials have been removed
          */

         /* setup data of linear expression */
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &lindata, polynomialdata->nmonomials + 1) );

         for( i = 0; i < polynomialdata->nmonomials; ++i )
         {
            assert(polynomialdata->monomials[i]->childidxs[0] == i);
            assert(polynomialdata->monomials[i]->exponents[0] == 1.0);
            lindata[i] = polynomialdata->monomials[i]->coef;  /*lint !e644*/
         }
         lindata[i] = polynomialdata->constant;

         polynomialdataFree(blkmem, &polynomialdata);
         *op = SCIP_EXPR_LINEAR;
         data->data = (void*)lindata;

         return SCIP_OKAY;
      }
   }

   if( maxdegree == 2 && (polynomialdata->nmonomials > 1 || polynomialdata->constant != 0.0 || polynomialdata->monomials[0]->coef != 1.0) )
   {
      /* polynomial is quadratic expression with more than one summand or with a constant or a square or bilinear term with coefficient != 1.0, so turn into SCIP_EXPR_QUADRATIC */
      SCIP_EXPRDATA_QUADRATIC* quaddata;
      int quadelemidx;

      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &quaddata) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &quaddata->quadelems, polynomialdata->nmonomials - nlinmonomials) );
      quaddata->nquadelems = polynomialdata->nmonomials - nlinmonomials;
      quaddata->constant = polynomialdata->constant;
      quaddata->sorted = FALSE; /* quadratic data is sorted different than polynomials */

      if( nlinmonomials > 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &quaddata->lincoefs, nchildren) );
         BMSclearMemoryArray(quaddata->lincoefs, nchildren);
      }
      else
         quaddata->lincoefs = NULL;

      quadelemidx = 0;
      for( i = 0; i < polynomialdata->nmonomials; ++i )
      {
         assert(polynomialdata->monomials[i]->nfactors == 1 || polynomialdata->monomials[i]->nfactors == 2);
         if( polynomialdata->monomials[i]->nfactors == 1 )
         {
            if( polynomialdata->monomials[i]->exponents[0] == 1.0 )
            {
               /* monomial is a linear term */
               assert(quaddata->lincoefs != NULL);
               quaddata->lincoefs[polynomialdata->monomials[i]->childidxs[0]] += polynomialdata->monomials[i]->coef;
            }
            else
            {
               /* monomial should be a square term */
               assert(polynomialdata->monomials[i]->exponents[0] == 2.0);
               assert(quadelemidx < quaddata->nquadelems);
               quaddata->quadelems[quadelemidx].idx1 = polynomialdata->monomials[i]->childidxs[0];
               quaddata->quadelems[quadelemidx].idx2 = polynomialdata->monomials[i]->childidxs[0];
               quaddata->quadelems[quadelemidx].coef = polynomialdata->monomials[i]->coef;
               ++quadelemidx;
            }
         }
         else
         {
            /* monomial should be a bilinear term */
            assert(polynomialdata->monomials[i]->exponents[0] == 1.0);
            assert(polynomialdata->monomials[i]->exponents[1] == 1.0);
            assert(quadelemidx < quaddata->nquadelems);
            quaddata->quadelems[quadelemidx].idx1 = MIN(polynomialdata->monomials[i]->childidxs[0], polynomialdata->monomials[i]->childidxs[1]);
            quaddata->quadelems[quadelemidx].idx2 = MAX(polynomialdata->monomials[i]->childidxs[0], polynomialdata->monomials[i]->childidxs[1]);
            quaddata->quadelems[quadelemidx].coef = polynomialdata->monomials[i]->coef;
            ++quadelemidx;
         }
      }
      assert(quadelemidx == quaddata->nquadelems);

      polynomialdataFree(blkmem, &polynomialdata);

      *op = SCIP_EXPR_QUADRATIC;
      data->data = (void*)quaddata;

      return SCIP_OKAY;
   }

   if( polynomialdata->constant == 0.0 && polynomialdata->nmonomials == 1 && polynomialdata->monomials[0]->coef == 1.0 )
   {
      /* polynomial is product of children */
      monomial = polynomialdata->monomials[0];
      assert(monomial->nfactors == nchildren);

      if( monomial->nfactors == 1 )
      {
         /* polynomial is x^k for some k */
         assert(monomial->exponents[0] != 1.0); /* should have been handled before */
         assert(monomial->childidxs[0] == 0);

         if( monomial->exponents[0] == 2.0 )
         {
            /* polynomial is x^2, so turn into SCIP_EXPR_SQUARE */

            polynomialdataFree(blkmem, &polynomialdata);
            data->data = NULL;

            *op = SCIP_EXPR_SQUARE;

            return SCIP_OKAY;
         }

         if( EPSISINT(monomial->exponents[0], 0.0) )  /*lint !e835*/
         {
            /* k is an integer, so turn into SCIP_EXPR_INTPOWER */
            int exponent;

            exponent = (int)EPSROUND(monomial->exponents[0], 0.0);  /*lint !e835*/

            polynomialdataFree(blkmem, &polynomialdata);

            *op = SCIP_EXPR_INTPOWER;
            data->intval = exponent;

            return SCIP_OKAY;
         }

         if( monomial->exponents[0] == 0.5 )
         {
            /* polynomial is sqrt(x), so turn into SCIP_EXPR_SQRT */

            polynomialdataFree(blkmem, &polynomialdata);
            data->data = NULL;

            *op = SCIP_EXPR_SQRT;

            return SCIP_OKAY;
         }

         {
            /* polynomial is x^a with a some real number, so turn into SCIP_EXPR_REALPOWER */
            SCIP_Real exponent;

            exponent = monomial->exponents[0];

            polynomialdataFree(blkmem, &polynomialdata);

            *op = SCIP_EXPR_REALPOWER;
            data->dbl = exponent;

            return SCIP_OKAY;
         }
      }

      if( maxdegree == 2 && monomial->nfactors == 2 )
      {
         /* polynomial is product of two children, so turn into SCIP_EXPR_MUL */
         assert(monomial->exponents[0] == 1.0);
         assert(monomial->exponents[1] == 1.0);

         polynomialdataFree(blkmem, &polynomialdata);
         data->data = NULL;

         *op = SCIP_EXPR_MUL;

         return SCIP_OKAY;
      }

      if( maxdegree == monomial->nfactors )
      {
         /* polynomial is a product of n children, so turn into SCIP_EXPR_PRODUCT */

         polynomialdataFree(blkmem, &polynomialdata);
         data->data = NULL;

         *op = SCIP_EXPR_PRODUCT;

         return SCIP_OKAY;
      }

      if( monomial->nfactors == 2 && monomial->exponents[0] == 1.0 && monomial->exponents[1] == -1.0 )
      {
         /* polynomial is x/y, so turn into SCIP_EXPR_DIV */
         assert(monomial->childidxs[0] == 0);
         assert(monomial->childidxs[1] == 1);

         polynomialdataFree(blkmem, &polynomialdata);
         data->data = NULL;

         *op = SCIP_EXPR_DIV;

         return SCIP_OKAY;
      }

      if( monomial->nfactors == 2 && monomial->exponents[0] == -1.0 && monomial->exponents[1] == 1.0 )
      {
         /* polynomial is y/x, so turn into SCIP_EXPR_DIV */
         void* tmp;

         assert(monomial->childidxs[0] == 0);
         assert(monomial->childidxs[1] == 1);

         polynomialdataFree(blkmem, &polynomialdata);
         data->data = NULL;

         /* swap children */
         tmp = children[1];           /*lint !e613*/
         children[1] = children[0];   /*lint !e613*/
         children[0] = tmp;           /*lint !e613*/

         *op = SCIP_EXPR_DIV;

         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** adds copies of expressions to the array of children of a sum, product, linear, quadratic, or polynomial expression
 *
 *  For a sum or product expression, this corresponds to add additional summands and factors, resp.
 *  For a linear expression, this corresponds to add each expression with coefficient 1.0.
 *  For a quadratic or polynomial expression, only the children array may be enlarged, the expression itself remains the same.
 */
static
SCIP_RETCODE exprsimplifyAddChildren(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< quadratic or polynomial expression */
   int                   nexprs,             /**< number of expressions to add */
   SCIP_EXPR**           exprs,              /**< expressions to add */
   SCIP_Bool             comparechildren,    /**< whether to compare expressions with already existing children (no effect for sum and product) */
   SCIP_Real             eps,                /**< which epsilon to use when comparing expressions */
   int*                  childmap            /**< array where to store mapping of indices from exprs to children array in expr, or NULL if not of interest */
   )
{
   int i;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_SUM || expr->op == SCIP_EXPR_PRODUCT || expr->op == SCIP_EXPR_LINEAR || expr->op == SCIP_EXPR_QUADRATIC || expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(exprs != NULL || nexprs == 0);

   if( nexprs == 0 )
      return SCIP_OKAY;

   switch( expr->op )
   {
   case SCIP_EXPR_SUM:
   case SCIP_EXPR_PRODUCT:
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &expr->children, expr->nchildren, expr->nchildren + nexprs) );
      for( i = 0; i < nexprs; ++i )
      {
         SCIP_CALL( SCIPexprCopyDeep(blkmem, &expr->children[expr->nchildren + i], exprs[i]) );  /*lint !e613*/
         if( childmap != NULL )
            childmap[i] = expr->nchildren + i;
      }
      expr->nchildren += nexprs;

      break;
   }

   case SCIP_EXPR_LINEAR:
   case SCIP_EXPR_QUADRATIC:
   case SCIP_EXPR_POLYNOMIAL:
   {
      int j;
      int orignchildren;
      SCIP_Bool existsalready;

      orignchildren = expr->nchildren;
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &expr->children, expr->nchildren, expr->nchildren + nexprs) );

      for( i = 0; i < nexprs; ++i )
      {
         existsalready = FALSE;
         if( comparechildren )
            for( j = 0; j < orignchildren; ++j )
               /* during simplification of polynomials, their may be NULL's in children array */
               if( expr->children[j] != NULL && SCIPexprAreEqual(expr->children[j], exprs[i], eps) )  /*lint !e613*/
               {
                  existsalready = TRUE;
                  break;
               }

         if( !existsalready )
         {
            /* add copy of exprs[j] to children array */
            SCIP_CALL( SCIPexprCopyDeep(blkmem, &expr->children[expr->nchildren], exprs[i]) );  /*lint !e613*/
            if( childmap != NULL )
               childmap[i] = expr->nchildren;
            ++expr->nchildren;
         }
         else
         {
            if( childmap != NULL )
               childmap[i] = j;  /*lint !e644*/
            if( expr->op == SCIP_EXPR_LINEAR )
            {
               /* if linear expression, increase coefficient by 1.0 */
               ((SCIP_Real*)expr->data.data)[j] += 1.0;
            }
         }
      }

      /* shrink children array to actually used size */
      assert(comparechildren || expr->nchildren == orignchildren + nexprs);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &expr->children, orignchildren + nexprs, expr->nchildren) );

      if( expr->op == SCIP_EXPR_LINEAR && expr->nchildren > orignchildren )
      {
         /* if linear expression, then add 1.0 coefficients for new expressions */
         SCIP_Real* data;

         data = (SCIP_Real*)expr->data.data;
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &data, orignchildren + 1, expr->nchildren + 1) );
         data[expr->nchildren] = data[orignchildren]; /* move constant from old end to new end */
         for( i = orignchildren; i < expr->nchildren; ++i )
            data[i] = 1.0;
         expr->data.data = (void*)data;
      }
      else if( expr->op == SCIP_EXPR_QUADRATIC && expr->nchildren > orignchildren )
      {
         /* if quadratic expression, then add 0.0 linear coefficients for new expressions */
         SCIP_EXPRDATA_QUADRATIC* data;

         data = (SCIP_EXPRDATA_QUADRATIC*)expr->data.data;
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &data->lincoefs, orignchildren, expr->nchildren) );
         BMSclearMemoryArray(&data->lincoefs[orignchildren], expr->nchildren - orignchildren);  /*lint !e866*/
      }

      break;
   }

   default:
      SCIPerrorMessage("exprsimplifyAddChildren cannot be called for operand %d\n", expr->op);
      return SCIP_INVALIDDATA;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** converts expressions into polynomials, where possible and obvious */
static
SCIP_RETCODE exprsimplifyConvertToPolynomials(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR*            expr                /**< expression to convert */
   )
{
   int i;

   assert(expr != NULL);

   for( i = 0; i < expr->nchildren; ++i )
   {
      SCIP_CALL( exprsimplifyConvertToPolynomials(blkmem, expr->children[i]) );
   }

   SCIP_CALL( exprConvertToPolynomial(blkmem, &expr->op, &expr->data, expr->nchildren) );

   return SCIP_OKAY;
}

/** removes duplicate children in a polynomial expression
 *
 *  Leaves NULL's in children array.
 */
static
SCIP_RETCODE exprsimplifyRemoveDuplicatePolynomialChildren(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             eps                 /**< threshold for zero */
   )
{
   SCIP_Bool foundduplicates;
   int* childmap;
   int i;
   int j;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetOperator(expr) == SCIP_EXPR_POLYNOMIAL);

   if( expr->nchildren == 0 )
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &childmap, expr->nchildren) );

   foundduplicates = FALSE;
   for( i = 0; i < expr->nchildren; ++i )
   {
      if( expr->children[i] == NULL )
         continue;
      childmap[i] = i;  /*lint !e644*/

      for( j = i+1; j < expr->nchildren; ++j )
      {
         if( expr->children[j] == NULL )
            continue;

         if( SCIPexprAreEqual(expr->children[i], expr->children[j], eps) )
         {
            /* forget about expr j and remember that is to be replaced by i */
            SCIPexprFreeDeep(blkmem, &expr->children[j]);
            childmap[j] = i;
            foundduplicates = TRUE;
         }
      }
   }

   /* apply childmap to monomials */
   if( foundduplicates )
      polynomialdataApplyChildmap((SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data, childmap);

   /* free childmap */
   BMSfreeBlockMemoryArray(blkmem, &childmap, expr->nchildren);

   return SCIP_OKAY;
}

/** eliminates NULL's in children array and shrinks it to actual size */
static
SCIP_RETCODE exprsimplifyRemovePolynomialNullChildren(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int* childmap;
   int lastnonnull;
   int i;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetOperator(expr) == SCIP_EXPR_POLYNOMIAL);

   if( expr->nchildren == 0 )
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &childmap, expr->nchildren) );

   /* close gaps in children array */
   lastnonnull = expr->nchildren-1;
   while( lastnonnull >= 0 && expr->children[lastnonnull] == NULL )
      --lastnonnull;
   for( i = 0; i <= lastnonnull; ++i )
   {
      if( expr->children[i] != NULL )
      {
         childmap[i] = i; /* child at index i is not moved */  /*lint !e644*/
         continue;
      }
      assert(expr->children[lastnonnull] != NULL);

      /* move child at lastnonnull to position i */
      expr->children[i] = expr->children[lastnonnull];
      expr->children[lastnonnull] = NULL;
      childmap[lastnonnull] = i;

      /* update lastnonnull */
      --lastnonnull;
      while( lastnonnull >= 0 && expr->children[lastnonnull] == NULL )
         --lastnonnull;
   }
   assert(i > lastnonnull);

   /* apply childmap to monomials */
   if( lastnonnull < expr->nchildren-1 )
      polynomialdataApplyChildmap((SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data, childmap);

   BMSfreeBlockMemoryArray(blkmem, &childmap, expr->nchildren);

   /* shrink children array */
   if( lastnonnull >= 0 )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &expr->children, expr->nchildren, lastnonnull+1) );
      expr->nchildren = lastnonnull+1;
   }
   else
   {
      BMSfreeBlockMemoryArray(blkmem, &expr->children, expr->nchildren);
      expr->nchildren = 0;
   }

   return SCIP_OKAY;
}

/** checks which children are still in use and frees those which are not */
static
SCIP_RETCODE exprsimplifyRemovePolynomialUnusedChildren(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR*            expr                /**< polynomial expression */
   )
{
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
   SCIP_EXPRDATA_MONOMIAL* monomial;
   SCIP_Bool* childinuse;
   int i;
   int j;

   assert(blkmem != NULL);
   assert(expr != NULL);

   if( expr->nchildren == 0 )
      return SCIP_OKAY;

   polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data;
   assert(polynomialdata != NULL);

   /* check which children are still in use */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &childinuse, expr->nchildren) );
   BMSclearMemoryArray(childinuse, expr->nchildren);  /*lint !e644*/
   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      monomial = polynomialdata->monomials[i];
      assert(monomial != NULL);

      for( j = 0; j < monomial->nfactors; ++j )
      {
         assert(monomial->childidxs[j] >= 0);
         assert(monomial->childidxs[j] < expr->nchildren);
         childinuse[monomial->childidxs[j]] = TRUE;
      }
   }

   /* free children that are not used in any monomial */
   for( i = 0; i < expr->nchildren; ++i )
      if( expr->children[i] != NULL && !childinuse[i] )
         SCIPexprFreeDeep(blkmem, &expr->children[i]);

   BMSfreeBlockMemoryArray(blkmem, &childinuse, expr->nchildren);

   return SCIP_OKAY;
}

/** flattens polynomials in polynomials, check for constants in non-polynomials expressions
 *
 *  exprsimplifyConvertToPolynomials should have been called before to eliminate simple polynomial operands.
 */
static
SCIP_RETCODE exprsimplifyFlattenPolynomials(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             eps,                /**< threshold, under which values are treat as 0 */
   int                   maxexpansionexponent/**< maximal exponent for which we still expand non-monomial polynomials */
   )
{
   int i;

   assert(expr != NULL);

   for( i = 0; i < expr->nchildren; ++i )
   {
      SCIP_CALL( exprsimplifyFlattenPolynomials(blkmem, messagehdlr, expr->children[i], eps, maxexpansionexponent) );
   }

   switch( SCIPexprGetOperator(expr) )
   {
   case SCIP_EXPR_VARIDX:
   case SCIP_EXPR_CONST:
   case SCIP_EXPR_PARAM:
   case SCIP_EXPR_PLUS:
   case SCIP_EXPR_MINUS:
   case SCIP_EXPR_MUL:
   case SCIP_EXPR_DIV:
   case SCIP_EXPR_SQUARE:
   case SCIP_EXPR_SQRT:
   case SCIP_EXPR_INTPOWER:
   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_SIGNPOWER:
      break;

   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   {
      /* check if argument is a constant */
      if( (expr->children[0]->op == SCIP_EXPR_POLYNOMIAL && SCIPexprGetNChildren(expr->children[0]) == 0) ||
         expr->children[0]->op == SCIP_EXPR_CONST )
      {
         SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
         SCIP_Real exprval;

         /* since child0 has no children and it's polynomial was flattened, it should have no monomials */
         assert(expr->children[0]->op != SCIP_EXPR_POLYNOMIAL || SCIPexprGetNMonomials(expr->children[0]) == 0);

         /* evaluate expression in constant polynomial */
         SCIP_CALL( SCIPexprEval(expr, NULL, NULL, &exprval) );

         /* create polynomial */
         SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 0, NULL, exprval, FALSE) );

         expr->op = SCIP_EXPR_POLYNOMIAL;
         expr->data.data = (void*)polynomialdata;

         /* forget child */
         SCIPexprFreeDeep(blkmem, &expr->children[0]);
         BMSfreeBlockMemoryArray(blkmem, &expr->children, 1);
         expr->nchildren = 0;
      }

      break;
   }

   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   {
      /* check if both arguments are constants */
      if( ((expr->children[0]->op == SCIP_EXPR_POLYNOMIAL && SCIPexprGetNChildren(expr->children[0]) == 0) || expr->children[0]->op == SCIP_EXPR_CONST) &&
         ((expr->children[1]->op == SCIP_EXPR_POLYNOMIAL && SCIPexprGetNChildren(expr->children[1]) == 0) || expr->children[1]->op == SCIP_EXPR_CONST) )
      {
         SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
         SCIP_Real exprval;

         /* since children have no children and it's polynomial was flattened, it should have no monomials */
         assert(expr->children[0]->op != SCIP_EXPR_POLYNOMIAL || SCIPexprGetNMonomials(expr->children[0]) == 0);
         assert(expr->children[1]->op != SCIP_EXPR_POLYNOMIAL || SCIPexprGetNMonomials(expr->children[1]) == 0);

         /* evaluate expression in constants */
         SCIP_CALL( SCIPexprEval(expr, NULL, NULL, &exprval) );

         /* create polynomial */
         SCIP_CALL( polynomialdataCreate(blkmem, &polynomialdata, 0, NULL, exprval, FALSE) );

         expr->op = SCIP_EXPR_POLYNOMIAL;
         expr->data.data = (void*)polynomialdata;

         /* forget children */
         SCIPexprFreeDeep(blkmem, &expr->children[0]);
         SCIPexprFreeDeep(blkmem, &expr->children[1]);
         BMSfreeBlockMemoryArray(blkmem, &expr->children, 2);
         expr->nchildren = 0;
      }

      break;
   }

   case SCIP_EXPR_SUM:
   case SCIP_EXPR_PRODUCT:
   case SCIP_EXPR_LINEAR:
   case SCIP_EXPR_QUADRATIC:
   case SCIP_EXPR_USER:
      break;

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_Bool removechild;
      int* childmap;
      int childmapsize;
      int j;

      /* simplify current polynomial */
      SCIP_CALL( exprsimplifyRemoveDuplicatePolynomialChildren(blkmem, expr, eps) );
      SCIPexprMergeMonomials(blkmem, expr, eps, TRUE);

      polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data;
      assert(polynomialdata != NULL);

      SCIPdebugMessage("expand factors in expression ");
      SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
      SCIPdebugPrintf("\n");

      childmap = NULL;
      childmapsize = 0;

      /* resolve children that are constants
       * we do this first, because it reduces the degree and number of factors in the monomials,
       *   thereby allowing some expansions of polynomials that may not be possible otherwise, e.g., turning c0*c1 with c0=quadratic and c1=constant into a single monomial
       */
      for( i = 0; i < expr->nchildren; ++i )
      {
         if( expr->children[i] == NULL )
            continue;

         if( SCIPexprGetOperator(expr->children[i]) != SCIP_EXPR_CONST )
            continue;

         removechild = TRUE; /* we intend to delete children[i] */

         if( childmapsize < expr->children[i]->nchildren )
         {
            int newsize;

            newsize = calcGrowSize(expr->children[i]->nchildren);
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &childmap, childmapsize, newsize) );
            childmapsize = newsize;
         }

         /* put constant of child i into every monomial where child i is used */
         for( j = 0; j < polynomialdata->nmonomials; ++j )
         {
            int factorpos;
            SCIP_Bool success;

            monomial = polynomialdata->monomials[j];
            /* if monomial is not sorted, then polynomial should not be sorted either, or have only one monomial */
            assert(monomial->sorted || !polynomialdata->sorted || polynomialdata->nmonomials <= 1);

            if( SCIPexprFindMonomialFactor(monomial, i, &factorpos) )
            {
               assert(factorpos >= 0);
               assert(factorpos < monomial->nfactors);
               /* assert that factors have been merged */
               assert(factorpos == 0 || monomial->childidxs[factorpos-1] != i);
               assert(factorpos == monomial->nfactors-1 || monomial->childidxs[factorpos+1] != i);

               /* SCIPdebugMessage("attempt expanding child %d at monomial %d factor %d\n", i, j, factorpos);
                  SCIPdebug( SCIPexprPrint(expr, NULL, NULL, NULL) ); SCIPdebugPrintf("\n");
                  SCIPdebug( SCIPexprPrint(expr->children[i], NULL, NULL, NULL) ); SCIPdebugPrintf("\n"); */

               if( !EPSISINT(monomial->exponents[factorpos], 0.0) && SCIPexprGetOpReal(expr->children[i]) < 0.0 )  /*lint !e835*/
               {
                  /* if constant is negative and our exponent is not integer, then cannot do expansion */
                  SCIPmessagePrintWarning(messagehdlr, "got negative constant %g to the power of a noninteger exponent %g\n",
                     SCIPexprGetOpReal(expr->children[i]), monomial->exponents[factorpos]);
                  success = FALSE;
               }
               else
               {
                  monomial->coef *= pow(SCIPexprGetOpReal(expr->children[i]), monomial->exponents[factorpos]);

                  /* move last factor to position factorpos */
                  if( factorpos < monomial->nfactors-1 )
                  {
                     monomial->exponents[factorpos] = monomial->exponents[monomial->nfactors-1];
                     monomial->childidxs[factorpos] = monomial->childidxs[monomial->nfactors-1];
                  }
                  --monomial->nfactors;
                  monomial->sorted = FALSE;
                  polynomialdata->sorted = FALSE;

                  success = TRUE;
               }

               if( !success )
                  removechild = FALSE;
            }
         }

         /* forget about child i, if it is not used anymore */
         if( removechild )
            SCIPexprFreeDeep(blkmem, &expr->children[i]);

         /* simplify current polynomial again */
         SCIPexprMergeMonomials(blkmem, expr, eps, TRUE);
      }

      /* try to resolve children that are polynomials itself */
      for( i = 0; i < expr->nchildren; ++i )
      {
         if( expr->children[i] == NULL )
            continue;

         if( SCIPexprGetOperator(expr->children[i]) != SCIP_EXPR_POLYNOMIAL )
            continue;

         removechild = TRUE; /* we intend to delete children[i] */

         if( childmapsize < expr->children[i]->nchildren )
         {
            int newsize;

            newsize = calcGrowSize(expr->children[i]->nchildren);
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &childmap, childmapsize, newsize) );
            childmapsize = newsize;
         }

         /* add children of child i */
         SCIP_CALL( exprsimplifyAddChildren(blkmem, expr, expr->children[i]->nchildren, expr->children[i]->children, TRUE, eps, childmap) );

         /* put polynomial of child i into every monomial where child i is used */
         j = 0;
         while( j < polynomialdata->nmonomials )
         {
            int factorpos;
            SCIP_Bool success;

            monomial = polynomialdata->monomials[j];
            /* if monomial is not sorted, then polynomial should not be sorted either, or have only one monomial */
            assert(monomial->sorted || !polynomialdata->sorted || polynomialdata->nmonomials <= 1);

            if( SCIPexprFindMonomialFactor(monomial, i, &factorpos) )
            {
               assert(factorpos >= 0);
               assert(factorpos < monomial->nfactors);
               /* assert that factors have been merged */
               assert(factorpos == 0 || monomial->childidxs[factorpos-1] != i);
               assert(factorpos == monomial->nfactors-1 || monomial->childidxs[factorpos+1] != i);

               /* SCIPdebugMessage("attempt expanding child %d at monomial %d factor %d\n", i, j, factorpos);
                  SCIPdebug( SCIPexprPrint(expr, NULL, NULL, NULL) ); SCIPdebugPrintf("\n");
                  SCIPdebug( SCIPexprPrint(expr->children[i], NULL, NULL, NULL) ); SCIPdebugPrintf("\n"); */

               SCIP_CALL( polynomialdataExpandMonomialFactor(blkmem, messagehdlr, polynomialdata, j, factorpos,
                     (SCIP_EXPRDATA_POLYNOMIAL*)expr->children[i]->data.data, childmap, maxexpansionexponent, &success) );

               if( !success )
               {
                  removechild = FALSE;
                  ++j;
               }
            }
            else
               ++j;

            /* expansion may remove monomials[j], move a monomial from the end to position j, or add new monomials to the end of polynomialdata
             * we thus repeat with index j, if a factor was successfully expanded
             */
         }

         /* forget about child i, if it is not used anymore */
         if( removechild )
            SCIPexprFreeDeep(blkmem, &expr->children[i]);

         /* simplify current polynomial again */
         SCIPexprMergeMonomials(blkmem, expr, eps, TRUE);
      }

      BMSfreeBlockMemoryArrayNull(blkmem, &childmap, childmapsize);

      /* free children that are not in use anymore */
      SCIP_CALL( exprsimplifyRemovePolynomialUnusedChildren(blkmem, expr) );

      /* remove NULLs from children array */
      SCIP_CALL( exprsimplifyRemovePolynomialNullChildren(blkmem, expr) );

      /* if no children left, then it's a constant polynomial -> change into EXPR_CONST */
      if( expr->nchildren == 0 )
      {
         SCIP_Real val;

         /* if no children, then it should also have no monomials */
         assert(polynomialdata->nmonomials == 0);

         val = polynomialdata->constant;
         polynomialdataFree(blkmem, &polynomialdata);

         expr->op = SCIP_EXPR_CONST;
         expr->data.dbl = val;
      }

      SCIPdebugMessage("-> ");
      SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
      SCIPdebugPrintf("\n");

      break;
   }

   case SCIP_EXPR_LAST:
      break;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** separates linear monomials from an expression, if it is a polynomial expression
 *
 *  Separates only those linear terms whose variable is not used otherwise in the expression.
 */
static
SCIP_RETCODE exprsimplifySeparateLinearFromPolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             eps,                /**< threshold, under which positive values are treat as 0 */
   int                   nvars,              /**< number of variables in expression */
   int*                  nlinvars,           /**< buffer to store number of linear variables in linear part */
   int*                  linidxs,            /**< array to store indices of variables in expression tree which belong to linear part */
   SCIP_Real*            lincoefs            /**< array to store coefficients of linear part */
   )
{
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
   SCIP_EXPRDATA_MONOMIAL* monomial;
   int* varsusage;
   int* childusage;
   int childidx;
   int i;
   int j;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(nlinvars != NULL);
   assert(linidxs != NULL);
   assert(lincoefs != NULL);

   *nlinvars = 0;

   if( SCIPexprGetOperator(expr) != SCIP_EXPR_POLYNOMIAL )
      return SCIP_OKAY;

   if( SCIPexprGetNChildren(expr) == 0 )
      return SCIP_OKAY;

   polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data;
   assert(polynomialdata != NULL);

   /* get variable usage */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &varsusage, nvars) );
   BMSclearMemoryArray(varsusage, nvars);  /*lint !e644*/
   SCIPexprGetVarsUsage(expr, varsusage);

   /* get child usage: how often each child is used in the polynomial */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &childusage, expr->nchildren) );
   BMSclearMemoryArray(childusage, expr->nchildren);  /*lint !e644*/
   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      monomial = polynomialdata->monomials[i];
      assert(monomial != NULL);
      for( j = 0; j < monomial->nfactors; ++j )
      {
         assert(monomial->childidxs[j] >= 0);
         assert(monomial->childidxs[j] < expr->nchildren);
         ++childusage[monomial->childidxs[j]];
      }
   }

   /* move linear monomials out of polynomial */
   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      monomial = polynomialdata->monomials[i];
      assert(monomial != NULL);
      if( monomial->nfactors != 1 )
         continue;
      if( monomial->exponents[0] != 1.0 )
         continue;
      childidx = monomial->childidxs[0];
      if( SCIPexprGetOperator(expr->children[childidx]) != SCIP_EXPR_VARIDX )
         continue;

      /* we are at a linear monomial in a variable */
      assert(SCIPexprGetOpIndex(expr->children[childidx]) < nvars);
      if( childusage[childidx] == 1 && varsusage[SCIPexprGetOpIndex(expr->children[childidx])] == 1 )
      {
         /* if the child expression is not used in another monomial (which would due to merging be not linear)
          * and if the variable is not used somewhere else in the tree,
          * then move this monomial into linear part and free child
          */
         linidxs[*nlinvars]  = SCIPexprGetOpIndex(expr->children[childidx]);
         lincoefs[*nlinvars] = monomial->coef;
         ++*nlinvars;

         SCIPexprFreeDeep(blkmem, &expr->children[childidx]);
         monomial->coef = 0.0;
         monomial->nfactors = 0;
      }
   }

   BMSfreeBlockMemoryArray(blkmem, &varsusage, nvars);
   BMSfreeBlockMemoryArray(blkmem, &childusage, expr->nchildren);

   if( *nlinvars > 0 )
   {
      /* if we did something, cleanup polynomial (e.g., remove monomials with coefficient 0.0) */
      polynomialdataMergeMonomials(blkmem, polynomialdata, eps, FALSE);
      SCIP_CALL( exprsimplifyRemovePolynomialNullChildren(blkmem, expr) );
   }

   return SCIP_OKAY;
}

/** converts polynomial expressions back into simpler expressions, where possible */
static
SCIP_RETCODE exprsimplifyUnconvertPolynomials(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR*            expr                /**< expression to convert back */
   )
{
   int i;

   assert(blkmem != NULL);
   assert(expr != NULL);

   for( i = 0; i < expr->nchildren; ++i )
   {
      SCIP_CALL( exprsimplifyUnconvertPolynomials(blkmem, expr->children[i]) );
   }

   if( expr->op != SCIP_EXPR_POLYNOMIAL )
      return SCIP_OKAY;

   SCIP_CALL( exprUnconvertPolynomial(blkmem, &expr->op, &expr->data, expr->nchildren, (void**)expr->children) );

   return SCIP_OKAY;
}

static
SCIP_DECL_HASHGETKEY( exprparseVarTableGetKey )
{  /*lint --e{715}*/
   return (void*)((char*)elem + sizeof(int));
}

/** parses a variable name from a string and creates corresponding expression
 *
 *  Creates a new variable index if variable not seen before, updates varnames and vartable structures.
 */
static
SCIP_RETCODE exprparseReadVariable(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   const char**          str,                /**< pointer to the string to be parsed */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to created expression */
   int*                  nvars,              /**< running number of encountered variables so far */
   int**                 varnames,           /**< pointer to buffer to store new variable names */
   int*                  varnameslength,     /**< pointer to length of the varnames buffer array */
   SCIP_HASHTABLE*       vartable,           /**< hash table for variable names and corresponding expression index */
   SCIP_Real             coefficient,        /**< coefficient to be used when creating the expression */
   const char*           varnameendptr       /**< if a \<varname\> should be parsed, set this to NULL. Then, str points to the '<'
                                                  else, str should point to the first letter of the varname, and varnameendptr should
                                                  point one char behind the last char of the variable name */
   )
{
   int namelength;
   int varidx;
   char varname[SCIP_MAXSTRLEN];
   void* element;

   assert(blkmem != NULL);
   assert(str != NULL);
   assert(expr != NULL);
   assert(nvars != NULL);
   assert(varnames != NULL);
   assert(vartable != NULL);

   if( varnameendptr == NULL )
   {
      ++*str;
      varnameendptr = *str;
      while( varnameendptr[0] != '>' )
         ++varnameendptr;
   }

   namelength = varnameendptr - *str; /*lint !e712*/
   if( namelength >= SCIP_MAXSTRLEN )
   {
      SCIPerrorMessage("Variable name %.*s is too long for buffer in exprparseReadVariable.\n", namelength, str);
      return SCIP_READERROR;
   }

   memcpy(varname, *str, namelength * sizeof(char));
   varname[namelength] = '\0';

   element = SCIPhashtableRetrieve(vartable, varname);
   if( element != NULL )
   {
      /* variable is old friend */
      assert(strcmp((char*)element + sizeof(int), varname) == 0);

      varidx = *(int*)element;
   }
   else
   {
      /* variable is new */
      varidx = *nvars;

      (*varnameslength) -= (int)(1 + (strlen(varname) + 1) / sizeof(int) + 1);
      if( *varnameslength < 0 )
      {
         SCIPerrorMessage("Buffer in exprparseReadVariable is too short for varaible name %.*s.\n", namelength, str);
         return SCIP_READERROR;
      }

      /* store index of variable and variable name in varnames buffer */
      **varnames = varidx;
      strcpy((char*)(*varnames + 1), varname);

      /* insert variable into hashtable */
      SCIP_CALL( SCIPhashtableInsert(vartable, (void*)*varnames) );

      ++*nvars;
      *varnames += 1 + (strlen(varname) + 1) / sizeof(int) + 1;
   }

   /* create VARIDX expression, put into LINEAR expression if we have coefficient != 1 */
   SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_VARIDX, varidx) );  /*lint !e613*/
   if( coefficient != 1.0 )
   {
      SCIP_CALL( SCIPexprCreateLinear(blkmem, expr, 1, expr, &coefficient, 0.0) );
   }

   /* Move pointer to char behind end of variable */
   *str = varnameendptr + 1;

   /* consprint sometimes prints a variable type identifier which we don't need */
   if( (*str)[0] == '[' && (*str)[2] == ']' &&
       ((*str)[1] == SCIP_VARTYPE_BINARY_CHAR  ||
        (*str)[1] == SCIP_VARTYPE_INTEGER_CHAR ||
        (*str)[1] == SCIP_VARTYPE_IMPLINT_CHAR ||
        (*str)[1] == SCIP_VARTYPE_CONTINUOUS_CHAR ) )
      *str += 3;

   return SCIP_OKAY;
}

/** if str[0] points to an opening parenthesis, this function sets endptr to point to the matching closing bracket in str
 *
 *  Searches for at most length characters.
 */
static
SCIP_RETCODE exprparseFindClosingParenthesis(
   const char*           str,                /**< pointer to the string to be parsed */
   const char**          endptr,             /**< pointer to point to the closing parenthesis */
   int                   length              /**< length of the string to be parsed */
   )
{
   int nopenbrackets;

   assert(str[0] == '(');

   *endptr = str;

   /* find the end of this expression */
   nopenbrackets = 0;
   while( (*endptr - str ) < length && !(nopenbrackets == 1 && *endptr[0] == ')') )
   {
      if( *endptr[0] == '(')
         ++nopenbrackets;
      if( *endptr[0] == ')')
         --nopenbrackets;
      ++*endptr;
   }

   if( *endptr[0] != ')' )
   {
      SCIPerrorMessage("unable to find closing parenthesis in unbalanced expression %.*s\n", length, str);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** this function sets endptr to point to the next separating comma in str
 *
 *  That is, for a given string like "x+f(x,y),z", endptr will point to the comma before "z"
 *
 *  Searches for at most length characters.
 */
static
SCIP_RETCODE exprparseFindSeparatingComma(
   const char*           str,                /**< pointer to the string to be parsed */
   const char**          endptr,             /**< pointer to point to the comma */
   int                   length              /**< length of the string to be parsed */
   )
{
   int nopenbrackets;

   *endptr = str;

   /* find a comma without open brackets */
   nopenbrackets = 0;
   while( (*endptr - str ) < length && !(nopenbrackets == 0 && *endptr[0] == ',') )
   {
      if( *endptr[0] == '(')
         ++nopenbrackets;
      if( *endptr[0] == ')')
         --nopenbrackets;
      ++*endptr;
   }

   if( *endptr[0] != ',' )
   {
      SCIPerrorMessage("unable to find separating comma in unbalanced expression %.*s\n", length, str);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** parses an expression from a string */
static
SCIP_RETCODE exprParse(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to created expression */
   const char*           str,                /**< pointer to the string to be parsed */
   int                   length,             /**< length of the string to be parsed */
   const char*           lastchar,           /**< pointer to the last char of str that should be parsed */
   int*                  nvars,              /**< running number of encountered variables so far */
   int**                 varnames,           /**< pointer to buffer to store new variable names */
   int*                  varnameslength,     /**< pointer to length of the varnames buffer array */
   SCIP_HASHTABLE*       vartable,           /**< hash table for variable names and corresponding expression index */
   int                   recursiondepth      /**< current recursion depth */
   )
{   /*lint --e{712,747}*/
   SCIP_EXPR* arg1;
   SCIP_EXPR* arg2;
   const char* subexpptr;
   const char* subexpendptr;
   const char* strstart;
   const char* endptr;
   char* nonconstendptr;
   SCIP_Real number;
   int subexplength;
   int nopenbrackets;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(str != NULL);
   assert(lastchar >= str);
   assert(nvars != NULL);
   assert(varnames != NULL);
   assert(vartable != NULL);

   assert(recursiondepth < 100);

   strstart = str; /* might be needed for error message... */

   SCIPdebugMessage("exprParse (%i): parsing %.*s\n", recursiondepth, (int) (lastchar-str + 1), str);

   /* ignore whitespace */
   while( isspace((unsigned char)*str) )
      ++str;

   /* look for a sum or difference not contained in brackets */
   subexpptr = str;
   nopenbrackets = 0;

   /* find the end of this expression
    * a '+' right at the beginning indicates a coefficient, not treated here, or a summation
    */
   while( subexpptr != lastchar && !(nopenbrackets == 0 && (subexpptr[0] == '+' || subexpptr[0] == '-') && subexpptr != str) )
   {
      if( subexpptr[0] == '(')
         ++nopenbrackets;
      if( subexpptr[0] == ')')
         --nopenbrackets;
      ++subexpptr;
   }

   if( subexpptr != lastchar )
   {
      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg1, str, (int) ((subexpptr - 1) - str + 1), subexpptr - 1, nvars,
            varnames, varnameslength, vartable, recursiondepth + 1) );

      if( subexpptr[0] == '+' )
         ++subexpptr;
      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg2, subexpptr , (int) (lastchar - (subexpptr ) + 1), lastchar, nvars,
            varnames, varnameslength, vartable, recursiondepth + 1) );

      /* make new expression from two arguments
       * we always use add, because we leave the operator between the found expressions in the second argument
       * this way, we do not have to worry about ''minus brackets'' in the case of more then two summands:
       *   a - b - c = a + (-b -c)
       */
      SCIP_CALL( SCIPexprAdd(blkmem, expr, 1.0, arg1, 1.0, arg2, 0.0) );

      SCIPdebugMessage("exprParse (%i): returns expression ", recursiondepth);
      SCIPdebug( SCIPexprPrint(*expr, messagehdlr, NULL, NULL, NULL, NULL) );
      SCIPdebug( SCIPmessagePrintInfo(messagehdlr, "\n") );

      return SCIP_OKAY;
   }

   /* check for a bracketed subexpression */
   if( str[0] == '(' )
   {
      nopenbrackets = 0;

      subexplength = -1;    /* we do not want the closing bracket in the string */
      subexpptr = str + 1;  /* leave out opening bracket */

      /* find the end of this expression */
      while( subexplength < length && !(nopenbrackets == 1 && str[0] == ')') )
      {
         if( str[0] == '(' )
            ++nopenbrackets;
         if( str[0] == ')' )
            --nopenbrackets;
         ++str;
         ++subexplength;
      }
      subexpendptr = str - 1; /* leave out closing bracket */

      SCIP_CALL( exprParse(blkmem, messagehdlr, expr, subexpptr, subexplength, subexpendptr, nvars, varnames,
            varnameslength, vartable, recursiondepth + 1) );
      ++str;
   }
   else if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+')
         && (isdigit((unsigned char)str[1]) || str[1] == ' ')) )
   {
      /* check if there is a lonely minus coming, indicating a -1.0 */
      if( str[0] == '-'  && str[1] == ' ' )
      {
         number = -1.0;
         nonconstendptr = (char*) str + 1;
      }
      /* check if there is a number coming */
      else if( !SCIPstrToRealValue(str, &number, &nonconstendptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_READERROR;
      }
      str = nonconstendptr;

      /* ignore whitespace */
      while( isspace((unsigned char)*str) && str != lastchar )
         ++str;

      if( str[0] != '*' && str[0] != '/' && str[0] != '+' && str[0] != '-' && str[0] != '^' )
      {
         if( str < lastchar )
         {
            SCIP_CALL( exprParse(blkmem, messagehdlr, expr, str, (int)(lastchar - str) + 1, lastchar, nvars, varnames,
                  varnameslength, vartable, recursiondepth + 1) );
            SCIP_CALL( SCIPexprMulConstant(blkmem, expr, *expr, number) );
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_CONST, number) );
         }
         str = lastchar + 1;
      }
      else
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_CONST, number) );
      }
   }
   else if( str[0] == '<' )
   {
      /* check if expressions begins with a variable */
      SCIP_CALL( exprparseReadVariable(blkmem, &str, expr, nvars, varnames, varnameslength, vartable, 1.0, NULL) );
   }
   /* four character operators */
   else if( strncmp(str, "sqrt", 4) == 0 )
   {
      str += 4;
      SCIP_CALL( exprparseFindClosingParenthesis(str, &endptr, length) );
      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg1, str + 1, endptr - str - 1, endptr -1, nvars, varnames,
            varnameslength, vartable, recursiondepth + 1) );
      str = endptr + 1;

      SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_SQRT, arg1) );
   }
   /* three character operators with 1 argument */
   else if(
      strncmp(str, "abs", 3) == 0 ||
      strncmp(str, "cos", 3) == 0 ||
      strncmp(str, "exp", 3) == 0 ||
      strncmp(str, "log", 3) == 0 ||
      strncmp(str, "sin", 3) == 0 ||
      strncmp(str, "sqr", 3) == 0 ||
      strncmp(str, "tan", 3) == 0 )
   {
      const char* opname = str;

      str += 3;
      SCIP_CALL( exprparseFindClosingParenthesis(str, &endptr, length) );
      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg1, str + 1, endptr - str - 1, endptr -1, nvars, varnames,
            varnameslength, vartable, recursiondepth + 1) );
      str = endptr + 1;

      if( strncmp(opname, "abs", 3) == 0 )
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_ABS, arg1) );
      }
      else if( strncmp(opname, "cos", 3) == 0 )
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_COS, arg1) );
      }
      else if( strncmp(opname, "exp", 3) == 0 )
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_EXP, arg1) );
      }
      else if( strncmp(opname, "log", 3) == 0 )
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_LOG, arg1) );
      }
      else if( strncmp(opname, "sin", 3) == 0 )
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_SIN, arg1) );
      }
      else if( strncmp(opname, "sqr", 3) == 0 )
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_SQUARE, arg1) );
      }
      else
      {
         assert(strncmp(opname, "tan", 3) == 0);
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_TAN, arg1) );
      }
   }
   /* three character operators with 2 arguments */
   else if(
      strncmp(str, "max", 3) == 0 ||
      strncmp(str, "min", 3) == 0 )
   {
      /* we have a string of the form "min(...,...)" or "max(...,...)"
       * first find the closing parenthesis, then the comma
       */
      const char* comma;
      SCIP_EXPROP op;

      op = (str[1] == 'a' ? SCIP_EXPR_MAX : SCIP_EXPR_MIN);

      str += 3;
      SCIP_CALL( exprparseFindClosingParenthesis(str, &endptr, length) );

      SCIP_CALL( exprparseFindSeparatingComma(str+1, &comma, endptr - str - 1) );

      /* parse first argument [str+1..comma-1] */
      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg1, str + 1, comma - str - 1, comma - 1, nvars, varnames,
            varnameslength, vartable, recursiondepth + 1) );

      /* parse second argument [comma+1..endptr] */
      ++comma;
      while( comma < endptr && *comma == ' ' )
         ++comma;

      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg2, comma, endptr - comma, endptr - 1, nvars, varnames,
            varnameslength, vartable, recursiondepth + 1) );

      SCIP_CALL( SCIPexprCreate(blkmem, expr, op, arg1, arg2) );

      str = endptr + 1;
   }
   else if( strncmp(str, "power", 5) == 0 )
   {
      /* we have a string of the form "power(...,integer)" (thus, intpower)
       * first find the closing parenthesis, then the comma
       */
      const char* comma;
      int exponent;

      str += 5;
      SCIP_CALL( exprparseFindClosingParenthesis(str, &endptr, length) );

      SCIP_CALL( exprparseFindSeparatingComma(str+1, &comma, endptr - str - 1) );

      /* parse first argument [str+1..comma-1] */
      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg1, str + 1, comma - str - 1, comma - 1, nvars, varnames,
            varnameslength, vartable, recursiondepth + 1) );

      ++comma;
      /* parse second argument [comma, endptr-1]: it needs to be an integer */
      while( comma < endptr && *comma == ' ' )
         ++comma;
      if( !isdigit((unsigned char)comma[0]) && !((comma[0] == '-' || comma[0] == '+') && isdigit((unsigned char)comma[1])) )
      {
         SCIPerrorMessage("error parsing integer exponent from <%s>\n", comma);
      }
      if( !SCIPstrToIntValue(comma, &exponent, &nonconstendptr) )
      {
         SCIPerrorMessage("error parsing integer from <%s>\n", comma);
         return SCIP_READERROR;
      }

      SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_INTPOWER, arg1, exponent) );

      str = endptr + 1;
   }
   else if( strncmp(str, "realpower", 9) == 0 || strncmp(str, "signpower", 9) == 0 )
   {
      /* we have a string of the form "realpower(...,double)" or "signpower(...,double)"
       * first find the closing parenthesis, then the comma
       */
      const char* opname = str;
      const char* comma;

      str += 9;
      SCIP_CALL( exprparseFindClosingParenthesis(str, &endptr, length) );

      SCIP_CALL( exprparseFindSeparatingComma(str+1, &comma, endptr - str - 1) );

      /* parse first argument [str+1..comma-1] */
      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg1, str + 1, comma - str - 1, comma - 1, nvars, varnames,
            varnameslength, vartable, recursiondepth + 1) );

      ++comma;
      /* parse second argument [comma, endptr-1]: it needs to be an number */
      while( comma < endptr && *comma == ' ' )
         ++comma;
      if( !isdigit((unsigned char)comma[0]) && !((comma[0] == '-' || comma[0] == '+') && isdigit((unsigned char)comma[1])) )
      {
         SCIPerrorMessage("error parsing number exponent from <%s>\n", comma);
      }
      if( !SCIPstrToRealValue(comma, &number, &nonconstendptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", comma);
         return SCIP_READERROR;
      }

      if( strncmp(opname, "realpower", 9) == 0 )
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_REALPOWER, arg1, number) );
      }
      else
      {
         assert(strncmp(opname, "signpower", 9) == 0);
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_SIGNPOWER, arg1, number) );
      }

      str = endptr + 1;
   }
   else if( isalpha(*str) || *str == '_' || *str == '#' )
   {
      /* check for a variable, that was not recognized earlier because somebody omitted the '<' and '>' we need for
       * SCIPparseVarName, making everyones life harder;
       * we allow only variable names starting with a character or underscore here
       */
      const char* varnamestartptr = str;

      /* allow only variable names containing characters, digits, hash marks, and underscores here */
      while( isalnum(str[0]) || str[0] == '_' || str[0] == '#' )
         ++str;

      SCIP_CALL( exprparseReadVariable(blkmem, &varnamestartptr, expr, nvars, varnames, varnameslength,
            vartable, 1.0, str) );
   }
   else
   {
      SCIPerrorMessage("parsing of invalid expression %.*s.\n", (int) (lastchar - str + 1), str);
      return SCIP_READERROR;
   }

   /* if we are one char behind lastchar, we are done */
   if( str == lastchar + 1)
   {
      SCIPdebugMessage("exprParse (%i): returns expression ", recursiondepth);
      SCIPdebug( SCIPexprPrint(*expr, messagehdlr, NULL, NULL, NULL, NULL) );
      SCIPdebug( SCIPmessagePrintInfo(messagehdlr, "\n") );

      return SCIP_OKAY;
   }

   /* check if we are still in bounds */
   if( str > lastchar + 1)
   {
      SCIPerrorMessage("error finding first expression in \"%.*s\" took us outside of given subexpression length\n", length, strstart);
      return SCIP_READERROR;
   }

   /* ignore whitespace */
   while( isspace((unsigned char)*str) && str != lastchar + 1 )
      ++str;

   /* maybe now we're done? */
   if( str >= lastchar + 1)
   {
      SCIPdebugMessage("exprParse (%i): returns expression ", recursiondepth);
      SCIPdebug( SCIPexprPrint(*expr, messagehdlr, NULL, NULL, NULL, NULL) );
      SCIPdebug( SCIPmessagePrintInfo(messagehdlr, "\n") );

      return SCIP_OKAY;
   }

   if( str[0] == '^' )
   {
      /* a '^' behind the found expression indicates a power */
      SCIP_Real constant;

      arg1 = *expr;
      ++str;
      while( isspace((unsigned char)*str) && str != lastchar + 1 )
         ++str;

      if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
      {
         /* there is a number coming */
         if( !SCIPstrToRealValue(str, &number, &nonconstendptr) )
         {
            SCIPerrorMessage("error parsing number from <%s>\n", str);
            return SCIP_READERROR;
         }

         SCIP_CALL( SCIPexprCreate(blkmem, &arg2, SCIP_EXPR_CONST, number) );
         str = nonconstendptr;

         constant = SCIPexprGetOpReal(arg2);
         SCIPexprFreeDeep(blkmem, &arg2);

         /* expr^number is intpower or realpower */
         if( EPSISINT(constant, 0.0) ) /*lint !e835*/
         {
            SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_INTPOWER, arg1, (int)constant) );
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_REALPOWER, arg1, constant) );
         }
      }
      else if( str[0] == '(' )
      {
         /* we use exprParse to evaluate the exponent */

         SCIP_CALL( exprparseFindClosingParenthesis(str, &endptr, length) );
         SCIP_CALL( exprParse(blkmem, messagehdlr, &arg2, str + 1, endptr - str - 1, endptr -1, nvars, varnames,
               varnameslength, vartable, recursiondepth + 1) );

         if( SCIPexprGetOperator(arg2) != SCIP_EXPR_CONST )
         {
            /* reformulate arg1^arg2 as exp(arg2 * log(arg1)) */
            SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_LOG, arg1) );
            SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_MUL, *expr, arg2) );
            SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_EXP, *expr) );
         }
         else
         {
            /* expr^number is intpower or realpower */
            constant = SCIPexprGetOpReal(arg2);
            SCIPexprFreeDeep(blkmem, &arg2);
            if( EPSISINT(constant, 0.0) ) /*lint !e835*/
            {
               SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_INTPOWER, arg1, (int)constant) );
            }
            else
            {
               SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_REALPOWER, arg1, constant) );
            }
         }
         str = endptr + 1;
      }
      else
      {
         SCIPerrorMessage("unexpected string following ^ in  %.*s\n", length, str);
         return SCIP_READERROR;
      }

      /* ignore whitespace */
      while( isspace((unsigned char)*str) && str != lastchar + 1 )
         ++str;
   }

   /* check for a two argument operator that is not a multiplication */
   if( str <= lastchar && (str[0] == '+' || str[0] == '-' || str[0] == '/') )
   {
      char op;

      op = str[0];
      arg1 = *expr;

      /* step forward over the operator to go to the beginning of the second argument */
      ++str;

      SCIP_CALL( exprParse(blkmem, messagehdlr, &arg2, str, (int) (lastchar - str + 1), lastchar, nvars, varnames,
            varnameslength, vartable, recursiondepth + 1) );
      str = lastchar + 1;

      /* make new expression from two arguments */
      if( op == '+')
      {
         SCIP_CALL( SCIPexprAdd(blkmem, expr, 1.0, arg1, 1.0, arg2, 0.0) );
      }
      else if( op == '-')
      {
         SCIP_CALL( SCIPexprAdd(blkmem, expr, 1.0, arg1, -1.0, arg2, 0.0) );
      }
      else if( op == '*' )
      {
         if( SCIPexprGetOperator(arg1) == SCIP_EXPR_CONST )
         {
            SCIP_CALL( SCIPexprMulConstant(blkmem, expr, arg2, SCIPexprGetOpReal(arg1)) );
         }
         else if( SCIPexprGetOperator(arg2) == SCIP_EXPR_CONST )
         {
            SCIP_CALL( SCIPexprMulConstant(blkmem, expr, arg1, SCIPexprGetOpReal(arg2)) );
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_MUL, arg1, arg2) );
         }
      }
      else
      {
         assert(op == '/');

         if( SCIPexprGetOperator(arg2) == SCIP_EXPR_CONST )
         {
            SCIP_CALL( SCIPexprMulConstant(blkmem, expr, arg1, 1.0 / SCIPexprGetOpReal(arg2)) );
            SCIPexprFreeShallow(blkmem, &arg2);
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_DIV, arg1, arg2) );
         }
      }
   }

   /* ignore whitespace */
   while( isspace((unsigned char)*str) )
      ++str;

   /* we are either done or we have a multiplication? */
   if( str >= lastchar + 1 )
   {
      SCIPdebugMessage("exprParse (%i): returns expression ", recursiondepth);
      SCIPdebug( SCIPexprPrint(*expr, messagehdlr, NULL, NULL, NULL, NULL) );
      SCIPdebug( SCIPmessagePrintInfo(messagehdlr, "\n") );

      return SCIP_OKAY;
   }

   /* if there is a part of the string left to be parsed, we assume that this as a multiplication */
   arg1 = *expr;

   /* stepping over multiplication operator if needed */
   if( str[0] == '*' )
   {
      ++str;
   }
   else if( str[0] != '(' )
   {
      SCIPdebugMessage("No operator found, assuming a multiplication before %.*s\n", (int) (lastchar - str + 1), str);
   }

   SCIP_CALL( exprParse(blkmem, messagehdlr, &arg2, str, (int) (lastchar - str + 1), lastchar, nvars, varnames,
         varnameslength, vartable, recursiondepth + 1) );

   if( SCIPexprGetOperator(arg1) == SCIP_EXPR_CONST )
   {
      SCIP_CALL( SCIPexprMulConstant(blkmem, expr, arg2, SCIPexprGetOpReal(arg1)) );
      SCIPexprFreeDeep(blkmem, &arg1);
   }
   else if( SCIPexprGetOperator(arg2) == SCIP_EXPR_CONST )
   {
      SCIP_CALL( SCIPexprMulConstant(blkmem, expr, arg1, SCIPexprGetOpReal(arg2)) );
      SCIPexprFreeDeep(blkmem, &arg2);
   }
   else
   {
      SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_MUL, arg1, arg2) );
   }

   SCIPdebugMessage("exprParse (%i): returns expression ", recursiondepth);
   SCIPdebug( SCIPexprPrint(*expr, messagehdlr, NULL, NULL, NULL, NULL) );
   SCIPdebug( SCIPmessagePrintInfo(messagehdlr, "\n") );

   return SCIP_OKAY;
}

/**@} */

/**@name Expression methods */
/**@{ */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPexprGetOperator
#undef SCIPexprGetNChildren
#undef SCIPexprGetChildren
#undef SCIPexprGetOpIndex
#undef SCIPexprGetOpReal
#undef SCIPexprGetOpData
#undef SCIPexprGetRealPowerExponent
#undef SCIPexprGetIntPowerExponent
#undef SCIPexprGetSignPowerExponent
#undef SCIPexprGetLinearCoefs
#undef SCIPexprGetLinearConstant
#undef SCIPexprGetQuadElements
#undef SCIPexprGetQuadConstant
#undef SCIPexprGetQuadLinearCoefs
#undef SCIPexprGetNQuadElements
#undef SCIPexprGetMonomials
#undef SCIPexprGetNMonomials
#undef SCIPexprGetPolynomialConstant
#undef SCIPexprGetMonomialCoef
#undef SCIPexprGetMonomialNFactors
#undef SCIPexprGetMonomialChildIndices
#undef SCIPexprGetMonomialExponents
#undef SCIPexprGetUserData
#undef SCIPexprHasUserEstimator
#undef SCIPexprGetUserEvalCapability

/** gives operator of expression */
SCIP_EXPROP SCIPexprGetOperator(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->op;
}

/** gives number of children of an expression */
int SCIPexprGetNChildren(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->nchildren;
}

/** gives pointer to array with children of an expression */
SCIP_EXPR** SCIPexprGetChildren(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->children;
}

/** gives index belonging to a SCIP_EXPR_VARIDX or SCIP_EXPR_PARAM operand */
int SCIPexprGetOpIndex(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_VARIDX || expr->op == SCIP_EXPR_PARAM);

   return expr->data.intval;
}

/** gives real belonging to a SCIP_EXPR_CONST operand */
SCIP_Real SCIPexprGetOpReal(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_CONST);

   return expr->data.dbl;
}

/** gives void* belonging to a complex operand */
void* SCIPexprGetOpData(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op >= SCIP_EXPR_SUM); /* only complex operands store their data as void* */

   return expr->data.data;
}

/** gives exponent belonging to a SCIP_EXPR_REALPOWER expression */
SCIP_Real SCIPexprGetRealPowerExponent(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_REALPOWER);

   return expr->data.dbl;
}

/** gives exponent belonging to a SCIP_EXPR_INTPOWER expression */
int SCIPexprGetIntPowerExponent(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_INTPOWER);

   return expr->data.intval;
}

/** gives exponent belonging to a SCIP_EXPR_SIGNPOWER expression */
SCIP_Real SCIPexprGetSignPowerExponent(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_SIGNPOWER);

   return expr->data.dbl;
}

/** gives linear coefficients belonging to a SCIP_EXPR_LINEAR expression */
SCIP_Real* SCIPexprGetLinearCoefs(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_LINEAR);
   assert(expr->data.data != NULL);

   /* the coefficients are stored in the first nchildren elements of the array stored as expression data */
   return (SCIP_Real*)expr->data.data;
}

/** gives constant belonging to a SCIP_EXPR_LINEAR expression */
SCIP_Real SCIPexprGetLinearConstant(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_LINEAR);
   assert(expr->data.data != NULL);

   /* the constant is stored in the nchildren's element of the array stored as expression data */
   return ((SCIP_Real*)expr->data.data)[expr->nchildren];
}

/** gives quadratic elements belonging to a SCIP_EXPR_QUADRATIC expression */
SCIP_QUADELEM* SCIPexprGetQuadElements(
   SCIP_EXPR*            expr                /**< quadratic expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_QUADRATIC);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_QUADRATIC*)expr->data.data)->quadelems;
}

/** gives constant belonging to a SCIP_EXPR_QUADRATIC expression */
SCIP_Real SCIPexprGetQuadConstant(
   SCIP_EXPR*            expr                /**< quadratic expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_QUADRATIC);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_QUADRATIC*)expr->data.data)->constant;
}

/** gives linear coefficients belonging to a SCIP_EXPR_QUADRATIC expression
 * can be NULL if all coefficients are 0.0 */
SCIP_Real* SCIPexprGetQuadLinearCoefs(
   SCIP_EXPR*            expr                /**< quadratic expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_QUADRATIC);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_QUADRATIC*)expr->data.data)->lincoefs;
}

/** gives number of quadratic elements belonging to a SCIP_EXPR_QUADRATIC expression */
int SCIPexprGetNQuadElements(
   SCIP_EXPR*            expr                /**< quadratic expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_QUADRATIC);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_QUADRATIC*)expr->data.data)->nquadelems;
}

/** gives the monomials belonging to a SCIP_EXPR_POLYNOMIAL expression */
SCIP_EXPRDATA_MONOMIAL** SCIPexprGetMonomials(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data)->monomials;
}

/** gives the number of monomials belonging to a SCIP_EXPR_POLYNOMIAL expression */
int SCIPexprGetNMonomials(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data)->nmonomials;
}

/** gives the constant belonging to a SCIP_EXPR_POLYNOMIAL expression */
SCIP_Real SCIPexprGetPolynomialConstant(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data)->constant;
}

/** gets coefficient of a monomial */
SCIP_Real SCIPexprGetMonomialCoef(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   )
{
   assert(monomial != NULL);

   return monomial->coef;
}

/** gets number of factors of a monomial */
int SCIPexprGetMonomialNFactors(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   )
{
   assert(monomial != NULL);

   return monomial->nfactors;
}

/** gets indices of children corresponding to factors of a monomial */
int* SCIPexprGetMonomialChildIndices(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   )
{
   assert(monomial != NULL);

   return monomial->childidxs;
}

/** gets exponents in factors of a monomial */
SCIP_Real* SCIPexprGetMonomialExponents(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   )
{
   assert(monomial != NULL);

   return monomial->exponents;
}

/** gets user data of a user expression */
SCIP_USEREXPRDATA* SCIPexprGetUserData(
   SCIP_EXPR*              expr
   )
{
   assert(expr != NULL);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_USER*)expr->data.data)->userdata;
}

/** indicates whether a user expression has the estimator callback defined */
SCIP_Bool SCIPexprHasUserEstimator(
   SCIP_EXPR*              expr
   )
{
   assert(expr != NULL);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_USER*)expr->data.data)->estimate != NULL;
}

/** gives the evaluation capability of a user expression */
SCIP_EXPRINTCAPABILITY SCIPexprGetUserEvalCapability(
   SCIP_EXPR*              expr
   )
{
   assert(expr != NULL);
   assert(expr->data.data != NULL);

   return ((SCIP_EXPRDATA_USER*)expr->data.data)->evalcapability;
}

/** creates a simple expression */
SCIP_RETCODE SCIPexprCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   SCIP_EXPROP           op,                 /**< operand of expression */
   ...                                       /**< arguments of operand */
   )
{
   va_list         ap;
   SCIP_EXPR**     children;
   SCIP_EXPROPDATA opdata;

   assert(blkmem != NULL);
   assert(expr   != NULL);

   switch( op )
   {
   case SCIP_EXPR_VARIDX:
   case SCIP_EXPR_PARAM:
   {
      va_start( ap, op );  /*lint !e838*/
      opdata.intval = va_arg( ap, int );  /*lint !e416 !e826*/
      va_end( ap );  /*lint !e826*/

      assert( opdata.intval >= 0 );

      SCIP_CALL( exprCreate( blkmem, expr, op, 0, NULL, opdata ) );
      break;
   }

   case SCIP_EXPR_CONST:
   {
      va_start(ap, op );  /*lint !e838*/
      opdata.dbl = va_arg( ap, SCIP_Real );  /*lint !e416 !e826*/
      va_end( ap );  /*lint !e826*/

      SCIP_CALL( exprCreate( blkmem, expr, op, 0, NULL, opdata ) );
      break;
   }

   /* operands with two children */
   case SCIP_EXPR_PLUS     :
   case SCIP_EXPR_MINUS    :
   case SCIP_EXPR_MUL      :
   case SCIP_EXPR_DIV      :
   case SCIP_EXPR_MIN      :
   case SCIP_EXPR_MAX      :
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &children, 2) );  /*lint !e506*/

      va_start(ap, op );  /*lint !e838*/
      children[0] = va_arg( ap, SCIP_EXPR* );  /*lint !e416 !e826*/
      children[1] = va_arg( ap, SCIP_EXPR* );  /*lint !e416 !e826*/
      assert(children[0] != NULL);
      assert(children[1] != NULL);
      va_end( ap );  /*lint !e826*/
      opdata.data = NULL; /* to avoid compiler warning about use of uninitialised value */

      SCIP_CALL( exprCreate( blkmem, expr, op, 2, children, opdata ) );
      break;
   }

   /* operands with one child */
   case SCIP_EXPR_SQUARE:
   case SCIP_EXPR_SQRT  :
   case SCIP_EXPR_EXP   :
   case SCIP_EXPR_LOG   :
   case SCIP_EXPR_SIN   :
   case SCIP_EXPR_COS   :
   case SCIP_EXPR_TAN   :
      /* case SCIP_EXPR_ERF   : */
      /* case SCIP_EXPR_ERFI  : */
   case SCIP_EXPR_ABS   :
   case SCIP_EXPR_SIGN  :
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &children, 1) );  /*lint !e506*/

      va_start(ap, op );  /*lint !e838*/
      children[0] = va_arg( ap, SCIP_EXPR* );  /*lint !e416 !e826*/
      assert(children[0] != NULL);
      va_end( ap );  /*lint !e826*/
      opdata.data = NULL; /* to avoid compiler warning about use of uninitialised value */

      SCIP_CALL( exprCreate( blkmem, expr, op, 1, children, opdata ) );
      break;
   }

   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_SIGNPOWER:
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &children, 1) );  /*lint !e506*/

      va_start(ap, op );  /*lint !e838*/
      children[0] = va_arg( ap, SCIP_EXPR* );  /*lint !e416 !e826*/
      assert(children[0] != NULL);
      opdata.dbl = va_arg( ap, SCIP_Real);  /*lint !e416 !e826*/
      va_end( ap );  /*lint !e826*/

      SCIP_CALL( exprCreate( blkmem, expr, op, 1, children, opdata ) );
      break;
   }

   case SCIP_EXPR_INTPOWER:
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &children, 1) );  /*lint !e506*/

      va_start(ap, op );  /*lint !e838*/
      children[0] = va_arg( ap, SCIP_EXPR* );  /*lint !e416 !e826*/
      assert(children[0] != NULL);
      opdata.intval = va_arg( ap, int);  /*lint !e416 !e826*/
      va_end( ap );  /*lint !e826*/

      SCIP_CALL( exprCreate( blkmem, expr, op, 1, children, opdata ) );
      break;
   }

   /* complex operands */
   case SCIP_EXPR_SUM    :
   case SCIP_EXPR_PRODUCT:
   {
      int nchildren;
      SCIP_EXPR** childrenarg;

      opdata.data = NULL; /* to avoid compiler warning about use of uninitialised value */

      va_start(ap, op );  /*lint !e838*/
      /* first argument should be number of children */
      nchildren = va_arg( ap, int );  /*lint !e416 !e826*/
      assert(nchildren >= 0);

      /* for a sum or product of 0 terms we can finish here */
      if( nchildren == 0 )
      {
         SCIP_RETCODE retcode;
         retcode = exprCreate( blkmem, expr, op, 0, NULL, opdata);
         va_end( ap );  /*lint !e826*/
         SCIP_CALL( retcode );
         break;
      }

      /* next argument should be array of children expressions */
      childrenarg = va_arg( ap, SCIP_EXPR** );  /*lint !e416 !e826*/
      assert(childrenarg != NULL);
      va_end( ap );  /*lint !e826*/

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &children, childrenarg, nchildren) );

      SCIP_CALL( exprCreate( blkmem, expr, op, nchildren, children, opdata) );
      break;
   }

   case SCIP_EXPR_LINEAR :
   case SCIP_EXPR_QUADRATIC:
   case SCIP_EXPR_POLYNOMIAL:
   case SCIP_EXPR_USER:
   {
      SCIPerrorMessage("cannot create complex expression linear, quadratic, polynomial, or user with SCIPexprCreate\n");
      return SCIP_INVALIDDATA;
   }

   case SCIP_EXPR_LAST:
      SCIPABORT();
      break;
   }

   return SCIP_OKAY;
}

/** copies an expression including its children */
SCIP_RETCODE SCIPexprCopyDeep(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to copied expression */
   SCIP_EXPR*            sourceexpr          /**< expression to copy */
   )
{
   assert(blkmem     != NULL);
   assert(targetexpr != NULL);
   assert(sourceexpr != NULL);

   SCIP_ALLOC( BMSduplicateBlockMemory(blkmem, targetexpr, sourceexpr) );

   if( sourceexpr->nchildren )
   {
      int i;

      /* alloc memory for children expressions */
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*targetexpr)->children, sourceexpr->nchildren) );

      /* copy children expressions */
      for( i = 0; i < sourceexpr->nchildren; ++i )
      {
         SCIP_CALL( SCIPexprCopyDeep(blkmem, &(*targetexpr)->children[i], sourceexpr->children[i]) );
      }
   }
   else
   {
      assert((*targetexpr)->children == NULL); /* otherwise, sourceexpr->children was not NULL, which is wrong */
   }

   /* call operands data copy callback for complex operands
    * for simple operands BMSduplicate above should have done the job
    */
   if( exprOpTable[sourceexpr->op].copydata != NULL )
   {
      SCIP_CALL( exprOpTable[sourceexpr->op].copydata(blkmem, sourceexpr->nchildren, sourceexpr->data, &(*targetexpr)->data) );
   }

   return SCIP_OKAY;
}

/** frees an expression including its children */
void SCIPexprFreeDeep(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr                /**< pointer to expression to free */
   )
{
   assert(blkmem != NULL);
   assert(expr   != NULL);
   assert(*expr  != NULL);

   /* call operands data free callback, if given */
   if( exprOpTable[(*expr)->op].freedata != NULL )
   {
      exprOpTable[(*expr)->op].freedata(blkmem, (*expr)->nchildren, (*expr)->data);
   }

   if( (*expr)->nchildren )
   {
      int i;

      assert( (*expr)->children != NULL );

      for( i = 0; i < (*expr)->nchildren; ++i )
      {
         SCIPexprFreeDeep(blkmem, &(*expr)->children[i]);
         assert((*expr)->children[i] == NULL);
      }

      BMSfreeBlockMemoryArray(blkmem, &(*expr)->children, (*expr)->nchildren);
   }
   else
   {
      assert( (*expr)->children == NULL );
   }

   BMSfreeBlockMemory(blkmem, expr);
}

/** frees an expression but not its children */
void SCIPexprFreeShallow(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr                /**< pointer to expression to free */
   )
{
   assert(blkmem != NULL);
   assert(expr   != NULL);
   assert(*expr  != NULL);

   /* call operands data free callback, if given */
   if( exprOpTable[(*expr)->op].freedata != NULL )
   {
      exprOpTable[(*expr)->op].freedata(blkmem, (*expr)->nchildren, (*expr)->data);
   }

   BMSfreeBlockMemoryArrayNull(blkmem, &(*expr)->children, (*expr)->nchildren);

   BMSfreeBlockMemory(blkmem, expr);
}

/** creates an expression from the addition of two given expression, with coefficients, and a constant
 *
 *  The given expressions may be modified or freed, otherwise it will be used a child expression.
 *  Favors creation and maintaining of SCIP_EXPR_LINEAR over SCIP_EXPR_PLUS or SCIP_EXPR_SUM.
 */
SCIP_RETCODE SCIPexprAdd(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to store pointer to created expression */
   SCIP_Real             coef1,              /**< coefficient of first term */
   SCIP_EXPR*            term1,              /**< expression of first term, or NULL */
   SCIP_Real             coef2,              /**< coefficient of second term */
   SCIP_EXPR*            term2,              /**< expression of second term, or NULL */
   SCIP_Real             constant            /**< constant term to add */
   )
{
   assert(blkmem != NULL);
   assert(expr != NULL);

   /* @todo could do something special with quadratic and polynomial expressions */

   if( term1 != NULL && SCIPexprGetOperator(term1) == SCIP_EXPR_CONST )
   {
      constant += coef1 * SCIPexprGetOpReal(term1);
      SCIPexprFreeDeep(blkmem, &term1);
   }

   if( term2 != NULL && SCIPexprGetOperator(term2) == SCIP_EXPR_CONST )
   {
      constant += coef2 * SCIPexprGetOpReal(term2);
      SCIPexprFreeDeep(blkmem, &term2);
   }

   if( term1 == NULL && term2 == NULL )
   {
      SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_CONST, constant) );
      return SCIP_OKAY;
   }

   if( term1 != NULL && SCIPexprGetOperator(term1) == SCIP_EXPR_LINEAR && coef1 != 1.0 )
   {
      /* multiply coefficients and constant of linear expression term1 by coef1 */
      SCIP_Real* data;
      int i;

      data = (SCIP_Real*)term1->data.data;
      assert(data != NULL);

      /* loop one more index to multiply also constant of linear expression */
      for( i = 0; i <= term1->nchildren; ++i )
         data[i] *= coef1;

      coef1 = 1.0;
   }

   if( term2 != NULL && SCIPexprGetOperator(term2) == SCIP_EXPR_LINEAR && coef2 != 1.0 )
   {
      /* multiply coefficients and constant of linear expression term2 by coef2 */
      SCIP_Real* data;
      int i;

      data = (SCIP_Real*)term2->data.data;
      assert(data != NULL);

      /* loop one more index to multiply also constant of linear expression */
      for( i = 0; i <= term2->nchildren; ++i )
         data[i] *= coef2;

      coef2 = 1.0;
   }

   if( term1 == NULL || term2 == NULL )
   {
      if( term1 == NULL )
      {
         term1 = term2;
         coef1 = coef2;
      }
      if( constant != 0.0 || coef1 != 1.0 )
      {
         if( SCIPexprGetOperator(term1) == SCIP_EXPR_LINEAR )
         {
            assert(coef1 == 1.0);

            /* add constant to existing linear expression */
            SCIP_CALL( SCIPexprAddToLinear(blkmem, term1, 0, NULL, NULL, constant) );
            *expr = term1;
         }
         else
         {
            /* create new linear expression for coef1 * term1 + constant */
            SCIP_CALL( SCIPexprCreateLinear(blkmem, expr, 1, &term1, &coef1, constant) );
         }
      }
      else
      {
         assert(constant == 0.0);
         assert(coef1 == 1.0);
         *expr = term1;
      }

      return SCIP_OKAY;
   }

   if( SCIPexprGetOperator(term1) == SCIP_EXPR_LINEAR && SCIPexprGetOperator(term2) == SCIP_EXPR_LINEAR )
   {
      /* add 2nd linear expression to first one */
      assert(coef1 == 1.0);
      assert(coef2 == 1.0);

      SCIP_CALL( SCIPexprAddToLinear(blkmem, term1, SCIPexprGetNChildren(term2), SCIPexprGetLinearCoefs(term2), SCIPexprGetChildren(term2), SCIPexprGetLinearConstant(term2) + constant) );
      SCIPexprFreeShallow(blkmem, &term2);

      *expr = term1;

      return SCIP_OKAY;
   }

   if( SCIPexprGetOperator(term2) == SCIP_EXPR_LINEAR )
   {
      /* if only term2 is linear, then swap */
      SCIP_EXPR* tmp;

      tmp = term2;
      assert(coef2 == 1.0);

      term2 = term1;
      coef2 = coef1;
      term1 = tmp;
      coef1 = 1.0;
   }

   if( SCIPexprGetOperator(term1) == SCIP_EXPR_LINEAR )
   {
      /* add coef2*term2 as extra child to linear expression term1 */
      assert(coef1 == 1.0);

      SCIP_CALL( SCIPexprAddToLinear(blkmem, term1, 1, &coef2, &term2, constant) );
      *expr = term1;

      return SCIP_OKAY;
   }

   /* both terms are not linear, then create new linear term for sum */
   {
      SCIP_Real coefs[2];
      SCIP_EXPR* children[2];

      coefs[0] = coef1;
      coefs[1] = coef2;
      children[0] = term1;
      children[1] = term2;

      SCIP_CALL( SCIPexprCreateLinear(blkmem, expr, 2, children, coefs, constant) );
   }

   return SCIP_OKAY;
}

/** creates an expression from the multiplication of an expression with a constant
 *
 *  The given expressions may be modified or freed, otherwise it will be used a child expression.
 *  Favors creation and maintaining SCIP_EXPR_LINEAR over SCIP_EXPR_PLUS or SCIP_EXPR_SUM.
 */
SCIP_RETCODE SCIPexprMulConstant(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to created expression */
   SCIP_EXPR*            term,               /**< term to multiply by factor */
   SCIP_Real             factor              /**< factor */
   )
{
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(term != NULL);

   if( factor == 0.0 )
   {
      SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_CONST, 0.0) );

      SCIPexprFreeDeep(blkmem, &term);

      return SCIP_OKAY;
   }
   if( factor == 1.0 )
   {
      *expr = term;
      return SCIP_OKAY;
   }

   switch( SCIPexprGetOperator(term) )
   {
      case SCIP_EXPR_CONST :
      {
         SCIP_CALL( SCIPexprCreate(blkmem, expr, SCIP_EXPR_CONST, factor * SCIPexprGetOpReal(term)) );
         SCIPexprFreeDeep(blkmem, &term);
         break;
      }

      case SCIP_EXPR_LINEAR :
      {
         SCIP_Real* data;
         int i;

         data = (SCIP_Real*)term->data.data;
         assert(data != NULL);

         /* loop one more index to multiply also constant of linear expression */
         for( i = 0; i <= SCIPexprGetNChildren(term); ++i )
            data[i] *= factor;

         *expr = term;
         break;
      }

      case SCIP_EXPR_QUADRATIC :
      {
         SCIP_EXPRDATA_QUADRATIC* data;
         int i;

         data = (SCIP_EXPRDATA_QUADRATIC*)term->data.data;

         data->constant *= factor;

         if( data->lincoefs != NULL )
            for( i = 0; i < term->nchildren; ++i )
               data->lincoefs[i] *= factor;

         for( i = 0; i < data->nquadelems; ++i )
            data->quadelems[i].coef *= factor;

         *expr = term;
         break;
      }

      case SCIP_EXPR_POLYNOMIAL :
      {
         SCIP_EXPRDATA_POLYNOMIAL* data;
         int i;

         data = (SCIP_EXPRDATA_POLYNOMIAL*)term->data.data;

         data->constant *= factor;

         for( i = 0; i < data->nmonomials; ++i )
            data->monomials[i]->coef *= factor;

         *expr = term;
         break;
      }

      default:
      {
         SCIP_CALL( SCIPexprCreateLinear(blkmem, expr, 1, &term, &factor, 0.0) );
         break;
      }

   } /*lint !e788 */

   return SCIP_OKAY;
}

/** creates a SCIP_EXPR_LINEAR expression that is (affine) linear in its children: constant + sum_i coef_i child_i */
SCIP_RETCODE SCIPexprCreateLinear(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children of expression */
   SCIP_Real*            coefs,              /**< coefficients of children */
   SCIP_Real             constant            /**< constant part */
   )
{
   SCIP_EXPROPDATA opdata;
   SCIP_EXPR**     childrencopy;
   SCIP_Real*      data;

   assert(nchildren >= 0);
   assert(children != NULL || nchildren == 0);
   assert(coefs    != NULL || nchildren == 0);

   if( nchildren > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &childrencopy, children, nchildren) );
   }
   else
      childrencopy = NULL;

   /* we store the coefficients and the constant in a single array and make this our operand data */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &data, nchildren + 1) );
   BMScopyMemoryArray(data, coefs, nchildren);  /*lint !e644*/
   data[nchildren] = constant;

   opdata.data = (void*)data;

   SCIP_CALL( exprCreate( blkmem, expr, SCIP_EXPR_LINEAR, nchildren, childrencopy, opdata) );

   return SCIP_OKAY;
}

/** adds new terms to a linear expression */
SCIP_RETCODE SCIPexprAddToLinear(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< linear expression */
   int                   nchildren,          /**< number of children to add */
   SCIP_Real*            coefs,              /**< coefficients of additional children */
   SCIP_EXPR**           children,           /**< additional children expressions */
   SCIP_Real             constant            /**< constant to add */
   )
{
   SCIP_Real* data;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_LINEAR);
   assert(nchildren >= 0);
   assert(coefs != NULL || nchildren == 0);
   assert(children != NULL || nchildren == 0);

   data = (SCIP_Real*)expr->data.data;
   assert(data != NULL);

   /* handle simple case of adding a constant */
   if( nchildren == 0 )
   {
      data[expr->nchildren] += constant;

      return SCIP_OKAY;
   }

   /* add new children to expr's children array */
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &expr->children, expr->nchildren, expr->nchildren + nchildren) );
   BMScopyMemoryArray(&expr->children[expr->nchildren], children, nchildren);  /*lint !e866*/

   /* add constant and new coefs to expr's data array */
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &data, expr->nchildren + 1, expr->nchildren + nchildren + 1) );
   data[expr->nchildren + nchildren] = data[expr->nchildren] + constant;
   BMScopyMemoryArray(&data[expr->nchildren], coefs, nchildren);  /*lint !e866*/
   expr->data.data = (void*)data;

   expr->nchildren += nchildren;

   return SCIP_OKAY;
}

/** creates a SCIP_EXPR_QUADRATIC expression: constant + sum_i coef_i child_i + sum_i coef_i child1_i child2_i */
SCIP_RETCODE SCIPexprCreateQuadratic(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children of expression */
   SCIP_Real             constant,           /**< constant */
   SCIP_Real*            lincoefs,           /**< linear coefficients of children, or NULL if all 0.0 */
   int                   nquadelems,         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems           /**< quadratic elements specifying coefficients and child indices */
   )
{
   SCIP_EXPROPDATA opdata;
   SCIP_EXPR**     childrencopy;
   SCIP_EXPRDATA_QUADRATIC* data;

   assert(nchildren >= 0);
   assert(children  != NULL || nchildren == 0);
   assert(quadelems != NULL || nquadelems == 0);

   if( nchildren > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &childrencopy, children, nchildren) );
   }
   else
      childrencopy = NULL;

   SCIP_CALL( quadraticdataCreate(blkmem, &data, constant, nchildren, lincoefs, nquadelems, quadelems) );

   opdata.data = (void*)data;

   SCIP_CALL( exprCreate( blkmem, expr, SCIP_EXPR_QUADRATIC, nchildren, childrencopy, opdata) );

   return SCIP_OKAY;
}

/** ensures that quadratic elements of a quadratic expression are sorted */
void SCIPexprSortQuadElems(
   SCIP_EXPR*            expr                /**< quadratic expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_QUADRATIC);
   assert(expr->data.data != NULL);

   quadraticdataSort((SCIP_EXPRDATA_QUADRATIC*)expr->data.data);
}

/** creates a SCIP_EXPR_POLYNOMIAL expression from an array of monomials: constant + sum_i monomial_i */
SCIP_RETCODE SCIPexprCreatePolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children of expression */
   int                   nmonomials,         /**< number of monomials */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< monomials */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool             copymonomials       /**< should monomials by copied or ownership be assumed? */
   )
{
   SCIP_EXPROPDATA opdata;
   SCIP_EXPR**     childrencopy;
   SCIP_EXPRDATA_POLYNOMIAL* data;

   assert(nchildren >= 0);
   assert(children != NULL || nchildren == 0);
   assert(monomials   != NULL || nmonomials   == 0);

   if( nchildren > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &childrencopy, children, nchildren) );
   }
   else
      childrencopy = NULL;

   SCIP_CALL( polynomialdataCreate(blkmem, &data, nmonomials, monomials, constant, copymonomials) );
   opdata.data = (void*)data;

   SCIP_CALL( exprCreate( blkmem, expr, SCIP_EXPR_POLYNOMIAL, nchildren, childrencopy, opdata) );

   return SCIP_OKAY;
}

/** adds an array of monomials to a SCIP_EXPR_POLYNOMIAL expression */
SCIP_RETCODE SCIPexprAddMonomials(
   BMS_BLKMEM*           blkmem,             /**< block memory of expression */
   SCIP_EXPR*            expr,               /**< expression */
   int                   nmonomials,         /**< number of monomials to add */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< the monomials to add */
   SCIP_Bool             copymonomials       /**< should monomials by copied or ownership be assumed? */
   )
{
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(monomials != NULL || nmonomials == 0);

   if( nmonomials == 0 )
      return SCIP_OKAY;

   SCIP_CALL( polynomialdataAddMonomials(blkmem, (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data, nmonomials, monomials, copymonomials) );

   return SCIP_OKAY;
}

/** changes the constant in a SCIP_EXPR_POLYNOMIAL expression */
void SCIPexprChgPolynomialConstant(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             constant            /**< new value for constant */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   ((SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data)->constant = constant;
}

/** multiplies each summand of a polynomial by a given constant */
void SCIPexprMultiplyPolynomialByConstant(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_Real             factor              /**< constant factor */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   polynomialdataMultiplyByConstant(blkmem, (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data, factor);
}

/** multiplies each summand of a polynomial by a given monomial */
SCIP_RETCODE SCIPexprMultiplyPolynomialByMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_EXPRDATA_MONOMIAL*  factor,          /**< monomial factor */
   int*                  childmap            /**< map children in factor to children in expr, or NULL for 1:1 */
   )
{
   assert(blkmem != NULL);
   assert(factor != NULL);
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   SCIP_CALL( polynomialdataMultiplyByMonomial(blkmem, (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data, factor, childmap) );

   return SCIP_OKAY;
}

/** multiplies this polynomial by a polynomial
 *
 *  Factor needs to be different from expr.
 *  Children of factor need to be children of expr already, w.r.t. an optional mapping of child indices.
 */
SCIP_RETCODE SCIPexprMultiplyPolynomialByPolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_EXPR*            factor,             /**< polynomial factor */
   int*                  childmap            /**< map children in factor to children in expr, or NULL for 1:1 */
   )
{
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);
   assert(factor != NULL);
   assert(factor->op == SCIP_EXPR_POLYNOMIAL);
   assert(factor->data.data != NULL);
   assert(expr != factor);

#ifndef NDEBUG
   if( childmap == NULL )
   {
      int i;
      assert(factor->nchildren == expr->nchildren);
      for( i = 0; i < factor->nchildren; ++i )
         assert(SCIPexprAreEqual(expr->children[i], factor->children[i], 0.0));
   }
   else
   {
      int i;
      for( i = 0; i < factor->nchildren; ++i )
      {
         assert(childmap[i] >= 0);
         assert(childmap[i] < expr->nchildren);
         assert(SCIPexprAreEqual(expr->children[childmap[i]], factor->children[i], 0.0));
      }
   }
#endif

   SCIP_CALL( polynomialdataMultiplyByPolynomial(blkmem, (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data, (SCIP_EXPRDATA_POLYNOMIAL*)factor->data.data, childmap) );

   return SCIP_OKAY;
}

/** takes a power of the polynomial
 *
 *  Exponent need to be an integer.
 *  Polynomial needs to be a monomial, if exponent is negative.
 */
SCIP_RETCODE SCIPexprPolynomialPower(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   int                   exponent            /**< exponent of power operation */
   )
{
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   SCIP_CALL( polynomialdataPower(blkmem, (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data, exponent) );

   return SCIP_OKAY;
}

/** merges monomials in a polynomial expression that differ only in coefficient into a single monomial
 *
 *  Eliminates monomials with coefficient between -eps and eps.
 */
void SCIPexprMergeMonomials(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_Real             eps,                /**< threshold under which numbers are treat as zero */
   SCIP_Bool             mergefactors        /**< whether to merge factors in monomials too */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   polynomialdataMergeMonomials(blkmem, (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data, eps, mergefactors);
}

/** checks if two monomials are equal */
SCIP_Bool SCIPexprAreMonomialsEqual(
   SCIP_EXPRDATA_MONOMIAL* monomial1,        /**< first monomial */
   SCIP_EXPRDATA_MONOMIAL* monomial2,        /**< second monomial */
   SCIP_Real             eps                 /**< threshold under which numbers are treated as 0.0 */
   )
{
   int i;

   assert(monomial1 != NULL);
   assert(monomial2 != NULL);

   if( monomial1->nfactors != monomial2->nfactors )
      return FALSE;

   if( !EPSEQ(monomial1->coef, monomial2->coef, eps) )
      return FALSE;

   SCIPexprSortMonomialFactors(monomial1);
   SCIPexprSortMonomialFactors(monomial2);

   for( i = 0; i < monomial1->nfactors; ++i )
   {
      if( monomial1->childidxs[i] != monomial2->childidxs[i] ||
         !EPSEQ(monomial1->exponents[i], monomial2->exponents[i], eps) )
         return FALSE;
   }

   return TRUE;
}

/** changes coefficient of monomial */
void SCIPexprChgMonomialCoef(
   SCIP_EXPRDATA_MONOMIAL*  monomial,              /**< monomial */
   SCIP_Real             newcoef             /**< new coefficient */
   )
{
   assert(monomial != NULL);

   monomial->coef = newcoef;
}

/** adds factors to a monomial */
SCIP_RETCODE SCIPexprAddMonomialFactors(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   int                   nfactors,           /**< number of factors to add */
   int*                  childidxs,          /**< indices of children corresponding to factors */
   SCIP_Real*            exponents           /**< exponent in each factor */
   )
{
   assert(monomial != NULL);
   assert(nfactors >= 0);
   assert(childidxs != NULL || nfactors == 0);
   assert(exponents != NULL || nfactors == 0);

   if( nfactors == 0 )
      return SCIP_OKAY;

   SCIP_CALL( monomialdataEnsureFactorsSize(blkmem, monomial, monomial->nfactors + nfactors) );
   assert(monomial->nfactors + nfactors <= monomial->factorssize);

   BMScopyMemoryArray(&monomial->childidxs[monomial->nfactors], childidxs, nfactors);  /*lint !e866*/
   BMScopyMemoryArray(&monomial->exponents[monomial->nfactors], exponents, nfactors);  /*lint !e866*/

   monomial->nfactors += nfactors;
   monomial->sorted = (monomial->nfactors <= 1);

   return SCIP_OKAY;
}

/** multiplies a monomial with a monomial */
SCIP_RETCODE SCIPexprMultiplyMonomialByMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   SCIP_EXPRDATA_MONOMIAL* factor,           /**< factor monomial */
   int*                  childmap            /**< map to apply to children in factor, or NULL for 1:1 */
   )
{
   assert(monomial != NULL);
   assert(factor != NULL);

   if( factor->coef == 0.0 )
   {
      monomial->nfactors = 0;
      monomial->coef = 0.0;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPexprAddMonomialFactors(blkmem, monomial, factor->nfactors, factor->childidxs, factor->exponents) );

   if( childmap != NULL )
   {
      int i;
      for( i = monomial->nfactors - factor->nfactors; i < monomial->nfactors; ++i )
         monomial->childidxs[i] = childmap[monomial->childidxs[i]];
   }

   monomial->coef *= factor->coef;

   return SCIP_OKAY;
}

/** replaces the monomial by a power of the monomial
 *
 *  Allows only integers as exponent.
 */
void SCIPexprMonomialPower(
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   int                   exponent            /**< integer exponent of power operation */
   )
{
   int i;

   assert(monomial != NULL);

   if( exponent == 1 )
      return;

   if( exponent == 0 )
   {
      /* x^0 = 1, unless x = 0; 0^0 = 0 */
      if( monomial->coef != 0.0 )
         monomial->coef = 1.0;
      monomial->nfactors = 0;
      return;
   }

   monomial->coef = pow(monomial->coef, (SCIP_Real)exponent);
   for( i = 0; i < monomial->nfactors; ++i )
      monomial->exponents[i] *= exponent;
}

/** merges factors that correspond to the same child by adding exponents
 *
 *  Eliminates factors with exponent between -eps and eps.
 */
void SCIPexprMergeMonomialFactors(
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   SCIP_Real             eps                 /**< threshold under which numbers are treated as 0.0 */
   )
{
   int i;
   int offset;

   assert(monomial != NULL);
   assert(eps >= 0.0);

   SCIPexprSortMonomialFactors(monomial);

   /* merge factors with same child index by adding up their exponents
    * delete factors with exponent 0.0 */
   offset = 0;
   i = 0;
   while( i + offset < monomial->nfactors )
   {
      if( offset > 0 )
      {
         assert(monomial->childidxs[i] == -1);
         assert(monomial->childidxs[i+offset] >= 0);
         monomial->childidxs[i] = monomial->childidxs[i+offset];
         monomial->exponents[i] = monomial->exponents[i+offset];
#ifndef NDEBUG
         monomial->childidxs[i+offset] = -1;
#endif
      }

      while( i+offset+1 < monomial->nfactors && monomial->childidxs[i] == monomial->childidxs[i+offset+1] )
      {
         monomial->exponents[i] += monomial->exponents[i+offset+1];
#ifndef NDEBUG
         monomial->childidxs[i+offset+1] = -1;
#endif
         ++offset;
      }

      if( EPSZ(monomial->exponents[i], eps) )
      {
#ifndef NDEBUG
         monomial->childidxs[i] = -1;
#endif
         ++offset;
         continue;
      }
      else if( EPSISINT(monomial->exponents[i], eps) )
         monomial->exponents[i] = EPSROUND(monomial->exponents[i], eps);

      ++i;
   }

#ifndef NDEBUG
   for( ; i < monomial->nfactors; ++i )
      assert(monomial->childidxs[i] == -1);
#endif

   monomial->nfactors -= offset;

   if( EPSEQ(monomial->coef, 1.0, eps) )
      monomial->coef = 1.0;
   else if( EPSEQ(monomial->coef, -1.0, eps) )
      monomial->coef = -1.0;
}

/** ensures that monomials of a polynomial are sorted */
void SCIPexprSortMonomials(
   SCIP_EXPR*            expr                /**< polynomial expression */
   )
{
   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_POLYNOMIAL);
   assert(expr->data.data != NULL);

   polynomialdataSortMonomials((SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data);
}

/** creates a monomial */
SCIP_RETCODE SCIPexprCreateMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_MONOMIAL** monomial,        /**< buffer where to store pointer to new monomial */
   SCIP_Real             coef,               /**< coefficient of monomial */
   int                   nfactors,           /**< number of factors in monomial */
   int*                  childidxs,          /**< indices of children corresponding to factors, or NULL if identity */
   SCIP_Real*            exponents           /**< exponent in each factor, or NULL if all 1.0 */
   )
{
   assert(blkmem != NULL);
   assert(monomial != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, monomial) );

   (*monomial)->coef     = coef;
   (*monomial)->nfactors = nfactors;
   (*monomial)->factorssize = nfactors;
   (*monomial)->sorted = (nfactors <= 1);

   if( nfactors > 0 )
   {
      if( childidxs != NULL )
      {
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*monomial)->childidxs, childidxs, nfactors) );
      }
      else
      {
         int i;

         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*monomial)->childidxs, nfactors) );
         for( i = 0; i < nfactors; ++i )
            (*monomial)->childidxs[i] = i;
      }

      if( exponents != NULL )
      {
         SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*monomial)->exponents, exponents, nfactors) );
      }
      else
      {
         int i;

         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*monomial)->exponents, nfactors) );
         for( i = 0; i < nfactors; ++i )
            (*monomial)->exponents[i] = 1.0;
      }
   }
   else
   {
      (*monomial)->childidxs = NULL;
      (*monomial)->exponents = NULL;
   }

   return SCIP_OKAY;
}

/** frees a monomial */
void SCIPexprFreeMonomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRDATA_MONOMIAL** monomial         /**< pointer to monomial that should be freed */
   )
{
   assert(blkmem != NULL);
   assert( monomial != NULL);
   assert(*monomial != NULL);

   if( (*monomial)->factorssize > 0 )
   {
      assert((*monomial)->childidxs != NULL);
      assert((*monomial)->exponents != NULL);

      BMSfreeBlockMemoryArray(blkmem, &(*monomial)->childidxs, (*monomial)->factorssize);
      BMSfreeBlockMemoryArray(blkmem, &(*monomial)->exponents, (*monomial)->factorssize);
   }
   assert((*monomial)->childidxs == NULL);
   assert((*monomial)->exponents == NULL);

   BMSfreeBlockMemory(blkmem, monomial);
}

/** ensures that factors in a monomial are sorted */
void SCIPexprSortMonomialFactors(
   SCIP_EXPRDATA_MONOMIAL* monomial          /**< monomial */
   )
{
   assert(monomial != NULL);

   if( monomial->sorted )
      return;

   if( monomial->nfactors > 0 )
      SCIPsortIntReal(monomial->childidxs, monomial->exponents, monomial->nfactors);

   monomial->sorted = TRUE;
}

/** finds a factor corresponding to a given child index in a monomial
 *
 *  Note that if the factors have not been merged, the position of some factor corresponding to a given child is given.
 *  Returns TRUE if a factor is found, FALSE if not.
 */
SCIP_Bool SCIPexprFindMonomialFactor(
   SCIP_EXPRDATA_MONOMIAL* monomial,         /**< monomial */
   int                   childidx,           /**< index of the child which factor to search for */
   int*                  pos                 /**< buffer to store position of factor */
   )
{
   assert(monomial != NULL);

   if( monomial->nfactors == 0 )
      return FALSE;

   SCIPexprSortMonomialFactors(monomial);

   return SCIPsortedvecFindInt(monomial->childidxs, childidx, monomial->nfactors, pos);
}

/** creates a user expression */
SCIP_RETCODE SCIPexprCreateUser(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR**           expr,               /**< pointer to buffer for expression address */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children of expression */
   SCIP_USEREXPRDATA*    data,               /**< user data for expression, expression assumes ownership */
   SCIP_EXPRINTCAPABILITY evalcapability,    /**< capability of evaluation functions (partially redundant, currently) */
   SCIP_DECL_USEREXPREVAL    ((*eval)),      /**< evaluation function */
   SCIP_DECL_USEREXPRINTEVAL ((*inteval)),   /**< interval evaluation function, or NULL if not implemented */
   SCIP_DECL_USEREXPRCURV    ((*curv)),      /**< curvature check function */
   SCIP_DECL_USEREXPRPROP    ((*prop)),      /**< interval propagation function, or NULL if not implemented */
   SCIP_DECL_USEREXPRESTIMATE ((*estimate)), /**< estimation function, or NULL if convex, concave, or not implemented */
   SCIP_DECL_USEREXPRCOPYDATA ((*copydata)), /**< expression data copy function, or NULL if nothing to copy */
   SCIP_DECL_USEREXPRFREEDATA ((*freedata)), /**< expression data free function, or NULL if nothing to free */
   SCIP_DECL_USEREXPRPRINT ((*print))        /**< expression print function, or NULL for default string "user" */
   )
{
   SCIP_EXPROPDATA opdata;
   SCIP_EXPRDATA_USER* userexprdata;
   SCIP_EXPR** childrencopy;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(eval != NULL);
   assert((evalcapability & SCIP_EXPRINTCAPABILITY_FUNCVALUE) != 0);  /* the function evaluation is not optional */
   assert(((evalcapability & SCIP_EXPRINTCAPABILITY_INTFUNCVALUE) == 0) || inteval != NULL);  /* if capability says it can do interval evaluation, then the corresponding callback needs to be provided */
   assert(curv != NULL);
   assert(copydata != NULL || data == NULL);
   assert(freedata != NULL || data == NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &userexprdata) );

   userexprdata->userdata = data;
   userexprdata->evalcapability = evalcapability;
   userexprdata->eval = eval;
   userexprdata->inteval = inteval;
   userexprdata->curv = curv;
   userexprdata->prop = prop;
   userexprdata->estimate = estimate;
   userexprdata->copydata = copydata;
   userexprdata->freedata = freedata;
   userexprdata->print = print;

   opdata.data = (void*) userexprdata;

   if( nchildren == 0 )
   {
      SCIP_CALL( exprCreate(blkmem, expr, SCIP_EXPR_USER, 0, NULL, opdata) );
      return SCIP_OKAY;
   }
   assert(children != NULL);

   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &childrencopy, children, nchildren) );

   SCIP_CALL( exprCreate( blkmem, expr, SCIP_EXPR_USER, nchildren, childrencopy, opdata) );

   return SCIP_OKAY;
}

/** indicates whether the expression contains a SCIP_EXPR_PARAM */
SCIP_Bool SCIPexprHasParam(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int i;

   assert(expr != NULL);

   if( expr->op == SCIP_EXPR_PARAM )
      return TRUE;

   for( i = 0; i < expr->nchildren; ++i )
      if( SCIPexprHasParam(expr->children[i]) )
         return TRUE;

   return FALSE;
}

/** gets maximal degree of expression, or SCIP_EXPR_DEGREEINFINITY if not a polynomial */
SCIP_RETCODE SCIPexprGetMaxDegree(
   SCIP_EXPR*            expr,               /**< expression */
   int*                  maxdegree           /**< buffer to store maximal degree */
   )
{
   int child1;
   int child2;

   assert(expr      != NULL);
   assert(maxdegree != NULL);

   switch( expr->op )
   {
   case SCIP_EXPR_VARIDX:
      *maxdegree = 1;
      break;

   case SCIP_EXPR_CONST:
   case SCIP_EXPR_PARAM:
      *maxdegree = 0;
      break;

   case SCIP_EXPR_PLUS:
   case SCIP_EXPR_MINUS:
   {
      assert(expr->children[0] != NULL);
      assert(expr->children[1] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );

      *maxdegree = MAX(child1, child2);
      break;
   }

   case SCIP_EXPR_MUL:
   {
      assert(expr->children[0] != NULL);
      assert(expr->children[1] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );

      *maxdegree = child1 + child2;
      break;
   }

   case SCIP_EXPR_DIV:
   {
      assert(expr->children[0] != NULL);
      assert(expr->children[1] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );

      /* if not division by constant, then it is not a polynomial */
      *maxdegree = (child2 != 0) ? SCIP_EXPR_DEGREEINFINITY : child1;
      break;
   }

   case SCIP_EXPR_SQUARE:
   {
      assert(expr->children[0] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );

      *maxdegree = 2 * child1;
      break;
   }

   case SCIP_EXPR_SQRT:
   {
      assert(expr->children[0] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );

      /* if not squareroot of constant, then no polynomial */
      *maxdegree = (child1 != 0) ? SCIP_EXPR_DEGREEINFINITY : 0;
      break;
   }

   case SCIP_EXPR_REALPOWER:
   {
      assert(expr->children[0] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );

      /* constant ^ constant has degree 0 */
      if( child1 == 0 )
      {
         *maxdegree = 0;
         break;
      }

      /* non-polynomial ^ constant is not a polynomial */
      if( child1 >= SCIP_EXPR_DEGREEINFINITY )
      {
         *maxdegree = SCIP_EXPR_DEGREEINFINITY;
         break;
      }

      /* so it is polynomial ^ constant
       * let's see whether the constant is integral */

      if( expr->data.dbl == 0.0 ) /* polynomial ^ 0 == 0 */
         *maxdegree = 0;
      else if( expr->data.dbl > 0.0 && (int)expr->data.dbl == expr->data.dbl ) /* natural exponent gives polynomial again */  /*lint !e777*/
         *maxdegree = child1 * (int)expr->data.dbl;
      else /* negative or nonintegral exponent does not give polynomial */
         *maxdegree = SCIP_EXPR_DEGREEINFINITY;

      break;
   }

   case SCIP_EXPR_INTPOWER:
   {
      assert(expr->children[0] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );

      /* constant ^ integer or something ^ 0 has degree 0 */
      if( child1 == 0 || expr->data.intval == 0 )
      {
         *maxdegree = 0;
         break;
      }

      /* non-polynomial ^ integer  or  something ^ negative  is not a polynomial */
      if( child1 >= SCIP_EXPR_DEGREEINFINITY || expr->data.intval < 0 )
      {
         *maxdegree = SCIP_EXPR_DEGREEINFINITY;
         break;
      }

      /* so it is polynomial ^ natural, which gives a polynomial again */
      *maxdegree = child1 * expr->data.intval;

      break;
   }

   case SCIP_EXPR_SIGNPOWER:
   {
      assert(expr->children[0] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );

      /* if child is not constant, then it is no polynomial */
      *maxdegree = child1 != 0 ? SCIP_EXPR_DEGREEINFINITY : 0;
      break;
   }

   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   case SCIP_EXPR_USER:
   {
      assert(expr->children[0] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );

      /* if argument is not a constant, then no polynomial, otherwise it is a constant */
      *maxdegree = (child1 != 0) ? SCIP_EXPR_DEGREEINFINITY : 0;
      break;
   }

   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   {
      assert(expr->children[0] != NULL);
      assert(expr->children[1] != NULL);

      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[0], &child1) );
      SCIP_CALL( SCIPexprGetMaxDegree(expr->children[1], &child2) );

      /* if any of the operands is not constant, then it is no polynomial */
      *maxdegree = (child1 != 0 || child2 != 0) ? SCIP_EXPR_DEGREEINFINITY : 0;
      break;
   }

   case SCIP_EXPR_SUM:
   case SCIP_EXPR_LINEAR:
   {
      int i;

      *maxdegree = 0;
      for( i = 0; i < expr->nchildren && *maxdegree < SCIP_EXPR_DEGREEINFINITY; ++i )
      {
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[i], &child1) );
         if( child1 > *maxdegree )
            *maxdegree = child1;
      }

      break;
   }

   case SCIP_EXPR_PRODUCT:
   {
      int i;

      *maxdegree = 0;
      for( i = 0; i < expr->nchildren; ++i )
      {
         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[i], &child1) );
         if( child1 >= SCIP_EXPR_DEGREEINFINITY )
         {
            *maxdegree = SCIP_EXPR_DEGREEINFINITY;
            break;
         }
         *maxdegree += child1;
      }

      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* quadraticdata;
      int childidx;
      int quadidx;

      quadraticdata = (SCIP_EXPRDATA_QUADRATIC*)expr->data.data;

      /* make sure quadratic elements are sorted */
      quadraticdataSort(quadraticdata);

      *maxdegree = 0;
      quadidx = 0;
      for( childidx = 0; childidx < expr->nchildren; ++childidx )
      {
         /* if no linear or no quadratic coefficient with current child on first position, then nothing to do */
         if( (quadraticdata->lincoefs == NULL || quadraticdata->lincoefs[childidx] == 0.0) &&
            (quadidx < quadraticdata->nquadelems && quadraticdata->quadelems[quadidx].idx1 > childidx) )
            continue;

         SCIP_CALL( SCIPexprGetMaxDegree(expr->children[childidx], &child1) );
         if( child1 == SCIP_EXPR_DEGREEINFINITY )
         {
            *maxdegree = SCIP_EXPR_DEGREEINFINITY;
            break;
         }

         while( quadidx < quadraticdata->nquadelems && quadraticdata->quadelems[quadidx].idx1 == childidx )
         {
            if( quadraticdata->quadelems[quadidx].idx2 == childidx )
            {
               /* square term */
               if( 2*child1 > *maxdegree )
                  *maxdegree = 2*child1;
            }
            else
            {
               /* bilinear term */
               SCIP_CALL( SCIPexprGetMaxDegree(expr->children[quadraticdata->quadelems[quadidx].idx2], &child2) );
               if( child2 == SCIP_EXPR_DEGREEINFINITY )
               {
                  *maxdegree = SCIP_EXPR_DEGREEINFINITY;
                  break;
               }
               if( child1 + child2 > *maxdegree )
                  *maxdegree = child1 + child2;
            }
            ++quadidx;
         }
         if( *maxdegree == SCIP_EXPR_DEGREEINFINITY )
            break;
      }

      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomialdata;
      int monomialdegree;
      int i;
      int j;

      polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data;

      *maxdegree = 0;
      for( i = 0; i < polynomialdata->nmonomials && *maxdegree < SCIP_EXPR_DEGREEINFINITY; ++i )
      {
         monomialdata = polynomialdata->monomials[i];
         assert(monomialdata != NULL);

         /* compute degree of monomial = sum of degree of factors */
         monomialdegree = 0;
         for( j = 0; j < monomialdata->nfactors; ++j )
         {
            SCIP_CALL( SCIPexprGetMaxDegree(expr->children[monomialdata->childidxs[j]], &child1) );

            /* if the exponent of the factor is not a natural number and the child is not constant (degree 0),
             * then we report that we are not really a polynomial */
            if( child1 != 0 && (monomialdata->exponents[j] < 0.0 || (int)monomialdata->exponents[j] != monomialdata->exponents[j]) )
            {
               *maxdegree = SCIP_EXPR_DEGREEINFINITY;
               break;
            }

            monomialdegree += child1 * (int)monomialdata->exponents[j];
         }

         if( monomialdegree > *maxdegree )
            *maxdegree = monomialdegree;
      }

      break;
   }

   case SCIP_EXPR_LAST:
      SCIPABORT();
      break;
   }

   return SCIP_OKAY;
}

/** counts usage of variables in expression */
void SCIPexprGetVarsUsage(
   SCIP_EXPR*            expr,               /**< expression to update */
   int*                  varsusage           /**< array with counters of variable usage */
   )
{
   int i;

   assert(expr != NULL);
   assert(varsusage != NULL);

   if( expr->op == SCIP_EXPR_VARIDX )
   {
      ++varsusage[expr->data.intval];
   }

   for( i = 0; i < expr->nchildren; ++i )
      SCIPexprGetVarsUsage(expr->children[i], varsusage);
}

/** compares whether two expressions are the same
 *
 *  Inconclusive, i.e., may give FALSE even if expressions are equivalent (x*y != y*x).
 */
SCIP_Bool SCIPexprAreEqual(
   SCIP_EXPR*            expr1,              /**< first expression */
   SCIP_EXPR*            expr2,              /**< second expression */
   SCIP_Real             eps                 /**< threshold under which numbers are assumed to be zero */
   )
{
   assert(expr1 != NULL);
   assert(expr2 != NULL);

   if( expr1 == expr2 )
      return TRUE;

   if( expr1->op != expr2->op )
      return FALSE;

   switch( expr1->op )
   {
   case SCIP_EXPR_VARIDX:
   case SCIP_EXPR_PARAM:
      return expr1->data.intval == expr2->data.intval;

   case SCIP_EXPR_CONST:
      return EPSEQ(expr1->data.dbl, expr2->data.dbl, eps);

      /* operands with two children */
   case SCIP_EXPR_PLUS     :
   case SCIP_EXPR_MINUS    :
   case SCIP_EXPR_MUL      :
   case SCIP_EXPR_DIV      :
   case SCIP_EXPR_MIN      :
   case SCIP_EXPR_MAX      :
      return SCIPexprAreEqual(expr1->children[0], expr2->children[0], eps) && SCIPexprAreEqual(expr1->children[1], expr2->children[1], eps);

      /* operands with one child */
   case SCIP_EXPR_SQUARE:
   case SCIP_EXPR_SQRT  :
   case SCIP_EXPR_EXP   :
   case SCIP_EXPR_LOG   :
   case SCIP_EXPR_SIN   :
   case SCIP_EXPR_COS   :
   case SCIP_EXPR_TAN   :
      /* case SCIP_EXPR_ERF   : */
      /* case SCIP_EXPR_ERFI  : */
   case SCIP_EXPR_ABS   :
   case SCIP_EXPR_SIGN  :
      return SCIPexprAreEqual(expr1->children[0], expr2->children[0], eps);

   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_SIGNPOWER:
      return EPSEQ(expr1->data.dbl, expr2->data.dbl, eps) && SCIPexprAreEqual(expr1->children[0], expr2->children[0], eps);

   case SCIP_EXPR_INTPOWER:
      return expr1->data.intval == expr2->data.intval && SCIPexprAreEqual(expr1->children[0], expr2->children[0], eps);

      /* complex operands */
   case SCIP_EXPR_SUM    :
   case SCIP_EXPR_PRODUCT:
   {
      int i;

      /* @todo sort children and have sorted flag in data? */

      if( expr1->nchildren != expr2->nchildren )
         return FALSE;

      for( i = 0; i < expr1->nchildren; ++i )
      {
         if( !SCIPexprAreEqual(expr1->children[i], expr2->children[i], eps) )
            return FALSE;
      }

      return TRUE;
   }

   case SCIP_EXPR_LINEAR :
   {
      SCIP_Real* data1;
      SCIP_Real* data2;
      int i;

      /* @todo sort children and have sorted flag in data? */

      if( expr1->nchildren != expr2->nchildren )
         return FALSE;

      data1 = (SCIP_Real*)expr1->data.data;
      data2 = (SCIP_Real*)expr2->data.data;

      /* check if constant and coefficients are equal */
      for( i = 0; i < expr1->nchildren + 1; ++i )
         if( !EPSEQ(data1[i], data2[i], eps) )
            return FALSE;

      /* check if children are equal */
      for( i = 0; i < expr1->nchildren; ++i )
      {
         if( !SCIPexprAreEqual(expr1->children[i], expr2->children[i], eps) )
            return FALSE;
      }

      return TRUE;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* data1;
      SCIP_EXPRDATA_QUADRATIC* data2;
      int i;

      if( expr1->nchildren != expr2->nchildren )
         return FALSE;

      data1 = (SCIP_EXPRDATA_QUADRATIC*)expr1->data.data;
      data2 = (SCIP_EXPRDATA_QUADRATIC*)expr2->data.data;

      if( data1->nquadelems != data2->nquadelems )
         return FALSE;

      if( !EPSEQ(data1->constant, data2->constant, eps) )
         return FALSE;

      /* check if linear part is equal */
      if( data1->lincoefs != NULL || data2->lincoefs != NULL )
         for( i = 0; i < expr1->nchildren; ++i )
         {
            if( data1->lincoefs == NULL )
            {
               if( !EPSZ(data2->lincoefs[i], eps) )
                  return FALSE;
            }
            else if( data2->lincoefs == NULL )
            {
               if( !EPSZ(data1->lincoefs[i], eps) )
                  return FALSE;
            }
            else if( !EPSEQ(data1->lincoefs[i], data2->lincoefs[i], eps) )
               return FALSE;
         }

      SCIPexprSortQuadElems(expr1);
      SCIPexprSortQuadElems(expr2);

      /* check if quadratic elements are equal */
      for( i = 0; i < data1->nquadelems; ++i )
         if( data1->quadelems[i].idx1 != data2->quadelems[i].idx1 ||
            data1->quadelems[i].idx2 != data2->quadelems[i].idx2 ||
            !EPSEQ(data1->quadelems[i].coef, data2->quadelems[i].coef, eps) )
            return FALSE;

      /* check if children are equal */
      for( i = 0; i < expr1->nchildren; ++i )
         if( !SCIPexprAreEqual(expr1->children[i], expr2->children[i], eps) )
            return FALSE;

      return TRUE;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* data1;
      SCIP_EXPRDATA_POLYNOMIAL* data2;
      int i;

      if( expr1->nchildren != expr2->nchildren )
         return FALSE;

      data1 = (SCIP_EXPRDATA_POLYNOMIAL*)expr1->data.data;
      data2 = (SCIP_EXPRDATA_POLYNOMIAL*)expr2->data.data;

      if( data1->nmonomials != data2->nmonomials )
         return FALSE;

      if( !EPSEQ(data1->constant, data2->constant, eps) )
         return FALSE;

      /* make sure polynomials are sorted */
      SCIPexprSortMonomials(expr1);
      SCIPexprSortMonomials(expr2);

      /* check if monomials are equal */
      for( i = 0; i < data1->nmonomials; ++i )
      {
         if( !SCIPexprAreMonomialsEqual(data1->monomials[i], data2->monomials[i], eps) )
            return FALSE;
      }

      /* check if children are equal */
      for( i = 0; i < expr1->nchildren; ++i )
      {
         if( !SCIPexprAreEqual(expr1->children[i], expr2->children[i], eps) )
            return FALSE;
      }

      return TRUE;
   }

   case SCIP_EXPR_USER:
   {
      /* @todo could implement this via another user callback */
      return FALSE;
   }

   case SCIP_EXPR_LAST:
      break;
   }

   SCIPerrorMessage("this should never happen\n");
   SCIPABORT();
   return FALSE;  /*lint !e527*/
}

/** aims at simplifying an expression and splitting of a linear expression
 *
 *  If linear variables are split off, expression interpreter data, if stored in the tree, is freed.
 */
SCIP_RETCODE SCIPexprSimplify(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             eps,                /**< threshold, under which positive values are treat as 0 */
   int                   maxexpansionexponent,/**< maximal exponent for which we still expand non-monomial polynomials */
   int                   nvars,              /**< number of variables in expression */
   int*                  nlinvars,           /**< buffer to store number of linear variables in linear part, or NULL if linear part should not be separated */
   int*                  linidxs,            /**< array to store indices of variables in expression tree which belong to linear part, or NULL */
   SCIP_Real*            lincoefs            /**< array to store coefficients of linear part, or NULL */
   )
{
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(eps >= 0.0);

   SCIPdebugMessage("simplify expression: ");
   SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
   SCIPdebugPrintf("\n");

   SCIP_CALL( exprsimplifyConvertToPolynomials(blkmem, expr) );

   SCIPdebugMessage("converted to polynomials: ");
   SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
   SCIPdebugPrintf("\n");

   SCIP_CALL( exprsimplifyFlattenPolynomials(blkmem, messagehdlr, expr, eps, maxexpansionexponent) );

   SCIPdebugMessage("polynomials flattened: ");
   SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
   SCIPdebugPrintf("\n");

   if( nlinvars != NULL )
   {
      /* separate linear part from root polynomial */
      SCIP_CALL( exprsimplifySeparateLinearFromPolynomial(blkmem, expr, eps, nvars, nlinvars, linidxs, lincoefs) );

      SCIPdebugMessage("separated linear part: ");
      SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
      SCIPdebugPrintf("\n");
   }

   SCIP_CALL( exprsimplifyUnconvertPolynomials(blkmem, expr) );

   SCIPdebugMessage("converted back from polynomials: ");
   SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
   SCIPdebugPrintf("\n");

   return SCIP_OKAY;
}

/** evaluates an expression w.r.t. given values for children expressions */
SCIP_RETCODE SCIPexprEvalShallow(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            argvals,            /**< values for children, can be NULL if the expression has no children */
   SCIP_Real*            varvals,            /**< values for variables, can be NULL if the expression operand is not a variable */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression operand is not a parameter */
   SCIP_Real*            val                 /**< buffer to store value */
   )
{
   assert(expr != NULL);
   assert(argvals != NULL || expr->nchildren == 0);

   /* evaluate this expression */
   assert( exprOpTable[expr->op].eval != NULL );
   SCIP_CALL( exprOpTable[expr->op].eval(expr->data, expr->nchildren, argvals, varvals, param, val) );

   return SCIP_OKAY;
}

/** evaluates an expression w.r.t. a point */
SCIP_RETCODE SCIPexprEval(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            varvals,            /**< values for variables, can be NULL if the expression is constant */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_Real*            val                 /**< buffer to store value */
   )
{
   int i;
   SCIP_Real  staticbuf[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_Real* buf;

   /* if many children, get large enough memory to store argument values */
   if( expr->nchildren > SCIP_EXPRESSION_MAXCHILDEST )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&buf, expr->nchildren) );
   }
   else
   {
      buf = staticbuf;
   }

   /* evaluate children */
   for( i = 0; i < expr->nchildren; ++i )
   {
      SCIP_CALL( SCIPexprEval(expr->children[i], varvals, param, &buf[i]) );  /*lint !e644*/
   }

   /* evaluate this expression */
   assert( exprOpTable[expr->op].eval != NULL );
   SCIP_CALL( exprOpTable[expr->op].eval(expr->data, expr->nchildren, buf, varvals, param, val) );

   /* free memory, if allocated before */
   if( staticbuf != buf )
   {
      BMSfreeMemoryArray(&buf);
   }

   return SCIP_OKAY;
}

/** evaluates an expression w.r.t. given interval values for children expressions */
SCIP_RETCODE SCIPexprEvalIntShallow(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        argvals,            /**< interval values for children, can be NULL if the expression has no children */
   SCIP_INTERVAL*        varvals,            /**< interval values for variables, can be NULL if the expression is constant */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_INTERVAL*        val                 /**< buffer to store value */
   )
{
   assert(expr != NULL);
   assert(argvals != NULL || expr->nchildren == 0);

   /* evaluate this expression */
   assert( exprOpTable[expr->op].inteval != NULL );
   SCIP_CALL( exprOpTable[expr->op].inteval(infinity, expr->data, expr->nchildren, argvals, varvals, param, val) );

   return SCIP_OKAY;
}

/** evaluates an expression w.r.t. an interval */
SCIP_RETCODE SCIPexprEvalInt(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values for variables, can be NULL if the expression is constant */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_INTERVAL*        val                 /**< buffer to store value */
   )
{
   int i;
   SCIP_INTERVAL  staticbuf[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_INTERVAL* buf;

   /* if many children, get large enough memory to store argument values */
   if( expr->nchildren > SCIP_EXPRESSION_MAXCHILDEST )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&buf, expr->nchildren) );
   }
   else
   {
      buf = staticbuf;
   }

   /* evaluate children */
   for( i = 0; i < expr->nchildren; ++i )
   {
      SCIP_CALL( SCIPexprEvalInt(expr->children[i], infinity, varvals, param, &buf[i]) );  /*lint !e644*/
   }

   /* evaluate this expression */
   assert( exprOpTable[expr->op].inteval != NULL );
   SCIP_CALL( exprOpTable[expr->op].inteval(infinity, expr->data, expr->nchildren, buf, varvals, param, val) );

   /* free memory, if allocated before */
   if( staticbuf != buf )
   {
      BMSfreeMemoryArray(&buf);
   }

   return SCIP_OKAY;
}

/** evaluates a user expression w.r.t. given values for children expressions */
SCIP_RETCODE SCIPexprEvalUser(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real*            argvals,            /**< values for children */
   SCIP_Real*            val,                /**< buffer to store function value */
   SCIP_Real*            gradient,           /**< buffer to store gradient values, or NULL if not requested */
   SCIP_Real*            hessian             /**< buffer to store values of full Hessian, or NULL if not requested */
   )
{
   SCIP_EXPRDATA_USER* exprdata;

   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_USER);
   assert(argvals != NULL || expr->nchildren == 0);

   exprdata = (SCIP_EXPRDATA_USER*) expr->data.data;
   assert(exprdata->eval != NULL);

   SCIP_CALL( exprdata->eval(exprdata->userdata, expr->nchildren, argvals, val, gradient, hessian) );

   return SCIP_OKAY;
}

/** evaluates a user expression w.r.t. an interval */
SCIP_RETCODE SCIPexprEvalIntUser(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        argvals,            /**< values for children */
   SCIP_INTERVAL*        val,                /**< buffer to store value */
   SCIP_INTERVAL*        gradient,           /**< buffer to store gradient values, or NULL if not requested */
   SCIP_INTERVAL*        hessian             /**< buffer to store values of full Hessian, or NULL if not requested */
   )
{
   SCIP_EXPRDATA_USER* exprdata;

   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_USER);
   assert(argvals != NULL || expr->nchildren == 0);

   exprdata = (SCIP_EXPRDATA_USER*) expr->data.data;

   if( exprdata->inteval == NULL )
   {
      int i;

      for( i = 0; i < expr->nchildren; ++i )
         SCIPintervalSetEntire(infinity, &argvals[i]); /*lint !e613*/
   }
   else
   {
      SCIP_CALL( exprdata->inteval(infinity, exprdata->userdata, expr->nchildren, argvals, val, gradient, hessian) );
   }

   return SCIP_OKAY;
}

/** internal curvature check method */
static
SCIP_RETCODE doCheckCurvature(
   SCIP_EXPR*            expr,               /**< expression to check */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        varbounds,          /**< domains of variables */
   SCIP_INTERVAL*        childbounds,        /**< child bounds buffer array */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_EXPRCURV*        curv,               /**< buffer to store curvature of expression */
   SCIP_EXPRCURV*        childcurv,          /**< buffer array for curvature of children */
   SCIP_INTERVAL*        bounds              /**< buffer to store bounds on expression */
   )
{
   int i;

   assert(childbounds != NULL);
   assert(childcurv != NULL);

   /* check curvature and compute bounds of children
    * constant children can be considered as always linear */
   for( i = 0; i < expr->nchildren; ++i )
   {
      SCIP_CALL( SCIPexprCheckCurvature(expr->children[i], infinity, varbounds, param, &childcurv[i], &childbounds[i]) );  /*lint !e644*/
      if( childbounds[i].inf == childbounds[i].sup )  /*lint !e777*/
         childcurv[i] = SCIP_EXPRCURV_LINEAR;
   }

   /* get curvature and bounds of expr */
   assert(exprOpTable[expr->op].curv != NULL);
   assert(exprOpTable[expr->op].inteval != NULL);

   SCIP_CALL( exprOpTable[expr->op].curv(infinity, expr->data, expr->nchildren, childbounds, childcurv, curv) );
   SCIP_CALL( exprOpTable[expr->op].inteval(infinity, expr->data, expr->nchildren, childbounds, varbounds, param, bounds) );

   return SCIP_OKAY;
}

/** tries to determine the curvature type of an expression w.r.t. given variable domains */
SCIP_RETCODE SCIPexprCheckCurvature(
   SCIP_EXPR*            expr,               /**< expression to check */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        varbounds,          /**< domains of variables */
   SCIP_Real*            param,              /**< values for parameters, can be NULL if the expression is not parameterized */
   SCIP_EXPRCURV*        curv,               /**< buffer to store curvature of expression */
   SCIP_INTERVAL*        bounds              /**< buffer to store bounds on expression */
   )
{
   SCIP_INTERVAL  childboundsbuf[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_INTERVAL* childbounds = NULL;
   SCIP_EXPRCURV  childcurvbuf[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_EXPRCURV* childcurv = NULL;
   SCIP_RETCODE retcode = SCIP_OKAY;

   assert(expr != NULL);
   assert(curv != NULL);
   assert(bounds != NULL);

   /* if many children, get large enough memory to store argument values */
   if( expr->nchildren > SCIP_EXPRESSION_MAXCHILDEST )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&childbounds, expr->nchildren) );
      SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&childcurv, expr->nchildren), TERMINATE );
   }
   else
   {
      childbounds = childboundsbuf;
      childcurv   = childcurvbuf;
   }

   retcode = doCheckCurvature(expr, infinity, varbounds, childbounds, param, curv, childcurv, bounds);

TERMINATE:
   /* free memory, if allocated before */
   if( childboundsbuf != childbounds )
   {
      BMSfreeMemoryArrayNull(&childcurv);
      BMSfreeMemoryArrayNull(&childbounds);
   }

   return retcode;
}

/** under-/overestimates a user expression w.r.t. to given values and bounds for children expressions */
SCIP_RETCODE SCIPexprEstimateUser(
   SCIP_EXPR*           expr,           /**< expression */
   SCIP_Real            infinity,       /**< value to use for infinity */
   SCIP_Real*           argvals,        /**< values for children */
   SCIP_INTERVAL*       argbounds,      /**< bounds for children */
   SCIP_Bool            overestimate,   /**< whether to overestimate the expression */
   SCIP_Real*           coeffs,         /**< buffer to store the linear coefficients for each child expression that gives a valid under-/overestimator */
   SCIP_Real*           constant,       /**< buffer to store the constant value of the linear under-/overestimator */
   SCIP_Bool*           success         /**< buffer to store whether an estimator was successfully computed */
   )
{
   SCIP_EXPRDATA_USER* exprdata;

   assert(expr != NULL);
   assert(expr->op == SCIP_EXPR_USER);
   assert(argvals != NULL || expr->nchildren == 0);
   assert(argbounds != NULL || expr->nchildren == 0);

   exprdata = (SCIP_EXPRDATA_USER*) expr->data.data;

   if( exprdata->estimate != NULL )
   {
      SCIP_CALL( exprdata->estimate(infinity, exprdata->userdata, expr->nchildren, argvals, argbounds, overestimate, coeffs, constant, success ) );
   }
   else
   {
      *success = FALSE;
   }

   return SCIP_OKAY;
}

/** substitutes variables (SCIP_EXPR_VARIDX) by expressions
 *
 *  Note that only the children of the given expr are checked!
 *  A variable with index i is replaced by a copy of substexprs[i], if the latter is not NULL.
 *  If substexprs[i] == NULL, then the variable expression i is not touched.
 */
SCIP_RETCODE SCIPexprSubstituteVars(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPR*            expr,               /**< expression, which of the children may be replaced */
   SCIP_EXPR**           substexprs          /**< array of substitute expressions; single entries can be NULL */
   )
{
   int i;

   assert(blkmem != NULL);
   assert(expr   != NULL);
   assert(substexprs != NULL);

   for( i = 0; i < expr->nchildren; ++i )
   {
      if( expr->children[i]->op == SCIP_EXPR_VARIDX )
      {
         int varidx;
         varidx = expr->children[i]->data.intval;

         assert(varidx >= 0);
         if( substexprs[varidx] != NULL )
         {
            /* replace child i by copy of substexprs[expr->children[i]->opdata.intval] */
            SCIPexprFreeDeep(blkmem, &expr->children[i]);
            SCIP_CALL( SCIPexprCopyDeep(blkmem, &expr->children[i], substexprs[varidx]) );
         }
      }
      else
      {
         /* call recursively */
         SCIP_CALL( SCIPexprSubstituteVars(blkmem, expr->children[i], substexprs) );
      }
   }

   return SCIP_OKAY;
}

/** updates variable indices in expression tree */
void SCIPexprReindexVars(
   SCIP_EXPR*            expr,               /**< expression to update */
   int*                  newindices          /**< new indices of variables */
   )
{
   int i;

   assert(expr != NULL);
   assert(newindices != NULL);

   if( expr->op == SCIP_EXPR_VARIDX )
   {
      expr->data.intval = newindices[expr->data.intval];
      assert(expr->data.intval >= 0);
   }

   for( i = 0; i < expr->nchildren; ++i )
      SCIPexprReindexVars(expr->children[i], newindices);
}

/** updates parameter indices in expression tree */
void SCIPexprReindexParams(
   SCIP_EXPR*            expr,               /**< expression to update */
   int*                  newindices          /**< new indices of variables */
   )
{
   int i;

   assert(expr != NULL);
   assert(newindices != NULL);

   if( expr->op == SCIP_EXPR_PARAM )
   {
      expr->data.intval = newindices[expr->data.intval];
      assert(expr->data.intval >= 0);
   }

   for( i = 0; i < expr->nchildren; ++i )
      SCIPexprReindexParams(expr->children[i], newindices);
}

/** prints an expression */
void SCIPexprPrint(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file for printing, or NULL for stdout */
   const char**          varnames,           /**< names of variables, or NULL for default names */
   const char**          paramnames,         /**< names of parameters, or NULL for default names */
   SCIP_Real*            paramvals           /**< values of parameters, or NULL for not printing */
   )
{
   assert( expr != NULL );

   switch( expr->op )
   {
      /* @Note: 'expr->data.intval' is either between 0 and number of variables-1, if it uses the varnames array, or
       *        between 0 and number of params in the expression tree, if it uses the paramnames array
       *        because, here, we cannot get the values above we cannot assert them
       */
   case SCIP_EXPR_VARIDX:
      if( varnames != NULL )
      {
         assert(varnames[expr->data.intval] != NULL);
         SCIPmessageFPrintInfo(messagehdlr, file, "<%s>", varnames[expr->data.intval]);
      }
      else
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "var%d", expr->data.intval);
      }
      break;

   case SCIP_EXPR_PARAM:
      if( paramnames != NULL )
      {
         assert(paramnames[expr->data.intval] != NULL);
         SCIPmessageFPrintInfo(messagehdlr, file, "%s", paramnames[expr->data.intval]);
      }
      else
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "param%d", expr->data.intval );
      }
      if( paramvals != NULL )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "[%g]", paramvals[expr->data.intval] );
      }
      break;

   case SCIP_EXPR_CONST:
      if (expr->data.dbl < 0.0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "(%g)", expr->data.dbl );
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "%g", expr->data.dbl );
      break;

   case SCIP_EXPR_PLUS:
      SCIPmessageFPrintInfo(messagehdlr, file, "(");
      SCIPexprPrint(expr->children[0], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, " + ");
      SCIPexprPrint(expr->children[1], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, ")");
      break;

   case SCIP_EXPR_MINUS:
      SCIPmessageFPrintInfo(messagehdlr, file, "(");
      SCIPexprPrint(expr->children[0], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, " - ");
      SCIPexprPrint(expr->children[1], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, ")");
      break;

   case SCIP_EXPR_MUL:
      SCIPmessageFPrintInfo(messagehdlr, file, "(");
      SCIPexprPrint(expr->children[0], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, " * ");
      SCIPexprPrint(expr->children[1], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, ")");
      break;

   case SCIP_EXPR_DIV:
      SCIPmessageFPrintInfo(messagehdlr, file, "(");
      SCIPexprPrint(expr->children[0], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, " / ");
      SCIPexprPrint(expr->children[1], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, ")");
      break;

   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_SIGNPOWER:
      SCIPmessageFPrintInfo(messagehdlr, file, "%s(", exprOpTable[expr->op].name);
      SCIPexprPrint(expr->children[0], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, ", %g)", expr->data.dbl);
      break;

   case SCIP_EXPR_INTPOWER:
      SCIPmessageFPrintInfo(messagehdlr, file, "power(");
      SCIPexprPrint(expr->children[0], messagehdlr, file, varnames, paramnames, paramvals);
      SCIPmessageFPrintInfo(messagehdlr, file, ", %d)", expr->data.intval);
      break;

   case SCIP_EXPR_SQUARE:
   case SCIP_EXPR_SQRT:
   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   {
      int i;

      SCIPmessageFPrintInfo(messagehdlr, file, "%s(", exprOpTable[expr->op].name);

      for( i = 0; i < expr->nchildren; ++i )
      {
         SCIPexprPrint(expr->children[i], messagehdlr, file, varnames, paramnames, paramvals);
         if( i + 1 < expr->nchildren )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, ", ");
         }
      }

      SCIPmessageFPrintInfo(messagehdlr, file, ")");
      break;
   }

   case SCIP_EXPR_SUM:
   case SCIP_EXPR_PRODUCT:
   {
      switch( expr->nchildren )
      {
      case 0:
         SCIPmessageFPrintInfo(messagehdlr, file, expr->op == SCIP_EXPR_SUM ? "0" : "1");
         break;
      case 1:
         SCIPexprPrint(expr->children[0], messagehdlr, file, varnames, paramnames, paramvals);
         break;
      default:
      {
         int i;
         char opstr[SCIP_MAXSTRLEN];

         SCIPmessageFPrintInfo(messagehdlr, file, "(");
         for( i = 0; i < expr->nchildren; ++i )
         {
            if( i > 0 )
            {
               (void) SCIPsnprintf(opstr, SCIP_MAXSTRLEN, "%s", expr->op == SCIP_EXPR_SUM ? " + " : " * ");
               SCIPmessageFPrintInfo(messagehdlr, file, opstr);
            }
            SCIPexprPrint(expr->children[i], messagehdlr, file, varnames, paramnames, paramvals);
         }
         SCIPmessageFPrintInfo(messagehdlr, file, ")");
      }
      }
      break;
   }

   case SCIP_EXPR_LINEAR:
   {
      SCIP_Real constant;
      int i;

      constant = ((SCIP_Real*)expr->data.data)[expr->nchildren];

      if( expr->nchildren == 0 )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%.20g", constant);
         break;
      }

      SCIPmessageFPrintInfo(messagehdlr, file, "(");

      if( constant != 0.0 )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%.20g", constant);
      }

      for( i = 0; i < expr->nchildren; ++i )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, " %+.20g ", ((SCIP_Real*)expr->data.data)[i]);
         SCIPexprPrint(expr->children[i], messagehdlr, file, varnames, paramnames, paramvals);
      }

      SCIPmessageFPrintInfo(messagehdlr, file, ")");
      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* quadraticdata;
      int i;

      quadraticdata = (SCIP_EXPRDATA_QUADRATIC*)expr->data.data;
      assert(quadraticdata != NULL);

      SCIPmessageFPrintInfo(messagehdlr, file, "(");

      if( quadraticdata->constant != 0.0 )
         SCIPmessageFPrintInfo(messagehdlr, file, " %+.20g ", quadraticdata->constant);

      if( quadraticdata->lincoefs != NULL )
         for( i = 0; i < expr->nchildren; ++i )
         {
            if( quadraticdata->lincoefs[i] == 0.0 )
               continue;
            SCIPmessageFPrintInfo(messagehdlr, file, " %+.20g ", quadraticdata->lincoefs[i]);
            SCIPexprPrint(expr->children[i], messagehdlr, file, varnames, paramnames, paramvals);
         }

      for( i = 0; i < quadraticdata->nquadelems; ++i )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, " %+.20g ", quadraticdata->quadelems[i].coef);
         SCIPexprPrint(expr->children[quadraticdata->quadelems[i].idx1], messagehdlr, file, varnames, paramnames, paramvals);
         if( quadraticdata->quadelems[i].idx1 == quadraticdata->quadelems[i].idx2 )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "^2");
         }
         else
         {
            SCIPmessageFPrintInfo(messagehdlr, file, " * ");
            SCIPexprPrint(expr->children[quadraticdata->quadelems[i].idx2], messagehdlr, file, varnames, paramnames, paramvals);
         }
      }

      SCIPmessageFPrintInfo(messagehdlr, file, ")");
      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL*   monomialdata;
      int i;
      int j;

      SCIPmessageFPrintInfo(messagehdlr, file, "(");

      polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)expr->data.data;
      assert(polynomialdata != NULL);

      if( polynomialdata->constant != 0.0 || polynomialdata->nmonomials == 0 )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%.20g", polynomialdata->constant);
      }

      for( i = 0; i < polynomialdata->nmonomials; ++i )
      {
         monomialdata = polynomialdata->monomials[i];
         SCIPmessageFPrintInfo(messagehdlr, file, " %+.20g", monomialdata->coef);

         for( j = 0; j < monomialdata->nfactors; ++j )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, " * ");

            SCIPexprPrint(expr->children[monomialdata->childidxs[j]], messagehdlr, file, varnames, paramnames, paramvals);
            if( monomialdata->exponents[j] < 0.0 )
            {
               SCIPmessageFPrintInfo(messagehdlr, file, "^(%.20g)", monomialdata->exponents[j]);
            }
            else if( monomialdata->exponents[j] != 1.0 )
            {
               SCIPmessageFPrintInfo(messagehdlr, file, "^%.20g", monomialdata->exponents[j]);
            }
         }
      }

      SCIPmessageFPrintInfo(messagehdlr, file, ")");
      break;
   }

   case SCIP_EXPR_USER:
   {
      SCIP_EXPRDATA_USER* exprdata;
      int i;

      exprdata = (SCIP_EXPRDATA_USER*)expr->data.data;
      assert(exprdata != NULL);

      if( exprdata->print != NULL )
      {
         exprdata->print(exprdata->userdata, messagehdlr, file);
      }
      else
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "user");
      }

      SCIPmessageFPrintInfo(messagehdlr, file, "(");
      for( i = 0; i < expr->nchildren; ++i )
      {
         if( i > 0 )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, ",");
         }
         SCIPexprPrint(expr->children[i], messagehdlr, file, varnames, paramnames, paramvals);
      }
      SCIPmessageFPrintInfo(messagehdlr, file, ")");

      break;
   }

   case SCIP_EXPR_LAST:
   {
      SCIPerrorMessage("invalid expression\n");
      SCIPABORT();
   }
   }
}

/** parses an expression from a string */
SCIP_RETCODE SCIPexprParse(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to created expression */
   const char*           str,                /**< pointer to the string to be parsed */
   const char*           lastchar,           /**< pointer to the last char of str that should be parsed */
   int*                  nvars,              /**< buffer to store number of variables */
   int*                  varnames,           /**< buffer to store variable names, prefixed by index (as int) */
   int                   varnameslength      /**< length of the varnames buffer array */
   )
{
   SCIP_HASHTABLE* vartable;
   SCIP_RETCODE retcode;

   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(str != NULL);
   assert(lastchar != NULL);
   assert(nvars != NULL);
   assert(varnames != NULL);

   *nvars = 0;

   /* create a hash table for variable names and corresponding expression index
    * for each variable, we store its name, prefixed with the assigned index in the first sizeof(int) bytes
    */
   SCIP_CALL( SCIPhashtableCreate(&vartable, blkmem, 10, exprparseVarTableGetKey, SCIPhashKeyEqString,
         SCIPhashKeyValString, NULL) );

   retcode = exprParse(blkmem, messagehdlr, expr, str, (int) (lastchar - str + 1), lastchar, nvars, &varnames,
      &varnameslength, vartable, 0);

   SCIPhashtableFree(&vartable);

   return retcode;
}


/**@} */

/**@name Expression tree methods */
/**@{ */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPexprtreeGetRoot
#undef SCIPexprtreeGetNVars
#undef SCIPexprtreeGetNParams
#undef SCIPexprtreeGetParamVals
#undef SCIPexprtreeSetParamVal
#undef SCIPexprtreeGetInterpreterData
#undef SCIPexprtreeSetInterpreterData
#undef SCIPexprtreeFreeInterpreterData
#undef SCIPexprtreeHasParam
#undef SCIPexprtreeGetMaxDegree
#undef SCIPexprtreeEval
#undef SCIPexprtreeEvalInt
#undef SCIPexprtreePrint

/** returns root expression of an expression tree */
SCIP_EXPR* SCIPexprtreeGetRoot(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree != NULL);

   return tree->root;
}

/** returns number of variables in expression tree */
int SCIPexprtreeGetNVars(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree != NULL);

   return tree->nvars;
}

/** returns number of parameters in expression tree */
int SCIPexprtreeGetNParams(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree != NULL);

   return tree->nparams;
}

/** returns values of parameters or NULL if none */
SCIP_Real* SCIPexprtreeGetParamVals(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree != NULL);

   return tree->params;
}

/** sets value of a single parameter in expression tree */
void SCIPexprtreeSetParamVal(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   paramidx,           /**< index of parameter */
   SCIP_Real             paramval            /**< new value of parameter */
   )
{
   assert(tree != NULL);
   assert(paramidx >= 0);
   assert(paramidx < tree->nparams);
   assert(tree->params != NULL);

   tree->params[paramidx] = paramval;
}

/** gets data of expression tree interpreter, or NULL if not set */
SCIP_EXPRINTDATA* SCIPexprtreeGetInterpreterData(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree != NULL);

   return tree->interpreterdata;
}

/** sets data of expression tree interpreter */
void SCIPexprtreeSetInterpreterData(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_EXPRINTDATA*     interpreterdata     /**< expression interpreter data */
   )
{
   assert(tree != NULL);
   assert(interpreterdata != NULL);
   assert(tree->interpreterdata == NULL);

   tree->interpreterdata = interpreterdata;
}

/** frees data of expression tree interpreter, if any */
SCIP_RETCODE SCIPexprtreeFreeInterpreterData(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   if( tree->interpreterdata != NULL )
   {
      SCIP_CALL( SCIPexprintFreeData(&tree->interpreterdata) );
      assert(tree->interpreterdata == NULL);
   }

   return SCIP_OKAY;
}

/** indicates whether there are parameterized constants (SCIP_EXPR_PARAM) in expression tree */
SCIP_Bool SCIPexprtreeHasParam(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree != NULL);

   return SCIPexprHasParam(tree->root);
}

/** Gives maximal degree of expression in expression tree.
 *
 *  If constant expression, gives 0,
 *  if linear expression, gives 1,
 *  if polynomial expression, gives its maximal degree,
 *  otherwise (nonpolynomial nonconstant expressions) gives at least SCIP_EXPR_DEGREEINFINITY.
 */
SCIP_RETCODE SCIPexprtreeGetMaxDegree(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int*                  maxdegree           /**< buffer to store maximal degree */
   )
{
   assert(tree != NULL);

   SCIP_CALL( SCIPexprGetMaxDegree(tree->root, maxdegree) );

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. a point */
SCIP_RETCODE SCIPexprtreeEval(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values for variables */
   SCIP_Real*            val                 /**< buffer to store expression tree value */
   )
{
   assert(tree    != NULL);
   assert(varvals != NULL || tree->nvars == 0);
   assert(val     != NULL);

   SCIP_CALL( SCIPexprEval(tree->root, varvals, tree->params, val) );

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. an interval */
SCIP_RETCODE SCIPexprtreeEvalInt(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< intervals for variables */
   SCIP_INTERVAL*        val                 /**< buffer to store expression tree value */
   )
{
   assert(tree    != NULL);
   assert(varvals != NULL || tree->nvars == 0);
   assert(val     != NULL);

   SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, varvals, tree->params, val) );

   return SCIP_OKAY;
}

/** prints an expression tree */
void SCIPexprtreePrint(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file for printing, or NULL for stdout */
   const char**          varnames,           /**< names of variables, or NULL for default names */
   const char**          paramnames          /**< names of parameters, or NULL for default names */
   )
{
   assert(tree != NULL);

   SCIPexprPrint(tree->root, messagehdlr, file, varnames, paramnames, tree->params);
}


/** creates an expression tree */
SCIP_RETCODE SCIPexprtreeCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRTREE**       tree,               /**< buffer to store address of created expression tree */
   SCIP_EXPR*            root,               /**< pointer to root expression, not copied deep !, can be NULL */
   int                   nvars,              /**< number of variables in variable mapping */
   int                   nparams,            /**< number of parameters in expression */
   SCIP_Real*            params              /**< values for parameters, or NULL (if NULL but nparams > 0, then params is initialized with zeros) */
   )
{
   assert(blkmem != NULL);
   assert(tree   != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, tree) );

   (*tree)->blkmem    = blkmem;
   (*tree)->root      = root;
   (*tree)->nvars     = nvars;
   (*tree)->vars      = NULL;
   (*tree)->nparams   = nparams;
   (*tree)->interpreterdata = NULL;

   if( params != NULL )
   {
      assert(nparams > 0);
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*tree)->params, params, nparams) );
   }
   else if( nparams > 0 )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*tree)->params, nparams) );
      BMSclearMemoryArray((*tree)->params, nparams);
   }
   else
   {
      assert(nparams == 0);
      (*tree)->params = NULL;
   }

   return SCIP_OKAY;
}

/** copies an expression tree */
SCIP_RETCODE SCIPexprtreeCopy(
   BMS_BLKMEM*           blkmem,             /**< block memory that should be used in new expression tree */
   SCIP_EXPRTREE**       targettree,         /**< buffer to store address of copied expression tree */
   SCIP_EXPRTREE*        sourcetree          /**< expression tree to copy */
   )
{
   assert(blkmem     != NULL);
   assert(targettree != NULL);
   assert(sourcetree != NULL);

   /* copy expression tree "header" */
   SCIP_ALLOC( BMSduplicateBlockMemory(blkmem, targettree, sourcetree) );

   /* we may have a new block memory; and we do not want to keep the others interpreter data */
   (*targettree)->blkmem          = blkmem;
   (*targettree)->interpreterdata = NULL;

   /* copy variables, if any */
   if( sourcetree->vars != NULL )
   {
      assert(sourcetree->nvars > 0);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*targettree)->vars, sourcetree->vars, sourcetree->nvars) );
   }

   /* copy parameters, if any */
   if( sourcetree->params != NULL )
   {
      assert(sourcetree->nparams > 0);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*targettree)->params, sourcetree->params, sourcetree->nparams) );
   }

   /* copy expression */
   SCIP_CALL( SCIPexprCopyDeep(blkmem, &(*targettree)->root, sourcetree->root) );

   return SCIP_OKAY;
}

/** frees an expression tree */
SCIP_RETCODE SCIPexprtreeFree(
   SCIP_EXPRTREE**       tree                /**< pointer to expression tree that is freed */
   )
{
   assert( tree != NULL);
   assert(*tree != NULL);

   SCIP_CALL( SCIPexprtreeFreeInterpreterData(*tree) );

   if( (*tree)->root != NULL )
   {
      SCIPexprFreeDeep((*tree)->blkmem, &(*tree)->root);
      assert((*tree)->root == NULL);
   }

   BMSfreeBlockMemoryArrayNull((*tree)->blkmem, &(*tree)->vars,   (*tree)->nvars  );
   BMSfreeBlockMemoryArrayNull((*tree)->blkmem, &(*tree)->params, (*tree)->nparams);

   BMSfreeBlockMemory((*tree)->blkmem, tree);

   return SCIP_OKAY;
}

/** sets number and values of all parameters in expression tree */
SCIP_RETCODE SCIPexprtreeSetParams(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nparams,            /**< number of parameters */
   SCIP_Real*            paramvals           /**< values of parameters, can be NULL if nparams == 0 */
   )
{
   assert(tree != NULL);
   assert(paramvals != NULL || nparams == 0);

   if( nparams == 0 )
   {
      BMSfreeBlockMemoryArrayNull(tree->blkmem, &tree->params, tree->nparams);
   }
   else if( tree->params != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(tree->blkmem, &tree->params, tree->nparams, nparams) );
      BMScopyMemoryArray(tree->params, paramvals, nparams);
   }
   else
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(tree->blkmem, &tree->params, paramvals, nparams) );
   }

   tree->nparams = nparams;
   assert(tree->params != NULL || tree->nparams == 0);

   return SCIP_OKAY;
}


/** gives the number of usages for each variable in the expression tree */
void SCIPexprtreeGetVarsUsage(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int*                  varsusage           /**< array where to store for each variable how often it is used in the tree */
   )
{
   assert(tree != NULL);
   assert(varsusage != NULL);

   if( tree->nvars == 0 )
      return;

   BMSclearMemoryArray(varsusage, tree->nvars);
   SCIPexprGetVarsUsage(tree->root, varsusage);
}

/** aims at simplifying an expression and splitting of a linear expression
 *
 *  If linear variables are split off, expression interpreter data, if stored in the tree, is freed.
 */
SCIP_RETCODE SCIPexprtreeSimplify(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Real             eps,                /**< threshold, under which positive values are treat as 0 */
   int                   maxexpansionexponent,/**< maximal exponent for which we still expand non-monomial polynomials */
   int*                  nlinvars,           /**< buffer to store number of linear variables in linear part, or NULL if linear part should not be separated */
   int*                  linidxs,            /**< array to store indices of variables in expression tree which belong to linear part, or NULL */
   SCIP_Real*            lincoefs            /**< array to store coefficients of linear part, or NULL */
   )
{
#ifndef NDEBUG
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real* testx;
   SCIP_Real testval_before;
   SCIP_Real testval_after;
   int i;
#endif

   assert(tree != NULL);

#ifndef NDEBUG
   SCIP_CALL( SCIPrandomCreate(&randnumgen, tree->blkmem, 42) );

   SCIP_ALLOC( BMSallocMemoryArray(&testx, SCIPexprtreeGetNVars(tree)) );  /*lint !e666*/
   for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
      testx[i] = SCIPrandomGetReal(randnumgen, -100.0, 100.0);  /*lint !e644*/
   SCIP_CALL( SCIPexprtreeEval(tree, testx, &testval_before) );

   SCIPrandomFree(&randnumgen, tree->blkmem);
#endif

   /* we should be careful about declaring numbers close to zero as zero, so take eps^2 as tolerance */
   SCIP_CALL( SCIPexprSimplify(tree->blkmem, messagehdlr, tree->root, eps*eps, maxexpansionexponent, tree->nvars, nlinvars, linidxs, lincoefs) );

#ifndef NDEBUG
   SCIP_CALL( SCIPexprtreeEval(tree, testx, &testval_after) );
   if( nlinvars != NULL && testval_before == testval_before )  /*lint !e777*/
      for( i = 0; i < *nlinvars; ++i )
         testval_after += lincoefs[i] * testx[linidxs[i]];
   assert(testval_before != testval_before || testval_before == testval_after || EPSZ(SCIPrelDiff(testval_before, testval_after), eps));  /*lint !e777*/
   BMSfreeMemoryArray(&testx);
#endif

   /* removing something from the the tree may invalidate the interpreter data */
   if( nlinvars != NULL && *nlinvars > 0 )
      SCIP_CALL( SCIPexprtreeFreeInterpreterData(tree) );

   return SCIP_OKAY;
}

/** adds an expression to the root expression of the tree
 *
 *  The root is replaced with an SCIP_EXPR_PLUS expression which has the previous root and the given expression (or a copy of it) as children.
 */
SCIP_RETCODE SCIPexprtreeAddExpr(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_EXPR*            expr,               /**< expression to add to tree */
   SCIP_Bool             copyexpr            /**< whether expression should be copied */
   )
{
   assert(tree != NULL);
   assert(tree->root != NULL);

   /* adding something to the tree may invalidate the interpreter data */
   SCIP_CALL( SCIPexprtreeFreeInterpreterData(tree) );

   if( copyexpr )
   {
      SCIP_CALL( SCIPexprCopyDeep(tree->blkmem, &expr, expr) );
   }

   SCIP_CALL( SCIPexprCreate(tree->blkmem, &tree->root, SCIP_EXPR_PLUS, tree->root, expr) );

   return SCIP_OKAY;
}

/** tries to determine the curvature type of an expression tree w.r.t. given variable domains */
SCIP_RETCODE SCIPexprtreeCheckCurvature(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varbounds,          /**< domains of variables */
   SCIP_EXPRCURV*        curv,               /**< buffer to store curvature of expression */
   SCIP_INTERVAL*        bounds              /**< buffer to store bounds on expression, or NULL if not needed */
   )
{
   SCIP_INTERVAL exprbounds;

   assert(tree != NULL);
   assert(tree->root != NULL);

   SCIP_CALL( SCIPexprCheckCurvature(tree->root, infinity, varbounds, tree->params, curv, &exprbounds) );

   if( bounds != NULL )
      *bounds = exprbounds;

   return SCIP_OKAY;
}

/** substitutes variables (SCIP_EXPR_VARIDX) in an expression tree by expressions
 *
 *  A variable with index i is replaced by a copy of substexprs[i], if that latter is not NULL.
 *  If substexprs[i] == NULL, then the variable expression i is not touched.
 */
SCIP_RETCODE SCIPexprtreeSubstituteVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_EXPR**           substexprs          /**< array of substitute expressions; single entries can be NULL */
   )
{
   assert(tree != NULL);
   assert(tree->root != NULL);

   if( tree->root->op == SCIP_EXPR_VARIDX )
   {
      int varidx;

      varidx = tree->root->data.intval;
      assert(varidx >= 0);
      if( substexprs[varidx] != NULL )
      {
         /* substitute root expression */
         SCIPexprFreeDeep(tree->blkmem, &tree->root);
         SCIP_CALL( SCIPexprCopyDeep(tree->blkmem, &tree->root, substexprs[varidx]) );
      }
   }
   else
   {
      /* check children (and grandchildren and so on...) of root expression */
      SCIP_CALL( SCIPexprSubstituteVars(tree->blkmem, tree->root, substexprs) );
   }

   /* substitution of variables should invalidate interpreter data */
   SCIP_CALL( SCIPexprtreeFreeInterpreterData(tree) );

   return SCIP_OKAY;
}

/**@} */

/**@name Quadratic element methods */
/**@{ */

/** comparing two quadratic elements
 *
 *  a is better than b if index1 of a is smaller than index1 of b or index1 of both is equal but index2 of a is smaller than index2 of b
 */
#define QUADELEMS_ISBETTER(a, b) ( ((a).idx1 < (b).idx1) || ((a).idx1 == (b).idx1 && (a).idx2 < (b).idx2) )

/** swaps two quadratic elements */
#define QUADELEMS_SWAP(x,y)                     \
   {                                            \
      SCIP_QUADELEM temp = x;                   \
      x = y;                                    \
      y = temp;                                 \
   }

/** quicksort an array of quadratic elements; pivot is the medial element (taken from scip/sorttpl.c) */
static
void quadelemsQuickSort(
   SCIP_QUADELEM*        elems,              /**< array to be sorted */
   int                   start,              /**< starting index */
   int                   end                 /**< ending index */
   )
{
   assert(start <= end);

   /* use quick sort for long lists */
   while( end - start >= 25 ) /* 25 was SORTTPL_SHELLSORTMAX in sorttpl.c */
   {
      SCIP_QUADELEM pivotkey;
      int lo;
      int hi;
      int mid;

      /* select pivot element */
      mid = (start+end)/2;
      pivotkey = elems[mid];

      /* partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end] */
      lo = start;
      hi = end;
      for( ;; )
      {
         while( lo < end   &&  QUADELEMS_ISBETTER(elems[lo], pivotkey) )
            lo++;
         while( hi > start && !QUADELEMS_ISBETTER(elems[hi], pivotkey) )
            hi--;

         if( lo >= hi )
            break;

         QUADELEMS_SWAP(elems[lo], elems[hi]);

         lo++;
         hi--;
      }
      assert(hi == lo-1 || hi == start);

      /* skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)*/
      while( lo < end && !QUADELEMS_ISBETTER(pivotkey, elems[lo]) )
         lo++;

      /* make sure that we have at least one element in the smaller partition */
      if( lo == start )
      {
         /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
         assert(!QUADELEMS_ISBETTER(elems[mid], pivotkey)); /* the pivot element did not change its position */
         assert(!QUADELEMS_ISBETTER(pivotkey, elems[mid]));
         QUADELEMS_SWAP(elems[lo], elems[mid]);
         lo++;
      }

      /* sort the smaller partition by a recursive call, sort the larger part without recursion */
      if( hi - start <= end - lo )
      {
         /* sort [start,hi] with a recursive call */
         if( start < hi )
            quadelemsQuickSort(elems, start, hi);

         /* now focus on the larger part [lo,end] */
         start = lo;
      }
      else
      {
         /* sort [lo,end] with a recursive call */
         if( lo < end )
            quadelemsQuickSort(elems, lo, end);

         /* now focus on the larger part [start,hi] */
         end = hi;
      }
   }

   /* use shell sort on the remaining small list */
   if( end - start >= 1 )
   {
      static const int incs[3] = {1, 5, 19}; /* sequence of increments */
      int k;

      for( k = 2; k >= 0; --k )
      {
         int h;
         int i;

         for( h = incs[k], i = h + start; i <= end; ++i )
         {
            int j;
            SCIP_QUADELEM tempkey = elems[i];

            j = i;
            while( j >= h && QUADELEMS_ISBETTER(tempkey, elems[j-h]) )
            {
               elems[j] = elems[j-h];
               j -= h;
            }

            elems[j] = tempkey;
         }
      }
   }
}

/** sorts an array of quadratic elements
 *
 *  The elements are sorted such that the first index is increasing and
 *  such that among elements with the same first index, the second index is increasing.
 *  For elements with same first and second index, the order is not defined.
 */
void SCIPquadelemSort(
   SCIP_QUADELEM*        quadelems,          /**< array of quadratic elements */
   int                   nquadelems          /**< number of quadratic elements */
   )
{
   if( nquadelems == 0 )
      return;

#ifndef NDEBUG
   {
      int i;
      for( i = 0; i < nquadelems; ++i )
         assert(quadelems[i].idx1 <= quadelems[i].idx2);
   }
#endif

   quadelemsQuickSort(quadelems, 0, nquadelems-1);
}

/** Finds an index pair in a sorted array of quadratic elements.
 *
 *  If (idx1,idx2) is found in quadelems, then returns TRUE and stores position of quadratic element in *pos.
 *  If (idx1,idx2) is not found in quadelems, then returns FALSE and stores position where a quadratic element with these indices would be inserted in *pos.
 *  Assumes that idx1 <= idx2.
 */
SCIP_Bool SCIPquadelemSortedFind(
   SCIP_QUADELEM*        quadelems,          /**< array of quadratic elements */
   int                   idx1,               /**< index of first  variable in element to search for */
   int                   idx2,               /**< index of second variable in element to search for */
   int                   nquadelems,         /**< number of quadratic elements in array */
   int*                  pos                 /**< buffer to store position of found quadratic element or position where it would be inserted, or NULL */
   )
{
   int left;
   int right;

   assert(quadelems != NULL || nquadelems == 0);
   assert(idx1 <= idx2);

   if( nquadelems == 0 )
   {
      if( pos != NULL )
         *pos = 0;
      return FALSE;
   }

   left = 0;
   right = nquadelems - 1;
   while( left <= right )
   {
      int middle;

      middle = (left+right)/2;
      assert(0 <= middle && middle < nquadelems);

      if( idx1 < quadelems[middle].idx1 || (idx1 == quadelems[middle].idx1 && idx2 < quadelems[middle].idx2) )  /*lint !e613*/
         right = middle - 1;
      else if( quadelems[middle].idx1 < idx1 || (quadelems[middle].idx1 == idx1 && quadelems[middle].idx2 < idx2) )  /*lint !e613*/
         left  = middle + 1;
      else
      {
         if( pos != NULL )
            *pos = middle;
         return TRUE;
      }
   }
   assert(left == right+1);

   if( pos != NULL )
      *pos = left;
   return FALSE;
}

/** Adds quadratic elements with same index and removes elements with coefficient 0.0.
 *
 *  Assumes that elements have been sorted before.
 */
void SCIPquadelemSqueeze(
   SCIP_QUADELEM*        quadelems,          /**< array of quadratic elements */
   int                   nquadelems,         /**< number of quadratic elements */
   int*                  nquadelemsnew       /**< pointer to store new (reduced) number of quadratic elements */
   )
{
   int i;
   int next;

   assert(quadelems     != NULL);
   assert(nquadelemsnew != NULL);
   assert(nquadelems    >= 0);

   i = 0;
   next = 0;
   while( next < nquadelems )
   {
      /* assert that array is sorted */
      assert(QUADELEMS_ISBETTER(quadelems[i], quadelems[next]) ||
         (quadelems[i].idx1 == quadelems[next].idx1 && quadelems[i].idx2 == quadelems[next].idx2));

      /* skip elements with coefficient 0.0 */
      if( quadelems[next].coef == 0.0 )
      {
         ++next;
         continue;
      }

      /* if next element has same index as previous one, add it to the previous one */
      if( i >= 1 &&
         quadelems[i-1].idx1 == quadelems[next].idx1 &&
         quadelems[i-1].idx2 == quadelems[next].idx2 )
      {
         quadelems[i-1].coef += quadelems[next].coef;
         ++next;
         continue;
      }

      /* otherwise, move next element to current position */
      quadelems[i] = quadelems[next];
      ++i;
      ++next;
   }
   assert(next == nquadelems);

   /* now i should point to the position after the last valid element, i.e., it is the remaining number of elements */
   *nquadelemsnew = i;
}

/**@} */

/**@name Expression graph node private methods */
/**@{ */

/** adds a parent to an expression graph node */
static
SCIP_RETCODE exprgraphNodeAddParent(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node where to add a parent */
   SCIP_EXPRGRAPHNODE*   parent              /**< parent node */
   )
{
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(node->depth >= 0);
   assert(node->pos >= 0);
   assert(parent != NULL);
   assert(parent->depth >= 0);
   assert(parent->pos >= 0);
   assert(parent->depth > node->depth); /* a parent node need to have larger depth */

   ensureBlockMemoryArraySize(blkmem, &node->parents, &node->parentssize, node->nparents + 1);
   assert(node->nparents < node->parentssize);

   node->parents[node->nparents] = parent;
   ++node->nparents;

   /* update sorted flag */
   node->parentssorted = (node->nparents <= 1) || (node->parentssorted && (exprgraphnodecomp((void*)node->parents[node->nparents-2], (void*)parent) <= 0));

   return SCIP_OKAY;
}

/** ensures that array of parents in a node is sorted */
static
void exprgraphNodeSortParents(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   if( node->parentssorted )
   {
#ifndef NDEBUG
      int i;
      for( i = 1; i < node->nparents; ++i )
         assert(exprgraphnodecomp((void*)node->parents[i-1], (void*)node->parents[i]) <= 0);
#endif
      return;
   }

   SCIPsortPtr((void**)node->parents, exprgraphnodecomp, node->nparents);

   node->parentssorted = TRUE;
}

/** removes a parent from an expression graph node
 *
 *  If the node is not used and has no other parents, then it is freed.
 */
static
SCIP_RETCODE exprgraphNodeRemoveParent(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  node,               /**< expression graph node where to remove a parent, *node will be set to NULL */
   SCIP_EXPRGRAPHNODE*   parent              /**< parent node to remove */
   )
{
   SCIP_EXPRGRAPHNODE* node_;
   int pos;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(*node != NULL);
   assert((*node)->depth >= 0);
   assert((*node)->pos >= 0);
   assert((*node)->nparents > 0);
   assert(parent != NULL);
   assert(parent->depth >= 0);
   assert(parent->pos >= 0);
   assert(parent->depth > (*node)->depth); /* a parent node need to have larger depth */

   /* find parent */
   exprgraphNodeSortParents(*node);
   (void) SCIPsortedvecFindPtr((void**)(*node)->parents, exprgraphnodecomp, (void*)parent, (*node)->nparents, &pos);
   assert(pos >= 0);
   assert(pos < (*node)->nparents);
   assert((*node)->parents[pos] == parent);

   /* move last parent to pos, if pos is before last
    * update sorted flag */
   if( pos < (*node)->nparents-1 )
   {
      (*node)->parents[pos]  = (*node)->parents[(*node)->nparents-1];
      (*node)->parentssorted = ((*node)->nparents <= 2);
   }
   --(*node)->nparents;

   /* keep pointer to *node in case it is still used */
   node_ = (*node)->nuses > 0 ? *node : NULL;

   /* capture and release node so it is freed if possible */
   SCIPexprgraphCaptureNode(*node);
   SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );

   /* restore pointer, if node still exists */
   *node = node_;

   return SCIP_OKAY;
}

/** checks if a node is parent of a node */
static
SCIP_Bool exprgraphNodeIsParent(
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_EXPRGRAPHNODE*   parent              /**< parent to look for */
   )
{
   int pos;

   assert(node   != NULL);
   assert(parent != NULL);

   /* if depth of node is at least as high as depth of parent, parent cannot be parent of node */
   if( node->depth >= parent->depth || node->nparents == 0 )
      return FALSE;
   assert(node->parents != NULL);

   /* ensure parents array is sorted */
   exprgraphNodeSortParents(node);

   return SCIPsortedvecFindPtr((void**)node->parents, exprgraphnodecomp, (void*)parent, node->nparents, &pos);
}

/** adds expression graph nodes to the array of children of a sum, product, linear, quadratic, or polynomial expression
 *
 *  For a sum or product expression, this corresponds to add additional summands and factors, resp.
 *  For a linear expression, this corresponds to add each expression with coefficient 1.0.
 *  For a quadratic or polynomial expression, only the children array may be enlarged, the expression itself remains the same.
 *
 *  It is assumed that node and all exprs are in the expression graph already.
 *  It is assumed that all expressions that are added have lower depth than node.
 */
static
SCIP_RETCODE exprgraphNodeAddChildren(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   int                   nexprs,             /**< number of children to add */
   SCIP_EXPRGRAPHNODE**  exprs,              /**< children nodes to add */
   int*                  childmap            /**< array where to store mapping of indices from exprs to children array in node, or NULL if not of interest */
   )
{
   int i;
   int j;
   int orignchildren;
   SCIP_Bool existsalready;

   assert(blkmem != NULL);
   assert(node != NULL);
   assert(node->depth > 0);
   assert(node->pos >= 0);
   assert(node->op == SCIP_EXPR_SUM || node->op == SCIP_EXPR_PRODUCT || node->op == SCIP_EXPR_LINEAR || node->op == SCIP_EXPR_QUADRATIC || node->op == SCIP_EXPR_POLYNOMIAL);
   assert(exprs != NULL || nexprs == 0);

   if( nexprs == 0 )
      return SCIP_OKAY;

   orignchildren = node->nchildren;
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &node->children, node->nchildren, node->nchildren + nexprs) );

   for( i = 0; i < nexprs; ++i )
   {
      assert(exprs[i]->depth >= 0);           /*lint !e613*/
      assert(exprs[i]->pos >= 0);             /*lint !e613*/
      assert(exprs[i]->depth < node->depth);  /*lint !e613*/

      /* check if exprs[i] is a child already, if not SUM or PRODUCT */
      existsalready = FALSE;
      if( node->op != SCIP_EXPR_SUM && node->op != SCIP_EXPR_PRODUCT )
         for( j = 0; j < orignchildren; ++j )
            /* during simplification of polynomials, their may be NULL's in children array */
            if( node->children[j] != NULL && node->children[j] == exprs[i] )  /*lint !e613*/
            {
               existsalready = TRUE;
               break;
            }

      if( !existsalready )
      {
         /* add exprs[i] to children array */
         node->children[node->nchildren] = exprs[i];  /*lint !e613*/
         SCIP_CALL( exprgraphNodeAddParent(blkmem, exprs[i], node) );  /*lint !e613*/
         if( childmap != NULL )
            childmap[i] = node->nchildren;
         ++node->nchildren;
      }
      else
      {
         if( childmap != NULL )
            childmap[i] = j;  /*lint !e644*/
         if( node->op == SCIP_EXPR_LINEAR )
         {
            /* if linear expression, increase coefficient by 1.0 */
            ((SCIP_Real*)node->data.data)[j] += 1.0;
         }
      }
   }

   /* shrink children array to actually used size */
   SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &node->children, orignchildren + nexprs, node->nchildren) );

   if( node->op == SCIP_EXPR_LINEAR && node->nchildren > orignchildren )
   {
      /* if linear expression, then add 1.0 coefficients for new expressions */
      SCIP_Real* data;

      data = (SCIP_Real*)node->data.data;
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &data, orignchildren + 1, node->nchildren + 1) );
      data[node->nchildren] = data[orignchildren]; /* move constant from old end to new end */
      for( i = orignchildren; i < node->nchildren; ++i )
         data[i] = 1.0;
      node->data.data = (void*)data;
   }
   else if( node->op == SCIP_EXPR_QUADRATIC && node->nchildren > orignchildren )
   {
      /* if quadratic expression, then add 0.0 linear coefficients for new expressions */
      SCIP_EXPRDATA_QUADRATIC* data;

      data = (SCIP_EXPRDATA_QUADRATIC*)node->data.data;
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &data->lincoefs, orignchildren, node->nchildren) );
      BMSclearMemoryArray(&data->lincoefs[orignchildren], node->nchildren - orignchildren);  /*lint !e866*/
   }

   node->simplified = FALSE;

   return SCIP_OKAY;
}

/** replaces a child node by another node
 *
 *  Assumes that both nodes represent the same expression.
 *  If this node was the last parent of oldchild and oldchild is not in use, then it is freed.
 *  newchild must have deeper depth than node.
 */
static
SCIP_RETCODE exprgraphNodeReplaceChild(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< pointer to expression graph node */
   SCIP_EXPRGRAPHNODE**  oldchild,           /**< child node that should be replaced, it may be freed */
   SCIP_EXPRGRAPHNODE*   newchild            /**< node that should take position of oldchild */
   )
{
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(oldchild != NULL);
   assert(*oldchild != NULL);
   assert(newchild != NULL);

   if( *oldchild == newchild )
      return SCIP_OKAY;

   SCIPdebugMessage("replace child %p in node %p by %p\n", (void*)*oldchild, (void*)node, (void*)newchild);

   /* search for oldchild in children array */
   for( i = 0; i < node->nchildren; ++i )
   {
      if( node->children[i] == *oldchild )
      {
         /* add as parent to newchild */
         SCIP_CALL( exprgraphNodeAddParent(exprgraph->blkmem, newchild, node) );

         /* remove as parent from oldchild */
         SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, oldchild, node) );

         /* set newchild as child i */
         node->children[i] = newchild;

         /* we're done */
         break;
      }
   }
   assert(i < node->nchildren); /* assert that oldchild has been found in children array */

   node->simplified = FALSE;

   return SCIP_OKAY;
}

/** comparison of SCIP_EXPRGRAPHNODE's that are of type SCIP_EXPR_CONST
 *
 *  A node is larger than another node, if their corresponding constants are related that way.
 */
static
SCIP_DECL_SORTPTRCOMP(exprgraphConstNodeComp)
{
   assert(elem1 != NULL);
   assert(elem2 != NULL);
   assert(((SCIP_EXPRGRAPHNODE*)elem1)->op == SCIP_EXPR_CONST);
   assert(((SCIP_EXPRGRAPHNODE*)elem2)->op == SCIP_EXPR_CONST);
   assert(((SCIP_EXPRGRAPHNODE*)elem1)->data.dbl == ((SCIP_EXPRGRAPHNODE*)elem1)->data.dbl); /* assert that const value is not nan */  /*lint !e777*/
   assert(((SCIP_EXPRGRAPHNODE*)elem2)->data.dbl == ((SCIP_EXPRGRAPHNODE*)elem2)->data.dbl); /* assert that const value is not nan */  /*lint !e777*/

   if( ((SCIP_EXPRGRAPHNODE*)elem1)->data.dbl > ((SCIP_EXPRGRAPHNODE*)elem2)->data.dbl )
      return 1;
   else if( ((SCIP_EXPRGRAPHNODE*)elem1)->data.dbl < ((SCIP_EXPRGRAPHNODE*)elem2)->data.dbl )
      return -1;
   else
      return 0;
}

/** sort array of nodes that holds constants */
static
void exprgraphSortConstNodes(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   )
{
   assert(exprgraph != NULL);

   if( exprgraph->constssorted )
      return;

   SCIPsortPtr((void**)exprgraph->constnodes, exprgraphConstNodeComp, exprgraph->nconsts);

   exprgraph->constssorted = TRUE;
}

/** finds position of expression graph node corresponding to a constant in constnodes array */
static
SCIP_Bool exprgraphFindConstNodePos(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< node to search for */
   int*                  pos                 /**< buffer to store position of node, if found */
   )
{
   int left;
   int right;
   int middle;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_CONST);
   assert(node->depth == 0);
   assert(node->pos >= 0);
   assert(pos != NULL);

   exprgraphSortConstNodes(exprgraph);
   assert(exprgraph->constssorted);

   /* find a node with constant node->data.dbl using binary search */
   left = 0;
   right = exprgraph->nconsts-1;
   *pos = -1;
   while( left <= right )
   {
      middle = (left+right)/2;
      assert(0 <= middle && middle < exprgraph->nconsts);

      if( node->data.dbl < exprgraph->constnodes[middle]->data.dbl )
         right = middle - 1;
      else if( node->data.dbl > exprgraph->constnodes[middle]->data.dbl )
         left  = middle + 1;
      else
      {
         *pos = middle;
         break;
      }
   }
   assert(left == right+1 || *pos >= 0);
   if( left == right+1 )
      return FALSE;

   /* search left of *pos to find node */
   while( exprgraph->constnodes[*pos] != node && *pos > 0 && exprgraph->constnodes[*pos-1]->data.dbl == node->data.dbl )  /*lint !e777*/
      --*pos;
   /* search right of *pos to find node */
   while( exprgraph->constnodes[*pos] != node && *pos < exprgraph->nconsts-1 && exprgraph->constnodes[*pos+1]->data.dbl == node->data.dbl )  /*lint !e777*/
      ++*pos;

   return exprgraph->constnodes[*pos] == node;
}

/** creates an expression graph node */
static
SCIP_RETCODE exprgraphCreateNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   SCIP_EXPROP           op,                 /**< operator type of expression */
   SCIP_EXPROPDATA       opdata              /**< operator data of expression */
   )
{
   assert(blkmem != NULL);
   assert(node   != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, node) );
   BMSclearMemory(*node);

   (*node)->op   = op;
   (*node)->data = opdata;

   /* mark graph position as not in graph yet */
   (*node)->depth = -1;
   (*node)->pos   = -1;

   /* arrays of length 0 are trivially sorted */
   (*node)->parentssorted  = TRUE;

   /* set bounds interval to entire */
   (*node)->boundstatus = SCIP_EXPRBOUNDSTATUS_VALID;
   SCIPintervalSetEntire(SCIP_REAL_MAX, &(*node)->bounds);

   /* set initial value to invalid */
   (*node)->value = SCIP_INVALID;

   /* set initial curvature to linear for variables, parameters, and constants and unknown otherwise */
   if( op == SCIP_EXPR_VARIDX || op == SCIP_EXPR_CONST || op == SCIP_EXPR_PARAM )
      (*node)->curv = SCIP_EXPRCURV_LINEAR;
   else
      (*node)->curv = SCIP_EXPRCURV_UNKNOWN;

   /* per default, a node is enabled */
   (*node)->enabled = TRUE;

   return SCIP_OKAY;
}

/** prints the expression corresponding to a node (not recursively) */
static
void exprgraphPrintNodeExpression(
   SCIP_EXPRGRAPHNODE*   node,               /**< node of expression graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   const char**          varnames,           /**< variable names, or NULL for generic names */
   SCIP_Bool             printchildrenbounds /**< whether to print bounds of children */
   )
{
   int i;

   assert(node != NULL);

   switch( node->op )
   {
   case SCIP_EXPR_VARIDX:
      if( varnames != NULL )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "<%s>", (const char*)varnames[node->data.intval]);
      }
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "x%d", node->data.intval);
      break;

   case SCIP_EXPR_CONST:
      SCIPmessageFPrintInfo(messagehdlr, file, "%g", node->data.dbl);
      break;

   case SCIP_EXPR_PARAM:
      SCIPmessageFPrintInfo(messagehdlr, file, "param%d", node->data.intval);
      break;

   case SCIP_EXPR_PLUS:
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c0[%10g,%10g]", node->children[0]->bounds.inf, node->children[0]->bounds.sup);
      SCIPmessageFPrintInfo(messagehdlr, file, "+");
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c1[%10g,%10g]", node->children[1]->bounds.inf, node->children[1]->bounds.sup);
      break;

   case SCIP_EXPR_MINUS:
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c0[%10g,%10g]", node->children[0]->bounds.inf, node->children[0]->bounds.sup);
      SCIPmessageFPrintInfo(messagehdlr, file, "-");
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c1[%10g,%10g]", node->children[1]->bounds.inf, node->children[1]->bounds.sup);
      break;

   case SCIP_EXPR_MUL:
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c0[%10g,%10g]", node->children[0]->bounds.inf, node->children[0]->bounds.sup);
      SCIPmessageFPrintInfo(messagehdlr, file, "*");
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c1[%10g,%10g]", node->children[1]->bounds.inf, node->children[1]->bounds.sup);
      break;

   case SCIP_EXPR_DIV:
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c0[%10g,%10g]", node->children[0]->bounds.inf, node->children[0]->bounds.sup);
      SCIPmessageFPrintInfo(messagehdlr, file, "/");
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c1[%10g,%10g]", node->children[1]->bounds.inf, node->children[1]->bounds.sup);
      break;

   case SCIP_EXPR_SQUARE:
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c0[%10g,%10g]", node->children[0]->bounds.inf, node->children[0]->bounds.sup);
      SCIPmessageFPrintInfo(messagehdlr, file, "^2");
      break;

   case SCIP_EXPR_REALPOWER:
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "c0[%10g,%10g]", node->children[0]->bounds.inf, node->children[0]->bounds.sup);
      SCIPmessageFPrintInfo(messagehdlr, file, "^%g", node->data.dbl);
      break;

   case SCIP_EXPR_SIGNPOWER:
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "sign(c0)|c0[%10g,%10g]|^%g",
            node->children[0]->bounds.inf, node->children[0]->bounds.sup, node->data.dbl);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "sign(c0)|c0|^%g", node->data.dbl);
      break;

   case SCIP_EXPR_INTPOWER:
      SCIPmessageFPrintInfo(messagehdlr, file, "c0");
      if( printchildrenbounds )
         SCIPmessageFPrintInfo(messagehdlr, file, "[%10g,%10g]", node->children[0]->bounds.inf, node->children[0]->bounds.sup);
      SCIPmessageFPrintInfo(messagehdlr, file, "^%d", node->data.intval);
      break;

   case SCIP_EXPR_SQRT:
   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* SCIP_EXPR_ERF       = 20, */  /**< gaussian error function (1 operand) */
      /* SCIP_EXPR_ERFI      = 21, */  /**< imaginary part of gaussian error function (1 operand) */
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
      SCIPmessageFPrintInfo(messagehdlr, file, "%s", (const char*)SCIPexpropGetName(node->op));
      if( printchildrenbounds )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "(c0[%10g,%10g]", node->children[0]->bounds.inf, node->children[0]->bounds.sup);
         if( node->nchildren == 2 )
            SCIPmessageFPrintInfo(messagehdlr, file, ",c1[%10g,%10g]", node->children[1]->bounds.inf, node->children[1]->bounds.sup);
         SCIPmessageFPrintInfo(messagehdlr, file, ")");
      }
      break;

   case SCIP_EXPR_SUM:
      if( printchildrenbounds )
         for( i = 0; i < node->nchildren; ++i )
         {
            if( i > 0 )
               SCIPmessageFPrintInfo(messagehdlr, file, "+");
            SCIPmessageFPrintInfo(messagehdlr, file, "c%d[%10g,%10g]", i, node->children[i]->bounds.inf, node->children[i]->bounds.sup);
         }
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "+");
      break;

   case SCIP_EXPR_PRODUCT:
      if( printchildrenbounds )
         for( i = 0; i < node->nchildren; ++i )
         {
            if( i > 0 )
               SCIPmessageFPrintInfo(messagehdlr, file, "*");
            SCIPmessageFPrintInfo(messagehdlr, file, "c%d[%10g,%10g]", i, node->children[i]->bounds.inf, node->children[i]->bounds.sup);
         }
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "*");
      break;

   case SCIP_EXPR_LINEAR:
   {
      SCIP_Real constant;

      constant = ((SCIP_Real*)node->data.data)[node->nchildren];

      if( constant != 0.0 || node->nchildren == 0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "%g", constant);

      for( i = 0; i < node->nchildren; ++i )
      {
         if( ((SCIP_Real*)node->data.data)[i] == 1.0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "+");
         else if( ((SCIP_Real*)node->data.data)[i] == -1.0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "-");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, "%+g*", ((SCIP_Real*)node->data.data)[i]);
         SCIPmessageFPrintInfo(messagehdlr, file, "c%d", i);
         if( printchildrenbounds )
            SCIPmessageFPrintInfo(messagehdlr, file, "[%10g,%10g]", node->children[i]->bounds.inf, node->children[i]->bounds.sup);
      }

      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* quadraticdata;

      quadraticdata = (SCIP_EXPRDATA_QUADRATIC*)node->data.data;
      assert(quadraticdata != NULL);

      if( quadraticdata->constant != 0.0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "%g", quadraticdata->constant);

      if( quadraticdata->lincoefs != NULL )
         for( i = 0; i < node->nchildren; ++i )
         {
            if( quadraticdata->lincoefs[i] == 0.0 )
               continue;
            SCIPmessageFPrintInfo(messagehdlr, file, "%+g*c%d", quadraticdata->lincoefs[i], i);
            if( printchildrenbounds )
               SCIPmessageFPrintInfo(messagehdlr, file, "[%10g,%10g]", node->children[i]->bounds.inf, node->children[i]->bounds.sup);
         }

      for( i = 0; i < quadraticdata->nquadelems; ++i )
      {
         if( quadraticdata->quadelems[i].coef == 1.0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "+");
         else if( quadraticdata->quadelems[i].coef == -1.0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "-");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, "%+g*", quadraticdata->quadelems[i].coef);
         SCIPmessageFPrintInfo(messagehdlr, file, "c%d", quadraticdata->quadelems[i].idx1);
         if( printchildrenbounds )
            SCIPmessageFPrintInfo(messagehdlr, file, "[%10g,%10g]", node->children[quadraticdata->quadelems[i].idx1]->bounds.inf, node->children[quadraticdata->quadelems[i].idx1]->bounds.sup);
         if( quadraticdata->quadelems[i].idx1 == quadraticdata->quadelems[i].idx2 )
            SCIPmessageFPrintInfo(messagehdlr, file, "^2");
         else
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "*c%d", quadraticdata->quadelems[i].idx2);
            if( printchildrenbounds )
               SCIPmessageFPrintInfo(messagehdlr, file, "[%10g,%10g]", node->children[quadraticdata->quadelems[i].idx2]->bounds.inf, node->children[quadraticdata->quadelems[i].idx2]->bounds.sup);
         }
      }

      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL*   monomialdata;
      int j;

      polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data;
      assert(polynomialdata != NULL);

      if( polynomialdata->constant != 0.0 || polynomialdata->nmonomials == 0 )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "%g", polynomialdata->constant);
      }

      for( i = 0; i < polynomialdata->nmonomials; ++i )
      {
         monomialdata = polynomialdata->monomials[i];
         if( monomialdata->coef == 1.0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "+");
         else if( monomialdata->coef == -1.0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "-");
         else
            SCIPmessageFPrintInfo(messagehdlr, file, "%+g", monomialdata->coef);

         for( j = 0; j < monomialdata->nfactors; ++j )
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "c%d", monomialdata->childidxs[j]);
            if( printchildrenbounds )
               SCIPmessageFPrintInfo(messagehdlr, file, "[%10g,%10g]", node->children[monomialdata->childidxs[j]]->bounds.inf, node->children[monomialdata->childidxs[j]]->bounds.sup);
            if( monomialdata->exponents[j] < 0.0 )
               SCIPmessageFPrintInfo(messagehdlr, file, "^(%g)", monomialdata->exponents[j]);
            else if( monomialdata->exponents[j] != 1.0 )
               SCIPmessageFPrintInfo(messagehdlr, file, "^%g", monomialdata->exponents[j]);
         }
      }

      break;
   }

   case SCIP_EXPR_LAST:
      SCIPABORT();
      break;

   default:
      SCIPmessageFPrintInfo(messagehdlr, file, "%s", SCIPexpropGetName(node->op));
      break;
   } /*lint !e788*/
}

/** prints a node of an expression graph */
static
void exprgraphPrintNodeDot(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< node of expression graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   const char**          varnames            /**< variable names, or NULL for generic names */
   )
{
   SCIP_Real color;
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(file != NULL);

   color = (SCIP_Real)node->op / (SCIP_Real)SCIP_EXPR_LAST;
   SCIPmessageFPrintInfo(messagehdlr, file, "n%d_%d [fillcolor=\"%g,%g,%g\", label=\"", node->depth, node->pos, color, color, color);

   exprgraphPrintNodeExpression(node, messagehdlr, file, varnames, FALSE);

   SCIPmessageFPrintInfo(messagehdlr, file, "\\n[%g,%g]", node->bounds.inf, node->bounds.sup);
   if( node->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDRELAXED )
      SCIPmessageFPrintInfo(messagehdlr, file, "!");
   if( node->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDTIGHTENED )
      SCIPmessageFPrintInfo(messagehdlr, file, "*");
   if( node->boundstatus & SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENT )
      SCIPmessageFPrintInfo(messagehdlr, file, "+");

   SCIPmessageFPrintInfo(messagehdlr, file, "\"");

   if( !node->enabled )
      SCIPmessageFPrintInfo(messagehdlr, file, ", style=dotted");

   SCIPmessageFPrintInfo(messagehdlr, file, "]\n");

   /* add edges from node to children */
   for( i = 0; i < node->nchildren; ++i )
      SCIPmessageFPrintInfo(messagehdlr, file, "n%d_%d -> n%d_%d [label=\"c%d\"]\n", node->depth, node->pos, node->children[i]->depth, node->children[i]->pos, i);
}

/** evaluate node of expression graph w.r.t. values stored in children */
static
SCIP_RETCODE exprgraphNodeEval(
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_Real*            varvals             /**< values for variables */
   )
{
   int i;
   SCIP_Real  staticbuf[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_Real* buf;

   assert(node != NULL);

   /* if many children, get large enough memory to store argument values */
   if( node->nchildren > SCIP_EXPRESSION_MAXCHILDEST )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&buf, node->nchildren) );
   }
   else
   {
      buf = staticbuf;
   }

   /* get values of children */
   for( i = 0; i < node->nchildren; ++i )
   {
      assert(node->children[i]->value != SCIP_INVALID);  /*lint !e777*/
      buf[i] = node->children[i]->value;  /*lint !e644*/
   }

   /* evaluate this expression */
   assert(exprOpTable[node->op].eval != NULL);
   SCIP_CALL( exprOpTable[node->op].eval(node->data, node->nchildren, buf, varvals, NULL, &node->value) );
   assert(node->value != SCIP_INVALID);  /*lint !e777*/

   /* free memory, if allocated before */
   if( staticbuf != buf )
   {
      BMSfreeMemoryArray(&buf);
   }

   return SCIP_OKAY;
}

/** evaluates node including subtree */
static
SCIP_RETCODE exprgraphNodeEvalWithChildren(
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_Real*            varvals             /**< values for variables */
   )
{
   int i;

   assert(node != NULL);

   for( i = 0; i < node->nchildren; ++i )
   {
      SCIP_CALL( exprgraphNodeEvalWithChildren(node->children[i], varvals) );
   }

   SCIP_CALL( exprgraphNodeEval(node, varvals) );

   return SCIP_OKAY;
}

/** updates bounds of a node if a children has changed its bounds */
static
SCIP_RETCODE exprgraphNodeUpdateBounds(
   SCIP_EXPRGRAPHNODE*   node,               /**< node of expression graph */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Real             minstrength,        /**< minimal required relative bound strengthening in a node to trigger a bound recalculation in parent nodes */
   SCIP_Bool             parenttightenisinvalid /**< whether to consider bounds that have been tightened by parents as invalid */
   )
{
   SCIP_INTERVAL  childboundsstatic[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_INTERVAL* childbounds;
   SCIP_INTERVAL newbounds;
   int i;

   assert(node != NULL);
   assert(node->depth >= 1); /* node should be in graph and not be at depth 0 (i.e., no variable, constant, or parameter) */
   assert(node->pos >= 0);   /* node should be in graph */
   assert(node->op != SCIP_EXPR_VARIDX);
   assert(node->op != SCIP_EXPR_PARAM);

   /* if we still have valid bounds and also no child got a bound tightening, then nothing to do
    * if node is disabled, then also do nothing */
   if( node->boundstatus == SCIP_EXPRBOUNDSTATUS_VALID || !node->enabled )
      return SCIP_OKAY;

   /* if many children, get large enough memory to store children bounds */
   if( node->nchildren > SCIP_EXPRESSION_MAXCHILDEST )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&childbounds, node->nchildren) );
   }
   else
   {
      childbounds = childboundsstatic;
   }

   /* assemble bounds of children */
   for( i = 0; i < node->nchildren; ++i )
   {
      /* child should have valid and non-empty bounds */
      assert(!(node->children[i]->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDRELAXED));
      assert(!SCIPintervalIsEmpty(infinity, node->children[i]->bounds));

      childbounds[i] = node->children[i]->bounds;  /*lint !e644*/
   }

   /* call interval evaluation function for this operand */
   assert( exprOpTable[node->op].inteval != NULL );
   SCIPintervalSet(&newbounds, 0.0);
   SCIP_CALL( exprOpTable[node->op].inteval(infinity, node->data, node->nchildren, childbounds, NULL, NULL, &newbounds) );

   /* free memory, if allocated before */
   if( childbounds != childboundsstatic )
   {
      BMSfreeMemoryArray(&childbounds);
   }

   /* NOTE: if you change code below, please make analog changes also in SCIPexprgraphUpdateNodeBoundsCurvature */

   /* if bounds of a children were relaxed or our bounds were tightened by a (now possibly invalid) reverse propagation from a parent
    * and now our bounds are relaxed, then we have to propagate this upwards to ensure valid bounds
    *
    * if bounds were tightened (considerably), then tell this to those parents which think that they have valid bounds
    *
    * finally, if there was only a little tightening, then keep this updated bounds, but don't notify parents
    */
   if( (newbounds.inf < node->bounds.inf || newbounds.sup > node->bounds.sup) &&
      ((node->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDRELAXED) || ((node->boundstatus & SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENT) && parenttightenisinvalid)) )
   {
      for( i = 0; i < node->nparents; ++i )
         node->parents[i]->boundstatus = SCIP_EXPRBOUNDSTATUS_CHILDRELAXED;

      node->bounds = newbounds;
   }
   else if( isLbBetter(minstrength, newbounds.inf, node->bounds.inf, node->bounds.sup) ||
      (     isUbBetter(minstrength, newbounds.sup, node->bounds.inf, node->bounds.sup)) )
   {
      for( i = 0; i < node->nparents; ++i )
         node->parents[i]->boundstatus |= SCIP_EXPRBOUNDSTATUS_CHILDTIGHTENED;

      node->bounds = newbounds;
   }
   else
   {
      SCIPintervalIntersect(&node->bounds, node->bounds, newbounds);
   }

   SCIPdebugMessage("updated bounds of node %p (%d,%d) op %s to [%g,%g]\n", (void*)node, node->depth, node->pos, SCIPexpropGetName(node->op), node->bounds.inf, node->bounds.sup);

   /* node now has valid bounds */
   node->boundstatus = SCIP_EXPRBOUNDSTATUS_VALID;

   return SCIP_OKAY;
}

/** propagate bounds of a node into children by reverting the nodes expression */
static
void exprgraphNodePropagateBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< node in expression graph with no parents */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Real             minstrength,        /**< minimal required relative bound strengthening in a node to trigger a propagation into children nodes */
   SCIP_Bool*            cutoff              /**< buffer to store whether a node's bounds were propagated to an empty interval */
   )
{
   SCIP_INTERVAL childbounds;
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->depth >= 0); /* node should be in graph */
   assert(node->pos >= 0);   /* node should be in graph */
   assert(minstrength >= 0.0);
   assert(cutoff != NULL);
   assert(!SCIPintervalIsEmpty(infinity, node->bounds)); /* should not call backward prop. for a node that yield a cutoff already */
   assert(!node->enabled || !(node->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDRELAXED)); /* there should be no unprocessed relaxations of children bounds, if node is enabled */

   /* if we have no recent bound tightening from a parent, then no use in reverse-propagating our bounds */
   if( (node->boundstatus & SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTRECENT) == 0 )
      return;

   /* if node is not enabled, then do nothing */
   if( !node->enabled )
      return;

   /* tell children that they should propagate their bounds even if not tightened */
   if( (node->boundstatus & SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTFORCE) == SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTFORCE )
      minstrength = -1.0;

   /* we will do something, so reset boundstatus to "tightened-by-parent, but not recently" */
   node->boundstatus = SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENT;

   /* SCIPdebugMessage("propagating node %p (%d,%d) op %s: [%10g,%10g] = ", (void*)node, node->depth, node->pos, SCIPexpropGetName(node->op), node->bounds.inf, node->bounds.sup);
    * SCIPdebug( exprgraphPrintNodeExpression(node, messagehdlr, NULL, NULL, TRUE) );
    * SCIPdebugPrintf("\n");
    */

   /* @todo add callback to exprOpTable for this */

   switch( node->op )
   {
   case SCIP_EXPR_VARIDX:
   case SCIP_EXPR_CONST:
   case SCIP_EXPR_PARAM:
      /* cannot propagate bound changes further */
      break;

   case SCIP_EXPR_PLUS:
   {
      assert(node->nchildren == 2);
      /* f = c0 + c1 -> c0 = f - c1, c1 = f - c0 */

      SCIPintervalSub(infinity, &childbounds, node->bounds, node->children[1]->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      if( *cutoff )
         break;

      SCIPintervalSub(infinity, &childbounds, node->bounds, node->children[0]->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[1], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_MINUS:
   {
      assert(node->nchildren == 2);
      /* f = c0 - c1 -> c0 = f + c1, c1 = c0 - f */

      SCIPintervalAdd(infinity, &childbounds, node->bounds, node->children[1]->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      if( *cutoff )
         break;

      SCIPintervalSub(infinity, &childbounds, node->children[0]->bounds, node->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[1], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_MUL:
   {
      assert(node->nchildren == 2);
      /* f = c0 * c1 -> c0 = f / c1, c1 = f / c0 */

      SCIPintervalDiv(infinity, &childbounds, node->bounds, node->children[1]->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      if( *cutoff )
         break;

      SCIPintervalDiv(infinity, &childbounds, node->bounds, node->children[0]->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[1], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_DIV:
   {
      assert(node->nchildren == 2);
      /* f = c0 / c1 -> c0 = f * c1, c1 = c0 / f */

      SCIPintervalMul(infinity, &childbounds, node->bounds, node->children[1]->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      if( *cutoff )
         break;

      SCIPintervalDiv(infinity, &childbounds, node->children[0]->bounds, node->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[1], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_SQUARE:
   {
      assert(node->nchildren == 1);
      /* f = c0^2 -> c0 = sqrt(f) union -sqrt(f) */

      if( node->bounds.sup < 0.0 )
      {
         *cutoff = TRUE;
         break;
      }

      SCIPintervalSquareRoot(infinity, &childbounds, node->bounds);
      if( node->children[0]->bounds.inf <= -childbounds.inf )
         SCIPintervalSetBounds(&childbounds, -childbounds.sup, childbounds.sup);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_SQRT:
   {
      assert(node->nchildren == 1);
      /* f = sqrt(c0) -> c0 = f^2 */

      SCIPintervalSquare(infinity, &childbounds, node->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_REALPOWER:
   {
      assert(node->nchildren == 1);

      SCIPintervalPowerScalarInverse(infinity, &childbounds, node->children[0]->bounds, node->data.dbl, node->bounds);

      if( SCIPintervalIsEmpty(infinity, childbounds) )
      {
         *cutoff = TRUE;
         break;
      }
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_SIGNPOWER:
   {
      assert(node->nchildren == 1);

      if( node->data.dbl != 0.0 )
      {
         SCIPintervalSignPowerScalar(infinity, &childbounds, node->bounds, 1.0/node->data.dbl);
      }
      else
      {
         /* behaves like SCIP_EXPR_SIGN */
         SCIPintervalSetBounds(&childbounds,
            (node->bounds.inf <= -1.0 && node->bounds.sup >= -1.0) ? -infinity : 0.0,
            (node->bounds.inf <=  1.0 && node->bounds.sup >=  1.0) ?  infinity : 0.0);
      }

      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_INTPOWER:
   {
      assert(node->nchildren == 1);

      SCIPintervalPowerScalarInverse(infinity, &childbounds, node->children[0]->bounds, (SCIP_Real)node->data.intval, node->bounds);

      if( SCIPintervalIsEmpty(infinity, childbounds) )
      {
         *cutoff = TRUE;
         break;
      }
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_EXP:
   {
      assert(node->nchildren == 1);
      /* f = exp(c0) -> c0 = log(f) */

      if( node->bounds.sup < 0.0 )
      {
         *cutoff = TRUE;
         break;
      }

      SCIPintervalLog(infinity, &childbounds, node->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_LOG:
   {
      assert(node->nchildren == 1);
      /* f = log(c0) -> c0 = exp(f) */

      SCIPintervalExp(infinity, &childbounds, node->bounds);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   {
      assert(node->nchildren == 1);

      /* @todo implement */

      break;
   }

   case SCIP_EXPR_ABS:
   {
      assert(node->nchildren == 1);

      /* use identity if child bounds are non-negative */
      if( node->children[0]->bounds.inf >= 0 )
      {
         SCIPintervalSetBounds(&childbounds, node->bounds.inf, node->bounds.sup);
      }
      /* use -identity if child bounds are non-positive */
      else if( node->children[0]->bounds.sup <= 0 )
      {
         assert(node->bounds.inf <= node->bounds.sup);
         SCIPintervalSetBounds(&childbounds, -node->bounds.sup, -node->bounds.inf);
      }
      /* f = |c0| -> c0 = -f union f = [-f.sup, f.sup] */
      else
      {
         SCIPintervalSetBounds(&childbounds, -node->bounds.sup, node->bounds.sup);
      }

      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_SIGN:
   {
      assert(node->nchildren == 1);
      /* f = sign(c0) -> c0 = ([-infty,0] if -1 in f) union ([0,infty] if 1 in f) */

      SCIPintervalSetBounds(&childbounds,
         (node->bounds.inf <= -1.0 && node->bounds.sup >= -1.0) ? -infinity : 0.0,
         (node->bounds.inf <=  1.0 && node->bounds.sup >=  1.0) ?  infinity : 0.0);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_MIN:
   {
      assert(node->nchildren == 2);
      /* f = min(c0,c1) -> f <= c0, f <= c1
       * if c1 > f -> c0 = f
       * if c0 > f -> c1 = f
       */

      SCIPintervalSetBounds(&childbounds, node->bounds.inf,
         node->children[1]->bounds.inf > node->bounds.sup ? node->bounds.sup : infinity);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      if( *cutoff )
         break;

      SCIPintervalSetBounds(&childbounds, node->bounds.inf,
         node->children[0]->bounds.inf > node->bounds.sup ? node->bounds.sup : infinity);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[1], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_MAX:
   {
      assert(node->nchildren == 2);
      /* f = max(c0, c1) -> f >= c0, f >= c1
       * if c1 < f -> c0 = f
       * if c0 < f -> c1 = f
       */

      SCIPintervalSetBounds(&childbounds,
         node->children[1]->bounds.sup < node->bounds.inf ? node->bounds.inf : -infinity,
         node->bounds.sup);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

      SCIPintervalSetBounds(&childbounds,
         node->children[0]->bounds.sup < node->bounds.inf ? node->bounds.inf : -infinity,
         node->bounds.sup);
      SCIPexprgraphTightenNodeBounds(exprgraph, node->children[1], childbounds, minstrength, infinity, cutoff);

      break;
   }

   case SCIP_EXPR_SUM:
   {
      SCIP_ROUNDMODE prevroundmode;

      /* f = sum_i c_i -> c_i = f - sum_{j,j!=i} c_j */

      SCIP_Real minlinactivity;
      SCIP_Real maxlinactivity;
      int minlinactivityinf;
      int maxlinactivityinf;

      if( node->nchildren == 0 )
         break;

      if( SCIPintervalIsEntire(infinity, node->bounds) )
         break;

      minlinactivity = 0.0;
      maxlinactivity = 0.0;
      minlinactivityinf = 0;
      maxlinactivityinf = 0;

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      for( i = 0; i < node->nchildren; ++i )
      {
         assert(!SCIPintervalIsEmpty(infinity, node->children[i]->bounds));

         /* minimal activity is only useful if node has a finite upper bound */
         if( node->bounds.sup < infinity )
         {
            if( node->children[i]->bounds.inf <= -infinity )
            {
               ++minlinactivityinf;
            }
            else
            {
               assert(node->children[i]->bounds.inf < infinity);
               minlinactivity += node->children[i]->bounds.inf;
            }
         }

         /* maximal activity is only useful if node has a finite lower bound
          * we compute negated maximal activity here so we can keep downward rounding
          */
         if( node->bounds.inf > -infinity )
         {
            if( node->children[i]->bounds.sup >= infinity )
            {
               ++maxlinactivityinf;
            }
            else
            {
               assert(node->children[i]->bounds.sup > -infinity);
               maxlinactivity -= node->children[i]->bounds.sup;
            }
         }
      }
      maxlinactivity = -maxlinactivity; /* correct sign */

      /* if there are too many unbounded bounds, then could only compute infinite bounds for children, so give up */
      if( (minlinactivityinf >= 2 || node->bounds.sup >=  infinity) &&
         ( maxlinactivityinf >= 2 || node->bounds.inf <= -infinity)
         )
      {
         SCIPintervalSetRoundingMode(prevroundmode);
         break;
      }

      for( i = 0; i < node->nchildren && !*cutoff; ++i )
      {
         /* upper bounds of c_i is
          *   node->bounds.sup - (minlinactivity - c_i.inf), if c_i.inf > -infinity and minlinactivityinf == 0
          *   node->bounds.sup - minlinactivity, if c_i.inf == -infinity and minlinactivityinf == 1
          */
         SCIPintervalSetEntire(infinity, &childbounds);
         if( node->bounds.sup <  infinity )
         {
            /* we are still in downward rounding mode, so negate and negate to get upward rounding */
            if( node->children[i]->bounds.inf <= -infinity && minlinactivityinf <= 1 )
            {
               assert(minlinactivityinf == 1);
               childbounds.sup = SCIPintervalNegateReal(minlinactivity - node->bounds.sup);
            }
            else if( minlinactivityinf == 0 )
            {
               childbounds.sup = SCIPintervalNegateReal(minlinactivity - node->bounds.sup - node->children[i]->bounds.inf);
            }
         }

         /* lower bounds of c_i is
          *   node->bounds.inf - (maxlinactivity - c_i.sup), if c_i.sup < infinity and maxlinactivityinf == 0
          *   node->bounds.inf - maxlinactivity, if c_i.sup == infinity and maxlinactivityinf == 1
          */
         if( node->bounds.inf > -infinity )
         {
            if( node->children[i]->bounds.sup >= infinity && maxlinactivityinf <= 1 )
            {
               assert(maxlinactivityinf == 1);
               childbounds.inf = node->bounds.inf - maxlinactivity;
            }
            else if( maxlinactivityinf == 0 )
            {
               childbounds.inf = node->bounds.inf - maxlinactivity + node->children[i]->bounds.sup;
            }
         }

         SCIPexprgraphTightenNodeBounds(exprgraph, node->children[i], childbounds, minstrength, infinity, cutoff);
      }

      SCIPintervalSetRoundingMode(prevroundmode);

      break;
   }

   case SCIP_EXPR_PRODUCT:
   {
      int j;
      /* f = prod_i c_i -> c_i = f / prod_{j:j!=i} c_j */

      /* too expensive (runtime here is quadratic in number of children) */
      if( node->nchildren > 10 )
         break;

      /* useless */
      if( SCIPintervalIsEntire(infinity, node->bounds) )
         break;

      for( i = 0; i < node->nchildren && !*cutoff; ++i )
      {
         /* compute prod_{j:j!=i} c_j */
         SCIPintervalSet(&childbounds, 1.0);
         for( j = 0; j < node->nchildren; ++j )
         {
            if( i == j )
               continue;
            SCIPintervalMul(infinity, &childbounds, childbounds, node->children[j]->bounds);

            /* if there is 0.0 in the product, then later division will hardly give useful bounds, so giveup for this i */
            if( childbounds.inf <= 0.0 && childbounds.sup >= 0.0 )
               break;
         }

         if( j == node->nchildren )
         {
            SCIPintervalDiv(infinity, &childbounds, node->bounds, childbounds); /* f / prod_{j:j!=i} c_j */
            SCIPexprgraphTightenNodeBounds(exprgraph, node->children[i], childbounds, minstrength, infinity, cutoff);
         }
      }

      break;
   }

   case SCIP_EXPR_LINEAR:
   {
      SCIP_ROUNDMODE prevroundmode;
      SCIP_Real* coefs;

      /* f = constant + sum_i a_ic_i -> c_i = (f - constant - sum_{j,j!=i} a_jc_j) / c_i */

      SCIP_Real minlinactivity;
      SCIP_Real maxlinactivity;
      int minlinactivityinf;
      int maxlinactivityinf;

      if( node->nchildren == 0 )
         break;

      if( SCIPintervalIsEntire(infinity, node->bounds) )
         break;

      coefs = (SCIP_Real*)node->data.data;

      minlinactivity =  coefs[node->nchildren];
      maxlinactivity = -coefs[node->nchildren];
      minlinactivityinf = 0;
      maxlinactivityinf = 0;

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      for( i = 0; i < node->nchildren; ++i )
      {
         assert(!SCIPintervalIsEmpty(infinity, node->children[i]->bounds));

         /* minimal activity is only useful if node has a finite upper bound */
         if( node->bounds.sup < infinity )
         {
            if( coefs[i] >= 0.0 )
            {
               if( node->children[i]->bounds.inf <= -infinity )
               {
                  ++minlinactivityinf;
               }
               else
               {
                  assert(node->children[i]->bounds.inf < infinity);
                  minlinactivity += coefs[i] * node->children[i]->bounds.inf;
               }
            }
            else
            {
               if( node->children[i]->bounds.sup >= infinity )
               {
                  ++minlinactivityinf;
               }
               else
               {
                  assert(node->children[i]->bounds.sup > -infinity);
                  minlinactivity += coefs[i] * node->children[i]->bounds.sup;
               }
            }
         }

         /* maximal activity is only useful if node has a finite lower bound
          * we compute negated maximal activity here so we can keep downward rounding
          */
         if( node->bounds.inf > -infinity )
         {
            if( coefs[i] >= 0.0 )
            {
               if( node->children[i]->bounds.sup >= infinity )
               {
                  ++maxlinactivityinf;
               }
               else
               {
                  assert(node->children[i]->bounds.sup > -infinity);
                  maxlinactivity += SCIPintervalNegateReal(coefs[i]) * node->children[i]->bounds.sup;
               }
            }
            else
            {
               if( node->children[i]->bounds.inf <= -infinity )
               {
                  ++maxlinactivityinf;
               }
               else
               {
                  assert(node->children[i]->bounds.inf < infinity);
                  maxlinactivity += SCIPintervalNegateReal(coefs[i]) * node->children[i]->bounds.inf;
               }
            }
         }
      }
      maxlinactivity = SCIPintervalNegateReal(maxlinactivity); /* correct sign */

      /* SCIPdebugMessage("activity = [%10g,%10g] ninf = [%d,%d]; bounds = [%10g,%10g]\n", minlinactivity, maxlinactivity, minlinactivityinf, maxlinactivityinf, node->bounds.inf, node->bounds.sup); */

      /* if there are too many unbounded bounds, then could only compute infinite bounds for children, so give up */
      if( (minlinactivityinf >= 2 || node->bounds.sup >=  infinity) &&
         (maxlinactivityinf >= 2 || node->bounds.inf <= -infinity)
         )
      {
         SCIPintervalSetRoundingMode(prevroundmode);
         break;
      }

      for( i = 0; i < node->nchildren && !*cutoff; ++i )
      {
         SCIP_INTERVAL ac;

         if( coefs[i] == 0.0 )
            continue;

         /* contribution of child i to activity (coefs[i] * node->children[i]->bounds) */
         SCIPintervalSet(&ac, 0.0);
         if( coefs[i] >= 0.0 )
         {
            if( node->children[i]->bounds.inf > -infinity )
               ac.inf = coefs[i] * node->children[i]->bounds.inf;
            if( node->children[i]->bounds.sup < infinity )
               ac.sup = SCIPintervalNegateReal(SCIPintervalNegateReal(coefs[i]) * node->children[i]->bounds.sup);
         }
         else
         {
            if( node->children[i]->bounds.sup < infinity )
               ac.inf = coefs[i] * node->children[i]->bounds.sup;
            if( node->children[i]->bounds.inf > -infinity )
               ac.sup = -SCIPintervalNegateReal(coefs[i] * node->children[i]->bounds.inf);
         }

         SCIPintervalSetEntire(infinity, &childbounds);
         if( coefs[i] > 0.0 )
         {
            /* upper bounds of c_i is
             *   (node->bounds.sup - minlinactivity)/coefs[i] + c_i.inf, if c_i.inf > -infinity and minlinactivityinf == 0
             *   (node->bounds.sup - minlinactivity)/coefs[i], if c_i.inf == -infinity and minlinactivityinf == 1
             */
            if( node->bounds.sup <  infinity )
            {
               /* we are still in downward rounding mode, so negate to get upward rounding */
               if( node->children[i]->bounds.inf <= -infinity && minlinactivityinf <= 1 )
               {
                  assert(minlinactivityinf == 1);
                  childbounds.sup = SCIPintervalNegateReal((minlinactivity - node->bounds.sup)/coefs[i]);
               }
               else if( minlinactivityinf == 0 )
               {
                  childbounds.sup = SCIPintervalNegateReal((minlinactivity - ac.inf - node->bounds.sup)/coefs[i]);
               }
            }

            /* lower bounds of c_i is
             *   (node->bounds.inf - maxlinactivity)/coefs[i] + c_i.sup, if c_i.sup < infinity and maxlinactivityinf == 0
             *   (node->bounds.inf - maxlinactivity)/coefs[i], if c_i.sup == infinity and maxlinactivityinf == 1
             */
            if( node->bounds.inf > -infinity )
            {
               if( node->children[i]->bounds.sup >= infinity && maxlinactivityinf <= 1 )
               {
                  assert(maxlinactivityinf == 1);
                  childbounds.inf = (node->bounds.inf - maxlinactivity)/coefs[i];
               }
               else if( maxlinactivityinf == 0 )
               {
                  childbounds.inf = (node->bounds.inf - maxlinactivity + ac.sup)/coefs[i];
               }
            }
         }
         else
         {
            /* (a-b)/c in downward rounding may not result in a lower bound on (a-b)/c if c is negative
             * thus, we do (b-a)/(-c) in downward rounding
             */
            /* lower bounds of c_i is
             *   (node->bounds.sup - minlinactivity)/coefs[i] + c_i.sup, if c_i.sup < infinity and minlinactivityinf == 0
             *   (node->bounds.sup - minlinactivity)/coefs[i], if c_i.sup == infinity and minlinactivityinf == 1
             */
            if( node->bounds.sup <  infinity )
            {
               if( node->children[i]->bounds.sup >= infinity && minlinactivityinf <= 1 )
               {
                  assert(minlinactivityinf == 1);
                  childbounds.inf = (minlinactivity - node->bounds.sup)/SCIPintervalNegateReal(coefs[i]);
               }
               else if( minlinactivityinf == 0 )
               {
                  childbounds.inf = (minlinactivity - ac.inf - node->bounds.sup)/SCIPintervalNegateReal(coefs[i]);
               }
            }

            /* upper bounds of c_i is
             *   (node->bounds.inf - maxlinactivity)/coefs[i] + c_i.inf, if c_i.inf > -infinity and maxlinactivityinf == 0
             *   (node->bounds.inf - maxlinactivity)/coefs[i], if c_i.inf == -infinity and maxlinactivityinf == 1
             */
            if( node->bounds.inf > -infinity )
            {
               /* we are still in downward rounding mode, so negate to get upward rounding */
               if( node->children[i]->bounds.inf <= -infinity && maxlinactivityinf <= 1 )
               {
                  assert(maxlinactivityinf == 1);
                  childbounds.sup = SCIPintervalNegateReal((node->bounds.inf - maxlinactivity)/SCIPintervalNegateReal(coefs[i]));
               }
               else if( maxlinactivityinf == 0 )
               {
                  childbounds.sup = SCIPintervalNegateReal((node->bounds.inf - maxlinactivity + ac.sup)/SCIPintervalNegateReal(coefs[i]));
               }
            }
         }

         SCIPexprgraphTightenNodeBounds(exprgraph, node->children[i], childbounds, minstrength, infinity, cutoff);
      }

      SCIPintervalSetRoundingMode(prevroundmode);

      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* quaddata;
      SCIP_INTERVAL tmp;
      SCIP_INTERVAL a;
      SCIP_INTERVAL b;
      SCIP_INTERVAL c;
      SCIP_QUADELEM* quadelems;
      int nquadelems;
      SCIP_Real* lincoefs;
      int k;

      /* f = constant + sum_i c_i + sum_k a_k c_(i_k) c_(j_k)
       * turn into quadratic univariate equation a*c_i^2 + b*c_i = c for each child
       */

      quaddata   = (SCIP_EXPRDATA_QUADRATIC*)node->data.data;
      quadelems  = quaddata->quadelems;
      nquadelems = quaddata->nquadelems;
      lincoefs   = quaddata->lincoefs;

      /* too expensive, runtime here is O(nchildren * nquadelems) = O(nquadelems^2) since nchildren <= 2*nquadelems usually */
      if( nquadelems > 10 )
         break;

      if( SCIPintervalIsEntire(infinity, node->bounds) )
         break;

      if( node->nchildren == 2 && nquadelems > 0 )
      {
         /* if it's a bivariate quadratic expression with bilinear term, do something special */
         SCIP_Real ax;  /* square coefficient of first  child */
         SCIP_Real ay;  /* square coefficient of second child */
         SCIP_Real axy; /* bilinear coefficient */

         ax = 0.0;
         ay = 0.0;
         axy = 0.0;
         for( i = 0; i < nquadelems; ++i )
            if( quadelems[i].idx1 == 0 && quadelems[i].idx2 == 0 )
               ax += quadelems[i].coef;
            else if( quadelems[i].idx1 == 1 && quadelems[i].idx2 == 1 )
               ay += quadelems[i].coef;
            else
               axy += quadelems[i].coef;

         c = node->bounds;
         SCIPintervalSubScalar(infinity, &c, c, quaddata->constant);

         /* compute bounds for x */
         SCIPintervalSolveBivariateQuadExpressionAllScalar(
            infinity, &childbounds, ax, ay, axy,
            lincoefs != NULL ? lincoefs[0] : 0.0, lincoefs != NULL ? lincoefs[1] : 0.0,
            c, node->children[0]->bounds, node->children[1]->bounds
            );
         if( (childbounds.inf > node->children[0]->bounds.inf + 1e-9 || childbounds.sup + 1e-9 < node->children[0]->bounds.sup) )
         {
            SCIPdebugMessage("%g x^2 + %g y^2 + %g xy + %g x + %g y in [%g,%g], x = [%g,%g], y = [%g,%g] -> x in [%g,%g], cutoff = %d\n",
               ax, ay, axy, lincoefs != NULL ? lincoefs[0] : 0.0, lincoefs != NULL ? lincoefs[1] : 0.0,
               c.inf, c.sup, node->children[0]->bounds.inf, node->children[0]->bounds.sup,
               node->children[1]->bounds.inf, node->children[1]->bounds.sup, childbounds.inf, childbounds.sup, (int)SCIPintervalIsEmpty(infinity, childbounds)
               );
         }

         if( SCIPintervalIsEmpty(infinity, childbounds) )
            *cutoff = TRUE;
         else
            SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);
         if( *cutoff )
            break;

         /* compute bounds for y */
         SCIPintervalSolveBivariateQuadExpressionAllScalar(
            infinity, &childbounds, ay, ax, axy,
            lincoefs != NULL ? lincoefs[1] : 0.0, lincoefs != NULL ? lincoefs[0] : 0.0,
            c, node->children[1]->bounds, node->children[0]->bounds
            );

         if( (childbounds.inf > node->children[1]->bounds.inf + 1e-9 || childbounds.sup + 1e-9 < node->children[1]->bounds.sup) )
         {
            SCIPdebugMessage("%g x^2 + %g y^2 + %g xy + %g x + %g y in [%g,%g], x = [%g,%g], y = [%g,%g] -> y in [%g,%g], cutoff = %d\n",
               ax, ay, axy, lincoefs != NULL ? lincoefs[0] : 0.0, lincoefs != NULL ? lincoefs[1] : 0.0,
               c.inf, c.sup, node->children[0]->bounds.inf, node->children[0]->bounds.sup,
               node->children[1]->bounds.inf, node->children[1]->bounds.sup, childbounds.inf, childbounds.sup, (int)SCIPintervalIsEmpty(infinity, childbounds)
               );
         }

         if( SCIPintervalIsEmpty(infinity, childbounds) )
            *cutoff = TRUE;
         else
            SCIPexprgraphTightenNodeBounds(exprgraph, node->children[1], childbounds, minstrength, infinity, cutoff);
         if( *cutoff )
            break;

         break;
      }

      for( i = 0; i < node->nchildren && !*cutoff; ++i )
      {
         SCIPintervalSet(&a, 0.0);
         SCIPintervalSet(&b, lincoefs != NULL ? lincoefs[i] : 0.0);
         c = node->bounds;
         SCIPintervalSubScalar(infinity, &c, c, quaddata->constant);

         /* move linear terms not corresponding to i into c
          * @todo do this faster, see EXPR_LINEAR
          */
         if( lincoefs != NULL )
            for( k = 0; k < node->nchildren; ++k )
               if( i != k && lincoefs[k] != 0.0 )
               {
                  SCIPintervalMulScalar(infinity, &tmp, node->children[k]->bounds, lincoefs[k]);
                  SCIPintervalSub(infinity, &c, c, tmp);
               }

         for( k = 0; k < nquadelems; ++k )
         {
            if( quadelems[k].idx1 == i && quadelems[k].idx2 == i )
            {
               SCIPintervalAddScalar(infinity, &a, a, quadelems[k].coef);
            }
            else if( quadelems[k].idx1 == i )
            {
               SCIPintervalMulScalar(infinity, &tmp, node->children[quadelems[k].idx2]->bounds, quadelems[k].coef);
               SCIPintervalAdd(infinity, &b, b, tmp);
            }
            else if( quadelems[k].idx2 == i )
            {
               SCIPintervalMulScalar(infinity, &tmp, node->children[quadelems[k].idx1]->bounds, quadelems[k].coef);
               SCIPintervalAdd(infinity, &b, b, tmp);
            }
            else if( quadelems[k].idx1 == quadelems[k].idx2 )
            {
               SCIPintervalSquare(infinity, &tmp, node->children[quadelems[k].idx1]->bounds);
               SCIPintervalMulScalar(infinity, &tmp, tmp, quadelems[k].coef);
               SCIPintervalSub(infinity, &c, c, tmp);
            }
            else
            {
               SCIPintervalMul(infinity, &tmp, node->children[quadelems[k].idx1]->bounds, node->children[quadelems[k].idx2]->bounds);
               SCIPintervalMulScalar(infinity, &tmp, tmp, quadelems[k].coef);
               SCIPintervalSub(infinity, &c, c, tmp);
            }
         }

         SCIPdebugMessage("solve %gc%d^2 + [%10g,%10g]c%d = [%10g,%10g]\n",
            a.inf, i, b.inf, b.sup, i, c.inf, c.sup);
         SCIPintervalSolveUnivariateQuadExpression(infinity, &childbounds, a, b, c);
         if( SCIPintervalIsEmpty(infinity, childbounds) )
            *cutoff = TRUE;
         else
            SCIPexprgraphTightenNodeBounds(exprgraph, node->children[i], childbounds, minstrength, infinity, cutoff);
      }

      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL**  monomials;
      SCIP_EXPRDATA_MONOMIAL*   monomial;
      int nmonomials;
      int j;
      int k;
      SCIP_Real n;
      int nexpisdoublen;
      int nexpishalfn;
      char abc_flag;

      SCIP_INTERVAL monomialcoef;
      SCIP_INTERVAL tmp;
      SCIP_INTERVAL a;
      SCIP_INTERVAL b;
      SCIP_INTERVAL c;

      /* f = constant + sum_i coef_i prod_j c_{i_j}^e_{i_j}
       * for each child x, write as a*x^(2n) + b*x^n = c for some n!=0
       *
       * we determine n by setting n to the first exponent of x that we see
       * then we count how often we see x^(2n) and x^(n/2)
       * if the number of x^(n/2) exceeds the number of x^(2n), then we half n
       */

      polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data;
      monomials  = polynomialdata->monomials;
      nmonomials = polynomialdata->nmonomials;

      if( SCIPintervalIsEntire(infinity, node->bounds) )
         break;

      for( i = 0; i < node->nchildren && !*cutoff; ++i )
      {
         n = 0.0;
         nexpisdoublen = 0;
         nexpishalfn = 0;
         for( j = 0; j < nmonomials; ++j )
         {
            monomial = monomials[j];
            for( k = 0; k < monomial->nfactors; ++k )
            {
               if( monomial->childidxs[k] == i )
               {
                  if( n == 0.0 )
                     n = monomial->exponents[k];
                  else if( n == 2*monomial->exponents[k] )  /*lint !e777*/
                     ++nexpishalfn;
                  else if( 2*n == monomial->exponents[k] )  /*lint !e777*/
                     ++nexpisdoublen;
               }
            }
         }

         if( n == 0.0 )
         {
            /* child does not appear in polynomial -> cannot deduce bound */
            continue;
         }

         /* half n if there are more monomials with x^(n/2) than monomials with x^(2n) */
         if( nexpishalfn > nexpisdoublen )
            n /= 2.0;

         SCIPintervalSet(&a, 0.0);
         SCIPintervalSet(&b, 0.0);
         SCIPintervalSubScalar(infinity, &c, node->bounds, polynomialdata->constant);

         for( j = 0; j < nmonomials; ++j )
         {
            monomial = monomials[j];
            SCIPintervalSet(&monomialcoef, monomial->coef);
            abc_flag = 'c';
            for( k = 0; k < monomial->nfactors; ++k )
            {
               if( monomial->childidxs[k] == i )
               {
                  assert(abc_flag == 'c'); /* child should appear only once per monom */
                  if( n > 0.0 )
                  {
                     if( monomial->exponents[k] > 2.0*n )
                     {
                        abc_flag = 'a';
                        SCIPintervalPowerScalar(infinity, &tmp, node->children[i]->bounds, monomial->exponents[k] - 2.0*n);
                        SCIPintervalMul(infinity, &monomialcoef, monomialcoef, tmp);
                     }
                     else if( monomial->exponents[k] == 2*n )  /*lint !e777*/
                     {
                        abc_flag = 'a';
                     }
                     else if( monomial->exponents[k] > n )
                     {
                        abc_flag = 'b';
                        SCIPintervalPowerScalar(infinity, &tmp, node->children[i]->bounds, monomial->exponents[k] - n);
                        SCIPintervalMul(infinity, &monomialcoef, monomialcoef, tmp);
                     }
                     else if( monomial->exponents[k] == n )  /*lint !e777*/
                     {
                        abc_flag = 'b';
                     }
                     else
                     {
                        SCIPintervalPowerScalar(infinity, &tmp, node->children[i]->bounds, monomial->exponents[k]);
                        SCIPintervalMul(infinity, &monomialcoef, monomialcoef, tmp);
                     }
                  }
                  else
                  {
                     assert(n < 0.0);
                     if( monomial->exponents[k] < 2.0*n )
                     {
                        abc_flag = 'a';
                        SCIPintervalPowerScalar(infinity, &tmp, node->children[i]->bounds, monomial->exponents[k] - 2.0*n);
                        SCIPintervalMul(infinity, &monomialcoef, monomialcoef, tmp);
                     }
                     else if( monomial->exponents[k] == 2*n )  /*lint !e777*/
                     {
                        abc_flag = 'a';
                     }
                     else if( monomial->exponents[k] < n )
                     {
                        abc_flag = 'b';
                        SCIPintervalPowerScalar(infinity, &tmp, node->children[i]->bounds, monomial->exponents[k] - n);
                        SCIPintervalMul(infinity, &monomialcoef, monomialcoef, tmp);
                     }
                     else if( monomial->exponents[k] == n )  /*lint !e777*/
                     {
                        abc_flag = 'b';
                     }
                     else
                     {
                        SCIPintervalPowerScalar(infinity, &tmp, node->children[i]->bounds, monomial->exponents[k]);
                        SCIPintervalMul(infinity, &monomialcoef, monomialcoef, tmp);
                     }
                  }
               }
               else
               {
                  SCIPintervalPowerScalar(infinity, &tmp, node->children[monomial->childidxs[k]]->bounds, monomial->exponents[k]);
                  SCIPintervalMul(infinity, &monomialcoef, monomialcoef, tmp);
               }
            }

            if( abc_flag == 'a' )
            {
               SCIPintervalAdd(infinity, &a, a, monomialcoef);
               /* if monomialcoef is such that a exceeds value for infinity, then stop */
               if( a.inf >= infinity || a.sup <= -infinity )
                  break;
            }
            else if( abc_flag == 'b' )
            {
               SCIPintervalAdd(infinity, &b, b, monomialcoef);
               /* if monomialcoef is such that b exceeds value for infinity, then stop */
               if( b.inf >= infinity || b.sup <= -infinity )
                  break;
            }
            else
            {
               SCIPintervalSub(infinity, &c, c, monomialcoef);
               /* if monomialcoef is such that c exceeds value for infinity, then stop */
               if( c.inf >= infinity || c.sup <= -infinity )
                  break;
            }
         }

         /* if we run out of numbers (within -infinity,infinity) above, then stop */
         if( j < nmonomials )
            continue;

         /* now have equation a*child^(2n) + b*child^n = c
          * solve a*y^2 + b*y = c, then child^n = y
          */
         SCIPdebugMessage("solve [%10g,%10g]c%d^%g + [%10g,%10g]c%d^%g = [%10g,%10g]",
            a.inf, a.sup, i, 2*n, b.inf, b.sup, i, n, c.inf, c.sup);
         SCIPintervalSolveUnivariateQuadExpression(infinity, &tmp, a, b, c);
         SCIPdebugPrintf(" -> c%d^%g = [%10g, %10g]", i, n, tmp.inf, tmp.sup);

         if( SCIPintervalIsEmpty(infinity, tmp) )
         {
            *cutoff = TRUE;
            break;
         }

         SCIPintervalPowerScalarInverse(infinity, &childbounds, node->children[i]->bounds, n, tmp);
         SCIPdebugPrintf(" -> c%d = [%10g, %10g]\n", i, childbounds.inf, childbounds.sup);
         if( SCIPintervalIsEmpty(infinity, childbounds) )
         {
            SCIPdebugMessage(" -> cutoff\n");
            *cutoff = TRUE;
            break;
         }

         SCIPexprgraphTightenNodeBounds(exprgraph, node->children[i], childbounds, minstrength, infinity, cutoff);

         /* SCIPdebugMessage("-> node %p (%d,%d): [%10g,%10g] = ", (void*)node, node->depth, node->pos, node->bounds.inf, node->bounds.sup);
            SCIPdebug( exprgraphPrintNodeExpression(node, messagehdlr, NULL, NULL, TRUE) );
            SCIPdebugPrintf("\n"); */
      }

      break;
   }

   case SCIP_EXPR_USER:
   {
      SCIP_INTERVAL* childrenbounds;
      SCIP_EXPRDATA_USER* exprdata;
      int c;

      exprdata = (SCIP_EXPRDATA_USER*)node->data.data;

      /* do nothing if callback not implemented */
      if( exprdata->prop == NULL )
         break;

      /* if only one child, do faster */
      if( node->nchildren == 1 )
      {
         childbounds = node->children[0]->bounds;
         SCIP_CALL_ABORT( exprdata->prop(infinity, exprdata->userdata, 1, &childbounds, node->bounds, cutoff) );

         if( !*cutoff )
            SCIPexprgraphTightenNodeBounds(exprgraph, node->children[0], childbounds, minstrength, infinity, cutoff);

         break;
      }

      SCIP_ALLOC_ABORT( BMSallocBlockMemoryArray(exprgraph->blkmem, &childrenbounds, node->nchildren) );
      for( c = 0; c < node->nchildren; ++c )
         childrenbounds[c] = node->children[c]->bounds;

      SCIP_CALL_ABORT( exprdata->prop(infinity, exprdata->userdata, node->nchildren, childrenbounds, node->bounds, cutoff) );

      for( c = 0; !*cutoff && c < node->nchildren; ++c )
      {
         SCIPexprgraphTightenNodeBounds(exprgraph, node->children[c], childrenbounds[c], minstrength, infinity, cutoff);
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &childrenbounds, node->nchildren);

      break;
   }

   case SCIP_EXPR_LAST:
      SCIPABORT();
      break;
   }
}

/** removes duplicate children in a polynomial expression node
 *
 *  Leaves NULL's in children array.
 */
static
SCIP_RETCODE exprgraphNodeRemovePolynomialDuplicateChildren(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   SCIP_Bool foundduplicates;
   int* childmap;
   int i;
   int j;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_POLYNOMIAL);

   if( node->nchildren == 0 )
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childmap, node->nchildren) );

   foundduplicates = FALSE;
   for( i = 0; i < node->nchildren; ++i )
   {
      if( node->children[i] == NULL )
         continue;
      childmap[i] = i;  /*lint !e644*/

      for( j = i+1; j < node->nchildren; ++j )
      {
         if( node->children[j] == NULL )
            continue;

         if( node->children[i] == node->children[j] )
         {
            /* node should be parent of children[j] at least twice,
             * so we remove it once
             */
            SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &node->children[j], node) );
            node->children[j] = NULL;
            assert(exprgraphNodeIsParent(node->children[i], node));

            childmap[j] = i;
            foundduplicates = TRUE;
         }
      }
   }

   /* apply childmap to monomials */
   if( foundduplicates )
      polynomialdataApplyChildmap((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data, childmap);

   /* free childmap */
   BMSfreeBlockMemoryArray(exprgraph->blkmem, &childmap, node->nchildren);

   return SCIP_OKAY;
}

/** eliminates NULL's in children array and shrinks it to actual size */
static
SCIP_RETCODE exprgraphNodeRemovePolynomialNullChildren(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   int* childmap;
   int lastnonnull;
   int i;

   assert(blkmem != NULL);
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_POLYNOMIAL);

   if( node->nchildren == 0 )
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &childmap, node->nchildren) );

   /* close gaps in children array */
   lastnonnull = node->nchildren-1;
   while( lastnonnull >= 0 && node->children[lastnonnull] == NULL )
      --lastnonnull;
   for( i = 0; i <= lastnonnull; ++i )
   {
      if( node->children[i] != NULL )
      {
         childmap[i] = i; /* child at index i is not moved */  /*lint !e644*/
         continue;
      }
      assert(node->children[lastnonnull] != NULL);

      /* move child at lastnonnull to position i */
      node->children[i] = node->children[lastnonnull];
      node->children[lastnonnull] = NULL;
      childmap[lastnonnull] = i;

      /* update lastnonnull */
      --lastnonnull;
      while( lastnonnull >= 0 && node->children[lastnonnull] == NULL )
         --lastnonnull;
   }
   assert(i > lastnonnull);

   /* apply childmap to monomials */
   if( lastnonnull < node->nchildren-1 )
      polynomialdataApplyChildmap((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data, childmap);

   BMSfreeBlockMemoryArray(blkmem, &childmap, node->nchildren);

   /* shrink children array */
   if( lastnonnull >= 0 )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &node->children, node->nchildren, lastnonnull+1) );
      node->nchildren = lastnonnull+1;
   }
   else
   {
      BMSfreeBlockMemoryArray(blkmem, &node->children, node->nchildren);
      node->nchildren = 0;
   }

   return SCIP_OKAY;
}

/** aims at simplifying a node in an expression graph, assuming all children have been simplified
 *
 *  Converts node into polynomial, if possible and not constant.
 */
static
SCIP_RETCODE exprgraphNodeSimplify(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Real             eps,                /**< threshold, under which positive values are treat as 0 */
   int                   maxexpansionexponent,/**< maximal exponent for which we still expand non-monomial polynomials */
   SCIP_Bool*            havechange          /**< flag to set if the node has been changed */
   )
{
   SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
   SCIP_EXPRDATA_MONOMIAL* monomial;
   BMS_BLKMEM* blkmem;
   SCIP_Bool removechild;
   SCIP_Bool* childinuse;
   int* childmap;
   int childmapsize;
   int i;
   int j;
   int orignchildren;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->depth > 0); /* simplifier is not thought for nodes at depth 0 */
   assert(havechange != NULL);

   blkmem = exprgraph->blkmem;
   assert(blkmem != NULL);

   SCIPdebugMessage("attempt simplification of node %p (%d,%d)\n", (void*)node, node->depth, node->pos);

   /* if all children are constants, then turn this node into constant */
   for( i = 0; i < node->nchildren; ++i )
      if( node->children[i]->op != SCIP_EXPR_CONST )
         break;
   if( node->nchildren > 0 && i == node->nchildren )
   {
      /* get value of node */
      SCIP_CALL( exprgraphNodeEvalWithChildren(node, NULL) );
      assert(node->value != SCIP_INVALID);  /*lint !e777*/

      SCIPdebugMessage("turn node %p (%d,%d) into constant %g\n", (void*)node, node->depth, node->pos, node->value);
      SCIPdebug( exprgraphPrintNodeExpression(node, messagehdlr, NULL, NULL, TRUE) );
      SCIPdebugPrintf("\n");

      /* free expression data */
      if( exprOpTable[node->op].freedata != NULL )
         exprOpTable[node->op].freedata(blkmem, node->nchildren, node->data);

      /* disconnect from children */
      for( i = 0; i < node->nchildren; ++i )
      {
         SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &node->children[i], node) );
      }
      BMSfreeBlockMemoryArray(blkmem, &node->children, node->nchildren);
      node->nchildren = 0;

      /* turn into constant expression */
      node->op = SCIP_EXPR_CONST;
      node->data.dbl = node->value;

      *havechange = TRUE;
      node->simplified = TRUE;

      return SCIP_OKAY;
   }

   /* @todo for sign, min, max, abs, knowing bounds on children may allow simplification
    * @todo log(product) -> sum(log)
    * @todo product(exp) -> exp(sum)
    * @todo exp(x)^p -> exp(p*x)
    * @todo exp(const*log(x)) -> x^const
    */

   SCIP_CALL( exprConvertToPolynomial(blkmem, &node->op, &node->data, node->nchildren) );

   if( node->op != SCIP_EXPR_POLYNOMIAL )
   {
      node->simplified = TRUE;
      return SCIP_OKAY;
   }

   polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data;
   assert(polynomialdata != NULL);

   orignchildren = node->nchildren;

   /* check if we have duplicate children and merge */
   SCIP_CALL( exprgraphNodeRemovePolynomialDuplicateChildren(exprgraph, node) );
   polynomialdataMergeMonomials(blkmem, polynomialdata, eps, TRUE);

   SCIPdebugMessage("expand factors in expression node ");
   SCIPdebug( exprgraphPrintNodeExpression(node, messagehdlr, NULL, NULL, FALSE) );
   SCIPdebugPrintf("\n");

   childmap = NULL;
   childmapsize = 0;

   /* resolve children that are constants
    * we do this first, because it reduces the degree and number of factors in the monomials,
    *   thereby allowing some expansions of polynomials that may not be possible otherwise, e.g., turning c0*c1 with c0=quadratic and c1=constant into a single monomial
    */
   for( i = 0; i < node->nchildren; ++i )
   {
      if( node->children[i] == NULL )
         continue;

      /* convert children to polynomial, if not constant or polynomial
       * if child was simplified in this round, it may have already been converted, and then nothing happens
       * but if child was already simplified, then it was not converted, and thus we try it here
       */
      if( node->children[i]->op != SCIP_EXPR_CONST )
         continue;

      SCIPdebugMessage("expand child %d in expression node ", i);
      SCIPdebug( exprgraphPrintNodeExpression(node, messagehdlr, NULL, NULL, FALSE) );
      SCIPdebugPrintf("\n\tchild = ");
      SCIPdebug( exprgraphPrintNodeExpression(node->children[i], messagehdlr, NULL, NULL, FALSE) );
      SCIPdebugPrintf("\n");

      removechild = TRUE; /* we intend to release children[i] */

      ensureBlockMemoryArraySize(blkmem, &childmap, &childmapsize, node->children[i]->nchildren);

      /* put constant of child i into every monomial where child i is used */
      for( j = 0; j < polynomialdata->nmonomials; ++j )
      {
         int factorpos;

         monomial = polynomialdata->monomials[j];
         /* if monomial is not sorted, then polynomial should not be sorted either, or have only one monomial */
         assert(monomial->sorted || !polynomialdata->sorted || polynomialdata->nmonomials <= 1);

         if( SCIPexprFindMonomialFactor(monomial, i, &factorpos) )
         {
            assert(factorpos >= 0);
            assert(factorpos < monomial->nfactors);
            /* assert that factors have been merged */
            assert(factorpos == 0 || monomial->childidxs[factorpos-1] != i);
            assert(factorpos == monomial->nfactors-1 || monomial->childidxs[factorpos+1] != i);

            SCIPdebugMessage("attempt expanding child %d at monomial %d factor %d\n", i, j, factorpos);

            if( !EPSISINT(monomial->exponents[factorpos], 0.0) && node->children[i]->data.dbl < 0.0 )  /*lint !e835*/
            {
               /* if constant is negative and our exponent is not integer, then cannot do expansion */
               SCIPmessagePrintWarning(messagehdlr, "got negative constant %g to the power of a noninteger exponent %g\n", node->children[i]->data.dbl, monomial->exponents[factorpos]);
               removechild = FALSE;
            }
            else
            {
               monomial->coef *= pow(node->children[i]->data.dbl, monomial->exponents[factorpos]);

               /* move last factor to position factorpos */
               if( factorpos < monomial->nfactors-1 )
               {
                  monomial->exponents[factorpos] = monomial->exponents[monomial->nfactors-1];
                  monomial->childidxs[factorpos] = monomial->childidxs[monomial->nfactors-1];
               }
               --monomial->nfactors;
               monomial->sorted = FALSE;
               polynomialdata->sorted = FALSE;

               *havechange = TRUE;
            }
         }
      }

      /* forget about child i, if it is not used anymore */
      if( removechild )
      {
         /* remove node from list of parents of child i */
         SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &node->children[i], node) );
         node->children[i] = NULL;
      }

      /* simplify current polynomial again */
      polynomialdataMergeMonomials(blkmem, polynomialdata, eps, TRUE);
   }

   /* resolve children that are polynomials itself */
   for( i = 0; i < node->nchildren; ++i )
   {
      if( node->children[i] == NULL )
         continue;

      /* convert children to polynomial, if not constant or polynomial
       * if child was simplified in this round, it may have already been converted, and then nothing happens
       * but if child was already simplified, then it was not converted, and thus we try it here
       */
      SCIP_CALL( exprConvertToPolynomial(blkmem, &node->children[i]->op, &node->children[i]->data, node->children[i]->nchildren) );

      if( node->children[i]->op != SCIP_EXPR_POLYNOMIAL )
         continue;

      SCIPdebugMessage("expand child %d in expression node %p = ", i, (void*)node);
      SCIPdebug( exprgraphPrintNodeExpression(node, messagehdlr, NULL, NULL, FALSE) );
      SCIPdebug( SCIPmessagePrintInfo(messagehdlr, "\n\tchild = ") );
      SCIPdebug( exprgraphPrintNodeExpression(node->children[i], messagehdlr, NULL, NULL, FALSE) );
      SCIPdebug( SCIPmessagePrintInfo(messagehdlr, "\n") );

      removechild = TRUE; /* we intend to release children[i] */

      ensureBlockMemoryArraySize(blkmem, &childmap, &childmapsize, node->children[i]->nchildren);

      /* add children of child i to node */
      SCIP_CALL( exprgraphNodeAddChildren(blkmem, node, node->children[i]->nchildren, node->children[i]->children, childmap) );

      /* put polynomial of child i into every monomial where child i is used */
      j = 0;
      while( j < polynomialdata->nmonomials )
      {
         int factorpos;
         SCIP_Bool success;

         monomial = polynomialdata->monomials[j];
         /* if monomial is not sorted, then polynomial should not be sorted either, or have only one monomial */
         assert(monomial->sorted || !polynomialdata->sorted || polynomialdata->nmonomials <= 1);

         /* make sure factors are merged, should only be potentially necessary if not sorted, see also #1848 */
         if( !monomial->sorted )
            SCIPexprMergeMonomialFactors(monomial, eps);

         if( !SCIPexprFindMonomialFactor(monomial, i, &factorpos) )
         {
            ++j;
            continue;
         }

         assert(factorpos >= 0);
         assert(factorpos < monomial->nfactors);
         /* assert that factors have been merged */
         assert(factorpos == 0 || monomial->childidxs[factorpos-1] != i);
         assert(factorpos == monomial->nfactors-1 || monomial->childidxs[factorpos+1] != i);

         SCIPdebugMessage("attempt expanding child %d at monomial %d factor %d\n", i, j, factorpos);

         SCIP_CALL( polynomialdataExpandMonomialFactor(blkmem, messagehdlr, polynomialdata, j, factorpos,
               (SCIP_EXPRDATA_POLYNOMIAL*)node->children[i]->data.data, childmap, maxexpansionexponent, &success) );

         if( !success )
         {
            removechild = FALSE;
            ++j;
         }
         else
            *havechange = TRUE;

         /* expansion may remove monomials[j], move a monomial from the end to position j, or add new monomials to the end of polynomialdata
          * we thus repeat with index j, if a factor was successfully expanded
          */
      }

      /* forget about child i, if it is not used anymore */
      if( removechild )
      {
         /* remove node from list of parents of child i */
         SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &node->children[i], node) );
         node->children[i] = NULL;
      }
   }

   /* simplify current polynomial again */
   polynomialdataMergeMonomials(blkmem, polynomialdata, eps, TRUE);

   BMSfreeBlockMemoryArrayNull(blkmem, &childmap, childmapsize);

   /* check which children are still in use */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &childinuse, node->nchildren) );
   BMSclearMemoryArray(childinuse, node->nchildren);  /*lint !e644*/
   for( i = 0; i < polynomialdata->nmonomials; ++i )
   {
      monomial = polynomialdata->monomials[i];
      assert(monomial != NULL);

      for( j = 0; j < monomial->nfactors; ++j )
      {
         assert(monomial->childidxs[j] >= 0);
         assert(monomial->childidxs[j] < node->nchildren);
         childinuse[monomial->childidxs[j]] = TRUE;
      }
   }

   /* free children that are not used in any monomial */
   for( i = 0; i < node->nchildren; ++i )
      if( node->children[i] != NULL && !childinuse[i] )
      {
         SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &node->children[i], node) );
         node->children[i] = NULL;
      }

   BMSfreeBlockMemoryArray(blkmem, &childinuse, node->nchildren);

   /* remove NULLs from children array */
   SCIP_CALL( exprgraphNodeRemovePolynomialNullChildren(blkmem, node) );

   /* if no children, then it's a constant polynomial -> change into EXPR_CONST */
   if( node->nchildren == 0 )
   {
      SCIP_Real val;

      /* if no children, then it should also have no monomials */
      assert(polynomialdata->nmonomials == 0);

      val = polynomialdata->constant;
      polynomialdataFree(blkmem, &polynomialdata);

      node->op = SCIP_EXPR_CONST;
      node->data.dbl = val;
      node->value = val;
   }

   /* if no factor in a monomial was replaced, the number of children should not have changed
    * but if we found duplicates in the children array, then it should be reduced, and we want to count this as a change too
    */
   *havechange |= (node->nchildren < orignchildren);  /*lint !e514*/

   node->simplified = TRUE;

   SCIPdebugMessage("-> %p = ", (void*)node);
   SCIPdebug( exprgraphPrintNodeExpression(node, messagehdlr, NULL, NULL, FALSE) );
   SCIPdebugPrintf("\n");

   return SCIP_OKAY;
}

/** creates an expression from a given node in an expression graph
 *
 *  Assembles mapping of variables from graph to tree.
 */
static
SCIP_RETCODE exprgraphNodeCreateExpr(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node which expression should be created */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to created expression */
   int*                  nexprvars,          /**< current number of variables in expression */
   int*                  varidx              /**< current mapping of variable indices from graph to expression */
   )
{
   SCIP_EXPR** childexprs;
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(expr != NULL);
   assert(nexprvars != NULL);
   assert(*nexprvars >= 0);
   assert(varidx != NULL);

   childexprs = NULL;
   if( node->nchildren > 0 )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childexprs, node->nchildren) );
      for( i = 0; i < node->nchildren; ++i )
      {
         SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, node->children[i], &childexprs[i], nexprvars, varidx) );  /*lint !e613*/
      }
   }

   switch( node->op )
   {
   case SCIP_EXPR_VARIDX:
   {
      /* check if the variable already has an index assigned in the expression tree
       * if not, create one and increase nexprvars
       */
      assert(node->data.intval >= 0);
      assert(node->data.intval < exprgraph->nvars);
      assert(varidx[node->data.intval] >= -1);
      assert(varidx[node->data.intval] < *nexprvars);
      if( varidx[node->data.intval] == -1 )
      {
         varidx[node->data.intval] = *nexprvars;
         ++*nexprvars;
      }

      SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, expr, SCIP_EXPR_VARIDX, varidx[node->data.intval]) );
      break;
   }

   case SCIP_EXPR_CONST:
   {
      SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, expr, node->op, node->data.dbl) );
      break;
   }

   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_SIGNPOWER:
   {
      assert(node->nchildren == 1);
      assert(childexprs != NULL);
      SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, expr, node->op, childexprs[0], node->data.dbl) );  /*lint !e613*/
      break;
   }

   case SCIP_EXPR_INTPOWER:
   {
      assert(node->nchildren == 1);
      assert(childexprs != NULL);
      SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, expr, node->op, childexprs[0], node->data.intval) );  /*lint !e613*/
      break;
   }

   case SCIP_EXPR_PLUS:
   case SCIP_EXPR_MINUS:
   case SCIP_EXPR_MUL:
   case SCIP_EXPR_DIV:
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   {
      assert(node->nchildren == 2);
      assert(childexprs != NULL);
      SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, expr, node->op, childexprs[0], childexprs[1]) );  /*lint !e613*/
      break;
   }

   case SCIP_EXPR_SQUARE:
   case SCIP_EXPR_SQRT:
   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   {
      assert(node->nchildren == 1);
      assert(childexprs != NULL);
      SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, expr, node->op, childexprs[0]) );  /*lint !e613*/
      break;
   }

   case SCIP_EXPR_SUM:
   case SCIP_EXPR_PRODUCT:
   {
      SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, expr, node->op, node->nchildren, childexprs) );
      break;
   }

   case SCIP_EXPR_LINEAR:
   {
      assert(node->data.data != NULL);

      SCIP_CALL( SCIPexprCreateLinear(exprgraph->blkmem, expr, node->nchildren, childexprs, (SCIP_Real*)node->data.data, ((SCIP_Real*)node->data.data)[node->nchildren]) );
      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* quaddata;

      quaddata = (SCIP_EXPRDATA_QUADRATIC*)node->data.data;
      assert(quaddata != NULL);

      SCIP_CALL( SCIPexprCreateQuadratic(exprgraph->blkmem, expr, node->nchildren, childexprs,
            quaddata->constant, quaddata->lincoefs, quaddata->nquadelems, quaddata->quadelems) );
      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;

      polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data;
      assert(polynomialdata != NULL);

      SCIP_CALL( SCIPexprCreatePolynomial(exprgraph->blkmem, expr, node->nchildren, childexprs,
            polynomialdata->nmonomials, polynomialdata->monomials, polynomialdata->constant, TRUE) );

      break;
   }

   case SCIP_EXPR_USER:
   {
      SCIP_EXPRDATA_USER* exprdata;
      SCIP_USEREXPRDATA* userdata;

      exprdata = (SCIP_EXPRDATA_USER*)node->data.data;
      assert(exprdata != NULL);

      if( exprdata->copydata != NULL )
      {
         SCIP_CALL( exprdata->copydata(exprgraph->blkmem, node->nchildren, exprdata->userdata, &userdata) );
      }
      else
         userdata = exprdata->userdata;

      SCIP_CALL( SCIPexprCreateUser(exprgraph->blkmem, expr, node->nchildren, childexprs,
         userdata, exprdata->evalcapability, exprdata->eval, exprdata->inteval, exprdata->curv, exprdata->prop, exprdata->estimate, exprdata->copydata, exprdata->freedata, exprdata->print) );

      break;
   }

   case SCIP_EXPR_LAST:
   case SCIP_EXPR_PARAM:
   {
      SCIPerrorMessage("expression operand %d not supported here\n", node->op);
      return SCIP_ERROR;
   }
   }

   BMSfreeBlockMemoryArrayNull(exprgraph->blkmem, &childexprs, node->nchildren);

   return SCIP_OKAY;
}

/** counts how often expression graph variables are used in a subtree of the expression graph
 *
 *  @note The function does not clear the array first, but only increases already existing counts.
 */
static
void exprgraphNodeGetVarsUsage(
   SCIP_EXPRGRAPHNODE*   node,               /**< root node of expression graph subtree */
   int*                  varsusage           /**< array where to count usage of variables, length must be at least the number of variables in the graph */
   )
{
   int i;

   assert(node != NULL);
   assert(varsusage != NULL);

   if( node->op == SCIP_EXPR_VARIDX )
   {
      ++varsusage[node->data.intval];
      return;
   }

   for( i = 0; i < node->nchildren; ++i )
      exprgraphNodeGetVarsUsage(node->children[i], varsusage);
}

/** checks whether a node can be put into a component when checking block separability of an expression
 *
 *  If a variable used by node is already in another component, components are merged and component number is updated.
 */
static
void exprgraphNodeCheckSeparabilityComponent(
   SCIP_EXPRGRAPHNODE*   node,               /**< node to which we assign a component */
   int*                  compnr,             /**< component number to assign, may be reduced if variables overlap */
   int                   nchildcomps,        /**< number of entries for which childcomps have been set already */
   int*                  childcomps,         /**< component numbers of children */
   int                   nvars,              /**< number of variables */
   int*                  varcomps            /**< component numbers of variables */
   )
{
   int varidx;
   int i;

   assert(node != NULL);
   assert(compnr != NULL);
   assert(*compnr >= 0);
   assert(childcomps != NULL);
   assert(varcomps != NULL);

   if( node->op != SCIP_EXPR_VARIDX )
   {
      for( i = 0; i < node->nchildren; ++i )
         exprgraphNodeCheckSeparabilityComponent(node->children[i], compnr, nchildcomps, childcomps, nvars, varcomps);
      return;
   }

   varidx = node->data.intval;
   assert(varidx >= 0);
   assert(varidx < nvars);

   if( varcomps[varidx] == -1 )
   {
      /* first time we get to this variable, so set it's component to compnr and we are done */
      varcomps[varidx] = *compnr;
      return;
   }

   if( varcomps[varidx] == *compnr )
   {
      /* variable is already in current component, that's also good and we are done */
      return;
   }

   /* variable is already in another component, so have to merge component compnr into that component
    * do this by updating varcomps and childcomps */
   for( i = 0; i < nvars; ++i )
      if( varcomps[i] == *compnr )
         varcomps[i] = varcomps[varidx];
   for( i = 0; i < nchildcomps; ++i )
      if( childcomps[i] == *compnr )
         childcomps[i] = varcomps[varidx];
   *compnr = varcomps[varidx];
}

/**@} */

/**@name Expression graph private methods */
/**@{ */

/** assert that expression graph has at least a given depth */
static
SCIP_RETCODE exprgraphEnsureDepth(
   SCIP_EXPRGRAPH*       exprgraph,          /**< buffer to store pointer to expression graph */
   int                   mindepth            /**< minimal depth that should be ensured */
   )
{
   int olddepth;

   assert(exprgraph != NULL);
   assert(exprgraph->blkmem != NULL);

   if( mindepth <= exprgraph->depth )
      return SCIP_OKAY;

   olddepth = exprgraph->depth;
   ensureBlockMemoryArraySize3(exprgraph->blkmem, &exprgraph->nodessize, &exprgraph->nnodes, &exprgraph->nodes, &exprgraph->depth, mindepth);
   assert(exprgraph->depth >= mindepth);

   /* initialize new array entries to 0 and NULL, resp. */
   BMSclearMemoryArray(&exprgraph->nodessize[olddepth], exprgraph->depth - olddepth);  /*lint !e866*/
   BMSclearMemoryArray(&exprgraph->nnodes[olddepth],    exprgraph->depth - olddepth);  /*lint !e866*/
   BMSclearMemoryArray(&exprgraph->nodes[olddepth],     exprgraph->depth - olddepth);  /*lint !e866*/

   return SCIP_OKAY;
}

/** remove a variable from the variables arrays, assuming that its node will be removed or converted next */
static
SCIP_RETCODE exprgraphRemoveVar(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   int                   varidx              /**< variable index */
   )
{
   SCIP_EXPRGRAPHNODE* varnode;
   void* var;

   assert(exprgraph != NULL);
   assert(varidx >= 0);
   assert(varidx < exprgraph->nvars);

   varnode = exprgraph->varnodes[varidx];
   assert(varnode->data.intval == varidx);

   var = exprgraph->vars[varidx];

   /* call varremove callback method, if set */
   if( exprgraph->exprgraphvarremove != NULL )
   {
      SCIP_CALL( exprgraph->exprgraphvarremove(exprgraph, exprgraph->userdata, var, varnode) );
   }

   /* remove variable from hashmap */
   SCIP_CALL( SCIPhashmapRemove(exprgraph->varidxs, var) );

   /* move last variable to position varidx and give it the new index */
   if( varidx < exprgraph->nvars-1 )
   {
      /* call callback method, if set */
      if( exprgraph->exprgraphvarchgidx != NULL )
      {
         SCIP_CALL( exprgraph->exprgraphvarchgidx(exprgraph, exprgraph->userdata, exprgraph->vars[exprgraph->nvars-1], exprgraph->varnodes[exprgraph->nvars-1], exprgraph->nvars-1, varidx) );
      }

      exprgraph->vars[varidx]      = exprgraph->vars[exprgraph->nvars-1];
      exprgraph->varbounds[varidx] = exprgraph->varbounds[exprgraph->nvars-1];
      exprgraph->varnodes[varidx]  = exprgraph->varnodes[exprgraph->nvars-1];
      exprgraph->varnodes[varidx]->data.intval = varidx;
      SCIP_CALL( SCIPhashmapSetImage(exprgraph->varidxs, exprgraph->vars[varidx], (void*)(size_t)(varidx)) );
   }
   --exprgraph->nvars;

   return SCIP_OKAY;
}

/** moves a node in an expression graph to a different depth
 *
 *  New depth must be larger than children depth.
 *  Moves parent nodes to higher depth, if needed.
 *  Variable nodes cannot be moved.
 */
static
SCIP_RETCODE exprgraphMoveNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< node that shall be moved */
   int                   newdepth            /**< new depth to which to move node */
   )
{
   int olddepth;
   int oldpos;
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->depth >= 0); /* node should be in graph */
   assert(newdepth >= 0);

   /* if already on aimed depth, then don't need to move */
   if( node->depth == newdepth )
      return SCIP_OKAY;

   SCIPdebugMessage("move node %p (%d,%d) to depth %d\n", (void*)node, node->depth, node->pos, newdepth);

#ifndef NDEBUG
   /* assert that children are at lower depth than new depth */
   for( i = 0; i < node->nchildren; ++i )
      assert(node->children[i]->depth < newdepth);
#endif

   /* move parents to higher depth, if needed */
   for( i = 0; i < node->nparents; ++i )
   {
      if( node->parents[i]->depth <= newdepth )
      {
         /* move parent to depth+1 */
         SCIP_CALL( exprgraphMoveNode(exprgraph, node->parents[i], newdepth+1) );
         assert(node->parents[i]->depth > newdepth);
      }
   }

   /* ensure that graph is deep enough */
   SCIP_CALL( exprgraphEnsureDepth(exprgraph, newdepth+1) );
   assert(exprgraph->depth > newdepth);

   olddepth = node->depth;
   oldpos   = node->pos;

   /* add node to new depth */
   ensureBlockMemoryArraySize(exprgraph->blkmem, &exprgraph->nodes[newdepth], &exprgraph->nodessize[newdepth], exprgraph->nnodes[newdepth]+1);  /*lint !e866*/
   node->depth = newdepth;
   node->pos   = exprgraph->nnodes[newdepth];
   exprgraph->nodes[newdepth][node->pos] = node;
   ++exprgraph->nnodes[newdepth];

   /* by moving the node to a new depth, the parents array in all its childrens may not be sorted anymore (parents order depends on depth) */
   for( i = 0; i < node->nchildren; ++i )
      node->children[i]->parentssorted = FALSE;

   /* move last node at previous depth to previous position, if it wasn't last */
   if( oldpos < exprgraph->nnodes[olddepth]-1 )
   {
      exprgraph->nodes[olddepth][oldpos] = exprgraph->nodes[olddepth][exprgraph->nnodes[olddepth]-1];
      exprgraph->nodes[olddepth][oldpos]->pos = oldpos;

      /* by moving the node to a new position, the parents array in all its children may not be sorted anymore (parents order depends on depth) */
      for( i = 0; i < exprgraph->nodes[olddepth][oldpos]->nchildren; ++i )
         exprgraph->nodes[olddepth][oldpos]->children[i]->parentssorted = FALSE;
   }
   --exprgraph->nnodes[olddepth];

   if( node->depth == 0 )
   {
      /* if at depth 0, then it need to be a node for either a constant or a variable */
      assert(node->op == SCIP_EXPR_CONST || node->op == SCIP_EXPR_VARIDX);
      if( node->op == SCIP_EXPR_CONST )
      {
         /* add node to constnodes array of exprgraph @todo should use SCIPsortedvecInsertPtr? */
         ensureBlockMemoryArraySize(exprgraph->blkmem, &exprgraph->constnodes, &exprgraph->constssize, exprgraph->nconsts + 1);
         exprgraph->constnodes[exprgraph->nconsts] = node;
         ++exprgraph->nconsts;
         exprgraph->constssorted = exprgraph->nconsts <= 1 || (exprgraph->constssorted && exprgraphConstNodeComp(exprgraph->constnodes[exprgraph->nconsts-2], node) < 0);
      }
      else
      {
         /* adding a variable by moving it from a higher depth seems awkward, how did the variable get there in the first place? */
         SCIPerrorMessage("cannot move variable nodes to depth 0\n");
         return SCIP_ERROR;
      }

      /* nodes at depth 0 always have curvature linear, even before any curvature check was running */
      node->curv = SCIP_EXPRCURV_LINEAR;
   }

   return SCIP_OKAY;
}

/** given a list of children, tries to find a common parent that represents a given operator with the same given data */
static
SCIP_RETCODE exprgraphFindParentByOperator(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   int                   nchildren,          /**< number of children */
   SCIP_EXPRGRAPHNODE**  children,           /**< children which parents to inspect */
   SCIP_EXPROP           op,                 /**< operator */
   SCIP_EXPROPDATA       opdata,             /**< operator data */
   SCIP_EXPR**           exprchildren,       /**< children of expression to consider when modifying (reordering) operator data, or NULL */
   SCIP_EXPRGRAPHNODE**  parent              /**< buffer to store parent node if any is found, or NULL if none found */
   )
{
   SCIP_EXPRGRAPHNODE** parentcands;
   int nparentcands;
   int parentcandssize;
   int i;
   int p;

   assert(exprgraph != NULL);
   assert(nchildren > 0);
   assert(children != NULL);
   assert(parent != NULL);

   *parent = NULL;

   /* create initial set of parent candidates as
    * all parents of first child that have the same operator type and the same number of children
    * additionally, some easy conditions for complex expression types:
    * if expression type is int/real/signpower, then compare also exponent,
    * if expression type is linear, then compare also constant part,
    * if expression type is quadratic, then compare also number of quadratic elements,
    * if expression type is polynomial, then compare also number of monmials and constant part
    */
   parentcandssize = children[0]->nparents;
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &parentcands, parentcandssize) );
   nparentcands = 0;
   for( p = 0; p < children[0]->nparents; ++p )
      if( children[0]->parents[p]->op == op &&
         children[0]->parents[p]->nchildren == nchildren &&
         (op != SCIP_EXPR_INTPOWER   || opdata.intval == children[0]->parents[p]->data.intval) &&
         (op != SCIP_EXPR_REALPOWER  || opdata.dbl == children[0]->parents[p]->data.dbl) &&  /*lint !e777*/
         (op != SCIP_EXPR_SIGNPOWER  || opdata.dbl == children[0]->parents[p]->data.dbl) &&  /*lint !e777*/
         (op != SCIP_EXPR_LINEAR     || ((SCIP_Real*)opdata.data)[nchildren] == ((SCIP_Real*)children[0]->parents[p]->data.data)[nchildren]) &&  /*lint !e777*/
         (op != SCIP_EXPR_QUADRATIC  || ((SCIP_EXPRDATA_QUADRATIC*)opdata.data)->nquadelems == ((SCIP_EXPRDATA_QUADRATIC*)children[0]->parents[p]->data.data)->nquadelems) &&
         (op != SCIP_EXPR_QUADRATIC  || ((SCIP_EXPRDATA_QUADRATIC*)opdata.data)->constant == ((SCIP_EXPRDATA_QUADRATIC*)children[0]->parents[p]->data.data)->constant) &&  /*lint !e777*/
         (op != SCIP_EXPR_POLYNOMIAL || ((SCIP_EXPRDATA_POLYNOMIAL*)opdata.data)->nmonomials == ((SCIP_EXPRDATA_POLYNOMIAL*)children[0]->parents[p]->data.data)->nmonomials) &&
         (op != SCIP_EXPR_POLYNOMIAL || ((SCIP_EXPRDATA_POLYNOMIAL*)opdata.data)->constant == ((SCIP_EXPRDATA_POLYNOMIAL*)children[0]->parents[p]->data.data)->constant)  /*lint !e777*/
         )
      {
         parentcands[nparentcands++] = children[0]->parents[p];  /*lint !e644*/
      }

   /* for all remaining children, remove parent candidates, that are not in their list of parents */
   for( i = 1; i < nchildren && nparentcands > 0; ++i )
   {
      p = 0;
      while( p < nparentcands )
      {
         /* if parentcands[p] is a parent of childnodes[i], then move last parent candidate to position p,
          * otherwise keep candidate and check next one
          */
         if( !exprgraphNodeIsParent(children[i], parentcands[p]) )
         {
            parentcands[p] = parentcands[nparentcands-1];
            --nparentcands;
         }
         else
            ++p;
      }
   }

   SCIPdebugMessage("check %d parent candidates for expr with operator %d and %d children\n", nparentcands, op, nchildren);

   if( nparentcands == 0 )
   {
      BMSfreeBlockMemoryArray(exprgraph->blkmem, &parentcands, children[0]->nparents);
      return SCIP_OKAY;
   }

   /* at this point, all parents in parentcands have the nodes in children as children and are of the same operator type
    * check if there is also one which corresponds to same expression and store that one in *parent
    */
   switch( op )
   {
      /* commutative operands with no data */
   case SCIP_EXPR_PLUS   :
   case SCIP_EXPR_MUL    :
   case SCIP_EXPR_MIN    :
   case SCIP_EXPR_MAX    :
   case SCIP_EXPR_SUM    :
   case SCIP_EXPR_PRODUCT:
   case SCIP_EXPR_SQUARE :
   case SCIP_EXPR_SQRT   :
   case SCIP_EXPR_EXP    :
   case SCIP_EXPR_LOG    :
   case SCIP_EXPR_SIN    :
   case SCIP_EXPR_COS    :
   case SCIP_EXPR_TAN    :
      /* case SCIP_EXPR_ERF    : */
      /* case SCIP_EXPR_ERFI   : */
   case SCIP_EXPR_ABS    :
   case SCIP_EXPR_SIGN   :
   {
      /* sort childnodes, if needed for later */
      if( nchildren > 2 )
         SCIPsortPtr((void**)children, exprgraphnodecomp, nchildren);
      for( p = 0; p < nparentcands; ++p )
      {
         assert(parentcands[p]->op        == op);        /* that was the first  criterium for adding a node to parentcands */
         assert(parentcands[p]->nchildren == nchildren); /* that was the second criterium for adding a node to parentcands */

         if( nchildren == 1 )
         {
            assert(parentcands[p]->children[0] == children[0]);
            /* same operand, same child, so same expression */
            *parent = parentcands[p];
            break;
         }
         else if( nchildren == 2 )
         {
            /* We know that every node in children is also a child of parentcands[p].
             * However, if there are duplicates in children, then it can happen that not every child of parentcands[p] is also one of the children.
             * So only if children equals parentcands[p]->children in some permutation, the expressions are the same.
             */
            if( (parentcands[p]->children[0] == children[0] && parentcands[p]->children[1] == children[1]) ||
               ( parentcands[p]->children[0] == children[1] && parentcands[p]->children[1] == children[0]) )
            {
               *parent = parentcands[p];
               break;
            }
         }
         else
         {
            /* as in the case for two nodes, we need to check whether parentcands[p]->children and children are equal up to permutation */

            /* sort children of parent candidate */
            SCIPsortPtr((void**)parentcands[p]->children, exprgraphnodecomp, nchildren);

            /* check if childnodes and parentcands[p]->children are the same */
            for( i = 0; i < nchildren; ++i )
               if( children[i] != parentcands[p]->children[i] )
                  break;
            if( i == nchildren )
            {
               /* yeah, found an exact match */
               *parent = parentcands[p];
               break;
            }
         }
      }

      break;
   }

   /* non-commutative operands with two children */
   case SCIP_EXPR_MINUS    :
   case SCIP_EXPR_DIV      :
   {
      for( p = 0; p < nparentcands; ++p )
      {
         assert(parentcands[p]->op == op); /* that was the first  criterium for adding a node to parentcands */
         assert(parentcands[p]->nchildren == 2); /* that was the second criterium for adding a node to parentcands */
         /* order of operands matters, so check if childnodes have same order as children of parent candidate (and are the same nodes too) */
         if( parentcands[p]->children[0] == children[0] && parentcands[p]->children[1] == children[1] )
         {
            /* yeah, found one */
            *parent = parentcands[p];
            break;
         }
      }

      break;
   }

   /* operands with one child and data */
   case SCIP_EXPR_INTPOWER:
   {
      assert(parentcands[0]->op == op); /* that was the first  criterium for adding a node to parentcands */
      assert(parentcands[0]->nchildren == 1); /* that was the second criterium for adding a node to parentcands */
      assert(parentcands[0]->children[0] == children[0]); /* that's what exprgraphNodeIsParent should have ensured */
      assert(parentcands[0]->data.intval == opdata.intval); /* that was another criterium for adding a node to parentcands */

      /* yeah, have one with same exponent */
      *parent = parentcands[0];

      break;
   }

   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_SIGNPOWER:
   {
      assert(parentcands[0]->op == op); /* that was the first  criterium for adding a node to parentcands */
      assert(parentcands[0]->nchildren == 1); /* that was the second criterium for adding a node to parentcands */
      assert(parentcands[0]->children[0] == children[0]); /* that's what exprgraphNodeIsParent should have ensured */
      assert(parentcands[0]->data.dbl == opdata.dbl); /* that was another criterium for adding a node to parentcands */  /*lint !e777*/

      /* yeah, have one with same exponent */
      *parent = parentcands[0];

      break;
   }

   /* commutative operands with n children and data */
   case SCIP_EXPR_LINEAR:
   {
      SCIP_Real* exprcoef;
      SCIP_Real* candcoef;

      exprcoef = (SCIP_Real*)opdata.data;
      /* sort childnodes, take care that children in expression are sorted the same way if given (so we don't mess up assignment of coefficients) */
      if( exprchildren != NULL )
         SCIPsortPtrPtrReal((void**)children, (void**)exprchildren, exprcoef, exprgraphnodecomp, nchildren);
      else
         SCIPsortPtrReal((void**)children, exprcoef, exprgraphnodecomp, nchildren);
      for( p = 0; p < nparentcands; ++p )
      {
         assert(parentcands[p]->op        == op);        /* that was the first  criterium for adding a node to parentcands */
         assert(parentcands[p]->nchildren == nchildren); /* that was the second criterium for adding a node to parentcands */

         candcoef = (SCIP_Real*)parentcands[p]->data.data;
         assert(exprcoef[nchildren] == candcoef[nchildren]); /* that was a criterium for adding a node to parentcands */  /*lint !e777*/

         /* sort children of parent candidate */
         SCIPsortPtrReal((void**)parentcands[p]->children, candcoef, exprgraphnodecomp, nchildren);

         /* check if children and coefficients in parent candidate and expression are the same */
         for( i = 0; i < nchildren; ++i )
         {
            if( children[i] != parentcands[p]->children[i] )
               break;
            if( exprcoef[i] != candcoef[i] )  /*lint !e777*/
               break;
         }
         if( i < nchildren )
            continue;

         /* yeah, found an exact match */
         *parent = parentcands[p];
         break;
      }

      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* exprdata;
      SCIP_Real* exprlincoef;
      SCIP_Real* candlincoef;
      SCIP_EXPRDATA_QUADRATIC* canddata;
      int* perm;
      int* invperm;

      exprdata = (SCIP_EXPRDATA_QUADRATIC*)opdata.data;
      exprlincoef = exprdata->lincoefs;

      /* sort children in expr and parentcands and update indices in quadelems accordingly, then sort quadelems again and compare */

      /* sort expr->children and childnodes and store inverse permutation in invperm */
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &invperm, nchildren) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &perm,    nchildren) );
      for( i = 0; i < nchildren; ++i )
         invperm[i] = i;  /*lint !e644*/

      if( exprlincoef != NULL )
         if( exprchildren != NULL )
            SCIPsortPtrPtrRealInt((void**)children, (void**)exprchildren, exprlincoef, invperm, exprgraphnodecomp, nchildren);
         else
            SCIPsortPtrRealInt((void**)children, exprlincoef, invperm, exprgraphnodecomp, nchildren);
      else
         if( exprchildren != NULL )
            SCIPsortPtrPtrInt((void**)children, (void**)exprchildren, invperm, exprgraphnodecomp, nchildren);
         else
            SCIPsortPtrInt((void**)children, invperm, exprgraphnodecomp, nchildren);

      /* compute permutation from its inverse */
      for( i = 0; i < nchildren; ++i )
         perm[invperm[i]] = i;  /*lint !e644*/

      /* apply permuation to exprdata->quadelems and sort again */
      for( i = 0; i < exprdata->nquadelems; ++i )
      {
         exprdata->quadelems[i].idx1 = perm[exprdata->quadelems[i].idx1];
         exprdata->quadelems[i].idx2 = perm[exprdata->quadelems[i].idx2];
         if( exprdata->quadelems[i].idx1 > exprdata->quadelems[i].idx2 )
         {
            int tmp;
            tmp = exprdata->quadelems[i].idx1;
            exprdata->quadelems[i].idx1 = exprdata->quadelems[i].idx2;
            exprdata->quadelems[i].idx2 = tmp;
         }
      }
      SCIPquadelemSort(exprdata->quadelems, exprdata->nquadelems);
      exprdata->sorted = TRUE;

      for( p = 0; p < nparentcands; ++p )
      {
         assert(parentcands[p]->op        == op);        /* that was the first  criterium for adding a node to parentcands */
         assert(parentcands[p]->nchildren == nchildren); /* that was the second criterium for adding a node to parentcands */

         canddata = (SCIP_EXPRDATA_QUADRATIC*)parentcands[p]->data.data;
         candlincoef = canddata->lincoefs;
         assert(canddata->nquadelems == exprdata->nquadelems); /* that was a criterium for adding a node to parentcands */
         assert(canddata->constant   == exprdata->constant);   /* that was a criterium for adding a node to parentcands */  /*lint !e777*/

         /* sort parentcands[p]->children and store inverse permutation in invperm */
         for( i = 0; i < nchildren; ++i )
            invperm[i] = i;

         if( candlincoef != NULL )
            SCIPsortPtrRealInt((void**)parentcands[p]->children, candlincoef, invperm, exprgraphnodecomp, parentcands[p]->nchildren);
         else
            SCIPsortPtrInt((void**)parentcands[p]->children, invperm, exprgraphnodecomp, nchildren);

         /* compute permutation from its inverse */
         for( i = 0; i < nchildren; ++i )
            perm[invperm[i]] = i;

         /* apply permutation to canddata->quadelems */
         for( i = 0; i < canddata->nquadelems; ++i )
         {
            canddata->quadelems[i].idx1 = perm[canddata->quadelems[i].idx1];
            canddata->quadelems[i].idx2 = perm[canddata->quadelems[i].idx2];
            if( canddata->quadelems[i].idx1 > canddata->quadelems[i].idx2 )
            {
               int tmp;
               tmp = canddata->quadelems[i].idx1;
               canddata->quadelems[i].idx1 = canddata->quadelems[i].idx2;
               canddata->quadelems[i].idx2 = tmp;
            }
         }
         SCIPquadelemSort(canddata->quadelems, canddata->nquadelems);
         canddata->sorted = TRUE;

         /* check if children and linear coefficients in parent candidate and expression are the same */
         for( i = 0; i < nchildren; ++i )
         {
            if( children[i] != parentcands[p]->children[i] )
               break;
            if( (exprlincoef == NULL ? 0.0 : exprlincoef[i]) != (candlincoef == NULL ? 0.0 : candlincoef[i]) )  /*lint !e777*/
               break;
         }
         if( i < nchildren )
            continue;

         assert(exprdata->nquadelems == canddata->nquadelems);
         for( i = 0; i < exprdata->nquadelems; ++i )
         {
            if( exprdata->quadelems[i].idx1 != canddata->quadelems[i].idx1 ||
               exprdata->quadelems[i].idx2 != canddata->quadelems[i].idx2 ||
               exprdata->quadelems[i].coef != canddata->quadelems[i].coef )  /*lint !e777*/
               break;
         }
         if( i == exprdata->nquadelems )
         {
            /* yeah, parentcands[p] is same quadratic expression as expr */
            *parent = parentcands[p];
            break;
         }
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &perm,    nchildren);
      BMSfreeBlockMemoryArray(exprgraph->blkmem, &invperm, nchildren);

      break;
   }

   /* @todo in one GlobalLib instance, two polynoms differ only in the sign of all coefficients, it would be nice to recognize this somehow */
   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* exprdata;
      SCIP_EXPRDATA_POLYNOMIAL* canddata;
      int* perm;
      int* invperm;

      exprdata = (SCIP_EXPRDATA_POLYNOMIAL*)opdata.data;

      /* sort children in expr and parentcands and update child indices in polynomialdata, then sort monomials again and compare */

      /* sort exprchildren and childnodes and store inverse permutation in invperm */
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &invperm, nchildren) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &perm,    nchildren) );
      for( i = 0; i < nchildren; ++i )
         invperm[i] = i;  /*lint !e644*/

      if( exprchildren != NULL )
         SCIPsortPtrPtrInt((void**)children, (void**)exprchildren, invperm, exprgraphnodecomp, nchildren);
      else
         SCIPsortPtrInt((void**)children, invperm, exprgraphnodecomp, nchildren);

      /* compute permutation from its inverse */
      for( i = 0; i < nchildren; ++i )
         perm[invperm[i]] = i;  /*lint !e644*/

      /* apply permutation to exprdata and sort again */
      polynomialdataApplyChildmap(exprdata, perm);
      polynomialdataSortMonomials(exprdata);

      for( p = 0; p < nparentcands; ++p )
      {
         assert(parentcands[p]->op        == op);        /* that was the first  criterium for adding a node to parentcands */
         assert(parentcands[p]->nchildren == nchildren); /* that was the second criterium for adding a node to parentcands */

         canddata = (SCIP_EXPRDATA_POLYNOMIAL*)parentcands[p]->data.data;
         assert(canddata->nmonomials == exprdata->nmonomials); /* that was a criterium for adding a node to parentcands */
         assert(canddata->constant   == exprdata->constant);   /* that was a criterium for adding a node to parentcands */  /*lint !e777*/

         /* sort parentcands[p]->children and store inverse permutation in invperm */
         for( i = 0; i < nchildren; ++i )
            invperm[i] = i;

         SCIPsortPtrInt((void**)parentcands[p]->children, invperm, exprgraphnodecomp, nchildren);

         /* compute permutation from its inverse */
         for( i = 0; i < nchildren; ++i )
            perm[invperm[i]] = i;

         /* apply permutation to canddata and sort again */
         polynomialdataApplyChildmap(canddata, perm);
         polynomialdataSortMonomials(canddata);

         /* check if children are equal */
         for( i = 0; i < nchildren; ++i )
            if( children[i] != parentcands[p]->children[i] )
               break;
         if( i < nchildren )
            continue;

         /* check if monomials are equal */
         for( i = 0; i < exprdata->nmonomials; ++i )
            if( !SCIPexprAreMonomialsEqual(exprdata->monomials[i], canddata->monomials[i], 0.0) )
               break;
         if( i == exprdata->nmonomials )
         {
            /* yeah, parentcands[p] is same polynomial expression as expr */
            *parent = parentcands[p];
            break;
         }
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &perm,    nchildren);
      BMSfreeBlockMemoryArray(exprgraph->blkmem, &invperm, nchildren);

      break;
   }

   case SCIP_EXPR_USER:
   {
      /* @todo need comparison function on user data to decide whether a parent candidate fits */
      break;
   }

   case SCIP_EXPR_VARIDX:
   case SCIP_EXPR_PARAM:
   case SCIP_EXPR_CONST:
   case SCIP_EXPR_LAST:
      SCIPerrorMessage("expression operand %d unexpected here\n", op);
      return SCIP_ERROR;
   }

   BMSfreeBlockMemoryArray(exprgraph->blkmem, &parentcands, parentcandssize);

   return SCIP_OKAY;
}

/** adds an expression into an expression graph
 *
 *  Enables corresponding nodes.
 */
static
SCIP_RETCODE exprgraphAddExpr(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPR*            expr,               /**< expression to add */
   void**                vars,               /**< variables corresponding to VARIDX expressions */
   SCIP_Real*            params,             /**< parameter values */
   SCIP_EXPRGRAPHNODE**  exprnode,           /**< buffer to store expression graph node corresponding to root of this expression */
   SCIP_Bool*            exprnodeisnew       /**< buffer to indicate whether the node in *exprnode has been newly created for this expression (otherwise, expression was already in graph) */
   )
{
   SCIP_EXPRGRAPHNODE** childnodes;
   SCIP_Bool childisnew;
   SCIP_Bool nochildisnew;
   SCIP_EXPROPDATA opdata;
   int i;

   assert(exprgraph != NULL);
   assert(expr != NULL);
   assert(exprnode != NULL);
   assert(exprnodeisnew != NULL);

   if( expr->op == SCIP_EXPR_VARIDX )
   {
      /* find node corresponding to variable and add if not existing yet */
      assert(expr->nchildren == 0);

      SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, &vars[expr->data.intval], exprnode) );
      assert(*exprnode != NULL);
      assert((*exprnode)->op == SCIP_EXPR_VARIDX);
      assert((*exprnode)->data.intval >= 0);
      assert((*exprnode)->data.intval < exprgraph->nvars);
      assert(exprgraph->vars[(*exprnode)->data.intval] == vars[expr->data.intval]);

      *exprnodeisnew = (*exprnode)->nuses == 0 && (*exprnode)->nparents == 0;

      return SCIP_OKAY;
   }

   if( expr->op == SCIP_EXPR_CONST )
   {
      /* find node corresponding to constant and add if not existing yet */
      assert(expr->nchildren == 0);

      SCIP_CALL( SCIPexprgraphAddConst(exprgraph, expr->data.dbl, exprnode) );
      assert(*exprnode != NULL);
      assert((*exprnode)->op == SCIP_EXPR_CONST);
      assert((*exprnode)->data.dbl == expr->data.dbl);  /*lint !e777*/

      *exprnodeisnew = (*exprnode)->nuses == 0 && (*exprnode)->nparents == 0;

      return SCIP_OKAY;
   }

   if( expr->op == SCIP_EXPR_PARAM )
   {
      /* find node corresponding to constant corresponding to parameter and add if not existing yet */
      assert(expr->nchildren == 0);
      assert(params != NULL);

      SCIP_CALL( SCIPexprgraphAddConst(exprgraph, params[expr->data.intval], exprnode) );
      assert(*exprnode != NULL);
      assert((*exprnode)->op == SCIP_EXPR_CONST);
      assert((*exprnode)->data.dbl == params[expr->data.intval]);  /*lint !e777*/

      *exprnodeisnew = (*exprnode)->nuses == 0 && (*exprnode)->nparents == 0;

      return SCIP_OKAY;
   }

   /* expression should be variable or constant or have children */
   assert(expr->nchildren > 0);

   /* add children expressions into expression graph
    * check if we can find a common parent
    */
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childnodes, expr->nchildren) );
   nochildisnew = TRUE;
   for( i = 0; i < expr->nchildren; ++i )
   {
      SCIP_CALL( exprgraphAddExpr(exprgraph, expr->children[i], vars, params, &childnodes[i], &childisnew) );  /*lint !e644*/
      assert(childnodes[i] != NULL);
      nochildisnew &= !childisnew;  /*lint !e514*/
   }

   /* if all children were known already, check if there is also already a node for the expression that we aim to add */
   if( nochildisnew )
   {
      SCIP_CALL( exprgraphFindParentByOperator(exprgraph, expr->nchildren, childnodes, expr->op, expr->data, expr->children, exprnode) );

      if( *exprnode != NULL )
      {
         /* node already existing, make sure it is enabled */
         (*exprnode)->enabled = TRUE;
         *exprnodeisnew = FALSE;

         /* SCIPdebugMessage("reused node %p (%d,%d) for expr ", (void*)*exprnode, (*exprnode)->depth, (*exprnode)->pos);
          * SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
          * SCIPdebugPrintf("\n");
          */

         BMSfreeBlockMemoryArray(exprgraph->blkmem, &childnodes, expr->nchildren);
         return SCIP_OKAY;
      }
   }

   SCIPdebugMessage("add expr with operator %d and %d children\n", expr->op, expr->nchildren);

   /* copy expression data */
   if( exprOpTable[expr->op].copydata != NULL )
   {
      SCIP_CALL( exprOpTable[expr->op].copydata(exprgraph->blkmem, expr->nchildren, expr->data, &opdata) );
   }
   else
   {
      opdata = expr->data;
   }

   SCIP_CALL( exprgraphCreateNode(exprgraph->blkmem, exprnode, expr->op, opdata) );
   SCIP_CALL( SCIPexprgraphAddNode(exprgraph, *exprnode, -1, expr->nchildren, childnodes) );
   *exprnodeisnew = TRUE;

   BMSfreeBlockMemoryArray(exprgraph->blkmem, &childnodes, expr->nchildren);

   /* SCIPdebugMessage("created new node %p (%d,%d) for expr ", (void*)*exprnode, (*exprnode)->depth, (*exprnode)->pos);
    * SCIPdebug( SCIPexprPrint(expr, messagehdlr, NULL, NULL, NULL, NULL) );
    * SCIPdebugPrintf("\n");
    */

   return SCIP_OKAY;
}

/** sets bounds in variable nodes to those stored in exprgraph's varbounds array */
static
void exprgraphUpdateVarNodeBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Bool*            clearreverseprop,   /**< flag to set if we had reset bound tightenings from reverse propagation */
   SCIP_Bool*            boundchanged        /**< buffer to store whether a variables bound has changes, compared to those stored in nodes */
   )
{
   SCIP_EXPRGRAPHNODE* node;
   int i;
   int p;

   assert(exprgraph != NULL);
   assert(clearreverseprop != NULL);
   assert(boundchanged != NULL);

   *boundchanged = FALSE;
   for( i = 0; i < exprgraph->nvars; ++i )
   {
      node = exprgraph->varnodes[i];

      if( node->bounds.inf == exprgraph->varbounds[i].inf &&  /*lint !e777*/
         +node->bounds.sup == exprgraph->varbounds[i].sup )   /*lint !e777*/
      {
         node->boundstatus = SCIP_EXPRBOUNDSTATUS_VALID;
         continue;
      }

      if( exprgraph->varbounds[i].inf > exprgraph->varbounds[i].sup )
      {
         /* hmm, may happen due to numerics, let's be conservative and relax bounds to something that seems reasonable */
         SCIP_Real tmp;

         tmp = exprgraph->varbounds[i].inf;
         exprgraph->varbounds[i].inf = MIN(tmp, exprgraph->varbounds[i].sup);
         exprgraph->varbounds[i].sup = MAX(tmp, exprgraph->varbounds[i].sup);
      }

      if( exprgraph->varbounds[i].inf < node->bounds.inf ||
         +exprgraph->varbounds[i].sup > node->bounds.sup )
      {
         for( p = 0; p < node->nparents; ++p )
            node->parents[p]->boundstatus = SCIP_EXPRBOUNDSTATUS_CHILDRELAXED;

         node->bounds = exprgraph->varbounds[i];
         SCIPdebugMessage("registered relaxed bound [%g,%g] of var %d for propagation\n", node->bounds.inf, node->bounds.sup, i);

         *boundchanged = TRUE;

         /* if a childs bounds are relaxed, then a previous reverse propagation may be invalid, so we have to clear its remainings */
         *clearreverseprop = TRUE;
      }
      else if( isLbBetter(1e-9, exprgraph->varbounds[i].inf, node->bounds.inf, node->bounds.sup) ||
         (     isUbBetter(1e-9, exprgraph->varbounds[i].sup, node->bounds.inf, node->bounds.sup)) )
      {
         for( p = 0; p < node->nparents; ++p )
            node->parents[p]->boundstatus |= SCIP_EXPRBOUNDSTATUS_CHILDTIGHTENED;

         node->bounds = exprgraph->varbounds[i];
         SCIPdebugMessage("registered tightened bound [%g,%g] of var %d for propagation\n", node->bounds.inf, node->bounds.sup, i);

         *boundchanged = TRUE;
      }
      else
      {
         node->bounds = exprgraph->varbounds[i];
         SCIPdebugMessage("registered slightly tightened bound [%g,%g] of var %d for propagation\n", node->bounds.inf, node->bounds.sup, i);
      }

      node->boundstatus = SCIP_EXPRBOUNDSTATUS_VALID;
   }
}

/**@} */

/**@name Expression graph node methods */
/**@{ */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPexprgraphCaptureNode
#undef SCIPexprgraphIsNodeEnabled
#undef SCIPexprgraphGetNodeNChildren
#undef SCIPexprgraphGetNodeChildren
#undef SCIPexprgraphGetNodeNParents
#undef SCIPexprgraphGetNodeParents
#undef SCIPexprgraphGetNodeDepth
#undef SCIPexprgraphGetNodePosition
#undef SCIPexprgraphGetNodeOperator
#undef SCIPexprgraphGetNodeOperatorIndex
#undef SCIPexprgraphGetNodeOperatorReal
#undef SCIPexprgraphGetNodeVar
#undef SCIPexprgraphGetNodeRealPowerExponent
#undef SCIPexprgraphGetNodeIntPowerExponent
#undef SCIPexprgraphGetNodeSignPowerExponent
#undef SCIPexprgraphGetNodeLinearCoefs
#undef SCIPexprgraphGetNodeLinearConstant
#undef SCIPexprgraphGetNodeQuadraticConstant
#undef SCIPexprgraphGetNodeQuadraticLinearCoefs
#undef SCIPexprgraphGetNodeQuadraticQuadElements
#undef SCIPexprgraphGetNodeQuadraticNQuadElements
#undef SCIPexprgraphGetNodePolynomialMonomials
#undef SCIPexprgraphGetNodePolynomialNMonomials
#undef SCIPexprgraphGetNodePolynomialConstant
#undef SCIPexprgraphGetNodeUserData
#undef SCIPexprgraphHasNodeUserEstimator
#undef SCIPexprgraphGetNodeBounds
#undef SCIPexprgraphGetNodeVal
#undef SCIPexprgraphGetNodeCurvature

/** captures node, i.e., increases number of uses */
void SCIPexprgraphCaptureNode(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node to capture */
   )
{
   assert(node->nuses >= 0);

   SCIPdebugMessage("capture node %p\n", (void*)node);

   ++node->nuses;
}

/** returns whether a node is currently enabled */
SCIP_Bool SCIPexprgraphIsNodeEnabled(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node to enable */
   )
{
   assert(node != NULL);

   return node->enabled;
}

/** gets number of children of a node in an expression graph */
int SCIPexprgraphGetNodeNChildren(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->nchildren;
}

/** gets children of a node in an expression graph */
SCIP_EXPRGRAPHNODE** SCIPexprgraphGetNodeChildren(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->children;
}

/** gets number of parents of a node in an expression graph */
int SCIPexprgraphGetNodeNParents(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->nparents;
}

/** gets parents of a node in an expression graph */
SCIP_EXPRGRAPHNODE** SCIPexprgraphGetNodeParents(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->parents;
}

/** gets depth of node in expression graph */
int SCIPexprgraphGetNodeDepth(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->depth;
}

/** gets position of node in expression graph at its depth level */
int SCIPexprgraphGetNodePosition(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->pos;
}

/** gets operator of a node in an expression graph */
SCIP_EXPROP SCIPexprgraphGetNodeOperator(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->op;
}

/** gives index belonging to a SCIP_EXPR_VARIDX or SCIP_EXPR_PARAM operand */
int SCIPexprgraphGetNodeOperatorIndex(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_VARIDX || node->op == SCIP_EXPR_PARAM);

   return node->data.intval;
}

/** gives real belonging to a SCIP_EXPR_CONST operand */
SCIP_Real SCIPexprgraphGetNodeOperatorReal(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_CONST);

   return node->data.dbl;
}

/** gives variable belonging to a SCIP_EXPR_VARIDX expression */
void* SCIPexprgraphGetNodeVar(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_VARIDX);
   assert(node->data.intval >= 0);
   assert(node->data.intval < exprgraph->nvars);

   return exprgraph->vars[node->data.intval];
}

/** gives exponent belonging to a SCIP_EXPR_REALPOWER expression */
SCIP_Real SCIPexprgraphGetNodeRealPowerExponent(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_REALPOWER);

   return node->data.dbl;
}

/** gives exponent belonging to a SCIP_EXPR_INTPOWER expression */
int SCIPexprgraphGetNodeIntPowerExponent(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_INTPOWER);

   return node->data.intval;
}

/** gives exponent belonging to a SCIP_EXPR_SIGNPOWER expression */
SCIP_Real SCIPexprgraphGetNodeSignPowerExponent(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_SIGNPOWER);

   return node->data.dbl;
}

/** gives linear coefficients belonging to a SCIP_EXPR_LINEAR expression */
SCIP_Real* SCIPexprgraphGetNodeLinearCoefs(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_LINEAR);

   return (SCIP_Real*)node->data.data;
}

/** gives constant belonging to a SCIP_EXPR_LINEAR expression  */
SCIP_Real SCIPexprgraphGetNodeLinearConstant(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_LINEAR);
   assert(node->data.data != NULL);

   return ((SCIP_Real*)node->data.data)[node->nchildren];
}

/** gives constant belonging to a SCIP_EXPR_QUADRATIC expression */
SCIP_Real SCIPexprgraphGetNodeQuadraticConstant(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_QUADRATIC);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_QUADRATIC*)node->data.data)->constant;
}

/** gives linear coefficients belonging to a SCIP_EXPR_QUADRATIC expression, or NULL if all coefficients are 0.0 */
SCIP_Real* SCIPexprgraphGetNodeQuadraticLinearCoefs(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_QUADRATIC);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_QUADRATIC*)node->data.data)->lincoefs;
}

/** gives quadratic elements belonging to a SCIP_EXPR_QUADRATIC expression */
SCIP_QUADELEM* SCIPexprgraphGetNodeQuadraticQuadElements(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_QUADRATIC);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_QUADRATIC*)node->data.data)->quadelems;
}

/** gives number of quadratic elements belonging to a SCIP_EXPR_QUADRATIC expression */
int SCIPexprgraphGetNodeQuadraticNQuadElements(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_QUADRATIC);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_QUADRATIC*)node->data.data)->nquadelems;
}

/** gives the monomials belonging to a SCIP_EXPR_POLYNOMIAL expression */
SCIP_EXPRDATA_MONOMIAL** SCIPexprgraphGetNodePolynomialMonomials(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_POLYNOMIAL);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data)->monomials;
}

/** gives the number of monomials belonging to a SCIP_EXPR_POLYNOMIAL expression */
int SCIPexprgraphGetNodePolynomialNMonomials(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_POLYNOMIAL);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data)->nmonomials;
}

/** gives the constant belonging to a SCIP_EXPR_POLYNOMIAL expression */
SCIP_Real SCIPexprgraphGetNodePolynomialConstant(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_POLYNOMIAL);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data)->constant;
}

/** gives the curvature of a single monomial belonging to a SCIP_EXPR_POLYNOMIAL expression
 *
 *  Assumes that curvature of children and bounds of children and node itself are valid.
 */
SCIP_RETCODE SCIPexprgraphGetNodePolynomialMonomialCurvature(
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   int                   monomialidx,        /**< index of monomial */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_EXPRCURV*        curv                /**< buffer to store monomial curvature */
   )
{
   SCIP_EXPRDATA_MONOMIAL* monomial;
   SCIP_INTERVAL  childboundsstatic[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_EXPRCURV  childcurvstatic[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_INTERVAL* childbounds = NULL;
   SCIP_EXPRCURV* childcurv = NULL;
   SCIP_EXPRGRAPHNODE* child;
   SCIP_RETCODE retcode = SCIP_OKAY;
   int i;

   assert(node != NULL);
   assert(node->depth >= 0); /* node should be in graph */
   assert(node->pos >= 0);   /* node should be in graph */
   assert(node->enabled);    /* node should be enabled, otherwise we may not have uptodate bounds and curvatures in children */
   assert(node->boundstatus == SCIP_EXPRBOUNDSTATUS_VALID);  /* we assume node bounds to be valid */
   assert(node->op == SCIP_EXPR_POLYNOMIAL);
   assert(node->data.data != NULL);
   assert(monomialidx >= 0);
   assert(monomialidx < ((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data)->nmonomials);
   assert(curv != NULL);

   if( SCIPintervalIsEmpty(infinity, node->bounds) )
   {
      *curv = SCIP_EXPRCURV_LINEAR;
      return SCIP_OKAY;
   }

   monomial = ((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data)->monomials[monomialidx];
   assert(monomial != NULL);

   /* if many children, get large enough memory to store children bounds */
   if( monomial->nfactors > SCIP_EXPRESSION_MAXCHILDEST )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&childbounds, monomial->nfactors) );
      SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&childcurv, monomial->nfactors), TERMINATE );
   }
   else
   {
      childbounds = childboundsstatic;
      childcurv   = childcurvstatic;
   }

   /* assemble bounds and curvature of children */
   for( i = 0; i < monomial->nfactors; ++i )
   {
      child = node->children[monomial->childidxs[i]];
      assert(child != NULL);

      /* child should have valid and non-empty bounds */
      assert(!(child->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDRELAXED));
      assert(!SCIPintervalIsEmpty(infinity, child->bounds));
      /* nodes at depth 0 are always linear */
      assert(child->depth > 0 || child->curv == SCIP_EXPRCURV_LINEAR);

      childbounds[i] = child->bounds;  /*lint !e644*/
      childcurv[i]   = child->curv;    /*lint !e644*/
   }

   /* check curvature */
   *curv = SCIPexprcurvMonomial(monomial->nfactors, monomial->exponents, NULL, childcurv, childbounds);
   *curv = SCIPexprcurvMultiply(monomial->coef, *curv);

   /* free memory, if allocated before */
TERMINATE:
   if( childbounds != childboundsstatic )
   {
      BMSfreeMemoryArrayNull(&childbounds);
      BMSfreeMemoryArrayNull(&childcurv);
   }

   return retcode;
}

/** gives the user data belonging to a SCIP_EXPR_USER expression */
SCIP_USEREXPRDATA* SCIPexprgraphGetNodeUserData(
   SCIP_EXPRGRAPHNODE*   node
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_USER);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_USER*)node->data.data)->userdata;
}

/** indicates whether a user expression has the estimator callback defined */
SCIP_Bool SCIPexprgraphHasNodeUserEstimator(
   SCIP_EXPRGRAPHNODE*   node
   )
{
   assert(node != NULL);
   assert(node->op == SCIP_EXPR_USER);
   assert(node->data.data != NULL);

   return ((SCIP_EXPRDATA_USER*)node->data.data)->estimate != NULL;
}

/** gets bounds of a node in an expression graph */
SCIP_INTERVAL SCIPexprgraphGetNodeBounds(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->bounds;
}

/** gets value of expression associated to node from last evaluation call */
SCIP_Real SCIPexprgraphGetNodeVal(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->value;
}

/** gets curvature of expression associated to node from last curvature check call */
SCIP_EXPRCURV SCIPexprgraphGetNodeCurvature(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   assert(node != NULL);

   return node->curv;
}

/** creates an expression graph node */
SCIP_RETCODE SCIPexprgraphCreateNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   SCIP_EXPROP           op,                 /**< operator type of expression */
   ...
   )
{
   va_list         ap;
   SCIP_EXPROPDATA opdata;

   assert(blkmem != NULL);
   assert(node   != NULL);

   *node = NULL;

   switch( op )
   {
   case SCIP_EXPR_VARIDX    :
   case SCIP_EXPR_PARAM     :
   case SCIP_EXPR_CONST     :
   case SCIP_EXPR_LINEAR    :
   case SCIP_EXPR_QUADRATIC :
   case SCIP_EXPR_POLYNOMIAL:
   case SCIP_EXPR_USER      :
   {
      SCIPerrorMessage("cannot create node with operand %d via SCIPexprgraphCreateNode\n");
      SCIPABORT();
      return SCIP_ERROR;  /*lint !e527*/
   }

   /* operands without data */
   case SCIP_EXPR_PLUS   :
   case SCIP_EXPR_MINUS  :
   case SCIP_EXPR_MUL    :
   case SCIP_EXPR_DIV    :
   case SCIP_EXPR_MIN    :
   case SCIP_EXPR_MAX    :
   case SCIP_EXPR_SQUARE :
   case SCIP_EXPR_SQRT   :
   case SCIP_EXPR_EXP    :
   case SCIP_EXPR_LOG    :
   case SCIP_EXPR_SIN    :
   case SCIP_EXPR_COS    :
   case SCIP_EXPR_TAN    :
      /* case SCIP_EXPR_ERF : */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_ABS    :
   case SCIP_EXPR_SIGN   :
   case SCIP_EXPR_SUM    :
   case SCIP_EXPR_PRODUCT:
      opdata.data = NULL;
      break;

   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_SIGNPOWER:
   {
      va_start(ap, op );  /*lint !e838*/
      opdata.dbl = va_arg( ap, SCIP_Real);  /*lint !e416 !e826*/
      va_end( ap );  /*lint !e826*/

      break;
   }

   case SCIP_EXPR_INTPOWER:
   {
      va_start(ap, op );  /*lint !e838*/
      opdata.intval = va_arg( ap, int);  /*lint !e416 !e826*/
      va_end( ap );  /*lint !e826*/

      break;
   }

   case SCIP_EXPR_LAST:
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   SCIP_CALL( exprgraphCreateNode(blkmem, node, op, opdata) ); /*lint !e644*/

   return SCIP_OKAY;
}

/** creates an expression graph node for a linear expression */
SCIP_RETCODE SCIPexprgraphCreateNodeLinear(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   int                   ncoefs,             /**< number of coefficients */
   SCIP_Real*            coefs,              /**< coefficients of linear expression */
   SCIP_Real             constant            /**< constant of linear expression */
   )
{
   SCIP_EXPROPDATA opdata;
   SCIP_Real* data;

   assert(blkmem != NULL);
   assert(node   != NULL);

   /* we store the coefficients and the constant in a single array and make this our operand data */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &data, ncoefs + 1) );
   BMScopyMemoryArray(data, coefs, ncoefs);  /*lint !e644*/
   data[ncoefs] = constant;

   opdata.data = data;
   SCIP_CALL( exprgraphCreateNode(blkmem, node, SCIP_EXPR_LINEAR, opdata) );

   return SCIP_OKAY;
}

/** creates an expression graph node for a quadratic expression */
SCIP_RETCODE SCIPexprgraphCreateNodeQuadratic(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   int                   nchildren,          /**< number of children */
   SCIP_Real*            lincoefs,           /**< linear coefficients for children, or NULL */
   int                   nquadelems,         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems,          /**< quadratic elements, or NULL if nquadelems == 0 */
   SCIP_Real             constant            /**< constant */
   )
{
   SCIP_EXPROPDATA opdata;
   SCIP_EXPRDATA_QUADRATIC* data;

   assert(blkmem != NULL);
   assert(node   != NULL);
   assert(quadelems != NULL || nquadelems == 0);

   SCIP_CALL( quadraticdataCreate(blkmem, &data, constant, nchildren, lincoefs, nquadelems, quadelems) );

   opdata.data = data;
   SCIP_CALL( exprgraphCreateNode(blkmem, node, SCIP_EXPR_QUADRATIC, opdata) );

   return SCIP_OKAY;
}

/** creates an expression graph node for a polynomial expression */
SCIP_RETCODE SCIPexprgraphCreateNodePolynomial(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   int                   nmonomials,         /**< number of monomials */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< monomials */
   SCIP_Real             constant,           /**< constant of polynomial */
   SCIP_Bool             copymonomials       /**< whether to copy monomials or to assume ownership */
   )
{
   SCIP_EXPROPDATA opdata;
   SCIP_EXPRDATA_POLYNOMIAL* data;

   assert(blkmem != NULL);
   assert(node   != NULL);
   assert(monomials != NULL || nmonomials == 0);

   SCIP_CALL( polynomialdataCreate(blkmem, &data, nmonomials, monomials, constant, copymonomials) );

   opdata.data = data;
   SCIP_CALL( exprgraphCreateNode(blkmem, node, SCIP_EXPR_POLYNOMIAL, opdata) );

   return SCIP_OKAY;
}

/** adds monomials to an expression graph node that is a polynomial expression */
SCIP_RETCODE SCIPexprgraphNodePolynomialAddMonomials(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE*   node,               /**< store expression graph node with polynomial operator */
   int                   nmonomials,         /**< number of monomials */
   SCIP_EXPRDATA_MONOMIAL** monomials,       /**< monomials */
   SCIP_Bool             copymonomials       /**< whether to copy monomials or to assume ownership */
   )
{
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_POLYNOMIAL);
   assert(monomials != NULL || nmonomials == 0);

   SCIP_CALL( polynomialdataAddMonomials(blkmem, (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data, nmonomials, monomials, copymonomials) );

   return SCIP_OKAY;
}

/** creates an expression graph node for a user expression */
SCIP_RETCODE SCIPexprgraphCreateNodeUser(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node,               /**< buffer to store expression graph node */
   SCIP_USEREXPRDATA*    data,               /**< user data for expression, node assumes ownership */
   SCIP_EXPRINTCAPABILITY evalcapability,    /**< evaluation capability */
   SCIP_DECL_USEREXPREVAL    ((*eval)),      /**< evaluation function */
   SCIP_DECL_USEREXPRINTEVAL ((*inteval)),   /**< interval evaluation function */
   SCIP_DECL_USEREXPRCURV    ((*curv)),      /**< curvature check function */
   SCIP_DECL_USEREXPRPROP    ((*prop)),      /**< interval propagation function */
   SCIP_DECL_USEREXPRESTIMATE ((*estimate)), /**< estimation function, or NULL if convex, concave, or not implemented */
   SCIP_DECL_USEREXPRCOPYDATA ((*copydata)), /**< expression data copy function, or NULL if nothing to copy */
   SCIP_DECL_USEREXPRFREEDATA ((*freedata)), /**< expression data free function, or NULL if nothing to free */
   SCIP_DECL_USEREXPRPRINT ((*print))        /**< expression print function, or NULL for default string "user" */
   )
{
   SCIP_EXPROPDATA opdata;
   SCIP_EXPRDATA_USER* exprdata;

   assert(blkmem != NULL);
   assert(node   != NULL);
   assert(eval != NULL);
   assert((evalcapability & SCIP_EXPRINTCAPABILITY_FUNCVALUE) != 0);  /* the function evaluation is not optional */
   assert(((evalcapability & SCIP_EXPRINTCAPABILITY_INTFUNCVALUE) == 0) || inteval != NULL);  /* if capability says it can do interval evaluation, then the corresponding callback needs to be provided */
   assert(copydata != NULL || data == NULL);
   assert(freedata != NULL || data == NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &exprdata) );

   exprdata->userdata = data;
   exprdata->evalcapability = evalcapability;
   exprdata->eval = eval;
   exprdata->estimate = estimate;
   exprdata->inteval = inteval;
   exprdata->curv = curv;
   exprdata->prop = prop;
   exprdata->copydata = copydata;
   exprdata->freedata = freedata;
   exprdata->print = print;

   opdata.data = (void*) exprdata;

   SCIP_CALL( exprgraphCreateNode(blkmem, node, SCIP_EXPR_USER, opdata) );

   return SCIP_OKAY;
}

/** given a node of an expression graph, splitup a linear part which variables are not used somewhere else in the same expression
 *
 *  E.g., if the expression is 1 + x + y + y^2, one gets 1 + x and the node remains at y + y^2.
 *  If the node is a linear expression, it may be freed.
 *  If it is not linear, the node may change, i.e., the remaining nonlinear part may be stored in a new node.
 *  It is assumed that the user had captured the node.
 *  It is assumed that the expression graph has been simplified before.
 */
SCIP_RETCODE SCIPexprgraphNodeSplitOffLinear(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  node,               /**< expression graph node where to splitup linear part */
   int                   linvarssize,        /**< length of linvars and lincoefs arrays */
   int*                  nlinvars,           /**< buffer to store length of linear term that have been splitup */
   void**                linvars,            /**< buffer to store variables of linear part */
   SCIP_Real*            lincoefs,           /**< buffer to store coefficients of linear part */
   SCIP_Real*            constant            /**< buffer to store constant part */
   )
{
   int orignvars;
   int* varsusage;
   SCIP_EXPRGRAPHNODE* orignode;
   SCIP_Bool havechange;
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(*node != NULL);
   assert((*node)->nuses > 0);
   assert(nlinvars != NULL);
   assert(linvars  != NULL || linvarssize == 0);
   assert(lincoefs != NULL || linvarssize == 0);
   assert(constant != NULL);

   *constant = 0.0;
   *nlinvars = 0;

   SCIPdebugMessage("split off linear part for %s node %p (%d,%d)\n", SCIPexpropGetName((*node)->op), (void*)*node, (*node)->depth, (*node)->pos);

   /* do some obvious and easy cases */
   switch( (*node)->op )
   {
   case SCIP_EXPR_VARIDX:
   {
      if( linvarssize >= 1 )
      {
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0;                                    /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
      }
      return SCIP_OKAY;
   }

   case SCIP_EXPR_CONST:
   {
      *constant = (*node)->data.dbl;
      SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );

      return SCIP_OKAY;
   }

   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_SIGNPOWER:
   {
      if( (*node)->data.dbl == 1.0 && (*node)->children[0]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[0]->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0;                                                 /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
      }
      return SCIP_OKAY;
   }

   case SCIP_EXPR_INTPOWER:
   {
      if( (*node)->data.intval == 1 && (*node)->children[0]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[0]->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0;                                                 /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
      }
      return SCIP_OKAY;
   }

   case SCIP_EXPR_PLUS:
   {
      if( (*node)->children[0]->op == SCIP_EXPR_CONST && (*node)->children[1]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *constant = (*node)->children[0]->data.dbl;
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[1]->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0;                                                 /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );

         return SCIP_OKAY;
      }
      else if( (*node)->children[1]->op == SCIP_EXPR_CONST && (*node)->children[0]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *constant = (*node)->children[1]->data.dbl;
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[0]->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0;                                                 /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );

         return SCIP_OKAY;
      }
      else if( (*node)->children[0]->op == SCIP_EXPR_VARIDX && (*node)->children[1]->op == SCIP_EXPR_VARIDX && linvarssize >= 2 )
      {
         *nlinvars = 2;
         linvars[0]  = exprgraph->vars[(*node)->children[0]->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0;                                                 /*lint !e613*/
         linvars[1]  = exprgraph->vars[(*node)->children[1]->data.intval];  /*lint !e613*/
         lincoefs[1] = 1.0;                                                 /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );

         return SCIP_OKAY;
      }
      else if( ((*node)->children[0]->op == SCIP_EXPR_VARIDX || (*node)->children[1]->op == SCIP_EXPR_VARIDX) && linvarssize >= 1 )
      {
         /* handle this one later */
         break;
      }
      return SCIP_OKAY;
   }

   case SCIP_EXPR_MINUS:
   {
      if( (*node)->children[0]->op == SCIP_EXPR_CONST && (*node)->children[1]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *constant = (*node)->children[0]->data.dbl;
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[1]->data.intval];  /*lint !e613*/
         lincoefs[0] = -1.0;                                                /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );

         return SCIP_OKAY;
      }
      else if( (*node)->children[1]->op == SCIP_EXPR_CONST && (*node)->children[0]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *constant = -(*node)->children[1]->data.dbl;
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[0]->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0;                                                 /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );

         return SCIP_OKAY;
      }
      else if( (*node)->children[0]->op == SCIP_EXPR_VARIDX && (*node)->children[1]->op == SCIP_EXPR_VARIDX && linvarssize >= 2 )
      {
         *nlinvars = 2;
         linvars[0]  = exprgraph->vars[(*node)->children[0]->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0;                                                 /*lint !e613*/
         linvars[1]  = exprgraph->vars[(*node)->children[1]->data.intval];  /*lint !e613*/
         lincoefs[1] = -1.0;                                                /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );

         return SCIP_OKAY;
      }
      else if( ((*node)->children[0]->op == SCIP_EXPR_VARIDX || (*node)->children[1]->op == SCIP_EXPR_VARIDX) && linvarssize >= 1 )
      {
         /* handle this one later */
         break;
      }
      return SCIP_OKAY;
   }

   case SCIP_EXPR_MUL:
   {
      if( (*node)->children[0]->op == SCIP_EXPR_CONST && (*node)->children[1]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[1]->data.intval];  /*lint !e613*/
         lincoefs[0] = (*node)->children[0]->data.dbl;                      /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
      }
      else if( (*node)->children[1]->op == SCIP_EXPR_CONST && (*node)->children[0]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[0]->data.intval];  /*lint !e613*/
         lincoefs[0] = (*node)->children[1]->data.dbl;                      /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
      }
      return SCIP_OKAY;
   }

   case SCIP_EXPR_DIV:
   {
      if( (*node)->children[1]->op != SCIP_EXPR_CONST )
         return SCIP_OKAY;

      if( (*node)->children[0]->op == SCIP_EXPR_VARIDX && linvarssize >= 1 )
      {
         *nlinvars = 1;
         linvars[0]  = exprgraph->vars[(*node)->children[0]->data.intval];  /*lint !e613*/
         lincoefs[0] = 1.0/(*node)->children[1]->data.dbl;                  /*lint !e613*/

         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
      }
      return SCIP_OKAY;
   }

   case SCIP_EXPR_SQUARE:
   case SCIP_EXPR_SQRT:
   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
      return SCIP_OKAY;

   case SCIP_EXPR_PRODUCT:
   case SCIP_EXPR_USER:
      return SCIP_OKAY;

   case SCIP_EXPR_SUM:
   case SCIP_EXPR_LINEAR:
   case SCIP_EXPR_QUADRATIC:
   case SCIP_EXPR_POLYNOMIAL:
   default:
   {
      /* check if there is a child that is a variable */
      for( i = 0; i < (*node)->nchildren; ++i )
      {
         if( (*node)->children[i]->op == SCIP_EXPR_VARIDX )
            break;
      }

      if( i == (*node)->nchildren )
         return SCIP_OKAY;

      break;
   }
   }  /*lint !e788*/

   /* count how often variables are used in this expression */
   assert(exprgraph->nvars > 0); /* in a simplified expr graph with no variables, there can only be const nodes, but these were handled above */
   orignvars = exprgraph->nvars;
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &varsusage, exprgraph->nvars) );
   BMSclearMemoryArray(varsusage, exprgraph->nvars);  /*lint !e644*/

   exprgraphNodeGetVarsUsage(*node, varsusage);

   /* duplicate node if it has parents or more than one user */
   orignode = NULL;
   if( (*node)->nparents > 0 || (*node)->nuses > 1 )
   {
      SCIP_EXPROPDATA data;

      orignode = *node;

      if( exprOpTable[orignode->op].copydata != NULL )
      {
         SCIP_CALL( exprOpTable[orignode->op].copydata(exprgraph->blkmem, orignode->nchildren, orignode->data, &data) );
      }
      else
         data = orignode->data;

      SCIP_CALL( exprgraphCreateNode(exprgraph->blkmem, node, orignode->op, data) );
      SCIP_CALL( SCIPexprgraphAddNode(exprgraph, *node, -1, orignode->nchildren, orignode->children) );
      SCIPexprgraphCaptureNode(*node);
   }

   havechange = FALSE;
   /* split up constant and linear part */
   switch( (*node)->op )
   {
   case SCIP_EXPR_PLUS:
   case SCIP_EXPR_MINUS:
   {
      SCIP_EXPRGRAPHNODE* varchild;
      SCIP_EXPRGRAPHNODE* otherchild;
      int varidx;

      /* we had looked at this above already and only continued if exactly one node is still a child and linvarssize is >= 1 */
      assert((*node)->children[0]->op == SCIP_EXPR_VARIDX || (*node)->children[1]->op == SCIP_EXPR_VARIDX);
      assert((*node)->children[0]->op != SCIP_EXPR_VARIDX || (*node)->children[1]->op != SCIP_EXPR_VARIDX);
      assert(linvarssize >= 1);

      varchild   = (*node)->children[0]->op == SCIP_EXPR_VARIDX ? (*node)->children[0] : (*node)->children[1];
      otherchild = (*node)->children[0]->op == SCIP_EXPR_VARIDX ? (*node)->children[1] : (*node)->children[0];
      varidx = varchild->data.intval;
      /* if variable is used in other child (which should be nonlinear), we don't take it */
      if( varsusage[varidx] > 1 )
         break;

      /* add to linear variables */
      *nlinvars = 1;
      linvars[0]  = exprgraph->vars[varidx];  /*lint !e613*/
      if( (*node)->op == SCIP_EXPR_MINUS && varchild == (*node)->children[1] )
         lincoefs[0] = -1.0;  /*lint !e613*/
      else
         lincoefs[0] =  1.0;  /*lint !e613*/

      if( (*node)->op == SCIP_EXPR_PLUS || (*node)->children[0] == otherchild )
      {
         /* replace *node by otherchild */
         SCIPexprgraphCaptureNode(otherchild);
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
         *node = otherchild;
      }
      else
      {
         SCIP_Real* lindata;

         /* turn *node into linear expression -1.0 * otherchild */

         /* reduce to one child */
         SCIP_ALLOC( BMSreallocBlockMemoryArray(exprgraph->blkmem, &(*node)->children, 2, 1) );  /*lint !e506*/
         (*node)->children[0] = otherchild;
         (*node)->nchildren = 1;
         (*node)->op = SCIP_EXPR_LINEAR;

         /* setup linear data -1.0 * child0 + 0.0 */
         SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &lindata, 2) );  /*lint !e506*/
         lindata[0] = -1.0;
         lindata[1] =  0.0;
         (*node)->data.data = (void*)lindata;

         /* remove *node as parent of varchild */
         SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &varchild, *node) );
      }

      havechange = TRUE;

      break;
   }

   case SCIP_EXPR_SUM:
   {
      int nchildren;

      i = 0;
      nchildren = (*node)->nchildren;
      while( i < nchildren )
      {
         /* sort out constants */
         if( (*node)->children[i]->op == SCIP_EXPR_CONST )
         {
            *constant += (*node)->children[i]->data.dbl;
            SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &(*node)->children[i], *node) );

            if( i < nchildren-1 )
            {
               (*node)->children[i] = (*node)->children[nchildren-1];
               (*node)->children[nchildren-1] = NULL;
            }
            --nchildren;

            continue;
         }

         /* keep every child that is not a constant or variable */
         if( (*node)->children[i]->op != SCIP_EXPR_VARIDX )
         {
            ++i;
            continue;
         }

         /* skip variables that are used in other parts of the expression */
         if( varsusage[(*node)->children[i]->data.intval] > 1 )
         {
            ++i;
            continue;
         }

         /* move variable into linear part, if still space */
         if( *nlinvars < linvarssize )
         {
            linvars[*nlinvars]  = exprgraph->vars[(*node)->children[i]->data.intval];  /*lint !e613*/
            lincoefs[*nlinvars] = 1.0;                                                 /*lint !e613*/
            ++*nlinvars;

            SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &(*node)->children[i], *node) );
            if( i < nchildren-1 )
            {
               (*node)->children[i] = (*node)->children[nchildren-1];
               (*node)->children[nchildren-1] = NULL;
            }
            --nchildren;

            continue;
         }
      }
      assert(i == nchildren);

      if( nchildren == 0 )
      {
         /* all children were removed */
         havechange = TRUE;
         BMSfreeBlockMemoryArray(exprgraph->blkmem, &(*node)->children, (*node)->nchildren);
         (*node)->nchildren = 0;
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
         break;
      }

      if( nchildren < (*node)->nchildren )
      {
         /* some children were removed */
         havechange = TRUE;
         SCIP_ALLOC( BMSreallocBlockMemoryArray(exprgraph->blkmem, &(*node)->children, (*node)->nchildren, nchildren) );
         (*node)->nchildren = nchildren;
      }

      if( havechange && (*node)->nchildren == 1 )
      {
         /* replace node by its child */
         SCIP_EXPRGRAPHNODE* child;

         child = (*node)->children[0];
         SCIPexprgraphCaptureNode(child);
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
         *node = child;

         break;
      }

      break;
   }

   case SCIP_EXPR_LINEAR:
   {
      int nchildren;
      SCIP_Real* coefs;

      coefs = (SCIP_Real*)(*node)->data.data;
      assert(coefs != NULL);

      /* remove constant, if nonzero */
      if( coefs[(*node)->nchildren] != 0.0 )
      {
         *constant = coefs[(*node)->nchildren];
         coefs[(*node)->nchildren] = 0.0;
         havechange = TRUE;
      }

      i = 0;
      nchildren = (*node)->nchildren;
      while( i < nchildren )
      {
         /* sort out constants */
         if( (*node)->children[i]->op == SCIP_EXPR_CONST )
         {
            *constant += coefs[i] * (*node)->children[i]->data.dbl;
            SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &(*node)->children[i], *node) );

            if( i < nchildren-1 )
            {
               (*node)->children[i] = (*node)->children[nchildren-1];
               (*node)->children[nchildren-1] = NULL;
               coefs[i] = coefs[nchildren-1];
               coefs[nchildren-1] = 0.0;
            }
            --nchildren;

            continue;
         }

         /* keep everything that is not a constant or variable */
         if( (*node)->children[i]->op != SCIP_EXPR_VARIDX )
         {
            ++i;
            continue;
         }

         /* skip variables that are used in other parts of the expression */
         if( varsusage[(*node)->children[i]->data.intval] > 1 )
         {
            ++i;
            continue;
         }

         /* move variable into linear part, if still space */
         if( *nlinvars < linvarssize )
         {
            linvars[*nlinvars]  = exprgraph->vars[(*node)->children[i]->data.intval];  /*lint !e613*/
            lincoefs[*nlinvars] = coefs[i];                                            /*lint !e613*/
            ++*nlinvars;

            SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &(*node)->children[i], *node) );
            if( i < nchildren-1 )
            {
               (*node)->children[i] = (*node)->children[nchildren-1];
               (*node)->children[nchildren-1] = NULL;
               coefs[i] = coefs[nchildren-1];
               coefs[nchildren-1] = 0.0;
            }
            --nchildren;

            continue;
         }
      }
      assert(i == nchildren);

      if( nchildren == 0 )
      {
         /* all children were removed */
         havechange = TRUE;
         BMSfreeBlockMemoryArray(exprgraph->blkmem, &(*node)->children, (*node)->nchildren);
         BMSfreeBlockMemoryArray(exprgraph->blkmem, &coefs, (*node)->nchildren+1);
         (*node)->data.data = NULL;
         (*node)->nchildren = 0;
         (*node)->op = SCIP_EXPR_SUM; /* because we freed the constraint data already */
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
         break;
      }

      if( nchildren < (*node)->nchildren )
      {
         /* some children were removed */
         havechange = TRUE;
         SCIP_ALLOC( BMSreallocBlockMemoryArray(exprgraph->blkmem, &(*node)->children, (*node)->nchildren, nchildren) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(exprgraph->blkmem, &coefs, (*node)->nchildren+1, nchildren+1) );
         coefs[nchildren] = 0.0;
         (*node)->data.data = (void*)coefs;
         (*node)->nchildren = nchildren;
      }

      if( havechange && (*node)->nchildren == 1 && coefs[0] == 1.0 )
      {
         /* replace node by its child */
         SCIP_EXPRGRAPHNODE* child;

         child = (*node)->children[0];
         SCIPexprgraphCaptureNode(child);
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
         *node = child;

         break;
      }

      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* quaddata;
      SCIP_Bool* childused;
      int* childmap;
      int nchildren;

      quaddata = (SCIP_EXPRDATA_QUADRATIC*)(*node)->data.data;
      assert(quaddata != NULL);

      /* remove constant, if nonzero */
      if( quaddata->constant != 0.0 )
      {
         *constant = quaddata->constant;
         quaddata->constant = 0.0;
         havechange = TRUE;
      }

      /* if there is no linear part or no space left for linear variables, then stop */
      if( quaddata->lincoefs != NULL || linvarssize == 0 )
         break;

      /* check which childs are used in quadratic terms */
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childused, (*node)->nchildren) );
      BMSclearMemoryArray(childused, (*node)->nchildren);  /*lint !e644*/

      for( i = 0; i < quaddata->nquadelems; ++i )
      {
         childused[quaddata->quadelems[i].idx1] = TRUE;
         childused[quaddata->quadelems[i].idx2] = TRUE;
      }

      /* alloc space for mapping of children indices */
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childmap, (*node)->nchildren) );

      nchildren = (*node)->nchildren;
      for( i = 0; i < nchildren; ++i )
      {
         childmap[i] = i;  /*lint !e644*/
         if( *nlinvars >= linvarssize )
            continue;
         /* skip child if not variable or also used in quadratic part or other parts of expression */
         if( (*node)->children[i]->op != SCIP_EXPR_VARIDX )
            continue;
         if( childused[i] )
            continue;
         if( varsusage[(*node)->children[i]->data.intval] > 1 )
            continue;

         /* put variable into linear part */
         linvars[*nlinvars]  = exprgraph->vars[(*node)->children[i]->data.intval];  /*lint !e613*/
         lincoefs[*nlinvars] = quaddata->lincoefs[i];                               /*lint !e613*/
         quaddata->lincoefs[i] = 0.0;
         ++*nlinvars;

         SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &(*node)->children[i], *node) );

         /* move last child to position i */
         if( i < nchildren-1 )
         {
            (*node)->children[i] = (*node)->children[nchildren-1];
            quaddata->lincoefs[i] = quaddata->lincoefs[nchildren-1];
            childused[i] = childused[nchildren-1];
            childmap[nchildren-1] = i;
         }
         --nchildren;
         childmap[i] = -1;

         havechange = TRUE;
         --i; /* look at i again */
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &childused, (*node)->nchildren);  /*lint !e850*/

      if( nchildren < (*node)->nchildren )
      {
         /* apply childmap to quadratic term */
         for( i = 0; i < quaddata->nquadelems; ++i )
         {
            quaddata->quadelems[i].idx1 = childmap[quaddata->quadelems[i].idx1];
            quaddata->quadelems[i].idx2 = childmap[quaddata->quadelems[i].idx2];
            if( quaddata->quadelems[i].idx1 > quaddata->quadelems[i].idx2 )
            {
               int tmp;
               tmp = quaddata->quadelems[i].idx1;
               quaddata->quadelems[i].idx1 = quaddata->quadelems[i].idx2;
               quaddata->quadelems[i].idx2 = tmp;
            }
         }
         quaddata->sorted = FALSE;
      }
      BMSfreeBlockMemoryArray(exprgraph->blkmem, &childmap, (*node)->nchildren);

      if( nchildren == 0 )
      {
         /* all children were removed (so it was actually a linear expression) */
         havechange = TRUE;
         BMSfreeBlockMemoryArray(exprgraph->blkmem, &(*node)->children, (*node)->nchildren);
         exprFreeDataQuadratic(exprgraph->blkmem, (*node)->nchildren, (*node)->data);
         (*node)->data.data = NULL;
         (*node)->nchildren = 0;
         (*node)->op = SCIP_EXPR_SUM;
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
         break;
      }

      if( nchildren < (*node)->nchildren )
      {
         /* reduce number of children */
         SCIP_ALLOC( BMSreallocBlockMemoryArray(exprgraph->blkmem, &(*node)->children, (*node)->nchildren, nchildren) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(exprgraph->blkmem, &quaddata->lincoefs, (*node)->nchildren, nchildren) );
         (*node)->nchildren = nchildren;
      }

      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* polynomialdata;
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_Bool* childused;
      int childidx;
      int j;

      polynomialdata = (SCIP_EXPRDATA_POLYNOMIAL*)(*node)->data.data;
      assert(polynomialdata != NULL);

      /* make sure linear monomials are merged */
      polynomialdataMergeMonomials(exprgraph->blkmem, polynomialdata, 0.0, FALSE);

      /* remove constant, if nonzero */
      if( polynomialdata->constant != 0.0 )
      {
         *constant = polynomialdata->constant;
         polynomialdata->constant = 0.0;
         havechange = TRUE;
      }

      /* if there is no space for linear variables, then stop */
      if( linvarssize == 0 )
         break;

      /* get nonlinear child usage: how often each child is used in the polynomial in a nonlinear monomial */
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childused, (*node)->nchildren) );
      BMSclearMemoryArray(childused, (*node)->nchildren);  /*lint !e644*/
      for( i = 0; i < polynomialdata->nmonomials; ++i )
      {
         monomial = polynomialdata->monomials[i];
         assert(monomial != NULL);
         if( monomial->nfactors == 0 || (monomial->nfactors == 1 && monomial->exponents[0] == 1.0) )
            continue;
         for( j = 0; j < monomial->nfactors; ++j )
         {
            assert(monomial->childidxs[j] >= 0);
            assert(monomial->childidxs[j] < (*node)->nchildren);
            childused[monomial->childidxs[j]] = TRUE;
         }
      }

      /* move linear monomials out of polynomial */
      for( i = 0; i < polynomialdata->nmonomials && *nlinvars < linvarssize; ++i )
      {
         monomial = polynomialdata->monomials[i];
         assert(monomial != NULL);

         /* sort out constants */
         if( monomial->nfactors == 0 )
         {
            if( monomial->coef != 0.0 )
            {
               *constant += monomial->coef;
               havechange = TRUE;
            }
            continue;
         }

         if( monomial->nfactors != 1 )
            continue;
         if( monomial->exponents[0] != 1.0 )
            continue;
         childidx = monomial->childidxs[0];
         assert((*node)->children[childidx] != NULL); /* should be due to merge in the beginning */
         if( (*node)->children[childidx]->op != SCIP_EXPR_VARIDX )
            continue;
         if( childused[childidx] || varsusage[(*node)->children[childidx]->data.intval] > 1 )
            continue;

         /* we are at a linear monomial in a variable that is not used somewhere else in nonlinear form */

         /* put variable into linear part */
         linvars[*nlinvars]  = exprgraph->vars[(*node)->children[childidx]->data.intval];  /*lint !e613*/
         lincoefs[*nlinvars] = monomial->coef;                                             /*lint !e613*/
         ++*nlinvars;

         monomial->coef = 0.0;
         monomial->nfactors = 0;
         polynomialdata->sorted = FALSE;

         SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &(*node)->children[childidx], *node) );
         (*node)->children[childidx] = NULL;

         havechange = TRUE;
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &childused, (*node)->nchildren);

      if( *nlinvars > 0 )
      {
         /* if we did something, cleanup polynomial (e.g., remove monomials with coefficient 0.0) */
         polynomialdataMergeMonomials(exprgraph->blkmem, polynomialdata, 0.0, FALSE);
         SCIP_CALL( exprgraphNodeRemovePolynomialNullChildren(exprgraph->blkmem, *node) );
      }

      if( (*node)->nchildren == 0 )
      {
         assert(polynomialdata->nmonomials == 0);
         assert(polynomialdata->constant == 0.0);
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
         havechange = TRUE;
         break;
      }

      break;
   }

   default: ;
   }  /*lint !e788*/

   BMSfreeBlockMemoryArray(exprgraph->blkmem, &varsusage, orignvars);

   if( orignode != NULL )
   {
      /* if node was duplicated, we need to forget about original or duplicate */
      if( !havechange )
      {
         /* if nothing has changed, then forget about duplicate */
         assert(*constant == 0.0);
         assert(*nlinvars == 0);
         assert(*node != NULL);
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, node) );
         *node = orignode;
      }
      else
      {
         /* if something changed, then release original node */
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, &orignode) );
      }
   }
   else if( havechange && *node != NULL )
   {
      /* if node was not duplicated and not removed but changed, then invalidate value, bounds, and simplified status */
      (*node)->value = SCIP_INVALID;
      (*node)->simplified = FALSE;
      (*node)->boundstatus = SCIP_EXPRBOUNDSTATUS_CHILDTIGHTENED;
      SCIPintervalSetEntire(SCIP_REAL_MAX, &(*node)->bounds);
      exprgraph->needvarboundprop = TRUE;
   }

   return SCIP_OKAY;
}

/** moves parents from a one node to another node
 *
 *  In other words, replaces the child srcnode by targetnode in all parents of srcnode.
 *  srcnode may be freed, if not captured.
 *  It is assumed that targetnode represents the same expression as srcnode.
 */
SCIP_RETCODE SCIPexprgraphMoveNodeParents(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  srcnode,            /**< node which parents to move */
   SCIP_EXPRGRAPHNODE*   targetnode          /**< node where to move parents to */
   )
{
   assert(exprgraph != NULL);
   assert(srcnode != NULL);
   assert(*srcnode != NULL);
   assert(targetnode != NULL);

   while( *srcnode != NULL && (*srcnode)->nparents > 0 )
   {
      if( (*srcnode)->parents[0]->depth <= targetnode->depth )
      {
         SCIP_CALL( exprgraphMoveNode(exprgraph, (*srcnode)->parents[0], targetnode->depth+1) );
      }
      SCIP_CALL( exprgraphNodeReplaceChild(exprgraph, (*srcnode)->parents[0], srcnode, targetnode) );
   }
   assert(*srcnode == NULL || (*srcnode)->nuses > 0);

   return SCIP_OKAY;
}

/** releases node, i.e., decreases number of uses
 *
 *  node is freed if no parents and no other uses.
 *  Children are recursively released if they have no other parents.
 *  Nodes that are removed are also freed.
 *  If node correspond to a variable, then the variable is removed from the expression graph;
 *  similarly for constants.
 */
SCIP_RETCODE SCIPexprgraphReleaseNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  node                /**< expression graph node to release */
   )
{
   int i;

   assert(exprgraph != NULL);
   assert(node  != NULL);
   assert(*node != NULL);
   assert((*node)->depth >= 0); /* node should be in graph */
   assert((*node)->pos   >= 0); /* node should be in graph */
   assert((*node)->depth < exprgraph->depth);
   assert((*node)->pos < exprgraph->nnodes[(*node)->depth]);
   assert((*node)->nuses >= 1);
   assert(exprgraph->nodes[(*node)->depth][(*node)->pos] == *node);

   SCIPdebugMessage("release node %p\n", (void*)*node);

   --(*node)->nuses;

   /* do nothing if node still has parents or is still in use */
   if( (*node)->nparents > 0 || (*node)->nuses > 0 )
   {
      SCIPdebugMessage("skip removing node %p (%d, %d) with %d parents and %d uses from expression graph\n", (void*)*node, (*node)->depth, (*node)->pos, (*node)->nparents, (*node)->nuses);
      *node = NULL;
      return SCIP_OKAY;
   }

   SCIPdebugMessage("remove node %p (%d, %d) with op %s from expression graph\n", (void*)*node, (*node)->depth, (*node)->pos, SCIPexpropGetName((*node)->op));

   /* notify children about removal of its parent
    * they are also freed, if possible */
   for( i = 0; i < (*node)->nchildren; ++i )
   {
      SCIP_CALL( exprgraphNodeRemoveParent(exprgraph, &(*node)->children[i], *node) );
      (*node)->children[i] = NULL;
   }

   if( (*node)->op == SCIP_EXPR_VARIDX )
   {
      assert((*node)->depth == 0);
      SCIP_CALL( exprgraphRemoveVar(exprgraph, (*node)->data.intval) );
   }
   else if( (*node)->op == SCIP_EXPR_CONST && (*node)->depth == 0 )
   {
      int constidx;

      (void) exprgraphFindConstNodePos(exprgraph, *node, &constidx);
      assert(constidx >= 0);
      assert(constidx < exprgraph->nconsts);
      assert(exprgraph->constnodes[constidx] == *node);

      /* move last constant to position constidx */
      if( constidx < exprgraph->nconsts-1 )
      {
         exprgraph->constnodes[constidx] = exprgraph->constnodes[exprgraph->nconsts-1];
         exprgraph->constssorted = (exprgraph->nconsts <= 2);
      }
      --exprgraph->nconsts;
   }
   else
   {
      /* only variables and constants are allowed at depth 0 */
      assert((*node)->depth > 0);
   }

   /* remove node from nodes array in expression graph */
   if( (*node)->pos < exprgraph->nnodes[(*node)->depth]-1 )
   {
      /* move last node at depth of *node to position of *node */
      exprgraph->nodes[(*node)->depth][(*node)->pos] = exprgraph->nodes[(*node)->depth][exprgraph->nnodes[(*node)->depth]-1];
      exprgraph->nodes[(*node)->depth][(*node)->pos]->pos = (*node)->pos;

      /* moving the node may change the order in the parents array of each child */
      for( i = 0; i < exprgraph->nodes[(*node)->depth][(*node)->pos]->nchildren; ++i )
         exprgraph->nodes[(*node)->depth][(*node)->pos]->children[i]->parentssorted = FALSE;
   }
   --exprgraph->nnodes[(*node)->depth];

   /* node is now not in graph anymore */
   (*node)->depth = -1;
   (*node)->pos   = -1;

   /* free node */
   SCIPexprgraphFreeNode(exprgraph->blkmem, node);

   *node = NULL;

   return SCIP_OKAY;
}

/* @todo should be a private method and node creation should already capture a node instead of waiting that it's added to the graph */
/** frees a node of an expression graph */
void SCIPexprgraphFreeNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPHNODE**  node                /**< pointer to expression graph node that should be freed */
   )
{
   assert(blkmem != NULL);
   assert( node != NULL);
   assert(*node != NULL);
   assert((*node)->depth == -1); /* node should not be in graph anymore */
   assert((*node)->pos   == -1); /* node should not be in graph anymore */
   assert((*node)->nuses == 0); /* node should not be in use */

   /* free operator data, if needed */
   if( exprOpTable[(*node)->op].freedata != NULL )
      exprOpTable[(*node)->op].freedata(blkmem, (*node)->nchildren, (*node)->data);

   /* free arrays of children and parent nodes */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*node)->children, (*node)->nchildren);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*node)->parents,  (*node)->parentssize);

   /* free node struct */
   BMSfreeBlockMemory(blkmem, node);
}

/** enables a node and recursively all its children in an expression graph */
void SCIPexprgraphEnableNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node to enable */
   )
{
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->depth >= 0);
   assert(node->pos >= 0);

   if( node->enabled )
      return;

   SCIPdebugMessage("enable node %p (%d,%d)\n", (void*)node, node->depth, node->pos);

   node->enabled = TRUE;
   for( i = 0; i < node->nchildren; ++i )
      SCIPexprgraphEnableNode(exprgraph, node->children[i]);

   /* make sure bounds are updated in next bound propagation round */
   SCIPintervalSetEntire(SCIP_REAL_MAX, &node->bounds);
   exprgraph->needvarboundprop = TRUE;
}

/** disables a node and recursively all children which have no enabled parents in an expression graph */
void SCIPexprgraphDisableNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node to enable */
   )
{
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->depth >= 0);
   assert(node->pos >= 0);

   if( !node->enabled )
      return;

   /* workaround: don't disable nodes if there could be more users than the one who is disabling the node
    * otherwise, if there are several nonlinear constraints using the same expression graph node as root node,
    * we might get enabled constraints with disabled node
    */
   if( node->nuses > 1 )
      return;

   /* if all parents of node are disabled, then also node can be disabled */
   node->enabled = FALSE;
   for( i = 0; i < node->nparents; ++i )
      if( node->parents[i]->enabled )
      {
         node->enabled = TRUE;
         return;
      }

   SCIPdebugMessage("disabled node %p (%d,%d), nuses = %d\n", (void*)node, node->depth, node->pos, node->nuses);

   for( i = 0; i < node->nchildren; ++i )
      SCIPexprgraphDisableNode(exprgraph, node->children[i]);
}

/** returns whether the node has siblings in the expression graph */
SCIP_Bool SCIPexprgraphHasNodeSibling(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   int p;

   assert(node != NULL);

   for( p = 0; p < node->nparents; ++p )
      if( node->parents[p]->nchildren > 1 )
         return TRUE;

   return FALSE;
}

/** returns whether all children of an expression graph node are variable nodes
 *
 *  Returns TRUE for nodes without children.
 */
SCIP_Bool SCIPexprgraphAreAllNodeChildrenVars(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   int i;

   assert(node != NULL);

   for( i = 0; i < node->nchildren; ++i )
      if( node->children[i]->op != SCIP_EXPR_VARIDX )
         return FALSE;

   return TRUE;
}

/** returns whether the node has an ancestor which has a nonlinear expression operand */
SCIP_Bool SCIPexprgraphHasNodeNonlinearAncestor(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   int p;

   for( p = 0; p < node->nparents; ++p )
   {
      assert(node->parents[p]->depth > node->depth);
      switch( node->parents[p]->op )
      {
      case SCIP_EXPR_PLUS:
      case SCIP_EXPR_MINUS:
      case SCIP_EXPR_SUM:
      case SCIP_EXPR_LINEAR:
         if( SCIPexprgraphHasNodeNonlinearAncestor(node->parents[p]) )
            return TRUE;
         break;

#ifndef NDEBUG
      case SCIP_EXPR_VARIDX:
      case SCIP_EXPR_CONST:
      case SCIP_EXPR_PARAM:
         assert(0); /* these expressions cannot have children */
         break;
#endif

      default:
         /* parent has nonlinear expression operand */
         return TRUE;
      }/*lint !e788*/
   }

   return FALSE;
}

/** prints an expression graph node */
void SCIPexprgraphPrintNode(
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   )
{
   assert(node != NULL);

   exprgraphPrintNodeExpression(node, messagehdlr, file, NULL, FALSE);
}

/** tightens the bounds in a node of the graph
 *
 *  Preparation for reverse propagation.
 *  Sets bound status to SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTRECENT if tightening is strong enough and not cutoff.
 */
void SCIPexprgraphTightenNodeBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< node in expression graph with no parents */
   SCIP_INTERVAL         nodebounds,         /**< new bounds for node */
   SCIP_Real             minstrength,        /**< minimal required relative bound strengthening in a node to trigger a propagation into children nodes (set to negative value if propagation should always be triggered) */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Bool*            cutoff              /**< buffer to store whether a node's bounds were propagated to an empty interval */
   )
{
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->depth >= 0);
   assert(node->pos >= 0);
   assert(!SCIPintervalIsEmpty(infinity, nodebounds));
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* if node is disabled, then ignore new bounds */
   if( !node->enabled )
   {
      SCIPdebugMessage("ignore bound tightening for node %p (%d,%d)\n", (void*)node, node->depth, node->pos);
      return;
   }

   SCIPdebugMessage("tighten bounds of node %p (%d,%d) from [%10g, %10g] by [%10g, %10g]",
      (void*)node, node->depth, node->pos,
      node->bounds.inf, node->bounds.sup, nodebounds.inf, nodebounds.sup);

   /* bounds in node should be valid */
   assert(!(node->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDRELAXED));

   if( nodebounds.inf > node->bounds.sup || nodebounds.sup < node->bounds.inf )
   {
      *cutoff = TRUE;
      SCIPdebugPrintf(" -> cutoff\n");
      return;
   }

   /* if minstrength is negative, always mark that node has recently tightened bounds,
    * if bounds are considerably improved or tightening leads to an empty interval,
    * mark that node has recently tightened bounds
    * if bounds are only slightly improved, set the status to tightened by parent,
    * so next propagateVarBound round will reset the bounds
    */
   if( minstrength < 0.0 )
      node->boundstatus |= SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTFORCE;
   else if(
      isLbBetter(minstrength, nodebounds.inf, node->bounds.inf, node->bounds.sup) ||
      isUbBetter(minstrength, nodebounds.sup, node->bounds.inf, node->bounds.sup) )
      node->boundstatus |= SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENTRECENT;
   else if( nodebounds.inf > node->bounds.inf || nodebounds.sup < node->bounds.sup )
      node->boundstatus |= SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENT;

   SCIPintervalIntersect(&node->bounds, node->bounds, nodebounds);
   SCIPdebugPrintf(" -> [%10g, %10g] status %d\n", node->bounds.inf, node->bounds.sup, node->boundstatus);
}

/** ensures that bounds and curvature information in a node is uptodate
 *
 *  Assumes that bounds and curvature in children are uptodate.
 */
SCIP_RETCODE SCIPexprgraphUpdateNodeBoundsCurvature(
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Real             minstrength,        /**< minimal required relative bound strengthening to trigger a bound recalculation in parent nodes */
   SCIP_Bool             clearreverseprop    /**< whether to reset bound tightenings from reverse propagation */
   )
{
   SCIP_INTERVAL  childboundsstatic[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_EXPRCURV  childcurvstatic[SCIP_EXPRESSION_MAXCHILDEST];
   SCIP_INTERVAL* childbounds = NULL;
   SCIP_EXPRCURV* childcurv = NULL;
   SCIP_RETCODE retcode = SCIP_OKAY;
   int i;

   assert(node != NULL);
   assert(node->depth >= 0); /* node should be in graph */
   assert(node->pos >= 0);   /* node should be in graph */
   assert(node->enabled);    /* node should be enabled, otherwise we may not have uptodate bounds and curvatures in children */

   if( node->depth == 0 )
   {
      /* we cannot update bound tightenings in variable nodes here */
      assert(!clearreverseprop || !(node->boundstatus & SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENT));
      return SCIP_OKAY;
   }

   assert(node->op != SCIP_EXPR_VARIDX);
   assert(node->op != SCIP_EXPR_PARAM);

   /* if many children, get large enough memory to store children bounds */
   if( node->nchildren > SCIP_EXPRESSION_MAXCHILDEST )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&childbounds, node->nchildren) );
      SCIP_ALLOC_TERMINATE(retcode, BMSallocMemoryArray(&childcurv, node->nchildren), TERMINATE);
   }
   else
   {
      childbounds = childboundsstatic;
      childcurv   = childcurvstatic;
   }

   /* assemble bounds and curvature of children */
   for( i = 0; i < node->nchildren; ++i )
   {
      /* child should have valid and non-empty bounds */
      assert(!(node->children[i]->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDRELAXED));
      assert(!SCIPintervalIsEmpty(infinity, node->children[i]->bounds));
      /* nodes at depth 0 are always linear */
      assert(node->children[i]->depth > 0 || node->children[i]->curv == SCIP_EXPRCURV_LINEAR);

      childbounds[i] = node->children[i]->bounds;  /*lint !e644*/
      childcurv[i]   = node->children[i]->curv;    /*lint !e644*/
   }

   /* if we do not have valid bounds, then update
    * code below is copied from exprgraphNodeUpdateBounds */
   if( node->boundstatus != SCIP_EXPRBOUNDSTATUS_VALID )
   {
      SCIP_INTERVAL newbounds;

      /* calling interval evaluation function for this operand */
      assert( exprOpTable[node->op].inteval != NULL );
      SCIP_CALL_TERMINATE( retcode, exprOpTable[node->op].inteval(infinity, node->data, node->nchildren, childbounds, NULL, NULL, &newbounds), TERMINATE );

      /* if bounds of a children were relaxed or our bounds were tightened by a (now possibly invalid) reverse propagation from a parent
       * and now our bounds are relaxed, then we have to propagate this upwards to ensure valid bounds
       *
       * if bounds were tightened (considerably), then tell this to those parents which think that they have valid bounds
       *
       * finally, if there was only a little tightening, then keep this updated bounds, but don't notify parents
       */
      if( (newbounds.inf < node->bounds.inf || newbounds.sup > node->bounds.sup) &&
         ((node->boundstatus & SCIP_EXPRBOUNDSTATUS_CHILDRELAXED) || ((node->boundstatus & SCIP_EXPRBOUNDSTATUS_TIGHTENEDBYPARENT) && clearreverseprop)) )
      {
         for( i = 0; i < node->nparents; ++i )
            node->parents[i]->boundstatus = SCIP_EXPRBOUNDSTATUS_CHILDRELAXED;

         node->bounds = newbounds;
      }
      else if( isLbBetter(minstrength, newbounds.inf, node->bounds.inf, node->bounds.sup) ||
         (     isUbBetter(minstrength, newbounds.sup, node->bounds.inf, node->bounds.sup)) )
      {
         for( i = 0; i < node->nparents; ++i )
            node->parents[i]->boundstatus |= SCIP_EXPRBOUNDSTATUS_CHILDTIGHTENED;

         node->bounds = newbounds;
      }
      else
      {
         SCIPintervalIntersect(&node->bounds, node->bounds, newbounds);
      }

      SCIPdebugMessage("updated bounds of node %p (%d,%d) op %s to [%g,%g]\n", (void*)node, node->depth, node->pos, SCIPexpropGetName(node->op), node->bounds.inf, node->bounds.sup);

      /* node now has valid bounds */
      node->boundstatus = SCIP_EXPRBOUNDSTATUS_VALID;
   }

   /* update curvature */
   if( SCIPintervalIsEmpty(infinity, node->bounds) )
   {
      node->curv = SCIP_EXPRCURV_LINEAR;

      SCIPdebugMessage("node %p(%d,%d) has empty domain in SCIPexprgraphUpdateNodeBoundsCurvature\n", (void*)node, node->depth, node->pos);
   }
   else
   {
      SCIP_CALL_TERMINATE( retcode, exprOpTable[node->op].curv(infinity, node->data, node->nchildren, childbounds, childcurv, &node->curv), TERMINATE );

      /* SCIPdebugMessage("curvature %s for %s = ", SCIPexprcurvGetName(node->curv), SCIPexpropGetName(node->op));
       * SCIPdebug( exprgraphPrintNodeExpression(node, NULL, NULL, TRUE) );
       * SCIPdebugPrintf("\n");
       */
   }
TERMINATE:
   /* free memory, if allocated before */
   if( childbounds != childboundsstatic )
   {
      BMSfreeMemoryArrayNull(&childbounds);
      BMSfreeMemoryArrayNull(&childcurv);
   }

   return retcode;
}

/**@} */

/**@name Expression graph methods */
/**@{ */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPexprgraphGetDepth
#undef SCIPexprgraphGetNNodes
#undef SCIPexprgraphGetNodes
#undef SCIPexprgraphGetNVars
#undef SCIPexprgraphGetVars
#undef SCIPexprgraphGetVarNodes
#undef SCIPexprgraphSetVarNodeValue
#undef SCIPexprgraphSetVarsBounds
#undef SCIPexprgraphSetVarBounds
#undef SCIPexprgraphSetVarNodeBounds
#undef SCIPexprgraphSetVarNodeLb
#undef SCIPexprgraphSetVarNodeUb
#undef SCIPexprgraphGetVarsBounds

/** get current maximal depth of expression graph */
int SCIPexprgraphGetDepth(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   )
{
   assert(exprgraph != NULL);

   return exprgraph->depth;
}

/** gets array with number of nodes at each depth of expression graph */
int* SCIPexprgraphGetNNodes(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   )
{
   assert(exprgraph != NULL);

   return exprgraph->nnodes;
}

/** gets nodes of expression graph, one array per depth */
SCIP_EXPRGRAPHNODE*** SCIPexprgraphGetNodes(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   )
{
   assert(exprgraph != NULL);

   return exprgraph->nodes;
}

/** gets number of variables in expression graph */
int SCIPexprgraphGetNVars(
   SCIP_EXPRGRAPH*       exprgraph           /**< pointer to expression graph that should be freed */
   )
{
   assert(exprgraph != NULL);

   return exprgraph->nvars;
}

/** gets array of variables in expression graph */
void** SCIPexprgraphGetVars(
   SCIP_EXPRGRAPH*       exprgraph           /**< pointer to expression graph that should be freed */
   )
{
   assert(exprgraph != NULL);

   return exprgraph->vars;
}

/** gets array of expression graph nodes corresponding to variables */
SCIP_EXPRGRAPHNODE** SCIPexprgraphGetVarNodes(
   SCIP_EXPRGRAPH*       exprgraph           /**< pointer to expression graph that should be freed */
   )
{
   assert(exprgraph != NULL);

   return exprgraph->varnodes;
}

/** sets value for a single variable given as expression graph node */
void SCIPexprgraphSetVarNodeValue(
   SCIP_EXPRGRAPHNODE*   varnode,            /**< expression graph node corresponding to variable */
   SCIP_Real             value               /**< new value for variable */
   )
{
   assert(varnode != NULL);
   assert(varnode->op == SCIP_EXPR_VARIDX);

   varnode->value = value;
}

/** sets bounds for variables */
void SCIPexprgraphSetVarsBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_INTERVAL*        varbounds           /**< new bounds for variables */
   )
{
   assert(exprgraph != NULL);
   assert(varbounds != NULL || exprgraph->nvars == 0);

   BMScopyMemoryArray(exprgraph->varbounds, varbounds, exprgraph->nvars);
}

/** sets bounds for a single variable */
void SCIPexprgraphSetVarBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   void*                 var,                /**< variable */
   SCIP_INTERVAL         varbounds           /**< new bounds of variable */
   )
{
   int pos;

   assert(exprgraph != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(exprgraph->varidxs, var));

   pos = (int)(size_t)SCIPhashmapGetImage(exprgraph->varidxs, var);
   assert(pos < exprgraph->nvars);
   assert(exprgraph->vars[pos] == var);

   exprgraph->varbounds[pos] = varbounds;
}

/** sets bounds for a single variable given as expression graph node */
void SCIPexprgraphSetVarNodeBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   varnode,            /**< expression graph node corresponding to variable */
   SCIP_INTERVAL         varbounds           /**< new bounds of variable */
   )
{
   int pos;

   assert(exprgraph != NULL);
   assert(varnode != NULL);

   pos = varnode->data.intval;
   assert(pos >= 0);
   assert(pos < exprgraph->nvars);
   assert(exprgraph->varnodes[pos] == varnode);

   exprgraph->varbounds[pos] = varbounds;
}

/** sets lower bound for a single variable given as expression graph node */
void SCIPexprgraphSetVarNodeLb(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   varnode,            /**< expression graph node corresponding to variable */
   SCIP_Real             lb                  /**< new lower bound for variable */
   )
{
   int pos;

   assert(exprgraph != NULL);
   assert(varnode != NULL);

   pos = varnode->data.intval;
   assert(pos >= 0);
   assert(pos < exprgraph->nvars);
   assert(exprgraph->varnodes[pos] == varnode);

   exprgraph->varbounds[pos].inf = lb;
}

/** sets upper bound for a single variable given as expression graph node */
void SCIPexprgraphSetVarNodeUb(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   varnode,            /**< expression graph node corresponding to variable */
   SCIP_Real             ub                  /**< new upper bound for variable */
   )
{
   int pos;

   assert(exprgraph != NULL);
   assert(varnode != NULL);

   pos = varnode->data.intval;
   assert(pos >= 0);
   assert(pos < exprgraph->nvars);
   assert(exprgraph->varnodes[pos] == varnode);

   exprgraph->varbounds[pos].sup = ub;
}

/** gets bounds that are stored for all variables */
SCIP_INTERVAL* SCIPexprgraphGetVarsBounds(
   SCIP_EXPRGRAPH*       exprgraph           /**< expression graph */
   )
{
   return exprgraph->varbounds;
}

/** creates an empty expression graph */
SCIP_RETCODE SCIPexprgraphCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRGRAPH**      exprgraph,          /**< buffer to store pointer to expression graph */
   int                   varssizeinit,       /**< minimal initial size for variables array, or -1 to choose automatically */
   int                   depthinit,          /**< minimal initial depth of expression graph, or -1 to choose automatically */
   SCIP_DECL_EXPRGRAPHVARADDED((*exprgraphvaradded)), /**< callback method to invoke when a variable has been added to the expression graph, or NULL if not needed */
   SCIP_DECL_EXPRGRAPHVARREMOVE((*exprgraphvarremove)), /**< callback method to invoke when a variable will be removed from the expression graph, or NULL if not needed */
   SCIP_DECL_EXPRGRAPHVARCHGIDX((*exprgraphvarchgidx)), /**< callback method to invoke when a variable changes its index in the expression graph, or NULL if not needed */
   void*                 userdata            /**< user data to pass to callback functions */
   )
{
   assert(blkmem != NULL);
   assert(exprgraph != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, exprgraph) );
   BMSclearMemory(*exprgraph);
   (*exprgraph)->blkmem = blkmem;

   /* create nodes's arrays */
   SCIP_CALL( exprgraphEnsureDepth(*exprgraph, MAX(1, depthinit)) );
   assert((*exprgraph)->depth >= 1);

   /* create var's arrays and hashmap */
   ensureBlockMemoryArraySize3((*exprgraph)->blkmem, &(*exprgraph)->varnodes, &(*exprgraph)->vars, &(*exprgraph)->varbounds, &(*exprgraph)->varssize, varssizeinit);
   SCIP_CALL( SCIPhashmapCreate(&(*exprgraph)->varidxs, (*exprgraph)->blkmem, (*exprgraph)->varssize) );

   /* empty array of constants is sorted */
   (*exprgraph)->constssorted = TRUE;

   /* store callback functions and user data */
   (*exprgraph)->exprgraphvaradded = exprgraphvaradded;
   (*exprgraph)->exprgraphvarremove = exprgraphvarremove;
   (*exprgraph)->exprgraphvarchgidx = exprgraphvarchgidx;
   (*exprgraph)->userdata = userdata;

   return SCIP_OKAY;
}

/** frees an expression graph */
SCIP_RETCODE SCIPexprgraphFree(
   SCIP_EXPRGRAPH**      exprgraph           /**< pointer to expression graph that should be freed */
   )
{
   BMS_BLKMEM* blkmem;
   int d;

   assert( exprgraph != NULL);
   assert(*exprgraph != NULL);
   assert((*exprgraph)->nvars == 0);
   assert((*exprgraph)->nconsts == 0);

   blkmem = (*exprgraph)->blkmem;
   assert(blkmem != NULL);

   /* free nodes arrays */
   for( d = 0; d < (*exprgraph)->depth; ++d )
   {
      assert((*exprgraph)->nnodes[d] == 0);
      BMSfreeBlockMemoryArrayNull(blkmem, &(*exprgraph)->nodes[d], (*exprgraph)->nodessize[d]);  /*lint !e866*/
   }
   assert((*exprgraph)->nodes     != NULL);
   assert((*exprgraph)->nnodes    != NULL);
   assert((*exprgraph)->nodessize != NULL);
   BMSfreeBlockMemoryArray(blkmem, &(*exprgraph)->nodes,     (*exprgraph)->depth);
   BMSfreeBlockMemoryArray(blkmem, &(*exprgraph)->nnodes,    (*exprgraph)->depth);
   BMSfreeBlockMemoryArray(blkmem, &(*exprgraph)->nodessize, (*exprgraph)->depth);

   /* free variables arrays and hashmap */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*exprgraph)->vars,      (*exprgraph)->varssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*exprgraph)->varnodes,  (*exprgraph)->varssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*exprgraph)->varbounds, (*exprgraph)->varssize);
   SCIPhashmapFree(&(*exprgraph)->varidxs);

   /* free constants array */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*exprgraph)->constnodes, (*exprgraph)->constssize);

   /* free graph struct */
   BMSfreeBlockMemory(blkmem, exprgraph);

   return SCIP_OKAY;
}

/** adds an expression graph node to an expression graph
 *
 *  Expression graph assumes ownership of node.
 *  Children are notified about new parent.
 *  Depth will be chosen to be the maximum of mindepth and the depth of all children plus one.
 */
SCIP_RETCODE SCIPexprgraphAddNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node to add */
   int                   mindepth,           /**< minimal depth in expression graph where to add node, e.g., 0 or smaller to choose automatically */
   int                   nchildren,          /**< number of children */
   SCIP_EXPRGRAPHNODE**  children            /**< children nodes, or NULL if no children */
   )
{
   SCIP_Bool childvalsvalid;
   int depth;
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->pos < 0);          /* node should have no position in graph yet */
   assert(node->depth < 0);        /* node should have no position in graph yet */
   assert(node->nchildren == 0);   /* node should not have stored children yet */
   assert(node->children == NULL); /* node should not have stored children yet */
   assert(node->nparents == 0);    /* node should not have parents stored yet */
   assert(children != NULL || nchildren == 0);

   /* choose depth as maximal depth of children + 1, and at least mindepth */
   depth = MAX(0, mindepth);
   for( i = 0; i < nchildren; ++i )
   {
      if( children[i]->depth >= depth )  /*lint !e613*/
         depth = children[i]->depth + 1;  /*lint !e613*/
   }

   /* ensure that expression graph is deep enough */
   SCIP_CALL( exprgraphEnsureDepth(exprgraph, depth+1) );
   assert(exprgraph->depth > depth);

   /* ensure enough space for nodes at depth depth */
   ensureBlockMemoryArraySize(exprgraph->blkmem, &exprgraph->nodes[depth], &exprgraph->nodessize[depth], exprgraph->nnodes[depth]+1);  /*lint !e866*/

   /* add node to graph */
   node->depth = depth;
   node->pos   = exprgraph->nnodes[depth];
   exprgraph->nodes[depth][node->pos] = node;
   ++exprgraph->nnodes[depth];

   /* add as parent to children
    * and check if children has valid values */
   childvalsvalid = TRUE;
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( exprgraphNodeAddParent(exprgraph->blkmem, children[i], node) );  /*lint !e613*/
      childvalsvalid &= (children[i]->value != SCIP_INVALID);  /*lint !e777 !e514 !e613*/
   }
   /* store children */
   if( nchildren > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(exprgraph->blkmem, &node->children, children, nchildren) );
      node->nchildren = nchildren;
   }

   if( node->op == SCIP_EXPR_CONST )
   {
      /* set bounds to constant value of node */
      node->boundstatus = SCIP_EXPRBOUNDSTATUS_VALID;
      SCIPintervalSet(&node->bounds, node->data.dbl);
   }
   else
   {
      /* set bounds to entire, set status to "should recompute", and note that we have to update something */
      node->boundstatus = SCIP_EXPRBOUNDSTATUS_CHILDTIGHTENED;
      SCIPintervalSetEntire(SCIP_REAL_MAX, &node->bounds);
      exprgraph->needvarboundprop = TRUE;
   }

   /* if not a variable, set value of node according to values of children (if all have valid values) */
   if( node->op != SCIP_EXPR_VARIDX && childvalsvalid )
   {
      SCIP_CALL( exprgraphNodeEval(node, NULL) );
   }

   return SCIP_OKAY;
}

/** adds variables to an expression graph, if not existing yet
 *
 *  Also already existing nodes are enabled.
 */
SCIP_RETCODE SCIPexprgraphAddVars(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   int                   nvars,              /**< number of variables to add */
   void**                vars,               /**< variables to add */
   SCIP_EXPRGRAPHNODE**  varnodes            /**< array to store nodes corresponding to variables, or NULL if not of interest */
   )
{
   SCIP_EXPRGRAPHNODE* node;
   SCIP_EXPROPDATA opdata;
   int i;

   assert(exprgraph != NULL);
   assert(exprgraph->depth >= 1);
   assert(vars != NULL || nvars == 0);

   /* if there are no variables yet, then it's quite likely that we will create new nodes for all vars, so can easily estimate how much space we will need in variables array and nodes at depth 0 arrays */
   if( exprgraph->nvars == 0 )
   {
      ensureBlockMemoryArraySize3(exprgraph->blkmem, &exprgraph->vars, &exprgraph->varnodes, &exprgraph->varbounds, &exprgraph->varssize, exprgraph->nvars + nvars);
      ensureBlockMemoryArraySize(exprgraph->blkmem, &exprgraph->nodes[0], &exprgraph->nodessize[0], exprgraph->nnodes[0] + nvars);
   }

   for( i = 0; i < nvars; ++i )
   {
      /* skip variables that exist already */
      if( SCIPhashmapExists(exprgraph->varidxs, vars[i]) )/*lint !e613*/
      {
         (void) SCIPexprgraphFindVarNode(exprgraph, vars[i], &node);  /*lint !e613*/
         assert(node != NULL);

         /* enable node */
         node->enabled = TRUE;

         if( varnodes != NULL )
            varnodes[i] = node;

         continue;
      }

      /* create new variable expression */
      opdata.intval = exprgraph->nvars;
      SCIP_CALL( exprgraphCreateNode(exprgraph->blkmem, &node, SCIP_EXPR_VARIDX, opdata) );

      /* add expression node to expression graph at depth 0 */
      SCIP_CALL( SCIPexprgraphAddNode(exprgraph, node, 0, 0, NULL) );

      /* add variable node to vars arrays and hashmap */
      ensureBlockMemoryArraySize3(exprgraph->blkmem, &exprgraph->vars, &exprgraph->varnodes, &exprgraph->varbounds, &exprgraph->varssize, exprgraph->nvars + 1);
      exprgraph->vars[exprgraph->nvars] = vars[i];  /*lint !e613*/
      exprgraph->varnodes[exprgraph->nvars] = node;
      SCIPintervalSetEntire(SCIP_REAL_MAX, &exprgraph->varbounds[exprgraph->nvars]);
      SCIP_CALL( SCIPhashmapInsert(exprgraph->varidxs, vars[i], (void*)(size_t)exprgraph->nvars) );  /*lint !e613*/
      ++exprgraph->nvars;

      if( varnodes != NULL )
         varnodes[i] = node;

      SCIPdebugMessage("added node %p (%d, %d) for new variable %d\n", (void*)node, node->depth, node->pos, node->data.intval);

      /* call callback method, if set */
      if( exprgraph->exprgraphvaradded != NULL )
      {
         SCIP_CALL( exprgraph->exprgraphvaradded(exprgraph, exprgraph->userdata, vars[i], node) );  /*lint !e613*/
      }
   }

   return SCIP_OKAY;
}

/** adds a constant to an expression graph, if not existing yet
 *
 *  Also already existing nodes are enabled.
 */
SCIP_RETCODE SCIPexprgraphAddConst(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             constant,           /**< constant to add */
   SCIP_EXPRGRAPHNODE**  constnode           /**< buffer to store pointer to expression graph node corresponding to constant */
   )
{
   SCIP_EXPROPDATA opdata;

   assert(exprgraph != NULL);
   assert(constnode != NULL);

   /* check if there is already an expression for this constant */
   if( SCIPexprgraphFindConstNode(exprgraph, constant, constnode) )
   {
      assert(*constnode != NULL);
      assert((*constnode)->op == SCIP_EXPR_CONST);
      assert((*constnode)->data.dbl == constant);  /*lint !e777*/
      (*constnode)->enabled = TRUE;
      return SCIP_OKAY;
   }

   /* create new node for constant */
   opdata.dbl = constant;
   SCIP_CALL( exprgraphCreateNode(exprgraph->blkmem, constnode, SCIP_EXPR_CONST, opdata) );

   /* add node to expression graph at depth 0 */
   SCIP_CALL( SCIPexprgraphAddNode(exprgraph, *constnode, 0, 0, NULL) );
   assert((*constnode)->depth == 0);
   assert((*constnode)->pos   == exprgraph->nnodes[0]-1);

   /* add node to constnodes arrays; @todo should use SCIPsortedvecInsertPtr? */
   ensureBlockMemoryArraySize(exprgraph->blkmem, &exprgraph->constnodes, &exprgraph->constssize, exprgraph->nconsts + 1);
   exprgraph->constnodes[exprgraph->nconsts] = *constnode;
   ++exprgraph->nconsts;
   exprgraph->constssorted = exprgraph->nconsts <= 1 || (exprgraph->constssorted && exprgraphConstNodeComp(exprgraph->constnodes[exprgraph->nconsts-2], *constnode) < 0);

   SCIPdebugMessage("added node %p (%d, %d) for new constant %g\n", (void*)constnode, (*constnode)->depth, (*constnode)->pos, (*constnode)->data.dbl);

   return SCIP_OKAY;
}

/** adds sum of expression trees into expression graph
 *
 *  node will also be captured.
 *
 *  @note Parameters will be converted into constants
 */
SCIP_RETCODE SCIPexprgraphAddExprtreeSum(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   int                   nexprtrees,         /**< number of expression trees to add */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees that should be added */
   SCIP_Real*            coefs,              /**< coefficients of expression trees, or NULL if all 1.0 */
   SCIP_EXPRGRAPHNODE**  rootnode,           /**< buffer to store expression graph node corresponding to root of expression tree */
   SCIP_Bool*            rootnodeisnew       /**< buffer to indicate whether the node in *rootnode has been newly created for this expression tree (otherwise, expression tree was already in graph) */
   )
{
   SCIP_Bool allone;

   assert(exprgraph != NULL);
   assert(nexprtrees > 0);
   assert(exprtrees != NULL);
   assert(rootnode != NULL);
   assert(rootnodeisnew != NULL);

   *rootnode = NULL;

   if( nexprtrees == 1 && (coefs == NULL || coefs[0] == 1.0) )
   {
      assert(exprtrees[0] != NULL);
      assert(exprtrees[0]->vars != NULL || exprtrees[0]->nvars == 0);

      SCIP_CALL( exprgraphAddExpr(exprgraph, exprtrees[0]->root, exprtrees[0]->vars, exprtrees[0]->params, rootnode, rootnodeisnew) );
   }
   else
   {
      SCIP_EXPROP op;
      SCIP_EXPRGRAPHNODE** rootnodes;
      SCIP_Bool rootnodeisnew_;
      int i;

      *rootnodeisnew = TRUE;
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &rootnodes, nexprtrees) );

      allone = TRUE;
      for( i = 0; i < nexprtrees; ++i )
      {
         assert(exprtrees[i] != NULL);
         assert(exprtrees[i]->vars != NULL || exprtrees[i]->nvars == 0);

         SCIP_CALL( exprgraphAddExpr(exprgraph, exprtrees[i]->root, exprtrees[i]->vars, exprtrees[i]->params, &rootnodes[i], &rootnodeisnew_) );  /*lint !e644*/
         assert(rootnodes[i] != NULL);
         *rootnodeisnew &= rootnodeisnew_;

         allone &= (coefs == NULL || coefs[i] == 1.0);  /*lint !e514*/
      }

      /* decide which operand we want to use for the root node */
      if( coefs == NULL || allone )
         op = nexprtrees == 2 ? SCIP_EXPR_PLUS : SCIP_EXPR_SUM;
      else if( nexprtrees == 2 && coefs[0] == 1.0 && coefs[1] == -1.0 )
         op = SCIP_EXPR_MINUS;
      else if( nexprtrees == 2 && coefs[1] == 1.0 && coefs[0] == -1.0 )
      {
         SCIP_EXPRGRAPHNODE* tmp;

         tmp = rootnodes[0];
         rootnodes[0] = rootnodes[1];
         rootnodes[1] = tmp;
         op = SCIP_EXPR_MINUS;
      }
      else
         op = SCIP_EXPR_LINEAR;

      if( op != SCIP_EXPR_LINEAR )
      {
         SCIP_EXPROPDATA data;
         data.data = NULL;

         if( !*rootnodeisnew )
         {
            /* all exprtrees were in the graph already, check if we also already have a node for the sum of rootnodes */
            SCIP_CALL( exprgraphFindParentByOperator(exprgraph, nexprtrees, rootnodes, op, data, NULL, rootnode) );
         }

         if( *rootnode == NULL )
         {
            /* create new node for sum of rootnodes and add to exprgraph */
            SCIP_CALL( exprgraphCreateNode(exprgraph->blkmem, rootnode, op, data) );
            SCIP_CALL( SCIPexprgraphAddNode(exprgraph, *rootnode, -1, nexprtrees, rootnodes) );
            *rootnodeisnew = TRUE;
         }
         else
         {
            /* apparently we already had a node for the sum of rootnodes, so signal this to the caller */
            *rootnodeisnew = FALSE;
         }
      }
      else
      {
         SCIP_EXPROPDATA data;
         SCIP_Real* lindata;

         assert(op == SCIP_EXPR_LINEAR);

         /* setup data for linear expression: copy of coefficients and 0.0 as constant term */
         SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &lindata, nexprtrees+1) );
         BMScopyMemoryArray(lindata, coefs, nexprtrees);  /*lint !e644*/
         lindata[nexprtrees] = 0.0;
         data.data = lindata;

         if( !*rootnodeisnew )
         {
            /* all exprtrees were in the graph already, check if we also already have a node for the linear combination of rootnodes */
            SCIP_CALL( exprgraphFindParentByOperator(exprgraph, nexprtrees, rootnodes, SCIP_EXPR_LINEAR, data, NULL, rootnode) );
         }

         if( *rootnode == NULL )
         {
            /* create new node for linear combination of rootnodes and add to exprgraph */
            SCIP_CALL( exprgraphCreateNode(exprgraph->blkmem, rootnode, SCIP_EXPR_LINEAR, data) );
            SCIP_CALL( SCIPexprgraphAddNode(exprgraph, *rootnode, -1, nexprtrees, rootnodes) );
            *rootnodeisnew = TRUE;
         }
         else
         {
            /* apparently we already had a node for the sum of rootnodes, so signal this to the caller */
            *rootnodeisnew = FALSE;
            BMSfreeBlockMemoryArray(exprgraph->blkmem, &lindata, nexprtrees+1);
         }
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &rootnodes, nexprtrees);
   }
   assert(*rootnode != NULL);

   SCIPexprgraphCaptureNode(*rootnode);

   SCIPdebugMessage("%s node %p (%d,%d) is root for %d expression trees\n",
      rootnodeisnew ? "new" : "old", (void*)*rootnode, (*rootnode)->depth, (*rootnode)->pos, nexprtrees);

   return SCIP_OKAY;
}

/** replaces variable in expression graph by a linear sum of variables
 *
 *  Variables will be added if not in the graph yet.
 */
SCIP_RETCODE SCIPexprgraphReplaceVarByLinearSum(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   void*                 var,                /**< variable to replace */
   int                   ncoefs,             /**< number of coefficients in linear term */
   SCIP_Real*            coefs,              /**< coefficients in linear term, or NULL if ncoefs == 0 */
   void**                vars,               /**< variables in linear term */
   SCIP_Real             constant            /**< constant offset */
   )
{
   SCIP_EXPRGRAPHNODE* varnode;
   SCIP_Real* lindata;
   int varidx;
   int i;

   assert(exprgraph != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(exprgraph->varidxs, var));
   assert(coefs != NULL || ncoefs == 0);
   assert(vars  != NULL || ncoefs == 0);

   varidx = (int)(size_t)SCIPhashmapGetImage(exprgraph->varidxs, var);
   assert(varidx < exprgraph->nvars);
   assert(exprgraph->vars[varidx] == var);
   varnode = exprgraph->varnodes[varidx];
   assert(varnode != NULL);
   assert(varnode->data.intval == varidx);

   if( ncoefs == 0 || (ncoefs == 1 && constant == 0.0 && coefs[0] == 1.0) )  /*lint !e613*/
   {
      /* variable is replaced by constant or variable */
      SCIP_EXPRGRAPHNODE* node;

      /* check if there is already a node for this constant or variable */
      node = NULL;
      if( ncoefs == 0 )
      {
         (void)SCIPexprgraphFindConstNode(exprgraph, constant, &node);
         assert(node == NULL || node->data.dbl == constant);  /*lint !e777*/
      }
      else
      {
         (void)SCIPexprgraphFindVarNode(exprgraph, vars[0], &node);  /*lint !e613*/
         assert(node == NULL || exprgraph->vars[node->data.intval] == vars[0]);  /*lint !e613*/
      }

      if( node != NULL )
      {
         SCIPdebugMessage("try to replace varnode %p (%d uses, %d parents) by %p\n", (void*)varnode, varnode->nuses, varnode->nparents, (void*)node);

         /* tell parents of varnode to replace child varnode by node, this may free varnode */
         SCIP_CALL( SCIPexprgraphMoveNodeParents(exprgraph, &varnode, node) );

         /* if varnode is in use, turn it into SCIP_EXPR_SUM with node as only child */
         if( varnode != NULL )
         {
            assert(varnode->nuses > 0);
            assert(varnode->nparents == 0);

            /* remove variable (but don't free it's node) from graph */
            SCIP_CALL( exprgraphRemoveVar(exprgraph, varidx) );

            /* move varnode up to depth 1 */
            SCIP_CALL( exprgraphMoveNode(exprgraph, varnode, 1) );

            /* turn into EXPR_SUM expression */
            varnode->op = SCIP_EXPR_SUM;
            varnode->data.data = NULL;
            SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &varnode->children, 1) );  /*lint !e506*/
            varnode->children[0] = node;
            varnode->nchildren = 1;
            SCIP_CALL( exprgraphNodeAddParent(exprgraph->blkmem, node, varnode) );

            varnode->value = node->value;
            varnode->bounds = node->bounds;
            varnode->boundstatus = (node->boundstatus == SCIP_EXPRBOUNDSTATUS_VALID) ? SCIP_EXPRBOUNDSTATUS_VALID : SCIP_EXPRBOUNDSTATUS_CHILDRELAXED;
         }
      }
      else if( ncoefs == 0 )
      {
         /* turn node into EXPR_CONST node */

         /* remove variable (but don't free it's node) from graph */
         SCIP_CALL( exprgraphRemoveVar(exprgraph, varidx) );

         /* convert into EXPR_CONST node */
         varnode->op = SCIP_EXPR_CONST;
         varnode->data.dbl = constant;

         varnode->value = constant;
         SCIPintervalSet(&varnode->bounds, constant);
         varnode->boundstatus = SCIP_EXPRBOUNDSTATUS_VALID;

         /* add to constnodes arrays; @todo should use SCIPsortedvecInsertPtr? */
         ensureBlockMemoryArraySize(exprgraph->blkmem, &exprgraph->constnodes, &exprgraph->constssize, exprgraph->nconsts + 1);
         exprgraph->constnodes[exprgraph->nconsts] = varnode;
         ++exprgraph->nconsts;
         exprgraph->constssorted = exprgraph->nconsts <= 1 || (exprgraph->constssorted && exprgraphConstNodeComp(exprgraph->constnodes[exprgraph->nconsts-2], varnode) < 0);
      }
      else
      {
         /* turn node into EXPR_VARIDX node for new variable */

         /* remove variable (but don't free it's node) from graph */
         SCIP_CALL( exprgraphRemoveVar(exprgraph, varidx) );

         varnode->data.intval = exprgraph->nvars;

         /* add variable node to vars arrays and hashmap */
         ensureBlockMemoryArraySize3(exprgraph->blkmem, &exprgraph->vars, &exprgraph->varnodes, &exprgraph->varbounds, &exprgraph->varssize, exprgraph->nvars + 1);
         exprgraph->vars[exprgraph->nvars] = vars[0];  /*lint !e613*/
         exprgraph->varnodes[exprgraph->nvars] = varnode;
         SCIPintervalSetEntire(SCIP_REAL_MAX, &exprgraph->varbounds[exprgraph->nvars]);
         SCIP_CALL( SCIPhashmapInsert(exprgraph->varidxs, vars[0], (void*)(size_t)exprgraph->nvars) );  /*lint !e613*/
         ++exprgraph->nvars;

         /* call callback method, if set */
         if( exprgraph->exprgraphvaradded != NULL )
         {
            SCIP_CALL( exprgraph->exprgraphvaradded(exprgraph, exprgraph->userdata, vars[0], varnode) );  /*lint !e613*/
         }
      }

      /* mark varnode and its parents as not simplified */
      if( varnode != NULL )
      {
         varnode->simplified = FALSE;
         for( i = 0; i < varnode->nparents; ++i )
            varnode->parents[i]->simplified = FALSE;
      }

      return SCIP_OKAY;
   }

   /* turn varnode into EXPR_LINEAR */

   /* remove variable (but don't free it's node) from graph */
   SCIP_CALL( exprgraphRemoveVar(exprgraph, varidx) );

   /* move varnode up to depth 1 */
   SCIP_CALL( exprgraphMoveNode(exprgraph, varnode, 1) );

   /* convert into EXPR_LINEAR node */
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &lindata, ncoefs + 1) );
   BMScopyMemoryArray(lindata, coefs, ncoefs);  /*lint !e644*/
   lindata[ncoefs] = constant;
   varnode->data.data = (void*)lindata;
   varnode->op = SCIP_EXPR_LINEAR;

   /* add nodes corresponding to vars to expression graph, if not existing yet */
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &varnode->children, ncoefs) );
   SCIP_CALL( SCIPexprgraphAddVars(exprgraph, ncoefs, vars, varnode->children) );
   varnode->nchildren = ncoefs;

   /* notify vars about new parent varnode */
   for( i = 0; i < ncoefs; ++i )
   {
      SCIP_CALL( exprgraphNodeAddParent(exprgraph->blkmem, varnode->children[i], varnode) );
   }

   /* set value and bounds to invalid, curvature can remain (still linear) */
   varnode->value = SCIP_INVALID;
   varnode->boundstatus = SCIP_EXPRBOUNDSTATUS_CHILDRELAXED;

   /* mark varnode and its parents as not simplified */
   varnode->simplified = FALSE;
   for( i = 0; i < varnode->nparents; ++i )
      varnode->parents[i]->simplified = FALSE;

   return SCIP_OKAY;
}

/** finds expression graph node corresponding to a variable */
SCIP_Bool SCIPexprgraphFindVarNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   void*                 var,                /**< variable to search for */
   SCIP_EXPRGRAPHNODE**  varnode             /**< buffer to store node corresponding to variable, if found, or NULL if not found */
   )
{
   int pos;

   assert(exprgraph != NULL);
   assert(var != NULL);
   assert(varnode != NULL);

   if( !SCIPhashmapExists(exprgraph->varidxs, var) )
   {
      *varnode = NULL;
      return FALSE;
   }

   pos = (int)(size_t)SCIPhashmapGetImage(exprgraph->varidxs, var);
   assert(pos < exprgraph->nvars);

   *varnode = exprgraph->varnodes[pos];
   assert(*varnode != NULL);
   assert((*varnode)->op == SCIP_EXPR_VARIDX);

   return TRUE;
}

/** finds expression graph node corresponding to a constant */
SCIP_Bool SCIPexprgraphFindConstNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             constant,           /**< constant to search for */
   SCIP_EXPRGRAPHNODE**  constnode           /**< buffer to store node corresponding to constant, if found, or NULL if not found */
   )
{
   int left;
   int right;
   int middle;

   assert(exprgraph != NULL);
   assert(constnode != NULL);
   assert(constant == constant); /* cannot search for nan */  /*lint !e777*/

   exprgraphSortConstNodes(exprgraph);
   assert(exprgraph->constssorted);

   /* find node using binary search */
   left = 0;
   right = exprgraph->nconsts-1;
   *constnode = NULL;

   while( left <= right )
   {
      middle = (left+right)/2;
      assert(0 <= middle && middle < exprgraph->nconsts);

      if( constant < exprgraph->constnodes[middle]->data.dbl )
         right = middle - 1;
      else if( constant > exprgraph->constnodes[middle]->data.dbl )
         left  = middle + 1;
      else
      {
         *constnode = exprgraph->constnodes[middle];
         break;
      }
   }
   if( left == right+1 )
      return FALSE;

   assert(*constnode != NULL);
   assert((*constnode)->op == SCIP_EXPR_CONST);
   assert((*constnode)->data.dbl == constant);  /*lint !e777*/

   return TRUE;
}

/** prints an expression graph in dot format */
SCIP_RETCODE SCIPexprgraphPrintDot(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   const char**          varnames            /**< variable names, or NULL for generic names */
   )
{
   int d;
   int i;

   assert(exprgraph != NULL);

   if( file == NULL )
      file = stdout;

   SCIPmessageFPrintInfo(messagehdlr, file, "strict digraph exprgraph {\n");
   SCIPmessageFPrintInfo(messagehdlr, file, "node [fontcolor=white, style=filled, rankdir=LR]\n");

   for( d = 0; d < exprgraph->depth; ++d )
   {
      if( exprgraph->nnodes[d] == 0 )
         continue;

      for( i = 0; i < exprgraph->nnodes[d]; ++i )
      {
         exprgraphPrintNodeDot(exprgraph, exprgraph->nodes[d][i], messagehdlr, file, varnames);
      }
   }

   /* tell dot that all nodes of depth 0 have the same rank */
   SCIPmessageFPrintInfo(messagehdlr, file, "{rank=same;");
   for( i = 0; i < exprgraph->nnodes[0]; ++i )
      SCIPmessageFPrintInfo(messagehdlr, file, " n0_%d", i);
   SCIPmessageFPrintInfo(messagehdlr, file, "}\n");

   /* tell dot that all nodes without parent have the same rank */
   SCIPmessageFPrintInfo(messagehdlr, file, "{rank=same;");
   for( d = 0; d < exprgraph->depth; ++d )
      for( i = 0; i < exprgraph->nnodes[d]; ++i )
         if( exprgraph->nodes[d][i]->nparents == 0 )
            SCIPmessageFPrintInfo(messagehdlr, file, " n%d_%d", d, i);
   SCIPmessageFPrintInfo(messagehdlr, file, "}\n");

   SCIPmessageFPrintInfo(messagehdlr, file, "}\n");

   return SCIP_OKAY;
}

/** evaluates nodes of expression graph for given values of variables */
SCIP_RETCODE SCIPexprgraphEval(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real*            varvals             /**< values for variables */
   )
{
   int d;
   int i;

   assert(exprgraph != NULL);
   assert(varvals != NULL || exprgraph->nvars == 0);

   for( d = 0; d < exprgraph->depth; ++d )
      for( i = 0; i < exprgraph->nnodes[d]; ++i )
      {
         SCIP_CALL( exprgraphNodeEval(exprgraph->nodes[d][i], varvals) );
      }

   return SCIP_OKAY;
}

/** propagates bound changes in variables forward through the expression graph */
SCIP_RETCODE SCIPexprgraphPropagateVarBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Bool             clearreverseprop,   /**< whether to reset bound tightenings from reverse propagation */
   SCIP_Bool*            domainerror         /**< buffer to store whether a node with empty bounds has been found, propagation is interrupted in this case */
   )
{
   SCIP_EXPRGRAPHNODE* node;
   SCIP_Bool boundchanged;
   int d;
   int i;

   assert(exprgraph != NULL);
   assert(domainerror != NULL);

   *domainerror = FALSE;

   /* update bounds in varnodes of expression graph */
   exprgraphUpdateVarNodeBounds(exprgraph, &clearreverseprop, &boundchanged);

   /* if variable bounds have not changed and we do not have to clear a previous backward propagation, we can just return */
   if( !boundchanged && !clearreverseprop && !exprgraph->needvarboundprop )
   {
      SCIPdebugMessage("no bounds changed and clearreverseprop is FALSE -> skip propagation of variable bounds\n");
      return SCIP_OKAY;
   }

   /* propagate bound changes, interrupt if we get to a node with empty bounds */
   for( d = 1; d < exprgraph->depth; ++d )
   {
      for( i = 0; i < exprgraph->nnodes[d]; ++i )
      {
         node = exprgraph->nodes[d][i];
         SCIP_CALL( exprgraphNodeUpdateBounds(node, infinity, 1e-9, clearreverseprop) );
         if( SCIPintervalIsEmpty(infinity, node->bounds) )
         {
            SCIPdebugMessage("bounds of node %p(%d,%d) empty, stop bounds propagation\n", (void*)node, node->depth, node->pos);
            /* we keep exprgraph->needvarboundprop at TRUE, since we interrupt propagation */
            *domainerror = TRUE;
            return SCIP_OKAY;
         }
      }
   }

   exprgraph->needvarboundprop = FALSE;

   return SCIP_OKAY;
}

/** propagates bound changes in nodes backward through the graph
 *
 *  New bounds are not stored in varbounds, but only in nodes corresponding to variables.
 *  NOTE: it is assumed that SCIPexprgraphPropagateVarBounds was called before if variable bounds were relaxed.
 */
void SCIPexprgraphPropagateNodeBounds(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Real             minstrength,        /**< minimal required relative bound strengthening in a node to trigger a propagation into children nodes */
   SCIP_Bool*            cutoff              /**< buffer to store whether a node's bounds were propagated to an empty interval */
   )
{
   SCIP_EXPRGRAPHNODE* node;
   int d;
   int i;

   assert(exprgraph != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   for( d = exprgraph->depth-1; d >= 0 && !*cutoff; --d )
   {
      for( i = 0; i < exprgraph->nnodes[d] && !*cutoff; ++i )
      {
         node = exprgraph->nodes[d][i];
         exprgraphNodePropagateBounds(exprgraph, node, infinity, minstrength, cutoff);
      }
   }
   if( *cutoff )
      return;
}

/** updates curvature information in expression graph nodes w.r.t. currently stored variable bounds
 *
 *  Implies update of bounds in expression graph.
 */
SCIP_RETCODE SCIPexprgraphCheckCurvature(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_Bool             clearreverseprop    /**< whether to reset bound tightenings from reverse propagation */
   )
{
   SCIP_EXPRGRAPHNODE* node;
   SCIP_Bool boundchanged;
   int d;
   int i;

   assert(exprgraph != NULL);

   /* update bounds in varnodes of expression graph */
   exprgraphUpdateVarNodeBounds(exprgraph, &clearreverseprop, &boundchanged);

#ifndef NDEBUG
   for( i = 0; i < exprgraph->nnodes[0]; ++i )
      assert(exprgraph->nodes[0][i]->curv == SCIP_EXPRCURV_LINEAR);
#endif

   for( d = 1; d < exprgraph->depth; ++d )
      for( i = 0; i < exprgraph->nnodes[d]; ++i )
      {
         node = exprgraph->nodes[d][i];
         assert(node != NULL);

         SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, infinity, 1e-9, clearreverseprop) );

         if( SCIPintervalIsEmpty(infinity, node->bounds) )
         {
            SCIPerrorMessage("SCIPexprgraphCheckCurvature gets domain error while propagating variables bounds, ignoring...\n");
            return SCIP_OKAY;
         }
      }

   return SCIP_OKAY;
}

/** aims at simplifying an expression graph
 *
 *  A domain error can occur when variables were fixed to values for which a parent expression is not defined (e.g., 0^(-1) or log(-1)).
 */
SCIP_RETCODE SCIPexprgraphSimplify(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Real             eps,                /**< threshold, under which positive values are treat as 0 */
   int                   maxexpansionexponent,/**< maximal exponent for which we still expand non-monomial polynomials */
   SCIP_Bool*            havechange,         /**< buffer to indicate whether the graph has been modified */
   SCIP_Bool*            domainerror         /**< buffer to indicate whether a domain error has been encountered, i.e., some expressions turned into NaN */
   )
{
   SCIP_EXPRGRAPHNODE* node;
   SCIP_Bool havechangenode;
   SCIP_Bool allsimplified;
   int d;
   int i;
   int j;

#ifndef NDEBUG
   SCIP_Real* testx;
   SCIP_HASHMAP* testvalidx;
   SCIP_Real* testvals;
   SCIP_RANDNUMGEN* randnumgen;
   int testvalssize;
   int ntestvals;
#endif

   assert(exprgraph != NULL);
   assert(eps >= 0.0);
   assert(havechange != NULL);
   assert(domainerror != NULL);

#ifndef NDEBUG
   SCIP_CALL( SCIPrandomCreate(&randnumgen, exprgraph->blkmem, 862) ); /* see also #1848 */
   SCIP_CALL( SCIPhashmapCreate(&testvalidx, exprgraph->blkmem, 1000) );
   testvals = NULL;
   ntestvals = 0;
   testvalssize = 0;

   SCIP_ALLOC( BMSallocMemoryArray(&testx, exprgraph->nvars) );
   for( i = 0; i < exprgraph->nvars; ++i )
      testx[i] = SCIPrandomGetReal(randnumgen, -100.0, 100.0);  /*lint !e644*/
   SCIP_CALL( SCIPexprgraphEval(exprgraph, testx) );
   for( d = 1; d < exprgraph->depth; ++d )
      for( i = 0; i < exprgraph->nnodes[d]; ++i )
      {
         node = exprgraph->nodes[d][i];
         assert(node != NULL);

         /* nodes that are in use should not be removed by simplifier, so for those we store their value and check if it remains the same after simplifier was run */
         if( node->nuses > 0 )
         {
            ensureBlockMemoryArraySize(exprgraph->blkmem, &testvals, &testvalssize, ntestvals+1);
            SCIP_CALL( SCIPhashmapInsert(testvalidx, (void*)node, (void*)(size_t)ntestvals) );
            testvals[ntestvals] = SCIPexprgraphGetNodeVal(node);  /*lint !e613 !e794*/
            ++ntestvals;
         }
      }

   SCIPrandomFree(&randnumgen, exprgraph->blkmem);
#endif

#ifdef SCIP_OUTPUT
   {
      FILE* file;
      file = fopen("exprgraph_beforesimplify.dot", "w");
      if( file != NULL )
      {
         SCIP_CALL( SCIPexprgraphPrintDot(exprgraph, messagehdlr, file, NULL) );
         fclose(file);
      }
   }
#endif

   *havechange = FALSE;  /* we have not changed any node yet */
   *domainerror = FALSE; /* no domain errors encountered so far */
   allsimplified = TRUE; /* all nodes we looked at are simplified */

   /* call node simplifier from bottom up
    * for each node, convert to polynomials and merge in child nodes that are not in use otherwise, if possible
    */
   for( d = 1; d < exprgraph->depth && !*domainerror; ++d )
   {
      for( i = 0; i < exprgraph->nnodes[d]; ++i )
      {
         node = exprgraph->nodes[d][i];
         assert(node != NULL);

         havechangenode = FALSE; /* node did not change yet */

         if( node->op != SCIP_EXPR_CONST )
         {
            /* skip nodes that are already simplified */
            if( node->simplified )
               continue;

            allsimplified = FALSE;  /* looks like we found a node that has not been simplified */

            /* we should be careful about declaring numbers close to zero as zero, so take eps^2 as tolerance */
            SCIP_CALL( exprgraphNodeSimplify(exprgraph, node, messagehdlr, eps*eps, maxexpansionexponent, &havechangenode) );
            assert(node->simplified == TRUE);
            *havechange |= havechangenode;
         }

         /* if node was or has been converted into constant, may move to depth 0 */
         if( node->op == SCIP_EXPR_CONST )
         {
            SCIP_EXPRGRAPHNODE* constnode;

            if( !SCIPisFinite(node->value) )  /*lint !e777*/
            {
               SCIPdebugMessage("Expression graph simplify turned node into NaN or inf.\n");
               *domainerror = TRUE;
               break;
            }

            /* check if there is already a node for this constant */
            if( SCIPexprgraphFindConstNode(exprgraph, node->value, &constnode) )
            {
               assert(constnode->op == SCIP_EXPR_CONST);
               assert(constnode->data.dbl == node->value);  /*lint !e777*/

               if( node->nparents > 0 )
               {
                  /* move parents of this node to constnode, node may be freed if not in use */
                  SCIP_CALL( SCIPexprgraphMoveNodeParents(exprgraph, &node, constnode) );
                  /* node should have no parents anymore, so it should have been freed if not in use */
                  assert(node == NULL || node->nuses > 0);
                  havechangenode = TRUE;

                  /* if node was freed, exprgraph->nodes[i] points to the next node that need to be simplified */
                  if( node == NULL )
                  {
                     --i;
                     continue;
                  }
               }
               assert(node != NULL);
               assert(node->nuses > 0);

               if( constnode->nuses == 0 )
               {
                  /* move node to depth 0, adding it to constnodes */
                  SCIP_CALL( exprgraphMoveNode(exprgraph, node, 0) );

                  /* move parents of constnode to node, so constnode is freed */
                  SCIP_CALL( SCIPexprgraphMoveNodeParents(exprgraph, &constnode, node) );
                  assert(constnode == NULL);
                  havechangenode = TRUE;

                  /* node moved to depth 0, so exprgraph->nodes[i] points to the next node that need to be simplified */
                  --i;
                  continue;
               }
            }
            else
            {
               /* move to depth 0, adding it to constnodes */
               SCIP_CALL( exprgraphMoveNode(exprgraph, node, 0) );

               /* node moved to depth 0, so exprgraph->nodes[i] points to the next node that need to be simplified */
               --i;
            }
         }

         /* if there was a change, mark parents as not simplified */
         if( havechangenode )
            for( j = 0; j < node->nparents; ++j )
               node->parents[j]->simplified = FALSE;
      }
   }  /*lint !e850*/

   /* if we did nothing, clean up and escape from here */
   if( allsimplified || *domainerror )
      goto EXPRGRAPHSIMPLIFY_CLEANUP;

   /* @todo find duplicate subexpressions in expression graph */

   /* unconvert polynomials into simpler expressions, where possible */
   for( d = 1; d < exprgraph->depth; ++d )
   {
      for( i = 0; i < exprgraph->nnodes[d]; ++i )
      {
         node = exprgraph->nodes[d][i];
         assert(node != NULL);

         if( node->op != SCIP_EXPR_POLYNOMIAL )
            continue;

         SCIP_CALL( exprUnconvertPolynomial(exprgraph->blkmem, &node->op, &node->data, node->nchildren, (void**)node->children) );

         if( node->op == SCIP_EXPR_SUM && node->nchildren == 1 )
         {
            /* node is identity w.r.t only child
             * replace node as child of parents by child of node
             */

            for( j = 0; node != NULL && j < node->nparents; ++j )
            {
               SCIP_CALL( exprgraphNodeReplaceChild(exprgraph, node->parents[j], &node, node->children[0]) );
            }
            /* node should have no parents anymore, so it should have been freed if not in use */
            assert(node == NULL || node->nuses > 0);

            /* if node was freed, exprgraph->nodes[i] points to the next node that need to be unconverted */
            if( node == NULL )
               --i;
         }
      }
   }  /*lint !e850*/

#ifdef SCIP_OUTPUT
   {
      FILE* file;
      file = fopen("exprgraph_aftersimplify.dot", "w");
      if( file != NULL )
      {
         SCIP_CALL( SCIPexprgraphPrintDot(exprgraph, messagehdlr, file, NULL) );
         fclose(file);
      }
   }
#endif

#ifndef NDEBUG
   for( d = 1; d < exprgraph->depth; ++d )
      for( i = 0; i < exprgraph->nnodes[d]; ++i )
      {
         int idx;
         SCIP_Real testval_before;
         SCIP_Real testval_after;

         node = exprgraph->nodes[d][i];
         assert(node != NULL);

         SCIP_CALL( exprgraphNodeEval(node, NULL) );

         /* nodes that are in use should not have been removed by simplifier, check if they still have the same value in our testpoint */
         if( node->nuses > 0 )
         {
            assert(SCIPhashmapExists(testvalidx, (void*)node));

            idx = (int)(size_t)(void*)SCIPhashmapGetImage(testvalidx, (void*)node);
            assert(idx < ntestvals);

            testval_before = testvals[idx];  /*lint !e613*/
            testval_after = SCIPexprgraphGetNodeVal(node);

            assert(!SCIPisFinite(testval_before) || EPSZ(SCIPrelDiff(testval_before, testval_after), eps));  /*lint !e777*/
         }
      }
#endif

 EXPRGRAPHSIMPLIFY_CLEANUP:
#ifndef NDEBUG
   BMSfreeMemoryArray(&testx);
   BMSfreeBlockMemoryArrayNull(exprgraph->blkmem, &testvals, testvalssize);
   SCIPhashmapFree(&testvalidx);
#endif

   return SCIP_OKAY;
}

/** creates an expression tree from a given node in an expression graph */
SCIP_RETCODE SCIPexprgraphGetTree(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   rootnode,           /**< expression graph node that should represent root of expression tree */
   SCIP_EXPRTREE**       exprtree            /**< buffer to store pointer to created expression tree */
   )
{
   SCIP_EXPR* root;
   int nexprvars;
   int* varidx;
   int i;

   assert(exprgraph != NULL);
   assert(rootnode  != NULL);
   assert(rootnode->depth >= 0);
   assert(rootnode->pos >= 0);
   assert(exprtree != NULL);

   /* buffer where to store mapping of expression graph variable indices to expression tree variable indices */
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &varidx, exprgraph->nvars) );

   /* initially, no variable appears in the expression tree */
   for( i = 0; i < exprgraph->nvars; ++i )
      varidx[i] = -1;  /*lint !e644*/
   nexprvars = 0;

   /* create expression from the subgraph that has rootnode as root */
   SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, rootnode, &root, &nexprvars, varidx) );

   /* create expression tree for this expression */
   SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, exprtree, root, nexprvars, 0, NULL) );

   /* copy variables into expression tree */
   if( nexprvars > 0 )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &(*exprtree)->vars, nexprvars) );
      for( i = 0; i < exprgraph->nvars; ++i )
      {
         assert(varidx[i] >= -1);
         assert(varidx[i] < nexprvars);
         if( varidx[i] >= 0 )
            (*exprtree)->vars[varidx[i]] = exprgraph->vars[i];
      }
   }

   BMSfreeBlockMemoryArray(exprgraph->blkmem, &varidx, exprgraph->nvars);

   return SCIP_OKAY;
}

/** creates a sum of expression trees with pairwise disjoint variables from a given node in an expression graph
 *
 *  Giving SCIPexprgraphGetNodeNChildren() for exprtreesize is always sufficient.
 */
SCIP_RETCODE SCIPexprgraphGetSeparableTrees(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node which represents expression to get */
   int                   exprtreessize,      /**< length of exprtrees and exprtreecoefs arrays, need to be at least one */
   int*                  nexprtrees,         /**< buffer to store number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< array where to store expression trees */
   SCIP_Real*            exprtreecoefs       /**< array where to store coefficients of expression trees */
   )
{
   int ncomponents;
   int* childcomp;
   int* varcomp;
   int compnr;
   SCIP_Bool haveoverlap;
   int i;
   int j;
   int k;

   SCIP_EXPR** exprs;
   int nexprs;
   int* childmap;
   int* childmapinv;
   int* varidx;
   int nexprvars;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->depth >= 0);
   assert(node->pos >= 0);
   assert(exprtreessize > 0);
   assert(nexprtrees != NULL);
   assert(exprtrees != NULL);
   assert(exprtreecoefs != NULL);

   /* easy cases: if we have space for only one tree or there is only one child or only one variable in the graph,
    * or the node operator is not separable, fallback to SCIPexprgraphGetTree */
   if( exprtreessize == 1 || node->nchildren <= 1 || exprgraph->nvars <= 1 ||
      (  node->op  != SCIP_EXPR_PLUS   &&
         node->op  != SCIP_EXPR_MINUS  &&
         node->op  != SCIP_EXPR_SUM    &&
         node->op  != SCIP_EXPR_LINEAR &&
         (node->op != SCIP_EXPR_QUADRATIC || ((SCIP_EXPRDATA_QUADRATIC*)node->data.data)->nquadelems <= 1) &&
         (node->op != SCIP_EXPR_POLYNOMIAL || ((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data)->nmonomials <= 1)) )
   {
      *nexprtrees = 1;
      exprtreecoefs[0] = 1.0;
      SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node, exprtrees) );

      return SCIP_OKAY;
   }

   /* find components in node->children <-> variables graph */
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childcomp, node->nchildren) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &varcomp, exprgraph->nvars) );
   for( i = 0; i < exprgraph->nvars; ++i )
      varcomp[i] = -1;  /*lint !e644*/

   haveoverlap = FALSE;
   for( i = 0; i < node->nchildren; ++i )
   {
      compnr = i;
      exprgraphNodeCheckSeparabilityComponent(node->children[i], &compnr, i-1, childcomp, exprgraph->nvars, varcomp);  /*lint !e644*/
      assert(compnr >= 0);
      assert(compnr < node->nchildren);
      childcomp[i] = compnr;

      /* remember if component number was changed by CheckComponent */
      if( compnr != i )
         haveoverlap = TRUE;
   }

   BMSfreeBlockMemoryArray(exprgraph->blkmem, &varcomp, exprgraph->nvars);

   if( node->op == SCIP_EXPR_QUADRATIC )
   {
      /* merge components for products of children from different components */
      SCIP_EXPRDATA_QUADRATIC* data;

      data = (SCIP_EXPRDATA_QUADRATIC*)node->data.data;
      assert(data != NULL);

      for( i = 0; i < data->nquadelems; ++i )
         if( childcomp[data->quadelems[i].idx1] != childcomp[data->quadelems[i].idx2] )
         {
            /* reassign all children in component childcomp[data->quadelems[i].idx2] to component childcomp[data->quadelems[i].idx1] */
            compnr = childcomp[data->quadelems[i].idx2];
            for( j = 0; j < node->nchildren; ++j )
               if( childcomp[j] == compnr )
                  childcomp[j] = childcomp[data->quadelems[i].idx1];
            assert(childcomp[data->quadelems[i].idx1] == childcomp[data->quadelems[i].idx2]);
            haveoverlap = TRUE;
         }
   }
   else if( node->op == SCIP_EXPR_POLYNOMIAL )
   {
      /* merge components for monomials of children from different components */
      SCIP_EXPRDATA_POLYNOMIAL* data;

      data = (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data;
      assert(data != NULL);

      for( i = 0; i < data->nmonomials; ++i )
         for( j = 1; j < data->monomials[i]->nfactors; ++j )
            if( childcomp[data->monomials[i]->childidxs[j]] != childcomp[data->monomials[i]->childidxs[0]] )
            {
               /* reassign all children in component childcomp[data->monomials[i]->childidxs[j]] to component childcomp[data->monomials[i]->childidxs[0]] */
               compnr = childcomp[data->monomials[i]->childidxs[j]];
               for( k = 0; k < node->nchildren; ++k )
                  if( childcomp[k] == compnr )
                     childcomp[k] = childcomp[data->monomials[i]->childidxs[0]];
               assert(childcomp[data->monomials[i]->childidxs[j]] == childcomp[data->monomials[i]->childidxs[0]]);
               haveoverlap = TRUE;
            }
   }

   if( haveoverlap )
   {
      /* some component numbers are unused, thus relabel and count final number of components */
      int* compmap;

      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &compmap, node->nchildren) );
      for( i = 0; i < node->nchildren; ++i )
         compmap[i] = -1;  /*lint !e644*/

      ncomponents = 0;
      for( i = 0; i < node->nchildren; ++i )
      {
         if( compmap[childcomp[i]] == -1 )
            compmap[childcomp[i]] = ncomponents++;
         childcomp[i] = compmap[childcomp[i]];
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &compmap, node->nchildren);
   }
   else
   {
      ncomponents = node->nchildren;
   }

   if( ncomponents == 1 )
   {
      /* it turned out that expression is not block separable, so fallback to SCIPexprgraphGetTree */
      BMSfreeBlockMemoryArray(exprgraph->blkmem, &childcomp, node->nchildren);

      *nexprtrees = 1;
      exprtreecoefs[0] = 1.0;
      SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node, exprtrees) );

      return SCIP_OKAY;
   }

   if( ncomponents > exprtreessize )
   {
      /* if we have not enough space for all expressions, merge components with number > exprtreessize into component exprtreessize */
      for( i = 0; i < node->nchildren; ++i )
         if( childcomp[i] >= exprtreessize )
            childcomp[i] = exprtreessize-1;
      ncomponents = exprtreessize;
   }

   assert(ncomponents >= 2);

   /* setup expression trees for each component */
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &exprs, node->nchildren) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &varidx, exprgraph->nvars) ); /* mapping of expression graph variable indices to expression tree variable indices */
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childmap, node->nchildren) ); /* mapping of child indices from node to expressions belonging to a single component */
   SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childmapinv, node->nchildren) ); /* mapping of child indices from expressions belonging to a single component to node */
   for( i = 0; i < ncomponents; ++i )
   {
      /* initially, no variable appears in the expression tree */
      for( j = 0; j < exprgraph->nvars; ++j )
         varidx[j] = -1;  /*lint !e644*/
      nexprvars = 0;

      /* collect expressions from children belonging to component i */
      nexprs = 0;
      for( j = 0; j < node->nchildren; ++j )
      {
         assert(childcomp[j] >= 0);
         assert(childcomp[j] < ncomponents);
         if( childcomp[j] != i )
            continue;

         /* create expression from the subgraph that has child j as root */
         SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, node->children[j], &exprs[nexprs], &nexprvars, varidx) );  /*lint !e644*/
         childmap[j] = nexprs;     /*lint !e644*/
         childmapinv[nexprs] = j;  /*lint !e644*/
         ++nexprs;
      }

      /* setup expression tree for component i */
      switch( node->op )
      {
      case SCIP_EXPR_PLUS:
      {
         assert(ncomponents == 2);
         assert(nexprs == 1);

         SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], exprs[0], nexprvars, 0, NULL) );
         exprtreecoefs[i] = 1.0;

         break;
      }

      case SCIP_EXPR_MINUS:
      {
         assert(ncomponents == 2);
         assert(nexprs == 1);

         SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], exprs[0], nexprvars, 0, NULL) );
         /* if component i consists of first child, then it has coefficient 1.0, otherwise it has coefficient -1 */
         assert(childmapinv[0] == 0 || childmapinv[0] == 1);
         exprtreecoefs[i] = (childmapinv[0] == 0 ? 1.0 : -1.0);

         break;
      }

      case SCIP_EXPR_SUM:
      {
         if( nexprs == 1 )
         {
            /* component corresponds to exactly one child of node */
            SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], exprs[0], nexprvars, 0, NULL) );
         }
         else
         {
            /* component corresponds to a sum of children of node */
            SCIP_EXPR* sumexpr;

            SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &sumexpr, SCIP_EXPR_SUM, nexprs, exprs) );
            SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], sumexpr, nexprvars, 0, NULL) );
         }
         exprtreecoefs[i] = 1.0;

         break;
      }

      case SCIP_EXPR_LINEAR:
      {
         SCIP_Real* nodecoefs;
         SCIP_EXPR* sumexpr;

         nodecoefs = (SCIP_Real*)node->data.data;

         /* if there is a constant, then we put it into the expression of the first component */
         if( nexprs == 1 && (i > 0 || nodecoefs[node->nchildren] == 0.0) )
         {
            /* component corresponds to exactly one child of node */
            SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], exprs[0], nexprvars, 0, NULL) );
            exprtreecoefs[i] = nodecoefs[childmapinv[0]];
         }
         else if( nexprs == 1 )
         {
            /* component corresponds to a sum of one child and a constant */
            assert(i == 0);
            assert(nodecoefs[node->nchildren] != 0.0);
            assert(nodecoefs[childmapinv[0]] != 0.0);
            SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &sumexpr, SCIP_EXPR_CONST, nodecoefs[node->nchildren]/nodecoefs[childmapinv[0]]) );
            SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &sumexpr, SCIP_EXPR_PLUS, sumexpr, exprs[0]) );
            SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], sumexpr, nexprvars, 0, NULL) );
            exprtreecoefs[i] = nodecoefs[childmapinv[0]];
         }
         else
         {
            /* component corresponds to a linear combination of children of node */

            if( nexprs == 2 && nodecoefs[childmapinv[0]] == nodecoefs[childmapinv[1]] && (i > 0 || nodecoefs[node->nchildren] == 0.0) )  /*lint !e777*/
            {
               /* if two expressions with equal sign, then create PLUS expression */
               SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &sumexpr, SCIP_EXPR_PLUS, exprs[0], exprs[1]) );
               exprtreecoefs[i] = nodecoefs[childmapinv[0]];
            }
            else if( nexprs == 2 && nodecoefs[childmapinv[0]] == -nodecoefs[childmapinv[1]] && (i > 0 || nodecoefs[node->nchildren] == 0.0) )  /*lint !e777*/
            {
               /* if two expressions with opposite sign, then create MINUS expression */
               SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &sumexpr, SCIP_EXPR_MINUS, exprs[0], exprs[1]) );
               exprtreecoefs[i] = nodecoefs[childmapinv[0]];
            }
            else
            {
               /* assemble coefficents and create SUM or LINEAR expression */
               SCIP_Real* coefs;
               SCIP_Bool allcoefsequal;

               SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &coefs, nexprs) );
               allcoefsequal = TRUE;
               coefs[0] = nodecoefs[childmapinv[0]];  /*lint !e644*/
               for( j = 0; j < nexprs; ++j )
               {
                  coefs[j] = nodecoefs[childmapinv[j]];
                  allcoefsequal &= (coefs[j] == coefs[0]);  /*lint !e777 !e514*/
               }

               /* if all coefficients are equal and no constant, create SUM expression, otherwise LINEAR expression */
               if( allcoefsequal && (i > 0 || nodecoefs[node->nchildren] == 0.0) )
               {
                  SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &sumexpr, SCIP_EXPR_SUM, nexprs, exprs) );
                  exprtreecoefs[i] = coefs[0];
               }
               else
               {
                  SCIP_CALL( SCIPexprCreateLinear(exprgraph->blkmem, &sumexpr, nexprs, exprs, coefs, i == 0 ? nodecoefs[node->nchildren] : 0.0) );
                  exprtreecoefs[i] = 1.0;
               }

               BMSfreeBlockMemoryArray(exprgraph->blkmem, &coefs, nexprs);
            }

            SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], sumexpr, nexprvars, 0, NULL) );
         }

         break;
      }

      case SCIP_EXPR_QUADRATIC:
      {
         SCIP_EXPR* quadexpr;
         SCIP_EXPRDATA_QUADRATIC* nodedata;
         SCIP_Real* lincoefs;
         SCIP_QUADELEM* quadelems;
         int nquadelems;

         nodedata = (SCIP_EXPRDATA_QUADRATIC*)node->data.data;

         exprtreecoefs[i] = 1.0;

         /* assemble coefficients corresponding to component i */
         if( nodedata->lincoefs != NULL )
         {
            SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &lincoefs, nexprs) );
            for( j = 0; j < nexprs; ++j )
               lincoefs[j] = nodedata->lincoefs[childmapinv[j]];  /*lint !e771*/
         }
         else
            lincoefs = NULL;

         SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &quadelems, nodedata->nquadelems) );
         nquadelems = 0;
         for( j = 0; j < nodedata->nquadelems; ++j )
         {
            assert(childcomp[nodedata->quadelems[j].idx1] == childcomp[nodedata->quadelems[j].idx2]);
            if( childcomp[nodedata->quadelems[j].idx1] != i )
               continue;
            quadelems[nquadelems].idx1 = MIN(childmap[nodedata->quadelems[j].idx1], childmap[nodedata->quadelems[j].idx2]);  /*lint !e644*/
            quadelems[nquadelems].idx2 = MAX(childmap[nodedata->quadelems[j].idx1], childmap[nodedata->quadelems[j].idx2]);
            quadelems[nquadelems].coef = nodedata->quadelems[j].coef;
            ++nquadelems;
         }

         /* put constant into first component */
         SCIP_CALL( SCIPexprCreateQuadratic(exprgraph->blkmem, &quadexpr, nexprs, exprs, i == 0 ? nodedata->constant : 0.0, lincoefs, nquadelems, quadelems) );
         SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], quadexpr, nexprvars, 0, NULL) );

         BMSfreeBlockMemoryArray(exprgraph->blkmem, &quadelems, nodedata->nquadelems);
         BMSfreeBlockMemoryArrayNull(exprgraph->blkmem, &lincoefs, nexprs);

         break;
      }

      case SCIP_EXPR_POLYNOMIAL:
      {
         SCIP_EXPR* polyexpr;
         SCIP_EXPRDATA_POLYNOMIAL* nodedata;
         SCIP_EXPRDATA_MONOMIAL** monomials;
         SCIP_Real constant;
         int nmonomials;

         nodedata = (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data;

         constant = nodedata->constant;
         exprtreecoefs[i] = 1.0;

         /* collect monomials belonging to component i */
         SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &monomials, nodedata->nmonomials) );
         nmonomials = 0;
         for( j = 0; j < nodedata->nmonomials; ++j )
         {
            if( nodedata->monomials[j]->nfactors == 0 )
            {
               constant += nodedata->monomials[j]->coef;
               continue;
            }
            if( childcomp[nodedata->monomials[j]->childidxs[0]] != i )
               continue;

            SCIP_CALL( SCIPexprCreateMonomial(exprgraph->blkmem, &monomials[nmonomials], nodedata->monomials[j]->coef, nodedata->monomials[j]->nfactors,
                  nodedata->monomials[j]->childidxs, nodedata->monomials[j]->exponents) );  /*lint !e644*/
            for( k = 0; k < monomials[nmonomials]->nfactors; ++k )
            {
               assert(childcomp[nodedata->monomials[j]->childidxs[k]] == i);
               monomials[nmonomials]->childidxs[k] = childmap[monomials[nmonomials]->childidxs[k]];
            }
            ++nmonomials;
         }

         SCIP_CALL( SCIPexprCreatePolynomial(exprgraph->blkmem, &polyexpr, nexprs, exprs, nmonomials, monomials, i == 0 ? constant : 0.0, FALSE) );
         SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[i], polyexpr, nexprvars, 0, NULL) );

         BMSfreeBlockMemoryArray(exprgraph->blkmem, &monomials, nodedata->nmonomials);

         break;
      }

      default:
         SCIPerrorMessage("unexpected operator type %d\n", node->op);
         return SCIP_ERROR;
      }  /*lint !e788*/

      /* copy variables into expression tree */
      if( nexprvars > 0 )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &exprtrees[i]->vars, nexprvars) );
         for( j = 0; j < exprgraph->nvars; ++j )
         {
            assert(varidx[j] >= -1);
            assert(varidx[j] < nexprvars);
            if( varidx[j] >= 0 )
               exprtrees[i]->vars[varidx[j]] = exprgraph->vars[j];
         }
      }
   }

   BMSfreeBlockMemoryArray(exprgraph->blkmem, &exprs, node->nchildren);
   BMSfreeBlockMemoryArray(exprgraph->blkmem, &varidx, exprgraph->nvars);
   BMSfreeBlockMemoryArray(exprgraph->blkmem, &childmap, node->nchildren);
   BMSfreeBlockMemoryArray(exprgraph->blkmem, &childmapinv, node->nchildren);
   BMSfreeBlockMemoryArray(exprgraph->blkmem, &childcomp, node->nchildren);

   *nexprtrees = ncomponents;

   return SCIP_OKAY;
}

/** returns how often expression graph variables are used in a subtree of the expression graph */
void SCIPexprgraphGetSubtreeVarsUsage(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< root node of expression graph subtree */
   int*                  varsusage           /**< array where to count usage of variables, length must be at least the number of variables in the graph */
   )
{
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(varsusage != NULL);

   BMSclearMemoryArray(varsusage, exprgraph->nvars);

   exprgraphNodeGetVarsUsage(node, varsusage);
}

/** gives the number of summands which the expression of an expression graph node consists of */
int SCIPexprgraphGetSumTreesNSummands(
   SCIP_EXPRGRAPHNODE*   node                /**< expression graph node */
   )
{
   switch( node->op )
   {
   case SCIP_EXPR_PLUS:
   case SCIP_EXPR_MINUS:
      return 2;

   case SCIP_EXPR_SUM:
   case SCIP_EXPR_LINEAR:
      return node->nchildren;

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* nodedata;

      nodedata = (SCIP_EXPRDATA_QUADRATIC*)node->data.data;
      return (nodedata->lincoefs != NULL ? node->nchildren : 0) + nodedata->nquadelems;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* nodedata;

      nodedata = (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data;
      return nodedata->nmonomials;
   }

   default:
      return 1;
   }  /*lint !e788*/
}

/** creates a sum of expression trees, possibly sharing variables, from a given node in an expression graph */
SCIP_RETCODE SCIPexprgraphGetSumTrees(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node which represents expression to get */
   int                   exprtreessize,      /**< length of exprtrees and exptreecoefs arrays, should be at least SCIPexprgraphGetSumTreesNSummands() */
   int*                  nexprtrees,         /**< buffer to store number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< array where to store expression trees */
   SCIP_Real*            exprtreecoefs       /**< array where to store coefficients of expression trees */
   )
{
   int* varidx;
   int nexprvars;
   int i;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(node->depth >= 0);
   assert(node->pos >= 0);
   assert(exprtreessize > 0);
   assert(nexprtrees != NULL);
   assert(exprtrees != NULL);
   assert(exprtreecoefs != NULL);

   /* if node is not separable, fallback to SCIPexprgraphGetTree */
   if( node->op != SCIP_EXPR_PLUS  &&
      node->op  != SCIP_EXPR_MINUS &&
      node->op  != SCIP_EXPR_SUM   &&
      (node->op != SCIP_EXPR_LINEAR || node->nchildren <= 1) &&
      (node->op != SCIP_EXPR_QUADRATIC || (((SCIP_EXPRDATA_QUADRATIC*)node->data.data)->lincoefs == NULL && ((SCIP_EXPRDATA_QUADRATIC*)node->data.data)->nquadelems <= 1)) &&
      (node->op != SCIP_EXPR_POLYNOMIAL || ((SCIP_EXPRDATA_POLYNOMIAL*)node->data.data)->nmonomials <= 1) )
   {
      *nexprtrees = 1;
      exprtreecoefs[0] = 1.0;
      SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node, exprtrees) );

      return SCIP_OKAY;
   }

   switch( node->op )
   {
   case SCIP_EXPR_PLUS:
   {
      assert(exprtreessize >= 2);

      SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node->children[0], &exprtrees[0]) );
      SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node->children[1], &exprtrees[1]) );

      exprtreecoefs[0] = 1.0;
      exprtreecoefs[1] = 1.0;

      *nexprtrees = 2;
      break;
   }

   case SCIP_EXPR_MINUS:
   {
      assert(exprtreessize >= 2);

      SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node->children[0], &exprtrees[0]) );
      SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node->children[1], &exprtrees[1]) );

      exprtreecoefs[0] =  1.0;
      exprtreecoefs[1] = -1.0;

      *nexprtrees = 2;
      break;
   }

   case SCIP_EXPR_SUM:
   {
      assert(exprtreessize >= node->nchildren);

      for( i = 0; i < node->nchildren; ++i )
      {
         SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node->children[i], &exprtrees[i]) );
         exprtreecoefs[i] = 1.0;
      }

      *nexprtrees = node->nchildren;
      break;
   }

   case SCIP_EXPR_LINEAR:
   {
      SCIP_Real* nodecoefs;

      assert(exprtreessize >= node->nchildren);
      assert(node->nchildren > 0);

      nodecoefs = (SCIP_Real*)node->data.data;
      assert(nodecoefs != NULL);

      for( i = 0; i < node->nchildren; ++i )
      {
         SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node->children[i], &exprtrees[i]) );
         exprtreecoefs[i] = nodecoefs[i];
      }

      /* add constant to first summand, if nonzero; need to divide by coef of this exprtree */
      if( nodecoefs[node->nchildren] != 0.0 )
      {
         SCIP_EXPR* constexpr_;

         SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &constexpr_, SCIP_EXPR_CONST, nodecoefs[node->nchildren] / exprtreecoefs[0]) );
         SCIP_CALL( SCIPexprtreeAddExpr(exprtrees[0], constexpr_, FALSE) );
      }

      *nexprtrees = node->nchildren;
      break;
   }

   case SCIP_EXPR_QUADRATIC:
   {
      SCIP_EXPRDATA_QUADRATIC* nodedata;
      SCIP_Real* lincoefs;
      SCIP_QUADELEM* quadelems;
      int nquadelems;
      SCIP_EXPR* expr;
      int j;

      nodedata = (SCIP_EXPRDATA_QUADRATIC*)node->data.data;
      lincoefs = nodedata->lincoefs;
      quadelems = nodedata->quadelems;
      nquadelems = nodedata->nquadelems;

      assert(exprtreessize >= (lincoefs != NULL ? node->nchildren : 0) + nquadelems);
      assert(node->nchildren > 0);

      *nexprtrees = 0;
      if( lincoefs != NULL )
      {
         for( i = 0; i < node->nchildren; ++i )
         {
            if( lincoefs[i] == 0.0 )
               continue;
            SCIP_CALL( SCIPexprgraphGetTree(exprgraph, node->children[i], &exprtrees[*nexprtrees]) );
            exprtreecoefs[*nexprtrees] = lincoefs[i];
            ++*nexprtrees;
         }
      }

      /* buffer where to store mapping of expression graph variable indices to expression tree variable indices */
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &varidx, exprgraph->nvars) );

      for( i = 0; i < nquadelems; ++i )
      {
         /* initially, no variable appears in the expression tree */
         for( j = 0; j < exprgraph->nvars; ++j )
            varidx[j] = -1;  /*lint !e644*/
         nexprvars = 0;

         /* create expression from the subgraph at quadelems[i].idx1 */
         SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, node->children[quadelems[i].idx1], &expr, &nexprvars, varidx) );

         if( quadelems[i].idx1 == quadelems[i].idx2 )
         {
            /* create expression for square of expr */
            SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &expr, SCIP_EXPR_SQUARE, expr) );
         }
         else
         {
            SCIP_EXPR* expr2;

            /* create expression from the subgraph at quadelems[i].idx2, may add more variables into varidx */
            SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, node->children[quadelems[i].idx2], &expr2, &nexprvars, varidx) );
            /* create expression for product */
            SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &expr, SCIP_EXPR_MUL, expr, expr2) );
         }

         /* create expression tree for expr */
         SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[*nexprtrees], expr, nexprvars, 0, NULL) );

         /* copy variables into expression tree */
         if( nexprvars > 0 )
         {
            SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &exprtrees[*nexprtrees]->vars, nexprvars) );
            for( j = 0; j < exprgraph->nvars; ++j )
            {
               assert(varidx[j] >= -1);
               assert(varidx[j] < nexprvars);
               if( varidx[j] >= 0 )
                  exprtrees[*nexprtrees]->vars[varidx[j]] = exprgraph->vars[j];
            }
         }

         exprtreecoefs[*nexprtrees] = quadelems[i].coef;

         ++*nexprtrees;
      }

      /* add constant to first summand, if nonzero; need to divide by coef of this exprtree */
      if( nodedata->constant != 0.0 )
      {
         SCIP_EXPR* constexpr_;

         assert(*nexprtrees > 0);
         SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &constexpr_, SCIP_EXPR_CONST, nodedata->constant / exprtreecoefs[0]) );
         SCIP_CALL( SCIPexprtreeAddExpr(exprtrees[0], constexpr_, FALSE) );
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &varidx, exprgraph->nvars);

      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_POLYNOMIAL* nodedata;
      SCIP_EXPRDATA_MONOMIAL** monomials;
      SCIP_Real constant;
      int nmonomials;
      SCIP_EXPR* expr;
      int* childidxs;
      int j;

      nodedata = (SCIP_EXPRDATA_POLYNOMIAL*)node->data.data;
      monomials = nodedata->monomials;
      nmonomials = nodedata->nmonomials;
      constant = nodedata->constant;

      assert(exprtreessize >= nmonomials);
      assert(node->nchildren > 0);

      *nexprtrees = 0;

      /* buffer where to store mapping of expression graph variable indices to expression tree variable indices */
      SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &varidx, exprgraph->nvars) );

      for( i = 0; i < nmonomials; ++i )
      {
         /* initially, no variable appears in the expression tree */
         for( j = 0; j < exprgraph->nvars; ++j )
            varidx[j] = -1;
         nexprvars = 0;

         if( monomials[i]->nfactors == 1 )
         {
            /* create expression from the subgraph at only factor */
            SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, node->children[monomials[i]->childidxs[0]], &expr, &nexprvars, varidx) );

            /* put exponent in, if not 1.0 */
            if( monomials[i]->exponents[0] == 1.0 )
               ;
            else if( monomials[i]->exponents[0] == 2.0 )
            {
               SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &expr, SCIP_EXPR_SQUARE, expr) );
            }
            else if( EPSISINT(monomials[i]->exponents[0], 0.0) )  /*lint !e835*/
            {
               SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &expr, SCIP_EXPR_INTPOWER, expr, (int)monomials[i]->exponents[0]) );
            }
            else
            {
               SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &expr, SCIP_EXPR_REALPOWER, expr, monomials[i]->exponents[0]) );
            }
         }
         else if( monomials[i]->nfactors == 2 && monomials[i]->exponents[0] == 1.0 && monomials[i]->exponents[1] == 1.0 )
         {
            SCIP_EXPR* expr2;

            /* create expressions for both factors */
            SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, node->children[monomials[i]->childidxs[0]], &expr,  &nexprvars, varidx) );
            SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, node->children[monomials[i]->childidxs[1]], &expr2, &nexprvars, varidx) );

            /* create expression for product of factors */
            SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &expr, SCIP_EXPR_MUL, expr, expr2) );
         }
         else
         {
            SCIP_EXPRDATA_MONOMIAL* monomial;
            SCIP_EXPR** exprs;
            int f;

            /* create expression for each factor, assemble varidx and nexprvars
             * create child indices (= identity) */
            SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &exprs,     monomials[i]->nfactors) );
            SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &childidxs, monomials[i]->nfactors) );
            for( f = 0; f < monomials[i]->nfactors; ++f )
            {
               SCIP_CALL( exprgraphNodeCreateExpr(exprgraph, node->children[monomials[i]->childidxs[f]], &exprs[f], &nexprvars, varidx) );  /*lint !e644*/
               childidxs[f] = f;  /*lint !e644*/
            }

            /* create monomial and polynomial expression for this monomial
             * add also constant here, but need to divide by monomial coefficient, since we set the exprtreecoefs to monomial coef
             */
            SCIP_CALL( SCIPexprCreateMonomial(exprgraph->blkmem, &monomial, 1.0, monomials[i]->nfactors, childidxs, monomials[i]->exponents) );
            SCIP_CALL( SCIPexprCreatePolynomial(exprgraph->blkmem, &expr, monomials[i]->nfactors, exprs, 1, &monomial, constant / monomials[i]->coef, FALSE) );
            constant = 0.0;

            BMSfreeBlockMemoryArray(exprgraph->blkmem, &exprs,     monomials[i]->nfactors);
            BMSfreeBlockMemoryArray(exprgraph->blkmem, &childidxs, monomials[i]->nfactors);
         }

         /* create expression tree for expr */
         SCIP_CALL( SCIPexprtreeCreate(exprgraph->blkmem, &exprtrees[*nexprtrees], expr, nexprvars, 0, NULL) );

         /* copy variables into expression tree */
         if( nexprvars > 0 )
         {
            SCIP_ALLOC( BMSallocBlockMemoryArray(exprgraph->blkmem, &exprtrees[*nexprtrees]->vars, nexprvars) );
            for( j = 0; j < exprgraph->nvars; ++j )
            {
               assert(varidx[j] >= -1);
               assert(varidx[j] < nexprvars);
               if( varidx[j] >= 0 )
                  exprtrees[*nexprtrees]->vars[varidx[j]] = exprgraph->vars[j];
            }
         }

         exprtreecoefs[*nexprtrees] = monomials[i]->coef;

         ++*nexprtrees;
      }

      /* add constant to first summand, if still nonzero; need to divide by coefficient of the this exprtree */
      if( constant != 0.0 )
      {
         SCIP_EXPR* constexpr_;

         assert(*nexprtrees > 0);
         SCIP_CALL( SCIPexprCreate(exprgraph->blkmem, &constexpr_, SCIP_EXPR_CONST, constant / exprtreecoefs[0]) );
         SCIP_CALL( SCIPexprtreeAddExpr(exprtrees[0], constexpr_, FALSE) );
      }

      BMSfreeBlockMemoryArray(exprgraph->blkmem, &varidx, exprgraph->nvars);

      break;
   }

   default:
      SCIPerrorMessage("unexpected operator type %d\n", node->op);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/**@} */
