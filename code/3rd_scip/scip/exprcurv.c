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

/**@file   exprcurv.c
 * @ingroup OTHER_CFILES
 * @brief  functions to work with curvature (convex, concave, etc)
 * @author Stefan Vigerske
 *
 * Declarations are in pub_expr.h
 */

#include "scip/pub_expr.h"

/** curvature names as strings */
static
const char* curvnames[4] =
   {
      "unknown",
      "convex",
      "concave",
      "linear"
   };

#ifdef NDEBUG
#undef SCIPexprcurvAdd
#undef SCIPexprcurvNegate
#undef SCIPexprcurvMultiply
#endif

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

      /* something like x^(-2) may look convex on each side of zero, but is not convex on the whole interval
       * due to the singularity at 0.0 */
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

/** gives required curvature for base so that base^exponent has given curvature under given bounds on base and constant exponent
 *
 * returns curvature unknown if expected curvature cannot be obtained
 */
SCIP_EXPRCURV SCIPexprcurvPowerInv(
   SCIP_INTERVAL         basebounds,         /**< bounds on base function */
   SCIP_Real             exponent,           /**< exponent, must not be 0 */
   SCIP_EXPRCURV         powercurv           /**< expected curvature for power */
   )
{
   SCIP_Bool expisint;

   assert(basebounds.inf <= basebounds.sup);
   assert(exponent != 0.0);
   assert(powercurv != SCIP_EXPRCURV_UNKNOWN);

   if( exponent == 1.0 )
      return powercurv;

   /* power is usually never linear, now that exponent != 1 */
   if( powercurv == SCIP_EXPRCURV_LINEAR )
      return SCIP_EXPRCURV_UNKNOWN;

   expisint = EPSISINT(exponent, 0.0); /*lint !e835*/

   /* if exponent is fractional, then power is only defined for a non-negative base
    * boundtightening should have ensured this before calling this function,
    * but sometimes this does not work and so we correct this here for us
    */
   if( !expisint && basebounds.inf < 0.0 )
   {
      basebounds.inf = 0.0;
      if( basebounds.sup < 0.0 )
         return SCIP_EXPRCURV_UNKNOWN;
   }

   /* if basebounds contains 0.0, consider negative and positive interval separately, if possible */
   if( basebounds.inf < 0.0 && basebounds.sup > 0.0 )
   {
      SCIP_INTERVAL leftbounds;
      SCIP_INTERVAL rightbounds;
      SCIP_EXPRCURV leftcurv;
      SCIP_EXPRCURV rightcurv;

      /* something like x^(-2) may look convex on each side of zero, but is not convex on the whole
       * interval due to the singularity at 0.0 */
      if( exponent < 0.0 )
         return SCIP_EXPRCURV_UNKNOWN;

      SCIPintervalSetBounds(&leftbounds,  basebounds.inf, 0.0);
      SCIPintervalSetBounds(&rightbounds, 0.0, basebounds.sup);

      leftcurv = SCIPexprcurvPowerInv(leftbounds, exponent, powercurv);
      rightcurv = SCIPexprcurvPowerInv(rightbounds, exponent, powercurv);

      /* now need to intersect */
      if( leftcurv == SCIP_EXPRCURV_LINEAR )
         return rightcurv;
      if( rightcurv == SCIP_EXPRCURV_LINEAR )
         return leftcurv;
      if( leftcurv == SCIP_EXPRCURV_UNKNOWN || rightcurv == SCIP_EXPRCURV_UNKNOWN )
         return SCIP_EXPRCURV_UNKNOWN;
      assert(leftcurv == SCIP_EXPRCURV_CONVEX || leftcurv == SCIP_EXPRCURV_CONCAVE);
      assert(rightcurv == SCIP_EXPRCURV_CONVEX || rightcurv == SCIP_EXPRCURV_CONCAVE);
      return SCIP_EXPRCURV_LINEAR;
   }
   assert(basebounds.inf >= 0.0 || basebounds.sup <= 0.0);

   /* inverting the logic from SCIPexprcurvPower here */
   if( powercurv == SCIP_EXPRCURV_CONVEX )
   {
      SCIP_Real sign;

      if( basebounds.sup <= 0.0 && exponent < 0.0 && expisint && ((int)exponent)%2 == 0 )
         return SCIP_EXPRCURV_CONVEX;
      if( basebounds.inf >= 0.0 && exponent > 1.0 )
         return SCIP_EXPRCURV_CONVEX;
      if( basebounds.sup <= 0.0 && exponent > 1.0 && expisint && ((int)exponent)%2 == 0 )
         return SCIP_EXPRCURV_CONCAVE;
      if( basebounds.inf >= 0.0 && exponent < 0.0 )
         return SCIP_EXPRCURV_CONCAVE;

      /* base^(exponent-2) is negative, if base < 0.0 and exponent is odd */
      sign = exponent * (exponent - 1.0);
      assert(basebounds.inf >= 0.0 || expisint);
      if( basebounds.inf < 0.0 && ((int)exponent)%2 != 0 )
         sign *= -1.0;
      assert(sign != 0.0);

      if( sign > 0.0 )
         return SCIP_EXPRCURV_LINEAR;
   }
   else
   {
      SCIP_Real sign;

      assert(powercurv == SCIP_EXPRCURV_CONCAVE);  /* linear handled at top, unknown should not be the case */

      if( basebounds.sup <= 0.0 && exponent < 0.0 && expisint && ((int)exponent)%2 != 0 )
         return SCIP_EXPRCURV_CONVEX;
      if( basebounds.sup <= 0.0 && exponent > 1.0 && expisint && ((int)exponent)%2 != 0 )
         return SCIP_EXPRCURV_CONCAVE;
      if( basebounds.inf >= 0.0 && exponent < 1.0 && exponent >= 0.0 )
         return SCIP_EXPRCURV_CONCAVE;

      /* base^(exponent-2) is negative, if base < 0.0 and exponent is odd */
      sign = exponent * (exponent - 1.0);
      assert(basebounds.inf >= 0.0 || expisint);
      if( basebounds.inf < 0.0 && ((int)exponent)%2 != 0 )
         sign *= -1.0;
      assert(sign != 0.0);

      if( sign < 0.0 )
         return SCIP_EXPRCURV_LINEAR;
   }

   return SCIP_EXPRCURV_UNKNOWN;
}

/** gives curvature for a monomial with given curvatures and bounds for each factor
 *
 *  See Maranas and Floudas, Finding All Solutions of Nonlinearly Constrained Systems of Equations, JOGO 7, 1995
 *  for the categorization in the case that all factors are linear.
 *
 *  Exponents can also be negative or rational.
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
   SCIP_INTERVAL bounds;
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

      e = exponents != NULL ? exponents[j] : 1.0;
      bounds = factorbounds[f];  /*lint !e613*/

      /* if argument is negative, then exponent should be integer; correct bounds if that doesn't hold */
      if( !EPSISINT(e, 0.0) && bounds.inf < 0.0 )  /*lint !e835*/
      {
         bounds.inf = 0.0;
         if( bounds.sup < 0.0 )
            return SCIP_EXPRCURV_UNKNOWN;
      }

      if( bounds.inf < 0.0 && bounds.sup > 0.0 )
         return SCIP_EXPRCURV_UNKNOWN;

      if( e < 0.0 )
         ++nnegative;
      else
         ++npositive;
      sum += e;

      if( bounds.inf < 0.0 )
      {
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

/** for a monomial with given bounds for each factor, gives condition on the curvature of each factor,
 * so that monomial has a requested curvature, if possible
 *
 * @return whether `monomialcurv` can be achieved
 */
SCIP_Bool SCIPexprcurvMonomialInv(
   SCIP_EXPRCURV         monomialcurv,       /**< desired curvature */
   int                   nfactors,           /**< number of factors in monomial */
   SCIP_Real*            exponents,          /**< exponents in monomial, or NULL if all 1.0 */
   SCIP_INTERVAL*        factorbounds,       /**< bounds of each factor */
   SCIP_EXPRCURV*        factorcurv          /**< buffer to store required curvature of each factor */
   )
{
   int nnegative;
   int npositive;
   SCIP_INTERVAL bounds;
   SCIP_Real e;
   SCIP_Real sum;
   int j;

   assert(monomialcurv != SCIP_EXPRCURV_UNKNOWN);
   assert(nfactors >= 1);
   assert(factorbounds != NULL);
   assert(factorcurv != NULL);

   if( nfactors == 1 )
   {
      factorcurv[0] = SCIPexprcurvPowerInv(factorbounds[0], exponents != NULL ? exponents[0] : 1.0, monomialcurv);
      return factorcurv[0] != SCIP_EXPRCURV_UNKNOWN;
   }

   /* any decent monomial with at least 2 factors is not linear */
   if( monomialcurv == SCIP_EXPRCURV_LINEAR )
      return FALSE;

   /* count positive and negative exponents, sum of exponents; flip negative factors */
   nnegative = 0; /* number of negative exponents */
   npositive = 0; /* number of positive exponents */
   sum = 0.0;     /* sum of exponents */
   for( j = 0; j < nfactors; ++j )
   {
      e = exponents != NULL ? exponents[j] : 1.0;
      assert(e != 0.0);  /* should have been simplified away */

      bounds = factorbounds[j];

      /* if argument is negative, then exponent should be integer
       * if that didn't happen, consider argument as if non-negative
       */
      if( !EPSISINT(e, 0.0) && bounds.inf < 0.0 )  /*lint !e835*/
      {
         bounds.inf = 0.0;
         if( bounds.sup < 0.0 )
            return FALSE;
      }

      /* mixed signs are bad */
      if( bounds.inf < 0.0 && bounds.sup > 0.0 )
         return FALSE;

      if( e < 0.0 )
         ++nnegative;
      else
         ++npositive;
      sum += e;

      if( bounds.inf < 0.0 )
      {
         /* flip j'th argument: (f_j)^(exp_j) = (-1)^(exp_j) (-f_j)^(exp_j)
          * thus, negate monomial, if exponent is odd, i.e., (-1)^(exp_j) = -1
          */
         if( (int)e % 2 != 0 )
            monomialcurv = SCIPexprcurvNegate(monomialcurv);
      }
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

   if( monomialcurv == SCIP_EXPRCURV_CONVEX )
   {
      if( nnegative < nfactors-1 )  /* at least two positive exponents */
         return FALSE;
      if( nnegative < nfactors && !EPSGE(sum, 1.0, 1e-9) )  /* one negative exponent, but sum is not >= 1 */
         return FALSE;

      /* monomial will be convex, if each factor is convex if exp_j >= 0, or concave if exp_j <= 0, i.e., exp_j*f_j'' >= 0 */
      for( j = 0; j < nfactors; ++j )
      {
         e = exponents != NULL ? exponents[j] : 1.0;

         /* if factor is negative, then factorcurv[j] need to be flipped, which we can also get by flipping e */
         if( factorbounds[j].inf < 0.0 && EPSISINT(e, 0.0) )  /*lint !e835*/
            e = -e;
         if( e >= 0.0 )
            factorcurv[j] = SCIP_EXPRCURV_CONVEX;
         else
            factorcurv[j] = SCIP_EXPRCURV_CONCAVE;
      }
   }
   else
   {
      assert(monomialcurv == SCIP_EXPRCURV_CONCAVE);
      if( npositive < nfactors )  /* at least one negative exponent */
         return FALSE;
      if( !EPSLE(sum, 1.0, 1e-9) )  /* sum is not <= 1 */
         return FALSE;

      /* monomial will be concave, if each factor is concave */
      for( j = 0; j < nfactors; ++j )
      {
         e = exponents != NULL ? exponents[j] : 1.0;

         /* if factor is negative, then factorcurv[j] need to be flipped, i.e. convex */
         if( factorbounds[j].inf < 0.0 && EPSISINT(e, 0.0) )  /*lint !e835*/
            factorcurv[j] = SCIP_EXPRCURV_CONVEX;
         else
            factorcurv[j] = SCIP_EXPRCURV_CONCAVE;
      }
   }

   return TRUE;
}

/** gives name as string for a curvature */
const char* SCIPexprcurvGetName(
   SCIP_EXPRCURV         curv                /**< curvature */
   )
{
   assert(0 <= curv && curv <= SCIP_EXPRCURV_LINEAR);  /*lint !e685 !e2650 !e587 !e831 !e641 !e568*/

   return curvnames[curv];
}
