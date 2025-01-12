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

/**@file   expr_pow.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  power expression handler
 * @author Benjamin Mueller
 * @author Ksenia Bestuzheva
 *
 * @todo signpower for exponent < 1 ?
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint --e{835}*/
/*lint -e777*/

#include <string.h>

#include "scip/expr_pow.h"
#include "scip/pub_expr.h"
#include "scip/expr_value.h"
#include "scip/expr_product.h"
#include "scip/expr_sum.h"
#include "scip/expr_exp.h"
#include "scip/expr_abs.h"

#define POWEXPRHDLR_NAME         "pow"
#define POWEXPRHDLR_DESC         "power expression"
#define POWEXPRHDLR_PRECEDENCE   55000
#define POWEXPRHDLR_HASHKEY      SCIPcalcFibHash(21163.0)

#define SIGNPOWEXPRHDLR_NAME       "signpower"
#define SIGNPOWEXPRHDLR_DESC       "signed power expression"
#define SIGNPOWEXPRHDLR_PRECEDENCE 56000
#define SIGNPOWEXPRHDLR_HASHKEY    SCIPcalcFibHash(21163.1)

#define INITLPMAXPOWVAL          1e+06 /**< maximal allowed absolute value of power expression at bound,
                                        *   used for adjusting bounds in the convex case in initestimates */

/*
 * Data structures
 */

/** sign of a value (-1 or +1)
 *
 * 0.0 has sign +1 here (shouldn't matter, though)
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

#define SIGNPOW_ROOTS_KNOWN 10   /**< up to which (integer) exponents precomputed roots have been stored */

/** The positive root of the polynomial (n-1) y^n + n y^(n-1) - 1 is needed in separation.
 *  Here we store these roots for small integer values of n.
 */
static
SCIP_Real signpow_roots[SIGNPOW_ROOTS_KNOWN+1] = {
   -1.0,                     /* no root for n=0 */
   -1.0,                     /* no root for n=1 */
   0.41421356237309504880,   /* root for n=2 (-1+sqrt(2)) */
   0.5,                      /* root for n=3 */
   0.56042566045031785945,   /* root for n=4 */
   0.60582958618826802099,   /* root for n=5 */
   0.64146546982884663257,   /* root for n=6 */
   0.67033204760309682774,   /* root for n=7 */
   0.69428385661425826738,   /* root for n=8 */
   0.71453772716733489700,   /* root for n=9 */
   0.73192937842370733350    /* root for n=10 */
};

/** expression handler data */
struct SCIP_ExprhdlrData
{
   SCIP_Real             minzerodistance;    /**< minimal distance from zero to enforce for child in bound tightening */
   SCIP_Bool             warnedonpole;       /**< whether we warned on enforcing a minimal distance from zero for child */
};

/** expression data */
struct SCIP_ExprData
{
   SCIP_Real             exponent;           /**< exponent */
   SCIP_Real             root;               /**< positive root of (n-1) y^n + n y^(n-1) - 1, or
                                                  SCIP_INVALID if not computed yet */
};

/*
 * Local methods
 */

/** computes positive root of the polynomial (n-1) y^n + n y^(n-1) - 1 for n > 1 */
static
SCIP_RETCODE computeSignpowerRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            root,               /**< buffer where to store computed root */
   SCIP_Real             exponent            /**< exponent n */
   )
{
   SCIP_Real polyval;
   SCIP_Real gradval;
   int iter;

   assert(scip != NULL);
   assert(exponent > 1.0);
   assert(root != NULL);

   /* lookup for popular integer exponent */
   if( SCIPisIntegral(scip, exponent) && exponent-0.5 < SIGNPOW_ROOTS_KNOWN )
   {
      *root = signpow_roots[(int)SCIPfloor(scip, exponent+0.5)];
      return SCIP_OKAY;
   }

   /* lookup for weymouth exponent */
   if( SCIPisEQ(scip, exponent, 1.852) )
   {
      *root = 0.39821689389382575186;
      return SCIP_OKAY;
   }

   /* search for a positive root of (n-1) y^n + n y^(n-1) - 1
    * use the closest precomputed root as starting value
    */
   if( exponent >= SIGNPOW_ROOTS_KNOWN )
      *root = signpow_roots[SIGNPOW_ROOTS_KNOWN];
   else if( exponent <= 2.0 )
      *root = signpow_roots[2];
   else
      *root = signpow_roots[(int)SCIPfloor(scip, exponent)];

   for( iter = 0; iter < 1000; ++iter )
   {
      polyval = (exponent - 1.0) * pow(*root, exponent) + exponent * pow(*root, exponent - 1.0) - 1.0;
      if( fabs(polyval) < 1e-12 && SCIPisZero(scip, polyval) )
         break;

      /* gradient of (n-1) y^n + n y^(n-1) - 1 is n(n-1)y^(n-1) + n(n-1)y^(n-2) */
      gradval = (exponent - 1.0) * exponent * (pow(*root, exponent - 1.0) + pow(*root, exponent - 2.0));
      if( SCIPisZero(scip, gradval) )
         break;

      /* update root by adding -polyval/gradval (Newton's method) */
      *root -= polyval / gradval;
      if( *root < 0.0 )
         *root = 0.0;
   }

   if( !SCIPisZero(scip, polyval) )
   {
      SCIPerrorMessage("failed to compute root for exponent %g\n", exponent);
      return SCIP_ERROR;
   }
   SCIPdebugMsg(scip, "root for %g is %.20g, certainty = %g\n", exponent, *root, polyval);
   /* @todo cache root value for other expressions (an exponent seldom comes alone)?? (they are actually really fast to compute...) */

   return SCIP_OKAY;
}

/** computes negative root of the polynomial (n-1) y^n - n y^(n-1) + 1 for n < -1 */
static
SCIP_RETCODE computeHyperbolaRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            root,               /**< buffer where to store computed root */
   SCIP_Real             exponent            /**< exponent n */
   )
{
   SCIP_Real polyval;
   SCIP_Real gradval;
   int iter;

   assert(scip != NULL);
   assert(exponent < -1.0);
   assert(root != NULL);

   *root = -2.0;  /* that's the solution for n=-2 */

   for( iter = 0; iter < 1000; ++iter )
   {
      polyval = (exponent - 1.0) * pow(*root, exponent) - exponent * pow(*root, exponent - 1.0) + 1.0;
      if( fabs(polyval) < 1e-12 && SCIPisZero(scip, polyval) )
         break;

      /* gradient of (n-1) y^n - n y^(n-1) + 1 is n(n-1)y^(n-1) - n(n-1)y^(n-2) */
      gradval = (exponent - 1.0) * exponent * (pow(*root, exponent - 1.0) - pow(*root, exponent - 2.0));
      if( SCIPisZero(scip, gradval) )
         break;

      /* update root by adding -polyval/gradval (Newton's method) */
      *root -= polyval / gradval;
      if( *root >= 0.0 )
         *root = -1;
   }

   if( !SCIPisZero(scip, polyval) )
   {
      SCIPerrorMessage("failed to compute root for exponent %g\n", exponent);
      return SCIP_ERROR;
   }
   SCIPdebugMsg(scip, "root for %g is %.20g, certainty = %g\n", exponent, *root, polyval);
   /* @todo cache root value for other expressions (an exponent seldom comes alone)?? (they are actually really fast to compute...) */

   return SCIP_OKAY;
}

/** creates expression data */
static
SCIP_RETCODE createData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRDATA**       exprdata,           /**< pointer where to store expression data */
   SCIP_Real             exponent            /**< exponent of the power expression */
   )
{
   assert(exprdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, exprdata) );

   (*exprdata)->exponent = exponent;
   (*exprdata)->root = SCIP_INVALID;

   return SCIP_OKAY;
}

/** computes a tangent at a reference point by linearization
 *
 * for a normal power, linearization in xref is xref^exponent + exponent * xref^(exponent-1) (x - xref)
 * = (1-exponent) * xref^exponent + exponent * xref^(exponent-1) * x
 *
 * for a signpower, linearization is the same if xref is positive
 * for xref negative it is -(-xref)^exponent + exponent * (-xref)^(exponent-1) (x-xref)
 * = (1-exponent) * (-xref)^(exponent-1) * xref + exponent * (-xref)^(exponent-1) * x
 */
static
void computeTangent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             signpower,          /**< are we signpower or normal power */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_Real             xref,               /**< reference point where to linearize */
   SCIP_Real*            constant,           /**< buffer to store constant term of secant */
   SCIP_Real*            slope,              /**< buffer to store slope of secant */
   SCIP_Bool*            success             /**< buffer to store whether secant could be computed */
   )
{
   SCIP_Real xrefpow;

   assert(scip != NULL);
   assert(constant != NULL);
   assert(slope != NULL);
   assert(success != NULL);
   assert(xref != 0.0 || exponent > 0.0);
   /* non-integral exponent -> reference point must be >= 0 or we do signpower */
   assert(EPSISINT(exponent, 0.0) || signpower || !SCIPisNegative(scip, xref));

   /* TODO power is not differentiable at 0.0 for exponent < 0
    * should we forbid here that xref > 0, do something smart here, or just return success=FALSE?
    */
   /* assert(exponent >= 1.0 || xref > 0.0); */

   if( !EPSISINT(exponent, 0.0) && !signpower && xref < 0.0 )
      xref = 0.0;

   xrefpow = pow(signpower ? REALABS(xref) : xref, exponent - 1.0);

   /* if huge xref and/or exponent too large, then pow may overflow */
   if( !SCIPisFinite(xrefpow) )
   {
      *success = FALSE;
      return;
   }

   *constant = (1.0 - exponent) * xrefpow * xref;
   *slope = exponent * xrefpow;
   *success = TRUE;
}

/** computes a secant between lower and upper bound
 *
 * secant is xlb^exponent + (xub^exponent - xlb^exponent) / (xub - xlb) * (x - xlb)
 * = xlb^exponent - slope * xlb + slope * x  with slope = (xub^exponent - xlb^exponent) / (xub - xlb)
 * same if signpower
 */
static
void computeSecant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             signpower,          /**< are we signpower or normal power */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_Real             xlb,                /**< lower bound on x */
   SCIP_Real             xub,                /**< upper bound on x */
   SCIP_Real*            constant,           /**< buffer to store constant term of secant */
   SCIP_Real*            slope,              /**< buffer to store slope of secant */
   SCIP_Bool*            success             /**< buffer to store whether secant could be computed */
   )
{
   assert(scip != NULL);
   assert(constant != NULL);
   assert(slope != NULL);
   assert(success != NULL);
   assert(xlb >= 0.0 || EPSISINT(exponent, 0.0) || signpower);
   assert(xub >= 0.0 || EPSISINT(exponent, 0.0) || signpower);
   assert(exponent != 1.0);

   *success = FALSE;

   /* infinite bounds will not work */
   if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
      return;

   /* first handle some special cases */
   if( xlb == xub )
   {
      /* usually taken care of in separatePointPow already, but we might be called with different bounds here,
       * e.g., when handling odd or signed power
       */
      *slope = 0.0;
      *constant = pow(xlb, exponent);
   }
   else if( EPSISINT(exponent / 2.0, 0.0) && !signpower && xub > 0.1 && SCIPisFeasEQ(scip, xlb, -xub) )
   {
      /* for normal power with even exponents with xlb ~ -xub the slope would be very close to 0
       * since xub^n - xlb^n is prone to cancellation here, we omit computing this secant (it's probably useless)
       * unless the bounds are close to 0 as well (xub <= 0.1 in the "if" above)
       * or we have exactly xlb=-xub, where we can return a clean 0.0 (though it's probably useless)
       */
      if( xlb == -xub )
      {
         *slope = 0.0;
         *constant = pow(xlb, exponent);
      }
      else
      {
         return;
      }
   }
   else if( xlb == 0.0 && exponent > 0.0 )
   {
      assert(xub >= 0.0);
      *slope = pow(xub, exponent-1.0);
      *constant = 0.0;
   }
   else if( xub == 0.0 && exponent > 0.0 )
   {
      /* normal pow: slope = - xlb^exponent / (-xlb) = xlb^(exponent-1)
       * signpower:  slope = (-xlb)^exponent / (-xlb) = (-xlb)^(exponent-1)
       */
      assert(xlb <= 0.0);  /* so signpower or exponent is integral */
      if( signpower )
         *slope = pow(-xlb, exponent-1.0);
      else
         *slope = pow(xlb, exponent-1.0);
      *constant = 0.0;
   }
   else if( SCIPisEQ(scip, xlb, xub) && (!signpower || xlb >= 0.0 || xub <= 0.0) )
   {
      /* Computing the slope as (xub^n - xlb^n)/(xub-xlb) can lead to cancellation.
       * To avoid this, we replace xub^n by a Taylor expansion of pow at xlb:
       *   xub^n = xlb^n + n xlb^(n-1) (xub-xlb) + 0.5 n*(n-1) xlb^(n-2) (xub-xlb)^2 + 1/6 n*(n-1)*(n-2) xi^(n-3) (xub-xlb)^3  for some xlb < xi < xub
       * Dropping the last term, the slope is  (with an error of O((xub-xlb)^2) = 1e-18)
       *   n*xlb^(n-1) + 0.5 n*(n-1) xlb^(n-2)*(xub-xlb)
       * = n*xlb^(n-1) (1 - 0.5*(n-1)) + 0.5 n*(n-1) xlb^(n-2)*xub
       * = 0.5*n*((3-n)*xlb^(n-1) + (n-1) xlb^(n-2)*xub)
       *
       * test n=2: 0.5*2*((3-2)*xlb + (2-1) 1*xub) = xlb + xub    ok
       *      n=3: 0.5*3*((3-3)*xlb + (3-1) xlb*xub) = 3*xlb*xub ~ xlb^2 + xlb*xub + xub^2  ok
       *
       * The constant is
       *   xlb^n - 0.5*n*((3-n) xlb^(n-1) + (n-1) xlb^(n-2)*xub) * xlb
       * = xlb^n - 0.5*n*(3-n) xlb^n - 0.5*n*(n-1) xlb^(n-1)*xub
       * = (1-0.5*n*(3-n)) xlb^n - 0.5 n*(n-1) xlb^(n-1) xub
       *
       * test n=2: (1-0.5*2*(3-2)) xlb^2 - 0.5 2*(2-1) xlb xub = -xlb*xub
       *           old formula: xlb^2 - (xlb+xub) * xlb = -xlb*xub  ok
       *
       * For signpower with xub <= 0, we can negate xlb and xub:
       *   slope: (sign(xub)|xub|^n - sign(xlb)*|xlb|^n) / (xub-xlb) = -((-xub)^n - (-xlb)^n) / (xub - xlb) = ((-xub)^n - (-xlb)^n) / (-xub - (-xlb))
       *   constant: sign(xlb)|xlb|^n + slope * (xub - xlb) = -((-xlb)^n - slope * (xub - xlb)) = -((-xlb)^n + slope * ((-xub) - (-xlb)))
       */
      SCIP_Real xlb_n;   /* xlb^n */
      SCIP_Real xlb_n1;  /* xlb^(n-1) */
      SCIP_Real xlb_n2;  /* xlb^(n-2) */

      if( signpower && xub <= 0.0 )
      {
         xlb *= -1.0;
         xub *= -1.0;
      }

      xlb_n = pow(xlb, exponent);
      xlb_n1 = pow(xlb, exponent - 1.0);
      xlb_n2 = pow(xlb, exponent - 2.0);

      *slope = 0.5*exponent * ((3.0-exponent) * xlb_n1 + (exponent-1.0) * xlb_n2 * xub);
      *constant = (1.0 - 0.5*exponent*(3.0-exponent)) * xlb_n - 0.5*exponent*(exponent-1.0) * xlb_n1 * xub;

      if( signpower && xub <= 0.0 )
         *constant *= -1.0;
   }
   else
   {
      SCIP_Real lbval;
      SCIP_Real ubval;

      if( signpower )
         lbval = SIGN(xlb) * pow(REALABS(xlb), exponent);
      else
         lbval = pow(xlb, exponent);
      if( !SCIPisFinite(lbval) )
         return;

      if( signpower )
         ubval = SIGN(xub) * pow(REALABS(xub), exponent);
      else
         ubval = pow(xub, exponent);
      if( !SCIPisFinite(ubval) )
         return;

      /* we still can have bad numerics when xlb^exponent and xub^exponent are very close, but xlb and xub are not
       * for now, only check that things did not cancel out completely
       */
      if( lbval == ubval )
         return;

      *slope = (ubval - lbval) / (xub - xlb);
      *constant = lbval - *slope * xlb;
   }

   /* check whether we had overflows */
   if( !SCIPisFinite(*slope) || !SCIPisFinite(*constant) )
      return;

   *success = TRUE;
}

/** Separation for parabola
 *
 * - even positive powers: x^2, x^4, x^6 with x arbitrary, or
 * - positive powers > 1: x^1.5, x^2.5 with x >= 0
 <pre>
  100 +--------------------------------------------------------------------+
      |*               +                 +                +               *|
   90 |**                                                     x**2 ********|
      |  *                                                              *  |
   80 |-+*                                                              *+-|
      |   **                                                          **   |
   70 |-+   *                                                        *   +-|
      |     **                                                      **     |
   60 |-+     *                                                    *     +-|
      |       **                                                  **       |
   50 |-+       *                                                *       +-|
      |          **                                            **          |
   40 |-+          *                                          *          +-|
      |             **                                      **             |
   30 |-+            **                                    **            +-|
      |                **                                **                |
   20 |-+                **                            **                +-|
      |                   ***                        ***                   |
   10 |-+                   ***                    ***                   +-|
      |                +       *****     +    *****       +                |
    0 +--------------------------------------------------------------------+
     -10              -5                 0                5                10
 </pre>
 */
static
void estimateParabola(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_Bool             overestimate,       /**< should the power be overestimated? */
   SCIP_Real             xlb,                /**< lower bound on x */
   SCIP_Real             xub,                /**< upper bound on x */
   SCIP_Real             xref,               /**< reference point (where to linearize) */
   SCIP_Real*            constant,           /**< buffer to store constant term of estimator */
   SCIP_Real*            slope,              /**< buffer to store slope of estimator */
   SCIP_Bool*            islocal,            /**< buffer to store whether estimator only locally valid, that is,
                                                  it depends on given bounds */
   SCIP_Bool*            success             /**< buffer to store whether estimator could be computed */
   )
{
   assert(scip != NULL);
   assert(constant != NULL);
   assert(slope != NULL);
   assert(islocal != NULL);
   assert(success != NULL);
   assert((exponent >= 0.0 && EPSISINT(exponent/2.0, 0.0)) || (exponent > 1.0 && xlb >= 0.0));

   if( !overestimate )
   {
      computeTangent(scip, FALSE, exponent, xref, constant, slope, success);
      *islocal = FALSE;
   }
   else
   {
      /* overestimation -> secant */
      computeSecant(scip, FALSE, exponent, xlb, xub, constant, slope, success);
      *islocal = TRUE;
   }
}


/** Separation for signpower
 *
 * - odd positive powers, x^3, x^5, x^7
 * - sign(x)|x|^n for n > 1
 * - lower bound on x is negative (otherwise one should use separation for parabola)
 <pre>
  100 +--------------------------------------------------------------------+
      |                +                 +                +              **|
      |                                                   x*abs(x) ******* |
      |                                                              **    |
      |                                                            **      |
   50 |-+                                                       ***      +-|
      |                                                      ***           |
      |                                                   ***              |
      |                                               *****                |
      |                                          *****                     |
    0 |-+                        ****************                        +-|
      |                     *****                                          |
      |                *****                                               |
      |              ***                                                   |
      |           ***                                                      |
  -50 |-+      ***                                                       +-|
      |      **                                                            |
      |    **                                                              |
      |  **                                                                |
      |**              +                 +                +                |
 -100 +--------------------------------------------------------------------+
     -10              -5                 0                5                10
 </pre>
 */
static
void estimateSignedpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_Real             root,               /**< positive root of the polynomial (n-1) y^n + n y^(n-1) - 1,
                                                  if xubglobal > 0 */
   SCIP_Bool             overestimate,       /**< should the power be overestimated? */
   SCIP_Real             xlb,                /**< lower bound on x, assumed to be non-positive */
   SCIP_Real             xub,                /**< upper bound on x */
   SCIP_Real             xref,               /**< reference point (where to linearize) */
   SCIP_Real             xlbglobal,          /**< global lower bound on x */
   SCIP_Real             xubglobal,          /**< global upper bound on x */
   SCIP_Real*            constant,           /**< buffer to store constant term of estimator */
   SCIP_Real*            slope,              /**< buffer to store slope of estimator */
   SCIP_Bool*            islocal,            /**< buffer to store whether estimator only locally valid, that is,
                                                  it depends on given bounds */
   SCIP_Bool*            branchcand,         /**< buffer to indicate whether estimator would improve by branching
                                                  on it */
   SCIP_Bool*            success             /**< buffer to store whether estimator could be computed */
   )
{
   assert(scip != NULL);
   assert(constant != NULL);
   assert(slope != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);  /* the default */
   assert(success != NULL);
   assert(exponent >= 1.0);
   assert(xlb < 0.0); /* otherwise estimateParabola should have been called */
   assert(xubglobal <= 0.0 || (root > 0.0 && root < 1.0));

   *success = FALSE;

   if( !SCIPisPositive(scip, xub) )
   {
      /* easy case */
      if( !overestimate )
      {
         /* underestimator is secant */
         computeSecant(scip, TRUE, exponent, xlb, xub, constant, slope, success);
         *islocal = TRUE;
      }
      else
      {
         /* overestimator is tangent */

         /* we must linearize left of 0 */
         if( xref > 0.0 )
            xref = 0.0;

         computeTangent(scip, TRUE, exponent, xref, constant, slope, success);

         /* if global upper bound is > 0, then the tangent is only valid locally if the reference point is right of
          * -root*xubglobal
          */
         *islocal = SCIPisPositive(scip, xubglobal) && xref > -root * xubglobal;

         /* tangent doesn't move after branching */
         *branchcand = FALSE;
      }
   }
   else
   {
      SCIP_Real c;

      if( !overestimate )
      {
         /* compute the special point which decides between secant and tangent */
         c = -xlb * root;

         if( xref < c )
         {
            /* underestimator is secant between xlb and c */
            computeSecant(scip, TRUE, exponent, xlb, c, constant, slope, success);
            *islocal = TRUE;
         }
         else
         {
            /* underestimator is tangent */
            computeTangent(scip, TRUE, exponent, xref, constant, slope, success);

            /* if reference point is left of -root*xlbglobal (c w.r.t. global bounds),
             * then tangent is not valid w.r.t. global bounds
             */
            *islocal = xref < -root * xlbglobal;

            /* tangent doesn't move after branching */
            *branchcand = FALSE;
         }
      }
      else
      {
         /* compute the special point which decides between secant and tangent */
         c = -xub * root;

         if( xref <= c )
         {
            /* overestimator is tangent */
            computeTangent(scip, TRUE, exponent, xref, constant, slope, success);

            /* if reference point is right of -root*xubglobal (c w.r.t. global bounds),
             * then tangent is not valid w.r.t. global bounds
             */
            *islocal = xref > -root * xubglobal;

            /* tangent doesn't move after branching */
            *branchcand = FALSE;
         }
         else
         {
            /* overestimator is secant */
            computeSecant(scip, TRUE, exponent, c, xub, constant, slope, success);
            *islocal = TRUE;
         }
      }
   }
}

/** Separation for positive hyperbola
 *
 * - x^-2, x^-4 with x arbitrary
 * - x^-0.5, x^-1, x^-1.5, x^-3, x^-5 with x >= 0
 <pre>
  5 +----------------------------------------------------------------------+
    |                 +               * +*               +                 |
    |                                 *  *                 x**(-2) ******* |
  4 |-+                               *  *                               +-|
    |                                 *  *                                 |
    |                                 *  *                                 |
    |                                 *  *                                 |
  3 |-+                               *   *                              +-|
    |                                *    *                                |
    |                                *    *                                |
  2 |-+                              *    *                              +-|
    |                                *    *                                |
    |                               *      *                               |
  1 |-+                             *      *                             +-|
    |                               *      *                               |
    |                             **        **                             |
    |                   **********            **********                   |
  0 |*******************                                *******************|
    |                                                                      |
    |                 +                 +                +                 |
 -1 +----------------------------------------------------------------------+
   -10               -5                 0                5                 10
 </pre>
 */
static
void estimateHyperbolaPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_Real             root,               /**< negative root of the polynomial (n-1) y^n - n y^(n-1) + 1,
                                                  if x has mixed sign (w.r.t. global bounds?) and underestimating */
   SCIP_Bool             overestimate,       /**< should the power be overestimated? */
   SCIP_Real             xlb,                /**< lower bound on x */
   SCIP_Real             xub,                /**< upper bound on x */
   SCIP_Real             xref,               /**< reference point (where to linearize) */
   SCIP_Real             xlbglobal,          /**< global lower bound on x */
   SCIP_Real             xubglobal,          /**< global upper bound on x */
   SCIP_Real*            constant,           /**< buffer to store constant term of estimator */
   SCIP_Real*            slope,              /**< buffer to store slope of estimator */
   SCIP_Bool*            islocal,            /**< buffer to store whether estimator only locally valid, that is,
                                                  it depends on given bounds */
   SCIP_Bool*            branchcand,         /**< buffer to indicate whether estimator would improve by branching
                                                  on it */
   SCIP_Bool*            success             /**< buffer to store whether estimator could be computed */
   )
{
   assert(scip != NULL);
   assert(constant != NULL);
   assert(slope != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);  /* the default */
   assert(success != NULL);
   assert(exponent < 0.0);
   assert(EPSISINT(exponent/2.0, 0.0) || xlb >= 0.0);

   *success = FALSE;

   if( !overestimate )
   {
      if( xlb >= 0.0 || xub <= 0.0 )
      {
         /* underestimate and fixed sign -> tangent */

         /* make sure xref has the same sign as xlb,xub */
         if( xref < 0.0 && xlb >= 0.0 )
            xref = xlb;
         else if( xref > 0.0 && xub <= 0.0 )
            xref = xub;

         if( SCIPisZero(scip, xref) )
         {
            /* estimator would need to have an (essentially) infinite scope
             * first try to make up a better refpoint
             */
            if( xub > 0.0 )
            {
               /* thus xlb >= 0.0; stay close to xlb (probably = 0) */
               if( !SCIPisInfinity(scip, xub) )
                  xref = 0.9 * xlb + 0.1 * xub;
               else
                  xref = 0.1;
            }
            else
            {
               /* xub <= 0.0; stay close to xub (probably = 0) */
               if( !SCIPisInfinity(scip, -xlb) )
                  xref = 0.1 * xlb + 0.9 * xub;
               else
                  xref = 0.1;
            }

            /* if still close to 0, then also bounds are close to 0, then just give up */
            if( SCIPisZero(scip, xref) )
               return;
         }

         computeTangent(scip, FALSE, exponent, xref, constant, slope, success);

         /* tangent will not change if branching on x (even if only locally valid, see checks below) */
         *branchcand = FALSE;

         if( EPSISINT(exponent/2.0, 0.0) )
         {
            /* for even exponents (as in the picture):
             * if x has fixed sign globally, then our tangent is also globally valid
             * however, if x has mixed sign, then it depends on the constellation between reference point and global
             * bounds, whether the tangent is globally valid (see also the longer discussion for the mixed-sign
             * underestimator below )
             */
            if( xref > 0.0 && xlbglobal < 0.0 )
            {
               assert(xubglobal > 0.0);  /* since xref > 0.0 */
               assert(root < 0.0); /* root needs to be given */
               /* if on right side, then tangent is only locally valid if xref is too much to the left */
               *islocal = xref < xlbglobal * root;
            }
            else if( xref < 0.0 && xubglobal > 0.0 )
            {
               assert(xlbglobal < 0.0);  /* since xref < 0.0 */
               assert(root < 0.0); /* root needs to be given */
               /* if on left side, then tangent is only locally valid if xref is too much to the right */
               *islocal = xref > xubglobal * root;
            }
            else
               *islocal = FALSE;
         }
         else
         {
            /* for odd exponents, the tangent is only locally valid if the sign of x is not fixed globally */
            *islocal = xlbglobal * xubglobal < 0.0;
         }
      }
      else
      {
         /* underestimate but mixed sign */
         if( SCIPisInfinity(scip, -xlb) )
         {
            if( SCIPisInfinity(scip, xub) )
            {
               /* underestimator is constant 0, but that is globally valid */
               *constant = 0.0;
               *slope = 0.0;
               *islocal = FALSE;
               *success = TRUE;
               return;
            }

            /* switch sign of x (mirror on ordinate) to make left bound finite and use its estimator */
            estimateHyperbolaPositive(scip, exponent, root, overestimate, -xub, -xlb, -xref, -xubglobal, -xlbglobal,
                                      constant, slope, islocal, branchcand, success);
            if( *success )
               *slope = -*slope;
         }
         else
         {
            /* The convex envelope of x^exponent for x in [xlb, infinity] is a line (secant) between xlb and some positive
             * coordinate xhat, and x^exponent for x > xhat.
             * Further, on [xlb,xub] with xub < xhat, the convex envelope is the secant between xlb and xub.
             *
             * To find xhat, consider the affine-linear function  l(x) = xlb^n + c * (x - xlb)   where n = exponent
             * we look for a value of x such that f(x) and l(x) coincide and such that l(x) will be tangent to f(x) on that
             * point, that is
             * xhat > 0 such that f(xhat) = l(xhat) and f'(xhat) = l'(xhat)
             * => xhat^n = xlb^n + c * (xhat - xlb)   and   n * xhat^(n-1) = c
             * => xhat^n = xlb^n + n * xhat^n - n * xhat^(n-1) * xlb
             * => 0 = xlb^n + (n-1) * xhat^n - n * xhat^(n-1) * xlb
             *
             * Divide by xlb^n, one gets a polynomial that looks very much like the one for signpower, but a sign is
             * different (since this is *not signed* power):
             * 0 = 1 + (n-1) * y^n - n * y^(n-1)  where y = xhat/xlb
             *
             * The solution y < 0 (because xlb < 0 and we want xhat > 0) is what we expect to be given as "root".
             */
            assert(root < 0.0); /* root needs to be given */
            if( xref <= xlb * root )
            {
               /* If the reference point is left of xhat (=xlb*root), then we can take the
                * secant between xlb and root*xlb (= tangent at root*xlb).
                * However, if xub < root*xlb, then we can tilt the estimator to be the secant between xlb and xub.
                */
               computeSecant(scip, FALSE, exponent, xlb, MIN(xlb * root, xub), constant, slope, success);
               *islocal = TRUE;
            }
            else
            {
               /* If reference point is right of xhat, then take the tangent at xref.
                * This will still be underestimating for x in [xlb,0], too.
                * The tangent is globally valid, if we had also generated w.r.t. global bounds.
                */
               computeTangent(scip, FALSE, exponent, xref, constant, slope, success);
               *islocal = xref < xlbglobal * root;
               *branchcand = FALSE;
            }
         }
      }
   }
   else
   {
      /* overestimate and mixed sign -> pole is within domain -> cannot overestimate */
      if( xlb < 0.0 && xub > 0.0 )
         return;

      /* overestimate and fixed sign -> secant */
      computeSecant(scip, FALSE, exponent, xlb, xub, constant, slope, success);
      *islocal = TRUE;
   }

}

/** Separation for mixed-sign hyperbola
 *
 * - x^-1, x^-3, x^-5 without x >= 0 (either x arbitrary or x negative)
 <pre>
    +----------------------------------------------------------------------+
    |                 +                 *                +                 |
  4 |-+                                  *                 x**(-1) *******-|
    |                                    *                                 |
    |                                    *                                 |
    |                                    *                                 |
  2 |-+                                  *                               +-|
    |                                     *                                |
    |                                      **                              |
    |                                        *********                     |
  0 |*********************                            *********************|
    |                     *********                                        |
    |                              **                                      |
    |                                *                                     |
 -2 |-+                               *                                  +-|
    |                                 *                                    |
    |                                 *                                    |
    |                                 *                                    |
 -4 |-+                               *                                  +-|
    |                 +                *+                +                 |
    +----------------------------------------------------------------------+
   -10               -5                 0                5                 10
 </pre>
 */
static
void estimateHyperbolaMixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_Bool             overestimate,       /**< should the power be overestimated? */
   SCIP_Real             xlb,                /**< lower bound on x */
   SCIP_Real             xub,                /**< upper bound on x */
   SCIP_Real             xref,               /**< reference point (where to linearize) */
   SCIP_Real             xlbglobal,          /**< global lower bound on x */
   SCIP_Real             xubglobal,          /**< global upper bound on x */
   SCIP_Real*            constant,           /**< buffer to store constant term of estimator */
   SCIP_Real*            slope,              /**< buffer to store slope of estimator */
   SCIP_Bool*            islocal,            /**< buffer to store whether estimator only locally valid, that is,
                                                  it depends on given bounds */
   SCIP_Bool*            branchcand,         /**< buffer to indicate whether estimator would improve by branching
                                                  on it */
   SCIP_Bool*            success             /**< buffer to store whether estimator could be computed */
   )
{
   assert(scip != NULL);
   assert(constant != NULL);
   assert(slope != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);  /* the default */
   assert(success != NULL);
   assert(exponent < 0.0);
   assert(EPSISINT((exponent-1.0)/2.0, 0.0));
   assert(xlb < 0.0);

   *success = FALSE;

   if( xub <= 0.0 )
   {
      /* x is negative */
      if( !overestimate )
      {
         /* underestimation -> secant */
         computeSecant(scip, FALSE, exponent, xlb, xub, constant, slope, success);
         *islocal = TRUE;
      }
      else if( !SCIPisZero(scip, xlb/10.0) )
      {
         /* overestimation -> tangent */

         /* need to linearize left of 0 */
         if( xref > 0.0 )
            xref = 0.0;

         if( SCIPisZero(scip, xref) )
         {
            /* if xref is very close to 0.0, then slope would be infinite
             * try to move closer to lower bound (if xlb < -10*eps)
             */
            if( !SCIPisInfinity(scip, -xlb) )
               xref = 0.1*xlb + 0.9*xub;
            else
               xref = 0.1;
         }

         computeTangent(scip, FALSE, exponent, xref, constant, slope, success);

         /* if x does not have a fixed sign globally, then our tangent is not globally valid
          * (power is not convex on global domain)
          */
         *islocal = xlbglobal * xubglobal < 0.0;

         /* tangent doesn't move by branching */
         *branchcand = FALSE;
      }
      /* else: xlb is very close to zero, xub is <= 0, so slope would be infinite
       * (for any reference point in [xlb, xub]) -> do not estimate
       */
   }
   /* else: x has mixed sign -> pole is within domain -> cannot estimate */
}

/** Separation for roots with exponent in [0,1]
 *
 * - x^0.5 with x >= 0
 <pre>
  8 +----------------------------------------------------------------------+
    |             +             +              +             +             |
  7 |-+                                                     x**0.5 ********|
    |                                                             *********|
    |                                                      ********        |
  6 |-+                                             ********             +-|
    |                                         ******                       |
  5 |-+                                 ******                           +-|
    |                             ******                                   |
    |                        *****                                         |
  4 |-+                  ****                                            +-|
    |               *****                                                  |
  3 |-+          ****                                                    +-|
    |         ***                                                          |
    |      ***                                                             |
  2 |-+  **                                                              +-|
    |  **                                                                  |
  1 |**                                                                  +-|
    |*                                                                     |
    |*            +             +              +             +             |
  0 +----------------------------------------------------------------------+
    0             10            20             30            40            50
 </pre>
 */
static
void estimateRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_Bool             overestimate,       /**< should the power be overestimated? */
   SCIP_Real             xlb,                /**< lower bound on x */
   SCIP_Real             xub,                /**< upper bound on x */
   SCIP_Real             xref,               /**< reference point (where to linearize) */
   SCIP_Real*            constant,           /**< buffer to store constant term of estimator */
   SCIP_Real*            slope,              /**< buffer to store slope of estimator */
   SCIP_Bool*            islocal,            /**< buffer to store whether estimator only locally valid, that is,
                                                  it depends on given bounds */
   SCIP_Bool*            success             /**< buffer to store whether estimator could be computed */
   )
{
   assert(scip != NULL);
   assert(constant != NULL);
   assert(slope != NULL);
   assert(islocal != NULL);
   assert(success != NULL);
   assert(exponent > 0.0);
   assert(exponent < 1.0);
   assert(xlb >= 0.0);

   if( !overestimate )
   {
      /* underestimate -> secant */
      computeSecant(scip, FALSE, exponent, xlb, xub, constant, slope, success);
      *islocal = TRUE;
   }
   else
   {
      /* overestimate -> tangent */

      /* need to linearize right of 0 */
      if( xref < 0.0 )
         xref = 0.0;

      if( SCIPisZero(scip, xref) )
      {
         if( SCIPisZero(scip, xub) )
         {
            *success = FALSE;
            *islocal = FALSE;
            return;
         }

         /* if xref is 0 (then xlb=0 probably), then slope is infinite, then try to move away from 0 */
         if( xub < 0.2 )
            xref = 0.5 * xlb + 0.5 * xub;
         else
            xref = 0.1;
      }

      computeTangent(scip, FALSE, exponent, xref, constant, slope, success);
      *islocal = FALSE;
   }
}

/** builds an estimator for a power function */
static
SCIP_RETCODE buildPowEstimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRDATA*        exprdata,           /**< expression data */
   SCIP_Bool             overestimate,       /**< is this an overestimator? */
   SCIP_Real             childlb,            /**< local lower bound on the child */
   SCIP_Real             childub,            /**< local upper bound on the child */
   SCIP_Real             childglb,           /**< global lower bound on the child */
   SCIP_Real             childgub,           /**< global upper bound on the child */
   SCIP_Bool             childintegral,      /**< whether child is integral */
   SCIP_Real             refpoint,           /**< reference point */
   SCIP_Real             exponent,           /**< esponent */
   SCIP_Real*            coef,               /**< pointer to store the coefficient of the estimator */
   SCIP_Real*            constant,           /**< pointer to store the constant of the estimator */
   SCIP_Bool*            success,            /**< pointer to store whether the estimator was built successfully */
   SCIP_Bool*            islocal,            /**< pointer to store whether the estimator is valid w.r.t. local bounds
                                                  only */
   SCIP_Bool*            branchcand          /**< pointer to indicate whether to consider child for branching
                                                  (initialized to TRUE) */
   )
{
   SCIP_Bool isinteger;
   SCIP_Bool iseven;

   assert(scip != NULL);
   assert(exprdata != NULL);
   assert(coef != NULL);
   assert(constant != NULL);
   assert(success != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);

   isinteger = EPSISINT(exponent, 0.0);
   iseven = isinteger && EPSISINT(exponent / 2.0, 0.0);

   if( exponent == 2.0 )
   {
      /* important special case: quadratic case */
      /* initialize, because SCIPaddSquareXyz only adds to existing values */
      *success = TRUE;
      *coef = 0.0;
      *constant = 0.0;

      if( overestimate )
      {
         SCIPaddSquareSecant(scip, 1.0, childlb, childub, coef, constant, success);
         *islocal = TRUE; /* secants are only valid locally */
      }
      else
      {
         SCIPaddSquareLinearization(scip, 1.0, refpoint, childintegral, coef, constant, success);
         *islocal = FALSE; /* linearizations are globally valid */
         *branchcand = FALSE;  /* there is no improvement due to branching */
      }
   }
   else if( exponent > 0.0 && iseven )
   {
      estimateParabola(scip, exponent, overestimate, childlb, childub, refpoint, constant, coef, islocal, success);
      /* if estimate is locally valid, then we computed a secant and so branching can improve it */
      *branchcand = *islocal;
   }
   else if( exponent > 1.0 && childlb >= 0.0 )
   {
      /* make sure we linearize in convex region (if we will linearize) */
      if( refpoint < 0.0 )
         refpoint = 0.0;

      estimateParabola(scip, exponent, overestimate, childlb, childub, refpoint, constant, coef, islocal, success);

      /* if estimate is locally valid, then we computed a secant and so branching can improve it */
      *branchcand = *islocal;

      /* if odd power, then check whether tangent on parabola is also globally valid, that is reference point is
       * right of -root*global-lower-bound
       */
      if( !*islocal && !iseven && childglb < 0.0 )
      {
         if( SCIPisInfinity(scip, -childglb) )
            *islocal = TRUE;
         else
         {
            if( exprdata->root == SCIP_INVALID )
            {
               SCIP_CALL( computeSignpowerRoot(scip, &exprdata->root, exponent) );
            }
            *islocal = refpoint < exprdata->root * (-childglb);
         }
      }
   }
   else if( exponent > 1.0 )  /* and !iseven && childlb < 0.0 due to previous if */
   {
      /* compute root if not known yet; only needed if mixed sign (global child ub > 0) */
      if( exprdata->root == SCIP_INVALID && childgub > 0.0 )
      {
         SCIP_CALL( computeSignpowerRoot(scip, &exprdata->root, exponent) );
      }
      estimateSignedpower(scip, exponent, exprdata->root, overestimate, childlb, childub, refpoint,
                          -childglb, childgub, constant, coef, islocal, branchcand, success);
   }
   else if( exponent < 0.0 && (iseven || childlb >= 0.0) )
   {
      /* compute root if not known yet; only needed if mixed sign (globally) and iseven */
      if( exprdata->root == SCIP_INVALID && iseven )
      {
         SCIP_CALL( computeHyperbolaRoot(scip, &exprdata->root, exponent) );
      }
      estimateHyperbolaPositive(scip, exponent, exprdata->root, overestimate, childlb, childub, refpoint,
                                childglb, childgub, constant, coef, islocal, branchcand, success);
   }
   else if( exponent < 0.0 )
   {
      assert(!iseven); /* should hold due to previous if */
      assert(childlb < 0.0); /* should hold due to previous if */
      assert(isinteger); /* should hold because childlb < 0.0 (same as assert above) */

      estimateHyperbolaMixed(scip, exponent, overestimate, childlb, childub, refpoint, childglb, childgub,
                             constant, coef, islocal, branchcand, success);
   }
   else
   {
      assert(exponent < 1.0); /* the only case that should be left */
      assert(exponent > 0.0); /* should hold due to previous if */

      estimateRoot(scip, exponent, overestimate, childlb, childub, refpoint, constant, coef, islocal, success);

      /* if estimate is locally valid, then we computed a secant, and so branching can improve it */
      *branchcand = *islocal;
   }

   return SCIP_OKAY;
}

/** fills an array of reference points for estimating on the convex side */
static
void addTangentRefpoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             exponent,           /**< exponent of the power expression */
   SCIP_Real             lb,                 /**< lower bound on the child variable */
   SCIP_Real             ub,                 /**< upper bound on the child variable */
   SCIP_Real*            refpoints           /**< array to store the reference points */
   )
{
   SCIP_Real maxabsbnd;

   assert(refpoints != NULL);

   maxabsbnd = pow(INITLPMAXPOWVAL, 1 / exponent);

   /* make sure the absolute values of bounds are not too large */
   if( ub > -maxabsbnd )
      lb = MAX(lb, -maxabsbnd);
   if( lb < maxabsbnd )
      ub = MIN(ub,  maxabsbnd);

   /* in the case when ub < -maxabsbnd or lb > maxabsbnd, we still want to at least make bounds finite */
   if( SCIPisInfinity(scip, -lb) )
      lb = MIN(-10.0, ub - 0.1*REALABS(ub));  /*lint !e666 */
   if( SCIPisInfinity(scip,  ub) )
      ub = MAX( 10.0, lb + 0.1*REALABS(lb));  /*lint !e666 */

   refpoints[0] = (7.0 * lb + ub) / 8.0;
   refpoints[1] = (lb + ub) / 2.0;
   refpoints[2] = (lb + 7.0 * ub) / 8.0;
}

/** fills an array of reference points for sign(x)*abs(x)^n or x^n (n odd), where x has mixed signs
 *
 *  The reference points are: the lower and upper bounds (one for secant and one for tangent);
 *  and for the second tangent, the point on the convex part of the function between the point
 *  deciding between tangent and secant, and the corresponding bound
 */
static
SCIP_RETCODE addSignpowerRefpoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRDATA*        exprdata,           /**< expression data */
   SCIP_Real             lb,                 /**< lower bound on the child variable */
   SCIP_Real             ub,                 /**< upper bound on the child variable */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_Bool             underestimate,      /**< are the refpoints for an underestimator */
   SCIP_Real*            refpoints           /**< array to store the reference points */
   )
{
   assert(refpoints != NULL);

   if( (underestimate && SCIPisInfinity(scip, -lb)) || (!underestimate && SCIPisInfinity(scip, ub)) )
      return SCIP_OKAY;

   if( exprdata->root == SCIP_INVALID )
   {
      SCIP_CALL( computeSignpowerRoot(scip, &exprdata->root, exponent) );
   }

   /* make bounds finite (due to a previous if, only one can be infinite here) */
   if( SCIPisInfinity(scip, -lb) )
      lb = -ub * exprdata->root - 1.0;
   else if( SCIPisInfinity(scip,  ub) )
      ub = -lb * exprdata->root + 1.0;

   if( underestimate )
   {
      /* secant point */
      refpoints[0] = lb;

      /* tangent points, depending on the special point */
      if( -lb * exprdata->root < ub - 2.0 )
         refpoints[2] = ub;
      if( -lb * exprdata->root < ub - 4.0 )
         refpoints[1] = (-lb * exprdata->root + ub) / 2.0;
   }

   if( !underestimate )
   {
      /* secant point */
      refpoints[2] = ub;

      /* tangent points, depending on the special point */
      if( -ub * exprdata->root > lb + 2.0 )
         refpoints[0] = lb;
      if( -ub * exprdata->root > lb + 4.0 )
         refpoints[1] = (lb - ub * exprdata->root) / 2.0;
   }

   return SCIP_OKAY;
}

/** choose reference points for adding initestimates cuts for a power expression */
static
SCIP_RETCODE chooseRefpointsPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRDATA*        exprdata,           /**< expression data */
   SCIP_Real             lb,                 /**< lower bound on the child variable */
   SCIP_Real             ub,                 /**< upper bound on the child variable */
   SCIP_Real*            refpointsunder,     /**< array to store reference points for underestimators */
   SCIP_Real*            refpointsover,      /**< array to store reference points for overestimators */
   SCIP_Bool             underestimate,      /**< whether refpoints for underestimation are needed */
   SCIP_Bool             overestimate        /**< whether refpoints for overestimation are needed */
   )
{
   SCIP_Bool convex;
   SCIP_Bool concave;
   SCIP_Bool mixedsign;
   SCIP_Bool even;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(exprdata != NULL);
   assert(refpointsunder != NULL && refpointsover != NULL);

   exponent = exprdata->exponent;
   even = EPSISINT(exponent, 0.0) && EPSISINT(exponent / 2.0, 0.0);

   convex = FALSE;
   concave = FALSE;
   mixedsign = lb < 0.0 && ub > 0.0;

   /* convex case:
    * - parabola with an even degree or positive domain
    * - hyperbola with a positive domain
    * - even hyperbola with a negative domain
    */
   if( (exponent > 1.0 && (lb >= 0 || even)) || (exponent < 0.0 && lb >= 0) || (exponent < 0.0 && even && ub <= 0.0) )
      convex = TRUE;
   /* concave case:
    * - parabola or hyperbola with a negative domain and (due to previous if) an uneven degree
    * - root
    */
   else if( ub <= 0 || (exponent > 0.0 && exponent < 1.0) )
      concave = TRUE;

   if( underestimate )
   {
      if( convex )
         addTangentRefpoints(scip, exponent, lb, ub, refpointsunder);
      else if( (concave && !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub)) ||
               (exponent < 0.0 && even && mixedsign) ) /* concave with finite bounds or mixed even hyperbola */
      {
         /* for secant, refpoint doesn't matter, but we add it to signal that the corresponding cut should be created */
         refpointsunder[0] = (lb + ub) / 2.0;
      }
      else if( exponent > 1.0 && !even && mixedsign ) /* mixed signpower */
      {
         SCIP_CALL( addSignpowerRefpoints(scip, exprdata, lb, ub, exponent, TRUE, refpointsunder) );
      }
      else /* mixed odd hyperbola or an infinite bound */
         assert((exponent < 0.0 && !even && mixedsign) || SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub));
   }

   if( overestimate )
   {
      if( convex && !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
         refpointsover[0] = (lb + ub) / 2.0;
      else if( concave )
         addTangentRefpoints(scip, exponent, lb, ub, refpointsover);
      else if( exponent > 1.0 && !even && mixedsign ) /* mixed signpower */
      {
         SCIP_CALL( addSignpowerRefpoints(scip, exprdata, lb, ub, exponent, FALSE, refpointsover) );
      }
      else /* mixed hyperbola or an infinite bound */
         assert((exponent < 0.0 && mixedsign) || SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub));
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of expression handler
 */

/** compares two power expressions
 *
 * the order of two power (normal or signed) is base_1^expo_1 < base_2^expo_2 if and only if
 * base_1 < base2 or, base_1 = base_2 and expo_1 < expo_2
 */
static
SCIP_DECL_EXPRCOMPARE(comparePow)
{  /*lint --e{715}*/
   SCIP_Real expo1;
   SCIP_Real expo2;
   int compareresult;

   /**! [SnippetExprComparePow] */
   compareresult = SCIPcompareExpr(scip, SCIPexprGetChildren(expr1)[0], SCIPexprGetChildren(expr2)[0]);
   if( compareresult != 0 )
      return compareresult;

   expo1 = SCIPgetExponentExprPow(expr1);
   expo2 = SCIPgetExponentExprPow(expr2);

   return expo1 == expo2 ? 0 : expo1 < expo2 ? -1 : 1;
   /**! [SnippetExprComparePow] */
}

/** simplifies a pow expression
 *
 * Evaluates the power function when its child is a value expression
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyPow)
{  /*lint --e{715}*/
   SCIP_EXPR* base;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   base = SCIPexprGetChildren(expr)[0];
   assert(base != NULL);

   exponent = SCIPgetExponentExprPow(expr);

   SCIPdebugPrintf("[simplifyPow] simplifying power with expo %g\n", exponent);

   /* enforces POW1 */
   if( exponent == 0.0 )
   {
      SCIPdebugPrintf("[simplifyPow] POW1\n");
      /* TODO: more checks? */
      assert(!SCIPisExprValue(scip, base) || SCIPgetValueExprValue(base) != 0.0);
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, 1.0, ownercreate, ownercreatedata) );
      return SCIP_OKAY;
   }

   /* enforces POW2 */
   if( exponent == 1.0 )
   {
      SCIPdebugPrintf("[simplifyPow] POW2\n");
      *simplifiedexpr = base;
      SCIPcaptureExpr(*simplifiedexpr);
      return SCIP_OKAY;
   }

   /* enforces POW3 */
   if( SCIPisExprValue(scip, base) )
   {
      SCIP_Real baseval;

      SCIPdebugPrintf("[simplifyPow] POW3\n");
      baseval = SCIPgetValueExprValue(base);

      /* the assert below was failing on st_e35 for baseval=-1e-15 and fractional exponent
       * in the subNLP heuristic; I assume that this was because baseval was evaluated after
       * variable fixings and that there were just floating-point inaccuracies and 0 was meant,
       * so I treat -1e-15 as 0 here
       */
      if( baseval < 0.0 && fmod(exponent, 1.0) != 0.0 && baseval > -SCIPepsilon(scip) )
         baseval = 0.0;

      /* TODO check if those are all important asserts */
      assert(baseval >= 0.0 || fmod(exponent, 1.0) == 0.0);
      assert(baseval != 0.0 || exponent != 0.0);

      if( baseval != 0.0 || exponent > 0.0 )
      {
         SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, pow(baseval, exponent), ownercreate, ownercreatedata) );
         return SCIP_OKAY;
      }
   }

   /* enforces POW11 (exp(x)^n = exp(n*x)) */
   if( SCIPisExprExp(scip, base) )
   {
      SCIP_EXPR* child;
      SCIP_EXPR* prod;
      SCIP_EXPR* exponential;
      SCIP_EXPR* simplifiedprod;

      SCIPdebugPrintf("[simplifyPow] POW11\n");
      child = SCIPexprGetChildren(base)[0];

      /* multiply child of exponential with exponent */
      SCIP_CALL( SCIPcreateExprProduct(scip, &prod, 1, &child, exponent, ownercreate, ownercreatedata) );

      /* simplify product */
      SCIP_CALL( SCIPcallExprSimplify(scip, prod, &simplifiedprod, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &prod) );

      /* create exponential with new child */
      SCIP_CALL( SCIPcreateExprExp(scip, &exponential, simplifiedprod, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedprod) );

      /* the final simplified expression is the simplification of the just created exponential */
      SCIP_CALL( SCIPcallExprSimplify(scip, exponential, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exponential) );

      return SCIP_OKAY;
   }

   /* enforces POW10 */
   if( SCIPisExprVar(scip, base) )
   {
      SCIP_VAR* basevar;

      SCIPdebugPrintf("[simplifyPow] POW10\n");
      basevar = SCIPgetVarExprVar(base);

      assert(basevar != NULL);

      /* TODO: if exponent is negative, we could fix the binary variable to 1. However, this is a bit tricky because
       * variables can not be tighten in EXITPRE, where the simplify is also called
       */
      if( SCIPvarIsBinary(basevar) && exponent > 0.0 )
      {
         *simplifiedexpr = base;
         SCIPcaptureExpr(*simplifiedexpr);
         return SCIP_OKAY;
      }
   }

   if( EPSISINT(exponent, 0.0) )
   {
      SCIP_EXPR* aux;
      SCIP_EXPR* simplifiedaux;

      /* enforces POW5
       * given (pow n (prod 1.0 expr_1 ... expr_k)) we distribute the exponent:
       * -> (prod 1.0 (pow n expr_1) ... (pow n expr_k))
       * notes: - since base is simplified and its coefficient is 1.0 (SP8)
       *        - n is an integer (excluding 1 and 0; see POW1-2 above)
       */
      if( SCIPisExprProduct(scip, base) )
      {
         SCIP_EXPR* auxproduct;
         int i;

         /* create empty product */
         SCIP_CALL( SCIPcreateExprProduct(scip, &auxproduct, 0, NULL, 1.0, ownercreate, ownercreatedata) );

         for( i = 0; i < SCIPexprGetNChildren(base); ++i )
         {
            /* create (pow n expr_i) and simplify */
            SCIP_CALL( SCIPcreateExprPow(scip, &aux, SCIPexprGetChildren(base)[i], exponent, ownercreate, ownercreatedata) );
            SCIP_CALL( simplifyPow(scip, aux, &simplifiedaux, ownercreate, ownercreatedata) );
            SCIP_CALL( SCIPreleaseExpr(scip, &aux) );

            /* append (pow n expr_i) to product */
            SCIP_CALL( SCIPappendExprChild(scip, auxproduct, simplifiedaux) );
            SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedaux) );
         }

         /* simplify (prod 1.0 (pow n expr_1) ... (pow n expr_k))
          * this calls simplifyProduct directly, since we know its children are simplified */
         SCIP_CALL( SCIPcallExprSimplify(scip, auxproduct, simplifiedexpr, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &auxproduct) );
         return SCIP_OKAY;
      }

      /* enforces POW6
       * given (pow n (sum 0.0 coef expr)) we can move `pow` inside `sum`:
       * (pow n (sum 0.0 coef expr) ) -> (sum 0.0 coef^n (pow n expr))
       * notes: - since base is simplified and its constant is 0, then coef != 1.0 (SS7)
       *        - n is an integer (excluding 1 and 0; see POW1-2 above)
       */
      if( SCIPisExprSum(scip, base) && SCIPexprGetNChildren(base) == 1 && SCIPgetConstantExprSum(base) == 0.0 )
      {
         SCIP_Real newcoef;

         SCIPdebugPrintf("[simplifyPow] seeing a sum with one term, exponent %g\n", exponent);

         /* assert SS7 holds */
         assert(SCIPgetCoefsExprSum(base)[0] != 1.0);

         /* create (pow n expr) and simplify it
          * note: we call simplifyPow directly, since we know that `expr` is simplified */
         newcoef = pow(SCIPgetCoefsExprSum(base)[0], exponent);
         SCIP_CALL( SCIPcreateExprPow(scip, &aux, SCIPexprGetChildren(base)[0], exponent, ownercreate, ownercreatedata) );
         SCIP_CALL( simplifyPow(scip, aux, &simplifiedaux, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &aux) );

         /* create (sum (pow n expr)) and simplify it
          * this calls simplifySum directly, since we know its children are simplified */
         SCIP_CALL( SCIPcreateExprSum(scip, &aux, 1, &simplifiedaux, &newcoef, 0.0, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPcallExprSimplify(scip, aux, simplifiedexpr, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &aux) );
         SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedaux) );
         return SCIP_OKAY;
      }

      /* enforces POW7
       * (const + sum alpha_i expr_i)^2 = sum alpha_i^2 expr_i^2
       * + sum_{j < i} 2*alpha_i alpha_j expr_i expr_j
       * + sum const alpha_i expr_i
       * TODO: put some limits on the number of children of the sum being expanded
       */
      if( SCIPisExprSum(scip, base) && exponent == 2.0 )
      {
         int i;
         int nchildren;
         int nexpandedchildren;
         SCIP_EXPR* expansion;
         SCIP_EXPR** expandedchildren;
         SCIP_Real* coefs;
         SCIP_Real constant;

         SCIPdebugPrintf("[simplifyPow] expanding sum^%g\n", exponent);

         nchildren = SCIPexprGetNChildren(base);
         nexpandedchildren = nchildren * (nchildren + 1) / 2 + nchildren;
         SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nexpandedchildren) );
         SCIP_CALL( SCIPallocBufferArray(scip, &expandedchildren, nexpandedchildren) );

         for( i = 0; i < nchildren; ++i )
         {
            int j;
            SCIP_EXPR* expansionchild;
            SCIP_EXPR* prodchildren[2];
            prodchildren[0] = SCIPexprGetChildren(base)[i];

            /* create and simplify expr_i * expr_j */
            for( j = 0; j < i; ++j )
            {
               prodchildren[1] = SCIPexprGetChildren(base)[j];
               coefs[i*(i+1)/2 + j] = 2 * SCIPgetCoefsExprSum(base)[i] * SCIPgetCoefsExprSum(base)[j];

               SCIP_CALL( SCIPcreateExprProduct(scip, &expansionchild, 2, prodchildren, 1.0, ownercreate,
                          ownercreatedata) );
               SCIP_CALL( SCIPcallExprSimplify(scip, expansionchild, &expandedchildren[i*(i+1)/2 + j],
                          ownercreate, ownercreatedata) ); /* this calls simplifyProduct */
               SCIP_CALL( SCIPreleaseExpr(scip, &expansionchild) );
            }
            /* create and simplify expr_i * expr_i */
            prodchildren[1] = SCIPexprGetChildren(base)[i];
            coefs[i*(i+1)/2 + i] = SCIPgetCoefsExprSum(base)[i] * SCIPgetCoefsExprSum(base)[i];

            SCIP_CALL( SCIPcreateExprProduct(scip, &expansionchild, 2, prodchildren, 1.0, ownercreate,
                       ownercreatedata) );
            SCIP_CALL( SCIPcallExprSimplify(scip, expansionchild, &expandedchildren[i*(i+1)/2 + i], ownercreate,
                       ownercreatedata) ); /* this calls simplifyProduct */
            SCIP_CALL( SCIPreleaseExpr(scip, &expansionchild) );
         }
         /* create const * alpha_i expr_i */
         for( i = 0; i < nchildren; ++i )
         {
            coefs[i + nexpandedchildren - nchildren] = 2 * SCIPgetConstantExprSum(base) * SCIPgetCoefsExprSum(base)[i];
            expandedchildren[i + nexpandedchildren - nchildren] = SCIPexprGetChildren(base)[i];
         }

         constant = SCIPgetConstantExprSum(base);
         constant *= constant;
         /* create sum of all the above and simplify it with simplifySum since all of its children are simplified! */
         SCIP_CALL( SCIPcreateExprSum(scip, &expansion, nexpandedchildren, expandedchildren, coefs, constant,
                    ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPcallExprSimplify(scip, expansion, simplifiedexpr, ownercreate,
                    ownercreatedata) ); /* this calls simplifySum */

         /* release everything */
         SCIP_CALL( SCIPreleaseExpr(scip, &expansion) );
         /* release the *created* expanded children */
         for( i = 0; i < nexpandedchildren - nchildren; ++i )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &expandedchildren[i]) );
         }
         SCIPfreeBufferArray(scip, &expandedchildren);
         SCIPfreeBufferArray(scip, &coefs);

         return SCIP_OKAY;
      }
   }
   else
   {
      /* enforces POW9
       *
       * FIXME code of POW6 is very similar
       */
      if( SCIPexprGetNChildren(base) == 1
         && SCIPisExprSum(scip, base)
         && SCIPgetConstantExprSum(base) == 0.0
         && SCIPgetCoefsExprSum(base)[0] >= 0.0 )
      {
         SCIP_EXPR* simplifiedaux;
         SCIP_EXPR* aux;
         SCIP_Real newcoef;

         SCIPdebugPrintf("[simplifyPow] seeing a sum with one term, exponent %g\n", exponent);
         /* assert SS7 holds */
         assert(SCIPgetCoefsExprSum(base)[0] != 1.0);

         /* create (pow n expr) and simplify it
          * note: we call simplifyPow directly, since we know that `expr` is simplified */
         SCIP_CALL( SCIPcreateExprPow(scip, &aux, SCIPexprGetChildren(base)[0], exponent, ownercreate,
                    ownercreatedata) );
         SCIP_CALL( simplifyPow(scip, aux, &simplifiedaux, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &aux) );

         /* create (sum (pow n expr)) and simplify it
          * this calls simplifySum directly, since we know its child is simplified! */
         newcoef = pow(SCIPgetCoefsExprSum(base)[0], exponent);
         SCIP_CALL( SCIPcreateExprSum(scip, &aux, 1, &simplifiedaux, &newcoef, 0.0, ownercreate,
                    ownercreatedata) );
         SCIP_CALL( SCIPcallExprSimplify(scip, aux, simplifiedexpr, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &aux) );
         SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedaux) );

         return SCIP_OKAY;
      }
   }

   /* enforces POW8
    * given (pow n (pow expo expr)) we distribute the exponent:
    * -> (pow n*expo expr)
    * notes: n is not 1 or 0, see POW1-2 above
    */
   if( SCIPisExprPower(scip, base) )
   {
      SCIP_Real newexponent;
      SCIP_Real baseexponent;

      baseexponent = SCIPgetExponentExprPow(base);
      newexponent = baseexponent * exponent;

      /* some checks (see POW8 definition in scip_expr.h) to make sure we don't loose an
       * implicit SCIPexprGetChildren(base)[0] >= 0 constraint
       *
       * if newexponent is fractional, then we will still need expr >= 0
       * if both exponents were integer, then we never required and will not require expr >= 0
       * if base exponent was an even integer, then we did not require expr >= 0
       * (but may need to use |expr|^newexponent)
       */
      if( !EPSISINT(newexponent, 0.0) ||
          (EPSISINT(baseexponent, 0.0) && EPSISINT(exponent, 0.0)) ||
          (EPSISINT(baseexponent, 0.0) && ((int)baseexponent) % 2 == 0) )
      {
         SCIP_EXPR* aux;

         if( EPSISINT(baseexponent, 0.0) && ((int)baseexponent) % 2 == 0 &&
             (!EPSISINT(newexponent, 0.0) || ((int)newexponent) % 2 == 1) )
         {
            /* If base exponent was even integer and new exponent will be fractional,
             * then simplify to |expr|^newexponent to allow eval for expr < 0.
             * If base exponent was even integer and new exponent will be odd integer,
             * then simplify to |expr|^newexponent to preserve value for expr < 0.
             */
            SCIP_EXPR* simplifiedaux;

            SCIP_CALL( SCIPcreateExprAbs(scip, &aux, SCIPexprGetChildren(base)[0], ownercreate, ownercreatedata) );
            SCIP_CALL( SCIPcallExprSimplify(scip, aux, &simplifiedaux, ownercreate, ownercreatedata) );
            SCIP_CALL( SCIPreleaseExpr(scip, &aux) );
            SCIP_CALL( SCIPcreateExprPow(scip, &aux, simplifiedaux, newexponent, ownercreate, ownercreatedata) );
            SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedaux) );
         }
         else
         {
            SCIP_CALL( SCIPcreateExprPow(scip, &aux, SCIPexprGetChildren(base)[0], newexponent, ownercreate,
                       ownercreatedata) );
         }

         SCIP_CALL( simplifyPow(scip, aux, simplifiedexpr, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &aux) );

         return SCIP_OKAY;
      }
   }

   SCIPdebugPrintf("[simplifyPow] power is simplified\n");
   *simplifiedexpr = expr;

   /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
   SCIPcaptureExpr(*simplifiedexpr);

   return SCIP_OKAY;
}

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrPow)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrPow(scip) );

   return SCIP_OKAY;
}

/** expression handler free callback */
static
SCIP_DECL_EXPRFREEHDLR(freehdlrPow)
{  /*lint --e{715}*/
   assert(exprhdlrdata != NULL);
   assert(*exprhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, exprhdlrdata);

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataPow)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* sourceexprdata;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   sourceexprdata = SCIPexprGetData(sourceexpr);
   assert(sourceexprdata != NULL);

   *targetexprdata = NULL;

   SCIP_CALL( createData(targetscip, targetexprdata, sourceexprdata->exponent) );

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataPow)
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
/** @todo: use precedence for better printing */
static
SCIP_DECL_EXPRPRINT(printPow)
{  /*lint --e{715}*/
   assert(expr != NULL);

   /**! [SnippetExprPrintPow] */
   switch( stage )
   {
      case SCIP_EXPRITER_ENTEREXPR :
      {
         /* print function with opening parenthesis */
         SCIPinfoMessage(scip, file, "(");
         break;
      }

      case SCIP_EXPRITER_VISITINGCHILD :
      {
         assert(currentchild == 0);
         break;
      }

      case SCIP_EXPRITER_LEAVEEXPR :
      {
         SCIP_Real exponent = SCIPgetExponentExprPow(expr);

         /* print closing parenthesis */
         if( exponent >= 0.0 )
            SCIPinfoMessage(scip, file, ")^%g", exponent);
         else
            SCIPinfoMessage(scip, file, ")^(%g)", exponent);

         break;
      }

      case SCIP_EXPRITER_VISITEDCHILD :
      default:
         break;
   }
   /**! [SnippetExprPrintPow] */

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalPow)
{  /*lint --e{715}*/
   SCIP_Real exponent;
   SCIP_Real base;

   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID);

   exponent = SCIPgetExponentExprPow(expr);
   base = SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]);

   *val = pow(base, exponent);

   /* if there is a domain, pole, or range error, pow() should return some kind of NaN, infinity, or HUGE_VAL
    * we could also work with floating point exceptions or errno, but I am not sure this would be thread-safe
    */
   if( !SCIPisFinite(*val) || *val == HUGE_VAL || *val == -HUGE_VAL )
      *val = SCIP_INVALID;

   return SCIP_OKAY;
}

/** derivative evaluation callback
 *
 * computes <gradient, children.dot>
 * if expr is child^p, then computes
 * p child^(p-1) dot(child)
 */
static
SCIP_DECL_EXPRFWDIFF(fwdiffPow)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_Real exponent;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(!SCIPisExprValue(scip, child));
   assert(SCIPexprGetDot(child) != SCIP_INVALID);

   exponent = SCIPgetExponentExprPow(expr);
   assert(exponent != 0.0);

   /* x^exponent is not differentiable for x = 0 and exponent in ]0,1[ */
   if( exponent > 0.0 && exponent < 1.0 && SCIPexprGetEvalValue(child) == 0.0 )
      *dot = SCIP_INVALID;
   else
      *dot = exponent * pow(SCIPexprGetEvalValue(child), exponent - 1.0) * SCIPexprGetDot(child);

   return SCIP_OKAY;
}

/** expression backward forward derivative evaluation callback
 *
 * computes partial/partial child ( <gradient, children.dot> )
 * if expr is child^n, then computes
 * n * (n - 1) child^(n-2) dot(child)
 */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffPow)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_Real exponent;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID);
   assert(childidx == 0);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(!SCIPisExprValue(scip, child));
   assert(SCIPexprGetDot(child) != SCIP_INVALID);

   exponent = SCIPgetExponentExprPow(expr);
   assert(exponent != 0.0);

   /* x^exponent is not twice differentiable for x = 0 and exponent in ]0,1[ u ]1,2[ */
   if( exponent > 0.0 && exponent < 2.0 && SCIPexprGetEvalValue(child) == 0.0 && exponent != 1.0 )
      *bardot = SCIP_INVALID;
   else
      *bardot = exponent * (exponent - 1.0) * pow(SCIPexprGetEvalValue(child), exponent - 2.0) * SCIPexprGetDot(child);

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffPow)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_Real childval;
   SCIP_Real exponent;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);
   assert(childidx == 0);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(!SCIPisExprValue(scip, child));

   childval = SCIPexprGetEvalValue(child);
   assert(childval != SCIP_INVALID);

   exponent = SCIPgetExponentExprPow(expr);
   assert(exponent != 0.0);

   /* x^exponent is not differentiable for x = 0 and exponent in ]0,1[ */
   if( exponent > 0.0 && exponent < 1.0 && childval == 0.0 )
      *val = SCIP_INVALID;
   else
      *val = exponent * pow(childval, exponent - 1.0);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalPow)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;
   SCIP_Real exponent;

   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   exponent = SCIPgetExponentExprPow(expr);

   if( exponent < 0.0 )
   {
      SCIP_EXPRHDLRDATA* exprhdlrdata;
      exprhdlrdata = SCIPexprhdlrGetData(SCIPexprGetHdlr(expr));
      assert(exprhdlrdata != NULL);

      if( exprhdlrdata->minzerodistance > 0.0 )
      {
         /* avoid small interval around 0 if possible, see also reversepropPow */
         if( childinterval.inf > -exprhdlrdata->minzerodistance && childinterval.inf < exprhdlrdata->minzerodistance )
         {
            if( !exprhdlrdata->warnedonpole && SCIPgetVerbLevel(scip) > SCIP_VERBLEVEL_NONE )
            {
               SCIPinfoMessage(scip, NULL, "Changing lower bound for child of pow(.,%g) from %g to %g.\n"
                  "Check your model formulation or use option expr/" POWEXPRHDLR_NAME "/minzerodistance to avoid this warning.\n",
                  exponent, childinterval.inf, exprhdlrdata->minzerodistance);
               SCIPinfoMessage(scip, NULL, "Expression: ");
               SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
               SCIPinfoMessage(scip, NULL, "\n");
               exprhdlrdata->warnedonpole = TRUE;
            }
            childinterval.inf = exprhdlrdata->minzerodistance;
         }
         else if( childinterval.sup < exprhdlrdata->minzerodistance
                  && childinterval.sup > -exprhdlrdata->minzerodistance )
         {
            if( !exprhdlrdata->warnedonpole && SCIPgetVerbLevel(scip) > SCIP_VERBLEVEL_NONE )
            {
               SCIPinfoMessage(scip, NULL, "Changing upper bound for child of pow(.,%g) from %g to %g.\n"
                  "Check your model formulation or use option expr/" POWEXPRHDLR_NAME "/minzerodistance to avoid this warning.\n",
                  exponent, childinterval.sup, -exprhdlrdata->minzerodistance);
               SCIPinfoMessage(scip, NULL, "Expression: ");
               SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
               SCIPinfoMessage(scip, NULL, "\n");
               exprhdlrdata->warnedonpole = TRUE;
            }
            childinterval.sup = -exprhdlrdata->minzerodistance;
         }
      }
   }

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
   {
      SCIPintervalSetEmpty(interval);
      return SCIP_OKAY;
   }

   SCIPintervalPowerScalar(SCIP_INTERVAL_INFINITY, interval, childinterval, exponent);

   /* make sure 0^negative is an empty interval, as some other codes do not handle intervals like [inf,inf] well
    * TODO maybe change SCIPintervalPowerScalar?
    */
   if( exponent < 0.0 && childinterval.inf == 0.0 && childinterval.sup == 0.0 )
      SCIPintervalSetEmpty(interval);

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_EXPRESTIMATE(estimatePow)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   SCIP_EXPR* child;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real exponent;
   SCIP_Bool isinteger;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), POWEXPRHDLR_NAME) == 0);
   assert(refpoint != NULL);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);  /* the default */
   assert(success != NULL);

   *success = FALSE;

   /* get aux variables: we over- or underestimate childvar^exponent  */
   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   SCIPdebugMsg(scip, "%sestimation of x^%g at x=%.15g\n",
      overestimate ? "over" : "under", SCIPgetExponentExprPow(expr), *refpoint);

   /* we can not generate a cut at +/- infinity */
   if( SCIPisInfinity(scip, REALABS(*refpoint)) )
      return SCIP_OKAY;

   childlb = localbounds[0].inf;
   childub = localbounds[0].sup;

   exprdata = SCIPexprGetData(expr);
   exponent = exprdata->exponent;
   assert(exponent != 1.0 && exponent != 0.0); /* this should have been simplified */

   /* if child is constant, then return a constant estimator
    * this can help with small infeasibilities if boundtightening is relaxing bounds too much
    */
   if( childlb == childub )
   {
      *coefs = 0.0;
      *constant = pow(childlb, exponent);
      *success = TRUE;
      *islocal = globalbounds[0].inf != globalbounds[0].sup;
      *branchcand = FALSE;
      return SCIP_OKAY;
   }

   isinteger = EPSISINT(exponent, 0.0);

   /* if exponent is not integral, then child must be non-negative */
   if( !isinteger && childlb < 0.0 )
   {
      /* somewhere we should have tightened the bound on x, but small tightening are not always applied by SCIP
       * it is ok to do this tightening here, but let's assert that we were close to 0.0 already
       */
      assert(SCIPisFeasZero(scip, childlb));
      childlb = 0.0;
   }
   assert(isinteger || childlb >= 0.0);

   SCIP_CALL( buildPowEstimator(scip, exprdata, overestimate, childlb, childub, globalbounds[0].inf,
         globalbounds[0].sup, SCIPexprIsIntegral(child), MAX(childlb, *refpoint), exponent, coefs,
         constant, success, islocal, branchcand) );

   return SCIP_OKAY;
}

/** expression reverse propagaton callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropPow)
{  /*lint --e{715}*/
   SCIP_INTERVAL child;
   SCIP_INTERVAL interval;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   exponent = SCIPgetExponentExprPow(expr);
   child = childrenbounds[0];

   SCIPdebugMsg(scip, "reverseprop x^%g in [%.15g,%.15g], x = [%.15g,%.15g]", exponent, bounds.inf, bounds.sup,
         child.inf, child.sup);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, child) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, bounds) )
   {
      /* if exponent is not integral, then make sure that child is non-negative */
      if( !EPSISINT(exponent, 0.0) && child.inf < 0.0 )
      {
         SCIPintervalSetBounds(&interval, 0.0, child.sup);
      }
      else
      {
         SCIPdebugMsgPrint(scip, "-> no improvement\n");
         return SCIP_OKAY;
      }
   }
   else
   {
      /* f = pow(c0, alpha) -> c0 = pow(f, 1/alpha) */
      SCIPintervalPowerScalarInverse(SCIP_INTERVAL_INFINITY, &interval, child, exponent, bounds);
   }

   if( exponent < 0.0 )
   {
      SCIP_EXPRHDLRDATA* exprhdlrdata;

      exprhdlrdata = SCIPexprhdlrGetData(SCIPexprGetHdlr(expr));
      assert(exprhdlrdata != NULL);

      if( exprhdlrdata->minzerodistance > 0.0 )
      {
         /* push lower bound from >= -epsilon to >=  epsilon to avoid pole at 0 (domain error)
          * push upper bound from <=  epsilon to <= -epsilon to avoid pole at 0 (domain error)
          * this can lead to a cutoff if domain would otherwise be very close around 0
          */
         if( interval.inf > -exprhdlrdata->minzerodistance && interval.inf < exprhdlrdata->minzerodistance )
         {
            if( !exprhdlrdata->warnedonpole && SCIPgetVerbLevel(scip) > SCIP_VERBLEVEL_NONE )
            {
               SCIPinfoMessage(scip, NULL, "Changing lower bound for child of pow(.,%g) from %g to %g.\n"
                  "Check your model formulation or use option expr/" POWEXPRHDLR_NAME "/minzerodistance to avoid this warning.\n",
                  exponent, interval.inf, exprhdlrdata->minzerodistance);
               SCIPinfoMessage(scip, NULL, "Expression: ");
               SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
               SCIPinfoMessage(scip, NULL, "\n");
               exprhdlrdata->warnedonpole = TRUE;
            }
            interval.inf = exprhdlrdata->minzerodistance;
         }
         else if( interval.sup < exprhdlrdata->minzerodistance && interval.sup > -exprhdlrdata->minzerodistance )
         {
            if( !exprhdlrdata->warnedonpole && SCIPgetVerbLevel(scip) > SCIP_VERBLEVEL_NONE )
            {
               SCIPinfoMessage(scip, NULL, "Changing lower bound for child of pow(.,%g) from %g to %g.\n"
                  "Check your model formulation or use option expr/" POWEXPRHDLR_NAME "/minzerodistance to avoid this warning.\n",
                  exponent, interval.sup, -exprhdlrdata->minzerodistance);
               SCIPinfoMessage(scip, NULL, "Expression: ");
               SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
               SCIPinfoMessage(scip, NULL, "\n");
               exprhdlrdata->warnedonpole = TRUE;
            }
            interval.sup = -exprhdlrdata->minzerodistance;
         }
      }
   }

   SCIPdebugMsgPrint(scip, " -> [%.15g,%.15g]\n", interval.inf, interval.sup);

   childrenbounds[0] = interval;

   return SCIP_OKAY;
}

/** initial estimates callback for a power expression */
static
SCIP_DECL_EXPRINITESTIMATES(initestimatesPow)
{
   SCIP_EXPRDATA* exprdata;
   SCIP_EXPR* child;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real exponent;
   SCIP_Bool isinteger;
   SCIP_Bool branchcand;
   SCIP_Bool success;
   SCIP_Bool islocal;
   SCIP_Real refpointsunder[3] = {SCIP_INVALID, SCIP_INVALID, SCIP_INVALID};
   SCIP_Real refpointsover[3] = {SCIP_INVALID, SCIP_INVALID, SCIP_INVALID};
   SCIP_Bool overest[6] = {FALSE, FALSE, FALSE, TRUE, TRUE, TRUE};
   int i;

   assert(scip != NULL);
   assert(expr != NULL);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   childlb = bounds[0].inf;
   childub = bounds[0].sup;

   /* if child is essentially constant, then there should be no point in separation */
   if( SCIPisEQ(scip, childlb, childub) )
   {
      SCIPdebugMsg(scip, "skip initestimates as child seems essentially fixed [%.15g,%.15g]\n", childlb, childub);
      return SCIP_OKAY;
   }

   exprdata = SCIPexprGetData(expr);
   exponent = exprdata->exponent;
   assert(exponent != 1.0 && exponent != 0.0); /* this should have been simplified */

   isinteger = EPSISINT(exponent, 0.0);

   /* if exponent is not integral, then child must be non-negative */
   if( !isinteger && childlb < 0.0 )
   {
      /* somewhere we should have tightened the bound on x, but small tightening are not always applied by SCIP
       * it is ok to do this tightening here, but let's assert that we were close to 0.0 already
       */
      assert(SCIPisFeasZero(scip, childlb));
      childlb = 0.0;
   }
   assert(isinteger || childlb >= 0.0);

   /* TODO simplify to get 3 refpoints for either under- or overestimate */
   SCIP_CALL( chooseRefpointsPow(scip, exprdata, childlb, childub, refpointsunder, refpointsover, !overestimate,
         overestimate) );

   for( i = 0; i < 6 && *nreturned < SCIP_EXPR_MAXINITESTIMATES; ++i )
   {
      SCIP_Real refpoint;

      if( (overest[i] && !overestimate) || (!overest[i] && overestimate) )
         continue;

      assert(overest[i] || i < 3); /* make sure that no out-of-bounds array access will be attempted */
      refpoint = overest[i] ? refpointsover[i % 3] : refpointsunder[i % 3];

      if( refpoint == SCIP_INVALID )
         continue;

      assert(SCIPisLE(scip, refpoint, childub) && SCIPisGE(scip, refpoint, childlb));

      branchcand = TRUE;
      SCIP_CALL( buildPowEstimator(scip, exprdata, overest[i], childlb, childub, childlb, childub,
            SCIPexprIsIntegral(child), refpoint, exponent, coefs[*nreturned], &constant[*nreturned],
            &success, &islocal, &branchcand) );

      if( success )
      {
         SCIPdebugMsg(scip, "initestimate x^%g for base in [%g,%g] at ref=%g, over:%u -> %g*x+%g\n", exponent,
               childlb, childub, refpoint, overest[i], coefs[*nreturned][0], constant[*nreturned]);
         ++*nreturned;
      }
   }

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_EXPRHASH(hashPow)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   /* TODO include exponent into hashkey */
   *hashkey = POWEXPRHDLR_HASHKEY;
   *hashkey ^= childrenhashes[0];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvaturePow)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_INTERVAL childinterval;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(exprcurvature != SCIP_EXPRCURV_UNKNOWN);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   exponent = SCIPgetExponentExprPow(expr);
   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   SCIP_CALL( SCIPevalExprActivity(scip, child) );
   childinterval = SCIPexprGetActivity(child);

   *childcurv = SCIPexprcurvPowerInv(childinterval, exponent, exprcurvature);
   /* SCIPexprcurvPowerInv return unknown actually means that curv cannot be obtained */
   *success = *childcurv != SCIP_EXPRCURV_UNKNOWN;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityPow)
{  /*lint --e{715}*/
   SCIP_INTERVAL interval;
   SCIP_Real exponent;
   SCIP_Real inf;
   SCIP_Real sup;
   SCIP_Bool expisint;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(childidx == 0);

   assert(SCIPexprGetChildren(expr)[0] != NULL);
   SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(expr)[0]) );
   interval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   *result = SCIP_MONOTONE_UNKNOWN;
   inf = SCIPintervalGetInf(interval);
   sup = SCIPintervalGetSup(interval);
   exponent = SCIPgetExponentExprPow(expr);
   expisint = EPSISINT(exponent, 0.0); /*lint !e835*/

   if( expisint )
   {
      SCIP_Bool expisodd = ceil(exponent/2) != exponent/2;

      if( expisodd )
      {
         /* x^1, x^3, ... */
         if( exponent >= 0.0 )
            *result = SCIP_MONOTONE_INC;

         /* ..., x^-3, x^-1 are decreasing if 0 is not in ]inf,sup[ */
         else if( inf >= 0.0 || sup <= 0.0 )
            *result = SCIP_MONOTONE_DEC;
      }
      /* ..., x^-4, x^-2, x^2, x^4, ... */
      else
      {
         /* function is not monotone if 0 is in ]inf,sup[ */
         if( inf >= 0.0 )
            *result = exponent >= 0.0 ? SCIP_MONOTONE_INC : SCIP_MONOTONE_DEC;
         else if( sup <= 0.0 )
            *result = exponent >= 0.0 ? SCIP_MONOTONE_DEC : SCIP_MONOTONE_INC;
      }
   }
   else
   {
      /* note that the expression is not defined for negative input values
       *
       * - increasing iff exponent >= 0
       * - decreasing iff exponent <= 0
       */
      *result = exponent >= 0.0 ? SCIP_MONOTONE_INC : SCIP_MONOTONE_DEC;
   }

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityPow)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_Real exponent;
   SCIP_Bool expisint;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   *isintegral = FALSE;

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /* expression can not be integral if child is not */
   if( !SCIPexprIsIntegral(child) )
      return SCIP_OKAY;

   exponent = SCIPgetExponentExprPow(expr);
   assert(exponent != 0.0);
   expisint = EPSISINT(exponent, 0.0); /*lint !e835*/

   /* expression is integral if and only if exponent non-negative and integral */
   *isintegral = expisint && exponent >= 0.0;

   return SCIP_OKAY;
}

/** simplifies a signpower expression
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifySignpower)
{  /*lint --e{715}*/
   SCIP_EXPR* base;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   base = SCIPexprGetChildren(expr)[0];
   assert(base != NULL);

   exponent = SCIPgetExponentExprPow(expr);
   SCIPdebugPrintf("[simplifySignpower] simplifying power with expo %g\n", exponent);
   assert(exponent >= 1.0);

   /* enforces SPOW2 */
   if( exponent == 1.0 )
   {
      SCIPdebugPrintf("[simplifySignpower] POW2\n");
      *simplifiedexpr = base;
      SCIPcaptureExpr(*simplifiedexpr);
      return SCIP_OKAY;
   }

   /* enforces SPOW3 */
   if( SCIPisExprValue(scip, base) )
   {
      SCIP_Real baseval;

      SCIPdebugPrintf("[simplifySignpower] POW3\n");
      baseval = SCIPgetValueExprValue(base);

      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, SIGN(baseval) * pow(REALABS(baseval), exponent),
                 ownercreate, ownercreatedata) );

      return SCIP_OKAY;
   }

   /* enforces SPOW11 (exp(x)^n = exp(n*x))
    * since exp() is always nonnegative, we can treat signpower as normal power here
    */
   if( SCIPisExprExp(scip, base) )
   {
      SCIP_EXPR* child;
      SCIP_EXPR* prod;
      SCIP_EXPR* exponential;
      SCIP_EXPR* simplifiedprod;

      SCIPdebugPrintf("[simplifySignpower] POW11\n");
      child = SCIPexprGetChildren(base)[0];

      /* multiply child of exponential with exponent */
      SCIP_CALL( SCIPcreateExprProduct(scip, &prod, 1, &child, exponent, ownercreate, ownercreatedata) );

      /* simplify product */
      SCIP_CALL( SCIPcallExprSimplify(scip, prod, &simplifiedprod, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &prod) );

      /* create exponential with new child */
      SCIP_CALL( SCIPcreateExprExp(scip, &exponential, simplifiedprod, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedprod) );

      /* the final simplified expression is the simplification of the just created exponential */
      SCIP_CALL( SCIPcallExprSimplify(scip, exponential, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exponential) );

      return SCIP_OKAY;
   }

   /* enforces SPOW6 */
   if( EPSISINT(exponent, 0.0) && ((int)exponent) % 2 == 1 )
   {
      SCIP_EXPR* aux;

      /* we do not just change the expression data of expression to say it is a normal power, since, at the moment,
       * simplify identifies that expressions changed by checking that the pointer of the input expression is
       * different from the returned (simplified) expression
      */
      SCIP_CALL( SCIPcreateExprPow(scip, &aux, base, exponent, ownercreate, ownercreatedata) );

      SCIP_CALL( simplifyPow(scip, aux, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &aux) );

      return SCIP_OKAY;
   }

   /* enforces SPOW10 */
   if( SCIPisExprVar(scip, base) )
   {
      SCIP_VAR* basevar;

      SCIPdebugPrintf("[simplifySignpower] POW10\n");
      basevar = SCIPgetVarExprVar(base);

      assert(basevar != NULL);

      if( SCIPvarIsBinary(basevar) )
      {
         *simplifiedexpr = base;
         SCIPcaptureExpr(*simplifiedexpr);
         return SCIP_OKAY;
      }
   }

   /* TODO if( SCIPisExprSignpower(scip, base) ... */

   /* enforces SPOW8
    * given (signpow n (pow expo expr)) we distribute the exponent:
    * -> (signpow n*expo expr) for even n  (i.e., sign(x^n) * |x|^n = 1 * x^n)
    * notes: n is an even integer (see SPOW6 above)
    * FIXME: doesn't this extend to any exponent?
    * If (pow expo expr) can be negative, it should mean that (-1)^expo = -1
    * then (signpow n (pow expo expr)) = sign(expr^expo) * |expr^expo|^n
    * then sign(expr^expo) = sign(expr) and |expr^expo| = |expr|^expo and so
    * (signpow n (pow expo expr)) = sign(expr^expo) * |expr^expo|^n = sign(expr) * |expr|^(expo*n) = signpow n*expo expr
    */
   if( EPSISINT(exponent, 0.0) && SCIPisExprPower(scip, base) )
   {
      SCIP_EXPR* aux;
      SCIP_Real newexponent;

      assert(((int)exponent) % 2 == 0 ); /* odd case should have been handled by SPOW6 */

      newexponent = SCIPgetExponentExprPow(base) * exponent;
      SCIP_CALL( SCIPcreateExprSignpower(scip, &aux, SCIPexprGetChildren(base)[0], newexponent,
                 ownercreate, ownercreatedata) );
      SCIP_CALL( simplifySignpower(scip, aux, simplifiedexpr, ownercreate, ownercreatedata) );

      SCIP_CALL( SCIPreleaseExpr(scip, &aux) );

      return SCIP_OKAY;
   }

   /* enforces SPOW9 */
   if( SCIPisExprSum(scip, base)
      && SCIPexprGetNChildren(base) == 1
      && SCIPgetConstantExprSum(base) == 0.0 )
   {
      SCIP_EXPR* simplifiedaux;
      SCIP_EXPR* aux;
      SCIP_Real newcoef;

      SCIPdebugPrintf("[simplifySignpower] seeing a sum with one term, exponent %g\n", exponent);
      /* assert SS7 holds */
      assert(SCIPgetCoefsExprSum(base)[0] != 1.0);

      /* create (signpow n expr) and simplify it
       * note: we call simplifySignpower directly, since we know that `expr` is simplified */
      SCIP_CALL( SCIPcreateExprSignpower(scip, &aux, SCIPexprGetChildren(base)[0], exponent,
                 ownercreate, ownercreatedata) );
      newcoef = SIGN(SCIPgetCoefsExprSum(base)[0]) * pow(REALABS(SCIPgetCoefsExprSum(base)[0]), exponent);
      SCIP_CALL( simplifySignpower(scip, aux, &simplifiedaux, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &aux) );

      /* create (sum (signpow n expr)) and simplify it
       * this calls simplifySum directly, since we know its child is simplified */
      SCIP_CALL( SCIPcreateExprSum(scip, &aux, 1, &simplifiedaux, &newcoef, 0.0, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPcallExprSimplify(scip, aux, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &aux) );
      SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedaux) );
      return SCIP_OKAY;
   }

   SCIPdebugPrintf("[simplifySignpower] signpower is simplified\n");
   *simplifiedexpr = expr;

   /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
   SCIPcaptureExpr(*simplifiedexpr);

   return SCIP_OKAY;
}

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrSignpower)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrSignpower(scip) );

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_EXPRPRINT(printSignpower)
{  /*lint --e{715}*/
   assert(expr != NULL);

   switch( stage )
   {
      case SCIP_EXPRITER_ENTEREXPR :
      {
         SCIPinfoMessage(scip, file, "signpower(");
         break;
      }

      case SCIP_EXPRITER_VISITINGCHILD :
      {
         assert(currentchild == 0);
         break;
      }

      case SCIP_EXPRITER_LEAVEEXPR :
      {
         SCIPinfoMessage(scip, file, ",%g)", SCIPgetExponentExprPow(expr));
         break;
      }

      case SCIP_EXPRITER_VISITEDCHILD :
      default:
         break;
   }

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_EXPRPARSE(parseSignpower)
{  /*lint --e{715}*/
   SCIP_EXPR* childexpr;
   SCIP_Real exponent;

   assert(expr != NULL);

   /**! [SnippetExprParseSignpower] */
   /* parse child expression string */
   SCIP_CALL( SCIPparseExpr(scip, &childexpr, string, endstring, ownercreate, ownercreatedata) );
   assert(childexpr != NULL);

   string = *endstring;
   while( *string == ' ' )
      ++string;

   if( *string != ',' )
   {
      SCIPerrorMessage("Expected comma after first argument of signpower().\n");
      return SCIP_READERROR;
   }
   ++string;

   if( !SCIPparseReal(scip, string, &exponent, (char**)endstring) )
   {
      SCIPerrorMessage("Expected numeric exponent for second argument of signpower().\n");
      return SCIP_READERROR;
   }

   if( exponent <= 1.0 || SCIPisInfinity(scip, exponent) )
   {
      SCIPerrorMessage("Expected finite exponent >= 1.0 for signpower().\n");
      return SCIP_READERROR;
   }

   /* create signpower expression */
   SCIP_CALL( SCIPcreateExprSignpower(scip, expr, childexpr, exponent, ownercreate, ownercreatedata) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the signpower expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &childexpr) );

   *success = TRUE;
   /**! [SnippetExprParseSignpower] */

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalSignpower)
{  /*lint --e{715}*/
   SCIP_Real exponent;
   SCIP_Real base;

   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID);

   exponent = SCIPgetExponentExprPow(expr);
   base = SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]);

   *val = SIGN(base) * pow(REALABS(base), exponent);

   /* if there is a range error, pow() should return some kind of infinity, or HUGE_VAL
    * we could also work with floating point exceptions or errno, but I am not sure this would be thread-safe
    */
   if( !SCIPisFinite(*val) || *val == HUGE_VAL || *val == -HUGE_VAL )
      *val = SCIP_INVALID;

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffSignpower)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_Real childval;
   SCIP_Real exponent;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);
   assert(childidx == 0);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(child)), "val") != 0);

   childval = SCIPexprGetEvalValue(child);
   assert(childval != SCIP_INVALID);

   exponent = SCIPgetExponentExprPow(expr);
   assert(exponent >= 1.0);

   *val = exponent * pow(REALABS(childval), exponent - 1.0);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalSignpower)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
   {
      SCIPintervalSetEmpty(interval);
      return SCIP_OKAY;
   }

   SCIPintervalSignPowerScalar(SCIP_INTERVAL_INFINITY, interval, childinterval, SCIPgetExponentExprPow(expr));

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_EXPRESTIMATE(estimateSignpower)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real childglb;
   SCIP_Real childgub;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), SIGNPOWEXPRHDLR_NAME) == 0);
   assert(refpoint != NULL);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);  /* the default */
   assert(success != NULL);

   *success = FALSE;

   SCIPdebugMsg(scip, "%sestimation of signed x^%g at x=%g\n", overestimate ? "over" : "under",
         SCIPgetExponentExprPow(expr), *refpoint);

   /* we can not generate a cut at +/- infinity */
   if( SCIPisInfinity(scip, REALABS(*refpoint)) )
      return SCIP_OKAY;

   childlb = localbounds[0].inf;
   childub = localbounds[0].sup;

   childglb = globalbounds[0].inf;
   childgub = globalbounds[0].sup;

   exprdata = SCIPexprGetData(expr);
   exponent = exprdata->exponent;
   assert(exponent > 1.0); /* exponent == 1 should have been simplified */

   /* if child is constant, then return a constant estimator
    * this can help with small infeasibilities if boundtightening is relaxing bounds too much
    */
   if( childlb == childub )
   {
      *coefs = 0.0;
      *constant = SIGN(childlb)*pow(REALABS(childlb), exponent);
      *success = TRUE;
      *islocal = childglb != childgub;
      *branchcand = FALSE;
      return SCIP_OKAY;
   }

   if( childlb >= 0.0 )
   {
      estimateParabola(scip, exponent, overestimate, childlb, childub, MAX(0.0, *refpoint), constant, coefs,
            islocal, success);

      *branchcand = *islocal;

      /* if odd or signed power, then check whether tangent on parabola is also globally valid, that is
       * reference point is right of -root*global-lower-bound
       */
      if( !*islocal && childglb < 0.0 )
      {
         if( SCIPisInfinity(scip, -childglb) )
            *islocal = TRUE;
         else
         {
            if( exprdata->root == SCIP_INVALID )
            {
               SCIP_CALL( computeSignpowerRoot(scip, &exprdata->root, exponent) );
            }
            *islocal = *refpoint < exprdata->root * (-childglb);
         }
      }
   }
   else  /* and childlb < 0.0 due to previous if */
   {
      /* compute root if not known yet; only needed if mixed sign (global child ub > 0) */
      if( exprdata->root == SCIP_INVALID && childgub > 0.0 )
      {
         SCIP_CALL( computeSignpowerRoot(scip, &exprdata->root, exponent) );
      }
      estimateSignedpower(scip, exponent, exprdata->root, overestimate, childlb, childub, *refpoint,
         childglb, childgub, constant, coefs, islocal, branchcand, success);
   }

   return SCIP_OKAY;
}

/** initial estimates callback for a signpower expression */
static
SCIP_DECL_EXPRINITESTIMATES(initestimatesSignpower)
{
   SCIP_EXPRDATA* exprdata;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real exponent;
   SCIP_Bool branchcand;
   SCIP_Bool success;
   SCIP_Bool islocal;
   SCIP_Real refpointsunder[3] = {SCIP_INVALID, SCIP_INVALID, SCIP_INVALID};
   SCIP_Real refpointsover[3] = {SCIP_INVALID, SCIP_INVALID, SCIP_INVALID};
   SCIP_Bool overest[6] = {FALSE, FALSE, FALSE, TRUE, TRUE, TRUE};
   SCIP_Real refpoint;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), SIGNPOWEXPRHDLR_NAME) == 0);

   childlb = bounds[0].inf;
   childub = bounds[0].sup;

   /* if child is essentially constant, then there should be no point in separation */
   if( SCIPisEQ(scip, childlb, childub) )
   {
      SCIPdebugMsg(scip, "skip initestimates as child seems essentially fixed [%.15g,%.15g]\n", childlb, childub);
      return SCIP_OKAY;
   }

   exprdata = SCIPexprGetData(expr);
   exponent = exprdata->exponent;
   assert(exponent > 1.0); /* this should have been simplified */

   if( childlb >= 0.0 )
   {
      if( !overestimate )
         addTangentRefpoints(scip, exponent, childlb, childub, refpointsunder);
      if( overestimate && !SCIPisInfinity(scip, childub) )
         refpointsover[0] = (childlb + childub) / 2.0;
   }
   else if( childub <= 0.0 )
   {
      if( !overestimate && !SCIPisInfinity(scip, -childlb) )
         refpointsunder[0] = (childlb + childub) / 2.0;
      if( overestimate )
         addTangentRefpoints(scip, exponent, childlb, childub, refpointsunder);
   }
   else
   {
      SCIP_CALL( addSignpowerRefpoints(scip, exprdata, childlb, childub, exponent, !overestimate, refpointsunder) );
   }

   /* add cuts for all refpoints */
   for( i = 0; i < 6 && *nreturned < SCIP_EXPR_MAXINITESTIMATES; ++i )
   {
      if( (overest[i] && !overestimate) || (!overest[i] && overestimate) )
         continue;

      assert(overest[i] || i < 3); /* make sure that no out-of-bounds array access will be attempted */
      refpoint = overest[i] ? refpointsover[i % 3] : refpointsunder[i % 3];
      if( refpoint == SCIP_INVALID )
         continue;
      assert(SCIPisLE(scip, refpoint, childub) && SCIPisGE(scip, refpoint, childlb));

      if( childlb >= 0 )
      {
         estimateParabola(scip, exponent, overest[i], childlb, childub, refpoint, &constant[*nreturned], coefs[*nreturned],
               &islocal, &success);
      }
      else
      {
         /* compute root if not known yet; only needed if mixed sign (global child ub > 0) */
         if( exprdata->root == SCIP_INVALID && childub > 0.0 )
         {
            SCIP_CALL( computeSignpowerRoot(scip, &exprdata->root, exponent) );
         }
         branchcand = TRUE;
         estimateSignedpower(scip, exponent, exprdata->root, overest[i], childlb, childub, refpoint,
               childlb, childub, &constant[*nreturned], coefs[*nreturned], &islocal,
               &branchcand, &success);
      }

      if( success )
         ++*nreturned;
   }

   return SCIP_OKAY;
}

/** expression reverse propagaton callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropSignpower)
{  /*lint --e{715}*/
   SCIP_INTERVAL interval;
   SCIP_INTERVAL exprecip;
   SCIP_Real exponent;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   exponent = SCIPgetExponentExprPow(expr);

   SCIPdebugMsg(scip, "reverseprop signpow(x,%g) in [%.15g,%.15g]", exponent, bounds.inf, bounds.sup);

   if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, bounds) )
   {
      SCIPdebugMsgPrint(scip, "-> no improvement\n");
      return SCIP_OKAY;
   }

   /* f = pow(c0, alpha) -> c0 = pow(f, 1/alpha) */
   SCIPintervalSet(&exprecip, exponent);
   SCIPintervalReciprocal(SCIP_INTERVAL_INFINITY, &exprecip, exprecip);
   if( exprecip.inf == exprecip.sup )
   {
      SCIPintervalSignPowerScalar(SCIP_INTERVAL_INFINITY, &interval, bounds, exprecip.inf);
   }
   else
   {
      SCIP_INTERVAL interval1, interval2;
      SCIPintervalSignPowerScalar(SCIP_INTERVAL_INFINITY, &interval1, bounds, exprecip.inf);
      SCIPintervalSignPowerScalar(SCIP_INTERVAL_INFINITY, &interval2, bounds, exprecip.sup);
      SCIPintervalUnify(&interval, interval1, interval2);
   }

   SCIPdebugMsgPrint(scip, " -> [%.15g,%.15g]\n", interval.inf, interval.sup);

   childrenbounds[0] = interval;

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_EXPRHASH(hashSignpower)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   /* TODO include exponent into hashkey */
   *hashkey = SIGNPOWEXPRHDLR_HASHKEY;
   *hashkey ^= childrenhashes[0];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureSignpower)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_INTERVAL childinterval;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(exprcurvature != SCIP_EXPRCURV_UNKNOWN);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   SCIP_CALL( SCIPevalExprActivity(scip, child) );
   childinterval = SCIPexprGetActivity(child);

   if( exprcurvature == SCIP_EXPRCURV_CONVEX )
   {
      /* signpower is only convex if argument is convex and non-negative */
      *childcurv = SCIP_EXPRCURV_CONVEX;
      *success = childinterval.inf >= 0.0;
   }
   else if( exprcurvature == SCIP_EXPRCURV_CONCAVE )
   {
      /* signpower is only concave if argument is concave and non-positive */
      *childcurv = SCIP_EXPRCURV_CONCAVE;
      *success = childinterval.sup <= 0.0;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicitySignpower)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);

   *result = SCIP_MONOTONE_INC;
   return SCIP_OKAY;
}

/** creates the handler for power expression and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrPow(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_EXPRHDLRDATA* exprhdlrdata;

   SCIP_CALL( SCIPallocClearBlockMemory(scip, &exprhdlrdata) );

   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, POWEXPRHDLR_NAME, POWEXPRHDLR_DESC, POWEXPRHDLR_PRECEDENCE,
              evalPow, exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrPow, freehdlrPow);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataPow, freedataPow);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyPow);
   SCIPexprhdlrSetPrint(exprhdlr, printPow);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalPow);
   SCIPexprhdlrSetEstimate(exprhdlr, initestimatesPow, estimatePow);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropPow);
   SCIPexprhdlrSetHash(exprhdlr, hashPow);
   SCIPexprhdlrSetCompare(exprhdlr, comparePow);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffPow, fwdiffPow, bwfwdiffPow);
   SCIPexprhdlrSetCurvature(exprhdlr, curvaturePow);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityPow);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityPow);

   SCIP_CALL( SCIPaddRealParam(scip, "expr/" POWEXPRHDLR_NAME "/minzerodistance",
      "minimal distance from zero to enforce for child in bound tightening",
      &exprhdlrdata->minzerodistance, FALSE, SCIPepsilon(scip), 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates the handler for signed power expression and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrSignpower(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, SIGNPOWEXPRHDLR_NAME, SIGNPOWEXPRHDLR_DESC,
              SIGNPOWEXPRHDLR_PRECEDENCE, evalSignpower, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrSignpower, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataPow, freedataPow);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifySignpower);
   SCIPexprhdlrSetPrint(exprhdlr, printSignpower);
   SCIPexprhdlrSetParse(exprhdlr, parseSignpower);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalSignpower);
   SCIPexprhdlrSetEstimate(exprhdlr, initestimatesSignpower, estimateSignpower);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropSignpower);
   SCIPexprhdlrSetHash(exprhdlr, hashSignpower);
   SCIPexprhdlrSetCompare(exprhdlr, comparePow);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffSignpower, NULL, NULL);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureSignpower);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicitySignpower);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityPow);

   return SCIP_OKAY;
}

/** creates a power expression */
SCIP_RETCODE SCIPcreateExprPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_Real             exponent,           /**< exponent of the power expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(child != NULL);

   SCIP_CALL( createData(scip, &exprdata, exponent) );
   assert(exprdata != NULL);

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPgetExprhdlrPower(scip), exprdata, 1, &child, ownercreate,
              ownercreatedata) );

   return SCIP_OKAY;
}

/** creates a signpower expression */
SCIP_RETCODE SCIPcreateExprSignpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_Real             exponent,           /**< exponent of the power expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindExprhdlr(scip, SIGNPOWEXPRHDLR_NAME) != NULL);

   SCIP_CALL( createData(scip, &exprdata, exponent) );
   assert(exprdata != NULL);

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPfindExprhdlr(scip, SIGNPOWEXPRHDLR_NAME), exprdata, 1, &child,
              ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/** indicates whether expression is of signpower-type */  /*lint -e{715}*/
SCIP_Bool SCIPisExprSignpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{ /*lint --e{715}*/
   assert(expr != NULL);

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), SIGNPOWEXPRHDLR_NAME) == 0;
}

/** computes coefficients of linearization of a square term in a reference point */
void SCIPaddSquareLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             refpoint,           /**< point where to linearize */
   SCIP_Bool             isint,              /**< whether corresponding variable is a discrete variable, and thus linearization could be moved */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( sqrcoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, REALABS(refpoint)) )
   {
      *success = FALSE;
      return;
   }

   if( !isint || SCIPisIntegral(scip, refpoint) )
   {
      SCIP_Real tmp;

      /* sqrcoef * x^2  ->  tangent in refpoint = sqrcoef * 2 * refpoint * (x - refpoint) */

      tmp = sqrcoef * refpoint;

      if( SCIPisInfinity(scip, 2.0 * REALABS(tmp)) )
      {
         *success = FALSE;
         return;
      }

      *lincoef += 2.0 * tmp;
      tmp *= refpoint;
      *linconstant -= tmp;
   }
   else
   {
      /* sqrcoef * x^2 -> secant between f=floor(refpoint) and f+1 = sqrcoef * (f^2 + ((f+1)^2 - f^2) * (x-f))
       * = sqrcoef * (-f*(f+1) + (2*f+1)*x)
       */
      SCIP_Real f;
      SCIP_Real coef;
      SCIP_Real constant;

      f = SCIPfloor(scip, refpoint);

      coef     =  sqrcoef * (2.0 * f + 1.0);
      constant = -sqrcoef * f * (f + 1.0);

      if( SCIPisInfinity(scip, REALABS(coef)) || SCIPisInfinity(scip, REALABS(constant)) )
      {
         *success = FALSE;
         return;
      }

      *lincoef     += coef;
      *linconstant += constant;
   }
}

/** computes coefficients of secant of a square term */
void SCIPaddSquareSecant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             lb,                 /**< lower bound on variable */
   SCIP_Real             ub,                 /**< upper bound on variable */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of secant */
   SCIP_Real*            linconstant,        /**< buffer to add constant of secant */
   SCIP_Bool*            success             /**< buffer to set to FALSE if secant has failed due to large numbers or unboundedness */
   )
{
   SCIP_Real coef;
   SCIP_Real constant;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lb));
   assert(!SCIPisInfinity(scip, -ub));
   assert(SCIPisLE(scip, lb, ub));
   assert(lincoef != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( sqrcoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
   {
      /* unboundedness */
      *success = FALSE;
      return;
   }

   /* sqrcoef * x^2 -> sqrcoef * (lb * lb + (ub*ub - lb*lb)/(ub-lb) * (x-lb)) = sqrcoef * (lb*lb + (ub+lb)*(x-lb))
    *  = sqrcoef * ((lb+ub)*x - lb*ub)
    */
   coef     =  sqrcoef * (lb + ub);
   constant = -sqrcoef * lb * ub;
   if( SCIPisInfinity(scip, REALABS(coef)) || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoef     += coef;
   *linconstant += constant;
}

/* from pub_expr.h */

/** gets the exponent of a power or signed power expression */  /*lint -e{715}*/
SCIP_Real SCIPgetExponentExprPow(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   return exprdata->exponent;
}
