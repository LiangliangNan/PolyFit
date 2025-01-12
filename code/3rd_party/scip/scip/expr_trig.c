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

/**@file   expr_trig.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  handler for sine and cosine expressions
 * @author Fabian Wegscheider
 *
 * The estimator/separator code always computes underestimators for sin(x).
 * For overestimators of cos(x), we first reduce to underestimators of sin(x).
 *
 * Overestimator for sin(x):
 *   Assume that a*y+b <= sin(y) for y in [-ub,-lb].
 *   Then we have a*(-y)-b >= -sin(y) = sin(-y) for y in [-ub,-lb].
 *   Thus, a*x-b >= sin(x) for x in [lb,ub].
 *
 * Underestimator for cos(x):
 *   Assume that a*y+b <= sin(y) for y in [lb+pi/2,ub+pi/2].
 *   Then we have a*(x+pi/2) + b <= sin(x+pi/2) = cos(x) for x in [lb,ub].
 *   Thus, a*x + (b+a*pi/2) <= cos(x) for x in [lb,ub].
 *
 * Overestimator for cos(x):
 *   Assume that a*z+b <= sin(z) for z in [-(ub+pi/2),-(lb+pi/2)].
 *   Then, a*y-b >= sin(y) for y in [lb+pi/2,ub+pi/2].
 *   Then, a*x-b+a*pi/2 >= cos(x) for x in [lb,ub].
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <string.h>
#include <math.h>
#include "scip/expr_trig.h"
#include "scip/expr_value.h"

/* fundamental expression handler properties */
#define SINEXPRHDLR_NAME         "sin"
#define SINEXPRHDLR_DESC         "sine expression"
#define SINEXPRHDLR_PRECEDENCE   91000
#define SINEXPRHDLR_HASHKEY      SCIPcalcFibHash(82457.0)

#define COSEXPRHDLR_NAME         "cos"
#define COSEXPRHDLR_DESC         "cosine expression"
#define COSEXPRHDLR_PRECEDENCE   92000
#define COSEXPRHDLR_HASHKEY      SCIPcalcFibHash(82463.0)

#define MAXCHILDABSVAL        1e+6                       /**< maximum absolute value that is accepted for propagation */
#define NEWTON_NITERATIONS    100
#define NEWTON_PRECISION      1e-12

/*
 * Local methods
 */

/** evaluates the function a*x + b - sin(x) for some coefficient a and constant b at a given point p
 *
 *  the constants a and b are expected to be stored in that order in params
 */
static
SCIP_DECL_NEWTONEVAL(function1)
{  /*lint --e{715}*/
   assert(params != NULL);
   assert(nparams == 2);

   return params[0]*point + params[1] - sin(point);
}

/** evaluates the derivative of a*x + b - sin(x) for some coefficient a and constant b at a given point p
 *
 *  the constants a and b are expected to be stored in that order in params
 */
static
SCIP_DECL_NEWTONEVAL(derivative1)
{  /*lint --e{715}*/
   assert(params != NULL);
   assert(nparams == 2);

   return params[0] - cos(point);
}

/** evaluates the function sin(x) + (alpha - x)*cos(x) - sin(alpha) for some constant alpha at a given point p
 *
 *  the constant alpha is expected to be stored in params
 */
static
SCIP_DECL_NEWTONEVAL(function2)
{  /*lint --e{715}*/
   assert(params != NULL);
   assert(nparams == 1);

   return sin(point) + (params[0] - point) * cos(point) - sin(params[0]);
}

/** evaluates the derivative of sin(x) + (alpha - x)*cos(x) - sin(alpha) for some constant alpha at a given point p
 *
 *  the constant alpha is expected to be stored in params
 */
static
SCIP_DECL_NEWTONEVAL(derivative2)
{  /*lint --e{715}*/
   assert(params != NULL);
   assert(nparams == 1);

   return (point - params[0]) * sin(point);
}

/** helper function to compute the secant if it is a valid underestimator
 *
 *  returns true if the estimator was computed successfully
 */
static
SCIP_Bool computeSecantSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   )
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);
   assert(lb < ub);

   /* if range is too big, secant is not underestimating */
   if( ub - lb >= M_PI )
      return FALSE;

   /* if bounds are not within positive bay, secant is not underestimating */
   if( sin(lb) < 0.0 || sin(ub) < 0.0  || (sin(lb) == 0.0 && cos(lb) < 0.0) )
      return FALSE;

   *lincoef = (sin(ub) - sin(lb)) / (ub - lb);
   *linconst = sin(ub) - (*lincoef) * ub;

   return TRUE;
}

/** helper function to compute the tangent at lower bound if it is underestimating
 *
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeLeftTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Real             lb                  /**< lower bound of argument variable */
   )
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);

   if( SCIPisInfinity(scip, -lb) )
      return FALSE;

   /* left tangent is only underestimating in [pi, 1.5*pi) *2kpi */
   if( sin(lb) > 0.0 || cos(lb) >= 0.0 )
      return FALSE;

   *lincoef = cos(lb);
   *linconst = sin(lb) - (*lincoef) * lb;

   return TRUE;
}

/* TODO: fix this, more cases can be considered, see at unit test
 * the underestimating of the tangents depends not only on the ub but also on the lower bound.
 * right now, this function is only checking whether the tangent underestimates independently of the lower bound!
 */
/** helper function to compute the tangent at upper bound if it is an underestimator
 *
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeRightTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   )
{
   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);

   if( SCIPisInfinity(scip, ub) )
      return FALSE;

   /* right tangent is only underestimating in (1.5*pi, 2*pi] *2kpi */
   if( sin(ub) > 0.0 || cos(ub) <= 0.0 )
      return FALSE;

   *lincoef = cos(ub);
   *linconst = sin(ub) - (*lincoef) * ub;

   return TRUE;
}

/** helper function to compute the tangent at solution point if it is an underestimator
 *
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeSolTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub,                 /**< upper bound of argument variable */
   SCIP_Real             solpoint            /**< solution point to be separated */
   )
{
   SCIP_Real params[2];
   SCIP_Real startingpoints[3];
   SCIP_Real solpointmodpi;
   SCIP_Real intersection;
   int i;

   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);

   /* tangent is only underestimating in negative bay */
   if( sin(solpoint) > 0.0 )
      return FALSE;

   /* compute solution point mod pi */
   solpointmodpi = fmod(solpoint, M_PI);
   if( solpoint < 0.0 )
      solpointmodpi += M_PI;

   /* if the point is too far away from the bounds or is at a multiple of pi, then tangent is not underestimating */
   if( SCIPisGE(scip, solpoint - lb, 2*M_PI) || SCIPisGE(scip, ub - solpoint, 2*M_PI)
      || SCIPisZero(scip, solpointmodpi) )
      return FALSE;

   params[0] = cos(solpoint);
   params[1] = sin(solpoint) - params[0] * solpoint;

   /* choose starting points for Newton procedure */
   if( SCIPisGT(scip, solpointmodpi, M_PI_2) )
   {
      startingpoints[0] = solpoint + (M_PI - solpointmodpi) + M_PI_2;
      startingpoints[1] = startingpoints[0] + M_PI_2;
      startingpoints[2] = startingpoints[1] + M_PI_2;
   }
   else
   {
      startingpoints[0] = solpoint - solpointmodpi - M_PI_2;
      startingpoints[1] = startingpoints[0] - M_PI_2;
      startingpoints[2] = startingpoints[1] - M_PI_2;
   }

   /* use Newton procedure to test if cut is valid */
   for( i = 0; i < 3; ++i )
   {
      intersection = SCIPcalcRootNewton(function1, derivative1, params, 2, startingpoints[i], NEWTON_PRECISION,
         NEWTON_NITERATIONS);

      if( intersection != SCIP_INVALID && !SCIPisEQ(scip, intersection, solpoint) ) /*lint !e777*/
         break;
   }

   /* if Newton failed or intersection point lies within bounds, underestimator is not valid */
   if( intersection == SCIP_INVALID || (intersection >= lb && intersection <= ub) ) /*lint !e777*/
      return FALSE;

   *lincoef = params[0];
   *linconst = params[1];

   return TRUE;
}

/** helper function to compute the secant between lower bound and some point of the graph such that it underestimates
 *
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeLeftSecantSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   )
{
   SCIP_Real lbmodpi;
   SCIP_Real tangentpoint;
   SCIP_Real startingpoint;

   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);
   assert(lb < ub);

   if( SCIPisInfinity(scip, -lb) )
      return FALSE;

   /* compute shifted bounds for case evaluation */
   lbmodpi = fmod(lb, M_PI);
   if( lb < 0.0 )
      lbmodpi += M_PI;

   /* choose starting point for Newton procedure */
   if( cos(lb) < 0.0 )
   {
      /* in [pi/2,pi] underestimating doesn't work; otherwise, take the midpoint of possible area */
      if( SCIPisLE(scip, sin(lb), 0.0) )
         return FALSE;
      else
         startingpoint = lb + 1.25*M_PI - lbmodpi;
   }
   else
   {
      /* in ascending area, take the midpoint of the possible area in descending part */
      /* for lb < 0 but close to zero, we may have sin(lb) = 0 but lbmodpi = pi, which gives a starting point too close to lb
       * but for sin(lb) around 0 we know that the tangent point needs to be in [lb+pi,lb+pi+pi/2]
       */
      if( SCIPisZero(scip, sin(lb)) )
         startingpoint = lb + 1.25*M_PI;
      else if( sin(lb) < 0.0 )
         startingpoint = lb + 2.25*M_PI - lbmodpi;
      else
         startingpoint = lb + 1.25*M_PI - lbmodpi;
   }

   /* use Newton procedure to find the point where the tangent intersects sine at lower bound */
   tangentpoint = SCIPcalcRootNewton(function2, derivative2, &lb, 1, startingpoint, NEWTON_PRECISION,
      NEWTON_NITERATIONS);

   /* if Newton procedure failed, no cut is added */
   if( tangentpoint == SCIP_INVALID ) /*lint !e777*/
      return FALSE;

   /* if the computed point lies outside the bounds, it is shifted to upper bound */
   if( SCIPisGE(scip, tangentpoint, ub) )
   {
      tangentpoint = ub;

      /* check whether affine function is still underestimating */
      if( SCIPisLE(scip, sin(0.5 * (ub + lb)), sin(lb) + 0.5*(sin(ub) - sin(lb))) )
         return FALSE;
   }

   if( SCIPisEQ(scip, tangentpoint, lb) )  /*lint !e777 */
      return FALSE;

   /* compute secant between lower bound and connection point */
   *lincoef = (sin(tangentpoint) - sin(lb)) / (tangentpoint - lb);
   *linconst = sin(lb) - (*lincoef) * lb;

   /* if the bounds are too close to each other, it's possible that the underestimator is not valid */
   if( *lincoef >= cos(lb) )
      return FALSE;

   SCIPdebugMsg(scip, "left secant: %g + %g*x <= sin(x) on [%g,%g]\n", *linconst, *lincoef, lb, ub);

   return TRUE;
}

/** helper function to compute the secant between upper bound and some point of the graph such that it underestimates
 *
 *  returns true if the underestimator was computed successfully
 */
static
SCIP_Bool computeRightSecantSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of tangent */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of tangent */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   )
{
   SCIP_Real ubmodpi;
   SCIP_Real tangentpoint;
   SCIP_Real startingpoint;

   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconst != NULL);
   assert(lb < ub);

   if( SCIPisInfinity(scip, ub) )
      return FALSE;

   /* compute shifted bounds for case evaluation */
   ubmodpi = fmod(ub, M_PI);
   if( ub < 0.0 )
      ubmodpi += M_PI;

   /* choose starting point for Newton procedure */
   if( cos(ub) > 0.0 )
   {
      /* in [3*pi/2,2*pi] underestimating doesn't work; otherwise, take the midpoint of possible area */
      if( SCIPisLE(scip, sin(ub), 0.0) )
         return FALSE;
      else
         startingpoint = ub - M_PI_4 - ubmodpi;
   }
   else
   {
      /* in descending area, take the midpoint of the possible area in ascending part */
      /* for ub < 0 but close to zero, we may have sin(ub) = 0 but ubmodpi = pi, which gives a starting point too close to ub
       * but for sin(ub) around 0 we know that the tangent point needs to be in [ub-(pi+pi/2),ub-pi]
       */
      if( SCIPisZero(scip, sin(ub)) )
         startingpoint = ub - 1.25*M_PI;
      else if( sin(ub) < 0.0 )
         startingpoint = ub - 1.25*M_PI - ubmodpi;
      else
         startingpoint = ub - M_PI_4 - ubmodpi;
   }

   /* use Newton procedure to find the point where the tangent intersects sine at lower bound */
   tangentpoint = SCIPcalcRootNewton(function2, derivative2, &ub, 1, startingpoint, NEWTON_PRECISION,
      NEWTON_NITERATIONS);

   /* if Newton procedure failed, no underestimator is found */
   if( tangentpoint == SCIP_INVALID ) /*lint !e777*/
      return FALSE;

   /* if the computed point lies outside the bounds, it is shifted to upper bound */
   if( SCIPisLE(scip, tangentpoint, lb) )
   {
      tangentpoint = lb;

      /* check whether affine function is still underestimating */
      if( SCIPisLE(scip, sin(0.5 * (ub + lb)), sin(lb) + 0.5*(sin(ub) - sin(lb))) )
         return FALSE;
   }

   if( SCIPisEQ(scip, tangentpoint, ub) )  /*lint !e777 */
      return FALSE;

   /* compute secant between lower bound and connection point */
   *lincoef = (sin(tangentpoint) - sin(ub)) / (tangentpoint - ub);
   *linconst = sin(ub) - (*lincoef) * ub;

   /* if the bounds are to close to each other, it's possible that the underestimator is not valid */
   if( *lincoef <= cos(lb) )
      return FALSE;

   return TRUE;
}

/** helper function to compute the new interval for child in reverse propagation */
static
SCIP_RETCODE computeRevPropIntervalSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTERVAL         parentbounds,       /**< bounds for sine expression */
   SCIP_INTERVAL         childbounds,        /**< bounds for child expression */
   SCIP_INTERVAL*        newbounds           /**< buffer to store new child bounds */
   )
{
   SCIP_Real newinf = childbounds.inf;
   SCIP_Real newsup = childbounds.sup;

   /* if the absolute values of the bounds are too large, skip reverse propagation
    * TODO: if bounds are close but too large, shift them to [0,2pi] and do the computation there
    */
   if( ABS(newinf) > MAXCHILDABSVAL || ABS(newsup) > MAXCHILDABSVAL )
   {
      SCIPintervalSetBounds(newbounds, newinf, newsup);
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -newinf) )
   {
      /* l(x) and u(x) are lower/upper bound of child, l(s) and u(s) are lower/upper bound of sin expr
       *
       * if sin(l(x)) < l(s), we are looking for k minimal s.t. a + 2k*pi > l(x) where a = asin(l(s))
       * then the new lower bound is a + 2k*pi
       */
      if( SCIPisLT(scip, sin(newinf), parentbounds.inf) )
      {
         SCIP_Real a = asin(parentbounds.inf);
         int k = (int) ceil((newinf - a) / (2.0*M_PI));
         newinf = a + 2.0*M_PI * k;
      }

      /* if sin(l(x)) > u(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) where a = asin(u(s))
       * then the new lower bound is pi - a + 2k*pi
       */
      else if( SCIPisGT(scip, sin(newinf), parentbounds.sup) )
      {
         SCIP_Real a = asin(parentbounds.sup);
         int k = (int) ceil((newinf + a) / (2.0*M_PI) - 0.5);
         newinf = M_PI * (2.0*k + 1.0) - a;
      }

      assert(newinf >= childbounds.inf);
      assert(SCIPisFeasGE(scip, sin(newinf), parentbounds.inf));
      assert(SCIPisFeasLE(scip, sin(newinf), parentbounds.sup));
   }

   if( !SCIPisInfinity(scip, newsup) )
   {
      /* if sin(u(x)) > u(s), we are looking for k minimal s.t. a + 2k*pi > u(x) - 2*pi where a = asin(u(s))
       * then the new upper bound is a + 2k*pi
       */
      if ( SCIPisGT(scip, sin(newsup), parentbounds.sup) )
      {
         SCIP_Real a = asin(parentbounds.sup);
         int k = (int) ceil((newsup - a ) / (2.0*M_PI)) - 1;
         newsup = a + 2.0*M_PI * k;
      }

      /* if sin(u(x)) < l(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) - 2*pi where a = asin(l(s))
       * then the new upper bound is pi - a + 2k*pi
       */
      if( SCIPisLT(scip, sin(newsup), parentbounds.inf) )
      {
         SCIP_Real a = asin(parentbounds.inf);
         int k = (int) ceil((newsup + a) / (2.0*M_PI) - 0.5) - 1;
         newsup = M_PI * (2.0*k + 1.0) - a;
      }

      assert(newsup <= childbounds.sup);
      assert(SCIPisFeasGE(scip, sin(newsup), parentbounds.inf));
      assert(SCIPisFeasLE(scip, sin(newsup), parentbounds.sup));
   }

   /* if the new interval is invalid, the old one was already invalid */
   if( newinf <= newsup )
      SCIPintervalSetBounds(newbounds, newinf, newsup);
   else
      SCIPintervalSetEmpty(newbounds);

   return SCIP_OKAY;
}

/** helper function to compute coefficients and constant term of a linear estimator at a given point
 *
 *  The function will try to compute the following estimators in that order:
 *  - soltangent: tangent at specified refpoint
 *  - secant: secant between the points (lb,sin(lb)) and (ub,sin(ub))
 *  - left secant: secant between lower bound and some point of the graph
 *  - right secant: secant between upper bound and some point of the graph
 *
 *  They are ordered such that a successful computation for one of them cannot be improved by following ones in terms
 *  of value at the reference point.
 */
static
SCIP_Bool computeEstimatorsTrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< sin or cos expression */
   SCIP_Real*            lincoef,            /**< buffer to store the linear coefficient */
   SCIP_Real*            linconst,           /**< buffer to store the constant term */
   SCIP_Real             refpoint,           /**< point at which to underestimate (can be SCIP_INVALID) */
   SCIP_Real             childlb,            /**< lower bound of child variable */
   SCIP_Real             childub,            /**< upper bound of child variable */
   SCIP_Bool             underestimate       /**< whether the estimator should be underestimating */
   )
{
   SCIP_Bool success;
   SCIP_Bool iscos;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "sin") == 0
         || strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "cos") == 0);
   assert(SCIPisLE(scip, childlb, childub));

   /* if child is essentially constant, then there should be no point in estimation */
   if( SCIPisEQ(scip, childlb, childub) ) /* @todo maybe return a constant estimator? */
      return FALSE;

   iscos = strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "cos") == 0;

   /* for cos expressions, the bounds have to be shifted before and after computation */
   if( iscos )
   {
      childlb += M_PI_2;
      childub += M_PI_2;
      refpoint += M_PI_2;
   }

   if( !underestimate )
   {
      SCIP_Real tmp = childlb;
      childlb = -childub;
      childub = -tmp;
      refpoint *= -1;
   }

   /* try out tangent at solution point */
   success = computeSolTangentSin(scip, lincoef, linconst, childlb, childub, refpoint);

   /* otherwise, try out secant */
   if( !success )
      success = computeSecantSin(scip, lincoef, linconst, childlb, childub);

   /* otherwise, try left secant */
   if( !success )
      success = computeLeftSecantSin(scip, lincoef, linconst, childlb, childub);

   /* otherwise, try right secant */
   if( !success )
      success = computeRightSecantSin(scip, lincoef, linconst, childlb, childub);

   if( !success )
      return FALSE;

   /* for overestimators, mirror back */
   if( !underestimate )
      (*linconst) *= -1.0;

   /* for cos expressions, shift back */
   if( iscos )
      (*linconst) += (*lincoef) * M_PI_2;

   return TRUE;
}

/** helper function to create initial cuts for sine and cosine separation
 *
 *  The following 5 cuts can be generated:
 *  - secant: secant between the bounds (lb,sin(lb)) and (ub,sin(ub))
 *  - left/right secant: secant between lower/upper bound and some point of the graph
 *  - left/right tangent: tangents at the lower/upper bounds
 */
static
SCIP_RETCODE computeInitialCutsTrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< sin or cos expression */
   SCIP_Real             childlb,            /**< lower bound of child variable */
   SCIP_Real             childub,            /**< upper bound of child variable */
   SCIP_Bool             underestimate,      /**< whether the cuts should be underestimating */
   SCIP_Real**           coefs,              /**< buffer to store coefficients of computed estimators */
   SCIP_Real*            constant,           /**< buffer to store constant of computed estimators */
   int*                  nreturned           /**< buffer to store number of estimators that have been computed */
   )
{
   SCIP_Bool iscos;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "sin") == 0 || strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "cos") == 0);
   assert(SCIPisLE(scip, childlb, childub));

   /* caller must ensure that variable is not already fixed */
   assert(!SCIPisEQ(scip, childlb, childub));

   *nreturned = 0;

   /* for cos expressions, the bounds have to be shifted before and after computation */
   iscos = strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), "cos") == 0;
   if( iscos )
   {
      childlb += M_PI_2;
      childub += M_PI_2;
   }

   /*
    * Compute all initial cuts
    * For each linear equation z = a*x + b with bounds [lb,ub] the parameters can be computed by:
    *
    * a = cos(x^)    and     b = sin(x^) - a * x^        where x^ is any known point in [lb,ub]
    *
    * and the resulting cut is       a*x + b <=/>= z           depending on over-/underestimation
    */

   if( ! underestimate )
   {
      SCIP_Real aux;
      aux = childlb;
      childlb = -childub;
      childub = -aux;
   }

   /* if we can generate a secant between the bounds, then we have convex (concave) hull */
   if( computeSecantSin(scip, coefs[*nreturned], &constant[*nreturned], childlb, childub) )
      (*nreturned)++;
   else
   {
      /* try generating a secant between lb (ub) and some point < ub (> lb); otherwise try with tangent at lb (ub)*/
      if( computeLeftSecantSin(scip, coefs[*nreturned], &constant[*nreturned], childlb, childub) )
         (*nreturned)++;
      else if( computeLeftTangentSin(scip, coefs[*nreturned], &constant[*nreturned], childlb) )
         (*nreturned)++;

      /* try generating a secant between ub (lb) and some point > lb (< ub); otherwise try with tangent at ub (lb)*/
      if( computeRightSecantSin(scip, coefs[*nreturned], &constant[*nreturned], childlb, childub) )
         (*nreturned)++;
      else if( computeRightTangentSin(scip, coefs[*nreturned], &constant[*nreturned], childub) )
         (*nreturned)++;
   }

   /* for cos expressions, the estimator needs to be shifted back to match original bounds */
   for( i = 0; i < *nreturned; ++i )
   {
      if( ! underestimate )
         constant[i] *= -1.0;

      if( iscos)
      {
         constant[i] += coefs[i][0] * M_PI_2;
      }
   }

   return SCIP_OKAY;
}

/* helper function that computes the curvature of a sine expression for given bounds and curvature of child */
static
SCIP_EXPRCURV computeCurvatureSin(
   SCIP_EXPRCURV         childcurvature,     /**< curvature of child */
   SCIP_Real             lb,                 /**< lower bound of child */
   SCIP_Real             ub                  /**< upper bound of child */
   )
{
   SCIP_Real lbsin = sin(lb);
   SCIP_Real ubsin = sin(ub);
   SCIP_Real lbcos = cos(lb);
   SCIP_Real ubcos = cos(ub);

   /* curvature can only be determined if bounds lie within one bay*/
   if( (ub - lb <= M_PI) && (lbsin * ubsin >= 0.0) )
   {
      /* special case that both sin(ub) and sin(lb) are 0 (i.e. ub - lb = pi) */
      if( lbsin == 0.0 && ubsin == 0.0 )
      {
         if( childcurvature == SCIP_EXPRCURV_LINEAR )
            return (fmod(lb, 2.0*M_PI) == 0.0) ? SCIP_EXPRCURV_CONCAVE : SCIP_EXPRCURV_CONVEX;
      }

      /* if sine is monotone on the interval, the curvature depends on the child curvature and on the segment */
      else if( lbcos * ubcos >= 0.0 )
      {
         /* on [0, pi/2], sine is concave iff child is concave */
         if( lbsin >= 0.0 && lbcos >= 0.0 && ((int)(childcurvature & SCIP_EXPRCURV_CONCAVE) != 0))
            return SCIP_EXPRCURV_CONCAVE;

         /* on [pi/2, pi], sine is concave iff child is convex */
         if( lbsin >= 0.0 && lbcos <= 0.0 && ((int)(childcurvature & SCIP_EXPRCURV_CONVEX) != 0))
            return SCIP_EXPRCURV_CONCAVE;

         /* on [pi, 3pi/2], sine is convex iff child is concave */
         if( lbsin <= 0.0 && lbcos <= 0.0 && ((int)(childcurvature & SCIP_EXPRCURV_CONCAVE) != 0))
            return SCIP_EXPRCURV_CONVEX;

         /* on [3pi/2, 2pi], sine is convex iff child is convex */
         if( lbsin <= 0.0 && lbcos >= 0.0 && ((int)(childcurvature & SCIP_EXPRCURV_CONVEX) != 0))
            return SCIP_EXPRCURV_CONVEX;
      }

      /* otherwise, we can only say something if the child is linear */
      else if( childcurvature == SCIP_EXPRCURV_LINEAR )
         return (lbsin >= 0.0 && ubsin >= 0.0) ? SCIP_EXPRCURV_CONCAVE : SCIP_EXPRCURV_CONVEX;
   }

   return SCIP_EXPRCURV_UNKNOWN;
}

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrSin)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrSin(scip) );

   return SCIP_OKAY;
}

/** simplifies a sine expression
 *
 * Evaluates the sine value function when its child is a value expression.
 *
 * TODO: add further simplifications
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifySin)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /* check for value expression */
   if( SCIPisExprValue(scip, child) )
   {
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, sin(SCIPgetValueExprValue(child)), ownercreate,
            ownercreatedata) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_EXPRPARSE(parseSin)
{  /*lint --e{715}*/
   SCIP_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseExpr(scip, &childexpr, string, endstring, ownercreate, ownercreatedata) );
   assert(childexpr != NULL);

   /* create sine expression */
   SCIP_CALL( SCIPcreateExprSin(scip, expr, childexpr, ownercreate, ownercreatedata) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the sine expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_EXPREVAL(evalSin)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = sin(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffSin)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(expr != NULL);
   assert(childidx == 0);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(child)), "val") != 0);

   *val = cos(SCIPexprGetEvalValue(child));

   return SCIP_OKAY;
}

/** derivative evaluation callback
 *
 * Computes <gradient, children.dot>, that is, cos(child) dot(child).
 */
static
SCIP_DECL_EXPRFWDIFF(fwdiffSin)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(expr != NULL);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(child)), "val") != 0);
   assert(SCIPexprGetDot(child) != SCIP_INVALID); /*lint !e777*/

   *dot = cos(SCIPexprGetEvalValue(child)) * SCIPexprGetDot(child);

   return SCIP_OKAY;
}

/** expression backward forward derivative evaluation callback
 *
 * Computes partial/partial child ( <gradient, children.dot> ), that is, -sin(child) dot(child).
 */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffSin)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(expr != NULL);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID); /*lint !e777*/
   assert(childidx == 0);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(child)), "val") != 0);
   assert(SCIPexprGetDot(child) != SCIP_INVALID); /*lint !e777*/

   *bardot = -sin(SCIPexprGetEvalValue(child)) * SCIPexprGetDot(child);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalSin)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalSin(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_EXPRINITESTIMATES(initEstimatesSin)
{  /*lint --e{715}*/
   SCIP_Real childlb;
   SCIP_Real childub;

   childlb = bounds[0].inf;
   childub = bounds[0].sup;

   /* no need for cut if child is fixed */
   if( SCIPisRelEQ(scip, childlb, childub) )
      return SCIP_OKAY;

   /* compute cuts */
   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, ! overestimate, coefs, constant, nreturned) );

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_EXPRESTIMATE(estimateSin)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), SINEXPRHDLR_NAME) == 0);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);

   *success = computeEstimatorsTrig(scip, expr, coefs, constant, refpoint[0], localbounds[0].inf,
         localbounds[0].sup, ! overestimate);
   *islocal = TRUE;  /* TODO there are cases where cuts would be globally valid */

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropSin)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPintervalGetInf(bounds) >= -1.0);
   assert(SCIPintervalGetSup(bounds) <= 1.0);

   /* compute the new child interval */
   SCIP_CALL( computeRevPropIntervalSin(scip, bounds, childrenbounds[0], childrenbounds) );

   return SCIP_OKAY;
}

/** sine hash callback */
static
SCIP_DECL_EXPRHASH(hashSin)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   *hashkey = SINEXPRHDLR_HASHKEY;
   *hashkey ^= childrenhashes[0];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureSin)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_INTERVAL childinterval;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   SCIP_CALL( SCIPevalExprActivity(scip, child) );
   childinterval = SCIPexprGetActivity(child);

   /* TODO rewrite SCIPcomputeCurvatureSin so it provides the reverse operation */
   *success = TRUE;
   if( computeCurvatureSin(SCIP_EXPRCURV_CONVEX, childinterval.inf, childinterval.sup) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_CONVEX;
   else if( computeCurvatureSin(SCIP_EXPRCURV_CONCAVE, childinterval.inf, childinterval.sup) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_CONCAVE;
   if( computeCurvatureSin(SCIP_EXPRCURV_LINEAR, childinterval.inf, childinterval.sup) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_LINEAR;
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicitySin)
{  /*lint --e{715}*/
   SCIP_INTERVAL interval;
   SCIP_Real inf;
   SCIP_Real sup;
   int k;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   assert(SCIPexprGetChildren(expr)[0] != NULL);
   SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(expr)[0]) );
   interval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   *result = SCIP_MONOTONE_UNKNOWN;
   inf = SCIPintervalGetInf(interval);
   sup = SCIPintervalGetSup(interval);

   /* expression is not monotone because the interval is too large */
   if( sup - inf > M_PI )
      return SCIP_OKAY;

   /* compute k s.t. PI * (2k+1) / 2 <= interval.inf <= PI * (2k+3) / 2 */
   k = (int)floor(inf/M_PI - 0.5);
   assert(M_PI * (2.0*k + 1.0) / 2.0 <= inf);
   assert(M_PI * (2.0*k + 3.0) / 2.0 >= inf);

   /* check whether [inf,sup] are in containing in an interval for which the sine function is monotone */
   if( M_PI * (2.0*k + 3.0) / 2.0 <= sup )
      *result = ((k % 2 + 2) % 2) == 1 ? SCIP_MONOTONE_INC : SCIP_MONOTONE_DEC;

   return SCIP_OKAY;
}


/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrCos)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrCos(scip) );

   return SCIP_OKAY;
}

/** simplifies a cosine expression
 *
 * Evaluates the cosine value function when its child is a value expression.
 *
 * TODO: add further simplifications
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyCos)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /* check for value expression */
   if( SCIPisExprValue(scip, child) )
   {
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, cos(SCIPgetValueExprValue(child)), ownercreate,
               ownercreatedata) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_EXPRPARSE(parseCos)
{  /*lint --e{715}*/
   SCIP_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseExpr(scip, &childexpr, string, endstring, ownercreate, ownercreatedata) );
   assert(childexpr != NULL);

   /* create cosine expression */
   SCIP_CALL( SCIPcreateExprCos(scip, expr, childexpr, ownercreate, ownercreatedata) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the cosine expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_EXPREVAL(evalCos)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = cos(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffCos)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(expr != NULL);
   assert(childidx == 0);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(child)), "val") != 0);

   *val = -sin(SCIPexprGetEvalValue(child));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalCos(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_EXPRINITESTIMATES(initEstimatesCos)
{
   SCIP_Real childlb;
   SCIP_Real childub;

   childlb = bounds[0].inf;
   childub = bounds[0].sup;

   /* no need for cut if child is fixed */
   if( SCIPisRelEQ(scip, childlb, childub) )
      return SCIP_OKAY;

   /* compute cuts */
   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, ! overestimate, coefs, constant, nreturned) );

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_EXPRESTIMATE(estimateCos)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), COSEXPRHDLR_NAME) == 0);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);

   *success = computeEstimatorsTrig(scip, expr, coefs, constant, refpoint[0], localbounds[0].inf,
         localbounds[0].sup, ! overestimate);
   *islocal = TRUE;  /* TODO there are cases where cuts would be globally valid */

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL newbounds;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   /* bounds should have been intersected with activity, which is within [-1,1] */
   assert(SCIPintervalGetInf(bounds) >= -1.0);
   assert(SCIPintervalGetSup(bounds) <= 1.0);

   /* get the child interval */
   newbounds = childrenbounds[0];

   /* shift child interval to match sine */
   SCIPintervalAddScalar(SCIP_INTERVAL_INFINITY, &newbounds, newbounds, M_PI_2);  /* TODO use bounds on Pi/2 instead of approximation of Pi/2 */

   /* compute the new child interval */
   SCIP_CALL( computeRevPropIntervalSin(scip, bounds, newbounds, &newbounds) );

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newbounds) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* shift the new interval back */
   SCIPintervalAddScalar(SCIP_INTERVAL_INFINITY, &childrenbounds[0], newbounds, -M_PI_2);  /* TODO use bounds on Pi/2 instead of approximation of Pi/2 */

   return SCIP_OKAY;
}

/** cosine hash callback */
static
SCIP_DECL_EXPRHASH(hashCos)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   *hashkey = COSEXPRHDLR_HASHKEY;
   *hashkey ^= childrenhashes[0];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureCos)
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

   /* TODO rewrite SCIPcomputeCurvatureSin so it provides the reverse operation */
   *success = TRUE;
   if( computeCurvatureSin(SCIP_EXPRCURV_CONCAVE, childinterval.inf + M_PI_2, childinterval.sup + M_PI_2) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_CONCAVE;
   else if( computeCurvatureSin(SCIP_EXPRCURV_CONVEX, childinterval.inf + M_PI_2, childinterval.sup + M_PI_2) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_CONVEX;
   else if( computeCurvatureSin(SCIP_EXPRCURV_LINEAR, childinterval.inf + M_PI_2, childinterval.sup + M_PI_2) == exprcurvature )
      *childcurv = SCIP_EXPRCURV_LINEAR;
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL interval;
   SCIP_Real inf;
   SCIP_Real sup;
   int k;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   assert(SCIPexprGetChildren(expr)[0] != NULL);
   SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(expr)[0]) );
   interval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   *result = SCIP_MONOTONE_UNKNOWN;
   inf = SCIPintervalGetInf(interval);
   sup = SCIPintervalGetSup(interval);

   /* expression is not monotone because the interval is too large */
   if( sup - inf > M_PI )
      return SCIP_OKAY;

   /* compute k s.t. PI * k <= interval.inf <= PI * (k+1) */
   k = (int)floor(inf/M_PI);
   assert(M_PI * k <= inf);
   assert(M_PI * (k+1) >= inf);

   /* check whether [inf,sup] are contained in an interval for which the cosine function is monotone */
   if( sup <= M_PI * (k+1) )
      *result = ((k % 2 + 2) % 2) == 0 ? SCIP_MONOTONE_DEC : SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** creates the handler for sin expressions and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrSin(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   /* include expression handler */
   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, SINEXPRHDLR_NAME, SINEXPRHDLR_DESC, SINEXPRHDLR_PRECEDENCE, evalSin, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrSin, NULL);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifySin);
   SCIPexprhdlrSetParse(exprhdlr, parseSin);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalSin);
   SCIPexprhdlrSetEstimate(exprhdlr, initEstimatesSin, estimateSin);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropSin);
   SCIPexprhdlrSetHash(exprhdlr, hashSin);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffSin, fwdiffSin, bwfwdiffSin);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureSin);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicitySin);

   return SCIP_OKAY;
}

/** creates the handler for cos expressions and includes it SCIP */
SCIP_RETCODE SCIPincludeExprhdlrCos(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   /* include expression handler */
   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, COSEXPRHDLR_NAME, COSEXPRHDLR_DESC, COSEXPRHDLR_PRECEDENCE, evalCos, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrCos, NULL);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyCos);
   SCIPexprhdlrSetParse(exprhdlr, parseCos);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalCos);
   SCIPexprhdlrSetEstimate(exprhdlr, initEstimatesCos, estimateCos);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropCos);
   SCIPexprhdlrSetHash(exprhdlr, hashCos);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffCos, NULL, NULL);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureCos);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityCos);

   return SCIP_OKAY;
}

/** creates a sin expression */
SCIP_RETCODE SCIPcreateExprSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindExprhdlr(scip, SINEXPRHDLR_NAME) != NULL);

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPfindExprhdlr(scip, SINEXPRHDLR_NAME), NULL, 1, &child, ownercreate,
            ownercreatedata) );

   return SCIP_OKAY;
}


/** creates a cos expression */
SCIP_RETCODE SCIPcreateExprCos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindExprhdlr(scip, COSEXPRHDLR_NAME) != NULL);

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPfindExprhdlr(scip, COSEXPRHDLR_NAME), NULL, 1, &child, ownercreate,
            ownercreatedata) );

   return SCIP_OKAY;
}

/** indicates whether expression is of sine-type */  /*lint -e{715}*/
SCIP_Bool SCIPisExprSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), SINEXPRHDLR_NAME) == 0;
}

/** indicates whether expression is of cosine-type */  /*lint -e{715}*/
SCIP_Bool SCIPisExprCos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), COSEXPRHDLR_NAME) == 0;
}
