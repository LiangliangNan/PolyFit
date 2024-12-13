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

/**@file   nlhdlr_bilinear.c
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  bilinear nonlinear handler
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/nlhdlr_bilinear.h"
#include "scip/cons_nonlinear.h"
#include "scip/expr_product.h"
#include "scip/expr_var.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "bilinear"
#define NLHDLR_DESC               "bilinear handler for expressions"
#define NLHDLR_DETECTPRIORITY     -10  /**< it is important that the nlhdlr runs after the default nldhlr */
#define NLHDLR_ENFOPRIORITY       -10

#define MIN_INTERIORITY           0.01 /**< minimum interiority for a reference point for applying separation */
#define MIN_ABSBOUNDSIZE          0.1  /**< minimum size of variable bounds for applying separation */

/* properties of the bilinear nlhdlr statistics table */
#define TABLE_NAME_BILINEAR                 "nlhdlr_bilinear"
#define TABLE_DESC_BILINEAR                 "bilinear nlhdlr statistics table"
#define TABLE_POSITION_BILINEAR             14800                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_BILINEAR       SCIP_STAGE_INITSOLVE   /**< output of the statistics table is only printed from this stage onwards */


/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_NlhdlrExprData
{
   SCIP_Real             underineqs[6];      /**< inequalities for underestimation */
   int                   nunderineqs;        /**< total number of inequalities for underestimation */
   SCIP_Real             overineqs[6];       /**< inequalities for overestimation */
   int                   noverineqs;         /**< total number of inequalities for overestimation */
   SCIP_Longint          lastnodeid;         /**< id of the last node that has been used for separation */
   int                   nseparoundslastnode; /**< number of separation calls of the last node */
};

/** nonlinear handler data */
struct SCIP_NlhdlrData
{
   SCIP_EXPR**           exprs;             /**< expressions that have been detected by the nlhdlr */
   int                   nexprs;            /**< total number of expression that have been detected */
   int                   exprsize;          /**< size of exprs array */
   SCIP_HASHMAP*         exprmap;           /**< hashmap to store the position of each expression in the exprs array */

   /* parameter */
   SCIP_Bool             useinteval;        /**< whether to use the interval evaluation callback of the nlhdlr */
   SCIP_Bool             usereverseprop;    /**< whether to use the reverse propagation callback of the nlhdlr */
   int                   maxseparoundsroot; /**< maximum number of separation rounds in the root node */
   int                   maxseparounds;     /**< maximum number of separation rounds in a local node */
   int                   maxsepadepth;      /**< maximum depth to apply separation */
};

/*
 * Local methods
 */

/** helper function to compute the violation of an inequality of the form xcoef * x <= ycoef * y + constant for two
 *  corner points of the domain [lbx,ubx] x [lby,uby]
 */
static
void getIneqViol(
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             xcoef,              /**< x-coefficient */
   SCIP_Real             ycoef,              /**< y-coefficient */
   SCIP_Real             constant,           /**< constant */
   SCIP_Real*            viol1,              /**< buffer to store the violation of the first corner point */
   SCIP_Real*            viol2               /**< buffer to store the violation of the second corner point */
   )
{
   SCIP_Real norm;
   assert(viol1 != NULL);
   assert(viol2 != NULL);

   norm = SQRT(SQR(xcoef) + SQR(ycoef));

   /* inequality can be used for underestimating xy if and only if xcoef * ycoef > 0 */
   if( xcoef * ycoef >= 0 )
   {
      /* violation for top-left and bottom-right corner */
      *viol1 = MAX(0, (xcoef * SCIPvarGetLbLocal(x)  - ycoef * SCIPvarGetUbLocal(y) - constant) / norm); /*lint !e666*/ /*lint !e661*/
      *viol2 = MAX(0, (xcoef * SCIPvarGetUbLocal(x)  - ycoef * SCIPvarGetLbLocal(y) - constant) / norm); /*lint !e666*/ /*lint !e661*/
   }
   else
   {
      /* violation for top-right and bottom-left corner */
      *viol1 = MAX(0, (xcoef * SCIPvarGetUbLocal(x)  - ycoef * SCIPvarGetUbLocal(y) - constant) / norm); /*lint !e666*/ /*lint !e661*/
      *viol2 = MAX(0, (xcoef * SCIPvarGetLbLocal(x)  - ycoef * SCIPvarGetLbLocal(y) - constant) / norm); /*lint !e666*/ /*lint !e661*/
   }
}

/** auxiliary function to decide whether to use inequalities for a strong relaxation of bilinear terms or not */
static
SCIP_Bool useBilinIneqs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< x variable */
   SCIP_VAR*             y,                  /**< y variable */
   SCIP_Real             refx,               /**< reference point for x */
   SCIP_Real             refy                /**< reference point for y */
   )
{
   SCIP_Real lbx;
   SCIP_Real ubx;
   SCIP_Real lby;
   SCIP_Real uby;
   SCIP_Real interiorityx;
   SCIP_Real interiorityy;
   SCIP_Real interiority;

   assert(x != NULL);
   assert(y != NULL);
   assert(x != y);

   /* get variable bounds */
   lbx = SCIPvarGetLbLocal(x);
   ubx = SCIPvarGetUbLocal(x);
   lby = SCIPvarGetLbLocal(y);
   uby = SCIPvarGetUbLocal(y);

   /* compute interiority */
   interiorityx = MIN(refx-lbx, ubx-refx) / MAX(ubx-lbx, SCIPepsilon(scip)); /*lint !e666*/
   interiorityy = MIN(refy-lby, uby-refy) / MAX(uby-lby, SCIPepsilon(scip)); /*lint !e666*/
   interiority = 2.0*MIN(interiorityx, interiorityy);

   return ubx - lbx >= MIN_ABSBOUNDSIZE && uby - lby >= MIN_ABSBOUNDSIZE && interiority >= MIN_INTERIORITY;
}

/** helper function to update the best relaxation for a bilinear term when using valid linear inequalities */
static
void updateBilinearRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR* RESTRICT    x,                  /**< first variable */
   SCIP_VAR* RESTRICT    y,                  /**< second variable */
   SCIP_Real             bilincoef,          /**< coefficient of the bilinear term */
   SCIP_SIDETYPE         violside,           /**< side of quadratic constraint that is violated */
   SCIP_Real             refx,               /**< reference point for the x variable */
   SCIP_Real             refy,               /**< reference point for the y variable */
   SCIP_Real* RESTRICT   ineqs,              /**< coefficients of each linear inequality; stored as triple (xcoef,ycoef,constant) */
   int                   nineqs,             /**< total number of inequalities */
   SCIP_Real             mccormickval,       /**< value of the McCormick relaxation at the reference point */
   SCIP_Real* RESTRICT   bestcoefx,          /**< pointer to update the x coefficient */
   SCIP_Real* RESTRICT   bestcoefy,          /**< pointer to update the y coefficient */
   SCIP_Real* RESTRICT   bestconst,          /**< pointer to update the constant */
   SCIP_Real* RESTRICT   bestval,            /**< value of the best relaxation that have been found so far */
   SCIP_Bool*            success             /**< buffer to store whether we found a better relaxation */
   )
{
   SCIP_Real constshift[2] = {0.0, 0.0};
   SCIP_Real constant;
   SCIP_Real xcoef;
   SCIP_Real ycoef;
   SCIP_Real lbx;
   SCIP_Real ubx;
   SCIP_Real lby;
   SCIP_Real uby;
   SCIP_Bool update;
   SCIP_Bool overestimate;
   int i;

   assert(x != y);
   assert(!SCIPisZero(scip, bilincoef));
   assert(nineqs >= 0 && nineqs <= 2);
   assert(bestcoefx != NULL);
   assert(bestcoefy != NULL);
   assert(bestconst != NULL);
   assert(bestval != NULL);

   /* no inequalities available */
   if( nineqs == 0 )
      return;
   assert(ineqs != NULL);

   lbx = SCIPvarGetLbLocal(x);
   ubx = SCIPvarGetUbLocal(x);
   lby = SCIPvarGetLbLocal(y);
   uby = SCIPvarGetUbLocal(y);
   overestimate = (violside == SCIP_SIDETYPE_LEFT);

   /* check cases for which we can't compute a tighter relaxation */
   if( SCIPisFeasLE(scip, refx, lbx) || SCIPisFeasGE(scip, refx, ubx)
      || SCIPisFeasLE(scip, refy, lby) || SCIPisFeasGE(scip, refy, uby) )
      return;

   /* due to the feasibility tolerances of the LP and NLP solver, it might possible that the reference point is
    * violating the linear inequalities; to ensure that we compute a valid underestimate, we relax the linear
    * inequality by changing its constant part
    */
   for( i = 0; i < nineqs; ++i )
   {
      constshift[i] = MAX(0.0, ineqs[3*i] * refx - ineqs[3*i+1] * refy - ineqs[3*i+2]);
      SCIPdebugMsg(scip, "constant shift of inequality %d = %.16f\n", i, constshift[i]);
   }

   /* try to use both inequalities */
   if( nineqs == 2 )
   {
      SCIPcomputeBilinEnvelope2(scip, bilincoef, lbx, ubx, refx, lby, uby, refy, overestimate, ineqs[0], ineqs[1],
         ineqs[2] + constshift[0], ineqs[3], ineqs[4], ineqs[5] + constshift[1], &xcoef, &ycoef, &constant, &update);

      if( update )
      {
         SCIP_Real val = xcoef * refx + ycoef * refy + constant;
         SCIP_Real relimpr = 1.0 - (REALABS(val - bilincoef * refx * refy) + 1e-4) / (REALABS(*bestval - bilincoef * refx * refy) + 1e-4);
         SCIP_Real absimpr = REALABS(val - (*bestval));

         /* update relaxation if possible */
         if( relimpr > 0.05 && absimpr > 1e-3 && ((overestimate && SCIPisRelLT(scip, val, *bestval))
             || (!overestimate && SCIPisRelGT(scip, val, *bestval))) )
         {
            *bestcoefx = xcoef;
            *bestcoefy = ycoef;
            *bestconst = constant;
            *bestval = val;
            *success = TRUE;
         }
      }
   }

   /* use inequalities individually */
   for( i = 0; i < nineqs; ++i )
   {
      SCIPcomputeBilinEnvelope1(scip, bilincoef, lbx, ubx, refx, lby, uby, refy, overestimate, ineqs[3*i], ineqs[3*i+1],
         ineqs[3*i+2] + constshift[i], &xcoef, &ycoef, &constant, &update);

      if( update )
      {
         SCIP_Real val = xcoef * refx + ycoef * refy + constant;
         SCIP_Real relimpr = 1.0 - (REALABS(val - bilincoef * refx * refy) + 1e-4)
               / (REALABS(mccormickval - bilincoef * refx * refy) + 1e-4);
         SCIP_Real absimpr = REALABS(val - (*bestval));

         /* update relaxation if possible */
         if( relimpr > 0.05 && absimpr > 1e-3 && ((overestimate && SCIPisRelLT(scip, val, *bestval))
             || (!overestimate && SCIPisRelGT(scip, val, *bestval))) )
         {
            *bestcoefx = xcoef;
            *bestcoefy = ycoef;
            *bestconst = constant;
            *bestval = val;
            *success = TRUE;
         }
      }
   }
}

/** helper function to determine whether a given point satisfy given inequalities */
static
SCIP_Bool isPointFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x,                  /**< x-coordinate */
   SCIP_Real             y,                  /**< y-coordinate */
   SCIP_Real             lbx,                /**< lower bound of x */
   SCIP_Real             ubx,                /**< upper bound of x */
   SCIP_Real             lby,                /**< lower bound of y */
   SCIP_Real             uby,                /**< upper bound of y */
   SCIP_Real*            ineqs,              /**< inequalities of the form coefx x <= coefy y + constant */
   int                   nineqs              /**< total number of inequalities */
   )
{
   int i;

   assert(ineqs != NULL);
   assert(nineqs > 0);

   /* check whether point satisfies the bounds */
   if( SCIPisLT(scip, x, lbx) || SCIPisGT(scip, x, ubx)
      || SCIPisLT(scip, y, lby) || SCIPisGT(scip, y, uby) )
      return FALSE;

   /* check whether point satisfy the linear inequalities */
   for( i = 0; i < nineqs; ++i )
   {
      SCIP_Real coefx = ineqs[3*i];
      SCIP_Real coefy = ineqs[3*i+1];
      SCIP_Real constant = ineqs[3*i+2];

      /* TODO check with an absolute comparison? */
      if( SCIPisGT(scip, coefx*x - coefy*y - constant, 0.0) )
         return FALSE;
   }

   return TRUE;
}

/** helper function for computing all vertices of the polytope described by the linear inequalities and the local
 *  extrema of the bilinear term along each inequality
 *
 * @note there are at most 22 points where the min/max can be achieved (given that there are at most 4 inequalities)
 *          - corners of [lbx,ubx]x[lby,uby] (4)
 *          - two intersection points for each inequality with the box (8)
 *          - global maximum / minimum on each inequality (4)
 *          - intersection between two inequalities (6)
 */
static
void getFeasiblePointsBilinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler, if levelset == TRUE, otherwise can be NULL */
   SCIP_EXPR*            expr,               /**< product expression */
   SCIP_INTERVAL         exprbounds,         /**< bounds on product expression, only used if levelset == TRUE */
   SCIP_Real*            underineqs,         /**< inequalities for underestimation */
   int                   nunderineqs,        /**< total number of inequalities for underestimation */
   SCIP_Real*            overineqs,          /**< inequalities for overestimation */
   int                   noverineqs,         /**< total number of inequalities for overestimation */
   SCIP_Bool             levelset,           /**< should the level set be considered? */
   SCIP_Real*            xs,                 /**< array to store x-coordinates of computed points */
   SCIP_Real*            ys,                 /**< array to store y-coordinates of computed points */
   int*                  npoints             /**< buffer to store the total number of computed points */
   )
{
   SCIP_EXPR* child1;
   SCIP_EXPR* child2;
   SCIP_Real ineqs[12];
   SCIP_INTERVAL boundsx;
   SCIP_INTERVAL boundsy;
   SCIP_Real lbx;
   SCIP_Real ubx;
   SCIP_Real lby;
   SCIP_Real uby;
   int nineqs = 0;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL || !levelset);
   assert(expr != NULL);
   assert(xs != NULL);
   assert(ys != NULL);
   assert(SCIPexprGetNChildren(expr) == 2);
   assert(noverineqs + nunderineqs > 0);
   assert(noverineqs + nunderineqs <= 4);

   *npoints = 0;

   /* collect inequalities */
   for( i = 0; i < noverineqs; ++i )
   {
      SCIPdebugMsg(scip, "over-inequality  %d: %g*x <= %g*y + %g\n", i, overineqs[3*i], overineqs[3*i+1], overineqs[3*i+2]);
      ineqs[3*nineqs] = overineqs[3*i];
      ineqs[3*nineqs+1] = overineqs[3*i+1];
      ineqs[3*nineqs+2] = overineqs[3*i+2];
      ++nineqs;
   }
   for( i = 0; i < nunderineqs; ++i )
   {
      SCIPdebugMsg(scip, "under-inequality %d: %g*x <= %g*y + %g 0\n", i, underineqs[3*i], underineqs[3*i+1], underineqs[3*i+2]);
      ineqs[3*nineqs] = underineqs[3*i];
      ineqs[3*nineqs+1] = underineqs[3*i+1];
      ineqs[3*nineqs+2] = underineqs[3*i+2];
      ++nineqs;
   }
   assert(nineqs == noverineqs + nunderineqs);

   /* collect children */
   child1 = SCIPexprGetChildren(expr)[0];
   child2 = SCIPexprGetChildren(expr)[1];
   assert(child1 != NULL && child2 != NULL);
   assert(child1 != child2);

   /* collect bounds of children */
   if( !levelset )
   {
      /* if called from inteval, then use activity */
      boundsx = SCIPexprGetActivity(child1);
      boundsy = SCIPexprGetActivity(child2);
   }
   else
   {
      /* if called from reverseprop, then use bounds */
      boundsx = SCIPgetExprBoundsNonlinear(scip, child1);
      boundsy = SCIPgetExprBoundsNonlinear(scip, child2);

      /* if children bounds are empty, then returning with *npoints==0 is the way to go */
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, boundsx)
          || SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, boundsy) )
         return;
   }
   lbx = boundsx.inf;
   ubx = boundsx.sup;
   lby = boundsy.inf;
   uby = boundsy.sup;
   SCIPdebugMsg(scip, "x = [%g,%g], y=[%g,%g]\n", lbx, ubx, lby, uby);

   /* corner points that satisfy all inequalities */
   for( i = 0; i < 4; ++i )
   {
      SCIP_Real cx = i < 2 ? lbx : ubx;
      SCIP_Real cy = (i % 2) == 0 ? lby : uby;

      SCIPdebugMsg(scip, "corner point (%g,%g) feasible? %u\n", cx, cy, isPointFeasible(scip, cx, cy, lbx, ubx, lby, uby, ineqs, nineqs));

      if( isPointFeasible(scip, cx, cy, lbx, ubx, lby, uby, ineqs, nineqs) )
      {
         xs[*npoints] = cx;
         ys[*npoints] = cy;
         ++(*npoints);
      }
   }

   /* intersection point of inequalities with [lbx,ubx] x [lby,uby] and extremum of xy on each inequality */
   for( i = 0; i < nineqs; ++i )
   {
      SCIP_Real coefx = ineqs[3*i];
      SCIP_Real coefy = ineqs[3*i+1];
      SCIP_Real constant = ineqs[3*i+2];
      SCIP_Real px[5] = {lbx, ubx, (coefy*lby + constant)/coefx, (coefy*uby + constant)/coefx, 0.0};
      SCIP_Real py[5] = {(coefx*lbx - constant)/coefy, (coefx*ubx - constant)/coefy, lby, uby, 0.0};
      int j;

      /* the last entry corresponds to the extremum of xy on the line */
      py[4] = (-constant) / (2.0 * coefy);
      px[4] = constant / (2.0 * coefx);

      for( j = 0; j < 5; ++j )
      {
         SCIPdebugMsg(scip, "intersection point (%g,%g) feasible? %u\n", px[j], py[j], isPointFeasible(scip, px[j], py[j], lbx, ubx, lby, uby, ineqs, nineqs));
         if( isPointFeasible(scip, px[j], py[j], lbx, ubx, lby, uby, ineqs, nineqs) )
         {
            xs[*npoints] = px[j];
            ys[*npoints] = py[j];
            ++(*npoints);
         }
      }
   }

   /* intersection point between two inequalities */
   for( i = 0; i < nineqs - 1; ++i )
   {
      SCIP_Real coefx1 = ineqs[3*i];
      SCIP_Real coefy1 = ineqs[3*i+1];
      SCIP_Real constant1 = ineqs[3*i+2];
      int j;

      for( j = i + 1; j < nineqs; ++j )
      {
         SCIP_Real coefx2 = ineqs[3*j];
         SCIP_Real coefy2 = ineqs[3*j+1];
         SCIP_Real constant2 = ineqs[3*j+2];
         SCIP_Real px;
         SCIP_Real py;

         /* no intersection point -> skip */
         if( SCIPisZero(scip, coefx2*coefy1 - coefx1 * coefy2) )
            continue;

         py = (constant2 * coefx1 - constant1 * coefx2)/ (coefx2 * coefy1 - coefx1 * coefy2);
         px = (coefy1 * py + constant1) / coefx1;
         assert(SCIPisRelEQ(scip, px, (coefy2 * py + constant2) / coefx2));

         if( isPointFeasible(scip, px, py, lbx, ubx, lby, uby, ineqs, nineqs) )
         {
            xs[*npoints] = px;
            ys[*npoints] = py;
            ++(*npoints);
         }
      }
   }

   assert(*npoints <= 22);

   /* consider the intersection of the level set with
    *
    * 1. the boundary of the box
    * 2. the linear inequalities
    *
    * this adds at most for 4 (level set curves) * 4 (inequalities) * 2 (intersection points) for all linear
    * inequalities and 4 (level set curves) * 2 (intersection points) with the boundary of the box
    */
   if( !levelset )
      return;

   /* compute intersection of level sets with the boundary */
   for( i = 0; i < 2; ++i )
   {
      SCIP_Real vals[4] = {lbx, ubx, lby, uby};
      SCIP_Real val;
      int k;

      /* fix auxiliary variable to its lower or upper bound and consider the coefficient of the product */
      val = (i == 0) ? exprbounds.inf : exprbounds.sup;
      val /= SCIPgetCoefExprProduct(expr);

      for( k = 0; k < 4; ++k )
      {
         if( !SCIPisZero(scip, vals[k]) )
         {
            SCIP_Real res = val / vals[k];

            assert(SCIPisRelGE(scip, SCIPgetCoefExprProduct(expr)*res*vals[k], exprbounds.inf));
            assert(SCIPisRelLE(scip, SCIPgetCoefExprProduct(expr)*res*vals[k], exprbounds.sup));

            /* fix x to lbx or ubx */
            if( k < 2 && isPointFeasible(scip, vals[k], res, lbx, ubx, lby, uby, ineqs, nineqs) )
            {
               xs[*npoints] = vals[k];
               ys[*npoints] = res;
               ++(*npoints);
            }
            /* fix y to lby or uby */
            else if( k >= 2 &&  isPointFeasible(scip, res, vals[k], lbx, ubx, lby, uby, ineqs, nineqs) )
            {
               xs[*npoints] = res;
               ys[*npoints] = vals[k];
               ++(*npoints);
            }
         }
      }
   }

   /* compute intersection points of level sets with the linear inequalities */
   for( i = 0; i < nineqs; ++i )
   {
      SCIP_INTERVAL result;
      SCIP_Real coefx = ineqs[3*i];
      SCIP_Real coefy = ineqs[3*i+1];
      SCIP_Real constant = ineqs[3*i+2];
      SCIP_INTERVAL sqrcoef;
      SCIP_INTERVAL lincoef;
      SCIP_Real px;
      SCIP_Real py;
      int k;

      /* solve system of coefx x = coefy y + constant and X = xy which is the same as computing the solutions of
       *
       *    (coefy / coefx) y^2 + (constant / coefx) y = inf(X) or sup(X)
       */
      SCIPintervalSet(&sqrcoef, coefy / coefx);
      SCIPintervalSet(&lincoef, constant / coefx);

      for( k = 0; k < 2; ++k )
      {
         SCIP_INTERVAL rhs;
         SCIP_INTERVAL ybnds;

         /* set right-hand side */
         if( k == 0 )
            SCIPintervalSet(&rhs, exprbounds.inf);
         else
            SCIPintervalSet(&rhs, exprbounds.sup);

         SCIPintervalSetBounds(&ybnds, lby, uby);
         SCIPintervalSolveUnivariateQuadExpression(SCIP_INTERVAL_INFINITY, &result, sqrcoef, lincoef, rhs, ybnds);

         /* interval is empty -> no solution available */
         if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, result) )
            continue;

         /* compute and check point */
         py = SCIPintervalGetInf(result);
         px = (coefy * py + constant) / coefx;

         if( isPointFeasible(scip, px, py, lbx, ubx, lby, uby, ineqs, nineqs) )
         {
            xs[*npoints] = px;
            ys[*npoints] = py;
            ++(*npoints);
         }

         /* check for a second solution */
         if( SCIPintervalGetInf(result) != SCIPintervalGetSup(result) ) /*lint !e777*/
         {
            py = SCIPintervalGetSup(result);
            px = (coefy * py + constant) / coefx;

            if( isPointFeasible(scip, px, py, lbx, ubx, lby, uby, ineqs, nineqs) )
            {
               xs[*npoints] = px;
               ys[*npoints] = py;
               ++(*npoints);
            }
         }
      }
   }

   assert(*npoints <= 62);
}

/** computes interval for a bilinear term when using at least one inequality */
static
SCIP_INTERVAL intevalBilinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< product expression */
   SCIP_Real*            underineqs,         /**< inequalities for underestimation */
   int                   nunderineqs,        /**< total number of inequalities for underestimation */
   SCIP_Real*            overineqs,          /**< inequalities for overestimation */
   int                   noverineqs          /**< total number of inequalities for overestimation */
   )
{
   SCIP_INTERVAL interval = {0., 0.};
   SCIP_Real xs[22];
   SCIP_Real ys[22];
   SCIP_Real inf;
   SCIP_Real sup;
   int npoints;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 2);
   assert(noverineqs + nunderineqs <= 4);

   /* no inequalities available -> skip computation */
   if( noverineqs == 0 && nunderineqs == 0 )
   {
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &interval);
      return interval;
   }

   /* x or y has empty interval -> empty */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(SCIPexprGetChildren(expr)[0])) ||
       SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(SCIPexprGetChildren(expr)[1])) )
   {
      SCIPintervalSetEmpty(&interval);
      return interval;
   }

   /* compute all feasible points (since we use levelset == FALSE, the value of interval doesn't matter) */
   getFeasiblePointsBilinear(scip, NULL, expr, interval, underineqs, nunderineqs, overineqs,
         noverineqs, FALSE, xs, ys, &npoints);

   /* no feasible point left -> return an empty interval */
   if( npoints == 0 )
   {
      SCIPintervalSetEmpty(&interval);
      return interval;
   }

   /* compute the minimum and maximum over all computed points */
   inf = xs[0] * ys[0];
   sup = inf;
   SCIPdebugMsg(scip, "point 0: (%g,%g) -> inf = sup = %g\n", xs[0], ys[0], inf);
   for( i = 1; i < npoints; ++i )
   {
      inf = MIN(inf, xs[i] * ys[i]);
      sup = MAX(sup, xs[i] * ys[i]);
      SCIPdebugMsg(scip, "point %d: (%g,%g) -> inf = %g, sup = %g\n", i, xs[i], ys[i], inf, sup);
   }
   assert(inf <= sup);

   /* adjust infinite values */
   inf = MAX(inf, -SCIP_INTERVAL_INFINITY);
   sup = MIN(sup, SCIP_INTERVAL_INFINITY);

   /* multiply resulting interval with coefficient of the product expression */
   SCIPintervalSetBounds(&interval, inf, sup);
   if( SCIPgetCoefExprProduct(expr) != 1.0 )
      SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &interval, interval, SCIPgetCoefExprProduct(expr));

   return interval;
}

/** uses inequalities for bilinear terms to get stronger bounds during reverse propagation */
static
void reversePropBilinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_EXPR*            expr,               /**< product expression */
   SCIP_INTERVAL         exprbounds,         /**< bounds on product expression */
   SCIP_Real*            underineqs,         /**< inequalities for underestimation */
   int                   nunderineqs,        /**< total number of inequalities for underestimation */
   SCIP_Real*            overineqs,          /**< inequalities for overestimation */
   int                   noverineqs,         /**< total number of inequalities for overestimation */
   SCIP_INTERVAL*        intervalx,          /**< buffer to store the new interval for x */
   SCIP_INTERVAL*        intervaly           /**< buffer to store the new interval for y */
   )
{
   SCIP_Real xs[62];
   SCIP_Real ys[62];
   SCIP_Real exprinf;
   SCIP_Real exprsup;
   SCIP_Bool first = TRUE;
   int npoints;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(intervalx != NULL);
   assert(intervaly != NULL);
   assert(SCIPexprGetNChildren(expr) == 2);

   assert(noverineqs + nunderineqs > 0);

   /* set intervals to be empty */
   SCIPintervalSetEmpty(intervalx);
   SCIPintervalSetEmpty(intervaly);

   /* compute feasible points */
   getFeasiblePointsBilinear(scip, conshdlr, expr, exprbounds, underineqs, nunderineqs, overineqs,
         noverineqs, TRUE, xs, ys, &npoints);

   /* no feasible points left -> problem is infeasible */
   if( npoints == 0 )
      return;

   /* get bounds of the product expression */
   exprinf = exprbounds.inf;
   exprsup = exprbounds.sup;

   /* update intervals with the computed points */
   for( i = 0; i < npoints; ++i )
   {
      SCIP_Real val = SCIPgetCoefExprProduct(expr) * xs[i] * ys[i];

#ifndef NDEBUG
      {
         SCIP_Real lbx = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]).inf;
         SCIP_Real ubx = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]).sup;
         SCIP_Real lby = SCIPexprGetActivity(SCIPexprGetChildren(expr)[1]).inf;
         SCIP_Real uby = SCIPexprGetActivity(SCIPexprGetChildren(expr)[1]).sup;

         assert(nunderineqs == 0 || isPointFeasible(scip, xs[i], ys[i], lbx, ubx, lby, uby, underineqs, nunderineqs));
         assert(noverineqs == 0 || isPointFeasible(scip, xs[i], ys[i], lbx, ubx, lby, uby, overineqs, noverineqs));
      }
#endif

      /* only accept points for which the value of x*y is in the interval of the product expression
       *
       * NOTE: in order to consider all relevant points, we are a bit conservative here and relax the interval of
       *      the expression by SCIPfeastol()
       */
      if( SCIPisRelGE(scip, val, exprinf - SCIPfeastol(scip)) && SCIPisRelLE(scip, val, exprsup + SCIPfeastol(scip)) )
      {
         if( first )
         {
            SCIPintervalSet(intervalx, xs[i]);
            SCIPintervalSet(intervaly, ys[i]);
            first = FALSE;
         }
         else
         {
            (*intervalx).inf = MIN((*intervalx).inf, xs[i]);
            (*intervalx).sup = MAX((*intervalx).sup, xs[i]);
            (*intervaly).inf = MIN((*intervaly).inf, ys[i]);
            (*intervaly).sup = MAX((*intervaly).sup, ys[i]);
         }

         SCIPdebugMsg(scip, "consider points (%g,%g)=%g for reverse propagation\n", xs[i], ys[i], val);
      }
   }
}

/** helper function to compute the convex envelope of a bilinear term when two linear inequalities are given; we
 *  use the same notation and formulas as in Locatelli 2016
 */
static
void computeBilinEnvelope2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x,                  /**< reference point for x */
   SCIP_Real             y,                  /**< reference point for y */
   SCIP_Real             mi,                 /**< coefficient of x in the first linear inequality */
   SCIP_Real             qi,                 /**< constant in the first linear inequality */
   SCIP_Real             mj,                 /**< coefficient of x in the second linear inequality */
   SCIP_Real             qj,                 /**< constant in the second linear inequality */
   SCIP_Real* RESTRICT   xi,                 /**< buffer to store x coordinate of the first point */
   SCIP_Real* RESTRICT   yi,                 /**< buffer to store y coordinate of the first point */
   SCIP_Real* RESTRICT   xj,                 /**< buffer to store x coordinate of the second point */
   SCIP_Real* RESTRICT   yj,                 /**< buffer to store y coordinate of the second point */
   SCIP_Real* RESTRICT   xcoef,              /**< buffer to store the x coefficient of the envelope */
   SCIP_Real* RESTRICT   ycoef,              /**< buffer to store the y coefficient of the envelope */
   SCIP_Real* RESTRICT   constant            /**< buffer to store the constant of the envelope */
   )
{
   SCIP_Real QUAD(xiq);
   SCIP_Real QUAD(yiq);
   SCIP_Real QUAD(xjq);
   SCIP_Real QUAD(yjq);
   SCIP_Real QUAD(xcoefq);
   SCIP_Real QUAD(ycoefq);
   SCIP_Real QUAD(constantq);
   SCIP_Real QUAD(tmpq);

   assert(xi != NULL);
   assert(yi != NULL);
   assert(xj != NULL);
   assert(yj != NULL);
   assert(xcoef != NULL);
   assert(ycoef != NULL);
   assert(constant != NULL);

   if( SCIPisEQ(scip, mi, mj) )
   {
      /* xi = (x + mi * y - qi) / (2.0*mi) */
      SCIPquadprecProdDD(xiq, mi, y);
      SCIPquadprecSumQD(xiq, xiq, x);
      SCIPquadprecSumQD(xiq, xiq, -qi);
      SCIPquadprecDivQD(xiq, xiq, 2.0 * mi);
      assert(EPSEQ((x + mi * y - qi) / (2.0*mi), QUAD_TO_DBL(xiq), 1e-3));

      /* yi = mi*(*xi) + qi */
      SCIPquadprecProdQD(yiq, xiq, mi);
      SCIPquadprecSumQD(yiq, yiq, qi);
      assert(EPSEQ(mi*QUAD_TO_DBL(xiq) + qi, QUAD_TO_DBL(yiq), 1e-3));

      /* xj = (*xi) + (qi - qj)/ (2.0*mi) */
      SCIPquadprecSumDD(xjq, qi, -qj);
      SCIPquadprecDivQD(xjq, xjq, 2.0 * mi);
      SCIPquadprecSumQQ(xjq, xjq, xiq);
      assert(EPSEQ(QUAD_TO_DBL(xiq) + (qi - qj)/ (2.0*mi), QUAD_TO_DBL(xjq), 1e-3));

      /* yj = mj * (*xj) + qj */
      SCIPquadprecProdQD(yjq, xjq, mj);
      SCIPquadprecSumQD(yjq, yjq, qj);
      assert(EPSEQ(mj * QUAD_TO_DBL(xjq) + qj, QUAD_TO_DBL(yjq), 1e-3));

      /* ycoef = (*xi) + (qi - qj) / (4.0*mi) note that this is wrong in Locatelli 2016 */
      SCIPquadprecSumDD(ycoefq, qi, -qj);
      SCIPquadprecDivQD(ycoefq, ycoefq, 4.0 * mi);
      SCIPquadprecSumQQ(ycoefq, ycoefq, xiq);
      assert(EPSEQ(QUAD_TO_DBL(xiq) + (qi - qj) / (4.0*mi), QUAD_TO_DBL(ycoefq), 1e-3));

      /* xcoef = 2.0*mi*(*xi) - mi * (*ycoef) + qi */
      SCIPquadprecProdQD(xcoefq, xiq, 2.0 * mi);
      SCIPquadprecProdQD(tmpq, ycoefq, -mi);
      SCIPquadprecSumQQ(xcoefq, xcoefq, tmpq);
      SCIPquadprecSumQD(xcoefq, xcoefq, qi);
      assert(EPSEQ(2.0*mi*QUAD_TO_DBL(xiq) - mi * QUAD_TO_DBL(ycoefq) + qi, QUAD_TO_DBL(xcoefq), 1e-3));

      /* constant = -mj*SQR(*xj) - (*ycoef) * qj */
      SCIPquadprecSquareQ(constantq, xjq);
      SCIPquadprecProdQD(constantq, constantq, -mj);
      SCIPquadprecProdQD(tmpq, ycoefq, -qj);
      SCIPquadprecSumQQ(constantq, constantq, tmpq);
      /* assert(EPSEQ(-mj*SQR(QUAD_TO_DBL(xjq)) - QUAD_TO_DBL(ycoefq) * qj, QUAD_TO_DBL(constantq), 1e-3)); */

      *xi = QUAD_TO_DBL(xiq);
      *yi = QUAD_TO_DBL(yiq);
      *xj = QUAD_TO_DBL(xjq);
      *yj = QUAD_TO_DBL(yjq);
      *ycoef = QUAD_TO_DBL(ycoefq);
      *xcoef = QUAD_TO_DBL(xcoefq);
      *constant = QUAD_TO_DBL(constantq);
   }
   else if( mi > 0.0 )
   {
      assert(mj > 0.0);

      /* xi = (y + SQRT(mi*mj)*x - qi) / (REALABS(mi) + SQRT(mi*mj)) */
      SCIPquadprecProdDD(xiq, mi, mj);
      SCIPquadprecSqrtQ(xiq, xiq);
      SCIPquadprecProdQD(xiq, xiq, x);
      SCIPquadprecSumQD(xiq, xiq, y);
      SCIPquadprecSumQD(xiq, xiq, -qi); /* (y + SQRT(mi*mj)*x - qi) */
      SCIPquadprecProdDD(tmpq, mi, mj);
      SCIPquadprecSqrtQ(tmpq, tmpq);
      SCIPquadprecSumQD(tmpq, tmpq, REALABS(mi)); /* REALABS(mi) + SQRT(mi*mj) */
      SCIPquadprecDivQQ(xiq, xiq, tmpq);
      assert(EPSEQ((y + SQRT(mi*mj)*x - qi) / (REALABS(mi) + SQRT(mi*mj)), QUAD_TO_DBL(xiq), 1e-3));

      /* yi = mi*(*xi) + qi */
      SCIPquadprecProdQD(yiq, xiq, mi);
      SCIPquadprecSumQD(yiq, yiq, qi);
      assert(EPSEQ(mi*(QUAD_TO_DBL(xiq)) + qi, QUAD_TO_DBL(yiq), 1e-3));

      /* xj = (y + SQRT(mi*mj)*x - qj) / (REALABS(mj) + SQRT(mi*mj)) */
      SCIPquadprecProdDD(xjq, mi, mj);
      SCIPquadprecSqrtQ(xjq, xjq);
      SCIPquadprecProdQD(xjq, xjq, x);
      SCIPquadprecSumQD(xjq, xjq, y);
      SCIPquadprecSumQD(xjq, xjq, -qj); /* (y + SQRT(mi*mj)*x - qj) */
      SCIPquadprecProdDD(tmpq, mi, mj);
      SCIPquadprecSqrtQ(tmpq, tmpq);
      SCIPquadprecSumQD(tmpq, tmpq, REALABS(mj)); /* REALABS(mj) + SQRT(mi*mj) */
      SCIPquadprecDivQQ(xjq, xjq, tmpq);
      assert(EPSEQ((y + SQRT(mi*mj)*x - qj) / (REALABS(mj) + SQRT(mi*mj)), QUAD_TO_DBL(xjq), 1e-3));

      /* yj = mj*(*xj) + qj */
      SCIPquadprecProdQD(yjq, xjq, mj);
      SCIPquadprecSumQD(yjq, yjq, qj);
      assert(EPSEQ(mj*QUAD_TO_DBL(xjq) + qj, QUAD_TO_DBL(yjq), 1e-3));

      /* ycoef = (2.0*mj*(*xj) + qj - 2.0*mi*(*xi) - qi) / (mj - mi) */
      SCIPquadprecProdQD(ycoefq, xjq, 2.0 * mj);
      SCIPquadprecSumQD(ycoefq, ycoefq, qj);
      SCIPquadprecProdQD(tmpq, xiq, -2.0 * mi);
      SCIPquadprecSumQQ(ycoefq, ycoefq, tmpq);
      SCIPquadprecSumQD(ycoefq, ycoefq, -qi);
      SCIPquadprecSumDD(tmpq, mj, -mi);
      SCIPquadprecDivQQ(ycoefq, ycoefq, tmpq);
      assert(EPSEQ((2.0*mj*QUAD_TO_DBL(xjq) + qj - 2.0*mi*QUAD_TO_DBL(xiq) - qi) / (mj - mi), QUAD_TO_DBL(ycoefq), 1e-3));

      /* xcoef = 2.0*mj*(*xj) + qj - mj*(*ycoef) */
      SCIPquadprecProdQD(xcoefq, xjq, 2.0 * mj);
      SCIPquadprecSumQD(xcoefq, xcoefq, qj);
      SCIPquadprecProdQD(tmpq, ycoefq, -mj);
      SCIPquadprecSumQQ(xcoefq, xcoefq, tmpq);
      assert(EPSEQ(2.0*mj*QUAD_TO_DBL(xjq) + qj - mj*QUAD_TO_DBL(ycoefq), QUAD_TO_DBL(xcoefq), 1e-3));

      /* constant = -mj*SQR(*xj) - (*ycoef) * qj */
      SCIPquadprecSquareQ(constantq, xjq);
      SCIPquadprecProdQD(constantq, constantq, -mj);
      SCIPquadprecProdQD(tmpq, ycoefq, -qj);
      SCIPquadprecSumQQ(constantq, constantq, tmpq);
      /* assert(EPSEQ(-mj*SQR(QUAD_TO_DBL(xjq)) - QUAD_TO_DBL(ycoefq) * qj, QUAD_TO_DBL(constantq), 1e-3)); */

      *xi = QUAD_TO_DBL(xiq);
      *yi = QUAD_TO_DBL(yiq);
      *xj = QUAD_TO_DBL(xjq);
      *yj = QUAD_TO_DBL(yjq);
      *ycoef = QUAD_TO_DBL(ycoefq);
      *xcoef = QUAD_TO_DBL(xcoefq);
      *constant = QUAD_TO_DBL(constantq);
   }
   else
   {
      assert(mi < 0.0 && mj < 0.0);

      /* apply variable transformation x = -x in case for overestimation */
      computeBilinEnvelope2(scip, -x, y, -mi, qi, -mj, qj, xi, yi, xj, yj, xcoef, ycoef, constant);

      /* revert transformation; multiply cut by -1 and change -x by x */
      *xi = -(*xi);
      *xj = -(*xj);
      *ycoef = -(*ycoef);
      *constant = -(*constant);
   }
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputBilinear)
{ /*lint --e{715}*/
   SCIP_NLHDLR* nlhdlr;
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_HASHMAP* hashmap;
   SCIP_EXPRITER* it;
   int resfound = 0;
   int restotal = 0;
   int c;

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   assert(conshdlr != NULL);
   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, NLHDLR_NAME);
   assert(nlhdlr != NULL);
   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPhashmapCreate(&hashmap, SCIPblkmem(scip), nlhdlrdata->nexprs) );
   SCIP_CALL( SCIPcreateExpriter(scip, &it) );

   for( c = 0; c < nlhdlrdata->nexprs; ++c )
   {
      assert(!SCIPhashmapExists(hashmap, nlhdlrdata->exprs[c]));
      SCIP_CALL( SCIPhashmapInsertInt(hashmap, nlhdlrdata->exprs[c], 0) );
   }

   /* count in how many constraints each expression is contained */
   for( c = 0; c < SCIPconshdlrGetNConss(conshdlr); ++c )
   {
      SCIP_CONS* cons = SCIPconshdlrGetConss(conshdlr)[c];
      SCIP_EXPR* expr;

      SCIP_CALL( SCIPexpriterInit(it, SCIPgetExprNonlinear(cons), SCIP_EXPRITER_DFS, FALSE) );

      for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*//*lint !e440*/
      {
         if( SCIPhashmapExists(hashmap, expr) )
         {
            int nuses = SCIPhashmapGetImageInt(hashmap, expr);
            SCIP_CALL( SCIPhashmapSetImageInt(hashmap, expr, nuses + 1) );
         }
      }
   }

   /* compute success ratio */
   for( c = 0; c < nlhdlrdata->nexprs; ++c )
   {
      SCIP_NLHDLREXPRDATA* nlhdlrexprdata;
      int nuses;

      nuses = SCIPhashmapGetImageInt(hashmap, nlhdlrdata->exprs[c]);
      assert(nuses > 0);

      nlhdlrexprdata = SCIPgetNlhdlrExprDataNonlinear(nlhdlr, nlhdlrdata->exprs[c]);
      assert(nlhdlrexprdata != NULL);

      if( nlhdlrexprdata->nunderineqs > 0 || nlhdlrexprdata->noverineqs > 0 )
         resfound += nuses;
      restotal += nuses;
   }

   /* print statistics */
   SCIPinfoMessage(scip, file, "Bilinear Nlhdlr    : %10s %10s\n", "#found", "#total");
   SCIPinfoMessage(scip, file, "  %-17s:", "");
   SCIPinfoMessage(scip, file, " %10d", resfound);
   SCIPinfoMessage(scip, file, " %10d", restotal);
   SCIPinfoMessage(scip, file, "\n");

   /* free memory */
   SCIPfreeExpriter(&it);
   SCIPhashmapFree(&hashmap);

   return SCIP_OKAY;
}


/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrBilinear)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrBilinear(targetscip) );

   return SCIP_OKAY;
}

/** callback to free data of handler */
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataBilinear)
{ /*lint --e{715}*/
   assert(nlhdlrdata != NULL);
   assert((*nlhdlrdata)->nexprs == 0);

   if( (*nlhdlrdata)->exprmap != NULL )
   {
      assert(SCIPhashmapGetNElements((*nlhdlrdata)->exprmap) == 0);
      SCIPhashmapFree(&(*nlhdlrdata)->exprmap);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrdata)->exprs, (*nlhdlrdata)->exprsize);
   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataBilinear)
{  /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   int pos;

   assert(expr != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);
   assert(nlhdlrdata->nexprs > 0);
   assert(nlhdlrdata->exprs != NULL);
   assert(nlhdlrdata->exprmap != NULL);
   assert(SCIPhashmapExists(nlhdlrdata->exprmap, (void*)expr));

   pos = SCIPhashmapGetImageInt(nlhdlrdata->exprmap, (void*)expr);
   assert(pos >= 0 && pos < nlhdlrdata->nexprs);
   assert(nlhdlrdata->exprs[pos] == expr);

   /* move the last expression to the free position */
   if( nlhdlrdata->nexprs > 0 && pos != nlhdlrdata->nexprs - 1 )
   {
      SCIP_EXPR* lastexpr = nlhdlrdata->exprs[nlhdlrdata->nexprs - 1];
      assert(expr != lastexpr);
      assert(SCIPhashmapExists(nlhdlrdata->exprmap, (void*)lastexpr));

      nlhdlrdata->exprs[pos] = lastexpr;
      nlhdlrdata->exprs[nlhdlrdata->nexprs - 1] = NULL;
      SCIP_CALL( SCIPhashmapSetImageInt(nlhdlrdata->exprmap, (void*)lastexpr, pos) );
   }

   /* remove expression from the nonlinear handler data */
   SCIP_CALL( SCIPhashmapRemove(nlhdlrdata->exprmap, (void*)expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   --nlhdlrdata->nexprs;

   /* free nonlinear handler expression data */
   SCIPfreeBlockMemoryNull(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** callback to be called in initialization */
#define nlhdlrInitBilinear NULL

/** callback to be called in deinitialization */
static
SCIP_DECL_NLHDLREXIT(nlhdlrExitBilinear)
{  /*lint --e{715}*/
   assert(SCIPnlhdlrGetData(nlhdlr) != NULL);
   assert(SCIPnlhdlrGetData(nlhdlr)->nexprs == 0);

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectBilinear)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(expr != NULL);
   assert(participating != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata);

   /* only during solving will we have the extra inequalities that we rely on so much here */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE )
      return SCIP_OKAY;

   /* check for product expressions with two children */
   if( SCIPisExprProduct(scip, expr) && SCIPexprGetNChildren(expr) == 2
      && (nlhdlrdata->exprmap == NULL || !SCIPhashmapExists(nlhdlrdata->exprmap, (void*)expr)) )
   {
      SCIP_EXPR** children;
      SCIP_Bool valid;
      int c;

      children = SCIPexprGetChildren(expr);
      assert(children != NULL);

      /* detection is only successful if both children will have auxiliary variable or are variables
       * that are not binary variables */
      valid = TRUE;
      for( c = 0; c < 2; ++c )
      {
         assert(children[c] != NULL);
         if( SCIPgetExprNAuxvarUsesNonlinear(children[c]) == 0 &&
            (!SCIPisExprVar(scip, children[c]) || SCIPvarIsBinary(SCIPgetVarExprVar(children[c]))) )
         {
            valid = FALSE;
            break;
         }
      }

      if( valid )
      {
         /* create expression data for the nonlinear handler */
         SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
         (*nlhdlrexprdata)->lastnodeid = -1;

         /* ensure that there is enough memory to store the detected expression */
         if( nlhdlrdata->exprsize < nlhdlrdata->nexprs + 1 )
         {
            int newsize = SCIPcalcMemGrowSize(scip, nlhdlrdata->nexprs + 1);
            assert(newsize > nlhdlrdata->exprsize);

            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrdata->exprs, nlhdlrdata->exprsize, newsize) );
            nlhdlrdata->exprsize = newsize;
         }

         /* create expression map, if not done so far */
         if( nlhdlrdata->exprmap == NULL )
         {
            SCIP_CALL( SCIPhashmapCreate(&nlhdlrdata->exprmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
         }

#ifndef NDEBUG
         {
            int i;

            for( i = 0; i < nlhdlrdata->nexprs; ++i )
               assert(nlhdlrdata->exprs[i] != expr);
         }
#endif

         /* add expression to nlhdlrdata and capture it */
         nlhdlrdata->exprs[nlhdlrdata->nexprs] = expr;
         SCIPcaptureExpr(expr);
         SCIP_CALL( SCIPhashmapInsertInt(nlhdlrdata->exprmap, (void*)expr, nlhdlrdata->nexprs) );
         ++nlhdlrdata->nexprs;

         /* tell children that we will use their auxvar and use its activity for both estimate and domain propagation */
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, children[0], TRUE, nlhdlrdata->useinteval
               || nlhdlrdata->usereverseprop, TRUE, TRUE) );
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, children[1], TRUE, nlhdlrdata->useinteval
               || nlhdlrdata->usereverseprop, TRUE, TRUE) );
      }
   }

   if( *nlhdlrexprdata != NULL )
   {
      /* we want to join separation and domain propagation, if not disabled by parameter */
      *participating = SCIP_NLHDLR_METHOD_SEPABOTH;
      if( nlhdlrdata->useinteval || nlhdlrdata->usereverseprop )
         *participating |= SCIP_NLHDLR_METHOD_ACTIVITY;
   }

#ifdef SCIP_DEBUG
   if( *participating )
   {
      SCIPdebugMsg(scip, "detected expr ");
      SCIPprintExpr(scip, expr, NULL);
      SCIPinfoMessage(scip, NULL, " participating: %d\n", *participating);
   }
#endif

   return SCIP_OKAY;
}

/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxBilinear)
{ /*lint --e{715}*/
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_Real coef;

   assert(SCIPisExprProduct(scip, expr));
   assert(SCIPexprGetNChildren(expr) == 2);

   var1 = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[0]);
   assert(var1 != NULL);
   var2 = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[1]);
   assert(var2 != NULL);
   coef = SCIPgetCoefExprProduct(expr);

   *auxvalue = coef * SCIPgetSolVal(scip, sol, var1) * SCIPgetSolVal(scip, sol, var2);

   return SCIP_OKAY;
}

/** separation initialization method of a nonlinear handler (called during CONSINITLP) */
#define nlhdlrInitSepaBilinear NULL

/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
#define nlhdlrExitSepaBilinear NULL

/** nonlinear handler separation callback */
#define nlhdlrEnfoBilinear NULL

/** nonlinear handler under/overestimation callback */
static
SCIP_DECL_NLHDLRESTIMATE(nlhdlrEstimateBilinear)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_VAR* auxvar;
   SCIP_Real lincoefx = 0.0;
   SCIP_Real lincoefy = 0.0;
   SCIP_Real linconstant = 0.0;
   SCIP_Real refpointx;
   SCIP_Real refpointy;
   SCIP_Real violation;
   SCIP_Longint nodeid;
   SCIP_Bool mccsuccess = TRUE;
   SCIP_ROWPREP* rowprep;

   assert(rowpreps != NULL);

   *success = FALSE;
   *addedbranchscores = FALSE;

   /* check whether an inequality is available */
   if( nlhdlrexprdata->noverineqs == 0 && nlhdlrexprdata->nunderineqs == 0 )
      return SCIP_OKAY;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   nodeid = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   /* update last node */
   if( nlhdlrexprdata->lastnodeid != nodeid )
   {
      nlhdlrexprdata->lastnodeid = nodeid;
      nlhdlrexprdata->nseparoundslastnode = 0;
   }

   /* update separation round */
   ++nlhdlrexprdata->nseparoundslastnode;

   /* check working limits */
   if( (SCIPgetDepth(scip) == 0 && nlhdlrexprdata->nseparoundslastnode > nlhdlrdata->maxseparoundsroot)
         || (SCIPgetDepth(scip) > 0 && nlhdlrexprdata->nseparoundslastnode > nlhdlrdata->maxseparounds)
         || SCIPgetDepth(scip) > nlhdlrdata->maxsepadepth )
      return SCIP_OKAY;

   /* collect variables */
   x = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[0]);
   assert(x != NULL);
   y = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[1]);
   assert(y != NULL);
   auxvar = SCIPgetExprAuxVarNonlinear(expr);
   assert(auxvar != NULL);

   /* get and adjust the reference points */
   refpointx = MIN(MAX(SCIPgetSolVal(scip, sol, x), SCIPvarGetLbLocal(x)),SCIPvarGetUbLocal(x)); /*lint !e666*/
   refpointy = MIN(MAX(SCIPgetSolVal(scip, sol, y), SCIPvarGetLbLocal(y)),SCIPvarGetUbLocal(y)); /*lint !e666*/
   assert(SCIPisLE(scip, refpointx, SCIPvarGetUbLocal(x)) && SCIPisGE(scip, refpointx, SCIPvarGetLbLocal(x)));
   assert(SCIPisLE(scip, refpointy, SCIPvarGetUbLocal(y)) && SCIPisGE(scip, refpointy, SCIPvarGetLbLocal(y)));

   /* use McCormick inequalities to decide whether we want to separate or not */
   SCIPaddBilinMcCormick(scip, SCIPgetCoefExprProduct(expr), SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), refpointx,
         SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y), refpointy, overestimate, &lincoefx, &lincoefy, &linconstant,
         &mccsuccess);

   /* too large values in McCormick inequalities -> skip */
   if( !mccsuccess )
      return SCIP_OKAY;

   /* compute violation for the McCormick relaxation */
   violation = lincoefx * refpointx + lincoefy * refpointy + linconstant - SCIPgetSolVal(scip, sol, auxvar);
   if( overestimate )
      violation = -violation;

   /* only use a tighter relaxations if McCormick does not separate the reference point */
   if( SCIPisFeasLE(scip, violation, 0.0) && useBilinIneqs(scip, x, y, refpointx, refpointy) )
   {
      SCIP_Bool useoverestineq = SCIPgetCoefExprProduct(expr) > 0.0 ? overestimate : !overestimate;
      SCIP_Real mccormickval = lincoefx * refpointx + lincoefy * refpointy + linconstant;
      SCIP_Real* ineqs;
      SCIP_Real bestval;
      int nineqs;

      /* McCormick relaxation is too weak */
      bestval = mccormickval;

      /* get the inequalities that might lead to a tighter relaxation */
      if( useoverestineq )
      {
         ineqs = nlhdlrexprdata->overineqs;
         nineqs = nlhdlrexprdata->noverineqs;
      }
      else
      {
         ineqs = nlhdlrexprdata->underineqs;
         nineqs = nlhdlrexprdata->nunderineqs;
      }

      /* use linear inequalities to update relaxation */
      updateBilinearRelaxation(scip, x, y, SCIPgetCoefExprProduct(expr),
         overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT,
         refpointx, refpointy, ineqs, nineqs, mccormickval,
         &lincoefx, &lincoefy, &linconstant, &bestval,
         success);

#ifndef NDEBUG
      /* check whether cut is really valid */
      if( *success )
      {
         assert(!overestimate || SCIPisLE(scip, bestval, mccormickval));
         assert(overestimate || SCIPisGE(scip, bestval, mccormickval));
      }
#endif
   }

   if( *success )
   {
      SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );
      SCIProwprepAddConstant(rowprep, linconstant);
      SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, 2) );
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, x, lincoefx) );
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, y, lincoefy) );
      SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep) );
   }

   return SCIP_OKAY;
}

/** nonlinear handler interval evaluation callback */
static
SCIP_DECL_NLHDLRINTEVAL(nlhdlrIntevalBilinear)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   assert(nlhdlrexprdata != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( nlhdlrdata->useinteval && nlhdlrexprdata->nunderineqs + nlhdlrexprdata->noverineqs > 0 )
   {
      SCIP_INTERVAL tmp = intevalBilinear(scip, expr, nlhdlrexprdata->underineqs, nlhdlrexprdata->nunderineqs,
         nlhdlrexprdata->overineqs, nlhdlrexprdata->noverineqs);

      /* intersect intervals if we have learned a tighter interval */
      if( SCIPisGT(scip, tmp.inf, (*interval).inf) || SCIPisLT(scip, tmp.sup, (*interval).sup) )
         SCIPintervalIntersect(interval, *interval, tmp);
   }

   return SCIP_OKAY;
}

/** nonlinear handler callback for reverse propagation */
static
SCIP_DECL_NLHDLRREVERSEPROP(nlhdlrReversepropBilinear)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(nlhdlrexprdata != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( nlhdlrdata->usereverseprop && nlhdlrexprdata->nunderineqs + nlhdlrexprdata->noverineqs > 0 )
   {
      SCIP_EXPR* childx;
      SCIP_EXPR* childy;
      SCIP_INTERVAL intervalx;
      SCIP_INTERVAL intervaly;

      assert(SCIPexprGetNChildren(expr) == 2);
      childx = SCIPexprGetChildren(expr)[0];
      childy = SCIPexprGetChildren(expr)[1];
      assert(childx != NULL && childy != NULL);

      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &intervalx);
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &intervaly);

      /* compute bounds on x and y */
      reversePropBilinear(scip, conshdlr, expr, bounds, nlhdlrexprdata->underineqs, nlhdlrexprdata->nunderineqs,
         nlhdlrexprdata->overineqs, nlhdlrexprdata->noverineqs, &intervalx, &intervaly);

      /* tighten bounds of x */
      SCIPdebugMsg(scip, "try to tighten bounds of x: [%g,%g] -> [%g,%g]\n",
         SCIPgetExprBoundsNonlinear(scip, childx).inf, SCIPgetExprBoundsNonlinear(scip, childx).sup,
         intervalx.inf, intervalx.sup);

      SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, SCIPexprGetChildren(expr)[0], intervalx, infeasible,
         nreductions) );

      if( !(*infeasible) )
      {
         /* tighten bounds of y */
         SCIPdebugMsg(scip, "try to tighten bounds of y: [%g,%g] -> [%g,%g]\n",
            SCIPgetExprBoundsNonlinear(scip, childy).inf, SCIPgetExprBoundsNonlinear(scip, childy).sup,
            intervaly.inf, intervaly.sup);
         SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, SCIPexprGetChildren(expr)[1], intervaly,
            infeasible, nreductions) );
      }
   }

   return SCIP_OKAY;
}

/*
 * nonlinear handler specific interface methods
 */

/** includes bilinear nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrBilinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   /**! [SnippetIncludeNlhdlrBilinear] */
   /* create nonlinear handler specific data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   BMSclearMemory(nlhdlrdata);

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectBilinear, nlhdlrEvalauxBilinear, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrBilinear);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrdataBilinear);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataBilinear);
   SCIPnlhdlrSetInitExit(nlhdlr, nlhdlrInitBilinear, nlhdlrExitBilinear);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaBilinear, nlhdlrEnfoBilinear, nlhdlrEstimateBilinear, nlhdlrExitSepaBilinear);
   SCIPnlhdlrSetProp(nlhdlr, nlhdlrIntevalBilinear, nlhdlrReversepropBilinear);

   /* parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/useinteval",
         "whether to use the interval evaluation callback of the nlhdlr",
         &nlhdlrdata->useinteval, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/usereverseprop",
         "whether to use the reverse propagation callback of the nlhdlr",
         &nlhdlrdata->usereverseprop, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/maxseparoundsroot",
         "maximum number of separation rounds in the root node",
         &nlhdlrdata->maxseparoundsroot, FALSE, 10, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/maxseparounds",
         "maximum number of separation rounds in a local node",
         &nlhdlrdata->maxseparounds, FALSE, 1, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/maxsepadepth",
         "maximum depth to apply separation",
         &nlhdlrdata->maxsepadepth, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   /* statistic table */
   assert(SCIPfindTable(scip, TABLE_NAME_BILINEAR) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_BILINEAR, TABLE_DESC_BILINEAR, FALSE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputBilinear,
         NULL, TABLE_POSITION_BILINEAR, TABLE_EARLIEST_STAGE_BILINEAR) );
   /**! [SnippetIncludeNlhdlrBilinear] */

   return SCIP_OKAY;
}

/** returns an array of expressions that have been detected by the bilinear nonlinear handler */
SCIP_EXPR** SCIPgetExprsBilinear(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(nlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(nlhdlr), NLHDLR_NAME) == 0);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata);

   return nlhdlrdata->exprs;
}

/** returns the total number of expressions that have been detected by the bilinear nonlinear handler */
int SCIPgetNExprsBilinear(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(nlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(nlhdlr), NLHDLR_NAME) == 0);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata);

   return nlhdlrdata->nexprs;
}

/** adds a globally valid inequality of the form \f$\text{xcoef}\cdot x \leq \text{ycoef} \cdot y + \text{constant}\f$ to a product expression of the form \f$x\cdot y\f$ */
SCIP_RETCODE SCIPaddIneqBilinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_EXPR*            expr,               /**< product expression */
   SCIP_Real             xcoef,              /**< x coefficient */
   SCIP_Real             ycoef,              /**< y coefficient */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool*            success             /**< buffer to store whether inequality has been accepted */
   )
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_Real* ineqs;
   SCIP_Real viol1;
   SCIP_Real viol2;
   SCIP_Bool underestimate;
   int nineqs;
   int i;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(nlhdlr), NLHDLR_NAME) == 0);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 2);
   assert(xcoef != SCIP_INVALID); /*lint !e777 */
   assert(ycoef != SCIP_INVALID); /*lint !e777 */
   assert(constant != SCIP_INVALID); /*lint !e777 */
   assert(success != NULL);

   *success = FALSE;

   /* find nonlinear handler expression handler data */
   nlhdlrexprdata = SCIPgetNlhdlrExprDataNonlinear(nlhdlr, expr);

   if( nlhdlrexprdata == NULL )
   {
      SCIPwarningMessage(scip, "nonlinear expression data has not been found. "
            "Skip SCIPaddConsExprExprProductBilinearIneq()\n");
      return SCIP_OKAY;
   }

   /* ignore inequalities that only yield to a (possible) bound tightening */
   if( SCIPisFeasZero(scip, xcoef) || SCIPisFeasZero(scip, ycoef) )
      return SCIP_OKAY;

   /* collect variables */
   x = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[0]);
   y = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[1]);
   assert(x != NULL);
   assert(y != NULL);
   assert(x != y);

   /* normalize inequality s.t. xcoef in {-1,1} */
   if( !SCIPisEQ(scip, REALABS(xcoef), 1.0) )
   {
      constant /= REALABS(xcoef);
      ycoef /= REALABS(xcoef);
      xcoef /= REALABS(xcoef);
   }

   /* coefficients of the inequality determine whether the inequality can be used for under- or overestimation */
   underestimate = xcoef * ycoef > 0;

   SCIPdebugMsg(scip, "add inequality for a bilinear term: %g %s <= %g %s + %g (underestimate=%u)\n", xcoef,
      SCIPvarGetName(x), ycoef, SCIPvarGetName(y), constant, underestimate);

   /* compute violation of the inequality of the important corner points */
   getIneqViol(x, y, xcoef, ycoef, constant, &viol1, &viol2);
   SCIPdebugMsg(scip, "violations of inequality = (%g,%g)\n", viol1, viol2);

   /* inequality does not cutoff one of the important corner points -> skip */
   if( SCIPisFeasLE(scip, MAX(viol1, viol2), 0.0) )
      return SCIP_OKAY;

   if( underestimate )
   {
      ineqs = nlhdlrexprdata->underineqs;
      nineqs = nlhdlrexprdata->nunderineqs;
   }
   else
   {
      ineqs = nlhdlrexprdata->overineqs;
      nineqs = nlhdlrexprdata->noverineqs;
   }
   assert( nineqs >= 0 );

   /* check for a duplicate */
   for( i = 0; i < nineqs; ++i )
   {
      if( SCIPisFeasEQ(scip, xcoef, ineqs[3*i]) && SCIPisFeasEQ(scip, ycoef, ineqs[3*i+1]) /*lint !e661*/
         && SCIPisFeasEQ(scip, constant, ineqs[3*i+2]) )
      {
         SCIPdebugMsg(scip, "inequality already found -> skip\n");
         return SCIP_OKAY;
      }
   }

   {
      SCIP_Real viols1[2] = {0.0, 0.0};
      SCIP_Real viols2[2] = {0.0, 0.0};

      /* compute violations of existing inequalities */
      for( i = 0; i < nineqs; ++i )
      {
         getIneqViol(x, y, ineqs[3*i], ineqs[3*i+1], ineqs[3*i+2], &viols1[i], &viols2[i]); /*lint !e661*/

         /* check whether an existing inequality is dominating the candidate */
         if( SCIPisGE(scip, viols1[i], viol1) && SCIPisGE(scip, viols2[i], viol2) ) /*lint !e661*/
         {
            SCIPdebugMsg(scip, "inequality is dominated by %d -> skip\n", i);
            return SCIP_OKAY;
         }

         /* replace inequality if candidate is dominating it */
         if( SCIPisLT(scip, viols1[i], viol1) && SCIPisLT(scip, viols2[i], viol2) ) /*lint !e661*/
         {
            SCIPdebugMsg(scip, "inequality dominates %d -> replace\n", i);
            ineqs[3*i] = xcoef; /*lint !e661*/
            ineqs[3*i+1] = ycoef; /*lint !e661*/
            ineqs[3*i+2] = constant; /*lint !e661*/
            *success = TRUE;
         }
      }

      /* inequality is not dominated by other inequalities -> add if we have less than 2 inequalities */
      if( nineqs < 2 )
      {
         ineqs[3*nineqs] = xcoef;
         ineqs[3*nineqs + 1] = ycoef;
         ineqs[3*nineqs + 2] = constant;
         *success = TRUE;
         SCIPdebugMsg(scip, "add inequality\n");

         /* increase number of inequalities */
         if( underestimate )
            ++(nlhdlrexprdata->nunderineqs);
         else
            ++(nlhdlrexprdata->noverineqs);
      }
   }

   if( *success )
   {
      /* With the added inequalities, we can potentially compute tighter activities for the expression,
       * so constraints that contain this expression should be propagated again.
       * We don't have a direct expression to constraint mapping, though. This call marks all expr-constraints
       * which include any of the variables that this expression depends on for propagation.
       */
      SCIP_CALL( SCIPmarkExprPropagateNonlinear(scip, expr) );
   }

   return SCIP_OKAY;
}

/** computes coefficients of linearization of a bilinear term in a reference point */
void SCIPaddBilinLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             refpointx,          /**< point where to linearize first  variable */
   SCIP_Real             refpointy,          /**< point where to linearize second variable */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   SCIP_Real constant;

   assert(scip != NULL);
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( bilincoef == 0.0 )
      return;

   if( SCIPisInfinity(scip, REALABS(refpointx)) || SCIPisInfinity(scip, REALABS(refpointy)) )
   {
      *success = FALSE;
      return;
   }

   /* bilincoef * x * y ->  bilincoef * (refpointx * refpointy + refpointy * (x - refpointx) + refpointx * (y - refpointy))
    *                    = -bilincoef * refpointx * refpointy + bilincoef * refpointy * x + bilincoef * refpointx * y
    */

   constant = -bilincoef * refpointx * refpointy;

   if( SCIPisInfinity(scip, REALABS(bilincoef * refpointx)) || SCIPisInfinity(scip, REALABS(bilincoef * refpointy))
      || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoefx    += bilincoef * refpointy;
   *lincoefy    += bilincoef * refpointx;
   *linconstant += constant;
}

/** computes coefficients of McCormick under- or overestimation of a bilinear term */
void SCIPaddBilinMcCormick(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   )
{
   SCIP_Real constant;
   SCIP_Real coefx;
   SCIP_Real coefy;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lbx));
   assert(!SCIPisInfinity(scip, -ubx));
   assert(!SCIPisInfinity(scip,  lby));
   assert(!SCIPisInfinity(scip, -uby));
   assert(SCIPisInfinity(scip,  -lbx) || SCIPisLE(scip, lbx, ubx));
   assert(SCIPisInfinity(scip,  -lby) || SCIPisLE(scip, lby, uby));
   assert(SCIPisInfinity(scip,  -lbx) || SCIPisLE(scip, lbx, refpointx));
   assert(SCIPisInfinity(scip,  -lby) || SCIPisLE(scip, lby, refpointy));
   assert(SCIPisInfinity(scip,  ubx) || SCIPisGE(scip, ubx, refpointx));
   assert(SCIPisInfinity(scip,  uby) || SCIPisGE(scip, uby, refpointy));
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( bilincoef == 0.0 )
      return;

   if( overestimate )
      bilincoef = -bilincoef;

   if( SCIPisRelEQ(scip, lbx, ubx) && SCIPisRelEQ(scip, lby, uby) )
   {
      /* both x and y are mostly fixed */
      SCIP_Real cand1;
      SCIP_Real cand2;
      SCIP_Real cand3;
      SCIP_Real cand4;

      coefx = 0.0;
      coefy = 0.0;

      /* estimate x * y by constant */
      cand1 = lbx * lby;
      cand2 = lbx * uby;
      cand3 = ubx * lby;
      cand4 = ubx * uby;

      /* take most conservative value for underestimator */
      if( bilincoef < 0.0 )
         constant = bilincoef * MAX( MAX(cand1, cand2), MAX(cand3, cand4) );
      else
         constant = bilincoef * MIN( MIN(cand1, cand2), MIN(cand3, cand4) );
   }
   else if( bilincoef > 0.0 )
   {
      /* either x or y is not fixed and coef > 0.0 */
      if( !SCIPisInfinity(scip, -lbx) && !SCIPisInfinity(scip, -lby) &&
         (SCIPisInfinity(scip,  ubx) || SCIPisInfinity(scip,  uby)
            || (uby - refpointy) * (ubx - refpointx) >= (refpointy - lby) * (refpointx - lbx)) )
      {
         if( SCIPisRelEQ(scip, lbx, ubx) )
         {
            /* x*y = lbx * y + (x-lbx) * y >= lbx * y + (x-lbx) * lby >= lbx * y + min{(ubx-lbx) * lby, 0 * lby} */
            coefx    =  0.0;
            coefy    =  bilincoef * lbx;
            constant =  bilincoef * (lby < 0.0 ? (ubx-lbx) * lby : 0.0);
         }
         else if( SCIPisRelEQ(scip, lby, uby) )
         {
            /* x*y = lby * x + (y-lby) * x >= lby * x + (y-lby) * lbx >= lby * x + min{(uby-lby) * lbx, 0 * lbx} */
            coefx    =  bilincoef * lby;
            coefy    =  0.0;
            constant =  bilincoef * (lbx < 0.0 ? (uby-lby) * lbx : 0.0);
         }
         else
         {
            coefx    =  bilincoef * lby;
            coefy    =  bilincoef * lbx;
            constant = -bilincoef * lbx * lby;
         }
      }
      else if( !SCIPisInfinity(scip, ubx) && !SCIPisInfinity(scip, uby) )
      {
         if( SCIPisRelEQ(scip, lbx, ubx) )
         {
            /* x*y = ubx * y + (x-ubx) * y >= ubx * y + (x-ubx) * uby >= ubx * y + min{(lbx-ubx) * uby, 0 * uby} */
            coefx    =  0.0;
            coefy    =  bilincoef * ubx;
            constant =  bilincoef * (uby > 0.0 ? (lbx - ubx) * uby : 0.0);
         }
         else if( SCIPisRelEQ(scip, lby, uby) )
         {
            /* x*y = uby * x + (y-uby) * x >= uby * x + (y-uby) * ubx >= uby * x + min{(lby-uby) * ubx, 0 * ubx} */
            coefx    =  bilincoef * uby;
            coefy    =  0.0;
            constant =  bilincoef * (ubx > 0.0 ? (lby - uby) * ubx : 0.0);
         }
         else
         {
            coefx    =  bilincoef * uby;
            coefy    =  bilincoef * ubx;
            constant = -bilincoef * ubx * uby;
         }
      }
      else
      {
         *success = FALSE;
         return;
      }
   }
   else
   {
      /* either x or y is not fixed and coef < 0.0 */
      if( !SCIPisInfinity(scip,  ubx) && !SCIPisInfinity(scip, -lby) &&
         (SCIPisInfinity(scip, -lbx) || SCIPisInfinity(scip,  uby)
            || (ubx - lbx) * (refpointy - lby) <= (uby - lby) * (refpointx - lbx)) )
      {
         if( SCIPisRelEQ(scip, lbx, ubx) )
         {
            /* x*y = ubx * y + (x-ubx) * y <= ubx * y + (x-ubx) * lby <= ubx * y + max{(lbx-ubx) * lby, 0 * lby} */
            coefx    =  0.0;
            coefy    =  bilincoef * ubx;
            constant =  bilincoef * (lby < 0.0 ? (lbx - ubx) * lby : 0.0);
         }
         else if( SCIPisRelEQ(scip, lby, uby) )
         {
            /* x*y = lby * x + (y-lby) * x <= lby * x + (y-lby) * ubx <= lby * x + max{(uby-lby) * ubx, 0 * ubx} */
            coefx    =  bilincoef * lby;
            coefy    =  0.0;
            constant =  bilincoef * (ubx > 0.0 ? (uby - lby) * ubx : 0.0);
         }
         else
         {
            coefx    =  bilincoef * lby;
            coefy    =  bilincoef * ubx;
            constant = -bilincoef * ubx * lby;
         }
      }
      else if( !SCIPisInfinity(scip, -lbx) && !SCIPisInfinity(scip, uby) )
      {
         if( SCIPisRelEQ(scip, lbx, ubx) )
         {
            /* x*y = lbx * y + (x-lbx) * y <= lbx * y + (x-lbx) * uby <= lbx * y + max{(ubx-lbx) * uby, 0 * uby} */
            coefx    =  0.0;
            coefy    =  bilincoef * lbx;
            constant =  bilincoef * (uby > 0.0 ? (ubx - lbx) * uby : 0.0);
         }
         else if( SCIPisRelEQ(scip, lby, uby) )
         {
            /* x*y = uby * x + (y-uby) * x <= uby * x + (y-uby) * lbx <= uby * x + max{(lby-uby) * lbx, 0 * lbx} */
            coefx    =  bilincoef * uby;
            coefy    =  0.0;
            constant =  bilincoef * (lbx < 0.0 ? (lby - uby) * lbx : 0.0);
         }
         else
         {
            coefx    =  bilincoef * uby;
            coefy    =  bilincoef * lbx;
            constant = -bilincoef * lbx * uby;
         }
      }
      else
      {
         *success = FALSE;
         return;
      }
   }

   if( SCIPisInfinity(scip, REALABS(coefx)) || SCIPisInfinity(scip, REALABS(coefy))
      || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   if( overestimate )
   {
      coefx    = -coefx;
      coefy    = -coefy;
      constant = -constant;
   }

   SCIPdebugMsg(scip, "%.15g * x[%.15g,%.15g] * y[%.15g,%.15g] %c= %.15g * x %+.15g * y %+.15g\n", bilincoef, lbx, ubx,
      lby, uby, overestimate ? '<' : '>', coefx, coefy, constant);

   *lincoefx    += coefx;
   *lincoefy    += coefy;
   *linconstant += constant;
}

/** computes coefficients of linearization of a bilinear term in a reference point when given a linear inequality
 *  involving only the variables of the bilinear term
 *
 *  @note the formulas are extracted from "Convex envelopes of bivariate functions through the solution of KKT systems"
 *        by Marco Locatelli
 */
void SCIPcomputeBilinEnvelope1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real             xcoef,              /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             ycoef,              /**< y coefficient of linear inequality */
   SCIP_Real             constant,           /**< constant of linear inequality */
   SCIP_Real* RESTRICT   lincoefx,           /**< buffer to store coefficient of first  variable in linearization */
   SCIP_Real* RESTRICT   lincoefy,           /**< buffer to store coefficient of second variable in linearization */
   SCIP_Real* RESTRICT   linconstant,        /**< buffer to store constant of linearization */
   SCIP_Bool* RESTRICT   success             /**< buffer to store whether linearization was successful */
   )
{
   SCIP_Real xs[2] = {lbx, ubx};
   SCIP_Real ys[2] = {lby, uby};
   SCIP_Real minx;
   SCIP_Real maxx;
   SCIP_Real miny;
   SCIP_Real maxy;
   SCIP_Real QUAD(lincoefyq);
   SCIP_Real QUAD(lincoefxq);
   SCIP_Real QUAD(linconstantq);
   SCIP_Real QUAD(denomq);
   SCIP_Real QUAD(mjq);
   SCIP_Real QUAD(qjq);
   SCIP_Real QUAD(xjq);
   SCIP_Real QUAD(yjq);
   SCIP_Real QUAD(tmpq);
   SCIP_Real vx;
   SCIP_Real vy;
   int n;
   int i;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lbx));
   assert(!SCIPisInfinity(scip, -ubx));
   assert(!SCIPisInfinity(scip,  lby));
   assert(!SCIPisInfinity(scip, -uby));
   assert(SCIPisLE(scip, lbx, ubx));
   assert(SCIPisLE(scip, lby, uby));
   assert(SCIPisLE(scip, lbx, refpointx));
   assert(SCIPisGE(scip, ubx, refpointx));
   assert(SCIPisLE(scip, lby, refpointy));
   assert(SCIPisGE(scip, uby, refpointy));
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);
   assert(xcoef == 0.0 || xcoef == -1.0 || xcoef == 1.0); /*lint !e777*/
   assert(ycoef != SCIP_INVALID && ycoef != 0.0); /*lint !e777*/
   assert(constant != SCIP_INVALID); /*lint !e777*/

   *success = FALSE;
   *lincoefx = SCIP_INVALID;
   *lincoefy = SCIP_INVALID;
   *linconstant = SCIP_INVALID;

   /* reference point does not satisfy linear inequality */
   if( SCIPisFeasGT(scip, xcoef * refpointx - ycoef * refpointy - constant, 0.0) )
      return;

   /* compute minimal and maximal bounds on x and y for accepting the reference point */
   minx = lbx + 0.01 * (ubx-lbx);
   maxx = ubx - 0.01 * (ubx-lbx);
   miny = lby + 0.01 * (uby-lby);
   maxy = uby - 0.01 * (uby-lby);

   /* check whether the reference point is in [minx,maxx]x[miny,maxy] */
   if( SCIPisLE(scip, refpointx, minx) || SCIPisGE(scip, refpointx, maxx)
      || SCIPisLE(scip, refpointy, miny) || SCIPisGE(scip, refpointy, maxy) )
      return;

   /* always consider xy without the bilinear coefficient */
   if( bilincoef < 0.0 )
      overestimate = !overestimate;

   /* we use same notation as in "Convex envelopes of bivariate functions through the solution of KKT systems", 2016 */
   /* mj = xcoef / ycoef */
   SCIPquadprecDivDD(mjq, xcoef, ycoef);

   /* qj = -constant / ycoef */
   SCIPquadprecDivDD(qjq, -constant, ycoef);

   /* mj > 0 => underestimate; mj < 0 => overestimate */
   if( SCIPisNegative(scip, QUAD_TO_DBL(mjq)) != overestimate )
      return;

   /* get the corner point that satisfies the linear inequality xcoef*x <= ycoef*y + constant */
   if( !overestimate )
   {
      ys[0] = uby;
      ys[1] = lby;
   }

   vx = SCIP_INVALID;
   vy = SCIP_INVALID;
   n = 0;
   for( i = 0; i < 2; ++i )
   {
      SCIP_Real activity = xcoef * xs[i] - ycoef * ys[i] - constant;
      if( SCIPisLE(scip, activity, 0.0) )
      {
         /* corner point is satisfies inequality */
         vx = xs[i];
         vy = ys[i];
      }
      else if( SCIPisFeasGT(scip, activity, 0.0) )
         /* corner point is clearly cut off */
         ++n;
   }

   /* skip if no corner point satisfies the inequality or if no corner point is cut off
    * (that is, all corner points satisfy the inequality almost [1e-9..1e-6]) */
   if( n != 1 || vx == SCIP_INVALID || vy == SCIP_INVALID ) /*lint !e777*/
      return;

   /* denom = mj*(refpointx - vx) + vy - refpointy */
   SCIPquadprecSumDD(denomq, refpointx, -vx); /* refpoint - vx */
   SCIPquadprecProdQQ(denomq, denomq, mjq); /* mj * (refpoint - vx) */
   SCIPquadprecSumQD(denomq, denomq, vy); /* mj * (refpoint - vx) + vy */
   SCIPquadprecSumQD(denomq, denomq, -refpointy); /* mj * (refpoint - vx) + vy - refpointy */

   if( SCIPisZero(scip, QUAD_TO_DBL(denomq)) )
      return;

   /* (xj,yj) is the projection onto the line xcoef*x = ycoef*y + constant */
   /* xj = (refpointx*(vy - qj) - vx*(refpointy - qj)) / denom */
   SCIPquadprecProdQD(xjq, qjq, -1.0); /* - qj */
   SCIPquadprecSumQD(xjq, xjq, vy); /* vy - qj */
   SCIPquadprecProdQD(xjq, xjq, refpointx); /* refpointx * (vy - qj) */
   SCIPquadprecProdQD(tmpq, qjq, -1.0); /* - qj */
   SCIPquadprecSumQD(tmpq, tmpq, refpointy); /* refpointy - qj */
   SCIPquadprecProdQD(tmpq, tmpq, -vx); /* - vx * (refpointy - qj) */
   SCIPquadprecSumQQ(xjq, xjq, tmpq); /* refpointx * (vy - qj) - vx * (refpointy - qj) */
   SCIPquadprecDivQQ(xjq, xjq, denomq); /* (refpointx * (vy - qj) - vx * (refpointy - qj)) / denom */

   /* yj = mj * xj + qj */
   SCIPquadprecProdQQ(yjq, mjq, xjq);
   SCIPquadprecSumQQ(yjq, yjq, qjq);

   assert(SCIPisFeasEQ(scip, xcoef*QUAD_TO_DBL(xjq) - ycoef*QUAD_TO_DBL(yjq) - constant, 0.0));

   /* check whether the projection is in [minx,maxx] x [miny,maxy]; this avoids numerical difficulties when the
    * projection is close to the variable bounds
    */
   if( SCIPisLE(scip, QUAD_TO_DBL(xjq), minx) || SCIPisGE(scip, QUAD_TO_DBL(xjq), maxx)
      || SCIPisLE(scip, QUAD_TO_DBL(yjq), miny) || SCIPisGE(scip, QUAD_TO_DBL(yjq), maxy) )
      return;

   assert(vy - QUAD_TO_DBL(mjq)*vx - QUAD_TO_DBL(qjq) != 0.0);

   /* lincoefy = (mj*SQR(xj) - 2.0*mj*vx*xj - qj*vx + vx*vy) / (vy - mj*vx - qj) */
   SCIPquadprecSquareQ(lincoefyq, xjq); /* xj^2 */
   SCIPquadprecProdQQ(lincoefyq, lincoefyq, mjq); /* mj * xj^2 */
   SCIPquadprecProdQQ(tmpq, mjq, xjq); /* mj * xj */
   SCIPquadprecProdQD(tmpq, tmpq, -2.0 * vx); /* -2 * vx * mj * xj */
   SCIPquadprecSumQQ(lincoefyq, lincoefyq, tmpq); /* mj * xj^2 -2 * vx * mj * xj */
   SCIPquadprecProdQD(tmpq, qjq, -vx); /* -qj * vx */
   SCIPquadprecSumQQ(lincoefyq, lincoefyq, tmpq); /* mj * xj^2 -2 * vx * mj * xj -qj * vx */
   SCIPquadprecProdDD(tmpq, vx, vy); /* vx * vy */
   SCIPquadprecSumQQ(lincoefyq, lincoefyq, tmpq); /* mj * xj^2 -2 * vx * mj * xj -qj * vx + vx * vy */
   SCIPquadprecProdQD(tmpq, mjq, vx); /* mj * vx */
   SCIPquadprecSumQD(tmpq, tmpq, -vy); /* -vy + mj * vx */
   SCIPquadprecSumQQ(tmpq, tmpq, qjq); /* -vy + mj * vx + qj */
   QUAD_SCALE(tmpq, -1.0); /* vy - mj * vx - qj */
   SCIPquadprecDivQQ(lincoefyq, lincoefyq, tmpq); /* (mj * xj^2 -2 * vx * mj * xj -qj * vx + vx * vy) / (vy - mj * vx - qj) */

   /* lincoefx = 2.0*mj*xj + qj - mj*(*lincoefy) */
   SCIPquadprecProdQQ(lincoefxq, mjq, xjq); /* mj * xj */
   QUAD_SCALE(lincoefxq, 2.0); /* 2 * mj * xj */
   SCIPquadprecSumQQ(lincoefxq, lincoefxq, qjq); /* 2 * mj * xj + qj */
   SCIPquadprecProdQQ(tmpq, mjq, lincoefyq); /* mj * lincoefy */
   QUAD_SCALE(tmpq, -1.0); /* - mj * lincoefy */
   SCIPquadprecSumQQ(lincoefxq, lincoefxq, tmpq); /* 2 * mj * xj + qj - mj * lincoefy */

   /* linconstant = -mj*SQR(xj) - (*lincoefy)*qj */
   SCIPquadprecSquareQ(linconstantq, xjq); /* xj^2 */
   SCIPquadprecProdQQ(linconstantq, linconstantq, mjq); /* mj * xj^2 */
   QUAD_SCALE(linconstantq, -1.0); /* - mj * xj^2 */
   SCIPquadprecProdQQ(tmpq, lincoefyq, qjq); /* lincoefy * qj */
   QUAD_SCALE(tmpq, -1.0); /* - lincoefy * qj */
   SCIPquadprecSumQQ(linconstantq, linconstantq, tmpq); /* - mj * xj^2 - lincoefy * qj */

   /* consider the bilinear coefficient */
   SCIPquadprecProdQD(lincoefxq, lincoefxq, bilincoef);
   SCIPquadprecProdQD(lincoefyq, lincoefyq, bilincoef);
   SCIPquadprecProdQD(linconstantq, linconstantq, bilincoef);
   *lincoefx = QUAD_TO_DBL(lincoefxq);
   *lincoefy = QUAD_TO_DBL(lincoefyq);
   *linconstant = QUAD_TO_DBL(linconstantq);

   /* cut needs to be tight at (vx,vy) and (xj,yj); otherwise we consider the cut to be numerically bad */
   *success = SCIPisFeasEQ(scip, (*lincoefx)*vx + (*lincoefy)*vy + (*linconstant), bilincoef*vx*vy)
      && SCIPisFeasEQ(scip, (*lincoefx)*QUAD_TO_DBL(xjq) + (*lincoefy)*QUAD_TO_DBL(yjq) + (*linconstant),
      bilincoef*QUAD_TO_DBL(xjq)*QUAD_TO_DBL(yjq));

#ifndef NDEBUG
   {
      SCIP_Real activity = (*lincoefx)*refpointx + (*lincoefy)*refpointy + (*linconstant);

      /* cut needs to under- or overestimate the bilinear term at the reference point */
      if( bilincoef < 0.0 )
         overestimate = !overestimate;

      if( overestimate )
         assert(SCIPisFeasGE(scip, activity, bilincoef*refpointx*refpointy));
      else
         assert(SCIPisFeasLE(scip, activity, bilincoef*refpointx*refpointy));
   }
#endif
}

/** computes coefficients of linearization of a bilinear term in a reference point when given two linear inequalities
 *  involving only the variables of the bilinear term
 *
 *  @note the formulas are extracted from "Convex envelopes of bivariate functions through the solution of KKT systems"
 *        by Marco Locatelli
 */
void SCIPcomputeBilinEnvelope2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real             xcoef1,             /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             ycoef1,             /**< y coefficient of linear inequality */
   SCIP_Real             constant1,          /**< constant of linear inequality */
   SCIP_Real             xcoef2,             /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             ycoef2,             /**< y coefficient of linear inequality */
   SCIP_Real             constant2,          /**< constant of linear inequality */
   SCIP_Real* RESTRICT   lincoefx,           /**< buffer to store coefficient of first  variable in linearization */
   SCIP_Real* RESTRICT   lincoefy,           /**< buffer to store coefficient of second variable in linearization */
   SCIP_Real* RESTRICT   linconstant,        /**< buffer to store constant of linearization */
   SCIP_Bool* RESTRICT   success             /**< buffer to store whether linearization was successful */
   )
{
   SCIP_Real mi, mj, qi, qj, xi, xj, yi, yj;
   SCIP_Real xcoef, ycoef, constant;
   SCIP_Real minx, maxx, miny, maxy;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lbx));
   assert(!SCIPisInfinity(scip, -ubx));
   assert(!SCIPisInfinity(scip,  lby));
   assert(!SCIPisInfinity(scip, -uby));
   assert(SCIPisLE(scip, lbx, ubx));
   assert(SCIPisLE(scip, lby, uby));
   assert(SCIPisLE(scip, lbx, refpointx));
   assert(SCIPisGE(scip, ubx, refpointx));
   assert(SCIPisLE(scip, lby, refpointy));
   assert(SCIPisGE(scip, uby, refpointy));
   assert(lincoefx != NULL);
   assert(lincoefy != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);
   assert(xcoef1 != 0.0 && xcoef1 != SCIP_INVALID); /*lint !e777*/
   assert(ycoef1 != SCIP_INVALID && ycoef1 != 0.0); /*lint !e777*/
   assert(constant1 != SCIP_INVALID); /*lint !e777*/
   assert(xcoef2 != 0.0 && xcoef2 != SCIP_INVALID); /*lint !e777*/
   assert(ycoef2 != SCIP_INVALID && ycoef2 != 0.0); /*lint !e777*/
   assert(constant2 != SCIP_INVALID); /*lint !e777*/

   *success = FALSE;
   *lincoefx = SCIP_INVALID;
   *lincoefy = SCIP_INVALID;
   *linconstant = SCIP_INVALID;

   /* reference point does not satisfy linear inequalities */
   if( SCIPisFeasGT(scip, xcoef1 * refpointx - ycoef1 * refpointy - constant1, 0.0)
      || SCIPisFeasGT(scip, xcoef2 * refpointx - ycoef2 * refpointy - constant2, 0.0) )
      return;

   /* compute minimal and maximal bounds on x and y for accepting the reference point */
   minx = lbx + 0.01 * (ubx-lbx);
   maxx = ubx - 0.01 * (ubx-lbx);
   miny = lby + 0.01 * (uby-lby);
   maxy = uby - 0.01 * (uby-lby);

   /* check the reference point is in the interior of the domain */
   if( SCIPisLE(scip, refpointx, minx) || SCIPisGE(scip, refpointx, maxx)
      || SCIPisLE(scip, refpointy, miny) || SCIPisFeasGE(scip, refpointy, maxy) )
      return;

   /* the sign of the x-coefficients of the two inequalities must be different; otherwise the convex or concave
    * envelope can be computed via SCIPcomputeBilinEnvelope1 for each inequality separately
    */
   if( (xcoef1 > 0) == (xcoef2 > 0) )
      return;

   /* always consider xy without the bilinear coefficient */
   if( bilincoef < 0.0 )
      overestimate = !overestimate;

   /* we use same notation as in "Convex envelopes of bivariate functions through the solution of KKT systems", 2016 */
   mi = xcoef1 / ycoef1;
   qi = -constant1 / ycoef1;
   mj = xcoef2 / ycoef2;
   qj = -constant2 / ycoef2;

   /* mi, mj > 0 => underestimate; mi, mj < 0 => overestimate */
   if( SCIPisNegative(scip, mi) != overestimate || SCIPisNegative(scip, mj) != overestimate )
      return;

   /* compute cut according to Locatelli 2016 */
   computeBilinEnvelope2(scip, refpointx, refpointy, mi, qi, mj, qj, &xi, &yi, &xj, &yj, &xcoef, &ycoef, &constant);
   assert(SCIPisRelEQ(scip, mi*xi + qi, yi));
   assert(SCIPisRelEQ(scip, mj*xj + qj, yj));

   /* it might happen that (xi,yi) = (xj,yj) if the two lines intersect */
   if( SCIPisEQ(scip, xi, xj) && SCIPisEQ(scip, yi, yj) )
      return;

   /* check whether projected points are in the interior */
   if( SCIPisLE(scip, xi, minx) || SCIPisGE(scip, xi, maxx) || SCIPisLE(scip, yi, miny) || SCIPisGE(scip, yi, maxy) )
      return;
   if( SCIPisLE(scip, xj, minx) || SCIPisGE(scip, xj, maxx) || SCIPisLE(scip, yj, miny) || SCIPisGE(scip, yj, maxy) )
      return;

   *lincoefx = bilincoef * xcoef;
   *lincoefy = bilincoef * ycoef;
   *linconstant = bilincoef * constant;

   /* cut needs to be tight at (vx,vy) and (xj,yj) */
   *success = SCIPisFeasEQ(scip, (*lincoefx)*xi + (*lincoefy)*yi + (*linconstant), bilincoef*xi*yi)
      && SCIPisFeasEQ(scip, (*lincoefx)*xj + (*lincoefy)*yj + (*linconstant), bilincoef*xj*yj);

#ifndef NDEBUG
   {
      SCIP_Real activity = (*lincoefx)*refpointx + (*lincoefy)*refpointy + (*linconstant);

      /* cut needs to under- or overestimate the bilinear term at the reference point */
      if( bilincoef < 0.0 )
         overestimate = !overestimate;

      if( overestimate )
         assert(SCIPisFeasGE(scip, activity, bilincoef*refpointx*refpointy));
      else
         assert(SCIPisFeasLE(scip, activity, bilincoef*refpointx*refpointy));
   }
#endif
}
