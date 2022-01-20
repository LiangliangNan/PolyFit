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

/**@file   presol_implfree.c
 * @brief  perform multi-aggregation on implied free variables
 * @author Dieter Weninger
 *
 * This presolver tries to find implied free variables within equalities
 * which are convenient for multi-aggregation.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>


#include "scip/pub_matrix.h"
#include "presol_implfree.h"

#define PRESOL_NAME            "implfree"
#define PRESOL_DESC            "exploit implied free variables for multi-aggregation"
#define PRESOL_PRIORITY            -1000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               0     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define MAXABSRATIO     ((SCIP_Real)1000.0)     /**< max abs coefficients ratio */
#define SIDECHANGERATIO   ((SCIP_Real)10.0)     /**< max side change ratio */


/*
 * Local methods
 */

/** calculate max activity of one row without one column */
static
SCIP_Real getMaxActSingleRowWithoutCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col                 /**< column index */
   )
{
   int c;
   SCIP_Real val;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   SCIP_Real maxactivity;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(matrix != NULL);

   maxactivity = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      c = *rowpnt;
      val = *valpnt;

      if( c == col )
         continue;

      if( val > 0.0 )
      {
         ub = SCIPmatrixGetColUb(matrix, c);
         assert(!SCIPisInfinity(scip, ub));
         maxactivity += val * ub;
      }
      else if( val < 0.0 )
      {
         lb = SCIPmatrixGetColLb(matrix, c);
         assert(!SCIPisInfinity(scip, -lb));
         maxactivity += val * lb;
      }
   }

   return maxactivity;
}

/** calculate min activity of one row without one column */
static
SCIP_Real getMinActSingleRowWithoutCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col                 /**< column index */
   )
{
   int c;
   SCIP_Real val;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   SCIP_Real minactivity;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(matrix != NULL);

   minactivity = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      c = *rowpnt;
      val = *valpnt;

      if( c == col )
         continue;

      if( val > 0.0 )
      {
         lb = SCIPmatrixGetColLb(matrix, c);
         assert(!SCIPisInfinity(scip, -lb));
         minactivity += val * lb;
      }
      else if( val < 0.0 )
      {
         ub = SCIPmatrixGetColUb(matrix, c);
         assert(!SCIPisInfinity(scip, ub));
         minactivity += val * ub;
      }
   }

   return minactivity;
}

/** get min/max residual activity without the specified column */
static
void getActivityResiduals(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this variable in this row */
   SCIP_Real*            minresactivity,     /**< minimum residual activity of this row */
   SCIP_Real*            maxresactivity,     /**< maximum residual activity of this row */
   SCIP_Bool*            isminsettoinfinity, /**< flag indicating if minresactiviy is set to infinity */
   SCIP_Bool*            ismaxsettoinfinity  /**< flag indicating if maxresactiviy is set to infinity */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   int nmaxactneginf;
   int nmaxactposinf;
   int nminactneginf;
   int nminactposinf;
   SCIP_Real maxactivity;
   SCIP_Real minactivity;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);
   assert(isminsettoinfinity != NULL);
   assert(ismaxsettoinfinity != NULL);

   lb = SCIPmatrixGetColLb(matrix, col);
   ub = SCIPmatrixGetColUb(matrix, col);

   assert(!SCIPisInfinity(scip, lb));
   assert(!SCIPisInfinity(scip, -ub));

   *isminsettoinfinity = FALSE;
   *ismaxsettoinfinity = FALSE;

   nmaxactneginf = SCIPmatrixGetRowNMaxActNegInf(matrix, row);
   nmaxactposinf = SCIPmatrixGetRowNMaxActPosInf(matrix, row);
   nminactneginf = SCIPmatrixGetRowNMinActNegInf(matrix, row);
   nminactposinf = SCIPmatrixGetRowNMinActPosInf(matrix, row);

   maxactivity = SCIPmatrixGetRowMaxActivity(matrix, row);
   minactivity = SCIPmatrixGetRowMinActivity(matrix, row);

   if( val >= 0.0 )
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(nmaxactposinf >= 1);
         if( nmaxactposinf == 1 && nmaxactneginf == 0 )
            *maxresactivity = getMaxActSingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nmaxactneginf + nmaxactposinf) > 0 )
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
         else
            *maxresactivity = maxactivity - val * ub;
      }

      if( SCIPisInfinity(scip, -lb) )
      {
         assert(nminactneginf >= 1);
         if( nminactneginf == 1 && nminactposinf == 0 )
            *minresactivity = getMinActSingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nminactneginf + nminactposinf) > 0 )
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
         else
            *minresactivity = minactivity - val * lb;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(nmaxactneginf >= 1);
         if( nmaxactneginf == 1 && nmaxactposinf == 0 )
            *maxresactivity = getMaxActSingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nmaxactneginf + nmaxactposinf) > 0 )
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
         else
            *maxresactivity = maxactivity - val * lb;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         assert(nminactposinf >= 1);
         if( nminactposinf == 1 && nminactneginf == 0 )
            *minresactivity = getMinActSingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nminactneginf + nminactposinf) > 0 )
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
         else
            *minresactivity = minactivity - val * ub;
      }
   }
}

/** calculate the bounds of one variable from one row */
static
void getVarBoundsOfRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index of variable */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this column in this row */
   SCIP_Real*            rowlb,              /**< lower bound of row */
   SCIP_Bool*            lbfound,            /**< flag indicating that a lower bound was calculated */
   SCIP_Real*            rowub,              /**< upper bound of row */
   SCIP_Bool*            ubfound             /**< flag indicating that an upper bound was calculated */
   )
{
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(rowlb != NULL);
   assert(lbfound != NULL);
   assert(rowub != NULL);
   assert(ubfound != NULL);

   *rowlb = -SCIPinfinity(scip);
   *rowub = SCIPinfinity(scip);
   *lbfound = FALSE;
   *ubfound = FALSE;

   getActivityResiduals(scip, matrix, col, row, val, &minresactivity, &maxresactivity,
         &isminsettoinfinity, &ismaxsettoinfinity);

   lhs = SCIPmatrixGetRowLhs(matrix, row);
   rhs = SCIPmatrixGetRowRhs(matrix, row);

   if( val > 0.0 )
   {
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) )
      {
         *rowub = (rhs - minresactivity)/val;
         *ubfound = TRUE;
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) )
      {
         *rowlb = (lhs - maxresactivity)/val;
         *lbfound = TRUE;
      }
   }
   else
   {
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) )
      {
         *rowlb = (rhs - minresactivity)/val;
         *lbfound = TRUE;
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) )
      {
         *rowub = (lhs - maxresactivity)/val;
         *ubfound = TRUE;
      }
   }
}

/** verify if the variable is implied free */
static
SCIP_Bool isVarImpliedFree(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index for implied free test */
   int                   row,                /**< constraint planned for multi-aggregation */
   SCIP_Bool*            lockedcons,         /**< constraint not suitable for bound implication */
   int*                  impllbrowidx,       /**< row which implies the lower bound */
   int*                  implubrowidx        /**< row which implied the upper bound */
   )
{
   SCIP_Real varub;
   SCIP_Real varlb;
   SCIP_Real impliedub;
   SCIP_Real impliedlb;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;
   SCIP_Bool impliedfree;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lockedcons != NULL);
   assert(impllbrowidx != NULL);
   assert(implubrowidx != NULL);

   impliedfree = FALSE;

   varub = SCIPmatrixGetColUb(matrix, col);
   varlb = SCIPmatrixGetColLb(matrix, col);

   impliedub = SCIPinfinity(scip);
   impliedlb = -SCIPinfinity(scip);

   *impllbrowidx = -1;
   *implubrowidx = -1;

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for( ; (colpnt < colend); colpnt++, valpnt++ )
   {
      SCIP_Real rowlb;
      SCIP_Real rowub;
      SCIP_Bool lbfound;
      SCIP_Bool ubfound;

      if( lockedcons[*colpnt] )
      {
         impliedfree = FALSE;
         break;
      }

      if(*colpnt == row)
         continue;

      getVarBoundsOfRow(scip, matrix, col, *colpnt, *valpnt,
            &rowlb, &lbfound, &rowub, &ubfound);

      if( lbfound && rowlb > impliedlb )
      {
         impliedlb = rowlb;
         *impllbrowidx = *colpnt;
      }

      if( ubfound && rowub < impliedub )
      {
         impliedub = rowub;
         *implubrowidx = *colpnt;
      }
   }

   if( !SCIPisInfinity(scip, -varlb) && !SCIPisInfinity(scip, varub) &&
       SCIPisFeasLE(scip, impliedub, varub) && SCIPisFeasLE(scip, varlb, impliedlb) )
   {
      impliedfree = TRUE;
   }

   return impliedfree;
}

/** calculate the amount of fill-in getting from multi-aggregation */
static
SCIP_Real getFillIn(
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   col,                /**< column index */
   int                   row                 /**< row index */
   )
{
   SCIP_Real fillin;
   int rowcnt;
   int colcnt;
   int nonzerosnew;
   int nonzerosold;

   assert(matrix != NULL);

   rowcnt = SCIPmatrixGetRowNNonzs(matrix, row);
   colcnt = SCIPmatrixGetColNNonzs(matrix, col);
   assert(rowcnt > 0 && colcnt > 0);

   nonzerosold = rowcnt + colcnt - 1;
   nonzerosnew = (colcnt - 1) * (rowcnt - 1);

   fillin = (SCIP_Real)nonzerosnew / (SCIP_Real)nonzerosold;

   return fillin;
}

/** verify if the multi-aggregation is numerical stable concerning the side change */
static
SCIP_Real sideChangeNumericalStable(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_Real             oldside,            /**< old rhs or lhs of constraint */
   SCIP_Real             aggrconst,          /**< multiaggrconstant */
   SCIP_Real             val                 /**< coefficient of multiaggvariable not within the aggregated constraint */
   )
{
   SCIP_Real enumerator;
   SCIP_Real denominator;

   enumerator = REALABS( -(val * aggrconst) - oldside );
   denominator = MAX(1.0,REALABS(oldside)); /*lint !e666*/

   return enumerator/denominator <= SIDECHANGERATIO;
}

/** calculate the huge contribution counters */
static
void getNumHugeActivities(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int*                  maxactposhuge,      /**< max activity positive contribution counter */
   int*                  maxactneghuge,      /**< max activity negative contribution counter */
   int*                  minactposhuge,      /**< min activity positive contribution counter */
   int*                  minactneghuge       /**< min activity negative contribution counter */
   )
{
   SCIP_Real val;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   int col;
   int row;
   SCIP_Real lb;
   SCIP_Real ub;
   int nrows;

   assert(scip != NULL);
   assert(matrix != NULL);

   nrows = SCIPmatrixGetNRows(matrix);

   for( row = 0; row < nrows; row++ )
   {
      maxactposhuge[row] = 0;
      maxactneghuge[row] = 0;
      minactposhuge[row] = 0;
      minactneghuge[row] = 0;

      rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
      rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
      valpnt = SCIPmatrixGetRowValPtr(matrix, row);

      for( ; rowpnt < rowend; rowpnt++, valpnt++ )
      {
         col = *rowpnt;
         val = *valpnt;

         lb = SCIPmatrixGetColLb(matrix, col);
         ub = SCIPmatrixGetColUb(matrix, col);

         if( val > 0.0 )
         {
            if( !SCIPisInfinity(scip, ub) )
            {
               if( SCIPisHugeValue(scip, val * ub) )
                  maxactposhuge[row]++;
               else if( SCIPisHugeValue(scip, -val * ub) )
                  maxactneghuge[row]++;
            }

            if( !SCIPisInfinity(scip, -lb) )
            {
               if( SCIPisHugeValue(scip, val * lb) )
                  minactposhuge[row]++;
               else if( SCIPisHugeValue(scip, -val * lb) )
                  minactneghuge[row]++;
            }
         }
         else
         {
            if( !SCIPisInfinity(scip, -lb) )
            {
               if( SCIPisHugeValue(scip, val * lb) )
                  maxactposhuge[row]++;
               else if( SCIPisHugeValue(scip, -val * lb) )
                  maxactneghuge[row]++;
            }

            if( !SCIPisInfinity(scip, ub) )
            {
               if( SCIPisHugeValue(scip, val * ub) )
                  minactposhuge[row]++;
               else if( SCIPisHugeValue(scip, -val * ub) )
                  minactneghuge[row]++;
            }
         }
      }
   }
}

/** verify if activities could be inexact */
static
void getActivityRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   row,                /**< row index */
   int                   col,                /**< column index */
   SCIP_Real             val,                /**< coefficient value of variable in linear constraint */
   int*                  maxactposhuge,      /**< max activity positive contribution counter */
   int*                  maxactneghuge,      /**< max activity negative contribution counter */
   int*                  minactposhuge,      /**< min activity positive contribution counter */
   int*                  minactneghuge,      /**< min activity negative contribution counter */
   SCIP_Bool*            minisrelax,         /**< flag indicating min activity could be inexact */
   SCIP_Bool*            maxisrelax          /**< flag indicating max activity could be inexact */
   )
{
   int numinf;
   int numhuge;
   SCIP_Real lb;
   SCIP_Real ub;
   int nmaxactposinf;
   int nmaxactneginf;
   int nminactposinf;
   int nminactneginf;

   *minisrelax = FALSE;
   *maxisrelax = FALSE;

   nmaxactposinf = SCIPmatrixGetRowNMaxActPosInf(matrix, row);
   nmaxactneginf = SCIPmatrixGetRowNMaxActNegInf(matrix, row);
   nminactposinf = SCIPmatrixGetRowNMinActPosInf(matrix, row);
   nminactneginf = SCIPmatrixGetRowNMinActNegInf(matrix, row);

   assert(nmaxactposinf >= 0);
   assert(nmaxactneginf >= 0);
   assert(nminactposinf >= 0);
   assert(nminactneginf >= 0);
   assert(maxactposhuge[row] >= 0);
   assert(maxactneghuge[row] >= 0);
   assert(minactposhuge[row] >= 0);
   assert(minactneghuge[row] >= 0);

   numinf = nmaxactposinf + nmaxactneginf + nminactposinf + nminactneginf;
   numhuge = maxactposhuge[row] + maxactneghuge[row] + minactposhuge[row] + minactneghuge[row];

   lb = SCIPmatrixGetColLb(matrix, col);
   ub = SCIPmatrixGetColUb(matrix, col);

   if( val > 0 )
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(nmaxactposinf >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *maxisrelax = TRUE;
      }
      else if( SCIPisHugeValue(scip, ub * val) )
      {
         assert(maxactposhuge[row] >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *maxisrelax = TRUE;
      }
      else if(SCIPisHugeValue(scip, ub * -val) )
      {
         assert(maxactneghuge[row] >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *maxisrelax = TRUE;
      }
      else
      {
         if( numinf + numhuge >= 1 )
            *maxisrelax = TRUE;
      }

      if( SCIPisInfinity(scip, -lb) )
      {
         assert(nminactneginf >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *minisrelax = TRUE;
      }
      else if( SCIPisHugeValue(scip, lb * val) )
      {
         assert(minactposhuge[row] >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *minisrelax = TRUE;
      }
      else if(SCIPisHugeValue(scip, lb * -val) )
      {
         assert(minactneghuge[row] >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *minisrelax = TRUE;
      }
      else
      {
         if( numinf + numhuge >= 1 )
            *minisrelax = TRUE;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(nmaxactneginf >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *maxisrelax = TRUE;
      }
      else if( SCIPisHugeValue(scip, lb * val) )
      {
         assert(maxactposhuge[row] >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *maxisrelax = TRUE;
      }
      else if(SCIPisHugeValue(scip, lb * -val) )
      {
         assert(maxactneghuge[row] >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *maxisrelax = TRUE;
      }
      else
      {
         if( numinf + numhuge >= 1 )
            *maxisrelax = TRUE;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         assert(nminactposinf >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *minisrelax = TRUE;
      }
      else if( SCIPisHugeValue(scip, ub * val) )
      {
         assert(minactposhuge[row] >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *minisrelax = TRUE;
      }
      else if(SCIPisHugeValue(scip, ub * -val) )
      {
         assert(minactneghuge[row] >= 1);
         if( numinf + numhuge - 1 >= 1 )
            *minisrelax = TRUE;
      }
      else
      {
         if( numinf + numhuge >= 1 )
            *minisrelax = TRUE;
      }
   }
}


/** verify if the multi-aggregation is numerical stable */
static
SCIP_Bool numericalStable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             aggrconst,          /**< aggregation constant */
   int*                  maxactposhuge,      /**< max activity positive contribution counter */
   int*                  maxactneghuge,      /**< max activity negative contribution counter */
   int*                  minactposhuge,      /**< min activity positive contribution counter */
   int*                  minactneghuge       /**< min activity negative contribution counter */
   )
{
   int* rowpnt;
   int* rowend;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;
   SCIP_Bool numericalstable;
   SCIP_Bool minisrelax;
   SCIP_Bool maxisrelax;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real minabsval;
   SCIP_Real maxabsval;

   numericalstable = TRUE;

   minabsval = SCIPinfinity(scip);
   maxabsval = -1.0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for( ; rowpnt < rowend; rowpnt++, valpnt++ )
   {
      SCIP_Real absval;

      absval = REALABS(*valpnt);
      if( absval < minabsval )
         minabsval = absval;
      if( absval > maxabsval )
         maxabsval = absval;
   }

   /* 1. criterion:
    * do not allow multi-aggregation if the abs coefficients differ too much
    */
   if( maxabsval / minabsval > MAXABSRATIO )
   {
      numericalstable = FALSE;
      return numericalstable;
   }

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for( ; colpnt < colend; colpnt++, valpnt++ )
   {
      int r;
      SCIP_Real val;

      r = *colpnt;
      val = *valpnt;

      /* 2. criterion:
       * do not allow huge contributions to the activities
       */
      getActivityRelax(scip, matrix, r, col, val,
         maxactposhuge, maxactneghuge, minactposhuge, minactneghuge,
         &minisrelax, &maxisrelax);

      if( minisrelax || maxisrelax )
      {
         numericalstable = FALSE;
         break;
      }

      if( r == row )
         continue;

      /* 3. criterion:
       * do not allow multi-aggregation if great sides changes appear
       */
      lhs = SCIPmatrixGetRowLhs(matrix, r);
      if( !SCIPisInfinity(scip, -lhs) )
      {
         if( !sideChangeNumericalStable(scip, lhs, aggrconst, val) )
         {
            numericalstable = FALSE;
            break;
         }
      }
      rhs = SCIPmatrixGetRowRhs(matrix, r);
      if( !SCIPisInfinity(scip, rhs) )
      {
         if( !sideChangeNumericalStable(scip, rhs, aggrconst, val) )
         {
            numericalstable = FALSE;
            break;
         }
      }
   }

   return numericalstable;
}

/** calculate sides after aggregation */
static
void calcNewSidesAfterAggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             coef,               /**< coefficient of aggregated variable */
   SCIP_Real*            newlhs,             /**< new calculated left hand side */
   SCIP_Real*            newrhs              /**< new calculated right hand side */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(newlhs != NULL);
   assert(newrhs != NULL);

   lhs = SCIPmatrixGetRowLhs(matrix, row);
   rhs = SCIPmatrixGetRowRhs(matrix, row);

   lb = SCIPmatrixGetColLb(matrix, col);
   ub = SCIPmatrixGetColUb(matrix, col);

   if( coef > 0.0 )
   {
      if( SCIPisInfinity(scip, -lb) )
         *newrhs = SCIPinfinity(scip);
      else
         *newrhs = rhs - coef * lb;

      if( SCIPisInfinity(scip, ub) )
         *newlhs = -SCIPinfinity(scip);
      else
         *newlhs = lhs - coef * ub;
   }
   else
   {
      if( SCIPisInfinity(scip, -lb) )
         *newlhs = -SCIPinfinity(scip);
      else
         *newlhs = rhs - coef * lb;

      if( SCIPisInfinity(scip, ub) )
         *newrhs = SCIPinfinity(scip);
      else
         *newrhs = lhs - coef * ub;
   }

   assert(SCIPisLE(scip, *newlhs, *newrhs));
}

/** verify if the constraint becomes redundant after aggregation */
static
SCIP_Bool isConsRedundant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             coef                /**< coefficient of aggregated variable */
   )
{
   SCIP_Real newlhs;
   SCIP_Real newrhs;
   SCIP_Bool consredundant;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;

   consredundant = FALSE;

   calcNewSidesAfterAggregation(scip, matrix, col, row, coef, &newlhs, &newrhs);
   getActivityResiduals(scip, matrix, col, row, coef, &minresactivity, &maxresactivity, &isminsettoinfinity, &ismaxsettoinfinity);

   if( !isminsettoinfinity && !ismaxsettoinfinity )
      consredundant = (SCIPisFeasLE(scip, newlhs, minresactivity) && SCIPisFeasLE(scip, maxresactivity, newrhs));

   return consredundant;
}

/** identify candidates for multi-aggregation */
static
void getMultiaggVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   SCIP_Bool*            multiaggvars,       /**< array indicating multi-aggregatable variables */
   int*                  nummultiaggvars,    /**< number of multi-aggregatable variables */
   int*                  multiaggequalities, /**< array holding aggregation equality row indices */
   SCIP_Bool*            consredundant,      /**< array indicating which constraint became redundant */
   SCIP_Bool*            lockedcons,         /**< constraints which could not be used for bound implication */
   SCIP_Bool*            skipvars,           /**< array holding which variables should be investigated */
   int*                  maxactposhuge,      /**< max activity positive contribution counter */
   int*                  maxactneghuge,      /**< max activity negative contribution counter */
   int*                  minactposhuge,      /**< min activity positive contribution counter */
   int*                  minactneghuge       /**< min activity negative contribution counter */
   )
{
   int r;
   SCIP_Real bestfillin;
   int bestvaridx;
   SCIP_Bool bestconsredundant;
   int tmpimpllbrowidx;
   int tmpimplubrowidx;
   int bestimpllbrowidx;
   int bestimplubrowidx;
   int nrows;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(multiaggvars != NULL);
   assert(nummultiaggvars != NULL);
   assert(multiaggequalities != NULL);
   assert(consredundant != NULL);
   assert(lockedcons != NULL);
   assert(skipvars != NULL);

   nrows = SCIPmatrixGetNRows(matrix);

   for( r = 0; r < nrows; r++ )
   {
      /* consider only equalities with min three variables */
      if( SCIPmatrixGetRowNNonzs(matrix, r) > 2 &&
         SCIPisEQ(scip, SCIPmatrixGetRowLhs(matrix, r), SCIPmatrixGetRowRhs(matrix, r)) )
      {
         int* rowpnt;
         int* rowend;
         SCIP_Real* valpnt;

         rowpnt = SCIPmatrixGetRowIdxPtr(matrix, r);
         rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, r);
         valpnt = SCIPmatrixGetRowValPtr(matrix, r);

         bestfillin = 1.0;
         bestvaridx = -1;
         bestimpllbrowidx = -1; /* only for lint */
         bestimplubrowidx = -1; /* only for lint */
         bestconsredundant = FALSE;

         for( ; rowpnt < rowend; rowpnt++, valpnt++ )
         {
            SCIP_VAR* var;
            SCIP_Real aggrconst;
            SCIP_Real fillin;

            var = SCIPmatrixGetVar(matrix, *rowpnt);

            if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
               continue;

            /* search for a continuous variable for aggregation which is implied free,
               produces less fill-in and leads to numerical stability */
            if( !skipvars[*rowpnt] && isVarImpliedFree(scip, matrix, *rowpnt, r, lockedcons, &tmpimpllbrowidx, &tmpimplubrowidx) )
            {
               assert(tmpimpllbrowidx >= 0 && tmpimplubrowidx >= 0);
               aggrconst = SCIPmatrixGetRowRhs(matrix, r) / (*valpnt);

               if( numericalStable(scip, matrix, *rowpnt, r, aggrconst,
                     maxactposhuge, maxactneghuge, minactposhuge, minactneghuge) )
               {
                  fillin = getFillIn(matrix, *rowpnt, r);
                  if( isConsRedundant(scip, matrix, *rowpnt, r, *valpnt) )
                  {
                     if( (fillin < bestfillin) || (!bestconsredundant && fillin <= 1.0) )
                     {
                        bestfillin = fillin;
                        bestvaridx = *rowpnt;
                        bestimpllbrowidx = tmpimpllbrowidx;
                        bestimplubrowidx = tmpimplubrowidx;
                        bestconsredundant = TRUE;
                     }
                  }
                  else
                  {
                     /* prefer redundant constraints */
                     if( bestconsredundant )
                        continue;

                     if( fillin < bestfillin )
                     {
                        bestfillin = fillin;
                        bestvaridx = *rowpnt;
                        bestimpllbrowidx = tmpimpllbrowidx;
                        bestimplubrowidx = tmpimplubrowidx;
                     }
                  }
               }
            }
            else
            {
               skipvars[*rowpnt] = TRUE;
            }
         }

         if( bestvaridx > -1 && multiaggvars[bestvaridx] != TRUE )
         {
            assert(bestvaridx < SCIPmatrixGetNColumns(matrix));
            multiaggvars[bestvaridx] = TRUE;
            multiaggequalities[bestvaridx] = r;
            lockedcons[bestimpllbrowidx] = TRUE;
            lockedcons[bestimplubrowidx] = TRUE;
            (*nummultiaggvars)++;

            if( bestconsredundant )
            {
               consredundant[r] = TRUE;
               lockedcons[r] = TRUE;
            }
         }
      }
   }
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyImplfree)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolImplfree(scip) );

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecImplfree)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPgetNVars(scip) == 0 || SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   if( SCIPgetNContVars(scip) == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized && complete )
   {
      SCIP_Bool* isvartomultiagg;
      int* multiaggequalities;
      int nummultiaggvars;
      SCIP_Bool* consredundant;
      SCIP_Bool* lockedcons;
      SCIP_Bool* skipvars;
      int nrows;
      int ncols;
      int* maxactposhuge;
      int* maxactneghuge;
      int* minactposhuge;
      int* minactneghuge;

      nrows = SCIPmatrixGetNRows(matrix);
      ncols = SCIPmatrixGetNColumns(matrix);
      assert(SCIPgetNVars(scip) == ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &isvartomultiagg, ncols) );
      BMSclearMemoryArray(isvartomultiagg, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &multiaggequalities, ncols) );
      BMSclearMemoryArray(multiaggequalities, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &consredundant, nrows) );
      BMSclearMemoryArray(consredundant, nrows);

      SCIP_CALL( SCIPallocBufferArray(scip, &lockedcons, nrows) );
      BMSclearMemoryArray(lockedcons, nrows);

      SCIP_CALL( SCIPallocBufferArray(scip, &skipvars, ncols) );
      BMSclearMemoryArray(skipvars, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &maxactposhuge, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &maxactneghuge, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &minactposhuge, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &minactneghuge, nrows) );

      getNumHugeActivities(scip, matrix, maxactposhuge, maxactneghuge, minactposhuge, minactneghuge);

      nummultiaggvars = 0;
      getMultiaggVars(scip, matrix, isvartomultiagg,&nummultiaggvars, multiaggequalities,
         consredundant, lockedcons, skipvars, maxactposhuge, maxactneghuge, minactposhuge, minactneghuge);

      if( nummultiaggvars > 0 )
      {
         int v;
         for( v = 0; v < ncols; v++ )
         {
            if( isvartomultiagg[v] )
            {
               SCIP_VAR* multiaggvar;
               SCIP_VAR** vars;
               SCIP_Real* scalars;
               SCIP_Real multiaggcoef;
               SCIP_Real aggrconst;
               SCIP_Bool infeasible;
               SCIP_Bool aggregated;
               SCIP_CONS* multiaggcons;
               int row;
               int cnt;
               int* rowpnt;
               int* rowend;
               SCIP_Real* valpnt;
               int nvars;
               SCIP_Real rhs;

               multiaggvar = SCIPmatrixGetVar(matrix, v);
               assert(SCIPvarGetType(multiaggvar) == SCIP_VARTYPE_CONTINUOUS);

               row = multiaggequalities[v];
               assert(row < nrows);

               multiaggcons = SCIPmatrixGetCons(matrix, row);
               rhs = SCIPmatrixGetRowRhs(matrix, row);

               /* get the number of variables without the multi-agg variable itself */
               nvars = SCIPmatrixGetRowNNonzs(matrix, row) - 1;

               SCIP_CALL( SCIPallocBufferArray(scip, &scalars, nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

               /* get multi-agg variable coefficient and vars */
               rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
               rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
               valpnt = SCIPmatrixGetRowValPtr(matrix, row);
               cnt = 0;
               multiaggcoef = 0.0;
               for( ; rowpnt < rowend; rowpnt++, valpnt++ )
               {
                  if( *rowpnt == v )
                  {
                     multiaggcoef = *valpnt;
                     continue;
                  }

                  vars[cnt++] = SCIPmatrixGetVar(matrix, *rowpnt);
               }

               /* avoid division by zero */
               if( SCIPisEQ(scip, multiaggcoef, 0.0) )
                  continue;

               assert(SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, row)) ||
                  SCIPisInfinity(scip, rhs) ||
                  SCIPisEQ(scip, SCIPmatrixGetRowLhs(matrix, row), rhs));
               assert(multiaggcoef != 0.0);

               /* we have to distinguished two cases */
               if( !SCIPisInfinity(scip, rhs) )
                  aggrconst = rhs / multiaggcoef;
               else
                  aggrconst = SCIPmatrixGetRowLhs(matrix, row) / multiaggcoef;

               /* calculate scalars */
               rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
               rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
               valpnt = SCIPmatrixGetRowValPtr(matrix, row);
               cnt = 0;
               SCIPdebugMsg(scip, "constraint <%s>: multi-aggregate <%s> ==", SCIPconsGetName(multiaggcons), SCIPvarGetName(multiaggvar));
               for( ; rowpnt < rowend; rowpnt++, valpnt++ )
               {
                  if( *rowpnt == v )
                     continue;

                  scalars[cnt] = -(*valpnt) / multiaggcoef; /*lint !e414*/
                  SCIPdebugMsgPrint(scip, " %+.15g<%s>", scalars[cnt], SCIPvarGetName(SCIPmatrixGetVar(matrix, *rowpnt)));

                  cnt++;
               }

               SCIPdebugMsgPrint(scip, " %+.15g, bounds of <%s>: [%.15g,%.15g]\n",
                  aggrconst, SCIPvarGetName(multiaggvar), SCIPvarGetLbGlobal(multiaggvar), SCIPvarGetUbGlobal(multiaggvar));

               /* perform multi-aggregation */
               SCIP_CALL( SCIPmultiaggregateVar(scip, multiaggvar, nvars, vars, scalars, aggrconst, &infeasible, &aggregated) );
               assert(aggregated);

               SCIPfreeBufferArray(scip, &vars);
               SCIPfreeBufferArray(scip, &scalars);

               /* check for infeasible aggregation */
               if( infeasible )
               {
                  SCIPdebugMsg(scip, "constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(multiaggcons));
                  return SCIP_OKAY;
               }

               (*naggrvars)++;

               /* remove constraint if it is redundant */
               if( consredundant[row] )
               {
                  SCIP_CALL( SCIPdelCons(scip, multiaggcons) );
                  (*ndelconss)++;
               }
            }
         }
      }

      SCIPfreeBufferArray(scip, &minactneghuge);
      SCIPfreeBufferArray(scip, &minactposhuge);
      SCIPfreeBufferArray(scip, &maxactneghuge);
      SCIPfreeBufferArray(scip, &maxactposhuge);
      SCIPfreeBufferArray(scip, &skipvars);
      SCIPfreeBufferArray(scip, &lockedcons);
      SCIPfreeBufferArray(scip, &consredundant);
      SCIPfreeBufferArray(scip, &multiaggequalities);
      SCIPfreeBufferArray(scip, &isvartomultiagg);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the implied free presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolImplfree(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecImplfree, NULL) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyImplfree) );

   return SCIP_OKAY;
}
