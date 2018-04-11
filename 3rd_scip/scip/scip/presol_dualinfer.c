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

/**@file    presol_dualinfer.c
 * @ingroup PRESOLVERS
 * @brief   dual inference presolver
 * @author  Dieter Weninger
 *
 * This presolver does bound strengthening on continuous variables
 * for getting better bounds on dual variables y. The bounds on the dual
 * variables are then used to derive variable fixings or side changes.
 * We distinguish two cases:
 * i)  positive reduced costs:    c_j - sup{A_{.j}^T y} > 0 => x_j = l_j
 * ii) positive dual lower bound: y_i > 0 =>  A_{i.}x = b_i
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/pub_matrix.h"
#include "scip/cons_linear.h"
#include "presol_dualinfer.h"

#define PRESOL_NAME             "dualinfer"
#define PRESOL_DESC             "exploit dual informations for fixings and side changes"
#define PRESOL_PRIORITY             -2000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS                0    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */
#define MAX_LOOPS                       3    /**< maximal number of dual bound strengthening loops */


/*
 * Data structures
 */


/** type of fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,
   NOFIX   =  0
};
typedef enum Fixingdirection FIXINGDIRECTION;

/** type of side change */
enum SideChange
{
   RHSTOLHS = -1,
   NOCHANGE = 0,
   LHSTORHS = 1
};
typedef enum SideChange SIDECHANGE;


/*
 * Local methods
 */

/** calculate max activity of one row without one column */
static
SCIP_Real getMaxActivitySingleRowWithoutCol(
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
SCIP_Real getMinActivitySingleRowWithoutCol(
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
void getMinMaxActivityResiduals(
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
            *maxresactivity = getMaxActivitySingleRowWithoutCol(scip, matrix, row, col);
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
            *minresactivity = getMinActivitySingleRowWithoutCol(scip, matrix, row, col);
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
            *maxresactivity = getMaxActivitySingleRowWithoutCol(scip, matrix, row, col);
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
            *minresactivity = getMinActivitySingleRowWithoutCol(scip, matrix, row, col);
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

/** calculate the upper bound of one variable from one row */
static
void getVarUpperBoundOfRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index of variable */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this column in this row */
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

   assert(rowub != NULL);
   assert(ubfound != NULL);

   *rowub = SCIPinfinity(scip);
   *ubfound = FALSE;

   getMinMaxActivityResiduals(scip, matrix, col, row, val, &minresactivity, &maxresactivity,
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
   }
   else
   {
      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) )
      {
         *rowub = (lhs - maxresactivity)/val;
         *ubfound = TRUE;
      }
   }
}


/** verify which variable upper bounds are implied */
static
SCIP_Bool isUpperBoundImplied(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col                 /**< column index for implied free test */
   )
{
   SCIP_Real varub;
   SCIP_Real impliedub;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;
   SCIP_Bool ubimplied;

   assert(scip != NULL);
   assert(matrix != NULL);

   ubimplied = FALSE;
   impliedub = SCIPinfinity(scip);

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for( ; (colpnt < colend); colpnt++, valpnt++ )
   {
      SCIP_Real rowub;
      SCIP_Bool ubfound;

      getVarUpperBoundOfRow(scip, matrix, col, *colpnt, *valpnt, &rowub, &ubfound);

      if( ubfound && (rowub < impliedub) )
         impliedub = rowub;
   }

   varub = SCIPmatrixGetColUb(matrix, col);

   if( !SCIPisInfinity(scip, varub) && SCIPisFeasLE(scip, impliedub, varub) )
      ubimplied = TRUE;

   return ubimplied;
}


/** calculate minimal column activity from one continuous variable without one row */
static
SCIP_Real getMinColActWithoutRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   withoutrow,         /**< exclude row index */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   int                   part                /**< which of part of the dual variables is active */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   SCIP_Real mincolactivity;
   int row;
   int p;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);

   mincolactivity = 0;
   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for( ; colpnt < colend; colpnt++, valpnt++ )
   {
      row = *colpnt;
      val = *valpnt;

      if( row == withoutrow )
      {
         if( SCIPmatrixIsRowRhsInfinity(matrix, row) )
            continue;

         /* treat remaining part of equalities or ranged rows */
         if( part == 0 )
         {
            /* consider second part */
            val = -val;
            p = 1;
         }
         else
         {
            /* consider first part */
            p = 0;
         }

         if( val > 0.0 )
         {
            assert(!SCIPisInfinity(scip, -lbdual[p][row]));
            mincolactivity += val * lbdual[p][row];
         }
         else if( val < 0.0 )
         {
            assert(!SCIPisInfinity(scip, ubdual[p][row]));
            mincolactivity += val * ubdual[p][row];
         }
      }
      else
      {
         p = 0;

         do
         {
            if( val > 0.0 )
            {
               assert(!SCIPisInfinity(scip, -lbdual[p][row]));
               mincolactivity += val * lbdual[p][row];
            }
            else if( val < 0.0 )
            {
               assert(!SCIPisInfinity(scip, ubdual[p][row]));
               mincolactivity += val * ubdual[p][row];
            }

            val = -val;
            p++;
         }
         while ( !SCIPmatrixIsRowRhsInfinity(matrix, row) && p < 2 );
      }
   }

   return mincolactivity;
}


/** calculate minimal/maximal column residual activities */
static
void calcColActResidual(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< matrix coefficient */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   int                   part,               /**< which of part of the dual variables is active */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   int*                  mincolactinf,       /**< number of (negative) infinite contributions to minimal column activity */
   SCIP_Real*            mincolresact        /**< minimal column residual activity */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(mincolact != NULL);
   assert(mincolactinf != NULL);
   assert(mincolresact != NULL);

   if( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part == 1 )
      val = -val;

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, -lbdual[part][row]) )
      {
         assert(mincolactinf[col] >= 1);
         if( mincolactinf[col] == 1 )
            *mincolresact = getMinColActWithoutRow(scip, matrix, col, row, lbdual, ubdual, part);
         else
            *mincolresact = -SCIPinfinity(scip);
      }
      else
      {
         if( mincolactinf[col] > 0 )
            *mincolresact = -SCIPinfinity(scip);
         else
            *mincolresact = mincolact[col] - val * lbdual[part][row];
      }
   }
   else if( val < 0.0 )
   {
      if( SCIPisInfinity(scip, ubdual[part][row]) )
      {
         assert(mincolactinf[col] >= 1);
         if( mincolactinf[col] == 1 )
            *mincolresact = getMinColActWithoutRow(scip, matrix, col, row, lbdual, ubdual, part);
         else
            *mincolresact = -SCIPinfinity(scip);
      }
      else
      {
         if( mincolactinf[col] > 0 )
            *mincolresact = -SCIPinfinity(scip);
         else
            *mincolresact = mincolact[col] - val * ubdual[part][row];
      }
   }
}


/** calculate minimal/maximal column activities */
static
void calcColActivity(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Bool             onlyconvars,        /**< consider only continuous variables for activity calculation */
   int                   startcol,           /**< start column for activity calculations */
   int                   stopcol,            /**< stop column for activity calculations */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactinf,       /**< number of (positive) infinite contributions to maximal column activity */
   int*                  mincolactinf,       /**< number of (negative) infinite contributions to minimal column activity */
   SCIP_Bool*            isimplfree          /**< is upper bound of variable implied */
   )
{
   SCIP_VAR* var;
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   int row;
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactinf != NULL);
   assert(mincolactinf != NULL);

   for( c = startcol; c < stopcol; c++ )
   {
      var = SCIPmatrixGetVar(matrix, c);

      if( onlyconvars )
      {
         /* exclude non-continuous variables for bound calculation */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            mincolact[c] = -SCIPinfinity(scip);
            maxcolact[c] = SCIPinfinity(scip);
            maxcolactinf[c] = 1;
            mincolactinf[c] = 1;
            continue;
         }
      }

      /* exclude some bound cases for activity calculation */
      if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) ||
         !(isimplfree[c] || SCIPisInfinity(scip, SCIPvarGetUbGlobal(var))) )
      {
         mincolact[c] = -SCIPinfinity(scip);
         maxcolact[c] = SCIPinfinity(scip);
         maxcolactinf[c] = 1;
         mincolactinf[c] = 1;
         continue;
      }

      /* init activities */
      mincolact[c] = 0;
      maxcolact[c] = 0;
      maxcolactinf[c] = 0;
      mincolactinf[c] = 0;

      colpnt = SCIPmatrixGetColIdxPtr(matrix, c);
      colend = colpnt + SCIPmatrixGetColNNonzs(matrix, c);
      valpnt = SCIPmatrixGetColValPtr(matrix, c);

      /* calculate column activities */
      for( ; colpnt < colend; colpnt++, valpnt++ )
      {
         row = *colpnt;
         val = *valpnt;

         /* consider >= inequalities */
         if( val > 0 )
         {
            if(SCIPisInfinity(scip, ubdual[0][row]))
               maxcolactinf[c]++;
            else
               maxcolact[c] += val * ubdual[0][row];

            if(SCIPisInfinity(scip, -lbdual[0][row]))
               mincolactinf[c]++;
            else
               mincolact[c] += val * lbdual[0][row];
         }
         else if( val < 0.0 )
         {
            if(SCIPisInfinity(scip, -lbdual[0][row]))
               maxcolactinf[c]++;
            else
               maxcolact[c] += val * lbdual[0][row];

            if(SCIPisInfinity(scip, ubdual[0][row]))
               mincolactinf[c]++;
            else
               mincolact[c] += val * ubdual[0][row];
         }

         /* consider equations and ranged rows */
         if( !SCIPmatrixIsRowRhsInfinity(matrix, row) )
         {
            val = -val;
            if( val > 0 )
            {
               if(SCIPisInfinity(scip, ubdual[1][row]))
                  maxcolactinf[c]++;
               else
                  maxcolact[c] += val * ubdual[1][row];

               if(SCIPisInfinity(scip, -lbdual[1][row]))
                  mincolactinf[c]++;
               else
                  mincolact[c] += val * lbdual[1][row];
            }
            else if( val < 0.0 )
            {
               if(SCIPisInfinity(scip, -lbdual[1][row]))
                  maxcolactinf[c]++;
               else
                  maxcolact[c] += val * lbdual[1][row];

               if(SCIPisInfinity(scip, ubdual[1][row]))
                  mincolactinf[c]++;
               else
                  mincolact[c] += val * ubdual[1][row];
            }
         }
      }

      /* update column activities if infinity counters are greater 0 */
      if( mincolactinf[c] > 0 )
         mincolact[c] = -SCIPinfinity(scip);
      if( maxcolactinf[c] > 0 )
         maxcolact[c] = SCIPinfinity(scip);
   }
}


/** update bounds on dual variables */
static
void updateDualBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Real             objval,             /**< objective function value */
   SCIP_Real             val,                /**< matrix coefficient */
   int                   row,                /**< row index */
   SCIP_Real             mincolresact,       /**< minimal column residual activity */
   SCIP_Real*            lbdual[2],          /**< dual lower bounds (ranged, equality) */
   SCIP_Real*            ubdual[2],          /**< dual upper bounds (ranged, equatity)*/
   int                   part,               /**< which of part of the dual variables is active */
   int*                  boundchanges,       /**< number of bound changes */
   SCIP_Bool*            updateinfcnt        /**< flag indicating to update infinity counters */
   )
{
   SCIP_Real newlbdual;
   SCIP_Real newubdual;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(boundchanges != NULL);
   assert(updateinfcnt != NULL);

   *updateinfcnt = FALSE;

   /* return if minimal residual activity is infinite */
   if( SCIPisInfinity(scip, mincolresact) || SCIPisInfinity(scip, -mincolresact) )
      return;

   if( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part == 1 )
      val = -val;

   if( val > 0 )
   {
      newubdual = (objval - mincolresact) / val;
      assert(SCIPisGE(scip, newubdual, 0.0));

      if( newubdual < ubdual[part][row] )
      {
         if( SCIPisInfinity(scip, ubdual[part][row]) )
            *updateinfcnt = TRUE;

         ubdual[part][row] = newubdual;
         (*boundchanges)++;
      }
   }
   else if( val < 0 )
   {
      newlbdual = (objval - mincolresact) / val;
      assert(SCIPisGE(scip, ubdual[part][row], newlbdual));

      if( newlbdual > lbdual[part][row] )
      {
         if( SCIPisInfinity(scip, -lbdual[part][row]) )
            *updateinfcnt = TRUE;

         lbdual[part][row] = newlbdual;
         (*boundchanges)++;
      }
   }
}


/** update minimal/maximal column activity infinity counters */
static
void infCntUpdate(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Real             val,                /**< matrix coefficient */
   int                   row,                /**< row index */
   SCIP_Real*            lbdual[2],          /**< lower bounds of dual variables */
   SCIP_Real*            ubdual[2],          /**< upper bounds of dual variables */
   int                   part,               /**< which part of the dual variables is active */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactinf,       /**< number of (positive) infinity contributions to maximal column activity */
   int*                  mincolactinf,       /**< number of (negative) infinity contributions to minimal column activity */
   SCIP_Bool*            isimplfree          /**< is upper bound of variable implied */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   SCIP_Real aij;
   int c;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual[0] != NULL);
   assert(ubdual[0] != NULL);
   assert(lbdual[1] != NULL);
   assert(ubdual[1] != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactinf != NULL);
   assert(mincolactinf != NULL);

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   if( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part == 1 )
      val = -val;

   /* look at all column entries present within row and update the
    * corresponding infinity counters. if one counter gets to zero,
    * then calculate this column activity new.
    */
   if( val > 0 )
   {
      /* finite upper bound change */
      for(; (rowpnt < rowend); rowpnt++, valpnt++ )
      {
         c = *rowpnt;
         var = SCIPmatrixGetVar(matrix, c);

         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS  ||
            SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)))
            continue;

         aij = *valpnt;

         if( aij < 0 )
         {
            assert(mincolactinf[c] > 0);
            mincolactinf[c]--;

            if( mincolactinf[c] == 0 )
               calcColActivity(scip, matrix, TRUE, c, c+1, lbdual, ubdual,
                  mincolact, maxcolact,
                  maxcolactinf, mincolactinf, isimplfree);
         }
      }
   }
   else if( val < 0 )
   {
      /* finite lower bound change */
      for(; (rowpnt < rowend); rowpnt++, valpnt++ )
      {
         c = *rowpnt;
         var = SCIPmatrixGetVar(matrix, c);

         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS  ||
            SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)))
            continue;

         aij = *valpnt;

         if( aij > 0 )
         {
            assert(mincolactinf[c] > 0);
            mincolactinf[c]--;

            if( mincolactinf[c] == 0 )
               calcColActivity(scip, matrix, TRUE, c, c+1, lbdual, ubdual,
                  mincolact, maxcolact,
                  maxcolactinf, mincolactinf, isimplfree);
         }
      }
   }
}


/** dual bound strengthening */
static
SCIP_RETCODE dualBoundStrengthening(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   FIXINGDIRECTION*      varstofix,          /**< array holding information for later upper/lower bound fixing */
   int*                  npossiblefixings,   /**< number of possible fixings */
   SIDECHANGE*           sidestochange,      /**< array holding if this is an implied equality */
   int*                  npossiblesidechanges/**< number of possible equality changes */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* lbdual[2];
   SCIP_Real* ubdual[2];
   SCIP_Real* mincolact;
   SCIP_Real* maxcolact;
   int* maxcolactinf;
   int* mincolactinf;
   int* colpnt;
   int* colend;
   int boundchanges;
   int loops;
   int c;
   int r;
   int nrows;
   int ncols;
   SCIP_Bool* isubimplied;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(varstofix != NULL);
   assert(npossiblefixings != NULL);
   assert(sidestochange != NULL);
   assert(npossiblesidechanges != NULL);

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &mincolact, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolact, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolactinf, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincolactinf, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(lbdual[0]), nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(ubdual[0]), nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(lbdual[1]), nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(ubdual[1]), nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &isubimplied, ncols) );

   for( c = 0; c < ncols; c++ )
   {
      if( isUpperBoundImplied(scip,matrix,c) )
         isubimplied[c] = TRUE;
      else
         isubimplied[c] = FALSE;
   }

   /* initialize bounds of dual variables */
   for( r = 0; r < nrows; r++ )
   {
      /* >= inequalities */
      lbdual[0][r] = 0.0;
      ubdual[0][r] = SCIPinfinity(scip);

      /* additional dual variable for equalities and ranged rows */
      lbdual[1][r] = 0.0;
      ubdual[1][r] = SCIPinfinity(scip);
   }

   loops = 0;
   boundchanges = 1;

   while( boundchanges && loops < MAX_LOOPS )
   {
      loops++;
      boundchanges = 0;

      calcColActivity(scip, matrix, TRUE, 0, ncols, lbdual, ubdual,
         mincolact, maxcolact, maxcolactinf, mincolactinf, isubimplied);

      for( c = 0 ; c < ncols; c++ )
      {
         SCIP_Real objval;
         SCIP_Real mincolresact;
         SCIP_Bool updateinfcnt;
         SCIP_VAR* var;

         var = SCIPmatrixGetVar(matrix, c);

         /* consider only continuous variables and certain bound cases */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS ||
            SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) ||
            !(isubimplied[c] || SCIPisInfinity(scip, SCIPvarGetUbGlobal(var))))
            continue;

         objval = SCIPvarGetObj(var);
         colpnt = SCIPmatrixGetColIdxPtr(matrix, c);
         colend = colpnt + SCIPmatrixGetColNNonzs(matrix, c);
         valpnt = SCIPmatrixGetColValPtr(matrix, c);
         mincolresact = -SCIPinfinity(scip);

         for( ; colpnt < colend; colpnt++, valpnt++ )
         {
            int row;
            SCIP_Real val;
            int part;

            row = *colpnt;
            val = *valpnt;
            mincolresact = -SCIPinfinity(scip);

            part = 0;

            do
            {
               /* calulate column activity residuals */
               calcColActResidual(scip, matrix, c, row, val, lbdual, ubdual,
                  part, mincolact, mincolactinf, &mincolresact);

               /* update dual bounds */
               updateDualBounds(scip, matrix, objval, val, row, mincolresact,
                  lbdual, ubdual, part, &boundchanges, &updateinfcnt);

               /* update infinity counters if bound changed properly */
               if( updateinfcnt )
                  infCntUpdate(scip, matrix, val, row, lbdual, ubdual, part,
                     mincolact, maxcolact, maxcolactinf, mincolactinf, isubimplied);

               part++;
            }
            while( !SCIPmatrixIsRowRhsInfinity(matrix, row) && part < 2 );
         }
      }
   }

   /* calculate minimal/maximal column activity over all columns */
   calcColActivity(scip, matrix, FALSE, 0, ncols, lbdual, ubdual,
      mincolact, maxcolact, maxcolactinf, mincolactinf, isubimplied);

   for( c = 0; c < ncols; c++ )
   {
      SCIP_Real objval;
      SCIP_VAR* var;

      var = SCIPmatrixGetVar(matrix, c);
      objval = SCIPvarGetObj(var);

      /* positive reduced costs: c_j - sup{(A_{.j})^T}y > 0 => x_j = 0 */
      if( SCIPisLT(scip, maxcolact[c], objval) && varstofix[c] == NOFIX )
      {
         varstofix[c] = FIXATLB;
         (*npossiblefixings)++;
      }
   }

   for( r = 0; r < nrows; r++ )
   {
      /* implied equality: y_i > 0 =>  A_{.i}x - b_i = 0 */
      if( SCIPmatrixIsRowRhsInfinity(matrix, r) )
      {
         if( SCIPisGT(scip, lbdual[0][r], 0.0) )
         {
            /* change >= inequality to equality */
            sidestochange[r] = RHSTOLHS;
            (*npossiblesidechanges)++;
         }
      }
      else
      {
         if( !SCIPmatrixIsRowRhsInfinity(matrix, r) &&
             !SCIPisEQ(scip,SCIPmatrixGetRowLhs(matrix, r),SCIPmatrixGetRowRhs(matrix, r)) )
         {
            /* for ranged rows we have to decide which side determines the equality */
            if( SCIPisGT(scip, lbdual[0][r], 0.0) && !SCIPisGT(scip, lbdual[1][r], 0.0) && sidestochange[r]==NOCHANGE )
            {
               sidestochange[r] = RHSTOLHS;
               (*npossiblesidechanges)++;
            }

            if( !SCIPisGT(scip, lbdual[0][r], 0.0) && SCIPisGT(scip, lbdual[1][r], 0.0) && sidestochange[r]==NOCHANGE)
            {
               sidestochange[r] = LHSTORHS;
               (*npossiblesidechanges)++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &isubimplied);
   SCIPfreeBufferArray(scip, &ubdual[1]);
   SCIPfreeBufferArray(scip, &lbdual[1]);
   SCIPfreeBufferArray(scip, &ubdual[0]);
   SCIPfreeBufferArray(scip, &lbdual[0]);
   SCIPfreeBufferArray(scip, &mincolactinf);
   SCIPfreeBufferArray(scip, &maxcolactinf);
   SCIPfreeBufferArray(scip, &maxcolact);
   SCIPfreeBufferArray(scip, &mincolact);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDualinfer)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDualinfer(scip) );

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualinfer)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) )
      return SCIP_OKAY;

   if( SCIPgetNContVars(scip)==0 )
      return SCIP_OKAY;

   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized && complete )
   {
      FIXINGDIRECTION* varstofix;
      int npossiblefixings;
      int nconvarsfixed;
      int nintvarsfixed;
      int nbinvarsfixed;
      SIDECHANGE* sidestochange;
      int npossiblesidechanges;
      int nsideschanged;
      int i;
      int nrows;
      int ncols;
      SCIP_Bool locksconsistent;
      SCIP_VAR* var;

      npossiblefixings = 0;
      nconvarsfixed = 0;
      nintvarsfixed = 0;
      nbinvarsfixed = 0;
      npossiblesidechanges = 0;
      nsideschanged = 0;
      locksconsistent = TRUE;

      nrows = SCIPmatrixGetNRows(matrix);
      ncols = SCIPmatrixGetNColumns(matrix);

      /* the locks of the continuous variable must be consistent */
      for(i = 0; i < ncols; i++)
      {
         var = SCIPmatrixGetVar(matrix, i);
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            if( SCIPvarGetNLocksUp(var) != SCIPmatrixGetColNUplocks(matrix, i) ||
               SCIPvarGetNLocksDown(var) != SCIPmatrixGetColNDownlocks(matrix, i) )
            {
               locksconsistent = FALSE;
               break;
            }
         }
      }

      if( locksconsistent )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, ncols) );
         SCIP_CALL( SCIPallocBufferArray(scip, &sidestochange, nrows) );

         BMSclearMemoryArray(varstofix, ncols);
         BMSclearMemoryArray(sidestochange, nrows);

         SCIP_CALL( dualBoundStrengthening(scip, matrix, varstofix, &npossiblefixings, sidestochange, &npossiblesidechanges) );

         if( npossiblefixings > 0 )
         {
            for( i = ncols - 1; i >= 0; --i )
            {
               SCIP_Bool infeasible;
               SCIP_Bool fixed;

               if( varstofix[i] == FIXATLB )
               {
                  SCIP_Real lb;

                  var = SCIPmatrixGetVar(matrix, i);

                  if( SCIPvarGetNLocksUp(var) != SCIPmatrixGetColNUplocks(matrix, i) ||
                     SCIPvarGetNLocksDown(var) != SCIPmatrixGetColNDownlocks(matrix, i) )
                  {
                     /* no fixing, locks for this variable not consistent */
                     continue;
                  }

                  lb = SCIPvarGetLbLocal(var);

                  /* fix at lower bound */
                  SCIP_CALL( SCIPfixVar(scip, var, lb, &infeasible, &fixed) );
                  if( infeasible )
                  {
                     SCIPdebugMsg(scip, " -> infeasible fixing\n");
                     *result = SCIP_CUTOFF;
                     break;
                  }
                  assert(fixed);
                  (*nfixedvars)++;
                  *result = SCIP_SUCCESS;

                  if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                     nconvarsfixed++;
                  else if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
                     nbinvarsfixed++;
                  else
                     nintvarsfixed++;
               }
            }
         }

         if( npossiblesidechanges > 0 )
         {
            for( i = 0; i < nrows; i++ )
            {
               SCIP_CONS* cons;
               SCIP_CONSHDLR* conshdlr;
               const char* conshdlrname;

               if( sidestochange[i] == NOCHANGE )
                  continue;

               cons = SCIPmatrixGetCons(matrix,i);
               conshdlr = SCIPconsGetHdlr(cons);
               conshdlrname = SCIPconshdlrGetName(conshdlr);

               /* TODO: implement further constraint types */
               if( strcmp(conshdlrname, "linear") == 0 )
               {
                  SCIP_Real lhs;
                  SCIP_Real rhs;
                  SCIP_Real matrixlhs;
                  SCIP_Real matrixrhs;

                  lhs = SCIPgetLhsLinear(scip, cons);
                  rhs = SCIPgetRhsLinear(scip, cons);
                  matrixlhs = SCIPmatrixGetRowLhs(matrix, i);
                  matrixrhs = SCIPmatrixGetRowRhs(matrix, i);

                  if( sidestochange[i] == RHSTOLHS )
                  {
                     assert(!SCIPisEQ(scip, matrixlhs, matrixrhs));

                     if( SCIPisEQ(scip, matrixlhs, lhs) )
                        SCIP_CALL( SCIPchgRhsLinear(scip, cons, matrixlhs) );
                     else
                        SCIP_CALL( SCIPchgLhsLinear(scip, cons, -matrixlhs) );

                     nsideschanged++;
                     (*nchgsides)++;
                  }
                  else
                  {
                     assert(!SCIPisEQ(scip, matrixlhs, matrixrhs));

                     if( SCIPisEQ(scip, matrixrhs, rhs) )
                        SCIP_CALL( SCIPchgLhsLinear(scip, cons, matrixrhs) );
                     else
                        SCIP_CALL( SCIPchgRhsLinear(scip, cons, -matrixrhs) );

                     nsideschanged++;
                     (*nchgsides)++;
                  }
               }
            }
         }

         SCIPfreeBufferArray(scip, &sidestochange);
         SCIPfreeBufferArray(scip, &varstofix);

         if( (nconvarsfixed + nintvarsfixed + nbinvarsfixed) > 0 || npossiblesidechanges > 0)
         {
            SCIPdebugMsg(scip, "### fixed vars [cont: %d, int: %d, bin: %d], changed sides [%d]\n",
               nconvarsfixed, nintvarsfixed, nbinvarsfixed, nsideschanged);
         }
      }
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the dual inference presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualinfer(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecDualinfer, NULL) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDualinfer) );

   return SCIP_OKAY;
}
