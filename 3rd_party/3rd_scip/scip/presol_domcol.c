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

/**@file    presol_domcol.c
 * @ingroup PRESOLVERS
 * @brief   dominated column presolver
 * @author  Dieter Weninger
 * @author  Gerald Gamrath
 *
 * This presolver looks for dominance relations between variable pairs.
 * From a dominance relation and certain bound/clique-constellations
 * variable fixings mostly at the lower bound of the dominated variable can be derived.
 * Additionally it is possible to improve bounds by predictive bound strengthening.
 *
 * @todo Also run on general CIPs, if the number of locks of the investigated variables comes only from (upgraded)
 * linear constraints.
 *
 * @todo Instead of choosing variables for comparison out of one row, we should try to use 'hashvalues' for columns that
 *       indicate in which constraint type (<=, >=, or ranged row / ==) they are existing. Then sort the variables (and
 *       corresponding data) after the ranged row/equation hashvalue and only try to derive dominance on variables with
 *       the same hashvalue on ranged row/equation and fitting hashvalues for the other constraint types.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/pub_matrix.h"
#include "presol_domcol.h"

#define PRESOL_NAME            "domcol"
#define PRESOL_DESC            "dominated column presolver"
#define PRESOL_PRIORITY            -1000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_NUMMINPAIRS         1024     /**< minimal number of pair comparisons */
#define DEFAULT_NUMMAXPAIRS      1048576     /**< maximal number of pair comparisons */

#define DEFAULT_PREDBNDSTR         FALSE     /**< should predictive bound strengthening be applied? */
#define DEFAULT_CONTINUOUS_RED      TRUE     /**< should reductions for continuous variables be carried out? */



/*
 * Data structures
  */

/** control parameters */
struct SCIP_PresolData
{
   int                   numminpairs;        /**< minimal number of pair comparisons */
   int                   nummaxpairs;        /**< maximal number of pair comparisons */
   int                   numcurrentpairs;    /**< current number of pair comparisons */
   SCIP_Bool             predbndstr;         /**< flag indicating if predictive bound strengthening should be applied */
   SCIP_Bool             continuousred;      /**< flag indicating if reductions for continuous variables should be performed */
};

/** type of fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,         /**< fix variable at lower bound */
   NOFIX   =  0,         /**< do not fix variable */
   FIXATUB =  1          /**< fix variable at upper bound */
};
typedef enum Fixingdirection FIXINGDIRECTION;


/*
 * Local methods
 */
#ifdef SCIP_DEBUG
/** print a row from the constraint matrix */
static
void printRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row                 /**< row index for printing */
   )
{
   int* rowpnt;
   int* rowend;
   int c;
   SCIP_Real val;
   SCIP_Real* valpnt;
   char relation;
   SCIP_VAR* var;

   relation='-';
   if( !SCIPmatrixIsRowRhsInfinity(matrix, row) &&
      SCIPisEQ(scip, SCIPmatrixGetRowLhs(matrix, row), SCIPmatrixGetRowRhs(matrix, row)) )
   {
      relation='e';
   }
   else if( !SCIPmatrixIsRowRhsInfinity(matrix, row) &&
      !SCIPisEQ(scip, SCIPmatrixGetRowLhs(matrix, row), SCIPmatrixGetRowRhs(matrix, row)) )
   {
      relation='r';
   }
   else
   {
      relation='g';
   }

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   SCIPdebugMsgPrint(scip, "\n(L:%g) [%c] %g  <=",
      (SCIPmatrixGetRowNMinActPosInf(matrix, row) + SCIPmatrixGetRowNMinActNegInf(matrix,row) > 0) ?
      -SCIPinfinity(scip) :
      SCIPmatrixGetRowMinActivity(matrix, row), relation, SCIPmatrixGetRowLhs(matrix, row));
   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      c = *rowpnt;
      val = *valpnt;
      var = SCIPmatrixGetVar(matrix, c);
      SCIPdebugMsgPrint(scip, "  %g{%s[idx:%d][bnd:%g,%g]}",
         val, SCIPvarGetName(var), c, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }
   SCIPdebugMsgPrint(scip, " <=  %g (U:%g)", (SCIPmatrixGetRowNMaxActPosInf(matrix, row) + SCIPmatrixGetRowNMaxActNegInf(matrix, row) > 0) ?
      SCIPinfinity(scip) :
      SCIPmatrixGetRowRhs(matrix, row), SCIPmatrixGetRowMaxActivity(matrix , row));
}

/** print all rows from the constraint matrix containing a variable */
static
SCIP_RETCODE printRowsOfCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col                 /**< column index for printing */
   )
{
   int numrows;
   int* rows;
   int i;
   int* colpnt;
   int* colend;

   numrows = SCIPmatrixGetColNNonzs(matrix, col);

   SCIP_CALL( SCIPallocBufferArray(scip, &rows, numrows) );

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   for( i = 0; (colpnt < colend); colpnt++, i++ )
   {
      rows[i] = *colpnt;
   }

   SCIPdebugMsgPrint(scip, "\n-------");
   SCIPdebugMsgPrint(scip, "\ncol %d number rows: %d",col,numrows);
   for( i = 0; i < numrows; i++ )
   {
      printRow(scip, matrix, rows[i]);
   }
   SCIPdebugMsgPrint(scip, "\n-------");

   SCIPfreeBufferArray(scip, &rows);

   return SCIP_OKAY;
}

/** print information about a dominance relation */
static
SCIP_RETCODE printDomRelInfo(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_VAR*             dominatingvar,      /**< dominating variable */
   int                   dominatingidx,      /**< index of dominating variable */
   SCIP_VAR*             dominatedvar,       /**< dominated variable */
   int                   dominatedidx,       /**< index of dominated variable */
   SCIP_Real             dominatingub,       /**< predicted upper bound of dominating variable */
   SCIP_Real             dominatingwclb      /**< worst case lower bound of dominating variable */
   )
{
   char type;

   assert(SCIPvarGetType(dominatingvar)==SCIPvarGetType(dominatedvar));

   switch(SCIPvarGetType(dominatingvar))
   {
   case SCIP_VARTYPE_CONTINUOUS:
      type='C';
      break;
   case SCIP_VARTYPE_BINARY:
      type='B';
      break;
   case SCIP_VARTYPE_INTEGER:
   case SCIP_VARTYPE_IMPLINT:
      type='I';
      break;
   default:
      SCIPerrorMessage("unknown variable type\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   SCIPdebugMsgPrint(scip, "\n\n### [%c], obj:%g->%g,\t%s[idx:%d](nrows:%d)->%s[idx:%d](nrows:%d)\twclb=%g, ub'=%g, ub=%g",
      type, SCIPvarGetObj(dominatingvar), SCIPvarGetObj(dominatedvar),
      SCIPvarGetName(dominatingvar), dominatingidx, SCIPmatrixGetColNNonzs(matrix, dominatingidx),
      SCIPvarGetName(dominatedvar), dominatedidx, SCIPmatrixGetColNNonzs(matrix, dominatedidx),
      dominatingwclb, dominatingub, SCIPvarGetUbGlobal(dominatingvar));

   SCIP_CALL( printRowsOfCol(scip, matrix, dominatingidx) );

   return SCIP_OKAY;
}
#endif


/** get minimum/maximum residual activity for the specified variable and with another variable set to its upper bound */
static
void getActivityResidualsUpperBound(
   SCIP*                 scip,
   SCIP_MATRIX*          matrix,
   int                   row,
   int                   col,
   SCIP_Real             coef,
   int                   upperboundcol,
   SCIP_Real             upperboundcoef,
   SCIP_Real*            minresactivity,
   SCIP_Real*            maxresactivity,
   SCIP_Bool*            success
   )
{
   SCIP_VAR* var;
   SCIP_VAR* upperboundvar;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real lbupperboundvar;
   SCIP_Real ubupperboundvar;
   SCIP_Real tmpmaxact;
   SCIP_Real tmpminact;
   int maxactinf;
   int minactinf;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(row < SCIPmatrixGetNRows(matrix));
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);

   var = SCIPmatrixGetVar(matrix, col);
   upperboundvar = SCIPmatrixGetVar(matrix, upperboundcol);
   assert(var != NULL);
   assert(upperboundvar != NULL);

   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);

   maxactinf = SCIPmatrixGetRowNMaxActPosInf(matrix, row) + SCIPmatrixGetRowNMaxActNegInf(matrix, row);
   minactinf = SCIPmatrixGetRowNMinActPosInf(matrix, row) + SCIPmatrixGetRowNMinActNegInf(matrix, row);

   assert(!SCIPisInfinity(scip, lb));
   assert(!SCIPisInfinity(scip, -ub));

   lbupperboundvar = SCIPvarGetLbGlobal(upperboundvar);
   ubupperboundvar = SCIPvarGetUbGlobal(upperboundvar);

   assert(!SCIPisInfinity(scip, lbupperboundvar));
   assert(!SCIPisInfinity(scip, -ubupperboundvar));

   tmpminact = SCIPmatrixGetRowMinActivity(matrix, row);
   tmpmaxact = SCIPmatrixGetRowMaxActivity(matrix, row);

   /* first, adjust min and max activity such that upperboundvar is always set to its upper bound */

   /* abort if the upper bound is infty */
   if( SCIPisInfinity(scip, ubupperboundvar) )
   {
      *success = FALSE;
      return;
   }
   else
   {
      /* coef > 0 --> the lower bound is used for the minactivity --> update the minactivity */
      if( upperboundcoef > 0.0 )
      {
         if( SCIPisInfinity(scip, -lbupperboundvar) )
         {
            assert(minactinf > 0);
            minactinf--;
            tmpminact += (upperboundcoef * ubupperboundvar);
         }
         else
         {
            tmpminact = tmpminact - (upperboundcoef * lbupperboundvar) + (upperboundcoef * ubupperboundvar);
         }
      }
      /* coef < 0 --> the lower bound is used for the maxactivity --> update the maxactivity */
      else
      {
         if( SCIPisInfinity(scip, -lbupperboundvar) )
         {
            assert(maxactinf > 0);
            maxactinf--;
            tmpmaxact += (upperboundcoef * ubupperboundvar);
         }
         else
         {
            tmpmaxact = tmpmaxact - (upperboundcoef * lbupperboundvar) + (upperboundcoef * ubupperboundvar);
         }
      }
   }

   *success = TRUE;

   /* the coefficient is positive, so the upper bound contributed to the maxactivity and the lower bound to the minactivity */
   if( coef >= 0.0 )
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(maxactinf >= 1);
         if( maxactinf == 1 )
         {
            *maxresactivity = tmpmaxact;
         }
         else
            *maxresactivity = SCIPinfinity(scip);
      }
      else
      {
         if( maxactinf > 0 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = tmpmaxact - coef * ub;
      }

      if( SCIPisInfinity(scip, -lb) )
      {
         assert(minactinf >= 1);
         if( minactinf == 1 )
         {
            *minresactivity = tmpminact;
         }
         else
            *minresactivity = -SCIPinfinity(scip);
      }
      else
      {
         if( minactinf > 0 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = tmpminact - coef * lb;
      }
   }
   /* the coefficient is negative, so the lower bound contributed to the maxactivity and the upper bound to the minactivity */
   else
   {
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(maxactinf >= 1);
         if( maxactinf == 1 )
         {
            *maxresactivity = tmpmaxact;
         }
         else
            *maxresactivity = SCIPinfinity(scip);
      }
      else
      {
         if( maxactinf > 0 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = tmpmaxact - coef * lb;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         assert(minactinf >= 1);
         if( minactinf == 1 )
         {
            *minresactivity = tmpminact;
         }
         else
            *minresactivity = -SCIPinfinity(scip);
      }
      else
      {
         if( minactinf > 0 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = tmpminact - coef * ub;
      }
   }
}

/** get minimum/maximum residual activity for the specified variable and with another variable set to its lower bound */
static
void getActivityResidualsLowerBound(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col,                /**< column index */
   SCIP_Real             coef,               /**< coefficient of the column in this row */
   int                   lowerboundcol,      /**< column index of variable to set to its lower bound */
   SCIP_Real             lowerboundcoef,     /**< coefficient in this row of the column to be set to its lower bound */
   SCIP_Real*            minresactivity,     /**< minimum residual activity of this row */
   SCIP_Real*            maxresactivity,     /**< maximum residual activity of this row */
   SCIP_Bool*            success             /**< pointer to store whether the computation was successful */
   )
{
   SCIP_VAR* var;
   SCIP_VAR* lowerboundvar;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real lblowerboundvar;
   SCIP_Real ublowerboundvar;
   SCIP_Real tmpmaxact;
   SCIP_Real tmpminact;
   int maxactinf;
   int minactinf;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(row < SCIPmatrixGetNRows(matrix));
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);

   var = SCIPmatrixGetVar(matrix, col);;
   lowerboundvar = SCIPmatrixGetVar(matrix, lowerboundcol);
   assert(var != NULL);
   assert(lowerboundvar != NULL);

   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);

   maxactinf = SCIPmatrixGetRowNMaxActPosInf(matrix, row) + SCIPmatrixGetRowNMaxActNegInf(matrix, row);
   minactinf = SCIPmatrixGetRowNMinActPosInf(matrix, row) + SCIPmatrixGetRowNMinActNegInf(matrix, row);

   assert(!SCIPisInfinity(scip, lb));
   assert(!SCIPisInfinity(scip, -ub));

   lblowerboundvar = SCIPvarGetLbGlobal(lowerboundvar);
   ublowerboundvar = SCIPvarGetUbGlobal(lowerboundvar);

   assert(!SCIPisInfinity(scip, lblowerboundvar));
   assert(!SCIPisInfinity(scip, -ublowerboundvar));

   tmpminact = SCIPmatrixGetRowMinActivity(matrix, row);
   tmpmaxact = SCIPmatrixGetRowMaxActivity(matrix, row);

   /* first, adjust min and max activity such that lowerboundvar is always set to its lower bound */

   /* abort if the lower bound is -infty */
   if( SCIPisInfinity(scip, -lblowerboundvar) )
   {
      *success = FALSE;
      return;
   }
   else
   {
      /* coef > 0 --> the upper bound is used for the maxactivity --> update the maxactivity */
      if( lowerboundcoef > 0.0 )
      {
         if( SCIPisInfinity(scip, ublowerboundvar) )
         {
            assert(maxactinf > 0);
            maxactinf--;
            tmpmaxact += (lowerboundcoef * lblowerboundvar);
         }
         else
         {
            tmpmaxact = tmpmaxact - (lowerboundcoef * ublowerboundvar) + (lowerboundcoef * lblowerboundvar);
         }
      }
      /* coef < 0 --> the upper bound is used for the minactivity --> update the minactivity */
      else
      {
         if( SCIPisInfinity(scip, ublowerboundvar) )
         {
            assert(minactinf > 0);
            minactinf--;
            tmpminact += (lowerboundcoef * lblowerboundvar);
         }
         else
         {
            tmpminact = tmpminact - (lowerboundcoef * ublowerboundvar) + (lowerboundcoef * lblowerboundvar);
         }
      }
   }

   *success = TRUE;

   /* the coefficient is positive, so the upper bound contributed to the maxactivity and the lower bound to the minactivity */
   if( coef >= 0.0 )
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(maxactinf >= 1);
         if( maxactinf == 1 )
         {
            *maxresactivity = tmpmaxact;
         }
         else
            *maxresactivity = SCIPinfinity(scip);
      }
      else
      {
         if( maxactinf > 0 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = tmpmaxact - coef * ub;
      }

      if( SCIPisInfinity(scip, -lb) )
      {
         assert(minactinf >= 1);
         if( minactinf == 1 )
         {
            *minresactivity = tmpminact;
         }
         else
            *minresactivity = -SCIPinfinity(scip);
      }
      else
      {
         if( minactinf > 0 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = tmpminact - coef * lb;
      }
   }
   /* the coefficient is negative, so the lower bound contributed to the maxactivity and the upper bound to the minactivity */
   else
   {
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(maxactinf >= 1);
         if( maxactinf == 1 )
         {
            *maxresactivity = tmpmaxact;
         }
         else
            *maxresactivity = SCIPinfinity(scip);
      }
      else
      {
         if( maxactinf > 0 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = tmpmaxact - coef * lb;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         assert(minactinf >= 1);
         if( minactinf == 1 )
         {
            *minresactivity = tmpminact;
         }
         else
            *minresactivity = -SCIPinfinity(scip);
      }
      else
      {
         if( minactinf > 0 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = tmpminact - coef * ub;
      }
   }
}

/** Calculate bounds of the dominated variable by rowbound analysis.
 *  We use a special kind of predictive rowbound analysis by first setting the dominating variable to its upper bound.
 */
static
SCIP_RETCODE calcVarBoundsDominated(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< current row index */
   int                   coldominating,      /**< column index of dominating variable */
   SCIP_Real             valdominating,      /**< row coefficient of dominating variable */
   int                   coldominated,       /**< column index of dominated variable */
   SCIP_Real             valdominated,       /**< row coefficient of dominated variable */
   SCIP_Bool*            ubcalculated,       /**< was a upper bound calculated? */
   SCIP_Real*            calculatedub,       /**< predicted upper bound */
   SCIP_Bool*            wclbcalculated,     /**< was a lower worst case bound calculated? */
   SCIP_Real*            calculatedwclb,     /**< predicted worst case lower bound */
   SCIP_Bool*            lbcalculated,       /**< was a lower bound calculated? */
   SCIP_Real*            calculatedlb,       /**< predicted lower bound */
   SCIP_Bool*            wcubcalculated,     /**< was a worst case upper bound calculated? */
   SCIP_Real*            calculatedwcub      /**< calculated worst case upper bound */
   )
{
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < SCIPmatrixGetNRows(matrix));
   assert(0 <= coldominating && coldominating < SCIPmatrixGetNColumns(matrix));
   assert(0 <= coldominated && coldominated < SCIPmatrixGetNColumns(matrix));

   assert(ubcalculated != NULL);
   assert(calculatedub != NULL);
   assert(wclbcalculated != NULL);
   assert(calculatedwclb != NULL);
   assert(lbcalculated != NULL);
   assert(calculatedlb != NULL);
   assert(wcubcalculated != NULL);
   assert(calculatedwcub != NULL);

   assert(!SCIPisZero(scip, valdominated));
   assert(SCIPmatrixGetVar(matrix, coldominated) != NULL);

   *ubcalculated = FALSE;
   *wclbcalculated = FALSE;
   *lbcalculated = FALSE;
   *wcubcalculated = FALSE;

   /* no rowbound analysis for multiaggregated variables, which should not exist, because the matrix only consists of
    * active variables
    */
   assert(SCIPvarGetStatus(SCIPmatrixGetVar(matrix, coldominating)) != SCIP_VARSTATUS_MULTAGGR);
   assert(SCIPvarGetStatus(SCIPmatrixGetVar(matrix, coldominated)) != SCIP_VARSTATUS_MULTAGGR);

   lhs = SCIPmatrixGetRowLhs(matrix, row);
   rhs = SCIPmatrixGetRowRhs(matrix, row);
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));

   /* get minimum/maximum activity of this row without the dominated variable */
   getActivityResidualsUpperBound(scip, matrix, row, coldominated, valdominated, coldominating, valdominating,
      &minresactivity, &maxresactivity, &success);

   if( !success )
      return SCIP_OKAY;

   assert(!SCIPisInfinity(scip, minresactivity));
   assert(!SCIPisInfinity(scip, -maxresactivity));

   *calculatedub = SCIPinfinity(scip);
   *calculatedwclb = -SCIPinfinity(scip);
   *calculatedlb = -SCIPinfinity(scip);
   *calculatedwcub = SCIPinfinity(scip);

   /* predictive rowbound analysis */

   if( valdominated > 0.0 )
   {
      /* lower bound calculation */
      if( !SCIPisInfinity(scip, maxresactivity) )
      {
         *calculatedlb = (lhs - maxresactivity)/valdominated;
         *lbcalculated = TRUE;
      }

      /* worst case calculation of lower bound */
      if( !SCIPisInfinity(scip, -minresactivity) )
      {
         *calculatedwclb = (lhs - minresactivity)/valdominated;
         *wclbcalculated = TRUE;
      }
      else
      {
         /* worst case lower bound is infinity */
         *calculatedwclb = SCIPinfinity(scip);
         *wclbcalculated = TRUE;
      }

      /* we can only calculate upper bounds, of the right hand side is finite */
      if( !SCIPmatrixIsRowRhsInfinity(matrix, row) )
      {
         /* upper bound calculation */
         if( !SCIPisInfinity(scip, -minresactivity) )
         {
            *calculatedub = (rhs - minresactivity)/valdominated;
            *ubcalculated = TRUE;
         }

         /* worst case calculation of upper bound */
         if( !SCIPisInfinity(scip, maxresactivity) )
         {
            *calculatedwcub = (rhs - maxresactivity)/valdominated;
            *wcubcalculated = TRUE;
         }
         else
         {
            /* worst case upper bound is -infinity */
            *calculatedwcub = -SCIPinfinity(scip);
            *wcubcalculated = TRUE;
         }
      }
   }
   else
   {
      /* upper bound calculation */
      if( !SCIPisInfinity(scip, maxresactivity) )
      {
         *calculatedub = (lhs - maxresactivity)/valdominated;
         *ubcalculated = TRUE;
      }

      /* worst case calculation of upper bound */
      if( !SCIPisInfinity(scip, -minresactivity) )
      {
         *calculatedwcub = (lhs - minresactivity)/valdominated;
         *wcubcalculated = TRUE;
      }
      else
      {
         /* worst case upper bound is -infinity */
         *calculatedwcub = -SCIPinfinity(scip);
         *wcubcalculated = TRUE;
      }

      /* we can only calculate lower bounds, of the right hand side is finite */
      if( !SCIPmatrixIsRowRhsInfinity(matrix, row) )
      {
         /* lower bound calculation */
         if( !SCIPisInfinity(scip, -minresactivity) )
         {
            *calculatedlb = (rhs - minresactivity)/valdominated;
            *lbcalculated = TRUE;
         }

         /* worst case calculation of lower bound */
         if( !SCIPisInfinity(scip, maxresactivity) )
         {
            *calculatedwclb = (rhs - maxresactivity)/valdominated;
            *wclbcalculated = TRUE;
         }
         else
         {
            /* worst case lower bound is infinity */
            *calculatedwclb = SCIPinfinity(scip);
            *wclbcalculated = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/** Calculate bounds of the dominating variable by rowbound analysis.
 *  We use a special kind of predictive rowbound analysis by first setting the dominated variable to its lower bound.
 */
static
SCIP_RETCODE calcVarBoundsDominating(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< current row index */
   int                   coldominating,      /**< column index of dominating variable */
   SCIP_Real             valdominating,      /**< row coefficient of dominating variable */
   int                   coldominated,       /**< column index of dominated variable */
   SCIP_Real             valdominated,       /**< row coefficient of dominated variable */
   SCIP_Bool*            ubcalculated,       /**< was a upper bound calculated? */
   SCIP_Real*            calculatedub,       /**< predicted upper bound */
   SCIP_Bool*            wclbcalculated,     /**< was a lower worst case bound calculated? */
   SCIP_Real*            calculatedwclb,     /**< predicted worst case lower bound */
   SCIP_Bool*            lbcalculated,       /**< was a lower bound calculated? */
   SCIP_Real*            calculatedlb,       /**< predicted lower bound */
   SCIP_Bool*            wcubcalculated,     /**< was a worst case upper bound calculated? */
   SCIP_Real*            calculatedwcub      /**< calculated worst case upper bound */
   )
{
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < SCIPmatrixGetNRows(matrix) );
   assert(0 <= coldominating && coldominating < SCIPmatrixGetNColumns(matrix));
   assert(0 <= coldominated && coldominated < SCIPmatrixGetNColumns(matrix));

   assert(ubcalculated != NULL);
   assert(calculatedub != NULL);
   assert(wclbcalculated != NULL);
   assert(calculatedwclb != NULL);
   assert(lbcalculated != NULL);
   assert(calculatedlb != NULL);
   assert(wcubcalculated != NULL);
   assert(calculatedwcub != NULL);

   assert(!SCIPisZero(scip, valdominating));
   assert(SCIPmatrixGetVar(matrix, coldominating) != NULL);

   *ubcalculated = FALSE;
   *wclbcalculated = FALSE;
   *lbcalculated = FALSE;
   *wcubcalculated = FALSE;

   /* no rowbound analysis for multiaggregated variables, which should not exist, because the matrix only consists of
    * active variables
    */
   assert(SCIPvarGetStatus(SCIPmatrixGetVar(matrix, coldominating)) != SCIP_VARSTATUS_MULTAGGR);
   assert(SCIPvarGetStatus(SCIPmatrixGetVar(matrix, coldominated)) != SCIP_VARSTATUS_MULTAGGR);

   lhs = SCIPmatrixGetRowLhs(matrix, row);
   rhs = SCIPmatrixGetRowRhs(matrix, row);
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));

   /* get minimum/maximum activity of this row without the dominating variable */
   getActivityResidualsLowerBound(scip, matrix, row, coldominating, valdominating, coldominated, valdominated,
      &minresactivity, &maxresactivity, &success);

   if( !success )
      return SCIP_OKAY;

   assert(!SCIPisInfinity(scip, minresactivity));
   assert(!SCIPisInfinity(scip, -maxresactivity));

   *calculatedub = SCIPinfinity(scip);
   *calculatedwclb = -SCIPinfinity(scip);
   *calculatedlb = -SCIPinfinity(scip);
   *calculatedwcub = SCIPinfinity(scip);

   /* predictive rowbound analysis */

   if( valdominating > 0.0 )
   {
      /* lower bound calculation */
      if( !SCIPisInfinity(scip, maxresactivity) )
      {
         *calculatedlb = (lhs - maxresactivity)/valdominating;
         *lbcalculated = TRUE;
      }

      /* worst case calculation of lower bound */
      if( !SCIPisInfinity(scip, -minresactivity) )
      {
         *calculatedwclb = (lhs - minresactivity)/valdominating;
         *wclbcalculated = TRUE;
      }
      else
      {
         /* worst case lower bound is infinity */
         *calculatedwclb = SCIPinfinity(scip);
         *wclbcalculated = TRUE;
      }

      /* we can only calculate upper bounds, of the right hand side is finite */
      if( !SCIPmatrixIsRowRhsInfinity(matrix, row) )
      {
         /* upper bound calculation */
         if( !SCIPisInfinity(scip, -minresactivity) )
         {
            *calculatedub = (rhs - minresactivity)/valdominating;
            *ubcalculated = TRUE;
         }

         /* worst case calculation of upper bound */
         if( !SCIPisInfinity(scip, maxresactivity) )
         {
            *calculatedwcub = (rhs - maxresactivity)/valdominating;
            *wcubcalculated = TRUE;
         }
         else
         {
            /* worst case upper bound is -infinity */
            *calculatedwcub = -SCIPinfinity(scip);
            *wcubcalculated = TRUE;
         }
      }
   }
   else
   {
      /* upper bound calculation */
      if( !SCIPisInfinity(scip, maxresactivity) )
      {
         *calculatedub = (lhs - maxresactivity)/valdominating;
         *ubcalculated = TRUE;
      }

      /* worst case calculation of upper bound */
      if( !SCIPisInfinity(scip, -minresactivity) )
      {
         *calculatedwcub = (lhs - minresactivity)/valdominating;
         *wcubcalculated = TRUE;
      }
      else
      {
         /* worst case upper bound is -infinity */
         *calculatedwcub = -SCIPinfinity(scip);
         *wcubcalculated = TRUE;
      }

      /* we can only calculate lower bounds, of the right hand side is finite */
      if( !SCIPmatrixIsRowRhsInfinity(matrix, row) )
      {
         /* lower bound calculation */
         if( !SCIPisInfinity(scip, -minresactivity) )
         {
            *calculatedlb = (rhs - minresactivity)/valdominating;
            *lbcalculated = TRUE;
         }

         /* worst case calculation of lower bound */
         if( !SCIPisInfinity(scip, maxresactivity) )
         {
            *calculatedwclb = (rhs - maxresactivity)/valdominating;
            *wclbcalculated = TRUE;
         }
         else
         {
            /* worst case lower bound is infinity */
            *calculatedwclb = SCIPinfinity(scip);
            *wclbcalculated = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/** try to find new variable bounds and update them when they are better then the old bounds */
static
SCIP_RETCODE updateBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col1,               /**< dominating variable index */
   SCIP_Real             val1,               /**< dominating variable coefficient */
   int                   col2,               /**< dominated variable index */
   SCIP_Real             val2,               /**< dominated variable coefficient */
   SCIP_Bool             predictdominating,  /**< flag indicating if bounds of dominating or dominated variable are predicted */
   SCIP_Real*            upperbound,         /**< predicted upper bound */
   SCIP_Real*            wclowerbound,       /**< predicted worst case lower bound */
   SCIP_Real*            lowerbound,         /**< predicted lower bound */
   SCIP_Real*            wcupperbound        /**< predicted worst case upper bound */
   )
{
   SCIP_Bool ubcalculated;
   SCIP_Bool wclbcalculated;
   SCIP_Bool lbcalculated;
   SCIP_Bool wcubcalculated;
   SCIP_Real newub;
   SCIP_Real newwclb;
   SCIP_Real newlb;
   SCIP_Real newwcub;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(upperbound != NULL);
   assert(wclowerbound != NULL);
   assert(lowerbound != NULL);
   assert(wcupperbound != NULL);

   newub = SCIPinfinity(scip);
   newlb = -SCIPinfinity(scip);
   newwcub = newub;
   newwclb = newlb;

   if( predictdominating )
   {
      /* do predictive rowbound analysis for the dominating variable */
      SCIP_CALL( calcVarBoundsDominating(scip, matrix, row, col1, val1, col2, val2,
            &ubcalculated, &newub, &wclbcalculated, &newwclb,
            &lbcalculated, &newlb, &wcubcalculated, &newwcub) );
   }
   else
   {
      /* do predictive rowbound analysis for the dominated variable */
      SCIP_CALL( calcVarBoundsDominated(scip, matrix, row, col1, val1, col2, val2,
            &ubcalculated, &newub, &wclbcalculated, &newwclb,
            &lbcalculated, &newlb, &wcubcalculated, &newwcub) );
   }

   /* update bounds in case if they are better */
   if( ubcalculated )
   {
      if( newub < *upperbound )
         *upperbound = newub;
   }
   if( wclbcalculated )
   {
      if( newwclb > *wclowerbound )
         *wclowerbound = newwclb;
   }
   if( lbcalculated )
   {
      if( newlb > *lowerbound )
         *lowerbound = newlb;
   }
   if( wcubcalculated )
   {
      if( newwcub < *wcupperbound )
         *wcupperbound = newwcub;
   }

   return SCIP_OKAY;
}

/** detect parallel columns by using the algorithm of Bixby and Wagner
 *  see paper: "A note on Detecting Simple Redundancies in Linear Systems", June 1986
 */
static
SCIP_RETCODE detectParallelCols(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int*                  pclass,             /**< parallel column classes */
   SCIP_Bool*            varineq             /**< indicating if variable is within an equation */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* values;
   SCIP_Real* scale;
   int* classsizes;
   int* pcset;
   int* rowpnt;
   int* rowend;
   int* colindices;
   int* pcs;
   SCIP_Real startval;
   SCIP_Real aij;
   int startpc;
   int startk;
   int startt;
   int pcsetfill;
   int colidx;
   int k;
   int t;
   int m;
   int i;
   int r;
   int newpclass;
   int pc;
   int nrows;
   int ncols;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(pclass != NULL);
   assert(varineq != NULL);

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &classsizes, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scale, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcset, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &values, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colindices, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcs, ncols) );

   BMSclearMemoryArray(scale, ncols);
   BMSclearMemoryArray(pclass, ncols);
   BMSclearMemoryArray(classsizes, ncols);

   classsizes[0] = ncols;
   pcsetfill = 0;
   for( t = 1; t < ncols; ++t )
      pcset[pcsetfill++] = t;

   /* loop over all rows */
   for( r = 0; r < nrows; ++r )
   {
      /* we consider only non-empty equations or ranged rows */
      if( !SCIPmatrixIsRowRhsInfinity(matrix, r) && SCIPmatrixGetRowNNonzs(matrix, r) > 0 )
      {
         rowpnt = SCIPmatrixGetRowIdxPtr(matrix, r);
         rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, r);
         valpnt = SCIPmatrixGetRowValPtr(matrix, r);

         i = 0;
         for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
         {
            aij = *valpnt;
            colidx = *rowpnt;

            /* remember variable was part of an equation or ranged row */
            varineq[colidx] = TRUE;

            if( scale[colidx] == 0.0 )
               scale[colidx] = aij;
            assert(scale[colidx] != 0.0);

            colindices[i] = colidx;
            values[i] = aij / scale[colidx];
            pc = pclass[colidx];
            assert(pc < ncols);

            /* update class sizes and pclass set */
            assert(classsizes[pc] > 0);
            classsizes[pc]--;
            if( classsizes[pc] == 0 )
            {
               assert(pcsetfill < ncols);
               pcset[pcsetfill++] = pc;
            }
            pcs[i] = pc;

            i++;
         }

         assert(i > 0);

         /* sort on the pclass values */
         if( i > 1 )
         {
            SCIPsortIntIntReal(pcs, colindices, values, i);
         }

         k = 0;
         while( TRUE ) /*lint !e716*/
         {
            assert(k < i);
            startpc = pcs[k];
            startk = k;

            /* find pclass-sets */
            while( k < i && pcs[k] == startpc )
               k++;

            /* sort on the A values which have equal pclass values */
            if( k - startk > 1 )
               SCIPsortRealInt(&(values[startk]), &(colindices[startk]), k - startk);

            t = 0;
            while( TRUE ) /*lint !e716*/
            {
               assert(startk + t < i);
               startval = values[startk + t];
               startt = t;

               /* find A-sets */
               while( t < k - startk && SCIPisEQ(scip, startval, values[startk + t]) )
                  t++;

               /* get new pclass */
               newpclass = pcset[0];
               assert(pcsetfill > 0);
               pcset[0] = pcset[--pcsetfill];

               /* renumbering */
               for( m = startk + startt; m < startk + t; m++ )
               {
                  assert(m < i);
                  assert(colindices[m] < ncols);
                  assert(newpclass < ncols);

                  pclass[colindices[m]] = newpclass;
                  classsizes[newpclass]++;
               }

               if( t == k - startk )
                  break;
            }

            if( k == SCIPmatrixGetRowNNonzs(matrix, r) )
               break;
         }
      }
   }

   SCIPfreeBufferArray(scip, &pcs);
   SCIPfreeBufferArray(scip, &colindices);
   SCIPfreeBufferArray(scip, &values);
   SCIPfreeBufferArray(scip, &pcset);
   SCIPfreeBufferArray(scip, &scale);
   SCIPfreeBufferArray(scip, &classsizes);

   return SCIP_OKAY;
}


/** try to improve variable bounds by predictive bound strengthening */
static
SCIP_RETCODE predBndStr(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_VAR*             dominatingvar,      /**< dominating variable */
   int                   dominatingidx,      /**< column index of the dominating variable */
   SCIP_Real             dominatingub,       /**< predicted upper bound of the dominating variable */
   SCIP_Real             dominatinglb,       /**< predicted lower bound of the dominating variable */
   SCIP_Real             dominatingwcub,     /**< predicted worst case upper bound of the dominating variable */
   SCIP_VAR*             dominatedvar,       /**< dominated variable */
   int                   dominatedidx,       /**< column index of the dominated variable */
   SCIP_Real             dominatedub,        /**< predicted upper bound of the dominated variable */
   SCIP_Real             dominatedwclb,      /**< predicted worst case lower bound of the dominated variable */
   SCIP_Real             dominatedlb,        /**< predicted lower bound of the dominated variable */
   FIXINGDIRECTION*      varstofix,          /**< array holding fixing information */
   int*                  nchgbds             /**< count number of bound changes */
   )
{
   /* we compare only variables from the same type */
   if( !(SCIPvarGetType(dominatingvar) == SCIPvarGetType(dominatedvar) ||
         SCIPvarIsBinary(dominatingvar) == SCIPvarIsBinary(dominatedvar) ||
         (SCIPvarGetType(dominatingvar) == SCIP_VARTYPE_INTEGER && SCIPvarGetType(dominatedvar) == SCIP_VARTYPE_IMPLINT) ||
         (SCIPvarGetType(dominatedvar) == SCIP_VARTYPE_INTEGER && SCIPvarGetType(dominatingvar) == SCIP_VARTYPE_IMPLINT)) )
   {
      return SCIP_OKAY;
   }

   if( varstofix[dominatingidx] == NOFIX )
   {
      /* assume x dominates y (x->y). we get this bound from a positive coefficient
       * of the dominating variable. because x->y the dominated variable y has
       * a positive coefficient too. thus y contributes to the minactivity with its
       * lower bound. but this case is considered within predictive bound analysis.
       * thus the dominating upper bound is a common upper bound.
       */
      if( !SCIPisInfinity(scip, dominatingub) )
      {
         SCIP_Real newub;
         SCIP_Real oldub;
         SCIP_Real lb;

         newub = dominatingub;
         oldub = SCIPvarGetUbGlobal(dominatingvar);
         lb = SCIPvarGetLbGlobal(dominatingvar);

         /* if( SCIPvarGetType(dominatingvar) != SCIP_VARTYPE_CONTINUOUS )
            newub = SCIPceil(scip, newub); */

         if( SCIPisLE(scip, lb, newub) && SCIPisLT(scip, newub, oldub) )
         {
            SCIPdebugMsg(scip, "[ub]\tupper bound for dominating variable <%s> changed: [%.17f,%.17f] -> [%.17f,%.17f]\n",
               SCIPvarGetName(dominatingvar), lb, oldub, lb, newub);
            SCIP_CALL( SCIPchgVarUb(scip, dominatingvar, newub) );
            (*nchgbds)++;
         }
      }

      /* assume x dominates y (x->y). we get this lower bound of the dominating variable
       * from a negative coefficient within a <= relation. if y has a positive coefficient
       * we get a common lower bound of x from predictive bound analysis. if y has a
       * negative coefficient we get an improved lower bound of x because the minactivity
       * is greater. for discrete variables we have to round down the lower bound.
       */
      if( !SCIPisInfinity(scip, -dominatinglb) )
      {
         SCIP_Real newlb;
         SCIP_Real oldlb;
         SCIP_Real ub;

         newlb = dominatinglb;
         oldlb = SCIPvarGetLbGlobal(dominatingvar);
         ub = SCIPvarGetUbGlobal(dominatingvar);

         if( SCIPvarGetType(dominatingvar) != SCIP_VARTYPE_CONTINUOUS )
            newlb = SCIPfloor(scip, newlb);

         if( SCIPisLT(scip, oldlb, newlb) && SCIPisLE(scip, newlb, ub) )
         {
            SCIPdebugMsg(scip, "[lb]\tlower bound of dominating variable <%s> changed: [%.17f,%.17f] -> [%.17f,%.17f]\n",
               SCIPvarGetName(dominatingvar), oldlb, ub, newlb, ub);
            SCIP_CALL( SCIPchgVarLb(scip, dominatingvar, newlb) );
            (*nchgbds)++;
         }
      }

      /* assume x dominates y (x->y). we get this bound from a positive coefficient
       * of x within a <= relation. from x->y it follows, that y has a positive
       * coefficient in this row too. the worst case upper bound of x is an estimation
       * of how great x can be in every case. if the objective coefficient of x is
       * negative we get thus a lower bound of x. for discrete variables we have
       * to round down the lower bound before.
       */
      if( !SCIPisInfinity(scip, dominatingwcub) && SCIPisNegative(scip, SCIPvarGetObj(dominatingvar)))
      {
         SCIP_Real newlb;
         SCIP_Real oldlb;
         SCIP_Real ub;

         newlb = dominatingwcub;
         oldlb = SCIPvarGetLbGlobal(dominatingvar);
         ub = SCIPvarGetUbGlobal(dominatingvar);

         if( SCIPvarGetType(dominatingvar) != SCIP_VARTYPE_CONTINUOUS )
            newlb = SCIPfloor(scip, newlb);

         if( SCIPisLT(scip, oldlb, newlb) && SCIPisLE(scip, newlb, ub) )
         {
            SCIPdebugMsg(scip, "[wcub]\tlower bound of dominating variable <%s> changed: [%.17f,%.17f] -> [%.17f,%.17f]\n",
               SCIPvarGetName(dominatingvar), oldlb, ub, newlb, ub);
            SCIP_CALL( SCIPchgVarLb(scip, dominatingvar, newlb) );
            (*nchgbds)++;
         }
      }
   }

   if( varstofix[dominatedidx] == NOFIX )
   {
      /* assume x dominates y (x->y). we get this bound for a positive coefficient of y
       * within a <= relation. if x has a negative coefficient we get a common upper
       * bound of y. if x has a positive coefficient we get an improved upper bound
       * of y because the minactivity is greater.
       */
      if( !SCIPisInfinity(scip, dominatedub) )
      {
         SCIP_Real newub;
         SCIP_Real oldub;
         SCIP_Real lb;

         newub = dominatedub;
         oldub = SCIPvarGetUbGlobal(dominatedvar);
         lb = SCIPvarGetLbGlobal(dominatedvar);

         if( SCIPisLE(scip, lb, newub) && SCIPisLT(scip, newub, oldub) )
         {
            SCIPdebugMsg(scip, "[ub]\tupper bound of dominated variable <%s> changed: [%.17f,%.17f] -> [%.17f,%.17f]\n",
               SCIPvarGetName(dominatedvar), lb, oldub, lb, newub);
            SCIP_CALL( SCIPchgVarUb(scip, dominatedvar, newub) );
            (*nchgbds)++;
         }
      }

      /* assume x dominates y (x->y). we get this bound only from a negative
       * coefficient of y within a <= relation. because of x->y then x has a negative
       * coefficient too. the worst case lower bound is an estimation what property
       * the dominated variable must have if the dominating variable is at its upper bound.
       * to get an upper bound of the dominated variable we need to consider a positive
       * objective coefficient. for discrete variables we have to round up the upper bound.
       */
      if( !SCIPisInfinity(scip, -dominatedwclb) && SCIPisPositive(scip, SCIPvarGetObj(dominatedvar)) )
      {
         SCIP_Real newub;
         SCIP_Real oldub;
         SCIP_Real lb;

         newub = dominatedwclb;
         oldub = SCIPvarGetUbGlobal(dominatedvar);
         lb = SCIPvarGetLbGlobal(dominatedvar);

         if( SCIPvarGetType(dominatedvar) != SCIP_VARTYPE_CONTINUOUS )
            newub = SCIPceil(scip, newub);

         if( SCIPisLE(scip, lb, newub) && SCIPisLT(scip, newub, oldub) )
         {
            SCIPdebugMsg(scip, "[wclb]\tupper bound of dominated variable <%s> changed: [%.17f,%.17f] -> [%.17f,%.17f]\n",
               SCIPvarGetName(dominatedvar), lb, oldub, lb, newub);
            SCIP_CALL( SCIPchgVarUb(scip, dominatedvar, newub) );
            (*nchgbds)++;
         }
      }

      /* assume x dominates y (x->y). we get a lower bound of y from a negative
       * coefficient within a <= relation. but if x->y then x has a negative
       * coefficient too and contributes with its upper bound to the minactivity.
       * thus in all we have a common lower bound calculation and no rounding is necessary.
       */
      if( !SCIPisInfinity(scip, -dominatedlb) )
      {
         SCIP_Real newlb;
         SCIP_Real oldlb;
         SCIP_Real ub;

         newlb = dominatedlb;
         oldlb = SCIPvarGetLbGlobal(dominatedvar);
         ub = SCIPvarGetUbGlobal(dominatedvar);

         if( SCIPisLT(scip, oldlb, newlb) && SCIPisLE(scip, newlb, ub) )
         {
            SCIPdebugMsg(scip, "[lb]\tlower bound of dominated variable <%s> changed: [%.17f,%.17f] -> [%.17f,%.17f]\n",
               SCIPvarGetName(dominatedvar), oldlb, ub, newlb, ub);
            SCIP_CALL( SCIPchgVarLb(scip, dominatedvar, newlb) );
            (*nchgbds)++;
         }
      }
   }

   return SCIP_OKAY;
}

/** try to find variable fixings */
static
SCIP_RETCODE findFixings(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix structure */
   SCIP_VAR*             dominatingvar,      /**< dominating variable */
   int                   dominatingidx,      /**< column index of the dominating variable */
   SCIP_Real             dominatingub,       /**< predicted upper bound of the dominating variable */
   SCIP_Real             dominatingwclb,     /**< predicted worst case lower bound of the dominating variable */
   SCIP_Real             dominatinglb,       /**< predicted lower bound of the dominating variable */
   SCIP_Real             dominatingwcub,     /**< predicted worst case upper bound of the dominating variable */
   SCIP_VAR*             dominatedvar,       /**< dominated variable */
   int                   dominatedidx,       /**< column index of the dominated variable */
   FIXINGDIRECTION*      varstofix,          /**< array holding fixing information */
   SCIP_Bool             onlybinvars,        /**< flag indicating only binary variables are present */
   SCIP_Bool             onlyoneone,         /**< when onlybinvars is TRUE, flag indicates if both binary variables are in clique */
   int*                  nfixings            /**< counter for possible fixings */
   )
{
   /* we compare only variables from the same type */
   if( !(SCIPvarGetType(dominatingvar) == SCIPvarGetType(dominatedvar) ||
         SCIPvarIsBinary(dominatingvar) == SCIPvarIsBinary(dominatedvar) ||
         (SCIPvarGetType(dominatingvar) == SCIP_VARTYPE_INTEGER && SCIPvarGetType(dominatedvar) == SCIP_VARTYPE_IMPLINT) ||
         (SCIPvarGetType(dominatedvar) == SCIP_VARTYPE_INTEGER && SCIPvarGetType(dominatingvar) == SCIP_VARTYPE_IMPLINT)) )
   {
      return SCIP_OKAY;
   }

   if( varstofix[dominatedidx] == NOFIX && SCIPmatrixGetColNNonzs(matrix, dominatingidx) == 1
      && SCIPmatrixGetColNNonzs(matrix, dominatedidx) == 1 )
   {
      /* We have a x->y dominance relation and only one equality constraint
       * where the dominating variable x with an infinity upper bound and the
       * dominated variable y are present. Then the dominated variable y
       * can be fixed at its lower bound.
       */
      int row;
      row = *(SCIPmatrixGetColIdxPtr(matrix, dominatedidx));

      if( SCIPisEQ(scip, SCIPmatrixGetRowLhs(matrix, row), SCIPmatrixGetRowRhs(matrix, row)) &&
         SCIPisInfinity(scip, SCIPvarGetUbGlobal(dominatingvar)) )
      {
         varstofix[dominatedidx] = FIXATLB;
         (*nfixings)++;

         return SCIP_OKAY;
      }
   }

   if( varstofix[dominatedidx] == NOFIX && !SCIPisNegative(scip, SCIPvarGetObj(dominatedvar)) )
   {
      if( !SCIPisInfinity(scip, -dominatingwclb) &&
         SCIPisLE(scip, dominatingwclb, SCIPvarGetUbGlobal(dominatingvar)) )
      {
         /* we have a x->y dominance relation with a positive obj coefficient
          * of the dominated variable y. we need to secure feasibility
          * by testing if the predicted lower worst case bound is less equal the
          * current upper bound. it is possible, that the lower worst case bound
          * is infinity and the upper bound of the dominating variable x is
          * infinity too.
          */
         varstofix[dominatedidx] = FIXATLB;
         (*nfixings)++;
      }
   }

   if( varstofix[dominatedidx] == NOFIX && !SCIPisInfinity(scip, dominatingub) &&
      SCIPisLE(scip, dominatingub, SCIPvarGetUbGlobal(dominatingvar)) )
   {
      /* we have a x->y dominance relation with an arbitrary obj coefficient
       * of the dominating variable x. in all cases we have to look
       * if the predicted upper bound of the dominating variable is great enough.
       * by testing, that the predicted upper bound is not infinity we avoid problems
       * with x->y e.g.
       *    min  -x -y
       *    s.t. -x -y <= -1
       *    0<=x<=1, 0<=y<=1
       * where y is not at their lower bound.
       */
      varstofix[dominatedidx] = FIXATLB;
      (*nfixings)++;
   }

   if( varstofix[dominatingidx] == NOFIX && !SCIPisPositive(scip, SCIPvarGetObj(dominatingvar)) )
   {
      /* we have a x->y dominance relation with a negative obj coefficient
       * of the dominating variable x. if the worst case upper bound is
       * greater equal than upper bound, we fix x at the upper bound
       */
      if( !SCIPisInfinity(scip, dominatingwcub) &&
         SCIPisGE(scip, dominatingwcub, SCIPvarGetUbGlobal(dominatingvar)) )
      {
         varstofix[dominatingidx] = FIXATUB;
         (*nfixings)++;
      }
   }

   if( varstofix[dominatingidx] == NOFIX && !SCIPisInfinity(scip, -dominatinglb) &&
      SCIPisGE(scip, dominatinglb, SCIPvarGetUbGlobal(dominatingvar)) )
   {
       /* we have a x->y dominance relation with an arbitrary obj coefficient
        * of the dominating variable x. if the predicted lower bound is greater
        * equal than upper bound, we fix x at the upper bound.
        */
      varstofix[dominatingidx] = FIXATUB;
      (*nfixings)++;
   }

   if( onlybinvars )
   {
      if( varstofix[dominatedidx] == NOFIX && (onlyoneone || SCIPvarsHaveCommonClique(dominatingvar, TRUE, dominatedvar, TRUE, TRUE)) )
      {
         /* We have a (1->1)-clique with dominance relation (x->y) (x dominates y).
          * From this dominance relation, we know (1->0) is possible and not worse than (0->1)
          * concerning the objective function. It follows that only (1->0) or (0->0) are possible,
          * but in both cases y has the value 0 => y=0.
          */
         varstofix[dominatedidx] = FIXATLB;
         (*nfixings)++;
      }

      if( varstofix[dominatingidx] == NOFIX && SCIPvarsHaveCommonClique(dominatingvar, FALSE, dominatedvar, FALSE, TRUE) )
      {
         /* We have a (0->0)-clique with dominance relation x->y (x dominates y).
          * From this dominance relation, we know (1->0) is possible and not worse than (0->1)
          * concerning the objective function. It follows that only (1->0) or (1->1) are possible,
          * but in both cases x has the value 1 => x=1
          */
         varstofix[dominatingidx] = FIXATUB;
         (*nfixings)++;
      }
   }
   else
      assert(!onlyoneone);

   return SCIP_OKAY;
}

/** find dominance relation between variable pairs */
static
SCIP_RETCODE findDominancePairs(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   int*                  searchcols,         /**< indexes of variables for pair comparisons */
   int                   searchsize,         /**< number of variables for pair comparisons */
   SCIP_Bool             onlybinvars,        /**< flag indicating searchcols contains only binary variable indexes */
   FIXINGDIRECTION*      varstofix,          /**< array holding information for later upper/lower bound fixing */
   int*                  nfixings,           /**< found number of possible fixings */
   SCIP_Longint*         ndomrelations,      /**< found number of dominance relations */
   int*                  nchgbds             /**< number of changed bounds */
   )
{
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   SCIP_Real tmpupperbounddominatingcol1;
   SCIP_Real tmpupperbounddominatingcol2;
   SCIP_Real tmpwclowerbounddominatingcol1;
   SCIP_Real tmpwclowerbounddominatingcol2;
   SCIP_Real tmplowerbounddominatingcol1;
   SCIP_Real tmplowerbounddominatingcol2;
   SCIP_Real tmpwcupperbounddominatingcol1;
   SCIP_Real tmpwcupperbounddominatingcol2;
   int* rows1;
   int* rows2;
   int nrows1;
   int nrows2;
   SCIP_Real tmpupperbounddominatedcol1;
   SCIP_Real tmpupperbounddominatedcol2;
   SCIP_Real tmpwclowerbounddominatedcol1;
   SCIP_Real tmpwclowerbounddominatedcol2;
   SCIP_Real tmplowerbounddominatedcol1;
   SCIP_Real tmplowerbounddominatedcol2;
   SCIP_Real tmpwcupperbounddominatedcol1;
   SCIP_Real tmpwcupperbounddominatedcol2;
   SCIP_Real obj1;
   SCIP_Bool col1domcol2;
   SCIP_Bool col2domcol1;
   SCIP_Bool onlyoneone;
   int cnt1;
   int cnt2;
   int col1;
   int col2;
   int r1;
   int r2;
   int paircnt;
   int oldnfixings;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(presoldata != NULL);
   assert(searchcols != NULL);
   assert(varstofix != NULL);
   assert(nfixings != NULL);
   assert(ndomrelations != NULL);
   assert(nchgbds != NULL);

   paircnt = 0;
   oldnfixings = *nfixings;

   /* pair comparisons */
   for( cnt1 = 0; cnt1 < searchsize; cnt1++ )
   {
      SCIP_VAR* varcol1;
      SCIP_VAR* varcol2;

      /* get index of the first variable */
      col1 = searchcols[cnt1];

      if( varstofix[col1] == FIXATLB )
         continue;

      varcol1 = SCIPmatrixGetVar(matrix, col1);
      obj1 = SCIPvarGetObj(varcol1);

      for( cnt2 = cnt1+1; cnt2 < searchsize; cnt2++ )
      {
         /* get index of the second variable */
         col2 = searchcols[cnt2];
         varcol2 = SCIPmatrixGetVar(matrix, col2);
         onlyoneone = FALSE;

         /* we always have minimize as obj sense */

         /* column 1 dominating column 2 */
         col1domcol2 = (obj1 <= SCIPvarGetObj(varcol2));

         /* column 2 dominating column 1 */
         col2domcol1 = (SCIPvarGetObj(varcol2) <= obj1);

         /* search only if nothing was found yet */
         col1domcol2 = col1domcol2 && (varstofix[col2] == NOFIX);
         col2domcol1 = col2domcol1 && (varstofix[col1] == NOFIX);

         /* we only search for a dominance relation if the lower bounds are not negative */
         if( !onlybinvars )
         {
            if( SCIPisLT(scip, SCIPvarGetLbGlobal(varcol1), 0.0) ||
               SCIPisLT(scip, SCIPvarGetLbGlobal(varcol2), 0.0) )
            {
               col1domcol2 = FALSE;
               col2domcol1 = FALSE;
            }
         }

         /* pair comparison control */
         if( paircnt == presoldata->numcurrentpairs )
         {
            assert(*nfixings >= oldnfixings);
            if( *nfixings == oldnfixings )
            {
               /* not enough fixings found, decrement number of comparisons */
               presoldata->numcurrentpairs >>= 1; /*lint !e702*/
               if( presoldata->numcurrentpairs < presoldata->numminpairs )
                  presoldata->numcurrentpairs = presoldata->numminpairs;

               /* stop searching in this row */
               return SCIP_OKAY;
            }
            oldnfixings = *nfixings;
            paircnt = 0;

            /* increment number of comparisons */
            presoldata->numcurrentpairs <<= 1; /*lint !e701*/
            if( presoldata->numcurrentpairs > presoldata->nummaxpairs )
               presoldata->numcurrentpairs = presoldata->nummaxpairs;
         }
         paircnt++;

         if( !col1domcol2 && !col2domcol1 )
            continue;

         /* get the data for both columns */
         vals1 = SCIPmatrixGetColValPtr(matrix, col1);
         rows1 = SCIPmatrixGetColIdxPtr(matrix, col1);
         nrows1 = SCIPmatrixGetColNNonzs(matrix, col1);
         r1 = 0;
         vals2 = SCIPmatrixGetColValPtr(matrix, col2);
         rows2 = SCIPmatrixGetColIdxPtr(matrix, col2);
         nrows2 = SCIPmatrixGetColNNonzs(matrix, col2);
         r2 = 0;

         /* do we have a obj constant? */
         if( nrows1 == 0 || nrows2 == 0 )
            continue;

         /* initialize temporary bounds of dominating variable */
         tmpupperbounddominatingcol1 = SCIPinfinity(scip);
         tmpupperbounddominatingcol2 = tmpupperbounddominatingcol1;
         tmpwclowerbounddominatingcol1 = -SCIPinfinity(scip);
         tmpwclowerbounddominatingcol2 = tmpwclowerbounddominatingcol1;
         tmplowerbounddominatingcol1 = -SCIPinfinity(scip);
         tmplowerbounddominatingcol2 = tmplowerbounddominatingcol1;
         tmpwcupperbounddominatingcol1 = SCIPinfinity(scip);
         tmpwcupperbounddominatingcol2 = tmpwcupperbounddominatingcol1;

         /* initialize temporary bounds of dominated variable */
         tmpupperbounddominatedcol1 = SCIPinfinity(scip);
         tmpupperbounddominatedcol2 = tmpupperbounddominatedcol1;
         tmpwclowerbounddominatedcol1 = -SCIPinfinity(scip);
         tmpwclowerbounddominatedcol2 = tmpwclowerbounddominatedcol1;
         tmplowerbounddominatedcol1 = -SCIPinfinity(scip);
         tmplowerbounddominatedcol2 = tmplowerbounddominatedcol1;
         tmpwcupperbounddominatedcol1 = SCIPinfinity(scip);
         tmpwcupperbounddominatedcol2 = tmpwcupperbounddominatedcol1;

         /* compare rows of this column pair */
         while( (col1domcol2 || col2domcol1) && (r1 < nrows1 || r2 < nrows2) )
         {
            assert((r1 >= nrows1-1) || (rows1[r1] < rows1[r1+1]));
            assert((r2 >= nrows2-1) || (rows2[r2] < rows2[r2+1]));

            /* there is a nonredundant row containing column 1 but not column 2 */
            if( r1 < nrows1 && (r2 == nrows2 || rows1[r1] < rows2[r2]) )
            {
               /* dominance depends on the type of relation */
               if( !SCIPmatrixIsRowRhsInfinity(matrix, rows1[r1]) )
               {
                  /* no dominance relation for equations or ranged rows */
                  col2domcol1 = FALSE;
                  col1domcol2 = FALSE;
               }
               else
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals1[r1] > 0.0 )
                     col2domcol1 = FALSE;
                  else
                     col1domcol2 = FALSE;
               }

               r1++;
            }
            /* there is a nonredundant row containing column 2, but not column 1 */
            else if( r2 < nrows2 && (r1 == nrows1 || rows1[r1] > rows2[r2]) )
            {
               /* dominance depends on the type of relation */
               if( !SCIPmatrixIsRowRhsInfinity(matrix, rows2[r2]) )
               {
                  /* no dominance relation for equations or ranged rows */
                  col2domcol1 = FALSE;
                  col1domcol2 = FALSE;
               }
               else
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals2[r2] < 0.0 )
                     col2domcol1 = FALSE;
                  else
                     col1domcol2 = FALSE;
               }

               r2++;
            }
            /* if both columns appear in a common row, compare the coefficients */
            else
            {
               assert(r1 < nrows1 && r2 < nrows2);
               assert(rows1[r1] == rows2[r2]);

               /* if both columns are binary variables we check if they have a common clique
                * and do not calculate any bounds
                */
               if( onlybinvars && !onlyoneone )
               {
                  if( vals1[r1] < 0 && vals2[r2] < 0 )
                  {
                     if( (SCIPmatrixGetRowNMaxActPosInf(matrix, rows1[r1]) + SCIPmatrixGetRowNMaxActNegInf(matrix, rows1[r1]) == 0)
                        && SCIPisFeasLE(scip, SCIPmatrixGetRowMaxActivity(matrix, rows1[r1]) + MAX(vals1[r1], vals2[r2]),
                           SCIPmatrixGetRowLhs(matrix, rows1[r1])) )
                     {
                        onlyoneone = TRUE;
                     }
                  }

                  if( !onlyoneone && !SCIPmatrixIsRowRhsInfinity(matrix, rows1[r1]) )
                  {
                     if ( vals1[r1] > 0 && vals2[r2] > 0 )
                     {
                        if( (SCIPmatrixGetRowNMinActPosInf(matrix, rows1[r1]) + SCIPmatrixGetRowNMinActNegInf(matrix, rows1[r1]) == 0)
                           && SCIPisFeasGE(scip, SCIPmatrixGetRowMinActivity(matrix, rows1[r1]) + MIN(vals1[r1], vals2[r2]),
                              SCIPmatrixGetRowRhs(matrix, rows1[r1])) )
                        {
                           onlyoneone = TRUE;
                        }
                     }
                  }

                  if( onlyoneone )
                  {
                     /* reset bounds */
                     tmpupperbounddominatingcol1 = SCIPinfinity(scip);
                     tmpupperbounddominatingcol2 = tmpupperbounddominatingcol1;
                     tmpwclowerbounddominatingcol1 = -SCIPinfinity(scip);
                     tmpwclowerbounddominatingcol2 = tmpwclowerbounddominatingcol1;
                     tmplowerbounddominatingcol1 = -SCIPinfinity(scip);
                     tmplowerbounddominatingcol2 = tmplowerbounddominatingcol1;
                     tmpwcupperbounddominatingcol1 = SCIPinfinity(scip);
                     tmpwcupperbounddominatingcol2 = tmpwcupperbounddominatingcol1;

                     tmpupperbounddominatedcol1 = SCIPinfinity(scip);
                     tmpupperbounddominatedcol2 = tmpupperbounddominatedcol1;
                     tmpwclowerbounddominatedcol1 = -SCIPinfinity(scip);
                     tmpwclowerbounddominatedcol2 = tmpwclowerbounddominatedcol1;
                     tmplowerbounddominatedcol1 = -SCIPinfinity(scip);
                     tmplowerbounddominatedcol2 = tmplowerbounddominatedcol1;
                     tmpwcupperbounddominatedcol1 = SCIPinfinity(scip);
                     tmpwcupperbounddominatedcol2 = tmpwcupperbounddominatedcol1;
                  }
               }

               /* dominance depends on the type of inequality */
               if( !SCIPmatrixIsRowRhsInfinity(matrix, rows1[r1]) )
               {
                  /* no dominance relation if coefficients differ for equations or ranged rows */
                  if( !SCIPisEQ(scip, vals1[r1], vals2[r2]) )
                  {
                     col2domcol1 = FALSE;
                     col1domcol2 = FALSE;
                  }
               }
               else
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals1[r1] > vals2[r2] )
                     col2domcol1 = FALSE;
                  else if( vals1[r1] < vals2[r2] )
                     col1domcol2 = FALSE;
               }

               /* we do not use bound calulations if two binary variable are in one common clique.
                * for the other cases we claim the same sign for the coefficients to
                * achieve monotonically decreasing predictive bound functions.
                */
               if( !onlyoneone &&
                  ((vals1[r1] < 0 && vals2[r2] < 0) || (vals1[r1] > 0 && vals2[r2] > 0)) )
               {
                  if( col1domcol2 )
                  {
                     /* update bounds of dominating variable for column 1 */
                     SCIP_CALL( updateBounds(scip, matrix, rows1[r1],
                           col1, vals1[r1], col2, vals2[r2], TRUE,
                           &tmpupperbounddominatingcol1, &tmpwclowerbounddominatingcol1,
                           &tmplowerbounddominatingcol1, &tmpwcupperbounddominatingcol1) );

                     /* update bounds of dominated variable for column 1 */
                     SCIP_CALL( updateBounds(scip, matrix, rows1[r1],
                           col1, vals1[r1], col2, vals2[r2], FALSE,
                           &tmpupperbounddominatedcol1, &tmpwclowerbounddominatedcol1,
                           &tmplowerbounddominatedcol1, &tmpwcupperbounddominatedcol1) );
                  }

                  if( col2domcol1 )
                  {
                     /* update bounds of dominating variable for column 2 */
                     SCIP_CALL( updateBounds(scip, matrix, rows2[r2],
                           col2, vals2[r2], col1, vals1[r1], TRUE,
                           &tmpupperbounddominatingcol2, &tmpwclowerbounddominatingcol2,
                           &tmplowerbounddominatingcol2, &tmpwcupperbounddominatingcol2) );

                     /* update bounds of dominated variable for column 2 */
                     SCIP_CALL( updateBounds(scip, matrix, rows2[r2],
                           col2, vals2[r2], col1, vals1[r1], FALSE,
                           &tmpupperbounddominatedcol2, &tmpwclowerbounddominatedcol2,
                           &tmplowerbounddominatedcol2, &tmpwcupperbounddominatedcol2) );
                  }
               }

               r1++;
               r2++;
            }
         }

         /* a column is only dominated, if there are no more rows in which it is contained */
         col1domcol2 = col1domcol2 && r2 == nrows2;
         col2domcol1 = col2domcol1 && r1 == nrows1;

         if( !col1domcol2 && !col2domcol1 )
            continue;

         /* no dominance relation for left equations or ranged rows */
         while( r1 < nrows1 )
         {
            if( !SCIPmatrixIsRowRhsInfinity(matrix, rows1[r1]) )
            {
               col2domcol1 = FALSE;
               col1domcol2 = FALSE;
               break;
            }
            r1++;
         }
         if( !col1domcol2 && !col2domcol1 )
            continue;
         while( r2 < nrows2 )
         {
            if( !SCIPmatrixIsRowRhsInfinity(matrix, rows2[r2]) )
            {
               col2domcol1 = FALSE;
               col1domcol2 = FALSE;
               break;
            }
            r2++;
         }

         if( col1domcol2 || col2domcol1 )
            (*ndomrelations)++;

         if( col1domcol2 && col2domcol1 )
         {
            /* prefer the variable as dominating variable with the greater upper bound */
            if( SCIPisGE(scip, SCIPvarGetUbGlobal(varcol1), SCIPvarGetUbGlobal(varcol2)) )
               col2domcol1 = FALSE;
            else
               col1domcol2 = FALSE;
         }

         /* use dominance relation and clique/bound-information
          * to find variable fixings. additionally try to strengthen
          * variable bounds by predictive bound strengthening.
          */
         if( col1domcol2 )
         {
            SCIP_CALL( findFixings(scip, matrix, varcol1, col1,
                  tmpupperbounddominatingcol1, tmpwclowerbounddominatingcol1,
                  tmplowerbounddominatingcol1, tmpwcupperbounddominatingcol1,
                  varcol2, col2,
                  varstofix, onlybinvars, onlyoneone, nfixings) );

            if( presoldata->predbndstr )
            {
               SCIP_CALL( predBndStr(scip, varcol1, col1,
                     tmpupperbounddominatingcol1,
                     tmplowerbounddominatingcol1, tmpwcupperbounddominatingcol1,
                     varcol2, col2,
                     tmpupperbounddominatedcol1, tmpwclowerbounddominatedcol1,
                     tmplowerbounddominatedcol1,
                     varstofix, nchgbds) );
            }
         }
         else if( col2domcol1 )
         {
            SCIP_CALL( findFixings(scip, matrix, varcol2, col2,
                  tmpupperbounddominatingcol2, tmpwclowerbounddominatingcol2,
                  tmplowerbounddominatingcol2, tmpwcupperbounddominatingcol2,
                  varcol1, col1,
                  varstofix, onlybinvars, onlyoneone, nfixings) );

            if( presoldata->predbndstr )
            {
               SCIP_CALL( predBndStr(scip, varcol2, col2,
                     tmpupperbounddominatingcol2,
                     tmplowerbounddominatingcol2, tmpwcupperbounddominatingcol2,
                     varcol1, col1,
                     tmpupperbounddominatedcol2, tmpwclowerbounddominatedcol2,
                     tmplowerbounddominatedcol2,
                     varstofix, nchgbds) );
            }
         }
         if( varstofix[col1] == FIXATLB )
            break;
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDomcol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDomcol(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeDomcol)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDomcol)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* don't run for pure LPs */
   if( !presoldata->continuousred && (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == 0) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized && complete )
   {
      int nfixings;
      SCIP_Longint ndomrelations;
      int v;
      int r;
      FIXINGDIRECTION* varstofix;
      SCIP_Bool* varsprocessed;
      int nrows;
      int ncols;
      int* rowidxsorted;
      int* rowsparsity;
      int varcount;
      int* consearchcols;
      int* intsearchcols;
      int* binsearchcols;
      int nconfill;
      int nintfill;
      int nbinfill;
#ifdef SCIP_DEBUG
      int nconvarsfixed = 0;
      int nintvarsfixed = 0;
      int nbinvarsfixed = 0;
#endif
      int* pclass;
      int* colidx;
      int pclassstart;
      int pc;
      SCIP_Bool* varineq;

      nfixings = 0;
      ndomrelations = 0;
      nrows = SCIPmatrixGetNRows(matrix);
      ncols = SCIPmatrixGetNColumns(matrix);
      assert(SCIPgetNVars(scip) == ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, ncols) );
      BMSclearMemoryArray(varstofix, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &varsprocessed, ncols) );
      BMSclearMemoryArray(varsprocessed, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &consearchcols, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &intsearchcols, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &binsearchcols, ncols) );

      SCIP_CALL( SCIPallocBufferArray(scip, &rowidxsorted, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowsparsity, nrows) );
      for( r = 0; r < nrows; ++r )
      {
         rowidxsorted[r] = r;
         rowsparsity[r] = SCIPmatrixGetRowNNonzs(matrix, r);
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &pclass, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colidx, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varineq, ncols) );
      for( v = 0; v < ncols; v++ )
      {
         colidx[v] = v;
         varineq[v] = FALSE;
      }

      /* init pair comparision control */
      presoldata->numcurrentpairs = presoldata->nummaxpairs;

      varcount = 0;

      /* 1.stage: search dominance relations of parallel columns
       *          within equalities and ranged rows
       */
      if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
      {
         SCIP_CALL( detectParallelCols(scip, matrix, pclass, varineq) );
         SCIPsortIntInt(pclass, colidx, ncols);

         pc = 0;
         while( pc < ncols )
         {
            int varidx;

            varidx = 0;
            nconfill = 0;
            nintfill = 0;
            nbinfill = 0;

            pclassstart = pclass[pc];
            while( pc < ncols && pclassstart == pclass[pc] )
            {
               SCIP_VAR* var;

               varidx = colidx[pc];
               var = SCIPmatrixGetVar(matrix, varidx);

               /* we only regard variables which were not processed yet and
                * are present within equalities or ranged rows
                */
               if( !varsprocessed[varidx] && varineq[varidx] )
               {
                  /* we search only for dominance relations between the same variable type */
                  if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  {
                     consearchcols[nconfill++] = varidx;
                  }
                  else if( SCIPvarIsBinary(var) )
                  {
                     binsearchcols[nbinfill++] = varidx;
                  }
                  else
                  {
                     assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
                     intsearchcols[nintfill++] = varidx;
                  }
               }
               ++pc;
            }

            /* continuous variables */
            if( nconfill > 1 && presoldata->continuousred )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, presoldata, consearchcols, nconfill, FALSE,
                     varstofix, &nfixings, &ndomrelations, nchgbds) );

               for( v = 0; v < nconfill; ++v )
                  varsprocessed[consearchcols[v]] = TRUE;

               varcount += nconfill;
            }
            else if( nconfill == 1 )
            {
               if( varineq[varidx] )
                  varsprocessed[consearchcols[0]] = TRUE;
            }

            /* integer and impl-integer variables */
            if( nintfill > 1 )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, presoldata, intsearchcols, nintfill, FALSE,
                     varstofix, &nfixings, &ndomrelations, nchgbds) );

               for( v = 0; v < nintfill; ++v )
                  varsprocessed[intsearchcols[v]] = TRUE;

               varcount += nintfill;
            }
            else if( nintfill == 1 )
            {
               if( varineq[varidx] )
                  varsprocessed[intsearchcols[0]] = TRUE;
            }

            /* binary variables */
            if( nbinfill > 1 )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, presoldata, binsearchcols, nbinfill, TRUE,
                     varstofix, &nfixings, &ndomrelations, nchgbds) );

               for( v = 0; v < nbinfill; ++v )
                  varsprocessed[binsearchcols[v]] = TRUE;

               varcount += nbinfill;
            }
            else if( nbinfill == 1 )
            {
               if( varineq[varidx] )
                  varsprocessed[binsearchcols[0]] = TRUE;
            }

            if( varcount >= ncols )
               break;
         }
      }

      /* 2.stage: search dominance relations for the remaining columns
       *          by increasing row-sparsity
       */
      if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
      {
         SCIPsortIntInt(rowsparsity, rowidxsorted, nrows);

         for( r = 0; r < nrows; ++r )
         {
            int rowidx;
            int* rowpnt;
            int* rowend;

            /* break if the time limit was reached; since the check is expensive,
             * we only check all 1000 constraints
             */
            if( (r % 1000 == 0) && SCIPisStopped(scip) )
               break;

            rowidx = rowidxsorted[r];
            rowpnt = SCIPmatrixGetRowIdxPtr(matrix, rowidx);
            rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, rowidx);

            if( SCIPmatrixGetRowNNonzs(matrix, rowidx) == 1 )
               continue;

            nconfill = 0;
            nintfill = 0;
            nbinfill = 0;

            for( ; rowpnt < rowend; rowpnt++ )
            {
               if( !(varsprocessed[*rowpnt]) )
               {
                  int varidx;
                  SCIP_VAR* var;

                  varidx = *rowpnt;
                  var = SCIPmatrixGetVar(matrix, varidx);

                  /* we search only for dominance relations between the same variable type */
                  if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  {
                     consearchcols[nconfill++] = varidx;
                  }
                  else if( SCIPvarIsBinary(var) )
                  {
                     binsearchcols[nbinfill++] = varidx;
                  }
                  else
                  {
                     assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
                     intsearchcols[nintfill++] = varidx;
                  }
               }
            }

            /* continuous variables */
            if( nconfill > 1 && presoldata->continuousred )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, presoldata, consearchcols, nconfill, FALSE,
                     varstofix, &nfixings, &ndomrelations, nchgbds) );

               for( v = 0; v < nconfill; ++v )
                  varsprocessed[consearchcols[v]] = TRUE;

               varcount += nconfill;
            }

            /* integer and impl-integer variables */
            if( nintfill > 1 )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, presoldata, intsearchcols, nintfill, FALSE,
                     varstofix, &nfixings, &ndomrelations, nchgbds) );

               for( v = 0; v < nintfill; ++v )
                  varsprocessed[intsearchcols[v]] = TRUE;

               varcount += nintfill;
            }

            /* binary variables */
            if( nbinfill > 1 )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, presoldata, binsearchcols, nbinfill, TRUE,
                     varstofix, &nfixings, &ndomrelations, nchgbds) );

               for( v = 0; v < nbinfill; ++v )
                  varsprocessed[binsearchcols[v]] = TRUE;

               varcount += nbinfill;
            }

            if( varcount >= ncols )
               break;
         }
      }

      if( nfixings > 0 )
      {
         int oldnfixedvars;

         oldnfixedvars = *nfixedvars;

         for( v = ncols - 1; v >= 0; --v )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            SCIP_VAR* var;

            var = SCIPmatrixGetVar(matrix,v);

            if( SCIPvarGetNLocksUp(var) != SCIPmatrixGetColNUplocks(matrix, v) ||
               SCIPvarGetNLocksDown(var) != SCIPmatrixGetColNDownlocks(matrix, v) )
            {
               /* no fixing, locks not consistent */
               continue;
            }

            if( varstofix[v] == FIXATLB )
            {
               SCIP_Real lb;

               lb = SCIPvarGetLbGlobal(var);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPisFeasIntegral(scip, lb));

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

#ifdef SCIP_DEBUG
               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  nconvarsfixed++;
               else if( SCIPvarIsBinary(var) )
                  nbinvarsfixed++;
               else
                  nintvarsfixed++;
#endif
            }
            else if( varstofix[v] == FIXATUB )
            {
               SCIP_Real ub;

               ub = SCIPvarGetUbGlobal(var);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPisFeasIntegral(scip, ub));

               /* fix at upper bound */
               SCIP_CALL( SCIPfixVar(scip, var, ub, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMsg(scip, " -> infeasible fixing\n");
                  *result = SCIP_CUTOFF;

                  break;
               }
               assert(fixed);
               (*nfixedvars)++;

#ifdef SCIP_DEBUG
               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  nconvarsfixed++;
               else if( SCIPvarIsBinary(var) )
                  nbinvarsfixed++;
               else
                  nintvarsfixed++;
#endif
            }
         }

         if( *result != SCIP_CUTOFF && *nfixedvars > oldnfixedvars )
            *result = SCIP_SUCCESS;
      }

      SCIPfreeBufferArray(scip, &varineq);
      SCIPfreeBufferArray(scip, &colidx);
      SCIPfreeBufferArray(scip, &pclass);
      SCIPfreeBufferArray(scip, &rowsparsity);
      SCIPfreeBufferArray(scip, &rowidxsorted);
      SCIPfreeBufferArray(scip, &binsearchcols);
      SCIPfreeBufferArray(scip, &intsearchcols);
      SCIPfreeBufferArray(scip, &consearchcols);
      SCIPfreeBufferArray(scip, &varsprocessed);
      SCIPfreeBufferArray(scip, &varstofix);

#ifdef SCIP_DEBUG
      if( (nconvarsfixed + nintvarsfixed + nbinvarsfixed) > 0 )
      {
         SCIPdebugMsg(scip, "### %d vars [%" SCIP_LONGINT_FORMAT " dom] => fixed [cont: %d, int: %d, bin: %d], %scutoff detected\n",
            ncols, ndomrelations, nconvarsfixed, nintvarsfixed, nbinvarsfixed, (*result != SCIP_CUTOFF) ? "no " : "");
      }
#endif
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the domcol presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDomcol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create domcol presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecDomcol, presoldata) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDomcol) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeDomcol) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/domcol/numminpairs",
         "minimal number of pair comparisons",
         &presoldata->numminpairs, FALSE, DEFAULT_NUMMINPAIRS, 100, DEFAULT_NUMMAXPAIRS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/domcol/nummaxpairs",
         "maximal number of pair comparisons",
         &presoldata->nummaxpairs, FALSE, DEFAULT_NUMMAXPAIRS, DEFAULT_NUMMINPAIRS, 1000000000, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/domcol/predbndstr",
         "should predictive bound strengthening be applied?",
         &presoldata->predbndstr, FALSE, DEFAULT_PREDBNDSTR, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/domcol/continuousred",
         "should reductions for continuous variables be performed?",
         &presoldata->continuousred, FALSE, DEFAULT_CONTINUOUS_RED, NULL, NULL) );

   return SCIP_OKAY;
}
