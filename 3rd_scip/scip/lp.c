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

/**@file   lp.c
 * @brief  LP management methods and data structures
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gerald Gamrath
 *
 *  In LP management, we have to differ between the current LP and the SCIP_LP
 *  stored in the LP solver. All LP methods affect the current LP only.
 *  Before solving the current LP with the LP solver or setting an LP state,
 *  the LP solvers data has to be updated to the current LP with a call to
 *  lpFlush().
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/intervalarith.h"
#include "scip/clock.h"
#include "scip/misc.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/event.h"
#include "scip/pub_message.h"
#include "lpi/lpi.h"



/*
 * debug messages
 */

#ifdef SCIP_DEBUG
/** method is to print in row in case SCIP_DEBUG is defined */
static
void debugRowPrint(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   int i;

   assert(row != NULL);

   /* print row name */
   if( row->name != NULL && row->name[0] != '\0' )
   {
      SCIPsetDebugMsgPrint(set, "%s: ", row->name);
   }

   /* print left hand side */
   SCIPsetDebugMsgPrint(set, "%.15g <= ", row->lhs);

   /* print coefficients */
   if( row->len == 0 )
   {
      SCIPsetDebugMsgPrint(set, "0 ");
   }
   for( i = 0; i < row->len; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(row->cols[i]->var != NULL);
      assert(SCIPvarGetName(row->cols[i]->var) != NULL);
      assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
      SCIPsetDebugMsgPrint(set, "%+.15g<%s> ", row->vals[i], SCIPvarGetName(row->cols[i]->var));
   }

   /* print constant */
   if( REALABS(row->constant) > SCIP_DEFAULT_EPSILON )
   {
      SCIPsetDebugMsgPrint(set, "%+.15g ", row->constant);
   }

   /* print right hand side */
   SCIPsetDebugMsgPrint(set, "<= %.15g\n", row->rhs);
}
#else
#define debugRowPrint(x,y) /**/
#endif

#ifdef SCIP_DEBUG
/** method to output column if SCIP_DEBUG is define */
static
void debugColPrint(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col                 /**< LP column */
   )
{
   int r;

   assert(col != NULL);
   assert(col->var != NULL);

   /* print bounds */
   SCIPsetDebugMsgPrint(set, "(obj: %.15g) [%.15g,%.15g], ", col->obj, col->lb, col->ub);

   /* print coefficients */
   if( col->len == 0 )
   {
      SCIPsetDebugMsgPrint(set, "<empty>");
   }
   for( r = 0; r < col->len; ++r )
   {
      assert(col->rows[r] != NULL);
      assert(col->rows[r]->name != NULL);
      SCIPsetDebugMsgPrint(set, "%+.15g<%s> ", col->vals[r], col->rows[r]->name);
   }
   SCIPsetDebugMsgPrint(set, "\n");
}
#else
#define debugColPrint(x,y) /**/
#endif

/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that chgcols array can store at least num entries */
static
SCIP_RETCODE ensureChgcolsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgcols <= lp->chgcolssize);

   if( num > lp->chgcolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->chgcols, newsize) );
      lp->chgcolssize = newsize;
   }
   assert(num <= lp->chgcolssize);

   return SCIP_OKAY;
}

/** ensures, that chgrows array can store at least num entries */
static
SCIP_RETCODE ensureChgrowsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nchgrows <= lp->chgrowssize);

   if( num > lp->chgrowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->chgrows, newsize) );
      lp->chgrowssize = newsize;
   }
   assert(num <= lp->chgrowssize);

   return SCIP_OKAY;
}

/** ensures, that lpicols array can store at least num entries */
static
SCIP_RETCODE ensureLpicolsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlpicols <= lp->lpicolssize);

   if( num > lp->lpicolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->lpicols, newsize) );
      lp->lpicolssize = newsize;
   }
   assert(num <= lp->lpicolssize);

   return SCIP_OKAY;
}

/** ensures, that lpirows array can store at least num entries */
static
SCIP_RETCODE ensureLpirowsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlpirows <= lp->lpirowssize);

   if( num > lp->lpirowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->lpirows, newsize) );
      lp->lpirowssize = newsize;
   }
   assert(num <= lp->lpirowssize);

   return SCIP_OKAY;
}

/** ensures, that cols array can store at least num entries */
static
SCIP_RETCODE ensureColsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->ncols <= lp->colssize);

   if( num > lp->colssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->cols, newsize) );
      lp->colssize = newsize;
   }
   assert(num <= lp->colssize);

   return SCIP_OKAY;
}

/** ensures, that lazy cols array can store at least num entries */
static
SCIP_RETCODE ensureLazycolsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nlazycols <= lp->lazycolssize);

   if( num > lp->lazycolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->lazycols, newsize) );
      lp->lazycolssize = newsize;
   }
   assert(num <= lp->lazycolssize);

   return SCIP_OKAY;
}

/** ensures, that rows array can store at least num entries */
static
SCIP_RETCODE ensureRowsSize(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(lp->nrows <= lp->rowssize);

   if( num > lp->rowssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&lp->rows, newsize) );
      lp->rowssize = newsize;
   }
   assert(num <= lp->rowssize);

   return SCIP_OKAY;
}

/** ensures, that row array of column can store at least num entries */
static
SCIP_RETCODE colEnsureSize(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(col != NULL);
   assert(col->len <= col->size);

   if( num > col->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->rows, col->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->vals, col->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &col->linkpos, col->size, newsize) );
      col->size = newsize;
   }
   assert(num <= col->size);

   return SCIP_OKAY;
}

/** save current LP values dependent on the solution */
static
SCIP_RETCODE lpStoreSolVals(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_STAT*            stat,               /**< problem statistics */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_LPSOLVALS* storedsolvals;

   assert(lp != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);

   /* allocate memory for storage */
   if( lp->storedsolvals == NULL )
   {
      SCIP_ALLOC( BMSallocMemory(&lp->storedsolvals) );
   }
   storedsolvals = lp->storedsolvals;

   /* store values */
   storedsolvals->lpsolstat = lp->lpsolstat;
   storedsolvals->lpobjval = lp->lpobjval;
   storedsolvals->primalfeasible = lp->primalfeasible;
   storedsolvals->primalchecked = lp->primalchecked;
   storedsolvals->dualfeasible = lp->dualfeasible;
   storedsolvals->dualchecked = lp->dualchecked;
   storedsolvals->solisbasic = lp->solisbasic;
   storedsolvals->lpissolved = lp->solved;

   return SCIP_OKAY;
}

/** restore LP solution values in column */
static
SCIP_RETCODE lpRestoreSolVals(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Longint          validlp             /**< number of lp for which restored values are valid */
   )
{
   SCIP_LPSOLVALS* storedsolvals;

   assert(lp != NULL);
   assert(blkmem != NULL);

   /* if stored values are available, restore them */
   storedsolvals = lp->storedsolvals;
   if( storedsolvals != NULL )
   {
      lp->solved = storedsolvals->lpissolved;
      lp->validsollp = validlp;

      lp->lpsolstat = storedsolvals->lpsolstat;
      lp->lpobjval = storedsolvals->lpobjval;
      lp->primalfeasible = storedsolvals->primalfeasible;
      lp->primalchecked = storedsolvals->primalchecked;
      lp->dualfeasible = storedsolvals->dualfeasible;
      lp->dualchecked = storedsolvals->dualchecked;
      lp->solisbasic = storedsolvals->solisbasic;

      /* solution values are stored only for LPs solved to optimality or unboundedness */
      assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL ||
         lp->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY ||
         lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT ||
         lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_ITERLIMIT ||
         lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_TIMELIMIT ||
         lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE ||
         lp->validsollp == -1);
   }
   /* no values available, mark LP as unsolved */
   else
   {
      lp->solved = FALSE;
      lp->validsollp = -1;

      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      lp->lpobjval = SCIP_INVALID;
      lp->primalfeasible = FALSE;
      lp->primalchecked = FALSE;
      lp->dualfeasible = FALSE;
      lp->dualchecked = FALSE;
      lp->solisbasic = FALSE;
      lp->validfarkaslp = -1;
   }

   /* intentionally keep storage space allocated */

   return SCIP_OKAY;
}

/** save current LP solution values stored in each column */
static
SCIP_RETCODE colStoreSolVals(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_COLSOLVALS* storedsolvals;

   assert(col != NULL);
   assert(blkmem != NULL);

   /* allocate memory for storage */
   if( col->storedsolvals == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &col->storedsolvals) );
   }
   storedsolvals = col->storedsolvals;

   /* store values */
   storedsolvals->primsol = col->primsol;
   storedsolvals->redcost = col->redcost;
   storedsolvals->basisstatus = col->basisstatus; /*lint !e641 !e732*/

   return SCIP_OKAY;
}

/** restore LP solution values in column */
static
SCIP_RETCODE colRestoreSolVals(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Longint          validlp,            /**< number of lp for which restored values are valid */
   SCIP_Bool             freebuffer          /**< should buffer for LP solution values be freed? */
   )
{
   SCIP_COLSOLVALS* storedsolvals;

   assert(col != NULL);
   assert(blkmem != NULL);

   /* if stored values are available, restore them */
   storedsolvals = col->storedsolvals;
   if( storedsolvals != NULL )
   {
      col->primsol = storedsolvals->primsol;
      col->redcost = storedsolvals->redcost;
      col->validredcostlp = validlp;
      col->basisstatus = storedsolvals->basisstatus; /*lint !e641 !e732*/

      /* we do not save the farkas coefficient, since it can be recomputed; thus, we invalidate it here */
      col->validfarkaslp = -1;
   }
   /* if the column was created after performing the storage (possibly during probing), we treat it as implicitly zero;
    * we make sure to invalidate the reduced cost and farkas coefficient, which are not available
    */
   else
   {
      col->primsol = 0.0;
      col->validredcostlp = -1;
      col->validfarkaslp = -1;
      col->basisstatus = SCIP_BASESTAT_ZERO; /*lint !e641*/
   }

   /* free memory */
   if( freebuffer )
   {
      BMSfreeBlockMemoryNull(blkmem, &col->storedsolvals);
      assert(col->storedsolvals == NULL);
   }

   return SCIP_OKAY;
}

/** save current LP solution values stored in each column */
static
SCIP_RETCODE rowStoreSolVals(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             infeasible          /**< is the solution infeasible? */
   )
{
   SCIP_ROWSOLVALS* storedsolvals;

   assert(row != NULL);
   assert(blkmem != NULL);

   /* allocate memory for storage */
   if( row->storedsolvals == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &row->storedsolvals) );
   }
   storedsolvals = row->storedsolvals;

   /* store values */
   if ( infeasible )
   {
      storedsolvals->dualsol = row->dualfarkas;
      storedsolvals->activity = SCIP_INVALID;
      storedsolvals->basisstatus = SCIP_BASESTAT_BASIC;  /*lint !e641*/
   }
   else
   {
      storedsolvals->dualsol = row->dualsol;
      storedsolvals->activity = row->activity;
      storedsolvals->basisstatus = row->basisstatus; /*lint !e641 !e732*/
   }

   return SCIP_OKAY;
}

/** restore LP solution values in row */
static
SCIP_RETCODE rowRestoreSolVals(
   SCIP_ROW*             row,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Longint          validlp,            /**< number of lp for which restored values are valid */
   SCIP_Bool             freebuffer,         /**< should buffer for LP solution values be freed? */
   SCIP_Bool             infeasible          /**< is the solution infeasible? */
   )
{
   SCIP_ROWSOLVALS* storedsolvals;

   assert(row != NULL);
   assert(blkmem != NULL);

   /* if stored values are available, restore them */
   storedsolvals = row->storedsolvals;
   if( storedsolvals != NULL )
   {
      if ( infeasible )
         row->dualfarkas = storedsolvals->dualsol;
      else
         row->dualsol = storedsolvals->dualsol;
      row->activity = storedsolvals->activity;
      row->validactivitylp = validlp;
      row->basisstatus = storedsolvals->basisstatus; /*lint !e641 !e732*/
   }
   /* if the row was created after performing the storage (possibly during probing), we treat it as basic;
    * we make sure to invalidate the reduced cost and farkas coefficient, which are not available
    */
   else
   {
      row->dualsol = 0.0;
      row->dualfarkas = 0.0;
      row->activity = SCIP_INVALID;
      row->validactivitylp = -1;
      row->basisstatus = SCIP_BASESTAT_BASIC; /*lint !e641*/
   }

   /* free memory */
   if( freebuffer )
   {
      BMSfreeBlockMemoryNull(blkmem, &row->storedsolvals);
      assert(row->storedsolvals == NULL);
   }

   return SCIP_OKAY;
}

/** ensures, that column array of row can store at least num entries */
SCIP_RETCODE SCIProwEnsureSize(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(row != NULL);
   assert(row->len <= row->size);

   if( num > row->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->cols, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->cols_index, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->vals, row->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &row->linkpos, row->size, newsize) );
      row->size = newsize;
   }
   assert(num <= row->size);

   return SCIP_OKAY;
}


#if 0 /* enable this to check the sortings within rows (for debugging, very slow!) */
static SCIP_Bool msgdisp_checkrow = FALSE;

static
void checkRow(
   SCIP_ROW*             row
   )
{
   int i;

   if( !msgdisp_checkrow )
   {
      printf("LP ROW CHECKING ACTIVATED! THIS IS VERY SLOW!\n");
      msgdisp_checkrow = TRUE;
   }

   /* validate sorting of LP part of row */
   if( row->lpcolssorted && row->nlpcols > 0)
   {
      assert(row->cols_index[0] == row->cols[0]->index);
      for( i = 1; i < row->nlpcols; ++i )
      {
         assert(row->cols_index[i] == row->cols[i]->index);
         assert(row->cols_index[i] >= row->cols_index[i-1]);
      }
   }

   /* validate sorting of non-LP part of row */
   if( row->nonlpcolssorted && row->len > row->nlpcols )
   {
      assert(row->cols_index[row->nlpcols] == row->cols[row->nlpcols]->index);
      for( i = row->nlpcols + 1; i < row->len; ++i )
      {
         assert(row->cols_index[i] == row->cols[i]->index);
         assert(row->cols_index[i] >= row->cols_index[i-1]);
      }
   }
}
#else
#define checkRow(row) /**/
#endif

#if 0 /* enable this to check norms of rows (for debugging, very slow!) */
static
void checkRowSqrnorm(
   SCIP_ROW*             row
   )
{
   SCIP_COL** cols;
   SCIP_Real sqrnorm;
   int c;

   cols = row->cols;
   assert(cols != NULL || row->len == 0);

   sqrnorm = 0.0;

   for( c = row->len - 1; c >= 0; --c )
   {
      if( cols[c]->lppos >= 0 )
         sqrnorm += SQR(row->vals[c]);
   }

   assert(ABS(sqrnorm - row->sqrnorm) < 1e-06 * MAX(1.0,sqrnorm));
}

static
void checkRowSumnorm(
   SCIP_ROW*             row
   )
{
   SCIP_COL** cols;
   SCIP_Real sumnorm;
   int c;

   cols = row->cols;
   assert(cols != NULL || row->len == 0);

   sumnorm = 0.0;

   for( c = row->len - 1; c >= 0; --c )
   {
      if( cols[c]->lppos >= 0 )
         sumnorm += REALABS(row->vals[c]);
   }

   assert(ABS(sumnorm - row->sumnorm) < 1e-06 * MAX(1.0,sumnorm));
}

static
void checkRowObjprod(
   SCIP_ROW*             row
   )
{
   SCIP_COL** cols;
   SCIP_Real objprod;
   int c;

   cols = row->cols;
   assert(cols != NULL || row->len == 0);

   objprod = 0.0;

   for( c = row->len - 1; c >= 0; --c )
   {
      if( cols[c]->lppos >= 0 )
         objprod += row->vals[c] * cols[c]->unchangedobj;
   }

   assert(ABS(objprod - row->objprod) < 1e-06 * MAX(1.0,objprod));
}
#else
#define checkRowSqrnorm(row) /**/
#define checkRowSumnorm(row) /**/
#define checkRowObjprod(row) /**/
#endif

/*
 * Local methods for pseudo and loose objective values
 */

/* recompute the loose objective value from scratch, if it was marked to be unreliable before */
static
void recomputeLooseObjectiveValue(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   SCIP_VAR** vars;
   SCIP_Real obj;
   int nvars;
   int v;

   assert(lp != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(!lp->looseobjvalid);

   vars = prob->vars;
   nvars = prob->nvars;
   lp->looseobjval = 0.0;

   /* iterate over all variables in the problem */
   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_LOOSE )
      {
         obj = SCIPvarGetObj(vars[v]);

         /* we are only interested in variables with a finite impact, because the infinity counters should be correct */
         if( SCIPsetIsPositive(set, obj) && !SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(vars[v])) )
            lp->looseobjval += obj * SCIPvarGetLbLocal(vars[v]);
         else if( SCIPsetIsNegative(set, obj) && !SCIPsetIsInfinity(set, SCIPvarGetUbLocal(vars[v])) )
            lp->looseobjval += obj * SCIPvarGetUbLocal(vars[v]);
      }
   }

   /* the recomputed value is reliable */
   lp->rellooseobjval = lp->looseobjval;
   lp->looseobjvalid = TRUE;
}

/* recompute the pseudo solution value from scratch, if it was marked to be unreliable before */
static
void recomputePseudoObjectiveValue(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(lp != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(!lp->pseudoobjvalid);

   vars = prob->vars;
   nvars = prob->nvars;
   lp->pseudoobjval = 0.0;

   /* iterate over all variables in the problem */
   for( v = 0; v < nvars; ++v )
   {
      /* we are only interested in variables with a finite impact, because the infinity counters should be correct */
      if( SCIPsetIsPositive(set, SCIPvarGetObj(vars[v])) &&
         !SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(vars[v])) )
      {
         lp->pseudoobjval += SCIPvarGetObj(vars[v]) * SCIPvarGetLbLocal(vars[v]);
      }
      else if( SCIPsetIsNegative(set, SCIPvarGetObj(vars[v])) &&
         !SCIPsetIsInfinity(set, SCIPvarGetUbLocal(vars[v])) )
      {
         lp->pseudoobjval += SCIPvarGetObj(vars[v]) * SCIPvarGetUbLocal(vars[v]);
      }
   }

   /* the recomputed value is reliable */
   lp->relpseudoobjval = lp->pseudoobjval;
   lp->pseudoobjvalid = TRUE;
}

/* recompute the global pseudo solution value from scratch, if it was marked to be unreliable before */
static
void recomputeGlbPseudoObjectiveValue(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(lp != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(!lp->glbpseudoobjvalid);

   vars = prob->vars;
   nvars = prob->nvars;
   lp->glbpseudoobjval = 0.0;

   /* iterate over all variables in the problem */
   for( v = 0; v < nvars; ++v )
   {
      /* we are only interested in variables with a finite impact, because the infinity counters should be correct */
      if( SCIPsetIsPositive(set, SCIPvarGetObj(vars[v])) &&
         !SCIPsetIsInfinity(set, -SCIPvarGetLbGlobal(vars[v])) )
      {
         lp->glbpseudoobjval += SCIPvarGetObj(vars[v]) * SCIPvarGetLbGlobal(vars[v]);
      }
      else if( SCIPsetIsNegative(set, SCIPvarGetObj(vars[v])) &&
         !SCIPsetIsInfinity(set, SCIPvarGetUbGlobal(vars[v])) )
      {
         lp->glbpseudoobjval += SCIPvarGetObj(vars[v]) * SCIPvarGetUbGlobal(vars[v]);
      }
   }

   /* the recomputed value is reliable */
   lp->relglbpseudoobjval = lp->glbpseudoobjval;
   lp->glbpseudoobjvalid = TRUE;
}

/** gets finite part of objective value of current LP that results from LOOSE variables only */
static
SCIP_Real getFiniteLooseObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert((lp->nloosevars > 0) || (lp->looseobjvalinf == 0 && lp->looseobjval == 0.0));
   assert(lp->flushed);
   assert(lp->looseobjvalinf == 0);

   /* recalculate the loose objective value, if needed */
   if( !lp->looseobjvalid )
      recomputeLooseObjectiveValue(lp, set, prob);

   return lp->looseobjval;
}

/** gets finite part of pseudo objective value of current LP */
static
SCIP_Real getFinitePseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   /* recalculate the pseudo objective value, if needed */
   if( !lp->pseudoobjvalid )
      recomputePseudoObjectiveValue(lp, set, prob);

   return lp->pseudoobjval;
}

/*
 * Sorting and searching rows and columns
 */


/** comparison method for sorting rows by non-decreasing index */
SCIP_DECL_SORTPTRCOMP(SCIProwComp)
{
   assert(elem1 != NULL);
   assert(elem2 != NULL);

   if( ((SCIP_ROW*)elem1)->index < ((SCIP_ROW*)elem2)->index )
      return -1;
   else if( ((SCIP_ROW*)elem1)->index > ((SCIP_ROW*)elem2)->index )
      return +1;
   else
   {
      assert(SCIProwGetIndex((SCIP_ROW*)(elem1)) == SCIProwGetIndex(((SCIP_ROW*)elem2)));
      return 0;
   }
}


/** sorts column entries of linked rows currently in the LP such that lower row indices precede higher ones */
static
void colSortLP(
   SCIP_COL*             col                 /**< column to be sorted */
   )
{
   int i;

   assert(col != NULL);

   /* check, if column is already sorted in the LP part */
   if( col->lprowssorted )
      return;

   /* sort coefficients */
   SCIPsortPtrRealInt((void**)col->rows, col->vals, col->linkpos, SCIProwComp, col->nlprows );

   /* update links */
   for( i = 0; i < col->nlprows; ++i )
   {
      if( col->linkpos[i] >= 0 )
      {
         assert(col->rows[i]->cols[col->linkpos[i]] == col);
         assert(col->rows[i]->linkpos[col->linkpos[i]] >= 0);
         col->rows[i]->linkpos[col->linkpos[i]] = i;
      }
   }

   col->lprowssorted = TRUE;
}

/** sorts column entries of unlinked rows or rows currently not in the LP such that lower row indices precede higher
 *  ones
 */
static
void colSortNonLP(
   SCIP_COL*             col                 /**< column to be sorted */
   )
{
   int i;

   assert(col != NULL);

   /* check, if column is already sorted in the non-LP part */
   if( col->nonlprowssorted )
      return;

   /* sort coefficients */
   SCIPsortPtrRealInt((void**)(&(col->rows[col->nlprows])), &(col->vals[col->nlprows]), &(col->linkpos[col->nlprows]), SCIProwComp, col->len - col->nlprows);

   /* update links */
   for( i = col->nlprows; i < col->len; ++i )
   {
      if( col->linkpos[i] >= 0 )
      {
         assert(col->rows[i]->cols[col->linkpos[i]] == col);
         assert(col->rows[i]->linkpos[col->linkpos[i]] >= 0);
         col->rows[i]->linkpos[col->linkpos[i]] = i;
      }
   }

   col->nonlprowssorted = TRUE;
}

/** sorts row entries of linked columns currently in the LP such that lower column indices precede higher ones */
static
void rowSortLP(
   SCIP_ROW*             row                 /**< row to be sorted */
   )
{
   int i;

   assert(row != NULL);

   /* check, if row is already sorted in the LP part, or if the sorting should be delayed */
   if( row->lpcolssorted || row->delaysort )
      return;

   /* sort coefficients */
   SCIPsortIntPtrIntReal(row->cols_index, (void**)row->cols, row->linkpos, row->vals, row->nlpcols);

   /* update links */
   for( i = 0; i < row->nlpcols; ++i )
   {
      if( row->linkpos[i] >= 0 )
      {
         assert(row->cols[i]->rows[row->linkpos[i]] == row);
         assert(row->cols[i]->linkpos[row->linkpos[i]] >= 0);
         row->cols[i]->linkpos[row->linkpos[i]] = i;
      }
   }

   row->lpcolssorted = TRUE;
}

/** sorts row entries of unlinked columns or columns currently not in the LP such that lower column indices precede
 *  higher ones
 */
static
void rowSortNonLP(
   SCIP_ROW*             row                 /**< row to be sorted */
   )
{
   int i;

   assert(row != NULL);

   checkRow(row);

   /* check, if row is already sorted in the non-LP part, or if the sorting should be delayed */
   if( row->nonlpcolssorted || row->delaysort )
      return;

   /* sort coefficients */
   SCIPsortIntPtrIntReal(&(row->cols_index[row->nlpcols]), (void**)(&(row->cols[row->nlpcols])), &(row->linkpos[row->nlpcols]), &(row->vals[row->nlpcols]), row->len - row->nlpcols);

   /* update links */
   for( i = row->nlpcols; i < row->len; ++i )
   {
      if( row->linkpos[i] >= 0 )
      {
         assert(row->cols[i]->rows[row->linkpos[i]] == row);
         assert(row->cols[i]->linkpos[row->linkpos[i]] >= 0);
         row->cols[i]->linkpos[row->linkpos[i]] = i;
      }
   }

   checkRow(row);

   row->nonlpcolssorted = TRUE;
}

/** searches coefficient in part of the column, returns position in col vector or -1 if not found */
static
int colSearchCoefPart(
   SCIP_COL*             col,                /**< column to be searched in */
   const SCIP_ROW*       row,                /**< coefficient to be searched for */
   int                   minpos,             /**< first position of search range */
   int                   maxpos              /**< last position of search range */
   )
{
   int pos;
   int idx;
   int searchidx;

   assert(col != NULL);
   assert(row != NULL);

   /* binary search */
   searchidx = row->index;
   while(minpos <= maxpos)
   {
      pos = (minpos + maxpos)/2;
      assert(0 <= pos && pos < col->len);
      assert(col->rows[pos] != NULL);
      assert((pos < col->nlprows) == (col->rows[pos]->lppos >= 0 && col->linkpos[pos] >= 0));
      idx = col->rows[pos]->index;
      if( searchidx == idx )
         return pos;
      else if( searchidx < idx )
         maxpos = pos-1;
      else
         minpos = pos+1;
   }

   return -1;
}

/** searches coefficient in column, returns position in col vector or -1 if not found */
static
int colSearchCoef(
   SCIP_COL*             col,                /**< column to be searched in */
   const SCIP_ROW*       row                 /**< coefficient to be searched for */
   )
{
   int pos;

   assert(col != NULL);
   assert(row != NULL);

   pos = -1;

   /* search in the linked LP rows */
   if( row->lppos >= 0 )
   {
      /* column has to be sorted, such that binary search works */
      colSortLP(col);
      assert(col->lprowssorted);

      pos = colSearchCoefPart(col, row, 0, col->nlprows-1);
      if( pos >= 0 )
         return pos;
   }

   /* search in the non-LP/unlinked rows */
   if( row->lppos == -1 || col->nunlinked > 0 )
   {
      /* column has to be sorted, such that binary search works */
      colSortNonLP(col);
      assert(col->nonlprowssorted);

      pos = colSearchCoefPart(col, row, col->nlprows, col->len-1);
   }

   return pos;
}

/** searches coefficient in part of the row, returns position in col vector or -1 if not found */
static
int rowSearchCoefPart(
   SCIP_ROW*             row,                /**< row to be searched in */
   const SCIP_COL*       col,                /**< coefficient to be searched for */
   int                   minpos,             /**< first position of search range */
   int                   maxpos              /**< last position of search range */
   )
{
   int pos;
   int idx;
   int searchidx;

   assert(row != NULL);
   assert(col != NULL);

   /* binary search */
   searchidx = col->index;
   while(minpos <= maxpos)
   {
      pos = (minpos + maxpos)/2;
      assert(0 <= pos && pos < row->len);
      assert(row->cols[pos] != NULL);
      assert((pos < row->nlpcols) == (row->cols[pos]->lppos >= 0 && row->linkpos[pos] >= 0));
      assert(row->cols_index[pos] == row->cols[pos]->index);
      idx = row->cols_index[pos];
      if( searchidx == idx )
         return pos;
      else if( searchidx < idx )
         maxpos = pos-1;
      else
         minpos = pos+1;
   }

   return -1;
}

/** searches coefficient in row, returns position in row vector or -1 if not found;
 *  if the sorting of the row is delayed, returns -1
 */
static
int rowSearchCoef(
   SCIP_ROW*             row,                /**< row to be searched in */
   const SCIP_COL*       col                 /**< coefficient to be searched for */
   )
{
   int pos;

   assert(row != NULL);
   assert(col != NULL);

   if( row->delaysort )
      return -1;

   pos = -1;

   /* search in the linked LP columns */
   if( col->lppos >= 0 )
   {
      /* row has to be sorted, such that binary search works */
      rowSortLP(row);
      assert(row->lpcolssorted);

      pos = rowSearchCoefPart(row, col, 0, row->nlpcols-1);
   }

   /* search in the non-LP/unlinked columns */
   if( pos == -1 && (col->lppos == -1 || row->nunlinked > 0) )
   {
      /* row has to be sorted, such that binary search works */
      rowSortNonLP(row);
      assert(row->nonlpcolssorted);

      pos = rowSearchCoefPart(row, col, row->nlpcols, row->len-1);
   }

#ifndef NDEBUG
   /* validate result */
   assert(-1 <= pos && pos < row->len);
   if( pos >= 0 )
      assert(row->cols[pos] == col);
   else
   {
      int i;
      for( i = 0; i < row->len; ++i )
         assert(row->cols[i] != col);
   }
#endif

   return pos;
}

/** moves a coefficient in a column to a different place, and updates all corresponding data structures */
static
void colMoveCoef(
   SCIP_COL*             col,                /**< LP column */
   int                   oldpos,             /**< old position of coefficient */
   int                   newpos              /**< new position of coefficient */
   )
{
   assert(col != NULL);
   assert(0 <= oldpos && oldpos < col->len);
   assert(0 <= newpos && newpos < col->len);
   assert(col->rows[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   col->rows[newpos] = col->rows[oldpos];
   col->vals[newpos] = col->vals[oldpos];
   col->linkpos[newpos] = col->linkpos[oldpos];

   /* update link position in row */
   if( col->linkpos[newpos] >= 0 )
   {
      assert(col->rows[newpos]->cols[col->linkpos[newpos]] == col);
      assert(col->rows[newpos]->linkpos[col->linkpos[newpos]] == oldpos);

      col->rows[newpos]->linkpos[col->linkpos[newpos]] = newpos;
   }

   /* update sorted flags */
   if( col->rows[newpos]->lppos >= 0 && col->linkpos[newpos] >= 0 )
      col->lprowssorted = FALSE;
   else
      col->nonlprowssorted = FALSE;
}

/** swaps two coefficients in a column, and updates all corresponding data structures */
static
void colSwapCoefs(
   SCIP_COL*             col,                /**< LP column */
   int                   pos1,               /**< position of first coefficient */
   int                   pos2                /**< position of second coefficient */
   )
{
   SCIP_ROW* tmprow;
   SCIP_Real tmpval;
   int tmplinkpos;

   assert(col != NULL);
   assert(0 <= pos1 && pos1 < col->len);
   assert(0 <= pos2 && pos2 < col->len);
   assert(col->rows[pos1] != NULL);

   if( pos1 == pos2 )
      return;

   /* swap coefficients */
   tmprow = col->rows[pos2];
   tmpval = col->vals[pos2];
   tmplinkpos = col->linkpos[pos2];

   col->rows[pos2] = col->rows[pos1];
   col->vals[pos2] = col->vals[pos1];
   col->linkpos[pos2] = col->linkpos[pos1];

   col->rows[pos1] = tmprow;
   col->vals[pos1] = tmpval;
   col->linkpos[pos1] = tmplinkpos;

   /* update link position in rows */
   if( col->linkpos[pos1] >= 0 )
   {
      assert(col->rows[pos1]->cols[col->linkpos[pos1]] == col);
      assert(col->rows[pos1]->linkpos[col->linkpos[pos1]] == pos2);

      col->rows[pos1]->linkpos[col->linkpos[pos1]] = pos1;
   }
   if( col->linkpos[pos2] >= 0 )
   {
      assert(col->rows[pos2]->cols[col->linkpos[pos2]] == col);
      assert(col->rows[pos2]->linkpos[col->linkpos[pos2]] == pos1);

      col->rows[pos2]->linkpos[col->linkpos[pos2]] = pos2;
   }

   /* update sorted flags */
   if( col->rows[pos1]->lppos >= 0 && col->linkpos[pos1] >= 0 )
      col->lprowssorted = FALSE;
   else
      col->nonlprowssorted = FALSE;
   if( col->rows[pos2]->lppos >= 0 && col->linkpos[pos2] >= 0 )
      col->lprowssorted = FALSE;
   else
      col->nonlprowssorted = FALSE;
}

/** moves a coefficient in a row to a different place, and updates all corresponding data structures */
static
void rowMoveCoef(
   SCIP_ROW*             row,                /**< LP row */
   int                   oldpos,             /**< old position of coefficient */
   int                   newpos              /**< new position of coefficient */
   )
{
   assert(row != NULL);
   assert(0 <= oldpos && oldpos < row->len);
   assert(0 <= newpos && newpos < row->len);
   assert(row->cols[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   row->cols[newpos] = row->cols[oldpos];
   row->cols_index[newpos] = row->cols_index[oldpos];
   row->vals[newpos] = row->vals[oldpos];
   row->linkpos[newpos] = row->linkpos[oldpos];

   /* update link position in column */
   if( row->linkpos[newpos] >= 0 )
   {
      assert(row->cols[newpos]->rows[row->linkpos[newpos]] == row);
      assert(row->cols[newpos]->linkpos[row->linkpos[newpos]] == oldpos);

      row->cols[newpos]->linkpos[row->linkpos[newpos]] = newpos;
   }

   /* update sorted flags */
   if( row->cols[newpos]->lppos >= 0 && row->linkpos[newpos] >= 0 )
      row->lpcolssorted = FALSE;
   else
      row->nonlpcolssorted = FALSE;
}

/** swaps two coefficients in a row, and updates all corresponding data structures */
static
void rowSwapCoefs(
   SCIP_ROW*             row,                /**< LP row */
   int                   pos1,               /**< position of first coefficient */
   int                   pos2                /**< position of second coefficient */
   )
{
   SCIP_COL* tmpcol;
   SCIP_Real tmpval;
   int tmpindex;
   int tmplinkpos;

   assert(row != NULL);
   assert(0 <= pos1 && pos1 < row->len);
   assert(0 <= pos2 && pos2 < row->len);
   assert(row->cols[pos1] != NULL);
   assert(row->cols[pos1]->index == row->cols_index[pos1]);

   if( pos1 == pos2 )
      return;

   /* swap coefficients */
   tmpcol = row->cols[pos2];
   tmpindex = row->cols_index[pos2];
   tmpval = row->vals[pos2];
   tmplinkpos = row->linkpos[pos2];

   row->cols[pos2] = row->cols[pos1];
   row->cols_index[pos2] = row->cols_index[pos1];
   row->vals[pos2] = row->vals[pos1];
   row->linkpos[pos2] = row->linkpos[pos1];

   row->cols[pos1] = tmpcol;
   row->cols_index[pos1] = tmpindex;
   row->vals[pos1] = tmpval;
   row->linkpos[pos1] = tmplinkpos;

   /* update link position in columns */
   if( row->linkpos[pos1] >= 0 )
   {
      assert(row->cols[pos1]->rows[row->linkpos[pos1]] == row);
      assert(row->cols[pos1]->linkpos[row->linkpos[pos1]] == pos2);

      row->cols[pos1]->linkpos[row->linkpos[pos1]] = pos1;
   }
   if( row->linkpos[pos2] >= 0 )
   {
      assert(row->cols[pos2]->rows[row->linkpos[pos2]] == row);
      assert(row->cols[pos2]->linkpos[row->linkpos[pos2]] == pos1);

      row->cols[pos2]->linkpos[row->linkpos[pos2]] = pos2;
   }

   /* update sorted flags */
   if( row->cols[pos1]->lppos >= 0 && row->linkpos[pos1] >= 0 )
      row->lpcolssorted = FALSE;
   else
      row->nonlpcolssorted = FALSE;
   if( row->cols[pos2]->lppos >= 0 && row->linkpos[pos2] >= 0 )
      row->lpcolssorted = FALSE;
   else
      row->nonlpcolssorted = FALSE;
}

/** issues a ROWCOEFCHANGED event on the given row */
static
SCIP_RETCODE rowEventCoefChanged(
   SCIP_ROW*             row,                /**< row which coefficient has changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_COL*             col,                /**< the column which coefficient has changed */
   SCIP_Real             oldval,             /**< old value of the coefficient */
   SCIP_Real             newval              /**< new value of the coefficient */
   )
{
   assert(row != NULL);
   assert(row->eventfilter != NULL);
   assert(col != NULL);

   /* check, if the row is being tracked for coefficient changes
    * if so, issue ROWCOEFCHANGED event
    */
   if( (row->eventfilter->len > 0 && (row->eventfilter->eventmask & SCIP_EVENTTYPE_ROWCOEFCHANGED) != 0) )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowCoefChanged(&event, blkmem, row, col, oldval, newval) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, row->eventfilter, &event) );
   }

   return SCIP_OKAY;
}

/** issues a ROWCONSTCHANGED event on the given row */
static
SCIP_RETCODE rowEventConstantChanged(
   SCIP_ROW*             row,                /**< row which coefficient has changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             oldval,             /**< old value of the constant */
   SCIP_Real             newval              /**< new value of the constant */
   )
{
   assert(row != NULL);
   assert(row->eventfilter != NULL);

   /* check, if the row is being tracked for coefficient changes
    * if so, issue ROWCONSTCHANGED event
    */
   if( (row->eventfilter->len > 0 && (row->eventfilter->eventmask & SCIP_EVENTTYPE_ROWCONSTCHANGED) != 0) )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowConstChanged(&event, blkmem, row, oldval, newval) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, row->eventfilter, &event) );
   }

   return SCIP_OKAY;
}

/** issues a ROWSIDECHANGED event on the given row */
static
SCIP_RETCODE rowEventSideChanged(
   SCIP_ROW*             row,                /**< row which coefficient has changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_SIDETYPE         side,               /**< the side that has changed */
   SCIP_Real             oldval,             /**< old value of side */
   SCIP_Real             newval              /**< new value of side */
   )
{
   assert(row != NULL);
   assert(row->eventfilter != NULL);

   /* check, if the row is being tracked for coefficient changes
    * if so, issue ROWSIDECHANGED event
    */
   if( (row->eventfilter->len > 0 && (row->eventfilter->eventmask & SCIP_EVENTTYPE_ROWSIDECHANGED) != 0) )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowSideChanged(&event, blkmem, row, side, oldval, newval) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, row->eventfilter, &event) );
   }

   return SCIP_OKAY;
}

#if 0 /* enable this to check links between columns and rows in LP data structure (for debugging, very slow!) */

#ifdef NDEBUG
#define ASSERT(x) do { if( !(x) ) abort(); } while( FALSE )
#else
#define ASSERT(x) assert(x)
#endif

static SCIP_Bool msgdisp_checklinks = FALSE;


static
void checkLinks(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_COL* col;
   SCIP_ROW* row;
   int i;
   int j;

   ASSERT(lp != NULL);

   if( !msgdisp_checklinks )
   {
      printf("LP LINK CHECKING ACTIVATED! THIS IS VERY SLOW!\n");
      msgdisp_checklinks = TRUE;
   }

   for( i = 0; i < lp->ncols; ++i )
   {
      col = lp->cols[i];
      ASSERT(col != NULL);
      ASSERT(!lp->flushed || col->lppos >= 0 || col->primsol == 0.0);
      ASSERT(!lp->flushed || col->lppos >= 0 || col->farkascoef == 0.0);
      ASSERT(col->nlprows <= col->len);
      ASSERT(col->lppos == -1 || col->lppos >= lp->lpifirstchgcol || col->nunlinked == 0);

      for( j = 0; j < col->len; ++j )
      {
         row = col->rows[j];
         ASSERT(row != NULL);
         ASSERT(!lp->flushed || col->lppos == -1 || col->linkpos[j] >= 0);
         ASSERT(col->linkpos[j] == -1 || row->cols[col->linkpos[j]] == col);
         ASSERT(col->linkpos[j] == -1 || EPSEQ(row->vals[col->linkpos[j]], col->vals[j], 1e-6));
         ASSERT((j < col->nlprows) == (col->linkpos[j] >= 0 && row->lppos >= 0));
      }
   }

   for( i = 0; i < lp->nrows; ++i )
   {
      row = lp->rows[i];
      ASSERT(row != NULL);
      ASSERT(!lp->flushed || row->lppos >= 0 || row->dualsol == 0.0);
      ASSERT(!lp->flushed || row->lppos >= 0 || row->dualfarkas == 0.0);
      ASSERT(row->nlpcols <= row->len);
      ASSERT(row->lppos == -1 || row->lppos >= lp->lpifirstchgrow || row->nunlinked == 0);

      for( j = 0; j < row->len; ++j )
      {
         col = row->cols[j];
         ASSERT(col != NULL);
         ASSERT(!lp->flushed || row->lppos == -1 || row->linkpos[j] >= 0);
         ASSERT(row->linkpos[j] == -1 || col->rows[row->linkpos[j]] == row);
         ASSERT(row->linkpos[j] == -1 || EPSEQ(col->vals[row->linkpos[j]], row->vals[j], 1e-6));
         ASSERT((j < row->nlpcols) == (row->linkpos[j] >= 0 && col->lppos >= 0));
      }
   }
}

#undef ASSERT

#else
#define checkLinks(lp) /**/
#endif

/*
 * Changing announcements
 */

/** announces, that the given coefficient in the constraint matrix changed */
static
void coefChanged(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_COL*             col,                /**< LP col */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(row != NULL);
   assert(col != NULL);
   assert(lp != NULL);

   if( row->lpipos >= 0 && col->lpipos >= 0 )
   {
      assert(row->lpipos < lp->nlpirows);
      assert(col->lpipos < lp->nlpicols);

      /* we have to remember the change only in the row or in the column, 
       * because the readdition of one vector would change the other automatically.
       */
      if( row->lpipos >= lp->lpifirstchgrow )
         row->coefchanged = TRUE;
      else if( col->lpipos >= lp->lpifirstchgcol )
         col->coefchanged = TRUE;
      else if( lp->lpifirstchgrow - row->lpipos <= lp->lpifirstchgcol - col->lpipos )
      {
         row->coefchanged = TRUE;
         lp->lpifirstchgrow = row->lpipos;
      }
      else
      {
         col->coefchanged = TRUE;
         lp->lpifirstchgcol = col->lpipos;
      }

      /* mark the current LP unflushed */
      lp->flushed = FALSE;
   }

   row->pseudoactivity = SCIP_INVALID;
   row->minactivity = SCIP_INVALID;
   row->maxactivity = SCIP_INVALID;
   row->validpsactivitydomchg = -1;
   row->validactivitybdsdomchg = -1;
}



/*
 * local column changing methods
 */

/* forward declaration for colAddCoef() */
static
SCIP_RETCODE rowAddCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val,                /**< value of coefficient */
   int                   linkpos             /**< position of row in the column's row array, or -1 */
   );

/** adds a previously non existing coefficient to an LP column */
static
SCIP_RETCODE colAddCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val,                /**< value of coefficient */
   int                   linkpos             /**< position of column in the row's col array, or -1 */
   )
{
   int pos;

   assert(blkmem != NULL);
   assert(col != NULL);
   assert(col->nlprows <= col->len);
   assert(col->var != NULL);
   assert(row != NULL);
   assert(!SCIPsetIsZero(set, val));
   /*assert(colSearchCoef(col, row) == -1);*/ /* this assert would lead to slight differences in the solution process */

   SCIP_CALL( colEnsureSize(col, blkmem, set, col->len+1) );
   assert(col->rows != NULL);
   assert(col->vals != NULL);
   assert(col->linkpos != NULL);

   pos = col->len;
   col->len++;

   /* if the row is in current LP and is linked to the column, we have to insert it at the end of the linked LP rows
    * part of the column's arrays
    */
   if( row->lppos >= 0 && linkpos >= 0 )
   {
      /* move the first non-LP/not linked row to the end */
      if( col->nlprows < pos )
      {
         colMoveCoef(col, col->nlprows, pos);
         pos = col->nlprows;
      }
      col->nlprows++;
   }

   /* in case the coefficient is integral w.r.t. numerics we explicitly round the coefficient to an integral value */
   val = SCIPsetIsIntegral(set, val) ? SCIPsetRound(set, val) : val;

   /* insert the row at the correct position and update the links */
   col->rows[pos] = row;
   col->vals[pos] = val;
   col->linkpos[pos] = linkpos;
   if( linkpos == -1 )
   {
      col->nunlinked++;

      /* if the column is in current LP, we have to link it to the row, because otherwise, the primal information
       * of the row is not complete
       */
      if( col->lppos >= 0 )
      {
         /* this call might swap the current row with the first non-LP/not linked row, s.t. insertion position
          * has to be updated
          */
         SCIP_CALL( rowAddCoef(row, blkmem, set, eventqueue, lp, col, val, pos) );
         if( row->lppos >= 0 )
            pos = col->nlprows-1;
         linkpos = col->linkpos[pos];

         assert(0 <= linkpos && linkpos < row->len);
         assert(row->cols[linkpos] == col);
         assert(col->rows[pos] == row);
         assert(col->rows[pos]->cols[col->linkpos[pos]] == col);
         assert(col->rows[pos]->linkpos[col->linkpos[pos]] == pos);
      }
   }
   else
   {
      assert(row->linkpos[linkpos] == -1);
      assert(row->nunlinked > 0);
      row->linkpos[linkpos] = pos;
      row->nunlinked--;

      /* if the column is in current LP, now both conditions, row->cols[linkpos]->lppos >= 0 and row->linkpos[linkpos] >= 0
       * hold, so we have to move the column to the linked LP-cols part of the row's cols array
       */
      if( col->lppos >= 0 )
      {
         row->nlpcols++;
         rowSwapCoefs(row, linkpos, row->nlpcols-1);

         /* if no swap was necessary, mark nonlpcols to be unsorted */
         if( linkpos == row->nlpcols-1 )
            row->lpcolssorted = FALSE;
      }
   }

   /* update the sorted flags */
   if( row->lppos >= 0 && linkpos >= 0 )
   {
      assert(col->nlprows >= 1);
      assert(col->rows[col->nlprows-1] == row);
      if( col->nlprows > 1 )
         col->lprowssorted = col->lprowssorted && (col->rows[col->nlprows-2]->index < row->index);
   }
   else
   {
      assert(col->len - col->nlprows >= 1);
      assert(col->rows[col->len-1] == row);
      if( col->len - col->nlprows > 1 )
         col->nonlprowssorted = col->nonlprowssorted && (col->rows[col->len-2]->index < row->index);
   }

   coefChanged(row, col, lp);

   SCIPsetDebugMsg(set, "added coefficient %g * <%s> at position %d (%d/%d) to column <%s> (nunlinked=%d)\n",
      val, row->name, pos, col->nlprows, col->len, SCIPvarGetName(col->var), col->nunlinked);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from column */
static
SCIP_RETCODE colDelCoefPos(
   SCIP_COL*             col,                /**< column to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   pos                 /**< position in column vector to delete */
   )
{
   SCIP_ROW* row;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] != NULL);
   assert(col->linkpos[pos] == -1 || col->rows[pos]->cols[col->linkpos[pos]] == col);
   assert((pos < col->nlprows) == (col->linkpos[pos] >= 0 && col->rows[pos]->lppos >= 0));

   row = col->rows[pos];
   assert((row->lppos >= 0) == (pos < col->nlprows));

   /*SCIPsetDebugMsg(set, "deleting coefficient %g * <%s> at position %d from column <%s>\n",
     col->vals[pos], row->name, pos, SCIPvarGetName(col->var));*/

   if( col->linkpos[pos] == -1 )
      col->nunlinked--;

   /* if row is a linked LP row, move last linked LP coefficient to position of empty slot (deleted coefficient) */
   if( pos < col->nlprows )
   {
      colMoveCoef(col, col->nlprows-1, pos);
      col->nlprows--;
      pos = col->nlprows;
   }

   /* move last coefficient to position of empty slot */
   colMoveCoef(col, col->len-1, pos);
   col->len--;

   coefChanged(row, col, lp);

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP column */
static
SCIP_RETCODE colChgCoefPos(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   pos,                /**< position in column vector to change */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] != NULL);
   assert(col->linkpos[pos] == -1 || col->rows[pos]->cols[col->linkpos[pos]] == col);

   /*debugMsg(scip, "changing coefficient %g * <%s> at position %d of column <%s> to %g\n",
     col->vals[pos], col->rows[pos]->name, pos, SCIPvarGetName(col->var), val);*/

   /* in case the coefficient is integral w.r.t. numerics we explicitly round the coefficient to an integral value */
   val = SCIPsetIsIntegral(set, val) ? SCIPsetRound(set, val) : val;

   if( SCIPsetIsZero(set, val) )
   {
      /* delete existing coefficient */
      SCIP_CALL( colDelCoefPos(col, set, lp, pos) );
   }
   else if( !SCIPsetIsEQ(set, col->vals[pos], val) )
   {
      /* change existing coefficient */
      col->vals[pos] = val;
      coefChanged(col->rows[pos], col, lp);
   }

   return SCIP_OKAY;
}




/*
 * local row changing methods
 */

/** update row norms after addition of coefficient */
static
void rowAddNorms(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col,                /**< column of added coefficient */
   SCIP_Real             val,                /**< value of added coefficient */
   SCIP_Bool             updateidxvals       /**< update min/max idx and min/max val? */
   )
{
   SCIP_Real absval;

   assert(row != NULL);
   assert(row->nummaxval >= 0);
   assert(row->numminval >= 0);
   assert(set != NULL);
   assert(col != NULL);

   absval = REALABS(val);
   assert(!SCIPsetIsZero(set, absval));

   /* Euclidean norm, sum norm, and objective function scalar product only take LP columns into account */
   if( col->lppos >= 0 )
   {
      /* update squared Euclidean norm and sum norm */
      row->sqrnorm += SQR(absval);
      row->sumnorm += absval;

      /* update objective function scalar product */
      row->objprod += val * col->unchangedobj;
   }

   if( updateidxvals )
   {
      /* update min/maxidx */
      row->minidx = MIN(row->minidx, col->index);
      row->maxidx = MAX(row->maxidx, col->index);

      /* update maximal and minimal non-zero value */
      if( row->nummaxval > 0 )
      {
         if( SCIPsetIsGT(set, absval, row->maxval) )
         {
            row->maxval = absval;
            row->nummaxval = 1;
         }
         else if( SCIPsetIsGE(set, absval, row->maxval) )
         {
            /* make sure the maxval is always exactly the same */
            row->maxval = MAX(absval, row->maxval);
            row->nummaxval++;
         }
      }
      if( row->numminval > 0 )
      {
         if( SCIPsetIsLT(set, absval, row->minval) )
         {
            row->minval = absval;
            row->numminval = 1;
         }
         else if( SCIPsetIsLE(set, absval, row->minval) )
         {
            /* make sure the minval is always exactly the same */
            row->minval = MIN(absval, row->minval);
            row->numminval++;
         }
      }
   }
   else
   {
      assert(row->minidx <= col->index);
      assert(row->maxidx >= col->index);
      assert(row->numminval <= 0 || absval >= row->minval);
      assert(row->nummaxval <= 0 || absval <= row->maxval);
   }
}

/** update row norms after deletion of coefficient */
static
void rowDelNorms(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col,                /**< column of deleted coefficient */
   SCIP_Real             val,                /**< value of deleted coefficient */
   SCIP_Bool             forcenormupdate,    /**< should the norms be updated even if lppos of column is -1? */
   SCIP_Bool             updateindex,        /**< should the minimal/maximal column index of row be updated? */
   SCIP_Bool             updateval           /**< should the minimal/maximal value of row be updated? */
   )
{
   SCIP_Real absval;

   assert(row != NULL);
   assert(row->nummaxval >= 0);
   assert(row->numminval >= 0);
   assert(set != NULL);
   assert(col != NULL);

   absval = REALABS(val);
   assert(!SCIPsetIsZero(set, absval));
   assert(row->nummaxval == 0 || row->maxval >= absval);
   assert(row->numminval == 0 || row->minval <= absval);

   /* update min/maxidx validity */
   if( updateindex && (col->index == row->minidx || col->index == row->maxidx) )
      row->validminmaxidx = FALSE;

   /* Euclidean norm, sum norm, and objective function scalar product only take LP columns into account */
   if( forcenormupdate || col->lppos >= 0 )
   {
      /* update squared Euclidean norm and sum norm */
      row->sqrnorm -= SQR(absval);
      row->sqrnorm = MAX(row->sqrnorm, 0.0);
      row->sumnorm -= absval;
      row->sumnorm = MAX(row->sumnorm, 0.0);

      /* update objective function scalar product */
      row->objprod -= val * col->unchangedobj;
   }

   if( updateval )
   {
      /* update maximal and minimal non-zero value */
      if( row->nummaxval > 0 )
      {
         if( SCIPsetIsGE(set, absval, row->maxval) )
            row->nummaxval--;
      }
      if( row->numminval > 0 )
      {
         if( SCIPsetIsLE(set, absval, row->minval) )
            row->numminval--;
      }
   }
}

/** adds a previously non existing coefficient to an LP row */
static
SCIP_RETCODE rowAddCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val,                /**< value of coefficient */
   int                   linkpos             /**< position of row in the column's row array, or -1 */
   )
{
   int pos;

   assert(row != NULL);
   assert(row->nlpcols <= row->len);
   assert(blkmem != NULL);
   assert(col != NULL);
   assert(col->var != NULL);
   assert(col->var_probindex == SCIPvarGetProbindex(col->var));
   assert(!SCIPsetIsZero(set, val));
   /*assert(rowSearchCoef(row, col) == -1);*/ /* this assert would lead to slight differences in the solution process */

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot add a coefficient to the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIProwEnsureSize(row, blkmem, set, row->len+1) );
   assert(row->cols != NULL);
   assert(row->vals != NULL);

   pos = row->len;
   row->len++;

   /* if the column is in current LP and is linked to the row, we have to insert it at the end of the linked LP columns
    * part of the row's arrays
    */
   if( col->lppos >= 0 && linkpos >= 0 )
   {
      /* move the first non-LP/not linked column to the end */
      if( row->nlpcols < pos )
      {
         rowMoveCoef(row, row->nlpcols, pos);
         pos = row->nlpcols;
      }
      row->nlpcols++;
   }

   /* in case the coefficient is integral w.r.t. numerics we explicitly round the coefficient to an integral value */
   val = SCIPsetIsIntegral(set, val) ? SCIPsetRound(set, val) : val;

   /* insert the column at the correct position and update the links */
   row->cols[pos] = col;
   row->cols_index[pos] = col->index;
   row->vals[pos] = val;
   row->linkpos[pos] = linkpos;
   row->integral = row->integral && SCIPcolIsIntegral(col) && SCIPsetIsIntegral(set, val);
   if( linkpos == -1 )
   {
      row->nunlinked++;

      /* if the row is in current LP, we have to link it to the column, because otherwise, the dual information
       * of the column is not complete
       */
      if( row->lppos >= 0 )
      {
         /* this call might swap the current column with the first non-LP/not linked column, s.t. insertion position
          * has to be updated
          */
         SCIP_CALL( colAddCoef(col, blkmem, set, eventqueue, lp, row, val, pos) );
         if( col->lppos >= 0 )
            pos = row->nlpcols-1;
         linkpos = row->linkpos[pos];

         assert(0 <= linkpos && linkpos < col->len);
         assert(col->rows[linkpos] == row);
         assert(row->cols[pos] == col);
         assert(row->cols[pos]->rows[row->linkpos[pos]] == row);
         assert(row->cols[pos]->linkpos[row->linkpos[pos]] == pos);
      }
   }
   else
   {
      assert(col->linkpos[linkpos] == -1);
      assert(col->nunlinked > 0);
      col->linkpos[linkpos] = pos;
      col->nunlinked--;

      /* if the row is in current LP, now both conditions, col->rows[linkpos]->lppos >= 0 and col->linkpos[linkpos] >= 0
       * hold, so we have to move the row to the linked LP-rows part of the column's rows array
       */
      if( row->lppos >= 0 )
      {
         col->nlprows++;
         colSwapCoefs(col, linkpos, col->nlprows-1);

         /* if no swap was necessary, mark lprows to be unsorted */
         if( linkpos == col->nlprows-1 )
            col->lprowssorted = FALSE;
      }
   }

   /* update the sorted flags */
   if( col->lppos >= 0 && linkpos >= 0 )
   {
      assert(row->nlpcols >= 1);
      assert(row->cols[row->nlpcols-1] == col);
      if( row->nlpcols > 1 )
      {
         assert(row->cols_index[row->nlpcols-2] == row->cols[row->nlpcols-2]->index);
         row->lpcolssorted = row->lpcolssorted && (row->cols_index[row->nlpcols-2] < col->index);
      }
   }
   else
   {
      assert(row->len - row->nlpcols >= 1);
      assert(row->cols[row->len-1] == col);
      if( row->len - row->nlpcols > 1 )
      {
         assert(row->cols_index[row->len-2] == row->cols[row->len-2]->index);
         row->nonlpcolssorted = row->nonlpcolssorted && (row->cols_index[row->len-2] < col->index);
      }
   }

   /* update row norm */
   rowAddNorms(row, set, col, val, TRUE);

   coefChanged(row, col, lp);

   SCIPsetDebugMsg(set, "added coefficient %g * <%s> at position %d (%d/%d) to row <%s> (nunlinked=%d)\n",
      val, SCIPvarGetName(col->var), pos, row->nlpcols, row->len, row->name, row->nunlinked);

   /* issue row coefficient changed event */
   SCIP_CALL( rowEventCoefChanged(row, blkmem, set, eventqueue, col, 0.0, val) );

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
SCIP_RETCODE rowDelCoefPos(
   SCIP_ROW*             row,                /**< row to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   pos                 /**< position in row vector to delete */
   )
{
   SCIP_COL* col;
   SCIP_Real val;

   assert(row != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < row->len);
   assert(row->cols[pos] != NULL);
   assert((pos < row->nlpcols) == (row->linkpos[pos] >= 0 && row->cols[pos]->lppos >= 0));

   col = row->cols[pos];
   val = row->vals[pos];
   assert((pos < row->nlpcols) == (col->lppos >= 0 && row->linkpos[pos] >= 0));

   /*SCIPsetDebugMsg(set, "deleting coefficient %g * <%s> at position %d from row <%s>\n",
     val, SCIPvarGetName(col->var), pos, row->name);*/

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot delete a coefficient from the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   if( row->linkpos[pos] == -1 )
      row->nunlinked--;

   /* if column is a linked LP column, move last linked LP coefficient to position of empty slot (deleted coefficient) */
   if( pos < row->nlpcols )
   {
      rowMoveCoef(row, row->nlpcols-1, pos);
      assert(!row->lpcolssorted);
      row->nlpcols--;
      pos = row->nlpcols;
   }

   /* move last coefficient to position of empty slot */
   rowMoveCoef(row, row->len-1, pos);
   row->len--;

   /* update norms */
   rowDelNorms(row, set, col, val, FALSE, TRUE, TRUE);

   coefChanged(row, col, lp);

   /* issue row coefficient changed event */
   SCIP_CALL( rowEventCoefChanged(row, blkmem, set, eventqueue, col, val, 0.0) );

   return SCIP_OKAY;
}

/** changes a coefficient at given position of an LP row */
static
SCIP_RETCODE rowChgCoefPos(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   int                   pos,                /**< position in row vector to change */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   SCIP_COL* col;

   assert(row != NULL);
   assert(0 <= pos && pos < row->len);

   /*SCIPsetDebugMsg(set, "changing coefficient %g * <%s> at position %d of row <%s> to %g\n",
     row->vals[pos], SCIPvarGetName(row->cols[pos]->var), pos, row->name, val);*/

   if( row->nlocks > 0 )
   {
      SCIPerrorMessage("cannot change a coefficient of the locked unmodifiable row <%s>\n", row->name);
      return SCIP_INVALIDDATA;
   }

   /* in case the coefficient is integral w.r.t. numerics we explicitly round the coefficient to an integral value */
   val = SCIPsetIsIntegral(set, val) ? SCIPsetRound(set, val) : val;
   col = row->cols[pos];
   assert(row->cols[pos] != NULL);

   if( SCIPsetIsZero(set, val) )
   {
      /* delete existing coefficient */
      SCIP_CALL( rowDelCoefPos(row, blkmem, set, eventqueue, lp, pos) );
   }
   else if( !SCIPsetIsEQ(set, row->vals[pos], val) )
   {
      SCIP_Real oldval;

      oldval = row->vals[pos];

      /* change existing coefficient */
      rowDelNorms(row, set, col, row->vals[pos], FALSE, FALSE, TRUE);
      row->vals[pos] = val;
      row->integral = row->integral && SCIPcolIsIntegral(col) && SCIPsetIsIntegral(set, val);
      rowAddNorms(row, set, col, row->vals[pos], TRUE);
      coefChanged(row, col, lp);

      /* issue row coefficient changed event */
      SCIP_CALL( rowEventCoefChanged(row, blkmem, set, eventqueue, col, oldval, val) );
   }

   return SCIP_OKAY;
}

/** notifies LP row, that its sides were changed */
static
SCIP_RETCODE rowSideChanged(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SIDETYPE         sidetype            /**< type of side: left or right hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);

   if( row->lpipos >= 0 )
   {
      /* insert row in the chgrows list (if not already there) */
      if( !row->lhschanged && !row->rhschanged )
      {
         SCIP_CALL( ensureChgrowsSize(lp, set, lp->nchgrows+1) );
         lp->chgrows[lp->nchgrows] = row;
         lp->nchgrows++;
      }

      /* mark side change in the row */
      switch( sidetype )
      {
      case SCIP_SIDETYPE_LEFT:
         row->lhschanged = TRUE;
         break;
      case SCIP_SIDETYPE_RIGHT:
         row->rhschanged = TRUE;
         break;
      default:
         SCIPerrorMessage("unknown row side type\n");
         SCIPABORT();
         return SCIP_INVALIDDATA;  /*lint !e527*/
      }

      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      assert(lp->nchgrows > 0);
   }

   return SCIP_OKAY;
}




/*
 * double linked coefficient matrix methods 
 */

/** insert column coefficients in corresponding rows */
static
SCIP_RETCODE colLink(
   SCIP_COL*             col,                /**< column data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( col->nunlinked > 0 )
   {
      SCIPsetDebugMsg(set, "linking column <%s>\n", SCIPvarGetName(col->var));

      /* unlinked rows can only be in the non-LP/unlinked rows part of the rows array */
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(!SCIPsetIsZero(set, col->vals[i]));
         if( col->linkpos[i] == -1 )
         {
            /* this call might swap the current row with the first non-LP/not linked row, but this is of no harm */
            SCIP_CALL( rowAddCoef(col->rows[i], blkmem, set, eventqueue, lp, col, col->vals[i], i) );
         }
         assert(col->rows[i]->cols[col->linkpos[i]] == col);
         assert(col->rows[i]->linkpos[col->linkpos[i]] == i);
         assert(col->nlprows == 0 || col->rows[col->nlprows-1]->cols[col->linkpos[col->nlprows-1]] == col);
         assert(col->nlprows == 0 || col->rows[col->nlprows-1]->linkpos[col->linkpos[col->nlprows-1]] == col->nlprows-1);
      }
   }
   assert(col->nunlinked == 0);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** removes column coefficients from corresponding rows */
static
SCIP_RETCODE colUnlink(
   SCIP_COL*             col,                /**< column data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( col->nunlinked < col->len )
   {
      SCIPsetDebugMsg(set, "unlinking column <%s>\n", SCIPvarGetName(col->var));
      for( i = 0; i < col->len; ++i )
      {
         if( col->linkpos[i] >= 0 )
         {
            assert(col->rows[i]->cols[col->linkpos[i]] == col);
            SCIP_CALL( rowDelCoefPos(col->rows[i], blkmem, set, eventqueue, lp, col->linkpos[i]) );
            col->linkpos[i] = -1;
            col->nunlinked++;
         }
      }
   }
   assert(col->nunlinked == col->len);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** insert row coefficients in corresponding columns */
static
SCIP_RETCODE rowLink(
   SCIP_ROW*             row,                /**< row data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->nunlinked > 0 )
   {
      SCIPsetDebugMsg(set, "linking row <%s>\n", row->name);

      /* unlinked columns can only be in the non-LP/unlinked columns part of the cols array */
      for( i = row->nlpcols; i < row->len; ++i )
      {
         assert(!SCIPsetIsZero(set, row->vals[i]));
         if( row->linkpos[i] == -1 )
         {
            /* this call might swap the current column with the first non-LP/not linked column, but this is of no harm */
            SCIP_CALL( colAddCoef(row->cols[i], blkmem, set, eventqueue, lp, row, row->vals[i], i) );
         }
         assert(row->cols[i]->rows[row->linkpos[i]] == row);
         assert(row->cols[i]->linkpos[row->linkpos[i]] == i);
         assert(row->nlpcols == 0 || row->cols[row->nlpcols-1]->rows[row->linkpos[row->nlpcols-1]] == row);
         assert(row->nlpcols == 0 || row->cols[row->nlpcols-1]->linkpos[row->linkpos[row->nlpcols-1]] == row->nlpcols-1);
      }
   }
   assert(row->nunlinked == 0);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** removes row coefficients from corresponding columns */
static
SCIP_RETCODE rowUnlink(
   SCIP_ROW*             row,                /**< row data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(row != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( row->nunlinked < row->len )
   {
      SCIPsetDebugMsg(set, "unlinking row <%s>\n", row->name);
      for( i = 0; i < row->len; ++i )
      {
         if( row->linkpos[i] >= 0 )
         {
            assert(row->cols[i]->rows[row->linkpos[i]] == row);
            SCIP_CALL( colDelCoefPos(row->cols[i], set, lp, row->linkpos[i]) );
            row->nunlinked++;
         }
      }
   }
   assert(row->nunlinked == row->len);

   return SCIP_OKAY;
}




/*
 * local LP parameter methods
 */

/** sets parameter of type int in LP solver, ignoring unknown parameters */
static
SCIP_RETCODE lpSetIntpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   int                   value,              /**< value to set parameter to */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   SCIP_RETCODE retcode;

   assert(lp != NULL);
   assert(success != NULL);

   retcode = SCIPlpiSetIntpar(lp->lpi, lpparam, value);

   /* check, if parameter is unknown */
   if( retcode == SCIP_PARAMETERUNKNOWN )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   *success = TRUE;

   return retcode;
}

/** sets parameter of type SCIP_Bool in LP solver, ignoring unknown parameters */
static
SCIP_RETCODE lpSetBoolpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Bool             value,              /**< value to set parameter to */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   return lpSetIntpar(lp, lpparam, (int)value, success);
}

/** sets parameter of type SCIP_Real in LP solver, ignoring unknown parameters */
static
SCIP_RETCODE lpSetRealpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Real             value,              /**< value to set parameter to */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   SCIP_RETCODE retcode;

   assert(lp != NULL);
   assert(success != NULL);

   retcode = SCIPlpiSetRealpar(lp->lpi, lpparam, value);

   /* check, if parameter is unknown */
   if( retcode == SCIP_PARAMETERUNKNOWN )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   *success = TRUE;

   return retcode;
}

#ifndef NDEBUG
/** checks, that parameter of type int in LP solver has the given value, ignoring unknown parameters */
static
SCIP_RETCODE lpCheckIntpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   int                   value               /**< value parameter should have */
   )
{
   SCIP_RETCODE retcode;
   int lpivalue;

   assert(lp != NULL);

   retcode = SCIPlpiGetIntpar(lp->lpi, lpparam, &lpivalue);

   /* ignore unknown parameter error */
   if( retcode == SCIP_PARAMETERUNKNOWN )
      return SCIP_OKAY;

   /* check value */
   assert(lpivalue == value);

   return retcode;
}

/** checks, that parameter of type SCIP_Bool in LP solver has the given value, ignoring unknown parameters */
static
SCIP_RETCODE lpCheckBoolpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Bool             value               /**< value parameter should have */
   )
{
   return lpCheckIntpar(lp, lpparam, (int)value);
}

/** checks, that parameter of type SCIP_Real in LP solver has the given value, ignoring unknown parameters */
static
SCIP_RETCODE lpCheckRealpar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPPARAM          lpparam,            /**< LP parameter */
   SCIP_Real             value               /**< value parameter should have */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Real lpivalue;

   assert(lp != NULL);

   retcode = SCIPlpiGetRealpar(lp->lpi, lpparam, &lpivalue);

   /* ignore unknown parameter error */
   if( retcode == SCIP_PARAMETERUNKNOWN )
      return SCIP_OKAY;

   /* check value */
   assert(lpivalue == value); /*lint !e777*/

   return retcode;
}
#else
#define lpCheckIntpar(lp, lpparam, value) SCIP_OKAY
#define lpCheckBoolpar(lp, lpparam, value) SCIP_OKAY
#define lpCheckRealpar(lp, lpparam, value) SCIP_OKAY
#endif

/** should the objective limit of the LP solver be disabled */
#define lpCutoffDisabled(set) (set->lp_disablecutoff == 1 || (set->nactivepricers > 0 && set->lp_disablecutoff == 2))

/** sets the objective limit of the LP solver
 *
 *  Note that we are always minimizing.
 */
static
SCIP_RETCODE lpSetObjlim(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             objlim              /**< new objective limit */
   )
{
   assert(lp != NULL);
   assert(set != NULL);

   /* We disabled the objective limit in the LP solver or we want so solve exactly and thus cannot rely on the LP
    * solver's objective limit handling, so we return here and do not apply the objective limit. */
   if( lpCutoffDisabled(set) || set->misc_exactsolve )
      return SCIP_OKAY;

   /* convert SCIP infinity value to lp-solver infinity value if necessary */
   if( SCIPsetIsInfinity(set, objlim) )
      objlim = SCIPlpiInfinity(lp->lpi);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_OBJLIM, lp->lpiobjlim) );

   if( objlim != lp->lpiobjlim ) /*lint !e777*/
   {
      SCIP_Bool success;

      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_OBJLIM, objlim, &success) );
      if( success )
      {
         /* mark the current solution invalid */
         lp->solved = FALSE;
         lp->primalfeasible = FALSE;
         lp->primalchecked = FALSE;
         lp->lpobjval = SCIP_INVALID;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         lp->lpiobjlim = objlim;
      }
   }

   return SCIP_OKAY;
}

/** sets the feasibility tolerance of the LP solver */
static
SCIP_RETCODE lpSetFeastol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             feastol,            /**< new feasibility tolerance */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(feastol >= 0.0);
   assert(success != NULL);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_FEASTOL, lp->lpifeastol) );

   if( feastol != lp->lpifeastol ) /*lint !e777*/
   {
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_FEASTOL, feastol, success) );
      if( *success )
      {
         if( lp->nrows > 0 && feastol < lp->lpifeastol )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->primalfeasible = FALSE;
            lp->primalchecked = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpifeastol = feastol;
      }
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the reduced costs feasibility tolerance of the LP solver */
static
SCIP_RETCODE lpSetDualfeastol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             dualfeastol,        /**< new reduced costs feasibility tolerance */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(dualfeastol >= 0.0);
   assert(success != NULL);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_DUALFEASTOL, lp->lpidualfeastol) );

   if( dualfeastol != lp->lpidualfeastol ) /*lint !e777*/
   {
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_DUALFEASTOL, dualfeastol, success) );
      if( *success )
      {
         if( lp->nrows > 0 && dualfeastol < lp->lpidualfeastol )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->dualfeasible = FALSE;
            lp->dualchecked = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpidualfeastol = dualfeastol;
      }
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the convergence tolerance used in barrier algorithm of the LP solver */
static
SCIP_RETCODE lpSetBarrierconvtol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             barrierconvtol,     /**< new convergence tolerance used in barrier algorithm */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(barrierconvtol >= 0.0);
   assert(success != NULL);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_BARRIERCONVTOL, lp->lpibarrierconvtol) );

   if( barrierconvtol != lp->lpibarrierconvtol ) /*lint !e777*/
   {
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_BARRIERCONVTOL, barrierconvtol, success) );
      if( *success )
      {
         if( lp->nrows > 0 && barrierconvtol < lp->lpibarrierconvtol
            && (lp->lastlpalgo == SCIP_LPALGO_BARRIER || lp->lastlpalgo == SCIP_LPALGO_BARRIERCROSSOVER) )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->dualfeasible = FALSE;
            lp->dualchecked = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpibarrierconvtol = barrierconvtol;
      }
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the FROMSCRATCH setting of the LP solver */
static
SCIP_RETCODE lpSetFromscratch(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             fromscratch,        /**< new FROMSCRATCH setting */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckBoolpar(lp, SCIP_LPPAR_FROMSCRATCH, lp->lpifromscratch) );

   if( fromscratch != lp->lpifromscratch )
   {
      SCIP_CALL( lpSetBoolpar(lp, SCIP_LPPAR_FROMSCRATCH, fromscratch, success) );
      if( *success )
         lp->lpifromscratch = fromscratch;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the FASTMIP setting of the LP solver */
static
SCIP_RETCODE lpSetFastmip(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   fastmip,            /**< new FASTMIP setting */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);
   assert(0 <= fastmip && fastmip <= 1);

   SCIP_CALL( lpCheckIntpar(lp, SCIP_LPPAR_FASTMIP, lp->lpifastmip) );

   if( fastmip != lp->lpifastmip )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_FASTMIP, fastmip, success) );
      if( *success )
         lp->lpifastmip = fastmip;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the SCALING setting of the LP solver */
static
SCIP_RETCODE lpSetScaling(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   scaling,            /**< new SCALING setting */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckIntpar(lp, SCIP_LPPAR_SCALING, lp->lpiscaling) );

   if( scaling != lp->lpiscaling )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_SCALING, scaling, success) );
      if( *success )
         lp->lpiscaling = scaling;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the number of THREADS  of the LP solver */
static
SCIP_RETCODE lpSetThreads(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   threads,            /**< new number of threads used to solve the LP */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckIntpar(lp, SCIP_LPPAR_THREADS, lp->lpithreads) );

   if( threads != lp->lpithreads )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_THREADS, threads, success) );
      if( *success )
         lp->lpithreads = threads;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the PRESOLVING setting of the LP solver */
static
SCIP_RETCODE lpSetPresolving(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             presolving,         /**< new PRESOLVING setting */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckBoolpar(lp, SCIP_LPPAR_PRESOLVING, lp->lpipresolving) );

   if( presolving != lp->lpipresolving )
   {
      SCIP_CALL( lpSetBoolpar(lp, SCIP_LPPAR_PRESOLVING, presolving, success) );
      if( *success )
         lp->lpipresolving = presolving;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the ROWREPSWITCH setting of the LP solver */
static
SCIP_RETCODE lpSetRowrepswitch(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             rowrepswitch,       /**< new ROWREPSWITCH value */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_ROWREPSWITCH, lp->lpirowrepswitch) );

   if( rowrepswitch != lp->lpirowrepswitch )  /*lint !e777*/
   {
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_ROWREPSWITCH, rowrepswitch, success) );
      if( *success )
         lp->lpirowrepswitch = rowrepswitch;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the iteration limit of the LP solver */
static
SCIP_RETCODE lpSetIterationLimit(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   itlim               /**< maximal number of LP iterations to perform, or -1 for no limit */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);
   assert(itlim >= -1);

   if( itlim == -1 )
      itlim = INT_MAX;

   SCIP_CALL( lpCheckIntpar(lp, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

   if( itlim != lp->lpiitlim )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_LPITLIM, itlim, &success) );
      if( success )
      {
         if( itlim > lp->lpiitlim )
         {
            /* mark the current solution invalid */
            lp->solved = FALSE;
            lp->lpobjval = SCIP_INVALID;
            lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
         }
         lp->lpiitlim = itlim;
      }
   }

   return SCIP_OKAY;
}

/** sets the pricing strategy of the LP solver */
static
SCIP_RETCODE lpSetPricing(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_PRICING          pricing             /**< pricing strategy */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);

   SCIP_CALL( lpCheckIntpar(lp, SCIP_LPPAR_PRICING, (int)lp->lpipricing) );

   if( pricing != lp->lpipricing )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_PRICING, (int)pricing, &success) );
      if( success )
         lp->lpipricing = pricing;
   }

   return SCIP_OKAY;
}

/** sets the pricing strategy of the LP solver (given the character representation of the strategy) */
static
SCIP_RETCODE lpSetPricingChar(
   SCIP_LP*              lp,                 /**< current LP data */
   char                  pricingchar         /**< character representing the pricing strategy */
   )
{
   SCIP_PRICING pricing;

   switch( pricingchar )
   {
   case 'l':
      pricing = SCIP_PRICING_LPIDEFAULT;
      break;
   case 'a':
      pricing = SCIP_PRICING_AUTO;
      break;
   case 'f':
      pricing = SCIP_PRICING_FULL;
      break;
   case 'p':
      pricing = SCIP_PRICING_PARTIAL;
      break;
   case 's':
      pricing = SCIP_PRICING_STEEP;
      break;
   case 'q':
      pricing = SCIP_PRICING_STEEPQSTART;
      break;
   case 'd':
      pricing = SCIP_PRICING_DEVEX;
      break;
   default:
      SCIPerrorMessage("invalid LP pricing parameter <%c>\n", pricingchar);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( lpSetPricing(lp, pricing) );

   return SCIP_OKAY;
}

/** sets the verbosity of the LP solver */
static
SCIP_RETCODE lpSetLPInfo(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             lpinfo              /**< should the LP solver display status messages? */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);

   SCIP_CALL( lpCheckBoolpar(lp, SCIP_LPPAR_LPINFO, lp->lpilpinfo) );

   if( lpinfo != lp->lpilpinfo )
   {
      SCIP_CALL( lpSetBoolpar(lp, SCIP_LPPAR_LPINFO, lpinfo, &success) );
      if( success )
         lp->lpilpinfo = lpinfo;
   }

   return SCIP_OKAY;
}

/** sets the CONDITIONLIMIT setting of the LP solver */
static
SCIP_RETCODE lpSetConditionLimit(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             condlimit,          /**< new CONDITIONLIMIT value */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   SCIP_CALL( lpCheckRealpar(lp, SCIP_LPPAR_CONDITIONLIMIT, lp->lpiconditionlimit) );

   if( condlimit != lp->lpiconditionlimit )  /*lint !e777*/
   {
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_CONDITIONLIMIT, condlimit, success) );
      if( *success )
         lp->lpiconditionlimit = condlimit;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the type of timer of the LP solver */
static
SCIP_RETCODE lpSetTiming(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLOCKTYPE        timing,             /**< new timing value */
   SCIP_Bool             enabled,            /**< is timing enabled? */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   int lptiming;

   assert(lp != NULL);
   assert(success != NULL);
   assert((int) SCIP_CLOCKTYPE_CPU == 1 && (int) SCIP_CLOCKTYPE_WALL == 2);

   SCIP_CALL( lpCheckIntpar(lp, SCIP_LPPAR_TIMING, lp->lpitiming) );

   if( !enabled )
      lptiming = 0;
   else
      lptiming = (int) timing;

   if( lptiming != lp->lpitiming )  /*lint !e777*/
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_TIMING, lptiming, success) );
      if( *success )
         lp->lpitiming = lptiming;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the initial random seed of the LP solver */
static
SCIP_RETCODE lpSetRandomseed(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   randomseed,         /**< new initial random seed */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   /* we don't check this parameter because SoPlex will always return its current random seed, not the initial one */

   if( randomseed == 0 )
   {
      lp->lpirandomseed = randomseed;
      *success = TRUE;
   }
   else if( randomseed != lp->lpirandomseed )  /*lint !e777*/
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_RANDOMSEED, randomseed, success) );
      if( *success )
         lp->lpirandomseed = randomseed;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the LP solution polishing method */
static
SCIP_RETCODE lpSetSolutionPolishing(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             polishing,          /**< LP solution polishing activated (0: disabled, 1: enabled) */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   if( polishing != lp->lpisolutionpolishing )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_POLISHING, (polishing ? 1 : 0), success) );
      if( *success )
         lp->lpisolutionpolishing = polishing;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** sets the LP refactorization interval */
static
SCIP_RETCODE lpSetRefactorInterval(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   refactor,           /**< LP refactorization interval (0: automatic) */
   SCIP_Bool*            success             /**< pointer to store whether the parameter was successfully changed */
   )
{
   assert(lp != NULL);
   assert(success != NULL);

   if( refactor != lp->lpirefactorinterval )
   {
      SCIP_CALL( lpSetIntpar(lp, SCIP_LPPAR_REFACTOR, refactor, success) );
      if( *success )
         lp->lpirefactorinterval = refactor;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}


/*
 * Column methods
 */

/** creates an LP column */
SCIP_RETCODE SCIPcolCreate(
   SCIP_COL**            col,                /**< pointer to column data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable, this column represents */
   int                   len,                /**< number of nonzeros in the column */
   SCIP_ROW**            rows,               /**< array with rows of column entries */
   SCIP_Real*            vals,               /**< array with coefficients of column entries */
   SCIP_Bool             removable           /**< should the column be removed from the LP due to aging or cleanup? */
   )
{
   int i;

   assert(col != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(var != NULL);
   assert(len >= 0);
   assert(len == 0 || (rows != NULL && vals != NULL));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, col) );

   if( len > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*col)->rows, rows, len) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*col)->vals, vals, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*col)->linkpos, len) );

      for( i = 0; i < len; ++i )
      {
         assert(rows[i] != NULL);
         assert(!SCIPsetIsZero(set, vals[i]));
         (*col)->linkpos[i] = -1;
      }
   }
   else
   {
      (*col)->rows = NULL;
      (*col)->vals = NULL;
      (*col)->linkpos = NULL;
   }

   (*col)->var = var;
   (*col)->obj = SCIPvarGetObj(var);
   (*col)->unchangedobj = SCIPvarGetUnchangedObj(var);
   (*col)->lb = SCIPvarGetLbLocal(var);
   (*col)->ub = SCIPvarGetUbLocal(var);
   (*col)->flushedobj = 0.0;
   (*col)->flushedlb = 0.0;
   (*col)->flushedub = 0.0;
   (*col)->index = stat->ncolidx;
   SCIPstatIncrement(stat, set, ncolidx);
   (*col)->size = len;
   (*col)->len = len;
   (*col)->nlprows = 0;
   (*col)->nunlinked = len;
   (*col)->lppos = -1;
   (*col)->lpipos = -1;
   (*col)->lpdepth = -1;
   (*col)->primsol = 0.0;
   (*col)->redcost = SCIP_INVALID;
   (*col)->farkascoef = SCIP_INVALID;
   (*col)->minprimsol = (*col)->ub;
   (*col)->maxprimsol = (*col)->lb;
   (*col)->sbdown = SCIP_INVALID;
   (*col)->sbup = SCIP_INVALID;
   (*col)->sbsolval  = SCIP_INVALID;
   (*col)->sblpobjval = SCIP_INVALID;
   (*col)->sbnode = -1;
   (*col)->validredcostlp = -1;
   (*col)->validfarkaslp = -1;
   (*col)->validsblp = -1;
   (*col)->sbitlim = -1;
   (*col)->nsbcalls = 0;
   (*col)->age = 0;
   (*col)->obsoletenode = -1;
   (*col)->var_probindex = SCIPvarGetProbindex(var);
   (*col)->basisstatus = SCIP_BASESTAT_ZERO; /*lint !e641*/
   (*col)->lprowssorted = TRUE;
   (*col)->nonlprowssorted = (len <= 1);
   (*col)->objchanged = FALSE;
   (*col)->lbchanged = FALSE;
   (*col)->ubchanged = FALSE;
   (*col)->coefchanged = FALSE;
   (*col)->integral = SCIPvarIsIntegral(var);
   (*col)->removable = removable;
   (*col)->sbdownvalid = FALSE;
   (*col)->sbupvalid = FALSE;
   (*col)->lazylb = SCIPvarGetLbLazy(var);
   (*col)->lazyub = SCIPvarGetUbLazy(var);
   (*col)->storedsolvals = NULL;

   return SCIP_OKAY;
}

/** frees an LP column */
SCIP_RETCODE SCIPcolFree(
   SCIP_COL**            col,                /**< pointer to LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(col != NULL);
   assert(*col != NULL);
   assert((*col)->var != NULL);
   assert(SCIPvarGetStatus((*col)->var) == SCIP_VARSTATUS_COLUMN);
   assert(&(*col)->var->data.col == col); /* SCIPcolFree() has to be called from SCIPvarFree() */
   assert((*col)->lppos == -1);
   assert((*col)->lpipos == -1);

   /* remove column indices from corresponding rows */
   SCIP_CALL( colUnlink(*col, blkmem, set, eventqueue, lp) );

   BMSfreeBlockMemoryNull(blkmem, &(*col)->storedsolvals);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*col)->rows, (*col)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*col)->vals, (*col)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*col)->linkpos, (*col)->size);
   BMSfreeBlockMemory(blkmem, col);

   return SCIP_OKAY;
}

/** output column to file stream */
void SCIPcolPrint(
   SCIP_COL*             col,                /**< LP column */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int r;

   assert(col != NULL);
   assert(col->var != NULL);

   /* print bounds */
   SCIPmessageFPrintInfo(messagehdlr, file, "(obj: %.15g) [%.15g,%.15g], ", col->obj, col->lb, col->ub);

   /* print coefficients */
   if( col->len == 0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "<empty>");
   for( r = 0; r < col->len; ++r )
   {
      assert(col->rows[r] != NULL);
      assert(col->rows[r]->name != NULL);
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.15g<%s> ", col->vals[r], col->rows[r]->name);
   }
   SCIPmessageFPrintInfo(messagehdlr, file, "\n");
}

/** sorts column entries such that LP rows precede non-LP rows and inside both parts lower row indices precede higher ones
 */
void SCIPcolSort(
   SCIP_COL*             col                 /**< column to be sorted */
   )
{
   /* sort LP rows */
   colSortLP(col);

   /* sort non-LP rows */
   colSortNonLP(col);
}

/** adds a previously non existing coefficient to an LP column */
SCIP_RETCODE SCIPcolAddCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   SCIP_CALL( colAddCoef(col, blkmem, set, eventqueue, lp, row, val, -1) );

   checkLinks(lp);

   return SCIP_OKAY;
}

/** deletes existing coefficient from column */
SCIP_RETCODE SCIPcolDelCoef(
   SCIP_COL*             col,                /**< column to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(col != NULL);
   assert(col->var != NULL);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);

   /* search the position of the row in the column's row vector */
   pos = colSearchCoef(col, row);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for row <%s> doesn't exist in column <%s>\n", row->name, SCIPvarGetName(col->var));
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < col->len);
   assert(col->rows[pos] == row);

   /* if row knows of the column, remove the column from the row's col vector */
   if( col->linkpos[pos] >= 0 )
   {
      assert(row->cols[col->linkpos[pos]] == col);
      assert(row->cols_index[col->linkpos[pos]] == col->index);
      assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
      SCIP_CALL( rowDelCoefPos(row, blkmem, set, eventqueue, lp, col->linkpos[pos]) );
   }

   /* delete the row from the column's row vector */
   SCIP_CALL( colDelCoefPos(col, set, lp, pos) );

   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP column */
SCIP_RETCODE SCIPcolChgCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   int pos;

   assert(col != NULL);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);

   /* search the position of the row in the column's row vector */
   pos = colSearchCoef(col, row);

   /* check, if row already exists in the column's row vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( colAddCoef(col, blkmem, set, eventqueue, lp, row, val, -1) );
   }
   else
   {
      /* modify already existing coefficient */
      assert(0 <= pos && pos < col->len);
      assert(col->rows[pos] == row);

      /* if row knows of the column, change the corresponding coefficient in the row */
      if( col->linkpos[pos] >= 0 )
      {
         assert(row->cols[col->linkpos[pos]] == col);
         assert(row->cols_index[col->linkpos[pos]] == col->index);
         assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
         SCIP_CALL( rowChgCoefPos(row, blkmem, set, eventqueue, lp, col->linkpos[pos], val) );
      }

      /* change the coefficient in the column */
      SCIP_CALL( colChgCoefPos(col, set, lp, pos, val) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** increases value of an existing or non-existing coefficient in an LP column */
SCIP_RETCODE SCIPcolIncCoef(
   SCIP_COL*             col,                /**< LP column */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             incval              /**< value to add to the coefficient */
   )
{
   int pos;

   assert(col != NULL);
   assert(lp != NULL);
   assert(!lp->diving);
   assert(row != NULL);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   /* search the position of the row in the column's row vector */
   pos = colSearchCoef(col, row);

   /* check, if row already exists in the column's row vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( colAddCoef(col, blkmem, set, eventqueue, lp, row, incval, -1) );
   }
   else
   {
      /* modify already existing coefficient */
      assert(0 <= pos && pos < col->len);
      assert(col->rows[pos] == row);

      /* if row knows of the column, change the corresponding coefficient in the row */
      if( col->linkpos[pos] >= 0 )
      {
         assert(row->cols[col->linkpos[pos]] == col);
         assert(row->cols_index[col->linkpos[pos]] == col->index);
         assert(SCIPsetIsEQ(set, row->vals[col->linkpos[pos]], col->vals[pos]));
         SCIP_CALL( rowChgCoefPos(row, blkmem, set, eventqueue, lp, col->linkpos[pos], col->vals[pos] + incval) );
      }

      /* change the coefficient in the column */
      SCIP_CALL( colChgCoefPos(col, set, lp, pos, col->vals[pos] + incval) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** insert column in the chgcols list (if not already there) */
static
SCIP_RETCODE insertColChgcols(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   if( !col->objchanged && !col->lbchanged && !col->ubchanged )
   {
      SCIP_CALL( ensureChgcolsSize(lp, set, lp->nchgcols+1) );
      lp->chgcols[lp->nchgcols] = col;
      lp->nchgcols++;
   }

   /* mark the current LP unflushed */
   lp->flushed = FALSE;

   return SCIP_OKAY;
}

/** Is the new value reliable or may we have cancellation?
 *
 *  @note: Here we only consider cancellations which can occur during decreasing the oldvalue to newvalue; not the
 *  cancellations which can occur during increasing the oldvalue to the newvalue
 */
static
SCIP_Bool isNewValueUnreliable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newvalue,           /**< new value */
   SCIP_Real             oldvalue            /**< old reliable value */
   )
{
   SCIP_Real quotient;

   assert(set != NULL);
   assert(oldvalue < SCIP_INVALID);

   quotient = (REALABS(newvalue)+1.0) / (REALABS(oldvalue) + 1.0);

   return SCIPsetIsZero(set, quotient);
}

/** update norms of objective function vector */
static
void lpUpdateObjNorms(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             oldobj,             /**< old objective value of variable */
   SCIP_Real             newobj              /**< new objective value of variable */
   )
{
   if( REALABS(newobj) != REALABS(oldobj) )   /*lint !e777*/
   {
      if( !lp->objsqrnormunreliable )
      {
         SCIP_Real oldvalue;

         oldvalue = lp->objsqrnorm;
         lp->objsqrnorm += SQR(newobj) - SQR(oldobj);

         /* due to numerical cancellations, we recalculate lp->objsqrnorm using all variables */
         if( SCIPsetIsLT(set, lp->objsqrnorm, 0.0) || isNewValueUnreliable(set, lp->objsqrnorm, oldvalue) )
            lp->objsqrnormunreliable = TRUE;
         else
         {
            assert(SCIPsetIsGE(set, lp->objsqrnorm, 0.0));

            /* due to numerical troubles it still can appear that lp->objsqrnorm is a little bit smaller than 0 */
            lp->objsqrnorm = MAX(lp->objsqrnorm, 0.0);

            assert(lp->objsqrnorm >= 0.0);
         }
      }

      lp->objsumnorm += REALABS(newobj) - REALABS(oldobj);
      lp->objsumnorm = MAX(lp->objsumnorm, 0.0);
   }
}

/** changes objective value of column */
SCIP_RETCODE SCIPcolChgObj(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "changing objective value of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->obj, newobj);

   /* only add actual changes */
   if( !SCIPsetIsEQ(set, col->obj, newobj) )
   {
      /* only variables with a real position in the LPI can be inserted */
      if( col->lpipos >= 0 )
      {
         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lp) );

         /* mark objective value change in the column */
         col->objchanged = TRUE;

         assert(lp->nchgcols > 0);
      }
      /* in any case, when the sign of the objective (and thereby the best bound) changes, the variable has to enter the
       * LP and the LP has to be flushed
       */
      else if( (col->obj < 0.0 && newobj >= 0.0 && SCIPsetIsZero(set, col->ub))
         || (col->obj >= 0.0 && newobj < 0.0 && SCIPsetIsZero(set, col->lb)) )
      {
         /* mark the LP unflushed */
         lp->flushed = FALSE;
      }
   }

   /* store new objective function value */
   col->obj = newobj;

   /* update original objective value, as long as we are not in diving or probing and changed objective values */
   if( !lp->divingobjchg )
   {
      SCIP_Real oldobj = col->unchangedobj;

      assert(SCIPsetIsEQ(set, newobj, SCIPvarGetUnchangedObj(col->var)));
      col->unchangedobj = newobj;

      /* update the objective function vector norms */
      lpUpdateObjNorms(lp, set, oldobj, newobj);
   }

   return SCIP_OKAY;
}

/** changes lower bound of column */
SCIP_RETCODE SCIPcolChgLb(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newlb               /**< new lower bound value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "changing lower bound of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->lb, newlb);

   /* only add actual changes */
   if( !SCIPsetIsEQ(set, col->lb, newlb) )
   {
      /* only variables with a real position in the LPI can be inserted */
      if( col->lpipos >= 0 )
      {
         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lp) );

         /* mark bound change in the column */
         col->lbchanged = TRUE;

         assert(lp->nchgcols > 0);
      }
      /* in any case, when the best bound is zero and gets changed, the variable has to enter the LP and the LP has to be
       * flushed
       */
      else if( col->obj >= 0.0 && SCIPsetIsZero(set, col->lb) )
      {
         /* mark the LP unflushed */
         lp->flushed = FALSE;
      }
   }

   col->lb = newlb;

   return SCIP_OKAY;
}

/** changes upper bound of column */
SCIP_RETCODE SCIPcolChgUb(
   SCIP_COL*             col,                /**< LP column to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newub               /**< new upper bound value */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "changing upper bound of column <%s> from %f to %f\n", SCIPvarGetName(col->var), col->ub, newub);

   /* only add actual changes */
   if( !SCIPsetIsEQ(set, col->ub, newub) )
   {
      /* only variables with a real position in the LPI can be inserted */
      if( col->lpipos >= 0 )
      {
         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lp) );

         /* mark bound change in the column */
         col->ubchanged = TRUE;

         assert(lp->nchgcols > 0);
      }
      /* in any case, when the best bound is zero and gets changed, the variable has to enter the LP and the LP has to be
       * flushed
       */
      else if( col->obj < 0.0 && SCIPsetIsZero(set, col->ub) )
      {
         /* mark the LP unflushed */
         lp->flushed = FALSE;
      }
   }

   col->ub = newub;

   return SCIP_OKAY;
}

/** calculates the reduced costs of a column using the given dual solution vector */
SCIP_Real SCIPcolCalcRedcost(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            dualsol             /**< dual solution vector for current LP rows */
   )
{
   SCIP_ROW* row;
   SCIP_Real redcost;
   int i;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(dualsol != NULL);

   redcost = col->obj;
   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->lppos >= 0);
      redcost -= col->vals[i] * dualsol[row->lppos];
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
            redcost -= col->vals[i] * dualsol[row->lppos];
      }
   }
#ifndef NDEBUG
   else
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos == -1);
         assert(col->linkpos[i] >= 0);
      }
   }
#endif

   return redcost;
}

/** calculates the reduced costs of a column using the dual solution stored in the rows */
static
SCIP_Real colCalcInternalRedcost(
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_ROW* row;
   SCIP_Real redcost;
   int i;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);

   redcost = col->obj;
   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->dualsol < SCIP_INVALID);
      assert(row->lppos >= 0);
      assert(col->linkpos[i] >= 0);
      redcost -= col->vals[i] * row->dualsol;
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos >= 0 || row->dualsol == 0.0);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
            redcost -= col->vals[i] * row->dualsol;
      }
   }
#ifndef NDEBUG
   else
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->dualsol == 0.0);
         assert(row->lppos == -1);
         assert(col->linkpos[i] >= 0);
      }
   }
#endif

   return redcost;
}

/** gets the reduced costs of a column in last LP or after recalculation */
SCIP_Real SCIPcolGetRedcost(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(col != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(col->validredcostlp <= stat->lpcount);
   assert(lp->validsollp == stat->lpcount);

   if( col->validredcostlp < stat->lpcount )
   {
      col->redcost = colCalcInternalRedcost(col);
      col->validredcostlp = stat->lpcount;
   }
   assert(col->validredcostlp == stat->lpcount);
   assert(col->redcost < SCIP_INVALID);

   return col->redcost;
}

/** gets the feasibility of (the dual row of) a column in last LP or after recalculation */
SCIP_Real SCIPcolGetFeasibility(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(col != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->validsollp == stat->lpcount);

   /* A column's reduced cost is defined as
    *   redcost  = obj - activity,  activity = y^T * col.   (activity = obj - redcost)
    * The activity is equal to the activity of the corresponding row in the dual LP.
    * The column's feasibility is the feasibility of the corresponding row in the dual LP.
    * The sides of the dual row depend on the bounds of the column:
    *  - lb == ub      :  dual row is a free row with infinite sides
    *  -  0 <= lb <  ub:         activity <= obj  =>  0 <= redcost
    *  - lb <   0 <  ub:  obj <= activity <= obj  =>  0 <= redcost <= 0
    *  - lb <  ub <=  0:  obj <= activity         =>       redcost <= 0
    */
   if( SCIPsetIsEQ(set, col->lb, col->ub) )
   {
      /* dual row is free */
      return SCIPsetInfinity(set);
   }
   else
   {
      SCIP_Real redcost;

      /* calculate reduced costs */
      redcost = SCIPcolGetRedcost(col, stat, lp);

      if( !SCIPsetIsNegative(set, col->lb) )
      {
         /* dual row is  activity <= obj  <=>  redcost >= 0 */
         return redcost;
      }
      else if( SCIPsetIsPositive(set, col->ub) )
      {
         /* dual row is  activity == obj  <=>  redcost == 0 */
         return -REALABS(redcost);
      }
      else
      {
         /* dual row is  activity >= obj  <=>  redcost <= 0 */
         return -redcost;
      }
   }
}

/** calculates the Farkas coefficient y^T A_i of a column i using the given dual Farkas vector y */
SCIP_Real SCIPcolCalcFarkasCoef(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            dualfarkas          /**< dense dual Farkas vector for current LP rows */
   )
{
   SCIP_ROW* row;
   SCIP_Real farkas;
   int i;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(dualfarkas != NULL);

   farkas = 0.0;
   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->lppos >= 0);
      farkas += col->vals[i] * dualfarkas[row->lppos];
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
            farkas += col->vals[i] * dualfarkas[row->lppos];
      }
   }
#ifndef NDEBUG
   else
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos == -1);
         assert(col->linkpos[i] >= 0);
      }
   }
#endif

   return farkas;
}

/** gets the Farkas coefficient y^T A_i of a column i in last LP (which must be infeasible) */
static
SCIP_Real colCalcInternalFarkasCoef(
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_ROW* row;
   SCIP_Real farkas;
   int i;

   assert(col != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);

   farkas = 0.0;
   for( i = 0; i < col->nlprows; ++i )
   {
      row = col->rows[i];
      assert(row != NULL);
      assert(row->dualfarkas < SCIP_INVALID);
      assert(row->lppos >= 0);
      assert(col->linkpos[i] >= 0);
      farkas += col->vals[i] * row->dualfarkas;
   }

   if( col->nunlinked > 0 )
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->lppos >= 0 || row->dualfarkas == 0.0);
         assert(row->lppos == -1 || col->linkpos[i] == -1);
         if( row->lppos >= 0 )
            farkas += col->vals[i] * row->dualfarkas;
      }
   }
#ifndef NDEBUG
   else
   {
      for( i = col->nlprows; i < col->len; ++i )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->dualfarkas == 0.0);
         assert(row->lppos == -1);
         assert(col->linkpos[i] >= 0);
      }
   }
#endif

   return farkas;
}

/** gets the Farkas coefficient of a column in last LP (which must be infeasible) */
SCIP_Real SCIPcolGetFarkasCoef(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(col != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(col->validfarkaslp <= stat->lpcount);
   assert(lp->validfarkaslp == stat->lpcount);

   if( col->validfarkaslp < stat->lpcount )
   {
      col->farkascoef = colCalcInternalFarkasCoef(col);
      col->validfarkaslp = stat->lpcount;
   }
   assert(col->validfarkaslp == stat->lpcount);
   assert(col->farkascoef < SCIP_INVALID);

   return col->farkascoef;
}

/** gets the Farkas value of a column in last LP (which must be infeasible), i.e. the Farkas coefficient y^T A_i times
 *  the best bound for this coefficient, i.e. max{y^T A_i x_i | lb <= x_i <= ub}
 */
SCIP_Real SCIPcolGetFarkasValue(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real farkascoef;

   assert(col != NULL);

   farkascoef = SCIPcolGetFarkasCoef(col, stat, lp);

   if( farkascoef > 0.0 )
      return col->ub * farkascoef;
   else
      return col->lb * farkascoef;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPlpStartStrongbranch(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);
   assert(!lp->strongbranching);

   lp->strongbranching = TRUE;
   SCIPdebugMessage("starting strong branching ...\n");
   SCIP_CALL( SCIPlpiStartStrongbranch(lp->lpi) );

   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPlpEndStrongbranch(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);
   assert(lp->strongbranching);

   lp->strongbranching = FALSE;
   SCIPdebugMessage("ending strong branching ...\n");
   SCIP_CALL( SCIPlpiEndStrongbranch(lp->lpi) );

   return SCIP_OKAY;
}

/** sets strong branching information for a column variable */
void SCIPcolSetStrongbranchData(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Real             lpobjval,           /**< objective value of the current LP */
   SCIP_Real             primsol,            /**< primal solution value of the column in the current LP */
   SCIP_Real             sbdown,             /**< dual bound after branching column down */
   SCIP_Real             sbup,               /**< dual bound after branching column up */
   SCIP_Bool             sbdownvalid,        /**< is the returned down value a valid dual bound? */
   SCIP_Bool             sbupvalid,          /**< is the returned up value a valid dual bound? */
   SCIP_Longint          iter,               /**< total number of strong branching iterations */
   int                   itlim               /**< iteration limit applied to the strong branching call */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPcolIsIntegral(col));
   assert(SCIPvarIsIntegral(col->var));
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(col->lpipos >= 0);
   assert(col->lppos >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->strongbranchprobing);
   assert(col->lppos < lp->ncols);
   assert(lp->cols[col->lppos] == col);
   assert(itlim >= 1);

   col->sblpobjval = lpobjval;
   col->sbsolval = primsol;
   col->validsblp = stat->lpcount - stat->nsbdivinglps;
   col->sbnode = stat->nnodes;

   col->sbitlim = itlim;
   col->nsbcalls++;

   col->sbdown = MIN(sbdown, lp->cutoffbound);
   col->sbup = MIN(sbup, lp->cutoffbound);
   col->sbdownvalid = sbdownvalid;
   col->sbupvalid = sbupvalid;

   SCIPstatIncrement(stat, set, nstrongbranchs);
   SCIPstatAdd(stat, set, nsblpiterations, iter);
   if( stat->nnodes == 1 )
   {
      SCIPstatIncrement(stat, set, nrootstrongbranchs);
      SCIPstatAdd(stat, set, nrootsblpiterations, iter);
   }
}

/** invalidates strong branching information for a column variable */
void SCIPcolInvalidateStrongbranchData(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPcolIsIntegral(col));
   assert(SCIPvarIsIntegral(col->var));
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(col->lpipos >= 0);
   assert(col->lppos >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->strongbranchprobing);
   assert(col->lppos < lp->ncols);
   assert(lp->cols[col->lppos] == col);

   col->sbdown = SCIP_INVALID;
   col->sbup = SCIP_INVALID;
   col->sbdownvalid = FALSE;
   col->sbupvalid = FALSE;
   col->validsblp = -1;
   col->sbsolval = SCIP_INVALID;
   col->sblpobjval = SCIP_INVALID;
   col->sbnode = -1;
   col->sbitlim = -1;
}


/** gets strong branching information on a column variable */
SCIP_RETCODE SCIPcolGetStrongbranch(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Bool             integral,           /**< should integral strong branching be performed? */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   assert(col != NULL);
   assert(col->var != NULL);
   assert(SCIPcolIsIntegral(col));
   assert(SCIPvarIsIntegral(col->var));
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(col->primsol < SCIP_INVALID);
   assert(col->lpipos >= 0);
   assert(col->lppos >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->strongbranching);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);
   assert(lp->validsollp == stat->lpcount);
   assert(col->lppos < lp->ncols);
   assert(lp->cols[col->lppos] == col);
   assert(itlim >= 1);
   /*   assert(down != NULL);
    *  assert(up != NULL); temporary hack for cloud branching
    */
   assert(lperror != NULL);

   *lperror = FALSE;

   if( col->validsblp != stat->lpcount - stat->nsbdivinglps || itlim > col->sbitlim )
   {
      col->validsblp = stat->lpcount - stat->nsbdivinglps;
      col->sbsolval = col->primsol;
      col->sblpobjval = SCIPlpGetObjval(lp, set, prob);
      col->sbnode = stat->nnodes;
      assert(integral || !SCIPsetIsFeasIntegral(set, col->primsol));

      /* if a loose variables has an infinite best bound, the LP bound is -infinity and no gain can be achieved */
      if( lp->looseobjvalinf > 0 )
      {
         col->sbdown = -SCIPsetInfinity(set);
         col->sbup = -SCIPsetInfinity(set);
         col->sbdownvalid = FALSE;
         col->sbupvalid = FALSE;
      }
      else
      {
         SCIP_RETCODE retcode;
         SCIP_Real sbdown;
         SCIP_Real sbup;
         SCIP_Bool sbdownvalid;
         SCIP_Bool sbupvalid;
         int iter;

         SCIPsetDebugMsg(set, "performing strong branching on variable <%s>(%g) with %d iterations\n",
            SCIPvarGetName(col->var), col->primsol, itlim);

         /* start timing */
         SCIPclockStart(stat->strongbranchtime, set);

         /* call LPI strong branching */
         col->sbitlim = itlim;
         col->nsbcalls++;

         sbdown = lp->lpobjval;
         sbup = lp->lpobjval;

         if( integral )
            retcode = SCIPlpiStrongbranchInt(lp->lpi, col->lpipos, col->primsol, itlim, down  == NULL ? NULL : &sbdown, up  == NULL ? NULL : &sbup, &sbdownvalid, &sbupvalid, &iter);
         else
         {
            assert( ! SCIPsetIsIntegral(set, col->primsol) );
            retcode = SCIPlpiStrongbranchFrac(lp->lpi, col->lpipos, col->primsol, itlim, down == NULL ? NULL : &sbdown, up == NULL ? NULL :  &sbup, &sbdownvalid, &sbupvalid, &iter);
         }


         /* check return code for errors */
         if( retcode == SCIP_LPERROR )
         {
            *lperror = TRUE;
            col->sbdown = SCIP_INVALID;
            col->sbup = SCIP_INVALID;
            col->sbdownvalid = FALSE;
            col->sbupvalid = FALSE;
            col->validsblp = -1;
            col->sbsolval = SCIP_INVALID;
            col->sblpobjval = SCIP_INVALID;
            col->sbnode = -1;
         }
         else
         {
            SCIP_Real looseobjval;

            *lperror = FALSE;
            SCIP_CALL( retcode );

            looseobjval = getFiniteLooseObjval(lp, set, prob);
            col->sbdown = MIN(sbdown + looseobjval, lp->cutoffbound);
            col->sbup = MIN(sbup + looseobjval, lp->cutoffbound);

            col->sbdownvalid = sbdownvalid;
            col->sbupvalid = sbupvalid;

            /* update strong branching statistics */
            if( iter == -1 )
            {
               /* calculate average iteration number */
               iter = stat->ndualresolvelps > 0 ? (int)(2*stat->ndualresolvelpiterations / stat->ndualresolvelps)
                  : stat->nduallps > 0 ? (int)((stat->nduallpiterations / stat->nduallps) / 5)
                  : stat->nprimalresolvelps > 0 ? (int)(2*stat->nprimalresolvelpiterations / stat->nprimalresolvelps)
                  : stat->nprimallps > 0 ? (int)((stat->nprimallpiterations / stat->nprimallps) / 5)
                  : 0;
               if( iter/2 >= itlim )
                  iter = 2*itlim;
            }
            SCIPstatIncrement(stat, set, nstrongbranchs);
            SCIPstatAdd(stat, set, nsblpiterations, iter);
            if( stat->nnodes == 1 )
            {
               SCIPstatIncrement(stat, set, nrootstrongbranchs);
               SCIPstatAdd(stat, set, nrootsblpiterations, iter);
            }
         }

         /* stop timing */
         SCIPclockStop(stat->strongbranchtime, set);
      }
   }
   assert(*lperror || col->sbdown < SCIP_INVALID);
   assert(*lperror || col->sbup < SCIP_INVALID);

   if( down != NULL)
      *down = col->sbdown;
   if( up != NULL )
      *up = col->sbup;
   if( downvalid != NULL )
      *downvalid = col->sbdownvalid;
   if( upvalid != NULL )
      *upvalid = col->sbupvalid;

   return SCIP_OKAY;
}

/** gets strong branching information on column variables */
SCIP_RETCODE SCIPcolGetStrongbranches(
   SCIP_COL**            cols,               /**< LP columns */
   int                   ncols,              /**< number of columns */
   SCIP_Bool             integral,           /**< should integral strong branching be performed? */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Real* sbdown;
   SCIP_Real* sbup;
   SCIP_Bool* sbdownvalid;
   SCIP_Bool* sbupvalid;
   SCIP_Real* primsols;
   SCIP_COL** subcols;
   int* lpipos;
   int* subidx;
   int nsubcols;
   int iter;
   int j;

   assert(cols != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL);
   assert(lp->validsollp == stat->lpcount);
   assert(itlim >= 1);
   assert(down != NULL);
   assert(up != NULL);
   assert(lperror != NULL);

   *lperror = FALSE;

   if ( ncols <= 0 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(stat->strongbranchtime, set);

   /* initialize storage */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &subcols, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &subidx, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lpipos, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsols, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &sbdown, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &sbup, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &sbdownvalid, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &sbupvalid, ncols) );

   nsubcols = 0;
   for( j = 0; j < ncols; ++j )
   {
      SCIP_COL* col;
      col = cols[j];

      assert(col->lppos < lp->ncols);
      assert(lp->cols[col->lppos] == col);
      assert(SCIPcolIsIntegral(col));
      assert(SCIPvarIsIntegral(col->var));
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);
      assert(col->primsol < SCIP_INVALID);
      assert(col->lpipos >= 0);
      assert(col->lppos >= 0);

      if( col->validsblp != stat->lpcount - stat->nsbdivinglps || itlim > col->sbitlim )
      {
         col->validsblp = stat->lpcount - stat->nsbdivinglps;
         col->sbsolval = col->primsol;
         col->sblpobjval = SCIPlpGetObjval(lp, set, prob);
         col->sbnode = stat->nnodes;
         assert(!SCIPsetIsFeasIntegral(set, col->primsol));

         /* if a loose variables has an infinite best bound, the LP bound is -infinity and no gain can be achieved */
         if( lp->looseobjvalinf > 0 )
         {
            /* directly set up column and result vectors*/
            col->sbdown = -SCIPsetInfinity(set);
            col->sbup = -SCIPsetInfinity(set);
            col->sbdownvalid = FALSE;
            col->sbupvalid = FALSE;
            down[j] = col->sbdown;
            up[j] = col->sbup;
            if( downvalid != NULL )
               downvalid[j] = col->sbdownvalid;
            if( upvalid != NULL )
               upvalid[j] = col->sbupvalid;
         }
         else
         {
            col->sbitlim = itlim;
            col->nsbcalls++;

            lpipos[nsubcols] = col->lpipos;
            primsols[nsubcols] = col->primsol;
            assert( integral || ! SCIPsetIsFeasIntegral(set, col->primsol) );
            subidx[nsubcols] = j;
            subcols[nsubcols++] = col;
         }
      }
      else
      {
         /* directly set up resulting values (use stored values) */
         down[j] = col->sbdown;
         up[j] = col->sbup;
         if( downvalid != NULL )
            downvalid[j] = col->sbdownvalid;
         if( upvalid != NULL )
            upvalid[j] = col->sbupvalid;
      }
   }

   SCIPsetDebugMsg(set, "performing strong branching on %d variables with %d iterations\n", ncols, itlim);

   /* call LPI strong branching */
   if ( integral )
      retcode = SCIPlpiStrongbranchesInt(lp->lpi, lpipos, nsubcols, primsols, itlim, sbdown, sbup, sbdownvalid, sbupvalid, &iter);
   else
      retcode = SCIPlpiStrongbranchesFrac(lp->lpi, lpipos, nsubcols, primsols, itlim, sbdown, sbup, sbdownvalid, sbupvalid, &iter);

   /* check return code for errors */
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;

      for( j = 0; j < nsubcols; ++j )
      {
         SCIP_COL* col;
         int idx;

         col = subcols[j];
         idx = subidx[j];

         col->sbdown = SCIP_INVALID;
         col->sbup = SCIP_INVALID;
         col->sbdownvalid = FALSE;
         col->sbupvalid = FALSE;
         col->validsblp = -1;
         col->sbsolval = SCIP_INVALID;
         col->sblpobjval = SCIP_INVALID;
         col->sbnode = -1;

         down[idx] = col->sbdown;
         up[idx] = col->sbup;
         if( downvalid != NULL )
            downvalid[idx] = col->sbdownvalid;
         if( upvalid != NULL )
            upvalid[idx] = col->sbupvalid;
      }
   }
   else
   {
      SCIP_Real looseobjval;

      *lperror = FALSE;
      SCIP_CALL( retcode );

      looseobjval = getFiniteLooseObjval(lp, set, prob);

      for( j = 0; j < nsubcols; ++j )
      {
         SCIP_COL* col;
         int idx;

         col = subcols[j];
         idx = subidx[j];

         assert( col->sbdown < SCIP_INVALID);
         assert( col->sbup < SCIP_INVALID);

         col->sbdown = MIN(sbdown[j] + looseobjval, lp->cutoffbound);
         col->sbup = MIN(sbup[j] + looseobjval, lp->cutoffbound);
         col->sbdownvalid = sbdownvalid[j];
         col->sbupvalid = sbupvalid[j];

         down[idx] = col->sbdown;
         up[idx] = col->sbup;
         if( downvalid != NULL )
            downvalid[idx] = col->sbdownvalid;
         if( upvalid != NULL )
            upvalid[idx] = col->sbupvalid;
      }

      /* update strong branching statistics */
      if( iter == -1 )
      {
         /* calculate average iteration number */
         iter = stat->ndualresolvelps > 0 ? (int)(2*stat->ndualresolvelpiterations / stat->ndualresolvelps)
            : stat->nduallps > 0 ? (int)((stat->nduallpiterations / stat->nduallps) / 5)
            : stat->nprimalresolvelps > 0 ? (int)(2*stat->nprimalresolvelpiterations / stat->nprimalresolvelps)
            : stat->nprimallps > 0 ? (int)((stat->nprimallpiterations / stat->nprimallps) / 5)
            : 0;
         if( iter/2 >= itlim )
            iter = 2*itlim;
      }
      SCIPstatAdd(stat, set, nstrongbranchs, ncols);
      SCIPstatAdd(stat, set, nsblpiterations, iter);
      if( stat->nnodes == 1 )
      {
         SCIPstatAdd(stat, set, nrootstrongbranchs, ncols);
         SCIPstatAdd(stat, set, nrootsblpiterations, iter);
      }
   }

   SCIPsetFreeBufferArray(set, &sbupvalid);
   SCIPsetFreeBufferArray(set, &sbdownvalid);
   SCIPsetFreeBufferArray(set, &sbup);
   SCIPsetFreeBufferArray(set, &sbdown);
   SCIPsetFreeBufferArray(set, &primsols);
   SCIPsetFreeBufferArray(set, &lpipos);
   SCIPsetFreeBufferArray(set, &subidx);
   SCIPsetFreeBufferArray(set, &subcols);

   /* stop timing */
   SCIPclockStop(stat->strongbranchtime, set);

   return SCIP_OKAY;
}

/** gets last strong branching information available for a column variable;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given column;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
void SCIPcolGetStrongbranchLast(
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real*            down,               /**< stores dual bound after branching column down, or NULL */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up, or NULL */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Real*            solval,             /**< stores LP solution value of column at last strong branching call, or NULL */
   SCIP_Real*            lpobjval            /**< stores LP objective value at last strong branching call, or NULL */
   )
{
   assert(col != NULL);

   if( down != NULL )
      *down = col->sbdown;
   if( up != NULL )
      *up = col->sbup;
   if( downvalid != NULL )
      *downvalid = col->sbdownvalid;
   if( upvalid != NULL )
      *upvalid = col->sbupvalid;
   if( solval != NULL )
      *solval = col->sbsolval;
   if( lpobjval != NULL )
      *lpobjval = col->sblpobjval;
}

/** if strong branching was already applied on the column at the current node, returns the number of LPs solved after
 *  the LP where the strong branching on this column was applied;
 *  if strong branching was not yet applied on the column at the current node, returns INT_MAX
 */
SCIP_Longint SCIPcolGetStrongbranchLPAge(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(col != NULL);
   assert(stat != NULL);

   return (col->sbnode != stat->nnodes ? SCIP_LONGINT_MAX : stat->lpcount - stat->nsbdivinglps - col->validsblp);
}

/** marks a column to be not removable from the LP in the current node because it became obsolete */
void SCIPcolMarkNotRemovableLocal(
   SCIP_COL*             col,                /**< LP column */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(col  != NULL);
   assert(stat != NULL);
   assert(stat->nnodes > 0);

   /* lpRemoveObsoleteCols() does not remove a column if the node number stored in obsoletenode equals the current node number */
   col->obsoletenode = stat->nnodes;
}


/*
 * Row methods
 */

/** calculates row norms and min/maxidx from scratch, and checks for sorting */
static
void rowCalcNorms(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(row != NULL);
   assert(set != NULL);

   row->sqrnorm = 0.0;
   row->sumnorm = 0.0;
   row->objprod = 0.0;
   row->maxval = 0.0;
   row->nummaxval = 1;
   row->minval = SCIPsetInfinity(set);
   row->numminval = 1;
   row->minidx = INT_MAX;
   row->maxidx = INT_MIN;
   row->validminmaxidx = TRUE;
   row->lpcolssorted = TRUE;
   row->nonlpcolssorted = TRUE;

   /* check, if row is sorted
    * calculate sqrnorm, sumnorm, maxval, minval, minidx, and maxidx
    */
   for( i = 0; i < row->nlpcols; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(!SCIPsetIsZero(set, row->vals[i]));
      assert(row->cols[i]->lppos >= 0);
      assert(row->linkpos[i] >= 0);
      assert(row->cols[i]->index == row->cols_index[i]);

      rowAddNorms(row, set, row->cols[i], row->vals[i], TRUE);
      if( i > 0 )
      {
         assert(row->cols[i-1]->index == row->cols_index[i-1]);
         row->lpcolssorted = row->lpcolssorted && (row->cols_index[i-1] < row->cols_index[i]);
      }
   }
   for( i = row->nlpcols; i < row->len; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(!SCIPsetIsZero(set, row->vals[i]));
      assert(row->cols[i]->lppos == -1 || row->linkpos[i] == -1);
      assert(row->cols[i]->index == row->cols_index[i]);

      rowAddNorms(row, set, row->cols[i], row->vals[i], TRUE);
      if( i > row->nlpcols )
      {
         assert(row->cols[i-1]->index == row->cols_index[i-1]);
         row->nonlpcolssorted = row->nonlpcolssorted && (row->cols_index[i-1] < row->cols_index[i]);
      }
   }
}

/** calculates min/maxval and min/maxidx from scratch */
static
void rowCalcIdxsAndVals(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_COL* col;
   SCIP_Real absval;
   int i;

   assert(row != NULL);
   assert(set != NULL);

   row->maxval = 0.0;
   row->nummaxval = 1;
   row->numintcols = 0;
   row->minval = SCIPsetInfinity(set);
   row->numminval = 1;
   row->minidx = INT_MAX;
   row->maxidx = INT_MIN;
   row->validminmaxidx = TRUE;

   /* calculate maxval, minval, minidx, and maxidx */
   for( i = 0; i < row->len; ++i )
   {
      col = row->cols[i];
      assert(col != NULL);
      assert(!SCIPsetIsZero(set, row->vals[i]));

      absval = REALABS(row->vals[i]);
      assert(!SCIPsetIsZero(set, absval));

      /* update min/maxidx */
      row->minidx = MIN(row->minidx, col->index);
      row->maxidx = MAX(row->maxidx, col->index);
      row->numintcols += SCIPcolIsIntegral(col); /*lint !e713*/

      /* update maximal and minimal non-zero value */
      if( row->nummaxval > 0 )
      {
         if( SCIPsetIsGT(set, absval, row->maxval) )
         {
            row->maxval = absval;
            row->nummaxval = 1;
         }
         else if( SCIPsetIsGE(set, absval, row->maxval) )
         {
            /* make sure the maxval is always exactly the same */
            row->maxval = MAX(absval, row->maxval);
            row->nummaxval++;
         }
      }
      if( row->numminval > 0 )
      {
         if( SCIPsetIsLT(set, absval, row->minval) )
         {
            row->minval = absval;
            row->numminval = 1;
         }
         else if( SCIPsetIsLE(set, absval, row->minval) )
         {
            /* make sure the minval is always exactly the same */
            row->minval = MIN(absval, row->minval);
            row->numminval++;
         }
      }
   }
}

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real*            intval              /**< pointer to store the scaled integral value, or NULL */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = floor(sval);
   upval = ceil(sval);

   if( SCIPrelDiff(sval, downval) <= maxdelta )
   {
      if( intval != NULL )
         *intval = downval;
      return TRUE;
   }
   else if( SCIPrelDiff(sval, upval) >= mindelta )
   {
      if( intval != NULL )
         *intval = upval;
      return TRUE;
   }

   return FALSE;
}

/** scales row with given factor, and rounds coefficients to integers if close enough;
 *  the constant is automatically moved to the sides;
 *  if the row's activity is proven to be integral, the sides are automatically rounded to the next integer
 */
static
SCIP_RETCODE rowScale(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             scaleval,           /**< value to scale row with */
   SCIP_Bool             integralcontvars,   /**< should the coefficients of the continuous variables also be made integral,
                                              *   if they are close to integral values? */
   SCIP_Real             minrounddelta,      /**< minimal relative difference of scaled coefficient s*c and integral i,
                                              *   upto which the integral is used instead of the scaled real coefficient */
   SCIP_Real             maxrounddelta       /**< maximal relative difference of scaled coefficient s*c and integral i
                                              *   upto which the integral is used instead of the scaled real coefficient */
   )
{
   SCIP_COL* col;
   SCIP_Real val;
   SCIP_Real newval;
   SCIP_Real intval;
   SCIP_Real mindelta;
   SCIP_Real maxdelta;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool mindeltainf;
   SCIP_Bool maxdeltainf;
   int oldlen; 
   int c;

   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->len == 0 || row->vals != NULL);
   assert(SCIPsetIsPositive(set, scaleval));
   assert(-1.0 < minrounddelta && minrounddelta <= 0.0);
   assert(0.0 <= maxrounddelta && maxrounddelta < 1.0);

   SCIPsetDebugMsg(set, "scale row <%s> with %g (tolerance=[%g,%g])\n", row->name, scaleval, minrounddelta, maxrounddelta);

   mindelta = 0.0;
   maxdelta = 0.0;
   mindeltainf = FALSE;
   maxdeltainf = FALSE;
   oldlen = row->len;

   /* scale the row coefficients, thereby recalculating whether the row's activity is always integral;
    * if the row coefficients are rounded to the nearest integer value, calculate the maximal activity difference,
    * this rounding can lead to
    */
   row->integral = TRUE;

   c = 0;
   while( c < row->len )
   {
      col = row->cols[c];
      val = row->vals[c];
      assert(!SCIPsetIsZero(set, val));

      /* get local or global bounds for column, depending on the local or global feasibility of the row */
      if( row->local )
      {
         lb = col->lb;
         ub = col->ub;
      }
      else
      {
         lb = SCIPvarGetLbGlobal(col->var);
         ub = SCIPvarGetUbGlobal(col->var);
      }

      /* calculate scaled coefficient */
      newval = val * scaleval;
      if( (integralcontvars || SCIPcolIsIntegral(col) || SCIPsetIsIntegral(set, newval))
         && isIntegralScalar(val, scaleval, minrounddelta, maxrounddelta, &intval) )
      {
         if( !SCIPsetIsEQ(set, intval, newval) )
         {
            if( intval < newval )
            {
               mindelta += (intval - newval)*ub;
               maxdelta += (intval - newval)*lb;
               mindeltainf = mindeltainf || SCIPsetIsInfinity(set, ub);
               maxdeltainf = maxdeltainf || SCIPsetIsInfinity(set, -lb);
            }
            else
            {
               mindelta += (intval - newval)*lb;
               maxdelta += (intval - newval)*ub;
               mindeltainf = mindeltainf || SCIPsetIsInfinity(set, -lb);
               maxdeltainf = maxdeltainf || SCIPsetIsInfinity(set, ub);
            }
         }
         newval = intval;
      }

      if( !SCIPsetIsEQ(set, val, newval) )
      {
         /* if column knows of the row, change the corresponding coefficient in the column */
         if( row->linkpos[c] >= 0 )
         {
            assert(col->rows[row->linkpos[c]] == row);
            assert(SCIPsetIsEQ(set, col->vals[row->linkpos[c]], row->vals[c]));
            SCIP_CALL( colChgCoefPos(col, set, lp, row->linkpos[c], newval) );
         }

         /* change the coefficient in the row, and update the norms and integrality status */
         SCIP_CALL( rowChgCoefPos(row, blkmem, set, eventqueue, lp, c, newval) );

         /* current coefficient has been deleted from the row because it was almost zero */
         if( oldlen != row->len ) 
         {  
            assert(row->len == oldlen - 1);
            c--;
            oldlen = row->len;
         }

      }
      else
         row->integral = row->integral && SCIPcolIsIntegral(col) && SCIPsetIsIntegral(set, val);

      ++c;
   }

   /* scale the row sides, and move the constant to the sides; relax the sides with accumulated delta in order
    * to not destroy feasibility due to rounding
    */
   /**@todo ensure that returned cut does not have infinite lhs and rhs */
   if( !SCIPsetIsInfinity(set, -row->lhs) )
   {
      if( mindeltainf )
         newval = -SCIPsetInfinity(set);
      else
      {
         newval = (row->lhs - row->constant) * scaleval + mindelta;
         if( SCIPsetIsIntegral(set, newval) || (row->integral && !row->modifiable) )
            newval = SCIPsetSumCeil(set, newval);
      }
      SCIP_CALL( SCIProwChgLhs(row, blkmem, set, eventqueue, lp, newval) );
   }
   if( !SCIPsetIsInfinity(set, row->rhs) )
   {
      if( maxdeltainf )
         newval = SCIPsetInfinity(set);
      else
      {
         newval = (row->rhs - row->constant) * scaleval + maxdelta;
         if( SCIPsetIsIntegral(set, newval) || (row->integral && !row->modifiable) )
            newval = SCIPsetSumFloor(set, newval);
      }
      SCIP_CALL( SCIProwChgRhs(row, blkmem, set, eventqueue, lp, newval) );
   }

   /* clear the row constant */
   SCIP_CALL( SCIProwChgConstant(row, blkmem, set, stat, eventqueue, lp, 0.0) );

   SCIPsetDebugMsg(set, "scaled row <%s> (integral: %u)\n", row->name, row->integral);
   debugRowPrint(set, row);

#ifdef SCIP_DEBUG
   /* check integrality status of row */
   for( c = 0; c < row->len && SCIPcolIsIntegral(row->cols[c]) && SCIPsetIsIntegral(set, row->vals[c]); ++c )
   {}
   assert(row->integral == (c == row->len));
#endif

   /* invalid the activity */
   row->validactivitylp = -1;

   return SCIP_OKAY;
}

/** creates and captures an LP row */
SCIP_RETCODE SCIProwCreate(
   SCIP_ROW**            row,                /**< pointer to LP row data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_ROWORIGINTYPE    origintype,         /**< type of origin of row */
   void*                 origin,             /**< pointer to constraint handler or separator who created the row (NULL if unkown) */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   assert(row != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(len >= 0);
   assert(len == 0 || (cols != NULL && vals != NULL));
   /* note, that the assert tries to avoid numerical troubles in the LP solver.
    * in case, for example, lhs > rhs but they are equal with tolerances, one could pass lhs=rhs=lhs+rhs/2 to
    * SCIProwCreate() (see cons_linear.c: detectRedundantConstraints())
    */
   assert(lhs <= rhs);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, row) );

   (*row)->integral = TRUE;
   if( len > 0 )
   {
      SCIP_VAR* var;
      int i;

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->cols, cols, len) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->vals, vals, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->cols_index, len) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*row)->linkpos, len) );

      for( i = 0; i < len; ++i )
      {
         assert(cols[i] != NULL);
         assert(!SCIPsetIsZero(set, vals[i]));

         var = cols[i]->var;
         (*row)->cols_index[i] = cols[i]->index;
         (*row)->linkpos[i] = -1;
         if( SCIPsetIsIntegral(set, (*row)->vals[i]) )
         {
            (*row)->vals[i] = SCIPsetRound(set, (*row)->vals[i]);
            (*row)->integral = (*row)->integral && SCIPvarIsIntegral(var);
         }
         else
         {
            (*row)->integral = FALSE;
         }
      }
   }
   else
   {
      (*row)->cols = NULL;
      (*row)->cols_index = NULL;
      (*row)->vals = NULL;
      (*row)->linkpos = NULL;
   }

   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*row)->name, name, strlen(name)+1) );
   (*row)->constant = 0.0;
   (*row)->lhs = lhs;
   (*row)->rhs = rhs;
   (*row)->flushedlhs = -SCIPsetInfinity(set);
   (*row)->flushedrhs = SCIPsetInfinity(set);
   (*row)->sqrnorm = 0.0;
   (*row)->sumnorm = 0.0;
   (*row)->objprod = 0.0;
   (*row)->maxval = 0.0;
   (*row)->minval = SCIPsetInfinity(set);
   (*row)->dualsol = 0.0;
   (*row)->activity = SCIP_INVALID;
   (*row)->dualfarkas = 0.0;
   (*row)->pseudoactivity = SCIP_INVALID;
   (*row)->minactivity = SCIP_INVALID;
   (*row)->maxactivity = SCIP_INVALID;
   (*row)->origin = origin;
   (*row)->eventfilter = NULL;
   (*row)->index = stat->nrowidx;
   SCIPstatIncrement(stat, set, nrowidx);
   (*row)->size = len;
   (*row)->len = len;
   (*row)->nlpcols = 0;
   (*row)->nunlinked = len;
   (*row)->nuses = 0;
   (*row)->lppos = -1;
   (*row)->lpipos = -1;
   (*row)->lpdepth = -1;
   (*row)->minidx = INT_MAX;
   (*row)->maxidx = INT_MIN;
   (*row)->nummaxval = 0;
   (*row)->numminval = 0;
   (*row)->numintcols = -1;
   (*row)->validactivitylp = -1;
   (*row)->validpsactivitydomchg = -1;
   (*row)->validactivitybdsdomchg = -1;
   (*row)->nlpsaftercreation = 0L;
   (*row)->activeinlpcounter = 0L;
   (*row)->age = 0;
   (*row)->rank = 0;
   (*row)->obsoletenode = -1;
   (*row)->basisstatus = SCIP_BASESTAT_BASIC; /*lint !e641*/
   (*row)->lpcolssorted = TRUE;
   (*row)->nonlpcolssorted = (len <= 1);
   (*row)->delaysort = FALSE;
   (*row)->validminmaxidx = FALSE;
   (*row)->lhschanged = FALSE;
   (*row)->rhschanged = FALSE;
   (*row)->coefchanged = FALSE;
   (*row)->local = local;
   (*row)->modifiable = modifiable;
   (*row)->nlocks = 0;
   (*row)->origintype = origintype; /*lint !e641*/
   (*row)->removable = removable;
   (*row)->inglobalcutpool = FALSE;
   (*row)->storedsolvals = NULL;

   /* calculate row norms and min/maxidx, and check if row is sorted */
   rowCalcNorms(*row, set);

   /* capture the row */
   SCIProwCapture(*row);

   /* create event filter */
   SCIP_CALL( SCIPeventfilterCreate(&(*row)->eventfilter, blkmem) );

   return SCIP_OKAY;
} /*lint !e715*/

/** frees an LP row */
SCIP_RETCODE SCIProwFree(
   SCIP_ROW**            row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses == 0);
   assert((*row)->lppos == -1);
   assert((*row)->eventfilter != NULL);

   /* remove column indices from corresponding rows */
   SCIP_CALL( rowUnlink(*row, set, lp) );

   /* free event filter */
   SCIP_CALL( SCIPeventfilterFree(&(*row)->eventfilter, blkmem, set) );

   BMSfreeBlockMemoryNull(blkmem, &(*row)->storedsolvals);
   BMSfreeBlockMemoryArray(blkmem, &(*row)->name, strlen((*row)->name)+1);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->cols, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->cols_index, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->vals, (*row)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*row)->linkpos, (*row)->size);
   BMSfreeBlockMemory(blkmem, row);

   return SCIP_OKAY;
}

/** output row to file stream */
void SCIProwPrint(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   assert(row != NULL);

   /* print row name */
   if( row->name != NULL && row->name[0] != '\0' )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, "%s: ", row->name);
   }

   /* print left hand side */
   SCIPmessageFPrintInfo(messagehdlr, file, "%.15g <= ", row->lhs);

   /* print coefficients */
   if( row->len == 0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "0 ");
   for( i = 0; i < row->len; ++i )
   {
      assert(row->cols[i] != NULL);
      assert(row->cols[i]->var != NULL);
      assert(SCIPvarGetName(row->cols[i]->var) != NULL);
      assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.15g<%s> ", row->vals[i], SCIPvarGetName(row->cols[i]->var));
   }

   /* print constant */
   if( REALABS(row->constant) > SCIP_DEFAULT_EPSILON )
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.15g ", row->constant);

   /* print right hand side */
   SCIPmessageFPrintInfo(messagehdlr, file, "<= %.15g\n", row->rhs);
}

/** increases usage counter of LP row */
void SCIProwCapture(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->nuses >= 0);
   assert(row->nlocks <= (unsigned int)(row->nuses)); /*lint !e574*/

   SCIPdebugMessage("capture row <%s> with nuses=%d and nlocks=%u\n", row->name, row->nuses, row->nlocks);
   row->nuses++;
}

/** decreases usage counter of LP row, and frees memory if necessary */
SCIP_RETCODE SCIProwRelease(
   SCIP_ROW**            row,                /**< pointer to LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(blkmem != NULL);
   assert(row != NULL);
   assert(*row != NULL);
   assert((*row)->nuses >= 1);
   assert((*row)->nlocks < (unsigned int)((*row)->nuses)); /*lint !e574*/

   SCIPsetDebugMsg(set, "release row <%s> with nuses=%d and nlocks=%u\n", (*row)->name, (*row)->nuses, (*row)->nlocks);
   (*row)->nuses--;
   if( (*row)->nuses == 0 )
   {
      SCIP_CALL( SCIProwFree(row, blkmem, set, lp) );
   }

   *row = NULL;

   return SCIP_OKAY;
}

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
void SCIProwLock(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   /* check, if row is modifiable */
   if( !row->modifiable )
   {
      SCIPdebugMessage("lock row <%s> with nuses=%d and nlocks=%u\n", row->name, row->nuses, row->nlocks);
      row->nlocks++;
   }
}

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
void SCIProwUnlock(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   /* check, if row is modifiable */
   if( !row->modifiable )
   {
      SCIPdebugMessage("unlock row <%s> with nuses=%d and nlocks=%u\n", row->name, row->nuses, row->nlocks);
      assert(row->nlocks > 0);
      row->nlocks--;
   }
}

/** adds a previously non existing coefficient to an LP row */
SCIP_RETCODE SCIProwAddCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   assert(lp != NULL);
   assert(!lp->diving || row->lppos == -1);

   SCIP_CALL( rowAddCoef(row, blkmem, set, eventqueue, lp, col, val, -1) );

   checkLinks(lp);

   return SCIP_OKAY;
}

/** deletes coefficient from row */
SCIP_RETCODE SCIProwDelCoef(
   SCIP_ROW*             row,                /**< row to be changed */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(row != NULL);
   assert(!row->delaysort);
   assert(lp != NULL);
   assert(!lp->diving || row->lppos == -1);
   assert(col != NULL);
   assert(col->var != NULL);

   /* search the position of the column in the row's col vector */
   pos = rowSearchCoef(row, col);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for column <%s> doesn't exist in row <%s>\n", SCIPvarGetName(col->var), row->name);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < row->len);
   assert(row->cols[pos] == col);
   assert(row->cols_index[pos] == col->index);

   /* if column knows of the row, remove the row from the column's row vector */
   if( row->linkpos[pos] >= 0 )
   {
      assert(col->rows[row->linkpos[pos]] == row);
      assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
      SCIP_CALL( colDelCoefPos(col, set, lp, row->linkpos[pos]) );
   }

   /* delete the column from the row's col vector */
   SCIP_CALL( rowDelCoefPos(row, blkmem, set, eventqueue, lp, pos) );

   checkLinks(lp);

   return SCIP_OKAY;
}

/** changes or adds a coefficient to an LP row */
SCIP_RETCODE SCIProwChgCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   int pos;

   assert(row != NULL);
   assert(!row->delaysort);
   assert(lp != NULL);
   assert(!lp->diving || row->lppos == -1);
   assert(col != NULL);

   /* search the position of the column in the row's col vector */
   pos = rowSearchCoef(row, col);

   /* check, if column already exists in the row's col vector */
   if( pos == -1 )
   {
      /* add previously not existing coefficient */
      SCIP_CALL( rowAddCoef(row, blkmem, set, eventqueue, lp, col, val, -1) );
   }
   else
   {
      /* modify already existing coefficient */
      assert(0 <= pos && pos < row->len);
      assert(row->cols[pos] == col);
      assert(row->cols_index[pos] == col->index);

      /* if column knows of the row, change the corresponding coefficient in the column */
      if( row->linkpos[pos] >= 0 )
      {
         assert(col->rows[row->linkpos[pos]] == row);
         assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
         SCIP_CALL( colChgCoefPos(col, set, lp, row->linkpos[pos], val) );
      }

      /* change the coefficient in the row */
      SCIP_CALL( rowChgCoefPos(row, blkmem, set, eventqueue, lp, pos, val) );
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** increases value of an existing or non-existing coefficient in an LP row */
SCIP_RETCODE SCIProwIncCoef(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_COL*             col,                /**< LP column */
   SCIP_Real             incval              /**< value to add to the coefficient */
   )
{
   int pos;

   assert(row != NULL);
   assert(lp != NULL);
   assert(!lp->diving || row->lppos == -1);
   assert(col != NULL);

   if( SCIPsetIsZero(set, incval) )
      return SCIP_OKAY;

   /* search the position of the column in the row's col vector */
   pos = rowSearchCoef(row, col);

   /* check, if column already exists in the row's col vector */
   if( pos == -1 )
   {
      /* coefficient doesn't exist, or sorting is delayed: add coefficient to the end of the row's arrays */
      SCIP_CALL( rowAddCoef(row, blkmem, set, eventqueue, lp, col, incval, -1) );
   }
   else
   {
      /* modify already existing coefficient */
      assert(0 <= pos && pos < row->len);
      assert(row->cols[pos] == col);
      assert(row->cols_index[pos] == col->index);

      /* if column knows of the row, change the corresponding coefficient in the column */
      if( row->linkpos[pos] >= 0 )
      {
         assert(col->rows[row->linkpos[pos]] == row);
         assert(SCIPsetIsEQ(set, col->vals[row->linkpos[pos]], row->vals[pos]));
         SCIP_CALL( colChgCoefPos(col, set, lp, row->linkpos[pos], row->vals[pos] + incval) );
      }

      /* change the coefficient in the row */
      SCIP_CALL( rowChgCoefPos(row, blkmem, set, eventqueue, lp, pos, row->vals[pos] + incval) );
   }

   checkLinks(lp);

   /* invalid the activity */
   row->validactivitylp = -1;

   return SCIP_OKAY;
}

/** changes constant value of a row */
SCIP_RETCODE SCIProwChgConstant(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             constant            /**< new constant value */
   )
{
   assert(row != NULL);
   assert(row->lhs <= row->rhs);
   assert(!SCIPsetIsInfinity(set, REALABS(constant)));
   assert(stat != NULL);
   assert(lp != NULL);
   assert(!lp->diving || row->lppos == -1);

   if( !SCIPsetIsEQ(set, constant, row->constant) )
   {
      SCIP_Real oldconstant;

      if( row->validpsactivitydomchg == stat->domchgcount )
      {
         assert(row->pseudoactivity < SCIP_INVALID);
         row->pseudoactivity += constant - row->constant;
      }
      if( row->validactivitybdsdomchg == stat->domchgcount )
      {
         assert(row->minactivity < SCIP_INVALID);
         assert(row->maxactivity < SCIP_INVALID);
         row->minactivity += constant - row->constant;
         row->maxactivity += constant - row->constant;
      }

      if( !SCIPsetIsInfinity(set, -row->lhs) )
      {
         SCIP_CALL( rowSideChanged(row, set, lp, SCIP_SIDETYPE_LEFT) );
      }
      if( !SCIPsetIsInfinity(set, row->rhs) )
      {
         SCIP_CALL( rowSideChanged(row, set, lp, SCIP_SIDETYPE_RIGHT) );
      }

      oldconstant = row->constant;

      row->constant = constant;

      /* issue row constant changed event */
      SCIP_CALL( rowEventConstantChanged(row, blkmem, set, eventqueue, oldconstant, constant) );
   }

   return SCIP_OKAY;
}

/** add constant value to a row */
SCIP_RETCODE SCIProwAddConstant(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             addval              /**< constant value to add to the row */
   )
{
   assert(row != NULL);
   assert(row->lhs <= row->rhs);
   assert(!SCIPsetIsInfinity(set, REALABS(addval)));
   assert(stat != NULL);
   assert(lp != NULL);
   assert(!lp->diving || row->lppos == -1);

   if( !SCIPsetIsZero(set, addval) )
   {
      SCIP_CALL( SCIProwChgConstant(row, blkmem, set, stat, eventqueue, lp, row->constant + addval) );
   }

   return SCIP_OKAY;
}

/** changes left hand side of LP row */
SCIP_RETCODE SCIProwChgLhs(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);

   if( !SCIPsetIsEQ(set, row->lhs, lhs) )
   {
      SCIP_Real oldlhs;

      oldlhs = row->lhs;

      row->lhs = lhs;
      SCIP_CALL( rowSideChanged(row, set, lp, SCIP_SIDETYPE_LEFT) );

      if( !lp->diving )
      {
         /* issue row side changed event */
         SCIP_CALL( rowEventSideChanged(row, blkmem, set, eventqueue, SCIP_SIDETYPE_LEFT, oldlhs, lhs) );
      }
   }

   return SCIP_OKAY;
}

/** changes right hand side of LP row */
SCIP_RETCODE SCIProwChgRhs(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   assert(row != NULL);
   assert(lp != NULL);

   if( !SCIPsetIsEQ(set, row->rhs, rhs) )
   {
      SCIP_Real oldrhs;

      oldrhs = row->rhs;

      row->rhs = rhs;
      SCIP_CALL( rowSideChanged(row, set, lp, SCIP_SIDETYPE_RIGHT) );

      if( !lp->diving )
      {
         /* issue row side changed event */
         SCIP_CALL( rowEventSideChanged(row, blkmem, set, eventqueue, SCIP_SIDETYPE_RIGHT, oldrhs, rhs) );
      }
   }

   return SCIP_OKAY;
}

/** changes the local flag of LP row */
SCIP_RETCODE SCIProwChgLocal(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Bool             local               /**< new value for local flag */
   )
{
   assert(row != NULL);

   row->local = local;

   return SCIP_OKAY;
}

/** additional scalars that are tried in integrality scaling */
static const SCIP_Real scalars[] = {3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};
static const int nscalars = 9;

/** tries to find a value, such that all row coefficients, if scaled with this value become integral */
SCIP_RETCODE SCIProwCalcIntegralScalar(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   )
{
#ifndef NDEBUG
   SCIP_COL* col;
#endif
   SCIP_Longint gcd;
   SCIP_Longint scm;
   SCIP_Longint nominator;
   SCIP_Longint denominator;
   SCIP_Real val;
   SCIP_Real absval;
   SCIP_Real minval;
   SCIP_Real scaleval;
   SCIP_Real twomultval;
   SCIP_Bool scalable;
   SCIP_Bool twomult;
   SCIP_Bool rational;
   int c;
   int s;

   /**@todo call misc.c:SCIPcalcIntegralScalar() instead - if usecontvars == FALSE, filter the integer variables first */
   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->len == 0 || row->cols_index != NULL);
   assert(row->len == 0 || row->vals != NULL);
   assert(maxdnom >= 1);
   assert(mindelta < 0.0);
   assert(maxdelta > 0.0);
   assert(success != NULL);

   SCIPsetDebugMsg(set, "trying to find rational representation for row <%s> (contvars: %u)\n", SCIProwGetName(row), usecontvars);
   SCIPdebug( val = 0; ); /* avoid warning "val might be used uninitialized; see SCIPdebugMessage lastval=%g below */

   if( intscalar != NULL )
      *intscalar = SCIP_INVALID;
   *success = FALSE;

   /* get minimal absolute non-zero value */
   minval = SCIP_REAL_MAX;
   for( c = 0; c < row->len; ++c )
   {
#ifndef NDEBUG
      col = row->cols[c];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);
#endif
      val = row->vals[c];
      assert(!SCIPsetIsZero(set, val));

      if( val < mindelta || val > maxdelta )
      {
         absval = REALABS(val);
         minval = MIN(minval, absval);
      }
   }
   if( minval == SCIP_REAL_MAX ) /*lint !e777*/
   {
      /* all coefficients are zero (inside tolerances) */
      if( intscalar != NULL )
         *intscalar = 1.0;
      *success = TRUE;
      SCIPsetDebugMsg(set, " -> all values are zero (inside tolerances)\n");

      return SCIP_OKAY;
   }
   assert(minval > MIN(-mindelta, maxdelta)); 
   assert(SCIPsetIsPositive(set, minval));
   assert(!SCIPsetIsInfinity(set, minval));

   /* try, if row coefficients can be made integral by multiplying them with the reciprocal of the smallest coefficient
    * and a power of 2
    */
   scaleval = 1.0/minval;
   scalable = (scaleval <= maxscale);
   for( c = 0; c < row->len && scalable; ++c )
   {
      /* don't look at continuous variables, if we don't have to */
      if( !usecontvars && !SCIPcolIsIntegral(row->cols[c]) )
         continue;

      /* check, if the coefficient can be scaled with a simple scalar */
      val = row->vals[c];
      absval = REALABS(val);
      while( scaleval <= maxscale
         && (absval * scaleval < 0.5 || !isIntegralScalar(val, scaleval, mindelta, maxdelta, NULL)) )
      {
         for( s = 0; s < nscalars; ++s )
         {
            if( isIntegralScalar(val, scaleval * scalars[s], mindelta, maxdelta, NULL) )
            {
               scaleval *= scalars[s];
               break;
            }
         }
         if( s >= nscalars )
            scaleval *= 2.0;
      }
      scalable = (scaleval <= maxscale);
      SCIPsetDebugMsg(set, " -> val=%g, scaleval=%g, val*scaleval=%g, scalable=%u\n", val, scaleval, val*scaleval, scalable);
   }
   if( scalable )
   {
      /* make row coefficients integral by dividing them by the smallest coefficient
       * (and multiplying them with a power of 2)
       */
      assert(scaleval <= maxscale);
      if( intscalar != NULL )
         *intscalar = scaleval;
      *success = TRUE;
      SCIPsetDebugMsg(set, " -> integrality can be achieved by scaling with %g (minval=%g)\n", scaleval, minval);

      return SCIP_OKAY;
   }

   /* try, if row coefficients can be made integral by multiplying them by a power of 2 */
   twomultval = 1.0;
   twomult = (twomultval <= maxscale);
   for( c = 0; c < row->len && twomult; ++c )
   {
      /* don't look at continuous variables, if we don't have to */
      if( !usecontvars && !SCIPcolIsIntegral(row->cols[c]) )
         continue;

      /* check, if the coefficient can be scaled with a simple scalar */
      val = row->vals[c];
      absval = REALABS(val);
      while( twomultval <= maxscale
         && (absval * twomultval < 0.5 || !isIntegralScalar(val, twomultval, mindelta, maxdelta, NULL)) )
      {
         for( s = 0; s < nscalars; ++s )
         {
            if( isIntegralScalar(val, twomultval * scalars[s], mindelta, maxdelta, NULL) )
            {
               twomultval *= scalars[s];
               break;
            }
         }
         if( s >= nscalars )
            twomultval *= 2.0;
      }
      twomult = (twomultval <= maxscale);
      SCIPsetDebugMsg(set, " -> val=%g, twomult=%g, val*twomult=%g, twomultable=%u\n",
         val, twomultval, val*twomultval, twomult);
   }
   if( twomult )
   {
      /* make row coefficients integral by multiplying them with a power of 2 */
      assert(twomultval <= maxscale);
      if( intscalar != NULL )
         *intscalar = twomultval;
      *success = TRUE;
      SCIPsetDebugMsg(set, " -> integrality can be achieved by scaling with %g (power of 2)\n", twomultval);

      return SCIP_OKAY;
   }

   /* convert each coefficient into a rational number, calculate the greatest common divisor of the numerators
    * and the smallest common multiple of the denominators
    */
   gcd = 1;
   scm = 1;
   rational = (maxdnom > 1);

   /* first coefficient (to initialize gcd) */
   for( c = 0; c < row->len && rational; ++c )
   {
      if( usecontvars || SCIPcolIsIntegral(row->cols[c]) )
      {
         val = row->vals[c];
         rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
         if( rational && nominator != 0 )
         {
            assert(denominator > 0);
            gcd = ABS(nominator);
            scm = denominator;
            rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
            SCIPsetDebugMsg(set, " -> first rational: val: %g == %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ", gcd=%" SCIP_LONGINT_FORMAT ", scm=%" SCIP_LONGINT_FORMAT ", rational=%u\n",
               val, nominator, denominator, gcd, scm, rational);
            break;
         }
      }
   }

   /* remaining coefficients */
   for( ++c; c < row->len && rational; ++c )
   {
      if( usecontvars || SCIPcolIsIntegral(row->cols[c]) )
      {
         val = row->vals[c];
         rational = SCIPrealToRational(val, mindelta, maxdelta, maxdnom, &nominator, &denominator);
         if( rational && nominator != 0 )
         {
            assert(denominator > 0);
            gcd = SCIPcalcGreComDiv(gcd, ABS(nominator));
            scm *= denominator / SCIPcalcGreComDiv(scm, denominator);
            rational = ((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
            SCIPsetDebugMsg(set, " -> next rational : val: %g == %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ", gcd=%" SCIP_LONGINT_FORMAT ", scm=%" SCIP_LONGINT_FORMAT ", rational=%u\n",
               val, nominator, denominator, gcd, scm, rational);
         }
      }
   }

   if( rational )
   {
      /* make row coefficients integral by multiplying them with the smallest common multiple of the denominators */
      assert((SCIP_Real)scm/(SCIP_Real)gcd <= maxscale);
      if( intscalar != NULL )
         *intscalar = (SCIP_Real)scm/(SCIP_Real)gcd;
      *success = TRUE;
      SCIPsetDebugMsg(set, " -> integrality can be achieved by scaling with %g (rational:%" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ")\n",
         (SCIP_Real)scm/(SCIP_Real)gcd, scm, gcd);
   }
   else
   {
      assert(!(*success));
      SCIPsetDebugMsg(set, " -> rationalizing failed: gcd=%" SCIP_LONGINT_FORMAT ", scm=%" SCIP_LONGINT_FORMAT ", lastval=%g\n", gcd, scm, val); /*lint !e771*/
   }

   return SCIP_OKAY;
}

/** tries to scale row, s.t. all coefficients become integral */
SCIP_RETCODE SCIProwMakeIntegral(
   SCIP_ROW*             row,                /**< LP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal value to scale row with */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Bool*            success             /**< stores whether row could be made rational */
   )
{
   SCIP_Real intscalar;

   assert(success != NULL);

   /* calculate scalar to make coefficients integral */
   SCIP_CALL( SCIProwCalcIntegralScalar(row, set, mindelta, maxdelta, maxdnom, maxscale, usecontvars,
         &intscalar, success) );

   if( *success )
   {
      /* scale the row */
      SCIP_CALL( rowScale(row, blkmem, set, eventqueue, stat, lp, intscalar, usecontvars, mindelta, maxdelta) );
   }

   return SCIP_OKAY;
}

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
void SCIProwSort(
   SCIP_ROW*             row                 /**< row to be sorted */
   )
{
   assert(row != NULL);

   /* sort LP columns */
   rowSortLP(row);

   /* sort non-LP columns */
   rowSortNonLP(row);

#ifdef SCIP_MORE_DEBUG
   /* check the sorting */
   {
      int c;
      if( !row->delaysort )
      {
         for( c = 1; c < row->nlpcols; ++c )
            assert(row->cols[c]->index >= row->cols[c-1]->index);
         for( c = row->nlpcols + 1; c < row->len; ++c )
            assert(row->cols[c]->index >= row->cols[c-1]->index);
      }
   }
#endif
}

/** sorts row, and merges equal column entries (resulting from lazy sorting and adding) into a single entry; removes
 *  zero entries from row 
 *  the row must not be linked to the columns; otherwise, we would need to update the columns as
 *  well, which is too expensive
 */
static
void rowMerge(
   SCIP_ROW*             row,                /**< row to be sorted */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   assert(!row->delaysort);
   assert(row->nunlinked == row->len);
   assert(row->nlpcols == 0);

   SCIPsetDebugMsg(set, "merging row <%s>\n", row->name);

   /* do nothing on empty rows; if row is sorted, nothing has to be done */
   if( row->len > 0 && (!row->lpcolssorted || !row->nonlpcolssorted) )
   {
      SCIP_COL** cols;
      int* cols_index;
      SCIP_Real* vals;
      int s;
      int t;

      /* make sure, the row is sorted */
      SCIProwSort(row);
      assert(row->lpcolssorted);
      assert(row->nonlpcolssorted);

      /* merge equal columns, thereby recalculating whether the row's activity is always integral */
      cols = row->cols;
      cols_index = row->cols_index;
      vals = row->vals;
      assert(cols != NULL);
      assert(cols_index != NULL);
      assert(vals != NULL);

      t = 0;
      row->integral = TRUE;
      assert(!SCIPsetIsZero(set, vals[0]));
      assert(row->linkpos[0] == -1);

      for( s = 1; s < row->len; ++s )
      {
         assert(!SCIPsetIsZero(set, vals[s]));
         assert(row->linkpos[s] == -1);

         if( cols[s] == cols[t] )
         {
            /* merge entries with equal column */
            vals[t] += vals[s];
         }
         else
         {
            /* go to the next entry, overwriting current entry if coefficient is zero */
            if( !SCIPsetIsZero(set, vals[t]) )
            {
               /* in case the coefficient is integral w.r.t. numerics we explicitly round the coefficient to an integral value */
               vals[t] = SCIPsetIsIntegral(set, vals[t]) ? SCIPsetRound(set, vals[t]) : vals[t];

               row->integral = row->integral && SCIPcolIsIntegral(cols[t]) && SCIPsetIsIntegral(set, vals[t]);
               t++;
            }
            cols[t] = cols[s];
            cols_index[t] = cols_index[s];
            vals[t] = vals[s];
         }
      }
      if( !SCIPsetIsZero(set, vals[t]) )
      {
         row->integral = row->integral && SCIPcolIsIntegral(cols[t]) && SCIPsetIsIntegral(set, vals[t]);
         t++;
      }
      assert(s == row->len);
      assert(t <= row->len);

      row->len = t;
      row->nunlinked = t;

      /* if equal entries were merged, we have to recalculate the norms, since the squared Euclidean norm is wrong */
      if( t < s )
         rowCalcNorms(row, set);
   }

#ifndef NDEBUG
   /* check for double entries */
   {
      int i;
      int j;

      for( i = 0; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->index == row->cols_index[i]);
         for( j = i+1; j < row->len; ++j )
            assert(row->cols[i] != row->cols[j]);
      }
   }
#endif
}

/** enables delaying of row sorting */
void SCIProwDelaySort(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(!row->delaysort);

   row->delaysort = TRUE;
}

/** disables delaying of row sorting, sorts row and merges coefficients with equal columns */
void SCIProwForceSort(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);
   assert(row->delaysort);

   row->delaysort = FALSE;
   rowMerge(row, set);
}

/** recalculates the current activity of a row */
void SCIProwRecalcLPActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COL* col;
   int c;

   assert(row != NULL);
   assert(stat != NULL);

   row->activity = row->constant;
   for( c = 0; c < row->nlpcols; ++c )
   {
      col = row->cols[c];
      assert(col != NULL);
      assert(col->primsol < SCIP_INVALID);
      assert(col->lppos >= 0);
      assert(row->linkpos[c] >= 0);
      row->activity += row->vals[c] * col->primsol;
   }

   if( row->nunlinked > 0 )
   {
      for( c = row->nlpcols; c < row->len; ++c )
      {
         col = row->cols[c];
         assert(col != NULL);
         assert(col->lppos >= 0 || col->primsol == 0.0);
         assert(col->lppos == -1 || row->linkpos[c] == -1);
         if( col->lppos >= 0 )
            row->activity += row->vals[c] * col->primsol;
      }
   }
#ifndef NDEBUG
   else
   {
      for( c = row->nlpcols; c < row->len; ++c )
      {
         col = row->cols[c];
         assert(col != NULL);
         assert(col->primsol == 0.0);
         assert(col->lppos == -1);
         assert(row->linkpos[c] >= 0);
      }
   }
#endif

   row->validactivitylp = stat->lpcount;
}

/** returns the activity of a row in the current LP solution */
SCIP_Real SCIProwGetLPActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real inf;
   SCIP_Real activity;

   assert(row != NULL);
   assert(stat != NULL);
   assert(lp != NULL);
   assert(row->validactivitylp <= stat->lpcount);
   assert(lp->validsollp == stat->lpcount);

   if( row->validactivitylp != stat->lpcount )
      SCIProwRecalcLPActivity(row, stat);
   assert(row->validactivitylp == stat->lpcount);
   assert(row->activity < SCIP_INVALID);

   activity = row->activity;
   inf = SCIPsetInfinity(set);
   activity = MAX(activity, -inf);
   activity = MIN(activity, +inf);

   return activity;
}

/** returns the feasibility of a row in the current LP solution: negative value means infeasibility */
SCIP_Real SCIProwGetLPFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real activity;

   assert(row != NULL);

   activity = SCIProwGetLPActivity(row, set, stat, lp);

   return MIN(row->rhs - activity, activity - row->lhs);
}

/** returns the feasibility of a row in the relaxed solution solution: negative value means infeasibility
 *
 *  @todo Implement calculation of activities similar to LPs.
 */
SCIP_Real SCIProwGetRelaxFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Real inf;
   SCIP_Real activity;
   SCIP_COL* col;
   int c;

   assert( row != NULL );
   assert( stat != NULL );

   activity = row->constant;
   for (c = 0; c < row->nlpcols; ++c)
   {
      col = row->cols[c];
      assert( col != NULL );
      assert( col->lppos >= 0 );
      assert( col->var != NULL );
      assert( row->linkpos[c] >= 0 );
      activity += row->vals[c] * SCIPvarGetRelaxSol(col->var, set);
   }

   if ( row->nunlinked > 0 )
   {
      for (c = row->nlpcols; c < row->len; ++c)
      {
         col = row->cols[c];
         assert( col != NULL );
         assert( col->lppos == -1 || row->linkpos[c] == -1 );
         if ( col->lppos >= 0 )
         {
            assert( col->var != NULL );
            activity += row->vals[c] * SCIPvarGetRelaxSol(col->var, set);
         }
      }
   }
#ifndef NDEBUG
   else
   {
      for (c = row->nlpcols; c < row->len; ++c)
      {
         col = row->cols[c];
         assert( col != NULL );
         assert( col->lppos == -1 );
         assert( row->linkpos[c] >= 0 );
      }
   }
#endif
   inf = SCIPsetInfinity(set);
   activity = MAX(activity, -inf);
   activity = MIN(activity, +inf);

   return MIN(row->rhs - activity, activity - row->lhs);
}

/** returns the feasibility of a row in the current NLP solution: negative value means infeasibility
 *
 *  @todo Implement calculation of activities similar to LPs.
 */
SCIP_Real SCIProwGetNLPFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Real inf;
   SCIP_Real activity;
   SCIP_COL* col;
   int c;

   assert( row != NULL );
   assert( stat != NULL );

   activity = row->constant;
   for (c = 0; c < row->nlpcols; ++c)
   {
      col = row->cols[c];
      assert( col != NULL );
      assert( col->lppos >= 0 );
      assert( col->var != NULL );
      assert( row->linkpos[c] >= 0 );
      activity += row->vals[c] * SCIPvarGetNLPSol(col->var);
   }

   if ( row->nunlinked > 0 )
   {
      for (c = row->nlpcols; c < row->len; ++c)
      {
         col = row->cols[c];
         assert( col != NULL );
         assert( col->lppos == -1 || row->linkpos[c] == -1 );
         if ( col->lppos >= 0 )
         {
            assert( col->var != NULL );
            activity += row->vals[c] * SCIPvarGetNLPSol(col->var);
         }
      }
   }
#ifndef NDEBUG
   else
   {
      for (c = row->nlpcols; c < row->len; ++c)
      {
         col = row->cols[c];
         assert( col != NULL );
         assert( col->lppos == -1 );
         assert( row->linkpos[c] >= 0 );
      }
   }
#endif
   inf = SCIPsetInfinity(set);
   activity = MAX(activity, -inf);
   activity = MIN(activity, +inf);

   return MIN(row->rhs - activity, activity - row->lhs);
}

/** calculates the current pseudo activity of a row */
void SCIProwRecalcPseudoActivity(
   SCIP_ROW*             row,                /**< row data */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COL* col;
   int i;

   assert(row != NULL);
   assert(stat != NULL);

   row->pseudoactivity = row->constant;
   for( i = 0; i < row->len; ++i )
   {
      col = row->cols[i];
      assert(col != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0 && col->lppos >= 0));
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);

      row->pseudoactivity += SCIPcolGetBestBound(col) * row->vals[i];
   }
   row->validpsactivitydomchg = stat->domchgcount;
   assert(!row->integral || EPSISINT(row->pseudoactivity - row->constant, SCIP_DEFAULT_SUMEPSILON));
}

/** returns the pseudo activity of a row in the current pseudo solution */
SCIP_Real SCIProwGetPseudoActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Real inf;
   SCIP_Real activity;

   assert(row != NULL);
   assert(stat != NULL);
   assert(row->validpsactivitydomchg <= stat->domchgcount);

   /* check, if pseudo activity has to be calculated */
   if( row->validpsactivitydomchg != stat->domchgcount )
      SCIProwRecalcPseudoActivity(row, stat);
   assert(row->validpsactivitydomchg == stat->domchgcount);
   assert(row->pseudoactivity < SCIP_INVALID);

   activity = row->pseudoactivity;
   inf = SCIPsetInfinity(set);
   activity = MAX(activity, -inf);
   activity = MIN(activity, +inf);

   return activity;
}

/** returns the pseudo feasibility of a row in the current pseudo solution: negative value means infeasibility */
SCIP_Real SCIProwGetPseudoFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Real pseudoactivity;

   assert(row != NULL);

   pseudoactivity = SCIProwGetPseudoActivity(row, set, stat);

   return MIN(row->rhs - pseudoactivity, pseudoactivity - row->lhs);
}

/** returns the activity of a row for a given solution */
SCIP_Real SCIProwGetSolActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_COL* col;
   SCIP_Real inf;
   SCIP_Real activity;
   SCIP_Real solval;
   int i;

   assert(row != NULL);

   activity = row->constant;
   for( i = 0; i < row->len; ++i )
   {
      col = row->cols[i];
      assert(col != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0 && col->lppos >= 0));
      solval = SCIPsolGetVal(sol, set, stat, col->var);
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
      {
         if( SCIPsetIsInfinity(set, -row->lhs) )
            solval = (row->vals[i] >= 0.0 ? col->lb : col->ub);
         else if( SCIPsetIsInfinity(set, row->rhs) )
            solval = (row->vals[i] >= 0.0 ? col->ub : col->lb);
         else
            solval = (col->lb + col->ub)/2.0;
      }
      activity += row->vals[i] * solval;
   }

   inf = SCIPsetInfinity(set);
   activity = MAX(activity, -inf);
   activity = MIN(activity, +inf);

   return activity;
}

/** returns the feasibility of a row for the given solution */
SCIP_Real SCIProwGetSolFeasibility(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_Real activity;

   assert(row != NULL);

   activity = SCIProwGetSolActivity(row, set, stat, sol);

   return MIN(row->rhs - activity, activity - row->lhs);
}

/** calculates minimal and maximal activity of row w.r.t. the column's bounds */
static
void rowCalcActivityBounds(
   SCIP_ROW*             row,                /**< row data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   SCIP_COL* col;
   SCIP_Real val;
   SCIP_Bool mininfinite;
   SCIP_Bool maxinfinite;
   int i;

   assert(row != NULL);
   assert(!SCIPsetIsInfinity(set, REALABS(row->constant)));
   assert(stat != NULL);

   /* calculate activity bounds */
   mininfinite = FALSE;
   maxinfinite = FALSE;
   row->minactivity = row->constant;
   row->maxactivity = row->constant;
   for( i = 0; i < row->len && (!mininfinite || !maxinfinite); ++i )
   {
      col = row->cols[i];
      assert(col != NULL);
      assert((i < row->nlpcols) == (row->linkpos[i] >= 0 && col->lppos >= 0));
      val = row->vals[i];
      if( val >= 0.0 )
      {
         mininfinite = mininfinite || SCIPsetIsInfinity(set, -col->lb);
         maxinfinite = maxinfinite || SCIPsetIsInfinity(set, col->ub);
         if( !mininfinite )
            row->minactivity += val * col->lb;
         if( !maxinfinite )
            row->maxactivity += val * col->ub;
      }
      else
      {
         mininfinite = mininfinite || SCIPsetIsInfinity(set, col->ub);
         maxinfinite = maxinfinite || SCIPsetIsInfinity(set, -col->lb);
         if( !mininfinite )
            row->minactivity += val * col->ub;
         if( !maxinfinite )
            row->maxactivity += val * col->lb;
      }
   }

   if( mininfinite )
      row->minactivity = -SCIPsetInfinity(set);
   if( maxinfinite )
      row->maxactivity = SCIPsetInfinity(set);
   row->validactivitybdsdomchg = stat->domchgcount;

   assert(!row->integral || mininfinite || REALABS(row->minactivity - row->constant) > 1.0/SCIP_DEFAULT_SUMEPSILON
      || EPSISINT(row->minactivity - row->constant, SCIP_DEFAULT_SUMEPSILON));
   assert(!row->integral || maxinfinite || REALABS(row->maxactivity - row->constant) > 1.0/SCIP_DEFAULT_SUMEPSILON
      || EPSISINT(row->maxactivity - row->constant, SCIP_DEFAULT_SUMEPSILON));
}

/** returns the minimal activity of a row w.r.t. the columns' bounds */
SCIP_Real SCIProwGetMinActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(row->validactivitybdsdomchg <= stat->domchgcount);

   /* check, if activity bounds has to be calculated */
   if( row->validactivitybdsdomchg != stat->domchgcount )
      rowCalcActivityBounds(row, set, stat);
   assert(row->validactivitybdsdomchg == stat->domchgcount);
   assert(row->minactivity < SCIP_INVALID);
   assert(row->maxactivity < SCIP_INVALID);

   return row->minactivity;
}

/** returns the maximal activity of a row w.r.t. the columns' bounds */
SCIP_Real SCIProwGetMaxActivity(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(row != NULL);
   assert(stat != NULL);
   assert(row->validactivitybdsdomchg <= stat->domchgcount);

   /* check, if activity bounds has to be calculated */
   if( row->validactivitybdsdomchg != stat->domchgcount )
      rowCalcActivityBounds(row, set, stat);
   assert(row->validactivitybdsdomchg == stat->domchgcount);
   assert(row->minactivity < SCIP_INVALID);
   assert(row->maxactivity < SCIP_INVALID);

   return row->maxactivity;
}

/** returns whether the row is unmodifiable and redundant w.r.t. the columns' bounds */
SCIP_Bool SCIProwIsRedundant(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(row != NULL);

   if( row->modifiable )
      return FALSE;
   if( !SCIPsetIsInfinity(set, -row->lhs) )
   {
      SCIP_Real minactivity;

      minactivity = SCIProwGetMinActivity(row, set, stat);
      if( SCIPsetIsFeasLT(set, minactivity, row->lhs) )
         return FALSE;
   }
   if( !SCIPsetIsInfinity(set, row->rhs) )
   {
      SCIP_Real maxactivity;

      maxactivity = SCIProwGetMaxActivity(row, set, stat);
      if( SCIPsetIsFeasGT(set, maxactivity, row->rhs) )
         return FALSE;
   }

   return TRUE;
}

/** gets maximal absolute value of row vector coefficients */
SCIP_Real SCIProwGetMaxval(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);

   if( row->nummaxval == 0 )
      rowCalcIdxsAndVals(row, set);
   assert(row->nummaxval > 0);
   assert(row->maxval >= 0.0 || row->len == 0);

   return row->maxval;
}

/** gets minimal absolute value of row vector's non-zero coefficients */
SCIP_Real SCIProwGetMinval(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);

   if( row->numminval == 0 )
      rowCalcIdxsAndVals(row, set);
   assert(row->numminval > 0);
   assert(row->minval >= 0.0 || row->len == 0);

   return row->minval;
}

/** gets maximal column index of row entries */
int SCIProwGetMaxidx(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);

   if( row->validminmaxidx == 0 )
      rowCalcIdxsAndVals(row, set);
   assert(row->maxidx >= 0 || row->len == 0);
   assert(row->validminmaxidx);

   return row->maxidx;
}

/** gets minimal column index of row entries */
int SCIProwGetMinidx(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);

   if( row->validminmaxidx == 0 )
      rowCalcIdxsAndVals(row, set);
   assert(row->minidx >= 0 || row->len == 0);
   assert(row->validminmaxidx);

   return row->minidx;
}

/** gets number of integral columns in row */
int SCIProwGetNumIntCols(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(row != NULL);

   if( row->numintcols == -1 )
      rowCalcIdxsAndVals(row, set);

   assert(row->numintcols <= row->len && row->numintcols >= 0);

   return row->numintcols;
}

/** returns row's efficacy with respect to the current LP solution: e = -feasibility/norm */
SCIP_Real SCIProwGetLPEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real norm;
   SCIP_Real feasibility;
   SCIP_Real eps;

   assert(set != NULL);

   switch( set->sepa_efficacynorm )
   {
   case 'e':
      norm = SCIProwGetNorm(row);
      break;
   case 'm':
      norm = SCIProwGetMaxval(row, set);
      break;
   case 's':
      norm = SCIProwGetSumNorm(row);
      break;
   case 'd':
      norm = (row->len == 0 ? 0.0 : 1.0);
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", set->sepa_efficacynorm);
      SCIPABORT();
      norm = 0.0; /*lint !e527*/
   }

   eps = SCIPsetSumepsilon(set);
   norm = MAX(norm, eps);
   feasibility = SCIProwGetLPFeasibility(row, set, stat, lp);

   return -feasibility / norm;
}

/** returns whether the row's efficacy with respect to the current LP solution is greater than the minimal cut efficacy */
SCIP_Bool SCIProwIsLPEfficacious(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             root                /**< should the root's minimal cut efficacy be used? */
   )
{
   SCIP_Real efficacy;

   efficacy = SCIProwGetLPEfficacy(row, set, stat, lp);

   return SCIPsetIsEfficacious(set, root, efficacy);
}

/** returns row's efficacy with respect to the given primal solution: e = -feasibility/norm */
SCIP_Real SCIProwGetSolEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_Real norm;
   SCIP_Real feasibility;
   SCIP_Real eps;

   assert(set != NULL);

   switch( set->sepa_efficacynorm )
   {
   case 'e':
      norm = SCIProwGetNorm(row);
      break;
   case 'm':
      norm = SCIProwGetMaxval(row, set);
      break;
   case 's':
      norm = SCIProwGetSumNorm(row);
      break;
   case 'd':
      norm = (row->len == 0 ? 0.0 : 1.0);
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", set->sepa_efficacynorm);
      SCIPABORT();
      norm = 0.0; /*lint !e527*/
   }

   eps = SCIPsetSumepsilon(set);
   norm = MAX(norm, eps);
   feasibility = SCIProwGetSolFeasibility(row, set, stat, sol);

   return -feasibility / norm;
}

/** returns whether the row's efficacy with respect to the given primal solution is greater than the minimal cut
 *  efficacy
 */
SCIP_Bool SCIProwIsSolEfficacious(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             root                /**< should the root's minimal cut efficacy be used? */
   )
{
   SCIP_Real efficacy;

   efficacy = SCIProwGetSolEfficacy(row, set, stat, sol);

   return SCIPsetIsEfficacious(set, root, efficacy);
}

/** returns row's efficacy with respect to the relaxed solution: e = -feasibility/norm */
SCIP_Real SCIProwGetRelaxEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   SCIP_Real norm;
   SCIP_Real feasibility;
   SCIP_Real eps;

   assert(set != NULL);

   switch( set->sepa_efficacynorm )
   {
   case 'e':
      norm = SCIProwGetNorm(row);
      break;
   case 'm':
      norm = SCIProwGetMaxval(row, set);
      break;
   case 's':
      norm = SCIProwGetSumNorm(row);
      break;
   case 'd':
      norm = (row->len == 0 ? 0.0 : 1.0);
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", set->sepa_efficacynorm);
      SCIPABORT();
      norm = 0.0; /*lint !e527*/
   }

   eps = SCIPsetSumepsilon(set);
   norm = MAX(norm, eps);
   feasibility = SCIProwGetRelaxFeasibility(row, set, stat);

   return -feasibility / norm;
}

/** returns row's efficacy with respect to the NLP solution: e = -feasibility/norm */
SCIP_Real SCIProwGetNLPEfficacy(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   SCIP_Real norm;
   SCIP_Real feasibility;
   SCIP_Real eps;

   assert(set != NULL);

   switch( set->sepa_efficacynorm )
   {
   case 'e':
      norm = SCIProwGetNorm(row);
      break;
   case 'm':
      norm = SCIProwGetMaxval(row, set);
      break;
   case 's':
      norm = SCIProwGetSumNorm(row);
      break;
   case 'd':
      norm = (row->len == 0 ? 0.0 : 1.0);
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", set->sepa_efficacynorm);
      SCIPABORT();
      norm = 0.0; /*lint !e527*/
   }

   eps = SCIPsetSumepsilon(set);
   norm = MAX(norm, eps);
   feasibility = SCIProwGetNLPFeasibility(row, set, stat);

   return -feasibility / norm;
}

/** returns the scalar product of the coefficient vectors of the two given rows
 *
 *  @note the scalar product is computed w.r.t. the current LP columns only
 *  @todo also consider non-LP columns for the computation?
 */
SCIP_Real SCIProwGetScalarProduct(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2                /**< second LP row */
   )
{
   SCIP_Real scalarprod;
   int* row1colsidx;
   int* row2colsidx;
   int i1;
   int i2;

   assert(row1 != NULL);
   assert(row2 != NULL);

   /* Sort the column indices of both rows.
    *
    * The columns in a row are divided into two parts: LP columns, which are currently in the LP and non-LP columns;
    * we sort the rows, but that only ensures that within these two parts, columns are sorted w.r.t. their index.
    * Normally, this should be suficient, because a column contained in both rows should either be one of the LP columns
    * for both or one of the non-LP columns for both.
    * However, directly after a row was created, before a row is added to the LP, the row is not linked to all its
    * columns and all columns are treated as non-LP columns. Moreover, for example when doing column generation,
    * columns can be added later and remain unlinked while all previously added columns might already be linked.
    * Therefore, we have to be very careful about whether we can rely on the partitioning of the variables.
    *
    * We distinguish the following cases:
    *
    * 1) both rows have no unlinked columns
    *    -> we just check the LP partitions
    *
    * 2) exactly one row is completely unlinked, the other one is completely linked
    *    -> we compare the non-LP (unlinked) partition with the LP partition of the other row
    *       (thus all common LP columns are regarded)
    *
    * 3) we have unlinked and LP columns in both rows
    *    -> we need to compare four partitions at once
    *
    * 4a) we have one row with unlinked and LP columns and the other without any unlinked columns
    *     -> we need to compare three partitions: the LP part of the completely linked row and both partitions of the
    *        other row
    *
    * 4b) we have one row with unlinked and LP columns and the other is completely unlinked
    *     -> we need to compare three partitions: the complete unlinked row and both partitions of the other row
    *
    * 5) both rows are completely unlinked
    *    -> we need to compare two partitions: both complete rows
    */
   SCIProwSort(row1);
   assert(row1->lpcolssorted);
   assert(row1->nonlpcolssorted);
   SCIProwSort(row2);
   assert(row2->lpcolssorted);
   assert(row2->nonlpcolssorted);

   assert(row1->nunlinked <= row1->len - row1->nlpcols);
   assert(row2->nunlinked <= row2->len - row2->nlpcols);

   row1colsidx = row1->cols_index;
   row2colsidx = row2->cols_index;

#ifndef NDEBUG
   /* check that we can rely on the partition into LP columns and non-LP columns if the rows are completely linked */
   if( row1->nunlinked == 0 && row2->nunlinked == 0 )
   {
      i1 = 0;
      i2 = row2->nlpcols;
      while( i1 < row1->nlpcols && i2 < row2->len )
      {
         assert(row1->cols[i1] != row2->cols[i2]);
         if( row1->cols[i1]->index < row2->cols[i2]->index )
            ++i1;
         else
         {
            assert(row1->cols[i1]->index > row2->cols[i2]->index);
            ++i2;
         }
      }
      assert(i1 == row1->nlpcols || i2 == row2->len);

      i1 = row1->nlpcols;
      i2 = 0;
      while( i1 < row1->len && i2 < row2->nlpcols )
      {
         assert(row1->cols[i1] != row2->cols[i2]);
         if( row1->cols[i1]->index < row2->cols[i2]->index )
            ++i1;
         else
         {
            assert(row1->cols[i1]->index > row2->cols[i2]->index);
            ++i2;
         }
      }
      assert(i1 == row1->len || i2 == row2->nlpcols);
   }
#endif

   /* The "easy" cases 1) and 2) */
   if( (row1->nunlinked == 0 && row2->nunlinked == 0) ||
      ((row1->nlpcols == row1->len || row1->nunlinked == row1->len)
         && (row2->nlpcols == row2->len || row2->nunlinked == row2->len)
         && (row1->nunlinked == 0 || row2->nunlinked == 0)) )
   {
      assert(row1->nunlinked == 0 || row1->nunlinked == row1->len);
      assert(row2->nunlinked == 0 || row2->nunlinked == row2->len);

      /* set the iterators to the last column we want to regard in the row: nunlinked is either 0 or row->len,
       * therefore, we get nlpcols if nunlinked is 0 and row->len if the row is completely unlinked
       */
      i1 = MAX(row1->nlpcols, row1->nunlinked) - 1;
      i2 = MAX(row2->nlpcols, row2->nunlinked) - 1;
      scalarprod = 0.0;

      /* calculate the scalar product */
      while( i1 >= 0 && i2 >= 0 )
      {
         assert(row1->cols[i1]->index == row1colsidx[i1]);
         assert(row2->cols[i2]->index == row2colsidx[i2]);
         assert((row1->cols[i1] == row2->cols[i2]) == (row1colsidx[i1] == row2colsidx[i2]));
         if( row1colsidx[i1] < row2colsidx[i2] )
            --i2;
         else if( row1colsidx[i1] > row2colsidx[i2] )
            --i1;
         else
         {
            scalarprod += row1->vals[i1] * row2->vals[i2];
            --i1;
            --i2;
         }
      }
   }
   /* the "harder" cases 3) - 5): start with four partitions and reduce their number iteratively */
   else
   {
      SCIP_Bool lpcols;
      int ilp1;
      int inlp1;
      int ilp2;
      int inlp2;
      int end1;
      int end2;

      scalarprod = 0;
      ilp1 = 0;
      ilp2 = 0;

      /* if a row is completely linked (case 4a), we do not have to consider its non-LP columns */
      inlp1 = (row1->nunlinked > 0 ? row1->nlpcols : row1->len);
      inlp2 = (row2->nunlinked > 0 ? row2->nlpcols : row2->len);

      /* handle the case of four partitions (case 3) until one partition is finished;
       * cases 4a), 4b), and 5) will fail the while-condition
       */
      while( ilp1 < row1->nlpcols && inlp1 < row1->len && ilp2 < row2->nlpcols && inlp2 < row2->len )
      {
         assert(row1->cols[ilp1]->index == row1colsidx[ilp1]);
         assert(row1->cols[inlp1]->index == row1colsidx[inlp1]);
         assert(row2->cols[ilp2]->index == row2colsidx[ilp2]);
         assert(row2->cols[inlp2]->index == row2colsidx[inlp2]);
         assert((row1->cols[ilp1] == row2->cols[ilp2]) == (row1colsidx[ilp1] == row2colsidx[ilp2]));
         assert((row1->cols[ilp1] == row2->cols[inlp2]) == (row1colsidx[ilp1] == row2colsidx[inlp2]));
         assert((row1->cols[inlp1] == row2->cols[ilp2]) == (row1colsidx[inlp1] == row2colsidx[ilp2]));
         assert((row1->cols[inlp1] == row2->cols[inlp2]) == (row1colsidx[inlp1] == row2colsidx[inlp2]));

         /* rows have the same linked LP columns */
         if( row1colsidx[ilp1] == row2colsidx[ilp2] )
         {
            scalarprod += row1->vals[ilp1] * row2->vals[ilp2];
            ++ilp1;
            ++ilp2;
         }
         /* LP column of row1 is the same as unlinked column of row2 */
         else if( row1colsidx[ilp1] == row2colsidx[inlp2] )
         {
            scalarprod += row1->vals[ilp1] * row2->vals[inlp2];
            ++ilp1;
            ++inlp2;
         }
         /* unlinked column of row1 is the same as LP column of row2 */
         else if( row1colsidx[inlp1] == row2colsidx[ilp2] )
         {
            scalarprod += row1->vals[inlp1] * row2->vals[ilp2];
            ++inlp1;
            ++ilp2;
         }
         /* two unlinked LP columns are the same */
         else if( row1colsidx[inlp1] == row2colsidx[inlp2] && row1->cols[inlp1]->lppos >= 0 )
         {
            scalarprod += row1->vals[inlp1] * row2->vals[inlp2];
            ++inlp1;
            ++inlp2;
         }
         /* increase smallest counter */
         else if( row1colsidx[ilp1] < row1colsidx[inlp1] )
         {
            if( row2colsidx[ilp2] < row2colsidx[inlp2] )
            {
               if( row1colsidx[ilp1] < row2colsidx[ilp2] )
                  ++ilp1;
               else
                  ++ilp2;
            }
            else
            {
               if( row1colsidx[ilp1] < row2colsidx[inlp2] )
                  ++ilp1;
               else
                  ++inlp2;
            }
         }
         else
         {
            if( row2colsidx[ilp2] < row2colsidx[inlp2] )
            {
               if( row1colsidx[inlp1] < row2colsidx[ilp2] )
                  ++inlp1;
               else
                  ++ilp2;
            }
            else
            {
               if( row1colsidx[inlp1] < row2colsidx[inlp2] )
                  ++inlp1;
               else
                  ++inlp2;
            }
         }
      }

      /* One partition was completely handled, we just have to handle the three remaining partitions:
       * the remaining partition of this row and the two partitions of the other row.
       * If necessary, we swap the partitions to ensure that row1 is the row with only one remaining partition.
       */
      if( ilp1 != row1->nlpcols && inlp1 != row1->len )
      {
         int tmpilp;
         int tmpinlp;

         assert(ilp2 == row2->nlpcols || inlp2 == row2->len);

         SCIPswapPointers((void**) &row1, (void**) &row2);
         SCIPswapPointers((void**) &row1colsidx, (void**) &row2colsidx);
         tmpilp = ilp1;
         tmpinlp = inlp1;
         ilp1 = ilp2;
         inlp1 = inlp2;
         ilp2 = tmpilp;
         inlp2 = tmpinlp;
      }

      /* determine section of row 1 that we want to look at (current iterator = begin, end, LP-columns?)
       * -> this merges cases 4a) and 4b)
       */
      if( ilp1 == row1->nlpcols )
      {
         i1 = inlp1;
         end1 = row1->len;
         lpcols = FALSE;
      }
      else
      {
         assert(inlp1 == row1->len);

         i1 = ilp1;
         end1 = row1->nlpcols;
         lpcols = TRUE;
      }

      /* handle the case of three partitions (case 4) until one partition is finished, this reduces our problem to case 1), 2), or 5);
       * case 5) will fail the while-condition
       */
      while( i1 < end1 && ilp2 < row2->nlpcols && inlp2 < row2->len )
      {
         assert(row1->cols[i1]->index == row1colsidx[i1]);
         assert(row2->cols[ilp2]->index == row2colsidx[ilp2]);
         assert(row2->cols[inlp2]->index == row2colsidx[inlp2]);
         assert((row1->cols[i1] == row2->cols[ilp2]) == (row1colsidx[i1] == row2colsidx[ilp2]));
         assert((row1->cols[i1] == row2->cols[inlp2]) == (row1colsidx[i1] == row2colsidx[inlp2]));

         /* current column in row 1 is the same as the current LP column in row 2 */
         if( row1colsidx[i1] == row2colsidx[ilp2] )
         {
            scalarprod += row1->vals[i1] * row2->vals[ilp2];
            ++i1;
            ++ilp2;
         }
         /* linked or unlinked LP column of row1 is the same as unlinked column of row2 */
         else if( row1colsidx[i1] == row2colsidx[inlp2] && (lpcols || row1->cols[i1]->lppos >= 0) )
         {
            scalarprod += row1->vals[i1] * row2->vals[inlp2];
            ++i1;
            ++inlp2;
         }
         /* increase smallest counter */
         else if( row2colsidx[ilp2] < row2colsidx[inlp2] )
         {
            if( row1colsidx[i1] < row2colsidx[ilp2] )
               ++i1;
            else
               ++ilp2;
         }
         else
         {
            if( row1colsidx[i1] < row2colsidx[inlp2] )
               ++i1;
            else
               ++inlp2;
         }
      }

      /* if the second section of row 1 was finished, we can stop; otherwise, we have to consider the remaining parts of
       * the two rows
       */
      if( i1 < end1 )
      {
         /* determine section of row 2 that we want to look at (current iterator = begin, end, LP-columns?) */
         if( ilp2 == row2->nlpcols )
         {
            i2 = inlp2;
            end2 = row2->len;
            lpcols = FALSE;
         }
         else
         {
            assert(inlp2 == row2->len);

            i2 = ilp2;
            end2 = row2->nlpcols;
         }

         /* handle the case of two partitions (standard case 5, or case 1 or 2 due to partition reduction) */
         while( i1 < end1 && i2 < end2 )
         {
            assert(row1->cols[i1]->index == row1colsidx[i1]);
            assert(row2->cols[i2]->index == row2colsidx[i2]);
            assert((row1->cols[i1] == row2->cols[i2]) == (row1colsidx[i1] == row2colsidx[i2]));

            /* linked or unlinked LP column of row1 is the same as linked or unlinked LP column of row2 */
            if( row1colsidx[i1] == row2colsidx[i2] && (lpcols || row1->cols[i1]->lppos >= 0) )
            {
               scalarprod += row1->vals[i1] * row2->vals[i2];
               ++i1;
               ++i2;
            }
            /* increase smallest counter */
            else if( row1colsidx[i1] < row2colsidx[i2] )
               ++i1;
            else
               ++i2;
         }
      }
   }

   return scalarprod;
}

/** returns the discrete scalar product of the coefficient vectors of the two given rows */
static
int SCIProwGetDiscreteScalarProduct(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2                /**< second LP row */
   )
{
   int prod;
   int* row1colsidx;
   int* row2colsidx;
   int i1;
   int i2;

   assert(row1 != NULL);
   assert(row2 != NULL);

   /* Sort the column indices of both rows.
    *
    * The columns in a row are divided into two parts: LP columns, which are currently in the LP and non-LP columns;
    * we sort the rows, but that only ensures that within these two parts, columns are sorted w.r.t. their index.
    * Normally, this should be suficient, because a column contained in both rows should either be one of the LP columns
    * for both or one of the non-LP columns for both.
    * However, directly after a row was created, before a row is added to the LP, the row is not linked to all its
    * columns and all columns are treated as non-LP columns. Moreover, for example when doing column generation,
    * columns can be added later and remain unlinked while all previously added columns might already be linked.
    * Therefore, we have to be very careful about whether we can rely on the partitioning of the variables.
    *
    * We distinguish the following cases:
    *
    * 1) both rows have no unlinked columns
    *    -> we just check the LP partitions
    *
    * 2) exactly one row is completely unlinked, the other one is completely linked
    *    -> we compare the non-LP (unlinked) partition with the LP partition of the other row
    *       (thus all common LP columns are regarded)
    *
    * 3) we have unlinked and LP columns in both rows
    *    -> we need to compare four partitions at once
    *
    * 4a) we have one row with unlinked and LP columns and the other without any unlinked columns
    *     -> we need to compare three partitions: the LP part of the completely linked row and both partitions of the
    *        other row
    *
    * 4b) we have one row with unlinked and LP columns and the other is completely unlinked
    *     -> we need to compare three partitions: the complete unlinked row and both partitions of the other row
    *
    * 5) both rows are completely unlinked
    *    -> we need to compare two partitions: both complete rows
    */
   SCIProwSort(row1);
   assert(row1->lpcolssorted);
   assert(row1->nonlpcolssorted);
   SCIProwSort(row2);
   assert(row2->lpcolssorted);
   assert(row2->nonlpcolssorted);

   assert(row1->nunlinked <= row1->len - row1->nlpcols);
   assert(row2->nunlinked <= row2->len - row2->nlpcols);

   row1colsidx = row1->cols_index;
   row2colsidx = row2->cols_index;

#ifndef NDEBUG
   /* check that we can rely on the partition into LP columns and non-LP columns if the rows are completely linked */
   if( row1->nunlinked == 0 && row2->nunlinked == 0 )
   {
      i1 = 0;
      i2 = row2->nlpcols;
      while( i1 < row1->nlpcols && i2 < row2->len )
      {
         assert(row1->cols[i1] != row2->cols[i2]);
         if( row1->cols[i1]->index < row2->cols[i2]->index )
            ++i1;
         else
         {
            assert(row1->cols[i1]->index > row2->cols[i2]->index);
            ++i2;
         }
      }
      assert(i1 == row1->nlpcols || i2 == row2->len);

      i1 = row1->nlpcols;
      i2 = 0;
      while( i1 < row1->len && i2 < row2->nlpcols )
      {
         assert(row1->cols[i1] != row2->cols[i2]);
         if( row1->cols[i1]->index < row2->cols[i2]->index )
            ++i1;
         else
         {
            assert(row1->cols[i1]->index > row2->cols[i2]->index);
            ++i2;
         }
      }
      assert(i1 == row1->len || i2 == row2->nlpcols);
   }
#endif

   /* The "easy" cases 1) and 2) */
   if( (row1->nunlinked == 0 && row2->nunlinked == 0) ||
      ((row1->nlpcols == row1->len || row1->nunlinked == row1->len)
         && (row2->nlpcols == row2->len || row2->nunlinked == row2->len)
         && (row1->nunlinked == 0 || row2->nunlinked == 0)) )
   {
      assert(row1->nunlinked == 0 || row1->nunlinked == row1->len);
      assert(row2->nunlinked == 0 || row2->nunlinked == row2->len);

      /* set the iterators to the last column we want to regard in the row: nunlinked is either 0 or row->len,
       * therefore, we get nlpcols if nunlinked is 0 and row->len if the row is completely unlinked
       */
      i1 = MAX(row1->nlpcols, row1->nunlinked) - 1;
      i2 = MAX(row2->nlpcols, row2->nunlinked) - 1;
      prod = 0;

      /* calculate the scalar product */
      while( i1 >= 0 && i2 >= 0 )
      {
         assert(row1->cols[i1]->index == row1colsidx[i1]);
         assert(row2->cols[i2]->index == row2colsidx[i2]);
         assert((row1->cols[i1] == row2->cols[i2]) == (row1colsidx[i1] == row2colsidx[i2]));
         if( row1colsidx[i1] < row2colsidx[i2] )
            --i2;
         else if( row1colsidx[i1] > row2colsidx[i2] )
            --i1;
         else
         {
            ++prod;
            --i1;
            --i2;
         }
      }
   }
   /* the "harder" cases 3) - 5): start with four partitions and reduce their number iteratively */
   else
   {
      SCIP_Bool lpcols;
      int ilp1;
      int inlp1;
      int ilp2;
      int inlp2;
      int end1;
      int end2;

      prod = 0;
      ilp1 = 0;
      ilp2 = 0;

      /* if a row is completely linked (case 4a), we do not have to consider its non-LP columns */
      inlp1 = (row1->nunlinked > 0 ? row1->nlpcols : row1->len);
      inlp2 = (row2->nunlinked > 0 ? row2->nlpcols : row2->len);

      /* handle the case of four partitions (case 3) until one partition is finished;
       * cases 4a), 4b), and 5) will fail the while-condition
       */
      while( ilp1 < row1->nlpcols && inlp1 < row1->len && ilp2 < row2->nlpcols && inlp2 < row2->len )
      {
         assert(row1->cols[ilp1]->index == row1colsidx[ilp1]);
         assert(row1->cols[inlp1]->index == row1colsidx[inlp1]);
         assert(row2->cols[ilp2]->index == row2colsidx[ilp2]);
         assert(row2->cols[inlp2]->index == row2colsidx[inlp2]);
         assert((row1->cols[ilp1] == row2->cols[ilp2]) == (row1colsidx[ilp1] == row2colsidx[ilp2]));
         assert((row1->cols[ilp1] == row2->cols[inlp2]) == (row1colsidx[ilp1] == row2colsidx[inlp2]));
         assert((row1->cols[inlp1] == row2->cols[ilp2]) == (row1colsidx[inlp1] == row2colsidx[ilp2]));
         assert((row1->cols[inlp1] == row2->cols[inlp2]) == (row1colsidx[inlp1] == row2colsidx[inlp2]));

         /* rows have the same linked LP columns */
         if( row1colsidx[ilp1] == row2colsidx[ilp2] )
         {
            ++prod;
            ++ilp1;
            ++ilp2;
         }
         /* LP column of row1 is the same as unlinked column of row2 */
         else if( row1colsidx[ilp1] == row2colsidx[inlp2] )
         {
            ++prod;
            ++ilp1;
            ++inlp2;
         }
         /* unlinked column of row1 is the same as LP column of row2 */
         else if( row1colsidx[inlp1] == row2colsidx[ilp2] )
         {
            ++prod;
            ++inlp1;
            ++ilp2;
         }
         /* two unlinked LP columns are the same */
         else if( row1colsidx[inlp1] == row2colsidx[inlp2] && row1->cols[inlp1]->lppos >= 0 )
         {
            ++prod;
            ++inlp1;
            ++inlp2;
         }
         /* increase smallest counter */
         else if( row1colsidx[ilp1] < row1colsidx[inlp1] )
         {
            if( row2colsidx[ilp2] < row2colsidx[inlp2] )
            {
               if( row1colsidx[ilp1] < row2colsidx[ilp2] )
                  ++ilp1;
               else
                  ++ilp2;
            }
            else
            {
               if( row1colsidx[ilp1] < row2colsidx[inlp2] )
                  ++ilp1;
               else
                  ++inlp2;
            }
         }
         else
         {
            if( row2colsidx[ilp2] < row2colsidx[inlp2] )
            {
               if( row1colsidx[inlp1] < row2colsidx[ilp2] )
                  ++inlp1;
               else
                  ++ilp2;
            }
            else
            {
               if( row1colsidx[inlp1] < row2colsidx[inlp2] )
                  ++inlp1;
               else
                  ++inlp2;
            }
         }
      }

      /* One partition was completely handled, we just have to handle the three remaining partitions:
       * the remaining partition of this row and the two partitions of the other row.
       * If necessary, we swap the partitions to ensure that row1 is the row with only one remaining partition.
       */
      if( ilp1 != row1->nlpcols && inlp1 != row1->len )
      {
         int tmpilp;
         int tmpinlp;

         assert(ilp2 == row2->nlpcols || inlp2 == row2->len);

         SCIPswapPointers((void**) &row1, (void**) &row2);
         SCIPswapPointers((void**) &row1colsidx, (void**) &row2colsidx);
         tmpilp = ilp1;
         tmpinlp = inlp1;
         ilp1 = ilp2;
         inlp1 = inlp2;
         ilp2 = tmpilp;
         inlp2 = tmpinlp;
      }

      /* determine section of row 1 that we want to look at (current iterator = begin, end, LP-columns?)
       * -> this merges cases 4a) and 4b)
       */
      if( ilp1 == row1->nlpcols )
      {
         i1 = inlp1;
         end1 = row1->len;
         lpcols = FALSE;
      }
      else
      {
         assert(inlp1 == row1->len);

         i1 = ilp1;
         end1 = row1->nlpcols;
         lpcols = TRUE;
      }

      /* handle the case of three partitions (case 4) until one partition is finished, this reduces our problem to case 1), 2), or 5);
       * case 5) will fail the while-condition
       */
      while( i1 < end1 && ilp2 < row2->nlpcols && inlp2 < row2->len )
      {
         assert(row1->cols[i1]->index == row1colsidx[i1]);
         assert(row2->cols[ilp2]->index == row2colsidx[ilp2]);
         assert(row2->cols[inlp2]->index == row2colsidx[inlp2]);
         assert((row1->cols[i1] == row2->cols[ilp2]) == (row1colsidx[i1] == row2colsidx[ilp2]));
         assert((row1->cols[i1] == row2->cols[inlp2]) == (row1colsidx[i1] == row2colsidx[inlp2]));

         /* current column in row 1 is the same as the current LP column in row 2 */
         if( row1colsidx[i1] == row2colsidx[ilp2] )
         {
            ++prod;
            ++i1;
            ++ilp2;
         }
         /* linked or unlinked LP column of row1 is the same as unlinked column of row2 */
         else if( row1colsidx[i1] == row2colsidx[inlp2] && (lpcols || row1->cols[i1]->lppos >= 0) )
         {
            ++prod;
            ++i1;
            ++inlp2;
         }
         /* increase smallest counter */
         else if( row2colsidx[ilp2] < row2colsidx[inlp2] )
         {
            if( row1colsidx[i1] < row2colsidx[ilp2] )
               ++i1;
            else
               ++ilp2;
         }
         else
         {
            if( row1colsidx[i1] < row2colsidx[inlp2] )
               ++i1;
            else
               ++inlp2;
         }
      }

      /* if the second section of row 1 was finished, we can stop; otherwise, we have to consider the remaining parts of
       * the two rows
       */
      if( i1 < end1 )
      {
         /* determine section of row 2 that we want to look at (current iterator = begin, end, LP-columns?) */
         if( ilp2 == row2->nlpcols )
         {
            i2 = inlp2;
            end2 = row2->len;
            lpcols = FALSE;
         }
         else
         {
            assert(inlp2 == row2->len);

            i2 = ilp2;
            end2 = row2->nlpcols;
         }

         /* handle the case of two partitions (standard case 5, or case 1 or 2 due to partition reduction) */
         while( i1 < end1 && i2 < end2 )
         {
            assert(row1->cols[i1]->index == row1colsidx[i1]);
            assert(row2->cols[i2]->index == row2colsidx[i2]);
            assert((row1->cols[i1] == row2->cols[i2]) == (row1colsidx[i1] == row2colsidx[i2]));

            /* linked or unlinked LP column of row1 is the same as linked or unlinked LP column of row2 */
            if( row1colsidx[i1] == row2colsidx[i2] && (lpcols || row1->cols[i1]->lppos >= 0) )
            {
               ++prod;
               ++i1;
               ++i2;
            }
            /* increase smallest counter */
            else if( row1colsidx[i1] < row2colsidx[i2] )
               ++i1;
            else
               ++i2;
         }
      }
   }

   return prod;
}

/** returns the degree of parallelism between the hyperplanes defined by the two row vectors v, w:
 *  p = |v*w|/(|v|*|w|);
 *  the hyperplanes are parallel, iff p = 1, they are orthogonal, iff p = 0
 */
SCIP_Real SCIProwGetParallelism(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2,               /**< second LP row */
   char                  orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   )
{
   SCIP_Real parallelism;
   SCIP_Real scalarprod;

   switch( orthofunc )
   {
   case 'e':
      scalarprod = SCIProwGetScalarProduct(row1, row2);
      if( scalarprod == 0.0 )
      {
         parallelism = 0.0;
         break;
      }

      if( SCIProwGetNorm(row1) == 0.0 )
      {
         /* In theory, this should not happen if the scalarproduct is not zero
          * But due to bug 520 (also issue 44), it is possible that norms are not correct.
          * Thus, if the norm is so bad that it is even 0, then reevaluate it here.
          * But as we don't have set available here, we cannot call rowCalcNorms, so do it by hand.
          */
         int i;
         for( i = 0; i < row1->len; ++i )
            if( row1->cols[i]->lppos >= 0 )
               row1->sqrnorm += SQR(row1->vals[i]);
         assert(SCIProwGetNorm(row1) != 0.0);
      }

      if( SCIProwGetNorm(row2) == 0.0 )
      {
         /* same as for row1 above: reeval norms if it is 0, which is wrong */
         int i;
         for( i = 0; i < row2->len; ++i )
            if( row2->cols[i]->lppos >= 0 )
               row2->sqrnorm += SQR(row2->vals[i]);
         assert(SCIProwGetNorm(row2) != 0.0);
      }

      parallelism = REALABS(scalarprod) / (SCIProwGetNorm(row1) * SCIProwGetNorm(row2));
      break;

   case 'd':
      scalarprod = (SCIP_Real) SCIProwGetDiscreteScalarProduct(row1, row2);
      parallelism = scalarprod / (sqrt((SCIP_Real) SCIProwGetNNonz(row1)) * sqrt((SCIP_Real) SCIProwGetNNonz(row2)));
      break;

   default:
      SCIPerrorMessage("invalid orthogonality function parameter '%c'\n", orthofunc);
      SCIPABORT();
      parallelism = 0.0; /*lint !e527*/
   }

   return parallelism;
}

/** returns the degree of orthogonality between the hyperplanes defined by the two row vectors v, w:
 *  o = 1 - |v*w|/(|v|*|w|);
 *  the hyperplanes are orthogonal, iff p = 1, they are parallel, iff p = 0
 */
SCIP_Real SCIProwGetOrthogonality(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2,               /**< second LP row */
   char                  orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   )
{
   return 1.0 - SCIProwGetParallelism(row1, row2, orthofunc);
}

/** gets parallelism of row with objective function: if the returned value is 1, the row is parallel to the objective
 *  function, if the value is 0, it is orthogonal to the objective function
 */
SCIP_Real SCIProwGetObjParallelism(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real prod;
   SCIP_Real parallelism;

   assert(row != NULL);
   assert(lp != NULL);

   if( lp->objsqrnormunreliable )
      SCIPlpRecalculateObjSqrNorm(set, lp);

   assert(!lp->objsqrnormunreliable);
   assert(lp->objsqrnorm >= 0.0);

   checkRowSqrnorm(row);
   checkRowObjprod(row);

   prod = row->sqrnorm * lp->objsqrnorm;

   parallelism = SCIPsetIsPositive(set, prod) ? REALABS(row->objprod) / SQRT(prod) : 0.0;
   assert(SCIPsetIsSumGE(set, parallelism, 0.0));
   assert(SCIPsetIsSumLE(set, parallelism, 1.0));
   parallelism = MIN(parallelism, 1.0);
   parallelism = MAX(parallelism, 0.0);

   return parallelism;
}

/** includes event handler with given data in row's event filter */
SCIP_RETCODE SCIProwCatchEvent(
   SCIP_ROW*             row,                /**< row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   assert(row != NULL);
   assert(row->eventfilter != NULL);
   assert((eventtype & ~SCIP_EVENTTYPE_ROWCHANGED) == 0);
   assert((eventtype &  SCIP_EVENTTYPE_ROWCHANGED) != 0);

   SCIPsetDebugMsg(set, "catch event of type 0x%" SCIP_EVENTTYPE_FORMAT " of row <%s> with handler %p and data %p\n",
      eventtype, row->name, (void*)eventhdlr, (void*)eventdata);

   SCIP_CALL( SCIPeventfilterAdd(row->eventfilter, blkmem, set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** deletes event handler with given data from row's event filter */
SCIP_RETCODE SCIProwDropEvent(
   SCIP_ROW*             row,                /**< row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int                   filterpos           /**< position of event filter entry returned by SCIPvarCatchEvent(), or -1 */
   )
{
   assert(row != NULL);
   assert(row->eventfilter != NULL);

   SCIPsetDebugMsg(set, "drop event of row <%s> with handler %p and data %p\n", row->name, (void*)eventhdlr, (void*)eventdata);

   SCIP_CALL( SCIPeventfilterDel(row->eventfilter, blkmem, set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** marks a row to be not removable from the LP in the current node because it became obsolete */
void SCIProwMarkNotRemovableLocal(
   SCIP_ROW*             row,                /**< LP row */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(row  != NULL);
   assert(stat != NULL);
   assert(stat->nnodes > 0);

   /* lpRemoveObsoleteRows() does not remove a row if the node number stored in obsoletenode equals the current node number */
   row->obsoletenode = stat->nnodes;
}

/*
 * LP solver data update
 */

/** resets column data to represent a column not in the LP solver */
static
void markColDeleted(
   SCIP_COL*             col                 /**< column to be marked deleted */
   )
{
   assert(col != NULL);

   col->lpipos = -1;
   col->primsol = 0.0;
   col->redcost = SCIP_INVALID;
   col->farkascoef = SCIP_INVALID;
   col->sbdown = SCIP_INVALID;
   col->sbup = SCIP_INVALID;
   col->sbdownvalid = FALSE;
   col->sbupvalid = FALSE;
   col->validredcostlp = -1;
   col->validfarkaslp = -1;
   col->sbitlim = -1;
   col->basisstatus = SCIP_BASESTAT_ZERO; /*lint !e641*/
}

/** applies all cached column removals to the LP solver */
static
SCIP_RETCODE lpFlushDelCols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->lpifirstchgcol <= lp->nlpicols);
   assert(lp->lpifirstchgcol <= lp->ncols);

   /* find the first column to change */
   while( lp->lpifirstchgcol < lp->nlpicols
      && lp->lpifirstchgcol < lp->ncols
      && lp->cols[lp->lpifirstchgcol]->lpipos == lp->lpifirstchgcol
      && !lp->cols[lp->lpifirstchgcol]->coefchanged )
   {
      assert(lp->cols[lp->lpifirstchgcol] == lp->lpicols[lp->lpifirstchgcol]);
      lp->lpifirstchgcol++;
   }

   /* shrink LP to the part which didn't change */
   if( lp->lpifirstchgcol < lp->nlpicols )
   {
      int i;

      assert(!lp->diving);
      SCIPdebugMessage("flushing col deletions: shrink LP from %d to %d columns\n", lp->nlpicols, lp->lpifirstchgcol);
      SCIP_CALL( SCIPlpiDelCols(lp->lpi, lp->lpifirstchgcol, lp->nlpicols-1) );
      for( i = lp->lpifirstchgcol; i < lp->nlpicols; ++i )
      {
         markColDeleted(lp->lpicols[i]);
      }
      lp->nlpicols = lp->lpifirstchgcol;
      lp->flushdeletedcols = TRUE;
      lp->updateintegrality = TRUE;

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->primalchecked = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   assert(lp->nlpicols == lp->lpifirstchgcol);

   return SCIP_OKAY;
}

/** computes for the given column the lower and upper bound that should be flushed into the LP
 *  depending on lazy bounds and diving mode; in diving mode, lazy bounds are ignored, i.e.,
 *  the bounds are explicitly added to the LP in any case
 */
static
void computeLPBounds(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col,                /**< column to compute bounds for */
   SCIP_Real             lpiinf,             /**< infinity value if the LP solver */
   SCIP_Real*            lb,                 /**< pointer to store the new lower bound */
   SCIP_Real*            ub                  /**< pointer to store the new upper bound */
   )
{
   assert(lp != NULL);
   assert(set != NULL);
   assert(col != NULL);
   assert(lb != NULL);
   assert(ub != NULL);

   /* get the correct new lower bound:
    * if lazy lower bound exists and is larger than lower bound, set lower bound to infinity;
    * if we are in diving mode, ignore lazy bounds and always take the lower bound
    */
   if( SCIPsetIsInfinity(set, -col->lb) || (SCIPsetIsLE(set, col->lb, col->lazylb) && !SCIPlpDiving(lp)) )
      (*lb) = -lpiinf;
   else
      (*lb) = col->lb;
   /* get the correct new upper bound:
    * if lazy upper bound exists and is larger than upper bound, set upper bound to infinity;
    * if we are in diving mode, ignore lazy bounds and always take the upper bound
    */
   if( SCIPsetIsInfinity(set, col->ub) || (SCIPsetIsGE(set, col->ub, col->lazyub) && !SCIPlpDiving(lp)) )
      (*ub) = lpiinf;
   else
      (*ub) = col->ub;
}

/** applies all cached column additions to the LP solver */
static
SCIP_RETCODE lpFlushAddCols(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int* beg;
   int* ind;
   SCIP_Real* val;
   char** name;
   SCIP_COL* col;
   SCIP_Real lpiinf;
   int c;
   int pos;
   int nnonz;
   int naddcols;
   int naddcoefs;
   int i;
   int lpipos;

   assert(lp != NULL);
   assert(lp->lpifirstchgcol == lp->nlpicols);
   assert(blkmem != NULL);
   assert(set != NULL);

   /* if there are no columns to add, we are ready */
   if( lp->ncols == lp->nlpicols )
      return SCIP_OKAY;

   /* add the additional columns */
   assert(!lp->diving);
   assert(lp->ncols > lp->nlpicols);
   SCIP_CALL( ensureLpicolsSize(lp, set, lp->ncols) );

   /* get the solver's infinity value */
   lpiinf = SCIPlpiInfinity(lp->lpi);

   /* count the (maximal) number of added coefficients, calculate the number of added columns */
   naddcols = lp->ncols - lp->nlpicols;
   naddcoefs = 0;
   for( c = lp->nlpicols; c < lp->ncols; ++c )
      naddcoefs += lp->cols[c]->len;
   assert(naddcols > 0);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &obj, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lb, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ub, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &beg, naddcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &val, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, naddcols) );

   /* fill temporary memory with column data */
   nnonz = 0;
   for( pos = 0, c = lp->nlpicols; c < lp->ncols; ++pos, ++c )
   {
      col = lp->cols[c];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);
      assert(col->lppos == c);
      assert(nnonz + col->nlprows <= naddcoefs);

      SCIPsetDebugMsg(set, "flushing added column <%s>: ", SCIPvarGetName(col->var));
      debugColPrint(set, col);

      /* Because the column becomes a member of the LP solver, it now can take values
       * different from zero. That means, we have to include the column in the corresponding
       * row vectors.
       */
      SCIP_CALL( colLink(col, blkmem, set, eventqueue, lp) );

      lp->lpicols[c] = col;
      col->lpipos = c;
      col->primsol = SCIP_INVALID;
      col->redcost = SCIP_INVALID;
      col->farkascoef = SCIP_INVALID;
      col->sbdown = SCIP_INVALID;
      col->sbup = SCIP_INVALID;
      col->sbdownvalid = FALSE;
      col->sbupvalid = FALSE;
      col->validredcostlp = -1;
      col->validfarkaslp = -1;
      col->sbitlim = -1;
      col->objchanged = FALSE;
      col->lbchanged = FALSE;
      col->ubchanged = FALSE;
      col->coefchanged = FALSE;
      obj[pos] = col->obj;

      /* compute bounds that should be flushed into the LP (taking into account lazy bounds) */
      computeLPBounds(lp, set, col, lpiinf, &(lb[pos]), &(ub[pos]));

      beg[pos] = nnonz;
      name[pos] = (char*)SCIPvarGetName(col->var);

      col->flushedobj = obj[pos];
      col->flushedlb = lb[pos];
      col->flushedub = ub[pos];

      for( i = 0; i < col->nlprows; ++i )
      {
         assert(col->rows[i] != NULL);
         lpipos = col->rows[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->nrows);
            assert(nnonz < naddcoefs);
            ind[nnonz] = lpipos;
            val[nnonz] = col->vals[i];
            nnonz++;
         }
      }
#ifndef NDEBUG
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lpipos == -1); /* because the row deletions are already performed */
      }
#endif
   }

   /* call LP interface */
   SCIPsetDebugMsg(set, "flushing col additions: enlarge LP from %d to %d columns\n", lp->nlpicols, lp->ncols);
   SCIP_CALL( SCIPlpiAddCols(lp->lpi, naddcols, obj, lb, ub, name, nnonz, beg, ind, val) );
   lp->nlpicols = lp->ncols;
   lp->lpifirstchgcol = lp->nlpicols;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &name);
   SCIPsetFreeBufferArray(set, &val);
   SCIPsetFreeBufferArray(set, &ind);
   SCIPsetFreeBufferArray(set, &beg);
   SCIPsetFreeBufferArray(set, &ub);
   SCIPsetFreeBufferArray(set, &lb);
   SCIPsetFreeBufferArray(set, &obj);

   lp->flushaddedcols = TRUE;
   lp->updateintegrality = TRUE;

   /* mark the LP unsolved */
   lp->solved = FALSE;
   lp->dualfeasible = FALSE;
   lp->dualchecked = FALSE;
   lp->lpobjval = SCIP_INVALID;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

   return SCIP_OKAY;
}

/** resets row data to represent a row not in the LP solver */
static
void markRowDeleted(
   SCIP_ROW*             row                 /**< row to be marked deleted */
   )
{
   assert(row != NULL);

   row->lpipos = -1;
   row->dualsol = 0.0;
   row->activity = SCIP_INVALID;
   row->dualfarkas = 0.0;
   row->basisstatus = SCIP_BASESTAT_BASIC; /*lint !e641*/
   row->validactivitylp = -1;
}

/** applies all cached row removals to the LP solver */
static
SCIP_RETCODE lpFlushDelRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(lp != NULL);
   assert(lp->lpifirstchgrow <= lp->nlpirows);
   assert(lp->lpifirstchgrow <= lp->nrows);

   /* find the first row to change */
   while( lp->lpifirstchgrow < lp->nlpirows
      && lp->lpifirstchgrow < lp->nrows
      && lp->rows[lp->lpifirstchgrow]->lpipos == lp->lpifirstchgrow
      && !lp->rows[lp->lpifirstchgrow]->coefchanged )
   {
      assert(lp->rows[lp->lpifirstchgrow] == lp->lpirows[lp->lpifirstchgrow]);
      lp->lpifirstchgrow++;
   }

   /* shrink LP to the part which didn't change */
   if( lp->lpifirstchgrow < lp->nlpirows )
   {
      int i;

      SCIPsetDebugMsg(set, "flushing row deletions: shrink LP from %d to %d rows\n", lp->nlpirows, lp->lpifirstchgrow);
      SCIP_CALL( SCIPlpiDelRows(lp->lpi, lp->lpifirstchgrow, lp->nlpirows-1) );
      for( i = lp->lpifirstchgrow; i < lp->nlpirows; ++i )
      {
         markRowDeleted(lp->lpirows[i]);
         SCIP_CALL( SCIProwRelease(&lp->lpirows[i], blkmem, set, lp) );
      }
      lp->nlpirows = lp->lpifirstchgrow;
      lp->flushdeletedrows = TRUE;

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->dualchecked = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   assert(lp->nlpirows == lp->lpifirstchgrow);

   return SCIP_OKAY;
}

/** applies all cached row additions and removals to the LP solver */
static
SCIP_RETCODE lpFlushAddRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   int* beg;
   int* ind;
   SCIP_Real* val;
   char** name;
   SCIP_ROW* row;
   SCIP_Real lpiinf;
   int r;
   int pos;
   int nnonz;
   int naddrows;
   int naddcoefs;
   int i;
   int lpipos;

   assert(lp != NULL);
   assert(lp->lpifirstchgrow == lp->nlpirows);
   assert(blkmem != NULL);

   /* if there are no rows to add, we are ready */
   if( lp->nrows == lp->nlpirows )
      return SCIP_OKAY;

   /* add the additional rows */
   assert(lp->nrows > lp->nlpirows);
   SCIP_CALL( ensureLpirowsSize(lp, set, lp->nrows) );

   /* get the solver's infinity value */
   lpiinf = SCIPlpiInfinity(lp->lpi);

   /* count the (maximal) number of added coefficients, calculate the number of added rows */
   naddrows = lp->nrows - lp->nlpirows;
   naddcoefs = 0;
   for( r = lp->nlpirows; r < lp->nrows; ++r )
      naddcoefs += lp->rows[r]->len;
   assert(naddrows > 0);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lhs, naddrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rhs, naddrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &beg, naddrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &val, naddcoefs) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, naddrows) );

   /* fill temporary memory with row data */
   nnonz = 0;
   for( pos = 0, r = lp->nlpirows; r < lp->nrows; ++pos, ++r )
   {
      row = lp->rows[r];
      assert(row != NULL);
      assert(row->lppos == r);
      assert(nnonz + row->nlpcols <= naddcoefs);

      SCIPsetDebugMsg(set, "flushing added row <%s>: ", row->name);
      debugRowPrint(set, row);

      /* Because the row becomes a member of the LP solver, its dual variable now can take values
       * different from zero. That means, we have to include the row in the corresponding
       * column vectors.
       */
      SCIP_CALL( rowLink(row, blkmem, set, eventqueue, lp) );

      SCIProwCapture(row);
      lp->lpirows[r] = row;
      row->lpipos = r;
      row->dualsol = SCIP_INVALID;
      row->activity = SCIP_INVALID;
      row->dualfarkas = SCIP_INVALID;
      row->validactivitylp = -1;
      row->lhschanged = FALSE;
      row->rhschanged = FALSE;
      row->coefchanged = FALSE;
      if( SCIPsetIsInfinity(set, -row->lhs) )
         lhs[pos] = -lpiinf;
      else
         lhs[pos] = row->lhs - row->constant;
      if( SCIPsetIsInfinity(set, row->rhs) )
         rhs[pos] = lpiinf;
      else
         rhs[pos] = row->rhs - row->constant;
      beg[pos] = nnonz;
      name[pos] = row->name;

      row->flushedlhs = lhs[pos];
      row->flushedrhs = rhs[pos];

      SCIPsetDebugMsg(set, "flushing added row (SCIP_LPI): %+g <=", lhs[pos]);
      for( i = 0; i < row->nlpcols; ++i )
      {
         assert(row->cols[i] != NULL);
         lpipos = row->cols[i]->lpipos;
         if( lpipos >= 0 )
         {
            assert(lpipos < lp->ncols);
            assert(nnonz < naddcoefs);
            SCIPsetDebugMsgPrint(set, " %+gx%d(<%s>)", row->vals[i], lpipos+1, SCIPvarGetName(row->cols[i]->var));
            ind[nnonz] = lpipos;
            val[nnonz] = row->vals[i];
            nnonz++;
         }
      }
      SCIPsetDebugMsgPrint(set, " <= %+g\n", rhs[pos]);
#ifndef NDEBUG
      for( i = row->nlpcols; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->lpipos == -1); /* because the column deletions are already performed */
      }
#endif
   }

   /* call LP interface */
   SCIPsetDebugMsg(set, "flushing row additions: enlarge LP from %d to %d rows\n", lp->nlpirows, lp->nrows);
   SCIP_CALL( SCIPlpiAddRows(lp->lpi, naddrows, lhs, rhs, name, nnonz, beg, ind, val) );
   lp->nlpirows = lp->nrows;
   lp->lpifirstchgrow = lp->nlpirows;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &name);
   SCIPsetFreeBufferArray(set, &val);
   SCIPsetFreeBufferArray(set, &ind);
   SCIPsetFreeBufferArray(set, &beg);
   SCIPsetFreeBufferArray(set, &rhs);
   SCIPsetFreeBufferArray(set, &lhs);

   lp->flushaddedrows = TRUE;

   /* mark the LP unsolved */
   lp->solved = FALSE;
   lp->primalfeasible = FALSE;
   lp->primalchecked = FALSE;
   lp->lpobjval = SCIP_INVALID;
   lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

   return SCIP_OKAY;
}

/** applies all cached column bound and objective changes to the LP */
static
SCIP_RETCODE lpFlushChgCols(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
#ifndef NDEBUG
   SCIP_Bool lpinone = (strcmp( SCIPlpiGetSolverName(), "NONE") == 0);
#endif
   SCIP_COL* col;
   int* objind;
   int* bdind;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real lpiinf;
   int nobjchg;
   int nbdchg;
   int i;

   assert(lp != NULL);

   if( lp->nchgcols == 0 )
      return SCIP_OKAY;

   /* get the solver's infinity value */
   lpiinf = SCIPlpiInfinity(lp->lpi);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &objind, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &obj, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bdind, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lb, lp->ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ub, lp->ncols) );

   /* collect all cached bound and objective changes */
   nobjchg = 0;
   nbdchg = 0;
   for( i = 0; i < lp->nchgcols; ++i )
   {
      col = lp->chgcols[i];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);

      if( col->lpipos >= 0 )
      {
#ifndef NDEBUG
         /* do not check consistency of data with LPI in case of LPI=none */
         if( !lpinone )
         {
            SCIP_Real lpiobj;
            SCIP_Real lpilb;
            SCIP_Real lpiub;

            SCIP_CALL( SCIPlpiGetObj(lp->lpi, col->lpipos, col->lpipos, &lpiobj) );
            SCIP_CALL( SCIPlpiGetBounds(lp->lpi, col->lpipos, col->lpipos, &lpilb, &lpiub) );
            assert(SCIPsetIsFeasEQ(set, lpiobj, col->flushedobj));
            assert((SCIPsetIsInfinity(set, -lpilb) && SCIPsetIsInfinity(set, -col->flushedlb))
                  || (!SCIPsetIsInfinity(set, -lpilb) && !SCIPsetIsInfinity(set, -col->flushedlb) && SCIPsetIsFeasEQ(set, lpilb, col->flushedlb)));
            assert((SCIPsetIsInfinity(set, lpiub) && SCIPsetIsInfinity(set, col->flushedub))
                  || (!SCIPsetIsInfinity(set, lpiub) && !SCIPsetIsInfinity(set, col->flushedub) && SCIPsetIsFeasEQ(set, lpiub, col->flushedub)));
         }
#endif

         if( col->objchanged )
         {
            SCIP_Real newobj;

            newobj = col->obj;
            if( col->flushedobj != newobj ) /*lint !e777*/
            {
               assert(nobjchg < lp->ncols);
               objind[nobjchg] = col->lpipos;
               obj[nobjchg] = newobj;
               nobjchg++;
               col->flushedobj = newobj;
            }
            col->objchanged = FALSE;
         }

         if( col->lbchanged || col->ubchanged )
         {
            SCIP_Real newlb;
            SCIP_Real newub;

            /* compute bounds that should be flushed into the LP (taking into account lazy bounds) */
            computeLPBounds(lp, set, col, lpiinf, &newlb, &newub);

            if( col->flushedlb != newlb || col->flushedub != newub ) /*lint !e777*/
            {
               assert(nbdchg < lp->ncols);
               bdind[nbdchg] = col->lpipos;
               lb[nbdchg] = newlb;
               ub[nbdchg] = newub;
               nbdchg++;
               col->flushedlb = newlb;
               col->flushedub = newub;
            }
            col->lbchanged = FALSE;
            col->ubchanged = FALSE;
         }
      }
      /* maybe lb/ub/objchanged should all be set to false when lpipos is -1 */
   }

   /* change objective values in LP */
   if( nobjchg > 0 )
   {
      SCIPsetDebugMsg(set, "flushing objective changes: change %d objective values of %d changed columns\n", nobjchg, lp->nchgcols);
      SCIP_CALL( SCIPlpiChgObj(lp->lpi, nobjchg, objind, obj) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->dualchecked = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   /* change bounds in LP */
   if( nbdchg > 0 )
   {
      SCIPsetDebugMsg(set, "flushing bound changes: change %d bounds of %d changed columns\n", nbdchg, lp->nchgcols);
      SCIP_CALL( SCIPlpiChgBounds(lp->lpi, nbdchg, bdind, lb, ub) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->primalchecked = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   lp->nchgcols = 0;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &ub);
   SCIPsetFreeBufferArray(set, &lb);
   SCIPsetFreeBufferArray(set, &bdind);
   SCIPsetFreeBufferArray(set, &obj);
   SCIPsetFreeBufferArray(set, &objind);

   return SCIP_OKAY;
}

/** applies all cached row side changes to the LP */
static
SCIP_RETCODE lpFlushChgRows(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
#ifndef NDEBUG
   SCIP_Bool lpinone = (strcmp( SCIPlpiGetSolverName(), "NONE") == 0);
#endif
   SCIP_ROW* row;
   int* ind;
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   SCIP_Real lpiinf;
   int i;
   int nchg;

   assert(lp != NULL);

   if( lp->nchgrows == 0 )
      return SCIP_OKAY;

   /* get the solver's infinity value */
   lpiinf = SCIPlpiInfinity(lp->lpi);

   /* get temporary memory for changes */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lhs, lp->nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rhs, lp->nrows) );

   /* collect all cached left and right hand side changes */
   nchg = 0;
   for( i = 0; i < lp->nchgrows; ++i )
   {
      row = lp->chgrows[i];
      assert(row != NULL);

      if( row->lpipos >= 0 )
      {
#ifndef NDEBUG
         /* do not check consistency of data with LPI in case of LPI=none */
         if( !lpinone )
         {
            SCIP_Real lpilhs;
            SCIP_Real lpirhs;

            SCIP_CALL( SCIPlpiGetSides(lp->lpi, row->lpipos, row->lpipos, &lpilhs, &lpirhs) );
            assert(SCIPsetIsSumEQ(set, lpilhs, row->flushedlhs));
            assert(SCIPsetIsSumEQ(set, lpirhs, row->flushedrhs));
         }
#endif
         if( row->lhschanged || row->rhschanged )
         {
            SCIP_Real newlhs;
            SCIP_Real newrhs;

            newlhs = (SCIPsetIsInfinity(set, -row->lhs) ? -lpiinf : row->lhs - row->constant);
            newrhs = (SCIPsetIsInfinity(set, row->rhs) ? lpiinf : row->rhs - row->constant);
            if( row->flushedlhs != newlhs || row->flushedrhs != newrhs ) /*lint !e777*/
            {
               assert(nchg < lp->nrows);
               ind[nchg] = row->lpipos;
               lhs[nchg] = newlhs;
               rhs[nchg] = newrhs;
               nchg++;
               row->flushedlhs = newlhs;
               row->flushedrhs = newrhs;
            }
            row->lhschanged = FALSE;
            row->rhschanged = FALSE;
         }
      }
   }

   /* change left and right hand sides in LP */
   if( nchg > 0 )
   {
      SCIPsetDebugMsg(set, "flushing side changes: change %d sides of %d rows\n", nchg, lp->nchgrows);
      SCIP_CALL( SCIPlpiChgSides(lp->lpi, nchg, ind, lhs, rhs) );

      /* mark the LP unsolved */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->primalchecked = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   lp->nchgrows = 0;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &rhs);
   SCIPsetFreeBufferArray(set, &lhs);
   SCIPsetFreeBufferArray(set, &ind);

   return SCIP_OKAY;
}

/** copy integrality information to the LP */
static
SCIP_RETCODE lpCopyIntegrality(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   int nintegers;
   int* integerInfo;
   SCIP_VAR* var;

   assert(lp != NULL);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &integerInfo, lp->ncols) );

   /* count total number of integralities */
   nintegers = 0;

   for( i = 0; i < lp->ncols; ++i )
   {
      var = SCIPcolGetVar(lp->cols[i]);
      if( SCIPvarIsIntegral(var) || SCIPvarIsBinary(var) )
      {
         integerInfo[i] = 1;
         ++nintegers;
      }
      else
         integerInfo[i] = 0;
   }

   /* only pass integrality information if integer variables are present */
   if( nintegers > 0 )
   {
      SCIP_CALL( SCIPlpiSetIntegralityInformation(lp->lpi, lp->ncols, integerInfo) );
   }
   else
   {
      SCIP_CALL( SCIPlpiSetIntegralityInformation(lp->lpi, 0, NULL) );
   }

   SCIPsetFreeBufferArray(set, &integerInfo);

   /* mark integralities to be updated */
   lp->updateintegrality = FALSE;

   return SCIP_OKAY;
}

/** applies all cached changes to the LP solver */
SCIP_RETCODE SCIPlpFlush(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(lp != NULL);
   assert(blkmem != NULL);

   SCIPsetDebugMsg(set, "flushing LP changes: old (%d cols, %d rows), nchgcols=%d, nchgrows=%d, firstchgcol=%d, firstchgrow=%d, new (%d cols, %d rows), flushed=%u\n",
      lp->nlpicols, lp->nlpirows, lp->nchgcols, lp->nchgrows, lp->lpifirstchgcol, lp->lpifirstchgrow, lp->ncols, lp->nrows, lp->flushed);

   if( !lp->flushed )
   {
      lp->flushdeletedcols = FALSE;
      lp->flushaddedcols = FALSE;
      lp->flushdeletedrows = FALSE;
      lp->flushaddedrows = FALSE;

      SCIP_CALL( lpFlushDelCols(lp) );
      SCIP_CALL( lpFlushDelRows(lp, blkmem, set) );
      SCIP_CALL( lpFlushChgCols(lp, set) );
      SCIP_CALL( lpFlushChgRows(lp, set) );
      SCIP_CALL( lpFlushAddCols(lp, blkmem, set, eventqueue) );
      SCIP_CALL( lpFlushAddRows(lp, blkmem, set, eventqueue) );

      lp->flushed = TRUE;

      checkLinks(lp);
   }

   /* if the cutoff bound was changed in between, we want to re-optimize the LP even if nothing else has changed */
   if( lp->cutoffbound != lp->lpiobjlim && lp->ncols > 0 ) /*lint !e777*/
   {
      lp->solved = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   assert(lp->nlpicols == lp->ncols);
   assert(lp->lpifirstchgcol == lp->nlpicols);
   assert(lp->nlpirows == lp->nrows);
   assert(lp->lpifirstchgrow == lp->nlpirows);
   assert(lp->nchgcols == 0);
   assert(lp->nchgrows == 0);
#ifndef NDEBUG
   {
      int ncols;
      int nrows;

      SCIP_CALL( SCIPlpiGetNCols(lp->lpi, &ncols) );
      SCIP_CALL( SCIPlpiGetNRows(lp->lpi, &nrows) );
      assert(ncols == lp->ncols);
      assert(nrows == lp->nrows);
   }
#endif

   return SCIP_OKAY;
}

/** marks the LP to be flushed, even if the LP thinks it is not flushed */
SCIP_RETCODE SCIPlpMarkFlushed(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
#ifndef NDEBUG
   SCIP_Bool lpinone = (strcmp( SCIPlpiGetSolverName(), "NONE") == 0);
#endif
   int i;

   assert(lp != NULL);

#ifndef NDEBUG
   /* check, if there are really no column or row deletions or coefficient changes left */
   while( lp->lpifirstchgcol < lp->nlpicols
      && lp->lpifirstchgcol < lp->ncols
      && lp->cols[lp->lpifirstchgcol]->lpipos == lp->lpifirstchgcol
      && !lp->cols[lp->lpifirstchgcol]->coefchanged )
   {
      assert(lp->cols[lp->lpifirstchgcol] == lp->lpicols[lp->lpifirstchgcol]);
      lp->lpifirstchgcol++;
   }
   assert(lp->nlpicols == lp->lpifirstchgcol);

   while( lp->lpifirstchgrow < lp->nlpirows
      && lp->lpifirstchgrow < lp->nrows
      && lp->rows[lp->lpifirstchgrow]->lpipos == lp->lpifirstchgrow
      && !lp->rows[lp->lpifirstchgrow]->coefchanged )
   {
      assert(lp->rows[lp->lpifirstchgrow] == lp->lpirows[lp->lpifirstchgrow]);
      lp->lpifirstchgrow++;
   }
   assert(lp->nlpirows == lp->lpifirstchgrow);
#endif

   lp->lpifirstchgcol = lp->nlpicols;
   lp->lpifirstchgrow = lp->nlpirows;

   /* check, if there are really no column or row additions left */
   assert(lp->ncols == lp->nlpicols);
   assert(lp->nrows == lp->nlpirows);

   /* mark the changed columns to be unchanged, and check, if this is really correct */
   for( i = 0; i < lp->nchgcols; ++i )
   {
      SCIP_COL* col;

      col = lp->chgcols[i];
      assert(col != NULL);
      assert(col->var != NULL);
      assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(col->var) == col);

      if( col->lpipos >= 0 )
      {
#ifndef NDEBUG
         /* do not check consistency of data with LPI in case of LPI=none */
         if( !lpinone )
         {
            SCIP_Real lpiobj;
            SCIP_Real lpilb;
            SCIP_Real lpiub;

            SCIP_CALL( SCIPlpiGetObj(lp->lpi, col->lpipos, col->lpipos, &lpiobj) );
            SCIP_CALL( SCIPlpiGetBounds(lp->lpi, col->lpipos, col->lpipos, &lpilb, &lpiub) );
            assert(SCIPsetIsSumEQ(set, lpiobj, col->flushedobj));
            assert(SCIPsetIsSumEQ(set, lpilb, col->flushedlb));
            assert(SCIPsetIsSumEQ(set, lpiub, col->flushedub));
            assert(col->flushedobj == col->obj); /*lint !e777*/
            assert(col->flushedlb == (SCIPsetIsInfinity(set, -col->lb) ? -SCIPlpiInfinity(lp->lpi) : col->lb)); /*lint !e777*/
            assert(col->flushedub == (SCIPsetIsInfinity(set, col->ub) ? SCIPlpiInfinity(lp->lpi) : col->ub)); /*lint !e777*/
         }
#endif
         col->objchanged = FALSE;
         col->lbchanged = FALSE;
         col->ubchanged = FALSE;
      }
      /* maybe lb/ub/objchanged should  be set to false also when lpipos is -1 */
   }
   lp->nchgcols = 0;

   /* mark the changed rows to be unchanged, and check, if this is really correct */
   for( i = 0; i < lp->nchgrows; ++i )
   {
      SCIP_ROW* row;

      row = lp->chgrows[i];
      assert(row != NULL);

      if( row->lpipos >= 0 )
      {
#ifndef NDEBUG
         /* do not check consistency of data with LPI in case of LPI=none */
         if( !lpinone )
         {
            SCIP_Real lpilhs;
            SCIP_Real lpirhs;

            SCIP_CALL( SCIPlpiGetSides(lp->lpi, row->lpipos, row->lpipos, &lpilhs, &lpirhs) );
            assert(SCIPsetIsSumEQ(set, lpilhs, row->flushedlhs));
            assert(SCIPsetIsSumEQ(set, lpirhs, row->flushedrhs));
            assert(row->flushedlhs == (SCIPsetIsInfinity(set, -row->lhs) ? -SCIPlpiInfinity(lp->lpi) : row->lhs - row->constant)); /*lint !e777*/
            assert(row->flushedrhs == (SCIPsetIsInfinity(set, row->rhs) ? SCIPlpiInfinity(lp->lpi) : row->rhs - row->constant)); /*lint !e777*/
         }
#endif
         row->lhschanged = FALSE;
         row->rhschanged = FALSE;
      }
   }
   lp->nchgrows = 0;

   /* mark the LP to be flushed */
   lp->flushed = TRUE;

   checkLinks(lp);

   return SCIP_OKAY;
}




/*
 * LP methods
 */

/** updates link data after addition of column */
static
void colUpdateAddLP(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_ROW* row;
   int i;
   int pos;

   assert(col != NULL);
   assert(col->lppos >= 0);

   /* update column arrays of all linked rows */
   for( i = 0; i < col->len; ++i )
   {
      pos = col->linkpos[i];
      if( pos >= 0 )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->linkpos[pos] == i);
         assert(row->cols[pos] == col);
         assert(row->nlpcols <= pos && pos < row->len);

         row->nlpcols++;
         rowSwapCoefs(row, pos, row->nlpcols-1);
         assert(row->cols[row->nlpcols-1] == col);

         /* if no swap was necessary, mark lpcols to be unsorted */
         if( pos == row->nlpcols-1 )
            row->lpcolssorted = FALSE;

         /* update norms */
         rowAddNorms(row, set, col, row->vals[row->nlpcols-1], FALSE);
      }
   }
}

/** updates link data after addition of row */
static
void rowUpdateAddLP(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_COL* col;
   int i;
   int pos;

   assert(row != NULL);
   assert(row->lppos >= 0);

   /* update row arrays of all linked columns */
   for( i = 0; i < row->len; ++i )
   {
      pos = row->linkpos[i];
      if( pos >= 0 )
      {
         col = row->cols[i];
         assert(col != NULL);
         assert(col->linkpos[pos] == i);
         assert(col->rows[pos] == row);
         assert(col->nlprows <= pos && pos < col->len);

         col->nlprows++;
         colSwapCoefs(col, pos, col->nlprows-1);

         /* if no swap was necessary, mark lprows to be unsorted */
         if( pos == col->nlprows-1 )
            col->lprowssorted = FALSE;
      }
   }
}

/** updates link data after removal of column */
static
void colUpdateDelLP(
   SCIP_COL*             col,                /**< LP column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_ROW* row;
   int i;
   int pos;

   assert(col != NULL);
   assert(col->lppos == -1);

   /* update column arrays of all linked rows */
   for( i = 0; i < col->len; ++i )
   {
      pos = col->linkpos[i];
      if( pos >= 0 )
      {
         row = col->rows[i];
         assert(row != NULL);
         assert(row->linkpos[pos] == i);
         assert(row->cols[pos] == col);
         assert(0 <= pos && pos < row->nlpcols);

         /* update norms */
         rowDelNorms(row, set, col, row->vals[pos], TRUE, FALSE, FALSE);

         row->nlpcols--;
         rowSwapCoefs(row, pos, row->nlpcols);

         /* if no swap was necessary, mark nonlpcols to be unsorted */
         if( pos == row->nlpcols )
            row->nonlpcolssorted = FALSE;
      }
   }
}

/** updates link data after removal of row */
static
void rowUpdateDelLP(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_COL* col;
   int i;
   int pos;

   assert(row != NULL);
   assert(row->lppos == -1);

   /* update row arrays of all linked columns */
   for( i = 0; i < row->len; ++i )
   {
      pos = row->linkpos[i];
      if( pos >= 0 )
      {
         col = row->cols[i];
         assert(col != NULL);
         assert(0 <= pos && pos < col->nlprows);
         assert(col->linkpos[pos] == i);
         assert(col->rows[pos] == row);

         col->nlprows--;
         colSwapCoefs(col, pos, col->nlprows);

         /* if no swap was necessary, mark lprows to be unsorted */
         if( pos == col->nlprows )
            col->nonlprowssorted = FALSE;
      }
   }
}

static
SCIP_RETCODE allocDiveChgSideArrays(
   SCIP_LP*              lp,                 /**< LP data object */
   int                   initsize            /**< initial size of the arrays */
   )
{
   assert(lp != NULL);
   assert(lp->divechgsides == NULL);
   assert(lp->divechgsidetypes == NULL);
   assert(lp->divechgrows == NULL);
   assert(lp->ndivechgsides == 0);
   assert(lp->divechgsidessize == 0);
   assert(initsize > 0);

   lp->divechgsidessize = initsize;
   SCIP_ALLOC( BMSallocMemoryArray(&lp->divechgsides, lp->divechgsidessize) );
   SCIP_ALLOC( BMSallocMemoryArray(&lp->divechgsidetypes, lp->divechgsidessize) );
   SCIP_ALLOC( BMSallocMemoryArray(&lp->divechgrows, lp->divechgsidessize) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE reallocDiveChgSideArrays(
   SCIP_LP*              lp,                 /**< LP data object */
   int                   minsize,            /**< minimal number of elements */
   SCIP_Real             growfact            /**< growing factor */
   )
{
   assert(lp != NULL);
   assert(lp->divechgsides != NULL);
   assert(lp->divechgsidetypes != NULL);
   assert(lp->divechgrows != NULL);
   assert(lp->ndivechgsides > 0);
   assert(lp->divechgsidessize > 0);
   assert(minsize > 0);

   if( minsize <= lp->divechgsidessize )
      return SCIP_OKAY;

   lp->divechgsidessize = MAX(minsize, (int)(lp->divechgsidessize * growfact));
   SCIP_ALLOC( BMSreallocMemoryArray(&lp->divechgsides, lp->divechgsidessize) );
   SCIP_ALLOC( BMSreallocMemoryArray(&lp->divechgsidetypes, lp->divechgsidessize) );
   SCIP_ALLOC( BMSreallocMemoryArray(&lp->divechgrows, lp->divechgsidessize) );

   return SCIP_OKAY;
}

static
void freeDiveChgSideArrays(
   SCIP_LP*              lp                  /**< LP data object */
   )
{
   assert(lp != NULL);
   assert(lp->divechgsides != NULL);
   assert(lp->divechgsidetypes != NULL);
   assert(lp->divechgrows != NULL);
   assert(lp->ndivechgsides == 0);
   assert(lp->divechgsidessize > 0);

   BMSfreeMemoryArrayNull(&lp->divechgsides);
   BMSfreeMemoryArrayNull(&lp->divechgsidetypes);
   BMSfreeMemoryArrayNull(&lp->divechgrows);
   lp->divechgsidessize = 0;
}

#define DIVESTACKINITSIZE 100

/** creates empty LP data object */
SCIP_RETCODE SCIPlpCreate(
   SCIP_LP**             lp,                 /**< pointer to LP data object */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name                /**< problem name */
   )
{
   SCIP_Bool success;

   assert(lp != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(name != NULL);

   SCIP_ALLOC( BMSallocMemory(lp) );

   /* open LP Solver interface */
   SCIP_CALL( SCIPlpiCreate(&(*lp)->lpi, messagehdlr, name, SCIP_OBJSEN_MINIMIZE) );

   (*lp)->lpicols = NULL;
   (*lp)->lpirows = NULL;
   (*lp)->chgcols = NULL;
   (*lp)->chgrows = NULL;
   (*lp)->cols = NULL;
   (*lp)->lazycols = NULL;
   (*lp)->rows = NULL;
   (*lp)->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   (*lp)->lpobjval = 0.0;
   (*lp)->glbpseudoobjval = 0.0;
   (*lp)->relglbpseudoobjval = 0.0;
   (*lp)->glbpseudoobjvalid = TRUE;
   (*lp)->glbpseudoobjvalinf = 0;
   (*lp)->pseudoobjval = 0.0;
   (*lp)->relpseudoobjval = 0.0;
   (*lp)->pseudoobjvalid = TRUE;
   (*lp)->pseudoobjvalinf = 0;
   (*lp)->looseobjval = 0.0;
   (*lp)->rellooseobjval = 0.0;
   (*lp)->looseobjvalid = TRUE;
   (*lp)->looseobjvalinf = 0;
   (*lp)->nloosevars = 0;
   (*lp)->rootlpobjval = SCIP_INVALID;
   (*lp)->rootlooseobjval = SCIP_INVALID;
   (*lp)->cutoffbound = SCIPsetInfinity(set);
   (*lp)->objsqrnorm = 0.0;
   (*lp)->objsumnorm = 0.0;
   (*lp)->lpicolssize = 0;
   (*lp)->nlpicols = 0;
   (*lp)->lpirowssize = 0;
   (*lp)->nlpirows = 0;
   (*lp)->lpifirstchgcol = 0;
   (*lp)->lpifirstchgrow = 0;
   (*lp)->colssize = 0;
   (*lp)->ncols = 0;
   (*lp)->lazycolssize = 0;
   (*lp)->nlazycols = 0;
   (*lp)->rowssize = 0;
   (*lp)->nrows = 0;
   (*lp)->chgcolssize = 0;
   (*lp)->nchgcols = 0;
   (*lp)->chgrowssize = 0;
   (*lp)->nchgrows = 0;
   (*lp)->firstnewcol = 0;
   (*lp)->firstnewrow = 0;
   (*lp)->nremovablecols = 0;
   (*lp)->nremovablerows = 0;
   (*lp)->validsollp = stat->lpcount; /* the initial (empty) SCIP_LP is solved with primal and dual solution of zero */
   (*lp)->validfarkaslp = -1;
   (*lp)->objsqrnormunreliable = FALSE;
   (*lp)->flushdeletedcols = FALSE;
   (*lp)->flushaddedcols = FALSE;
   (*lp)->flushdeletedrows = FALSE;
   (*lp)->flushaddedrows = FALSE;
   (*lp)->updateintegrality = TRUE;
   (*lp)->flushed = TRUE;
   (*lp)->solved = TRUE;
   (*lp)->primalfeasible = TRUE;
   (*lp)->primalchecked = TRUE;
   (*lp)->dualfeasible = TRUE;
   (*lp)->dualchecked = TRUE;
   (*lp)->solisbasic = FALSE;
   (*lp)->rootlpisrelax = TRUE;
   (*lp)->isrelax = TRUE;
   (*lp)->installing = FALSE;
   (*lp)->strongbranching = FALSE;
   (*lp)->strongbranchprobing = FALSE;
   (*lp)->probing = FALSE;
   (*lp)->diving = FALSE;
   (*lp)->divingobjchg = FALSE;
   (*lp)->divinglazyapplied = FALSE;
   (*lp)->divelpistate = NULL;
   (*lp)->divelpwasprimfeas = TRUE;
   (*lp)->divelpwasprimchecked = TRUE;
   (*lp)->divelpwasdualfeas = TRUE;
   (*lp)->divelpwasdualchecked = TRUE;
   (*lp)->divechgsides = NULL;
   (*lp)->divechgsidetypes = NULL;
   (*lp)->divechgrows = NULL;
   (*lp)->ndivechgsides = 0;
   (*lp)->divechgsidessize = 0;
   (*lp)->ndivingrows = 0;
   (*lp)->divinglpiitlim = INT_MAX;
   (*lp)->resolvelperror = FALSE;
   (*lp)->divenolddomchgs = 0;
   (*lp)->adjustlpval = FALSE;
   (*lp)->lpiobjlim = SCIPlpiInfinity((*lp)->lpi);
   (*lp)->lpifeastol = SCIPsetLpfeastol(set);
   (*lp)->lpidualfeastol = SCIPsetDualfeastol(set);
   (*lp)->lpibarrierconvtol = SCIPsetBarrierconvtol(set);
   (*lp)->lpifromscratch = FALSE;
   (*lp)->lpifastmip = set->lp_fastmip;
   (*lp)->lpiscaling = set->lp_scaling;
   (*lp)->lpipresolving = set->lp_presolving;
   (*lp)->lpilpinfo = set->disp_lpinfo;
   (*lp)->lpirowrepswitch = set->lp_rowrepswitch;
   (*lp)->lpisolutionpolishing = (set->lp_solutionpolishing > 0);
   (*lp)->lpirefactorinterval = set->lp_refactorinterval;
   (*lp)->lpiconditionlimit = set->lp_conditionlimit;
   (*lp)->lpiitlim = INT_MAX;
   (*lp)->lpipricing = SCIP_PRICING_AUTO;
   (*lp)->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;
   (*lp)->lpithreads = set->lp_threads;
   (*lp)->lpitiming = (int) set->time_clocktype;
   (*lp)->lpirandomseed = set->random_randomseed;
   (*lp)->storedsolvals = NULL;

   /* allocate arrays for diving */
   SCIP_CALL( allocDiveChgSideArrays(*lp, DIVESTACKINITSIZE) );

   /* set default parameters in LP solver */
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_OBJLIM, (*lp)->lpiobjlim, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: objective limit cannot be set -- can lead to unnecessary simplex iterations\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_FEASTOL, (*lp)->lpifeastol, &success) );
   (*lp)->lpihasfeastol = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: primal feasibility tolerance cannot be set -- tolerance of SCIP and LP solver may differ\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_DUALFEASTOL, (*lp)->lpidualfeastol, &success) );
   (*lp)->lpihasdualfeastol = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: dual feasibility tolerance cannot be set -- tolerance of SCIP and LP solver may differ\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_BARRIERCONVTOL, (*lp)->lpibarrierconvtol, &success) );
   (*lp)->lpihasbarrierconvtol = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: barrier convergence tolerance cannot be set -- tolerance of SCIP and LP solver may differ\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetBoolpar(*lp, SCIP_LPPAR_FROMSCRATCH, (*lp)->lpifromscratch, &success) );
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_FASTMIP, (*lp)->lpifastmip, &success) );
   (*lp)->lpihasfastmip = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: fastmip setting not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_SCALING, (*lp)->lpiscaling, &success) );
   (*lp)->lpihasscaling = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: scaling not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetBoolpar(*lp, SCIP_LPPAR_PRESOLVING, (*lp)->lpipresolving, &success) );
   (*lp)->lpihaspresolving = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: presolving not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_TIMING, (*lp)->lpitiming, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: clock type cannot be set\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_LPITLIM, (*lp)->lpiitlim, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: iteration limit cannot be set -- can lead to unnecessary simplex iterations\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_PRICING, (int)(*lp)->lpipricing, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: pricing strategy cannot be set -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetBoolpar(*lp, SCIP_LPPAR_LPINFO, (*lp)->lpilpinfo, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: lpinfo setting not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_ROWREPSWITCH, (*lp)->lpirowrepswitch, &success) );
   (*lp)->lpihasrowrep = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: row representation of the basis not available -- SCIP parameter lp/rowrepswitch has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_POLISHING, ((*lp)->lpisolutionpolishing ? 1 : 0), &success) );
   (*lp)->lpihaspolishing = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: solution polishing not available -- SCIP parameter lp/solutionpolishing has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_REFACTOR, (*lp)->lpirefactorinterval, &success) );
   (*lp)->lpihasrefactor = success;
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: refactorization interval not available -- SCIP parameter lp/refactorinterval has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetRealpar(*lp, SCIP_LPPAR_CONDITIONLIMIT, (*lp)->lpiconditionlimit, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: condition number limit for the basis not available -- SCIP parameter lp/conditionlimit has no effect\n",
         SCIPlpiGetSolverName());
   }
   SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_THREADS, (*lp)->lpithreads, &success) );
   if( !success )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "LP Solver <%s>: number of threads settings not available -- SCIP parameter has no effect\n",
         SCIPlpiGetSolverName());
   }
   /* keep the default LP random seed if this parameter is set to 0 (current default) */
   if( (*lp)->lpirandomseed != 0 )
   {
      SCIP_CALL( lpSetIntpar(*lp, SCIP_LPPAR_RANDOMSEED, (*lp)->lpirandomseed, &success) );
      if( !success )
      {
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "LP Solver <%s>: random seed parameter not available -- SCIP parameter has no effect\n",
            SCIPlpiGetSolverName());
      }
   }

   /* Check that infinity value of LP-solver is at least as large as the one used in SCIP. This is necessary, because we
    * transfer SCIP infinity values to the ones by the LPI, but not the converse. */
   if ( set->num_infinity > SCIPlpiInfinity((*lp)->lpi) )
   {
      SCIPerrorMessage("The infinity value of the LP solver has to be at least as large as the one of SCIP.\n");
      return SCIP_PARAMETERWRONGVAL;
   }

   return SCIP_OKAY;
}

/** frees LP data object */
SCIP_RETCODE SCIPlpFree(
   SCIP_LP**             lp,                 /**< pointer to LP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   int i;

   assert(lp != NULL);
   assert(*lp != NULL);

   SCIP_CALL( SCIPlpClear(*lp, blkmem, set, eventqueue, eventfilter) );

   freeDiveChgSideArrays(*lp);

   /* release LPI rows */
   for( i = 0; i < (*lp)->nlpirows; ++i )
   {
      SCIP_CALL( SCIProwRelease(&(*lp)->lpirows[i], blkmem, set, *lp) );
   }

   if( (*lp)->lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&(*lp)->lpi) );
   }

   BMSfreeMemoryNull(&(*lp)->storedsolvals);
   BMSfreeMemoryArrayNull(&(*lp)->lpicols);
   BMSfreeMemoryArrayNull(&(*lp)->lpirows);
   BMSfreeMemoryArrayNull(&(*lp)->chgcols);
   BMSfreeMemoryArrayNull(&(*lp)->chgrows);
   BMSfreeMemoryArrayNull(&(*lp)->lazycols);
   BMSfreeMemoryArrayNull(&(*lp)->cols);
   BMSfreeMemoryArrayNull(&(*lp)->rows);
   BMSfreeMemory(lp);

   return SCIP_OKAY;
}

/** resets the LP to the empty LP by removing all columns and rows from LP, releasing all rows, and flushing the
 *  changes to the LP solver
 */
SCIP_RETCODE SCIPlpReset(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   assert(stat != NULL);

   SCIP_CALL( SCIPlpClear(lp, blkmem, set, eventqueue, eventfilter) );
   SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );

   /* mark the empty LP to be solved */
   lp->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   lp->lpobjval = 0.0;
   lp->validsollp = stat->lpcount; /* the initial (empty) SCIP_LP is solved with primal and dual solution of zero */
   lp->validfarkaslp = -1;
   lp->solved = TRUE;
   lp->primalfeasible = TRUE;
   lp->primalchecked = TRUE;
   lp->dualfeasible = TRUE;
   lp->dualchecked = TRUE;
   lp->solisbasic = FALSE;
   lp->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;

   return SCIP_OKAY;
}

/** adds a column to the LP */
SCIP_RETCODE SCIPlpAddCol(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COL*             col,                /**< LP column */
   int                   depth               /**< depth in the tree where the column addition is performed */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);
   assert(col != NULL);
   assert(col->len == 0 || col->rows != NULL);
   assert(col->lppos == -1);
   assert(col->var != NULL);
   assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetCol(col->var) == col);
   assert(SCIPvarIsIntegral(col->var) == col->integral);

   SCIPsetDebugMsg(set, "adding column <%s> to LP (%d rows, %d cols)\n", SCIPvarGetName(col->var), lp->nrows, lp->ncols);
#ifdef SCIP_DEBUG
   {
      int i;
      SCIPsetDebugMsgPrint(set, "  (obj: %g) [%g,%g]", col->obj, col->lb, col->ub);
      for( i = 0; i < col->len; ++i )
         SCIPsetDebugMsgPrint(set, " %+g<%s>", col->vals[i], col->rows[i]->name);
      SCIPsetDebugMsgPrint(set, "\n");
   }
#endif

   SCIP_CALL( ensureColsSize(lp, set, lp->ncols+1) );
   lp->cols[lp->ncols] = col;
   col->lppos = lp->ncols;
   col->lpdepth = depth;
   col->age = 0;
   lp->ncols++;
   if( col->removable )
      lp->nremovablecols++;

   if( !SCIPsetIsInfinity(set, -col->lazylb) || !SCIPsetIsInfinity(set, col->lazyub) )
   {
      SCIP_CALL( ensureLazycolsSize(lp, set, lp->nlazycols+1) );
      lp->lazycols[lp->nlazycols] = col;
      lp->nlazycols++;
   }

   /* mark the current LP unflushed */
   lp->flushed = FALSE;

   /* update column arrays of all linked rows */
   colUpdateAddLP(col, set);

   /* update the objective function vector norms */
   lpUpdateObjNorms(lp, set, 0.0, col->unchangedobj);

   checkLinks(lp);

   return SCIP_OKAY;
}

/** adds a row to the LP and captures it */
SCIP_RETCODE SCIPlpAddRow(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_ROW*             row,                /**< LP row */
   int                   depth               /**< depth in the tree where the row addition is performed */
   )
{
   assert(lp != NULL);
   assert(row != NULL);
   assert(row->len == 0 || row->cols != NULL);
   assert(row->lppos == -1);

   SCIProwCapture(row);
   SCIProwLock(row);

   SCIPsetDebugMsg(set, "adding row <%s> to LP (%d rows, %d cols)\n", row->name, lp->nrows, lp->ncols);
#ifdef SCIP_DEBUG
   {
      int i;
      SCIPsetDebugMsgPrint(set, "  %g <=", row->lhs);
      for( i = 0; i < row->len; ++i )
         SCIPsetDebugMsgPrint(set, " %+g<%s>", row->vals[i], SCIPvarGetName(row->cols[i]->var));
      if( !SCIPsetIsZero(set, row->constant) )
         SCIPsetDebugMsgPrint(set, " %+g", row->constant);
      SCIPsetDebugMsgPrint(set, " <= %g\n", row->rhs);
   }
#endif

   SCIP_CALL( ensureRowsSize(lp, set, lp->nrows+1) );
   lp->rows[lp->nrows] = row;
   row->lppos = lp->nrows;
   row->lpdepth = depth;
   row->age = 0;
   lp->nrows++;
   if( row->removable )
      lp->nremovablerows++;

   /* mark the current LP unflushed */
   lp->flushed = FALSE;

   /* update row arrays of all linked columns */
   rowUpdateAddLP(row);

   checkLinks(lp);

   rowCalcNorms(row, set);

   /* check, if row addition to LP events are tracked
    * if so, issue ROWADDEDLP event
    */
   if( (eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWADDEDLP) != 0) )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowAddedLP(&event, blkmem, row) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
   }

   return SCIP_OKAY;
}


#ifndef NDEBUG
/** method checks if all columns in the lazycols array have at least one lazy bound and also have a counter part in the
 *  cols array; furthermore, it is checked if columns in the cols array which have a lazy bound have a counter part in
 *  the lazycols array
 */
static
void checkLazyColArray(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Bool contained;
   int c;
   int i;

   assert(lp != NULL);

   /* check if each column in the lazy column array has a counter part in the column array */
   for( i = 0; i < lp->nlazycols; ++i )
   {
      /* check if each lazy column has at least on lazy bound */
      assert(lp->lazycols[i] != NULL);
      assert(!SCIPsetIsInfinity(set, lp->lazycols[i]->lazyub) || !SCIPsetIsInfinity(set, -lp->lazycols[i]->lazylb));

      contained = FALSE;
      for( c = 0; c < lp->ncols; ++c )
      {
         if( lp->lazycols[i] == lp->cols[c] )
         {
            assert(!SCIPsetIsInfinity(set, lp->cols[c]->lazyub) || !SCIPsetIsInfinity(set, -lp->cols[c]->lazylb));
            contained = TRUE;
         }
      }
      assert(contained);
   }

   /* check if each column in the column array which has at least one lazy bound has a counter part in the lazy column *
    * array */
   for( c = 0; c < lp->ncols; ++c )
   {
      contained = FALSE;
      assert(lp->cols[c] != NULL);

      for( i = 0; i < lp->nlazycols; ++i )
      {
         if( lp->lazycols[i] == lp->cols[c] )
         {
            contained = TRUE;
         }
      }

      assert(contained == (!SCIPsetIsInfinity(set, lp->cols[c]->lazyub) || !SCIPsetIsInfinity(set, -lp->cols[c]->lazylb)));
   }
}
#else
#define checkLazyColArray(lp, set) /**/
#endif

/** removes all columns after the given number of cols from the LP */
SCIP_RETCODE SCIPlpShrinkCols(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   newncols            /**< new number of columns in the LP */
   )
{
   SCIP_COL* col;
   int c;

   assert(lp != NULL);

   SCIPsetDebugMsg(set, "shrinking LP from %d to %d columns\n", lp->ncols, newncols);
   assert(0 <= newncols);
   assert(newncols <= lp->ncols);

   if( newncols < lp->ncols )
   {
      assert(!lp->diving);

      for( c = lp->ncols-1; c >= newncols; --c )
      {
         col = lp->cols[c];
         assert(col != NULL);
         assert(col->len == 0 || col->rows != NULL);
         assert(col->var != NULL);
         assert(SCIPvarGetStatus(col->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(col->var) == lp->cols[c]);
         assert(col->lppos == c);

         /* mark column to be removed from the LP */
         col->lppos = -1;
         col->lpdepth = -1;
         lp->ncols--;

         /* count removable columns */
         if( col->removable )
            lp->nremovablecols--;

         /* update column arrays of all linked rows */
         colUpdateDelLP(col, set);

         /* update the objective function vector norms */
         lpUpdateObjNorms(lp, set, col->unchangedobj, 0.0);
      }
      assert(lp->ncols == newncols);
      lp->lpifirstchgcol = MIN(lp->lpifirstchgcol, newncols);

      /* remove columns which are deleted from the lazy column array */
      c = 0;
      while( c < lp->nlazycols )
      {
         if( lp->lazycols[c]->lppos < 0 )
         {
            lp->lazycols[c] = lp->lazycols[lp->nlazycols-1];
            lp->nlazycols--;
         }
         else
            c++;
      }

      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      checkLazyColArray(lp, set);
      checkLinks(lp);
   }
   assert(lp->nremovablecols <= lp->ncols);

   return SCIP_OKAY;
}

/** removes and releases all rows after the given number of rows from the LP */
SCIP_RETCODE SCIPlpShrinkRows(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   int                   newnrows            /**< new number of rows in the LP */
   )
{
   SCIP_ROW* row;
   int r;

   assert(lp != NULL);
   assert(0 <= newnrows && newnrows <= lp->nrows);

   SCIPsetDebugMsg(set, "shrinking LP from %d to %d rows\n", lp->nrows, newnrows);
   if( newnrows < lp->nrows )
   {
      for( r = lp->nrows-1; r >= newnrows; --r )
      {
         row = lp->rows[r];
         assert(row != NULL);
         assert(row->len == 0 || row->cols != NULL);
         assert(row->lppos == r);

         /* mark row to be removed from the LP */
         row->lppos = -1;
         row->lpdepth = -1;
         lp->nrows--;

         /* count removable rows */
         if( row->removable )
            lp->nremovablerows--;

         /* update row arrays of all linked columns */
         rowUpdateDelLP(row);

         SCIProwUnlock(lp->rows[r]);

         /* check, if row deletion events are tracked
          * if so, issue ROWDELETEDLP event
          */
         if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDLP) != 0 )
         {
            SCIP_EVENT* event;

            SCIP_CALL( SCIPeventCreateRowDeletedLP(&event, blkmem, lp->rows[r]) );
            SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
         }

         SCIP_CALL( SCIProwRelease(&lp->rows[r], blkmem, set, lp) );
      }
      assert(lp->nrows == newnrows);
      lp->lpifirstchgrow = MIN(lp->lpifirstchgrow, newnrows);

      /* mark the current LP unflushed */
      lp->flushed = FALSE;

      checkLinks(lp);
   }
   assert(lp->nremovablerows <= lp->nrows);

   return SCIP_OKAY;
}

/** removes all columns and rows from LP, releases all rows */
SCIP_RETCODE SCIPlpClear(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   SCIPsetDebugMsg(set, "clearing LP\n");
   SCIP_CALL( SCIPlpShrinkCols(lp, set, 0) );
   SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, eventqueue, eventfilter, 0) );

   return SCIP_OKAY;
}

/** remembers number of columns and rows to track the newly added ones */
void SCIPlpMarkSize(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   lp->firstnewrow = lp->nrows;
   lp->firstnewcol = lp->ncols;
}

/** sets the remembered number of columns and rows to the given values */
void SCIPlpSetSizeMark(
   SCIP_LP*              lp,                 /**< current LP data */
   int                   nrows,              /**< number of rows to set the size marker to */
   int                   ncols               /**< number of columns to set the size marker to */
   )
{
   assert(lp != NULL);
   assert(!lp->diving);

   lp->firstnewrow = nrows;
   lp->firstnewcol = ncols;
}

/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
SCIP_RETCODE SCIPlpGetBasisInd(
   SCIP_LP*              lp,                 /**< LP data */
   int*                  basisind            /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);
   assert(basisind != NULL);

   SCIP_CALL( SCIPlpiGetBasisInd(lp->lpi, basisind) );

   return SCIP_OKAY;
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpGetBase(
   SCIP_LP*              lp,                 /**< LP data */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);

   SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

   return SCIP_OKAY;
}

/** gets a row from the inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpGetBInvRow(
   SCIP_LP*              lp,                 /**< LP data */
   int                   r,                  /**< row number */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                               *  (-1: if we do not store sparsity informations) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);
   assert(0 <= r && r < lp->nrows);  /* the basis matrix is nrows x nrows */
   assert(coef != NULL);

   SCIP_CALL( SCIPlpiGetBInvRow(lp->lpi, r, coef, inds, ninds) );

   return SCIP_OKAY;
}

/** gets a column from the inverse basis matrix B^-1 */
SCIP_RETCODE SCIPlpGetBInvCol(
   SCIP_LP*              lp,                 /**< LP data */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP
                                              *   returned by SCIPcolGetLPPos(); you have to call SCIPgetBasisInd()
                                              *   to get the array which links the B^-1 column numbers to the row and
                                              *   column numbers of the LP! c must be between 0 and nrows-1, since the
                                              *   basis has the size nrows * nrows */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                               *  (-1: if we do not store sparsity informations) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);
   assert(0 <= c && c < lp->nrows);  /* the basis matrix is nrows x nrows */
   assert(coef != NULL);

   SCIP_CALL( SCIPlpiGetBInvCol(lp->lpi, c, coef, inds, ninds) );

   return SCIP_OKAY;
}

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
SCIP_RETCODE SCIPlpGetBInvARow(
   SCIP_LP*              lp,                 /**< LP data */
   int                   r,                  /**< row number */
   SCIP_Real*            binvrow,            /**< row in B^-1 from prior call to SCIPlpGetBInvRow(), or NULL */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);
   assert(0 <= r && r < lp->nrows);  /* the basis matrix is nrows x nrows */
   assert(coef != NULL);

   SCIP_CALL( SCIPlpiGetBInvARow(lp->lpi, r, binvrow, coef, inds, ninds) );

   return SCIP_OKAY;
}

/** gets a column from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A),
 *  i.e., it computes B^-1 * A_c with A_c being the c'th column of A
 */
SCIP_RETCODE SCIPlpGetBInvACol(
   SCIP_LP*              lp,                 /**< LP data */
   int                   c,                  /**< column number which can be accessed by SCIPcolGetLPPos() */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *  (-1: if we do not store sparsity informations) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->solisbasic);
   assert(0 <= c && c < lp->ncols);
   assert(coef != NULL);

   SCIP_CALL( SCIPlpiGetBInvACol(lp->lpi, c, coef, inds, ninds) );

   return SCIP_OKAY;
}

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
SCIP_RETCODE SCIPlpSumRows(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   SCIP_Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   SCIP_Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   )
{
   SCIP_ROW* row;
   int r;
   int i;
   int idx;
   SCIP_Bool lhsinfinite;
   SCIP_Bool rhsinfinite;

   assert(lp != NULL);
   assert(prob != NULL);
   assert(weights != NULL);
   assert(sumcoef != NULL);
   assert(sumlhs != NULL);
   assert(sumrhs != NULL);

   /**@todo test, if a column based summation is faster */

   SCIP_CALL( SCIPrealarrayClear(sumcoef) );
   SCIP_CALL( SCIPrealarrayExtend(sumcoef, set->mem_arraygrowinit, set->mem_arraygrowfac, 0, prob->nvars-1) );
   *sumlhs = 0.0;
   *sumrhs = 0.0;
   lhsinfinite = FALSE;
   rhsinfinite = FALSE;
   for( r = 0; r < lp->nrows; ++r )
   {
      if( !SCIPsetIsZero(set, weights[r]) )
      {
         row = lp->rows[r];
         assert(row != NULL);
         assert(row->len == 0 || row->cols != NULL);
         assert(row->len == 0 || row->cols_index != NULL);
         assert(row->len == 0 || row->vals != NULL);

         /* add the row coefficients to the sum */
         for( i = 0; i < row->len; ++i )
         {
            assert(row->cols[i] != NULL);
            assert(row->cols[i]->var != NULL);
            assert(SCIPvarGetStatus(row->cols[i]->var) == SCIP_VARSTATUS_COLUMN);
            assert(SCIPvarGetCol(row->cols[i]->var) == row->cols[i]);
            assert(SCIPvarGetProbindex(row->cols[i]->var) == row->cols[i]->var_probindex);
            idx = row->cols[i]->var_probindex;
            assert(0 <= idx && idx < prob->nvars);
            SCIP_CALL( SCIPrealarrayIncVal(sumcoef, set->mem_arraygrowinit, set->mem_arraygrowfac, idx, weights[r] * row->vals[i]) );
         }

         /* add the row sides to the sum, depending on the sign of the weight */
         if( weights[r] > 0.0 )
         {
            lhsinfinite = lhsinfinite || SCIPsetIsInfinity(set, -row->lhs);
            if( !lhsinfinite )
               (*sumlhs) += weights[r] * (row->lhs - row->constant);
            rhsinfinite = rhsinfinite || SCIPsetIsInfinity(set, row->rhs);
            if( !rhsinfinite )
               (*sumrhs) += weights[r] * (row->rhs - row->constant);
         }
         else
         {
            lhsinfinite = lhsinfinite || SCIPsetIsInfinity(set, row->rhs);
            if( !lhsinfinite )
               (*sumlhs) += weights[r] * (row->rhs - row->constant);
            rhsinfinite = rhsinfinite || SCIPsetIsInfinity(set, -row->lhs);
            if( !rhsinfinite )
               (*sumrhs) += weights[r] * (row->lhs - row->constant);
         }
      }
   }

   if( lhsinfinite )
      *sumlhs = -SCIPsetInfinity(set);
   if( rhsinfinite )
      *sumrhs = SCIPsetInfinity(set);

   return SCIP_OKAY;
}

/** stores LP state (like basis information) into LP state object */
SCIP_RETCODE SCIPlpGetState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(blkmem != NULL);
   assert(lpistate != NULL);

   /* check whether there is no lp */
   if( lp->nlpicols == 0 && lp->nlpirows == 0 )
      *lpistate = NULL;
   else
   {
      SCIP_CALL( SCIPlpiGetState(lp->lpi, blkmem, lpistate) );
   }

   return SCIP_OKAY;
}

/** loads LP state (like basis information) into solver */
SCIP_RETCODE SCIPlpSetState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LPISTATE*        lpistate,           /**< LP state information (like basis information) */
   SCIP_Bool             wasprimfeas,        /**< primal feasibility when LP state information was stored */
   SCIP_Bool             wasprimchecked,     /**< true if the LP solution has passed the primal feasibility check */
   SCIP_Bool             wasdualfeas,        /**< dual feasibility when LP state information was stored */
   SCIP_Bool             wasdualchecked      /**< true if the LP solution has passed the dual feasibility check */
   )
{
   assert(lp != NULL);
   assert(blkmem != NULL);

   /* flush changes to the LP solver */
   SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );
   assert(lp->flushed);

   if( lp->solved && lp->solisbasic )
      return SCIP_OKAY;

   /* set LPI state in the LP solver */
   if( lpistate == NULL )
      lp->solisbasic = FALSE;
   else
   {
      SCIP_CALL( SCIPlpiSetState(lp->lpi, blkmem, lpistate) );
      lp->solisbasic = SCIPlpiHasStateBasis(lp->lpi, lpistate);
   }
   /* @todo: setting feasibility to TRUE might be wrong because in probing mode, the state is even saved when the LP was
    *        flushed and solved, also, e.g., when we hit the iteration limit
    */
   lp->primalfeasible = wasprimfeas;
   lp->primalchecked = wasprimchecked;
   lp->dualfeasible = wasdualfeas;
   lp->dualchecked = wasdualchecked;

   return SCIP_OKAY;
}

/** frees LP state information */
SCIP_RETCODE SCIPlpFreeState(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   )
{
   assert(lp != NULL);

   if( *lpistate != NULL )
   {
      SCIP_CALL( SCIPlpiFreeState(lp->lpi, blkmem, lpistate) );
   }

   return SCIP_OKAY;
}

/** stores pricing norms into LP norms object */
SCIP_RETCODE SCIPlpGetNorms(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LP pricing norms information */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(blkmem != NULL);
   assert(lpinorms != NULL);

   /* check whether there is no lp */
   if( lp->nlpicols == 0 && lp->nlpirows == 0 )
      *lpinorms = NULL;
   else
   {
      SCIP_CALL( SCIPlpiGetNorms(lp->lpi, blkmem, lpinorms) );
   }

   return SCIP_OKAY;
}

/** loads pricing norms from LP norms object into solver */
SCIP_RETCODE SCIPlpSetNorms(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS*        lpinorms            /**< LP pricing norms information */
   )
{
   assert(lp != NULL);
   assert(blkmem != NULL);
   assert(lp->flushed);

   /* set LPI norms in the LP solver */
   if( lpinorms != NULL )
   {
      SCIP_CALL( SCIPlpiSetNorms(lp->lpi, blkmem, lpinorms) );
   }

   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpFreeNorms(
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LP pricing norms information */
   )
{
   assert(lp != NULL);

   SCIP_CALL( SCIPlpiFreeNorms(lp->lpi, blkmem, lpinorms) );

   return SCIP_OKAY;
}

/** return the current cutoff bound of the lp */
SCIP_Real SCIPlpGetCutoffbound(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->cutoffbound;
}

/** sets the upper objective limit of the LP solver */
SCIP_RETCODE SCIPlpSetCutoffbound(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             cutoffbound         /**< new upper objective limit */
   )
{
   assert(lp != NULL);

   SCIPsetDebugMsg(set, "setting LP upper objective limit from %g to %g\n", lp->cutoffbound, cutoffbound);

   /* if the objective function was changed in diving, the cutoff bound has no meaning (it will be set correctly
    * in SCIPendDive())
    */
   if( SCIPlpDivingObjChanged(lp) )
   {
      assert(SCIPsetIsInfinity(set, lp->cutoffbound));
      return SCIP_OKAY;
   }

   /* if the cutoff bound is increased, and the LP was proved to exceed the old cutoff, it is no longer solved */
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT && cutoffbound > lp->cutoffbound )
   {
      /* mark the current solution invalid */
      lp->solved = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }
   /* if the cutoff bound is decreased below the current optimal value, the LP now exceeds the objective limit;
    * if the objective limit in the LP solver was disabled, the solution status of the LP is not changed
    */
   else if( !lpCutoffDisabled(set) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL
      && SCIPlpGetObjval(lp, set, prob) >= cutoffbound )
   {
      assert(lp->flushed);
      assert(lp->solved);
      lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
   }

   lp->cutoffbound = cutoffbound;

   return SCIP_OKAY;
}

/** returns the name of the given LP algorithm */
static
const char* lpalgoName(
   SCIP_LPALGO           lpalgo              /**< LP algorithm */
   )
{
   switch( lpalgo )
   {
   case SCIP_LPALGO_PRIMALSIMPLEX:
      return "primal simplex";
   case SCIP_LPALGO_DUALSIMPLEX:
      return "dual simplex";
   case SCIP_LPALGO_BARRIER:
      return "barrier";
   case SCIP_LPALGO_BARRIERCROSSOVER:
      return "barrier/crossover";
   default:
      SCIPerrorMessage("invalid LP algorithm\n");
      SCIPABORT();
      return "invalid"; /*lint !e527*/
   }
}

/** calls LPI to perform primal simplex, measures time and counts iterations, gets basis feasibility status */
static
SCIP_RETCODE lpPrimalSimplex(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_Real timedelta;
   SCIP_RETCODE retcode;
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   SCIPsetDebugMsg(set, "solving LP %" SCIP_LONGINT_FORMAT " (%d cols, %d rows) with primal simplex (diving=%d, nprimallps=%" SCIP_LONGINT_FORMAT ", ndivinglps=%" SCIP_LONGINT_FORMAT ")\n",
      stat->lpcount+1, lp->ncols, lp->nrows, lp->diving || lp->probing, stat->nprimallps, stat->ndivinglps);

   *lperror = FALSE;

#if 0 /* for debugging: write all root node LP's */
   if( stat->nnodes == 1 && !lp->diving && !lp->probing )
   {
      char fname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(fname, SCIP_MAXSTRLEN, "lp%" SCIP_LONGINT_FORMAT "_%" SCIP_LONGINT_FORMAT ".lp", stat->nnodes, stat->lpcount);
      SCIP_CALL( SCIPlpWrite(lp, fname) );
      SCIPsetDebugMsg("wrote LP to file <%s> (primal simplex, objlim=%.15g, feastol=%.15g/%.15g, fromscratch=%d, fastmip=%d, scaling=%d, presolving=%d)\n",
         fname, lp->lpiobjlim, lp->lpifeastol, lp->lpidualfeastol,
         lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling, lp->lpipresolving);
   }
#endif

   /* start timing */
   if( lp->diving || lp->probing )
   {
      if( lp->strongbranchprobing )
         SCIPclockStart(stat->strongbranchtime, set);
      else
         SCIPclockStart(stat->divinglptime, set);

      timedelta = 0.0;   /* unused for diving or probing */
   }
   else
   {
      SCIPclockStart(stat->primallptime, set);
      timedelta = -SCIPclockGetTime(stat->primallptime);
   }

   /* call primal simplex */
   retcode = SCIPlpiSolvePrimal(lp->lpi);
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      SCIPsetDebugMsg(set, "(node %" SCIP_LONGINT_FORMAT ") primal simplex solving error in LP %" SCIP_LONGINT_FORMAT "\n", stat->nnodes, stat->nlps);
   }
   else
   {
      SCIP_CALL( retcode );
   }
   lp->lastlpalgo = SCIP_LPALGO_PRIMALSIMPLEX;
   lp->solisbasic = TRUE;

   /* stop timing */
   if( lp->diving || lp->probing )
   {
      if( lp->strongbranchprobing )
         SCIPclockStop(stat->strongbranchtime, set);
      else
         SCIPclockStop(stat->divinglptime, set);
   }
   else
   {
      timedelta += SCIPclockGetTime(stat->primallptime);
      SCIPclockStop(stat->primallptime, set);
   }

   /* count number of iterations */
   SCIPstatIncrement(stat, set, lpcount);
   SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      if( !lp->strongbranchprobing )
      {
         SCIPstatIncrement(stat, set, nlps);
         SCIPstatAdd( stat, set, nlpiterations, iterations );
      }
      if( resolve && !lp->lpifromscratch && stat->nlps > 1 )
      {
         SCIPstatIncrement(stat, set, nprimalresolvelps );
         SCIPstatAdd(stat, set, nprimalresolvelpiterations, iterations);
      }
      if( lp->diving || lp->probing )
      {
         if( lp->strongbranchprobing )
         {
            SCIPstatIncrement(stat, set, nsbdivinglps);
            SCIPstatAdd(stat, set, nsbdivinglpiterations, iterations);
         }
         else
         {
            SCIPstatUpdate(stat, set, lastdivenode, stat->nnodes);
            SCIPstatIncrement(stat, set, ndivinglps);
            SCIPstatAdd(stat, set, ndivinglpiterations, iterations);
         }
      }
      else
      {
         SCIPstatIncrement(stat, set, nprimallps);
         SCIPstatAdd(stat, set, nprimallpiterations, iterations);
      }
   }
   else
   {
      if ( ! lp->diving && ! lp->probing )
      {
         SCIPstatIncrement(stat, set, nprimalzeroitlps);
         SCIPstatAdd(stat, set, primalzeroittime, timedelta);
      }

      if ( keepsol && !(*lperror) )
      {
         /* the solution didn't change: if the solution was valid before resolve, it is still valid */
         if( lp->validsollp == stat->lpcount-1 )
            lp->validsollp = stat->lpcount;
         if( lp->validfarkaslp == stat->lpcount-1 )
            lp->validfarkaslp = stat->lpcount;
      }
   }

   SCIPsetDebugMsg(set, "solved LP %" SCIP_LONGINT_FORMAT " with primal simplex (diving=%d, nprimallps=%" SCIP_LONGINT_FORMAT ", ndivinglps=%" SCIP_LONGINT_FORMAT ")\n",
      stat->lpcount, lp->diving || lp->probing, stat->nprimallps, stat->ndivinglps);

   return SCIP_OKAY;
}

/** calls LPI to perform dual simplex, measures time and counts iterations */
static
SCIP_RETCODE lpDualSimplex(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_Real timedelta;
   SCIP_RETCODE retcode;
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   SCIPsetDebugMsg(set, "solving LP %" SCIP_LONGINT_FORMAT " (%d cols, %d rows) with dual simplex (diving=%d, nduallps=%" SCIP_LONGINT_FORMAT ", ndivinglps=%" SCIP_LONGINT_FORMAT ")\n",
      stat->lpcount+1, lp->ncols, lp->nrows, lp->diving || lp->probing, stat->nduallps, stat->ndivinglps);

   *lperror = FALSE;

#if 0 /* for debugging: write all root node LP's */
   if( stat->nnodes == 1 && !lp->diving && !lp->probing )
   {
      char fname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(fname, SCIP_MAXSTRLEN, "lp%" SCIP_LONGINT_FORMAT "_%" SCIP_LONGINT_FORMAT ".lp", stat->nnodes, stat->lpcount);
      SCIP_CALL( SCIPlpWrite(lp, fname) );
      SCIPsetDebugMsg("wrote LP to file <%s> (dual simplex, objlim=%.15g, feastol=%.15g/%.15g, fromscratch=%d, fastmip=%d, scaling=%d, presolving=%d)\n",
         fname, lp->lpiobjlim, lp->lpifeastol, lp->lpidualfeastol,
         lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling, lp->lpipresolving);
   }
#endif

   /* start timing */
   if( lp->diving || lp->probing )
   {
      if( lp->strongbranchprobing )
         SCIPclockStart(stat->strongbranchtime, set);
      else
         SCIPclockStart(stat->divinglptime, set);

      timedelta = 0.0;   /* unused for diving or probing */
   }
   else
   {
      SCIPclockStart(stat->duallptime, set);
      timedelta = -SCIPclockGetTime(stat->duallptime);      
   }

   /* call dual simplex */
   retcode = SCIPlpiSolveDual(lp->lpi);
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      SCIPsetDebugMsg(set, "(node %" SCIP_LONGINT_FORMAT ") dual simplex solving error in LP %" SCIP_LONGINT_FORMAT "\n", stat->nnodes, stat->nlps);
   }
   else
   {
      SCIP_CALL( retcode );
   }
   lp->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;
   lp->solisbasic = TRUE;

   /* stop timing */
   if( lp->diving || lp->probing )
   {
      if( lp->strongbranchprobing )
         SCIPclockStop(stat->strongbranchtime, set);
      else
         SCIPclockStop(stat->divinglptime, set);
   }
   else
   {
      timedelta += SCIPclockGetTime(stat->duallptime);
      SCIPclockStop(stat->duallptime, set);
   }

   /* count number of iterations */
   SCIPstatIncrement(stat, set, lpcount);
   SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      if( !lp->strongbranchprobing )
      {
         SCIPstatIncrement(stat, set, nlps);
         SCIPstatAdd(stat, set, nlpiterations, iterations);
      }
      if( resolve && !lp->lpifromscratch && stat->nlps > 1  )
      {
         SCIPstatIncrement(stat, set, ndualresolvelps);
         SCIPstatAdd(stat, set, ndualresolvelpiterations, iterations);
      }
      if( lp->diving || lp->probing )
      {
         if( lp->strongbranchprobing )
         {
            SCIPstatIncrement(stat, set, nsbdivinglps);
            SCIPstatAdd(stat, set, nsbdivinglpiterations, iterations);
         }
         else
         {
            SCIPstatUpdate(stat, set, lastdivenode, stat->nnodes);
            SCIPstatIncrement(stat, set, ndivinglps);
            SCIPstatAdd(stat, set, ndivinglpiterations, iterations);
         }
      }
      else
      {
         SCIPstatIncrement(stat, set, nduallps);
         SCIPstatAdd(stat, set, nduallpiterations, iterations);
      }
   }
   else
   {
      if ( ! lp->diving && ! lp->probing )
      {
         SCIPstatIncrement(stat, set, ndualzeroitlps);
         SCIPstatAdd(stat, set, dualzeroittime, timedelta);
      }

      if( keepsol && !(*lperror) )
      {
         /* the solution didn't change: if the solution was valid before resolve, it is still valid */
         if( lp->validsollp == stat->lpcount-1 )
            lp->validsollp = stat->lpcount;
         if( lp->validfarkaslp == stat->lpcount-1 )
            lp->validfarkaslp = stat->lpcount;
      }
   }

   SCIPsetDebugMsg(set, "solved LP %" SCIP_LONGINT_FORMAT " with dual simplex (diving=%d, nduallps=%" SCIP_LONGINT_FORMAT ", ndivinglps=%" SCIP_LONGINT_FORMAT ")\n",
      stat->lpcount, lp->diving || lp->probing, stat->nduallps, stat->ndivinglps);

   return SCIP_OKAY;
}

/** calls LPI to perform lexicographic dual simplex to find a lexicographically minimal optimal solution, measures time and counts iterations
 *
 *  We follow the approach of the following paper to find a lexicographically minimal optimal
 *  solution:
 *
 *  Zanette, Fischetti, Balas@n
 *  Can pure cutting plane algorithms work?@n
 *  IPCO 2008, Bertinoro, Italy.
 *
 *  We do, however, not aim for the exact lexicographically minimal optimal solutions, but perform a
 *  heuristic, i.e., we limit the number of components which are minimized.
 *
 *  More precisely, we first solve the problem with the dual simplex algorithm. Then we fix those
 *  nonbasic variables to their current value (i.e., one of the bounds except maybe for free
 *  variables) that have nonzero reduced cost. This fixes the objective function value, because only
 *  pivots that will not change the objective are allowed afterwards.
 *
 *  Then the not yet fixed variables are considered in turn. If they are at their lower bounds and
 *  nonbasic, they are fixed to this bound, since their value cannot be decreased further. Once a
 *  candidate is found, we set the objective to minimize this variable. We run the primal simplex
 *  algorithm (since the objective is changed the solution is not dual feasible anymore; if
 *  variables out of the basis have been fixed to their lower bound, the basis is also not primal
 *  feasible anymore). After the optimization, we again fix nonbasic variables that have nonzero
 *  reduced cost. We then choose the next variable and iterate.
 *
 *  We stop the process once we do not find candidates or have performed a maximum number of
 *  iterations.
 *
 *  @todo Does this really produce a lexicographically minimal solution? 
 *  @todo Can we skip the consideration of basic variables that are at their lower bound? How can we
 *    guarantee that these variables will not be changed in later stages? We can fix these variables
 *    to their lower bound, but this destroys the basis.
 *  @todo Should we use lexicographical minimization in diving/probing or not?
 */
static
SCIP_RETCODE lpLexDualSimplex(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_Real timedelta;
   SCIP_RETCODE retcode;
   int totalIterations;
   int lexIterations;
   int iterations;
   int rounds;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   SCIPsetDebugMsg(set, "solving LP %" SCIP_LONGINT_FORMAT " (%d cols, %d rows) with lex dual simplex (diving=%d, nduallps=%" SCIP_LONGINT_FORMAT ", ndivinglps=%" SCIP_LONGINT_FORMAT ")\n",
      stat->lpcount+1, lp->ncols, lp->nrows, lp->diving || lp->probing, stat->nduallps, stat->ndivinglps);

   *lperror = FALSE;

   /* start timing */
   if( lp->diving || lp->probing )
   {
      if( lp->strongbranchprobing )
         SCIPclockStart(stat->strongbranchtime, set);
      else
         SCIPclockStart(stat->divinglptime, set);

      timedelta = 0.0;   /* unused for diving or probing */
   }
   else
   {
      SCIPclockStart(stat->duallptime, set);
      timedelta = -SCIPclockGetTime(stat->duallptime);      
   }

   /* call dual simplex for first lp */
   retcode = SCIPlpiSolveDual(lp->lpi);
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      SCIPsetDebugMsg(set, "(node %" SCIP_LONGINT_FORMAT ") dual simplex solving error in LP %" SCIP_LONGINT_FORMAT "\n", stat->nnodes, stat->nlps);
   }
   else
   {
      SCIP_CALL( retcode );
   }
   SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
   totalIterations = iterations;

   /* stop timing */
   if( lp->diving || lp->probing )
   {
      if( lp->strongbranchprobing )
         SCIPclockStop(stat->strongbranchtime, set);
      else
         SCIPclockStop(stat->divinglptime, set);
   }
   else
   {
      timedelta += SCIPclockGetTime(stat->duallptime);
      SCIPclockStop(stat->duallptime, set);
   }

   /* count number of iterations */
   SCIPstatIncrement(stat, set, lpcount);
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      if( lp->strongbranchprobing )
      {
         SCIPstatAdd(stat, set, nlpiterations, iterations);
      }
      if( resolve && !lp->lpifromscratch && stat->nlps > 1  )
      {
         SCIPstatIncrement(stat, set, ndualresolvelps);
         SCIPstatAdd(stat, set, ndualresolvelpiterations, iterations);
      }
      if( lp->diving || lp->probing )
      {
         if( lp->strongbranchprobing )
         {
            SCIPstatIncrement(stat, set, nsbdivinglps);
            SCIPstatAdd(stat, set, nsbdivinglpiterations, iterations);
         }
         else
         {
            SCIPstatUpdate(stat, set, lastdivenode, stat->nnodes);
            SCIPstatIncrement(stat, set, ndivinglps);
            SCIPstatAdd(stat, set, ndivinglpiterations, iterations);
         }
      }
      else
      {
         SCIPstatIncrement(stat, set, nduallps);
         SCIPstatAdd(stat, set, nduallpiterations, iterations);
      }
   }
   else
   {
      if ( ! lp->diving && ! lp->probing )
      {
         SCIPstatIncrement(stat, set, ndualzeroitlps);
         SCIPstatAdd(stat, set, dualzeroittime, timedelta);
      }
   }
   lexIterations = 0;

   /* search for lexicographically minimal optimal solution */
   if( !lp->diving && !lp->probing && SCIPlpiIsOptimal(lp->lpi) )
   {
      SCIP_Bool chooseBasic;
      SCIP_Real* primsol;
      SCIP_Real* dualsol;
      SCIP_Real* redcost;
      int* cstat;
      int* rstat;
      SCIP_Real* newobj;
      SCIP_Real* newlb;
      SCIP_Real* newub;
      SCIP_Real* newlhs;
      SCIP_Real* newrhs;
      SCIP_Real* oldlb;
      SCIP_Real* oldub;
      SCIP_Real* oldlhs;
      SCIP_Real* oldrhs;
      SCIP_Real* oldobj;
      SCIP_Bool* fixedc;
      SCIP_Bool* fixedr;
      int* indcol;
      int* indrow;
      int* indallcol;
      int* indallrow;
      int nDualDeg;
      int r, c;
      int cntcol;
      int cntrow;
      int nruns;
      int pos;

      chooseBasic = set->lp_lexdualbasic;

      /* start timing */
      SCIPclockStart(stat->lexduallptime, set);

      /* get all solution information */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &dualsol, lp->nlpirows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &redcost, lp->nlpicols) );
      if( chooseBasic )
      {
         SCIP_CALL( SCIPsetAllocBufferArray(set, &primsol, lp->nlpicols) );
      }
      else
         primsol = NULL;

      /* get basic and nonbasic information */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, lp->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, lp->nlpirows) );

      /* save bounds, lhs/rhs, and objective */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &oldobj, lp->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &oldlb, lp->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &oldub, lp->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &oldlhs, lp->nlpirows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &oldrhs, lp->nlpirows) );
      SCIP_CALL( SCIPlpiGetBounds(lp->lpi, 0, lp->nlpicols-1, oldlb, oldub) );
      SCIP_CALL( SCIPlpiGetSides(lp->lpi, 0, lp->nlpirows-1, oldlhs, oldrhs) );
      SCIP_CALL( SCIPlpiGetObj(lp->lpi, 0, lp->nlpicols-1, oldobj) );

      /* get storage for several arrays */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &newlb, lp->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &newub, lp->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &indcol, lp->nlpicols) );

      SCIP_CALL( SCIPsetAllocBufferArray(set, &newlhs, lp->nlpirows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &newrhs, lp->nlpirows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &indrow, lp->nlpirows) );

      SCIP_CALL( SCIPsetAllocBufferArray(set, &indallcol, lp->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &indallrow, lp->nlpirows) );

      SCIP_CALL( SCIPsetAllocBufferArray(set, &fixedc, lp->nlpicols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &fixedr, lp->nlpirows) );

      /* initialize: set objective to 0, get fixed variables */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &newobj, lp->nlpicols) );
      for( c = 0; c < lp->nlpicols; ++c )
      {
         newobj[c] = 0.0;
         indallcol[c] = c;
         if( SCIPsetIsFeasEQ(set, oldlb[c], oldub[c]) )
            fixedc[c] = TRUE;
         else
            fixedc[c] = FALSE;
      }

      /* initialize: get fixed slack variables */
      for( r = 0; r < lp->nlpirows; ++r )
      {
         indallrow[r] = r;
         if( SCIPsetIsFeasEQ(set, oldlhs[r], oldrhs[r]) )
            fixedr[r] = TRUE;
         else
            fixedr[r] = FALSE;
      }

#ifdef DEBUG_LEXDUAL
      {
         int j;

         if( !chooseBasic )
         {
            assert(primsol == NULL);
            SCIP_CALL( SCIPsetAllocBufferArray(set, &primsol, lp->nlpicols) );
         }
         assert(primsol != NULL);
         SCIP_CALL( SCIPlpiGetSol(lp->lpi, NULL, primsol, NULL, NULL, NULL) );
         SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

         for( j = 0; j < lp->nlpicols; ++j )
         {
            if( fixedc[j] )
            {
               SCIPsetDebugMsg(set, "%f (%d) [f] ", primsol[j], j);
            }
            else
            {
               char type;
               switch( (SCIP_BASESTAT) cstat[j] )
               {
               case SCIP_BASESTAT_LOWER:
                  type = 'l'; 
                  break;
               case SCIP_BASESTAT_UPPER:
                  type = 'u';
                  break;
               case SCIP_BASESTAT_ZERO:
                  type = 'z';
                  break;
               case SCIP_BASESTAT_BASIC:
                  type = 'b';
                  break;
               default:
                  type = '?';
                  SCIPerrorMessage("unknown base stat %d\n", cstat[j]);
                  SCIPABORT();
               }
               SCIPsetDebugMsg(set, "%f (%d) [%c] ", primsol[j], j, type);
            }
         }
         SCIPsetDebugMsg(set, "\n\n");

         if( !chooseBasic )
         {
            SCIPsetFreeBufferArray(set, &primsol);
            assert(primsol == NULL);
         }
      }
#endif

      /* perform lexicographic rounds */
      pos = -1;
      nruns = 0;
      rounds = 0;
      /* SCIP_CALL( lpSetLPInfo(lp, TRUE) ); */
      do
      {
         int oldpos;

         /* get current solution */
         if( chooseBasic )
            SCIP_CALL( SCIPlpiGetSol(lp->lpi, NULL, primsol, dualsol, NULL, redcost) );
         else
         {
            SCIP_CALL( SCIPlpiGetSol(lp->lpi, NULL, NULL, dualsol, NULL, redcost) );
            assert(primsol == NULL);
         }

         /* get current basis */
         SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

         /* check columns: find first candidate (either basic or nonbasic and zero reduced cost) and fix variables */
         nDualDeg = 0;
         cntcol = 0;
         oldpos = pos;
         pos = -1;
         for( c = 0; c < lp->nlpicols; ++c )
         {
            if( !fixedc[c] )
            {
               /* check whether variable is in basis */
               if( (SCIP_BASESTAT) cstat[c] == SCIP_BASESTAT_BASIC )
               {
                  /* store first candidate */
                  if( pos == -1 && c > oldpos )
                  {
                     if( !chooseBasic || !SCIPsetIsIntegral(set, primsol[c]) ) /*lint !e613*/
                        pos = c;
                  }
               }
               else
               {
                  /* reduced cost == 0 -> possible candidate */
                  if( SCIPsetIsDualfeasZero(set, redcost[c]) )
                  {
                     ++nDualDeg;
                     /* only if we have not yet found a candidate */
                     if( pos == -1 && c > oldpos )
                     {
                        /* if the variable is at its lower bound - fix it, because its value cannot be reduced */
                        if( (SCIP_BASESTAT) cstat[c] == SCIP_BASESTAT_LOWER )
                        {
                           newlb[cntcol] = oldlb[c];
                           newub[cntcol] = oldlb[c];
                           indcol[cntcol++] = c;
                           fixedc[c] = TRUE;
                        }
                        else /* found a non-fixed candidate */
                        {
                           if( !chooseBasic )
                              pos = c;
                        }
                     }
                  }
                  else
                  {
                     /* nonzero reduced cost -> variable can be fixed */
                     if( (SCIP_BASESTAT) cstat[c] == SCIP_BASESTAT_LOWER )
                     {
                        newlb[cntcol] = oldlb[c];
                        newub[cntcol] = oldlb[c];
                     }
                     else
                     {
                        if( (SCIP_BASESTAT) cstat[c] == SCIP_BASESTAT_UPPER )
                        {
                           newlb[cntcol] = oldub[c];
                           newub[cntcol] = oldub[c];
                        }
                        else
                        {
                           assert((SCIP_BASESTAT) cstat[c] == SCIP_BASESTAT_ZERO);
                           newlb[cntcol] = 0.0;
                           newub[cntcol] = 0.0;
                        }
                     }
                     indcol[cntcol++] = c;
                     fixedc[c] = TRUE;
                  }
               }
            }
         }

         /* check rows */
         cntrow = 0;
         for( r = 0; r < lp->nlpirows; ++r )
         {
            if( !fixedr[r] )
            {
               /* consider only nonbasic rows */
               if( (SCIP_BASESTAT) rstat[r] != SCIP_BASESTAT_BASIC )
               {
                  assert((SCIP_BASESTAT) rstat[r] != SCIP_BASESTAT_ZERO);
                  if( SCIPsetIsFeasZero(set, dualsol[r]) )
                     ++nDualDeg;
                  else
                  {
                     if( SCIPsetIsFeasPositive(set, dualsol[r]) ) 
                     {
                        assert(!SCIPsetIsInfinity(set, -oldlhs[r]));
                        newlhs[cntrow] = oldlhs[r];
                        newrhs[cntrow] = oldlhs[r];
                     }
                     else
                     {
                        assert(!SCIPsetIsInfinity(set, oldrhs[r]));
                        newlhs[cntrow] = oldrhs[r];
                        newrhs[cntrow] = oldrhs[r];
                     }
                     indrow[cntrow++] = r;
                     fixedr[r] = TRUE;
                  }
               }
            }
         }

         if( nDualDeg > 0 && pos >= 0 )
         {
            assert(0 <= pos && pos < lp->nlpicols && pos > oldpos);

            /* change objective */
            if( nruns == 0 )
            {
               /* set objective to appropriate unit vector for first run */
               newobj[pos] = 1.0;
               SCIP_CALL( SCIPlpiChgObj(lp->lpi, lp->nlpicols, indallcol, newobj) );
            }
            else
            {
               /* set obj. coef. to 1 for other runs (ones remain in previous positions) */
               SCIP_Real obj = 1.0;
               SCIP_CALL( SCIPlpiChgObj(lp->lpi, 1, &pos, &obj) );
            }

            /* fix variables */
            SCIP_CALL( SCIPlpiChgBounds(lp->lpi, cntcol, indcol, newlb, newub) );
            SCIP_CALL( SCIPlpiChgSides(lp->lpi, cntrow, indrow, newlhs, newrhs) );

            /* solve with primal simplex, because we are primal feasible, but not necessarily dual feasible */
            retcode = SCIPlpiSolvePrimal(lp->lpi);
            if( retcode == SCIP_LPERROR )
            {
               *lperror = TRUE;
               SCIPsetDebugMsg(set, "(node %" SCIP_LONGINT_FORMAT ") in lex-dual: primal simplex solving error in LP %" SCIP_LONGINT_FORMAT "\n", stat->nnodes, stat->nlps);
            }
            else
            {
               SCIP_CALL( retcode );
            }
            SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
            lexIterations += iterations;

#ifdef DEBUG_LEXDUAL
            if( iterations > 0 )
            {
               int j;

               if( !chooseBasic )
               {
                  assert(primsol == NULL);
                  SCIP_CALL( SCIPsetAllocBufferArray(set, &primsol, lp->nlpicols) );
               }
               assert(primsol != NULL);
               SCIP_CALL( SCIPlpiGetSol(lp->lpi, NULL, primsol, NULL, NULL, NULL) );

               for( j = 0; j < lp->nlpicols; ++j )
               {
                  if( fixedc[j] )
                  {
                     SCIPsetDebugMsg(set, "%f (%d) [f] ", primsol[j], j);
                  }
                  else
                  {
                     char cstart = '[';
                     char cend = ']';
                     char type;

                     if(j == pos)
                     {
                        cstart = '*';
                        cend = '*';
                     }

                     switch( (SCIP_BASESTAT) cstat[j] )
                     {
                     case SCIP_BASESTAT_LOWER:
                        type = 'l'; 
                        break;
                     case SCIP_BASESTAT_UPPER:
                        type = 'u';
                        break;
                     case SCIP_BASESTAT_ZERO:
                        type = 'z'; 
                        break;
                     case SCIP_BASESTAT_BASIC:
                        type = 'b'; 
                        break;
                     default: 
                        type = '?';
                        SCIPerrorMessage("unknown base state %d\n", cstat[j]);
                        SCIPABORT();
                     }
                     SCIPsetDebugMsg(set, "%f (%d) %c%c%c ", primsol[j], j, cstart, type, cend);
                  }
               }
               SCIPsetDebugMsg(set, "\n\n");

               if( !chooseBasic )
               {
                  SCIPsetFreeBufferArray(set, &primsol);
                  assert(primsol == NULL);
               }
            }
#endif

            /* count only as round if iterations have been performed */
            if( iterations > 0 )
               ++rounds;
            ++nruns;
         }
      }
      while( pos >= 0 && nDualDeg > 0 && (set->lp_lexdualmaxrounds == -1 || rounds < set->lp_lexdualmaxrounds) );

      /* reset bounds, lhs/rhs, and obj */
      SCIP_CALL( SCIPlpiChgBounds(lp->lpi, lp->nlpicols, indallcol, oldlb, oldub) );
      SCIP_CALL( SCIPlpiChgSides(lp->lpi, lp->nlpirows, indallrow, oldlhs, oldrhs) );
      SCIP_CALL( SCIPlpiChgObj(lp->lpi, lp->nlpicols, indallcol, oldobj) );

      /* resolve to update solvers internal data structures - should only produce few pivots - is this needed? */
      retcode = SCIPlpiSolveDual(lp->lpi);
      if( retcode == SCIP_LPERROR )
      {
         *lperror = TRUE;
         SCIPsetDebugMsg(set, "(node %" SCIP_LONGINT_FORMAT ") dual simplex solving error in LP %" SCIP_LONGINT_FORMAT "\n", stat->nnodes, stat->nlps);
      }
      else
      {
         SCIP_CALL( retcode );
      }
      assert(SCIPlpiIsOptimal(lp->lpi));
      SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
      lexIterations += iterations;

      /* SCIP_CALL( lpSetLPInfo(lp, set->disp_lpinfo) ); */

      /* count number of iterations */
      if( totalIterations == 0 && lexIterations > 0 && !lp->strongbranchprobing )
         SCIPstatIncrement(stat, set, nlps);

      if( lexIterations > 0 ) /* don't count the resolves after removing unused columns/rows */
      {
         SCIPstatAdd(stat, set, nlpiterations, lexIterations);
         if( resolve && !lp->lpifromscratch && stat->nlps > 1  )
         {
            SCIPstatIncrement(stat, set, nlexdualresolvelps);
            SCIPstatAdd(stat, set, nlexdualresolvelpiterations, lexIterations);
         }
         SCIPstatIncrement(stat, set, nlexduallps);
         SCIPstatAdd(stat, set, nlexduallpiterations, lexIterations);

         totalIterations += lexIterations;
      }

      /* free space */
      SCIPsetFreeBufferArray(set, &newobj);

      SCIPsetFreeBufferArray(set, &fixedr);
      SCIPsetFreeBufferArray(set, &fixedc);

      SCIPsetFreeBufferArray(set, &indallrow);
      SCIPsetFreeBufferArray(set, &indallcol);

      SCIPsetFreeBufferArray(set, &indrow);
      SCIPsetFreeBufferArray(set, &newrhs);
      SCIPsetFreeBufferArray(set, &newlhs);

      SCIPsetFreeBufferArray(set, &indcol);
      SCIPsetFreeBufferArray(set, &newub);
      SCIPsetFreeBufferArray(set, &newlb);

      SCIPsetFreeBufferArray(set, &oldobj);
      SCIPsetFreeBufferArray(set, &oldrhs);
      SCIPsetFreeBufferArray(set, &oldlhs);
      SCIPsetFreeBufferArray(set, &oldub);
      SCIPsetFreeBufferArray(set, &oldlb);

      SCIPsetFreeBufferArray(set, &rstat);
      SCIPsetFreeBufferArray(set, &cstat);

      SCIPsetFreeBufferArray(set, &redcost);
      SCIPsetFreeBufferArray(set, &dualsol);
      if( chooseBasic )
         SCIPsetFreeBufferArray(set, &primsol);

      /* stop timing */
      SCIPclockStop(stat->lexduallptime, set);

      SCIPsetDebugMsg(set, "solved LP %" SCIP_LONGINT_FORMAT " with lex dual simplex (diving=%d, nduallps=%" SCIP_LONGINT_FORMAT ", ndivinglps=%" SCIP_LONGINT_FORMAT ")\n",
         stat->lpcount, lp->diving || lp->probing, stat->nduallps, stat->ndivinglps);
   }
   lp->lastlpalgo = SCIP_LPALGO_DUALSIMPLEX;
   lp->solisbasic = TRUE;

   if( totalIterations > 0 && !lp->strongbranchprobing )
      SCIPstatIncrement(stat, set, nlps);
   else
   {
      if( keepsol && !(*lperror) )
      {
         /* the solution didn't change: if the solution was valid before resolve, it is still valid */
         if( lp->validsollp == stat->lpcount-1 )
            lp->validsollp = stat->lpcount;
         if( lp->validfarkaslp == stat->lpcount-1 )
            lp->validfarkaslp = stat->lpcount;
      }
   }

   return SCIP_OKAY;
}

/** calls LPI to perform barrier, measures time and counts iterations, gets basis feasibility status */
static
SCIP_RETCODE lpBarrier(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool             crossover,          /**< should crossover be performed? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_Real timedelta;
   SCIP_RETCODE retcode;
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   SCIPsetDebugMsg(set, "solving LP %" SCIP_LONGINT_FORMAT " (%d cols, %d rows) with barrier%s (diving=%d, nbarrierlps=%" SCIP_LONGINT_FORMAT ", ndivinglps=%" SCIP_LONGINT_FORMAT ")\n",
      stat->lpcount+1, lp->ncols, lp->nrows, crossover ? "/crossover" : "", lp->diving || lp->probing,
      stat->nbarrierlps, stat->ndivinglps);

   *lperror = FALSE;

#if 0 /* for debugging: write all root node LP's */
   if( stat->nnodes == 1 && !lp->diving && !lp->probing )
   {
      char fname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(fname, SCIP_MAXSTRLEN, "lp%" SCIP_LONGINT_FORMAT "_%" SCIP_LONGINT_FORMAT ".lp", stat->nnodes, stat->lpcount);
      SCIP_CALL( SCIPlpWrite(lp, fname) );
      SCIPsetDebugMsg("wrote LP to file <%s> (barrier, objlim=%.15g, feastol=%.15g/%.15g, convtol=%.15g, fromscratch=%d, fastmip=%d, scaling=%d, presolving=%d)\n",
         fname, lp->lpiobjlim, lp->lpifeastol, lp->lpidualfeastol, lp->lpibarrierconvtol,
         lp->lpifromscratch, lp->lpifastmip, lp->lpiscaling, lp->lpipresolving);
   }
#endif

   /* start timing */
   if( lp->diving || lp->probing )
   {
      if( lp->strongbranchprobing )
         SCIPclockStart(stat->strongbranchtime, set);
      else
         SCIPclockStart(stat->divinglptime, set);

      timedelta = 0.0;   /* unused for diving or probing */
   }
   else
   {
      SCIPclockStart(stat->barrierlptime, set);
      timedelta = -SCIPclockGetTime(stat->duallptime);      
   }

   /* call barrier algorithm */
   retcode = SCIPlpiSolveBarrier(lp->lpi, crossover);
   if( retcode == SCIP_LPERROR )
   {
      *lperror = TRUE;
      SCIPsetDebugMsg(set, "(node %" SCIP_LONGINT_FORMAT ") barrier solving error in LP %" SCIP_LONGINT_FORMAT "\n", stat->nnodes, stat->nlps);
   }
   else
   {
      SCIP_CALL( retcode );
   }
   lp->lastlpalgo = (crossover ? SCIP_LPALGO_BARRIERCROSSOVER : SCIP_LPALGO_BARRIER);
   lp->solisbasic = crossover;

   /* stop timing */
   if( lp->diving || lp->probing )
   {
      if( lp->strongbranchprobing )
         SCIPclockStop(stat->strongbranchtime, set);
      else
         SCIPclockStop(stat->divinglptime, set);
   }
   else
   {
      SCIPclockStop(stat->barrierlptime, set);
      timedelta = -SCIPclockGetTime(stat->duallptime);      
   }

   /* count number of iterations */
   SCIPstatIncrement(stat, set, lpcount);
   SCIP_CALL( SCIPlpGetIterations(lp, &iterations) );
   if( iterations > 0 ) /* don't count the resolves after removing unused columns/rows */
   {
      if( !lp->strongbranchprobing )
      {
         SCIPstatIncrement(stat, set, nlps);
         SCIPstatAdd(stat, set, nlpiterations, iterations);
      }
      if( lp->diving || lp->probing )
      {
         if( lp->strongbranchprobing )
         {
            SCIPstatIncrement(stat, set, nsbdivinglps);
            SCIPstatAdd(stat, set, nsbdivinglpiterations, iterations);
         }
         else
         {
            SCIPstatUpdate(stat, set, lastdivenode, stat->nnodes);
            SCIPstatIncrement(stat, set, ndivinglps);
            SCIPstatAdd(stat, set, ndivinglpiterations, iterations);
         }
      }
      else
      {
         SCIPstatIncrement(stat, set, nbarrierlps);
         SCIPstatAdd(stat, set, nbarrierlpiterations, iterations);
      }
   }
   else
   {
      if ( ! lp->diving && ! lp->probing )
      {
         SCIPstatIncrement(stat, set, nbarrierzeroitlps);
         SCIPstatAdd(stat, set, barrierzeroittime, timedelta);
      }

      if( keepsol && !(*lperror) )
      {
         /* the solution didn't change: if the solution was valid before resolve, it is still valid */
         if( lp->validsollp == stat->lpcount-1 )
            lp->validsollp = stat->lpcount;
         if( lp->validfarkaslp == stat->lpcount-1 )
            lp->validfarkaslp = stat->lpcount;
      }
   }

   SCIPsetDebugMsg(set, "solved LP %" SCIP_LONGINT_FORMAT " with barrier%s (diving=%d, nduallps=%" SCIP_LONGINT_FORMAT ", nbarrierlps=%" SCIP_LONGINT_FORMAT ")\n",
      stat->lpcount, crossover ? "/crossover" : "", lp->diving || lp->probing, stat->nbarrierlps, stat->ndivinglps);

   return SCIP_OKAY;
}

/** solves the LP with the given algorithm */
static
SCIP_RETCODE lpAlgorithm(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            timelimit,          /**< pointer to store whether the time limit was hit */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_Real lptimelimit;
   SCIP_Bool success;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lperror != NULL);

   /* check if a time limit is set, and set time limit for LP solver accordingly */
   lptimelimit = SCIPlpiInfinity(lp->lpi);
   if( set->istimelimitfinite )
      lptimelimit = set->limit_time - SCIPclockGetTime(stat->solvingtime);

   success = FALSE;
   if( lptimelimit > 0.0 )
      SCIP_CALL( lpSetRealpar(lp, SCIP_LPPAR_LPTILIM, lptimelimit, &success) );

   if( lptimelimit <= 0.0 || !success )
   {
      SCIPsetDebugMsg(set, "time limit of %f seconds could not be set\n", lptimelimit);
      *lperror = ((lptimelimit > 0.0) ? TRUE : FALSE);
      *timelimit = TRUE;
      return SCIP_OKAY;
   }
   SCIPsetDebugMsg(set, "calling LP algorithm <%s> with a time limit of %f seconds\n", lpalgoName(lpalgo), lptimelimit);

   /* call appropriate LP algorithm */
   switch( lpalgo )
   {
   case SCIP_LPALGO_PRIMALSIMPLEX:
      SCIP_CALL( lpPrimalSimplex(lp, set, stat, resolve, keepsol, lperror) );
      break;

   case SCIP_LPALGO_DUALSIMPLEX:
      /* run dual lexicographic simplex if required */
      if( set->lp_lexdualalgo && (!set->lp_lexdualrootonly || stat->maxdepth == 0) && (!set->lp_lexdualstalling || lp->installing) )
      {
         SCIP_CALL( lpLexDualSimplex(lp, set, stat, resolve, keepsol, lperror) );
      }
      else
      {
         SCIP_CALL( lpDualSimplex(lp, set, stat, resolve, keepsol, lperror) );
      }
      break;

   case SCIP_LPALGO_BARRIER:
      SCIP_CALL( lpBarrier(lp, set, stat, FALSE, keepsol, lperror) );
      break;

   case SCIP_LPALGO_BARRIERCROSSOVER:
      SCIP_CALL( lpBarrier(lp, set, stat, TRUE, keepsol, lperror) );
      break;

   default:
      SCIPerrorMessage("invalid LP algorithm\n");
      return SCIP_INVALIDDATA;
   }

   if( !(*lperror) )
   {
      /* check for primal and dual feasibility */
      SCIP_CALL( SCIPlpiGetSolFeasibility(lp->lpi, &lp->primalfeasible, &lp->dualfeasible) );

      SCIPsetDebugMsg(set, "LP feasibility: primalfeasible=%u, dualfeasible=%u\n", lp->primalfeasible, lp->dualfeasible);
   }

   return SCIP_OKAY;
}

/** maximal number of verblevel-high messages about numerical trouble in LP that will be printed
 * when this number is reached and display/verblevel is not full, then further messages are suppressed in this run
 */
#define MAXNUMTROUBLELPMSGS 10

/** prints message about numerical trouble
 *
 * If message has verblevel at most high and display/verblevel is not full,
 * then the message is not printed if already MAXNUMTROUBLELPMSGS messages
 * were printed before in the current run.
 */
static
void lpNumericalTroubleMessage(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VERBLEVEL        verblevel,          /**< verbosity level of message */
   const char*           formatstr,          /**< message format string */
   ...                                       /**< arguments to format string */
   )
{
   va_list ap;

   assert(verblevel > SCIP_VERBLEVEL_NONE);
   assert(verblevel <= SCIP_VERBLEVEL_FULL);
   assert(set->disp_verblevel <= SCIP_VERBLEVEL_FULL);

   if( set->disp_verblevel < SCIP_VERBLEVEL_FULL )
   {
      if( verblevel <= SCIP_VERBLEVEL_HIGH )
      {
         /* if already max number of messages about numerical trouble in LP on verblevel at most high, then skip message */
         if( stat->nnumtroublelpmsgs > MAXNUMTROUBLELPMSGS )
            return;

         /* increase count on messages with verblevel high */
         ++stat->nnumtroublelpmsgs ;
      }

      /* if messages wouldn't be printed, then return already */
      if( verblevel > set->disp_verblevel )
         return;
   }

   /* print common begin of message */
   SCIPmessagePrintInfo(messagehdlr,
      "(node %" SCIP_LONGINT_FORMAT ") numerical troubles in LP %" SCIP_LONGINT_FORMAT " -- ",
      stat->nnodes, stat->nlps);

   /* print individual part of message */
   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(messagehdlr, NULL, formatstr, ap);
   va_end(ap);

   /* warn that further messages will be suppressed */
   if( set->disp_verblevel < SCIP_VERBLEVEL_FULL && verblevel <= SCIP_VERBLEVEL_HIGH && stat->nnumtroublelpmsgs > MAXNUMTROUBLELPMSGS )
   {
      SCIPmessagePrintInfo(messagehdlr, " -- further messages will be suppressed (use display/verblevel=5 to see all)");
   }

   /* print closing new-line */
   SCIPmessagePrintInfo(messagehdlr, "\n");
}

#define FEASTOLTIGHTFAC 0.001
/** solves the LP with the given LP algorithm, and tries to resolve numerical problems */
static
SCIP_RETCODE lpSolveStable(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
   int                   itlim,              /**< maximal number of LP iterations to perform in first LP calls (before solving from scratch), or -1 for no limit */
   int                   harditlim,          /**< maximal number of LP iterations to perform (hard limit for all LP calls), or -1 for no limit */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   int                   fastmip,            /**< which FASTMIP setting of LP solver should be used? */
   SCIP_Bool             tightprimfeastol,   /**< should a tighter primal feasibility tolerance be used? */
   SCIP_Bool             tightdualfeastol,   /**< should a tighter dual feasibility tolerance be used? */
   SCIP_Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            timelimit,          /**< pointer to store whether the time limit was hit */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_Bool success;
   SCIP_Bool success2;
   SCIP_Bool success3;
   SCIP_Bool simplex;
   SCIP_Bool itlimishard;
   SCIP_Bool usepolishing;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);
   assert(timelimit != NULL);

   *lperror = FALSE;

   /**@todo implement solving the LP when loose variables with infinite best bound are present; for this, we need to
    *       solve with deactivated objective limit in order to determine whether we are (a) infeasible or (b) feasible
    *       and hence unbounded; to handle case (b) we need to store an array of loose variables with best bound in
    *       SCIP_LP such that we can return a primal ray
    */
   if( lp->looseobjvalinf > 0 )
   {
      SCIPerrorMessage("cannot solve LP when loose variable with infinite best bound is present\n");
      return SCIP_ERROR;
   }

   /* check, whether we solve with a simplex algorithm */
   simplex = (lpalgo == SCIP_LPALGO_PRIMALSIMPLEX || lpalgo == SCIP_LPALGO_DUALSIMPLEX);

   /* check whether the iteration limit is a hard one */
   itlimishard = (itlim == harditlim);

   /* check whether solution polishing should be used */
   if( lp->lpihaspolishing && (set->lp_solutionpolishing == 2 || (set->lp_solutionpolishing == 1 && stat->nnodes == 1 && !lp->probing)
         || (set->lp_solutionpolishing == 3 && ((lp->probing && !lp->strongbranchprobing) || lp->diving))) )
   {
      usepolishing = TRUE;
      if( lp->updateintegrality )
      {
         SCIP_CALL( lpCopyIntegrality(lp, set) );
      }
   }
   else
      usepolishing = FALSE;

   /* solve with given settings (usually fast but imprecise) */
   if( SCIPsetIsInfinity(set, lp->cutoffbound) )
   {
      SCIP_CALL( lpSetObjlim(lp, set, lp->cutoffbound) );
   }
   else
   {
      SCIP_CALL( lpSetObjlim(lp, set, lp->cutoffbound - getFiniteLooseObjval(lp, set, prob)) );
   }
   SCIP_CALL( lpSetIterationLimit(lp, itlim) );
   SCIP_CALL( lpSetFeastol(lp, tightprimfeastol ? FEASTOLTIGHTFAC * SCIPsetLpfeastol(set) : SCIPsetLpfeastol(set), &success) );
   SCIP_CALL( lpSetDualfeastol(lp, tightdualfeastol ? FEASTOLTIGHTFAC * SCIPsetDualfeastol(set) : SCIPsetDualfeastol(set),
         &success) );
   SCIP_CALL( lpSetBarrierconvtol(lp, (tightprimfeastol || tightdualfeastol) ? FEASTOLTIGHTFAC * SCIPsetBarrierconvtol(set)
         : SCIPsetBarrierconvtol(set), &success) );
   SCIP_CALL( lpSetFromscratch(lp, fromscratch, &success) );
   SCIP_CALL( lpSetFastmip(lp, fastmip, &success) );
   SCIP_CALL( lpSetScaling(lp, set->lp_scaling, &success) );
   SCIP_CALL( lpSetPresolving(lp, set->lp_presolving, &success) );
   SCIP_CALL( lpSetRowrepswitch(lp, set->lp_rowrepswitch, &success) );
   SCIP_CALL( lpSetPricingChar(lp, set->lp_pricing) );
   SCIP_CALL( lpSetThreads(lp, set->lp_threads, &success) );
   SCIP_CALL( lpSetLPInfo(lp, set->disp_lpinfo) );
   SCIP_CALL( lpSetConditionLimit(lp, set->lp_conditionlimit, &success) );
   SCIP_CALL( lpSetTiming(lp, set->time_clocktype, set->time_enabled, &success) );
   SCIP_CALL( lpSetRandomseed(lp, SCIPsetInitializeRandomSeed(set, set->random_randomseed), &success) );
   SCIP_CALL( lpSetSolutionPolishing(lp, usepolishing, &success) );
   SCIP_CALL( lpSetRefactorInterval(lp, set->lp_refactorinterval, &success) );
   SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );
   resolve = FALSE; /* only the first solve should be counted as resolving call */

   /* check for stability; iteration limit exceeded is also treated like instability if the iteration limit is soft */
   if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi) && (itlimishard || !SCIPlpiIsIterlimExc(lp->lpi))) )
      return SCIP_OKAY;
   else if( !set->lp_checkstability )
   {
      SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
      if( success )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
         return SCIP_OKAY;
      }
   }

   /* In the following, whenever the LP iteration limit is exceeded in an LP solving call, we leave out the 
    * remaining resolving calls with changed settings and go directly to solving the LP from scratch.
    */

   /* if FASTMIP is turned on, solve again without FASTMIP (starts from the solution of the last LP solving call);
    * do this only if the iteration limit was not exceeded in the last LP solving call 
    */
   if( fastmip > 0 && simplex && ((*lperror) || !SCIPlpiIsIterlimExc(lp->lpi)) )
   {
      SCIP_CALL( lpSetFastmip(lp, 0, &success) );
      if( success )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again with %s without FASTMIP", lpalgoName(lpalgo));
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

         /* check for stability */
         if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi) && (itlimishard || !SCIPlpiIsIterlimExc(lp->lpi))) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }
      }
   }

   /* if the iteration limit was exceeded in the last LP solving call, we leave out the remaining resolving calls with changed settings
    * and go directly to solving the LP from scratch 
    */
   if( (*lperror) || !SCIPlpiIsIterlimExc(lp->lpi) )
   {
      /* solve again with opposite scaling setting (starts from the solution of the last LP solving call) */
      SCIP_CALL( lpSetScaling(lp, (set->lp_scaling > 0) ? 0 : 1, &success) );
      if( success )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again with %s %s scaling", lpalgoName(lpalgo), (set->lp_scaling == 0) ? "with" : "without");
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

         /* check for stability */
         if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi) && (itlimishard || !SCIPlpiIsIterlimExc(lp->lpi))) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }

         /* reset scaling */
         SCIP_CALL( lpSetScaling(lp, set->lp_scaling, &success) );
         assert(success);
      }
   }

   /* if the iteration limit was exceeded in the last LP solving call, we leave out the remaining resolving calls with changed settings
    * and go directly to solving the LP from scratch
    */
   if( (*lperror) || !SCIPlpiIsIterlimExc(lp->lpi) )
   {
      /* solve again with opposite presolving setting (starts from the solution of the last LP solving call) */
      SCIP_CALL( lpSetPresolving(lp, !set->lp_presolving, &success) );
      if( success )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again with %s %s presolving",
            lpalgoName(lpalgo), !set->lp_presolving ? "with" : "without");
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

         /* check for stability */
         if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi) && (itlimishard || !SCIPlpiIsIterlimExc(lp->lpi))) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }

         /* reset presolving */
         SCIP_CALL( lpSetPresolving(lp, set->lp_presolving, &success) );
         assert(success);
      }
   }

   /* solve again with a tighter feasibility tolerance (starts from the solution of the last LP solving call);
    * do this only if the iteration limit was not exceeded in the last LP solving call 
    */
   if( ((simplex && (!tightprimfeastol || !tightdualfeastol)) || (!tightprimfeastol && !tightdualfeastol)) &&
      ((*lperror) || !SCIPlpiIsIterlimExc(lp->lpi)) )
   {
      success = FALSE;
      if( !tightprimfeastol )
      {
         SCIP_CALL( lpSetFeastol(lp, FEASTOLTIGHTFAC * SCIPsetLpfeastol(set), &success) );
      }

      success2 = FALSE;
      if( !tightdualfeastol )
      {
         SCIP_CALL( lpSetDualfeastol(lp, FEASTOLTIGHTFAC * SCIPsetDualfeastol(set), &success2) );
      }

      success3 = FALSE;
      if( !simplex && !tightprimfeastol && !tightdualfeastol )
      {
         SCIP_CALL( lpSetBarrierconvtol(lp, FEASTOLTIGHTFAC * SCIPsetBarrierconvtol(set), &success3) );
      }

      if( success || success2 || success3 )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again with %s with tighter primal and dual feasibility tolerance", lpalgoName(lpalgo));
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

         /* check for stability */
         if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi) && (itlimishard || !SCIPlpiIsIterlimExc(lp->lpi))) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }

         /* reset feasibility tolerance */
         if( !tightprimfeastol )
         {
            SCIP_CALL( lpSetFeastol(lp, SCIPsetLpfeastol(set), &success) );
         }
         if( !tightdualfeastol )
         {
            SCIP_CALL( lpSetDualfeastol(lp, SCIPsetDualfeastol(set), &success) );
         }
         if( !simplex && !tightprimfeastol && !tightdualfeastol )
         {
            SCIP_CALL( lpSetBarrierconvtol(lp, SCIPsetBarrierconvtol(set), &success) );
         }
      }
   }

   /* all LPs solved after this point are solved from scratch, so set the LP iteration limit to the hard limit;
    * the given iteration limit might be a soft one to restrict resolving calls only */
   SCIP_CALL( lpSetIterationLimit(lp, harditlim) );

   /* if not already done, solve again from scratch */
   if( !fromscratch && simplex )
   {
      SCIP_CALL( lpSetFromscratch(lp, TRUE, &success) );
      if( success )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again from scratch with %s", lpalgoName(lpalgo));
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

         /* check for stability */
         if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi)) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }
      }
   }

   /* solve again, use other simplex this time */
   if( simplex )
   {
      lpalgo = (lpalgo == SCIP_LPALGO_PRIMALSIMPLEX ? SCIP_LPALGO_DUALSIMPLEX : SCIP_LPALGO_PRIMALSIMPLEX);
      lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again from scratch with %s", lpalgoName(lpalgo));
      SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

      /* check for stability */
      if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi)) )
         return SCIP_OKAY;
      else if( !set->lp_checkstability )
      {
         SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
         if( success )
         {
            lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
            return SCIP_OKAY;
         }
      }

      /* solve again with opposite scaling and other simplex */
      SCIP_CALL( lpSetScaling(lp, (set->lp_scaling > 0) ? 0 : 1, &success) );
      if( success )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again from scratch with %s %s scaling",
            lpalgoName(lpalgo), (set->lp_scaling == 0) ? "with" : "without");
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

         /* check for stability */
         if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi)) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }

         /* reset scaling */
         SCIP_CALL( lpSetScaling(lp, set->lp_scaling, &success) );
         assert(success);
      }

      /* solve again with opposite presolving and other simplex */
      SCIP_CALL( lpSetPresolving(lp, !set->lp_presolving, &success) );
      if( success )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again from scratch with %s %s presolving",
            lpalgoName(lpalgo), !set->lp_presolving ? "with" : "without");
         SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

         /* check for stability */
         if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi)) )
            return SCIP_OKAY;
         else if( !set->lp_checkstability )
         {
            SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
            if( success )
            {
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
               return SCIP_OKAY;
            }
         }

         /* reset presolving */
         SCIP_CALL( lpSetPresolving(lp, set->lp_presolving, &success) );
         assert(success);
      }

      /* solve again with tighter feasibility tolerance, use other simplex this time */
      if( !tightprimfeastol || !tightdualfeastol )
      {
         success = FALSE;
         if( !tightprimfeastol )
         {
            SCIP_CALL( lpSetFeastol(lp, FEASTOLTIGHTFAC * SCIPsetLpfeastol(set), &success) );
         }

         success2 = FALSE;
         if( !tightdualfeastol )
         {
            SCIP_CALL( lpSetDualfeastol(lp, FEASTOLTIGHTFAC * SCIPsetDualfeastol(set), &success2) );
         }

         if( success || success2 )
         {
            lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "solve again from scratch with %s with tighter feasibility tolerance", lpalgoName(lpalgo));
            SCIP_CALL( lpAlgorithm(lp, set, stat, lpalgo, resolve, keepsol, timelimit, lperror) );

            /* check for stability */
            if( *timelimit || (!(*lperror) && SCIPlpiIsStable(lp->lpi)) )
               return SCIP_OKAY;
            else if( !set->lp_checkstability )
            {
               SCIP_CALL( SCIPlpiIgnoreInstability(lp->lpi, &success) );
               if( success )
               {
                  lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "ignoring instability of %s", lpalgoName(lpalgo));
                  return SCIP_OKAY;
               }
            }

            /* reset feasibility tolerance */
            if( !tightprimfeastol )
            {
               SCIP_CALL( lpSetFeastol(lp, SCIPsetLpfeastol(set), &success) );
            }
            if( !tightdualfeastol )
            {
               SCIP_CALL( lpSetDualfeastol(lp, SCIPsetDualfeastol(set), &success) );
            }
         }
      }
   }

   /* nothing worked -- exit with an LPERROR */
   lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_HIGH, "unresolved");
   *lperror = TRUE;

   return SCIP_OKAY;
}

/** adjust the LP objective value if its greater/less than +/- SCIPsetInfinity() */
static
void adjustLPobjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert(lp != NULL);
   assert(set != NULL);

   if( SCIPsetIsInfinity(set, lp->lpobjval) && lp->lpobjval != SCIPsetInfinity(set) ) /*lint !e777*/
   {
      if( !lp->adjustlpval )
      {
         SCIPmessagePrintWarning(messagehdlr, "LP solution value is above SCIP's infinity value\n");
         lp->adjustlpval = TRUE;
      }
      lp->lpobjval = SCIPsetInfinity(set);
   }
   else if( SCIPsetIsInfinity(set, -lp->lpobjval) && lp->lpobjval != -SCIPsetInfinity(set) ) /*lint !e777*/
   {
      if( !lp->adjustlpval )
      {
         SCIPmessagePrintWarning(messagehdlr, "LP solution value is below SCIP's -infinity value\n");
         lp->adjustlpval = TRUE;
      }
      lp->lpobjval = -SCIPsetInfinity(set);
   }
}

/** solves the LP with the given algorithm and evaluates return status */
static
SCIP_RETCODE lpSolve(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LPALGO           lpalgo,             /**< LP algorithm that should be applied */
   int                   resolveitlim,       /**< maximal number of LP iterations to perform in resolving calls, or -1 for no limit */
   int                   harditlim,          /**< maximal number of LP iterations to perform (hard limit for all LP calls), or -1 for no limit */
   SCIP_Bool             needprimalray,      /**< if the LP is unbounded, do we need a primal ray? */
   SCIP_Bool             needdualray,        /**< if the LP is infeasible, do we need a dual ray? */
   SCIP_Bool             resolve,            /**< is this a resolving call (starting with feasible basis)? */
   int                   fastmip,            /**< which FASTMIP setting of LP solver should be used? */
   SCIP_Bool             tightprimfeastol,   /**< should a tighter primal feasibility tolerance be used? */
   SCIP_Bool             tightdualfeastol,   /**< should a tighter dual feasibility tolerance be used? */
   SCIP_Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_Bool solvedprimal;
   SCIP_Bool solveddual;
   SCIP_Bool timelimit;
   int itlim;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lperror != NULL);

   checkLinks(lp);

   solvedprimal = FALSE;
   solveddual = FALSE;
   timelimit = FALSE;

   /* select the basic iteration limit depending on whether this is a resolving call or not */
   itlim = ( resolve ? resolveitlim : harditlim );

 SOLVEAGAIN:
   /* call simplex */
   SCIP_CALL( lpSolveStable(lp, set, messagehdlr, stat, prob, lpalgo, itlim, harditlim, resolve, fastmip, tightprimfeastol, tightdualfeastol, fromscratch,
         keepsol, &timelimit, lperror) );
   resolve = FALSE; /* only the first solve should be counted as resolving call */
   solvedprimal = solvedprimal || (lp->lastlpalgo == SCIP_LPALGO_PRIMALSIMPLEX);
   solveddual = solveddual || (lp->lastlpalgo == SCIP_LPALGO_DUALSIMPLEX);

   /* check, if an error occurred */
   if( *lperror )
   {
      SCIPsetDebugMsg(set, "unresolved error while solving LP with %s\n", lpalgoName(lp->lastlpalgo));
      lp->solved = FALSE;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      return SCIP_OKAY;
   }

   /* check, if a time limit was exceeded */
   if( timelimit )
   {
      SCIPsetDebugMsg(set, "time limit exceeded before solving LP\n");
      lp->solved = TRUE;
      lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
      lp->lpobjval = -SCIPsetInfinity(set);
      return SCIP_OKAY;
   }

   /* only one should return true */
   assert(!(SCIPlpiIsOptimal(lp->lpi) && SCIPlpiIsObjlimExc(lp->lpi) && SCIPlpiIsPrimalInfeasible(lp->lpi) && 
         SCIPlpiExistsPrimalRay(lp->lpi) && SCIPlpiIsIterlimExc(lp->lpi) && SCIPlpiIsTimelimExc(lp->lpi)));

   /* evaluate solution status */
   if( SCIPlpiIsOptimal(lp->lpi) )
   {
      assert(lp->primalfeasible);
      assert(lp->dualfeasible);
      lp->lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
      SCIP_CALL( SCIPlpiGetObjval(lp->lpi, &lp->lpobjval) );
      adjustLPobjval(lp, set, messagehdlr);

      if( !SCIPsetIsInfinity(set, lp->lpiobjlim) && SCIPsetIsGE(set, lp->lpobjval, lp->lpiobjlim) )
      {
         /* the solver may return the optimal value, even if this is greater or equal than the upper bound */
         SCIPsetDebugMsg(set, "optimal solution %.15g exceeds objective limit %.15g\n", lp->lpobjval, lp->lpiobjlim);
         lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
         lp->lpobjval = SCIPsetInfinity(set);
      }
      /* if we did not disable the cutoff bound in the LP solver, the LP solution status should be objective limit
       * reached if the LP objective value is greater than the cutoff bound
       */
      assert(lpCutoffDisabled(set) || lp->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || SCIPsetIsInfinity(set, lp->cutoffbound)
         || SCIPsetIsLE(set, lp->lpobjval + getFiniteLooseObjval(lp, set, prob), lp->cutoffbound));
   }
   else if( SCIPlpiIsObjlimExc(lp->lpi) )
   {
      assert(!lpCutoffDisabled(set));
      lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
      lp->lpobjval = SCIPsetInfinity(set);
   }
   else if( SCIPlpiIsPrimalInfeasible(lp->lpi) )
   {
      /* because of numerical instability lpalgo != lp->lastlpalgo might happen - hence, we have to check both */
      if( needdualray && !SCIPlpiHasDualRay(lp->lpi) && !solveddual && lpalgo != SCIP_LPALGO_DUALSIMPLEX )
      {
         assert(lp->lastlpalgo != SCIP_LPALGO_DUALSIMPLEX);
         lpalgo = SCIP_LPALGO_DUALSIMPLEX;
         goto SOLVEAGAIN;
      }
      lp->lpsolstat = SCIP_LPSOLSTAT_INFEASIBLE;
      lp->lpobjval = SCIPsetInfinity(set);
   }
   else if( SCIPlpiExistsPrimalRay(lp->lpi) )
   {
      /* because of numerical instability lpalgo != lp->lastlpalgo might happen - hence, we have to check both */
      if( needprimalray && !SCIPlpiHasPrimalRay(lp->lpi) && !solvedprimal && lpalgo != SCIP_LPALGO_PRIMALSIMPLEX )
      {
         assert(lp->lastlpalgo != SCIP_LPALGO_PRIMALSIMPLEX);
         lpalgo = SCIP_LPALGO_PRIMALSIMPLEX;
         goto SOLVEAGAIN;
      }
      lp->lpsolstat = SCIP_LPSOLSTAT_UNBOUNDEDRAY;
      lp->lpobjval = -SCIPsetInfinity(set);
   }
   else if( SCIPlpiIsIterlimExc(lp->lpi) )
   {
      SCIP_CALL( SCIPlpiGetObjval(lp->lpi, &lp->lpobjval) );
      adjustLPobjval(lp, set, messagehdlr);
      lp->lpsolstat = SCIP_LPSOLSTAT_ITERLIMIT;
   }
   else if( SCIPlpiIsTimelimExc(lp->lpi) )
   {
      lp->lpobjval = -SCIPsetInfinity(set);
      lp->lpsolstat = SCIP_LPSOLSTAT_TIMELIMIT;
   }
   else if( !solveddual && lpalgo != SCIP_LPALGO_DUALSIMPLEX)
   {
      assert(lp->lastlpalgo != SCIP_LPALGO_DUALSIMPLEX);
      lpalgo = SCIP_LPALGO_DUALSIMPLEX;
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %" SCIP_LONGINT_FORMAT ") solution status of LP %" SCIP_LONGINT_FORMAT " could not be proven (internal status:%d) -- solve again with %s\n",
         stat->nnodes, stat->nlps, SCIPlpiGetInternalStatus(lp->lpi), lpalgoName(lpalgo));
      goto SOLVEAGAIN;
   }
   else if( !solvedprimal && lpalgo != SCIP_LPALGO_PRIMALSIMPLEX)
   {
      assert(lp->lastlpalgo != SCIP_LPALGO_PRIMALSIMPLEX);
      lpalgo = SCIP_LPALGO_PRIMALSIMPLEX;
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "(node %" SCIP_LONGINT_FORMAT ") solution status of LP %" SCIP_LONGINT_FORMAT " could not be proven (internal status:%d) -- solve again with %s\n",
         stat->nnodes, stat->nlps, SCIPlpiGetInternalStatus(lp->lpi), lpalgoName(lpalgo));
      goto SOLVEAGAIN;
   }
   else
   {
      SCIPerrorMessage("(node %" SCIP_LONGINT_FORMAT ") error or unknown return status of %s in LP %" SCIP_LONGINT_FORMAT " (internal status: %d)\n",
         stat->nnodes, lpalgoName(lp->lastlpalgo), stat->nlps, SCIPlpiGetInternalStatus(lp->lpi));
      lp->lpsolstat = SCIP_LPSOLSTAT_ERROR;
      return SCIP_LPERROR;
   }

   lp->solved = TRUE;

   SCIPsetDebugMsg(set, "solving LP with %s returned solstat=%d (internal status: %d, primalfeasible=%u, dualfeasible=%u)\n",
      lpalgoName(lp->lastlpalgo), lp->lpsolstat, SCIPlpiGetInternalStatus(lp->lpi),
      SCIPlpiIsPrimalFeasible(lp->lpi), SCIPlpiIsDualFeasible(lp->lpi));

   return SCIP_OKAY;
}

/** flushes the LP and solves it with the primal or dual simplex algorithm, depending on the current basis feasibility */
static
SCIP_RETCODE lpFlushAndSolve(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   resolveitlim,       /**< maximal number of LP iterations to perform in resolving calls, or -1 for no limit */
   int                   harditlim,          /**< maximal number of LP iterations to perform (hard limit for all LP calls), or -1 for no limit */
   SCIP_Bool             needprimalray,      /**< if the LP is unbounded, do we need a primal ray? */
   SCIP_Bool             needdualray,        /**< if the LP is infeasible, do we need a dual ray? */
   int                   fastmip,            /**< which FASTMIP setting of LP solver should be used? */
   SCIP_Bool             tightprimfeastol,   /**< should a tighter primal feasibility tolerance be used? */
   SCIP_Bool             tightdualfeastol,   /**< should a tighter dual feasibility tolerance be used? */
   SCIP_Bool             fromscratch,        /**< should the LP be solved from scratch without using current basis? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_Bool resolve;
   char algo;

   assert(lp != NULL);
   assert(set != NULL);
   assert(lperror != NULL);

   /* flush changes to the LP solver */
   SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );
   fastmip = ((!lp->flushaddedcols && !lp->flushdeletedcols) ? fastmip : 0); /* turn off FASTMIP if columns were changed */

   /* select LP algorithm to apply */
   resolve = lp->solisbasic && (lp->dualfeasible || lp->primalfeasible) && !fromscratch;
   algo = resolve ? set->lp_resolvealgorithm : set->lp_initalgorithm;

   switch( algo )
   {
   case 's':
      /* select simplex method */
      if( lp->dualfeasible || !lp->primalfeasible )
      {
         SCIPsetDebugMsg(set, "solving dual LP\n");
         SCIP_CALL( lpSolve(lp, set, messagehdlr, stat, prob, SCIP_LPALGO_DUALSIMPLEX, resolveitlim, harditlim, needprimalray,
               needdualray, resolve, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );
      }
      else
      {
         SCIPsetDebugMsg(set, "solving primal LP\n");
         SCIP_CALL( lpSolve(lp, set, messagehdlr, stat, prob, SCIP_LPALGO_PRIMALSIMPLEX, resolveitlim, harditlim, needprimalray,
               needdualray, resolve, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );
      }
      break;

   case 'p':
      SCIPsetDebugMsg(set, "solving primal LP\n");
      SCIP_CALL( lpSolve(lp, set, messagehdlr, stat, prob, SCIP_LPALGO_PRIMALSIMPLEX, resolveitlim, harditlim, needprimalray,
            needdualray, resolve, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );
      break;

   case 'd':
      SCIPsetDebugMsg(set, "solving dual LP\n");
      SCIP_CALL( lpSolve(lp, set, messagehdlr, stat, prob, SCIP_LPALGO_DUALSIMPLEX, resolveitlim, harditlim, needprimalray,
            needdualray, resolve, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );
      break;

   case 'b':
      SCIPsetDebugMsg(set, "solving barrier LP\n");
      SCIP_CALL( lpSolve(lp, set, messagehdlr, stat, prob, SCIP_LPALGO_BARRIER, resolveitlim, harditlim, needprimalray,
            needdualray, resolve, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );
      break;

   case 'c':
      SCIPsetDebugMsg(set, "solving barrier LP with crossover\n");
      SCIP_CALL( lpSolve(lp, set, messagehdlr, stat, prob, SCIP_LPALGO_BARRIERCROSSOVER, resolveitlim, harditlim, needprimalray,
            needdualray, resolve, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );
      break;

   default:
      SCIPerrorMessage("invalid parameter setting <%c> for LP algorithm\n", algo);
      return SCIP_PARAMETERWRONGVAL;
   }
   assert(!(*lperror) || !lp->solved);

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** checks if the lazy bounds are valid */
static
void checkLazyBounds(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_COL* col;
   int c;

   assert(lp->flushed);

   for( c = 0; c < lp->nlazycols; ++c )
   {
      col = lp->lazycols[c];

      /* in case lazy bounds are given, check that the primal solution satisfies them */
      assert(SCIPsetIsInfinity(set, -col->lazylb) || SCIPsetIsFeasGE(set, col->primsol, col->lazylb));
      assert(SCIPsetIsInfinity(set, col->lazyub) || SCIPsetIsFeasLE(set, col->primsol, col->lazyub));
   }
}
#else
#define checkLazyBounds(lp, set) /**/
#endif

/** marks all lazy columns to be changed; this is needed for reloading/removing bounds of these columns before and after
 *  diving
 */
static
SCIP_RETCODE updateLazyBounds(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_COL* col;
   int c;

   assert(lp->nlazycols > 0);

   /* return, if we are in diving, and bounds were already applied
    * or if we are not in diving and bounds were not applied
    */
   if( lp->diving == lp->divinglazyapplied )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "mark all lazy columns as changed in order to reload bounds (diving=%u, applied=%u)\n",
      lp->diving, lp->divinglazyapplied);

   for( c = 0; c < lp->nlazycols; ++c )
   {
      col = lp->lazycols[c];

      /* if the column has a lazy lower bound, mark its lower bounds as changed */
      if( !SCIPsetIsInfinity(set, -col->lazylb) )
      {
         assert((!(lp->divinglazyapplied)) || (col->flushedlb == col->lb)); /*lint !e777*/
         assert(lp->divinglazyapplied || SCIPsetIsGT(set, col->lb, col->lazylb)
            || (col->flushedlb == -SCIPlpiInfinity(lp->lpi))); /*lint !e777*/

         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lp) );

         /* mark bound change in the column */
         col->lbchanged = TRUE;
      }

      /* if the column has a lazy upper bound, mark its upper bounds as changed */
      if( !SCIPsetIsInfinity(set, col->lazyub) )
      {
         assert((!(lp->divinglazyapplied)) || (col->flushedub == col->ub)); /*lint !e777*/
         assert(lp->divinglazyapplied || SCIPsetIsLT(set, col->ub, col->lazyub)
            || (col->flushedub == SCIPlpiInfinity(lp->lpi))); /*lint !e777*/

         /* insert column in the chgcols list (if not already there) */
         SCIP_CALL( insertColChgcols(col, set, lp) );

         /* mark bound change in the column */
         col->ubchanged = TRUE;
      }
   }

   /* update lp->divinglazyapplied flag: if we are in diving mode, we just applied the lazy bounds,
    * if not, we just removed them
    */
   lp->divinglazyapplied = lp->diving;

   return SCIP_OKAY;
}

/** returns the iteration limit for an LP resolving call */
static
int lpGetResolveItlim(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   itlim               /**< hard iteration limit */
   )
{
   /* no limit set or average not yet reliable */
   if( (set->lp_resolveiterfac == -1) || stat->nlps - stat->nrootlps < 5 )
      return itlim;
   /* set itlim to INT_MAX if it is -1 to reduce the number of cases to be regarded in the following */
   if( itlim == -1 )
      itlim = INT_MAX;
   /* return resolveiterfac * average iteration number per call after root, but at least resolveitermin and at most the hard iteration limit */
   return (int) MIN(itlim, MAX(set->lp_resolveitermin, \
         (set->lp_resolveiterfac * (stat->nlpiterations - stat->nrootlpiterations) / (SCIP_Real)(stat->nlps - stat->nrootlps))));
}



/** solves the LP with simplex algorithm, and copy the solution into the column's data */
SCIP_RETCODE SCIPlpSolveAndEval(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Longint          itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool             limitresolveiters,  /**< should LP iterations for resolving calls be limited?
                                              *   (limit is computed within the method w.r.t. the average LP iterations) */
   SCIP_Bool             aging,              /**< should aging and removal of obsolete cols/rows be applied? */
   SCIP_Bool             keepsol,            /**< should the old LP solution be kept if no iterations were performed? */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Bool needprimalray;
   SCIP_Bool needdualray;
   int harditlim;
   int resolveitlim;

   assert(lp != NULL);
   assert(prob != NULL);
   assert(prob->nvars >= lp->ncols);
   assert(lperror != NULL);

   SCIPsetDebugMsg(set, "solving LP: %d rows, %d cols, primalfeasible=%u, dualfeasible=%u, solved=%u, diving=%u, probing=%u, cutoffbnd=%g\n",
      lp->nrows, lp->ncols, lp->primalfeasible, lp->dualfeasible, lp->solved, lp->diving, lp->probing, lp->cutoffbound);

   retcode = SCIP_OKAY;
   *lperror = FALSE;

   /* check whether we need a proof of unboundedness or infeasibility by a primal or dual ray */
   needprimalray = TRUE;
   needdualray = (!SCIPprobAllColsInLP(prob, set, lp) || set->misc_exactsolve
      || (set->conf_enable && set->conf_useinflp != 'o'));

   /* compute the limit for the number of LP resolving iterations, if needed (i.e. if limitresolveiters == TRUE) */
   harditlim = (int) MIN(itlim, INT_MAX);
   resolveitlim = ( limitresolveiters ? lpGetResolveItlim(set, stat, harditlim) : harditlim );
   assert(harditlim == -1 || (resolveitlim <= harditlim));

   /* if there are lazy bounds, check whether the bounds should explicitly be put into the LP (diving was started)
    * or removed from the LP (diving was ended)
    */
   if( lp->nlazycols > 0 )
   {
      /* @todo avoid loosing primal feasibility here after changing the objective already did destroy dual feasibility;
       * first resolve LP?
       */
      SCIP_CALL( updateLazyBounds(lp, set) );
      assert(lp->diving == lp->divinglazyapplied);
   }

   /* flush changes to the LP solver */
   SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );
   assert(lp->flushed);

   /* if the time limit was reached in the last call and the LP did not change, lp->solved is set to TRUE, but we want
    * to run again anyway, since there seems to be some time left / the time limit was increased
    */
   if( !lp->solved || (lp->lpsolstat == SCIP_LPSOLSTAT_TIMELIMIT && stat->status != SCIP_STATUS_TIMELIMIT) )
   {
      SCIP_Bool* primalfeaspointer;
      SCIP_Bool* dualfeaspointer;
      SCIP_Bool primalfeasible;
      SCIP_Bool dualfeasible;
      SCIP_Bool rayfeasible;
      SCIP_Bool tightprimfeastol;
      SCIP_Bool tightdualfeastol;
      SCIP_Bool fromscratch;
      SCIP_Bool wasfromscratch;
      SCIP_Longint oldnlps;
      int fastmip;

      /* set initial LP solver settings */
      fastmip = ((lp->lpihasfastmip && !lp->flushaddedcols && !lp->flushdeletedcols && stat->nnodes > 1) ? set->lp_fastmip : 0);
      tightprimfeastol = FALSE;
      tightdualfeastol = FALSE;
      fromscratch = FALSE;
      primalfeasible = FALSE;
      dualfeasible = FALSE;
      wasfromscratch = (stat->nlps == 0);

   SOLVEAGAIN:
      /* solve the LP */
      oldnlps = stat->nlps;
      SCIP_CALL( lpFlushAndSolve(lp, blkmem, set, messagehdlr, stat, prob, eventqueue, resolveitlim, harditlim, needprimalray,
            needdualray, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );
      SCIPsetDebugMsg(set, "lpFlushAndSolve() returned solstat %d (error=%u)\n", SCIPlpGetSolstat(lp), *lperror);
      assert(!(*lperror) || !lp->solved);

      /* check for error */
      if( *lperror )
      {
         retcode = SCIP_OKAY;
         goto TERMINATE;
      }

      /* evaluate solution status */
      switch( SCIPlpGetSolstat(lp) )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         /* get LP solution and possibly check the solution's feasibility again */
         if( set->lp_checkprimfeas )
         {
            primalfeaspointer = &primalfeasible;
            lp->primalchecked = TRUE;
         }
         else
         {
            /* believe in the primal feasibility of the LP solution */
            primalfeasible = TRUE;
            primalfeaspointer = NULL;
            lp->primalchecked = FALSE;
         }
         if( set->lp_checkdualfeas )
         {
            dualfeaspointer = &dualfeasible;
            lp->dualchecked = TRUE;
         }
         else
         {
            /* believe in the dual feasibility of the LP solution */
            dualfeasible = TRUE;
            dualfeaspointer = NULL;
            lp->dualchecked = FALSE;
         }

         SCIP_CALL( SCIPlpGetSol(lp, set, stat, primalfeaspointer, dualfeaspointer) );

         /* in debug mode, check that lazy bounds (if present) are not violated */
         checkLazyBounds(lp, set);

         if( primalfeasible && dualfeasible && aging && !lp->diving && stat->nlps > oldnlps )
         {
            /* update ages and remove obsolete columns and rows from LP */
            SCIP_CALL( SCIPlpUpdateAges(lp, stat) );
            if( stat->nlps % ((set->lp_rowagelimit+1)/2 + 1) == 0 ) /*lint !e776*/
            {
               SCIP_CALL( SCIPlpRemoveNewObsoletes(lp, blkmem, set, stat, eventqueue, eventfilter) );
            }

            if( !lp->solved )
            {
               /* resolve LP after removing obsolete columns and rows */
               SCIPsetDebugMsg(set, "removed obsoletes - resolve LP again: %d rows, %d cols\n", lp->nrows, lp->ncols);
               aging = FALSE; /* to prevent infinite loops */
               goto SOLVEAGAIN;
            }
         }
         if( !primalfeasible || !dualfeasible )
         {
            SCIP_Bool simplex = (lp->lastlpalgo == SCIP_LPALGO_PRIMALSIMPLEX || lp->lastlpalgo == SCIP_LPALGO_DUALSIMPLEX);

            if( (fastmip > 0) && simplex )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again without FASTMIP */
               SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %" SCIP_LONGINT_FORMAT ") solution of LP %" SCIP_LONGINT_FORMAT " not optimal (pfeas=%d, dfeas=%d) -- solving again without FASTMIP\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               fastmip = 0;
               goto SOLVEAGAIN;
            }
            else if( (!primalfeasible && !tightprimfeastol) || (!dualfeasible && !tightdualfeastol) )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again with tighter feasibility
                * tolerance
                */
               SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %" SCIP_LONGINT_FORMAT ") solution of LP %" SCIP_LONGINT_FORMAT " not optimal (pfeas=%d, dfeas=%d) -- solving again with tighter feasibility tolerance\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               tightprimfeastol = tightprimfeastol || !primalfeasible;
               tightdualfeastol = tightdualfeastol || !dualfeasible;
               goto SOLVEAGAIN;
            }
            else if( !fromscratch && !wasfromscratch && simplex )
            {
               /* solution is infeasible (this can happen due to numerical problems): solve again from scratch */
               SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %" SCIP_LONGINT_FORMAT ") solution of LP %" SCIP_LONGINT_FORMAT " not optimal (pfeas=%d, dfeas=%d) -- solving again from scratch\n",
                  stat->nnodes, stat->nlps, primalfeasible, dualfeasible);
               fromscratch = TRUE;
               goto SOLVEAGAIN;
            }
            else
            {
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "unresolved");
               lp->solved = FALSE;
               lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
               *lperror = TRUE;
            }
         }
         SCIPsetDebugMsg(set, " -> LP objective value: %g + %g = %g (solstat=%d, cutoff=%g)\n",
            lp->lpobjval, getFiniteLooseObjval(lp, set, prob), lp->lpobjval + getFiniteLooseObjval(lp, set, prob),
            lp->lpsolstat, lp->cutoffbound);
         break;

      case SCIP_LPSOLSTAT_INFEASIBLE:
         SCIPsetDebugMsg(set, " -> LP infeasible\n");
         if( !SCIPprobAllColsInLP(prob, set, lp) || set->misc_exactsolve )
         {
            if( SCIPlpiHasDualRay(lp->lpi) )
            {
               SCIP_CALL( SCIPlpGetDualfarkas(lp, set, stat) );
            }
            /* it might happen that we have no infeasibility proof for the current LP (e.g. if the LP was always solved
             * with the primal simplex due to numerical problems) - treat this case like an LP error
             */
            else
            {
               SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %" SCIP_LONGINT_FORMAT ") infeasibility of LP %" SCIP_LONGINT_FORMAT " could not be proven by dual ray\n", stat->nnodes, stat->nlps);
               lp->solved = FALSE;
               lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
               *lperror = TRUE;
            }
         }
         break;

      case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
         if( set->lp_checkprimfeas )
         {
            /* get unbounded LP solution and check the solution's feasibility again */
            SCIP_CALL( SCIPlpGetUnboundedSol(lp, set, stat, &primalfeasible, &rayfeasible) );

            lp->primalchecked = TRUE;
         }
         else
         {
            /* get unbounded LP solution believing in the feasibility of the LP solution */
            SCIP_CALL( SCIPlpGetUnboundedSol(lp, set, stat, NULL, NULL) );

            primalfeasible = TRUE;
            rayfeasible = TRUE;
            lp->primalchecked = FALSE;
         }

         /* in debug mode, check that lazy bounds (if present) are not violated */
         checkLazyBounds(lp, set);

         SCIPsetDebugMsg(set, " -> LP has unbounded primal ray (primalfeas=%u, rayfeas=%u)\n",
            primalfeasible, rayfeasible);

         if( !primalfeasible || !rayfeasible )
         {
            SCIP_Bool simplex = (lp->lastlpalgo == SCIP_LPALGO_PRIMALSIMPLEX || lp->lastlpalgo == SCIP_LPALGO_DUALSIMPLEX);

            if( (fastmip > 0) && simplex )
            {
               /* unbounded solution is infeasible (this can happen due to numerical problems): solve again without FASTMIP */
               SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %" SCIP_LONGINT_FORMAT ") solution of unbounded LP %" SCIP_LONGINT_FORMAT " not optimal (pfeas=%d, rfeas=%d) -- solving again without FASTMIP\n",
                  stat->nnodes, stat->nlps, primalfeasible, rayfeasible);
               fastmip = 0;
               goto SOLVEAGAIN;
            }
            else if( !tightprimfeastol )
            {
               /* unbounded solution is infeasible (this can happen due to numerical problems): solve again with tighter feasibility
                * tolerance
                */
               SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %" SCIP_LONGINT_FORMAT ") solution of unbounded LP %" SCIP_LONGINT_FORMAT " not optimal (pfeas=%d, rfeas=%d) -- solving again with tighter primal feasibility tolerance\n",
                  stat->nnodes, stat->nlps, primalfeasible, rayfeasible);
               tightprimfeastol = TRUE;
               goto SOLVEAGAIN;
            }
            else if( !fromscratch && simplex )
            {
               /* unbounded solution is infeasible (this can happen due to numerical problems): solve again from scratch */
               SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
                  "(node %" SCIP_LONGINT_FORMAT ") solution of unbounded LP %" SCIP_LONGINT_FORMAT " not optimal (pfeas=%d, rfeas=%d) -- solving again from scratch\n",
                  stat->nnodes, stat->nlps, primalfeasible, rayfeasible);
               fromscratch = TRUE;
               goto SOLVEAGAIN;
            }
            else
            {
               /* unbounded solution is infeasible (this can happen due to numerical problems) and nothing helped:
                * forget about the LP at this node and mark it to be unsolved
                */
               lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "unresolved, LP unbounded");
               lp->solved = FALSE;
               lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
               *lperror = TRUE;
            }
         }

         break;

      case SCIP_LPSOLSTAT_OBJLIMIT:
         assert(!lpCutoffDisabled(set));
         /* if we do branch-and-price, make sure that a dual feasible solution exists, that exceeds the objective limit;
          * With FASTMIP setting, CPLEX does not apply the final pivot to reach the dual solution exceeding the objective
          * limit. Therefore, we have to either turn off FASTMIP and resolve the problem or continue solving it without
          * objective limit for at least one iteration. We first try to continue with FASTMIP for one additional simplex
          * iteration using the steepest edge pricing rule. If this does not fix the problem, we temporarily disable
          * FASTMIP and solve again.
          */
         if( !SCIPprobAllColsInLP(prob, set, lp) && fastmip )
         {
            SCIP_LPI* lpi;
            SCIP_Real objval;

            lpi = SCIPlpGetLPI(lp);

            assert(lpi != NULL);
            /* actually, SCIPsetIsGE(set, lp->lpobjval, lp->lpiuobjlim) should hold, but we are a bit less strict in
             * the assert by using !SCIPsetIsFeasNegative()
             */
            assert(SCIPlpiIsObjlimExc(lpi) || !SCIPsetIsFeasNegative(set, lp->lpobjval - lp->lpiobjlim));

            SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );

            /* do one additional simplex step if the computed dual solution doesn't exceed the objective limit */
            if( SCIPsetIsLT(set, objval, lp->lpiobjlim) )
            {
               SCIP_Real tmpcutoff;
               char tmppricingchar;
               SCIP_LPSOLSTAT solstat;

               SCIPsetDebugMsg(set, "objval = %f < %f = lp->lpiobjlim, but status objlimit\n", objval, lp->lpiobjlim);

               /* we want to resolve from the current basis (also if the LP had to be solved from scratch) */
               fromscratch = FALSE;

               /* temporarily disable cutoffbound, which also disables the objective limit */
               tmpcutoff = lp->cutoffbound;
               lp->cutoffbound = SCIPlpiInfinity(lpi);

               /* set lp pricing strategy to steepest edge */
               SCIP_CALL( SCIPsetGetCharParam(set, "lp/pricing", &tmppricingchar) );
               SCIP_CALL( SCIPsetSetCharParam(set, messagehdlr, "lp/pricing", 's') );

               /* resolve LP with an iteration limit of 1 */
               SCIP_CALL( lpSolve(lp, set, messagehdlr, stat, prob, SCIP_LPALGO_DUALSIMPLEX, 1, 1,
                     FALSE, FALSE, TRUE, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );

               /* reinstall old cutoff bound and lp pricing strategy */
               lp->cutoffbound = tmpcutoff;
               SCIP_CALL( SCIPsetSetCharParam(set, messagehdlr, "lp/pricing", tmppricingchar) );

               /* get objective value */
               SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );

               /* get solution status for the lp */
               solstat = SCIPlpGetSolstat(lp);
               assert(solstat != SCIP_LPSOLSTAT_OBJLIMIT);

               if( !(*lperror) && solstat != SCIP_LPSOLSTAT_ERROR && solstat != SCIP_LPSOLSTAT_NOTSOLVED )
               {
                  SCIPsetDebugMsg(set, " ---> new objval = %f (solstat: %d, 1 add. step)\n", objval, solstat);
               }

               /* disable fastmip for subsequent LP calls (if objective limit is not yet exceeded or LP solution is infeasible) */
               fastmip = 0;

               /* the solution is still not exceeding the objective limit and the solving process
                * was stopped due to time or iteration limit, solve again with fastmip turned off
                */
               if( solstat == SCIP_LPSOLSTAT_ITERLIMIT &&
                  SCIPsetIsLT(set, objval, lp->cutoffbound - getFiniteLooseObjval(lp, set, prob)) )
               {
                  assert(!(*lperror));

                  SCIP_CALL( lpSolve(lp, set, messagehdlr, stat, prob, SCIP_LPALGO_DUALSIMPLEX, -1, -1,
                        FALSE, FALSE, TRUE, fastmip, tightprimfeastol, tightdualfeastol, fromscratch, keepsol, lperror) );

                  /* get objective value */
                  SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );

                  /* get solution status for the lp */
                  solstat = SCIPlpGetSolstat(lp);

                  SCIPsetDebugMsg(set, " ---> new objval = %f (solstat: %d, without fastmip)\n", objval, solstat);
               }

               /* check for lp errors */
               if( *lperror || solstat == SCIP_LPSOLSTAT_ERROR || solstat == SCIP_LPSOLSTAT_NOTSOLVED )
               {
                  SCIPsetDebugMsg(set, "unresolved error while resolving LP in order to exceed the objlimit\n");
                  lp->solved = FALSE;
                  lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;

                  retcode = *lperror ? SCIP_OKAY : SCIP_LPERROR;
                  goto TERMINATE;
               }

               lp->solved = TRUE;

               /* optimal solution / objlimit with fastmip turned off / itlimit or timelimit, but objlimit exceeded */
               if( solstat == SCIP_LPSOLSTAT_OPTIMAL || solstat == SCIP_LPSOLSTAT_OBJLIMIT
                  || ( (solstat == SCIP_LPSOLSTAT_ITERLIMIT || solstat == SCIP_LPSOLSTAT_TIMELIMIT)
                     && SCIPsetIsGE(set, objval, lp->cutoffbound - getFiniteLooseObjval(lp, set, prob)) ) )
               {
                  /* get LP solution and possibly check the solution's feasibility again */
                  if( set->lp_checkprimfeas )
                  {
                     primalfeaspointer = &primalfeasible;
                     lp->primalchecked = TRUE;
                  }
                  else
                  {
                     /* believe in the primal feasibility of the LP solution */
                     primalfeasible = TRUE;
                     primalfeaspointer = NULL;
                     lp->primalchecked = FALSE;
                  }
                  if( set->lp_checkdualfeas )
                  {
                     dualfeaspointer = &dualfeasible;
                     lp->dualchecked = TRUE;
                  }
                  else
                  {
                     /* believe in the dual feasibility of the LP solution */
                     dualfeasible = TRUE;
                     dualfeaspointer = NULL;
                     lp->dualchecked = FALSE;
                  }

                  SCIP_CALL( SCIPlpGetSol(lp, set, stat, primalfeaspointer, dualfeaspointer) );

                  /* in debug mode, check that lazy bounds (if present) are not violated by an optimal LP solution */
                  if( solstat == SCIP_LPSOLSTAT_OPTIMAL )
                  {
                     checkLazyBounds(lp, set);
                  }

                  /* if objective value is larger than the cutoff bound, set solution status to objective
                   * limit reached and objective value to infinity, in case solstat = SCIP_LPSOLSTAT_OBJLIMIT,
                   * this was already done in the lpSolve() method
                   */
                  if( SCIPsetIsGE(set, objval, lp->cutoffbound - getFiniteLooseObjval(lp, set, prob)) )
                  {
                     lp->lpsolstat = SCIP_LPSOLSTAT_OBJLIMIT;
                     lp->lpobjval = SCIPsetInfinity(set);
                  }

                  /* LP solution is not feasible or objective limit was reached without the LP value really exceeding
                   * the cutoffbound; mark the LP to be unsolved
                   */
                  if( !primalfeasible || !dualfeasible
                     || (solstat == SCIP_LPSOLSTAT_OBJLIMIT &&
                        !SCIPsetIsGE(set, objval, lp->cutoffbound - getFiniteLooseObjval(lp, set, prob))) )
                  {
                     lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_HIGH, "unresolved");
                     lp->solved = FALSE;
                     lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
                     *lperror = TRUE;
                  }

                  SCIPsetDebugMsg(set, " -> LP objective value: %g + %g = %g (solstat=%d, cutoff=%g)\n",
                     lp->lpobjval, getFiniteLooseObjval(lp, set, prob), lp->lpobjval + getFiniteLooseObjval(lp, set, prob),
                     lp->lpsolstat, lp->cutoffbound);
               }
               /* infeasible solution */
               else if( solstat == SCIP_LPSOLSTAT_INFEASIBLE )
               {
                  SCIP_CALL( SCIPlpGetDualfarkas(lp, set, stat) );
                  SCIPsetDebugMsg(set, " -> LP infeasible\n");
               }
               /* unbounded solution */
               else if( solstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
               {
                  if( set->lp_checkprimfeas )
                  {
                     /* get unbounded LP solution and check the solution's feasibility again */
                     SCIP_CALL( SCIPlpGetUnboundedSol(lp, set, stat, &primalfeasible, &rayfeasible) );

                     lp->primalchecked = TRUE;
                  }
                  else
                  {
                     /* get unbounded LP solution believing in its feasibility */
                     SCIP_CALL( SCIPlpGetUnboundedSol(lp, set, stat, NULL, NULL) );

                     primalfeasible = TRUE;
                     rayfeasible = TRUE;
                     lp->primalchecked = FALSE;
                  }

                  SCIPsetDebugMsg(set, " -> LP has unbounded primal ray\n");

                  /* in debug mode, check that lazy bounds (if present) are not violated */
                  checkLazyBounds(lp, set);

                  if( !primalfeasible || !rayfeasible )
                  {
                     /* unbounded solution is infeasible (this can happen due to numerical problems):
                      * forget about the LP at this node and mark it to be unsolved

                      * @todo: like in the default LP solving evaluation, solve without fastmip,
                      * with tighter feasibility tolerance and from scratch
                      */
                     lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "unresolved, unbounded LP");
                     lp->solved = FALSE;
                     lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
                     *lperror = TRUE;
                  }

               }

               assert(lp->lpsolstat != SCIP_LPSOLSTAT_ITERLIMIT);
               assert(SCIPsetIsGE(set, objval, lp->cutoffbound - getFiniteLooseObjval(lp, set, prob))
                  || lp->lpsolstat != SCIP_LPSOLSTAT_OBJLIMIT);
            }
            else
            {
               SCIP_CALL( SCIPlpGetSol(lp, set, stat, NULL, NULL) );
            }
         }
         else if( !SCIPprobAllColsInLP(prob, set, lp) || set->misc_exactsolve )
         {
            SCIP_CALL( SCIPlpGetSol(lp, set, stat, NULL, NULL) );
         }
         SCIPsetDebugMsg(set, " -> LP objective limit reached\n");
         break;

      case SCIP_LPSOLSTAT_ITERLIMIT:
         SCIPsetDebugMsg(set, " -> LP iteration limit exceeded\n");
         break;

      case SCIP_LPSOLSTAT_TIMELIMIT:
         SCIPsetDebugMsg(set, " -> LP time limit exceeded\n");

         /* make sure that we evaluate the time limit exactly in order to avoid erroneous warning */
         stat->nclockskipsleft = 0;
         if( !SCIPsolveIsStopped(set, stat, FALSE) )
         {
            SCIPmessagePrintWarning(messagehdlr, "LP solver reached time limit, but SCIP time limit is not exceeded yet; "
               "you might consider switching the clock type of SCIP\n");
            stat->status = SCIP_STATUS_TIMELIMIT;
         }
         break;

      case SCIP_LPSOLSTAT_ERROR:
      case SCIP_LPSOLSTAT_NOTSOLVED:
         SCIPerrorMessage("error in LP solver\n");
         retcode = SCIP_LPERROR;
         goto TERMINATE;

      default:
         SCIPerrorMessage("unknown LP solution status\n");
         retcode = SCIP_ERROR;
         goto TERMINATE;
      }
   }
   assert(!(*lperror) || !lp->solved);

 TERMINATE:
   /* if the LP had to be solved from scratch, we have to reset this flag since it is stored in the LPI; otherwise it
    * may happen that we continue to solve from scratch during strong branching */
   if( lp->lpifromscratch )
   {
      SCIP_Bool success;
      (void) lpSetFromscratch(lp, FALSE, &success);
      SCIPsetDebugMsg(set, "resetting parameter SCIP_LPPARAM_FROMSCRATCH to FALSE %s\n", success ? "" : "failed");
   }

   return retcode;
}

/** gets solution status of current LP */
SCIP_LPSOLSTAT SCIPlpGetSolstat(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved || lp->lpsolstat == SCIP_LPSOLSTAT_NOTSOLVED);

   return (lp->flushed ? lp->lpsolstat : SCIP_LPSOLSTAT_NOTSOLVED);
}

/** gets objective value of current LP
 *
 *  @note This method returns the objective value of the current LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the LP solution status is
 *        SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 */
SCIP_Real SCIPlpGetObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert((lp->nloosevars > 0) || (lp->looseobjvalinf == 0 && lp->looseobjval == 0.0));
   assert(set != NULL);

   if( !lp->flushed )
      return SCIP_INVALID;
   else if( SCIPsetIsInfinity(set, lp->lpobjval) || SCIPsetIsInfinity(set, -lp->lpobjval))
      return lp->lpobjval;
   else if( lp->looseobjvalinf > 0 )
      return -SCIPsetInfinity(set);
   else
   {
      /* recalculate the loose objective value, if needed */
      if( !lp->looseobjvalid )
         recomputeLooseObjectiveValue(lp, set, prob);

      return lp->lpobjval + lp->looseobjval;
   }
}

/** gets part of objective value of current LP that results from COLUMN variables only */
SCIP_Real SCIPlpGetColumnObjval(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->solved);

   return (lp->flushed ? lp->lpobjval : SCIP_INVALID);
}

/** gets part of objective value of current LP that results from LOOSE variables only */
SCIP_Real SCIPlpGetLooseObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert((lp->nloosevars > 0) || (lp->looseobjvalinf == 0 && lp->looseobjval == 0.0));
   assert(set != NULL);

   if( !lp->flushed )
      return SCIP_INVALID;
   else if( lp->looseobjvalinf > 0 )
      return -SCIPsetInfinity(set);
   else
      return getFiniteLooseObjval(lp, set, prob);
}

/** remembers the current LP objective value as root solution value */
void SCIPlpStoreRootObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(lp != NULL);

   lp->rootlpobjval = SCIPlpGetColumnObjval(lp);
   lp->rootlooseobjval = SCIPlpGetLooseObjval(lp, set, prob);
}

/** invalidates the root LP solution value */
void SCIPlpInvalidateRootObjval(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   lp->rootlpobjval = SCIP_INVALID;
   lp->rootlooseobjval = SCIP_INVALID;
}

/** recomputes local and global pseudo objective values */
void SCIPlpRecomputeLocalAndGlobalPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(lp != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   vars = prob->vars;
   nvars = prob->nvars;

   lp->glbpseudoobjvalinf = 0;
   lp->glbpseudoobjval = 0.0;

   lp->pseudoobjvalinf = 0;
   lp->pseudoobjval = 0.0;

   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real obj = SCIPvarGetObj(vars[v]);

      if( SCIPsetIsPositive(set, obj) )
      {
         /* update the global pseudo objective value */
         if( SCIPsetIsInfinity(set, -SCIPvarGetLbGlobal(vars[v])) )
            ++(lp->glbpseudoobjvalinf);
         else
            lp->glbpseudoobjval += obj * SCIPvarGetLbGlobal(vars[v]);

         /* update the local pseudo objective value */
         if( SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(vars[v])) )
            ++(lp->pseudoobjvalinf);
         else
            lp->pseudoobjval += obj * SCIPvarGetLbLocal(vars[v]);
      }

      if( SCIPsetIsNegative(set, obj) )
      {
         /* update the global pseudo objective value */
         if( SCIPsetIsInfinity(set, SCIPvarGetUbGlobal(vars[v])) )
            ++(lp->glbpseudoobjvalinf);
         else
            lp->glbpseudoobjval += obj * SCIPvarGetUbGlobal(vars[v]);

         /* update the local pseudo objective value */
         if( SCIPsetIsInfinity(set, SCIPvarGetUbLocal(vars[v])) )
            ++(lp->pseudoobjvalinf);
         else
            lp->pseudoobjval += obj * SCIPvarGetUbLocal(vars[v]);
      }
   }

   /* the recomputed values are reliable */
   lp->relglbpseudoobjval = lp->glbpseudoobjval;
   lp->glbpseudoobjvalid = TRUE;
   lp->relpseudoobjval = lp->pseudoobjval;
   lp->pseudoobjvalid = TRUE;
}

/** gets the global pseudo objective value; that is all variables set to their best (w.r.t. the objective function)
 *  global bound
 */
SCIP_Real SCIPlpGetGlobalPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(lp != NULL);
   assert(lp->glbpseudoobjvalinf >= 0);
   assert(set != NULL);

   if( lp->glbpseudoobjvalinf > 0 || set->nactivepricers > 0 )
      return -SCIPsetInfinity(set);
   else
   {
      /* recalculate the global pseudo solution value, if needed */
      if( !lp->glbpseudoobjvalid )
         recomputeGlbPseudoObjectiveValue(lp, set, prob);

      /* if the global pseudo objective value is smaller than -infinity, we just return -infinity */
      if( SCIPsetIsInfinity(set, -lp->glbpseudoobjval) )
         return -SCIPsetInfinity(set);

      if( SCIPsetIsInfinity(set, lp->glbpseudoobjval) )
         return SCIPsetInfinity(set);

      return lp->glbpseudoobjval;
   }
}

/** gets the pseudo objective value for the current search node; that is all variables set to their best (w.r.t. the
 *  objective function) local bound
 */
SCIP_Real SCIPlpGetPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(lp != NULL);
   assert(lp->pseudoobjvalinf >= 0);
   assert(set != NULL);

   if( lp->pseudoobjvalinf > 0 || set->nactivepricers > 0 )
      return -SCIPsetInfinity(set);
   else
   {
      /* recalculate the pseudo solution value, if needed */
      if( !lp->pseudoobjvalid )
         recomputePseudoObjectiveValue(lp, set, prob);

      /* if the pseudo objective value is smaller than -infinity, we just return -infinity */
      if( SCIPsetIsInfinity(set, -lp->pseudoobjval) )
         return -SCIPsetInfinity(set);

      if( SCIPsetIsInfinity(set, lp->pseudoobjval) )
         return SCIPsetInfinity(set);

      return lp->pseudoobjval;
   }
}

/** gets pseudo objective value, if a bound of the given variable would be modified in the given way */
SCIP_Real SCIPlpGetModifiedPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   SCIP_Real pseudoobjval;
   int pseudoobjvalinf;
   SCIP_Real obj;

   pseudoobjval = getFinitePseudoObjval(lp, set, prob);
   pseudoobjvalinf = lp->pseudoobjvalinf;
   obj = SCIPvarGetObj(var);
   if( !SCIPsetIsZero(set, obj) && boundtype == SCIPvarGetBestBoundType(var) )
   {
      if( SCIPsetIsInfinity(set, REALABS(oldbound)) )
         pseudoobjvalinf--;
      else
         pseudoobjval -= oldbound * obj;
      assert(pseudoobjvalinf >= 0);
      if( SCIPsetIsInfinity(set, REALABS(newbound)) )
         pseudoobjvalinf++;
      else
         pseudoobjval += newbound * obj;
   }
   assert(pseudoobjvalinf >= 0);

   if( pseudoobjvalinf > 0 || set->nactivepricers > 0 )
      return -SCIPsetInfinity(set);
   else
      return pseudoobjval;
}

/** gets pseudo objective value, if a bound of the given variable would be modified in the given way;
 *  perform calculations with interval arithmetic to get an exact lower bound
 */
SCIP_Real SCIPlpGetModifiedProvedPseudoObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   SCIP_Real pseudoobjval;
   int pseudoobjvalinf;
   SCIP_Real obj;

   assert(lp->pseudoobjvalid);

   pseudoobjval = lp->pseudoobjval;
   pseudoobjvalinf = lp->pseudoobjvalinf;
   obj = SCIPvarGetObj(var);
   if( !SCIPsetIsZero(set, obj) && boundtype == SCIPvarGetBestBoundType(var) )
   {
      SCIP_INTERVAL objint;
      SCIP_INTERVAL bd;
      SCIP_INTERVAL prod;
      SCIP_INTERVAL psval;

      SCIPintervalSet(&psval, pseudoobjval);
      SCIPintervalSet(&objint, SCIPvarGetObj(var));

      if( SCIPsetIsInfinity(set, REALABS(oldbound)) )
         pseudoobjvalinf--;
      else
      {
         SCIPintervalSet(&bd, oldbound);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, objint);
         SCIPintervalSub(SCIPsetInfinity(set), &psval, psval, prod);
      }
      assert(pseudoobjvalinf >= 0);
      if( SCIPsetIsInfinity(set, REALABS(newbound)) )
         pseudoobjvalinf++;
      else
      {
         SCIPintervalSet(&bd, newbound);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, objint);
         SCIPintervalAdd(SCIPsetInfinity(set), &psval, psval, prod);
      }

      pseudoobjval = SCIPintervalGetInf(psval);
   }
   assert(pseudoobjvalinf >= 0);

   if( pseudoobjvalinf > 0 || set->nactivepricers > 0 )
      return -SCIPsetInfinity(set);
   else
      return pseudoobjval;
}

/** compute the objective delta due the new objective coefficient */
static
void getObjvalDeltaObj(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             oldobj,             /**< old objective value of variable */
   SCIP_Real             newobj,             /**< new objective value of variable */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real*            deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   assert(!SCIPsetIsInfinity(set, REALABS(oldobj)));
   assert(!SCIPsetIsInfinity(set, REALABS(newobj)));
   assert(!SCIPsetIsInfinity(set, lb));
   assert(!SCIPsetIsInfinity(set, -ub));
   assert(!SCIPsetIsEQ(set, oldobj, newobj));

   (*deltaval) = 0.0;
   (*deltainf) = 0;

   if( SCIPsetIsPositive(set, oldobj) )
   {
      /* sign of objective did not change */
      if( SCIPsetIsPositive(set, newobj) )
      {
         /* if the bound is finite, calculate the deltaval */
         if( !SCIPsetIsInfinity(set, -lb) )
            (*deltaval) = lb * (newobj - oldobj);
      }
      /* sign of objective did change, so the best bound does change */
      else if( SCIPsetIsNegative(set, newobj) )
      {
         if( SCIPsetIsInfinity(set, -lb) )
         {
            /* old best bound was infinite while new one is not */
            if( !SCIPsetIsInfinity(set, ub) )
            {
               (*deltainf) = -1;
               (*deltaval) = ub * newobj;
            }
         }
         else
         {
            /* new best bound is infinite while old one was not */
            if( SCIPsetIsInfinity(set, ub) )
            {
               (*deltainf) = 1;
               (*deltaval) = -lb * oldobj;
            }
            /* neither old nor new best bound is infinite, so just calculate the deltaval */
            else
            {
               (*deltaval) = (ub * newobj) - (lb * oldobj);
            }
         }
      }
      /* new objective is 0.0 */
      else
      {
         if( SCIPsetIsInfinity(set, -lb) )
            (*deltainf) = -1;
         else
            (*deltaval) = -lb * oldobj;
      }

   }
   else if( SCIPsetIsNegative(set, oldobj) )
   {
      /* sign of objective did not change */
      if( SCIPsetIsNegative(set, newobj) )
      {
         /* if the bound is finite, calculate the deltaval */
         if( !SCIPsetIsInfinity(set, ub) )
            (*deltaval) = ub * (newobj - oldobj);
      }
      /* sign of objective did change, so the best bound does change */
      else if( SCIPsetIsPositive(set, newobj) )
      {
         if( SCIPsetIsInfinity(set, ub) )
         {
            /* old best bound was infinite while new one is not */
            if( !SCIPsetIsInfinity(set, -lb) )
            {
               (*deltainf) = -1;
               (*deltaval) = lb * newobj;
            }
         }
         else
         {
            /* new best bound is infinite while old one was not */
            if( SCIPsetIsInfinity(set, -lb) )
            {
               (*deltainf) = 1;
               (*deltaval) = -ub * oldobj;
            }
            /* neither old nor new best bound is infinite, so just calculate the deltaval */
            else
            {
               (*deltaval) = (lb * newobj) - (ub * oldobj);
            }
         }
      }
      /* new objective is 0.0 */
      else
      {
         if( SCIPsetIsInfinity(set, ub) )
            (*deltainf) = -1;
         else
            (*deltaval) = -ub * oldobj;
      }
   }
   /* old objective was 0.0 */
   else
   {
      if( SCIPsetIsNegative(set, newobj) )
      {
         if( SCIPsetIsInfinity(set, ub) )
            (*deltainf) = 1;
         else
            (*deltaval) = ub * newobj;
      }
      else if( SCIPsetIsPositive(set, newobj) )
      {
         if( SCIPsetIsInfinity(set, -lb) )
            (*deltainf) = 1;
         else
            (*deltaval) = lb * newobj;
      }
   }
}

/** compute the objective delta due the new lower bound */
static
void getObjvalDeltaLb(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             obj,                /**< objective value of variable */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb,              /**< new lower bound of variable */
   SCIP_Real*            deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   assert(!SCIPsetIsInfinity(set, REALABS(obj)));
   assert(!SCIPsetIsInfinity(set, oldlb));
   assert(!SCIPsetIsInfinity(set, -oldlb) || !SCIPsetIsInfinity(set, -newlb));
   assert(SCIPsetIsPositive(set, obj)); /* we only need to update if the objective is positive */

   if( SCIPsetIsInfinity(set, -oldlb) )
   {
      if( !SCIPsetIsInfinity(set, newlb) )
      {
         (*deltainf) = -1;
         (*deltaval) = newlb * obj;
      }
      else
      {
         (*deltainf) = 0;
         (*deltaval) = 0.0;
      }
   }
   else if( SCIPsetIsInfinity(set, REALABS(newlb)) )
   {
      (*deltainf) = 1;
      (*deltaval) = -oldlb * obj;
   }
   else
   {
      (*deltainf) = 0;
      (*deltaval) = obj * (newlb - oldlb);
   }
}

/** compute the objective delta due the new upper bound */
static
void getObjvalDeltaUb(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             obj,                /**< objective value of variable */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub,              /**< new upper bound of variable */
   SCIP_Real*            deltaval,           /**< pointer to store the delta value */
   int*                  deltainf            /**< pointer to store the number of variables with infinite best bound */
   )
{
   assert(!SCIPsetIsInfinity(set, REALABS(obj)));
   assert(!SCIPsetIsInfinity(set, -oldub));
   assert(!SCIPsetIsInfinity(set, oldub) || !SCIPsetIsInfinity(set, newub));
   assert(SCIPsetIsNegative(set, obj)); /* we only need to update if the objective is negative */

   if( SCIPsetIsInfinity(set, oldub) )
   {
      if( !SCIPsetIsInfinity(set, -newub) )
      {
         (*deltainf) = -1;
         (*deltaval) = newub * obj;
      }
      else
      {
         (*deltainf) = 0;
         (*deltaval) = 0.0;
      }
   }
   else if( SCIPsetIsInfinity(set, REALABS(newub)) )
   {
      (*deltainf) = 1;
      (*deltaval) = -oldub * obj;
   }
   else
   {
      (*deltainf) = 0;
      (*deltaval) = obj * (newub - oldub);
   }
}

/** updates current pseudo and loose objective values for a change in a variable's objective value or bounds */
static
void lpUpdateObjval(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             deltaval,           /**< delta value in the objective function */
   int                   deltainf,           /**< delta value for the number of variables with infinite best bound */
   SCIP_Bool             local,              /**< should the local pseudo objective value be updated? */
   SCIP_Bool             loose,              /**< should the loose objective value be updated? */
   SCIP_Bool             global              /**< should the global pseudo objective value be updated? */
   )
{
   assert(lp != NULL);
   assert(lp->looseobjvalinf >= 0);
   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->glbpseudoobjvalinf >= 0);

   /* update the pseudo objective value */
   if( local )
   {
      lp->pseudoobjvalinf += deltainf;
      if( lp->pseudoobjvalid )
      {
         lp->pseudoobjval += deltaval;

         /* if the absolute value was increased, this is regarded as reliable,
          * otherwise, we check whether we can still trust the updated value
          */
         if( REALABS(lp->relpseudoobjval) < REALABS(lp->pseudoobjval) )
            lp->relpseudoobjval = lp->pseudoobjval;
         else if( SCIPsetIsUpdateUnreliable(set, lp->pseudoobjval, lp->relpseudoobjval) )
            lp->pseudoobjvalid = FALSE;
      }

      /* after changing a local bound on a LOOSE variable, we have to update the loose objective value, too */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
         loose = TRUE;
   }
   /* update the loose objective value */
   if( loose )
   {
      lp->looseobjvalinf += deltainf;

      if( deltaval != 0.0 && lp->looseobjvalid )
      {
         lp->looseobjval += deltaval;

         /* if the absolute value was increased, this is regarded as reliable,
          * otherwise, we check whether we can still trust the updated value
          */
         if( REALABS(lp->rellooseobjval) < REALABS(lp->looseobjval) )
            lp->rellooseobjval = lp->looseobjval;
         else if( SCIPsetIsUpdateUnreliable(set, lp->looseobjval, lp->rellooseobjval) )
            lp->looseobjvalid = FALSE;
      }
   }
   /* update the root pseudo objective values */
   if( global )
   {
      lp->glbpseudoobjvalinf += deltainf;
      if( lp->glbpseudoobjvalid )
      {
         lp->glbpseudoobjval += deltaval;

         /* if the absolute value was increased, this is regarded as reliable,
          * otherwise, we check whether we can still trust the updated value
          */
         if( REALABS(lp->relglbpseudoobjval) < REALABS(lp->glbpseudoobjval) )
            lp->relglbpseudoobjval = lp->glbpseudoobjval;
         else if( SCIPsetIsUpdateUnreliable(set, lp->glbpseudoobjval, lp->relglbpseudoobjval) )
            lp->glbpseudoobjvalid = FALSE;
      }
   }

   assert(lp->looseobjvalinf >= 0);
   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->glbpseudoobjvalinf >= 0);
}

/** updates current pseudo and loose objective values for a change in a variable's objective value or bounds;
 *  pseudo objective value is calculated with interval arithmetics to get a proved lower bound
 */
static
SCIP_RETCODE lpUpdateVarProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldobj,             /**< old objective value of variable */
   SCIP_Real             oldlb,              /**< old objective value of variable */
   SCIP_Real             oldub,              /**< old objective value of variable */
   SCIP_Real             newobj,             /**< new objective value of variable */
   SCIP_Real             newlb,              /**< new objective value of variable */
   SCIP_Real             newub               /**< new objective value of variable */
   )
{
   SCIP_INTERVAL deltaval;
   SCIP_INTERVAL bd;
   SCIP_INTERVAL obj;
   SCIP_INTERVAL prod;
   SCIP_INTERVAL psval;
   int deltainf;

   assert(lp != NULL);
   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->looseobjvalinf >= 0);
   assert(!SCIPsetIsInfinity(set, REALABS(oldobj)));
   assert(!SCIPsetIsInfinity(set, oldlb));
   assert(!SCIPsetIsInfinity(set, -oldub));
   assert(!SCIPsetIsInfinity(set, REALABS(newobj)));
   assert(!SCIPsetIsInfinity(set, newlb));
   assert(!SCIPsetIsInfinity(set, -newub));
   assert(var != NULL);

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("LP was informed of an objective change of a non-active variable\n");
      return SCIP_INVALIDDATA;
   }

   assert(SCIPvarGetProbindex(var) >= 0);

   SCIPintervalSet(&deltaval, 0.0);
   deltainf = 0;

   /* subtract old pseudo objective value */
   if( oldobj > 0.0 )
   {
      if( SCIPsetIsInfinity(set, -oldlb) )
         deltainf--;
      else
      {
         SCIPintervalSet(&bd, oldlb);
         SCIPintervalSet(&obj, oldobj);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, obj);
         SCIPintervalSub(SCIPsetInfinity(set), &deltaval, deltaval, prod);  /* deltaval -= oldlb * oldobj; */
      }
   }
   else if( oldobj < 0.0 )
   {
      if( SCIPsetIsInfinity(set, oldub) )
         deltainf--;
      else
      {
         SCIPintervalSet(&bd, oldub);
         SCIPintervalSet(&obj, oldobj);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, obj);
         SCIPintervalSub(SCIPsetInfinity(set), &deltaval, deltaval, prod);  /* deltaval -= oldub * oldobj; */
      }
   }

   /* add new pseudo objective value */
   if( newobj > 0.0 )
   {
      if( SCIPsetIsInfinity(set, -newlb) )
         deltainf++;
      else
      {
         SCIPintervalSet(&bd, newlb);
         SCIPintervalSet(&obj, newobj);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, obj);
         SCIPintervalAdd(SCIPsetInfinity(set), &deltaval, deltaval, prod);  /* deltaval += newlb * newobj; */
      }
   }
   else if( newobj < 0.0 )
   {
      if( SCIPsetIsInfinity(set, newub) )
         deltainf++;
      else
      {
         SCIPintervalSet(&bd, newub);
         SCIPintervalSet(&obj, newobj);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, obj);
         SCIPintervalAdd(SCIPsetInfinity(set), &deltaval, deltaval, prod);  /* deltaval += newub * newobj; */
      }
   }

   /* update the pseudo and loose objective values */
   SCIPintervalSet(&psval, lp->pseudoobjval);
   SCIPintervalAdd(SCIPsetInfinity(set), &psval, psval, deltaval);
   lp->pseudoobjval = SCIPintervalGetInf(psval);
   lp->pseudoobjvalinf += deltainf;
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIPintervalSet(&psval, lp->looseobjval);
      SCIPintervalAdd(SCIPsetInfinity(set), &psval, psval, deltaval);
      lp->looseobjval = SCIPintervalGetInf(psval);
      lp->looseobjvalinf += deltainf;
   }

   assert(lp->pseudoobjvalinf >= 0);
   assert(lp->looseobjvalinf >= 0);

   return SCIP_OKAY;
}

/** updates current pseudo and loose objective value for a change in a variable's objective value */
SCIP_RETCODE SCIPlpUpdateVarObj(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldobj,             /**< old objective value of variable */
   SCIP_Real             newobj              /**< new objective value of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldobj != newobj ) /*lint !e777*/
      {
         SCIP_CALL( lpUpdateVarProved(lp, set, var, oldobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               newobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldobj, newobj) )
      {
         SCIP_Real deltaval;
         int deltainf;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetProbindex(var) >= 0);

         /* the objective coefficient can only be changed during presolving, that implies that the global and local
          * domain of the variable are the same
          */
         assert(lp->probing || SCIPsetIsEQ(set, SCIPvarGetLbGlobal(var), SCIPvarGetLbLocal(var)));
         assert(lp->probing || SCIPsetIsEQ(set, SCIPvarGetUbGlobal(var), SCIPvarGetUbLocal(var)));

         /* compute the pseudo objective delta due the new objective coefficient */
         getObjvalDeltaObj(set, oldobj, newobj, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), &deltaval, &deltainf);

         /* update the local pseudo objective value */
         lpUpdateObjval(lp, set, var, deltaval, deltainf, TRUE, FALSE, FALSE);

         /* compute the pseudo objective delta due the new objective coefficient */
         getObjvalDeltaObj(set, oldobj, newobj, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), &deltaval, &deltainf);

         /* update the global pseudo objective value */
         lpUpdateObjval(lp, set, var, deltaval, deltainf, FALSE, FALSE, TRUE);
      }
   }

   return SCIP_OKAY;
}


/** updates current root pseudo objective value for a global change in a variable's lower bound */
SCIP_RETCODE SCIPlpUpdateVarLbGlobal(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb               /**< new lower bound of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( !SCIPsetIsEQ(set, oldlb, newlb) && SCIPsetIsPositive(set, SCIPvarGetObj(var)) )
   {
      SCIP_Real deltaval;
      int deltainf;

      /* compute the pseudo objective delta due the new lower bound */
      getObjvalDeltaLb(set, SCIPvarGetObj(var), oldlb, newlb, &deltaval, &deltainf);

      /* update the root pseudo objective values */
      lpUpdateObjval(lp, set, var, deltaval, deltainf, FALSE, FALSE, TRUE);

   }

   return SCIP_OKAY;
}

/** updates current pseudo and loose objective value for a change in a variable's lower bound */
SCIP_RETCODE SCIPlpUpdateVarLb(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb               /**< new lower bound of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldlb != newlb && SCIPvarGetObj(var) > 0.0 ) /*lint !e777*/
      {
         SCIP_CALL( lpUpdateVarProved(lp, set, var, SCIPvarGetObj(var), oldlb, SCIPvarGetUbLocal(var), 
               SCIPvarGetObj(var), newlb, SCIPvarGetUbLocal(var)) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldlb, newlb) && SCIPsetIsPositive(set, SCIPvarGetObj(var)) )
      {
         SCIP_Real deltaval;
         int deltainf;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetProbindex(var) >= 0);

         /* compute the pseudo objective delta due the new lower bound */
         getObjvalDeltaLb(set, SCIPvarGetObj(var), oldlb, newlb, &deltaval, &deltainf);

         /* update the pseudo and loose objective values */
         lpUpdateObjval(lp, set, var, deltaval, deltainf, TRUE, FALSE, FALSE);
      }
   }

   return SCIP_OKAY;
}

/** updates current root pseudo objective value for a global change in a variable's upper bound */
SCIP_RETCODE SCIPlpUpdateVarUbGlobal(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub               /**< new upper bound of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( !SCIPsetIsEQ(set, oldub, newub) && SCIPsetIsNegative(set, SCIPvarGetObj(var)) )
   {
      SCIP_Real deltaval;
      int deltainf;

      /* compute the pseudo objective delta due the new upper bound */
      getObjvalDeltaUb(set, SCIPvarGetObj(var), oldub, newub, &deltaval, &deltainf);

      /* update the root pseudo objective values */
      lpUpdateObjval(lp, set, var, deltaval, deltainf, FALSE, FALSE, TRUE);
   }

   return SCIP_OKAY;
}

/** updates current pseudo objective value for a change in a variable's upper bound */
SCIP_RETCODE SCIPlpUpdateVarUb(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable that changed */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub               /**< new upper bound of variable */
   )
{
   assert(set != NULL);
   assert(var != NULL);

   if( set->misc_exactsolve )
   {
      if( oldub != newub && SCIPvarGetObj(var) < 0.0 ) /*lint !e777*/
      {
         SCIP_CALL( lpUpdateVarProved(lp, set, var, SCIPvarGetObj(var), SCIPvarGetLbLocal(var), oldub, 
               SCIPvarGetObj(var), SCIPvarGetLbLocal(var), newub) );
      }
   }
   else
   {
      if( !SCIPsetIsEQ(set, oldub, newub) && SCIPsetIsNegative(set, SCIPvarGetObj(var)) )
      {
         SCIP_Real deltaval;
         int deltainf;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetProbindex(var) >= 0);

         /* compute the pseudo objective delta due the new upper bound */
         getObjvalDeltaUb(set, SCIPvarGetObj(var), oldub, newub, &deltaval, &deltainf);

         /* update the pseudo and loose objective values */
         lpUpdateObjval(lp, set, var, deltaval, deltainf, TRUE, FALSE, FALSE);
      }
   }

   return SCIP_OKAY;
}

/** informs LP, that given variable was added to the problem */
SCIP_RETCODE SCIPlpUpdateAddVar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that is now a LOOSE problem variable */
   )
{
   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   /* add the variable to the loose objective value sum */
   SCIP_CALL( SCIPlpUpdateVarObj(lp, set, var, 0.0, SCIPvarGetObj(var)) );

   /* update the loose variables counter */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      lp->nloosevars++;

   return SCIP_OKAY;
}

/** informs LP, that given variable is to be deleted from the problem */
SCIP_RETCODE SCIPlpUpdateDelVar(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable that will be deleted from the problem */
   )
{
   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   /* subtract the variable from the loose objective value sum */
   SCIP_CALL( SCIPlpUpdateVarObj(lp, set, var, SCIPvarGetObj(var), 0.0) );

   /* update the loose variables counter */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
   {
      SCIPlpDecNLoosevars(lp);
   }

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable */
static
SCIP_RETCODE lpUpdateVarColumn(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(lp != NULL);
   assert(lp->nloosevars > 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(lp->looseobjvalinf >= 0);

   obj = SCIPvarGetObj(var);

   /* update loose objective value */
   if( SCIPsetIsPositive(set, obj) )
   {
      lb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsInfinity(set, -lb) )
         lp->looseobjvalinf--;
      else
         lpUpdateObjval(lp, set, var, -lb * obj, 0, FALSE, TRUE, FALSE);
   }
   else if( SCIPsetIsNegative(set, obj) )
   {
      ub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsInfinity(set, ub) )
         lp->looseobjvalinf--;
      else
         lpUpdateObjval(lp, set, var, -ub * obj, 0, FALSE, TRUE, FALSE);
   }

   SCIPlpDecNLoosevars(lp);

   assert(lp->looseobjvalinf >= 0);

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable
 *  pseudo objective value is calculated with interval arithmetics to get a proved lower bound
 */
static
SCIP_RETCODE lpUpdateVarColumnProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   SCIP_INTERVAL bd;
   SCIP_INTERVAL ob;
   SCIP_INTERVAL prod;
   SCIP_INTERVAL loose;
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(lp != NULL);
   assert(lp->nloosevars > 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(SCIPvarGetProbindex(var) >= 0);

   obj = SCIPvarGetObj(var);

   SCIPintervalSet(&loose, lp->looseobjval);

   /* update loose objective value corresponding to the deletion of variable */
   if( obj > 0.0 )
   {
      lb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsInfinity(set, -lb) )
         lp->looseobjvalinf--;
      else
      {
         SCIPintervalSet(&bd, lb);
         SCIPintervalSet(&ob, obj);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, ob);
         SCIPintervalSub(SCIPsetInfinity(set), &loose, loose, prod);  /* lp->looseobjval -= lb * obj; */
      }
   }
   else if( SCIPsetIsNegative(set, obj) )
   {
      ub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsInfinity(set, ub) )
         lp->looseobjvalinf--;
      else
      {
         SCIPintervalSet(&bd, ub);
         SCIPintervalSet(&ob, obj);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, ob);
         SCIPintervalSub(SCIPsetInfinity(set), &loose, loose, prod);  /* lp->looseobjval -= ub * obj; */
      }
   }
   lp->nloosevars--;

   /* get rid of numerical problems: set loose objective value explicitly to zero, if no loose variables remain */
   if( lp->nloosevars == 0 )
   {
      assert(lp->looseobjvalinf == 0);
      lp->looseobjval = 0.0;
   }
   else
      lp->looseobjval = SCIPintervalGetInf(loose);

   return SCIP_OKAY;
}

/** informs LP, that given formerly loose problem variable is now a column variable */
SCIP_RETCODE SCIPlpUpdateVarColumn(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from LOOSE to COLUMN */
   )
{
   assert(set != NULL);

   if( set->misc_exactsolve )
   {
      SCIP_CALL( lpUpdateVarColumnProved(lp, set, var) );
   }
   else
   {
      SCIP_CALL( lpUpdateVarColumn(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable */
static
SCIP_RETCODE lpUpdateVarLoose(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(lp->looseobjvalinf >= 0);

   obj = SCIPvarGetObj(var);

   /* update loose objective value corresponding to the addition of variable */
   if( SCIPsetIsPositive(set, obj) )
   {
      lb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsInfinity(set, -lb) )
         lp->looseobjvalinf++;
      else
         lpUpdateObjval(lp, set, var, lb * obj, 0, FALSE, TRUE, FALSE);
   }
   else if( SCIPsetIsNegative(set, obj) )
   {
      ub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsInfinity(set, ub) )
         lp->looseobjvalinf++;
      else
         lpUpdateObjval(lp, set, var, ub * obj, 0, FALSE, TRUE, FALSE);
   }
   lp->nloosevars++;

   assert(lp->looseobjvalinf >= 0);

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable
 *  pseudo objective value is calculated with interval arithmetics to get a proved lower bound
 */
static
SCIP_RETCODE lpUpdateVarLooseProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   SCIP_INTERVAL bd;
   SCIP_INTERVAL ob;
   SCIP_INTERVAL prod;
   SCIP_INTERVAL loose;
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(lp != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);
   assert(SCIPvarGetProbindex(var) >= 0);

   obj = SCIPvarGetObj(var);

   SCIPintervalSet(&loose, lp->looseobjval);

   /* update loose objective value corresponding to the deletion of variable */
   if( obj > 0.0 )
   {
      lb = SCIPvarGetLbLocal(var);
      if( SCIPsetIsInfinity(set, -lb) )
         lp->looseobjvalinf++;
      else
      {
         SCIPintervalSet(&bd, lb);
         SCIPintervalSet(&ob, obj);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, ob);
         SCIPintervalAdd(SCIPsetInfinity(set), &loose, loose, prod);  /* lp->looseobjval += lb * obj; */
      }
   }
   else if( SCIPsetIsNegative(set, obj) )
   {
      ub = SCIPvarGetUbLocal(var);
      if( SCIPsetIsInfinity(set, ub) )
         lp->looseobjvalinf++;
      else
      {
         SCIPintervalSet(&bd, ub);
         SCIPintervalSet(&ob, obj);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, bd, ob);
         SCIPintervalAdd(SCIPsetInfinity(set), &loose, loose, prod);  /* lp->looseobjval += ub * obj; */
      }
   }
   lp->nloosevars++;

   lp->looseobjval = SCIPintervalGetInf(loose);

   return SCIP_OKAY;
}

/** informs LP, that given formerly column problem variable is now again a loose variable */
SCIP_RETCODE SCIPlpUpdateVarLoose(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< problem variable that changed from COLUMN to LOOSE */
   )
{
   assert(set != NULL);

   if( set->misc_exactsolve )
   {
      SCIP_CALL( lpUpdateVarLooseProved(lp, set, var) );
   }
   else
   {
      SCIP_CALL( lpUpdateVarLoose(lp, set, var) );
   }

   return SCIP_OKAY;
}

/** decrease the number of loose variables by one */
void SCIPlpDecNLoosevars(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->nloosevars > 0);

   lp->nloosevars--;

   /* get rid of numerical problems: set loose objective value explicitly to zero, if no loose variables remain */
   if( lp->nloosevars == 0 )
   {
      assert(lp->looseobjvalinf == 0);
      lp->looseobjval = 0.0;
   }
}

/** stores the LP solution in the columns and rows */
SCIP_RETCODE SCIPlpGetSol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            dualfeasible        /**< pointer to store whether the solution is dual feasible, or NULL */
   )
{
   SCIP_COL** lpicols;
   SCIP_ROW** lpirows;
   SCIP_Real* primsol;
   SCIP_Real* dualsol;
   SCIP_Real* activity;
   SCIP_Real* redcost;
   SCIP_Real primalbound;
   SCIP_Real dualbound;
   SCIP_Bool stillprimalfeasible;
   SCIP_Bool stilldualfeasible;
   int* cstat;
   int* rstat;
   SCIP_Longint lpcount;
   int nlpicols;
   int nlpirows;
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp->validsollp <= stat->lpcount);

   /* initialize return and feasibility flags; if primal oder dual feasibility shall not be checked, we set the
    * corresponding flag immediately to FALSE to skip all checks
    */
   if( primalfeasible == NULL )
      stillprimalfeasible = FALSE;
   else
   {
      *primalfeasible = TRUE;
      stillprimalfeasible = TRUE;
   }
   if( dualfeasible == NULL )
      stilldualfeasible = FALSE;
   else
   {
      *dualfeasible = TRUE;
      stilldualfeasible = TRUE;
   }

   /* check if the values are already calculated */
   if( lp->validsollp == stat->lpcount )
      return SCIP_OKAY;
   lp->validsollp = stat->lpcount;

   SCIPsetDebugMsg(set, "getting new LP solution %" SCIP_LONGINT_FORMAT " for solstat %d\n",
      stat->lpcount, SCIPlpGetSolstat(lp));

   lpicols = lp->lpicols;
   lpirows = lp->lpirows;
   nlpicols = lp->nlpicols;
   nlpirows = lp->nlpirows;
   lpcount = stat->lpcount;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsol, nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualsol, nlpirows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activity, nlpirows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redcost, nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, nlpirows) );

   SCIP_CALL( SCIPlpiGetSol(lp->lpi, NULL, primsol, dualsol, activity, redcost) );
   if( lp->solisbasic )
   {
      SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );
   }
   else
   {
      BMSclearMemoryArray(cstat, nlpicols);
      BMSclearMemoryArray(rstat, nlpirows);
   }

   primalbound = 0.0;
   dualbound = 0.0;

   /* copy primal solution and reduced costs into columns */
   for( c = 0; c < nlpicols; ++c )
   {
      assert( 0 <= cstat[c] && cstat[c] < 4 );
      lpicols[c]->primsol = primsol[c];
      lpicols[c]->minprimsol = MIN(lpicols[c]->minprimsol, primsol[c]);
      lpicols[c]->maxprimsol = MAX(lpicols[c]->maxprimsol, primsol[c]);
      lpicols[c]->redcost = redcost[c];
      lpicols[c]->basisstatus = (unsigned int) cstat[c];
      lpicols[c]->validredcostlp = lpcount;
      if( stillprimalfeasible )
      {
         stillprimalfeasible =
            (SCIPsetIsInfinity(set, -lpicols[c]->lb) || !SCIPsetIsFeasNegative(set, lpicols[c]->primsol - lpicols[c]->lb))
            && (SCIPsetIsInfinity(set, lpicols[c]->ub) || !SCIPsetIsFeasPositive(set, lpicols[c]->primsol - lpicols[c]->ub));
         primalbound += (lpicols[c]->primsol * lpicols[c]->obj);
      }
      if( lp->lastlpalgo == SCIP_LPALGO_BARRIER )
      {
         double compslack;

         /* complementary slackness in barrier solutions is measured as product of primal slack and dual multiplier;
          * we use a slack of at most 1, because otherwise we multiply by something like SCIPinfinty() for unbounded
          * variables, which would magnify even the tiniest violation in the dual multiplier
          */
         if( stilldualfeasible )
         {
            compslack = MIN((lpicols[c]->primsol - lpicols[c]->lb), 1.0) * lpicols[c]->redcost;
            stilldualfeasible = !SCIPsetIsDualfeasPositive(set, compslack);
         }
         if( stilldualfeasible )
         {
            compslack = MIN((lpicols[c]->ub - lpicols[c]->primsol), 1.0) * lpicols[c]->redcost;
            stilldualfeasible = !SCIPsetIsDualfeasNegative(set, compslack);
         }

         SCIPsetDebugMsg(set, " col <%s> [%.9g,%.9g]: primsol=%.9f, redcost=%.9f, pfeas=%u/%u(%u), dfeas=%d/%d(%u)\n",
            SCIPvarGetName(lpicols[c]->var), lpicols[c]->lb, lpicols[c]->ub, lpicols[c]->primsol, lpicols[c]->redcost,
            SCIPsetIsFeasGE(set, lpicols[c]->primsol, lpicols[c]->lb),
            SCIPsetIsFeasLE(set, lpicols[c]->primsol, lpicols[c]->ub),
            primalfeasible != NULL ? stillprimalfeasible : TRUE,
            !SCIPsetIsDualfeasPositive(set, MIN((lpicols[c]->primsol - lpicols[c]->lb), 1.0) * lpicols[c]->redcost),
            !SCIPsetIsDualfeasNegative(set, MIN((lpicols[c]->ub - lpicols[c]->primsol), 1.0) * lpicols[c]->redcost),
            dualfeasible != NULL ? stilldualfeasible : TRUE);
      }
      else
      {
         /* if dual feasibility check is disabled, set reduced costs of basic variables to 0 */
         if( dualfeasible == NULL && lpicols[c]->basisstatus == (unsigned int) SCIP_BASESTAT_BASIC )
         {
            lpicols[c]->redcost = 0.0;
         }

         /* complementary slackness means that if a variable is not at its lower or upper bound, its reduced costs
          * must be non-positive or non-negative, respectively; in particular, if a variable is strictly within its
          * bounds, its reduced cost must be zero
          */
         if( stilldualfeasible
            && (SCIPsetIsInfinity(set, -lpicols[c]->lb) || SCIPsetIsFeasGT(set, lpicols[c]->primsol, lpicols[c]->lb)) )
            stilldualfeasible = !SCIPsetIsDualfeasPositive(set, lpicols[c]->redcost);
         if( stilldualfeasible
            && (SCIPsetIsInfinity(set, lpicols[c]->ub) || SCIPsetIsFeasLT(set, lpicols[c]->primsol, lpicols[c]->ub)) )
            stilldualfeasible = !SCIPsetIsDualfeasNegative(set, lpicols[c]->redcost);

         SCIPsetDebugMsg(set, " col <%s> [%.9g,%.9g]: primsol=%.9f, redcost=%.9f, pfeas=%u/%u(%u), dfeas=%d/%d(%u)\n",
            SCIPvarGetName(lpicols[c]->var), lpicols[c]->lb, lpicols[c]->ub, lpicols[c]->primsol, lpicols[c]->redcost,
            SCIPsetIsFeasGE(set, lpicols[c]->primsol, lpicols[c]->lb),
            SCIPsetIsFeasLE(set, lpicols[c]->primsol, lpicols[c]->ub),
            primalfeasible != NULL ? stillprimalfeasible : TRUE,
            !SCIPsetIsFeasGT(set, lpicols[c]->primsol, lpicols[c]->lb) || !SCIPsetIsDualfeasPositive(set, lpicols[c]->redcost),
            !SCIPsetIsFeasLT(set, lpicols[c]->primsol, lpicols[c]->ub) || !SCIPsetIsDualfeasNegative(set, lpicols[c]->redcost),
            dualfeasible != NULL ? stilldualfeasible : TRUE);
      }

      /* we intentionally use an exact positive/negative check because ignoring small reduced cost values may lead to a
       * wrong bound value; if the corresponding bound is +/-infinity, we use zero reduced cost (if stilldualfeasible is
       * TRUE, we are in the case that the reduced cost is tiny with wrong sign)
       */
      if( stilldualfeasible )
      {
         if( lpicols[c]->redcost > 0.0 && !SCIPsetIsInfinity(set, -lpicols[c]->lb) )
            dualbound += (lpicols[c]->redcost * lpicols[c]->lb);
         else if( lpicols[c]->redcost < 0.0 && !SCIPsetIsInfinity(set, lpicols[c]->ub) )
            dualbound += (lpicols[c]->redcost * lpicols[c]->ub);
      }
   }

   /* copy dual solution and activities into rows */
   for( r = 0; r < nlpirows; ++r )
   {
      assert( 0 <= rstat[r] && rstat[r] < 4 );
      lpirows[r]->dualsol = dualsol[r];
      lpirows[r]->activity = activity[r] + lpirows[r]->constant;
      lpirows[r]->basisstatus = (unsigned int) rstat[r]; /*lint !e732*/
      lpirows[r]->validactivitylp = lpcount;
      if( stillprimalfeasible )
      {
         stillprimalfeasible =
            (SCIPsetIsInfinity(set,-lpirows[r]->lhs) ||SCIPsetIsFeasGE(set, lpirows[r]->activity, lpirows[r]->lhs))
            && (SCIPsetIsInfinity(set, lpirows[r]->rhs) || SCIPsetIsFeasLE(set, lpirows[r]->activity, lpirows[r]->rhs));
      }
      if( lp->lastlpalgo == SCIP_LPALGO_BARRIER )
      {
         double compslack;

         /* complementary slackness in barrier solutions is measured as product of primal slack and dual multiplier;
          * we use a slack of at most 1, because otherwise we multiply by something like SCIPinfinity() for unbounded
          * variables, which would magnify even the tiniest violation in the dual multiplier
          */
         if( stilldualfeasible )
         {
            compslack = MIN((lpirows[r]->activity - lpirows[r]->lhs), 1.0) * lpirows[r]->dualsol;
            stilldualfeasible = !SCIPsetIsDualfeasPositive(set, compslack);
         }
         if( stilldualfeasible )
         {
            compslack = MIN((lpirows[r]->rhs - lpirows[r]->activity), 1.0) * lpirows[r]->dualsol;
            stilldualfeasible = !SCIPsetIsDualfeasNegative(set, compslack);
         }

         SCIPsetDebugMsg(set, " row <%s> [%.9g,%.9g]: activity=%.9f, dualsol=%.9f, pfeas=%u/%u(%u), dfeas=%d/%d(%u)\n",
            lpirows[r]->name, lpirows[r]->lhs, lpirows[r]->rhs, lpirows[r]->activity, lpirows[r]->dualsol,
            SCIPsetIsFeasGE(set, lpirows[r]->activity, lpirows[r]->lhs),
            SCIPsetIsFeasLE(set, lpirows[r]->activity, lpirows[r]->rhs),
            primalfeasible != NULL ? stillprimalfeasible : TRUE,
            !SCIPsetIsDualfeasPositive(set, MIN((lpirows[r]->activity - lpirows[r]->lhs), 1.0) * lpirows[r]->dualsol),
            !SCIPsetIsDualfeasNegative(set, MIN((lpirows[r]->rhs - lpirows[r]->activity), 1.0) * lpirows[r]->dualsol),
            dualfeasible != NULL ? stilldualfeasible : TRUE);
      }
      else
      {
         /* complementary slackness means that if the activity of a row is not at its left-hand or right-hand side,
          * its dual multiplier must be non-positive or non-negative, respectively; in particular, if the activity is
          * strictly within left-hand and right-hand side, its dual multiplier must be zero
          */
         if( stilldualfeasible &&
               (SCIPsetIsInfinity(set, -lpirows[r]->lhs) || SCIPsetIsFeasGT(set, lpirows[r]->activity, lpirows[r]->lhs)) )
            stilldualfeasible = !SCIPsetIsDualfeasPositive(set, lpirows[r]->dualsol);
         if( stilldualfeasible &&
               (SCIPsetIsInfinity(set,lpirows[r]->rhs) || SCIPsetIsFeasLT(set, lpirows[r]->activity, lpirows[r]->rhs)) )
            stilldualfeasible = !SCIPsetIsDualfeasNegative(set, lpirows[r]->dualsol);

         SCIPsetDebugMsg(set, " row <%s> [%.9g,%.9g] + %.9g: activity=%.9f, dualsol=%.9f, pfeas=%u/%u(%u), dfeas=%d/%d(%u)\n",
            lpirows[r]->name, lpirows[r]->lhs, lpirows[r]->rhs, lpirows[r]->constant, lpirows[r]->activity, lpirows[r]->dualsol,
            SCIPsetIsFeasGE(set, lpirows[r]->activity, lpirows[r]->lhs),
            SCIPsetIsFeasLE(set, lpirows[r]->activity, lpirows[r]->rhs),
            primalfeasible != NULL ? stillprimalfeasible : TRUE,
            !SCIPsetIsFeasGT(set, lpirows[r]->activity, lpirows[r]->lhs) || !SCIPsetIsDualfeasPositive(set, lpirows[r]->dualsol),
            !SCIPsetIsFeasLT(set, lpirows[r]->activity, lpirows[r]->rhs) || !SCIPsetIsDualfeasNegative(set, lpirows[r]->dualsol),
            dualfeasible != NULL ? stilldualfeasible : TRUE);
      }

      /* we intentionally use an exact positive/negative check because ignoring small dual multipliers may lead to a
       * wrong bound value; if the corresponding side is +/-infinity, we use a zero dual multiplier (if
       * stilldualfeasible is TRUE, we are in the case that the dual multiplier is tiny with wrong sign)
       */
      if( stilldualfeasible )
      {
         if( lpirows[r]->dualsol > 0.0 && !SCIPsetIsInfinity(set, -lpirows[r]->lhs) )
            dualbound += (lpirows[r]->dualsol * (lpirows[r]->lhs - lpirows[r]->constant));
         else if( lpirows[r]->dualsol < 0.0 && !SCIPsetIsInfinity(set, lpirows[r]->rhs) )
            dualbound += (lpirows[r]->dualsol * (lpirows[r]->rhs - lpirows[r]->constant));
      }
   }

   /* if the objective value returned by the LP solver is smaller than the internally computed primal bound, then we
    * declare the solution primal infeasible; we assume primalbound and lp->lpobjval to be equal if they are both +/-
    * infinity
    */
   /**@todo alternatively, if otherwise the LP solution is feasible, we could simply update the objective value */
   if( stillprimalfeasible && !(SCIPsetIsInfinity(set, primalbound) && SCIPsetIsInfinity(set, lp->lpobjval))
      && !(SCIPsetIsInfinity(set, -primalbound) && SCIPsetIsInfinity(set, -lp->lpobjval)) )
   {
      stillprimalfeasible = SCIPsetIsFeasLE(set, primalbound, lp->lpobjval);
      SCIPsetDebugMsg(set, " primalbound=%.9f, lpbound=%.9g, pfeas=%u(%u)\n", primalbound, lp->lpobjval,
         SCIPsetIsFeasLE(set, primalbound, lp->lpobjval), primalfeasible != NULL ? stillprimalfeasible : TRUE);
   }

   /* if the objective value returned by the LP solver is smaller than the internally computed dual bound, we declare
    * the solution dual infeasible; we assume dualbound and lp->lpobjval to be equal if they are both +/- infinity
    */
   /**@todo alternatively, if otherwise the LP solution is feasible, we could simply update the objective value */
   if( stilldualfeasible && !(SCIPsetIsInfinity(set, dualbound) && SCIPsetIsInfinity(set, lp->lpobjval))
      && !(SCIPsetIsInfinity(set, -dualbound) && SCIPsetIsInfinity(set, -lp->lpobjval)) )
   {
      stilldualfeasible =  SCIPsetIsFeasGE(set, dualbound, lp->lpobjval);
      SCIPsetDebugMsg(set, " dualbound=%.9f, lpbound=%.9g, dfeas=%u(%u)\n", dualbound, lp->lpobjval,
         SCIPsetIsFeasGE(set, dualbound, lp->lpobjval), dualfeasible != NULL ? stilldualfeasible : TRUE);
   }

   if( primalfeasible != NULL )
      *primalfeasible = stillprimalfeasible;
   if( dualfeasible != NULL )
      *dualfeasible = stilldualfeasible;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);
   SCIPsetFreeBufferArray(set, &redcost);
   SCIPsetFreeBufferArray(set, &activity);
   SCIPsetFreeBufferArray(set, &dualsol);
   SCIPsetFreeBufferArray(set, &primsol);

   return SCIP_OKAY;
}

/** stores LP solution with infinite objective value in the columns and rows */
SCIP_RETCODE SCIPlpGetUnboundedSol(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            primalfeasible,     /**< pointer to store whether the solution is primal feasible, or NULL */
   SCIP_Bool*            rayfeasible         /**< pointer to store whether the primal ray is a feasible unboundedness proof, or NULL */
   )
{
   SCIP_COL** lpicols;
   SCIP_ROW** lpirows;
   SCIP_Real* primsol;
   SCIP_Real* activity;
   SCIP_Real* ray;
   SCIP_Real rayobjval;
   SCIP_Real rayscale;
   SCIP_Longint lpcount;
   int nlpicols;
   int nlpirows;
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
   assert(SCIPsetIsInfinity(set, -lp->lpobjval));
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp->validsollp <= stat->lpcount);

   if( primalfeasible != NULL )
      *primalfeasible = TRUE;
   if( rayfeasible != NULL )
      *rayfeasible = TRUE;

   /* check if the values are already calculated */
   if( lp->validsollp == stat->lpcount )
      return SCIP_OKAY;

   /* check if the LP solver is able to provide a primal unbounded ray */
   if( !SCIPlpiHasPrimalRay(lp->lpi) )
   {
      SCIPerrorMessage("LP solver has no primal ray to prove unboundedness\n");
      return SCIP_LPERROR;
   }

   lp->validsollp = stat->lpcount;

   SCIPsetDebugMsg(set, "getting new unbounded LP solution %" SCIP_LONGINT_FORMAT "\n", stat->lpcount);

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsol, lp->nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &activity, lp->nlpirows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ray, lp->nlpicols) );

   /* get primal feasible point */
   SCIP_CALL( SCIPlpiGetSol(lp->lpi, NULL, primsol, NULL, activity, NULL) );

   /* get primal unbounded ray */
   SCIP_CALL( SCIPlpiGetPrimalRay(lp->lpi, ray) );

   lpicols = lp->lpicols;
   lpirows = lp->lpirows;
   nlpicols = lp->nlpicols;
   nlpirows = lp->nlpirows;
   lpcount = stat->lpcount;

   /* calculate the objective value decrease of the ray */
   rayobjval = 0.0;
   for( c = 0; c < nlpicols; ++c )
   {
      assert(lpicols[c] != NULL);
      assert(lpicols[c]->var != NULL);

      /* there should only be a nonzero value in the ray if there is no finite bound in this direction */
      if( rayfeasible != NULL )
         *rayfeasible = *rayfeasible
            && (!SCIPsetIsNegative(set, ray[c]) || SCIPsetIsInfinity(set, -lpicols[c]->lb))
            && (!SCIPsetIsPositive(set, ray[c]) || SCIPsetIsInfinity(set,  lpicols[c]->ub));

      /* check primal feasibility of (finite) primal solution; note that the comparisons ensure that the primal
       * solution is within SCIP's infinity bounds; otherwise the rayscale below is not well-defined
       */
      if( primalfeasible != NULL )
         *primalfeasible = *primalfeasible
            && !SCIPsetIsFeasNegative(set, primsol[c] - lpicols[c]->lb)
            && !SCIPsetIsFeasPositive(set, primsol[c] - lpicols[c]->ub);

      if( !SCIPsetIsZero(set, ray[c]) )
         rayobjval += ray[c] * lpicols[c]->obj;
   }

   /* if the finite point is already infeasible, we do not have to add the ray */
   if( primalfeasible != NULL && !(*primalfeasible) )
   {
      rayscale = 0.0;
   }
   /* if the ray is already infeasible (due to numerics), we do not want to add the ray */
   else if( rayfeasible != NULL && !(*rayfeasible) )
   {
      rayscale = 0.0;
   }
   /* due to numerical problems, the objective of the ray might be nonnegative,
    *
    * @todo How to check for negative objective value here?
    */
   else if( !SCIPsetIsNegative(set, rayobjval) )
   {
      if( rayfeasible != NULL )
      {
         *rayfeasible = FALSE;
      }

      rayscale = 0.0;
   }
   else
   {
      assert(rayobjval != 0.0);

      /* scale the ray, such that the resulting point has infinite objective value */
      rayscale = -2*SCIPsetInfinity(set)/rayobjval;
      assert(SCIPsetIsFeasPositive(set, rayscale));

      /* ensure that unbounded point does not violate the bounds of the variables */
      for( c = 0; c < nlpicols; ++c )
      {
         if( SCIPsetIsPositive(set, ray[c]) )
            rayscale = MIN(rayscale, (lpicols[c]->ub - primsol[c])/ray[c]);
         else if( SCIPsetIsNegative(set, ray[c]) )
            rayscale = MIN(rayscale, (lpicols[c]->lb - primsol[c])/ray[c]);

         assert(SCIPsetIsFeasPositive(set, rayscale));
      }
   }

   SCIPsetDebugMsg(set, "unbounded LP solution: rayobjval=%f, rayscale=%f\n", rayobjval, rayscale);

   /* calculate the unbounded point: x' = x + rayscale * ray */
   for( c = 0; c < nlpicols; ++c )
   {
      if( SCIPsetIsZero(set, ray[c]) )
         lpicols[c]->primsol = primsol[c];
      else
      {
         SCIP_Real primsolval = primsol[c] + rayscale * ray[c];
         lpicols[c]->primsol = MAX(-SCIPsetInfinity(set), MIN(SCIPsetInfinity(set), primsolval)); /*lint !e666*/
      }
      lpicols[c]->redcost = SCIP_INVALID;
      lpicols[c]->validredcostlp = -1;
   }

   for( r = 0; r < nlpirows; ++r )
   {
      lpirows[r]->dualsol = SCIP_INVALID;
      lpirows[r]->activity = activity[r] + lpirows[r]->constant;
      lpirows[r]->validactivitylp = lpcount;

      /* check for feasibility of the rows */
      if( primalfeasible != NULL )
         *primalfeasible = *primalfeasible
            && (SCIPsetIsInfinity(set, -lpirows[r]->lhs) || SCIPsetIsFeasGE(set, lpirows[r]->activity, lpirows[r]->lhs))
            && (SCIPsetIsInfinity(set, lpirows[r]->rhs) || SCIPsetIsFeasLE(set, lpirows[r]->activity, lpirows[r]->rhs));
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &ray);
   SCIPsetFreeBufferArray(set, &activity);
   SCIPsetFreeBufferArray(set, &primsol);

   return SCIP_OKAY;
}

/** returns primal ray proving the unboundedness of the current LP */
SCIP_RETCODE SCIPlpGetPrimalRay(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            ray                 /**< array for storing primal ray values, they are stored w.r.t. the problem index of the variables,
                                              *   so the size of this array should be at least number of active variables
                                              *   (all entries have to be initialized to 0 before) */
   )
{
   SCIP_COL** lpicols;
   SCIP_Real* lpiray;
   SCIP_VAR* var;
   int nlpicols;
   int c;

   assert(lp != NULL);
   assert(set != NULL);
   assert(ray != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
   assert(SCIPsetIsInfinity(set, -lp->lpobjval));

   /* check if the LP solver is able to provide a primal unbounded ray */
   if( !SCIPlpiHasPrimalRay(lp->lpi) )
   {
      SCIPerrorMessage("LP solver has no primal ray for unbounded LP\n");
      return SCIP_LPERROR;
   }

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lpiray, lp->nlpicols) );

   SCIPsetDebugMsg(set, "getting primal ray values\n");

   /* get primal unbounded ray */
   SCIP_CALL( SCIPlpiGetPrimalRay(lp->lpi, lpiray) );

   lpicols = lp->lpicols;
   nlpicols = lp->nlpicols;

   /* store the ray values of active problem variables */
   for( c = 0; c < nlpicols; c++ )
   {
      assert(lpicols[c] != NULL);

      var = lpicols[c]->var;
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) != -1);
      ray[SCIPvarGetProbindex(var)] = lpiray[c];
   }

   SCIPsetFreeBufferArray(set, &lpiray);

   return SCIP_OKAY;
}

/** stores the dual Farkas multipliers for infeasibility proof in rows */
SCIP_RETCODE SCIPlpGetDualfarkas(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COL** lpicols;
   SCIP_ROW** lpirows;
   SCIP_Real* dualfarkas;
   int nlpicols;
   int nlpirows;
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
   assert(set != NULL);
   assert(stat != NULL);
   assert(lp->validfarkaslp <= stat->lpcount);

   /* check if the values are already calculated */
   if( lp->validfarkaslp == stat->lpcount )
      return SCIP_OKAY;
   lp->validfarkaslp = stat->lpcount;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualfarkas, lp->nlpirows) );

   /* get dual Farkas infeasibility proof */
   SCIP_CALL( SCIPlpiGetDualfarkas(lp->lpi, dualfarkas) );

   lpicols = lp->lpicols;
   lpirows = lp->lpirows;
   nlpicols = lp->nlpicols;
   nlpirows = lp->nlpirows;

   /* store infeasibility proof in rows */
   SCIPsetDebugMsg(set, "LP is infeasible:\n");
   for( r = 0; r < nlpirows; ++r )
   {
      SCIPsetDebugMsg(set, " row <%s>: dualfarkas=%f\n", lpirows[r]->name, dualfarkas[r]);
      lpirows[r]->dualfarkas = dualfarkas[r];
      lpirows[r]->dualsol = SCIP_INVALID;
      lpirows[r]->activity = 0.0;
      lpirows[r]->validactivitylp = -1L;
      lpirows[r]->basisstatus = (unsigned int) SCIP_BASESTAT_BASIC;
   }

   /* set columns as invalid */
   for( c = 0; c < nlpicols; ++c )
   {
      lpicols[c]->primsol = SCIP_INVALID;
      lpicols[c]->redcost = SCIP_INVALID;
      lpicols[c]->validredcostlp = -1L;
      lpicols[c]->validfarkaslp = -1L;
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &dualfarkas);

   return SCIP_OKAY;
}

/** get number of iterations used in last LP solve */
SCIP_RETCODE SCIPlpGetIterations(
   SCIP_LP*              lp,                 /**< current LP data */
   int*                  iterations          /**< pointer to store the iteration count */
   )
{
   assert(lp != NULL);

   SCIP_CALL( SCIPlpiGetIterations(lp->lpi, iterations) );

   return SCIP_OKAY;
}

/** increases age of columns with solution value 0.0 and basic rows with activity not at its bounds,
 *  resets age of non-zero columns and sharp rows
 */
SCIP_RETCODE SCIPlpUpdateAges(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_COL** lpicols;
   SCIP_ROW** lpirows;
   int nlpicols;
   int nlpirows;
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(lp->nlpicols == lp->ncols);
   assert(lp->nlpirows == lp->nrows);
   assert(stat != NULL);
   assert(lp->validsollp == stat->lpcount);

   SCIPdebugMessage("updating LP ages\n");

   lpicols = lp->lpicols;
   lpirows = lp->lpirows;
   nlpicols = lp->nlpicols;
   nlpirows = lp->nlpirows;

   for( c = 0; c < nlpicols; ++c )
   {
      assert(lpicols[c] == lp->cols[c]);
      if( lpicols[c]->primsol == 0.0 )  /* non-basic columns to remove are exactly at 0.0 */
         lpicols[c]->age++;
      else
         lpicols[c]->age = 0;
      /*SCIPstatDebugMsg(stat, " -> col <%s>: primsol=%f, age=%d\n",
        SCIPvarGetName(lpicols[c]->var), lpicols[c]->primsol, lpicols[c]->age);*/
   }

   for( r = 0; r < nlpirows; ++r )
   {
      lpirows[r]->nlpsaftercreation++;
      assert(lpirows[r] == lp->rows[r]);

      if( lpirows[r]->dualsol == 0.0 ) /* basic rows to remove are exactly at 0.0 */
      {
         lpirows[r]->age++;
      }
      else
      {
         lpirows[r]->activeinlpcounter++;
         lpirows[r]->age = 0;
      }
      /*debugMsg(scip, " -> row <%s>: activity=%f, age=%d\n", lpirows[r]->name, lpirows[r]->activity, lpirows[r]->age);*/
   }

   return SCIP_OKAY;
}

/* deletes the marked columns from the LP and the LP interface */
static
SCIP_RETCODE lpDelColset(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int*                  coldstat            /**< deletion status of columns:  1 if column should be deleted, 0 if not */
   )
{
   SCIP_COL* col;
   int ncols;
   int c;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(!lp->diving);
   assert(coldstat != NULL);
   assert(lp->nlazycols <= lp->ncols); 

   ncols = lp->ncols;

   /* delete columns in LP solver */
   SCIP_CALL( SCIPlpiDelColset(lp->lpi, coldstat) );

   /* update LP data respectively */
   for( c = 0; c < ncols; ++c )
   {
      col = lp->cols[c];
      assert(col != NULL);
      assert(col == lp->lpicols[c]);
      assert(coldstat[c] <= c);
      col->lppos = coldstat[c];
      if( coldstat[c] == -1 )
      {
         assert(col->removable);

         /* mark column to be deleted from the LPI, update column arrays of all linked rows, and update the objective
          * function vector norms
          */
         markColDeleted(col);
         colUpdateDelLP(col, set);
         lpUpdateObjNorms(lp, set, col->unchangedobj, 0.0);
         col->lpdepth = -1;

         lp->cols[c] = NULL;
         lp->lpicols[c] = NULL;
         lp->ncols--;
         lp->nremovablecols--;
         lp->nlpicols--;
      }
      else if( coldstat[c] < c )
      {
         assert(lp->cols[coldstat[c]] == NULL);
         assert(lp->lpicols[coldstat[c]] == NULL);
         lp->cols[coldstat[c]] = col;
         lp->lpicols[coldstat[c]] = col;
         lp->cols[coldstat[c]]->lppos = coldstat[c];
         lp->cols[coldstat[c]]->lpipos = coldstat[c];
         lp->cols[c] = NULL;
         lp->lpicols[c] = NULL;
      }
   }

   /* remove columns which are deleted from the lazy column array */
   c = 0;
   while( c < lp->nlazycols )
   {
      if( lp->lazycols[c]->lpipos < 0 )
      {
         lp->lazycols[c] = lp->lazycols[lp->nlazycols-1];
         lp->nlazycols--;
      }
      else
         c++;
   }

   /* mark LP to be unsolved */
   if( lp->ncols < ncols )
   {
      assert(lp->ncols == lp->nlpicols);
      assert(lp->nchgcols == 0);
      assert(lp->flushed);

      lp->lpifirstchgcol = lp->nlpicols;

      /* mark the current solution invalid */
      lp->solved = FALSE;
      lp->primalfeasible = FALSE;
      lp->primalchecked = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   checkLazyColArray(lp, set);
   checkLinks(lp);

   return SCIP_OKAY;
}

/* deletes the marked rows from the LP and the LP interface */
static
SCIP_RETCODE lpDelRowset(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   int*                  rowdstat            /**< deletion status of rows:  1 if row should be deleted, 0 if not */
   )
{
   SCIP_ROW* row;
   int nrows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->nrows == lp->nlpirows);
   assert(!lp->diving);
   assert(rowdstat != NULL);

   nrows = lp->nrows;

   /* delete rows in LP solver */
   SCIP_CALL( SCIPlpiDelRowset(lp->lpi, rowdstat) );

   /* update LP data respectively */
   for( r = 0; r < nrows; ++r )
   {
      row = lp->rows[r];
      assert(row == lp->lpirows[r]);
      assert(rowdstat[r] <= r);
      assert(row != NULL);
      row->lppos = rowdstat[r];
      if( rowdstat[r] == -1 )
      {
         if( row->removable )
            lp->nremovablerows--;

         /* mark row to be deleted from the LPI and update row arrays of all linked columns */
         markRowDeleted(row);
         rowUpdateDelLP(row);
         row->lpdepth = -1;

         /* check, if row deletion events are tracked
          * if so, issue ROWDELETEDLP event
          */
         if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDLP) != 0 )
         {
            SCIP_EVENT* event;

            SCIP_CALL( SCIPeventCreateRowDeletedLP(&event, blkmem, row) );
            SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
         }

         SCIP_CALL( SCIProwRelease(&lp->lpirows[r], blkmem, set, lp) );
         SCIProwUnlock(lp->rows[r]);
         SCIP_CALL( SCIProwRelease(&lp->rows[r], blkmem, set, lp) );
         assert(lp->lpirows[r] == NULL);
         assert(lp->rows[r] == NULL);
         lp->nrows--;
         lp->nlpirows--;
      }
      else if( rowdstat[r] < r )
      {
         assert(lp->rows[rowdstat[r]] == NULL);
         assert(lp->lpirows[rowdstat[r]] == NULL);
         lp->rows[rowdstat[r]] = row;
         lp->lpirows[rowdstat[r]] = row;
         lp->rows[rowdstat[r]]->lppos = rowdstat[r];
         lp->rows[rowdstat[r]]->lpipos = rowdstat[r];
         lp->rows[r] = NULL;
         lp->lpirows[r] = NULL;
      }
   }

   /* mark LP to be unsolved */
   if( lp->nrows < nrows )
   {
      assert(lp->nrows == lp->nlpirows);
      assert(lp->nchgrows == 0);
      assert(lp->flushed);

      lp->lpifirstchgrow = lp->nlpirows;

      /* mark the current solution invalid */
      lp->solved = FALSE;
      lp->dualfeasible = FALSE;
      lp->dualchecked = FALSE;
      lp->lpobjval = SCIP_INVALID;
      lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   checkLinks(lp);

   return SCIP_OKAY;
}

/** removes all non-basic columns, that are too old, beginning with the given firstcol */
static
SCIP_RETCODE lpRemoveObsoleteCols(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   int                   firstcol            /**< first column to check for clean up */
   )
{
   SCIP_COL** cols;
#ifndef NDEBUG
   SCIP_COL** lpicols;
#endif
   int* coldstat;
   int ncols;
   int ndelcols;
   int c;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(lp->nremovablecols <= lp->ncols);
   assert(!lp->diving);
   assert(set != NULL);
   assert(stat != NULL);

   if( lp->nremovablecols == 0 || set->lp_colagelimit == -1 || !lp->solisbasic )
      return SCIP_OKAY;

   ncols = lp->ncols;
   cols = lp->cols;
#ifndef NDEBUG
   lpicols = lp->lpicols;
#endif

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &coldstat, ncols) );

   /* mark obsolete columns to be deleted */
   ndelcols = 0;
   BMSclearMemoryArray(coldstat, ncols);
   for( c = firstcol; c < ncols; ++c )
   {
      assert(cols[c] == lpicols[c]);
      assert(cols[c]->lppos == c);
      assert(cols[c]->lpipos == c);
      if( cols[c]->removable
         && cols[c]->obsoletenode != stat->nnodes /* don't remove column a second time from same node (avoid cycling), or a first time if marked nonremovable locally */
         && cols[c]->age > set->lp_colagelimit
         && (SCIP_BASESTAT)cols[c]->basisstatus != SCIP_BASESTAT_BASIC
         && SCIPsetIsZero(set, SCIPcolGetBestBound(cols[c])) ) /* bestbd != 0 -> column would be priced in next time */
      {
         assert(cols[c]->primsol == 0.0);
         coldstat[c] = 1;
         ndelcols++;
         cols[c]->obsoletenode = stat->nnodes;
         SCIPsetDebugMsg(set, "removing obsolete col <%s>: primsol=%f, bounds=[%g,%g]\n",
            SCIPvarGetName(cols[c]->var), cols[c]->primsol, cols[c]->lb, cols[c]->ub);
      }
   }

   SCIPsetDebugMsg(set, "removing %d/%d obsolete columns from LP\n", ndelcols, ncols);

   /* delete the marked columns in the LP solver interface, update the LP respectively */
   if( ndelcols > 0 )
   {
      SCIP_CALL( lpDelColset(lp, set, coldstat) );
   }
   assert(lp->ncols == ncols - ndelcols);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &coldstat);

   return SCIP_OKAY;
}

/** removes all basic rows, that are too old, beginning with the given firstrow */
static
SCIP_RETCODE lpRemoveObsoleteRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   int                   firstrow            /**< first row to check for clean up */
   )
{
   SCIP_ROW** rows;
#ifndef NDEBUG
   SCIP_ROW** lpirows;
#endif
   int* rowdstat;
   int nrows;
   int ndelrows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->nrows == lp->nlpirows);
   assert(lp->nremovablerows <= lp->nrows);
   assert(!lp->diving);
   assert(set != NULL);
   assert(stat != NULL);

   if( lp->nremovablerows == 0 || set->lp_rowagelimit == -1 || !lp->solisbasic )
      return SCIP_OKAY;

   nrows = lp->nrows;
   rows = lp->rows;
#ifndef NDEBUG
   lpirows = lp->lpirows;
#endif

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rowdstat, nrows) );

   /* mark obsolete rows to be deleted */
   ndelrows = 0;
   BMSclearMemoryArray(rowdstat, nrows);
   for( r = firstrow; r < nrows; ++r )
   {
      assert(rows[r] == lpirows[r]);
      assert(rows[r]->lppos == r);
      assert(rows[r]->lpipos == r);
      if( rows[r]->removable
         && rows[r]->obsoletenode != stat->nnodes  /* don't remove row a second time from same node (avoid cycling), or a first time if marked nonremovable locally */
         && rows[r]->age > set->lp_rowagelimit
         && (SCIP_BASESTAT)rows[r]->basisstatus == SCIP_BASESTAT_BASIC )
      {
         rowdstat[r] = 1;
         ndelrows++;
         rows[r]->obsoletenode = stat->nnodes;
         SCIPsetDebugMsg(set, "removing obsolete row <%s>: activity=%f, sides=[%g,%g]\n",
            rows[r]->name, rows[r]->activity, rows[r]->lhs, rows[r]->rhs);
      }
   }

   SCIPsetDebugMsg(set, "removing %d/%d obsolete rows from LP\n", ndelrows, nrows);

   /* delete the marked rows in the LP solver interface, update the LP respectively */
   if( ndelrows > 0 )
   {
      SCIP_CALL( lpDelRowset(lp, blkmem, set, eventqueue, eventfilter, rowdstat) );
   }
   assert(lp->nrows == nrows - ndelrows);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &rowdstat);

   return SCIP_OKAY;
}

/** removes all non-basic columns and basic rows in the part of the LP created at the current node, that are too old */
SCIP_RETCODE SCIPlpRemoveNewObsoletes(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(set != NULL);

   SCIPsetDebugMsg(set, "removing obsolete columns starting with %d/%d, obsolete rows starting with %d/%d\n",
      lp->firstnewcol, lp->ncols, lp->firstnewrow, lp->nrows);

   if( lp->firstnewcol < lp->ncols )
   {
      SCIP_CALL( lpRemoveObsoleteCols(lp, set, stat, lp->firstnewcol) );
   }
   if( lp->firstnewrow < lp->nrows )
   {
      SCIP_CALL( lpRemoveObsoleteRows(lp, blkmem, set, stat, eventqueue, eventfilter, lp->firstnewrow) );
   }

   return SCIP_OKAY;
}

/** removes all non-basic columns and basic rows in whole LP, that are too old */
SCIP_RETCODE SCIPlpRemoveAllObsoletes(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(set != NULL);

   SCIPsetDebugMsg(set, "removing all obsolete columns and rows\n");

   if( 0 < lp->ncols )
   {
      SCIP_CALL( lpRemoveObsoleteCols(lp, set, stat, 0) );
   }
   if( 0 < lp->nrows )
   {
      SCIP_CALL( lpRemoveObsoleteRows(lp, blkmem, set, stat, eventqueue, eventfilter, 0) );
   }

   return SCIP_OKAY;
}

/** removes all non-basic columns at 0.0 beginning with the given firstcol */
static
SCIP_RETCODE lpCleanupCols(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   int                   firstcol            /**< first column to check for clean up */
   )
{
   SCIP_COL** cols;
   SCIP_COL** lpicols;
   int* coldstat;
   int ncols;
   int ndelcols;
   int c;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(!lp->diving);
   assert(stat != NULL);
   assert(lp->validsollp == stat->lpcount);
   assert(0 <= firstcol && firstcol < lp->ncols);

   if( lp->nremovablecols == 0 || !lp->solisbasic )
      return SCIP_OKAY;

   ncols = lp->ncols;
   cols = lp->cols;
   lpicols = lp->lpicols;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &coldstat, ncols) );

   /* mark unused columns to be deleted */
   ndelcols = 0;
   BMSclearMemoryArray(coldstat, ncols);
   for( c = firstcol; c < ncols; ++c )
   {
      assert(cols[c] == lpicols[c]);
      assert(cols[c]->lppos == c);
      assert(cols[c]->lpipos == c);
      if( lpicols[c]->removable
         && (SCIP_BASESTAT)lpicols[c]->basisstatus != SCIP_BASESTAT_BASIC
         && lpicols[c]->primsol == 0.0 /* non-basic columns to remove are exactly at 0.0 */
         && SCIPsetIsZero(set, SCIPcolGetBestBound(cols[c])) ) /* bestbd != 0 -> column would be priced in next time */
      {
         coldstat[c] = 1;
         ndelcols++;
      }
   }

   SCIPsetDebugMsg(set, "removing %d/%d unused columns from LP\n", ndelcols, ncols);

   /* delete the marked columns in the LP solver interface, update the LP respectively */
   if( ndelcols > 0 )
   {
      SCIP_CALL( lpDelColset(lp, set, coldstat) );
   }
   assert(lp->ncols == ncols - ndelcols);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &coldstat);

   return SCIP_OKAY;
}

/** removes all basic rows beginning with the given firstrow */
static
SCIP_RETCODE lpCleanupRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   int                   firstrow            /**< first row to check for clean up */
   )
{
#ifndef NDEBUG
   SCIP_ROW** rows;
#endif
   SCIP_ROW** lpirows;
   int* rowdstat;
   int nrows;
   int ndelrows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(lp->nrows == lp->nlpirows);
   assert(!lp->diving);
   assert(stat != NULL);
   assert(lp->validsollp == stat->lpcount);
   assert(0 <= firstrow && firstrow < lp->nrows);

   if( lp->nremovablerows == 0 || !lp->solisbasic  )
      return SCIP_OKAY;

#ifndef NDEBUG
   rows = lp->rows;
#endif
   nrows = lp->nrows;
   lpirows = lp->lpirows;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rowdstat, nrows) );

   /* mark unused rows to be deleted */
   ndelrows = 0;
   BMSclearMemoryArray(rowdstat, nrows);
   for( r = firstrow; r < nrows; ++r )
   {
      assert(rows[r] == lpirows[r]);
      assert(rows[r]->lppos == r);
      assert(rows[r]->lpipos == r);
      if( lpirows[r]->removable && (SCIP_BASESTAT)lpirows[r]->basisstatus == SCIP_BASESTAT_BASIC )
      {
         rowdstat[r] = 1;
         ndelrows++;
      }
   }

   SCIPsetDebugMsg(set, "removing %d/%d unused rows from LP\n", ndelrows, nrows);

   /* delete the marked rows in the LP solver interface, update the LP respectively */
   if( ndelrows > 0 )
   {
      SCIP_CALL( lpDelRowset(lp, blkmem, set, eventqueue, eventfilter, rowdstat) );
   }
   assert(lp->nrows == nrows - ndelrows);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &rowdstat);

   return SCIP_OKAY;
}

/** removes all non-basic columns at 0.0 and basic rows in the part of the LP created at the current node */
SCIP_RETCODE SCIPlpCleanupNew(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   SCIP_Bool cleanupcols;
   SCIP_Bool cleanuprows;

   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(set != NULL);

   /* check, if we want to clean up the columns and rows */
   cleanupcols = (root ? set->lp_cleanupcolsroot : set->lp_cleanupcols);
   cleanuprows = (root ? set->lp_cleanuprowsroot : set->lp_cleanuprows);

   SCIPsetDebugMsg(set, "removing unused columns starting with %d/%d (%u), unused rows starting with %d/%d (%u), LP algo: %d, basic sol: %u\n",
      lp->firstnewcol, lp->ncols, cleanupcols, lp->firstnewrow, lp->nrows, cleanuprows, lp->lastlpalgo, lp->solisbasic);

   if( cleanupcols && lp->firstnewcol < lp->ncols )
   {
      SCIP_CALL( lpCleanupCols(lp, set, stat, lp->firstnewcol) );
   }
   if( cleanuprows && lp->firstnewrow < lp->nrows )
   {
      SCIP_CALL( lpCleanupRows(lp, blkmem, set, stat, eventqueue, eventfilter, lp->firstnewrow) );
   }

   return SCIP_OKAY;
}

/** removes all non-basic columns at 0.0 and basic rows in the whole LP */
SCIP_RETCODE SCIPlpCleanupAll(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   SCIP_Bool cleanupcols;
   SCIP_Bool cleanuprows;

   assert(lp != NULL);
   assert(lp->solved);
   assert(!lp->diving);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
   assert(set != NULL);

   /* check, if we want to clean up the columns and rows */
   cleanupcols = (root ? set->lp_cleanupcolsroot : set->lp_cleanupcols);
   cleanuprows = (root ? set->lp_cleanuprowsroot : set->lp_cleanuprows);

   SCIPsetDebugMsg(set, "removing all unused columns (%u) and rows (%u), LP algo: %d, basic sol: %u\n",
      cleanupcols, cleanuprows, lp->lastlpalgo, lp->solisbasic);

   if( cleanupcols && 0 < lp->ncols )
   {
      SCIP_CALL( lpCleanupCols(lp, set, stat, 0) );
   }
   if( cleanuprows && 0 < lp->nrows )
   {
      SCIP_CALL( lpCleanupRows(lp, blkmem, set, stat, eventqueue, eventfilter, 0) );
   }

   return SCIP_OKAY;
}

/** removes all redundant rows that were added at the current node */
SCIP_RETCODE SCIPlpRemoveRedundantRows(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   )
{
#ifndef NDEBUG
   SCIP_ROW** rows;
#endif
   SCIP_ROW** lpirows;
   int* rowdstat;
   int nrows;
   int ndelrows;
   int r;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->ncols == lp->nlpicols);
   assert(lp->nrows == lp->nlpirows);
   assert(!lp->diving);
   assert(stat != NULL);
   assert(lp->validsollp == stat->lpcount);
   assert(lp->firstnewrow <= lp->nrows);

   if( lp->firstnewrow == lp->nrows )
      return SCIP_OKAY;

#ifndef NDEBUG
   rows = lp->rows;
#endif
   nrows = lp->nrows;
   lpirows = lp->lpirows;

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rowdstat, nrows) );

   /* mark redundant rows to be deleted (only delete basic rows!) */
   ndelrows = 0;
   BMSclearMemoryArray(rowdstat, nrows);
   for( r = lp->firstnewrow; r < nrows; ++r )
   {
      assert(rows[r] == lpirows[r]);
      assert(rows[r]->lppos == r);
      assert(rows[r]->lpipos == r);
      if( (!lp->solisbasic || (SCIP_BASESTAT)lpirows[r]->basisstatus == SCIP_BASESTAT_BASIC)
         && SCIProwIsRedundant(lpirows[r], set, stat) )
      {
         SCIPsetDebugMsg(set, "basic row <%s> is redundant: sides=[%g,%g], act=[%g,%g]\n",
            SCIProwGetName(lpirows[r]), SCIProwGetLhs(lpirows[r]), SCIProwGetRhs(lpirows[r]),
            SCIProwGetMinActivity(lpirows[r], set, stat), SCIProwGetMaxActivity(lpirows[r], set, stat));
         rowdstat[r] = 1;
         ndelrows++;
      }
   }

   SCIPsetDebugMsg(set, "removing %d/%d redundant basic rows from LP\n", ndelrows, nrows);

   /* delete the marked rows in the LP solver interface, update the LP respectively */
   if( ndelrows > 0 )
   {
      SCIP_CALL( lpDelRowset(lp, blkmem, set, eventqueue, eventfilter, rowdstat) );
   }
   assert(lp->nrows == nrows - ndelrows);

   /* release temporary memory */
   SCIPsetFreeBufferArray(set, &rowdstat);

   return SCIP_OKAY;
}

/** initiates LP diving */
SCIP_RETCODE SCIPlpStartDive(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   int c;
   int r;

   assert(lp != NULL);
   assert(lp->flushed || !lp->solved);
   assert(!lp->diving);
   assert(!lp->probing);
   assert(lp->divelpistate == NULL);
   assert(lp->divelpwasprimfeas);
   assert(lp->divelpwasdualfeas);
   assert(lp->validsollp <= stat->lpcount);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(lp->ndivechgsides == 0);

   SCIPsetDebugMsg(set, "diving started (LP flushed: %u, LP solved: %u, solstat: %d)\n",
      lp->flushed, lp->solved, SCIPlpGetSolstat(lp));

#ifndef NDEBUG
   for( c = 0; c < lp->ncols; ++c )
   {
      assert(lp->cols[c] != NULL);
      assert(lp->cols[c]->var != NULL);
      assert(SCIPvarGetStatus(lp->cols[c]->var) == SCIP_VARSTATUS_COLUMN);
      assert(SCIPvarGetCol(lp->cols[c]->var) == lp->cols[c]);
      assert(SCIPsetIsFeasEQ(set, SCIPvarGetObj(lp->cols[c]->var), lp->cols[c]->obj));
      assert(SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(lp->cols[c]->var), lp->cols[c]->lb));
      assert(SCIPsetIsFeasEQ(set, SCIPvarGetUbLocal(lp->cols[c]->var), lp->cols[c]->ub));
   }
#endif

   /* save current LPI state (basis information) */
   SCIP_CALL( SCIPlpiGetState(lp->lpi, blkmem, &lp->divelpistate) );
   lp->divelpwasprimfeas = lp->primalfeasible;
   lp->divelpwasdualfeas = lp->dualfeasible;
   lp->divelpwasprimchecked = lp->primalchecked;
   lp->divelpwasdualchecked = lp->dualchecked;

   /* save current LP values dependent on the solution */
   SCIP_CALL( lpStoreSolVals(lp, stat, blkmem) );
   assert(lp->storedsolvals != NULL);
   if( !set->lp_resolverestore && lp->solved )
   {
      SCIP_Bool store = TRUE;

      switch ( lp->lpsolstat )
      {
      case SCIP_LPSOLSTAT_OPTIMAL:
         SCIP_CALL( SCIPlpGetSol(lp, set, stat, NULL, NULL) );
         assert(lp->validsollp == stat->lpcount);
         break;
      case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
         SCIP_CALL( SCIPlpGetUnboundedSol(lp, set, stat, NULL, NULL) );
         assert(lp->validsollp == stat->lpcount);
         break;
      case SCIP_LPSOLSTAT_OBJLIMIT:
      case SCIP_LPSOLSTAT_ITERLIMIT:
      case SCIP_LPSOLSTAT_TIMELIMIT:
         SCIP_CALL( SCIPlpGetSol(lp, set, stat, NULL, NULL) );
         assert(lp->validsollp == stat->lpcount);
         break;
      case SCIP_LPSOLSTAT_INFEASIBLE:
         SCIP_CALL( SCIPlpGetDualfarkas(lp, set, stat) );
         break;
      case SCIP_LPSOLSTAT_NOTSOLVED:
      case SCIP_LPSOLSTAT_ERROR:
      default:
         store = FALSE;
      }

      if ( store )
      {
         for( c = 0; c < lp->ncols; ++c )
         {
            SCIP_CALL( colStoreSolVals(lp->cols[c], blkmem) );
         }
         for( r = 0; r < lp->nrows; ++r )
         {
            SCIP_CALL( rowStoreSolVals(lp->rows[r], blkmem, lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE) );
         }
      }
   }

   /* store LPI iteration limit */
   SCIP_CALL( SCIPlpiGetIntpar(lp->lpi, SCIP_LPPAR_LPITLIM, &lp->divinglpiitlim) );

   /* remember the number of domain changes */
   lp->divenolddomchgs = stat->domchgcount;

   /* store current number of rows */
   lp->ndivingrows = lp->nrows;

   /* switch to diving mode */
   lp->diving = TRUE;

   return SCIP_OKAY;
}

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
SCIP_RETCODE SCIPlpEndDive(
   SCIP_LP*              lp,                 /**< current LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_VAR**            vars,               /**< array with all active variables */
   int                   nvars               /**< number of active variables */
   )
{
   SCIP_VAR* var;
   int v;

   assert(lp != NULL);
   assert(lp->diving);
   assert(blkmem != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIPsetDebugMsg(set, "diving ended (LP flushed: %u, solstat: %d)\n", lp->flushed, SCIPlpGetSolstat(lp));

   /* reset all columns' objective values and bounds to its original values */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPcolChgObj(SCIPvarGetCol(var), set, lp, SCIPvarGetObj(var)) );
         SCIP_CALL( SCIPcolChgLb(SCIPvarGetCol(var), set, lp, SCIPvarGetLbLocal(var)) );
         SCIP_CALL( SCIPcolChgUb(SCIPvarGetCol(var), set, lp, SCIPvarGetUbLocal(var)) );
      }
   }

   /* remove rows which were added in diving mode */
   SCIP_CALL( SCIPlpShrinkRows(lp, blkmem, set, eventqueue, eventfilter, lp->ndivingrows) );

   /* undo changes to left hand sides and right hand sides */
   while( lp->ndivechgsides > 0 )
   {
      SCIP_Real oldside;
      SCIP_SIDETYPE sidetype;
      SCIP_ROW* row;

      lp->ndivechgsides--;
      oldside = lp->divechgsides[lp->ndivechgsides];
      sidetype = lp->divechgsidetypes[lp->ndivechgsides];
      row = lp->divechgrows[lp->ndivechgsides];

      if( sidetype == SCIP_SIDETYPE_LEFT )
      {
         SCIP_CALL( SCIProwChgLhs(row, blkmem, set, eventqueue, lp, oldside) );
      }
      else
      {
         SCIP_CALL( SCIProwChgRhs(row, blkmem, set, eventqueue, lp, oldside) );
      }
   }

   /* restore LPI iteration limit */
   SCIP_CALL( lpSetIterationLimit(lp, lp->divinglpiitlim) );

   /* reload LPI state saved at start of diving and free it afterwards; it may be NULL, in which case simply nothing
    * happens
    */
   SCIP_CALL( SCIPlpSetState(lp, blkmem, set, eventqueue, lp->divelpistate,
         lp->divelpwasprimfeas, lp->divelpwasprimchecked, lp->divelpwasdualfeas, lp->divelpwasdualchecked) );
   SCIP_CALL( SCIPlpFreeState(lp, blkmem, &lp->divelpistate) );
   lp->divelpwasprimfeas = TRUE;
   lp->divelpwasdualfeas = TRUE;
   lp->divelpwasprimchecked = TRUE;
   lp->divelpwasdualchecked = TRUE;
   assert(lp->divelpistate == NULL);

   /* switch to standard (non-diving) mode */
   lp->diving = FALSE;
   lp->divingobjchg = FALSE;

   /* if the LP was solved before starting the dive, but not to optimality (or unboundedness), then we need to solve the
    * LP again to reset the solution (e.g. we do not save the Farkas proof for infeasible LPs, because we assume that we
    * are not called in this case, anyway); restoring by solving the LP again in either case can be forced by setting
    * the parameter resolverestore to TRUE
    * restoring an unbounded ray after solve does not seem to work currently (bug 631), so we resolve also in this case
    */
   assert(lp->storedsolvals != NULL);
   if( lp->storedsolvals->lpissolved
      && (set->lp_resolverestore || lp->storedsolvals->lpsolstat != SCIP_LPSOLSTAT_OPTIMAL || lp->divenolddomchgs < stat->domchgcount) )
   {
      SCIP_Bool lperror;

      SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr,  blkmem, stat, eventqueue, eventfilter, prob, -1LL, FALSE, FALSE, FALSE, &lperror) );
      if( lperror )
      {
         lpNumericalTroubleMessage(messagehdlr, set, stat, SCIP_VERBLEVEL_FULL, "unresolved when resolving LP after diving");
         lp->resolvelperror = TRUE;
      }
      else if( SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL
         && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_INFEASIBLE
         && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY
         && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OBJLIMIT )
      {
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "LP was not resolved to a sufficient status after diving\n");
         lp->resolvelperror = TRUE;
      }
   }
   /* otherwise, we can just reload the buffered LP solution values at start of diving; this has the advantage that we
    * are guaranteed to continue with the same LP status as before diving, while in numerically difficult cases, a
    * re-solve as above can lead to a different LP status
    */
   else
   {
      int c;
      int r;

      /* if there are lazy bounds, remove them from the LP */
      if( lp->nlazycols > 0 )
      {
         /* @todo avoid loosing primal feasibility here after changing the objective already did destroy dual feasibility;
          * first resolve LP?
          */
         SCIP_CALL( updateLazyBounds(lp, set) );
         assert(lp->diving == lp->divinglazyapplied);

         /* flush changes to the LP solver */
         SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );
      }

      /* increment lp counter to ensure that we do not use solution values from the last solved diving lp */
      SCIPstatIncrement(stat, set, lpcount);

      /* restore LP solution values in lp data, columns and rows */
      if( lp->storedsolvals->lpissolved &&
         (lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL ||
            lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDEDRAY ||
            lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT ||
            lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_ITERLIMIT ||
            lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_TIMELIMIT ||
            lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE)
         )
      {
         SCIP_CALL( lpRestoreSolVals(lp, blkmem, stat->lpcount) );

         for( c = 0; c < lp->ncols; ++c )
         {
            SCIP_CALL( colRestoreSolVals(lp->cols[c], blkmem, stat->lpcount, set->lp_freesolvalbuffers) );
         }
         for( r = 0; r < lp->nrows; ++r )
         {
            SCIP_CALL( rowRestoreSolVals(lp->rows[r], blkmem, stat->lpcount, set->lp_freesolvalbuffers, lp->storedsolvals->lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE) );
         }
      }
      else
      {
         SCIP_CALL( lpRestoreSolVals(lp, blkmem, -1LL) );
      }
   }

#ifndef NDEBUG
   {
      int c;
      for( c = 0; c < lp->ncols; ++c )
      {
         assert(lp->cols[c] != NULL);
         assert(lp->cols[c]->var != NULL);
         assert(SCIPvarGetStatus(lp->cols[c]->var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(lp->cols[c]->var) == lp->cols[c]);
         assert(SCIPsetIsEQ(set, SCIPvarGetObj(lp->cols[c]->var), lp->cols[c]->obj));
         assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(lp->cols[c]->var), lp->cols[c]->lb));
         assert(SCIPsetIsEQ(set, SCIPvarGetUbLocal(lp->cols[c]->var), lp->cols[c]->ub));
      }
   }
#endif

   return SCIP_OKAY;
}

#define DIVESTACKGROWFACT 1.5

/** records a current row side such that any change will be undone after diving */
SCIP_RETCODE SCIPlpRecordOldRowSideDive(
   SCIP_LP*              lp,                 /**< LP data object */
   SCIP_ROW*             row,                /**< row affected by the change */
   SCIP_SIDETYPE         sidetype            /**< side type */
   )
{
   assert(lp != NULL);
   assert(row != NULL);

   if( lp->ndivechgsides == lp->divechgsidessize )
   {
      SCIP_CALL( reallocDiveChgSideArrays(lp, lp->divechgsidessize + 1, DIVESTACKGROWFACT) );
   }
   assert(lp->ndivechgsides < lp->divechgsidessize);

   lp->divechgsides[lp->ndivechgsides] = (sidetype == SCIP_SIDETYPE_LEFT) ? row->lhs : row->rhs;
   lp->divechgsidetypes[lp->ndivechgsides] = sidetype;
   lp->divechgrows[lp->ndivechgsides] = row;
   lp->ndivechgsides++;

   return SCIP_OKAY;
}

/** informs the LP that probing mode was initiated */
SCIP_RETCODE SCIPlpStartProbing(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(!lp->probing);
   assert(!lp->strongbranching);
   assert(!lp->strongbranchprobing);

   lp->probing = TRUE;

   return SCIP_OKAY;
}

/** informs the LP that probing mode was finished */
SCIP_RETCODE SCIPlpEndProbing(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->probing);
   assert(!lp->strongbranching);
   assert(!lp->strongbranchprobing);

   lp->probing = FALSE;

   return SCIP_OKAY;
}

/** informs the LP that the probing mode is now used for strongbranching */
void SCIPlpStartStrongbranchProbing(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->probing);
   assert(!lp->strongbranching);
   assert(!lp->strongbranchprobing);

   lp->strongbranchprobing = TRUE;
}

/** informs the LP that the probing mode is not used for strongbranching anymore */
void SCIPlpEndStrongbranchProbing(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->probing);
   assert(!lp->strongbranching);
   assert(lp->strongbranchprobing);

   lp->strongbranchprobing = FALSE;
}

/** calculates y*b + min{(c - y*A)*x | lb <= x <= ub} for given vectors y and c;
 *  the vector b is defined with b[i] = lhs[i] if y[i] >= 0, b[i] = rhs[i] if y[i] < 0
 *  Calculating this value in interval arithmetics gives a proved lower LP bound for the following reason (assuming,
 *  we have only left hand sides):
 *           min{cx       |  b <=  Ax, lb <= x <= ub}
 *   >=      min{cx       | yb <= yAx, lb <= x <= ub}   (restriction in minimum is relaxed)
 *   == yb + min{cx - yb  | yb <= yAx, lb <= x <= ub}   (added yb - yb == 0)
 *   >= yb + min{cx - yAx | yb <= yAx, lb <= x <= ub}   (because yAx >= yb inside minimum)
 *   >= yb + min{cx - yAx |            lb <= x <= ub}   (restriction in minimum is relaxed)
 */
static
SCIP_RETCODE provedBound(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             usefarkas,          /**< use y = dual Farkas and c = 0 instead of y = dual solution and c = obj? */
   SCIP_Real*            bound               /**< result of interval arithmetic minimization */
   )
{
   SCIP_INTERVAL* yinter;
   SCIP_INTERVAL b;
   SCIP_INTERVAL ytb;
   SCIP_INTERVAL prod;
   SCIP_INTERVAL diff;
   SCIP_INTERVAL x;
   SCIP_INTERVAL minprod;
   SCIP_INTERVAL a;
   SCIP_ROW* row;
   SCIP_COL* col;
   SCIP_Real y;
   SCIP_Real c;
   int i;
   int j;

   assert(lp != NULL);
   assert(lp->solved);
   assert(set != NULL);
   assert(bound != NULL);

   /* allocate buffer for storing y in interval arithmetic */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &yinter, lp->nrows) );

   /* create y vector in interval arithmetic, setting near zeros to zero; calculate y^Tb */
   SCIPintervalSet(&ytb, 0.0);
   for( j = 0; j < lp->nrows; ++j )
   {
      row = lp->rows[j];
      assert(row != NULL);

      y = (usefarkas ? row->dualfarkas : row->dualsol);

      if( SCIPsetIsFeasPositive(set, y) )
      {
         SCIPintervalSet(&yinter[j], y);
         SCIPintervalSet(&b, row->lhs - row->constant);
      }
      else if( SCIPsetIsFeasNegative(set, y) )
      {
         SCIPintervalSet(&yinter[j], y);
         SCIPintervalSet(&b, row->rhs - row->constant);
      }
      else
      {
         SCIPintervalSet(&yinter[j], 0.0);
         SCIPintervalSet(&b, 0.0);
      }

      SCIPintervalMul(SCIPsetInfinity(set), &prod, yinter[j], b);
      SCIPintervalAdd(SCIPsetInfinity(set), &ytb, ytb, prod);
   }

   /* calculate min{(c^T - y^TA)x} */
   SCIPintervalSet(&minprod, 0.0);
   for( j = 0; j < lp->ncols; ++j )
   {
      col = lp->cols[j];
      assert(col != NULL);
      assert(col->nunlinked == 0);

      SCIPintervalSetBounds(&x, SCIPcolGetLb(col), SCIPcolGetUb(col));

      c = usefarkas ? 0.0 : col->obj;
      SCIPintervalSet(&diff, c);

      for( i = 0; i < col->nlprows; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos >= 0);
         assert(col->linkpos[i] >= 0);
         SCIPintervalSet(&a, col->vals[i]);
         SCIPintervalMul(SCIPsetInfinity(set), &prod, yinter[col->rows[i]->lppos], a);
         SCIPintervalSub(SCIPsetInfinity(set), &diff, diff, prod);
      }

#ifndef NDEBUG
      for( i = col->nlprows; i < col->len; ++i )
      {
         assert(col->rows[i] != NULL);
         assert(col->rows[i]->lppos == -1);
         assert(col->rows[i]->dualsol == 0.0);
         assert(col->rows[i]->dualfarkas == 0.0);
         assert(col->linkpos[i] >= 0);
      }
#endif

      SCIPintervalSetBounds(&x, col->lb, col->ub);
      SCIPintervalMul(SCIPsetInfinity(set), &diff, diff, x);
      SCIPintervalAdd(SCIPsetInfinity(set), &minprod, minprod, diff);
   }

   /* add y^Tb */
   SCIPintervalAdd(SCIPsetInfinity(set), &minprod, minprod, ytb);

   /* free buffer for storing y in interval arithmetic */
   SCIPsetFreeBufferArray(set, &yinter);

   *bound = SCIPintervalGetInf(minprod);

   return SCIP_OKAY;
}

/** gets proven lower (dual) bound of last LP solution */
SCIP_RETCODE SCIPlpGetProvedLowerbound(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            bound               /**< pointer to store proven dual bound */
   )
{
   SCIP_CALL( provedBound(lp, set, FALSE, bound) );

   SCIPsetDebugMsg(set, "proved lower bound of LP: %.15g\n", *bound);

   return SCIP_OKAY;
}

/** gets proven dual bound of last LP solution */
SCIP_RETCODE SCIPlpIsInfeasibilityProved(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool*            proved              /**< pointer to store whether infeasibility is proven */
   )
{
   SCIP_Real bound;

   assert(proved != NULL);

   SCIP_CALL( provedBound(lp, set, TRUE, &bound) );

   *proved = (bound > 0.0);

   SCIPsetDebugMsg(set, "proved Farkas value of LP: %g -> infeasibility %sproved\n", bound, *proved ? "" : "not ");

   return SCIP_OKAY;
}



/** writes LP to a file */
SCIP_RETCODE SCIPlpWrite(
   SCIP_LP*              lp,                 /**< current LP data */
   const char*           fname               /**< file name */
   )
{
   assert(lp != NULL);
   assert(lp->flushed);
   assert(fname != NULL);

   SCIP_CALL( SCIPlpiWriteLP(lp->lpi, fname) );

   return SCIP_OKAY;
}

/** writes MIP relaxation of the current B&B node to a file */
SCIP_RETCODE SCIPlpWriteMip(
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           fname,              /**< file name */
   SCIP_Bool             genericnames,       /**< should generic names like x_i and row_j be used in order to avoid
                                              *   troubles with reserved symbols? */
   SCIP_Bool             origobj,            /**< should the original objective function be used? */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< objective scaling factor */
   SCIP_Real             objoffset,          /**< objective offset, e.g., caused by variable fixings in presolving */
   SCIP_Bool             lazyconss           /**< output removable rows as lazy constraints? */
   )
{
   FILE* file;
   int i;
   int j;
   char rowname[SCIP_MAXSTRLEN];
   SCIP_Real coeff;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(fname != NULL);

   SCIPsetDebugMsg(set, "Start to write MIP to file <%s>\n", fname);
   file = fopen(fname, "w");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for writing\n", fname);
      SCIPprintSysError(fname);
      return SCIP_FILECREATEERROR;
   }

   /* print comments */
   if( genericnames )
      SCIPmessageFPrintInfo(messagehdlr, file, "\\ Original Variable and Constraint Names have been replaced by generic names.\n"); 
   else
   {
      SCIPmessageFPrintInfo(messagehdlr, file, "\\ Warning: Variable and Constraint Names should not contain special characters like '+', '=' etc.\n");
      SCIPmessageFPrintInfo(messagehdlr, file, "\\ If this is the case, the model may be corrupted!\n");
   }

   if( origobj && objoffset != 0.0 )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, "\\ An artificial variable 'objoffset' has been added and fixed to 1.\n"); 
      SCIPmessageFPrintInfo(messagehdlr, file, "\\ Switching this variable to 0 will disable the offset in the objective.\n\n"); 
   }

   /* print objective function */
   /**@note the transformed problem in SCIP is always a minimization problem */
   if( !origobj || objsense == SCIP_OBJSENSE_MINIMIZE )
      SCIPmessageFPrintInfo(messagehdlr, file, "Minimize");
   else
      SCIPmessageFPrintInfo(messagehdlr, file, "Maximize");

   /* print objective */
   SCIPmessageFPrintInfo(messagehdlr, file, "\nObj:");
   j = 0;
   for( i = 0; i < lp->ncols; ++i )
   {
      if( lp->cols[i]->obj != 0.0 )
      {
         coeff = lp->cols[i]->obj;
         if( origobj )
         {
            coeff *= (SCIP_Real) objsense;
            coeff *= objscale;
         }

         if( genericnames )
            SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g x_%d", coeff, lp->cols[i]->lppos);
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g %s", coeff, lp->cols[i]->var->name);

         ++j;
         if( j % 10 == 0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "\n     ");
      }
   }
   /* add artificial variable 'objoffset' to transfer objective offset */
   if( origobj && objoffset != 0.0 )
      SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g objoffset", objoffset * (SCIP_Real) objsense * objscale);

   /* print constraint section */
   SCIPmessageFPrintInfo(messagehdlr, file, "\nSubject to\n");
   for( i = 0; i < lp->nrows; i++ )
   {
      char type = 'i';

      /* skip removable rows if we want to write them as lazy constraints */
      if ( lazyconss && SCIProwIsRemovable(lp->rows[i]) )
         continue;

      /* constraint types: 'l' means: only lhs exists, 'r' means: only rhs exists, 'e' means: both sides exist and are
       * equal, 'b' and 'B' mean: both sides exist, if the type is 'b', the lhs will be written, if the type is 'B', 
       * the rhs will be written. Ergo: set type to b first, change it to 'B' afterwards and go back to WRITEROW.
       * type 'i' means: lhs and rhs are both infinite */      
      if( SCIPsetIsInfinity(set, REALABS(lp->rows[i]->lhs)) && !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->rhs)) )
         type = 'r';
      else if( !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->lhs)) && SCIPsetIsInfinity(set, REALABS(lp->rows[i]->rhs)) )
         type = 'l';
      else if( !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->lhs)) && SCIPsetIsEQ(set, lp->rows[i]->lhs, lp->rows[i]->rhs) )
         type = 'e';
      else if( !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->lhs)) && !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->rhs)) )
         type = 'b';

      /* print name of row */
      if( genericnames )
         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "row_%d", lp->rows[i]->lppos);
      else
         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s", lp->rows[i]->name);

   WRITEROW:
      switch( type )
      {
      case 'r':
      case 'l':
      case 'e':
         if( strlen(rowname) > 0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "%s: ", rowname);
         break;
      case 'i':
         SCIPmessageFPrintInfo(messagehdlr, file, "\\\\ WARNING: The lhs and the rhs of the row with original name <%s>", lp->rows[i]->name);
         SCIPmessageFPrintInfo(messagehdlr, file, "are not in a valid range. The following two constraints may be corrupted!\n");
         SCIPmessagePrintWarning(messagehdlr, "The lhs and rhs of row <%s> are not in a valid range.\n", lp->rows[i]->name);
         type = 'b';
         /*lint -fallthrough*/
      case 'b':
         SCIPmessageFPrintInfo(messagehdlr, file, "%s_lhs: ", rowname);
         break;
      default:
         assert(type == 'B');
         SCIPmessageFPrintInfo(messagehdlr, file, "%s_rhs: ", rowname);
         break;
      }

      /* print coefficients and variables */
      for( j = 0; j < lp->rows[i]->nlpcols; ++j )
      {
         if( genericnames )
            SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g x_%d", lp->rows[i]->vals[j], lp->rows[i]->cols[j]->lppos);
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g %s", lp->rows[i]->vals[j], lp->rows[i]->cols[j]->var->name);

         if( (j+1) % 10 == 0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "\n          ");
      }

      /* print right hand side */
      switch( type )
      {
      case 'b':
         SCIPmessageFPrintInfo(messagehdlr, file, " >= %.15g\n", lp->rows[i]->lhs - lp->rows[i]->constant);
         type = 'B';
         goto WRITEROW;
      case 'l':
         SCIPmessageFPrintInfo(messagehdlr, file, " >= %.15g\n", lp->rows[i]->lhs - lp->rows[i]->constant);
         break;
      case 'B':
      case 'r':
         SCIPmessageFPrintInfo(messagehdlr, file, " <= %.15g\n", lp->rows[i]->rhs - lp->rows[i]->constant);
         break;
      case 'e':
         SCIPmessageFPrintInfo(messagehdlr, file, " = %.15g\n", lp->rows[i]->lhs - lp->rows[i]->constant);
         break;
      default:
         SCIPerrorMessage("Undefined row type!\n");
         return SCIP_ERROR;
      }
   }

   if ( lazyconss )
   {
      /* print lazy constraint section */
      SCIPmessageFPrintInfo(messagehdlr, file, "lazy constraints\n");
      for( i = 0; i < lp->nrows; i++ )
      {
         char type = 'i';

         /* skip non-removable rows if we want to write lazy constraints */
         if ( ! SCIProwIsRemovable(lp->rows[i]) )
            continue;

         /* constraint types: 'l' means: only lhs exists, 'r' means: only rhs exists, 'e' means: both sides exist and are
          * equal, 'b' and 'B' mean: both sides exist, if the type is 'b', the lhs will be written, if the type is 'B', 
          * the rhs will be written. Ergo: set type to b first, change it to 'B' afterwards and go back to WRITEROW.
          * type 'i' means: lhs and rhs are both infinite */      
         if( SCIPsetIsInfinity(set, REALABS(lp->rows[i]->lhs)) && !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->rhs)) )
            type = 'r';
         else if( !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->lhs)) && SCIPsetIsInfinity(set, REALABS(lp->rows[i]->rhs)) )
            type = 'l';
         else if( !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->lhs)) && SCIPsetIsEQ(set, lp->rows[i]->lhs, lp->rows[i]->rhs) )
            type = 'e';
         else if( !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->lhs)) && !SCIPsetIsInfinity(set, REALABS(lp->rows[i]->rhs)) )
            type = 'b';

         /* print name of row */
         if( genericnames )
            (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "row_%d", lp->rows[i]->lppos);
         else
            (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s", lp->rows[i]->name);

      WRITELAZYROW:
         switch( type )
         {
         case 'r':
         case 'l':
         case 'e':
            if( strlen(rowname) > 0 )
               SCIPmessageFPrintInfo(messagehdlr, file, "%s: ", rowname);
            break;
         case 'i':
            SCIPmessageFPrintInfo(messagehdlr, file, "\\\\ WARNING: The lhs and the rhs of the row with original name <%s>", lp->rows[i]->name);
            SCIPmessageFPrintInfo(messagehdlr, file, "are not in a valid range. The following two constraints may be corrupted!\n");
            SCIPmessagePrintWarning(messagehdlr, "The lhs and rhs of row <%s> are not in a valid range.\n",lp->rows[i]->name);
            type = 'b';
            /*lint -fallthrough*/
         case 'b':
            SCIPmessageFPrintInfo(messagehdlr, file, "%s_lhs: ", rowname);
            break;
         default:
            assert(type == 'B');
            SCIPmessageFPrintInfo(messagehdlr, file, "%s_rhs: ", rowname);
            break;
         }

         /* print coefficients and variables */
         for( j = 0; j < lp->rows[i]->nlpcols; ++j )
         {
            if( genericnames )
               SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g x_%d", lp->rows[i]->vals[j], lp->rows[i]->cols[j]->lppos);
            else
               SCIPmessageFPrintInfo(messagehdlr, file, " %+.15g %s", lp->rows[i]->vals[j], lp->rows[i]->cols[j]->var->name);

            if( (j+1) % 10 == 0 )
               SCIPmessageFPrintInfo(messagehdlr, file, "\n          ");
         }

         /* print right hand side */
         switch( type )
         {
         case 'b':
            SCIPmessageFPrintInfo(messagehdlr, file, " >= %.15g\n", lp->rows[i]->lhs - lp->rows[i]->constant);
            type = 'B';
            goto WRITELAZYROW;
         case 'l':
            SCIPmessageFPrintInfo(messagehdlr, file, " >= %.15g\n", lp->rows[i]->lhs - lp->rows[i]->constant);
            break;
         case 'B':
         case 'r':
            SCIPmessageFPrintInfo(messagehdlr, file, " <= %.15g\n", lp->rows[i]->rhs - lp->rows[i]->constant);
            break;
         case 'e':
            SCIPmessageFPrintInfo(messagehdlr, file, " = %.15g\n", lp->rows[i]->lhs - lp->rows[i]->constant);
            break;
         default:
            SCIPerrorMessage("Undefined row type!\n");
            return SCIP_ERROR;
         }
      }
   }

   /* print variable bounds */
   SCIPmessageFPrintInfo(messagehdlr, file, "Bounds\n");
   for( i = 0; i < lp->ncols; ++i )
   {
      if( !SCIPsetIsInfinity(set,-lp->cols[i]->lb) || !SCIPsetIsInfinity(set,lp->cols[i]->ub) )
      {
         /* print lower bound as far this one is not infinity */
         if( !SCIPsetIsInfinity(set,-lp->cols[i]->lb) )
            SCIPmessageFPrintInfo(messagehdlr, file, " %.15g <=", lp->cols[i]->lb);

         /* print variable name */
         if( genericnames )
            SCIPmessageFPrintInfo(messagehdlr, file, " x_%d ", lp->cols[i]->lppos);
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %s ", lp->cols[i]->var->name);

         /* print upper bound as far this one is not infinity */
         if( !SCIPsetIsInfinity(set,lp->cols[i]->ub) )
            SCIPmessageFPrintInfo(messagehdlr, file, "<= %.15g", lp->cols[i]->ub);
         SCIPmessageFPrintInfo(messagehdlr, file, "\n");
      }
   }
   if( origobj && objoffset != 0.0 )
      SCIPmessageFPrintInfo(messagehdlr, file, " objoffset = 1\n");

   /* print integer variables */
   SCIPmessageFPrintInfo(messagehdlr, file, "Generals\n");
   j = 0;
   for( i = 0; i < lp->ncols; ++i )
   {
      if( SCIPvarIsIntegral(lp->cols[i]->var) )
      {
         /* print variable name */
         if( genericnames )
            SCIPmessageFPrintInfo(messagehdlr, file, " x_%d ", lp->cols[i]->lppos);
         else
            SCIPmessageFPrintInfo(messagehdlr, file, " %s ", lp->cols[i]->var->name);

         j++;
         if( j % 10 == 0 )
            SCIPmessageFPrintInfo(messagehdlr, file, "\n");
      }
   }

   SCIPmessageFPrintInfo(messagehdlr, file, "\nEnd");
   fclose(file);

   return SCIP_OKAY;
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPcolGetObj
#undef SCIPcolGetLb
#undef SCIPcolGetUb
#undef SCIPcolGetBestBound
#undef SCIPcolGetPrimsol
#undef SCIPcolGetMinPrimsol
#undef SCIPcolGetMaxPrimsol
#undef SCIPcolGetBasisStatus
#undef SCIPcolGetVar
#undef SCIPcolGetIndex
#undef SCIPcolIsIntegral
#undef SCIPcolIsRemovable
#undef SCIPcolGetLPPos
#undef SCIPcolGetLPDepth
#undef SCIPcolIsInLP
#undef SCIPcolGetNNonz
#undef SCIPcolGetNLPNonz
#undef SCIPcolGetRows
#undef SCIPcolGetVals
#undef SCIPcolGetStrongbranchNode
#undef SCIPcolGetNStrongbranchs
#undef SCIPboundtypeOpposite
#undef SCIProwGetNNonz
#undef SCIProwGetNLPNonz
#undef SCIProwGetCols
#undef SCIProwGetVals
#undef SCIProwGetConstant
#undef SCIProwGetNorm
#undef SCIProwGetSumNorm
#undef SCIProwGetLhs
#undef SCIProwGetRhs
#undef SCIProwGetDualsol
#undef SCIProwGetDualfarkas
#undef SCIProwGetBasisStatus
#undef SCIProwGetName
#undef SCIProwGetIndex
#undef SCIProwGetAge
#undef SCIProwGetRank
#undef SCIProwIsIntegral
#undef SCIProwIsLocal
#undef SCIProwIsModifiable
#undef SCIProwIsRemovable
#undef SCIProwGetOrigintype
#undef SCIProwGetOriginCons
#undef SCIProwGetOriginSepa
#undef SCIProwIsInGlobalCutpool
#undef SCIProwGetLPPos
#undef SCIProwGetLPDepth
#undef SCIProwIsInLP
#undef SCIProwGetActiveLPCount
#undef SCIProwGetNLPsAfterCreation
#undef SCIProwChgRank
#undef SCIPlpGetCols
#undef SCIPlpGetNCols
#undef SCIPlpGetRows
#undef SCIPlpGetNRows
#undef SCIPlpGetNewcols
#undef SCIPlpGetNNewcols
#undef SCIPlpGetNewrows
#undef SCIPlpGetNNewrows
#undef SCIPlpGetObjNorm
#undef SCIPlpGetRootObjval
#undef SCIPlpGetRootColumnObjval
#undef SCIPlpGetRootLooseObjval
#undef SCIPlpGetLPI
#undef SCIPlpSetIsRelax
#undef SCIPlpIsRelax
#undef SCIPlpIsSolved
#undef SCIPlpIsSolBasic
#undef SCIPlpDiving
#undef SCIPlpDivingObjChanged
#undef SCIPlpMarkDivingObjChanged
#undef SCIPlpUnmarkDivingObjChanged
#undef SCIPlpDivingRowsChanged

/** gets objective value of column */
SCIP_Real SCIPcolGetObj(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->obj;
}

/** gets lower bound of column */
SCIP_Real SCIPcolGetLb(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->lb;
}

/** gets upper bound of column */
SCIP_Real SCIPcolGetUb(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->ub;
}

/** gets best bound of column with respect to the objective function */
SCIP_Real SCIPcolGetBestBound(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->obj >= 0.0 )
      return col->lb;
   else
      return col->ub;
}

/** gets the primal LP solution of a column */
SCIP_Real SCIPcolGetPrimsol(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   if( col->lppos >= 0 )
      return col->primsol;
   else
      return 0.0;
}

/** gets the minimal LP solution value, this column ever assumed */
SCIP_Real SCIPcolGetMinPrimsol(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->minprimsol;
}

/** gets the maximal LP solution value, this column ever assumed */
SCIP_Real SCIPcolGetMaxPrimsol(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->maxprimsol;
}

/** gets the basis status of a column in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_ZERO for columns not in the current SCIP_LP
 */
SCIP_BASESTAT SCIPcolGetBasisStatus(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(col->lppos >= 0 || (SCIP_BASESTAT)col->basisstatus == SCIP_BASESTAT_ZERO);

   return (SCIP_BASESTAT)col->basisstatus;
}

/** gets variable this column represents */
SCIP_VAR* SCIPcolGetVar(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->var;
}

/** gets unique index of col */
int SCIPcolGetIndex(
   SCIP_COL*             col                 /**< LP col */
   )
{
   assert(col != NULL);

   return col->index;
}

/** returns whether the associated variable is of integral type (binary, integer, implicit integer) */
SCIP_Bool SCIPcolIsIntegral(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(SCIPvarIsIntegral(col->var) == col->integral);

   return col->integral;
}

/** returns TRUE iff column is removable from the LP (due to aging or cleanup) */
SCIP_Bool SCIPcolIsRemovable(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->removable;
}

/** gets position of column in current LP, or -1 if it is not in LP */
int SCIPcolGetLPPos(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return col->lppos;
}

/** gets depth in the tree where the column entered the LP, or -1 if it is not in LP */
int SCIPcolGetLPDepth(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return col->lpdepth;
}

/** returns TRUE iff column is member of current LP */
SCIP_Bool SCIPcolIsInLP(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert((col->lppos == -1) == (col->lpdepth == -1));

   return (col->lppos >= 0);
}

/** get number of nonzero entries in column vector */
int SCIPcolGetNNonz(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->len;
}

/** get number of nonzero entries in column vector, that correspond to rows currently in the SCIP_LP;
 *
 *  @warning This method is only applicable on columns, that are completely linked to their rows (e.g. a column
 *  that is in the current LP and the LP was solved, or a column that was in a solved LP and didn't change afterwards
 */
int SCIPcolGetNLPNonz(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);
   assert(col->nunlinked == 0);

   return col->nlprows;
}

/** gets array with rows of nonzero entries */
SCIP_ROW** SCIPcolGetRows(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->rows;
}

/** gets array with coefficients of nonzero entries */
SCIP_Real* SCIPcolGetVals(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->vals;
}

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
SCIP_Longint SCIPcolGetStrongbranchNode(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->sbnode;
}

/** gets number of times, strong branching was applied in current run on the given column */
int SCIPcolGetNStrongbranchs(
   SCIP_COL*             col                 /**< LP column */
   )
{
   assert(col != NULL);

   return col->nsbcalls;
}

/** gets opposite bound type of given bound type */
SCIP_BOUNDTYPE SCIPboundtypeOpposite(
   SCIP_BOUNDTYPE        boundtype           /**< type of bound (lower or upper) */
   )
{
   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

   return (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);
}

/** get number of nonzero entries in row vector */
int SCIProwGetNNonz(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->len;
}

/** get number of nonzero entries in row vector, that correspond to columns currently in the SCIP_LP;
 *
 *  @warning This method is only applicable on rows, that are completely linked to their columns (e.g. a row
 *  that is in the current LP and the LP was solved, or a row that was in a solved LP and didn't change afterwards
 */
int SCIProwGetNLPNonz(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->nunlinked == 0);

   return row->nlpcols;
}

/** gets array with columns of nonzero entries */
SCIP_COL** SCIProwGetCols(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->cols;
}

/** gets array with coefficients of nonzero entries */
SCIP_Real* SCIProwGetVals(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->vals;
}

/** gets constant shift of row */
SCIP_Real SCIProwGetConstant(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->constant;
}

/** gets Euclidean norm of row vector */
SCIP_Real SCIProwGetNorm(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   checkRowSqrnorm(row);

   return sqrt(row->sqrnorm);
}

/** gets sum norm of row vector (sum of absolute values of coefficients) */
SCIP_Real SCIProwGetSumNorm(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   checkRowSumnorm(row);

   return row->sumnorm;
}

/** returns the left hand side of the row */
SCIP_Real SCIProwGetLhs(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->lhs;
}

/** returns the right hand side of the row */
SCIP_Real SCIProwGetRhs(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->rhs;
}

/** gets the dual LP solution of a row */
SCIP_Real SCIProwGetDualsol(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   if( row->lppos >= 0 )
      return row->dualsol;
   else
      return 0.0;
}

/** gets the dual Farkas coefficient of a row in an infeasible LP */
SCIP_Real SCIProwGetDualfarkas(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   if( row->lppos >= 0 )
      return row->dualfarkas;
   else
      return 0.0;
}

/** gets the basis status of a row in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_BASIC for rows not in the current SCIP_LP
 */
SCIP_BASESTAT SCIProwGetBasisStatus(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert(row->lppos >= 0 || (SCIP_BASESTAT)row->basisstatus == SCIP_BASESTAT_BASIC);

   return (SCIP_BASESTAT)row->basisstatus;
}

/** returns the name of the row */
const char* SCIProwGetName(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->name;
}

/** gets unique index of row */
int SCIProwGetIndex(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->index;
}

/** gets age of row */
int SCIProwGetAge(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->age;
}

/** gets rank of row */
int SCIProwGetRank(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->rank;
}

/** returns TRUE iff the activity of the row (without the row's constant) is always integral in a feasible solution */
SCIP_Bool SCIProwIsIntegral(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->integral;
}

/** returns TRUE iff row is only valid locally */
SCIP_Bool SCIProwIsLocal(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->local;
}

/** returns TRUE iff row is modifiable during node processing (subject to column generation) */
SCIP_Bool SCIProwIsModifiable(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->modifiable;
}

/** returns TRUE iff row is removable from the LP (due to aging or cleanup) */
SCIP_Bool SCIProwIsRemovable(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->removable;
}

/** returns type of origin that created the row */
SCIP_ROWORIGINTYPE SCIProwGetOrigintype(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert( row != NULL );

   return (SCIP_ROWORIGINTYPE) row->origintype;
}

/** returns origin constraint handler that created the row (NULL if not available) */
SCIP_CONSHDLR* SCIProwGetOriginCons(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert( row != NULL );

   if ( (SCIP_ROWORIGINTYPE) row->origintype == SCIP_ROWORIGINTYPE_CONS )
   {
      assert( row->origin != NULL );
      return (SCIP_CONSHDLR*) row->origin;
   }
   return NULL;
}

/** returns origin separator that created the row (NULL if not available) */
SCIP_SEPA* SCIProwGetOriginSepa(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert( row != NULL );

   if ( (SCIP_ROWORIGINTYPE) row->origintype == SCIP_ROWORIGINTYPE_SEPA )
   {
      assert( row->origin != NULL );
      return (SCIP_SEPA*) row->origin;
   }
   return NULL;
}

/** returns TRUE iff row is member of the global cut pool */
SCIP_Bool SCIProwIsInGlobalCutpool(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);

   return row->inglobalcutpool;
}

/** gets position of row in current LP, or -1 if it is not in LP */
int SCIProwGetLPPos(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return row->lppos;
}

/** gets depth in the tree where the row entered the LP, or -1 if it is not in LP */
int SCIProwGetLPDepth(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return row->lpdepth;
}

/** returns TRUE iff row is member of current LP */
SCIP_Bool SCIProwIsInLP(
   SCIP_ROW*             row                 /**< LP row */
   )
{
   assert(row != NULL);
   assert((row->lppos == -1) == (row->lpdepth == -1));

   return (row->lppos >= 0);
}

/** changes the rank of LP row */
void SCIProwChgRank(
   SCIP_ROW*             row,                /**< LP row */
   int                   rank                /**< new value for rank */
   )
{
   assert(row != NULL);

   row->rank = rank;
}

/** returns the number of times that this row has been sharp in an optimal LP solution */
SCIP_Longint SCIProwGetActiveLPCount(
   SCIP_ROW*             row                 /**< row */
   )
{
   assert(row != NULL);

   return row->activeinlpcounter;
}

/** returns the number of LPs since this row has been created */
SCIP_Longint SCIProwGetNLPsAfterCreation(
   SCIP_ROW*             row                 /**< row */
   )
{
   assert(row != NULL);

   return row->nlpsaftercreation;
}

/** gets array with columns of the LP */
SCIP_COL** SCIPlpGetCols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->cols;
}

/** gets current number of columns in LP */
int SCIPlpGetNCols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->ncols;
}

/** gets array with rows of the LP */
SCIP_ROW** SCIPlpGetRows(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->rows;
}

/** gets current number of rows in LP */
int SCIPlpGetNRows(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->nrows;
}

/** gets array with newly added columns after the last mark */
SCIP_COL** SCIPlpGetNewcols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return &(lp->cols[lp->firstnewcol]);
}

/** gets number of newly added columns after the last mark */
int SCIPlpGetNNewcols(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewcol && lp->firstnewcol <= lp->ncols);

   return lp->ncols - lp->firstnewcol;
}

/** gets array with newly added rows after the last mark */
SCIP_ROW** SCIPlpGetNewrows(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return &(lp->rows[lp->firstnewrow]);
}

/** gets number of newly added rows after the last mark */
int SCIPlpGetNNewrows(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(0 <= lp->firstnewrow && lp->firstnewrow <= lp->nrows);

   return lp->nrows - lp->firstnewrow;
}

/** recalculates Euclidean norm of objective function vector of column variables if it have gotten unreliable during calculation */
void SCIPlpRecalculateObjSqrNorm(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< LP data */
   )
{
   if( lp->objsqrnormunreliable )
   {
      SCIP_COL** cols;
      int c;

      cols = lp->cols;
      assert(cols != NULL || lp->ncols == 0);

      lp->objsqrnorm = 0.0;

      for( c = lp->ncols - 1; c >= 0; --c )
      {
         lp->objsqrnorm += SQR(cols[c]->unchangedobj);  /*lint !e613*/
      }
      assert(SCIPsetIsGE(set, lp->objsqrnorm, 0.0));

      /* due to numerical troubles it still can appear that lp->objsqrnorm is a little bit smaller than 0 */
      lp->objsqrnorm = MAX(lp->objsqrnorm, 0.0);

      lp->objsqrnormunreliable = FALSE;
   }
   return;
}

/** gets Euclidean norm of objective function vector of column variables, only use this method if
 *  lp->objsqrnormunreliable == FALSE, so probably you have to call SCIPlpRecalculateObjSqrNorm before */
SCIP_Real SCIPlpGetObjNorm(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);
   assert(!lp->objsqrnormunreliable);
   assert(lp->objsqrnorm >= 0.0);

   return SQRT(lp->objsqrnorm);
}

/** sets whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound */
void SCIPlpSetRootLPIsRelax(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             isrelax             /**< is the root lp a relaxation of the problem? */
   )
{
   assert(lp != NULL);

   lp->rootlpisrelax = isrelax;
}

/** returns whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound */
SCIP_Bool SCIPlpIsRootLPRelax(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   return lp->rootlpisrelax;
}

/** gets the objective value of the root node LP; returns SCIP_INVALID if the root node LP was not (yet) solved */
SCIP_Real SCIPlpGetRootObjval(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   return MIN(lp->rootlpobjval + lp->rootlooseobjval, SCIP_INVALID);
}

/** gets part of the objective value of the root node LP that results from COLUMN variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 */
SCIP_Real SCIPlpGetRootColumnObjval(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   return lp->rootlpobjval;
}

/** gets part of the objective value of the root node LP that results from LOOSE variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 */
SCIP_Real SCIPlpGetRootLooseObjval(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   return lp->rootlooseobjval;
}

/** gets the LP solver interface */
SCIP_LPI* SCIPlpGetLPI(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->lpi;
}

/** sets whether the current LP is a relaxation of the current problem and its optimal objective value is a local lower bound */
void SCIPlpSetIsRelax(
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             relax               /**< is the current lp a relaxation? */
   )
{
   assert(lp != NULL);

   lp->isrelax = relax;
}

/** returns whether the current LP is a relaxation of the problem for which it has been solved and its 
 *  solution value a valid local lower bound? 
 */
SCIP_Bool SCIPlpIsRelax(
   SCIP_LP*              lp                  /**< LP data */
   )
{
   assert(lp != NULL);

   return lp->isrelax;
}

/** returns whether the current LP is flushed and solved */
SCIP_Bool SCIPlpIsSolved(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->flushed && lp->solved;
}

/** return whether the current LP solution passed the primal feasibility check */
SCIP_Bool SCIPlpIsPrimalReliable(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return (lp->primalchecked && lp->primalfeasible);
}

/** return whether the current LP solution passed the dual feasibility check */
SCIP_Bool SCIPlpIsDualReliable(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return (lp->dualchecked && lp->dualfeasible);
}

/** returns whether the current LP solution is a basic solution */
SCIP_Bool SCIPlpIsSolBasic(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->solisbasic;
}

/** returns whether the LP is in diving mode */
SCIP_Bool SCIPlpDiving(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->diving;
}

/** returns whether the LP is in diving mode and the objective value of at least one column was changed */
SCIP_Bool SCIPlpDivingObjChanged(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);

   return lp->divingobjchg;
}

/** marks the diving LP to have a changed objective function */
void SCIPlpMarkDivingObjChanged(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->diving || lp->probing);

   lp->divingobjchg = TRUE;
}

/** marks the diving LP to not have a changed objective function anymore */
void SCIPlpUnmarkDivingObjChanged(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->diving || lp->probing);

   lp->divingobjchg = FALSE;
}

/* returns TRUE if at least one left/right hand side of an LP row was changed during diving mode */
SCIP_Bool SCIPlpDivingRowsChanged(
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(lp != NULL);
   assert(lp->diving || lp->ndivechgsides == 0);

   return (lp->ndivechgsides > 0);
}

/** compute relative interior point with auxiliary lpi, see SCIPlpComputeRelIntPoint() */
static
SCIP_RETCODE computeRelIntPoint(
   SCIP_LPI*             lpi,                /**< auxiliary LP interface */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             relaxrows,          /**< should the rows be relaxed */
   SCIP_Bool             inclobjcutoff,      /**< should a row for the objective cutoff be included */
   SCIP_Real             timelimit,          /**< time limit for LP solver */
   int                   iterlimit,          /**< iteration limit for LP solver */
   SCIP_Real*            point,              /**< array to store relative interior point on exit */
   SCIP_Bool*            success             /**< buffer to indicate whether interior point was successfully computed */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Real* primal;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* matvals;
   SCIP_Real* matlhs;
   SCIP_Real* matrhs;
   SCIP_Real objval;
   SCIP_Real alpha;
   int* matinds;
   int* matbeg;
#ifndef NDEBUG
   int nslacks;
#endif
   int nnewcols;
   int ntotnonz = 0;
   int ntotrows = 0;
   int matrowidx;
   int matidx;
   int cnt;
   int j;
   int i;

   assert(lpi != NULL);

   retcode = SCIPlpiSetRealpar(lpi, SCIP_LPPAR_FEASTOL, SCIPsetLpfeastol(set));
   if( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERUNKNOWN )
   {
      SCIP_CALL( retcode );
   }

   retcode = SCIPlpiSetRealpar(lpi, SCIP_LPPAR_DUALFEASTOL, SCIPsetDualfeastol(set));
   if( retcode != SCIP_OKAY && retcode != SCIP_PARAMETERUNKNOWN )
   {
      SCIP_CALL( retcode );
   }

   /* get storage */
   nnewcols = 3*lp->ncols + 2*lp->nrows + (inclobjcutoff ? 1 : 0) + 1;
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lb, nnewcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ub, nnewcols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &obj, nnewcols) );

   /* create original columns (bounds are relaxed below, unless the variable is fixed) */
   for( j = 0; j < lp->ncols; ++j )
   {
      /* note: if the variable is fixed we cannot simply fix the variables (because alpha scales the problem) */
      obj[j] = 0.0;
      lb[j] = -SCIPlpiInfinity(lpi);
      ub[j] =  SCIPlpiInfinity(lpi);
      /* note: we could also use the original bounds - free variables seem to be faster. */
   }

   /* add artificial alpha variable */
   nnewcols = lp->ncols;
   obj[nnewcols] = 0.0;
   lb[nnewcols] = 1.0;
   ub[nnewcols] = SCIPlpiInfinity(lpi);
   ++nnewcols;

   /* create slacks for rows */
   for( i = 0; i < lp->nrows; ++i )
   {
      SCIP_ROW* row;

      row = lp->rows[i];
      assert( row != NULL );

      if( SCIProwIsModifiable(row) )
         continue;

      /* make sure row is sorted */
      rowSortLP(row);
      assert( row->lpcolssorted );

      /* check whether we have an equation */
      if( SCIPsetIsEQ(set, row->lhs, row->rhs) )
      {
         assert( !SCIPsetIsInfinity(set, REALABS(row->lhs)) );
         assert( !SCIPsetIsInfinity(set, REALABS(row->rhs)) );
         ntotnonz += row->nlpcols + 1;
         ++ntotrows;
      }
      else
      {
         /* otherwise add slacks for each side if necessary */
         if ( ! SCIPsetIsInfinity(set, REALABS(row->lhs)) )
         {
            if ( relaxrows )
            {
               lb[nnewcols] = 0.0;
               ub[nnewcols] = 1.0;
               obj[nnewcols++] = 1.0;
               ntotnonz += row->nlpcols + 2;
            }
            else
               ntotnonz += row->nlpcols + 1;
            ++ntotrows;
         }
         if ( ! SCIPsetIsInfinity(set, REALABS(row->rhs)) )
         {
            if ( relaxrows )
            {
               lb[nnewcols] = 0.0;
               ub[nnewcols] = 1.0;
               obj[nnewcols++] = 1.0;
               ntotnonz += row->nlpcols + 2;
            }
            else
               ntotnonz += row->nlpcols + 1;
            ++ntotrows;
         }
      }
   }

   /* create slacks for objective cutoff row */
   if( inclobjcutoff && relaxrows )
   {
      /* add slacks for right hand side */
      lb[nnewcols] = 0.0;
      ub[nnewcols] = 1.0;
      obj[nnewcols++] = 1.0;
      ntotnonz += lp->ncols + 2;
      ++ntotrows;
   }

   /* create slacks for bounds */
   for( j = 0; j < lp->ncols; ++j )
   {
      SCIP_COL* col;

      col = lp->cols[j];
      assert( col != NULL );

      /* no slacks for fixed variables */
      if( SCIPsetIsEQ(set, col->lb, col->ub) )
      {
         ++ntotrows;
         ntotnonz += 2;
      }
      else
      {
         /* add slacks for each bound if necessary */
         if ( ! SCIPsetIsInfinity(set, REALABS(col->lb)) )
         {
            lb[nnewcols] = 0.0;
            ub[nnewcols] = 1.0;
            obj[nnewcols++] = 1.0;
            ntotnonz += 3;
            ++ntotrows;
         }
         if( ! SCIPsetIsInfinity(set, REALABS(col->ub)) )
         {
            lb[nnewcols] = 0.0;
            ub[nnewcols] = 1.0;
            obj[nnewcols++] = 1.0;
            ntotnonz += 3;
            ++ntotrows;
         }
      }
   }
#ifndef NDEBUG
   nslacks = nnewcols - lp->ncols - 1;
   assert( nslacks >= 0 );
   assert( nnewcols <= 3*lp->ncols + 2*lp->nrows + (inclobjcutoff ? 1 : 0) + 1 );
#endif

   /* add columns */
   SCIP_CALL( SCIPlpiAddCols(lpi, nnewcols, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );

   /* free storage */
   SCIPsetFreeBufferArray(set, &obj);
   SCIPsetFreeBufferArray(set, &ub);
   SCIPsetFreeBufferArray(set, &lb);

   /* prepare storage for rows */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &matinds, ntotnonz) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &matvals, ntotnonz) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &matbeg, ntotrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &matlhs, ntotrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &matrhs, ntotrows) );

   /* create rows arising from original rows */
   cnt = 0;
   matrowidx = 0;
   matidx = 0;
   for( i = 0; i < lp->nrows; ++i )
   {
      SCIP_ROW* row;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nnonz;

      row = lp->rows[i];
      assert( row != NULL );

      if( SCIProwIsModifiable(row) )
         continue;
      assert( row->lpcolssorted );

      /* get row data */
      lhs = row->lhs - (SCIPsetIsInfinity(set, -row->lhs) ? 0.0 : row->constant);
      rhs = row->rhs - (SCIPsetIsInfinity(set,  row->rhs) ? 0.0 : row->constant);
      nnonz = row->nlpcols;
      assert( nnonz <= lp->ncols );
      rowcols = row->cols;
      rowvals = row->vals;

      /* if we have an equation */
      if( SCIPsetIsEQ(set, lhs, rhs) )
      {
         /* set up indices */
         matbeg[matrowidx] = matidx;
         for( j = 0; j < nnonz; ++j )
         {
            assert( rowcols[j] != NULL );
            assert( 0 <= rowcols[j]->lppos && rowcols[j]->lppos < lp->ncols );
            assert( lp->cols[rowcols[j]->lppos] == rowcols[j] );
            assert( ! SCIPsetIsZero(set, rowvals[j]) );
            matinds[matidx] = rowcols[j]->lppos;
            matvals[matidx++] = rowvals[j];
            assert( matidx <= ntotnonz );
         }

         /* add artificial variable */
         if ( ! SCIPsetIsZero(set, rhs) )
         {
            matinds[matidx] = lp->ncols;
            matvals[matidx++] = -rhs;
            assert( matidx <= ntotnonz );
         }

         matlhs[matrowidx] = 0.0;
         matrhs[matrowidx++] = 0.0;
         assert( matrowidx <= ntotrows );
      }
      else
      {
         SCIP_Real abslhs = REALABS(lhs);
         SCIP_Real absrhs = REALABS(rhs);

         assert(!SCIPsetIsEQ(set, lhs, rhs));

         /* treat lhs */
         if( !SCIPsetIsInfinity(set, abslhs) )
         {
            /* set up indices */
            matbeg[matrowidx] = matidx;
            for( j = 0; j < nnonz; ++j )
            {
               assert( rowcols[j] != NULL );
               assert( 0 <= rowcols[j]->lppos && rowcols[j]->lppos < lp->ncols );
               assert( lp->cols[rowcols[j]->lppos] == rowcols[j] );
               assert( ! SCIPsetIsZero(set, rowvals[j]) );
               matinds[matidx] = rowcols[j]->lppos;
               matvals[matidx++] = rowvals[j];
               assert( matidx <= ntotnonz );
            }

            /* add artificial variable */
            if ( ! SCIPsetIsZero(set, lhs) )
            {
               matinds[matidx] = lp->ncols;
               matvals[matidx++] = -lhs;
               assert( matidx <= ntotnonz );
            }

            if( relaxrows )
            {
               /* add slack variable */
               matvals[matidx] = -MAX(1.0, lhs);      /*lint !e679*/
               matinds[matidx++] = lp->ncols + 1 + cnt; /*lint !e679*/
               assert( matidx <= ntotnonz );
               ++cnt;
            }

            matlhs[matrowidx] = 0.0;
            matrhs[matrowidx++] = SCIPlpiInfinity(lpi);
            assert( matrowidx <= ntotrows );
         }

         /* treat rhs */
         if( !SCIPsetIsInfinity(set, absrhs) )
         {
            /* set up indices */
            matbeg[matrowidx] = matidx;
            for( j = 0; j < nnonz; ++j )
            {
               assert( rowcols[j] != NULL );
               assert( 0 <= rowcols[j]->lppos && rowcols[j]->lppos < lp->ncols );
               assert( lp->cols[rowcols[j]->lppos] == rowcols[j] );
               assert( ! SCIPsetIsZero(set, rowvals[j]) );
               matinds[matidx] = rowcols[j]->lppos;
               matvals[matidx++] = rowvals[j];
               assert( matidx <= ntotnonz );
            }

            /* add artificial variable */
            if ( ! SCIPsetIsZero(set, rhs) )
            {
               matinds[matidx] = lp->ncols;
               matvals[matidx++] = -rhs;
               assert( matidx <= ntotnonz );
            }

            if( relaxrows )
            {
               /* add slack variable */
               matvals[matidx] = MAX(1.0, absrhs);    /*lint !e679*/
               matinds[matidx++] = lp->ncols + 1 + cnt; /*lint !e679*/
               ++cnt;
            }

            matlhs[matrowidx] = -SCIPlpiInfinity(lpi);
            matrhs[matrowidx++] = 0.0;
            assert( matrowidx <= ntotrows );
         }
      }
   }

   /* create row arising from objective cutoff */
   if( inclobjcutoff )
   {
      SCIP_Real rhs;

      /* get row data */
      assert(lp->looseobjvalinf == 0);
      rhs = lp->cutoffbound - getFiniteLooseObjval(lp, set, prob);

      /* set up indices and coefficients */
      matbeg[matrowidx] = matidx;
      for( j = 0; j < lp->ncols; ++j )
      {
         assert( lp->cols[j] != NULL );
         assert( 0 <= lp->cols[j]->lppos && lp->cols[j]->lppos < lp->ncols );
         assert( lp->cols[lp->cols[j]->lppos] == lp->cols[j] );

         if( ! SCIPsetIsZero(set, lp->cols[j]->obj) )
         {
            matinds[matidx] = lp->cols[j]->lppos;
            matvals[matidx++] = lp->cols[j]->obj;
            assert( matidx <= ntotnonz );
         }
      }

      /* treat rhs */

      /* add artificial variable */
      if ( ! SCIPsetIsZero(set, rhs) )
      {
         matinds[matidx] = lp->ncols;
         matvals[matidx++] = -rhs;
         assert( matidx <= ntotnonz );
      }

      if( relaxrows )
      {
         SCIP_Real absrhs = REALABS(rhs);

         /* add slack variable */
         matvals[matidx] = MAX(1.0, absrhs);
         matinds[matidx++] = lp->ncols + 1 + cnt;
         assert( matidx <= ntotnonz );
         ++cnt;
      }
      matlhs[matrowidx] = -SCIPsetInfinity(set);
      matrhs[matrowidx++] = 0.0;
      assert( matrowidx <= ntotrows );
   }

   /* create rows arising from bounds */
   for( j = 0; j < lp->ncols; ++j )
   {
      SCIP_COL* col;
      SCIP_Real abscollb;
      SCIP_Real abscolub;

      col = lp->cols[j];
      assert( col != NULL );
      assert( col->lppos == j );

      /* fixed variable */
      if( SCIPsetIsEQ(set, col->lb, col->ub) )
      {
         /* set up index of column */
         matbeg[matrowidx] = matidx;

         matinds[matidx] = j;
         matvals[matidx++] = 1.0;
         assert( matidx <= ntotnonz );

         /* add artificial variable */
         if ( ! SCIPsetIsZero(set, col->ub) )
         {
            matinds[matidx] = lp->ncols;
            matvals[matidx++] = -col->lb;
            assert( matidx <= ntotnonz );
         }

         matlhs[matrowidx] = 0.0;
         matrhs[matrowidx++] = 0.0;
         assert( matrowidx <= ntotrows );

         continue;
      }

      abscollb = REALABS(col->lb);
      abscolub = REALABS(col->ub);

      /* lower bound */
      if ( ! SCIPsetIsInfinity(set, abscollb) )
      {
         /* set up index of column */
         matbeg[matrowidx] = matidx;

         matinds[matidx] = j;
         matvals[matidx++] = 1.0;
         assert( matidx <= ntotnonz );

         /* add artificial variable */
         if ( ! SCIPsetIsZero(set, col->lb) )
         {
            matinds[matidx] = lp->ncols;
            matvals[matidx++] = -col->lb;
            assert( matidx <= ntotnonz );
         }

         /* add slack variable */
         matvals[matidx] = -MAX(1.0, abscollb);
         matinds[matidx++] = lp->ncols + 1 + cnt;
         assert( matidx <= ntotnonz );
         ++cnt;

         matlhs[matrowidx] = 0.0;
         matrhs[matrowidx++] = SCIPsetInfinity(set);
         assert( matrowidx <= ntotrows );
      }

      /* upper bound */
      if ( ! SCIPsetIsInfinity(set, abscolub) )
      {
         /* set up index of column */
         matbeg[matrowidx] = matidx;

         matinds[matidx] = j;
         matvals[matidx++] = 1.0;
         assert( matidx <= ntotnonz );

         /* add artificial variable */
         if ( ! SCIPsetIsZero(set, col->ub) )
         {
            matinds[matidx] = lp->ncols;
            matvals[matidx++] = -col->ub;
            assert( matidx <= ntotnonz );
         }

         /* add slack variable */
         matvals[matidx] = MAX(1.0, abscolub);
         matinds[matidx++] = lp->ncols + 1 + cnt;
         assert( matidx <= ntotnonz );
         ++cnt;

         matlhs[matrowidx] = -SCIPsetInfinity(set);
         matrhs[matrowidx++] = 0.0;
         assert( matrowidx <= ntotrows );
      }
   }
   assert( cnt == nslacks );
   assert( matrowidx == ntotrows );

   /* add rows */
   SCIP_CALL( SCIPlpiAddRows(lpi, ntotrows, matlhs, matrhs, NULL, matidx, matbeg, matinds, matvals) );

   SCIPsetFreeBufferArray(set, &matrhs);
   SCIPsetFreeBufferArray(set, &matlhs);
   SCIPsetFreeBufferArray(set, &matbeg);
   SCIPsetFreeBufferArray(set, &matvals);
   SCIPsetFreeBufferArray(set, &matinds);

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPlpiWriteLP(lpi, "relativeInterior.lp") );
#endif

#ifndef NDEBUG
   {
      int ncols;
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
      assert( ncols == nnewcols );
   }
#endif

   /* set time limit */
   if( SCIPsetIsInfinity(set, timelimit) )
      timelimit = SCIPlpiInfinity(lpi);
   retcode = SCIPlpiSetRealpar(lpi, SCIP_LPPAR_LPTILIM, timelimit);

   /* check, if parameter is unknown */
   if ( retcode == SCIP_PARAMETERUNKNOWN )
      SCIPmessagePrintWarning(messagehdlr, "Could not set time limit of LP solver for relative interior point computation.\n");

   /* set iteration limit */
   retcode = SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, iterlimit);

   /* check, if parameter is unknown */
   if ( retcode == SCIP_PARAMETERUNKNOWN )
      SCIPmessagePrintWarning(messagehdlr, "Could not set iteration limit of LP solver for relative interior point computation.\n");

   /* solve and store point */
   /* SCIP_CALL( SCIPlpiSolvePrimal(lpi) ); */
   SCIP_CALL( SCIPlpiSolveDual(lpi) );   /* dual is usually faster */

#ifndef NDEBUG
   if ( SCIPlpiIsIterlimExc(lpi) )
      SCIPmessagePrintWarning(messagehdlr, "Iteration limit exceeded in relative interior point computation.\n");
   if ( SCIPlpiIsTimelimExc(lpi) )
      SCIPmessagePrintWarning(messagehdlr, "Time limit exceeded in relative interior point computation.\n");
#endif

   if( SCIPlpiIsOptimal(lpi) )
   {
      /* get primal solution */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &primal, nnewcols) );
      SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primal, NULL, NULL, NULL) );
      alpha = primal[lp->ncols];
      assert( SCIPsetIsFeasGE(set, alpha, 1.0) );

      SCIPsetDebugMsg(set, "Solved relative interior lp with objective %g.\n", objval);

      /* construct relative interior point */
      for( j = 0; j < lp->ncols; ++j )
         point[j] = primal[j]/alpha;

#ifdef SCIP_DEBUG
      /* check whether the point is a relative interior point */
      cnt = 0;
      if( relaxrows )
      {
         for( i = 0; i < lp->nrows; ++i )
         {
            SCIP_ROW* row;
            SCIP_COL** rowcols;
            SCIP_Real* rowvals;
            SCIP_Real lhs;
            SCIP_Real rhs;
            SCIP_Real sum;
            int nnonz;

            row = lp->rows[i];
            assert( row != NULL );

            /* get row data */
            lhs = row->lhs - (SCIPsetIsInfinity(set, -row->lhs) ? 0.0 : row->constant);
            rhs = row->rhs - (SCIPsetIsInfinity(set,  row->rhs) ? 0.0 : row->constant);
            nnonz = row->nlpcols;
            assert( nnonz <= lp->ncols );
            rowcols = row->cols;
            rowvals = row->vals;

            sum = 0.0;
            for( j = 0; j < nnonz; ++j )
               sum += rowvals[j] * primal[rowcols[j]->lppos];
            sum /= alpha;

            /* if we have an equation */
            if( SCIPsetIsEQ(set, lhs, rhs) )
            {
               assert( SCIPsetIsFeasEQ(set, sum, lhs) );
            }
            else
            {
               /* treat lhs */
               if( !SCIPsetIsInfinity(set, REALABS(lhs)) )
               {
                  assert( SCIPsetIsFeasZero(set, primal[lp->ncols+1+cnt]) || SCIPsetIsFeasGT(set, sum, lhs) );
                  ++cnt;
               }
               /* treat rhs */
               if( !SCIPsetIsInfinity(set, REALABS(rhs)) )
               {
                  assert( SCIPsetIsFeasZero(set, primal[lp->ncols+1+cnt]) || SCIPsetIsFeasLT(set, sum, rhs) );
                  ++cnt;
               }
            }
         }
         if( inclobjcutoff )
         {
            SCIP_Real sum;
#ifndef NDEBUG
            SCIP_Real rhs;

            rhs = lp->cutoffbound - getFiniteLooseObjval(lp, set, prob);
#endif
            sum = 0.0;
            for( j = 0; j < lp->ncols; ++j )
               sum += lp->cols[j]->obj * primal[lp->cols[j]->lppos];
            sum /= alpha;

            assert( SCIPsetIsFeasZero(set, primal[lp->ncols+1+cnt]) || SCIPsetIsFeasLT(set, sum, rhs) );
            ++cnt;
         }
      }
      /* check bounds */
      for( j = 0; j < lp->ncols; ++j )
      {
         SCIP_COL* col;
#ifndef NDEBUG
         SCIP_Real val;
#endif

         col = lp->cols[j];
         assert( col != NULL );
#ifndef NDEBUG
         val = primal[col->lppos] / alpha;
#endif
         /* if the variable is not fixed */
         if( !SCIPsetIsEQ(set, col->lb, col->ub) )
         {
            /* treat lb */
            if( !SCIPsetIsInfinity(set, REALABS(col->lb)) )
            {
               assert( SCIPsetIsFeasZero(set, primal[lp->ncols+1+cnt]) || SCIPsetIsFeasGT(set, val, col->lb) );
               ++cnt;
            }
            /* treat rhs */
            if( !SCIPsetIsInfinity(set, REALABS(col->ub)) )
            {
               assert( SCIPsetIsFeasZero(set, primal[lp->ncols+1+cnt]) || SCIPsetIsFeasLT(set, val, col->ub) );
               ++cnt;
            }
         }
      }
#endif

      /* free */
      SCIPsetFreeBufferArray(set, &primal);

      *success = TRUE;
   }

   return SCIP_OKAY;
}

/** compute relative interior point
 *
 *  We use the approach of@par
 *  R. Freund, R. Roundy, M. J. Todd@par
 *  "Identifying the Set of Always-Active Constraints in a System of Linear Inequalities by a Single Linear Program"@par
 *  Tech. Rep, No. 1674-85, Sloan School, M.I.T., 1985
 *
 *  to compute a relative interior point for the current LP.
 *
 *  Assume the original LP looks as follows:
 *  \f[
 *     \begin{array}{rrl}
 *        \min & c^T x &\\
 *             & A x & \geq a\\
 *             & B x & \leq b\\
 *             & D x & = d.
 *     \end{array}
 *  \f]
 *  Note that bounds should be included in the system.
 *
 *  To find an interior point the following LP does the job:
 *  \f[
 *     \begin{array}{rrl}
 *        \max & 1^T y &\\
 *             & A x - y - \alpha a & \geq 0\\
 *             & B x + y - \alpha b & \leq 0\\
 *             & D x - \alpha d & = 0\\
 *             & 0 \leq y & \leq 1\\
 *             & \alpha & \geq 1.
 *     \end{array}
 *  \f]
 *  If the original LP is feasible, this LP is feasible as well. Any optimal solution yields the relative interior point
 *  \f$x^*_j/\alpha^*\f$. Note that this will just produce some relative interior point. It does not produce a
 *  particular relative interior point, e.g., one that maximizes the distance to the boundary in some norm.
 */
SCIP_RETCODE SCIPlpComputeRelIntPoint(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             relaxrows,          /**< should the rows be relaxed */
   SCIP_Bool             inclobjcutoff,      /**< should a row for the objective cutoff be included */
   SCIP_Real             timelimit,          /**< time limit for LP solver */
   int                   iterlimit,          /**< iteration limit for LP solver */
   SCIP_Real*            point,              /**< array to store relative interior point on exit */
   SCIP_Bool*            success             /**< buffer to indicate whether interior point was successfully computed */
   )
{
   SCIP_LPI* lpi;

   SCIP_RETCODE retcode;

   assert(set != NULL);
   assert(lp  != NULL);
   assert(point != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check time and iteration limits */
   if ( timelimit <= 0.0 || iterlimit <= 0 )
      return SCIP_OKAY;

   /* exit if there are no columns */
   assert(lp->nrows >= 0);
   assert(lp->ncols >= 0);
   if( lp->ncols == 0 )
      return SCIP_OKAY;

   /* disable objective cutoff if we have none */
   if( inclobjcutoff && (SCIPsetIsInfinity(set, lp->cutoffbound) || lp->looseobjvalinf > 0 || lp->looseobjval == SCIP_INVALID) ) /*lint !e777 */
      inclobjcutoff = FALSE;

   SCIPsetDebugMsg(set, "Computing relative interior point to current LP.\n");

   /* if there are no rows, we return the zero point */
   if( lp->nrows == 0 && !inclobjcutoff )
   {
      /* create zero point */
      BMSclearMemoryArray(point, lp->ncols);
      *success = TRUE;

      return SCIP_OKAY;
   }

   /* create auxiliary LP */
   SCIP_CALL( SCIPlpiCreate(&lpi, messagehdlr, "relativeInterior", SCIP_OBJSEN_MAXIMIZE) );

   /* catch return code and ensure that lpi is freed, anyway */
   retcode = computeRelIntPoint(lpi, set, messagehdlr, lp, prob, relaxrows, inclobjcutoff, timelimit, iterlimit, point, success);

   SCIP_CALL( SCIPlpiFree(&lpi) );

   return retcode;
}
