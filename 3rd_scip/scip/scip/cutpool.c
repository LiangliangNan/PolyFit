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

/**@file   cutpool.c
 * @brief  methods for storing cuts in a cut pool
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Gerald Gamrath
 * @author Marc Pfetsch
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/lp.h"
#include "scip/cons.h"
#include "scip/sepa.h"
#include "scip/sepastore.h"
#include "scip/cutpool.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_cutpool.h"



/*
 * Hash functions
 */

/** gets the hash key of a cut */
static
SCIP_DECL_HASHGETKEY(hashGetKeyCut)
{  /*lint --e{715}*/
   SCIP_CUT* cut;

   cut = (SCIP_CUT*)elem;
   assert(cut != NULL);
   assert(cut->row != NULL);

   /* the key of a cut is the row */
   return cut->row;
}

/** returns TRUE iff both cuts are identical */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqCut)
{  /*lint --e{715}*/
   SCIP_ROW* row1;
   SCIP_ROW* row2;
   SCIP_Real row1scale;
   SCIP_Real row2scale;
   SCIP_SET* set;

   row1 = (SCIP_ROW*)key1;
   row2 = (SCIP_ROW*)key2;
   assert(row1 != NULL);
   assert(row2 != NULL);

   /* return true if the row is the same */
   if( row1 == row2 )
      return TRUE;

   assert(row1->validminmaxidx);
   assert(row2->validminmaxidx);

   /* compare the trivial characteristics of the rows */
   if( row1->len != row2->len
      || row1->minidx != row2->minidx
      || row1->maxidx != row2->maxidx
      )
      return FALSE;

   set = (SCIP_SET*) userptr;

   /* set scale for the rows such that the largest absolute coefficient is 1.0 */
   row1scale = 1.0 / SCIProwGetMaxval(row1, set);
   row2scale = 1.0 / SCIProwGetMaxval(row2, set);

   /* check if scaled min value is feas equal first */
   if( !SCIPsetIsFeasEQ(set, row1scale * SCIProwGetMinval(row1, set),
                             row2scale * SCIProwGetMinval(row2, set)) )
      return FALSE;

   SCIProwSort(row1);
   assert(row1->lpcolssorted);
   assert(row1->nonlpcolssorted);

   SCIProwSort(row2);
   assert(row2->lpcolssorted);
   assert(row2->nonlpcolssorted);

   /* currently we are only handling rows which are completely linked or not linked at all */
   assert(row1->nunlinked == 0 || row1->nlpcols == 0);
   assert(row2->nunlinked == 0 || row2->nlpcols == 0);

   /* set scale sign such that the rows are of the form ax <= b */
   if( SCIPsetIsInfinity(set, row1->rhs) )
      row1scale = -row1scale;
   if( SCIPsetIsInfinity(set, row2->rhs) )
      row2scale = -row2scale;

   /* both rows have LP columns, or none of them has, or one has only LP colums and the other only non-LP columns,
    * so we can rely on the sorting of the columns
    */
   if( (row1->nlpcols == 0) == (row2->nlpcols == 0)
      || (row1->nlpcols == 0 && row2->nlpcols == row2->len)
      || (row1->nlpcols == row1->len && row2->nlpcols == 0) )
   {
      int i;

      if( (row1->nlpcols == 0) == (row2->nlpcols == 0) )
      {
#ifndef NDEBUG
         /* in debug mode, we check that we can rely on the partition into LP columns and non-LP columns */
         int i2;

         i = 0;
         i2 = row2->nlpcols;
         while( i < row1->nlpcols && i2 < row2->len )
         {
            assert(row1->cols[i] != row2->cols[i2]);
            if( row1->cols_index[i] < row2->cols_index[i2] )
               ++i;
            else
            {
               assert(row1->cols_index[i] > row2->cols_index[i2]);
               ++i2;
            }
         }
         assert(i == row1->nlpcols || i2 == row2->len);

         i = row1->nlpcols;
         i2 = 0;
         while( i < row1->len && i2 < row2->nlpcols )
         {
            assert(row1->cols[i] != row2->cols[i2]);
            if( row1->cols_index[i] < row2->cols_index[i2] )
               ++i;
            else
            {
               assert(row1->cols_index[i] > row2->cols_index[i2]);
               ++i2;
            }
         }
         assert(i == row1->len || i2 == row2->nlpcols);
#endif

         /* both rows are linked and the number of lpcolumns is not equal so they cannot be equal */
         if( row1->nlpcols != row2->nlpcols )
            return FALSE;
      }

      /* compare the columns of the rows */
      for( i = 0; i < row1->len; ++i )
      {
         if( row1->cols_index[i] != row2->cols_index[i] )
            return FALSE;
      }

      /* compare the coefficients of the rows */
      for( i = 0; i < row1->len; ++i )
      {
         if( !SCIPsetIsFeasEQ(set, row1scale * row1->vals[i], row2scale * row2->vals[i]) )
            return FALSE;
      }
   }
   /* one row has LP columns, but the other not, that could be because the one without was just created and isn't
    * linked yet; in this case, one column could be an LP column in one row and a non-LP column in the other row, so we
    * cannot rely on the partition; thus, we iteratively check whether the next column of row1 is either the next LP
    * column of row2 or the next non-LP column of row2 and the coefficients are equal
    */
   else
   {
      int i1;
      int ilp;
      int inlp;

      /* ensure that row1 is the row without LP columns, switch the rows, if neccessary */
      if( row2->nlpcols == 0 )
      {
         SCIP_ROW* tmprow;
         SCIP_Real tmpscale;

         tmprow = row2;
         row2 = row1;
         row1 = tmprow;

         tmpscale = row2scale;
         row2scale = row1scale;
         row1scale = tmpscale;
      }
      assert(row1->nlpcols == 0 && row2->nlpcols > 0);

      ilp = 0;
      inlp = row2->nlpcols;

      /* compare the columns and coefficients of the rows */
      for( i1 = 0; i1 < row1->len; ++i1 )
      {
         /* current column of row1 is the current LP column of row2, check the coefficient */
         if( ilp < row2->nlpcols && row1->cols[i1] == row2->cols[ilp] )
         {
            if( !SCIPsetIsFeasEQ(set, row1scale * row1->vals[i1], row2scale * row2->vals[ilp]) )
               return FALSE;
            else
               ++ilp;
         }
         /* current column of row1 is the current non-LP column of row2, check the coefficient */
         else if( inlp < row2->len && row1->cols[i1] == row2->cols[inlp] )
         {
            if( !SCIPsetIsFeasEQ(set, row1scale * row1->vals[i1], row2scale * row2->vals[inlp]) )
               return FALSE;
            else
               ++inlp;
         }
         /* current column of row1 is neither the current LP column of row2, nor the current non-LP column of row 2 */
         else
            return FALSE;
      }
   }

   return TRUE;
}

static
SCIP_DECL_HASHKEYVAL(hashKeyValCut)
{  /*lint --e{715}*/
   SCIP_ROW* row;
   int i;
   SCIP_Real scale;
   SCIP_SET* set;
   uint64_t hash;

   set = (SCIP_SET*) userptr;
   row = (SCIP_ROW*)key;
   assert(row != NULL);
   assert(row->len > 0);

   scale = 1.0 / SCIProwGetMaxval(row, set);
   if( SCIPsetIsInfinity(set, row->rhs) )
      scale = -scale;

   hash = (uint64_t) (long) row->len;

   for( i = 0; i < row->len; ++i )
   {
      SCIP_Real val = scale * row->vals[i];

      hash += SCIPhashTwo(SCIPrealHashCode(val), row->cols_index[i]);
   }

   return hash;
}


/*
 * dynamic memory arrays
 */

/** resizes cuts array to be able to store at least num entries */
static
SCIP_RETCODE cutpoolEnsureCutsMem(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(cutpool != NULL);
   assert(set != NULL);

   if( num > cutpool->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&cutpool->cuts, newsize) );
      cutpool->cutssize = newsize;
   }
   assert(num <= cutpool->cutssize);

   return SCIP_OKAY;
}



/*
 * Cut methods
 */

/** creates a cut and captures the row */
static
SCIP_RETCODE cutCreate(
   SCIP_CUT**            cut,                /**< pointer to store the cut */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row this cut represents */
   )
{
   assert(cut != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* allocate cut memory */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, cut) );
   (*cut)->row = row;
   (*cut)->age = 0;
   (*cut)->processedlp = -1;
   (*cut)->processedlpsol = -1;
   (*cut)->pos = -1;

   /* capture row */
   SCIProwCapture(row);

   return SCIP_OKAY;
}

/** frees a cut and releases the row */
static
SCIP_RETCODE cutFree(
   SCIP_CUT**            cut,                /**< pointer to store the cut */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(cut != NULL);
   assert(*cut != NULL);
   assert((*cut)->row != NULL);
   assert(blkmem != NULL);

   /* release row */
   SCIP_CALL( SCIProwRelease(&(*cut)->row, blkmem, set, lp) );

   /* free cut memory */
   BMSfreeBlockMemory(blkmem, cut);

   return SCIP_OKAY;
}

/** returns whether the cut's age exceeds the age limit */
static
SCIP_Bool cutIsAged(
   SCIP_CUT*             cut,                /**< cut to check */
   int                   agelimit            /**< maximum age a cut can reach before it is deleted from the pool, or -1 */
   )
{
   assert(cut != NULL);

   /* since agelimit can be -1 cast to unsigned before comparison, then it is the maximum unsigned value in that case */
   return (unsigned int)cut->age > (unsigned int)agelimit;
}

/** gets the row of the cut */
SCIP_ROW* SCIPcutGetRow(
   SCIP_CUT*             cut                 /**< cut */
   )
{
   assert(cut != NULL);

   return cut->row;
}

/** gets the age of the cut: the number of consecutive cut pool separation rounds where the cut was neither in the LP nor violated */
int SCIPcutGetAge(
   SCIP_CUT*             cut                 /**< cut */
   )
{
   assert(cut != NULL);

   return cut->age;
}

/** returns the ratio of LPs where the row belonging to this cut was active in an LP solution, i.e.
 *  where the age of its row has not been increased
 *
 *  @see SCIPcutGetAge() to get the age of a cut
 */
SCIP_Real SCIPcutGetLPActivityQuot(
   SCIP_CUT*             cut                 /**< cut */
   )
{
   SCIP_Longint nlpsaftercreation;
   SCIP_Longint activeinlpcounter;

   assert(cut != NULL);
   assert(cut->row != NULL);

   nlpsaftercreation = SCIProwGetNLPsAfterCreation(cut->row);
   activeinlpcounter = SCIProwGetActiveLPCount(cut->row);

   return (nlpsaftercreation > 0 ? activeinlpcounter / (SCIP_Real)nlpsaftercreation : 0.0);
}

/*
 * Cutpool methods
 */

/** creates cut pool */
SCIP_RETCODE SCIPcutpoolCreate(
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   agelimit,           /**< maximum age a cut can reach before it is deleted from the pool */
   SCIP_Bool             globalcutpool       /**< is this the global cut pool of SCIP? */
   )
{
   assert(cutpool != NULL);
   assert(agelimit >= -1);

   SCIP_ALLOC( BMSallocMemory(cutpool) );

   SCIP_CALL( SCIPclockCreate(&(*cutpool)->poolclock, SCIP_CLOCKTYPE_DEFAULT) );

   SCIP_CALL( SCIPhashtableCreate(&(*cutpool)->hashtable, blkmem, 
         (set->misc_usesmalltables ? SCIP_HASHSIZE_CUTPOOLS_SMALL : SCIP_HASHSIZE_CUTPOOLS),
         hashGetKeyCut, hashKeyEqCut, hashKeyValCut, (void*) set) );

   (*cutpool)->cuts = NULL;
   (*cutpool)->cutssize = 0;
   (*cutpool)->ncuts = 0;
   (*cutpool)->nremovablecuts = 0;
   (*cutpool)->agelimit = agelimit;
   (*cutpool)->processedlp = -1;
   (*cutpool)->processedlpsol = -1;
   (*cutpool)->processedlpefficacy = SCIP_INVALID;
   (*cutpool)->processedlpsolefficacy = SCIP_INVALID;
   (*cutpool)->firstunprocessed = 0;
   (*cutpool)->firstunprocessedsol = 0;
   (*cutpool)->maxncuts = 0;
   (*cutpool)->ncalls = 0;
   (*cutpool)->ncutsfound = 0;
   (*cutpool)->globalcutpool = globalcutpool;

   return SCIP_OKAY;
}

/** frees cut pool */
SCIP_RETCODE SCIPcutpoolFree(
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(cutpool != NULL);
   assert(*cutpool != NULL);

   /* remove all cuts from the pool */
   SCIP_CALL( SCIPcutpoolClear(*cutpool, blkmem, set, lp) );

   /* free clock */
   SCIPclockFree(&(*cutpool)->poolclock);

   /* free hash table */
   SCIPhashtableFree(&(*cutpool)->hashtable);

   BMSfreeMemoryArrayNull(&(*cutpool)->cuts);
   BMSfreeMemory(cutpool);

   return SCIP_OKAY;
}

/** removes all rows from the cut pool */
SCIP_RETCODE SCIPcutpoolClear(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   int i;

   assert(cutpool != NULL);

   /* free cuts */
   for( i = 0; i < cutpool->ncuts; ++i )
   {
      if( cutpool->globalcutpool )
         cutpool->cuts[i]->row->inglobalcutpool = FALSE;
      SCIProwUnlock(cutpool->cuts[i]->row);
      SCIP_CALL( cutFree(&cutpool->cuts[i], blkmem, set, lp) );
   }

   cutpool->ncuts = 0;
   cutpool->nremovablecuts = 0;

   return SCIP_OKAY;
}

/** removes the cut from the cut pool */
static
SCIP_RETCODE cutpoolDelCut(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CUT*             cut                 /**< cut to remove */
   )
{
   int pos;

   assert(cutpool != NULL);
   assert(cutpool->firstunprocessed <= cutpool->ncuts);
   assert(cutpool->firstunprocessedsol <= cutpool->ncuts);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(cutpool->processedlp <= stat->lpcount);
   assert(cutpool->processedlpsol <= stat->lpcount);
   assert(cut != NULL);
   assert(cut->row != NULL);

   pos = cut->pos;
   assert(0 <= pos && pos < cutpool->ncuts);
   assert(cutpool->cuts[pos] == cut);

   /* decrease the number of removable cuts counter (row might have changed its removable status -> counting might not
    * be correct
    */
   if( SCIProwIsRemovable(cut->row) && cutpool->nremovablecuts > 0 )
      cutpool->nremovablecuts--;

   /* if this is the global cut pool of SCIP, mark the row to not be member anymore */
   if( cutpool->globalcutpool )
   {
      assert(cut->row->inglobalcutpool);
      cut->row->inglobalcutpool = FALSE;
   }

   /* remove the cut from the hash table */
   assert(SCIPhashtableExists(cutpool->hashtable, (void*)cut));
   SCIP_CALL( SCIPhashtableRemove(cutpool->hashtable, (void*)cut) );
   assert(! SCIPhashtableExists(cutpool->hashtable, (void*)cut));

   /* unlock the row */
   SCIProwUnlock(cut->row);

   /* free the cut */
   SCIP_CALL( cutFree(&cutpool->cuts[pos], blkmem, set, lp) );

   /* move the last cut of the pool to the free position */
   if( pos < cutpool->ncuts-1 )
   {
      cutpool->cuts[pos] = cutpool->cuts[cutpool->ncuts-1];
      cutpool->cuts[pos]->pos = pos;
      assert(cutpool->cuts[pos]->processedlp <= stat->lpcount);
      assert(cutpool->cuts[pos]->processedlpsol <= stat->lpcount);
      if( cutpool->cuts[pos]->processedlp < stat->lpcount )
         cutpool->firstunprocessed = MIN(cutpool->firstunprocessed, pos);
      if( cutpool->cuts[pos]->processedlpsol < stat->lpcount )
         cutpool->firstunprocessedsol = MIN(cutpool->firstunprocessedsol, pos);
   }
   else
   {
      cutpool->firstunprocessed = MIN(cutpool->firstunprocessed, cutpool->ncuts-1);
      cutpool->firstunprocessedsol = MIN(cutpool->firstunprocessedsol, cutpool->ncuts-1);
   }

   cutpool->ncuts--;

   return SCIP_OKAY;
}

/** checks if cut is already existing */
SCIP_Bool SCIPcutpoolIsCutNew(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CUT* othercut;
   assert(cutpool != NULL);
   assert(row != NULL);

   if( row->len == 0 )
   {
      /* trivial cut is only new if it proves infeasibility */
      return SCIPsetIsFeasLT(set, row->constant, row->lhs) || SCIPsetIsFeasGT(set, row->constant, row->rhs);
   }

   othercut = (SCIP_CUT*)SCIPhashtableRetrieve(cutpool->hashtable, (void*)row);
   /* check in hash table, if cut already exists in the pool */
   if( othercut == NULL )
   {
      return TRUE;
   }
   else if( othercut->row != row )
   {
      SCIP_ROW* otherrow = othercut->row;
      SCIP_Real otherrhs;
      SCIP_Real rhs;
      SCIP_Real scale;
      SCIP_Real otherscale;

      /* since we are comparing the improvement with an absolute value, we apply a
       * scale to both rows such that the max absolute value is 1.0.
       * Then bring the cut into the form ax <= b
       */
      scale = 1.0 / SCIProwGetMaxval(row, set);
      otherscale = 1.0 / SCIProwGetMaxval(otherrow, set);

      if( SCIPsetIsInfinity(set, otherrow->rhs) )
      {
         otherrhs = otherscale * (otherrow->constant - otherrow->lhs);
      }
      else
      {
         otherrhs = otherscale * (otherrow->rhs - otherrow->constant);
      }

      if( SCIPsetIsInfinity(set, row->rhs) )
      {
         rhs = scale * (row->constant - row->lhs);
      }
      else
      {
         rhs = scale * (row->rhs - row->constant);
      }

      if( SCIPsetIsFeasLT(set, rhs, otherrhs) )
         return TRUE;
   }

   return FALSE;
}

/** if not already existing, adds row to cut pool and captures it */
SCIP_RETCODE SCIPcutpoolAddRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CUT* othercut;
   assert(cutpool != NULL);
   assert(row != NULL);

   if( row->len == 0 )
      return SCIP_OKAY;

   othercut = (SCIP_CUT*)SCIPhashtableRetrieve(cutpool->hashtable, (void*)row);
   /* check in hash table, if cut already exists in the pool */
   if( othercut == NULL )
   {
      SCIP_CALL( SCIPcutpoolAddNewRow(cutpool, blkmem, set, stat, lp, row) );
   }
   else
   {
      SCIP_ROW* otherrow = othercut->row;
      SCIP_Real otherrhs;
      SCIP_Real rhs;
      SCIP_Real scale;
      SCIP_Real otherscale;

      /* since we are comparing the improvement with an absolute value, we apply a
       * scale to both rows such that the max absolute value is 1.0.
       * Then bring the cut into the form ax <= b
       */
      scale = 1.0 / SCIProwGetMaxval(row, set);
      otherscale = 1.0 / SCIProwGetMaxval(otherrow, set);

      if( SCIPsetIsInfinity(set, otherrow->rhs) )
      {
         otherrhs = otherscale * (otherrow->constant - otherrow->lhs);
      }
      else
      {
         otherrhs = otherscale * (otherrow->rhs - otherrow->constant);
      }

      if( SCIPsetIsInfinity(set, row->rhs) )
      {
         rhs = scale * (row->constant - row->lhs);
      }
      else
      {
         rhs = scale * (row->rhs - row->constant);
      }

      if( SCIPsetIsFeasLT(set, rhs, otherrhs) )
      {
         SCIP_CALL( cutpoolDelCut(cutpool, blkmem, set, stat, lp, othercut) );

         /* use recursion, since in rare cases new cut might compare equal to multiple other cuts
          * that do not compare equal themselve due to non-transitivity of epsilon comparisons
          */
         SCIP_CALL( SCIPcutpoolAddRow(cutpool, blkmem, set, stat, lp, row) );
      }
   }

   return SCIP_OKAY;
}

/** adds row to cut pool and captures it; doesn't check for multiple cuts */
SCIP_RETCODE SCIPcutpoolAddNewRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_Real thisefficacy;
   SCIP_CUT* cut;

   assert(cutpool != NULL);
   assert(row != NULL);

   /* check, if row is modifiable or local */
   if( SCIProwIsModifiable(row) )
   {
      SCIPerrorMessage("cannot store modifiable row <%s> in a cut pool\n", SCIProwGetName(row));
      return SCIP_INVALIDDATA;
   }
   if( SCIProwIsLocal(row) )
   {
      SCIPerrorMessage("cannot store locally valid row <%s> in a cut pool\n", SCIProwGetName(row));
      return SCIP_INVALIDDATA;
   }

   assert(! row->inglobalcutpool);

   /* only called to ensure that minidx and maxidx are up-to-date */
   (void) SCIProwGetMaxidx(row, set);
   assert(row->validminmaxidx);

   /* create the cut */
   SCIP_CALL( cutCreate(&cut, blkmem, row) );
   cut->pos = cutpool->ncuts;

   /* add cut to the pool */
   SCIP_CALL( cutpoolEnsureCutsMem(cutpool, set, cutpool->ncuts+1) );
   cutpool->cuts[cutpool->ncuts] = cut;
   cutpool->ncuts++;
   cutpool->maxncuts = MAX(cutpool->maxncuts, cutpool->ncuts);
   if( SCIProwIsRemovable(row) )
      cutpool->nremovablecuts++;

   assert(!SCIPhashtableExists(cutpool->hashtable, (void*)cut));

   /* insert cut in the hash table */
   SCIP_CALL( SCIPhashtableInsert(cutpool->hashtable, (void*)cut) );

   assert(SCIPhashtableExists(cutpool->hashtable, (void*)cut));

   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      thisefficacy = SCIProwGetLPEfficacy(row, set, stat, lp);
      stat->bestefficacy = MAX(thisefficacy, stat->bestefficacy);
   }

   /* if this is the global cut pool of SCIP, mark the row to be member of the pool */
   if( cutpool->globalcutpool )
      row->inglobalcutpool = TRUE;

   /* lock the row */
   SCIProwLock(row);

   return SCIP_OKAY;
}

/** removes the LP row from the cut pool */
SCIP_RETCODE SCIPcutpoolDelRow(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row                 /**< row to remove */
   )
{
   SCIP_CUT* cut;

   assert(cutpool != NULL);
   assert(row != NULL);

   /* find the cut in hash table */
   cut = (SCIP_CUT*)SCIPhashtableRetrieve(cutpool->hashtable, (void*)row);
   if( cut == NULL )
   {
      SCIPerrorMessage("row <%s> is not existing in cutpool %p\n", SCIProwGetName(row), cutpool);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( cutpoolDelCut(cutpool, blkmem, set, stat, lp, cut) );

   return SCIP_OKAY;
}


/** separates cuts of the cut pool */
SCIP_RETCODE SCIPcutpoolSeparate(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL for LP-solution) */
   SCIP_Bool             cutpoolisdelayed,   /**< is the cutpool delayed (count cuts found)? */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CUT* cut;
   SCIP_Bool found;
   SCIP_Bool cutoff;
   SCIP_Real minefficacy;
   SCIP_Bool retest;
   int firstunproc;
   int oldncuts;
   int nefficaciouscuts;
   int c;

   assert(cutpool != NULL);
   assert(stat != NULL);
   assert(cutpool->processedlp <= stat->lpcount);
   assert(cutpool->processedlpsol <= stat->lpcount);
   assert(cutpool->firstunprocessed <= cutpool->ncuts);
   assert(cutpool->firstunprocessedsol <= cutpool->ncuts);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* don't separate cut pool in the root node, if there are no removable cuts */
   if( root && cutpool->nremovablecuts == 0 )
      return SCIP_OKAY;

   if ( sol == NULL )
   {
      if( cutpool->processedlp < stat->lpcount )
         cutpool->firstunprocessed = 0;
      if( cutpool->firstunprocessed == cutpool->ncuts )
         return SCIP_OKAY;
      firstunproc = cutpool->firstunprocessed;
   }
   else
   {
      if( cutpool->processedlpsol < stat->lpcount )
         cutpool->firstunprocessedsol = 0;
      if( cutpool->firstunprocessedsol == cutpool->ncuts )
         return SCIP_OKAY;
      firstunproc = cutpool->firstunprocessedsol;
   }

   *result = SCIP_DIDNOTFIND;
   cutpool->ncalls++;
   found = FALSE;
   minefficacy = stat->bestefficacy * stat->minefficacyfac;

   if( sol == NULL )
   {
      retest = cutpool->processedlpefficacy > minefficacy;
      cutpool->processedlpefficacy = minefficacy;
   }
   else
   {
      retest = cutpool->processedlpsolefficacy > minefficacy;
      cutpool->processedlpsolefficacy = minefficacy;
   }

   SCIPsetDebugMsg(set, "separating%s cut pool %p with %d cuts, beginning with cut %d\n", ( sol == NULL ) ? "" : " solution from", (void*)cutpool, cutpool->ncuts, firstunproc);

   /* start timing */
   SCIPclockStart(cutpool->poolclock, set);

   /* remember the current total number of found cuts */
   oldncuts = SCIPsepastoreGetNCuts(sepastore);
   nefficaciouscuts = 0;

   /* process all unprocessed cuts in the pool */
   cutoff = FALSE;
   for( c = firstunproc; c < cutpool->ncuts; ++c )
   {
      SCIP_Longint proclp;

      cut = cutpool->cuts[c];
      assert(cut != NULL);
      assert(cut->processedlp <= stat->lpcount);
      assert(cut->processedlpsol <= stat->lpcount);
      assert(cut->pos == c);

      proclp = ( sol == NULL ) ? cut->processedlp : cut->processedlpsol;

      if( retest || proclp < stat->lpcount )
      {
         SCIP_ROW* row;

         if ( sol == NULL )
            cut->processedlp = stat->lpcount;
         else
            cut->processedlpsol = stat->lpcount;

         row = cut->row;
         if( !SCIProwIsInLP(row) )
         {
            SCIP_Real efficacy;

            /* if the cut is a bound change (i.e. a row with only one variable), add it as bound change instead of LP
             * row; hence, we want to remove the bound change cut from the SCIP cut pool
             */
            if( !SCIProwIsModifiable(row) && SCIProwGetNNonz(row) == 1 )
            {
               /* insert bound change cut into separation store which will force that cut */
               SCIP_CALL( SCIPsepastoreAddCut(sepastore, blkmem, set, stat, eventqueue, eventfilter, lp, row, FALSE, root, &cutoff) );
               SCIP_CALL( cutpoolDelCut(cutpool, blkmem, set, stat, lp, cut) );

               if ( cutoff )
                  break;

               continue;
            }

            efficacy = sol == NULL ? SCIProwGetLPEfficacy(row, set, stat, lp) : SCIProwGetSolEfficacy(row, set, stat, sol);
            if( SCIPsetIsFeasPositive(set, efficacy) )
               ++nefficaciouscuts;

            if( efficacy >= minefficacy )
            {
               /* insert cut in separation storage */
               SCIPsetDebugMsg(set, " -> separated cut <%s> from the cut pool (feasibility: %g)\n",
                  SCIProwGetName(row), ( sol == NULL ) ? SCIProwGetLPFeasibility(row, set, stat, lp) : SCIProwGetSolFeasibility(row, set, stat, sol) );
               SCIP_CALL( SCIPsepastoreAddCut(sepastore, blkmem, set, stat, eventqueue, eventfilter, lp, row, FALSE, root, &cutoff) );

               /* count cuts */
               if ( cutpoolisdelayed )
               {
                  if ( SCIProwGetOriginSepa(row) != NULL )
                  {
                     SCIP_SEPA* sepa;

                     sepa = SCIProwGetOriginSepa(row);
                     SCIPsepaIncNCutsFound(sepa);
                     SCIPsepaIncNCutsFoundAtNode(sepa);
                  }
                  else if ( SCIProwGetOriginCons(row) != NULL )
                  {
                     SCIP_CONSHDLR* conshdlr;

                     conshdlr = SCIProwGetOriginCons(row);
                     SCIPconshdlrIncNCutsFound(conshdlr);
                  }
               }

               found = TRUE;
               cut->age = 0;

               if ( cutoff )
                  break;
            }
            else
            {
               cut->age++;
               if( cutIsAged(cut, cutpool->agelimit) )
               {
                  SCIP_CALL( cutpoolDelCut(cutpool, blkmem, set, stat, lp, cut) );
               }
            }
         }
      }
   }

   if ( sol == NULL )
   {
      cutpool->processedlp = stat->lpcount;
      cutpool->firstunprocessed = cutpool->ncuts;
   }
   else
   {
      cutpool->processedlpsol = stat->lpcount;
      cutpool->firstunprocessedsol = cutpool->ncuts;
   }

   if( nefficaciouscuts > 0 )
   {
      int maxncuts = SCIPsetGetSepaMaxcuts(set, root);
      int ncuts = SCIPsepastoreGetNCuts(sepastore) - oldncuts;

      maxncuts = MIN(maxncuts, nefficaciouscuts);

      /* update the number of found cuts */
      cutpool->ncutsfound += ncuts;

      if( ncuts > (0.5 * maxncuts) )
      {
         stat->ncutpoolfails = MIN(stat->ncutpoolfails - 1, -1);
      }
      else if( ncuts == 0 || (ncuts < (0.05 * maxncuts)) )
      {
         stat->ncutpoolfails = MAX(stat->ncutpoolfails + 1, 1);
      }
   }

   if( stat->ncutpoolfails == (root ? 2 : 10) )
   {
      cutpool->firstunprocessed = 0;
      cutpool->firstunprocessedsol = 0;
      stat->minefficacyfac *= 0.5;
      stat->ncutpoolfails = 0;
   }
   else if( stat->ncutpoolfails == -2 )
   {
      stat->minefficacyfac *= 1.2;
      stat->ncutpoolfails = 0;
   }

   /* stop timing */
   SCIPclockStop(cutpool->poolclock, set);

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if( found )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/** gets array of cuts in the cut pool */
SCIP_CUT** SCIPcutpoolGetCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->cuts;
}

/** gets number of cuts in the cut pool */
int SCIPcutpoolGetNCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->ncuts;
}

/** gets maximum number of cuts that were stored in the cut pool at the same time */
int SCIPcutpoolGetMaxNCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->maxncuts;
}

/** gets time in seconds used for separating cuts from the pool */
SCIP_Real SCIPcutpoolGetTime(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return SCIPclockGetTime(cutpool->poolclock);
}

/** get number of times, the cut pool was separated */
SCIP_Longint SCIPcutpoolGetNCalls(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->ncalls;
}

/** get total number of cuts that were separated from the cut pool */
SCIP_Longint SCIPcutpoolGetNCutsFound(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   )
{
   assert(cutpool != NULL);

   return cutpool->ncutsfound;
}

