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

/**@file    nlpioracle.c
 * @brief   implementation of NLPI oracle interface
 * @author  Stefan Vigerske
 *
 * @todo jacobi evaluation should be sparse
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "nlpi/nlpioracle.h"
#include "nlpi/pub_expr.h"
#include "nlpi/exprinterpret.h"

#include <string.h> /* for strlen */

/**@name NLPI Oracle data structures */
/**@{ */

struct SCIP_NlpiOracleCons
{
   SCIP_Real             lhs;                /**< left hand side (for constraint) or constant (for objective) */
   SCIP_Real             rhs;                /**< right hand side (for constraint) or constant (for objective) */

   int                   linsize;            /**< length of linidxs and lincoefs arrays */
   int                   nlinidxs;           /**< number of linear variable indices and coefficients */
   int*                  linidxs;            /**< variable indices in linear part, or NULL if none */
   SCIP_Real*            lincoefs;           /**< variable coefficients in linear part, of NULL if none */

   int                   quadsize;           /**< length of quadelems array */
   int                   nquadelems;         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems;          /**< quadratic elements, or NULL if none */

   int*                  exprvaridxs;        /**< indices of variables in expression tree, or NULL if no exprtree */
   SCIP_EXPRTREE*        exprtree;           /**< expression tree for nonlinear part, or NULL if none */

   char*                 name;               /**< name of constraint */
};
typedef struct SCIP_NlpiOracleCons SCIP_NLPIORACLECONS;

struct SCIP_NlpiOracle
{
   BMS_BLKMEM*           blkmem;             /**< block memory */
   SCIP_Real             infinity;           /**< value for infinity */
   char*                 name;               /**< name of problem */

   int                   varssize;           /**< length of variables related arrays */
   int                   nvars;              /**< number of variables */
   SCIP_Real*            varlbs;             /**< array with variable lower bounds */
   SCIP_Real*            varubs;             /**< array with variable upper bounds */
   char**                varnames;           /**< array with variable names */
   int*                  vardegrees;         /**< array with maximal degree of variable over objective and all constraints */
   SCIP_Bool             vardegreesuptodate; /**< whether the variable degrees are up to date */

   int                   consssize;          /**< length of constraints related arrays */
   int                   nconss;             /**< number of constraints */
   SCIP_NLPIORACLECONS** conss;              /**< constraints, or NULL if none */

   SCIP_NLPIORACLECONS*  objective;          /**< objective */

   int*                  jacoffsets;         /**< rowwise jacobi sparsity pattern: constraint offsets in jaccols */
   int*                  jaccols;            /**< rowwise jacobi sparsity pattern: indices of variables appearing in constraints */

   int*                  heslagoffsets;      /**< rowwise sparsity pattern of hessian matrix of Lagrangian: row offsets in heslagcol */
   int*                  heslagcols;         /**< rowwise sparsity pattern of hessian matrix of Lagrangian: column indices; sorted for each row */


   SCIP_EXPRINT*         exprinterpreter;    /**< interpreter for expression trees: evaluation and derivatives */
};

/**@} */

/**@name Local functions */
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

/** ensures that variables related arrays in oracle have at least a given length */
static
SCIP_RETCODE ensureVarsSize(
   SCIP_NLPIORACLE*      oracle,             /**< NLPIORACLE data structure */
   int                   minsize             /**< minimal required size */
   )
{
   assert(oracle != NULL);

   if( minsize > oracle->varssize )
   {
      int newsize;

      newsize = calcGrowSize(minsize);
      assert(newsize >= minsize);

      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varlbs, oracle->varssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varubs, oracle->varssize, newsize) );
      if( oracle->varnames != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->varnames, oracle->varssize, newsize) );
      }
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->vardegrees, oracle->varssize, newsize) );

      oracle->varssize = newsize;
   }
   assert(oracle->varssize >= minsize);

   return SCIP_OKAY;
}

/** ensures that constraints array in oracle has at least a given length */
static
SCIP_RETCODE ensureConssSize(
   SCIP_NLPIORACLE*      oracle,             /**< NLPIORACLE data structure */
   int                   minsize             /**< minimal required size */
   )
{
   assert(oracle != NULL);

   if( minsize > oracle->consssize )
   {
      int newsize;

      newsize = calcGrowSize(minsize);
      assert(newsize >= minsize);

      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->conss, oracle->consssize, newsize) );
      oracle->consssize = newsize;
   }
   assert(oracle->consssize >= minsize);

   return SCIP_OKAY;
}

/** ensures that arrays for linear part in a oracle constraints have at least a given length */
static
SCIP_RETCODE ensureConsLinSize(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPIORACLECONS*  cons,               /**< oracle constraint */
   int                   minsize             /**< minimal required size */
   )
{
   assert(blkmem != NULL);
   assert(cons != NULL);

   if( minsize > cons->linsize )
   {
      int newsize;

      newsize = calcGrowSize(minsize);
      assert(newsize >= minsize);

      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &cons->linidxs,  cons->linsize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &cons->lincoefs, cons->linsize, newsize) );
      cons->linsize = newsize;
   }
   assert(cons->linsize >= minsize);

   return SCIP_OKAY;
}

/** ensures that arrays for quadratic part in a oracle constraints have at least a given length */
static
SCIP_RETCODE ensureConsQuadSize(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPIORACLECONS*  cons,               /**< oracle constraint */
   int                   minsize             /**< minimal required size */
   )
{
   assert(blkmem != NULL);
   assert(cons != NULL);

   if( minsize > cons->quadsize )
   {
      int newsize;

      newsize = calcGrowSize(minsize);
      assert(newsize >= minsize);

      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &cons->quadelems, cons->quadsize, newsize) );
      cons->quadsize = newsize;
   }
   assert(cons->quadsize >= minsize);

   return SCIP_OKAY;
}

/** ensures that a given array of integers has at least a given length */
static
SCIP_RETCODE ensureIntArraySize(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int**                 intarray,           /**< array of integers */
   int*                  len,                /**< length of array (modified if reallocated) */
   int                   minsize             /**< minimal required array length */
   )
{
   assert(blkmem != NULL);
   assert(intarray != NULL);
   assert(len != NULL);

   if( minsize > *len )
   {
      int newsize;

      newsize = calcGrowSize(minsize);
      assert(newsize >= minsize);

      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, intarray, *len, newsize) );
      *len = newsize;
   }
   assert(*len >= minsize);

   return SCIP_OKAY;
}

/** Invalidates the sparsity pattern of the Jacobian.
 *  Should be called when constraints are added or deleted.
 */
static
void invalidateJacobiSparsity(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p invalidate jacobian sparsity\n", (void*)oracle);

   if( oracle->jacoffsets == NULL )
   { /* nothing to do */
      assert(oracle->jaccols == NULL);
      return;
   }

   assert(oracle->jaccols != NULL);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->jaccols,    oracle->jacoffsets[oracle->nconss]);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->jacoffsets, oracle->nconss + 1);
}

/** Invalidates the sparsity pattern of the Hessian of the Lagragian.
 *  Should be called when the objective is set or constraints are added or deleted.
 */
static
void invalidateHessianLagSparsity(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p invalidate hessian lag sparsity\n", (void*)oracle);

   if( oracle->heslagoffsets == NULL )
   { /* nothing to do */
      assert(oracle->heslagcols == NULL);
      return;
   }

   assert(oracle->heslagcols != NULL);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->heslagcols,    oracle->heslagoffsets[oracle->nvars]);
   BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->heslagoffsets, oracle->nvars + 1);
}

/** sorts a linear term, merges duplicate entries and removes entries with coefficient 0.0 */
static
void sortLinearCoefficients(
   int*                  nidxs,              /**< number of variables */
   int*                  idxs,               /**< indices of variables */
   SCIP_Real*            coefs               /**< coefficients of variables */
   )
{
   int offset;
   int j;

   assert(nidxs != NULL);
   assert(idxs  != NULL || *nidxs == 0);
   assert(coefs != NULL || *nidxs == 0);

   if( *nidxs == 0 )
      return;

   SCIPsortIntReal(idxs, coefs, *nidxs);

   offset = 0;
   j = 0;
   while( j+offset < *nidxs )
   {
      assert(idxs[j] >= 0);  /*lint !e613*/

      /* move j+offset to j, if different */
      if( offset > 0 )
      {
         idxs[j]  = idxs[j+offset];   /*lint !e613*/
         coefs[j] = coefs[j+offset];  /*lint !e613*/
      }

      /* add up coefs for j+offset+1... as long as they have the same index */
      while( j+offset+1 < *nidxs && idxs[j] == idxs[j+offset+1] )  /*lint !e613*/
      {
         coefs[j] += coefs[j+offset+1];  /*lint !e613*/
         ++offset;
      }

      /* if j'th element is 0, increase offset, otherwise increase j */
      if( coefs[j] == 0.0 )  /*lint !e613*/
         ++offset;
      else
         ++j;
   }
   *nidxs -= offset;
}

/** creates a NLPI constraint from given constraint data */
static
SCIP_RETCODE createConstraint(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPIORACLECONS** cons,               /**< buffer where to store pointer to constraint */
   int                   nlinidxs,           /**< length of linear part */
   const int*            linidxs,            /**< indices of linear part, or NULL if nlinidxs == 0 */
   const SCIP_Real*      lincoefs,           /**< coefficients of linear part, or NULL if nlinidxs == 0 */
   int                   nquadelems,         /**< lenght of quadratic part */
   const SCIP_QUADELEM*  quadelems,          /**< quadratic elements, or NULL if nquadelems == 0 */
   const int*            exprvaridxs,        /**< indicies of variables in expression tree, or NULL if exprtree == NULL */
   const SCIP_EXPRTREE*  exprtree,           /**< expression tree, or NULL */
   SCIP_Real             lhs,                /**< left-hand-side of constraint */
   SCIP_Real             rhs,                /**< right-hand-side of constraint */
   const char*           name                /**< name of constraint, or NULL */
   )
{
   assert(blkmem != NULL);
   assert(cons != NULL);
   assert(nlinidxs >= 0);
   assert(linidxs != NULL  || nlinidxs == 0);
   assert(lincoefs != NULL || nlinidxs == 0);
   assert(nquadelems >= 0);
   assert(quadelems != NULL || nquadelems == 0);
   assert(exprvaridxs != NULL || exprtree == NULL);
   assert(EPSLE(lhs, rhs, SCIP_DEFAULT_EPSILON));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, cons) );
   assert(*cons != NULL);
   BMSclearMemory(*cons);

   if( nlinidxs > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*cons)->linidxs,  linidxs,  nlinidxs) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*cons)->lincoefs, lincoefs, nlinidxs) );
      (*cons)->linsize  = nlinidxs;
      (*cons)->nlinidxs = nlinidxs;

      /* sort, merge duplicates, remove zero's */
      sortLinearCoefficients(&(*cons)->nlinidxs, (*cons)->linidxs, (*cons)->lincoefs);
      assert((*cons)->linidxs[0] >= 0);
   }

   if( nquadelems > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*cons)->quadelems, quadelems, nquadelems) );
      (*cons)->nquadelems = nquadelems;
      (*cons)->quadsize   = nquadelems;

      /* sort and squeeze quadratic part */
      SCIPquadelemSort((*cons)->quadelems, nquadelems);
      SCIPquadelemSqueeze((*cons)->quadelems, nquadelems, &(*cons)->nquadelems);
      assert((*cons)->nquadelems == 0 || (*cons)->quadelems[0].idx1 >= 0);
      assert((*cons)->nquadelems == 0 || (*cons)->quadelems[0].idx2 >= 0);
   }

   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeCopy(blkmem, &(*cons)->exprtree, (SCIP_EXPRTREE*)exprtree) );

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*cons)->exprvaridxs, exprvaridxs, SCIPexprtreeGetNVars((SCIP_EXPRTREE*)exprtree)) );
   }

   if( lhs > rhs )
   {
      assert(EPSEQ(lhs, rhs, SCIP_DEFAULT_EPSILON));
      lhs = rhs;
   }
   (*cons)->lhs = lhs;
   (*cons)->rhs = rhs;

   if( name != NULL )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*cons)->name, name, strlen(name)+1) );
   }

   return SCIP_OKAY;
}

/** frees a constraint */
static
void freeConstraint(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPIORACLECONS** cons                /**< pointer to constraint that should be freed */
   )
{
   assert(blkmem != NULL);
   assert(cons   != NULL);
   assert(*cons  != NULL);

   SCIPdebugMessage("free constraint %p\n", (void*)*cons);

   BMSfreeBlockMemoryArrayNull(blkmem, &(*cons)->linidxs, (*cons)->linsize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*cons)->lincoefs, (*cons)->linsize);

   BMSfreeBlockMemoryArrayNull(blkmem, &(*cons)->quadelems, (*cons)->quadsize);

   if( (*cons)->exprtree != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*cons)->exprvaridxs, SCIPexprtreeGetNVars((*cons)->exprtree));
      SCIP_CALL_ABORT( SCIPexprtreeFree(&(*cons)->exprtree) );
   }

   if( (*cons)->name != NULL )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &(*cons)->name, strlen((*cons)->name)+1);
   }

   BMSfreeBlockMemory(blkmem, cons);
   assert(*cons == NULL);
}

/** frees all constraints */
static
void freeConstraints(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   int i;

   assert(oracle != NULL);

   SCIPdebugMessage("%p free constraints\n", (void*)oracle);

   for( i = 0; i < oracle->nconss; ++i )
   {
      freeConstraint(oracle->blkmem, &oracle->conss[i]);
      assert(oracle->conss[i] == NULL);
   }
   oracle->nconss = 0;

   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->conss, oracle->consssize);
   oracle->consssize = 0;
}

/** moves one variable
 * The place where it moves to need to be empty (all NULL) but allocated.
 * Note that this function does not update the variable indices in the constraints!
 */
static
SCIP_RETCODE moveVariable(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to store NLPIORACLE data structure */
   int                   fromidx,            /**< index of variable to move */
   int                   toidx               /**< index of place where to move variable to */
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p move variable\n", (void*)oracle);

   assert(0 <= fromidx);
   assert(0 <= toidx);
   assert(fromidx < oracle->nvars);
   assert(toidx   < oracle->nvars);

   assert(oracle->varlbs[toidx] <= -oracle->infinity);
   assert(oracle->varubs[toidx] >=  oracle->infinity);
   assert(oracle->varnames == NULL || oracle->varnames[toidx] == NULL);
   assert(!oracle->vardegreesuptodate || oracle->vardegrees[toidx] == -1);

   oracle->varlbs[toidx] = oracle->varlbs[fromidx];
   oracle->varubs[toidx] = oracle->varubs[fromidx];

   oracle->varlbs[fromidx] = -oracle->infinity;
   oracle->varubs[fromidx] =  oracle->infinity;

   oracle->vardegrees[toidx]   = oracle->vardegrees[fromidx];
   oracle->vardegrees[fromidx] = -1;

   if( oracle->varnames != NULL )
   {
      oracle->varnames[toidx]   = oracle->varnames[fromidx];
      oracle->varnames[fromidx] = NULL;
   }

   return SCIP_OKAY;
}

/** frees all variables */
static
void freeVariables(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   int i;

   assert(oracle != NULL);

   SCIPdebugMessage("%p free variables\n", (void*)oracle);

   if( oracle->varnames != NULL )
   {
      for( i = 0; i < oracle->nvars; ++i )
      {
         if( oracle->varnames[i] != NULL )
         {
            BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->varnames[i], strlen(oracle->varnames[i])+1);  /*lint !e866*/
         }
      }
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->varnames, oracle->varssize);
   }
   oracle->nvars = 0;
   oracle->vardegreesuptodate = TRUE; 

   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->varlbs,     oracle->varssize);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->varubs,     oracle->varssize);
   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &oracle->vardegrees, oracle->varssize);

   oracle->varssize = 0;
}

/** increases variable degrees in oracle w.r.t. variables occuring in a single constraint */
static
void updateVariableDegreesCons(
   SCIP_NLPIORACLE*      oracle,             /**< oracle data structure */
   SCIP_NLPIORACLECONS*  cons                /**< oracle constraint */
   )
{
   int j;

   assert(oracle != NULL);
   assert(oracle->vardegrees != NULL);
   assert(cons != NULL);

   for( j = 0; j < cons->nlinidxs; ++j )
      if( oracle->vardegrees[cons->linidxs[j]] < 1 )
         oracle->vardegrees[cons->linidxs[j]] = 1;

   for( j = 0; j < cons->nquadelems; ++j )
   {
      if( oracle->vardegrees[cons->quadelems[j].idx1] < 2 )
         oracle->vardegrees[cons->quadelems[j].idx1] = 2;

      if( oracle->vardegrees[cons->quadelems[j].idx2] < 2 )
         oracle->vardegrees[cons->quadelems[j].idx2] = 2;
   }

   /* we could use exprtreeGetDegree to get actual degree of a variable in tree,
    * but so far no solver could make use of this information */
   if( cons->exprtree != NULL )
      for( j = SCIPexprtreeGetNVars(cons->exprtree)-1; j >= 0; --j )
         oracle->vardegrees[cons->exprvaridxs[j]] = INT_MAX;
}

/** Updates the degrees of all variables. */
static
void updateVariableDegrees(
   SCIP_NLPIORACLE*      oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   int c;

   assert(oracle != NULL);
   assert(oracle->nvars == 0 || oracle->vardegrees != NULL);
   assert(oracle->objective != NULL);

   SCIPdebugMessage("%p update variable degrees\n", (void*)oracle);

   if( oracle->vardegreesuptodate || oracle->nvars == 0 )
      return;

   /* assume all variables do not appear in NLP */
   BMSclearMemoryArray(oracle->vardegrees, oracle->nvars);

   updateVariableDegreesCons(oracle, oracle->objective);
   for( c = 0; c < oracle->nconss; ++c )
      updateVariableDegreesCons(oracle, oracle->conss[c]);

   oracle->vardegreesuptodate = TRUE;
}

/** applies a mapping of indices to one array of indices */
static
void mapIndices(
   int*                  indexmap,           /**< mapping from old variable indices to new indices */
   int                   nindices,           /**< number of indices in indices1 and indices2 */
   int*                  indices             /**< array of indices to adjust */
   )
{
   assert(indexmap != NULL);
   assert(nindices == 0 || indices != NULL);

   for( ; nindices ; --nindices, ++indices )
   {
      assert(indexmap[*indices] >= 0);
      *indices = indexmap[*indices];
   }
}

/** removes entries with index -1 (marked as deleted) from array of linear elements
 * assumes that array is sorted by index, i.e., all -1 are at the beginning
 */
static
void clearDeletedLinearElements(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int**                 linidxs,            /**< variable indices */
   SCIP_Real**           coefs,              /**< variable coefficients */
   int*                  nidxs               /**< number of indices */
   )
{
   int i;
   int offset;

   SCIPdebugMessage("clear deleted linear elements\n");

   assert(blkmem   != NULL);
   assert(linidxs  != NULL);
   assert(*linidxs != NULL);
   assert(coefs    != NULL);
   assert(*coefs   != NULL);
   assert(nidxs    != NULL);
   assert(*nidxs   > 0);

   /* search for beginning of non-delete entries @todo binary search? */
   for( offset = 0; offset < *nidxs; ++offset )
      if( (*linidxs)[offset] >= 0 )
         break;

   /* nothing was deleted */
   if( offset == 0 )
      return;

   /* some or all elements were deleted -> move remaining ones front */
   for( i = 0; i < *nidxs - offset; ++i )
   {
      (*linidxs)[i] = (*linidxs)[i+offset];
      (*coefs)[i]   = (*coefs)  [i+offset];
   }
   *nidxs -= offset;
}

/** removes entries with index pair (-1,-1) (marked as deleted) from array of quadratic elements
 * assumes that array is sorted, i.e., all deleted elements are at the beginning
 */
static
void clearDeletedQuadElements(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_QUADELEM**       quadelems,          /**< quadratic elements */
   int*                  nquadelems          /**< number of quadratic elements */
   )
{
   int i;
   int offset;

   SCIPdebugMessage("clear deleted quad elements\n");

   assert(blkmem      != NULL);
   assert(quadelems   != NULL);
   assert(*quadelems  != NULL);
   assert(nquadelems  != NULL);
   assert(*nquadelems > 0);

   /* search for beginning of non-delete entries @todo binary search? */
   for( offset = 0; offset < *nquadelems; ++offset )
   {
      /* either both variables are marked as deleted or none of them */
      assert(((*quadelems)[offset].idx1 >= 0) == ((*quadelems)[offset].idx2 >= 0));
      if( (*quadelems)[offset].idx1 >= 0 )
         break;
   }

   /* nothing was deleted */
   if( offset == 0 )
      return;

   /* some or all elements were deleted -> move remaining ones front */
   for( i = 0; i < *nquadelems - offset; ++i )
      (*quadelems)[i] = (*quadelems)[i+offset];
   *nquadelems -= offset;
}

/** applies a mapping of indices to an array of quadratic elements */
static
void mapIndicesQuad(
   int*                  indexmap,           /**< mapping from old variable indices to new indices */
   int                   nelems,             /**< number of quadratic elements */
   SCIP_QUADELEM*        elems               /**< array of quadratic elements to adjust */
   )
{
   assert(indexmap != NULL);
   assert(nelems == 0 || elems != NULL);

   for( ; nelems ; --nelems, ++elems )
   {
      assert(indexmap[elems->idx1] >= 0);
      assert(indexmap[elems->idx2] >= 0);
      elems->idx1 = indexmap[elems->idx1];
      elems->idx2 = indexmap[elems->idx2];
      /* swap indices if not idx1 <= idx2 */
      if( elems->idx1 > elems->idx2 )
      {
         int tmp = elems->idx2;
         elems->idx2 = elems->idx1;
         elems->idx1 = tmp;
      }
   }
}

/** computes the value of a function */
static
SCIP_RETCODE evalFunctionValue(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_NLPIORACLECONS*  cons,               /**< oracle constraint */
   const SCIP_Real*      x,                  /**< the point where to evaluate */
   SCIP_Real*            val                 /**< pointer to store function value */
   )
{  /*lint --e{715}*/
   assert(oracle != NULL);
   assert(cons != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(val != NULL);

   SCIPdebugMessage("%p eval function value\n", (void*)oracle);

   *val = 0.0;

   if( cons->nlinidxs > 0 )
   {
      int*       linidxs;
      SCIP_Real* lincoefs;
      int        nlin;

      nlin     = cons->nlinidxs;
      linidxs  = cons->linidxs;
      lincoefs = cons->lincoefs;
      assert(linidxs  != NULL);
      assert(lincoefs != NULL);
      assert(x != NULL);

      for( ; nlin > 0; --nlin, ++linidxs, ++lincoefs )
         *val += *lincoefs * x[*linidxs];
   }

   if( cons->nquadelems > 0 )
   {
      SCIP_QUADELEM* quadelems;
      int nquadelems;

      quadelems  = cons->quadelems;
      nquadelems = cons->nquadelems;
      assert(quadelems != NULL);
      assert(x != NULL);

      for( ; nquadelems > 0; --nquadelems, ++quadelems )
         *val += quadelems->coef * x[quadelems->idx1] * x[quadelems->idx2];
   }

   if( cons->exprtree != NULL )
   {
      SCIP_Real* xx;
      int        i;
      SCIP_Real  nlval;
      int        nvars;

      nvars = SCIPexprtreeGetNVars(cons->exprtree);

      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &xx, nvars) );
      for( i = 0; i < nvars; ++i )
      {
         assert(cons->exprvaridxs[i] >= 0);
         assert(cons->exprvaridxs[i] < oracle->nvars);
         xx[i] = x[cons->exprvaridxs[i]];  /*lint !e613 !e644*/
      }

      SCIP_CALL( SCIPexprintEval(oracle->exprinterpreter, cons->exprtree, xx, &nlval) );
      if( nlval != nlval || ABS(nlval) >= oracle->infinity )  /*lint !e777*/
         *val  = nlval;
      else
         *val += nlval;

      BMSfreeBlockMemoryArray(oracle->blkmem, &xx, nvars);
   }

   return SCIP_OKAY;
}

/** computes the value and gradient of a function
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
static
SCIP_RETCODE evalFunctionGradient(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_NLPIORACLECONS*  cons,               /**< oracle constraint */
   const SCIP_Real*      x,                  /**< the point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real* RESTRICT   val,                /**< pointer to store function value */
   SCIP_Real* RESTRICT   grad                /**< pointer to store function gradient */
   )
{  /*lint --e{715}*/
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(val != NULL);
   assert(grad != NULL);

   SCIPdebugMessage("%p eval function gradient\n", (void*)oracle);

   *val = 0.0;
   BMSclearMemoryArray(grad, oracle->nvars);

   if( cons->nlinidxs > 0 )
   {
      int*       linidxs;
      SCIP_Real* lincoefs;
      int        nlin;

      nlin     = cons->nlinidxs;
      linidxs  = cons->linidxs;
      lincoefs = cons->lincoefs;
      assert(linidxs  != NULL);
      assert(lincoefs != NULL);
      assert(x != NULL);

      for( ; nlin > 0; --nlin, ++linidxs, ++lincoefs )
      {
         *val += *lincoefs * x[*linidxs];
         assert(grad[*linidxs] == 0.0);   /* we do not like duplicate indices */
         grad[*linidxs] = *lincoefs;
      }
   }

   if( cons->nquadelems > 0 )
   {
      SCIP_Real tmp;
      SCIP_QUADELEM* quadelems;
      int nquadelems;

      quadelems  = cons->quadelems;
      nquadelems = cons->nquadelems;
      assert(quadelems != NULL);
      assert(x != NULL);

      for( ; nquadelems > 0; --nquadelems, ++quadelems )
      {
         tmp = quadelems->coef * x[quadelems->idx1];
         *val += tmp * x[quadelems->idx2];
         grad[quadelems->idx2] += tmp;
         grad[quadelems->idx1] += quadelems->coef * x[quadelems->idx2];
      }
   }

   if( cons->exprtree != NULL )
   {
      SCIP_Real* xx;
      SCIP_Real* g;
      int        i;
      SCIP_Real  nlval;
      int        nvars;

      xx = NULL;
      nvars = SCIPexprtreeGetNVars(cons->exprtree);

      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &g, nvars) );

      if( isnewx )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &xx, nvars) );
         for( i = 0; i < nvars; ++i )
         {
            assert(cons->exprvaridxs[i] >= 0);
            assert(cons->exprvaridxs[i] < oracle->nvars);
            xx[i] = x[cons->exprvaridxs[i]];  /*lint !e613*/
         }
      }

      SCIPdebugMessage("eval gradient of ");
      SCIPdebug( if( isnewx ) {printf("\nx ="); for( i = 0; i < nvars; ++i) printf(" %g", xx[i]); /*lint !e613*/ printf("\n");} )

      SCIP_CALL( SCIPexprintGrad(oracle->exprinterpreter, cons->exprtree, xx, isnewx, &nlval, g) );  /*lint !e644*/

      SCIPdebug( printf("g ="); for( i = 0; i < nvars; ++i) printf(" %g", g[i]); printf("\n"); )

      if( nlval != nlval || ABS(nlval) >= oracle->infinity )  /*lint !e777*/
      {
         BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
         BMSfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
         SCIPdebugMessage("gradient evaluation yield invalid function value %g\n", nlval);
         return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
      }
      else
      {
         *val += nlval;
         for( i = 0; i < nvars; ++i )
            if( !SCIPisFinite(g[i]) )  /*lint !e777*/
            {
               SCIPdebugMessage("gradient evaluation yield invalid gradient value %g\n", g[i]);
               BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
               BMSfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
               return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
            }
            else
            {
               grad[cons->exprvaridxs[i]] += g[i];
            }
      }

      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
      BMSfreeBlockMemoryArray(oracle->blkmem, &g, nvars);
   }

   return SCIP_OKAY;
}

/** collects nonzeros entries in colnz and increases the nzcount given indices of quadratic terms */
static
SCIP_RETCODE hessLagSparsitySetNzFlagForQuad(
   SCIP_NLPIORACLE*      oracle,             /**< NLPI oracle */
   int**                 colnz,              /**< indices of nonzero entries for each column */
   int*                  collen,             /**< space allocated to store indices of nonzeros for each column */
   int*                  colnnz,             /**< number of nonzero entries for each column */
   int*                  nzcount,            /**< counter for total number of nonzeros; should be increased whenever some colnnz is increased */
   int                   length,             /**< length of quadratic part */
   SCIP_QUADELEM*        quadelems           /**< quadratic elements */
   )
{
   int pos;

   SCIPdebugMessage("%p hess lag sparsity set nzflag for quad\n", (void*)oracle);

   assert(oracle != NULL);
   assert(colnz  != NULL);
   assert(collen != NULL);
   assert(colnnz != NULL);
   assert(nzcount != NULL);
   assert(quadelems != NULL);
   assert(length >= 0);

   for( ; length > 0; --length, ++quadelems )
   {
      assert(quadelems->idx1 <= quadelems->idx2);

      if( colnz[quadelems->idx2] == NULL || !SCIPsortedvecFindInt(colnz[quadelems->idx2], quadelems->idx1, colnnz[quadelems->idx2], &pos) )
      {
         SCIP_CALL( ensureIntArraySize(oracle->blkmem, &colnz[quadelems->idx2], &collen[quadelems->idx2], colnnz[quadelems->idx2]+1) );
         SCIPsortedvecInsertInt(colnz[quadelems->idx2], quadelems->idx1, &colnnz[quadelems->idx2], NULL);
         ++(*nzcount);
      }
   }

   return SCIP_OKAY;
}

/** collects indices of nonzero entries in the lower-left part of the hessian matrix of an expression
 * adds the indices to a given set of indices, avoiding duplicates */
static
SCIP_RETCODE hessLagSparsitySetNzFlagForExprtree(
   SCIP_NLPIORACLE*      oracle,             /**< NLPI oracle */
   int**                 colnz,              /**< indices of nonzero entries for each column */
   int*                  collen,             /**< space allocated to store indices of nonzeros for each column */
   int*                  colnnz,             /**< number of nonzero entries for each column */
   int*                  nzcount,            /**< counter for total number of nonzeros; should be increased when nzflag is set to 1 the first time */
   int*                  exprvaridx,         /**< indices of variables from expression tree in NLP */
   SCIP_EXPRTREE*        exprtree,           /**< expression tree */
   int                   dim                 /**< dimension of matrix */
   )
{
   SCIP_Real*  x;
   SCIP_Bool*  hesnz;
   int         i;
   int         j;
   int         nvars;
   int         nn;
   int         row;
   int         col;
   int         pos;

   assert(oracle != NULL);
   assert(colnz  != NULL);
   assert(collen != NULL);
   assert(colnnz != NULL);
   assert(nzcount != NULL);
   assert(exprvaridx != NULL);
   assert(exprtree != NULL);
   assert(dim >= 0);

   SCIPdebugMessage("%p hess lag sparsity set nzflag for exprtree\n", (void*)oracle);

   nvars = SCIPexprtreeGetNVars(exprtree);
   nn = nvars * nvars;

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &x,     nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &hesnz, nn) );

   for( i = 0; i < nvars; ++i )
      x[i] = 2.0; /* hope that this value does not make much trouble for the evaluation routines */  /*lint !e644*/

   SCIP_CALL( SCIPexprintHessianSparsityDense(oracle->exprinterpreter, exprtree, x, hesnz) );  /*lint !e644*/

   for( i = 0; i < nvars; ++i ) /* rows */
      for( j = 0; j <= i; ++j ) /* cols */
      {
         if( !hesnz[i*nvars + j] )
            continue;

         row = MAX(exprvaridx[i], exprvaridx[j]);
         col = MIN(exprvaridx[i], exprvaridx[j]);

         assert(row <  dim);
         assert(col <= row);

         if( colnz[row] == NULL || !SCIPsortedvecFindInt(colnz[row], col, colnnz[row], &pos) )
         {
            SCIP_CALL( ensureIntArraySize(oracle->blkmem, &colnz[row], &collen[row], colnnz[row]+1) );
            SCIPsortedvecInsertInt(colnz[row], col, &colnnz[row], NULL);
            ++(*nzcount);
         }
      }

   BMSfreeBlockMemoryArray(oracle->blkmem, &x, nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &hesnz, nn);

   return SCIP_OKAY;
}

/** adds quadratic part into hessian structure */
static
SCIP_RETCODE hessLagAddQuad(
   SCIP_Real             weight,             /**< weight of quadratic part */
   int                   length,             /**< number of elements in matrix of quadratic part */
   SCIP_QUADELEM*        quadelems,          /**< elements in matrix of quadratic part */
   int*                  hesoffset,          /**< row offsets in sparse matrix that is to be filled */ 
   int*                  hescol,             /**< column indices in sparse matrix that is to be filled */
   SCIP_Real*            values              /**< buffer for values of sparse matrix that is to be filled */
   )
{
   int idx;

   SCIPdebugMessage("hess lag add quad\n");

   assert(length >= 0);
   assert(quadelems != NULL || length == 0);
   assert(hesoffset != NULL);
   assert(hescol != NULL);
   assert(values != NULL);

   for( ; length > 0; --length, ++quadelems )  /*lint !e613*/
   {
      assert(quadelems->idx1 <= quadelems->idx2);  /*lint !e613*/
      if( !SCIPsortedvecFindInt(&hescol[hesoffset[quadelems->idx2]], quadelems->idx1, hesoffset[quadelems->idx2 + 1] - hesoffset[quadelems->idx2], &idx) )  /*lint !e613*/
      {
         SCIPerrorMessage("Could not find entry in hessian sparsity\n");
         return SCIP_ERROR;
      }
      values[hesoffset[quadelems->idx2] + idx] += weight * ((quadelems->idx1 == quadelems->idx2) ? 2 * quadelems->coef : quadelems->coef);  /*lint !e613*/
   }

   return SCIP_OKAY;
}

/** adds hessian of an expression into hessian structure */
static
SCIP_RETCODE hessLagAddExprtree(
   SCIP_NLPIORACLE*      oracle,             /**< oracle */
   SCIP_Real             weight,             /**< weight of quadratic part */
   const SCIP_Real*      x,                  /**< point for which hessian should be returned */
   SCIP_Bool             new_x,              /**< whether point has been evaluated before */
   int*                  exprvaridx,         /**< NLP indices for variables in expression tree */
   SCIP_EXPRTREE*        exprtree,           /**< expression tree */
   int*                  hesoffset,          /**< row offsets in sparse matrix that is to be filled */
   int*                  hescol,             /**< column indices in sparse matrix that is to be filled */
   SCIP_Real*            values              /**< buffer for values of sparse matrix that is to be filled */
   )
{
   SCIP_Real* xx;
   SCIP_Real* h;
   SCIP_Real* hh;
   int        i;
   int        j;
   int        nvars;
   int        nn;
   int        row;
   int        col;
   int        idx;
   SCIP_Real  val;

   SCIPdebugMessage("%p hess lag add exprtree\n", (void*)oracle);

   assert(oracle != NULL);
   assert(x != NULL || new_x == FALSE);

   nvars = exprtree != NULL ? SCIPexprtreeGetNVars(exprtree) : 0;
   if( nvars == 0 )
      return SCIP_OKAY;

   assert(exprtree != NULL);
   assert(exprvaridx != NULL);
   assert(hesoffset != NULL);
   assert(hescol != NULL);
   assert(values != NULL);

   nn = nvars * nvars;

   xx = NULL;
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &h, nn) );

   if( new_x )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &xx, nvars) );
      for( i = 0; i < nvars; ++i )
      {
         assert(exprvaridx[i] >= 0);
         xx[i] = x[exprvaridx[i]];  /*lint !e613*/
      }
   }

   SCIP_CALL( SCIPexprintHessianDense(oracle->exprinterpreter, exprtree, xx, new_x, &val, h) );  /*lint !e644*/
   if( val != val )  /*lint !e777*/
   {
      SCIPdebugMessage("hessian evaluation yield invalid function value %g\n", val);
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
      BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
      return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
   }

   hh = h;
   for( i = 0; i < nvars; ++i ) /* rows */
   {
      for( j = 0; j <= i; ++j, ++hh ) /* cols */
      {
         if( !*hh )
            continue;

         if( !SCIPisFinite(*hh) )  /*lint !e777*/
         {
            SCIPdebugMessage("hessian evaluation yield invalid hessian value %g\n", *hh);
            BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
            BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
            return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
         }

         row = MAX(exprvaridx[i], exprvaridx[j]);
         col = MIN(exprvaridx[i], exprvaridx[j]);

         if( !SCIPsortedvecFindInt(&hescol[hesoffset[row]], col, hesoffset[row+1] - hesoffset[row], &idx) )
         {
            SCIPerrorMessage("Could not find entry (%d, %d) in hessian sparsity\n", row, col);
            BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
            BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);
            return SCIP_ERROR;
         }

         values[hesoffset[row] + idx] += weight * *hh;
      }
      hh += nvars - j;
   }

   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &h, nn);

   return SCIP_OKAY;
}

/** prints a name, if available, makes sure it has not more than 64 characters, and adds a unique prefix if the longnames flag is set */
static
void printName(
   char*                 buffer,             /**< buffer to print to, has to be not NULL */
   char*                 name,               /**< name, or NULL */
   int                   idx,                /**< index of var or cons which the name corresponds to */
   char                  prefix,             /**< a letter (typically 'x' or 'e') to distinguish variable and equation names, if names[idx] is not available */
   const char*           suffix,             /**< a suffer to add to the name, or NULL */
   SCIP_Bool             longnames           /**< whether prefixes for long names should be added */
   )
{
   assert(idx >= 0 && idx < 100000); /* to ensure that we do not exceed the size of the buffer */

   if( longnames )
   {
      if( name != NULL )
         sprintf(buffer, "%c%05d%.*s%s", prefix, idx, suffix ? (int)(57-strlen(suffix)) : 57, name, suffix ? suffix : "");
      else
         sprintf(buffer, "%c%05d", prefix, idx);
   }
   else
   {
      if( name != NULL )
      {
         assert(strlen(name) + (suffix ? strlen(suffix) : 0) <= 64);
         sprintf(buffer, "%s%s", name, suffix ? suffix : "");
      }
      else
         sprintf(buffer, "%c%d%s", prefix, idx, suffix ? suffix : "");
   }
}

/** prints a function */
static
SCIP_RETCODE printFunction(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, has to be not NULL */
   SCIP_NLPIORACLECONS*  cons,               /**< constraint which function to print */
   SCIP_Bool             longvarnames,       /**< whether variable names need to be shorten to 64 characters */
   SCIP_Bool             longequnames        /**< whether equation names need to be shorten to 64 characters */
   )
{  /*lint --e{715}*/
   int i;
   char namebuf[70];

   SCIPdebugMessage("%p print function\n", (void*)oracle);

   assert(oracle != NULL);
   assert(file != NULL);
   assert(cons != NULL);

   for( i = 0; i < cons->nlinidxs; ++i )
   {
      printName(namebuf, oracle->varnames != NULL ? oracle->varnames[cons->linidxs[i]] : NULL, cons->linidxs[i], 'x', NULL, longvarnames);
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.20g*%s", cons->lincoefs[i], namebuf);
      if( i % 10 == 9 )
         SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }

   for( i = 0; i < cons->nquadelems; ++i )
   {
      printName(namebuf, oracle->varnames != NULL ? oracle->varnames[cons->quadelems[i].idx1] : NULL, cons->quadelems[i].idx1, 'x', NULL, longvarnames);
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.20g*%s", cons->quadelems[i].coef, namebuf);
      printName(namebuf, oracle->varnames != NULL ? oracle->varnames[cons->quadelems[i].idx2] : NULL, cons->quadelems[i].idx2, 'x', NULL, longvarnames);
      SCIPmessageFPrintInfo(messagehdlr, file, "*%s", namebuf);
      if( i % 10 == 9 )
         SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }

   if( cons->exprtree != NULL )
   {
      char** varnames;
      SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &varnames, SCIPexprtreeGetNVars(cons->exprtree)) );  /*lint !e666*/

      /* setup variable names */
      for( i = 0; i < SCIPexprtreeGetNVars(cons->exprtree); ++i )
      {
         assert(cons->exprvaridxs[i] < 1e+20);
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &varnames[i], 70) );  /*lint !e866 !e506 !e644*/
         printName(varnames[i], oracle->varnames != NULL ? oracle->varnames[cons->exprvaridxs[i]] : NULL, cons->exprvaridxs[i], 'x', NULL, longvarnames);
      }

      SCIPmessageFPrintInfo(messagehdlr, file, " +");
      SCIPexprtreePrint(cons->exprtree, messagehdlr, file, (const char**)varnames, NULL);

      for( i = 0; i < SCIPexprtreeGetNVars(cons->exprtree); ++i )
      {
         BMSfreeBlockMemoryArray(oracle->blkmem, &varnames[i], 70);  /*lint !e866*/
      }
      BMSfreeBlockMemoryArray(oracle->blkmem, &varnames, SCIPexprtreeGetNVars(cons->exprtree));
   }

   return SCIP_OKAY;
}

/** returns whether an expression is contains nonsmooth operands (min, max, abs, ...) */
static
SCIP_Bool exprIsNonSmooth(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int i;

   assert(expr != NULL);
   assert(SCIPexprGetChildren(expr) != NULL || SCIPexprGetNChildren(expr) == 0);

   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      if( exprIsNonSmooth(SCIPexprGetChildren(expr)[i]) )
         return TRUE;
   }

   switch( SCIPexprGetOperator(expr) )
   {
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   case SCIP_EXPR_SIGNPOWER:
      return TRUE;

   default: ;
   } /*lint !e788*/

   return FALSE;
}

/**@} */

/**@name public function */
/**@{ */

/** creates an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPIORACLE**     oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   assert(blkmem != NULL);
   assert(oracle != NULL);

   SCIPdebugMessage("%p oracle create\n", (void*)oracle);

   SCIP_ALLOC( BMSallocMemory(oracle) );
   BMSclearMemory(*oracle);

   (*oracle)->blkmem   = blkmem;
   (*oracle)->infinity = SCIP_DEFAULT_INFINITY;
   (*oracle)->vardegreesuptodate = TRUE;

   SCIPdebugMessage("Oracle initializes expression interpreter %s\n", SCIPexprintGetName());
   SCIP_CALL( SCIPexprintCreate(blkmem, &(*oracle)->exprinterpreter) );

   /* create zero objective function */
   SCIP_CALL( createConstraint((*oracle)->blkmem, &(*oracle)->objective, 0, NULL, NULL, 0, NULL, NULL, NULL, 0.0, 0.0, NULL) );

   return SCIP_OKAY;
}

/** frees an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleFree(
   SCIP_NLPIORACLE**     oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle  != NULL);
   assert(*oracle != NULL);

   SCIPdebugMessage("%p oracle free\n", (void*)oracle);

   invalidateJacobiSparsity(*oracle);
   invalidateHessianLagSparsity(*oracle);

   freeConstraint((*oracle)->blkmem, &(*oracle)->objective);
   freeConstraints(*oracle);
   freeVariables(*oracle);

   SCIP_CALL( SCIPexprintFree(&(*oracle)->exprinterpreter) );

   if( (*oracle)->name != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleSetProblemName(*oracle, NULL) );
   }

   BMSfreeMemory(oracle);

   return SCIP_OKAY;
}

/** sets the value for infinity */
SCIP_RETCODE SCIPnlpiOracleSetInfinity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real             infinity            /**< value to use for infinity */
   )
{
   assert(oracle != NULL);
   assert(infinity > 0.0);

   SCIPdebugMessage("%p set infinity\n", (void*)oracle);

   oracle->infinity = infinity;

   return SCIP_OKAY;
}

/** gets the value for infinity */
SCIP_Real SCIPnlpiOracleGetInfinity(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p get infinity\n", (void*)oracle);

   return oracle->infinity;
}

/** sets the problem name (used for printing) */
SCIP_RETCODE SCIPnlpiOracleSetProblemName(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const char*           name                /**< name of problem */
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p set problem name\n", (void*)oracle);

   if( oracle->name != NULL )
   {
      BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->name, strlen(oracle->name)+1);
   }

   if( name != NULL )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->name, name, strlen(name)+1) );
   }

   return SCIP_OKAY;
}

/** gets the problem name, or NULL if none set */
const char* SCIPnlpiOracleGetProblemName(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p get problem name\n", (void*)oracle);

   return oracle->name;
}

/** adds variables */
SCIP_RETCODE SCIPnlpiOracleAddVars(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to add */
   const SCIP_Real*      lbs,                /**< array with lower bounds of new variables, or NULL if all -infinity */
   const SCIP_Real*      ubs,                /**< array with upper bounds of new variables, or NULL if all +infinity */
   const char**          varnames            /**< array with names of new variables, or NULL if no names should be stored */
   )
{
   int i;

   assert(oracle != NULL);

   SCIPdebugMessage("%p add vars\n", (void*)oracle);

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(nvars > 0);

   SCIP_CALL( ensureVarsSize(oracle, oracle->nvars + nvars) );

   if( lbs != NULL )
   {
      BMScopyMemoryArray(&oracle->varlbs[oracle->nvars], lbs, nvars);  /*lint !e866*/
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varlbs[oracle->nvars+i] = -oracle->infinity;

   if( ubs != NULL )
   {
      BMScopyMemoryArray(&oracle->varubs[oracle->nvars], ubs, nvars);  /*lint !e866*/

      /* ensure variable bounds are consistent */
      for( i = oracle->nvars; i < oracle->nvars + nvars; ++i )
      {
         if( oracle->varlbs[i] > oracle->varubs[i] )
         {
            assert(EPSEQ(oracle->varlbs[i], oracle->varubs[i], SCIP_DEFAULT_EPSILON));
            oracle->varlbs[i] = oracle->varubs[i];
         }
      }
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varubs[oracle->nvars+i] =  oracle->infinity;

   if( varnames != NULL )
   {
      if( oracle->varnames == NULL )
      {
         SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->varnames, oracle->varssize) );
         BMSclearMemoryArray(oracle->varnames, oracle->nvars);
      }

      for( i = 0; i < nvars; ++i )
      {
         if( varnames[i] != NULL )
         {
            SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &oracle->varnames[oracle->nvars+i], varnames[i], strlen(varnames[i])+1) );  /*lint !e866*/
         }
         else
            oracle->varnames[oracle->nvars+i] = NULL;
      }
   }
   else if( oracle->varnames != NULL )
   {
      BMSclearMemoryArray(&oracle->varnames[oracle->nvars], nvars);  /*lint !e866*/
   }

   BMSclearMemoryArray(&oracle->vardegrees[oracle->nvars], nvars);  /*lint !e866*/

   /* @TODO update sparsity pattern by extending heslagoffsets */
   invalidateHessianLagSparsity(oracle);

   oracle->nvars += nvars;

   return SCIP_OKAY;
}

/** adds constraints 
 * 
 *  linear coefficients: row(=constraint) oriented matrix;
 *  quadratic coefficients: row oriented matrix for each constraint
 */
SCIP_RETCODE SCIPnlpiOracleAddConstraints(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to add */
   const SCIP_Real*      lhss,               /**< array with left-hand sides of constraints, or NULL if all -infinity */
   const SCIP_Real*      rhss,               /**< array with right-hand sides of constraints, or NULL if all +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   const int*            nquadelems,         /**< number of elements in matrix of quadratic part for each constraint,
                                              * may be NULL in case of no quadratic part in any constraint */
   SCIP_QUADELEM* const* quadelems,          /**< quadratic elements specifying quadratic part for each constraint, entry of array may be NULL in case of no quadratic part,
                                              * may be NULL in case of no quadratic part in any constraint */
   int* const*           exprvaridxs,        /**< NULL if no nonquadratic parts, otherwise epxrvaridxs[.] maps variable indices in expression tree to indices in nlp */
   SCIP_EXPRTREE* const* exprtrees,          /**< NULL if no nonquadratic parts, otherwise exprtrees[.] gives nonquadratic part, 
                                              *   or NULL if no nonquadratic part in this constraint */
   const char**          consnames           /**< names of new constraints, or NULL if no names should be stored */
   )
{  /*lint --e{715}*/
   SCIP_NLPIORACLECONS* cons;
   SCIP_Bool addednlcon;  /* whether a nonlinear constraint was added */
   int c;

   assert(oracle != NULL);

   SCIPdebugMessage("%p add constraints\n", (void*)oracle);

   if( nconss == 0 )
      return SCIP_OKAY;

   assert(nconss > 0);

   addednlcon = FALSE;

   invalidateJacobiSparsity(oracle); /* @TODO we could also update (extend) the sparsity pattern */

   SCIP_CALL( ensureConssSize(oracle, oracle->nconss + nconss) );
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( createConstraint(oracle->blkmem, &cons,
            nlininds != NULL ? nlininds[c] : 0,
            lininds != NULL ? lininds[c] : NULL,
            linvals != NULL ? linvals[c] : NULL,
            nquadelems != NULL ? nquadelems[c] : 0,
            quadelems != NULL ? quadelems[c] : NULL,
            exprvaridxs != NULL ? exprvaridxs[c] : NULL,
            exprtrees != NULL ? exprtrees[c] : NULL,
            lhss != NULL ? lhss[c] : -oracle->infinity,
            rhss != NULL ? rhss[c] :  oracle->infinity,
            consnames != NULL ? consnames[c] : NULL
            ) );

      if( cons->nquadelems > 0 )
         addednlcon = TRUE;

      if( cons->exprtree != NULL )
      {
         addednlcon = TRUE;
         SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, cons->exprtree) );
      }

      /* keep variable degrees updated */
      if( oracle->vardegreesuptodate )
         updateVariableDegreesCons(oracle, cons);

      oracle->conss[oracle->nconss+c] = cons;
   }
   oracle->nconss += nconss;

   if( addednlcon == TRUE )
      invalidateHessianLagSparsity(oracle);

   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected
 * 
 *  May change sparsity pattern.
 */
SCIP_RETCODE SCIPnlpiOracleSetObjective(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real       constant,           /**< constant part of objective */
   int                   nlin,               /**< number of linear variable coefficients */ 
   const int*            lininds,            /**< indices of linear variables, or NULL if no linear part */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables, or NULL if no linear part */
   int                   nquadelems,         /**< number of entries in matrix of quadratic part */
   const SCIP_QUADELEM*  quadelems,          /**< entries in matrix of quadratic part, may be NULL in case of no quadratic part */
   const int*            exprvaridxs,        /**< maps variable indices in expression tree to indices in nlp, or NULL if no nonquadratic part */
   const SCIP_EXPRTREE*  exprtree            /**< expression tree of nonquadratic part, or NULL if no nonquadratic part */
   )
{  /*lint --e{715}*/
   assert(oracle != NULL);
   assert(REALABS(constant) < oracle->infinity);

   SCIPdebugMessage("%p set objective\n", (void*)oracle);

   if( nquadelems > 0 || oracle->objective->quadsize > 0 || exprtree != NULL || oracle->objective->exprtree != NULL )
      invalidateHessianLagSparsity(oracle);

   /* clear previous objective */
   freeConstraint(oracle->blkmem, &oracle->objective);

   SCIP_CALL( createConstraint(oracle->blkmem, &oracle->objective,
         nlin, lininds, linvals, nquadelems, quadelems, exprvaridxs, exprtree, constant, constant, NULL) );

   if( oracle->objective->exprtree != NULL )
   {
      SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, oracle->objective->exprtree) );
   }

   oracle->vardegreesuptodate = FALSE;

   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiOracleChgVarBounds(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds, or NULL if all should be -infty */
   const SCIP_Real*      ubs                 /**< new upper bounds, or NULL if all should be +infty */
   )
{
   int i;

   assert(oracle != NULL);
   assert(indices != NULL || nvars == 0);

   SCIPdebugMessage("%p chg var bounds\n", (void*)oracle);

   for( i = 0; i < nvars; ++i )
   {
      assert(indices != NULL);
      assert(indices[i] >= 0);
      assert(indices[i] < oracle->nvars);

      oracle->varlbs[indices[i]] = (lbs != NULL ? lbs[i] : -oracle->infinity);
      oracle->varubs[indices[i]] = (ubs != NULL ? ubs[i] :  oracle->infinity);

      if( oracle->varlbs[indices[i]] > oracle->varubs[indices[i]] )
      {
         /* inconsistent bounds; let's assume it's due to rounding and make them equal */
         assert(EPSEQ(oracle->varlbs[indices[i]], oracle->varubs[indices[i]], SCIP_DEFAULT_EPSILON));
         oracle->varlbs[indices[i]] = oracle->varubs[indices[i]];
      }
   }

   return SCIP_OKAY;
}

/** change constraint bounds */
SCIP_RETCODE SCIPnlpiOracleChgConsSides(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to change bounds */
   const int*            indices,            /**< indices of constraints to change bounds */
   const SCIP_Real*      lhss,               /**< new left-hand sides, or NULL if all should be -infty */
   const SCIP_Real*      rhss                /**< new right-hand sides, or NULL if all should be +infty */
   )
{
   int i;

   assert(oracle != NULL);
   assert(indices != NULL || nconss == 0);

   SCIPdebugMessage("%p chg cons sides\n", (void*)oracle);

   for( i = 0; i < nconss; ++i )
   {
      assert(indices != NULL);
      assert(indices[i] >= 0);
      assert(indices[i] < oracle->nconss);

      oracle->conss[indices[i]]->lhs = (lhss != NULL ? lhss[i] : -oracle->infinity);
      oracle->conss[indices[i]]->rhs = (rhss != NULL ? rhss[i] :  oracle->infinity);
      if( oracle->conss[indices[i]]->lhs > oracle->conss[indices[i]]->rhs )
      {
         assert(EPSEQ(oracle->conss[indices[i]]->lhs, oracle->conss[indices[i]]->rhs, SCIP_DEFAULT_EPSILON));
         oracle->conss[indices[i]]->lhs = oracle->conss[indices[i]]->rhs;
      }
   }

   return SCIP_OKAY;
}

/** deletes a set of variables */
SCIP_RETCODE SCIPnlpiOracleDelVarSet(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< deletion status of vars in input (1 if var should be deleted, 0 if not); 
                                              *   new position of var in output (-1 if var was deleted) */
   )
{  /*lint --e{715}*/
   int c;
   int lastgood; /* index of the last variable that should be kept */
   SCIP_NLPIORACLECONS* cons;

   assert(oracle != NULL);

   SCIPdebugMessage("%p del var set\n", (void*)oracle);

   invalidateJacobiSparsity(oracle);
   invalidateHessianLagSparsity(oracle);

   lastgood = oracle->nvars - 1;
   while( lastgood >= 0 && delstats[lastgood] == 1 )
      --lastgood;
   if( lastgood < 0 )
   {
      /* all variables should be deleted */
      assert(oracle->nconss == 0); /* we could relax this by checking that all constraints are constant */
      assert(oracle->objective->exprtree == NULL || SCIPexprtreeGetNVars(oracle->objective->exprtree) == 0);
      oracle->objective->nquadelems = 0;
      oracle->objective->nlinidxs = 0;
      for( c = 0; c < oracle->nvars; ++c )
         delstats[c] = -1;
      freeVariables(oracle);
      return SCIP_OKAY;
   }

   /* delete variables at the end */
   for( c = oracle->nvars - 1; c > lastgood; --c )
   {
      if( oracle->varnames && oracle->varnames[c] != NULL )
      {
         BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->varnames[c], strlen(oracle->varnames[c])+1);  /*lint !e866*/
      }
      delstats[c] = -1;
   }

   /* go through variables from the beginning on
    * if variable should be deleted, free it and move lastgood variable to this position
    * then update lastgood */
   for( c = 0; c <= lastgood; ++c )
   {
      if( delstats[c] == 0 )
      { /* variable should not be deleted and is kept on position c */
         delstats[c] = c;
         continue;
      }
      assert(delstats[c] == 1); /* variable should be deleted */

      if( oracle->varnames && oracle->varnames[c] != NULL )
      {
         BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->varnames[c], strlen(oracle->varnames[c])+1);  /*lint !e866*/
      }
      delstats[c] = -1;

      /* move constraint at position lastgood to position c */
      SCIP_CALL( moveVariable(oracle, lastgood, c) );
      delstats[lastgood] = c; /* mark that lastgood variable is now at position c */

      /* move lastgood forward, delete variables on the way */
      --lastgood;
      while( lastgood > c && delstats[lastgood] == 1)
      {
         if( oracle->varnames && oracle->varnames[lastgood] != NULL )
         {
            BMSfreeBlockMemoryArray(oracle->blkmem, &oracle->varnames[lastgood], strlen(oracle->varnames[lastgood])+1);  /*lint !e866*/
         }
         delstats[lastgood] = -1;
         --lastgood;
      }
   }
   assert(c == lastgood);

   for( c = -1; c < oracle->nconss; ++c )
   {
      cons = c < 0 ? oracle->objective : oracle->conss[c];
      assert(cons != NULL);

      /* update indices in linear part, sort indices, and then clear elements that are marked as deleted */
      mapIndices(delstats, cons->nlinidxs, cons->linidxs);
      SCIPsortIntReal(cons->linidxs, cons->lincoefs, cons->nlinidxs);
      clearDeletedLinearElements(oracle->blkmem, &cons->linidxs, &cons->lincoefs, &cons->nlinidxs);

      /* update indices in quadratic part, sort elements, and then clear elements that are marked as deleted */
      mapIndicesQuad(delstats, cons->quadsize, cons->quadelems);
      SCIPquadelemSort(cons->quadelems, cons->quadsize);
      clearDeletedQuadElements(oracle->blkmem, &cons->quadelems, &cons->quadsize);

      if( cons->exprtree != NULL )
      {
         mapIndices(delstats, SCIPexprtreeGetNVars(cons->exprtree), cons->exprvaridxs);
         /* assert that all variables from this expression have been deleted */
         assert(SCIPexprtreeGetNVars(cons->exprtree) == 0 || cons->exprvaridxs[SCIPexprtreeGetNVars(cons->exprtree)-1] == -1);
         BMSfreeBlockMemoryArrayNull(oracle->blkmem, &cons->exprvaridxs, SCIPexprtreeGetNVars(cons->exprtree));
         SCIP_CALL( SCIPexprtreeFree(&cons->exprtree) );
      }
   }

   oracle->nvars = lastgood+1;

   return SCIP_OKAY;
}

/** deletes a set of constraints */
SCIP_RETCODE SCIPnlpiOracleDelConsSet(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of rows in input (1 if row should be deleted, 0 if not); 
                                              *   new position of row in output (-1 if row was deleted) */
   )
{  /*lint --e{715}*/
   int c;
   int lastgood; /* index of the last constraint that should be kept */ 

   assert(oracle != NULL);

   SCIPdebugMessage("%p del cons set\n", (void*)oracle);

   invalidateJacobiSparsity(oracle);
   invalidateHessianLagSparsity(oracle);
   oracle->vardegreesuptodate = FALSE;

   lastgood = oracle->nconss - 1;
   while( lastgood >= 0 && delstats[lastgood] == 1)
      --lastgood;
   if( lastgood < 0 )
   {
      /* all constraints should be deleted */
      for( c = 0; c < oracle->nconss; ++c )
         delstats[c] = -1;
      freeConstraints(oracle);
      return SCIP_OKAY;
   }

   /* delete constraints at the end */
   for( c = oracle->nconss - 1; c > lastgood; --c )
   {
      freeConstraint(oracle->blkmem, &oracle->conss[c]);
      assert(oracle->conss[c] == NULL);
      delstats[c] = -1;
   }

   /* go through constraint from the beginning on
    * if constraint should be deleted, free it and move lastgood constraint to this position
    * then update lastgood */
   for( c = 0; c <= lastgood; ++c )
   {
      if( delstats[c] == 0 )
      {
         /* constraint should not be deleted and is kept on position c */
         delstats[c] = c;
         continue;
      }
      assert(delstats[c] == 1); /* constraint should be deleted */

      freeConstraint(oracle->blkmem, &oracle->conss[c]);
      assert(oracle->conss[c] == NULL);
      delstats[c] = -1;

      /* move constraint at position lastgood to position c */
      oracle->conss[c] = oracle->conss[lastgood];
      assert(oracle->conss[c] != NULL);
      delstats[lastgood] = c; /* mark that lastgood constraint is now at position c */
      oracle->conss[lastgood] = NULL;

      /* move lastgood forward, delete constraints on the way */
      --lastgood;
      while( lastgood > c && delstats[lastgood] == 1)
      {
         freeConstraint(oracle->blkmem, &oracle->conss[lastgood]);
         assert(oracle->conss[lastgood] == NULL);
         delstats[lastgood] = -1;
         --lastgood;
      }
   }
   assert(c == lastgood);

   oracle->nconss = lastgood+1;

   return SCIP_OKAY;
}

/** changes linear coefficients in one constraint or objective */
SCIP_RETCODE SCIPnlpiOracleChgLinearCoefs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where linear coefficients should be changed, or -1 for objective */
   int                   nentries,           /**< number of coefficients to change */
   const int*            varidxs,            /**< array with indices of variables which coefficients should be changed */
   const SCIP_Real*      newcoefs            /**< array with new coefficients of variables */
   )
{  /*lint --e{715}*/
   SCIP_NLPIORACLECONS* cons;
   SCIP_Bool needsort;
   int       i;

   SCIPdebugMessage("%p chg linear coefs\n", (void*)oracle);

   assert(oracle != NULL);
   assert(varidxs != NULL || nentries == 0);
   assert(newcoefs != NULL || nentries == 0);
   assert(considx >= -1);
   assert(considx < oracle->nconss);

   if( nentries == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("change %d linear coefficients in cons %d\n", nentries, considx);

   needsort = FALSE;

   cons = considx < 0 ? oracle->objective : oracle->conss[considx];

   if( cons->linsize == 0 )
   {
      /* first time we have linear coefficients in this constraint (or objective) */
      assert(cons->linidxs  == NULL);
      assert(cons->lincoefs == NULL);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &cons->linidxs,  varidxs,  nentries) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &cons->lincoefs, newcoefs, nentries) );
      cons->linsize  = nentries;
      cons->nlinidxs = nentries;

      needsort = TRUE;
   }
   else
   {
      int pos;

      for( i = 0; i < nentries; ++i )
      {
         assert(varidxs[i] >= 0);             /*lint !e613*/
         assert(varidxs[i] < oracle->nvars);  /*lint !e613*/

         if( SCIPsortedvecFindInt(cons->linidxs, varidxs[i], cons->nlinidxs, &pos) )  /*lint !e613*/
         {
            SCIPdebugMessage("replace coefficient of var %d at pos %d by %g\n", varidxs[i], pos, newcoefs[i]);  /*lint !e613*/

            cons->lincoefs[pos] = newcoefs[i];  /*lint !e613*/

            /* remember that we need to sort/merge/squeeze array if coefficient became zero here */
            needsort |= (newcoefs[i] == 0.0);  /*lint !e613 !e514*/
         }
         else if( newcoefs[i] != 0.0 )  /*lint !e613*/
         {
            /* append new entry */
            SCIPdebugMessage("add coefficient of var %d at pos %d, value %g\n", varidxs[i], cons->nlinidxs, newcoefs[i]);  /*lint !e613*/

            SCIP_CALL( ensureConsLinSize(oracle->blkmem, cons, cons->nlinidxs + (nentries-i)) );
            cons->linidxs[cons->nlinidxs]  = varidxs[i];   /*lint !e613*/
            cons->lincoefs[cons->nlinidxs] = newcoefs[i];  /*lint !e613*/
            ++cons->nlinidxs;

            needsort = TRUE;
         }
      }
   }

   if( needsort )
   {
      int oldlen;

      invalidateJacobiSparsity(oracle);

      oldlen = cons->nlinidxs;
      sortLinearCoefficients(&cons->nlinidxs, cons->linidxs, cons->lincoefs);

      /* if sorting removed an entry, then the var degrees are not uptodate anymore */
      oracle->vardegreesuptodate &= (cons->nlinidxs == oldlen);  /*lint !e514*/

      /* increase variable degrees of variables to 1 */
      if( oracle->vardegreesuptodate )
         for( i = 0; i < cons->nlinidxs; ++i )
            oracle->vardegrees[varidxs[i]] = MAX(1, oracle->vardegrees[varidxs[i]]);  /*lint !e613*/
   }

   return SCIP_OKAY;
}

/** changes (or adds) coefficients in the quadratic part of one constraint or objective */
SCIP_RETCODE SCIPnlpiOracleChgQuadCoefs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where quadratic coefficients should be changed, or -1 for objective */
   int                   nquadelems,         /**< number of entries in quadratic constraint to change */
   const SCIP_QUADELEM*  quadelems           /**< new elements in quadratic matrix (replacing already existing ones or adding new ones) */
   )
{  /*lint --e{715}*/
   SCIP_NLPIORACLECONS* cons;
   SCIP_Bool needsort;
   int       i;

   SCIPdebugMessage("%p chg quad coefs\n", (void*)oracle);

   assert(oracle != NULL);
   assert(quadelems != NULL || nquadelems == 0);
   assert(considx >= -1);
   assert(considx < oracle->nconss);

   if( nquadelems == 0 )
      return SCIP_OKAY;

   needsort = FALSE;

   cons = considx < 0 ? oracle->objective : oracle->conss[considx];

   if( cons->quadsize == 0 )
   {
      /* first time we have quadratic coefficients in this constraint (or objective) */
      assert(cons->quadelems == NULL);

      SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &cons->quadelems, quadelems, nquadelems) );
      cons->quadsize  = nquadelems;
      cons->nquadelems = nquadelems;

      needsort = TRUE;
   }
   else
   {
      int pos;

      for( i = 0; i < nquadelems; ++i )
      {
         assert(quadelems[i].idx1 >= 0);  /*lint !e613*/
         assert(quadelems[i].idx2 >= 0);  /*lint !e613*/
         assert(quadelems[i].idx1 < oracle->nvars);  /*lint !e613*/
         assert(quadelems[i].idx2 < oracle->nvars);  /*lint !e613*/

         /* if we already have an entry for quadelems[i], then just replace the coefficient, otherwise append new entry */
         if( SCIPquadelemSortedFind(cons->quadelems, quadelems[i].idx1, quadelems[i].idx2, cons->nquadelems, &pos) )  /*lint !e613*/
         {
            SCIPdebugMessage("replace coefficient of var%d*var%d at pos %d by %g\n", quadelems[i].idx1, quadelems[i].idx2, pos, quadelems[i].coef);  /*lint !e613*/

            cons->quadelems[pos].coef = quadelems[i].coef;  /*lint !e613*/

            /* remember that we need to sort/merge/squeeze array if coefficient became zero here */
            needsort |= (quadelems[i].coef == 0.0);  /*lint !e613 !e514*/
         }
         else
         {
            /* append new entry */
            SCIPdebugMessage("add coefficient of var%d*var%d at pos %d, value %g\n", quadelems[i].idx1, quadelems[i].idx2, cons->nquadelems, quadelems[i].coef);  /*lint !e613*/

            SCIP_CALL( ensureConsQuadSize(oracle->blkmem, cons, cons->nquadelems + (nquadelems-i)) );
            cons->quadelems[cons->nquadelems] = quadelems[i];  /*lint !e613*/
            ++cons->nquadelems;

            needsort = TRUE;
         }
      }
   }

   if( needsort )
   {
      int oldsize;

      invalidateJacobiSparsity(oracle);
      invalidateHessianLagSparsity(oracle);

      oldsize = cons->nquadelems;
      SCIPquadelemSort(cons->quadelems, cons->nquadelems);
      SCIPquadelemSqueeze(cons->quadelems, cons->nquadelems, &cons->nquadelems);

      /* if sorting removed an entry, then the var degrees are not uptodate anymore */
      oracle->vardegreesuptodate &= (cons->nquadelems == oldsize);  /*lint !e514*/

      /* increase variable degrees of variables to 2 */
      if( oracle->vardegreesuptodate )
         for( i = 0; i < cons->nquadelems; ++i )
         {
            oracle->vardegrees[cons->quadelems[i].idx1] = MAX(2, oracle->vardegrees[cons->quadelems[i].idx1]);
            oracle->vardegrees[cons->quadelems[i].idx2] = MAX(2, oracle->vardegrees[cons->quadelems[i].idx2]);
         }
   }

   return SCIP_OKAY;
}

/** replaces expression tree of one constraint or objective  */
SCIP_RETCODE SCIPnlpiOracleChgExprtree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where expression tree should be changed, or -1 for objective */
   const int*            exprvaridxs,        /**< problem indices of variables in expression tree */
   const SCIP_EXPRTREE*  exprtree            /**< new expression tree, or NULL */
   )
{
   SCIP_NLPIORACLECONS* cons;
   int j;

   SCIPdebugMessage("%p chg exprtree\n", (void*)oracle);

   assert(oracle != NULL);
   assert(considx >= -1);
   assert(considx < oracle->nconss);
   assert((exprvaridxs != NULL) == (exprtree != NULL));

   invalidateHessianLagSparsity(oracle);
   invalidateJacobiSparsity(oracle);

   cons = considx < 0 ? oracle->objective : oracle->conss[considx];

   /* free previous expression tree */
   if( cons->exprtree != NULL )
   {
      BMSfreeBlockMemoryArray(oracle->blkmem, &cons->exprvaridxs, SCIPexprtreeGetNVars(cons->exprtree));
      SCIP_CALL( SCIPexprtreeFree(&cons->exprtree));
      oracle->vardegreesuptodate = FALSE;
   }

   /* if user did not want to set new tree, then we are done */
   if( exprtree == NULL )
      return SCIP_OKAY;

   assert(oracle->exprinterpreter != NULL);

   /* install new expression tree */
   SCIP_CALL( SCIPexprtreeCopy(oracle->blkmem, &cons->exprtree, (SCIP_EXPRTREE*)exprtree) );
   SCIP_CALL( SCIPexprintCompile(oracle->exprinterpreter, cons->exprtree) );
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(oracle->blkmem, &cons->exprvaridxs, exprvaridxs, SCIPexprtreeGetNVars(cons->exprtree)) );

   /* increase variable degree to keep them up to date
    * could get more accurate degree via getMaxDegree function in exprtree, but no solver would use this information so far
    */
   if( oracle->vardegreesuptodate )
      for( j = 0; j < SCIPexprtreeGetNVars(cons->exprtree); ++j )
      {
         assert(cons->exprvaridxs[j] >= 0);
         assert(cons->exprvaridxs[j] <  oracle->nvars);
         oracle->vardegrees[cons->exprvaridxs[j]] = INT_MAX;
      }

   return SCIP_OKAY;
}

/** changes one parameter of expression tree of one constraint or objective
 */
SCIP_RETCODE SCIPnlpiOracleChgExprParam(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where parameter should be changed in expression tree, or -1 for objective */
   int                   paramidx,           /**< index of parameter */
   SCIP_Real             paramval            /**< new value of parameter */
   )
{
   SCIPdebugMessage("%p chg expr param\n", (void*)oracle);

   assert(oracle != NULL);
   assert(considx >= -1);
   assert(considx < oracle->nconss);
   assert(paramidx >= 0);
   assert(considx >= 0  || oracle->objective->exprtree != NULL);
   assert(considx >= 0  || paramidx < SCIPexprtreeGetNParams(oracle->objective->exprtree));
   assert(considx == -1 || oracle->conss[considx]->exprtree != NULL);
   assert(considx == -1 || paramidx < SCIPexprtreeGetNParams(oracle->conss[considx]->exprtree));

   SCIPexprtreeSetParamVal(considx >= 0 ? oracle->conss[considx]->exprtree : oracle->objective->exprtree, paramidx, paramval);

   return SCIP_OKAY;
}

/** changes the constant value in the objective function
 */
SCIP_RETCODE SCIPnlpiOracleChgObjConstant(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real             objconstant         /**< new value for objective constant */
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p chg obj constant\n", (void*)oracle);

   oracle->objective->lhs = objconstant;
   oracle->objective->rhs = objconstant;

   return SCIP_OKAY;
}

/** gives the current number of variables */
int SCIPnlpiOracleGetNVars(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   return oracle->nvars;
}

/** gives the current number of constraints */
int SCIPnlpiOracleGetNConstraints(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   return oracle->nconss;
}

/** gives the variables lower bounds */
const SCIP_Real* SCIPnlpiOracleGetVarLbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   return oracle->varlbs;
}

/** gives the variables upper bounds */
const SCIP_Real* SCIPnlpiOracleGetVarUbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   return oracle->varubs;
}

/** gives the variables names, or NULL if not set */
char** SCIPnlpiOracleGetVarNames(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   return oracle->varnames;
}

/** Gives maximum degree of a variable w.r.t. objective and all constraints.
 *  The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
int SCIPnlpiOracleGetVarDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   varidx
   )
{
   assert(oracle != NULL);
   assert(varidx >= 0);
   assert(varidx < oracle->nvars);

   updateVariableDegrees(oracle);

   return oracle->vardegrees[varidx];
}

/** Gives maximum degree of all variables w.r.t. objective and all constraints.
 *  The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
int* SCIPnlpiOracleGetVarDegrees(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   updateVariableDegrees(oracle);

   return oracle->vardegrees;
}

/** gives left-hand side of a constraint */
SCIP_Real SCIPnlpiOracleGetConstraintLhs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   )
{
   assert(oracle != NULL);
   assert(considx >= 0);
   assert(considx < oracle->nconss);

   return oracle->conss[considx]->lhs;
}

/** gives right-hand side of a constraint */
SCIP_Real SCIPnlpiOracleGetConstraintRhs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   )
{
   assert(oracle != NULL);
   assert(considx >= 0);
   assert(considx < oracle->nconss);

   return oracle->conss[considx]->rhs;
}

/** gives name of a constraint, may be NULL */
char* SCIPnlpiOracleGetConstraintName(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   )
{
   assert(oracle != NULL);
   assert(considx >= 0);
   assert(considx < oracle->nconss);

   return oracle->conss[considx]->name;
}

/** gives maximum degree of a constraint or objective
 *  The degree is the maximal degree of all summands,, and is infinity for nonpolynomial terms.
 */ 
int SCIPnlpiOracleGetConstraintDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< index of constraint for which the degree is requested, or -1 for objective */
   )
{
   SCIP_NLPIORACLECONS* cons;

   assert(oracle != NULL);
   assert(considx >= -1);
   assert(considx < oracle->nconss);

   cons = considx < 0 ? oracle->objective : oracle->conss[considx];

   /* could do something more clever like exprtreeGetMaxDegree, but no solver uses this so far */
   if( cons->exprtree != NULL )
      return INT_MAX;

   if( cons->nquadelems > 0 )
      return 2;

   if( cons->nlinidxs > 0 )
      return 1;

   return 0;
}

/** Gives maximum degree over all constraints and the objective (or over all variables, resp.).
 * Thus, if this function returns 0, then the objective and all constraints are constant.
 * If it returns 1, then the problem in linear.
 * If it returns 2, then its a QP, QCP, or QCQP.
 * And if it returns > 2, then it is an NLP.
 */
int SCIPnlpiOracleGetMaxDegree(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   int i;
   int maxdegree;

   assert(oracle != NULL);

   SCIPdebugMessage("%p get max degree\n", (void*)oracle);

   updateVariableDegrees(oracle);

   maxdegree = 0;
   for( i = 0; i < oracle->nvars; ++i )
      if( oracle->vardegrees[i] > maxdegree )
      {
         maxdegree = oracle->vardegrees[i];
         if( maxdegree == INT_MAX )
            break;
      }

   return maxdegree;
}

/** Gives the evaluation capabilities that are shared among all expression trees in the problem. */
SCIP_EXPRINTCAPABILITY SCIPnlpiOracleGetEvalCapability(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   int c;
   SCIP_EXPRINTCAPABILITY evalcapability;

   assert(oracle != NULL);

   if( oracle->objective->exprtree != NULL )
      evalcapability = SCIPexprintGetExprtreeCapability(oracle->exprinterpreter, oracle->objective->exprtree);
   else
      evalcapability = SCIP_EXPRINTCAPABILITY_ALL;

   for( c = 0; c < oracle->nconss; ++c )
      if( oracle->conss[c]->exprtree != NULL )
         evalcapability &= SCIPexprintGetExprtreeCapability(oracle->exprinterpreter, oracle->conss[c]->exprtree);

   return evalcapability;
}

/** evaluates the objective function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveValue(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            objval              /**< pointer to store objective value */  
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p eval obj value\n", (void*)oracle);

   SCIP_CALL_QUIET( evalFunctionValue(oracle, oracle->objective, x, objval) );

   assert(oracle->objective->lhs == oracle->objective->rhs);  /*lint !e777*/
   *objval += oracle->objective->lhs;

   return SCIP_OKAY;
}

/** evaluates one constraint function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValue(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint to evaluate */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            conval              /**< pointer to store constraint value */  
   )
{
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);

   SCIPdebugMessage("%p eval cons value\n", (void*)oracle);

   SCIP_CALL_QUIET( evalFunctionValue(oracle, oracle->conss[considx], x, conval) );

   return SCIP_OKAY;
}

/** evaluates all constraint functions in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValues(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            convals             /**< buffer to store constraint values */  
   )
{
   int i;

   SCIPdebugMessage("%p eval cons values\n", (void*)oracle);

   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(convals != NULL);

   for( i = 0; i < oracle->nconss; ++i )
   {
      SCIP_CALL_QUIET( evalFunctionValue(oracle, oracle->conss[i], x, &convals[i]) );
   }

   return SCIP_OKAY;
}

/** computes the objective gradient in a given point
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveGradient(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Real*            objgrad             /**< pointer to store (dense) objective gradient */  
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p eval obj grad\n", (void*)oracle);

   SCIP_CALL_QUIET( evalFunctionGradient(oracle, oracle->objective, x, isnewx, objval, objgrad) );

   assert(oracle->objective->lhs == oracle->objective->rhs);  /*lint !e777*/
   *objval += oracle->objective->lhs;

   return SCIP_OKAY;
}

/** computes a constraints gradient in a given point
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to compute gradient for */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            conval,             /**< pointer to store constraint value */
   SCIP_Real*            congrad             /**< pointer to store (dense) constraint gradient */  
   )
{
   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);

   SCIPdebugMessage("%p eval cons grad\n", (void*)oracle);

   SCIP_CALL_QUIET( evalFunctionGradient(oracle, oracle->conss[considx], x, isnewx, conval, congrad) );

   return SCIP_OKAY;
}

/** gets sparsity pattern (rowwise) of Jacobian matrix
 * 
 *  Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 *  Adding or deleting constraints destroys the sparsity structure and make another call to this function necessary. 
 */
SCIP_RETCODE SCIPnlpiOracleGetJacobianSparsity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, offset[nconss] gives length of col, can be NULL */
   )
{
   SCIP_NLPIORACLECONS* cons;
   SCIP_Bool* nzflag;
   int nnz;
   int maxnnz;
   int i;
   int j;

   assert(oracle != NULL);

   SCIPdebugMessage("%p get jacobian sparsity\n", (void*)oracle);

   if( oracle->jacoffsets != NULL )
   {
      assert(oracle->jaccols != NULL);
      if( offset != NULL )
         *offset = oracle->jacoffsets;
      if( col != NULL )
         *col = oracle->jaccols;
      return SCIP_OKAY;
   }

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->jacoffsets, oracle->nconss + 1) );

   maxnnz = MIN(oracle->nvars, 10) * oracle->nconss;  /* initial guess */
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, maxnnz) );

   if( maxnnz == 0 )
   {
      /* no variables */
      BMSclearMemoryArray(oracle->jacoffsets, oracle->nconss + 1);
      if( offset != NULL )
         *offset = oracle->jacoffsets;
      if( col != NULL )
         *col = oracle->jaccols;
      return SCIP_OKAY;
   }
   nnz = 0;

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &nzflag, oracle->nvars) );

   for( i = 0; i < oracle->nconss; ++i )
   {
      oracle->jacoffsets[i] = nnz;

      cons = oracle->conss[i];
      assert(cons != NULL);

      if( cons->nquadelems == 0 && cons->exprtree == NULL )
      {
         /* for a linear constraint, we can just copy the linear indices from the constraint into the sparsity pattern */
         if( cons->nlinidxs > 0 )
         {
            SCIP_CALL( ensureIntArraySize(oracle->blkmem, &oracle->jaccols, &maxnnz, nnz + cons->nlinidxs) );
            BMScopyMemoryArray(&oracle->jaccols[nnz], cons->linidxs, cons->nlinidxs);  /*lint !e866*/
            nnz += cons->nlinidxs;
         }
         continue;
      }
      else if( cons->nlinidxs == 0 && cons->nquadelems == 0 )
      {
         /* for a constraint with exprtree only, we can just copy the exprvaridxs from the constraint into the sparsity pattern */
         int nvars;

         assert(cons->exprtree != NULL); /* this had been the first case */

         nvars = SCIPexprtreeGetNVars(cons->exprtree);
         assert(cons->exprvaridxs != NULL || nvars == 0);

         if( nvars > 0 )
         {
            SCIP_CALL( ensureIntArraySize(oracle->blkmem, &oracle->jaccols, &maxnnz, nnz + nvars) );
            BMScopyMemoryArray(&oracle->jaccols[nnz], cons->exprvaridxs, nvars);  /*lint !e866*/
            nnz += nvars;
         }
         continue;
      }

      /* check which variables appear in constraint i
       * @todo this could be done faster for very sparse constraint by assembling all appearing variables, sorting, and removing duplicates
       */
      BMSclearMemoryArray(nzflag, oracle->nvars);  /*lint !e644*/

      for( j = 0; j < cons->nlinidxs; ++j )
         nzflag[cons->linidxs[j]] = TRUE;

      for( j = 0; j < cons->nquadelems; ++j )
      {
         nzflag[cons->quadelems[j].idx1] = TRUE;
         nzflag[cons->quadelems[j].idx2] = TRUE;
      }

      if( cons->exprvaridxs != NULL )
      {
         assert(cons->exprtree != NULL);
         for( j = SCIPexprtreeGetNVars(cons->exprtree)-1; j >= 0; --j )
            nzflag[cons->exprvaridxs[j]] = TRUE;
      }

      /* store variables indices in jaccols */
      for( j = 0; j < oracle->nvars; ++j )
      {
         if( nzflag[j] == FALSE )
            continue;

         SCIP_CALL( ensureIntArraySize(oracle->blkmem, &oracle->jaccols, &maxnnz, nnz + 1) );
         oracle->jaccols[nnz] = j;
         ++nnz;
      }
   }

   oracle->jacoffsets[oracle->nconss] = nnz;

   /* shrink jaccols array to nnz */
   if( nnz < maxnnz )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(oracle->blkmem, &oracle->jaccols, maxnnz, nnz) );
   }

   BMSfreeBlockMemoryArray(oracle->blkmem, &nzflag, oracle->nvars);

   if( offset != NULL )
      *offset = oracle->jacoffsets;
   if( col != NULL )
      *col = oracle->jaccols;

   return SCIP_OKAY;
}

/** evaluates the Jacobi matrix in a given point
 * 
 *  The values in the Jacobi matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetJacobianSparsity.
 *  The user need to call SCIPnlpiOracleGetJacobianSparsity at least ones before using this function. 
 *
 * @return SCIP_INVALIDDATA, if the Jacobian could not be evaluated (domain error, etc.)
 */
SCIP_RETCODE SCIPnlpiOracleEvalJacobian(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            convals,            /**< pointer to store constraint values, can be NULL */ 
   SCIP_Real*            jacobi              /**< pointer to store sparse jacobian values */  
   )
{
   SCIP_NLPIORACLECONS* cons;
   SCIP_RETCODE retcode;
   SCIP_Real* grad;
   SCIP_Real* xx;
   SCIP_Real nlval;
   int i;
   int j;
   int k;
   int l;

   SCIPdebugMessage("%p eval jacobian\n", (void*)oracle);

   assert(oracle != NULL);
   assert(jacobi != NULL);

   assert(oracle->jacoffsets != NULL);
   assert(oracle->jaccols    != NULL);

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &grad, oracle->nvars) );
   xx = NULL;

   retcode = SCIP_OKAY;

   j = oracle->jacoffsets[0];
   k = 0;
   for( i = 0; i < oracle->nconss; ++i )
   {
      cons = oracle->conss[i];
      assert(cons != NULL);

      if( cons->nquadelems == 0 && cons->exprtree == NULL )
      {
         /* for a linear constraint, we can just copy the linear coefs from the constraint into the jacobian */
         if( cons->nlinidxs > 0 )
         {
            BMScopyMemoryArray(&jacobi[k], cons->lincoefs, cons->nlinidxs);  /*lint !e866*/
            j += cons->nlinidxs;
            k += cons->nlinidxs;
         }
         assert(j == oracle->jacoffsets[i+1]);
         continue;
      }

      if( cons->nlinidxs == 0 && cons->nquadelems == 0 )
      {
         /* for a constraint with exprtree only, we can just copy gradient of the exprtree from the constraint into jacobian */
         int nvars;

         assert(cons->exprtree != NULL); /* this had been the first case */

         nvars = SCIPexprtreeGetNVars(cons->exprtree);
         assert(nvars <= oracle->nvars);
         assert(cons->exprvaridxs != NULL || nvars == 0);

         if( nvars > 0 )
         {
            if( isnewx )
            {
               if( xx == NULL )
               {
                  /* if no xx yet, alloc it; make it large enough in case we need it again */
                  SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &xx, oracle->nvars) );
               }
               for( l = 0; l < nvars; ++l )
                  xx[l] = x[cons->exprvaridxs[l]];  /*lint !e613*/
            }

            SCIPdebugMessage("eval gradient of ");
            SCIPdebug( if( isnewx ) {printf("\nx ="); for( l = 0; l < nvars; ++l) printf(" %g", xx[l]); /*lint !e613*/ printf("\n");} )

            SCIP_CALL( SCIPexprintGrad(oracle->exprinterpreter, cons->exprtree, xx, isnewx, &nlval, grad) );  /*lint !e644*/

            SCIPdebug( printf("g ="); for( l = 0; l < nvars; ++l) printf(" %g", grad[l]); printf("\n"); )

            if( nlval != nlval || ABS(nlval) >= oracle->infinity )  /*lint !e777*/
            {
               SCIPdebugMessage("gradient evaluation yield invalid function value %g\n", nlval);
               retcode = SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
               break;
            }
            else
            {
               if( convals != NULL )
                  convals[i] = nlval;
               for( l = 0; l < nvars; ++l )
               {
                  assert(oracle->jaccols[j+l] == cons->exprvaridxs[l]);
                  if( !SCIPisFinite(grad[l]) )  /*lint !e777*/
                  {
                     SCIPdebugMessage("gradient evaluation yield invalid gradient value %g\n", grad[l]);
                     retcode = SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
                     break;
                  }
                  else
                     jacobi[k++] = grad[l];
               }
               if( l < nvars )
                  break;
               j += nvars;
            }
         }
         else if( convals != NULL )
         {
            SCIPdebugMessage("eval value of constant ");

            SCIP_CALL( SCIPexprintEval(oracle->exprinterpreter, cons->exprtree, NULL, &convals[i]) );
         }
         continue;
      }

      /* do dense eval @todo could do it sparse */
      retcode = SCIPnlpiOracleEvalConstraintGradient(oracle, i, x, isnewx, (convals ? &convals[i] : &nlval), grad);
      if( retcode != SCIP_OKAY )
         break;

      while( j < oracle->jacoffsets[i+1] )
         jacobi[k++] = grad[oracle->jaccols[j++]];
   }

   BMSfreeBlockMemoryArrayNull(oracle->blkmem, &xx, oracle->nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &grad, oracle->nvars);

   return retcode;
}

/** gets sparsity pattern of the Hessian matrix of the Lagrangian
 * 
 *  Note that internal data is returned in *offset and *col, thus the user must not allocate memory there.
 *  Adding or deleting variables, objective, or constraints may destroy the sparsity structure and make another call to this function necessary.
 *  Only elements of the lower left triangle and the diagonal are counted.
 */
SCIP_RETCODE SCIPnlpiOracleGetHessianLagSparsity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, offset[nconss] gives length of col, can be NULL */
   )
{
   int** colnz;   /* nonzeros in Hessian corresponding to one column */
   int*  collen;  /* collen[i] is length of array colnz[i] */
   int*  colnnz;  /* colnnz[i] is number of entries in colnz[i] (<= collen[i]) */
   int   nnz;
   int   i;
   int   j;
   int   cnt;

   assert(oracle != NULL);

   SCIPdebugMessage("%p get hessian lag sparsity\n", (void*)oracle);

   if( oracle->heslagoffsets != NULL )
   {
      assert(oracle->heslagcols != NULL);
      if( offset != NULL )
         *offset = oracle->heslagoffsets;
      if( col != NULL )
         *col = oracle->heslagcols;
      return SCIP_OKAY;
   }

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->heslagoffsets, oracle->nvars + 1) );

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &colnz,  oracle->nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &collen, oracle->nvars) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &colnnz, oracle->nvars) );
   BMSclearMemoryArray(colnz,  oracle->nvars);  /*lint !e644*/
   BMSclearMemoryArray(collen, oracle->nvars);  /*lint !e644*/
   BMSclearMemoryArray(colnnz, oracle->nvars);  /*lint !e644*/
   nnz = 0;

   if( oracle->objective->nquadelems != 0 )
   {
      SCIP_CALL( hessLagSparsitySetNzFlagForQuad(oracle, colnz, collen, colnnz, &nnz, oracle->objective->nquadelems, oracle->objective->quadelems) );
   }

   if( oracle->objective->exprtree != NULL )
   {
      SCIP_CALL( hessLagSparsitySetNzFlagForExprtree(oracle, colnz, collen, colnnz, &nnz, oracle->objective->exprvaridxs, oracle->objective->exprtree, oracle->nvars) );
   }

   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->conss[i]->nquadelems != 0 )
      {
         SCIP_CALL( hessLagSparsitySetNzFlagForQuad(oracle, colnz, collen, colnnz, &nnz, oracle->conss[i]->nquadelems, oracle->conss[i]->quadelems) );
      }

      if( oracle->conss[i]->exprtree != NULL )
      {
         SCIP_CALL( hessLagSparsitySetNzFlagForExprtree(oracle, colnz, collen, colnnz, &nnz, oracle->conss[i]->exprvaridxs, oracle->conss[i]->exprtree, oracle->nvars) );
      }
   }

   SCIP_ALLOC( BMSallocBlockMemoryArray(oracle->blkmem, &oracle->heslagcols, nnz) );

   /* set hessian sparsity from colnz, colnnz */
   cnt = 0;
   for( i = 0; i < oracle->nvars; ++i )
   {
      oracle->heslagoffsets[i] = cnt;
      for( j = 0; j < colnnz[i]; ++j )
      {
         assert(cnt < nnz);
         oracle->heslagcols[cnt++] = colnz[i][j];
      }
      BMSfreeBlockMemoryArrayNull(oracle->blkmem, &colnz[i], collen[i]);  /*lint !e866*/
      collen[i] = 0;
   }
   oracle->heslagoffsets[oracle->nvars] = cnt;
   assert(cnt == nnz);

   BMSfreeBlockMemoryArray(oracle->blkmem, &colnz,  oracle->nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &colnnz, oracle->nvars);
   BMSfreeBlockMemoryArray(oracle->blkmem, &collen, oracle->nvars);

   if( offset != NULL )
      *offset = oracle->heslagoffsets;
   if( col != NULL )
      *col = oracle->heslagcols;

   return SCIP_OKAY;
}

/** evaluates the Hessian matrix of the Lagrangian in a given point
 * 
 *  The values in the Hessian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetHessianLagSparsity.
 *  The user must call SCIPnlpiOracleGetHessianLagSparsity at least ones before using this function. 
 *  Only elements of the lower left triangle and the diagonal are computed.
 *
 * @return SCIP_INVALIDDATA, if the Hessian could not be evaluated (domain error, etc.)
 */
SCIP_RETCODE SCIPnlpiOracleEvalHessianLag(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real             objfactor,          /**< weight for objective function */
   const SCIP_Real*      lambda,             /**< weights (Lagrangian multipliers) for the constraints */ 
   SCIP_Real*            hessian             /**< pointer to store sparse hessian values */  
   )
{  /*lint --e{715}*/
   int i;

   assert(oracle != NULL);
   assert(x != NULL);
   assert(lambda != NULL || oracle->nconss == 0);
   assert(hessian != NULL);

   assert(oracle->heslagoffsets != NULL);
   assert(oracle->heslagcols != NULL);

   SCIPdebugMessage("%p eval hessian lag\n", (void*)oracle);

   for( i = oracle->heslagoffsets[oracle->nvars] - 1; i >= 0; --i )
      hessian[i] = 0.0;

   if( objfactor != 0.0 )
   {
      SCIP_CALL( hessLagAddQuad(objfactor, oracle->objective->nquadelems, oracle->objective->quadelems, oracle->heslagoffsets, oracle->heslagcols, hessian) );
      SCIP_CALL_QUIET( hessLagAddExprtree(oracle, objfactor, x, isnewx, oracle->objective->exprvaridxs, oracle->objective->exprtree, oracle->heslagoffsets, oracle->heslagcols, hessian) );
   }

   for( i = 0; i < oracle->nconss; ++i )
   {
      assert( lambda != NULL ); /* for lint */
      if( lambda[i] == 0.0 )
         continue;
      SCIP_CALL( hessLagAddQuad(lambda[i], oracle->conss[i]->nquadelems, oracle->conss[i]->quadelems, oracle->heslagoffsets, oracle->heslagcols, hessian) );
      SCIP_CALL_QUIET( hessLagAddExprtree(oracle, lambda[i], x, isnewx, oracle->conss[i]->exprvaridxs, oracle->conss[i]->exprtree, oracle->heslagoffsets, oracle->heslagcols, hessian) );
   }

   return SCIP_OKAY;
}

/** prints the problem to a file. */
SCIP_RETCODE SCIPnlpiOraclePrintProblem(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   )
{  /*lint --e{777} */
   int i;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(oracle != NULL);

   SCIPdebugMessage("%p print problem\n", (void*)oracle);

   if( file == NULL )
      file = stdout;

   SCIPmessageFPrintInfo(messagehdlr, file, "NLPI Oracle %s: %d variables and %d constraints\n", oracle->name ? oracle->name : "", oracle->nvars, oracle->nconss);
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
         SCIPmessageFPrintInfo(messagehdlr, file, "%10s", oracle->varnames[i]);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "x%09d", i);
      SCIPmessageFPrintInfo(messagehdlr, file, ": [%8g, %8g]", oracle->varlbs[i], oracle->varubs[i]);
      if( oracle->vardegreesuptodate )
         SCIPmessageFPrintInfo(messagehdlr, file, "\t degree: %d", oracle->vardegrees[i]);
      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }

   SCIPmessageFPrintInfo(messagehdlr, file, "objective: ");
   SCIP_CALL( printFunction(oracle, messagehdlr, file, oracle->objective, FALSE, FALSE) );
   if( oracle->objective->lhs != 0.0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.20g", oracle->objective->lhs);
   SCIPmessageFPrintInfo(messagehdlr, file, "\n");

   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->conss[i]->name != NULL )
         SCIPmessageFPrintInfo(messagehdlr, file, "%10s", oracle->conss[i]->name);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "con%07d", i);

      lhs = oracle->conss[i]->lhs;
      rhs = oracle->conss[i]->rhs;
      SCIPmessageFPrintInfo(messagehdlr, file, ": ");
      if( lhs > -oracle->infinity && rhs < oracle->infinity && lhs != rhs )
         SCIPmessageFPrintInfo(messagehdlr, file, "%.20g <= ", lhs);

      SCIP_CALL( printFunction(oracle, messagehdlr, file, oracle->conss[i], FALSE, FALSE) );

      if( lhs == rhs )
         SCIPmessageFPrintInfo(messagehdlr, file, " = %.20g", rhs);
      else if( rhs <  oracle->infinity )
         SCIPmessageFPrintInfo(messagehdlr, file, " <= %.20g", rhs);
      else if( lhs > -oracle->infinity )
         SCIPmessageFPrintInfo(messagehdlr, file, " >= %.20g", lhs);

      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }

   return SCIP_OKAY;
}

/** prints the problem to a file in GAMS format
 * If there are variable (equation, resp.) names with more than 9 characters, then variable (equation, resp.) names are prefixed with an unique identifier.
 * This is to make it easier to identify variables solution output in the listing file.
 * Names with more than 64 characters are shorten to 64 letters due to GAMS limits.
 */
SCIP_RETCODE SCIPnlpiOraclePrintProblemGams(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real*            initval,            /**< starting point values for variables or NULL */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   )
{  /*lint --e{777} */
   int i;
   int nllevel; /* level of nonlinearity of problem: linear = 0, quadratic, smooth nonlinear, nonsmooth */
   static const char* nllevelname[4] = { "LP", "QCP", "NLP", "DNLP" };
   char problemname[SCIP_MAXSTRLEN];
   char namebuf[70];
   SCIP_Bool havelongvarnames;
   SCIP_Bool havelongequnames;

   SCIPdebugMessage("%p print problem gams\n", (void*)oracle);

   assert(oracle != NULL);

   if( file == NULL )
      file = stdout;

   nllevel = 0;

   havelongvarnames = FALSE;
   for( i = 0; i < oracle->nvars; ++i )
      if( oracle->varnames != NULL && oracle->varnames[i] != NULL && strlen(oracle->varnames[i]) > 9 )
      {
         havelongvarnames = TRUE;
         break;
      }

   havelongequnames = FALSE;
   for( i = 0; i < oracle->nconss; ++i )
      if( oracle->conss[i]->name && strlen(oracle->conss[i]->name) > 9 )
      {
         havelongequnames = TRUE;
         break;
      }

   SCIPmessageFPrintInfo(messagehdlr, file, "$offlisting\n");
   SCIPmessageFPrintInfo(messagehdlr, file, "$offdigit\n");
   SCIPmessageFPrintInfo(messagehdlr, file, "* NLPI Oracle Problem %s\n", oracle->name ? oracle->name : "");
   SCIPmessageFPrintInfo(messagehdlr, file, "Variables ");
   for( i = 0; i < oracle->nvars; ++i )
   {
      printName(namebuf, oracle->varnames != NULL ? oracle->varnames[i] : NULL, i, 'x', NULL, havelongvarnames);
      SCIPmessageFPrintInfo(messagehdlr, file, "%s, ", namebuf);
      if( i % 10 == 9 )
         SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }
   SCIPmessageFPrintInfo(messagehdlr, file, "NLPIORACLEOBJVAR;\n\n");
   for( i = 0; i < oracle->nvars; ++i )
   {
      char* name;
      name = oracle->varnames != NULL ? oracle->varnames[i] : NULL;
      if( oracle->varlbs[i] == oracle->varubs[i] )
      {
         printName(namebuf, name, i, 'x', NULL, havelongvarnames);
         SCIPmessageFPrintInfo(messagehdlr, file, "%s.fx = %.20g;\t", namebuf, oracle->varlbs[i]);
      }
      else
      {
         if( oracle->varlbs[i] > -oracle->infinity )
         {
            printName(namebuf, name, i, 'x', NULL, havelongvarnames);
            SCIPmessageFPrintInfo(messagehdlr, file, "%s.lo = %.20g;\t", namebuf, oracle->varlbs[i]);
         }
         if( oracle->varubs[i] <  oracle->infinity )
         {
            printName(namebuf, name, i, 'x', NULL, havelongvarnames);
            SCIPmessageFPrintInfo(messagehdlr, file, "%s.up = %.20g;\t", namebuf, oracle->varubs[i]);
         }
      }
      if( initval != NULL )
      {
         printName(namebuf, name, i, 'x', NULL, havelongvarnames);
         SCIPmessageFPrintInfo(messagehdlr, file, "%s.l = %.20g;\t", namebuf, initval[i]);
      }
      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }
   SCIPmessageFPrintInfo(messagehdlr, file, "\n");

   SCIPmessageFPrintInfo(messagehdlr, file, "Equations ");
   for( i = 0; i < oracle->nconss; ++i )
   {
      printName(namebuf, oracle->conss[i]->name, i, 'e', NULL, havelongequnames);
      SCIPmessageFPrintInfo(messagehdlr, file, "%s, ", namebuf);

      if( oracle->conss[i]->lhs > -oracle->infinity && oracle->conss[i]->rhs < oracle->infinity && oracle->conss[i]->lhs != oracle->conss[i]->rhs )
      {
         /* ranged row: add second constraint */
         printName(namebuf, oracle->conss[i]->name, i, 'e', "_RNG", havelongequnames);
         SCIPmessageFPrintInfo(messagehdlr, file, "%s, ", namebuf);
      }
      if( i % 10 == 9 )
         SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }
   SCIPmessageFPrintInfo(messagehdlr, file, "NLPIORACLEOBJ;\n\n");

   SCIPmessageFPrintInfo(messagehdlr, file, "NLPIORACLEOBJ.. NLPIORACLEOBJVAR =E= ");
   SCIP_CALL( printFunction(oracle, messagehdlr, file, oracle->objective, havelongvarnames, havelongequnames) );
   if( oracle->objective->lhs != 0.0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.20g", oracle->objective->lhs);
   SCIPmessageFPrintInfo(messagehdlr, file, ";\n");

   for( i = 0; i < oracle->nconss; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;

      printName(namebuf, oracle->conss[i]->name, i, 'e', NULL, havelongequnames);
      SCIPmessageFPrintInfo(messagehdlr, file, "%s.. ", namebuf);

      SCIP_CALL( printFunction(oracle, messagehdlr, file, oracle->conss[i], havelongvarnames, havelongequnames) );

      lhs = oracle->conss[i]->lhs;
      rhs = oracle->conss[i]->rhs;

      if( lhs == rhs )
         SCIPmessageFPrintInfo(messagehdlr, file, " =E= %.20g", rhs);
      else if( rhs <  oracle->infinity )
         SCIPmessageFPrintInfo(messagehdlr, file, " =L= %.20g", rhs);
      else if( lhs > -oracle->infinity )
         SCIPmessageFPrintInfo(messagehdlr, file, " =G= %.20g", lhs);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, " =N= 0");
      SCIPmessageFPrintInfo(messagehdlr, file, ";\n");

      if( lhs > -oracle->infinity && rhs < oracle->infinity && lhs != rhs )
      {
         printName(namebuf, oracle->conss[i]->name, i, 'e', "_RNG", havelongequnames);
         SCIPmessageFPrintInfo(messagehdlr, file, "%s.. ", namebuf);

         SCIP_CALL( printFunction(oracle, messagehdlr, file, oracle->conss[i], havelongvarnames, havelongequnames) );

         SCIPmessageFPrintInfo(messagehdlr, file, " =G= %.20g;\n", lhs);
      }

      if( nllevel <= 0 && oracle->conss[i]->nquadelems > 0 )
         nllevel = 1;
      if( nllevel <= 1 && oracle->conss[i]->exprtree != NULL )
         nllevel = 2;
      if( nllevel <= 2 && oracle->conss[i]->exprtree != NULL && exprIsNonSmooth(SCIPexprtreeGetRoot(oracle->conss[i]->exprtree)) )
         nllevel = 3;
   }

   (void) SCIPsnprintf(problemname, SCIP_MAXSTRLEN, "%s", oracle->name ? oracle->name : "m");

   SCIPmessageFPrintInfo(messagehdlr, file, "Model %s / all /;\n", problemname);
   SCIPmessageFPrintInfo(messagehdlr, file, "option limrow = 0;\n");
   SCIPmessageFPrintInfo(messagehdlr, file, "option limcol = 0;\n");
   SCIPmessageFPrintInfo(messagehdlr, file, "Solve %s minimizing NLPIORACLEOBJVAR using %s;\n", problemname, nllevelname[nllevel]);

   return SCIP_OKAY;
}

/**@} */
