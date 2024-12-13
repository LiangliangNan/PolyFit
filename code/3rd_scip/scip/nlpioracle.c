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

/**@file    nlpioracle.c
 * @ingroup OTHER_CFILES
 * @brief   implementation of NLPI oracle
 * @author  Stefan Vigerske
 *
 * @todo jacobi evaluation should be sparse
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/nlpioracle.h"
#include "scip/exprinterpret.h"
#include "scip/expr_pow.h"
#include "scip/expr_varidx.h"

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

   SCIP_EXPR*            expr;               /**< expression for nonlinear part, or NULL if none */
   SCIP_EXPRINTDATA*     exprintdata;        /**< expression interpret data for expression, or NULL if no expr or not compiled yet */

   char*                 name;               /**< name of constraint */
};
typedef struct SCIP_NlpiOracleCons SCIP_NLPIORACLECONS;

struct SCIP_NlpiOracle
{
   char*                 name;               /**< name of problem */

   int                   varssize;           /**< length of variables related arrays */
   int                   nvars;              /**< number of variables */
   SCIP_Real*            varlbs;             /**< array with variable lower bounds */
   SCIP_Real*            varubs;             /**< array with variable upper bounds */
   char**                varnames;           /**< array with variable names */
   int*                  varlincount;        /**< array with number of appearances of variable in linear part of objective or constraints */
   int*                  varnlcount;         /**< array with number of appearances of variable in nonlinear part of objective or constraints */

   int                   consssize;          /**< length of constraints related arrays */
   int                   nconss;             /**< number of constraints */
   SCIP_NLPIORACLECONS** conss;              /**< constraints, or NULL if none */

   SCIP_NLPIORACLECONS*  objective;          /**< objective */

   int*                  jacoffsets;         /**< rowwise jacobi sparsity pattern: constraint offsets in jaccols */
   int*                  jaccols;            /**< rowwise jacobi sparsity pattern: indices of variables appearing in constraints */

   int*                  heslagoffsets;      /**< rowwise sparsity pattern of hessian matrix of Lagrangian: row offsets in heslagcol */
   int*                  heslagcols;         /**< rowwise sparsity pattern of hessian matrix of Lagrangian: column indices; sorted for each row */

   SCIP_EXPRINT*         exprinterpreter;    /**< interpreter for expressions: evaluation and derivatives */
   SCIP_CLOCK*           evalclock;          /**< clock measuring evaluation time */
};

/**@} */

/*lint -e440*/
/*lint -e441*/
/*lint -e866*/

/**@name Local functions */
/**@{ */

/** ensures that those arrays in oracle that store information on variables have at least a given length */
static
SCIP_RETCODE ensureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< NLPIORACLE data structure */
   int                   minsize             /**< minimal required size */
   )
{
   assert(oracle != NULL);

   if( minsize > oracle->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, minsize);
      assert(newsize >= minsize);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &oracle->varlbs, oracle->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &oracle->varubs, oracle->varssize, newsize) );
      if( oracle->varnames != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &oracle->varnames, oracle->varssize, newsize) );
      }
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &oracle->varlincount, oracle->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &oracle->varnlcount, oracle->varssize, newsize) );

      oracle->varssize = newsize;
   }
   assert(oracle->varssize >= minsize);

   return SCIP_OKAY;
}

/** ensures that constraints array in oracle has at least a given length */
static
SCIP_RETCODE ensureConssSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< NLPIORACLE data structure */
   int                   minsize             /**< minimal required size */
   )
{
   assert(oracle != NULL);

   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &oracle->conss, &oracle->consssize, minsize) );
   assert(oracle->consssize >= minsize);

   return SCIP_OKAY;
}

/** ensures that arrays for linear part in a oracle constraints have at least a given length */
static
SCIP_RETCODE ensureConsLinSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLECONS*  cons,               /**< oracle constraint */
   int                   minsize             /**< minimal required size */
   )
{
   assert(cons != NULL);

   if( minsize > cons->linsize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, minsize);
      assert(newsize >= minsize);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &cons->linidxs,  cons->linsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &cons->lincoefs, cons->linsize, newsize) );
      cons->linsize = newsize;
   }
   assert(cons->linsize >= minsize);

   return SCIP_OKAY;
}

/** ensures that a given array of integers has at least a given length */
static
SCIP_RETCODE ensureIntArraySize(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 intarray,           /**< array of integers */
   int*                  len,                /**< length of array (modified if reallocated) */
   int                   minsize             /**< minimal required array length */
   )
{
   assert(intarray != NULL);
   assert(len != NULL);

   SCIP_CALL( SCIPensureBlockMemoryArray(scip, intarray, len, minsize) );
   assert(*len >= minsize);

   return SCIP_OKAY;
}

/** Invalidates the sparsity pattern of the Jacobian.
 *  Should be called when constraints are added or deleted.
 */
static
void invalidateJacobiSparsity(
   SCIP*                 scip,               /**< SCIP data structure */
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
   SCIPfreeBlockMemoryArray(scip, &oracle->jaccols,    oracle->jacoffsets[oracle->nconss]);
   SCIPfreeBlockMemoryArray(scip, &oracle->jacoffsets, oracle->nconss + 1);
}

/** Invalidates the sparsity pattern of the Hessian of the Lagragian.
 *  Should be called when the objective is set or constraints are added or deleted.
 */
static
void invalidateHessianLagSparsity(
   SCIP*                 scip,               /**< SCIP data structure */
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
   SCIPfreeBlockMemoryArray(scip, &oracle->heslagcols,    oracle->heslagoffsets[oracle->nvars]);
   SCIPfreeBlockMemoryArray(scip, &oracle->heslagoffsets, oracle->nvars + 1);
}

/** increases or decreases variable counts in oracle w.r.t. linear and nonlinear appearance */
static
SCIP_RETCODE updateVariableCounts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< oracle data structure */
   int                   factor,             /**< whether to add (factor=1) or remove (factor=-1) variable counts */
   int                   nlinidxs,           /**< number of linear indices */
   const int*            linidxs,            /**< indices of variables in linear part */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int j;

   assert(oracle != NULL);
   assert(oracle->varlincount != NULL || (nlinidxs == 0 && expr == NULL));
   assert(oracle->varnlcount != NULL || (nlinidxs == 0 && expr == NULL));
   assert(factor == 1 || factor == -1);
   assert(nlinidxs == 0 || linidxs != NULL);

   for( j = 0; j < nlinidxs; ++j )
   {
      oracle->varlincount[linidxs[j]] += factor;
      assert(oracle->varlincount[linidxs[j]] >= 0);
   }

   if( expr != NULL )
   {
      SCIP_EXPRITER* it;

      SCIP_CALL( SCIPcreateExpriter(scip, &it) );
      SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );

      for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
         if( SCIPisExprVaridx(scip, expr) )
         {
            oracle->varnlcount[SCIPgetIndexExprVaridx(expr)] += factor;
            assert(oracle->varnlcount[SCIPgetIndexExprVaridx(expr)] >= 0);
         }

      SCIPfreeExpriter(&it);
   }

   return SCIP_OKAY;
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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_NLPIORACLECONS** cons,               /**< buffer where to store pointer to constraint */
   int                   nlinidxs,           /**< length of linear part */
   const int*            linidxs,            /**< indices of linear part, or NULL if nlinidxs == 0 */
   const SCIP_Real*      lincoefs,           /**< coefficients of linear part, or NULL if nlinidxs == 0 */
   SCIP_EXPR*            expr,               /**< expression, or NULL */
   SCIP_Real             lhs,                /**< left-hand-side of constraint */
   SCIP_Real             rhs,                /**< right-hand-side of constraint */
   const char*           name                /**< name of constraint, or NULL */
   )
{
   assert(cons != NULL);
   assert(nlinidxs >= 0);
   assert(linidxs != NULL  || nlinidxs == 0);
   assert(lincoefs != NULL || nlinidxs == 0);
   assert(EPSLE(lhs, rhs, SCIP_DEFAULT_EPSILON));

   SCIP_CALL( SCIPallocClearBlockMemory(scip, cons) );
   assert(*cons != NULL);

   if( nlinidxs > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*cons)->linidxs,  linidxs,  nlinidxs) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*cons)->lincoefs, lincoefs, nlinidxs) );
      (*cons)->linsize  = nlinidxs;
      (*cons)->nlinidxs = nlinidxs;

      /* sort, merge duplicates, remove zero's */
      sortLinearCoefficients(&(*cons)->nlinidxs, (*cons)->linidxs, (*cons)->lincoefs);
      assert((*cons)->linidxs[0] >= 0);
   }

   if( expr != NULL )
   {
      (*cons)->expr = expr;
      SCIPcaptureExpr(expr);

      SCIP_CALL( SCIPexprintCompile(scip, oracle->exprinterpreter, (*cons)->expr, &(*cons)->exprintdata) );
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
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*cons)->name, name, strlen(name)+1) );
   }

   /* add variable counts */
   SCIP_CALL( updateVariableCounts(scip, oracle, 1, (*cons)->nlinidxs, (*cons)->linidxs, (*cons)->expr) );

   return SCIP_OKAY;
}

/** frees a constraint */
static
SCIP_RETCODE freeConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_NLPIORACLECONS** cons,               /**< pointer to constraint that should be freed */
   SCIP_Bool             updatevarcount      /**< whether the update variable counts (typically TRUE) */
   )
{
   assert(oracle != NULL);
   assert(cons   != NULL);
   assert(*cons  != NULL);

   SCIPdebugMessage("free constraint %p\n", (void*)*cons);

   /* remove variable counts */
   if( updatevarcount )
   {
      SCIP_CALL( updateVariableCounts(scip, oracle, -1, (*cons)->nlinidxs, (*cons)->linidxs, (*cons)->expr) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*cons)->linidxs, (*cons)->linsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*cons)->lincoefs, (*cons)->linsize);

   if( (*cons)->expr != NULL )
   {
      SCIP_CALL( SCIPexprintFreeData(scip, oracle->exprinterpreter, (*cons)->expr, &(*cons)->exprintdata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &(*cons)->expr) );
   }

   if( (*cons)->name != NULL )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*cons)->name, strlen((*cons)->name)+1);
   }

   SCIPfreeBlockMemory(scip, cons);
   assert(*cons == NULL);

   return SCIP_OKAY;
}

/** frees all constraints
 *
 * \attention This omits updating the variable counts in the oracle.
 */
static
SCIP_RETCODE freeConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   int i;

   assert(oracle != NULL);

   SCIPdebugMessage("%p free constraints\n", (void*)oracle);

   for( i = 0; i < oracle->nconss; ++i )
   {
      SCIP_CALL( freeConstraint(scip, oracle, &oracle->conss[i], FALSE) );
      assert(oracle->conss[i] == NULL);
   }
   oracle->nconss = 0;

   SCIPfreeBlockMemoryArrayNull(scip, &oracle->conss, oracle->consssize);
   oracle->consssize = 0;

   return SCIP_OKAY;
}

/** moves one variable
 * The place where it moves to need to be empty (all NULL) but allocated.
 * Note that this function does not update the variable indices in the constraints!
 */
static
SCIP_RETCODE moveVariable(
   SCIP*                 scip,               /**< SCIP data structure */
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

   assert(oracle->varnames == NULL || oracle->varnames[toidx] == NULL);

   oracle->varlbs[toidx] = oracle->varlbs[fromidx];
   oracle->varubs[toidx] = oracle->varubs[fromidx];
   oracle->varlbs[fromidx] = -SCIPinfinity(scip);
   oracle->varubs[fromidx] =  SCIPinfinity(scip);

   oracle->varlincount[toidx] = oracle->varlincount[fromidx];
   oracle->varnlcount[toidx] = oracle->varnlcount[fromidx];
   oracle->varlincount[fromidx] = 0;
   oracle->varnlcount[fromidx] = 0;

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
   SCIP*                 scip,               /**< SCIP data structure */
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
            SCIPfreeBlockMemoryArray(scip, &oracle->varnames[i], strlen(oracle->varnames[i])+1);  /*lint !e866*/
         }
      }
      SCIPfreeBlockMemoryArrayNull(scip, &oracle->varnames, oracle->varssize);
   }
   oracle->nvars = 0;

   SCIPfreeBlockMemoryArrayNull(scip, &oracle->varlbs, oracle->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &oracle->varubs, oracle->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &oracle->varlincount, oracle->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &oracle->varnlcount, oracle->varssize);

   oracle->varssize = 0;
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
   int**                 linidxs,            /**< variable indices */
   SCIP_Real**           coefs,              /**< variable coefficients */
   int*                  nidxs               /**< number of indices */
   )
{
   int i;
   int offset;

   SCIPdebugMessage("clear deleted linear elements\n");

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

/** computes the value of a function */
static
SCIP_RETCODE evalFunctionValue(
   SCIP*                 scip,               /**< SCIP data structure */
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

   if( cons->expr != NULL )
   {
      SCIP_Real  nlval;

      SCIP_CALL( SCIPexprintEval(scip, oracle->exprinterpreter, cons->expr, cons->exprintdata, (SCIP_Real*)x, &nlval) );
      if( !SCIPisFinite(nlval) || SCIPisInfinity(scip, ABS(nlval)) )
         *val  = nlval;
      else
         *val += nlval;
   }

   return SCIP_OKAY;
}

/** computes the value and gradient of a function
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
static
SCIP_RETCODE evalFunctionGradient(
   SCIP*                 scip,               /**< SCIP data structure */
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

   if( cons->expr != NULL )
   {
      SCIP_Real nlval;
      int       i;

      SCIPdebugMsg(scip, "eval gradient of ");
      SCIPdebug( if( isnewx ) {printf("\nx ="); for( i = 0; i < oracle->nvars; ++i) printf(" %g", x[i]); printf("\n");} )

      SCIP_CALL( SCIPexprintGrad(scip, oracle->exprinterpreter, cons->expr, cons->exprintdata, (SCIP_Real*)x, isnewx, &nlval, grad) );

      SCIPdebug( printf("g ="); for( i = 0; i < oracle->nvars; ++i) printf(" %g", grad[i]); printf("\n"); )

      /* check for eval error */
      if( !SCIPisFinite(nlval) || SCIPisInfinity(scip, ABS(nlval)) )
      {
         SCIPdebugMessage("gradient evaluation yield invalid function value %g\n", nlval);
         return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
      }
      for( i = 0; i < oracle->nvars; ++i )
         if( !SCIPisFinite(grad[i]) )
         {
            SCIPdebugMessage("gradient evaluation yield invalid gradient value %g\n", grad[i]);
            return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
         }

      *val += nlval;
   }

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
         grad[*linidxs] += *lincoefs;
      }
   }

   return SCIP_OKAY;
}

/** collects indices of nonzero entries in the lower-left part of the hessian matrix of an expression
 * adds the indices to a given set of indices, avoiding duplicates */
static
SCIP_RETCODE hessLagSparsitySetNzFlagForExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< NLPI oracle */
   int**                 colnz,              /**< indices of nonzero entries for each column */
   int*                  collen,             /**< space allocated to store indices of nonzeros for each column */
   int*                  colnnz,             /**< number of nonzero entries for each column */
   int*                  nzcount,            /**< counter for total number of nonzeros; should be increased when nzflag is set to 1 the first time */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< expression interpreter data for expression */
   int                   dim                 /**< dimension of matrix */
   )
{
   SCIP_Real* x;
   int* rowidxs;
   int* colidxs;
   int nnz;
   int row;
   int col;
   int pos;
   int i;

   assert(oracle != NULL);
   assert(colnz  != NULL);
   assert(collen != NULL);
   assert(colnnz != NULL);
   assert(nzcount != NULL);
   assert(expr != NULL);
   assert(dim >= 0);

   SCIPdebugMessage("%p hess lag sparsity set nzflag for expr\n", (void*)oracle);

   SCIP_CALL( SCIPallocBufferArray(scip, &x, oracle->nvars) );
   for( i = 0; i < oracle->nvars; ++i )
      x[i] = 2.0; /* hope that this value does not make much trouble for the evaluation routines */

   SCIP_CALL( SCIPexprintHessianSparsity(scip, oracle->exprinterpreter, expr, exprintdata, x, &rowidxs, &colidxs, &nnz) );

   for( i = 0; i < nnz; ++i )
   {
      row = rowidxs[i];
      col = colidxs[i];

      assert(row < oracle->nvars);
      assert(col <= row);

      if( colnz[row] == NULL || !SCIPsortedvecFindInt(colnz[row], col, colnnz[row], &pos) )
      {
         SCIP_CALL( ensureIntArraySize(scip, &colnz[row], &collen[row], colnnz[row]+1) );
         SCIPsortedvecInsertInt(colnz[row], col, &colnnz[row], NULL);
         ++*nzcount;
      }
   }

   SCIPfreeBufferArray(scip, &x);

   return SCIP_OKAY;
}

/** adds hessian of an expression into hessian structure */
static
SCIP_RETCODE hessLagAddExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< oracle */
   SCIP_Real             weight,             /**< weight of quadratic part */
   const SCIP_Real*      x,                  /**< point for which hessian should be returned */
   SCIP_Bool             new_x,              /**< whether point has been evaluated before */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< expression interpreter data for expression */
   int*                  hesoffset,          /**< row offsets in sparse matrix that is to be filled */
   int*                  hescol,             /**< column indices in sparse matrix that is to be filled */
   SCIP_Real*            values              /**< buffer for values of sparse matrix that is to be filled */
   )
{
   SCIP_Real val;
   SCIP_Real* h;
   int* rowidxs;
   int* colidxs;
   int nnz;
   int row;
   int col;
   int pos;
   int i;

   SCIPdebugMessage("%p hess lag add expr\n", (void*)oracle);

   assert(oracle != NULL);
   assert(x != NULL || new_x == FALSE);
   assert(expr != NULL);
   assert(hesoffset != NULL);
   assert(hescol != NULL);
   assert(values != NULL);

   SCIP_CALL( SCIPexprintHessian(scip, oracle->exprinterpreter, expr, exprintdata, (SCIP_Real*)x, new_x, &val, &rowidxs, &colidxs, &h, &nnz) );
   if( !SCIPisFinite(val) )
   {
      SCIPdebugMessage("hessian evaluation yield invalid function value %g\n", val);
      return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
   }

   for( i = 0; i < nnz; ++i )
   {
      if( !SCIPisFinite(h[i]) )
      {
         SCIPdebugMessage("hessian evaluation yield invalid hessian value %g\n", *h);
         return SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
      }

      if( h[i] == 0.0 )
         continue;

      row = rowidxs[i];
      col = colidxs[i];

      if( !SCIPsortedvecFindInt(&hescol[hesoffset[row]], col, hesoffset[row+1] - hesoffset[row], &pos) )
      {
         SCIPerrorMessage("Could not find entry (%d, %d) in hessian sparsity\n", row, col);
         return SCIP_ERROR;
      }

      values[hesoffset[row] + pos] += weight * h[i];
   }

   return SCIP_OKAY;
}

/** prints a name, if available, makes sure it has not more than 64 characters, and adds a unique prefix if the longnames flag is set */
static
void printName(
   char*                 buffer,             /**< buffer to print to, has to be not NULL and should be at least 65 bytes */
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
         (void) SCIPsnprintf(buffer, 64, "%c%05d%.*s%s", prefix, idx, suffix != NULL ? (int)(57-strlen(suffix)) : 57, name, suffix ? suffix : "");
      else
         (void) SCIPsnprintf(buffer, 64, "%c%05d", prefix, idx);
   }
   else
   {
      if( name != NULL )
      {
         assert(strlen(name) + (suffix != NULL ? strlen(suffix) : 0) <= 64);
         (void) SCIPsnprintf(buffer, 64, "%s%s", name, suffix != NULL ? suffix : "");
      }
      else
      {
         assert(1 + 5 + (suffix != NULL ? strlen(suffix) : 0) <= 64);
         (void) SCIPsnprintf(buffer, 64, "%c%d%s", prefix, idx, suffix != NULL ? suffix : "");
      }
   }
}

/** prints a function */
static
SCIP_RETCODE printFunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   FILE*                 file,               /**< file to print to, has to be not NULL */
   SCIP_NLPIORACLECONS*  cons,               /**< constraint which function to print */
   SCIP_Bool             longvarnames        /**< whether variable names need to be shorten to 64 characters */
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
      SCIPinfoMessage(scip, file, "%+.15g*%s", cons->lincoefs[i], namebuf);
      if( i % 10 == 9 )
         SCIPinfoMessage(scip, file, "\n");
   }

   if( cons->expr != NULL )
   {
      /* TODO SCIPprintExpr does not use the variable names in oracle->varnames, probably that should be changed */
      SCIPinfoMessage(scip, file, " +");
      SCIP_CALL( SCIPprintExpr(scip, cons->expr, file) );
   }

   return SCIP_OKAY;
}

/** returns whether an expression contains nonsmooth operands (min, max, abs, ...) */
static
SCIP_RETCODE exprIsNonSmooth(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool*            nonsmooth           /**< buffer to store whether expression seems nonsmooth */
   )
{
   SCIP_EXPRITER* it;

   assert(expr != NULL);
   assert(nonsmooth != NULL);

   *nonsmooth = FALSE;

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );

   for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      const char* hdlrname;
      if( SCIPisExprSignpower(scip, expr) )
      {
         *nonsmooth = TRUE;
         break;
      }
      hdlrname = SCIPexprhdlrGetName(SCIPexprGetHdlr(expr));
      if( strcmp(hdlrname, "abs") == 0 )
      {
         *nonsmooth = TRUE;
         break;
      }
      if( strcmp(hdlrname, "min") == 0 )
      {
         *nonsmooth = TRUE;
         break;
      }
      if( strcmp(hdlrname, "max") == 0 )
      {
         *nonsmooth = TRUE;
         break;
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/**@} */

/**@name public function */
/**@{ */

/** creates an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE**     oracle              /**< pointer to store NLPIORACLE data structure */
   )
{
   SCIP_Bool nlpieval;

   assert(oracle != NULL);

   SCIPdebugMessage("%p oracle create\n", (void*)oracle);

   SCIP_CALL( SCIPallocMemory(scip, oracle) );
   BMSclearMemory(*oracle);

   SCIPdebugMessage("Oracle initializes expression interpreter %s\n", SCIPexprintGetName());
   SCIP_CALL( SCIPexprintCreate(scip, &(*oracle)->exprinterpreter) );

   SCIP_CALL( SCIPcreateClock(scip, &(*oracle)->evalclock) );

   SCIP_CALL( SCIPgetBoolParam(scip, "timing/nlpieval", &nlpieval) );
   if( !nlpieval )
      SCIPsetClockEnabled((*oracle)->evalclock, FALSE);

   /* create zero objective function */
   SCIP_CALL( createConstraint(scip, *oracle, &(*oracle)->objective, 0, NULL, NULL, NULL, 0.0, 0.0, NULL) );

   return SCIP_OKAY;
}

/** frees an NLPIORACLE data structure */
SCIP_RETCODE SCIPnlpiOracleFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE**     oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle  != NULL);
   assert(*oracle != NULL);

   SCIPdebugMessage("%p oracle free\n", (void*)oracle);

   invalidateJacobiSparsity(scip, *oracle);
   invalidateHessianLagSparsity(scip, *oracle);

   SCIP_CALL( freeConstraint(scip, *oracle, &(*oracle)->objective, FALSE) );
   SCIP_CALL( freeConstraints(scip, *oracle) );
   freeVariables(scip, *oracle);

   SCIP_CALL( SCIPfreeClock(scip, &(*oracle)->evalclock) );

   SCIP_CALL( SCIPexprintFree(scip, &(*oracle)->exprinterpreter) );

   if( (*oracle)->name != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleSetProblemName(scip, *oracle, NULL) );
   }

   BMSfreeMemory(oracle);

   return SCIP_OKAY;
}

/** sets the problem name (used for printing) */
SCIP_RETCODE SCIPnlpiOracleSetProblemName(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const char*           name                /**< name of problem */
   )
{
   assert(oracle != NULL);

   SCIPdebugMessage("%p set problem name\n", (void*)oracle);

   if( oracle->name != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &oracle->name, strlen(oracle->name)+1);
   }

   if( name != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &oracle->name, name, strlen(name)+1) );
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
   SCIP*                 scip,               /**< SCIP data structure */
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

   SCIP_CALL( ensureVarsSize(scip, oracle, oracle->nvars + nvars) );

   if( lbs != NULL )
   {
      BMScopyMemoryArray(&oracle->varlbs[oracle->nvars], lbs, nvars);
   }
   else
      for( i = 0; i < nvars; ++i )
         oracle->varlbs[oracle->nvars+i] = -SCIPinfinity(scip);

   if( ubs != NULL )
   {
      BMScopyMemoryArray(&oracle->varubs[oracle->nvars], ubs, nvars);

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
         oracle->varubs[oracle->nvars+i] =  SCIPinfinity(scip);

   if( varnames != NULL )
   {
      if( oracle->varnames == NULL )
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &oracle->varnames, oracle->varssize) );
      }

      for( i = 0; i < nvars; ++i )
      {
         if( varnames[i] != NULL )
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &oracle->varnames[oracle->nvars+i], varnames[i], strlen(varnames[i])+1) );
         }
         else
            oracle->varnames[oracle->nvars+i] = NULL;
      }
   }
   else if( oracle->varnames != NULL )
   {
      BMSclearMemoryArray(&oracle->varnames[oracle->nvars], nvars);
   }

   BMSclearMemoryArray(&oracle->varlincount[oracle->nvars], nvars);
   BMSclearMemoryArray(&oracle->varnlcount[oracle->nvars], nvars);

   /* @TODO update sparsity pattern by extending heslagoffsets */
   invalidateHessianLagSparsity(scip, oracle);

   oracle->nvars += nvars;

   return SCIP_OKAY;
}

/** adds constraints
 *
 *  linear coefficients: row(=constraint) oriented matrix;
 *  quadratic coefficients: row oriented matrix for each constraint
 */
SCIP_RETCODE SCIPnlpiOracleAddConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to add */
   const SCIP_Real*      lhss,               /**< array with left-hand sides of constraints, or NULL if all -infinity */
   const SCIP_Real*      rhss,               /**< array with right-hand sides of constraints, or NULL if all +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   SCIP_EXPR**           exprs,              /**< NULL if no nonlinear parts, otherwise exprs[.] gives nonlinear part,
                                              *   or NULL if no nonlinear part in this constraint */
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

   invalidateJacobiSparsity(scip, oracle); /* @TODO we could also update (extend) the sparsity pattern */

   SCIP_CALL( ensureConssSize(scip, oracle, oracle->nconss + nconss) );
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( createConstraint(scip, oracle, &cons,
            nlininds != NULL ? nlininds[c] : 0,
            lininds != NULL ? lininds[c] : NULL,
            linvals != NULL ? linvals[c] : NULL,
            exprs != NULL ? exprs[c] : NULL,
            lhss != NULL ? lhss[c] : -SCIPinfinity(scip),
            rhss != NULL ? rhss[c] :  SCIPinfinity(scip),
            consnames != NULL ? consnames[c] : NULL
            ) );

      if( cons->expr != NULL )
         addednlcon = TRUE;

      oracle->conss[oracle->nconss+c] = cons;
   }
   oracle->nconss += nconss;

   if( addednlcon == TRUE )
      invalidateHessianLagSparsity(scip, oracle);

   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected
 *
 *  May change sparsity pattern.
 */
SCIP_RETCODE SCIPnlpiOracleSetObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real       constant,           /**< constant part of objective */
   int                   nlin,               /**< number of linear variable coefficients */
   const int*            lininds,            /**< indices of linear variables, or NULL if no linear part */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables, or NULL if no linear part */
   SCIP_EXPR*            expr                /**< expression of nonlinear part, or NULL if no nonlinear part */
   )
{  /*lint --e{715}*/
   assert(oracle != NULL);
   assert(!SCIPisInfinity(scip, REALABS(constant)));

   SCIPdebugMessage("%p set objective\n", (void*)oracle);

   if( expr != NULL || oracle->objective->expr != NULL )
      invalidateHessianLagSparsity(scip, oracle);

   /* clear previous objective */
   SCIP_CALL( freeConstraint(scip, oracle, &oracle->objective, TRUE) );

   /* create new objective */
   SCIP_CALL( createConstraint(scip, oracle, &oracle->objective,
         nlin, lininds, linvals, expr, constant, constant, NULL) );

   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiOracleChgVarBounds(
   SCIP*                 scip,               /**< SCIP data structure */
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

      oracle->varlbs[indices[i]] = (lbs != NULL ? lbs[i] : -SCIPinfinity(scip));
      oracle->varubs[indices[i]] = (ubs != NULL ? ubs[i] :  SCIPinfinity(scip));

      if( oracle->varlbs[indices[i]] > oracle->varubs[indices[i]] )
      {
         /* inconsistent bounds; let's assume it's due to rounding and make them equal */
         assert(EPSEQ(oracle->varlbs[indices[i]], oracle->varubs[indices[i]], SCIP_DEFAULT_EPSILON));
         oracle->varlbs[indices[i]] = oracle->varubs[indices[i]];
      }
   }

   return SCIP_OKAY;
}

/** change constraint sides */
SCIP_RETCODE SCIPnlpiOracleChgConsSides(
   SCIP*                 scip,               /**< SCIP data structure */
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

      oracle->conss[indices[i]]->lhs = (lhss != NULL ? lhss[i] : -SCIPinfinity(scip));
      oracle->conss[indices[i]]->rhs = (rhss != NULL ? rhss[i] :  SCIPinfinity(scip));
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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< deletion status of vars in input (1 if var should be deleted, 0 if not);
                                              *   new position of var in output (-1 if var was deleted) */
   )
{  /*lint --e{715}*/
   int c;
   int lastgood; /* index of the last variable that should be kept */
   SCIP_NLPIORACLECONS* cons;
   SCIP_EXPRITER* it;

   assert(oracle != NULL);

   SCIPdebugMessage("%p del var set\n", (void*)oracle);

   invalidateJacobiSparsity(scip, oracle);
   invalidateHessianLagSparsity(scip, oracle);

   lastgood = oracle->nvars - 1;
   while( lastgood >= 0 && delstats[lastgood] == 1 )
      --lastgood;
   if( lastgood < 0 )
   {
      /* all variables should be deleted */
      assert(oracle->nconss == 0); /* we could relax this by checking that all constraints are constant */
      oracle->objective->nlinidxs = 0;
      for( c = 0; c < oracle->nvars; ++c )
         delstats[c] = -1;
      freeVariables(scip, oracle);
      return SCIP_OKAY;
   }

   /* delete variables at the end */
   for( c = oracle->nvars - 1; c > lastgood; --c )
   {
      if( oracle->varnames && oracle->varnames[c] != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &oracle->varnames[c], strlen(oracle->varnames[c])+1);
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

      if( oracle->varnames != NULL && oracle->varnames[c] != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &oracle->varnames[c], strlen(oracle->varnames[c])+1);
      }
      delstats[c] = -1;

      /* move variable at position lastgood to position c */
      SCIP_CALL( moveVariable(scip, oracle, lastgood, c) );
      delstats[lastgood] = c; /* mark that lastgood variable is now at position c */

      /* move lastgood forward, delete variables on the way */
      --lastgood;
      while( lastgood > c && delstats[lastgood] == 1)
      {
         if( oracle->varnames && oracle->varnames[lastgood] != NULL )
         {
            SCIPfreeBlockMemoryArray(scip, &oracle->varnames[lastgood], strlen(oracle->varnames[lastgood])+1);
         }
         delstats[lastgood] = -1;
         --lastgood;
      }
   }
   assert(c == lastgood);

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );

   for( c = -1; c < oracle->nconss; ++c )
   {
      cons = c < 0 ? oracle->objective : oracle->conss[c];
      assert(cons != NULL);

      /* update indices in linear part, sort indices, and then clear elements that are marked as deleted */
      mapIndices(delstats, cons->nlinidxs, cons->linidxs);
      SCIPsortIntReal(cons->linidxs, cons->lincoefs, cons->nlinidxs);
      clearDeletedLinearElements(&cons->linidxs, &cons->lincoefs, &cons->nlinidxs);

      if( cons->expr != NULL )
      {
         /* update variable indices in varidx expressions */
         SCIP_EXPR* expr;
         SCIP_Bool keptvar = FALSE;  /* whether any of the variables in expr was not deleted */
#ifndef NDEBUG
         SCIP_Bool delvar = FALSE; /* whether any of the variables in expr was deleted */
#endif

         SCIP_CALL( SCIPexpriterInit(it, cons->expr, SCIP_EXPRITER_DFS, FALSE) );
         for( expr = cons->expr; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
         {
            if( !SCIPisExprVaridx(scip, expr) )
               continue;

            if( delstats[SCIPgetIndexExprVaridx(expr)] >= 0 )
            {
               /* if variable is not deleted, then set its new index */
               keptvar = TRUE;
               SCIPsetIndexExprVaridx(expr, delstats[SCIPgetIndexExprVaridx(expr)]);

               /* if variable is kept, then there must not have been any variable that was deleted */
               assert(!delvar);
            }
            else
            {
#ifndef NDEBUG
               delvar = TRUE;
#endif
               /* if variable is deleted, then there must not have been any variable that was kept
                * (either all variables are deleted, which removes the expr, or none)
                */
               assert(!keptvar);
            }
         }
         if( !keptvar )
         {
            SCIP_CALL( SCIPexprintFreeData(scip, oracle->exprinterpreter, cons->expr, &cons->exprintdata) );
            SCIP_CALL( SCIPreleaseExpr(scip, &cons->expr) );
         }
      }
   }

   SCIPfreeExpriter(&it);

   oracle->nvars = lastgood+1;

   return SCIP_OKAY;
}

/** deletes a set of constraints */
SCIP_RETCODE SCIPnlpiOracleDelConsSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of rows in input (1 if row should be deleted, 0 if not);
                                              *   new position of row in output (-1 if row was deleted) */
   )
{  /*lint --e{715}*/
   int c;
   int lastgood; /* index of the last constraint that should be kept */

   assert(oracle != NULL);

   SCIPdebugMessage("%p del cons set\n", (void*)oracle);

   invalidateJacobiSparsity(scip, oracle);
   invalidateHessianLagSparsity(scip, oracle);

   lastgood = oracle->nconss - 1;
   while( lastgood >= 0 && delstats[lastgood] == 1)
      --lastgood;
   if( lastgood < 0 )
   {
      /* all constraints should be deleted */
      for( c = 0; c < oracle->nconss; ++c )
         delstats[c] = -1;
      SCIP_CALL( freeConstraints(scip, oracle) );

      /* the previous call did not keep variable counts uptodate
       * since we only have an objective function left, we reset the counts to the ones of the objective
       */
      BMSclearMemoryArray(oracle->varlincount, oracle->nvars);
      BMSclearMemoryArray(oracle->varnlcount, oracle->nvars);
      SCIP_CALL( updateVariableCounts(scip, oracle, 1, oracle->objective->nlinidxs, oracle->objective->linidxs, oracle->objective->expr) );

      return SCIP_OKAY;
   }

   /* delete constraints at the end */
   for( c = oracle->nconss - 1; c > lastgood; --c )
   {
      SCIP_CALL( freeConstraint(scip, oracle, &oracle->conss[c], TRUE) );
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

      SCIP_CALL( freeConstraint(scip, oracle, &oracle->conss[c], TRUE) );
      assert(oracle->conss[c] == NULL);
      delstats[c] = -1;

      /* move constraint at position lastgood to position c */
      oracle->conss[c] = oracle->conss[lastgood];
      assert(oracle->conss[c] != NULL);
      delstats[lastgood] = c; /* mark that lastgood constraint is now at position c */
      oracle->conss[lastgood] = NULL;
      --lastgood;

      /* move lastgood forward, delete constraints on the way */
      while( lastgood > c && delstats[lastgood] == 1)
      {
         SCIP_CALL( freeConstraint(scip, oracle, &oracle->conss[lastgood], TRUE) );
         assert(oracle->conss[lastgood] == NULL);
         delstats[lastgood] = -1;
         --lastgood;
      }
   }
   assert(c == lastgood+1);

   oracle->nconss = lastgood+1;

   return SCIP_OKAY;
}

/** changes (or adds) linear coefficients in one constraint or objective */
SCIP_RETCODE SCIPnlpiOracleChgLinearCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
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

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &cons->linidxs,  varidxs,  nentries) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &cons->lincoefs, newcoefs, nentries) );
      cons->linsize  = nentries;
      cons->nlinidxs = nentries;

      SCIP_CALL( updateVariableCounts(scip, oracle, 1, nentries, varidxs, NULL) );

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

            if( newcoefs[i] == 0.0 )
            {
               --oracle->varlincount[varidxs[i]];
               assert(oracle->varlincount[varidxs[i]] >= 0);
            }
         }
         else if( newcoefs[i] != 0.0 )  /*lint !e613*/
         {
            /* append new entry */
            SCIPdebugMessage("add coefficient of var %d at pos %d, value %g\n", varidxs[i], cons->nlinidxs, newcoefs[i]);  /*lint !e613*/

            SCIP_CALL( ensureConsLinSize(scip, cons, cons->nlinidxs + (nentries-i)) );
            cons->linidxs[cons->nlinidxs]  = varidxs[i];   /*lint !e613*/
            cons->lincoefs[cons->nlinidxs] = newcoefs[i];  /*lint !e613*/
            ++cons->nlinidxs;

            ++oracle->varlincount[varidxs[i]];

            needsort = TRUE;
         }
      }
   }

   if( needsort )
   {
      invalidateJacobiSparsity(scip, oracle);
      sortLinearCoefficients(&cons->nlinidxs, cons->linidxs, cons->lincoefs);
   }

   return SCIP_OKAY;
}

/** replaces expression of one constraint or objective */
SCIP_RETCODE SCIPnlpiOracleChgExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where expression should be changed, or -1 for objective */
   SCIP_EXPR*            expr                /**< new expression, or NULL */
   )
{
   SCIP_NLPIORACLECONS* cons;

   SCIPdebugMessage("%p chg expr\n", (void*)oracle);

   assert(oracle != NULL);
   assert(considx >= -1);
   assert(considx < oracle->nconss);

   invalidateHessianLagSparsity(scip, oracle);
   invalidateJacobiSparsity(scip, oracle);

   cons = considx < 0 ? oracle->objective : oracle->conss[considx];

   /* free previous expression */
   if( cons->expr != NULL )
   {
      SCIP_CALL( updateVariableCounts(scip, oracle, -1, 0, NULL, cons->expr) );
      SCIP_CALL( SCIPexprintFreeData(scip, oracle->exprinterpreter, cons->expr, &cons->exprintdata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &cons->expr) );
   }

   /* if user did not want to set new expr, then we are done */
   if( expr == NULL )
      return SCIP_OKAY;

   assert(oracle->exprinterpreter != NULL);

   /* install new expression */
   cons->expr = expr;
   SCIPcaptureExpr(cons->expr);
   SCIP_CALL( SCIPexprintCompile(scip, oracle->exprinterpreter, cons->expr, &cons->exprintdata) );

   /* keep variable counts up to date */
   SCIP_CALL( updateVariableCounts(scip, oracle, 1, 0, NULL, cons->expr) );

   return SCIP_OKAY;
}

/** changes the constant value in the objective function */ /*lint -e{715}*/
SCIP_RETCODE SCIPnlpiOracleChgObjConstant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real             objconstant         /**< new value for objective constant */
   )
{ /*lint --e{715}*/
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

/** indicates whether variable appears nonlinear in any objective or constraint */  /*lint --e{715}*/
SCIP_Bool SCIPnlpiOracleIsVarNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   varidx              /**< the variable to check */
   )
{
   assert(oracle != NULL);
   assert(varidx >= 0);
   assert(varidx < oracle->nvars);
   assert(oracle->varnlcount != NULL);

   return oracle->varnlcount[varidx] > 0;
}

/** returns number of linear and nonlinear appearances of variables in objective and constraints */  /*lint --e{715}*/
void SCIPnlpiOracleGetVarCounts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           lincounts,          /**< buffer to return pointer to array of counts of linear appearances */
   const int**           nlcounts            /**< buffer to return pointer to array of counts of nonlinear appearances */
   )
{
   assert(oracle != NULL);
   assert(lincounts != NULL);
   assert(nlcounts != NULL);

   *lincounts = oracle->varlincount;
   *nlcounts = oracle->varnlcount;
}

/** gives constant term of objective */
SCIP_Real SCIPnlpiOracleGetObjectiveConstant(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);
   assert(oracle->objective->lhs == oracle->objective->rhs);  /*lint !e777*/

   return oracle->objective->lhs;
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

/** indicates whether constraint is nonlinear */
SCIP_Bool SCIPnlpiOracleIsConstraintNonlinear(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< index of constraint for which nonlinearity status is returned, or -1 for objective */
   )
{
   SCIP_NLPIORACLECONS* cons;

   assert(oracle != NULL);
   assert(considx >= -1);
   assert(considx < oracle->nconss);

   cons = considx < 0 ? oracle->objective : oracle->conss[considx];

   return cons->expr != NULL;
}

/** gives the evaluation capabilities that are shared among all expressions in the problem */
SCIP_EXPRINTCAPABILITY SCIPnlpiOracleGetEvalCapability(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   int c;
   SCIP_EXPRINTCAPABILITY evalcapability;

   assert(oracle != NULL);

   if( oracle->objective->expr != NULL )
      evalcapability = SCIPexprintGetExprCapability(scip, oracle->exprinterpreter, oracle->objective->expr, oracle->objective->exprintdata);
   else
      evalcapability = SCIP_EXPRINTCAPABILITY_ALL;

   for( c = 0; c < oracle->nconss; ++c )
      if( oracle->conss[c]->expr != NULL )
         evalcapability &= SCIPexprintGetExprCapability(scip, oracle->exprinterpreter, oracle->conss[c]->expr, oracle->conss[c]->exprintdata);

   return evalcapability;
}

/** evaluates the objective function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            objval              /**< pointer to store objective value */
   )
{
   SCIP_RETCODE retcode;

   assert(oracle != NULL);

   SCIPdebugMessage("%p eval obj value\n", (void*)oracle);

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );
   retcode = evalFunctionValue(scip, oracle, oracle->objective, x, objval);
   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   assert(oracle->objective->lhs == oracle->objective->rhs);  /*lint !e777*/
   if( retcode == SCIP_OKAY )
      *objval += oracle->objective->lhs;

   return retcode;
}

/** evaluates one constraint function in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint to evaluate */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            conval              /**< pointer to store constraint value */
   )
{
   SCIP_RETCODE retcode;

   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);

   SCIPdebugMessage("%p eval cons value\n", (void*)oracle);

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );
   retcode = evalFunctionValue(scip, oracle, oracle->conss[considx], x, conval);
   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   return retcode;
}

/** evaluates all constraint functions in a given point */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            convals             /**< buffer to store constraint values */
   )
{
   SCIP_RETCODE retcode = SCIP_OKAY;
   int i;

   SCIPdebugMessage("%p eval cons values\n", (void*)oracle);

   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(convals != NULL);

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );
   for( i = 0; i < oracle->nconss; ++i )
   {
      retcode = evalFunctionValue(scip, oracle, oracle->conss[i], x, &convals[i]);
      if( retcode != SCIP_OKAY )
         break;
   }
   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   return retcode;
}

/** computes the objective gradient in a given point
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            objval,             /**< pointer to store objective value */
   SCIP_Real*            objgrad             /**< pointer to store (dense) objective gradient */
   )
{
   SCIP_RETCODE retcode;
   assert(oracle != NULL);

   SCIPdebugMessage("%p eval obj grad\n", (void*)oracle);

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );
   retcode = evalFunctionGradient(scip, oracle, oracle->objective, x, isnewx, objval, objgrad);
   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   assert(oracle->objective->lhs == oracle->objective->rhs);  /*lint !e777*/
   if( retcode == SCIP_OKAY )
      *objval += oracle->objective->lhs;

   return retcode;
}

/** computes a constraints gradient in a given point
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to compute gradient for */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            conval,             /**< pointer to store constraint value */
   SCIP_Real*            congrad             /**< pointer to store (dense) constraint gradient */
   )
{
   SCIP_RETCODE retcode;

   assert(oracle != NULL);
   assert(x != NULL || oracle->nvars == 0);
   assert(conval != NULL);

   SCIPdebugMessage("%p eval cons grad\n", (void*)oracle);

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );
   retcode = evalFunctionGradient(scip, oracle, oracle->conss[considx], x, isnewx, conval, congrad);
   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   return retcode;
}

/** gets sparsity pattern (rowwise) of Jacobian matrix
 *
 *  Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 *  Adding or deleting constraints destroys the sparsity structure and make another call to this function necessary.
 */
SCIP_RETCODE SCIPnlpiOracleGetJacobianSparsity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, offset[nconss] gives length of col, can be NULL */
   )
{
   SCIP_NLPIORACLECONS* cons;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
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

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &oracle->jacoffsets, oracle->nconss + 1) );

   maxnnz = MIN(oracle->nvars, 10) * oracle->nconss;  /* initial guess */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &oracle->jaccols, maxnnz) );

   if( maxnnz == 0 )
   {
      /* no variables */
      BMSclearMemoryArray(oracle->jacoffsets, oracle->nconss + 1);
      if( offset != NULL )
         *offset = oracle->jacoffsets;
      if( col != NULL )
         *col = oracle->jaccols;

      SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

      return SCIP_OKAY;
   }
   nnz = 0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nzflag, oracle->nvars) );

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   for( i = 0; i < oracle->nconss; ++i )
   {
      oracle->jacoffsets[i] = nnz;

      cons = oracle->conss[i];
      assert(cons != NULL);

      if( cons->expr == NULL )
      {
         /* for a linear constraint, we can just copy the linear indices from the constraint into the sparsity pattern */
         if( cons->nlinidxs > 0 )
         {
            SCIP_CALL( ensureIntArraySize(scip, &oracle->jaccols, &maxnnz, nnz + cons->nlinidxs) );
            BMScopyMemoryArray(&oracle->jaccols[nnz], cons->linidxs, cons->nlinidxs);
            nnz += cons->nlinidxs;
         }
         continue;
      }

      /* check which variables appear in constraint i
       * @todo this could be done faster for very sparse constraint by assembling all appearing variables, sorting, and removing duplicates
       */
      BMSclearMemoryArray(nzflag, oracle->nvars);

      for( j = 0; j < cons->nlinidxs; ++j )
         nzflag[cons->linidxs[j]] = TRUE;

      for( expr = SCIPexpriterRestartDFS(it, cons->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
         if( SCIPisExprVaridx(scip, expr) )
         {
            assert(SCIPgetIndexExprVaridx(expr) < oracle->nvars);
            nzflag[SCIPgetIndexExprVaridx(expr)] = TRUE;
         }

      /* store variables indices in jaccols */
      for( j = 0; j < oracle->nvars; ++j )
      {
         if( nzflag[j] == FALSE )
            continue;

         SCIP_CALL( ensureIntArraySize(scip, &oracle->jaccols, &maxnnz, nnz + 1) );
         oracle->jaccols[nnz] = j;
         ++nnz;
      }
   }

   SCIPfreeExpriter(&it);

   oracle->jacoffsets[oracle->nconss] = nnz;

   /* shrink jaccols array to nnz */
   if( nnz < maxnnz )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &oracle->jaccols, maxnnz, nnz) );
   }

   SCIPfreeBlockMemoryArray(scip, &nzflag, oracle->nvars);

   if( offset != NULL )
      *offset = oracle->jacoffsets;
   if( col != NULL )
      *col = oracle->jaccols;

   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   return SCIP_OKAY;
}

/** evaluates the Jacobian matrix in a given point
 *
 *  The values in the Jacobian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetJacobianSparsity().
 *  The user need to call SCIPnlpiOracleGetJacobianSparsity() at least ones before using this function.
 *
 * @return SCIP_INVALIDDATA, if the Jacobian could not be evaluated (domain error, etc.)
 */
SCIP_RETCODE SCIPnlpiOracleEvalJacobian(
   SCIP*                 scip,               /**< SCIP data structure */
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

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &grad, oracle->nvars) );

   retcode = SCIP_OKAY;

   j = oracle->jacoffsets[0];  /* TODO isn't oracle->jacoffsets[0] == 0 and thus always j == k ? */
   k = 0;
   for( i = 0; i < oracle->nconss; ++i )
   {
      cons = oracle->conss[i];
      assert(cons != NULL);

      if( cons->expr == NULL )
      {
         if( convals != NULL )
            convals[i] = 0.0;

         /* for a linear constraint, we can just copy the linear coefs from the constraint into the jacobian */
         if( cons->nlinidxs > 0 )
         {
            BMScopyMemoryArray(&jacobi[k], cons->lincoefs, cons->nlinidxs);
            j += cons->nlinidxs;
            k += cons->nlinidxs;
            if( convals != NULL )
               for( l = 0; l < cons->nlinidxs; ++l )
                  convals[i] += cons->lincoefs[l] * x[cons->linidxs[l]];
         }
         assert(j == oracle->jacoffsets[i+1]);
         continue;
      }

      /* eval grad for nonlinear and add to jacobi */
      SCIPdebugMsg(scip, "eval gradient of ");
      SCIPdebug( if( isnewx ) {printf("\nx ="); for( l = 0; l < oracle->nvars; ++l) printf(" %g", x[l]); printf("\n");} )

      SCIP_CALL( SCIPexprintGrad(scip, oracle->exprinterpreter, cons->expr, cons->exprintdata, (SCIP_Real*)x, isnewx, &nlval, grad) );

      SCIPdebug( printf("g ="); for( l = oracle->jacoffsets[i]; l < oracle->jacoffsets[i+1]; ++l) printf(" %g", grad[oracle->jaccols[l]]); printf("\n"); )

      if( !SCIPisFinite(nlval) || SCIPisInfinity(scip, ABS(nlval)) )
      {
         SCIPdebugMessage("gradient evaluation yield invalid function value %g\n", nlval);
         retcode = SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
         goto TERMINATE;
      }
      if( convals != NULL )
         convals[i] = nlval;

      /* add linear part to grad */
      for( l = 0; l < cons->nlinidxs; ++l )
      {
         if( convals != NULL )
            convals[i] += cons->lincoefs[l] * x[cons->linidxs[l]];
         /* if grad[cons->linidxs[l]] is not finite, then adding a finite value doesn't change that, so don't check that here */
         grad[cons->linidxs[l]] += cons->lincoefs[l];
      }

      /* store complete gradient (linear + nonlinear) in jacobi
       * use the already evaluated sparsity pattern to pick only elements from grad that could have been set
       */
      assert(j == oracle->jacoffsets[i]);
      for( ; j < oracle->jacoffsets[i+1]; ++j )
      {
         if( !SCIPisFinite(grad[oracle->jaccols[j]]) )
         {
            SCIPdebugMessage("gradient evaluation yield invalid gradient value %g\n", grad[l]);
            retcode = SCIP_INVALIDDATA; /* indicate that the function could not be evaluated at given point */
            goto TERMINATE;
         }
         jacobi[k++] = grad[oracle->jaccols[j]];
         /* reset to 0 for next constraint */
         grad[oracle->jaccols[j]] = 0.0;
      }

#ifndef NDEBUG
      /* check that exprint really wrote only into expected elements of grad
       * TODO remove after some testing for better performance of debug runs */
      for( l = 0; l < oracle->nvars; ++l )
         assert(grad[l] == 0.0);
#endif
   }

TERMINATE:
   /* if there was an eval error, then we may have interrupted before cleaning up the grad buffer */
   if( retcode == SCIP_INVALIDDATA )
      BMSclearMemoryArray(grad, oracle->nvars);

   SCIPfreeCleanBufferArray(scip, &grad);

   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   return retcode;
}

/** gets sparsity pattern of the Hessian matrix of the Lagrangian
 *
 *  Note that internal data is returned in *offset and *col, thus the user must not to allocate memory there.
 *  Adding or deleting variables, objective, or constraints may destroy the sparsity structure and make another call to this function necessary.
 *  Only elements of the lower left triangle and the diagonal are counted.
 */
SCIP_RETCODE SCIPnlpiOracleGetHessianLagSparsity(
   SCIP*                 scip,               /**< SCIP data structure */
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

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &oracle->heslagoffsets, oracle->nvars + 1) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &colnz,  oracle->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &collen, oracle->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &colnnz, oracle->nvars) );
   BMSclearMemoryArray(colnz,  oracle->nvars);
   BMSclearMemoryArray(collen, oracle->nvars);
   BMSclearMemoryArray(colnnz, oracle->nvars);
   nnz = 0;

   if( oracle->objective->expr != NULL )
   {
      SCIP_CALL( hessLagSparsitySetNzFlagForExpr(scip, oracle, colnz, collen, colnnz, &nnz, oracle->objective->expr, oracle->objective->exprintdata, oracle->nvars) );
   }

   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->conss[i]->expr != NULL )
      {
         SCIP_CALL( hessLagSparsitySetNzFlagForExpr(scip, oracle, colnz, collen, colnnz, &nnz, oracle->conss[i]->expr, oracle->conss[i]->exprintdata, oracle->nvars) );
      }
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &oracle->heslagcols, nnz) );

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
      SCIPfreeBlockMemoryArrayNull(scip, &colnz[i], collen[i]);
      collen[i] = 0;
   }
   oracle->heslagoffsets[oracle->nvars] = cnt;
   assert(cnt == nnz);

   SCIPfreeBlockMemoryArray(scip, &colnz,  oracle->nvars);
   SCIPfreeBlockMemoryArray(scip, &colnnz, oracle->nvars);
   SCIPfreeBlockMemoryArray(scip, &collen, oracle->nvars);

   if( offset != NULL )
      *offset = oracle->heslagoffsets;
   if( col != NULL )
      *col = oracle->heslagcols;

   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   return SCIP_OKAY;
}

/** evaluates the Hessian matrix of the Lagrangian in a given point
 *
 *  The values in the Hessian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetHessianLagSparsity().
 *  The user must call SCIPnlpiOracleGetHessianLagSparsity() at least ones before using this function.
 *  Only elements of the lower left triangle and the diagonal are computed.
 *
 * @return SCIP_INVALIDDATA, if the Hessian could not be evaluated (domain error, etc.)
 */
SCIP_RETCODE SCIPnlpiOracleEvalHessianLag(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx_obj,         /**< has the point x changed since the last call to an objective evaluation function? */
   SCIP_Bool             isnewx_cons,        /**< has the point x changed since the last call to the constraint evaluation function? */
   SCIP_Real             objfactor,          /**< weight for objective function */
   const SCIP_Real*      lambda,             /**< weights (Lagrangian multipliers) for the constraints */
   SCIP_Real*            hessian             /**< pointer to store sparse hessian values */
   )
{  /*lint --e{715}*/
   SCIP_RETCODE retcode = SCIP_OKAY;
   int i;

   assert(oracle != NULL);
   assert(x != NULL);
   assert(lambda != NULL || oracle->nconss == 0);
   assert(hessian != NULL);

   assert(oracle->heslagoffsets != NULL);
   assert(oracle->heslagcols != NULL);

   SCIPdebugMessage("%p eval hessian lag\n", (void*)oracle);

   SCIP_CALL( SCIPstartClock(scip, oracle->evalclock) );

   BMSclearMemoryArray(hessian, oracle->heslagoffsets[oracle->nvars]);

   if( objfactor != 0.0 && oracle->objective->expr != NULL )
   {
      retcode = hessLagAddExpr(scip, oracle, objfactor, x, isnewx_obj, oracle->objective->expr, oracle->objective->exprintdata, oracle->heslagoffsets, oracle->heslagcols, hessian);
   }

   for( i = 0; i < oracle->nconss && retcode == SCIP_OKAY; ++i )
   {
      assert( lambda != NULL ); /* for lint */
      if( lambda[i] == 0.0 || oracle->conss[i]->expr == NULL )
         continue;
      retcode = hessLagAddExpr(scip, oracle, lambda[i], x, isnewx_cons, oracle->conss[i]->expr, oracle->conss[i]->exprintdata, oracle->heslagoffsets, oracle->heslagcols, hessian);
   }

   SCIP_CALL( SCIPstopClock(scip, oracle->evalclock) );

   return retcode;
}

/** resets clock that measures evaluation time */
SCIP_RETCODE SCIPnlpiOracleResetEvalTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   SCIP_CALL( SCIPresetClock(scip, oracle->evalclock) );

   return SCIP_OKAY;
}

/** gives time spend in evaluation since last reset of clock
 *
 * Gives 0 if the eval clock is disabled.
 */
SCIP_Real SCIPnlpiOracleGetEvalTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   )
{
   assert(oracle != NULL);

   return SCIPgetClockTime(scip, oracle->evalclock);
}

/** prints the problem to a file. */
SCIP_RETCODE SCIPnlpiOraclePrintProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
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

   SCIPinfoMessage(scip, file, "NLPI Oracle %s: %d variables and %d constraints\n", oracle->name ? oracle->name : "", oracle->nvars, oracle->nconss);
   for( i = 0; i < oracle->nvars; ++i )
   {
      if( oracle->varnames != NULL && oracle->varnames[i] != NULL )
         SCIPinfoMessage(scip, file, "%10s (x%d)", oracle->varnames[i], i);  /* give also name x%d as it will be by expression-print (printFunction) */
      else
         SCIPinfoMessage(scip, file, "x%09d", i);
      SCIPinfoMessage(scip, file, ": [%8g, %8g]", oracle->varlbs[i], oracle->varubs[i]);
      SCIPinfoMessage(scip, file, "\t #linear: %d #nonlinear: %d\n", oracle->varlincount[i], oracle->varnlcount[i]);
   }

   SCIPinfoMessage(scip, file, "objective: ");
   SCIP_CALL( printFunction(scip, oracle, file, oracle->objective, FALSE) );
   if( oracle->objective->lhs != 0.0 )
      SCIPinfoMessage(scip, file, "%+.15g", oracle->objective->lhs);
   SCIPinfoMessage(scip, file, "\n");

   for( i = 0; i < oracle->nconss; ++i )
   {
      if( oracle->conss[i]->name != NULL )
         SCIPinfoMessage(scip, file, "%10s", oracle->conss[i]->name);
      else
         SCIPinfoMessage(scip, file, "con%07d", i);

      lhs = oracle->conss[i]->lhs;
      rhs = oracle->conss[i]->rhs;
      SCIPinfoMessage(scip, file, ": ");
      if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && lhs != rhs )
         SCIPinfoMessage(scip, file, "%.15g <= ", lhs);

      SCIP_CALL( printFunction(scip, oracle, file, oracle->conss[i], FALSE) );

      if( lhs == rhs )
         SCIPinfoMessage(scip, file, " = %.15g", rhs);
      else if( !SCIPisInfinity(scip, rhs) )
         SCIPinfoMessage(scip, file, " <= %.15g", rhs);
      else if( !SCIPisInfinity(scip, -lhs) )
         SCIPinfoMessage(scip, file, " >= %.15g", lhs);

      SCIPinfoMessage(scip, file, "\n");
   }

   return SCIP_OKAY;
}

/** prints the problem to a file in GAMS format
 *
 * If there are variable (equation, resp.) names with more than 9 characters, then variable (equation, resp.) names are prefixed with an unique identifier.
 * This is to make it easier to identify variables solution output in the listing file.
 * Names with more than 64 characters are shorten to 64 letters due to GAMS limits.
 */
SCIP_RETCODE SCIPnlpiOraclePrintProblemGams(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real*            initval,            /**< starting point values for variables or NULL */
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

   SCIPinfoMessage(scip, file, "$offlisting\n");
   SCIPinfoMessage(scip, file, "$offdigit\n");
   SCIPinfoMessage(scip, file, "* NLPI Oracle Problem %s\n", oracle->name ? oracle->name : "");
   SCIPinfoMessage(scip, file, "Variables ");
   for( i = 0; i < oracle->nvars; ++i )
   {
      printName(namebuf, oracle->varnames != NULL ? oracle->varnames[i] : NULL, i, 'x', NULL, havelongvarnames);
      SCIPinfoMessage(scip, file, "%s, ", namebuf);
      if( i % 10 == 9 )
         SCIPinfoMessage(scip, file, "\n");
   }
   SCIPinfoMessage(scip, file, "NLPIORACLEOBJVAR;\n\n");
   for( i = 0; i < oracle->nvars; ++i )
   {
      char* name;
      name = oracle->varnames != NULL ? oracle->varnames[i] : NULL;
      if( oracle->varlbs[i] == oracle->varubs[i] )
      {
         printName(namebuf, name, i, 'x', NULL, havelongvarnames);
         SCIPinfoMessage(scip, file, "%s.fx = %.15g;\t", namebuf, oracle->varlbs[i]);
      }
      else
      {
         if( !SCIPisInfinity(scip, -oracle->varlbs[i]) )
         {
            printName(namebuf, name, i, 'x', NULL, havelongvarnames);
            SCIPinfoMessage(scip, file, "%s.lo = %.15g;\t", namebuf, oracle->varlbs[i]);
         }
         if( !SCIPisInfinity(scip, oracle->varubs[i]) )
         {
            printName(namebuf, name, i, 'x', NULL, havelongvarnames);
            SCIPinfoMessage(scip, file, "%s.up = %.15g;\t", namebuf, oracle->varubs[i]);
         }
      }
      if( initval != NULL )
      {
         printName(namebuf, name, i, 'x', NULL, havelongvarnames);
         SCIPinfoMessage(scip, file, "%s.l = %.15g;\t", namebuf, initval[i]);
      }
      SCIPinfoMessage(scip, file, "\n");
   }
   SCIPinfoMessage(scip, file, "\n");

   SCIPinfoMessage(scip, file, "Equations ");
   for( i = 0; i < oracle->nconss; ++i )
   {
      printName(namebuf, oracle->conss[i]->name, i, 'e', NULL, havelongequnames);
      SCIPinfoMessage(scip, file, "%s, ", namebuf);

      if( !SCIPisInfinity(scip, -oracle->conss[i]->lhs) && !SCIPisInfinity(scip, oracle->conss[i]->rhs) && oracle->conss[i]->lhs != oracle->conss[i]->rhs )
      {
         /* ranged row: add second constraint */
         printName(namebuf, oracle->conss[i]->name, i, 'e', "_RNG", havelongequnames);
         SCIPinfoMessage(scip, file, "%s, ", namebuf);
      }
      if( i % 10 == 9 )
         SCIPinfoMessage(scip, file, "\n");
   }
   SCIPinfoMessage(scip, file, "NLPIORACLEOBJ;\n\n");

   SCIPinfoMessage(scip, file, "NLPIORACLEOBJ.. NLPIORACLEOBJVAR =E= ");
   SCIP_CALL( printFunction(scip, oracle, file, oracle->objective, havelongvarnames) );
   if( oracle->objective->lhs != 0.0 )
      SCIPinfoMessage(scip, file, "%+.15g", oracle->objective->lhs);
   SCIPinfoMessage(scip, file, ";\n");

   for( i = 0; i < oracle->nconss; ++i )
   {
      SCIP_Real lhs;
      SCIP_Real rhs;

      printName(namebuf, oracle->conss[i]->name, i, 'e', NULL, havelongequnames);
      SCIPinfoMessage(scip, file, "%s.. ", namebuf);

      SCIP_CALL( printFunction(scip, oracle, file, oracle->conss[i], havelongvarnames) );

      lhs = oracle->conss[i]->lhs;
      rhs = oracle->conss[i]->rhs;

      if( lhs == rhs )
         SCIPinfoMessage(scip, file, " =E= %.15g", rhs);
      else if( !SCIPisInfinity(scip, rhs) )
         SCIPinfoMessage(scip, file, " =L= %.15g", rhs);
      else if( !SCIPisInfinity(scip, -lhs) )
         SCIPinfoMessage(scip, file, " =G= %.15g", lhs);
      else
         SCIPinfoMessage(scip, file, " =N= 0");
      SCIPinfoMessage(scip, file, ";\n");

      if( !SCIPisInfinity(scip, lhs) && !SCIPisInfinity(scip, rhs) && lhs != rhs )
      {
         printName(namebuf, oracle->conss[i]->name, i, 'e', "_RNG", havelongequnames);
         SCIPinfoMessage(scip, file, "%s.. ", namebuf);

         SCIP_CALL( printFunction(scip, oracle, file, oracle->conss[i], havelongvarnames) );

         SCIPinfoMessage(scip, file, " =G= %.15g;\n", lhs);
      }

      if( nllevel <= 1 && oracle->conss[i]->expr != NULL )
         nllevel = 2;
      if( nllevel <= 2 && oracle->conss[i]->expr != NULL )
      {
         SCIP_Bool nonsmooth;
         SCIP_CALL( exprIsNonSmooth(scip, oracle->conss[i]->expr, &nonsmooth) );
         if( nonsmooth )
            nllevel = 3;
      }
   }

   (void) SCIPsnprintf(problemname, SCIP_MAXSTRLEN, "%s", oracle->name ? oracle->name : "m");

   SCIPinfoMessage(scip, file, "Model %s / all /;\n", problemname);
   SCIPinfoMessage(scip, file, "option limrow = 0;\n");
   SCIPinfoMessage(scip, file, "option limcol = 0;\n");
   SCIPinfoMessage(scip, file, "Solve %s minimizing NLPIORACLEOBJVAR using %s;\n", problemname, nllevelname[nllevel]);

   return SCIP_OKAY;
}

/**@} */
