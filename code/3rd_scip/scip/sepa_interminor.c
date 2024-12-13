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

/**@file   sepa_interminor.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  minor separator with intersection cuts
 * @author Felipe Serrano
 * @author Antonia Chmiela
 *
 * Let X be the matrix of auxiliary variables added for bilinear terms, X_{ij} = x_ix_j.
 * The separator enforces quadratic constraints det(2x2 minor of X) = 0 via intersection cuts.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_interminor.h"
#include "scip/expr.h"
#include "scip/expr_var.h"
#include "scip/expr_pow.h"
#include "scip/expr_product.h"
#include "scip/nlpi_ipopt.h"
#include "scip/cons_nonlinear.h"



#define SEPA_NAME              "interminor"
#define SEPA_DESC              "intersection cuts separator to ensure that 2x2 minors of X (= xx') have determinant 0"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MINCUTVIOL         1e-4 /**< default minimum required violation of a cut */
#define DEFAULT_RANDSEED            157 /**< default random seed */
#define DEFAULT_MAXROUNDS            10 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define BINSEARCH_MAXITERS          120 /**< default iteration limit for binary search */
#define DEFAULT_USESTRENGTHENING  FALSE /**< default for using strengthend intersection cuts to separate */
#define DEFAULT_USEBOUNDS         FALSE /**< default for using nonnegativity bounds when separating */

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_VAR**            minors;             /**< variables of 2x2 minors; each minor is stored like (auxvar_x^2,auxvar_y^2,auxvar_xy) */
   SCIP_Bool*            isdiagonal;         /**< bool array determining if the variables appearing in the minor are diagonal */
   int                   nminors;            /**< total number of minors */
   int                   minorssize;         /**< size of minors array */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   SCIP_Bool             detectedminors;     /**< has minor detection be called? */
   SCIP_Real             mincutviol;         /**< minimum required violation of a cut */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generation */
   SCIP_Bool             usestrengthening;   /**< whether to use strengthened intersection cuts to separate minors */
   SCIP_Bool             usebounds;          /**< whether to also enforce nonegativity bounds of principle minors */
};

/* these represent a row */
struct rowdata
{
   int*                  vals;               /**< index of the column */
   int                   rowidx;             /**< index corresponding to variable of that row */
   int                   nvals;              /**< number of nonzero entries in column */
   int                   valssize;           /**< size of the array that is currently allocated */
   SCIP_HASHMAP*         auxvars;            /**< entry of the matrix */
};

/*
 * Local methods
 */

/** helper method to store a 2x2 minor in the separation data */
static
SCIP_RETCODE sepadataAddMinor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR*             auxvarxik,          /**< auxiliary variable X_ik = x_i * x_k */
   SCIP_VAR*             auxvarxil,          /**< auxiliary variable X_il = x_i * x_l */
   SCIP_VAR*             auxvarxjk,          /**< auxiliary variable X_jk = x_j * x_k */
   SCIP_VAR*             auxvarxjl,          /**< auxiliary variable X_jl = x_j * x_l */
   SCIP_Bool             isauxvarxikdiag,    /**< is X_ik diagonal? (i.e. i = k) */
   SCIP_Bool             isauxvarxildiag,    /**< is X_il diagonal? (i.e. i = l) */
   SCIP_Bool             isauxvarxjkdiag,    /**< is X_jk diagonal? (i.e. j = k) */
   SCIP_Bool             isauxvarxjldiag     /**< is X_jl diagonal? (i.e. j = l) */
   )
{
   assert(sepadata != NULL);
   assert(auxvarxik != NULL);
   assert(auxvarxil != NULL);
   assert(auxvarxjk != NULL);
   assert(auxvarxjl != NULL);
   assert(auxvarxik != auxvarxil);
   assert(auxvarxjk != auxvarxjl);

   SCIPdebugMsg(scip, "store 2x2 minor: [%s %s, %s %s]\n", SCIPvarGetName(auxvarxik), SCIPvarGetName(auxvarxil),
         SCIPvarGetName(auxvarxjk), SCIPvarGetName(auxvarxjl));

   /* reallocate if necessary */
   if( sepadata->minorssize < 4 * (sepadata->nminors + 1) )
   {
      int newsize = SCIPcalcMemGrowSize(scip, 4 * (sepadata->nminors + 1));
      assert(newsize >= 4 * (sepadata->nminors + 1));

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->minors), sepadata->minorssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->isdiagonal), sepadata->minorssize, newsize) );
      sepadata->minorssize = newsize;
   }

   /* store minor */
   sepadata->minors[4 * sepadata->nminors] = auxvarxik;
   sepadata->minors[4 * sepadata->nminors + 1] = auxvarxil;
   sepadata->minors[4 * sepadata->nminors + 2] = auxvarxjk;
   sepadata->minors[4 * sepadata->nminors + 3] = auxvarxjl;
   sepadata->isdiagonal[4 * sepadata->nminors] = isauxvarxikdiag;
   sepadata->isdiagonal[4 * sepadata->nminors + 1] = isauxvarxildiag;
   sepadata->isdiagonal[4 * sepadata->nminors + 2] = isauxvarxjkdiag;
   sepadata->isdiagonal[4 * sepadata->nminors + 3] = isauxvarxjldiag;
   ++(sepadata->nminors);

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, auxvarxik) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvarxil) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvarxjk) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvarxjl) );

   return SCIP_OKAY;
}

/** helper method to clear separation data */
static
SCIP_RETCODE sepadataClear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   int i;

   assert(sepadata != NULL);

   SCIPdebugMsg(scip, "clear separation data\n");

   /* release captured variables */
   for( i = 0; i < 4 * sepadata->nminors; ++i )
   {
      assert(sepadata->minors[i] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &sepadata->minors[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &sepadata->minors, sepadata->minorssize);
   SCIPfreeBlockMemoryArrayNull(scip, &sepadata->isdiagonal, sepadata->minorssize);

   /* reset counters */
   sepadata->nminors = 0;
   sepadata->minorssize = 0;

   return SCIP_OKAY;
}

/** helper method to get the variables associated to a minor */
static
SCIP_RETCODE getMinorVars(
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   int                   idx,                /**< index of the stored minor */
   SCIP_VAR**            auxvarxik,          /**< auxiliary variable X_ik = x_i * x_k */
   SCIP_VAR**            auxvarxil,          /**< auxiliary variable X_il = x_i * x_l */
   SCIP_VAR**            auxvarxjk,          /**< auxiliary variable X_jk = x_j * x_k */
   SCIP_VAR**            auxvarxjl,          /**< auxiliary variable X_jl = x_j * x_l */
   SCIP_Bool*            isauxvarxikdiag,    /**< is X_ik diagonal? (i.e. i = k) */
   SCIP_Bool*            isauxvarxildiag,    /**< is X_il diagonal? (i.e. i = l) */
   SCIP_Bool*            isauxvarxjkdiag,    /**< is X_jk diagonal? (i.e. j = k) */
   SCIP_Bool*            isauxvarxjldiag     /**< is X_jl diagonal? (i.e. j = l) */
   )
{
   assert(auxvarxik != NULL);
   assert(auxvarxil != NULL);
   assert(auxvarxjk != NULL);
   assert(auxvarxjl != NULL);

   *auxvarxik = sepadata->minors[4 * idx];
   *auxvarxil = sepadata->minors[4 * idx + 1];
   *auxvarxjk = sepadata->minors[4 * idx + 2];
   *auxvarxjl = sepadata->minors[4 * idx + 3];

   *isauxvarxikdiag = sepadata->isdiagonal[4 * idx];
   *isauxvarxildiag = sepadata->isdiagonal[4 * idx + 1];
   *isauxvarxjkdiag = sepadata->isdiagonal[4 * idx + 2];
   *isauxvarxjldiag = sepadata->isdiagonal[4 * idx + 3];

   return SCIP_OKAY;
}


/** adds a new entry (i.e., auxvar) of in (row, col) of matrix M.
 *
 * we have a matrix, M, indexed by the variables
 * M(xi, xk) is the auxiliary variable of xi * xk if it exists
 * We store, for each row of the matrix, the indices of the nonzero column entries (assoc with the given row) and the auxiliary variable for xi * xk
 * The nonzero column entries are stored as an array (struct rowdata)
 * So we have a hashmap mapping each variable (row of the matrix) with its array representing the nonzero entries of the row.
 */
static
SCIP_RETCODE insertIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         rowmap,             /**< hashmap of the rows of the matrix */
   SCIP_VAR*             row,                /**< variable corresponding to row of new entry */
   SCIP_VAR*             col,                /**< variable corresponding to column of new entry */
   SCIP_VAR*             auxvar,             /**< auxvar to insert into the matrix */
   int*                  rowindices,         /**< array of indices of all variables corresponding to a row */
   int*                  nrows               /**< number of rows */
   )
{
   SCIPdebugMsg(scip, "inserting %s in row %s and col %s \n", SCIPvarGetName(auxvar), SCIPvarGetName(row), SCIPvarGetName(col));

   /* check whether variable has an array associated to it */
   if( SCIPhashmapExists(rowmap, (void*)row) )
   {
      struct rowdata* arr;

      arr = (struct rowdata*)SCIPhashmapGetImage(rowmap, (void *)row);

      /* reallocate if necessary */
      if( arr->valssize < arr->nvals + 1 )
      {
         int newsize = SCIPcalcMemGrowSize(scip, arr->nvals + 1);
         assert(newsize > arr->nvals + 1);

         SCIP_CALL( SCIPreallocBufferArray(scip, &(arr->vals), newsize) );
         arr->valssize = newsize;
      }

      /* insert */
      arr->vals[arr->nvals] = SCIPvarGetProbindex(col);
      SCIP_CALL( SCIPhashmapInsert(arr->auxvars, (void*)col, (void *)auxvar) );
      arr->nvals += 1;
   }
   else
   {
      struct rowdata* arr;

      /* create index array */
      SCIP_CALL( SCIPallocBuffer(scip, &arr) );
      arr->valssize = 10;
      arr->nvals = 0;
      SCIP_CALL( SCIPallocBufferArray(scip, &arr->vals, arr->valssize) );
      SCIP_CALL( SCIPhashmapCreate(&arr->auxvars, SCIPblkmem(scip), arr->valssize) );

      /* insert */
      arr->rowidx = SCIPvarGetProbindex(row);
      arr->vals[arr->nvals] = SCIPvarGetProbindex(col);
      SCIP_CALL( SCIPhashmapInsert(arr->auxvars, (void*)col, (void *)auxvar) );
      arr->nvals += 1;

      /* store in hashmap */
      SCIP_CALL( SCIPhashmapInsert(rowmap, (void*)row, (void *)arr) );

      /* remember the new row */
      rowindices[*nrows] = SCIPvarGetProbindex(row);
      *nrows += 1;
   }

   return SCIP_OKAY;
}

/** method to detect and store principal minors */
static
SCIP_RETCODE detectMinors(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPRITER* it;
   SCIP_HASHMAP* rowmap;
   int* rowvars = NULL;
   int* intersection;
   int nrowvars = 0;
   int c;
   int i;

#ifdef SCIP_STATISTIC
   SCIP_Real totaltime = -SCIPgetTotalTime(scip);
#endif

   assert(sepadata != NULL);

   /* check whether minor detection has been called already */
   if( sepadata->detectedminors )
      return SCIP_OKAY;

   assert(sepadata->minors == NULL);
   assert(sepadata->nminors == 0);

   /* we assume that the auxiliary variables in the nonlinear constraint handler have been already generated */
   sepadata->detectedminors = TRUE;

   /* check whether there are nonlinear constraints available */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   if( conshdlr == NULL || SCIPconshdlrGetNConss(conshdlr) == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "call detectMinors()\n");

   /* allocate memory */
   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPhashmapCreate(&rowmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowvars, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &intersection, SCIPgetNVars(scip)) );

   /* initialize iterator */
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   for( c = 0; c < SCIPconshdlrGetNConss(conshdlr); ++c )
   {
      SCIP_CONS* cons;
      SCIP_EXPR* expr;
      SCIP_EXPR* root;

      cons = SCIPconshdlrGetConss(conshdlr)[c];
      assert(cons != NULL);
      root = SCIPgetExprNonlinear(cons);
      assert(root != NULL);

      for( expr = SCIPexpriterRestartDFS(it, root); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*//*lint !e440*/
      {
         SCIP_EXPR** children;
         SCIP_VAR* auxvar;

         SCIPdebugMsg(scip, "visit expression %p in constraint %s\n", (void*)expr, SCIPconsGetName(cons));

         /* check whether the expression has an auxiliary variable */
         auxvar = SCIPgetExprAuxVarNonlinear(expr);
         if( auxvar == NULL )
         {
            SCIPdebugMsg(scip, "expression has no auxiliary variable -> skip\n");
            continue;
         }

         children = SCIPexprGetChildren(expr);

         /* check for expr = (x)^2 */
         if( SCIPexprGetNChildren(expr) == 1 && SCIPisExprPower(scip, expr)
            && SCIPgetExponentExprPow(expr) == 2.0
            && SCIPgetExprAuxVarNonlinear(children[0]) != NULL )
         {
            SCIP_VAR* quadvar;

            assert(children[0] != NULL);

            quadvar = SCIPgetExprAuxVarNonlinear(children[0]);
            assert(quadvar != NULL);

            SCIP_CALL( insertIndex(scip, rowmap, quadvar, quadvar, auxvar, rowvars, &nrowvars) );
         }
         /* check for expr = x_i * x_k */
         else if( SCIPexprGetNChildren(expr) == 2 && SCIPisExprProduct(scip, expr)
            && SCIPgetExprAuxVarNonlinear(children[0]) != NULL && SCIPgetExprAuxVarNonlinear(children[1]) != NULL )
         {
            SCIP_VAR* xi;
            SCIP_VAR* xk;

            assert(children[0] != NULL);
            assert(children[1] != NULL);

            xi = SCIPgetExprAuxVarNonlinear(children[0]);
            xk = SCIPgetExprAuxVarNonlinear(children[1]);

            SCIP_CALL( insertIndex(scip, rowmap, xk, xi, auxvar, rowvars, &nrowvars) );
            SCIP_CALL( insertIndex(scip, rowmap, xi, xk, auxvar, rowvars, &nrowvars) );
         }
      }
   }

   /* sort the column entries */
   for( i = 0; i < nrowvars; ++i )
   {
      struct rowdata* row;

      row = (struct rowdata*)SCIPhashmapGetImage(rowmap, (void *)SCIPgetVars(scip)[rowvars[i]]);
      SCIPsortInt(row->vals, row->nvals);
   }

   /* store 2x2 minors */
   /* TODO: we might store some minors twice since the matrix is symmetric. Handle that! (see unit test for example) */
   for( i = 0; i < nrowvars; ++i )
   {
      int j;
      struct rowdata* rowi;

      rowi = (struct rowdata*)SCIPhashmapGetImage(rowmap, (void *)SCIPgetVars(scip)[rowvars[i]]);

      for( j = i + 1; j < nrowvars; ++j )
      {
         struct rowdata* rowj;
         int ninter;

         rowj = (struct rowdata*)SCIPhashmapGetImage(rowmap, (void *)SCIPgetVars(scip)[rowvars[j]]);

         SCIPcomputeArraysIntersectionInt(rowi->vals, rowi->nvals, rowj->vals, rowj->nvals, intersection, &ninter);

         if( ninter > 1)
         {
            int p;

            for( p = 0; p < ninter - 1; ++p )
            {
               int q;

               for( q = p + 1; q < ninter; ++q )
               {
                  SCIP_HASHMAP* rowicols;
                  SCIP_HASHMAP* rowjcols;
                  SCIP_VAR* colk;
                  SCIP_VAR* coll;
                  SCIP_VAR* auxvarik;
                  SCIP_VAR* auxvaril;
                  SCIP_VAR* auxvarjk;
                  SCIP_VAR* auxvarjl;
                  int ii;
                  int jj;
                  int k;
                  int l;
                  SCIP_Bool isauxvarikdiag = FALSE;
                  SCIP_Bool isauxvarildiag = FALSE;
                  SCIP_Bool isauxvarjkdiag = FALSE;
                  SCIP_Bool isauxvarjldiag = FALSE;

                  ii = rowi->rowidx;
                  jj = rowj->rowidx;
                  k = intersection[p];
                  l = intersection[q];

                  rowicols = rowi->auxvars;
                  rowjcols = rowj->auxvars;

                  colk = SCIPgetVars(scip)[k];
                  coll = SCIPgetVars(scip)[l];

                  auxvarik = (SCIP_VAR*) SCIPhashmapGetImage(rowicols, colk);
                  auxvaril = (SCIP_VAR*) SCIPhashmapGetImage(rowicols, coll);
                  auxvarjk = (SCIP_VAR*) SCIPhashmapGetImage(rowjcols, colk);
                  auxvarjl = (SCIP_VAR*) SCIPhashmapGetImage(rowjcols, coll);

                  if( ii == k )
                     isauxvarikdiag = TRUE;
                  else if( ii == l )
                     isauxvarildiag = TRUE;
                  if( jj == k )
                     isauxvarjkdiag = TRUE;
                  else if( jj == l )
                     isauxvarjldiag = TRUE;

                  SCIP_CALL( sepadataAddMinor(scip, sepadata, auxvarik, auxvaril, auxvarjk, auxvarjl,
                        isauxvarikdiag, isauxvarildiag, isauxvarjkdiag, isauxvarjldiag) );
               }
            }
         }
      }
      SCIPfreeBufferArrayNull(scip, &rowi->vals);
      SCIPhashmapFree(&rowi->auxvars);
      SCIPfreeBufferArrayNull(scip, &rowi);
   }


   SCIPdebugMsg(scip, "found %d principal minors in total\n", sepadata->nminors);

   /* free memory */
   SCIPfreeBufferArray(scip, &intersection);
   SCIPfreeBufferArray(scip, &rowvars);
   SCIPhashmapFree(&rowmap);
   SCIPfreeExpriter(&it);

#ifdef SCIP_STATISTIC
   totaltime += SCIPgetTotalTime(scip);
   SCIPstatisticMessage("MINOR DETECT %s %f %d %d\n", SCIPgetProbName(scip), totaltime, sepadata->nminors, maxminors);
#endif

   return SCIP_OKAY;
}

/** constructs map between lp position of a basic variable and its row in the tableau */
static
SCIP_RETCODE constructBasicVars2TableauRowMap(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  map                 /**< buffer to store the map */
   )
{
   int* basisind;
   int nrows;
   int i;

   nrows = SCIPgetNLPRows(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );

   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );
   for( i = 0; i < nrows; ++i )
   {
      if( basisind[i] >= 0 )
         map[basisind[i]] = i;
   }

   SCIPfreeBufferArray(scip, &basisind);

   return SCIP_OKAY;
}

/** The restriction of the function representing the maximal S-free set to zlp + t * ray has the form
 * SQRT(A t^2 + B t + C) - (D t + E).
 * This function computes the coefficients A, B, C, D, E for the given ray.
 */
static
SCIP_RETCODE computeRestrictionToRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            ray,                /**< coefficients of ray */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real*            coefs,              /**< buffer to store A, B, C, D, and E of cases 1, 2, 3, or 4a*/
   SCIP_Real*            coefs4b,            /**< buffer to store A, B, C, D, and E of case 4b */
   SCIP_Real*            coefscondition,     /**< buffer to store coefs for checking whether we are in case 4a or 4b */
   SCIP_Bool             usebounds,          /**< TRUE if we want to separate non-negative bound */
   SCIP_Real*            ad,                 /**< coefs a and d for the hyperplane aTx + dTy <= 0 */
   SCIP_Bool*            success             /**< FALSE if we need to abort generation because of numerics */
   )
{
   SCIP_Real eigenvectors[16] = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
   SCIP_Real eigenvalues[4] = {0.5, 0.5, -0.5, -0.5};
   SCIP_Real eigencoef = 0.7071067811865475244008443621048490;
   SCIP_Real* a;
   SCIP_Real* b;
   SCIP_Real* c;
   SCIP_Real* d;
   SCIP_Real* e;
   SCIP_Real min;
   SCIP_Real max;
   SCIP_Real norm1;
   SCIP_Real norm2;
   int negidx;
   int posidx;
   int i;

   *success = TRUE;

   /* set all coefficients to zero */
   memset(coefs, 0, 5 * sizeof(SCIP_Real));
   memset(coefs4b, 0, 5 * sizeof(SCIP_Real));
   norm1 = 0.0;
   norm2 = 0.0;

   a = coefs;
   b = coefs + 1;
   c = coefs + 2;
   d = coefs + 3;
   e = coefs + 4;

   negidx = 2;
   posidx = 0;
   for( i = 0; i < 4; ++i )
   {
      int j;
      SCIP_Real vzlp;
      SCIP_Real vdotray;

      vzlp = 0;
      vdotray = 0;

      /* compute eigenvec * ray and eigenvec * solution */
      for( j = 0; j < 4; ++j )
      {
         vdotray += eigencoef * eigenvectors[4 * i + j] * ray[j];
         vzlp += eigencoef * eigenvectors[4 * i + j] * SCIPvarGetLPSol(vars[j]);
      }

      if( eigenvalues[i] > 0 )
      {
         /* positive eigenvalue: compute D and E */
         *d += eigenvalues[i] * vzlp * vdotray;
         *e += eigenvalues[i] * SQR( vzlp );

         if( usebounds )
         {
            norm1 += eigenvalues[i] * (1 - SQR( ad[posidx] )) * SQR( vzlp );
            norm2 += SQRT( eigenvalues[i] ) * ad[posidx] * vzlp;
            ++posidx;
         }

      }
      else
      {
         /* negative eigenvalue: compute A, B, and C */
         *a -= eigenvalues[i] * SQR( vdotray );
         *b -= 2.0 * eigenvalues[i] * vzlp * vdotray;
         *c -= eigenvalues[i] * SQR( vzlp );

         if( usebounds )
         {
            coefs4b[0] -= eigenvalues[i] * (1 - SQR( ad[negidx] )) * SQR( vdotray );
            coefs4b[1] -= 2.0 * eigenvalues[i] * (1 - SQR( ad[negidx] )) * vzlp * vdotray;
            coefs4b[2] -= eigenvalues[i] * (1 - SQR( ad[negidx] )) * SQR( vzlp );
            coefs4b[3] += SQRT( -eigenvalues[i] ) * ad[negidx] * vdotray;
            coefs4b[4] += SQRT( -eigenvalues[i] ) * ad[negidx] * vzlp;
            ++negidx;
         }
      }
   }

   assert(*e > 0);

   if( SQRT( *c ) - SQRT( *e ) >= 0.0 )
   {
      assert(SQRT( *c ) - SQRT( *e ) < 1e-6);
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* finish computation of coefficients when using bounds */
   if( usebounds )
   {
      coefscondition[0] = norm2 / SQRT( *e );
      coefscondition[1] = coefs4b[3];
      coefscondition[2] = coefs4b[4];

      coefs4b[0] *= norm1 / *e;
      coefs4b[1] *= norm1 / *e;
      coefs4b[2] *= norm1 / *e;
      coefs4b[3] *= norm2 / SQRT( *e );
      coefs4b[4] *= norm2 / SQRT( *e );

      coefs4b[3] += *d / SQRT( *e );
      coefs4b[4] += SQRT( *e );

      assert( SQRT( coefs4b[2] ) - coefs4b[4] < 0.0 );
   }

   /* finish computation of D and E */
   *e = SQRT( *e );
   *d /= *e;

   /* maybe we want to avoid a large dynamism between A, B and C */
   max = 0.0;
   min = SCIPinfinity(scip);
   for( i = 0; i < 3; ++i )
   {
      SCIP_Real absval;

      absval = ABS(coefs[i]);
      if( max < absval )
         max = absval;
      if( absval != 0.0 && absval < min )
         min = absval;
   }

   if( SCIPisHugeValue(scip, max / min) )
   {
#ifdef DEBUG_INTERSECTIONCUT
      printf("Bad numerics: max(A,B,C)/min(A,B,C) is too large (%g)\n", max / min);
#endif
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* some sanity checks */
   assert(*c >= 0); /* radicand at zero */
   assert(SQRT( *c ) - *e < 0); /* the function at 0 must be negative */
   assert(*a >= 0); /* the function inside the root is convex */

#ifdef DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "Restriction yields: a,b,c,d,e %g %g %g %g %g\n", coefs[0], coefs[1], coefs[2], coefs[3], coefs[4]);
#endif

   return SCIP_OKAY;
}

/** returns phi(zlp + t * ray) = SQRT(A t^2 + B t + C) - (D t + E) */  /*lint -e{715}*/
static
SCIP_Real evalPhiAtRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             t,                  /**< argument of phi restricted to ray */
   SCIP_Real             a,                  /**< value of A */
   SCIP_Real             b,                  /**< value of B */
   SCIP_Real             c,                  /**< value of C */
   SCIP_Real             d,                  /**< value of D */
   SCIP_Real             e                   /**< value of E */
   )
{
#ifdef INTERCUTS_DBLDBL
   SCIP_Real QUAD(lin);
   SCIP_Real QUAD(disc);
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(root);

   /* d * t + e */
   SCIPquadprecProdDD(lin, d, t);
   SCIPquadprecSumQD(lin, lin, e);

   /* a * t * t */
   SCIPquadprecSquareD(disc, t);
   SCIPquadprecProdQD(disc, disc, a);

   /* b * t */
   SCIPquadprecProdDD(tmp, b, t);

   /* a * t * t + b * t */
   SCIPquadprecSumQQ(disc, disc, tmp);

   /* a * t * t + b * t + c */
   SCIPquadprecSumQD(disc, disc, c);

   /* sqrt(above): can't take sqrt of 0! */
   if( QUAD_TO_DBL(disc) == 0 )
   {
      QUAD_ASSIGN(root, 0.0);
   }
   else
   {
      SCIPquadprecSqrtQ(root, disc);
   }

   /* final result */
   QUAD_SCALE(lin, -1.0);
   SCIPquadprecSumQQ(tmp, root, lin);

   assert(!SCIPisInfinity(scip, t) || QUAD_TO_DBL(tmp) <= 0);

   return  QUAD_TO_DBL(tmp);
#else
   return SQRT( a * t * t + b * t + c ) - ( d * t + e );
#endif
}

/** helper function of computeRoot: we want phi to be <= 0 */
static
void doBinarySearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             a,                  /**< value of A */
   SCIP_Real             b,                  /**< value of B */
   SCIP_Real             c,                  /**< value of C */
   SCIP_Real             d,                  /**< value of D */
   SCIP_Real             e,                  /**< value of E */
   SCIP_Real*            sol                 /**< buffer to store solution; also gives initial point */
   )
{
   SCIP_Real lb = 0.0;
   SCIP_Real ub = *sol;
   SCIP_Real curr;
   int i;

   for( i = 0; i < BINSEARCH_MAXITERS; ++i )
   {
      SCIP_Real phival;

      curr = (lb + ub) / 2.0;
      phival = evalPhiAtRay(scip, curr, a, b, c, d, e);
#ifdef INTERCUT_MOREDEBUG
      printf("%d: lb,ub %.10f, %.10f. curr = %g -> phi at curr %g -> phi at lb %g \n", i, lb, ub, curr, phival, evalPhiAtRay(scip, lb, a, b, c, d, e));
#endif

      if( phival <= 0.0 )
      {
         lb = curr;
         if( SCIPisFeasZero(scip, phival) || SCIPisFeasEQ(scip, ub, lb) )
            break;
      }
      else
         ub = curr;
   }

   *sol = lb;

}

/** checks if we are in case 4a, i.e., if
 * (num(xhat_{r+1}(zlp)) / E) * SQRT(A * tsol^2 + B * tsol + C) + w(ray) * tsol + num(yhat_{s+1}(zlp)) <= 0
 */
static
SCIP_Real isCase4a(
   SCIP_Real             tsol,               /**< t in the above formula */
   SCIP_Real*            coefs,              /**< coefficients A, B, C, D, and E of case 4a */
   SCIP_Real*            coefscondition      /**< extra coefficients needed for the evaluation of the condition:
                                                num(xhat_{r+1}(zlp)) / E; w(ray); num(yhat_{s+1}(zlp)) */
   )
{
   return (coefscondition[0] * SQRT( coefs[0] * SQR( tsol ) + coefs[1] * tsol + coefs[2] ) + coefscondition[1] *
      tsol + coefscondition[2]) <= 0.0;
}

/**  finds smallest positive root phi by finding the smallest positive root of
 * (A - D^2) t^2 + (B - 2 D*E) t + (C - E^2) = 0
 *
 * However, we are conservative and want a solution such that phi is negative, but close to 0;
 * thus we correct the result with a binary search
 */
static
SCIP_Real computeRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            coefs               /**< value of A */
   )
{
   SCIP_Real sol;
   SCIP_INTERVAL bounds;
   SCIP_INTERVAL result;
   SCIP_Real a = coefs[0];
   SCIP_Real b = coefs[1];
   SCIP_Real c = coefs[2];
   SCIP_Real d = coefs[3];
   SCIP_Real e = coefs[4];

   /* there is an intersection point if and only if SQRT(A) > D: here we are beliving in math, this might cause
    * numerical issues
    */
   if( SQRT( a ) <= d )
   {
      sol = SCIPinfinity(scip);

      return sol;
   }

   SCIPintervalSetBounds(&bounds, 0.0, SCIPinfinity(scip));

   /* SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar finds all x such that a x^2 + b x >= c and x in bounds.
    * it is known that if tsol is the root we are looking for, then gamma(zlp + t * ray) <= 0 between 0 and tsol, thus
    * tsol is the smallest t such that (A - D^2) t^2 + (B - 2 D*E) t + (C - E^2) >= 0
    */
   SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(SCIP_INTERVAL_INFINITY, &result, a - d * d, b - 2.0 * d *
         e, -(c - e * e), bounds);

   /* it can still be empty because of our infinity, I guess... */
   sol = SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, result) ? SCIPinfinity(scip) : SCIPintervalGetInf(result);

   /* check that solution is acceptable, ideally it should be <= 0, however when it is positive, we trigger a binary
    * search to make it negative. This binary search might return a solution point that is not at accurately 0 as the
    * one obtained from the function above. Thus, it might fail to satisfy the condition of case 4b in some cases, e.g.,
    * ex8_3_1, bchoco05, etc
    */
   if( evalPhiAtRay(scip, sol, a, b, c, d, e) <= 1e-10 )
   {
#ifdef INTERCUT_MOREDEBUG
      printf("interval solution returned %g -> phival = %g, believe it\n", sol, evalPhiAtRay(sol, a, b, c, d, e));
      printf("don't do bin search\n");
#endif

      return sol;
   }
   else
   {
      /* perform a binary search to make it negative: this might correct a wrong infinity (e.g. crudeoil_lee1_05) */
#ifdef INTERCUT_MOREDEBUG
      printf("do bin search because phival is %g\n", evalPhiAtRay(scip, sol, a, b, c, d, e));
#endif
      doBinarySearch(scip, a, b, c, d, e, &sol);
   }

   return sol;
}

/** The maximal S-free set is gamma(z) <= 0; we find the intersection point of the ray `ray` starting from zlp with the
 * boundary of the S-free set.
 * That is, we find t >= 0 such that gamma(zlp + t * ray) = 0.
 *
 * In cases 1,2, and 3, gamma is of the form
 *    gamma(zlp + t * ray) = SQRT(A t^2 + B t + C) - (D t + E)
 *
 * In the case 4 gamma is of the form
 *    gamma(zlp + t * ray) = SQRT(A t^2 + B t + C) - (D t + E)          if some condition holds
 *                           SQRT(A' t^2 + B' t + C') - (D' t + E')     otherwise
 *
 * It can be shown (given the special properties of gamma) that the smallest positive root of each function of the form
 * SQRT(a t^2 + b t + c) - (d t + e)
 * is the same as the smallest positive root of the quadratic equation:
 *       (SQRT(a t^2 + b t + c) - (d t + e)) * (SQRT(a t^2 + b t + c) + (d t + e)) = 0
 *  <==> (a - d^2) t^2 + (b - 2 d*e) t + (c - e^2) = 0
 *
 * So, in cases 1, 2, and 3, this function just returns the solution of the above equation.
 * In case 4, it first solves the equation assuming we are in the first piece.
 * If there is no solution, then the second piece can't have a solution (first piece >= second piece for all t)
 * Then we check if the solution satisfies the condition.
 * If it doesn't then we solve the equation for the second piece.
 * If it has a solution, then it _has_ to be the solution.
 */
static
SCIP_Real computeIntersectionPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             usebounds,          /**< whether we are in case 4 or not */
   SCIP_Real*            coefs,              /**< values of A, B, C, D, and E of cases 1, 2, 3, or 4a */
   SCIP_Real*            coefs4b,            /**< values of A, B, C, D, and E of case 4b */
   SCIP_Real*            coefscondition      /**< values needed to evaluate condition of case 4 */
   )
{
   SCIP_Real sol;
   SCIP_Real sol4b;

   assert(coefs != NULL);

   if( ! usebounds )
      return computeRoot(scip, coefs);

   assert(coefs4b != NULL);
   assert(coefscondition != NULL);

   /* compute solution of first piece */
   sol = computeRoot(scip, coefs);

   /* if there is no solution --> second piece doesn't have solution */
   if( SCIPisInfinity(scip, sol) )
   {
      /* this assert fails on multiplants_mtg5 the problem is that sqrt(A) <= D in 4a but not in 4b,
       * now, this is impossible since the phi4a >= phi4b, so actually sqrt(A) is 10e-15 away from
       * D in 4b
       */
      /* assert(SCIPisInfinity(scip, computeRoot(scip, coefs4b))); */
      return sol;
   }

   /* if solution of 4a is in 4a, then return */
   if( isCase4a(sol, coefs, coefscondition) )
      return sol;

   /* not on 4a --> then the intersection point is whatever 4b says: as phi4a >= phi4b, the solution of phi4b should
    * always be larger (but shouldn't be equal at this point given that isCase4a failed, and the condition function
    * evaluates to 0 when phi4a == phi4b) than the solution of phi4a; However, because of numerics (or limits in the
    * binary search) we can find a slightly smaller solution; thus, we just keep the larger one
    */
   sol4b = computeRoot(scip, coefs4b);

   return MAX(sol, sol4b);
}

/** adds cutcoef * (col - col*) to rowprep */
static
SCIP_RETCODE addColToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store intersection cut */
   SCIP_Real             cutcoef,            /**< cut coefficient */
   SCIP_COL*             col                 /**< column to add to rowprep */
   )
{
   assert(col != NULL);

#ifdef DEBUG_INTERCUTS_NUMERICS
   SCIPinfoMessage(scip, NULL, "adding col %s to cut. %g <= col <= %g\n", SCIPvarGetName(SCIPcolGetVar(col)),
      SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)));
   SCIPinfoMessage(scip, NULL, "col is active at %s. Value %.15f\n", SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER ? "lower bound" :
      "upper bound" , SCIPcolGetPrimsol(col));
#endif

   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPcolGetVar(col), cutcoef) );
   SCIProwprepAddConstant(rowprep, -cutcoef * SCIPcolGetPrimsol(col) );

   return SCIP_OKAY;
}

/** adds cutcoef * (slack - slack*) to rowprep
  *
  * row is lhs <= <coefs, vars> + constant <= rhs, thus slack is defined by
  * slack + <coefs, vars> + constant = side
  * If row (slack) is at upper, it means that <coefs,vars*> + constant = rhs, and so
  * slack* = side - rhs --> slack - slack* = rhs - <coefs, vars> - constant.
  * If row (slack) is at lower, then <coefs,vars*> + constant = lhs, and so
  * slack* = side - lhs --> slack - slack* = lhs - <coefs, vars> - constant.
  */
static
SCIP_RETCODE addRowToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store intersection cut */
   SCIP_Real             cutcoef,            /**< cut coefficient */
   SCIP_ROW*             row,                /**< row, whose slack we are adding to rowprep */
   SCIP_Bool*            success             /**< buffer to store whether the row is nonbasic enough */
   )
{
   int i;
   SCIP_COL** rowcols;
   SCIP_Real* rowcoefs;
   int nnonz;

   assert(row != NULL);

   rowcols = SCIProwGetCols(row);
   rowcoefs = SCIProwGetVals(row);
   nnonz = SCIProwGetNLPNonz(row);

#ifdef DEBUG_INTERCUTS_NUMERICS
   SCIPinfoMessage(scip, NULL, "adding slack var row_%d to cut. %g <= row <= %g\n", SCIProwGetLPPos(row), SCIProwGetLhs(row), SCIProwGetRhs(row));
   SCIPinfoMessage(scip, NULL, "row is active at %s = %.15f Activity %.15f\n", SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER ? "lhs" :
   "rhs" , SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER ? SCIProwGetLhs(row) : SCIProwGetRhs(row),
   SCIPgetRowActivity(scip, row));
#endif

   if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER )
   {
      assert(!SCIPisInfinity(scip, -SCIProwGetLhs(row)));
      if( ! SCIPisFeasEQ(scip, SCIProwGetLhs(row), SCIPgetRowActivity(scip, row)) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      SCIProwprepAddConstant(rowprep, SCIProwGetLhs(row) * cutcoef);
   }
   else
   {
      assert(!SCIPisInfinity(scip, SCIProwGetRhs(row)));
      if( ! SCIPisFeasEQ(scip, SCIProwGetRhs(row), SCIPgetRowActivity(scip, row)) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      SCIProwprepAddConstant(rowprep, SCIProwGetRhs(row) * cutcoef);
   }

   for( i = 0; i < nnonz; i++ )
   {
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPcolGetVar(rowcols[i]), -rowcoefs[i] * cutcoef) );
   }

   SCIProwprepAddConstant(rowprep, -SCIProwGetConstant(row) * cutcoef);

   return SCIP_OKAY;
}

/** get the tableau rows of the variables in vars */
static
SCIP_RETCODE getTableauRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in the minor */
   int*                  basicvarpos2tableaurow,/**< map between basic var and its tableau row */
   SCIP_HASHMAP*         tableau,            /**< map between var an its tableau row */
   SCIP_Real**           tableaurows,        /**< buffer to store tableau row */
   SCIP_Bool*            success             /**< set to TRUE if no variable had basisstat = ZERO */
   )
{
   int v;
   int nrows;
   int ncols;

   *success = TRUE;

   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   /* check if we have the tableau row of the variable and if not compute it */
   for( v = 0; v < 4; ++v )
   {
      if( ! SCIPhashmapExists(tableau, (void*)vars[v]) )
      {
         SCIP_COL* col;

         /* get column of variable */
         col = SCIPvarGetCol(vars[v]);

         /* if variable is basic, then get its tableau row and insert it in the hashmap */
         if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
         {
            int lppos;
            SCIP_Real* densetableaurow;

            lppos = SCIPcolGetLPPos(col);
            SCIP_CALL( SCIPallocBufferArray(scip, &densetableaurow, ncols + nrows) );

            SCIP_CALL( SCIPgetLPBInvRow(scip, basicvarpos2tableaurow[lppos], &densetableaurow[ncols], NULL, NULL) );
            SCIP_CALL( SCIPgetLPBInvARow(scip, basicvarpos2tableaurow[lppos], &densetableaurow[ncols], densetableaurow, NULL, NULL) );

            /* insert tableau row in hashmap*/
            SCIP_CALL( SCIPhashmapInsert(tableau, (void*)vars[v], (void *)densetableaurow) );

         }
         else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
         {
            *success = FALSE;
            return SCIP_OKAY; /* don't even bother */
         }
         else
         {
            SCIP_CALL( SCIPhashmapInsert(tableau, (void*)vars[v], (void *)NULL) );
         }

      }

      /* get tableau row of var */
      tableaurows[v] = (SCIP_Real *)SCIPhashmapGetImage(tableau, (void*)vars[v]);
   }
   return SCIP_OKAY;
}

/** computes the cut coefs of the  non-basic (non-slack) variables (correspond to cols) and adds them to the
 * intersection cut
 */
static
SCIP_RETCODE addCols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real**           tableaurows,        /**< tableau rows corresponding to the variables in vars */
   SCIP_ROWPREP*         rowprep,            /**< store cut */
   SCIP_Real*            rays,               /**< buffer to store rays */
   int*                  nrays,              /**< pointer to store number of nonzero rays */
   int*                  rayslppos,          /**< buffer to store lppos of nonzero rays */
   SCIP_Real*            interpoints,        /**< buffer to store intersection points or NULL if not needed */
   SCIP_Bool             usebounds,          /**< TRUE if we want to separate non-negative bound */
   SCIP_Real*            ad,                 /**< coefs a and d for the hyperplane aTx + dTy <= 0 */
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   int i;
   int ncols;
   SCIP_COL** cols;

   *success = TRUE;

   /* loop over non-basic (non-slack) variables */
   cols = SCIPgetLPCols(scip);
   ncols = SCIPgetNLPCols(scip);
   for( i = 0; i < ncols; ++i )
   {
      SCIP_COL* col;
      SCIP_Real coefs[5];
      SCIP_Real coefs4b[5];
      SCIP_Real coefscondition[3];
      SCIP_Real factor;
      SCIP_Bool israynonzero;
      SCIP_Real cutcoef;
      SCIP_Real interpoint;
      int v;

      col = cols[i];

      /* set factor to store entries of ray as = [-BinvL, BinvU] */
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER )
         factor = -1.0;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER )
         factor = 1.0;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }
      else
         continue;

      /* build the ray */
      israynonzero = FALSE;
      for( v = 0; v < 4; ++v )
      {
         if( tableaurows[v] != NULL )
            rays[(*nrays) * 4 + v] = factor * (SCIPisZero(scip, tableaurows[v][i]) ? 0.0 : tableaurows[v][i]);
         else
         {
            if( col == SCIPvarGetCol(vars[v]) )
               rays[(*nrays) * 4 + v] = -factor;
            else
               rays[(*nrays) * 4 + v] = 0.0;
         }

         israynonzero = israynonzero || (rays[(*nrays) * 4 + v] != 0.0);
      }

      /* do nothing if ray is 0 */
      if( ! israynonzero )
         continue;

      /* compute the cut */
      SCIP_CALL( computeRestrictionToRay(scip, &rays[(*nrays) * 4], vars, coefs, coefs4b, coefscondition, usebounds,
            ad, success) );

      if( *success == FALSE )
         return SCIP_OKAY;

      /* compute intersection point */
      interpoint = computeIntersectionPoint(scip, usebounds, coefs, coefs4b, coefscondition);

      /* store intersection points */
      interpoints[*nrays] = interpoint;

      /* remember lppos */
      rayslppos[*nrays] = i;

      /* count nonzero rays */
      *nrays += 1;

      /* compute cut coef */
      cutcoef = SCIPisInfinity(scip, interpoint) ? 0.0 : 1.0 / interpoint;

      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      assert(SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER || SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER);
      SCIP_CALL( addColToCut(scip, rowprep, SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER ? -cutcoef :
            cutcoef, col) );
   }

   return SCIP_OKAY;
}

/** computes the cut coefs of the non-basic slack variables (correspond to rows) and adds them to the
 * intersection cut
 */
static
SCIP_RETCODE addRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real**           tableaurows,        /**< tableau rows corresponding to the variables in vars */
   SCIP_ROWPREP*         rowprep,            /**< store cut */
   SCIP_Real*            rays,               /**< buffer to store rays */
   int*                  nrays,              /**< pointer to store number of nonzero rays */
   int*                  rayslppos,          /**< buffer to store lppos of nonzero rays */
   SCIP_Real*            interpoints,        /**< buffer to store intersection points or NULL if not needed */
   SCIP_Bool             usebounds,          /**< TRUE if we want to separate non-negative bound */
   SCIP_Real*            ad,                 /**< coefs a and d for the hyperplane aTx + dTy <= 0 */
   SCIP_Bool*            success             /**< pointer to store whether the generation of cutcoefs was successful */
   )
{
   int i;
   int nrows;
   int ncols;
   SCIP_ROW** rows;

   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   *success = TRUE;

   /* loop over non-basic slack variables */
   rows = SCIPgetLPRows(scip);
   for( i = 0; i < nrows; ++i )
   {
      SCIP_ROW* row;
      SCIP_Real coefs[5];
      SCIP_Real coefs4b[5];
      SCIP_Real coefscondition[3];
      SCIP_Real factor;
      SCIP_Bool israynonzero;
      SCIP_Real cutcoef;
      SCIP_Real interpoint;
      int v;

      row = rows[i];

      /* set factor to store entries of ray as = [BinvL, -BinvU] */
      if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER )
         factor = 1.0;
      else if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER )
         factor = -1.0;
      else if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_ZERO )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }
      else
         continue;

      /* build the ray */
      israynonzero = FALSE;
      for( v = 0; v < 4; ++v )
      {
         int idx;

         idx = ncols + i;

         if( tableaurows[v] != NULL )
            rays[(*nrays) * 4 + v] = factor * (SCIPisZero(scip, tableaurows[v][idx]) ? 0.0 : tableaurows[v][idx]);
         else
         {
            /* TODO: We assume that slack variables can never occure in the minor. This is correct, right? */
            rays[(*nrays) * 4 + v] = 0.0;
         }

         israynonzero = israynonzero || (rays[(*nrays) * 4 + v] != 0.0);
      }

      /* do nothing if ray is 0 */
      if( ! israynonzero )
         continue;

      /* compute the cut */
      SCIP_CALL( computeRestrictionToRay(scip, &rays[(*nrays) * 4], vars, coefs, coefs4b, coefscondition, usebounds,
            ad, success) );

      if( *success == FALSE )
         return SCIP_OKAY;

      /* compute intersection point */
      interpoint = computeIntersectionPoint(scip, usebounds, coefs, coefs4b, coefscondition);

      /* store intersection points */
      interpoints[*nrays] = interpoint;

      /* store lppos of ray, make it negative so we can differentiate between cols and rows */
      rayslppos[*nrays] = -i - 1;

      /* count nonzero rays */
      *nrays += 1;

      /* compute cut coef */
      cutcoef = SCIPisInfinity(scip, interpoint) ? 0.0 : 1.0 / interpoint;

      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      assert(SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER || SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER);

      SCIP_CALL( addRowToCut(scip, rowprep, SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER ? cutcoef :
            -cutcoef, row, success) ); /* rows have flipper base status! */
   }

   return SCIP_OKAY;
}

/* checks if two rays are linearly dependent */
static
SCIP_Bool raysAreDependent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            ray1,               /**< coefficients of ray 1 */
   SCIP_Real*            ray2,               /**< coefficients of ray 2 */
   SCIP_Real*            coef                /**< pointer to store coef (s.t. r1 = coef * r2) in case rays are
                                                  dependent */
   )
{
   int i;

   *coef = 0.0;

   for( i = 0; i < 4; ++i )
   {
      /* rays cannot be dependent if one ray has zero entry and the other one doesn't */
      if( (SCIPisZero(scip, ray1[i]) && ! SCIPisZero(scip, ray2[i])) ||
         (! SCIPisZero(scip, ray1[i]) && SCIPisZero(scip, ray2[i])) )
      {
         return FALSE;
      }

      if( *coef != 0.0 )
      {
         /* cannot be dependent if the coefs aren't equal for all entries */
         if( ! SCIPisFeasEQ(scip, *coef, ray1[i] / ray2[i]) )
            return FALSE;
      }
      else
         *coef = ray1[i] / ray2[i];
   }

   return TRUE;
}

/** finds the smallest negative steplength for the current ray r_idx such that the combination
 * of r_idx with all rays not in the recession cone is in the recession cone
 */
static
SCIP_RETCODE findRho(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            rays,               /**< rays */
   int                   nrays,              /**< number of nonzero rays */
   int                   idx,                /**< index of current ray we want to find rho for */
   SCIP_Real*            interpoints,        /**< intersection points of nonzero rays */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real*            rho,                /**< pointer to store the optimal rho */
   SCIP_Bool             usebounds,          /**< TRUE if we want to separate non-negative bound */
   SCIP_Real*            ad,                 /**< coefs a and d for the hyperplane aTx + dTy <= 0 */
   SCIP_Bool*            success             /**< TRUE if computation of rho was successful */
   )
{
   int i;

   *success = TRUE;

   /* go through all rays not in the recession cone and compute the largest negative steplength possible. The
    * smallest of them is then the steplength rho we use for the current ray */
   *rho = 0;
   for( i = 0; i < nrays; ++i )
   {
      SCIP_Real currentrho;
      SCIP_Real coef;

      if( SCIPisInfinity(scip, interpoints[i]) )
         continue;

      /* if the rays are linearly independent, we don't need to search for rho */
      if( raysAreDependent(scip, &rays[4 * i], &rays[4 * idx], &coef) )
         currentrho = coef * interpoints[i];
      else
      {
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real alpha;
         int j;

         /* do binary search by lookig at the convex combinations of r_i and r_j */
         lb = 0.0;
         ub = 1.0;

         for( j = 0; j < BINSEARCH_MAXITERS; ++j )
         {
            SCIP_Real coefs[5];
            SCIP_Real coefs4b[5];
            SCIP_Real coefscondition[3];
            SCIP_Real newray[4];
            SCIP_Real interpoint;
            int k;

            alpha = (lb + ub) / 2.0;

            /* build the ray alpha * ray_i + (1 - alpha) * ray_idx */
            for( k = 0; k < 4; ++k )
               newray[k] = alpha * rays[4 * i + k] - (1 - alpha) * rays[4 * idx + k];

            /* restrict phi to the "new" ray */
            SCIP_CALL( computeRestrictionToRay(scip, newray, vars, coefs, coefs4b, coefscondition, usebounds,
                  ad, success) );

            if( ! *success )
               return SCIP_OKAY;

            /* check if restriction to "new" ray is numerically nasty. If so, treat the corresponding rho as if phi is
             * positive
             */

            /* compute intersection point */
            interpoint = computeIntersectionPoint(scip, usebounds, coefs, coefs4b, coefscondition);

            /* no root exists */
            if( SCIPisInfinity(scip, interpoint) )
            {
               lb = alpha;
               if( SCIPisEQ(scip, ub, lb) )
                  break;
            }
            else
               ub = alpha;
         }

         /* now we found the best convex combination which we use to derive the corresponding coef. If alpha = 0, we
          * cannot move the ray in the recession cone, i.e. strengthening is not possible */
         if( SCIPisZero(scip, alpha) )
         {
            *rho = -SCIPinfinity(scip);
            return SCIP_OKAY;
         }
         else
            currentrho = (alpha - 1) * interpoints[i] / alpha;
      }

      if( currentrho < *rho )
         *rho = currentrho;
   }

   return SCIP_OKAY;
}

/** computes negative steplengths for the rays that are in the recession cone of the S-free set, i.e.,
 * which have an infinite intersection point.
 */
static
SCIP_RETCODE computeNegCutcoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real*            rays,               /**< rays */
   int                   nrays,              /**< number of nonzero rays */
   int*                  rayslppos,          /**< lppos of nonzero rays */
   SCIP_Real*            interpoints,        /**< intersection points */
   SCIP_ROWPREP*         rowprep,            /**< rowprep for the generated cut */
   SCIP_Bool             usebounds,          /**< TRUE if we want to separate non-negative bound */
   SCIP_Real*            ad,                 /**< coefs a and d for the hyperplane aTx + dTy <= 0 */
   SCIP_Bool*            success             /**< if a cut candidate could be computed */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int i;

   *success = TRUE;

   cols = SCIPgetLPCols(scip);
   rows = SCIPgetLPRows(scip);

   /* go through all intersection points that are equal to infinity -> these correspond to the rays which are in the
    * recession cone of C, i.e. the rays for which we (possibly) can compute a negative steplength */
   for( i = 0; i < nrays ; ++i )
   {
      SCIP_Real rho;
      SCIP_Real cutcoef;
      int lppos;

      if( !SCIPisInfinity(scip, interpoints[i]) )
         continue;

      /* compute the smallest rho */
      SCIP_CALL( findRho(scip, rays, nrays, i, interpoints, vars, &rho, usebounds, ad, success) );

      if( ! *success )
         continue;

      /* compute cut coef */
      cutcoef = SCIPisInfinity(scip, -rho) ? 0.0 : 1.0 / rho;

      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      lppos = rayslppos[i];
      if( lppos < 0 )
      {
         lppos = -lppos - 1;

         assert(SCIProwGetBasisStatus(rows[lppos]) == SCIP_BASESTAT_LOWER || SCIProwGetBasisStatus(rows[lppos]) ==
               SCIP_BASESTAT_UPPER);

         SCIP_CALL( addRowToCut(scip, rowprep, SCIProwGetBasisStatus(rows[lppos]) == SCIP_BASESTAT_UPPER ? cutcoef :
                  -cutcoef, rows[lppos], success) ); /* rows have flipped base status! */

         if( ! *success )
            return SCIP_OKAY;
      }
      else
      {
         assert(SCIPcolGetBasisStatus(cols[lppos]) == SCIP_BASESTAT_UPPER || SCIPcolGetBasisStatus(cols[lppos]) ==
               SCIP_BASESTAT_LOWER);
         SCIP_CALL( addColToCut(scip, rowprep, SCIPcolGetBasisStatus(cols[lppos]) == SCIP_BASESTAT_UPPER ? -cutcoef :
                  cutcoef, cols[lppos]) );
      }
   }

   return SCIP_OKAY;
}

/** separates cuts for stored principal minors */
static
SCIP_RETCODE separateDeterminant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR*             xik,                /**< variable X_ik = x_i * x_k */
   SCIP_VAR*             xil,                /**< variable X_il = x_i * x_l */
   SCIP_VAR*             xjk,                /**< variable X_jk = x_j * x_k */
   SCIP_VAR*             xjl,                /**< variable X_jl = x_j * x_l */
   SCIP_Bool*            isxikdiag,          /**< is X_ik diagonal? (i.e. i = k) */
   SCIP_Bool*            isxildiag,          /**< is X_il diagonal? (i.e. i = l) */
   SCIP_Bool*            isxjkdiag,          /**< is X_jk diagonal? (i.e. j = k) */
   SCIP_Bool*            isxjldiag,          /**< is X_jl diagonal? (i.e. j = l) */
   int*                  basicvarpos2tableaurow,/**< map from basic var to its tableau row */
   SCIP_HASHMAP*         tableau,            /**< map from var to its tableau row */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* vars[4] = {xik, xjl, xil, xjk};
   SCIP_Real* tableaurows[4];
   SCIP_Real* interpoints;
   SCIP_Real* rays;
   int nrays;
   int* rayslppos;
   int ncols;
   int nrows;
   SCIP_Bool success;
   SCIP_Real ad[4] = {0.0, 0.0, 0.0, 0.0};
   SCIP_Real solxik;
   SCIP_Real solxil;
   SCIP_Real solxjk;
   SCIP_Real solxjl;

   ncols = SCIPgetNLPCols(scip);
   nrows = SCIPgetNLPRows(scip);

   /* allocate memory for intersection points */
   SCIP_CALL( SCIPallocBufferArray(scip, &interpoints, ncols + nrows) );

   /* allocate memory for rays */
   SCIP_CALL( SCIPallocBufferArray(scip, &rays, 4 * (ncols + nrows)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rayslppos, ncols + nrows) );

   /* cut (in the nonbasic space) is of the form alpha^T x >= 1 */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, TRUE) );
   SCIProwprepAddSide(rowprep, 1.0);

   /* check if we have the tableau row of the variable and if not compute it */
   SCIP_CALL( getTableauRows(scip, vars, basicvarpos2tableaurow, tableau, tableaurows, &success) );

   if( ! success )
      goto CLEANUP;

   /* if we want to enforce bounds, set the right a and d to enforce aTx + dTy <= 0 */
   if( sepadata->usebounds )
   {
      solxik = SCIPvarGetLPSol(xik);
      solxil = SCIPvarGetLPSol(xil);
      solxjk = SCIPvarGetLPSol(xjk);
      solxjl = SCIPvarGetLPSol(xjl);

      if( isxikdiag && SCIPisFeasNegative(scip, solxik) )
      {
         ad[0] = -1.0;
         ad[2] = 1.0;
      }
      else if( isxjldiag && SCIPisFeasNegative(scip, solxjl) )
      {
         ad[0] = -1.0;
         ad[2] = -1.0;
      }
      else if( isxildiag && SCIPisFeasNegative(scip, solxil) )
      {
         ad[1] = 1.0;
         ad[3] = -1.0;
      }
      else if( isxjkdiag && SCIPisFeasNegative(scip, solxjk) )
      {
         ad[1] = -1.0;
         ad[3] = -1.0;
      }
   }

   nrays = 0;
   /* loop over each non-basic var; get the ray; compute cut coefficient */
   SCIP_CALL( addCols(scip, vars, tableaurows, rowprep, rays, &nrays, rayslppos, interpoints, sepadata->usebounds, ad, &success) );

   if( ! success )
      goto CLEANUP;

   /* loop over non-basic slack variables */
   SCIP_CALL( addRows(scip, vars, tableaurows, rowprep, rays, &nrays, rayslppos, interpoints, sepadata->usebounds, ad, &success) );

   if( ! success )
      goto CLEANUP;

   /* do strengthening */
   if( sepadata->usestrengthening )
   {
      SCIP_CALL( computeNegCutcoefs(scip, vars, rays, nrays, rayslppos, interpoints, rowprep, sepadata->usebounds, ad, &success) );

      if( ! success )
         goto CLEANUP;
   }

   /* merge coefficients that belong to same variable */
   SCIPmergeRowprepTerms(scip, rowprep);

   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, sepadata->mincutviol, NULL, &success) );

   /* if cleanup was successfull, create row out of rowprep and add it */
   if( success )
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;

      /* create row */
      SCIP_CALL( SCIPgetRowprepRowSepa(scip, &row, rowprep, sepa) );

      assert(SCIPgetCutEfficacy(scip, NULL, row) > 0.0);

      /* add row */
      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

      if( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

CLEANUP:
   SCIPfreeRowprep(scip, &rowprep);
   SCIPfreeBuffer(scip, &rayslppos);
   SCIPfreeBuffer(scip, &rays);
   SCIPfreeBuffer(scip, &interpoints);

   return SCIP_OKAY;
}


/** separates cuts for stored principal minors */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_HASHMAP* tableau = NULL;
   int* basicvarpos2tableaurow = NULL; /* map between basic var and its tableau row */
   int i;

   assert(sepa != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* check whether there are some minors available */
   if( sepadata->nminors == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* loop over the minors and if they are violated build cut */
   for( i = 0; i < sepadata->nminors && (*result != SCIP_CUTOFF); ++i )
   {
      SCIP_VAR* auxvarxik;
      SCIP_VAR* auxvarxil;
      SCIP_VAR* auxvarxjk;
      SCIP_VAR* auxvarxjl;
      SCIP_Bool isauxvarxikdiag;
      SCIP_Bool isauxvarxildiag;
      SCIP_Bool isauxvarxjkdiag;
      SCIP_Bool isauxvarxjldiag;
      SCIP_Real solxik;
      SCIP_Real solxil;
      SCIP_Real solxjk;
      SCIP_Real solxjl;
      SCIP_Real det;

      /* get variables of the i-th minor */
      SCIP_CALL( getMinorVars(sepadata, i, &auxvarxik, &auxvarxil, &auxvarxjk, &auxvarxjl, &isauxvarxikdiag,
            &isauxvarxildiag, &isauxvarxjkdiag, &isauxvarxjldiag) );

      /* get current solution values */
      solxik = SCIPvarGetLPSol(auxvarxik);
      solxil = SCIPvarGetLPSol(auxvarxil);
      solxjk = SCIPvarGetLPSol(auxvarxjk);
      solxjl = SCIPvarGetLPSol(auxvarxjl);

      det = solxik * solxjl - solxil * solxjk;

      if( SCIPisFeasZero(scip, det) )
         continue;

      if( basicvarpos2tableaurow == NULL )
      {
         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, SCIPgetNLPCols(scip)) );
         SCIP_CALL( SCIPhashmapCreate(&tableau, SCIPblkmem(scip), SCIPgetNVars(scip)) );

         /* construct basicvar to tableau row map */
         SCIP_CALL( constructBasicVars2TableauRowMap(scip, basicvarpos2tableaurow) );
      }
      assert(tableau != NULL);

      if( SCIPisFeasPositive(scip, det) )
      {
         SCIP_CALL( separateDeterminant(scip, sepa, sepadata, auxvarxik, auxvarxil, auxvarxjk, auxvarxjl, &isauxvarxikdiag,
                  &isauxvarxildiag, &isauxvarxjkdiag, &isauxvarxjldiag, basicvarpos2tableaurow, tableau, result) );
      }
      else
      {
         assert(SCIPisFeasNegative(scip, det));
         SCIP_CALL( separateDeterminant(scip, sepa, sepadata, auxvarxil, auxvarxik, auxvarxjl, auxvarxjk, &isauxvarxildiag,
                  &isauxvarxikdiag, &isauxvarxjldiag, &isauxvarxjkdiag, basicvarpos2tableaurow, tableau, result) );
      }
   }

   /* all minors were feasible, so no memory to free */
   if( basicvarpos2tableaurow == NULL )
      return SCIP_OKAY;

   /* free memory */
   for( i = 0; i < SCIPhashmapGetNEntries(tableau); ++i )
   {
      SCIP_HASHMAPENTRY* entry;

      entry = SCIPhashmapGetEntry(tableau, i);

      if( entry != NULL )
      {
         SCIP_Real* tableaurow;

         tableaurow = (SCIP_Real *) SCIPhashmapEntryGetImage(entry);

         SCIPfreeBufferArrayNull(scip, &tableaurow);
      }
   }
   SCIPhashmapFree(&tableau);
   SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyMinor)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaInterminor(scip) );

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->minors == NULL);
   assert(sepadata->nminors == 0);
   assert(sepadata->minorssize == 0);

   /* free separator data */
   SCIPfreeBlockMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->randnumgen == NULL);

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &sepadata->randnumgen, DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;
}


/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   assert(sepadata->randnumgen != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &sepadata->randnumgen);

   return SCIP_OKAY;
}


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
static
SCIP_DECL_SEPAINITSOL(sepaInitsolMinor)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* clear separation data */
   SCIP_CALL( sepadataClear(scip, sepadata) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   int ncalls;
   int currentdepth;

   /* need routine to compute eigenvalues/eigenvectors */
   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   currentdepth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call the separator a given number of times at each node */
   if( (currentdepth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (currentdepth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
   {
      SCIPdebugMsg(scip, "reached round limit for node\n");
      return SCIP_OKAY;
   }

   /* try to detect minors */
   SCIP_CALL( detectMinors(scip, sepadata) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, result) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the minor separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaInterminor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata = NULL;
   SCIP_SEPA* sepa = NULL;

   /* create minor separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   BMSclearMemory(sepadata);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpMinor, NULL,
         sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMinor) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMinor) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitMinor) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitMinor) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolMinor) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolMinor) );

   /* add minor separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/usestrengthening",
         "whether to use strengthened intersection cuts to separate minors",
         &sepadata->usestrengthening, FALSE, DEFAULT_USESTRENGTHENING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/usebounds",
         "whether to also enforce nonegativity bounds of principle minors",
         &sepadata->usebounds, FALSE, DEFAULT_USEBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/mincutviol",
         "minimum required violation of a cut",
         &sepadata->mincutviol, FALSE, DEFAULT_MINCUTVIOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxroundsroot",
         "maximal number of separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
