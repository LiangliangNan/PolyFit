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

/**@file   sepa_minor.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  principal minor separator
 * @author Benjamin Mueller
 *
 * @todo detect non-principal minors and use them to derive split cuts
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_minor.h"
#include "scip/cons_nonlinear.h"
#include "scip/nlpi_ipopt.h"

#define SEPA_NAME              "minor"
#define SEPA_DESC              "separator to ensure that 2x2 principal minors of X - xx' are positive semi-definite"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXMINORSCONST     3000 /**< default constant for the maximum number of minors, i.e., max(const, fac * # quadratic terms) */
#define DEFAULT_MAXMINORSFAC       10.0 /**< default factor for the maximum number of minors, i.e., max(const, fac * # quadratic terms) */
#define DEFAULT_MINCUTVIOL         1e-4 /**< default minimum required violation of a cut */
#define DEFAULT_RANDSEED            157 /**< default random seed */
#define DEFAULT_MAXROUNDS            10 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_IGNOREPACKINGCONSS TRUE /**< default for ignoring circle packing constraints during minor detection */

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_VAR**            minors;             /**< variables of 2x2 minors; each minor is stored like (auxvar_x^2,auxvar_y^2,auxvar_xy) */
   int                   nminors;            /**< total number of minors */
   int                   minorssize;         /**< size of minors array */
   int                   maxminorsconst;     /**< constant for the maximum number of minors, i.e., max(const, fac * # quadratic terms) */
   SCIP_Real             maxminorsfac;       /**< factor for the maximum number of minors, i.e., max(const, fac * # quadratic terms) */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   SCIP_Bool             detectedminors;     /**< has minor detection be called? */
   SCIP_Real             mincutviol;         /**< minimum required violation of a cut */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generation */
   SCIP_Bool             ignorepackingconss; /**< whether to ignore circle packing constraints during minor detection */
};

/*
 * Local methods
 */

/** helper method to store a 2x2 minor in the separation data */
static
SCIP_RETCODE sepadataAddMinor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR*             x,                  /**< x variable */
   SCIP_VAR*             y,                  /**< y variable */
   SCIP_VAR*             auxvarxx,           /**< auxiliary variable for x*x */
   SCIP_VAR*             auxvaryy,           /**< auxiliary variable for y*y */
   SCIP_VAR*             auxvarxy            /**< auxiliary variable for x*y */
   )
{
   assert(sepadata != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(x != y);
   assert(auxvarxx != NULL);
   assert(auxvaryy != NULL);
   assert(auxvarxy != NULL);
   assert(auxvarxx != auxvaryy);
   assert(auxvarxx != auxvarxy);
   assert(auxvaryy != auxvarxy);

   SCIPdebugMsg(scip, "store 2x2 minor: %s %s %s for x=%s y=%s\n", SCIPvarGetName(auxvarxx), SCIPvarGetName(auxvaryy),
      SCIPvarGetName(auxvarxy), SCIPvarGetName(x), SCIPvarGetName(y));

   /* reallocate if necessary */
   if( sepadata->minorssize < 5 * (sepadata->nminors + 1) )
   {
      int newsize = SCIPcalcMemGrowSize(scip, 5 * (sepadata->nminors + 1));
      assert(newsize > 5 * (sepadata->nminors + 1));

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->minors), sepadata->minorssize, newsize) );
      sepadata->minorssize = newsize;
   }

   /* store minor */
   sepadata->minors[5 * sepadata->nminors] = x;
   sepadata->minors[5 * sepadata->nminors + 1] = y;
   sepadata->minors[5 * sepadata->nminors + 2] = auxvarxx;
   sepadata->minors[5 * sepadata->nminors + 3] = auxvaryy;
   sepadata->minors[5 * sepadata->nminors + 4] = auxvarxy;
   ++(sepadata->nminors);

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, x) );
   SCIP_CALL( SCIPcaptureVar(scip, y) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvarxx) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvaryy) );
   SCIP_CALL( SCIPcaptureVar(scip, auxvarxy) );

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
   for( i = 0; i < 5 * sepadata->nminors; ++i )
   {
      assert(sepadata->minors[i] != NULL);
      SCIP_CALL( SCIPreleaseVar(scip, &sepadata->minors[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &sepadata->minors, sepadata->minorssize);

   /* reset counters */
   sepadata->nminors = 0;
   sepadata->minorssize = 0;

   return SCIP_OKAY;
}

/** helper method to identify non-overlapping constraints in circle packing */
static
SCIP_Bool isPackingCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< nonlinear constraint */
   )
{
   SCIP_EXPR* root;
   SCIP_VAR* quadvars[4] = {NULL, NULL, NULL, NULL};
   SCIP_VAR* bilinvars[4] = {NULL, NULL, NULL, NULL};
   int nbilinvars = 0;
   int nquadvars = 0;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   root = SCIPgetExprNonlinear(cons);
   assert(root != NULL);
   nchildren = SCIPexprGetNChildren(root);

   /* non-overlapping constraint has 6 terms (2 bilinear + 4 quadratic) */
   if( nchildren != 6 || !SCIPisExprSum(scip, root) )
      return FALSE;

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_EXPR* expr;
      SCIP_EXPR** children;

      /* get child */
      expr = SCIPexprGetChildren(root)[i];
      assert(expr != NULL);
      children = SCIPexprGetChildren(expr);

      /* case: expr = x^2; x is no auxiliary variable */
      if( SCIPisExprPower(scip, expr) && SCIPgetExponentExprPow(expr) == 2.0
         && SCIPisExprVar(scip, children[0]) )
      {
         SCIP_VAR* x;

         /* too many quadratic variables -> stop */
         if( nquadvars > 3 )
            return FALSE;

         x = SCIPgetVarExprVar(children[0]);
         assert(x != NULL);

         quadvars[nquadvars++] = x;
      }
      /* case: expr = x * y; x and y are no auxiliary variables */
      else if( SCIPisExprProduct(scip, expr) && SCIPexprGetNChildren(expr) == 2
         && SCIPisExprVar(scip, children[0]) && SCIPisExprVar(scip, children[1]) )
      {
         SCIP_VAR* x;
         SCIP_VAR* y;

         /* too many bilinear variables -> stop */
         if( nbilinvars > 2 )
            return FALSE;

         x = SCIPgetVarExprVar(children[0]);
         assert(x != NULL);
         y = SCIPgetVarExprVar(children[1]);
         assert(y != NULL);
         assert(x != y);

         bilinvars[nbilinvars++] = x;
         bilinvars[nbilinvars++] = y;
      }
      else
      {
         return FALSE;
      }
   }

   /* number of bilinear and quadratic terms do not fit */
   if( nbilinvars != 4 || nquadvars != 4 )
      return FALSE;

   /* each quadratic variable has to appear in exactly one bilinear terms */
   for( i = 0; i < nquadvars; ++i )
   {
      int counter = 0;
      int j;

      for( j = 0; j < nbilinvars; ++j )
      {
         if( quadvars[i] == bilinvars[j] )
            ++counter;
      }

      if( counter != 1 )
         return FALSE;
   }

   return TRUE;
}

/** helper method to get the variables associated to a minor */
static
SCIP_RETCODE getMinorVars(
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   int                   idx,                /**< index of the stored minor */
   SCIP_VAR**            x,                  /**< pointer to store x variable */
   SCIP_VAR**            y,                  /**< pointer to store x variable */
   SCIP_VAR**            auxvarxx,           /**< pointer to store auxiliary variable for x*x */
   SCIP_VAR**            auxvaryy,           /**< pointer to store auxiliary variable for y*y */
   SCIP_VAR**            auxvarxy            /**< pointer to store auxiliary variable for x*y */
   )
{
   assert(sepadata != NULL);
   assert(idx >= 0 && idx < sepadata->nminors);
   assert(auxvarxx != NULL);
   assert(auxvaryy != NULL);
   assert(auxvarxy != NULL);

   *x = sepadata->minors[5 * idx];
   *y = sepadata->minors[5 * idx + 1];
   *auxvarxx = sepadata->minors[5 * idx + 2];
   *auxvaryy = sepadata->minors[5 * idx + 3];
   *auxvarxy = sepadata->minors[5 * idx + 4];

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
   SCIP_HASHMAP* quadmap;
   SCIP_VAR** xs;
   SCIP_VAR** ys;
   SCIP_VAR** auxvars;
   int* perm = NULL;
   int nbilinterms = 0;
   int nquadterms = 0;
   int maxminors;
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
   SCIP_CALL( SCIPhashmapCreate(&quadmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &xs, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ys, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxvars, SCIPgetNVars(scip)) );

   /* initialize iterator */
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR);

   for( c = 0; c < SCIPconshdlrGetNConss(conshdlr); ++c )
   {
      SCIP_CONS* cons;
      SCIP_EXPR* expr;
      SCIP_EXPR* root;

      cons = SCIPconshdlrGetConss(conshdlr)[c];
      assert(cons != NULL);
      root = SCIPgetExprNonlinear(cons);
      assert(root != NULL);

      /* ignore circle packing constraints; the motivation for this is that in circle packing instance not only the SDP
       * relaxation is weak (see "Packing circles in a square: a theoretical comparison of various convexification
       * techniques", http://www.optimization-online.org/DB_HTML/2017/03/5911.html), but it also hurts performance
       */
      if( sepadata->ignorepackingconss && isPackingCons(scip, cons) )
      {
         SCIPdebugMsg(scip, "ignore packing constraints %s\n", SCIPconsGetName(cons));
         continue;
      }

      for( expr = SCIPexpriterRestartDFS(it, root); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/ /*lint !e440*/
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
            assert(!SCIPhashmapExists(quadmap, (void*)quadvar));
            SCIPdebugMsg(scip, "found %s = (%s)^2\n", SCIPvarGetName(auxvar), SCIPvarGetName(quadvar));

            /* hash the quadratic variable to its corresponding auxiliary variable */
            SCIP_CALL( SCIPhashmapInsert(quadmap, (void*)quadvar, auxvar) );
            ++nquadterms;
         }
         /* check for expr = x * y */
         else if( SCIPexprGetNChildren(expr) == 2 && SCIPisExprProduct(scip, expr)
            && SCIPgetExprAuxVarNonlinear(children[0]) != NULL && SCIPgetExprAuxVarNonlinear(children[1]) != NULL )
         {
            SCIP_VAR* x;
            SCIP_VAR* y;

            assert(children[0] != NULL);
            assert(children[1] != NULL);

            x = SCIPgetExprAuxVarNonlinear(children[0]);
            y = SCIPgetExprAuxVarNonlinear(children[1]);

            /* ignore binary variables */
            if( !SCIPvarIsBinary(x) && !SCIPvarIsBinary(y) )
            {
               xs[nbilinterms] = SCIPgetExprAuxVarNonlinear(children[0]);
               ys[nbilinterms] = SCIPgetExprAuxVarNonlinear(children[1]);
               auxvars[nbilinterms] = auxvar;
               SCIPdebugMsg(scip, "found %s = %s * %s\n", SCIPvarGetName(auxvar), SCIPvarGetName(xs[nbilinterms]), SCIPvarGetName(ys[nbilinterms]));
               ++nbilinterms;
            }
         }
      }
   }
   assert(nbilinterms < SCIPgetNVars(scip));
   SCIPdebugMsg(scip, "stored %d bilinear terms in total\n", nbilinterms);

   /* use max(maxminorsconst, maxminorsfac * # quadratic terms) as a limit for the maximum number of minors */
   maxminors = (int) MAX(sepadata->maxminorsconst, sepadata->maxminorsfac * nquadterms);
   SCIPdebugMsg(scip, "maximum number of minors = %d\n", maxminors);

   /* permute bilinear terms if there are too many of them; the motivation for this is that we don't want to
    * prioritize variables because of the order in the bilinear terms where they appear; however, variables that
    * appear more often in bilinear terms might be more important than others so the corresponding bilinear terms
    * are more likely to be chosen
    */
   if( maxminors < nbilinterms && maxminors < SQR(nquadterms) )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, nbilinterms) );

      for( i = 0; i < nbilinterms; ++i )
         perm[i] = i;

      /* permute array */
      SCIPrandomPermuteIntArray(sepadata->randnumgen, perm, 0, nbilinterms);
   }

   /* store 2x2 principal minors */
   for( i = 0; i < nbilinterms && sepadata->nminors < maxminors; ++i )
   {
      SCIP_VAR* x;
      SCIP_VAR* y;
      SCIP_VAR* auxvarxy;

      if( perm == NULL )
      {
         x = xs[i];
         y = ys[i];
         auxvarxy = auxvars[i];
      }
      else
      {
         x = xs[perm[i]];
         y = ys[perm[i]];
         auxvarxy = auxvars[perm[i]];
      }

      assert(x != NULL);
      assert(y != NULL);
      assert(auxvarxy != NULL);
      assert(x != y);

      if( SCIPhashmapExists(quadmap, (void*)x) && SCIPhashmapExists(quadmap, (void*)y) )
      {
         SCIP_VAR* auxvarxx;
         SCIP_VAR* auxvaryy;

         auxvarxx = (SCIP_VAR*)SCIPhashmapGetImage(quadmap, (void*)x);
         assert(auxvarxx != NULL);
         auxvaryy = (SCIP_VAR*)SCIPhashmapGetImage(quadmap, (void*)y);
         assert(auxvaryy != NULL);

         /* store minor into the separation data */
         SCIP_CALL( sepadataAddMinor(scip, sepadata, x, y, auxvarxx, auxvaryy, auxvarxy) );
      }
   }
   SCIPdebugMsg(scip, "found %d principal minors in total\n", sepadata->nminors);

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &perm);
   SCIPfreeBufferArray(scip, &auxvars);
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);
   SCIPhashmapFree(&quadmap);
   SCIPfreeExpriter(&it);

#ifdef SCIP_STATISTIC
   totaltime += SCIPgetTotalTime(scip);
   SCIPstatisticMessage("MINOR DETECT %s %f %d %d\n", SCIPgetProbName(scip), totaltime, sepadata->nminors, maxminors);
#endif

   return SCIP_OKAY;
}

/** helper method to compute eigenvectors and eigenvalues */
static
SCIP_RETCODE getEigenValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x,                  /**< solution value of x */
   SCIP_Real             y,                  /**< solution value of y */
   SCIP_Real             xx,                 /**< solution value of x*x */
   SCIP_Real             yy,                 /**< solution value of y*y */
   SCIP_Real             xy,                 /**< solution value of x*y */
   SCIP_Real*            eigenvals,          /**< array to store eigenvalues (at least of size 3) */
   SCIP_Real*            eigenvecs,          /**< array to store eigenvalues (at least of size 9) */
   SCIP_Bool*            success             /**< pointer to store whether eigenvalue computation was successful */
   )
{
   assert(eigenvals != NULL);
   assert(eigenvecs != NULL);
   assert(success != NULL);

   *success = TRUE;

   /* construct matrix */
   eigenvecs[0] = 1.0;
   eigenvecs[1] = x;
   eigenvecs[2] = y;
   eigenvecs[3] = x;
   eigenvecs[4] = xx;
   eigenvecs[5] = xy;
   eigenvecs[6] = y;
   eigenvecs[7] = xy;
   eigenvecs[8] = yy;

   /* use LAPACK to compute the eigenvalues and eigenvectors */
   if( SCIPcallLapackDsyevIpopt(TRUE, 3, eigenvecs, eigenvals) != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Failed to compute eigenvalues and eigenvectors of augmented quadratic form matrix.\n");
      *success = FALSE;
   }

   return SCIP_OKAY;
}

/** generate and add a cut */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< solution to separate (might be NULL) */
   SCIP_VAR*             x,                  /**< x variable */
   SCIP_VAR*             y,                  /**< y variable */
   SCIP_VAR*             xx,                 /**< auxiliary variable for x*x */
   SCIP_VAR*             yy,                 /**< auxiliary variable for y*y */
   SCIP_VAR*             xy,                 /**< auxiliary variable for x*y */
   SCIP_Real*            eigenvec,           /**< array containing an eigenvector */
   SCIP_Real             eigenval,           /**< eigenvalue */
   SCIP_Real             mincutviol,         /**< minimal required violation */
   SCIP_RESULT*          result              /**< pointer to update the result */
   )
{
   SCIP_VAR* vars[5] = {x, y, xx, yy, xy};
   SCIP_Real coefs[5];
   SCIP_Real constant;
   SCIP_ROWPREP* rowprep;
   SCIP_Bool success;

   assert(x != NULL);
   assert(y != NULL);
   assert(xx != NULL);
   assert(yy != NULL);
   assert(xy != NULL);
   assert(eigenvec != NULL);
   assert(mincutviol >= 0.0);
   assert(result != NULL);

   /* check whether the resulting cut is violated enough */
   if( !SCIPisFeasLT(scip, eigenval, -mincutviol) )
      return SCIP_OKAY;

   /* the resulting cut reads as
    *              (1 x  y )  (v0)
    *  (v0 v1 v2)  (x xx xy)  (v1)  >= 0
    *              (y xy yy)  (v2)
    *  where v is the eigenvector corresponding to a negative eigenvalue
    *  that is,
    *  v0^2 + 2 v0 v1 * x + 2 v0 v2 * y + v1^2 * xx + v2^2 * yy + 2 v1 v2 * xy >= 0
    */
   constant = SQR(eigenvec[0]);
   coefs[0] = 2.0 * eigenvec[0] * eigenvec[1];
   coefs[1] = 2.0 * eigenvec[0] * eigenvec[2];
   coefs[2] = SQR(eigenvec[1]);
   coefs[3] = SQR(eigenvec[2]);
   coefs[4] = 2.0 * eigenvec[1] * eigenvec[2];

   /* create rowprep */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, FALSE) );
   SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, 5, vars, coefs) );
   SCIProwprepAddConstant(rowprep, constant);
   SCIPdebug( SCIPprintRowprep(scip, rowprep, NULL) );
   SCIPdebugMsg(scip, "cut violation %g mincutviol = %g\n", SCIPgetRowprepViolation(scip, rowprep, sol, NULL), mincutviol);

   /* cleanup coefficient and side, esp treat epsilon to integral values; don't consider scaling up here */
   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, 0.0, NULL, &success) );

   /* check cut violation */
   if( success && SCIPgetRowprepViolation(scip, rowprep, sol, NULL) > mincutviol )
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;

      /* set name of rowprep */
      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "minor_%s_%s_%s_%lld", SCIPvarGetName(xx), SCIPvarGetName(yy),
         SCIPvarGetName(xy), SCIPgetNLPs(scip));

      /* create, add, and release row */
      SCIP_CALL( SCIPgetRowprepRowSepa(scip, &row, rowprep, sepa) );
      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );

      /* update result pointer */
      *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
   }

   /* free rowprep */
   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** separates cuts for stored principal minors */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_SEPADATA* sepadata;
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

   for( i = 0; i < sepadata->nminors && (*result != SCIP_CUTOFF); ++i )
   {
      SCIP_Real eigenvals[3];
      SCIP_Real eigenvecs[9];
      SCIP_VAR* x;
      SCIP_VAR* y;
      SCIP_VAR* xx;
      SCIP_VAR* yy;
      SCIP_VAR* xy;
      SCIP_Real solx;
      SCIP_Real soly;
      SCIP_Real solxx;
      SCIP_Real solyy;
      SCIP_Real solxy;
      SCIP_Bool success;
      int k;

      /* get variables of the i-th minor */
      SCIP_CALL( getMinorVars(sepadata, i, &x, &y, &xx, &yy, &xy) );
      assert(x != NULL);
      assert(y != NULL);
      assert(xx != NULL);
      assert(yy != NULL);
      assert(xy != NULL);

      /* get current solution values */
      solx = SCIPgetSolVal(scip, sol, x);
      soly = SCIPgetSolVal(scip, sol, y);
      solxx = SCIPgetSolVal(scip, sol, xx);
      solyy = SCIPgetSolVal(scip, sol, yy);
      solxy = SCIPgetSolVal(scip, sol, xy);
      SCIPdebugMsg(scip, "solution values (x,y,xx,yy,xy)=(%g,%g,%g,%g,%g)\n", solx, soly, solxx, solyy, solxy);

      /* compute eigenvalues and eigenvectors */
      SCIP_CALL( getEigenValues(scip, solx, soly, solxx, solyy, solxy, eigenvals, eigenvecs, &success) );
      if( !success )
         continue;

      /* try to generate a cut for each negative eigenvalue */
      for( k = 0; k < 3 && (*result != SCIP_CUTOFF); ++k )
      {
         SCIPdebugMsg(scip, "eigenvalue = %g  eigenvector = (%g,%g,%g)\n", eigenvals[k], eigenvecs[3*k], eigenvecs[3*k + 1], eigenvecs[3*k + 2]);
         SCIP_CALL( addCut(scip, sepa, sol, x, y, xx, yy, xy, &eigenvecs[3*k], eigenvals[k], sepadata->mincutviol, result) );
         SCIPdebugMsg(scip, "result: %d\n", *result);
      }
   }

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
   SCIP_CALL( SCIPincludeSepaMinor(scip) );

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

   /* need routine to compute eigenvalues/eigenvectors */
   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
   {
      SCIPdebugMsg(scip, "reached round limit for node\n");
      return SCIP_OKAY;
   }

   /* try to detect minors */
   SCIP_CALL( detectMinors(scip, sepadata) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, NULL, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMinor)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   int ncalls;

   /* need routine to compute eigenvalues/eigenvectors */
   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
   {
      SCIPdebugMsg(scip, "reached round limit for node\n");
      return SCIP_OKAY;
   }

   /* try to detect minors */
   SCIP_CALL( detectMinors(scip, SCIPsepaGetData(sepa)) );

   /* call separation method */
   SCIP_CALL( separatePoint(scip, sepa, sol, result) );

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the minor separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMinor(
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
         sepaExeclpMinor, sepaExecsolMinor,
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
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxminorsconst",
         "constant for the maximum number of minors, i.e., max(const, fac * # quadratic terms)",
         &sepadata->maxminorsconst, FALSE, DEFAULT_MAXMINORSCONST, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/maxminorsfac",
         "factor for the maximum number of minors, i.e., max(const, fac * # quadratic terms)",
         &sepadata->maxminorsfac, FALSE, DEFAULT_MAXMINORSFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );

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

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/ignorepackingconss",
         "whether to ignore circle packing constraints during minor detection",
         &sepadata->ignorepackingconss, FALSE, DEFAULT_IGNOREPACKINGCONSS, NULL, NULL) );

   return SCIP_OKAY;
}
