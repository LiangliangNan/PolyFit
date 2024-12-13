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

/**@file   nlhdlr_soc.c
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  nonlinear handler for second order cone constraints

 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Fabian Wegscheider
 *
 * This is a nonlinear handler for second order cone constraints of the form
 *
 * \f[\sqrt{\sum_{i=1}^{n} (v_i^T x + \beta_i)^2} \leq v_{n+1}^T x + \beta_{n+1}.\f]
 *
 * Note that \f$v_i\f$, for \f$i \leq n\f$, could be 0, thus allowing a positive constant term inside the root.
 *
 * @todo test if it makes sense to only disaggregate when nterms > some parameter
 *
 */

#include <string.h>

#include "scip/nlhdlr_soc.h"
#include "scip/cons_nonlinear.h"
#include "scip/expr_pow.h"
#include "scip/expr_sum.h"
#include "scip/expr_var.h"
#include "scip/debug.h"
#include "scip/pub_nlhdlr.h"
#include "scip/nlpi_ipopt.h"


/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "soc"
#define NLHDLR_DESC               "nonlinear handler for second-order cone structures"
#define NLHDLR_DETECTPRIORITY       100 /**< priority of the nonlinear handler for detection */
#define NLHDLR_ENFOPRIORITY         100 /**< priority of the nonlinear handler for enforcement */
#define DEFAULT_MINCUTEFFICACY     1e-5 /**< default value for parameter mincutefficacy */
#define DEFAULT_COMPEIGENVALUES    TRUE /**< default value for parameter compeigenvalues */

/*
 * Data structures
 */

/** nonlinear handler expression data. The data is structured in the following way:
 *
 *  A 'term' is one of the arguments of the quadratic terms, i.e. \f$v_i^T x + beta_i\f$.
 *  The last term is always the one on the right-hand side. This means that nterms is
 *  equal to n+1 in the above description.
 *
 *  - vars contains a list of all expressions which are treated as variables (no duplicates)
 *  - offsets contains the constants beta_i of each term
 *  - transcoefs contains the non-zero values of the transformation vectors v_i of each term
 *  - transcoefsidx contains for each entry of transcoefs the position of the respective variable in vars
 *  - termbegins contains the index at which the transcoefs of each term start, with a sentinel value
 *  - nterms is the total number of terms appearing on both sides
 *  - nvars is the total number of unique variables appearing (length of vars)
 *
 *  Note that the numbers of nonzeroes in v_i is termbegins[i+1] - termbegins[i] and that
 *  the total number of entries in transcoefs and transcoefsidx is termbegins[nterms]
 *
 *  The disaggregation is implicitly stored in the variables disvars and disrow. An SOC as
 *  described above is replaced by n smaller SOCs
 *
 *              (v_i^T x + beta_i)^2 <= disvar_i     * (v_{n+1}^T x + beta_{n+1})
 *
 *  and the row       sum_i disvar_i <= v_{n+1}^T x + beta_{n+1}.
 *
 *  The disaggregation only happens if we have more than 3 terms.
 *
 *  Example: The constraint SQRT(5 + (3x - 4y + 2)^2 + y^2 + 7z^2) <= 5x - y - 1
 *           results in the following nlhdlrexprdata:
 *
 *           vars = {x, y, z}
 *           offsets = {2, 0, 0, SQRT(5), -1}
 *           transcoefs = {3, -4, 1, SQRT(7), 5, -1}
 *           transcoefsidx = {0, 1, 1, 2, 0, 1}
 *           termbegins = {0, 2, 3, 4, 4, 6}
 *           nvars = 3
 *           nterms = 5
 *
 * @note: due to the current implementation, the constant term is the second to last term, except when the SOC was a rotated
 * SOC, e.g., 1 + x^2 - y*z, i.e., when detected by detectSocQuadraticSimple. In that case, the constant is third to
 * last term.
 */
struct SCIP_NlhdlrExprData
{
   SCIP_EXPR**           vars;               /**< expressions which (aux)variables appear on both sides (x) */
   SCIP_Real*            offsets;            /**< offsets of both sides (beta_i) */
   SCIP_Real*            transcoefs;         /**< non-zeros of linear transformation vectors (v_i) */
   int*                  transcoefsidx;      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  termbegins;         /**< starting indices of transcoefs for each term */
   int                   nvars;              /**< total number of variables appearing */
   int                   nterms;             /**< number of summands in the SQRT +1 for RHS (n+1) */

   /* variables for cone disaggregation */
   SCIP_VAR**            disvars;            /**< disaggregation variables for each term in lhs */
   SCIP_ROW*             disrow;             /**< disaggregation row */
};

struct SCIP_NlhdlrData
{
   SCIP_Real             mincutefficacy;     /**< minimum efficacy a cut need to be added */
   SCIP_Bool             compeigenvalues;    /**< whether Eigenvalue computations should be done to detect complex cases */
};

/*
 * Local methods
 */

#ifdef SCIP_DEBUG
/** prints the nlhdlr expression data */
static
void printNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< pointer to store nonlinear handler expression data */
   )
{
   int nterms;
   int i;
   int j;

   nterms = nlhdlrexprdata->nterms;

   SCIPinfoMessage(scip, NULL, "SQRT( ");

   for( i = 0; i < nterms - 1; ++i )
   {
      int startidx;

      startidx = nlhdlrexprdata->termbegins[i];

      /* v_i is 0 */
      if( startidx == nlhdlrexprdata->termbegins[i + 1] )
      {
         assert(nlhdlrexprdata->offsets[i] != 0.0);

         SCIPinfoMessage(scip, NULL, "%f", SQR(nlhdlrexprdata->offsets[i]));
         continue;
      }

      /* v_i is not 0 */
      SCIPinfoMessage(scip, NULL, "(");

      for( j = startidx; j < nlhdlrexprdata->termbegins[i + 1]; ++j )
      {
         if( nlhdlrexprdata->transcoefs[j] != 1.0 )
            SCIPinfoMessage(scip, NULL, "%f*", nlhdlrexprdata->transcoefs[j]);
         if( SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]) != NULL )
         {
            SCIPinfoMessage(scip, NULL, "%s", SCIPvarGetName(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]])));
            SCIPinfoMessage(scip, NULL, "(%p)", (void*)nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]);
         }
         else
            SCIPinfoMessage(scip, NULL, "%p", (void*)nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]);

         if( j < nlhdlrexprdata->termbegins[i + 1] - 1 )
            SCIPinfoMessage(scip, NULL, " + ");
         else if( nlhdlrexprdata->offsets[i] != 0.0 )
            SCIPinfoMessage(scip, NULL, " + %f", nlhdlrexprdata->offsets[i]);
      }

      SCIPinfoMessage(scip, NULL, ")^2");

      if( i < nterms - 2 )
         SCIPinfoMessage(scip, NULL, " + ");
   }

   SCIPinfoMessage(scip, NULL, " ) <= ");

   for( j = nlhdlrexprdata->termbegins[nterms-1]; j < nlhdlrexprdata->termbegins[nterms]; ++j )
   {
      if( nlhdlrexprdata->transcoefs[j] != 1.0 )
         SCIPinfoMessage(scip, NULL, "%f*", nlhdlrexprdata->transcoefs[j]);
      if( SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]) != NULL )
         SCIPinfoMessage(scip, NULL, "%s", SCIPvarGetName(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]])));
      else
         SCIPinfoMessage(scip, NULL, "%p", (void*)nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]);

      if( j < nlhdlrexprdata->termbegins[nterms] - 1 )
         SCIPinfoMessage(scip, NULL, " + ");
      else if( nlhdlrexprdata->offsets[nterms-1] != 0.0 )
         SCIPinfoMessage(scip, NULL, " + %f", nlhdlrexprdata->offsets[nterms-1]);
   }

   SCIPinfoMessage(scip, NULL, "\n");
}
#endif

/** helper method to create variables for the cone disaggregation */
static
SCIP_RETCODE createDisaggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< nonlinear handler expression data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int ndisvars;
   int i;

   assert(nlhdlrexprdata != NULL);

   ndisvars = nlhdlrexprdata->nterms - 1;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexprdata->disvars, ndisvars) );

   /* create disaggregation variables representing the epigraph of (v_i^T x + beta_i)^2 / (v_{n+1}^T x + beta_{n+1}) */
   for( i = 0; i < ndisvars; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_%d", (void*) expr, i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &nlhdlrexprdata->disvars[i], name, 0.0, SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS) );
      SCIPvarMarkRelaxationOnly(nlhdlrexprdata->disvars[i]);

      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[i]) );
      SCIP_CALL( SCIPaddVarLocksType(scip, nlhdlrexprdata->disvars[i], SCIP_LOCKTYPE_MODEL, 1, 1) );
   }

   return SCIP_OKAY;
}

/** helper method to free variables for the cone disaggregation */
static
SCIP_RETCODE freeDisaggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< nonlinear handler expression data */
   )
{
   int ndisvars;
   int i;

   assert(nlhdlrexprdata != NULL);

   if( nlhdlrexprdata->disvars == NULL )
      return SCIP_OKAY;

   ndisvars = nlhdlrexprdata->nterms - 1;

   /* release variables */
   for( i = 0; i < ndisvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocksType(scip, nlhdlrexprdata->disvars[i], SCIP_LOCKTYPE_MODEL, -1, -1) );
      SCIP_CALL( SCIPreleaseVar(scip, &nlhdlrexprdata->disvars[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &nlhdlrexprdata->disvars, ndisvars);

   return SCIP_OKAY;
}

/** helper method to create the disaggregation row \f$\text{disvars}_i \leq v_{n+1}^T x + \beta_{n+1}\f$ */
static
SCIP_RETCODE createDisaggrRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< nonlinear handler expression data */
   )
{
   SCIP_Real beta;
   char name[SCIP_MAXSTRLEN];
   int ndisvars;
   int nterms;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->disrow == NULL);

   nterms = nlhdlrexprdata->nterms;
   beta = nlhdlrexprdata->offsets[nterms - 1];

   ndisvars = nterms - 1;

   /* create row 0 <= beta_{n+1} */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_row", (void*) expr);
   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &nlhdlrexprdata->disrow, conshdlr, name,
         -SCIPinfinity(scip), beta, FALSE, FALSE, TRUE) );

   /* add disvars to row */
   for( i = 0; i < ndisvars; ++i )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, nlhdlrexprdata->disrow, nlhdlrexprdata->disvars[i], 1.0) );
   }

   /* add rhs vars to row */
   for( i = nlhdlrexprdata->termbegins[nterms - 1]; i < nlhdlrexprdata->termbegins[nterms]; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[i]]);
      assert(var != NULL);

      coef = -nlhdlrexprdata->transcoefs[i];

      SCIP_CALL( SCIPaddVarToRow(scip, nlhdlrexprdata->disrow, var, coef) );
   }

   return SCIP_OKAY;
}

/** helper method to create nonlinear handler expression data */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           vars,               /**< expressions which variables appear on both sides (\f$x\f$) */
   SCIP_Real*            offsets,            /**< offsets of bot sides (\f$beta_i\f$) */
   SCIP_Real*            transcoefs,         /**< non-zeroes of linear transformation vectors (\f$v_i\f$) */
   int*                  transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  termbegins,         /**< starting indices of transcoefs for each term */
   int                   nvars,              /**< total number of variables appearing */
   int                   nterms,             /**< number of summands in the SQRT, +1 for RHS */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata      /**< pointer to store nonlinear handler expression data */
   )
{
   int ntranscoefs;

   assert(vars != NULL);
   assert(offsets != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(termbegins != NULL);
   assert(nlhdlrexprdata != NULL);

   ntranscoefs = termbegins[nterms];

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexprdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, vars, nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->offsets, offsets, nterms) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefs, transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefsidx, transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->termbegins, termbegins, nterms + 1) );
   (*nlhdlrexprdata)->nvars = nvars;
   (*nlhdlrexprdata)->nterms = nterms;

   (*nlhdlrexprdata)->disrow = NULL;
   (*nlhdlrexprdata)->disvars = NULL;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "created nlhdlr data for the following soc expression:\n");
   printNlhdlrExprData(scip, *nlhdlrexprdata);
   printf("x is %p\n", (void *)vars[0]);
#endif

   return SCIP_OKAY;
}

/** helper method to free nonlinear handler expression data */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata      /**< pointer to free nonlinear handler expression data */
   )
{
   int ntranscoefs;

   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   /* free variables and row for cone disaggregation */
   SCIP_CALL( freeDisaggrVars(scip, *nlhdlrexprdata) );

   ntranscoefs = (*nlhdlrexprdata)->termbegins[(*nlhdlrexprdata)->nterms];

   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->termbegins, (*nlhdlrexprdata)->nterms + 1);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefsidx, ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefs, ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->offsets, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, (*nlhdlrexprdata)->nvars);
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** evaluate a single term of the form \f$v_i^T x + \beta_i\f$ */
static
SCIP_Real evalSingleTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata,      /**< nonlinear handler expression data */
   SCIP_SOL*             sol,                /**< solution */
   int                   k                   /**< term to be evaluated */
   )
{
   SCIP_Real result;
   int i;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(0 <= k && k < nlhdlrexprdata->nterms);

   result = nlhdlrexprdata->offsets[k];

   for( i = nlhdlrexprdata->termbegins[k]; i < nlhdlrexprdata->termbegins[k + 1]; ++i )
   {
      SCIP_Real varval = SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[i]]));
      result += nlhdlrexprdata->transcoefs[i] * varval;
   }

   return result;
}

/** computes gradient cut for a 2D or 3D SOC
 *
 *  A 3D SOC looks like
 *  \f[
 *    \sqrt{ (v_1^T x + \beta_1)^2 + (v_2^T x + \beta_2)^2 } \leq v_3^T x + \beta_3
 *  \f]
 *
 *  Let \f$f(x)\f$ be the left-hand-side. The partial derivatives of \f$f\f$ are given by
 *  \f[
 *    \frac{\delta f}{\delta x_j} = \frac{(v_1)_j(v_1^T x + \beta_1) + (v_2)_j (v_2^T x + \beta_2)}{f(x)}
 *  \f]
 *
 *  and the gradient cut is then \f$f(x^*) + \nabla f(x^*)(x - x^*) \leq v_3^T x + \beta_3\f$.
 *
 *  A 2D SOC is
 *  \f[
 *    |v_1^T x + \beta_1| \leq v_2^T x + \beta_2
 *  \f]
 *  but we build the cut using the same procedure as for 3D.
 */
static
SCIP_RETCODE generateCutSolSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_CONS*            cons,               /**< the constraint that expr is part of */
   SCIP_SOL*             sol,                /**< solution to separate or NULL for the LP solution */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_Real             mincutviolation,    /**< minimal required cut violation */
   SCIP_Real             rhsval,             /**< value of last term at sol */
   SCIP_ROW**            cut                 /**< pointer to store a cut */
   )
{
   SCIP_ROWPREP* rowprep;
   SCIP_Real* transcoefs;
   SCIP_Real cutcoef;
   SCIP_Real fvalue;
   SCIP_Real valterms[2] = {0.0, 0.0}; /* for lint */
   SCIP_Real cutrhs;
   SCIP_EXPR** vars;
   SCIP_VAR* cutvar;
   int* transcoefsidx;
   int* termbegins;
   int nterms;
   int i;
   int j;

   assert(expr != NULL);
   assert(cons != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(mincutviolation >= 0.0);
   assert(cut != NULL);

   vars = nlhdlrexprdata->vars;
   transcoefs = nlhdlrexprdata->transcoefs;
   transcoefsidx = nlhdlrexprdata->transcoefsidx;
   termbegins = nlhdlrexprdata->termbegins;
   nterms = nlhdlrexprdata->nterms;

   *cut = NULL;

   /* evaluate lhs terms and compute f(x*) */
   fvalue = 0.0;
   for( i = 0; i < nterms - 1; ++i )
   {
      valterms[i] = evalSingleTerm(scip, nlhdlrexprdata, sol, i);
      fvalue += SQR( valterms[i] );
   }
   fvalue = SQRT( fvalue );

   /* don't generate cut if we are not violated @todo: remove this once core detects better when a nlhdlr's cons is
    * violated
    */
   if( fvalue - rhsval <= mincutviolation )
   {
      SCIPdebugMsg(scip, "do not generate cut: rhsval %g, fvalue %g violation is %g\n", rhsval, fvalue, fvalue - rhsval);
      return SCIP_OKAY;
   }

   /* if f(x*) = 0 then SOC can't be violated and we shouldn't be here */
   assert(fvalue > 0.0);

   /* create cut */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, termbegins[nterms]) );

   /* cut is f(x*) + \nabla f(x*)^T (x - x*) \leq v_n^T x + \beta_n, i.e.,
    * \nabla f(x*)^T x - v_n^T x \leq \beta_n + \nabla f(x*)^T x* - f(x*)
    * thus cutrhs is \beta_n - f(x*) + \nabla f(x*)^T x*
    */
   cutrhs = nlhdlrexprdata->offsets[nterms - 1] - fvalue;

   /* add cut coefficients from lhs terms and compute cut's rhs */
   for( j = 0; j < nterms - 1; ++j )
   {
      for( i = termbegins[j]; i < termbegins[j + 1]; ++i )
      {
         cutvar = SCIPgetExprAuxVarNonlinear(vars[transcoefsidx[i]]);

         /* cutcoef is (the first part of) the partial derivative w.r.t cutvar */
         cutcoef = transcoefs[i] * valterms[j] / fvalue;

         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

         cutrhs += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
      }
   }

   /* add terms for v_n */
   for( i = termbegins[nterms - 1]; i < termbegins[nterms]; ++i )
   {
      cutvar = SCIPgetExprAuxVarNonlinear(vars[transcoefsidx[i]]);
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, -transcoefs[i]) );
   }

   /* add side */
   SCIProwprepAddSide(rowprep, cutrhs);

   SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIPinfinity(scip), NULL) );

   if( SCIPgetRowprepViolation(scip, rowprep, sol, NULL) >= mincutviolation )
   {
      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "soc%d_%p_%" SCIP_LONGINT_FORMAT, nterms, (void*) expr, SCIPgetNLPs(scip));
      SCIP_CALL( SCIPgetRowprepRowCons(scip, cut, rowprep, cons) );
   }
   else
   {
      SCIPdebugMsg(scip, "%d-SOC rowprep violation %g below mincutviolation %g\n", nterms, SCIPgetRowprepViolation(scip,
               rowprep, sol, NULL), mincutviolation);
      /* SCIPprintRowprep(scip, rowprep, NULL); */
   }

   /* free memory */
   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** helper method to compute and add a gradient cut for the k-th cone disaggregation
 *
 *  After the SOC constraint \f$\sqrt{\sum_{i = 0}^{n-1} (v_i^T x + \beta_i)^2} \leq v_n^T x + \beta_n\f$
 *  has been disaggregated into the row \f$\sum_{i = 0}^{n-1} y_i \leq v_n^T x + \beta_n\f$ and the smaller SOC constraints
 *  \f[
 *    (v_i^T x + \beta_i)^2 \leq (v_n^T x + \beta_n) y_i \text{ for } i \in \{0, \ldots, n -1\},
 *  \f]
 *  we want to separate one of the small rotated cones.
 *  We first transform it into standard form:
 *  \f[
 *    \sqrt{4(v_i^T x + \beta_i)^2 + (v_n^T x + \beta_n - y_i)^2} - v_n^T x - \beta_n - y_i \leq 0.
 *  \f]
 *  Let \f$f(x,y)\f$ be the left-hand-side. We now compute the gradient by
 *  \f{align*}{
 *    \frac{\delta f}{\delta x_j} &= \frac{(v_i)_j(4v_i^T x + 4\beta_i) + (v_n)_j(v_n^T x + \beta_n - y_i)}{\sqrt{4(v_i^T x + \beta_i)^2 + (v_n^T x + \beta_n - y_i)^2}} - (v_n)_j \\
 *    \frac{\delta f}{\delta y_i} &= \frac{y_i - v_n^T x -\beta_n}{\sqrt{4(v_i^T x + \beta_i)^2 + (v_n^T x + \beta_n - y_i)^2}} - 1
 *  \f}
 *  and the gradient cut is then \f$f(x^*, y^*) + \nabla f(x^*,y^*)((x,y) - (x^*, y^*)) \leq 0\f$.
 */
static
SCIP_RETCODE generateCutSolDisagg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_CONS*            cons,               /**< the constraint that expr is part of */
   SCIP_SOL*             sol,                /**< solution to separate or NULL for the LP solution */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   int                   disaggidx,          /**< index of disaggregation to separate */
   SCIP_Real             mincutviolation,    /**< minimal required cut violation */
   SCIP_Real             rhsval,             /**< value of the rhs term */
   SCIP_ROW**            cut                 /**< pointer to store a cut */
   )
{
   SCIP_EXPR** vars;
   SCIP_VAR** disvars;
   SCIP_Real* transcoefs;
   int* transcoefsidx;
   int* termbegins;
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* cutvar;
   SCIP_Real cutcoef;
   SCIP_Real fvalue;
   SCIP_Real disvarval;
   SCIP_Real lhsval;
   SCIP_Real constant;
   SCIP_Real denominator;
   int ncutvars;
   int nterms;
   int i;

   assert(expr != NULL);
   assert(cons != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(disaggidx < nlhdlrexprdata->nterms);
   assert(mincutviolation >= 0.0);
   assert(cut != NULL);

   vars = nlhdlrexprdata->vars;
   disvars = nlhdlrexprdata->disvars;
   transcoefs = nlhdlrexprdata->transcoefs;
   transcoefsidx = nlhdlrexprdata->transcoefsidx;
   termbegins = nlhdlrexprdata->termbegins;
   nterms = nlhdlrexprdata->nterms;

   /* nterms is equal to n in the description and disaggidx is in {0, ..., n - 1} */

   *cut = NULL;

   disvarval = SCIPgetSolVal(scip, sol, disvars[disaggidx]);

   lhsval = evalSingleTerm(scip, nlhdlrexprdata, sol, disaggidx);

   denominator = SQRT(4.0 * SQR(lhsval) + SQR(rhsval - disvarval));

   /* compute value of function to be separated (f(x*,y*)) */
   fvalue = denominator - rhsval - disvarval;

   /* if the disagg soc is not violated don't compute cut */
   if( fvalue <= mincutviolation )
   {
      SCIPdebugMsg(scip, "skip cut on disaggregation index %d as violation=%g below minviolation %g\n", disaggidx,
            fvalue, mincutviolation);
      return SCIP_OKAY;
   }

   /* if the denominator is 0 -> the constraint can't be violated */
   assert(!SCIPisZero(scip, denominator));
   /* if v_disaggidx^T x + beta_disaggidx is 0 -> the constraint can't be violated */
   assert(!SCIPisZero(scip, lhsval));

   /* compute upper bound on the number of variables in cut: vars in rhs + vars in term + disagg var */
   ncutvars = (termbegins[nterms] - termbegins[nterms-1]) + (termbegins[disaggidx + 1] - termbegins[disaggidx]) + 1;

   /* create cut */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, ncutvars) );

   /* constant will be grad_f(x*,y*)^T  (x*, y*) */
   constant = 0.0;

   /* a variable could appear on the lhs and rhs, but we add the coefficients separately  */

   /* add terms for v_disaggidx */
   for( i = termbegins[disaggidx]; i < termbegins[disaggidx + 1]; ++i )
   {
      cutvar = SCIPgetExprAuxVarNonlinear(vars[transcoefsidx[i]]);
      assert(cutvar != NULL);

      /* cutcoef is (the first part of) the partial derivative w.r.t cutvar */
      cutcoef = 4.0 * lhsval * transcoefs[i] / denominator;

      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

      constant += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
   }

   /* add terms for v_n */
   for( i = termbegins[nterms - 1]; i < termbegins[nterms]; ++i )
   {
      cutvar = SCIPgetExprAuxVarNonlinear(vars[transcoefsidx[i]]);
      assert(cutvar != NULL);

      /* cutcoef is the (second part of) the partial derivative w.r.t cutvar */
      cutcoef = (rhsval - disvarval) * transcoefs[i] / denominator - transcoefs[i];

      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

      constant += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
   }

   /* add term for disvar: cutcoef is the the partial derivative w.r.t. the disaggregation variable */
   cutcoef = (disvarval - rhsval) / denominator - 1.0;
   cutvar = disvars[disaggidx];

   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

   constant += cutcoef * SCIPgetSolVal(scip, sol, cutvar);

   /* add side */
   SCIProwprepAddSide(rowprep, constant - fvalue);

   SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIPinfinity(scip), NULL) );

   if( SCIPgetRowprepViolation(scip, rowprep, sol, NULL) >= mincutviolation )
   {
      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "soc_%p_%d_%" SCIP_LONGINT_FORMAT, (void*) expr, disaggidx, SCIPgetNLPs(scip));
      SCIP_CALL( SCIPgetRowprepRowCons(scip, cut, rowprep, cons) );
   }
   else
   {
      SCIPdebugMsg(scip, "rowprep violation %g below mincutviolation %g\n", SCIPgetRowprepViolation(scip, rowprep, sol,
               NULL), mincutviolation);
      /* SCIPprintRowprep(scip, rowprep, NULL); */
   }

   /* free memory */
   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** checks if an expression is quadratic and collects all occurring expressions
 *
 * @pre `expr2idx` and `occurringexprs` need to be initialized with capacity 2 * nchildren
 *
 * @note We assume that a linear term always appears before its corresponding
 * quadratic term in quadexpr; this should be ensured by canonicalize
 */
static
SCIP_RETCODE checkAndCollectQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            quadexpr,           /**< candidate for a quadratic expression */
   SCIP_HASHMAP*         expr2idx,           /**< hashmap to store expressions */
   SCIP_EXPR**           occurringexprs,     /**< array to store expressions */
   int*                  nexprs,             /**< buffer to store number of expressions */
   SCIP_Bool*            success             /**< buffer to store whether the check was successful */
   )
{
   SCIP_EXPR** children;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(quadexpr != NULL);
   assert(expr2idx != NULL);
   assert(occurringexprs != NULL);
   assert(nexprs != NULL);
   assert(success != NULL);

   *nexprs = 0;
   *success = FALSE;
   children = SCIPexprGetChildren(quadexpr);
   nchildren = SCIPexprGetNChildren(quadexpr);

   /* iterate in reverse order to ensure that quadratic terms are found before linear terms */
   for( i = nchildren - 1; i >= 0; --i )
   {
      SCIP_EXPR* child;

      child = children[i];
      if( SCIPisExprPower(scip, child) )
      {
         SCIP_EXPR* childarg;

         if( SCIPgetExponentExprPow(child) != 2.0 )
            return SCIP_OKAY;

         childarg = SCIPexprGetChildren(child)[0];

         if( !SCIPhashmapExists(expr2idx, (void*) childarg) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) childarg, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = childarg;

            ++(*nexprs);
         }
      }
      else if( SCIPisExprVar(scip, child) && SCIPvarIsBinary(SCIPgetVarExprVar(child)) )
      {
         if( !SCIPhashmapExists(expr2idx, (void*) child) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) child, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = child;

            ++(*nexprs);
         }
      }
      else if( SCIPisExprProduct(scip, child) )
      {
         SCIP_EXPR* childarg1;
         SCIP_EXPR* childarg2;

         if( SCIPexprGetNChildren(child) != 2 )
            return SCIP_OKAY;

         childarg1 = SCIPexprGetChildren(child)[0];
         childarg2 = SCIPexprGetChildren(child)[1];

         if( !SCIPhashmapExists(expr2idx, (void*) childarg1) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) childarg1, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = childarg1;

            ++(*nexprs);
         }

         if( !SCIPhashmapExists(expr2idx, (void*) childarg2) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) childarg2, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = childarg2;

            ++(*nexprs);
         }
      }
      else
      {
         /* if there is a linear term without corresponding quadratic term, it is not a SOC */
         if( !SCIPhashmapExists(expr2idx, (void*) child) )
            return SCIP_OKAY;
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/* builds the constraint defining matrix and vector of a quadratic expression
 *
 * @pre `quadmatrix` and `linvector` need to be initialized with size `nexprs`^2 and `nexprs`, resp.
 */
static
void buildQuadExprMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            quadexpr,           /**< the quadratic expression */
   SCIP_HASHMAP*         expr2idx,           /**< hashmap mapping the occurring expressions to their index */
   int                   nexprs,             /**< number of occurring expressions */
   SCIP_Real*            quadmatrix,         /**< pointer to store (the lower-left triangle of) the quadratic matrix */
   SCIP_Real*            linvector           /**< pointer to store the linear vector */
   )
{
   SCIP_EXPR** children;
   SCIP_Real* childcoefs;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(quadexpr != NULL);
   assert(expr2idx != NULL);
   assert(quadmatrix != NULL);
   assert(linvector != NULL);

   children = SCIPexprGetChildren(quadexpr);
   nchildren = SCIPexprGetNChildren(quadexpr);
   childcoefs = SCIPgetCoefsExprSum(quadexpr);

   /* iterate over children to build the constraint defining matrix and vector */
   for( i = 0; i < nchildren; ++i )
   {
      int varpos;

      if( SCIPisExprPower(scip, children[i]) )
      {
         assert(SCIPgetExponentExprPow(children[i]) == 2.0);
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPexprGetChildren(children[i])[0]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPexprGetChildren(children[i])[0]);
         assert(0 <= varpos && varpos < nexprs);

         quadmatrix[varpos * nexprs + varpos] = childcoefs[i];
      }
      else if( SCIPisExprVar(scip, children[i]) && SCIPvarIsBinary(SCIPgetVarExprVar(children[i])) )
      {
         assert(SCIPhashmapExists(expr2idx, (void*) children[i]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);
         assert(0 <= varpos && varpos < nexprs);

         quadmatrix[varpos * nexprs + varpos] = childcoefs[i];
      }
      else if( SCIPisExprProduct(scip, children[i]) )
      {
         int varpos2;

         assert(SCIPexprGetNChildren(children[i]) == 2);
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPexprGetChildren(children[i])[0]));
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPexprGetChildren(children[i])[1]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPexprGetChildren(children[i])[0]);
         assert(0 <= varpos && varpos < nexprs);

         varpos2 = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPexprGetChildren(children[i])[1]);
         assert(0 <= varpos2 && varpos2 < nexprs);
         assert(varpos != varpos2);

         /* Lapack uses only the lower left triangle of the symmetric matrix */
         quadmatrix[MIN(varpos, varpos2) * nexprs + MAX(varpos, varpos2)] = childcoefs[i] / 2.0;
      }
      else
      {
         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);
         assert(0 <= varpos && varpos < nexprs);

         linvector[varpos] = childcoefs[i];
      }
   }
}

/** tries to fill the nlhdlrexprdata for a potential quadratic SOC expression
 *
 * We say "try" because the expression might still turn out not to be a SOC at this point.
 */
static
SCIP_RETCODE tryFillNlhdlrExprDataQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           occurringexprs,     /**< array of all occurring expressions (nvars many) */
   SCIP_Real*            eigvecmatrix,       /**< array containing the Eigenvectors */
   SCIP_Real*            eigvals,            /**< array containing the Eigenvalues */
   SCIP_Real*            bp,                 /**< product of linear vector b * P (eigvecmatrix^t) */
   int                   nvars,              /**< number of variables */
   int*                  termbegins,         /**< pointer to store the termbegins */
   SCIP_Real*            transcoefs,         /**< pointer to store the transcoefs */
   int*                  transcoefsidx,      /**< pointer to store the transcoefsidx */
   SCIP_Real*            offsets,            /**< pointer to store the offsets */
   SCIP_Real*            lhsconstant,        /**< pointer to store the lhsconstant */
   int*                  nterms,             /**< pointer to store the total number of terms */
   SCIP_Bool*            success             /**< whether the expression is indeed a SOC */
   )
{
   SCIP_Real sqrteigval;
   int nextterm = 0;
   int nexttranscoef = 0;
   int specialtermidx;
   int i;
   int j;

   assert(scip != NULL);
   assert(occurringexprs != NULL);
   assert(eigvecmatrix != NULL);
   assert(eigvals != NULL);
   assert(bp != NULL);
   assert(termbegins != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(offsets != NULL);
   assert(lhsconstant != NULL);
   assert(success != NULL);

   *success = FALSE;
   *nterms = 0;

   /* we have lhsconstant + x^t A x + b x <= 0 and A has a single negative eigenvalue; try to build soc;
    * we now store all the v_i^T x + beta_i on the lhs, and compute the constant
    */
   specialtermidx = -1;
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPisZero(scip, eigvals[i]) )
         continue;

      if( eigvals[i] < 0.0 )
      {
         assert(specialtermidx == -1); /* there should only be one negative eigenvalue */

         specialtermidx = i;

         *lhsconstant -= bp[i] * bp[i] / (4.0 * eigvals[i]);

         continue;
      }

      assert(eigvals[i] > 0.0);
      sqrteigval = SQRT(eigvals[i]);

      termbegins[nextterm] = nexttranscoef;
      offsets[nextterm] = bp[i] / (2.0 * sqrteigval);
      *lhsconstant -= bp[i] * bp[i] / (4.0 * eigvals[i]);

      /* set transcoefs */
      for( j = 0; j < nvars; ++j )
      {
         if( !SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
         {
            transcoefs[nexttranscoef] = sqrteigval * eigvecmatrix[i * nvars + j];
            transcoefsidx[nexttranscoef] = j;

            ++nexttranscoef;
         }
      }
      ++nextterm;
   }
   assert(specialtermidx > -1);

   /* process constant; if constant is negative -> no soc */
   if( SCIPisNegative(scip, *lhsconstant) )
      return SCIP_OKAY;

   /* we need lhsconstant to be >= 0 */
   if( *lhsconstant < 0.0 )
      *lhsconstant = 0.0;

   /* store constant term */
   if( *lhsconstant > 0.0 )
   {
      termbegins[nextterm] = nexttranscoef;
      offsets[nextterm] = SQRT( *lhsconstant );
      ++nextterm;
   }

   /* now process rhs */
   {
      SCIP_Real rhstermlb;
      SCIP_Real rhstermub;
      SCIP_Real signfactor;

      assert(-eigvals[specialtermidx] > 0.0);
      sqrteigval = SQRT(-eigvals[specialtermidx]);

      termbegins[nextterm] = nexttranscoef;
      offsets[nextterm] = -bp[specialtermidx] / (2.0 * sqrteigval);

      /* the expression can only be an soc if the resulting rhs term does not change sign;
       * the rhs term is a linear combination of variables, so estimate its bounds
       */
      rhstermlb = offsets[nextterm];
      for( j = 0; j < nvars; ++j )
      {
         SCIP_INTERVAL activity;
         SCIP_Real aux;

         if( SCIPisZero(scip, eigvecmatrix[specialtermidx * nvars + j]) )
            continue;

         SCIP_CALL( SCIPevalExprActivity(scip, occurringexprs[j]) );
         activity = SCIPexprGetActivity(occurringexprs[j]);

         if( eigvecmatrix[specialtermidx * nvars + j] > 0.0 )
         {
            aux = activity.inf;
            assert(!SCIPisInfinity(scip, aux));
         }
         else
         {
            aux = activity.sup;
            assert(!SCIPisInfinity(scip, -aux));
         }

         if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
         {
            rhstermlb = -SCIPinfinity(scip);
            break;
         }
         else
            rhstermlb += sqrteigval * eigvecmatrix[specialtermidx * nvars + j] * aux;
      }

      rhstermub = offsets[nextterm];
      for( j = 0; j < nvars; ++j )
      {
         SCIP_INTERVAL activity;
         SCIP_Real aux;

         if( SCIPisZero(scip, eigvecmatrix[specialtermidx * nvars + j]) )
            continue;

         SCIP_CALL( SCIPevalExprActivity(scip, occurringexprs[j]) );
         activity = SCIPexprGetActivity(occurringexprs[j]);

         if( eigvecmatrix[specialtermidx * nvars + j] > 0.0 )
         {
            aux = activity.sup;
            assert(!SCIPisInfinity(scip, -aux));
         }
         else
         {
            aux = activity.inf;
            assert(!SCIPisInfinity(scip, aux));
         }

         if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
         {
            rhstermub = SCIPinfinity(scip);
            break;
         }
         else
            rhstermub += sqrteigval * eigvecmatrix[specialtermidx * nvars + j] * aux;
      }

      /* since we are just interested in obtaining an interval that contains the real bounds
       * and is tight enough so that we can identify that the rhsvar does not change sign,
       * we swap the bounds in case of numerical troubles
       */
      if( rhstermub < rhstermlb )
      {
         assert(SCIPisEQ(scip, rhstermub, rhstermlb));
         SCIPswapReals(&rhstermub, &rhstermlb);
      }

      /* if rhs changes sign -> not a SOC */
      if( SCIPisLT(scip, rhstermlb, 0.0) && SCIPisGT(scip, rhstermub, 0.0) )
         return SCIP_OKAY;

      signfactor = SCIPisLE(scip, rhstermub, 0.0) ? -1.0 : 1.0;

      offsets[nextterm] *= signfactor;

      /* set transcoefs for rhs term */
      for( j = 0; j < nvars; ++j )
      {
         if( SCIPisZero(scip, eigvecmatrix[specialtermidx * nvars + j]) )
            continue;

         transcoefs[nexttranscoef] = signfactor * sqrteigval * eigvecmatrix[specialtermidx * nvars + j];
         transcoefsidx[nexttranscoef] = j;

         ++nexttranscoef;
      }

      /* if rhs is a constant this method shouldn't have been called */
      assert(nexttranscoef > termbegins[nextterm]);

      /* finish processing term */
      ++nextterm;
   }

   *nterms = nextterm;

   /* sentinel value */
   termbegins[nextterm] = nexttranscoef;

   *success = TRUE;

   return SCIP_OKAY;
}

/** detects if expr &le; auxvar is of the form SQRT(sum_i coef_i (expr_i + shift_i)^2 + const) &le; auxvar
 *
 * @note if a user inputs the above expression with `const` = -epsilon, then `const` is going to be set to 0.
 */
static
SCIP_RETCODE detectSocNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_EXPR** children;
   SCIP_EXPR* child;
   SCIP_EXPR** vars;
   SCIP_HASHMAP* expr2idx;
   SCIP_HASHSET* linexprs;
   SCIP_Real* childcoefs;
   SCIP_Real* offsets;
   SCIP_Real* transcoefs;
   SCIP_Real constant;
   SCIP_Bool issoc;
   int* transcoefsidx;
   int* termbegins;
   int nchildren;
   int nterms;
   int nvars;
   int nextentry;
   int i;

   assert(expr != NULL);
   assert(success != NULL);

   *success = FALSE;
   issoc = TRUE;

   /* relation is not "<=" -> skip */
   if( SCIPgetExprNLocksPosNonlinear(expr) == 0 )
      return SCIP_OKAY;

   assert(SCIPexprGetNChildren(expr) > 0);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /* check whether expression is a SQRT and has a sum as child with at least 2 children and a non-negative constant */
   if( ! SCIPisExprPower(scip, expr)
      || SCIPgetExponentExprPow(expr) != 0.5
      || !SCIPisExprSum(scip, child)
      || SCIPexprGetNChildren(child) < 2
      || SCIPgetConstantExprSum(child) < 0.0)
   {
      return SCIP_OKAY;
   }

   /* assert(SCIPvarGetLbLocal(auxvar) >= 0.0); */

   /* get children of the sum */
   children = SCIPexprGetChildren(child);
   nchildren = SCIPexprGetNChildren(child);
   childcoefs = SCIPgetCoefsExprSum(child);

   /* TODO: should we initialize the hashmap with size SCIPgetNVars() so that it never has to be resized? */
   SCIP_CALL( SCIPhashmapCreate(&expr2idx, SCIPblkmem(scip), nchildren) );
   SCIP_CALL( SCIPhashsetCreate(&linexprs, SCIPblkmem(scip), nchildren) );

   /* we create coefs array here already, since we have to fill it in first loop in case of success
    * +1 for auxvar
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, nchildren+1) );

   nterms = 0;

   /* check if all children are squares or linear terms with matching square term:
    * if the i-th child is (pow, expr, 2) we store the association <|expr -> i|> in expr2idx and if expr was in
    * linexprs, we remove it from there.
    * if the i-th child is expr' (different from (pow, expr, 2)) and expr' is not a key of expr2idx, we add it
    * to linexprs.
    * if at the end there is any expr in linexpr -> we do not have a separable quadratic function.
    */
   for( i = 0; i < nchildren; ++i )
   {
      /* handle quadratic expressions children */
      if( SCIPisExprPower(scip, children[i]) && SCIPgetExponentExprPow(children[i]) == 2.0 )
      {
         SCIP_EXPR* squarearg = SCIPexprGetChildren(children[i])[0];

         if( !SCIPhashmapExists(expr2idx, (void*) squarearg) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void *) squarearg, nterms) );
         }

         if( childcoefs[i] < 0.0 )
         {
            issoc = FALSE;
            break;
         }
         transcoefs[nterms] = SQRT(childcoefs[i]);

         SCIP_CALL( SCIPhashsetRemove(linexprs, (void*) squarearg) );
         ++nterms;
      }
      /* handle binary variable children */
      else if( SCIPisExprVar(scip, children[i]) && SCIPvarIsBinary(SCIPgetVarExprVar(children[i])) )
      {
         assert(!SCIPhashmapExists(expr2idx, (void*) children[i]));
         assert(!SCIPhashsetExists(linexprs, (void*) children[i]));

         SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void *) children[i], nterms) );

         if( childcoefs[i] < 0.0 )
         {
            issoc = FALSE;
            break;
         }
         transcoefs[nterms] = SQRT(childcoefs[i]);

         ++nterms;
      }
      else
      {
         if( !SCIPhashmapExists(expr2idx, (void*) children[i]) )
         {
            SCIP_CALL( SCIPhashsetInsert(linexprs, SCIPblkmem(scip), (void*) children[i]) );
         }
      }
   }

   /* there are linear terms without corresponding quadratic terms or it was detected not to be soc */
   if( SCIPhashsetGetNElements(linexprs) > 0 || ! issoc )
   {
      SCIPfreeBufferArray(scip, &transcoefs);
      SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   /* add one to terms counter for auxvar */
   ++nterms;

   constant = SCIPgetConstantExprSum(child);

   /* compute constant of possible soc expression to check its sign */
   for( i = 0; i < nchildren; ++i )
   {
      if( ! SCIPisExprPower(scip, children[i]) || SCIPgetExponentExprPow(children[i]) != 2.0 )
      {
         int auxvarpos;

         assert(SCIPhashmapExists(expr2idx, (void*) children[i]) );
         auxvarpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);

         constant -= SQR(0.5 * childcoefs[i] / transcoefs[auxvarpos]);
      }
   }

   /* if the constant is negative -> no SOC */
   if( SCIPisNegative(scip, constant) )
   {
      SCIPfreeBufferArray(scip, &transcoefs);
      SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }
   else if( SCIPisZero(scip, constant) )
      constant = 0.0;
   assert(constant >= 0.0);

   /* at this point, we have found an SOC structure */
   *success = TRUE;

   nvars = nterms;

   /* add one to terms counter for constant term */
   if( constant > 0.0 )
      ++nterms;

   /* allocate temporary memory to collect data */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, nterms + 1) );

   /* fill in data for non constant terms of lhs; initialize their offsets */
   for( i = 0; i < nvars - 1; ++i )
   {
      transcoefsidx[i] = i;
      termbegins[i] = i;
      offsets[i] = 0.0;
   }

   /* add constant term and rhs */
   vars[nvars - 1] = expr;
   if( constant > 0.0 )
   {
      /* constant term */
      termbegins[nterms - 2] = nterms - 2;
      offsets[nterms - 2] = SQRT(constant);

      /* rhs */
      termbegins[nterms - 1] = nterms - 2;
      offsets[nterms - 1] = 0.0;
      transcoefsidx[nterms - 2] = nvars - 1;
      transcoefs[nterms - 2] = 1.0;

      /* sentinel value */
      termbegins[nterms] = nterms - 1;
   }
   else
   {
      /* rhs */
      termbegins[nterms - 1] = nterms - 1;
      offsets[nterms - 1] = 0.0;
      transcoefsidx[nterms - 1] = nvars - 1;
      transcoefs[nterms - 1] = 1.0;

      /* sentinel value */
      termbegins[nterms] = nterms;
   }

   /* request required auxiliary variables and fill vars and offsets array */
   nextentry = 0;
   for( i = 0; i < nchildren; ++i )
   {
      if( SCIPisExprPower(scip, children[i]) && SCIPgetExponentExprPow(children[i]) == 2.0 )
      {
         SCIP_EXPR* squarearg;

         squarearg = SCIPexprGetChildren(children[i])[0];
         assert(SCIPhashmapGetImageInt(expr2idx, (void*) squarearg) == nextentry);

         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, squarearg, TRUE, FALSE, FALSE, FALSE) );

         vars[nextentry] = squarearg;
         ++nextentry;
      }
      else if( SCIPisExprVar(scip, children[i]) && SCIPvarIsBinary(SCIPgetVarExprVar(children[i])) )
      {
         /* handle binary variable children: no need to request auxvar */
         assert(SCIPhashmapGetImageInt(expr2idx, (void*) children[i]) == nextentry);
         vars[nextentry] = children[i];
         ++nextentry;
      }
      else
      {
         int auxvarpos;

         assert(SCIPhashmapExists(expr2idx, (void*) children[i]));
         auxvarpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);

         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, children[i], TRUE, FALSE, FALSE, FALSE) );

         offsets[auxvarpos] = 0.5 * childcoefs[i] / transcoefs[auxvarpos];
      }
   }
   assert(nextentry == nvars - 1);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n", (void*)expr);
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, " <= auxvar\n");
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, offsets, transcoefs, transcoefsidx, termbegins,
            nvars, nterms, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

   /* free memory */
   SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
   SCIPhashmapFree(&expr2idx);
   SCIPfreeBufferArray(scip, &termbegins);
   SCIPfreeBufferArray(scip, &transcoefsidx);
   SCIPfreeBufferArray(scip, &offsets);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &transcoefs);

   return SCIP_OKAY;
}

/** helper method to detect c + sum_i coef_i expr_i^2 - coef_k expr_k^2 &le; 0
 *  and c + sum_i coef_i expr_i^2 - coef_k expr_k expr_l &le; 0
 *
 *  binary linear variables are interpreted as quadratic terms
 *
 *  @todo: extend this function to detect  c + sum_i coef_i (expr_i + const_i)^2 - ...
 *  this would probably share a lot of code with detectSocNorm
 */
static
SCIP_RETCODE detectSocQuadraticSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            enforcebelow,       /**< pointer to store whether we enforce <= (TRUE) or >= (FALSE); only valid when success is TRUE */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_EXPR** children;
   SCIP_EXPR** vars = NULL;
   SCIP_Real* childcoefs;
   SCIP_Real* offsets = NULL;
   SCIP_Real* transcoefs = NULL;
   int* transcoefsidx = NULL;
   int* termbegins = NULL;
   SCIP_Real constant;
   SCIP_Real lhsconstant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real rhssign;
   SCIP_INTERVAL expractivity;
   int ntranscoefs;
   int nposquadterms;
   int nnegquadterms;
   int nposbilinterms;
   int nnegbilinterms;
   int rhsidx;
   int lhsidx;
   int specialtermidx;
   int nchildren;
   int nnzinterms;
   int nterms;
   int nvars;
   int nextentry;
   int i;
   SCIP_Bool ishyperbolic;

   assert(expr != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is a sum with at least 2 children */
   if( ! SCIPisExprSum(scip, expr) || SCIPexprGetNChildren(expr) < 2 )
      return SCIP_OKAY;

   /* get children of the sum */
   children = SCIPexprGetChildren(expr);
   nchildren = SCIPexprGetNChildren(expr);
   constant = SCIPgetConstantExprSum(expr);

   /* we duplicate the child coefficients since we have to manipulate them */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &childcoefs, SCIPgetCoefsExprSum(expr), nchildren) ); /*lint !e666*/

   /* initialize data */
   lhsidx = -1;
   rhsidx = -1;
   nposquadterms = 0;
   nnegquadterms = 0;
   nposbilinterms = 0;
   nnegbilinterms = 0;

   /* check if all children are quadratic or binary linear and count number of positive and negative terms */
   for( i = 0; i < nchildren; ++i )
   {
      if( SCIPisExprPower(scip, children[i]) && SCIPgetExponentExprPow(children[i]) == 2.0 )
      {
         if( childcoefs[i] > 0.0 )
         {
            ++nposquadterms;
            lhsidx = i;
         }
         else
         {
            ++nnegquadterms;
            rhsidx = i;
         }
      }
      else if( SCIPisExprVar(scip, children[i]) && SCIPvarIsBinary(SCIPgetVarExprVar(children[i])) )
      {
         if( childcoefs[i] > 0.0 )
         {
            ++nposquadterms;
            lhsidx = i;
         }
         else
         {
            ++nnegquadterms;
            rhsidx = i;
         }
      }
      else if( SCIPisExprProduct(scip, children[i]) && SCIPexprGetNChildren(children[i]) == 2 )
      {
         if( childcoefs[i] > 0.0 )
         {
            ++nposbilinterms;
            lhsidx = i;
         }
         else
         {
            ++nnegbilinterms;
            rhsidx = i;
         }
      }
      else
      {
         goto CLEANUP;
      }

      /* more than one positive eigenvalue and more than one negative eigenvalue -> can't be convex */
      if( nposquadterms > 1 && nnegquadterms > 1 )
         goto CLEANUP;

      /* more than one bilinear term -> can't be handled by this method */
      if( nposbilinterms + nnegbilinterms > 1 )
         goto CLEANUP;

      /* one positive bilinear term and also at least one positive quadratic term -> not a simple SOC */
      if( nposbilinterms > 0 && nposquadterms > 0 )
         goto CLEANUP;

      /* one negative bilinear term and also at least one negative quadratic term -> not a simple SOC */
      if( nnegbilinterms > 0 && nnegquadterms > 0 )
         goto CLEANUP;
   }

   if( nposquadterms == nchildren || nnegquadterms == nchildren )
      goto CLEANUP;

   assert(nposquadterms <= 1 || nnegquadterms <= 1);
   assert(nposbilinterms + nnegbilinterms <= 1);
   assert(nposbilinterms == 0 || nposquadterms == 0);
   assert(nnegbilinterms == 0 || nnegquadterms == 0);

   /* if a bilinear term is involved, it is a hyperbolic expression */
   ishyperbolic = (nposbilinterms + nnegbilinterms > 0);

   if( conslhs == SCIP_INVALID || consrhs == SCIP_INVALID )  /*lint !e777*/
   {
      SCIP_CALL( SCIPevalExprActivity(scip, expr) );
      expractivity = SCIPexprGetActivity(expr);

      lhs = (conslhs == SCIP_INVALID ? expractivity.inf : conslhs); /*lint !e777*/
      rhs = (consrhs == SCIP_INVALID ? expractivity.sup : consrhs); /*lint !e777*/
   }
   else
   {
      lhs = conslhs;
      rhs = consrhs;
   }

   /* detect case and store lhs/rhs information */
   if( (ishyperbolic && nnegbilinterms > 0) || (!ishyperbolic && nnegquadterms < 2) )
   {
      /* we have -x*y + z^2 ... -> we want to write  z^2 ... <= x*y;
       * or we have -x^2 + y^2  ... -> we want to write y^2 ... <= x^2;
       * in any case, we need a finite rhs
       */
      assert(nnegbilinterms == 1 || nnegquadterms == 1);
      assert(rhsidx != -1);

      /* if rhs is infinity, it can't be soc
       * TODO: if it can't be soc, then we should enforce the caller so that we do not try the more complex quadratic
       * method
       */
      if( SCIPisInfinity(scip, rhs) )
         goto CLEANUP;

      specialtermidx = rhsidx;
      lhsconstant = constant - rhs;
      *enforcebelow = TRUE; /* enforce expr <= rhs */
   }
   else
   {
      assert(lhsidx != -1);

      /* if lhs is infinity, it can't be soc */
      if( SCIPisInfinity(scip, -lhs) )
         goto CLEANUP;

      specialtermidx = lhsidx;
      lhsconstant = lhs - constant;

      /* negate all coefficients */
      for( i = 0; i < nchildren; ++i )
         childcoefs[i] = -childcoefs[i];
      *enforcebelow = FALSE; /* enforce lhs <= expr */
   }
   assert(childcoefs[specialtermidx] != 0.0);

   if( ishyperbolic )
   {
      SCIP_INTERVAL yactivity;
      SCIP_INTERVAL zactivity;

      assert(SCIPexprGetNChildren(children[specialtermidx]) == 2);

      SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(children[specialtermidx])[0]) );
      yactivity = SCIPexprGetActivity(SCIPexprGetChildren(children[specialtermidx])[0]);

      SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(children[specialtermidx])[1]) );
      zactivity = SCIPexprGetActivity(SCIPexprGetChildren(children[specialtermidx])[1]);

      if( SCIPisNegative(scip, yactivity.inf + zactivity.inf) )
      {
         /* the sum of the expressions in the bilinear term changes sign -> no SOC */
         if( SCIPisPositive(scip, yactivity.sup + zactivity.sup) )
            goto CLEANUP;

         rhssign = -1.0;
      }
      else
         rhssign = 1.0;

      lhsconstant *= 4.0 / -childcoefs[specialtermidx];
   }
   else if( SCIPisExprVar(scip, children[specialtermidx]) )
   {
      /* children[specialtermidx] can be a variable, in which case we treat it as if it is squared */
      rhssign = 1.0;
   }
   else
   {
      SCIP_INTERVAL rhsactivity;

      assert(SCIPexprGetNChildren(children[specialtermidx]) == 1);
      SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(children[specialtermidx])[0]) );
      rhsactivity = SCIPexprGetActivity(SCIPexprGetChildren(children[specialtermidx])[0]);

      if( rhsactivity.inf < 0.0 )
      {
         /* rhs variable changes sign -> no SOC */
         if( rhsactivity.sup > 0.0 )
            goto CLEANUP;

         rhssign = -1.0;
      }
      else
         rhssign = 1.0;
   }

   if( SCIPisNegative(scip, lhsconstant) )
      goto CLEANUP;

   if( SCIPisZero(scip, lhsconstant) )
      lhsconstant = 0.0;

   /*
    * we have found an SOC-representable expression. Now build the nlhdlrexprdata
    *
    * in the non-hyperbolic case, c + sum_i coef_i expr_i^2 - coef_k expr_k^2 <= 0 is transformed to
    * SQRT( c + sum_i coef_i expr_i^2 ) <= coef_k expr_k
    * so there are nchildren many vars, nchildren (+ 1 if c != 0) many terms, nchildren many coefficients in the vs
    * in SOC representation
    *
    * in the hyperbolic case, c + sum_i coef_i expr_i^2 - coef_k expr_k expr_l <= 0 is transformed to
    * SQRT( 4(c + sum_i coef_i expr_i^2) + (expr_k - expr_l)^2 ) <= expr_k + expr_l
    * so there are nchildren + 1many vars, nchildren + 1(+ 1 if c != 0) many terms, nchildren + 3 many coefficients in
    * the vs in SOC representation
    */

   ntranscoefs = ishyperbolic ? nchildren + 3 : nchildren;
   nvars = ishyperbolic ? nchildren + 1 : nchildren;
   nterms = nvars;

   /* constant term */
   if( lhsconstant > 0.0 )
      nterms++;

   /* SOC was detected, allocate temporary memory for data to collect */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, nterms + 1) );

   *success = TRUE;
   nextentry = 0;

   /* collect all the v_i and beta_i */
   nnzinterms = 0;
   for( i = 0; i < nchildren; ++i )
   {
      /* variable and coef for rhs have to be set to the last entry */
      if( i == specialtermidx )
         continue;

      /* extract (unique) variable appearing in term */
      if( SCIPisExprVar(scip, children[i]) )
      {
         vars[nextentry] = children[i];

         assert(SCIPvarIsBinary(SCIPgetVarExprVar(vars[nextentry])));
      }
      else
      {
         assert(SCIPisExprPower(scip, children[i]));

         /* notify that we will require auxiliary variable */
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(children[i])[0], TRUE, FALSE, FALSE, FALSE) );
         vars[nextentry] = SCIPexprGetChildren(children[i])[0];
      }
      assert(vars[nextentry] != NULL);

      /* store v_i and beta_i */
      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = 0.0;

      transcoefsidx[nnzinterms] = nextentry;
      if( ishyperbolic )
      {
         /* we eliminate the coefficient of the bilinear term to arrive at standard form */
         assert(4.0 * childcoefs[i] / -childcoefs[specialtermidx] > 0.0);
         transcoefs[nnzinterms] = SQRT(4.0 * childcoefs[i] / -childcoefs[specialtermidx]);
      }
      else
      {
         assert(childcoefs[i] > 0.0);
         transcoefs[nnzinterms] = SQRT(childcoefs[i]);
      }

      /* finish adding nonzeros */
      ++nnzinterms;

      /* finish processing term */
      ++nextentry;
   }
   assert(nextentry == nchildren - 1);

   /* store term for constant (v_i = 0) */
   if( lhsconstant > 0.0 )
   {
      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = SQRT(lhsconstant);

      /* finish processing term; this term has 0 nonzero thus we do not increase nnzinterms */
      ++nextentry;
   }

   if( !ishyperbolic )
   {
      /* store rhs term */
      if( SCIPisExprVar(scip, children[specialtermidx]) )
      {
         /* this should be the "children[specialtermidx] can be a variable, in which case we treat it as if it is squared" case */
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, children[specialtermidx], TRUE, FALSE, FALSE, FALSE) );
         vars[nvars - 1] = children[specialtermidx];
      }
      else
      {
         assert(SCIPisExprPower(scip, children[specialtermidx]));
         assert(SCIPexprGetChildren(children[specialtermidx]) != NULL);
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(children[specialtermidx])[0], TRUE, FALSE, FALSE, FALSE) );
         vars[nvars - 1] = SCIPexprGetChildren(children[specialtermidx])[0];
      }

      assert(childcoefs[specialtermidx] < 0.0);

      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = 0.0;
      transcoefs[nnzinterms] = rhssign * SQRT(-childcoefs[specialtermidx]);
      transcoefsidx[nnzinterms] = nvars - 1;

      /* finish adding nonzeros */
      ++nnzinterms;

      /* finish processing term */
      ++nextentry;
   }
   else
   {
      /* store last lhs term and rhs term coming from the bilinear term */
      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(children[specialtermidx])[0], TRUE, FALSE, FALSE, FALSE) );
      vars[nvars - 2] = SCIPexprGetChildren(children[specialtermidx])[0];

      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(children[specialtermidx])[1], TRUE, FALSE, FALSE, FALSE) );
      vars[nvars - 1] = SCIPexprGetChildren(children[specialtermidx])[1];

      /* at this point, vars[nvars - 2] = expr_k and vars[nvars - 1] = expr_l;
       * on the lhs we have the term (expr_k - expr_l)^2
       */
      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = 0.0;

      /* expr_k */
      transcoefsidx[nnzinterms] = nvars - 2;
      transcoefs[nnzinterms] = 1.0;
      ++nnzinterms;

      /* - expr_l */
      transcoefsidx[nnzinterms] = nvars - 1;
      transcoefs[nnzinterms] = -1.0;
      ++nnzinterms;

      /* finish processing term */
      ++nextentry;

      /* on rhs we have +/-(expr_k + expr_l) */
      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = 0.0;

      /* rhssing * expr_k */
      transcoefsidx[nnzinterms] = nvars - 2;
      transcoefs[nnzinterms] = rhssign;
      ++nnzinterms;

      /* rhssing * expr_l */
      transcoefsidx[nnzinterms] = nvars - 1;
      transcoefs[nnzinterms] = rhssign;
      ++nnzinterms;

      /* finish processing term */
      ++nextentry;
   }
   assert(nextentry == nterms);
   assert(nnzinterms == ntranscoefs);

   /* sentinel value */
   termbegins[nextentry] = nnzinterms;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n%f <= ", (void*)expr, lhs);
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, "<= %f\n", rhs);
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, offsets, transcoefs, transcoefsidx, termbegins, nvars, nterms,
            nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

CLEANUP:
   SCIPfreeBufferArrayNull(scip, &termbegins);
   SCIPfreeBufferArrayNull(scip, &transcoefsidx);
   SCIPfreeBufferArrayNull(scip, &transcoefs);
   SCIPfreeBufferArrayNull(scip, &offsets);
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &childcoefs);

   return SCIP_OKAY;
}

/** detects complex quadratic expressions that can be represented as SOC constraints
 *
 *  These are quadratic expressions with either exactly one positive or exactly one negative eigenvalue,
 *  in addition to some extra conditions. One needs to write the quadratic as
 *  sum eigval_i (eigvec_i . x)^2 + c &le; -eigval_k (eigvec_k . x)^2, where eigval_k is the negative eigenvalue,
 *  and c must be positive and (eigvec_k . x) must not change sign.
 *  This is described in more details in
 *  Mahajan, Ashutosh & Munson, Todd, Exploiting Second-Order Cone Structure for Global Optimization, 2010.
 *
 *  The eigen-decomposition is computed using Lapack.
 *  Binary linear variables are interpreted as quadratic terms.
 *
 * @todo: In the case -b <= a + x^2 - y^2 <= b, it is possible to represent both sides by SOC. Currently, the
 * datastructure can only handle one SOC. If this should appear more often, it could be worth to extend it,
 * such that both sides can be handled (see e.g. instance chp_partload).
 * FS: this shouldn't be possible. For a <= b + x^2 - y^2 <= c to be SOC representable on both sides, we would need
 * that a - b >= 0 and b -c >= 0, but this implies that a >= c and assuming the constraint is not trivially infeasible,
 * a <= b. Thus, a = b = c and the constraint is x^2 == y^2.
 *
 * @todo: Since cons_nonlinear multiplies as many terms out as possible during presolving, some SOC-representable
 * structures cannot be detected, (see e.g. instances bearing or wager). There is currently no obvious way
 * to handle this.
 */
static
SCIP_RETCODE detectSocQuadraticComplex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            enforcebelow,       /**< pointer to store whether we enforce <= (TRUE) or >= (FALSE); only
                                              *   valid when success is TRUE */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_EXPR** occurringexprs;
   SCIP_HASHMAP* expr2idx;
   SCIP_Real* offsets;
   SCIP_Real* transcoefs;
   SCIP_Real* eigvecmatrix;
   SCIP_Real* eigvals;
   SCIP_Real* lincoefs;
   SCIP_Real* bp;
   int* transcoefsidx;
   int* termbegins;
   SCIP_Real constant;
   SCIP_Real lhsconstant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_INTERVAL expractivity;
   int nvars;
   int nterms;
   int nchildren;
   int npos;
   int nneg;
   int ntranscoefs;
   int i;
   int j;
   SCIP_Bool rhsissoc;
   SCIP_Bool lhsissoc;
   SCIP_Bool isquadratic;

   assert(expr != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is a sum with at least 2 children */
   if( ! SCIPisExprSum(scip, expr) || SCIPexprGetNChildren(expr) < 2 )
   {
      return SCIP_OKAY;
   }

   /* we need Lapack to compute eigenvalues/vectors below */
   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   /* get children of the sum */
   nchildren = SCIPexprGetNChildren(expr);
   constant = SCIPgetConstantExprSum(expr);

   /* initialize data */
   offsets = NULL;
   transcoefs = NULL;
   transcoefsidx = NULL;
   termbegins = NULL;
   bp = NULL;

   SCIP_CALL( SCIPhashmapCreate(&expr2idx, SCIPblkmem(scip), 2 * nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &occurringexprs, 2 * nchildren) );

   /* check if the expression is quadratic and collect all occurring expressions */
   SCIP_CALL( checkAndCollectQuadratic(scip, expr, expr2idx, occurringexprs, &nvars, &isquadratic) );

   if( !isquadratic )
   {
      SCIPfreeBufferArray(scip, &occurringexprs);
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   /* check that nvars*nvars doesn't get too large, see also SCIPcomputeExprQuadraticCurvature() */
   if( nvars > 7000 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "nlhdlr_soc - number of quadratic variables is too large (%d) to check the curvature\n", nvars);
      SCIPfreeBufferArray(scip, &occurringexprs);
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   assert(SCIPhashmapGetNElements(expr2idx) == nvars);

   /* create datastructures for constaint defining matrix and vector */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &eigvecmatrix, nvars * nvars) ); /*lint !e647*/
   SCIP_CALL( SCIPallocClearBufferArray(scip, &lincoefs, nvars) );

   /* build constraint defining matrix (stored in eigvecmatrix) and vector (stored in lincoefs) */
   buildQuadExprMatrix(scip, expr, expr2idx, nvars, eigvecmatrix, lincoefs);

   SCIP_CALL( SCIPallocBufferArray(scip, &eigvals, nvars) );

   /* compute eigenvalues and vectors, A = PDP^t
    * note: eigvecmatrix stores P^t, i.e., P^t_{i,j} = eigvecmatrix[i*nvars+j]
    */
   if( SCIPcallLapackDsyevIpopt(TRUE, nvars, eigvecmatrix, eigvals) != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Failed to compute eigenvalues and eigenvectors for expression:\n");

#ifdef SCIP_DEBUG
      SCIPdismantleExpr(scip, NULL, expr);
#endif

      goto CLEANUP;
   }

   SCIP_CALL( SCIPallocClearBufferArray(scip, &bp, nvars) );

   nneg = 0;
   npos = 0;
   ntranscoefs = 0;

   /* set small eigenvalues to 0 and compute b*P */
   for( i = 0; i < nvars; ++i )
   {
      for( j = 0; j < nvars; ++j )
      {
         bp[i] += lincoefs[j] * eigvecmatrix[i * nvars + j];

         /* count the number of transcoefs to be used later */
         if( !SCIPisZero(scip, eigvals[i]) && !SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
            ++ntranscoefs;
      }

      if( SCIPisZero(scip, eigvals[i]) )
      {
         /* if there is a purely linear variable, the constraint can't be written as a SOC */
         if( !SCIPisZero(scip, bp[i]) )
            goto CLEANUP;

         bp[i] = 0.0;
         eigvals[i] = 0.0;
      }
      else if( eigvals[i] > 0.0 )
         npos++;
      else
         nneg++;
   }

   /* a proper SOC constraint needs at least 2 variables */
   if( npos + nneg < 2 )
      goto CLEANUP;

   /* determine whether rhs or lhs of cons is potentially SOC, if any */
   rhsissoc = (nneg == 1 && SCIPgetExprNLocksPosNonlinear(expr) > 0);
   lhsissoc = (npos == 1 && SCIPgetExprNLocksNegNonlinear(expr) > 0);

   if( rhsissoc || lhsissoc )
   {
      if( conslhs == SCIP_INVALID || consrhs == SCIP_INVALID ) /*lint !e777*/
      {
         SCIP_CALL( SCIPevalExprActivity(scip, expr) );
         expractivity = SCIPexprGetActivity(expr);
         lhs = (conslhs == SCIP_INVALID ? expractivity.inf : conslhs); /*lint !e777*/
         rhs = (consrhs == SCIP_INVALID ? expractivity.sup : consrhs); /*lint !e777*/
      }
      else
      {
         lhs = conslhs;
         rhs = consrhs;
      }
   }
   else
   {
      /* if none of the sides is potentially SOC, stop */
      goto CLEANUP;
   }

   /* @TODO: what do we do if both sides are possible? */
   if( !rhsissoc )
   {
      assert(lhsissoc);

      /* lhs is potentially SOC, change signs */
      lhsconstant = lhs - constant;  /*lint !e644*/

      for( i = 0; i < nvars; ++i )
      {
         eigvals[i] = -eigvals[i];
         bp[i] = -bp[i];
      }
      *enforcebelow = FALSE; /* enforce lhs <= expr */
   }
   else
   {
      lhsconstant = constant - rhs;  /*lint !e644*/
      *enforcebelow = TRUE; /* enforce expr <= rhs */
   }

   /* initialize remaining datastructures for nonlinear handler */
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, npos + nneg + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, npos + nneg + 2) );

   /* try to fill the nlhdlrexprdata (at this point, it can still fail) */
   SCIP_CALL( tryFillNlhdlrExprDataQuad(scip, occurringexprs, eigvecmatrix, eigvals, bp, nvars, termbegins, transcoefs,
         transcoefsidx, offsets, &lhsconstant, &nterms, success) );

   if( !(*success) )
      goto CLEANUP;

   assert(0 < nterms && nterms <= npos + nneg + 1);
   assert(ntranscoefs == termbegins[nterms]);

   /*
    * at this point, the expression passed all checks and is SOC-representable
    */

   /* register all requests for auxiliary variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, occurringexprs[i], TRUE, FALSE, FALSE, FALSE) );
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n%f <= ", (void*)expr, lhs);
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, "<= %f\n", rhs);
#endif

   /* finally, create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, occurringexprs, offsets, transcoefs, transcoefsidx, termbegins, nvars, nterms,
            nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

CLEANUP:
   SCIPfreeBufferArrayNull(scip, &termbegins);
   SCIPfreeBufferArrayNull(scip, &transcoefsidx);
   SCIPfreeBufferArrayNull(scip, &transcoefs);
   SCIPfreeBufferArrayNull(scip, &offsets);
   SCIPfreeBufferArrayNull(scip, &bp);
   SCIPfreeBufferArray(scip, &eigvals);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &eigvecmatrix);
   SCIPfreeBufferArray(scip, &occurringexprs);
   SCIPhashmapFree(&expr2idx);

   return SCIP_OKAY;
}

/** helper method to detect SOC structures
 *
 * The detection runs in 3 steps:
 *  1. check if expression is a norm of the form \f$\sqrt{\sum_i (\text{sqrcoef}_i\, \text{expr}_i^2 + \text{lincoef}_i\, \text{expr}_i) + c}\f$
 *  which can be transformed to the form \f$\sqrt{\sum_i (\text{coef}_i \text{expr}_i + \text{const}_i)^2 + c^*}\f$ with \f$c^* \geq 0\f$.\n
 *    -> this results in the SOC     expr &le; auxvar(expr)
 *
 *    TODO we should generalize and check for sqrt(positive-semidefinite-quadratic)
 *
 *  2. check if expression represents a quadratic function of one of the following forms (all coefs > 0)
 *     1. \f$(\sum_i   \text{coef}_i \text{expr}_i^2) - \text{coef}_k \text{expr}_k^2 \leq \text{RHS}\f$ or
 *     2. \f$(\sum_i - \text{coef}_i \text{expr}_i^2) + \text{coef}_k \text{expr}_k^2 \geq \text{LHS}\f$ or
 *     3. \f$(\sum_i   \text{coef}_i \text{expr}_i^2) - \text{coef}_k \text{expr}_k \text{expr}_l \leq \text{RHS}\f$ or
 *     4. \f$(\sum_i - \text{coef}_i \text{expr}_i^2) + \text{coef}_k \text{expr}_k \text{expr}_l \geq \text{LHS}\f$,
 *
 *     where RHS &ge; 0 or LHS &le; 0, respectively. For LHS and RHS we use the constraint sides if it is a root expr
 *     and the bounds of the auxiliary variable otherwise.
 *     The last two cases are called hyperbolic or rotated second order cone.\n
 *     -> this results in the SOC \f$\sqrt{(\sum_i \text{coef}_i \text{expr}_i^2) - \text{RHS}} \leq \sqrt{\text{coef}_k} \text{expr}_k\f$
 *                            or  \f$\sqrt{4(\sum_i \text{coef}_i \text{expr}_i^2) - 4\text{RHS} + (\text{expr}_k - \text{expr}_l)^2)} \leq \text{expr}_k + \text{expr}_l\f$.
 *                            (analogously for the LHS cases)
 *
 *  3. check if expression represents a quadratic inequality of the form \f$f(x) = x^TAx + b^Tx + c \leq 0\f$ such that \f$f(x)\f$
 *  has exactly one negative eigenvalue plus some extra conditions, see detectSocQuadraticComplex().
 *
 *  Note that step 3 is only performed if parameter `compeigenvalues` is set to TRUE.
 */
static
SCIP_RETCODE detectSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            enforcebelow,       /**< pointer to store whether we enforce <= (TRUE) or >= (FALSE); only
                                              *   valid when success is TRUE */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   assert(expr != NULL);
   assert(nlhdlrdata != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is given as norm as described in case 1 above: if we have a constraint
    * sqrt(sum x_i^2) <= constant, then it might be better not to handle this here; thus, we only call detectSocNorm
    * when the expr is _not_ the root of a constraint
    */
   if( conslhs == SCIP_INVALID && consrhs == SCIP_INVALID ) /*lint !e777*/
   {
      SCIP_CALL( detectSocNorm(scip, expr, nlhdlrexprdata, success) );
      *enforcebelow = *success;
   }

   if( !(*success) )
   {
      /* check whether expression is a simple soc-respresentable quadratic expression as described in case 2 above */
      SCIP_CALL( detectSocQuadraticSimple(scip, expr, conslhs, consrhs, nlhdlrexprdata, enforcebelow, success) );
   }

   if( !(*success) && nlhdlrdata->compeigenvalues )
   {
      /* check whether expression is a more complex soc-respresentable quadratic expression as described in case 3 */
      SCIP_CALL( detectSocQuadraticComplex(scip, expr, conslhs, consrhs, nlhdlrexprdata, enforcebelow, success) );
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrSoc)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrSoc(targetscip) );

   return SCIP_OKAY;
}

/** callback to free data of handler */
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataSoc)
{ /*lint --e{715}*/
   assert(nlhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataSoc)
{  /*lint --e{715}*/
   assert(*nlhdlrexprdata != NULL);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );

   return SCIP_OKAY;
}


/** callback to be called in initialization */
#if 0
static
SCIP_DECL_NLHDLRINIT(nlhdlrInitSoc)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrInitSoc NULL
#endif


/** callback to be called in deinitialization */
#if 0
static
SCIP_DECL_NLHDLREXIT(nlhdlrExitSoc)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrExitSoc NULL
#endif


/** callback to detect structure in expression tree */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectSoc)
{ /*lint --e{715}*/
   SCIP_Real conslhs;
   SCIP_Real consrhs;
   SCIP_Bool enforcebelow;
   SCIP_Bool success;
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(expr != NULL);

   /* don't try if no sepa is required
    * TODO implement some bound strengthening
    */
   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) == SCIP_NLHDLR_METHOD_SEPABOTH )
      return SCIP_OKAY;

   assert(SCIPgetExprNAuxvarUsesNonlinear(expr) > 0);  /* since some sepa is required, there should have been demand for it */

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   conslhs = (cons == NULL ? SCIP_INVALID : SCIPgetLhsNonlinear(cons));
   consrhs = (cons == NULL ? SCIP_INVALID : SCIPgetRhsNonlinear(cons));

   SCIP_CALL( detectSOC(scip, nlhdlrdata, expr, conslhs, consrhs, nlhdlrexprdata, &enforcebelow, &success) );

   if( !success )
      return SCIP_OKAY;

   /* inform what we can do */
   *participating = enforcebelow ? SCIP_NLHDLR_METHOD_SEPABELOW : SCIP_NLHDLR_METHOD_SEPAABOVE;

   /* if we have been successful on sqrt(...) <= auxvar, then we enforce
    * otherwise, expr is quadratic and we separate for expr <= ub(auxvar) only
    * in that case, we enforce only if expr is the root of a constraint, since then replacing auxvar by up(auxvar) does not relax anything (auxvar <= ub(auxvar) is the only constraint on auxvar)
    */
   if( (SCIPisExprPower(scip, expr) && SCIPgetExponentExprPow(expr) == 0.5) || (cons != NULL) )
      *enforcing |= *participating;

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler
 * @todo: remember if we are in the original variables and avoid reevaluating
 */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxSoc)
{ /*lint --e{715}*/
   int i;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->vars != NULL);
   assert(nlhdlrexprdata->transcoefs != NULL);
   assert(nlhdlrexprdata->transcoefsidx != NULL);
   assert(nlhdlrexprdata->nterms > 1);

   /* if the original expression is a norm, evaluate w.r.t. the auxiliary variables */
   if( SCIPisExprPower(scip, expr) )
   {
      assert(SCIPgetExponentExprPow(expr) == 0.5);

      /* compute sum_i coef_i expr_i^2 */
      *auxvalue = 0.0;
      for( i = 0; i < nlhdlrexprdata->nterms - 1; ++i )
      {
         SCIP_Real termval;

         termval = evalSingleTerm(scip, nlhdlrexprdata, sol, i);
         *auxvalue += SQR(termval);
      }

      assert(*auxvalue >= 0.0);

      /* compute SQRT(sum_i coef_i expr_i^2) */
      *auxvalue = SQRT(*auxvalue);
   }
   /* otherwise, evaluate the original quadratic expression w.r.t. the created auxvars of the children */
   else
   {
      SCIP_EXPR** children;
      SCIP_Real* childcoefs;
      int nchildren;

      assert(SCIPisExprSum(scip, expr));

      children = SCIPexprGetChildren(expr);
      childcoefs = SCIPgetCoefsExprSum(expr);
      nchildren = SCIPexprGetNChildren(expr);

      *auxvalue = SCIPgetConstantExprSum(expr);

      for( i = 0; i < nchildren; ++i )
      {
         if( SCIPisExprPower(scip, children[i]) )
         {
            SCIP_VAR* argauxvar;
            SCIP_Real solval;

            assert(SCIPgetExponentExprPow(children[i]) == 2.0);

            argauxvar = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(children[i])[0]);
            assert(argauxvar != NULL);

            solval = SCIPgetSolVal(scip, sol, argauxvar);
            *auxvalue += childcoefs[i] * SQR( solval );
         }
         else if( SCIPisExprProduct(scip, children[i]) )
         {
            SCIP_VAR* argauxvar1;
            SCIP_VAR* argauxvar2;

            assert(SCIPexprGetNChildren(children[i]) == 2);

            argauxvar1 = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(children[i])[0]);
            argauxvar2 = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(children[i])[1]);
            assert(argauxvar1 != NULL);
            assert(argauxvar2 != NULL);

            *auxvalue += childcoefs[i] * SCIPgetSolVal(scip, sol, argauxvar1) * SCIPgetSolVal(scip, sol, argauxvar2);
         }
         else
         {
            SCIP_VAR* argauxvar;

            argauxvar = SCIPgetExprAuxVarNonlinear(children[i]);
            assert(argauxvar != NULL);

            *auxvalue += childcoefs[i] * SCIPgetSolVal(scip, sol, argauxvar);
         }
      }
   }

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSINITLP) */
static
SCIP_DECL_NLHDLRINITSEPA(nlhdlrInitSepaSoc)
{ /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   /* if we have 3 or more terms in lhs create variable and row for disaggregation */
   if( nlhdlrexprdata->nterms > 3 )
   {
      /* create variables for cone disaggregation */
      SCIP_CALL( createDisaggrVars(scip, expr, nlhdlrexprdata) );

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real lhsval;
         SCIP_Real rhsval;
         SCIP_Real disvarval;
         int ndisvars;
         int nterms;
         int i;
         int k;

         /*  the debug solution value of the disaggregation variables is set to
          *      (v_i^T x + beta_i)^2 / (v_{n+1}^T x + beta_{n+1})
          *  if (v_{n+1}^T x + beta_{n+1}) is different from 0.
          *  Otherwise, the debug solution value is set to 0.
          */

         nterms = nlhdlrexprdata->nterms;

         /* find value of rhs */
         rhsval = nlhdlrexprdata->offsets[nterms - 1];
         for( i = nlhdlrexprdata->termbegins[nterms - 1]; i < nlhdlrexprdata->termbegins[nterms]; ++i )
         {
            SCIP_VAR* var;
            SCIP_Real varval;

            var = SCIPgetVarExprVar(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[i]]);

            SCIP_CALL( SCIPdebugGetSolVal(scip, var, &varval) );
            rhsval += nlhdlrexprdata->transcoefs[i] * varval;
         }

         /* set value of disaggregation vars */
         ndisvars = nlhdlrexprdata->nterms - 1;

         if( SCIPisZero(scip, rhsval) )
         {
            for( i = 0; i < ndisvars; ++i )
            {
               SCIP_CALL( SCIPdebugAddSolVal(scip, nlhdlrexprdata->disvars[i], 0.0) );
            }
         }
         else
         {
            /* set value for each disaggregation variable corresponding to quadratic term */
            for( k = 0; k < ndisvars; ++k )
            {
               lhsval = nlhdlrexprdata->offsets[k];

               for( i = nlhdlrexprdata->termbegins[k]; i < nlhdlrexprdata->termbegins[k + 1]; ++i )
               {
                  SCIP_VAR* var;
                  SCIP_Real varval;

                  var = SCIPgetVarExprVar(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[i]]);

                  SCIP_CALL( SCIPdebugGetSolVal(scip, var, &varval) );
                  lhsval += nlhdlrexprdata->transcoefs[i] * varval;
               }

               disvarval = SQR(lhsval) / rhsval;

               SCIP_CALL( SCIPdebugAddSolVal(scip, nlhdlrexprdata->disvars[k], disvarval) );
            }
         }
      }
#endif

      /* create the disaggregation row and store it in nlhdlrexprdata */
      SCIP_CALL( createDisaggrRow(scip, conshdlr, expr, nlhdlrexprdata) );
   }

   /* TODO add something to the LP as well, at least the disaggregation row */

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
static
SCIP_DECL_NLHDLREXITSEPA(nlhdlrExitSepaSoc)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   /* free disaggreagation row */
   if( nlhdlrexprdata->disrow != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &nlhdlrexprdata->disrow) );
   }

   return SCIP_OKAY;
}


/** nonlinear handler separation callback */
static
SCIP_DECL_NLHDLRENFO(nlhdlrEnfoSoc)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_Real rhsval;
   int ndisaggrs;
   int k;
   SCIP_Bool infeasible;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nterms < 4 || nlhdlrexprdata->disrow != NULL);
   assert(nlhdlrexprdata->nterms > 1);

   *result = SCIP_DIDNOTFIND;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   rhsval = evalSingleTerm(scip, nlhdlrexprdata, sol, nlhdlrexprdata->nterms - 1);

   /* if there are three or two terms just compute gradient cut */
   if( nlhdlrexprdata->nterms < 4 )
   {
      SCIP_ROW* row;

      /* compute gradient cut */
      SCIP_CALL( generateCutSolSOC(scip, expr, cons, sol, nlhdlrexprdata, SCIPgetLPFeastol(scip), rhsval, &row) );

      /* TODO this code repeats below, factorize out */
      if( row != NULL )
      {
         SCIP_Real cutefficacy;

         cutefficacy = SCIPgetCutEfficacy(scip, sol, row);

         SCIPdebugMsg(scip, "generated row for %d-SOC, efficacy=%g, minefficacy=%g, allowweakcuts=%u\n",
            nlhdlrexprdata->nterms, cutefficacy, nlhdlrdata->mincutefficacy, allowweakcuts);

         /* check whether cut is applicable */
         if( SCIPisCutApplicable(scip, row) && (allowweakcuts || cutefficacy >= nlhdlrdata->mincutefficacy) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

#ifdef SCIP_CONSNONLINEAR_ROWNOTREMOVABLE
            /* mark row as not removable from LP for current node, if in enforcement (==addbranchscores) (this can prevent some cycling) */
            if( addbranchscores )
               SCIPmarkRowNotRemovableLocal(scip, row);
#endif

            if( infeasible )
               *result = SCIP_CUTOFF;
            else
               *result = SCIP_SEPARATED;
         }

         /* release row */
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
#ifdef SCIP_DEBUG
      else
      {
         SCIPdebugMsg(scip, "failed to generate SOC\n");
      }
#endif

      return SCIP_OKAY;
   }

   ndisaggrs = nlhdlrexprdata->nterms - 1;

   /* check whether the aggregation row is in the LP */
   if( !SCIProwIsInLP(nlhdlrexprdata->disrow) && -SCIPgetRowSolFeasibility(scip, nlhdlrexprdata->disrow, sol) > SCIPgetLPFeastol(scip) )
   {
      SCIP_CALL( SCIPaddRow(scip, nlhdlrexprdata->disrow, TRUE, &infeasible) );
      SCIPdebugMsg(scip, "added disaggregation row to LP, cutoff=%u\n", infeasible);

      if( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      *result = SCIP_SEPARATED;
   }

   for( k = 0; k < ndisaggrs && *result != SCIP_CUTOFF; ++k )
   {
      SCIP_ROW* row;

      /* compute gradient cut */
      SCIP_CALL( generateCutSolDisagg(scip, expr, cons, sol, nlhdlrexprdata, k, SCIPgetLPFeastol(scip), rhsval, &row) );

      if( row != NULL )
      {
         SCIP_Real cutefficacy;

         cutefficacy = SCIPgetCutEfficacy(scip, sol, row);

         SCIPdebugMsg(scip, "generated row for disaggregation %d, efficacy=%g, minefficacy=%g, allowweakcuts=%u\n",
            k, cutefficacy, nlhdlrdata->mincutefficacy, allowweakcuts);

         /* check whether cut is applicable */
         if( SCIPisCutApplicable(scip, row) && (allowweakcuts || cutefficacy >= nlhdlrdata->mincutefficacy) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

#ifdef SCIP_CONSNONLINEAR_ROWNOTREMOVABLE
            /* mark row as not removable from LP for current node, if in enforcement (==addbranchscores) (this can prevent some cycling) */
            if( addbranchscores )
               SCIPmarkRowNotRemovableLocal(scip, row);
#endif

            if( infeasible )
               *result = SCIP_CUTOFF;
            else
               *result = SCIP_SEPARATED;
         }

         /* release row */
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }

   return SCIP_OKAY;
}

/*
 * nonlinear handler specific interface methods
 */

/** includes SOC nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrSoc(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   /* create nonlinear handler data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY, NLHDLR_ENFOPRIORITY, nlhdlrDetectSoc, nlhdlrEvalauxSoc, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrSoc);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrdataSoc);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataSoc);
   SCIPnlhdlrSetInitExit(nlhdlr, nlhdlrInitSoc, nlhdlrExitSoc);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaSoc, nlhdlrEnfoSoc, NULL, nlhdlrExitSepaSoc);

   /* add soc nlhdlr parameters */
   /* TODO should we get rid of this and use separating/mineffiacy(root) instead, which is 1e-4? */
   SCIP_CALL( SCIPaddRealParam(scip, "nlhdlr/" NLHDLR_NAME "/mincutefficacy",
         "Minimum efficacy which a cut needs in order to be added.",
         &nlhdlrdata->mincutefficacy, FALSE, DEFAULT_MINCUTEFFICACY, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/compeigenvalues",
         "Should Eigenvalue computations be done to detect complex cases in quadratic constraints?",
         &nlhdlrdata->compeigenvalues, FALSE, DEFAULT_COMPEIGENVALUES, NULL, NULL) );

   return SCIP_OKAY;
}

/** checks whether constraint is SOC representable in original variables and returns the SOC representation
 *
 * The SOC representation has the form:
 * \f$\sqrt{\sum_{i=1}^{n} (v_i^T x + \beta_i)^2} - v_{n+1}^T x - \beta_{n+1} \lessgtr 0\f$,
 * where \f$n+1 = \text{nterms}\f$ and the inequality type is given by sidetype (`SCIP_SIDETYPE_RIGHT` if inequality
 * is \f$\leq\f$, `SCIP_SIDETYPE_LEFT` if \f$\geq\f$).
 *
 * For each term (i.e. for each \f$i\f$ in the above notation as well as \f$n+1\f$), the constant \f$\beta_i\f$ is given by the
 * corresponding element `offsets[i-1]` and `termbegins[i-1]` is the starting position of the term in arrays
 * `transcoefs` and `transcoefsidx`. The overall number of nonzeros is `termbegins[nterms]`.
 *
 * Arrays `transcoefs` and `transcoefsidx` have size `termbegins[nterms]` and define the linear expressions \f$v_i^T x\f$
 * for each term. For a term \f$i\f$ in the above notation, the nonzeroes are given by elements
 * `termbegins[i-1]...termbegins[i]` of `transcoefs` and `transcoefsidx`. There may be no nonzeroes for some term (i.e.,
 * constant terms are possible). `transcoefs` contains the coefficients \f$v_i\f$ and `transcoefsidx` contains positions of
 * variables in the `vars` array.
 *
 * The `vars` array has size `nvars` and contains \f$x\f$ variables; each variable is included at most once.
 *
 * The arrays should be freed by calling SCIPfreeSOCArraysNonlinear().
 *
 * This function uses the methods that are used in the detection algorithm of the SOC nonlinear handler.
 */
SCIP_RETCODE SCIPisSOCNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_Bool             compeigenvalues,    /**< whether eigenvalues should be computed to detect complex cases */
   SCIP_Bool*            success,            /**< pointer to store whether SOC structure has been detected */
   SCIP_SIDETYPE*        sidetype,           /**< pointer to store which side of cons is SOC representable; only
                                              *   valid when success is TRUE */
   SCIP_VAR***           vars,               /**< variables (x) that appear on both sides; no duplicates are allowed */
   SCIP_Real**           offsets,            /**< offsets of both sides (beta_i) */
   SCIP_Real**           transcoefs,         /**< non-zeros of linear transformation vectors (v_i) */
   int**                 transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int**                 termbegins,         /**< starting indices of transcoefs for each term */
   int*                  nvars,              /**< total number of variables appearing (i.e. size of vars) */
   int*                  nterms              /**< number of summands in the SQRT +1 for RHS (n+1) */
   )
{
   SCIP_NLHDLRDATA nlhdlrdata;
   SCIP_NLHDLREXPRDATA *nlhdlrexprdata;
   SCIP_Real conslhs;
   SCIP_Real consrhs;
   SCIP_EXPR* expr;
   SCIP_Bool enforcebelow;
   int i;

   assert(cons != NULL);

   expr = SCIPgetExprNonlinear(cons);
   assert(expr != NULL);

   nlhdlrdata.mincutefficacy = 0.0;
   nlhdlrdata.compeigenvalues = compeigenvalues;

   conslhs = SCIPgetLhsNonlinear(cons);
   consrhs = SCIPgetRhsNonlinear(cons);

   SCIP_CALL( detectSOC(scip, &nlhdlrdata, expr, conslhs, consrhs, &nlhdlrexprdata, &enforcebelow, success) );

   /* the constraint must be SOC representable in original variables */
   if( *success )
   {
      assert(nlhdlrexprdata != NULL);

      for( i = 0; i < nlhdlrexprdata->nvars; ++i )
      {
         if( !SCIPisExprVar(scip, nlhdlrexprdata->vars[i]) )
         {
            *success = FALSE;
            break;
         }
      }
   }

   if( *success )
   {
      *sidetype = enforcebelow ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, vars, nlhdlrexprdata->nvars) );

      for( i = 0; i < nlhdlrexprdata->nvars; ++i )
      {
         (*vars)[i] = SCIPgetVarExprVar(nlhdlrexprdata->vars[i]);
         assert((*vars)[i] != NULL);
      }
      SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->vars, nlhdlrexprdata->nvars);
      *offsets = nlhdlrexprdata->offsets;
      *transcoefs = nlhdlrexprdata->transcoefs;
      *transcoefsidx = nlhdlrexprdata->transcoefsidx;
      *termbegins = nlhdlrexprdata->termbegins;
      *nvars = nlhdlrexprdata->nvars;
      *nterms = nlhdlrexprdata->nterms;
      SCIPfreeBlockMemory(scip, &nlhdlrexprdata);
   }
   else
   {
      if( nlhdlrexprdata != NULL )
      {
         SCIP_CALL( freeNlhdlrExprData(scip, &nlhdlrexprdata) );
      }
      *vars = NULL;
      *offsets = NULL;
      *transcoefs = NULL;
      *transcoefsidx = NULL;
      *termbegins = NULL;
      *nvars = 0;
      *nterms = 0;
   }

   return SCIP_OKAY;
}

/** frees arrays created by SCIPisSOCNonlinear() */
void SCIPfreeSOCArraysNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variables that appear on both sides (x) */
   SCIP_Real**           offsets,            /**< offsets of both sides (beta_i) */
   SCIP_Real**           transcoefs,         /**< non-zeros of linear transformation vectors (v_i) */
   int**                 transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int**                 termbegins,         /**< starting indices of transcoefs for each term */
   int                   nvars,              /**< total number of variables appearing */
   int                   nterms              /**< number of summands in the SQRT +1 for RHS (n+1) */
   )
{
   int ntranscoefs;

   if( nvars == 0 )
      return;

   assert(vars != NULL);
   assert(offsets != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(termbegins != NULL);

   ntranscoefs = (*termbegins)[nterms];

   SCIPfreeBlockMemoryArray(scip, termbegins, nterms + 1);
   SCIPfreeBlockMemoryArray(scip, transcoefsidx, ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, transcoefs, ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, offsets, nterms);
   SCIPfreeBlockMemoryArray(scip, vars, nvars);
}
