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

/**@file   nlhdlr_quadratic.c
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  nonlinear handler to handle quadratic expressions
 * @author Felipe Serrano
 * @author Antonia Chmiela
 *
 * Some definitions:
 * - a `BILINEXPRTERM` is a product of two expressions
 * - a `QUADEXPRTERM` stores an expression `expr` that is known to appear in a nonlinear, quadratic term, that is
 *   `expr^2` or `expr*other_expr`. It stores its `sqrcoef` (that can be 0), its linear coef and all the bilinear expression
 *   terms in which expr appears.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define DEBUG_INTERSECTIONCUT */
/* #define INTERCUT_MOREDEBUG */
/* #define INTERCUTS_VERBOSE */

#ifdef INTERCUTS_VERBOSE
#define INTER_LOG
#endif

#ifdef INTER_LOG
#define INTERLOG(x) if( SCIPgetSubscipDepth(scip) == 0 && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL ) { x }
#else
#define INTERLOG(x)
#endif

#include <string.h>

#include "scip/cons_nonlinear.h"
#include "scip/pub_nlhdlr.h"
#include "scip/nlhdlr_quadratic.h"
#include "scip/expr_pow.h"
#include "scip/expr_sum.h"
#include "scip/expr_var.h"
#include "scip/expr_product.h"
#include "scip/pub_misc_rowprep.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME                    "quadratic"
#define NLHDLR_DESC                    "handler for quadratic expressions"
#define NLHDLR_DETECTPRIORITY          1
#define NLHDLR_ENFOPRIORITY            100

/* properties of the quadratic nlhdlr statistics table */
#define TABLE_NAME_QUADRATIC           "nlhdlr_quadratic"
#define TABLE_DESC_QUADRATIC           "quadratic nlhdlr statistics table"
#define TABLE_POSITION_QUADRATIC       14700                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_QUADRATIC SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

/* some default values */
#define INTERCUTS_MINVIOL              1e-4
#define DEFAULT_USEINTERCUTS           FALSE
#define DEFAULT_USESTRENGTH            FALSE
#define DEFAULT_USEBOUNDS              FALSE
#define BINSEARCH_MAXITERS             120
#define DEFAULT_NCUTSROOT              20
#define DEFAULT_NCUTS                  2

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_NlhdlrExprData
{
   SCIP_EXPR*            qexpr;              /**< quadratic expression (stored here again for convenient access) */

   SCIP_EXPRCURV         curvature;          /**< curvature of the quadratic representation of the expression */

   SCIP_INTERVAL         linactivity;        /**< activity of linear part */

   /* activities of quadratic parts as defined in nlhdlrIntevalQuadratic */
   SCIP_Real             minquadfiniteact;   /**< minimum activity of quadratic part where only terms with finite min
                                                  activity contribute */
   SCIP_Real             maxquadfiniteact;   /**< maximum activity of quadratic part where only terms with finite max
                                                  activity contribute */
   int                   nneginfinityquadact;/**< number of quadratic terms contributing -infinity to activity */
   int                   nposinfinityquadact;/**< number of quadratic terms contributing +infinity to activity */
   SCIP_INTERVAL*        quadactivities;     /**< activity of each quadratic term as defined in nlhdlrIntevalQuadratic */
   SCIP_INTERVAL         quadactivity;       /**< activity of quadratic part (sum of quadactivities) */
   SCIP_Longint          activitiestag;      /**< value of activities tag when activities were computed */

   SCIP_CONS*            cons;               /**< if expr is the root of constraint cons, store cons; otherwise NULL */
   SCIP_Bool             separating;         /**< whether we are using the nlhdlr also for separation */
   SCIP_Bool             origvars;           /**< whether the quad expr in qexpr is in original (non-aux) variables */

   int                   ncutsadded;         /**< number of intersection cuts added for this quadratic */
};

/** nonlinear handler data */
struct SCIP_NlhdlrData
{
   int                   ncutsgenerated;     /**< total number of cuts that where generated by separateQuadratic */
   int                   ncutsadded;         /**< total number of cuts that where generated by separateQuadratic and actually added */
   SCIP_Longint          lastnodenumber;     /**< number of last node for which cuts were (allowed to be) generated */
   int                   lastncuts;          /**< number of cuts already generated */

   /* parameter */
   SCIP_Bool             useintersectioncuts; /**< whether to use intersection cuts for quadratic constraints or not */
   SCIP_Bool             usestrengthening;   /**< whether the strengthening should be used */
   SCIP_Bool             useboundsasrays;    /**< use bounds of variables in quadratic as rays for intersection cuts */
   int                   ncutslimit;         /**< limit for number of cuts generated consecutively */
   int                   ncutslimitroot;     /**< limit for number of cuts generated at root node */
   int                   maxrank;            /**< maximal rank a slackvar can have */
   SCIP_Real             mincutviolation;    /**< minimal cut violation the generated cuts must fulfill to be added to the LP */
   SCIP_Real             minviolation;       /**< minimal violation the constraint must fulfill such that a cut can be generated */
   int                   atwhichnodes;       /**< determines at which nodes cut is used (if it's -1, it's used only at the root node,
                                                  if it's n >= 0, it's used at every multiple of n) */
   int                   nstrengthlimit;     /**< limit for number of rays we do the strengthening for */
   SCIP_Real             cutcoefsum;         /**< sum of average cutcoefs of a cut */
   SCIP_Bool             ignorebadrayrestriction; /**< should cut be generated even with bad numerics when restricting to ray? */
   SCIP_Bool             ignorehighre;       /**< should cut be added even when range / efficacy is large? */

   /* statistics */
   int                   ncouldimprovedcoef; /**< number of times a coefficient could improve but didn't because of numerics */
   int                   nbadrayrestriction; /**< number of times a cut was aborted because of numerics when restricting to ray */
   int                   nbadnonbasic;       /**< number of times a cut was aborted because the nonbasic row was not nonbasic enough */
   int                   nhighre;            /**< number of times a cut was not added because range / efficacy was too large */
   int                   nphinonneg;         /**< number of times a cut was aborted because phi is nonnegative at 0 */
   int                   nstrengthenings;    /**< number of successful strengthenings */
   int                   nboundcuts;         /**< number of successful bound cuts */
   SCIP_Real             ncalls;             /**< number of calls to separation */
   SCIP_Real             densitysum;         /**< sum of density of cuts */
};

/* structure to store rays. note that for a given ray, the entries in raysidx are sorted. */
struct Rays
{
   SCIP_Real*            rays;               /**< coefficients of rays */
   int*                  raysidx;            /**< to which var the coef belongs; vars are [quadvars, linvars, auxvar] */
   int*                  raysbegin;          /**< positions of rays: coefs of i-th ray [raybegin[i], raybegin[i+1]) */
   int*                  lpposray;           /**< lp pos of var associated with the ray;
                                               >= 0 -> lppos of var; < 0 -> var is row -lpposray -1 */

   int                   rayssize;           /**< size of rays and rays idx */
   int                   nrays;              /**< size of raysbegin is nrays + 1; size of lpposray */
};
typedef struct Rays RAYS;


/*
 * Callback methods of the table
 */

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputQuadratic)
{ /*lint --e{715}*/
   SCIP_NLHDLR* nlhdlr;
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSHDLR* conshdlr;

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   assert(conshdlr != NULL);
   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, NLHDLR_NAME);
   assert(nlhdlr != NULL);
   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata);

   /* print statistics */
   SCIPinfoMessage(scip, file, "Quadratic Nlhdlr   : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n", "GenCuts", "AddCuts", "CouldImpr", "NLargeRE",
         "AbrtBadRay", "AbrtPosPhi", "AbrtNonBas", "NStrength", "AveCutcoef", "AveDensity", "AveBCutsFrac");
   SCIPinfoMessage(scip, file, "  %-17s:", "Quadratic Nlhdlr");
   SCIPinfoMessage(scip, file, " %10d", nlhdlrdata->ncutsgenerated);
   SCIPinfoMessage(scip, file, " %10d", nlhdlrdata->ncutsadded);
   SCIPinfoMessage(scip, file, " %10d", nlhdlrdata->ncouldimprovedcoef);
   SCIPinfoMessage(scip, file, " %10d", nlhdlrdata->nhighre);
   SCIPinfoMessage(scip, file, " %10d", nlhdlrdata->nbadrayrestriction);
   SCIPinfoMessage(scip, file, " %10d", nlhdlrdata->nphinonneg);
   SCIPinfoMessage(scip, file, " %10d", nlhdlrdata->nbadnonbasic);
   SCIPinfoMessage(scip, file, " %10d", nlhdlrdata->nstrengthenings);
   SCIPinfoMessage(scip, file, " %10g", nlhdlrdata->nstrengthenings > 0 ? nlhdlrdata->cutcoefsum / nlhdlrdata->nstrengthenings : 0.0);
   SCIPinfoMessage(scip, file, " %10g", nlhdlrdata->ncutsadded > 0 ? nlhdlrdata->densitysum / nlhdlrdata->ncutsadded : 0.0);
   SCIPinfoMessage(scip, file, " %10g", nlhdlrdata->ncalls > 0 ? nlhdlrdata->nboundcuts / nlhdlrdata->ncalls : 0.0);
   SCIPinfoMessage(scip, file, "\n");

   return SCIP_OKAY;
}


/*
 * static methods
 */

/** adds cutcoef * (col - col*) to rowprep */
static
SCIP_RETCODE addColToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store intersection cut */
   SCIP_SOL*             sol,                /**< solution to separate */
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
   SCIProwprepAddConstant(rowprep, -cutcoef * SCIPgetSolVal(scip, sol, SCIPcolGetVar(col)) );

   return SCIP_OKAY;
}

/** adds cutcoef * (slack - slack*) to rowprep
 *
  * row is lhs &le; <coefs, vars> + constant &le; rhs, thus slack is defined by
  * slack + <coefs.vars> + constant = side
  *
  * If row (slack) is at upper, it means that <coefs,vars*> + constant = rhs, and so
  * slack* = side - rhs --> slack - slack* = rhs - <coefs, vars> - constant.
  *
  * If row (slack) is at lower, then <coefs,vars*> + constant = lhs, and so
  * slack* = side - lhs --> slack - slack* = lhs - <coefs, vars> - constant.
  */
static
SCIP_RETCODE addRowToCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store intersection cut */
   SCIP_Real             cutcoef,            /**< cut coefficient */
   SCIP_ROW*             row,                /**< row, whose slack we are ading to rowprep */
   SCIP_Bool*            success             /**< if the row is nonbasic enough */
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

/** counts the number of basic variables in the quadratic expr */
static
int countBasicVars(
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   SCIP_Bool*            nozerostat          /**< whether there is no variable with basis status zero */
   )
{
   SCIP_EXPR* qexpr;
   SCIP_EXPR** linexprs;
   SCIP_COL* col;
   int i;
   int nbasicvars = 0;
   int nquadexprs;
   int nlinexprs;

   *nozerostat = TRUE;

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, NULL, NULL, NULL);

   /* loop over quadratic vars */
   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_EXPR* expr;

      SCIPexprGetQuadraticQuadTerm(qexpr, i, &expr, NULL, NULL, NULL, NULL, NULL);

      col = SCIPvarGetCol(SCIPgetExprAuxVarNonlinear(expr));
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
         nbasicvars += 1;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
      {
         *nozerostat = FALSE;
         return 0;
      }
   }

   /* loop over linear vars */
   for( i = 0; i < nlinexprs; ++i )
   {
      col = SCIPvarGetCol(SCIPgetExprAuxVarNonlinear(linexprs[i]));
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
         nbasicvars += 1;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
      {
         *nozerostat = FALSE;
         return 0;
      }
   }

   /* finally consider the aux var (if it exists) */
   if( auxvar != NULL )
   {
      col = SCIPvarGetCol(auxvar);
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
         nbasicvars += 1;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
      {
         *nozerostat = FALSE;
         return 0;
      }
   }

   return nbasicvars;
}

/** stores the row of the tableau where `col` is basic
 *
 *  In general, we will have
 *
 *      basicvar1 = tableaurow var1
 *      basicvar2 = tableaurow var2
 *      ...
 *      basicvarn = tableaurow varn
 *
 *  However, we want to store the the tableau row by columns. Thus, we need to know which of the basic vars `col` is.
 *
 *  Note we only store the entries of the nonbasic variables
 */
static
SCIP_RETCODE storeDenseTableauRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col,                /**< basic columns to store its tableau row */
   int*                  basicvarpos2tableaurow,/**< map between basic var and its tableau row */
   int                   nbasiccol,          /**< which basic var this is */
   int                   raylength,          /**< the length of a ray (the total number of basic vars) */
   SCIP_Real*            binvrow,            /**< buffer to store row of Binv */
   SCIP_Real*            binvarow,           /**< buffer to store row of Binv A */
   SCIP_Real*            tableaurows         /**< pointer to store the tableau rows */
   )
{
   int nrows;
   int ncols;
   int lppos;
   int i;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int nray;

   assert(nbasiccol < raylength);
   assert(col != NULL);
   assert(binvrow != NULL);
   assert(binvarow != NULL);
   assert(tableaurows != NULL);
   assert(basicvarpos2tableaurow != NULL);
   assert(SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC);

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   lppos = SCIPcolGetLPPos(col);

   assert(basicvarpos2tableaurow[lppos] >= 0);

   SCIP_CALL( SCIPgetLPBInvRow(scip, basicvarpos2tableaurow[lppos], binvrow, NULL, NULL) );
   SCIP_CALL( SCIPgetLPBInvARow(scip, basicvarpos2tableaurow[lppos], binvrow, binvarow, NULL, NULL) );

   nray = 0;
   for( i = 0; i < ncols; ++i )
      if( SCIPcolGetBasisStatus(cols[i]) != SCIP_BASESTAT_BASIC )
      {
         tableaurows[nbasiccol + nray * raylength] = binvarow[i];
         nray++;
      }
   for( ; i < ncols + nrows; ++i )
      if( SCIProwGetBasisStatus(rows[i - ncols]) != SCIP_BASESTAT_BASIC )
      {
         tableaurows[nbasiccol + nray * raylength] = binvrow[i - ncols];
         nray++;
      }

   return SCIP_OKAY;
}

/** stores the rows of the tableau corresponding to the basic variables in the quadratic expression
 *
 * Also return a map storing to which var the entry of a ray corresponds, i.e., if the tableau is
 *
 *     basicvar_1 = ray1_1 nonbasicvar_1 + ...
 *     basicvar_2 = ray1_2 nonbasicvar_1 + ...
 *     ...
 *     basicvar_n = ray1_n nonbasicvar_1 + ...
 *
 * The map maps k to the position of basicvar_k in the variables of the constraint assuming the variables are sorted as
 * [quadratic vars, linear vars, auxvar].
 */
static
SCIP_RETCODE storeDenseTableauRowsByColumns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   int                   raylength,          /**< length of a ray of the tableau */
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   SCIP_Real*            tableaurows,        /**< buffer to store the tableau rows */
   int*                  rayentry2conspos    /**< buffer to store the map */
   )
{
   SCIP_EXPR* qexpr;
   SCIP_EXPR** linexprs;
   SCIP_Real* binvarow;
   SCIP_Real* binvrow;
   SCIP_COL* col;
   int* basicvarpos2tableaurow; /* map between basic var and its tableau row */
   int nrayentries;
   int nquadexprs;
   int nlinexprs;
   int nrows;
   int ncols;
   int i;

   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &basicvarpos2tableaurow, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvarow, ncols) );

   for( i = 0; i < ncols; ++i )
      basicvarpos2tableaurow[i] = -1;
   SCIP_CALL( constructBasicVars2TableauRowMap(scip, basicvarpos2tableaurow) );

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, NULL, NULL, NULL);

   /* entries of quadratic basic vars */
   nrayentries = 0;
   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_EXPR* expr;
      SCIPexprGetQuadraticQuadTerm(qexpr, i, &expr, NULL, NULL, NULL, NULL, NULL);

      col = SCIPvarGetCol(SCIPgetExprAuxVarNonlinear(expr));
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
      {
         SCIP_CALL( storeDenseTableauRow(scip, col, basicvarpos2tableaurow, nrayentries, raylength, binvrow, binvarow,
                  tableaurows) );

         rayentry2conspos[nrayentries] = i;
         nrayentries++;
      }
   }
   /* entries of linear vars */
   for( i = 0; i < nlinexprs; ++i )
   {
      col = SCIPvarGetCol(SCIPgetExprAuxVarNonlinear(linexprs[i]));
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
      {
         SCIP_CALL( storeDenseTableauRow(scip, col, basicvarpos2tableaurow, nrayentries, raylength, binvrow, binvarow,
                  tableaurows) );

         rayentry2conspos[nrayentries] = nquadexprs + i;
         nrayentries++;
      }
   }
   /* entry of aux var (if it exists) */
   if( auxvar != NULL )
   {
      col = SCIPvarGetCol(auxvar);
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC )
      {
         SCIP_CALL( storeDenseTableauRow(scip, col, basicvarpos2tableaurow, nrayentries, raylength, binvrow, binvarow,
                  tableaurows) );

         rayentry2conspos[nrayentries] = nquadexprs + nlinexprs;
         nrayentries++;
      }
   }
   assert(nrayentries == raylength);

#ifdef  DEBUG_INTERSECTIONCUT
   for( i = 0; i < ncols; ++i )
   {
      SCIPinfoMessage(scip, NULL, "%d column of tableau is: ",i);
      for( int j = 0; j < raylength; ++j )
         SCIPinfoMessage(scip, NULL, "%g ", tableaurows[i * raylength + j]);
      SCIPinfoMessage(scip, NULL, "\n");
   }
#endif

   SCIPfreeBufferArray(scip, &binvarow);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basicvarpos2tableaurow);

   return SCIP_OKAY;
}

/** initializes rays data structure */
static
SCIP_RETCODE createRays(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS**                rays                /**< rays data structure */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, rays) );
   BMSclearMemory(*rays);

   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->rays, SCIPgetNLPCols(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->raysidx, SCIPgetNLPCols(scip)) );

   /* overestimate raysbegin and lpposray */
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->raysbegin, SCIPgetNLPCols(scip) + SCIPgetNLPRows(scip) + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->lpposray, SCIPgetNLPCols(scip) + SCIPgetNLPRows(scip)) );
   (*rays)->raysbegin[0] = 0;

   (*rays)->rayssize = SCIPgetNLPCols(scip);

   return SCIP_OKAY;
}

/** initializes rays data structure for bound rays */
static
SCIP_RETCODE createBoundRays(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS**                rays,               /**< rays data structure */
   int                   size                /**< number of rays to allocate */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, rays) );
   BMSclearMemory(*rays);

   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->rays, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->raysidx, size) );

   /* overestimate raysbegin and lpposray */
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->raysbegin, size + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*rays)->lpposray, size) );
   (*rays)->raysbegin[0] = 0;

   (*rays)->rayssize = size;

   return SCIP_OKAY;
}

/** frees rays data structure */
static
void freeRays(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS**                rays                /**< rays data structure */
   )
{
   if( *rays == NULL )
      return;

   SCIPfreeBufferArray(scip, &(*rays)->lpposray);
   SCIPfreeBufferArray(scip, &(*rays)->raysbegin);
   SCIPfreeBufferArray(scip, &(*rays)->raysidx);
   SCIPfreeBufferArray(scip, &(*rays)->rays);

   SCIPfreeBuffer(scip, rays);
}

/** inserts entry to rays data structure; checks and resizes if more space is needed */
static
SCIP_RETCODE insertRayEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS*                 rays,               /**< rays data structure */
   SCIP_Real             coef,               /**< coefficient to insert */
   int                   coefidx,            /**< index of coefficient (conspos of var it corresponds to) */
   int                   coefpos             /**< where to insert coefficient */
   )
{
   /* check for size */
   if( rays->rayssize <= coefpos + 1 )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, coefpos + 1);
      SCIP_CALL( SCIPreallocBufferArray(scip, &(rays->rays), newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &(rays->raysidx), newsize) );
      rays->rayssize = newsize;
   }

   /* insert entry */
   rays->rays[coefpos] = coef;
   rays->raysidx[coefpos] = coefidx;

   return SCIP_OKAY;
}

/** constructs map between the lppos of a variables and its position in the constraint assuming the constraint variables
 * are sorted as [quad vars, lin vars, aux var (if it exists)]
 *
 * If a variable doesn't appear in the constraint, then its position is -1.
 */
static
void constructLPPos2ConsPosMap(
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_VAR*             auxvar,             /**< aux var of the expr */
   int*                  map                 /**< buffer to store the mapping */
   )
{
   SCIP_EXPR* qexpr;
   SCIP_EXPR** linexprs;
   int nquadexprs;
   int nlinexprs;
   int lppos;
   int i;

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, NULL, NULL, NULL);

   /* set pos of quadratic vars */
   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_EXPR* expr;
      SCIPexprGetQuadraticQuadTerm(qexpr, i, &expr, NULL, NULL, NULL, NULL, NULL);

      lppos = SCIPcolGetLPPos(SCIPvarGetCol(SCIPgetExprAuxVarNonlinear(expr)));
      map[lppos] = i;
   }
   /* set pos of lin vars */
   for( i = 0; i < nlinexprs; ++i )
   {
      lppos = SCIPcolGetLPPos(SCIPvarGetCol(SCIPgetExprAuxVarNonlinear(linexprs[i])));
      map[lppos] = nquadexprs + i;
   }
   /* set pos of aux var (if it exists) */
   if( auxvar != NULL )
   {
      lppos = SCIPcolGetLPPos(SCIPvarGetCol(auxvar));
      map[lppos] = nquadexprs + nlinexprs;
   }

   return;
}

/** inserts entries of factor * nray-th column of densetableaucols into rays data structure */
static
SCIP_RETCODE insertRayEntries(
   SCIP*                 scip,               /**< SCIP data structure */
   RAYS*                 rays,               /**< rays data structure */
   SCIP_Real*            densetableaucols,   /**< column of the tableau in dense format */
   int*                  rayentry2conspos,   /**< map between rayentry and conspos of associated var */
   int                   raylength,          /**< length of a tableau column */
   int                   nray,               /**< which tableau column to insert */
   int                   conspos,            /**< conspos of ray's nonbasic var in the cons; -1 if not in the cons */
   SCIP_Real             factor,             /**< factor to multiply each tableau col */
   int*                  nnonz,              /**< position to start adding the ray in rays and buffer to store nnonz */
   SCIP_Bool*            success             /**< we can't separate if there is a nonzero ray with basis status ZERO */
   )
{
   int i;

   *success = TRUE;

   for( i = 0; i < raylength; ++i )
   {
      SCIP_Real coef;

      /* we have a nonzero ray with base stat zero -> can't generate cut */
      if( factor == 0.0 && ! SCIPisZero(scip, densetableaucols[nray * raylength + i]) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      coef = factor * densetableaucols[nray * raylength + i];

      /* this might be a source of numerical issues
       * TODO: in case of problems, an idea would be to scale the ray entries; compute the cut coef and scale it back down
       * another idea would be to check against a smaller epsilion.
       * The problem is that if the cut coefficient is 1/t where lpsol + t*ray intersects the S-free set.
       * Now if t is super big, then a super small coefficient would have had an impact...
       */
      if( SCIPisZero(scip, coef) )
         continue;

      /* check if nonbasic var entry should come before this one */
      if( conspos > -1 && conspos < rayentry2conspos[i] )
      {
         /* add nonbasic entry */
         assert(factor != 0.0);

#ifdef  DEBUG_INTERSECTIONCUT
         SCIPinfoMessage(scip, NULL, "ray belongs to nonbasic var pos %d in cons\n", conspos);
#endif

         SCIP_CALL( insertRayEntry(scip, rays, -factor, conspos, *nnonz) );
         (*nnonz)++;

         /* we are done with nonbasic entry */
         conspos = -1;
      }

      SCIP_CALL( insertRayEntry(scip, rays, coef, rayentry2conspos[i], *nnonz) );
      (*nnonz)++;
   }

   /* if nonbasic entry was not added and should still be added, then it should go at the end */
   if( conspos > -1 )
   {
      /* add nonbasic entry */
      assert(factor != 0.0);

#ifdef  DEBUG_INTERSECTIONCUT
      SCIPinfoMessage(scip, NULL, "ray belongs to nonbasic var pos %d in cons\n", conspos);
#endif

      SCIP_CALL( insertRayEntry(scip, rays, -factor, conspos, *nnonz) );
      (*nnonz)++;
   }

   /* finished ray entry; store its end */
   rays->raysbegin[rays->nrays + 1] = *nnonz;

#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "entries of ray %d are between [%d, %d):\n", rays->nrays, rays->raysbegin[rays->nrays], *nnonz);
   for( i = rays->raysbegin[rays->nrays]; i < *nnonz; ++i )
      SCIPinfoMessage(scip, NULL, "(%d, %g), ", rays->raysidx[i], rays->rays[i]);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   return SCIP_OKAY;
}

/** stores rays in sparse form
 *
 * The first rays correspond to the nonbasic variables
 * and the last rays to the nonbasic slack variables.
 *
 * More details: The LP tableau is of the form
 *
 *     basicvar_1 = ray1_1 nonbasicvar_1 + ... + raym_1 nonbasicvar_m
 *     basicvar_2 = ray1_2 nonbasicvar_1 + ... + raym_2 nonbasicvar_m
 *     ...
 *     basicvar_n = ray1_n nonbasicvar_1 + ... + raym_n nonbasicvar_m
 *     nonbasicvar_1 = 1.0 nonbasicvar_1 + ... +    0.0 nonbasicvar_m
 *     ...
 *     nonbasicvar_m = 0.0 nonbasicvar_1 + ... +    1.0 nonbasicvar_m
 *
 *  so rayk = (rayk_1, ... rayk_n, e_k)
 *  We store the entries of the rays associated to the variables present in the quadratic expr.
 *  We do not store zero rays.
 *
 *  Also, we store the rays as if every nonbasic variable was at lower (so that all rays moves to infinity)
 *  Since the tableau is:
 *
 *      basicvar + Binv L (nonbasic_lower - lb) + Binv U (nonbasic_upper - ub) = basicvar_sol
 *
 *  then:
 *
 *      basicvar = basicvar_sol - Binv L (nonbasic_lower - lb) + Binv U (ub - nonbasic_upper)
 *
 *  and so the entries of the rays associated with the basic variables are:
 *  rays_basicvars = [-BinvL, BinvU].
 *
 *  So we flip the sign of the rays associated to nonbasic vars at lower.
 *  In constrast, the nonbasic part of the ray has a 1.0 for nonbasic at lower and a -1.0 for nonbasic at upper, i.e.
 *  nonbasic_lower = lb + 1.0(nonbasic_lower - lb) and
 *  nonbasic_upper = ub - 1.0(ub - nonbasic_upper)
 */
static
SCIP_RETCODE createAndStoreSparseRays(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   RAYS**                raysptr,            /**< buffer to store rays datastructure */
   SCIP_Bool*            success             /**< we can't separate if there is a var with basis status ZERO */
   )
{
   SCIP_Real* densetableaucols;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   RAYS* rays;
   int* rayentry2conspos;
   int* lppos2conspos;
   int nnonbasic;
   int nrows;
   int ncols;
   int nnonz;
   int raylength;
   int i;

   nrows = SCIPgetNLPRows(scip);
   ncols = SCIPgetNLPCols(scip);

   *success = TRUE;

   raylength = countBasicVars(nlhdlrexprdata, auxvar, success);
   if( ! *success )
   {
      SCIPdebugMsg(scip, "failed to store sparse rays: there is a var with base status zero\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &densetableaucols, raylength * (ncols + nrows)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rayentry2conspos, raylength) );

   /* construct dense tableau and map between ray entries and position of corresponding var in quad cons */
   SCIP_CALL( storeDenseTableauRowsByColumns(scip, nlhdlrexprdata, raylength, auxvar,
            densetableaucols, rayentry2conspos) );

   /* build rays sparsely now */
   SCIP_CALL( SCIPallocBufferArray(scip, &lppos2conspos, ncols) );
   for( i = 0; i < ncols; ++i )
      lppos2conspos[i] = -1;

   constructLPPos2ConsPosMap(nlhdlrexprdata, auxvar, lppos2conspos);

   /* store sparse rays */
   SCIP_CALL( createRays(scip, raysptr) );
   rays = *raysptr;

   nnonz = 0;
   nnonbasic = 0;

   /* go through nonbasic variables */
   cols = SCIPgetLPCols(scip);
   for( i = 0; i < ncols; ++i )
   {
      int oldnnonz = nnonz;
      SCIP_COL* col;
      SCIP_Real factor;

      col = cols[i];

      /* set factor to store basic entries of ray as = [-BinvL, BinvU] */
      if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER )
         factor = -1.0;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER )
         factor = 1.0;
      else if( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_ZERO )
         factor = 0.0;
      else
         continue;

      SCIP_CALL( insertRayEntries(scip, rays, densetableaucols, rayentry2conspos, raylength, nnonbasic,
               lppos2conspos[SCIPcolGetLPPos(col)], factor, &nnonz, success) );

#ifdef  DEBUG_INTERSECTIONCUT
      SCIPinfoMessage(scip, NULL, "looked at ray of var %s with basestat %d, it has %d nonzeros\n-----------------\n",
            SCIPvarGetName(SCIPcolGetVar(col)), SCIPcolGetBasisStatus(col), nnonz - oldnnonz);
#endif
      if( ! (*success) )
      {
#ifdef  DEBUG_INTERSECTIONCUT
         SCIPdebugMsg(scip, "nonzero ray associated with variable <%s> has base status zero -> abort storing rays\n",
               SCIPvarGetName(SCIPcolGetVar(col)));
#endif
         goto CLEANUP;
      }

      /* if ray is non zero remember who it belongs to */
      assert(oldnnonz <= nnonz);
      if( oldnnonz < nnonz )
      {
         rays->lpposray[rays->nrays] = SCIPcolGetLPPos(col);
         (rays->nrays)++;
      }
      nnonbasic++;
   }

   /* go through nonbasic slack variables */
   rows = SCIPgetLPRows(scip);
   for( i = 0; i < nrows; ++i )
   {
      int oldnnonz = nnonz;
      SCIP_ROW* row;
      SCIP_Real factor;

      row = rows[i];

      /* set factor to store basic entries of ray as = [-BinvL, BinvU]; basic status of rows are flipped! See lpi.h! */
      assert(SCIProwGetBasisStatus(row) != SCIP_BASESTAT_ZERO);
      if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER )
         factor = 1.0;
      else if( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER )
         factor =-1.0;
      else
         continue;

      SCIP_CALL( insertRayEntries(scip, rays, densetableaucols, rayentry2conspos, raylength, nnonbasic, -1, factor,
               &nnonz, success) );
      assert(*success);

      /* if ray is non zero remember who it belongs to */
#ifdef  DEBUG_INTERSECTIONCUT
      SCIPinfoMessage(scip, NULL, "looked at ray of row %d, it has %d nonzeros\n-----------------\n", i, nnonz - oldnnonz);
#endif
      assert(oldnnonz <= nnonz);
      if( oldnnonz < nnonz )
      {
         rays->lpposray[rays->nrays] = -SCIProwGetLPPos(row) - 1;
         (rays->nrays)++;
      }
      nnonbasic++;
   }

CLEANUP:
   SCIPfreeBufferArray(scip, &lppos2conspos);
   SCIPfreeBufferArray(scip, &rayentry2conspos);
   SCIPfreeBufferArray(scip, &densetableaucols);

   if( ! *success )
   {
      freeRays(scip, &rays);
   }

   return SCIP_OKAY;
}

/* TODO: which function this comment belongs to? */
/* this function determines how the maximal S-free set is going to look like
 *
 * There are 4 possibilities: after writing the quadratic constraint
 * \f$q(z) \leq 0\f$
 * as
 * \f$\Vert x(z)\Vert^2 - \Vert y\Vert^2 + w(z) + kappa \leq 0\f$,
 * the cases are determined according to the following:
 * - Case 1: w = 0 and kappa = 0
 * - Case 2: w = 0 and kappa > 0
 * - Case 3: w = 0 and kappa < 0
 * - Case 4: w != 0
 */

/** compute quantities for intersection cuts
 *
 * Assume the quadratic is stored as
 * \f[ q(z) = z_q^T Q z_q + b_q^T z_q + b_l z_l + c - z_a \f]
 * where:
 *  - \f$z_q\f$ are the quadratic vars
 *  - \f$z_l\f$ are the linear vars
 *  - \f$z_a\f$ is the aux var if it exists
 *
 * We can rewrite it as
 * \f[ \Vert x(z)\Vert^2 - \Vert y\Vert^2 + w(z) + \kappa \leq 0. \f]
 * To do this transformation and later to compute the actual cut we need to compute and store some quantities.
 * Let
 *    - \f$I_0\f$, \f$I_+\f$, and \f$I_-\f$ be the index set of zero, positive, and negative eigenvalues, respectively
 *    - \f$v_i\f$ be the i-th eigenvector of \f$Q\f$
 *    - \f$zlp\f$ be the lp value of the variables \f$z\f$
 *
 * The quantities we need are:
 *    - \f$vb_i = v_i^T b\f$ for \f$i \in I_+ \cup I_-\f$
 *    - \f$vzlp_i = v_i^T zlp_q\f$ for \f$i \in I_+ \cup I_-\f$
 *    - \f$\kappa = c - 1/4 \sum_{i \in I_+ \cup I_-} (v_i^T b_q)^2 / eigval_i\f$
 *    - \f$w(z) = (\sum_{i \in I_0} v_i^T b_q v_i^T) z_q + b_l^T z_l - z_a\f$
 *    - \f$w(zlp)\f$
 *
 * @return \f$\kappa\f$ and the vector \f$\sum_{i \in I_0} v_i^T b_q v_i^T\f$
 *
 * @note if the constraint is q(z) &le; rhs, then the constant when writing the constraint as quad &le; 0 is c - rhs.
 * @note if the quadratic constraint we are separating is q(z) &ge; lhs, then we multiply by -1.
 * In practice, what changes is
 *    - the sign of the eigenvalues
 *    - the sign of \f$b_q\f$ and \f$b_l\f$
 *    - the sign of the coefficient of the auxvar (if it exists)
 *    - the constant of the quadratic written as quad &le; 0 is lhs - c
 * @note The eigenvectors _do not_ change sign!
 */
static
SCIP_RETCODE intercutsComputeCommonQuantities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   SCIP_Real             sidefactor,         /**< 1.0 if the violated constraint is q &le; rhs, -1.0 otherwise */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_Real*            vb,                 /**< buffer to store \f$v_i^T b\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            vzlp,               /**< buffer to store \f$v_i^T zlp_q\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            wcoefs,             /**< buffer to store the coefs of quad vars of w */
   SCIP_Real*            wzlp,               /**< pointer to store the value of w at zlp */
   SCIP_Real*            kappa               /**< pointer to store the value of kappa */
   )
{
   SCIP_EXPR* qexpr;
   SCIP_EXPR** linexprs;
   SCIP_Real* eigenvectors;
   SCIP_Real* eigenvalues;
   SCIP_Real* lincoefs;
   SCIP_Real constant; /* constant of the quadratic when written as <= 0 */
   int nquadexprs;
   int nlinexprs;
   int i;
   int j;

   assert(sidefactor == 1.0 || sidefactor == -1.0);
   assert(nlhdlrexprdata->cons != NULL || auxvar != NULL );

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, NULL, &eigenvalues,
         &eigenvectors);

   assert( eigenvalues != NULL );

   /* first get constant of quadratic when written as quad <= 0 */
   if( nlhdlrexprdata->cons != NULL )
      constant = (sidefactor == 1.0) ? constant - SCIPgetRhsNonlinear(nlhdlrexprdata->cons) :
         SCIPgetLhsNonlinear(nlhdlrexprdata->cons) - constant;
   else
      constant = (sidefactor * constant);

   *kappa = 0.0;
   *wzlp = 0.0;
   BMSclearMemoryArray(wcoefs, nquadexprs);

   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_Real vdotb;
      SCIP_Real vdotzlp;
      int offset;

      offset = i * nquadexprs;

      /* compute v_i^T b and v_i^T zlp */
      vdotb = 0;
      vdotzlp = 0;
      for( j = 0; j < nquadexprs; ++j )
      {
         SCIP_EXPR* expr;
         SCIP_Real lincoef;

         SCIPexprGetQuadraticQuadTerm(qexpr, j, &expr, &lincoef, NULL, NULL, NULL, NULL);

         vdotb += (sidefactor * lincoef) * eigenvectors[offset + j];
#ifdef INTERCUT_MOREDEBUG
         printf("vdotb: offset %d, eigenvector %d = %g, lincoef quad %g\n", offset, j,
               eigenvectors[offset + j], lincoef);
#endif
         vdotzlp += SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr)) * eigenvectors[offset + j];
      }
      vb[i] = vdotb;
      vzlp[i] = vdotzlp;

      if( ! SCIPisZero(scip, eigenvalues[i]) )
      {
         /* nonzero eigenvalue: compute kappa */
         *kappa += SQR(vdotb) / (sidefactor * eigenvalues[i]);
      }
      else
      {
         /* compute coefficients of w and compute w at zlp */
         for( j = 0; j < nquadexprs; ++j )
            wcoefs[j] += vdotb * eigenvectors[offset + j];

         *wzlp += vdotb * vdotzlp;
      }
   }

   /* finish kappa computation */
   *kappa *= -0.25;
   *kappa += constant;

   /* finish w(zlp) computation: linear part (including auxvar, if applicable) */
   for( i = 0; i < nlinexprs; ++i )
   {
      *wzlp += (sidefactor * lincoefs[i]) * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(linexprs[i]));
   }
   if( auxvar != NULL )
   {
      *wzlp += (sidefactor * -1.0) * SCIPgetSolVal(scip, sol, auxvar);
   }

#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "Computed common quantities needed for intercuts:\n");
   SCIPinfoMessage(scip, NULL, "   kappa = %g\n   quad part w = ", *kappa);
   for( i = 0; i < nquadexprs; ++i )
   {
      SCIPinfoMessage(scip, NULL, "%g ", wcoefs[i]);
   }
   SCIPinfoMessage(scip, NULL, "\n");
#endif
   return SCIP_OKAY;
}

/** computes eigenvec^T ray */
static
SCIP_Real computeEigenvecDotRay(
   SCIP_Real*            eigenvec,           /**< eigenvector */
   int                   nquadvars,          /**< number of quadratic vars (length of eigenvec) */
   SCIP_Real*            raycoefs,           /**< coefficients of ray */
   int*                  rayidx,             /**< index of consvar the ray coef is associated to */
   int                   raynnonz            /**< length of raycoefs and rayidx */
   )
{
   SCIP_Real retval;
   int i;

   retval = 0.0;
   for( i = 0; i < raynnonz; ++i )
   {
      /* rays are sorted; the first entries correspond to the quad vars, so we stop after first nonquad var entry */
      if( rayidx[i] >= nquadvars )
         break;

      retval += eigenvec[rayidx[i]] * raycoefs[i];
   }

   return retval;
}

/** computes linear part of evaluation of w(ray): \f$b_l^T ray_l - ray_a\f$
 *
 * @note we can know whether the auxiliary variable appears by the entries of the ray
 */
static
SCIP_Real computeWRayLinear(
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_Real             sidefactor,         /**< 1.0 if the violated constraint is q &le; rhs, -1.0 otherwise */
   SCIP_Real*            raycoefs,           /**< coefficients of ray */
   int*                  rayidx,             /**< ray coef[i] affects var at pos rayidx[i] in consvar */
   int                   raynnonz            /**< length of raycoefs and rayidx */
   )
{
   SCIP_EXPR* qexpr;
   SCIP_Real* lincoefs;
   SCIP_Real retval;
   int nquadexprs;
   int nlinexprs;

   int i;
   int start;

#ifdef INTERCUT_MOREDEBUG
   printf("Computing w(ray) \n");
#endif
   retval = 0.0;
   start = raynnonz - 1;

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, NULL, &nlinexprs, NULL, &lincoefs, &nquadexprs, NULL, NULL, NULL);

   /* process ray entry associated to the auxvar if it applies */
   if( rayidx[raynnonz - 1] == nquadexprs + nlinexprs )
   {
#ifdef INTERCUT_MOREDEBUG
      printf("wray auxvar term %g \n", (sidefactor * -1.0) * raycoefs[raynnonz - 1]);
#endif
      retval += (sidefactor * -1.0) * raycoefs[raynnonz - 1];
      start--;
   }

   /* process the rest of the entries */
   for( i = start; i >= 0; --i )
   {
      /* rays are sorted; last entries correspond to the lin vars, so we stop after first quad var entry */
      if( rayidx[i] < nquadexprs )
         break;

#ifdef INTERCUT_MOREDEBUG
      printf("wray var in pos %d term %g:: lincoef %g raycoef %g\n", rayidx[i], (sidefactor *
            lincoefs[rayidx[i] - nquadexprs]) * raycoefs[i], lincoefs[rayidx[i] - nquadexprs] ,raycoefs[i]);
#endif
      retval += (sidefactor * lincoefs[rayidx[i] - nquadexprs]) * raycoefs[i] ;
   }

   return retval;
}

/** calculate coefficients of restriction of the function to given ray.
 *
 * The restriction of the function representing the maximal S-free set to zlp + t * ray has the form
 * SQRT(A t^2 + B t + C) - (D t + E) for cases 1, 2, and 3.
 * For case 4 it is a piecewise defined function and each piece is of the aforementioned form.
 *
 * This function computes the coefficients A, B, C, D, E for the given ray.
 * In case 4, it computes the coefficients for both pieces, in addition to coefficients needed to evaluate the condition
 * in the piecewise definition of the function.
 *
 * The parameter iscase4 tells the function if it is case 4 or not.
 */
static
SCIP_RETCODE computeRestrictionToRay(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_Real             sidefactor,         /**< 1.0 if the violated constraint is q &le; rhs, -1.0 otherwise */
   SCIP_Bool             iscase4,            /**< whether we are in case 4 */
   SCIP_Real*            raycoefs,           /**< coefficients of ray */
   int*                  rayidx,             /**< index of consvar the ray coef is associated to */
   int                   raynnonz,           /**< length of raycoefs and rayidx */
   SCIP_Real*            vb,                 /**< array containing \f$v_i^T b\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            vzlp,               /**< array containing \f$v_i^T zlp_q\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            wcoefs,             /**< coefficients of w for the qud vars or NULL if w is 0 */
   SCIP_Real             wzlp,               /**< value of w at zlp */
   SCIP_Real             kappa,              /**< value of kappa */
   SCIP_Real*            coefs1234a,         /**< buffer to store A, B, C, D, and E of cases 1, 2, 3, or 4a */
   SCIP_Real*            coefs4b,            /**< buffer to store A, B, C, D, and E of case 4b (or NULL if not needed) */
   SCIP_Real*            coefscondition,     /**< buffer to store data to evaluate condition to decide case 4a or 4b */
   SCIP_Bool*            success             /**< did we successfully compute the coefficients? */
   )
{
   SCIP_EXPR* qexpr;
   int nquadexprs;
   SCIP_Real* eigenvectors;
   SCIP_Real* eigenvalues;
   SCIP_Real* a;
   SCIP_Real* b;
   SCIP_Real* c;
   SCIP_Real* d;
   SCIP_Real* e;
   SCIP_Real wray;
   int i;

   *success = TRUE;

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, NULL, NULL, NULL, NULL, &nquadexprs, NULL, &eigenvalues, &eigenvectors);

#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "\n############################################\n");
   SCIPinfoMessage(scip, NULL, "Restricting to ray:\n");
   for( i = 0; i < raynnonz; ++i )
   {
      SCIPinfoMessage(scip, NULL, "(%d, %g), ", rayidx[i], raycoefs[i]);
   }
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   assert(coefs1234a != NULL);

   /* set all coefficients to zero */
   memset(coefs1234a, 0, 5 * sizeof(SCIP_Real));
   if( iscase4 )
   {
      assert(coefs4b != NULL);
      assert(coefscondition != NULL);
      assert(wcoefs != NULL);

      memset(coefs4b, 0, 5 * sizeof(SCIP_Real));
      memset(coefscondition, 0, 3 * sizeof(SCIP_Real));
   }

   a = coefs1234a;
   b = coefs1234a + 1;
   c = coefs1234a + 2;
   d = coefs1234a + 3;
   e = coefs1234a + 4;
   wray = 0;

   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_Real dot = 0.0;
      SCIP_Real vdotray;

      if( SCIPisZero(scip, eigenvalues[i]) )
      {
         if( wcoefs == NULL )
            continue;
      }
      else
      {
         dot = vzlp[i] + vb[i] / (2.0 * (sidefactor * eigenvalues[i]));
      }
      vdotray = computeEigenvecDotRay(&eigenvectors[i * nquadexprs], nquadexprs, raycoefs, rayidx, raynnonz);

      if( SCIPisZero(scip, eigenvalues[i]) )
      {
         /* zero eigenvalue (and wcoefs not null) -> case 4: compute w(r) */
         assert(wcoefs != NULL);
         wray += vb[i] * vdotray;
#ifdef INTERCUT_MOREDEBUG
         printf(" wray += %g, vb[i] %g and vdotray %g\n", vb[i] * vdotray,vb[i],vdotray);
#endif
      }
      else if( sidefactor * eigenvalues[i] > 0 )
      {
         /* positive eigenvalue: compute common part of D and E */
         *d += (sidefactor * eigenvalues[i]) * dot * vdotray;
         *e += (sidefactor * eigenvalues[i]) * SQR( dot );

#ifdef INTERCUT_MOREDEBUG
         printf("Positive eigenvalue: computing D: v^T ray %g, v^T( zlp + b/theta ) %g and theta %g \n", vdotray, dot, (sidefactor * eigenvalues[i]));
#endif
      }
      else
      {
         /* negative eigenvalue: compute common part of A, B, and C */
         *a -= (sidefactor * eigenvalues[i]) * SQR( vdotray );
         *b -= 2.0 * (sidefactor * eigenvalues[i]) *  dot * vdotray;
         *c -= (sidefactor * eigenvalues[i]) * SQR( dot );

#ifdef INTERCUT_MOREDEBUG
         printf("Negative eigenvalue: computing A: v^T ray %g, and theta %g \n", vdotray, (sidefactor * eigenvalues[i]));
#endif
      }
   }

   if( ! iscase4 )
   {
      /* We are in one of the first 3 cases */
      *e += MAX(kappa, 0.0);
      *c -= MIN(kappa, 0.0);

      /* finish computation of D and E */
      assert(*e > 0);
      *e = SQRT( *e );
      *d /= *e;

      /* some sanity checks only applicable to these cases (more at the end) */
      assert(*c >= 0);

      /* In theory, the function at 0 must be negative. Because of bad numerics this might not always hold, so we abort
       * the generation of the cut in this case.
       */
      if( SQRT( *c ) - *e >= 0 )
      {
         /* check if it's really a numerical problem */
         assert(SQRT( *c ) > 10e+15 || *e > 10e+15 || SQRT( *c ) - *e < 10e+9);

         INTERLOG(printf("Bad numerics: phi(0) >= 0\n"); )
         *success = FALSE;
         return SCIP_OKAY;
      }
   }
   else
   {
      SCIP_Real norm;
      SCIP_Real xextra;
      SCIP_Real yextra;

      norm = SQRT( 1 + SQR( kappa ) );
      xextra = wzlp + kappa + norm;
      yextra = wzlp + kappa - norm;

      /* finish computing w(ray), the linear part is missing */
      wray += computeWRayLinear(nlhdlrexprdata, sidefactor, raycoefs, rayidx, raynnonz);

      /*
       * coefficients of case 4b
       */
      /* at this point E is \|x(zlp)\|^2, so we can finish A, B, and C */
      coefs4b[0] = (*a) * (*e);
      coefs4b[1] = (*b) * (*e);
      coefs4b[2] = (*c) * (*e);

      /* finish D and E */
      coefs4b[3] = *d;
      coefs4b[4] = (*e) + xextra / 2.0;

      /* when \|x(zlp)\|^2 is too large, we can divide everything by \|x(zlp)\| */
      if( *e > 100 )
      {
         coefs4b[0] = (*a);
         coefs4b[1] = (*b);
         coefs4b[2] = (*c);

         /* finish D and E */
         coefs4b[3] = *d / SQRT( *e );
         coefs4b[4] = SQRT( *e ) + (xextra / (2.0 * SQRT( *e )));
      }

      /*
       * coefficients of case 4a
       */
      /* finish A, B, and C */
      *a += SQR( wray ) / (4.0 * norm);
      *b += 2.0 * yextra * (wray) / (4.0 * norm);
      *c += SQR( yextra ) / (4.0 * norm);

      /* finish D and E */
      *e +=  SQR( xextra ) / (4.0 * norm);
      *e = SQRT( *e );

      *d += xextra * (wray) / (4.0 * norm);
      *d /= *e;

      /*
       * coefficients of condition: stores -numerator of x_{r+1}/ norm xhat, w(ray), and numerator of y_{s+1} at zlp
       *
       */
      /* at this point E is \| \hat{x} (zlp)\| */
      coefscondition[0] = - xextra / (*e);
      coefscondition[1] = wray;
      coefscondition[2] = yextra;
   }

#ifdef  DEBUG_INTERSECTIONCUT
   if( ! iscase4 )
   {
      SCIPinfoMessage(scip, NULL, "Restriction yields case 1,2 or 3: a,b,c,d,e %g %g %g %g %g\n", coefs1234a[0], coefs1234a[1], coefs1234a[2],
            coefs1234a[3], coefs1234a[4]);
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "Restriction yields\n   Case 4a: a,b,c,d,e %g %g %g %g %g\n", coefs1234a[0],
            coefs1234a[1], coefs1234a[2], coefs1234a[3], coefs1234a[4]);
      SCIPinfoMessage(scip, NULL, "   Case 4b: a,b,c,d,e %g %g %g %g %g\n", coefs4b[0], coefs4b[1], coefs4b[2],
            coefs4b[3], coefs4b[4]);
      SCIPinfoMessage(scip, NULL, "   Condition: xextra/e, wray, yextra %g %g %g g\n", coefscondition[0],
            coefscondition[1], coefscondition[2]);
   }
#endif

   /* some sanity check applicable to all cases */
   assert(*a >= 0); /* the function inside the root is convex */
   assert(*c >= 0); /* radicand at zero */

   if( iscase4 )
   {
      assert(coefs4b[0] >= 0); /* the function inside the root is convex */
      assert(coefs4b[2] >= 0); /* radicand at zero */
   }
   /*assert(4.0 * (*a) * (*c) >= SQR( *b ) ); *//* the function is defined everywhere, so minimum of radicand must be nonnegative */

   return SCIP_OKAY;
}

/** returns phi(zlp + t * ray) = SQRT(A t^2 + B t + C) - (D t + E) */
static
SCIP_Real evalPhiAtRay(
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

   assert(t != 1e20 || QUAD_TO_DBL(tmp) <= 0);

   return  QUAD_TO_DBL(tmp);
#else
   return SQRT( a * t * t + b * t + c ) - ( d * t + e );
#endif
}

/** checks whether case 4a applies
 *
 * The condition for being in case 4a is
 * \f[ -\lambda_{r+1} \Vert \hat y(zlp + tsol\, ray)\Vert + \hat y_{s+1}(zlp + tsol\, ray) \leq 0\f]
 *
 * This reduces to
 * \f[ -num(\hat x_{r+1}(zlp)) \sqrt{A t^2 + B t + C} / E  + w(ray) \cdot t + num(\hat y_{s+1}(zlp)) \leq 0\f]
 * where num is the numerator.
 */
static
SCIP_Real isCase4a(
   SCIP_Real             tsol,               /**< t in the above formula */
   SCIP_Real*            coefs4a,            /**< coefficients A, B, C, D, and E of case 4a */
   SCIP_Real*            coefscondition      /**< extra coefficients needed for the evaluation of the condition:
                                              *   \f$num(\hat x_{r+1}(zlp)) / E\f$; \f$w(ray)\f$; \f$num(\hat y_{s+1}(zlp))\f$ */
   )
{
   return (coefscondition[0] * SQRT( coefs4a[0] * SQR( tsol ) + coefs4a[1] * tsol + coefs4a[2] ) + coefscondition[1] *
         tsol + coefscondition[2]) <= 0.0;
}

/** helper function of computeRoot: we want phi to be &le; 0 */
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
      phival = evalPhiAtRay(curr, a, b, c, d, e);
#ifdef INTERCUT_MOREDEBUG
      printf("%d: lb,ub %.10f, %.10f. curr = %g -> phi at curr %g -> phi at lb %g \n", i, lb, ub, curr, phival, evalPhiAtRay(lb, a, b, c, d, e));
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

/** finds smallest positive root phi by finding the smallest positive root of
 * (A - D^2) t^2 + (B - 2 D*E) t + (C - E^2) = 0
 *
 * However, we are conservative and want a solution such that phi is negative, but close to 0.
 * Thus, we correct the result with a binary search.
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

#ifdef INTERCUT_MOREDEBUG
   {
      SCIP_Real binsol;
      binsol = SCIPinfinity(scip);
      doBinarySearch(scip, a, b, c, d, e, &binsol);
      printf("got root %g: with binsearch get %g\n", sol, binsol);
   }
#endif

   /* check that solution is acceptable, ideally it should be <= 0, however when it is positive, we trigger a binary
    * search to make it negative. This binary search might return a solution point that is not at accurately 0 as the
    * one obtained from the function above. Thus, it might fail to satisfy the condition of case 4b in some cases, e.g.,
    * ex8_3_1, bchoco05, etc
    */
   if( evalPhiAtRay(sol, a, b, c, d, e) <= 1e-10 )
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
      printf("do bin search because phival is %g\n", evalPhiAtRay(sol, a, b, c, d, e));
#endif
      doBinarySearch(scip, a, b, c, d, e, &sol);
   }

   return sol;
}

/** The maximal S-free set is \f$\gamma(z) \leq 0\f$; we find the intersection point of the ray `ray` starting from zlp with the
 * boundary of the S-free set.
 * That is, we find \f$t \geq 0\f$ such that \f$\gamma(zlp + t \cdot \text{ray}) = 0\f$.
 *
 * In cases 1,2, and 3, gamma is of the form
 *    \f[ \gamma(zlp + t \cdot \text{ray}) = \sqrt{A t^2 + B t + C} - (D t + E) \f]
 *
 * In the case 4 gamma is of the form
 *    \f[ \gamma(zlp + t \cdot \text{ray}) =
 *      \begin{cases}
 *        \sqrt{A t^2 + B t + C} - (D t + E), & \text{if some condition holds}, \\
 *        \sqrt{A' t^2 + B' t + C'} - (D' t + E'), & \text{otherwise.}
 *      \end{cases}
 *    \f]
 *
 * It can be shown (given the special properties of \f$\gamma\f$) that the smallest positive root of each function of the form
 * \f$\sqrt{a t^2 + b t + c} - (d t + e)\f$
 * is the same as the smallest positive root of the quadratic equation:
 * \f{align}{
 *     & \sqrt{a t^2 + b t + c} - (d t + e)) (\sqrt{a t^2 + b t + c} + (d t + e)) = 0 \\  \Leftrightarrow
 *     & (a - d^2) t^2 + (b - 2 d\,e) t + (c - e^2) = 0
 * \f}
 *
 * So, in cases 1, 2, and 3, this function just returns the solution of the above equation.
 * In case 4, it first solves the equation assuming we are in the first piece.
 * If there is no solution, then the second piece can't have a solution (first piece &ge; second piece for all t)
 * Then we check if the solution satisfies the condition.
 * If it doesn't then we solve the equation for the second piece.
 * If it has a solution, then it _has_ to be the solution.
 */
static
SCIP_Real computeIntersectionPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nlhdlr data */
   SCIP_Bool             iscase4,            /**< whether we are in case 4 or not */
   SCIP_Real*            coefs1234a,         /**< values of A, B, C, D, and E of cases 1, 2, 3, or 4a */
   SCIP_Real*            coefs4b,            /**< values of A, B, C, D, and E of case 4b */
   SCIP_Real*            coefscondition      /**< values needed to evaluate condition of case 4 */
   )
{
   SCIP_Real sol1234a;
   SCIP_Real sol4b;

   assert(coefs1234a != NULL);

#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "Computing intersection point for case 4? %d\n", iscase4);
#endif
   if( ! iscase4 )
      return computeRoot(scip, coefs1234a);

   assert(coefs4b != NULL);
   assert(coefscondition != NULL);

   /* compute solution of first piece */
#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "Compute root in 4a\n");
#endif
   sol1234a = computeRoot(scip, coefs1234a);

   /* if there is no solution --> second piece doesn't have solution */
   if( SCIPisInfinity(scip, sol1234a) )
   {
      /* this assert fails on multiplants_mtg5 the problem is that sqrt(A) <= D in 4a but not in 4b,
       * now, this is impossible since the phi4a >= phi4b, so actually sqrt(A) is 10e-15 away from
       * D in 4b
       */
      /* assert(SCIPisInfinity(scip, computeRoot(scip, coefs4b))); */
      return sol1234a;
   }

   /* if solution of 4a is in 4a, then return */
   if( isCase4a(sol1234a, coefs1234a, coefscondition) )
      return sol1234a;

#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "Root not in 4a -> Compute root in 4b\n");
#endif

   /* not on 4a --> then the intersection point is whatever 4b says: as phi4a >= phi4b, the solution of phi4b should
    * always be larger (but shouldn't be equal at this point given that isCase4a failed, and the condition function
    * evaluates to 0 when phi4a == phi4b) than the solution of phi4a; However, because of numerics (or limits in the
    * binary search) we can find a slightly smaller solution; thus, we just keep the larger one
    */
   sol4b = computeRoot(scip, coefs4b);

   /* this assert fails in many instances, e.g. water, because sol4b < sol1234a  */
   /* assert(SCIPisInfinity(scip, sol4b) || !isCase4a(sol4b, coefs1234a, coefscondition)); */
   /* count number of times we could have improved the coefficient by 10% */
   if( sol4b < sol1234a && evalPhiAtRay(1.1 * sol1234a, coefs4b[0], coefs4b[1], coefs4b[2], coefs4b[3], coefs4b[4]) <=
         0.0 )
      nlhdlrdata->ncouldimprovedcoef++;

   return MAX(sol1234a, sol4b);
}

/** checks if numerics of the coefficients are not too bad */
static
SCIP_Bool areCoefsNumericsGood(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nlhdlr data */
   SCIP_Real*            coefs1234a,         /**< coefficients for case 1-3 and 4a */
   SCIP_Real*            coefs4b,            /**< coefficients for case 4b */
   SCIP_Bool             iscase4             /**< whether we are in case 4 */
   )
{
   SCIP_Real max;
   SCIP_Real min;
   int j;

   /* check at phi at 0 is negative (note; this could be checked before restricting to the ray) also, if this
    * succeeds for one ray, it should suceed for every ray
    */
   if( SQRT( coefs1234a[2] ) - coefs1234a[4] >= 0.0 )
   {
      INTERLOG(printf("Bad numerics: phi(0) >= 0\n"); )
      nlhdlrdata->nphinonneg++;
      return FALSE;
   }

   /* maybe we want to avoid a large dynamism between A, B and C */
   if( nlhdlrdata->ignorebadrayrestriction )
   {
      max = 0.0;
      min = SCIPinfinity(scip);
      for( j = 0; j < 3; ++j )
      {
         SCIP_Real absval;

         absval = REALABS(coefs1234a[j]);
         if( max < absval )
            max = absval;
         if( absval != 0.0 && absval < min )
            min = absval;
      }

      if( SCIPisHugeValue(scip, max / min) )
      {
         INTERLOG(printf("Bad numerics 1 2 3 or 4a: max(A,B,C)/min(A,B,C) is too large (%g)\n", max / min); )
         nlhdlrdata->nbadrayrestriction++;
         return FALSE;
      }

      if( iscase4 )
      {
         max = 0.0;
         min = SCIPinfinity(scip);
         for( j = 0; j < 3; ++j )
         {
            SCIP_Real absval;

            absval = ABS(coefs4b[j]);
            if( max < absval )
               max = absval;
            if( absval != 0.0 && absval < min )
               min = absval;
         }

         if( SCIPisHugeValue(scip, max / min) )
         {
            INTERLOG(printf("Bad numeric 4b: max(A,B,C)/min(A,B,C) is too large (%g)\n", max / min); )
            nlhdlrdata->nbadrayrestriction++;
            return FALSE;
         }
      }
   }

   return TRUE;
}

/** computes intersection cut cuts off sol (because solution sol violates the quadratic constraint cons)
 * and stores it in rowprep. Here, we don't use any strengthening.
 */
static
SCIP_RETCODE computeIntercut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nlhdlr data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   RAYS*                 rays,               /**< rays */
   SCIP_Real             sidefactor,         /**< 1.0 if the violated constraint is q &le; rhs, -1.0 otherwise */
   SCIP_Bool             iscase4,            /**< whether we are in case 4 */
   SCIP_Real*            vb,                 /**< array containing \f$v_i^T b\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            vzlp,               /**< array containing \f$v_i^T zlp_q\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            wcoefs,             /**< coefficients of w for the qud vars or NULL if w is 0 */
   SCIP_Real             wzlp,               /**< value of w at zlp */
   SCIP_Real             kappa,              /**< value of kappa */
   SCIP_ROWPREP*         rowprep,            /**< rowprep for the generated cut */
   SCIP_Real*            interpoints,        /**< array to store intersection points for all rays or NULL if nothing
                                                  needs to be stored */
   SCIP_SOL*             sol,                /**< solution we want to separate */
   SCIP_Bool*            success             /**< if a cut candidate could be computed */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   int i;

   cols = SCIPgetLPCols(scip);
   rows = SCIPgetLPRows(scip);

   /* for every ray: compute cut coefficient and add var associated to ray into cut */
   for( i = 0; i < rays->nrays; ++i )
   {
      SCIP_Real interpoint;
      SCIP_Real cutcoef;
      int lppos;
      SCIP_Real coefs1234a[5];
      SCIP_Real coefs4b[5];
      SCIP_Real coefscondition[3];

      /* restrict phi to ray */
      SCIP_CALL( computeRestrictionToRay(scip, nlhdlrexprdata, sidefactor, iscase4,
               &rays->rays[rays->raysbegin[i]], &rays->raysidx[rays->raysbegin[i]], rays->raysbegin[i + 1] -
               rays->raysbegin[i], vb, vzlp, wcoefs, wzlp, kappa, coefs1234a, coefs4b, coefscondition, success) );

      if( ! *success )
         return SCIP_OKAY;

      /* if restriction to ray is numerically nasty -> abort cut separation */
      *success = areCoefsNumericsGood(scip, nlhdlrdata, coefs1234a, coefs4b, iscase4);

      if( ! *success )
         return SCIP_OKAY;

      /* compute intersection point */
      interpoint = computeIntersectionPoint(scip, nlhdlrdata, iscase4, coefs1234a, coefs4b, coefscondition);

#ifdef  DEBUG_INTERSECTIONCUT
      SCIPinfoMessage(scip, NULL, "interpoint for ray %d is %g\n", i, interpoint);
#endif

      /* store intersection point */
      if( interpoints != NULL )
         interpoints[i] = interpoint;

      /* compute cut coef */
      cutcoef = SCIPisInfinity(scip, interpoint) ? 0.0 : 1.0 / interpoint;

      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      lppos = rays->lpposray[i];
      if( lppos < 0 )
      {
         lppos = -lppos - 1;

         assert(SCIProwGetBasisStatus(rows[lppos]) == SCIP_BASESTAT_LOWER || SCIProwGetBasisStatus(rows[lppos]) ==
               SCIP_BASESTAT_UPPER);

         SCIP_CALL( addRowToCut(scip, rowprep, SCIProwGetBasisStatus(rows[lppos]) == SCIP_BASESTAT_UPPER ? cutcoef :
                  -cutcoef, rows[lppos], success) ); /* rows have flipper base status! */

         if( ! *success )
         {
            INTERLOG(printf("Bad numeric: now not nonbasic enough\n");)
            nlhdlrdata->nbadnonbasic++;
            return SCIP_OKAY;
         }
      }
      else
      {
         if( ! nlhdlrdata->useboundsasrays )
         {
            assert(SCIPcolGetBasisStatus(cols[lppos]) == SCIP_BASESTAT_UPPER || SCIPcolGetBasisStatus(cols[lppos]) ==
                  SCIP_BASESTAT_LOWER);
            SCIP_CALL( addColToCut(scip, rowprep, sol, SCIPcolGetBasisStatus(cols[lppos]) == SCIP_BASESTAT_UPPER ? -cutcoef :
                  cutcoef, cols[lppos]) );
         }
         else
         {
            SCIP_CALL( addColToCut(scip, rowprep, sol, rays->rays[i] == -1 ? -cutcoef : cutcoef, cols[lppos]) );
         }
      }
   }

   return SCIP_OKAY;
}

/** combine ray 1 and 2 to obtain new ray coef1 * ray1 - coef2 * ray2 in sparse format */
static
void combineRays(
   SCIP_Real*            raycoefs1,          /**< coefficients of ray 1 */
   int*                  rayidx1,            /**< index of consvar of the ray 1 coef is associated to */
   int                   raynnonz1,          /**< length of raycoefs1 and rayidx1 */
   SCIP_Real*            raycoefs2,          /**< coefficients of ray 2 */
   int*                  rayidx2,            /**< index of consvar of the ray 2 coef is associated to */
   int                   raynnonz2,          /**< length of raycoefs2 and rayidx2 */
   SCIP_Real*            newraycoefs,        /**< coefficients of combined ray */
   int*                  newrayidx,          /**< index of consvar of the combined ray coef is associated to */
   int*                  newraynnonz,        /**< pointer to length of newraycoefs and newrayidx */
   SCIP_Real             coef1,              /**< coef of ray 1 */
   SCIP_Real             coef2               /**< coef of ray 2 */
   )
{
   int idx1;
   int idx2;

   idx1 = 0;
   idx2 = 0;
   *newraynnonz = 0;

   while( idx1 < raynnonz1 || idx2 < raynnonz2 )
   {
      /* if the pointers look at different variables (or one already arrieved at the end), only one pointer can move
       * on
       */
      if( idx1 >= raynnonz1 || (idx2 < raynnonz2 && rayidx1[idx1] > rayidx2[idx2]) )
      {
         /*printf("case 1 \n"); */
         newraycoefs[*newraynnonz] = - coef2 * raycoefs2[idx2];
         newrayidx[*newraynnonz] = rayidx2[idx2];
         ++(*newraynnonz);
         ++idx2;
      }
      else if( idx2 >= raynnonz2 || rayidx1[idx1] < rayidx2[idx2] )
      {
         /*printf("case 2 \n"); */
         newraycoefs[*newraynnonz] = coef1 * raycoefs1[idx1];
         newrayidx[*newraynnonz] = rayidx1[idx1];
         ++(*newraynnonz);
         ++idx1;
      }
      /* if both pointers look at the same variable, just compute the difference and move both pointers */
      else if( rayidx1[idx1] == rayidx2[idx2] )
      {
         /*printf("case 3 \n"); */
         newraycoefs[*newraynnonz] = coef1 * raycoefs1[idx1] - coef2 * raycoefs2[idx2];
         newrayidx[*newraynnonz] = rayidx1[idx1];
         ++(*newraynnonz);
         ++idx1;
         ++idx2;
      }
   }
}

/** checks if two rays are linearly dependent */
static
SCIP_Bool raysAreDependent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            raycoefs1,          /**< coefficients of ray 1 */
   int*                  rayidx1,            /**< index of consvar of the ray 1 coef is associated to */
   int                   raynnonz1,          /**< length of raycoefs1 and rayidx1 */
   SCIP_Real*            raycoefs2,          /**< coefficients of ray 2 */
   int*                  rayidx2,            /**< index of consvar of the ray 2 coef is associated to */
   int                   raynnonz2,          /**< length of raycoefs2 and rayidx2 */
   SCIP_Real*            coef                /**< pointer to store coef (s.t. r1 = coef * r2) in case rays are
                                              *   dependent */
   )
{
   int i;

   /* cannot be dependent if they have different number of non-zero entries */
   if( raynnonz1 != raynnonz2 )
      return FALSE;

   *coef = 0.0;

   for( i = 0; i < raynnonz1; ++i )
   {
      /* cannot be dependent if different variables have non-zero entries */
      if( rayidx1[i] != rayidx2[i] ||
         (SCIPisZero(scip, raycoefs1[i]) && !SCIPisZero(scip, raycoefs2[i])) ||
         (!SCIPisZero(scip, raycoefs1[i]) && SCIPisZero(scip, raycoefs2[i])) )
      {
         return FALSE;
      }

      if( *coef != 0.0 )
      {
         /* cannot be dependent if the coefs aren't equal for all entries */
         if( ! SCIPisEQ(scip, *coef, raycoefs1[i] / raycoefs2[i]) )
         {
            return FALSE;
         }
      }
      else
         *coef = raycoefs1[i] / raycoefs2[i];
   }

   return TRUE;
}

/** checks if the ray alpha * ray_i + (1 - alpha) * ray_j is in the recession cone of the S-free set. To do so,
  * we check if phi restricted to the ray has a positive root.
  */
static
SCIP_RETCODE rayInRecessionCone(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nlhdlr data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   RAYS*                 rays,               /**< rays */
   int                   j,                  /**< index of current ray in recession cone */
   int                   i,                  /**< index of current ray not in recession cone */
   SCIP_Real             sidefactor,         /**< 1.0 if the violated constraint is q &le; rhs, -1.0 otherwise */
   SCIP_Bool             iscase4,            /**< whether we are in case 4 */
   SCIP_Real*            vb,                 /**< array containing \f$v_i^T b\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            vzlp,               /**< array containing \f$v_i^T zlp_q\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            wcoefs,             /**< coefficients of w for the quad vars or NULL if w is 0 */
   SCIP_Real             wzlp,               /**< value of w at zlp */
   SCIP_Real             kappa,              /**< value of kappa */
   SCIP_Real             alpha,              /**< coef for combining the two rays */
   SCIP_Bool*            inreccone,          /**< pointer to store whether the ray is in the recession cone or not */
   SCIP_Bool*            success             /**< Did numerical troubles occur? */
   )
{
   SCIP_Real coefs1234a[5];
   SCIP_Real coefs4b[5];
   SCIP_Real coefscondition[3];
   SCIP_Real interpoint;
   SCIP_Real* newraycoefs;
   int* newrayidx;
   int newraynnonz;

  *inreccone = FALSE;

   /* allocate memory for new ray */
   newraynnonz = (rays->raysbegin[i + 1] - rays->raysbegin[i]) + (rays->raysbegin[j + 1] - rays->raysbegin[j]);
   SCIP_CALL( SCIPallocBufferArray(scip, &newraycoefs, newraynnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newrayidx, newraynnonz) );

   /* build the ray alpha * ray_i + (1 - alpha) * ray_j */
   combineRays(&rays->rays[rays->raysbegin[i]], &rays->raysidx[rays->raysbegin[i]], rays->raysbegin[i + 1] -
           rays->raysbegin[i], &rays->rays[rays->raysbegin[j]], &rays->raysidx[rays->raysbegin[j]],
           rays->raysbegin[j + 1] - rays->raysbegin[j], newraycoefs, newrayidx, &newraynnonz, alpha,
           -1 + alpha);

   /* restrict phi to the "new" ray */
   SCIP_CALL( computeRestrictionToRay(scip, nlhdlrexprdata, sidefactor, iscase4, newraycoefs, newrayidx,
           newraynnonz, vb, vzlp, wcoefs, wzlp, kappa, coefs1234a, coefs4b, coefscondition, success) );

   if( ! *success )
      goto CLEANUP;

   /* check if restriction to "new" ray is numerically nasty. If so, treat the corresponding rho as if phi is
    * positive
    */

   /* compute intersection point */
   interpoint = computeIntersectionPoint(scip, nlhdlrdata, iscase4, coefs1234a, coefs4b, coefscondition);

   /* no root exists */
   if( SCIPisInfinity(scip, interpoint) )
      *inreccone = TRUE;

CLEANUP:
   SCIPfreeBufferArray(scip, &newrayidx);
   SCIPfreeBufferArray(scip, &newraycoefs);

   return SCIP_OKAY;
}

/** finds the smallest negative steplength for the current ray r_idx such that the combination
 * of r_idx with all rays not in the recession cone is in the recession cone
 */
static
SCIP_RETCODE findRho(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nlhdlr data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   RAYS*                 rays,               /**< rays */
   int                   idx,                /**< index of current ray we want to find rho for */
   SCIP_Real             sidefactor,         /**< 1.0 if the violated constraint is q &le; rhs, -1.0 otherwise */
   SCIP_Bool             iscase4,            /**< whether we are in case 4 */
   SCIP_Real*            vb,                 /**< array containing \f$v_i^T b\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            vzlp,               /**< array containing \f$v_i^T zlp_q\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            wcoefs,             /**< coefficients of w for the quad vars or NULL if w is 0 */
   SCIP_Real             wzlp,               /**< value of w at zlp */
   SCIP_Real             kappa,              /**< value of kappa */
   SCIP_Real*            interpoints,        /**< array to store intersection points for all rays or NULL if nothing
                                              *   needs to be stored */
   SCIP_Real*            rho,                /**< pointer to store the optimal rho */
   SCIP_Bool*            success             /**< could we successfully find the right rho? */
   )
{
   int i;

   /* go through all rays not in the recession cone and compute the largest negative steplength possible. The
    * smallest of them is then the steplength rho we use for the current ray */
   *rho = 0.0;
   for( i = 0; i < rays->nrays; ++i )
   {
      SCIP_Real currentrho;
      SCIP_Real coef;

      if( SCIPisInfinity(scip, interpoints[i]) )
         continue;

      /* if we cannot strengthen enough, we don't strengthen at all */
      if( SCIPisInfinity(scip, -*rho) )
      {
         *rho = -SCIPinfinity(scip);
         return SCIP_OKAY;
      }

      /* if the rays are linearly independent, we don't need to search for rho */
      if( raysAreDependent(scip, &rays->rays[rays->raysbegin[i]], &rays->raysidx[rays->raysbegin[i]],
          rays->raysbegin[i + 1] - rays->raysbegin[i], &rays->rays[rays->raysbegin[idx]],
          &rays->raysidx[rays->raysbegin[idx]], rays->raysbegin[idx + 1] - rays->raysbegin[idx], &coef) )
      {
         currentrho = coef * interpoints[i];
      }
      else
      {
         /*  since the two rays are linearly independent, we need to find the biggest alpha such that
          *  alpha * ray_i + (1 - alpha) * ray_idx in the recession cone is. For every ray i, we compute
          *  such a alpha but take the smallest one of them. We use "maxalpha" to keep track of this.
          *  Since we know that we can only use alpha < maxalpha, we don't need to do the whole binary search
          *  for every ray i. We only need to search the intervall [0, maxalpha]. Thereby, we start by checking
          *  if alpha = maxalpha is already feasable */

         SCIP_Bool inreccone;
         SCIP_Real alpha;
         SCIP_Real lb;
         SCIP_Real ub;
         int j;

         lb = 0.0;
         ub = 1.0;
         for( j = 0; j < BINSEARCH_MAXITERS; ++j )
         {
            alpha = (lb + ub) / 2.0;

            if( SCIPisZero(scip, alpha) )
            {
               alpha = 0.0;
               break;
            }

            SCIP_CALL( rayInRecessionCone(scip, nlhdlrdata, nlhdlrexprdata, rays, idx, i, sidefactor, iscase4, vb,
                  vzlp, wcoefs, wzlp, kappa, alpha, &inreccone, success) );

            if( ! *success )
               return SCIP_OKAY;

            /* no root exists */
            if( inreccone )
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
            currentrho = (alpha - 1) * interpoints[i] / alpha;  /*lint !e795*/
      }

      if( currentrho < *rho )
         *rho = currentrho;

      if( *rho < -10e+06 )
         *rho = -SCIPinfinity(scip);

      /* if rho is too small, don't add it */
      if( SCIPisZero(scip, *rho) )
         *success = FALSE;
   }

   return SCIP_OKAY;
}

/** computes intersection cut using negative edge extension to strengthen rays that do not intersect
 * (i.e., rays in the recession cone)
 */
static
SCIP_RETCODE computeStrengthenedIntercut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nlhdlr data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   RAYS*                 rays,               /**< rays */
   SCIP_Real             sidefactor,         /**< 1.0 if the violated constraint is q &le; rhs, -1.0 otherwise */
   SCIP_Bool             iscase4,            /**< whether we are in case 4 */
   SCIP_Real*            vb,                 /**< array containing \f$v_i^T b\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            vzlp,               /**< array containing \f$v_i^T zlp_q\f$ for \f$i \in I_+ \cup I_-\f$ */
   SCIP_Real*            wcoefs,             /**< coefficients of w for the qud vars or NULL if w is 0 */
   SCIP_Real             wzlp,               /**< value of w at zlp */
   SCIP_Real             kappa,              /**< value of kappa */
   SCIP_ROWPREP*         rowprep,            /**< rowprep for the generated cut */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_Bool*            success,            /**< if a cut candidate could be computed */
   SCIP_Bool*            strengthsuccess     /**< if strengthening was successfully applied */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_Real* interpoints;
   SCIP_Real avecutcoef;
   int counter;
   int i;

   *success = TRUE;
   *strengthsuccess = FALSE;

   cols = SCIPgetLPCols(scip);
   rows = SCIPgetLPRows(scip);

   /* allocate memory for intersection points */
   SCIP_CALL( SCIPallocBufferArray(scip, &interpoints, rays->nrays) );

   /* compute all intersection points and store them in interpoints; build not-stregthened intersection cut */
   SCIP_CALL( computeIntercut(scip, nlhdlrdata, nlhdlrexprdata, rays, sidefactor, iscase4, vb, vzlp, wcoefs, wzlp, kappa,
            rowprep, interpoints, sol, success) );

   if( ! *success )
      goto CLEANUP;

   /* keep track of the number of attempted strengthenings and average cutcoef */
   counter = 0;
   avecutcoef = 0.0;

   /* go through all intersection points that are equal to infinity -> these correspond to the rays which are in the
    * recession cone of C, i.e. the rays for which we (possibly) can compute a negative steplength */
   for( i = 0; i < rays->nrays; ++i )
   {
      SCIP_Real rho;
      SCIP_Real cutcoef;
      int lppos;

      if( !SCIPisInfinity(scip, interpoints[i]) )
         continue;

      /* if we reached the limit of strengthenings, we stop */
      if( counter >= nlhdlrdata->nstrengthlimit )
         break;

      /* compute the smallest rho */
      SCIP_CALL( findRho(scip, nlhdlrdata, nlhdlrexprdata, rays, i, sidefactor, iscase4, vb, vzlp, wcoefs, wzlp, kappa,
               interpoints, &rho, success));

      /* compute cut coef */
      if( ! *success  || SCIPisInfinity(scip, -rho) )
         cutcoef = 0.0;
      else
         cutcoef = 1.0 / rho;

      /* track average cut coef */
      counter += 1;
      avecutcoef += cutcoef;

      if( ! SCIPisZero(scip, cutcoef) )
         *strengthsuccess = TRUE;

      /* add var to cut: if variable is nonbasic at upper we have to flip sign of cutcoef */
      lppos = rays->lpposray[i];
      if( lppos < 0 )
      {
         lppos = -lppos - 1;

         assert(SCIProwGetBasisStatus(rows[lppos]) == SCIP_BASESTAT_LOWER || SCIProwGetBasisStatus(rows[lppos]) ==
               SCIP_BASESTAT_UPPER);

         SCIP_CALL( addRowToCut(scip, rowprep, SCIProwGetBasisStatus(rows[lppos]) == SCIP_BASESTAT_UPPER ? cutcoef :
                  -cutcoef, rows[lppos], success) ); /* rows have flipper base status! */

         if( ! *success )
         {
            INTERLOG(printf("Bad numeric: row not nonbasic enough\n");)
            nlhdlrdata->nbadnonbasic++;
            return SCIP_OKAY;
         }
      }
      else
      {
         assert(SCIPcolGetBasisStatus(cols[lppos]) == SCIP_BASESTAT_UPPER || SCIPcolGetBasisStatus(cols[lppos]) ==
               SCIP_BASESTAT_LOWER);
         SCIP_CALL( addColToCut(scip, rowprep, sol, SCIPcolGetBasisStatus(cols[lppos]) == SCIP_BASESTAT_UPPER ? -cutcoef :
                  cutcoef, cols[lppos]) );
      }
   }

   if( counter > 0 )
      nlhdlrdata->cutcoefsum += avecutcoef / counter;

CLEANUP:
   SCIPfreeBufferArray(scip, &interpoints);

   return SCIP_OKAY;
}

/** sets variable in solution "vertex" to its nearest bound */
static
SCIP_RETCODE setVarToNearestBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_SOL*             vertex,             /**< new solution to separate */
   SCIP_VAR*             var,                /**< var we want to find nearest bound to */
   SCIP_Real*            factor,             /**< is vertex for current var at lower or upper? */
   SCIP_Bool*            success             /**< TRUE if no variable is bounded */
   )
{
   SCIP_Real solval;
   SCIP_Real bound;

   solval = SCIPgetSolVal(scip, sol, var);
   *success = TRUE;

   /* find nearest bound */
   if( SCIPisInfinity(scip, SCIPvarGetLbLocal(var)) && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   else if( solval - SCIPvarGetLbLocal(var) < SCIPvarGetUbLocal(var) - solval )
   {
      bound = SCIPvarGetLbLocal(var);
      *factor = 1.0;
   }
   else
   {
      bound = SCIPvarGetUbLocal(var);
      *factor = -1.0;
   }

   /* set val to bound in solution */
   SCIP_CALL( SCIPsetSolVal(scip, vertex, var, bound) );
   return SCIP_OKAY;
}

/** This function finds vertex (w.r.t. bounds of variables appearing in the quadratic) that is closest to the current
 * solution we want to separate.
 *
 * Furthermore, we store the rays corresponding to the unit vectors, i.e.,
 *    - if \f$x_i\f$ is at its lower bound in vertex --> \f$r_i =  e_i\f$
 *    - if \f$x_i\f$ is at its upper bound in vertex --> \f$r_i = -e_i\f$
 */
static
SCIP_RETCODE findVertexAndGetRays(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,      /**< nlhdlr expression data */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_SOL*             vertex,             /**< new 'vertex' (w.r.t. bounds) solution to separate */
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   RAYS**                raysptr,            /**< pointer to new bound rays */
   SCIP_Bool*            success             /**< TRUE if no variable is unbounded */
   )
{
   SCIP_EXPR* qexpr;
   SCIP_EXPR** linexprs;
   RAYS* rays;
   int nquadexprs;
   int nlinexprs;
   int raylength;
   int i;

   *success = TRUE;

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, NULL, NULL, NULL);

   raylength = (auxvar == NULL) ? nquadexprs + nlinexprs : nquadexprs + nlinexprs + 1;

   /* create rays */
   SCIP_CALL( createBoundRays(scip, raysptr, raylength) );
   rays = *raysptr;

   rays->rayssize = raylength;
   rays->nrays = raylength;

   /* go through quadratic variables */
   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_EXPR* expr;
      SCIPexprGetQuadraticQuadTerm(qexpr, i, &expr, NULL, NULL, NULL, NULL, NULL);

      rays->raysbegin[i] = i;
      rays->raysidx[i] = i;
      rays->lpposray[i] = SCIPcolGetLPPos(SCIPvarGetCol(SCIPgetExprAuxVarNonlinear(expr)));

      SCIP_CALL( setVarToNearestBound(scip, sol, vertex, SCIPgetExprAuxVarNonlinear(expr),
         &rays->rays[i], success) );

      if( ! *success )
         return SCIP_OKAY;
   }

   /* go through linear variables */
   for( i = 0; i < nlinexprs; ++i )
   {
      rays->raysbegin[i + nquadexprs] = i + nquadexprs;
      rays->raysidx[i + nquadexprs] = i + nquadexprs;
      rays->lpposray[i + nquadexprs] = SCIPcolGetLPPos(SCIPvarGetCol(SCIPgetExprAuxVarNonlinear(linexprs[i])));

      SCIP_CALL( setVarToNearestBound(scip, sol, vertex, SCIPgetExprAuxVarNonlinear(linexprs[i]),
         &rays->rays[i + nquadexprs], success) );

      if( ! *success )
         return SCIP_OKAY;
   }

   /* consider auxvar if it exists */
   if( auxvar != NULL )
   {
      rays->raysbegin[nquadexprs + nlinexprs] = nquadexprs + nlinexprs;
      rays->raysidx[nquadexprs + nlinexprs] = nquadexprs + nlinexprs;
      rays->lpposray[nquadexprs + nlinexprs] = SCIPcolGetLPPos(SCIPvarGetCol(auxvar));

      SCIP_CALL( setVarToNearestBound(scip, sol, vertex, auxvar, &rays->rays[nquadexprs + nlinexprs], success) );

      if( ! *success )
         return SCIP_OKAY;
   }

   rays->raysbegin[raylength] = raylength;

   return SCIP_OKAY;
}

/** checks if the quadratic constraint is violated by sol */
static
SCIP_Bool isQuadConsViolated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_VAR*             auxvar,             /**< aux var of expr or NULL if not needed (e.g. separating real cons) */
   SCIP_SOL*             sol,                /**< solution to check feasibility for */
   SCIP_Real             sidefactor          /**< 1.0 if the violated constraint is q &le; rhs, -1.0 otherwise */
   )
{
   SCIP_EXPR* qexpr;
   SCIP_EXPR** linexprs;
   SCIP_Real* lincoefs;
   SCIP_Real constant;
   SCIP_Real val;
   int nquadexprs;
   int nlinexprs;
   int nbilinexprs;
   int i;

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs,
      &nbilinexprs, NULL, NULL);

   val = 0.0;

   /* go through quadratic terms */
   for( i = 0; i < nquadexprs; i++ )
   {
      SCIP_EXPR* expr;
      SCIP_Real quadlincoef;
      SCIP_Real sqrcoef;
      SCIP_Real solval;

      SCIPexprGetQuadraticQuadTerm(qexpr, i, &expr, &quadlincoef, &sqrcoef, NULL, NULL, NULL);

      solval = SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr));

      /* add square term */
      val += sqrcoef * SQR(solval);

      /* add linear term */
      val += quadlincoef * solval;
   }

   /* go through bilinear terms */
   for( i = 0; i < nbilinexprs; i++ )
   {
      SCIP_EXPR* expr1;
      SCIP_EXPR* expr2;
      SCIP_Real bilincoef;

      SCIPexprGetQuadraticBilinTerm(qexpr, i, &expr1, &expr2, &bilincoef, NULL, NULL);

      val += bilincoef * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr1))
         * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr2));
   }

   /* go through linear terms */
   for( i = 0; i < nlinexprs; i++ )
   {
      val += lincoefs[i] * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(linexprs[i]));
   }

   /* add auxvar if exists and get constant */
   if( auxvar != NULL )
   {
      val -= SCIPgetSolVal(scip, sol, auxvar);

      constant = (sidefactor == 1.0) ? constant - SCIPgetRhsNonlinear(nlhdlrexprdata->cons) :
         SCIPgetLhsNonlinear(nlhdlrexprdata->cons) - constant;
   }
   else
      constant = (sidefactor * constant);

   val = (sidefactor * val);

   /* now constraint is q(z) <= const */
   if( val <= constant )
      return FALSE;
   else
      return TRUE;
}

/** generates intersection cut that cuts off sol (which violates the quadratic constraint cons) */
static
SCIP_RETCODE generateIntercut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expr */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nlhdlr data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_CONS*            cons,               /**< violated constraint that contains expr */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_ROWPREP*         rowprep,            /**< rowprep for the generated cut */
   SCIP_Bool             overestimate,       /**< TRUE if viol cons is q(z) &ge; lhs; FALSE if q(z) &le; rhs */
   SCIP_Bool*            success             /**< whether separation was successfull or not */
   )
{
   SCIP_EXPR* qexpr;
   RAYS* rays;
   SCIP_VAR* auxvar;
   SCIP_Real sidefactor;
   SCIP_Real* vb;      /* eigenvectors * b */
   SCIP_Real* vzlp;    /* eigenvectors * lpsol */
   SCIP_Real* wcoefs;  /* coefficients affecting quadterms in w */
   SCIP_Real wzlp;     /* w(lpsol) */
   SCIP_Real kappa;
   SCIP_Bool iscase4;
   SCIP_SOL* vertex;
   SCIP_SOL* soltoseparate;
   int nquadexprs;
   int nlinexprs;
   int i;

   /* count number of calls */
   (nlhdlrdata->ncalls++);

   qexpr = nlhdlrexprdata->qexpr;
   SCIPexprGetQuadraticData(qexpr, NULL, &nlinexprs, NULL, NULL, &nquadexprs, NULL, NULL, NULL);

#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "Generating intersection cut for quadratic expr %p aka", (void*)expr);
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   *success = TRUE;
   iscase4 = TRUE;

   /* in nonbasic space cut is >= 1 */
   assert(SCIProwprepGetSide(rowprep) == 0.0);
   SCIProwprepAddSide(rowprep, 1.0);
   SCIProwprepSetSidetype(rowprep, SCIP_SIDETYPE_LEFT);
   assert(SCIProwprepGetSide(rowprep) == 1.0);

   auxvar = (nlhdlrexprdata->cons != cons) ? SCIPgetExprAuxVarNonlinear(expr) : NULL;
   sidefactor = overestimate ? -1.0 : 1.0;

   rays = NULL;

   /* check if we use tableau or bounds as rays */
   if( ! nlhdlrdata->useboundsasrays )
   {
      SCIP_CALL( createAndStoreSparseRays(scip, nlhdlrexprdata, auxvar, &rays, success) );

      if( ! *success )
      {
         INTERLOG(printf("Failed to get rays: there is a var with base status ZERO!\n"); )
         return SCIP_OKAY;
      }

      soltoseparate = sol;
   }
   else
   {
      SCIP_Bool violated;

      if( auxvar != NULL )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      /* create new solution */
      SCIP_CALL( SCIPcreateSol(scip, &vertex, NULL) );

      /* find nearest vertex of the box to separate and compute rays */
      SCIP_CALL( findVertexAndGetRays(scip, nlhdlrexprdata, sol, vertex, auxvar, &rays, success) );

      if( ! *success )
      {
         INTERLOG(printf("Failed to use bounds as rays: variable is unbounded!\n"); )
         freeRays(scip, &rays);
         SCIP_CALL( SCIPfreeSol(scip, &vertex) );
         return SCIP_OKAY;
      }

      /* check if vertex is violated */
      violated = isQuadConsViolated(scip, nlhdlrexprdata, auxvar, vertex, sidefactor);

      if( ! violated )
      {
         INTERLOG(printf("Failed to use bounds as rays: nearest vertex is not violated!\n"); )
         freeRays(scip, &rays);
         SCIP_CALL( SCIPfreeSol(scip, &vertex) );
         *success = FALSE;
         return SCIP_OKAY;
      }

      soltoseparate = vertex;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vb, nquadexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vzlp, nquadexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &wcoefs, nquadexprs) );

   SCIP_CALL( intercutsComputeCommonQuantities(scip, nlhdlrexprdata, auxvar, sidefactor, soltoseparate, vb, vzlp, wcoefs, &wzlp, &kappa) );

   /* check if we are in case 4 */
   if( nlinexprs == 0 && auxvar == NULL )
   {
      for( i = 0; i < nquadexprs; ++i )
         if( wcoefs[i] != 0.0 )
            break;

      if( i == nquadexprs )
      {
         /* from now on wcoefs is going to be NULL --> case 1, 2 or 3 */
         SCIPfreeBufferArray(scip, &wcoefs);
         iscase4 = FALSE;
      }
   }

   /* compute (strengthened) intersection cut */
   if( nlhdlrdata->usestrengthening )
   {
      SCIP_Bool strengthsuccess;

      SCIP_CALL( computeStrengthenedIntercut(scip, nlhdlrdata, nlhdlrexprdata, rays, sidefactor, iscase4, vb, vzlp, wcoefs,
            wzlp, kappa, rowprep, soltoseparate, success, &strengthsuccess) );

      if( *success && strengthsuccess )
         nlhdlrdata->nstrengthenings++;
   }
   else
   {
      SCIP_CALL( computeIntercut(scip, nlhdlrdata, nlhdlrexprdata, rays, sidefactor, iscase4, vb, vzlp, wcoefs, wzlp, kappa,
            rowprep, NULL, soltoseparate, success) );
   }

   SCIPfreeBufferArrayNull(scip, &wcoefs);
   SCIPfreeBufferArray(scip, &vzlp);
   SCIPfreeBufferArray(scip, &vb);
   freeRays(scip, &rays);

   if( nlhdlrdata->useboundsasrays )
   {
      SCIP_CALL( SCIPfreeSol(scip, &vertex) );
   }

   return SCIP_OKAY;
}

/** returns whether a quadratic form is "propagable"
 *
 * It is propagable, if a variable (aka child expr) appears at least twice, which is the case if at least two of the following hold:
 * - it appears as a linear term (coef*expr)
 * - it appears as a square term (coef*expr^2)
 * - it appears in a bilinear term
 * - it appears in another bilinear term
 */
static
SCIP_Bool isPropagable(
   SCIP_EXPR*            qexpr               /**< quadratic representation data */
   )
{
   int nquadexprs;
   int i;

   SCIPexprGetQuadraticData(qexpr, NULL, NULL, NULL, NULL, &nquadexprs, NULL, NULL, NULL);

   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_Real lincoef;
      SCIP_Real sqrcoef;
      int nadjbilin;

      SCIPexprGetQuadraticQuadTerm(qexpr, i, NULL, &lincoef, &sqrcoef, &nadjbilin, NULL, NULL);

      if( (lincoef != 0.0) + (sqrcoef != 0.0) + nadjbilin >= 2 )  /*lint !e514*/ /* actually MIN(2, nadjbilin), but we check >= 2 */
         return TRUE;
   }

   return FALSE;
}

/** returns whether a quadratic term is "propagable"
 *
 * A term is propagable, if its variable (aka child expr) appears at least twice, which is the case if at least two of the following hold:
 * - it appears as a linear term (coef*expr)
 * - it appears as a square term (coef*expr^2)
 * - it appears in a bilinear term
 * - it appears in another bilinear term
 */
static
SCIP_Bool isPropagableTerm(
   SCIP_EXPR*            qexpr,              /**< quadratic representation data */
   int                   idx                 /**< index of quadratic term to consider */
   )
{
   SCIP_Real lincoef;
   SCIP_Real sqrcoef;
   int nadjbilin;

   SCIPexprGetQuadraticQuadTerm(qexpr, idx, NULL, &lincoef, &sqrcoef, &nadjbilin, NULL,  NULL);

   return (lincoef != 0.0) + (sqrcoef != 0.0) + nadjbilin >= 2;  /*lint !e514*/ /* actually MIN(2, nadjbilin), but we check >= 2 */
}

/** solves a quadratic equation \f$ a\, \text{expr}^2 + b\, \text{expr} \in \text{rhs} \f$ (with \f$b\f$ an interval)
 * and reduces bounds on `expr` or deduces infeasibility if possible
 */
static
SCIP_RETCODE propagateBoundsQuadExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression for which to solve */
   SCIP_Real             sqrcoef,            /**< square coefficient */
   SCIP_INTERVAL         b,                  /**< interval acting as linear coefficient */
   SCIP_INTERVAL         rhs,                /**< interval acting as rhs */
   SCIP_Bool*            infeasible,         /**< buffer to store if propagation produced infeasibility */
   int*                  nreductions         /**< buffer to store the number of interval reductions */
   )
{
   SCIP_INTERVAL a;
   SCIP_INTERVAL exprbounds;
   SCIP_INTERVAL newrange;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

#ifdef DEBUG_PROP
   SCIPinfoMessage(scip, NULL, "Propagating <expr> by solving a <expr>^2 + b <expr> in rhs, where <expr> is: ");
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "expr in [%g, %g], a = %g, b = [%g, %g] and rhs = [%g, %g]\n",
         SCIPintervalGetInf(SCIPgetExprBoundsNonlinear(scip, expr)),
         SCIPintervalGetSup(SCIPgetExprBoundsNonlinear(scip, expr)), sqrcoef, b.inf, b.sup,
         rhs.inf, rhs.sup);
#endif

   exprbounds = SCIPgetExprBoundsNonlinear(scip, expr);
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, exprbounds) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* compute solution of a*x^2 + b*x in rhs */
   SCIPintervalSet(&a, sqrcoef);
   SCIPintervalSolveUnivariateQuadExpression(SCIP_INTERVAL_INFINITY, &newrange, a, b, rhs, exprbounds);

#ifdef DEBUG_PROP
   SCIPinfoMessage(scip, NULL, "Solution [%g, %g]\n", newrange.inf, newrange.sup);
#endif

   SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, expr, newrange, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** solves a linear equation \f$ b\, \text{expr} \in \text{rhs} \f$ (with \f$b\f$ a scalar) and reduces bounds on `expr` or deduces infeasibility if possible */
static
SCIP_RETCODE propagateBoundsLinExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression for which to solve */
   SCIP_Real             b,                  /**< linear coefficient */
   SCIP_INTERVAL         rhs,                /**< interval acting as rhs */
   SCIP_Bool*            infeasible,         /**< buffer to store if propagation produced infeasibility */
   int*                  nreductions         /**< buffer to store the number of interval reductions */
   )
{
   SCIP_INTERVAL newrange;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

#ifdef DEBUG_PROP
   SCIPinfoMessage(scip, NULL, "Propagating <expr> by solving %g <expr> in [%g, %g], where <expr> is: ", b, rhs.inf, rhs.sup);
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   /* compute solution of b*x in rhs */
   SCIPintervalDivScalar(SCIP_INTERVAL_INFINITY, &newrange, rhs, b);

#ifdef DEBUG_PROP
   SCIPinfoMessage(scip, NULL, "Solution [%g, %g]\n", newrange.inf, newrange.sup);
#endif

   SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, expr, newrange, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** returns max of a/x - c*x for x in {x1, x2} with x1, x2 > 0 */
static
SCIP_Real computeMaxBoundaryForBilinearProp(
   SCIP_Real             a,                  /**< coefficient a */
   SCIP_Real             c,                  /**< coefficient c */
   SCIP_Real             x1,                 /**< coefficient x1 > 0 */
   SCIP_Real             x2                  /**< coefficient x2 > 0 */
   )
{
   SCIP_Real cneg;
   SCIP_Real cand1;
   SCIP_Real cand2;
   SCIP_ROUNDMODE roundmode;

   assert(x1 > 0.0);
   assert(x2 > 0.0);

   cneg = SCIPintervalNegateReal(c);

   roundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();
   cand1 = a/x1 + cneg*x1;
   cand2 = a/x2 + cneg*x2;
   SCIPintervalSetRoundingMode(roundmode);

   return MAX(cand1, cand2);
}

/** returns max of a/x - c*x for x in dom; it assumes that dom is contained in (0, +inf) */
static
SCIP_Real computeMaxForBilinearProp(
   SCIP_Real             a,                  /**< coefficient a */
   SCIP_Real             c,                  /**< coefficient c */
   SCIP_INTERVAL         dom                 /**< domain of x */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_INTERVAL argmax;
   SCIP_Real negunresmax;
   SCIP_Real boundarymax;
   assert(dom.inf > 0);

   /* if a >= 0, then the function is convex which means the maximum is at one of the boundaries
    *
    * if c = 0, then the function is monotone which means the maximum is also at one of the boundaries
    *
    * if a < 0, then the function is concave. The function then has a maximum if and only if there is a point with derivative 0,
    * that is, iff -a/x^2 - c = 0 has a solution; i.e. if -a/c >= 0, i.e. (using a<0 and c != 0), c > 0.
    * Otherwise (that is, c<0), the maximum is at one of the boundaries.
    */
   if( a >= 0.0 || c <= 0.0 )
      return computeMaxBoundaryForBilinearProp(a, c, dom.inf, dom.sup);

   /* now, the (unrestricted) maximum is at sqrt(-a/c).
    * if the argmax is not in the interior of dom then the solution is at a boundary, too
    * we check this by computing an interval that contains sqrt(-a/c) first
    */
   SCIPintervalSet(&argmax, -a);
   SCIPintervalDivScalar(SCIP_INTERVAL_INFINITY, &argmax, argmax, c);
   SCIPintervalSquareRoot(SCIP_INTERVAL_INFINITY, &argmax, argmax);

   /* if the interval containing sqrt(-a/c) does not intersect with the interior of dom, then
    * the (restricted) maximum is at a boundary (we could even say at which boundary, but that doesn't save much)
    */
   if( argmax.sup <= dom.inf || argmax.inf >= dom.sup )
      return computeMaxBoundaryForBilinearProp(a, c, dom.inf, dom.sup);

   /* the maximum at sqrt(-a/c) is -2*sqrt(-a*c), so we compute an upper bound for that by computing a lower bound for 2*sqrt(-a*c) */
   roundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();
   negunresmax = 2.0*SCIPnextafter(sqrt(SCIPintervalNegateReal(a)*c), 0.0);
   SCIPintervalSetRoundingMode(roundmode);

   /* if the interval containing sqrt(-a/c) is contained in dom, then we can return -negunresmax */
   if( argmax.inf >= dom.inf && argmax.sup <= dom.sup )
      return -negunresmax;

   /* now what is left is the case where we cannot say for sure whether sqrt(-a/c) is contained in dom or not
    * so we are conservative and return the max of both cases, i.e.,
    * the max of the upper bounds on -2*sqrt(-a*c), a/dom.inf-c*dom.inf, a/dom.sup-c*dom.sup.
    */
   boundarymax = computeMaxBoundaryForBilinearProp(a, c, dom.inf, dom.sup);
   return MAX(boundarymax, -negunresmax);
}

/** computes the range of rhs/x - coef * x for x in exprdom; this is used for the propagation of bilinear terms
 *
 * If 0 is in the exprdom, we set range to \f$\mathbb{R}\f$ (even though this is not quite correct, it is correct for the
 * intended use of the function).
 * TODO: maybe check before calling it whether 0 is in the domain and then just avoid calling it
 *
 * If rhs is [A,B] and x > 0, then we want the min of A/x - coef*x and max of B/x - coef*x for x in [exprdom].
 * If rhs is [A,B] and x < 0, then we want the min of B/x - coef*x and max of A/x - coef*x for x in [exprdom].
 * However, this is the same as min of -B/x + coef*x and max of -A/x + coef*x for x in -[exprdom].
 * Thus, we can always reduce to x > 0 by multiplying [exprdom], rhs, and coef by -1.
 */
static
void computeRangeForBilinearProp(
   SCIP_INTERVAL         exprdom,            /**< expression for which to solve */
   SCIP_Real             coef,               /**< expression for which to solve */
   SCIP_INTERVAL         rhs,                /**< rhs used for computation */
   SCIP_INTERVAL*        range               /**< storage for the resulting range */
   )
{
   SCIP_Real max;
   SCIP_Real min;

   if( exprdom.inf <= 0.0 && 0.0 <= exprdom.sup )
   {
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, range);
      return;
   }

   /* reduce to positive case */
   if( exprdom.sup < 0 )
   {
      SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &exprdom, exprdom, -1.0);
      SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &rhs, rhs, -1.0);
      coef *= -1.0;
   }
   assert(exprdom.inf > 0.0);

   /* compute maximum and minimum */
   max = computeMaxForBilinearProp(rhs.sup, coef, exprdom);
   min = -computeMaxForBilinearProp(-rhs.inf, -coef, exprdom);

   /* set interval */
   SCIPintervalSetBounds(range, min, max);
}

/** reverse propagates coef_i expr_i + constant in rhs */
static
SCIP_RETCODE reversePropagateLinearExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           linexprs,           /**< linear expressions */
   int                   nlinexprs,          /**< number of linear expressions */
   SCIP_Real*            lincoefs,           /**< coefficients of linear expressions */
   SCIP_Real             constant,           /**< constant */
   SCIP_INTERVAL         rhs,                /**< rhs */
   SCIP_Bool*            infeasible,         /**< buffer to store whether an exps' bounds were propagated to an empty interval */
   int*                  nreductions         /**< buffer to store the number of interval reductions of all exprs */
   )
{
   SCIP_INTERVAL* oldboundslin;
   SCIP_INTERVAL* newboundslin;
   int i;

   if( nlinexprs == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &oldboundslin, nlinexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newboundslin, nlinexprs) );

   for( i = 0; i < nlinexprs; ++i )
      oldboundslin[i] = SCIPexprGetActivity(linexprs[i]);  /* TODO use SCIPgetExprBoundsNonlinear(scip, linexprs[i]) ? */

   *nreductions = SCIPintervalPropagateWeightedSum(SCIP_INTERVAL_INFINITY, nlinexprs,
            oldboundslin, lincoefs, constant, rhs, newboundslin, infeasible);

   if( *nreductions > 0 && !*infeasible )
   {
      /* SCIP is more conservative with what constitutes a reduction than interval arithmetic so we follow SCIP */
      *nreductions = 0;
      for( i = 0; i < nlinexprs && ! (*infeasible); ++i )
      {
         SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, linexprs[i], newboundslin[i], infeasible, nreductions) );
      }
   }

   SCIPfreeBufferArray(scip, &newboundslin);
   SCIPfreeBufferArray(scip, &oldboundslin);

   return SCIP_OKAY;
}


/*
 * Callback methods of nonlinear handler
 */

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeexprdataQuadratic)
{  /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   if( (*nlhdlrexprdata)->quadactivities != NULL )
   {
      int nquadexprs;
      SCIPexprGetQuadraticData((*nlhdlrexprdata)->qexpr, NULL, NULL, NULL, NULL, &nquadexprs, NULL, NULL, NULL);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->quadactivities, nquadexprs);
   }

   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree
 *
 * A term is quadratic if
 * - it is a product expression of two expressions, or
 * - it is power expression of an expression with exponent 2.0.
 *
 * We define a _propagable_ quadratic expression as a quadratic expression whose termwise propagation does not yield the
 * best propagation. In other words, is a quadratic expression that suffers from the dependency problem.
 *
 * Specifically, a propagable quadratic expression is a sum expression such that there is at least one expr that appears
 * at least twice (because of simplification, this means it appears in a quadratic terms and somewhere else).
 * For example: \f$x^2 + y^2\f$ is not a propagable quadratic expression; \f$x^2 + x\f$ is a propagable quadratic expression;
 * \f$x^2 + x y\f$ is also a propagable quadratic expression
 *
 * Furthermore, we distinguish between propagable and non-propagable terms. A term is propagable if any of the expressions
 * involved in it appear somewhere else. For example, \f$xy + z^2 + z\f$ is a propagable quadratic, the term \f$xy\f$ is
 * non-propagable, and \f$z^2\f$ is propagable. For propagation, non-propagable terms are handled as if they were linear
 * terms, that is, we do not use the activity of \f$x\f$ and \f$y\f$ to compute the activity of \f$xy\f$ but rather we use directly
 * the activity of \f$xy\f$. Similarly, we do not backward propagate to \f$x\f$ and \f$y\f$ (the product expr handler will do this),
 * but we backward propagate to \f$x*y\f$. More technically, we register \f$xy\f$ for its activity usage, rather than\f$x\f$ and \f$y\f$.
 *
 * For propagation, we store the quadratic in our data structure in the following way: We count how often a variable
 * appears. Then, a bilinear product expr_i * expr_j is stored as expr_i * expr_j if # expr_i appears > # expr_j
 * appears. When # expr_i appears = # expr_j appears, it then it will be stored as expr_i * expr_j if and only if
 * expr_i < expr_j, where '<' is the expression order (see \ref EXPR_ORDER "Ordering Rules" in \ref scip_expr.h).
 * Heuristically, this should be useful for propagation. The intuition is that by factoring out the variable that
 * appears most often we should be able to take care of the dependency problem better.
 *
 * Simple convex quadratics like \f$x^2 + y^2\f$ are ignored since the default nlhdlr will take care of them.
 *
 * @note The expression needs to be simplified (in particular, it is assumed to be sorted).
 * @note Common subexpressions are also assumed to have been identified, the hashing will fail otherwise!
 *
 * Sorted implies that:
 *  - expr < expr^2: bases are the same, but exponent 1 < 2
 *  - expr < expr * other_expr: u*v < w holds if and only if v < w (OR8), but here w = u < v, since expr comes before
 *  other_expr in the product
 *  - expr < other_expr * expr: u*v < w holds if and only if v < w (OR8), but here v = w
 *
 *  Thus, if we see somebody twice, it is a propagable quadratic.
 *
 * It also implies that
 *  - expr^2 < expr * other_expr
 *  - other_expr * expr < expr^2
 *
 * It also implies that x^-2 < x^-1, but since, so far, we do not interpret x^-2 as (x^-1)^2, it is not a problem.
 */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectQuadratic)
{  /*lint --e{715,774}*/
   SCIP_NLHDLREXPRDATA* nlexprdata;
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_Real* eigenvalues;
   SCIP_Bool isquadratic;
   SCIP_Bool propagable;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);
   assert(nlhdlrexprdata != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* don't check if all enforcement methods are already ensured */
   if( (*enforcing & SCIP_NLHDLR_METHOD_ALL) == SCIP_NLHDLR_METHOD_ALL )
      return SCIP_OKAY;

   /* if it is not a sum of at least two terms, it is not interesting */
   /* TODO: constraints of the form l<= x*y <= r ? */
   if( ! SCIPisExprSum(scip, expr) || SCIPexprGetNChildren(expr) < 2 )
      return SCIP_OKAY;

   /* If we are in a subSCIP we don't want to separate intersection cuts */
   if( SCIPgetSubscipDepth(scip) > 0 )
      nlhdlrdata->useintersectioncuts = FALSE;

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "Nlhdlr quadratic detecting expr %p aka ", (void*)expr);
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "Have to enforce %d\n", *enforcing);
#endif

   /* check whether expression is quadratic (a sum with at least one square or bilinear term) */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, expr, &isquadratic) );

   /* not quadratic -> nothing for us */
   if( !isquadratic )
   {
      SCIPdebugMsg(scip, "expr %p is not quadratic -> abort detect\n", (void*)expr);
      return SCIP_OKAY;
   }

   propagable = isPropagable(expr);

   /* if we are not propagable and are in presolving, return */
   if( !propagable && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIPdebugMsg(scip, "expr %p is not propagable and in presolving -> abort detect\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* if we do not use intersection cuts and are not propagable, then we do not want to handle it at all;
    * if not propagable, then we need to check the curvature to decide if we want to generate intersection cuts
    */
   if( !propagable && !nlhdlrdata->useintersectioncuts )
   {
      SCIPdebugMsg(scip, "expr %p is not propagable -> abort detect\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* store quadratic in nlhdlrexprdata */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   nlexprdata = *nlhdlrexprdata;
   nlexprdata->qexpr = expr;
   nlexprdata->cons = cons;

#ifdef DEBUG_DETECT
   SCIPinfoMessage(scip, NULL, "Nlhdlr quadratic detected:\n");
   SCIP_CALL( SCIPprintExprQuadratic(scip, conshdlr, qexpr) );
#endif

   /* every propagable quadratic expression will be handled since we can propagate */
   if( propagable )
   {
      SCIP_EXPR** linexprs;
      int nlinexprs;
      int nquadexprs;
      int nbilin;
      int i;

      *participating |= SCIP_NLHDLR_METHOD_ACTIVITY;
      *enforcing |= SCIP_NLHDLR_METHOD_ACTIVITY;

      SCIPexprGetQuadraticData(expr, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, &nbilin, NULL, NULL);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlexprdata->quadactivities, nquadexprs) );

      /* notify children of quadratic that we will need their activity for propagation */
      for( i = 0; i < nlinexprs; ++i )
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, linexprs[i], FALSE, TRUE, FALSE, FALSE) );

      for( i = 0; i < nquadexprs; ++i )
      {
         SCIP_EXPR* argexpr;
         if( isPropagableTerm(expr, i) )
         {
            SCIPexprGetQuadraticQuadTerm(expr, i, &argexpr, NULL, NULL, &nbilin, NULL, NULL);
            SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, argexpr, FALSE, TRUE, FALSE, FALSE) );

#ifdef DEBUG_DETECT
            SCIPinfoMessage(scip, NULL, "quadterm %d propagable, using %p, unbounded=%d\n", i, (void*)argexpr, nbilin >
                  0 && SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(scip, argexpr)));
#endif
         }
         else
         {
            /* non-propagable quadratic is either a single square term or a single bilinear term
             * we should make use nlhdlrs in pow or product for this term, so we register usage of the square or product
             * expr instead of argexpr
             */
            SCIP_EXPR* sqrexpr;
            int* adjbilin;

            SCIPexprGetQuadraticQuadTerm(expr, i, &argexpr, NULL, NULL, &nbilin, &adjbilin, &sqrexpr);

            if( sqrexpr != NULL )
            {
               SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, sqrexpr, FALSE, TRUE, FALSE, FALSE) );
               assert(nbilin == 0);

#ifdef DEBUG_DETECT
               SCIPinfoMessage(scip, NULL, "quadterm %d non-propagable square, using %p\n", i, (void*)sqrexpr);
#endif
            }
            else
            {
               /* we have expr1 * other_expr or other_expr * expr1; know that expr1 is non propagable, but to decide if
                * we want the bounds of expr1 or of the product expr1 * other_expr (or other_expr * expr1), we have to
                * decide whether other_expr is also non propagable; due to the way we sort bilinear terms (by
                * frequency), we can deduce that other_expr doesn't appear anywhere else (i.e. is non propagable) if the
                * product is of the form expr1 * other_expr; however, if we see other_expr * expr1 we need to find
                * other_expr and check whether it is propagable
                */
               SCIP_EXPR* expr1;
               SCIP_EXPR* prodexpr;

               assert(nbilin == 1);
               SCIPexprGetQuadraticBilinTerm(expr, adjbilin[0], &expr1, NULL, NULL, NULL, &prodexpr);

               if( expr1 == argexpr )
               {
                  SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, prodexpr, FALSE, TRUE, FALSE, FALSE) );

#ifdef DEBUG_DETECT
                  SCIPinfoMessage(scip, NULL, "quadterm %d non-propagable product, using %p\n", i, (void*)prodexpr);
#endif
               }
               else
               {
                  int j;
                  /* check if other_expr is propagable in which case we need the bounds of expr1; otherwise we just need
                   * the bounds of the product and this will be (or was) registered when the loop takes us to the
                   * quadexpr other_expr.
                   * TODO this should be done faster, maybe store pos1 in bilinexprterm or store quadexprterm's in bilinexprterm
                   */
                  for( j = 0; j < nquadexprs; ++j )
                  {
                     SCIP_EXPR* exprj;
                     SCIPexprGetQuadraticQuadTerm(expr, j, &exprj, NULL, NULL, NULL, NULL, NULL);
                     if( expr1 == exprj )
                     {
                        if( isPropagableTerm(expr, j) )
                        {
                           SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, argexpr, FALSE, TRUE, FALSE, FALSE) );
#ifdef DEBUG_DETECT
                           SCIPinfoMessage(scip, NULL, "quadterm %d non-propagable alien product, using %p\n", i, (void*)argexpr);
#endif
                        }
                        break;
                     }
                  }
               }
            }
         }
      }
   }

   /* check if we are going to separate or not */
   nlexprdata->curvature = SCIP_EXPRCURV_UNKNOWN;

   /* for now, we do not care about separation if it is not required */
   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) == SCIP_NLHDLR_METHOD_SEPABOTH )
   {
      /* if nobody can do anything, remove data */
      if( *participating == SCIP_NLHDLR_METHOD_NONE ) /*lint !e845*/
      {
         SCIP_CALL( nlhdlrFreeexprdataQuadratic(scip, nlhdlr, expr, nlhdlrexprdata) );
      }
      else
      {
         SCIPdebugMsg(scip, "expr %p is quadratic and propagable -> propagate\n", (void*)expr);
      }
      return SCIP_OKAY;
   }

   assert(SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE);  /* separation should only be required in (init)solving stage */

   /* check if we can do something more: check curvature of quadratic function stored in nlexprdata
    * this is currently only used to decide whether we want to separate, so it can be skipped if in presolve
    */
   SCIPdebugMsg(scip, "checking curvature of expr %p\n", (void*)expr);
   SCIP_CALL( SCIPcomputeExprQuadraticCurvature(scip, expr, &nlexprdata->curvature, NULL, nlhdlrdata->useintersectioncuts) );

   /* get eigenvalues to be able to check whether they were computed */
   SCIPexprGetQuadraticData(expr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);

   /* if we use intersection cuts then we can handle any non-convex quadratic */
   if( nlhdlrdata->useintersectioncuts && eigenvalues != NULL && (*enforcing & SCIP_NLHDLR_METHOD_SEPABELOW) ==
         FALSE && nlexprdata->curvature != SCIP_EXPRCURV_CONVEX )
   {
      *participating |= SCIP_NLHDLR_METHOD_SEPABELOW;
   }

   if( nlhdlrdata->useintersectioncuts && eigenvalues != NULL && (*enforcing & SCIP_NLHDLR_METHOD_SEPAABOVE) == FALSE &&
         nlexprdata->curvature != SCIP_EXPRCURV_CONCAVE )
   {
      *participating |= SCIP_NLHDLR_METHOD_SEPAABOVE;
   }

   /* if nobody can do anything, remove data */
   if( *participating == SCIP_NLHDLR_METHOD_NONE ) /*lint !e845*/
   {
      SCIP_CALL( nlhdlrFreeexprdataQuadratic(scip, nlhdlr, expr, nlhdlrexprdata) );
      return SCIP_OKAY;
   }

   /* we only need auxiliary variables if we are going to separate */
   if( *participating & SCIP_NLHDLR_METHOD_SEPABOTH )
   {
      SCIP_EXPR** linexprs;
      int nquadexprs;
      int nlinexprs;
      int i;

      SCIPexprGetQuadraticData(expr, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, NULL, NULL, NULL);

      for( i = 0; i < nlinexprs; ++i ) /* expressions appearing linearly */
      {
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, linexprs[i], TRUE, FALSE, FALSE, FALSE) );
      }
      for( i = 0; i < nquadexprs; ++i ) /* expressions appearing quadratically */
      {
         SCIP_EXPR* quadexpr;
         SCIPexprGetQuadraticQuadTerm(expr, i, &quadexpr, NULL, NULL, NULL, NULL, NULL);
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, quadexpr, TRUE, FALSE, FALSE, FALSE) );
      }

      SCIPdebugMsg(scip, "expr %p is quadratic and propagable -> propagate and separate\n", (void*)expr);

      nlexprdata->separating = TRUE;
   }
   else
   {
      SCIPdebugMsg(scip, "expr %p is quadratic and propagable -> propagate only\n", (void*)expr);
   }

   if( SCIPexprAreQuadraticExprsVariables(expr) )
   {
      SCIPexprSetCurvature(expr, nlexprdata->curvature);
      SCIPdebugMsg(scip, "expr is %s in the original variables\n", nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE ? "concave" : "convex");
      nlexprdata->origvars = TRUE;
   }

   return SCIP_OKAY;
}

/** nonlinear handler auxiliary evaluation callback */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxQuadratic)
{  /*lint --e{715}*/
   int i;
   int nlinexprs;
   int nquadexprs;
   int nbilinexprs;
   SCIP_Real constant;
   SCIP_Real* lincoefs;
   SCIP_EXPR** linexprs;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(auxvalue != NULL);
   assert(nlhdlrexprdata->separating);
   assert(nlhdlrexprdata->qexpr == expr);

   /* if the quadratic is in the original variable we can just evaluate the expression */
   if( nlhdlrexprdata->origvars )
   {
      *auxvalue = SCIPexprGetEvalValue(expr);
      return SCIP_OKAY;
   }

   /* TODO there was a
     *auxvalue = SCIPevalExprQuadratic(scip, nlhdlrexprdata->qexpr, sol);
     here; any reason why not using this anymore?
   */

   SCIPexprGetQuadraticData(expr, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, &nbilinexprs, NULL, NULL);

   *auxvalue = constant;

   for( i = 0; i < nlinexprs; ++i ) /* linear exprs */
      *auxvalue += lincoefs[i] * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(linexprs[i]));

   for( i = 0; i < nquadexprs; ++i ) /* quadratic terms */
   {
      SCIP_Real solval;
      SCIP_Real lincoef;
      SCIP_Real sqrcoef;
      SCIP_EXPR* qexpr;

      SCIPexprGetQuadraticQuadTerm(expr, i, &qexpr, &lincoef, &sqrcoef, NULL, NULL, NULL);

      solval = SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(qexpr));
      *auxvalue += (lincoef + sqrcoef * solval) * solval;
   }

   for( i = 0; i < nbilinexprs; ++i ) /* bilinear terms */
   {
      SCIP_EXPR* expr1;
      SCIP_EXPR* expr2;
      SCIP_Real coef;

      SCIPexprGetQuadraticBilinTerm(expr, i, &expr1, &expr2, &coef, NULL, NULL);

      *auxvalue += coef * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr1)) * SCIPgetSolVal(scip, sol,
            SCIPgetExprAuxVarNonlinear(expr2));
   }

   return SCIP_OKAY;
}

/** nonlinear handler enforcement callback */
static
SCIP_DECL_NLHDLRENFO(nlhdlrEnfoQuadratic)
{  /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_ROWPREP* rowprep;
   SCIP_Bool success = FALSE;
   SCIP_NODE* node;
   int depth;
   SCIP_Longint nodenumber;
   SCIP_Real* eigenvalues;
   SCIP_Real violation;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->qexpr == expr);

   INTERLOG(printf("Starting interesection cuts!\n");)

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   /* estimate should take care of convex quadratics */
   if( ( overestimate && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONCAVE) ||
       (!overestimate && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX) )
   {
      INTERLOG(printf("Convex, no need of interesection cuts!\n");)
      return SCIP_OKAY;
   }

   /* nothing to do if we can't use intersection cuts */
   if( ! nlhdlrdata->useintersectioncuts )
   {
      INTERLOG(printf("We don't use intersection cuts!\n");)
      return SCIP_OKAY;
   }

   /* right now can use interesction cuts only if a basic LP solution is at hand; TODO: in principle we can do something
    * even if it is not optimal
    */
   if( sol != NULL || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL || !SCIPisLPSolBasic(scip) )
   {
      INTERLOG(printf("LP solutoin not good!\n");)
      return SCIP_OKAY;
   }

   /* only separate at selected nodes */
   node = SCIPgetCurrentNode(scip);
   depth = SCIPnodeGetDepth(node);
   if( (nlhdlrdata->atwhichnodes == -1 && depth != 0) || (nlhdlrdata->atwhichnodes != -1 && depth % nlhdlrdata->atwhichnodes != 0) )
   {
      INTERLOG(printf("Don't separate at this node\n");)
      return SCIP_OKAY;
   }

   /* do not add more than ncutslimitroot cuts in root node and ncutslimit cuts in the non-root nodes */
   nodenumber = SCIPnodeGetNumber(node);
   if( nlhdlrdata->lastnodenumber != nodenumber )
   {
      nlhdlrdata->lastnodenumber = nodenumber;
      nlhdlrdata->lastncuts = nlhdlrdata->ncutsadded;
   }
   /*else if( (depth > 0 && nlhdlrdata->ncutsadded - nlhdlrdata->lastncuts >= nlhdlrdata->ncutslimit) || (depth == 0 &&
            nlhdlrdata->ncutsadded - nlhdlrdata->lastncuts >= nlhdlrdata->ncutslimitroot)) */
   /* allow the addition of a certain number of cuts per quadratic */
   if( (depth > 0 && nlhdlrexprdata->ncutsadded >= nlhdlrdata->ncutslimit) || (depth == 0 &&
      nlhdlrexprdata->ncutsadded >= nlhdlrdata->ncutslimitroot) )
   {
      INTERLOG(printf("Too many cuts added already\n");)
      return SCIP_OKAY;
   }

   /* can't separate if we do not have an eigendecomposition */
   SCIPexprGetQuadraticData(expr, NULL, NULL, NULL, NULL, NULL, NULL, &eigenvalues, NULL);
   if( eigenvalues == NULL )
   {
      INTERLOG(printf("No known eigenvalues!\n");)
      return SCIP_OKAY;
   }

   /* if constraint is not sufficiently violated -> do nothing */
   if( cons != nlhdlrexprdata->cons )
   {
      /* constraint is w.r.t auxvar */
      violation = auxvalue - SCIPgetSolVal(scip, NULL, SCIPgetExprAuxVarNonlinear(expr));
      violation = ABS( violation );
   }
   else
      /* quadratic is a constraint */
      violation = MAX( SCIPgetLhsNonlinear(nlhdlrexprdata->cons) - auxvalue, auxvalue -
            SCIPgetRhsNonlinear(nlhdlrexprdata->cons)); /*lint !e666*/

   if( violation < nlhdlrdata->minviolation )
   {
      INTERLOG(printf("Violation %g is just too small\n", violation); )
      return SCIP_OKAY;
   }

   /* we can't build an intersection cut when the expr is the root of some constraint and also a subexpression of
    * another constraint because we initialize data differently TODO: how differently? */
   /* TODO: I don't think this is needed */
   if( nlhdlrexprdata->cons != NULL && cons != nlhdlrexprdata->cons )
   {
      INTERLOG(printf("WARNING!! expr is root of one constraint and subexpr of another!\n"); )
      return SCIP_OKAY;
   }

   /* if we are the root of a constraint and we are feasible w.r.t our auxiliary variables, that is, auxvalue is
    * actually feasible for the sides of the constraint, then do not separate
    */
   if( cons == nlhdlrexprdata->cons && ((overestimate && (SCIPgetLhsNonlinear(cons)) - auxvalue < SCIPfeastol(scip)) ||
            (! overestimate && (auxvalue - SCIPgetRhsNonlinear(cons) < SCIPfeastol(scip)))) )
   {
      INTERLOG(printf("We are actually feasible for the sides of the constraint\n"); )
      return SCIP_OKAY;
   }

#ifdef  DEBUG_INTERSECTIONCUT
   SCIPinfoMessage(scip, NULL, "Build intersection cut for \n");
   if( cons == nlhdlrexprdata->cons )
   {
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
   }
   else
   {
      SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, " == %s\n", SCIPvarGetName(SCIPgetExprAuxVarNonlinear(expr)));
   }
   SCIPinfoMessage(scip, NULL, "We need to %sestimate\n", overestimate ? "over" : "under" );
   SCIPinfoMessage(scip, NULL, "LP sol: \n");
   SCIP_CALL( SCIPprintTransSol(scip, NULL, NULL, FALSE) );
#endif
   *result = SCIP_DIDNOTFIND;

   /* cut (in the nonbasic space) is of the form alpha^T x >= 1 */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_LEFT, TRUE) );
   INTERLOG(printf("Generating inter cut\n"); )

   SCIP_CALL( generateIntercut(scip, expr, nlhdlrdata, nlhdlrexprdata, cons, sol, rowprep, overestimate, &success) );
   INTERLOG(if( !success) printf("Generation failed\n"); )

   /* we generated something, let us see if it survives the clean up */
   if( success )
   {
      assert(sol == NULL);
      nlhdlrdata->ncutsgenerated += 1;
      nlhdlrexprdata->ncutsadded += 1;

      /* merge coefficients that belong to same variable */
      SCIPmergeRowprepTerms(scip, rowprep);

      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, nlhdlrdata->mincutviolation, &violation, &success) );
      INTERLOG(if( !success) printf("Clean up failed\n"); )
   }

   /* if cut looks good (numerics ok and cutting off solution), then turn into row and add to sepastore */
   if( success )
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;

      /* count number of bound cuts */
      if( nlhdlrdata->useboundsasrays )
         nlhdlrdata->nboundcuts += 1;

      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "%s_intersection_quadratic%p_lp%" SCIP_LONGINT_FORMAT,
         overestimate ? "over" : "under",
         (void*)expr,
         SCIPgetNLPs(scip));

      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

      /*printf("## New cut\n");
      printf(" -> found maxquad-free cut <%s>: act=%f, lhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n\n",
            SCIProwGetName(row), SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row), SCIProwGetNorm(row),
            SCIPgetCutEfficacy(scip, NULL, row),
            SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
            SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row)); */

      /* intersection cuts can be numerically nasty; we do some extra numerical checks here */
      /*printf("SCIP DEPTH %d got a cut with violation %g, efficacy %g and r/e %g\n", SCIPgetSubscipDepth(scip),
       * violation, SCIPgetCutEfficacy(scip, NULL, row), SCIPgetRowMaxCoef(scip, row) / SCIPgetRowMinCoef(scip, row) /
       * SCIPgetCutEfficacy(scip, NULL, row));
       */
      assert(SCIPgetCutEfficacy(scip, NULL, row) > 0.0);
      if( ! nlhdlrdata->ignorehighre || SCIPgetRowMaxCoef(scip, row) / SCIPgetRowMinCoef(scip, row) / SCIPgetCutEfficacy(scip, NULL, row) < 1e9 )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "adding cut ");
         SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif

         SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
         }
         else
         {
            *result = SCIP_SEPARATED;
            nlhdlrdata->ncutsadded += 1;
            nlhdlrdata->densitysum += (SCIP_Real) SCIProwprepGetNVars(rowprep) / (SCIP_Real) SCIPgetNVars(scip);
         }
      }
      else
      {
         nlhdlrdata->nhighre++;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** nonlinear handler forward propagation callback
 *
 * This method should solve the problem
 * <pre>
 *    max/min quad expression over box constraints
 * </pre>
 * However, this problem is difficult so we are satisfied with a proxy.
 * Interval arithmetic suffices when no variable appears twice, however this is seldom the case, so we try
 * to take care of the dependency problem to some extent:
 * Let \f$P_l = \{i : \text{expr}_l \text{expr}_i \,\text{is a bilinear expr}\}\f$.
 * 1. partition the quadratic expression as sum of quadratic functions \f$\sum_l q_l\f$
 *    where \f$q_l = a_l \text{expr}_l^2 + c_l \text{expr}_l + \sum_{i \in P_l} b_{il} \text{expr}_i \text{expr}_l\f$
 * 2. build interval quadratic functions, i.e., \f$a x^2 + b x\f$ where \f$b\f$ is an interval, i.e.,
 *    \f$a_l \text{expr}_l^2 + [\sum_{i \in P_l} b_{il} \text{expr}_i + c_l] \text{expr}_l\f$
 * 3. compute \f$\min/\max \{ a x^2 + b x : x \in [x] \}\f$ for each interval quadratic, i.e.,
 *    \f$\min/\max a_l \text{expr}_l^2 + \text{expr}_l [\sum_{i \in P_l} b_{il} \text{expr}_i + c_l] : \text{expr}_l \in [\text{expr}_l]\f$
 *
 * Notes:
 * 1. The \f$l\f$-th quadratic expr (expressions that appear quadratically) is associated with \f$q_l\f$.
 * 2. `nlhdlrdata->quadactivities[l]` is the activity of \f$q_l\f$ as computed in the description above.
 * 3. The \f$q_l\f$ of a quadratic term might be empty, in which case `nlhdlrdata->quadactivities[l]` is [0,0].\n
 *    For example, consider \f$x^2 + xy\f$. There are two quadratic expressions, \f$x\f$ and \f$y\f$.
 *    The \f$q\f$ associated to \f$x\f$ is \f$x^2 + xy\f$, while the \f$q\f$ associated to \f$y\f$ is empty.
 *    Thus, `nlhdlrdata->quadactivities[1]` is [0,0] in this case.
 *    The logic is to avoid considering the term \f$xy\f$ twice.
 *
 * @note The order matters! If \f$\text{expr}_i\, \text{expr}_l\f$ is a term in the quadratic, then \f$i\f$ is *not* in \f$P_l\f$
 */
static
SCIP_DECL_NLHDLRINTEVAL(nlhdlrIntevalQuadratic)
{ /*lint --e{715}*/
   SCIP_EXPR** linexprs;
   SCIP_Real* lincoefs;
   SCIP_Real constant;
   int nquadexprs;
   int nlinexprs;

   assert(scip != NULL);
   assert(expr != NULL);

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->quadactivities != NULL);
   assert(nlhdlrexprdata->qexpr == expr);

   SCIPdebugMsg(scip, "Interval evaluation of quadratic expr\n");

   SCIPexprGetQuadraticData(expr, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, NULL, NULL, NULL);

   /*
    * compute activity of linear part, if some linear term has changed
    */
   {
      int i;

      SCIPdebugMsg(scip, "Computing activity of linear part\n");

      SCIPintervalSet(&nlhdlrexprdata->linactivity, constant);
      for( i = 0; i < nlinexprs; ++i )
      {
         SCIP_INTERVAL linterminterval;

         linterminterval = SCIPexprGetActivity(linexprs[i]);
         if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, linterminterval) )
         {
            SCIPdebugMsg(scip, "Activity of linear part is empty due to child %d\n", i);
            SCIPintervalSetEmpty(interval);
            return SCIP_OKAY;
         }
         SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &linterminterval, linterminterval, lincoefs[i]);
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &nlhdlrexprdata->linactivity, nlhdlrexprdata->linactivity, linterminterval);
      }

      SCIPdebugMsg(scip, "Activity of linear part is [%g, %g]\n", nlhdlrexprdata->linactivity.inf,
            nlhdlrexprdata->linactivity.sup);
   }

   /*
    * compute activity of quadratic part
    */
   {
      int i;

      SCIPdebugMsg(scip, "Computing activity of quadratic part\n");

      nlhdlrexprdata->nneginfinityquadact = 0;
      nlhdlrexprdata->nposinfinityquadact = 0;
      nlhdlrexprdata->minquadfiniteact = 0.0;
      nlhdlrexprdata->maxquadfiniteact = 0.0;
      SCIPintervalSet(&nlhdlrexprdata->quadactivity, 0.0);

      for( i = 0; i < nquadexprs; ++i )
      {
         SCIP_Real quadlb;
         SCIP_Real quadub;
         SCIP_EXPR* qexpr;
         SCIP_Real lincoef;
         SCIP_Real sqrcoef;
         int nadjbilin;
         int* adjbilin;
         SCIP_EXPR* sqrexpr;

         SCIPexprGetQuadraticQuadTerm(expr, i, &qexpr, &lincoef, &sqrcoef, &nadjbilin, &adjbilin, &sqrexpr);

         if( !isPropagableTerm(expr, i) )
         {
            /* term is not propagable, i.e., the exprs involved in term only appear once; thus use the activity of the
             * quadratic term directly and not the activity of the exprs involed in the term. See also documentation of
             * DETECT
             */
            SCIP_INTERVAL tmp;

            assert(lincoef == 0.0);

            if( sqrcoef != 0.0 )
            {
               assert(sqrexpr != NULL);
               assert(nadjbilin == 0);

               tmp = SCIPexprGetActivity(sqrexpr);
               if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, tmp) )
               {
                  SCIPintervalSetEmpty(interval);
                  return SCIP_OKAY;
               }

               SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &tmp, tmp, sqrcoef);
               quadlb = tmp.inf;
               quadub = tmp.sup;

#ifdef DEBUG_PROP
               SCIPinfoMessage(scip, NULL, "Computing activity for quadratic term %g <expr>, where <expr> is: ", sqrcoef);
               SCIP_CALL( SCIPprintExpr(scip, sqrexpr, NULL) );
#endif
            }
            else
            {
               SCIP_EXPR* expr1;
               SCIP_EXPR* prodexpr;
               SCIP_Real prodcoef;

               assert(nadjbilin == 1);
               SCIPexprGetQuadraticBilinTerm(expr, adjbilin[0], &expr1, NULL, &prodcoef, NULL, &prodexpr);

               if( expr1 == qexpr )
               {
                  /* the quadratic expression expr1 appears only as expr1 * expr2, so its 'q' is expr1 * expr2 */
                  tmp = SCIPexprGetActivity(prodexpr);
                  if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, tmp) )
                  {
                     SCIPintervalSetEmpty(interval);
                     return SCIP_OKAY;
                  }

                  SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &tmp, tmp, prodcoef);
                  quadlb = tmp.inf;
                  quadub = tmp.sup;

#ifdef DEBUG_PROP
                  SCIPinfoMessage(scip, NULL, "Computing activity for quadratic term %g <expr>, where <expr> is: ", prodcoef);
                  SCIP_CALL( SCIPprintExpr(scip, prodexpr, NULL) );
#endif
               }
               else
               {
                  /* the quadratic expression expr1 appears as expr2 * expr1, thus its 'q' is empty, see also the Notes
                   * in the documentation of the function
                   */
                  SCIPintervalSet(&nlhdlrexprdata->quadactivities[i], 0.0);
                  continue;
               }
            }
         }
         else
         {
            int j;
            SCIP_INTERVAL b;

            SCIPexprGetQuadraticQuadTerm(expr, i, &qexpr, &lincoef, &sqrcoef, &nadjbilin, &adjbilin, NULL);

            if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(qexpr)) )
            {
               SCIPintervalSetEmpty(interval);
               return SCIP_OKAY;
            }

            /* b = [c_l] */
            SCIPintervalSet(&b, lincoef);
#ifdef DEBUG_PROP
            SCIPinfoMessage(scip, NULL, "b := %g\n", lincoef);
#endif
            for( j = 0; j < nadjbilin; ++j )
            {
               SCIP_INTERVAL bterm;
               SCIP_EXPR* expr1;
               SCIP_EXPR* expr2;
               SCIP_Real bilincoef;

               SCIPexprGetQuadraticBilinTerm(expr, adjbilin[j], &expr1, &expr2, &bilincoef, NULL, NULL);

               if( expr1 != qexpr )
                  continue;

               bterm = SCIPexprGetActivity(expr2);
               if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bterm) )
               {
                  SCIPintervalSetEmpty(interval);
                  return SCIP_OKAY;
               }

               /* b += [b_jl * expr_j] for j \in P_l */
               SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &bterm, bterm, bilincoef);
               SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &b, b, bterm);

#ifdef DEBUG_PROP
               SCIPinfoMessage(scip, NULL, "b += %g * [expr2], where <expr2> is: ", bilincoef);
               SCIP_CALL( SCIPprintExpr(scip, expr2, NULL) );
               SCIPinfoMessage(scip, NULL, " [%g,%g]\n", SCIPexprGetActivity(expr2).inf, SCIPexprGetActivity(expr2).sup);
#endif
            }

            /* TODO: under which assumptions do we know that we just need to compute min or max? its probably the locks that give some information here */
            quadub = SCIPintervalQuadUpperBound(SCIP_INTERVAL_INFINITY, sqrcoef, b,
               SCIPexprGetActivity(qexpr));

            /* TODO: implement SCIPintervalQuadLowerBound */
            {
               SCIP_INTERVAL minusb;
               SCIPintervalSetBounds(&minusb, -SCIPintervalGetSup(b), -SCIPintervalGetInf(b));

               quadlb = -SCIPintervalQuadUpperBound(SCIP_INTERVAL_INFINITY, -sqrcoef, minusb,
                  SCIPexprGetActivity(qexpr));
            }

#ifdef DEBUG_PROP
            SCIPinfoMessage(scip, NULL, "Computing activity for quadratic term %g <expr>^2 + [%g,%g] <expr>, where <expr> is: ", sqrcoef, b.inf, b.sup);
            SCIP_CALL( SCIPprintExpr(scip, qexpr, NULL) );
#endif
         }
#ifdef DEBUG_PROP
         SCIPinfoMessage(scip, NULL, " -> [%g, %g]\n", quadlb, quadub);
#endif

         SCIPintervalSetBounds(&nlhdlrexprdata->quadactivities[i], quadlb, quadub);
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &nlhdlrexprdata->quadactivity, nlhdlrexprdata->quadactivity, nlhdlrexprdata->quadactivities[i]);

         /* get number of +/-infinity contributions and compute finite activity */
         if( quadlb <= -SCIP_INTERVAL_INFINITY )
            nlhdlrexprdata->nneginfinityquadact++;
         else
         {
            SCIP_ROUNDMODE roundmode;

            roundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeDownwards();

            nlhdlrexprdata->minquadfiniteact += quadlb;

            SCIPintervalSetRoundingMode(roundmode);
         }
         if( quadub >= SCIP_INTERVAL_INFINITY )
            nlhdlrexprdata->nposinfinityquadact++;
         else
         {
            SCIP_ROUNDMODE roundmode;

            roundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeUpwards();

            nlhdlrexprdata->maxquadfiniteact += quadub;

            SCIPintervalSetRoundingMode(roundmode);
         }
      }

      SCIPdebugMsg(scip, "Activity of quadratic part is [%g, %g]\n", nlhdlrexprdata->quadactivity.inf, nlhdlrexprdata->quadactivity.sup);
   }

   /* interval evaluation is linear activity + quadactivity */
   SCIPintervalAdd(SCIP_INTERVAL_INFINITY, interval, nlhdlrexprdata->linactivity,  nlhdlrexprdata->quadactivity);

   nlhdlrexprdata->activitiestag = SCIPgetCurBoundsTagNonlinear(SCIPfindConshdlr(scip, "nonlinear"));

   return SCIP_OKAY;
}

/** nonlinear handler reverse propagation callback
 *
 * @note the implemented technique is a proxy for solving the problem min/max{ x_i : quad expr in [quad expr] }
 * and as such can be improved.
 */
static
SCIP_DECL_NLHDLRREVERSEPROP(nlhdlrReversepropQuadratic)
{ /*lint --e{715}*/
   SCIP_EXPR** linexprs;
   SCIP_EXPR** bilinexprs; /* TODO: should this be stored in the nlhdlr expr data? */
   SCIP_Real* bilincoefs;
   SCIP_Real* lincoefs;
   SCIP_Real constant;
   int nquadexprs;
   int nlinexprs;

   SCIP_INTERVAL rhs;
   SCIP_INTERVAL quadactivity;
   int i;

   SCIPdebugMsg(scip, "Reverse propagation of quadratic expr given bounds = [%g,%g]\n", bounds.inf, bounds.sup);

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->quadactivities != NULL);
   assert(nlhdlrexprdata->qexpr == expr);

   *nreductions = 0;

   /* not possible to conclude finite bounds if the interval of the expression is [-inf,inf] */
   if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, bounds) )
   {
      SCIPdebugMsg(scip, "expr's range is R -> cannot reverse propagate\n");
      return SCIP_OKAY;
   }

   /* ensure that partial activities as stored in nlhdlrexprdata are uptodate
    * if the activity stored in expr is more recent than the partial activities stored in this nlhdlrexprdata,
    * then we should reevaluate the partial activities
    */
   if( SCIPexprGetActivityTag(expr) > nlhdlrexprdata->activitiestag )
   {
      SCIP_CALL( nlhdlrIntevalQuadratic(scip, nlhdlr, expr, nlhdlrexprdata, &quadactivity, NULL, NULL) );
   }

   SCIPexprGetQuadraticData(expr, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, NULL, NULL, NULL);

   /* propagate linear part in rhs = expr's interval - quadratic activity; first, reconstruct the quadratic activity */
   SCIPintervalSetBounds(&quadactivity,
         nlhdlrexprdata->nneginfinityquadact > 0 ? -SCIP_INTERVAL_INFINITY : nlhdlrexprdata->minquadfiniteact,
         nlhdlrexprdata->nposinfinityquadact > 0 ?  SCIP_INTERVAL_INFINITY : nlhdlrexprdata->maxquadfiniteact);

   SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs, bounds, quadactivity);

   SCIP_CALL( reversePropagateLinearExpr(scip, linexprs, nlinexprs, lincoefs, constant, rhs, infeasible, nreductions) );

   /* stop if we find infeasibility */
   if( *infeasible )
      return SCIP_OKAY;

   /* propagate quadratic part in expr's interval - linear activity, where linear activity was computed in INTEVAL.
    * The idea is basically to write interval quadratics for each expr and then solve for expr.
    *
    * One way of achieving this is:
    * - for each expression expr_i, write the quadratic expression as a_i expr^2_i + expr_i ( \sum_{j \in J_i} b_ij
    *   expr_j + c_i ) + quadratic expression in expr_k for k \neq i
    * - compute the interval b = [\sum_{j \in J_i} b_ij expr_j + c_i], where J_i are all the indices j such that the
    *   bilinear expression expr_i expr_j appears
    * - use some technique (like the one in nlhdlrIntevalQuadratic), to evaluate the activity of rest_i = [quadratic
    *   expression in expr_k for k \neq i].
    * - solve a_i expr_i^2 + b expr_i \in rhs_i := [expr activity] - rest_i
    *
    * However, this might be expensive, especially computing rest_i. Hence, we implement a simpler version.
    * - we use the same partition as in nlhdlrIntevalQuadratic for the bilinear terms. This way, b = [\sum_{j \in P_i}
    *   b_ij expr_j + c_i], where P_i is the set of indices j such that expr_i * expr_j appears in that order
    * - we evaluate the activity of rest_i as sum_{k \neq i} [\min q_k, \max q_k] where q_k = a_k expr_k^2 + [\sum_{j
    *   \in P_k} b_jk expr_j + c_k] expr_k. The intervals [\min q_k, \max q_k] were already computed in
    *   nlhdlrIntevalQuadratic, so we just reuse them.
    *
    * A downside of the above is that we might not deduce any bounds for variables that appear less often. For example,
    * consider x^2 + x * y + x * z + y * z + z. This quadratic gets partitioned as (x^2 + x*y + x*z) + (z*y + z). The
    * first parenthesis is interpreted as a function of x, while the second one as a function of z.
    * To also get bounds on y, after reverse propagating x in x^2 + x*y + x*z \in rhs, we rewrite this as y + z \in rhs/x -
    * x and propagate the y + z).
    * In general, after reverse propagating expr_i, we consider
    *   \sum_{j \in J_i} b_ij expr_j in ([expr activity] - quadratic expression in expr_k for k \neq i - c_i) / expr_i - a_i expr_i,
    * compute an interval for the right hand side (see computeRangeForBilinearProp) and use that to propagate the
    * linear sum on the left hand side.
    *
    * Note: this last step generalizes a technique that appeared in the classic cons_quadratic.
    * The idea of that technique was to borrow a bilinear term expr_k expr_l when propagating expr_l and the quadratic
    * function for expr_k was simple enough.
    * Since in P_l we only consider the indices of expressions that appear multiplying expr_l as _second_ factor, we
    * would lose the bilinear terms expr_k * expr_l, which contributes to the dependency problem.
    * The problem is that the contribution of b_kl * expr_k * expr_l to rest_i is not just [b_kl * expr_k * expr_l], but
    * rather quadactivities[k] (= max/min of a_k expr_k^2 + expr_k * [c_k + sum_i \in P_k b_ki expr_i]).
    * Thus, we _cannot_ just substract [b_kl * expr_k * expr_l] from rest_i.
    * But, if expr_k only appears as expr_k * expr_l, then  quadactivities[k] = [b_kl * expr_k * expr_l]. So this
    * case was handled in old cons_quadratic.
    *
    *
    * TODO: handle simple cases
    * TODO: identify early when there is nothing to be gain
    */
   SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs, bounds, nlhdlrexprdata->linactivity);
   SCIP_CALL( SCIPallocBufferArray(scip, &bilinexprs, nquadexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bilincoefs, nquadexprs) );

   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_INTERVAL rhs_i;
      SCIP_INTERVAL rest_i;
      SCIP_EXPR* qexpr;
      SCIP_Real lincoef;
      SCIP_Real sqrcoef;
      int nadjbilin;
      int* adjbilin;
      SCIP_EXPR* sqrexpr;

      SCIPexprGetQuadraticQuadTerm(expr, i, &qexpr, &lincoef, &sqrcoef, &nadjbilin, &adjbilin, &sqrexpr);

      /* rhs_i = rhs - rest_i.
       * to compute rest_i = [\sum_{k \neq i} q_k] we just have to substract
       * the activity of q_i from quadactivity; however, care must be taken about infinities;
       * if [q_i].sup = +infinity and there is = 1 contributing +infinity -> rest_i.sup = maxquadfiniteact
       * if [q_i].sup = +infinity and there is > 1 contributing +infinity -> rest_i.sup = +infinity
       * if [q_i].sup = finite and there is > 0 contributing +infinity -> rest_i.sup = +infinity
       * if [q_i].sup = finite and there is = 0 contributing +infinity -> rest_i.sup = maxquadfiniteact - [q_i].sup
       *
       * the same holds when replacing sup with inf, + with - and max(quadfiniteact) with min(...)
       */
      /* compute rest_i.sup */
      if( SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]) < SCIP_INTERVAL_INFINITY &&
         nlhdlrexprdata->nposinfinityquadact == 0 )
      {
         SCIP_ROUNDMODE roundmode;

         roundmode = SCIPintervalGetRoundingMode();
         SCIPintervalSetRoundingModeUpwards();
         rest_i.sup = nlhdlrexprdata->maxquadfiniteact - SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]);

         SCIPintervalSetRoundingMode(roundmode);
      }
      else if( SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]) >= SCIP_INTERVAL_INFINITY &&
         nlhdlrexprdata->nposinfinityquadact == 1 )
         rest_i.sup = nlhdlrexprdata->maxquadfiniteact;
      else
         rest_i.sup = SCIP_INTERVAL_INFINITY;

      /* compute rest_i.inf */
      if( SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]) > -SCIP_INTERVAL_INFINITY &&
         nlhdlrexprdata->nneginfinityquadact == 0 )
      {
         SCIP_ROUNDMODE roundmode;

         roundmode = SCIPintervalGetRoundingMode();
         SCIPintervalSetRoundingModeDownwards();
         rest_i.inf = nlhdlrexprdata->minquadfiniteact - SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]);

         SCIPintervalSetRoundingMode(roundmode);
      }
      else if( SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]) <= -SCIP_INTERVAL_INFINITY &&
         nlhdlrexprdata->nneginfinityquadact == 1 )
         rest_i.inf = nlhdlrexprdata->minquadfiniteact;
      else
         rest_i.inf = -SCIP_INTERVAL_INFINITY;

#ifdef SCIP_DISABLED_CODE  /* I (SV) added the following in cons_quadratic to fix/workaround some bug. Maybe we'll need this here, too? */
      /* FIXME in theory, rest_i should not be empty here
       * what we tried to do here is to remove the contribution of the i'th bilinear term (=bilinterm) to [minquadactivity,maxquadactivity] from rhs
       * however, quadactivity is computed differently (as x*(a1*y1+...+an*yn)) than q_i (a*ak*yk) and since interval arithmetics do overestimation,
       * it can happen that q_i is actually slightly larger than quadactivity, which results in rest_i being (slightly) empty
       * a proper fix could be to compute the quadactivity also as x*a1*y1+...+x*an*yn if sqrcoef=0, but due to taking
       * also infinite bounds into account, this complicates the code even further
       * instead, I'll just work around this by turning an empty rest_i into a small non-empty one
       */
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, rest_i) )
      {
         assert(SCIPisSumRelEQ(scip, rest_i.inf, rest_i.sup));
         SCIPswapReals(&rest_i.inf, &rest_i.sup);
      }
#endif
      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, rest_i));

      /* compute rhs_i */
      SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs_i, rhs, rest_i);

      if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, rhs_i) )
         continue;

      /* try to propagate */
      if( !isPropagableTerm(expr, i) )
      {
         assert(lincoef == 0.0);

         if( sqrcoef != 0.0 )
         {
            assert(sqrexpr != NULL);
            assert(nadjbilin == 0);

            /* solve sqrcoef sqrexpr in rhs_i */
            SCIP_CALL( propagateBoundsLinExpr(scip, sqrexpr, sqrcoef, rhs_i, infeasible, nreductions) );
         }
         else
         {
            /* qexpr only appears in a term of the form qexpr * other_expr (or other_expr * qexpr); we only care about
             * getting bounds for the product, thus we will compute these bounds when qexpr appears as qexpr *
             * other_expr; note that if it appears as other_expr * qexpr, then when we process other_expr bounds for the
             * product will be computed
             * TODO: we can actually avoid computing rhs_i in the case that qexpr is not propagable and it appears as
             * other_expr * qexpr
             */
            SCIP_EXPR* expr1;
            SCIP_EXPR* prodexpr;
            SCIP_Real prodcoef;

            assert(nadjbilin == 1);
            SCIPexprGetQuadraticBilinTerm(expr, adjbilin[0], &expr1, NULL, &prodcoef, NULL, &prodexpr);

            if( expr1 == qexpr )
            {
               /* solve prodcoef prodexpr in rhs_i */
               SCIP_CALL( propagateBoundsLinExpr(scip, prodexpr, prodcoef, rhs_i, infeasible, nreductions) );
            }
         }
      }
      else
      {
         SCIP_INTERVAL b;
         SCIP_EXPR* expr1 = NULL;
         SCIP_EXPR* expr2 = NULL;
         SCIP_Real bilincoef = 0.0;
         int nbilin = 0;
         int pos2 = 0;
         int j;

         /* set b to [c_l] */
         SCIPintervalSet(&b, lincoef);

         /* add [\sum_{j \in P_l} b_lj expr_j + c_l] into b */
         for( j = 0; j < nadjbilin; ++j )
         {
            SCIP_INTERVAL bterm;
            SCIP_INTERVAL expr2bounds;

            SCIPexprGetQuadraticBilinTerm(expr, adjbilin[j], &expr1, &expr2, &bilincoef, &pos2, NULL);

            if( expr1 != qexpr )
               continue;

            expr2bounds = SCIPgetExprBoundsNonlinear(scip, expr2);
            if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, expr2bounds) )
            {
               *infeasible = TRUE;
               break;
            }

            /* b += [b_lj * expr_j] for j \in P_l */
            SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &bterm, expr2bounds, bilincoef);
            SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &b, b, bterm);

            /* remember b_lj and expr_j to propagate them too */
            bilinexprs[nbilin] = expr2;
            bilincoefs[nbilin] = bilincoef;
            nbilin++;
         }

         if( !*infeasible )
         {
            /* solve a_i expr_i^2 + b expr_i in rhs_i */
            SCIP_CALL( propagateBoundsQuadExpr(scip, qexpr, sqrcoef, b, rhs_i, infeasible, nreductions) );
         }

         if( nbilin > 0 && !*infeasible )
         {
            /* if 0 is not in [expr_i], then propagate bilincoefs^T bilinexpr in rhs_i/expr_i - a_i expr_i - c_i */
            SCIP_INTERVAL bilinrhs;
            SCIP_INTERVAL qexprbounds;

            qexprbounds = SCIPgetExprBoundsNonlinear(scip, qexpr);
            if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, qexprbounds) )
            {
               *infeasible = TRUE;
            }
            else
            {
               /* compute bilinrhs := [rhs_i/expr_i - a_i expr_i] */
               computeRangeForBilinearProp(qexprbounds, sqrcoef, rhs_i, &bilinrhs);

               if( !SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, bilinrhs) )
               {
                  int nreds;

                  /* propagate \sum_{j \in P_i} b_ij expr_j + c_i in bilinrhs */
                  SCIP_CALL( reversePropagateLinearExpr(scip, bilinexprs, nbilin, bilincoefs, lincoef, bilinrhs,
                           infeasible, &nreds) );

                  /* TODO FIXME: we are overestimating the number of reductions: an expr might be tightened many times! */
                  *nreductions += nreds;
               }
            }
         }
      }

      /* stop if we find infeasibility */
      if( *infeasible )
         break;
   }

   SCIPfreeBufferArray(scip, &bilincoefs);
   SCIPfreeBufferArray(scip, &bilinexprs);

   return SCIP_OKAY;
}

/** callback to free data of handler */
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataQuadratic)
{ /*lint --e{715}*/
   assert(nlhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrQuadratic)
{  /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrQuadratic(targetscip) );

   return SCIP_OKAY;
}

/** includes quadratic nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   /* create nonlinear handler specific data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   BMSclearMemory(nlhdlrdata);

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectQuadratic, nlhdlrEvalauxQuadratic, nlhdlrdata) );

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrQuadratic);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrdataQuadratic);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeexprdataQuadratic);
   SCIPnlhdlrSetSepa(nlhdlr, NULL, nlhdlrEnfoQuadratic, NULL, NULL);
   SCIPnlhdlrSetProp(nlhdlr, nlhdlrIntevalQuadratic, nlhdlrReversepropQuadratic);

   /* parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/useintersectioncuts",
         "whether to use intersection cuts for quadratic constraints to separate",
         &nlhdlrdata->useintersectioncuts, FALSE, DEFAULT_USEINTERCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/usestrengthening",
         "whether the strengthening should be used",
         &nlhdlrdata->usestrengthening, FALSE, DEFAULT_USESTRENGTH, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/useboundsasrays",
         "use bounds of variables in quadratic as rays for intersection cuts",
         &nlhdlrdata->useboundsasrays, FALSE, DEFAULT_USEBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/ncutslimit",
         "limit for number of cuts generated consecutively",
         &nlhdlrdata->ncutslimit, FALSE, DEFAULT_NCUTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/ncutslimitroot",
         "limit for number of cuts generated at root node",
         &nlhdlrdata->ncutslimitroot, FALSE, DEFAULT_NCUTSROOT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/maxrank",
         "maximal rank a slackvar can have",
         &nlhdlrdata->maxrank, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "nlhdlr/" NLHDLR_NAME "/mincutviolation",
         "minimal cut violation the generated cuts must fulfill to be added to the LP",
         &nlhdlrdata->mincutviolation, FALSE, 1e-4, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "nlhdlr/" NLHDLR_NAME "/minviolation",
         "minimal violation the constraint must fulfill such that a cut is generated",
         &nlhdlrdata->mincutviolation, FALSE, INTERCUTS_MINVIOL, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/atwhichnodes",
         "determines at which nodes cut is used (if it's -1, it's used only at the root node, if it's n >= 0, it's used at every multiple of n",
         &nlhdlrdata->atwhichnodes, FALSE, 1, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/nstrengthlimit",
         "limit for number of rays we do the strengthening for",
         &nlhdlrdata->nstrengthlimit, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/ignorebadrayrestriction",
         "should cut be generated even with bad numerics when restricting to ray?",
         &nlhdlrdata->ignorebadrayrestriction, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/ignorenhighre",
         "should cut be added even when range / efficacy is large?",
         &nlhdlrdata->ignorehighre, FALSE, TRUE, NULL, NULL) );

   /* statistic table */
   assert(SCIPfindTable(scip, TABLE_NAME_QUADRATIC) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_QUADRATIC, TABLE_DESC_QUADRATIC, FALSE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputQuadratic,
         NULL, TABLE_POSITION_QUADRATIC, TABLE_EARLIEST_STAGE_QUADRATIC) );
   return SCIP_OKAY;
}
