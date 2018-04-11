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

/**@file   sepa_disjunctive.c
 * @brief  disjunctive cut separator
 * @author Tobias Fischer
 * @author Marc Pfetsch
 *
 * We separate disjunctive cuts for two term disjunctions of the form \f$x_1 = 0 \vee x_2 = 0\f$. They can be generated
 * directly from the simplex tableau. For further information, we refer to@n
 * "A complementarity-based partitioning and disjunctive cut algorithm for mathematical programming problems with
 * equilibrium constraints"@n
 * Júdice, J.J., Sherali, H.D., Ribeiro, I.M., Faustino, A.M., Journal of Global Optimization 36(1), 89–114 (2006)
 *
 * Cut coefficients belonging to integer variables can be strengthened by the 'monoidal cut strengthening' procedure, see@n
 * "Strengthening cuts for mixed integer programs"@n
 * Egon Balas, Robert G. Jeroslow, European Journal of Operational Research, Volume 4, Issue 4, 1980, Pages 224-234
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "scip/sepa_disjunctive.h"
#include "scip/cons_sos1.h"


#define SEPA_NAME              "disjunctive"
#define SEPA_DESC              "disjunctive cut separator"
#define SEPA_PRIORITY                10 /**< priority for separation */
#define SEPA_FREQ                     0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define SEPA_MAXBOUNDDIST           0.0 /**< maximal relative distance from the current node's dual bound to primal bound
                                         *   compared to best node's dual bound for applying separation.*/
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                 TRUE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXRANK              20 /**< maximal rank of a cut that could not be scaled to integral coefficients (-1: unlimited) */
#define DEFAULT_MAXRANKINTEGRAL      -1 /**< maximal rank of a cut that could be scaled to integral coefficients (-1: unlimited) */
#define DEFAULT_MAXWEIGHTRANGE    1e+03 /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
#define DEFAULT_STRENGTHEN         TRUE /**< strengthen cut if integer variables are present */

#define DEFAULT_MAXDEPTH             -1 /**< node depth of separating cuts (-1: no limit) */
#define DEFAULT_MAXROUNDS            25 /**< maximal number of separation rounds in a branching node (-1: no limit) */
#define DEFAULT_MAXROUNDSROOT       100 /**< maximal number of separation rounds in the root node (-1: no limit) */
#define DEFAULT_MAXINVCUTS           50 /**< maximal number of cuts investigated per iteration in a branching node */
#define DEFAULT_MAXINVCUTSROOT      250 /**< maximal number of cuts investigated per iteration in the root node */
#define DEFAULT_MAXCONFSDELAY    100000 /**< delay separation if number of conflict graph edges is larger than predefined value (-1: no limit) */

#define MAKECONTINTEGRAL          FALSE /**< convert continuous variable to integral variables in SCIPmakeRowIntegral() */


/** separator data */
struct SCIP_SepaData
{
   SCIP_Bool             strengthen;         /**< strengthen cut if integer variables are present */
   SCIP_CONSHDLR*        conshdlr;           /**< SOS1 constraint handler */
   SCIP_Real             maxweightrange;     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   int                   maxrank;            /**< maximal rank of a cut that could not be scaled to integral coefficients (-1: unlimited) */
   int                   maxrankintegral;    /**< maximal rank of a cut that could be scaled to integral coefficients (-1: unlimited) */
   int                   maxdepth;           /**< node depth of separating cuts (-1: no limit) */
   int                   maxrounds;          /**< maximal number of separation rounds in a branching node (-1: no limit) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: no limit) */
   int                   maxinvcuts;         /**< maximal number of cuts separated per iteration in a branching node */
   int                   maxinvcutsroot;     /**< maximal number of cuts separated per iteration in the root node */
   int                   maxconfsdelay;      /**< delay separation if number of conflict graph edges is larger than predefined value (-1: no limit) */
   int                   lastncutsfound;     /**< total number of cuts found after last call of separator */
};


/** gets rank of variable corresponding to row of \f$B^{-1}\f$ */
static
int getVarRank(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_Real*            binvrow,            /**< row of \f$B^{-1}\f$ */
   SCIP_Real*            rowsmaxval,         /**< maximal row multiplicator from nonbasic matrix A_N */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_ROW**            rows,               /**< rows of LP relaxation of scip */
   int                   nrows               /**< number of rows */
   )
{
   SCIP_Real maxweight = 0.0;
   int maxrank = 0;
   int r;

   assert( scip != NULL );
   assert( binvrow != NULL || nrows == 0 );
   assert( rowsmaxval != NULL || nrows == 0 );
   assert( rows != NULL || nrows == 0 );

   /* compute maximum row weights resulting from multiplication */
   for (r = 0; r < nrows; ++r)
   {
      SCIP_Real val;

      val = REALABS(binvrow[r] * rowsmaxval[r]);/*lint !e613*/
      if ( SCIPisGT(scip, val, maxweight) )
         maxweight = val;
   }

   /* compute rank */
   for (r = 0; r < nrows; ++r)
   {
      SCIP_Real val;
      int rank;

      val = REALABS(binvrow[r] * rowsmaxval[r]);/*lint !e613*/
      rank = SCIProwGetRank(rows[r]);/*lint !e613*/
      if ( rank > maxrank && SCIPisGT(scip, val * maxweightrange, maxweight) )
         maxrank = rank;
   }

   return maxrank;
}


/** gets the nonbasic coefficients of a simplex row */
static
SCIP_RETCODE getSimplexCoefficients(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_ROW**            rows,               /**< LP rows */
   int                   nrows,              /**< number LP rows */
   SCIP_COL**            cols,               /**< LP columns */
   int                   ncols,              /**< number of LP columns */
   SCIP_Real*            coef,               /**< row of \f$B^{-1} \cdot A\f$ */
   SCIP_Real*            binvrow,            /**< row of \f$B^{-1}\f$ */
   SCIP_Real*            simplexcoefs,       /**< pointer to store the nonbasic simplex-coefficients */
   int*                  nonbasicnumber      /**< pointer to store the number of nonbasic simplex-coefficients */
   )
{
   int r;
   int c;

   assert( scip != NULL );
   assert( rows != NULL );
   assert( nonbasicnumber != NULL );
   assert( simplexcoefs != NULL );
   assert( cols != NULL );

   *nonbasicnumber = 0;

   /* note: the non-slack variables have to be added first (see the function generateDisjCutSOS1()) */

   /* get simplex-coefficients of the non-basic non-slack variables */
   for (c = 0; c < ncols; ++c)
   {
      SCIP_COL* col;

      col = cols[c];
      assert( col != NULL );
      if ( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER  || SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER )
         simplexcoefs[(*nonbasicnumber)++] = coef[c];
   }

   /* get simplex-coefficients of the non-basic slack variables */
   for (r = 0; r < nrows; ++r)
   {
      SCIP_ROW* row;
      row = rows[r];
      assert( row != NULL );

      if ( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER || SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER )
      {
         assert( SCIPisFeasZero(scip, SCIPgetRowLPActivity(scip, row) - SCIProwGetRhs(row)) || SCIPisFeasZero(scip, SCIPgetRowLPActivity(scip, row) - SCIProwGetLhs(row)) );

         simplexcoefs[(*nonbasicnumber)++] = binvrow[r];
      }
   }

   return SCIP_OKAY;
}


/** computes a disjunctive cut inequality based on two simplex taubleau rows */
static
SCIP_RETCODE generateDisjCutSOS1(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_ROW**            rows,               /**< LP rows */
   int                   nrows,              /**< number of LP rows */
   SCIP_COL**            cols,               /**< LP columns */
   int                   ncols,              /**< number of LP columns */
   int                   ndisjcuts,          /**< number of disjunctive cuts found so far */
   SCIP_Bool             scale,              /**< should cut be scaled */
   SCIP_Bool             strengthen,         /**< should cut be strengthened if integer variables are present */
   SCIP_Real             cutlhs1,            /**< left hand side of the first simplex row */
   SCIP_Real             cutlhs2,            /**< left hand side of the second simplex row */
   SCIP_Real             bound1,             /**< bound of first simplex row */
   SCIP_Real             bound2,             /**< bound of second simplex row */
   SCIP_Real*            simplexcoefs1,      /**< simplex coefficients of first row */
   SCIP_Real*            simplexcoefs2,      /**< simplex coefficients of second row */
   SCIP_Real*            cutcoefs,           /**< pointer to store cut coefficients (length: nscipvars) */
   SCIP_ROW**            row,                /**< pointer to store disjunctive cut inequality */
   SCIP_Bool*            madeintegral        /**< pointer to store whether cut has been scaled to integral values */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_COL** rowcols;
   SCIP_COL* col;
   SCIP_Real* rowvals;
   SCIP_Real lhsrow;
   SCIP_Real rhsrow;
   SCIP_Real cutlhs;
   SCIP_Real sgn;
   SCIP_Real lb;
   SCIP_Real ub;
   int nonbasicnumber = 0;
   int rownnonz;
   int ind;
   int r;
   int c;

   assert( scip != NULL );
   assert( row != NULL );
   assert( rows != NULL );
   assert( cols != NULL );
   assert( simplexcoefs1 != NULL );
   assert( simplexcoefs2 != NULL );
   assert( cutcoefs != NULL );
   assert( sepa != NULL );
   assert( madeintegral != NULL );

   *madeintegral = FALSE;

   /* check signs */
   if ( SCIPisFeasPositive(scip, cutlhs1) == SCIPisFeasPositive(scip, cutlhs2) )
      sgn = 1.0;
   else
      sgn = -1.0;

   /* check bounds */
   if ( SCIPisInfinity(scip, REALABS(bound1)) || SCIPisInfinity(scip, REALABS(bound2)) )
      strengthen = FALSE;

   /* compute left hand side of row (a later update is possible, see below) */
   cutlhs = sgn * cutlhs1 * cutlhs2;

   /* add cut-coefficients of the non-basic non-slack variables */
   for (c = 0; c < ncols; ++c)
   {
      col = cols[c];
      assert( col != NULL );
      ind = SCIPcolGetLPPos(col);
      assert( ind >= 0 );

      if ( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_LOWER )
      {
         lb = SCIPcolGetLb(col);

         /* for integer variables we may obtain stronger coefficients */
         if ( strengthen && SCIPcolIsIntegral(col) )
         {
            SCIP_Real mval;
            SCIP_Real mvalfloor;
            SCIP_Real mvalceil;

            mval = (cutlhs2 * simplexcoefs1[nonbasicnumber] - cutlhs1 * simplexcoefs2[nonbasicnumber]) / (cutlhs2 * bound1 + cutlhs1 * bound2);
            mvalfloor = SCIPfloor(scip, mval);
            mvalceil = SCIPceil(scip, mval);

            cutcoefs[ind] = MIN(sgn * cutlhs2 * (simplexcoefs1[nonbasicnumber] - mvalfloor * bound1), sgn * cutlhs1 * (simplexcoefs2[nonbasicnumber] + mvalceil * bound2));
            assert( SCIPisFeasLE(scip, cutcoefs[ind], MAX(sgn * cutlhs2 * simplexcoefs1[nonbasicnumber], sgn * cutlhs1 * simplexcoefs2[nonbasicnumber])) );
         }
         else
            cutcoefs[ind] = MAX(sgn * cutlhs2 * simplexcoefs1[nonbasicnumber], sgn * cutlhs1 * simplexcoefs2[nonbasicnumber]);

         cutlhs += cutcoefs[ind] * lb;
         ++nonbasicnumber;
      }
      else if ( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_UPPER )
      {
         ub = SCIPcolGetUb(col);

         /* for integer variables we may obtain stronger coefficients */
         if ( strengthen && SCIPcolIsIntegral(col) )
         {
            SCIP_Real mval;
            SCIP_Real mvalfloor;
            SCIP_Real mvalceil;

            mval = (cutlhs2 * simplexcoefs1[nonbasicnumber] - cutlhs1 * simplexcoefs2[nonbasicnumber]) / (cutlhs2 * bound1 + cutlhs1 * bound2);
            mvalfloor = SCIPfloor(scip, -mval);
            mvalceil = SCIPceil(scip, -mval);

            cutcoefs[ind] = MAX(sgn * cutlhs2 * (simplexcoefs1[nonbasicnumber] + mvalfloor * bound1), sgn * cutlhs1 * (simplexcoefs2[nonbasicnumber] - mvalceil * bound2));
            assert( SCIPisFeasLE(scip, -cutcoefs[ind], -MIN(sgn * cutlhs2 * simplexcoefs1[nonbasicnumber], sgn * cutlhs1 * simplexcoefs2[nonbasicnumber])) );
         }
         else
            cutcoefs[ind] = MIN(sgn * cutlhs2 * simplexcoefs1[nonbasicnumber], sgn * cutlhs1 * simplexcoefs2[nonbasicnumber]);

         cutlhs += cutcoefs[ind] * ub;
         ++nonbasicnumber;
      }
      else
      {
         assert( SCIPcolGetBasisStatus(col) != SCIP_BASESTAT_ZERO );
         cutcoefs[ind] = 0.0;
      }
   }

   /* add cut-coefficients of the non-basic slack variables */
   for (r = 0; r < nrows; ++r)
   {
      rhsrow = SCIProwGetRhs(rows[r]) - SCIProwGetConstant(rows[r]);
      lhsrow = SCIProwGetLhs(rows[r]) - SCIProwGetConstant(rows[r]);

      assert( SCIProwGetBasisStatus(rows[r]) != SCIP_BASESTAT_ZERO );
      assert( SCIPisFeasZero(scip, lhsrow - rhsrow) || SCIPisNegative(scip, lhsrow - rhsrow) );
      assert( SCIProwIsInLP(rows[r]) );

      if ( SCIProwGetBasisStatus(rows[r]) != SCIP_BASESTAT_BASIC )
      {
         SCIP_Real cutcoef;

         if ( SCIProwGetBasisStatus(rows[r]) == SCIP_BASESTAT_UPPER )
         {
            assert( SCIPisFeasZero(scip, SCIPgetRowLPActivity(scip, rows[r]) - SCIProwGetRhs(rows[r])) );

            cutcoef = MAX(sgn * cutlhs2 * simplexcoefs1[nonbasicnumber], sgn * cutlhs1 * simplexcoefs2[nonbasicnumber]);
            cutlhs -= cutcoef * rhsrow;
            ++nonbasicnumber;
         }
         else /* SCIProwGetBasisStatus(rows[r]) == SCIP_BASESTAT_LOWER */
         {
            assert( SCIProwGetBasisStatus(rows[r]) == SCIP_BASESTAT_LOWER );
            assert( SCIPisFeasZero(scip, SCIPgetRowLPActivity(scip, rows[r]) - SCIProwGetLhs(rows[r])) );

            cutcoef = MIN(sgn * cutlhs2 * simplexcoefs1[nonbasicnumber], sgn * cutlhs1 * simplexcoefs2[nonbasicnumber]);
            cutlhs -= cutcoef * lhsrow;
            ++nonbasicnumber;
         }

         rownnonz = SCIProwGetNNonz(rows[r]);
         rowvals = SCIProwGetVals(rows[r]);
         rowcols = SCIProwGetCols(rows[r]);

         for (c = 0; c < rownnonz; ++c)
         {
            ind = SCIPcolGetLPPos(rowcols[c]);

            /* if column is not in LP, then return without generating cut */
            if ( ind < 0 )
            {
               *row = NULL;
               return SCIP_OKAY;
            }

            cutcoefs[ind] -= cutcoef * rowvals[c];
         }
      }
   }

   /* create cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_%d_%d", SCIPsepaGetName(sepa), SCIPgetNLPs(scip), ndisjcuts);
   if ( SCIPgetDepth(scip) == 0 )
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, row, sepa, cutname, cutlhs, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
   else
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, row, sepa, cutname, cutlhs, SCIPinfinity(scip), TRUE, FALSE, TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, *row) );
   for (c = 0; c < ncols; ++c)
   {
      ind = SCIPcolGetLPPos(cols[c]);
      assert( ind >= 0 );
      if ( ! SCIPisFeasZero(scip, cutcoefs[ind]) )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, *row, SCIPcolGetVar(cols[c]), cutcoefs[ind] ) );
      }
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, *row) );

   /* try to scale the cut to integral values
    * @todo find better but still stable disjunctive cut settings
    */
   if ( scale )
   {
      SCIP_Longint maxdnom;
      SCIP_Real maxscale;

      assert( SCIPgetDepth(scip) >= 0 );
      if( SCIPgetDepth(scip) == 0 )
      {
	 maxdnom = 100;
	 maxscale = 100.0;
      }
      else
      {
	 maxdnom = 10;
	 maxscale = 10.0;
      }

      SCIP_CALL( SCIPmakeRowIntegral(scip, *row, -SCIPepsilon(scip), SCIPsumepsilon(scip), maxdnom, maxscale, TRUE, madeintegral) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyDisjunctive)
{
   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaDisjunctive(scip) );

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeDisjunctive)/*lint --e{715}*/
{
   SCIP_SEPADATA* sepadata;

   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** solving process initialization method of separator */
static
SCIP_DECL_SEPAEXITSOL(sepaInitsolDisjunctive)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   sepadata->conshdlr = SCIPfindConshdlr(scip, "SOS1");

   return SCIP_OKAY;
}


/** LP solution separation method for disjunctive cuts */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpDisjunctive)
{
   SCIP_SEPADATA* sepadata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_DIGRAPH* conflictgraph;
   SCIP_ROW** rows;
   SCIP_COL** cols;
   SCIP_Real* cutcoefs = NULL;
   SCIP_Real* simplexcoefs1 = NULL;
   SCIP_Real* simplexcoefs2 = NULL;
   SCIP_Real* coef = NULL;
   SCIP_Real* binvrow = NULL;
   SCIP_Real* rowsmaxval = NULL;
   SCIP_Real* violationarray = NULL;
   int* fixings1 = NULL;
   int* fixings2 = NULL;
   int* basisind = NULL;
   int* basisrow = NULL;
   int* varrank = NULL;
   int* edgearray = NULL;
   int nedges;
   int ndisjcuts;
   int nrelevantedges;
   int nsos1vars;
   int nconss;
   int maxcuts;
   int ncalls;
   int depth;
   int ncols;
   int nrows;
   int ind;
   int j;
   int i;

   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   if( !allowlocal )
      return SCIP_OKAY;

   /* only generate disjunctive cuts if we are not close to terminating */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only generate disjunctive cuts if an optimal LP solution is at hand */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only generate disjunctive cuts if the LP solution is basic */
   if ( ! SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* get LP data */
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* return if LP has no columns or no rows */
   if ( ncols == 0 || nrows == 0 )
      return SCIP_OKAY;

   assert( cols != NULL );
   assert( rows != NULL );

   /* get sepa data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   /* get constraint handler */
   conshdlr = sepadata->conshdlr;
   if ( conshdlr == NULL )
      return SCIP_OKAY;

   /* get number of constraints */
   nconss = SCIPconshdlrGetNConss(conshdlr);
   if ( nconss == 0 )
      return SCIP_OKAY;

   /* check for maxdepth < depth, maxinvcutsroot = 0 and maxinvcuts = 0 */
   depth = SCIPgetDepth(scip);
   if ( ( sepadata->maxdepth >= 0 && sepadata->maxdepth < depth )
      || ( depth == 0 && sepadata->maxinvcutsroot == 0 )
      || ( depth > 0 && sepadata->maxinvcuts == 0 ) )
      return SCIP_OKAY;

   /* only call the cut separator a given number of times at each node */
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   if ( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* get conflict graph and number of conflict graph edges (note that the digraph arcs were added in both directions) */
   conflictgraph = SCIPgetConflictgraphSOS1(conshdlr);
   nedges = (int)SCIPceil(scip, (SCIP_Real)SCIPdigraphGetNArcs(conflictgraph)/2);

   /* if too many conflict graph edges, the separator can be slow: delay it until no other cuts have been found */
   if ( sepadata->maxconfsdelay >= 0 && nedges >= sepadata->maxconfsdelay )
   {
      int ncutsfound;

      ncutsfound = SCIPgetNCutsFound(scip);
      if ( ncutsfound > sepadata->lastncutsfound || ! SCIPsepaWasLPDelayed(sepa) )
      {
         sepadata->lastncutsfound = ncutsfound;
         *result = SCIP_DELAYED;
         return SCIP_OKAY;
      }
   }

   /* check basis status */
   for (j = 0; j < ncols; ++j)
   {
      if ( SCIPcolGetBasisStatus(cols[j]) == SCIP_BASESTAT_ZERO )
         return SCIP_OKAY;
   }

   /* check accuracy of rows */
   for (j = 0; j < nrows; ++j)
   {
      SCIP_ROW* row;

      row = rows[j];

      if ( ( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_UPPER  && ! SCIPisEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetRhs(row)) )
           || ( SCIProwGetBasisStatus(row) == SCIP_BASESTAT_LOWER && ! SCIPisEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row)) ) )
         return SCIP_OKAY;
   }

   /* get number of SOS1 variables */
   nsos1vars = SCIPgetNSOS1Vars(conshdlr);

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearray, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixings1, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixings2, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violationarray, nedges) );

   /* get all violated conflicts {i, j} in the conflict graph and sort them based on the degree of a violation value */
   nrelevantedges = 0;
   for (j = 0; j < nsos1vars; ++j)
   {
      SCIP_VAR* var;

      var = SCIPnodeGetVarSOS1(conflictgraph, j);

      if ( SCIPvarIsActive(var) && ! SCIPisFeasZero(scip, SCIPcolGetPrimsol(SCIPvarGetCol(var))) && SCIPcolGetBasisStatus(SCIPvarGetCol(var)) == SCIP_BASESTAT_BASIC )
      {
         int* succ;
         int nsucc;

         /* get successors and number of successors */
         nsucc = SCIPdigraphGetNSuccessors(conflictgraph, j);
         succ = SCIPdigraphGetSuccessors(conflictgraph, j);

         for (i = 0; i < nsucc; ++i)
         {
            SCIP_VAR* varsucc;
            int succind;

            succind = succ[i];
            varsucc = SCIPnodeGetVarSOS1(conflictgraph, succind);
            if ( SCIPvarIsActive(varsucc) && succind < j && ! SCIPisFeasZero(scip, SCIPgetSolVal(scip, NULL, varsucc) ) &&
                 SCIPcolGetBasisStatus(SCIPvarGetCol(varsucc)) == SCIP_BASESTAT_BASIC )
            {
               fixings1[nrelevantedges] = j;
               fixings2[nrelevantedges] = succind;
               edgearray[nrelevantedges] = nrelevantedges;
               violationarray[nrelevantedges++] = SCIPgetSolVal(scip, NULL, var) * SCIPgetSolVal(scip, NULL, varsucc);
            }
         }
      }
   }

   /* sort violation score values */
   if ( nrelevantedges > 0)
      SCIPsortDownRealInt(violationarray, edgearray, nrelevantedges);
   else
   {
      SCIPfreeBufferArrayNull(scip, &violationarray);
      SCIPfreeBufferArrayNull(scip, &fixings2);
      SCIPfreeBufferArrayNull(scip, &fixings1);
      SCIPfreeBufferArrayNull(scip, &edgearray);

      return SCIP_OKAY;
   }
   SCIPfreeBufferArrayNull(scip, &violationarray);

   /* compute maximal number of cuts */
   if ( SCIPgetDepth(scip) == 0 )
      maxcuts = MIN(sepadata->maxinvcutsroot, nrelevantedges);
   else
      maxcuts = MIN(sepadata->maxinvcuts, nrelevantedges);
   assert( maxcuts > 0 );

   /* allocate buffer arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &varrank, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowsmaxval, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisrow, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coef, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &simplexcoefs1, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &simplexcoefs2, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );

   /* get basis indices */
   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );

   /* create vector "basisrow" with basisrow[column of non-slack basis variable] = corresponding row of B^-1;
    * compute maximum absolute value of nonbasic row coefficients */
   for (j = 0; j < nrows; ++j)
   {
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_ROW* row;
      SCIP_Real val;
      SCIP_Real max = 0.0;
      int nnonz;

      /* fill basisrow vector */
      ind = basisind[j];
      if ( ind >= 0 )
         basisrow[ind] = j;

      /* compute maximum absolute value of nonbasic row coefficients */
      row = rows[j];
      assert( row != NULL );
      rowvals = SCIProwGetVals(row);
      nnonz = SCIProwGetNNonz(row);
      rowcols = SCIProwGetCols(row);

      for (i = 0; i < nnonz; ++i)
      {
         if ( SCIPcolGetBasisStatus(rowcols[i]) == SCIP_BASESTAT_LOWER  || SCIPcolGetBasisStatus(rowcols[i]) == SCIP_BASESTAT_UPPER )
         {
            val = REALABS(rowvals[i]);
            if ( SCIPisFeasGT(scip, val, max) )
               max = REALABS(val);
         }
      }

      /* handle slack variable coefficient and save maximum value */
      rowsmaxval[j] = MAX(max, 1.0);
   }

   /* initialize variable ranks with -1 */
   for (j = 0; j < ncols; ++j)
      varrank[j] = -1;

   /* free buffer array */
   SCIPfreeBufferArrayNull(scip, &basisind);

   /* for the most promising disjunctions: try to generate disjunctive cuts */
   ndisjcuts = 0;
   for (i = 0; i < maxcuts; ++i)
   {
      SCIP_Bool madeintegral;
      SCIP_Real cutlhs1;
      SCIP_Real cutlhs2;
      SCIP_Real bound1;
      SCIP_Real bound2;
      SCIP_ROW* row = NULL;
      SCIP_VAR* var;
      SCIP_COL* col;

      int nonbasicnumber;
      int cutrank = 0;
      int edgenumber;
      int rownnonz;

      edgenumber = edgearray[i];

      /* determine first simplex row */
      var = SCIPnodeGetVarSOS1(conflictgraph, fixings1[edgenumber]);
      col = SCIPvarGetCol(var);
      ind = SCIPcolGetLPPos(col);
      assert( ind >= 0 );
      assert( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC );

      /* get the 'ind'th row of B^-1 and B^-1 \cdot A */
      SCIP_CALL( SCIPgetLPBInvRow(scip, basisrow[ind], binvrow, NULL, NULL) );
      SCIP_CALL( SCIPgetLPBInvARow(scip, basisrow[ind], binvrow, coef, NULL, NULL) );

      /* get the simplex-coefficients of the non-basic variables */
      SCIP_CALL( getSimplexCoefficients(scip, rows, nrows, cols, ncols, coef, binvrow, simplexcoefs1, &nonbasicnumber) );

      /* get rank of variable if not known already */
      if ( varrank[ind] < 0 )
         varrank[ind] = getVarRank(scip, binvrow, rowsmaxval, sepadata->maxweightrange, rows, nrows);
      cutrank = MAX(cutrank, varrank[ind]);

      /* get right hand side and bound of simplex talbeau row */
      cutlhs1 = SCIPcolGetPrimsol(col);
      if ( SCIPisFeasPositive(scip, cutlhs1) )
         bound1 = SCIPcolGetUb(col);
      else
         bound1 = SCIPcolGetLb(col);


      /* determine second simplex row */
      var = SCIPnodeGetVarSOS1(conflictgraph, fixings2[edgenumber]);
      col = SCIPvarGetCol(var);
      ind = SCIPcolGetLPPos(col);
      assert( ind >= 0 );
      assert( SCIPcolGetBasisStatus(col) == SCIP_BASESTAT_BASIC );

      /* get the 'ind'th row of B^-1 and B^-1 \cdot A */
      SCIP_CALL( SCIPgetLPBInvRow(scip, basisrow[ind], binvrow, NULL, NULL) );
      SCIP_CALL( SCIPgetLPBInvARow(scip, basisrow[ind], binvrow, coef, NULL, NULL) );

      /* get the simplex-coefficients of the non-basic variables */
      SCIP_CALL( getSimplexCoefficients(scip, rows, nrows, cols, ncols, coef, binvrow, simplexcoefs2, &nonbasicnumber) );

      /* get rank of variable if not known already */
      if ( varrank[ind] < 0 )
         varrank[ind] = getVarRank(scip, binvrow, rowsmaxval, sepadata->maxweightrange, rows, nrows);
      cutrank = MAX(cutrank, varrank[ind]);

      /* get right hand side and bound of simplex talbeau row */
      cutlhs2 = SCIPcolGetPrimsol(col);
      if ( SCIPisFeasPositive(scip, cutlhs2) )
         bound2 = SCIPcolGetUb(col);
      else
         bound2 = SCIPcolGetLb(col);

      /* add coefficients to cut */
      SCIP_CALL( generateDisjCutSOS1(scip, sepa, rows, nrows, cols, ncols, ndisjcuts, MAKECONTINTEGRAL, sepadata->strengthen, cutlhs1, cutlhs2, bound1, bound2, simplexcoefs1, simplexcoefs2, cutcoefs, &row, &madeintegral) );
      if ( row == NULL )
         continue;

      /* raise cutrank for present cut */
      ++cutrank;

      /* check if there are numerical evidences */
      if ( ( madeintegral && ( sepadata->maxrankintegral == -1 || cutrank <= sepadata->maxrankintegral ) )
         || ( ! madeintegral && ( sepadata->maxrank == -1 || cutrank <= sepadata->maxrank ) ) )
      {
         /* possibly add cut to LP if it is useful; in case the lhs of the cut is minus infinity (due to scaling) the cut is useless */
         rownnonz = SCIProwGetNNonz(row);
         if ( rownnonz > 0 && ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) && ! SCIProwIsInLP(row) && SCIPisCutEfficacious(scip, NULL, row) )
         {
            SCIP_Bool infeasible;

            /* set cut rank */
            SCIProwChgRank(row, cutrank);

            /* add cut */
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            if ( infeasible )
            {
               *result = SCIP_CUTOFF;
               break;
            }
            ++ndisjcuts;
         }
      }

      /* release row */
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   /* save total number of cuts found so far */
   sepadata->lastncutsfound = SCIPgetNCutsFound(scip);

   /* evaluate the result of the separation */
   if ( *result != SCIP_CUTOFF )
   {
      if ( ndisjcuts > 0 )
         *result = SCIP_SEPARATED;
      else
         *result = SCIP_DIDNOTFIND;
   }

   SCIPdebugMsg(scip, "Number of found disjunctive cuts: %d.\n", ndisjcuts);

   /* free buffer arrays */
   SCIPfreeBufferArrayNull(scip, &cutcoefs);
   SCIPfreeBufferArrayNull(scip, &simplexcoefs2);
   SCIPfreeBufferArrayNull(scip, &simplexcoefs1);
   SCIPfreeBufferArrayNull(scip, &coef);
   SCIPfreeBufferArrayNull(scip, &binvrow);
   SCIPfreeBufferArrayNull(scip, &basisrow);
   SCIPfreeBufferArrayNull(scip, &fixings2);
   SCIPfreeBufferArrayNull(scip, &fixings1);
   SCIPfreeBufferArrayNull(scip, &edgearray);
   SCIPfreeBufferArrayNull(scip, &rowsmaxval);
   SCIPfreeBufferArrayNull(scip, &varrank);

   return SCIP_OKAY;
}


/** creates the disjunctive cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaDisjunctive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata = NULL;
   SCIP_SEPA* sepa = NULL;

   /* create separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepadata->conshdlr = NULL;
   sepadata->lastncutsfound = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpDisjunctive, NULL, sepadata) );

   assert( sepa != NULL );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyDisjunctive) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeDisjunctive) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolDisjunctive) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/strengthen",
         "strengthen cut if integer variables are present.",
         &sepadata->strengthen, TRUE, DEFAULT_STRENGTHEN, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxdepth",
         "node depth of separating bipartite disjunctive cuts (-1: no limit)",
         &sepadata->maxdepth, TRUE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxrounds",
         "maximal number of separation rounds per iteration in a branching node (-1: no limit)",
         &sepadata->maxrounds, TRUE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxroundsroot",
         "maximal number of separation rounds in the root node (-1: no limit)",
         &sepadata->maxroundsroot, TRUE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxinvcuts",
         "maximal number of cuts investigated per iteration in a branching node",
         &sepadata->maxinvcuts, TRUE, DEFAULT_MAXINVCUTS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxinvcutsroot",
         "maximal number of cuts investigated per iteration in the root node",
         &sepadata->maxinvcutsroot, TRUE, DEFAULT_MAXINVCUTSROOT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxconfsdelay",
         "delay separation if number of conflict graph edges is larger than predefined value (-1: no limit)",
         &sepadata->maxconfsdelay, TRUE, DEFAULT_MAXCONFSDELAY, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxrank",
         "maximal rank of a disj. cut that could not be scaled to integral coefficients (-1: unlimited)",
         &sepadata->maxrank, FALSE, DEFAULT_MAXRANK, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxrankintegral",
         "maximal rank of a disj. cut that could be scaled to integral coefficients (-1: unlimited)",
         &sepadata->maxrankintegral, FALSE, DEFAULT_MAXRANKINTEGRAL, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/maxweightrange",
         "maximal valid range max(|weights|)/min(|weights|) of row weights",
         &sepadata->maxweightrange, TRUE, DEFAULT_MAXWEIGHTRANGE, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
