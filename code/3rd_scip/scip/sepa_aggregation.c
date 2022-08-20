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

/**@file   sepa_aggregation.c
 * @brief  flow cover and complemented mixed integer rounding cuts separator (Marchand's version)
 * @author Robert Lion Gottwald
 * @author Kati Wolter
 * @author Tobias Achterberg
 *
 * For an overview see:
 *
 * Marchand, H., & Wolsey, L. A. (2001).@n
 * Aggregation and mixed integer rounding to solve MIPs.@n
 * Operations research, 49(3), 363-371.
 *
 * Some remarks:
 * - In general, continuous variables are less prefered than integer variables, since their cut
 *   coefficient is worse.
 * - We seek for aggregations that project out continuous variables that are far away from their bound,
 *   since if it is at its bound then it doesn't contribute to the violation
 * - These aggregations are also useful for the flowcover separation, so after building an aggregation
 *   we try to generate a MIR cut and a flowcover cut.
 * - We only keep the best cut.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_aggregation.h"
#include "scip/pub_misc.h"
#include "scip/cuts.h"


#define SEPA_NAME              "aggregation"
#define SEPA_DESC              "aggregation heuristic for complemented mixed integer rounding cuts and flowcover cuts"
#define SEPA_PRIORITY             -3000
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS            -1 /**< maximal number of cmir separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of cmir separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXTRIES            200 /**< maximal number of rows to start aggregation with per separation round
                                         *   (-1: unlimited) */
#define DEFAULT_MAXTRIESROOT         -1 /**< maximal number of rows to start aggregation with per round in the root node
                                         *   (-1: unlimited) */
#define DEFAULT_MAXFAILS             20 /**< maximal number of consecutive unsuccessful aggregation tries (-1: unlimited) */
#define DEFAULT_MAXFAILSROOT        100 /**< maximal number of consecutive unsuccessful aggregation tries in the root node
                                         *   (-1: unlimited) */
#define DEFAULT_MAXAGGRS              3 /**< maximal number of aggregations for each row per separation round */
#define DEFAULT_MAXAGGRSROOT          6 /**< maximal number of aggregations for each row per round in the root node */
#define DEFAULT_MAXSEPACUTS         100 /**< maximal number of cmir cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     500 /**< maximal number of cmir cuts separated per separation round in root node */
#define DEFAULT_MAXSLACK            0.0 /**< maximal slack of rows to be used in aggregation */
#define DEFAULT_MAXSLACKROOT        0.1 /**< maximal slack of rows to be used in aggregation in the root node */
#define DEFAULT_DENSITYSCORE       1e-4 /**< weight of row density in the aggregation scoring of the rows */
#define DEFAULT_SLACKSCORE         1e-3 /**< weight of slack in the aggregation scoring of the rows */
#define DEFAULT_MAXAGGDENSITY      0.20 /**< maximal density of aggregated row */
#define DEFAULT_MAXROWDENSITY      0.05 /**< maximal density of row to be used in aggregation */
#define DEFAULT_DENSITYOFFSET       100 /**< additional number of variables allowed in row on top of density */
#define DEFAULT_MAXROWFAC          1e+4 /**< maximal row aggregation factor */
#define DEFAULT_MAXTESTDELTA         -1 /**< maximal number of different deltas to try (-1: unlimited) */
#define DEFAULT_AGGRTOL            1e-2 /**< aggregation heuristic: we try to delete continuous variables from the current
                                         *   aggregation, whose distance to its tightest bound is >= L - DEFAULT_AGGRTOL,
                                         *   where L is the largest of the distances between a continuous variable's value
                                         *   and its tightest bound in the current aggregation */
#define DEFAULT_TRYNEGSCALING      TRUE /**< should negative values also be tested in scaling? */
#define DEFAULT_FIXINTEGRALRHS     TRUE /**< should an additional variable be complemented if f0 = 0? */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */

#define BOUNDSWITCH                 0.5
#define POSTPROCESS                TRUE
#define USEVBDS                    TRUE
#define MINFRAC                    0.05
#define MAXFRAC                    0.999
#define MAKECONTINTEGRAL          FALSE
#define IMPLINTSARECONT


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_Real             maxslack;           /**< maximal slack of rows to be used in aggregation */
   SCIP_Real             maxslackroot;       /**< maximal slack of rows to be used in aggregation in the root node */
   SCIP_Real             densityscore;       /**< weight of row density in the aggregation scoring of the rows */
   SCIP_Real             slackscore;         /**< weight of slack in the aggregation scoring of the rows */
   SCIP_Real             maxaggdensity;      /**< maximal density of aggregated row */
   SCIP_Real             maxrowdensity;      /**< maximal density of row to be used in aggregation */
   SCIP_Real             maxrowfac;          /**< maximal row aggregation factor */
   SCIP_Real             aggrtol;            /**< tolerance for bound distance used in aggregation heuristic */
   int                   maxrounds;          /**< maximal number of cmir separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of cmir separation rounds in the root node (-1: unlimited) */
   int                   maxtries;           /**< maximal number of rows to start aggregation with per separation round
                                              *   (-1: unlimited) */
   int                   maxtriesroot;       /**< maximal number of rows to start aggregation with per round in the root node
                                              *   (-1: unlimited) */
   int                   maxfails;           /**< maximal number of consecutive unsuccessful aggregation tries
                                              *   (-1: unlimited) */
   int                   maxfailsroot;       /**< maximal number of consecutive unsuccessful aggregation tries in the root
                                              *   node (-1: unlimited) */
   int                   maxaggrs;           /**< maximal number of aggregations for each row per separation round */
   int                   maxaggrsroot;       /**< maximal number of aggregations for each row per round in the root node */
   int                   maxsepacuts;        /**< maximal number of cmir cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cmir cuts separated per separation round in root node */
   int                   densityoffset;      /**< additional number of variables allowed in row on top of density */
   int                   maxtestdelta;       /**< maximal number of different deltas to try (-1: unlimited) */
   SCIP_Bool             trynegscaling;      /**< should negative values also be tested in scaling? */
   SCIP_Bool             fixintegralrhs;     /**< should an additional variable be complemented if f0 = 0? */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_SEPA*            cmir;               /**< separator for adding cmir cuts */
   SCIP_SEPA*            flowcover;          /**< separator for adding flowcover cuts */
};

typedef
struct AggregationData {
   SCIP_Real*            bounddist;          /**< bound distance of continuous variables */
   int*                  bounddistinds;      /**< problem indices of the continUous variables corresponding to the bounddistance value */
   int                   nbounddistvars;     /**< number of continuous variables that are not at their bounds */
   SCIP_ROW**            aggrrows;           /**< array of rows suitable for substitution of continuous variable */
   SCIP_Real*            aggrrowscoef;       /**< coefficient of continuous variable in row that is suitable for substitution of that variable */
   int                   aggrrowssize;       /**< size of aggrrows array */
   int                   naggrrows;          /**< occupied positions in aggrrows array */
   int*                  aggrrowsstart;      /**< array with start positions of suitable rows for substitution for each
                                              *   continuous variable with non-zero bound distance */
   int*                  ngoodaggrrows;      /**< array with number of rows suitable for substitution that only contain
                                              *   one continuous variable that is not at it's bound */
   int*                  nbadvarsinrow;      /**< number of continuous variables that are not at their bounds for each row */
   SCIP_AGGRROW*         aggrrow;            /**< store aggregation row here so that it can be reused */
} AGGREGATIONDATA;

/*
 * Local methods
 */

/** adds given cut to LP if violated */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_Bool             makeintegral,       /**< should cut be scaled to integral coefficients if possible? */
   SCIP_Real*            cutcoefs,           /**< coefficients of active variables in cut */
   int*                  cutinds,            /**< problem indices of variables in cut */
   int                   cutnnz,             /**< number of non-zeros in cut */
   SCIP_Real             cutrhs,             /**< right hand side of cut */
   SCIP_Real             cutefficacy,        /**< efficacy of cut */
   SCIP_Bool             cutislocal,         /**< is the cut only locally valid? */
   SCIP_Bool             cutremovable,       /**< should the cut be removed from the LP due to aging or cleanup? */
   int                   cutrank,            /**< rank of the cut */
   const char*           cutclassname,       /**< name of cut class to use for row names */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts,              /**< pointer to count the number of added cuts */
   SCIP_ROW**            thecut              /**< pointer to return cut if it was added */
   )
{
   assert(scip != NULL);
   assert(cutcoefs != NULL);
   assert(cutoff != NULL);
   assert(ncuts != NULL);

   *cutoff = FALSE;

   if( cutnnz > 0 && SCIPisEfficacious(scip, cutefficacy) )
   {
      SCIP_VAR** vars;
      int i;
      SCIP_ROW* cut;
      char cutname[SCIP_MAXSTRLEN];
      SCIP_Bool success;

      /* get active problem variables */
      vars = SCIPgetVars(scip);

      /* create the cut */
      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s%d_%d", cutclassname, SCIPgetNLPs(scip), *ncuts);

tryagain:
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs,
                                        cutislocal, FALSE, cutremovable) );

      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

      for( i = 0; i < cutnnz; ++i )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[i]], cutcoefs[i]) );
      }

      /* set cut rank */
      SCIProwChgRank(cut, cutrank);

      SCIPdebugMsg(scip, " -> found potential %s cut <%s>: rhs=%f, eff=%f\n",
                   cutclassname, cutname, cutrhs, cutefficacy);
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );


      /* if requested, try to scale the cut to integral values  but only if the scaling is small; otherwise keep the fractional cut */
      if( makeintegral && SCIPgetRowNumIntCols(scip, cut) == SCIProwGetNNonz(cut) )
      {
         SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
               1000LL, 1000.0, MAKECONTINTEGRAL, &success) );

         if( SCIPisInfinity(scip, SCIProwGetRhs(cut)) )
         {
            /* release the row */
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );

            /* the scaling destroyed the cut, so try to add it again but this time do not scale it */
            makeintegral = FALSE;
            goto tryagain;
         }
      }
      else
      {
         success = FALSE;
      }

      if( success && !SCIPisCutEfficacious(scip, sol, cut) )
      {
         SCIPdebugMsg(scip, " -> %s cut <%s> no longer efficacious: rhs=%f, eff=%f\n",
                      cutclassname, cutname, cutrhs, cutefficacy);
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );

         SCIP_CALL( SCIPreleaseRow(scip, &cut) );

         /* the cut is not efficacious anymore due to the scaling so do not add it */
         return SCIP_OKAY;
      }

      SCIPdebugMsg(scip, " -> found %s cut <%s>: rhs=%f, eff=%f, rank=%d, min=%f, max=%f (range=%g)\n",
                     cutclassname, cutname, cutrhs, cutefficacy, SCIProwGetRank(cut),
                     SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                     SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );

      SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

      if( SCIPisCutNew(scip, cut) )
      {
         (*ncuts)++;

         if( !cutislocal )
         {
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
         }
         else
         {
            SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
         }

         *thecut = cut;
      }
      else
      {
         /* release the row */
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      }
   }

   return SCIP_OKAY;
}

/** setup data for aggregating rows */
static
SCIP_RETCODE setupAggregationData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to separate, NULL for LP solution */
   SCIP_Bool             allowlocal,         /**< should local cuts be allowed */
   AGGREGATIONDATA*      aggrdata            /**< pointer to aggregation data to setup */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nintvars;
   int ncontvars;
   int firstcontvar;
   int nimplvars;
   SCIP_ROW** rows;
   int nrows;
   int i;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, &ncontvars) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   SCIP_CALL( SCIPallocBufferArray(scip, &aggrdata->bounddist, ncontvars + nimplvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrdata->bounddistinds, ncontvars + nimplvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrdata->ngoodaggrrows, ncontvars + nimplvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrdata->aggrrowsstart, ncontvars + nimplvars + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrdata->nbadvarsinrow, nrows) );
   SCIP_CALL( SCIPaggrRowCreate(scip, &aggrdata->aggrrow) );
   BMSclearMemoryArray(aggrdata->nbadvarsinrow, nrows);

   aggrdata->nbounddistvars = 0;
   aggrdata->aggrrows = NULL;
   aggrdata->aggrrowscoef = NULL;
   aggrdata->aggrrowssize = 0;
   aggrdata->naggrrows = 0;

   firstcontvar = nvars - ncontvars;

   for( i = nbinvars + nintvars; i < nvars; ++i )
   {
      SCIP_Real bounddist;

      /* compute the bound distance of the variable */
      {
         SCIP_Real primsol;
         SCIP_Real distlb;
         SCIP_Real distub;
         SCIP_Real bestlb;
         SCIP_Real bestub;
         SCIP_Real bestvlb;
         SCIP_Real bestvub;
         int bestvlbidx;
         int bestvubidx;

         if( allowlocal )
         {
            bestlb = SCIPvarGetLbLocal(vars[i]);
            bestub = SCIPvarGetUbLocal(vars[i]);
         }
         else
         {
            bestlb = SCIPvarGetLbGlobal(vars[i]);
            bestub = SCIPvarGetUbGlobal(vars[i]);
         }

         SCIP_CALL( SCIPgetVarClosestVlb(scip, vars[i], sol, &bestvlb, &bestvlbidx) );
         SCIP_CALL( SCIPgetVarClosestVub(scip, vars[i], sol, &bestvub, &bestvubidx) );
         if( bestvlbidx >= 0 )
            bestlb = MAX(bestlb, bestvlb);
         if( bestvubidx >= 0 )
            bestub = MIN(bestub, bestvub);

         primsol = SCIPgetSolVal(scip, sol, vars[i]);
         distlb = primsol - bestlb;
         distub = bestub - primsol;

         bounddist = MIN(distlb, distub);
         bounddist = MAX(bounddist, 0.0);

         /* prefer continuous variables over implicit integers to be aggregated out */
         if( i < firstcontvar )
            bounddist *= 0.1;
      }

      /* when variable is not at its bound, we want to project it out, so add it to the aggregation data */
      if( !SCIPisZero(scip, bounddist) )
      {
         int k = aggrdata->nbounddistvars++;
         aggrdata->bounddist[k] = bounddist;
         aggrdata->bounddistinds[k] = i;
         aggrdata->aggrrowsstart[k] = aggrdata->naggrrows;

         /* the current variable is a bad variable (continuous, not at its bound): increase the number of bad variable
          * count on each row this variables appears in; also each of these rows can be used to project the variable out
          * so store them.
          */
         if( SCIPvarIsInLP(vars[i]) )
         {
            SCIP_COL* col = SCIPvarGetCol(vars[i]);
            SCIP_ROW** colrows = SCIPcolGetRows(col);
            SCIP_Real* colrowvals = SCIPcolGetVals(col);
            int ncolnonzeros = SCIPcolGetNLPNonz(col);
            int aggrrowsminsize = aggrdata->naggrrows + ncolnonzeros;

            if( aggrrowsminsize > aggrdata->aggrrowssize )
            {
               SCIP_CALL( SCIPreallocBufferArray(scip, &aggrdata->aggrrows, aggrrowsminsize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &aggrdata->aggrrowscoef, aggrrowsminsize) );
            }

            for( k = 0; k < ncolnonzeros; ++k )
            {
               /* ignore modifiable rows */
               if( SCIProwIsModifiable(colrows[k]) )
                  continue;

               ++aggrdata->nbadvarsinrow[SCIProwGetLPPos(colrows[k])];
               aggrdata->aggrrows[aggrdata->naggrrows] = colrows[k];
               aggrdata->aggrrowscoef[aggrdata->naggrrows] = colrowvals[k];
               ++aggrdata->naggrrows;
            }
         }
      }
   }

   /* add sentinel entry at the end */
   aggrdata->aggrrowsstart[aggrdata->nbounddistvars] = aggrdata->naggrrows;

   /* for each continous variable that is not at its bounds check if there is a
    * row where it is the only such variable ("good" rows). In the array with the rows that are
    * suitable for substituting this variable move the good rows to the beginning
    * and store the number of good rows for each of the variables.
    * If a variable has at least one good row, then it is a "better" variable and we make
    * the value of the bounddistance for this variable negative, to mark it.
    * Note that better variables are continous variables that are not at their bounds
    * and can be projected out without introducing bad variables (by using a good row).
    */
   {
      int beg;
      int end;

      beg = aggrdata->aggrrowsstart[0];
      for( i = 0; i < aggrdata->nbounddistvars; ++i )
      {
         int k;
         int ngoodrows;

         end = aggrdata->aggrrowsstart[i + 1];
         ngoodrows = 0;
         for( k = beg; k < end; ++k )
         {
            int lppos = SCIProwGetLPPos(aggrdata->aggrrows[k]);
            if( aggrdata->nbadvarsinrow[lppos] == 1
               && SCIPisEQ(scip, SCIProwGetLhs(aggrdata->aggrrows[k]), SCIProwGetRhs(aggrdata->aggrrows[k])) )
            {
               int nextgoodrowpos = beg + ngoodrows;
               if( k > nextgoodrowpos )
               {
                  SCIPswapPointers((void**) (&aggrdata->aggrrows[k]), (void**) (&aggrdata->aggrrows[nextgoodrowpos]));
                  SCIPswapReals(&aggrdata->aggrrowscoef[k], &aggrdata->aggrrowscoef[nextgoodrowpos]);
               }
               ++ngoodrows;
            }
         }
         if( ngoodrows > 0 )
         {
            aggrdata->bounddist[i] = -aggrdata->bounddist[i];
         }
         aggrdata->ngoodaggrrows[i] = ngoodrows;
         beg = end;
      }
   }

   return SCIP_OKAY;
}

/** free resources held in aggregation data */
static
void destroyAggregationData(
   SCIP*                 scip,               /**< SCIP datastructure */
   AGGREGATIONDATA*      aggrdata            /**< pointer to ggregation data */
   )
{
   SCIPaggrRowFree(scip, &aggrdata->aggrrow);
   SCIPfreeBufferArrayNull(scip, &aggrdata->aggrrowscoef);
   SCIPfreeBufferArrayNull(scip, &aggrdata->aggrrows);
   SCIPfreeBufferArray(scip, &aggrdata->nbadvarsinrow);
   SCIPfreeBufferArray(scip, &aggrdata->aggrrowsstart);
   SCIPfreeBufferArray(scip, &aggrdata->ngoodaggrrows);
   SCIPfreeBufferArray(scip, &aggrdata->bounddistinds);
   SCIPfreeBufferArray(scip, &aggrdata->bounddist);
}

/** retrieves the candidate rows for canceling out the given variable, also returns the number of "good" rows which are the
 *  rows stored at the first ngoodrows positions. A row is good if its continuous variables are all at their bounds, except
 *  maybe the given continuous variable (in probvaridx)
 */
static
SCIP_Bool getRowAggregationCandidates(
   AGGREGATIONDATA*      aggrdata,           /**< pointer to ggregation data */
   int                   probvaridx,         /**< problem index of variables to retrieve candidates for */
   SCIP_ROW***           rows,               /**< pointer to store array to candidate rows */
   SCIP_Real**           rowvarcoefs,        /**< pointer to store array of coefficients of given variable in the corresponding rows */
   int*                  nrows,              /**< pointer to return number of rows in returned arrays */
   int*                  ngoodrows           /**< pointer to return number of "good" rows in the returned arrays */
   )
{
   int aggrdataidx;

   if( !SCIPsortedvecFindInt(aggrdata->bounddistinds, probvaridx, aggrdata->nbounddistvars, &aggrdataidx) )
      return FALSE;

   *rows = aggrdata->aggrrows + aggrdata->aggrrowsstart[aggrdataidx];
   *nrows = aggrdata->aggrrowsstart[aggrdataidx + 1] - aggrdata->aggrrowsstart[aggrdataidx];
   *rowvarcoefs = aggrdata->aggrrowscoef + aggrdata->aggrrowsstart[aggrdataidx];
   *ngoodrows = aggrdata->ngoodaggrrows[aggrdataidx];

   return TRUE;
}

/** find the bound distance value in the aggregation data struct for the given variable problem index */
static
SCIP_Real aggrdataGetBoundDist(
   AGGREGATIONDATA*      aggrdata,           /**< SCIP datastructure */
   int                   probvaridx          /**< problem index of variables to retrieve candidates for */
   )
{
   int aggrdataidx;

   if( !SCIPsortedvecFindInt(aggrdata->bounddistinds, probvaridx, aggrdata->nbounddistvars, &aggrdataidx) )
      return 0.0;

   return aggrdata->bounddist[aggrdataidx];
}

/** Aggregates the next row suitable for cancelling out an active continuous variable.
 *  Equality rows that contain no other active continuous variables are preffered and apart from that
 *  the scores for the rows are used to determine which row is aggregated next
 */
static
SCIP_RETCODE aggregateNextRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_Real*            rowlhsscores,       /**< aggregation scores for left hand sides of row */
   SCIP_Real*            rowrhsscores,       /**< aggregation scores for right hand sides of row */
   AGGREGATIONDATA*      aggrdata,           /**< aggregation data */
   SCIP_AGGRROW*         aggrrow,            /**< current aggregation row */
   int*                  naggrs,             /**< pointer to increase counter if real aggregation took place */
   SCIP_Bool*            success             /**< pointer to return whether another row was added to the aggregation row */
   )
{
   int i;
   int firstcontvar;
   int* badvarinds;
   SCIP_Real* badvarbddist;
   int nbadvars;
   SCIP_Real minbddist;
   SCIP_ROW* bestrow;
   SCIP_Real bestrowscore;
   SCIP_Real aggrfac;
   int bestrowside;

   int nnz = SCIPaggrRowGetNNz(aggrrow);
   int* inds = SCIPaggrRowGetInds(aggrrow);

   *success = FALSE;

   {
      int nbinvars;
      int nintvars;
      int ncontvars;

      nbinvars = SCIPgetNBinVars(scip);
      nintvars = SCIPgetNIntVars(scip);
      firstcontvar =  nbinvars + nintvars;
      ncontvars = SCIPgetNVars(scip) - firstcontvar;

      SCIP_CALL( SCIPallocBufferArray(scip, &badvarinds, MIN(ncontvars, nnz)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &badvarbddist, MIN(ncontvars, nnz)) );
   }

   nbadvars = 0;

   for( i = 0; i < nnz; ++i )
   {
      SCIP_Real bounddist;

      /* only consider continuous variables */
      if( inds[i] < firstcontvar )
         continue;

      bounddist = aggrdataGetBoundDist(aggrdata, inds[i]);

      if( bounddist == 0.0 )
         continue;

      badvarinds[nbadvars] = inds[i];
      badvarbddist[nbadvars] = bounddist;
      ++nbadvars;
   }

   if( nbadvars == 0 )
      goto TERMINATE;

   SCIPsortDownRealInt(badvarbddist, badvarinds, nbadvars);

   aggrfac = 0.0;
   bestrowscore = 0.0;
   bestrowside = 0;
   minbddist = 0.0;
   bestrow = NULL;

   /* because the "good" bad variables have a negative bound distance, they are at the end */
   for( i = nbadvars - 1; i >= 0; --i )
   {
      int probvaridx;
      SCIP_ROW** candrows;
      SCIP_Real* candrowcoefs;
      int nrows;
      int ngoodrows;
      int k;

      /* if the bound distance is not negative, there are no more good variables so stop */
      if( badvarbddist[i] > 0.0 )
         break;

      /* if no best row was found yet, this variable has the currently best bound distance */
      if( aggrfac == 0.0 )
         minbddist = -badvarbddist[i] * (1.0 - sepadata->aggrtol);

      /* if the bound distance of the current variable is smaller than the minimum bound distance stop looping */
      if( -badvarbddist[i] < minbddist )
         break;

      probvaridx = badvarinds[i];

      if( !getRowAggregationCandidates(aggrdata, probvaridx, &candrows, &candrowcoefs, &nrows, &ngoodrows) )
      {
         SCIPABORT();
      }
      assert(ngoodrows > 0); /* bounddistance was negative for this variable, so it should have good rows */

      for( k = 0; k < ngoodrows; ++k )
      {
         SCIP_Real rowaggrfac;
         int lppos;

         /* do not add rows twice */
         if( SCIPaggrRowHasRowBeenAdded(aggrrow, candrows[k]) )
            continue;

         rowaggrfac = - SCIPaggrRowGetProbvarValue(aggrrow, probvaridx) / candrowcoefs[k];

         /* if factor is too extreme skip this row */
         if( SCIPisFeasZero(scip, rowaggrfac) || REALABS(rowaggrfac) > sepadata->maxrowfac )
            continue;

         lppos = SCIProwGetLPPos(candrows[k]);
         /* row could be used and good rows are equalities, so ignore sidetype */
         {
            SCIP_Real rowscore = MAX(rowlhsscores[lppos], rowrhsscores[lppos]);

            /* if this rows score is better than the currently best score, remember it */
            if( aggrfac == 0.0 || rowscore > bestrowscore )
            {
               bestrow = candrows[k];
               aggrfac = rowaggrfac;
               bestrowscore = rowscore;
               bestrowside = 0;
            }
         }
      }
   }

   /* found a row among the good rows, so aggregate it and stop */
   if( aggrfac != 0.0 )
   {
      ++(*naggrs);
      SCIP_CALL( SCIPaggrRowAddRow(scip, aggrrow, bestrow, aggrfac, bestrowside) );
      SCIPaggrRowRemoveZeros(scip, aggrrow, success);
      goto TERMINATE;
   }

   for( i = 0; i < nbadvars; ++i )
   {
      int probvaridx;
      SCIP_ROW** candrows;
      SCIP_Real* candrowcoefs;
      int nrows;
      int ngoodrows;
      int k;

      /* if the bound distance is negative, there are no more variables to be tested, so stop */
      if( badvarbddist[i] < 0.0 )
         break;

      /* if no best row was found yet, this variable has the currently best bound distance */
      if( aggrfac == 0.0 )
         minbddist = badvarbddist[i] * (1.0 - sepadata->aggrtol);

      /* if the bound distance of the current variable is smaller than the minimum bound distance stop looping */
      if( badvarbddist[i] < minbddist )
         break;

      probvaridx = badvarinds[i];

      if( !getRowAggregationCandidates(aggrdata, probvaridx, &candrows, &candrowcoefs, &nrows, &ngoodrows) )
      {
         SCIPABORT();
      }

      /* bounddistance was positive for this variable, so it should not have good rows */
      assert(ngoodrows == 0);

      for( k = 0; k < nrows; ++k )
      {
         SCIP_Real rowaggrfac;
         int lppos;

         /* do not add rows twice */
         if( SCIPaggrRowHasRowBeenAdded(aggrrow, candrows[k]) )
            continue;

         rowaggrfac = - SCIPaggrRowGetProbvarValue(aggrrow, probvaridx) / candrowcoefs[k];

         /* if factor is too extreme skip this row */
         if( SCIPisFeasZero(scip, rowaggrfac) || REALABS(rowaggrfac) > sepadata->maxrowfac )
            continue;

         lppos = SCIProwGetLPPos(candrows[k]);

         /* row could be used, decide which side */
         {
            SCIP_Real rowscore;
            int rowside;

            /* either both or none of the rowscores are 0.0 so use the one which gives a positive slack */
            if( (rowaggrfac < 0.0 && !SCIPisInfinity(scip, -SCIProwGetLhs(candrows[k]))) ||
                  SCIPisInfinity(scip, SCIProwGetRhs(candrows[k])) )
            {
               rowscore = rowlhsscores[lppos];
               rowside = -1;
            }
            else
            {
               rowscore = rowrhsscores[lppos];
               rowside = 1;
            }

            /* if this rows score is better than the currently best score, remember it */
            if( aggrfac == 0.0 || SCIPisGT(scip, rowscore, bestrowscore) ||
                (SCIPisEQ(scip, rowscore, bestrowscore) && aggrdata->nbadvarsinrow[lppos] < aggrdata->nbadvarsinrow[SCIProwGetLPPos(bestrow)]) )
            {
               bestrow = candrows[k];
               aggrfac = rowaggrfac;
               bestrowscore = rowscore;
               bestrowside = rowside;
            }
         }
      }
   }

    /* found a row so aggregate it */
   if( aggrfac != 0.0 )
   {
      ++(*naggrs);
      SCIP_CALL( SCIPaggrRowAddRow(scip, aggrrow, bestrow, aggrfac, bestrowside) );
      SCIPaggrRowRemoveZeros(scip, aggrrow, success);
   }

TERMINATE:
   SCIPfreeBufferArray(scip, &badvarbddist);
   SCIPfreeBufferArray(scip, &badvarinds);

   return SCIP_OKAY;
}

/** aggregates different single mixed integer constraints by taking linear combinations of the rows of the LP  */
static
SCIP_RETCODE aggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   AGGREGATIONDATA*      aggrdata,           /**< pointer to aggregation data */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             allowlocal,         /**< should local cuts be allowed */
   SCIP_Real*            rowlhsscores,       /**< aggregation scores for left hand sides of row */
   SCIP_Real*            rowrhsscores,       /**< aggregation scores for right hand sides of row */
   int                   startrow,           /**< index of row to start aggregation */
   int                   maxaggrs,           /**< maximal number of aggregations */
   SCIP_Bool*            wastried,           /**< pointer to store whether the given startrow was actually tried */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  cutinds,            /**< buffer array to store temporarily cut */
   SCIP_Real*            cutcoefs,           /**< buffer array to store temporarily cut */
   SCIP_Bool             negate,             /**< should the start row be multiplied by -1 */
   int*                  ncuts               /**< pointer to count the number of generated cuts */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_ROW** rows;

   SCIP_Real startweight;
   SCIP_Real startrowact;
   int maxaggrnonzs;
   int naggrs;
   int nrows;
   int maxtestdelta;

   SCIP_Real cutrhs;
   SCIP_Real cutefficacy;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(rowlhsscores != NULL);
   assert(rowrhsscores != NULL);
   assert(wastried != NULL);
   assert(cutoff != NULL);
   assert(ncuts != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   *cutoff = FALSE;
   *wastried = FALSE;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert(nrows == 0 || rows != NULL);
   assert(0 <= startrow && startrow < nrows);

   SCIPdebugMsg(scip, "start c-MIR aggregation with row <%s> (%d/%d)\n", SCIProwGetName(rows[startrow]), startrow, nrows);

   /* calculate maximal number of non-zeros in aggregated row */
   maxaggrnonzs = (int)(sepadata->maxaggdensity * SCIPgetNLPCols(scip)) + sepadata->densityoffset;

   startrowact = SCIPgetRowSolActivity(scip, rows[startrow], sol);

   if( startrowact <= 0.5 * SCIProwGetLhs(rows[startrow]) + 0.5 * SCIProwGetRhs(rows[startrow]) )
      startweight = -1.0;
   else
      startweight = 1.0;

   maxtestdelta = sepadata->maxtestdelta == -1 ? INT_MAX : sepadata->maxtestdelta;

   /* add start row to the initially empty aggregation row (aggrrow) */
   SCIP_CALL( SCIPaggrRowAddRow(scip, aggrdata->aggrrow, rows[startrow], negate ? -startweight : startweight, 0) );

   /* try to generate cut from the current aggregated row; add cut if found, otherwise add another row to aggrrow
    * in order to get rid of a continuous variable
    */
   naggrs = 0;
   while( naggrs <= maxaggrs )
   {
      int cutrank;
      int cutnnz;
      SCIP_Bool aggrsuccess;
      SCIP_Bool cmirsuccess;
      SCIP_Bool cmircutislocal;
      SCIP_Bool flowcoversuccess;
      SCIP_Real flowcoverefficacy;
      SCIP_Bool flowcovercutislocal;
      SCIP_ROW* cut;

      *wastried = TRUE;

      /* Step 1:
       * try to generate a MIR cut out of the current aggregated row
       */

      flowcoverefficacy =  -SCIPinfinity(scip);
      SCIP_CALL( SCIPcalcFlowCover(scip, sol, POSTPROCESS, BOUNDSWITCH, allowlocal, aggrdata->aggrrow,
         cutcoefs, &cutrhs, cutinds, &cutnnz, &flowcoverefficacy, &cutrank, &flowcovercutislocal, &flowcoversuccess) );

      cutefficacy = flowcoverefficacy;
      SCIP_CALL( SCIPcutGenerationHeuristicCMIR(scip, sol, POSTPROCESS, BOUNDSWITCH, USEVBDS, allowlocal, maxtestdelta, NULL, NULL, MINFRAC, MAXFRAC,
         aggrdata->aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, &cutrank, &cmircutislocal, &cmirsuccess) );

      cut = NULL;

      if( cmirsuccess )
      {
         SCIP_CALL( addCut(scip, sol, sepadata->cmir, FALSE, cutcoefs, cutinds, cutnnz, cutrhs, cutefficacy, cmircutislocal,
               sepadata->dynamiccuts, cutrank, "cmir", cutoff, ncuts, &cut) );
      }
      else if ( flowcoversuccess )
      {
         SCIP_CALL( addCut(scip, sol, sepadata->flowcover, FALSE, cutcoefs, cutinds, cutnnz, cutrhs, cutefficacy, flowcovercutislocal,
               sepadata->dynamiccuts, cutrank, "flowcover", cutoff, ncuts, &cut) );
      }

      if ( *cutoff )
      {
         if( cut != NULL )
         {
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
         break;
      }

      /* if the cut was successfully added, decrease the score of the rows used in the aggregation and clean the aggregation
       * row (and call this function again with a different start row for aggregation)
       */
      if( cut != NULL )
      {
         int* rowinds;
         int i;

         rowinds = SCIPaggrRowGetRowInds(aggrdata->aggrrow);
         nrows = SCIPaggrRowGetNRows(aggrdata->aggrrow);

         /* decrease row score of used rows slightly */
         for( i = 0; i < nrows; ++i )
         {
            SCIP_Real fac = 1.0 - 0.999 * SCIProwGetParallelism(rows[rowinds[i]], cut, 'e');

            rowlhsscores[rowinds[i]] *= fac;
            rowrhsscores[rowinds[i]] *= fac;
         }

         SCIP_CALL( SCIPreleaseRow(scip, &cut) );

         SCIPdebugMsg(scip, " -> abort aggregation: cut found\n");
         break;
      }

      /* Step 2:
       * aggregate an additional row in order to remove a continuous variable
       */

      /* abort, if we reached the maximal number of aggregations */
      if( naggrs == maxaggrs )
      {
         SCIPdebugMsg(scip, " -> abort aggregation: maximal number of aggregations reached\n");
         break;
      }

      SCIP_CALL( aggregateNextRow(scip, sepadata, rowlhsscores, rowrhsscores, aggrdata, aggrdata->aggrrow,
            &naggrs, &aggrsuccess) );

      /* no suitable aggregation was found or number of non-zeros is now too large so abort */
      if( ! aggrsuccess || SCIPaggrRowGetNNz(aggrdata->aggrrow) > maxaggrnonzs || SCIPaggrRowGetNNz(aggrdata->aggrrow) == 0 )
      {
         break;
      }

      SCIPdebugMsg(scip, " ->  current aggregation has %d/%d nonzeros and consists of %d/%d rows\n",
          SCIPaggrRowGetNNz(aggrdata->aggrrow), maxaggrnonzs, naggrs, maxaggrs);
   }

   SCIPaggrRowClear(aggrdata->aggrrow);

   return SCIP_OKAY;
}

/** gives an estimate of how much the activity of this row is
 *  affected by fractionality in the current solution
 */
static
SCIP_Real getRowFracActivity(
   SCIP_ROW*             row,                /**< the LP row */
   SCIP_Real*            fractionalities     /**< array of fractionalities for each variable */
   )
{
   int nlpnonz;
   int i;
   SCIP_COL** cols;
   SCIP_Real* vals;
   SCIP_Real fracsum = 0.0;

   cols = SCIProwGetCols(row);
   vals = SCIProwGetVals(row);
   nlpnonz = SCIProwGetNLPNonz(row);

   for( i = 0; i < nlpnonz; ++i )
   {
      SCIP_VAR* var = SCIPcolGetVar(cols[i]);
      fracsum += REALABS(vals[i] * fractionalities[SCIPvarGetProbindex(var)]);
   }

   return fracsum;


}

/** searches and adds c-MIR cuts that separate the given primal solution */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the c-MIR separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             allowlocal,         /**< should local cuts be allowed */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   AGGREGATIONDATA aggrdata;
   SCIP_SEPADATA* sepadata;
   SCIP_VAR** vars;
   SCIP_Real* varsolvals;
   SCIP_Real* bestcontlbs;
   SCIP_Real* bestcontubs;
   SCIP_Real* fractionalities;
   SCIP_ROW** rows;
   SCIP_Real* rowlhsscores;
   SCIP_Real* rowrhsscores;
   SCIP_Real* rowscores;
   int* roworder;
   SCIP_Real maxslack;
   SCIP_Bool cutoff = FALSE;
   int nvars;
   int nintvars;
   int ncontvars;
   int nrows;
   int nnonzrows;
   int zerorows;
   int ntries;
   int nfails;
   int depth;
   int ncalls;
   int maxtries;
   int maxfails;
   int maxaggrs;
   int maxsepacuts;
   int ncuts;
   int r;
   int v;

   int* cutinds;
   SCIP_Real* cutcoefs;

   assert(result != NULL);
   assert(*result == SCIP_DIDNOTRUN);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call the cmir cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* get all rows and number of columns */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert(nrows == 0 || rows != NULL);

   /* nothing to do, if LP is empty */
   if( nrows == 0 )
      return SCIP_OKAY;

   /* check whether SCIP was stopped in the meantime */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* get active problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   ncontvars = SCIPgetNContVars(scip);
#ifdef IMPLINTSARECONT
   ncontvars += SCIPgetNImplVars(scip); /* also aggregate out implicit integers */
#endif
   nintvars = nvars - ncontvars;
   assert(nvars == 0 || vars != NULL);

   /* nothing to do, if problem has no variables */
   if( nvars == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "separating c-MIR cuts\n");

   *result = SCIP_DIDNOTFIND;

   /* get data structure */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowlhsscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowrhsscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowscores, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &roworder, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestcontlbs, ncontvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestcontubs, ncontvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fractionalities, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );

   /* get the solution values for all active variables */
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, varsolvals) );

   /* calculate the fractionality of the integer variables in the current solution */
   for( v = 0; v < nintvars; ++v )
   {
      fractionalities[v] = SCIPfeasFrac(scip, varsolvals[v]);
      fractionalities[v] = MIN(fractionalities[v], 1.0 - fractionalities[v]);
   }

   /* calculate the fractionality of the continuous variables in the current solution;
    * The fractionality of a continuous variable x is defined to be a * f_y,
    * if there is a variable bound x <= a * y + c where f_y is the fractionality of y
    * and in the current solution the variable bound has no slack.
    */
   for( ; v < nvars; ++v )
   {
      SCIP_VAR** vlbvars;
      SCIP_VAR** vubvars;
      SCIP_Real* vlbcoefs;
      SCIP_Real* vubcoefs;
      SCIP_Real closestvlb;
      SCIP_Real closestvub;
      int closestvlbidx;
      int closestvubidx;

      SCIP_CALL( SCIPgetVarClosestVlb(scip, vars[v], sol, &closestvlb, &closestvlbidx) );
      SCIP_CALL( SCIPgetVarClosestVub(scip, vars[v], sol, &closestvub, &closestvubidx) );

      vlbvars = SCIPvarGetVlbVars(vars[v]);
      vubvars = SCIPvarGetVubVars(vars[v]);
      vlbcoefs = SCIPvarGetVlbCoefs(vars[v]);
      vubcoefs = SCIPvarGetVubCoefs(vars[v]);

      fractionalities[v] = 0.0;
      if( closestvlbidx != -1 && SCIPisEQ(scip, varsolvals[v], closestvlb) )
      {
         int vlbvarprobidx = SCIPvarGetProbindex(vlbvars[closestvlbidx]);
         SCIP_Real frac = SCIPfeasFrac(scip, varsolvals[vlbvarprobidx]);
         if( frac < 0.0 )
            frac = 0.0;
         assert(frac >= 0.0 && frac < 1.0);
         frac = MIN(frac, 1.0 - frac) * vlbcoefs[closestvlbidx];
         fractionalities[v] += frac;
      }

      if( closestvubidx != -1 && SCIPisEQ(scip, varsolvals[v], closestvub) )
      {
         int vubvarprobidx = SCIPvarGetProbindex(vubvars[closestvubidx]);
         SCIP_Real frac = SCIPfeasFrac(scip, varsolvals[vubvarprobidx]);
         if( frac < 0.0 )
            frac = 0.0;
         assert(frac >= 0.0 && frac < 1.0);
         frac = MIN(frac, 1.0 - frac) * vubcoefs[closestvubidx];
         fractionalities[v] += frac;
      }
   }

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
   {
      maxtries = sepadata->maxtriesroot;
      maxfails = sepadata->maxfailsroot;
      maxaggrs = sepadata->maxaggrsroot;
      maxsepacuts = sepadata->maxsepacutsroot;
      maxslack = sepadata->maxslackroot;
   }
   else
   {
      maxtries = sepadata->maxtries;
      maxfails = sepadata->maxfails;
      maxaggrs = sepadata->maxaggrs;
      maxsepacuts = sepadata->maxsepacuts;
      maxslack = sepadata->maxslack;
   }

   /* calculate aggregation scores for both sides of all rows, and sort rows by decreasing maximal score
    * TODO: document score definition */

   /* count the number of non-zero rows and zero rows. these values are used for the sorting of the rowscores.
    * only the non-zero rows need to be sorted. */
   nnonzrows = 0;
   zerorows = 0;
   for( r = 0; r < nrows; r++ )
   {
      int nnonz;
      int i;

      assert(SCIProwGetLPPos(rows[r]) == r);

      nnonz = SCIProwGetNLPNonz(rows[r]);
      if( nnonz == 0 )
      {
         /* ignore empty rows */
         rowlhsscores[r] = 0.0;
         rowrhsscores[r] = 0.0;

         /* add the row number to the back of roworder for zero rows */
         zerorows++;
         rowscores[r] = 0.0;
         roworder[nrows - zerorows] = r;
      }
      else
      {
         SCIP_Real activity;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real dualsol;
         SCIP_Real dualscore;
         SCIP_Real rowdensity;
         SCIP_Real rownorm;
         SCIP_Real slack;
         SCIP_Real fracact;
         SCIP_Real fracscore;
         SCIP_Real objnorm;

         objnorm = SCIPgetObjNorm(scip);
         objnorm = MAX(objnorm, 1.0);

         fracact = getRowFracActivity(rows[r], fractionalities);
         dualsol = (sol == NULL ? SCIProwGetDualsol(rows[r]) : 1.0);
         activity = SCIPgetRowSolActivity(scip, rows[r], sol);
         lhs = SCIProwGetLhs(rows[r]);
         rhs = SCIProwGetRhs(rows[r]);
         rownorm = SCIProwGetNorm(rows[r]);
         rownorm = MAX(rownorm, 0.1);
         rowdensity = (SCIP_Real)(nnonz - sepadata->densityoffset)/(SCIP_Real)nvars;
         assert(SCIPisPositive(scip, rownorm));
         fracscore = fracact / rownorm;

         slack = (activity - lhs)/rownorm;
         dualscore = MAX(fracscore * dualsol/objnorm, 0.0001);
         if( !SCIPisInfinity(scip, -lhs) && SCIPisLE(scip, slack, maxslack)
            && (allowlocal || !SCIProwIsLocal(rows[r])) /*lint !e506 !e774*/
            && rowdensity <= sepadata->maxrowdensity
            && rowdensity <= sepadata->maxaggdensity )  /*lint !e774*/
         {
            rowlhsscores[r] = dualscore + sepadata->densityscore * (1.0-rowdensity)
               + sepadata->slackscore * MAX(1.0 - slack, 0.0);
            assert(rowlhsscores[r] > 0.0);
         }
         else
            rowlhsscores[r] = 0.0;

         slack = (rhs - activity)/rownorm;
         dualscore = MAX(-fracscore * dualsol/objnorm, 0.0001);
         if( !SCIPisInfinity(scip, rhs) && SCIPisLE(scip, slack, maxslack)
            && (allowlocal || !SCIProwIsLocal(rows[r])) /*lint !e506 !e774*/
            && rowdensity <= sepadata->maxrowdensity
            && rowdensity <= sepadata->maxaggdensity )  /*lint !e774*/
         {
            rowrhsscores[r] = dualscore + sepadata->densityscore * (1.0-rowdensity)
               + sepadata->slackscore * MAX(1.0 - slack, 0.0);
            assert(rowrhsscores[r] > 0.0);
         }
         else
            rowrhsscores[r] = 0.0;

         /* for the row order only use the fractionality score since it best indicates how likely it is to find a cut */
         rowscores[r] = fracscore;
         if( rowscores[r] == 0.0 )
         {
            /* add the row number to the back of roworder for zero rows */
            zerorows++;
            roworder[nrows - zerorows] = r;
         }
         else
         {
            /* insert the row number in the correct position of roworder */
            for( i = nnonzrows; i > 0 && rowscores[r] > rowscores[roworder[i - 1]]; --i )
               roworder[i] = roworder[i - 1];
            roworder[i] = r;

            nnonzrows++;
         }
      }

      SCIPdebugMsg(scip, " -> row %d <%s>: lhsscore=%g rhsscore=%g maxscore=%g\n", r, SCIProwGetName(rows[r]),
            rowlhsscores[r], rowrhsscores[r], rowscores[r]);
   }
   assert(nrows == nnonzrows + zerorows);

   /* calculate the data required for performing the row aggregation */
   SCIP_CALL( setupAggregationData(scip, sol, allowlocal, &aggrdata) );

   ncuts = 0;
   if( maxtries < 0 )
      maxtries = INT_MAX;
   if( maxfails < 0 )
      maxfails = INT_MAX;
   else if( depth == 0 && 2 * SCIPgetNSepaRounds(scip) < maxfails )
      maxfails += maxfails - 2 * SCIPgetNSepaRounds(scip); /* allow up to double as many fails in early separounds of root node */

   /* start aggregation heuristic for each row in the LP and generate resulting cuts */
   ntries = 0;
   nfails = 0;
   for( r = 0; r < nnonzrows && ntries < maxtries && ncuts < maxsepacuts && !SCIPisStopped(scip); r++ )
   {
      SCIP_Bool wastried;
      int oldncuts;

      oldncuts = ncuts;
      SCIP_CALL( aggregation(scip, &aggrdata, sepa, sol, allowlocal, rowlhsscores, rowrhsscores,
                             roworder[r], maxaggrs, &wastried, &cutoff, cutinds, cutcoefs, FALSE, &ncuts) );

      /* if trynegscaling is true we start the aggregation heuristic again for this row, but multiply it by -1 first.
       * This is done by calling the aggregation function with the parameter negate equal to TRUE
       */
      if( sepadata->trynegscaling && !cutoff )
      {
         SCIP_CALL( aggregation(scip, &aggrdata, sepa, sol, allowlocal, rowlhsscores, rowrhsscores,
                             roworder[r], maxaggrs, &wastried, &cutoff, cutinds, cutcoefs, TRUE, &ncuts) );
      }

      if ( cutoff )
         break;

      if( !wastried )
      {
         continue;
      }
      ntries++;

      if( ncuts == oldncuts )
      {
         nfails++;
         if( nfails >= maxfails )
         {
            break;
         }
      }
      else
      {
         nfails = 0;
      }
   }

   /* free data structure */
   destroyAggregationData(scip, &aggrdata);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &cutinds);
   SCIPfreeBufferArray(scip, &fractionalities);
   SCIPfreeBufferArray(scip, &bestcontubs);
   SCIPfreeBufferArray(scip, &bestcontlbs);
   SCIPfreeBufferArray(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &roworder);
   SCIPfreeBufferArray(scip, &rowscores);
   SCIPfreeBufferArray(scip, &rowrhsscores);
   SCIPfreeBufferArray(scip, &rowlhsscores);

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyAggregation)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaAggregation(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeAggregation)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpAggregation)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( separateCuts(scip, sepa, NULL, allowlocal, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolAggregation)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( separateCuts(scip, sepa, sol, allowlocal, result) );

   return SCIP_OKAY;
}

/** LP solution separation method of dummy separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpDummy)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of dummy separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolDummy)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the cmir separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaAggregation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create cmir separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );

   /* include dummy separators */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepadata->flowcover, "flowcover", "dummy separator for adding flowcover cuts", -100000, -1, SEPA_MAXBOUNDDIST,
                                   SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpDummy, sepaExecsolDummy, NULL) );

   assert(sepadata->flowcover != NULL);

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepadata->cmir, "cmir", "dummy separator for adding cmir cuts", -100000, -1, SEPA_MAXBOUNDDIST,
                                   SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpDummy, sepaExecsolDummy, NULL) );

   assert(sepadata->cmir != NULL);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpAggregation, sepaExecsolAggregation,
         sepadata) );

   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyAggregation) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeAggregation) );

   /* add cmir separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxrounds",
         "maximal number of cmir separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxroundsroot",
         "maximal number of cmir separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxtries",
         "maximal number of rows to start aggregation with per separation round (-1: unlimited)",
         &sepadata->maxtries, TRUE, DEFAULT_MAXTRIES, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxtriesroot",
         "maximal number of rows to start aggregation with per separation round in the root node (-1: unlimited)",
         &sepadata->maxtriesroot, TRUE, DEFAULT_MAXTRIESROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxfails",
         "maximal number of consecutive unsuccessful aggregation tries (-1: unlimited)",
         &sepadata->maxfails, TRUE, DEFAULT_MAXFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxfailsroot",
         "maximal number of consecutive unsuccessful aggregation tries in the root node (-1: unlimited)",
         &sepadata->maxfailsroot, TRUE, DEFAULT_MAXFAILSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxaggrs",
         "maximal number of aggregations for each row per separation round",
         &sepadata->maxaggrs, TRUE, DEFAULT_MAXAGGRS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxaggrsroot",
         "maximal number of aggregations for each row per separation round in the root node",
         &sepadata->maxaggrsroot, TRUE, DEFAULT_MAXAGGRSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsepacuts",
         "maximal number of cmir cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsepacutsroot",
         "maximal number of cmir cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/maxslack",
         "maximal slack of rows to be used in aggregation",
         &sepadata->maxslack, TRUE, DEFAULT_MAXSLACK, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/maxslackroot",
         "maximal slack of rows to be used in aggregation in the root node",
         &sepadata->maxslackroot, TRUE, DEFAULT_MAXSLACKROOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/densityscore",
         "weight of row density in the aggregation scoring of the rows",
         &sepadata->densityscore, TRUE, DEFAULT_DENSITYSCORE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/slackscore",
         "weight of slack in the aggregation scoring of the rows",
         &sepadata->slackscore, TRUE, DEFAULT_SLACKSCORE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/maxaggdensity",
         "maximal density of aggregated row",
         &sepadata->maxaggdensity, TRUE, DEFAULT_MAXAGGDENSITY, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/maxrowdensity",
         "maximal density of row to be used in aggregation",
         &sepadata->maxrowdensity, TRUE, DEFAULT_MAXROWDENSITY, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/densityoffset",
         "additional number of variables allowed in row on top of density",
         &sepadata->densityoffset, TRUE, DEFAULT_DENSITYOFFSET, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/maxrowfac",
         "maximal row aggregation factor",
         &sepadata->maxrowfac, TRUE, DEFAULT_MAXROWFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxtestdelta",
         "maximal number of different deltas to try (-1: unlimited)",
         &sepadata->maxtestdelta, TRUE, DEFAULT_MAXTESTDELTA, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/aggrtol",
         "tolerance for bound distances used to select continuous variable in current aggregated constraint to be eliminated",
         &sepadata->aggrtol, TRUE, DEFAULT_AGGRTOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/trynegscaling",
         "should negative values also be tested in scaling?",
         &sepadata->trynegscaling, TRUE, DEFAULT_TRYNEGSCALING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/fixintegralrhs",
         "should an additional variable be complemented if f0 = 0?",
         &sepadata->fixintegralrhs, TRUE, DEFAULT_FIXINTEGRALRHS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
