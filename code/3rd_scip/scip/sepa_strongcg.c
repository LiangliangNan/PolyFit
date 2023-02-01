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

/**@file   sepa_strongcg.c
 * @brief  Strong CG Cuts (Letchford & Lodi)
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_strongcg.h"
#include "scip/pub_misc.h"


#define SEPA_NAME              "strongcg"
#define SEPA_DESC              "Strong CG cuts separator (Letchford and Lodi)"
#define SEPA_PRIORITY             -2000
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of strong CG separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        20 /**< maximal number of strong CG separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          20 /**< maximal number of strong CG cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     500 /**< maximal number of strong CG cuts separated per separation round in root node */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_RANDSEED             54 /**< initial random seed */

#define SEPARATEROWS           /* separate rows with integral slack */

#define BOUNDSWITCH              0.9999
#define POSTPROCESS                TRUE
#define USEVBDS                    TRUE
#define MINFRAC                    0.05
#define MAXFRAC                    0.95

#define MAXAGGRLEN(nvars)          (0.1*(nvars)+1000) /**< maximal length of base inequality */


/** separator data */
struct SCIP_SepaData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   int                   maxrounds;          /**< maximal number of strong CG separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of strong CG separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of strong CG cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of strong CG cuts separated per separation round in root node */
   int                   lastncutsfound;     /**< total number of cuts found after last call of separator */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
};

/*
 * Callback methods
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyStrongcg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaStrongcg(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeStrongcg)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitStrongcg)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* create and initialize random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &sepadata->randnumgen, DEFAULT_RANDSEED) );

   return SCIP_OKAY;
}

/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitStrongcg)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeRandom(scip, &sepadata->randnumgen);

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpStrongcg)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_VAR** vars;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_AGGRROW* aggrrow;
   SCIP_Real* varsolvals;
   SCIP_Real* binvrow;
   SCIP_Real* cutcoefs;
   SCIP_Real* basisfrac;
   SCIP_Real cutrhs;
   int* basisind;
   int* basisperm;
   int* inds;
   int* cutinds;
   int cutnnz;
   int ninds;
   int nvars;
   int ncols;
   int nrows;
   int ncalls;
   int depth;
   int maxsepacuts;
   int ncuts;
   int c;
   int i;
   int cutrank;
   SCIP_Bool success;
   SCIP_Bool cutislocal;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) || !allowlocal )
      return SCIP_OKAY;

   /* only call the strong CG cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if the LP solution is basic */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* get variables data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get LP data */
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   if( ncols == 0 || nrows == 0 )
      return SCIP_OKAY;

#if 0 /* if too many columns, separator is usually very slow: delay it until no other cuts have been found */
   if( ncols >= 50*nrows )
      return SCIP_OKAY;
   if( ncols >= 5*nrows )
   {
      int ncutsfound;

      ncutsfound = SCIPgetNCutsFound(scip);
      if( ncutsfound > sepadata->lastncutsfound || !SCIPsepaWasLPDelayed(sepa) )
      {
         sepadata->lastncutsfound = ncutsfound;
         *result = SCIP_DELAYED;
         return SCIP_OKAY;
      }
   }
#endif

   *result = SCIP_DIDNOTFIND;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisperm, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisfrac, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nrows) );
   SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );

   varsolvals = NULL; /* allocate this later, if needed */

   /* get basis indices */
   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );

   for( i = 0; i < nrows; ++i )
   {
      SCIP_Real frac = 0.0;

      c = basisind[i];

      basisperm[i] = i;

      if( c >= 0 )
      {
         SCIP_VAR* var;

         assert(c < ncols);
         var = SCIPcolGetVar(cols[c]);
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            frac = SCIPfeasFrac(scip, SCIPcolGetPrimsol(cols[c]));
            frac = MIN(frac, 1.0 - frac);
         }
      }
#ifdef SEPARATEROWS
      else
      {
         SCIP_ROW* row;

         assert(0 <= -c-1 && -c-1 < nrows);
         row = rows[-c-1];
         if( SCIProwIsIntegral(row) && !SCIProwIsModifiable(row) )
         {
            frac = SCIPfeasFrac(scip, SCIPgetRowActivity(scip, row));
            frac = MIN(frac, 1.0 - frac);
         }
      }
#endif

      if( frac >= MINFRAC )
      {
         /* slightly change fractionality to have random order for equal fractions */
         basisfrac[i] = frac + SCIPrandomGetReal(sepadata->randnumgen, -1e-6, 1e-6);
      }
      else
      {
         basisfrac[i] = 0.0;
      }
   }

   /* sort basis indices by fractionality */
   SCIPsortDownRealInt(basisfrac, basisperm, nrows);

   /* get the maximal number of cuts allowed in a separation round */
   if( depth == 0 )
      maxsepacuts = sepadata->maxsepacutsroot;
   else
      maxsepacuts = sepadata->maxsepacuts;

   SCIPdebugMsg(scip, "searching strong CG cuts: %d cols, %d rows, maxcuts=%d\n",
      ncols, nrows, maxsepacuts);

   /* for all basic columns belonging to integer variables, try to generate a strong CG cut */
   ncuts = 0;
   for( i = 0; i < nrows && ncuts < maxsepacuts && !SCIPisStopped(scip) && *result != SCIP_CUTOFF; ++i )
   {
      int j;
      SCIP_Real cutefficacy;

      if( basisfrac[i] == 0.0 )
         break;

      j = basisperm[i];
      c = basisind[j];

      /* get the row of B^-1 for this basic integer variable with fractional solution value */
      SCIP_CALL( SCIPgetLPBInvRow(scip, j, binvrow, inds, &ninds) );

#ifdef SCIP_DEBUG
      /* initialize variables, that might not have been initialized in SCIPcalcMIR if success == FALSE */
      cutefficacy = 0.0;
      cutrhs = SCIPinfinity(scip);
#endif

      /* create the aggregation row using the B^-1 row as weights */
      SCIP_CALL( SCIPaggrRowSumRows(scip, aggrrow, binvrow, inds, ninds,
         FALSE, allowlocal, 1, (int) MAXAGGRLEN(nvars), &success) );

      if( !success )
         continue;

      /* create a strong CG cut out of the aggregation row */
      SCIP_CALL( SCIPcalcStrongCG(scip, NULL, POSTPROCESS, BOUNDSWITCH, USEVBDS, allowlocal, MINFRAC, MAXFRAC,
         1.0, aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, &cutrank, &cutislocal, &success) );

      assert(allowlocal || !cutislocal);
      SCIPdebugMsg(scip, " -> success=%u: rhs: %g, efficacy: %g\n", success, cutrhs, cutefficacy);

      if( !success )
         continue;

      /* if successful, convert dense cut into sparse row, and add the row as a cut */
      if( SCIPisEfficacious(scip, cutefficacy) )
      {
         SCIP_ROW* cut;
         char cutname[SCIP_MAXSTRLEN];
         int v;

         SCIPdebugMsg(scip, " -> strong CG cut for <%s>: act=%f, rhs=%f, norm=%f, eff=%f, rank=%d\n",
                        c >= 0 ? SCIPvarGetName(SCIPcolGetVar(cols[c])) : SCIProwGetName(rows[-c-1]),
                        cutefficacy * SCIPgetVectorEfficacyNorm(scip, cutcoefs, cutnnz) + cutrhs, cutrhs, SCIPgetVectorEfficacyNorm(scip, cutcoefs, cutnnz), cutefficacy, cutrank);

         /* create the cut */
         if( c >= 0 )
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "scg%d_x%d", SCIPgetNLPs(scip), c);
         else
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "scg%d_s%d", SCIPgetNLPs(scip), -c-1);

         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, sepadata->dynamiccuts) );

         /*SCIPdebug( SCIP_CALL(SCIPprintRow(scip, cut, NULL)) );*/
         SCIProwChgRank(cut, cutrank);

         /* cache the row extension and only flush them if the cut gets added */
         SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

         /* collect all non-zero coefficients */
         for( v = 0; v < cutnnz; ++v )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[v]], cutcoefs[v]) );
         }

         assert(success);

         if( !SCIPisCutEfficacious(scip, NULL, cut) )
         {
            SCIPdebugMsg(scip, " -> strong CG cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
                           cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                           SCIPgetCutEfficacy(scip, NULL, cut));
            /*SCIPdebug( SCIP_CALL(SCIPprintRow(scip, cut, NULL)) );*/
            success = FALSE;
         }
         else
         {
            SCIP_Bool infeasible = FALSE;

            /* flush all changes before adding the cut */
            SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

            SCIPdebugMsg(scip, " -> found strong CG cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                           cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                           SCIPgetCutEfficacy(scip, NULL, cut),
                           SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                           SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
            /*SCIPdebug( SCIP_CALL(SCIPprintRow(scip, cut, NULL)) );*/

            if( SCIPisCutNew(scip, cut) )
            {
               if( !cutislocal )
               {
                  SCIP_CALL( SCIPaddPoolCut(scip, cut) );
               }
               else
               {
                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );
               }

               ncuts++;

               if( infeasible )
               {
                  *result = SCIP_CUTOFF;
               }
               else
               {
                  *result = SCIP_SEPARATED;
               }
            }
         }

         /* release the row */
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      }
   }

   /* free temporary memory */
   SCIPaggrRowFree(scip, &aggrrow);
   SCIPfreeBufferArrayNull(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basisfrac);
   SCIPfreeBufferArray(scip, &basisperm);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &cutinds);
   SCIPfreeBufferArray(scip, &cutcoefs);

   SCIPdebugMsg(scip, "end searching strong CG cuts: found %d cuts\n", ncuts);

   sepadata->lastncutsfound = SCIPgetNCutsFound(scip);

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the Strong CG cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaStrongcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepadata->lastncutsfound = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpStrongcg, NULL,
         sepadata) );

   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyStrongcg) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeStrongcg) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitStrongcg) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitStrongcg) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/strongcg/maxrounds",
         "maximal number of strong CG separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/strongcg/maxroundsroot",
         "maximal number of strong CG separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/strongcg/maxsepacuts",
         "maximal number of strong CG cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/strongcg/maxsepacutsroot",
         "maximal number of strong CG cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/strongcg/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

