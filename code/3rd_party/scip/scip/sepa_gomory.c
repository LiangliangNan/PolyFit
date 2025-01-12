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

/**@file   sepa_gomory.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  Gomory MIR Cuts
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Domenico Salvagnin
 * @author Marc Pfetsch
 * @author Leona Gottwald
 */

/**@todo try k-Gomory-cuts (s. Cornuejols: K-Cuts: A Variation of Gomory Mixed Integer Cuts from the LP Tableau)
 *
 * @todo Try cuts on the objective tableau row.
 *
 * @todo Also try negative basis inverse row?
 *
 * @todo It happens that the SCIPcalcMIR() function returns with the same cut for different calls. Check if this is a
 *       bug or do not use it for the MIP below and turn off presolving and all heuristics:
 *
 *  Max y
 *  Subject to
 *  c1: -x + y <= 1
 *  c2: 2x + 3y <= 12
 *  c3: 3x + 2y <= 12
 *  Bounds
 *  0 <= x
 *  0 <= y
 *  General
 *  x
 *  y
 *  END
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cuts.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_sepa.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sepa.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/sepa_gomory.h"
#include <string.h>

#define SEPA_NAME              "gomory"
#define SEPA_DESC              "separator for Gomory mixed-integer and strong CG cuts from LP tableau rows"
#define SEPA_PRIORITY             -1000
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of gomory separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        10 /**< maximal number of gomory separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of gomory cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of gomory cuts separated per separation round in root node */
#define DEFAULT_MAXRANK              -1 /**< maximal rank of a gomory cut that could not be scaled to integral coefficients (-1: unlimited) */
#define DEFAULT_MAXRANKINTEGRAL      -1 /**< maximal rank of a gomory cut that could be scaled to integral coefficients (-1: unlimited) */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_AWAY               0.01 /**< minimal integrality violation of a basis variable in order to try Gomory cut */
#define DEFAULT_MAKEINTEGRAL      FALSE /**< try to scale all cuts to integral coefficients */
#define DEFAULT_FORCECUTS          TRUE /**< if conversion to integral coefficients failed still consider the cut */
#define DEFAULT_SEPARATEROWS       TRUE /**< separate rows with integral slack */
#define DEFAULT_DELAYEDCUTS       FALSE /**< should cuts be added to the delayed cut pool? */
#define DEFAULT_SIDETYPEBASIS      TRUE /**< choose side types of row (lhs/rhs) based on basis information? */
#define DEFAULT_TRYSTRONGCG        TRUE /**< try to generate strengthened Chvatal-Gomory cuts? */
#define DEFAULT_GENBOTHGOMSCG      TRUE /**< should both Gomory and strong CG cuts be generated (otherwise take best) */
#define DEFAULT_RANDSEED             53 /**< initial random seed */

#define BOUNDSWITCH              0.9999 /**< threshold for bound switching - see SCIPcalcMIR() */
#define POSTPROCESS                TRUE /**< apply postprocessing after MIR calculation - see SCIPcalcMIR() */
#define USEVBDS                    TRUE /**< use variable bounds - see SCIPcalcMIR() */
#define FIXINTEGRALRHS            FALSE /**< try to generate an integral rhs - see SCIPcalcMIR() */
#define MAKECONTINTEGRAL          FALSE /**< convert continuous variable to integral variables in SCIPmakeRowIntegral() */

#define MAXAGGRLEN(nvars)          (0.1*(nvars)+1000) /**< maximal length of base inequality */


/** separator data */
struct SCIP_SepaData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_SEPA*            strongcg;           /**< strong CG cut separator */
   SCIP_SEPA*            gomory;             /**< gomory cut separator */
   SCIP_Real             away;               /**< minimal integrality violation of a basis variable in order to try Gomory cut */
   int                   maxrounds;          /**< maximal number of gomory separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of gomory separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of gomory cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of gomory cuts separated per separation round in root node */
   int                   maxrank;            /**< maximal rank of a gomory cut that could not be scaled to integral coefficients (-1: unlimited) */
   int                   maxrankintegral;    /**< maximal rank of a gomory cut that could be scaled to integral coefficients (-1: unlimited) */
   int                   lastncutsfound;     /**< total number of cuts found after last call of separator */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_Bool             makeintegral;       /**< try to scale all cuts to integral coefficients */
   SCIP_Bool             forcecuts;          /**< if conversion to integral coefficients failed still consider the cut */
   SCIP_Bool             separaterows;       /**< separate rows with integral slack */
   SCIP_Bool             delayedcuts;        /**< should cuts be added to the delayed cut pool? */
   SCIP_Bool             sidetypebasis;      /**< choose side types of row (lhs/rhs) based on basis information? */
   SCIP_Bool             trystrongcg;        /**< try to generate strengthened Chvatal-Gomory cuts? */
   SCIP_Bool             genbothgomscg;      /**< should both Gomory and strong CG cuts be generated (otherwise take best) */
};


/** returns TRUE if the cut can be taken, otherwise FALSE if there some numerical evidences */
static
SCIP_RETCODE evaluateCutNumerics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< data of the separator */
   SCIP_ROW*             cut,                /**< cut to check */
   SCIP_Longint          maxdnom,            /**< maximal denominator to use for scaling */
   SCIP_Real             maxscale,           /**< maximal scaling factor */
   SCIP_Bool*            useful              /**< pointer to store if the cut is useful */
   )
{
   SCIP_Bool madeintegral = FALSE;

   assert(useful != NULL);

   *useful = FALSE;

   if( sepadata->makeintegral && SCIPgetRowNumIntCols(scip, cut) == SCIProwGetNNonz(cut) )
   {
      /* try to scale the cut to integral values */
      SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
            maxdnom, maxscale, MAKECONTINTEGRAL, &madeintegral) );

      if( !madeintegral && !sepadata->forcecuts )
         return SCIP_OKAY;

      /* in case the right hand side is plus infinity (due to scaling) the cut is useless so we are not taking it at all */
      if( madeintegral && SCIPisInfinity(scip, SCIProwGetRhs(cut)) )
         return SCIP_OKAY;
   }

   /* discard integral cut if the rank is too high */
   if( madeintegral && sepadata->maxrankintegral != -1 && (SCIProwGetRank(cut) > sepadata->maxrankintegral) )
      return SCIP_OKAY;

   /* discard cut if the rank is too high */
   if( !madeintegral && (sepadata->maxrank != -1) && (SCIProwGetRank(cut) > sepadata->maxrank) )
      return SCIP_OKAY;

   *useful = TRUE;

   return SCIP_OKAY;
}


/** add cut */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR**            vars,               /**< array of variables */
   int                   c,                  /**< index of basic variable (< 0 for slack variables) */
   SCIP_Longint          maxdnom,            /**< maximal denominator to use for scaling */
   SCIP_Real             maxscale,           /**< maximal scaling factor */
   int                   cutnnz,             /**< number of nonzeros in cut */
   int*                  cutinds,            /**< variable indices in cut */
   SCIP_Real*            cutcoefs,           /**< cut cofficients */
   SCIP_Real             cutefficacy,        /**< cut efficacy */
   SCIP_Real             cutrhs,             /**< rhs of cut */
   SCIP_Bool             cutislocal,         /**< whether cut is local */
   int                   cutrank,            /**< rank of cut */
   SCIP_Bool             strongcg,           /**< whether the cut arises from the strong-CG procedure */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff appeared */
   int*                  naddedcuts          /**< pointer to store number of added cuts */
   )
{
   assert(scip != NULL);
   assert(cutoff != NULL);
   assert(naddedcuts != NULL);

   if( cutnnz == 0 && SCIPisFeasNegative(scip, cutrhs) ) /*lint !e644*/
   {
      SCIPdebugMsg(scip, " -> gomory cut detected infeasibility with cut 0 <= %g.\n", cutrhs);
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* Only take efficient cuts, except for cuts with one non-zero coefficient (= bound
    * changes); the latter cuts will be handled internally in sepastore. */
   if( SCIPisEfficacious(scip, cutefficacy) || ( cutnnz == 1 && SCIPisFeasPositive(scip, cutefficacy) ) )
   {
      SCIP_ROW* cut;
      SCIP_SEPA* cutsepa;
      char cutname[SCIP_MAXSTRLEN];
      int v;

      /* construct cut name */
      if( strongcg )
      {
         cutsepa = sepadata->strongcg;

         if( c >= 0 )
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "scg%" SCIP_LONGINT_FORMAT "_x%d", SCIPgetNLPs(scip), c);
         else
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "scg%" SCIP_LONGINT_FORMAT "_s%d", SCIPgetNLPs(scip), -c-1);
      }
      else
      {
         cutsepa = sepadata->gomory;

         if( c >= 0 )
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "gom%" SCIP_LONGINT_FORMAT "_x%d", SCIPgetNLPs(scip), c);
         else
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "gom%" SCIP_LONGINT_FORMAT "_s%d", SCIPgetNLPs(scip), -c-1);
      }

      /* create empty cut */
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, cutsepa, cutname, -SCIPinfinity(scip), cutrhs,
            cutislocal, FALSE, sepadata->dynamiccuts) );

      /* set cut rank */
      SCIProwChgRank(cut, cutrank); /*lint !e644*/

      /* cache the row extension and only flush them if the cut gets added */
      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

      /* collect all non-zero coefficients */
      for( v = 0; v < cutnnz; ++v )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[v]], cutcoefs[v]) );
      }

      /* flush all changes before adding the cut */
      SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

      if( SCIProwGetNNonz(cut) == 0 )
      {
         assert( SCIPisFeasNegative(scip, cutrhs) );
         SCIPdebugMsg(scip, " -> gomory cut detected infeasibility with cut 0 <= %g.\n", cutrhs);
         *cutoff = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIProwGetNNonz(cut) == 1 )
      {
         /* Add the bound change as cut to avoid that the LP gets modified. This would mean that the LP is not flushed
          * and the method SCIPgetLPBInvRow() fails; SCIP internally will apply this bound change automatically. */
         SCIP_CALL( SCIPaddRow(scip, cut, TRUE, cutoff) );
         ++(*naddedcuts);
      }
      else
      {
         SCIP_Bool useful;

         assert(SCIPisInfinity(scip, -SCIProwGetLhs(cut)));
         assert(!SCIPisInfinity(scip, SCIProwGetRhs(cut)));

         SCIPdebugMsg(scip, " -> %s cut <%s>: rhs=%f, eff=%f\n", strongcg ? "strong-CG" : "gomory", cutname, cutrhs, cutefficacy);

         SCIP_CALL( evaluateCutNumerics(scip, sepadata, cut, maxdnom, maxscale, &useful) );

         if( useful )
         {
            SCIPdebugMsg(scip, " -> found %s cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
               strongcg ? "strong-CG" : "gomory", cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut),
               SCIProwGetNorm(cut), SCIPgetCutEfficacy(scip, NULL, cut),
               SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
               SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));

            if( SCIPisCutNew(scip, cut) )
            {
               /* add global cuts which are not implicit bound changes to the cut pool */
               if( !cutislocal )
               {
                  if( sepadata->delayedcuts )
                  {
                     SCIP_CALL( SCIPaddDelayedPoolCut(scip, cut) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                  }
               }
               else
               {
                  /* local cuts we add to the sepastore */
                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
               }

               ++(*naddedcuts);
            }
         }
      }
      /* release the row */
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyGomory)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of separator */
   SCIP_CALL( SCIPincludeSepaGomory(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
/**! [SnippetSepaFreeGomory] */
static
SCIP_DECL_SEPAFREE(sepaFreeGomory)
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
/**! [SnippetSepaFreeGomory] */

/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitGomory)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* create and initialize random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &sepadata->randnumgen, DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;
}

/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitGomory)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeRandom(scip, &sepadata->randnumgen);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpGomory)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_VAR** vars;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_AGGRROW* aggrrow;
   SCIP_Real* binvrow;
   SCIP_Real* cutcoefs;
   SCIP_Real* basisfrac;
   int* basisind;
   int* basisperm;
   int* inds;
   int* cutinds;
   SCIP_Real maxscale;
   SCIP_Real minfrac;
   SCIP_Real maxfrac;
   SCIP_Longint maxdnom;
   SCIP_Bool cutoff;
   SCIP_Bool separatescg;
   SCIP_Bool separategmi;
   int naddedcuts;
   int nvars;
   int ncols;
   int nrows;
   int ncalls;
   int maxdepth;
   int maxsepacuts;
   int freq;
   int c;
   int i;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   minfrac = sepadata->away;
   maxfrac = 1.0 - sepadata->away;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call the gomory cut separator a given number of times at each node */
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

   /* check whether strong CG cuts should be separated */
   freq = SCIPsepaGetFreq(sepadata->strongcg);
   if( freq > 0 )
      separatescg = (depth % freq == 0);
   else
      separatescg = (freq == depth);

   /* check whether Gomory MI cuts should be separated */
   freq = SCIPsepaGetFreq(sepadata->gomory);
   if( freq > 0 )
      separategmi = (depth % freq == 0);
   else
      separategmi = (freq == depth);

   if( !separatescg && !separategmi )
      return SCIP_OKAY;

   /* get variables data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get LP data */
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   if( ncols == 0 || nrows == 0 )
      return SCIP_OKAY;

   /* set the maximal denominator in rational representation of gomory cut and the maximal scale factor to
    * scale resulting cut to integral values to avoid numerical instabilities
    */
   /**@todo find better but still stable gomory cut settings: look at dcmulti, gesa3, khb0525, misc06, p2756 */
   maxdepth = SCIPgetMaxDepth(scip);
   if( depth == 0 )
   {
      maxdnom = 1000;
      maxscale = 1000.0;
   }
   else if( depth <= maxdepth/4 )
   {
      maxdnom = 1000;
      maxscale = 1000.0;
   }
   else if( depth <= maxdepth/2 )
   {
      maxdnom = 100;
      maxscale = 100.0;
   }
   else
   {
      maxdnom = 10;
      maxscale = 10.0;
   }

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisperm, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisfrac, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nrows) );
   SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );

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
      else if( sepadata->separaterows )
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

      if( frac >= minfrac )
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

   SCIPdebugMsg(scip, "searching gomory cuts: %d cols, %d rows, maxdnom=%" SCIP_LONGINT_FORMAT ", maxscale=%g, maxcuts=%d\n",
      ncols, nrows, maxdnom, maxscale, maxsepacuts);

   cutoff = FALSE;
   naddedcuts = 0;

   /* for all basic columns belonging to integer variables, try to generate a gomory cut */
   for( i = 0; i < nrows && naddedcuts < maxsepacuts && !SCIPisStopped(scip) && !cutoff; ++i )
   {
      SCIP_Real cutrhs;
      SCIP_Real cutefficacy = 0.0;
      SCIP_Bool success;
      SCIP_Bool cutislocal;
      SCIP_Bool strongcgsuccess = FALSE;
      int ninds = -1;
      int cutnnz;
      int cutrank;
      int j;

      if( basisfrac[i] == 0.0 )
         break;

      j = basisperm[i];
      c = basisind[j];

      /* get the row of B^-1 for this basic integer variable with fractional solution value */
      SCIP_CALL( SCIPgetLPBInvRow(scip, j, binvrow, inds, &ninds) );

      SCIP_CALL( SCIPaggrRowSumRows(scip, aggrrow, binvrow, inds, ninds,
         sepadata->sidetypebasis, allowlocal, 2, (int) MAXAGGRLEN(nvars), &success) );

      if( !success )
         continue;

      /* try to create a strong CG cut out of the aggregation row */
      if( separatescg )
      {
         SCIP_CALL( SCIPcalcStrongCG(scip, NULL, POSTPROCESS, BOUNDSWITCH, USEVBDS, allowlocal, minfrac, maxfrac,
            1.0, aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, &cutrank, &cutislocal, &strongcgsuccess) );

         /* if we want to generate both cuts, add cut and reset cutefficacy and strongcgsuccess */
         if( strongcgsuccess && sepadata->genbothgomscg )
         {
            assert(allowlocal || !cutislocal); /*lint !e644*/
            SCIP_CALL( addCut(scip, sepadata, vars, c, maxdnom, maxscale, cutnnz, cutinds, cutcoefs, cutefficacy, cutrhs,
                  cutislocal, cutrank, TRUE, &cutoff, &naddedcuts) );
            cutefficacy = 0.0;
            strongcgsuccess = FALSE;
            if( cutoff )
               break;
         }
      }

      /* @todo Currently we are using the SCIPcalcMIR() function to compute the coefficients of the Gomory
       *       cut. Alternatively, we could use the direct version (see thesis of Achterberg formula (8.4)) which
       *       leads to cut a of the form \sum a_i x_i \geq 1. Rumor has it that these cuts are better.
       */

      /* try to create Gomory cut out of the aggregation row */
      if( separategmi )
      {
         /* SCIPcalcMIR will only override the cut if its efficacy is larger than the one of the strongcg cut */
         SCIP_CALL( SCIPcalcMIR(scip, NULL, POSTPROCESS, BOUNDSWITCH, USEVBDS, allowlocal, FIXINTEGRALRHS, NULL, NULL,
            minfrac, maxfrac, 1.0, aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, &cutrank, &cutislocal, &success) );

         if( success || strongcgsuccess )
         {
            assert(allowlocal || !cutislocal); /*lint !e644*/
            if( success )
               strongcgsuccess = FALSE;   /* Set strongcgsuccess to FALSE, since the MIR cut has overriden the strongcg cut. */

            SCIP_CALL( addCut(scip, sepadata, vars, c, maxdnom, maxscale, cutnnz, cutinds, cutcoefs, cutefficacy, cutrhs,
                  cutislocal, cutrank, strongcgsuccess, &cutoff, &naddedcuts) );
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basisfrac);
   SCIPfreeBufferArray(scip, &basisperm);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &cutinds);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPaggrRowFree(scip, &aggrrow);

   SCIPdebugMsg(scip, "end searching gomory cuts: found %d cuts\n", naddedcuts);

   sepadata->lastncutsfound = SCIPgetNCutsFound(scip);

   /* evaluate the result of the separation */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if ( naddedcuts > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** LP solution separation method of dummy separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpDummy)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of dummy separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolDummy)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** creates the Gomory MIR cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaGomory(
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
         sepaExeclpGomory, NULL,
         sepadata) );
   assert(sepa != NULL);

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepadata->strongcg, "strongcg", "separator for strong CG cuts", -100000, SEPA_FREQ, 0.0,
      SEPA_USESSUBSCIP, FALSE, sepaExeclpDummy, sepaExecsolDummy, NULL) );
   assert(sepadata->strongcg != NULL);

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepadata->gomory, "gomorymi", "separator for Gomory mixed-integer cuts", -100000, SEPA_FREQ, 0.0,
      SEPA_USESSUBSCIP, FALSE, sepaExeclpDummy, sepaExecsolDummy, NULL) );
   assert(sepadata->gomory != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyGomory) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeGomory) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitGomory) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitGomory) );

   /* mark main separator as a parent */
   SCIPsetSepaIsParentsepa(scip, sepa);

   /* set pointer from child separators to main separator */
   SCIPsetSepaParentsepa(scip, sepadata->strongcg, sepa);
   SCIPsetSepaParentsepa(scip, sepadata->gomory, sepa);

   /* add separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gomory/maxrounds",
         "maximal number of gomory separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gomory/maxroundsroot",
         "maximal number of gomory separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gomory/maxsepacuts",
         "maximal number of gomory cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gomory/maxsepacutsroot",
         "maximal number of gomory cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gomory/maxrank",
         "maximal rank of a gomory cut that could not be scaled to integral coefficients (-1: unlimited)",
         &sepadata->maxrank, FALSE, DEFAULT_MAXRANK, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/gomory/maxrankintegral",
         "maximal rank of a gomory cut that could be scaled to integral coefficients (-1: unlimited)",
         &sepadata->maxrankintegral, FALSE, DEFAULT_MAXRANKINTEGRAL, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/gomory/away",
         "minimal integrality violation of a basis variable in order to try Gomory cut",
         &sepadata->away, FALSE, DEFAULT_AWAY, 1e-4, 0.5, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gomory/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gomory/makeintegral",
         "try to scale cuts to integral coefficients",
         &sepadata->makeintegral, TRUE, DEFAULT_MAKEINTEGRAL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gomory/forcecuts",
         "if conversion to integral coefficients failed still consider the cut",
         &sepadata->forcecuts, TRUE, DEFAULT_FORCECUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gomory/separaterows",
         "separate rows with integral slack",
         &sepadata->separaterows, TRUE, DEFAULT_SEPARATEROWS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gomory/delayedcuts",
         "should cuts be added to the delayed cut pool?",
         &sepadata->delayedcuts, TRUE, DEFAULT_DELAYEDCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gomory/sidetypebasis",
         "choose side types of row (lhs/rhs) based on basis information?",
         &sepadata->sidetypebasis, TRUE, DEFAULT_SIDETYPEBASIS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gomory/trystrongcg",
         "try to generate strengthened Chvatal-Gomory cuts?",
         &sepadata->trystrongcg, TRUE, DEFAULT_TRYSTRONGCG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/gomory/genbothgomscg",
         "Should both Gomory and strong CG cuts be generated (otherwise take best)?",
         &sepadata->genbothgomscg, TRUE, DEFAULT_GENBOTHGOMSCG, NULL, NULL) );

   return SCIP_OKAY;
}
