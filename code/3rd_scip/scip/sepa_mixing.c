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

/**@file   sepa_mixing.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  mixing/star inequality separator
 * @author Weikun Chen
 * @author Marc Pfetsch
 *
 * This separator generates cuts based on the mixing set
 * \f[
 * X = \{ (x,y) \in \{0,1\}^{N \cup M} \times \mathbb{R} \, : \,  y \geq a_i x_i, \, \textrm{for} \, i \in N, \,
 *                           y \leq u -  a_i x_i, \, \textrm{for} \, i \in M, \, 0 \leq y \leq u  \},
 * \f]
 * where \f$0 \leq a_i \leq u\f$ for all \f$i\f$. This information can be obtained directly from the variable bounds data
 * structure. The separator will generate three classes of cuts.
 *
 * VLB: Let \f$T\f$ be a subset of \f$N\f$, wlog, \f$T = \{1,\ldots,r\}\f$ with \f$a_1 \leq a_2 \leq \ldots \leq a_r\f$.
 * Let \f$a_0 = 0\f$. The mixing/star VLB cut is of the form \f$ y \geq \sum_{i=1}^r (a_i - a_{i-1})x_i \f$.
 *
 * VUB: Let \f$T\f$ be a subset of \f$M\f$, wlog, \f$T = \{1,\ldots,r\}\f$ with \f$a_1 \leq a_2 \leq \ldots \leq a_r\f$.
 * Let \f$a_0 = 0\f$. The mixing/star VUB cut is of the form \f$ y \leq u - \sum_{i=1}^r (a_i - a_{i-1})x_i \f$.
 *
 * CONFLICT: Consider \f$i \in N\f$ and \f$j \in M\f$ with \f$a_i + a_j > u\f$. The conflict cut is
 * \f$x_i + x_j \leq 1\f$.
 *
 * A small example is described in the following to see the generated cuts.
 * \f[
 * Y = \{ (x,y) \in \{0,1\}^{4} \times \mathbb{R} \, : \, y \geq 2x_1, \, y \geq 3x_2, \, y \leq 4 - x_3, \,
 *                           y \leq 4 - 2 x_4, \, 0 \leq y \leq 4 \}.
 * \f]
 * In this small example, the mixing/star cuts \f$y \geq 2x_1 + x_2\f$ (VLB) and \f$y \leq 4 - x_3 - x_4\f$ (VUB) will be
 * considered to be generated. Besides the mixing cuts, we also consider the conflict cut \f$x_1 + x_3 \leq 1\f$ (CONFLICT).
 *
 *
 * For an overview see:
 * Atamturk, A., Nemhauser, G.L. and Savelsbergh, M.W.,@n
 * The mixed vertex packing problem.@n
 * Mathematical Programming, 89(1), 35-53, 2000.
 *
 * Some remarks:
 * - Besides the mixing inequality, we also add the conflict inequality.
 * - Currently, the performance is bad on the neos-565672 instance.
 *   The reason is that, after adding the separator, SCIP spends a lot of time at the stage of cutting plane generation.
 * - We do not consider sparsity of the cuts as we aim to find a most violated cut.
 * - Besides the most violated cut we consider, we also add an additional variable to make the cut be the strongest one,
 *   even the additional variable does not contribute any to the violation.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sepa.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cut.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sepa.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"
#include "scip/sepa_mixing.h"
#include "scip/scip_tree.h"
#include "scip/sepa_mixing.h"
#include <string.h>


#define SEPA_NAME              "mixing"
#define SEPA_DESC              "mixing inequality separator"
#define DEFAULT_MAXROUNDS            -1 /**< maximal number of mixing separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of mixing separation rounds in the root node (-1: unlimited) */
#define SEPA_PRIORITY               -50
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_USELOCALBOUNDS    FALSE /**< should local bounds be used? */
#define DEFAULT_ISCUTSONINTS      FALSE /**< should general/implied integer variables be used to generate cuts? */
#define DEFAULT_MAXNUNSUCCESSFUL     10 /**< maximal number of consecutive unsuccessful iterations */

/** separator-specific data for the mixing separator */
struct SCIP_SepaData
{
   SCIP_Bool             uselocalbounds;     /**< should local bounds be used? */
   SCIP_Bool             iscutsonints;       /**< should general/implied integer variables be used to generate cuts? */
   int                   maxrounds;          /**< maximal number of mixing separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of mixing separation rounds in the root node (-1: unlimited) */
   int                   nunsuccessful;      /**< number of consecutive unsuccessful iterations */
   int                   maxnunsuccessful;   /**< maximal number of consecutive unsuccessful iterations */
};

/*
 * local methods
 */

/** adds the given cut */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            cutcoefs,           /**< coefficients of active variables in cut */
   int*                  cutinds,            /**< problem indices of variables in cut */
   int                   cutnnz,             /**< number of non-zeros in cut */
   SCIP_Real             cutrhs,             /**< right hand side of cut */
   SCIP_Bool             cutislocal,         /**< Is the cut only locally valid? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to update number of cuts added */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_VAR** vars;
   SCIP_ROW* cut;
   int v;

   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
   assert(ncuts != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* get active problem variables */
   vars = SCIPgetVars(scip);

   /* construct cut name */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "mix%" SCIP_LONGINT_FORMAT "_x%d", SCIPgetNLPs(scip), *ncuts);

   /* create empty cut */
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs,
         cutislocal, FALSE, TRUE) );

   /* cache the row extension and only flush them if the cut gets added */
   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

   /* collect all non-zero coefficients */
   for( v = 0; v < cutnnz; ++v )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[v]], cutcoefs[v]) );
   }

   /* flush all changes before adding the cut */
   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

   if( SCIPisCutEfficacious(scip, sol, cut) )
   {
      /* set cut rank */
      SCIProwChgRank(cut, 1);

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "-> found cut (eff: %f): ", SCIPgetCutEfficacy(scip, sol, cut));
      SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

      if( cutislocal )
      {
         /* local cuts are added to the sepastore */
         SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
      }
      else
      {
         SCIP_CALL( SCIPaddPoolCut(scip, cut) );
      }
      (*ncuts)++;
   }

   /* release the row */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** searches and adds mixing cuts that are violated by the given solution value array
 *
 *  This function implements a separation heuristic that runs in linear time in comparison to the quadratic exact
 *  algorithm in Atamturk et al.:
 *  - Lower and upper bounds are considered separately. Then possible conflict cuts.
 *  - Variable lower/upper bounds data are collected, i.e., the corresponding binary variables and coefficients.
 *  - These lists are sorted non-increasingly according to the solution values.
 *  - Then binary variables are added in turn as long as their coefficients increase in order to make the coefficients
 *    nonnegative.  This clearly makes the cuts heuristic, since it is order dependent, but also sparser.
 *  - The procedure stops as soon as it reaches 0 solution values.
 *  - If the cut is efficious it is added.
 *
 *  @todo Check whether one wants to sort according to the quotient of solution values and coefficients to increase the
 *  chance of having smaller coefficients in the front.
 */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to store the number of generated cuts */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* vlbmixcoefs;
   SCIP_Real* vlbmixsols;
   SCIP_Real* vubmixcoefs;
   SCIP_Real* vubmixsols;
   SCIP_Real* cutcoefs;
   SCIP_Real cutrhs;
   int* vlbmixinds;
   int* vubmixinds;
   int* cutinds;
   int* vlbmixsigns;
   int* vubmixsigns;
   int firstvar;
   int nmaxvars;
   int nvars;
   int i;
   int k;

   assert(sepa != NULL);
   assert(cutoff != NULL);
   assert(ncuts != NULL);

   *cutoff = FALSE;
   *ncuts = 0;

   /* exit if there are no binary variables - ignore integer variables that have 0/1 bounds */
   nmaxvars = SCIPgetNBinVars(scip);
   if( nmaxvars <= 0 )
      return SCIP_OKAY;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* get the index of the first considered variable */
   if( sepadata->iscutsonints )
   {
      /* generate cuts based on all nonbinary variabels */
      firstvar = SCIPgetNBinVars(scip);
   }
   else
   {
      /* only generate cuts based on continuous variables */
      firstvar = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);
   }
   nvars = SCIPgetNVars(scip);
   if( firstvar == nvars )
      return SCIP_OKAY;

   vars = SCIPgetVars(scip);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &vlbmixcoefs, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vlbmixsols, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vlbmixinds, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vlbmixsigns, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vubmixcoefs, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vubmixsols, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vubmixinds, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vubmixsigns, nmaxvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nmaxvars + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, nmaxvars + 1) );

   for( i = firstvar; i < nvars; i++ )
   {
      SCIP_VAR** vlbvars;
      SCIP_Real* vlbcoefs;
      SCIP_Real* vlbconsts;
      SCIP_VAR** vubvars;
      SCIP_Real* vubcoefs;
      SCIP_Real* vubconsts;
      SCIP_Real maxabscoef;
      SCIP_Real activity;
      SCIP_Real lastcoef;
      SCIP_Real varsolval;
      SCIP_Real lb = SCIP_INVALID;
      SCIP_Real ub = SCIP_INVALID;
      SCIP_Bool islocallb = FALSE;  /* Is it a local lower bound or global lower bound? */
      SCIP_Bool islocalub = FALSE;  /* Is it a local upper bound or global upper bound? */
      SCIP_Bool cutislocal; /* Is it a local cut or global cut? */
      int vlbmixsize = 0;
      int vubmixsize = 0;
      int cutnnz = 0;
      int maxabsind;
      int maxabssign;
      int nvlb;
      int nvub;
      int j;

      var = vars[i];
      assert( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY );

      if( SCIPvarGetProbindex(var) < 0 )
         continue;

      nvlb = SCIPvarGetNVlbs(var);
      nvub = SCIPvarGetNVubs(var);

      if( nvlb == 0 && nvub == 0 )
         continue;

      /* skip lower bound if the LP solution value is equal to the upper bound of the continuous variable */
      varsolval = SCIPgetSolVal(scip, sol, var);
      if( SCIPisFeasEQ(scip, SCIPvarGetUbLocal(var), varsolval) )
         goto VUB;

      if( nvlb == 0 )
         goto VUB;

      /* get variable lower variable bounds information */
      vlbvars = SCIPvarGetVlbVars(var);
      vlbcoefs = SCIPvarGetVlbCoefs(var);
      vlbconsts = SCIPvarGetVlbConstants(var);

      maxabscoef = 0.0;
      maxabsind = -1;
      maxabssign = 0;

      lb = SCIPvarGetLbGlobal(var);
      if( sepadata->uselocalbounds && SCIPisLT(scip, lb, SCIPvarGetLbLocal(var)) )
      {
         /* this is a lcoal cut */
         islocallb = TRUE;
         lb = SCIPvarGetLbLocal(var);
      }

#ifdef SCIP_DEBUG
      for( j = 0; j < nvlb; j++ )
      {
         SCIP_Real tmplb;

         if( SCIPvarIsBinary(vlbvars[j]) && SCIPvarGetProbindex(vlbvars[j]) >= 0 )
         {
            tmplb = (vlbcoefs[j] > 0) ? vlbconsts[j] : (vlbconsts[j] + vlbcoefs[j]);
            assert( SCIPisFeasLE(scip, tmplb, lb) );
         }
      }
#endif

      assert( SCIPisFeasLE(scip, lb, SCIPvarGetUbLocal(var)) );

      /* extract the useful variable bounds information (binary and nonredundant) */
      for( j = 0; j < nvlb; j++ )
      {
         /* consider only active and binary variables */
         if( SCIPvarIsBinary(vlbvars[j]) && SCIPvarGetProbindex(vlbvars[j]) >= 0 )
         {
            SCIP_Real maxactivity;
            SCIP_Real coef;

            maxactivity = (vlbcoefs[j] > 0) ? (vlbconsts[j] + vlbcoefs[j]) : vlbconsts[j];

            /* skip redundant variable bound constraints */
            if( SCIPisFeasLE(scip, maxactivity, lb) )
               continue;

            if( vlbcoefs[j] > 0 )
            {
               coef = maxactivity - lb;
               vlbmixsigns[vlbmixsize] = 0;
            }
            else
            {
               coef = lb - maxactivity;
               vlbmixsigns[vlbmixsize] = 1;
            }
            assert(vlbmixsize <= nvars);

            vlbmixcoefs[vlbmixsize] = REALABS(coef);
            vlbmixinds[vlbmixsize] = SCIPvarGetProbindex(vlbvars[j]);
            vlbmixsols[vlbmixsize] = (! vlbmixsigns[vlbmixsize]) ? SCIPgetSolVal(scip, sol, vlbvars[j]) : (1.0 - SCIPgetSolVal(scip, sol, vlbvars[j]));

            /* update the maximal coefficient if needed */
            if( maxabscoef < vlbmixcoefs[vlbmixsize] )
            {
               maxabscoef = vlbmixcoefs[vlbmixsize];
               maxabsind = vlbmixinds[vlbmixsize];
               maxabssign = vlbmixsigns[vlbmixsize];
            }

            ++vlbmixsize;

            /* stop if size is exceeded; possibly ignore redundant variable bounds */
            if( vlbmixsize >= nmaxvars )
               break;
         }
      }
      assert( vlbmixsize <= nmaxvars );

      /* stop if no variable lower bounds information exists */
      if( vlbmixsize == 0 )
         goto VUB;

      /* stop if the current solution value of the transformed continuous variable is larger than the maximal coefficient */
      if( SCIPisFeasGT(scip, varsolval - lb, maxabscoef) )
         goto VUB;

      /* sort the solution values in non-increasing order */
      SCIPsortDownRealRealIntInt(vlbmixsols, vlbmixcoefs, vlbmixinds, vlbmixsigns, vlbmixsize);

      /* add the continuous variable */
      cutcoefs[cutnnz] = -1;
      cutinds[cutnnz] = SCIPvarGetProbindex(var);
      cutrhs = -lb;
      cutnnz++;

      activity = -(varsolval - lb);
      lastcoef = 0.0;

      /* loop over the variables and add the variable to the cut if its coefficient is larger than that of the last variable */
      for( j = 0; j < vlbmixsize; j++ )
      {
         SCIP_Real solval;

         solval = vlbmixsols[j];

         /* stop if we can not find a violated cut or we reached 0 solution values */
         if( activity + solval * (maxabscoef - lastcoef) < 0.0 || SCIPisFeasZero(scip, solval) )
            break;
         else
         {
            /* skip if we have already added a variable with bigger coefficient */
            if( SCIPisLE(scip, vlbmixcoefs[j], lastcoef) )
               continue;
            else
            {
               activity += (vlbmixcoefs[j] - lastcoef) * solval;
               if( vlbmixsigns[j] )
               {
                  cutcoefs[cutnnz] = lastcoef - vlbmixcoefs[j];
                  cutrhs -= vlbmixcoefs[j] - lastcoef;
               }
               else
                  cutcoefs[cutnnz] = vlbmixcoefs[j] - lastcoef;
               cutinds[cutnnz++] = vlbmixinds[j];
               lastcoef = vlbmixcoefs[j];
            }
         }
      }

      /* add the variable with maximal coefficient to make sure the cut is strong enough */
      if( SCIPisGT(scip, maxabscoef, lastcoef) )
      {
         /* do not update activity: either the corresponding solval is 0.0 or the cut will not become efficious */
         if( maxabssign )
         {
            cutcoefs[cutnnz] = lastcoef - maxabscoef;
            cutrhs -= maxabscoef - lastcoef;
         }
         else
            cutcoefs[cutnnz] = maxabscoef - lastcoef;
         cutinds[cutnnz++] = maxabsind;
      }
      assert( cutnnz <= nmaxvars + 1 );

      /* add the cut if the violtion is good enough and the number of nonzero coefficients is larger than 2 */
      if( SCIPisEfficacious(scip, activity) && cutnnz > 2 )
      {
         SCIP_CALL( addCut(scip, sepa, sol, cutcoefs, cutinds, cutnnz, cutrhs, islocallb, cutoff, ncuts) );
      }

   VUB:
      if( nvub == 0 )
         goto CONFLICT;

      /* get variable upper bounds information */
      vubvars = SCIPvarGetVubVars(var);
      vubcoefs = SCIPvarGetVubCoefs(var);
      vubconsts = SCIPvarGetVubConstants(var);

      maxabscoef = 0.0;
      maxabsind = -1;
      maxabssign = 0;

      /* stop if the lower bound is equal to the solution value of the continuous variable */
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), varsolval) )
         goto CONFLICT;

      ub = SCIPvarGetUbGlobal(var);
      if( sepadata->uselocalbounds && SCIPisGT(scip, ub, SCIPvarGetUbLocal(var)) )
      {
         /* this is a lcoal cut */
         islocalub = TRUE;
         ub = SCIPvarGetUbLocal(var);
      }

#ifdef SCIP_DEBUG
      for( j = 0; j < nvub; j++ )
      {
         SCIP_Real tmpub;

         if( SCIPvarIsBinary(vubvars[j]) && SCIPvarGetProbindex(vubvars[j]) >= 0 )
         {
            tmpub = (vubcoefs[j] < 0) ? vubconsts[j] : (vubconsts[j] + vubcoefs[j]);
            assert( SCIPisFeasGE(scip, tmpub, ub) );
         }
      }
#endif

      assert( SCIPisFeasLE(scip, SCIPvarGetLbLocal(var), ub) );

      /* extract the useful variable bounds information (binary and nonredundant) */
      for( j = 0; j < nvub; j++ )
      {
         /* consider only active and binary variables */
         if( SCIPvarIsBinary(vubvars[j]) && SCIPvarGetProbindex(vubvars[j]) >= 0 )
         {
            SCIP_Real minactivity;
            SCIP_Real coef;

            minactivity = (vubcoefs[j] < 0) ? (vubconsts[j] + vubcoefs[j]) : vubconsts[j];

            /* skip redundant variable bound constraints */
            if( SCIPisFeasLE(scip, ub, minactivity) )
               continue;

            if( vubcoefs[j] > 0 )
            {
               coef = ub - minactivity;
               vubmixsigns[vubmixsize] = 1;
            }
            else
            {
               coef = minactivity - ub;
               vubmixsigns[vubmixsize] = 0;
            }

            vubmixcoefs[vubmixsize] = REALABS(coef);
            vubmixinds[vubmixsize] = SCIPvarGetProbindex(vubvars[j]);
            vubmixsols[vubmixsize] = (! vubmixsigns[vubmixsize]) ? SCIPgetSolVal(scip, sol, vubvars[j]): (1.0 - SCIPgetSolVal(scip, sol, vubvars[j]));

            /* update the maximal coefficient if needed */
            if( maxabscoef < vubmixcoefs[vubmixsize] )
            {
               maxabscoef = vubmixcoefs[vubmixsize];
               maxabsind = vubmixinds[vubmixsize];
               maxabssign = vubmixsigns[vubmixsize];
            }

            ++vubmixsize;

            /* stop if size is exceeded; possibly ignore redundant variable bounds */
            if( vubmixsize >= nmaxvars )
               break;
         }
      }
      assert( vubmixsize <= nmaxvars );

      /* stop if no variable upper bounds information exists */
      if( vubmixsize == 0 )
         goto CONFLICT;

      /* stop if the current solution value of transformed continuous variable is larger than the maximal coefficient */
      if( SCIPisFeasGT(scip, ub - varsolval, maxabscoef) )
         goto CONFLICT;

      /* sort the solution values in non-increasing order */
      SCIPsortDownRealRealIntInt(vubmixsols, vubmixcoefs, vubmixinds, vubmixsigns, vubmixsize);

      /* add the continuous variables */
      cutnnz = 0;
      cutcoefs[cutnnz] = 1;
      cutinds[cutnnz] = SCIPvarGetProbindex(var);
      cutrhs = ub;
      cutnnz++;

      activity = varsolval - ub;
      lastcoef = 0.0;

      for( j = 0; j < vubmixsize; j++ )
      {
         SCIP_Real solval;

         solval = vubmixsols[j];

         /* stop if we can not find a violated cut or we reached 0 solution values */
         if( activity + solval * (maxabscoef - lastcoef) < 0.0 || SCIPisFeasZero(scip, solval) )
            break;
         else
         {
            /* skip if we have already added a variable with bigger coefficient */
            if( SCIPisLE(scip, vubmixcoefs[j], lastcoef) )
               continue;
            else
            {
               activity += (vubmixcoefs[j] - lastcoef) * solval;
               if( vubmixsigns[j] )
               {
                  cutcoefs[cutnnz] = lastcoef - vubmixcoefs[j];
                  cutrhs -= vubmixcoefs[j] - lastcoef;
               }
               else
                  cutcoefs[cutnnz] = vubmixcoefs[j] - lastcoef;
               cutinds[cutnnz++] = vubmixinds[j];
               lastcoef = vubmixcoefs[j];
            }
         }
      }

      /* add the variable with maximal coefficient if needed */
      if( SCIPisGT(scip, maxabscoef, lastcoef) )
      {
         /* do not update activity: either the corresponding solval is 0.0 or the cut will not become efficious */
         if( maxabssign )
         {
            cutcoefs[cutnnz] = lastcoef - maxabscoef;
            cutrhs -= maxabscoef - lastcoef;
         }
         else
            cutcoefs[cutnnz] = maxabscoef - lastcoef;
         cutinds[cutnnz++] = maxabsind;
      }
      assert( cutnnz <= nmaxvars + 1 );

      /* add the cut if the violtion is good enough and the number of nonzero coefficients is larger than 2 */
      if( SCIPisEfficacious(scip, activity) && cutnnz > 2 )
      {
         SCIP_CALL( addCut(scip, sepa, sol, cutcoefs, cutinds, cutnnz, cutrhs, islocalub, cutoff, ncuts) );
      }

   CONFLICT:
      /* combine the variable lower bounds information and upper bounds information together to generate cuts */
      /* stop if no useful variable lower (or upper) bounds information exists */
      if( vlbmixsize == 0 || vubmixsize == 0 )
         continue;

      assert( lb != SCIP_INVALID ); /*lint !e777*/
      assert( ub != SCIP_INVALID ); /*lint !e777*/

      cutislocal = islocallb || islocalub;
      for( j = 0; j < vlbmixsize; j++ )
      {
         SCIP_Real solval;

         solval = vlbmixsols[j];

         /* stop if no violated cut exists */
         if( ! SCIPisEfficacious(scip, solval + vubmixsols[0] - 1.0) )
            break;

         for( k = 0; k < vubmixsize; k++ )
         {
            /* only consider the inequality if its violation is good enough */
            if( SCIPisEfficacious(scip, solval + vubmixsols[k] - 1.0) )
            {
               SCIP_Real tmp;

               tmp = lb + vlbmixcoefs[j] + vubmixcoefs[k] - ub;

               /* add the cut if it is valid */
               if( SCIPisEfficacious(scip, tmp) )
               {
                  cutnnz = 2;
                  cutrhs = 1.0;
                  cutcoefs[0] = vlbmixsigns[j] ? -1.0 : 1.0;
                  cutcoefs[1] = vubmixsigns[k] ? -1.0 : 1.0;
                  cutinds[0] = vlbmixinds[j];
                  cutinds[1] = vubmixinds[k];
                  cutrhs = vlbmixsigns[j] ? (cutrhs - 1.0) : cutrhs;
                  cutrhs = vubmixsigns[k] ? (cutrhs - 1.0) : cutrhs;
                  SCIP_CALL( addCut(scip, sepa, sol, cutcoefs, cutinds, cutnnz, cutrhs, cutislocal, cutoff, ncuts) );
               }
            }
            else
               break;
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cutinds);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &vubmixsigns);
   SCIPfreeBufferArray(scip, &vubmixinds);
   SCIPfreeBufferArray(scip, &vubmixsols);
   SCIPfreeBufferArray(scip, &vubmixcoefs);
   SCIPfreeBufferArray(scip, &vlbmixsigns);
   SCIPfreeBufferArray(scip, &vlbmixinds);
   SCIPfreeBufferArray(scip, &vlbmixsols);
   SCIPfreeBufferArray(scip, &vlbmixcoefs);

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyMixing)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of separator */
   SCIP_CALL( SCIPincludeSepaMixing(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMixing)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* get separation data and free it */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   SCIPfreeBlockMemory(scip, &sepadata);

   /* reset data pointer to NULL */
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMixing)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_Bool cutoff;
   int nbinvars;
   int nvars;
   int ncuts;
   int ncalls;

   assert(sepa != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* only call the mixing cut separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* gets numver of active problem variables and number of binary variables */
   SCIP_CALL( SCIPgetVarsData(scip, NULL, &nvars, &nbinvars, NULL, NULL, NULL) );

   /* if all the active problem variables are binary, stop */
   if( nvars == nbinvars )
      return SCIP_OKAY;

   /* call the cut separation */
   SCIP_CALL( separateCuts(scip, sepa, NULL, &cutoff, &ncuts) );

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( ncuts > 0 )
   {
      SCIPdebugMsg(scip, "mixing separator generated %d cuts.\n", ncuts);
      *result = SCIP_SEPARATED;
   }
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecSolMixing)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_Bool cutoff;
   int nbinvars;
   int nvars;
   int ncuts;
   int ncalls;

   assert(sepa != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* do not run if we have reached the maximal number of consecutive unsuccessful calls */
   if( sepadata->nunsuccessful >= sepadata->maxnunsuccessful )
      return SCIP_OKAY;

   /* only call the mixing cut separator a given number of times at each node */
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* gets numver of active problem variables and number of binary variables */
   nvars = SCIPgetNVars(scip);
   nbinvars = SCIPgetNBinVars(scip);

   /* if all the active problem variables are binary, stop */
   if( nvars == nbinvars )
      return SCIP_OKAY;

   /* call the cut separation */
   SCIP_CALL( separateCuts(scip, sepa, sol, &cutoff, &ncuts) );

   /* adjust result code */
   if( cutoff )
   {
      sepadata->nunsuccessful = 0;
      *result = SCIP_CUTOFF;
   }
   else if( ncuts > 0 )
   {
      SCIPdebugMsg(scip, "mixing separator generated %d cuts.\n", ncuts);
      sepadata->nunsuccessful = 0;
      *result = SCIP_SEPARATED;
   }
   else
   {
      ++sepadata->nunsuccessful;
      *result = SCIP_DIDNOTFIND;
   }

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the mixing separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMixing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create mixing separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   assert(sepadata != NULL);
   sepadata->nunsuccessful = 0;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpMixing, sepaExecSolMixing,
         sepadata) );
   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMixing) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMixing) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/mixing/uselocalbounds",
         "Should local bounds be used?",
         &sepadata->uselocalbounds, TRUE, DEFAULT_USELOCALBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/mixing/iscutsonints",
         "Should general integer variables be used to generate cuts?",
         &sepadata->iscutsonints, TRUE, DEFAULT_ISCUTSONINTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/mixing/maxrounds",
         "maximal number of mixing separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/mixing/maxroundsroot",
         "maximal number of mixing separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/mixing/maxnunsuccessful",
         "maximal number of consecutive unsuccessful iterations",
         &sepadata->maxnunsuccessful, FALSE, DEFAULT_MAXNUNSUCCESSFUL, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
