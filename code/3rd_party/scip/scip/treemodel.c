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

/**@file   treemodel.c
 * @brief  Branching rules based on the Single-Variable-Branching (SVB) model
 * @author Daniel Anderson
 * @author Pierre Le Bodic
 *
 * The Single-Variable-Branching (SVB) model is a simplified model of
 * Branch & Bound trees, from which several nontrivial variable selection
 * rules arise. The Treemodel branching rule complements SCIP's hybrid
 * branching by suggesting improved branching variables given the current
 * pseudocosts and the current dual gap.
 *
 * Given a variable with dual bound changes (l, r) (both positive)
 * and an absolute gap G, the SVB model describes the tree that needs to be
 * built by branching on that same variable at every node until the value G
 * is reached at every leaf, starting from 0 at the root node.
 * If we do so for every variable, we can select the variable that produces
 * the smallest tree.
 * In the case where the gap is not known, then we can compute the growth rate
 * of the tree, which we call the ratio.
 * The ratio of a variable (l, r) is the factor by which the size of the tree
 * built using (l, r) that closes a gap G must be multiplied by to close a gap
 * G+1. This ratio is not constant for all gaps, but when G tends to infinity,
 * it converges to a fixed value we can compute numerically using a root finding
 * algorithm (e.g. Laguerre).
 * The ratio is used when the gap is too large (e.g. no primal bound known) or
 * to help approximate the size of the SVB tree for that variable.
 *
 * See the following publication for more detail:
 *
 * @par
 * Pierre Le Bodic and George Nemhauser@n
 * An abstract model for branching and its application to mixed integer programming@n
 * Mathematical Programming, 2017@n
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/treemodel.h"

#include "scip/history.h"
#include "scip/var.h"

#include <limits.h>

#define LAGUERRE_THRESHOLD    100      /**< Maximum value of r/l at which Laguerre is the prefered FP method */

/* Default parameters for the Treemodel branching rules */
#define DEFAULT_ENABLE         FALSE    /**< should candidate branching variables be scored using the Treemodel rule? */
#define DEFAULT_HIGHRULE       'r'      /**< scoring function to use at nodes predicted to be high in the tree.
					  * ('d'efault, 's'vts, 'r'atio, 't'ree sample) */
#define DEFAULT_LOWRULE        'r'      /**< scoring function to use at nodes predicted to be low in the tree
					  * ('d'efault, 's'vts, 'r'atio, 't'ree sample) */
#define DEFAULT_HEIGHT         10       /**< estimated tree height at which we switch from using the low rule to
					  * the high rule */
#define DEFAULT_FILTERHIGH     'a'      /**< should dominated candidates be filtered before using the high scoring
					  * function? ('a'uto, 't'rue, 'f'alse) */
#define DEFAULT_FILTERLOW      'a'      /**< should dominated candidates be filtered before using the low scoring
					  * function? ('a'uto, 't'rue, 'f'alse) */
#define DEFAULT_MAXFPITER      24       /**< maximum number of fixed-point iterations when computing the ratio */
#define DEFAULT_MAXSVTSHEIGHT  100      /**< maximum height to compute the SVTS score exactly before approximating */
#define DEFAULT_FALLBACKINF    'r'      /**< which method should be used as a fallback if the tree size estimates are
					  * infinite? ('d'efault, 'r'atio) */
#define DEFAULT_FALLBACKNOPRIM 'r'      /**< which method should be used as a fallback if there is no primal bound
					  * available? ('d'efault, 'r'atio) */
#define DEFAULT_SMALLPSCOST    0.1      /**< threshold at which pseudocosts are considered small, making hybrid scores
					  * more likely to be the deciding factor in branching */

/** parameters required by the Treemodel branching rules */
struct SCIP_Treemodel
{
   SCIP_Bool            enabled;             /**< should candidate branching variables be scored using the Treemodel
					       * rule? */
   char                 highrule;            /**< scoring function to use at nodes predicted to be high in the tree.
					       * ('d'efault, 's'vts, 'r'atio, 't'ree sample) */
   char                 lowrule;             /**< scoring function to use at nodes predicted to be low in the tree
					       * ('d'efault, 's'vts, 'r'atio, 't'ree sample) */
   int                  height;              /**< estimated tree height at which we switch from using the low rule to
					       * the high rule */
   char                 filterhigh;          /**< should dominated candidates be filtered before using the high
					       * scoring function? ('a'uto, 't'rue, 'f'alse) [ADVANCED] */
   char                 filterlow;           /**< should dominated candidates be filtered before using the low
					       * scoring function? ('a'uto, 't'rue, 'f'alse) [ADVANCED] */
   int                  maxfpiter;           /**< maximum number of fixed-point iterations when computing the ratio
					       * [ADVANCED] */
   int                  maxsvtsheight;       /**< maximum height to compute the SVTS score exactly before approximating
					       * [ADVANCED] */
   char                 fallbackinf;         /**< which method should be used as a fallback if the tree size estimates
					       * are infinite? ('d'efault, 'r'atio) [ADVANCED] */
   char                 fallbacknoprim;      /**< which method should be used as a fallback if there is no primal bound
					       * available? ('d'efault, 'r'atio) [ADVANCED] */
   SCIP_Real            smallpscost;         /**< threshold at which pseudocosts are considered small, making hybrid
					       * scores more likely to be the deciding factor in branching [ADVANCED] */
};

/** branching encoding of a variable's ratio
 * A variable's ratio is defined based upon its left and right LP gains, as the unique root > 1 of
 * the polynomial x^r - x^(r-l) -1, where l and r are the left and right LP gains.
 * We store the root as upratio^(invleft), with invleft = 1/l. The value upratio is thus
 * the ratio of the variable (1, r/l).
 * An extra boolean stores whether the encoded ratio is valid,
 * i.e. there were no numerical problems when computing it */
struct SCIP_Ratio
{
   SCIP_Real             upratio;           /**< "UnPowered" ratio, i.e. the ratio of the characteristic polynomial
					      * with gains (1, rightgain/leftgain) */
   SCIP_Real             invleft;           /**< "INVerse left gain, i.e. 1/leftgain */
   SCIP_Bool             valid;             /**< True iff the ratio computed is valid */
};
typedef struct SCIP_Ratio SCIP_RATIO;

/** a comparison method for the next method. It simply compares two SCIP_Real */
static
SCIP_DECL_SORTINDCOMP(sciprealcomp)
{
   SCIP_Real* value = (SCIP_Real*) dataptr;
   SCIP_Real diffval;

   assert(value != NULL);
   assert(ind1 >= 0 && ind2 >= 0);

   diffval = value[ind1] - value[ind2];
   if( diffval < 0.0 )
      return -1;
   else if( diffval > 0.0)
      return 1;
   else
      return 0;
}

/** given a pair of arrays of real non-negative values (a,b), with a <= b, computes
 * the pairs that belong to the pareto front (with a tolerance).
 * In other words, we are looking for non-dominated pairs of values.
 * One value and one array are computed after this method.
 * The value is the number of non-dominated elements.
 * The array is a boolean array that indicates if an element is dominated.
 * In case of a draw, only one variable is considered as non-dominated.
 */
static
SCIP_RETCODE findNonDominatedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            a,                  /**< the first set of values */
   SCIP_Real*            b,                  /**< the second set of values */
   int                   size,               /**< the size of array a (and b) */
   int*                  ndominated,         /**< returns the number of dominated elements */
   SCIP_Bool*            dominated           /**< returns the array of booleans that determine if an element is
					       * dominated */
   )
{
   SCIP_Real bestcurrenta;
   SCIP_Real besta;
   SCIP_Real currentb;
   int* permb;
   int* bestcurrents;
   int nbestcurrent;
   int indexinpermb;
   int origindex;
   int iterbestcurrent;

   assert(scip != NULL);
   assert(a != NULL);
   assert(b != NULL);
   assert(ndominated != NULL);
   assert(dominated != NULL);
   assert(size > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &bestcurrents, size) );

   /* we first find the permutation of indices of array b that corresponds to
    * the array of a non-increasing sort of its values */
   SCIP_CALL( SCIPallocBufferArray(scip, &permb, size) );
   for( origindex=0; origindex<size; ++origindex )
      permb[origindex] = origindex;

   SCIPsortDownInd(permb, sciprealcomp, (void*)b, size);

   *ndominated = 0;
   /* Now we will traverse the pair of arrays a and b by non-decreasing order of values of b
    * and mark the (non) dominated pairs */

   /* The current max value of a for all pairs that (almost) have the same value b */
   bestcurrenta = a[permb[0]];

   /* the current value b */
   currentb = b[permb[0]];
   /* the best pair(s) for the current value b */
   bestcurrents[0] = permb[0];
   nbestcurrent = 1;
   /* the value a to "beat" to be non-dominated */
   besta = -1;
   for( indexinpermb = 1; indexinpermb < size; ++indexinpermb )
   {
      origindex = permb[indexinpermb];
      assert(b[origindex] <= currentb);
      if( SCIPisLT(scip, b[origindex], currentb) )
      {
         /* If the above is true, then we went through all the previous elements that had value currentb */
         /* Thus the best element for value currentb is non-dominated if its value bestcurrenta is better
          * than the previous best besta */
         if( bestcurrenta > besta )
         {
            for( iterbestcurrent=0; iterbestcurrent < nbestcurrent; ++iterbestcurrent )
               dominated[bestcurrents[iterbestcurrent]] = FALSE;

            besta = bestcurrenta;
         }
         else
         {
            for( iterbestcurrent = 0; iterbestcurrent < nbestcurrent; ++iterbestcurrent )
            {
               dominated[bestcurrents[iterbestcurrent]] = TRUE;
               ++(*ndominated);
            }
         }
         bestcurrenta = a[origindex];
         currentb = b[origindex];
         bestcurrents[0] = origindex;
         nbestcurrent = 1;
      }
      else
      {
         /* Then the b values are (almost) equal and we need to compare values a */
         if( SCIPisGT(scip, a[origindex], bestcurrenta) )
         {
            /* Then the new value is better than the old one(s) */
            for( iterbestcurrent = 0; iterbestcurrent < nbestcurrent; ++iterbestcurrent )
            {
               dominated[bestcurrents[iterbestcurrent]] = TRUE;
               ++(*ndominated);
            }

            bestcurrenta = a[origindex];
            bestcurrents[0] = origindex;
            nbestcurrent = 1;
         }
         else
         {
            /* Then the new value is equal or dominated */
            if( SCIPisEQ(scip, a[origindex], bestcurrenta) )
            {
               bestcurrents[nbestcurrent] = origindex;
               ++nbestcurrent;
            }
            else
            {
               dominated[origindex] = TRUE;
               ++(*ndominated);
            }
         }
      }
   }
   /* Finally, we have to look at the last best variable */
   if( bestcurrenta > besta )
   {
      for( iterbestcurrent = 0; iterbestcurrent < nbestcurrent; ++iterbestcurrent )
         dominated[bestcurrents[iterbestcurrent]] = FALSE;
   }
   else
   {
      for( iterbestcurrent = 0; iterbestcurrent < nbestcurrent; ++iterbestcurrent )
      {
         dominated[bestcurrents[iterbestcurrent]] = TRUE;
         ++(*ndominated);
      }
   }

   SCIPfreeBufferArray(scip, &permb);
   SCIPfreeBufferArray(scip, &bestcurrents);
   return SCIP_OKAY;
}

/** returns true iff the variable with given gains has a ratio better (i.e smaller) than the given one */
static
SCIP_Bool hasBetterRatio(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RATIO*           branchratio,        /**< The variable's ratio to compare against */
   SCIP_Real             leftgain,           /**< the left gain of a variable */
   SCIP_Real             rightgain           /**< the right gain of a variable */
   )
{
   SCIP_Real result;

   assert(branchratio != NULL);
   assert(branchratio->valid);
   assert(SCIPisLE(scip, leftgain, rightgain));

   /* We evaluate the characteristic polynomial of the variable on the given ratio. */
   result = -1;
   if( leftgain > 0.0 && rightgain > 0.0 )
   {
      result = pow(branchratio->upratio, rightgain * branchratio->invleft) - pow(branchratio->upratio, (rightgain - leftgain) * branchratio->invleft) - 1; /*lint !e644*/
   }

   /* If the result is positive, then it has a better ratio. */
   return (result > 0.0);
}

/** computes the variable ratio corresponding to the left and right gains */
static
void computeVarRatio(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR*             var,                /**< the candidate branching variable */
   SCIP_Real             leftgain,           /**< the left gain of the variable */
   SCIP_Real             rightgain,          /**< the right gain of the variable */
   SCIP_RATIO*           branchratio         /**< storage for the computed ratio */
   )
{
   SCIP_Real ratio;
   SCIP_Real newratio;
   SCIP_Real r;
   int iters;

   assert(SCIPisGE(scip, leftgain, 0.0));
   assert(SCIPisGE(scip, rightgain, leftgain));

   if( SCIPisZero(scip, leftgain) || SCIPisZero(scip, rightgain) )
   {
      branchratio->valid = FALSE;
      return;
   }

   /* We scale left and right gains by dividing both by left */
   r = rightgain / leftgain;

   /* In the case where l and r are very close r may become < 1 */
   if( r <= 1 )
   {
      branchratio->valid = TRUE;
      branchratio->upratio = 2.0;
      branchratio->invleft = 1.0 / leftgain;
      return;
   }

   /* Check if this ratio has already been computed */
   if( SCIPhistoryIsRatioValid(var->history) && SCIPisEQ(scip, SCIPhistoryGetLastBalance(var->history), r) )
   {
      branchratio->valid = TRUE;
      branchratio->upratio = SCIPhistoryGetLastRatio(var->history);
      branchratio->invleft = 1.0 / leftgain;
      return;
   }

   /* Initialise the ratio at the previously computed ratio (if applicable) otherwise
    * use the lower bound 2^(1/r) <= phi <= 2^(1/l).
    * Note that we only use the previous ratio if the previous value of r/l was larger,
    * ie. the previous ratio was smaller since we want to initialise at a lower bound.
    */
   ratio = 1.0;
   newratio = pow(2.0, 1.0/r);
   if( SCIPhistoryIsRatioValid(var->history) && SCIPhistoryGetLastBalance(var->history) > r
      && SCIPhistoryGetLastRatio(var->history) > newratio )
      newratio = SCIPhistoryGetLastRatio(var->history);

   /* Depending on the value of rightgain/leftgain, we have two different methods to compute the ratio
    * If this value is bigger than 100, we use a fixed-point method. Otherwise, we use Laguerre's method
    * This is strictly for numerical efficiency and based on experiments.
    */

   /* Use Laguerre's method */
   if( r <= LAGUERRE_THRESHOLD )
   {
      /* We relax the epsilon after 5 iterations since we may not have enough precision to achieve any better
       * convergence */
      for( iters = 0; ((iters <= 5 && !SCIPisEQ(scip, ratio, newratio)) ||
              (iters > 5 && !SCIPisSumEQ(scip, ratio, newratio)))
              && iters < treemodel->maxfpiter && newratio > 1.0; iters++ )
      {
         double G, H, a, p, p1, p2, phi_r;

         ratio = newratio;
         phi_r = pow(ratio, r);
         p = phi_r - phi_r / ratio - 1.0;
         if( p != 0 )
         {
            p1 = (r * phi_r - (r - 1.0) * phi_r / ratio) / ratio;
            p2 = (r * (r - 1.0) *  phi_r - (r - 1.0) * (r - 2.0) * phi_r / ratio) / ratio / ratio;
            G = p1 / p;
            H = G * G - (p2 / p);
            a = r / (G + (G >= 0 ? 1.0 : -1.0) * sqrt((r - 1.0) * (r * H - G * G)));
            newratio = ratio - a;
         }
      }
   }
   /* Use fixed point method */
   else
   {
      /* We relax the epsilon after 10 iterations since we may not have enough precision to achieve any better
       * convergence */
      for( iters = 0; ((iters <= 10 && !SCIPisEQ(scip, ratio, newratio)) ||
              (iters > 10 && !SCIPisSumEQ(scip, ratio, newratio)))
              && iters < treemodel->maxfpiter && newratio > 1; iters++ )
      {
         ratio = newratio;
         newratio = pow(1.0-1.0/ratio, -1.0/r);
      }
   }

   /* We think that everything worked.
    * Note that the fixed point method is not guaranteed to converge due to numerical precision issues.
    * In the case that the method fails to converge, a fallback strategy must be used.
    * For instance, if this method is used for branching, then this variable can be ignored,
    * or the scores of all variables could be recomputed using a different method. */
   if( iters < treemodel->maxfpiter && newratio > 1.0 )
   {
      branchratio->valid = TRUE;
      branchratio->upratio = (ratio + newratio) / 2.0;
      branchratio->invleft = 1.0 / leftgain;
   }
   /* We (hopefully) make finding bugs easier by setting these values */
   else
   {
      branchratio->valid = FALSE;
      branchratio->upratio = -1.0;
      branchratio->invleft = -1.0;
   }

   /* Update the history */
   SCIPhistorySetRatioHistory(var->history, branchratio->valid, branchratio->upratio, r);
}

/** use the Ratio scoring function to select a branching candidate */
static
SCIP_RETCODE selectCandidateUsingRatio(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR**            branchcands,        /**< branching candidate storage */
   SCIP_Real*            mingains,           /**< minimum gain of rounding downwards or upwards */
   SCIP_Real*            maxgains,           /**< maximum gain of rounding downwards or upwards */
   SCIP_Bool             filterdominated,    /**< whether dominated variables have been filtered */
   SCIP_Bool*            dominated,          /**< whether each variable is dominated or not */
   int                   nbranchcands,       /**< the number of branching candidates */
   int*                  bestcand            /**< the best branching candidate found before the call,
					          and the best candidate after the call (possibly the same) */
   )
{
   SCIP_RATIO branchratio;
   SCIP_RATIO bestbranchratio;
   int c;

   /* We initialize bestbranchratio at the default bestcand ratio, since it is likely to have
    * a very good ratio and save evaluations of the ratio for many variables */
   int referencevar = *bestcand;
   computeVarRatio(scip, treemodel, branchcands[referencevar], mingains[referencevar], maxgains[referencevar], &bestbranchratio);

   for( c = 0; c < nbranchcands; ++c )
   {
      if( (!filterdominated || !dominated[c]) && c != referencevar )
      {
         if( !bestbranchratio.valid || hasBetterRatio(scip, &bestbranchratio, mingains[c], maxgains[c]) ) /*lint !e644*/
         {
            computeVarRatio(scip, treemodel, branchcands[c], mingains[c], maxgains[c], &branchratio);
            if( branchratio.valid ) /*lint !e644*/
            {
               *bestcand = c;
               bestbranchratio = branchratio;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** Returns the predicted treesize for the gap and given up and down gains */
static
SCIP_Real computeSVTS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR*             var,                /**< the candidate branching variable */
   SCIP_Real             absgap,             /**< the absolute gap to close (typically the local gap at the current node) */
   SCIP_Real             mingain,            /**< prediction of smaller objective gain of downwards/upwards */
   SCIP_Real             maxgain             /**< prediction of larger objective gain of downwards/upwards */
   )
{
   SCIP_Real prediction = SCIP_REAL_MAX;

   if( SCIPisGT(scip, mingain, 0.0) && !SCIPisInfinity(scip, absgap) )
   {
      SCIP_Real treesize;
      SCIP_Real gaptoreach;
      SCIP_Real scaledgap;
      SCIP_Real scaledgain;
      int mindepth;
      int nr;
      int ir;

      /* We implicitly set the minimum gain to 1, and the maximum gain and gap accordingly,
       * as the treesize does not change if we scale the gains and gap by a scalar  */
      scaledgain = maxgain / mingain;
      scaledgap = absgap / mingain;

      mindepth = (int) SCIPceil(scip, scaledgap / scaledgain);

      /* In the following case we compute the treesize for a smaller gap
       * and we will deduce the treesize of the scaledgap using the ratio */
      if( mindepth > treemodel->maxsvtsheight )
      {
         gaptoreach = scaledgap * (treemodel->maxsvtsheight - 1) / mindepth;

         assert(!SCIPisInfinity(scip, gaptoreach));
         assert(gaptoreach > scaledgain);
      }
      else
      {
         gaptoreach = scaledgap;
      }

      mindepth = (int) ceil(gaptoreach / scaledgain);
      assert(mindepth <= treemodel->maxsvtsheight);
      treesize = 1;

      /* nr is the number of times we turn right to reach a leaf */
      for( nr = 1; nr <= mindepth; ++nr )
      {
         SCIP_Real binomcoeff = 1.0;
         for( ir = 1; ir <= nr; ++ir )
         {
            binomcoeff *= (nr + ceil((gaptoreach - (nr - 1) * scaledgain)) - ir) / ir;
         }
         treesize += binomcoeff;
      }

      treesize = 2.0 * treesize - 1.0;

      assert(SCIPisGE(scip, treesize, 3.0));

      if( !SCIPisEQ(scip, scaledgap, gaptoreach) )
      {
         /* If we have not computed the treesize for the scaled gap but for max gap instead,
          * we use the ratio between two iterations to come up with an estimate of the treesize
          * for the scaled gap */
         if( !SCIPisInfinity(scip,treesize) )
         {
            SCIP_RATIO branchratio;
            computeVarRatio(scip, treemodel, var, mingain, maxgain, &branchratio);

            if( branchratio.valid ) /*lint !e644*/
               prediction = treesize * pow(branchratio.upratio, (scaledgap - gaptoreach) * branchratio.invleft); /*lint !e644*/
         }
      }
      else
      {
         prediction = treesize;
      }
   }

   return prediction;
}

/** use the SVTS scoring function to select a branching candidate */
static
SCIP_RETCODE selectCandidateUsingSVTS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR**            branchcands,        /**< branching candidate storage */
   SCIP_Real*            mingains,           /**< minimum gain of rounding downwards or upwards */
   SCIP_Real*            maxgains,           /**< maximum gain of rounding downwards or upwards */
   SCIP_Real*            tiebreakerscore,    /**< scores to use for tie breaking */
   SCIP_Real             localabsgap,        /**< The dual gap at the current node */
   SCIP_Bool             filterdominated,    /**< whether dominated variables have been filtered */
   SCIP_Bool*            dominated,          /**< whether each variable is dominated or not */
   int                   nbranchcands,       /**< the number of branching candidates */
   int                   ndominated,         /**< the number of dominated candidates */
   int*                  bestcand            /**< the best branching candidate found before the call,
					          and the best candidate after the call (possibly the same) */
   )
{
   SCIP_Real* treesizes;
   SCIP_Real referencetreesize;
   SCIP_Real score;
   SCIP_Real bestscore;
   SCIP_Real avgtreesize;
   int besttscand;
   int referencevar;
   int c;

   /* We will first measure the treesize for scip's default variable. If it is infinite then we don't compute
    * the treesize for other variables (even though it might be finite) and go directly to the fallback strategy */
   besttscand = *bestcand;
   referencevar = *bestcand;

   treesizes = NULL;
   bestscore = 0.0;
   avgtreesize = 0.0;
   if( !SCIPisInfinity(scip, localabsgap) )
   {
      referencetreesize = computeSVTS(scip, treemodel, branchcands[referencevar], localabsgap, mingains[referencevar],
		             maxgains[referencevar]);
      if( !SCIPisInfinity(scip, referencetreesize) )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &treesizes, nbranchcands) );
         treesizes[referencevar] = referencetreesize;

         for( c = 0; c < nbranchcands; ++c )
         {
            if( !filterdominated || !dominated[c] )
            {
               if( c != referencevar )
                  treesizes[c] = computeSVTS(scip, treemodel, branchcands[c], localabsgap, mingains[c], maxgains[c]);
               else
                  treesizes[c] = referencetreesize;

               avgtreesize += treesizes[c];
            }
            else
               treesizes[c] = SCIP_REAL_MAX;
         }
         avgtreesize = avgtreesize / (nbranchcands - ndominated);

         for( c = 0; c < nbranchcands; ++c )
         {
            score = (1.0 - 1.0 / (1.0 + avgtreesize / treesizes[c])) + 0.01 * tiebreakerscore[c];
            if(score > bestscore)
            {
               bestscore = score;
               besttscand = c;
            }
         }

         *bestcand = besttscand;

         SCIPfreeBufferArray(scip, &treesizes);
      }
      /* Apply infinite treesize fallback strategy */
      else if( treemodel->fallbackinf == 'r' )
      {
         SCIP_CALL( selectCandidateUsingRatio(scip, treemodel, branchcands, mingains, maxgains, filterdominated, dominated,
               nbranchcands, bestcand) );
      }
   }
   /* Apply no primal bound fallback strategy */
   else if( treemodel->fallbacknoprim == 'r' )
   {
      SCIP_CALL( selectCandidateUsingRatio(scip, treemodel, branchcands, mingains, maxgains, filterdominated, dominated,
            nbranchcands, bestcand) );
   }

   return SCIP_OKAY;
}

/** computes a^b for integer b */
static
SCIP_Real integerpow(
   SCIP_Real             a,                  /**< the base */
   int                   b                   /**< the integer exponent */
   )
{  /*lint --e{644}*/
   SCIP_Real ans;

   ans = 1.0;
   for( ; b; b /= 2 )
   {
      if( b & 1 )
         ans *= a;
      a *= a;
   }
   return ans;
}

/** returns the sampled tree size for the given lp gains and dual gap */
static
SCIP_Real computeSampleTreesize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR*             var,                /**< the candidate branching variable */
   SCIP_Real             absgap,             /**< the absolute gap to close (typically the local at the current node) */
   SCIP_Real             leftgain,           /**< The minimum gain from branching on this variable */
   SCIP_Real             rightgain           /**< The maximum gain from branching on this variable */
   )
{
   SCIP_RATIO branchratio;
   SCIP_Real prediction;
   SCIP_Real leftsize;
   SCIP_Real rightsize;
   SCIP_Real midsize;

   computeVarRatio(scip, treemodel, var, leftgain, rightgain, &branchratio);

   if( branchratio.valid ) /*lint !e644*/
   {  /*lint --e{644}*/
      SCIP_Real phi_l = branchratio.upratio;
      SCIP_Real phi_r = pow(branchratio.upratio, rightgain * branchratio.invleft);
      int kl = (int)ceil(absgap / leftgain);
      int kr = (int)ceil(absgap / rightgain);
      int k = (int)ceil(absgap / (leftgain + rightgain));
      SCIP_Real phi_lr = phi_l * phi_r;
      SCIP_Real phi_klr = integerpow(phi_lr, k);

      /* We compute an estimate of the size of the tree using the left-most leaf,
       * right-most leaf, and the leaf obtained from alternating left and right. */
      leftsize = (integerpow(phi_l, kl + 1) - 1.0) / (phi_l - 1.0);
      rightsize = (integerpow(phi_r, kr + 1) - 1.0) / (phi_r - 1.0);

      if( k * (leftgain + rightgain) < absgap + rightgain )
         midsize = (1.0 + phi_l) * (phi_klr * phi_lr - 1.0) / (phi_lr - 1.0) - phi_klr * phi_l;
      else
         midsize = (1.0 + phi_l) * (phi_klr - 1.0) / (phi_lr - 1.0);

      prediction = (leftsize + rightsize + midsize) / 3.0;
   }
   else
   {
      prediction = SCIP_REAL_MAX;
   }

   return prediction;
}

/** use the Tree Sampling scoring function to select a branching candidate */
static
SCIP_RETCODE selectCandidateUsingSampling(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR**            branchcands,        /**< branching candidate storage */
   SCIP_Real*            mingains,           /**< minimum gain of rounding downwards or upwards */
   SCIP_Real*            maxgains,           /**< maximum gain of rounding downwards or upwards */
   SCIP_Real*            tiebreakerscore,    /**< scores to use for tie breaking */
   SCIP_Real             localabsgap,        /**< The dual gap at the current node */
   SCIP_Bool             filterdominated,    /**< whether dominated variables have been filtered */
   SCIP_Bool*            dominated,          /**< whether each variable is dominated or not */
   int                   nbranchcands,       /**< the number of branching candidates */
   int                   ndominated,         /**< the number of dominated candidates */
   int*                  bestcand            /**< the best branching candidate found before the call,
					          and the best candidate after the call (possibly the same) */
   )
{
   SCIP_Real* treesizes;
   SCIP_Real referencetreesize;
   SCIP_Real score;
   SCIP_Real bestscore;
   SCIP_Real avgtreesize;
   int besttscand;
   int referencevar;
   int c;

   /* We will first measure the treesize for scip's default variable. If it is infinite then we don't compute
    * the treesize for other variables (even though it might be finite) and go directly to the fallback strategy */
   besttscand = *bestcand;
   referencevar = *bestcand;

   treesizes = NULL;
   bestscore = 0.0;
   avgtreesize = 0.0;
   if( !SCIPisInfinity(scip, localabsgap) )
   {
      referencetreesize = computeSampleTreesize(scip, treemodel, branchcands[referencevar], localabsgap, mingains[referencevar],
		             maxgains[referencevar]);

      if( !SCIPisInfinity(scip, referencetreesize) )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &treesizes, nbranchcands) );
         treesizes[referencevar] = referencetreesize;

         for( c = 0; c < nbranchcands; ++c )
         {
            if( !filterdominated || !dominated[c] )
            {
               if( c != referencevar )
                  treesizes[c] = computeSampleTreesize(scip, treemodel, branchcands[c], localabsgap, mingains[c], maxgains[c]);
               else
                  treesizes[c] = referencetreesize;

               avgtreesize += treesizes[c];
            }
            else
               treesizes[c] = SCIP_REAL_MAX;
         }
         avgtreesize = avgtreesize / (nbranchcands - ndominated);

         for( c = 0; c < nbranchcands; ++c )
         {
            score = (1.0 - 1.0 / (1.0 + avgtreesize / treesizes[c])) + 0.01 * tiebreakerscore[c];
            if( score > bestscore )
            {
               bestscore = score;
               besttscand = c;
            }
         }

         *bestcand = besttscand;

         SCIPfreeBufferArray(scip, &treesizes);
      }
      /* Apply infinite treesize fallback strategy */
      else if( treemodel->fallbackinf == 'r' )
      {
         SCIP_CALL( selectCandidateUsingRatio(scip, treemodel, branchcands, mingains, maxgains, filterdominated, dominated,
               nbranchcands, bestcand) );
      }
   }
   /* Apply no primal bound fallback strategy */
   else if( treemodel->fallbacknoprim == 'r' )
   {
      SCIP_CALL( selectCandidateUsingRatio(scip, treemodel, branchcands, mingains, maxgains, filterdominated, dominated,
            nbranchcands, bestcand) );
   }

   return SCIP_OKAY;
}

/** initialises the Treemodel parameter data structure */
SCIP_RETCODE SCIPtreemodelInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL**      treemodel           /**< Treemodel parameter data structure */
   )
{
   assert(treemodel != NULL);
   SCIP_CALL( SCIPallocBlockMemory(scip, treemodel) );
   assert(*treemodel != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "branching/treemodel/enable",
         "should candidate branching variables be scored using the Treemodel branching rules?",
         &(*treemodel)->enabled, FALSE, DEFAULT_ENABLE,
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/highrule",
         "scoring function to use at nodes predicted to be high in the tree ('d'efault, 's'vts, 'r'atio, 't'ree sample)",
         &(*treemodel)->highrule, FALSE, DEFAULT_HIGHRULE, "dsrt",
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/lowrule",
         "scoring function to use at nodes predicted to be low in the tree ('d'efault, 's'vts, 'r'atio, 't'ree sample)",
         &(*treemodel)->lowrule, FALSE, DEFAULT_LOWRULE, "dsrt",
         NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/treemodel/height",
         "estimated tree height at which we switch from using the low rule to the high rule",
         &(*treemodel)->height, FALSE, DEFAULT_HEIGHT, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/filterhigh",
         "should dominated candidates be filtered before using the high scoring function? ('a'uto, 't'rue, 'f'alse)",
         &(*treemodel)->filterhigh, TRUE, DEFAULT_FILTERHIGH, "atf",
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/filterlow",
         "should dominated candidates be filtered before using the low scoring function? ('a'uto, 't'rue, 'f'alse)",
         &(*treemodel)->filterlow, TRUE, DEFAULT_FILTERLOW, "atf",
         NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/treemodel/maxfpiter",
         "maximum number of fixed-point iterations when computing the ratio",
         &(*treemodel)->maxfpiter, TRUE, DEFAULT_MAXFPITER, 1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "branching/treemodel/maxsvtsheight",
         "maximum height to compute the SVTS score exactly before approximating",
         &(*treemodel)->maxsvtsheight, TRUE, DEFAULT_MAXSVTSHEIGHT, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/fallbackinf",
         "which method should be used as a fallback if the tree size estimates are infinite? ('d'efault, 'r'atio)",
         &(*treemodel)->fallbackinf, TRUE, DEFAULT_FALLBACKINF, "dr",
         NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "branching/treemodel/fallbacknoprim",
         "which method should be used as a fallback if there is no primal bound available? ('d'efault, 'r'atio)",
         &(*treemodel)->fallbacknoprim, TRUE, DEFAULT_FALLBACKNOPRIM, "dr",
         NULL, NULL) );
   SCIP_CALL ( SCIPaddRealParam(scip, "branching/treemodel/smallpscost",
         "threshold at which pseudocosts are considered small, making hybrid scores more likely to be the deciding factor in branching",
         &(*treemodel)->smallpscost, TRUE, DEFAULT_SMALLPSCOST,
         0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** frees the Treemodel parameter data structure */
SCIP_RETCODE SCIPtreemodelFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL**      treemodel           /**< Treemodel parameter data structure */
   )
{
   assert(treemodel != NULL);
   assert(*treemodel != NULL);

   SCIPfreeBlockMemory(scip, treemodel);

   assert(*treemodel == NULL);

   return SCIP_OKAY;
}

/** returns TRUE if the Treemodel branching rules are enabled */
SCIP_Bool SCIPtreemodelIsEnabled(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel           /**< Treemodel parameter data structure */
   )
{
   assert(scip != NULL);
   return treemodel->enabled;
}

/** apply the Treemodel branching rules to attempt to select a better
 *  branching candidate than the one selected by pseudocost branching
 */
SCIP_RETCODE SCIPtreemodelSelectCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TREEMODEL*       treemodel,          /**< Treemodel parameter data structure */
   SCIP_VAR**            branchcands,        /**< branching candidate storage */
   SCIP_Real*            mingains,           /**< minimum gain of rounding downwards or upwards */
   SCIP_Real*            maxgains,           /**< maximum gain of rounding downwards or upwards */
   SCIP_Real*            tiebreakerscore,    /**< scores to use for tie breaking */
   int                   nbranchcands,       /**< the number of branching candidates */
   int*                  bestcand            /**< the best branching candidate found before the call,
					          and the best candidate after the call (possibly the same) */
   )
{
   SCIP_Real localabsgap;           /* The gap at the current node */
   int bestcandheight;              /* The height of the best candidate according to SCIP */
   char scoringfunction;            /* Scoring function to use (based on the estimated tree height) */
   char filtersetting;              /* Whether we should apply filtering of dominated variables */

   assert(treemodel != NULL);
   assert(treemodel->enabled);
   assert(*bestcand >= 0);

   /* Compute the dual gap at the current node */
   if( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) )
      localabsgap = SCIPgetUpperbound(scip) - SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip));
   else
      localabsgap = SCIPinfinity(scip);

   /* Compute an estimate of the height of the current node using the bestcand variable */
   if( !SCIPisInfinity(scip, localabsgap) && SCIPisGT(scip, mingains[*bestcand], 0.0)
      && SCIPisLT(scip, localabsgap/mingains[*bestcand], 1.0 * INT_MAX))
      bestcandheight = (int)(localabsgap/mingains[*bestcand]);
   else
      bestcandheight = INT_MAX;

   /* Decide which scoring function to use based on the estimated height of the tree */
   if( bestcandheight < treemodel->height )
   {
      scoringfunction = treemodel->lowrule;
      filtersetting = treemodel->filterlow;
   }
   else
   {
      scoringfunction = treemodel->highrule;
      filtersetting = treemodel->filterhigh;
   }

   /* We are going to apply a Treemodel variable selection rule */
   if( scoringfunction != 'd' )
   {
      SCIP_Bool* dominated;            /* Whether variables are dominated */
      SCIP_Bool autofilter;            /* If auto filtering is chosen, should variables be filtered? */
      SCIP_Bool filterdominated;       /* Whether we should filter dominated variables */
      int ndominated;                  /* Number of dominated variables */

      /* Filtering dominated variables is suggested for SVTS and Tree Sampling rules */
      autofilter = (filtersetting == 'a' && (scoringfunction == 's' || scoringfunction == 't'));
      filterdominated = (autofilter || filtersetting == 't');

      /* If selected, find the dominated variables */
      if( filterdominated )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &dominated, nbranchcands) );
         SCIP_CALL( findNonDominatedVars(scip, mingains, maxgains, nbranchcands, &ndominated, dominated) );
      }
      else
      {
         dominated = NULL;
         ndominated = 0;
      }

      /* Invoke the selected scoring function */
      switch( scoringfunction )
      {
      case 's':
         SCIP_CALL( selectCandidateUsingSVTS(scip, treemodel, branchcands, mingains, maxgains, tiebreakerscore,
               localabsgap, filterdominated, dominated, nbranchcands, ndominated, bestcand) );
         break;
      case 'r':
         SCIP_CALL( selectCandidateUsingRatio(scip, treemodel, branchcands, mingains, maxgains, filterdominated,
               dominated, nbranchcands, bestcand) );
         break;
      case 't':
         SCIP_CALL( selectCandidateUsingSampling(scip, treemodel, branchcands, mingains, maxgains, tiebreakerscore,
               localabsgap, filterdominated, dominated, nbranchcands, ndominated, bestcand) );
         break;
      default:
         return SCIP_PARAMETERWRONGVAL;
      }

      /* Free dominated variable buffer if it was used */
      if( filterdominated )
      {
         assert(dominated != NULL);
         SCIPfreeBufferArray(scip, &dominated);
      }
   }

   return SCIP_OKAY;
}
