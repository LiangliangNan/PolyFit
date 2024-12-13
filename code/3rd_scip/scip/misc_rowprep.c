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

/**@file   misc_rowprep.c
 * @ingroup OTHER_CFILES
 * @brief  linear inequalities in preparation
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_misc_rowprep.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_var.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_sepa.h"
#include "scip/scip_sol.h"
#include "scip/scip_tree.h"
#include "scip/struct_misc.h"
#include "scip/struct_scip.h"
#include "scip/set.h"

#define ROWPREP_SCALEUP_VIOLNONZERO    (10.0*SCIPepsilon(scip))    /**< minimal violation for considering up-scaling of rowprep (we want to avoid upscaling very small violations) */
#define ROWPREP_SCALEUP_MINVIOLFACTOR  2.0                         /**< scale up will target a violation of ~MINVIOLFACTOR*minviol, where minviol is given by caller */
#define ROWPREP_SCALEUP_MAXMINCOEF     (1.0 / SCIPfeastol(scip))   /**< scale up only if min. coef is below this number (before scaling) */
#define ROWPREP_SCALEUP_MAXMAXCOEF     SCIPgetHugeValue(scip)      /**< scale up only if max. coef will not exceed this number by scaling */
#define ROWPREP_SCALEUP_MAXSIDE        SCIPgetHugeValue(scip)      /**< scale up only if side will not exceed this number by scaling */
#define ROWPREP_SCALEDOWN_MINMAXCOEF   (1.0 / SCIPfeastol(scip))   /**< scale down if max. coef is at least this number (before scaling) */
#define ROWPREP_SCALEDOWN_MINCOEF      SCIPfeastol(scip)           /**< scale down only if min. coef does not drop below this number by scaling */

#ifndef M_SQRT2
#define M_SQRT2 sqrt(2.0)
#endif

/** adds a variable to the `rowprep->modifiedvars` array, if recording of modification has been enabled and the variable is not fixed */
static
SCIP_RETCODE rowprepRecordModifiedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_VAR*             var                 /**< variable to add */
   )
{
   int oldsize;

   if( !rowprep->recordmodifications )
      return SCIP_OKAY;

   /* do not record for fixed variables, as they are not suitable for branching */
   if( SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      SCIPdebugMsg(scip, "skip recording modification for fixed variable <%s>[%g,%g]\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      return SCIP_OKAY;
   }

   /* increase modifiedvars array size */
   if( rowprep->nmodifiedvars >= rowprep->modifiedvarssize )
   {
      oldsize = rowprep->modifiedvarssize;
      rowprep->modifiedvarssize = SCIPcalcMemGrowSize(scip, rowprep->nmodifiedvars + 1);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &rowprep->modifiedvars, oldsize, rowprep->modifiedvarssize) );
   }

   rowprep->modifiedvars[rowprep->nmodifiedvars] = var;
   ++rowprep->nmodifiedvars;

   return SCIP_OKAY;
}

/** sort terms by absolute value of coefficients, from largest to smallest */
static
SCIP_RETCODE rowprepCleanupSortTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep             /**< rowprep to be sorted */
   )
{
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);

   /* special treatment for cuts with few variables */
   switch( rowprep->nvars )
   {
      case 0:
      case 1:
         break;

      case 2:
      {
         if( REALABS(rowprep->coefs[0]) < REALABS(rowprep->coefs[1]) )
         {
            SCIP_Real tmp1;
            SCIP_VAR* tmp2;

            tmp1 = rowprep->coefs[0];
            rowprep->coefs[0] = rowprep->coefs[1];
            rowprep->coefs[1] = tmp1;

            tmp2 = rowprep->vars[0];
            rowprep->vars[0] = rowprep->vars[1];
            rowprep->vars[1] = tmp2;
         }
         break;
      }

      default :
      {
         SCIP_Real* abscoefs;

         SCIP_CALL( SCIPallocBufferArray(scip, &abscoefs, rowprep->nvars) );
         for( i = 0; i < rowprep->nvars; ++i )
            abscoefs[i] = REALABS(rowprep->coefs[i]);
         SCIPsortDownRealRealPtr(abscoefs, rowprep->coefs, (void**)rowprep->vars, rowprep->nvars);
         SCIPfreeBufferArray(scip, &abscoefs);
      }
   }

   /* forget about coefs that are exactly zero (unlikely to have some) */
   while( rowprep->nvars > 0 && rowprep->coefs[rowprep->nvars-1] == 0.0 )
      --rowprep->nvars;

   return SCIP_OKAY;
}

/** try to improve coef range by aggregating row with variable bounds
 *
 * Assumes terms have been sorted by rowprepCleanupSortTerms().
 */
static
SCIP_RETCODE rowprepCleanupImproveCoefrange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_SOL*             sol,                /**< solution that we try to cut off, or NULL for LP solution */
   SCIP_Real             maxcoefrange        /**< maximal allowed coefficients range */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real ref;
   SCIP_Real coef;
   SCIP_Real mincoef;
   SCIP_Real maxcoef;
   SCIP_Real loss[2];
   int maxcoefidx;
   int pos;

   maxcoefidx = 0;
   if( rowprep->nvars > 0 )
   {
      maxcoef = REALABS(rowprep->coefs[0]);
      mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
   }
   else
      mincoef = maxcoef = 1.0;

   /* eliminate minimal or maximal coefs as long as coef range is too large
    * this is likely going to eliminate coefs that are within eps of 0.0
    * if not, then we do so after scaling (or should we enforce this here?)
    */
   while( maxcoef / mincoef > maxcoefrange  )
   {
      SCIPdebugMsg(scip, "cut coefficients have very large range: mincoef = %g maxcoef = %g\n", mincoef, maxcoef);

      /* max/min can only be > 1 if there is more than one var
       * we need this below for updating the max/min coef after eliminating a term
       */
      assert(rowprep->nvars > 1);

      /* try to reduce coef range by aggregating with variable bounds
       * that is, eliminate a term like a*x from a*x + ... <= side by adding -a*x <= -a*lb(x)
       * with ref(x) the reference point we try to eliminate, this would weaken the cut by a*(lb(x)-ref(x))
       *
       * we consider eliminating either the term with maximal or the one with minimal coefficient,
       * taking the one that leads to the least weakening of the cut
       *
       * TODO (suggested by @bzfserra, see !496):
       * - Also one could think of not completely removing the coefficient but do an aggregation that makes the coefficient look better. For instance:
       *   say you have $`a x + 0.1 y \leq r`$ and $`y`$ has only an upper bound, $`y \leq b`$,
       *   then you can't really remove $`y`$. However, you could aggregate it with $`0.9 \cdot (y \leq b)`$ to get
       *   $`a x + y \leq r + 0.9 b`$, which has better numerics (and hopefully still cuts the point... actually, if for the point you want to separate, $`y^* = b`$, then the violation is the same)
       */

      for( pos = 0; pos < 2; ++pos )
      {
         var = rowprep->vars[pos ? rowprep->nvars-1 : maxcoefidx];
         coef = rowprep->coefs[pos ? rowprep->nvars-1 : maxcoefidx];
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         ref = SCIPgetSolVal(scip, sol, var);
         assert(coef != 0.0);

         /* make sure reference point is something reasonable within the bounds, preferable the value from the solution */
         if( SCIPisInfinity(scip, REALABS(ref)) )
            ref = 0.0;
         ref = MAX(lb, MIN(ub, ref));

         /* check whether we can eliminate coef*var from rowprep and how much we would loose w.r.t. ref(x) */
         if( ((coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT)) )
         {
            /* we would need to aggregate with -coef*var <= -coef*lb(x) */
            if( SCIPisInfinity(scip, -lb) )
               loss[pos] = SCIP_INVALID;
            else
               loss[pos] = REALABS(coef) * (ref - lb);
         }
         else
         {
            assert((coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT));
            /* we would need to aggregate with -coef*var >= -coef*ub(x) */
            if( SCIPisInfinity(scip, ub) )
               loss[pos] = SCIP_INVALID;
            else
               loss[pos] = REALABS(coef) * (ub - ref);
         }
         assert(loss[pos] >= 0.0);  /* assuming SCIP_INVALID >= 0 */

         SCIPdebugMsg(scip, "aggregating %g*<%s> %c= ... with <%s>[%g] %c= %g looses %g\n",
            coef, SCIPvarGetName(var), rowprep->sidetype == SCIP_SIDETYPE_RIGHT ? '<' : '>',
            SCIPvarGetName(var), ref,
            ((coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT)) ? '>' : '<',
            ((coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT)) ? lb : ub, loss[pos]);
      }

      /*lint --e{777} */
      if( loss[0] == SCIP_INVALID && loss[1] == SCIP_INVALID )
         break;  /* cannot eliminate coefficient */

      /* select position with smaller loss */
      pos = (loss[1] == SCIP_INVALID || loss[1] > loss[0]) ? 0 : 1;

      /* now do the actual elimination */
      var = rowprep->vars[pos ? rowprep->nvars-1 : maxcoefidx];
      coef = rowprep->coefs[pos ? rowprep->nvars-1 : maxcoefidx];

      /* eliminate coef*var from rowprep: increase side */
      if( ((coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT)) )
      {
         /* we aggregate with -coef*var <= -coef*lb(x) */
         assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
         SCIProwprepAddConstant(rowprep, coef * SCIPvarGetLbLocal(var));
         rowprep->local |= SCIPisGT(scip, SCIPvarGetLbLocal(var), SCIPvarGetLbGlobal(var));
      }
      else
      {
         assert((coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT));
         /* we aggregate with -coef*var >= -coef*ub(x) */
         assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
         SCIProwprepAddConstant(rowprep, coef * SCIPvarGetUbLocal(var));
         rowprep->local |= SCIPisLT(scip, SCIPvarGetUbLocal(var), SCIPvarGetUbGlobal(var));
      }

      /* eliminate coef*var from rowprep: remove coef */
      if( pos == 0 )
      {
         /* set first term to zero */
         rowprep->coefs[maxcoefidx] = 0.0;

         /* update index */
         ++maxcoefidx;

         /* update maxcoef */
         maxcoef = REALABS(rowprep->coefs[maxcoefidx]);
      }
      else
      {
         /* forget last term */
         --rowprep->nvars;

         /* update mincoef */
         mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
      }

      /* (potentially) remember the variable that has been removed here */
      SCIP_CALL( rowprepRecordModifiedVar(scip, rowprep, var) );
   }

   /* if maximal coefs were removed, then there are now 0's in the beginning of the coefs array
    * -> move all remaining coefs and vars up front
    */
   if( maxcoefidx > 0 )
   {
      int i;
      for( i = maxcoefidx; i < rowprep->nvars; ++i )
      {
         rowprep->vars[i-maxcoefidx] = rowprep->vars[i];
         rowprep->coefs[i-maxcoefidx] = rowprep->coefs[i];
      }
      rowprep->nvars -= maxcoefidx;
   }

   return SCIP_OKAY;
}


/** scales up rowprep if it seems useful */
static
void rowprepCleanupScaleup(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_Real*            viol,               /**< violation of cut in sol (input and output) */
   SCIP_Real             minviol             /**< minimal violation we try to achieve */
   )
{
   SCIP_Real scalefactor;
   SCIP_Real mincoef;
   SCIP_Real maxcoef;

   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(viol != NULL);

   /* if violation is very small than better don't scale up */
   if( *viol < ROWPREP_SCALEUP_VIOLNONZERO )
      return;

   /* if violation is already above minviol, then nothing to do */
   if( *viol >= minviol )
      return;
   assert(!SCIPisInfinity(scip, *viol));

   /* if violation is sufficiently positive (>10*eps), but has not reached minviol,
    * then consider scaling up to reach approx MINVIOLFACTOR*minviol
    */
   scalefactor = ROWPREP_SCALEUP_MINVIOLFACTOR * minviol / *viol;

   /* scale by approx. scalefactor, if minimal coef is not so large yet and maximal coef and rhs don't get huge by doing so (or have been so before) */
   mincoef = rowprep->nvars > 0 ? REALABS(rowprep->coefs[rowprep->nvars-1]) : 1.0;
   maxcoef = rowprep->nvars > 0 ? REALABS(rowprep->coefs[0]) : 1.0;
   if( mincoef < ROWPREP_SCALEUP_MAXMINCOEF && scalefactor * maxcoef < ROWPREP_SCALEUP_MAXMAXCOEF && scalefactor * REALABS(rowprep->side) < ROWPREP_SCALEUP_MAXSIDE )
   {
      int scaleexp;

      /* SCIPinfoMessage(scip, NULL, "scale up by ~%g, viol=%g: ", scalefactor, myviol);
         SCIPprintRowprep(scip, rowprep, NULL); */

      /* SCIPscaleRowprep returns the actually applied scale factor */
      scaleexp = SCIPscaleRowprep(rowprep, scalefactor);
      *viol = ldexp(*viol, scaleexp);

      /* SCIPinfoMessage(scip, NULL, "scaled up by %g, viol=%g: ", ldexp(1.0, scaleexp), myviol);
         SCIPprintRowprep(scip, rowprep, NULL); */
   }
}

/** scales down rowprep if it improves coefs and keeps rowprep violated */
static
void rowprepCleanupScaledown(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_Real*            viol,               /**< violation of cut in sol (input and output) */
   SCIP_Real             minviol             /**< minimal violation we try to keep */
   )
{
   SCIP_Real scalefactor;

   /* if maxcoef < ROWPREP_SCALEDOWN_MINMAXCOEF (or no terms), then don't consider scaling down */
   if( rowprep->nvars == 0 || REALABS(rowprep->coefs[0]) < ROWPREP_SCALEDOWN_MINMAXCOEF )
      return;

   /* consider scaling down so that maxcoef ~ 10 */
   scalefactor = 10.0 / REALABS(rowprep->coefs[0]);

   /* if minimal violation would be lost by scaling down, then increase scalefactor such that minviol is still reached */
   if( *viol > minviol && !SCIPisInfinity(scip, *viol) && scalefactor * *viol < minviol )
   {
      assert(minviol > 0.0);  /* since viol >= 0, the if-condition should ensure that minviol > 0 */
      assert(*viol > 0.0);    /* since minviol > 0, the if-condition ensures viol > 0 */
      scalefactor = ROWPREP_SCALEUP_MINVIOLFACTOR * minviol / *viol;
   }

   /* scale by approx. scalefactor if scaling down and minimal coef does not get too small
    * myviol < minviol (-> scalefactor > 1) or mincoef < feastol before scaling is possible, in which case we also don't scale down
    */
   if( scalefactor < 1.0 && scalefactor * REALABS(rowprep->coefs[rowprep->nvars-1]) > ROWPREP_SCALEDOWN_MINCOEF )
   {
      int scaleexp;

      /* SCIPinfoMessage(scip, NULL, "scale down by ~%g, viol=%g: ", scalefactor, myviol);
         SCIPprintRowprep(scip, rowprep, NULL); */

      scaleexp = SCIPscaleRowprep(rowprep, scalefactor);
      if( !SCIPisInfinity(scip, *viol) )
         *viol = ldexp(*viol, scaleexp);

      /* SCIPinfoMessage(scip, NULL, "scaled down by %g, viol=%g: ", ldexp(1.0, scaleexp), myviol);
         SCIPprintRowprep(scip, rowprep, NULL); */
   }
}

/** rounds almost integral coefs to integrals, thereby trying to relax the cut */
static
SCIP_RETCODE rowprepCleanupIntegralCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_Real*            viol                /**< NULL or violation of cut in sol (input), set to SCIP_INVALID if some coef changed */
   )
{
   SCIP_Real coef;
   SCIP_Real roundcoef;
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);

   /* Coefficients smaller than epsilon are rounded to 0.0 when added to row and
    * coefficients very close to integral values are rounded to integers when added to LP.
    * Both cases can be problematic if variable value is very large (bad numerics).
    * Thus, we anticipate by rounding coef here, but also modify constant so that cut is still valid (if possible),
    * i.e., bound coef[i]*x by round(coef[i])*x + (coef[i]-round(coef[i])) * bound(x).
    * Or in other words, we aggregate with the variable bound.
    *
    * If the required bound of x is not finite, then only round coef (introduces an error).
    * @TODO If only the opposite bound is available, then one could move the coefficient
    *   away from the closest integer so that the SCIP_ROW won't try to round it.
    */
   for( i = 0; i < rowprep->nvars; ++i )
   {
      coef = rowprep->coefs[i];
      roundcoef = SCIPround(scip, coef);
      if( coef != roundcoef && SCIPisEQ(scip, coef, roundcoef) ) /*lint !e777*/
      {
         SCIP_Real xbnd;
         SCIP_VAR* var;

         var = rowprep->vars[i];
         if( rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
            if( rowprep->local )
               xbnd = coef > roundcoef ? SCIPvarGetLbLocal(var)  : SCIPvarGetUbLocal(var);
            else
               xbnd = coef > roundcoef ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var);
         else
            if( rowprep->local )
               xbnd = coef > roundcoef ? SCIPvarGetUbLocal(var)  : SCIPvarGetLbLocal(var);
            else
               xbnd = coef > roundcoef ? SCIPvarGetUbGlobal(var) : SCIPvarGetLbGlobal(var);

         if( !SCIPisInfinity(scip, REALABS(xbnd)) )
         {
            /* if there is a bound, then relax row side so rounding coef will not introduce an error */
            SCIPdebugMsg(scip, "var <%s> [%g,%g] has almost integral coef %.15g, round coefficient to %g and add constant %g\n",
               SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), coef, roundcoef, (coef-roundcoef) * xbnd);
            SCIProwprepAddConstant(rowprep, (coef-roundcoef) * xbnd);
         }
         else
         {
            /* if there is no bound, then we make the coef integral, too, even though this will introduce an error
             * however, SCIP_ROW would do this anyway, but doing this here might eliminate some epsilon coefs (so they don't determine mincoef below)
             * and helps to get a more accurate row violation value
             */
            SCIPdebugMsg(scip, "var <%s> [%g,%g] has almost integral coef %.15g, round coefficient to %g without relaxing side (!)\n",
               SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), coef, roundcoef);
         }
         rowprep->coefs[i] = roundcoef;
         if( viol != NULL )
            *viol = SCIP_INVALID;

         /* (potentially) remember the variable which coef has been modified here */
         SCIP_CALL( rowprepRecordModifiedVar(scip, rowprep, var) );
      }
   }

   /* forget about coefs that became exactly zero by the above step */
   while( rowprep->nvars > 0 && rowprep->coefs[rowprep->nvars-1] == 0.0 )
      --rowprep->nvars;

   return SCIP_OKAY;
}

/** relaxes almost zero side */
static
void rowprepCleanupSide(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_Real*            viol                /**< NULL or violation of cut in sol (input), set to SCIP_INVALID if some coef changed */
   )
{
   /* SCIP_ROW handling will replace a side close to 0 by 0.0, even if that makes the row more restrictive
    * we thus relax the side here so that it will either be 0 now or will not be rounded to 0 later
    */
   if( rowprep->side == 0.0 || !SCIPisZero(scip, rowprep->side) )
      return;

   if( rowprep->side > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
      rowprep->side =  1.1*SCIPepsilon(scip);
   else if( rowprep->side < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT )
      rowprep->side = -1.1*SCIPepsilon(scip);
   else
      rowprep->side = 0.0;

   if( rowprep->recordmodifications )
      rowprep->modifiedside = TRUE;

   if( viol != NULL )
      *viol = SCIP_INVALID;
}

#ifdef NDEBUG
/* Undo the defines from pub_misc_rowprep.h, which exist if NDEBUG is defined. */
#undef SCIProwprepGetNVars
#undef SCIProwprepGetVars
#undef SCIProwprepGetCoefs
#undef SCIProwprepGetSide
#undef SCIProwprepGetSidetype
#undef SCIProwprepIsLocal
#undef SCIProwprepGetName
#undef SCIProwprepGetNModifiedVars
#undef SCIProwprepGetModifiedVars
#undef SCIProwprepAddSide
#undef SCIProwprepAddConstant
#undef SCIProwprepSetSidetype
#undef SCIProwprepSetLocal
#undef SCIProwprepRecordModifications
#endif

/** creates a SCIP_ROWPREP datastructure
 *
 * Initial row represents 0 &le; 0.
 */
SCIP_RETCODE SCIPcreateRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep,            /**< buffer to store pointer to rowprep */
   SCIP_SIDETYPE         sidetype,           /**< whether cut will be or lower-equal or larger-equal type */
   SCIP_Bool             local               /**< whether cut will be valid only locally */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, rowprep) );
   BMSclearMemory(*rowprep);

   (*rowprep)->sidetype = sidetype;
   (*rowprep)->local = local;

   return SCIP_OKAY;
}

/** frees a SCIP_ROWPREP datastructure */
void SCIPfreeRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep             /**< pointer that stores pointer to rowprep */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(*rowprep != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*rowprep)->vars, (*rowprep)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*rowprep)->coefs, (*rowprep)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*rowprep)->modifiedvars, (*rowprep)->modifiedvarssize);
   SCIPfreeBlockMemory(scip, rowprep);
}

/** creates a copy of a SCIP_ROWPREP datastructure */
SCIP_RETCODE SCIPcopyRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        target,             /**< buffer to store pointer of rowprep copy */
   SCIP_ROWPREP*         source              /**< rowprep to copy */
   )
{
   assert(scip != NULL);
   assert(target != NULL);
   assert(source != NULL);

   SCIP_CALL( SCIPduplicateBlockMemory(scip, target, source) );
   if( source->coefs != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*target)->coefs, source->coefs, source->varssize) );
   }
   if( source->vars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*target)->vars, source->vars, source->varssize) );
   }

   (*target)->recordmodifications = FALSE;
   (*target)->modifiedvars = NULL;
   (*target)->modifiedvarssize = 0;
   (*target)->nmodifiedvars = 0;
   (*target)->modifiedside = FALSE;

   return SCIP_OKAY;
}

/** gives number of terms in rowprep */
int SCIProwprepGetNVars(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->nvars;
}

/** gives variables of rowprep (feel free to modify) */
SCIP_VAR** SCIProwprepGetVars(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->vars;
}

/** gives coefficients of rowprep (feel free to modify) */
SCIP_Real* SCIProwprepGetCoefs(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->coefs;
}

/** gives side of rowprep */
SCIP_Real SCIProwprepGetSide(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->side;
}

/** gives kind of inequality of rowprep */
SCIP_SIDETYPE SCIProwprepGetSidetype(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->sidetype;
}

/** returns whether rowprep is locally valid only */
SCIP_Bool SCIProwprepIsLocal(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->local;
}

/** returns name of rowprep (feel free to modify) */
char* SCIProwprepGetName(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->name;
}

/** returns number of variables which coefficients were modified in cleanup */
int SCIProwprepGetNModifiedVars(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->nmodifiedvars;
}

/** returns variables which coefficients were modified in cleanup */
SCIP_VAR** SCIProwprepGetModifiedVars(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   return rowprep->modifiedvars;
}

/** resets rowprep to have 0 terms and side 0.0 */
void SCIProwprepReset(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   rowprep->nvars = 0;
   rowprep->side = 0.0;

   rowprep->recordmodifications = FALSE;
   rowprep->nmodifiedvars = 0;
   rowprep->modifiedside = FALSE;
}

#ifdef NDEBUG
#undef SCIProwprepAddSide
#undef SCIProwprepAddConstant
#endif

/** adds constant value to side of rowprep */
void SCIProwprepAddSide(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Real             side                /**< constant value to be added to side */
   )
{
   assert(rowprep != NULL);

   rowprep->side += side;
}

/** adds constant term to rowprep
 *
 * Substracts constant from side.
 */
void SCIProwprepAddConstant(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Real             constant            /**< constant value to be added */
   )
{
   SCIProwprepAddSide(rowprep, -constant);
}

/** sets side type of rowprep */
void SCIProwprepSetSidetype(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_SIDETYPE         sidetype            /**< new side type */
   )
{
   assert(rowprep != NULL);

   rowprep->sidetype = sidetype;
}

/** sets whether rowprep is local */
void SCIProwprepSetLocal(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Bool             islocal             /**< whether rowprep is local */
   )
{
   assert(rowprep != NULL);

   rowprep->local = islocal;
}

/** enables recording for where modifications were done in cleanup */
void SCIProwprepRecordModifications(
   SCIP_ROWPREP*         rowprep             /**< rowprep */
   )
{
   assert(rowprep != NULL);

   rowprep->recordmodifications = TRUE;
}

/** prints a rowprep */
void SCIPprintRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be printed */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   )
{
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);

   if( *rowprep->name != '\0' )
   {
      SCIPinfoMessage(scip, file, "[%s](%c) ", rowprep->name, rowprep->local ? 'l' : 'g');
   }

   for( i = 0; i < rowprep->nvars; ++i )
   {
      SCIPinfoMessage(scip, file, "%+.15g*<%s> ", rowprep->coefs[i], SCIPvarGetName(rowprep->vars[i]));
   }

   SCIPinfoMessage(scip, file, rowprep->sidetype == SCIP_SIDETYPE_LEFT ? ">= %.15g\n" : "<= %.15g\n", rowprep->side);
}

/** prints a rowprep and values in solution */
void SCIPprintRowprepSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be printed */
   SCIP_SOL*             sol,                /**< solution for activity */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   )
{
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_Real term;
   SCIP_Real maxterm;
   SCIP_Real activity;
   SCIP_Real violation;
   int maxtermidx;
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);

   if( *rowprep->name != '\0' )
   {
      SCIPinfoMessage(scip, file, "[%s](%c) ", rowprep->name, rowprep->local ? 'l' : 'g');
   }

   activity = 0.0;
   maxterm = REALABS(rowprep->side);
   maxtermidx = -1;
   for( i = 0; i < rowprep->nvars; ++i )
   {
      coef = rowprep->coefs[i];
      var = rowprep->vars[i];
      SCIPinfoMessage(scip, file, "%+.15g*<%s>(%.15g) ", coef, SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var));

      term = coef * SCIPgetSolVal(scip, sol, var);
      if( REALABS(term) > maxterm )
      {
         maxterm = term;
         maxtermidx = i;
      }

      activity += term;
   }

   SCIPinfoMessage(scip, file, rowprep->sidetype == SCIP_SIDETYPE_LEFT ? ">= %.15g" : "<= %.15g", rowprep->side);

   if( rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
      /* cut is activity <= side -> violation is activity - side (if positive) */
      violation = activity - rowprep->side;
   else
      /* cut is activity >= side -> violation is side - activity (if positive) */
      violation = rowprep->side - activity;

   SCIPinfoMessage(scip, file, "; activity %.15g", activity);
   SCIPinfoMessage(scip, file, "; violation %e", violation);
   SCIPinfoMessage(scip, file, "; maxterm %e at pos %d\n", maxterm, maxtermidx);
}

/** ensures that rowprep has space for at least given number of additional terms
 *
 * Useful when knowing in advance how many terms will be added.
 */
SCIP_RETCODE SCIPensureRowprepSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   int                   size                /**< number of additional terms for which to alloc space in rowprep */
   )
{
   int oldsize;

   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(size >= 0);

   if( rowprep->varssize >= rowprep->nvars + size )
      return SCIP_OKAY;  /* already enough space left */

   /* realloc vars and coefs array */
   oldsize = rowprep->varssize;
   rowprep->varssize = SCIPcalcMemGrowSize(scip, rowprep->nvars + size);

   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &rowprep->vars,  oldsize, rowprep->varssize) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &rowprep->coefs, oldsize, rowprep->varssize) );

   return SCIP_OKAY;
}

/** adds a term coef*var to a rowprep */
SCIP_RETCODE SCIPaddRowprepTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             coef                /**< coefficient to add */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(var != NULL);

   if( coef == 0.0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, 1) );
   assert(rowprep->varssize > rowprep->nvars);

   rowprep->vars[rowprep->nvars] = var;
   rowprep->coefs[rowprep->nvars] = coef;
   ++rowprep->nvars;

   return SCIP_OKAY;
}

/** adds several terms coef*var to a rowprep */
SCIP_RETCODE SCIPaddRowprepTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   int                   nvars,              /**< number of terms to add */
   SCIP_VAR**            vars,               /**< variables to add */
   SCIP_Real*            coefs               /**< coefficients to add */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(vars != NULL || nvars == 0);
   assert(coefs != NULL || nvars == 0);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nvars) );
   assert(rowprep->varssize >= rowprep->nvars + nvars);

   /*lint --e{866} */
   BMScopyMemoryArray(rowprep->vars + rowprep->nvars, vars, nvars);
   BMScopyMemoryArray(rowprep->coefs + rowprep->nvars, coefs, nvars);
   rowprep->nvars += nvars;

   return SCIP_OKAY;
}

/** computes violation of rowprep in a given solution
 *
 * Can return whether the violation value is reliable from a floating-point accuracy point of view.
 * The value will not be deemed reliable when its calculation involved the subtraction of large numbers.
 * To be precise, the violation of an inequality \f$ \sum_i a_ix_i \leq b \f$ in a solution \f$x^*\f$ is deemed
 * reliable if \f$ |\sum_i a_ix^*_i - b| \geq 2^{-50} \max (|b|, \max_i |a_ix^*_i|) \f$.
 */
SCIP_Real SCIPgetRowprepViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_SOL*             sol,                /**< solution or NULL for LP solution */
   SCIP_Bool*            reliable            /**< buffer to store whether computed violation is reliable (numerically), or NULL if not of interest */
   )
{
   SCIP_Real activity;
   SCIP_Real maxterm;
   SCIP_Real term;
   SCIP_Real violation;
   SCIP_Real val;
   int i;

   activity = 0.0;
   maxterm = REALABS(rowprep->side);
   for( i = 0; i < rowprep->nvars; ++i )
   {
      /* Loose variable have the best bound as LP solution value.
       * HOWEVER, they become column variables when they are added to a row (via SCIPaddVarsToRow below).
       * When this happens, their LP solution value changes to 0.0!
       * So when calculating the row activity for an LP solution, we treat loose variable as if they were already column variables.
       */
      if( sol != NULL || SCIPvarGetStatus(rowprep->vars[i]) != SCIP_VARSTATUS_LOOSE )
      {
         val = SCIPgetSolVal(scip, sol, rowprep->vars[i]);

         /* If a variable is at infinity, then this should lead to an immediate decision.
          * Having different contradicting infinities is something I would now know how to handle and am ignoring now.
          */
         if( SCIPisInfinity(scip, val * (rowprep->coefs[i] >= 0.0 ? 1.0 : -1.0)) )
         {
            /* activity = SCIPinfinity(scip); */
            if( reliable != NULL )
               *reliable = TRUE;
            if( rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
               return SCIPinfinity(scip);  /* infinity <= side -> always violated */
            else
               return 0.0;  /* infinity >= side -> never violated */
         }
         if( SCIPisInfinity(scip, val * (rowprep->coefs[i] >= 0.0 ? -1.0 : 1.0)) )
         {
            /* activity = -SCIPinfinity(scip); */
            if( reliable != NULL )
               *reliable = TRUE;
            if( rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
               return 0.0;  /* -infinity <= side -> never violated */
            else
               return SCIPinfinity(scip);  /* -infinity >= side -> always violated */
         }

         term = rowprep->coefs[i] * val;
         activity += term;

         if( reliable != NULL && REALABS(term) > maxterm )
            maxterm = REALABS(term);
      }
   }

   if( rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
      /* cut is activity <= side -> violation is activity - side (if positive) */
      violation = activity - rowprep->side;
   else
      /* cut is activity >= side -> violation is side - activity (if positive) */
      violation = rowprep->side - activity;

   /* In double precision, the mantissa (or significand) of a floating point number has 52 bit.
    * Therefore, if the exponent in the violation is 52 (or more) less than the one of maxterm,
    * then it is essentially random.
    * We require here that the exponents differ by at most 50.
    * To be more robust w.r.t. scaling of the row, we look at the exponent of the quotient maxterm/violation
    * instead of the difference of the exponents of maxterm and violation.
    */
   if( reliable != NULL )
   {
      if( violation != 0.0 )
      {
         int exponent;
         (void) frexp(maxterm / violation, &exponent);  /* difference in exponents for maxterm and violation */
         *reliable = exponent <= 50;
      }
      else
         *reliable = TRUE;  /* not clear how to evaluate reliability here, so think positive */
   }

   return MAX(violation, 0.0);
}

/** computes violation of rowprep in a given solution and reports whether that value seem numerically reliable
 *
 * @see SCIPgetRowprepViolation()
 */
SCIP_Bool SCIPisRowprepViolationReliable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_SOL*             sol                 /**< solution or NULL for LP solution */
   )
{
   SCIP_Bool reliable;

   assert(scip != NULL);
   assert(rowprep != NULL);

   (void) SCIPgetRowprepViolation(scip, rowprep, sol, &reliable);

   return reliable;
}

/** Merge terms that use same variable and eliminate zero coefficients.
 *
 * Removes a variable if its bounds have a relative difference of below epsilon.
 * Local bounds are checked for local rows, otherwise global bounds are used.
 * If the bounds are not absolute equal, the bound that relaxes the row is used.
 *
 * Terms are sorted by variable (see SCIPvarComp()) after return.
 */
void SCIPmergeRowprepTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep             /**< rowprep to be cleaned up */
   )
{
   int i;
   int j;

   assert(scip != NULL);
   assert(rowprep != NULL);

   if( rowprep->nvars <= 1 )
      return;

   /* sort terms by variable index */
   SCIPsortPtrReal((void**)rowprep->vars, rowprep->coefs, SCIPvarComp, rowprep->nvars);

   /* merge terms with same variable, drop 0 coefficients */
   i = 0;
   j = 1;
   while( j < rowprep->nvars )
   {
      if( rowprep->vars[i] == rowprep->vars[j] )
      {
         /* merge term j into term i */
         rowprep->coefs[i] += rowprep->coefs[j];
         ++j;
         continue;
      }

      /* move term i into side if fixed */
      if( rowprep->local && SCIPisRelEQ(scip, SCIPvarGetLbLocal(rowprep->vars[i]), SCIPvarGetUbLocal(rowprep->vars[i])) )
      {
         if( (rowprep->coefs[i] > 0.0) == (rowprep->sidetype == SCIP_SIDETYPE_RIGHT) )
            rowprep->side -= rowprep->coefs[i] * SCIPvarGetLbLocal(rowprep->vars[i]);
         else
            rowprep->side -= rowprep->coefs[i] * SCIPvarGetUbLocal(rowprep->vars[i]);
         rowprep->coefs[i] = 0.0;  /* so will be cleaned out below */
      }
      else if( !rowprep->local && SCIPisRelEQ(scip, SCIPvarGetLbGlobal(rowprep->vars[i]), SCIPvarGetUbGlobal(rowprep->vars[i])) )
      {
         if( (rowprep->coefs[i] > 0.0) == (rowprep->sidetype == SCIP_SIDETYPE_RIGHT) )
            rowprep->side -= rowprep->coefs[i] * SCIPvarGetLbGlobal(rowprep->vars[i]);
         else
            rowprep->side -= rowprep->coefs[i] * SCIPvarGetUbGlobal(rowprep->vars[i]);
         rowprep->coefs[i] = 0.0;  /* so will be cleaned out below */
      }

      if( rowprep->coefs[i] == 0.0 )
      {
         /* move term j to position i */
         rowprep->coefs[i] = rowprep->coefs[j];
         rowprep->vars[i] = rowprep->vars[j];
         ++j;
         continue;
      }

      /* move term j to position i+1 and move on */
      if( j != i+1 )
      {
         rowprep->vars[i+1] = rowprep->vars[j];
         rowprep->coefs[i+1] = rowprep->coefs[j];
      }
      ++i;
      ++j;
   }

   /* move term i into side if fixed */
   if( rowprep->local && SCIPisRelEQ(scip, SCIPvarGetLbLocal(rowprep->vars[i]), SCIPvarGetUbLocal(rowprep->vars[i])) )
   {
      if( (rowprep->coefs[i] > 0.0) == (rowprep->sidetype == SCIP_SIDETYPE_RIGHT) )
         rowprep->side -= rowprep->coefs[i] * SCIPvarGetLbLocal(rowprep->vars[i]);
      else
         rowprep->side -= rowprep->coefs[i] * SCIPvarGetUbLocal(rowprep->vars[i]);
      rowprep->coefs[i] = 0.0;  /* so will be cleaned out below */
   }
   else if( !rowprep->local && SCIPisRelEQ(scip, SCIPvarGetLbGlobal(rowprep->vars[i]), SCIPvarGetUbGlobal(rowprep->vars[i])) )
   {
      if( (rowprep->coefs[i] > 0.0) == (rowprep->sidetype == SCIP_SIDETYPE_RIGHT) )
         rowprep->side -= rowprep->coefs[i] * SCIPvarGetLbGlobal(rowprep->vars[i]);
      else
         rowprep->side -= rowprep->coefs[i] * SCIPvarGetUbGlobal(rowprep->vars[i]);
      rowprep->coefs[i] = 0.0;  /* so will be cleaned out below */
   }

   /* remaining term can have coef zero -> forget about it */
   if( rowprep->coefs[i] == 0.0 )
      --i;

   /* i points to last term */
   rowprep->nvars = i+1;
}

/** Cleans up and attempts to improve rowprep
 *
 * Drops small or large coefficients if coefrange is too large, if this can be done by relaxing the row.
 * Scales coefficients up to reach minimal violation, if possible.
 * Scaling is omitted if violation is very small (\ref ROWPREP_SCALEUP_VIOLNONZERO) or
 * maximal coefficient would become huge (\ref ROWPREP_SCALEUP_MAXMAXCOEF).
 * Scales coefficients and side down if they are large and if the minimal violation is still reached.
 * Rounds coefficients close to integral values to integrals, if this can be done by relaxing the row.
 * Rounds side within epsilon of 0 to 0.0 or +/-1.1*epsilon, whichever relaxes the row least.
 *
 * After return, the terms in the rowprep will be sorted by absolute value of coefficient, in decreasing order.
 * Thus, the coefrange can be obtained via `REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1])` (if nvars>0).
 *
 * `success` is set to TRUE if and only if the rowprep satisfies the following:
 * - the coefrange is below `maxcoefrange`
 * - the violation is at least `minviol`
 * - the violation is reliable or `minviol` = 0
 * - the absolute value of coefficients are below SCIPinfinity()
 * - the absolute value of the side is below SCIPinfinity()
 */
SCIP_RETCODE SCIPcleanupRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be cleaned */
   SCIP_SOL*             sol,                /**< solution that we try to cut off, or NULL for LP solution */
   SCIP_Real             minviol,            /**< minimal absolute violation the row should achieve (w.r.t. sol) */
   SCIP_Real*            viol,               /**< buffer to store absolute violation of cleaned up cut in sol, or NULL if not of interest */
   SCIP_Bool*            success             /**< buffer to store whether cut cleanup was successful, or NULL if not of interest */
   )
{
   SCIP_Real myviol;
   SCIP_Bool violreliable = TRUE;
   SCIP_Real maxcoefrange;
#ifdef SCIP_DEBUG
   SCIP_Real mincoef = 1.0;
   SCIP_Real maxcoef = 1.0;
#endif

   maxcoefrange = SCIPsetGetSepaMaxCoefRatioRowprep(scip->set);

   if( rowprep->recordmodifications )
   {
      /* forget about possible previous modifications */
      rowprep->nmodifiedvars = 0;
      rowprep->modifiedside = FALSE;
   }

   /* sort term by absolute value of coef. */
   SCIP_CALL( rowprepCleanupSortTerms(scip, rowprep) );

#ifdef SCIP_DEBUG
   if( rowprep->nvars > 0 )
   {
      maxcoef = REALABS(rowprep->coefs[0]);
      mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
   }

   SCIPinfoMessage(scip, NULL, "starting cleanup, coefrange %g: ", maxcoef/mincoef);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* improve coefficient range by aggregating out variables */
   SCIP_CALL( rowprepCleanupImproveCoefrange(scip, rowprep, sol, maxcoefrange) );

   /* get current violation in sol (reliability info only needed if success is not NULL) */
   myviol = SCIPgetRowprepViolation(scip, rowprep, sol, success != NULL ? &violreliable : NULL);  /*lint !e826*/
   assert(myviol >= 0.0);

#ifdef SCIP_DEBUG
   if( rowprep->nvars > 0 )
   {
      maxcoef = REALABS(rowprep->coefs[0]);
      mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
   }

   SCIPinfoMessage(scip, NULL, "improved coefrange to %g, viol %g: ", maxcoef / mincoef, myviol);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* if there is interest in achieving some minimal violation, then possibly scale up to increase violation
    * this updates myviol; since this is only scaling the cut, it doesn't change anything about the reliability of the violation value */
   if( minviol > 0.0 )
   {
      /* first, try to achieve scip's minefficacy (typically 1e-4) */
      if( SCIPgetSepaMinEfficacy(scip) > minviol )
         rowprepCleanupScaleup(scip, rowprep, &myviol, SCIPgetSepaMinEfficacy(scip));
      /* in case scip minefficacy could not be reached or was smaller than minviol, try with the given minviol */
      rowprepCleanupScaleup(scip, rowprep, &myviol, minviol);
   }

   /* scale down to improve numerics, updates myviol (reliability doesn't change) */
   rowprepCleanupScaledown(scip, rowprep, &myviol, MAX(SCIPgetSepaMinEfficacy(scip), minviol)); /*lint !e666*/

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "applied scaling, viol %g: ", myviol);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* turn almost-integral coefs to integral values, may set myviol to SCIP_INVALID */
   SCIP_CALL( rowprepCleanupIntegralCoefs(scip, rowprep, &myviol) );

   /* relax almost-zero side, may set myviol to SCIP_INVALID */
   rowprepCleanupSide(scip, rowprep, &myviol);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "adjusted almost-integral coefs and sides, viol %g: ", myviol);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

#if !1
   /* compute final coefrange, if requested by caller */
   if( coefrange != NULL )
   {
      if( rowprep->nvars > 0 )
         *coefrange = REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1]);
      else
         *coefrange = 1.0;
   }
#endif

   /* check whether rowprep could be turned into a reasonable row */
   if( success != NULL )
   {
      *success = TRUE;

      /* check whether the coef.range is below maxcoefrange */
      if( rowprep->nvars > 0 && REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1]) > maxcoefrange )
      {
         SCIPdebugMsg(scip, "rowprep coefrange %g is above the limit %g\n", REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1]), maxcoefrange);
         *success = FALSE;
      }

      /* check whether coefficients are below SCIPinfinity (terms are order by coef value) */
      if( *success && rowprep->nvars > 0 && SCIPisInfinity(scip, REALABS(rowprep->coefs[0])) )
      {
         SCIPdebugMsg(scip, "rowprep coefficient %g is beyond value for infinity\n", rowprep->coefs[0]);
         *success = FALSE;
      }

      /* check whether the absolute value of the side is below SCIPinfinity */
      if( *success && SCIPisInfinity(scip, REALABS(rowprep->side)) )
      {
         SCIPdebugMsg(scip, "rowprep side %g is beyond value for infinity\n", rowprep->side);
         *success = FALSE;
      }

      /* check if violation is at least minviol and reliable, if minviol > 0 */
      if( *success && minviol > 0.0 )
      {
         /* may need to recompute violation if coefs or side was modified above */
         if( myviol == SCIP_INVALID )  /*lint !e777 */
            myviol = SCIPgetRowprepViolation(scip, rowprep, sol, &violreliable);

         if( !violreliable )
         {
            SCIPdebugMsg(scip, "rowprep violation %g is not reliable\n", myviol);
            *success = FALSE;
         }
         else if( myviol < minviol )
         {
            SCIPdebugMsg(scip, "rowprep violation %g is below minimal violation %g\n", myviol, minviol);
            *success = FALSE;
         }
      }
   }

   /* If we updated myviol correctly, then it should coincide with freshly computed violation.
    * I leave this assert off for now, since getting the tolerance in the EQ correctly is not trivial. We recompute viol below anyway.
    */
   /* assert(myviol == SCIP_INVALID || SCIPisEQ(scip, myviol, SCIPgetRowprepViolation(scip, rowprep, sol, NULL))); */

   /* compute final violation, if requested by caller */
   if( viol != NULL )  /*lint --e{777} */
      *viol = myviol == SCIP_INVALID ? SCIPgetRowprepViolation(scip, rowprep, sol, NULL) : myviol;

   return SCIP_OKAY;
}

/** Cleans up and attempts to improve rowprep without regard for violation
 *
 * Drops small or large coefficients if their ratio is beyond separating/maxcoefratiofacrowprep / numerics/feastol,
 * if this can be done by relaxing the row.
 * Scales coefficients and side to have maximal coefficient in `[1/maxcoefbound,maxcoefbound]`.
 * Rounds coefficients close to integral values to integrals, if this can be done by relaxing the row.
 * Rounds side within epsilon of 0 to 0.0 or +/-1.1*epsilon, whichever relaxes the row least.
 *
 * After return, the terms in the rowprep will be sorted by absolute value of coefficient, in decreasing order.
 * Thus, the coefratio can be obtained via `REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1])` (if nvars>0).
 *
 * `success` is set to TRUE if and only if the rowprep satisfies the following:
 * - the coefratio is below separating/maxcoefratiofacrowprep / numerics/feastol
 * - the absolute value of coefficients are below SCIPinfinity()
 * - the absolute value of the side is below SCIPinfinity()
 *
 * In difference to SCIPcleanupRowprep(), this function does not scale up the row to increase the absolute violation.
 */
SCIP_RETCODE SCIPcleanupRowprep2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be cleaned */
   SCIP_SOL*             sol,                /**< solution that we try to cut off, or NULL for LP solution */
   SCIP_Real             maxcoefbound,       /**< bound on absolute value of largest coefficient */
   SCIP_Bool*            success             /**< buffer to store whether cut cleanup was successful, or NULL if not of interest */
   )
{
   SCIP_Real maxcoefrange;
#ifdef SCIP_DEBUG
   SCIP_Real mincoef = 1.0;
   SCIP_Real maxcoef = 1.0;
#endif

   assert(maxcoefbound >= 1.0);

   maxcoefrange = SCIPsetGetSepaMaxCoefRatioRowprep(scip->set);

   if( rowprep->recordmodifications )
   {
      /* forget about possible previous modifications */
      rowprep->nmodifiedvars = 0;
      rowprep->modifiedside = FALSE;
   }

   /* sort term by absolute value of coef. */
   SCIP_CALL( rowprepCleanupSortTerms(scip, rowprep) );

#ifdef SCIP_DEBUG
   if( rowprep->nvars > 0 )
   {
      maxcoef = REALABS(rowprep->coefs[0]);
      mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
   }

   SCIPinfoMessage(scip, NULL, "starting cleanup, coefrange %g: ", maxcoef/mincoef);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* improve coefficient range by aggregating out variables */
   SCIP_CALL( rowprepCleanupImproveCoefrange(scip, rowprep, sol, maxcoefrange) );

#ifdef SCIP_DEBUG
   if( rowprep->nvars > 0 )
   {
      maxcoef = REALABS(rowprep->coefs[0]);
      mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
   }

   SCIPinfoMessage(scip, NULL, "improved coefrange to %g: ", maxcoef / mincoef);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* scale up or down to improve numerics
    * if maximal coef is below 1.0/maxcoefbound, scale up to reach ~ 1.0/maxcoefbound
    * if maximal coef is above maxcoefbound, scale down to ~ maxcoefbound
    */
   if( rowprep->nvars > 0 && !SCIPisInfinity(scip, maxcoefbound) )
   {
      SCIP_Real expon = 0.0;
      if( REALABS(rowprep->coefs[0]) < 1.0/maxcoefbound )
         expon = SCIPscaleRowprep(rowprep, (1.0/maxcoefbound) / REALABS(rowprep->coefs[0]));
      else if( REALABS(rowprep->coefs[0]) > maxcoefbound )
         expon = SCIPscaleRowprep(rowprep, maxcoefbound / REALABS(rowprep->coefs[0]));

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "applied scaling by %g: ", pow(2.0, expon));
      SCIPprintRowprep(scip, rowprep, NULL);
#else
      (void) expon;
#endif
   }

   /* turn almost-integral coefs to integral values */
   SCIP_CALL( rowprepCleanupIntegralCoefs(scip, rowprep, NULL) );

   /* relax almost-zero side */
   rowprepCleanupSide(scip, rowprep, NULL);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "adjusted almost-integral coefs and sides: ");
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* check whether rowprep could be turned into a reasonable row */
   if( success != NULL )
   {
      *success = TRUE;

      /* check whether the coef.range is below maxcoefrange */
      if( rowprep->nvars > 0 && REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1]) > maxcoefrange )
      {
         SCIPdebugMsg(scip, "rowprep coefrange %g is above the limit %g\n", REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1]), maxcoefrange);
         *success = FALSE;
      }

      /* check whether coefficients are below SCIPinfinity (terms are order by coef value) */
      if( *success && rowprep->nvars > 0 && SCIPisInfinity(scip, REALABS(rowprep->coefs[0])) )
      {
         SCIPdebugMsg(scip, "rowprep coefficient %g is beyond value for infinity\n", rowprep->coefs[0]);
         *success = FALSE;
      }

      /* check whether the absolute value of the side is below SCIPinfinity */
      if( *success && SCIPisInfinity(scip, REALABS(rowprep->side)) )
      {
         SCIPdebugMsg(scip, "rowprep side %g is beyond value for infinity\n", rowprep->side);
         *success = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** Scales up a rowprep to increase coefficients/sides that are within epsilon to an integer value, if possible.
 *
 * Computes the minimal fractionality of all fractional coefficients and the side of the rowprep.
 * If this fractionality is below epsilon, the rowprep is scaled up such that the fractionality exceeds epsilon,
 * if this will not put any coefficient or side above SCIPhugeValue().
 *
 * This function does not relax the rowprep.
 *
 * `success` is set to TRUE if the resulting rowprep can be turned into a SCIP_ROW, that is,
 * all coefs and the side is below SCIPinfinity() and fractionalities are above epsilon.
 * If `success` is set to FALSE, then the rowprep will not have been modified.
 *
 * @return The applied scaling factor, if `success` is set to TRUE.
 */
SCIP_Real SCIPscaleupRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be cleaned */
   SCIP_Real             minscaleup,         /**< minimal factor by which to scale up row, or <= 1.0 if to be ignored */
   SCIP_Bool*            success             /**< buffer to store whether rowprep could be turned into SCIP_ROW without loss, or NULL if not of interest */
   )
{
   SCIP_Real minfrac = 0.5;
   SCIP_Real minfrac0 = 0.5;
   SCIP_Real frac;
   SCIP_Real maxval;
   SCIP_Real factor = 1.0;
   SCIP_Bool makeintegral = TRUE;
   int i;

   /* find the smallest fractionality in rowprep sides and coefficients and the largest absolute coefficient/side */
   frac = REALABS(floor(rowprep->side + 0.5) - rowprep->side);
   if( frac != 0.0 )
   {
      if( REALABS(rowprep->side) > 0.5 )
      {
         if( frac < minfrac )
            minfrac = frac;
      }
      else if( frac < minfrac0 )
         minfrac0 = frac;
   }
   maxval = REALABS(rowprep->side);

   for( i = 0; i < rowprep->nvars; ++i )
   {
      frac = REALABS(floor(rowprep->coefs[i] + 0.5) - rowprep->coefs[i]);
      if( frac != 0.0 )
      {
         if( REALABS(rowprep->coefs[i]) > 0.5 )
         {
            if( frac < minfrac )
               minfrac = frac;
         }
         else if( frac < minfrac0 )
            minfrac0 = frac;
      }
      if( REALABS(rowprep->coefs[i]) > maxval )
         maxval = REALABS(rowprep->coefs[i]);
   }

   SCIPdebugMsg(scip, "minimal fractional of rowprep coefs and side is %g, max coef/side is %g\n", MIN(minfrac, minfrac0), maxval);

   /* in order for SCIP_ROW to not modify the coefs and side, they need to be more than epsilon way from an integer value
    *
    * If the integer value is 0, then scaling up the rowprep by epsilon/minfrac will increase its minimal fractionality
    * above epsilon.
    * If the integer value is not zero, then scaling up the rowprep by a well-chosen fractional number alpha will increase
    * the minimal fractionality by about alpha*integer-value mod 1. To reduce the chance that alpha*integer-value is integral,
    * we use a "very fractional" value for alpha.
    *
    * If the scaling increases the maximal coef/value beyond SCIPinfinity, then the rowprep would be useless.
    * We even check that we don't increase beyond SCIPhugeValue here
    */
   if( minfrac0 <= SCIPepsilon(scip) )
   {
      factor = 1.1 * SCIPepsilon(scip) / minfrac0;

      if( factor < minscaleup )
         factor = minscaleup;
   }
   else if( minfrac <= SCIPepsilon(scip) )
   {
      factor = MAX(M_SQRT2, minscaleup);
      makeintegral = FALSE;
   }
   else if( minscaleup > 1.0 )
   {
      factor = minscaleup;
   }
   else
   {
      /* do not scale up, only check whether maxval is already below infinity */
      if( success != NULL )
         *success = !SCIPisInfinity(scip, maxval);

      return 1.0;
   }

   if( !SCIPisHugeValue(scip, factor * maxval) )
   {
      if( makeintegral)
      {
         factor = SCIPscaleRowprep(rowprep, factor);

#ifdef SCIP_DEBUG
         factor = pow(2.0, factor);  /* SCIPscaleRowprep() actually returned log2 of factor */
#endif
      }
      else
      {
         /* multiply each coefficient by factor */
         for( i = 0; i < rowprep->nvars; ++i )
            rowprep->coefs[i] *= factor;

         /* multiply side by factor */
         rowprep->side *= factor;
      }
#ifdef SCIP_DEBUG
      maxval *= factor;
      SCIPinfoMessage(scip, NULL, "scaled up rowprep by %g (minfrac=%g, minscaleup=%g), maxval is now %g\n", factor, minfrac, minscaleup, maxval);
      SCIPprintRowprep(scip, rowprep, NULL);
#endif

      if( success != NULL )
         *success = TRUE;
   }
   else if( success != NULL )
      *success = FALSE;

   return factor;
}

/** scales a rowprep by given factor (after some rounding)
 *
 * @return Exponent of actually applied scaling factor, if written as \f$2^x\f$.
 */
int SCIPscaleRowprep(
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be scaled */
   SCIP_Real             factor              /**< suggested scale factor */
   )
{
   double v;
   int expon;
   int i;

   assert(rowprep != NULL);
   assert(factor > 0.0);

   /* write factor as v*2^expon with v in [0.5,1) */
   v = frexp(factor, &expon);
   /* adjust to v'*2^expon with v' in (0.5,1] by v'=v if v > 0.5, v'=1 if v=0.5 */
   if( v == 0.5 )
      --expon;

   /* multiply each coefficient by 2^expon */
   for( i = 0; i < rowprep->nvars; ++i )
      rowprep->coefs[i] = ldexp(rowprep->coefs[i], expon);

   /* multiply side by 2^expon */
   rowprep->side = ldexp(rowprep->side, expon);

   return expon;
}

/** generates a SCIP_ROW from a rowprep, setting its origin to given constraint handler */
SCIP_RETCODE SCIPgetRowprepRowConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(scip != NULL);
   assert(row != NULL);
   assert(rowprep != NULL);
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, row, conshdlr, rowprep->name,
      rowprep->sidetype == SCIP_SIDETYPE_LEFT  ? rowprep->side : -SCIPinfinity(scip),
      rowprep->sidetype == SCIP_SIDETYPE_RIGHT ? rowprep->side :  SCIPinfinity(scip),
      rowprep->local && (SCIPgetDepth(scip) > 0), FALSE, TRUE) );

   SCIP_CALL( SCIPaddVarsToRow(scip, *row, rowprep->nvars, rowprep->vars, rowprep->coefs) );

   return SCIP_OKAY;
}

/** generates a SCIP_ROW from a rowprep, setting its origin to given constraint */
SCIP_RETCODE SCIPgetRowprepRowCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(scip != NULL);
   assert(row != NULL);
   assert(rowprep != NULL);
   assert(cons != NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, row, cons, rowprep->name,
      rowprep->sidetype == SCIP_SIDETYPE_LEFT  ? rowprep->side : -SCIPinfinity(scip),
      rowprep->sidetype == SCIP_SIDETYPE_RIGHT ? rowprep->side :  SCIPinfinity(scip),
      rowprep->local && (SCIPgetDepth(scip) > 0), FALSE, TRUE) );

   SCIP_CALL( SCIPaddVarsToRow(scip, *row, rowprep->nvars, rowprep->vars, rowprep->coefs) );

   return SCIP_OKAY;
}

/** generates a SCIP_ROW from a rowprep, setting its origin to given separator */
SCIP_RETCODE SCIPgetRowprepRowSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(scip != NULL);
   assert(row != NULL);
   assert(rowprep != NULL);

   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, row, sepa, rowprep->name,
      rowprep->sidetype == SCIP_SIDETYPE_LEFT  ? rowprep->side : -SCIPinfinity(scip),
      rowprep->sidetype == SCIP_SIDETYPE_RIGHT ? rowprep->side :  SCIPinfinity(scip),
      rowprep->local && (SCIPgetDepth(scip) > 0), FALSE, TRUE) );

   SCIP_CALL( SCIPaddVarsToRow(scip, *row, rowprep->nvars, rowprep->vars, rowprep->coefs) );

   return SCIP_OKAY;
}
