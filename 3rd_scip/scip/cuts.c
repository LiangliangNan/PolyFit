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

/**@file   cuts.c
 * @brief  methods for aggregation of rows
 *
 * @author Jakob Witzig
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cuts.h"
#include "scip/misc.h"
#include "scip/struct_lp.h"
#include "scip/lp.h"
#include "scip/struct_cuts.h"
#include "scip/cons_knapsack.h"
#include "scip/struct_scip.h"
#include "scip/dbldblarith.h"

/* =========================================== general static functions =========================================== */
#ifdef SCIP_DEBUG
static
void printCutQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            cutcoefs,           /**< non-zero coefficients of cut */
   QUAD(SCIP_Real        cutrhs),            /**< right hand side of the MIR row */
   int*                  cutinds,            /**< indices of problem variables for non-zero coefficients */
   int                   cutnnz,             /**< number of non-zeros in cut */
   SCIP_Bool             ignorsol,
   SCIP_Bool             islocal
   )
{
   SCIP_Real QUAD(activity);
   SCIP_VAR** vars;
   int i;

   assert(scip != NULL);
   vars = SCIPgetVars(scip);

   SCIPdebugMessage("CUT:");
   QUAD_ASSIGN(activity, 0.0);
   for( i = 0; i < cutnnz; ++i )
   {
      SCIP_Real QUAD(coef);

      QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]);

      SCIPdebugPrintf(" %+g<%s>", QUAD_TO_DBL(coef), SCIPvarGetName(vars[cutinds[i]]));

      if( !ignorsol )
      {
         SCIPquadprecProdQD(coef, coef, (sol == NULL ? SCIPvarGetLPSol(vars[cutinds[i]]) : SCIPgetSolVal(scip, sol, vars[cutinds[i]])));
      }
      else
      {
         if( cutcoefs[i] > 0.0 )
         {
            SCIPquadprecProdQD(coef, coef, (islocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]])));
         }
         else
         {
            SCIPquadprecProdQD(coef, coef, (islocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]])));
         }
      }

      SCIPquadprecSumQQ(activity, activity, coef);
   }
   SCIPdebugPrintf(" <= %.6f (activity: %g)\n", QUAD_TO_DBL(cutrhs), QUAD_TO_DBL(activity));
}
#endif

/* macro to make sure a value is not equal to zero, i.e. NONZERO(x) != 0.0
 * will be TRUE for every x including 0.0
 *
 * To avoid branches it will add 1e-100 with the same sign as x to x which will
 * be rounded away for any sane non-zero value but will make sure the value is
 * never exactly 0.0
 */
#define NONZERO(x)   (COPYSIGN(1e-100, (x)) + (x))

/** add a scaled row to a dense vector indexed over the problem variables and keep the
 *  index of non-zeros up-to-date
 */
static
SCIP_RETCODE varVecAddScaledRowCoefs(
   int*RESTRICT          inds,               /**< pointer to array with variable problem indices of non-zeros in variable vector */
   SCIP_Real*RESTRICT    vals,               /**< array with values of variable vector */
   int*RESTRICT          nnz,                /**< number of non-zeros coefficients of variable vector */
   SCIP_ROW*             row,                /**< row coefficients to add to variable vector */
   SCIP_Real             scale               /**< scale for adding given row to variable vector */
   )
{
   /* add up coefficients */
   int i;

   assert(inds != NULL);
   assert(vals != NULL);
   assert(nnz != NULL);
   assert(row != NULL);

   /* add the non-zeros to the aggregation row and keep non-zero index up to date */
   for( i = 0 ; i < row->len; ++i )
   {
      SCIP_Real val;
      int probindex = row->cols[i]->var_probindex;

      val = vals[probindex];

      if( val == 0.0 )
         inds[(*nnz)++] = probindex;

      val += row->vals[i] * scale;

      /* the value must not be exactly zero due to sparsity pattern */
      val = NONZERO(val);

      assert(val != 0.0);
      vals[probindex] = val;
   }

   return SCIP_OKAY;
}

/** add a scaled row to a dense vector indexed over the problem variables and keep the
 *  index of non-zeros up-to-date
 */
static
SCIP_RETCODE varVecAddScaledRowCoefsQuad(
   int*RESTRICT          inds,               /**< pointer to array with variable problem indices of non-zeros in variable vector */
   SCIP_Real*RESTRICT    vals,               /**< array with values of variable vector */
   int*RESTRICT          nnz,                /**< number of non-zeros coefficients of variable vector */
   SCIP_ROW*             row,                /**< row coefficients to add to variable vector */
   SCIP_Real             scale               /**< scale for adding given row to variable vector */
   )
{
   /* add up coefficients */
   int i;

   assert(inds != NULL);
   assert(vals != NULL);
   assert(nnz != NULL);
   assert(row != NULL);

   /* add the non-zeros to the aggregation row and keep non-zero index up to date */
   for( i = 0 ; i < row->len; ++i )
   {
      SCIP_Real QUAD(val);
      int probindex = row->cols[i]->var_probindex;

      QUAD_ARRAY_LOAD(val, vals, probindex);

      if( QUAD_HI(val) == 0.0 )
         inds[(*nnz)++] = probindex;

      SCIPquadprecSumQD(val, val, row->vals[i] * scale);

      /* the value must not be exactly zero due to sparsity pattern */
      QUAD_HI(val) = NONZERO(QUAD_HI(val));
      assert(QUAD_HI(val) != 0.0);

      QUAD_ARRAY_STORE(vals, probindex, val);
   }

   return SCIP_OKAY;
}

/* calculates the cuts efficacy for the given solution */
static
SCIP_Real calcEfficacy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to calculate the efficacy for (NULL for LP solution) */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   SCIP_Real             cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int                   cutnnz              /**< the number of non-zeros in the cut */
   )
{
   SCIP_VAR** vars;
   SCIP_Real norm;
   SCIP_Real activity;
   int i;

   assert(scip != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);

   vars = SCIPgetVars(scip);

   activity = 0.0;
   for( i = 0; i < cutnnz; ++i )
      activity += cutcoefs[i] * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);

   norm = SCIPgetVectorEfficacyNorm(scip, cutcoefs, cutnnz);
   return (activity - cutrhs) / MAX(1e-6, norm);
}

static
SCIP_Real calcEfficacyNormQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            vals,               /**< array of the non-zero coefficients in the vector; this is a quad precision array! */
   int*                  inds,               /**< array of the problem indices of variables with a non-zero coefficient in the vector */
   int                   nnz                 /**< the number of non-zeros in the vector */
   )
{
   SCIP_Real norm;
   SCIP_Real QUAD(coef);
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   norm = 0.0;
   switch( scip->set->sepa_efficacynorm )
   {
   case 'e':
      for( i = 0; i < nnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, vals, inds[i]);
         norm += SQR(QUAD_TO_DBL(coef));
      }
      norm = SQRT(norm);
      break;
   case 'm':
      for( i = 0; i < nnz; ++i )
      {
         SCIP_Real absval;
         QUAD_ARRAY_LOAD(coef, vals, inds[i]);

         absval = REALABS(QUAD_TO_DBL(coef));
         norm = MAX(norm, absval);
      }
      break;
   case 's':
      for( i = 0; i < nnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, vals, inds[i]);
         norm += REALABS(QUAD_TO_DBL(coef));
      }
      break;
   case 'd':
      for( i = 0; i < nnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, vals, inds[i]);
         if( !SCIPisZero(scip, QUAD_TO_DBL(coef)) )
         {
            norm = 1.0;
            break;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", scip->set->sepa_efficacynorm);
      assert(FALSE);
   }

   return norm;
}

/* calculates the cuts efficacy for the given solution; the cut coefs are stored densely and in quad precision */
static
SCIP_Real calcEfficacyDenseStorageQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to calculate the efficacy for (NULL for LP solution) */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut; this is a quad precision array! */
   SCIP_Real             cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int                   cutnnz              /**< the number of non-zeros in the cut */
   )
{
   SCIP_VAR** vars;
   SCIP_Real norm;
   SCIP_Real activity;
   SCIP_Real QUAD(coef);
   int i;

   assert(scip != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
   assert(scip->set != NULL);

   vars = SCIPgetVars(scip);

   activity = 0.0;
   for( i = 0; i < cutnnz; ++i )
   {
      QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]);
      activity += QUAD_TO_DBL(coef) * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
   }

   norm = calcEfficacyNormQuad(scip, cutcoefs, cutinds, cutnnz);
   return (activity - cutrhs) / MAX(1e-6, norm);
}

/** safely remove all coefficients below the given value; returns TRUE if the cut became redundant */
static
SCIP_Bool removeZerosQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             minval,             /**< minimal absolute value of coefficients that should not be removed */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   QUAD(SCIP_Real*       cutrhs),            /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz              /**< the number of non-zeros in the cut */
   )
{
   int i;
   SCIP_VAR** vars;

   vars = SCIPgetVars(scip);

   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real QUAD(val);

      int v = cutinds[i];
      QUAD_ARRAY_LOAD(val, cutcoefs, v);

      if( EPSZ(QUAD_TO_DBL(val), minval) )
      {
         if( REALABS(QUAD_TO_DBL(val)) > QUAD_EPSILON )
         {
            /* adjust left and right hand sides with max contribution */
            if( QUAD_TO_DBL(val) < 0.0 )
            {
               SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[v]) : SCIPvarGetUbGlobal(vars[v]);
               if( SCIPisInfinity(scip, ub) )
                  return TRUE;
               else
               {
                  SCIPquadprecProdQD(val, val, ub);
                  SCIPquadprecSumQQ(*cutrhs, *cutrhs, -val);
               }
            }
            else
            {
               SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[v]) : SCIPvarGetLbGlobal(vars[v]);
               if( SCIPisInfinity(scip, -lb) )
                  return TRUE;
               else
               {
                  SCIPquadprecProdQD(val, val, lb);
                  SCIPquadprecSumQQ(*cutrhs, *cutrhs, -val);
               }
            }
         }

         QUAD_ASSIGN(val, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, val);

         /* remove non-zero entry */
         --(*cutnnz);
         cutinds[i] = cutinds[*cutnnz];
      }
      else
         ++i;
   }

   /* relax rhs to zero, if it's very close to */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   return FALSE;
}

/** safely remove all coefficients below the given value; returns TRUE if the cut became redundant */
static
SCIP_Bool removeZeros(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             minval,             /**< minimal absolute value of coefficients that should not be removed */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   QUAD(SCIP_Real*       cutrhs),            /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz              /**< the number of non-zeros in the cut */
   )
{
   int i;
   SCIP_VAR** vars;

   vars = SCIPgetVars(scip);

   /* loop over non-zeros and remove values below minval; values above QUAD_EPSILON are cancelled with their bound
    * to avoid numerical rounding errors
    */
   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real val;

      int v = cutinds[i];
      val = cutcoefs[v];

      if( EPSZ(val, minval) )
      {
         if( REALABS(val) > QUAD_EPSILON )
         {
            /* adjust left and right hand sides with max contribution */
            if( val < 0.0 )
            {
               SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[v]) : SCIPvarGetUbGlobal(vars[v]);
               if( SCIPisInfinity(scip, ub) )
                  return TRUE;
               else
               {
                  SCIPquadprecSumQD(*cutrhs, *cutrhs, -val * ub);
               }
            }
            else
            {
               SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[v]) : SCIPvarGetLbGlobal(vars[v]);
               if( SCIPisInfinity(scip, -lb) )
                  return TRUE;
               else
               {
                  SCIPquadprecSumQD(*cutrhs, *cutrhs, -val * lb);
               }
            }
         }

         cutcoefs[v] = 0.0;

         /* remove non-zero entry */
         --(*cutnnz);
         cutinds[i] = cutinds[*cutnnz];
      }
      else
         ++i;
   }

   /* relax rhs to zero, if it's very close to */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   return FALSE;
}

static
SCIP_DECL_SORTINDCOMP(compareAbsCoefsQuad)
{
   SCIP_Real abscoef1;
   SCIP_Real abscoef2;
   SCIP_Real QUAD(coef1);
   SCIP_Real QUAD(coef2);
   SCIP_Real* coefs = (SCIP_Real*) dataptr;

   QUAD_ARRAY_LOAD(coef1, coefs, ind1);
   QUAD_ARRAY_LOAD(coef2, coefs, ind2);

   abscoef1 = REALABS(QUAD_TO_DBL(coef1));
   abscoef2 = REALABS(QUAD_TO_DBL(coef2));

   if( abscoef1 < abscoef2 )
      return -1;
   if( abscoef2 < abscoef1 )
      return 1;

   return 0;
}

static
SCIP_DECL_SORTINDCOMP(compareAbsCoefs)
{
   SCIP_Real abscoef1;
   SCIP_Real abscoef2;
   SCIP_Real* coefs = (SCIP_Real*) dataptr;

   abscoef1 = REALABS(coefs[ind1]);
   abscoef2 = REALABS(coefs[ind2]);

   if( abscoef1 < abscoef2 )
      return -1;
   if( abscoef2 < abscoef1 )
      return 1;

   return 0;
}

/** change given coefficient to new given value, adjust right hand side using the variables bound;
 *  returns TRUE if the right hand side would need to be changed to infinity and FALSE otherwise
 */
static
SCIP_Bool chgCoeffWithBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable the coefficient belongs to */
   SCIP_Real             oldcoeff,           /**< old coefficient value */
   SCIP_Real             newcoeff,           /**< new coefficient value */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   QUAD(SCIP_Real*       cutrhs)             /**< pointer to adjust right hand side of cut */
   )
{
   SCIP_Real QUAD(delta);
   SCIPquadprecSumDD(delta, newcoeff, -oldcoeff);

   if( QUAD_TO_DBL(delta) > QUAD_EPSILON )
   {
      SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
      if( SCIPisInfinity(scip, ub) )
         return TRUE;
      else
      {
         SCIPquadprecProdQD(delta, delta, ub);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, delta);
      }
   }
   else if( QUAD_TO_DBL(delta) < -QUAD_EPSILON )
   {
      SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
      if( SCIPisInfinity(scip, -lb) )
         return TRUE;
      else
      {
         SCIPquadprecProdQD(delta, delta, lb);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, delta);
      }
   }

   return FALSE;
}

/** change given (quad) coefficient to new given value, adjust right hand side using the variables bound;
 *  returns TRUE if the right hand side would need to be changed to infinity and FALSE otherwise
 */
static
SCIP_Bool chgQuadCoeffWithBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable the coefficient belongs to */
   QUAD(SCIP_Real        oldcoeff),          /**< old coefficient value */
   SCIP_Real             newcoeff,           /**< new coefficient value */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   QUAD(SCIP_Real*       cutrhs)             /**< pointer to adjust right hand side of cut */
   )
{
   SCIP_Real QUAD(delta);

   SCIPquadprecSumQD(delta, -oldcoeff, newcoeff);

   if( QUAD_TO_DBL(delta) > QUAD_EPSILON )
   {
      SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
      if( SCIPisInfinity(scip, ub) )
         return TRUE;
      else
      {
         SCIPquadprecProdQD(delta, delta, ub);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, delta);
      }
   }
   else if( QUAD_TO_DBL(delta) < -QUAD_EPSILON )
   {
      SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
      if( SCIPisInfinity(scip, -lb) )
         return TRUE;
      else
      {
         SCIPquadprecProdQD(delta, delta, lb);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, delta);
      }
   }

   return FALSE;
}


/** scales the cut and then  tighten the coefficients of the given cut based on the maximal activity;
 *  see cons_linear.c consdataTightenCoefs() for details; the cut is given in a semi-sparse quad precision array;
 */
static
SCIP_RETCODE cutTightenCoefsQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   QUAD(SCIP_Real*       cutrhs),            /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< the number of non-zeros in the cut */
   SCIP_Bool*            redundant           /**< whether the cut was detected to be redundant */
   )
{
   int i;
   int nintegralvars;
   SCIP_Bool isintegral;
   SCIP_VAR** vars;
   SCIP_Real QUAD(maxacttmp);
   SCIP_Real maxact;
   SCIP_Real maxabsintval;
   SCIP_Real maxabscontval;

   QUAD_ASSIGN(maxacttmp, 0.0);

   vars = SCIPgetVars(scip);
   maxabsintval = 0.0;
   maxabscontval = 0.0;
   nintegralvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   isintegral = TRUE;

   *redundant = FALSE;

   /* compute the maximum activity and maximum absolute coefficient values for all and for integral variables in the cut */
   for( i = 0; i < *cutnnz; ++i )
   {
      SCIP_Real QUAD(val);

      assert(cutinds[i] >= 0);
      assert(vars[cutinds[i]] != NULL);

      QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);

      if( QUAD_TO_DBL(val) < 0.0 )
      {
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, -lb) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, -QUAD_TO_DBL(val));
         else
         {
            maxabscontval = MAX(maxabscontval, -QUAD_TO_DBL(val));
            isintegral = FALSE;
         }

         SCIPquadprecProdQD(val, val, lb);

         SCIPquadprecSumQQ(maxacttmp, maxacttmp, val);
      }
      else
      {
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, ub) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, QUAD_TO_DBL(val));
         else
         {
            maxabscontval = MAX(maxabscontval, QUAD_TO_DBL(val));
            isintegral = FALSE;
         }

         SCIPquadprecProdQD(val, val, ub);

         SCIPquadprecSumQQ(maxacttmp, maxacttmp, val);
      }
   }

   maxact = QUAD_TO_DBL(maxacttmp);

   /* cut is redundant in activity bounds */
   if( SCIPisFeasLE(scip, maxact, QUAD_TO_DBL(*cutrhs)) )
   {
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   /* cut is only on integral variables, try to scale to integral coefficients */
   if( isintegral )
   {
      SCIP_Real equiscale;
      SCIP_Real intscalar;
      SCIP_Bool success;
      SCIP_Real* intcoeffs;

      SCIP_CALL( SCIPallocBufferArray(scip, &intcoeffs, *cutnnz) );

      equiscale = 1.0 / MIN((maxact - QUAD_TO_DBL(*cutrhs)), maxabsintval);

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(val);

         QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);
         SCIPquadprecProdQD(val, val, equiscale);

         intcoeffs[i] = QUAD_TO_DBL(val);
      }

      SCIP_CALL( SCIPcalcIntegralScalar(intcoeffs, *cutnnz, -SCIPsumepsilon(scip), SCIPepsilon(scip),
            (SCIP_Longint)scip->set->sepa_maxcoefratio, scip->set->sepa_maxcoefratio, &intscalar, &success) );

      SCIPfreeBufferArray(scip, &intcoeffs);

      if( success )
      {
         /* if successful, apply the scaling */
         intscalar *= equiscale;

         SCIPquadprecProdQD(*cutrhs, *cutrhs, intscalar);

         for( i = 0; i < *cutnnz; )
         {
            SCIP_Real QUAD(val);
            SCIP_Real intval;

            QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);
            SCIPquadprecProdQD(val, val, intscalar);

            intval = SCIPround(scip, QUAD_TO_DBL(val));

            if( chgQuadCoeffWithBound(scip, vars[cutinds[i]], QUAD(val), intval, cutislocal, QUAD(cutrhs)) )
            {
               /* TODO maybe change the coefficient to the other value instead of discarding the cut? */
               *redundant = TRUE;
               return SCIP_OKAY;
            }

            if( intval != 0.0 )
            {
               QUAD_ASSIGN(val, intval);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], val);
               ++i;
            }
            else
            {
               /* this must not be -0.0, otherwise the clean buffer memory is not cleared properly */
               QUAD_ASSIGN(val, 0.0);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], val);
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
            }
         }

         SCIPquadprecEpsFloorQ(*cutrhs, *cutrhs, SCIPfeastol(scip)); /*lint !e666*/

         /* recompute the maximal activity after scaling to integral values */
         QUAD_ASSIGN(maxacttmp, 0.0);
         maxabsintval = 0.0;

         for( i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real QUAD(val);

            assert(cutinds[i] >= 0);
            assert(vars[cutinds[i]] != NULL);

            QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);

            if( QUAD_TO_DBL(val) < 0.0 )
            {
               SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, -QUAD_TO_DBL(val));

               SCIPquadprecProdQD(val, val, lb);

               SCIPquadprecSumQQ(maxacttmp, maxacttmp, val);
            }
            else
            {
               SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, QUAD_TO_DBL(val));

               SCIPquadprecProdQD(val, val, ub);

               SCIPquadprecSumQQ(maxacttmp, maxacttmp, val);
            }
         }

         maxact = QUAD_TO_DBL(maxacttmp);

         assert(EPSISINT(maxact, 1e-4));
         maxact = SCIPround(scip, maxact);
         QUAD_ASSIGN(maxacttmp, maxact);

         /* check again for redundancy */
         if( SCIPisFeasLE(scip, maxact, QUAD_TO_DBL(*cutrhs)) )
         {
            *redundant = TRUE;
            return SCIP_OKAY;
         }
      }
      else
      {
         /* otherwise, apply the equilibrium scaling */
         isintegral = FALSE;

         /* perform the scaling */
         SCIPquadprecProdQD(maxacttmp, maxacttmp, equiscale);

         SCIPquadprecProdQD(*cutrhs, *cutrhs, equiscale);
         maxabsintval *= equiscale;

         for( i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real QUAD(val);

            QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);
            SCIPquadprecProdQD(val, val, equiscale);
            QUAD_ARRAY_STORE(cutcoefs, cutinds[i], val);
         }
      }
   }
   else
   {
      /* cut has integer and continuous variables, so scale it to equilibrium */
      SCIP_Real scale;
      SCIP_Real maxabsval;

      maxabsval = maxact - QUAD_TO_DBL(*cutrhs);
      maxabsval = MIN(maxabsval, maxabsintval);
      maxabsval = MAX(maxabsval, maxabscontval);

      scale = 1.0 / maxabsval; /*lint !e795*/

      /* perform the scaling */
      SCIPquadprecProdQD(maxacttmp, maxacttmp, scale);
      maxact = QUAD_TO_DBL(maxacttmp);

      SCIPquadprecProdQD(*cutrhs, *cutrhs, scale);
      maxabsintval *= scale;

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(val);

         QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);
         SCIPquadprecProdQD(val, val, scale);
         QUAD_ARRAY_STORE(cutcoefs, cutinds[i], val);
      }
   }

   /* no coefficient tightening can be performed since the precondition doesn't hold for any of the variables */
   if( SCIPisGT(scip, maxact - maxabsintval, QUAD_TO_DBL(*cutrhs)) )
      return SCIP_OKAY;

   SCIPsortDownInd(cutinds, compareAbsCoefsQuad, (void*) cutcoefs, *cutnnz);

   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real QUAD(val);

      if( cutinds[i] >= nintegralvars )
      {
         ++i;
         continue;
      }

      QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);

      assert(SCIPvarIsIntegral(vars[cutinds[i]]));

      if( QUAD_TO_DBL(val) < 0.0 && SCIPisLE(scip, maxact + QUAD_TO_DBL(val), QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, *cutrhs, -maxacttmp);

         if( isintegral )
         {
            /* if the cut is integral, the true coefficient must also be integral;
             * thus we round it to the exact integral value
             */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) > QUAD_TO_DBL(val) )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQQ(delta, -val, coef);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);

            SCIPdebugPrintf("tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
                   QUAD_TO_DBL(val), QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp), lb,
                   cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]));

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(!SCIPisPositive(scip, QUAD_TO_DBL(coef)));

            if( SCIPisNegative(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
            }
            else
            {
               QUAD_ASSIGN(coef, 0.0);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else if( QUAD_TO_DBL(val) > 0.0 && SCIPisLE(scip, maxact - QUAD_TO_DBL(val), QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, maxacttmp, -*cutrhs);

         if( isintegral )
         {
            /* if the cut is integral, the true coefficient must also be integral;
             * thus we round it to the exact integral value
             */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) < QUAD_TO_DBL(val) )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQQ(delta, -val, coef);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);

            SCIPdebugPrintf("tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
                   QUAD_TO_DBL(val), QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp),
                   cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]), ub);

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(SCIPisGE(scip, QUAD_TO_DBL(coef), 0.0));
            if( SCIPisPositive(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
            }
            else
            {
               QUAD_ASSIGN(coef, 0.0);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else /* due to sorting we can stop completely if the precondition was not fulfilled for this variable */
         break;

      ++i;
   }

   return SCIP_OKAY;
}

/** scales the cut and then  tighten the coefficients of the given cut based on the maximal activity;
 *  see cons_linear.c consdataTightenCoefs() for details; the cut is given in a semi-sparse array;
 */
static
SCIP_RETCODE cutTightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   QUAD(SCIP_Real*       cutrhs),            /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< the number of non-zeros in the cut */
   SCIP_Bool*            redundant           /**< pointer to return whtether the cut was detected to be redundant */
   )
{
   int i;
   int nintegralvars;
   SCIP_Bool isintegral;
   SCIP_VAR** vars;
   SCIP_Real QUAD(maxacttmp);
   SCIP_Real maxact;
   SCIP_Real maxabsintval;
   SCIP_Real maxabscontval;

   QUAD_ASSIGN(maxacttmp, 0.0);

   vars = SCIPgetVars(scip);
   maxabsintval = 0.0;
   maxabscontval = 0.0;
   isintegral = TRUE;
   *redundant = FALSE;
   nintegralvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

   for( i = 0; i < *cutnnz; ++i )
   {
      SCIP_Real val;

      assert(cutinds[i] >= 0);
      assert(vars[cutinds[i]] != NULL);

      val = cutcoefs[cutinds[i]];

      if( val < 0.0 )
      {
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, -lb) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, -val);
         else
         {
            maxabscontval = MAX(maxabscontval, -val);
            isintegral = FALSE;
         }

         SCIPquadprecSumQD(maxacttmp, maxacttmp, val * lb);
      }
      else
      {
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, ub) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, val);
         else
         {
            maxabscontval = MAX(maxabscontval, val);
            isintegral = FALSE;
         }

         SCIPquadprecSumQD(maxacttmp, maxacttmp, val * ub);
      }
   }

   maxact = QUAD_TO_DBL(maxacttmp);

   /* cut is redundant in activity bounds */
   if( SCIPisFeasLE(scip, maxact, QUAD_TO_DBL(*cutrhs)) )
   {
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   /* cut is only on integral variables, try to scale to integral coefficients */
   if( isintegral )
   {
      SCIP_Real equiscale;
      SCIP_Real intscalar;
      SCIP_Bool success;
      SCIP_Real* intcoeffs;

      SCIP_CALL( SCIPallocBufferArray(scip, &intcoeffs, *cutnnz) );

      equiscale = 1.0 / MIN((maxact - QUAD_TO_DBL(*cutrhs)), maxabsintval);

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real scaleval;
         SCIP_Real val;

         val = cutcoefs[cutinds[i]];

         scaleval = val * equiscale;

         intcoeffs[i] = scaleval;
      }

      SCIP_CALL( SCIPcalcIntegralScalar(intcoeffs, *cutnnz, -SCIPsumepsilon(scip), SCIPepsilon(scip),
            (SCIP_Longint)scip->set->sepa_maxcoefratio, scip->set->sepa_maxcoefratio, &intscalar, &success) );

      SCIPfreeBufferArray(scip, &intcoeffs);

      if( success )
      {
         /* if successful, apply the scaling */
         intscalar *= equiscale;

         SCIPquadprecProdQD(*cutrhs, *cutrhs, intscalar);

         for( i = 0; i < *cutnnz; )
         {
            SCIP_Real val;
            SCIP_Real intval;

            val = cutcoefs[cutinds[i]];
            val *= intscalar;

            intval = SCIPround(scip, val);

            if( chgCoeffWithBound(scip, vars[cutinds[i]], val, intval, cutislocal, QUAD(cutrhs)) )
            {
               /* TODO maybe change the coefficient to the other value instead of discarding the cut? */
               *redundant = TRUE;
               return SCIP_OKAY;
            }

            if( intval != 0.0 )
            {
               cutcoefs[cutinds[i]] = intval;
               ++i;
            }
            else
            {
               /* this must not be -0.0, otherwise the clean buffer memory is not cleared properly */
               cutcoefs[cutinds[i]] = 0.0;
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
            }
         }

         SCIPquadprecEpsFloorQ(*cutrhs, *cutrhs, SCIPfeastol(scip)); /*lint !e666*/

         /* recompute the maximal activity after scaling to integral values */
         QUAD_ASSIGN(maxacttmp, 0.0);
         maxabsintval = 0.0;

         for( i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real val;

            assert(cutinds[i] >= 0);
            assert(vars[cutinds[i]] != NULL);

            val = cutcoefs[cutinds[i]];

            if( val < 0.0 )
            {
               SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, -val);

               val *= lb;

               SCIPquadprecSumQD(maxacttmp, maxacttmp, val);
            }
            else
            {
               SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, val);

               val *= ub;

               SCIPquadprecSumQD(maxacttmp, maxacttmp, val);
            }
         }

         maxact = QUAD_TO_DBL(maxacttmp);

         assert(EPSISINT(maxact, 1e-4));
         maxact = SCIPround(scip, maxact);
         QUAD_ASSIGN(maxacttmp, maxact);

         /* check again for redundancy */
         if( SCIPisFeasLE(scip, maxact, QUAD_TO_DBL(*cutrhs)) )
         {
            *redundant = TRUE;
            return SCIP_OKAY;
         }
      }
      else
      {
         /* otherwise, apply the equilibrium scaling */
         isintegral = FALSE;

         /* perform the scaling */
         SCIPquadprecProdQD(maxacttmp, maxacttmp, equiscale);

         SCIPquadprecProdQD(*cutrhs, *cutrhs, equiscale);
         maxabsintval *= equiscale;

         for( i = 0; i < *cutnnz; ++i )
            cutcoefs[cutinds[i]] *= equiscale;
      }
   }
   else
   {
      /* cut has integer and continuous variables, so scale it to equilibrium */
      SCIP_Real scale;
      SCIP_Real maxabsval;

      maxabsval = maxact - QUAD_TO_DBL(*cutrhs);
      maxabsval = MIN(maxabsval, maxabsintval);
      maxabsval = MAX(maxabsval, maxabscontval);

      scale = 1.0 / maxabsval; /*lint !e795*/

      /* perform the scaling */
      SCIPquadprecProdQD(maxacttmp, maxacttmp, scale);
      maxact = QUAD_TO_DBL(maxacttmp);

      SCIPquadprecProdQD(*cutrhs, *cutrhs, scale);
      maxabsintval *= scale;

      for( i = 0; i < *cutnnz; ++i )
         cutcoefs[cutinds[i]] *= scale;
   }

   /* no coefficient tightening can be performed since the precondition doesn't hold for any of the variables */
   if( SCIPisGT(scip, maxact - maxabsintval, QUAD_TO_DBL(*cutrhs)) )
      return SCIP_OKAY;

   SCIPsortDownInd(cutinds, compareAbsCoefs, (void*) cutcoefs, *cutnnz);

   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real val;

      if( cutinds[i] >= nintegralvars )
      {
         ++i;
         continue;
      }

      val = cutcoefs[cutinds[i]];

      assert(SCIPvarIsIntegral(vars[cutinds[i]]));

      if( val < 0.0 && SCIPisLE(scip, maxact + val, QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, -maxacttmp, *cutrhs);

         if( isintegral )
         {
            /* if the cut is integral, the true coefficient must also be integral;
             * thus we round it to the exact integral value
             */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) > val )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQD(delta, coef, -val);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
                   val, QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp), lb,
                   cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]));

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(!SCIPisPositive(scip, QUAD_TO_DBL(coef)));

            if( SCIPisNegative(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               cutcoefs[cutinds[i]] = QUAD_TO_DBL(coef);
            }
            else
            {
               cutcoefs[cutinds[i]] = 0.0;
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else if( val > 0.0 && SCIPisLE(scip, maxact - val, QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, maxacttmp, -*cutrhs);

         if( isintegral )
         {
            /* if the cut is integral, the true coefficient must also be integral;
             * thus we round it to the exact integral value
             */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) < val )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQD(delta, coef, -val);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
                   val, QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp),
                   cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]), ub);

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(! SCIPisNegative(scip, QUAD_TO_DBL(coef)));

            if( SCIPisPositive(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               cutcoefs[cutinds[i]] = QUAD_TO_DBL(coef);
            }
            else
            {
               cutcoefs[cutinds[i]] = 0.0;
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else /* due to sorting we can stop completely if the precondition was not fulfilled for this variable */
         break;

      ++i;
   }

   return SCIP_OKAY;
}

/** perform activity based coefficient tightening on the given cut; returns TRUE if the cut was detected
 *  to be redundant due to activity bounds
 */
SCIP_Bool SCIPcutsTightenCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< the number of non-zeros in the cut */
   int*                  nchgcoefs           /**< number of changed coefficients */
   )
{
   int i;
   int nintegralvars;
   SCIP_VAR** vars;
   SCIP_Real* absvals;
   SCIP_Real QUAD(maxacttmp);
   SCIP_Real maxact;
   SCIP_Real maxabsval;
   SCIP_Bool redundant;

   assert(nchgcoefs != NULL);

   QUAD_ASSIGN(maxacttmp, 0.0);

   vars = SCIPgetVars(scip);
   nintegralvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   maxabsval = 0.0;
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &absvals, *cutnnz) );

   *nchgcoefs = 0;
   redundant = FALSE;

   for( i = 0; i < *cutnnz; ++i )
   {
      assert(cutinds[i] >= 0);
      assert(vars[cutinds[i]] != NULL);

      if( cutcoefs[i] < 0.0 )
      {
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, -lb) )
            goto TERMINATE;

         if( cutinds[i] < nintegralvars )
         {
            maxabsval = MAX(maxabsval, -cutcoefs[i]);
            absvals[i] = -cutcoefs[i];
         }
         else
         {
            absvals[i] = 0.0;
         }

         SCIPquadprecSumQD(maxacttmp, maxacttmp, lb * cutcoefs[i]);
      }
      else
      {
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, ub) )
            goto TERMINATE;

         if( cutinds[i] < nintegralvars )
         {
            maxabsval = MAX(maxabsval, cutcoefs[i]);
            absvals[i] = cutcoefs[i];
         }
         else
         {
            absvals[i] = 0.0;
         }

         SCIPquadprecSumQD(maxacttmp, maxacttmp, ub * cutcoefs[i]);
      }
   }

   maxact = QUAD_TO_DBL(maxacttmp);

   /* cut is redundant in activity bounds */
   if( SCIPisFeasLE(scip, maxact, *cutrhs) )
   {
      redundant = TRUE;
      goto TERMINATE;
   }

   /* no coefficient tightening can be performed since the precondition doesn't hold for any of the variables */
   if( SCIPisGT(scip, maxact - maxabsval, *cutrhs) )
      goto TERMINATE;

   SCIPsortDownRealRealInt(absvals, cutcoefs, cutinds, *cutnnz);
   SCIPfreeBufferArray(scip, &absvals);

   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *cutnnz;)
   {
      if( cutinds[i] >= nintegralvars )
      {
         ++i;
         continue;
      }

      assert(SCIPvarIsIntegral(vars[cutinds[i]]));

      if( cutcoefs[i] < 0.0 && SCIPisLE(scip, maxact + cutcoefs[i], *cutrhs) )
      {
         SCIP_Real coef = (*cutrhs) - maxact;
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         coef = floor(coef);

         if( coef > cutcoefs[i] )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumDD(delta, coef, -cutcoefs[i]);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQD(tmp, delta, *cutrhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
                   cutcoefs[i], coef, (*cutrhs), QUAD_TO_DBL(tmp), lb,
                   cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]));

            *cutrhs = QUAD_TO_DBL(tmp);

            assert(!SCIPisPositive(scip, coef));

            ++(*nchgcoefs);

            if( SCIPisNegative(scip, coef) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               cutcoefs[i] = coef;
            }
            else
            {
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               cutcoefs[i] = cutcoefs[*cutnnz];
               continue;
            }
         }
      }
      else if( cutcoefs[i] > 0.0 && SCIPisLE(scip, maxact - cutcoefs[i], *cutrhs) )
      {
         SCIP_Real coef = maxact - (*cutrhs);
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         coef = ceil(coef);

         if( coef < cutcoefs[i] )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumDD(delta, coef, -cutcoefs[i]);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQD(tmp, delta, *cutrhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
                   cutcoefs[i], coef, (*cutrhs), QUAD_TO_DBL(tmp),
                   cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]), ub);

            *cutrhs = QUAD_TO_DBL(tmp);

            assert(!SCIPisNegative(scip, coef));

            ++(*nchgcoefs);

            if( SCIPisPositive(scip, coef) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               cutcoefs[i] = coef;
            }
            else
            {
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               cutcoefs[i] = cutcoefs[*cutnnz];
               continue;
            }
         }
      }
      else /* due to sorting we can stop completely if the precondition was not fulfilled for this variable */
         break;

      ++i;
   }

  TERMINATE:
   SCIPfreeBufferArrayNull(scip, &absvals);

   return redundant;
}

/* =========================================== aggregation row =========================================== */


/** create an empty aggregation row */
SCIP_RETCODE SCIPaggrRowCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow             /**< pointer to return aggregation row */
   )
{
   int nvars;
   assert(scip != NULL);
   assert(aggrrow != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, aggrrow) );

   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*aggrrow)->vals, QUAD_ARRAY_SIZE(nvars)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*aggrrow)->inds, nvars) );

   BMSclearMemoryArray((*aggrrow)->vals, QUAD_ARRAY_SIZE(nvars));

   (*aggrrow)->local = FALSE;
   (*aggrrow)->nnz = 0;
   (*aggrrow)->rank = 0;
   QUAD_ASSIGN((*aggrrow)->rhs, 0.0);
   (*aggrrow)->rowsinds = NULL;
   (*aggrrow)->slacksign = NULL;
   (*aggrrow)->rowweights = NULL;
   (*aggrrow)->nrows = 0;
   (*aggrrow)->rowssize = 0;

   return SCIP_OKAY;
}

/** free a aggregation row */
void SCIPaggrRowFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow             /**< pointer to aggregation row that should be freed */
   )
{
   int nvars;
   assert(scip != NULL);
   assert(aggrrow != NULL);

   nvars = SCIPgetNVars(scip);

   SCIPfreeBlockMemoryArray(scip, &(*aggrrow)->inds, nvars);
   SCIPfreeBlockMemoryArray(scip, &(*aggrrow)->vals, QUAD_ARRAY_SIZE(nvars)); /*lint !e647*/
   SCIPfreeBlockMemoryArrayNull(scip, &(*aggrrow)->rowsinds, (*aggrrow)->rowssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*aggrrow)->slacksign, (*aggrrow)->rowssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*aggrrow)->rowweights, (*aggrrow)->rowssize);
   SCIPfreeBlockMemory(scip, aggrrow);
}

/** output aggregation row to file stream */
void SCIPaggrRowPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< pointer to return aggregation row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** vars;
   SCIP_MESSAGEHDLR* messagehdlr;
   int i;

   assert(scip != NULL);
   assert(aggrrow != NULL);

   vars = SCIPgetVars(scip);
   assert(vars != NULL);

   messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr);

   /* print coefficients */
   if( aggrrow->nnz == 0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "0 ");

   for( i = 0; i < aggrrow->nnz; ++i )
   {
      SCIP_Real QUAD(val);

      QUAD_ARRAY_LOAD(val, aggrrow->vals, aggrrow->inds[i]);
      assert(SCIPvarGetProbindex(vars[aggrrow->inds[i]]) == aggrrow->inds[i]);
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.15g<%s> ", QUAD_TO_DBL(val), SCIPvarGetName(vars[aggrrow->inds[i]]));
   }

   /* print right hand side */
   SCIPmessageFPrintInfo(messagehdlr, file, "<= %.15g\n", QUAD_TO_DBL(aggrrow->rhs));
}

/** copy a aggregation row */
SCIP_RETCODE SCIPaggrRowCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow,            /**< pointer to return aggregation row */
   SCIP_AGGRROW*         source              /**< source aggregation row */
   )
{
   int nvars;

   assert(scip != NULL);
   assert(aggrrow != NULL);
   assert(source != NULL);

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBlockMemory(scip, aggrrow) );

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->vals, source->vals, QUAD_ARRAY_SIZE(nvars)) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->inds, source->inds, nvars) );
   (*aggrrow)->nnz = source->nnz;
   QUAD_ASSIGN_Q((*aggrrow)->rhs, source->rhs);

   if( source->nrows > 0 )
   {
      assert(source->rowsinds != NULL);
      assert(source->slacksign != NULL);
      assert(source->rowweights != NULL);

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->rowsinds, source->rowsinds, source->nrows) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->slacksign, source->slacksign, source->nrows) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->rowweights, source->rowweights, source->nrows) );
   }
   else
   {
      (*aggrrow)->rowsinds = NULL;
      (*aggrrow)->slacksign = NULL;
      (*aggrrow)->rowweights = NULL;
   }

   (*aggrrow)->nrows = source->nrows;
   (*aggrrow)->rowssize = source->nrows;
   (*aggrrow)->rank = source->rank;
   (*aggrrow)->local = source->local;

   return SCIP_OKAY;
}

/** add weighted row to aggregation row */
SCIP_RETCODE SCIPaggrRowAddRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_ROW*             row,                /**< row to add to aggregation row */
   SCIP_Real             weight,             /**< scale for adding given row to aggregation row */
   int                   sidetype            /**< specify row side type (-1 = lhs, 0 = automatic, 1 = rhs) */
   )
{
   int i;

   assert(row->lppos >= 0);

   /* update local flag */
   aggrrow->local = aggrrow->local || row->local;

   /* update rank */
   aggrrow->rank = MAX(row->rank, aggrrow->rank);

   {
      SCIP_Real sideval;
      SCIP_Bool uselhs;

      i = aggrrow->nrows++;

      if( aggrrow->nrows > aggrrow->rowssize )
      {
         int newsize = SCIPcalcMemGrowSize(scip, aggrrow->nrows);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowsinds, aggrrow->rowssize, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->slacksign, aggrrow->rowssize, newsize) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowweights, aggrrow->rowssize, newsize) );
         aggrrow->rowssize = newsize;
      }
      aggrrow->rowsinds[i] = SCIProwGetLPPos(row);
      aggrrow->rowweights[i] = weight;

      if ( sidetype == -1 )
      {
         assert( ! SCIPisInfinity(scip, -row->lhs) );
         uselhs = TRUE;
      }
      else if ( sidetype == 1 )
      {
         assert( ! SCIPisInfinity(scip, row->rhs) );
         uselhs = FALSE;
      }
      else
      {
         /* Automatically decide, whether we want to use the left or the right hand side of the row in the summation.
          * If possible, use the side that leads to a positive slack value in the summation.
          */
         if( SCIPisInfinity(scip, row->rhs) || (!SCIPisInfinity(scip, -row->lhs) && weight < 0.0) )
            uselhs = TRUE;
         else
            uselhs = FALSE;
      }

      if( uselhs )
      {
         aggrrow->slacksign[i] = -1;
         sideval = row->lhs - row->constant;
         if( row->integral )
            sideval = SCIPceil(scip, sideval); /* row is integral: round left hand side up */
      }
      else
      {
         aggrrow->slacksign[i] = +1;
         sideval = row->rhs - row->constant;
         if( row->integral )
            sideval = SCIPfloor(scip, sideval); /* row is integral: round right hand side up */
      }

      SCIPquadprecSumQD(aggrrow->rhs, aggrrow->rhs, weight * sideval);
   }

   /* add up coefficients */
   SCIP_CALL( varVecAddScaledRowCoefsQuad(aggrrow->inds, aggrrow->vals, &aggrrow->nnz, row, weight) );

   return SCIP_OKAY;
}

/** Removes a given variable @p var from position @p pos the aggregation row and updates the right-hand side according
 *  to sign of the coefficient, i.e., rhs -= coef * bound, where bound = lb if coef >= 0 and bound = ub, otherwise.
 *
 *  @note: The choice of global or local bounds depend on the validity (global or local) of the aggregation row.
 *
 *  @note: The list of non-zero indices will be updated by swapping the last non-zero index to @p pos.
 */
void SCIPaggrRowCancelVarWithBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_VAR*             var,                /**< variable that should be removed */
   int                   pos,                /**< position of the variable in the aggregation row */
   SCIP_Bool*            valid               /**< pointer to return whether the aggregation row is still valid */
   )
{
   SCIP_Real QUAD(val);
   int v;

   assert(valid != NULL);
   assert(pos >= 0);

   v = aggrrow->inds[pos];
   assert(v == SCIPvarGetProbindex(var));

   QUAD_ARRAY_LOAD(val, aggrrow->vals, v);

   *valid = TRUE;

   /* adjust left and right hand sides with max contribution */
   if( QUAD_TO_DBL(val) < 0.0 )
   {
      SCIP_Real ub = aggrrow->local ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
      if( SCIPisInfinity(scip, ub) )
         QUAD_ASSIGN(aggrrow->rhs, SCIPinfinity(scip));
      else
      {
         SCIPquadprecProdQD(val, val, ub);
         SCIPquadprecSumQQ(aggrrow->rhs, aggrrow->rhs, -val);
      }
   }
   else
   {
      SCIP_Real lb = aggrrow->local ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
      if( SCIPisInfinity(scip, -lb) )
         QUAD_ASSIGN(aggrrow->rhs, SCIPinfinity(scip));
      else
      {
         SCIPquadprecProdQD(val, val, lb);
         SCIPquadprecSumQQ(aggrrow->rhs, aggrrow->rhs, -val);
      }
   }

   QUAD_ASSIGN(val, 0.0);
   QUAD_ARRAY_STORE(aggrrow->vals, v, val);

   /* remove non-zero entry */
   --(aggrrow->nnz);
   aggrrow->inds[pos] = aggrrow->inds[aggrrow->nnz];

   if( SCIPisInfinity(scip, QUAD_HI(aggrrow->rhs)) )
      *valid = FALSE;
}

/** add the objective function with right-hand side @p rhs and scaled by @p scale to the aggregation row */
SCIP_RETCODE SCIPaggrRowAddObjectiveFunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_Real             rhs,                /**< right-hand side of the artificial row */
   SCIP_Real             scale               /**< scalar */
   )
{
   SCIP_VAR** vars;
   SCIP_Real QUAD(val);
   int nvars;

   assert(scip != NULL);
   assert(aggrrow != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* add all variables straight forward if the aggregation row is empty */
   if( aggrrow->nnz == 0 )
   {
      int i;
      for( i = 0; i < nvars; ++i )
      {
         assert(SCIPvarGetProbindex(vars[i]) == i);

         /* skip all variables with zero objective coefficient */
         if( SCIPisZero(scip, scale * SCIPvarGetObj(vars[i])) )
            continue;

         QUAD_ASSIGN(val, scale * SCIPvarGetObj(vars[i]));
         QUAD_ARRAY_STORE(aggrrow->vals, i, val);
         aggrrow->inds[aggrrow->nnz++] = i;
      }

      /* add right-hand side value */
      QUAD_ASSIGN(aggrrow->rhs, scale * rhs);
   }
   else
   {
      int i;
      /* add the non-zeros to the aggregation row and keep non-zero index up to date */
      for( i = 0 ; i < nvars; ++i )
      {
         assert(SCIPvarGetProbindex(vars[i]) == i);

         /* skip all variables with zero objective coefficient */
         if( SCIPisZero(scip, scale * SCIPvarGetObj(vars[i])) )
            continue;

         QUAD_ARRAY_LOAD(val, aggrrow->vals, i); /* val = aggrrow->vals[i] */

         if( QUAD_HI(val) == 0.0 )
            aggrrow->inds[aggrrow->nnz++] = i;

         SCIPquadprecSumQD(val, val, scale * SCIPvarGetObj(vars[i]));

         /* the value must not be exactly zero due to sparsity pattern */
         QUAD_HI(val) = NONZERO(QUAD_HI(val));
         assert(QUAD_HI(val) != 0.0);

         QUAD_ARRAY_STORE(aggrrow->vals, i, val);
      }

      /* add right-hand side value */
      SCIPquadprecSumQD(aggrrow->rhs, aggrrow->rhs, scale * rhs);
   }

   return SCIP_OKAY;
}

/** add weighted constraint to the aggregation row */
SCIP_RETCODE SCIPaggrRowAddCustomCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   int*                  inds,               /**< variable problem indices in constraint to add to the aggregation row */
   SCIP_Real*            vals,               /**< values of constraint to add to the aggregation row */
   int                   len,                /**< length of constraint to add to the aggregation row */
   SCIP_Real             rhs,                /**< right hand side of constraint to add to the aggregation row */
   SCIP_Real             weight,             /**< (positive) scale for adding given constraint to the aggregation row */
   int                   rank,               /**< rank to use for given constraint */
   SCIP_Bool             local               /**< is constraint only valid locally */
   )
{
   int i;

   assert(weight >= 0.0);
   assert(!SCIPisInfinity(scip, REALABS(weight * rhs)));

   /* update local flag */
   aggrrow->local = aggrrow->local || local;

   /* update rank */
   aggrrow->rank = MAX(rank, aggrrow->rank);

   /* add right hand side value */
   SCIPquadprecSumQD(aggrrow->rhs, aggrrow->rhs, weight * rhs);

   /* add the non-zeros to the aggregation row and keep non-zero index up to date */
   for( i = 0 ; i < len; ++i )
   {
      SCIP_Real QUAD(val);
      int probindex = inds[i];

      QUAD_ARRAY_LOAD(val, aggrrow->vals, probindex); /* val = aggrrow->vals[probindex] */

      if( QUAD_HI(val) == 0.0 )
         aggrrow->inds[aggrrow->nnz++] = probindex;

      SCIPquadprecSumQD(val, val, vals[i] * weight);

      /* the value must not be exactly zero due to sparsity pattern */
      QUAD_HI(val) = NONZERO(QUAD_HI(val));
      assert(QUAD_HI(val) != 0.0);

      QUAD_ARRAY_STORE(aggrrow->vals, probindex, val);
   }

   return SCIP_OKAY;
}

/** clear all entries int the aggregation row but don't free memory */
void SCIPaggrRowClear(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   int i;
   SCIP_Real QUAD(tmp);

   QUAD_ASSIGN(tmp, 0.0);

   for( i = 0; i < aggrrow->nnz; ++i )
   {
      QUAD_ARRAY_STORE(aggrrow->vals, aggrrow->inds[i], tmp);
   }

   aggrrow->nnz = 0;
   aggrrow->nrows = 0;
   aggrrow->rank = 0;
   QUAD_ASSIGN(aggrrow->rhs, 0.0);
   aggrrow->local = FALSE;
}

/** calculates the efficacy norm of the given aggregation row, which depends on the "separating/efficacynorm" parameter
 *
 *  @return the efficacy norm of the given aggregation row, which depends on the "separating/efficacynorm" parameter
 */
SCIP_Real SCIPaggrRowCalcEfficacyNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   return calcEfficacyNormQuad(scip, aggrrow->vals, aggrrow->inds, aggrrow->nnz);
}

/** Adds one row to the aggregation row. Differs from SCIPaggrRowAddRow() by providing some additional
 *  parameters required for SCIPaggrRowSumRows()
 */
static
SCIP_RETCODE addOneRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_ROW*             row,                /**< the row to add */
   SCIP_Real             weight,             /**< weight of row to add */
   SCIP_Bool             sidetypebasis,      /**< choose sidetypes of row (lhs/rhs) based on basis information? */
   SCIP_Bool             allowlocal,         /**< should local rows allowed to be used? */
   int                   negslack,           /**< should negative slack variables allowed to be used? (0: no, 1: only for integral rows, 2: yes) */
   int                   maxaggrlen,         /**< maximal length of aggregation row */
   SCIP_Bool*            rowtoolong          /**< is the aggregated row too long */
   )
{
   SCIP_Real sideval;
   SCIP_Bool uselhs;
   int i;

   *rowtoolong = FALSE;

   if( SCIPisFeasZero(scip, weight) || SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !allowlocal) )
   {
      return SCIP_OKAY;
   }

   if( sidetypebasis && !SCIPisEQ(scip, SCIProwGetLhs(row), SCIProwGetRhs(row)) )
   {
      SCIP_BASESTAT stat = SCIProwGetBasisStatus(row);

      if( stat == SCIP_BASESTAT_LOWER )
      {
         assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );
         uselhs = TRUE;
      }
      else if( stat == SCIP_BASESTAT_UPPER )
      {
         assert( ! SCIPisInfinity(scip, SCIProwGetRhs(row)) );
         uselhs = FALSE;
      }
      else if( SCIPisInfinity(scip, SCIProwGetRhs(row)) ||
         (weight < 0.0 && ! SCIPisInfinity(scip, -SCIProwGetLhs(row))) )
      {
         uselhs = TRUE;
      }
      else
      {
         uselhs = FALSE;
      }
   }
   else if( (weight < 0.0 && !SCIPisInfinity(scip, -row->lhs)) || SCIPisInfinity(scip, row->rhs) )
   {
      uselhs = TRUE;
   }
   else
   {
      uselhs = FALSE;
   }

   if( uselhs )
   {
      assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );

      if( weight > 0.0 && ((negslack == 0) || (negslack == 1 && !row->integral)) )
         return SCIP_OKAY;

      sideval = row->lhs - row->constant;
      /* row is integral? round left hand side up */
      if( row->integral )
         sideval = SCIPceil(scip, sideval);
   }
   else
   {
      assert( ! SCIPisInfinity(scip, SCIProwGetRhs(row)) );

      if( weight < 0.0 && ((negslack == 0) || (negslack == 1 && !row->integral)) )
         return SCIP_OKAY;

      sideval = row->rhs - row->constant;
      /* row is integral? round right hand side down */
      if( row->integral )
         sideval = SCIPfloor(scip, sideval);
   }

   /* add right hand side, update rank and local flag */
   SCIPquadprecSumQD(aggrrow->rhs, aggrrow->rhs, sideval * weight);
   aggrrow->rank = MAX(aggrrow->rank, row->rank);
   aggrrow->local = aggrrow->local || row->local;

   /* ensure the array for storing the row information is large enough */
   i = aggrrow->nrows++;
   if( aggrrow->nrows > aggrrow->rowssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, aggrrow->nrows);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowsinds, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->slacksign, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowweights, aggrrow->rowssize, newsize) );
      aggrrow->rowssize = newsize;
   }

   /* add information of addditional row */
   aggrrow->rowsinds[i] = row->lppos;
   aggrrow->rowweights[i] = weight;
   aggrrow->slacksign[i] = uselhs ? -1 : 1;

   /* add up coefficients */
   SCIP_CALL( varVecAddScaledRowCoefsQuad(aggrrow->inds, aggrrow->vals, &aggrrow->nnz, row, weight) );

   /* check if row is too long now */
   if( aggrrow->nnz > maxaggrlen )
      *rowtoolong = TRUE;

   return SCIP_OKAY;
}

/** aggregate rows using the given weights; the current content of the aggregation
 *  row, \p aggrrow, gets overwritten
 */
SCIP_RETCODE SCIPaggrRowSumRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  rowinds,            /**< array to store indices of non-zero entries of the weights array, or NULL */
   int                   nrowinds,           /**< number of non-zero entries in weights array, -1 if rowinds is NULL */
   SCIP_Bool             sidetypebasis,      /**< choose sidetypes of row (lhs/rhs) based on basis information? */
   SCIP_Bool             allowlocal,         /**< should local rows allowed to be used? */
   int                   negslack,           /**< should negative slack variables allowed to be used? (0: no, 1: only for integral rows, 2: yes) */
   int                   maxaggrlen,         /**< maximal length of aggregation row */
   SCIP_Bool*            valid               /**< is the aggregation valid */
   )
{
   SCIP_ROW** rows;
   SCIP_VAR** vars;
   int nrows;
   int nvars;
   int k;
   SCIP_Bool rowtoolong;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   SCIPaggrRowClear(aggrrow);
   *valid = FALSE;

   if( rowinds != NULL && nrowinds > -1 )
   {
      for( k = 0; k < nrowinds; ++k )
      {
         SCIP_CALL( addOneRow(scip, aggrrow, rows[rowinds[k]], weights[rowinds[k]], sidetypebasis, allowlocal,
                              negslack, maxaggrlen, &rowtoolong) );

         if( rowtoolong )
            return SCIP_OKAY;
      }
   }
   else
   {
      for( k = 0; k < nrows; ++k )
      {
         SCIP_CALL( addOneRow(scip, aggrrow, rows[k], weights[k], sidetypebasis, allowlocal,
                              negslack, maxaggrlen, &rowtoolong) );

         if( rowtoolong )
            return SCIP_OKAY;
      }
   }

   SCIPaggrRowRemoveZeros(scip, aggrrow, valid);

   return SCIP_OKAY;
}

/** checks for cut redundancy and performs activity based coefficient tightening;
 *  removes coefficients that are zero with QUAD_EPSILON tolerance and uses variable bounds
 *  to remove small coefficients (relative to the maximum absolute coefficient)
 */
static
SCIP_RETCODE postprocessCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut a local cut */
   int*                  cutinds,            /**< variable problem indices of non-zeros in cut */
   SCIP_Real*            cutcoefs,           /**< non-zeros coefficients of cut */
   int*                  nnz,                /**< number non-zeros coefficients of cut */
   SCIP_Real*            cutrhs,             /**< right hand side of cut */
   SCIP_Bool*            success             /**< pointer to return whether post-processing was succesful or cut is redundant */
   )
{
   int i;
   SCIP_Bool redundant;
   SCIP_Real maxcoef;
   SCIP_Real minallowedcoef;
   SCIP_Real QUAD(rhs);

   assert(scip != NULL);
   assert(cutinds != NULL);
   assert(cutcoefs != NULL);
   assert(cutrhs != NULL);

   *success = FALSE;

   QUAD_ASSIGN(rhs, *cutrhs);

   if( removeZeros(scip, SCIPfeastol(scip), cutislocal, cutcoefs, QUAD(&rhs), cutinds, nnz) )
   {
      /* right hand side was changed to infinity -> cut is redundant */
      return SCIP_OKAY;
   }

   SCIP_CALL( cutTightenCoefs(scip, cutislocal, cutcoefs, QUAD(&rhs), cutinds, nnz, &redundant) );

   if( redundant )
   {
      /* cut is redundant */
      return SCIP_OKAY;
   }

   maxcoef = 0.0;
   for( i = 0; i < *nnz; ++i )
   {
      SCIP_Real absval = REALABS(cutcoefs[cutinds[i]]);
      maxcoef = MAX(absval, maxcoef);
   }

   maxcoef /= scip->set->sepa_maxcoefratio;
   minallowedcoef = SCIPsumepsilon(scip);
   minallowedcoef = MAX(minallowedcoef, maxcoef);

   *success = ! removeZeros(scip, minallowedcoef, cutislocal, cutcoefs, QUAD(&rhs), cutinds, nnz);
   *cutrhs = QUAD_TO_DBL(rhs);

   return SCIP_OKAY;
}


/** checks for cut redundancy and performs activity based coefficient tightening;
 *  removes coefficients that are zero with QUAD_EPSILON tolerance and uses variable bounds
 *  to remove small coefficients (relative to the maximum absolute coefficient).
 *  The cutcoefs must be a quad precision array, i.e. allocated with size
 *  QUAD_ARRAY_SIZE(nvars) and accessed with QUAD_ARRAY_LOAD and QUAD_ARRAY_STORE
 *  macros.
 */
static
SCIP_RETCODE postprocessCutQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut a local cut */
   int*                  cutinds,            /**< variable problem indices of non-zeros in cut */
   SCIP_Real*            cutcoefs,           /**< non-zeros coefficients of cut */
   int*                  nnz,                /**< number non-zeros coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< right hand side of cut */
   SCIP_Bool*            success             /**< pointer to return whether the cleanup was successful or if it is useless */
   )
{
   int i;
   SCIP_Bool redundant;
   SCIP_Real maxcoef;
   SCIP_Real minallowedcoef;

   assert(scip != NULL);
   assert(cutinds != NULL);
   assert(cutcoefs != NULL);
   assert(QUAD_HI(cutrhs) != NULL);

   *success = FALSE;

   if( removeZerosQuad(scip, SCIPfeastol(scip), cutislocal, cutcoefs, QUAD(cutrhs), cutinds, nnz) )
   {
      /* right hand side was changed to infinity -> cut is redundant */
      return SCIP_OKAY;
   }

   SCIP_CALL( cutTightenCoefsQuad(scip, cutislocal, cutcoefs, QUAD(cutrhs), cutinds, nnz, &redundant) );
   if( redundant )
   {
      /* cut is redundant */
      return SCIP_OKAY;
   }

   maxcoef = 0.0;
   for( i = 0; i < *nnz; ++i )
   {
      SCIP_Real abscoef;
      SCIP_Real QUAD(coef);
      QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]); /* coef = cutcoefs[cutinds[i]] */
      abscoef = REALABS(QUAD_TO_DBL(coef));
      maxcoef = MAX(abscoef, maxcoef);
   }

   maxcoef /= scip->set->sepa_maxcoefratio;
   minallowedcoef = SCIPsumepsilon(scip);
   minallowedcoef = MAX(minallowedcoef, maxcoef);

   *success = ! removeZerosQuad(scip, minallowedcoef, cutislocal, cutcoefs, QUAD(cutrhs), cutinds, nnz);

   return SCIP_OKAY;
}

/** removes almost zero entries from the aggregation row. */
void SCIPaggrRowRemoveZeros(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_Bool*            valid               /**< pointer to return whether the aggregation row is still valid */
   )
{
   assert(aggrrow != NULL);

   *valid = ! removeZerosQuad(scip, SCIPsumepsilon(scip), aggrrow->local, aggrrow->vals, QUAD(&aggrrow->rhs), aggrrow->inds, &aggrrow->nnz);
}

/** get number of aggregated rows */
int SCIPaggrRowGetNRows(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->nrows;
}

/** get array with lp positions of rows used in aggregation */
int* SCIPaggrRowGetRowInds(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   assert(aggrrow != NULL);
   assert(aggrrow->rowsinds != NULL || aggrrow->nrows == 0);

   return aggrrow->rowsinds;
}

/** get array with weights of aggregated rows */
SCIP_Real* SCIPaggrRowGetRowWeights(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   assert(aggrrow != NULL);
   assert(aggrrow->rowweights != NULL || aggrrow->nrows == 0);

   return aggrrow->rowweights;
}

/** checks whether a given row has been added to the aggregation row */
SCIP_Bool SCIPaggrRowHasRowBeenAdded(
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_ROW*             row                 /**< row for which it is checked whether it has been added to the aggregation */
   )
{
   int i;
   int rowind;

   assert(aggrrow != NULL);
   assert(row != NULL);

   rowind = SCIProwGetLPPos(row);

   for( i = 0; i < aggrrow->nrows; ++i )
      if( aggrrow->rowsinds[i] == rowind )
         return TRUE;

   return FALSE;
}

/** gets the array of corresponding variable problem indices for each non-zero in the aggregation row */
int* SCIPaggrRowGetInds(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->inds;
}

/** gets the number of non-zeros in the aggregation row */
int SCIPaggrRowGetNNz(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->nnz;
}

/** gets the rank of the aggregation row */
int SCIPaggrRowGetRank(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->rank;
}

/** checks if the aggregation row is only valid locally */
SCIP_Bool SCIPaggrRowIsLocal(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->local;
}

/** gets the right hand side of the aggregation row */
SCIP_Real SCIPaggrRowGetRhs(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return QUAD_TO_DBL(aggrrow->rhs);
}

/* =========================================== c-MIR =========================================== */

#define MAXCMIRSCALE               1e+6 /**< maximal scaling (scale/(1-f0)) allowed in c-MIR calculations */

/** finds the best lower bound of the variable to use for MIR transformation */
static
SCIP_RETCODE findBestLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestlb,             /**< pointer to store best bound value */
   int*                  bestlbtype          /**< pointer to store best bound type */
   )
{
   assert(bestlb != NULL);
   assert(bestlbtype != NULL);

   *bestlb = SCIPvarGetLbGlobal(var);
   *bestlbtype = -1;

   if( allowlocal )
   {
      SCIP_Real loclb;

      loclb = SCIPvarGetLbLocal(var);
      if( SCIPisGT(scip, loclb, *bestlb) )
      {
         *bestlb = loclb;
         *bestlbtype = -2;
      }
   }

   if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      SCIP_Real bestvlb;
      int bestvlbidx;

      SCIP_CALL( SCIPgetVarClosestVlb(scip, var, sol, &bestvlb, &bestvlbidx) );
      if( bestvlbidx >= 0
         && (bestvlb > *bestlb || (*bestlbtype < 0 && SCIPisGE(scip, bestvlb, *bestlb))) )
      {
         SCIP_VAR** vlbvars;

         /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
         /**@todo this check is not needed for continuous variables; but allowing all but binary variables
          *       to be replaced by variable bounds seems to be buggy (wrong result on gesa2)
          */
         vlbvars = SCIPvarGetVlbVars(var);
         assert(vlbvars != NULL);
         if( SCIPvarGetProbindex(vlbvars[bestvlbidx]) < SCIPvarGetProbindex(var) )
         {
            *bestlb = bestvlb;
            *bestlbtype = bestvlbidx;
         }
      }
   }

   return SCIP_OKAY;
}

/** finds the best upper bound of the variable to use for MIR transformation */
static
SCIP_RETCODE findBestUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestub,             /**< pointer to store best bound value */
   int*                  bestubtype          /**< pointer to store best bound type */
   )
{
   assert(bestub != NULL);
   assert(bestubtype != NULL);

   *bestub = SCIPvarGetUbGlobal(var);
   *bestubtype = -1;

   if( allowlocal )
   {
      SCIP_Real locub;

      locub = SCIPvarGetUbLocal(var);
      if( SCIPisLT(scip, locub, *bestub) )
      {
         *bestub = locub;
         *bestubtype = -2;
      }
   }

   if( usevbds && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      SCIP_Real bestvub;
      int bestvubidx;

      SCIP_CALL( SCIPgetVarClosestVub(scip, var, sol, &bestvub, &bestvubidx) );
      if( bestvubidx >= 0
         && (bestvub < *bestub || (*bestubtype < 0 && SCIPisLE(scip, bestvub, *bestub))) )
      {
         SCIP_VAR** vubvars;

         /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
         /**@todo this check is not needed for continuous variables; but allowing all but binary variables
          *       to be replaced by variable bounds seems to be buggy (wrong result on gesa2)
          */
         vubvars = SCIPvarGetVubVars(var);
         assert(vubvars != NULL);
         if( SCIPvarGetProbindex(vubvars[bestvubidx]) < SCIPvarGetProbindex(var) )
         {
            *bestub = bestvub;
            *bestubtype = bestvubidx;
         }
      }
   }

   return SCIP_OKAY;
}

/** determine the best bounds with respect to the given solution for complementing the given variable */
static
SCIP_RETCODE determineBestBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to determine best bound for */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   SCIP_Bool             ignoresol,          /**< should the LP solution be ignored? (eg, apply MIR to dualray) */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real*            bestlb,             /**< pointer to store best lower bound of variable */
   SCIP_Real*            bestub,             /**< pointer to store best upper bound of variable */
   int*                  bestlbtype,         /**< pointer to store type of best lower bound of variable */
   int*                  bestubtype,         /**< pointer to store type of best upper bound of variable */
   SCIP_BOUNDTYPE*       selectedbound,      /**< pointer to store whether the lower bound or the upper bound should be preferred */
   SCIP_Bool*            freevariable        /**< pointer to store if this is a free variable */
   )
{
   int v;

   v = SCIPvarGetProbindex(var);

   /* check if the user specified a bound to be used */
   if( boundsfortrans != NULL && boundsfortrans[v] > -3 )
   {
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || ( boundsfortrans[v] == -2 || boundsfortrans[v] == -1 ));

      /* user has explicitly specified a bound to be used */
      if( boundtypesfortrans[v] == SCIP_BOUNDTYPE_LOWER )
      {
         /* user wants to use lower bound */
         *bestlbtype = boundsfortrans[v];
         if( *bestlbtype == -1 )
            *bestlb = SCIPvarGetLbGlobal(var); /* use global standard lower bound */
         else if( *bestlbtype == -2 )
            *bestlb = SCIPvarGetLbLocal(var);  /* use local standard lower bound */
         else
         {
            SCIP_VAR** vlbvars;
            SCIP_Real* vlbcoefs;
            SCIP_Real* vlbconsts;
            int k;

            assert(!ignoresol);

            /* use the given variable lower bound */
            vlbvars = SCIPvarGetVlbVars(var);
            vlbcoefs = SCIPvarGetVlbCoefs(var);
            vlbconsts = SCIPvarGetVlbConstants(var);
            k = boundsfortrans[v];
            assert(k >= 0 && k < SCIPvarGetNVlbs(var));
            assert(vlbvars != NULL);
            assert(vlbcoefs != NULL);
            assert(vlbconsts != NULL);

            *bestlb = vlbcoefs[k] * (sol == NULL ? SCIPvarGetLPSol(vlbvars[k]) : SCIPgetSolVal(scip, sol, vlbvars[k])) + vlbconsts[k];
         }

         assert(!SCIPisInfinity(scip, - *bestlb));
         *selectedbound = SCIP_BOUNDTYPE_LOWER;

         /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
         SCIP_CALL( findBestUb(scip, var, sol, usevbds && fixintegralrhs, allowlocal && fixintegralrhs, bestub, bestubtype) );
      }
      else
      {
         assert(boundtypesfortrans[v] == SCIP_BOUNDTYPE_UPPER);

         /* user wants to use upper bound */
         *bestubtype = boundsfortrans[v];
         if( *bestubtype == -1 )
            *bestub = SCIPvarGetUbGlobal(var); /* use global standard upper bound */
         else if( *bestubtype == -2 )
            *bestub = SCIPvarGetUbLocal(var);  /* use local standard upper bound */
         else
         {
            SCIP_VAR** vubvars;
            SCIP_Real* vubcoefs;
            SCIP_Real* vubconsts;
            int k;

            assert(!ignoresol);

            /* use the given variable upper bound */
            vubvars = SCIPvarGetVubVars(var);
            vubcoefs = SCIPvarGetVubCoefs(var);
            vubconsts = SCIPvarGetVubConstants(var);
            k = boundsfortrans[v];
            assert(k >= 0 && k < SCIPvarGetNVubs(var));
            assert(vubvars != NULL);
            assert(vubcoefs != NULL);
            assert(vubconsts != NULL);

            /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
            *bestub = vubcoefs[k] * (sol == NULL ? SCIPvarGetLPSol(vubvars[k]) : SCIPgetSolVal(scip, sol, vubvars[k])) + vubconsts[k];
         }

         assert(!SCIPisInfinity(scip, *bestub));
         *selectedbound = SCIP_BOUNDTYPE_UPPER;

         /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
         SCIP_CALL( findBestLb(scip, var, sol, usevbds && fixintegralrhs, allowlocal && fixintegralrhs, bestlb, bestlbtype) );
      }
   }
   else
   {
      SCIP_Real varsol;

      /* bound selection should be done automatically */

      /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
      SCIP_CALL( findBestLb(scip, var, sol, usevbds, allowlocal, bestlb, bestlbtype) );

      /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
      SCIP_CALL( findBestUb(scip, var, sol, usevbds, allowlocal, bestub, bestubtype) );

      /* check, if variable is free variable */
      if( SCIPisInfinity(scip, - *bestlb) && SCIPisInfinity(scip, *bestub) )
      {
         /* we found a free variable in the row with non-zero coefficient
            *  -> MIR row can't be transformed in standard form
            */
         *freevariable = TRUE;
         return SCIP_OKAY;
      }

      if( !ignoresol )
      {
         /* select transformation bound */
         varsol = (sol == NULL ? SCIPvarGetLPSol(var) : SCIPgetSolVal(scip, sol, var));

         if( SCIPisInfinity(scip, *bestub) ) /* if there is no ub, use lb */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( SCIPisInfinity(scip, - *bestlb) ) /* if there is no lb, use ub */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( SCIPisLT(scip, varsol, (1.0 - boundswitch) * (*bestlb) + boundswitch * (*bestub)) )
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( SCIPisGT(scip, varsol, (1.0 - boundswitch) * (*bestlb) + boundswitch * (*bestub)) )
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( *bestlbtype == -1 )  /* prefer global standard bounds */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( *bestubtype == -1 )  /* prefer global standard bounds */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( *bestlbtype >= 0 )   /* prefer variable bounds over local bounds */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( *bestubtype >= 0 )   /* prefer variable bounds over local bounds */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else                         /* no decision yet? just use lower bound */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
      }
      else
      {
         SCIP_Real glbub = SCIPvarGetUbGlobal(var);
         SCIP_Real glblb = SCIPvarGetLbGlobal(var);
         SCIP_Real distlb = REALABS(glblb - *bestlb);
         SCIP_Real distub = REALABS(glbub - *bestub);

         assert(!SCIPisInfinity(scip, - *bestlb) || !SCIPisInfinity(scip, *bestub));

         if( SCIPisInfinity(scip, - *bestlb) )
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( !SCIPisNegative(scip, *bestlb) )
         {
            if( SCIPisInfinity(scip, *bestub) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else if( SCIPisZero(scip, glblb) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else if( SCIPisLE(scip, distlb, distub) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else
               *selectedbound = SCIP_BOUNDTYPE_UPPER;
         }
         else
         {
            assert(!SCIPisInfinity(scip, - *bestlb));
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         }
      }
   }

   return SCIP_OKAY;
}

/** Transform equation \f$ a \cdot x = b; lb \leq x \leq ub \f$ into standard form
 *    \f$ a^\prime \cdot x^\prime = b,\; 0 \leq x^\prime \leq ub' \f$.
 *
 *  Transform variables (lb or ub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \mbox{if lb is used in transformation}\\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if ub is used in transformation}
 *  \end{array}
 *  \f]
 *  and move the constant terms \f$ a_j\, lb_j \f$ or \f$ a_j\, ub_j \f$ to the rhs.
 *
 *  Transform variables (vlb or vub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - (bl_j\, zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \mbox{if vlb is used in transf.} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if vub is used in transf.}
 *  \end{array}
 *  \f]
 *  move the constant terms \f$ a_j\, dl_j \f$ or \f$ a_j\, du_j \f$ to the rhs, and update the coefficient of the VLB variable:
 *  \f[
 *  \begin{array}{ll}
 *    a_{zl_j} := a_{zl_j} + a_j\, bl_j,& \mbox{or} \\
 *    a_{zu_j} := a_{zu_j} + a_j\, bu_j &
 *  \end{array}
 *  \f]
 */
static
SCIP_RETCODE cutsTransformMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   SCIP_Bool             ignoresol,          /**< should the LP solution be ignored? (eg, apply MIR to dualray) */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in MIR row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_Real QUAD(tmp);
   SCIP_Real* bestlbs;
   SCIP_Real* bestubs;
   int* bestlbtypes;
   int* bestubtypes;
   SCIP_BOUNDTYPE* selectedbounds;
   int i;
   int aggrrowintstart;
   int nvars;
   int firstcontvar;
   SCIP_VAR** vars;

   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);

   *freevariable = FALSE;
   *localbdsused = FALSE;

   /* allocate temporary memory to store best bounds and bound types */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbs, 2*(*nnz)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubs, 2*(*nnz)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbtypes, 2*(*nnz)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubtypes, 2*(*nnz)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selectedbounds, 2*(*nnz)) );

   /* start with continuous variables, because using variable bounds can affect the untransformed integral
    * variables, and these changes have to be incorporated in the transformation of the integral variables
    * (continuous variables have largest problem indices!)
    */
   SCIPsortDownInt(cutinds, *nnz);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   firstcontvar = nvars - SCIPgetNContVars(scip);

   /* determine the best bounds for the continous variables */
   for( i = 0; i < *nnz && cutinds[i] >= firstcontvar; ++i )
   {
      SCIP_CALL( determineBestBounds(scip, vars[cutinds[i]], sol, boundswitch, usevbds, allowlocal, fixintegralrhs,
                                     ignoresol, boundsfortrans, boundtypesfortrans,
                                     bestlbs + i, bestubs + i, bestlbtypes + i, bestubtypes + i, selectedbounds + i, freevariable) );

      if( *freevariable )
         goto TERMINATE;
   }

   /* remember start of integer variables in the aggrrow */
   aggrrowintstart = i;

   /* perform bound substitution for continuous variables */
   for( i = 0; i < aggrrowintstart; ++i )
   {
      SCIP_Real QUAD(coef);
      int v = cutinds[i];
      SCIP_VAR* var = vars[v];

      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      if( selectedbounds[i] == SCIP_BOUNDTYPE_LOWER )
      {
         assert(!SCIPisInfinity(scip, -bestlbs[i]));

         /* use lower bound as transformation bound: x'_j := x_j - lb_j */
         boundtype[i] = bestlbtypes[i];
         varsign[i] = +1;

         /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
         if( bestlbtypes[i] < 0 )
         {
            SCIPquadprecProdQD(tmp, coef, bestlbs[i]);
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
            *localbdsused = *localbdsused || (bestlbtypes[i] == -2);
         }
         else
         {
            SCIP_VAR** vlbvars;
            SCIP_Real* vlbcoefs;
            SCIP_Real* vlbconsts;
            SCIP_Real QUAD(zcoef);
            int zidx;

            vlbvars = SCIPvarGetVlbVars(var);
            vlbcoefs = SCIPvarGetVlbCoefs(var);
            vlbconsts = SCIPvarGetVlbConstants(var);
            assert(vlbvars != NULL);
            assert(vlbcoefs != NULL);
            assert(vlbconsts != NULL);

            assert(0 <= bestlbtypes[i] && bestlbtypes[i] < SCIPvarGetNVlbs(var));
            assert(SCIPvarIsActive(vlbvars[bestlbtypes[i]]));
            zidx = SCIPvarGetProbindex(vlbvars[bestlbtypes[i]]);
            assert(0 <= zidx && zidx < firstcontvar);

            SCIPquadprecProdQD(tmp, coef, vlbconsts[bestlbtypes[i]]);
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);

            /* check if integral variable already exists in the row and update sparsity pattern */
            QUAD_ARRAY_LOAD(zcoef, cutcoefs, zidx);
            if( QUAD_HI(zcoef) == 0.0 )
               cutinds[(*nnz)++] = zidx;

            SCIPquadprecProdQD(coef, coef, vlbcoefs[bestlbtypes[i]]);
            SCIPquadprecSumQQ(zcoef, zcoef, coef);
            QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
            QUAD_ARRAY_STORE(cutcoefs, zidx, zcoef);
            assert(QUAD_HI(zcoef) != 0.0);
         }
      }
      else
      {
         assert(!SCIPisInfinity(scip, bestubs[i]));

         /* use upper bound as transformation bound: x'_j := ub_j - x_j */
         boundtype[i] = bestubtypes[i];
         varsign[i] = -1;

         /* standard (bestubtype < 0) or variable (bestubtype >= 0) upper bound? */
         if( bestubtypes[i] < 0 )
         {
            SCIPquadprecProdQD(tmp, coef, bestubs[i]);
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
            *localbdsused = *localbdsused || (bestubtypes[i] == -2);
         }
         else
         {
            SCIP_VAR** vubvars;
            SCIP_Real* vubcoefs;
            SCIP_Real* vubconsts;
            SCIP_Real QUAD(zcoef);
            int zidx;

            vubvars = SCIPvarGetVubVars(var);
            vubcoefs = SCIPvarGetVubCoefs(var);
            vubconsts = SCIPvarGetVubConstants(var);
            assert(vubvars != NULL);
            assert(vubcoefs != NULL);
            assert(vubconsts != NULL);

            assert(0 <= bestubtypes[i] && bestubtypes[i] < SCIPvarGetNVubs(var));
            assert(SCIPvarIsActive(vubvars[bestubtypes[i]]));
            zidx = SCIPvarGetProbindex(vubvars[bestubtypes[i]]);
            assert(zidx >= 0);

            SCIPquadprecProdQD(tmp, coef, vubconsts[bestubtypes[i]]);
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);

            /* check if integral variable already exists in the row and update sparsity pattern */
            QUAD_ARRAY_LOAD(zcoef, cutcoefs, zidx);
            if( QUAD_HI(zcoef) == 0.0 )
               cutinds[(*nnz)++] = zidx;

            SCIPquadprecProdQD(coef, coef, vubcoefs[bestubtypes[i]]);
            SCIPquadprecSumQQ(zcoef, zcoef, coef);
            QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
            QUAD_ARRAY_STORE(cutcoefs, zidx, zcoef);
            assert(QUAD_HI(zcoef) != 0.0);
         }
      }
   }

   /* remove integral variables that now have a zero coefficient due to variable bound usage of continuous variables
    * and determine the bound to use for the integer variables that are left
    */
   while( i < *nnz )
   {
      SCIP_Real QUAD(coef);
      int v = cutinds[i];
      assert(cutinds[i] < firstcontvar);

      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      /* due to variable bound usage for the continous variables cancellation may have occurred */
      if( EPSZ(QUAD_TO_DBL(coef), QUAD_EPSILON) )
      {
         QUAD_ASSIGN(coef, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, coef);
         --(*nnz);
         cutinds[i] = cutinds[*nnz];
         /* do not increase i, since last element is copied to the i-th position */
         continue;
      }

      /* determine the best bounds for the integral variable, usevbd can be set to FALSE here as vbds are only used for continous variables */
      SCIP_CALL( determineBestBounds(scip, vars[v], sol, boundswitch, FALSE, allowlocal, fixintegralrhs,
                                     ignoresol, boundsfortrans, boundtypesfortrans,
                                     bestlbs + i, bestubs + i, bestlbtypes + i, bestubtypes + i, selectedbounds + i, freevariable) );

      /* increase i */
      ++i;

      if( *freevariable )
         goto TERMINATE;
   }

   /* now perform the bound substitution on the remaining integral variables which only uses standard bounds */
   for( i = aggrrowintstart; i < *nnz; ++i )
   {
      SCIP_Real QUAD(coef);
      int v = cutinds[i];

      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      /* perform bound substitution */
      if( selectedbounds[i] == SCIP_BOUNDTYPE_LOWER )
      {
         assert(!SCIPisInfinity(scip, - bestlbs[i]));
         assert(bestlbtypes[i] < 0);

         /* use lower bound as transformation bound: x'_j := x_j - lb_j */
         boundtype[i] = bestlbtypes[i];
         varsign[i] = +1;

         /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
         SCIPquadprecProdQD(tmp, coef, bestlbs[i]);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
         *localbdsused = *localbdsused || (bestlbtypes[i] == -2);
      }
      else
      {
         assert(!SCIPisInfinity(scip, bestubs[i]));
         assert(bestubtypes[i] < 0);

         /* use upper bound as transformation bound: x'_j := ub_j - x_j */
         boundtype[i] = bestubtypes[i];
         varsign[i] = -1;

         /* standard (bestubtype < 0) or variable (bestubtype >= 0) upper bound? */
         SCIPquadprecProdQD(tmp, coef, bestubs[i]);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
         *localbdsused = *localbdsused || (bestubtypes[i] == -2);
      }
   }

   if( fixintegralrhs )
   {
      SCIP_Real f0;

      /* check if rhs is fractional */
      f0 = EPSFRAC(QUAD_TO_DBL(*cutrhs), SCIPsumepsilon(scip));
      if( f0 < minfrac || f0 > maxfrac )
      {
         SCIP_Real bestviolgain;
         SCIP_Real bestnewf0;
         int besti;

         /* choose complementation of one variable differently such that f0 is in correct range */
         besti = -1;
         bestviolgain = -1e+100;
         bestnewf0 = 1.0;
         for( i = 0; i < *nnz; i++ )
         {
            int v;
            SCIP_Real QUAD(coef);

            v = cutinds[i];
            assert(0 <= v && v < nvars);

            QUAD_ARRAY_LOAD(coef, cutcoefs, v);
            assert(!EPSZ(QUAD_TO_DBL(coef), QUAD_EPSILON));

            if( boundtype[i] < 0
               && ((varsign[i] == +1 && !SCIPisInfinity(scip, bestubs[i]) && bestubtypes[i] < 0)
                  || (varsign[i] == -1 && !SCIPisInfinity(scip, -bestlbs[i]) && bestlbtypes[i] < 0)) )
            {
               SCIP_Real fj;
               SCIP_Real newfj;
               SCIP_Real newrhs;
               SCIP_Real newf0;
               SCIP_Real solval;
               SCIP_Real viol;
               SCIP_Real newviol;
               SCIP_Real violgain;

               /* currently:              a'_j =  varsign * a_j  ->  f'_j =  a'_j - floor(a'_j)
                * after complementation: a''_j = -varsign * a_j  -> f''_j = a''_j - floor(a''_j) = 1 - f'_j
                *                        rhs'' = rhs' + varsign * a_j * (lb_j - ub_j)
                * cut violation from f0 and fj:   f'_0 -  f'_j *  x'_j
                * after complementation:         f''_0 - f''_j * x''_j
                *
                * for continuous variables, we just set f'_j = f''_j = |a'_j|
                */
               newrhs = QUAD_TO_DBL(*cutrhs) + varsign[i] * QUAD_TO_DBL(coef) * (bestlbs[i] - bestubs[i]);
               newf0 = EPSFRAC(newrhs, SCIPsumepsilon(scip));
               if( newf0 < minfrac || newf0 > maxfrac )
                  continue;
               if( v >= firstcontvar )
               {
                  fj = REALABS(QUAD_TO_DBL(coef));
                  newfj = fj;
               }
               else
               {
                  fj = SCIPfrac(scip, varsign[i] * QUAD_TO_DBL(coef));
                  newfj = SCIPfrac(scip, -varsign[i] * QUAD_TO_DBL(coef));
               }

               if( !ignoresol )
               {
                  solval = (sol == NULL ? SCIPvarGetLPSol(vars[v]) : SCIPgetSolVal(scip, sol, vars[v]));
                  viol = f0 - fj * (varsign[i] == +1 ? solval - bestlbs[i] : bestubs[i] - solval);
                  newviol = newf0 - newfj * (varsign[i] == -1 ? solval - bestlbs[i] : bestubs[i] - solval);
                  violgain = newviol - viol;
               }
               else
               {
                  /* todo: this should be done, this can improve the dualray significantly */
                  SCIPerrorMessage("Cannot handle closest bounds with ignoring the LP solution.\n");
                  return SCIP_INVALIDCALL;
               }

               /* prefer larger violations; for equal violations, prefer smaller f0 values since then the possibility that
                * we f_j > f_0 is larger and we may improve some coefficients in rounding
                */
               if( SCIPisGT(scip, violgain, bestviolgain)
                  || (SCIPisGE(scip, violgain, bestviolgain) && newf0 < bestnewf0) )
               {
                  besti = i;
                  bestviolgain = violgain;
                  bestnewf0 = newf0;
               }
            }
         }

         if( besti >= 0 )
         {
            SCIP_Real QUAD(coef);
            assert(besti < *nnz);
            assert(boundtype[besti] < 0);
            assert(!SCIPisInfinity(scip, -bestlbs[besti]));
            assert(!SCIPisInfinity(scip, bestubs[besti]));

            QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[besti]);
            QUAD_SCALE(coef, varsign[besti]);

            /* switch the complementation of this variable */
            SCIPquadprecSumDD(tmp, bestlbs[besti], - bestubs[besti]);
            SCIPquadprecProdQQ(tmp, tmp, coef);
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);

            if( varsign[besti] == +1 )
            {
               /* switch to upper bound */
               assert(bestubtypes[besti] < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
               boundtype[besti] = bestubtypes[besti];
               varsign[besti] = -1;
            }
            else
            {
               /* switch to lower bound */
               assert(bestlbtypes[besti] < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
               boundtype[besti] = bestlbtypes[besti];
               varsign[besti] = +1;
            }
            *localbdsused = *localbdsused || (boundtype[besti] == -2);
         }
      }
   }

  TERMINATE:

   /*free temporary memory */
   SCIPfreeBufferArray(scip, &selectedbounds);
   SCIPfreeBufferArray(scip, &bestubtypes);
   SCIPfreeBufferArray(scip, &bestlbtypes);
   SCIPfreeBufferArray(scip, &bestubs);
   SCIPfreeBufferArray(scip, &bestlbs);

   return SCIP_OKAY;
}

/** Calculate fractionalities \f$ f_0 := b - down(b), f_j := a^\prime_j - down(a^\prime_j) \f$, and derive MIR cut \f$ \tilde{a} \cdot x' \leq down(b) \f$
 * \f[
 * \begin{array}{rll}
 *  integers :&  \tilde{a}_j = down(a^\prime_j),                        & if \qquad f_j \leq f_0 \\
 *            &  \tilde{a}_j = down(a^\prime_j) + (f_j - f_0)/(1 - f_0),& if \qquad f_j >  f_0 \\
 *  continuous:& \tilde{a}_j = 0,                                       & if \qquad a^\prime_j \geq 0 \\
 *             & \tilde{a}_j = a^\prime_j/(1 - f_0),                    & if \qquad a^\prime_j <  0
 * \end{array}
 * \f]
 *
 *  Transform inequality back to \f$ \hat{a} \cdot x \leq rhs \f$:
 *
 *  (lb or ub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{if lb was used in transformation} \\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{if ub was used in transformation}
 * \end{array}
 * \f]
 *  and move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j \cdot lb_j = -\hat{a}_j \cdot lb_j,& \mbox{or} \\
 *     \tilde{a}_j \cdot ub_j = -\hat{a}_j \cdot ub_j &
 * \end{array}
 * \f]
 *  to the rhs.
 *
 *  (vlb or vub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - (bl_j \cdot zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{(vlb)} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{(vub)}
 * \end{array}
 * \f]
 *  move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j\, dl_j = -\hat{a}_j\, dl_j,& \mbox{or} \\
 *     \tilde{a}_j\, du_j = -\hat{a}_j\, du_j &
 * \end{array}
 * \f]
 *  to the rhs, and update the VB variable coefficients:
 * \f[
 * \begin{array}{ll}
 *    \hat{a}_{zl_j} := \hat{a}_{zl_j} - \tilde{a}_j\, bl_j = \hat{a}_{zl_j} - \hat{a}_j\, bl_j,& \mbox{or} \\
 *    \hat{a}_{zu_j} := \hat{a}_{zu_j} + \tilde{a}_j\, bu_j = \hat{a}_{zu_j} - \hat{a}_j\, bu_j &
 * \end{array}
 * \f]
 */
static
SCIP_RETCODE cutsRoundMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*RESTRICT    cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*RESTRICT cutrhs),          /**< pointer to right hand side of cut */
   int*RESTRICT          cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*RESTRICT          nnz,                /**< number of non-zeros in cut */
   int*RESTRICT          varsign,            /**< stores the sign of the transformed variable in summation */
   int*RESTRICT          boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub) */
   QUAD(SCIP_Real        f0)                 /**< fractional value of rhs */
   )
{
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(onedivoneminusf0);
   int i;
   int firstcontvar;
   SCIP_VAR** vars;
   int ndelcontvars;

   assert(QUAD_HI(cutrhs) != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(boundtype != NULL);
   assert(varsign != NULL);
   assert(0.0 < QUAD_TO_DBL(f0) && QUAD_TO_DBL(f0) < 1.0);

   SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0);
   SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0);

   /* Loop backwards to process integral variables first and be able to delete coefficients of integral variables
    * without destroying the ordering of the aggrrow's non-zeros.
    * (due to sorting in cutsTransformMIR the ordering is continuous before integral)
    */

   firstcontvar = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   vars = SCIPgetVars(scip);
#ifndef NDEBUG
   /*in debug mode check that all continuous variables of the aggrrow come before the integral variables */
   i = 0;
   while( i < *nnz && cutinds[i] >= firstcontvar )
      ++i;

   while( i < *nnz )
   {
      assert(cutinds[i] < firstcontvar);
      ++i;
   }
#endif

   for( i = *nnz - 1; i >= 0 && cutinds[i] < firstcontvar; --i )
   {
      SCIP_VAR* var;
      SCIP_Real QUAD(cutaj);
      int v;

      v = cutinds[i];
      assert(0 <= v && v < SCIPgetNVars(scip));

      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[i] == +1 || varsign[i] == -1);

      /* calculate the coefficient in the retransformed cut */
      {
         SCIP_Real QUAD(aj);
         SCIP_Real QUAD(downaj);
         SCIP_Real QUAD(fj);

         QUAD_ARRAY_LOAD(aj, cutcoefs, v);
         QUAD_SCALE(aj, varsign[i]);

         SCIPquadprecEpsFloorQ(downaj, aj, SCIPepsilon(scip)); /*lint !e666*/
         SCIPquadprecSumQQ(fj, aj, -downaj);
         assert(QUAD_TO_DBL(fj) >= -SCIPepsilon(scip) && QUAD_TO_DBL(fj) < 1.0);

         if( SCIPisLE(scip, QUAD_TO_DBL(fj), QUAD_TO_DBL(f0)) )
         {
            QUAD_ASSIGN_Q(cutaj, downaj);
         }
         else
         {
            SCIPquadprecSumQQ(tmp, fj, -f0);
            SCIPquadprecProdQQ(tmp, tmp, onedivoneminusf0);
            SCIPquadprecSumQQ(cutaj, tmp, downaj);
         }

         QUAD_SCALE(cutaj, varsign[i]);
      }

      /* remove zero cut coefficients from cut */
      if( EPSZ(QUAD_TO_DBL(cutaj), QUAD_EPSILON) )
      {
         QUAD_ASSIGN(cutaj, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, cutaj);
         --*nnz;
         cutinds[i] = cutinds[*nnz];
         continue;
      }

      QUAD_ARRAY_STORE(cutcoefs, v, cutaj);

      /* integral var uses standard bound */
      assert(boundtype[i] < 0);

      /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
      if( varsign[i] == +1 )
      {
         /* lower bound was used */
         if( boundtype[i] == -1 )
         {
            assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
            SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbGlobal(var));
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp); /* rhs += cutaj * SCIPvarGetLbGlobal(var) */
         }
         else
         {
            assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
            SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbLocal(var));
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp); /* rhs += cutaj * SCIPvarGetLbLocal(var) */
         }
      }
      else
      {
         /* upper bound was used */
         if( boundtype[i] == -1 )
         {
            assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
            SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbGlobal(var));
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp); /* rhs += cutaj * SCIPvarGetUbGlobal(var) */
         }
         else
         {
            assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
            SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbLocal(var));
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp); /* rhs += cutaj * SCIPvarGetUbLocal(var) */
         }
      }
   }

   /* now process the continuous variables; postpone deletetion of zeros till all continuous variables have been processed */
   ndelcontvars = 0;
   while( i >= ndelcontvars )
   {
      SCIP_VAR* var;
      SCIP_Real QUAD(cutaj);
      int v;

      v = cutinds[i];
      assert(0 <= v && v < SCIPgetNVars(scip));

      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[i] == +1 || varsign[i] == -1);
      assert( v >= firstcontvar );

      /* calculate the coefficient in the retransformed cut */
      {
         SCIP_Real QUAD(aj);

         QUAD_ARRAY_LOAD(aj, cutcoefs, v);

         if( QUAD_TO_DBL(aj) * varsign[i] >= 0.0 )
            QUAD_ASSIGN(cutaj, 0.0);
         else
            SCIPquadprecProdQQ(cutaj, onedivoneminusf0, aj); /* cutaj = varsign[i] * aj * onedivoneminusf0; // a^_j */
      }

      /* remove zero cut coefficients from cut; move a continuous var from the beginning
       * to the current position, so that all integral variables stay behind the continuous
       * variables
       */
      if( EPSZ(QUAD_TO_DBL(cutaj), QUAD_EPSILON) )
      {
         QUAD_ASSIGN(cutaj, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, cutaj);
         cutinds[i] = cutinds[ndelcontvars];
         varsign[i] = varsign[ndelcontvars];
         boundtype[i] = boundtype[ndelcontvars];
         ++ndelcontvars;
         continue;
      }

      QUAD_ARRAY_STORE(cutcoefs, v, cutaj);

      /* check for variable bound use */
      if( boundtype[i] < 0 )
      {
         /* standard bound */

         /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
         if( varsign[i] == +1 )
         {
            /* lower bound was used */
            if( boundtype[i] == -1 )
            {
               assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
               SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbGlobal(var));
               SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
            }
            else
            {
               assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
               SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbLocal(var));
               SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
            }
         }
         else
         {
            /* upper bound was used */
            if( boundtype[i] == -1 )
            {
               assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
               SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbGlobal(var));
               SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
            }
            else
            {
               assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
               SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbLocal(var));
               SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
            }
         }
      }
      else
      {
         SCIP_VAR** vbz;
         SCIP_Real* vbb;
         SCIP_Real* vbd;
         SCIP_Real QUAD(zcoef);
         int vbidx;
         int zidx;

         /* variable bound */
         vbidx = boundtype[i];

         /* change mirrhs and cutaj of integer variable z_j of variable bound */
         if( varsign[i] == +1 )
         {
            /* variable lower bound was used */
            assert(0 <= vbidx && vbidx < SCIPvarGetNVlbs(var));
            vbz = SCIPvarGetVlbVars(var);
            vbb = SCIPvarGetVlbCoefs(var);
            vbd = SCIPvarGetVlbConstants(var);
         }
         else
         {
            /* variable upper bound was used */
            assert(0 <= vbidx && vbidx < SCIPvarGetNVubs(var));
            vbz = SCIPvarGetVubVars(var);
            vbb = SCIPvarGetVubCoefs(var);
            vbd = SCIPvarGetVubConstants(var);
         }
         assert(SCIPvarIsActive(vbz[vbidx]));
         zidx = SCIPvarGetProbindex(vbz[vbidx]);
         assert(0 <= zidx && zidx < firstcontvar);

         SCIPquadprecProdQD(tmp, cutaj, vbd[vbidx]);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);

         SCIPquadprecProdQD(tmp, cutaj, vbb[vbidx]);
         QUAD_ARRAY_LOAD(zcoef, cutcoefs, zidx);

         /* update sparsity pattern */
         if( QUAD_HI(zcoef) == 0.0 )
            cutinds[(*nnz)++] = zidx;

         SCIPquadprecSumQQ(zcoef, zcoef, -tmp);
         QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
         QUAD_ARRAY_STORE(cutcoefs, zidx, zcoef);
         assert(QUAD_HI(zcoef) != 0.0);
      }

      /* advance to next variable */
      --i;
   }

   /* fill the empty position due to deleted continuous variables */
   if( ndelcontvars > 0 )
   {
      assert(ndelcontvars <= *nnz);
      *nnz -= ndelcontvars;
      if( *nnz < ndelcontvars )
      {
         BMScopyMemoryArray(cutinds, cutinds + ndelcontvars, *nnz);
      }
      else
      {
         BMScopyMemoryArray(cutinds, cutinds + *nnz, ndelcontvars);
      }
   }

   return SCIP_OKAY;
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row: \f$ a^\prime_r = scale * weight[r] * slacksign[r]. \f$
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 *  \f[
 *  \begin{array}{rll}
 *    integers : & \hat{a}_r = \tilde{a}_r = down(a^\prime_r),                      & \mbox{if}\qquad f_r <= f0 \\
 *               & \hat{a}_r = \tilde{a}_r = down(a^\prime_r) + (f_r - f0)/(1 - f0),& \mbox{if}\qquad f_r >  f0 \\
 *    continuous:& \hat{a}_r = \tilde{a}_r = 0,                                     & \mbox{if}\qquad a^\prime_r >= 0 \\
 *               & \hat{a}_r = \tilde{a}_r = a^\prime_r/(1 - f0),                   & \mbox{if}\qquad a^\prime_r <  0
 *  \end{array}
 *  \f]
 *
 *  Substitute \f$ \hat{a}_r \cdot s_r \f$ by adding \f$ \hat{a}_r \f$ times the slack's definition to the cut.
 */
static
SCIP_RETCODE cutsSubstituteMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   int*                  rowinds,            /**< sparsity pattern of used rows */
   int                   nrowinds,           /**< number of used rows */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   QUAD(SCIP_Real        f0)                 /**< fractional value of rhs */
   )
{  /*lint --e{715}*/
   SCIP_ROW** rows;
   SCIP_Real QUAD(onedivoneminusf0);
   int i;

   assert(scip != NULL);
   assert(weights != NULL || nrowinds == 0);
   assert(slacksign != NULL || nrowinds == 0);
   assert(rowinds != NULL || nrowinds == 0);
   assert(scale > 0.0);
   assert(cutcoefs != NULL);
   assert(QUAD_HI(cutrhs) != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(0.0 < QUAD_TO_DBL(f0) && QUAD_TO_DBL(f0) < 1.0);

   SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0);
   SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0);

   rows = SCIPgetLPRows(scip);
   for( i = 0; i < nrowinds; i++ )
   {
      SCIP_ROW* row;
      SCIP_Real ar;
      SCIP_Real downar;
      SCIP_Real QUAD(cutar);
      SCIP_Real QUAD(fr);
      SCIP_Real QUAD(tmp);
      SCIP_Real mul;
      int r;

      r = rowinds[i]; /*lint !e613*/
      assert(0 <= r && r < SCIPgetNLPRows(scip));
      assert(slacksign[i] == -1 || slacksign[i] == +1); /*lint !e613*/
      assert(!SCIPisZero(scip, weights[i])); /*lint !e613*/

      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* get the slack's coefficient a'_r in the aggregated row */
      ar = slacksign[i] * scale * weights[i]; /*lint !e613*/

      /* calculate slack variable's coefficient a^_r in the cut */
      if( row->integral
         && ((slacksign[i] == +1 && SCIPisFeasIntegral(scip, row->rhs - row->constant))
            || (slacksign[i] == -1 && SCIPisFeasIntegral(scip, row->lhs - row->constant))) ) /*lint !e613*/
      {
         /* slack variable is always integral:
          *    a^_r = a~_r = down(a'_r)                      , if f_r <= f0
          *    a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
          */
         downar = EPSFLOOR(ar, QUAD_EPSILON);
         SCIPquadprecSumDD(fr, ar, -downar);
         if( SCIPisLE(scip, QUAD_TO_DBL(fr), QUAD_TO_DBL(f0)) )
            QUAD_ASSIGN(cutar, downar);
         else
         {
            SCIPquadprecSumQQ(cutar, fr, -f0);
            SCIPquadprecProdQQ(cutar, cutar, onedivoneminusf0);
            SCIPquadprecSumQD(cutar, cutar, downar);
         }
      }
      else
      {
         /* slack variable is continuous:
          *    a^_r = a~_r = 0                               , if a'_r >= 0
          *    a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
          */
         if( ar >= 0.0 )
            continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
         else
            SCIPquadprecProdQD(cutar, onedivoneminusf0, ar);
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( EPSZ(QUAD_TO_DBL(cutar), QUAD_EPSILON) )
         continue;

      /* depending on the slack's sign, we have
       *   a*x + c + s == rhs  =>  s == - a*x - c + rhs,  or  a*x + c - s == lhs  =>  s == a*x + c - lhs
       * substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
       */
      mul = -slacksign[i] * QUAD_TO_DBL(cutar); /*lint !e613*/

      /* add the slack's definition multiplied with a^_j to the cut */
      SCIP_CALL( varVecAddScaledRowCoefsQuad(cutinds, cutcoefs, nnz, row, mul) );

      /* move slack's constant to the right hand side */
      if( slacksign[i] == +1 ) /*lint !e613*/
      {
         SCIP_Real QUAD(rowrhs);

         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
         assert(!SCIPisInfinity(scip, row->rhs));
         SCIPquadprecSumDD(rowrhs, row->rhs, -row->constant);
         if( row->integral )
         {
            /* the right hand side was implicitly rounded down in row aggregation */
            QUAD_ASSIGN(rowrhs, SCIPfloor(scip, QUAD_TO_DBL(rowrhs)));
         }
         SCIPquadprecProdQQ(tmp, cutar, rowrhs);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
      }
      else
      {
         SCIP_Real QUAD(rowlhs);

         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
         assert(!SCIPisInfinity(scip, -row->lhs));
         SCIPquadprecSumDD(rowlhs, row->lhs, -row->constant);
         if( row->integral )
         {
            /* the left hand side was implicitly rounded up in row aggregation */
            QUAD_ASSIGN(rowlhs, SCIPceil(scip, QUAD_TO_DBL(rowlhs)));
         }
         SCIPquadprecProdQQ(tmp, cutar, rowlhs);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
      }
   }

   /* relax rhs to zero, if it's very close to */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   return SCIP_OKAY;
}

/** calculates an MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because
 *  these rows cannot participate in an MIR cut.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcalcMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to the aggrrow; must be positive */
   SCIP_AGGRROW*         aggrrow,            /**< aggrrow to compute MIR cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store efficacy of cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut */
   )
{
   int i;
   int nvars;
   int* varsign;
   int* boundtype;
   SCIP_Real* tmpcoefs;

   SCIP_Real QUAD(rhs);
   SCIP_Real QUAD(downrhs);
   SCIP_Real QUAD(f0);
   SCIP_Bool freevariable;
   SCIP_Bool localbdsused;

   assert(aggrrow != NULL);
   assert(SCIPisPositive(scip, scale));
   assert(success != NULL);

   SCIPdebugMessage("calculating MIR cut (scale: %g)\n", scale);

   *success = FALSE;

   /* allocate temporary memory */
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtype, nvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &tmpcoefs, QUAD_ARRAY_SIZE(nvars)) );

   /* initialize cut with aggregation */
   *cutnnz = aggrrow->nnz;
   *cutislocal = aggrrow->local;

   SCIPquadprecProdQD(rhs, aggrrow->rhs, scale);

   if( *cutnnz > 0 )
   {
      BMScopyMemoryArray(cutinds, aggrrow->inds, *cutnnz);

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(coef);

         int k = aggrrow->inds[i];
         QUAD_ARRAY_LOAD(coef, aggrrow->vals, k);

         SCIPquadprecProdQD(coef, coef, scale);

         QUAD_ARRAY_STORE(tmpcoefs, k, coef);

         assert(QUAD_HI(coef) != 0.0);
      }

      /* Transform equation  a*x == b, lb <= x <= ub  into standard form
       *   a'*x' == b, 0 <= x' <= ub'.
       *
       * Transform variables (lb or ub):
       *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
       *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
       * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
       *
       * Transform variables (vlb or vub):
       *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
       *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
       * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
       *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
       *   a_{zu_j} := a_{zu_j} + a_j * bu_j
       */
      SCIP_CALL( cutsTransformMIR(scip, sol, boundswitch, usevbds, allowlocal, fixintegralrhs, FALSE,
            boundsfortrans, boundtypesfortrans, minfrac, maxfrac, tmpcoefs, QUAD(&rhs), cutinds, cutnnz, varsign, boundtype, &freevariable, &localbdsused) );
      assert(allowlocal || !localbdsused);
      *cutislocal = *cutislocal || localbdsused;

      if( freevariable )
         goto TERMINATE;
      SCIPdebug(printCutQuad(scip, sol, tmpcoefs, QUAD(rhs), cutinds, *cutnnz, FALSE, FALSE));
   }

   /* Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers :  a~_j = down(a'_j)                      , if f_j <= f_0
    *             a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
    * continuous: a~_j = 0                               , if a'_j >= 0
    *             a~_j = a'_j/(1 - f0)                   , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */
   SCIPquadprecEpsFloorQ(downrhs, rhs, SCIPepsilon(scip)); /*lint !e666*/

   SCIPquadprecSumQQ(f0, rhs, -downrhs);

   if( QUAD_TO_DBL(f0) < minfrac || QUAD_TO_DBL(f0) > maxfrac )
      goto TERMINATE;

   /* We multiply the coefficients of the base inequality roughly by scale/(1-f0).
    * If this gives a scalar that is very big, we better do not generate this cut.
    */
   if( REALABS(scale)/(1.0 - QUAD_TO_DBL(f0)) > MAXCMIRSCALE )
      goto TERMINATE;

   /* renormalize f0 value */
   SCIPquadprecSumDD(f0, QUAD_HI(f0), QUAD_LO(f0));

   QUAD_ASSIGN_Q(rhs, downrhs);

   if( *cutnnz > 0 )
   {
      SCIP_CALL( cutsRoundMIR(scip, tmpcoefs, QUAD(&rhs), cutinds, cutnnz, varsign, boundtype, QUAD(f0)) );
      SCIPdebug(printCutQuad(scip, sol, tmpcoefs, QUAD(rhs), cutinds, *cutnnz, FALSE, FALSE));
   }

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a^_r = a~_r = down(a'_r)                      , if f_r <= f0
    *               a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
    *   continuous: a^_r = a~_r = 0                               , if a'_r >= 0
    *               a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
    *
    * Substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
    */
   SCIP_CALL( cutsSubstituteMIR(scip, aggrrow->rowweights, aggrrow->slacksign, aggrrow->rowsinds,
                                aggrrow->nrows, scale, tmpcoefs, QUAD(&rhs), cutinds, cutnnz, QUAD(f0)) );
   SCIPdebug( printCutQuad(scip, sol, tmpcoefs, QUAD(rhs), cutinds, *cutnnz, FALSE, FALSE) );

   if( postprocess )
   {
      /* remove all nearly-zero coefficients from MIR row and relax the right hand side correspondingly in order to
       * prevent numerical rounding errors
       */
      SCIP_CALL( postprocessCutQuad(scip, *cutislocal, cutinds, tmpcoefs, cutnnz, QUAD(&rhs), success) );
   }
   else
   {
      *success = ! removeZerosQuad(scip, SCIPsumepsilon(scip), *cutislocal, tmpcoefs, QUAD(&rhs), cutnnz, cutinds);
   }

   SCIPdebug(printCutQuad(scip, sol, tmpcoefs, QUAD(rhs), cutinds, *cutnnz, FALSE, FALSE));

   if( *success )
   {
      *cutrhs = QUAD_TO_DBL(rhs);

      /* clean tmpcoefs and go back to double precision */
      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(coef);
         int j = cutinds[i];

         QUAD_ARRAY_LOAD(coef, tmpcoefs, j);

         cutcoefs[i] = QUAD_TO_DBL(coef);
         QUAD_ASSIGN(coef, 0.0);
         QUAD_ARRAY_STORE(tmpcoefs, j, coef);
      }

      if( cutefficacy != NULL )
         *cutefficacy = calcEfficacy(scip, sol, cutcoefs, *cutrhs, cutinds, *cutnnz);

      if( cutrank != NULL )
         *cutrank = aggrrow->rank + 1;
   }

  TERMINATE:
   if( !(*success) )
   {
      SCIP_Real QUAD(tmp);

      QUAD_ASSIGN(tmp, 0.0);
      for( i = 0; i < *cutnnz; ++i )
      {
         QUAD_ARRAY_STORE(tmpcoefs, cutinds[i], tmp);
      }
   }
   /* free temporary memory */
   SCIPfreeCleanBufferArray(scip, &tmpcoefs);
   SCIPfreeBufferArray(scip, &boundtype);
   SCIPfreeBufferArray(scip, &varsign);

   return SCIP_OKAY;
}

/** compute the efficacy of the MIR cut for the given values without computing the cut.
 *  This is used for the CMIR cut generation heuristic.
 */
static
SCIP_Real computeMIREfficacy(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*RESTRICT    coefs,              /**< array with coefficients in row */
   SCIP_Real*RESTRICT    solvals,            /**< solution values of variables in the row */
   SCIP_Real             rhs,                /**< right hand side of MIR cut */
   SCIP_Real             contactivity,       /**< aggregated activity of continuous variables in the row */
   SCIP_Real             contsqrnorm,        /**< squared norm of continuous variables */
   SCIP_Real             delta,              /**< delta value to compute the violation for */
   int                   nvars,              /**< number of variables in the row, i.e. the size of coefs and solvals arrays */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac             /**< maximal fractionality of rhs to produce MIR cut for */
   )
{
   int i;
   SCIP_Real f0pluseps;
   SCIP_Real f0;
   SCIP_Real onedivoneminusf0;
   SCIP_Real scale;
   SCIP_Real downrhs;
   SCIP_Real norm;
   SCIP_Real contscale;

   scale = 1.0 / delta;

   rhs *= scale;

   downrhs = SCIPfloor(scip, rhs);

   f0 = rhs - downrhs;

   if( f0 < minfrac || f0 > maxfrac )
      return 0.0;

   onedivoneminusf0 = 1.0 / (1.0 - f0);

   contscale = scale * onedivoneminusf0;

   /* We multiply the coefficients of the base inequality roughly by scale/(1-f0).
    * If this gives a scalar that is very big, we better do not generate this cut.
    */
   if( contscale > MAXCMIRSCALE )
      return 0.0;

   rhs = downrhs;
   rhs -= contscale * contactivity;
   norm = SQR(contscale) * contsqrnorm;

   assert(!SCIPisFeasZero(scip, f0));
   assert(!SCIPisFeasZero(scip, 1.0 - f0));

   f0pluseps = f0 + SCIPepsilon(scip);

   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real floorai = floor(scale * coefs[i]);
      SCIP_Real fi = (scale * coefs[i]) - floorai;

      if( fi > f0pluseps )
         floorai += (fi - f0) * onedivoneminusf0;

      rhs -= solvals[i] * floorai;
      norm += SQR(floorai);
   }

   norm = SQRT(norm);

   return - rhs / MAX(norm, 1e-6);
}

/** calculates an MIR cut out of an aggregation of LP rows
 *
 *  Given the aggregation, it is transformed to a mixed knapsack set via complementation (using bounds or variable bounds)
 *  Then, different scalings of the mkset are used to generate a MIR and the best is chosen.
 *  One of the steps of the MIR is to round the coefficients of the integer variables down,
 *  so one would prefer to have integer coefficients for integer variables which are far away from their bounds in the
 *  mkset.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcutGenerationHeuristicCMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int                   maxtestdelta,       /**< maximum number of deltas to test */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_AGGRROW*         aggrrow,            /**< aggrrow to compute MIR cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store efficacy of best cut; only cuts that are strictly better than the value of
                                              *   this efficacy on input to this function are returned */
   int*                  cutrank,            /**< pointer to return rank of generated cut (or NULL) */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid and efficacious cut was returned */
   )
{
   int i;
   int firstcontvar;
   int nvars;
   int intstart;
   int ntmpcoefs;
   int* varsign;
   int* boundtype;
   int* mksetinds;
   SCIP_Real* mksetcoefs;
   SCIP_Real QUAD(mksetrhs);
   int mksetnnz;
   SCIP_Real* bounddist;
   int* bounddistpos;
   int nbounddist;
   SCIP_Real* tmpcoefs;
   SCIP_Real* tmpvalues;
   SCIP_Real* deltacands;
   int ndeltacands;
   SCIP_Real bestdelta;
   SCIP_Real bestefficacy;
   SCIP_Real maxabsmksetcoef;
   SCIP_VAR** vars;
   SCIP_Bool freevariable;
   SCIP_Bool localbdsused;
   SCIP_Real contactivity;
   SCIP_Real contsqrnorm;

   assert(aggrrow != NULL);
   assert(aggrrow->nrows + aggrrow->nnz >= 1);
   assert(success != NULL);

   *success = FALSE;
   nvars = SCIPgetNVars(scip);
   firstcontvar = nvars - SCIPgetNContVars(scip);
   vars = SCIPgetVars(scip);

   /* allocate temporary memory */

   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtype, nvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &mksetcoefs, QUAD_ARRAY_SIZE(nvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mksetinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcoefs, nvars + aggrrow->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvalues, nvars + aggrrow->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &deltacands, aggrrow->nnz + 2) );
   /* we only compute bound distance for integer variables; we allocate an array of length aggrrow->nnz to store this, since
    * this is the largest number of integer variables. (in contrast to the number of total variables which can be 2 *
    * aggrrow->nnz variables: if all are continuous and we use variable bounds to completement, we introduce aggrrow->nnz
    * extra vars)
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &bounddist, aggrrow->nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounddistpos, aggrrow->nnz) );

   /* initialize mkset with aggregation */
   mksetnnz = aggrrow->nnz;
   QUAD_ASSIGN_Q(mksetrhs, aggrrow->rhs);

   BMScopyMemoryArray(mksetinds, aggrrow->inds, mksetnnz);

   for( i = 0; i < mksetnnz; ++i )
   {
      int j = mksetinds[i];
      SCIP_Real QUAD(coef);
      QUAD_ARRAY_LOAD(coef, aggrrow->vals, j);
      QUAD_ARRAY_STORE(mksetcoefs, j, coef);
      assert(QUAD_HI(coef) != 0.0);
   }

   *cutislocal = aggrrow->local;

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a'*x' == b, 0 <= x' <= ub'.
    *
    * Transform variables (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
    * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
    *
    * Transform variables (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
    * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
    *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
    *   a_{zu_j} := a_{zu_j} + a_j * bu_j
    */
   SCIP_CALL( cutsTransformMIR(scip, sol, boundswitch, usevbds, allowlocal, FALSE, FALSE,
         boundsfortrans, boundtypesfortrans, minfrac, maxfrac, mksetcoefs, QUAD(&mksetrhs), mksetinds, &mksetnnz, varsign, boundtype, &freevariable, &localbdsused) );

   assert(allowlocal || !localbdsused);

   if( freevariable )
      goto TERMINATE;

   SCIPdebugMessage("transformed aggrrow row:\n");
   SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(mksetrhs), mksetinds, mksetnnz, FALSE, FALSE));

   /* found positions of integral variables that are strictly between their bounds */
   maxabsmksetcoef = -1.0;
   nbounddist = 0;

   for( i = mksetnnz - 1; i >= 0 && mksetinds[i] < firstcontvar; --i )
   {
      SCIP_VAR* var = vars[mksetinds[i]];
      SCIP_Real primsol = SCIPgetSolVal(scip, sol, var);
      SCIP_Real lb = SCIPvarGetLbLocal(var);
      SCIP_Real ub = SCIPvarGetUbLocal(var);
      SCIP_Real QUAD(coef);
      SCIP_Real absmksetcoef;

      QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[i]);

      absmksetcoef = REALABS(QUAD_TO_DBL(coef));

      maxabsmksetcoef = MAX(absmksetcoef, maxabsmksetcoef);

      if( SCIPisEQ(scip, primsol, lb) || SCIPisEQ(scip, primsol, ub) )
         continue;

      bounddist[nbounddist] = MIN(ub - primsol, primsol - lb);
      bounddistpos[nbounddist] = i;
      deltacands[nbounddist] = absmksetcoef;
      ++nbounddist;
   }

   /* no fractional variable; so abort here */
   if( nbounddist == 0 )
      goto TERMINATE;

   intstart = i + 1;
   ndeltacands = nbounddist;

   SCIPsortDownRealRealInt(bounddist, deltacands, bounddistpos, nbounddist);

   /* also test 1.0 and maxabsmksetcoef + 1.0 as last delta values */
   if( maxabsmksetcoef != -1.0 )
   {
      deltacands[ndeltacands++] = maxabsmksetcoef + 1.0;
   }

   deltacands[ndeltacands++] = 1.0;

   maxtestdelta = MIN(ndeltacands, maxtestdelta);

   /* For each delta
    * Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers :  a~_j = down(a'_j)                      , if f_j <= f_0
    *             a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
    * continuous: a~_j = 0                               , if a'_j >= 0
    *             a~_j = a'_j/(1 - f0)                   , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */

   ntmpcoefs = 0;
   for( i = intstart; i < mksetnnz; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real solval;
      SCIP_Real QUAD(coef);

      var = vars[mksetinds[i]];

      /* get the soltion value of the continuous variable */
      solval = SCIPgetSolVal(scip, sol, var);

      /* now compute the solution value in the transform space considering complementation */
      if( boundtype[i] == -1 )
      {
         /* variable was complemented with global (simple) bound */
         if( varsign[i] == -1 )
            solval = SCIPvarGetUbGlobal(var) - solval;
         else
            solval = solval - SCIPvarGetLbGlobal(var);
      }
      else
      {
         assert(boundtype[i] == -2);

         /* variable was complemented with local (simple) bound */
         if( varsign[i] == -1 )
            solval = SCIPvarGetUbLocal(var) - solval;
         else
            solval = solval - SCIPvarGetLbLocal(var);
      }

      tmpvalues[ntmpcoefs] = solval;
      QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[i]);
      tmpcoefs[ntmpcoefs] = varsign[i] * QUAD_TO_DBL(coef);
      ++ntmpcoefs;
   }

   assert(ntmpcoefs == mksetnnz - intstart);

   contactivity = 0.0;
   contsqrnorm = 0.0;
   for( i = 0; i < intstart; ++i )
   {
      SCIP_Real solval;
      SCIP_Real QUAD(mksetcoef);

      QUAD_ARRAY_LOAD(mksetcoef, mksetcoefs, mksetinds[i]);

      if( varsign[i] * QUAD_TO_DBL(mksetcoef) >= 0.0 )
         continue;

      /* get the soltion value of the continuous variable */
      solval = SCIPgetSolVal(scip, sol, vars[mksetinds[i]]);

      /* now compute the solution value in the transform space considering complementation */
      switch( boundtype[i] )
      {
         case -1:
            /* variable was complemented with global (simple) bound */
            if( varsign[i] == -1 )
               solval = SCIPvarGetUbGlobal(vars[mksetinds[i]]) - solval;
            else
               solval = solval - SCIPvarGetLbGlobal(vars[mksetinds[i]]);
            break;
         case -2:
            /* variable was complemented with local (simple) bound */
            if( varsign[i] == -1 )
               solval = SCIPvarGetUbLocal(vars[mksetinds[i]]) - solval;
            else
               solval = solval - SCIPvarGetLbLocal(vars[mksetinds[i]]);
            break;
         default:
            /* variable was complemented with a variable bound */
            if( varsign[i] == -1 )
            {
               SCIP_Real coef;
               SCIP_Real constant;
               SCIP_Real vbdsolval;

               coef = SCIPvarGetVubCoefs(vars[mksetinds[i]])[boundtype[i]];
               constant = SCIPvarGetVubConstants(vars[mksetinds[i]])[boundtype[i]];
               vbdsolval = SCIPgetSolVal(scip, sol, SCIPvarGetVubVars(vars[mksetinds[i]])[boundtype[i]]);

               solval = (coef * vbdsolval + constant) - solval;
            }
            else
            {
               SCIP_Real coef;
               SCIP_Real constant;
               SCIP_Real vbdsolval;

               coef = SCIPvarGetVlbCoefs(vars[mksetinds[i]])[boundtype[i]];
               constant = SCIPvarGetVlbConstants(vars[mksetinds[i]])[boundtype[i]];
               vbdsolval = SCIPgetSolVal(scip, sol, SCIPvarGetVlbVars(vars[mksetinds[i]])[boundtype[i]]);

               solval = solval - (coef * vbdsolval + constant);
            }
      }

      contactivity += solval * (QUAD_TO_DBL(mksetcoef) * varsign[i]);
      contsqrnorm += QUAD_TO_DBL(mksetcoef) * QUAD_TO_DBL(mksetcoef);
   }

   {
      SCIP_ROW** rows;

      rows = SCIPgetLPRows(scip);

      for( i = 0; i < aggrrow->nrows; ++i )
      {
         SCIP_ROW* row;
         SCIP_Real slackval;

         row = rows[aggrrow->rowsinds[i]];

         if( (aggrrow->rowweights[i] * aggrrow->slacksign[i]) >= 0.0 && !row->integral )
            continue;

         /* compute solution value of slack variable */
         slackval = SCIPgetRowSolActivity(scip, row, sol);

         if( aggrrow->slacksign[i] == +1 )
         {
            /* right hand side */
            assert(!SCIPisInfinity(scip, row->rhs));

            slackval = row->rhs - slackval;
         }
         else
         {
            /* left hand side */
            assert(aggrrow->slacksign[i] == -1);
            assert(!SCIPisInfinity(scip, -row->lhs));

            slackval = slackval - row->lhs;
         }

         if( row->integral )
         {
            /* if row is integral add variable to tmp arrays */
            tmpvalues[ntmpcoefs] = slackval;
            tmpcoefs[ntmpcoefs] = aggrrow->rowweights[i] * aggrrow->slacksign[i];
            ++ntmpcoefs;
         }
         else
         {
            SCIP_Real slackcoeff = (aggrrow->rowweights[i] * aggrrow->slacksign[i]);

            /* otherwise add it to continuous activity */
            contactivity += slackval * slackcoeff;
            contsqrnorm += SQR(slackcoeff);
         }
      }
   }

   /* try all candidates for delta and remember best */
   bestdelta = SCIP_INVALID;
   bestefficacy = -SCIPinfinity(scip);

   for( i = 0; i < maxtestdelta; ++i )
   {
      int j;
      SCIP_Real efficacy;

      /* check if we have seen this value of delta before */
      SCIP_Bool deltaseenbefore = FALSE;
      for( j = 0; j < i; ++j )
      {
         if( SCIPisEQ(scip, deltacands[i], deltacands[j]) )
         {
            deltaseenbefore = TRUE;
            break;
         }
      }

      /* skip this delta value and allow one more delta value if available */
      if( deltaseenbefore )
      {
         maxtestdelta = MIN(maxtestdelta + 1, ndeltacands);
         continue;
      }

      efficacy = computeMIREfficacy(scip, tmpcoefs, tmpvalues, QUAD_TO_DBL(mksetrhs), contactivity, contsqrnorm, deltacands[i], ntmpcoefs, minfrac, maxfrac);

      if( efficacy > bestefficacy )
      {
         bestefficacy = efficacy;
         bestdelta = deltacands[i];
      }
   }

   /* no delta was found that yielded any cut */
   if( bestdelta == SCIP_INVALID ) /*lint !e777*/
      goto TERMINATE;

   /* try bestdelta divided by 2, 4 and 8 */
   for( i = 2; i <= 8 ; i *= 2 )
   {
      SCIP_Real efficacy;
      SCIP_Real delta;

      delta = bestdelta / i;

      efficacy = computeMIREfficacy(scip, tmpcoefs, tmpvalues, QUAD_TO_DBL(mksetrhs), contactivity, contsqrnorm, delta, ntmpcoefs, minfrac, maxfrac);

      if( efficacy >= bestefficacy )
      {
         bestefficacy = efficacy;
         bestdelta = delta;
      }
   }

   /* try to improve efficacy by switching complementation of integral variables that are not at their bounds
    * in order of non-increasing bound distance
    */
   for( i = 0; i < nbounddist; ++i )
   {
      int k;
      SCIP_Real newefficacy;
      SCIP_Real QUAD(newrhs);
      SCIP_Real bestlb;
      SCIP_Real bestub;
      SCIP_Real oldsolval;
      int bestlbtype;
      int bestubtype;

      k = bounddistpos[i];

      SCIP_CALL( findBestLb(scip, vars[mksetinds[k]], sol, FALSE, allowlocal, &bestlb, &bestlbtype) );

      if( SCIPisInfinity(scip, -bestlb) )
         continue;

      SCIP_CALL( findBestUb(scip, vars[mksetinds[k]], sol, FALSE, allowlocal, &bestub, &bestubtype) );

      if( SCIPisInfinity(scip, bestub) )
         continue;

      /* switch the complementation of this variable */
#ifndef NDEBUG
      {
         SCIP_Real QUAD(coef);
         QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[k]);
         assert(SCIPisEQ(scip, tmpcoefs[k - intstart], varsign[k] * QUAD_TO_DBL(coef)));
      }
#endif

      /* compute this: newrhs = mksetrhs + tmpcoefs[k - intstart] * (bestlb - bestub); */
      SCIPquadprecSumQD(newrhs, mksetrhs, tmpcoefs[k - intstart] * (bestlb - bestub));
      tmpcoefs[k - intstart] = -tmpcoefs[k - intstart];

      oldsolval = tmpvalues[k - intstart];
      tmpvalues[k - intstart] = varsign[k] == +1 ? bestub - SCIPgetSolVal(scip, sol, vars[mksetinds[k]]) : SCIPgetSolVal(scip, sol, vars[mksetinds[k]]) - bestlb;

      /* compute new violation */
      newefficacy = computeMIREfficacy(scip, tmpcoefs, tmpvalues, QUAD_TO_DBL(newrhs), contactivity, contsqrnorm, bestdelta, ntmpcoefs, minfrac, maxfrac);

      /* check if violaton was increased */
      if( newefficacy > bestefficacy )
      {
         /* keep change of complementation */
         bestefficacy = newefficacy;
         QUAD_ASSIGN_Q(mksetrhs, newrhs);

         if( varsign[k] == +1 )
         {
            /* switch to upper bound */
            assert(bestubtype < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
            boundtype[k] = bestubtype;
            varsign[k] = -1;
         }
         else
         {
            /* switch to lower bound */
            assert(bestlbtype < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
            boundtype[k] = bestlbtype;
            varsign[k] = +1;
         }

         localbdsused = localbdsused || (boundtype[k] == -2);
      }
      else
      {
         /* undo the change of the complementation */
         tmpcoefs[k - intstart] = -tmpcoefs[k - intstart];
         tmpvalues[k - intstart] = oldsolval;
      }
   }

   if( bestefficacy > 0.0 )
   {
      SCIP_Real mirefficacy;
      SCIP_Real QUAD(downrhs);
      SCIP_Real QUAD(f0);
      SCIP_Real scale;

      scale = 1.0 / bestdelta;
      SCIPquadprecProdQD(mksetrhs, mksetrhs, scale);

      SCIPquadprecEpsFloorQ(downrhs, mksetrhs, SCIPepsilon(scip)); /*lint !e666*/
      SCIPquadprecSumQQ(f0, mksetrhs, -downrhs);

      /* renormaliize f0 value */
      SCIPquadprecSumDD(f0, QUAD_HI(f0), QUAD_LO(f0));

      for( i = 0; i < mksetnnz; ++i )
      {
         SCIP_Real QUAD(coef);

         QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[i]);
         SCIPquadprecProdQD(coef, coef, scale);
         QUAD_ARRAY_STORE(mksetcoefs, mksetinds[i], coef);
      }
      SCIPdebugMessage("applied best scale (=%.13g):\n", scale);
      SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(mksetrhs), mksetinds, mksetnnz, FALSE, FALSE));

      QUAD_ASSIGN_Q(mksetrhs, downrhs);

      SCIP_CALL( cutsRoundMIR(scip, mksetcoefs, QUAD(&mksetrhs), mksetinds, &mksetnnz, varsign, boundtype, QUAD(f0)) );

      SCIPdebugMessage("rounded MIR cut:\n");
      SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(mksetrhs), mksetinds, mksetnnz, FALSE, FALSE));

      /* substitute aggregated slack variables:
       *
       * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
       * variable only appears in its own row:
       *    a'_r = scale * weight[r] * slacksign[r].
       *
       * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
       *   integers :  a^_r = a~_r = down(a'_r)                      , if f_r <= f0
       *               a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
       *   continuous: a^_r = a~_r = 0                               , if a'_r >= 0
       *               a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
       *
       * Substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
       */
      SCIP_CALL( cutsSubstituteMIR(scip, aggrrow->rowweights, aggrrow->slacksign, aggrrow->rowsinds,
                                   aggrrow->nrows, scale, mksetcoefs, QUAD(&mksetrhs), mksetinds, &mksetnnz, QUAD(f0)) );

      SCIPdebugMessage("substituted slacks in MIR cut:\n");
      SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(mksetrhs), mksetinds, mksetnnz, FALSE, FALSE));

#ifndef NDEBUG
      {
         SCIP_Real efficacy = -QUAD_TO_DBL(mksetrhs);
         for( i = 0; i < mksetnnz; ++i )
         {
            SCIP_Real QUAD(coef);
            QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[i]);
            efficacy += QUAD_TO_DBL(coef) * SCIPgetSolVal(scip, sol, vars[mksetinds[i]]);
         }

         if( !EPSZ(SCIPrelDiff(efficacy, bestefficacy), 1e-4) )
         {
            SCIPdebugMessage("efficacy of cmir cut is different than expected efficacy: %f != %f\n", efficacy, bestefficacy);
         }
      }
#endif

      *cutislocal = *cutislocal || localbdsused;

      /* remove all nearly-zero coefficients from MIR row and relax the right hand side correspondingly in order to
       * prevent numerical rounding errors
       */
      if( postprocess )
      {
         SCIP_CALL( postprocessCutQuad(scip, *cutislocal, mksetinds, mksetcoefs, &mksetnnz, QUAD(&mksetrhs), success) );
      }
      else
      {
         *success = ! removeZerosQuad(scip, SCIPsumepsilon(scip), *cutislocal, mksetcoefs, QUAD(&mksetrhs), mksetinds, &mksetnnz);
      }

      SCIPdebugMessage("post-processed cut (success = %s):\n", *success ? "TRUE" : "FALSE");
      SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(mksetrhs), mksetinds, mksetnnz, FALSE, FALSE));

      if( *success )
      {
         mirefficacy = calcEfficacyDenseStorageQuad(scip, sol, mksetcoefs, QUAD_TO_DBL(mksetrhs), mksetinds, mksetnnz);

         if( SCIPisEfficacious(scip, mirefficacy) && mirefficacy > *cutefficacy )
         {
            BMScopyMemoryArray(cutinds, mksetinds, mksetnnz);
            for( i = 0; i < mksetnnz; ++i )
            {
               SCIP_Real QUAD(coef);
               int j = cutinds[i];

               QUAD_ARRAY_LOAD(coef, mksetcoefs, j);

               cutcoefs[i] = QUAD_TO_DBL(coef);
               QUAD_ASSIGN(coef, 0.0);
               QUAD_ARRAY_STORE(mksetcoefs, j, coef);
            }
            *cutnnz = mksetnnz;
            *cutrhs = QUAD_TO_DBL(mksetrhs);
            *cutefficacy = mirefficacy;
            if( cutrank != NULL )
               *cutrank = aggrrow->rank + 1;
            *cutislocal = *cutislocal || localbdsused;
         }
         else
         {
            *success = FALSE;
         }
      }
   }

  TERMINATE:
   /* if we aborted early we need to clean the mksetcoefs */
   if( !(*success) )
   {
      SCIP_Real QUAD(tmp);
      QUAD_ASSIGN(tmp, 0.0);

      for( i = 0; i < mksetnnz; ++i )
      {
         QUAD_ARRAY_STORE(mksetcoefs, mksetinds[i], tmp);
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &bounddistpos);
   SCIPfreeBufferArray(scip, &bounddist);
   SCIPfreeBufferArray(scip, &deltacands);
   SCIPfreeBufferArray(scip, &tmpvalues);
   SCIPfreeBufferArray(scip, &tmpcoefs);
   SCIPfreeBufferArray(scip, &mksetinds);
   SCIPfreeCleanBufferArray(scip, &mksetcoefs);
   SCIPfreeBufferArray(scip, &boundtype);
   SCIPfreeBufferArray(scip, &varsign);

   return SCIP_OKAY;
}

/* =========================================== flow cover =========================================== */

#define NO_EXACT_KNAPSACK

#ifndef NO_EXACT_KNAPSACK
#define MAXDNOM                  1000LL
#define MINDELTA                  1e-03
#define MAXDELTA                  1e-09
#define MAXSCALE                 1000.0
#define MAXDYNPROGSPACE         1000000
#endif

#define MAXABSVBCOEF               1e+5 /**< maximal absolute coefficient in variable bounds used for snf relaxation */
#define MAXBOUND                  1e+10   /**< maximal value of normal bounds used for snf relaxation */

/** structure that contains all data required to perform the sequence independent lifting
 */
typedef
struct LiftingData
{
   SCIP_Real*            M;                  /**< \f$ M_0 := 0.0 \f$ and \f$ M_i := M_i-1 + m_i \f$ */
   SCIP_Real*            m;                  /**< non-increasing array of variable upper bound coefficients
                                              *   for all variables in \f$ C^{++} \f$  and \f$ L^- \f$,
                                              *   where \f$ C = C^+ \cup C^- \f$ is the flowcover and
                                              *   \f$ C^{++} := \{ j \in C^+ \mid u_j > \lambda \} \f$
                                              *   \f$ L^- := \{ j \in (N^- \setminus C^-) \mid u_j > \lambda \} \f$
                                              */
   int                   r;                  /**< size of array m */
   int                   t;                  /**< index of smallest value in m that comes from a variable in \f$ C^{++} \f$ */
   SCIP_Real             d1;                 /**< right hand side of single-node-flow set plus the sum of all \f$ u_j \f$ for \f$ j \in C^- \f$ */
   SCIP_Real             d2;                 /**< right hand side of single-node-flow set plus the sum of all \f$ u_j \f$ for \f$ j \in N^- \f$ */
   SCIP_Real             lambda;             /**< excess of the flowcover */
   SCIP_Real             mp;                 /**< smallest variable bound coefficient of variable in \f$ C^{++} (min_{j \in C++} u_j) \f$ */
   SCIP_Real             ml;                 /**< \f$ ml := min(\lambda, \sum_{j \in C^+ \setminus C^{++}} u_j) \f$ */
} LIFTINGDATA;

/** structure that contains all the data that defines the single-node-flow relaxation of an aggregation row */
typedef
struct SNF_Relaxation
{
   int*                  transvarcoefs;      /**< coefficients of all vars in relaxed set */
   SCIP_Real*            transbinvarsolvals; /**< sol val of bin var in vub of all vars in relaxed set */
   SCIP_Real*            transcontvarsolvals;/**< sol val of all real vars in relaxed set */
   SCIP_Real*            transvarvubcoefs;   /**< coefficient in vub of all vars in relaxed set */
   int                   ntransvars;         /**< number of vars in relaxed set */
   SCIP_Real             transrhs;           /**< rhs in relaxed set */
   int*                  origbinvars;        /**< associated original binary var for all vars in relaxed set */
   int*                  origcontvars;       /**< associated original continuous var for all vars in relaxed set */
   SCIP_Real*            aggrcoefsbin;       /**< aggregation coefficient of the original binary var used to define the
                                              *   continuous variable in the relaxed set */
   SCIP_Real*            aggrcoefscont;      /**< aggregation coefficient of the original continous var used to define the
                                              *   continuous variable in the relaxed set */
   SCIP_Real*            aggrconstants;      /**< aggregation constant used to define the continuous variable in the relaxed set */
} SNF_RELAXATION;

/** get solution value and index of variable lower bound (with binary variable) which is closest to the current LP
 *  solution value of a given variable; candidates have to meet certain criteria in order to ensure the nonnegativity
 *  of the variable upper bound imposed on the real variable in the 0-1 single node flow relaxation associated with the
 *  given variable
 */
static
SCIP_RETCODE getClosestVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< given active problem variable */
   SCIP_SOL*             sol,                /**< solution to use for variable bound; NULL for LP solution */
   SCIP_Real*            rowcoefs,           /**< (dense) array of coefficients of row */
   int8_t*               binvarused,         /**< array that stores if a binary variable was already used (+1)
                                              *   was not used (0) or was not used but is contained in the row (-1)
                                              */
   SCIP_Real             bestsub,            /**< closest simple upper bound of given variable */
   SCIP_Real             rowcoef,            /**< coefficient of given variable in current row */
   SCIP_Real*            closestvlb,         /**< pointer to store the LP sol value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest vlb; -1 if no vlb was found */
   )
{
   int nvlbs;
   int nbinvars;

   assert(scip != NULL);
   assert(var != NULL);
   assert(bestsub == SCIPvarGetUbGlobal(var) || bestsub == SCIPvarGetUbLocal(var)); /*lint !e777*/
   assert(!SCIPisInfinity(scip, bestsub));
   assert(!EPSZ(rowcoef, QUAD_EPSILON));
   assert(rowcoefs != NULL);
   assert(binvarused != NULL);
   assert(closestvlb != NULL);
   assert(closestvlbidx != NULL);

   nvlbs = SCIPvarGetNVlbs(var);
   nbinvars = SCIPgetNBinVars(scip);

   *closestvlbidx = -1;
   *closestvlb = -SCIPinfinity(scip);
   if( nvlbs > 0 )
   {
      SCIP_VAR** vlbvars;
      SCIP_Real* vlbcoefs;
      SCIP_Real* vlbconsts;
      int i;

      vlbvars = SCIPvarGetVlbVars(var);
      vlbcoefs = SCIPvarGetVlbCoefs(var);
      vlbconsts = SCIPvarGetVlbConstants(var);

      for( i = 0; i < nvlbs; i++ )
      {
         SCIP_Real rowcoefbinvar;
         SCIP_Real val1;
         SCIP_Real val2;
         SCIP_Real vlbsol;
         SCIP_Real rowcoefsign;
         int probidxbinvar;

         if( bestsub > vlbconsts[i] )
            continue;

         /* for numerical reasons, ignore variable bounds with large absolute coefficient and
          * those which lead to an infinite variable bound coefficient (val2) in snf relaxation
          */
         if( REALABS(vlbcoefs[i]) > MAXABSVBCOEF  )
            continue;

         /* use only variable lower bounds l~_i * x_i + d_i with x_i binary which are active */
         probidxbinvar = SCIPvarGetProbindex(vlbvars[i]);

         /* if the variable is not active the problem index is -1, so we cast to unsigned int before the comparison which
          * ensures that the problem index is between 0 and nbinvars - 1
          */
         if( (unsigned int)probidxbinvar >= (unsigned int)nbinvars )
            continue;

         assert(SCIPvarIsBinary(vlbvars[i]));

         /* check if current variable lower bound l~_i * x_i + d_i imposed on y_j meets the following criteria:
          * (let a_j  = coefficient of y_j in current row,
          *      u_j  = closest simple upper bound imposed on y_j,
          *      c_i  = coefficient of x_i in current row)
          *   0. no other non-binary variable y_k has used a variable bound with x_i to get transformed variable y'_k yet
          * if a_j > 0:
          *   1. u_j <= d_i
          *   2. a_j ( u_j - d_i ) + c_i <= 0
          *   3. a_j l~_i + c_i <= 0
          * if a_j < 0:
          *   1. u_j <= d_i
          *   2. a_j ( u_j - d_i ) + c_i >= 0
          *   3. a_j l~_i + c_i >= 0
          */

         /* has already been used in the SNF relaxation */
         if( binvarused[probidxbinvar] == 1 )
            continue;

         /* get the row coefficient */
         {
            SCIP_Real QUAD(tmp);
            QUAD_ARRAY_LOAD(tmp, rowcoefs, probidxbinvar);
            rowcoefbinvar = QUAD_TO_DBL(tmp);
         }
         rowcoefsign = COPYSIGN(1.0, rowcoef);

         val2 = rowcoefsign * ((rowcoef * vlbcoefs[i]) + rowcoefbinvar);

         /* variable lower bound does not meet criteria */
         if( val2 > 0.0 || SCIPisInfinity(scip, -val2) )
            continue;

         val1 = rowcoefsign * ((rowcoef * (bestsub - vlbconsts[i])) + rowcoefbinvar);

         /* variable lower bound does not meet criteria */
         if( val1 > 0.0 )
            continue;

         vlbsol = vlbcoefs[i] * SCIPgetSolVal(scip, sol, vlbvars[i]) + vlbconsts[i];
         if( vlbsol > *closestvlb )
         {
            *closestvlb = vlbsol;
            *closestvlbidx = i;
         }
         assert(*closestvlbidx >= 0);

      }
   }

   return SCIP_OKAY;
}

/** get LP solution value and index of variable upper bound (with binary variable) which is closest to the current LP
 *  solution value of a given variable; candidates have to meet certain criteria in order to ensure the nonnegativity
 *  of the variable upper bound imposed on the real variable in the 0-1 single node flow relaxation associated with the
 *  given variable
 */
static
SCIP_RETCODE getClosestVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< given active problem variable */
   SCIP_SOL*             sol,                /**< solution to use for variable bound; NULL for LP solution */
   SCIP_Real*            rowcoefs,           /**< (dense) array of coefficients of row */
   int8_t*               binvarused,         /**< array that stores if a binary variable was already used (+1)
                                              *   was not used (0) or was not used but is contained in the row (-1)
                                              */
   SCIP_Real             bestslb,            /**< closest simple lower bound of given variable */
   SCIP_Real             rowcoef,            /**< coefficient of given variable in current row */
   SCIP_Real*            closestvub,         /**< pointer to store the LP sol value of the closest variable upper bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest vub; -1 if no vub was found */
   )
{
   int nvubs;
   int nbinvars;

   assert(scip != NULL);
   assert(var != NULL);
   assert(bestslb == SCIPvarGetLbGlobal(var) || bestslb == SCIPvarGetLbLocal(var)); /*lint !e777*/
   assert(!SCIPisInfinity(scip, - bestslb));
   assert(!EPSZ(rowcoef, QUAD_EPSILON));
   assert(rowcoefs != NULL);
   assert(binvarused != NULL);
   assert(closestvub != NULL);
   assert(closestvubidx != NULL);

   nvubs = SCIPvarGetNVubs(var);
   nbinvars = SCIPgetNBinVars(scip);

   *closestvubidx = -1;
   *closestvub = SCIPinfinity(scip);
   if( nvubs > 0 )
   {
      SCIP_VAR** vubvars;
      SCIP_Real* vubcoefs;
      SCIP_Real* vubconsts;
      int i;

      vubvars = SCIPvarGetVubVars(var);
      vubcoefs = SCIPvarGetVubCoefs(var);
      vubconsts = SCIPvarGetVubConstants(var);

      for( i = 0; i < nvubs; i++ )
      {
         SCIP_Real rowcoefbinvar;
         SCIP_Real val1;
         SCIP_Real val2;
         SCIP_Real vubsol;
         SCIP_Real rowcoefsign;
         int probidxbinvar;

         if( bestslb < vubconsts[i] )
            continue;

         /* for numerical reasons, ignore variable bounds with large absolute coefficient and
          * those which lead to an infinite variable bound coefficient (val2) in snf relaxation
          */
         if( REALABS(vubcoefs[i]) > MAXABSVBCOEF  )
            continue;

         /* use only variable upper bound u~_i * x_i + d_i with x_i binary and which are active */
         probidxbinvar = SCIPvarGetProbindex(vubvars[i]);

         /* if the variable is not active the problem index is -1, so we cast to unsigned int before the comparison which
          * ensures that the problem index is between 0 and nbinvars - 1
          */
         if( (unsigned int)probidxbinvar >= (unsigned int)nbinvars )
            continue;

         assert(SCIPvarIsBinary(vubvars[i]));


         /* checks if current variable upper bound u~_i * x_i + d_i meets the following criteria
          * (let a_j  = coefficient of y_j in current row,
          *      l_j  = closest simple lower bound imposed on y_j,
          *      c_i  = coefficient of x_i in current row)
          *   0. no other non-binary variable y_k has used a variable bound with x_i to get transformed variable y'_k
          * if a > 0:
          *   1. l_j >= d_i
          *   2. a_j ( l_i - d_i ) + c_i >= 0
          *   3. a_j u~_i + c_i >= 0
          * if a < 0:
          *   1. l_j >= d_i
          *   2. a_j ( l_j - d_i ) + c_i <= 0
          *   3. a_j u~_i + c_i <= 0
          */

         /* has already been used in the SNF relaxation */
         if( binvarused[probidxbinvar] == 1 )
            continue;

         /* get the row coefficient */
         {
            SCIP_Real QUAD(tmp);
            QUAD_ARRAY_LOAD(tmp, rowcoefs, probidxbinvar);
            rowcoefbinvar = QUAD_TO_DBL(tmp);
         }
         rowcoefsign = COPYSIGN(1.0, rowcoef);

         val2 = rowcoefsign * ((rowcoef * vubcoefs[i]) + rowcoefbinvar);

         /* variable upper bound does not meet criteria */
         if( val2 < 0.0 || SCIPisInfinity(scip, val2) )
            continue;

         val1 = rowcoefsign * ((rowcoef * (bestslb - vubconsts[i])) + rowcoefbinvar);

         /* variable upper bound does not meet criteria */
         if( val1 < 0.0 )
            continue;

         vubsol = vubcoefs[i] * SCIPgetSolVal(scip, sol, vubvars[i]) + vubconsts[i];
         if( vubsol < *closestvub )
         {
            *closestvub = vubsol;
            *closestvubidx = i;
         }
         assert(*closestvubidx >= 0);
      }
   }

   return SCIP_OKAY;
}

/** determines the bounds to use for constructing the single-node-flow relaxation of a variable in
 *  the given row.
 */
static
SCIP_RETCODE determineBoundForSNF(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to use for variable bound; NULL for LP solution */
   SCIP_VAR**            vars,               /**< array of problem variables */
   SCIP_Real*            rowcoefs,           /**< (dense) array of variable coefficients in the row */
   int*                  rowinds,            /**< array with positions of non-zero values in the rowcoefs array */
   int                   varposinrow,        /**< position of variable in the rowinds array for which the bounds should be determined */
   int8_t*               binvarused,         /**< array that stores if a binary variable was already used (+1)
                                              *   was not used (0) or was not used but is contained in the row (-1)
                                              */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Real*            bestlb,             /**< pointer to store best lower bound for transformation */
   SCIP_Real*            bestub,             /**< pointer to store best upper bound for transformation */
   SCIP_Real*            bestslb,            /**< pointer to store best simple lower bound for transformation */
   SCIP_Real*            bestsub,            /**< pointer to store best simple upper bound for transformation */
   int*                  bestlbtype,         /**< pointer to store type of best lower bound */
   int*                  bestubtype,         /**< pointer to store type of best upper bound */
   int*                  bestslbtype,        /**< pointer to store type of best simple lower bound */
   int*                  bestsubtype,        /**< pointer to store type of best simple upper bound */
   SCIP_BOUNDTYPE*       selectedbounds,     /**< pointer to store the preferred bound for the transformation */
   SCIP_Bool*            freevariable        /**< pointer to store if variable is a free variable */
   )
{
   SCIP_VAR* var;

   SCIP_Real rowcoef;
   SCIP_Real solval;

   int probidx;

   bestlb[varposinrow] = -SCIPinfinity(scip);
   bestub[varposinrow] = SCIPinfinity(scip);
   bestlbtype[varposinrow] = -3;
   bestubtype[varposinrow] = -3;

   probidx = rowinds[varposinrow];
   var = vars[probidx];
   {
      SCIP_Real QUAD(tmp);
      QUAD_ARRAY_LOAD(tmp, rowcoefs, probidx);
      rowcoef = QUAD_TO_DBL(tmp);
   }

   assert(!EPSZ(rowcoef, QUAD_EPSILON));

   /* get closest simple lower bound and closest simple upper bound */
   SCIP_CALL( findBestLb(scip, var, sol, FALSE, allowlocal, &bestslb[varposinrow], &bestslbtype[varposinrow]) );
   SCIP_CALL( findBestUb(scip, var, sol, FALSE, allowlocal, &bestsub[varposinrow], &bestsubtype[varposinrow]) );

   /* do not use too large bounds */
   if( bestslb[varposinrow] <= -MAXBOUND )
      bestslb[varposinrow] = -SCIPinfinity(scip);

   if( bestsub[varposinrow] >= MAXBOUND )
      bestsub[varposinrow] = SCIPinfinity(scip);

   solval = SCIPgetSolVal(scip, sol, var);

   SCIPdebugMsg(scip, "  %d: %g <%s, idx=%d, lp=%g, [%g(%d),%g(%d)]>:\n", varposinrow, rowcoef, SCIPvarGetName(var), probidx,
      solval, bestslb[varposinrow], bestslbtype[varposinrow], bestsub[varposinrow], bestsubtype[varposinrow]);

   /* mixed integer set cannot be relaxed to 0-1 single node flow set because both simple bounds are -infinity
    * and infinity, respectively
    */
   if( SCIPisInfinity(scip, -bestslb[varposinrow]) && SCIPisInfinity(scip, bestsub[varposinrow]) )
   {
      *freevariable = TRUE;
      return SCIP_OKAY;
   }

   /* get closest lower bound that can be used to define the real variable y'_j in the 0-1 single node flow
    * relaxation
    */
   if( !SCIPisInfinity(scip, bestsub[varposinrow]) )
   {
      bestlb[varposinrow] = bestslb[varposinrow];
      bestlbtype[varposinrow] = bestslbtype[varposinrow];

      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvlb;
         int bestvlbidx;

         SCIP_CALL( getClosestVlb(scip, var, sol, rowcoefs, binvarused, bestsub[varposinrow], rowcoef, &bestvlb, &bestvlbidx) );
         if( SCIPisGT(scip, bestvlb, bestlb[varposinrow]) )
         {
            bestlb[varposinrow] = bestvlb;
            bestlbtype[varposinrow] = bestvlbidx;
         }
      }
   }

   /* get closest upper bound that can be used to define the real variable y'_j in the 0-1 single node flow
    * relaxation
    */
   if( !SCIPisInfinity(scip, -bestslb[varposinrow]) )
   {
      bestub[varposinrow] = bestsub[varposinrow];
      bestubtype[varposinrow] = bestsubtype[varposinrow];

      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvub;
         int bestvubidx;

         SCIP_CALL( getClosestVub(scip, var, sol, rowcoefs, binvarused, bestslb[varposinrow], rowcoef, &bestvub, &bestvubidx) );
         if( SCIPisLT(scip, bestvub, bestub[varposinrow]) )
         {
            bestub[varposinrow] = bestvub;
            bestubtype[varposinrow] = bestvubidx;
         }
      }
   }
   SCIPdebugMsg(scip, "        bestlb=%g(%d), bestub=%g(%d)\n", bestlb[varposinrow], bestlbtype[varposinrow], bestub[varposinrow], bestubtype[varposinrow]);

   /* mixed integer set cannot be relaxed to 0-1 single node flow set because there are no suitable bounds
    * to define the transformed variable y'_j
    */
   if( SCIPisInfinity(scip, -bestlb[varposinrow]) && SCIPisInfinity(scip, bestub[varposinrow]) )
   {
      *freevariable = TRUE;
      return SCIP_OKAY;
   }

   *freevariable = FALSE;

   /* select best upper bound if it is closer to the LP value of y_j and best lower bound otherwise and use this bound
    * to define the real variable y'_j with 0 <= y'_j <= u'_j x_j in the 0-1 single node flow relaxation;
    * prefer variable bounds
    */
   if( SCIPisEQ(scip, solval, (1.0 - boundswitch) * bestlb[varposinrow] + boundswitch * bestub[varposinrow]) && bestlbtype[varposinrow] >= 0 )
   {
      selectedbounds[varposinrow] = SCIP_BOUNDTYPE_LOWER;
   }
   else if( SCIPisEQ(scip, solval, (1.0 - boundswitch) * bestlb[varposinrow] + boundswitch * bestub[varposinrow])
      && bestubtype[varposinrow] >= 0 )
   {
      selectedbounds[varposinrow] = SCIP_BOUNDTYPE_UPPER;
   }
   else if( SCIPisLE(scip, solval, (1.0 - boundswitch) * bestlb[varposinrow] + boundswitch * bestub[varposinrow]) )
   {
      selectedbounds[varposinrow] = SCIP_BOUNDTYPE_LOWER;
   }
   else
   {
      assert(SCIPisGT(scip, solval, (1.0 - boundswitch) * bestlb[varposinrow] + boundswitch * bestub[varposinrow]));
      selectedbounds[varposinrow] = SCIP_BOUNDTYPE_UPPER;
   }

   if( selectedbounds[varposinrow] == SCIP_BOUNDTYPE_LOWER && bestlbtype[varposinrow] >= 0 )
   {
      int vlbvarprobidx;
      SCIP_VAR** vlbvars = SCIPvarGetVlbVars(var);

       /* mark binary variable of vlb so that it is not used for other continuous variables
       * by setting it's position in the aggrrow to a negative value
       */
      vlbvarprobidx = SCIPvarGetProbindex(vlbvars[bestlbtype[varposinrow]]);
      binvarused[vlbvarprobidx] = 1;
   }
   else if ( selectedbounds[varposinrow] == SCIP_BOUNDTYPE_UPPER && bestubtype[varposinrow] >= 0 )
   {
      int vubvarprobidx;
      SCIP_VAR** vubvars = SCIPvarGetVubVars(var);

       /* mark binary variable of vub so that it is not used for other continuous variables
       * by setting it's position in the aggrrow to a negative value
       */
      vubvarprobidx = SCIPvarGetProbindex(vubvars[bestubtype[varposinrow]]);
      binvarused[vubvarprobidx] = 1;
   }

   return SCIP_OKAY;
}

/** construct a 0-1 single node flow relaxation (with some additional simple constraints) of a mixed integer set
 *  corresponding to the given aggrrow a * x <= rhs
 */
static
SCIP_RETCODE constructSNFRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            rowcoefs,           /**< array of coefficients of row */
   QUAD(SCIP_Real        rowrhs),            /**< pointer to right hand side of row */
   int*                  rowinds,            /**< array of variables problem indices for non-zero coefficients in row */
   int                   nnz,                /**< number of non-zeros in row */
   SNF_RELAXATION*       snf,                /**< stores the sign of the transformed variable in summation */
   SCIP_Bool*            success,            /**< stores whether the transformation was valid */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_VAR** vars;
   int i;
   int nnonbinvarsrow;
   int8_t* binvarused;
   int nbinvars;
   SCIP_Real QUAD(transrhs);

   /* arrays to store the selected bound for each non-binary variable in the row */
   SCIP_Real* bestlb;
   SCIP_Real* bestub;
   SCIP_Real* bestslb;
   SCIP_Real* bestsub;
   int* bestlbtype;
   int* bestubtype;
   int* bestslbtype;
   int* bestsubtype;
   SCIP_BOUNDTYPE* selectedbounds;

   *success = FALSE;

   SCIPdebugMsg(scip, "--------------------- construction of SNF relaxation ------------------------------------\n");

   nbinvars = SCIPgetNBinVars(scip);
   vars = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &bestlb, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestub, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestslb, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestsub, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbtype, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubtype, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestslbtype, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestsubtype, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selectedbounds, nnz) );

   /* sort descending to have continuous variables first */
   SCIPsortDownInt(rowinds, nnz);

   /* array to store whether a binary variable is in the row (-1) or has been used (1) due to variable bound usage */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &binvarused, nbinvars) );

   for( i = nnz - 1; i >= 0 && rowinds[i] < nbinvars; --i )
   {
      int j = rowinds[i];
      binvarused[j] = -1;
   }

   nnonbinvarsrow = i + 1;
   /* determine the bounds to use for transforming the non-binary variables */
   for( i = 0; i < nnonbinvarsrow; ++i )
   {
      SCIP_Bool freevariable;

      assert(rowinds[i] >= nbinvars);

      SCIP_CALL( determineBoundForSNF(scip, sol, vars, rowcoefs, rowinds, i, binvarused, allowlocal, boundswitch,
            bestlb, bestub, bestslb, bestsub, bestlbtype, bestubtype, bestslbtype, bestsubtype, selectedbounds, &freevariable) );

      if( freevariable )
      {
         int j;

         /* clear binvarused at indices of binary variables of row */
         for( j = nnz - 1; j >= nnonbinvarsrow; --j )
            binvarused[rowinds[j]] = 0;

         /* clear binvarused at indices of selected variable bounds */
         for( j = 0; j < i; ++j )
         {
            if( selectedbounds[j] == SCIP_BOUNDTYPE_LOWER && bestlbtype[j] >= 0 )
            {
               SCIP_VAR** vlbvars = SCIPvarGetVlbVars(vars[rowinds[j]]);
               binvarused[SCIPvarGetProbindex(vlbvars[bestlbtype[j]])] = 0;
            }
            else if( selectedbounds[j] == SCIP_BOUNDTYPE_UPPER && bestubtype[j] >= 0 )
            {
               SCIP_VAR** vubvars = SCIPvarGetVubVars(vars[rowinds[j]]);
               binvarused[SCIPvarGetProbindex(vubvars[bestubtype[j]])] = 0;
            }
         }

         /* terminate */
         goto TERMINATE;
      }
   }

   *localbdsused = FALSE;
   QUAD_ASSIGN_Q(transrhs, rowrhs);
   snf->ntransvars = 0;

   /* transform non-binary variables */
   for( i = 0; i < nnonbinvarsrow; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real QUAD(rowcoef);
      SCIP_Real solval;
      int probidx;

      probidx = rowinds[i];
      var = vars[probidx];
      QUAD_ARRAY_LOAD(rowcoef, rowcoefs, probidx);
      solval = SCIPgetSolVal(scip, sol, var);

      assert(probidx >= nbinvars);

      if( selectedbounds[i] == SCIP_BOUNDTYPE_LOWER )
      {
         /* use bestlb to define y'_j */

         assert(!SCIPisInfinity(scip, bestsub[i]));
         assert(!SCIPisInfinity(scip, - bestlb[i]));
         assert(bestsubtype[i] == -1 || bestsubtype[i] == -2);
         assert(bestlbtype[i] > -3 && bestlbtype[i] < SCIPvarGetNVlbs(var));

         /* store for y_j that bestlb is the bound used to define y'_j and that y'_j is the associated real variable
          * in the relaxed set
          */
         snf->origcontvars[snf->ntransvars] = probidx;

         if( bestlbtype[i] < 0 )
         {
            SCIP_Real QUAD(val);
            SCIP_Real QUAD(contsolval);
            SCIP_Real QUAD(rowcoeftimesbestsub);

            /* use simple lower bound in bestlb = l_j <= y_j <= u_j = bestsub to define
             *   y'_j = - a_j ( y_j - u_j ) with 0 <= y'_j <=   a_j ( u_j - l_j ) x_j and x_j = 1    if a_j > 0
             *   y'_j =   a_j ( y_j - u_j ) with 0 <= y'_j <= - a_j ( u_j - l_j ) x_j and x_j = 1    if a_j < 0,
             * put j into the set
             *   N2   if a_j > 0
             *   N1   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j u_j
             */
            SCIPquadprecSumDD(val, bestsub[i], -bestlb[i]);
            SCIPquadprecProdQQ(val, val, rowcoef);
            SCIPquadprecSumDD(contsolval, solval, -bestsub[i]);
            SCIPquadprecProdQQ(contsolval, contsolval, rowcoef);

            if( bestlbtype[i] == -2 || bestsubtype[i] == -2 )
               *localbdsused = TRUE;

            SCIPquadprecProdQD(rowcoeftimesbestsub, rowcoef, bestsub[i]);

            /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
            snf->origbinvars[snf->ntransvars] = -1;
            snf->aggrcoefsbin[snf->ntransvars] = 0.0;

            if( QUAD_TO_DBL(rowcoef) > QUAD_EPSILON )
            {
               snf->transvarcoefs[snf->ntransvars] = - 1;
               snf->transvarvubcoefs[snf->ntransvars] = QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = 1.0;
               snf->transcontvarsolvals[snf->ntransvars] = - QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrconstants[snf->ntransvars] = QUAD_TO_DBL(rowcoeftimesbestsub);
               snf->aggrcoefscont[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
            }
            else
            {
               assert(QUAD_TO_DBL(rowcoef) < QUAD_EPSILON);
               snf->transvarcoefs[snf->ntransvars] = 1;
               snf->transvarvubcoefs[snf->ntransvars] = - QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = 1.0;
               snf->transcontvarsolvals[snf->ntransvars] = QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrconstants[snf->ntransvars] = - QUAD_TO_DBL(rowcoeftimesbestsub);
               snf->aggrcoefscont[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
            }
            SCIPquadprecSumQQ(transrhs, transrhs, -rowcoeftimesbestsub);

            SCIPdebugMsg(scip, "    --> bestlb used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=1), rhs=%g-(%g*%g)=%g\n",
                         snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars, snf->transvarvubcoefs[snf->ntransvars],
                         snf->ntransvars, QUAD_TO_DBL(transrhs) + QUAD_TO_DBL(rowcoeftimesbestsub), QUAD_TO_DBL(rowcoef), bestsub, QUAD_TO_DBL(transrhs));
         }
         else
         {
            SCIP_Real QUAD(rowcoefbinary);
            SCIP_Real varsolvalbinary;
            SCIP_Real QUAD(val);
            SCIP_Real QUAD(contsolval);
            SCIP_Real QUAD(rowcoeftimesvlbconst);
            int vlbvarprobidx;

            SCIP_VAR** vlbvars = SCIPvarGetVlbVars(var);
            SCIP_Real* vlbconsts = SCIPvarGetVlbConstants(var);
            SCIP_Real* vlbcoefs = SCIPvarGetVlbCoefs(var);

            /* use variable lower bound in bestlb = l~_j x_j + d_j <= y_j <= u_j = bestsub to define
             *   y'_j = - ( a_j ( y_j - d_j ) + c_j x_j ) with 0 <= y'_j <= - ( a_j l~_j + c_j ) x_j    if a_j > 0
             *   y'_j =     a_j ( y_j - d_j ) + c_j x_j   with 0 <= y'_j <=   ( a_j l~_j + c_j ) x_j    if a_j < 0,
             * where c_j is the coefficient of x_j in the row, put j into the set
             *   N2   if a_j > 0
             *   N1   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j d_j
             */

            vlbvarprobidx = SCIPvarGetProbindex(vlbvars[bestlbtype[i]]);
            assert(binvarused[vlbvarprobidx] == 1);
            assert(vlbvarprobidx < nbinvars);

            QUAD_ARRAY_LOAD(rowcoefbinary, rowcoefs, vlbvarprobidx);
            varsolvalbinary = SCIPgetSolVal(scip, sol, vlbvars[bestlbtype[i]]);

            SCIPquadprecProdQD(val, rowcoef, vlbcoefs[bestlbtype[i]]);
            SCIPquadprecSumQQ(val, val, rowcoefbinary);
            {
               SCIP_Real QUAD(tmp);

               SCIPquadprecProdQD(tmp, rowcoefbinary, varsolvalbinary);
               SCIPquadprecSumDD(contsolval, solval, - vlbconsts[bestlbtype[i]]);
               SCIPquadprecProdQQ(contsolval, contsolval, rowcoef);
               SCIPquadprecSumQQ(contsolval, contsolval, tmp);
            }

            SCIPquadprecProdQD(rowcoeftimesvlbconst, rowcoef, vlbconsts[bestlbtype[i]]);

            /* clear the binvarpos array, since the variable has been processed */
            binvarused[vlbvarprobidx] = 0;

            /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
            snf->origbinvars[snf->ntransvars] = vlbvarprobidx;

            if( QUAD_TO_DBL(rowcoef) > QUAD_EPSILON )
            {
               snf->transvarcoefs[snf->ntransvars] = - 1;
               snf->transvarvubcoefs[snf->ntransvars] = - QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = varsolvalbinary;
               snf->transcontvarsolvals[snf->ntransvars] = - QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefsbin[snf->ntransvars] = - QUAD_TO_DBL(rowcoefbinary);
               snf->aggrcoefscont[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = QUAD_TO_DBL(rowcoeftimesvlbconst);
            }
            else
            {
               assert(QUAD_TO_DBL(rowcoef) < QUAD_EPSILON);
               snf->transvarcoefs[snf->ntransvars] = 1;
               snf->transvarvubcoefs[snf->ntransvars] = QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = varsolvalbinary;
               snf->transcontvarsolvals[snf->ntransvars] = QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefsbin[snf->ntransvars] = QUAD_TO_DBL(rowcoefbinary);
               snf->aggrcoefscont[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = - QUAD_TO_DBL(rowcoeftimesvlbconst);
            }
            SCIPquadprecSumQQ(transrhs, transrhs, -rowcoeftimesvlbconst);

            SCIPdebugMsg(scip, "    --> bestlb used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s), rhs=%g-(%g*%g)=%g\n",
                         snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars, snf->transvarvubcoefs[snf->ntransvars],
                         snf->ntransvars, SCIPvarGetName(vlbvars[bestlbtype[i]]), QUAD_TO_DBL(transrhs) + QUAD_TO_DBL(rowcoeftimesvlbconst), QUAD_TO_DBL(rowcoef),
                         vlbconsts[bestlbtype[i]], snf->transrhs );
         }
      }
      else
      {
         /* use bestub to define y'_j */

         assert(!SCIPisInfinity(scip, bestub[i]));
         assert(!SCIPisInfinity(scip, - bestslb[i]));
         assert(bestslbtype[i] == -1 || bestslbtype[i] == -2);
         assert(bestubtype[i] > -3 && bestubtype[i] < SCIPvarGetNVubs(var));

         /* store for y_j that y'_j is the associated real variable
          * in the relaxed set
          */
         snf->origcontvars[snf->ntransvars] = probidx;

         if( bestubtype[i] < 0 )
         {
            SCIP_Real QUAD(val);
            SCIP_Real QUAD(contsolval);
            SCIP_Real QUAD(rowcoeftimesbestslb);

            /* use simple upper bound in bestslb = l_j <= y_j <= u_j = bestub to define
             *   y'_j =   a_j ( y_j - l_j ) with 0 <= y'_j <=   a_j ( u_j - l_j ) x_j and x_j = 1    if a_j > 0
             *   y'_j = - a_j ( y_j - l_j ) with 0 <= y'_j <= - a_j ( u_j - l_j ) x_j and x_j = 1    if a_j < 0,
             * put j into the set
             *   N1   if a_j > 0
             *   N2   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j l_j
             */
            SCIPquadprecSumDD(val, bestub[i], - bestslb[i]);
            SCIPquadprecProdQQ(val, val, rowcoef);
            SCIPquadprecSumDD(contsolval, solval, - bestslb[i]);
            SCIPquadprecProdQQ(contsolval, contsolval, rowcoef);

            if( bestubtype[i] == -2 || bestslbtype[i] == -2 )
               *localbdsused = TRUE;

            SCIPquadprecProdQD(rowcoeftimesbestslb, rowcoef, bestslb[i]);

            /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
            snf->origbinvars[snf->ntransvars] = -1;
            snf->aggrcoefsbin[snf->ntransvars] = 0.0;

            if( QUAD_TO_DBL(rowcoef) > QUAD_EPSILON )
            {
               snf->transvarcoefs[snf->ntransvars] = 1;
               snf->transvarvubcoefs[snf->ntransvars] = QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = 1.0;
               snf->transcontvarsolvals[snf->ntransvars] = QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefscont[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = - QUAD_TO_DBL(rowcoeftimesbestslb);
            }
            else
            {
               assert(QUAD_TO_DBL(rowcoef) < QUAD_EPSILON);
               snf->transvarcoefs[snf->ntransvars] = - 1;
               snf->transvarvubcoefs[snf->ntransvars] = - QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = 1.0;
               snf->transcontvarsolvals[snf->ntransvars] = - QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefscont[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = QUAD_TO_DBL(rowcoeftimesbestslb);
            }
            SCIPquadprecSumQQ(transrhs, transrhs, -rowcoeftimesbestslb);

            SCIPdebugMsg(scip, "    --> bestub used for trans: ... %s y'_%d + ..., Y'_%d <= %g x_%d (=1), rhs=%g-(%g*%g)=%g\n",
                         snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars, snf->transvarvubcoefs[snf->ntransvars],
                         snf->ntransvars, QUAD_TO_DBL(transrhs) + QUAD_TO_DBL(rowcoeftimesbestslb), QUAD_TO_DBL(rowcoef), bestslb[i], QUAD_TO_DBL(transrhs));
         }
         else
         {
            SCIP_Real QUAD(rowcoefbinary);
            SCIP_Real varsolvalbinary;
            SCIP_Real QUAD(val);
            SCIP_Real QUAD(contsolval);
            SCIP_Real QUAD(rowcoeftimesvubconst);
            int vubvarprobidx;

            SCIP_VAR** vubvars = SCIPvarGetVubVars(var);
            SCIP_Real* vubconsts = SCIPvarGetVubConstants(var);
            SCIP_Real* vubcoefs = SCIPvarGetVubCoefs(var);

            /* use variable upper bound in bestslb = l_j <= y_j <= u~_j x_j + d_j = bestub to define
             *   y'_j =     a_j ( y_j - d_j ) + c_j x_j   with 0 <= y'_j <=   ( a_j u~_j + c_j ) x_j    if a_j > 0
             *   y'_j = - ( a_j ( y_j - d_j ) + c_j x_j ) with 0 <= y'_j <= - ( a_j u~_j + c_j ) x_j    if a_j < 0,
             * where c_j is the coefficient of x_j in the row, put j into the set
             *   N1   if a_j > 0
             *   N2   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j d_j
             */

            vubvarprobidx = SCIPvarGetProbindex(vubvars[bestubtype[i]]);
            assert(binvarused[vubvarprobidx] == 1);
            assert(vubvarprobidx < nbinvars);

            QUAD_ARRAY_LOAD(rowcoefbinary, rowcoefs, vubvarprobidx);
            varsolvalbinary = SCIPgetSolVal(scip, sol, vubvars[bestubtype[i]]);

            /* clear the binvarpos array, since the variable has been processed */
            binvarused[vubvarprobidx] = 0;

            SCIPquadprecProdQD(val, rowcoef, vubcoefs[bestubtype[i]]);
            SCIPquadprecSumQQ(val, val, rowcoefbinary);
            {
               SCIP_Real QUAD(tmp);
               SCIPquadprecProdQD(tmp, rowcoefbinary, varsolvalbinary);
               SCIPquadprecSumDD(contsolval, solval, - vubconsts[bestubtype[i]]);
               SCIPquadprecProdQQ(contsolval, contsolval, rowcoef);
               SCIPquadprecSumQQ(contsolval, contsolval, tmp);
            }

            SCIPquadprecProdQD(rowcoeftimesvubconst, rowcoef, vubconsts[bestubtype[i]]);
            /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
            snf->origbinvars[snf->ntransvars] = vubvarprobidx;

            if( QUAD_TO_DBL(rowcoef) > QUAD_EPSILON )
            {
               snf->transvarcoefs[snf->ntransvars] = 1;
               snf->transvarvubcoefs[snf->ntransvars] = QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = varsolvalbinary;
               snf->transcontvarsolvals[snf->ntransvars] = QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefsbin[snf->ntransvars] = QUAD_TO_DBL(rowcoefbinary);
               snf->aggrcoefscont[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = - QUAD_TO_DBL(rowcoeftimesvubconst);
            }
            else
            {
               assert(QUAD_TO_DBL(rowcoef) < QUAD_EPSILON);
               snf->transvarcoefs[snf->ntransvars] = - 1;
               snf->transvarvubcoefs[snf->ntransvars] = - QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = varsolvalbinary;
               snf->transcontvarsolvals[snf->ntransvars] = - QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefsbin[snf->ntransvars] = - QUAD_TO_DBL(rowcoefbinary);
               snf->aggrcoefscont[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = QUAD_TO_DBL(rowcoeftimesvubconst);
            }
            SCIPquadprecSumQQ(transrhs, transrhs, -rowcoeftimesvubconst);

            /* store for x_j that y'_j is the associated real variable in the 0-1 single node flow relaxation */

            SCIPdebugMsg(scip, "    --> bestub used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s), rhs=%g-(%g*%g)=%g\n",
                         snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars, snf->transvarvubcoefs[snf->ntransvars],
                         snf->ntransvars, SCIPvarGetName(vubvars[bestubtype[i]]), QUAD_TO_DBL(transrhs) + QUAD_TO_DBL(rowcoeftimesvubconst), QUAD_TO_DBL(rowcoef),
                         vubconsts[bestubtype[i]], QUAD_TO_DBL(transrhs));
         }
      }

      /* make sure the coefficient is not negative due to small numerical rounding errors */
      assert(snf->transvarvubcoefs[snf->ntransvars] > -QUAD_EPSILON);
      snf->transvarvubcoefs[snf->ntransvars] = MAX(snf->transvarvubcoefs[snf->ntransvars], 0.0);

      ++snf->ntransvars;
   }

   snf->transrhs = QUAD_TO_DBL(transrhs);

   /* transform remaining binary variables of row */
   for( i = nnonbinvarsrow; i < nnz; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real QUAD(rowcoef);
      int probidx;
      SCIP_Real val;
      SCIP_Real contsolval;
      SCIP_Real varsolval;

      probidx = rowinds[i];
      /* variable should be binary */
      assert(probidx < nbinvars);

      /* binary variable was processed together with a non-binary variable */
      if( binvarused[probidx] == 0 )
         continue;

      /* binary variable was not processed yet, so the binvarused value sould be -1 */
      assert(binvarused[probidx] == -1);

      /* set binvarused to zero since it has been processed */
      binvarused[probidx] = 0;

      var = vars[probidx];
      QUAD_ARRAY_LOAD(rowcoef, rowcoefs, probidx);

      assert(!EPSZ(QUAD_TO_DBL(rowcoef), QUAD_EPSILON));

      varsolval = SCIPgetSolVal(scip, sol, var);
      SCIPdebugMsg(scip, "  %d: %g <%s, idx=%d, lp=%g, [%g, %g]>:\n", i, QUAD_TO_DBL(rowcoef), SCIPvarGetName(var), probidx, varsolval,
         SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));

      /* define
       *    y'_j =   c_j x_j with 0 <= y'_j <=   c_j x_j    if c_j > 0
       *    y'_j = - c_j x_j with 0 <= y'_j <= - c_j x_j    if c_j < 0,
       * where c_j is the coefficient of x_j in the row and put j into the set
       *    N1   if c_j > 0
       *    N2   if c_j < 0.
       */
      val = QUAD_TO_DBL(rowcoef);
      contsolval = QUAD_TO_DBL(rowcoef) * varsolval;

      /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
      snf->origbinvars[snf->ntransvars] = probidx;
      snf->origcontvars[snf->ntransvars] = -1;
      snf->aggrcoefscont[snf->ntransvars] = 0.0;
      snf->aggrconstants[snf->ntransvars] = 0.0;

      if( QUAD_TO_DBL(rowcoef) > QUAD_EPSILON )
      {
         snf->transvarcoefs[snf->ntransvars] = 1;
         snf->transvarvubcoefs[snf->ntransvars] = val;
         snf->transbinvarsolvals[snf->ntransvars] = varsolval;
         snf->transcontvarsolvals[snf->ntransvars] = contsolval;

         /* aggregation information for y'_j */
         snf->aggrcoefsbin[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
      }
      else
      {
         assert(QUAD_TO_DBL(rowcoef) < QUAD_EPSILON);
         snf->transvarcoefs[snf->ntransvars] = - 1;
         snf->transvarvubcoefs[snf->ntransvars] = - val;
         snf->transbinvarsolvals[snf->ntransvars] = varsolval;
         snf->transcontvarsolvals[snf->ntransvars] = - contsolval;

         /* aggregation information for y'_j */
         snf->aggrcoefsbin[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
      }

      assert(snf->transvarcoefs[snf->ntransvars] == 1 || snf->transvarcoefs[snf->ntransvars] == - 1 );
      assert(SCIPisFeasGE(scip, snf->transbinvarsolvals[snf->ntransvars], 0.0)
         && SCIPisFeasLE(scip, snf->transbinvarsolvals[snf->ntransvars], 1.0));
      assert(SCIPisFeasGE(scip, snf->transvarvubcoefs[snf->ntransvars], 0.0)
         && !SCIPisInfinity(scip, snf->transvarvubcoefs[snf->ntransvars]));

      SCIPdebugMsg(scip, "   --> ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s))\n", snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars,
         snf->transvarvubcoefs[snf->ntransvars], snf->ntransvars, SCIPvarGetName(var) );

      /* updates number of variables in transformed problem */
      snf->ntransvars++;
   }

   /* construction was successful */
   *success = TRUE;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "constraint in constructed 0-1 single node flow relaxation: ");
   for( i = 0; i < snf->ntransvars; i++ )
   {
      SCIPdebugMsgPrint(scip, "%s y'_%d ", snf->transvarcoefs[i] == 1 ? "+" : "-", i);
   }
   SCIPdebugMsgPrint(scip, "<= %g\n", snf->transrhs);
#endif

  TERMINATE:

   SCIPfreeCleanBufferArray(scip, &binvarused);
   SCIPfreeBufferArray(scip, &selectedbounds);
   SCIPfreeBufferArray(scip, &bestsubtype);
   SCIPfreeBufferArray(scip, &bestslbtype);
   SCIPfreeBufferArray(scip, &bestubtype);
   SCIPfreeBufferArray(scip, &bestlbtype);
   SCIPfreeBufferArray(scip, &bestsub);
   SCIPfreeBufferArray(scip, &bestslb);
   SCIPfreeBufferArray(scip, &bestub);
   SCIPfreeBufferArray(scip, &bestlb);

   return SCIP_OKAY;
}

/** allocate buffer arrays for storing the single-node-flow relaxation */
static
SCIP_RETCODE allocSNFRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< pointer to snf relaxation to be destroyed */
   int                   nvars               /**< number of active problem variables */
   )
{
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->transvarcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->transbinvarsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->transcontvarsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->transvarvubcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->origbinvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->origcontvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->aggrcoefsbin, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->aggrcoefscont, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->aggrconstants, nvars) );

   return SCIP_OKAY;
}

/** free buffer arrays for storing the single-node-flow relaxation */
static
void destroySNFRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf                 /**< pointer to snf relaxation to be destroyed */
   )
{
   SCIPfreeBufferArray(scip, &snf->aggrconstants);
   SCIPfreeBufferArray(scip, &snf->aggrcoefscont);
   SCIPfreeBufferArray(scip, &snf->aggrcoefsbin);
   SCIPfreeBufferArray(scip, &snf->origcontvars);
   SCIPfreeBufferArray(scip, &snf->origbinvars);
   SCIPfreeBufferArray(scip, &snf->transvarvubcoefs);
   SCIPfreeBufferArray(scip, &snf->transcontvarsolvals);
   SCIPfreeBufferArray(scip, &snf->transbinvarsolvals);
   SCIPfreeBufferArray(scip, &snf->transvarcoefs);
}

/** solve knapsack problem in maximization form with "<" constraint approximately by greedy; if needed, one can provide
 *  arrays to store all selected items and all not selected items
 */
static
SCIP_RETCODE SCIPsolveKnapsackApproximatelyLT(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Real*            weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Real             capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval              /**< pointer to store optimal solution value, or NULL */
   )
{
   SCIP_Real* tempsort;
   SCIP_Real solitemsweight;
   SCIP_Real mediancapacity;
   int j;
   int i;
   int criticalitem;

   assert(weights != NULL);
   assert(profits != NULL);
   assert(SCIPisFeasGE(scip, capacity, 0.0));
   assert(!SCIPisInfinity(scip, capacity));
   assert(items != NULL);
   assert(nitems >= 0);

   if( solitems != NULL )
   {
      *nsolitems = 0;
      *nnonsolitems = 0;
   }
   if( solval != NULL )
      *solval = 0.0;

   /* allocate memory for temporary array used for sorting; array should contain profits divided by corresponding weights (p_1 / w_1 ... p_n / w_n )*/
   SCIP_CALL( SCIPallocBufferArray(scip, &tempsort, nitems) );

   /* initialize temporary array */
   for( i = nitems - 1; i >= 0; --i )
      tempsort[i] = profits[i] / weights[i];

   /* decrease capacity slightly to make it tighter than the original capacity */
   mediancapacity = capacity * (1 - SCIPfeastol(scip));

   /* rearrange items around  */
   SCIPselectWeightedDownRealRealInt(tempsort, profits, items, weights, mediancapacity, nitems, &criticalitem);

   /* free temporary array */
   SCIPfreeBufferArray(scip, &tempsort);

   /* select items as long as they fit into the knapsack */
   solitemsweight = 0.0;
   for( j = 0; j < nitems && SCIPisFeasLT(scip, solitemsweight + weights[j], capacity); j++ )
   {
      if( solitems != NULL )
      {
         solitems[*nsolitems] = items[j];
         (*nsolitems)++;
      }
      if( solval != NULL )
         (*solval) += profits[j];
      solitemsweight += weights[j];
   }


   /* continue to put items into the knapsack if they entirely fit */
   for( ; j < nitems; j++ )
   {
      if( SCIPisFeasLT(scip, solitemsweight + weights[j], capacity) )
      {
         if( solitems != NULL )
         {
            solitems[*nsolitems] = items[j];
            (*nsolitems)++;
         }
         if( solval != NULL )
            (*solval) += profits[j];
         solitemsweight += weights[j];
      }
      else if( solitems != NULL )
      {
         nonsolitems[*nnonsolitems] = items[j];
         (*nnonsolitems)++;
      }
   }

   return SCIP_OKAY;
}


/** build the flow cover which corresponds to the given exact or approximate solution of KP^SNF; given unfinished
 *  flow cover contains variables which have been fixed in advance
 */
static
void buildFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  coefs,              /**< coefficient of all real variables in N1&N2 */
   SCIP_Real*            vubcoefs,           /**< coefficient in vub of all real variables in N1&N2 */
   SCIP_Real             rhs,                /**< right hand side of 0-1 single node flow constraint */
   int*                  solitems,           /**< items in knapsack */
   int*                  nonsolitems,        /**< items not in knapsack */
   int                   nsolitems,          /**< number of items in knapsack */
   int                   nnonsolitems,       /**< number of items not in knapsack */
   int*                  nflowcovervars,     /**< pointer to store number of variables in flow cover */
   int*                  nnonflowcovervars,  /**< pointer to store number of variables not in flow cover */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */
   QUAD(SCIP_Real*       flowcoverweight),   /**< pointer to store weight of flow cover */
   SCIP_Real*            lambda              /**< pointer to store lambda */
   )
{
   int j;
   SCIP_Real QUAD(tmp);

   assert(scip != NULL);
   assert(coefs != NULL);
   assert(vubcoefs != NULL);
   assert(solitems != NULL);
   assert(nonsolitems != NULL);
   assert(nsolitems >= 0);
   assert(nnonsolitems >= 0);
   assert(nflowcovervars != NULL && *nflowcovervars >= 0);
   assert(nnonflowcovervars != NULL && *nnonflowcovervars >= 0);
   assert(flowcoverstatus != NULL);
   assert(QUAD_HI(flowcoverweight) != NULL);
   assert(lambda != NULL);

   /* get flowcover status for each item */
   for( j = 0; j < nsolitems; j++ )
   {
      /* j in N1 with z°_j = 1 => j in N1\C1 */
      if( coefs[solitems[j]] == 1 )
      {
         flowcoverstatus[solitems[j]] = -1;
         (*nnonflowcovervars)++;
      }
      /* j in N2 with z_j = 1 => j in C2 */
      else
      {
         assert(coefs[solitems[j]] == -1);
         flowcoverstatus[solitems[j]] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(*flowcoverweight, *flowcoverweight, -vubcoefs[solitems[j]]);
      }
   }
   for( j = 0; j < nnonsolitems; j++ )
   {
      /* j in N1 with z°_j = 0 => j in C1 */
      if( coefs[nonsolitems[j]] == 1 )
      {
         flowcoverstatus[nonsolitems[j]] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(*flowcoverweight, *flowcoverweight, vubcoefs[nonsolitems[j]]);
      }
      /* j in N2 with z_j = 0 => j in N2\C2 */
      else
      {
         assert(coefs[nonsolitems[j]] == -1);
         flowcoverstatus[nonsolitems[j]] = -1;
         (*nnonflowcovervars)++;
      }
   }

   /* get lambda = sum_{j in C1} u_j - sum_{j in C2} u_j - rhs */
   SCIPquadprecSumQD(tmp, *flowcoverweight, -rhs);
   *lambda = QUAD_TO_DBL(tmp);
}

#ifndef NO_EXACT_KNAPSACK

/** checks, whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = floor(sval);
   upval = ceil(sval);

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}

/** get integral number with error in the bounds which corresponds to given value scaled by a given scalar;
 *  should be used in connection with isIntegralScalar()
 */
static
SCIP_Longint getIntegralVal(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real upval;
   SCIP_Longint intval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   upval = ceil(sval);

   if( SCIPrelDiff(sval, upval) >= mindelta )
      intval = (SCIP_Longint) upval;
   else
      intval = (SCIP_Longint) (floor(sval));

   return intval;
}

/** get a flow cover (C1, C2) for a given 0-1 single node flow set
 *    {(x,y) in {0,1}^n x R^n : sum_{j in N1} y_j - sum_{j in N2} y_j <= b, 0 <= y_j <= u_j x_j},
 *  i.e., get sets C1 subset N1 and C2 subset N2 with sum_{j in C1} u_j - sum_{j in C2} u_j = b + lambda and lambda > 0
 */
static
SCIP_RETCODE getFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< the single node flow relaxation */
   int*                  nflowcovervars,     /**< pointer to store number of variables in flow cover */
   int*                  nnonflowcovervars,  /**< pointer to store number of variables not in flow cover */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */
   SCIP_Real*            lambda,             /**< pointer to store lambda */
   SCIP_Bool*            found               /**< pointer to store whether a cover was found */
   )
{
   SCIP_Real* transprofitsint;
   SCIP_Real* transprofitsreal;
   SCIP_Real* transweightsreal;
   SCIP_Longint* transweightsint;
   int* items;
   int* itemsint;
   int* nonsolitems;
   int* solitems;
   SCIP_Real QUAD(flowcoverweight);
   SCIP_Real QUAD(flowcoverweightafterfix);
   SCIP_Real n1itemsweight;
   SCIP_Real n2itemsminweight;
   SCIP_Real scalar;
   SCIP_Real transcapacityreal;
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
   SCIP_Bool kpexact;
#endif
   SCIP_Bool scalesuccess;
   SCIP_Bool transweightsrealintegral;
   SCIP_Longint transcapacityint;
   int nflowcovervarsafterfix;
   int nitems;
   int nn1items;
   int nnonflowcovervarsafterfix;
   int nnonsolitems;
   int nsolitems;
   int j;

   assert(scip != NULL);
   assert(snf->transvarcoefs != NULL);
   assert(snf->transbinvarsolvals != NULL);
   assert(snf->transvarvubcoefs != NULL);
   assert(snf->ntransvars > 0);
   assert(nflowcovervars != NULL);
   assert(nnonflowcovervars != NULL);
   assert(flowcoverstatus != NULL);
   assert(lambda != NULL);
   assert(found != NULL);

   SCIPdebugMsg(scip, "--------------------- get flow cover ----------------------------------------------------\n");

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &items, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &itemsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofitsreal, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofitsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsreal, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, snf->ntransvars) );

   BMSclearMemoryArray(flowcoverstatus, snf->ntransvars);
   *found = FALSE;
   *nflowcovervars = 0;
   *nnonflowcovervars = 0;

   QUAD_ASSIGN(flowcoverweight, 0.0);
   nflowcovervarsafterfix = 0;
   nnonflowcovervarsafterfix = 0;
   QUAD_ASSIGN(flowcoverweightafterfix, 0.0);
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
   kpexact = FALSE;
#endif

   /* fix some variables in advance according to the following fixing strategy
    *   put j into N1\C1,          if j in N1 and x*_j = 0,
    *   put j into C1,             if j in N1 and x*_j = 1,
    *   put j into C2,             if j in N2 and x*_j = 1,
    *   put j into N2\C2,          if j in N2 and x*_j = 0
    * and get the set of the remaining variables
    */
   SCIPdebugMsg(scip, "0. Fix some variables in advance:\n");
   nitems = 0;
   nn1items = 0;
   n1itemsweight = 0.0;
   n2itemsminweight = SCIP_REAL_MAX;
   for( j = 0; j < snf->ntransvars; j++ )
   {
      assert(snf->transvarcoefs[j] == 1 || snf->transvarcoefs[j] == -1);
      assert(SCIPisFeasGE(scip, snf->transbinvarsolvals[j], 0.0) && SCIPisFeasLE(scip,  snf->transbinvarsolvals[j], 1.0));
      assert(SCIPisFeasGE(scip, snf->transvarvubcoefs[j], 0.0));

      /* if u_j = 0, put j into N1\C1 and N2\C2, respectively */
      if( SCIPisFeasZero(scip, snf->transvarvubcoefs[j]) )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         continue;
      }

      /* x*_j is fractional */
      if( !SCIPisFeasIntegral(scip,  snf->transbinvarsolvals[j]) )
      {
         items[nitems] = j;
         nitems++;
         if( snf->transvarcoefs[j] == 1 )
         {
            n1itemsweight += snf->transvarvubcoefs[j];
            nn1items++;
         }
         else
            n2itemsminweight = MIN(n2itemsminweight, snf->transvarvubcoefs[j]);
      }
      /* j is in N1 and x*_j = 0 */
      else if( snf->transvarcoefs[j] == 1 &&  snf->transbinvarsolvals[j] < 0.5 )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMsg(scip, "     <%d>: in N1-C1\n", j);
      }
      /* j is in N1 and x*_j = 1 */
      else if( snf->transvarcoefs[j] == 1 &&  snf->transbinvarsolvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(flowcoverweight, flowcoverweight, snf->transvarvubcoefs[j]);
         SCIPdebugMsg(scip, "     <%d>: in C1\n", j);
      }
      /* j is in N2 and x*_j = 1 */
      else if( snf->transvarcoefs[j] == -1 &&  snf->transbinvarsolvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(flowcoverweight, flowcoverweight, -snf->transvarvubcoefs[j]);
         SCIPdebugMsg(scip, "     <%d>: in C2\n", j);
      }
      /* j is in N2 and x*_j = 0 */
      else
      {
         assert(snf->transvarcoefs[j] == -1 &&  snf->transbinvarsolvals[j] < 0.5);
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMsg(scip, "     <%d>: in N2-C2\n", j);
      }
   }
   assert((*nflowcovervars) + (*nnonflowcovervars) + nitems == snf->ntransvars);
   assert(nn1items >= 0);

   /* to find a flow cover, transform the following knapsack problem
    *
    * (KP^SNF)      max sum_{j in N1} ( x*_j - 1 ) z_j + sum_{j in N2} x*_j z_j
    *                   sum_{j in N1}          u_j z_j - sum_{j in N2} u_j  z_j > b
    *                                         z_j in {0,1} for all j in N1 & N2
    *
    * 1. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have
    *    positive weights and the constraint is a "<" constraint, by complementing all variables in N1
    *
    *    (KP^SNF_rat)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}          u_j z°_j + sum_{j in N2} u_j  z_j < - b + sum_{j in N1} u_j
    *                                                 z°_j in {0,1} for all j in N1
    *                                                  z_j in {0,1} for all j in N2,
    *    and solve it approximately under consideration of the fixing,
    * or
    * 2. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have
    *    positive integer weights and the constraint is a "<=" constraint, by complementing all variables in N1
    *    and multiplying the constraint by a suitable scalar C
    *
    *    (KP^SNF_int)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}        C u_j z°_j + sum_{j in N2} C u_j  z_j <= c
    *                                                   z°_j in {0,1} for all j in N1
    *                                                    z_j in {0,1} for all j in N2,
    *    where
    *      c = floor[ C (- b + sum_{j in N1} u_j ) ]      if frac[ C (- b + sum_{j in N1} u_j ) ] > 0
    *      c =        C (- b + sum_{j in N1} u_j )   - 1  if frac[ C (- b + sum_{j in N1} u_j ) ] = 0
    *    and solve it exactly under consideration of the fixing.
    */
   SCIPdebugMsg(scip, "1. Transform KP^SNF to KP^SNF_rat:\n");

   /* get weight and profit of variables in KP^SNF_rat and check, whether all weights are already integral */
   transweightsrealintegral = TRUE;
   for( j = 0; j < nitems; j++ )
   {
      transweightsreal[j] = snf->transvarvubcoefs[items[j]];

      if( !isIntegralScalar(transweightsreal[j], 1.0, -MINDELTA, MAXDELTA) )
         transweightsrealintegral = FALSE;

      if( snf->transvarcoefs[items[j]] == 1 )
      {
         transprofitsreal[j] = 1.0 -  snf->transbinvarsolvals[items[j]];
         SCIPdebugMsg(scip, "     <%d>: j in N1:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j],
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
      else
      {
         transprofitsreal[j] =  snf->transbinvarsolvals[items[j]];
         SCIPdebugMsg(scip, "     <%d>: j in N2:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j],
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
   }
   /* get capacity of knapsack constraint in KP^SNF_rat */
   transcapacityreal = - snf->transrhs + QUAD_TO_DBL(flowcoverweight) + n1itemsweight;
   SCIPdebugMsg(scip, "     transcapacity = -rhs(%g) + flowcoverweight(%g) + n1itemsweight(%g) = %g\n",
      snf->transrhs, QUAD_TO_DBL(flowcoverweight), n1itemsweight, transcapacityreal);

   /* there exists no flow cover if the capacity of knapsack constraint in KP^SNF_rat after fixing
    * is less than or equal to zero
    */
   if( SCIPisFeasLE(scip, transcapacityreal/10, 0.0) )
   {
      assert(!(*found));
      goto TERMINATE;
   }

   /* KP^SNF_rat has been solved by fixing some variables in advance */
   assert(nitems >= 0);
   if( nitems == 0)
   {
      /* get lambda = sum_{j in C1} u_j - sum_{j in C2} u_j - rhs */
      SCIPquadprecSumQD(flowcoverweight, flowcoverweight, -snf->transrhs);
      *lambda = QUAD_TO_DBL(flowcoverweight);
      *found = TRUE;
      goto TERMINATE;
   }

   /* Use the following strategy
    *   solve KP^SNF_int exactly,          if a suitable factor C is found and (nitems*capacity) <= MAXDYNPROGSPACE,
    *   solve KP^SNF_rat approximately,    otherwise
    */

   /* find a scaling factor C */
   if( transweightsrealintegral )
   {
      /* weights are already integral */
      scalar = 1.0;
      scalesuccess = TRUE;
   }
   else
   {
      scalesuccess = FALSE;
      SCIP_CALL( SCIPcalcIntegralScalar(transweightsreal, nitems, -MINDELTA, MAXDELTA, MAXDNOM, MAXSCALE, &scalar,
            &scalesuccess) );
   }

   /* initialize number of (non-)solution items, should be changed to a nonnegative number in all possible paths below */
   nsolitems = -1;
   nnonsolitems = -1;

   /* suitable factor C was found*/
   if( scalesuccess )
   {
      SCIP_Real tmp1;
      SCIP_Real tmp2;

      /* transform KP^SNF to KP^SNF_int */
      for( j = 0; j < nitems; ++j )
      {
         transweightsint[j] = getIntegralVal(transweightsreal[j], scalar, -MINDELTA, MAXDELTA);
         transprofitsint[j] = transprofitsreal[j];
         itemsint[j] = items[j];
      }
      if( isIntegralScalar(transcapacityreal, scalar, -MINDELTA, MAXDELTA) )
      {
         transcapacityint = getIntegralVal(transcapacityreal, scalar, -MINDELTA, MAXDELTA);
         transcapacityint -= 1;
      }
      else
         transcapacityint = (SCIP_Longint) (transcapacityreal * scalar);
      nflowcovervarsafterfix = *nflowcovervars;
      nnonflowcovervarsafterfix = *nnonflowcovervars;
      QUAD_ASSIGN_Q(flowcoverweightafterfix, flowcoverweight);

      tmp1 = (SCIP_Real) (nitems + 1);
      tmp2 = (SCIP_Real) ((transcapacityint) + 1);
      if( transcapacityint * nitems <= MAXDYNPROGSPACE && tmp1 * tmp2 <= INT_MAX / 8.0)
      {
         SCIP_Bool success;

         /* solve KP^SNF_int by dynamic programming */
         SCIP_CALL(SCIPsolveKnapsackExactly(scip, nitems, transweightsint, transprofitsint, transcapacityint,
               itemsint, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL, &success));

         if( !success )
         {
            /* solve KP^SNF_rat approximately */
            SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal,
                  transcapacityreal, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));
         }
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
         else
            kpexact = TRUE;
#endif
      }
      else
      {
         /* solve KP^SNF_rat approximately */
         SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
               items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));
         assert(!kpexact);
      }
   }
   else
   {
      /* solve KP^SNF_rat approximately */
      SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
            items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));
      assert(!kpexact);
   }

   assert(nsolitems != -1);
   assert(nnonsolitems != -1);

   /* build the flow cover from the solution of KP^SNF_rat and KP^SNF_int, respectively and the fixing */
   assert(*nflowcovervars + *nnonflowcovervars + nsolitems + nnonsolitems == snf->ntransvars);
   buildFlowCover(scip, snf->transvarcoefs, snf->transvarvubcoefs, snf->transrhs, solitems, nonsolitems, nsolitems, nnonsolitems, nflowcovervars,
      nnonflowcovervars, flowcoverstatus, QUAD(&flowcoverweight), lambda);
   assert(*nflowcovervars + *nnonflowcovervars == snf->ntransvars);

   /* if the found structure is not a flow cover, because of scaling, solve KP^SNF_rat approximately */
   if( SCIPisFeasLE(scip, *lambda, 0.0) )
   {
      assert(kpexact);

      /* solve KP^SNF_rat approximately */
      SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
            items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));
#ifdef SCIP_DEBUG /* this time only for SCIP_DEBUG, because only then, the variable is used again  */
      kpexact = FALSE;
#endif

      /* build the flow cover from the solution of KP^SNF_rat and the fixing */
      *nflowcovervars = nflowcovervarsafterfix;
      *nnonflowcovervars = nnonflowcovervarsafterfix;
      QUAD_ASSIGN_Q(flowcoverweight, flowcoverweightafterfix);

      assert(*nflowcovervars + *nnonflowcovervars + nsolitems + nnonsolitems == snf->ntransvars);
      buildFlowCover(scip, snf->transvarcoefs, snf->transvarvubcoefs, snf->transrhs, solitems, nonsolitems, nsolitems, nnonsolitems, nflowcovervars,
         nnonflowcovervars, flowcoverstatus, QUAD(&flowcoverweight), lambda);
      assert(*nflowcovervars + *nnonflowcovervars == snf->ntransvars);
   }
   *found = SCIPisFeasGT(scip, *lambda, 0.0);

  TERMINATE:
   assert((!*found) || SCIPisFeasGT(scip, *lambda, 0.0));
#ifdef SCIP_DEBUG
   if( *found )
   {
      SCIPdebugMsg(scip, "2. %s solution:\n", kpexact ? "exact" : "approximate");
      for( j = 0; j < snf->ntransvars; j++ )
      {
         if( snf->transvarcoefs[j] == 1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMsg(scip, "     C1: + y_%d [u_%d = %g]\n", j, j, snf->transvarvubcoefs[j]);
         }
         else if( snf->transvarcoefs[j] == -1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMsg(scip, "     C2: - y_%d [u_%d = %g]\n", j, j, snf->transvarvubcoefs[j]);
         }
      }
      SCIPdebugMsg(scip, "     flowcoverweight(%g) = rhs(%g) + lambda(%g)\n", QUAD_TO_DBL(flowcoverweight), snf->transrhs, *lambda);
   }
#endif

   /* free data structures */
   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &transweightsint);
   SCIPfreeBufferArray(scip, &transweightsreal);
   SCIPfreeBufferArray(scip, &transprofitsint);
   SCIPfreeBufferArray(scip, &transprofitsreal);
   SCIPfreeBufferArray(scip, &itemsint);
   SCIPfreeBufferArray(scip, &items);

   return SCIP_OKAY;
}

#else

/** get a flow cover \f$(C1, C2)\f$ for a given 0-1 single node flow set
 *    \f${(x,y) in {0,1}^n x R^n : sum_{j in N1} y_j - sum_{j in N2} y_j <= b, 0 <= y_j <= u_j x_j}\f$,
 *  i.e., get sets \f$ C1 \subset N1 \f$ and \f$ C2 \subset N2 \f$ with
 *  \f$ \sum_{j in C1} u_j - sum_{j in C2} u_j = b + lambda \f$ and \f$ lambda > 0 \f$
 */
static
SCIP_RETCODE getFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< the 0-1 single node flow relaxation */
   int*                  nflowcovervars,     /**< pointer to store number of variables in flow cover */
   int*                  nnonflowcovervars,  /**< pointer to store number of variables not in flow cover */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */
   SCIP_Real*            lambda,             /**< pointer to store lambda */
   SCIP_Bool*            found               /**< pointer to store whether a cover was found */
   )
{
   SCIP_Real* transprofitsreal;
   SCIP_Real* transweightsreal;
   SCIP_Longint* transweightsint;
   int* items;
   int* itemsint;
   int* nonsolitems;
   int* solitems;
   SCIP_Real QUAD(flowcoverweight);
   SCIP_Real n1itemsweight;
   SCIP_Real n2itemsminweight;
   SCIP_Real transcapacityreal;
   int nitems;
   int nn1items;
   int nnonsolitems;
   int nsolitems;
   int j;

   assert(scip != NULL);
   assert(snf->transvarcoefs != NULL);
   assert(snf->transbinvarsolvals != NULL);
   assert(snf->transvarvubcoefs != NULL);
   assert(snf->ntransvars > 0);
   assert(nflowcovervars != NULL);
   assert(nnonflowcovervars != NULL);
   assert(flowcoverstatus != NULL);
   assert(lambda != NULL);
   assert(found != NULL);

   SCIPdebugMsg(scip, "--------------------- get flow cover ----------------------------------------------------\n");

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &items, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &itemsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofitsreal, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsreal, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, snf->ntransvars) );

   BMSclearMemoryArray(flowcoverstatus, snf->ntransvars);
   *found = FALSE;
   *nflowcovervars = 0;
   *nnonflowcovervars = 0;

   QUAD_ASSIGN(flowcoverweight, 0.0);

   /* fix some variables in advance according to the following fixing strategy
    *   put j into N1\C1,          if j in N1 and x*_j = 0,
    *   put j into C1,             if j in N1 and x*_j = 1,
    *   put j into C2,             if j in N2 and x*_j = 1,
    *   put j into N2\C2,          if j in N2 and x*_j = 0
    * and get the set of the remaining variables
    */
   SCIPdebugMsg(scip, "0. Fix some variables in advance:\n");
   nitems = 0;
   nn1items = 0;
   n1itemsweight = 0.0;
   n2itemsminweight = SCIP_REAL_MAX;
   for( j = 0; j < snf->ntransvars; j++ )
   {
      assert(snf->transvarcoefs[j] == 1 || snf->transvarcoefs[j] == -1);
      assert(SCIPisFeasGE(scip, snf->transbinvarsolvals[j], 0.0) && SCIPisFeasLE(scip,  snf->transbinvarsolvals[j], 1.0));
      assert(SCIPisFeasGE(scip, snf->transvarvubcoefs[j], 0.0));

      /* if u_j = 0, put j into N1\C1 and N2\C2, respectively */
      if( SCIPisFeasZero(scip, snf->transvarvubcoefs[j]) )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         continue;
      }

      /* x*_j is fractional */
      if( !SCIPisFeasIntegral(scip,  snf->transbinvarsolvals[j]) )
      {
         items[nitems] = j;
         nitems++;
         if( snf->transvarcoefs[j] == 1 )
         {
            n1itemsweight += snf->transvarvubcoefs[j];
            nn1items++;
         }
         else
            n2itemsminweight = MIN(n2itemsminweight, snf->transvarvubcoefs[j]);
      }
      /* j is in N1 and x*_j = 0 */
      else if( snf->transvarcoefs[j] == 1 &&  snf->transbinvarsolvals[j] < 0.5 )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMsg(scip, "     <%d>: in N1-C1\n", j);
      }
      /* j is in N1 and x*_j = 1 */
      else if( snf->transvarcoefs[j] == 1 &&  snf->transbinvarsolvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(flowcoverweight, flowcoverweight, snf->transvarvubcoefs[j]);
         SCIPdebugMsg(scip, "     <%d>: in C1\n", j);
      }
      /* j is in N2 and x*_j = 1 */
      else if( snf->transvarcoefs[j] == -1 &&  snf->transbinvarsolvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(flowcoverweight, flowcoverweight, -snf->transvarvubcoefs[j]);
         SCIPdebugMsg(scip, "     <%d>: in C2\n", j);
      }
      /* j is in N2 and x*_j = 0 */
      else
      {
         assert(snf->transvarcoefs[j] == -1 &&  snf->transbinvarsolvals[j] < 0.5);
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMsg(scip, "     <%d>: in N2-C2\n", j);
      }
   }
   assert((*nflowcovervars) + (*nnonflowcovervars) + nitems == snf->ntransvars);
   assert(nn1items >= 0);

   /* to find a flow cover, transform the following knapsack problem
    *
    * (KP^SNF)      max sum_{j in N1} ( x*_j - 1 ) z_j + sum_{j in N2} x*_j z_j
    *                   sum_{j in N1}          u_j z_j - sum_{j in N2} u_j  z_j > b
    *                                         z_j in {0,1} for all j in N1 & N2
    *
    * 1. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have
    *    positive weights and the constraint is a "<" constraint, by complementing all variables in N1
    *
    *    (KP^SNF_rat)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}          u_j z°_j + sum_{j in N2} u_j  z_j < - b + sum_{j in N1} u_j
    *                                                 z°_j in {0,1} for all j in N1
    *                                                  z_j in {0,1} for all j in N2,
    *    and solve it approximately under consideration of the fixing,
    * or
    * 2. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have
    *    positive integer weights and the constraint is a "<=" constraint, by complementing all variables in N1
    *    and multiplying the constraint by a suitable scalar C
    *
    *    (KP^SNF_int)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}        C u_j z°_j + sum_{j in N2} C u_j  z_j <= c
    *                                                   z°_j in {0,1} for all j in N1
    *                                                    z_j in {0,1} for all j in N2,
    *    where
    *      c = floor[ C (- b + sum_{j in N1} u_j ) ]      if frac[ C (- b + sum_{j in N1} u_j ) ] > 0
    *      c =        C (- b + sum_{j in N1} u_j )   - 1  if frac[ C (- b + sum_{j in N1} u_j ) ] = 0
    *    and solve it exactly under consideration of the fixing.
    */
   SCIPdebugMsg(scip, "1. Transform KP^SNF to KP^SNF_rat:\n");

   /* get weight and profit of variables in KP^SNF_rat and check, whether all weights are already integral */
   for( j = 0; j < nitems; j++ )
   {
      transweightsreal[j] = snf->transvarvubcoefs[items[j]];

      if( snf->transvarcoefs[items[j]] == 1 )
      {
         transprofitsreal[j] = 1.0 -  snf->transbinvarsolvals[items[j]];
         SCIPdebugMsg(scip, "     <%d>: j in N1:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j],
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
      else
      {
         transprofitsreal[j] =  snf->transbinvarsolvals[items[j]];
         SCIPdebugMsg(scip, "     <%d>: j in N2:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j],
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
   }
   /* get capacity of knapsack constraint in KP^SNF_rat */
   transcapacityreal = - snf->transrhs + QUAD_TO_DBL(flowcoverweight) + n1itemsweight;
   SCIPdebugMsg(scip, "     transcapacity = -rhs(%g) + flowcoverweight(%g) + n1itemsweight(%g) = %g\n",
      snf->transrhs, QUAD_TO_DBL(flowcoverweight), n1itemsweight, transcapacityreal);

   /* there exists no flow cover if the capacity of knapsack constraint in KP^SNF_rat after fixing
    * is less than or equal to zero
    */
   if( SCIPisFeasLE(scip, transcapacityreal/10, 0.0) )
   {
      assert(!(*found));
      goto TERMINATE;
   }

   /* KP^SNF_rat has been solved by fixing some variables in advance */
   assert(nitems >= 0);
   if( nitems == 0 )
   {
      /* get lambda = sum_{j in C1} u_j - sum_{j in C2} u_j - rhs */
      SCIPquadprecSumQD(flowcoverweight, flowcoverweight, -snf->transrhs);
      *lambda = QUAD_TO_DBL(flowcoverweight);
      *found = TRUE;
      goto TERMINATE;
   }

   /* Solve the KP^SNF_rat approximately */

   /* initialize number of (non-)solution items, should be changed to a nonnegative number in all possible paths below */
   nsolitems = -1;
   nnonsolitems = -1;

   /* suitable factor C was found*/
   /* solve KP^SNF_rat approximately */
   SCIP_CALL(SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
                                              items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL));

   assert(nsolitems != -1);
   assert(nnonsolitems != -1);

   /* build the flow cover from the solution of KP^SNF_rat and KP^SNF_int, respectively and the fixing */
   assert(*nflowcovervars + *nnonflowcovervars + nsolitems + nnonsolitems == snf->ntransvars);
   buildFlowCover(scip, snf->transvarcoefs, snf->transvarvubcoefs, snf->transrhs, solitems, nonsolitems, nsolitems, nnonsolitems, nflowcovervars,
      nnonflowcovervars, flowcoverstatus, QUAD(&flowcoverweight), lambda);
   assert(*nflowcovervars + *nnonflowcovervars == snf->ntransvars);

   *found = SCIPisFeasGT(scip, *lambda, 0.0);

  TERMINATE:
   assert((!*found) || SCIPisFeasGT(scip, *lambda, 0.0));
#ifdef SCIP_DEBUG
   if( *found )
   {
      SCIPdebugMsg(scip, "2. approximate solution:\n");
      for( j = 0; j < snf->ntransvars; j++ )
      {
         if( snf->transvarcoefs[j] == 1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMsg(scip, "     C1: + y_%d [u_%d = %g]\n", j, j, snf->transvarvubcoefs[j]);
         }
         else if( snf->transvarcoefs[j] == -1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMsg(scip, "     C2: - y_%d [u_%d = %g]\n", j, j, snf->transvarvubcoefs[j]);
         }
      }
      SCIPdebugMsg(scip, "     flowcoverweight(%g) = rhs(%g) + lambda(%g)\n", QUAD_TO_DBL(flowcoverweight), snf->transrhs, *lambda);
   }
#endif

   /* free data structures */
   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &transweightsint);
   SCIPfreeBufferArray(scip, &transweightsreal);
   SCIPfreeBufferArray(scip, &transprofitsreal);
   SCIPfreeBufferArray(scip, &itemsint);
   SCIPfreeBufferArray(scip, &items);

   return SCIP_OKAY;
}

#endif

/** evaluate the super-additive lifting function for the lifted simple generalized flowcover inequalities
 *  for a given value \f$ x \in \{ u_j \mid j \in C- \} \f$.
 */
static
SCIP_Real evaluateLiftingFunction(
   SCIP*                 scip,               /**< SCIP data structure */
   LIFTINGDATA*          liftingdata,        /**< lifting data to use */
   SCIP_Real             x                   /**< value where to evaluate lifting function */
   )
{
   SCIP_Real QUAD(tmp);
   SCIP_Real xpluslambda;
   int i;

   xpluslambda = x + liftingdata->lambda;

   i = 0;
   while( i < liftingdata->r && SCIPisGT(scip, xpluslambda, liftingdata->M[i+1]) )
      ++i;

   if( i < liftingdata->t )
   {
      if( SCIPisLE(scip, liftingdata->M[i], x) )
      {
         assert(SCIPisLE(scip, xpluslambda, liftingdata->M[i+1]));
         return i * liftingdata->lambda;
      }

      assert(i > 0 && SCIPisLE(scip, liftingdata->M[i], xpluslambda) && x <= liftingdata->M[i]);

      /* return x - liftingdata->M[i] + i * liftingdata->lambda */
      SCIPquadprecProdDD(tmp, i, liftingdata->lambda);
      SCIPquadprecSumQD(tmp, tmp, x);
      SCIPquadprecSumQD(tmp, tmp, -liftingdata->M[i]);
      return QUAD_TO_DBL(tmp);
   }

   if( i < liftingdata->r )
   {
      assert(!SCIPisInfinity(scip, liftingdata->mp));

      /* p = liftingdata->m[i] - (liftingdata->mp - liftingdata->lambda) - liftingdata->ml; */
      SCIPquadprecSumDD(tmp, liftingdata->m[i], -liftingdata->mp);
      SCIPquadprecSumQD(tmp, tmp, -liftingdata->ml);
      SCIPquadprecSumQD(tmp, tmp, liftingdata->lambda);

      /* p = MAX(0.0, p); */
      if( QUAD_HI(tmp) < 0.0 )
      {
         QUAD_ASSIGN(tmp, 0.0);
      }

      SCIPquadprecSumQD(tmp, tmp, liftingdata->M[i]);
      SCIPquadprecSumQD(tmp, tmp, liftingdata->ml);

      if( SCIPisLT(scip, QUAD_TO_DBL(tmp), xpluslambda) )
         return i * liftingdata->lambda;

      assert(SCIPisFeasLE(scip, liftingdata->M[i], xpluslambda) &&
             SCIPisFeasLE(scip, xpluslambda, liftingdata->M[i] + liftingdata->ml +
             MAX(0.0, liftingdata->m[i] - (liftingdata->mp - liftingdata->lambda) - liftingdata->ml)));

      SCIPquadprecProdDD(tmp, i, liftingdata->lambda);
      SCIPquadprecSumQD(tmp, tmp, x);
      SCIPquadprecSumQD(tmp, tmp, - liftingdata->M[i]);
      return QUAD_TO_DBL(tmp);
   }

   assert(i == liftingdata->r && SCIPisLE(scip, liftingdata->M[liftingdata->r], xpluslambda));

   SCIPquadprecProdDD(tmp, liftingdata->r, liftingdata->lambda);
   SCIPquadprecSumQD(tmp, tmp, x);
   SCIPquadprecSumQD(tmp, tmp, - liftingdata->M[liftingdata->r]);
   return QUAD_TO_DBL(tmp);
}

/** computes
 * \f[
 * (\alpha_j, \beta_j) =
 *    \begin{cases}
 *       (0, 0) &\quad\text{if} M_i \leq u_j \leq M_{i+1} - \lambda \\
 *       (1, M_i - i \lambda) &\quad\text{if} M_i − \lambda < u_j < M_i \\
 *    \end{cases}
 * \f]
 */
static
void getAlphaAndBeta(
   SCIP*                 scip,               /**< SCIP data structure */
   LIFTINGDATA*          liftingdata,        /**< pointer to lifting function struct */
   SCIP_Real             vubcoef,            /**< vub coefficient to get alpha and beta for */
   int*                  alpha,              /**< get alpha coefficient for lifting */
   SCIP_Real*            beta                /**< get beta coefficient for lifting */
   )
{
   SCIP_Real vubcoefpluslambda;
   int i;

   vubcoefpluslambda = vubcoef + liftingdata->lambda;

   i = 0;
   while( i < liftingdata->r && SCIPisGT(scip, vubcoefpluslambda, liftingdata->M[i+1]) )
      ++i;

   if( SCIPisLT(scip, vubcoef, liftingdata->M[i]) )
   {
      SCIP_Real QUAD(tmp);
      assert(liftingdata->M[i] < vubcoefpluslambda);
      *alpha = 1;
      SCIPquadprecProdDD(tmp, -i, liftingdata->lambda);
      SCIPquadprecSumQD(tmp, tmp, liftingdata->M[i]);
      *beta = QUAD_TO_DBL(tmp);
   }
   else
   {
      assert(SCIPisSumLE(scip, liftingdata->M[i], vubcoef));
      assert(i == liftingdata->r || SCIPisLE(scip, vubcoefpluslambda, liftingdata->M[i+1]));
      *alpha = 0;
      *beta = 0.0;
   }
}

/** compute relevant data for performing the sequence independent lifting */
static
SCIP_RETCODE computeLiftingData(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< pointer to SNF relaxation */
   int*                  transvarflowcoverstatus, /**< pointer to store whether non-binary var is in L2 (2) or not (-1 or 1) */
   SCIP_Real             lambda,             /**< lambda */
   LIFTINGDATA*          liftingdata,        /**< pointer to lifting function struct */
   SCIP_Bool*            valid               /**< is the lifting data valid */
   )
{
   int i;
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(sumN2mC2LE);
   SCIP_Real QUAD(sumN2mC2GT);
   SCIP_Real QUAD(sumC1LE);
   SCIP_Real QUAD(sumC2);

   SCIP_CALL( SCIPallocBufferArray(scip, &liftingdata->m, snf->ntransvars) );

   liftingdata->r = 0;
   QUAD_ASSIGN(sumN2mC2LE, 0.0);
   QUAD_ASSIGN(sumC1LE, 0.0);
   QUAD_ASSIGN(sumN2mC2GT, 0.0);
   QUAD_ASSIGN(sumC2, 0.0);

   liftingdata->mp = SCIPinfinity(scip);

   *valid = FALSE;

   for( i = 0; i < snf->ntransvars; ++i )
   {
      int s = (snf->transvarcoefs[i] + 1) + (transvarflowcoverstatus[i] + 1)/2;

      switch(s)
      {
         case 0: /* var is in N2 \ C2 */
            assert(snf->transvarvubcoefs[i] >= 0.0);
            assert(snf->transvarcoefs[i] == -1 && transvarflowcoverstatus[i] == -1);

            if( SCIPisGT(scip, snf->transvarvubcoefs[i], lambda) )
            {
               SCIPquadprecSumQD(sumN2mC2GT, sumN2mC2GT, snf->transvarvubcoefs[i]);
               liftingdata->m[liftingdata->r++] = snf->transvarvubcoefs[i];
            }
            else
            {
               SCIPquadprecSumQD(sumN2mC2LE, sumN2mC2LE, snf->transvarvubcoefs[i]);
            }
            break;
         case 1: /* var is in C2 */
            assert(snf->transvarvubcoefs[i] > 0.0);
            assert(snf->transvarcoefs[i] == -1 && transvarflowcoverstatus[i] == 1);

            SCIPquadprecSumQD(sumC2, sumC2, snf->transvarvubcoefs[i]);
            break;
         case 3: /* var is in C1 */
            assert(snf->transvarcoefs[i] == 1 && transvarflowcoverstatus[i] == 1);
            assert(snf->transvarvubcoefs[i] > 0.0);

            if( SCIPisGT(scip, snf->transvarvubcoefs[i], lambda) )
            {
               liftingdata->m[liftingdata->r++] = snf->transvarvubcoefs[i];
               liftingdata->mp = MIN(liftingdata->mp, snf->transvarvubcoefs[i]);
            }
            else
            {
               SCIPquadprecSumQD(sumC1LE, sumC1LE, snf->transvarvubcoefs[i]);
            }
            break;
         default:
            assert(s == 2);
            continue;
      }
   }

   if( SCIPisInfinity(scip, liftingdata->mp) )
   {
      SCIPfreeBufferArray(scip, &liftingdata->m);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &liftingdata->M, liftingdata->r + 1) );

   *valid = TRUE;

   SCIPquadprecSumQQ(tmp, sumC1LE, sumN2mC2LE);
   liftingdata->ml = MIN(lambda, QUAD_TO_DBL(tmp));
   SCIPquadprecSumQD(tmp, sumC2, snf->transrhs);
   liftingdata->d1 = QUAD_TO_DBL(tmp);
   SCIPquadprecSumQQ(tmp, tmp, sumN2mC2GT);
   SCIPquadprecSumQQ(tmp, tmp, sumN2mC2LE);
   liftingdata->d2 = QUAD_TO_DBL(tmp);

   SCIPsortDownReal(liftingdata->m, liftingdata->r);

   /* compute M[i] = sum_{i \in [1,r]} m[i] where m[*] is sorted decreasingly and M[0] = 0 */
   QUAD_ASSIGN(tmp, 0.0);
   for( i = 0; i < liftingdata->r; ++i)
   {
      liftingdata->M[i] = QUAD_TO_DBL(tmp);
      SCIPquadprecSumQD(tmp, tmp, liftingdata->m[i]);
   }

   liftingdata->M[liftingdata->r] = QUAD_TO_DBL(tmp);

   SCIP_UNUSED( SCIPsortedvecFindDownReal(liftingdata->m, liftingdata->mp, liftingdata->r, &liftingdata->t) );
   assert(liftingdata->m[liftingdata->t] == liftingdata->mp || SCIPisInfinity(scip, liftingdata->mp)); /*lint !e777*/

   /* compute t largest index sucht that m_t = mp
    * note that liftingdata->m[t-1] == mp due to zero based indexing of liftingdata->m
    */
   ++liftingdata->t;
   while( liftingdata->t < liftingdata->r && liftingdata->m[liftingdata->t] == liftingdata->mp ) /*lint !e777*/
      ++liftingdata->t;

   liftingdata->lambda = lambda;

   return SCIP_OKAY;
}

/** destroy data used for the sequence independent lifting */
static
void destroyLiftingData(
   SCIP*                 scip,               /**< SCIP data structure */
   LIFTINGDATA*          liftingdata         /**< pointer to lifting function struct */
   )
{
   SCIPfreeBufferArray(scip, &liftingdata->M);
   SCIPfreeBufferArray(scip, &liftingdata->m);
}

/** store the simple lifted flowcover cut defined by the given data in the given arrays
 *  the array for storing the cut coefficients must be all zeros
 */
static
SCIP_RETCODE generateLiftedFlowCoverCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< pointer to SNF relaxation */
   SCIP_AGGRROW*         aggrrow,            /**< aggrrow used to construct SNF relaxation */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */
   SCIP_Real             lambda,             /**< lambda */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   SCIP_Real*            cutrhs,             /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   SCIP_Bool*            success             /**< was the cut successfully generated */
   )
{
   SCIP_Real QUAD(rhs);
   LIFTINGDATA liftingdata;
   int i;

   SCIP_CALL( computeLiftingData(scip, snf, flowcoverstatus, lambda, &liftingdata, success) );
   if( ! *success )
      return SCIP_OKAY;

   QUAD_ASSIGN(rhs, liftingdata.d1);

   *nnz = 0;

   for( i = 0; i < snf->ntransvars; ++i )
   {
      int s = (snf->transvarcoefs[i] + 1) + (flowcoverstatus[i] + 1)/2;

      switch(s)
      {
         case 0: /* var is in N2 \ C2 */
            if( SCIPisGT(scip, snf->transvarvubcoefs[i], lambda) )
            {
               /* var is in L- */
               if( snf->origbinvars[i] != -1 )
               {
                  assert(cutcoefs[snf->origbinvars[i]] == 0.0);
                  cutinds[*nnz] = snf->origbinvars[i];
                  cutcoefs[snf->origbinvars[i]] = -lambda;
                  ++(*nnz);
               }
               else
               {
                  SCIPquadprecSumQD(rhs, rhs, lambda);
               }
            }
            else
            {
               /* var is in L-- */
               if( snf->origcontvars[i] != -1 && snf->aggrcoefscont[i] != 0.0 )
               {
                  assert(cutcoefs[snf->origcontvars[i]] == 0.0);
                  cutinds[*nnz] = snf->origcontvars[i];
                  cutcoefs[snf->origcontvars[i]] = -snf->aggrcoefscont[i];
                  ++(*nnz);
               }

               if( snf->origbinvars[i] != -1 && snf->aggrcoefsbin[i] != 0.0 )
               {
                  assert(cutcoefs[snf->origbinvars[i]] == 0.0);
                  cutinds[*nnz] = snf->origbinvars[i];
                  cutcoefs[snf->origbinvars[i]] = -snf->aggrcoefsbin[i];
                  ++(*nnz);
               }

               SCIPquadprecSumQD(rhs, rhs, snf->aggrconstants[i]);
            }
            break;
         case 1: /* var is in C2 */
         {
            assert(snf->transvarvubcoefs[i] > 0.0);
            assert(snf->transvarcoefs[i] == -1 && flowcoverstatus[i] == 1);

            if( snf->origbinvars[i] != -1 )
            {
               SCIP_Real liftedbincoef = evaluateLiftingFunction(scip, &liftingdata, snf->transvarvubcoefs[i]);
               assert(cutcoefs[snf->origbinvars[i]] == 0.0);
               if( liftedbincoef != 0.0 )
               {
                  cutinds[*nnz] = snf->origbinvars[i];
                  cutcoefs[snf->origbinvars[i]] = -liftedbincoef;
                  ++(*nnz);
                  SCIPquadprecSumQD(rhs, rhs, -liftedbincoef);
               }
            }
            break;
         }
         case 2: /* var is in N1 \ C1 */
         {
            int alpha;
            SCIP_Real beta;

            assert(snf->transvarcoefs[i] == 1 && flowcoverstatus[i] == -1);

            getAlphaAndBeta(scip, &liftingdata, snf->transvarvubcoefs[i], &alpha, &beta);
            assert(alpha == 0 || alpha == 1);

            if( alpha == 1 )
            {
               SCIP_Real QUAD(binvarcoef);
               assert(beta > 0.0);

               if( snf->origcontvars[i] != -1 && snf->aggrcoefscont[i] != 0.0 )
               {
                  assert(cutcoefs[snf->origcontvars[i]] == 0.0);
                  cutinds[*nnz] = snf->origcontvars[i];
                  cutcoefs[snf->origcontvars[i]] = snf->aggrcoefscont[i];
                  ++(*nnz);
               }

               SCIPquadprecSumDD(binvarcoef, snf->aggrcoefsbin[i], -beta);
               if( snf->origbinvars[i] != -1 )
               {
                  SCIP_Real tmp;

                  assert(cutcoefs[snf->origbinvars[i]] == 0.0);

                  tmp = QUAD_TO_DBL(binvarcoef);
                  if( tmp != 0.0 )
                  {
                     cutinds[*nnz] = snf->origbinvars[i];
                     cutcoefs[snf->origbinvars[i]] = tmp;
                     ++(*nnz);
                  }
               }
               else
               {
                  SCIPquadprecSumQQ(rhs, rhs, -binvarcoef);
               }

               SCIPquadprecSumQD(rhs, rhs, -snf->aggrconstants[i]);
            }
            break;
         }
         case 3: /* var is in C1 */
         {
            SCIP_Real bincoef = snf->aggrcoefsbin[i];
            SCIP_Real constant = snf->aggrconstants[i];

            if( snf->origbinvars[i] != -1 && SCIPisGT(scip, snf->transvarvubcoefs[i], lambda) )
            {
               /* var is in C++ */
               SCIP_Real QUAD(tmp);
               SCIP_Real QUAD(tmp2);

               SCIPquadprecSumDD(tmp, snf->transvarvubcoefs[i], -lambda);

               SCIPquadprecSumQD(tmp2, tmp, constant);
               constant = QUAD_TO_DBL(tmp2);

               SCIPquadprecSumQD(tmp2, tmp, -bincoef);
               bincoef = -QUAD_TO_DBL(tmp2);
            }

            if( snf->origbinvars[i] != -1 && bincoef != 0.0 )
            {
               assert(cutcoefs[snf->origbinvars[i]] == 0.0);
               cutinds[*nnz] = snf->origbinvars[i];
               cutcoefs[snf->origbinvars[i]] = bincoef;
               ++(*nnz);
            }

            if( snf->origcontvars[i] != -1 && snf->aggrcoefscont[i] != 0.0 )
            {
               assert(cutcoefs[snf->origcontvars[i]] == 0.0);
               cutinds[*nnz] = snf->origcontvars[i];
               cutcoefs[snf->origcontvars[i]] = snf->aggrcoefscont[i];
               ++(*nnz);
            }

            SCIPquadprecSumQD(rhs, rhs, -constant);
            break;
         }
         default:
            SCIPABORT();
      }
   }

   destroyLiftingData(scip, &liftingdata);

   {
      SCIP_ROW** rows = SCIPgetLPRows(scip);
      for( i = 0; i < aggrrow->nrows; ++i )
      {
         SCIP_ROW* row;
         SCIP_Real rowlhs;
         SCIP_Real rowrhs;
         SCIP_Real slackub;
         SCIP_Real slackcoef;

         slackcoef = aggrrow->rowweights[i] * aggrrow->slacksign[i];
         assert(slackcoef != 0.0);

         /* positive slack was implicitly handled in flow cover separation */
         if( slackcoef > 0.0 )
            continue;

         row = rows[aggrrow->rowsinds[i]];

         /* add the slack's definition multiplied with its coefficient to the cut */
         SCIP_CALL( varVecAddScaledRowCoefs(cutinds, cutcoefs, nnz, row, -aggrrow->rowweights[i]) );

         /* retrieve sides of row */
         rowlhs = row->lhs - row->constant;
         rowrhs = row->rhs - row->constant;

         if( row->integral )
         {
            rowrhs = SCIPfloor(scip, rowrhs);
            rowlhs = SCIPceil(scip, rowlhs);
         }

         slackub = rowrhs - rowlhs;

         /* move slack's constant to the right hand side, and add lambda to the right hand side if the
          * upper bound of the slack is larger than lambda, since then an artifical binary variable
          * for the slack would get coefficient -lambda
          */
         if( aggrrow->slacksign[i] == +1 )
         {
            SCIP_Real rhsslack;
            /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
            assert(!SCIPisInfinity(scip, row->rhs));

            rhsslack = rowrhs - SCIPgetRowMinActivity(scip, row);
            slackub = -aggrrow->rowweights[i] * MIN(rhsslack, slackub);

            if( SCIPisGE(scip, slackub, lambda) )
               SCIPquadprecSumQD(rhs, rhs, lambda);

            SCIPquadprecSumQD(rhs, rhs, -aggrrow->rowweights[i] * rowrhs);
         }
         else
         {
            SCIP_Real lhsslack;
            /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
            assert(!SCIPisInfinity(scip, -row->lhs));

            lhsslack = SCIPgetRowMaxActivity(scip, row) - rowlhs;
            slackub = aggrrow->rowweights[i] * MIN(lhsslack, slackub);

            if( SCIPisGE(scip, slackub, lambda) )
               SCIPquadprecSumQD(rhs, rhs, lambda);

            SCIPquadprecSumQD(rhs, rhs, -aggrrow->rowweights[i] * rowlhs);
         }
      }
   }

   *cutrhs = QUAD_TO_DBL(rhs);

   /* relax rhs to zero, if it's very close to */
   if( *cutrhs < 0.0 && *cutrhs >= SCIPepsilon(scip) )
      *cutrhs = 0.0;

   return SCIP_OKAY;
}

/** calculates a lifted simple generalized flow cover cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in the cut.
 *  For further details we refer to:
 *
 *  Gu, Z., Nemhauser, G. L., & Savelsbergh, M. W. (1999). Lifted flow cover inequalities for mixed 0-1 integer programs.
 *  Mathematical Programming, 85(3), 439-467.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcalcFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute flow cover cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store the efficacy of the cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid cut was returned */
   )
{
   int i;
   int nvars;
   SCIP_Bool localbdsused;
   SNF_RELAXATION snf;
   SCIP_Real lambda;
   SCIP_Real* tmpcoefs;
   int *transvarflowcoverstatus;
   int nflowcovervars;
   int nnonflowcovervars;

   nvars = SCIPgetNVars(scip);

   *success = FALSE;

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &transvarflowcoverstatus, nvars) );
   SCIP_CALL( allocSNFRelaxation(scip,  &snf, nvars) );

   SCIPdebug( printCutQuad(scip, sol, aggrrow->vals, QUAD(aggrrow->rhs), aggrrow->inds, aggrrow->nnz, FALSE, aggrrow->local) );

   SCIP_CALL( constructSNFRelaxation(scip, sol, boundswitch, allowlocal, aggrrow->vals, QUAD(aggrrow->rhs), aggrrow->inds, aggrrow->nnz, &snf, success, &localbdsused) );

   if( ! *success )
   {
      goto TERMINATE;
   }

   *cutislocal = aggrrow->local || localbdsused;

   /* initialize lambda because gcc issues a stupid warning */
   lambda = 0.0;
   SCIP_CALL( getFlowCover(scip, &snf, &nflowcovervars, &nnonflowcovervars, transvarflowcoverstatus, &lambda, success) );

   if( ! *success )
   {
      goto TERMINATE;
   }

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &tmpcoefs, nvars) );

   SCIP_CALL( generateLiftedFlowCoverCut(scip, &snf, aggrrow, transvarflowcoverstatus, lambda, tmpcoefs, cutrhs, cutinds, cutnnz, success) );
   SCIPdebugMsg(scip, "computed flowcover_%lli_%i:\n", SCIPgetNLPs(scip), SCIPgetNCuts(scip));

   /* if success is FALSE generateLiftedFlowCoverCut wont have touched the tmpcoefs array so we dont need to clean it then */
   if( *success )
   {
      if( postprocess )
      {
         SCIP_CALL( postprocessCut(scip, *cutislocal, cutinds, tmpcoefs, cutnnz, cutrhs, success) );
      }
      else
      {
         SCIP_Real QUAD(rhs);

         QUAD_ASSIGN(rhs, *cutrhs);
         *success = ! removeZeros(scip, SCIPsumepsilon(scip), *cutislocal, tmpcoefs, QUAD(&rhs), cutinds, cutnnz);
         *cutrhs = QUAD_TO_DBL(rhs);
      }

      if( *success )
      {
         /* store cut sparse and calculate efficacy */
         for( i = 0; i < *cutnnz; ++i )
         {
            int j = cutinds[i];
            assert(tmpcoefs[j] != 0.0);
            cutcoefs[i] = tmpcoefs[j];
            tmpcoefs[j] = 0.0;
         }

         if( cutefficacy != NULL )
            *cutefficacy = calcEfficacy(scip, sol, cutcoefs, *cutrhs, cutinds, *cutnnz);

         if( cutrank != NULL )
            *cutrank = aggrrow->rank + 1;
      }
      else
      {
         /* clean buffer array */
         for( i = 0; i < *cutnnz; ++i )
         {
            int j = cutinds[i];
            assert(tmpcoefs[j] != 0.0);
            tmpcoefs[j] = 0.0;
         }
      }
   }

   SCIPfreeCleanBufferArray(scip, &tmpcoefs);

  TERMINATE:
   destroySNFRelaxation(scip, &snf);
   SCIPfreeBufferArray(scip, &transvarflowcoverstatus);

   return SCIP_OKAY;
}


/* =========================================== strongcg =========================================== */

/** Transform equation \f$ a \cdot x = b; lb \leq x \leq ub \f$ into standard form
 *    \f$ a^\prime \cdot x^\prime = b,\; 0 \leq x^\prime \leq ub' \f$.
 *
 *  Transform variables (lb or ub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \mbox{if lb is used in transformation}\\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if ub is used in transformation}
 *  \end{array}
 *  \f]
 *  and move the constant terms \f$ a_j\, lb_j \f$ or \f$ a_j\, ub_j \f$ to the rhs.
 *
 *  Transform variables (vlb or vub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - (bl_j\, zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \mbox{if vlb is used in transf.} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if vub is used in transf.}
 *  \end{array}
 *  \f]
 *  move the constant terms \f$ a_j\, dl_j \f$ or \f$ a_j\, du_j \f$ to the rhs, and update the coefficient of the VLB variable:
 *  \f[
 *  \begin{array}{ll}
 *    a_{zl_j} := a_{zl_j} + a_j\, bl_j,& \mbox{or} \\
 *    a_{zu_j} := a_{zu_j} + a_j\, bu_j &
 *  \end{array}
 *  \f]
 */
static
SCIP_RETCODE cutsTransformStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in MIR row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_Real* bestbds;
   int i;
   int aggrrowintstart;
   int nvars;
   int firstcontvar;
   SCIP_VAR** vars;

   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);

   *freevariable = FALSE;
   *localbdsused = FALSE;

   /* allocate temporary memory to store best bounds and bound types */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestbds, 2*(*nnz)) );

   /* start with continuous variables, because using variable bounds can affect the untransformed integral
    * variables, and these changes have to be incorporated in the transformation of the integral variables
    * (continuous variables have largest problem indices!)
    */
   SCIPsortDownInt(cutinds, *nnz);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   firstcontvar = nvars - SCIPgetNContVars(scip);

   /* determine best bounds for the continous variables such that they will have a positive coefficient in the transformation */
   for( i = 0; i < *nnz && cutinds[i] >= firstcontvar; ++i )
   {
      SCIP_Real QUAD(coef);
      int v = cutinds[i];

      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      if( QUAD_TO_DBL(coef) > 0.0 )
      {
         /* find closest lower bound in standard lower bound or variable lower bound for continuous variable so that it will have a positive coefficient */
         SCIP_CALL( findBestLb(scip, vars[v], sol, usevbds, allowlocal, bestbds + i, boundtype + i) );

         /* cannot create transformation for strongcg cut */
         if( SCIPisInfinity(scip, -bestbds[i]) )
         {
            *freevariable = TRUE;
            goto TERMINATE;
         }

         varsign[i] = +1;
      }
      else if( QUAD_TO_DBL(coef) < 0.0 )
      {
         /* find closest upper bound in standard upper bound or variable upper bound for continuous variable so that it will have a positive coefficient */
         SCIP_CALL( findBestUb(scip, vars[cutinds[i]], sol, usevbds, allowlocal, bestbds + i, boundtype + i) );

          /* cannot create transformation for strongcg cut */
         if( SCIPisInfinity(scip, bestbds[i]) )
         {
            *freevariable = TRUE;
            goto TERMINATE;
         }

         varsign[i] = -1;
      }
   }

   /* remember start of integer variables in the aggrrow */
   aggrrowintstart = i;

   /* perform bound substitution for continuous variables */
   for( i = 0; i < aggrrowintstart; ++i )
   {
      SCIP_Real QUAD(coef);
      SCIP_Real QUAD(tmp);
      int v = cutinds[i];
      SCIP_VAR* var = vars[v];
      assert(!SCIPisInfinity(scip, -varsign[i] * bestbds[i]));

      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
      if( boundtype[i] < 0 )
      {
         SCIPquadprecProdQD(tmp, coef, bestbds[i]);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
         *localbdsused = *localbdsused || (boundtype[i] == -2);
      }
      else
      {
         SCIP_VAR** vbdvars;
         SCIP_Real* vbdcoefs;
         SCIP_Real* vbdconsts;
         SCIP_Real QUAD(zcoef);
         int zidx;

         if( varsign[i] == +1 )
         {
            vbdvars = SCIPvarGetVlbVars(var);
            vbdcoefs = SCIPvarGetVlbCoefs(var);
            vbdconsts = SCIPvarGetVlbConstants(var);
            assert(0 <= boundtype[i] && boundtype[i] < SCIPvarGetNVlbs(var));
         }
         else
         {
            vbdvars = SCIPvarGetVubVars(var);
            vbdcoefs = SCIPvarGetVubCoefs(var);
            vbdconsts = SCIPvarGetVubConstants(var);
            assert(0 <= boundtype[i] && boundtype[i] < SCIPvarGetNVubs(var));
         }

         assert(vbdvars != NULL);
         assert(vbdcoefs != NULL);
         assert(vbdconsts != NULL);
         assert(SCIPvarIsActive(vbdvars[boundtype[i]]));

         zidx = SCIPvarGetProbindex(vbdvars[boundtype[i]]);
         assert(0 <= zidx && zidx < firstcontvar);

         SCIPquadprecProdQD(tmp, coef, vbdconsts[boundtype[i]]);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);

         /* check if integral variable already exists in the row */
         QUAD_ARRAY_LOAD(zcoef, cutcoefs, zidx);

         if( QUAD_HI(zcoef) == 0.0 )
            cutinds[(*nnz)++] = zidx;

         SCIPquadprecProdQD(tmp, coef, vbdcoefs[boundtype[i]]);
         SCIPquadprecSumQQ(zcoef, zcoef, tmp);

         QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
         assert(QUAD_HI(zcoef) != 0.0);

         QUAD_ARRAY_STORE(cutcoefs, zidx, zcoef);
      }
   }

   assert(i == aggrrowintstart);

   /* remove integral variables that now have a zero coefficient due to variable bound usage of continuous variables
    * and perform the bound substitution for the integer variables that are left using simple bounds
    */
   while( i < *nnz )
   {
      SCIP_Real QUAD(coef);
      SCIP_Real QUAD(tmp);
      SCIP_Real bestlb;
      SCIP_Real bestub;
      int bestlbtype;
      int bestubtype;
      SCIP_BOUNDTYPE selectedbound;
      int v = cutinds[i];

      assert(v < firstcontvar);
      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      /* due to variable bound usage for the continous variables cancellation may have occurred */
      if( EPSZ(QUAD_TO_DBL(coef), QUAD_EPSILON) )
      {
         QUAD_ASSIGN(coef, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, coef);
         --(*nnz);
         cutinds[i] = cutinds[*nnz];

         /* do not increase i, since last element is copied to the i-th position */
         continue;
      }

      /* determine the best bounds for the integral variable, usevbd can be set to FALSE here as vbds are only used for continous variables */
      SCIP_CALL( determineBestBounds(scip, vars[v], sol, boundswitch, FALSE, allowlocal, FALSE, FALSE, NULL, NULL,
                                     &bestlb, &bestub, &bestlbtype, &bestubtype, &selectedbound, freevariable) );

      /* check if we have an unbounded integral variable */
      if( *freevariable )
      {
         goto TERMINATE;
      }

      /* perform bound substitution */
      if( selectedbound == SCIP_BOUNDTYPE_LOWER )
      {
         boundtype[i] = bestlbtype;
         varsign[i] = +1;
         SCIPquadprecProdQD(tmp, coef, bestlb);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
      }
      else
      {
         assert(selectedbound == SCIP_BOUNDTYPE_UPPER);
         boundtype[i] = bestubtype;
         varsign[i] = -1;
         SCIPquadprecProdQD(tmp, coef, bestub);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
      }

      assert(boundtype[i] == -1 || boundtype[i] == -2);
      *localbdsused = *localbdsused || (boundtype[i] == -2);

      /* increase i */
      ++i;
   }

   /* relax rhs to zero if it is close to */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

  TERMINATE:
   /*free temporary memory */
   SCIPfreeBufferArray(scip, &bestbds);

   return SCIP_OKAY;
}

/** Calculate fractionalities \f$ f_0 := b - down(b) \f$, \f$ f_j := a^\prime_j - down(a^\prime_j) \f$ and
 *   integer \f$ k >= 1 \f$ with \f$ 1/(k + 1) <= f_0 < 1/k \f$ and \f$ (=> k = up(1/f_0) + 1) \f$
 *   integer \f$ 1 <= p_j <= k \f$ with \f$ f_0 + ((p_j - 1) * (1 - f_0)/k) < f_j <= f_0 + (p_j * (1 - f_0)/k)\f$ \f$ (=> p_j = up( k*(f_j - f_0)/(1 - f_0) )) \f$
 * and derive strong CG cut \f$ \tilde{a}*x^\prime <= down(b) \f$
 * \f[
 * \begin{array}{rll}
 * integers : &  \tilde{a}_j = down(a^\prime_j)                &, if \qquad f_j <= f_0 \\
 *            &  \tilde{a}_j = down(a^\prime_j) + p_j/(k + 1)  &, if \qquad f_j >  f_0 \\
 * continuous:&  \tilde{a}_j = 0                               &, if \qquad a^\prime_j >= 0 \\
 *            &  \mbox{no strong CG cut found}                 &, if \qquad a^\prime_j <  0
 * \end{array}
 * \f]
 *
 * Transform inequality back to \f$ \hat{a}*x <= rhs \f$:
 *
 *  (lb or ub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - lb_j,&   x_j == x^\prime_j + lb_j,&   a^\prime_j ==  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{if lb was used in transformation} \\
 *    x^\prime_j := ub_j - x_j,&   x_j == ub_j - x^\prime_j,&   a^\prime_j == -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{if ub was used in transformation}
 * \end{array}
 * \f]
 * \f[
 *  and move the constant terms
 * \begin{array}{rl}
 *    -\tilde{a}_j * lb_j == -\hat{a}_j * lb_j, & \mbox{or} \\
 *     \tilde{a}_j * ub_j == -\hat{a}_j * ub_j &
 * \end{array}
 * \f]
 *  to the rhs.
 *
 *  (vlb or vub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - (bl_j * zl_j + dl_j),&   x_j == x^\prime_j + (bl_j * zl_j + dl_j),&   a^\prime_j ==  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{(vlb)} \\
 *    x^\prime_j := (bu_j * zu_j + du_j) - x_j,&   x_j == (bu_j * zu_j + du_j) - x^\prime_j,&   a^\prime_j == -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{(vub)}
 * \end{array}
 * \f]
 *  move the constant terms
 * \f[
 * \begin{array}{rl}
 *    -\tilde{a}_j * dl_j == -\hat{a}_j * dl_j,& \mbox{or} \\
 *     \tilde{a}_j * du_j == -\hat{a}_j * du_j &
 * \end{array}
 * \f]
 *  to the rhs, and update the VB variable coefficients:
 * \f[
 * \begin{array}{ll}
 *    \hat{a}_{zl_j} := \hat{a}_{zl_j} - \tilde{a}_j * bl_j == \hat{a}_{zl_j} - \hat{a}_j * bl_j,& \mbox{or} \\
 *    \hat{a}_{zu_j} := \hat{a}_{zu_j} + \tilde{a}_j * bu_j == \hat{a}_{zu_j} - \hat{a}_j * bu_j &
 * \end{array}
 * \f]
 */
static
SCIP_RETCODE cutsRoundStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub)*/
   QUAD(SCIP_Real        f0),                /**< fractional value of rhs */
   SCIP_Real             k                   /**< factor to strengthen strongcg cut */
   )
{
   SCIP_Real QUAD(onedivoneminusf0);
   int i;
   int firstcontvar;
   SCIP_VAR** vars;
   int aggrrowintstart;

   assert(QUAD_HI(cutrhs) != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(boundtype != NULL);
   assert(varsign != NULL);
   assert(0.0 < QUAD_TO_DBL(f0) && QUAD_TO_DBL(f0) < 1.0);

   SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0);
   SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0);

   /* Loop backwards to process integral variables first and be able to delete coefficients of integral variables
    * without destroying the ordering of the aggrrow's non-zeros.
    * (due to sorting in cutsTransformStrongCG the ordering is continuous before integral)
    */

   firstcontvar = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   vars = SCIPgetVars(scip);
#ifndef NDEBUG
   /*in debug mode check, that all continuous variables of the aggrrow come before the integral variables */
   i = 0;
   while( i < *nnz && cutinds[i] >= firstcontvar )
      ++i;

   while( i < *nnz )
   {
      assert(cutinds[i] < firstcontvar);
      ++i;
   }
#endif

   /* integer variables */
   for( i = *nnz - 1; i >= 0 && cutinds[i] < firstcontvar; --i )
   {
      SCIP_VAR* var;
      SCIP_Real QUAD(aj);
      SCIP_Real QUAD(downaj);
      SCIP_Real QUAD(cutaj);
      SCIP_Real QUAD(fj);
      int v;

      v = cutinds[i];
      assert(0 <= v && v < SCIPgetNVars(scip));

      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) == v);
      assert(boundtype[i] == -1 || boundtype[i] == -2);
      assert(varsign[i] == +1 || varsign[i] == -1);

      /* calculate the coefficient in the retransformed cut */
      QUAD_ARRAY_LOAD(aj, cutcoefs, v);
      QUAD_SCALE(aj, varsign[i]);

      SCIPquadprecEpsFloorQ(downaj, aj, SCIPepsilon(scip)); /*lint !e666*/
      SCIPquadprecSumQQ(fj, aj, -downaj);

      if( SCIPisLE(scip, QUAD_TO_DBL(fj), QUAD_TO_DBL(f0)) )
         QUAD_ASSIGN_Q(cutaj, downaj); /* a^_j */
      else
      {
         SCIP_Real pj;

         SCIPquadprecSumQQ(cutaj, fj, -f0);
         SCIPquadprecProdQD(cutaj, cutaj, k);
         SCIPquadprecProdQQ(cutaj, cutaj, onedivoneminusf0);
         pj = SCIPceil(scip, QUAD_TO_DBL(cutaj));
         assert(pj >= 0); /* should be >= 1, but due to rounding bias can be 0 if fj almost equal to f0 */
         assert(pj <= k);
         SCIPquadprecDivDD(cutaj, pj, k + 1.0);
         SCIPquadprecSumQQ(cutaj, cutaj, downaj);
      }

      QUAD_SCALE(cutaj, varsign[i]);

      /* remove zero cut coefficients from cut */
      if( EPSZ(QUAD_TO_DBL(cutaj), QUAD_EPSILON) )
      {
         QUAD_ASSIGN(cutaj, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, cutaj);
         --*nnz;
         cutinds[i] = cutinds[*nnz];
         continue;
      }

      QUAD_ARRAY_STORE(cutcoefs, v, cutaj);

       /* integral var uses standard bound */
      assert(boundtype[i] < 0);

      /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
      if( varsign[i] == +1 )
      {
         SCIP_Real QUAD(tmp);

         /* lower bound was used */
         if( boundtype[i] == -1 )
         {
            assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
            SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbGlobal(var));
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
         }
         else
         {
            assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
            SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbLocal(var));
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
         }
      }
      else
      {
         SCIP_Real QUAD(tmp);

         /* upper bound was used */
         if( boundtype[i] == -1 )
         {
            assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
            SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbGlobal(var));
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
         }
         else
         {
            assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
            SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbLocal(var));
            SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
         }
      }
   }

   /* now process the continuous variables; postpone deletetion of zeros till all continuous variables have been processed */
   aggrrowintstart = i + 1;

#ifndef NDEBUG
   /* in a strong CG cut, cut coefficients of continuous variables are always zero; check this in debug mode */
   for( i = 0; i < aggrrowintstart; ++i )
   {
      int v;

      v = cutinds[i];
      assert(firstcontvar <= v && v < SCIPgetNVars(scip));

      {
         SCIP_VAR* var;
         SCIP_Real QUAD(aj);

         var = vars[v];
         assert(var != NULL);
         assert(!SCIPvarIsIntegral(var));
         assert(SCIPvarGetProbindex(var) == v);
         assert(varsign[i] == +1 || varsign[i] == -1);

         /* calculate the coefficient in the retransformed cut */
         QUAD_ARRAY_LOAD(aj, cutcoefs, v);
         QUAD_SCALE(aj, varsign[i]);

         assert(QUAD_TO_DBL(aj) >= 0.0);
      }
   }
#endif

   /* move integer variables to the empty position of the continuous variables */
   if( aggrrowintstart > 0 )
   {
      SCIP_Real QUAD(tmp);
      assert(aggrrowintstart <= *nnz);

      QUAD_ASSIGN(tmp, 0.0);

      for( i = 0; i < aggrrowintstart; ++i )
      {
         QUAD_ARRAY_STORE(cutcoefs, cutinds[i], tmp);
      }

      *nnz -= aggrrowintstart;
      if( *nnz < aggrrowintstart )
      {
         BMScopyMemoryArray(cutinds, cutinds + aggrrowintstart, *nnz);
      }
      else
      {
         BMScopyMemoryArray(cutinds, cutinds + *nnz, aggrrowintstart);
      }
   }

   return SCIP_OKAY;
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row: \f$ a^\prime_r = scale * weight[r] * slacksign[r] \f$.
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 * \f[
 * \begin{array}{rll}
 *    integers:  & \hat{a}_r = \tilde{a}_r = down(a^\prime_r)                  &, if \qquad f_r <= f0 \\
 *               & \hat{a}_r = \tilde{a}_r = down(a^\prime_r) + p_r/(k + 1)    &, if \qquad f_r >  f0 \\
 *    continuous:& \hat{a}_r = \tilde{a}_r = 0                                 &, if \qquad a^\prime_r >= 0 \\
 *               & \mbox{no strong CG cut found}                               &, if \qquad a^\prime_r <  0
 * \end{array}
 * \f]
 *
 *  Substitute \f$ \hat{a}_r * s_r \f$ by adding \f$ \hat{a}_r \f$ times the slack's definition to the cut.
 */
static
SCIP_RETCODE cutsSubstituteStrongCG(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   int*                  rowinds,            /**< sparsity pattern of used rows */
   int                   nrowinds,           /**< number of used rows */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   QUAD(SCIP_Real        f0),                /**< fractional value of rhs */
   SCIP_Real             k                   /**< factor to strengthen strongcg cut */
   )
{  /*lint --e{715}*/
   SCIP_ROW** rows;
   SCIP_Real QUAD(onedivoneminusf0);
   int i;

   assert(scip != NULL);
   assert(weights != NULL);
   assert(slacksign != NULL);
   assert(rowinds != NULL);
   assert(SCIPisPositive(scip, scale));
   assert(cutcoefs != NULL);
   assert(QUAD_HI(cutrhs) != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(0.0 < QUAD_TO_DBL(f0) && QUAD_TO_DBL(f0) < 1.0);

   SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0);
   SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0);

   rows = SCIPgetLPRows(scip);
   for( i = 0; i < nrowinds; i++ )
   {
      SCIP_ROW* row;
      SCIP_Real pr;
      SCIP_Real QUAD(ar);
      SCIP_Real downar;
      SCIP_Real QUAD(cutar);
      SCIP_Real QUAD(fr);
      SCIP_Real mul;
      int r;

      r = rowinds[i];
      assert(0 <= r && r < SCIPgetNLPRows(scip));
      assert(slacksign[i] == -1 || slacksign[i] == +1);
      assert(!SCIPisZero(scip, weights[i]));

      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* get the slack's coefficient a'_r in the aggregated row */
      SCIPquadprecProdDD(ar, slacksign[i] * scale, weights[i]);

      /* calculate slack variable's coefficient a^_r in the cut */
      if( row->integral )
      {
         /* slack variable is always integral: */
         downar = EPSFLOOR(QUAD_TO_DBL(ar), QUAD_EPSILON);
         SCIPquadprecSumQD(fr, ar, -downar);

         if( SCIPisLE(scip, QUAD_TO_DBL(fr), QUAD_TO_DBL(f0)) )
            QUAD_ASSIGN(cutar, downar);
         else
         {
            SCIPquadprecSumQQ(cutar, fr, -f0);
            SCIPquadprecProdQQ(cutar, cutar, onedivoneminusf0);
            SCIPquadprecProdQD(cutar, cutar, k);
            pr = SCIPceil(scip, QUAD_TO_DBL(cutar));
            assert(pr >= 0); /* should be >= 1, but due to rounding bias can be 0 if fr almost equal to f0 */
            assert(pr <= k);
            SCIPquadprecDivDD(cutar, pr, k + 1.0);
            SCIPquadprecSumQD(cutar, cutar, downar);
         }
      }
      else
      {
         /* slack variable is continuous: */
         assert(QUAD_TO_DBL(ar) >= 0.0);
         continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( EPSZ(QUAD_TO_DBL(cutar), QUAD_EPSILON) )
         continue;

      /* depending on the slack's sign, we have
       *   a*x + c + s == rhs  =>  s == - a*x - c + rhs,  or  a*x + c - s == lhs  =>  s == a*x + c - lhs
       * substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
       */
      mul = -slacksign[i] * QUAD_TO_DBL(cutar);

      /* add the slack's definition multiplied with a^_j to the cut */
      SCIP_CALL( varVecAddScaledRowCoefsQuad(cutinds, cutcoefs, nnz, row, mul) );

      /* move slack's constant to the right hand side */
      if( slacksign[i] == +1 )
      {
         SCIP_Real rhs;

         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
         assert(!SCIPisInfinity(scip, row->rhs));
         rhs = row->rhs - row->constant;
         if( row->integral )
         {
            /* the right hand side was implicitly rounded down in row aggregation */
            rhs = SCIPfloor(scip, rhs);
         }

         SCIPquadprecProdQD(cutar, cutar, rhs);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -cutar);
      }
      else
      {
         SCIP_Real lhs;

         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
         assert(!SCIPisInfinity(scip, -row->lhs));
         lhs = row->lhs - row->constant;
         if( row->integral )
         {
            /* the left hand side was implicitly rounded up in row aggregation */
            lhs = SCIPceil(scip, lhs);
         }

         SCIPquadprecProdQD(cutar, cutar, lhs);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, cutar);
      }
   }

   /* relax rhs to zero, if it's very close to */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   return SCIP_OKAY;
}


/** calculates a strong CG cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in a strongcg cut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcalcStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute a strong CG cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store the efficacy of the cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid cut was returned */
   )
{
   int i;
   int nvars;
   int* varsign;
   int* boundtype;
   SCIP_Real* tmpcoefs;
   SCIP_Real QUAD(downrhs);
   SCIP_Real QUAD(f0);
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(rhs);
   SCIP_Real k;
   SCIP_Bool freevariable;
   SCIP_Bool localbdsused;

   assert(scip != NULL);
   assert(aggrrow != NULL);
   assert(SCIPisPositive(scip, scale));
   assert(cutcoefs != NULL);
   assert(cutrhs != NULL);
   assert(cutinds != NULL);
   assert(success != NULL);
   assert(cutislocal != NULL);

   SCIPdebugMessage("calculating strong CG cut (scale: %g)\n", scale);

   *success = FALSE;

   /* allocate temporary memory */
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtype, nvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &tmpcoefs, QUAD_ARRAY_SIZE(nvars)) );

   /* initialize cut with aggregation */
   *cutnnz = aggrrow->nnz;
   *cutislocal = aggrrow->local;
   SCIPquadprecProdQD(rhs, aggrrow->rhs, scale);

   if( *cutnnz > 0 )
   {
      BMScopyMemoryArray(cutinds, aggrrow->inds, *cutnnz);

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(coef);
         int j = cutinds[i];

         QUAD_ARRAY_LOAD(coef, aggrrow->vals, j);
         SCIPquadprecProdQD(coef, coef, scale);

         QUAD_HI(coef) = NONZERO(QUAD_HI(coef));
         assert(QUAD_HI(coef) != 0.0);

         QUAD_ARRAY_STORE(tmpcoefs, j, coef);
      }

      /* Transform equation  a*x == b, lb <= x <= ub  into standard form
       *   a'*x' == b, 0 <= x' <= ub'.
       *
       * Transform variables (lb or ub):
       *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
       *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
       * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
       *
       * Transform variables (vlb or vub):
       *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
       *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
       * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
       *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
       *   a_{zu_j} := a_{zu_j} + a_j * bu_j
       */
      SCIP_CALL( cutsTransformStrongCG(scip, sol, boundswitch, usevbds, allowlocal,
                                       tmpcoefs, QUAD(&rhs), cutinds, cutnnz, varsign, boundtype, &freevariable, &localbdsused) );

      assert(allowlocal || !localbdsused);
      *cutislocal = *cutislocal || localbdsused;

      if( freevariable )
         goto TERMINATE;

      SCIPdebug(printCutQuad(scip, NULL, cutcoefs, QUAD(rhs), cutinds, *cutnnz, FALSE, FALSE));
   }

   /* Calculate
    *  - fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j)
    *  - integer k >= 1 with 1/(k + 1) <= f_0 < 1/k
    *    (=> k = up(1/f_0) + 1)
    *  - integer 1 <= p_j <= k with f_0 + ((p_j - 1) * (1 - f_0)/k) < f_j <= f_0 + (p_j * (1 - f_0)/k)
    *    (=> p_j = up( (f_j - f_0)/((1 - f_0)/k) ))
    * and derive strong CG cut
    *   a~*x' <= (k+1) * down(b)
    * integers :  a~_j = down(a'_j)                , if f_j <= f_0
    *             a~_j = down(a'_j) + p_j/(k + 1)  , if f_j >  f_0
    * continuous: a~_j = 0                         , if a'_j >= 0
    *             no strong CG cut found          , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */
   SCIPquadprecEpsFloorQ(downrhs, rhs, SCIPepsilon(scip)); /*lint !e666*/

   SCIPquadprecSumQQ(f0, rhs, -downrhs);
   if( QUAD_TO_DBL(f0) < minfrac || QUAD_TO_DBL(f0) > maxfrac )
      goto TERMINATE;

   /* renormalize the f0 value */
   SCIPquadprecSumDD(f0, QUAD_HI(f0), QUAD_LO(f0));

   SCIPquadprecDivDQ(tmp, 1.0, f0);
   k = SCIPround(scip, ceil(QUAD_TO_DBL(tmp)) - 1.0);

   QUAD_ASSIGN_Q(rhs, downrhs);

   if( *cutnnz > 0 )
   {
      SCIP_CALL( cutsRoundStrongCG(scip, tmpcoefs, QUAD(&rhs), cutinds, cutnnz, varsign, boundtype, QUAD(f0), k) );
      SCIPdebug(printCutQuad(scip, sol, cutcoefs, QUAD(rhs), cutinds, *cutnnz, FALSE, FALSE));
   }

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a^_r = a~_r = (k + 1) * down(a'_r)        , if f_r <= f0
    *               a^_r = a~_r = (k + 1) * down(a'_r) + p_r  , if f_r >  f0
    *   continuous: a^_r = a~_r = 0                           , if a'_r >= 0
    *               a^_r = a~_r = a'_r/(1 - f0)               , if a'_r <  0
    *
    * Substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
    */
   SCIP_CALL( cutsSubstituteStrongCG(scip, aggrrow->rowweights, aggrrow->slacksign, aggrrow->rowsinds,
                          aggrrow->nrows, scale, tmpcoefs, QUAD(&rhs), cutinds, cutnnz, QUAD(f0), k) );
   SCIPdebug(printCutQuad(scip, sol, cutcoefs, QUAD(rhs), cutinds, *cutnnz, FALSE, FALSE));

   /* remove all nearly-zero coefficients from strong CG row and relax the right hand side correspondingly in order to
    * prevent numerical rounding errors
    */
   if( postprocess )
   {
      SCIP_CALL( postprocessCutQuad(scip, *cutislocal, cutinds, tmpcoefs, cutnnz, QUAD(&rhs), success) );
   }
   else
   {
      *success = ! removeZerosQuad(scip, SCIPsumepsilon(scip), *cutislocal, tmpcoefs, QUAD(&rhs), cutinds, cutnnz);
   }
   SCIPdebug(printCutQuad(scip, sol, cutcoefs, QUAD(rhs), cutinds, *cutnnz, FALSE, FALSE));

   if( *success )
   {
      *cutrhs = QUAD_TO_DBL(rhs);

      /* store cut in given array in sparse representation and clean buffer array */
      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(coef);
         int j = cutinds[i];

         QUAD_ARRAY_LOAD(coef, tmpcoefs, j);
         assert(QUAD_HI(coef) != 0.0);

         cutcoefs[i] = QUAD_TO_DBL(coef);
         QUAD_ASSIGN(coef, 0.0);
         QUAD_ARRAY_STORE(tmpcoefs, j, coef);
      }

      if( cutefficacy != NULL )
         *cutefficacy = calcEfficacy(scip, sol, cutcoefs, *cutrhs, cutinds, *cutnnz);

      if( cutrank != NULL )
         *cutrank = aggrrow->rank + 1;
   }

  TERMINATE:

   /* if we aborted early the tmpcoefs array needs to be cleaned */
   if( !(*success) )
   {
      QUAD_ASSIGN(tmp, 0.0);

      for( i = 0; i < *cutnnz; ++i )
      {
         QUAD_ARRAY_STORE(tmpcoefs, cutinds[i], tmp);
      }
   }

   /* free temporary memory */
   SCIPfreeCleanBufferArray(scip, &tmpcoefs);
   SCIPfreeBufferArray(scip, &boundtype);
   SCIPfreeBufferArray(scip, &varsign);

   return SCIP_OKAY;
}
