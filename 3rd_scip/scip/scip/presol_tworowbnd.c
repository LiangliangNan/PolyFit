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

/**@file   presol_tworowbnd.c
 * @brief  do bound tightening by using two rows
 * @author Dieter Weninger
 *
 * Perform bound tightening on two inequalities with some common variables.
 *
 * Let two constraints be given:
 * \f{eqnarray*}{
 *   A_{iR} x_R + A_{iS} x_S              \geq b_i\\
 *   A_{kR} x_R              + A_{kT} x_T \geq b_k
 * \f}
 * with \f$N\f$ the set of variable indexes, \f$R \subseteq N\f$, \f$S \subseteq N\f$, \f$T \subseteq N\f$,
 * \f$R \cap S = \emptyset\f$, \f$R \cap T = \emptyset\f$, \f$S \cap T = \emptyset\f$ and \f$i \not= k\f$.
 *
 * Solve the following two LPs
 * \f{eqnarray*}{
 *   L = \min \{ A_{kR} x_R : A_{iR} x_R + A_{iS} x_S \geq b_i \}\\
 *   U = \max \{ A_{kR} x_R : A_{iR} x_R + A_{iS} x_S \geq b_i \}
 * \f}
 * and use \f$L\f$ and \f$U\f$ for getting bounds on \f$x_T\f$.
 *
 * If \f$L + \mbox{infimum}(A_{kT}x_T) \geq b_k\f$, then the second constraint above is redundant.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/pub_matrix.h"
#include "presol_tworowbnd.h"

#define PRESOL_NAME            "tworowbnd"
#define PRESOL_DESC            "do bound tigthening by using two rows"
#define PRESOL_PRIORITY          -500000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               0     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define SUPPORT_THRESHOLD            0.5     /**< threshold for two constraints overlap */
#define FASTMODE_THRESHOLD          1000     /**< max number of baserows for switching to fast mode */

/* uncomment the following define to debug solving the small LPs */
/* #define DEBUG_WRITE_CHECK_LPS */

/** type of bound change */
enum Bndchgtype
{
   LOWERBOUND = 1,
   UPPERBOUND = 2,
   BOTHBOUNDS = 3
};
typedef enum Bndchgtype BNDCHGTYPE;

/*
 * Local methods
 */

#ifdef DEBUG_WRITE_CHECK_LPS
/** write min and max LP to file */
static
void writeLPs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   otherrow,           /**< other row index */
   int                   numoverlap,         /**< overlap-size */
   int*                  overlapidx,         /**< overlap column indexes */
   int*                  othernonoverlapidx, /**< other row non overlap indexes */
   SCIP_Real*            coefbaseoverlap,    /**< base row overlap coefficients */
   SCIP_Real*            coefotheroverlap,   /**< other row overlap coefficients */
   SCIP_Real*            coefothernonoverlap,/**< other row non overlap coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds            /**< upper bounds */
   )
{
   FILE* filemax;
   FILE* filemin;
   SCIP_Real lhs;
   int i;
   int nothernonolap;

   lhs = SCIPmatrixGetRowLhs(matrix, otherrow);

   filemax = fopen("max.lp", "wt");
   filemin = fopen("min.lp", "wt");
   if( filemax != NULL && filemin != NULL )
   {
      fprintf(filemax, "max\n\t");
      fprintf(filemin, "min\n\t");

      for( i = 0; i < numoverlap; i++ )
      {
         if( coefbaseoverlap[i] > 0.0 )
         {
            fprintf(filemax, "+%.24f %s ", coefbaseoverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
            fprintf(filemin, "+%.24f %s ", coefbaseoverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
         }
         else
         {
            fprintf(filemax, "%.24f %s ", coefbaseoverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
            fprintf(filemin, "%.24f %s ", coefbaseoverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
         }
      }

      fprintf(filemax, "\ns.t.\n\t");
      fprintf(filemin, "\ns.t.\n\t");

      for( i = 0; i < numoverlap; i++ )
      {
         if( coefotheroverlap[i] > 0.0 )
         {
            fprintf(filemax, "+%.24f %s ", coefotheroverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
            fprintf(filemin, "+%.24f %s ", coefotheroverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
         }
         else
         {
            fprintf(filemax, "%.24f %s ", coefotheroverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
            fprintf(filemin, "%.24f %s ", coefotheroverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
         }
      }

      nothernonolap = SCIPmatrixGetRowNNonzs(matrix, otherrow) - numoverlap;

      for( i = 0; i < nothernonolap; i++ )
      {
         if( coefothernonoverlap[i] > 0.0 )
         {
            fprintf(filemax, "+%.24f %s ", coefothernonoverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])));
            fprintf(filemin, "+%.24f %s ", coefothernonoverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])));
         }
         else
         {
            fprintf(filemax, "%.24f %s ", coefothernonoverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])));
            fprintf(filemin, "%.24f %s ", coefothernonoverlap[i], SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])));
         }
      }
      fprintf(filemax, " >= %.24f\n", lhs);
      fprintf(filemin, " >= %.24f\n", lhs);

      fprintf(filemax, "bounds\n");
      fprintf(filemin, "bounds\n");

      for( i = 0; i < numoverlap; i++ )
      {
         if( !SCIPisInfinity(scip, -lowerbds[overlapidx[i]]) && !SCIPisInfinity(scip, upperbds[overlapidx[i]]) )
         {
            fprintf(filemax, "\t%.24f <= %s <= %.24f\n", lowerbds[overlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])), upperbds[overlapidx[i]]);
            fprintf(filemin, "\t%.24f <= %s <= %.24f\n", lowerbds[overlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])), upperbds[overlapidx[i]]);
         }
         else if( !SCIPisInfinity(scip, -lowerbds[overlapidx[i]]) )
         {
            fprintf(filemax, "\t%.24f <= %s\n", lowerbds[overlapidx[i]], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
            fprintf(filemin, "\t%.24f <= %s\n", lowerbds[overlapidx[i]], SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])));
         }
         else if( !SCIPisInfinity(scip, upperbds[overlapidx[i]]) )
         {
            fprintf(filemax, "\t%s <= %.24f\n", SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])), upperbds[overlapidx[i]]);
            fprintf(filemin, "\t%s <= %.24f\n", SCIPvarGetName(SCIPmatrixGetVar(matrix, overlapidx[i])), upperbds[overlapidx[i]]);
         }
      }

      for( i = 0; i < nothernonolap; i++ )
      {
         if( !SCIPisInfinity(scip, -lowerbds[othernonoverlapidx[i]]) && !SCIPisInfinity(scip, upperbds[othernonoverlapidx[i]]) )
         {
            fprintf(filemax, "\t%.24f <= %s <= %.24f\n", lowerbds[othernonoverlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])), upperbds[othernonoverlapidx[i]]);
            fprintf(filemin, "\t%.24f <= %s <= %.24f\n", lowerbds[othernonoverlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])), upperbds[othernonoverlapidx[i]]);
         }
         else if( !SCIPisInfinity(scip, -lowerbds[othernonoverlapidx[i]]) )
         {
            fprintf(filemax, "\t%.24f <= %s\n", lowerbds[othernonoverlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])));
            fprintf(filemin, "\t%.24f <= %s\n", lowerbds[othernonoverlapidx[i]],
               SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])));
         }
         else if( !SCIPisInfinity(scip, upperbds[othernonoverlapidx[i]]) )
         {
            fprintf(filemax, "\t%s <= %.24f\n", SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])), upperbds[othernonoverlapidx[i]]);
            fprintf(filemin, "\t%s <= %.24f\n", SCIPvarGetName(SCIPmatrixGetVar(matrix, othernonoverlapidx[i])), upperbds[othernonoverlapidx[i]]);
         }
      }


      fprintf(filemax, "end\n");
      fprintf(filemin, "end\n");

      fclose(filemax);
      fclose(filemin);
   }
   else
      SCIPABORT();
}
#endif


/** solve two LPs with one row (single constraint) each
 *
 * a1x + a3y      >= b1  (other row)
 * a2x      + a4z >= b2  (base row)
 *
 * minact = min{a2x : a1x + a3y >= b1}
 * maxact = max{a2x : a1x + a3y >= b1}
 */
static
void getActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int                   numoverlap,         /**< overlap-size */
   int*                  overlapidx,         /**< overlap column indexes */
   int*                  othernonoverlapidx, /**< other row non overlap indexes */
   SCIP_Real*            coefbaseoverlap,    /**< base row overlap coefficients */
   SCIP_Real*            coefotheroverlap,   /**< other row overlap coefficients */
   SCIP_Real*            coefothernonoverlap,/**< other row non overlap coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds,           /**< upper bounds */
   SCIP_Real*            tmplowerbds,        /**< tmp lower bounds */
   SCIP_Real*            tmpupperbds,        /**< tmp upper bounds */
   SCIP_Real*            minratios,          /**< min LP ratios */
   SCIP_Real*            maxratios,          /**< max LP ratios */
   int*                  minsortedidx,       /**< min LP sorted indexes */
   int*                  maxsortedidx,       /**< max LP sorted indexes */
   SCIP_Real*            minact,             /**< calculated overlap minimal activity w.r.t. to the other row */
   SCIP_Real*            maxact              /**< calculated overlap maximal activity w.r.t. to the other row */
   )
{/*lint --e{715}*/
   SCIP_Real val;
   int nothernonoverlap;
   SCIP_Real lhs;
   SCIP_Real minlhs;
   SCIP_Real maxlhs;
   SCIP_Bool consred;
   int nminratios;
   int nmaxratios;
   int i;

   *minact = 0;
   *maxact = 0;

#ifdef DEBUG_WRITE_CHECK_LPS
   SCIPmatrixPrintRow(scip,matrix,baserow);
   SCIPmatrixPrintRow(scip,matrix,otherrow);
   writeLPs(scip, matrix, otherrow, numoverlap, overlapidx, othernonoverlapidx,
      coefbaseoverlap, coefotheroverlap, coefothernonoverlap, lowerbds, upperbds);
#endif

   lhs = SCIPmatrixGetRowLhs(matrix, otherrow);
   assert(!SCIPisInfinity(scip, -lhs));

   nothernonoverlap = SCIPmatrixGetRowNNonzs(matrix, otherrow) - numoverlap;
   val = 0;
   consred = FALSE;

   /* compute maximal contribution of non-overlap part of the
      single constraint to the activity.
      maybe the single constraint is redundant */
   for( i = 0; i < nothernonoverlap; i++ )
   {
      if( coefothernonoverlap[i] < 0.0 )
      {
         if( SCIPisInfinity(scip, -lowerbds[othernonoverlapidx[i]]) )
         {
            consred = TRUE;
            break;
         }
         else
         {
            val += coefothernonoverlap[i] * lowerbds[othernonoverlapidx[i]];
         }
      }
      else if( coefothernonoverlap[i] > 0.0 )
      {
         if( SCIPisInfinity(scip, upperbds[othernonoverlapidx[i]]) )
         {
            consred = TRUE;
            break;
         }
         else
         {
            val += coefothernonoverlap[i] * upperbds[othernonoverlapidx[i]];
         }
      }
   }

   if( !consred )
   {
      SCIP_Real minlowerbnd;
      minlowerbnd = SCIPinfinity(scip);

      /* we want that every coefficient in the single constraint
         has a negative sign and hence we need to multiply
         some columns by -1 */
      for( i = 0; i < numoverlap; i++ )
      {
         tmplowerbds[i] = lowerbds[overlapidx[i]];
         tmpupperbds[i] = upperbds[overlapidx[i]];

         if( coefotheroverlap[i] > 0.0 )
         {
            /* multiply column by -1 and swap bounds */
            double tmp;
            tmp = tmplowerbds[i];
            tmplowerbds[i] = -tmpupperbds[i];
            tmpupperbds[i] = -tmp;

            coefotheroverlap[i] = -coefotheroverlap[i];
            coefbaseoverlap[i] = -coefbaseoverlap[i];
         }

         if( tmplowerbds[i] < minlowerbnd )
         {
            if( SCIPisInfinity(scip, -tmplowerbds[i]) )
            {
               /* lower bounds have to be finite for later boundshift */
               *minact = -SCIPinfinity(scip);
               *maxact = SCIPinfinity(scip);
               return;
            }
            minlowerbnd = tmplowerbds[i];
         }
      }

      if( minlowerbnd < 0.0 )
      {
         SCIP_Real bndshift = -minlowerbnd;
         if( bndshift > (SCIP_Real)SCIP_LONGINT_MAX )
         {
            /* shift value is too large */
            *minact = -SCIPinfinity(scip);
            *maxact = SCIPinfinity(scip);
            return;
         }
      }

      /* init left hand side values and consider non-overlapping contribution */
      minlhs = lhs - val;
      nminratios = 0;
      maxlhs = lhs - val;
      nmaxratios = 0;

      if( minlowerbnd < 0.0 )
      {
         double bndshift = -minlowerbnd;
         if( bndshift > (double)SCIP_LONGINT_MAX )
         {
            /* shift value is too large */
            *minact = -SCIPinfinity(scip);
            *maxact = SCIPinfinity(scip);
            return;
         }

         /* shift polyhedra into the positive orthant */
         for( i = 0; i < numoverlap; i++ )
         {
            minlhs += coefotheroverlap[i] * bndshift;
            *minact -= coefbaseoverlap[i] * bndshift;

            maxlhs += coefotheroverlap[i] * bndshift;
            *maxact -= coefbaseoverlap[i] * bndshift;

            tmplowerbds[i] += bndshift;
            tmpupperbds[i] += bndshift;

            assert(tmplowerbds[i] >= 0.0);
         }
      }

      /*
       * solve minimization LP
       *
       * we distinguish two column cases:
       *
       * a)           b)
       * min  -       min  +
       * s.t. -       s.t. -
       */

      for( i = 0; i < numoverlap; i++ )
      {
         if( coefbaseoverlap[i] > 0.0 )
         {
            /* b): fix variable to its lower bound */
            minlhs -= coefotheroverlap[i] * tmplowerbds[i];
            *minact += coefbaseoverlap[i] * tmplowerbds[i];
         }
         else
         {
            /* a): save coefficient ratios for later sorting */
            minratios[nminratios] = coefbaseoverlap[i] / coefotheroverlap[i];
            minsortedidx[nminratios] = i;
            nminratios++;

            /* consider lower bounds for left hand side and obj value */
            minlhs -= coefotheroverlap[i] * tmplowerbds[i];
            *minact += coefbaseoverlap[i] * tmplowerbds[i];
         }
      }

      /* sort the ratios for case a) */
      if( nminratios > 1 )
         SCIPsortRealInt(minratios, minsortedidx, nminratios);

      /* pack every variable on the highest possible value as long as we are feasible */
      for( i = nminratios-1; 0 <= i; i-- )
      {
         SCIP_Real tmpval;

         /* consider contribution from lower bounds */
         if( tmplowerbds[minsortedidx[i]] > 0 )
         {
            *minact -= coefbaseoverlap[minsortedidx[i]] * tmplowerbds[minsortedidx[i]];
            minlhs += coefotheroverlap[minsortedidx[i]] * tmplowerbds[minsortedidx[i]];
         }

         /* calculate highest possible variable value */
         tmpval = minlhs / coefotheroverlap[minsortedidx[i]];
         if( tmpval < tmplowerbds[minsortedidx[i]] )
         {
            /* infeasible */
            *minact = -SCIPinfinity(scip);
            break;
         }

         /* if the upper bound is large enough we are ready
            otherwise we set the variable to its upper bound and iterate */
         if( tmpval <= tmpupperbds[minsortedidx[i]] )
         {
            *minact += coefbaseoverlap[minsortedidx[i]] * tmpval;
            break;
         }
         else
         {
            *minact += coefbaseoverlap[minsortedidx[i]] * tmpupperbds[minsortedidx[i]];
            minlhs -= coefotheroverlap[minsortedidx[i]] * tmpupperbds[minsortedidx[i]];
         }
      }


      /*
       * solve maximization LP
       *
       * we distinguish two column cases:
       *
       * c)           d)
       * max  +       max  -
       * s.t. -       s.t. -
       */
      for( i = 0; i < numoverlap; i++ )
      {
         if( coefbaseoverlap[i] < 0.0 )
         {
            /* d): fix variable to its lower bound */
            maxlhs -= coefotheroverlap[i] * tmplowerbds[i];
            *maxact += coefbaseoverlap[i] * tmplowerbds[i];
         }
         else
         {
            /* c): save coefficient ratios for later sorting */
            maxratios[nmaxratios] = coefbaseoverlap[i] / coefotheroverlap[i];
            maxsortedidx[nmaxratios] = i;
            nmaxratios++;

            /* consider lower bounds for left hand side and obj value */
            maxlhs -= coefotheroverlap[i] * tmplowerbds[i];
            *maxact += coefbaseoverlap[i] * tmplowerbds[i];
         }
      }

      /* sort the ratios for case a) */
      if( nmaxratios > 1 )
         SCIPsortRealInt(maxratios, maxsortedidx, nmaxratios);

      /* pack every variable on the highest possible value as long as we are feasible */
      for( i = 0; i < nmaxratios; i++ )
      {
         SCIP_Real tmpval;

         /* consider contribution from lower bounds */
         if( tmplowerbds[maxsortedidx[i]] > 0 )
         {
            *maxact -= coefbaseoverlap[maxsortedidx[i]] * tmplowerbds[maxsortedidx[i]];
            maxlhs += coefotheroverlap[maxsortedidx[i]] * tmplowerbds[maxsortedidx[i]];
         }

         /* calculate highest possible variable value */
         tmpval = maxlhs / coefotheroverlap[maxsortedidx[i]];
         if( tmpval < tmplowerbds[maxsortedidx[i]] )
         {
            /* infeasible */
            *maxact = SCIPinfinity(scip);
            break;
         }

         /* if the upper bound is large enough we are ready
            otherwise we set the variable to its upper bound and iterate */
         if( tmpval <= tmpupperbds[maxsortedidx[i]] )
         {
            *maxact += coefbaseoverlap[maxsortedidx[i]] * tmpval;
            break;
         }
         else
         {
            *maxact += coefbaseoverlap[maxsortedidx[i]] * tmpupperbds[maxsortedidx[i]];
            maxlhs -= coefotheroverlap[maxsortedidx[i]] * tmpupperbds[maxsortedidx[i]];
         }
      }
   }
   else
   {
      /* single constraint is redundant.
         we calculate the value of the objective function */

      /* minimization LP */
      for( i = 0; i < numoverlap; i++ )
      {
         if( coefbaseoverlap[i] > 0.0 )
         {
            if( !SCIPisInfinity(scip, -lowerbds[overlapidx[i]]) )
            {
               *minact += coefbaseoverlap[i] * lowerbds[overlapidx[i]];
            }
            else
            {
               *minact = -SCIPinfinity(scip);
               break;
            }
         }
         else if( coefbaseoverlap[i] < 0.0 )
         {
            if( !SCIPisInfinity(scip, upperbds[overlapidx[i]]) )
            {
               *minact += coefbaseoverlap[i] * upperbds[overlapidx[i]];
            }
            else
            {
               *minact = -SCIPinfinity(scip);
               break;
            }
         }
      }

      /* maximization LP */
      for( i = 0; i < numoverlap; i++ )
      {
         if( coefbaseoverlap[i] > 0.0 )
         {
            if( !SCIPisInfinity(scip, upperbds[overlapidx[i]]) )
            {
               *maxact += coefbaseoverlap[i] * upperbds[overlapidx[i]];
            }
            else
            {
               *maxact = SCIPinfinity(scip);
               break;
            }
         }
         else if( coefbaseoverlap[i] < 0.0 )
         {
            if( !SCIPisInfinity(scip, -lowerbds[overlapidx[i]]) )
            {
               *maxact += coefbaseoverlap[i] * lowerbds[overlapidx[i]];
            }
            else
            {
               *maxact = SCIPinfinity(scip);
               break;
            }
         }
      }
   }

#ifdef DEBUG_WRITE_CHECK_LPS
   {
      SCIP_Real minsolve = 0.0;
      SCIP* subscip;
      SCIP_SOL* sol;
      SCIP_STATUS status;
      SCIP_VAR** vars;
      int nvars;
      int objnonzeros = 0;
      SCIP_CALL_ABORT( SCIPcreate(&subscip) );
      SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(subscip) );
      SCIP_CALL_ABORT( SCIPreadProb(subscip, "min.lp", NULL) );
      SCIP_CALL_ABORT( SCIPsetIntParam(subscip,"presolving/maxrounds",0) );
      vars = SCIPgetVars(subscip);
      nvars = SCIPgetNVars(subscip);
      for(i=0; i< nvars; i++)
      {
         if(SCIPvarGetObj(vars[i]) != 0.0)
            objnonzeros++;
      }
      assert(numoverlap == objnonzeros);

      SCIP_CALL_ABORT( SCIPsolve(subscip) );
      status = SCIPgetStatus(subscip);
      if(SCIP_STATUS_OPTIMAL == status)
      {
         sol = SCIPgetBestSol(subscip);
         minsolve = SCIPgetSolOrigObj(subscip, sol);
         assert(SCIPisEQ(scip,minsolve,*minact));
      }
      else
      {
         assert(SCIPisEQ(scip,-SCIPinfinity(scip),*minact));
      }
      SCIP_CALL_ABORT( SCIPfree(&subscip) );
   }
   {
      SCIP_Real maxsolve = 0.0;
      SCIP* subscip;
      SCIP_SOL* sol;
      SCIP_STATUS status;
      SCIP_CALL_ABORT( SCIPcreate(&subscip) );
      SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(subscip) );
      SCIP_CALL_ABORT( SCIPreadProb(subscip, "max.lp", NULL) );
      SCIP_CALL_ABORT( SCIPsetIntParam(subscip,"presolving/maxrounds",0) );
      SCIP_CALL_ABORT( SCIPsolve(subscip) );
      status = SCIPgetStatus(subscip);
      if(SCIP_STATUS_OPTIMAL == status)
      {
         sol = SCIPgetBestSol(subscip);
         maxsolve = SCIPgetSolOrigObj(subscip, sol);
         assert(SCIPisEQ(scip,maxsolve,*maxact));
      }
      else
      {
         assert(SCIPisEQ(scip,SCIPinfinity(scip),*maxact));
      }
      SCIP_CALL_ABORT( SCIPfree(&subscip) );
   }
#endif
}/*lint !e438*/

/** calculate min activity */
static
SCIP_Real getMinActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   len,                /**< length */
   int*                  varidxs,            /**< variables indexes */
   SCIP_Real*            coefs,              /**< coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds            /**< upper bounds */
   )
{
   SCIP_Real infimum;
   int i;

   infimum = 0;

   for( i = 0; i < len; i++ )
   {
      if( coefs[i] > 0.0 )
      {
         if( SCIPisInfinity(scip, -lowerbds[varidxs[i]]) )
         {
            infimum = -SCIPinfinity(scip);
            break;
         }
         else
         {
            infimum += coefs[i] * lowerbds[varidxs[i]];
         }
      }
      else
      {
         if( SCIPisInfinity(scip, upperbds[varidxs[i]]) )
         {
            infimum = -SCIPinfinity(scip);
            break;
         }
         else
         {
            infimum += coefs[i] * upperbds[varidxs[i]];
         }
      }
   }

   return infimum;
}

/** calculate max activity */
static
SCIP_Real getMaxActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   len,                /**< length */
   int*                  varidxs,            /**< variable indexes */
   SCIP_Real*            coefs,              /**< coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds,           /**< upper bounds */
   int*                  infcnt              /**< infinity counter */
   )
{
   SCIP_Real supremum;
   int i;

   *infcnt = 0;
   supremum = 0;

   for( i = 0; i < len; i++ )
   {
      if( coefs[i] < 0.0 )
      {
         if( SCIPisInfinity(scip, -lowerbds[varidxs[i]]) )
            (*infcnt)++;
         else
            supremum += coefs[i] * lowerbds[varidxs[i]];
      }
      else
      {
         if( SCIPisInfinity(scip, upperbds[varidxs[i]]) )
            (*infcnt)++;
         else
            supremum += coefs[i] * upperbds[varidxs[i]];
      }
   }

   if(*infcnt > 0)
      supremum = SCIPinfinity(scip);

   return supremum;
}

/** get max activity without one column */
static
SCIP_Real getMaxResActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   len,                /**< length */
   int*                  varidxs,            /**< variable indexes */
   SCIP_Real*            coefs,              /**< coefficients */
   SCIP_Real*            lowerbds,           /**< upper bounds */
   SCIP_Real*            upperbds,           /**< lower bounds */
   int                   idx                 /**< omitting index */
   )
{
   SCIP_Real supremum;
   int i;

   supremum = 0;

   for( i = 0; i < len; i++ )
   {
      if( i == idx )
         continue;

      if( coefs[i] < 0.0 )
      {
         assert(!SCIPisInfinity(scip, -lowerbds[varidxs[i]]));
         supremum += coefs[i] * lowerbds[varidxs[i]];
      }
      else
      {
         assert(!SCIPisInfinity(scip, upperbds[varidxs[i]]));
         supremum += coefs[i] * upperbds[varidxs[i]];
      }
   }

   return supremum;
}

/** apply bound tightening on two overlapping constraints */
static
void applyTightening(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int                   numoverlap,         /**< overlap-size */
   int*                  overlapidx,         /**< overlap column indexes */
   int*                  othernonoverlapidx, /**< other row non overlap indexes */
   int*                  basenonoverlapidx,  /**< base row non overlap indexes */
   SCIP_Real*            coefbaseoverlap,    /**< base row overlap coefficients */
   SCIP_Real*            coefotheroverlap,   /**< other row overlap coefficients */
   SCIP_Real*            coefbasenonoverlap, /**< base row non overlap coefficients */
   SCIP_Real*            coefothernonoverlap,/**< other row non overlap coefficients */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds,           /**< upper bounds */
   SCIP_Real*            tmplowerbds,        /**< tmp lower bounds */
   SCIP_Real*            tmpupperbds,        /**< tmp upper bounds */
   SCIP_Real*            minratios,          /**< min LP ratios */
   SCIP_Real*            maxratios,          /**< max LP ratios */
   int*                  minsortedidx,       /**< min LP sorted indexes */
   int*                  maxsortedidx,       /**< max LP sorted indexes */
   int*                  ntightenbnds,       /**< number of tightened bounds */
   BNDCHGTYPE*           tighten,            /**< tightened bounds */
   int*                  ndeletecons,        /**< number of redundant constraints */
   SCIP_Bool*            deletecons          /**< redundant constraints */
   )
{
   SCIP_Real maxactoverlap;
   SCIP_Real minactoverlap;
   SCIP_Real minactnonoverlap;
   SCIP_Real maxactnonoverlap;
   int len;
   SCIP_Real lhs;
   int i;

   getActivities(scip, matrix, baserow, otherrow, numoverlap, overlapidx, othernonoverlapidx,
      coefbaseoverlap, coefotheroverlap, coefothernonoverlap,
      lowerbds, upperbds, tmplowerbds, tmpupperbds, minratios, maxratios,
      minsortedidx, maxsortedidx, &minactoverlap, &maxactoverlap);

   len = SCIPmatrixGetRowNNonzs(matrix, baserow) - numoverlap;
   lhs = SCIPmatrixGetRowLhs(matrix, baserow);

   if( !SCIPisInfinity(scip, -minactoverlap) )
   {
      /* detect redundant constraints */
      minactnonoverlap = getMinActivity(scip, len, basenonoverlapidx, coefbasenonoverlap, lowerbds, upperbds);
      if( !SCIPisInfinity(scip, -minactnonoverlap) )
      {
         if( SCIPisGE(scip, minactoverlap + minactnonoverlap, lhs) )
         {
            if( !deletecons[baserow] )
            {
               (*ndeletecons)++;
               deletecons[baserow] = TRUE;
            }
         }
      }
   }

   if( !SCIPisInfinity(scip, maxactoverlap) )
   {
      int infcnt;
      SCIP_Real bnd;
      SCIP_Real tmpsup;

      /* bound tightening */
      maxactnonoverlap = getMaxActivity(scip, len, basenonoverlapidx, coefbasenonoverlap, lowerbds, upperbds, &infcnt);
      if( !SCIPisInfinity(scip, maxactnonoverlap) )
      {
         for( i = 0; i < len; i++ )
         {
            if( coefbasenonoverlap[i] < 0.0 )
            {
               /* get ub */
               tmpsup = maxactnonoverlap - (coefbasenonoverlap[i] * lowerbds[basenonoverlapidx[i]]);
               bnd = (lhs - (tmpsup + maxactoverlap)) / coefbasenonoverlap[i];
               if( bnd < upperbds[basenonoverlapidx[i]] )
               {
                  upperbds[basenonoverlapidx[i]] = bnd;
                  if( tighten[basenonoverlapidx[i]] != UPPERBOUND && tighten[basenonoverlapidx[i]] != BOTHBOUNDS )
                  {
                     (*ntightenbnds)++;
                     if( tighten[basenonoverlapidx[i]] == LOWERBOUND )
                        tighten[basenonoverlapidx[i]] = BOTHBOUNDS;
                     else
                        tighten[basenonoverlapidx[i]] = UPPERBOUND;
                  }
               }
            }
            else
            {
               /* get lb */
               tmpsup = maxactnonoverlap - (coefbasenonoverlap[i] * upperbds[basenonoverlapidx[i]]);
               bnd = (lhs - (tmpsup + maxactoverlap)) / coefbasenonoverlap[i];
               if( bnd > lowerbds[basenonoverlapidx[i]] )
               {
                  lowerbds[basenonoverlapidx[i]] = bnd;
                  if( tighten[basenonoverlapidx[i]] != LOWERBOUND && tighten[basenonoverlapidx[i]] != BOTHBOUNDS )
                  {
                     (*ntightenbnds)++;
                     if( tighten[basenonoverlapidx[i]] == UPPERBOUND )
                        tighten[basenonoverlapidx[i]] = BOTHBOUNDS;
                     else
                        tighten[basenonoverlapidx[i]] = LOWERBOUND;
                  }
               }
            }
         }
      }
      /* maximal activity in non-overlapping variables is +infinity */
      else
      {
         /* we can only do bound tightening, if we have exactly one infinite contribution*/
         if( infcnt == 1 )
         {
            for( i = 0; i < len; i++ )
            {
               if( coefbasenonoverlap[i] < 0.0 )
               {
                  if( SCIPisInfinity(scip, -lowerbds[basenonoverlapidx[i]]) )
                  {
                     /* get ub */
                     tmpsup = getMaxResActivity(scip, len, basenonoverlapidx, coefbasenonoverlap, lowerbds, upperbds, i);
                     assert(!SCIPisInfinity(scip, tmpsup));
                     bnd = (lhs - (tmpsup + maxactoverlap)) / coefbasenonoverlap[i];
                     if( bnd < upperbds[basenonoverlapidx[i]] )
                     {
                        upperbds[basenonoverlapidx[i]] = bnd;
                        if( tighten[basenonoverlapidx[i]] != UPPERBOUND && tighten[basenonoverlapidx[i]] != BOTHBOUNDS )
                        {
                           (*ntightenbnds)++;
                           if( tighten[basenonoverlapidx[i]] == LOWERBOUND )
                              tighten[basenonoverlapidx[i]] = BOTHBOUNDS;
                           else
                              tighten[basenonoverlapidx[i]] = UPPERBOUND;
                        }
                     }
                  }
               }
               else
               {
                  if( infcnt == 1 && SCIPisInfinity(scip, upperbds[basenonoverlapidx[i]]) )
                  {
                     /* get lb */
                     tmpsup = getMaxResActivity(scip, len, basenonoverlapidx, coefbasenonoverlap, lowerbds, upperbds, i);
                     assert(!SCIPisInfinity(scip, tmpsup));
                     bnd = (lhs - (tmpsup + maxactoverlap)) / coefbasenonoverlap[i];
                     if( bnd > lowerbds[basenonoverlapidx[i]] )
                     {
                        lowerbds[basenonoverlapidx[i]] = bnd;
                        if( tighten[basenonoverlapidx[i]] != LOWERBOUND && tighten[basenonoverlapidx[i]] != BOTHBOUNDS )
                        {
                           (*ntightenbnds)++;
                           if( tighten[basenonoverlapidx[i]] == UPPERBOUND )
                              tighten[basenonoverlapidx[i]] = BOTHBOUNDS;
                           else
                              tighten[basenonoverlapidx[i]] = LOWERBOUND;
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

/** extract coefficients from matrix */
static
void getCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int                   numoverlap,         /**< overlap-size */
   int*                  olapidxbaseorder,   /**< overlap column indexes in baserow order */
   int*                  olapidxotherorder,  /**< overlap column indexes in otherrow order */
   int*                  othernonoverlapidx, /**< other row non overlap indexes */
   int*                  basenonoverlapidx,  /**< base row non overlap indexes */
   SCIP_Real*            coefbaseoverlap,    /**< base row overlap coefficients */
   SCIP_Real*            coefotheroverlap,   /**< other row overlap coefficients */
   SCIP_Real*            coefbasenonoverlap, /**< base row non overlap coefficients */
   SCIP_Real*            coefothernonoverlap /**< other row non overlap coefficients */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   int baserowcnt;
   int otherrowcnt;
   int olapcnt;
   int nonolapcnt;

   /* get number of columns in the rows */
   baserowcnt = SCIPmatrixGetRowNNonzs(matrix, baserow);
   otherrowcnt = SCIPmatrixGetRowNNonzs(matrix, otherrow);
   assert(baserowcnt != 0 && otherrowcnt != 0);

#if 1 /* @todo why do we need this? */
   /* set end marker */
   if( numoverlap < SCIPmatrixGetNColumns(matrix) )
   {
      olapidxbaseorder[numoverlap] = -1;
      olapidxotherorder[numoverlap] = -1;
   }
#endif

   olapcnt = 0;
   nonolapcnt = 0;

   /* partition columns of base row into overlapping columns and non-overlapping columns (w.r.t. other row) and store
    * the corresponding coefficients
    */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, baserow);
   rowend = rowpnt + baserowcnt;
   valpnt = SCIPmatrixGetRowValPtr(matrix, baserow);

   for( ; rowpnt < rowend; rowpnt++, valpnt++ )
   {
      if( olapidxbaseorder[olapcnt] == *rowpnt )
      {
         coefbaseoverlap[olapcnt] = *valpnt;
         olapcnt++;
      }
      else
      {
         basenonoverlapidx[nonolapcnt] = *rowpnt;
         coefbasenonoverlap[nonolapcnt] = *valpnt;
         nonolapcnt++;
      }
   }

   assert(olapcnt+nonolapcnt == baserowcnt);
   assert(olapcnt == numoverlap);
   assert(nonolapcnt > 0);

   olapcnt = 0;
   nonolapcnt = 0;

   /* partition columns of other row into overlapping columns and non-overlapping columns (w.r.t. base row) and store
    * the corresponding coefficients
    */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, otherrow);
   rowend = rowpnt + otherrowcnt;
   valpnt = SCIPmatrixGetRowValPtr(matrix, otherrow);

   for( ; rowpnt < rowend; rowpnt++, valpnt++ )
   {
      if( olapidxotherorder[olapcnt] == *rowpnt )
      {
         coefotheroverlap[olapcnt] = *valpnt;
         olapcnt++;
      }
      else
      {
         othernonoverlapidx[nonolapcnt] = *rowpnt;
         coefothernonoverlap[nonolapcnt] = *valpnt;
         nonolapcnt++;
      }
   }

   assert(olapcnt+nonolapcnt == otherrowcnt);
   assert(olapcnt == numoverlap);
}

/** calculate overlap-size */
static
void getNumOverlap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int*                  countings,          /**< overlap counting helper array */
   int*                  clearinfo,          /**< reset helper array */
   int*                  numoverlap,         /**< overlap-size */
   int*                  olapidxotherorder   /**< overlap column indexes in otherrow order */
   )
{
   int* rowpnt;
   int* rowend;
   int noverlap;
   int baserowcnt;
   int otherrowcnt;
   int nclear;
   int i;

   noverlap = 0;
   nclear = 0;

   baserowcnt = SCIPmatrixGetRowNNonzs(matrix, baserow);
   otherrowcnt = SCIPmatrixGetRowNNonzs(matrix, otherrow);
   if( baserowcnt == 0 || otherrowcnt == 0 )
   {
      *numoverlap = noverlap;
      return;
   }

   /* set flags corresponding to baserow non-zeros */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, baserow);
   rowend = rowpnt + baserowcnt;
   for( ; rowpnt < rowend; rowpnt++ )
   {
      countings[*rowpnt] = 1;
      clearinfo[nclear] = *rowpnt;
      nclear++;
   }

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, otherrow);
   rowend = rowpnt + otherrowcnt;
   for( ; rowpnt < rowend; rowpnt++ )
   {
      if( countings[*rowpnt] == 1 )
      {
         /* collect overlapping indexes in otherrow order */
         olapidxotherorder[noverlap] = *rowpnt;
         noverlap++;
      }
   }

   for( i = 0; i < nclear; i++ )
      countings[clearinfo[i]] = 0;

   *numoverlap = noverlap;
}

static
void getOverlapBaseOrdered(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   baserow,            /**< base row index */
   int                   otherrow,           /**< other row index */
   int*                  countings,          /**< overlap counting helper array */
   int*                  clearinfo,          /**< reset helper array */
   int                   numoverlap,         /**< just calculated overlap-size */
   int*                  olapidxbaseorder    /**< overlap column indexes in baserow order */
   )
{
   int* rowpnt;
   int* rowend;
   int noverlap;
   int baserowcnt;
   int otherrowcnt;
   int nclear;
   int i;

   noverlap = 0;
   nclear = 0;

   baserowcnt = SCIPmatrixGetRowNNonzs(matrix, baserow);
   otherrowcnt = SCIPmatrixGetRowNNonzs(matrix, otherrow);

   /* set flags corresponding to otherrow non-zeros */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, otherrow);
   rowend = rowpnt + otherrowcnt;
   for( ; rowpnt < rowend; rowpnt++ )
   {
      countings[*rowpnt] = 1;
      clearinfo[nclear] = *rowpnt;
      nclear++;
   }

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, baserow);
   rowend = rowpnt + baserowcnt;
   for( ; rowpnt < rowend; rowpnt++ )
   {
      if( countings[*rowpnt] == 1 )
      {
         /* collect overlapping indexes in baserow order */
         olapidxbaseorder[noverlap] = *rowpnt;
         noverlap++;
      }
   }

   for( i = 0; i < nclear; i++ )
      countings[clearinfo[i]] = 0;

   assert(noverlap == numoverlap);
}


/** perform bound tightening on two rows with a specific support intersection */
static
SCIP_RETCODE calcTwoRowBnds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   nbaserows,          /**< number of base rows */
   int*                  baserows,           /**< base rows indexes */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds,           /**< upper bounds */
   int*                  ntightenbnds,       /**< number of tightened bounds */
   BNDCHGTYPE*           tighten,            /**< bound tighten information */
   int*                  ndeletecons,        /**< number of redundant constraints */
   SCIP_Bool*            deletecons          /**< redundant constraints */
   )
{
   int* rowpnt;
   int* rowend;
   int rowcnt;
   int br;
   int col;
   int* colpnt;
   int* colend;
   int colcnt;
   int numoverlap;
   int* olapidxbaseorder;
   int* olapidxotherorder;
   SCIP_Real threshold;
   int rowcnt2;
   int nrows;
   int ncols;
   int* othernonoverlapidx;
   int* basenonoverlapidx;
   SCIP_Real* coefbaseoverlap;
   SCIP_Real* coefotheroverlap;
   SCIP_Real* coefbasenonoverlap;
   SCIP_Real* coefothernonoverlap;
   int* countings;
   int* clearinfo;
   SCIP_Real* tmplowerbds;
   SCIP_Real* tmpupperbds;
   SCIP_Real* minratios;
   SCIP_Real* maxratios;
   int* minsortedidx;
   int* maxsortedidx;
   SCIP_Bool usefastmode;
   int* ignorerowidx;
   SCIP_Bool* ignorerow;
   int ignorerowcnt;
   int i;

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &olapidxbaseorder, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &olapidxotherorder, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &othernonoverlapidx, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basenonoverlapidx, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefbaseoverlap, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefotheroverlap, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefbasenonoverlap, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefothernonoverlap, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &countings, ncols) );
   BMSclearMemoryArray(countings, ncols);
   SCIP_CALL( SCIPallocBufferArray(scip, &clearinfo, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmplowerbds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpupperbds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minratios, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxratios, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minsortedidx, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxsortedidx, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ignorerowidx, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ignorerow, nrows) );
   BMSclearMemoryArray(ignorerow, nrows);

   /* use fast mode if too many base rows are present */
   if( nbaserows > FASTMODE_THRESHOLD )
      usefastmode = TRUE;
   else
      usefastmode = FALSE;

   for( br = 0; br < nbaserows; br++ )
   {
      ignorerowcnt = 0;

      rowpnt = SCIPmatrixGetRowIdxPtr(matrix, baserows[br]);
      rowcnt = SCIPmatrixGetRowNNonzs(matrix, baserows[br]);
      if( rowcnt == 0 )
         continue;

      rowend = rowpnt + rowcnt;

      for( ; (rowpnt < rowend); rowpnt++ )
      {
         col = *rowpnt;
         colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
         colcnt = SCIPmatrixGetColNNonzs(matrix, col);
         colend = colpnt + colcnt;
         for( ; (colpnt < colend); colpnt++ )
         {
            if( *colpnt == baserows[br] || ignorerow[*colpnt] )
               continue;

            /* we consider only >= constraints */
            if( !SCIPmatrixIsRowRhsInfinity(matrix, *colpnt) )
               continue;

            /* determine overlap-size */
            getNumOverlap(scip, matrix, baserows[br], *colpnt,
               countings, clearinfo, &numoverlap, olapidxotherorder);

            if( numoverlap == 0 )
               continue;

            rowcnt2 = SCIPmatrixGetRowNNonzs(matrix, *colpnt);
            threshold = (SCIP_Real)numoverlap/(SCIP_Real)MIN(rowcnt, rowcnt2);

            /* verify if overlap-size is ok */
            if( SUPPORT_THRESHOLD <= threshold && numoverlap < rowcnt )
            {
               getOverlapBaseOrdered(scip, matrix, baserows[br], *colpnt,
                  countings, clearinfo, numoverlap, olapidxbaseorder);

               getCoefficients(scip, matrix, baserows[br], *colpnt, numoverlap,
                  olapidxbaseorder, olapidxotherorder, othernonoverlapidx, basenonoverlapidx,
                  coefbaseoverlap, coefotheroverlap, coefbasenonoverlap, coefothernonoverlap);

               applyTightening(scip, matrix, baserows[br], *colpnt, numoverlap, olapidxotherorder,
                  othernonoverlapidx, basenonoverlapidx,
                  coefbaseoverlap, coefotheroverlap, coefbasenonoverlap, coefothernonoverlap,
                  lowerbds, upperbds, tmplowerbds, tmpupperbds, minratios, maxratios,
                  minsortedidx, maxsortedidx, ntightenbnds, tighten, ndeletecons, deletecons);
            }

            ignorerow[*colpnt] = TRUE;
            ignorerowidx[ignorerowcnt] = *colpnt;
            ignorerowcnt++;

            if( usefastmode )
               break;
         }

         if( usefastmode )
            break;
      }

      for( i = 0; i < ignorerowcnt; i++ )
         ignorerow[ignorerowidx[i]] = FALSE;
   }

   SCIPfreeBufferArray(scip, &ignorerow);
   SCIPfreeBufferArray(scip, &ignorerowidx);
   SCIPfreeBufferArray(scip, &maxsortedidx);
   SCIPfreeBufferArray(scip, &minsortedidx);
   SCIPfreeBufferArray(scip, &maxratios);
   SCIPfreeBufferArray(scip, &minratios);
   SCIPfreeBufferArray(scip, &tmpupperbds);
   SCIPfreeBufferArray(scip, &tmplowerbds);
   SCIPfreeBufferArray(scip, &clearinfo);
   SCIPfreeBufferArray(scip, &countings);
   SCIPfreeBufferArray(scip, &coefothernonoverlap);
   SCIPfreeBufferArray(scip, &coefbasenonoverlap);
   SCIPfreeBufferArray(scip, &coefotheroverlap);
   SCIPfreeBufferArray(scip, &coefbaseoverlap);
   SCIPfreeBufferArray(scip, &basenonoverlapidx);
   SCIPfreeBufferArray(scip, &othernonoverlapidx);
   SCIPfreeBufferArray(scip, &olapidxotherorder);
   SCIPfreeBufferArray(scip, &olapidxbaseorder);

   return SCIP_OKAY;
}

/** determine base rows */
static
SCIP_RETCODE getBaseRows(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int*                  nbaserows,          /**< number of present base rows */
   int*                  baserows            /**< indexes of base rows */
   )
{
   int nrows;
   int fill;
   int r;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(nbaserows != NULL);
   assert(baserows != NULL);

   nrows = SCIPmatrixGetNRows(matrix);

   fill = 0;
   for( r = 0; r < nrows; r++ )
   {
      if( !SCIPmatrixIsRowRhsInfinity(matrix, r) )
         continue;

      baserows[fill] = r;
      fill++;
   }

   *nbaserows = fill;

   return SCIP_OKAY;
}

/** get bounds of variables */
static
void getBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   SCIP_Real*            lowerbds,           /**< lower bounds */
   SCIP_Real*            upperbds            /**< upper bounds */
   )
{
   int c;
   int ncols;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lowerbds != NULL);
   assert(upperbds != NULL);

   ncols = SCIPmatrixGetNColumns(matrix);

   for( c = 0; c < ncols; c++ )
   {
      SCIP_VAR* var;
      var = SCIPmatrixGetVar(matrix, c);
      lowerbds[c] = SCIPvarGetLbGlobal(var);
      upperbds[c] = SCIPvarGetUbGlobal(var);
   }
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyTworowbnd)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolTworowbnd(scip) );

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecTworowbnd)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPgetNVars(scip) == 0 || SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized && complete )
   {
      int* baserows;
      int nbaserows;
      int ntightenbnds;
      BNDCHGTYPE* tighten;
      int ndeletecons;
      SCIP_Bool* deletecons;
      int ncols;
      int nrows;
      SCIP_Real* lowerbds;
      SCIP_Real* upperbds;

      ncols = SCIPmatrixGetNColumns(matrix);
      nrows = SCIPmatrixGetNRows(matrix);

      ntightenbnds = 0;
      ndeletecons = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &baserows, nrows) );

      SCIP_CALL( SCIPallocBufferArray(scip, &tighten, ncols) );
      BMSclearMemoryArray(tighten, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &lowerbds, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &upperbds, ncols) );
      getBounds(scip, matrix, lowerbds, upperbds);

      SCIP_CALL( SCIPallocBufferArray(scip, &deletecons, nrows) );
      BMSclearMemoryArray(deletecons, nrows);

      SCIP_CALL( getBaseRows(scip, matrix, &nbaserows, baserows) );

      SCIP_CALL( calcTwoRowBnds(scip, matrix,
            nbaserows, baserows, lowerbds, upperbds,
            &ntightenbnds, tighten, &ndeletecons, deletecons) );

      if( ntightenbnds > 0 )
      {
         int c;
         SCIP_VAR* var;
         SCIP_Bool infeas;
         SCIP_Bool tightened;

         for( c = 0; c < ncols; c++ )
         {
            if( tighten[c] == LOWERBOUND )
            {
               var = SCIPmatrixGetVar(matrix, c);
               SCIP_CALL( SCIPtightenVarLb(scip, var, lowerbds[c], FALSE, &infeas, &tightened) );
               if( tightened )
                  ++(*nchgbds);
            }
            else if( tighten[c] == UPPERBOUND )
            {
               var = SCIPmatrixGetVar(matrix, c);
               SCIP_CALL( SCIPtightenVarUb(scip, var, upperbds[c], FALSE, &infeas, &tightened) );
               if( tightened )
                  ++(*nchgbds);
            }
            else if( tighten[c] == BOTHBOUNDS )
            {
               var = SCIPmatrixGetVar(matrix, c);
               SCIP_CALL( SCIPtightenVarLb(scip, var, lowerbds[c], FALSE, &infeas, &tightened) );
               if( tightened )
                  ++(*nchgbds);
               SCIP_CALL( SCIPtightenVarUb(scip, var, upperbds[c], FALSE, &infeas, &tightened) );
               if( tightened )
                  ++(*nchgbds);
            }
         }
      }

      if( ndeletecons > 0 )
      {
         int r;
         for( r = 0; r < nrows; r++ )
         {
            if( deletecons[r] )
            {
               SCIP_CONS* cons;
               cons = SCIPmatrixGetCons(matrix, r);
               SCIP_CALL( SCIPdelCons(scip, cons) );

               (*ndelconss)++;
            }
         }
      }

      SCIPfreeBufferArray(scip, &deletecons);
      SCIPfreeBufferArray(scip, &upperbds);
      SCIPfreeBufferArray(scip, &lowerbds);
      SCIPfreeBufferArray(scip, &tighten);
      SCIPfreeBufferArray(scip, &baserows);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the tworowbnd presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolTworowbnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecTworowbnd, NULL) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyTworowbnd) );

   return SCIP_OKAY;
}
