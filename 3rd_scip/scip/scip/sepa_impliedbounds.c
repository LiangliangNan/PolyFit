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

/**@file   sepa_impliedbounds.c
 * @brief  implied bounds separator
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_impliedbounds.h"
#include "scip/pub_misc.h"


#define SEPA_NAME              "impliedbounds"
#define SEPA_DESC              "implied bounds separator"
#define SEPA_PRIORITY               -50
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define RELCUTCOEFMAXRANGE          1.0 /**< maximal allowed range of cut coefficients, relative to 1/feastol */
#define DEFAULT_USETWOSIZECLIQUES  TRUE /**< should violated inequalities for cliques with 2 variables be separated? */

/** separator-specific data for the implied bounds separator */
struct SCIP_SepaData
{
   SCIP_Bool             usetwosizecliques;  /**< should violated inequalities for cliques with 2 variables be separated? */
};

/*
 * Local methods
 */

/** adds given cut with two variables, if it is violated */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_Real             val1,               /**< given coefficient of first variable */
   SCIP_VAR*             var1,               /**< given first variable */
   SCIP_Real             solval1,            /**< current LP solution value of first variable */
   SCIP_Real             val2,               /**< given coefficient of second variable */
   SCIP_VAR*             var2,               /**< given second variable */
   SCIP_Real             solval2,            /**< current LP solution value of second variable */
   SCIP_Real             rhs,                /**< given right hand side of the cut to add */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to update number of cuts added */
   )
{
   SCIP_Real activity;

   assert(ncuts != NULL);
   assert(cutoff != NULL);
   *cutoff = FALSE;

   /* calculate activity of cut */
   activity = val1 * solval1 + val2 * solval2;
   /*SCIPdebugMsg(scip, " -> %g<%s>[%g] + %g<%s>[%g] <= %g (act: %g)\n",
     val1, SCIPvarGetName(var1), solval1, val2, SCIPvarGetName(var2), solval2, rhs, activity);*/

   /* check, if cut is violated */
   if( SCIPisEfficacious(scip, activity - rhs) )
   {
      SCIP_ROW* cut;
      char cutname[SCIP_MAXSTRLEN];

      /* create cut */
      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "implbd%d_%d", SCIPgetNLPs(scip), *ncuts);
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
      SCIP_CALL( SCIPaddVarToRow(scip, cut, var1, val1) );
      SCIP_CALL( SCIPaddVarToRow(scip, cut, var2, val2) );
      SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
      /* set cut rank: for implied bounds we always set to 1 */
      SCIProwChgRank(cut, 1);

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, " -> found cut (activity = %g): ", activity);
      SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

      SCIP_CALL( SCIPaddPoolCut(scip, cut) );
      (*ncuts)++;

      /* release cut */
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }

   return SCIP_OKAY;
}

/** searches and adds implied bound cuts that are violated by the given solution value array */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            solvals,            /**< array with solution values of all problem variables */
   SCIP_VAR**            fracvars,           /**< array of fractional variables */
   SCIP_Real*            fracvals,           /**< solution values of fractional variables */
   int                   nfracs,             /**< number of fractional variables */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to store the number of generated cuts */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_SEPADATA* sepadata;
   int ncliques;
   int i;

   assert(solvals != NULL);
   assert(fracvars != NULL || nfracs == 0);
   assert(fracvals != NULL || nfracs == 0);
   assert(cutoff != NULL);
   assert(ncuts != NULL);

   *cutoff = FALSE;
   *ncuts = 0;
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPdebugMsg(scip, "searching for implied bound cuts\n");

   /* search binary variables for violated implications */
   for( i = 0; i < nfracs; i++ )
   {
      SCIP_BOUNDTYPE* impltypes; 
      SCIP_Real* implbounds; 
      SCIP_VAR** implvars;
      int nimpl;
      int j;

      assert(fracvars != NULL);
      assert(fracvals != NULL);

      /* only process binary variables */
      if( SCIPvarGetType(fracvars[i]) != SCIP_VARTYPE_BINARY )
         continue;

      /* get implications of x == 1 */
      nimpl = SCIPvarGetNImpls(fracvars[i], TRUE);
      implvars = SCIPvarGetImplVars(fracvars[i], TRUE);
      impltypes = SCIPvarGetImplTypes(fracvars[i], TRUE);
      implbounds = SCIPvarGetImplBounds(fracvars[i], TRUE);

      /*SCIPdebugMsg(scip, "%d implications for <%s>[%g] == 1\n", nimpl, SCIPvarGetName(fracvars[i]), fracvals[i]);*/

      /* try to add cuts for implications of x == 1
       *    x == 1 -> y <= p:  y <= ub + x * (p - ub)  <==>  y + (ub - p) * x <=  ub
       *    x == 1 -> y >= p:  y >= lb + x * (p - lb)  <==> -y + (p - lb) * x <= -lb
       * with lb (ub) global lower (upper) bound of y
       */
      for( j = 0; j < nimpl; j++ )
      {
         SCIP_Real solval;

         assert(implvars != NULL);
         assert(impltypes != NULL);
         assert(implbounds != NULL);

         /* consider only implications with active implvar */
         if( SCIPvarGetProbindex(implvars[j]) < 0 )
            continue;

         solval = solvals[SCIPvarGetProbindex(implvars[j])];
         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            SCIP_Real ub;

            /* implication x == 1 -> y <= p */
            ub = SCIPvarGetUbGlobal(implvars[j]);

            /* consider only nonredundant and numerical harmless implications */
            if( SCIPisLE(scip, implbounds[j], ub) && (ub - implbounds[j]) * SCIPfeastol(scip) <= RELCUTCOEFMAXRANGE )
            {
               /* add cut if violated */
               SCIP_CALL( addCut(scip, sepa, 1.0, implvars[j], solval, (ub - implbounds[j]), fracvars[i], fracvals[i],
                     ub, cutoff, ncuts) );
               if ( *cutoff )
                  return SCIP_OKAY;
            }
         }
         else
         {
            SCIP_Real lb;

            /* implication x == 1 -> y >= p */
            lb = SCIPvarGetLbGlobal(implvars[j]);
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);

            /* consider only nonredundant and numerical harmless implications */
            if( SCIPisGE(scip, implbounds[j], lb) && (implbounds[j] - lb) * SCIPfeastol(scip) <= RELCUTCOEFMAXRANGE )
            {
               /* add cut if violated */
               SCIP_CALL( addCut(scip, sepa, -1.0, implvars[j], solval, (implbounds[j] - lb), fracvars[i], fracvals[i],
                     -lb, cutoff, ncuts) );
               if ( *cutoff )
                  return SCIP_OKAY;
            }
         }
      }

      /* get implications of x == 0 */
      nimpl = SCIPvarGetNImpls(fracvars[i], FALSE);
      implvars = SCIPvarGetImplVars(fracvars[i], FALSE);
      impltypes = SCIPvarGetImplTypes(fracvars[i], FALSE);
      implbounds = SCIPvarGetImplBounds(fracvars[i], FALSE);

      /*SCIPdebugMsg(scip, "%d implications for <%s>[%g] == 0\n", nimpl, SCIPvarGetName(fracvars[i]), fracvals[i]);*/

      /* try to add cuts for implications of x == 0
       *    x == 0 -> y <= p:  y <= p + x * (ub - p)  <==>  y + (p - ub) * x <=  p
       *    x == 0 -> y >= p:  y >= p + x * (lb - p)  <==> -y + (lb - p) * x <= -p
       * with lb (ub) global lower (upper) bound of y
       */
      for( j = 0; j < nimpl; j++ )
      {
         SCIP_Real solval;

         /* consider only implications with active implvar */
         if( SCIPvarGetProbindex(implvars[j]) < 0 )
            continue;

         solval = solvals[SCIPvarGetProbindex(implvars[j])];
         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            SCIP_Real ub;

            /* implication x == 0 -> y <= p */
            ub = SCIPvarGetUbGlobal(implvars[j]);

            /* consider only nonredundant and numerical harmless implications */
            if( SCIPisLE(scip, implbounds[j], ub) && (ub - implbounds[j]) * SCIPfeastol(scip) < RELCUTCOEFMAXRANGE )
            {
               /* add cut if violated */
               SCIP_CALL( addCut(scip, sepa, 1.0, implvars[j], solval, (implbounds[j] - ub), fracvars[i], fracvals[i],
                     implbounds[j], cutoff, ncuts) );
               if ( *cutoff )
                  return SCIP_OKAY;
            }
         }
         else
         {
            SCIP_Real lb;

            /* implication x == 0 -> y >= p */
            lb = SCIPvarGetLbGlobal(implvars[j]);
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);

            /* consider only nonredundant and numerical harmless implications */
            if( SCIPisGE(scip, implbounds[j], lb) && (implbounds[j] - lb) * SCIPfeastol(scip) < RELCUTCOEFMAXRANGE )
            {
               /* add cut if violated */
               SCIP_CALL( addCut(scip, sepa, -1.0, implvars[j], solval, (lb - implbounds[j]), fracvars[i], fracvals[i],
                     -implbounds[j], cutoff, ncuts) );
               if ( *cutoff )
                  return SCIP_OKAY;
            }
         }
      }
   }

   /* stop separation here if cliques should not be separated */
   if( ! sepadata->usetwosizecliques )
      return SCIP_OKAY;

   /* prepare clean clique data */
   SCIP_CALL( SCIPcleanupCliques(scip, cutoff) );

   if( *cutoff )
      return SCIP_OKAY;

   cliques = SCIPgetCliques(scip);
   ncliques = SCIPgetNCliques(scip);

   /* loop over cliques of size 2 which are essentially implications and add cuts if they are violated */
   for( i = 0; i < ncliques; ++i )
   {
      SCIP_CLIQUE* clique;
      SCIP_VAR** clqvars;
      SCIP_Bool* clqvals;
      SCIP_Real rhs;

      clique = cliques[i];
      /* only consider inequality cliques of size 2 */
      if( SCIPcliqueGetNVars(clique) != 2 || SCIPcliqueIsEquation(clique) )
         continue;

      /* get variables and values of the clique */
      clqvars = SCIPcliqueGetVars(clique);
      clqvals = SCIPcliqueGetValues(clique);

      /* clique variables should never be equal after clean up */
      assert(clqvars[0] != clqvars[1]);

      /* calculate right hand side of clique inequality, which is initially 1 and decreased by 1 for every occurence of
       * a negated variable in the clique
       */
      rhs = 1.0;
      if( ! clqvals[0] )
         rhs -= 1.0;
      if( ! clqvals[1] )
         rhs -= 1.0;

      /* Basic clique inequality is
       *
       *       cx * x + (1-cx) (1-x) + cy * y + (1-cy) * (1-y) <= 1,
       *
       * where x and y are the two binary variables in the clique and cx and cy are their clique values, where a
       * clique value of 0 means that the negation of the variable should be part of the inequality.
       * Hence, exactly one of the two possible terms for x and y has a nonzero coefficient
       */
      SCIP_CALL( addCut(scip, sepa,
            clqvals[0] ? 1.0 : -1.0, clqvars[0], SCIPgetSolVal(scip, sol, clqvars[0]),
            clqvals[1] ? 1.0 : -1.0, clqvars[1], SCIPgetSolVal(scip, sol, clqvars[1]),
            rhs, cutoff, ncuts) );

      /* terminate if cutoff was found */
      if( *cutoff )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyImpliedbounds)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaImpliedbounds(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeImpliedbounds)
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
SCIP_DECL_SEPAEXECLP(sepaExeclpImpliedbounds)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR** fracvars;
   SCIP_Real* solvals;
   SCIP_Real* fracvals;
   SCIP_Bool cutoff;
   int nvars;
   int nbinvars;
   int nfracs;
   int ncuts;

   assert(sepa != NULL);
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;

   /* gets active problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* get fractional problem variables */
   /* todo try out also separating fractional implicit integer variables */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &fracvars, &fracvals, NULL, &nfracs, NULL, NULL) );
   if( nfracs == 0 )
      return SCIP_OKAY;

   /* get solution values for all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPgetVarSols(scip, nvars, vars, solvals) );

   /* call the cut separation */
   SCIP_CALL( separateCuts(scip, sepa, NULL, solvals, fracvars, fracvals, nfracs, &cutoff, &ncuts) );

   /* adjust result code */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolImpliedbounds)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR** fracvars;
   SCIP_Real* solvals;
   SCIP_Real* fracvals;
   SCIP_Bool cutoff;
   int nvars;
   int nbinvars;
   int nfracs;
   int ncuts;
   int i;

   assert(sepa != NULL);
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;

   /* gets active problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* get solution values for all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, solvals) );

   /* get binary problem variables that are fractional in given solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &fracvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fracvals, nbinvars) );
   nfracs = 0;
   for( i = 0; i < nbinvars; ++i )
   {
      if( !SCIPisFeasIntegral(scip, solvals[i]) )
      {
         fracvars[nfracs] = vars[i];
         fracvals[nfracs] = solvals[i];
         nfracs++;
      }
   }

   /* call the cut separation */
   ncuts = 0;
   cutoff = FALSE;

   if( nfracs > 0 )
   {
      SCIP_CALL( separateCuts(scip, sepa, sol, solvals, fracvars, fracvals, nfracs, &cutoff, &ncuts) );
   }

   /* adjust result code */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &fracvals);
   SCIPfreeBufferArray(scip, &fracvars);
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the impliedbounds separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaImpliedbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create impliedbounds separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   assert(sepadata != NULL);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpImpliedbounds, sepaExecsolImpliedbounds,
         sepadata) );
   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyImpliedbounds) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeImpliedbounds) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/impliedbounds/usetwosizecliques",
         "should violated inequalities for cliques with 2 variables be separated?",
         &sepadata->usetwosizecliques, TRUE, DEFAULT_USETWOSIZECLIQUES, NULL, NULL) );

   return SCIP_OKAY;
}
