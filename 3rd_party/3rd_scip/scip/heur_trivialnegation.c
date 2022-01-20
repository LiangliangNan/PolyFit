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

/**@file   heur_trivialnegation.c
 * @brief  trivialnegation primal heuristic
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_trivialnegation.h"
#include "scip/pub_reopt.h"
#include "scip/struct_scip.h"

#define HEUR_NAME             "trivialnegation"
#define HEUR_DESC             "negate solution entries if an objective coefficient changes the sign, enters or leaves the objective."
#define HEUR_DISPCHAR         'j'
#define HEUR_PRIORITY         40000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

/*
 * Local methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTrivialnegation)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTrivialnegation(scip) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTrivialnegation)
{  /*lint --e{715}*/
   SCIP_SOL* lastbestsol;         /* best solution from last run */
   SCIP_SOL* allchanged;          /* solution with all entries negated */
   SCIP_SOL* feasiblechanged;     /* solution with all feasible entries negated */
   SCIP_SOL* singlenegatedsol;    /* solution with exactly one negated entry */
   SCIP_VAR** vars;
   int nbinvars;
   int i;

   SCIP_Real solval;

   vars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);

   *result = SCIP_DIDNOTRUN;

   if( !SCIPisReoptEnabled(scip) )
      return SCIP_OKAY;

   if( nbinvars < SCIPgetNVars(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get best solution from the run */
   lastbestsol = SCIPgetReoptLastOptSol(scip);

   if( lastbestsol == NULL )
      return SCIP_OKAY;

   /* initialize data structure */
   SCIP_CALL( SCIPcreateSol(scip, &allchanged, heur) );
   SCIP_CALL( SCIPcreateSol(scip, &feasiblechanged, heur) );
   SCIP_CALL( SCIPcreateSol(scip, &singlenegatedsol, heur) );

   /* copy the solutions */
   for( i = 0; i < nbinvars; i++ )
   {
      solval = SCIPgetSolVal(scip, lastbestsol, vars[i]);
      SCIP_CALL( SCIPsetSolVal(scip, allchanged, vars[i], solval) );
      SCIP_CALL( SCIPsetSolVal(scip, feasiblechanged, vars[i], solval) );
      SCIP_CALL( SCIPsetSolVal(scip, singlenegatedsol, vars[i], solval) );
   }

   assert(SCIPsolGetHeur(allchanged) == heur);
   assert(SCIPsolGetHeur(feasiblechanged) == heur);
   assert(SCIPsolGetHeur(singlenegatedsol) == heur);

   /* change the entries */
   for( i = 0; i < nbinvars; i++ )
   {
      SCIP_VAR* transvar;

      assert(SCIPvarIsActive(vars[i]));

      transvar = vars[i];

      if( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_BINARY )
      {
         SCIP_Real obj;
         SCIP_Real newcoef;
         SCIP_Real oldcoef;
         SCIP_Bool changed;

         /* perform negation only on variables that are not globally fixed */
         if( SCIPvarGetLbGlobal(vars[i]) > 0.5 || SCIPvarGetUbGlobal(vars[i]) < 0.5 )
            continue;

         SCIP_CALL( SCIPgetReoptOldObjCoef(scip, transvar, SCIPgetNReoptRuns(scip), &oldcoef));
         SCIP_CALL( SCIPgetReoptOldObjCoef(scip, transvar, SCIPgetNReoptRuns(scip)-1, &newcoef));

         /* check if variable entered or left the objective, or if its objective coefficient changed sign */
         changed = FALSE;
         if( !SCIPisFeasEQ(scip, oldcoef, newcoef) )
         {
            changed = SCIPisZero(scip, oldcoef) != SCIPisZero(scip, newcoef);
            changed |= SCIPisPositive(scip, oldcoef) == SCIPisNegative(scip, newcoef); /*lint !e514*/
         }

         SCIPdebugMsg(scip, "check variable <%s> which has %schanged from %g to %g\n", SCIPvarGetName(transvar), changed ? "" : "not ", oldcoef, newcoef);

         if( changed )
         {
            SCIP_Bool success;

            solval = SCIPgetSolVal(scip, lastbestsol, vars[i]);

            /* change solution value */
            SCIP_CALL( SCIPsetSolVal(scip, allchanged, vars[i], 1 - solval) );
            SCIP_CALL( SCIPsetSolVal(scip, feasiblechanged, vars[i], 1 - solval) );
            SCIP_CALL( SCIPsetSolVal(scip, singlenegatedsol, vars[i], 1 - solval) );

            /* try solution with all changes */
            success = FALSE;
            obj = SCIPgetSolTransObj(scip, allchanged);
            if( SCIPisFeasLT(scip, obj, SCIPgetCutoffbound(scip)) )
            {
               SCIPdebugMsg(scip, "try solution with all negations\n");
#ifdef SCIP_DEBUG
               SCIP_CALL( SCIPtrySol(scip, allchanged, TRUE, FALSE, TRUE, FALSE, TRUE, &success) );
#else
               SCIP_CALL( SCIPtrySol(scip, allchanged, FALSE, FALSE, TRUE, FALSE, TRUE, &success) );
#endif

               if( success )
               {
                  SCIPdebugMsg(scip, "found feasible solution solution:\n");
                  SCIPdebug( SCIP_CALL( SCIPprintSol(scip, allchanged, NULL, FALSE) ) );

                  *result = SCIP_FOUNDSOL;
               }
            }

            /* try solution with feasible changes */
            success = FALSE;
            obj = SCIPgetSolTransObj(scip, feasiblechanged);
            if( SCIPisFeasLT(scip, obj, SCIPgetCutoffbound(scip)) )
            {
               SCIPdebugMsg(scip, "try solution with feasible negations\n");
#ifdef SCIP_DEBUG
               SCIP_CALL( SCIPtrySol(scip, feasiblechanged, TRUE, FALSE, TRUE, FALSE, TRUE, &success) );
#else
               SCIP_CALL( SCIPtrySol(scip, feasiblechanged, FALSE, FALSE, TRUE, FALSE, TRUE, &success) );
#endif
               if( success )
               {
                  SCIPdebugMsg(scip, "found feasible solution solution:\n");
                  SCIPdebug( SCIP_CALL( SCIPprintSol(scip, feasiblechanged, NULL, FALSE) ) );

                  *result = SCIP_FOUNDSOL;
               }
            }

            if( !success )
            {
               /* reset solution with feasible changes */
               SCIP_CALL( SCIPsetSolVal(scip, feasiblechanged, vars[i], solval) );
            }

            /* try solution with exactly one changed value */
            obj = SCIPgetSolTransObj(scip, singlenegatedsol);
            if( SCIPisFeasLT(scip, obj, SCIPgetCutoffbound(scip)) )
            {
               success = FALSE;
               SCIPdebugMsg(scip, "try solution with a single negation\n");
#ifdef SCIP_DEBUG
               SCIP_CALL( SCIPtrySol(scip, singlenegatedsol, TRUE, FALSE, TRUE, FALSE, TRUE, &success) );
#else
               SCIP_CALL( SCIPtrySol(scip, singlenegatedsol, FALSE, FALSE, TRUE, FALSE, TRUE, &success) );
#endif
               if( success )
               {
                  SCIPdebugMsg(scip, "found feasible solution:\n");
                  SCIPdebug( SCIP_CALL( SCIPprintSol(scip, singlenegatedsol, NULL, FALSE) ) );

                  *result = SCIP_FOUNDSOL;
               }
            }

            /* reset solution with exactly one changed value */
            SCIP_CALL( SCIPsetSolVal(scip, singlenegatedsol, vars[i], solval) );
         }
      }
   }

   /* free solutions */
   SCIP_CALL( SCIPfreeSol(scip, &allchanged) );
   SCIP_CALL( SCIPfreeSol(scip, &feasiblechanged) );
   SCIP_CALL( SCIPfreeSol(scip, &singlenegatedsol) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the trivialnegation primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTrivialnegation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEUR* heur;

   /* include heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTrivialnegation, NULL) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTrivialnegation) );

   return SCIP_OKAY;
}
