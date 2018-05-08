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

/**@file   heur_mpec.c
 * @brief  mpec primal heuristic
 * @author Felipe Serrano
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_mpec.h"
#include "scip/heur_subnlp.h"
#include "nlpi/nlpi.h"


#define HEUR_NAME             "mpec"
#define HEUR_DESC             "regularization heuristic for convex and nonconvex MINLPs"
#define HEUR_DISPCHAR         'W'
#define HEUR_PRIORITY         -2050000
#define HEUR_FREQ             50
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE           /**< disable the heuristic in sub-SCIPs, even though it does not use any */

#define DEFAULT_INITTHETA     0.125          /**< default initial regularization right-hand side value (< 0.25) */
#define DEFAULT_SIGMA         0.5            /**< default regularization update factor (< 1) */
#define DEFAULT_MAXITER       100            /**< default maximum number of iterations of the MPEC loop */
#define DEFAULT_MAXNLPITER    500            /**< default maximum number of NLP iterations per solve */
#define DEFAULT_MINGAPLEFT    0.05           /**< default minimum amount of gap left in order to call the heuristic */
#define DEFAULT_SUBNLPTRIGGER 1e-3           /**< default maximum integrality violation before triggering a sub-NLP call */
#define DEFAULT_MAXNLPCOST    1e+8           /**< default maximum cost available for solving NLPs per call of the heuristic */
#define DEFAULT_MINIMPROVE    0.01           /**< default factor by which heuristic should at least improve the incumbent */
#define DEFAULT_MAXNUNSUCC    10             /**< default maximum number of consecutive calls for which the heuristic did not find an improving solution */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_NLPI*            nlpi;               /**< nlpi used to create the nlpi problem */
   SCIP_NLPIPROBLEM*     nlpiprob;           /**< nlpi problem representing the NLP relaxation */
   SCIP_HASHMAP*         var2idx;            /**< mapping between variables and nlpi indices */
   SCIP_HEUR*            subnlp;             /**< sub-NLP heuristic */

   SCIP_Real             inittheta;          /**< initial regularization right-hand side value */
   SCIP_Real             sigma;              /**< regularization update factor */
   SCIP_Real             subnlptrigger;      /**< maximum number of NLP iterations per solve */
   SCIP_Real             maxnlpcost;         /**< maximum cost available for solving NLPs per call of the heuristic */
   SCIP_Real             minimprove;         /**< factor by which heuristic should at least improve the incumbent */
   SCIP_Real             mingapleft;         /**< minimum amount of gap left in order to call the heuristic */
   int                   maxiter;            /**< maximum number of iterations of the MPEC loop */
   int                   maxnlpiter;         /**< maximum number of NLP iterations per solve */
   int                   nunsucc;             /**< number of consecutive calls for which the heuristic did not find an
                                              * improving solution */
   int                   maxnunsucc;         /**< maximum number of consecutive calls for which the heuristic did not
                                              * find an improving solution */
};


/*
 * Local methods
 */

/** creates the data structure for generating the current NLP relaxation */
static
SCIP_RETCODE createNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   SCIP_Real cutoff = SCIPinfinity(scip);

   assert(heurdata != NULL);
   assert(heurdata->nlpi != NULL);

   /* NLP has been already created */
   if( heurdata->nlpiprob != NULL )
      return SCIP_OKAY;

   /* compute cutoff value to ensure minimum improvement */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      assert( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) );

      if( !SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) )
      {
         cutoff = (1.0 - heurdata->minimprove) * SCIPgetUpperbound(scip)
            + heurdata->minimprove * SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound(scip) >= 0.0 )
            cutoff = ( 1.0 - heurdata->minimprove ) * SCIPgetUpperbound(scip);
         else
            cutoff = ( 1.0 + heurdata->minimprove ) * SCIPgetUpperbound(scip);
      }
      cutoff = MIN(upperbound, cutoff);
      SCIPdebugMsg(scip, "set objective limit %g in [%g,%g]\n", cutoff, SCIPgetLowerbound(scip),
         SCIPgetUpperbound(scip));
   }

   SCIP_CALL( SCIPnlpiCreateProblem(heurdata->nlpi, &heurdata->nlpiprob, "MPEC-nlp") );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->var2idx, SCIPblkmem(scip), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPcreateNlpiProb(scip, heurdata->nlpi, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip),
         heurdata->nlpiprob, heurdata->var2idx, NULL, cutoff, TRUE, FALSE) );

   return SCIP_OKAY;
}

/** frees the data structures for the NLP relaxation */
static
SCIP_RETCODE freeNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   assert(heurdata != NULL);

   /* NLP has not been created yet */
   if( heurdata->nlpiprob == NULL )
      return SCIP_OKAY;

   assert(heurdata->nlpi != NULL);
   assert(heurdata->var2idx != NULL);

   SCIPhashmapFree(&heurdata->var2idx);
   SCIP_CALL( SCIPnlpiFreeProblem(heurdata->nlpi, &heurdata->nlpiprob) );

   return SCIP_OKAY;
}

/** add or updates the regularization constraints to the NLP; for a given parameter theta we add for each non-fixed
 *  binary variable z the constraint z*(1-z) <= theta; if these constraint are already present we update the theta on
 *  the right-hand side
 */
static
SCIP_RETCODE addRegularScholtes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VAR**            binvars,            /**< array containing all non-fixed binary variables */
   int                   nbinvars,           /**< total number of non-fixed binary variables */
   SCIP_Real             theta,              /**< regularization parameter */
   SCIP_Bool             update              /**< should the regularization constraints be added or updated? */
   )
{
   int i;

   assert(binvars != NULL);
   assert(nbinvars > 0);

   /* add or update regularization for each non-fixed binary variables */
   if( !update )
   {
      SCIP_QUADELEM* quadelems;
      SCIP_Real* linvals;
      int* lininds;

      SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds, 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals, 1) );

      for( i = 0; i < nbinvars; ++i )
      {
         SCIP_VAR* var = binvars[i];
         SCIP_Real lhs = -SCIPinfinity(scip);
         SCIP_Real rhs = theta;
         int nlininds = 1;
         int nquadelems = 1;
         int idx;

         assert(var != NULL);
         assert(heurdata->var2idx != NULL);
         assert(SCIPhashmapExists(heurdata->var2idx, (void*)var));
         idx = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)var);

         lininds[0] = idx;
         linvals[0] = 1.0;
         quadelems->idx1 = lininds[0];
         quadelems->idx2 = lininds[0];
         quadelems->coef = -1.0;

         SCIP_CALL( SCIPnlpiAddConstraints(heurdata->nlpi, heurdata->nlpiprob, 1, &lhs, &rhs, &nlininds,
               &lininds, &linvals, &nquadelems, &quadelems, NULL, NULL, NULL) );
      }

      SCIPfreeBufferArray(scip, &linvals);
      SCIPfreeBufferArray(scip, &lininds);
      SCIPfreeBufferArray(scip, &quadelems);
   }
   else
   {
      int startidx = SCIPgetNNLPNlRows(scip) + 1; /* the cutoff is a separate constraint */
      SCIP_Real* lhss;
      SCIP_Real* rhss;
      int* indices;

      SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nbinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nbinvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &indices, nbinvars) );

      for( i = 0; i < nbinvars; ++i )
      {
         lhss[i] = -SCIPinfinity(scip);
         rhss[i] = theta;
         indices[i] = startidx + i;
      }

      SCIP_CALL( SCIPnlpiChgConsSides(heurdata->nlpi, heurdata->nlpiprob, nbinvars, indices, lhss, rhss) );

      SCIPfreeBufferArray(scip, &indices);
      SCIPfreeBufferArray(scip, &rhss);
      SCIPfreeBufferArray(scip, &lhss);
   }

   return SCIP_OKAY;
}

/** recursive helper function to count the number of nodes in a sub-tree */
static
int getExprSize(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int sum;
   int i;

   assert(expr != NULL);

   sum = 0;
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
      sum += getExprSize(child);
   }
   return 1 + sum;
}

/** returns the number of nodes in an expression tree */
static
int getExprtreeSize(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   if( tree == NULL )
      return 0;
   return getExprSize(SCIPexprtreeGetRoot(tree));
}

/** returns the available time limit that is left */
static
SCIP_RETCODE getTimeLeft(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            timeleft            /**< pointer to store the remaining time limit */
   )
{
   SCIP_Real timelim;

   assert(timeleft != NULL);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelim) );

   if( SCIPisInfinity(scip, timelim) )
      *timeleft = timelim - SCIPgetSolvingTime(scip);
   else
      *timeleft = SCIPinfinity(scip);

   return SCIP_OKAY;
}

/** main execution function of the MPEC heuristic */
static
SCIP_RETCODE heurExec(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< MPEC heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_NLPSTATISTICS* nlpstatistics = NULL;
   SCIP_VAR** binvars = NULL;
   SCIP_Real* initguess = NULL;
   SCIP_Real* ubs = NULL;
   SCIP_Real* lbs = NULL;
   int* indices = NULL;
   SCIP_Real theta = heurdata->inittheta;
   SCIP_Real nlpcostperiter = 0.0;
   SCIP_Real nlpcostleft = heurdata->maxnlpcost;
   SCIP_Bool reinit = TRUE;
   SCIP_Bool fixed = FALSE;
   SCIP_Bool subnlpcalled = FALSE;
   int nbinvars = 0;
   int i;

   assert(heurdata->nlpiprob != NULL);
   assert(heurdata->var2idx != NULL);
   assert(heurdata->nlpi != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, SCIPgetNBinVars(scip)) );
   SCIP_CALL( SCIPnlpStatisticsCreate(SCIPblkmem(scip), &nlpstatistics) );

   /* collect all non-fixed binary variables */
   for( i = 0; i < SCIPgetNBinVars(scip); ++i )
   {
      SCIP_VAR* var = SCIPgetVars(scip)[i];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      if( !SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         binvars[nbinvars++] = var;
   }

   /* all binary variables are fixed */
   SCIPdebugMsg(scip, "nbinvars %d\n", nbinvars);
   if( nbinvars == 0 )
      goto TERMINATE;

   SCIP_CALL( SCIPallocBufferArray(scip, &initguess, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indices, nbinvars) );

   /* compute estimate cost for each NLP iteration */
   for( i = 0; i < SCIPgetNNLPNlRows(scip); ++i )
   {
      SCIP_NLROW* nlrow = SCIPgetNLPNlRows(scip)[i];
      assert(nlrow != NULL);

      nlpcostperiter += 1.0 * SCIPnlrowGetNLinearVars(nlrow)
               + 2.0 * SCIPnlrowGetNQuadElems(nlrow)
               + 3.0 * getExprtreeSize(SCIPnlrowGetExprtree(nlrow));
   }

   /* set initial guess */
   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      SCIP_VAR* var = SCIPgetVars(scip)[i];
      initguess[i] = SCIPgetSolVal(scip, NULL, var);
      /* SCIPdebugMsg(scip, "set initial value for %s to %g\n", SCIPvarGetName(var), initguess[i]); */
   }
   SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, initguess, NULL, NULL, NULL) );

   /* set parameters of NLP solver */
   SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_FEASTOL,
         SCIPfeastol(scip) / 10.0) );
   SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_RELOBJTOL,
         SCIPdualfeastol(scip) / 10.0) );
   SCIP_CALL( SCIPnlpiSetIntPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_VERBLEVEL, 0) );

   /* main loop */
   for( i = 0; i < heurdata->maxiter && *result != SCIP_FOUNDSOL && nlpcostleft > 0.0 && !SCIPisStopped(scip); ++i )
   {
      SCIP_Real* primal = NULL;
      SCIP_Real timeleft = SCIPinfinity(scip);
      SCIP_Bool binaryfeasible;
      SCIP_Bool regularfeasible;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Real maxviolbin = 0.0;
      SCIP_Real maxviolreg = 0.0;
      int j;

      /* add or update regularization */
      SCIP_CALL( addRegularScholtes(scip, heurdata, binvars, nbinvars, theta, i > 0) );

      /* set working limits */
      SCIP_CALL( getTimeLeft(scip, &timeleft) );
      if( timeleft <= 0.0 )
      {
         SCIPdebugMsg(scip, "skip NLP solve; no time left\n");
         break;
      }

      SCIP_CALL( SCIPnlpiSetRealPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_TILIM, timeleft) );
      SCIP_CALL( SCIPnlpiSetIntPar(heurdata->nlpi, heurdata->nlpiprob, SCIP_NLPPAR_ITLIM, heurdata->maxnlpiter) );

      /* solve NLP */
      SCIP_CALL( SCIPnlpiSolve(heurdata->nlpi, heurdata->nlpiprob) );
      solstat = SCIPnlpiGetSolstat(heurdata->nlpi, heurdata->nlpiprob);

      /* give up if an error occurred or no primal values are accessible */
      if( solstat > SCIP_NLPSOLSTAT_LOCINFEASIBLE )
      {
         SCIPdebugMsg(scip, "error occured during NLP solve -> stop!\n");
         break;
      }

      /* update nlpcostleft */
      SCIP_CALL( SCIPnlpiGetStatistics(heurdata->nlpi, heurdata->nlpiprob, nlpstatistics) );
      nlpcostleft -= SCIPnlpStatisticsGetNIterations(nlpstatistics) * nlpcostperiter * nbinvars;
      SCIPdebugMsg(scip, "nlpcostleft = %e\n", nlpcostleft);

      SCIP_CALL( SCIPnlpiGetSolution(heurdata->nlpi, heurdata->nlpiprob, &primal, NULL, NULL, NULL, NULL) );
      assert(primal != NULL);

      /* check for binary feasibility */
      binaryfeasible = TRUE;
      regularfeasible = TRUE;
      for( j = 0; j < nbinvars; ++j )
      {
         int idx = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)binvars[j]);
         binaryfeasible = binaryfeasible && SCIPisFeasIntegral(scip, primal[idx]);
         regularfeasible = regularfeasible && SCIPisLE(scip, primal[idx] - SQR(primal[idx]), theta);

         maxviolreg = MAX(maxviolreg, primal[idx] - SQR(primal[idx]) - theta);
         maxviolbin = MAX(maxviolbin, MIN(primal[idx], 1.0-primal[idx]));
      }
      SCIPdebugMsg(scip, "maxviol-regularization %g maxviol-integrality %g\n", maxviolreg, maxviolbin);

      /* call sub-NLP heuristic when the maximum binary infeasibility is small enough (or this is the last iteration
       * because we reached the nlpcost limit)
       */
      if( !subnlpcalled && heurdata->subnlp != NULL
         && (SCIPisLE(scip, maxviolbin, heurdata->subnlptrigger) || nlpcostleft <= 0.0)
         && !SCIPisStopped(scip) )
      {
         SCIP_SOL* refpoint;
         SCIP_RESULT subnlpresult;

         SCIPdebugMsg(scip, "call sub-NLP heuristic because binary infeasibility is small enough\n");
         SCIP_CALL( SCIPcreateSol(scip, &refpoint, heur) );

         for( j = 0; j < SCIPgetNVars(scip); ++j )
         {
            SCIP_VAR* var = SCIPgetVars(scip)[j];
            SCIP_Real val = SCIPvarIsBinary(var) ? SCIPfeasRound(scip, primal[j]) : primal[j];
            SCIP_CALL( SCIPsetSolVal(scip, refpoint, var, val) );
         }

         SCIP_CALL( getTimeLeft(scip, &timeleft) );
         SCIP_CALL( SCIPapplyHeurSubNlp(scip, heurdata->subnlp, &subnlpresult, refpoint, -1LL, timeleft,
               heurdata->minimprove, NULL,
               NULL) );
         SCIP_CALL( SCIPfreeSol(scip, &refpoint) );
         SCIPdebugMsg(scip, "result of sub-NLP call: %d\n", subnlpresult);

         /* stop MPEC heuristic when the sub-NLP heuristic has found a feasible solution */
         if( subnlpresult == SCIP_FOUNDSOL )
         {
            SCIPdebugMsg(scip, "sub-NLP found a feasible solution -> stop!\n");
            break;
         }

         subnlpcalled = TRUE;
      }

      /* NLP feasible + binary feasible -> add solution and stop */
      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE && binaryfeasible )
      {
         SCIP_SOL* sol;
         SCIP_Bool stored;

         SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

         for( j = 0; j < SCIPgetNVars(scip); ++j )
         {
            SCIP_VAR* var = SCIPgetVars(scip)[j];
            assert(j == (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)var));
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, primal[j]) );
         }

#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, FALSE, &stored) );
#else
         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, TRUE, FALSE, &stored) );
#endif
         SCIPdebugMsg(scip, "found a solution (stored = %u)\n", stored);

         if( stored )
            *result = SCIP_FOUNDSOL;
         break;
      }

      /* NLP feasible + binary infeasible -> reduce theta */
      else if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE && !binaryfeasible )
      {
         BMScopyMemoryArray(initguess, primal, SCIPgetNVars(scip));
         SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, primal, NULL, NULL, NULL) );
         SCIPdebugMsg(scip, "update theta from %g -> %g\n", theta, theta*heurdata->sigma);

         if( !reinit )
         {
            SCIPdebugMsg(scip, "reinit fixed the infeasibility\n");
            reinit = TRUE;
         }

         theta *= heurdata->sigma;

         /* unfix binary variables */
         if( fixed )
         {
            SCIPdebugMsg(scip, "unfixing binary variables\n");
            for( j = 0; j < nbinvars; ++j )
            {
               lbs[j] = 0.0;
               ubs[j] = 1.0;
               indices[j] = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)binvars[j]);
            }
            SCIP_CALL( SCIPnlpiChgVarBounds(heurdata->nlpi, heurdata->nlpiprob, nbinvars, indices, lbs, ubs) );
            fixed = FALSE;
         }
      }

      /* NLP infeasible + regularization feasible -> stop (give up) */
      else if( solstat > SCIP_NLPSOLSTAT_FEASIBLE && regularfeasible )
      {
         SCIPdebugMsg(scip, "NLP is infeasible but regularization constraints are satisfied -> stop!\n");
         break;
      }

      /* NLP infeasible + binary infeasible -> set initial point / fix binary variables */
      else
      {
         assert(solstat > SCIP_NLPSOLSTAT_FEASIBLE && !regularfeasible);

         SCIPdebugMsg(scip, "NLP solution is not feasible for the NLP and the binary variables\n");

         /* stop if fixing did not resolve the infeasibility */
         if( fixed )
         {
            SCIPdebugMsg(scip, "fixing variables did not resolve infeasibility -> stop!\n");
            break;
         }

         /* fix variables if reinit is FALSE; otherwise set another initial point */
         if( !reinit )
         {
            int nfixedvars = 0;

            /* fix binary variables */
            for( j = 0; j < nbinvars; ++j )
            {
               int idx = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)binvars[j]);
               indices[j] = idx;

               if( SCIPisFeasLE(scip, primal[idx] - SQR(primal[idx]), theta) )
               {
                  lbs[j] = 0.0;
                  ubs[j] = 1.0;
               }
               else
               {
                  lbs[j] = primal[idx] >= 0.5 ? 0.0 : 1.0;
                  ubs[j] = primal[idx] >= 0.5 ? 0.0 : 1.0;
                  ++nfixedvars;
                  /* SCIPdebugMsg(scip, "fix binary variable %s = %g\n", SCIPvarGetName(binvars[j]), ubs[j]); */
               }
            }
            SCIPdebugMsg(scip, "fixed %d binary variables\n", nfixedvars);
            SCIP_CALL( SCIPnlpiChgVarBounds(heurdata->nlpi, heurdata->nlpiprob, nbinvars, indices, lbs, ubs) );
            fixed = TRUE;
         }
         else
         {
            SCIPdebugMsg(scip, "update initial guess\n");

            /* set initial point */
            for( j = 0; j < nbinvars; ++j )
            {
               int idx = (int)(size_t)SCIPhashmapGetImage(heurdata->var2idx, (void*)binvars[j]);
               initguess[idx] = primal[idx] >= 0.5 ? 0.0 : 1.0;
               /* SCIPdebugMsg(scip, "update init guess for %s to %g\n", SCIPvarGetName(binvars[j]), initguess[idx]); */
            }
            SCIP_CALL( SCIPnlpiSetInitialGuess(heurdata->nlpi, heurdata->nlpiprob, initguess, NULL, NULL, NULL) );
            reinit = FALSE;
         }
      }
   }

TERMINATE:
   SCIPfreeBufferArrayNull(scip, &indices);
   SCIPfreeBufferArrayNull(scip, &ubs);
   SCIPfreeBufferArrayNull(scip, &lbs);
   SCIPfreeBufferArrayNull(scip, &initguess);
   SCIPnlpStatisticsFree(SCIPblkmem(scip), &nlpstatistics);
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyMpec)
{  /*lint --e{715}*/
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurMpec(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);
   assert(heurdata->nlpi == NULL);

   if( SCIPgetNNlpis(scip) > 0 )
   {
      heurdata->nlpi = SCIPgetNlpis(scip)[0];
      heurdata->subnlp = SCIPfindHeur(scip, "subnlp");
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);
   heurdata->nlpi = NULL;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMpec)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   SCIP_CONSHDLR* andhdlr = SCIPfindConshdlr(scip, "and");
   SCIP_CONSHDLR* sosonehdlr = SCIPfindConshdlr(scip, "SOS1");
   SCIP_CONSHDLR* sostwohdlr = SCIPfindConshdlr(scip, "SOS2");

   assert(heurdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNIntVars(scip) > 0 || SCIPgetNBinVars(scip) == 0
      || heurdata->nlpi == NULL || !SCIPisNLPConstructed(scip)
      || heurdata->mingapleft > SCIPgetGap(scip)
      || heurdata->nunsucc > heurdata->maxnunsucc )
      return SCIP_OKAY;

   /* skip heuristic if constraints without a nonlinear representation are present */
   if( (andhdlr != NULL && SCIPconshdlrGetNConss(andhdlr) > 0) ||
      (sosonehdlr != NULL && SCIPconshdlrGetNConss(sosonehdlr) > 0) ||
      (sostwohdlr != NULL && SCIPconshdlrGetNConss(sostwohdlr) > 0) )
   {
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* call MPEC method */
   SCIP_CALL( createNLP(scip, heurdata) );
   SCIP_CALL( heurExec(scip, heur, heurdata, result) );
   SCIP_CALL( freeNLP(scip, heurdata) );

   /* update number of unsuccessful calls */
   heurdata->nunsucc = (*result == SCIP_FOUNDSOL) ? 0 : heurdata->nunsucc + 1;

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the mpec primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMpec(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata = NULL;
   SCIP_HEUR* heur = NULL;

   /* create mpec primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMpec, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyMpec) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMpec) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolMpec) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolMpec) );

   /* add mpec primal heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/inittheta",
         "initial regularization right-hand side value",
         &heurdata->inittheta, FALSE, DEFAULT_INITTHETA, 0.0, 0.25, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/sigma",
         "regularization update factor",
         &heurdata->sigma, FALSE, DEFAULT_SIGMA, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/subnlptrigger",
         "maximum number of NLP iterations per solve",
         &heurdata->subnlptrigger, FALSE, DEFAULT_SUBNLPTRIGGER, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxnlpcost",
         "maximum cost available for solving NLPs per call of the heuristic",
         &heurdata->maxnlpcost, FALSE, DEFAULT_MAXNLPCOST, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which heuristic should at least improve the incumbent",
         &heurdata->minimprove, FALSE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/mingapleft",
         "minimum amount of gap left in order to call the heuristic",
         &heurdata->mingapleft, FALSE, DEFAULT_MINGAPLEFT, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxiter",
         "maximum number of iterations of the MPEC loop",
         &heurdata->maxiter, FALSE, DEFAULT_MAXITER, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxnlpiter",
         "maximum number of NLP iterations per solve",
         &heurdata->maxnlpiter, FALSE, DEFAULT_MAXNLPITER, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxnunsucc",
         "maximum number of consecutive calls for which the heuristic did not find an improving solution",
         &heurdata->maxnunsucc, FALSE, DEFAULT_MAXNUNSUCC, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
