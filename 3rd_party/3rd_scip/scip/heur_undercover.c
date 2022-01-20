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

/**@file   heur_undercover.c
 * @brief  Undercover primal heuristic for MINLPs
 * @author Timo Berthold
 * @author Ambros Gleixner
 *
 * The undercover heuristic is designed for mixed-integer nonlinear programs and tries to fix a subset of variables such
 * as to make each constraint linear or convex. For this purpose it solves a binary program to automatically determine
 * the minimum number of variable fixings necessary. As fixing values, we use values from the LP relaxation, the NLP
 * relaxation, or the incumbent solution.
 *
 * @todo use the conflict analysis to analyze the infeasibility which arise after the probing of the cover worked and
 *       solve returned infeasible, instead of adding the Nogood/Conflict by hand; that has the advantage that the SCIP
 *       takes care of creating the conflict and might shrink the initial reason
 *
 * @todo do not use LP and NLP fixing values in the same run, e.g., fixingalts = "lni", but start a second dive if LP
 *       values fail
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_undercover.h"
#include "scip/pub_misc.h"
#include "nlpi/exprinterpret.h"

#define HEUR_NAME               "undercover"
#define HEUR_DESC               "solves a sub-CIP determined by a set covering approach"
#define HEUR_DISPCHAR           'U'
#define HEUR_PRIORITY           -1110000
#define HEUR_FREQ               0
#define HEUR_FREQOFS            0
#define HEUR_MAXDEPTH           -1
#define HEUR_TIMING             SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP        TRUE         /**< does the heuristic use a secondary SCIP instance? */

/* default values for user parameters, grouped by parameter type */
#define DEFAULT_FIXINGALTS      "li"         /**< sequence of fixing values used: 'l'p relaxation, 'n'lp relaxation, 'i'ncumbent solution */

#define DEFAULT_MAXNODES        (SCIP_Longint)500/**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINNODES        (SCIP_Longint)500/**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS        (SCIP_Longint)500/**< number of nodes added to the contingent of the total nodes */

#define DEFAULT_CONFLICTWEIGHT  1000.0       /**< weight for conflict score in fixing order */
#define DEFAULT_CUTOFFWEIGHT    1.0          /**< weight for cutoff score in fixing order */
#define DEFAULT_INFERENCEWEIGHT 1.0          /**< weight for inference score in fixing order */
#define DEFAULT_MAXCOVERSIZEVARS  1.0           /**< maximum coversize (as fraction of total number of variables) */
#define DEFAULT_MAXCOVERSIZECONSS SCIP_REAL_MAX /**< maximum coversize (as ratio to the percentage of non-affected constraints) */
#define DEFAULT_MINCOVEREDREL   0.15         /**< minimum percentage of nonlinear constraints in the original problem */
#define DEFAULT_MINCOVEREDABS   5            /**< minimum number of nonlinear constraints in the original problem */
#define DEFAULT_MINIMPROVE      0.0          /**< factor by which heuristic should at least improve the incumbent */
#define DEFAULT_NODESQUOT       0.1          /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_RECOVERDIV      0.9          /**< fraction of covering variables in the last cover which need to change their value when recovering */

#define DEFAULT_MAXBACKTRACKS   6            /**< maximum number of backtracks */
#define DEFAULT_MAXRECOVERS     0            /**< maximum number of recoverings */
#define DEFAULT_MAXREORDERS     1            /**< maximum number of reorderings of the fixing order */

#define DEFAULT_COVERINGOBJ     'u'          /**< objective function of the covering problem */
#define DEFAULT_FIXINGORDER     'v'          /**< order in which variables should be fixed */

#define DEFAULT_BEFORECUTS      TRUE         /**< should undercover called at root node before cut separation? */
#define DEFAULT_FIXINTFIRST     FALSE        /**< should integer variables in the cover be fixed first? */
#define DEFAULT_LOCKSROUNDING   TRUE         /**< shall LP values for integer vars be rounded according to locks? */
#define DEFAULT_ONLYCONVEXIFY   FALSE        /**< should we only fix/dom.red. variables creating nonconvexity? */
#define DEFAULT_POSTNLP         TRUE         /**< should the NLP heuristic be called to polish a feasible solution? */
#define DEFAULT_COVERBD         FALSE        /**< should bounddisjunction constraints be covered (or just copied)? */
#define DEFAULT_REUSECOVER      FALSE        /**< shall the cover be re-used if a conflict was added after an infeasible subproblem? */
#define DEFAULT_COPYCUTS        TRUE         /**< should all active cuts from the cutpool of the original scip be copied
                                              *   to constraints of the subscip
                                              */
#define DEFAULT_RANDSEED        43           /* initial random seed */

/* local defines */
#define COVERINGOBJS            "cdlmtu"     /**< list of objective functions of the covering problem */
#define FIXINGORDERS            "CcVv"       /**< list of orders in which variables can be fixed */
#define MAXNLPFAILS             1            /**< maximum number of fails after which we give up solving the NLP relaxation */
#define MAXPOSTNLPFAILS         1            /**< maximum number of fails after which we give up calling NLP local search */
#define MINTIMELEFT             1.0          /**< don't start expensive parts of the heuristics if less than this amount of time left */
#define SUBMIPSETUPCOSTS        200          /**< number of nodes equivalent for the costs for setting up the sub-CIP */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_CONSHDLR**       nlconshdlrs;        /**< array of nonlinear constraint handlers */
   SCIP_HEUR*            nlpheur;            /**< pointer to NLP local search heuristics */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   char*                 fixingalts;         /**< sequence of fixing values used: 'l'p relaxation, 'n'lp relaxation, 'i'ncumbent solution */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Longint          nusednodes;         /**< nodes already used by heuristic in earlier calls */
   SCIP_Real             conflictweight;     /**< weight for conflict score in fixing order */
   SCIP_Real             cutoffweight;       /**< weight for cutoff score in fixing order */
   SCIP_Real             inferenceweight;    /**< weight for inference score in foxing order */
   SCIP_Real             maxcoversizevars;   /**< maximum coversize (as fraction of total number of variables) */
   SCIP_Real             maxcoversizeconss;  /**< maximum coversize maximum coversize (as ratio to the percentage of non-affected constraints) */
   SCIP_Real             mincoveredrel;      /**< minimum percentage of nonlinear constraints in the original problem */
   SCIP_Real             minimprove;         /**< factor by which heuristic should at least improve the incumbent */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             recoverdiv;         /**< fraction of covering variables in the last cover which need to change their value when recovering */
   int                   mincoveredabs;      /**< minimum number of nonlinear constraints in the original problem */
   int                   maxbacktracks;      /**< maximum number of backtracks */
   int                   maxrecovers;        /**< maximum number of recoverings */
   int                   maxreorders;        /**< maximum number of reorderings of the fixing order */
   int                   nfixingalts;        /**< number of fixing alternatives */
   int                   nnlpfails;          /**< number of fails when solving the NLP relaxation after last success */
   int                   npostnlpfails;      /**< number of fails of the NLP local search after last success */
   int                   nnlconshdlrs;       /**< number of nonlinear constraint handlers */
   char                  coveringobj;        /**< objective function of the covering problem */
   char                  fixingorder;        /**< order in which variables should be fixed */
   SCIP_Bool             beforecuts;         /**< should undercover be called at root node before cut separation? */
   SCIP_Bool             fixintfirst;        /**< should integer variables in the cover be fixed first? */
   SCIP_Bool             globalbounds;       /**< should global bounds on variables be used instead of local bounds at focus node? */
   SCIP_Bool             locksrounding;      /**< shall LP values for integer vars be rounded according to locks? */
   SCIP_Bool             nlpsolved;          /**< has current NLP relaxation already been solved successfully? */
   SCIP_Bool             nlpfailed;          /**< has solving the NLP relaxation failed? */
   SCIP_Bool             onlyconvexify;      /**< should we only fix/dom.red. variables creating nonconvexity? */
   SCIP_Bool             postnlp;            /**< should the NLP heuristic be called to polish a feasible solution? */
   SCIP_Bool             coverbd;            /**< should bounddisjunction constraints be covered (or just copied)? */
   SCIP_Bool             reusecover;         /**< shall the cover be re-used if a conflict was added after an infeasible subproblem? */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   subproblem? */
};

/** working memory for retrieving dense sparsity of Hessian matrices */
struct HessianData
{
   SCIP_SOL*             evalsol;
   SCIP_Real*            varvals;
   SCIP_Bool*            sparsity;
   int                   nvars;
};

/*
 * Local methods
 */


/** determines, whether a variable is fixed to the given value */
static
SCIP_Bool varIsFixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_Real             val,                /**< value to check */
   SCIP_Bool             global              /**< should global bounds be used? */
   )
{
   SCIP_Bool isfixed;

   if( global )
      isfixed = SCIPisFeasEQ(scip, val, SCIPvarGetLbGlobal(var)) && SCIPisFeasEQ(scip, val, SCIPvarGetUbGlobal(var));
   else
      isfixed = SCIPisFeasEQ(scip, val, SCIPvarGetLbLocal(var)) && SCIPisFeasEQ(scip, val, SCIPvarGetUbLocal(var));

   return isfixed;
}


/** determines, whether a term is already constant, because the variable is fixed or the coefficient is zero */
static
SCIP_Bool termIsConstant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   SCIP_Real             coeff,              /**< coefficient to check */
   SCIP_Bool             global              /**< should global bounds be used? */
   )
{
   /* if the variable has zero coefficient in the original problem, the term is linear */
   if( SCIPisZero(scip, coeff) )
      return TRUE;

   /* if the variable is fixed in the original problem, the term is linear */
   if( global )
      return SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   else
      return SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

}


/** determines, whether a term is convex */
static
SCIP_Bool termIsConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lhs,                /**< left hand side of the constraint */
   SCIP_Real             rhs,                /**< right hand side of the constraint */
   SCIP_Bool             sign                /**< signature of the term */
   )
{
   return sign ? SCIPisInfinity(scip, -lhs) : SCIPisInfinity(scip, rhs);
}


/** increases counters */
static
void  incCounters(
   int*                  termcounter,        /**< array to count in how many nonlinear terms a variable appears */
   int*                  conscounter,        /**< array to count in how many constraints a variable appears */
   SCIP_Bool*            consmarker,         /**< was this variable already counted for this constraint? */
   int                   idx                 /**< problem index of the variable */
   )
{
   termcounter[idx]++;
   if( !consmarker[idx] )
   {
      conscounter[idx]++;
      consmarker[idx] = TRUE;
   }
   return;
}


/** update time limit */
static
SCIP_RETCODE updateTimelimit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_Real*            timelimit           /**< time limit */
   )
{
   *timelimit -= SCIPgetClockTime(scip, clck);
   SCIP_CALL( SCIPresetClock(scip, clck) );
   SCIP_CALL( SCIPstartClock(scip, clck) );

   return SCIP_OKAY;
}


/** analyzes a nonlinear row and adds constraints and fixings to the covering problem */
static
SCIP_RETCODE processNlRow(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< nonlinear row representation of a nonlinear constraint */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter for computing sparsity pattern of the Hessian;
                                              *   if NULL, we will simply fix all variables in the expression tree */
   struct HessianData*   hessiandata,        /**< working memory for retrieving dense sparsity of Hessian matrices */
   SCIP*                 coveringscip,       /**< SCIP data structure for the covering problem */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            coveringvars,       /**< array to store the covering problem's variables */
   int*                  termcounter,        /**< counter array for number of nonlinear nonzeros per variable */
   int*                  conscounter,        /**< counter array for number of nonlinear constraints per variable */
   SCIP_Bool*            consmarker,         /**< marker array if constraint has been counted in conscounter */
   SCIP_Bool             globalbounds,       /**< should global bounds on variables be used instead of local bounds at focus node? */
   SCIP_Bool             onlyconvexify,      /**< should we only fix/dom.red. variables creating nonconvexity? */
   SCIP_Bool*            success             /**< pointer to store whether row was processed successfully */
   )
{
   SCIP_EXPRTREE* exprtree;
   SCIP_Bool infeas;
   SCIP_Bool fixed;
   int t;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(nlrow != NULL);
   assert(exprint == NULL || hessiandata != NULL);
   assert(coveringscip != NULL);
   assert(nvars >= 1);
   assert(coveringvars != NULL);
   assert(termcounter != NULL);
   assert(conscounter != NULL);
   assert(consmarker != NULL);
   assert(success != NULL);

   *success = FALSE;
   BMSclearMemoryArray(consmarker, nvars);

   /* go through expression tree */
   exprtree = SCIPnlrowGetExprtree(nlrow);
   if( exprtree != NULL )
   {
      SCIP_VAR** exprtreevars;
      int nexprtreevars;
      int probidx1;
      int probidx2;

      /* get variables in expression tree */
      nexprtreevars = SCIPexprtreeGetNVars(exprtree);
      exprtreevars = SCIPexprtreeGetVars(exprtree);

      if( exprtreevars != NULL && nexprtreevars > 0 )
      {
         SCIP_Bool usehessian;
         int i;

         /* is an expression interpreter available which can return the sparsity pattern of the Hessian? */
         usehessian = exprint != NULL && (SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_HESSIAN);
         if( usehessian )
         {
            int idx1;
            int idx2;

            assert(hessiandata != NULL);
            assert(hessiandata->nvars == 0 || hessiandata->varvals != NULL);
            assert(hessiandata->nvars == 0 || hessiandata->sparsity != NULL);

            /* compile expression tree */
            SCIP_CALL( SCIPexprintCompile(exprint, exprtree) );

            /* ensure memory */
            if( hessiandata->nvars < nexprtreevars )
            {
               SCIP_CALL( SCIPreallocBufferArray(scip, &hessiandata->varvals, nexprtreevars) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &hessiandata->sparsity, nexprtreevars*nexprtreevars) );
               hessiandata->nvars = nexprtreevars;
            }

            /* get point at which to evaluate the Hessian sparsity */
            if( hessiandata->evalsol == NULL )
            {
               SCIP_CALL( SCIPcreateSol(scip, &hessiandata->evalsol, NULL) );
               SCIP_CALL( SCIPlinkCurrentSol(scip, hessiandata->evalsol) );
            }
            SCIP_CALL( SCIPgetSolVals(scip, hessiandata->evalsol, nexprtreevars, exprtreevars, hessiandata->varvals) );

            /* get sparsity of the Hessian at current LP solution */
            SCIP_CALL( SCIPexprintHessianSparsityDense(exprint, exprtree, hessiandata->varvals, hessiandata->sparsity) );

            for( idx1 = nexprtreevars-1; idx1 >= 0; idx1-- )
            {
               /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
               probidx1 = SCIPvarGetProbindex(exprtreevars[idx1]);
               if( probidx1 == -1 )
               {
                  SCIPdebugMsg(scip, "strange: inactive variables detected in nonlinear row <%s>\n", SCIPnlrowGetName(nlrow));
                  return SCIP_OKAY;
               }

               /* nonzero diagonal element of the Hessian: fix */
               if( hessiandata->sparsity[idx1*nexprtreevars + idx1]
                  && !termIsConstant(scip, exprtreevars[idx1], 1.0, globalbounds) )
               {
                  SCIP_CALL( SCIPfixVar(coveringscip, coveringvars[probidx1], 1.0, &infeas, &fixed) );
                  assert(!infeas);
                  assert(fixed);

                  /* update counters */
                  incCounters(termcounter, conscounter, consmarker, probidx1);

                  SCIPdebugMsg(scip, "fixing var <%s> in covering problem to 1\n", SCIPvarGetName(coveringvars[probidx1]));

                  /* if covering variable is fixed, then no need to still check non-diagonal elements */
                  continue;
               }

               /* two different variables relate nonlinearly */
               for( idx2 = nexprtreevars-1; idx2 > idx1; idx2-- )
               {
                  SCIP_CONS* coveringcons;
                  SCIP_VAR* coveringconsvars[2];

                  /* do not assume symmetry */
                  if( !hessiandata->sparsity[idx1*nexprtreevars + idx2] && !hessiandata->sparsity[idx2*nexprtreevars + idx1] )
                     continue;

                  /* if diagonal has entry already, then covering constraint would always be satisfied, thus no need to add */
                  if( hessiandata->sparsity[idx2*nexprtreevars + idx2] && !termIsConstant(scip, exprtreevars[idx2], 1.0, globalbounds) )
                     continue;

                  /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
                  probidx2 = SCIPvarGetProbindex(exprtreevars[idx2]);
                  if( probidx2 == -1 )
                  {
                     SCIPdebugMsg(scip, "strange: inactive variables detected in nonlinear row <%s>\n", SCIPnlrowGetName(nlrow));
                     return SCIP_OKAY;
                  }

                  /* if the term is linear because one of the variables is fixed, nothing to do */
                  if( termIsConstant(scip, exprtreevars[idx1], 1.0, globalbounds)
                     || termIsConstant(scip, exprtreevars[idx2], 1.0, globalbounds) )
                     continue;

                  /* create covering constraint */
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering_%d_%d", SCIPnlrowGetName(nlrow), idx1, idx2);
                  coveringconsvars[0] = coveringvars[probidx1];
                  coveringconsvars[1] = coveringvars[probidx2];
                  SCIP_CALL( SCIPcreateConsSetcover(coveringscip, &coveringcons, name, 2, coveringconsvars,
                        TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

                  if( coveringcons == NULL )
                  {
                     SCIPdebugMsg(scip, "failed to create set covering constraint <%s>\n", name);
                     return SCIP_OKAY;
                  }

                  /* add and release covering constraint */
                  SCIP_CALL( SCIPaddCons(coveringscip, coveringcons) );
                  SCIP_CALL( SCIPreleaseCons(coveringscip, &coveringcons) );

                  SCIPdebugMsg(scip, "added covering constraint for vars <%s> and <%s> in covering problem\n", SCIPvarGetName(coveringvars[probidx1]), SCIPvarGetName(coveringvars[probidx2]));

                  /* update counters for both variables */
                  incCounters(termcounter, conscounter, consmarker, probidx1);
                  incCounters(termcounter, conscounter, consmarker, probidx2);
               }
            }
         }
         /* fix all variables contained in the expression tree */
         else
         {
            for( i = nexprtreevars-1; i >= 0; i-- )
            {
               assert(exprtreevars[i] != NULL);

               /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
               probidx1 = SCIPvarGetProbindex(exprtreevars[i]);
               if( probidx1 == -1 )
               {
                  SCIPdebugMsg(scip, "strange: inactive variable <%s> detected in nonlinear row <%s>\n",
                     SCIPvarGetName(exprtreevars[i]), SCIPnlrowGetName(nlrow));
                  return SCIP_OKAY;
               }

               /* term is constant, nothing to do */
               if( termIsConstant(scip, exprtreevars[i], 1.0, globalbounds) )
                  continue;

               /* otherwise fix variable */
               SCIP_CALL( SCIPfixVar(coveringscip, coveringvars[probidx1], 1.0, &infeas, &fixed) );
               assert(!infeas);
               assert(fixed);

               /* update counters */
               incCounters(termcounter, conscounter, consmarker, probidx1);

               SCIPdebugMsg(scip, "fixing var <%s> in covering problem to 1\n", SCIPvarGetName(coveringvars[probidx1]));
            }
         }
      }
   }

   /* go through all quadratic terms */
   for( t = SCIPnlrowGetNQuadElems(nlrow)-1; t >= 0; t-- )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_VAR* bilinvar1;
      SCIP_VAR* bilinvar2;
      int probidx1;
      int probidx2;

      /* get quadratic term */
      quadelem = &SCIPnlrowGetQuadElems(nlrow)[t];

      /* get involved variables */
      bilinvar1 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx1];
      bilinvar2 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx2];
      assert(bilinvar1 != NULL);
      assert(bilinvar2 != NULL);

      /* if constraints with inactive variables are present, we will have difficulties creating the sub-CIP later */
      probidx1 = SCIPvarGetProbindex(bilinvar1);
      probidx2 = SCIPvarGetProbindex(bilinvar2);
      if( probidx1 == -1 || probidx2 == -1 )
      {
         SCIPdebugMsg(scip, "inactive variables detected in nonlinear row <%s>\n", SCIPnlrowGetName(nlrow));
         return SCIP_OKAY;
      }

      /* we have a square term */
      if( bilinvar1 == bilinvar2 )
      {
         /* term is constant, nothing to do */
         if( termIsConstant(scip, bilinvar1, quadelem->coef, globalbounds) )
            continue;

         /* if we only convexify and term is convex considering the bounds of the nlrow, nothing to do */
         if( onlyconvexify && termIsConvex(scip, SCIPnlrowGetLhs(nlrow), SCIPnlrowGetRhs(nlrow), quadelem->coef >= 0) )
            continue;

         /* otherwise variable has to be in the cover */
         SCIP_CALL( SCIPfixVar(coveringscip, coveringvars[probidx1], 1.0, &infeas, &fixed) );
         assert(!infeas);
         assert(fixed);

         /* update counters */
         incCounters(termcounter, conscounter, consmarker, probidx1);

         SCIPdebugMsg(scip, "fixing var <%s> in covering problem to 1\n", SCIPvarGetName(coveringvars[probidx1]));
      }
      /* we have a bilinear term */
      else
      {
         SCIP_CONS* coveringcons;
         SCIP_VAR* coveringconsvars[2];

         /* if the term is linear because one of the variables is fixed or the coefficient is zero, nothing to do */
         if( termIsConstant(scip, bilinvar1, quadelem->coef, globalbounds)
            || termIsConstant(scip, bilinvar2, quadelem->coef, globalbounds) )
            continue;

         /* create covering constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering%d", SCIPnlrowGetName(nlrow), t);
         coveringconsvars[0] = coveringvars[probidx1];
         coveringconsvars[1] = coveringvars[probidx2];
         SCIP_CALL( SCIPcreateConsSetcover(coveringscip, &coveringcons, name, 2, coveringconsvars,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

         if( coveringcons == NULL )
         {
            SCIPdebugMsg(scip, "failed to create set covering constraint <%s>\n", name);
            return SCIP_OKAY;
         }

         /* add and release covering constraint */
         SCIP_CALL( SCIPaddCons(coveringscip, coveringcons) );
         SCIP_CALL( SCIPreleaseCons(coveringscip, &coveringcons) );

         /* update counters for both variables */
         incCounters(termcounter, conscounter, consmarker, probidx1);
         incCounters(termcounter, conscounter, consmarker, probidx2);
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}


/** creates the covering problem to determine a number of variables to be fixed */
static
SCIP_RETCODE createCoveringProblem(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 coveringscip,       /**< SCIP data structure for the covering problem */
   SCIP_VAR**            coveringvars,       /**< array to store the covering problem's variables */
   SCIP_Bool             globalbounds,       /**< should global bounds on variables be used instead of local bounds at focus node? */
   SCIP_Bool             onlyconvexify,      /**< should we only fix/dom.red. variables creating nonconvexity? */
   SCIP_Bool             coverbd,            /**< should bounddisjunction constraints be covered (or just copied)? */
   char                  coveringobj,        /**< objective function of the covering problem */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSHDLR* conshdlr;
   SCIP_HASHMAP* nlrowmap;
   SCIP_EXPRINT* exprint;
   SCIP_Bool* consmarker;
   int* conscounter;
   int* termcounter;

   struct HessianData hessiandata;
   int nlocksup;
   int nlocksdown;
   int nvars;
   int i;
   int probindex;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(coveringscip != NULL);
   assert(coveringvars != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create problem data structure */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPgetProbName(scip));
   SCIP_CALL( SCIPcreateProb(coveringscip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* allocate and initialize to zero counter arrays for weighted objectives */
   SCIP_CALL( SCIPallocBufferArray(scip, &consmarker, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conscounter, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termcounter, nvars) );
   BMSclearMemoryArray(conscounter, nvars);
   BMSclearMemoryArray(termcounter, nvars);

   /* create expression interpreter */
   SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &exprint) );
   assert(exprint != NULL);

   /* initialize empty hessiandata; memory will be allocated in method processNlRow() as required */
   hessiandata.evalsol = NULL;
   hessiandata.varvals = NULL;
   hessiandata.sparsity = NULL;
   hessiandata.nvars = 0;

   /* create covering variable for each variable in the original problem (fix it or not?) in the same order as in the
    * original problem
    */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real ub = 1.0;

      /* if the variable in the original problem is fixed, then the corresponding cover variable cannot be 1 in any
       * optimal solution of the covering problem (see special termIsConstant treatment below)
       * since some calling code may assume that no fixed variables will appear in the cover (see #1845), but we
       * might not compute an optimal cover here, we fix these variable to 0 here
       */
      if( globalbounds )
      {
         if( SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i])) )
            ub = 0.0;
      }
      else
      {
         if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) )
            ub = 0.0;
      }

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPvarGetName(vars[i]));
      SCIP_CALL( SCIPcreateVar(coveringscip, &coveringvars[i], name, 0.0, ub, 1.0, SCIP_VARTYPE_BINARY,
            TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      assert(coveringvars[i] != NULL);
      SCIP_CALL( SCIPaddVar(coveringscip, coveringvars[i]) );
   }

   /* first, go through some special constraint handlers which we do not want to treat by looking at their nlrow
    * representation; we store these in a hash map and afterwards process all nlrows which are not found in the hash map 
    */
   nlrowmap = NULL;
   if( SCIPisNLPConstructed(scip) )
   {
      int nnlprows;

      nnlprows = SCIPgetNNLPNlRows(scip);
      if( nnlprows > 0 )
      {
         int mapsize;

         /* calculate size of hash map */
         conshdlr = SCIPfindConshdlr(scip, "quadratic");
         mapsize = (conshdlr == NULL) ? 0 : SCIPconshdlrGetNActiveConss(conshdlr);
         conshdlr = SCIPfindConshdlr(scip, "soc");
         if( conshdlr != NULL )
            mapsize += SCIPconshdlrGetNActiveConss(conshdlr);
         mapsize = MAX(mapsize, nnlprows);
         assert(mapsize > 0);

         /* create hash map */
         SCIP_CALL( SCIPhashmapCreate(&nlrowmap, SCIPblkmem(scip), mapsize) );
         assert(nlrowmap != NULL);
      }
   }

   /* go through all AND constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "and");
   if( conshdlr != NULL )
   {
      int c;

      for( c = SCIPconshdlrGetNActiveConss(conshdlr)-1; c >= 0; c-- )
      {
         SCIP_CONS* andcons;
         SCIP_CONS* coveringcons;
         SCIP_VAR** andvars;
         SCIP_VAR* andres;
         SCIP_VAR** coveringconsvars;
         SCIP_Real* coveringconsvals;
         SCIP_Bool negated;

         int ntofix;
         int v;

         /* get constraint and variables */
         andcons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(andcons != NULL);
         andvars = SCIPgetVarsAnd(scip, andcons);
         assert(andvars != NULL);

         /* "and" constraints are not passed to the NLP, hence nothing to store in the hash map */

         /* allocate memory for covering constraint */
         SCIP_CALL( SCIPallocBufferArray(coveringscip, &coveringconsvars, SCIPgetNVarsAnd(scip, andcons)+1) );
         SCIP_CALL( SCIPallocBufferArray(coveringscip, &coveringconsvals, SCIPgetNVarsAnd(scip, andcons)+1) );

         /* collect unfixed variables */
         BMSclearMemoryArray(consmarker, nvars);
         ntofix = 0;
         for( v = SCIPgetNVarsAnd(scip, andcons)-1; v >= 0; v-- )
         {
            assert(andvars[v] != NULL);
            negated = FALSE;

            /* if variable is fixed to 0, entire constraint can be linearized */
            if( varIsFixed(scip, andvars[v], 0.0, globalbounds) )
            {
               ntofix = 0;
               break;
            }

            /* if variable is fixed, nothing to do */
            if( termIsConstant(scip, andvars[v], 1.0, globalbounds) )
            {
               continue;
            }

            /* if constraints with inactive variables are present, we have to find the corresponding active variable */
            probindex = SCIPvarGetProbindex(andvars[v]);
            if( probindex == -1 )
            {
               SCIP_VAR* repvar;

               /* get binary representative of variable */
               SCIP_CALL( SCIPgetBinvarRepresentative(scip, andvars[v], &repvar, &negated) );
               assert(repvar != NULL);
               assert(SCIPvarGetStatus(repvar) != SCIP_VARSTATUS_FIXED);

               if( SCIPvarGetStatus(repvar) == SCIP_VARSTATUS_MULTAGGR )
               {
                  SCIPdebugMsg(scip, "strange: multiaggregated variable found <%s>\n", SCIPvarGetName(andvars[v]));
                  SCIPdebugMsg(scip, "inactive variables detected in constraint <%s>\n", SCIPconsGetName(andcons));
                  SCIPfreeBufferArray(coveringscip, &coveringconsvals);
                  SCIPfreeBufferArray(coveringscip, &coveringconsvars);
                  goto TERMINATE;
               }

               /* check for negation */
               if( SCIPvarIsNegated(repvar) )
               {
                  probindex = SCIPvarGetProbindex(SCIPvarGetNegationVar(repvar));
                  negated = TRUE;
               }
               else
               {
                  assert(SCIPvarIsActive(repvar));
                  probindex = SCIPvarGetProbindex(repvar);
                  negated = FALSE;
               }
            }
            assert(probindex >= 0);

            /* add covering variable for unfixed original variable */
            if( negated )
            {
               SCIP_CALL( SCIPgetNegatedVar(coveringscip, coveringvars[probindex], &coveringconsvars[ntofix]) );
            }
            else
               coveringconsvars[ntofix] = coveringvars[probindex];
            coveringconsvals[ntofix] = 1.0;
            ntofix++;
         }
         negated = FALSE;

         /* if constraints with inactive variables are present, we have to find the corresponding active variable */
         andres = SCIPgetResultantAnd(scip, andcons);
         assert(andres != NULL);
         probindex = SCIPvarGetProbindex(andres);

         /* if resultant is fixed this constraint can be either linearized or is redundant because all operands can be fixed */
         if( termIsConstant(scip, andres, 1.0, globalbounds) )
         {
            /* free memory for covering constraint */
            SCIPfreeBufferArray(coveringscip, &coveringconsvals);
            SCIPfreeBufferArray(coveringscip, &coveringconsvars);

            continue;
         }

         if( probindex == -1 )
         {
            SCIP_VAR* repvar;

            /* get binary representative of variable */
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, SCIPgetResultantAnd(scip, andcons), &repvar, &negated) );
            assert(repvar != NULL);
            assert(SCIPvarGetStatus(repvar) != SCIP_VARSTATUS_FIXED);

            if( SCIPvarGetStatus(repvar) == SCIP_VARSTATUS_MULTAGGR )
            {
               SCIPdebugMsg(scip, "strange: multiaggregated variable found <%s>\n", SCIPvarGetName(SCIPgetResultantAnd(scip, andcons)));
               SCIPdebugMsg(scip, "inactive variables detected in constraint <%s>\n", SCIPconsGetName(andcons));
               SCIPfreeBufferArray(coveringscip, &coveringconsvals);
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* check for negation */
            if( SCIPvarIsNegated(repvar) )
            {
               probindex = SCIPvarGetProbindex(SCIPvarGetNegationVar(repvar));
               negated = TRUE;
            }
            else
            {
               assert(SCIPvarIsActive(repvar));
               probindex = SCIPvarGetProbindex(repvar);
               negated = FALSE;
            }
         }
         else if( SCIPvarGetLbGlobal(andres) > 0.5 || SCIPvarGetUbGlobal(andres) < 0.5 )
         {
            /* free memory for covering constraint */
            SCIPfreeBufferArray(coveringscip, &coveringconsvals);
            SCIPfreeBufferArray(coveringscip, &coveringconsvars);

            continue;
         }
         assert(probindex >= 0);
         assert(!termIsConstant(scip, (negated ? SCIPvarGetNegatedVar(vars[probindex]) : vars[probindex]), 1.0, globalbounds));

         /* if less than two variables are unfixed or the resultant variable is fixed, the entire constraint can be
	  * linearized anyway
	  */
         if( ntofix >= 2 )
         {
            assert(ntofix <= SCIPgetNVarsAnd(scip, andcons));

            /* add covering variable for unfixed resultant */
            if( negated )
            {
               SCIP_CALL( SCIPgetNegatedVar(coveringscip, coveringvars[probindex], &coveringconsvars[ntofix]) );
            }
            else
               coveringconsvars[ntofix] = coveringvars[probindex];
            coveringconsvals[ntofix] = (SCIP_Real)(ntofix - 1);
            ntofix++;

            /* create covering constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPconsGetName(andcons));
            SCIP_CALL( SCIPcreateConsLinear(coveringscip, &coveringcons, name, ntofix, coveringconsvars, coveringconsvals,
                  (SCIP_Real)(ntofix - 2), SCIPinfinity(coveringscip),
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( coveringcons == NULL )
            {
               SCIPdebugMsg(scip, "failed to create linear constraint <%s>\n", name);
               SCIPfreeBufferArray(coveringscip, &coveringconsvals);
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* add and release covering constraint */
            SCIP_CALL( SCIPaddCons(coveringscip, coveringcons) );
            SCIP_CALL( SCIPreleaseCons(coveringscip, &coveringcons) );

            /* update counters */
            for( v = ntofix-1; v >= 0; v-- )
               if( SCIPvarIsNegated(coveringconsvars[v]) )
                  incCounters(termcounter, conscounter, consmarker, SCIPvarGetProbindex(SCIPvarGetNegationVar(coveringconsvars[v])));
               else
                  incCounters(termcounter, conscounter, consmarker, SCIPvarGetProbindex(coveringconsvars[v]));
         }

         /* free memory for covering constraint */
         SCIPfreeBufferArray(coveringscip, &coveringconsvals);
         SCIPfreeBufferArray(coveringscip, &coveringconsvars);
      }
   }

   /* go through all bounddisjunction constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "bounddisjunction");
   if( conshdlr != NULL && coverbd )
   {
      int c;

      for( c = SCIPconshdlrGetNActiveConss(conshdlr)-1; c >= 0; c-- )
      {
         SCIP_CONS* bdcons;
         SCIP_CONS* coveringcons;
         SCIP_VAR** bdvars;
         SCIP_VAR** coveringconsvars;
         SCIP_Real* coveringconsvals;

         int nbdvars;
         int ntofix;
         int v;

         /* get constraint and variables */
         bdcons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(bdcons != NULL);
         bdvars = SCIPgetVarsBounddisjunction(scip, bdcons);
         assert(bdvars != NULL);
         nbdvars = SCIPgetNVarsBounddisjunction(scip, bdcons);

         /* bounddisjunction constraints are not passed to the NLP, hence nothing to store in the hash map */

         /* allocate memory for covering constraint */
         SCIP_CALL( SCIPallocBufferArray(coveringscip, &coveringconsvars, nbdvars) );
         SCIP_CALL( SCIPallocBufferArray(coveringscip, &coveringconsvals, nbdvars) );

         /* collect unfixed variables */
         BMSclearMemoryArray(consmarker, nvars);
         ntofix = 0;
         for( v = nbdvars-1; v >= 0; v-- )
         {
            SCIP_Bool negated;

            assert(bdvars[v] != NULL);
            negated = FALSE;

            /* if variable is fixed, nothing to do */
            if( varIsFixed(scip, bdvars[v], globalbounds ? SCIPvarGetLbGlobal(bdvars[v]) : SCIPvarGetLbLocal(bdvars[v]),
                  globalbounds) )
            {
               continue;
            }

            /* if constraints with inactive variables are present, we have to find the corresponding active variable */
            probindex = SCIPvarGetProbindex(bdvars[v]);
            if( probindex == -1 )
            {
               SCIP_VAR* repvar;

               /* get binary representative of variable */
               SCIP_CALL( SCIPgetBinvarRepresentative(scip, bdvars[v], &repvar, &negated) );
               assert(repvar != NULL);
               assert(SCIPvarGetStatus(repvar) != SCIP_VARSTATUS_FIXED);

               if( SCIPvarGetStatus(repvar) == SCIP_VARSTATUS_MULTAGGR )
               {
                  SCIPdebugMsg(scip, "strange: multiaggregated variable found <%s>\n", SCIPvarGetName(bdvars[v]));
                  SCIPdebugMsg(scip, "inactive variables detected in constraint <%s>\n", SCIPconsGetName(bdcons));
                  SCIPfreeBufferArray(coveringscip, &coveringconsvals);
                  SCIPfreeBufferArray(coveringscip, &coveringconsvars);
                  goto TERMINATE;
               }

               /* check for negation */
               if( SCIPvarIsNegated(repvar) )
               {
                  probindex = SCIPvarGetProbindex(SCIPvarGetNegationVar(repvar));
                  negated = TRUE;
               }
               else
               {
                  assert(SCIPvarIsActive(repvar));
                  probindex = SCIPvarGetProbindex(repvar);
                  negated = FALSE;
               }
            }
            assert(probindex >= 0);

            /* add covering variable for unfixed original variable */
            if( negated )
            {
               SCIP_CALL( SCIPgetNegatedVar(coveringscip, coveringvars[probindex], &coveringconsvars[ntofix]) );
            }
            else
               coveringconsvars[ntofix] = coveringvars[probindex];
            coveringconsvals[ntofix] = 1.0;
            ntofix++;
         }

         /* if less than two variables are unfixed, the entire constraint can be linearized anyway */
         if( ntofix >= 2 )
         {
            assert(ntofix <= nbdvars);

            /* create covering constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPconsGetName(bdcons));
            SCIP_CALL( SCIPcreateConsLinear(coveringscip, &coveringcons, name, ntofix, coveringconsvars, coveringconsvals,
                  (SCIP_Real)(ntofix - 1), SCIPinfinity(coveringscip),
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( coveringcons == NULL )
            {
               SCIPdebugMsg(scip, "failed to create linear constraint <%s>\n", name);
               SCIPfreeBufferArray(coveringscip, &coveringconsvals);
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* add and release covering constraint */
            SCIP_CALL( SCIPaddCons(coveringscip, coveringcons) );
            SCIP_CALL( SCIPreleaseCons(coveringscip, &coveringcons) );

            /* update counters */
            for( v = ntofix-1; v >= 0; v-- )
               if( SCIPvarIsNegated(coveringconsvars[v]) )
                  incCounters(termcounter, conscounter, consmarker, SCIPvarGetProbindex(SCIPvarGetNegationVar(coveringconsvars[v])));
               else
                  incCounters(termcounter, conscounter, consmarker, SCIPvarGetProbindex(coveringconsvars[v]));
         }

         /* free memory for covering constraint */
         SCIPfreeBufferArray(coveringscip, &coveringconsvals);
         SCIPfreeBufferArray(coveringscip, &coveringconsvars);
      }
   }

   /* go through all indicator constraints in the original problem; fix the binary variable */
   conshdlr = SCIPfindConshdlr(scip, "indicator");
   if( conshdlr != NULL )
   {
      int c;

      for( c = SCIPconshdlrGetNActiveConss(conshdlr)-1; c >= 0; c-- )
      {
         SCIP_CONS* indcons;
         SCIP_VAR* binvar;
         SCIP_VAR* coveringvar;

         /* get constraint and variables */
         indcons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(indcons != NULL);
         binvar = SCIPgetBinaryVarIndicator(indcons);
         assert(binvar != NULL);

         /* indicator constraints are not passed to the NLP, hence nothing to store in the hash map */

         /* if variable is fixed, nothing to do */
         if( varIsFixed(scip, binvar, globalbounds ? SCIPvarGetLbGlobal(binvar) : SCIPvarGetLbLocal(binvar), globalbounds) )
         {
            continue;
         }

         /* if constraints with inactive variables are present, we have to find the corresponding active variable */
         probindex = SCIPvarGetProbindex(binvar);
         if( probindex == -1 )
         {
            SCIP_VAR* repvar;
            SCIP_Bool negated;

            /* get binary representative of variable */
            negated = FALSE;
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, binvar, &repvar, &negated) );
            assert(repvar != NULL);
            assert(SCIPvarGetStatus(repvar) != SCIP_VARSTATUS_FIXED);

            if( SCIPvarGetStatus(repvar) == SCIP_VARSTATUS_MULTAGGR )
            {
               SCIPdebugMsg(scip, "strange: multiaggregated variable found <%s>\n", SCIPvarGetName(binvar));
               SCIPdebugMsg(scip, "inactive variables detected in constraint <%s>\n", SCIPconsGetName(indcons));
               goto TERMINATE;
            }

            /* check for negation */
            if( SCIPvarIsNegated(repvar) )
               probindex = SCIPvarGetProbindex(SCIPvarGetNegationVar(repvar));
            else
            {
               assert(SCIPvarIsActive(repvar));
               probindex = SCIPvarGetProbindex(repvar);
            }
         }
         assert(probindex >= 0);

         /* get covering variable for unfixed binary variable in indicator constraint */
         coveringvar = coveringvars[probindex];

         /* require covering variable to be fixed such that indicator is linearized */
         SCIP_CALL( SCIPchgVarLb(coveringscip, coveringvar, 1.0) );

         /* update counters */
         BMSclearMemoryArray(consmarker, nvars);
         incCounters(termcounter, conscounter, consmarker, probindex);
      }
   }

   /* go through all quadratic constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "quadratic");
   if( conshdlr != NULL )
   {
      int c;

      for( c = SCIPconshdlrGetNActiveConss(conshdlr)-1; c >= 0; c-- )
      {
         SCIP_CONS* quadcons;
         SCIP_NLROW* nlrow;

         /* get constraint */
         quadcons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(quadcons != NULL);

         /* get nlrow representation and store it in hash map */
         SCIP_CALL( SCIPgetNlRowQuadratic(scip, quadcons, &nlrow) );
         assert(nlrow != NULL);
         if( nlrowmap != NULL )
         {
            assert(!SCIPhashmapExists(nlrowmap, nlrow));
            SCIP_CALL( SCIPhashmapInsert(nlrowmap, nlrow, quadcons) );
         }

         /* if we only want to convexify and curvature and bounds prove already convexity, nothing to do */
         if( onlyconvexify
            && ((SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, quadcons)) && SCIPisConvexQuadratic(scip, quadcons))
               || (SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, quadcons)) && SCIPisConcaveQuadratic(scip, quadcons))) )
            continue;

         /* process nlrow */
         *success = FALSE;
         SCIP_CALL( processNlRow(scip, nlrow, exprint, &hessiandata, coveringscip, nvars, coveringvars,
               termcounter, conscounter, consmarker, globalbounds, onlyconvexify, success) );

         if( *success == FALSE )
            goto TERMINATE;
      }

      *success = FALSE;
   }

   /* go through all "soc" constraints in the original problem */
   conshdlr = SCIPfindConshdlr(scip, "soc");
   if( conshdlr != NULL && !onlyconvexify )
   {
      int c;

      for( c = SCIPconshdlrGetNActiveConss(conshdlr)-1; c >= 0; c-- )
      {
         SCIP_CONS* soccons;
         SCIP_CONS* coveringcons;
         SCIP_VAR** soclhsvars;
         SCIP_VAR* socrhsvar;
         SCIP_VAR** coveringconsvars;
         SCIP_NLROW* nlrow;

         int ntofix;
         int v;

         /* get constraints and variables */
         soccons = SCIPconshdlrGetConss(conshdlr)[c];
         assert(soccons != NULL);
         socrhsvar = SCIPgetRhsVarSOC(scip, soccons);
         assert(socrhsvar != NULL);
         soclhsvars = SCIPgetLhsVarsSOC(scip, soccons);
         assert(soclhsvars != NULL);

         /* get nlrow representation and store it in hash map */
         SCIP_CALL( SCIPgetNlRowSOC(scip, soccons, &nlrow) );
         assert(nlrow != NULL);
         if( nlrowmap != NULL )
         {
            assert(!SCIPhashmapExists(nlrowmap, nlrow));
            SCIP_CALL( SCIPhashmapInsert(nlrowmap, nlrow, soccons) );
         }

         /* allocate memory for covering constraint */
         SCIP_CALL( SCIPallocBufferArray(coveringscip, &coveringconsvars, SCIPgetNLhsVarsSOC(scip, soccons)+1) );

         /* collect unfixed variables */
         BMSclearMemoryArray(consmarker, nvars);
         ntofix = 0;

         /* soc constraints should contain only active and multi-aggregated variables; the latter we do not handle */
         probindex = SCIPvarGetProbindex(socrhsvar);
         if( probindex == -1 )
         {
            SCIPdebugMsg(scip, "inactive variables detected in constraint <%s>\n", SCIPconsGetName(soccons));
            SCIPfreeBufferArray(coveringscip, &coveringconsvars);
            goto TERMINATE;
         }

         /* add covering variable for unfixed rhs variable */
         if( !termIsConstant(scip, socrhsvar, SCIPgetRhsCoefSOC(scip, soccons), globalbounds) )
         {
            SCIP_CALL( SCIPgetNegatedVar(coveringscip, coveringvars[probindex], &coveringconsvars[ntofix]) );
            ntofix++;
         }

         /* go through lhs variables */
         for( v = SCIPgetNLhsVarsSOC(scip, soccons)-1; v >= 0; v-- )
         {
            assert(soclhsvars[v] != NULL);

            /* soc constraints should contain only active and multi-aggregated variables; the latter we do not handle */
            probindex = SCIPvarGetProbindex(soclhsvars[v]);
            if( probindex == -1 )
            {
               SCIPdebugMsg(scip, "inactive variables detected in constraint <%s>\n", SCIPconsGetName(soccons));
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* add covering variable for unfixed lhs variable */
            if( !termIsConstant(scip, soclhsvars[v], SCIPgetLhsCoefsSOC(scip, soccons)[v], globalbounds) )
            {
               SCIP_CALL( SCIPgetNegatedVar(coveringscip, coveringvars[probindex], &coveringconsvars[ntofix]) );
               ntofix++;
            }
         }

         if( ntofix >= 2 )
         {
            /* create covering constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_covering", SCIPconsGetName(soccons));
            SCIP_CALL( SCIPcreateConsSetpack(coveringscip, &coveringcons, name, ntofix, coveringconsvars,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

            if( coveringcons == NULL )
            {
               SCIPdebugMsg(scip, "failed to create set packing constraint <%s>\n", name);
               SCIPfreeBufferArray(coveringscip, &coveringconsvars);
               goto TERMINATE;
            }

            /* add and release covering constraint */
            SCIP_CALL( SCIPaddCons(coveringscip, coveringcons) );
            SCIP_CALL( SCIPreleaseCons(coveringscip, &coveringcons) );

            /* update counters */
            for( v = ntofix-1; v >= 0; v-- )
               incCounters(termcounter, conscounter, consmarker, SCIPvarGetProbindex(SCIPvarGetNegatedVar(coveringconsvars[v])));
         }

         /* free memory for covering constraint */
         SCIPfreeBufferArray(coveringscip, &coveringconsvars);
      }
   }

   /* go through all yet unprocessed nlrows */
   if( nlrowmap != NULL )
   {
      SCIP_NLROW** nlrows;
      int nnlrows;

      assert(SCIPisNLPConstructed(scip));

      /* get nlrows */
      nnlrows = SCIPgetNNLPNlRows(scip);
      nlrows = SCIPgetNLPNlRows(scip);

      for( i = nnlrows-1; i >= 0; i-- )
      {
         assert(nlrows[i] != NULL);

         /* nlrow or corresponding constraint already processed */
         if( SCIPhashmapExists(nlrowmap, nlrows[i]) )
            continue;

         /* process nlrow */
         *success = FALSE;
         SCIP_CALL( processNlRow(scip, nlrows[i], exprint, &hessiandata, coveringscip, nvars, coveringvars,
               termcounter, conscounter, consmarker, globalbounds, onlyconvexify, success) );

         if( *success == FALSE )
            goto TERMINATE;
      }
   }

   /* set objective function of covering problem */
   switch( coveringobj )
   {
   case 'c': /* number of influenced nonlinear constraints */
      for( i = nvars-1; i >= 0; i-- )
      {
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], (SCIP_Real) conscounter[i]) );
      }
      break;
   case 'd': /* domain size */
      for( i = nvars-1; i >= 0; i-- )
      {
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i],
               (globalbounds ? SCIPvarGetUbGlobal(vars[i]) - SCIPvarGetLbGlobal(vars[i]) : SCIPvarGetUbLocal(vars[i]) - SCIPvarGetLbLocal(vars[i]))) );
      }
      break;
   case 'l': /* number of locks */
      for( i = nvars-1; i >= 0; i-- )
      {
         nlocksup = SCIPvarGetNLocksUp(vars[i]);
         nlocksdown = SCIPvarGetNLocksDown(vars[i]);
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], (SCIP_Real) (nlocksup+nlocksdown+1)) );
      }
      break;
   case 'm': /* min(up locks, down locks)+1 */
      for( i = nvars-1; i >= 0; i-- )
      {
         nlocksup = SCIPvarGetNLocksUp(vars[i]);
         nlocksdown = SCIPvarGetNLocksDown(vars[i]);
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], (SCIP_Real) (MIN(nlocksup, nlocksdown)+1)) );
      }
      break;
   case 't': /* number of influenced nonlinear terms */
      for( i = nvars-1; i >= 0; i-- )
      {
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], (SCIP_Real) termcounter[i]) );
      }
      break;
   case 'u': /* unit penalties */
      for( i = nvars-1; i >= 0; i-- )
      {
         SCIP_CALL( SCIPchgVarObj(coveringscip, coveringvars[i], 1.0) );
      }
      break;
   default:
      SCIPerrorMessage("invalid choice <%c> for covering objective\n", coveringobj);
      goto TERMINATE;
   }

   /* covering problem successfully set up */
   *success = TRUE;

 TERMINATE:
   /* free nlrow hash map */
   if( nlrowmap != NULL )
   {
      SCIPhashmapFree(&nlrowmap);
   }

   /* free hessian data */
   if( hessiandata.nvars > 0 )
   {
      assert(hessiandata.evalsol != NULL);
      assert(hessiandata.varvals != NULL);
      assert(hessiandata.sparsity != NULL);

      SCIPfreeBufferArray(scip, &hessiandata.sparsity);
      SCIPfreeBufferArray(scip, &hessiandata.varvals);

      SCIP_CALL( SCIPfreeSol(scip, &hessiandata.evalsol) );
   }
   else
   {
      assert(hessiandata.evalsol == NULL);
      assert(hessiandata.varvals == NULL);
      assert(hessiandata.sparsity == NULL);
   }

   /* free expression interpreter */
   SCIP_CALL( SCIPexprintFree(&exprint) );

   SCIPstatistic(
      {
	 int nnonzs;
	 nnonzs = 0;
	 for( i = 0; i < nvars; ++i)
	    nnonzs += termcounter[i];
	 SCIPstatisticPrintf("UCstats nnz/var: %9.6f\n", nnonzs/(SCIP_Real)nvars);
	 nnonzs = 0;
	 for( i = 0; i < nvars; ++i)
	    if( conscounter[i] > 0 )
	       nnonzs++;
	 SCIPstatisticPrintf("UCstats nlvars: %6d\n", nnonzs);
      }
      );

   /* free counter arrays for weighted objectives */
   SCIPfreeBufferArray(scip, &termcounter);
   SCIPfreeBufferArray(scip, &conscounter);
   SCIPfreeBufferArray(scip, &consmarker);

   return SCIP_OKAY;
}


/** adds a constraint to the covering problem to forbid the given cover */
static
SCIP_RETCODE forbidCover(
   SCIP*                 scip,               /**< SCIP data structure of the covering problem */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< variable array */
   int                   coversize,          /**< size of the cover */
   int*                  cover,              /**< problem indices of the variables in the cover */
   int                   diversification,    /**< how many unfixed variables have to change their value? */
   SCIP_Bool*            success,            /**< pointer to store whether the cutoff constraint was created successfully */
   SCIP_Bool*            infeas              /**< pointer to store whether the cutoff proves (local or global) infeasibility */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR** consvars;

   char name[SCIP_MAXSTRLEN];
   int nconsvars;
   int i;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars >= 1);
   assert(cover != NULL);
   assert(coversize >= 1);
   assert(coversize <= nvars);
   assert(diversification >= 1);
   assert(success != NULL);
   assert(infeas != NULL);

   *success = FALSE;
   *infeas = FALSE;

   /* allocate memory for constraint variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, coversize) );
   nconsvars = 0;
   cons = NULL;

   /* create constraint name */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "forbid_cover_assignment");

   /* if all variables in the cover are binary and we require only one variable to change its value, then we create a
    * set covering constraint
    */
   if( diversification == 1 )
   {
      /* build up constraint */
      for( i = coversize-1; i >= 0; i-- )
      {
         if( !SCIPisFeasGE(scip, SCIPvarGetLbLocal(vars[cover[i]]), 1.0) )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[cover[i]], &consvars[nconsvars]) );
            nconsvars++;
         }
      }

      /* if all covering variables are fixed to one, the constraint cannot be satisfied */
      if( nconsvars == 0 )
      {
         *infeas = TRUE;
      }
      else
      {
         /* create constraint */
         SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, name, nconsvars, consvars,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
      }
   }
   /* if all variables in the cover are binary and we require more variables to change their value, then we create a
    * linear constraint
    */
   else
   {
      SCIP_Real* consvals;
      SCIP_Real rhs;

      /* build up constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, coversize) );
      for( i = coversize-1; i >= 0; i-- )
      {
         if( !SCIPisFeasGE(scip, SCIPvarGetLbLocal(vars[cover[i]]), 1.0) )
         {
            consvars[nconsvars] = vars[cover[i]];
            consvals[nconsvars] = 1.0;
            nconsvars++;
         }
      }
      rhs = (SCIP_Real) (nconsvars-diversification);

      /* if too many covering variables are fixed to 1, the constraint cannot be satisfied */
      if( rhs < 0 )
      {
         *infeas = TRUE;
      }
      else
      {
         /* create constraint */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name,
               nconsvars, consvars, consvals, -SCIPinfinity(scip), rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &consvals);
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &consvars);

   /* if proven infeasible, we do not even add the constraint; otherwise we add and release the constraint if created
    * successfully 
    */
   if( !(*infeas) && cons != NULL )
   {
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      *success = TRUE;
   }

   return SCIP_OKAY;
}


/** adds a set covering or bound disjunction constraint to the original problem */
static
SCIP_RETCODE createConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   bdlen,              /**< length of bound disjunction */
   SCIP_VAR**            bdvars,             /**< array of variables in bound disjunction */
   SCIP_BOUNDTYPE*       bdtypes,            /**< array of bound types in bound disjunction */
   SCIP_Real*            bdbounds,           /**< array of bounds in bound disjunction */
   SCIP_Bool             local,              /**< is constraint valid only locally? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool*            success             /**< pointer to store whether the cutoff constraint was created successfully */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   SCIP_Bool isbinary;
   char name[SCIP_MAXSTRLEN];
   int i;

   assert(scip != NULL);
   assert(bdlen >= 1);
   assert(bdvars != NULL);
   assert(bdtypes != NULL);
   assert(bdbounds != NULL);
   assert(success != NULL);

   /* initialize */
   *success = FALSE;
   cons = NULL;
   consvars = NULL;

   /* create constraint name */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "undercover_cutoff");

   /* check if all variables in the cover are binary */
   isbinary = TRUE;
   for( i = bdlen-1; i >= 0 && isbinary; i-- )
   {
      isbinary = SCIPvarIsBinary(bdvars[i]);
   }

   /* if all variables in the cover are binary, then we create a logicor constraint */
   if( isbinary )
   {
      /* allocate memory for constraint variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, bdlen) );

      /* build up constraint */
      for( i = bdlen-1; i >= 0; i-- )
      {
         assert(bdtypes[i] == SCIP_BOUNDTYPE_LOWER || SCIPisFeasZero(scip, bdbounds[i]));
         assert(bdtypes[i] == SCIP_BOUNDTYPE_UPPER || SCIPisFeasEQ(scip, bdbounds[i], 1.0));

         if( bdtypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            consvars[i] = bdvars[i];
         }
         else
         {
            assert(bdtypes[i] == SCIP_BOUNDTYPE_UPPER);
            SCIP_CALL( SCIPgetNegatedVar(scip, bdvars[i], &consvars[i]) );
         }
      }

      /* create conflict constraint */
      SCIP_CALL( SCIPcreateConsLogicor(scip, &cons, name, bdlen, consvars,
            FALSE, TRUE, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );
   }
   /* otherwise we create a bound disjunction constraint as given */
   else
   {
      /* create conflict constraint */
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, bdlen, bdvars, bdtypes, bdbounds,
            FALSE, TRUE, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );
   }

   /* add and release constraint if created successfully */
   if( cons != NULL )
   {
      if( local )
      {
         SCIP_CALL( SCIPaddConsLocal(scip, cons, NULL) );
      }
      else
      {
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      *success = TRUE;
   }

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &consvars);

   return SCIP_OKAY;
}


/** solve covering problem */
static
SCIP_RETCODE solveCoveringProblem(
   SCIP*                 coveringscip,       /**< SCIP data structure for the covering problem */
   int                   ncoveringvars,      /**< number of the covering problem's variables */
   SCIP_VAR**            coveringvars,       /**< array of the covering problem's variables */
   int*                  coversize,          /**< size of the computed cover */
   int*                  cover,              /**< array to store indices of the variables in the computed cover
                                              *   (should be ready to hold ncoveringvars entries) */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Real             objlimit,           /**< upper bound on the cover size */
   SCIP_Bool*            success             /**< feasible cover found? */
   )
{
   SCIP_Real* solvals;
   SCIP_Real totalpenalty;
   SCIP_RETCODE retcode;
   int i;

   assert(coveringscip != NULL);
   assert(coveringvars != NULL);
   assert(cover != NULL);
   assert(coversize != NULL);
   assert(timelimit > 0.0);
   assert(memorylimit > 0.0);
   assert(success != NULL);

   *success = FALSE;

   /* forbid call of heuristics and separators solving sub-CIPs */
   SCIP_CALL( SCIPsetSubscipsOff(coveringscip, TRUE) );

   /* set presolving and separation to fast */
   SCIP_CALL( SCIPsetSeparating(coveringscip, SCIP_PARAMSETTING_FAST, TRUE) );
   SCIP_CALL( SCIPsetPresolving(coveringscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use inference branching */
   if( SCIPfindBranchrule(coveringscip, "inference") != NULL && !SCIPisParamFixed(coveringscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(coveringscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* only solve root */
   SCIP_CALL( SCIPsetLongintParam(coveringscip, "limits/nodes", 1LL) );

   SCIPdebugMsg(coveringscip, "timelimit = %g, memlimit = %g\n", timelimit, memorylimit);

   /* set time, memory, and objective limit */
   SCIP_CALL( SCIPsetRealParam(coveringscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(coveringscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetObjlimit(coveringscip, objlimit) );

   /* do not abort on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(coveringscip, "misc/catchctrlc", FALSE) );

   /* disable output to console in optimized mode, enable in SCIP's debug mode */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(coveringscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(coveringscip, "display/freq", 100000) );
#else
   SCIP_CALL( SCIPsetIntParam(coveringscip, "display/verblevel", 0) );
#endif


   /* solve covering problem */
   retcode = SCIPsolve(coveringscip);

   /* errors in solving the covering problem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(coveringscip, "Error while solving covering problem in Undercover heuristic; sub-SCIP terminated with code <%d>\n",retcode);
      return SCIP_OKAY;
   }

   /* check, whether a feasible cover was found */
   if( SCIPgetNSols(coveringscip) == 0 )
      return SCIP_OKAY;

   /* store solution */
   SCIP_CALL( SCIPallocBufferArray(coveringscip, &solvals, ncoveringvars) );
   SCIP_CALL( SCIPgetSolVals(coveringscip, SCIPgetBestSol(coveringscip), ncoveringvars, coveringvars, solvals) );

   *coversize = 0;
   totalpenalty = 0.0;
   for( i = 0; i < ncoveringvars; i++ )
   {
      if( solvals[i] > 0.5 )
      {
         cover[*coversize] = i;
         (*coversize)++;
      }
      totalpenalty += SCIPvarGetObj(coveringvars[i]);
   }

   /* print solution if we are in SCIP's debug mode */
   assert(SCIPgetBestSol(coveringscip) != NULL);
   SCIPdebugMsg(coveringscip, "found a feasible cover: %d/%d variables fixed, normalized penalty=%g\n\n",
      *coversize, SCIPgetNOrigVars(coveringscip), SCIPgetSolOrigObj(coveringscip, SCIPgetBestSol(coveringscip))/(totalpenalty+SCIPsumepsilon(coveringscip)));
   SCIPdebug( SCIP_CALL( SCIPprintSol(coveringscip, SCIPgetBestSol(coveringscip), NULL, FALSE) ) );
   SCIPdebugMsg(coveringscip, "\r                                                  \n");

   *success = TRUE;

   /* free array of solution values */
   SCIPfreeBufferArray(coveringscip, &solvals);

   return SCIP_OKAY;
}


/** computes fixing order and returns whether order has really been changed */
static
SCIP_RETCODE computeFixingOrder(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   nvars,              /**< number of variables in the original problem */
   SCIP_VAR**            vars,               /**< variables in the original problem */
   int                   coversize,          /**< size of the cover */
   int*                  cover,              /**< problem indices of the variables in the cover */
   int                   lastfailed,         /**< position in cover array of the variable the fixing of which yielded
                                              *   infeasibility in last dive (or >= coversize, in which case *success
                                              *   is always TRUE) */
   SCIP_Bool*            success             /**< has order really been changed? */
   )
{
   SCIP_Real* scores;
   SCIP_Real bestscore;
   SCIP_Bool sortdown;
   int i;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(nvars >= 1);
   assert(vars != NULL);
   assert(coversize >= 1);
   assert(cover != NULL);
   assert(lastfailed >= 0);

   *success = FALSE;

   /* if fixing the first variable had failed, do not try with another order */
   if( lastfailed == 0 )
      return SCIP_OKAY;

   /* allocate buffer array for score values */
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, coversize) );

   /* initialize best score to infinite value */
   sortdown = (heurdata->fixingorder == 'c' || heurdata->fixingorder == 'v' );
   bestscore = sortdown ? -SCIPinfinity(scip) : +SCIPinfinity(scip);

   /* compute score values */
   for( i = coversize-1; i >= 0; i-- )
   {
      SCIP_VAR* var;

      /* get variable in the cover */
      assert(cover[i] >= 0);
      assert(cover[i] < nvars);
      var = vars[cover[i]];

      if( heurdata->fixingorder == 'C' || heurdata->fixingorder == 'c' )
      {
         /* add a small pertubation value to the score to reduce performance variability */
         scores[i] = heurdata->conflictweight * SCIPgetVarConflictScore(scip, var)
            + heurdata->inferenceweight * SCIPgetVarAvgInferenceCutoffScore(scip, var, heurdata->cutoffweight)
            + SCIPrandomGetReal(heurdata->randnumgen, 0.0, SCIPepsilon(scip));
      }
      else if( heurdata->fixingorder == 'V' || heurdata->fixingorder == 'v' )
         scores[i] = cover[i];
      else
         return SCIP_PARAMETERWRONGVAL;

      assert(scores[i] >= 0.0);

      /* update best score */
      if( sortdown )
         bestscore = MAX(bestscore, scores[i]);
      else
         bestscore = MIN(bestscore, scores[i]);

   }

   /* put integers to the front */
   if( heurdata->fixintfirst )
   {
      for( i = coversize-1; i >= 0; i-- )
      {
         if( SCIPvarIsIntegral(vars[cover[i]]) )
         {
            if( sortdown )
               scores[i] += bestscore+1.0;
            else
               scores[i] = bestscore - 1.0/(scores[i]+1.0);
         }
      }
   }

   /* put last failed variable to the front */
   if( lastfailed < coversize )
   {
      if( sortdown )
         scores[lastfailed] += bestscore+2.0;
      else
         scores[lastfailed] = bestscore - 2.0/(scores[lastfailed]+1.0);
      i = cover[lastfailed];
   }

   /* sort by non-decreasing (negative) score */
   if( sortdown )
      SCIPsortDownRealInt(scores, cover, coversize);
   else
      SCIPsortRealInt(scores, cover, coversize);

   assert(lastfailed >= coversize || cover[0] == i);

   /* free buffer memory */
   SCIPfreeBufferArray(scip, &scores);

   *success = TRUE;

   return SCIP_OKAY;
}


/** gets fixing value */
static
SCIP_RETCODE getFixingValue(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VAR*             var,                /**< variable in the original problem */
   SCIP_Real*            val,                /**< buffer for returning fixing value */
   char                  fixalt,             /**< source of the fixing value: 'l'p relaxation, 'n'lp relaxation, 'i'ncumbent solution */
   SCIP_Bool*            success,            /**< could value be retrieved successfully? */
   int                   bdlen,              /**< current length of probing path */
   SCIP_VAR**            bdvars,             /**< array of variables with bound changes along probing path */
   SCIP_BOUNDTYPE*       bdtypes,            /**< array of bound types in bound disjunction */
   SCIP_Real*            oldbounds           /**< array of bounds before fixing */
   )
{
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(var != NULL);
   assert(val != NULL);
   assert(success != NULL);

   *success = FALSE;

   switch( fixalt )
   {
   case 'l':
      /* get the last LP solution value */
      *val = SCIPvarGetLPSol(var);
      *success = TRUE;
      break;
   case 'n':
      /* only call this function if NLP relaxation is available */
      assert(SCIPisNLPConstructed(scip));

      /* the solution values are already available */
      if( heurdata->nlpsolved )
      {
         assert(!heurdata->nlpfailed);

         /* retrieve NLP solution value */
         *val = SCIPvarGetNLPSol(var);
         *success = TRUE;
      }
      /* solve NLP relaxation unless it has not failed too often before */
      else if( !heurdata->nlpfailed )
      {  /*lint --e{850}*/
         SCIP_NLPSOLSTAT stat;
         int i;

         /* restore bounds at start of probing, since otherwise, if in backtrack mode, NLP solver becomes most likely
          * locally infeasible 
          */
         SCIP_CALL( SCIPstartDiveNLP(scip) );

         for( i = bdlen-1; i >= 0; i-- )
         {
            SCIP_VAR* relaxvar;
            SCIP_Real lb;
            SCIP_Real ub;

            relaxvar = bdvars[i];

            /* both bounds were tightened */
            if( i > 0 && bdvars[i-1] == relaxvar )
            {
               assert(bdtypes[i] != bdtypes[i-1]);

               lb = bdtypes[i] == SCIP_BOUNDTYPE_UPPER ? oldbounds[i] : oldbounds[i-1];
               ub = bdtypes[i] == SCIP_BOUNDTYPE_UPPER ? oldbounds[i-1] : oldbounds[i];
               i--;
            }
            /* lower bound was tightened */
            else if( bdtypes[i] == SCIP_BOUNDTYPE_UPPER )
            {
               lb = oldbounds[i];
               ub = SCIPvarGetUbLocal(relaxvar);
            }
            /* upper bound was tightened */
            else
            {
               lb = SCIPvarGetLbLocal(relaxvar);
               ub = oldbounds[i];
            }

            assert(SCIPisLE(scip, lb, SCIPvarGetLbLocal(relaxvar)));
            assert(SCIPisGE(scip, ub, SCIPvarGetUbLocal(relaxvar)));

            /* relax bounds */
            SCIP_CALL( SCIPchgVarBoundsDiveNLP(scip, relaxvar, lb, ub) );
         }

         /* activate NLP solver output if we are in SCIP's debug mode */
         SCIPdebug( SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, 1) ) );

         /* set starting point to lp solution */
         SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );

         /* solve NLP relaxation */
         SCIP_CALL( SCIPsolveDiveNLP(scip) );
         stat = SCIPgetNLPSolstat(scip);
         *success = stat == SCIP_NLPSOLSTAT_GLOBOPT || stat == SCIP_NLPSOLSTAT_LOCOPT || stat == SCIP_NLPSOLSTAT_FEASIBLE;

         SCIPdebugMsg(scip, "solving NLP relaxation to obtain fixing values %s (stat=%d)\n", *success ? "successful" : "failed", stat);

         if( *success )
         {
            /* solving succeeded */
            heurdata->nnlpfails = 0;
            heurdata->nlpsolved = TRUE;

            /* retrieve NLP solution value */
            *val = SCIPvarGetNLPSol(var);
         }
         else
         {
            /* solving failed */
            heurdata->nnlpfails++;
            heurdata->nlpfailed = TRUE;
            heurdata->nlpsolved = FALSE;

            SCIPdebugMsg(scip, "solving NLP relaxation failed (%d time%s%s)\n",
               heurdata->nnlpfails, heurdata->nnlpfails > 1 ? "s" : "", heurdata->nnlpfails >= MAXNLPFAILS ? ", will not be called again" : "");
         }
      }
      break;
   case 'i':
      /* only call this function if a feasible solution is available */
      assert(SCIPgetBestSol(scip) != NULL);

      /* get value in the incumbent solution */
      *val = SCIPgetSolVal(scip, SCIPgetBestSol(scip), var);
      *success = TRUE;
      break;
   default:
      break;
   }

   /* due to propagation (during probing) it might happen that the LP and NLP solution value of var might be outside of
    * its bounds
    */
   *val = MAX(*val, SCIPvarGetLbLocal(var)); /*lint !e666*/
   *val = MIN(*val, SCIPvarGetUbLocal(var)); /*lint !e666*/

   return SCIP_OKAY;
}


/** calculates up to four alternative values for backtracking, if fixing the variable failed.
 * The alternatives are the two bounds of the variable, and the averages of the bounds and the fixing value.
 * For infinite bounds, fixval +/- abs(fixval) will be used instead.
 */
static
void calculateAlternatives(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR*             var,                /**< variable to calculate alternatives for */
   SCIP_Real             fixval,             /**< reference fixing value */
   int*                  nalternatives,      /**< number of fixing values computed */
   SCIP_Real*            alternatives        /**< array to store the alternative fixing values */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   /* for binary variables, there is only two possible fixing values */
   if( SCIPvarIsBinary(var) )
   {
      if( SCIPisFeasEQ(scip, fixval, 0.0) || SCIPisFeasEQ(scip, fixval, 1.0) )
      {
         alternatives[0] = 1.0 - fixval;
         *nalternatives = 1;
      }
      else
      {
         alternatives[0] = 0.0;
         alternatives[1] = 1.0;
         *nalternatives = 2;
      }
      return;
   }

   /* get bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   /* if lower bound is infinite, use x'-|x'|; if x' is zero, use -1.0 instead */
   if( SCIPisInfinity(scip, -lb) )
      lb = SCIPisFeasZero(scip, fixval) ? -1.0 : fixval - ABS(fixval);

   /* if upper bound is infinite, use x'+|x'|; if x' is zero, use 1.0 instead */
   if( SCIPisInfinity(scip, ub) )
      ub = SCIPisFeasZero(scip, fixval) ? 1.0 : fixval + ABS(fixval);

   assert(!SCIPisEQ(scip, lb, ub));

   /* collect alternatives */
   *nalternatives = 0;

   /* use lower bound if not equal to x' */
   if( !SCIPisFeasEQ(scip, lb, fixval) )
   {
      alternatives[*nalternatives] = lb;
      (*nalternatives)++;
   }

   /* use upper bound if not equal to x' */
   if( !SCIPisFeasEQ(scip, ub, fixval) )
   {
      alternatives[*nalternatives] = ub;
      (*nalternatives)++;
   }

   /* use the average of x' and lower bound as alternative value, if this is not equal to any of the other values */
   if( !SCIPisFeasEQ(scip, lb, fixval) && (!SCIPvarIsIntegral(var) || !SCIPisFeasEQ(scip, lb, fixval-1)) )
   {
      alternatives[*nalternatives] = (lb+fixval)/2.0;
      (*nalternatives)++;
   }

   /* use the average of x' and upper bound as alternative value, if this is not equal to any of the other values */
   if( !SCIPisFeasEQ(scip, ub, fixval) && (!SCIPvarIsIntegral(var) || !SCIPisFeasEQ(scip, ub, fixval+1)) )
   {
      alternatives[*nalternatives] = (ub+fixval)/2.0;
      (*nalternatives)++;
   }

   assert(*nalternatives <= 4);

   return;
}


/** rounds the given fixing value */
static
SCIP_RETCODE roundFixingValue(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_Real*            val,                /**< fixing value to be rounded */
   SCIP_VAR*             var,                /**< corresponding variable */
   SCIP_Bool             locksrounding       /**< shall we round according to locks? (otherwise to nearest integer) */
   )
{
   SCIP_Real x;

   x = *val;

   /* if integral within feasibility tolerance, only shift to nearest integer */
   if( SCIPisFeasIntegral(scip, x) )
      x = SCIPfeasFrac(scip, x) < 0.5 ? SCIPfeasFloor(scip, x) : SCIPfeasCeil(scip, x);

   /* round in the direction of least locks with fractionality as tie breaker */
   else if( locksrounding )
   {
      if( SCIPvarGetNLocksDown(var) < SCIPvarGetNLocksUp(var) )
         x = SCIPfeasFloor(scip, x);
      else if( SCIPvarGetNLocksDown(var) > SCIPvarGetNLocksUp(var) )
         x = SCIPfeasCeil(scip, x);
      else
         x = SCIPfeasFrac(scip, x) < 0.5 ? SCIPfeasFloor(scip, x) : SCIPfeasCeil(scip, x);
   }
   /* round in the direction of least fractionality with locks as tie breaker */
   else
   {
      if( SCIPfeasFrac(scip, x) < 0.5)
         x = SCIPfeasFloor(scip, x);
      else if( SCIPfeasFrac(scip, x) > 0.5 )
         x = SCIPfeasCeil(scip, x);
      else
         x = SCIPvarGetNLocksDown(var) < SCIPvarGetNLocksUp(var) ? SCIPfeasFloor(scip, x) : SCIPfeasCeil(scip, x);
   }

   /* return rounded fixing value */
   *val = x;

   return SCIP_OKAY;
}


/** copy the solution of the subproblem to newsol */
static
SCIP_RETCODE copySol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_SOL**            newsol              /**< solution to the original problem */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem */
   int        nvars;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);
   assert(newsol != NULL);
   assert(*newsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );
   SCIP_CALL( SCIPsetSolVals(scip, *newsol, nvars, vars, subsolvals) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/** solve subproblem and pass best feasible solution to original SCIP instance */
static
SCIP_RETCODE solveSubproblem(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   int                   coversize,          /**< size of the cover */
   int*                  cover,              /**< problem indices of the variables in the cover */
   SCIP_Real*            fixedvals,          /**< fixing values for the variables in the cover */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Longint          nodelimit,          /**< node limit */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_Bool*            validsolved,        /**< was problem constructed from a valid copy and solved (proven optimal or infeasible)? */
   SCIP_SOL**            sol,                /**< best solution found in subproblem (if feasible); *sol must be NULL, solution will be created */
   SCIP_Longint*         nusednodes          /**< number of nodes used for solving the subproblem */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP* subscip;
   SCIP_VAR** subvars;
   SCIP_VAR** vars;
   SCIP_HASHMAP* varmap;
   SCIP_VAR** fixedvars;
   int nfixedvars;

   SCIP_RETCODE retcode;

   int nvars;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(cover != NULL);
   assert(fixedvals != NULL);
   assert(coversize >= 1);
   assert(timelimit > 0.0);
   assert(memorylimit > 0.0);
   assert(nodelimit >= 1);
   assert(nstallnodes >= 1);
   assert(validsolved != NULL);
   assert(sol != NULL);
   assert(*sol == NULL);
   assert(nusednodes != NULL);

   *validsolved = FALSE;
   *nusednodes = 0;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, coversize) );
   nfixedvars = coversize;
   /* fix subproblem variables in the cover */
   SCIPdebugMsg(scip, "fixing variables\n");
   for( i = coversize-1; i >= 0; i-- )
   {
      assert(cover[i] >= 0);
      assert(cover[i] < nvars);

      fixedvars[i] = vars[cover[i]];
   }

   /* create subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), nvars) );

   /* copy original problem to subproblem; do not copy pricers */
   SCIP_CALL( SCIPcopyConsCompression(scip, subscip, varmap, NULL, "undercoversub", fixedvars, fixedvals, nfixedvars,
         heurdata->globalbounds, FALSE, TRUE, validsolved) );

   if( heurdata->copycuts )
   {
      /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
      SCIP_CALL( SCIPcopyCuts(scip, subscip, varmap, NULL, heurdata->globalbounds, NULL) );
   }

   SCIPdebugMsg(scip, "problem copied, copy %svalid\n", *validsolved ? "" : "in");

   /* store subproblem variables */
   for( i = nvars-1; i >= 0; i-- )
   {
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);
      assert(subvars[i] != NULL);
   }

   /* free variable mapping hash map */
   SCIPhashmapFree(&varmap);

   /* set the parameters such that good solutions are found fast */
   SCIPdebugMsg(scip, "setting subproblem parameters\n");
   SCIP_CALL( SCIPsetEmphasis(subscip, SCIP_PARAMEMPHASIS_FEASIBILITY, TRUE) );
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(subscip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );

   /* deactivate expensive pre-root heuristics, since it may happen that the lp relaxation of the subproblem is already
    * infeasible; in this case, we do not want to waste time on heuristics before solving the root lp */
   if( !SCIPisParamFixed(subscip, "heuristics/shiftandpropagate/freq") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/shiftandpropagate/freq", -1) );
   }

   /* forbid recursive call of undercover heuristic */
   if( SCIPisParamFixed(subscip, "heuristics/" HEUR_NAME "/freq") )
   {
      SCIPwarningMessage(scip, "unfixing parameter heuristics/" HEUR_NAME "/freq in subscip of undercover heuristic to avoid recursive calls\n");
      SCIP_CALL( SCIPunfixParam(subscip, "heuristics/" HEUR_NAME "/freq") );
   }
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/" HEUR_NAME "/freq", -1) );

   SCIPdebugMsg(scip, "timelimit = %g, memlimit = %g, nodelimit = %" SCIP_LONGINT_FORMAT ", nstallnodes = %" SCIP_LONGINT_FORMAT "\n", timelimit, memorylimit, nodelimit, nstallnodes);

   SCIPdebugMsg(scip, "timelimit = %g, memlimit = %g, nodelimit = %" SCIP_LONGINT_FORMAT ", nstallnodes = %" SCIP_LONGINT_FORMAT "\n", timelimit, memorylimit, nodelimit, nstallnodes);

   /* disable statistic timing inside sub SCIP */
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

   /* set time, memory and node limits */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console in optimized mode, enable in SCIP's debug mode */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000) );
#else
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

   /* if there is already a solution, add an objective cutoff; note: this does not affect the validity of the subproblem
    * if we find solutions later, thus we do not set *validsolved to FALSE */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real cutoff;
      SCIP_Real upperbound;

      assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));
      upperbound = SCIPgetUpperbound(scip);

      if( SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) )
         cutoff = (upperbound >= 0 ? 1.0 - heurdata->minimprove : 1.0 + heurdata->minimprove) * upperbound;
      else
         cutoff = (1.0 - heurdata->minimprove) * upperbound + heurdata->minimprove * SCIPgetLowerbound(scip);

      cutoff = MIN(upperbound, cutoff);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );

      SCIPdebugMsg(scip, "adding objective cutoff=%g (minimprove=%g)\n", cutoff, heurdata->minimprove);
   }

   /* solve subproblem */
   SCIPdebugMsg(scip, "solving subproblem started\n");
   retcode = SCIPsolve(subscip);

   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in Undercover heuristic; sub-SCIP terminated with code <%d>\n",retcode);
      /* free array of subproblem variables, and subproblem */
      SCIPfreeBufferArray(scip, &subvars);
      SCIPfreeBufferArray(scip, &fixedvars);
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_OKAY;
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   /* store solving status; note: if we proved infeasibility in presence of an objective cutoff beyond the primal bound,
    * the subproblem was not a valid copy */
   *validsolved = *validsolved && (SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL
      || (SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE && (SCIPgetNSols(scip) == 0 || heurdata->minimprove <= 0.0)));
   *nusednodes = SCIPgetNNodes(subscip);

   /* if a solution was found for the subproblem, create corresponding solution in the original problem */
   if( SCIPgetNSols(subscip) > 0 && (SCIPgetStatus(subscip) != SCIP_STATUS_INFEASIBLE || heurdata->minimprove > 0.0) )
   {
      SCIP_SOL** subsols;
      SCIP_Bool success;
      int nsubsols;

      /* create solution */
      SCIP_CALL( SCIPcreateSol(scip, sol, heur) );

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      assert(subsols != NULL);

      success = FALSE;
      for( i = 0; i < nsubsols && !success; i++ )
      {
         /* transform solution to original problem */
         SCIP_CALL( copySol(scip, subscip, subvars, subsols[i], sol) );

         /* try to add new solution to scip */
         SCIP_CALL( SCIPtrySol(scip, *sol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
      }

      if( success )
      {
         assert(i >= 1);
         SCIPdebugMsg(scip, "heuristic found %d solutions in subproblem; solution %d feasible in original problem\n", nsubsols, i);
      }
      else
      {
         /* free solution structure, since we found no feasible solution */
         SCIP_CALL( SCIPfreeSol(scip, sol) );
         *sol = NULL;
      }

      /* if the best subproblem solution was not accepted in the original problem, we do not trust the solving status */
      *validsolved = *validsolved && i == 1;
   }

   if( *validsolved )
   {
      SCIP_CALL( SCIPmergeVariableStatistics(subscip, scip, subvars, vars, nvars) );
   }

   /* free array of subproblem variables, and subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIPfreeBufferArray(scip, &fixedvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/** perform fixing of a variable and record bound disjunction information */
static
SCIP_RETCODE performFixing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to fix */
   SCIP_Real             val,                /**< fixing value */
   SCIP_Bool*            infeas,             /**< pointer to store whether the fixing lead to infeasibility */
   int*                  bdlen,              /**< current length of bound disjunction */
   SCIP_VAR**            bdvars,             /**< array of variables in bound disjunction */
   SCIP_BOUNDTYPE*       bdtypes,            /**< array of bound types in bound disjunction */
   SCIP_Real*            bdbounds,           /**< array of bounds in bound disjunction */
   SCIP_Real*            oldbounds           /**< array of bounds before fixing */
   )
{
   SCIP_Longint ndomredsfound;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   int oldbdlen;

   assert(scip != NULL);
   assert(var != NULL);
   assert(val >= SCIPvarGetLbLocal(var));
   assert(val <= SCIPvarGetUbLocal(var));
   assert(infeas != NULL);
   assert(bdlen != NULL);
   assert(*bdlen >= 0);
   assert(*bdlen < 2*SCIPgetNVars(scip)-1);
   assert(bdvars != NULL);
   assert(bdtypes != NULL);
   assert(bdbounds != NULL);

   assert(!SCIPvarIsIntegral(var) || SCIPisFeasIntegral(scip, val));

   /* remember length of probing path */
   oldbdlen = *bdlen;

   /* get bounds of the variable to fix */
   oldlb = SCIPvarGetLbLocal(var);
   oldub = SCIPvarGetUbLocal(var);

   /* decrease upper bound to fixing value */
   *infeas = FALSE;
   if( SCIPisUbBetter(scip, val, oldlb, oldub) )
   {
      /* we only want to open a new probing node if we do not exceed the maximal tree depth */
      if( SCIPgetDepth(scip) < SCIP_MAXTREEDEPTH )
      {
         /* create next probing node */
         SCIP_CALL( SCIPnewProbingNode(scip) );
      }
      SCIP_CALL( SCIPchgVarUbProbing(scip, var, val) );

      SCIPdebugMsg(scip, "tentatively decreasing upper bound of variable <%s> to %g for probing\n",
         SCIPvarGetName(var), val);

      /* store bound disjunction information */
      bdvars[*bdlen] = var;
      bdtypes[*bdlen] = SCIP_BOUNDTYPE_LOWER;
      bdbounds[*bdlen] = SCIPvarIsIntegral(var) ? SCIPfeasCeil(scip, val)+1.0 : val;
      oldbounds[*bdlen] = oldub;
      (*bdlen)++;

      /* propagate the bound change; conflict analysis is performed automatically */
      SCIP_CALL( SCIPpropagateProbing(scip, 0, infeas, &ndomredsfound) );
      SCIPdebugMsg(scip, "  --> propagation reduced %" SCIP_LONGINT_FORMAT " further domains\n", ndomredsfound);

      /* if propagation led to a cutoff, we backtrack immediately */
      if( *infeas )
      {
         *bdlen = oldbdlen;
         return SCIP_OKAY;
      }

      /* store bound before propagation */
      oldbounds[*bdlen] = oldlb;

      /* move fixing value into the new domain, since it may be outside due to numerical issues or previous propagation */
      oldlb = SCIPvarGetLbLocal(var);
      oldub = SCIPvarGetUbLocal(var);
      val = MIN(val, oldub);
      val = MAX(val, oldlb);

      assert(!SCIPvarIsIntegral(var) || SCIPisFeasIntegral(scip, val));
   }

   /* update lower bound to fixing value */
   *infeas = FALSE;
   if( SCIPisLbBetter(scip, val, oldlb, oldub) )
   {
      /* we only want to open a new probing node if we do not exceed the maximal tree depth */
      if( SCIPgetDepth(scip) < SCIP_MAXTREEDEPTH )
      {
         /* create next probing node */
         SCIP_CALL( SCIPnewProbingNode(scip) );
      }
      SCIP_CALL( SCIPchgVarLbProbing(scip, var, val) );

      SCIPdebugMsg(scip, "tentatively increasing lower bound of variable <%s> to %g for probing\n",
         SCIPvarGetName(var), val);

      /* store bound disjunction information */
      bdvars[*bdlen] = var;
      bdtypes[*bdlen] = SCIP_BOUNDTYPE_UPPER;
      bdbounds[*bdlen] = SCIPvarIsIntegral(var) ? SCIPfeasCeil(scip, val)-1.0 : val;
      (*bdlen)++;

      /* propagate the bound change */
      SCIP_CALL( SCIPpropagateProbing(scip, 0, infeas, &ndomredsfound) );
      SCIPdebugMsg(scip, "  --> propagation reduced %" SCIP_LONGINT_FORMAT " further domains\n", ndomredsfound);

      /* if propagation led to a cutoff, we backtrack immediately */
      if( *infeas )
      {
         *bdlen = oldbdlen;
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE fixAndPropagate(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   int*                  cover,              /**< array with indices of the variables in the computed cover */
   int                   coversize,          /**< size of the cover */
   SCIP_Real*            fixingvals,         /**< fixing values for the variables in the cover */
   int*                  bdlen,              /**< current length of bound disjunction along the probing path */
   SCIP_VAR**            bdvars,             /**< array of variables in bound disjunction */
   SCIP_BOUNDTYPE*       bdtypes,            /**< array of bound types in bound disjunction */
   SCIP_Real*            bdbounds,           /**< array of bounds in bound disjunction */
   SCIP_Real*            oldbounds,          /**< array of bounds before fixing */
   int*                  nfixedints,         /**< pointer to store number of fixed integer variables */
   int*                  nfixedconts,        /**< pointer to store number of fixed continuous variables */
   int*                  lastfailed,         /**< position in cover array of the variable the fixing of which yielded
                                              *   infeasibility */
   SCIP_Bool*            infeas              /**< pointer to store whether fix-and-propagate led to an infeasibility */
   )
{
   SCIP_VAR** vars;                          /* original problem's variables */

   int i;
   SCIP_Bool lpsolved;

   /* start probing in original problem */
   lpsolved = SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL;
   SCIP_CALL( SCIPstartProbing(scip) );

   /* initialize data */
   *nfixedints = 0;
   *nfixedconts = 0;
   *bdlen = 0;
   vars = SCIPgetVars(scip);

   /* round-fix-propagate-analyze-backtrack for each variable in the cover */
   for( i = 0; i < coversize && !(*infeas); i++ )
   {
      SCIP_Real* boundalts;
      SCIP_Real* usedvals;
      SCIP_Real val;
      int nbacktracks;
      int nboundalts;
      int nfailedvals;
      int nusedvals;
      int probingdepth;
      int idx;

      /* get probindex of next variable in the cover */
      idx = cover[i];

      /* nothing to do if the variable was already fixed, e.g., by propagation */
      if( SCIPisEQ(scip, SCIPvarGetLbLocal(vars[idx]), SCIPvarGetUbLocal(vars[idx])) )
      {
         fixingvals[i] = SCIPvarGetLbLocal(vars[idx]);
         continue;
      }

      /* we will store the fixing values already used to avoid try the same value twice */
      SCIP_CALL( SCIPallocBufferArray(scip, &usedvals, heurdata->maxbacktracks+1) );
      nusedvals = 0;

      /* backtracking loop */
      *infeas = TRUE;
      nfailedvals = 0;
      nboundalts = 0;
      boundalts = NULL;
      val = 0.0;
      for( nbacktracks = 0; nbacktracks <= heurdata->maxbacktracks+nfailedvals && *infeas; nbacktracks++ )
      {
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Bool usedbefore;
         int j;

         probingdepth = SCIPgetProbingDepth(scip);

         /* get fixing value */
         if( nbacktracks < heurdata->nfixingalts )
         {
            SCIP_Bool success;

            /* if the lp relaxation is not solved, we do not even try to retrieve the lp solution value;
             * if the NLP relaxation is not constructed, we do not even try to retrieve the NLP solution value;
             * if there is no feasible solution yet, we do not even try to obtain the value in the incumbent */
            success = FALSE;
            if( (heurdata->fixingalts[nbacktracks] != 'l' || lpsolved)
               && (heurdata->fixingalts[nbacktracks] != 'n' || !heurdata->nlpfailed)
               && (heurdata->fixingalts[nbacktracks] != 'i' || SCIPgetBestSol(scip) != NULL) )
            {
               SCIP_CALL( getFixingValue(scip, heurdata, vars[idx], &val, heurdata->fixingalts[nbacktracks], &success, *bdlen, bdvars, bdtypes, oldbounds) );
            }

            if( !success )
            {
               SCIPdebugMsg(scip, "retrieving fixing value '%c' for variable <%s> failed, trying next in the list\n",
                  heurdata->fixingalts[nbacktracks], SCIPvarGetName(vars[idx]));
               nfailedvals++;
               continue;
            }

            /* for the first (successfully retrieved) fixing value, compute (at most 4) bound dependent
             * alternative fixing values */
            if( boundalts == NULL )
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &boundalts, 4) );
               nboundalts = 0;
               calculateAlternatives(scip, vars[idx], val, &nboundalts, boundalts);
               assert(nboundalts >= 0);
               assert(nboundalts <= 4);
            }
         }
         /* get alternative fixing value */
         else if( boundalts != NULL && nbacktracks <  heurdata->nfixingalts+nboundalts )
         {
            assert(nbacktracks-heurdata->nfixingalts >= 0);
            val = boundalts[nbacktracks-heurdata->nfixingalts];
         }
         else
            break;

         /* round fixing value */
         if( SCIPvarIsIntegral(vars[idx]) && !SCIPisIntegral(scip, val) )
         {
            SCIP_CALL( roundFixingValue(scip, &val, vars[idx], heurdata->locksrounding) );
            assert(SCIPisIntegral(scip, val));
         }

         /* move value into the domain, since it may be outside due to numerical issues or previous propagation */
         oldlb = SCIPvarGetLbLocal(vars[idx]);
         oldub = SCIPvarGetUbLocal(vars[idx]);
         val = MIN(val, oldub);
         val = MAX(val, oldlb);

         assert(!SCIPvarIsIntegral(vars[idx]) || SCIPisFeasIntegral(scip, val));

         /* check if this fixing value was already used */
         usedbefore = FALSE;
         for( j = nusedvals-1; j >= 0 && !usedbefore; j-- )
            usedbefore = SCIPisFeasEQ(scip, val, usedvals[j]);

         if( usedbefore )
         {
            nfailedvals++;
            continue;
         }

         /* store fixing value */
         assert(nusedvals < heurdata->maxbacktracks);
         usedvals[nusedvals] = val;
         nusedvals++;

         /* fix-propagate-analyze */
         SCIP_CALL( performFixing(scip, vars[idx], val, infeas, bdlen, bdvars, bdtypes, bdbounds, oldbounds) );

         /* if infeasible, backtrack and try alternative fixing value */
         if( *infeas )
         {
            SCIPdebugMsg(scip, "  --> cutoff detected - backtracking\n");
            SCIP_CALL( SCIPbacktrackProbing(scip, probingdepth) );
         }
      }

      /* free array of alternative backtracking values */
      if( boundalts != NULL)
         SCIPfreeBufferArray(scip, &boundalts);
      SCIPfreeBufferArray(scip, &usedvals);

      /* backtracking loop unsuccessful */
      if( *infeas )
      {
         SCIPdebugMsg(scip, "no feasible fixing value found for variable <%s> in fixing order\n",
            SCIPvarGetName(vars[idx]));
         break;
      }
      /* fixing successful */
      else
      {
         /* store successful fixing value */
         fixingvals[i] = val;

         /* statistics */
         if( SCIPvarGetType(vars[idx]) == SCIP_VARTYPE_CONTINUOUS )
            (*nfixedconts)++;
         else
            (*nfixedints)++;
      }
   }
   assert(*infeas || i == coversize);
   assert(!(*infeas) || i < coversize);

   /* end of dive */
   SCIP_CALL( SCIPendProbing(scip) );

   *lastfailed = i;

   return SCIP_OKAY;
}

/** main procedure of the undercover heuristic */
static
SCIP_RETCODE SCIPapplyUndercover(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Longint          nstallnodes         /**< number of stalling nodes for the subproblem */
   )
{
   SCIP_HEURDATA* heurdata;                  /* heuristic data */
   SCIP_VAR** vars;                          /* original problem's variables */
   SCIP_CLOCK* clock;                        /* clock for updating time limit */

   SCIP* coveringscip;                       /* SCIP data structure for covering problem */
   SCIP_VAR** coveringvars;                  /* covering variables */
   SCIP_Real* fixingvals;                    /* array for storing fixing values used */
   int* cover;                               /* array to store problem indices of variables in the computed cover */

   SCIP_VAR** bdvars;                        /* array of variables in bound disjunction along the probing path */
   SCIP_BOUNDTYPE* bdtypes;                  /* array of bound types in bound disjunction along the probing path */
   SCIP_Real* bdbounds;                      /* array of bounds in bound disjunction along the probing path */
   SCIP_Real* oldbounds;                     /* array of bounds before fixing along the probing path */

   SCIP_Real maxcoversize;

   int coversize;
   int nvars;
   int ncovers;
   int nunfixeds;
   int nnlconss;
   int i;

   SCIP_Bool success;
   SCIP_Bool reusecover;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND);

   /* create and start timing */
   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   /* initialize */
   fixingvals = NULL;
   cover = NULL;
   bdvars = NULL;
   bdtypes = NULL;
   bdbounds = NULL;
   oldbounds = NULL;
   coversize = 0;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* NLP relaxation has not been solved yet (only solve once, not again for each cover or dive, because it is expensive) */
   heurdata->nlpsolved = FALSE;

   /* if solving the NLP relaxation has failed too often in previous runs, or NLP and NLP solver is not available, we do
    * not even try 
    */
   heurdata->nlpfailed = heurdata->nnlpfails >= MAXNLPFAILS || !SCIPisNLPConstructed(scip) || SCIPgetNNlpis(scip) == 0;

   /* get variable data of original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* get number of nonlinear constraints */
   nnlconss = 0;
   for( i = 0; i < heurdata->nnlconshdlrs; ++i )
      nnlconss += SCIPconshdlrGetNActiveConss(heurdata->nlconshdlrs[i]);
   assert(nnlconss >= 0);
   assert(nnlconss <= SCIPgetNConss(scip));

   /* run only if problem is sufficiently nonlinear */
   if( nnlconss < (SCIP_Real) SCIPgetNConss(scip) * heurdata->mincoveredrel || nnlconss < heurdata->mincoveredabs )
   {
      SCIPdebugMsg(scip, "too few nonlinear constraints present, not running\n");

      /* free clock */
      SCIP_CALL( SCIPfreeClock(scip, &clock) );

      return SCIP_OKAY;
   }

   /* calculate upper bound for cover size */
   maxcoversize = nvars*heurdata->maxcoversizevars;
   if( heurdata->maxcoversizeconss < SCIP_REAL_MAX )
   {
      SCIP_Real maxcoversizeconss;
      maxcoversizeconss = heurdata->maxcoversizeconss * nnlconss / ((SCIP_Real) SCIPgetNConss(scip));
      maxcoversize = MIN(maxcoversize, maxcoversizeconss);
   }

   /* create covering problem */
   success = FALSE;
   SCIP_CALL( SCIPcreate(&coveringscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(coveringscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coveringvars, nvars) );
   SCIP_CALL( createCoveringProblem(scip, coveringscip, coveringvars, heurdata->globalbounds, heurdata->onlyconvexify,
         heurdata->coverbd, heurdata->coveringobj, &success) );

   if( !success )
   {
      SCIPdebugMsg(scip, "creating covering problem failed, terminating\n");
      goto TERMINATE;
   }
   else
   {
      SCIPdebugMsg(scip, "covering problem created successfully\n");
   }

   /* count number of unfixed covering variables */
   nunfixeds = 0;
   for( i = nvars-1; i >= 0; i-- )
   {
      if( SCIPisFeasEQ(coveringscip, SCIPvarGetLbGlobal(coveringvars[i]), 1.0) )
         nunfixeds++;
   }

   /* update time limit */
   SCIP_CALL( updateTimelimit(scip, clock, &timelimit) );

   if( timelimit <= MINTIMELEFT )
   {
      SCIPdebugMsg(scip, "time limit hit, terminating\n");
      goto TERMINATE;
   }

   /* update memory left */
   memorylimit -= SCIPgetMemUsed(coveringscip)/1048576.0;
   memorylimit -= SCIPgetMemExternEstim(coveringscip)/1048576.0;

   /* allocate memory for storing bound disjunction information along probing path */
   SCIP_CALL( SCIPallocBufferArray(scip, &bdvars, 2*nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bdtypes, 2*nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bdbounds, 2*nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oldbounds, 2*nvars) );

   /* initialize data for recovering loop */
   SCIP_CALL( SCIPallocBufferArray(scip, &cover, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixingvals, nvars) );
   ncovers = 0;
   success = FALSE;
   reusecover = FALSE;

   heurdata->nfixingalts = (int) strlen(heurdata->fixingalts);
   assert(heurdata->nfixingalts >= 1);

   /* recovering loop */
   while( (ncovers <= heurdata->maxrecovers || reusecover) && !success )
   {
      int lastfailed;
      int ndives;
      int nfixedints;
      int nfixedconts;
      int bdlen;                                /* current length of bound disjunction along the probing path */

      SCIP_Bool conflictcreated;

      SCIPdebugMsg(scip, "solving covering problem\n\n");
      success = FALSE;
      bdlen = 0;
      conflictcreated = FALSE;

      /* solve covering problem */
      if( !reusecover )
      {
         SCIP_CALL( solveCoveringProblem(coveringscip, nvars, coveringvars, &coversize, cover,
               timelimit, memorylimit + (SCIPgetMemExternEstim(coveringscip)+SCIPgetMemUsed(coveringscip))/1048576.0, maxcoversize, &success) );

         SCIPstatistic(
            if( ncovers == 0 && success )
               SCIPstatisticPrintf("UCstats coversize abs: %6d rel: %9.6f\n", coversize, 100.0*coversize /(SCIP_Real)nvars);
            );

         assert(coversize >= 0);
         assert(coversize <= nvars);
         ncovers++;

         /* free transformed covering problem immediately */
         SCIP_CALL( SCIPfreeTransform(coveringscip) );

         /* terminate if no feasible cover was found */
         if( !success )
         {
            SCIPdebugMsg(scip, "no feasible cover found in covering problem %d, terminating\n", ncovers);
            goto TERMINATE;
         }

         /* terminate, if cover is empty or too large */
         if( coversize == 0 || coversize > maxcoversize )
         {
            SCIPdebugMsg(scip, "terminating due to coversize=%d\n", coversize);
            goto TERMINATE;
         }

         /* terminate, if cover too large for the ratio of nonlinear constraints */
         if( heurdata->maxcoversizeconss < SCIP_REAL_MAX && coversize > heurdata->maxcoversizeconss * nnlconss / (SCIP_Real) SCIPgetNConss(scip) )
         {
            SCIPdebugMsg(scip, "terminating due to coversize=%d\n", coversize);
            goto TERMINATE;
         }
      }

      /* data setup */
      ndives = 0;
      nfixedints = 0;
      nfixedconts = 0;
      success = FALSE;
      lastfailed = reusecover ? MAX(1, coversize-1) : coversize;

      /* round-fix-propagate-analyze-backtrack-reorder loop */
      while( ndives <= heurdata->maxreorders && !success )
      {
         SCIP_Bool reordered;
         SCIP_Bool infeas;

         /* compute fixing order */
         SCIP_CALL( computeFixingOrder(scip, heurdata, nvars, vars, coversize, cover, lastfailed, &reordered) );
         reordered = reordered || ndives == 0;
         SCIPdebugMsg(scip, "%sordering variables in cover %s\n", ndives == 0 ? "" : "re", reordered ? "" : "failed");

         /* stop if order has not changed */
         if( !reordered )
            break;

         infeas = FALSE;
         SCIP_CALL( fixAndPropagate(scip, heurdata, cover, coversize, fixingvals, &bdlen, bdvars, bdtypes, bdbounds, oldbounds,
               &nfixedints, &nfixedconts, &lastfailed, &infeas) );
         ndives++;
         success = !infeas;
      }

      /* update time limit */
      SCIPdebugMsg(scip, "%d dive%s of fix-and-propagate for cover %d took %.1f seconds\n", ndives, ndives > 1 ? "s" : "", ncovers, SCIPgetClockTime(scip, clock));
      SCIP_CALL( updateTimelimit(scip, clock, &timelimit) );

      if( timelimit <= MINTIMELEFT )
      {
         SCIPdebugMsg(scip, "time limit hit, terminating\n");
         goto TERMINATE;
      }

      /* no feasible fixing could be found for the current cover */
      if( !success )
      {
         SCIPdebugMsg(scip, "no feasible fixing values found for cover %d\n", ncovers);
      }
      else
      {
         SCIP_SOL* sol;
         SCIP_Longint nsubnodes;
         SCIP_Bool validsolved;

         SCIPdebugMsg(scip, "heuristic successfully fixed %d variables (%d integral, %d continuous) during probing\n",
            nfixedints+nfixedconts, nfixedints, nfixedconts); /*lint !e771*/

         /* solve sub-CIP and pass feasible solutions to original problem */
         success = FALSE;
         validsolved = FALSE;
         sol = NULL;
         nsubnodes = 0;

         SCIP_CALL( solveSubproblem(scip, heur, coversize, cover, fixingvals,
               timelimit, memorylimit, heurdata->maxnodes, nstallnodes, &validsolved, &sol, &nsubnodes) );

         /* update number of sub-CIP nodes used by heuristic so far */
         heurdata->nusednodes += nsubnodes;

         /* if the subproblem was constructed from a valid copy and solved, try to forbid the assignment of fixing
          * values to variables in the cover
          */
         if( validsolved )
         {
            SCIP_Real maxvarsfac;
            SCIP_Bool useconf;
            int minmaxvars;

            SCIP_CALL( SCIPgetIntParam(scip, "conflict/minmaxvars", &minmaxvars) );
            SCIP_CALL( SCIPgetRealParam(scip, "conflict/maxvarsfac", &maxvarsfac) );

            useconf = bdlen > 0 && (bdlen <= minmaxvars || bdlen < maxvarsfac*nvars);

            if( useconf )
            {
               /* even if we had reset the global bounds at start of probing, the constraint might be only locally valid due to local constraints/cuts */
               SCIP_CALL( createConflict(scip, bdlen, bdvars, bdtypes, bdbounds, SCIPgetDepth(scip) > 0, TRUE, TRUE, &success) );
               conflictcreated = success;
            }

            SCIPdebugMsg(scip, "subproblem solved (%s), forbidding assignment in original problem %s, %sconflict length=%d\n",
               sol == NULL ? "infeasible" : "optimal",
               success ? "successful" : "failed", useconf ? "" : "skipped due to ", bdlen);
         }

         /* heuristic succeeded */
         success = (sol != NULL);
         if( success )
         {
            *result = SCIP_FOUNDSOL;
            success = TRUE;

            /* update time limit */
            SCIP_CALL( updateTimelimit(scip, clock, &timelimit) );

            /* call NLP local search heuristic unless it has failed too often */
            if( heurdata->postnlp && heurdata->npostnlpfails < MAXPOSTNLPFAILS )
            {
               if( nfixedconts == 0 && validsolved )
               {
                  SCIPdebugMsg(scip, "subproblem solved to optimality while all covering variables are integral, hence skipping NLP local search\n");
               }
               else if( timelimit <= MINTIMELEFT )
               {
                  SCIPdebugMsg(scip, "time limit hit, skipping NLP local search\n");
               }
               else if( heurdata->nlpheur == NULL )
               {
                  SCIPdebugMsg(scip, "NLP heuristic not found, skipping NLP local search\n");
               }
               else
               {
                  SCIP_RESULT nlpresult;

                  SCIP_CALL( SCIPapplyHeurSubNlp(scip, heurdata->nlpheur, &nlpresult, sol, -1LL, timelimit, heurdata->minimprove, NULL, NULL) );
                  SCIPdebugMsg(scip, "NLP local search %s\n", nlpresult == SCIP_FOUNDSOL ? "successful" : "failed");

                  if( nlpresult == SCIP_FOUNDSOL )
                     heurdata->npostnlpfails = 0;
                  else
                     heurdata->npostnlpfails++;
               }
            }

            /* free solution */
            SCIP_CALL( SCIPfreeSol(scip, &sol) );
         }
      }

      /* heuristic failed but we have another recovering try, hence we forbid the current cover in the covering problem */
      if( !success && ncovers <= heurdata->maxrecovers )
      {
         SCIP_Bool infeas;
         int diversification;

         /* compute minimal number of unfixed covering variables (in the cover) which have to change their value */
         diversification = (int) SCIPfeasCeil(scip, (heurdata->recoverdiv) * (SCIP_Real) nunfixeds);
         diversification = MAX(diversification, 1);

         /* forbid unsuccessful cover globally in covering problem */
         SCIP_CALL( forbidCover(coveringscip, nvars, coveringvars, coversize, cover, diversification, &success, &infeas) );

         if( infeas )
         {
            SCIPdebugMsg(scip, "recovering problem infeasible (diversification=%d), terminating\n", diversification);
            goto TERMINATE;
         }
         else if( !success )
         {
            SCIPdebugMsg(scip, "failed to forbid current cover in the covering problem, terminating\n");
            goto TERMINATE;
         }
         else
         {
            SCIPdebugMsg(scip, "added constraint to the covering problem in order to forbid current cover\n");
            success = FALSE;
         }
      }

      /* try to re-use the same cover at most once */
      if( heurdata->reusecover && !reusecover && conflictcreated )
         reusecover = TRUE;
      else
         reusecover = FALSE;
   }

 TERMINATE:
   if( *result != SCIP_FOUNDSOL && *result != SCIP_DELAYED )
   {
      SCIPdebugMsg(scip, "heuristic terminating unsuccessfully\n");
   }

   /* we must remain in NLP diving mode until here to be able to retrieve NLP solution values easily */
   /* assert((SCIPisNLPConstructed(scip) == FALSE && heurdata->nlpsolved == FALSE) ||
    * (SCIPisNLPConstructed(scip) == TRUE && heurdata->nlpsolved == SCIPnlpIsDiving(SCIPgetNLP(scip))));
    */
   if( heurdata->nlpsolved )
   {
      SCIP_CALL( SCIPendDiveNLP(scip) );
   }

   /* free arrays for storing the cover */
   SCIPfreeBufferArrayNull(scip, &fixingvals);
   SCIPfreeBufferArrayNull(scip, &cover);

   /* free arrays for storing bound disjunction information along probing path */
   SCIPfreeBufferArrayNull(scip, &oldbounds);
   SCIPfreeBufferArrayNull(scip, &bdbounds);
   SCIPfreeBufferArrayNull(scip, &bdtypes);
   SCIPfreeBufferArrayNull(scip, &bdvars);

   /* free covering problem */
   for( i = nvars-1; i >= 0; i-- )
   {
      SCIP_CALL( SCIPreleaseVar(coveringscip, &coveringvars[i]) );
   }
   SCIPfreeBufferArray(scip, &coveringvars);
   SCIP_CALL( SCIPfree(&coveringscip) );

   /* free clock */
   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyUndercover)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurUndercover(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic */
static
SCIP_DECL_HEUREXIT(heurExitUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   int h;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize counters to zero */
   heurdata->nusednodes = 0;
   heurdata->npostnlpfails = 0;
   heurdata->nnlpfails = 0;

   /* if the heuristic is called at the root node, we may want to be called directly after the initial root LP solve */
   if( heurdata->beforecuts && SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP);

   /* find nonlinear constraint handlers */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->nlconshdlrs, 7) );/*lint !e506*/
   h = 0;

   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "and");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;

   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "quadratic");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;

   if( heurdata->coverbd )
   {
      heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "bounddisjunction");
      if( heurdata->nlconshdlrs[h] != NULL )
         h++;
   }

   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "indicator");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;

   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "soc");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;

   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "nonlinear");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;

   heurdata->nlconshdlrs[h] = SCIPfindConshdlr(scip, "abspower");
   if( heurdata->nlconshdlrs[h] != NULL )
      h++;

   heurdata->nnlconshdlrs = h;
   assert( heurdata->nnlconshdlrs <= 7 );

   /* find NLP local search heuristic */
   heurdata->nlpheur = SCIPfindHeur(scip, "subnlp");

   /* add global linear constraints to NLP relaxation */
   if( SCIPisNLPConstructed(scip) && heurdata->nlpheur != NULL )
   {
      SCIP_CALL( SCIPaddLinearConsToNlpHeurSubNlp(scip, heurdata->nlpheur, TRUE, TRUE) );
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolUndercover)
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free array of nonlinear constraint handlers */
   SCIPfreeBlockMemoryArray(scip, &heurdata->nlconshdlrs, 7);

   /* reset timing, if it was changed temporary (at the root node) */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecUndercover)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;                  /* heuristic data */
   SCIP_Real timelimit;                      /* time limit for the subproblem */
   SCIP_Real memorylimit;                    /* memory limit for the subproblem */
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */
   SCIP_Bool run;

   int h;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only call heuristic once at the root */
   if( SCIPgetDepth(scip) == 0 && SCIPheurGetNCalls(heur) > 0 )
      return SCIP_OKAY;

   /* if we want to use NLP fixing values exclusively and no NLP solver is available, we cannot run */
   if( strcmp(heurdata->fixingalts, "n") == 0 && SCIPgetNNlpis(scip) == 0 )
   {
      SCIPdebugMsg(scip, "skipping undercover heuristic: want to use NLP fixing values exclusively, but no NLP solver available\n");
      return SCIP_OKAY;
   }

   /* calculate stallnode limit */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward heuristic if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur) + 1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= SUBMIPSETUPCOSTS * SCIPheurGetNCalls(heur);  /* account for the setup costs of the sub-CIP */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->nusednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);
   nstallnodes = MAX(nstallnodes, 1);

   /* only call heuristics if we have enough nodes left to call sub-CIP solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping undercover heuristic: nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* only call heuristics if we have enough time left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   if( timelimit <= 2*MINTIMELEFT )
   {
      SCIPdebugMsg(scip, "skipping undercover heuristic: time left=%g\n", timelimit);
      return SCIP_OKAY;
   }

   /* only call heuristics if we have enough memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   if( memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
   {
      SCIPdebugMsg(scip, "skipping undercover heuristic: too little memory\n");
      return SCIP_OKAY;
   }

   /* only call heuristic if nonlinear constraints are present */
   run = FALSE;
   for( h = heurdata->nnlconshdlrs-1; h >= 0 && !run; h-- )
   {
      run = (SCIPconshdlrGetNActiveConss(heurdata->nlconshdlrs[h]) > 0);
   }

   /* go through all nlrows and check for general nonlinearities */
   if( SCIPisNLPConstructed(scip) )
   {
      SCIP_NLROW** nlrows;
      int nnlrows;
      int i;

      /* get nlrows */
      nnlrows = SCIPgetNNLPNlRows(scip);
      nlrows = SCIPgetNLPNlRows(scip);

      /* check for an nlrow with nontrivial expression tree or quadratic terms; start from 0 since we expect the linear
       * nlrows at the end */
      for( i = nnlrows-1; i >= 0 && !run; i-- )
      {
         assert(nlrows[i] != NULL);
         run = SCIPnlrowGetExprtree(nlrows[i]) != NULL && SCIPexprtreeGetNVars(SCIPnlrowGetExprtree(nlrows[i])) > 0;
         run = run || SCIPnlrowGetNQuadVars(nlrows[i]) > 0;
      }
   }

   if( !run )
   {
      SCIPdebugMsg(scip, "skipping undercover heuristic: no nonlinear constraints found\n");
      return SCIP_OKAY;
   }

   /* only call heuristics if solving has not stopped yet */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* reset timing, if it was changed temporary (at the root node) */
   if( heurtiming != HEUR_TIMING )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);

   /* call heuristic */
   *result = SCIP_DIDNOTFIND;
   SCIPdebugMsg(scip, "calling undercover heuristic for <%s> at depth %d\n", SCIPgetProbName(scip), SCIPgetDepth(scip));

   SCIP_CALL( SCIPapplyUndercover(scip, heur, result, timelimit, memorylimit, nstallnodes) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the undercover primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurUndercover(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create undercover primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* always use local bounds */
   heurdata->globalbounds = FALSE;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecUndercover, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyUndercover) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeUndercover) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitUndercover) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitUndercover) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolUndercover) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolUndercover) );

   /* add string parameters */
   heurdata->fixingalts = NULL;
   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/" HEUR_NAME "/fixingalts",
         "prioritized sequence of fixing values used ('l'p relaxation, 'n'lp relaxation, 'i'ncumbent solution)",
         &heurdata->fixingalts, FALSE, DEFAULT_FIXINGALTS, NULL, NULL) );

   /* add longint parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   /* add real parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/conflictweight",
         "weight for conflict score in fixing order",
         &heurdata->conflictweight, TRUE, DEFAULT_CONFLICTWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/cutoffweight",
         "weight for cutoff score in fixing order",
         &heurdata->cutoffweight, TRUE, DEFAULT_CUTOFFWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/inferenceweight",
         "weight for inference score in fixing order",
         &heurdata->inferenceweight, TRUE, DEFAULT_INFERENCEWEIGHT, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxcoversizevars",
         "maximum coversize (as fraction of total number of variables)",
         &heurdata->maxcoversizevars, TRUE, DEFAULT_MAXCOVERSIZEVARS, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxcoversizeconss",
         "maximum coversize maximum coversize (as ratio to the percentage of non-affected constraints)",
         &heurdata->maxcoversizeconss, TRUE, DEFAULT_MAXCOVERSIZECONSS, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/mincoveredrel",
         "minimum percentage of nonlinear constraints in the original problem",
         &heurdata->mincoveredrel, TRUE, DEFAULT_MINCOVEREDREL, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which the heuristic should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, -1.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/recoverdiv",
         "fraction of covering variables in the last cover which need to change their value when recovering",
         &heurdata->recoverdiv, TRUE, DEFAULT_RECOVERDIV, 0.0, 1.0, NULL, NULL) );

   /* add int parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/mincoveredabs",
         "minimum number of nonlinear constraints in the original problem",
         &heurdata->mincoveredabs, TRUE, DEFAULT_MINCOVEREDABS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxbacktracks",
         "maximum number of backtracks in fix-and-propagate",
         &heurdata->maxbacktracks, TRUE, DEFAULT_MAXBACKTRACKS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxrecovers",
         "maximum number of recoverings",
         &heurdata->maxrecovers, TRUE, DEFAULT_MAXRECOVERS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxreorders",
         "maximum number of reorderings of the fixing order",
         &heurdata->maxreorders, TRUE, DEFAULT_MAXREORDERS, 0, INT_MAX, NULL, NULL) );

   /* add char parameters */
   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/coveringobj",
         "objective function of the covering problem (influenced nonlinear 'c'onstraints/'t'erms, 'd'omain size, 'l'ocks, 'm'in of up/down locks, 'u'nit penalties)",
         &heurdata->coveringobj, TRUE, DEFAULT_COVERINGOBJ, COVERINGOBJS, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/fixingorder",
         "order in which variables should be fixed (increasing 'C'onflict score, decreasing 'c'onflict score, increasing 'V'ariable index, decreasing 'v'ariable index",
         &heurdata->fixingorder, TRUE, DEFAULT_FIXINGORDER, FIXINGORDERS, NULL, NULL) );

   /* add bool parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/beforecuts",
         "should the heuristic be called at root node before cut separation?",
         &heurdata->beforecuts, TRUE, DEFAULT_BEFORECUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/fixintfirst",
         "should integer variables in the cover be fixed first?",
         &heurdata->fixintfirst, TRUE, DEFAULT_FIXINTFIRST, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/locksrounding",
         "shall LP values for integer vars be rounded according to locks?",
         &heurdata->locksrounding, TRUE, DEFAULT_LOCKSROUNDING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/onlyconvexify",
         "should we only fix variables in order to obtain a convex problem?",
         &heurdata->onlyconvexify, FALSE, DEFAULT_ONLYCONVEXIFY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/postnlp",
         "should the NLP heuristic be called to polish a feasible solution?",
         &heurdata->postnlp, FALSE, DEFAULT_POSTNLP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/coverbd",
         "should bounddisjunction constraints be covered (or just copied)?",
         &heurdata->coverbd, TRUE, DEFAULT_COVERBD, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/reusecover",
         "shall the cover be reused if a conflict was added after an infeasible subproblem?",
         &heurdata->reusecover, TRUE, DEFAULT_REUSECOVER, NULL, NULL) );

   return SCIP_OKAY;
}

/** create and solve covering problem */
static
SCIP_RETCODE computeCoverUndercover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 coveringscip,       /**< SCIP instance for covering problem */
   int*                  coversize,          /**< buffer for the size of the computed cover */
   SCIP_VAR**            cover,              /**< pointer to store the variables (of the original SCIP) in the computed cover
                                              *   (should be ready to hold SCIPgetNVars(scip) entries) */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Real             objlimit,           /**< objective limit: upper bound on coversize */
   SCIP_Bool             globalbounds,       /**< should global bounds on variables be used instead of local bounds at focus node? */
   SCIP_Bool             onlyconvexify,      /**< should we only fix/dom.red. variables creating nonconvexity? */
   SCIP_Bool             coverbd,            /**< should bounddisjunction constraints be covered (or just copied)? */
   char                  coveringobj,        /**< objective function of the covering problem ('b'ranching status,
                                              *   influenced nonlinear 'c'onstraints/'t'erms, 'd'omain size, 'l'ocks,
                                              *   'm'in of up/down locks, 'u'nit penalties, constraint 'v'iolation) */
   SCIP_Bool*            success             /**< feasible cover found? */
   )
{
   SCIP_VAR** coveringvars;                  /* covering variables */
   SCIP_VAR** vars;                          /* original variables */
   int* coverinds;                           /* indices of variables in the cover */
   int nvars;                                /* number of original variables */
   int i;

   assert(scip != NULL);
   assert(coveringscip != NULL);

   SCIP_CALL( SCIPincludeDefaultPlugins(coveringscip) );

   /* allocate memory for variables of the covering problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coveringvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coverinds, nvars) );

   SCIP_CALL( createCoveringProblem(scip, coveringscip, coveringvars, globalbounds, onlyconvexify, coverbd, coveringobj, success) );

   if( *success )
   {
      /* solve covering problem */
      SCIPdebugMsg(scip, "solving covering problem\n\n");

      SCIP_CALL( solveCoveringProblem(coveringscip, nvars, coveringvars, coversize, coverinds,
            timelimit, memorylimit + (SCIPgetMemExternEstim(coveringscip)+SCIPgetMemUsed(coveringscip))/1048576.0, objlimit, success) );

      if( *success )
      {
         assert(*coversize >= 0);
         assert(*coversize <= nvars);

         /* return original variables in the cover */
         for( i = *coversize-1; i >= 0; i-- )
         {
            assert(coverinds[i] >= 0);
            assert(coverinds[i] < nvars);
            cover[i] = vars[coverinds[i]];
         }
      }
   }
   else
   {
      SCIPdebugMsg(scip, "failure: covering problem could not be created\n");
   }

   /* free covering problem */
   for( i = nvars-1; i >= 0; i-- )
   {
      SCIP_CALL( SCIPreleaseVar(coveringscip, &coveringvars[i]) );
   }
   SCIPfreeBufferArray(scip, &coverinds);
   SCIPfreeBufferArray(scip, &coveringvars);

   return SCIP_OKAY;
}

/** computes a minimal set of covering variables */
SCIP_RETCODE SCIPcomputeCoverUndercover(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  coversize,          /**< buffer for the size of the computed cover */
   SCIP_VAR**            cover,              /**< pointer to store the variables (of the original SCIP) in the computed cover
                                              *   (should be ready to hold SCIPgetNVars(scip) entries) */
   SCIP_Real             timelimit,          /**< time limit */
   SCIP_Real             memorylimit,        /**< memory limit */
   SCIP_Real             objlimit,           /**< objective limit: upper bound on coversize */
   SCIP_Bool             globalbounds,       /**< should global bounds on variables be used instead of local bounds at focus node? */
   SCIP_Bool             onlyconvexify,      /**< should we only fix/dom.red. variables creating nonconvexity? */
   SCIP_Bool             coverbd,            /**< should bounddisjunction constraints be covered (or just copied)? */
   char                  coveringobj,        /**< objective function of the covering problem ('b'ranching status,
                                              *   influenced nonlinear 'c'onstraints/'t'erms, 'd'omain size, 'l'ocks,
                                              *   'm'in of up/down locks, 'u'nit penalties, constraint 'v'iolation) */
   SCIP_Bool*            success             /**< feasible cover found? */
   )
{
   SCIP* coveringscip;                       /* SCIP instance for covering problem */
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(coversize != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* create covering problem */
   SCIP_CALL( SCIPcreate(&coveringscip) );

   retcode = computeCoverUndercover(scip, coveringscip, coversize, cover,
         timelimit, memorylimit, objlimit,
         globalbounds, onlyconvexify, coverbd, coveringobj, success);

   /* free the covering problem scip instance before reacting on potential errors */
   SCIP_CALL( SCIPfree(&coveringscip) );

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}
