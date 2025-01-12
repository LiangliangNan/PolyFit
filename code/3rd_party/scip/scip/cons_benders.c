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

/**@file   cons_benders.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for Benders' decomposition
 * @author Stephen J. Maher
 *
 * Two constraint handlers are implemented for the generation of Benders' decomposition cuts. When included in a
 * problem, the Benders' decomposition constraint handlers generate cuts during the enforcement of LP and relaxation
 * solutions. Additionally, Benders' decomposition cuts can be generated when checking the feasibility of solutions with
 * respect to the subproblem constraints.
 *
 * This constraint handler has an enforcement priority that is less than the integer constraint handler. This means that
 * only integer feasible solutions from the LP solver are enforced by this constraint handler. This is the traditional
 * behaviour of the branch-and-check approach to Benders' decomposition. Additionally, the check priority is set low,
 * such that this expensive constraint handler is only called as a final check on primal feasible solutions.
 *
 * This constraint handler in the standard constraint handler that should be added when using Benders' decomposition.
 * Additionally, there is a flag in SCIPincludeConshdlrBenders that permits the addition of the LP constraint handler,
 * cons_benderslp. The use of both cons_benders and cons_benderslp allows the user to perform a multiphase Benders'
 * decomposition algorithm.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/cons_benders.h"
#include "scip/heur_trysol.h"
#include "scip/heuristics.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "benders"
#define CONSHDLR_DESC          "constraint handler to execute Benders' Decomposition"
#define CONSHDLR_ENFOPRIORITY      -100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -5000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_FAST /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */


#define DEFAULT_CHECKEDSOLSSIZE      20 /**< initial size of the checked sols array */
#define DEFAULT_ACTIVE            FALSE /**< is the constraint handler active? */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int*                  checkedsols;        /**< an array of solutions that this constraint has already checked */
   int                   ncheckedsols;       /**< the number of checked solutions */
   int                   checkedsolssize;    /**< the size of the checked solutions array */
   SCIP_Bool             active;             /**< is the constraint handler active? */
};

/*
 * Local methods
 */

/** constructs a new solution based upon the solutions to the Benders' decomposition subproblems */
static
SCIP_RETCODE constructValidSolution(
   SCIP*                 scip,               /**< the SCIP instance */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_BENDERSENFOTYPE  type                /**< the type of solution being enforced */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_SOL* newsol;
   SCIP_HEUR* heurtrysol;
   SCIP_BENDERS** benders;
   SCIP_VAR** auxiliaryvars;
   int nactivebenders;
   int nsubproblems;
   int i;
   int j;
   SCIP_Bool success = TRUE;

   /* don't propose new solutions if not in presolve or solving */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   benders = SCIPgetBenders(scip);
   nactivebenders = SCIPgetNActiveBenders(scip);

   /* if the solution is NULL, then we create the solution from the LP sol */
   if( sol != NULL )
   {
      assert(type == SCIP_BENDERSENFOTYPE_CHECK);
      SCIP_CALL( SCIPcreateSolCopy(scip, &newsol, sol) );
   }
   else
   {
      switch( type )
      {
      case SCIP_BENDERSENFOTYPE_LP:
         SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
         break;
      case SCIP_BENDERSENFOTYPE_PSEUDO:
         SCIP_CALL( SCIPcreatePseudoSol(scip, &newsol, NULL) );
         break;
      case SCIP_BENDERSENFOTYPE_RELAX:
         SCIP_CALL( SCIPcreateRelaxSol(scip, &newsol, NULL) );
         break;
      default:
         SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
         break;
      }  /*lint !e788*/
   }
   SCIP_CALL( SCIPunlinkSol(scip, newsol) );

   /* looping through all Benders' decompositions to construct the new solution */
   for( i = 0; i < nactivebenders; i++ )
   {
      /* getting the auxiliary variables and the number of subproblems from the Benders' decomposition structure */
      auxiliaryvars = SCIPbendersGetAuxiliaryVars(benders[i]);
      nsubproblems = SCIPbendersGetNSubproblems(benders[i]);

      /* setting the auxiliary variable in the new solution */
      for( j = 0; j < nsubproblems; j++ )
      {
         SCIP_Real objval;

         objval = SCIPbendersGetSubproblemObjval(benders[i], j);

         if( SCIPvarGetStatus(auxiliaryvars[j]) == SCIP_VARSTATUS_FIXED
            && !SCIPisEQ(scip, SCIPgetSolVal(scip, newsol, auxiliaryvars[j]), objval) )
         {
            success = FALSE;
            break;
         }
         else if( SCIPisLT(scip, SCIPgetSolVal(scip, newsol, auxiliaryvars[j]), objval) )
         {
            SCIP_CALL( SCIPsetSolVal(scip, newsol, auxiliaryvars[j], objval) );
         }
      }

      if( !success )
         break;
   }

   /* if setting the variable values was successful, then we try to add the solution */
   if( success ) /*lint !e774*/
   {
      /* checking the size of the checkedsols array and extending it is there is not enough memory */
      assert(conshdlrdata->ncheckedsols <= conshdlrdata->checkedsolssize);
      if( conshdlrdata->ncheckedsols + 1 > conshdlrdata->checkedsolssize )
      {
         int newsize;

         newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->ncheckedsols + 1);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->checkedsols, conshdlrdata->checkedsolssize, newsize) );
         conshdlrdata->checkedsolssize = newsize;
      }
      assert(conshdlrdata->ncheckedsols + 1 <= conshdlrdata->checkedsolssize);

      /* recording the solution number to avoid checking the solution again */
      conshdlrdata->checkedsols[conshdlrdata->ncheckedsols] = SCIPsolGetIndex(newsol);
      conshdlrdata->ncheckedsols++;

      /* getting the try solution heuristic */
      heurtrysol = SCIPfindHeur(scip, "trysol");

      /* passing the new solution to the trysol heuristic  */
      SCIP_CALL( SCIPcheckSol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
      if ( success )
      {
         SCIP_CALL( SCIPheurPassSolAddSol(scip, heurtrysol, newsol) );
         SCIPdebugMsg(scip, "Creating solution was successful.\n");
      }
      else
      {
         /* the solution might not be feasible, because of additional constraints */
         SCIPdebugMsg(scip, "Creating solution was not successful.\n");
      }
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/** checks the Benders' decomposition auxiliary variables for unboundedness. */
static
SCIP_Bool unboundedAuxiliaryVariables(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   SCIP_SOL*             sol                 /**< the primal solution to enforce, or NULL for the current LP/pseudo sol */
   )
{
   int nsubproblems;
   SCIP_Bool unbounded = FALSE;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* checking the auxiliary variable values for unboundedness */
   for( i = 0; i < nsubproblems; i++ )
   {
      if( SCIPisInfinity(scip, SCIPgetBendersAuxiliaryVarVal(scip, benders, sol, i))
         || SCIPisInfinity(scip, -SCIPgetBendersAuxiliaryVarVal(scip, benders, sol, i)) )
      {
         unbounded = TRUE;
         break;
      }
   }

   return unbounded;
}

/** enforces Benders' constraints for given solution
 *
 *  This method is called from cons_benderslp and cons_benders. If the method is called from cons_benderslp, then the
 *  solutions are not guaranteed to be integer feasible. This is because the default priority is set greater than the
 *  integer constraint handler. If this method is called from cons_benders, then, because the default enforcement
 *  priority is set less than that of the integer constraint handler, then it can be assumed that the solutions are
 *  integer feasible.
 *
 *  The checkint flag indicates whether integer feasibility can be assumed. If it is not assumed, i.e. checkint ==
 *  FALSE, then only the convex relaxations of the subproblems are solved. If integer feasibility is assumed, i.e.
 *  checkint == TRUE, then the convex relaxations and the full CIP are solved to generate Benders' cuts and check
 *  solution feasibility.
 */
SCIP_RETCODE SCIPconsBendersEnforceSolution(
   SCIP*                 scip,               /**< the SCIP instance */
   SCIP_SOL*             sol,                /**< the primal solution to enforce, or NULL for the current LP/pseudo sol */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_RESULT*          result,             /**< the result of the enforcement */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should integrality be considered when checking the subproblems */
   )
{
   SCIP_BENDERS** benders;
   SCIP_Bool infeasible;
   SCIP_Bool auxviol;
   int nactivebenders;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   (*result) = SCIP_FEASIBLE;
   infeasible = FALSE;
   auxviol = FALSE;

   benders = SCIPgetBenders(scip);
   nactivebenders = SCIPgetNActiveBenders(scip);

   for( i = 0; i < nactivebenders; i++ )
   {
      switch( type )
      {
         case SCIP_BENDERSENFOTYPE_LP:
            if( SCIPbendersCutLP(benders[i]) )
            {
               SCIP_Bool unbounded = FALSE;

               /* if the solution is unbounded, then it may not be possible to generate any Benders' decomposition
                * cuts. If the unboundedness is from the auxiliary variables, then cuts are required. Otherwise, if
                * the unboundedness comes from original variables, then the unboundedness needs to be handled by other
                * constraint handlers or the problem is reported as unbounded
                * */
               if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
               {
                  if( !unboundedAuxiliaryVariables(scip, benders[i], NULL) )
                  {
                     (*result) = SCIP_FEASIBLE;
                     auxviol = FALSE;
                     unbounded = TRUE;
                  }
               }

               if( !unbounded )
               {
                  SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], NULL, result, &infeasible, &auxviol, type, checkint) );
               }
            }
            break;
         case SCIP_BENDERSENFOTYPE_RELAX:
            if( SCIPbendersCutRelaxation(benders[i]) )
            {
               SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], sol, result, &infeasible, &auxviol, type, checkint) );
            }
            break;
         case SCIP_BENDERSENFOTYPE_PSEUDO:
            if( SCIPbendersCutPseudo(benders[i]) )
            {
               SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], NULL, result, &infeasible, &auxviol, type, checkint) );
            }
            break;
         case SCIP_BENDERSENFOTYPE_CHECK:
            SCIPwarningMessage(scip, "The conscheck callback is not supported\n");
            break;
         default:
            break;
      }  /*lint !e788*/

      /* The decompositions are checked until one is found that is not feasible. Not being feasible could mean that
       * infeasibility of the original problem has been proven or a constraint has been added. If the result
       * SCIP_DIDNOTRUN is returned, then the next decomposition is checked */
      if( (*result) != SCIP_FEASIBLE && (*result) != SCIP_DIDNOTRUN )
         break;
   }

   /* if the constraint handler was called with an integer feasible solution, then a feasible solution can be proposed */
   if( checkint )
   {
      /* in the case that the problem is feasible, this means that all subproblems are feasible. The auxiliary variables
       * still need to be updated. This is done by constructing a valid solution. */
      if( (*result) == SCIP_FEASIBLE && auxviol )
      {
         SCIP_CALL( constructValidSolution(scip, conshdlr, sol, type) );

         (*result) = SCIP_INFEASIBLE;
      }
   }

   /* if no Benders' decomposition were run, then the result is returned as SCIP_FEASIBLE. The SCIP_DIDNOTRUN result
    * indicates that no subproblems were checked or that cuts were disabled, so that it is not guaranteed that this
    * solution is feasible.
    */
   if( (*result) == SCIP_DIDNOTRUN )
      (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBenders)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeConshdlrBenders(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* freeing the constraint handler data */
   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   conshdlrdata->checkedsolssize = DEFAULT_CHECKEDSOLSSIZE;
   conshdlrdata->ncheckedsols = 0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->checkedsols, conshdlrdata->checkedsolssize) );

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* freeing the checked sols array */
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->checkedsols, conshdlrdata->checkedsolssize);

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->active )
   {
      SCIP_CALL( SCIPconsBendersEnforceSolution(scip, NULL, conshdlr, result, SCIP_BENDERSENFOTYPE_LP, TRUE) );
   }
   else
      (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->active )
   {
      SCIP_CALL( SCIPconsBendersEnforceSolution(scip, sol, conshdlr, result, SCIP_BENDERSENFOTYPE_RELAX, TRUE) );
   }
   else
      (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->active )
   {
      SCIP_CALL( SCIPconsBendersEnforceSolution(scip, NULL, conshdlr, result, SCIP_BENDERSENFOTYPE_PSEUDO, TRUE) );
   }
   else
      (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions
 *
 *  This function checks the feasibility of the Benders' decomposition master problem. In the case that the problem is
 *  feasible, then the auxiliary variables must be updated with the subproblem objective function values. It is not
 *  possible to simply update the auxiliary variable values, so a new solution is created.
 */
static
SCIP_DECL_CONSCHECK(consCheckBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_BENDERS** benders;
   int nactivebenders;
   int solindex;
   int i;
   SCIP_Bool performcheck;
   SCIP_Bool infeasible;
   SCIP_Bool auxviol;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   (*result) = SCIP_FEASIBLE;
   performcheck = TRUE;
   infeasible = FALSE;
   auxviol = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* if the constraint handler is active, then the check must be performed.  */
   if( conshdlrdata->active )
   {
      benders = SCIPgetBenders(scip);
      nactivebenders = SCIPgetNActiveBenders(scip);

      /* checking if the solution was constructed by this constraint handler */
      solindex = SCIPsolGetIndex(sol);
      for( i = 0; i < conshdlrdata->ncheckedsols; i++ )
      {
         if( conshdlrdata->checkedsols[i] == solindex )
         {
            conshdlrdata->checkedsols[0] = conshdlrdata->checkedsols[conshdlrdata->ncheckedsols - 1];
            conshdlrdata->ncheckedsols--;

            performcheck = FALSE;
            break;
         }
      }

      /* if the solution has not been checked before, then we must perform the check */
      if( performcheck && nactivebenders > 0 )
      {
         for( i = 0; i < nactivebenders; i++ )
         {
            SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], sol, result, &infeasible, &auxviol,
                  SCIP_BENDERSENFOTYPE_CHECK, TRUE) );

            /* in the case of multiple Benders' decompositions, the subproblems are solved until a constriant is added or
             * infeasibility is proven. So if the result is not SCIP_FEASIBLE, then the loop is exited */
            if( (*result) != SCIP_FEASIBLE )
               break;
         }

         /* in the case that the problem is feasible, this means that all subproblems are feasible. The auxiliary variables
          * still need to be updated. This is done by constructing a valid solution. */
         if( (*result) == SCIP_FEASIBLE )
         {
            if( auxviol )
            {
               if( !SCIPsolIsOriginal(sol) )
               {
                  SCIP_CALL( constructValidSolution(scip, conshdlr, sol, SCIP_BENDERSENFOTYPE_CHECK) );
               }

               if( printreason )
                  SCIPmessagePrintInfo(SCIPgetMessagehdlr(scip), "all subproblems are feasible but there is a violation in the auxiliary variables\n");

               (*result) = SCIP_INFEASIBLE;
            }
         }

         /* if no Benders' decomposition were run, then the result is returned as SCIP_FEASIBLE. The SCIP_DIDNOTRUN result
          * indicates that no subproblems were checked or that cuts were disabled, so that it is not guaranteed that this
          * solution is feasible.
          */
         if( (*result) == SCIP_DIDNOTRUN )
            (*result) = SCIP_FEASIBLE;
      }
   }

   return SCIP_OKAY;
}


/** the presolving method for the Benders' decomposition constraint handler
 *
 *  This method is used to update the lower bounds of the auxiliary problem and to identify infeasibility before the
 *  subproblems are solved. When SCIP is copied, the Benders' decomposition subproblems from the source SCIP are
 *  transferred to the target SCIP. So there is no need to perform this presolving step in the copied SCIP, since the
 *  computed bounds would be identical.
 */
static
SCIP_DECL_CONSPRESOL(consPresolBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_BENDERS** benders;
   int nactivebenders;
   int nsubproblems;
   int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   (*result) = SCIP_DIDNOTFIND;

   /* this presolving step is only valid for the main SCIP instance. If the SCIP instance is copied, then the presolving
    * step is not performed.
    */
   if( SCIPgetSubscipDepth(scip) > 0 )
   {
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* it is only possible to compute the lower bound of the subproblems if the constraint handler is active */
   if( conshdlrdata->active )
   {
      benders = SCIPgetBenders(scip);
      nactivebenders = SCIPgetNActiveBenders(scip);

      /* need to compute the lower bound for all active Benders' decompositions */
      for( i = 0; i < nactivebenders; i++ )
      {
         nsubproblems = SCIPbendersGetNSubproblems(benders[i]);

         for( j = 0; j < nsubproblems; j++ )
         {
            SCIP_VAR* auxiliaryvar;
            SCIP_Real lowerbound;
            SCIP_Bool infeasible;

            infeasible = FALSE;

            /* computing the lower bound of the subproblem by solving it without any variable fixings */
            SCIP_CALL( SCIPcomputeBendersSubproblemLowerbound(scip, benders[i], j, &lowerbound, &infeasible) );

            if( infeasible )
            {
               (*result) = SCIP_CUTOFF;
               break;
            }

            /* retrieving the auxiliary variable */
            auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders[i], j);

            /* only update the lower bound if it is greater than the current lower bound */
            if( SCIPisGT(scip, lowerbound, SCIPvarGetLbLocal(auxiliaryvar)) )
            {
               SCIPdebugMsg(scip, "Tightened lower bound of <%s> to %g\n", SCIPvarGetName(auxiliaryvar), lowerbound);
               /* updating the lower bound of the auxiliary variable */
               SCIP_CALL( SCIPchgVarLb(scip, auxiliaryvar, lowerbound) );

               (*nchgbds)++;
               (*result) = SCIP_SUCCESS;
            }

            /* stores the lower bound for the subproblem */
            SCIPbendersUpdateSubproblemLowerbound(benders[i], j, lowerbound);
         }

         if( (*result) == SCIP_CUTOFF )
            break;
      }
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler
 *  The auxiliary variables and the master problem variables need to lock added by the Benders' decomposition
 *  constraint. The auxiliary variables require a down lock. The master problem variable need both up and down lock.
 *  The master problem variables require locks in both directions because the coefficients in all potential Benders'
 *  cuts are not known in general.
 */
static
SCIP_DECL_CONSLOCK(consLockBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_BENDERS** benders;
   SCIP_VAR** vars;
   int nactivebenders;
   int nsubproblems;
   int nvars;
   int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(locktype == SCIP_LOCKTYPE_MODEL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* the locks should only be added if the Benders' decomposition constraint handler has been activated */
   if( conshdlrdata->active )
   {
      benders = SCIPgetBenders(scip);
      nactivebenders = SCIPgetNActiveBenders(scip);

      /* retrieving the master problem variables */
      SCIP_CALL( SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

      /* need to compute the lower bound for all active Benders' decompositions */
      for( i = 0; i < nactivebenders; i++ )
      {
         nsubproblems = SCIPbendersGetNSubproblems(benders[i]);

         /* if the auxiliary variable exists, then we need to add a down lock. Initially, a down lock is added to all
          * auxiliary variables during creating. This is because the creation of auxiliary variable occurs after
          * CONS_LOCK is called. The inclusion of the auxiliary variables in this function is to cover the case if locks
          * are added or removed after presolving.
          */
         for( j = 0; j < nsubproblems; j++ )
         {
            SCIP_VAR* auxvar;

            auxvar = SCIPbendersGetAuxiliaryVar(benders[i], j);

            if( auxvar != NULL )
            {
               SCIP_CALL( SCIPaddVarLocksType(scip, auxvar, locktype, nlockspos, nlocksneg) );
            }
         }

         /* adding up and down locks for all master problem variables. Since the locks for all constraint handlers
          * without constraints, no auxiliary variables have been added. As such, all variables are master variables.
          */
         for( j = 0; j < nvars; j++ )
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, vars[j], locktype, (nlockspos + nlocksneg)*nsubproblems,
                  (nlockspos + nlocksneg)*nsubproblems) );
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBenders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create benders constraint handler data */
   conshdlrdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   conshdlr = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpBenders, consEnfopsBenders, consCheckBenders, consLockBenders,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitBenders) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitBenders) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyBenders, NULL) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeBenders) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxBenders) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolBenders, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );

   /* add Benders' decomposition constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/active", "is the Benders' decomposition constraint handler active?",
         &conshdlrdata->active, FALSE, DEFAULT_ACTIVE, NULL, NULL));

   return SCIP_OKAY;
}
