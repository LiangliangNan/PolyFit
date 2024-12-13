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

/**@file   benders.c
 * @ingroup OTHER_CFILES
 * @brief  methods for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/dcmp.h"
#include "scip/paramset.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pricestore.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/benders.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"

#include "scip/struct_benders.h"
#include "scip/struct_benderscut.h"

#include "scip/benderscut.h"

/* Defaults for parameters */
#define SCIP_DEFAULT_TRANSFERCUTS         FALSE  /** should Benders' cuts generated in LNS heuristics be transferred to the main SCIP instance? */
#define SCIP_DEFAULT_CUTSASCONSS           TRUE  /** should the transferred cuts be added as constraints? */
#define SCIP_DEFAULT_LNSCHECK              TRUE  /** should the Benders' decomposition be used in LNS heuristics */
#define SCIP_DEFAULT_LNSMAXDEPTH             -1  /** maximum depth at which the LNS check is performed */
#define SCIP_DEFAULT_LNSMAXCALLS             10  /** the maximum number of Benders' decomposition calls in LNS heuristics */
#define SCIP_DEFAULT_LNSMAXCALLSROOT          0  /** the maximum number of root node Benders' decomposition calls in LNS heuristics */
#define SCIP_DEFAULT_SUBPROBFRAC            1.0  /** fraction of subproblems that are solved in each iteration */
#define SCIP_DEFAULT_UPDATEAUXVARBOUND    FALSE  /** should the auxiliary variable lower bound be updated by solving the subproblem */
#define SCIP_DEFAULT_AUXVARSIMPLINT       FALSE  /** set the auxiliary variables as implint if the subproblem objective is integer */
#define SCIP_DEFAULT_CUTCHECK              TRUE  /** should cuts be generated during the checking of solutions? */
#define SCIP_DEFAULT_STRENGTHENMULT         0.5  /** the convex combination multiplier for the cut strengthening */
#define SCIP_DEFAULT_NOIMPROVELIMIT           5  /** the maximum number of cut strengthening without improvement */
#define SCIP_DEFAULT_STRENGTHENPERTURB    1e-06  /** the amount by which the cut strengthening solution is perturbed */
#define SCIP_DEFAULT_STRENGTHENENABLED    FALSE  /** enable the core point cut strengthening approach */
#define SCIP_DEFAULT_STRENGTHENINTPOINT     'r'  /** where should the strengthening interior point be sourced from ('l'p relaxation, 'f'irst solution, 'i'ncumbent solution, 'r'elative interior point, vector of 'o'nes, vector of 'z'eros) */
#define SCIP_DEFAULT_NUMTHREADS               1  /** the number of parallel threads to use when solving the subproblems */
#define SCIP_DEFAULT_EXECFEASPHASE        FALSE  /** should a feasibility phase be executed during the root node processing */
#define SCIP_DEFAULT_SLACKVARCOEF          1e+6  /** the initial objective coefficient of the slack variables in the subproblem */
#define SCIP_DEFAULT_MAXSLACKVARCOEF       1e+9  /** the maximal objective coefficient of the slack variables in the subproblem */
#define SCIP_DEFAULT_CHECKCONSCONVEXITY    TRUE  /** should the constraints of the subproblem be checked for convexity? */

#define BENDERS_MAXPSEUDOSOLS                 5  /** the maximum number of pseudo solutions checked before suggesting
                                                  *  merge candidates */

#define BENDERS_ARRAYSIZE        1000    /**< the initial size of the added constraints/cuts arrays */

#define AUXILIARYVAR_NAME     "##bendersauxiliaryvar" /** the name for the Benders' auxiliary variables in the master problem */
#define SLACKVAR_NAME         "##bendersslackvar"     /** the name for the Benders' slack variables added to each
                                                       *  constraints in the subproblems */
#define NLINEARCONSHDLRS 5

/* event handler properties */
#define NODEFOCUS_EVENTHDLR_NAME         "bendersnodefocus"
#define NODEFOCUS_EVENTHDLR_DESC         "node focus event handler for Benders' decomposition"

#define MIPNODEFOCUS_EVENTHDLR_NAME      "bendersmipsolvenodefocus"
#define MIPNODEFOCUS_EVENTHDLR_DESC      "node focus event handler for the MIP solve method for Benders' decomposition"

#define UPPERBOUND_EVENTHDLR_NAME        "bendersupperbound"
#define UPPERBOUND_EVENTHDLR_DESC        "found solution event handler to terminate subproblem solve for a given upper bound"

#define NODESOLVED_EVENTHDLR_NAME        "bendersnodesolved"
#define NODESOLVED_EVENTHDLR_DESC        "node solved event handler for the Benders' integer cuts"


/** event handler data */
struct SCIP_EventhdlrData
{
   int                   filterpos;          /**< the event filter entry */
   int                   numruns;            /**< the number of times that the problem has been solved */
   SCIP_Real             upperbound;         /**< an upper bound for the problem */
   SCIP_Bool             solvecip;           /**< is the event called from a MIP subproblem solve*/
};


/* ---------------- Local methods for event handlers ---------------- */

/** initialises the members of the eventhandler data */
static
SCIP_RETCODE initEventhandlerData(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< the event handler data */
   )
{
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->filterpos = -1;
   eventhdlrdata->numruns = 0;
   eventhdlrdata->upperbound = -SCIPinfinity(scip);
   eventhdlrdata->solvecip = FALSE;

   return SCIP_OKAY;
}

/** initsol method for the event handlers */
static
SCIP_RETCODE initsolEventhandler(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< the event handlers data structure */
   SCIP_EVENTTYPE        eventtype           /**< event type mask to select events to catch */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   SCIP_CALL(SCIPcatchEvent(scip, eventtype, eventhdlr, NULL, &eventhdlrdata->filterpos));

   return SCIP_OKAY;
}

/** the exit sol method for the event handlers */
static
SCIP_RETCODE exitsolEventhandler(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< the event handlers data structure */
   SCIP_EVENTTYPE        eventtype           /**< event type mask to select events to catch */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   if( eventhdlrdata->filterpos >= 0 )
   {
      SCIP_CALL(SCIPdropEvent(scip, eventtype, eventhdlr, NULL, eventhdlrdata->filterpos));
      eventhdlrdata->filterpos = -1;
   }

   return SCIP_OKAY;
}

/** the exit method for the event handlers */
static
SCIP_RETCODE exitEventhandler(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< the event handlers data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   /* reinitialise the event handler data */
   SCIP_CALL( initEventhandlerData(scip, eventhdlrdata) );

   return SCIP_OKAY;
}

/** free method for the event handler */
static
SCIP_RETCODE freeEventhandler(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< the event handlers data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);

   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}



/* ---------------- Callback methods of node focus event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersNodefocus)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   /* sending an interrupt solve signal to return the control back to the Benders' decomposition plugin.
    * This will ensure the SCIP stage is SCIP_STAGE_SOLVING, allowing the use of probing mode. */
   SCIP_CALL( SCIPinterruptSolve(scip) );

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, eventhdlrdata->filterpos));
   eventhdlrdata->filterpos = -1;

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( initsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_NODEFOCUSED) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_NODEFOCUSED) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTFREE(eventFreeBendersNodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( freeEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}


/* ---------------- Callback methods of MIP solve node focus event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersMipnodefocus)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   /* interrupting the solve so that the control is returned back to the Benders' core. */
   if( eventhdlrdata->numruns == 0 && !eventhdlrdata->solvecip )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, eventhdlrdata->filterpos));
   eventhdlrdata->filterpos = -1;

   eventhdlrdata->numruns++;

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( initsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_NODEFOCUSED) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_NODEFOCUSED) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTFREE(eventFreeBendersMipnodefocus)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), MIPNODEFOCUS_EVENTHDLR_NAME) == 0);

   SCIP_CALL( freeEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/* ---------------- Callback methods of solution found event handler ---------------- */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersUpperbound)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SOL* bestsol;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   bestsol = SCIPgetBestSol(scip);

   if( SCIPisLT(scip, SCIPgetSolOrigObj(scip, bestsol)*(int)SCIPgetObjsense(scip), eventhdlrdata->upperbound) )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL( initsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_BESTSOLFOUND) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitsolEventhandler(scip, eventhdlr, SCIP_EVENTTYPE_BESTSOLFOUND) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL( exitEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTFREE(eventFreeBendersUpperbound)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), UPPERBOUND_EVENTHDLR_NAME) == 0);

   SCIP_CALL( freeEventhandler(scip, eventhdlr) );

   return SCIP_OKAY;
}

/** updates the upper bound in the event handler data */
static
SCIP_RETCODE updateEventhdlrUpperbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             upperbound          /**< the upper bound value */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   eventhdlr = SCIPfindEventhdlr(SCIPbendersSubproblem(benders, probnumber), UPPERBOUND_EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->upperbound = upperbound;

   return SCIP_OKAY;
}

/* ---------------- Callback methods of the node solved event handler ---------------- */

/** Updates the cut constant of the Benders' cuts data.
 *  This function solves the master problem with only the auxiliary variables in the objective function.
 */
static
SCIP_RETCODE updateSubproblemLowerbound(
   SCIP*                 masterprob,         /**< the SCIP instance of the master problem */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nsubproblems;
   int i;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   assert(masterprob != NULL);
   assert(benders != NULL);

   /* don't run in probing or in repropagation */
   if( SCIPinProbing(masterprob) || SCIPinRepropagation(masterprob) || SCIPinDive(masterprob) )
      return SCIP_OKAY;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   SCIP_CALL( SCIPstartProbing(masterprob) );

   /* change the master problem variables to 0 */
   nvars = SCIPgetNVars(masterprob);
   vars = SCIPgetVars(masterprob);

   /* setting the objective function coefficient to 0 for all variables */
   for( i = 0; i < nvars; i++ )
   {
      if( SCIPvarGetStatus(vars[i]) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPchgVarObjProbing(masterprob, vars[i], 0.0) );
      }
   }

   /* solving an LP for all subproblems to find the lower bound */
   for( i = 0; i < nsubproblems; i++)
   {
      SCIP_VAR* auxiliaryvar;

      auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, i);

      if( SCIPvarGetStatus(auxiliaryvar) != SCIP_VARSTATUS_COLUMN )
         continue;

      SCIP_CALL( SCIPchgVarObjProbing(masterprob, auxiliaryvar, 1.0) );

      /* solving the probing LP to get a lower bound on the auxiliary variables */
      SCIP_CALL( SCIPsolveProbingLP(masterprob, -1, &lperror, &cutoff) );

      if( !SCIPisInfinity(masterprob, -SCIPgetSolTransObj(masterprob, NULL)) )
         SCIPbendersUpdateSubproblemLowerbound(benders, i, SCIPgetSolTransObj(masterprob, NULL));

      SCIPdebugMsg(masterprob, "Cut constant for subproblem %d: %g\n", i,
         SCIPbendersGetSubproblemLowerbound(benders, i));

      SCIP_CALL( SCIPchgVarObjProbing(masterprob, auxiliaryvar, 0.0) );
   }

   SCIP_CALL( SCIPendProbing(masterprob) );

   return SCIP_OKAY;
}

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBendersNodesolved)
{  /*lint --e{715}*/
   SCIP_BENDERS* benders;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODESOLVED_EVENTHDLR_NAME) == 0);

   benders = (SCIP_BENDERS*)SCIPeventhdlrGetData(eventhdlr);   /*lint !e826*/

   if( SCIPbendersGetNSubproblems(benders) > 0
      && SCIPbendersGetNSubproblems(benders) > SCIPbendersGetNConvexSubproblems(benders) )
   {
      SCIP_CALL( updateSubproblemLowerbound(scip, benders) );
   }

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolBendersNodesolved)
{
   SCIP_BENDERS* benders;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), NODESOLVED_EVENTHDLR_NAME) == 0);

   /* getting the Benders' decomposition data structure */
   benders = (SCIP_BENDERS*)SCIPeventhdlrGetData(eventhdlr);   /*lint !e826*/

   /* The event is only caught if there is an active Benders' decomposition, the integer subproblem are solved and
    * the Benders' decomposition has not been copied in thread safe mode
    */
   if( SCIPbendersIsActive(benders) && !SCIPbendersOnlyCheckConvexRelax(benders, SCIPgetSubscipsOff(scip))
      && !benders->threadsafe )
   {
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, NULL) );
   }

   return SCIP_OKAY;
}


/* ---------------- Methods for the parallelisation of Benders' decomposition ---------------- */

/** comparison method for sorting the subproblems.
 *  The subproblem that has been called the least is prioritised
 */
static
SCIP_DECL_SORTPTRCOMP(benderssubcompdefault)
{
   SCIP_SUBPROBLEMSOLVESTAT* solvestat1;
   SCIP_SUBPROBLEMSOLVESTAT* solvestat2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   solvestat1 = (SCIP_SUBPROBLEMSOLVESTAT*)elem1;
   solvestat2 = (SCIP_SUBPROBLEMSOLVESTAT*)elem2;

   /* prefer subproblems with fewer calls, using the index as tie breaker */
   if( MAX(solvestat1->ncalls, solvestat2->ncalls) == 0 )
      return solvestat1->idx - solvestat2->idx;
   else if( solvestat1->ncalls != solvestat2->ncalls )
      return solvestat1->ncalls - solvestat2->ncalls;
   else
   {
      /* prefer the harder problem (with more average iterations) */
      int avgiterdiff = (int)solvestat2->avgiter - (int)solvestat1->avgiter;

      if( avgiterdiff != 0 )
         return avgiterdiff;

      return solvestat1->idx - solvestat2->idx;
   }

/* the code below does not give a total order of the elements */
#ifdef SCIP_DISABLED_CODE
   if( solvestat1->ncalls == 0 )
      if( solvestat2->ncalls == 0 )
         if( solvestat1->idx < solvestat2->idx )
            return -1;
         else
            return 1;
      else
         return -1;
   else if( solvestat2->ncalls == 0 )
      return 1;
   else
   {
      if( solvestat1->ncalls < solvestat2->ncalls )
         return -1;
      else if( solvestat2->ncalls < solvestat1->ncalls )
         return 1;
      else
      {
         /* we want to execute the hard subproblems first */
         if( solvestat1->avgiter > solvestat2->avgiter )
            return 1;
         else
            return -1;
      }
   }
#endif
}

/* Local methods */

/** A workaround for GCG. This is a temp vardata that is set for the auxiliary variables */
struct SCIP_VarData
{
   int                   vartype;             /**< the variable type. In GCG this indicates whether the variable is a
                                               *   master problem or subproblem variable. */
};

/** adds the auxiliary variables to the Benders' decomposition master problem */
static
SCIP_RETCODE addAuxiliaryVariablesToMaster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition structure */
   )
{
   SCIP_BENDERS* topbenders;        /* the highest priority Benders' decomposition */
   SCIP_VAR* auxiliaryvar;
   SCIP_VARDATA* vardata;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */
   SCIP_Bool shareauxvars;
   int i;

   /* this is a workaround for GCG. GCG expects that the variable has vardata when added. So a dummy vardata is created */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = -1;

   /* getting the highest priority Benders' decomposition */
   topbenders = SCIPgetBenders(scip)[0];

   /* if the current Benders is the highest priority Benders, then we need to create the auxiliary variables.
    * Otherwise, if the shareauxvars flag is set, then the auxiliary variables from the highest priority Benders' are
    * stored with this Benders. */
   shareauxvars = FALSE;
   if( topbenders != benders && SCIPbendersShareAuxVars(benders) )
      shareauxvars = TRUE;

   for( i = 0; i < SCIPbendersGetNSubproblems(benders); i++ )
   {
      /* if the auxiliary variables are shared, then a pointer to the variable is retrieved from topbenders,
       * otherwise the auxiliaryvariable is created. */
      if( shareauxvars )
      {
         auxiliaryvar = SCIPbendersGetAuxiliaryVar(topbenders, i);

         SCIP_CALL( SCIPcaptureVar(scip, auxiliaryvar) );
      }
      else
      {
         SCIP_VARTYPE vartype;

         /* set the variable type of the auxiliary variables to implied integer if the objective function of the
          * subproblem is guaranteed to be integer. This behaviour is controlled through a user parameter.
          * NOTE: It is only possible to determine if the objective function is integral if the subproblem is defined as
          * a SCIP instance, i.e. not NULL.
          */
         if( benders->auxvarsimplint && SCIPbendersSubproblem(benders, i) != NULL
            && SCIPisObjIntegral(SCIPbendersSubproblem(benders, i)) )
            vartype = SCIP_VARTYPE_IMPLINT;
         else
            vartype = SCIP_VARTYPE_CONTINUOUS;

         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_%d_%s", AUXILIARYVAR_NAME, i, SCIPbendersGetName(benders) );
         SCIP_CALL( SCIPcreateVarBasic(scip, &auxiliaryvar, varname, benders->subproblowerbound[i], SCIPinfinity(scip),
               1.0, vartype) );

         SCIPvarSetData(auxiliaryvar, vardata);

         SCIP_CALL( SCIPaddVar(scip, auxiliaryvar) );

         /* adding the down lock for the Benders' decomposition constraint handler */
         SCIP_CALL( SCIPaddVarLocksType(scip, auxiliaryvar, SCIP_LOCKTYPE_MODEL, 1, 0) );
      }

      benders->auxiliaryvars[i] = auxiliaryvar;
   }

   SCIPfreeBlockMemory(scip, &vardata);

   return SCIP_OKAY;
}

/** assigns the copied auxiliary variables in the target SCIP to the target Benders' decomposition data */
static
SCIP_RETCODE assignAuxiliaryVariables(
   SCIP*                 scip,               /**< SCIP data structure, the target scip */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_BENDERS* topbenders;        /* the highest priority Benders' decomposition */
   SCIP_VAR* targetvar;
   SCIP_VARDATA* vardata;
   char varname[SCIP_MAXSTRLEN];    /* the name of the auxiliary variable */
   SCIP_Bool shareauxvars;
   int subscipdepth;
   int i;
   int j;

   assert(scip != NULL);
   assert(benders != NULL);

   /* this is a workaround for GCG. GCG expects that the variable has vardata when added. So a dummy vardata is created */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vardata) );
   vardata->vartype = -1;

   /* getting the highest priority Benders' decomposition */
   topbenders = SCIPgetBenders(scip)[0];

   /* if the auxiliary variable are shared, then the variable name will have a suffix of the highest priority Benders'
    * name. So the shareauxvars flag indicates how to search for the auxiliary variables */
   shareauxvars = FALSE;
   if( topbenders != benders && SCIPbendersShareAuxVars(benders) )
      shareauxvars = TRUE;

   subscipdepth = SCIPgetSubscipDepth(scip);

   for( i = 0; i < SCIPbendersGetNSubproblems(benders); i++ )
   {
      char prefix[SCIP_MAXSTRLEN];
      char tmpprefix[SCIP_MAXSTRLEN];
      int len = 1;

      j = 0;
      targetvar = NULL;

      /* the prefix for the variable names is required for UG, since we don't know how many copies have been made. To
       * find the target variable, we start with an empty prefix. Then t_ is prepended until the target variable is
       * found
       */
      prefix[0] = '\0';
      while( targetvar == NULL && j <= subscipdepth )
      {
         if( shareauxvars )
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s%s_%d_%s", prefix, AUXILIARYVAR_NAME, i, SCIPbendersGetName(topbenders));
         else
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s%s_%d_%s", prefix, AUXILIARYVAR_NAME, i, SCIPbendersGetName(benders));

         /* finding the variable in the copied problem that has the same name as the auxiliary variable */
         targetvar = SCIPfindVar(scip, varname);

         (void) SCIPsnprintf(tmpprefix, len, "t_%s", prefix);
         len += 2;
         (void) strncpy(prefix, tmpprefix, len); /*lint !e732*/

         j++;
      }

      if( targetvar != NULL )
      {
         SCIPvarSetData(targetvar, vardata);

         benders->auxiliaryvars[i] = SCIPvarGetTransVar(targetvar);

         SCIP_CALL( SCIPcaptureVar(scip, benders->auxiliaryvars[i]) );
      }
      else
      {
         SCIPABORT();
      }
   }

   SCIPfreeBlockMemory(scip, &vardata);

   return SCIP_OKAY;
}

/** sets the subproblem objective value array to -infinity */
static
void resetSubproblemObjectiveValue(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP* subproblem;
   SCIP_Real inf;
   int nsubproblems;
   int i;

   assert(benders != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   for( i = 0; i < nsubproblems; i++ )
   {
      subproblem = SCIPbendersSubproblem(benders, i);
      if( subproblem != NULL )
         inf = SCIPinfinity(subproblem);
      else
         inf = SCIPsetInfinity(set);

      SCIPbendersSetSubproblemObjval(benders, i, inf);
   }
}

/** compares two Benders' decompositions w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPbendersComp)
{  /*lint --e{715}*/
   return ((SCIP_BENDERS*)elem2)->priority - ((SCIP_BENDERS*)elem1)->priority;
}

/** comparison method for sorting Benders' decompositions w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPbendersCompName)
{
   return strcmp(SCIPbendersGetName((SCIP_BENDERS*)elem1), SCIPbendersGetName((SCIP_BENDERS*)elem2));
}

/** method to call, when the priority of a Benders' decomposition was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdBendersPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetBendersPriority() to mark the Benders' decompositions as unsorted */
   SCIPsetBendersPriority(scip, (SCIP_BENDERS*)paramdata, SCIPparamGetInt(param)); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a variable mapping between the master problem variables of the source scip and the sub scip */
static
SCIP_RETCODE createMasterVarMapping(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition of the target SCIP instance */
   SCIP_SET*             sourceset,          /**< global SCIP settings from the source SCIP */
   SCIP_HASHMAP*         varmap              /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables; must not be NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* targetvar;
   int nvars;
   int i;

   assert(benders != NULL);
   assert(sourceset != NULL);
   assert(benders->iscopy);
   assert(benders->mastervarsmap == NULL);

   /* getting the master problem variable data */
   vars = SCIPgetVars(sourceset->scip);
   nvars = SCIPgetNVars(sourceset->scip);

   /* creating the hashmap for the mapping between the master variable of the target and source scip */
   SCIP_CALL( SCIPhashmapCreate(&benders->mastervarsmap, SCIPblkmem(sourceset->scip), nvars) );

   for( i = 0; i < nvars; i++ )
   {
      /* getting the variable pointer for the target SCIP variables. The variable mapping returns the target SCIP
       * varibale for a given source SCIP variable. */
      targetvar = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);
      if( targetvar != NULL )
      {
         SCIP_CALL( SCIPhashmapInsert(benders->mastervarsmap, targetvar, vars[i]) );
         SCIP_CALL( SCIPcaptureVar(sourceset->scip, vars[i]) );
      }
   }

   return SCIP_OKAY;
}

/** copies the given Benders' decomposition to a new SCIP */
SCIP_RETCODE SCIPbendersCopyInclude(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             sourceset,          /**< SCIP_SET of SCIP to copy from */
   SCIP_SET*             targetset,          /**< SCIP_SET of SCIP to copy to */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables; if NULL, then the transfer of cuts is not possible */
   SCIP_Bool             threadsafe,         /**< must the Benders' decomposition copy be thread safe */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   SCIP_BENDERS* targetbenders;  /* the copy of the Benders' decomposition struct in the target set */
   int i;

   assert(benders != NULL);
   assert(targetset != NULL);
   assert(valid != NULL);
   assert(targetset->scip != NULL);

   (*valid) = FALSE;

   if( benders->benderscopy != NULL && targetset->benders_copybenders && SCIPbendersIsActive(benders) )
   {
      SCIPsetDebugMsg(targetset, "including Benders' decomposition %s in subscip %p\n", SCIPbendersGetName(benders), (void*)targetset->scip);
      SCIP_CALL( benders->benderscopy(targetset->scip, benders, threadsafe) );

      /* copying the Benders' cuts */
      targetbenders = SCIPsetFindBenders(targetset, SCIPbendersGetName(benders));

      /* storing the pointer to the source scip instance */
      targetbenders->sourcescip = sourceset->scip;

      /* the flag is set to indicate that the Benders' decomposition is a copy */
      targetbenders->iscopy = TRUE;

      /* storing whether the lnscheck should be performed */
      targetbenders->lnscheck = benders->lnscheck;
      targetbenders->lnsmaxdepth = benders->lnsmaxdepth;
      targetbenders->lnsmaxcalls = benders->lnsmaxcalls;
      targetbenders->lnsmaxcallsroot = benders->lnsmaxcallsroot;

      /* storing whether the Benders' copy required thread safety */
      targetbenders->threadsafe = threadsafe;

      /* calling the copy method for the Benders' cuts */
      SCIPbendersSortBenderscuts(benders);
      for( i = 0; i < benders->nbenderscuts; i++ )
      {
         SCIP_CALL( SCIPbenderscutCopyInclude(targetbenders, benders->benderscuts[i], targetset) );
      }

      /* When the Benders' decomposition is copied then a variable mapping between the master problem variables is
       * required. This variable mapping is used to transfer the cuts generated in the target SCIP to the source SCIP.
       * The variable map is stored in the target Benders' decomposition. This will be freed when the sub-SCIP is freed.
       */
      if( varmap != NULL )
      {
         SCIP_CALL( createMasterVarMapping(targetbenders, sourceset, varmap) );
      }

      assert((varmap != NULL && targetbenders->mastervarsmap != NULL)
         || (varmap == NULL && targetbenders->mastervarsmap == NULL));
   }

   /* if the Benders' decomposition is active, then copy is not valid. */
   (*valid) = !SCIPbendersIsActive(benders);

   return SCIP_OKAY;
}

/** internal method for creating a Benders' decomposition structure */
static
SCIP_RETCODE doBendersCreate(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_Bool             shareauxvars,       /**< should this Benders' use the highest priority Benders aux vars */
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre)),/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre)),/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve)),/**< called prior to the subproblem solving loop */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex)),/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve)),/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub)),/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(benders != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   /* Checking whether the benderssolvesub and the bendersfreesub are both implemented or both are not implemented */
   if( (benderssolvesubconvex == NULL && benderssolvesub == NULL && bendersfreesub != NULL)
      || ((benderssolvesubconvex != NULL || benderssolvesub != NULL) && bendersfreesub == NULL) )
   {
      SCIPerrorMessage("Benders' decomposition <%s> requires that if bendersFreesub%s is implemented, then at least "
         "one of bendersSolvesubconvex%s or bendersSolvesub%s are implemented.\n", name, name, name, name);
      return SCIP_INVALIDCALL;
   }

   SCIP_ALLOC( BMSallocMemory(benders) );
   BMSclearMemory(*benders);
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benders)->desc, desc, strlen(desc)+1) );
   (*benders)->priority = priority;
   (*benders)->cutlp = cutlp;
   (*benders)->cutpseudo = cutpseudo;
   (*benders)->cutrelax = cutrelax;
   (*benders)->shareauxvars = shareauxvars;
   (*benders)->benderscopy = benderscopy;
   (*benders)->bendersfree = bendersfree;
   (*benders)->bendersinit = bendersinit;
   (*benders)->bendersexit = bendersexit;
   (*benders)->bendersinitpre = bendersinitpre;
   (*benders)->bendersexitpre = bendersexitpre;
   (*benders)->bendersinitsol = bendersinitsol;
   (*benders)->bendersexitsol = bendersexitsol;
   (*benders)->bendersgetvar = bendersgetvar;
   (*benders)->benderscreatesub = benderscreatesub;
   (*benders)->benderspresubsolve = benderspresubsolve;
   (*benders)->benderssolvesubconvex = benderssolvesubconvex;
   (*benders)->benderssolvesub = benderssolvesub;
   (*benders)->benderspostsolve = benderspostsolve;
   (*benders)->bendersfreesub = bendersfreesub;
   (*benders)->bendersdata = bendersdata;
   SCIP_CALL( SCIPclockCreate(&(*benders)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*benders)->bendersclock, SCIP_CLOCKTYPE_DEFAULT) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of Benders' decomposition <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*benders)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
         paramChgdBendersPriority, (SCIP_PARAMDATA*)(*benders)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutlp", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts be generated for LP solutions?", &(*benders)->cutlp, FALSE, cutlp, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutpseudo", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts be generated for pseudo solutions?", &(*benders)->cutpseudo, FALSE, cutpseudo, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutrelax", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts be generated for relaxation solutions?", &(*benders)->cutrelax, FALSE, cutrelax, NULL, NULL) ); /*lint !e740*/

   /* These parameters are left for the user to decide in a settings file. This departs from the usual SCIP convention
    * where the settings available at the creation of the plugin can be set in the function call.
    */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/transfercuts", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts from LNS heuristics be transferred to the main SCIP instance?", &(*benders)->transfercuts,
         FALSE, SCIP_DEFAULT_TRANSFERCUTS, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnscheck", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' decomposition be used in LNS heurisics?", &(*benders)->lnscheck, FALSE, SCIP_DEFAULT_LNSCHECK,
         NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnsmaxdepth", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "maximum depth at which the LNS check is performed (-1: no limit)", &(*benders)->lnsmaxdepth, TRUE,
         SCIP_DEFAULT_LNSMAXDEPTH, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnsmaxcalls", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "the maximum number of Benders' decomposition calls in LNS heuristics (-1: no limit)", &(*benders)->lnsmaxcalls,
         TRUE, SCIP_DEFAULT_LNSMAXCALLS, -1, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/lnsmaxcallsroot", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "the maximum number of root node Benders' decomposition calls in LNS heuristics (-1: no limit)",
         &(*benders)->lnsmaxcallsroot, TRUE, SCIP_DEFAULT_LNSMAXCALLSROOT, -1, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutsasconss", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should the transferred cuts be added as constraints?", &(*benders)->cutsasconss, FALSE,
         SCIP_DEFAULT_CUTSASCONSS, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/subprobfrac", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
         "fraction of subproblems that are solved in each iteration", &(*benders)->subprobfrac, FALSE,
         SCIP_DEFAULT_SUBPROBFRAC, 0.0, 1.0, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/updateauxvarbound", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should the auxiliary variable bound be updated by solving the subproblem?", &(*benders)->updateauxvarbound,
         FALSE, SCIP_DEFAULT_UPDATEAUXVARBOUND, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/auxvarsimplint", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "if the subproblem objective is integer, then define the auxiliary variables as implied integers?",
         &(*benders)->auxvarsimplint, FALSE, SCIP_DEFAULT_AUXVARSIMPLINT, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutcheck", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should Benders' cuts be generated while checking solutions?",
         &(*benders)->cutcheck, FALSE, SCIP_DEFAULT_CUTCHECK, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutstrengthenmult", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
         "the convex combination multiplier for the cut strengthening", &(*benders)->convexmult, FALSE,
         SCIP_DEFAULT_STRENGTHENMULT, 0.0, 1.0, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/noimprovelimit", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "the maximum number of cut strengthening without improvement", &(*benders)->noimprovelimit, TRUE,
         SCIP_DEFAULT_NOIMPROVELIMIT, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/corepointperturb", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
         "the constant use to perturb the cut strengthening core point", &(*benders)->perturbeps, FALSE,
         SCIP_DEFAULT_STRENGTHENPERTURB, 0.0, 1.0, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutstrengthenenabled", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should the core point cut strengthening be employed (only applied to fractional solutions or continuous subproblems)?",
         &(*benders)->strengthenenabled, FALSE, SCIP_DEFAULT_STRENGTHENENABLED, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/cutstrengthenintpoint", name);
   SCIP_CALL( SCIPsetAddCharParam(set, messagehdlr, blkmem, paramname,
         "where should the strengthening interior point be sourced from ('l'p relaxation, 'f'irst solution, 'i'ncumbent solution, 'r'elative interior point, vector of 'o'nes, vector of 'z'eros)",
         &(*benders)->strengthenintpoint, FALSE, SCIP_DEFAULT_STRENGTHENINTPOINT, "lfiroz", NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/numthreads", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "the number of threads to use when solving the subproblems", &(*benders)->numthreads, TRUE,
         SCIP_DEFAULT_NUMTHREADS, 1, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/execfeasphase", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should a feasibility phase be executed during the root node, i.e. adding slack variables to constraints to ensure feasibility",
         &(*benders)->execfeasphase, FALSE, SCIP_DEFAULT_EXECFEASPHASE, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/slackvarcoef", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
         "the initial objective coefficient of the slack variables in the subproblem", &(*benders)->slackvarcoef, FALSE,
         SCIP_DEFAULT_SLACKVARCOEF, 0.0, SCIPsetInfinity(set), NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/maxslackvarcoef", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname,
         "the maximal objective coefficient of the slack variables in the subproblem", &(*benders)->maxslackvarcoef, FALSE,
         SCIP_DEFAULT_MAXSLACKVARCOEF, 0.0, SCIPsetInfinity(set), NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/checkconsconvexity", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should the constraints of the subproblems be checked for convexity?", &(*benders)->checkconsconvexity, FALSE,
         SCIP_DEFAULT_CHECKCONSCONVEXITY, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a Benders' decomposition structure
 *
 *  To use the Benders' decomposition for solving a problem, it first has to be activated with a call to SCIPactivateBenders().
 */
SCIP_RETCODE SCIPbendersCreate(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of Benders' decomposition */
   const char*           desc,               /**< description of Benders' decomposition */
   int                   priority,           /**< priority of the Benders' decomposition */
   SCIP_Bool             cutlp,              /**< should Benders' cuts be generated for LP solutions */
   SCIP_Bool             cutpseudo,          /**< should Benders' cuts be generated for pseudo solutions */
   SCIP_Bool             cutrelax,           /**< should Benders' cuts be generated for relaxation solutions */
   SCIP_Bool             shareauxvars,       /**< should this Benders' use the highest priority Benders aux vars */
   SCIP_DECL_BENDERSCOPY ((*benderscopy)),   /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree)),   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit)),   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit)),   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre)),/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre)),/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol)),/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol)),/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)),/**< returns the master variable for a given subproblem variable */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub)),/**< creates a Benders' decomposition subproblem */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve)),/**< called prior to the subproblem solving loop */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex)),/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub)),/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve)),/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub)),/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_BENDERSDATA*     bendersdata         /**< Benders' decomposition data */
   )
{
   assert(benders != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_CALL_FINALLY( doBendersCreate(benders, set, messagehdlr, blkmem, name, desc, priority, cutlp, cutpseudo,
         cutrelax, shareauxvars, benderscopy, bendersfree, bendersinit, bendersexit, bendersinitpre, bendersexitpre,
         bendersinitsol, bendersexitsol, bendersgetvar, benderscreatesub, benderspresubsolve, benderssolvesubconvex,
         benderssolvesub, benderspostsolve, bendersfreesub, bendersdata), (void) SCIPbendersFree(benders, set) );

   return SCIP_OKAY;
}


/** releases the variables that have been captured in the hashmap */
static
SCIP_RETCODE releaseVarMappingHashmapVars(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   int nentries;
   int i;

   assert(scip != NULL);
   assert(benders != NULL);

   assert(benders->mastervarsmap != NULL);

   nentries = SCIPhashmapGetNEntries(benders->mastervarsmap);

   for( i = 0; i < nentries; ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      entry = SCIPhashmapGetEntry(benders->mastervarsmap, i);

      if( entry != NULL )
      {
         SCIP_VAR* var;
         var = (SCIP_VAR*) SCIPhashmapEntryGetImage(entry);

         SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }
   }

   return SCIP_OKAY;
}


/** calls destructor and frees memory of Benders' decomposition */
SCIP_RETCODE SCIPbendersFree(
   SCIP_BENDERS**        benders,            /**< pointer to Benders' decomposition data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(*benders != NULL);
   assert(!(*benders)->initialized);
   assert(set != NULL);

   /* call destructor of Benders' decomposition */
   if( (*benders)->bendersfree != NULL )
   {
      SCIP_CALL( (*benders)->bendersfree(set->scip, *benders) );
   }

   /* if the Benders' decomposition is a copy and a varmap has been passed to SCIP_BENDERS, then the variable map
    * between the source and the target SCIP needs to be freed.
    */
   if( (*benders)->iscopy && (*benders)->mastervarsmap != NULL )
   {
      SCIP_CALL( releaseVarMappingHashmapVars((*benders)->sourcescip, (*benders)) );
      SCIPhashmapFree(&(*benders)->mastervarsmap);
   }

   /* freeing the Benders' cuts */
   for( i = 0; i < (*benders)->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutFree(&((*benders)->benderscuts[i]), set) );
   }
   BMSfreeMemoryArrayNull(&(*benders)->benderscuts);

   SCIPclockFree(&(*benders)->bendersclock);
   SCIPclockFree(&(*benders)->setuptime);
   BMSfreeMemoryArray(&(*benders)->name);
   BMSfreeMemoryArray(&(*benders)->desc);
   BMSfreeMemory(benders);

   return SCIP_OKAY;
}

/* adds a slack variable to the given constraint */
static
SCIP_RETCODE addSlackVars(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_CONS*            cons,               /**< constraint to which the slack variable(s) is added to */
   SCIP_CONSHDLR**       linearconshdlrs,    /**< an array storing the linear constraint handlers */
   SCIP_CONSHDLR*        nlconshdlr,         /**< pointer to the nonlinear constraint handler */
   int                   nlinearconshdlrs    /**< the number of linear constraint handlers */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* var;
   SCIP_Real rhs;
   SCIP_Real lhs;
   SCIP_Real objcoef;
   int i;
   SCIP_Bool linearcons;
   SCIP_Bool success;
   char name[SCIP_MAXSTRLEN];

   conshdlr = SCIPconsGetHdlr(cons);

   /* assume that the constraint is not linear, then we check whether it is linear */
   linearcons = FALSE;

   /* checking whether the constraint is a linear constraint. If so, we add a coefficient to the constraint */
   for( i = 0; i < nlinearconshdlrs; ++i )
   {
      if( conshdlr == linearconshdlrs[i] )
      {
         linearcons = TRUE;
         break;
      }
   }

   if( !linearcons && conshdlr != nlconshdlr )
   {
      SCIPwarningMessage(scip, "The subproblem includes constraint <%s>. "
         "This is not supported and the slack variable will not be added to the constraint. Feasibility cuts may be invalid.\n",
         SCIPconshdlrGetName(conshdlr));
   }

   if( linearcons )
   {
      rhs = SCIPconsGetRhs(scip, cons, &success);
      assert(success);
      lhs = SCIPconsGetLhs(scip, cons, &success);
      assert(success);
   }
   else
   {
      rhs = SCIPgetRhsNonlinear(cons);
      lhs = SCIPgetLhsNonlinear(cons);
   }

   /* getting the objective coefficient for the slack variables */
   objcoef = benders->slackvarcoef;

   /* if the right hand side is finite, then we need to add a slack variable with a negative coefficient */
   if( !SCIPisInfinity(scip, rhs) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s_neg", SLACKVAR_NAME, SCIPconsGetName(cons) );

      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, SCIPinfinity(scip), objcoef, SCIP_VARTYPE_CONTINUOUS) );

      /* adding the slack variable to the subproblem */
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* adds the slack variable to the constraint */
      if( linearcons )
      {
         SCIP_CALL( SCIPconsAddCoef(scip, cons, var, -1.0) );
      }
      else
      {
         SCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, var, -1.0) );
      }

      /* releasing the variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* if the left hand side if finite, then we need to add a slack variable with a positive coefficient */
   if( !SCIPisInfinity(scip, -lhs) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s_pos", SLACKVAR_NAME, SCIPconsGetName(cons) );

      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, SCIPinfinity(scip), objcoef, SCIP_VARTYPE_CONTINUOUS) );

      /* adding the slack variable to the subproblem */
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* adds the slack variable to the constraint */
      if( linearcons )
      {
         SCIP_CALL( SCIPconsAddCoef(scip, cons, var, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPaddLinearVarNonlinear(scip, cons, var, 1.0) );
      }

      /* releasing the variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   return SCIP_OKAY;
}

/** adds the slack variables to each of the constraints for the generation of feasibility cuts for the given non-linear
 * subproblem
 */
static
SCIP_RETCODE addSlackVarsToConstraints(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_CONSHDLR* linearconshdlrs[NLINEARCONSHDLRS];
   SCIP_CONSHDLR* nlconshdlr;
   SCIP_CONS* cons;
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* get pointers to linear constraints handlers, so can avoid string comparisons */
   linearconshdlrs[0] = SCIPfindConshdlr(subproblem, "knapsack");
   linearconshdlrs[1] = SCIPfindConshdlr(subproblem, "linear");
   linearconshdlrs[2] = SCIPfindConshdlr(subproblem, "logicor");
   linearconshdlrs[3] = SCIPfindConshdlr(subproblem, "setppc");
   linearconshdlrs[4] = SCIPfindConshdlr(subproblem, "varbound");

   nlconshdlr = SCIPfindConshdlr(subproblem, "nonlinear");

   for( i = 0; i < SCIPgetNOrigConss(subproblem); ++i )
   {
      cons = SCIPgetOrigConss(subproblem)[i];

      /* adding the slack variables to the constraint */
      SCIP_CALL( addSlackVars(subproblem, benders, cons, linearconshdlrs, nlconshdlr, NLINEARCONSHDLRS) );
   }

   return SCIP_OKAY;
}

/** initialises a MIP subproblem by putting the problem into SCIP_STAGE_SOLVING. This is achieved by calling SCIPsolve
 *  and then interrupting the solve in a node focus event handler.
 *  The LP subproblem is also initialised using this method; however, a different event handler is added. This event
 *  handler will put the LP subproblem into probing mode.
 *  The MIP solving function is called to initialise the subproblem because this function calls SCIPsolve with the
 *  appropriate parameter settings for Benders' decomposition.
 */
static
SCIP_RETCODE initialiseSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            success             /**< was the initialisation process successful */
   )
{
   SCIP* subproblem;
   SCIP_STATUS solvestatus;
   SCIP_Bool cutoff;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));
   assert(success != NULL);

   (*success) = FALSE;

   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* Getting the problem into the right SCIP stage for solving */
   SCIP_CALL( SCIPbendersSolveSubproblemCIP(set->scip, benders, probnumber, &solvestatus, FALSE) );

   /* Constructing the LP that can be solved in later iterations */
   if( solvestatus != SCIP_STATUS_BESTSOLLIMIT && solvestatus != SCIP_STATUS_TIMELIMIT
      && solvestatus != SCIP_STATUS_MEMLIMIT )
   {
      assert(SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING);

      SCIP_CALL( SCIPconstructLP(subproblem, &cutoff) );
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** initialises an LP subproblem by putting the problem into probing mode. The probing mode is invoked in a node focus
 *  event handler. This event handler is added just prior to calling the initialise subproblem function.
 */
static
SCIP_RETCODE initialiseLPSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_Bool success;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* include event handler into SCIP */
   SCIP_CALL( SCIPallocBlockMemory(subproblem, &eventhdlrdata) );

   SCIP_CALL( initEventhandlerData(subproblem, eventhdlrdata) );

   SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, NODEFOCUS_EVENTHDLR_NAME, NODEFOCUS_EVENTHDLR_DESC,
         eventExecBendersNodefocus, eventhdlrdata) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersNodefocus) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersNodefocus) );
   SCIP_CALL( SCIPsetEventhdlrExit(subproblem, eventhdlr, eventExitBendersNodefocus) );
   SCIP_CALL( SCIPsetEventhdlrFree(subproblem, eventhdlr, eventFreeBendersNodefocus) );
   assert(eventhdlr != NULL);

   /* calling an initial solve to put the problem into probing mode */
   SCIP_CALL( initialiseSubproblem(benders, set, probnumber, &success) );

   return SCIP_OKAY; /*lint !e438*/
}

/** checks whether the convex relaxation of the subproblem is sufficient to solve the original problem to optimality
 *
 * We check whether we can conclude that the CIP is actually an LP or a convex NLP.
 * To do this, we check that all variables are of continuous type and that every constraint is either handled by known
 * linear constraint handler (knapsack, linear, logicor, setppc, varbound) or the nonlinear constraint handler.
 * In the latter case, we also check whether the nonlinear constraint is convex.
 * Further, nonlinear constraints are only considered if an NLP solver interface is available, i.e., and NLP could
 * be solved.
 * If constraints are present that cannot be identified as linear or convex nonlinear, then we assume that the
 * problem is not convex, thus solving its LP or NLP relaxation will not be sufficient.
 */
static
SCIP_RETCODE checkSubproblemConvexity(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number, or -1 for the master problem */
   )
{
   SCIP* subproblem;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS* cons;
   SCIP_HASHMAP* assumevarfixed;
   SCIP_VAR** vars;
   int nvars;
   int nbinvars;
   int nintvars;
   int nimplintvars;
   int i;
   int j;
   SCIP_Bool convexcons;
   SCIP_Bool discretevar;
   SCIP_Bool isnonlinear;
   SCIP_CONSHDLR* linearconshdlrs[NLINEARCONSHDLRS];
   SCIP_CONSHDLR* nlconshdlr = NULL;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= -1 && probnumber < SCIPbendersGetNSubproblems(benders));

   assumevarfixed = NULL;
   if( probnumber >= 0 )
      subproblem = SCIPbendersSubproblem(benders, probnumber);
   else
      subproblem = set->scip;

   assert(subproblem != NULL);

   convexcons = FALSE;
   discretevar = FALSE;
   isnonlinear = FALSE;

   /* getting the number of integer and binary variables to determine the problem type */
   SCIP_CALL( SCIPgetVarsData(subproblem, &vars, &nvars, &nbinvars, &nintvars, &nimplintvars, NULL) );

   /* if there are any binary, integer or implied integer variables, then the subproblems is marked as non-convex */
   if( nbinvars != 0 || nintvars != 0 || nimplintvars != 0 )
   {
      discretevar = TRUE;
   }

   /* get pointers to linear constraints handlers, so can avoid string comparisons */
   linearconshdlrs[0] = SCIPfindConshdlr(subproblem, "knapsack");
   linearconshdlrs[1] = SCIPfindConshdlr(subproblem, "linear");
   linearconshdlrs[2] = SCIPfindConshdlr(subproblem, "logicor");
   linearconshdlrs[3] = SCIPfindConshdlr(subproblem, "setppc");
   linearconshdlrs[4] = SCIPfindConshdlr(subproblem, "varbound");

   /* Get pointer to the nonlinear constraint handler if we also have an NLP solver to solve NLPs.
    * If there is no NLP solver, but there are (convex) nonlinear constraints, then the LP relaxation of subproblems
    * will (currently) not be sufficient to solve subproblems to optimality. Thus, we also take the presence of convex
    * nonlinear constraints as signal for having to solve the CIP eventually, thus, by abuse of notation,
    * return not-convex here. In summary, we do not need to have a special look onto non-linear constraints
    * if no NLP solver is present, and can treat them as any other constraint that is not of linear type.
    */
   if( SCIPgetNNlpis(subproblem) > 0 )
   {
      nlconshdlr = SCIPfindConshdlr(subproblem, "nonlinear");
   }

   /* if the nonlinear constraint handler exists, then we create a hashmap of variables that can be assumed to be fixed.
    * These variables correspond to the copies of the master variables in the subproblem
    */
   if( probnumber >= 0 && nlconshdlr != NULL )
   {
      SCIP_VAR* mappedvar;

      SCIP_CALL( SCIPhashmapCreate(&assumevarfixed, SCIPblkmem(set->scip), SCIPgetNVars(subproblem)) );

      /* finding the subproblem variables that correspond to master variables */
      for( i = 0; i < nvars; i++ )
      {
         /* getting the corresponding master problem variable for the given variable */
         SCIP_CALL( SCIPbendersGetVar(benders, set, vars[i], &mappedvar, -1) );

         /* if the mapped variable is not NULL, then it must be stored as a possible fixed variable */
         if( mappedvar != NULL )
         {
            SCIP_CALL( SCIPhashmapInsert(assumevarfixed, vars[i], vars[i]) );
         }
      }
   }

   for( i = 0; i < SCIPgetNOrigConss(subproblem); ++i )
   {
      cons = SCIPgetOrigConss(subproblem)[i];
      conshdlr = SCIPconsGetHdlr(cons);

      for( j = 0; j < NLINEARCONSHDLRS; ++j )
         if( conshdlr == linearconshdlrs[j] )
            break;

      /* if linear constraint, then we are good */
      if( j < NLINEARCONSHDLRS )
      {
#ifdef SCIP_MOREDEBUG
         SCIPdebugMsg(subproblem, "subproblem <%s>: constraint <%s> is linear\n", SCIPgetProbName(subproblem), SCIPconsGetName(cons));
#endif
         continue;
      }

      /* if cons_nonlinear (and nlconshdlr != NULL), then check whether convex */
      if( conshdlr == nlconshdlr )
      {
         SCIP_Bool isconvex;
         SCIP_EXPRCURV curv;
         SCIP_Bool havelhs;
         SCIP_Bool haverhs;

         isnonlinear = TRUE;

         havelhs = !SCIPisInfinity(subproblem, -SCIPgetLhsNonlinear(cons));
         haverhs = !SCIPisInfinity(subproblem,  SCIPgetRhsNonlinear(cons));
         if( havelhs && haverhs )
         {
            isconvex = FALSE;
         }
         else
         {
            /* look at curvature stored in cons, though at this stage this will be unknown a.a. */
            curv = SCIPgetCurvatureNonlinear(cons);
            isconvex = ((!havelhs || (curv & SCIP_EXPRCURV_CONCAVE) == SCIP_EXPRCURV_CONCAVE)) &&
                ((!haverhs || (curv & SCIP_EXPRCURV_CONVEX) == SCIP_EXPRCURV_CONVEX));

            if( !isconvex )
            {
               /* if not found convex, compute curvature via nlhdlr_convex and decide again */

               /* make sure activities are up to date. SCIPhasExprCurvature currently assumes that this is already the case */
               SCIP_CALL( SCIPevalExprActivity(subproblem, SCIPgetExprNonlinear(cons)) );

               SCIP_CALL( SCIPhasExprCurvature(subproblem, SCIPgetExprNonlinear(cons), havelhs ? SCIP_EXPRCURV_CONCAVE : SCIP_EXPRCURV_CONVEX, &isconvex, assumevarfixed) );
            }
         }

         if( isconvex )
         {
#ifdef SCIP_MOREDEBUG
            SCIPdebugMsg(subproblem, "subproblem <%s>: nonlinear constraint <%s> is convex\n", SCIPgetProbName(subproblem), SCIPconsGetName(cons));
#endif
            continue;
         }
         else
         {
#ifdef SCIP_MOREDEBUG
            SCIPdebugMsg(subproblem, "subproblem <%s>: nonlinear constraint <%s> not convex\n", SCIPgetProbName(subproblem), SCIPconsGetName(cons));
#endif
            goto TERMINATE;
         }
      }

#ifdef SCIP_MOREDEBUG
      SCIPdebugMsg(subproblem, "subproblem <%s>: potentially nonconvex constraint <%s>\n", SCIPgetProbName(subproblem), SCIPconsGetName(cons));
#endif
      goto TERMINATE;
   }

   /* if we made it until here, then all constraints are known and convex */
   convexcons = TRUE;

TERMINATE:
   /* setting the flag for the convexity of the subproblem. If convexity doesn't need to be checked, then it is assumed
    * that the subproblems are convex. However, if there are discrete variables, then the problem must be set as
    * non-convex. The discrete master variables will be changed to continuous, but this will happen at the first call to
    * SCIPbendersSetupSubproblem
    */
   if( probnumber >= 0 )
   {
      convexcons = convexcons || !benders->checkconsconvexity;

      if( convexcons && !discretevar )
         SCIPbendersSetSubproblemType(benders, probnumber, SCIP_BENDERSSUBTYPE_CONVEXCONT);
      else if( convexcons && discretevar )
         SCIPbendersSetSubproblemType(benders, probnumber, SCIP_BENDERSSUBTYPE_CONVEXDIS);
      else if( !convexcons && !discretevar )
         SCIPbendersSetSubproblemType(benders, probnumber, SCIP_BENDERSSUBTYPE_NONCONVEXCONT);
      else if( !convexcons && discretevar )
         SCIPbendersSetSubproblemType(benders, probnumber, SCIP_BENDERSSUBTYPE_NONCONVEXDIS);
      else
         SCIPABORT();
   }

   /* setting the non-linear subproblem flag */
   if( probnumber >= 0 )
      SCIPbendersSetSubproblemIsNonlinear(benders, probnumber, isnonlinear);
   else
      SCIPbendersSetMasterIsNonlinear(benders, isnonlinear);

   if( probnumber >= 0 )
   {
      SCIPsetDebugMsg(set, "subproblem <%s> has been found to be of type %d\n", SCIPgetProbName(subproblem),
         SCIPbendersGetSubproblemType(benders, probnumber));
   }

   /* releasing the fixed variable hashmap */
   if( assumevarfixed != NULL )
      SCIPhashmapFree(&assumevarfixed);

   return SCIP_OKAY;
}

/** creates the subproblems and registers it with the Benders' decomposition struct */
static
SCIP_RETCODE createSubproblems(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP* subproblem;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_VAR* mastervar;
   SCIP_VAR** vars;
   int nvars;
   int nsubproblems;
   int i;
   int j;

   assert(benders != NULL);
   assert(set != NULL);

   /* if the subproblems have already been created, then they will not be created again. This is the case if the
    * transformed problem has been freed and then retransformed. The subproblems should only be created when the problem
    * is first transformed. */
   if( benders->subprobscreated )
      return SCIP_OKAY;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* creating all subproblems */
   for( i = 0; i < nsubproblems; i++ )
   {
      /* calling the create subproblem call back method */
      SCIP_CALL( benders->benderscreatesub(set->scip, benders, i) );

      subproblem = SCIPbendersSubproblem(benders, i);

      /* the subproblem SCIP instance could be set to NULL. This is because user defined subproblem solving methods
       * could be used that don't solve a SCIP instance. Thus, the following setup of the subproblem SCIP instance is
       * not required.
       *
       * NOTE: since the subproblems are supplied as NULL pointers, the internal convexity check can not be performed.
       * The user needs to explicitly specify the subproblem type.
       */
      if( subproblem != NULL )
      {
         /* setting global limits for the subproblems. This overwrites the limits set by the user */
         SCIP_CALL( SCIPsetIntParam(subproblem, "limits/maxorigsol", 0) );

         /* getting the number of integer and binary variables to determine the problem type */
         SCIP_CALL( SCIPgetVarsData(subproblem, &vars, &nvars, NULL, NULL, NULL, NULL) );

         /* The objective function coefficients of the master problem are set to zero. This is necessary for the Benders'
          * decomposition algorithm, since the cut methods and the objective function check assumes that the objective
          * coefficients of the master problem variables are zero.
          *
          * This only occurs if the Benders' decomposition is not a copy. It is assumed that the correct objective
          * coefficients are given during the first subproblem creation.
          *
          * If the subproblems were copied, then the master variables will be checked to ensure that they have a zero
          * objective value.
          */
         if( !benders->iscopy || benders->threadsafe )
         {
            SCIP_Bool objchanged = FALSE;

            assert(SCIPgetStage(subproblem) == SCIP_STAGE_PROBLEM);
            for( j = 0; j < nvars; j++ )
            {
               /* retrieving the master problem variable */
               SCIP_CALL( SCIPbendersGetVar(benders, set, vars[j], &mastervar, -1) );

               /* if mastervar is not NULL, then the subproblem variable has a corresponding master problem variable */
               if( mastervar != NULL && !SCIPisZero(subproblem, SCIPvarGetObj(vars[j])) )
               {
                  SCIPverbMessage(subproblem, SCIP_VERBLEVEL_FULL, NULL, "Benders' decomposition: Changing the objective "
                     "coefficient of copy of master problem variable <%s> in subproblem %d to zero.\n",
                     SCIPvarGetName(mastervar), i);
                  /* changing the subproblem variable objective coefficient to zero */
                  SCIP_CALL( SCIPchgVarObj(subproblem, vars[j], 0.0) );

                  objchanged = TRUE;
               }
            }

            if( objchanged )
            {
               SCIPverbMessage(subproblem, SCIP_VERBLEVEL_HIGH, NULL, "Benders' decomposition: Objective coefficients of "
                  "copied of master problem variables has been changed to zero.\n");
            }
         }

         /* changing all of the master problem variable to continuous. */
         SCIP_CALL( SCIPbendersChgMastervarsToCont(benders, set, i) );

         /* checking the convexity of the subproblem. The convexity of the subproblem indicates whether the convex
          * relaxation is a valid relaxation for the problem
          */
         SCIP_CALL( checkSubproblemConvexity(benders, set, i) );

         /* if the problem is convex and has nonlinear constraints, then slack variables must be added to each of the
          * constraints
             */
         if( benders->execfeasphase ||
            (SCIPbendersGetSubproblemType(benders, i) <= SCIP_BENDERSSUBTYPE_CONVEXDIS
             && SCIPbendersSubproblemIsNonlinear(benders, i)) )
         {
            /* the slack variables are only added to the subproblem once. If the initialisation methods are called from a
             * copy, then the slack variables are not re-added. Alternatively, if the copy must be threadsafe, then the
             * subproblems are created from scratch again, so the slack variables need to be added.
             */
            if( !benders->iscopy || benders->threadsafe )
            {
               SCIP_CALL( addSlackVarsToConstraints(benders, set, i) );
            }

            /* setting the flag to indicate that slack variables have been added to the subproblem constraints. This is only
             * set if the slack variables have been added at the request of the user.
             */
            if( benders->execfeasphase )
               benders->feasibilityphase = TRUE;
         }

         /* after checking the subproblem for convexity, if the subproblem has convex constraints and continuous variables,
          * then the problem is entered into probing mode. Otherwise, it is initialised as a CIP
          */
         if( SCIPbendersGetSubproblemType(benders, i) == SCIP_BENDERSSUBTYPE_CONVEXCONT )
         {
            /* if the user has not implemented a solve subproblem callback, then the subproblem solves are performed
             * internally. To be more efficient the subproblem is put into probing mode. */
            if( benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL
               && SCIPgetStage(subproblem) <= SCIP_STAGE_PROBLEM )
            {
               SCIP_CALL( initialiseLPSubproblem(benders, set, i) );
            }
         }
         else
         {
            SCIP_EVENTHDLRDATA* eventhdlrdata_mipnodefocus;
            SCIP_EVENTHDLRDATA* eventhdlrdata_upperbound;

            /* because the subproblems could be reused in the copy, the event handler is not created again. If the
             * threadsafe is TRUE, then it is assumed that the subproblems are not reused.
             * NOTE: This currently works with the benders_default implementation. It may not be very general. */
            if( benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL
               && (!benders->iscopy || benders->threadsafe) )
            {
               SCIP_CALL( SCIPallocBlockMemory(subproblem, &eventhdlrdata_mipnodefocus) );
               SCIP_CALL( SCIPallocBlockMemory(subproblem, &eventhdlrdata_upperbound) );

               SCIP_CALL( initEventhandlerData(subproblem, eventhdlrdata_mipnodefocus) );
               SCIP_CALL( initEventhandlerData(subproblem, eventhdlrdata_upperbound) );

               /* include the first LP solved event handler into the subproblem */
               SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, MIPNODEFOCUS_EVENTHDLR_NAME,
                     MIPNODEFOCUS_EVENTHDLR_DESC, eventExecBendersMipnodefocus, eventhdlrdata_mipnodefocus) );
               SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersMipnodefocus) );
               SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersMipnodefocus) );
               SCIP_CALL( SCIPsetEventhdlrExit(subproblem, eventhdlr, eventExitBendersMipnodefocus) );
               SCIP_CALL( SCIPsetEventhdlrFree(subproblem, eventhdlr, eventFreeBendersMipnodefocus) );
               assert(eventhdlr != NULL);

               /* include the upper bound interrupt event handler into the subproblem */
               SCIP_CALL( SCIPincludeEventhdlrBasic(subproblem, &eventhdlr, UPPERBOUND_EVENTHDLR_NAME,
                     UPPERBOUND_EVENTHDLR_DESC, eventExecBendersUpperbound, eventhdlrdata_upperbound) );
               SCIP_CALL( SCIPsetEventhdlrInitsol(subproblem, eventhdlr, eventInitsolBendersUpperbound) );
               SCIP_CALL( SCIPsetEventhdlrExitsol(subproblem, eventhdlr, eventExitsolBendersUpperbound) );
               SCIP_CALL( SCIPsetEventhdlrExit(subproblem, eventhdlr, eventExitBendersUpperbound) );
               SCIP_CALL( SCIPsetEventhdlrFree(subproblem, eventhdlr, eventFreeBendersUpperbound) );
               assert(eventhdlr != NULL);
            }
         }
      }
      else
      {
         /* a user must specify the subproblem type if they are not supplying a SCIP instance. */
         if( SCIPbendersGetSubproblemType(benders, i) == SCIP_BENDERSSUBTYPE_UNKNOWN )
         {
            SCIPerrorMessage("If the subproblem is set to NULL, then the subproblem type must be specified.\n");
            SCIPerrorMessage("In the subproblem creation callback, call SCIPbendersSetSubproblemType with the appropriate problem type.\n");

            return SCIP_ERROR;
         }
      }
   }

   /* checking the convexity of the master problem. This information is useful for the cut generation methods, such as
    * non-good and integer cuts
    */
   SCIP_CALL( checkSubproblemConvexity(benders, set, -1) );

   benders->subprobscreated = TRUE;

   return SCIP_OKAY;
}


/** initializes Benders' decomposition */
SCIP_RETCODE SCIPbendersInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   if( benders->initialized )
   {
      SCIPerrorMessage("Benders' decomposition <%s> already initialized\n", benders->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(benders->setuptime);
      SCIPclockReset(benders->bendersclock);

      benders->ncalls = 0;
      benders->ncutsfound = 0;
      benders->ntransferred = 0;
   }

   /* start timing */
   SCIPclockStart(benders->setuptime, set);

   if( benders->bendersinit != NULL )
   {
      SCIP_CALL( benders->bendersinit(set->scip, benders) );
   }

   benders->initialized = TRUE;

   /* if the Benders' decomposition is a copy, then the auxiliary variables already exist. So they are registered with
    * the Benders' decomposition struct during the init stage. If the Benders' decomposition is not a copy, then the
    * auxiliary variables need to be created, which occurs in the initpre stage
    */
   if( benders->iscopy )
   {
      /* the copied auxiliary variables must be assigned to the target Benders' decomposition */
      SCIP_CALL( assignAuxiliaryVariables(set->scip, benders) );
   }

   /* creates the subproblems and sets up the probing mode for LP subproblems. This function calls the benderscreatesub
    * callback. */
   SCIP_CALL( createSubproblems(benders, set) );

   /* storing the solution tolerance set by the SCIP parameters */
   SCIP_CALL( SCIPsetGetRealParam(set, "benders/solutiontol", &benders->solutiontol) );

   /* allocating memory for the stored constraints array */
   if( benders->storedcutssize == 0 )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(SCIPblkmem(set->scip), &benders->storedcuts, BENDERS_ARRAYSIZE) );
      benders->storedcutssize = BENDERS_ARRAYSIZE;
      benders->nstoredcuts = 0;
   }

   /* initialising the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutInit(benders->benderscuts[i], set) );
   }

   /* stop timing */
   SCIPclockStop(benders->setuptime, set);

   return SCIP_OKAY;
}


/** Transfers Benders' cuts that were generated while solving a sub-SCIP to the original SCIP instance. This involves
 *  creating a constraint/cut that is equivalent to the generated cut in the sub-SCIP. This new constraint/cut is then
 *  added to the original SCIP instance.
 */
static
SCIP_RETCODE createAndAddTransferredCut(
   SCIP*                 sourcescip,         /**< the source SCIP from when the Benders' decomposition was copied */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure of the sub SCIP */
   SCIP_VAR**            vars,               /**< the variables from the source constraint */
   SCIP_Real*            vals,               /**< the coefficients of the variables in the source constriant */
   SCIP_Real             lhs,                /**< the LHS of the source constraint */
   SCIP_Real             rhs,                /**< the RHS of the source constraint */
   int                   nvars               /**< the number of variables in the source constraint */
   )
{
   SCIP_BENDERS* sourcebenders;     /* the Benders' decomposition of the source SCIP */
   SCIP_CONSHDLR* consbenders;      /* a helper variable for the Benders' decomposition constraint handler */
   SCIP_CONS* transfercons = NULL;  /* the constraint that is generated to transfer the constraints/cuts */
   SCIP_ROW* transfercut = NULL;    /* the cut that is generated to transfer the constraints/cuts */
   SCIP_VAR* sourcevar;             /* the source variable that will be added to the transferred cut */
   SCIP_VAR* origvar;
   SCIP_Real scalar;
   SCIP_Real constant;
   char cutname[SCIP_MAXSTRLEN];    /* the name of the transferred cut */
   int i;
   SCIP_Bool fail;

   assert(sourcescip != NULL);
   assert(benders != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   /* retrieving the source Benders' decomposition structure */
   sourcebenders = SCIPfindBenders(sourcescip, SCIPbendersGetName(benders));

   /* retrieving the Benders' decomposition constraint handler */
   consbenders = SCIPfindConshdlr(sourcescip, "benders");

   /* setting the name of the transferred cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "transferredcut_%d",
      SCIPbendersGetNTransferredCuts(sourcebenders) );

   /* TODO: It could be more efficient to pass an updated vars array with the vals array to the
    * SCIPcreateConsBasicLinear/SCIPcreateEmptyRowConshdlr. This should be implemented to improve the performance of the
    * Large Neighbourhood Benders Search.
    */

   /* creating an empty row/constraint for the transferred cut */
   if( sourcebenders->cutsasconss )
   {
      SCIP_CALL( SCIPcreateConsBasicLinear(sourcescip, &transfercons, cutname, 0, NULL, NULL, lhs, rhs) );
      SCIP_CALL( SCIPsetConsRemovable(sourcescip, transfercons, TRUE) );
   }
   else
   {
      SCIP_CALL( SCIPcreateEmptyRowConshdlr(sourcescip, &transfercut, consbenders, cutname, lhs, rhs, FALSE,
            FALSE, TRUE) );
   }

   fail = FALSE;
   for( i = 0; i < nvars; i++ )
   {
      /* getting the original variable for the transformed variable */
      origvar = vars[i];
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      /* getting the source var from the hash map */
      sourcevar = (SCIP_VAR*) SCIPhashmapGetImage(benders->mastervarsmap, origvar);

      /* if the source variable is not found, then the mapping in incomplete. So the constraint can not be
       * transferred. */
      if( sourcevar == NULL )
      {
         fail = TRUE;
         break;
      }

      if( sourcebenders->cutsasconss )
      {
         assert( transfercons != NULL );
         SCIP_CALL( SCIPaddCoefLinear(sourcescip, transfercons, sourcevar, vals[i]) );    /*lint !e644*/
      }
      else
      {
         assert( transfercut != NULL );
         SCIP_CALL( SCIPaddVarToRow(sourcescip, transfercut, sourcevar, vals[i]) );       /*lint !e644*/
      }
   }

   /* if all of the source variables were found to generate the cut */
   if( !fail )
   {
      if( sourcebenders->cutsasconss )
      {
         SCIP_CALL( SCIPaddCons(sourcescip, transfercons) );
      }
      else
      {
         SCIP_CALL( SCIPaddPoolCut(sourcescip, transfercut) );
      }

      sourcebenders->ntransferred++;
   }

   /* release the row/constraint */
   if( sourcebenders->cutsasconss )
   {
      /* only release if the creation of the constraint failed. */
      SCIP_CALL( SCIPreleaseCons(sourcescip, &transfercons) );
   }
   else
   {
      SCIP_CALL( SCIPreleaseRow(sourcescip, &transfercut) );
   }

   return SCIP_OKAY;
}


/** transfers the cuts generated in a subscip to the source scip */
static
SCIP_RETCODE transferBendersCuts(
   SCIP*                 sourcescip,         /**< the source SCIP from when the Benders' decomposition was copied */
   SCIP*                 subscip,            /**< the sub SCIP where the Benders' cuts were generated */
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure of the sub SCIP */
   )
{
   SCIP_BENDERS* sourcebenders;     /* the Benders' decomposition of the source SCIP */
   SCIP_VAR** vars;                 /* the variables of the added constraint/row */
   SCIP_Real* vals;                 /* the values of the added constraint/row */
   SCIP_Real lhs;                   /* the LHS of the added constraint/row */
   SCIP_Real rhs;                   /* the RHS of the added constraint/row */
   int naddedcuts;
   int nvars;
   int i;

   assert(subscip != NULL);
   assert(benders != NULL);

   /* retrieving the source Benders' decomposition structure */
   sourcebenders = SCIPfindBenders(sourcescip, SCIPbendersGetName(benders));

   /* exit if the cuts should not be transferred from the sub SCIP to the source SCIP. */
   if( !sourcebenders->transfercuts || benders->mastervarsmap == NULL )
      return SCIP_OKAY;

   /* retrieving the number of stored Benders' cuts */
   naddedcuts =  SCIPbendersGetNStoredCuts(benders);

   /* looping over all added cuts to construct the cut for the source scip */
   for( i = 0; i < naddedcuts; i++ )
   {
      /* collecting the variable information from the constraint */
      SCIP_CALL( SCIPbendersGetStoredCutData(benders, i, &vars, &vals, &lhs, &rhs, &nvars) );

      if( nvars > 0 )
      {
         /* create and add the cut to be transferred from the sub SCIP to the source SCIP */
         SCIP_CALL( createAndAddTransferredCut(sourcescip, benders, vars, vals, lhs, rhs, nvars) );
      }
   }

   return SCIP_OKAY;
}


/** calls exit method of Benders' decomposition */
SCIP_RETCODE SCIPbendersExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int nsubproblems;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   if( !benders->initialized )
   {
      SCIPerrorMessage("Benders' decomposition <%s> not initialized\n", benders->name);
      return SCIP_INVALIDCALL;
   }

   /* start timing */
   SCIPclockStart(benders->setuptime, set);

   if( benders->bendersexit != NULL )
   {
      SCIP_CALL( benders->bendersexit(set->scip, benders) );
   }

   /* if the Benders' decomposition is a copy, then is a variable mapping was provided, then the generated cuts will
    * be transferred to the source scip
    */
   if( benders->iscopy && benders->mastervarsmap != NULL )
   {
      SCIP_CALL( transferBendersCuts(benders->sourcescip, set->scip, benders) );
   }

   /* releasing the stored constraints */
   for( i = benders->nstoredcuts - 1; i >= 0; i-- )
   {
      SCIPfreeBlockMemoryArray(set->scip, &benders->storedcuts[i]->vals, benders->storedcuts[i]->nvars);
      SCIPfreeBlockMemoryArray(set->scip, &benders->storedcuts[i]->vars, benders->storedcuts[i]->nvars);
      SCIPfreeBlockMemory(set->scip, &benders->storedcuts[i]); /*lint !e866*/
   }

   BMSfreeBlockMemoryArray(SCIPblkmem(set->scip), &benders->storedcuts, benders->storedcutssize);
   benders->storedcutssize = 0;
   benders->nstoredcuts = 0;

   /* releasing all of the auxiliary variables */
   nsubproblems = SCIPbendersGetNSubproblems(benders);
   for( i = 0; i < nsubproblems; i++ )
   {
      /* it is possible that the master problem is not solved. As such, the auxiliary variables will not be created. So
       * we don't need to release the variables
       */
      if( benders->auxiliaryvars[i] != NULL )
      {
         /* we need to remove the locks from the auxiliary variables. This will be called always for the highest priority
          * Benders' plugin and others if the auxiliary variables are not shared
          */
         if( !benders->iscopy && SCIPvarGetNLocksDown(benders->auxiliaryvars[i]) > 0 )
            SCIP_CALL( SCIPaddVarLocksType(set->scip, benders->auxiliaryvars[i], SCIP_LOCKTYPE_MODEL, -1, 0) );

         SCIP_CALL( SCIPreleaseVar(set->scip, &benders->auxiliaryvars[i]) );
      }
   }

   /* if a corepoint has been used for cut strengthening, then this needs to be freed */
   if( benders->corepoint != NULL )
   {
      SCIP_CALL( SCIPfreeSol(set->scip, &benders->corepoint) );
   }

   /* calling the exit method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutExit(benders->benderscuts[i], set) );
   }

   benders->initialized = FALSE;

   /* stop timing */
   SCIPclockStop(benders->setuptime, set);

   return SCIP_OKAY;
}

/** Checks whether a subproblem is independent. */
static
SCIP_RETCODE checkSubproblemIndependence(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nsubproblems;
   int i;
   int j;

   assert(scip != NULL);
   assert(benders != NULL);

   /* retrieving the master problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* looping over all subproblems to check whether there exists at least one master problem variable */
   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_Bool independent = FALSE;

      /* if there are user defined solving or freeing functions, then it is not possible to declare the independence of
       * the subproblems.
       */
      if( benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL
         && benders->bendersfreesub == NULL )
      {
         independent = TRUE;

         for( j = 0; j < nvars; j++ )
         {
            SCIP_VAR* subprobvar;

            /* getting the subproblem problem variable corresponding to the master problem variable */
            SCIP_CALL( SCIPgetBendersSubproblemVar(scip, benders, vars[j], &subprobvar, i) );

            /* if the subporblem variable is not NULL, then the subproblem depends on the master problem */
            if( subprobvar != NULL )
            {
               independent = FALSE;
               break;
            }
         }

         /* setting the independent flag */
         SCIPbendersSetSubproblemIsIndependent(benders, i, independent);
      }
   }

   return SCIP_OKAY;
}

/** informs the Benders' decomposition that the presolving process is being started */
SCIP_RETCODE SCIPbendersInitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   /* if the Benders' decomposition is the original, then the auxiliary variables need to be created. If the Benders'
    * decomposition is a copy, then the auxiliary variables already exist. The assignment of the auxiliary variables
    * occurs in bendersInit
    */
   if( !benders->iscopy )
   {
      /* check the subproblem independence. This check is only performed if the user has not implemented a solve
       * subproblem function.
       */
      if( benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL )
        SCIP_CALL( checkSubproblemIndependence(set->scip, benders) );

      /* adding the auxiliary variables to the master problem */
      SCIP_CALL( addAuxiliaryVariablesToMaster(set->scip, benders) );
   }

   /* call presolving initialization method of Benders' decomposition */
   if( benders->bendersinitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersinitpre(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   return SCIP_OKAY;
}


/** informs the Benders' decomposition that the presolving process has completed */
SCIP_RETCODE SCIPbendersExitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(stat != NULL);

   /* call presolving  deinitialization method of Benders' decomposition */
   if( benders->bendersexitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexitpre(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs Benders' decomposition that the branch and bound process is being started */
SCIP_RETCODE SCIPbendersInitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   /* call solving process initialization method of Benders' decomposition */
   if( benders->bendersinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersinitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* calling the initsol method for the Benders' cuts */
   SCIPbendersSortBenderscuts(benders);
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutInitsol(benders->benderscuts[i], set) );
   }

   return SCIP_OKAY;
}

/** informs Benders' decomposition that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbendersExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int nsubproblems;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   nsubproblems = SCIPbendersGetNSubproblems(benders);
   /* freeing all subproblems that are independent, this is because they have not bee freed during the subproblem
    * solving loop.
    */
   for( i = 0; i < nsubproblems; i++ )
   {
      if( SCIPbendersSubproblemIsIndependent(benders, i) )
      {
         /* disabling the independence of the subproblem so that it can be freed */
         SCIPbendersSetSubproblemIsIndependent(benders, i, FALSE);

         /* freeing the independent subproblem */
         SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, i) );
      }
   }

   /* call solving process deinitialization method of Benders' decomposition */
   if( benders->bendersexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benders->setuptime, set);

      SCIP_CALL( benders->bendersexitsol(set->scip, benders) );

      /* stop timing */
      SCIPclockStop(benders->setuptime, set);
   }

   /* sorting the Benders' decomposition cuts in order of priority. Only a single cut is generated for each subproblem
    * per solving iteration. This is particularly important in the case of the optimality and feasibility cuts. Since
    * these work on two different solutions to the subproblem, it is not necessary to generate both cuts. So, once the
    * feasibility cut is generated, then no other cuts will be generated.
    */
   SCIPbendersSortBenderscuts(benders);

   /* calling the exitsol method for the Benders' cuts */
   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      SCIP_CALL( SCIPbenderscutExitsol(benders->benderscuts[i], set) );
   }

   return SCIP_OKAY;
}

/** activates Benders' decomposition such that it is called in LP solving loop */
SCIP_RETCODE SCIPbendersActivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nsubproblems        /**< the number subproblems used in this decomposition */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT || set->stage == SCIP_STAGE_PROBLEM);

   if( !benders->active )
   {
      benders->active = TRUE;
      set->nactivebenders++;
      set->benderssorted = FALSE;

      benders->nsubproblems = nsubproblems;
      benders->nactivesubprobs = nsubproblems;
      benders->prevlowerbound = -SCIPsetInfinity(set);
      benders->strengthenround = FALSE;

      /* allocating memory for the subproblems arrays */
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subproblems, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->auxiliaryvars, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->solvestat, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobobjval, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->bestsubprobobjval, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subproblowerbound, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobtype, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobisconvex, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobisnonlinear, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobsetup, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->indepsubprob, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->subprobenabled, benders->nsubproblems) );
      SCIP_ALLOC( BMSallocMemoryArray(&benders->mastervarscont, benders->nsubproblems) );

      /* creating the priority queue for the subproblem solving status */
      SCIP_CALL( SCIPpqueueCreate(&benders->subprobqueue, benders->nsubproblems, 1.1,
            benders->benderssubcomp == NULL ? benderssubcompdefault : benders->benderssubcomp, NULL) );

      for( i = 0; i < benders->nsubproblems; i++ )
      {
         SCIP_SUBPROBLEMSOLVESTAT* solvestat;

         benders->subproblems[i] = NULL;
         benders->auxiliaryvars[i] = NULL;
         benders->subprobobjval[i] = SCIPsetInfinity(set);
         benders->bestsubprobobjval[i] = SCIPsetInfinity(set);
         benders->subproblowerbound[i] = -SCIPsetInfinity(set);
         benders->subprobtype[i] = SCIP_BENDERSSUBTYPE_UNKNOWN;
         benders->subprobisconvex[i] = FALSE;
         benders->subprobisnonlinear[i] = FALSE;
         benders->subprobsetup[i] = FALSE;
         benders->indepsubprob[i] = FALSE;
         benders->subprobenabled[i] = TRUE;
         benders->mastervarscont[i] = FALSE;

         /* initialising the subproblem solving status */
         SCIP_ALLOC( BMSallocMemory(&solvestat) );
         solvestat->idx = i;
         solvestat->ncalls = 0;
         solvestat->avgiter = 0;
         benders->solvestat[i] = solvestat;

         /* inserting the initial elements into the priority queue */
         SCIP_CALL( SCIPpqueueInsert(benders->subprobqueue, benders->solvestat[i]) );
      }

      /* adding an eventhandler for updating the lower bound when the root node is solved. */
      eventhdlrdata = (SCIP_EVENTHDLRDATA*)benders;

      /* include event handler into SCIP */
      SCIP_CALL( SCIPincludeEventhdlrBasic(set->scip, &eventhdlr, NODESOLVED_EVENTHDLR_NAME, NODESOLVED_EVENTHDLR_DESC,
            eventExecBendersNodesolved, eventhdlrdata) );
      SCIP_CALL( SCIPsetEventhdlrInitsol(set->scip, eventhdlr, eventInitsolBendersNodesolved) );
      assert(eventhdlr != NULL);
   }

   return SCIP_OKAY;
}

/** deactivates Benders' decomposition such that it is no longer called in LP solving loop */
SCIP_RETCODE SCIPbendersDeactivate(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT || set->stage == SCIP_STAGE_PROBLEM);

   if( benders->active )
   {
      int nsubproblems;

      nsubproblems = SCIPbendersGetNSubproblems(benders);

#ifndef NDEBUG
      /* checking whether the auxiliary variables and subproblems are all NULL */
      for( i = 0; i < nsubproblems; i++ )
         assert(benders->auxiliaryvars[i] == NULL);
#endif

      /* if the subproblems were created by the Benders' decomposition core, then they need to be freed */
      if( benders->freesubprobs )
      {
         for( i = SCIPbendersGetNSubproblems(benders) - 1; i >= 0; i-- )
         {
            SCIP* subproblem = SCIPbendersSubproblem(benders, i);
            SCIP_CALL( SCIPfree(&subproblem) );
         }
      }

      benders->active = FALSE;
      set->nactivebenders--;
      set->benderssorted = FALSE;

      /* freeing the priority queue memory */
      SCIPpqueueFree(&benders->subprobqueue);

      for( i = nsubproblems - 1; i >= 0; i-- )
         BMSfreeMemory(&benders->solvestat[i]);

      /* freeing the memory allocated during the activation of the Benders' decomposition */
      BMSfreeMemoryArray(&benders->mastervarscont);
      BMSfreeMemoryArray(&benders->subprobenabled);
      BMSfreeMemoryArray(&benders->indepsubprob);
      BMSfreeMemoryArray(&benders->subprobsetup);
      BMSfreeMemoryArray(&benders->subprobisnonlinear);
      BMSfreeMemoryArray(&benders->subprobisconvex);
      BMSfreeMemoryArray(&benders->subprobtype);
      BMSfreeMemoryArray(&benders->subproblowerbound);
      BMSfreeMemoryArray(&benders->bestsubprobobjval);
      BMSfreeMemoryArray(&benders->subprobobjval);
      BMSfreeMemoryArray(&benders->auxiliaryvars);
      BMSfreeMemoryArray(&benders->solvestat);
      BMSfreeMemoryArray(&benders->subproblems);
   }

   return SCIP_OKAY;
}

/** returns whether the given Benders' decomposition is in use in the current problem */
SCIP_Bool SCIPbendersIsActive(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   assert(benders != NULL);

   return benders->active;
}

/** updates the lower bound for all auxiliary variables. This is called if the first LP enforced is unbounded. */
static
SCIP_RETCODE updateAuxiliaryVarLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RESULT*          result              /**< the result from updating the auxiliary variable lower bound */
   )
{
   int nsubproblems;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   (*result) = SCIP_DIDNOTRUN;

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   for( i = 0; i < nsubproblems; i++ )
   {
      SCIP_VAR* auxiliaryvar;
      SCIP_Real lowerbound;
      SCIP_Bool infeasible;

      infeasible = FALSE;

      /* computing the lower bound of the subproblem by solving it without any variable fixings */
      SCIP_CALL( SCIPbendersComputeSubproblemLowerbound(benders, set, i, &lowerbound, &infeasible) );

      /* if the subproblem is infeasible, then the original problem is infeasible */
      if( infeasible )
      {
         (*result) = SCIP_INFEASIBLE;
         break;
      }

      /* retrieving the auxiliary variable */
      auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, i);

      /* only update the lower bound if it is greater than the current lower bound */
      if( SCIPsetIsGT(set, lowerbound, SCIPvarGetLbGlobal(auxiliaryvar)) )
      {
         SCIPsetDebugMsg(set, "Tightened lower bound of <%s> to %g\n", SCIPvarGetName(auxiliaryvar), lowerbound);
         /* updating the lower bound of the auxiliary variable */
         SCIP_CALL( SCIPchgVarLb(set->scip, auxiliaryvar, lowerbound) );
         (*result) = SCIP_REDUCEDDOM;
      }

      /* stores the lower bound for the subproblem */
      SCIPbendersUpdateSubproblemLowerbound(benders, i, lowerbound);
   }

   return SCIP_OKAY;
}

/** sets the core point used for cut strengthening. If the strenghtenintpoint is set to 'i', then the core point is
 *  reinitialised each time the incumbent is updated
 */
static
SCIP_RETCODE setAndUpdateCorePoint(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_SOL* bestsol;

   assert(scip != NULL);
   assert(benders != NULL);

   /* if the core point is not NULL and the interior point is not reinitialised, then nothing is done  */
   if( benders->corepoint != NULL && benders->strengthenintpoint != 'i' )
      return SCIP_OKAY;

   bestsol = SCIPgetBestSol(scip);

   /* if the core point should be updated, then this only happens if the incumbent solution has been updated */
   if( benders->strengthenintpoint == 'i' && benders->initcorepoint == bestsol )
      return SCIP_OKAY;

   /* if a corepoint has been used for cut strengthening, then this needs to be freed */
   if( benders->corepoint != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &benders->corepoint) );
   }

   switch( benders->strengthenintpoint )
   {
      SCIP_VAR** vars;
      SCIP_Real timelimit;
      int nvars;
      int i;

      case 'l':
         SCIP_CALL( SCIPcreateLPSol(scip, &benders->corepoint, NULL) );
         SCIP_CALL( SCIPunlinkSol(scip, benders->corepoint) );
         break;
      case 'f':
      case 'i':
         SCIP_CALL( SCIPcreateSolCopy(scip, &benders->corepoint, bestsol) );
         SCIP_CALL( SCIPunlinkSol(scip, benders->corepoint) );
         benders->initcorepoint = bestsol;
         break;
      case 'r':
         /* prepare time limit */
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if ( ! SCIPisInfinity(scip, timelimit) )
            timelimit -= SCIPgetSolvingTime(scip);

         /* if there is time remaining, then compute the relative interior point. Otherwise, return the LP solution */
         if ( timelimit > 0.0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Computing relative interior point (time limit: %g, iter limit: %d) ...\n", timelimit, INT_MAX);
            SCIP_CALL( SCIPcomputeLPRelIntPoint(scip, TRUE, FALSE, timelimit, INT_MAX, &benders->corepoint) );
         }
         else
         {
            SCIP_CALL( SCIPcreateLPSol(scip, &benders->corepoint, NULL) );
            SCIP_CALL( SCIPunlinkSol(scip, benders->corepoint) );
         }
         break;
      case 'z':
         SCIP_CALL( SCIPcreateSol(scip, &benders->corepoint, NULL) );
         break;
      case 'o':
         SCIP_CALL( SCIPcreateSol(scip, &benders->corepoint, NULL) );

         /* getting the variable data so that the  */
         SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

         /* setting all variable values to 1.0 */
         for( i = 0; i < nvars; i++ )
         {
            SCIP_CALL( SCIPsetSolVal(scip, benders->corepoint, vars[i], 1.0) );
         }
         break;
      default:
         SCIP_CALL( SCIPcreateLPSol(scip, &benders->corepoint, NULL) );
         SCIP_CALL( SCIPunlinkSol(scip, benders->corepoint) );
   }

   return SCIP_OKAY;
}

/** performs cut strengthening by using an interior solution to generate cuts */
static
SCIP_RETCODE performInteriorSolCutStrengthening(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint,           /**< are the subproblems called during a check/enforce of integer sols? */
   SCIP_Bool             perturbsol,         /**< should the solution be perturbed to escape infeasibility? */
   SCIP_Bool*            auxviol,            /**< set to TRUE only if the solution is feasible but the aux vars are violated */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_Bool*            skipsolve,          /**< should the main solve be skipped as a result of this strengthening? */
   SCIP_RESULT*          result              /**< result of the pricing process */
   )
{
   SCIP_SOL* sepapoint;
   SCIP_VAR** vars;
   int prevcutsfound;
   int nvars;
   int i;

   assert(benders != NULL);
   assert(set != NULL);

   (*result) = SCIP_DIDNOTRUN;
   (*skipsolve) = FALSE;

   /* the cut stabilisation is only performed when enforcing LP solutions. The solution is not NULL if the stabilisation
    * is currently being performed. It is important to avoid recursion
    */
   if( type != SCIP_BENDERSENFOTYPE_LP || sol != NULL )
      return SCIP_OKAY;

   /* checking if a change to the lower bound has occurred */
   if( SCIPsetIsGT(set, SCIPgetLowerbound(set->scip), benders->prevlowerbound)
      || SCIPgetCurrentNode(set->scip) != benders->prevnode )
   {
      benders->prevnode = SCIPgetCurrentNode(set->scip);
      benders->prevlowerbound = SCIPgetLowerbound(set->scip);
      benders->noimprovecount = 0;
   }
   else
      benders->noimprovecount++;

   /* if the number of iterations without improvement exceeds 3*noimprovelimit, then the no stabilisation is performed
    */
   if( benders->noimprovecount > 3*benders->noimprovelimit )
      return SCIP_OKAY;

   /* if there is no incumbent solution, then it is not possible to create the core point and hence the strengthening
    * can not be performed
    */
   if( SCIPgetBestSol(set->scip) == NULL )
      return SCIP_OKAY;

   /* if no LP iterations have been performed since the last call of the cut strenghtening, then the strengthening is
    * aborted
    */
   if( benders->prevnlpiter == SCIPgetNLPIterations(set->scip) )
      return SCIP_OKAY;

   benders->prevnlpiter = SCIPgetNLPIterations(set->scip);

   /* if the separation point solution is NULL, then we create the solution using the current LP relaxation. */
   SCIP_CALL( setAndUpdateCorePoint(set->scip, benders) );

   /* creating the separation point
    * TODO: This could be a little to memory heavy, it may be better just to create the separation point once and then
    * update it each time.
    */
   SCIP_CALL( SCIPcreateLPSol(set->scip, &sepapoint, NULL) );
   SCIP_CALL( SCIPunlinkSol(set->scip, sepapoint) );

   SCIP_CALL( SCIPgetVarsData(set->scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(vars != NULL);

   /* creating a solution that is a convex combination of the LP solution and the separation point */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_VAR* subvar;
      SCIP_Real corepointval;
      SCIP_Real lpsolval;
      SCIP_Real newsolval;
      int j;

      corepointval = SCIPgetSolVal(set->scip, benders->corepoint, vars[i]);
      lpsolval = SCIPgetSolVal(set->scip, sol, vars[i]);
      newsolval = lpsolval;

      /* checking whether the master variable is mapped to any subproblem variables */
      subvar = NULL;
      j = 0;
      while( subvar == NULL && j < SCIPgetBendersNSubproblems(set->scip, benders)  )
      {
         SCIP_CALL( SCIPgetBendersSubproblemVar(set->scip, benders, vars[i], &subvar, j) );
         j++;
      }

      /* if the variable is a linking variable and it is not fixed, then a convex combination with the corepoint is
       * computed.
       */
      if( subvar != NULL && SCIPvarGetStatus(vars[i]) != SCIP_VARSTATUS_FIXED )
      {
         /* if the number of iterations without improvement exceeds noimprovelimit, then no convex combination is
          * created
          */
         if( !perturbsol && benders->noimprovecount <= benders->noimprovelimit )
         {
            newsolval = lpsolval*benders->convexmult + corepointval*(1 - benders->convexmult);

            /* updating the core point */
            SCIP_CALL( SCIPsetSolVal(set->scip, benders->corepoint, vars[i], newsolval) );
         }

         /* if the number of iterations without improvement is less than 2*noimprovelimit, then perturbation is
          * performed
          * TODO: This should be a random vector!!!!
          */
         if( perturbsol || benders->noimprovecount <= 2*benders->noimprovelimit )
            newsolval += benders->perturbeps;
      }

      /* updating the separation point */
      SCIP_CALL( SCIPsetSolVal(set->scip, sepapoint, vars[i], newsolval) );
   }

   /* storing the number of cuts found */
   prevcutsfound = SCIPbendersGetNCutsFound(benders);

   SCIPsetDebugMsg(set, "solving Benders' decomposition subproblems with stabilised point.\n");

   /* calling the subproblem solving method to generate cuts from the separation solution */
   SCIP_CALL( SCIPsolveBendersSubproblems(set->scip, benders, sepapoint, result, infeasible, auxviol, type, checkint) );

   SCIPsetDebugMsg(set, "solved Benders' decomposition subproblems with stabilised point. noimprovecount %d result %d\n",
      benders->noimprovecount, (*result));

   /* if constraints were added, then the main Benders' solving loop is skipped. */
   if( !(*infeasible) && ((*result) == SCIP_CONSADDED || (*result) == SCIP_SEPARATED) )
      (*skipsolve) = TRUE;

   /* capturing cut strengthening statistics */
   benders->nstrengthencalls++;
   benders->nstrengthencuts += (SCIPbendersGetNCutsFound(benders) - prevcutsfound);

   /* if no cuts were added, then the strengthening round is marked as failed */
   if( SCIPbendersGetNCutsFound(benders) == prevcutsfound )
      benders->nstrengthenfails++;

   /* freeing the sepapoint solution */
   SCIP_CALL( SCIPfreeSol(set->scip, &sepapoint) );

   return SCIP_OKAY;
}


/** Returns whether only the convex relaxations will be checked in this solve loop
 *  when Benders' is used in the LNS heuristics, only the convex relaxations of the master/subproblems are checked,
 *  i.e. no integer cuts are generated. In this case, then Benders' decomposition is performed under the assumption
 *  that all subproblems are convex relaxations.
 */
SCIP_Bool SCIPbendersOnlyCheckConvexRelax(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_Bool             subscipsoff         /**< flag indicating whether plugins using sub-SCIPs are deactivated */
   )
{
   return benders->iscopy && benders->lnscheck && subscipsoff;
}

/** returns the number of subproblems that will be checked in this iteration */
static
int numSubproblemsToCheck(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERSENFOTYPE  type                /**< the type of solution being enforced */
   )
{
   if( benders->ncalls == 0 || type == SCIP_BENDERSENFOTYPE_CHECK
      || SCIPbendersOnlyCheckConvexRelax(benders, SCIPsetGetSubscipsOff(set)) )
      return SCIPbendersGetNSubproblems(benders);
   else
      return (int) SCIPsetCeil(set, (SCIP_Real) SCIPbendersGetNSubproblems(benders)*benders->subprobfrac);
}

/** returns whether the solving of the given subproblem needs to be executed */
static
SCIP_Bool subproblemIsActive(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem index */
   )
{
   return (!SCIPbendersSubproblemIsIndependent(benders, probnumber)
      && SCIPbendersSubproblemIsEnabled(benders, probnumber));
}

/** creates an ordered list of subproblem indices to be solved */
static
void createSolveSubproblemIndexList(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   int**                 solveidx,           /**< a list of subproblem indices to the solved in the current iteration */
   int*                  nsolveidx           /**< the number of subproblem indices in the list */
   )
{
   int nsubproblems;
   int numtocheck;
   int subproblemcount;

   assert(benders != NULL);
   assert(set != NULL);
   assert((*solveidx) != NULL);
   assert(nsolveidx != NULL);
   assert(SCIPpqueueNElems(benders->subprobqueue) <= SCIPbendersGetNSubproblems(benders));

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* it is possible to only solve a subset of subproblems. This is given by a parameter. */
   numtocheck = numSubproblemsToCheck(benders, set, type);

   (*nsolveidx) = 0;

   subproblemcount = 0;
   while( subproblemcount < nsubproblems && subproblemcount < numtocheck )
   {
      SCIP_SUBPROBLEMSOLVESTAT* solvestat;

      solvestat = (SCIP_SUBPROBLEMSOLVESTAT*)SCIPpqueueRemove(benders->subprobqueue);
      (*solveidx)[(*nsolveidx)] = solvestat->idx;
      (*nsolveidx)++;

      subproblemcount++;
   }
}

/** updates the subproblem solving statistics and inserts the indices into the queue */
static
SCIP_RETCODE updateSubproblemStatQueue(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int*                  solveidx,           /**< the list of indices of subproblems that were solved */
   int                   nsolveidx,          /**< the number of subproblem indices */
   SCIP_Bool             updatestat          /**< should the statistics be updated */
   )
{
   int i;

   assert(benders != NULL);
   assert(solveidx != NULL);

   for( i = 0; i < nsolveidx; i++ )
   {
      SCIP* subproblem;
      SCIP_SUBPROBLEMSOLVESTAT* solvestat;

      subproblem = SCIPbendersSubproblem(benders, solveidx[i]);
      solvestat = benders->solvestat[solveidx[i]];
      assert(solvestat->idx == solveidx[i]);

      /* updating the solving statistics */
      if( updatestat )
      {
         if( subproblem == NULL )
            solvestat->avgiter = 1;
         else
            solvestat->avgiter = (SCIP_Real)(solvestat->avgiter*solvestat->ncalls + SCIPgetNLPIterations(subproblem))
               /(SCIP_Real)(solvestat->ncalls + 1);
         solvestat->ncalls++;
      }

      /* inserting the solving statistics into the priority queue */
      SCIP_CALL( SCIPpqueueInsert(benders->subprobqueue, solvestat) );
   }

   assert(SCIPpqueueNElems(benders->subprobqueue) == SCIPbendersGetNSubproblems(benders));

   return SCIP_OKAY;
}

/** Solves each of the Benders' decomposition subproblems for the given solution. All, or a fraction, of subproblems are
 *  solved before the Benders' decomposition cuts are generated.
 *  Since a convex relaxation of the subproblem could be solved to generate cuts, a parameter nverified is used to
 *  identified the number of subproblems that have been solved in their "original" form. For example, if the subproblem
 *  is a MIP, then if the LP is solved to generate cuts, this does not constitute a verification. The verification is
 *  only performed when the MIP is solved.
 */
static
SCIP_RETCODE solveBendersSubproblems(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the current solve loop */
   SCIP_Bool             checkint,           /**< are the subproblems called during a check/enforce of integer sols? */
   int*                  nverified,          /**< the number of subproblems verified in the current loop */
   int*                  solveidx,           /**< the indices of subproblems to be solved in this loop */
   int                   nsolveidx,          /**< the number of subproblems to be solved in this loop */
   SCIP_Bool**           subprobsolved,      /**< an array indicating the subproblems that were solved in this loop. */
   SCIP_BENDERSSUBSTATUS** substatus,        /**< array to store the status of the subsystem */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_Bool*            optimal,            /**< is the current solution optimal? */
   SCIP_Bool*            stopped             /**< was the solving process stopped? */
   )
{
   SCIP_Bool onlyconvexcheck;
#ifdef _OPENMP
   int numthreads;
   int maxnthreads;
#endif
   int i;
   int j;

   /* local variables for parallelisation of the solving loop */
   int locnverified = *nverified;
   SCIP_Bool locinfeasible = *infeasible;
   SCIP_Bool locoptimal = *optimal;
   SCIP_Bool locstopped = *stopped;

   SCIP_RETCODE retcode = SCIP_OKAY;

   assert(benders != NULL);
   assert(set != NULL);

   /* getting the number of threads to use when solving the subproblems. This will be either be
    * min(numthreads, maxnthreads).
    * NOTE: This may not be correct. The Benders' decomposition parallelisation should not take all minimum threads if
    * they are specified. The number of threads should be specified with the Benders' decomposition parameters.
    */
#ifdef _OPENMP
   SCIP_CALL( SCIPsetGetIntParam(set, "parallel/maxnthreads", &maxnthreads) );
   numthreads = MIN(benders->numthreads, maxnthreads);
#endif

   /* in the case of an LNS check, only the convex relaxations of the subproblems will be solved. This is a performance
    * feature, since solving the convex relaxation is typically much faster than solving the corresponding CIP. While
    * the CIP is not solved during the LNS check, the solutions are still of higher quality than when Benders' is not
    * employed.
    */
   onlyconvexcheck = SCIPbendersOnlyCheckConvexRelax(benders, SCIPsetGetSubscipsOff(set));

   SCIPsetDebugMsg(set, "Performing the subproblem solving process. Number of subproblems to check %d\n", nsolveidx);

   SCIPsetDebugMsg(set, "Benders' decomposition - solve loop %d\n", solveloop);

   if( type == SCIP_BENDERSENFOTYPE_CHECK && sol == NULL )
   {
      /* TODO: Check whether this is absolutely necessary. I think that this if statment can be removed. */
      locinfeasible = TRUE;
   }
   else
   {
      /* solving each of the subproblems for Benders' decomposition */
      /* TODO: ensure that the each of the subproblems solve and update the parameters with the correct return values
       */
#ifndef __INTEL_COMPILER
      #pragma omp parallel for num_threads(numthreads) private(i) reduction(&&:locoptimal) reduction(||:locinfeasible) reduction(+:locnverified) reduction(||:locstopped) reduction(min:retcode)
#endif
      for( j = 0; j < nsolveidx; j++ )
      {
         SCIP_Bool subinfeas = FALSE;
         SCIP_Bool convexsub;
         SCIP_Bool solvesub = TRUE;
         SCIP_Bool solved;

         i = solveidx[j];
         convexsub = SCIPbendersGetSubproblemType(benders, i) == SCIP_BENDERSSUBTYPE_CONVEXCONT;;

         /* the subproblem is initially flagged as not solved for this solving loop */
         (*subprobsolved)[i] = FALSE;

         /* setting the subsystem status to UNKNOWN at the start of each solve loop */
         (*substatus)[i] = SCIP_BENDERSSUBSTATUS_UNKNOWN;

         /* for the second solving loop, if the problem is an LP, it is not solved again. If the problem is a MIP,
          * then the subproblem objective function value is set to infinity. However, if the subproblem is proven
          * infeasible from the LP, then the IP loop is not performed.
          * If the solve loop is SCIP_BENDERSSOLVELOOP_USERCIP, then nothing is done. It is assumed that the user will
          * correctly update the objective function within the user-defined solving function.
          */
         if( solveloop == SCIP_BENDERSSOLVELOOP_CIP )
         {
            if( convexsub || (*substatus)[i] == SCIP_BENDERSSUBSTATUS_INFEAS )
               solvesub = FALSE;
            else
            {
               SCIPbendersSetSubproblemObjval(benders, i, SCIPbendersSubproblem(benders, i) != NULL ?
                  SCIPinfinity(SCIPbendersSubproblem(benders, i)) : SCIPsetInfinity(set));
            }
         }

         /* if the subproblem is independent, then it does not need to be solved. In this case, the nverified flag will
          * increase by one. When the subproblem is not independent, then it needs to be checked.
          */
         if( !subproblemIsActive(benders, i) )
         {
            /* NOTE: There is no need to update the optimal flag. This is because optimal is always TRUE until a
             * non-optimal subproblem is found.
             */
            /* if the auxiliary variable value is infinity, then the subproblem has not been solved yet. Currently the
             * subproblem statue is unknown. */
            if( SCIPsetIsInfinity(set, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i))
               || SCIPsetIsInfinity(set, -SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i))
               || SCIPsetIsInfinity(set, -SCIPbendersGetSubproblemLowerbound(benders, i)) )
            {
               SCIPbendersSetSubproblemObjval(benders, i, SCIPbendersSubproblem(benders, i) != NULL ?
                  SCIPinfinity(SCIPbendersSubproblem(benders, i)) : SCIPsetInfinity(set));

               (*substatus)[i] = SCIP_BENDERSSUBSTATUS_UNKNOWN;
               locoptimal = FALSE;

               SCIPsetDebugMsg(set, "Benders' decomposition: subproblem %d is not active, but has not been solved."
                 " setting status to UNKNOWN\n", i);
            }
            else
            {
               if( SCIPrelDiff(SCIPbendersGetSubproblemLowerbound(benders, i),
                     SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i)) < benders->solutiontol )
               {
                  SCIPbendersSetSubproblemObjval(benders, i, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i));
                  (*substatus)[i] = SCIP_BENDERSSUBSTATUS_OPTIMAL;
               }
               else
               {
                  SCIPbendersSetSubproblemObjval(benders, i, SCIPbendersGetSubproblemLowerbound(benders, i));
                  (*substatus)[i] = SCIP_BENDERSSUBSTATUS_AUXVIOL;
               }

               SCIPsetDebugMsg(set, "Benders' decomposition: subproblem %d is not active, setting status to OPTIMAL\n", i);
            }

            (*subprobsolved)[i] = TRUE;

            /* the nverified counter is only increased in the convex solving loop */
            if( solveloop == SCIP_BENDERSSOLVELOOP_CONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX )
               locnverified++;
         }
         else if( solvesub )
         {
            retcode = SCIPbendersExecSubproblemSolve(benders, set, sol, i, solveloop, FALSE, &solved, &subinfeas, type);

            /* the solution for the subproblem is only processed if the return code is SCIP_OKAY */
            if( retcode == SCIP_OKAY )
            {
#ifdef SCIP_DEBUG
               if( type == SCIP_BENDERSENFOTYPE_LP )
               {
               SCIPsetDebugMsg(set, "Enfo LP: Subproblem %d Type %d (%f < %f)\n", i,
                  SCIPbendersGetSubproblemType(benders, i), SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i),
                     SCIPbendersGetSubproblemObjval(benders, i));
               }
#endif
               (*subprobsolved)[i] = solved;

               locinfeasible = locinfeasible || subinfeas;
               if( subinfeas )
                  (*substatus)[i] = SCIP_BENDERSSUBSTATUS_INFEAS;

               /* if the subproblems are solved to check integer feasibility, then the optimality check must be performed.
                * This will only be performed if checkint is TRUE and the subproblem was solved. The subproblem may not be
                * solved if the user has defined a solving function
                */
               if( checkint && (*subprobsolved)[i] )
               {
                  /* if the subproblem is feasible, then it is necessary to update the value of the auxiliary variable to the
                   * objective function value of the subproblem.
                   */
                  if( !subinfeas )
                  {
                     SCIP_Bool subproboptimal;

                     subproboptimal = SCIPbendersSubproblemIsOptimal(benders, set, sol, i);

                     if( subproboptimal )
                        (*substatus)[i] = SCIP_BENDERSSUBSTATUS_OPTIMAL;
                     else
                        (*substatus)[i] = SCIP_BENDERSSUBSTATUS_AUXVIOL;

                     /* It is only possible to determine the optimality of a solution within a given subproblem in four
                      * different cases:
                      * i) solveloop == SCIP_BENDERSSOLVELOOP_CONVEX or USERCONVEX and the subproblem is convex.
                      * ii) solveloop == SCIP_BENDERSOLVELOOP_CONVEX  and only the convex relaxations will be checked.
                      * iii) solveloop == SCIP_BENDERSSOLVELOOP_USERCIP and the subproblem was solved, since the user has
                      * defined a solve function, it is expected that the solving is correctly executed.
                      * iv) solveloop == SCIP_BENDERSSOLVELOOP_CIP and the MIP for the subproblem has been solved.
                      */
                     if( convexsub || onlyconvexcheck
                        || solveloop == SCIP_BENDERSSOLVELOOP_CIP
                        || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP )
                        locoptimal = locoptimal && subproboptimal;

#ifdef SCIP_DEBUG
                     if( convexsub || solveloop >= SCIP_BENDERSSOLVELOOP_CIP )
                     {
                        if( subproboptimal )
                        {
                           SCIPsetDebugMsg(set, "Subproblem %d is Optimal (%f >= %f)\n", i,
                              SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i), SCIPbendersGetSubproblemObjval(benders, i));
                        }
                        else
                        {
                           SCIPsetDebugMsg(set, "Subproblem %d is NOT Optimal (%f < %f)\n", i,
                              SCIPbendersGetAuxiliaryVarVal(benders, set, sol, i), SCIPbendersGetSubproblemObjval(benders, i));
                        }
                     }
#endif

                     /* the nverified variable is only incremented when the original form of the subproblem has been solved.
                      * What is meant by "original" is that the LP relaxation of CIPs are solved to generate valid cuts. So
                      * if the subproblem is defined as a CIP, then it is only classified as checked if the CIP is solved.
                      * There are three cases where the "original" form is solved are:
                      * i) solveloop == SCIP_BENDERSSOLVELOOP_CONVEX or USERCONVEX and the subproblem is an LP
                      *    - the original form has been solved.
                      * ii) solveloop == SCIP_BENDERSSOLVELOOP_CIP or USERCIP and the CIP for the subproblem has been
                      *    solved.
                      * iii) or, only a convex check is performed.
                      */
                     if( ((solveloop == SCIP_BENDERSSOLVELOOP_CONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX)
                           && convexsub)
                        || ((solveloop == SCIP_BENDERSSOLVELOOP_CIP || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP)
                           && !convexsub)
                        || onlyconvexcheck )
                        locnverified++;
                  }
               }
            }
         }

         /* checking whether the limits have been exceeded in the master problem */
         locstopped = SCIPisStopped(set->scip);
      }
   }

   /* setting the input parameters to the local variables */
   SCIPsetDebugMsg(set, "Local variable values: nverified %d infeasible %u optimal %u stopped %u\n", locnverified,
      locinfeasible, locoptimal, locstopped);
   *nverified = locnverified;
   *infeasible = locinfeasible;
   *optimal = locoptimal;
   *stopped = locstopped;

   return retcode;
}

/** Calls the Benders' decompsition cuts for the given solve loop. There are four cases:
 *  i) solveloop == SCIP_BENDERSSOLVELOOP_CONVEX - only the LP Benders' cuts are called
 *  ii) solveloop == SCIP_BENDERSSOLVELOOP_CIP - only the CIP Benders' cuts are called
 *  iii) solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX - only the LP Benders' cuts are called
 *  iv) solveloop == SCIP_BENDERSSOLVELOOP_USERCIP - only the CIP Benders' cuts are called
 *
 *  The priority of the results are: SCIP_CONSADDED (SCIP_SEPARATED), SCIP_DIDNOTFIND, SCIP_FEASIBLE, SCIP_DIDNOTRUN. In
 *  this function, there are four levels of results that need to be assessed. These are:
 *  i) The result from the individual cut for the subproblem
 *  ii) The overall result for the subproblem from all cuts
 *  iii) the overall result for the solve loop from all cuts
 *  iv) the over all result from all solve loops.
 *  In each level, the priority of results must be adhered to.
 */
static
SCIP_RETCODE generateBendersCuts(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the current solve loop */
   SCIP_Bool             checkint,           /**< are the subproblems called during a check/enforce of integer sols? */
   SCIP_Bool*            subprobsolved,      /**< an array indicating the subproblems that were solved in this loop. */
   SCIP_BENDERSSUBSTATUS* substatus,         /**< array to store the status of the subsystem */
   int*                  solveidx,           /**< the indices of subproblems to be solved in this loop */
   int                   nsolveidx,          /**< the number of subproblems to be solved in this loop */
   int**                 mergecands,         /**< the subproblems that are merge candidates */
   int*                  npriomergecands,    /**< the number of priority merge candidates. */
   int*                  nmergecands,        /**< the number of merge candidates. */
   int*                  nsolveloops         /**< the number of solve loops, is updated w.r.t added cuts */
   )
{
   SCIP_BENDERSCUT** benderscuts;
   SCIP_RESULT solveloopresult;
   int nbenderscuts;
   SCIP_Longint addedcuts = 0;
   int i;
   int j;
   int k;
   SCIP_Bool onlyconvexcheck;

   assert(benders != NULL);
   assert(set != NULL);

   /* getting the Benders' decomposition cuts */
   benderscuts = SCIPbendersGetBenderscuts(benders);
   nbenderscuts = SCIPbendersGetNBenderscuts(benders);

   solveloopresult = SCIP_DIDNOTRUN;

   /* in the case of an LNS check, only the convex relaxations of the subproblems will be solved. This is a performance
    * feature, since solving the convex relaxation is typically much faster than solving the corresponding CIP. While
    * the CIP is not solved during the LNS check, the solutions are still of higher quality than when Benders' is not
    * employed.
    */
   onlyconvexcheck = SCIPbendersOnlyCheckConvexRelax(benders, SCIPsetGetSubscipsOff(set));

   /* It is only possible to add cuts to the problem if it has not already been solved */
   if( SCIPsetGetStage(set) < SCIP_STAGE_SOLVED
      && SCIPsetGetStage(set) != SCIP_STAGE_TRANSFORMED && SCIPsetGetStage(set) != SCIP_STAGE_PRESOLVED
      && (benders->cutcheck || type != SCIP_BENDERSENFOTYPE_CHECK) )
   {
      /* This is done in two loops. The first is by subproblem and the second is by cut type. */
      for( k = 0; k < nsolveidx; k++ )
      {
         SCIP_RESULT subprobresult;
         SCIP_Bool convexsub;

         i = solveidx[k];

         convexsub = SCIPbendersGetSubproblemType(benders, i) == SCIP_BENDERSSUBTYPE_CONVEXCONT;

         /* cuts can only be generated if the subproblem is not independent and if it has been solved. Additionally, the
          * status of the subproblem solving must not be INFEASIBLE while in a cut strengthening round.
          * The subproblem solved flag is important for the user-defined subproblem solving methods
          */
         if( subproblemIsActive(benders, i) && subprobsolved[i]
            && !(substatus[i] == SCIP_BENDERSSUBSTATUS_INFEAS && benders->strengthenround) )
         {
            subprobresult = SCIP_DIDNOTRUN;
            for( j = 0; j < nbenderscuts; j++ )
            {
               SCIP_RESULT cutresult;
               SCIP_Longint prevaddedcuts;

               assert(benderscuts[j] != NULL);

               prevaddedcuts = SCIPbenderscutGetNFound(benderscuts[j]);
               cutresult = SCIP_DIDNOTRUN;

               /* the result is updated only if a Benders' cut is generated or one was not found. However, if a cut has
                * been found in a previous iteration, then the result is returned as SCIP_CONSADDED or SCIP_SEPARATED.
                * This result is permitted because if a constraint was added, the solution that caused the error in the cut
                * generation will be cutoff from the master problem.
                */
               if( (SCIPbenderscutIsLPCut(benderscuts[j]) && (solveloop == SCIP_BENDERSSOLVELOOP_CONVEX
                        || solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX))
                  || (!SCIPbenderscutIsLPCut(benderscuts[j]) && ((solveloop == SCIP_BENDERSSOLVELOOP_CIP && !convexsub)
                        || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP)) )
                  SCIP_CALL( SCIPbenderscutExec(benderscuts[j], set, benders, sol, i, type, &cutresult) );

               addedcuts += (SCIPbenderscutGetNFound(benderscuts[j]) - prevaddedcuts);

               /* the result is updated only if a Benders' cut is generated */
               if( cutresult == SCIP_CONSADDED || cutresult == SCIP_SEPARATED )
               {
                  subprobresult = cutresult;

                  benders->ncutsfound++;

                  /* at most a single cut is generated for each subproblem */
                  break;
               }
               else
               {
                  /* checking from lowest priority result */
                  if( subprobresult == SCIP_DIDNOTRUN )
                     subprobresult = cutresult;
                  else if( subprobresult == SCIP_FEASIBLE && cutresult == SCIP_DIDNOTFIND )
                     subprobresult = cutresult;
                  /* if the subprobresult is SCIP_DIDNOTFIND, then it can't be updated. */
               }
            }

            /* the highest priority for the results is CONSADDED and SEPARATED. The solveloopresult will always be
             * updated if the subprobresult is either of these.
             */
            if( subprobresult == SCIP_CONSADDED || subprobresult == SCIP_SEPARATED )
            {
               solveloopresult = subprobresult;
            }
            else if( subprobresult == SCIP_FEASIBLE )
            {
               /* updating the solve loop result based upon the priority */
               if( solveloopresult == SCIP_DIDNOTRUN )
                  solveloopresult = subprobresult;
            }
            else if( subprobresult == SCIP_DIDNOTFIND )
            {
               /* updating the solve loop result based upon the priority */
               if( solveloopresult == SCIP_DIDNOTRUN || solveloopresult == SCIP_FEASIBLE )
                  solveloopresult = subprobresult;

               /* since a cut was not found, then merging could be useful to avoid this in subsequent iterations. The
                * candidate is labelled as a non-priority merge candidate
                */
               if( substatus[i] != SCIP_BENDERSSUBSTATUS_OPTIMAL )
               {
                  (*mergecands)[(*nmergecands)] = i;
                  (*nmergecands)++;
               }
            }
            else if( subprobresult == SCIP_DIDNOTRUN )
            {
               /* if the subproblem is infeasible and no cut generation methods were run, then the infeasibility will
                * never be resolved. As such, the subproblem will be merged into the master problem. If the subproblem
                * was not infeasible, then it is added as a possible merge candidate
                */
               if( substatus[i] == SCIP_BENDERSSUBSTATUS_INFEAS )
               {
                  (*mergecands)[(*nmergecands)] = (*mergecands)[(*npriomergecands)];
                  (*mergecands)[(*npriomergecands)] = i;
                  (*npriomergecands)++;
                  (*nmergecands)++;
               }
               else if( substatus[i] != SCIP_BENDERSSUBSTATUS_OPTIMAL )
               {
                  (*mergecands)[(*nmergecands)] = i;
                  (*nmergecands)++;
               }
            }
         }
      }
   }

   /* updating the overall result based upon the priorities */
   if( solveloopresult == SCIP_CONSADDED || solveloopresult == SCIP_SEPARATED )
   {
      (*result) = solveloopresult;
   }
   else if( solveloopresult == SCIP_FEASIBLE )
   {
      /* updating the solve loop result based upon the priority */
      if( (*result) == SCIP_DIDNOTRUN )
         (*result) = solveloopresult;
   }
   else if( solveloopresult == SCIP_DIDNOTFIND )
   {
      /* updating the solve loop result based upon the priority */
      if( (*result) == SCIP_DIDNOTRUN || (*result) == SCIP_FEASIBLE )
         (*result) = solveloopresult;
   }

   /* if no cuts were added, then the number of solve loops is increased */
   if( addedcuts == 0 && SCIPbendersGetNConvexSubproblems(benders) < SCIPbendersGetNSubproblems(benders)
      && checkint && !onlyconvexcheck )
      (*nsolveloops) = 2;

   return SCIP_OKAY;
}

/** Solves the subproblem using the current master problem solution.
 *
 *  The checkint flag indicates whether integer feasibility can be assumed. If it is not assumed, i.e. checkint ==
 *  FALSE, then only the convex relaxations of the subproblems are solved. If integer feasibility is assumed, i.e.
 *  checkint == TRUE, then the convex relaxations and the full CIP are solved to generate Benders' cuts and check
 *  solution feasibility.
 *
 *  TODO: consider allowing the possibility to pass solution information back from the subproblems instead of the scip
 *  instance. This would allow the use of different solvers for the subproblems, more importantly allowing the use of an
 *  LP solver for LP subproblems.
 */
SCIP_RETCODE SCIPbendersExec(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_RESULT*          result,             /**< result of the pricing process */
   SCIP_Bool*            infeasible,         /**< is the master problem infeasible with respect to the Benders' cuts? */
   SCIP_Bool*            auxviol,            /**< set to TRUE only if the solution is feasible but the aux vars are violated */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should the integer solution be checked by the subproblems */
   )
{
   int nsubproblems;
   int subproblemcount;
   int nsolveloops;
   int nverified;
   int nsolved;
   int* mergecands;
   int npriomergecands;
   int nmergecands;
   int* solveidx;
   int* executedidx;
   int nsolveidx;
   int nexecutedidx;
   int nfree;
   SCIP_Bool* subprobsolved;
   SCIP_BENDERSSUBSTATUS* substatus;
   SCIP_Bool optimal;
   SCIP_Bool allverified;
   SCIP_Bool success;
   SCIP_Bool stopped;
   int i;
   int l;

   success = TRUE;
   stopped = FALSE;

   SCIPsetDebugMsg(set, "Starting Benders' decomposition subproblem solving. type %d checkint %u\n", type, checkint);

#ifdef SCIP_MOREDEBUG
   SCIP_CALL( SCIPprintSol(set->scip, sol, NULL, FALSE) );
#endif

   /* start timing */
   SCIPclockStart(benders->bendersclock, set);

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   (*auxviol) = FALSE;
   (*infeasible) = FALSE;

   /* It is assumed that the problem is optimal, until a subproblem is found not to be optimal. However, not all
    * subproblems could be checked in each iteration. As such, it is not possible to state that the problem is optimal
    * if not all subproblems are checked. Situations where this may occur is when a subproblem is a MIP and only the LP
    * is solved. Also, in a distributed computation, then it may be advantageous to only solve some subproblems before
    * resolving the master problem. As such, for a problem to be optimal, then (optimal && allverified) == TRUE
    */
   optimal = TRUE;
   nverified = 0;
   nsolved = 0;

   assert(benders != NULL);
   assert(result != NULL);
   assert(infeasible != NULL);
   assert(auxviol != NULL);

   /* if the Benders' decomposition is called from a sub-SCIP and the sub-SCIPs have been deactivated, then it is
    * assumed that this is an LNS heuristic. As such, the check is not performed and the solution is assumed to be
    * feasible
    */
   if( benders->iscopy && set->subscipsoff
      && (!benders->lnscheck
         || (benders->lnsmaxdepth > -1 && SCIPgetDepth(benders->sourcescip) >= benders->lnsmaxdepth)
         || (benders->lnsmaxcalls > -1 && SCIPbendersGetNCalls(benders) >= benders->lnsmaxcalls)
         || (type != SCIP_BENDERSENFOTYPE_CHECK && SCIPgetDepth(set->scip) == 0 && benders->lnsmaxcallsroot > -1
            && SCIPbendersGetNCalls(benders) >= benders->lnsmaxcallsroot)) )
   {
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* it is not necessary to check all primal solutions by solving the Benders' decomposition subproblems.
    * Only the improving solutions are checked to improve efficiency of the algorithm.
    * If the solution is non-improving, the result FEASIBLE is returned. While this may be incorrect w.r.t to the
    * Benders' subproblems, this solution will never be the optimal solution. A non-improving solution may be used
    * within LNS primal heuristics. If this occurs, the improving solution, if found, will be checked by the solving
    * the Benders' decomposition subproblems.
    * TODO: Add a parameter to control this behaviour.
    */
   if( checkint && SCIPsetIsLE(set, SCIPgetPrimalbound(set->scip)*(int)SCIPgetObjsense(set->scip),
         SCIPgetSolOrigObj(set->scip, sol)*(int)SCIPgetObjsense(set->scip)) )
   {
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* if the enforcement type is SCIP_BENDERSENFOTYPE_LP and the LP is currently unbounded. This could mean that there
    * is no lower bound on the auxiliary variables. In this case, we try to update the lower bound for the auxiliary
    * variables.
    */
   if( type == SCIP_BENDERSENFOTYPE_LP && SCIPgetLPSolstat(set->scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY
      && benders->updateauxvarbound )
   {
      SCIP_CALL( updateAuxiliaryVarLowerbound(benders, set, result) );

      /* the auxiliary variable bound will only be updated once. */
      benders->updateauxvarbound = FALSE;
   }

   /* sets the stored objective function values of the subproblems to infinity */
   resetSubproblemObjectiveValue(benders, set);

   *result = SCIP_DIDNOTRUN;

   if( benders->benderspresubsolve != NULL && !benders->strengthenround )
   {
      SCIP_Bool skipsolve;

      skipsolve = FALSE;
      SCIP_CALL( benders->benderspresubsolve(set->scip, benders, sol, type, checkint, infeasible, auxviol, &skipsolve,
            result) );

      /* evaluate result */
      if( (*result) != SCIP_DIDNOTRUN
         && (*result) != SCIP_FEASIBLE
         && (*result) != SCIP_INFEASIBLE
         && (*result) != SCIP_CONSADDED
         && (*result) != SCIP_SEPARATED )
      {
         SCIPerrorMessage("the user-defined pre subproblem solving method for the Benders' decomposition <%s> returned "
            "invalid result <%d>\n", benders->name, *result);
         return SCIP_INVALIDRESULT;
      }

      /* if the solve must be skipped, then the solving loop is exited and the user defined result is returned */
      if( skipsolve )
      {
         SCIPsetDebugMsg(set, "skipping the subproblem solving for Benders' decomposition <%s>. "
            "returning result <%d>\n", benders->name, *result);
         return SCIP_OKAY;
      }
   }

   /* the cut strengthening is performed before the regular subproblem solve is called. To avoid recursion, the flag
    * strengthenround is set to TRUE when the cut strengthening is performed. The cut strengthening is not performed as
    * part of the large neighbourhood Benders' search.
    *
    * NOTE: cut strengthening is only applied for fractional solutions and integer solutions if there are no CIP
    * subproblems.
    */
   if( benders->strengthenenabled && !benders->strengthenround && !benders->iscopy
      && (!checkint || SCIPbendersGetNConvexSubproblems(benders) == SCIPbendersGetNSubproblems(benders)) )
   {
      SCIP_Bool skipsolve;

      benders->strengthenround = TRUE;
      /* if the user has not requested the solve to be skipped, then the cut strengthening is performed */
      SCIP_CALL( performInteriorSolCutStrengthening(benders, set, sol, type, checkint, FALSE, infeasible, auxviol,
            &skipsolve, result) );
      benders->strengthenround = FALSE;

      /* if the solve must be skipped, then the solving loop is exited and the user defined result is returned */
      if( skipsolve )
      {
         SCIPsetDebugMsg(set, "skipping the subproblem solving because cut strengthening found a cut "
            "for Benders' decomposition <%s>. Returning result <%d>\n", benders->name, *result);
         return SCIP_OKAY;
      }

      /* the result flag need to be reset to DIDNOTRUN for the main subproblem solve */
      (*result) = SCIP_DIDNOTRUN;
   }

   /* allocating memory for the infeasible subproblem array */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(set->scip, &subprobsolved, nsubproblems) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(set->scip, &substatus, nsubproblems) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(set->scip, &mergecands, nsubproblems) );
   npriomergecands = 0;
   nmergecands = 0;

   /* allocating the memory for the subproblem solving and cut generation indices */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(set->scip, &solveidx, nsubproblems) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(set->scip, &executedidx, nsubproblems) );
   nsolveidx = 0;
   nexecutedidx = 0;

   /* only a subset of the subproblems are initially solved. Both solving loops are executed for the subproblems to
    * check whether any cuts are generated. If a cut is generated, then no further subproblems are solved. If a cut is
    * not generated, then an additional set of subproblems are solved.
    */
   while( nsolved < nsubproblems )
   {
      /* getting the indices for the subproblems that will be solved */
      createSolveSubproblemIndexList(benders, set, type, &solveidx, &nsolveidx);

      /* by default the number of solve loops is 1. This is the case if all subproblems are LP or the user has defined a
       * benderssolvesub callback. If there is a subproblem that is not an LP, then 2 solve loops are performed. The first
       * loop is the LP solving loop, the second solves the subproblem to integer optimality.
       */
      nsolveloops = 1;

      for( l = 0; l < nsolveloops; l++ )
      {
         SCIP_BENDERSSOLVELOOP solveloop;    /* identifies what problem type is solve in this solve loop */

         /* if either benderssolvesubconvex or benderssolvesub are implemented, then the user callbacks are invoked */
         if( benders->benderssolvesubconvex != NULL || benders->benderssolvesub != NULL )
         {
            if( l == 0 )
               solveloop = SCIP_BENDERSSOLVELOOP_USERCONVEX;
            else
               solveloop = SCIP_BENDERSSOLVELOOP_USERCIP;
         }
         else
            solveloop = (SCIP_BENDERSSOLVELOOP) l;

         /* solving the subproblems for this round of enforcement/checking. */
         SCIP_CALL( solveBendersSubproblems(benders, set, sol, type, solveloop, checkint, &nverified,
               solveidx, nsolveidx, &subprobsolved, &substatus, infeasible, &optimal, &stopped) );

         /* if the solving has been stopped, then the subproblem solving and cut generation must terminate */
         if( stopped )
            break;

         /* Generating cuts for the subproblems. Cuts are only generated when the solution is from primal heuristics,
          * relaxations or the LP
          */
         if( type != SCIP_BENDERSENFOTYPE_PSEUDO )
         {
            SCIP_CALL( generateBendersCuts(benders, set, sol, result, type, solveloop, checkint, subprobsolved,
                  substatus, solveidx, nsolveidx, &mergecands, &npriomergecands, &nmergecands, &nsolveloops) );
         }
         else
         {
            /* The first solving loop solves the convex subproblems and the convex relaxations of the CIP subproblems. The
             * second solving loop solves the CIP subproblems. The second solving loop is only called if the integer
             * feasibility is being checked and if the convex subproblems and convex relaxations are not infeasible.
             */
            if( !(*infeasible) && checkint && !SCIPbendersOnlyCheckConvexRelax(benders, SCIPsetGetSubscipsOff(set))
               && SCIPbendersGetNConvexSubproblems(benders) < SCIPbendersGetNSubproblems(benders))
               nsolveloops = 2;
         }
      }

      nsolved += nsolveidx;

      /* storing the indices of the subproblems for which the solving loop was executed */
      for( i = 0; i < nsolveidx; i++ )
         executedidx[nexecutedidx++] = solveidx[i];

      /* if the result is CONSADDED or SEPARATED, then a cut is generated and no further subproblem processing is
       * required
       */
      if( (*result) == SCIP_CONSADDED || (*result) == SCIP_SEPARATED )
         break;
   }

   /* inserting the subproblems into the priority queue for the next solve call */
   SCIP_CALL( updateSubproblemStatQueue(benders, executedidx, nexecutedidx, TRUE) );

   if( stopped )
      goto TERMINATE;

   allverified = (nverified == nsubproblems);

   SCIPsetDebugMsg(set, "End Benders' decomposition subproblem solve. result %d infeasible %u auxviol %u nverified %d\n",
      *result, *infeasible, *auxviol, nverified);

#ifdef SCIP_DEBUG
   if( (*result) == SCIP_CONSADDED )
   {
      SCIPsetDebugMsg(set, "Benders' decomposition: Cut added\n");
   }
#endif

   /* if the number of checked pseudo solutions exceeds a set limit, then all subproblems are passed as merge
    * candidates. Currently, merging subproblems into the master problem is the only method for resolving numerical
    * troubles.
    *
    * We are only interested in the pseudo solutions that have been checked completely for integrality. This is
    * identified by checkint == TRUE. This means that the Benders' decomposition constraint is one of the last
    * constraint handlers that must resolve the infeasibility. If the Benders' decomposition framework can't resolve the
    * infeasibility, then this will result in an error.
    */
   if( type == SCIP_BENDERSENFOTYPE_PSEUDO && checkint )
   {
      benders->npseudosols++;

      if( benders->npseudosols > BENDERS_MAXPSEUDOSOLS )
      {
         /* if a priority merge candidate already exists, then no other merge candidates need to be added.*/
         if( npriomergecands == 0 )
         {
            /* all subproblems are added to the merge candidate list. The first active subproblem is added as a
             * priority merge candidate
             */
            nmergecands = 0;
            npriomergecands = 1;
            for( i = 0; i < nsubproblems; i++ )
            {
               /* only active subproblems are added to the merge candidate list */
               if( subproblemIsActive(benders, i) )
               {
                  mergecands[nmergecands] = i;
                  nmergecands++;
               }
            }

            SCIPverbMessage(set->scip, SCIP_VERBLEVEL_HIGH, NULL, "   The number of checked pseudo solutions exceeds the "
              "limit of %d. All active subproblems are merge candidates, with subproblem %d a priority candidate.\n",
              BENDERS_MAXPSEUDOSOLS, mergecands[0]);
         }
      }
   }
   else
      benders->npseudosols = 0;

   /* if the result is SCIP_DIDNOTFIND, then there was a error in generating cuts in all subproblems that are not
    * optimal. This result does not cutoff any solution, so the Benders' decomposition algorithm will fail.
    *
    * It could happen that the cut strengthening approach causes an error the cut generation. In this case, an error
    * should not be thrown. So, this check will be skipped when in a strengthening round.
    * TODO: Work out a way to ensure Benders' decomposition does not terminate due to a SCIP_DIDNOTFIND result.
    */
   if( (*result) == SCIP_DIDNOTFIND && !benders->strengthenround )
   {
      if( type == SCIP_BENDERSENFOTYPE_PSEUDO )
         (*result) = SCIP_SOLVELP;
      else
         (*result) = SCIP_INFEASIBLE;

      SCIPerrorMessage("An error was found when generating cuts for non-optimal subproblems of Benders' "
         "decomposition <%s>. Consider merging the infeasible subproblems into the master problem.\n", SCIPbendersGetName(benders));

      /* since no other cuts are generated, then this error will result in a crash. It is possible to avoid the error,
       * by merging the affected subproblem into the master problem.
       *
       * NOTE: If the error occurs while checking solutions, i.e. SCIP_BENDERSENFOTYPE_CHECK, then it is valid to set
       * the result to SCIP_INFEASIBLE and the success flag to TRUE
       */
      if( type != SCIP_BENDERSENFOTYPE_CHECK )
         success = FALSE;

      goto POSTSOLVE;
   }

   if( type == SCIP_BENDERSENFOTYPE_PSEUDO )
   {
      if( (*infeasible) || !allverified )
         (*result) = SCIP_SOLVELP;
      else
      {
         (*result) = SCIP_FEASIBLE;

         /* if the subproblems are not infeasible, but they are also not optimal. This means that there is a violation
          * in the auxiliary variable values. In this case, a feasible result is returned with the auxviol flag set to
          * TRUE.
          */
         (*auxviol) = !optimal;
      }
   }
   else if( checkint && (type == SCIP_BENDERSENFOTYPE_CHECK
         || ((*result) != SCIP_CONSADDED && (*result) != SCIP_SEPARATED)) )
   {
      /* if the subproblems are being solved as part of conscheck, then the results flag must be returned after the solving
       * has completed.
       */
      if( (*infeasible) || !allverified )
         (*result) = SCIP_INFEASIBLE;
      else
      {
         (*result) = SCIP_FEASIBLE;

         /* if the subproblems are not infeasible, but they are also not optimal. This means that there is a violation
          * in the auxiliary variable values. In this case, a feasible result is returned with the auxviol flag set to
          * TRUE.
          */
         (*auxviol) = !optimal;
      }
   }

POSTSOLVE:
   /* calling the post-solve call back for the Benders' decomposition algorithm. This allows the user to work directly
    * with the solved subproblems and the master problem */
   if( benders->benderspostsolve != NULL )
   {
      SCIP_Bool merged;

      merged = FALSE;

      SCIP_CALL( benders->benderspostsolve(set->scip, benders, sol, type, mergecands, npriomergecands, nmergecands,
            checkint, (*infeasible), &merged) );

      if( merged )
      {
         (*result) = SCIP_CONSADDED;

         /* since subproblems have been merged, then constraints have been added. This could resolve the unresolved
          * infeasibility, so the error has been corrected.
          */
         success = TRUE;
      }
      else if( !success )
      {
         SCIPerrorMessage("An error occurred during Benders' decomposition cut generations and no merging had been "
            "performed. It is not possible to continue solving the problem by Benders' decomposition\n");
      }
   }

TERMINATE:
   /* if the solving process has stopped, then all subproblems need to be freed */
   if( stopped )
      nfree = nsubproblems;
   else
      nfree = nexecutedidx;

   /* freeing the subproblems after the cuts are generated */
   subproblemcount = 0;
   while( subproblemcount < nfree )
   {
      int subidx;

      if( stopped )
         subidx = subproblemcount;
      else
         subidx = executedidx[subproblemcount];

      SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, subidx) );

      subproblemcount++;
   }

#ifndef NDEBUG
   for( i = 0; i < nsubproblems; i++ )
      assert(SCIPbendersSubproblem(benders, i) == NULL
         || SCIPgetStage(SCIPbendersSubproblem(benders, i)) < SCIP_STAGE_TRANSFORMED
         || !SCIPinProbing(SCIPbendersSubproblem(benders, i))
         || !subproblemIsActive(benders, i));
#endif

   /* increment the number of calls to the Benders' decomposition subproblem solve */
   benders->ncalls++;

   SCIPsetDebugMsg(set, "End Benders' decomposition execution method. result %d infeasible %u auxviol %u\n", *result,
      *infeasible, *auxviol);

   /* end timing */
   SCIPclockStop(benders->bendersclock, set);

   /* freeing memory */
   SCIPfreeBlockMemoryArray(set->scip, &executedidx, nsubproblems);
   SCIPfreeBlockMemoryArray(set->scip, &solveidx, nsubproblems);
   SCIPfreeBlockMemoryArray(set->scip, &mergecands, nsubproblems);
   SCIPfreeBlockMemoryArray(set->scip, &substatus, nsubproblems);
   SCIPfreeBlockMemoryArray(set->scip, &subprobsolved, nsubproblems);

   /* if there was an error in generating cuts and merging was not performed, then the solution is perturbed in an
    * attempt to generate a cut and correct the infeasibility
    */
   if( !success && !stopped )
   {
      SCIP_Bool skipsolve;
      SCIP_RESULT perturbresult;

      skipsolve = FALSE;

      benders->strengthenround = TRUE;
      /* if the user has not requested the solve to be skipped, then the cut strengthening is performed */
      SCIP_CALL( performInteriorSolCutStrengthening(benders, set, sol, type, checkint, TRUE, infeasible, auxviol,
            &skipsolve, &perturbresult) );
      benders->strengthenround = FALSE;

      if( perturbresult == SCIP_CONSADDED || perturbresult == SCIP_SEPARATED )
         (*result) = perturbresult;

      success = skipsolve;
   }

   /* if the Benders' decomposition subproblem check stopped, then we don't have a valid result. In this case, the
    * safest thing to do is report INFEASIBLE.
    */
   if( stopped )
      (*result) = SCIP_INFEASIBLE;

   /* if the subproblem verification identifies the solution as feasible, then a check whether slack variables have been
    * used is necessary. If any slack variables are non-zero, then the solution is reverified after the objective
    * coefficient for the slack variables is increased.
    */
   if( (*result) == SCIP_FEASIBLE )
   {
      SCIP_Bool activeslack;

      SCIP_CALL( SCIPbendersSolSlackVarsActive(benders, &activeslack) );
      SCIPsetDebugMsg(set, "Type: %d Active slack: %u Feasibility Phase: %u\n", type, activeslack,
         benders->feasibilityphase);
      if( activeslack )
      {
         if( type == SCIP_BENDERSENFOTYPE_CHECK )
            (*result) = SCIP_INFEASIBLE;
         else
         {
            /* increasing the value of the slack variable by a factor of 10 */
            benders->slackvarcoef *= 10.0;

            if( benders->slackvarcoef <= benders->maxslackvarcoef )
            {
               SCIPmessagePrintVerbInfo(SCIPgetMessagehdlr(set->scip), set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "Increasing the slack variable coefficient to %g.\n", benders->slackvarcoef);
            }
            else
            {
               SCIPmessagePrintVerbInfo(SCIPgetMessagehdlr(set->scip), set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "Fixing the slack variables to zero.\n");
            }

            /* resolving the subproblems with an increased slack variable */
            SCIP_CALL( SCIPsolveBendersSubproblems(set->scip, benders, sol, result, infeasible, auxviol, type, checkint) );
         }
      }
      else if( benders->feasibilityphase )
      {
         if( type != SCIP_BENDERSENFOTYPE_CHECK )
         {
            /* disabling the feasibility phase */
            benders->feasibilityphase = FALSE;

            /* resolving the subproblems with the slack variables set to zero */
            SCIP_CALL( SCIPsolveBendersSubproblems(set->scip, benders, sol, result, infeasible, auxviol, type, checkint) );
         }
      }
   }

   if( !success )
      return SCIP_ERROR;
   else
      return SCIP_OKAY;
}

/** solves the user-defined subproblem solving function */
static
SCIP_RETCODE executeUserDefinedSolvesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the solve loop iteration. The first iter is for LP, the second for IP */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_Real*            objective,          /**< the objective function value of the subproblem */
   SCIP_RESULT*          result              /**< the result from solving the subproblem */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);
   assert(benders->benderssolvesubconvex != NULL || benders->benderssolvesub != NULL);

   assert(solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP);

   (*objective) = -SCIPsetInfinity(set);

   /* calls the user defined subproblem solving method. Only the convex relaxations are solved during the Large
    * Neighbourhood Benders' Search. */
   if( solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX )
   {
      if( benders->benderssolvesubconvex != NULL )
      {
         SCIP_CALL( benders->benderssolvesubconvex(set->scip, benders, sol, probnumber,
               SCIPbendersOnlyCheckConvexRelax(benders, SCIPsetGetSubscipsOff(set)), objective, result) );
      }
      else
         (*result) = SCIP_DIDNOTRUN;
   }
   else if( solveloop == SCIP_BENDERSSOLVELOOP_USERCIP )
   {
      if( benders->benderssolvesub != NULL )
      {
         SCIP_CALL( benders->benderssolvesub(set->scip, benders, sol, probnumber, objective, result) );
      }
      else
         (*result) = SCIP_DIDNOTRUN;
   }

   /* evaluate result */
   if( (*result) != SCIP_DIDNOTRUN
      && (*result) != SCIP_FEASIBLE
      && (*result) != SCIP_INFEASIBLE
      && (*result) != SCIP_UNBOUNDED )
   {
      SCIPerrorMessage("the user-defined solving method for the Benders' decomposition <%s> returned invalid result <%d>\n",
         benders->name, *result);
      return SCIP_INVALIDRESULT;
   }

   if( (*result) == SCIP_INFEASIBLE )
      (*infeasible) = TRUE;

   if( (*result) == SCIP_FEASIBLE
      && (SCIPsetIsInfinity(set, -(*objective)) || SCIPsetIsInfinity(set, (*objective))) )
   {
      SCIPerrorMessage("the user-defined solving method for the Benders' decomposition <%s> returned objective value %g\n",
         benders->name, (*objective));
      return SCIP_ERROR;
   }

   /* if the result is SCIP_DIDNOTFIND, then an error is returned and SCIP will terminate. */
   if( (*result) == SCIP_DIDNOTFIND )
      return SCIP_ERROR;
   else
      return SCIP_OKAY;
}

/** executes the subproblem solving process */
SCIP_RETCODE SCIPbendersExecSubproblemSolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSSOLVELOOP solveloop,          /**< the solve loop iteration. The first iter is for LP, the second for IP */
   SCIP_Bool             enhancement,        /**< is the solve performed as part of and enhancement? */
   SCIP_Bool*            solved,             /**< flag to indicate whether the subproblem was solved */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   )
{  /*lint --e{715}*/
   SCIP* subproblem;
   SCIP_RESULT result;
   SCIP_Real objective;
   SCIP_STATUS solvestatus = SCIP_STATUS_UNKNOWN;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   SCIPsetDebugMsg(set, "Benders' decomposition: solving subproblem %d\n", probnumber);

   result = SCIP_DIDNOTRUN;
   objective = SCIPsetInfinity(set);

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   if( subproblem == NULL && (benders->benderssolvesubconvex == NULL || benders->benderssolvesub == NULL) )
   {
      SCIPerrorMessage("The subproblem %d is set to NULL, but both bendersSolvesubconvex%s and bendersSolvesub%s "
         "are not defined.\n", probnumber, benders->name, benders->name);
      SCIPABORT();
      return SCIP_ERROR;
   }

   /* initially setting the solved flag to FALSE */
   (*solved) = FALSE;

   /* if the subproblem solve callback is implemented, then that is used instead of the default setup */
   if( solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP )
   {
      /* calls the user defined subproblem solving method. Only the convex relaxations are solved during the Large
       * Neighbourhood Benders' Search. */
      SCIP_CALL( executeUserDefinedSolvesub(benders, set, sol, probnumber, solveloop, infeasible, &objective, &result) );

      /* if the result is DIDNOTRUN, then the subproblem was not solved */
      (*solved) = (result != SCIP_DIDNOTRUN);
   }
   else if( subproblem != NULL )
   {
      /* setting up the subproblem */
      if( solveloop == SCIP_BENDERSSOLVELOOP_CONVEX )
      {
         SCIP_CALL( SCIPbendersSetupSubproblem(benders, set, sol, probnumber, type) );

         /* if the limits of the master problem were hit during the setup process, then the subproblem will not have
          * been setup. In this case, the solving function must be exited.
          */
         if( !SCIPbendersSubproblemIsSetup(benders, probnumber) )
         {
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
            (*solved) = FALSE;
            return SCIP_OKAY;
         }
      }
      else
      {
         SCIP_CALL( updateEventhdlrUpperbound(benders, probnumber, SCIPbendersGetAuxiliaryVarVal(benders, set, sol, probnumber)) );
      }

      /* solving the subproblem
       * the LP of the subproblem is solved in the first solveloop.
       * In the second solve loop, the MIP problem is solved */
      if( solveloop == SCIP_BENDERSSOLVELOOP_CONVEX
         || SCIPbendersGetSubproblemType(benders, probnumber) == SCIP_BENDERSSUBTYPE_CONVEXCONT )
      {
         SCIP_CALL( SCIPbendersSolveSubproblemLP(set->scip, benders, probnumber, &solvestatus, &objective) );

         /* if the (N)LP was solved without error, then the subproblem is labelled as solved */
         if( solvestatus == SCIP_STATUS_OPTIMAL || solvestatus == SCIP_STATUS_INFEASIBLE )
            (*solved) = TRUE;

         if( solvestatus == SCIP_STATUS_INFEASIBLE )
            (*infeasible) = TRUE;
      }
      else
      {
         SCIP_SOL* bestsol;

         SCIP_CALL( SCIPbendersSolveSubproblemCIP(set->scip, benders, probnumber, &solvestatus, FALSE) );

         if( solvestatus == SCIP_STATUS_INFEASIBLE )
            (*infeasible) = TRUE;

         /* if the generic subproblem solving methods are used, then the CIP subproblems are always solved. */
         (*solved) = TRUE;

         bestsol = SCIPgetBestSol(subproblem);
         if( bestsol != NULL )
            objective = SCIPgetSolOrigObj(subproblem, bestsol)*(int)SCIPgetObjsense(set->scip);
         else
            objective = SCIPsetInfinity(set);
      }
   }
   else
   {
      SCIPABORT();
   }

   if( !enhancement )
   {
      /* The following handles the cases when the subproblem is OPTIMAL, INFEASIBLE and UNBOUNDED.
       * If a subproblem is unbounded, then the auxiliary variables are set to -infinity and the unbounded flag is
       * returned as TRUE. No cut will be generated, but the result will be set to SCIP_FEASIBLE.
       */
      if( solveloop == SCIP_BENDERSSOLVELOOP_CONVEX || solveloop == SCIP_BENDERSSOLVELOOP_CIP )
      {
         /* TODO: Consider whether other solutions status should be handled */
         if( solvestatus == SCIP_STATUS_OPTIMAL )
            SCIPbendersSetSubproblemObjval(benders, probnumber, objective);
         else if( solvestatus == SCIP_STATUS_INFEASIBLE )
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
         else if( solvestatus == SCIP_STATUS_USERINTERRUPT || solvestatus == SCIP_STATUS_BESTSOLLIMIT )
            SCIPbendersSetSubproblemObjval(benders, probnumber, objective);
         else if( solvestatus == SCIP_STATUS_MEMLIMIT || solvestatus == SCIP_STATUS_TIMELIMIT
            || solvestatus == SCIP_STATUS_UNKNOWN )
         {
            SCIPverbMessage(set->scip, SCIP_VERBLEVEL_FULL, NULL, "   Benders' decomposition: Error solving "
               "subproblem %d. No cut will be generated for this subproblem.\n", probnumber);
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
         }
         else if( solvestatus == SCIP_STATUS_UNBOUNDED )
         {
            SCIPerrorMessage("The Benders' decomposition subproblem %d is unbounded. This should not happen.\n",
               probnumber);
            SCIPABORT();
         }
         else
         {
            SCIPerrorMessage("Invalid status returned from solving Benders' decomposition subproblem %d. Solution status: %d\n",
               probnumber, solvestatus);
            SCIPABORT();
         }
      }
      else
      {
         assert(solveloop == SCIP_BENDERSSOLVELOOP_USERCONVEX || solveloop == SCIP_BENDERSSOLVELOOP_USERCIP);
         if( result == SCIP_FEASIBLE )
            SCIPbendersSetSubproblemObjval(benders, probnumber, objective);
         else if( result == SCIP_INFEASIBLE )
            SCIPbendersSetSubproblemObjval(benders, probnumber, SCIPsetInfinity(set));
         else if( result == SCIP_UNBOUNDED )
         {
            SCIPerrorMessage("The Benders' decomposition subproblem %d is unbounded. This should not happen.\n",
               probnumber);
            SCIPABORT();
         }
         else if( result != SCIP_DIDNOTRUN )
         {
            SCIPerrorMessage("Invalid result <%d> from user-defined subproblem solving method. This should not happen.\n",
               result);
         }
      }
   }

   return SCIP_OKAY;
}

/** sets up the subproblem using the solution to the master problem  */
SCIP_RETCODE SCIPbendersSetupSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSENFOTYPE  type                /**< the enforcement type calling this function */
   )
{
   SCIP* subproblem;
   SCIP_VAR** vars;
   SCIP_VAR* mastervar;
   SCIP_Real solval;
   int nvars;
   int i;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   /* the subproblem setup can only be performed if the subproblem is not NULL */
   if( subproblem == NULL )
   {
      SCIPerrorMessage("The subproblem %d is NULL. Thus, the subproblem setup must be performed manually in either "
         "bendersSolvesubconvex%s or bendersSolvesub%s.\n", probnumber, benders->name, benders->name);
      return SCIP_ERROR;
   }
   assert(subproblem != NULL);

   /* changing all of the master problem variable to continuous. */
   SCIP_CALL( SCIPbendersChgMastervarsToCont(benders, set, probnumber) );

   /* if the Benders' decomposition subproblem is convex and has continuous variables, then probing mode
    * must be started.
    * If the subproblem contains non-convex constraints or discrete variables, then the problem must be initialised,
    * and then put into SCIP_STAGE_SOLVING to be able to change the variable bounds. The probing mode is entered once
    * the variable bounds are set.
    * In the latter case, the transformed problem is freed after each subproblem solve round. */
   if( SCIPbendersGetSubproblemType(benders, probnumber) == SCIP_BENDERSSUBTYPE_CONVEXCONT )
   {
      SCIP_CALL( SCIPstartProbing(subproblem) );
   }
   else
   {
      SCIP_Bool success;

      SCIP_CALL( initialiseSubproblem(benders, set, probnumber, &success) );

      if( !success )
      {
         /* set the flag to indicate that the subproblems have been set up */
         SCIPbendersSetSubproblemIsSetup(benders, probnumber, FALSE);

         return SCIP_OKAY;
      }
   }

   vars = SCIPgetVars(subproblem);
   nvars = SCIPgetNVars(subproblem);

   /* looping over all variables in the subproblem to find those corresponding to the master problem variables. */
   /* TODO: It should be possible to store the pointers to the master variables to speed up the subproblem setup */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPbendersGetVar(benders, set, vars[i], &mastervar, -1) );

      if( mastervar != NULL )
      {
         /* It is possible due to numerics that the solution value exceeds the upper or lower bounds. When this
          * happens, it causes an error in the LP solver as a result of inconsistent bounds. So the following statements
          * are used to ensure that the bounds are not exceeded when applying the fixings for the Benders'
          * decomposition subproblems
          */
         solval = SCIPgetSolVal(set->scip, sol, mastervar);
         if( !SCIPisLT(set->scip, solval, SCIPvarGetUbLocal(vars[i])) )
            solval = SCIPvarGetUbLocal(vars[i]);
         else if( !SCIPisGT(set->scip, solval, SCIPvarGetLbLocal(vars[i])) )
            solval = SCIPvarGetLbLocal(vars[i]);

         /* fixing the variable in the subproblem */
         if( !SCIPisEQ(subproblem, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) )
         {
            if( SCIPisGT(subproblem, solval, SCIPvarGetLbLocal(vars[i])) )
            {
               SCIP_CALL( SCIPchgVarLb(subproblem, vars[i], solval) );
            }
            if( SCIPisLT(subproblem, solval, SCIPvarGetUbLocal(vars[i])) )
            {
               SCIP_CALL( SCIPchgVarUb(subproblem, vars[i], solval) );
            }
         }

         assert(SCIPisEQ(subproblem, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])));
      }
      else if( strstr(SCIPvarGetName(vars[i]), SLACKVAR_NAME) != NULL )
      {
         /* if the slack variables have been added to help improve feasibility, then they remain unfixed with a large
          * objective coefficient. Once the root node has been solved to optimality, then the slack variables are
          * fixed to zero.
          */
         if( benders->feasibilityphase && SCIPgetDepth(set->scip) == 0 && type != SCIP_BENDERSENFOTYPE_CHECK )
         {
            /* The coefficient update or variable fixing can only be performed if the subproblem is in probing mode.
             * If the slack var coef gets very large, then we fix the slack variable to 0 instead.
             */
            if( SCIPinProbing(subproblem) )
            {
               if( benders->slackvarcoef <= benders->maxslackvarcoef )
               {
                  SCIP_CALL( SCIPchgVarObjProbing(subproblem, vars[i], benders->slackvarcoef) );
               }
               else
               {
                  SCIP_CALL( SCIPchgVarUbProbing(subproblem, vars[i], 0.0) );
               }
            }
         }
         else
         {
            /* if the subproblem is non-linear and convex, then slack variables have been added to the subproblem. These
             * need to be fixed to zero when first solving the subproblem. However, if the slack variables have been added
             * by setting the execfeasphase runtime parameter, then they must not get fixed to zero
             */
            assert( !SCIPisEQ(subproblem, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) );
            assert( SCIPisZero(subproblem, SCIPvarGetLbLocal(vars[i])) );

            if( SCIPisLT(subproblem, 0.0, SCIPvarGetUbLocal(vars[i])) )
            {
               SCIP_CALL( SCIPchgVarUb(subproblem, vars[i], 0.0) );
            }
         }
      }
   }

   /* if the subproblem contain non-convex constraints or discrete variables, then the probing mode is entered after
    * setting up the subproblem
    */
   if( SCIPbendersGetSubproblemType(benders, probnumber) != SCIP_BENDERSSUBTYPE_CONVEXCONT )
   {
      SCIP_CALL( SCIPstartProbing(subproblem) );
   }

   /* set the flag to indicate that the subproblems have been set up */
   SCIPbendersSetSubproblemIsSetup(benders, probnumber, TRUE);

   return SCIP_OKAY;
}

/** Solve a Benders' decomposition subproblems. This will either call the user defined method or the generic solving
 *  methods. If the generic method is called, then the subproblem must be set up before calling this method. */
SCIP_RETCODE SCIPbendersSolveSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution, can be NULL */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool*            infeasible,         /**< returns whether the current subproblem is infeasible */
   SCIP_Bool             solvecip,           /**< directly solve the CIP subproblem */
   SCIP_Real*            objective           /**< the objective function value of the subproblem, can be NULL */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   (*infeasible) = FALSE;

   /* the subproblem must be set up before this function is called. */
   if( SCIPbendersSubproblem(benders, probnumber) != NULL && !SCIPbendersSubproblemIsSetup(benders, probnumber)
      && !SCIPbendersSubproblemIsIndependent(benders, probnumber) )
   {
      SCIPerrorMessage("Benders' decomposition subproblem %d must be set up before calling SCIPbendersSolveSubproblem(). Call SCIPsetupSubproblem() first.\n", probnumber);
      return SCIP_ERROR;
   }

   /* if the subproblem solve callback is implemented, then that is used instead of the default setup */
   if( benders->benderssolvesubconvex != NULL ||  benders->benderssolvesub != NULL)
   {
      SCIP_BENDERSSOLVELOOP solveloop;
      SCIP_RESULT result;
      SCIP_Real subobj;

      if( solvecip )
         solveloop = SCIP_BENDERSSOLVELOOP_USERCIP;
      else
         solveloop = SCIP_BENDERSSOLVELOOP_USERCONVEX;

      SCIP_CALL( executeUserDefinedSolvesub(benders, set, sol, probnumber, solveloop, infeasible, &subobj, &result) );

      if( objective != NULL )
         (*objective) = subobj;
   }
   else
   {
      SCIP* subproblem;

      subproblem = SCIPbendersSubproblem(benders, probnumber);
      assert(subproblem != NULL);

      /* solving the subproblem */
      if( solvecip && SCIPbendersGetSubproblemType(benders, probnumber) != SCIP_BENDERSSUBTYPE_CONVEXCONT )
      {
         SCIP_STATUS solvestatus;

         SCIP_CALL( SCIPbendersSolveSubproblemCIP(set->scip, benders, probnumber, &solvestatus, solvecip) );

         if( solvestatus == SCIP_STATUS_INFEASIBLE )
            (*infeasible) = TRUE;
         if( objective != NULL )
            (*objective) = SCIPgetSolOrigObj(subproblem, SCIPgetBestSol(subproblem))*(int)SCIPgetObjsense(subproblem);
      }
      else
      {
         SCIP_Bool success;

         /* if the subproblem has convex constraints and continuous variables, then it should have been initialised and
          * in SCIP_STAGE_SOLVING. In this case, the subproblem only needs to be put into probing mode.
          */
         if( SCIPbendersGetSubproblemType(benders, probnumber) == SCIP_BENDERSSUBTYPE_CONVEXCONT )
         {
            /* if the subproblem is not in probing mode, then it must be put into that mode for the LP solve. */
            if( !SCIPinProbing(subproblem) )
            {
               SCIP_CALL( SCIPstartProbing(subproblem) );
            }

            success = TRUE;
         }
         else
         {
            SCIP_CALL( initialiseSubproblem(benders, set, probnumber, &success) );
         }

         /* if setting up the subproblem was successful */
         if( success )
         {
            SCIP_STATUS solvestatus;
            SCIP_Real lpobjective;

            SCIP_CALL( SCIPbendersSolveSubproblemLP(set->scip, benders, probnumber, &solvestatus, &lpobjective) );

            if( solvestatus == SCIP_STATUS_INFEASIBLE )
               (*infeasible) = TRUE;
            else if( objective != NULL )
               (*objective) = lpobjective;
         }
         else
         {
            if( objective != NULL )
               (*objective) = SCIPinfinity(subproblem);
         }
      }
   }

   return SCIP_OKAY;
}

/** copies the time and memory limit from the master problem to the subproblem */
static
SCIP_RETCODE copyMemoryAndTimeLimits(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 subproblem          /**< the Benders' decomposition subproblem */
   )
{
   SCIP_Real mastertimelimit;
   SCIP_Real subtimelimit;
   SCIP_Real maxsubtimelimit;
   SCIP_Real mastermemorylimit;
   SCIP_Real submemorylimit;
   SCIP_Real maxsubmemorylimit;

   assert(scip != NULL);

   /* setting the time limit for the Benders' decomposition subproblems. It is set to 102% of the remaining time. */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &mastertimelimit) );
   maxsubtimelimit = SCIPparamGetRealMax(SCIPgetParam(subproblem, "limits/time"));
   subtimelimit = (mastertimelimit - SCIPgetSolvingTime(scip)) * 1.02;
   subtimelimit = MIN(subtimelimit, maxsubtimelimit);
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/time", MAX(0.0, subtimelimit)) );

   /* setting the memory limit for the Benders' decomposition subproblems. */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &mastermemorylimit) );
   maxsubmemorylimit = SCIPparamGetRealMax(SCIPgetParam(subproblem, "limits/memory"));
   submemorylimit = mastermemorylimit - (SCIPgetMemUsed(scip) + SCIPgetMemExternEstim(scip))/1048576.0;
   submemorylimit = MIN(submemorylimit, maxsubmemorylimit);
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/memory", MAX(0.0, submemorylimit)) );

   return SCIP_OKAY;
}

/** stores the original parameters from the subproblem */
static
SCIP_RETCODE storeOrigSubproblemParams(
   SCIP*                 subproblem,         /**< the SCIP data structure */
   SCIP_SUBPROBPARAMS*   origparams          /**< the original subproblem parameters */
   )
{
   assert(subproblem != NULL);
   assert(origparams != NULL);

   SCIP_CALL( SCIPgetRealParam(subproblem, "limits/memory", &origparams->limits_memory) );
   SCIP_CALL( SCIPgetRealParam(subproblem, "limits/time", &origparams->limits_time) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "conflict/enable", &origparams->conflict_enable) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &origparams->lp_disablecutoff) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/scaling", &origparams->lp_scaling) );
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/initalgorithm", &origparams->lp_initalg) );
   SCIP_CALL( SCIPgetCharParam(subproblem, "lp/resolvealgorithm", &origparams->lp_resolvealg) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "lp/alwaysgetduals", &origparams->lp_alwaysgetduals) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/scaleobj", &origparams->misc_scaleobj) );
   SCIP_CALL( SCIPgetBoolParam(subproblem, "misc/catchctrlc", &origparams->misc_catchctrlc) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxrounds", &origparams->prop_maxrounds) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "propagating/maxroundsroot", &origparams->prop_maxroundsroot) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "constraints/linear/propfreq", &origparams->cons_linear_propfreq) );

   return SCIP_OKAY;
}

/** sets the parameters for the subproblem */
static
SCIP_RETCODE setSubproblemParams(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 subproblem          /**< the subproblem SCIP instance */
   )
{
   assert(scip != NULL);
   assert(subproblem != NULL);

   /* copying memory and time limits */
   SCIP_CALL( copyMemoryAndTimeLimits(scip, subproblem) );

   /* Do we have to disable presolving? If yes, we have to store all presolving parameters. */
   SCIP_CALL( SCIPsetPresolving(subproblem, SCIP_PARAMSETTING_OFF, TRUE) );

   /* Disabling heuristics so that the problem is not trivially solved */
   SCIP_CALL( SCIPsetHeuristics(subproblem, SCIP_PARAMSETTING_OFF, TRUE) );

   /* store parameters that are changed for the generation of the subproblem cuts */
   SCIP_CALL( SCIPsetBoolParam(subproblem, "conflict/enable", FALSE) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", 1) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/scaling", 0) );

   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/initalgorithm", 'd') );
   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/resolvealgorithm", 'd') );

   SCIP_CALL( SCIPsetBoolParam(subproblem, "lp/alwaysgetduals", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/scaleobj", FALSE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/catchctrlc", FALSE) );

#ifndef SCIP_MOREDEBUG
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
#endif

   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", 0) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "constraints/linear/propfreq", -1) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "heuristics/alns/freq", -1) );

   SCIP_CALL( SCIPsetIntParam(subproblem, "separating/aggregation/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "separating/gomory/freq", -1) );

   return SCIP_OKAY;
}

/** resets the original parameters from the subproblem */
static
SCIP_RETCODE resetOrigSubproblemParams(
   SCIP*                 subproblem,         /**< the SCIP data structure */
   SCIP_SUBPROBPARAMS*   origparams          /**< the original subproblem parameters */
   )
{
   assert(subproblem != NULL);
   assert(origparams != NULL);

   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/memory", origparams->limits_memory) );
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/time", origparams->limits_time) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "conflict/enable", origparams->conflict_enable) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", origparams->lp_disablecutoff) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "lp/scaling", origparams->lp_scaling) );
   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/initalgorithm", origparams->lp_initalg) );
   SCIP_CALL( SCIPsetCharParam(subproblem, "lp/resolvealgorithm", origparams->lp_resolvealg) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "lp/alwaysgetduals", origparams->lp_alwaysgetduals) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/scaleobj", origparams->misc_scaleobj) );
   SCIP_CALL( SCIPsetBoolParam(subproblem, "misc/catchctrlc", origparams->misc_catchctrlc) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxrounds", origparams->prop_maxrounds) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "propagating/maxroundsroot", origparams->prop_maxroundsroot) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "constraints/linear/propfreq", origparams->cons_linear_propfreq) );

   return SCIP_OKAY;
}

/** solves the LP of the Benders' decomposition subproblem
 *
 *  This requires that the subproblem is in probing mode.
 */
SCIP_RETCODE SCIPbendersSolveSubproblemLP(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_STATUS*          solvestatus,        /**< status of subproblem solve */
   SCIP_Real*            objective           /**< optimal value of subproblem, if solved to optimality */
   )
{
   SCIP* subproblem;
   SCIP_SUBPROBPARAMS* origparams;
   SCIP_Bool solvenlp;

   assert(benders != NULL);
   assert(solvestatus != NULL);
   assert(objective != NULL);
   assert(SCIPbendersSubproblemIsSetup(benders, probnumber));

   /* TODO: This should be solved just as an LP, so as a MIP. There is too much overhead with the MIP.
    * Need to change status check for checking the LP. */
   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* only solve the NLP relaxation if the NLP has been constructed and there exists an NLPI. If it is not possible to
    * solve the NLP relaxation, then the LP relaxation is used to generate Benders' cuts
    */
   solvenlp = FALSE;
   if( SCIPisNLPConstructed(subproblem) && SCIPgetNNlpis(subproblem) > 0
      && SCIPbendersGetSubproblemType(benders, probnumber) <= SCIP_BENDERSSUBTYPE_CONVEXDIS )
      solvenlp = TRUE;

   *objective = SCIPinfinity(subproblem);

   assert(SCIPisNLPConstructed(subproblem) || SCIPisLPConstructed(subproblem));
   assert(SCIPinProbing(subproblem));

   /* allocating memory for the parameter storage */
   SCIP_CALL( SCIPallocBlockMemory(subproblem, &origparams) );

   /* store the original parameters of the subproblem */
   SCIP_CALL( storeOrigSubproblemParams(subproblem, origparams) );

   /* setting the subproblem parameters */
   SCIP_CALL( setSubproblemParams(scip, subproblem) );

   if( solvenlp )
   {
      SCIP_NLPSOLSTAT nlpsolstat;
      SCIP_NLPTERMSTAT nlptermstat;
#ifdef SCIP_MOREDEBUG
      SCIP_SOL* nlpsol;
#endif

      SCIP_CALL( SCIPsolveNLP(subproblem) );  /*lint !e666*/

      nlpsolstat = SCIPgetNLPSolstat(subproblem);
      nlptermstat = SCIPgetNLPTermstat(subproblem);
      SCIPdebugMsg(scip, "NLP solstat %d termstat %d\n", nlpsolstat, nlptermstat);

      if( nlptermstat == SCIP_NLPTERMSTAT_OKAY && (nlpsolstat == SCIP_NLPSOLSTAT_LOCINFEASIBLE || nlpsolstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE) )
      {
         /* trust infeasible only if terminated "okay" */
         (*solvestatus) = SCIP_STATUS_INFEASIBLE;
      }
      else if( nlpsolstat == SCIP_NLPSOLSTAT_LOCOPT || nlpsolstat == SCIP_NLPSOLSTAT_GLOBOPT
         || nlpsolstat == SCIP_NLPSOLSTAT_FEASIBLE )
      {
#ifdef SCIP_MOREDEBUG
         SCIP_CALL( SCIPcreateNLPSol(subproblem, &nlpsol, NULL) );
         SCIP_CALL( SCIPprintSol(subproblem, nlpsol, NULL, FALSE) );
         SCIP_CALL( SCIPfreeSol(subproblem, &nlpsol) );
#endif

         (*solvestatus) = SCIP_STATUS_OPTIMAL;
         (*objective) = SCIPretransformObj(subproblem, SCIPgetNLPObjval(subproblem));
      }
      else if( nlpsolstat == SCIP_NLPSOLSTAT_UNBOUNDED )
      {
         (*solvestatus) = SCIP_STATUS_UNBOUNDED;
         SCIPerrorMessage("The NLP of Benders' decomposition subproblem %d is unbounded. This should not happen.\n",
            probnumber);
         SCIPABORT();
      }
      else if( nlptermstat == SCIP_NLPTERMSTAT_TIMELIMIT )
      {
         (*solvestatus) = SCIP_STATUS_TIMELIMIT;
      }
      else if( nlptermstat == SCIP_NLPTERMSTAT_INTERRUPT )
      {
         (*solvestatus) = SCIP_STATUS_USERINTERRUPT;
      }
      else
      {
         SCIPerrorMessage("Invalid solution status: %d. Termination status: %d. Solving the NLP relaxation of Benders' decomposition subproblem %d.\n",
            nlpsolstat, nlptermstat, probnumber);
         SCIPABORT();
      }
   }
   else
   {
      SCIP_Bool lperror;
      SCIP_Bool cutoff;

      SCIP_CALL( SCIPsolveProbingLP(subproblem, -1, &lperror, &cutoff) );

      switch( SCIPgetLPSolstat(subproblem) )
      {
         case SCIP_LPSOLSTAT_INFEASIBLE:
         {
            (*solvestatus) = SCIP_STATUS_INFEASIBLE;
            break;
         }

         case SCIP_LPSOLSTAT_OPTIMAL :
         {
            (*solvestatus) = SCIP_STATUS_OPTIMAL;
            (*objective) = SCIPgetSolOrigObj(subproblem, NULL)*(int)SCIPgetObjsense(scip);
            break;
         }

         case SCIP_LPSOLSTAT_UNBOUNDEDRAY :
         {
            (*solvestatus) = SCIP_STATUS_UNBOUNDED;
            SCIPerrorMessage("The LP of Benders' decomposition subproblem %d is unbounded. This should not happen.\n",
               probnumber);
            SCIPABORT();
            break;
         }

         case SCIP_LPSOLSTAT_ERROR :
         case SCIP_LPSOLSTAT_NOTSOLVED :
         case SCIP_LPSOLSTAT_TIMELIMIT :
         {
            if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_TIMELIMIT )
               (*solvestatus) = SCIP_STATUS_TIMELIMIT;
            else
               (*solvestatus) = SCIP_STATUS_UNKNOWN;

            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "   Benders' decomposition: Error solving LP "
               "relaxation of subproblem %d. No cut will be generated for this subproblem.\n", probnumber);
            break;
         }

         case SCIP_LPSOLSTAT_OBJLIMIT:
         case SCIP_LPSOLSTAT_ITERLIMIT:
         default:
         {
            SCIPerrorMessage("Invalid status: %d. Solving the LP relaxation of Benders' decomposition subproblem %d.\n",
               SCIPgetLPSolstat(subproblem), probnumber);
            SCIPABORT();
            break;
         }
      }
   }

   /* resetting the subproblem parameters */
   SCIP_CALL( resetOrigSubproblemParams(subproblem, origparams) );

   /* freeing the parameter storage */
   SCIPfreeBlockMemory(subproblem, &origparams);

   return SCIP_OKAY;
}

/** solves the Benders' decomposition subproblem */
SCIP_RETCODE SCIPbendersSolveSubproblemCIP(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_STATUS*          solvestatus,        /**< status of subproblem solve */
   SCIP_Bool             solvecip            /**< directly solve the CIP subproblem */
   )
{
   SCIP* subproblem;
   SCIP_SUBPROBPARAMS* origparams;

   assert(benders != NULL);
   assert(solvestatus != NULL);

   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* allocating memory for the parameter storage */
   SCIP_CALL( SCIPallocBlockMemory(subproblem, &origparams) );

   /* store the original parameters of the subproblem */
   SCIP_CALL( storeOrigSubproblemParams(subproblem, origparams) );

   /* If the solve has been stopped for the subproblem, then we need to restart it to complete the solve. The subproblem
    * is stopped when it is a MIP so that LP cuts and IP cuts can be generated. */
   if( SCIPgetStage(subproblem) == SCIP_STAGE_SOLVING )
   {
      /* the subproblem should be in probing mode. Otherwise, the event handler did not work correctly */
      assert( SCIPinProbing(subproblem) );

      /* the probing mode needs to be stopped so that the MIP can be solved */
      SCIP_CALL( SCIPendProbing(subproblem) );

      /* the problem was interrupted in the event handler, so SCIP needs to be informed that the problem is to be restarted */
      SCIP_CALL( SCIPrestartSolve(subproblem) );
   }
   else if( solvecip )
   {
      /* if the MIP will be solved directly, then the probing mode needs to be skipped.
       * This is achieved by setting the solvecip flag in the event handler data to TRUE
       */
      SCIP_EVENTHDLR* eventhdlr;
      SCIP_EVENTHDLRDATA* eventhdlrdata;

      eventhdlr = SCIPfindEventhdlr(subproblem, MIPNODEFOCUS_EVENTHDLR_NAME);
      eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

      eventhdlrdata->solvecip = TRUE;
   }
   else
   {
      /* if the problem is not in probing mode, then we need to solve the LP. That requires all methods that will
       * modify the structure of the problem need to be deactivated */

      /* setting the subproblem parameters */
      SCIP_CALL( setSubproblemParams(scip, subproblem) );

#ifdef SCIP_EVENMOREDEBUG
      SCIP_CALL( SCIPsetBoolParam(subproblem, "display/lpinfo", TRUE) );
#endif
   }

#ifdef SCIP_MOREDEBUG
      SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_FULL) );
      SCIP_CALL( SCIPsetIntParam(subproblem, "display/freq", 1) );
#endif

   SCIP_CALL( SCIPsolve(subproblem) );

   *solvestatus = SCIPgetStatus(subproblem);

   if( *solvestatus != SCIP_STATUS_OPTIMAL && *solvestatus != SCIP_STATUS_UNBOUNDED
      && *solvestatus != SCIP_STATUS_INFEASIBLE && *solvestatus != SCIP_STATUS_USERINTERRUPT
      && *solvestatus != SCIP_STATUS_BESTSOLLIMIT && *solvestatus != SCIP_STATUS_TIMELIMIT
      && *solvestatus != SCIP_STATUS_MEMLIMIT )
   {
      SCIPerrorMessage("Invalid status: %d. Solving the CIP of Benders' decomposition subproblem %d.\n",
         *solvestatus, probnumber);
      SCIPABORT();
   }

   /* resetting the subproblem parameters */
   SCIP_CALL( resetOrigSubproblemParams(subproblem, origparams) );

   /* freeing the parameter storage */
   SCIPfreeBlockMemory(subproblem, &origparams);

   return SCIP_OKAY;
}

/** frees the subproblems */
SCIP_RETCODE SCIPbendersFreeSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(benders->bendersfreesub != NULL
      || (benders->bendersfreesub == NULL && benders->benderssolvesubconvex == NULL && benders->benderssolvesub == NULL));
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   if( benders->bendersfreesub != NULL )
   {
      SCIP_CALL( benders->bendersfreesub(set->scip, benders, probnumber) );
   }
   else
   {
      /* the subproblem is only freed if it is not independent */
      if( subproblemIsActive(benders, probnumber) )
      {
         SCIP* subproblem = SCIPbendersSubproblem(benders, probnumber);

         if( SCIPbendersGetSubproblemType(benders, probnumber) == SCIP_BENDERSSUBTYPE_CONVEXCONT )
         {
            /* ending probing mode to reset the current node. The probing mode will be restarted at the next solve */
            if( SCIPinProbing(subproblem) )
            {
               SCIP_CALL( SCIPendProbing(subproblem) );
            }
         }
         else
         {
            /* if the subproblems were solved as part of an enforcement stage, then they will still be in probing mode. The
             * probing mode must first be finished and then the problem can be freed */
            if( SCIPgetStage(subproblem) >= SCIP_STAGE_TRANSFORMED && SCIPinProbing(subproblem) )
            {
               SCIP_CALL( SCIPendProbing(subproblem) );
            }

            SCIP_CALL( SCIPfreeTransform(subproblem) );
         }
      }
   }

   /* setting the setup flag for the subproblem to FALSE */
   SCIPbendersSetSubproblemIsSetup(benders, probnumber, FALSE);
   return SCIP_OKAY;
}

/** compares the subproblem objective value with the auxiliary variable value for optimality */
SCIP_Bool SCIPbendersSubproblemIsOptimal(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP_Real auxiliaryvarval;
   SCIP_Bool optimal;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   optimal = FALSE;

   auxiliaryvarval = SCIPbendersGetAuxiliaryVarVal(benders, set, sol, probnumber);

   SCIPsetDebugMsg(set, "Subproblem %d - Auxiliary Variable: %g Subproblem Objective: %g Reldiff: %g Soltol: %g\n",
      probnumber, auxiliaryvarval, SCIPbendersGetSubproblemObjval(benders, probnumber),
      SCIPrelDiff(SCIPbendersGetSubproblemObjval(benders, probnumber), auxiliaryvarval), benders->solutiontol);

   if( SCIPrelDiff(SCIPbendersGetSubproblemObjval(benders, probnumber), auxiliaryvarval) < benders->solutiontol )
      optimal = TRUE;

   return optimal;
}

/** returns the value of the auxiliary variable value in a master problem solution */
SCIP_Real SCIPbendersGetAuxiliaryVarVal(
   SCIP_BENDERS*         benders,            /**< the benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP_VAR* auxiliaryvar;

   assert(benders != NULL);
   assert(set != NULL);

   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);
   assert(auxiliaryvar != NULL);

   return SCIPgetSolVal(set->scip, sol, auxiliaryvar);
}

/** Solves an independent subproblem to identify its lower bound. The lower bound is then used to update the bound on
 *  the auxiliary variable.
 */
SCIP_RETCODE SCIPbendersComputeSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber,         /**< the subproblem to be evaluated */
   SCIP_Real*            lowerbound,         /**< the lowerbound for the subproblem */
   SCIP_Bool*            infeasible          /**< was the subproblem found to be infeasible? */
   )
{
   SCIP* subproblem;
   SCIP_Real dualbound;
   SCIP_Real memorylimit;
   SCIP_Real timelimit;
   SCIP_Longint totalnodes;
   int disablecutoff;
   int verblevel;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;

   assert(benders != NULL);
   assert(set != NULL);

   if( benders->benderssolvesub != NULL || benders->benderssolvesubconvex != NULL )
   {
      (*lowerbound) = SCIPvarGetLbGlobal(SCIPbendersGetAuxiliaryVar(benders, probnumber));
      (*infeasible) = FALSE;

      SCIPinfoMessage(set->scip, NULL, "Benders' decomposition: a bendersSolvesub or bendersSolvesubconvex has been "
         "implemented. SCIPbendersComputeSubproblemLowerbound can not be executed.\n");
      SCIPinfoMessage(set->scip, NULL, "Set the auxiliary variable lower bound by calling "
         "SCIPbendersUpdateSubproblemLowerbound in bendersCreatesub. The auxiliary variable %d will remain as %g\n",
         probnumber, (*lowerbound));

      return SCIP_OKAY;
   }
   else
   {
      SCIPverbMessage(set->scip, SCIP_VERBLEVEL_FULL, NULL, "Benders' decomposition: Computing a lower bound for"
         " subproblem %d\n", probnumber);
   }

   /* getting the subproblem to evaluate */
   subproblem = SCIPbendersSubproblem(benders, probnumber);

   (*lowerbound) = -SCIPinfinity(subproblem);
   (*infeasible) = FALSE;

   SCIP_CALL( SCIPgetIntParam(subproblem, "display/verblevel", &verblevel) );
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );
#ifdef SCIP_MOREDEBUG
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", (int)SCIP_VERBLEVEL_HIGH) );
#endif

   /* copying memory and time limits */
   SCIP_CALL( SCIPgetRealParam(subproblem, "limits/time", &timelimit) );
   SCIP_CALL( SCIPgetRealParam(subproblem, "limits/memory", &memorylimit) );
   SCIP_CALL( copyMemoryAndTimeLimits(set->scip, subproblem) );

   /* if the subproblem is independent, then the default SCIP settings are used. Otherwise, only the root node is solved
    * to compute a lower bound on the subproblem
    */
   SCIP_CALL( SCIPgetLongintParam(subproblem, "limits/totalnodes", &totalnodes) );
   SCIP_CALL( SCIPgetIntParam(subproblem, "lp/disablecutoff", &disablecutoff) );
   if( !SCIPbendersSubproblemIsIndependent(benders, probnumber) )
   {
      SCIP_CALL( SCIPsetLongintParam(subproblem, "limits/totalnodes", 1LL) );
      SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", 1) );
   }

   /* if the subproblem not independent and is convex, then the probing LP is solved. Otherwise, the MIP is solved */
   dualbound = -SCIPinfinity(subproblem);
   if( SCIPbendersGetSubproblemType(benders, probnumber) == SCIP_BENDERSSUBTYPE_CONVEXCONT )
   {
      SCIP_Bool solvenlp = FALSE;

      assert(SCIPisLPConstructed(subproblem) || SCIPisNLPConstructed(subproblem));

      if( SCIPisNLPConstructed(subproblem) && SCIPgetNNlpis(subproblem) > 0
         && SCIPbendersGetSubproblemType(benders, probnumber) <= SCIP_BENDERSSUBTYPE_CONVEXDIS )
         solvenlp = TRUE;

      SCIP_CALL( SCIPstartProbing(subproblem) );
      if( solvenlp )
      {
         SCIP_NLPSOLSTAT nlpsolstat;
         SCIP_NLPTERMSTAT nlptermstat;

         SCIP_CALL( SCIPsolveNLP(subproblem) );  /*lint !e666*/

         nlpsolstat = SCIPgetNLPSolstat(subproblem);
         nlptermstat = SCIPgetNLPTermstat(subproblem);
         SCIPdebugMsg(set->scip, "NLP solstat %d termstat %d\n", nlpsolstat, nlptermstat);

         if( nlptermstat == SCIP_NLPTERMSTAT_OKAY && (nlpsolstat == SCIP_NLPSOLSTAT_LOCINFEASIBLE || nlpsolstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE) )
         {
            /* trust infeasible only if terminated "okay" */
            (*infeasible) = TRUE;
         }
         else if( nlpsolstat == SCIP_NLPSOLSTAT_LOCOPT || nlpsolstat == SCIP_NLPSOLSTAT_GLOBOPT
            || nlpsolstat == SCIP_NLPSOLSTAT_FEASIBLE )
         {
            dualbound = SCIPretransformObj(subproblem, SCIPgetNLPObjval(subproblem));
         }
      }
      else
      {
         SCIP_CALL( SCIPsolveProbingLP(subproblem, -1, &lperror, &cutoff) );

         if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_INFEASIBLE )
            (*infeasible) = TRUE;
         else if( SCIPgetLPSolstat(subproblem) == SCIP_LPSOLSTAT_OPTIMAL )
            dualbound = SCIPgetSolOrigObj(subproblem, NULL)*(int)SCIPgetObjsense(set->scip);
      }
   }
   else
   {
      SCIP_EVENTHDLRDATA* eventhdlrdata;

      /* if the subproblem is not convex, then event handlers have been added to interrupt the solve. These must be
       * disabled
       */
      eventhdlrdata = SCIPeventhdlrGetData(SCIPfindEventhdlr(subproblem, MIPNODEFOCUS_EVENTHDLR_NAME));
      eventhdlrdata->solvecip = TRUE;

      SCIP_CALL( SCIPsolve(subproblem) );

      if( SCIPgetStatus(subproblem) == SCIP_STATUS_INFEASIBLE )
         (*infeasible) = TRUE;
      else
         dualbound = SCIPgetDualbound(subproblem);
   }

   /* getting the lower bound value */
   (*lowerbound) = dualbound;

   if( !SCIPbendersSubproblemIsIndependent(benders, probnumber) )
   {
      SCIP_CALL( SCIPsetLongintParam(subproblem, "limits/totalnodes", totalnodes) );
      SCIP_CALL( SCIPsetIntParam(subproblem, "lp/disablecutoff", disablecutoff) );
   }
   SCIP_CALL( SCIPsetIntParam(subproblem, "display/verblevel", verblevel) );
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetRealParam(subproblem, "limits/time", timelimit) );

   /* the subproblem must be freed so that it is reset for the subsequent Benders' decomposition solves. If the
    * subproblems are independent, they are not freed. SCIPfreeBendersSubproblem must still be called, but in this
    * function the independent subproblems are not freed. However, they will still be freed at the end of the
    * solving process for the master problem.
    */
   SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, probnumber) );

   return SCIP_OKAY;
}

/** Merges a subproblem into the master problem. This process just adds a copy of the subproblem variables and
 *  constraints to the master problem, but keeps the subproblem stored in the Benders' decomposition data structure. The reason for
 *  keeping the subproblem available is for when it is queried for solutions after the problem is solved.
 *
 *  Once the subproblem is merged into the master problem, then the subproblem is flagged as disabled. This means that
 *  it will not be solved in the subsequent subproblem solving loops.
 *
 *  The associated auxiliary variables are kept in the master problem. The objective function of the merged subproblem
 *  is added as an underestimator constraint.
 */
SCIP_RETCODE SCIPbendersMergeSubproblemIntoMaster(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of subproblem variables corresponding
                                              *   to the newly created master variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of subproblem constraints to the
                                              *   corresponding newly created constraints, or NULL */
   int                   probnumber          /**< the number of the subproblem that will be merged into the master problem*/
   )
{
   SCIP* subproblem;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_VAR** vars;
   SCIP_VAR* auxiliaryvar;
   SCIP_CONS** conss;
   SCIP_CONS* objcons;
   int nvars;
   int nconss;
   int i;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   char varname[SCIP_MAXSTRLEN];
   char consname[SCIP_MAXSTRLEN];
   const char* origvarname;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   SCIPverbMessage(set->scip, SCIP_VERBLEVEL_HIGH, NULL, "   Benders' decomposition: Infeasibility of subproblem %d can't "
      "be resolved. Subproblem %d is being merged into the master problem.\n", probnumber, probnumber);

   /* freeing the subproblem because it will be flagged as independent. Since the subproblem is flagged as independent,
    * it will no longer be solved or freed within the solving loop.
    */
   SCIP_CALL( SCIPbendersFreeSubproblem(benders, set, probnumber) );

   subproblem = SCIPbendersSubproblem(benders, probnumber);

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(set->scip), SCIPgetNVars(subproblem)) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(set->scip), SCIPgetNConss(subproblem)) );
   }
   else
      localconsmap = consmap;

   /* retrieving the subproblem variable to build a subproblem mapping */
   vars = SCIPgetVars(subproblem);
   nvars = SCIPgetNVars(subproblem);

   /* creating the objective function constraint that will be added to the master problem */
   /* setting the name of the transferred cut */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "objectivecons_%d", probnumber );
   SCIP_CALL( SCIPcreateConsBasicLinear(set->scip, &objcons, consname, 0, NULL, NULL, -SCIPsetInfinity(set), 0.0) );
   SCIP_CALL( SCIPsetConsRemovable(set->scip, objcons, TRUE) );

   for( i = 0; i < nvars; i++ )
   {
      SCIP_VAR* mastervar = NULL;
      SCIP_Bool releasevar = FALSE;

      SCIP_CALL( SCIPgetBendersMasterVar(set->scip, benders, vars[i], &mastervar) );

      /* if the master problem variable is not NULL, then there is a corresponding variable in the master problem for
       * the given subproblem variable. In this case, the variable is added to the hashmap.
       */
      if( mastervar == NULL )
      {
         SCIP_VAR* origvar;
         SCIP_Real scalar;
         SCIP_Real constant;

         /* This is following the same process as in createVariableMappings. The original variable is used to map
          * between the subproblem and the master problem
          */
         origvar = vars[i];
         scalar = 1.0;
         constant = 0.0;
         SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

         /* retrieving the var name */
         origvarname = SCIPvarGetName(origvar);
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s", origvarname);

         /* creating and adding the variable to the Benders' decomposition master problem */
         SCIP_CALL( SCIPcreateVarBasic(set->scip, &mastervar, varname, SCIPvarGetLbOriginal(origvar),
            SCIPvarGetUbOriginal(origvar), 0.0, SCIPvarGetType(origvar)) );

         /* adding the variable to the master problem */
         SCIP_CALL( SCIPaddVar(set->scip, mastervar) );

         /* adds the variable to the objective function constraint */
         SCIP_CALL( SCIPaddCoefLinear(set->scip, objcons, mastervar, SCIPvarGetObj(origvar)) );

         /* the variable must be released */
         releasevar = TRUE;
      }

      /* creating the mapping betwen the subproblem var and the master var for the constraint copying */
      SCIP_CALL( SCIPhashmapInsert(localvarmap, vars[i], mastervar) );

      /* releasing the variable */
      if( releasevar )
      {
         SCIP_CALL( SCIPreleaseVar(set->scip, &mastervar) );
      }
   }

   /* getting the constraints from the subproblem that will be added to the master problem */
   conss = SCIPgetConss(subproblem);
   nconss = SCIPgetNConss(subproblem);

   /* getting a copy of all constraints and adding it to the master problem */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CONS* targetcons;
      SCIP_Bool initial;
      SCIP_Bool valid;

      /* NOTE: adding all subproblem constraints appears to cause an error when resolving the LP, which results in the
       * current incumbent being reported as optimal. To avoid this, only half of the subproblem constraints are added
       * the master problem. The remaining half are marked as lazy and are separated as required.
       */
      initial = (i < nconss/2);

      SCIP_CALL( SCIPgetConsCopy(subproblem, set->scip, conss[i], &targetcons, SCIPconsGetHdlr(conss[i]),
         localvarmap, localconsmap, NULL, initial, SCIPconsIsSeparated(conss[i]),
         SCIPconsIsEnforced(conss[i]), SCIPconsIsChecked(conss[i]), SCIPconsIsPropagated(conss[i]), FALSE,
         SCIPconsIsModifiable(conss[i]), SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]),
         FALSE, TRUE, &valid) );
      assert(SCIPconsIsInitial(conss[i]));
      assert(valid);

      SCIP_CALL( SCIPaddCons(set->scip, targetcons) );

      SCIP_CALL( SCIPreleaseCons(set->scip, &targetcons) );
   }

   /* freeing the hashmaps */
   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   /* adding the auxiliary variable to the objective constraint */
   auxiliaryvar = SCIPbendersGetAuxiliaryVar(benders, probnumber);
   SCIP_CALL( SCIPaddCoefLinear(set->scip, objcons, auxiliaryvar, -1.0) );

   /* adding the objective function constraint to the master problem */
   SCIP_CALL( SCIPaddCons(set->scip, objcons) );

   SCIP_CALL( SCIPreleaseCons(set->scip, &objcons) );

   /* the merged subproblem is no longer solved. This is indicated by setting the subproblem as disabled. The
    * subproblem still exists, but it is not solved in the solving loop.
    */
   SCIPbendersSetSubproblemEnabled(benders, probnumber, FALSE);

   return SCIP_OKAY;
}

/** when applying a decomposition from a supplied format, constraints must be transferred from the master problem to the
 *  subproblem. This is achieved by adding new constraints to the subproblem
 */
static
SCIP_RETCODE addConstraintToBendersSubproblem(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP*                 subproblem,         /**< the SCIP instance for the subproblem */
   SCIP_HASHMAP*         varmap,             /**< the variable hash map mapping the source variables to the target variables */
   SCIP_CONS*            sourcecons          /**< the constraint that being added to the subproblem */
   )
{
   SCIP* scip;
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   int nconsvars;
   int i;
   SCIP_Bool success;

   assert(set != NULL);
   assert(subproblem != NULL);
   assert(varmap != NULL);
   assert(sourcecons != NULL);

   SCIPdebugMessage("Adding constraint <%s> to Benders' decomposition subproblem\n", SCIPconsGetName(sourcecons));

   scip = set->scip;

   /* getting the variables that are in the constraint */
   SCIP_CALL( SCIPgetConsNVars(scip, sourcecons, &nconsvars, &success) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );

   SCIP_CALL( SCIPgetConsVars(scip, sourcecons, consvars, nconsvars, &success) );
   assert(success);

   /* checking all variables to see whether they already exist in the subproblem. If they don't exist, then the variable
    * is created
    */
   for( i = 0; i < nconsvars; i++ )
   {
      /* if the variable is not in the hashmap, then it doesn't exist in the subproblem */
      if( !SCIPhashmapExists(varmap, consvars[i]) )
      {
         SCIP_VAR* var;

         /* creating a variable as a copy of the original variable. */
         SCIP_CALL( SCIPcreateVar(subproblem, &var, SCIPvarGetName(consvars[i]), SCIPvarGetLbGlobal(consvars[i]),
               SCIPvarGetUbGlobal(consvars[i]), SCIPvarGetObj(consvars[i]), SCIPvarGetType(consvars[i]),
               SCIPvarIsInitial(consvars[i]), SCIPvarIsRemovable(consvars[i]), NULL, NULL, NULL, NULL, NULL) );

         /* adding the variable to the subproblem */
         SCIP_CALL( SCIPaddVar(subproblem, var) );

         /* adding the variable to the hash map so that it is copied correctly in the constraint */
         SCIP_CALL( SCIPhashmapInsert(varmap, consvars[i], var) );

         /* releasing the variable */
         SCIP_CALL( SCIPreleaseVar(subproblem, &var) );
      }
   }

   /* freeing the buffer memory for the consvars */
   SCIPfreeBufferArray(scip, &consvars);

   /* copying the constraint from the master scip to the subproblem */
   SCIP_CALL( SCIPgetConsCopy(scip, subproblem, sourcecons, &cons, SCIPconsGetHdlr(sourcecons), varmap, NULL,
         SCIPconsGetName(sourcecons), SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons), SCIPconsIsMarkedPropagate(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons), TRUE, &success) );

   /* if the copy failed, then the subproblem for the decomposition could not be performed. */
   if( !success )
   {
      SCIPerrorMessage("It is not possible to copy constraint <%s>. Benders' decomposition could not be applied.\n",
         SCIPconsGetName(sourcecons));
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPaddCons(subproblem, cons) );
   SCIP_CALL( SCIPreleaseCons(subproblem, &cons) );

   return SCIP_OKAY;
}

/** removes the variables and constraints from the master problem that have been transferred to a subproblem when the
 *  decomposition was applied.
 */
static
SCIP_RETCODE removeVariablesAndConstraintsFromMaster(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP_CONS**           conss,              /**< the master problem constraints */
   SCIP_VAR**            vars,               /**< the master problem variables, can be NULL */
   int*                  conslabels,         /**< the labels indicating the block for each constraint */
   int*                  varslabels,         /**< the labels indicating the block for each variable, can be NULL */
   int                   nconss,             /**< the number of constraints */
   int                   nvars               /**< the number of variables */
   )
{
   int i;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(conslabels != NULL);
   assert((vars != NULL && varslabels != NULL) || (vars == NULL && varslabels == NULL));

   /* removing constraints */
   for( i = nconss - 1; i >= 0; i-- )
   {
      if( conslabels[i] >= 0 && !SCIPconsIsDeleted(conss[i]) )
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );
   }

   /* removing variables */
   if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM && vars != NULL && varslabels != NULL )
   {
      for( i = nvars - 1; i >= 0; i-- )
      {
         if( varslabels[i] >= 0 && !SCIPvarIsDeleted(vars[i]) )
         {
            SCIP_Bool deleted;

            SCIP_CALL( SCIPdelVar(scip, vars[i], &deleted) );
            assert(deleted);
         }
      }
   }

   return SCIP_OKAY;
}

/** Applies a Benders' decomposition to the problem based upon the decomposition selected from the storage */
SCIP_RETCODE SCIPbendersApplyDecomposition(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DECOMP*          decomp              /**< the decomposition to apply to the problem */
   )
{
   SCIP** subproblems;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_HASHMAP** varmaps;
   int* varslabels;
   int* conslabels;
   int nvars;
   int nconss;
   int nblocks;
   int i;
   char subprobname[SCIP_MAXSTRLEN];

   assert(benders != NULL);
   assert(set != NULL);
   assert(decomp != NULL);

   SCIPdebugMessage("Applying a Benders' decomposition to <%s>\n", SCIPgetProbName(set->scip));

   /* retrieving the number of blocks for this decomposition */
   nblocks = SCIPdecompGetNBlocks(decomp);
   assert(nblocks > 0);

   /* initialising the subproblems for the Benders' decomposition */
   SCIP_CALL( SCIPallocBufferArray(set->scip, &subproblems, nblocks) );

   /* creating the subproblems before adding the constraints */
   for( i = 0; i < nblocks; i++ )
   {
      SCIP_Bool valid;

      SCIP_CALL( SCIPcreate(&subproblems[i]) );

      /* copying the plugins from the original SCIP instance to the subproblem SCIP */
      SCIP_CALL( SCIPcopyPlugins(set->scip, subproblems[i], TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, &valid) );

      (void) SCIPsnprintf(subprobname, SCIP_MAXSTRLEN, "sub_%s_%d", SCIPgetProbName(set->scip), i);
      SCIP_CALL( SCIPcreateProbBasic(subproblems[i], subprobname) );
   }

   /* TODO: Need to work out whether a check for original and transformed problem is necessary */

   /* getting the variables and constraints from the problem */
   SCIP_CALL( SCIPgetVarsData(set->scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   conss = SCIPgetConss(set->scip);
   nconss = SCIPgetNConss(set->scip);

   /* allocating buffer memory for the labels arrays */
   SCIP_CALL( SCIPallocBufferArray(set->scip, &varslabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(set->scip, &conslabels, nconss) );

   /* getting the labels for the variables and constraints from the decomposition */
   SCIPdecompGetVarsLabels(decomp, vars, varslabels, nvars);
   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

   /* creating the variable maps for adding the constraints to the subproblems */
   SCIP_CALL( SCIPallocBufferArray(set->scip, &varmaps, nblocks) );

   for( i = 0; i < nblocks; i++ )
   {
      SCIP_CALL( SCIPhashmapCreate(&varmaps[i], SCIPblkmem(subproblems[i]), nvars) );
   }

   /* copying the constraints to the appropriate subproblems */
   for( i = 0; i < nconss; i++ )
   {
      /* we are only interested in the constraints that are in the blocks. These are identified by a label >= 0 */
      if( conslabels[i] >= 0 )
      {
         SCIP_CALL( addConstraintToBendersSubproblem(set, subproblems[conslabels[i]], varmaps[conslabels[i]],
               conss[i]) );
      }
   }

   /* removing the variables and constraints from the master problem that have been added to the subproblem */
   SCIP_CALL( removeVariablesAndConstraintsFromMaster(set->scip, conss, vars, conslabels, varslabels, nconss, nvars) );

   /* creating the Benders' decomposition my calling the default plugin */
   SCIP_CALL( SCIPcreateBendersDefault(set->scip, subproblems, nblocks) );

   /* flag to the Benders' decomposition core that the subproblems need to be freed */
   benders->freesubprobs = TRUE;

   /* activating the Benders' constraint handler for the scenario stages.
    * TODO: consider whether the two-phase method should be activated by default in the scenario stages.
    */
   SCIP_CALL( SCIPsetBoolParam(set->scip, "constraints/benders/active", TRUE) );

   /* changing settings that are required for Benders' decomposition */
   SCIP_CALL( SCIPsetPresolving(set->scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(set->scip, "propagating/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(set->scip, "propagating/maxroundsroot", 0) );
   SCIP_CALL( SCIPsetIntParam(set->scip, "heuristics/trysol/freq", 1) );

   /* disabling aggregation since it can affect the mapping between the master and subproblem variables */
   SCIP_CALL( SCIPsetBoolParam(set->scip, "presolving/donotaggr", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(set->scip, "presolving/donotmultaggr", TRUE) );

   /* freeing the allocated memory */
   for( i = nblocks - 1; i >= 0; i-- )
   {
      SCIPhashmapFree(&varmaps[i]);
   }

   SCIPfreeBufferArray(set->scip, &varmaps);
   SCIPfreeBufferArray(set->scip, &conslabels);
   SCIPfreeBufferArray(set->scip, &varslabels);
   SCIPfreeBufferArray(set->scip, &subproblems);

   return SCIP_OKAY;
}

/** Returns the corresponding master or subproblem variable for the given variable.
 *  This provides a call back for the variable mapping between the master and subproblems. */
SCIP_RETCODE SCIPbendersGetVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< the variable for which the corresponding variable is desired */
   SCIP_VAR**            mappedvar,          /**< the variable that is mapped to var */
   int                   probnumber          /**< the problem number for the desired variable, -1 for the master problem */
   )
{
   assert(benders != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(mappedvar != NULL);
   assert(benders->bendersgetvar != NULL);

   (*mappedvar) = NULL;

   /* if the variable name matches the auxiliary variable, then the master variable is returned as NULL */
   if( strstr(SCIPvarGetName(var), AUXILIARYVAR_NAME) != NULL )
      return SCIP_OKAY;

   SCIP_CALL( benders->bendersgetvar(set->scip, benders, var, mappedvar, probnumber) );

   return SCIP_OKAY;
}

/** gets user data of Benders' decomposition */
SCIP_BENDERSDATA* SCIPbendersGetData(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->bendersdata;
}

/** sets user data of Benders' decomposition; user has to free old data in advance! */
void SCIPbendersSetData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSDATA*     bendersdata         /**< new Benders' decomposition user data */
   )
{
   assert(benders != NULL);

   benders->bendersdata = bendersdata;
}

/** sets copy callback of Benders' decomposition */
void SCIPbendersSetCopy(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSCOPY ((*benderscopy))    /**< copy callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->benderscopy = benderscopy;
}

/** sets destructor callback of Benders' decomposition */
void SCIPbendersSetFree(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREE ((*bendersfree))    /**< destructor of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersfree = bendersfree;
}

/** sets initialization callback of Benders' decomposition */
void SCIPbendersSetInit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINIT((*bendersinit))     /**< initialize the Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinit = bendersinit;
}

/** sets deinitialization callback of Benders' decomposition */
void SCIPbendersSetExit(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXIT((*bendersexit))     /**< deinitialize the Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexit = bendersexit;
}

/** sets presolving initialization callback of Benders' decomposition */
void SCIPbendersSetInitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre))/**< initialize presolving for Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinitpre = bendersinitpre;
}

/** sets presolving deinitialization callback of Benders' decomposition */
void SCIPbendersSetExitpre(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre))/**< deinitialize presolving for Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexitpre = bendersexitpre;
}

/** sets solving process initialization callback of Benders' decomposition */
void SCIPbendersSetInitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol))/**< solving process initialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersinitsol = bendersinitsol;
}

/** sets solving process deinitialization callback of Benders' decomposition */
void SCIPbendersSetExitsol(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol))/**< solving process deinitialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->bendersexitsol = bendersexitsol;
}

/** sets the pre subproblem solve callback of Benders' decomposition */
void SCIPbendersSetPresubsolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve))/**< called prior to the subproblem solving loop */
   )
{
   assert(benders != NULL);

   benders->benderspresubsolve = benderspresubsolve;
}

/** sets convex solve callback of Benders' decomposition */
void SCIPbendersSetSolvesubconvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex))/**< solving method for the convex Benders' decomposition subproblem */
   )
{
   assert(benders != NULL);

   benders->benderssolvesubconvex = benderssolvesubconvex;
}

/** sets solve callback of Benders' decomposition */
void SCIPbendersSetSolvesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub))/**< solving method for a Benders' decomposition subproblem */
   )
{
   assert(benders != NULL);

   benders->benderssolvesub = benderssolvesub;
}

/** sets post-solve callback of Benders' decomposition */
void SCIPbendersSetPostsolve(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve))/**< solving process deinitialization callback of Benders' decomposition */
   )
{
   assert(benders != NULL);

   benders->benderspostsolve = benderspostsolve;
}

/** sets post-solve callback of Benders' decomposition */
void SCIPbendersSetSubproblemComp(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_SORTPTRCOMP((*benderssubcomp))  /**< a comparator for defining the solving order of the subproblems */
   )
{
   assert(benders != NULL);

   benders->benderssubcomp = benderssubcomp;
}

/** sets free subproblem callback of Benders' decomposition */
void SCIPbendersSetFreesub(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub))/**< the freeing callback for the subproblem */
   )
{
   assert(benders != NULL);

   benders->bendersfreesub = bendersfreesub;
}

/** gets name of Benders' decomposition */
const char* SCIPbendersGetName(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->name;
}

/** gets description of Benders' decomposition */
const char* SCIPbendersGetDesc(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->desc;
}

/** gets priority of Benders' decomposition */
int SCIPbendersGetPriority(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->priority;
}

/** sets priority of Benders' decomposition */
void SCIPbendersSetPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(set != NULL);

   benders->priority = priority;
   set->benderssorted = FALSE;
}

/** gets the number of subproblems for the Benders' decomposition */
int SCIPbendersGetNSubproblems(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   )
{
   assert(benders != NULL);

   return benders->nsubproblems;
}

/** returns the SCIP instance for a given subproblem */
SCIP* SCIPbendersSubproblem(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition data structure */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < benders->nsubproblems);

   return benders->subproblems[probnumber];
}

/** gets the number of times, the Benders' decomposition was called and tried to find a variable with negative reduced costs */
int SCIPbendersGetNCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->ncalls;
}

/** gets the number of optimality cuts found by the collection of Benders' decomposition subproblems */
int SCIPbendersGetNCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->ncutsfound;
}

/** gets the number of cuts found from the strengthening round */
int SCIPbendersGetNStrengthenCutsFound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nstrengthencuts;
}

/** gets the number of calls to the strengthening round */
int SCIPbendersGetNStrengthenCalls(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nstrengthencalls;
}

/** gets the number of calls to the strengthening round that fail */
int SCIPbendersGetNStrengthenFails(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nstrengthenfails;
}

/** gets time in seconds used in this Benders' decomposition for setting up for next stages */
SCIP_Real SCIPbendersGetSetupTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->setuptime);
}

/** gets time in seconds used in this Benders' decomposition */
SCIP_Real SCIPbendersGetTime(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return SCIPclockGetTime(benders->bendersclock);
}

/** enables or disables all clocks of the Benders' decomposition, depending on the value of the flag */
void SCIPbendersEnableOrDisableClocks(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the Benders' decomposition be enabled? */
   )
{
   assert(benders != NULL);

   SCIPclockEnableOrDisable(benders->setuptime, enable);
   SCIPclockEnableOrDisable(benders->bendersclock, enable);
}

/** is Benders' decomposition initialized? */
SCIP_Bool SCIPbendersIsInitialized(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->initialized;
}

/** Are Benders' cuts generated from the LP solutions? */
SCIP_Bool SCIPbendersCutLP(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutlp;
}

/** Are Benders' cuts generated from the pseudo solutions? */
SCIP_Bool SCIPbendersCutPseudo(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutpseudo;
}

/** Are Benders' cuts generated from the relaxation solutions? */
SCIP_Bool SCIPbendersCutRelaxation(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->cutrelax;
}

/** should this Benders' use the auxiliary variables from the highest priority Benders' */
SCIP_Bool SCIPbendersShareAuxVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->shareauxvars;
}

/** adds a subproblem to the Benders' decomposition data. If a custom subproblem solving method is used, then the
 *  subproblem pointer can be set to NULL
 */
SCIP_RETCODE SCIPbendersAddSubproblem(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP*                 subproblem          /**< subproblem to be added to the data storage, can be NULL */
   )
{
   assert(benders != NULL);
   assert(benders->naddedsubprobs + 1 <= benders->nsubproblems);

   /* if the subproblem pointer is NULL, then the subproblem solving callback functions must be set. */
   if( subproblem == NULL && (!benders->benderssolvesubconvex || !benders->benderssolvesub) )
   {
      SCIPerrorMessage("The subproblem can only be set to NULL if both bendersSolvesubconvex%s and bendersSolvesub%s "
         "are defined.\n", benders->name, benders->name);
      return SCIP_ERROR;
   }

   benders->subproblems[benders->naddedsubprobs] = subproblem;

   benders->naddedsubprobs++;

   return SCIP_OKAY;
}

/** removes the subproblems from the Benders' decomposition data */
void SCIPbendersRemoveSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(benders->subproblems != NULL);

   BMSclearMemoryArray(&benders->subproblems, benders->naddedsubprobs);
   benders->naddedsubprobs = 0;
}

/** returns the auxiliary variable for the given subproblem */
SCIP_VAR* SCIPbendersGetAuxiliaryVar(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->auxiliaryvars[probnumber];
}

/** returns all auxiliary variables */
SCIP_VAR** SCIPbendersGetAuxiliaryVars(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->auxiliaryvars;
}

/** stores the objective function value of the subproblem for use in cut generation */
void SCIPbendersSetSubproblemObjval(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition structure */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             objval              /**< the objective function value for the subproblem */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* updating the best objval */
   if( objval < benders->bestsubprobobjval[probnumber] )
      benders->bestsubprobobjval[probnumber] = objval;

   benders->subprobobjval[probnumber] = objval;
}

/** returns the objective function value of the subproblem for use in cut generation */
SCIP_Real SCIPbendersGetSubproblemObjval(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobobjval[probnumber];
}

/** returns whether the solution has non-zero slack variables */
SCIP_RETCODE SCIPbendersSolSlackVarsActive(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_Bool*            activeslack         /**< flag to indicate whether a slack variable is active */
   )
{
   SCIP* subproblem;
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   int nsubproblems;
   int nvars;
   int ncontvars;
   int i;
   int j;
   SCIP_Bool freesol = FALSE;

   assert(benders != NULL);
   assert(activeslack != NULL);

   (*activeslack) = FALSE;

   /* if the slack variables have not been added, then we can immediately state that no slack variables are active */
   if( !benders->feasibilityphase )
   {
      return SCIP_OKAY;
   }

   nsubproblems = SCIPbendersGetNSubproblems(benders);

   /* checking all subproblems for active slack variables */
   for( i = 0; i < nsubproblems && !(*activeslack); i++ )
   {
      subproblem = SCIPbendersSubproblem(benders, i);

      /* if the subproblem is convex and an NLP, then we need to create the NLP solution. Otherwise, the solution can be
       * retrieved from the LP or CIP.
       */
      if( SCIPbendersGetSubproblemType(benders, i) == SCIP_BENDERSSUBTYPE_CONVEXCONT )
      {
         if( SCIPisNLPConstructed(subproblem) && SCIPgetNNlpis(subproblem) > 0 )
         {
            SCIP_CALL( SCIPcreateNLPSol(subproblem, &sol, NULL) );
         }
         else
         {
            SCIP_CALL( SCIPcreateCurrentSol(subproblem, &sol, NULL) );
         }
         freesol = TRUE;
      }
      else
         sol = SCIPgetBestSol(subproblem);

      /* getting the variable data. Only the continuous variables are important. */
      SCIP_CALL( SCIPgetVarsData(subproblem, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );

      /* checking all slack variables for non-zero solution values */
      for( j = nvars - 1; j >= nvars - ncontvars; j-- )
      {
         if( strstr(SCIPvarGetName(vars[j]), SLACKVAR_NAME) != NULL )
         {
            if( !SCIPisZero(subproblem, SCIPgetSolVal(subproblem, sol, vars[j])) )
            {
               (*activeslack) = TRUE;
               break;
            }
         }
      }

      /* freeing the LP and NLP solutions */
      if( freesol )
      {
         SCIP_CALL( SCIPfreeSol(subproblem, &sol) );
      }
   }

   return SCIP_OKAY;
}

/** sets the subproblem type
 *
 * The subproblem types are:
 *    - Convex constraints with continuous variables
 *    - Convex constraints with discrete variables
 *    - Non-convex constraints with continuous variables
 *    - Non-convex constraints with discrete variables
 */
void SCIPbendersSetSubproblemType(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_BENDERSSUBTYPE   subprobtype         /**< the subproblem type */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   if( subprobtype == SCIP_BENDERSSUBTYPE_CONVEXCONT
      && benders->subprobtype[probnumber] != SCIP_BENDERSSUBTYPE_CONVEXCONT )
      benders->nconvexsubprobs++;
   else if( subprobtype != SCIP_BENDERSSUBTYPE_CONVEXCONT
      && benders->subprobtype[probnumber] == SCIP_BENDERSSUBTYPE_CONVEXCONT )
      benders->nconvexsubprobs--;

   benders->subprobtype[probnumber] = subprobtype;

   assert(benders->nconvexsubprobs >= 0 && benders->nconvexsubprobs <= benders->nsubproblems);
}

/** returns the type of the subproblem
 *
 *  This type is used to determine whether the duals of the problem can be used to generate cuts
 */
SCIP_BENDERSSUBTYPE SCIPbendersGetSubproblemType(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobtype[probnumber];
}

/** sets the flag indicating whether a subproblem is convex
 *
 *  It is possible that this can change during the solving process. One example is when the three-phase method is
 *  employed, where the first phase solves the convex relaxation of both the master and subproblems, the second phase
 *  reintroduces the integrality constraints to the master problem and the third phase then reintroduces integrality
 *  constraints to the subproblems.
 */
void SCIPbendersSetSubproblemIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isconvex            /**< flag to indicate whether the subproblem is convex */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   if( isconvex && !benders->subprobisconvex[probnumber] )
      benders->nconvexsubprobs++;
   else if( !isconvex && benders->subprobisconvex[probnumber] )
      benders->nconvexsubprobs--;

   benders->subprobisconvex[probnumber] = isconvex;

   assert(benders->nconvexsubprobs >= 0 && benders->nconvexsubprobs <= benders->nsubproblems);
}

/** returns whether the subproblem is convex
 *
 *  This means that the dual solution can be used to generate cuts.
 */
SCIP_Bool SCIPbendersSubproblemIsConvex(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobisconvex[probnumber];
}

/** returns the number of subproblems that are convex */
int SCIPbendersGetNConvexSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nconvexsubprobs;
}

/** sets the flag indicating whether a subproblem contains non-linear constraints */
void SCIPbendersSetSubproblemIsNonlinear(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isnonlinear         /**< flag to indicate whether the subproblem contains non-linear constraints */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   if( isnonlinear && !benders->subprobisnonlinear[probnumber] )
      benders->nnonlinearsubprobs++;
   else if( !isnonlinear && benders->subprobisnonlinear[probnumber] )
      benders->nnonlinearsubprobs--;

   benders->subprobisnonlinear[probnumber] = isnonlinear;

   assert(benders->nnonlinearsubprobs >= 0 && benders->nnonlinearsubprobs <= benders->nsubproblems);
}

/** returns whether the subproblem contains non-linear constraints. */
SCIP_Bool SCIPbendersSubproblemIsNonlinear(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobisnonlinear[probnumber];
}

/** returns the number of subproblems that contain non-linear constraints  */
int SCIPbendersGetNNonlinearSubproblems(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nnonlinearsubprobs;
}

/** sets the flag indicating whether the master problem contains non-linear constraints */
void SCIPbendersSetMasterIsNonlinear(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_Bool             isnonlinear         /**< flag to indicate whether the subproblem contains non-linear constraints */
   )
{
   assert(benders != NULL);

   benders->masterisnonlinear = isnonlinear;
}

/** returns whether the master problem contains non-linear constraints. */
SCIP_Bool SCIPbendersMasterIsNonlinear(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->masterisnonlinear;
}

/** returns the flag indicating that Benders' decomposition is in a cut strengthening round */
SCIP_Bool SCIPbendersInStrengthenRound(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->strengthenround;
}

/** changes all of the master problem variables in the given subproblem to continuous. */
SCIP_RETCODE SCIPbendersChgMastervarsToCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   probnumber          /**< the subproblem number */
   )
{
   SCIP* subproblem;
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int chgvarscount;
   int origintvars;
   int i;
   SCIP_Bool infeasible;

   assert(benders != NULL);
   assert(set != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   subproblem = SCIPbendersSubproblem(benders, probnumber);
   assert(subproblem != NULL);

   /* only set the master problem variable to continuous if they have not already been changed. */
   if( !SCIPbendersGetMastervarsCont(benders, probnumber) )
   {
      SCIP_VAR* mastervar;

      /* retrieving the variable data */
      SCIP_CALL( SCIPgetVarsData(subproblem, &vars, NULL, &nbinvars, &nintvars, &nimplvars, NULL) );

      origintvars = nbinvars + nintvars + nimplvars;

      chgvarscount = 0;

      /* looping over all integer variables to change the master variables to continuous */
      i = 0;
      while( i < nbinvars + nintvars + nimplvars )
      {
         SCIP_CALL( SCIPbendersGetVar(benders, set, vars[i], &mastervar, -1) );

         if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS && mastervar != NULL )
         {
            /* changing the type of the subproblem variable corresponding to mastervar to CONTINUOUS */
            SCIP_CALL( SCIPchgVarType(subproblem, vars[i], SCIP_VARTYPE_CONTINUOUS, &infeasible) );

            assert(!infeasible);

            chgvarscount++;
            SCIP_CALL( SCIPgetVarsData(subproblem, NULL, NULL, &nbinvars, &nintvars, &nimplvars, NULL) );
         }
         else
            i++;
      }

      /* if all of the integer variables have been changed to continuous, then the subproblem could now be a convex
       * problem. This must be checked and if TRUE, then the LP subproblem is initialised and then put into probing
       * mode
       */
      if( chgvarscount > 0 && chgvarscount == origintvars )
      {
         /* checking the convexity of the subproblem */
         SCIP_CALL( checkSubproblemConvexity(benders, set, probnumber) );

         /* if the subproblem has convex constraints and continuous variables, then it is initialised and put into
          * probing mode
          */
         if( SCIPbendersGetSubproblemType(benders, probnumber) == SCIP_BENDERSSUBTYPE_CONVEXCONT )
         {
            SCIP_CALL( initialiseLPSubproblem(benders, set, probnumber) );
         }
      }

      SCIP_CALL( SCIPbendersSetMastervarsCont(benders, probnumber, TRUE) );
   }

   return SCIP_OKAY;
}

/** sets the subproblem setup flag */
void SCIPbendersSetSubproblemIsSetup(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             issetup             /**< flag to indicate whether the subproblem has been setup */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   benders->subprobsetup[probnumber] = issetup;
}

/** returns the subproblem setup flag */
SCIP_Bool SCIPbendersSubproblemIsSetup(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobsetup[probnumber];
}

/** sets the independent subproblem flag */
void SCIPbendersSetSubproblemIsIndependent(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             isindep             /**< flag to indicate whether the subproblem is independent */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* if the user has defined solving or freeing functions, then it is not possible to declare a subproblem as
    * independent. This is because declaring a subproblem as independent changes the solving loop, so it would change
    * the expected behaviour of the user defined plugin. If a user calls this function, then an error will be returned.
    */
   if( benders->benderssolvesubconvex != NULL || benders->benderssolvesub != NULL || benders->bendersfreesub != NULL )
   {
      SCIPerrorMessage("The user has defined either bendersSolvesubconvex%s, bendersSolvesub%s or bendersFreesub%s. "
         "Thus, it is not possible to declare the independence of a subproblem.\n", benders->name, benders->name,
         benders->name);
      SCIPABORT();
   }
   else
   {
      SCIP_Bool activesubprob;

      /* if the active status of the subproblem changes, then we must update the activesubprobs counter */
      activesubprob = subproblemIsActive(benders, probnumber);

      benders->indepsubprob[probnumber] = isindep;

      /* updating the activesubprobs counter */
      if( activesubprob && !subproblemIsActive(benders, probnumber) )
         benders->nactivesubprobs--;
      else if( !activesubprob && subproblemIsActive(benders, probnumber) )
         benders->nactivesubprobs++;

      assert(benders->nactivesubprobs >= 0 && benders->nactivesubprobs <= SCIPbendersGetNSubproblems(benders));
   }
}

/** returns whether the subproblem is independent */
SCIP_Bool SCIPbendersSubproblemIsIndependent(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->indepsubprob[probnumber];
}

/** Sets whether the subproblem is enabled or disabled. A subproblem is disabled if it has been merged into the master
 *  problem.
 */
void SCIPbendersSetSubproblemEnabled(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             enabled             /**< flag to indicate whether the subproblem is enabled */
   )
{
   SCIP_Bool activesubprob;

   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* if the active status of the subproblem changes, then we must update the activesubprobs counter */
   activesubprob = subproblemIsActive(benders, probnumber);

   benders->subprobenabled[probnumber] = enabled;

   /* updating the activesubprobs counter */
   if( activesubprob && !subproblemIsActive(benders, probnumber) )
      benders->nactivesubprobs--;
   else if( !activesubprob && subproblemIsActive(benders, probnumber) )
      benders->nactivesubprobs++;

   assert(benders->nactivesubprobs >= 0 && benders->nactivesubprobs <= SCIPbendersGetNSubproblems(benders));
}

/** returns whether the subproblem is enabled, i.e. the subproblem is still solved in the solving loop. */
SCIP_Bool SCIPbendersSubproblemIsEnabled(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subprobenabled[probnumber];
}

/** sets a flag to indicate whether the master variables are all set to continuous */
SCIP_RETCODE SCIPbendersSetMastervarsCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Bool             arecont             /**< flag to indicate whether the master problem variables are continuous */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   /* if the master variables were all continuous and now are not, then the subproblem must exit probing mode and be
    * changed to non-LP subproblem */
   if( benders->mastervarscont[probnumber] && !arecont )
   {
      SCIP_BENDERSSUBTYPE subtype;

      if( SCIPinProbing(SCIPbendersSubproblem(benders, probnumber)) )
      {
         SCIP_CALL( SCIPendProbing(SCIPbendersSubproblem(benders, probnumber)) );
      }

      subtype = SCIPbendersGetSubproblemType(benders, probnumber);
      assert(subtype == SCIP_BENDERSSUBTYPE_CONVEXCONT || subtype == SCIP_BENDERSSUBTYPE_NONCONVEXCONT);

      if( subtype == SCIP_BENDERSSUBTYPE_CONVEXCONT )
         SCIPbendersSetSubproblemType(benders, probnumber, SCIP_BENDERSSUBTYPE_CONVEXDIS);
      else if( subtype == SCIP_BENDERSSUBTYPE_NONCONVEXCONT )
         SCIPbendersSetSubproblemType(benders, probnumber, SCIP_BENDERSSUBTYPE_NONCONVEXDIS);
   }

   benders->mastervarscont[probnumber] = arecont;

   return SCIP_OKAY;
}

/** returns whether the master variables are all set to continuous */
SCIP_Bool SCIPbendersGetMastervarsCont(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->mastervarscont[probnumber];
}

/** returns the number of cuts that have been transferred from sub SCIPs to the master SCIP */
int SCIPbendersGetNTransferredCuts(
   SCIP_BENDERS*         benders             /**< the Benders' decomposition data structure */
   )
{
   assert(benders != NULL);

   return benders->ntransferred;
}

/** updates the lower bound for the subproblem. If the lower bound is not greater than the previously stored lowerbound,
 *  then no update occurs.
 */
void SCIPbendersUpdateSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber,         /**< the subproblem number */
   SCIP_Real             lowerbound          /**< the lower bound */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   if( EPSGE(lowerbound, benders->subproblowerbound[probnumber], 1e-06) )
      benders->subproblowerbound[probnumber] = lowerbound;
   else
   {
      SCIPdebugMessage("The lowerbound %g for subproblem %d is less than the currently stored lower bound %g\n",
         lowerbound, probnumber, benders->subproblowerbound[probnumber]);
   }
}

/** returns the stored lower bound for the given subproblem */
SCIP_Real SCIPbendersGetSubproblemLowerbound(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   probnumber          /**< the subproblem number */
   )
{
   assert(benders != NULL);
   assert(probnumber >= 0 && probnumber < SCIPbendersGetNSubproblems(benders));

   return benders->subproblowerbound[probnumber];
}

/** returns the number of cuts that have been added for storage */
int SCIPbendersGetNStoredCuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nstoredcuts;
}

/** returns the cuts that have been stored for transfer */
SCIP_RETCODE SCIPbendersGetStoredCutData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   int                   cutidx,             /**< the index for the cut data that is requested */
   SCIP_VAR***           vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real**           vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   SCIP_Real*            rhs,                /**< the right hand side of the cut */
   int*                  nvars               /**< the number of variables with non-zero coefficients in the cut */
   )
{
   assert(benders != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nvars != NULL);
   assert(cutidx >= 0 && cutidx < benders->nstoredcuts);

   (*vars) = benders->storedcuts[cutidx]->vars;
   (*vals) = benders->storedcuts[cutidx]->vals;
   (*lhs) = benders->storedcuts[cutidx]->lhs;
   (*rhs) = benders->storedcuts[cutidx]->rhs;
   (*nvars) = benders->storedcuts[cutidx]->nvars;

   return SCIP_OKAY;
}

/** returns the original problem data for the cuts that have been added by the Benders' cut plugin. The stored
 *  variables and values will populate the input vars and vals arrays. Thus, memory must be allocated for the vars and
 *  vals arrays
 */
SCIP_RETCODE SCIPbendersGetStoredCutOrigData(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition cut */
   int                   cutidx,             /**< the index for the cut data that is requested */
   SCIP_VAR***           vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real**           vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real*            lhs,                /**< the left hand side of the cut */
   SCIP_Real*            rhs,                /**< the right hand side of the cut */
   int*                  nvars,              /**< the number of variables with non-zero coefficients in the cut */
   int                   varssize            /**< the available slots in the array */
   )
{
   SCIP_VAR* origvar;
   SCIP_Real scalar;
   SCIP_Real constant;
   int i;

   assert(benders != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(lhs != NULL);
   assert(rhs != NULL);
   assert(nvars != NULL);
   assert(cutidx >= 0 && cutidx < benders->nstoredcuts);

   (*lhs) = benders->storedcuts[cutidx]->lhs;
   (*rhs) = benders->storedcuts[cutidx]->rhs;
   (*nvars) = benders->storedcuts[cutidx]->nvars;

   /* if there are enough slots, then store the cut variables and values */
   if( varssize >= *nvars )
   {
      for( i = 0; i < *nvars; i++ )
      {
         /* getting the original variable for the transformed variable */
         origvar = benders->storedcuts[cutidx]->vars[i];
         scalar = 1.0;
         constant = 0.0;
         SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

         (*vars)[i] = origvar;
         (*vals)[i] = benders->storedcuts[cutidx]->vals[i];
      }
   }

   return SCIP_OKAY;
}

/** adds the data for the generated cuts to the Benders' cut storage */
SCIP_RETCODE SCIPbendersStoreCut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition cut */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< the variables that have non-zero coefficients in the cut */
   SCIP_Real*            vals,               /**< the coefficients of the variables in the cut */
   SCIP_Real             lhs,                /**< the left hand side of the cut */
   SCIP_Real             rhs,                /**< the right hand side of the cut */
   int                   nvars               /**< the number of variables with non-zero coefficients in the cut */
   )
{
   SCIP_BENDERSCUTCUT* cut;

   assert(benders != NULL);
   assert(set != NULL);
   assert(vars != NULL);
   assert(vals != NULL);

   /* allocating the block memory for the cut storage */
   SCIP_CALL( SCIPallocBlockMemory(set->scip, &cut) );

   /* storing the cut data */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(set->scip, &cut->vars, vars, nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(set->scip, &cut->vals, vals, nvars) );
   cut->lhs = lhs;
   cut->rhs = rhs;
   cut->nvars = nvars;

   /* ensuring the required memory is available for the stored cuts array */
   if( benders->storedcutssize < benders->nstoredcuts + 1 )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, benders->nstoredcuts + 1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(SCIPblkmem(set->scip), &benders->storedcuts,
            benders->storedcutssize, newsize) );

      benders->storedcutssize = newsize;
   }
   assert(benders->storedcutssize >= benders->nstoredcuts + 1);

   /* adding the cuts to the Benders' cut storage */
   benders->storedcuts[benders->nstoredcuts] = cut;
   benders->nstoredcuts++;

   return SCIP_OKAY;
}

/** sets the sorted flags in the Benders' decomposition */
void SCIPbendersSetBenderscutsSorted(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_Bool             sorted              /**< the value to set the sorted flag to */
   )
{
   assert(benders != NULL);

   benders->benderscutssorted = sorted;
   benders->benderscutsnamessorted = sorted;
}

/** inserts a Benders' cut into the Benders' cuts list */
SCIP_RETCODE SCIPbendersIncludeBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERSCUT*      benderscut          /**< Benders' cut */
   )
{
   assert(benders != NULL);
   assert(benderscut != NULL);

   if( benders->nbenderscuts >= benders->benderscutssize )
   {
      benders->benderscutssize = SCIPsetCalcMemGrowSize(set, benders->nbenderscuts+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&benders->benderscuts, benders->benderscutssize) );
   }
   assert(benders->nbenderscuts < benders->benderscutssize);

   benders->benderscuts[benders->nbenderscuts] = benderscut;
   benders->nbenderscuts++;
   benders->benderscutssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the Benders' cut of the given name, or NULL if not existing */
SCIP_BENDERSCUT* SCIPfindBenderscut(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   const char*           name                /**< name of Benderscut' decomposition */
   )
{
   int i;

   assert(benders != NULL);
   assert(name != NULL);

   for( i = 0; i < benders->nbenderscuts; i++ )
   {
      if( strcmp(SCIPbenderscutGetName(benders->benderscuts[i]), name) == 0 )
         return benders->benderscuts[i];
   }

   return NULL;
}

/** returns the array of currently available Benders' cuts; active Benders' decomposition are in the first slots of
 * the array
 */
SCIP_BENDERSCUT** SCIPbendersGetBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutssorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutComp, benders->nbenderscuts);
      benders->benderscutssorted = TRUE;
      benders->benderscutsnamessorted = FALSE;
   }

   return benders->benderscuts;
}

/** returns the number of currently available Benders' cuts */
int SCIPbendersGetNBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   return benders->nbenderscuts;
}

/** sets the priority of a Benders' decomposition */
SCIP_RETCODE SCIPbendersSetBenderscutPriority(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' cut */
   int                   priority            /**< new priority of the Benders' decomposition */
   )
{
   assert(benders != NULL);
   assert(benderscut != NULL);

   benderscut->priority = priority;
   benders->benderscutssorted = FALSE;

   return SCIP_OKAY;
}

/** sorts Benders' decomposition cuts by priorities */
void SCIPbendersSortBenderscuts(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutssorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutComp, benders->nbenderscuts);
      benders->benderscutssorted = TRUE;
      benders->benderscutsnamessorted = FALSE;
   }
}

/** sorts Benders' decomposition cuts by name */
void SCIPbendersSortBenderscutsName(
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   assert(benders != NULL);

   if( !benders->benderscutsnamessorted )
   {
      SCIPsortPtr((void**)benders->benderscuts, SCIPbenderscutCompName, benders->nbenderscuts);
      benders->benderscutssorted = FALSE;
      benders->benderscutsnamessorted = TRUE;
   }
}
