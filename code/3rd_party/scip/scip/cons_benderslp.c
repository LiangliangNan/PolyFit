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

/**@file   cons_benderslp.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for benderslp decomposition
 * @author Stephen J. Maher
 *
 * Two constraint handlers are implemented for the generation of Benders' decomposition cuts. When included in a
 * problem, the Benders' decomposition constraint handlers generate cuts during the enforcement of LP and relaxation
 * solutions. Additionally, Benders' decomposition cuts can be generated when checking the feasibility of solutions with
 * respect to the subproblem constraints.
 *
 * This constraint handler has an enforcement priority that is greater than the integer constraint handler. This means
 * that all LP solutions will be first checked for feasibility with respect to the Benders' decomposition second stage
 * constraints before performing an integrality check. This is part of a multi-phase approach for solving mixed integer
 * programs by Benders' decomposition.
 *
 * A parameter is available to control the depth at which the non-integer LP solution are enforced by solving the
 * Benders' decomposition subproblems. This parameter is set to 0 by default, indicating that non-integer LP solutions
 * are enforced only at the root node.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/cons_benderslp.h"
#include "scip/cons_benders.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "benderslp"
#define CONSHDLR_DESC          "constraint handler for Benders' Decomposition to separate LP solutions"
#define CONSHDLR_ENFOPRIORITY  10000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY 10000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */


#define DEFAULT_CONSBENDERSLP_MAXDEPTH 0/**< depth at which Benders' decomposition cuts are generated from the LP
                                         *   solution (-1: always, 0: only at root) */
#define DEFAULT_CONSBENDERSLP_FREQ     0/**< the depth frequency for generating LP cuts after the max depth is reached (0: never, 1: all nodes, ...) */
#define DEFAULT_CONSBENDERSLP_STALLLIMIT 100/**< the number of nodes processed without a dual bound improvement before enforcing the LP relaxation, 0: no stall count applied */
#define DEFAULT_CONSBENDERSLP_ITERLIMIT 100 /**< after the root node, only iterlimit fractional LP solutions are used at each node to generate Benders' decomposition cuts. */
#define DEFAULT_ACTIVE            FALSE /**< is the constraint handler active? */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   /* parameters for controlling the two-phase method for Benders' decomposition */
   int                   maxdepth;           /**< the maximum depth at which Benders' cuts are generated from the LP */
   int                   freq;               /**< the depth frequency of generating LP cuts after the max depth is reached */
   SCIP_Bool             active;             /**< is the constraint handler active? */

   /* variable used to control the behaviour of the two-phase method for Benders' decomposition */
   SCIP_Longint          ncallsnode;         /**< the number of calls at the current node */
   SCIP_NODE*            currnode;           /**< the current node */
   SCIP_Real             prevbound;          /**< the previous dual bound */
   int                   iterlimit;          /**< the iteration limit for the first phase of the two-phase method at a node lower than the root. */
   int                   stallcount;         /**< the number of nodes processed since the last lower bound increase */
   int                   stalllimit;         /**< the number of nodes processed without bound improvement before enforcing the LP relaxation */
};


/*
 * Local methods
 */


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBenderslp)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeConshdlrBenderslp(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeBenderslp)
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



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBenderslp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* updating the stall count. If the bound has improved since the last call, then the stall count is set to zero */
   conshdlrdata->stallcount++;
   if( SCIPisLT(scip, conshdlrdata->prevbound, SCIPgetLowerbound(scip)) )
      conshdlrdata->stallcount = 0;

   conshdlrdata->prevbound = SCIPgetLowerbound(scip);
   conshdlrdata->ncallsnode++;

   /* if a new node is being processed, then the iteration counts are reset. */
   if( conshdlrdata->currnode != SCIPgetCurrentNode(scip) )
   {
      conshdlrdata->currnode = SCIPgetCurrentNode(scip);
      conshdlrdata->ncallsnode = 0;
   }

   /* the result is initially set to FEASIBLE. If the two-phase method is not executed, then the result will remain as
    * FEASIBLE. The actual feasibility of the Benders' decomposition subproblems is checked in cons_benders.
   */
   (*result) = SCIP_FEASIBLE;

   /* only check the Benders' decomposition subproblems for fractional LP solutions if the two-phase method is activated */
   if( !conshdlrdata->active )
      return SCIP_OKAY;

   /* checking whether the two-phase method is performed.
    * - If a maxdepth is specified
    *   - current depth is checked
    *   - the frequency is checked, if a frequency is specified
    *   - the stalllimit is checked if a stalllimit is specified.
    */
   if( conshdlrdata->maxdepth >= 0 && SCIPgetDepth(scip) > conshdlrdata->maxdepth
         && (conshdlrdata->freq == 0 || SCIPgetDepth(scip) % conshdlrdata->freq != 0)
         && (conshdlrdata->stalllimit == 0 || conshdlrdata->stallcount < conshdlrdata->stalllimit) )
      return SCIP_OKAY;

   /* if an iteration limit is specified, then this is imposed at nodes after the root node */
   if( SCIPgetDepth(scip) > 0 && conshdlrdata->ncallsnode >= conshdlrdata->iterlimit )
      return SCIP_OKAY;

   /* the two phase method is only performed at the root node for sub-SCIPs */
   if( SCIPgetSubscipDepth(scip) > 0 && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   /* checking the Benders' decomposition subproblems for feasibility. */
   SCIP_CALL( SCIPconsBendersEnforceSolution(scip, NULL, conshdlr, result, SCIP_BENDERSENFOTYPE_LP, FALSE) );

   /* if the stalllimit is exceeded and the subproblem were checked, then the stall count is reset to zero */
   if( conshdlrdata->stallcount >= conshdlrdata->stalllimit )
      conshdlrdata->stallcount = 0;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxBenderslp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   if( !conshdlrdata->active || (conshdlrdata->maxdepth >= 0 && SCIPgetDepth(scip) > conshdlrdata->maxdepth) )
      (*result) = SCIP_FEASIBLE;
   else
      SCIP_CALL( SCIPconsBendersEnforceSolution(scip, sol, conshdlr, result, SCIP_BENDERSENFOTYPE_RELAX, FALSE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBenderslp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   if( !conshdlrdata->active || (conshdlrdata->maxdepth >= 0 && SCIPgetDepth(scip) > conshdlrdata->maxdepth) )
      (*result) = SCIP_FEASIBLE;
   else
      SCIP_CALL( SCIPconsBendersEnforceSolution(scip, NULL, conshdlr, result, SCIP_BENDERSENFOTYPE_PSEUDO, FALSE) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions.
 *  The feasibility check for Benders' decomposition is performed in cons_benders. As such, it is redundant to perform
 *  the feasibility check here. As such, the solution is flagged as feasible, which will then be corrected in
 *  cons_benders if the solution is infeasible with respect to the second stage constraints
 */
static
SCIP_DECL_CONSCHECK(consCheckBenderslp)
{  /*lint --e{715}*/
   (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBenderslp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}



/*
 * constraint specific interface methods
 */

/** creates the handler for executing the Benders' decomposition subproblem solve on fractional LP solution and
 * includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBenderslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   /* create benderslp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);
   conshdlrdata->prevbound = -SCIPinfinity(scip);

   conshdlr = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpBenderslp, consEnfopsBenderslp, consCheckBenderslp, consLockBenderslp,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyBenderslp, NULL) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeBenderslp) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxBenderslp) );

   /* add Benders' decomposition LP constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxdepth",
         "depth at which Benders' decomposition cuts are generated from the LP solution (-1: always, 0: only at root)",
         &conshdlrdata->maxdepth, TRUE, DEFAULT_CONSBENDERSLP_MAXDEPTH, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/depthfreq",
         "the depth frequency for generating LP cuts after the max depth is reached (0: never, 1: all nodes, ...)",
         &conshdlrdata->freq, TRUE, DEFAULT_CONSBENDERSLP_FREQ, 0, SCIP_MAXTREEDEPTH, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/stalllimit",
         "the number of nodes processed without a dual bound improvement before enforcing the LP relaxation, 0: no stall count applied",
         &conshdlrdata->stalllimit, TRUE, DEFAULT_CONSBENDERSLP_STALLLIMIT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/iterlimit",
         "after the root node, only iterlimit fractional LP solutions are used at each node to generate Benders' decomposition cuts.",
         &conshdlrdata->iterlimit, TRUE, DEFAULT_CONSBENDERSLP_ITERLIMIT, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/active", "is the Benders' decomposition LP cut constraint handler active?",
         &conshdlrdata->active, FALSE, DEFAULT_ACTIVE, NULL, NULL));

   conshdlrdata->stallcount = 0;
   return SCIP_OKAY;
}
