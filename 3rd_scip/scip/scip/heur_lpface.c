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

/**@file   heur_lpface.c
 * @brief  lpface primal heuristic that searches the optimal LP face inside a sub-MIP
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/heur_lpface.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "lpface"
#define HEUR_DESC             "LNS heuristic that searches the optimal LP face inside a sub-MIP"
#define HEUR_DISPCHAR         '_'
#define HEUR_PRIORITY         -1104000
#define HEUR_FREQ             15
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL         /**< maximum number of nodes to regard in the subproblem                   */
#define DEFAULT_MINNODES      50LL           /**< minimum number of nodes to regard in the subproblem                   */
#define DEFAULT_MINFIXINGRATE 0.1            /**< required percentage of fixed integer variables in sub-MIP to run */
#define DEFAULT_NODESOFS      200LL          /**< number of nodes added to the contingent of the total nodes            */
#define DEFAULT_NODESQUOT     0.1            /**< subproblem nodes in relation to nodes of the original problem         */
#define DEFAULT_LPLIMFAC      2.0            /**< factor by which the limit on the number of LP depends on the node limit */
#define DEFAULT_USELPROWS     TRUE           /**< should subproblem be created out of the rows in the LP rows,
                                              *   otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?                                     */
#define DEFAULT_DUALBASISEQUATIONS FALSE     /**< should the dually nonbasic rows be turned into equations?        */
#define DEFAULT_KEEPSUBSCIP   FALSE          /**< should the heuristic continue solving the same sub-SCIP? */
#define DEFAULT_MINPATHLEN        5          /**< the minimum active search tree path length along which the lower bound
                                              *   hasn't changed before heuristic becomes active */
/* event handler properties */
#define EVENTHDLR_NAME         "Lpface"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"

/*
 * Data structures
 */

/** data structure to keep sub-SCIP across runs */
struct SubscipData
{
   SCIP*                 subscip;            /**< pointer to store sub-SCIP data structure */
   SCIP_VAR**            subvars;            /**< array of variables of the sub-problem */
   int                   nsubvars;           /**< number of sub-problem variables */
   SCIP_Real             objbound;           /**< lower bound on objective for when sub SCIP was created */
};
typedef struct SubscipData SUBSCIPDATA;

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem               */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem               */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes        */
   SCIP_Longint          usednodes;          /**< nodes already used by lpface in earlier calls                  */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem     */

   unsigned int          nfailures;          /**< number of failures since last successful call                     */
   SCIP_Longint          nextnodenumber;     /**< number of nodes at which lpface should be called the next time */
   SCIP_Real             lastlpobjinfeas;    /**< last LP objective where the sub-MIP was run to proven infeasibility */
   SCIP_Real             minfixingrate;      /**< required percentage of fixed integer variables in sub-MIP to run     */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?      */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?                                     */
   SCIP_Bool             dualbasisequations; /**< should the dually nonbasic rows be turned into equations?        */
   SCIP_Bool             keepsubscip;        /**< should the heuristic continue solving the same sub-SCIP? */
   char                  subscipobjective;   /**< objective function in the sub-SCIP: (z)ero, (r)oot-LP-difference,
                                              *   (i)nference, LP (f)ractionality, (o)riginal */

   SCIP_STATUS           submipstatus;       /**< return status of the sub-MIP */
   SCIP_Longint          submipnlpiters;     /**< number of LP iterations of the sub-MIP */
   SCIP_Real             submippresoltime;   /**< time required to presolve the sub-MIP */
   int                   nvarsfixed;         /**< the number of fixed variables by the heuristic */
   int                   minpathlen;         /**< the minimum active search tree path length along which the lower bound
                                              *   hasn't changed before heuristic becomes active */
   SUBSCIPDATA*          subscipdata;        /**< sub-SCIP data structure */
};

/*
 * Local methods
 */

/** fixes variables of the subproblem considering their reduced costs */
static
SCIP_RETCODE fixVariables(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            success             /**< pointer to store whether enough integer variables were fixed */
   )
{
   SCIP_VAR** vars;                          /* original scip variables                */
   SCIP_Real fixingrate;                     /* percentage of variables that are fixed */
   int nvars;
   int nbinvars;
   int nintvars;
   int i;
   int fixingcounter;

   /* get required data of the main scip problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   fixingcounter = 0;

   assert(nvars >= nbinvars + nintvars);

   /* loop over problem variables and fix all with nonzero reduced costs to their solution value */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real solval;
      SCIP_COL* col;
      SCIP_Real redcost;
      SCIP_Real lbglobal;
      SCIP_Real ubglobal;
      SCIP_VAR* var;

      var = vars[i];

      /* skip non-column variables */
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
         continue;

      solval = SCIPgetSolVal(scip, NULL, var);
      col = SCIPvarGetCol(vars[i]);
      assert(col != NULL);
      redcost = SCIPgetColRedcost(scip, col);
      lbglobal = SCIPvarGetLbGlobal(var);
      ubglobal = SCIPvarGetUbGlobal(var);

      /* fix the variable to its solution value if variable is nonbasic (i.e., at one of its bounds)
       *  with nonzero reduced costs
       */
      if( ! SCIPisDualfeasZero(scip, redcost) )
      {
         /* fix variable based on reduced cost information, respecting global bounds */
         if( (redcost > 0 && SCIPisFeasEQ(scip, solval, lbglobal)) ||
             (redcost < 0 && SCIPisFeasEQ(scip, solval, ubglobal)) )
         {
            SCIPdebugMsg(scip, "Fixing variable <%s> (obj: %g), local bounds [%.1g, %.1g], redcost %9.5g, LP sol val %9.5g\n",
                  SCIPvarGetName(var), SCIPvarGetObj(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), redcost, solval);
            assert(! SCIPisInfinity(scip, solval));
            assert(! SCIPisInfinity(scip, -solval));
            SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], solval) );
            SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], solval) );
            if( SCIPvarIsIntegral(var) )
               ++fixingcounter;
         }

      }
   }

   fixingrate = (SCIP_Real)fixingcounter / (SCIP_Real)(MAX(nbinvars + nintvars, 1));
   heurdata->nvarsfixed = fixingcounter;

   /* if all variables were fixed or amount of fixed variables is insufficient, skip residual part of
    * subproblem creation and abort immediately
    */
   *success = (fixingcounter < nvars && fixingrate >= heurdata->minfixingrate);

   SCIPdebugMsg(scip, " LP face heuristic fixed %senough variables (%d out of %d)\n",
      *success ? "": "not ", fixingcounter, nvars);

   return SCIP_OKAY;
}

/** creates the rows of the subproblem */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_Bool             dualbasisequations  /**< should the dually nonbasic rows be turned into equations? */
   )
{
   SCIP_ROW** rows;                          /* original scip rows                       */

   int nrows;
   int i;

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
      SCIP_VAR** consvars;                   /* new constraint's variables               */
      SCIP_COL** cols;                       /* original row's columns                   */
      SCIP_CONS* cons;                       /* new constraint                           */

      SCIP_Real* vals;                       /* variables' coefficient values of the row */
      SCIP_Real constant;                    /* constant added to the row                */
      SCIP_Real lhs;                         /* left hand side of the row                */
      SCIP_Real rhs;                         /* left right side of the row               */
      SCIP_Real dualsol;
      SCIP_Real rowsolactivity;
      int j;
      int nnonz;

      /* ignore rows that are only locally valid */
      if( SCIProwIsLocal(rows[i]) )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      /* only subtract constant if left hand side is not infinite */
      lhs = SCIProwGetLhs(rows[i]);
      if( ! SCIPisInfinity(scip, -lhs) )
         lhs -= constant;

      /* only subtract constant if right hand side is not infinite */
      rhs = SCIProwGetRhs(rows[i]);
      if( ! SCIPisInfinity(scip, rhs) )
         rhs -= constant;

      assert(lhs <= rhs);

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      dualsol = SCIProwGetDualsol(rows[i]);
      rowsolactivity = SCIPgetRowActivity(scip, rows[i]);

      /* transform into equation if the row is sharp and has a nonzero dual solution */
      if( dualbasisequations && ! SCIPisDualfeasZero(scip, dualsol) )
      {
         if( dualsol > 0.0 && SCIPisFeasEQ(scip, rowsolactivity, lhs) )
            rhs = lhs;
         else if( dualsol < 0.0 && SCIPisFeasEQ(scip, rowsolactivity, rhs) )
            lhs = rhs;
      }

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** creates the LP face subproblem by fixing nonbasic variables with nonzero reduced costs */
static
SCIP_RETCODE setupSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully */
   )
{
   SCIP_VAR** vars = SCIPgetVars(scip);
   int nvars = SCIPgetNVars(scip);
   SCIP_Real lowerbound;
   SCIP_CONS* origobjcons;
   int i;
#ifndef NDEBUG
   int nobjvars = 0;
#endif

   /* fix variables in subproblem with nonzero reduced costs */
   SCIP_CALL( fixVariables(scip, subscip, subvars, heurdata, success) );

   if( ! (*success) )
      return SCIP_OKAY;

   /* we copy the rows of the LP, if enough variables could be fixed and we work on the MIP relaxation of the problem */
   if( *success && heurdata->uselprows )
   {
      SCIP_CALL( createRows(scip, subscip, subvars, heurdata->dualbasisequations) );
   }

   /* add an equation that the objective function must be equal to the lower bound */
   lowerbound = SCIPgetLowerbound(scip);

   SCIP_CALL( SCIPcreateConsLinear(subscip, &origobjcons, "objbound_of_origscip", 0, NULL, NULL, lowerbound, lowerbound,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   for( i = 0; i < nvars; ++i)
   {
      if( ! SCIPisZero(subscip, SCIPvarGetObj(vars[i])) )
      {
         SCIP_CALL( SCIPaddCoefLinear(subscip, origobjcons, subvars[i], SCIPvarGetObj(vars[i])) );
#ifndef NDEBUG
         nobjvars++;
#endif
      }
   }
   assert(nobjvars == SCIPgetNObjVars(scip));

   SCIP_CALL( SCIPaddCons(subscip, origobjcons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &origobjcons) );

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEUR*            heur,               /**< lpface heuristic structure */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   int*                  solindex,           /**< pointer to store index of the solution */
   SCIP_Bool*            success             /**< pointer to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_Bool printreason;
   SCIP_Bool completely;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );
   *solindex = SCIPsolGetIndex(newsol);

#ifdef SCIP_DEBUG
   printreason = TRUE;
   completely = TRUE;
   SCIPdebugMsg(scip, "trying to transfer LP face solution with solution value %16.9g to main problem\n",
      SCIPretransformObj(scip, SCIPgetSolTransObj(scip, newsol)));
#else
   printreason = FALSE;
   completely = FALSE;
#endif

   /* try to add new solution to scip and free it immediately */
   *success = FALSE;
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, printreason, completely, TRUE, TRUE, TRUE, success) );

   SCIPdebugMsg(scip, "Transfer was %s successful\n", *success ? "" : "not");

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** updates heurdata after an unsuccessful run of lpface */
static
void updateFailureStatistic(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   /* increase number of failures, calculate next node at which lpface should be called and update actual solutions */
   heurdata->nfailures++;
   heurdata->nextnodenumber = (heurdata->nfailures <= 25
      ? SCIPgetNNodes(scip) + 100*(2LL << heurdata->nfailures) /*lint !e703*/
      : SCIP_LONGINT_MAX);
}

/** calculate a node limit based on node limiting parameters of the heuristic */
static
SCIP_Longint calcNodeLimit(
   SCIP*                 scip,               /**< (original) SCIP data structure */
   SCIP_HEUR*            heur,               /**< LP face heuristic */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   SCIP_Longint nodelimit;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nodelimit = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* count the setup costs for the sub-MIP as 100 nodes */
   nodelimit -= 100 * SCIPheurGetNCalls(heur);

   /* add the offset */
   nodelimit += heurdata->nodesofs;

   /* subtract previously used nodes */
   nodelimit -= heurdata->usednodes;

   /* do not use more than the maximum number of allowed nodes in one run */
   nodelimit = MIN(nodelimit, heurdata->maxnodes);

   /* if the subscip has been kept from a previous run, add the number of already processed nodes */
   if( heurdata->subscipdata->subscip != NULL )
      nodelimit += SCIPgetNNodes(heurdata->subscipdata->subscip);

   return nodelimit;
}

/** sets node, time, and memory limit according to the parameter settings of the heuristic */
static
SCIP_RETCODE setSubscipLimits(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< data structure of the sub-problem */
   SCIP_HEUR*            heur,               /**< LP face heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_Bool*            success             /**< did we successfully set all parameters up? */
   )
{
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Longint nodelimit;

   *success = TRUE;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   if( ! SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( ! SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= 2.0 * SCIPgetMemExternEstim(scip) / 1048576.0 )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* calculate node limit for the subproblem */
   nodelimit = calcNodeLimit(scip, heur, heurdata);

   /* we should have aborted the sub-SCIP procedure earlier if no additional nodes are allowed
    * with the current parameter settings
    */
   assert(nodelimit > SCIPgetNNodes(subscip));

   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );
   heurdata->nodelimit = nodelimit;

   /* set also the other two limits */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   return SCIP_OKAY;
}

/** sets all one-time parameter settings like search strategy, but no limits */
static
SCIP_RETCODE setSubscipParameters(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip             /**< data structure of the sub-problem */
   )
{

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* for debugging lpface, enable MIP output */
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1) );
#endif

   /* disable statistic timing inside sub SCIP */
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable expensive cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use restart depth first node selection */
   if( SCIPfindNodesel(subscip, "restartdfs") != NULL && ! SCIPisParamFixed(subscip, "nodeselection/restartdfs/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/restartdfs/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && ! SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && ! SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 500) );
   }

   /* enable conflict analysis, disable analysis of boundexceeding LPs, and restrict conflict pool */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetCharParam(subscip, "conflict/useboundlp", 'o') );
   }
   if( !SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }

   return SCIP_OKAY;
}

/** reset the sub-SCIP data to its default values */
static
void subscipdataReset(
   SUBSCIPDATA*          subscipdata         /**< data structure of the sub-problem */
   )
{
   subscipdata->subscip = NULL;
   subscipdata->subvars = NULL;
   subscipdata->nsubvars = 0;
   subscipdata->objbound = SCIP_INVALID;
}

/** free the stored sub-SCIP information */
static
SCIP_RETCODE subscipdataFreeSubscip(
   SCIP*                 scip,               /**< original SCIP data structure */
   SUBSCIPDATA*          subscipdata         /**< data structure of the sub-problem */
   )
{
   assert(subscipdata != NULL);

   /* free the subscipdata's scip */
   if( subscipdata->subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&subscipdata->subscip) );
   }

   subscipdata->subscip = NULL;

   /* free the subscip variables */
   if( subscipdata->subvars != NULL )
   {
      assert(subscipdata->nsubvars > 0);
      SCIPfreeBlockMemoryArray(scip, &subscipdata->subvars, subscipdata->nsubvars);
   }

   subscipdataReset(subscipdata);

   return SCIP_OKAY;
}

/** store the sub-SCIP to the data structure */
static
SCIP_RETCODE subscipdataCopySubscip(
   SCIP*                 scip,               /**< original SCIP data structure */
   SUBSCIPDATA*          subscipdata,        /**< data structure of the sub-problem */
   SCIP*                 subscip,            /**< sub scip data structure to keep */
   SCIP_VAR**            subvars,            /**< sub scip variable array in the order of the main SCIP variables */
   int                   nvars               /**< number of sub SCIP variables */
   )
{
   assert(scip != NULL);
   assert(subscipdata != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(nvars == SCIPgetNVars(scip));

   assert(subscipdata->subscip == NULL);
   assert(subscipdata->subvars == NULL);

   subscipdata->subscip = subscip;
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &subscipdata->subvars, subvars, nvars) );
   subscipdata->nsubvars = nvars;

   subscipdata->objbound = SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip));

   return SCIP_OKAY;
}

#ifdef SCIP_DEBUG
/** print debug message listing solving time, nodes, and status of sub-SCIP */
static
SCIP_RETCODE subscipGetInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip             /**< sub SCIP data */
   )
{
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Longint nodelimit;
   SCIP_Real time;
   SCIP_Longint nodes;
   SCIP_STATUS status;

   SCIP_CALL( SCIPgetRealParam(subscip, "limits/time", &timelimit) );
   SCIP_CALL( SCIPgetRealParam(subscip, "limits/memory", &memorylimit) );
   SCIP_CALL( SCIPgetLongintParam(subscip, "limits/nodes", &nodelimit) );

   time = SCIPgetSolvingTime(subscip);
   nodes = SCIPgetNNodes(subscip);
   status = SCIPgetStatus(subscip);

   SCIPdebugMsg(scip, "SCIP info: Time: %.1f (Limit: %.1f) Nodes: %"SCIP_LONGINT_FORMAT" (Limit: %"SCIP_LONGINT_FORMAT") Status: %d\n",
         time, timelimit, nodes, nodelimit, status);

   return SCIP_OKAY;
}
#endif

/** create the objective function based on the user selection */
static
SCIP_RETCODE changeSubvariableObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_VAR*             var,                /**< SCIP variable */
   SCIP_VAR*             subvar,             /**< sub-SCIP variable whose objective coefficient is changed */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure to control how the objective is changed */
   )
{
   SCIP_Real objcoeff;
   SCIP_Real upfrac;
   SCIP_Real downfrac;
   SCIP_Real lpsolval;
   SCIP_Real rootlpsolval;

   /* create the objective value based on the choice of the sub-SCIP objective */
   switch( heurdata->subscipobjective )
   {
      /* use zero as objective function */
      case 'z':
         objcoeff = 0.0;
         break;

      /* use current LP fractionality as objective */
      case 'f':
         lpsolval = SCIPvarGetLPSol(var);
         downfrac = SCIPfrac(scip, lpsolval);
         upfrac = 1.0 - downfrac;

         objcoeff = downfrac - upfrac;
         break;

      /* use root LP solution difference */
      case 'r':
         lpsolval = SCIPvarGetLPSol(var);
         rootlpsolval = SCIPvarGetRootSol(var);
         objcoeff = rootlpsolval - lpsolval;
         break;

      /* use average inferences */
      case 'i':
         objcoeff = SCIPgetVarAvgInferences(scip, var, SCIP_BRANCHDIR_DOWNWARDS)
            - SCIPgetVarAvgInferences(scip, var, SCIP_BRANCHDIR_UPWARDS);
         break;

      /* use original objective function */
      case 'o':
         objcoeff = SCIPvarGetObj(var);
         break;
      default:
         objcoeff = 0.0;
         break;
   }

   SCIP_CALL( SCIPchgVarObj(subscip, subvar, objcoeff) );

   return SCIP_OKAY;
}

/* ---------------- Callback methods of event handler ---------------- */

/** execution callback of the event handler for Lpface sub-SCIP
 *
 * we interrupt the solution process if we hit the LP iteration limit per node
 */
static
SCIP_DECL_EVENTEXEC(eventExecLpface)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetNLPs(scip) > heurdata->lplimfac * heurdata->nodelimit )
   {
      SCIPdebugMsg(scip, "interrupt after  %" SCIP_LONGINT_FORMAT " LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/** setup and solve the subproblem and catch the return code */
static
SCIP_RETCODE setupSubscipLpface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristics data */
   SCIP_VAR**            subvars,            /**< subproblem's variables */
   SCIP_VAR**            vars,               /**< original problem's variables */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   SCIP_Bool*            keepthisscip,       /**< should the subscip be kept or deleted? */
   int                   nvars               /**< number of original problem's variables */
   )
{
   SCIP_HASHMAP* varmapfw = NULL;            /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_Bool success;
   int i;

   assert( subscip != NULL );
   assert( heurdata!= NULL );
   assert( vars != NULL );

   /* create the variable hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
   success = FALSE;

   if( heurdata->uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_lpfacesub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_lpfacesub", SCIPgetProbName(scip));

      /* create the subproblem */
      SCIP_CALL( SCIPcreateProbBasic(subscip, probname) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, NULL, NULL, 0, TRUE) );

      /* copy parameter settings */
      SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );
   }
   else
   {
      SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "lpface", TRUE, FALSE, TRUE, &success) );

      if( heurdata->copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE, NULL) );
      }
   }

   /* fill subvars array with mapping from original variables and set the objective coefficient to the desired value */
   for( i = 0; i < nvars; i++ )
   {
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

      SCIP_CALL( changeSubvariableObjective(scip, subscip, vars[i], subvars[i], heurdata) );

   }
   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   success = FALSE;

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* fix variables that are at their bounds and have nonzero reduced costs  */
   SCIP_CALL( setupSubproblem(scip, subscip, subvars, heurdata, &success) );

   /* if creation of sub-SCIP was aborted (e.g. due to number of fixings), free sub-SCIP and abort */
   if( ! success )
   {
      *result = SCIP_DIDNOTRUN;

      /* this run will be counted as a failure since no new solution tuple could be generated or the neighborhood of the
       * solution was not fruitful in the sense that it was too big
       */
      updateFailureStatistic(scip, heurdata);

      /* we do not want to keep this SCIP */
      *keepthisscip = FALSE;

      return SCIP_OKAY;
   }
   /* set up sub-SCIP parameters */
   SCIP_CALL( setSubscipParameters(scip, subscip) );

   return SCIP_OKAY;
}

/** setup and solve the subproblem and catch the return code */
static
SCIP_RETCODE solveSubscipLpface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_HEUR*            heur,               /**< mutation heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristics data */
   SCIP_VAR**            subvars,            /**< subproblem's variables */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   SCIP_Real             focusnodelb,        /**< lower bound of the focus node */
   SCIP_Bool*            keepthisscip        /**< should the subscip be kept or deleted? */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_Bool success;
   int i;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( heur != NULL );
   assert( heurdata != NULL );
   assert( subvars != NULL );

   /* create event handler for LP events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecLpface, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* determine node, memory, and time limits for the sub-SCIP. Both node and time limit change with every call to
    * the heuristic
    */
   SCIP_CALL( setSubscipLimits(scip, subscip, heur, heurdata, &success) );

   /* if we did not succeed to set the limits of the subscip to let it run, we won't keep it any longer */
   if( !success )
   {
      *keepthisscip = FALSE;

      return SCIP_OKAY;
   }

   /* catch LP events of sub-SCIP */
   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

#ifdef WRITELPFACEPROB
   {
      char probfilename[] = "./lpface_prob.mps";
      char paramfilename[] = "./lpface_prob.set";
      SCIPinfoMessage(scip, NULL, "Writing problem and parameters to file: <%s> <%s>\n", probfilename, paramfilename);
      SCIP_CALL( SCIPwriteOrigProblem(subscip, probfilename, NULL, FALSE) );
      SCIP_CALL( SCIPwriteParams(subscip, paramfilename, TRUE, TRUE) );
   }
#endif

   /* we must not be infeasible at this stage */
   assert(SCIPgetStatus(subscip) != SCIP_STATUS_INFEASIBLE);

   /* solve the subproblem */
   SCIPdebugMsg(scip, "Solve Lpface subMIP\n");
   SCIPdebug(
      SCIP_CALL( subscipGetInfo(scip, subscip) );
   )

   /* Errors in solving the subproblem should not kill the overall solving process.
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   /* save useful information regarding the subscip runs */
   heurdata->usednodes += SCIPgetNNodes(subscip);
   heurdata->submipnlpiters += SCIPgetNLPIterations(subscip);
   heurdata->submippresoltime += SCIPgetPresolvingTime(subscip);
   heurdata->submipstatus = SCIPgetStatus(subscip);

   /* store the focus node lower bound as infeasible to avoid running on this face again */
   if( heurdata->submipstatus == SCIP_STATUS_INFEASIBLE )
   {
      heurdata->lastlpobjinfeas = focusnodelb;
      *keepthisscip = FALSE;
   }
   else if( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;
      int solindex;

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one is accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      solindex = -1;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &solindex, &success) );
      }

      /* we found an optimal solution and are done. Thus, we free the subscip immediately */
      if( success )
      {
         *keepthisscip = FALSE;
         *result = SCIP_FOUNDSOL;
      }

      /* if solution could not be added to problem => run is counted as a failure */
      if( ! success || solindex != SCIPsolGetIndex(SCIPgetBestSol(scip)) )
         updateFailureStatistic(scip, heurdata);
   }
   else
   {
      /* if no new solution was found, run was a failure */
      updateFailureStatistic(scip, heurdata);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLpface)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurLpface(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLpface)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

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
SCIP_DECL_HEURINIT(heurInitLpface)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->nfailures = 0;
   heurdata->nextnodenumber = 0;

   heurdata->submipstatus = SCIP_STATUS_UNKNOWN;
   heurdata->submipnlpiters = -1;
   heurdata->submippresoltime = 0.0;
   heurdata->nvarsfixed = -1;

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolLpface)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* reset the last infeasible objective because it lives in transformed space and must be invalidated at every restart */
   heurdata->lastlpobjinfeas = -SCIPinfinity(scip);

   assert(heurdata->subscipdata == NULL);

   /* variable order might have changed since the last run, reinitialize sub-SCIP data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata->subscipdata) );
   subscipdataReset(heurdata->subscipdata);

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process is exiting) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolLpface)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* variable order might change after restart, free the heuristic subscipdata */
   assert(heurdata->keepsubscip || heurdata->subscipdata->subscip == NULL);
   if( heurdata->subscipdata->subscip != NULL )
   {
      /* free kept data structures first */
      SCIP_CALL( subscipdataFreeSubscip(scip, heurdata->subscipdata) );

   }

   /* free the sub-SCIP data structure */
   SCIPfreeBlockMemory(scip, &heurdata->subscipdata);

   return SCIP_OKAY;
}

#ifdef SCIP_STATISTIC
/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitLpface)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
      "LP Face heuristic stats: Status: %d Nodes: %d LP iters: %d Fixed: %d Presolving time: %.2f\n",
      heurdata->submipstatus, heurdata->usednodes, heurdata->submipnlpiters, heurdata->nvarsfixed, heurdata->submippresoltime);

   return SCIP_OKAY;
}
#else
#define heurExitLpface NULL
#endif

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLpface)
{  /*lint --e{715}*/
   SCIP* subscip;                            /* the subproblem created by lpface       */
   SCIP_HEURDATA* heurdata;                  /* primal heuristic data                  */
   SCIP_VAR** vars;                          /* original problem's variables           */
   SCIP_VAR** subvars;                       /* subproblem's variables                 */
   SCIP_RETCODE retcode;
   SCIP_Bool keepthisscip;
   SCIP_Real focusnodelb;
   SCIP_Real rootlb;
   SCIP_Longint nodelimit;                   /* node limit for the subproblem          */
   int nvars;                                /* number of original problem's variables */
   int nbinvars;
   int nintvars;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DELAYED;

   /* we skip infeasible nodes */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* the node number to run the heuristic again was not yet reached */
   if( SCIPgetNNodes(scip) < heurdata->nextnodenumber )
      return SCIP_OKAY;

   /* do not run heuristic on nodes that were not solved to optimality */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* LP face requires that the LP defines a valid lower bound for the current node */
   if( ! SCIPisLPRelax(scip) || ! SCIPallColsInLP(scip) )
      return SCIP_OKAY;

   assert(SCIPgetCurrentNode(scip) != NULL);
   focusnodelb = SCIPgetNodeLowerbound(scip, SCIPgetCurrentNode(scip));

   /* from the checked conditions, the LP objective should be a valid lower bound for the current node */
   assert(SCIPisGE(scip, focusnodelb, SCIPgetLPObjval(scip)));

   /* do not run if the current focus node already has a lower bound higher than the LP value at the node,
    * for example, due to strong branching
    */
   if( SCIPisGT(scip, focusnodelb, SCIPgetLPObjval(scip)) )
      return SCIP_OKAY;

   /* delay heuristic if the active search tree path is not deep enough */
   if( SCIPgetDepth(scip) < heurdata->minpathlen - 1 )
      return SCIP_OKAY;

   /* only run at lower bound defining nodes */
   if( SCIPisGT(scip, focusnodelb, SCIPgetLowerbound(scip)) )
      return SCIP_OKAY;

   /* only run if lower bound has increased since last LP objective where the sub-MIP was solved to infeasibility */
   if( SCIPisEQ(scip, heurdata->lastlpobjinfeas, focusnodelb) )
      return SCIP_OKAY;

   /* make the reasoning stronger if the objective value must be integral */
   if( SCIPisObjIntegral(scip)
         && (! SCIPisIntegral(scip, focusnodelb) || SCIPisLT(scip, focusnodelb, heurdata->lastlpobjinfeas + 1.0)) )
      return SCIP_OKAY;

   rootlb = SCIPgetLowerboundRoot(scip);
   assert(SCIPisLE(scip, rootlb, focusnodelb));

   /* if the lower bound hasn't changed since the root node, we want to run anyway, otherwise we base our decision on the
    * total path length of the active search tree along which the lower bound did not change anymore.
    */
   if( SCIPisLT(scip, rootlb, focusnodelb) )
   {
      SCIP_NODE* parent;
      int nonimprovingpathlen = 0; /* the length of the current path (in edges) along which the lower bound stayed the same */

      parent = SCIPnodeGetParent(SCIPgetCurrentNode(scip));

      /* count the path length along which the dual bound has not changed */
      while( SCIPisEQ(scip, SCIPnodeGetLowerbound(parent), focusnodelb) && nonimprovingpathlen < heurdata->minpathlen )
      {
         ++nonimprovingpathlen;

         /* we cannot hit the root node because the root lower bound is strictly smaller */
         assert(SCIPnodeGetParent(parent) != NULL);
         parent = SCIPnodeGetParent(parent);
      }

      /* we return if the nonimproving path is too short measured by the heuristic frequency */
      if( nonimprovingpathlen < heurdata->minpathlen )
      {
         /* we do not delay the heuristic if the path has length zero, otherwise it may be called at children so that
          * the path length is sufficient
          */
         if( nonimprovingpathlen == 0 )
            *result = SCIP_DIDNOTRUN;

         return SCIP_OKAY;
      }
   }

   *result = SCIP_DIDNOTRUN;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nodelimit = calcNodeLimit(scip, heur, heurdata);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nodelimit < heurdata->minnodes )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
     return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   assert(nvars > 0);

   /* check whether discrete variables are present */
   if( nbinvars == 0 && nintvars == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   keepthisscip = heurdata->keepsubscip;

   /* check if variable number increased since last call to the sub-SCIP */
   if( heurdata->subscipdata->subscip != NULL && heurdata->subscipdata->nsubvars != nvars )
   {
      SCIPdebugMsg(scip, "Free subscip of LP face heuristic because variable number %d changed since last call (was %d)\n",
         nvars, heurdata->subscipdata->nsubvars);

      SCIP_CALL( subscipdataFreeSubscip(scip, heurdata->subscipdata) );
   }
   else if( heurdata->subscipdata->subscip != NULL && SCIPisGT(scip, focusnodelb, heurdata->subscipdata->objbound) )
   {
      SCIPdebugMsg(scip, "Free subscip of LP face heuristic because of different dual bound: %16.9g > %16.9g\n",
         SCIPretransformObj(scip, focusnodelb), SCIPretransformObj(scip, heurdata->subscipdata->objbound));

      SCIP_CALL( subscipdataFreeSubscip(scip, heurdata->subscipdata) );
   }


   /* retrieve the sub-SCIP from the heuristic data structure */
   if( heurdata->subscipdata->subscip != NULL )
   {
      subscip = heurdata->subscipdata->subscip;
      subvars = heurdata->subscipdata->subvars;
      nvars = heurdata->subscipdata->nsubvars;

      SCIPdebug(
         SCIPdebugMsg(scip, "Loaded sub-SCIP from previous run:\n");
         SCIP_CALL( subscipGetInfo(scip, subscip) );
         )
   }
   else
   {
      assert(heurdata->subscipdata->subscip == NULL);
      SCIPdebugMsg(scip, "Creating new sub-Problem for LP face heuristic\n");

      /* allocate memory to hold sub-SCIP variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

      /* initialize the subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      retcode = setupSubscipLpface(scip, subscip, heurdata, subvars, vars, result, &keepthisscip, nvars);

      SCIP_CALL( retcode );

      if( *result == SCIP_DIDNOTRUN )
         goto TERMINATE;
   }

   retcode = solveSubscipLpface(scip, subscip, heur, heurdata, subvars, result, focusnodelb, &keepthisscip);

   SCIP_CALL( retcode );

TERMINATE:
   /* free subproblem or store it for the next run of the heuristic */
   if( ! keepthisscip )
   {
      /* we only allocated buffer memory if no previous subscip was reinstalled */
      if( heurdata->subscipdata->subscip == NULL )
      {
         SCIPfreeBufferArray(scip, &subvars);
         SCIP_CALL( SCIPfree(&subscip) );
      }
      else
         SCIP_CALL( subscipdataFreeSubscip(scip, heurdata->subscipdata) );

      subscipdataReset(heurdata->subscipdata);
   }
   else
   {
      /* if the subscip has not yet been stored, we copy the subscip into the heuristic data to keep it for the next run */
      if( heurdata->subscipdata->subscip == NULL )
      {
         SCIP_CALL( subscipdataCopySubscip(scip, heurdata->subscipdata, subscip, subvars, nvars) );
         SCIPfreeBufferArray(scip, &subvars);
      }
      else
      {
         assert(heurdata->subscipdata->subscip == subscip);
         assert(heurdata->subscipdata->subvars == subvars);
         assert(heurdata->subscipdata->nsubvars == nvars);
      }
   }

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the lpface primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLpface(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Lpface primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   heurdata->subscipdata = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLpface, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLpface) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLpface) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLpface) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolLpface) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolLpface) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitLpface) );

   /* add lpface primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "required percentage of fixed integer variables in sub-MIP to run",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/dualbasisequations",
         "should dually nonbasic rows be turned into equations?",
         &heurdata->dualbasisequations, TRUE, DEFAULT_DUALBASISEQUATIONS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/keepsubscip",
         "should the heuristic continue solving the same sub-SCIP?",
         &heurdata->keepsubscip, TRUE, DEFAULT_KEEPSUBSCIP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/subscipobjective",
         "objective function in the sub-SCIP: (z)ero, (r)oot-LP-difference, (i)nference, LP (f)ractionality, (o)riginal",
         &heurdata->subscipobjective, TRUE, 'z', "forzi", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minpathlen",
         "the minimum active search tree path length along which lower bound hasn't changed before heuristic becomes active",
         &heurdata->minpathlen, TRUE, DEFAULT_MINPATHLEN, 0, 65531, NULL, NULL) );

   return SCIP_OKAY;
}
