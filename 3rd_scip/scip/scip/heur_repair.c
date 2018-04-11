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

/**@file   heur_repair.c
 * @brief  repair primal heuristic
 * @author Gregor Hendel
 * @author Thomas Nagel
 *
 */

/* This heuristic takes an infeasible solution and tries to repair it.
 * This can happen by variable fixing as long as the sum of all potential possible shiftings
 * is higher than alpha*slack or slack variables with a strong penalty on the objective function.
 * This heuristic cannot run if variable fixing and slack variables are turned off.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_repair.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_varbound.h"

#define HEUR_NAME             "repair"
#define HEUR_DESC             "tries to repair a primal infeasible solution"
#define HEUR_DISPCHAR         '!'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */
#define DEFAULT_MINFIXINGRATE 0.3       /* minimum percentage of integer variables that have to be fixed */

#define DEFAULT_NODESOFS      500       /* number of nodes added to the contingent of the total nodes */
#define DEFAULT_MAXNODES      5000      /* maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINNODES      50        /* minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem */

#define DEFAULT_FILENAME      "-"       /**< file name of a solution to be used as infeasible starting point */
#define DEFAULT_ROUNDIT       TRUE      /**< if it is TRUE : fractional variables which are not fractional in the given
                                         *   solution are rounded, if it is FALSE : solving process of this heuristic
                                         *   is stopped
                                         */
#define DEFAULT_USEOBJFACTOR  FALSE     /**< should a scaled objective function for original variables be used in repair
                                         *   subproblem?
                                         */
#define DEFAULT_USEVARFIX     TRUE      /**< should variable fixings be used in repair subproblem? */
#define DEFAULT_USESLACKVARS  FALSE     /**< should slack variables be used in repair subproblem? */
#define DEFAULT_ALPHA         2.0       /**< how many times the potential should be bigger than the slack? */

/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             infsol;             /**< infeasible solution to start with */
   char*                 filename;           /**< file name of a solution to be used as infeasible starting point */
   SCIP_Longint          usednodes;          /**< number of already used nodes by repair */
   SCIP_Longint          subnodes;           /**< number of nodes which were necessary to solve the sub-SCIP */
   SCIP_Longint          subiters;           /**< contains total number of iterations used in primal and dual simplex
                                              *   and barrier algorithm to solve the sub-SCIP
                                              */
   SCIP_Real             relvarfixed;        /**< relative number of fixed variables */
   SCIP_Real             alpha;              /**< how many times the potential should be bigger than the slack? */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed */
#ifdef SCIP_STATISTIC
   SCIP_Real             relviolatedvars;    /**< relative number of violated variables */
   SCIP_Real             subpresoltime;      /**< time for presolving the sub-SCIP */
   SCIP_Real             relviolatedcons;    /**< relative number of violated cons */
   SCIP_Real             originalsolval;     /**< value of the solution find by repair, in the original Problem*/
   SCIP_Real             improvedoldsol;     /**< value of the given solution after being improved by SCIP */
   int                   nviolatedvars;      /**< number of violated variables in the given solution */
   int                   norigvars;          /**< number of all variables in the given problem */
   int                   nviolatedcons;      /**< number of violated cons in the given solution */
   int                   norcons;            /**< number of all cons in the given problem */
#endif
   int                   nvarfixed;          /**< number of all variables fixed in the sub problem */
   int                   runs;               /**< number of branch and bound runs performed to solve the sub-SCIP */
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   int                   minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Bool             roundit;            /**< if it is TRUE : fractional variables which are not fractional in the
                                              *   given solution are rounded, if it is FALSE : solving process of this
                                              *   heuristic is stopped.
                                              */
   SCIP_Bool             useobjfactor;       /**< should a scaled objective function for original variables be used in
                                              *   repair subproblem?
                                              */
   SCIP_Bool             usevarfix;          /**< should variable fixings be used in repair subproblem? */
   SCIP_Bool             useslackvars;       /**< should slack variables be used in repair subproblem? */
};


/*
 * Local methods
 */

/** computes a factor, so that (factor) * (original objective upper bound) <= 1.*/
static
SCIP_RETCODE getObjectiveFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure */
   SCIP_Real*            factor,             /**< SCIP_Real to save the factor for the old objective function*/
   SCIP_Bool*            success             /**< SCIP_Bool: Is the factor real?*/
   )
{
   SCIP_VAR** vars;
   SCIP_Real lprelaxobj;
   SCIP_Real upperbound;
   SCIP_Real objoffset;
   int nvars;
   int i;

   *success = TRUE;
   *factor = 0.0;
   upperbound = 0.0;

   lprelaxobj = SCIPgetLowerbound(scip);

   if( SCIPisInfinity(scip, -lprelaxobj) )
   {
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) )
   {
      upperbound = SCIPgetUpperbound(scip);
   }
   else
   {
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

      /* tries to find an upper bound for the original objective function, by compute the worst objective value of the
       * LP-relaxation, which holds all variable bounds
       */
      for (i = 0; i < nvars; ++i)
      {
         upperbound = SCIPvarGetObj(vars[i]);
         if( SCIPisInfinity(scip, upperbound) || SCIPisInfinity(scip, -upperbound) )
         {
            /* TODO fancy diving function to find a solution for the max problem */
            *factor = 1 / SCIPinfinity(scip);
            return SCIP_OKAY;
         }
         else if( SCIPisZero(scip, upperbound) )
         {
            continue;
         }
         else if( SCIPisGT(scip, 0.0, upperbound) )
         {
            *factor += upperbound * SCIPvarGetLbGlobal(vars[i]);
         }
         else
         {
            *factor += upperbound * SCIPvarGetUbGlobal(vars[i]);
         }
      }
   }

   /* Ending-sequence */
   *factor = upperbound - lprelaxobj;
   if( !SCIPisZero(scip, *factor) )
   {
      *factor = 1.0 / *factor;
   }

   /* set an offset which guarantees positive objective values */
   objoffset = -lprelaxobj * (*factor);
   SCIP_CALL( SCIPaddOrigObjoffset(subscip, -objoffset) );

   return SCIP_OKAY;
}

/** returns the contributed potential for a variable */
static
SCIP_Real getPotentialContributed(
   SCIP*             scip,                   /**< SCIP data structure */
   SCIP_SOL*         sol,                    /**< infeasible solution */
   SCIP_VAR*         var,                    /**< variable, which potential should be returned */
   SCIP_Real         coefficient,            /**< variables coefficient in corresponding row */
   int               sgn                     /**< sign of the slack */
   )
{
   SCIP_Real potential;

   assert(NULL != scip);
   assert(NULL != var);

   if( 0 > sgn * coefficient )
   {
      if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
      {
         potential = SCIPinfinity(scip);
      }
      else
      {
         potential = coefficient * (SCIPgetSolVal(scip, sol, var) - SCIPvarGetLbGlobal(var));
      }
   }
   else
   {
      if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
      {
         potential = -SCIPinfinity(scip);
      }
      else
      {
         potential = coefficient * (SCIPgetSolVal(scip, sol, var) - SCIPvarGetUbGlobal(var));
      }
   }

   if( SCIPisZero(scip, potential) )
   {
      potential = 0.0;
   }
   return potential;
}

/** finds out if a variable can be fixed with respect to the potentials of all rows, if it is possible, the potentials
 *  of rows are adapted and TRUE is returned.
 */
static
SCIP_RETCODE tryFixVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_SOL*             sol,                /**< solution data structure */
   SCIP_Real*            potential,          /**< array with all potential values */
   SCIP_Real*            slack,              /**< array with all slack values */
   SCIP_VAR*             var,                /**< variable to be fixed? */
   SCIP_VAR*             subvar,             /**< representative variable for var in the sub-SCIP */
   int*                  inftycounter,       /**< counters how many variables have an infinity potential in a row */
   SCIP_HEURDATA*        heurdata,           /**< repairs heuristic data */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   )
{
   SCIP_ROW** rows;
   SCIP_COL* col;
   SCIP_Real* vals;
   SCIP_Real alpha;
   SCIP_Real solval;
   int nrows;
   int i;
   int sgn;
   int rowindex;

   assert(NULL != scip);
   assert(NULL != potential);
   assert(NULL != slack);
   assert(NULL != var);
   assert(NULL != inftycounter);
   assert(NULL != heurdata);

   alpha = heurdata->alpha;
   *infeasible = TRUE;
   *fixed = FALSE;

   solval = SCIPgetSolVal(scip, sol, var);

   if( SCIPisFeasLT(scip, solval, SCIPvarGetLbGlobal(var)) )
   {
      return SCIP_OKAY;
   }
   if( SCIPisFeasGT(scip, solval, SCIPvarGetUbGlobal(var)) )
   {
      return SCIP_OKAY;
   }

   col = SCIPvarGetCol(var);
   rows = SCIPcolGetRows(col);
   nrows = SCIPcolGetNLPNonz(col);
   vals = SCIPcolGetVals(col);

   if( NULL == rows )
   {
      SCIP_CALL( SCIPfixVar(subscip, subvar, solval,
               infeasible, fixed) );
      assert(!*infeasible && *fixed);
      heurdata->nvarfixed++;
      SCIPdebugMsg(scip,"Variable %s is fixed to %g\n",SCIPvarGetName(var), solval);
      return SCIP_OKAY;
   }
   assert(NULL != rows);

   /* iterate over rows, where the variable coefficient is nonzero */
   for( i = 0; i < nrows; ++i )
   {
      SCIP_Real contribution;
      rowindex = SCIProwGetLPPos(rows[i]);
      assert(rowindex >= 0);

      sgn = 1;

      if( SCIPisFeasZero(scip, slack[rowindex]) )
      {
         continue;
      }
      else if( SCIPisFeasGT(scip, 0.0 , slack[rowindex]) )
      {
         sgn = -1;
      }

      contribution = getPotentialContributed(scip, sol, var, vals[i], sgn);

      if( !SCIPisInfinity(scip, REALABS(contribution)) )
      {
         potential[rowindex] -= contribution;
      }
      else
      {
         inftycounter[rowindex]--;
      }

      assert(0 <= inftycounter[rowindex]);
      if( 0 == inftycounter[rowindex] && REALABS(potential[rowindex]) < alpha * REALABS(slack[rowindex]) )
      {
         /* revert the changes before */
         int j = i;
         for( ; j >= 0; --j )
         {
            sgn = 1;
            if( 0 == slack[rowindex] )
            {
               continue;
            }
            rowindex = SCIProwGetLPPos(rows[j]);
            if( 0 > slack[rowindex])
            {
               sgn = -1;
            }
            contribution = getPotentialContributed(scip, sol, var, vals[j], sgn);
            if( !SCIPisInfinity(scip, REALABS(contribution)) )
            {
               potential[rowindex] += contribution;
            }
            else
            {
               inftycounter[rowindex]++;
            }
         }
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPfixVar(subscip, subvar, solval, infeasible, fixed) );
   assert(!*infeasible && *fixed);
   heurdata->nvarfixed++;
   SCIPdebugMsg(scip,"Variable %s is fixed to %g\n",SCIPvarGetName(var),
         SCIPgetSolVal(scip, sol, var));

   return SCIP_OKAY;
}

/** checks if all integral variables in the given solution are integral. */
static
SCIP_RETCODE checkCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution pointer to the to be checked solution */
   SCIP_Bool             roundit,            /**< round fractional solution values of integer variables */
   SCIP_Bool*            success             /**< pointer to store if all integral variables are integral or could
                                              *   be rounded
                                              */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int nfracvars;
   int nbinvars;
   int nintvars;
   int i;

   assert(NULL != success);
   assert(NULL != sol);

   *success = TRUE;

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* check if the candidates are fractional and round them if necessary */
   nfracvars = nbinvars + nintvars;
   for( i = 0; i < nfracvars; ++i)
   {
      SCIP_Real value = SCIPgetSolVal(scip, sol, vars[i]);

      if( SCIPisInfinity(scip, REALABS(value)) )
      {
         *success = FALSE;
         SCIPdebugMsg(scip, "Variable with infinite solution value");

         return SCIP_OKAY;
      }
      if( !SCIPisFeasIntegral(scip, value) )
      {
         if( roundit )
         {
            SCIP_Real roundedvalue;

            if( SCIPvarGetNLocksUp(vars[i]) > SCIPvarGetNLocksDown(vars[i]) )
            {
               roundedvalue = SCIPceil(scip, value - 1.0);
            }
            else
            {
               roundedvalue = SCIPfloor(scip, value + 1.0);
            }

            SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], roundedvalue) );
         }
         else
         {
            *success = FALSE;
            SCIPdebugMsg(scip, "Repair: All variables are integral.\n");
            return SCIP_OKAY;
         }
      }
   }

   /* ensure that no other variables have infinite LP solution values */
   for( ; i < nvars; ++i )
   {
      if( SCIPisInfinity(scip, REALABS(SCIPgetSolVal(scip, sol, vars[i]))) )
      {
         *success = FALSE;
         SCIPdebugMsg(scip, "Variable with infinite solution value");

         return SCIP_OKAY;
      }
   }

   SCIPdebugMsg(scip, "All variables rounded.\n");
   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< Repair heuristic structure                          */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;                         /* the original problem's number of variables      */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

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
   /* try to add new solution to SCIP and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

#ifdef SCIP_STATISTICS
   {
      SCIP_HEURDATA* heurdata;
      heurdata = SCIPheurGetData(heur);

      if( *success )
      {
         heurdata->originalsolval = SCIPgetSolOrigObj(scip, newsol);
      }
   }
#endif

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** tries to fix variables as an approach to repair a solution. */
static
SCIP_RETCODE applyRepair(
   SCIP*                 scip,               /**< SCIP data structure of the problem */
   SCIP_HEUR*            heur,               /**< pointer to this heuristic instance */
   SCIP_RESULT*          result,             /**< pointer to return the result status */
   SCIP_Longint          nnodes              /**< nodelimit for sub-SCIP */
   )
{
   SCIP* subscip = NULL;
   SCIP_VAR** vars = NULL;
   SCIP_VAR** subvars = NULL;
   SCIP_ROW** rows;
   SCIP_CONS** subcons = NULL;
   int* nviolatedrows = NULL;
   int* permutation = NULL;
   int* inftycounter = NULL;
   SCIP_SOL* sol;
   SCIP_SOL* subsol = NULL;
   SCIP_HEURDATA* heurdata;
   SCIP_Real* potential = NULL;
   SCIP_Real* slacks = NULL;
   SCIP_RETCODE retcode = SCIP_OKAY;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Real factor;
   char probname[SCIP_MAXSTRLEN];
   int i;
   int nbinvars;
   int nintvars;
   int nvars;
   int nrows;
   int ndiscvars;
   int nfixeddiscvars;
   SCIP_Bool success;

   heurdata = SCIPheurGetData(heur);
   sol = heurdata->infsol;

   /* initializes the sub-SCIP */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* get name of the original problem and add the string "_repairsub" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_repairsub", SCIPgetProbName(scip));

   SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* a trivial feasible solution can be constructed if violations are modeled with slack variables */
   if( heurdata->useslackvars )
   {
      SCIP_CALL( SCIPcreateSol(subscip, &subsol, heur) );
   }

   /* gets all original variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nviolatedrows, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permutation, nvars) );

   SCIPdebugMsg(scip,"\n\n Calling objective factor calculation \n\n");
   if( heurdata->useobjfactor )
   {
      SCIP_CALL( getObjectiveFactor(scip, subscip, &factor, &success) );
   }
   else
   {
      factor = 0.0;
   }

   /* adds all original variables */
   ndiscvars = 0;
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real lborig;
      SCIP_Real uborig;
      SCIP_Real varslack;
      SCIP_Real objval;
      SCIP_Real value;
      SCIP_VARTYPE vartype;
      char varname[SCIP_MAXSTRLEN];
      char slackvarname[SCIP_MAXSTRLEN];
      char consvarname[SCIP_MAXSTRLEN];

#ifdef SCIP_STATISTIC
      heurdata->norigvars++;
#endif

      varslack = 0.0;
      lborig = SCIPvarGetLbGlobal(vars[i]);
      uborig = SCIPvarGetUbGlobal(vars[i]);
      value = SCIPgetSolVal(scip, sol, vars[i]);
      vartype = SCIPvarGetType(vars[i]);

      nviolatedrows[i] = 0;

      /* if the value of x is lower than the variables lower bound, sets the slack to a correcting value */
      if( heurdata->useslackvars && SCIPisFeasLT(scip, value, lborig) )
      {
         lb = value;
         varslack = lborig - value;
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], lb) );
      }
      else
      {
         lb = lborig;
      }

      /* if the value of x is bigger than the variables upper bound, sets the slack to a correcting value */
      if( heurdata->useslackvars && SCIPisFeasGT(scip, value, uborig) )
      {
         ub = value;
         varslack = uborig - value;
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], ub) );
      }
      else
      {
         ub = uborig;
      }

      if( heurdata->useobjfactor )
      {
         objval = SCIPvarGetObj(vars[i])*factor;

         if( SCIPisZero(scip, objval) )
         {
            objval = 0.0;
         }
      }
      else
      {
         objval = SCIPvarGetObj(vars[i]);
      }

      /* if a binary variable is out of bound, generalize it to an integer variable */
      if( !SCIPisFeasZero(scip, varslack) && SCIP_VARTYPE_BINARY == vartype )
      {
         vartype = SCIP_VARTYPE_INTEGER;
         SCIP_CALL( SCIPchgVarType(subscip, subvars[i], vartype, &success) );
      }

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "sub_%s", SCIPvarGetName(vars[i]));

      /* Adds the sub representing variable to the sub-SCIP. */
      SCIP_CALL( SCIPcreateVarBasic(subscip, &subvars[i], varname, lb, ub, objval, vartype) );
      SCIP_CALL( SCIPaddVar(subscip, subvars[i]) );

      /* a trivial feasible solution can be constructed if violations are modeled with slack variables */
      if( heurdata->useslackvars )
      {
         SCIP_CALL( SCIPsetSolVal(subscip, subsol, subvars[i], value) );
      }

      /* if necessary adds a constraint to represent the original bounds of x.*/
      if( !SCIPisFeasEQ(scip, varslack, 0.0) )
      {
         SCIP_VAR* newvar;
         (void) SCIPsnprintf(slackvarname, SCIP_MAXSTRLEN, "artificialslack_%s", SCIPvarGetName(vars[i]));
         (void) SCIPsnprintf(consvarname, SCIP_MAXSTRLEN, "boundcons_%s", SCIPvarGetName(vars[i]));

         /* initialize and add an artificial slack variable */
         if( heurdata->useobjfactor )
         {
            SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, slackvarname, 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS));
         }
         else
         {
            SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, slackvarname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY));
         }
         SCIP_CALL( SCIPaddVar(subscip, newvar) );

         /* set the value of the slack variable to 1 to punish the use of it.
          * note that a trivial feasible solution can be only constructed if violations are modeled with slack variables
          */
         if( heurdata->useslackvars )
         {
            SCIP_CALL( SCIPsetSolVal(subscip, subsol, newvar, 1.0) );
         }

         /* adds a linear constraint to represent the old bounds */
         SCIP_CALL( SCIPcreateConsBasicVarbound(subscip, &cons, consvarname, subvars[i], newvar, varslack, lb, ub) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseVar(subscip, &newvar) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

         /* increases the counter for violated vars */
#ifdef SCIP_STATISTIC
         heurdata->nviolatedvars++;
#endif
      }


#ifdef SCIP_STATISTIC
      if( SCIPisFeasLT(scip, value, lb) || SCIPisFeasGT(scip, value, ub) )
      {
         heurdata->nviolatedvars++;
      }
#endif
      if( SCIP_VARTYPE_BINARY == vartype || SCIP_VARTYPE_INTEGER == vartype )
      {
         ndiscvars++;
      }
   }

   /* check solution for feasibility regarding the LP rows (SCIPgetRowSolActivity()) */
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &potential, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &slacks, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subcons, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inftycounter, nrows) );

   /* Adds all original constraints and computes potentials and slacks */
   for (i = 0; i < nrows; ++i)
   {
      SCIP_COL** cols;
      SCIP_VAR** consvars;
      SCIP_Real* vals;
      SCIP_Real constant;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real rowsolact;
      int nnonz;
      int j;

#ifdef SCIP_STATISTIC
      heurdata->norcons++;
#endif

      /* gets the values to check the constraint */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIPisInfinity(scip, -SCIProwGetLhs(rows[i])) ? SCIProwGetLhs(rows[i]) : SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIPisInfinity(scip, SCIProwGetRhs(rows[i])) ? SCIProwGetRhs(rows[i]) : SCIProwGetRhs(rows[i]) - constant;
      rowsolact = SCIPgetRowSolActivity(scip, rows[i], sol) - constant;
      vals = SCIProwGetVals(rows[i]);
      potential[i] = 0.0;
      inftycounter[i] = 0;

      assert(SCIPisFeasLE(scip, lhs, rhs));

      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);
      SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, nnonz) );

      /* sets the slack if its necessary */
      if( SCIPisFeasLT(scip, rowsolact, lhs) )
      {
         slacks[i] = lhs - rowsolact;
#ifdef SCIP_STATISTIC
         heurdata->nviolatedcons++;
#endif
      }
      else if( SCIPisFeasGT(scip, rowsolact, rhs) )
      {
         slacks[i] = rhs - rowsolact;
#ifdef SCIP_STATISTIC
         heurdata->nviolatedcons++;
#endif
      }
      else
      {
         slacks[i] = 0.0;
      }

      /* translate all variables from the original SCIP to the sub-SCIP with sub-SCIP variables. */
      for( j = 0; j < nnonz; ++j )
      {
         SCIP_Real contribution;
         int pos;
         int sgn = 1;

         /* negative slack represents a right hand side violation */
         if( SCIPisFeasGT(scip, 0.0, slacks[i]) )
         {
            assert(!SCIPisInfinity(scip, rhs));
            sgn = -1;
         }
      #ifndef NDEBUG
         else
            assert(!SCIPisInfinity(scip, lhs));
      #endif

         pos = SCIPvarGetProbindex(SCIPcolGetVar(cols[j]));
         consvars[j] = subvars[pos];
         assert(pos >= 0);

         /* compute potentials */
         contribution = getPotentialContributed(scip, sol, vars[pos], vals[j], sgn);
         if( !SCIPisInfinity(scip, REALABS(contribution)) )
         {
            potential[i] += contribution;
         }
         else
         {
            inftycounter[i]++;
         }

         if( !SCIPisZero(scip, slacks[i]) )
         {
            nviolatedrows[pos]++;
         }
      }


      /* create a new linear constraint, representing the old one */
      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &subcons[i], SCIProwGetName(rows[i]),
               nnonz, consvars, vals, lhs, rhs) );

      if( heurdata->useslackvars )
      {
         SCIP_VAR* newvar;
         char varname[SCIP_MAXSTRLEN];

         /*if necessary adds a new artificial slack variable*/
         if( !SCIPisFeasEQ(subscip, slacks[i], 0.0) )
         {
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "artificialslack_%s", SCIProwGetName(rows[i]));
            SCIP_CALL( SCIPcreateVarBasic(subscip, &newvar, varname, 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(subscip, newvar) );

            /* a trivial feasible solution can be constructed if violations are modeled with slack variables */
            SCIP_CALL( SCIPsetSolVal(subscip, subsol, newvar, 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, subcons[i], newvar, slacks[i]) );
            SCIP_CALL( SCIPreleaseVar(subscip, &newvar) );
         }
      }

      /*Adds the Constraint and release it.*/
      SCIP_CALL( SCIPaddCons(subscip, subcons[i]) );
      SCIP_CALL( SCIPreleaseCons(subscip, &subcons[i]) );
      SCIPfreeBufferArray(subscip, &consvars);
   }

   if( heurdata->usevarfix )
   {
      /* get the greedy order */
      for( i = 0; i < nvars; ++i )
      {
         permutation[i] = i;
      }
      SCIPsortIntInt(nviolatedrows, permutation, nvars);

      /* loops over variables and greedily fix variables, but preserve the cover property that enough slack is given to
       * violated rows
       */
      nfixeddiscvars = 0;
      heurdata->nvarfixed = 0;
      for( i = 0; i < nvars; ++i )
      {
         SCIP_Bool infeasible = FALSE;
         SCIP_Bool fixed = TRUE;

         SCIP_CALL( tryFixVar(scip, subscip, sol, potential, slacks, vars[permutation[i]], subvars[permutation[i]], inftycounter, heurdata, &infeasible, &fixed) );

         if( fixed && (SCIP_VARTYPE_BINARY == SCIPvarGetType(subvars[permutation[i]])
            || SCIP_VARTYPE_INTEGER == SCIPvarGetType(subvars[permutation[i]])) )
         {
            nfixeddiscvars++;
         }
       }
      SCIPdebugMsg(scip,"fixings finished\n\n");
      if( heurdata->minfixingrate > ((SCIP_Real)nfixeddiscvars/MAX((SCIP_Real)ndiscvars,1.0)) )
      {
         goto TERMINATE;
      }
   }

   /* a trivial feasible solution can be constructed if violations are modeled with slack variables */
   if( heurdata->useslackvars )
   {
      SCIP_CALL( SCIPaddSolFree(subscip, &subsol, &success) );

      if( !success )
      {
         SCIPdebugMsg(scip, "Initial repair solution was not accepted.\n");
      }
   }

#ifdef SCIP_STATISTIC
   if( heurdata->useslackvars )
      heurdata->improvedoldsol = SCIPgetSolOrigObj(subscip, subsol);
#endif

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   /* subtract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip) / 1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip) / 1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= 2.0 * SCIPgetMemExternEstim(scip) / 1048576.0 )
      goto TERMINATE;

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nnodes) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetObjlimit(subscip,1.0) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", (int)SCIP_VERBLEVEL_NONE) );

#ifdef SCIP_DEBUG
   /* for debugging Repair, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", (int)SCIP_VERBLEVEL_FULL) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", -1) );
#endif

   /* solve the subproblem */
   retcode = SCIPsolve(subscip);

   /* errors in sub-SCIPs should not kill the overall solving process. Hence, we print a warning message. Only
    * in debug mode, SCIP will stop
    */
   if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while solving subproblem in REPAIR heuristic; sub-SCIP terminated with code <%d>\n", retcode);
      SCIPABORT();  /*lint --e{527}*/
      goto TERMINATE;
   }

   success = FALSE;

   /* if a solution is found, save its value and create a new solution instance for the original SCIP */
   if( SCIPgetBestSol(subscip) != NULL )
   {
#ifdef SCIP_STATISTIC
      heurdata->improvedoldsol = SCIPgetSolOrigObj(subscip, SCIPgetBestSol(subscip));
#endif
      /* print solving statistics of subproblem if we are in SCIP's debug mode */
      SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

      assert(SCIPgetNSols(subscip) > 0);
      SCIP_CALL( createNewSol(scip, subscip, subvars, heur, SCIPgetBestSol(subscip), &success) );

      if( success )
      {
         *result = SCIP_FOUNDSOL;
      }
   }
   else
   {
      SCIPdebugMsg(scip,"No solution found!\n");
   }

   if( SCIPgetStage(subscip) >= SCIP_STAGE_SOLVED )
   {
      heurdata->subiters = SCIPgetNLPIterations(subscip);
      heurdata->subnodes = SCIPgetNTotalNodes(subscip);
#ifdef SCIP_STATISTIC
      heurdata->subpresoltime = SCIPgetPresolvingTime(subscip);
#endif
      heurdata->runs = SCIPgetNRuns(subscip);
   }

   /* terminates the solving process  */
TERMINATE:
   if( NULL != sol )
   {
      SCIP_CALL( SCIPfreeSol(scip, &sol) );
   }
   SCIPfreeBufferArrayNull(scip, &nviolatedrows);
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &subvars[i]) );
   }
   SCIPfreeBufferArrayNull(scip, &inftycounter);
   SCIPfreeBufferArrayNull(scip, &subcons);
   SCIPfreeBufferArrayNull(scip, &slacks);
   SCIPfreeBufferArrayNull(scip, &potential);
   SCIPfreeBufferArrayNull(scip, &permutation);
   SCIPfreeBufferArrayNull(scip, &subvars);

   if( NULL != subsol )
   {
      SCIP_CALL( SCIPfreeSol(subscip, &subsol) );
   }

   SCIP_CALL( SCIPfree(&subscip) );

   SCIPdebugMsg(scip, "repair finished\n");
   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRepair)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);

   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRepair)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   heurdata->subiters = -1;
   heurdata->subnodes = -1;
   heurdata->runs = 0;

   heurdata->nvarfixed = 0;
   heurdata->relvarfixed = -1;

#ifdef SCIP_STATISTIC
   heurdata->subpresoltime = 0;

   heurdata->nviolatedvars = 0;
   heurdata->norigvars = 0;
   heurdata->relviolatedvars = 0;
   heurdata->nviolatedcons = 0;
   heurdata->norcons = 0;
   heurdata->relviolatedcons = 0;

   heurdata->originalsolval = SCIP_INVALID;

   heurdata->improvedoldsol = SCIP_UNKNOWN;
#endif

   heurdata->usednodes = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRepair)
{  /*lint --e{715}*/
#ifdef SCIP_STATISTIC
   SCIP_HEURDATA* heurdata;
   SCIP_Real time;
   SCIP_Real relvars;
   SCIP_Real relcons;
   SCIP_Real relfixed;
   char solval[SCIP_MAXSTRLEN];
   int violateds;
   int ninvars;
   int ninvcons;
   int nvars;
   int ncons;
   int iterations;
   int nodes;
   int runs;

   heurdata = SCIPheurGetData(heur);
   violateds = heurdata->nviolatedvars+heurdata->nviolatedcons;
   ninvars = heurdata->nviolatedvars;
   ninvcons = heurdata->nviolatedcons;
   nvars = heurdata->norigvars;
   ncons = heurdata->norcons;
   iterations = heurdata->subiters;
   nodes = heurdata->subnodes;
   time = heurdata->subpresoltime;
   runs = heurdata->runs;

   if( SCIP_INVALID == heurdata->originalsolval )
   {
      (void) SCIPsnprintf(solval, SCIP_MAXSTRLEN ,"--");
   }
   else
   {
      (void) SCIPsnprintf(solval, SCIP_MAXSTRLEN, "%15.9g", heurdata->originalsolval);
   }

   heurdata->relviolatedvars = MAX((SCIP_Real)heurdata->norigvars, 1.0);
   heurdata->relviolatedvars = heurdata->nviolatedvars/heurdata->relviolatedvars;
   heurdata->relviolatedcons = MAX((SCIP_Real)heurdata->norcons, 1.0);
   heurdata->relviolatedcons = heurdata->nviolatedcons/heurdata->relviolatedcons;

   heurdata->relvarfixed = MAX((SCIP_Real)heurdata->norigvars, 1.0);
   heurdata->relvarfixed = heurdata->nvarfixed/heurdata->relvarfixed;
   relvars = heurdata->relviolatedvars;
   relcons = heurdata->relviolatedcons;
   relfixed = heurdata->relvarfixed;

   /* prints all statistic data for a user*/
   SCIPstatistic(
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "<repair> \n total violations             : %10d\n", violateds);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " violated variables           : %10d\n", ninvars);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " total variables              : %10d\n", nvars)
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " relative violated variables  : %10.2f%%\n", 100 * relvars);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " violated constraints         : %10d\n", ninvcons);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " total constraints            : %10d\n", ncons);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " relative violated constraints: %10.2f%%\n", 100* relcons);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " fixed variables              : %10d\n", heurdata->nvarfixed);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " relative fixed variables     : %10.2f%%\n\n", 100* relfixed);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " iterations                   : %10d\n", iterations);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " nodes                        : %10d\n", nodes);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " number of runs               : %10d\n", runs);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " presolve time                : %10.2f\n", time);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "</repair>\n\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " value of best solution       : %10g\n", solval);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " improved orig. solval        : %10g\n", heurdata->improvedoldsol);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, message);
   )

#endif
   return SCIP_OKAY;
}

/** execution method of primal heuristic. Repair needs an incorrect solution, in which all variables are in their bound. */
static
SCIP_DECL_HEUREXEC(heurExecRepair)
{ /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_RETCODE retcode;
   SCIP_Bool success;
   SCIP_Bool error;
   SCIP_Longint nnodes;

   heurdata = SCIPheurGetData(heur);
   SCIPdebugMsg(scip, "%s\n", heurdata->filename);

   /* checks the result pointer */
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   /* if repair already ran or neither variable fixing nor slack variables are enabled, stop */
   if( 0 < SCIPheurGetNCalls(heur) || !(heurdata->usevarfix || heurdata->useslackvars) )
      return SCIP_OKAY;

   /* do not run if the neither the LP is constructed nor a user given solution exists */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && strcmp(heurdata->filename, DEFAULT_FILENAME) == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward REPAIR if it succeeded often */
   nnodes = (SCIP_Longint)(nnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nnodes -= (SCIP_Longint)(100.0 * SCIPheurGetNCalls(heur));  /* count the setup costs for the sub-MIP as 100 nodes */
   nnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nnodes -= heurdata->usednodes;
   nnodes = MIN(nnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nnodes < heurdata->minnodes )
      return SCIP_OKAY;

   if( !SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_Bool cutoff;

      SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

      /* manually cut off the node if the LP construction detected infeasibility (heuristics cannot return such a result) */
      if( cutoff )
      {
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetCurrentNode(scip)) );
         return SCIP_OKAY;
      }
   }

   /* create original solution */
   SCIP_CALL( SCIPcreateOrigSol(scip, &(heurdata->infsol), heur) );

   /* use read method to enter solution from a file */
   if( strcmp(heurdata->filename, DEFAULT_FILENAME) == 0 )
   {
      retcode = SCIPlinkLPSol(scip, heurdata->infsol);
   }
   else
   {
      error = FALSE;
      retcode = SCIPreadSolFile(scip, heurdata->filename, heurdata->infsol, FALSE, NULL, &error);
      assert(error || SCIP_OKAY == retcode);
   }

   if( SCIP_NOFILE == retcode )
   {
      assert(strcmp(heurdata->filename, DEFAULT_FILENAME) != 0);
      SCIPwarningMessage(scip, "cannot open file <%s> for reading\n", heurdata->filename);

      SCIP_CALL( SCIPfreeSol(scip, &(heurdata->infsol)) );
      return SCIP_OKAY;
   }
   else if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "cannot run repair, unknown return status <%d>\n", retcode);
      SCIP_CALL( SCIPfreeSol(scip, &(heurdata->infsol)) );
      return SCIP_OKAY;
   }
   SCIPdebugMsg(scip, "Repair: Solution file read.\n");

   /* checks the integrality of all discrete variable */
   SCIP_CALL( checkCands(scip, heurdata->infsol, heurdata->roundit, &success) );
   if( !success )
   {
      SCIPdebugMsg(scip,"Given solution is not integral, repair terminates.\n");
      SCIP_CALL( SCIPfreeSol(scip, &(heurdata->infsol)) );
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPtrySol(scip, heurdata->infsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );

   /* the solution is not feasible for the original problem; we will try to repair it */
   if( !success )
   {
      assert(NULL != heurdata->infsol);
      assert(heurdata->usevarfix || heurdata->useslackvars);
      SCIP_CALL( applyRepair(scip, heur, result, nnodes) );
   }
   else
   {
      SCIP_CALL( SCIPfreeSol(scip, &(heurdata->infsol)) );
   }

   return SCIP_OKAY;
}

/* primal heuristic specific interface methods */

/** creates the repair primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRepair(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create repair primal heuristic data */
   heurdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip ,&heurdata) );

   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRepair, heurdata) );

   assert(heur != NULL);
   assert(heurdata != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRepair) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRepair) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRepair) );

   /* add repair primal heuristic parameters */

   heurdata->filename = NULL;
   /* add string parameter for filename containing a solution */
   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/" HEUR_NAME "/filename",
         "file name of a solution to be used as infeasible starting point, [-] if not available",
         &heurdata->filename, FALSE, DEFAULT_FILENAME, NULL, NULL) );

   /* add bool parameter for decision how to deal with unfractional cands */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/roundit",
         "True : fractional variables which are not fractional in the given solution are rounded, "
         "FALSE : solving process of this heuristic is stopped. ",
         &heurdata->roundit, FALSE, DEFAULT_ROUNDIT, NULL, NULL));

   /* add bool parameter for decision how the objective function should be */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useobjfactor",
         "should a scaled objective function for original variables be used in repair subproblem?",
         &heurdata->useobjfactor, FALSE, DEFAULT_USEOBJFACTOR, NULL, NULL));

   /* add bool parameter for decision if variable fixings should be used */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usevarfix",
         "should variable fixings be used in repair subproblem?",
         &heurdata->usevarfix, FALSE, DEFAULT_USEVARFIX, NULL, NULL));

   /* add bool parameter for decision how the objective function should be */
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useslackvars",
         "should slack variables be used in repair subproblem?",
         &heurdata->useslackvars, FALSE, DEFAULT_USESLACKVARS, NULL, NULL));

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/alpha", "factor for the potential of var fixings",
         &heurdata->alpha, TRUE, DEFAULT_ALPHA, 0.0, 100.00, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );


   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "minimum percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
