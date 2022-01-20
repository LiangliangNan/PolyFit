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

/**@file   heur_completesol.c
 * @brief  COMPLETESOL - primal heuristic trying to complete given partial solutions
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "scip/heur_completesol.h"
#include "scip/scipdefplugins.h"       /* needed for the secondary SCIP instance */
#include "scip/pub_misc.h"
#include "scip/def.h"

#define HEUR_NAME             "completesol"
#define HEUR_DESC             "primal heuristic trying to complete given partial solutions"
#define HEUR_DISPCHAR         'h'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

/* default values for heuristic plugins */
#define DEFAULT_MAXNODES      5000LL    /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MAXUNKRATE    0.85      /**< maximum percentage of unknown solution values */
#define DEFAULT_ADDALLSOLS   FALSE      /**< should all subproblem solutions be added to the original SCIP? */
#define DEFAULT_MINNODES      50LL      /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS      500LL     /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT     0.1       /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_LPLIMFAC      2.0       /**< factor by which the limit on the number of LP depends on the node limit */
#define DEFAULT_OBJWEIGHT     1.0       /**< weight of the original objective function (1: only original objective) */
#define DEFAULT_BOUNDWIDENING 0.1       /**< bound widening factor applied to continuous variables
                                         *   (0: round bounds to next integer, 1: relax to global bounds)
                                         */
#define DEFAULT_MINIMPROVE    0.01      /**< factor by which the incumbent should be improved at least */
#define DEFAULT_MINOBJWEIGHT 1e-3       /**< minimal weight for original objective function (zero could lead to infinite solutions) */
#define DEFAULT_IGNORECONT  FALSE       /**< should solution values for continuous variables be ignored? */
#define DEFAULT_BESTSOLS        5       /**< heuristic stops, if the given number of improving solutions were found (-1: no limit) */
#define DEFAULT_MAXPROPROUNDS  10       /**< maximal number of iterations in propagation (-1: no limit) */
#define DEFAULT_MAXLPITER      -1LL     /**< maximal number of LP iterations (-1: no limit) */
#define DEFAULT_MAXCONTVARS    -1       /**< maximal number of continuous variables after presolving (-1: no limit) */
#define DEFAULT_BEFOREPRESOL  TRUE      /**< should the heuristic run before presolving? */

/* event handler properties */
#define EVENTHDLR_NAME         "Completesol"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Longint          maxlpiter;          /**< maximal number of LP iterations (-1: no limit) */
   SCIP_Real             maxunknownrate;     /**< maximal rate of changed coefficients in the objective function */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   SCIP_Real             objweight;          /**< weight of the original objective function (1: only original obj, 0: try to keep to given solution) */
   SCIP_Real             boundwidening;      /**< bound widening factor applied to continuous variables
                                              *   (0: fix variables to given solution values, 1: relax to global bounds)
                                              */
   SCIP_Real             minimprove;         /**< factor by which the incumbent should be improved at least */
   SCIP_Bool             addallsols;         /**< should all subproblem solutions be added to the original SCIP? */
   SCIP_Bool             ignorecont;         /**< should solution values for continuous variables be ignored? */
   SCIP_Bool             beforepresol;       /**< should the heuristic run before presolving? */
   int                   bestsols;           /**< heuristic stops, if the given number of improving solutions were found (-1: no limit) */
   int                   maxcontvars;        /**< maximal number of continuous variables after presolving (-1: no limit) */
   int                   maxproprounds;      /**< maximal number of iterations in propagation (-1: no limit) */
};

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecCompletesol)
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
      SCIPdebugMsg(scip, "interrupt after %" SCIP_LONGINT_FORMAT " LPs\n",SCIPgetNLPs(scip));
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

/** creates a subproblem by fixing a number of variables */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_HEURDATA*        heurdata,           /**< heuristic's private data structure */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_SOL*             partialsol,         /**< partial solution */
   SCIP_Bool*            tightened,          /**< array to store for which variables we have found bound tightenings */
   SCIP_Bool*            success             /**< pointer to store whether the creation was successful */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS* objcons;
   SCIP_Real epsobj;
   SCIP_Real cutoff;
   SCIP_Real upperbound;
   char consobjname[SCIP_MAXSTRLEN];
   int nvars;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(heurdata != NULL);

   *success = TRUE;

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip) + heurdata->minimprove * SCIPgetLowerbound(scip);
      else
      {
         if( SCIPgetUpperbound(scip) >= 0 )
            cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
         else
            cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
      }
      cutoff = MIN(upperbound, cutoff);
      SCIPdebugMsg(scip, "set cutoff=%g for sub-SCIP\n", cutoff);
   }
   else
      cutoff = SCIPinfinity(scip);

   /* calculate objective coefficients for all potential epsilons */
   if( SCIPisEQ(scip, heurdata->objweight, 1.0) )
      return SCIP_OKAY;
   else if( !SCIPisInfinity(scip, cutoff) )
      epsobj = 1.0;
   else
   {
      /* divide by objweight to avoid changing objective coefficient of original problem variables */
      epsobj = (1.0 - heurdata->objweight)/heurdata->objweight;

      /* scale with -1 if we have a maximization problem */
      if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE )
         epsobj *= -1.0;
   }

   /* get active variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   objcons = NULL;

   /* add constraints to measure the distance to the given partial solution */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Real solval;
      int idx;

      assert(SCIPvarIsActive(vars[i]));

      /* add objective function as a constraint, if a primal bound exists */
      if( SCIPisInfinity(scip, cutoff) )
      {
         /* create the constraints */
         if( objcons == NULL )
         {
            SCIP_Real lhs;
            SCIP_Real rhs;

            if( SCIPgetObjsense(subscip) == SCIP_OBJSENSE_MINIMIZE )
            {
               lhs = -SCIPinfinity(subscip);
               rhs = cutoff;
            }
            else
            {
               lhs = cutoff;
               rhs = SCIPinfinity(subscip);
            }

            (void)SCIPsnprintf(consobjname, SCIP_MAXSTRLEN, "obj");
            SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &objcons, consobjname, 0, NULL, NULL, lhs, rhs) );
         }

         /* add the variable to the constraints */
         SCIP_CALL( SCIPaddCoefLinear(subscip, objcons, subvars[i], SCIPvarGetObj(subvars[i])) );

         /* set objective coefficient to 0.0 */
         SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], 0.0) );
      }

      solval = SCIPgetSolVal(scip, partialsol, vars[i]);

      /* skip variables with unknown solution value */
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         continue;

      idx = SCIPvarGetProbindex(vars[i]);
      assert(idx >= 0);

      /* skip variables where we already found some bound tightenings */
      if( tightened[idx] == FALSE )
      {
         /* special case: vars[i] is binary; we do not add an extra variable, but we mimic the behavior we would get with it.
          * E.g., if the solval is 0.3, setting the variable to 0 would give a cost of 0.3 * epsobj, setting it to 1 gives
          * 0.7 * epsobj. Thus, 0.3 * epsobj can be treated as a constant in the objective function and the variable gets
          * an objective coefficient of 0.4 * epsobj.
          */
         if( SCIPvarIsBinary(vars[i]) )
         {
            SCIP_Real frac = SCIPfeasFrac(scip, solval);
            SCIP_Real objcoef;

            frac = MIN(frac, 1-frac);
            objcoef = (1 - 2*frac) * epsobj * (int)SCIPgetObjsense(scip);

            if( solval > 0.5 )
            {
               SCIP_CALL( SCIPchgVarObj(scip, vars[i], -objcoef) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarObj(scip, vars[i], objcoef) );
            }
         }
         else
         {
            SCIP_CONS* conspos;
            SCIP_CONS* consneg;
            SCIP_VAR* eps;
            char consnamepos[SCIP_MAXSTRLEN];
            char consnameneg[SCIP_MAXSTRLEN];
            char epsname[SCIP_MAXSTRLEN];

            /* create two new variables */
            (void)SCIPsnprintf(epsname, SCIP_MAXSTRLEN, "eps_%s", SCIPvarGetName(subvars[i]));

            SCIP_CALL( SCIPcreateVarBasic(subscip, &eps, epsname, 0.0, SCIPinfinity(scip), epsobj, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(subscip, eps) );

            /* create two constraints */
            (void)SCIPsnprintf(consnamepos, SCIP_MAXSTRLEN, "cons_%s_pos", SCIPvarGetName(subvars[i]));
            (void)SCIPsnprintf(consnameneg, SCIP_MAXSTRLEN, "cons_%s_neq", SCIPvarGetName(subvars[i]));

            /* x_{i} - s_{i} <= e_{i}   <==>   x_{i} - e_{i} <= s_{i} */
            SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &conspos, consnamepos, 0, NULL, NULL, -SCIPinfinity(scip), solval) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, conspos, subvars[i], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, conspos, eps, -1.0) );
            SCIP_CALL( SCIPaddCons(subscip, conspos) );
            SCIP_CALL( SCIPreleaseCons(subscip, &conspos) );

            /* s_{i} - x_{i} <= e_{i}   <==>   e_{i} - x_{i} >= s_{i} */
            SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &consneg, consnameneg, 0, NULL, NULL, solval, SCIPinfinity(scip)) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, consneg, subvars[i], -1.0) );
            SCIP_CALL( SCIPaddCoefLinear(subscip, consneg, eps, 1.0) );
            SCIP_CALL( SCIPaddCons(subscip, consneg) );
            SCIP_CALL( SCIPreleaseCons(subscip, &consneg) );

            /* release the variables */
            SCIP_CALL( SCIPreleaseVar(subscip, &eps) );
         }
      }
   }

   /* add and release the constraint representing the original objective function */
   if( objcons != NULL )
   {
      SCIP_CALL( SCIPaddCons(subscip, objcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &objcons) );
   }

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEUR*            heur,               /**< Completesol heuristic structure */
   SCIP_SOL*             subsol,             /**< solution of the subproblem or the partial */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   int nvars;                                /* the original problem's number of variables */
   SCIP_SOL* newsol;                         /* solution to be created for the original problem */
   int v;

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

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real solval = SCIPgetSolVal(subscip, subsol, subvars[v]);

      assert(!SCIPisInfinity(subscip, solval) && !SCIPisInfinity(subscip, -solval));
      assert(solval != SCIP_UNKNOWN); /*lint !e777*/

      SCIP_CALL( SCIPsetSolVal(scip, newsol, vars[v], solval) );
   }

   /* try to add new solution to SCIP and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   return SCIP_OKAY;
}

/** perform a probing bound change or fixes the variable */
static
SCIP_RETCODE chgProbingBound(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             newval,             /**< new bound */
   SCIP_BRANCHDIR        branchdir,          /**< bound change direction */
   SCIP_Bool*            success             /**< pointer to store whether the bound could be tightened */
   )
{
   SCIP_Real ub;
   SCIP_Real lb;

   assert(scip != NULL);
   assert(var != NULL);

   (*success) = FALSE;

   ub = SCIPvarGetUbLocal(var);
   lb = SCIPvarGetLbLocal(var);

   switch (branchdir) {
   case SCIP_BRANCHDIR_DOWNWARDS:
      if( SCIPisLT(scip, newval, ub) && SCIPisGE(scip, newval, lb) )
      {
         SCIP_CALL( SCIPchgVarUbProbing(scip, var, newval) );
         (*success) = TRUE;
      }
      break;
   case SCIP_BRANCHDIR_UPWARDS:
      if( SCIPisLE(scip, newval, ub) && SCIPisGT(scip, newval, lb) )
      {
         SCIP_CALL( SCIPchgVarLbProbing(scip, var, newval) );
         (*success) = TRUE;
      }
      break;
   case SCIP_BRANCHDIR_FIXED:
      if( SCIPisLE(scip, newval, ub) && SCIPisGE(scip, newval, lb) )
      {
         SCIP_CALL( SCIPfixVarProbing(scip, var, newval) );
         (*success) = TRUE;
      }
      break;
   default:
      return SCIP_INVALIDDATA;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** tries variables bound changes guided by the given solution */
static
SCIP_RETCODE tightenVariables(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic's private data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   SCIP_SOL*             sol,                /**< solution to guide the bound changes */
   SCIP_Bool*            tightened,          /**< array to store if variable bound could be tightened */
   SCIP_Bool*            success             /**< pointer to store the success */
   )
{
#ifndef NDEBUG
   SCIP_Bool incontsection;
#endif
   SCIP_Bool abortearly;
   SCIP_Bool cutoff;
   SCIP_Bool probingsuccess;
   SCIP_Longint ndomreds;
   SCIP_Longint ndomredssum;
   int nbndtightenings;
   int v;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(vars != NULL);
   assert(nvars >= 0);
   assert(sol != NULL);
   assert(tightened != NULL);

   assert(SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_PARTIAL);

   SCIPdebugMsg(scip, "> start probing along the solution values\n");

   *success = TRUE;
   abortearly = FALSE;
   nbndtightenings = 0;
   ndomredssum = 0;
#ifndef NDEBUG
   incontsection = FALSE;
#endif

   /* there is at least one integral variable; open one probing node for all non-continuous variables */
   if( nvars - SCIPgetNContVars(scip) > 0 )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
   }

   for( v = 0; v < nvars && !abortearly; v++ )
   {
      SCIP_Real solval;

      assert(SCIPvarIsActive(vars[v]));

      cutoff = FALSE;
      ndomreds = 0;

#ifndef NDEBUG
      incontsection |= (!SCIPvarIsIntegral(vars[v])); /*lint !e514*/
      assert(!incontsection || !SCIPvarIsIntegral(vars[v]));
#endif

      /* return if we have found enough domain reductions tightenings */
      if( ndomredssum > 0.3*nvars )
         break;

      solval = SCIPgetSolVal(scip, sol, vars[v]);

      /* skip unknown variables */
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         continue;
      assert(!SCIPisInfinity(scip, solval) && !SCIPisInfinity(scip, -solval));

      /* variable is binary or integer */
      if( SCIPvarIsIntegral(vars[v]) )
      {
         /* the solution value is integral, try to fix them */
         if( SCIPisIntegral(scip, solval) )
         {
            SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_FIXED, &probingsuccess) );
            tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
            ++nbndtightenings;

#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "> fix variable <%s> = [%g,%g] to %g \n", SCIPvarGetName(vars[v]),
                     SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]), solval);
#endif
         }
         else
         {
            SCIP_Real ub = SCIPceil(scip, solval) + 1.0;
            SCIP_Real lb = SCIPfloor(scip, solval) - 1.0;

            /* try tightening of upper bound */
            if( SCIPisLT(scip, ub, SCIPvarGetUbLocal(vars[v])) )
            {
               SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_DOWNWARDS, &probingsuccess) );
               tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
               ++nbndtightenings;

#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "> tighten upper bound of variable <%s>: %g to %g\n", SCIPvarGetName(vars[v]),
                     SCIPvarGetUbGlobal(vars[v]), ub);
#endif
            }

            /* try tightening of lower bound */
            if( SCIPisGT(scip, lb, SCIPvarGetLbLocal(vars[v])) )
            {
               SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_UPWARDS, &probingsuccess) );
               tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
               ++nbndtightenings;

#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "> tighten lower bound of variable <%s>: %g to %g\n", SCIPvarGetName(vars[v]),
                     SCIPvarGetLbGlobal(vars[v]), ub);
#endif
            }
         }
      }
      /* variable is continuous */
      else
      {
         /* fix to lb or ub */
         if( SCIPisEQ(scip, solval, SCIPvarGetLbLocal(vars[v])) || SCIPisEQ(scip, solval, SCIPvarGetUbLocal(vars[v])) )
         {
            /* open a new probing node */
            if( SCIPgetProbingDepth(scip) < SCIP_MAXTREEDEPTH-10 )
            {
               SCIP_CALL( SCIPnewProbingNode(scip) );

               SCIP_CALL( chgProbingBound(scip, vars[v], solval, SCIP_BRANCHDIR_FIXED, &probingsuccess) );

               /* skip propagation if the bound could not be changed, e.g., already tightened due to previous
                * domain propagation
                */
               if( probingsuccess )
               {
                  SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, &ndomreds) );
               }

               if( cutoff )
               {
                  ndomreds = 0;
                  SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
               }
               else
               {
                  assert(SCIPvarGetProbindex(vars[v]) >= 0);
                  tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
                  ++nbndtightenings;
#ifdef SCIP_MORE_DEBUG
                  SCIPdebugMsg(scip, "> fix variable <%s> = [%g,%g] to %g (ndomreds=%lld)\n", SCIPvarGetName(vars[v]),
                        SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]), solval, ndomreds);
#endif
               }
            }
            else
               /* abort probing */
               abortearly = TRUE;
         }
         else
         {
            SCIP_Real offset;
            SCIP_Real newub = SCIPvarGetUbGlobal(vars[v]);
            SCIP_Real newlb = SCIPvarGetLbGlobal(vars[v]);

            /* both bound are finite */
            if( !SCIPisInfinity(scip, -newlb) && !SCIPisInfinity(scip, newub) )
               offset = REALABS(heurdata->boundwidening * (newub-newlb));
            else
            {
               /* if one bound is finite, widen bound w.r.t. solution value and finite bound */
               if( !SCIPisInfinity(scip, -newlb) )
                  offset = REALABS(heurdata->boundwidening * (solval-newlb));
               else
               {
                  assert(!SCIPisInfinity(scip, newub));
                  offset = REALABS(heurdata->boundwidening * (newub-solval));
               }
            }

            /* update bounds */
            newub = SCIPceil(scip, solval) + offset;
            newlb = SCIPfloor(scip, solval) - offset;

            /* try tightening of upper bound */
            if( SCIPisLT(scip, newub, SCIPvarGetUbLocal(vars[v])) )
            {
               /* open a new probing node */
               if( SCIPgetProbingDepth(scip) < SCIP_MAXTREEDEPTH-10 )
               {
                  SCIP_CALL( SCIPnewProbingNode(scip) );
                  SCIP_CALL( chgProbingBound(scip, vars[v], newub, SCIP_BRANCHDIR_DOWNWARDS, &probingsuccess) );

                  /* skip propagation if the bound could not be changed, e.g., already tightened due to previous
                   * domain propagation
                   */
                  if( probingsuccess )
                  {
                     SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, &ndomreds) );
                  }

                  if( cutoff )
                  {
                     ndomreds = 0;

                     /* backtrack to last feasible probing node */
                     SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );

                     /* we can tighten the lower bound by newub */
                     SCIP_CALL( chgProbingBound(scip, vars[v], newub, SCIP_BRANCHDIR_UPWARDS, &probingsuccess) );

                     /* propagate the new bound */
                     SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, &ndomreds) );

                     /* there is no feasible solution w.r.t. the current bounds */
                     if( cutoff )
                     {
                        SCIPdebugMsg(scip, "> subproblem is infeasible within the local bounds\n");
                        *success = FALSE;
                        return SCIP_OKAY;
                     }
#ifdef SCIP_MORE_DEBUG
                     SCIPdebugMsg(scip, "> tighten lower bound of variable <%s>: %g to %g\n",
                           SCIPvarGetName(vars[v]), SCIPvarGetLbGlobal(vars[v]), newub);
#endif
                  }
                  else
                  {
                     assert(SCIPvarGetProbindex(vars[v]) >= 0);
                     tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
                     ++nbndtightenings;
#ifdef SCIP_MORE_DEBUG
                     SCIPdebugMsg(scip, "> tighten upper bound of variable <%s>: %g to %g (ndomreds=%lld)\n",
                           SCIPvarGetName(vars[v]), SCIPvarGetUbGlobal(vars[v]), newub, ndomreds);
#endif
                  }
               }
               else
                  /* abort probing */
                  abortearly = TRUE;
            }

            /* try tightening of lower bound */
            if( SCIPisGT(scip, newlb, SCIPvarGetLbLocal(vars[v])) )
            {
               /* open a new probing node */
               if( SCIPgetProbingDepth(scip) < SCIP_MAXTREEDEPTH-10 )
               {
                  SCIP_CALL( SCIPnewProbingNode(scip) );
                  SCIP_CALL( chgProbingBound(scip, vars[v], newlb, SCIP_BRANCHDIR_UPWARDS, &probingsuccess) );

                  /* skip propagation if the bound could not be changed, e.g., already tightened due to previous
                   * domain propagation
                   */
                  if( probingsuccess )
                  {
                     SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, &ndomreds) );
                  }

                  if( cutoff )
                  {
                     ndomreds = 0;

                     /* backtrack to last feasible probing node */
                     SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );

                     /* we can tighten the upper bound by newlb */
                     SCIP_CALL( chgProbingBound(scip, vars[v], newlb, SCIP_BRANCHDIR_DOWNWARDS, &probingsuccess) );

                     /* propagate the new bound */
                     SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, &ndomreds) );

                     /* there is no feasible solution w.r.t. the current bounds */
                     if( cutoff )
                     {
                        SCIPdebugMsg(scip, "> subproblem is infeasible within the local bounds\n");
                        *success = FALSE;
                        return SCIP_OKAY;
                     }
#ifdef SCIP_MORE_DEBUG
                     SCIPdebugMsg(scip, "> tighten upper bound of variable <%s>: %g to %g\n",
                           SCIPvarGetName(vars[v]), SCIPvarGetUbGlobal(vars[v]), newlb);
#endif
                  }
                  else
                  {
                     assert(SCIPvarGetProbindex(vars[v]) >= 0);
                     tightened[SCIPvarGetProbindex(vars[v])] = TRUE;
                     ++nbndtightenings;
#ifdef SCIP_MORE_DEBUG
                     SCIPdebugMsg(scip, "> tighten lower bound of variable <%s>: %g to %g (ndomreds=%lld)\n",
                           SCIPvarGetName(vars[v]), SCIPvarGetLbGlobal(vars[v]), newlb, ndomreds);
#endif
                  }
               }
               else
                  /* abort probing */
                  abortearly = TRUE;
            }
         }
      }

      ndomredssum += ndomreds;
   }

   SCIPdebugMsg(scip, "> found %d bound tightenings and %lld induced domain reductions (abort=%u).\n", nbndtightenings,
         ndomredssum, abortearly);

   return SCIP_OKAY;
}

/* setup and solve the sub-SCIP */
static
SCIP_RETCODE setupAndSolve(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic's private data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_SOL*             partialsol,         /**< partial solution */
   SCIP_Bool*            tightened           /**< array to store whether a variable was already tightened */
   )
{
   SCIP_HASHMAP* varmapf;
   SCIP_VAR** vars;
   SCIP_VAR** subvars = NULL;
   SCIP_EVENTHDLR* eventhdlr;
   int nvars;
   int i;

   SCIP_SOL** subsols;
   int nsubsols;

   SCIP_Bool valid;
   SCIP_Bool success;
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);
   assert(partialsol != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapf, SCIPblkmem(subscip), nvars) );

   eventhdlr = NULL;
   valid = FALSE;

   /* copy complete SCIP instance */
   SCIP_CALL( SCIPcopyConsCompression(scip, subscip, varmapf, NULL, "completesol", NULL, NULL, 0, FALSE, FALSE, TRUE, &valid) );
   SCIPdebugMsg(scip, "Copying the SCIP instance returned with valid=%d.\n", valid);

   /* create event handler for LP events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecCompletesol, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* allocate memory to align the SCIP and the sub-SCIP variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* map all variables */
   for( i = 0; i < nvars; i++ )
   {
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapf, vars[i]);
     assert(subvars[i] != NULL);
   }

   /* free hash map */
   SCIPhashmapFree(&varmapf);

   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   SCIP_CALL( createSubproblem(scip, subscip, heurdata, subvars, partialsol, tightened, &success) );
   if( !success )
   {
      SCIPdebugMsg(scip, "Error while creating completesol subproblem w.r.t. partial solution <%p>.\n", (void*)partialsol);
      goto TERMINATE;
   }
   SCIPdebugMsg(scip, "Completesol subproblem: %d vars, %d cons\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip));

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", SCIP_VERBLEVEL_FULL) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", -1) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", (int) SCIP_VERBLEVEL_NONE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* set limits for the subproblem */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   heurdata->nodelimit = heurdata->maxnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsols) );

   /* limit the number of LP iterations */
   SCIP_CALL( SCIPsetLongintParam(subscip, "lp/iterlim", heurdata->maxlpiter) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "lp/rootiterlim", heurdata->maxlpiter) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", FALSE) );
   }

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   /* solve the subproblem */
   SCIPdebugMsg(scip, "solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->maxnodes);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */

   retcode = SCIPpresolve(subscip);

   /* errors in presolving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while presolving subproblem in %s heuristic; sub-SCIP terminated with code <%d>\n", HEUR_NAME, retcode);

      SCIPABORT(); /*lint --e{527}*/

      goto TERMINATE;
   }

   if( SCIPgetStage(subscip) == SCIP_STAGE_PRESOLVED )
   {
      SCIPdebugMsg(scip, "presolved instance has bin=%d, int=%d, cont=%d variables\n",
            SCIPgetNBinVars(subscip), SCIPgetNIntVars(subscip), SCIPgetNContVars(subscip));

      /* check whether the presolved instance is small enough */
      if( heurdata->maxcontvars >= 0 && SCIPgetNContVars(subscip) > heurdata->maxcontvars )
      {
         SCIPdebugMsg(scip, "presolved instance has too many continuous variables (maxcontvars: %d)\n", heurdata->maxcontvars);
         goto TERMINATE;
      }

      /* set node limit of 1 if the presolved problem is an LP, otherwise we would start branching if an LP iteration
       * limit was set by the user.
       */
      if( !SCIPisNLPEnabled(subscip) && SCIPgetNContVars(subscip) == SCIPgetNVars(subscip) )
      {
         SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", 1LL) );
      }

      retcode = SCIPsolve(subscip);

      /* errors in solving the subproblem should not kill the overall solving process;
       * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
      if( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving subproblem in %s heuristic; sub-SCIP terminated with code <%d>\n", HEUR_NAME, retcode);

         SCIPABORT(); /*lint --e{527}*/

         goto TERMINATE;
      }
   }

   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   success = FALSE;
   for( i = 0; i < nsubsols && (!success || heurdata->addallsols); i++ )
   {
      SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      if( success )
         *result = SCIP_FOUNDSOL;
   }

   SCIPstatisticPrintf("%s statistic: fixed %6.3f integer variables, needed %6.1f seconds, %" SCIP_LONGINT_FORMAT " nodes, solution %10.4f found at node %" SCIP_LONGINT_FORMAT "\n",
      HEUR_NAME, 0.0, SCIPgetSolvingTime(subscip), SCIPgetNNodes(subscip), success ? SCIPgetPrimalbound(scip) : SCIPinfinity(scip),
      nsubsols > 0 ? SCIPsolGetNodenum(SCIPgetBestSol(subscip)) : -1 );

  TERMINATE:
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}

/** main procedure of the completesol heuristic, creates and solves a sub-SCIP */
static
SCIP_RETCODE applyCompletesol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic's private data structure */
   SCIP_RESULT*          result,             /**< result data structure */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem */
   SCIP_SOL*             partialsol          /**< partial solution */
   )
{
   SCIP* subscip;
   SCIP_VAR** vars;
   SCIP_Bool* tightened;
   SCIP_Bool success;
   SCIP_RETCODE retcode;
   int nvars;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);
   assert(result != NULL);
   assert(partialsol != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "+---+ Start Completesol heuristic +---+\n");

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if( !success )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get variable data */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* get buffer memory and initialize it to FALSE */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &tightened, nvars) );

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_CALL( tightenVariables(scip, heurdata, vars, nvars, partialsol, tightened, &success) );

   if( !success )
      goto ENDPROBING;

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   retcode = setupAndSolve(scip, subscip, heur, heurdata, result, nstallnodes, partialsol, tightened);

   /* free subproblem */
   SCIP_CALL( SCIPfree(&subscip) );

   SCIP_CALL( retcode );

  ENDPROBING:
   SCIPfreeBufferArray(scip, &tightened);
   SCIP_CALL( SCIPendProbing(scip) );



   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyCompletesol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurCompletesol(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeCompletesol)
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

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecCompletesol)
{/*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   SCIP_SOL** partialsols;
   SCIP_Longint nstallnodes;
   int npartialsols;
   int nunknown;
   int nfracints;
   int nvars;
   int s;
   int v;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DIDNOTRUN;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* do not run after restart */
   if( SCIPgetNRuns(scip) > 1 )
      return SCIP_OKAY;

   /* check whether we want to run before presolving */
   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL && !heurdata->beforepresol )
      return SCIP_OKAY;

   /* only run before root node */
   if( heurtiming == SCIP_HEURTIMING_BEFORENODE && SCIPgetCurrentNode(scip) != SCIPgetRootNode(scip) )
      return SCIP_OKAY;

   /* get variable data and return of no variables are left in the problem */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   if( nvars == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward Completesol if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMsg(scip, "skipping Complete: nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n",
         nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* check the number of variables with unknown value and continuous variables with fractional value */
   nfracints = 0;

   /* get all partial sols */
   npartialsols = SCIPgetNPartialSols(scip);
   partialsols = SCIPgetPartialSols(scip);

   /* loop over all partial solutions */
   for( s = 0; s < npartialsols; s++ )
   {
      SCIP_SOL* sol;
      SCIP_Real solval;
      SCIP_Real unknownrate;

      sol = partialsols[s];
      assert(sol != NULL);
      assert(SCIPsolIsPartial(sol));

      nunknown = 0;
      /* loop over all variables */
      for( v = 0; v < nvars; v++ )
      {
         assert(SCIPvarIsActive(vars[v]));

         /* skip continuous variables if they should ignored */
         if( !SCIPvarIsIntegral(vars[v]) && heurdata->ignorecont )
            continue;

         solval = SCIPgetSolVal(scip, sol, vars[v]);

         /* we only want to count variables that are unfixed after the presolving */
         if( solval == SCIP_UNKNOWN ) /*lint !e777*/
            ++nunknown;
         else if( SCIPvarIsIntegral(vars[v]) && !SCIPisIntegral(scip, solval) )
            ++nfracints;
      }

      if( heurdata->ignorecont )
         unknownrate = nunknown/((SCIP_Real)SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip));
      else
         unknownrate = nunknown/((SCIP_Real)nvars);
      SCIPdebugMsg(scip, "%d (rate %.4f) unknown solution values\n", nunknown, unknownrate);

      /* run the heuristic, if not too many unknown variables exist */
      if( unknownrate > heurdata->maxunknownrate )
         continue;

      /* all variables have a finite/known solution value all integer variables have an integral solution value,
       * and there are no continuous variables
       * in the sub-SCIP, all variables would be fixed, so create a new solution without solving a sub-SCIP
       */
      if( nunknown == 0 && nfracints == 0 && SCIPgetNContVars(scip) == 0 && SCIPgetNImplVars(scip) == 0 )
      {
         SCIP_VAR** origvars;
         SCIP_SOL* newsol;
         SCIP_Bool stored;
         int norigvars;

         origvars = SCIPgetOrigVars(scip);
         norigvars = SCIPgetNOrigVars(scip);

         SCIP_CALL( SCIPcreateOrigSol(scip, &newsol, heur) );

         for( v = 0; v < norigvars; v++ )
         {
            solval = SCIPgetSolVal(scip, sol, origvars[v]);
            assert(solval != SCIP_UNKNOWN); /*lint !e777*/

            SCIP_CALL( SCIPsetSolVal(scip, newsol, origvars[v], solval) );
         }

         SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );
         if( stored )
            *result = SCIP_FOUNDSOL;
      }
      else
      {
         /* run the heuristic */
         SCIP_CALL( applyCompletesol(scip, heur, heurdata, result, nstallnodes, sol) );
      }
   }

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the completesol primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCompletesol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create completesol primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   assert(heurdata != NULL);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecCompletesol, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyCompletesol) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeCompletesol) );

   /* add completesol primal heuristic parameters */

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxunknownrate",
         "maximal rate of unknown solution values",
         &heurdata->maxunknownrate, FALSE, DEFAULT_MAXUNKRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/addallsols",
         "should all subproblem solutions be added to the original SCIP?",
         &heurdata->addallsols, TRUE, DEFAULT_ADDALLSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/objweight",
         "weight of the original objective function (1: only original objective)",
         &heurdata->objweight, TRUE, DEFAULT_OBJWEIGHT, DEFAULT_MINOBJWEIGHT, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/boundwidening",
         "bound widening factor applied to continuous variables (0: fix variables to given solution values, 1: relax to global bounds)",
         &heurdata->boundwidening, TRUE, DEFAULT_BOUNDWIDENING, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which the incumbent should be improved at least",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/ignorecont",
         "should number of continuous variables be ignored?",
         &heurdata->ignorecont, FALSE, DEFAULT_IGNORECONT, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/solutions",
         "heuristic stops, if the given number of improving solutions were found (-1: no limit)",
         &heurdata->bestsols, FALSE, DEFAULT_BESTSOLS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxproprounds",
         "maximal number of iterations in propagation (-1: no limit)",
         &heurdata->maxproprounds, FALSE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/beforepresol",
         "should the heuristic run before presolving?",
         &heurdata->beforepresol, FALSE, DEFAULT_BEFOREPRESOL, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxlpiter",
         "maximal number of LP iterations (-1: no limit)",
         &heurdata->maxlpiter, FALSE, DEFAULT_MAXLPITER, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxcontvars",
         "maximal number of continuous variables after presolving",
         &heurdata->maxcontvars, FALSE, DEFAULT_MAXCONTVARS, -1, INT_MAX, NULL, NULL) );


   return SCIP_OKAY;
}
