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

/**@file   heur_dins.c
 * @brief  DINS primal heuristic (according to Ghosh)
 * @author Timo Berthold
 * @author Robert Waniek
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/heur_dins.h"

#define HEUR_NAME             "dins"
#define HEUR_DESC             "distance induced neighborhood search by Ghosh"
#define HEUR_DISPCHAR         'D'
#define HEUR_PRIORITY         -1105000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NODESOFS      5000LL    /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINNODES      50LL      /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which DINS should at least improve the incumbent          */
#define DEFAULT_NODESQUOT     0.05      /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_LPLIMFAC      1.5       /* factor by which the limit on the number of LP depends on the node limit  */
#define DEFAULT_MINFIXINGRATE 0.3       /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_NWAITINGNODES 200LL       /* number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_NEIGHBORHOODSIZE  18    /* radius of the incumbents neighborhood to be searched                */
#define DEFAULT_SOLNUM        5         /* number of pool-solutions to be checked for flag array update        */
#define DEFAULT_USELPROWS     FALSE     /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip        */

#define DEFAULT_BESTSOLLIMIT    3       /* limit on number of improving incumbent solutions in sub-CIP            */
#define DEFAULT_USEUCT        FALSE     /* should uct node selection be used at the beginning of the search?     */


/* event handler properties */
#define EVENTHDLR_NAME         "Dins"
#define EVENTHDLR_DESC         "LP event handler for " HEUR_NAME " heuristic"

/*
 * Data structures
 */

/** DINS primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Longint          nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   SCIP_Real             minimprove;         /**< factor by which DINS should at least improve the incumbent          */
   SCIP_Longint          usednodes;          /**< nodes already used by DINS in earlier calls                         */
   SCIP_Longint          lastnsolsfound;     /**< total number of found solutions at previous execution of DINS       */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_Real             nodelimit;          /**< the nodelimit employed in the current sub-SCIP, for the event handler*/
   SCIP_Real             lplimfac;           /**< factor by which the limit on the number of LP depends on the node limit */
   int                   neighborhoodsize;   /**< radius of the incumbent's neighborhood to be searched               */
   SCIP_Bool*            delta;              /**< stores whether a variable kept its value from root LP all the time  */
   int                   deltalength;        /**< if there are no binary variables, we need no flag array             */
   int                   solnum;             /**< number of pool-solutions to be checked for flag array update        */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?        */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?
                                              */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP            */
   SCIP_Bool             useuct;             /**< should uct node selection be used at the beginning of the search?  */
};


/*
 * Local methods
 */

/** compute tightened bounds for integer variables depending on how much the LP and the incumbent solution values differ */
static
void computeIntegerVariableBounds(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP_VAR*             var,                /**< the variable for which bounds should be computed */
   SCIP_Real*            lbptr,              /**< pointer to store the lower bound in the DINS sub-SCIP */
   SCIP_Real*            ubptr               /**< pointer to store the upper bound in the DINS sub-SCIP */
   )
{
   SCIP_Real mipsol;
   SCIP_Real lpsol;

   SCIP_Real lbglobal;
   SCIP_Real ubglobal;
   SCIP_SOL* bestsol;

   /* get the bounds for each variable */
   lbglobal = SCIPvarGetLbGlobal(var);
   ubglobal = SCIPvarGetUbGlobal(var);

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);
   /* get the current LP solution for each variable */
   lpsol = SCIPvarGetLPSol(var);

   /* get the current MIP solution for each variable */
   bestsol = SCIPgetBestSol(scip);
   mipsol = SCIPgetSolVal(scip, bestsol, var);


   /* if the solution values differ by 0.5 or more, the variable is rebounded, otherwise it is just copied */
   if( REALABS(lpsol - mipsol) >= 0.5 )
   {
      SCIP_Real range;

      *lbptr = lbglobal;
      *ubptr = ubglobal;

      /* create a equally sized range around lpsol for general integers: bounds are lpsol +- (mipsol-lpsol) */
      range = 2*lpsol-mipsol;

      if( mipsol >= lpsol )
      {
         range = SCIPfeasCeil(scip, range);
         *lbptr = MAX(*lbptr, range);

         /* when the bound new upper bound is equal to the current MIP solution, we set both bounds to the integral bound (without eps) */
         if( SCIPisFeasEQ(scip, mipsol, *lbptr) )
            *ubptr = *lbptr;
         else
            *ubptr = mipsol;
      }
      else
      {
         range = SCIPfeasFloor(scip, range);
         *ubptr = MIN(*ubptr, range);

         /* when the bound new upper bound is equal to the current MIP solution, we set both bounds to the integral bound (without eps) */
         if( SCIPisFeasEQ(scip, mipsol, *ubptr) )
            *lbptr = *ubptr;
         else
            *lbptr = mipsol;
      }

      /* the global domain of variables might have been reduced since incumbent was found: adjust lb and ub accordingly */
      *lbptr = MAX(*lbptr, lbglobal);
      *ubptr = MIN(*ubptr, ubglobal);
   }
   else
   {
      /* the global domain of variables might have been reduced since incumbent was found: adjust it accordingly */
      *lbptr = MAX(mipsol, lbglobal);
      *ubptr = MIN(mipsol, ubglobal);
   }
}

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE determineVariableFixings(
   SCIP*                 scip,               /**< SCIP data structure of the original problem                    */
   SCIP_HEUR*            heur,               /**< DINS heuristic structure                                       */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure                                       */
   SCIP_VAR**            vars,               /**< variables of the original problem                              */
   SCIP_VAR**            fixedvars,          /**< array to store variables that should be fixed in the sub-SCIP  */
   SCIP_Real*            fixedvals,          /**< array to store fixing values for fixed variables               */
   int                   nbinvars,           /**< number of binary variables of problem and subproblem           */
   int                   nintvars,           /**< number of general integer variables of problem and subproblem  */
   int*                  binfixings,         /**< pointer to store number of binary variables that should be fixed */
   int*                  intfixings          /**< pointer to store number of integer variables that should be fixed */
   )
{
   SCIP_SOL* bestsol;
   SCIP_SOL** sols;
   SCIP_Bool* delta;
   int i;
   int nsols;
   SCIP_Longint nsolsfound;
   int checklength;
   int nfixedvars;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(fixedvars != NULL);
   assert(fixedvals != NULL);
   assert(binfixings != NULL);
   assert(intfixings != NULL);
   assert(heur != NULL);

   /* get the best MIP-solution known so far */
   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   /* get solution pool and number of solutions in pool */
   sols = SCIPgetSols(scip);
   nsols = SCIPgetNSols(scip);
   nsolsfound = SCIPgetNSolsFound(scip);
   checklength = MIN(nsols, heurdata->solnum);
   assert(sols != NULL);
   assert(nsols > 0);

   /* if new binary variables have been created, e.g., due to column generation, reallocate the delta array */
   if( heurdata->deltalength < nbinvars )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, nbinvars);
      assert(newsize >= nbinvars);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &heurdata->delta, heurdata->deltalength, newsize) );

      /* initialize new part of delta array */
      for( i = heurdata->deltalength; i < newsize; i++ )
         heurdata->delta[i] = TRUE;

      heurdata->deltalength = newsize;
   }

   delta = heurdata->delta;
   /* fixing for binary variables */
   /* hard fixing for some with mipsol(s)=lpsolval=rootlpsolval and preparation for soft fixing for the remaining */
   nfixedvars = 0;
   *intfixings = *binfixings = 0;
   for( i = 0; i < nbinvars; i++ )
   {
      SCIP_Real lpsolval;
      SCIP_Real mipsolval;
      SCIP_Real rootlpsolval;
      int j;

      /* get the current LP solution for each variable */
      lpsolval = SCIPvarGetLPSol(vars[i]);
      /* get the current MIP solution for each variable */
      mipsolval = SCIPgetSolVal(scip, bestsol, vars[i]);
      /* get the root LP solution for each variable */
      rootlpsolval = SCIPvarGetRootSol(vars[i]);

      if( SCIPisFeasEQ(scip, lpsolval, mipsolval) && SCIPisFeasEQ(scip, mipsolval, rootlpsolval) )
      {
         /* update delta */
         if( nsols > 1 && heurdata->lastnsolsfound != nsolsfound && delta[i] ) /* no need to update delta[i] if already FALSE */
         {
            /* no need to update delta[i] if already FALSE or sols[i] already checked on previous run or worse than DINS-solution of last run */
            for( j = 1; delta[i] && j < checklength && SCIPgetSolHeur(scip, sols[j]) != heur ; j++ )
            {
               SCIP_Real solval;
               solval = SCIPgetSolVal(scip, sols[j], vars[i]);
               delta[i] = delta[i] && SCIPisFeasEQ(scip, mipsolval, solval);
            }
         }

         /* hard fixing if rootlpsolval=nodelpsolval=mipsolval(s) and delta (is TRUE) */
         if( delta[i] )
         {
            fixedvars[nfixedvars] = vars[i];
            fixedvals[nfixedvars] = mipsolval;
            ++nfixedvars;
         }
      }
   }

   *binfixings = nfixedvars;

   /* store the number of found solutions for next run */
   heurdata->lastnsolsfound = nsolsfound;

   /* compute a tighter bound for the variable in the subproblem;  */
   for( i = nbinvars; i < nbinvars + nintvars; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      computeIntegerVariableBounds(scip, vars[i], &lb, &ub);

      /* hard fixing if heuristic bounds coincide */
      if( ub - lb < 0.5 )
      {
         fixedvars[(nfixedvars)] = vars[i];
         fixedvals[(nfixedvars)] = lb;
         ++nfixedvars;
      }
   }

   *intfixings = nfixedvars - *binfixings;

   return SCIP_OKAY;
}

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE reboundIntegerVariables(
   SCIP*                 scip,               /**< SCIP data structure of the original problem                    */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem       */
   SCIP_VAR**            vars,               /**< variables of the original problem                              */
   SCIP_VAR**            subvars,            /**< variables of the DINS sub-SCIP                                 */
   int                   nbinvars,           /**< number of binary variables of problem and subproblem           */
   int                   nintvars            /**< number of general integer variables of problem and subproblem  */
   )
{
   int i;
   /* compute a tighter bound for the variable in the subproblem;  */
   for( i = nbinvars; i < nbinvars + nintvars; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      computeIntegerVariableBounds(scip, vars[i], &lb, &ub);

      /* change variable bounds in the DINS subproblem; if bounds coincide, variable should already be fixed */
      if( ub - lb >= 0.5 )
      {
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], lb) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], ub) );
      }
      else
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(subvars[i]), SCIPvarGetUbGlobal(subvars[i])));
      }
   }

   return SCIP_OKAY;
}

/** create the extra constraint of local branching and add it to subscip */
static
SCIP_RETCODE addLocalBranchingConstraint(
   SCIP*                 scip,               /**< SCIP data structure of the original problem */
   SCIP*                 subscip,            /**< SCIP data structure of the subproblem       */
   SCIP_VAR**            subvars,            /**< variables of the subproblem                 */
   SCIP_HEURDATA*        heurdata            /**< heuristic's data structure                  */
   )
{
   SCIP_CONS* cons;                     /* local branching constraint to create          */
   SCIP_VAR** vars;
   SCIP_SOL* bestsol;

   SCIP_Real* consvals;
   SCIP_Real solval;
   SCIP_Real lhs;
   SCIP_Real rhs;

   char consname[SCIP_MAXSTRLEN];

   int nbinvars;
   int i;

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_dinsLBcons", SCIPgetProbName(scip));

   /* get the data of the variables and the best solution */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, NULL, NULL, NULL) );
   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nbinvars) );

   /* set initial left and right hand sides of local branching constraint */
   lhs = 0.0;
   rhs = (SCIP_Real) heurdata->neighborhoodsize;

   /* create the distance function of the binary variables (to incumbent solution) */
   for( i = 0; i < nbinvars; i++ )
   {
      assert(SCIPvarGetType(subvars[i]) == SCIP_VARTYPE_BINARY);
      if( SCIPvarGetUbGlobal(subvars[i]) - SCIPvarGetLbGlobal(subvars[i]) < 0.5 )
      {
         consvals[i]=0.0;
         continue;
      }

      solval = SCIPgetSolVal(scip, bestsol, vars[i]);
      assert(SCIPisFeasIntegral(scip, solval));

      /* is variable i part of the binary support of the current solution? */
      if( SCIPisFeasEQ(scip, solval, 1.0) )
      {
         consvals[i] = -1.0;
         rhs -= 1.0;
         lhs -= 1.0;
      }
      else
         consvals[i] = 1.0;
   }

   /* creates local branching constraint and adds it to subscip */
   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, consname, nbinvars, subvars, consvals,
         lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   /* free local memory */
   SCIPfreeBufferArray(scip, &consvals);

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< DINS heuristic structure                            */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(heur != NULL);
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

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );
   if( *success )
   {
      SCIPdebugMsg(scip, "DINS successfully found new solution\n");
   }

   SCIPfreeBufferArray(scip, &subsolvals);
   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTEXEC(eventExecDins);

/** wrapper for the part of heuristic that runs a subscip. Wrapper is needed to avoid possible ressource leaks */
static
SCIP_RETCODE wrapperDins(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_HEUR*            heur,               /**< Heuristic pointer                                   */
   SCIP_HEURDATA*        heurdata,           /**< Heuristic's data                                    */
   SCIP_VAR**            vars,               /**< original problem's variables                        */
   SCIP_VAR**            fixedvars,          /**< Fixed variables of original SCIP                    */
   SCIP_Real*            fixedvals,          /**< Fixed values of original SCIP                       */
   SCIP_RESULT*          result,             /**< Result pointer                                      */
   int                   nvars,              /**< Number of variables                                 */
   int                   nbinvars,           /**< Number of binary variables in original SCIP         */
   int                   nintvars,           /**< Number of integer variables in original SCIP        */
   int                   binfixings,         /**< Number of binary fixing in original SCIP            */
   int                   intfixings,         /**< Number of integer fixings in original SCIP          */
   SCIP_Longint          nsubnodes           /**< Number of nodes in the subscip                      */
   )
{
   SCIP_VAR** subvars;                       /* variables of the subscip */
   SCIP_HASHMAP*  varmapfw;                  /* hashmap for mapping between vars of scip and subscip */
   SCIP_EVENTHDLR* eventhdlr;                /* event handler for LP events  */
   SCIP_Real upperbound;                     /* upperbound of the original SCIP */
   SCIP_Real cutoff;                         /* objective cutoff for the subproblem */

   SCIP_Bool success;

   int i;
   int nsubsols;

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

   success = FALSE;
   eventhdlr = NULL;

   /* create a problem copy as sub SCIP */
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "dins", fixedvars, fixedvals, binfixings + intfixings,
      heurdata->uselprows, heurdata->copycuts, &success, NULL) );

   /* create event handler for LP events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecDins, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMsg(scip, "Copying the SCIP instance was %ssuccessful.\n", success ? "" : "not ");

   SCIPdebugMsg(scip, "DINS subproblem: %d vars (%d binvars & %d intvars), %d cons\n",
      SCIPgetNVars(subscip), SCIPgetNBinVars(subscip) , SCIPgetNIntVars(subscip) , SCIPgetNConss(subscip));

   /* store subproblem variables that correspond to original variables */
   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* perform prepared softfixing for all unfixed vars if the number of unfixed vars is larger than the neighborhoodsize (otherwise it will be useless) */
   if( nbinvars - binfixings > heurdata->neighborhoodsize )
   {
      SCIP_CALL( addLocalBranchingConstraint(scip, subscip, subvars, heurdata) );
   }

   /* rebound integer variables if not all were fixed */
   if( intfixings < nintvars )
   {
      assert(nintvars > 0);
      SCIP_CALL( reboundIntegerVariables(scip, subscip, vars, subvars, nbinvars, nintvars) );
   }

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* set limits for the subproblem */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   heurdata->nodelimit = nsubnodes;
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nsubnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", MAX(10, nsubnodes/10)) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsollimit) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
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

   /* activate uct node selection at the top of the tree */
   if( heurdata->useuct && SCIPfindNodesel(subscip, "uct") != NULL && !SCIPisParamFixed(subscip, "nodeselection/uct/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/uct/stdpriority", INT_MAX/2) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
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

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 500) );
   }

   /* add an objective cutoff */
   assert(!SCIPisInfinity(scip, SCIPgetUpperbound(scip)));

   if( !SCIPisInfinity(scip, -1.0*SCIPgetLowerbound(scip)) )
   {
      cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip) + heurdata->minimprove * SCIPgetLowerbound(scip);
      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      cutoff = MIN(upperbound, cutoff);
   }
   else
   {
      if( SCIPgetUpperbound(scip) >= 0 )
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
      else
         cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      cutoff = MIN(upperbound, cutoff);
   }
   SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );

   /* catch LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPtransformProb(subscip) );
      SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );
   }

   /* solve the subproblem */
   SCIPdebugMsg(scip, "solving DINS sub-MIP with neighborhoodsize %d and maxnodes %" SCIP_LONGINT_FORMAT "\n", heurdata->neighborhoodsize, nsubnodes);

   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   SCIP_CALL_ABORT( SCIPsolve(subscip) );

   /* drop LP events of sub-SCIP */
   if( !heurdata->uselprows )
   {
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   heurdata->usednodes += SCIPgetNNodes(subscip);
   nsubsols = SCIPgetNSols(subscip);
   SCIPdebugMsg(scip, "DINS used %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT " nodes and found %d solutions\n", SCIPgetNNodes(subscip), nsubnodes, nsubsols);

   /* check, whether a  (new) solution was found */
   if( nsubsols > 0 )
   {
      SCIP_SOL** subsols;

      /* check, whether a solution was found; due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted */
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      }
      if( success )
         *result = SCIP_FOUNDSOL;
   }

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);

   return SCIP_OKAY;
}


/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecDins)
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


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyDins)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurDins(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeDins)
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


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolDins)
{
   SCIP_HEURDATA* heurdata;
   int i;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   heurdata->lastnsolsfound = 0;

   /* create flag array */
   heurdata->deltalength = SCIPgetNBinVars(scip);

   /* no binvars => no flag array needed */
   if( heurdata->deltalength > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->delta), heurdata->deltalength) );
      for( i = 0; i < heurdata->deltalength; i++ )
         heurdata->delta[i] = TRUE;
   }
   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolDins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free flag array if exist */
   if( heurdata->deltalength > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &(heurdata->delta), heurdata->deltalength);
   }
   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecDins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP* subscip;                            /* the subproblem created by DINS                               */
   SCIP_VAR** vars;                          /* variables of the original problem                            */
   SCIP_VAR** fixedvars;
   SCIP_Real* fixedvals;

   SCIP_Longint maxnnodes;                   /* maximum number of subnodes                                   */
   SCIP_Longint nsubnodes;                   /* nodelimit for subscip                                        */

   SCIP_RETCODE retcode;

   int nvars;                                /* number of variables in original SCIP                         */
   int nbinvars;                             /* number of binary variables in original SCIP                  */
   int nintvars;                             /* number of general integer variables in original SCIP         */
   int binfixings;
   int intfixings;

   SCIP_Bool success;                        /* used to store whether new solution was found or not          */

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* only call heuristic, if a CIP solution is at hand */
   if( SCIPgetNSols(scip) <= 0 )
      return SCIP_OKAY;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip, SCIPgetBestSol(scip)) < heurdata->nwaitingnodes )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* determine the node limit for the current process */
   maxnnodes = (SCIP_Longint) (heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward DINS if it succeeded often */
   maxnnodes = (SCIP_Longint) (maxnnodes * (1.0 + 2.0 * (SCIPheurGetNBestSolsFound(heur)+1.0) / (SCIPheurGetNCalls(heur) + 1.0)));

   /* count the setup costs for the sub-MIP as 100 nodes */
   maxnnodes -= 100 * SCIPheurGetNCalls(heur);
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nsubnodes = maxnnodes - heurdata->usednodes;
   nsubnodes = MIN(nsubnodes , heurdata->maxnodes);

   /* check whether we have enough nodes left to call sub problem solving */
   if( nsubnodes < heurdata->minnodes )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
     return SCIP_OKAY;

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   assert(nbinvars <= nvars);

   /* do not run heuristic if only continuous variables are present */
   if( nbinvars == 0 && nintvars == 0 )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   /* abort if no time is left or not enough memory to create a copy of SCIP */
   if( !success )
      return SCIP_OKAY;

   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, nbinvars + nintvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, nbinvars + nintvars) );

   /* collect variables that should be fixed in the DINS subproblem */
   binfixings = intfixings = 0;
   SCIP_CALL( determineVariableFixings(scip, heur, heurdata, vars, fixedvars, fixedvals, nbinvars, nintvars, &binfixings, &intfixings) );

   /* abort, if all integer variables were fixed (which should not happen for MIP),
    * but frequently happens for MINLPs using an LP relaxation */
   if( binfixings + intfixings == nbinvars + nintvars )
      goto TERMINATE;

   /* abort, if the amount of fixed variables is insufficient */
   if( (binfixings + intfixings) / (SCIP_Real)(MAX(nbinvars + nintvars, 1)) < heurdata->minfixingrate )
      goto TERMINATE;

   *result = SCIP_DIDNOTFIND;

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );


   retcode = wrapperDins(scip, subscip, heur, heurdata, vars, fixedvars, fixedvals, result, nvars, nbinvars, nintvars, binfixings, intfixings, nsubnodes);
   SCIP_CALL( SCIPfree(&subscip) );

   SCIP_CALL( retcode );

 TERMINATE:
    SCIPfreeBufferArray(scip, &fixedvals);
    SCIPfreeBufferArray(scip, &fixedvars);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the DINS primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurDins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Dins primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecDins, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyDins) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeDins) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolDins) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolDins) );

   /* add DINS primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, FALSE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/solnum",
         "number of pool-solutions to be checked for flag array update (for hard fixing of binary variables)",
         &heurdata->solnum, FALSE, DEFAULT_SOLNUM, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/neighborhoodsize",
         "radius (using Manhattan metric) of the incumbent's neighborhood to be searched",
         &heurdata->neighborhoodsize, FALSE, DEFAULT_NEIGHBORHOODSIZE, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which " HEUR_NAME " should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the limit on the number of LP depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "minimum percentage of integer variables that have to be fixable",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useuct",
         "should uct node selection be used at the beginning of the search?",
         &heurdata->useuct, TRUE, DEFAULT_USEUCT, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/bestsollimit",
         "limit on number of improving incumbent solutions in sub-CIP",
         &heurdata->bestsollimit, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
