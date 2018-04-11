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

/**@file   cons_components.c
 * @brief  constraint handler for handling independent components
 * @author Gerald Gamrath
 *
 * This constraint handler looks for independent components.
 */
/*#define DETAILED_OUTPUT*/
/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG*/
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_components.h"
#include "scip/debug.h"

#define CONSHDLR_NAME          "components"
#define CONSHDLR_DESC          "independent components constraint handler"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYPROP         TRUE /**< should propagation method be delayed, if other propagators found reductions? */

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_FINAL  /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define DEFAULT_MAXDEPTH             -1      /**< maximum depth of a node to run components detection (-1: disable component detection during solving) */
#define DEFAULT_MAXINTVARS          500      /**< maximum number of integer (or binary) variables to solve a subproblem directly in presolving (-1: no solving) */
#define DEFAULT_MINSIZE              50      /**< minimum absolute size (in terms of variables) to solve a component individually during branch-and-bound */
#define DEFAULT_MINRELSIZE          0.1      /**< minimum relative size (in terms of variables) to solve a component individually during branch-and-bound */
#define DEFAULT_NODELIMIT       10000LL      /**< maximum number of nodes to be solved in subproblems during presolving */
#define DEFAULT_INTFACTOR           1.0      /**< the weight of an integer variable compared to binary variables */
#define DEFAULT_FEASTOLFACTOR       1.0      /**< default value for parameter to increase the feasibility tolerance in all sub-SCIPs */

/*
 * Data structures
 */

/** data related to one problem (see below) */
typedef struct Problem PROBLEM;

/** data related to one component */
typedef struct Component
{
   PROBLEM*              problem;            /**< the problem this component belongs to */
   SCIP*                 subscip;            /**< sub-SCIP representing the component */
   SCIP_SOL*             workingsol;         /**< working solution for transferring solutions into the sub-SCIP */
   SCIP_VAR**            vars;               /**< variables belonging to this component (in complete problem) */
   SCIP_VAR**            subvars;            /**< variables belonging to this component (in subscip) */
   SCIP_VAR**            fixedvars;          /**< variables in the original SCIP which were copied while copying the component's
                                              *   constraints, but do not count to the subvars, because they were locally fixed */
   SCIP_VAR**            fixedsubvars;       /**< variables in the sub-SCIP which were copied while copying the component's
                                              *   constraints, but do not count to the subvars, because they were locally fixed */
   SCIP_Real             fixedvarsobjsum;    /**< objective contribution of all locally fixed variables */
   SCIP_Real             lastdualbound;      /**< dual bound after last optimization call for this component */
   SCIP_Real             lastprimalbound;    /**< primal bound after last optimization call for this component */
   SCIP_STATUS           laststatus;         /**< solution status of last optimization call for the sub-SCIP of this component */
   SCIP_Bool             solved;             /**< was this component solved already? */
   int                   ncalls;             /**< number of optimization calls for this component */
   int                   lastsolindex;       /**< index of best solution after last optimization call for this component */
   int                   lastbestsolindex;   /**< index of last best solution transferred to this component from the main problem */
   int                   nvars;              /**< number of variables belonging to this component */
   int                   nfixedvars;         /**< number of fixed variables copied during constraint copying */
   int                   fixedvarssize;      /**< size of fixedvars and fixedsubvars arrays */
   int                   number;             /**< component number */
} COMPONENT;

/** data related to one problem
 *  (corresponding to one node in the branch-and-bound tree and consisting of multiple components)
 */
struct Problem
{
   SCIP*                 scip;               /**< the SCIP instance this problem belongs to */
   COMPONENT*            components;         /**< independent components into which the problem can be divided */
   SCIP_PQUEUE*          compqueue;          /**< priority queue for components */
   SCIP_SOL*             bestsol;            /**< best solution found so far for the problem */
   char*                 name;               /**< name of the problem */
   SCIP_Real             fixedvarsobjsum;    /**< objective contribution of all locally fixed variables */
   SCIP_Real             lowerbound;         /**< lower bound of the problem */
   int                   ncomponents;        /**< number of independent components into which the problem can be divided */
   int                   componentssize;     /**< size of components array */
   int                   nfeascomps;         /**< number of components for which a feasible solution was found */
   int                   nsolvedcomps;       /**< number of components solved to optimality */
   int                   nlowerboundinf;     /**< number of components with lower bound equal to -infinity */
};


/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Longint          nodelimit;          /**< maximum number of nodes to be solved in subproblems */
   SCIP_Real             intfactor;          /**< the weight of an integer variable compared to binary variables */
   SCIP_Real             feastolfactor;      /**< parameter to increase the feasibility tolerance in all sub-SCIPs */
   int                   maxintvars;         /**< maximum number of integer (or binary) variables to solve a subproblem
                                              *   directly (-1: no solving) */
   int                   maxdepth;           /**< maximum depth of a node to run components detection (-1: disable
                                              *   component detection during solving) */
   int                   minsize;            /**< minimum absolute size (in terms of variables) to solve a component
                                              *   individually during branch-and-bound */
   SCIP_Real             minrelsize;         /**< minimum relative size (in terms of variables) to solve a component
                                              *   individually during branch-and-bound */
   int                   subscipdepth;       /**< depth offset of the current (sub-)problem compared to the original
                                              *   problem */
};


/** comparison method for sorting components */
static
SCIP_DECL_SORTPTRCOMP(componentSort)
{
   SCIP* scip;
   COMPONENT* comp1;
   COMPONENT* comp2;
   SCIP_Real gap1;
   SCIP_Real gap2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   comp1 = (COMPONENT*)elem1;
   comp2 = (COMPONENT*)elem2;

   if( comp1->ncalls == 0 )
      if( comp2->ncalls == 0 )
         return comp1->number - comp2->number;
      else
         return -1;
   else if( comp2->ncalls == 0 )
      return 1;

   /* the main sorting criterion is the absolute gap; however, we devide it by the number of solving calls for this
    * component to diversify the search if one component does not improve
    * @todo investigate other sorting criteria
    */
   gap1 = SQR(comp1->lastprimalbound - comp1->lastdualbound) / comp1->ncalls;
   gap2 = SQR(comp2->lastprimalbound - comp2->lastdualbound) / comp2->ncalls;

   assert(comp1->problem != NULL);
   assert(comp1->problem == comp2->problem);
   assert(comp1->problem->scip == comp2->problem->scip);

   scip = comp1->problem->scip;
   assert(scip != NULL);

   if( SCIPisFeasGT(scip, gap1, gap2) )
      return -1;
   else if( SCIPisFeasLT(scip, gap1, gap2) )
      return +1;
   else
      return comp1->number - comp2->number;
}

/** returns minimum size of components to be solved individually during the branch-and-bound search */
static
int getMinsize(
   SCIP*                 scip,               /**< main SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   int minsize;

   assert(conshdlrdata != NULL);

   minsize = (int)(conshdlrdata->minrelsize * SCIPgetNVars(scip));
   minsize = MAX(minsize, conshdlrdata->minsize);

   return minsize;
}

/** initialize component structure */
static
SCIP_RETCODE initComponent(
   PROBLEM*              problem             /**< subproblem structure */
   )
{
   COMPONENT* component;
   SCIP* scip;

   assert(problem != NULL);
   assert(problem->ncomponents < problem->componentssize);

   scip = problem->scip;
   assert(scip != NULL);

   component = &problem->components[problem->ncomponents];

   component->problem = problem;
   component->subscip = NULL;
   component->workingsol = NULL;
   component->vars = NULL;
   component->subvars = NULL;
   component->fixedvars = NULL;
   component->fixedsubvars = NULL;
   component->fixedvarsobjsum = 0.0;
   component->lastdualbound = -SCIPinfinity(scip);
   component->lastprimalbound = SCIPinfinity(scip);
   component->laststatus = SCIP_STATUS_UNKNOWN;
   component->solved = FALSE;
   component->ncalls = 0;
   component->lastsolindex = -1;
   component->lastbestsolindex = -1;
   component->nvars = 0;
   component->nfixedvars = 0;
   component->fixedvarssize = 0;
   component->number = problem->ncomponents;

   ++problem->ncomponents;

   return SCIP_OKAY;
}

/** free component structure */
static
SCIP_RETCODE freeComponent(
   COMPONENT*            component           /**< pointer to component structure */
   )
{
   PROBLEM* problem;
   SCIP* scip;

   assert(component != NULL);

   problem = component->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   SCIPdebugMsg(scip, "freeing component %d of problem <%s>\n", component->number, component->problem->name);

   assert((component->vars != NULL) == (component->subvars != NULL));
   if( component->vars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &component->vars, component->nvars);
      SCIPfreeBlockMemoryArray(scip, &component->subvars, component->nvars);
   }

   assert((component->fixedvars != NULL) == (component->fixedsubvars != NULL));
   if( component->fixedvars != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &component->fixedsubvars, component->fixedvarssize);
      SCIPfreeBlockMemoryArray(scip, &component->fixedvars, component->fixedvarssize);
   }

   /* free sub-SCIP belonging to this component and the working solution */
   if( component->subscip != NULL )
   {
      if( component->workingsol != NULL )
      {
         SCIP_CALL( SCIPfreeSol(component->subscip, &component->workingsol) );
      }

      SCIP_CALL( SCIPfree(&component->subscip) );
   }

   return SCIP_OKAY;
}


/** create the working solution for a given component, store fixed variables and the corresponding objective offset */
static
SCIP_RETCODE componentSetupWorkingSol(
   COMPONENT*            component,          /**< component structure */
   SCIP_HASHMAP*         varmap              /**< variable hashmap */
   )
{
   SCIP* subscip;

   assert(component != NULL);

   subscip = component->subscip;
   assert(subscip != NULL);
   assert(SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM);

   /* the solution should live in the primal, not the origprimal, of the sub-SCIP, so we need to transform it first */
   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcreateOrigSol(subscip, &(component->workingsol), NULL) );

   /* the number of variables was increased by copying the constraints */
   if( SCIPgetNOrigVars(subscip) > component->nvars )
   {
      PROBLEM* problem;
      SCIP* scip;
      SCIP_VAR** sourcevars;
      SCIP_VAR* subvar;
      int nsourcevars;
      int nnewvars;
      int idx = 0;
      int nvars;
      int v;

      problem = component->problem;
      assert(problem != NULL);

      scip = problem->scip;
      assert(scip != NULL);

      sourcevars = SCIPgetVars(scip);
      nsourcevars = SCIPgetNVars(scip);
      nnewvars = SCIPgetNOrigVars(subscip);
      nvars = component->nvars;

      component->fixedvarssize = nnewvars - nvars;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &component->fixedvars, component->fixedvarssize) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &component->fixedsubvars, component->fixedvarssize) );

      for( v = 0; v < nsourcevars; ++v )
      {
         subvar = (SCIP_VAR*)SCIPhashmapGetImage(varmap, sourcevars[v]);
         if( subvar != NULL && SCIPvarGetIndex(subvar) >= nvars )
         {
            /* the variable is either locally fixed or could be an inactive variable present in a constraint
             * for which an aggregation constraint linking it to the active variable was created in the subscip
             */
            assert(SCIPisZero(subscip, SCIPvarGetObj(subvar)) ||
               SCIPisEQ(subscip, SCIPvarGetLbGlobal(subvar), SCIPvarGetUbGlobal(subvar)));

            /* variable is gloablly fixed in sub-SCIP, so it was locally fixed in the main-SCIP */
            if( SCIPisEQ(subscip, SCIPvarGetLbGlobal(subvar), SCIPvarGetUbGlobal(subvar)) )
            {
               assert(SCIPisEQ(scip, SCIPvarGetLbLocal(sourcevars[v]), SCIPvarGetUbLocal(sourcevars[v])));

               component->fixedvarsobjsum += SCIPvarGetLbGlobal(subvar) * SCIPvarGetObj(subvar);
               component->fixedvars[idx] = sourcevars[v];
               component->fixedsubvars[idx] = subvar;
               ++idx;

               SCIP_CALL( SCIPsetSolVal(subscip, component->workingsol, subvar, SCIPvarGetLbGlobal(subvar)) );
            }
            /* inactive variable */
            else
            {
               assert(SCIPisZero(subscip, SCIPvarGetObj(subvar)));
            }
         }
         else
         {
            assert(subvar == NULL || SCIPisLT(scip, SCIPvarGetLbGlobal(sourcevars[v]), SCIPvarGetUbGlobal(sourcevars[v])));
            assert(subvar == NULL || SCIPisLT(subscip, SCIPvarGetLbGlobal(subvar), SCIPvarGetUbGlobal(subvar)));
         }
      }
      component->nfixedvars = idx;
      assert(component->nfixedvars <= component->fixedvarssize);
      SCIPdebugMsg(scip, "%d locally fixed variables have been copied, objective contribution: %g\n",
         component->nfixedvars, component->fixedvarsobjsum);
   }

   /* set up debug solution */
#ifdef WITH_DEBUG_SOLUTION
   {
      PROBLEM* problem;
      SCIP* scip;
      SCIP_Bool isvalid = FALSE;

      problem = component->problem;
      assert(problem != NULL);

      scip = problem->scip;
      assert(scip != NULL);

      SCIP_CALL( SCIPdebugSolIsValidInSubtree(scip, &isvalid) );

      if( isvalid )
      {
         SCIP_Real val;
         int i;

         SCIPdebugSolEnable(component->subscip);

         for( i = 0; i < component->nvars; ++i )
         {
            SCIP_CALL( SCIPdebugGetSolVal(scip, component->vars[i], &val) );
            SCIP_CALL( SCIPdebugAddSolVal(component->subscip, component->subvars[i], val) );
         }
         for( i = 0; i < component->nfixedvars; ++i )
         {
            SCIP_CALL( SCIPdebugGetSolVal(scip, component->fixedvars[i], &val) );
            SCIP_CALL( SCIPdebugAddSolVal(component->subscip, component->fixedsubvars[i], val) );
         }
      }
   }
#endif

   return SCIP_OKAY;
}

/** create a sub-SCIP for the given variables and constraints */
static
SCIP_RETCODE createSubscip(
   SCIP*                 scip,               /**< main SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP**                subscip             /**< pointer to store created sub-SCIP */
   )
{
   SCIP_Bool success;

   assert(conshdlrdata != NULL);

   /* create a new SCIP instance */
   SCIP_CALL( SCIPcreate(subscip) );

   /* copy plugins, we omit pricers (because we do not run if there are active pricers) and dialogs */
#ifdef SCIP_MORE_DEBUG /* we print statistics later, so we need to copy statistics tables */
   SCIP_CALL( SCIPcopyPlugins(scip, *subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, &success) );
#else
   SCIP_CALL( SCIPcopyPlugins(scip, *subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, &success) );
#endif

   /* the plugins were successfully copied */
   if( success )
   {
      SCIP_CONSHDLR* newconshdlr;
      SCIP_CONSHDLRDATA* newconshdlrdata;
#ifdef WITH_DEBUG_SOLUTION
      SCIP_Bool isvalid = FALSE;
#endif

      /* copy parameter settings */
      SCIP_CALL( SCIPcopyParamSettings(scip, *subscip) );

      /* some general settings should not be fixed */
      assert(!SCIPisParamFixed(*subscip, "limits/solutions"));
      assert(!SCIPisParamFixed(*subscip, "limits/bestsol"));
      assert(!SCIPisParamFixed(*subscip, "misc/usevartable"));
      assert(!SCIPisParamFixed(*subscip, "misc/useconstable"));
      assert(!SCIPisParamFixed(*subscip, "numerics/feastol"));
      assert(!SCIPisParamFixed(*subscip, "misc/usesmalltables"));

      /* disable solution limits */
      SCIP_CALL( SCIPsetIntParam(*subscip, "limits/solutions", -1) );
      SCIP_CALL( SCIPsetIntParam(*subscip, "limits/bestsol", -1) );

      /* reduce the effort spent for hash tables; however, if the debug solution is enabled and valid in this subtree,
       * hash tables are needed for installing the debug solution
       */
#ifdef WITH_DEBUG_SOLUTION
      SCIP_CALL( SCIPdebugSolIsValidInSubtree(scip, &isvalid) );
      if( !isvalid && SCIPgetStage(scip) > SCIP_STAGE_PRESOLVING )
#endif
      {
         SCIP_CALL( SCIPsetBoolParam(*subscip, "misc/usevartable", FALSE) );
         SCIP_CALL( SCIPsetBoolParam(*subscip, "misc/useconstable", FALSE) );
      }

      /* disable presolving */
      SCIP_CALL( SCIPsetPresolving(*subscip, SCIP_PARAMSETTING_OFF, TRUE) );

      /* disable component presolving and fix the parameter */
      SCIP_CALL( SCIPsetIntParam(*subscip, "constraints/" CONSHDLR_NAME "/maxprerounds", 0) );
      SCIP_CALL( SCIPfixParam(*subscip, "constraints/" CONSHDLR_NAME "/maxprerounds") );

      /* find the components constraint handler in the sub-SCIP and inform it about the actual depth in the tree */
      newconshdlr = SCIPfindConshdlr(*subscip, CONSHDLR_NAME);
      assert(newconshdlr != NULL);

      newconshdlrdata = SCIPconshdlrGetData(newconshdlr);
      assert(newconshdlrdata != NULL);
      newconshdlrdata->subscipdepth = conshdlrdata->subscipdepth + SCIPgetDepth(scip);

      /* disable output, unless in extended debug mode */
#ifndef SCIP_MORE_DEBUG
      SCIP_CALL( SCIPsetIntParam(*subscip, "display/verblevel", 0) );
#endif
   }
   else
   {
      SCIP_CALL( SCIPfree(subscip) );
      *subscip = NULL;
   }

   return SCIP_OKAY;
}

/** copies the given variables and constraints to the given sub-SCIP */
static
SCIP_RETCODE copyToSubscip(
   SCIP*                 scip,               /**< source SCIP */
   SCIP*                 subscip,            /**< target SCIP */
   const char*           name,               /**< name for copied problem */
   SCIP_VAR**            vars,               /**< array of variables to copy */
   SCIP_VAR**            subvars,            /**< array to fill with copied vars */
   SCIP_CONS**           conss,              /**< constraint to copy */
   SCIP_HASHMAP*         varmap,             /**< hashmap used for the copy process of variables */
   SCIP_HASHMAP*         consmap,            /**< hashmap used for the copy process of constraints */
   int                   nvars,              /**< number of variables to copy */
   int                   nconss,             /**< number of constraints to copy */
   SCIP_Bool*            success             /**< pointer to store whether copying was successful */
   )
{
   SCIP_CONS* newcons;
   int i;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(vars != NULL);
   assert(subvars != NULL);
   assert(conss != NULL);
   assert(varmap != NULL);
   assert(consmap != NULL);
   assert(success != NULL);

   *success = TRUE;

   /* create problem in sub-SCIP */
   SCIP_CALL( SCIPcopyProb(scip, subscip, varmap, consmap, FALSE, name) );

   /* copy variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPgetVarCopy(scip, subscip, vars[i], &subvars[i], varmap, consmap, FALSE, success) );

      /* abort if variable was not successfully copied */
      if( !(*success) )
         return SCIP_OKAY;
   }
   assert(nvars == SCIPgetNOrigVars(subscip));

   /* copy constraints */
   for( i = 0; i < nconss; ++i )
   {
      assert(!SCIPconsIsModifiable(conss[i]));

      /* copy the constraint */
      SCIP_CALL( SCIPgetConsCopy(scip, subscip, conss[i], &newcons, SCIPconsGetHdlr(conss[i]), varmap, consmap, NULL,
            SCIPconsIsInitial(conss[i]), SCIPconsIsSeparated(conss[i]), SCIPconsIsEnforced(conss[i]),
            SCIPconsIsChecked(conss[i]), SCIPconsIsPropagated(conss[i]), FALSE, FALSE,
            SCIPconsIsDynamic(conss[i]), SCIPconsIsRemovable(conss[i]), FALSE, FALSE, success) );

      /* abort if constraint was not successfully copied */
      if( !(*success) )
         return SCIP_OKAY;

      SCIP_CALL( SCIPaddCons(subscip, newcons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &newcons) );
   }

   return SCIP_OKAY;
}

/** create the sub-SCIP for a given component */
static
SCIP_RETCODE componentCreateSubscip(
   COMPONENT*            component,          /**< component structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_HASHMAP*         varmap,             /**< variable hashmap used to improve performance */
   SCIP_HASHMAP*         consmap,            /**< constraint hashmap used to improve performance */
   SCIP_CONS**           conss,              /**< constraints contained in this component */
   int                   nconss,             /**< number of constraints contained in this component */
   SCIP_Bool*            success             /**< pointer to store whether the copying process was successful */
   )
{
   char name[SCIP_MAXSTRLEN];
   PROBLEM* problem;
   SCIP* scip;
   int minsize;

   assert(component != NULL);
   assert(consmap != NULL);
   assert(conss != NULL);
   assert(success != NULL);
   assert(component->nvars > 0);

   problem = component->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   (*success) = TRUE;

   SCIP_CALL( createSubscip(scip, conshdlrdata, &component->subscip) );

   if( component->subscip != NULL )
   {
      /* get minimum size of components to solve individually and set the parameter in the sub-SCIP */
      minsize = getMinsize(scip, conshdlrdata);

      SCIP_CALL( SCIPsetIntParam(component->subscip, "constraints/" CONSHDLR_NAME "/minsize", minsize) );

      /* get name of the original problem and add "comp_nr" */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", problem->name, component->number);

      SCIP_CALL( copyToSubscip(scip, component->subscip, name, component->vars, component->subvars,
            conss, varmap, consmap, component->nvars, nconss, success) );

      if( !(*success) )
      {
         SCIP_CALL( SCIPfree(&component->subscip) );
         component->subscip = NULL;
      }
   }
   else
      (*success) = FALSE;

   return SCIP_OKAY;
}

/** solve a given sub-SCIP up to the given limits */
static
SCIP_RETCODE solveSubscip(
   SCIP*                 scip,               /**< main SCIP */
   SCIP*                 subscip,            /**< sub-SCIP to solve */
   SCIP_Longint          nodelimit,          /**< node limit */
   SCIP_Real             gaplimit            /**< gap limit */
   )
{
   SCIP_Real timelimit;
   SCIP_Real softtimelimit;
   SCIP_Real memorylimit;

   assert(scip != NULL);
   assert(subscip != NULL);

   /* set time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      timelimit += SCIPgetSolvingTime(subscip);
   }

   /* set soft time limit, if specified in main SCIP */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/softtime", &softtimelimit) );
   if( softtimelimit > -0.5 )
   {
      softtimelimit -= SCIPgetSolvingTime(scip);
      softtimelimit += SCIPgetSolvingTime(subscip);
      softtimelimit = MAX(softtimelimit, 0.0);
   }

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   /* @todo count memory of other components */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= 0.0)
   {
      SCIPdebugMessage("--> not solved (not enough memory or time left)\n");
      return SCIP_OKAY;
   }

   /* SCIP copy limits will set wrong time limits since it does not take into account time spent already in the
    * sub-SCIP; nevertheless, we call it to set the memory limit and unset all other limits, if set in the main SCIP
    */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );

   /* set time and memory limit for the subproblem */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/softtime", softtimelimit) );

   /* set gap limit */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/gap", gaplimit) );

   /* set node limit */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );

   /* solve the subproblem */
   SCIP_CALL( SCIPsolve(subscip) );

#ifdef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPprintBestSol(subscip, NULL, FALSE) );
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

   return SCIP_OKAY;
}

/** solve a connected component during presolving and evaluate the result */
static
SCIP_RETCODE solveAndEvalSubscip(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< the components constraint handler data */
   SCIP*                 subscip,            /**< sub-SCIP to be solved */
   SCIP_VAR**            vars,               /**< array of variables copied to this component */
   SCIP_VAR**            subvars,            /**< array of sub-SCIP variables corresponding to the vars array */
   SCIP_CONS**           conss,              /**< array of constraints copied to this component */
   int                   nvars,              /**< number of variables copied to this component */
   int                   nconss,             /**< number of constraints copied to this component */
   int*                  ndeletedconss,      /**< pointer to store the number of deleted constraints */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int*                  ntightenedbounds,   /**< pointer to store the number of bound tightenings */
   SCIP_RESULT*          result,             /**< pointer to store the result of the component solving */
   SCIP_Bool*            solved              /**< pointer to store if the problem was solved to optimality */
   )
{
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(subscip != NULL);
   assert(vars != NULL);
   assert(conss != NULL);
   assert(ndeletedconss != NULL);
   assert(nfixedvars != NULL);
   assert(ntightenedbounds != NULL);
   assert(result != NULL);

   *solved  = FALSE;

   SCIP_CALL( solveSubscip(scip, subscip, conshdlrdata->nodelimit, 0.0) );

   if( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
   {
      SCIP_SOL* sol;
      SCIP_VAR* var;
      SCIP_VAR* subvar;
      SCIP_Real* fixvals;
      SCIP_Bool feasible;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      sol = SCIPgetBestSol(subscip);

#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, TRUE, TRUE) );
#else
      SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, FALSE, FALSE) );
#endif

      SCIPdebugMessage("--> solved to optimality: time=%.2f, solution is%s feasible\n", SCIPgetSolvingTime(subscip), feasible ? "" : " not");

      SCIP_CALL( SCIPallocBufferArray(scip, &fixvals, nvars) );

      if( feasible )
      {
         SCIP_Real glb;
         SCIP_Real gub;

         /* get values of variables in the optimal solution */
         for( i = 0; i < nvars; ++i )
         {
            var = vars[i];
            subvar = subvars[i];

            /* get global bounds */
            glb = SCIPvarGetLbGlobal(var);
            gub = SCIPvarGetUbGlobal(var);

            if( subvar != NULL )
            {
               /* get solution value from optimal solution of the sub-SCIP */
               fixvals[i] = SCIPgetSolVal(subscip, sol, subvar);

               assert(SCIPisFeasLE(scip, fixvals[i], SCIPvarGetUbLocal(var)));
               assert(SCIPisFeasGE(scip, fixvals[i], SCIPvarGetLbLocal(var)));

               /* checking a solution is done with a relative tolerance of feasibility epsilon, if we really want to
                * change the bounds of the variables by fixing them, the old bounds must not be violated by more than
                * the absolute epsilon; therefore, we change the fixing values, if needed, and mark that the solution
                * has to be checked again
                */
               if( SCIPisGT(scip, fixvals[i], gub) )
               {
                  SCIPdebugMessage("variable <%s> fixval: %f violates global upperbound: %f\n",
                     SCIPvarGetName(var), fixvals[i], gub);
                  fixvals[i] = gub;
                  feasible = FALSE;
               }
               else if( SCIPisLT(scip, fixvals[i], glb) )
               {
                  SCIPdebugMessage("variable <%s> fixval: %f violates global lowerbound: %f\n",
                     SCIPvarGetName(var), fixvals[i], glb);
                  fixvals[i] = glb;
                  feasible = FALSE;
               }
               assert(SCIPisLE(scip, fixvals[i], SCIPvarGetUbLocal(var)));
               assert(SCIPisGE(scip, fixvals[i], SCIPvarGetLbLocal(var)));
            }
            else
            {
               /* the variable was not copied, so it was cancelled out of constraints during copying;
                * thus, the variable is not constrained and we fix it to its best bound
                */
               if( SCIPisPositive(scip, SCIPvarGetObj(var)) )
                  fixvals[i] = glb;
               else if( SCIPisNegative(scip, SCIPvarGetObj(var)) )
                  fixvals[i] = gub;
               else
               {
                  fixvals[i] = 0.0;
                  fixvals[i] = MIN(fixvals[i], gub);
                  fixvals[i] = MAX(fixvals[i], glb);
               }
            }
         }

         /* the solution value of at least one variable is feasible with a relative tolerance of feasibility epsilon,
          * but infeasible with an absolute tolerance of epsilon; try to set the variables to the bounds and check
          * solution again in the original space (changing the values might now introduce infeasibilities of constraints)
          */
         if( !feasible )
         {
            SCIP_Real origobj;

            SCIPdebugMessage("solution violates bounds by more than epsilon, check the corrected solution...\n");

            origobj = SCIPgetSolOrigObj(subscip, SCIPgetBestSol(subscip));

            SCIP_CALL( SCIPfreeTransform(subscip) );

            SCIP_CALL( SCIPcreateOrigSol(subscip, &sol, NULL) );

            /* set solution values of variables */
            for( i = 0; i < nvars; ++i )
            {
               SCIP_CALL( SCIPsetSolVal(subscip, sol, subvars[i], fixvals[i]) );
            }

            /* check the solution; integrality and bounds should be fulfilled and do not have to be checked */
            SCIP_CALL( SCIPcheckSol(subscip, sol, FALSE, FALSE, FALSE, FALSE, TRUE, &feasible) );

#ifndef NDEBUG
            /* in debug mode, we additionally check integrality and bounds */
            if( feasible )
            {
               SCIP_CALL( SCIPcheckSol(subscip, sol, FALSE, FALSE, TRUE, TRUE, FALSE, &feasible) );
               assert(feasible);
            }
#endif

            SCIPdebugMessage("--> corrected solution is%s feasible\n", feasible ? "" : " not");

            if( !SCIPisFeasEQ(subscip, SCIPsolGetOrigObj(sol), origobj) )
            {
               SCIPdebugMessage("--> corrected solution has a different objective value (old=%16.9g, corrected=%16.9g)\n",
                  origobj, SCIPsolGetOrigObj(sol));

               feasible = FALSE;
            }

            SCIP_CALL( SCIPfreeSol(subscip, &sol) );
         }

         /* if the solution is feasible, fix variables and delete constraints of the component */
         if( feasible )
         {
            /* fix variables */
            for( i = 0; i < nvars; ++i )
            {
               assert(SCIPisLE(scip, fixvals[i], SCIPvarGetUbLocal(vars[i])));
               assert(SCIPisGE(scip, fixvals[i], SCIPvarGetLbLocal(vars[i])));
               assert(SCIPisLE(scip, fixvals[i], SCIPvarGetUbGlobal(vars[i])));
               assert(SCIPisGE(scip, fixvals[i], SCIPvarGetLbGlobal(vars[i])));

               SCIP_CALL( SCIPfixVar(scip, vars[i], fixvals[i], &infeasible, &fixed) );
               SCIPvarMarkDeleteGlobalStructures(vars[i]);
               assert(!infeasible);
               assert(fixed);
               (*nfixedvars)++;
            }

            /* delete constraints */
            for( i = 0; i < nconss; ++i )
            {
               SCIP_CALL( SCIPdelCons(scip, conss[i]) );
               (*ndeletedconss)++;
            }

            *result = SCIP_SUCCESS;
            *solved = TRUE;
         }
      }

      SCIPfreeBufferArray(scip, &fixvals);
   }
   else if( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE )
   {
      *result = SCIP_CUTOFF;
   }
   else if( SCIPgetStatus(subscip) == SCIP_STATUS_UNBOUNDED || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD )
   {
      /* TODO: store unbounded ray in original SCIP data structure */
      *result = SCIP_UNBOUNDED;
   }
   else
   {
      SCIPdebugMessage("--> solving interrupted (status=%d, time=%.2f)\n",
         SCIPgetStatus(subscip), SCIPgetSolvingTime(subscip));

      /* transfer global fixings to the original problem; we can only do this, if we did not find a solution in the
       * subproblem, because otherwise, the primal bound might lead to dual reductions that cannot be transferred to
       * the original problem without also transferring the possibly suboptimal solution (which is currently not
       * possible)
       */
      if( SCIPgetNSols(subscip) == 0 )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         int ntightened;

         ntightened = 0;

         for( i = 0; i < nvars; ++i )
         {
            assert(subvars[i] != NULL);

            SCIP_CALL( SCIPtightenVarLb(scip, vars[i], SCIPvarGetLbGlobal(subvars[i]), FALSE,
                  &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened )
               ntightened++;

            SCIP_CALL( SCIPtightenVarUb(scip, vars[i], SCIPvarGetUbGlobal(subvars[i]), FALSE,
                  &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened )
               ntightened++;
         }

         *result = SCIP_SUCCESS;

         *ntightenedbounds += ntightened;

         SCIPdebugMessage("--> tightened %d bounds of variables due to global bounds in the sub-SCIP\n", ntightened);
      }
   }

   return SCIP_OKAY;
}

/** (continues) solving a connected component */
static
SCIP_RETCODE solveComponent(
   COMPONENT*            component,          /**< component structure */
   SCIP_Bool             lastcomponent,      /**< is this the last component to be solved? */
   SCIP_RESULT*          result              /**< pointer to store the result of the solving process */
   )
{
   PROBLEM* problem;
   SCIP* scip;
   SCIP* subscip;
   SCIP_SOL* bestsol;
   SCIP_Longint nodelimit;
   SCIP_Longint lastnnodes;
   SCIP_Real gaplimit;
   SCIP_STATUS status;

   assert(component != NULL);

   problem = component->problem;
   assert(problem != NULL);

   scip = problem->scip;
   assert(scip != NULL);

   subscip = component->subscip;
   assert(subscip != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("solve component <%s> (ncalls=%d, absgap=%.9g)\n",
      SCIPgetProbName(subscip), component->ncalls, component->lastprimalbound - component->lastdualbound);

   bestsol = SCIPgetBestSol(scip);

   /* update best solution of component */
   if( bestsol != NULL && SCIPsolGetIndex(bestsol) != component->lastbestsolindex )
   {
      SCIP_SOL* compsol = component->workingsol;
      SCIP_VAR** vars = component->vars;
      SCIP_VAR** subvars = component->subvars;
      int nvars = component->nvars;
      int v;

      component->lastbestsolindex = SCIPsolGetIndex(bestsol);

      /* set solution values of component variables */
      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPsetSolVal(subscip, compsol, subvars[v], SCIPgetSolVal(scip, bestsol, vars[v])) );
      }
#ifndef NDEBUG
      for( v = 0; v < component->nfixedvars; ++v )
      {
         assert(SCIPisEQ(scip, SCIPgetSolVal(subscip, compsol, component->fixedsubvars[v]),
               SCIPvarGetLbGlobal(component->fixedsubvars[v])));
      }
#endif

      if( SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM
         || SCIPisLT(subscip, SCIPgetSolOrigObj(subscip, compsol), SCIPgetPrimalbound(subscip)) )
      {
         SCIP_Bool feasible;

         SCIPdebugMessage("checking new solution in component <%s> inherited from problem <%s>: primal bound %.9g --> %.9g\n",
            SCIPgetProbName(subscip), problem->name,
            SCIPgetStage(subscip) == SCIP_STAGE_PROBLEM ? SCIPinfinity(subscip) : SCIPgetPrimalbound(subscip),
            SCIPgetSolOrigObj(subscip, compsol));

         SCIP_CALL( SCIPcheckSolOrig(subscip, compsol, &feasible, FALSE, FALSE) );
         if( feasible )
         {
            SCIPdebugMessage("... feasible, adding solution.\n");

            SCIP_CALL( SCIPaddSol(subscip, compsol, &feasible) );
         }

         /* We cannot take the value of compsol as a cutoff bound if it was not feasible; some of the fixed connecting
          * variables are different and might not allow for a better solution in this component, but still for far
          * better solutions in other components. Therefore, the only cutoffbound we can apply is the cutoffbound
          * of the problem reduced by the dual bounds of the other components
          */
         if( problem->nlowerboundinf == 0 || (problem->nlowerboundinf == 1
               && SCIPisInfinity(scip, -component->lastdualbound)) )
         {
            SCIP_Real newcutoffbound = SCIPgetSolTransObj(scip, bestsol);

            assert(problem->nlowerboundinf > 0 || SCIPisGE(scip, newcutoffbound, problem->lowerbound));

            newcutoffbound = newcutoffbound - problem->lowerbound + component->fixedvarsobjsum;

            if( problem->nlowerboundinf == 0 )
               newcutoffbound += component->lastdualbound;

            if( SCIPisSumLT(subscip, newcutoffbound, SCIPgetCutoffbound(subscip)) )
            {
               SCIPdebugMessage("update cutoff bound to %16.9g\n", newcutoffbound);

               SCIP_CALL( SCIPupdateCutoffbound(subscip, newcutoffbound) );
            }
         }
      }
   }

   assert(component->laststatus != SCIP_STATUS_OPTIMAL);

   SCIPdebugMsg(scip, "solve sub-SCIP for component <%s> (ncalls=%d, absgap=%16.9g)\n",
      SCIPgetProbName(component->subscip), component->ncalls, component->lastprimalbound - component->lastdualbound);

   if( component->ncalls == 0 )
   {
      nodelimit = 1LL;
      gaplimit = 0.0;

      lastnnodes = 0;
   }
   else
   {
      SCIP_Longint mainnodelimit;

      lastnnodes = SCIPgetNNodes(component->subscip);

      SCIP_CALL( SCIPgetLongintParam(scip, "limits/nodes", &mainnodelimit) );

      nodelimit = 2 * lastnnodes;
      nodelimit = MAX(nodelimit, 10LL);

      if( mainnodelimit != -1 )
      {
         assert(mainnodelimit >= lastnnodes);
         nodelimit = MIN(nodelimit, mainnodelimit - lastnnodes);
      }

      /* set a gap limit of half the current gap (at most 10%) */
      if( SCIPgetGap(component->subscip) < 0.2 )
         gaplimit = 0.5 * SCIPgetGap(component->subscip);
      else
         gaplimit = 0.1;

      if( lastcomponent )
         gaplimit = 0.0;
   }

   SCIP_CALL( solveSubscip(scip, subscip, nodelimit, gaplimit) );

   SCIPaddNNodes(scip, SCIPgetNNodes(subscip) - lastnnodes);

   SCIP_CALL( SCIPprintDisplayLine(scip, NULL, SCIP_VERBLEVEL_NORMAL, TRUE) );

   status = SCIPgetStatus(subscip);

   component->laststatus = status;
   ++component->ncalls;

   SCIPdebugMsg(scip, "--> (status=%d, nodes=%lld, time=%.2f): gap: %12.5g%% absgap: %16.9g\n",
      status, SCIPgetNNodes(subscip), SCIPgetSolvingTime(subscip), 100.0*SCIPgetGap(subscip),
      SCIPgetPrimalbound(subscip) - SCIPgetDualbound(subscip));

   *result = SCIP_SUCCESS;

   switch( status )
   {
   case SCIP_STATUS_OPTIMAL:
      component->solved = TRUE;
      break;
   case SCIP_STATUS_INFEASIBLE:
      component->solved = TRUE;

      /* the problem is really infeasible */
      if( SCIPisInfinity(subscip, SCIPgetPrimalbound(subscip)) )
      {
         *result = SCIP_CUTOFF;
      }
      /* the cutoff bound was reached; no solution better than the cutoff bound exists */
      else
      {
         *result = SCIP_SUCCESS;
         component->laststatus = SCIP_STATUS_OPTIMAL;
         assert(SCIPisLE(subscip, SCIPgetDualbound(subscip), SCIPgetPrimalbound(subscip)));
      }
      break;
   case SCIP_STATUS_UNBOUNDED:
   case SCIP_STATUS_INFORUNBD:
      /* TODO: store unbounded ray in original SCIP data structure */
      *result = SCIP_UNBOUNDED;
      component->solved = TRUE;
      break;
   case SCIP_STATUS_USERINTERRUPT:
      SCIP_CALL( SCIPinterruptSolve(scip) );
      break;
   case SCIP_STATUS_UNKNOWN:
   case SCIP_STATUS_NODELIMIT:
   case SCIP_STATUS_TOTALNODELIMIT:
   case SCIP_STATUS_STALLNODELIMIT:
   case SCIP_STATUS_TIMELIMIT:
   case SCIP_STATUS_MEMLIMIT:
   case SCIP_STATUS_GAPLIMIT:
   case SCIP_STATUS_SOLLIMIT:
   case SCIP_STATUS_BESTSOLLIMIT:
   case SCIP_STATUS_RESTARTLIMIT:
   default:
      break;
   }

   /* evaluate call */
   if( *result == SCIP_SUCCESS )
   {
      SCIP_SOL* sol = SCIPgetBestSol(subscip);
      SCIP_VAR* var;
      SCIP_VAR* subvar;
      SCIP_Real newdualbound;
      int v;

      /* get dual bound as the minimum of SCIP dual bound and sub-problems dual bound */
      newdualbound = SCIPgetDualbound(subscip) - component->fixedvarsobjsum;

      /* update dual bound of problem */
      if( !SCIPisEQ(scip, component->lastdualbound, newdualbound) )
      {
         assert(!SCIPisInfinity(scip, -newdualbound));

         /* first finite dual bound: decrease inf counter and add dual bound to problem dualbound */
         if( SCIPisInfinity(scip, -component->lastdualbound) )
         {
            --problem->nlowerboundinf;
            problem->lowerbound += newdualbound;
         }
         /* increase problem dual bound by dual bound delta */
         else
         {
            problem->lowerbound += (newdualbound - component->lastdualbound);
         }

         /* update problem dual bound if all problem components have a finite dual bound */
         if( problem->nlowerboundinf == 0 )
         {
            SCIPdebugMessage("component <%s>: dual bound increased from %16.9g to %16.9g, new dual bound of problem <%s>: %16.9g (gap: %16.9g, absgap: %16.9g)\n",
               SCIPgetProbName(subscip), component->lastdualbound, newdualbound, problem->name,
               SCIPretransformObj(scip, problem->lowerbound),
               problem->nfeascomps == problem->ncomponents ?
               (SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound)) /
               MAX( ABS( SCIPretransformObj(scip, problem->lowerbound) ), SCIPgetSolOrigObj(scip, problem->bestsol) ) /*lint !e666*/
               : SCIPinfinity(scip),
               problem->nfeascomps == problem->ncomponents ?
               SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound) : SCIPinfinity(scip));
            SCIP_CALL( SCIPupdateLocalLowerbound(scip, problem->lowerbound) );
         }

         /* store dual bound of this call */
         component->lastdualbound = newdualbound;
      }

      /* update primal solution of problem */
      if( sol != NULL && component->lastsolindex != SCIPsolGetIndex(sol) )
      {
         component->lastsolindex = SCIPsolGetIndex(sol);

         if( SCIPsolGetHeur(sol) != NULL )
            SCIPsolSetHeur(problem->bestsol, SCIPfindHeur(scip, SCIPheurGetName(SCIPsolGetHeur(sol))));
         else
            SCIPsolSetHeur(problem->bestsol, NULL);

         /* increase counter for feasible problems if no solution was known before */
         if( SCIPisInfinity(scip, component->lastprimalbound) )
            ++(problem->nfeascomps);

         /* update working best solution in problem */
         for( v = 0; v < component->nvars; ++v )
         {
            var = component->vars[v];
            subvar = component->subvars[v];
            assert(var != NULL);
            assert(subvar != NULL);
            assert(SCIPvarIsActive(var));

            SCIP_CALL( SCIPsetSolVal(scip, problem->bestsol, var, SCIPgetSolVal(subscip, sol, subvar)) );
         }

         /* if we have a feasible solution for each component, add the working solution to the main problem */
         if( problem->nfeascomps == problem->ncomponents )
         {
            SCIP_Bool feasible;
#ifdef SCIP_MORE_DEBUG
            SCIP_CALL( SCIPcheckSol(scip, problem->bestsol, TRUE, FALSE, TRUE, TRUE, TRUE, &feasible) );
            assert(feasible);
#endif
            SCIP_CALL( SCIPaddSol(scip, problem->bestsol, &feasible) );

            SCIPdebugMessage("component <%s>: primal bound decreased from %16.9g to %16.9g, new primal bound of problem <%s>: %16.9g (gap: %16.9g, absgap: %16.9g)\n",
               SCIPgetProbName(subscip), component->lastprimalbound, SCIPgetPrimalbound(subscip), problem->name,
               SCIPgetSolOrigObj(scip, problem->bestsol),
               problem->nfeascomps == problem->ncomponents ?
               (SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound)) /
               MAX( ABS( SCIPretransformObj(scip, problem->lowerbound) ),SCIPgetSolOrigObj(scip, problem->bestsol) ) /*lint !e666*/
               : SCIPinfinity(scip),
               problem->nfeascomps == problem->ncomponents ?
               SCIPgetSolOrigObj(scip, problem->bestsol) - SCIPretransformObj(scip, problem->lowerbound) : SCIPinfinity(scip));
         }

         /* store primal bound of this call */
         component->lastprimalbound = SCIPgetPrimalbound(subscip) - component->fixedvarsobjsum;
      }

      /* if the component was solved to optimality, we increase the respective counter and free the subscip */
      if( component->laststatus == SCIP_STATUS_OPTIMAL || component->laststatus == SCIP_STATUS_INFEASIBLE ||
            component->laststatus == SCIP_STATUS_UNBOUNDED || component->laststatus == SCIP_STATUS_INFORUNBD )
      {
         ++(problem->nsolvedcomps);
         component->solved = TRUE;

         /* free working solution and component */
         SCIP_CALL( SCIPfreeSol(subscip, &component->workingsol) );

         SCIP_CALL( SCIPfree(&subscip) );
         component->subscip = NULL;
      }
   }

   return SCIP_OKAY;
}

/** initialize subproblem structure */
static
SCIP_RETCODE initProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   PROBLEM**             problem,            /**< pointer to subproblem structure */
   SCIP_Real             fixedvarsobjsum,    /**< objective contribution of all locally fixed variables */
   int                   ncomponents         /**< number of independent components */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(problem != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBlockMemory(scip, problem) );
   assert(*problem != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*problem)->components, ncomponents) );

   /* create a priority queue for the components: we need exactly ncomponents slots in the queue so it should never be
    * resized
    */
   SCIP_CALL( SCIPpqueueCreate(&(*problem)->compqueue, ncomponents, 1.2, componentSort) );

   (*problem)->scip = scip;
   (*problem)->lowerbound = fixedvarsobjsum;
   (*problem)->fixedvarsobjsum = fixedvarsobjsum;
   (*problem)->ncomponents = 0;
   (*problem)->componentssize = ncomponents;
   (*problem)->nlowerboundinf = ncomponents;
   (*problem)->nfeascomps = 0;
   (*problem)->nsolvedcomps = 0;

   if( SCIPgetDepth(scip) == 0 )
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   else
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_node_%d", SCIPgetProbName(scip), SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*problem)->name, name, strlen(name)+1) );

   SCIP_CALL( SCIPcreateSol(scip, &(*problem)->bestsol, NULL) );

   for( v = 0; v < nvars; v++ )
   {
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v])) )
      {
         SCIP_CALL( SCIPsetSolVal(scip, (*problem)->bestsol, vars[v],
               (SCIPvarGetUbLocal(vars[v]) + SCIPvarGetLbLocal(vars[v]))/2) );
      }
   }

   SCIPdebugMessage("initialized problem <%s>\n", (*problem)->name);

   return SCIP_OKAY;
}

/** free subproblem structure */
static
SCIP_RETCODE freeProblem(
   PROBLEM**             problem             /**< pointer to problem to free */
   )
{
   SCIP* scip;
   int c;

   assert(problem != NULL);
   assert(*problem != NULL);

   scip = (*problem)->scip;
   assert(scip != NULL);

   /* free best solution */
   if( (*problem)->bestsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &(*problem)->bestsol) );
   }

   /* free all components */
   for( c = (*problem)->ncomponents - 1; c >= 0; --c )
   {
      SCIP_CALL( freeComponent(&(*problem)->components[c]) );
   }
   if( (*problem)->components != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*problem)->components, (*problem)->componentssize);
   }

   /* free priority queue */
   SCIPpqueueFree(&(*problem)->compqueue);

   /* free problem name */
   SCIPfreeMemoryArray(scip, &(*problem)->name);

   /* free PROBLEM struct and set the pointer to NULL */
   SCIPfreeBlockMemory(scip, problem);
   *problem = NULL;

   return SCIP_OKAY;
}

/** creates and captures a components constraint */
static
SCIP_RETCODE createConsComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   PROBLEM*              problem             /**< problem to be stored in the constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* find the samediff constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("components constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, (SCIP_CONSDATA*)problem,
         FALSE, FALSE, FALSE, FALSE, TRUE,
         TRUE, FALSE, FALSE, FALSE, TRUE) );

   return SCIP_OKAY;
}


/** sort the components by size and sort vars and conss arrays by component numbers */
static
SCIP_RETCODE sortComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_CONS**           conss,              /**< constraints */
   SCIP_VAR**            vars,               /**< variables */
   int*                  varcomponent,       /**< component numbers for the variables */
   int*                  conscomponent,      /**< array to store component numbers for the constraints */
   int                   nconss,             /**< number of constraints */
   int                   nvars,              /**< number of variables */
   int*                  firstvaridxpercons, /**< array with index of first variable in vars array for each constraint */
   int*                  ncompsminsize,      /**< pointer to store the number of components not exceeding the minimum size */
   int*                  ncompsmaxsize       /**< pointer to store the number of components not exceeding the maximum size */
   )
{
   SCIP_Real* compsize;
   int* permu;
   int ncomponents;
   int nbinvars;
   int nintvars;
   int ndiscvars;
   int ncontvars;
   int minsize;
   int v;
   int c;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(digraph != NULL);
   assert(conss != NULL);
   assert(vars != NULL);
   assert(firstvaridxpercons != NULL);

   /* compute minimum size of components to solve individually */
   minsize = getMinsize(scip, conshdlrdata);

   ncomponents = SCIPdigraphGetNComponents(digraph);
   *ncompsminsize = 0;
   *ncompsmaxsize = 0;

   /* We want to sort the components in increasing complexity (number of discrete variables,
    * integer weighted with factor intfactor, continuous used as tie-breaker).
    * Therefore, we now get the variables for each component, count the different variable types
    * and compute a size as described above. Then, we rename the components
    * such that for i < j, component i has no higher complexity than component j.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &compsize, ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &permu, ncomponents) );

   /* get number of variables in the components */
   for( c = 0; c < ncomponents; ++c )
   {
      int* cvars;
      int ncvars;

      SCIPdigraphGetComponent(digraph, c, &cvars, &ncvars);
      permu[c] = c;
      nbinvars = 0;
      nintvars = 0;

      for( v = 0; v < ncvars; ++v )
      {
         /* check whether variable is of binary or integer type */
         if( SCIPvarGetType(vars[cvars[v]]) == SCIP_VARTYPE_BINARY )
            nbinvars++;
         else if( SCIPvarGetType(vars[cvars[v]]) == SCIP_VARTYPE_INTEGER )
            nintvars++;
      }
      ncontvars = ncvars - nintvars - nbinvars;
      ndiscvars = (int)(nbinvars + conshdlrdata->intfactor * nintvars);
      compsize[c] = ((1000.0 * ndiscvars + (950.0 * ncontvars)/nvars));

      /* component fulfills the maxsize requirement */
      if( ndiscvars <= conshdlrdata->maxintvars )
         ++(*ncompsmaxsize);

      /* component fulfills the minsize requirement */
      if( ncvars >= minsize )
         ++(*ncompsminsize);
   }

   /* get permutation of component numbers such that the size of the components is increasing */
   SCIPsortRealInt(compsize, permu, ncomponents);

   /* now, we need the reverse direction, i.e., for each component number, we store its new number
    * such that the components are sorted; for this, we abuse the conscomponent array
    */
   for( c = 0; c < ncomponents; ++c )
      conscomponent[permu[c]] = c;

   /* for each variable, replace the old component number by the new one */
   for( c = 0; c < nvars; ++c )
      varcomponent[c] = conscomponent[varcomponent[c]];

   SCIPfreeBufferArray(scip, &permu);
   SCIPfreeBufferArray(scip, &compsize);

   /* do the mapping from calculated components per variable to corresponding
    * constraints and sort the component-arrays for faster finding the
    * actual variables and constraints belonging to one component
    */
   for( c = 0; c < nconss; c++ )
      conscomponent[c] = (firstvaridxpercons[c] == -1 ? -1 : varcomponent[firstvaridxpercons[c]]);

   SCIPsortIntPtr(varcomponent, (void**)vars, nvars);
   SCIPsortIntPtr(conscomponent, (void**)conss, nconss);

   return SCIP_OKAY;
}



/** create PROBLEM structure for the current node and split it into components */
static
SCIP_RETCODE createAndSplitProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_Real             fixedvarsobjsum,    /**< objective contribution of all locally fixed variables */
   SCIP_VAR**            sortedvars,         /**< array of unfixed variables sorted by components */
   SCIP_CONS**           sortedconss,        /**< array of (checked) constraints sorted by components */
   int*                  compstartsvars,     /**< start points of components in sortedvars array */
   int*                  compstartsconss,    /**< start points of components in sortedconss array */
   int                   ncomponents,        /**< number of components */
   PROBLEM**             problem             /**< pointer to store problem structure */
   )
{
   COMPONENT* component;
   SCIP_HASHMAP* consmap;
   SCIP_HASHMAP* varmap;
   SCIP_VAR** compvars;
   SCIP_CONS** compconss;
   SCIP_Bool success = TRUE;
   int nfixedvars = SCIPgetNVars(scip) - compstartsvars[ncomponents];
   int ncompconss;
   int comp;

   /* init subproblem data structure */
   SCIP_CALL( initProblem(scip, problem, fixedvarsobjsum, ncomponents) );
   assert((*problem)->components != NULL);

   /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
   SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), compstartsconss[ncomponents]) );

   /* loop over all components */
   for( comp = 0; comp < ncomponents; comp++ )
   {
      SCIP_CALL( initComponent(*problem) );
      assert((*problem)->ncomponents == comp+1);

      component = &(*problem)->components[comp];

      /* get component variables and store them in component structure */
      compvars = &(sortedvars[compstartsvars[comp]]);
      component->nvars = compstartsvars[comp + 1 ] - compstartsvars[comp];
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &component->vars, compvars, component->nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &component->subvars, component->nvars) );
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), component->nvars + nfixedvars) );

      /* get component constraints */
      compconss = &(sortedconss[compstartsconss[comp]]);
      ncompconss = compstartsconss[comp + 1] - compstartsconss[comp];

#ifdef DETAILED_OUTPUT
      /* print details about the component including its size */
      if( component->nvars > 1 && ncompconss > 1 )
      {
         int nbinvars = 0;
         int nintvars = 0;
         int ncontvars = 0;
         int i;

         for( i = 0; i < component->nvars; ++i )
         {
            if( SCIPvarGetType(compvars[i]) == SCIP_VARTYPE_BINARY )
               ++nbinvars;
            else if( SCIPvarGetType(compvars[i]) == SCIP_VARTYPE_INTEGER )
               ++nintvars;
            else
               ++ncontvars;
         }
         SCIPdebugMsg(scip, "component %d at node %lld, depth %d (%d): %d vars (%d bin, %d int, %d cont), %d conss\n",
            comp, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), SCIPgetDepth(scip), SCIPgetDepth(scip) + conshdlrdata->subscipdepth,
            component->nvars, nbinvars, nintvars, ncontvars, ncompconss);
      }
#endif
      assert(ncompconss > 0 || component->nvars == 1);

      SCIPdebugMsg(scip, "build sub-SCIP for component %d of problem <%s>: %d vars, %d conss\n",
         component->number, (*problem)->name, component->nvars, ncompconss);

#ifndef NDEBUG
      {
         int i;
         for( i = 0; i < component->nvars; ++i )
            assert(SCIPvarIsActive(component->vars[i]));
      }
#endif

      /* build subscip for component */
      SCIP_CALL( componentCreateSubscip(component, conshdlrdata, varmap, consmap, compconss, ncompconss, &success) );

      if( success )
      {
         SCIP_CALL( componentSetupWorkingSol(component, varmap) );

         /* add component to the priority queue of the problem structure */
         SCIP_CALL( SCIPpqueueInsert((*problem)->compqueue, component) );
      }

      SCIPhashmapFree(&varmap);

      if( !success )
         break;
   }

   SCIPhashmapFree(&consmap);

   if( !success )
   {
      /* free subproblem data structure since not all component could be copied */
      SCIP_CALL( freeProblem(problem) );
   }

   return SCIP_OKAY;
}

/** continue solving a problem  */
static
SCIP_RETCODE solveProblem(
   PROBLEM*              problem,            /**< problem structure */
   SCIP_RESULT*          result              /**< result pointer for the problem solve */
   )
{
   COMPONENT* component;
   SCIP_RESULT subscipresult;

   assert(problem != NULL);

   *result = SCIP_SUCCESS;

   component = (COMPONENT*)SCIPpqueueRemove(problem->compqueue);

   /* continue solving the component */
   SCIP_CALL( solveComponent(component, SCIPpqueueNElems(problem->compqueue) == 0, &subscipresult) );

   /* if infeasibility or unboundedness was detected, return this */
   if( subscipresult == SCIP_CUTOFF || subscipresult == SCIP_UNBOUNDED )
   {
      *result = subscipresult;
   }
   /* the component was not solved to optimality, so we need to re-insert it in the components queue */
   else if( !component->solved )
   {
      SCIP_CALL( SCIPpqueueInsert(problem->compqueue, component) );
      *result = SCIP_DELAYNODE;
   }
   /* no unsolved components are left, so this problem has be completely evaluated and the node can be pruned */
   else if( SCIPpqueueNElems(problem->compqueue) == 0 )
      *result = SCIP_CUTOFF;
   /* there are some unsolved components left, so we delay this node */
   else
      *result = SCIP_DELAYNODE;

   return SCIP_OKAY;
}

/*
 * Local methods
 */

/** loop over constraints, get active variables and fill directed graph */
static
SCIP_RETCODE fillDigraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  unfixedvarpos,      /**< mapping from variable problem index to unfixed var index */
   int                   nunfixedvars,       /**< number of unfixed variables */
   int*                  firstvaridxpercons, /**< array to store for each constraint the index in the local vars array
                                              *   of the first variable of the constraint */
   SCIP_Bool*            success             /**< flag indicating successful directed graph filling */
   )
{
   SCIP_VAR** consvars;
   int requiredsize;
   int nconsvars;
   int nvars;
   int idx1;
   int idx2;
   int c;
   int v;

   assert(scip != NULL);
   assert(digraph != NULL);
   assert(conss != NULL);
   assert(firstvaridxpercons != NULL);
   assert(success != NULL);

   *success = TRUE;

   nconsvars = 0;
   requiredsize = 0;
   nvars = SCIPgetNVars(scip);

   /* allocate buffer for storing active variables per constraint; size = nvars ensures that it will be big enough */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );

   for( c = 0; c < nconss; ++c )
   {
      /* check for reached timelimit */
      if( (c % 1000 == 0) && SCIPisStopped(scip) )
      {
         *success = FALSE;
         break;
      }

      /* get number of variables for this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, success) );

      if( !(*success) )
         break;

      /* reallocate consvars array, if needed */
      if( nconsvars > nvars )
      {
         nvars = nconsvars;
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, nvars) );
      }

#ifndef NDEBUG
      /* clearing variables array to check for consistency */
      if( nconsvars == nvars )
      {
	 BMSclearMemoryArray(consvars, nconsvars);
      }
      else
      {
	 assert(nconsvars < nvars);
	 BMSclearMemoryArray(consvars, nconsvars + 1);
      }
#endif

      /* get variables for this constraint */
      SCIP_CALL( SCIPgetConsVars(scip, conss[c], consvars, nvars, success) );

      if( !(*success) )
      {
#ifndef NDEBUG
	 /* it looks strange if returning the number of variables was successful but not returning the variables */
	 SCIPwarningMessage(scip, "constraint <%s> returned number of variables but returning variables failed\n", SCIPconsGetName(conss[c]));
#endif
         break;
      }

#ifndef NDEBUG
      /* check if returned variables are consistent with the number of variables that were returned */
      for( v = nconsvars - 1; v >= 0; --v )
	 assert(consvars[v] != NULL);
      if( nconsvars < nvars )
	 assert(consvars[nconsvars] == NULL);
#endif

      /* transform given variables to active variables */
      SCIP_CALL( SCIPgetActiveVars(scip, consvars, &nconsvars, nvars, &requiredsize) );
      assert(requiredsize <= nvars);

      firstvaridxpercons[c] = -1;

      /* store the index of the first unfixed variable and add edges to the directed graph */
      if( nconsvars > 0 )
      {
         v = 0;
         idx1 = -1;

         /* go through variables until the first unfixed one is reached (which has unfixedvarpos >= 0) */
         while( idx1 == -1 && v < nconsvars )
         {
            idx1 = SCIPvarGetProbindex(consvars[v]);
            assert(idx1 >= 0);
            idx1 = unfixedvarpos[idx1];
            assert(idx1 < nunfixedvars);
            ++v;
         }

         if( idx1 >= 0 )
         {
            /* save index of the first variable for later component assignment */
            firstvaridxpercons[c] = idx1;

            /* create sparse directed graph; sparse means to add only those edges necessary for component calculation,
             * i.e., add edges from the first variable to all others
             */
            for(; v < nconsvars; ++v )
            {
               idx2 = SCIPvarGetProbindex(consvars[v]);
               assert(idx2 >= 0);
               idx2 = unfixedvarpos[idx2];
               assert(idx2 < nunfixedvars);

               /* variable is fixed */
               if( idx2 < 0 )
                  continue;

               /* we add only one directed edge, because the other direction is automatically added for component computation */
               SCIP_CALL( SCIPdigraphAddArc(digraph, idx1, idx2, NULL) );
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

/** search for components in the problem */
static
SCIP_RETCODE findComponents(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< the components constraint handler data */
   SCIP_Real*            fixedvarsobjsum,    /**< objective contribution of all locally fixed variables, or NULL if
                                              *   fixed variables should not be disregarded */
   SCIP_VAR**            sortedvars,         /**< array to store variables sorted by components, should have enough size
                                              *   for all variables */
   SCIP_CONS**           sortedconss,        /**< array to store (checked) constraints sorted by components, should have
                                              *   enough size for all constraints */
   int*                  compstartsvars,     /**< start points of components in sortedvars array */
   int*                  compstartsconss,    /**< start points of components in sortedconss array */
   int*                  nsortedvars,        /**< pointer to store the number of variables belonging to any component */
   int*                  nsortedconss,       /**< pointer to store the number of (checked) constraints in components */
   int*                  ncomponents,        /**< pointer to store the number of components */
   int*                  ncompsminsize,      /**< pointer to store the number of components not exceeding the minimum size */
   int*                  ncompsmaxsize       /**< pointer to store the number of components not exceeding the maximum size */

   )
{
   SCIP_CONS** tmpconss;
   SCIP_VAR** vars;
   SCIP_Bool success;
   int ntmpconss;
   int nvars;
   int c;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(sortedvars != NULL);
   assert(sortedconss != NULL);
   assert(compstartsvars != NULL);
   assert(compstartsconss != NULL);
   assert(nsortedvars != NULL);
   assert(nsortedconss != NULL);
   assert(ncomponents != NULL);
   assert(ncompsminsize != NULL);
   assert(ncompsmaxsize != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   if( fixedvarsobjsum != NULL )
      *fixedvarsobjsum = 0.0;

   *ncomponents = 0;
   *ncompsminsize = 0;
   *ncompsmaxsize = 0;

   /* collect checked constraints for component detection */
   ntmpconss = SCIPgetNConss(scip);
   tmpconss = SCIPgetConss(scip);
   (*nsortedconss) = 0;
   for( c = 0; c < ntmpconss; c++ )
   {
      sortedconss[(*nsortedconss)] = tmpconss[c];
      (*nsortedconss)++;
   }

   if( nvars > 1 && *nsortedconss > 1 )
   {
      int* unfixedvarpos;
      int* firstvaridxpercons;
      int* varlocks;
      int nunfixedvars = 0;
      int v;

      /* arrays for storing the first variable in each constraint (for later component assignment), the number of
       * variable locks, and the positions in the sortedvars array for all unfixed variables
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &firstvaridxpercons, *nsortedconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varlocks, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &unfixedvarpos, nvars) );

      /* count number of varlocks for each variable (up + down locks) and multiply it by 2;
       * that value is used as an estimate of the number of arcs incident to the variable's node in the digraph
       * to be safe, we double this value
       */
      for( v = 0; v < nvars; ++v )
      {
         /* variable is not fixed or we do not want to disregard fixed variables */
         if( (fixedvarsobjsum == NULL) || SCIPisLT(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetUbLocal(vars[v])) )
         {
            assert(nunfixedvars <= v);
            sortedvars[nunfixedvars] = vars[v];
            varlocks[nunfixedvars] = 4 * (SCIPvarGetNLocksDown(vars[v]) + SCIPvarGetNLocksUp(vars[v]));
            unfixedvarpos[v] = nunfixedvars;
            ++nunfixedvars;
         }
         /* variable is fixed; update the objective sum of all fixed variables */
         else
         {
            unfixedvarpos[v] = -1;
            (*fixedvarsobjsum) += SCIPvarGetObj(vars[v]) * SCIPvarGetLbLocal(vars[v]);
         }
      }
      *nsortedvars = nunfixedvars;

      if( nunfixedvars > 0 )
      {
         SCIP_DIGRAPH* digraph;

         /* create and fill directed graph */
         SCIP_CALL( SCIPcreateDigraph(scip, &digraph, nunfixedvars) );
         SCIP_CALL( SCIPdigraphSetSizes(digraph, varlocks) );
         SCIP_CALL( fillDigraph(scip, digraph, sortedconss, *nsortedconss, unfixedvarpos, nunfixedvars, firstvaridxpercons, &success) );

         if( success )
         {
            int* varcomponent;
            int* conscomponent;

            SCIP_CALL( SCIPallocBufferArray(scip, &varcomponent, nunfixedvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &conscomponent, MAX(nunfixedvars,*nsortedconss)) );

            /* compute independent components */
            SCIP_CALL( SCIPdigraphComputeUndirectedComponents(digraph, 1, varcomponent, ncomponents) );

            if( *ncomponents > 1 )
            {
               int nconss = *nsortedconss;
               int i;

               nvars = *nsortedvars;

               SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
                  "cons components found %d undirected components at node %lld, depth %d (%d)\n",
                  *ncomponents, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), SCIPgetDepth(scip), SCIPgetDepth(scip) + conshdlrdata->subscipdepth);

               /* sort components by size and sort variables and constraints by component number */
               SCIP_CALL( sortComponents(scip, conshdlrdata, digraph, sortedconss, sortedvars, varcomponent, conscomponent, nconss, *nsortedvars,
                     firstvaridxpercons, ncompsminsize, ncompsmaxsize) );

               /* determine start indices of components in sortedvars and sortedconss array */
               i = 0;

               while( i < nconss && conscomponent[i] == -1 )
                  ++i;

               for( c = 0; c < *ncomponents + 1; ++c )
               {
                  assert(i == nconss || conscomponent[i] >= c);

                  compstartsconss[c] = i;

                  while( i < nconss && conscomponent[i] == c )
                     ++i;
               }

               for( c = 0, i = 0; c < *ncomponents + 1; ++c )
               {
                  assert(i == nvars || varcomponent[i] >= c);

                  compstartsvars[c] = i;

                  while( i < nvars && varcomponent[i] == c )
                     ++i;
               }

#ifndef NDEBUG
               for( c = 0; c < *ncomponents; ++c )
               {
                  for( i = compstartsconss[c]; i < compstartsconss[c+1]; ++i )
                     assert(conscomponent[i] == c);
                  for( i = compstartsvars[c]; i < compstartsvars[c+1]; ++i )
                     assert(varcomponent[i] == c);
               }
#endif
            }

            SCIPfreeBufferArray(scip, &conscomponent);
            SCIPfreeBufferArray(scip, &varcomponent);
         }

         SCIPdigraphFree(&digraph);
      }

      SCIPfreeBufferArray(scip, &unfixedvarpos);
      SCIPfreeBufferArray(scip, &varlocks);
      SCIPfreeBufferArray(scip, &firstvaridxpercons);
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyComponents)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrComponents(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(conshdlrFreeComponents)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropComponents)
{  /*lint --e{715}*/
   PROBLEM* problem;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Longint nodelimit;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPconshdlrGetNActiveConss(conshdlr) >= 0);
   assert(SCIPconshdlrGetNActiveConss(conshdlr) <= 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* do not try to detect independent components if the depth is too high */
   if( SCIPgetDepth(scip) + conshdlrdata->subscipdepth > conshdlrdata->maxdepth
      && SCIPconshdlrGetNActiveConss(conshdlr) == 0 )
      return SCIP_OKAY;

   /* don't run in probing or in repropagation */
   if( SCIPinProbing(scip) || SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   /* do not run, if not all variables are explicitly known */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* we do not want to run, if there are no variables left */
   if( SCIPgetNVars(scip) == 0 )
      return SCIP_OKAY;

   /* check for a reached timelimit */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* the components constraint handler does kind of dual reductions */
   if( !SCIPallowDualReds(scip) || !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   problem = NULL;
   *result = SCIP_DIDNOTFIND;

   /* the current node already has a components constraint storing a problem split into individual components */
   if( SCIPconshdlrGetNActiveConss(conshdlr) >= 1 )
   {
      assert(SCIPconshdlrGetNActiveConss(conshdlr) == 1);

      problem = (PROBLEM*)SCIPconsGetData(SCIPconshdlrGetConss(conshdlr)[0]);
   }
   /* no components constraint at the current node, search for components */
   else
   {
      SCIP_Real fixedvarsobjsum;
      SCIP_VAR** sortedvars;
      SCIP_CONS** sortedconss;
      int* compstartsvars;
      int* compstartsconss;
      int nsortedvars;
      int nsortedconss;
      int ncomponents;
      int ncompsminsize;
      int ncompsmaxsize;

      assert(SCIPconshdlrGetNActiveConss(conshdlr) == 0);

      /* allocate memory for sorted components */
      SCIP_CALL( SCIPallocBufferArray(scip, &sortedvars, SCIPgetNVars(scip)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sortedconss, SCIPgetNConss(scip)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &compstartsvars, SCIPgetNVars(scip) + 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &compstartsconss, SCIPgetNVars(scip) + 1) );

      /* search for components */
      SCIP_CALL( findComponents(scip, conshdlrdata, &fixedvarsobjsum, sortedvars, sortedconss, compstartsvars,
            compstartsconss, &nsortedvars, &nsortedconss, &ncomponents, &ncompsminsize, &ncompsmaxsize) );

      if( ncompsminsize > 1 )
      {
         SCIP_CONS* cons;

         SCIPdebugMsg(scip, "found %d components (%d fulfulling the minsize requirement) at node %lld at depth %d (%d)\n",
            ncomponents, ncompsminsize, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), SCIPgetDepth(scip),
            SCIPgetDepth(scip) + conshdlrdata->subscipdepth);

         /* if there are components with size smaller than the limit, we merge them with the smallest component */
         if( ncomponents > ncompsminsize )
         {
            int minsize;
            int size;
            int c;
            int m = 0;

            /* compute minimum size of components to solve individually */
            minsize = getMinsize(scip, conshdlrdata);

            for( c = 0; c < ncomponents; ++c )
            {
               size = compstartsvars[c+1] - compstartsvars[c];

               if( size >= minsize )
               {
                  ++m;
                  compstartsvars[m] = compstartsvars[c+1];
                  compstartsconss[m] = compstartsconss[c+1];
               }
               /* the last component is too small */
               else if( c == ncomponents - 1 )
               {
                  assert(m == ncompsminsize);
                  compstartsvars[m] = compstartsvars[c+1];
                  compstartsconss[m] = compstartsconss[c+1];
               }
            }
            assert(m == ncompsminsize);
            assert(compstartsvars[m] == nsortedvars);
            assert(compstartsconss[m] == nsortedconss);

            ncomponents = m;
         }

         SCIP_CALL( createAndSplitProblem(scip, conshdlrdata, fixedvarsobjsum, sortedvars, sortedconss, compstartsvars,
               compstartsconss, ncomponents, &problem) );

         /* if the problem is not NULL, copying worked fine */
         if( problem != NULL )
         {
            SCIP_CALL( createConsComponents(scip, &cons, problem->name, problem) );
            SCIP_CALL( SCIPaddConsNode(scip, SCIPgetCurrentNode(scip), cons, NULL) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }

      SCIPfreeBufferArray(scip, &compstartsconss);
      SCIPfreeBufferArray(scip, &compstartsvars);
      SCIPfreeBufferArray(scip, &sortedconss);
      SCIPfreeBufferArray(scip, &sortedvars);
   }

   /* (continue to) solve the problem
    *
    * If the problem was not solved to optimality yet, the result code is set to SCIP_DELAYNODE, so that after the
    * propagation is finished, the node is put back into the queue of open nodes and solving the components of the
    * problem will be continued when the node is focused and propagated the next time.
    * However, if we are at the root node, we continue solving the problem until it is solved or some limit is reached
    * since there are no other nodes to process and we want to avoid calling other propagation methods or heuristics
    * again and again
    */
   SCIP_CALL( SCIPgetLongintParam(scip, "limits/nodes", &nodelimit) );
   if( nodelimit == -1 )
      nodelimit = SCIP_LONGINT_MAX;

   do
   {
      if( problem != NULL )
      {
         SCIP_CALL( solveProblem(problem, result) );
      }
   } while( *result == SCIP_DELAYNODE && SCIPgetDepth(scip) == 0 && !SCIPisStopped(scip) && SCIPgetNNodes(scip) < nodelimit);

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolComponents)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** sortedvars;
   SCIP_CONS** sortedconss;
   int* compstartsvars;
   int* compstartsconss;
   int nsortedvars;
   int nsortedconss;
   int ncomponents;
   int ncompsminsize;
   int ncompsmaxsize;
   int nvars;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPconshdlrGetNActiveConss(conshdlr) >= 0);
   assert(SCIPconshdlrGetNActiveConss(conshdlr) <= 1);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* do not run, if not all variables are explicitly known */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   nvars = SCIPgetNVars(scip);

   /* we do not want to run, if there are no variables left */
   if( nvars == 0 )
      return SCIP_OKAY;

   /* only call the components presolving, if presolving would be stopped otherwise */
   if( !SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   /* the components constraint handler does kind of dual reductions */
   if( !SCIPallowDualReds(scip) || !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   /* check for a reached timelimit */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   assert(SCIPconshdlrGetNActiveConss(conshdlr) == 0);

   /* allocate memory for sorted components */
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedvars, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedconss, SCIPgetNConss(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &compstartsvars, SCIPgetNVars(scip) + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &compstartsconss, SCIPgetNVars(scip) + 1) );

   /* search for components */
   SCIP_CALL( findComponents(scip, conshdlrdata, NULL, sortedvars, sortedconss, compstartsvars,
         compstartsconss, &nsortedvars, &nsortedconss, &ncomponents, &ncompsminsize, &ncompsmaxsize) );

   if( ncompsmaxsize > 0 )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP* subscip;
      SCIP_HASHMAP* consmap;
      SCIP_HASHMAP* varmap;
      SCIP_VAR** compvars;
      SCIP_VAR** subvars;
      SCIP_CONS** compconss;
      SCIP_Bool success;
      SCIP_Bool solved;
      int nsolved = 0;
      int ncompvars;
      int ncompconss;
      int comp;

      SCIPdebugMsg(scip, "found %d components (%d with small size) during presolving; overall problem size: %d vars (%d int, %d bin, %d cont), %d conss\n",
         ncomponents, ncompsmaxsize, SCIPgetNVars(scip), SCIPgetNBinVars(scip), SCIPgetNIntVars(scip), SCIPgetNContVars(scip) + SCIPgetNImplVars(scip), SCIPgetNConss(scip));

      /* build subscip */
      SCIP_CALL( createSubscip(scip, conshdlrdata, &subscip) );

      if( subscip == NULL )
         goto TERMINATE;

      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/usesmalltables", TRUE) );
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/" CONSHDLR_NAME "/propfreq", -1) );

      /* hashmap mapping from original constraints to constraints in the sub-SCIPs (for performance reasons) */
      SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), nsortedconss) );

      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nsortedvars) );

      /* loop over all components */
      for( comp = 0; comp < ncompsmaxsize && !SCIPisStopped(scip); comp++ )
      {
#ifdef WITH_DEBUG_SOLUTION
         if( SCIPgetStage(subscip) > SCIP_STAGE_INIT )
         {
            SCIP_CALL( SCIPfree(&subscip) );
            SCIP_CALL( createSubscip(scip, conshdlrdata, &subscip) );
         }
#endif
         /* get component variables */
         compvars = &(sortedvars[compstartsvars[comp]]);
         ncompvars = compstartsvars[comp + 1 ] - compstartsvars[comp];

         /* get component constraints */
         compconss = &(sortedconss[compstartsconss[comp]]);
         ncompconss = compstartsconss[comp + 1] - compstartsconss[comp];

         /* if we have an unlocked variable, let duality fixing do the job! */
         if( ncompconss == 0 )
         {
            assert(ncompvars == 1);
            continue;
         }

         SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), ncompvars) );
#ifdef DETAILED_OUTPUT
         {
            int nbinvars = 0;
            int nintvars = 0;
            int ncontvars = 0;
            int i;

            for( i = 0; i < ncompvars; ++i )
            {
               if( SCIPvarGetType(compvars[i]) == SCIP_VARTYPE_BINARY )
                  ++nbinvars;
               else if( SCIPvarGetType(compvars[i]) == SCIP_VARTYPE_INTEGER )
                  ++nintvars;
               else
                  ++ncontvars;
            }
            SCIPdebugMsg(scip, "solve component %d: %d vars (%d bin, %d int, %d cont), %d conss\n",
               comp, ncompvars, nbinvars, nintvars, ncontvars, ncompconss);
         }
#endif
#ifndef NDEBUG
         {
            int i;
            for( i = 0; i < ncompvars; ++i )
               assert(SCIPvarIsActive(compvars[i]));
         }
#endif

         /* get name of the original problem and add "comp_nr" */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp_%d", SCIPgetProbName(scip), comp);

         SCIP_CALL( copyToSubscip(scip, subscip, name, compvars, subvars,
               compconss, varmap, consmap, ncompvars, ncompconss, &success) );

         if( !success )
         {
            SCIPhashmapFree(&varmap);
            continue;
         }

            /* set up debug solution */
#ifdef WITH_DEBUG_SOLUTION
         {
            SCIP_SOL* debugsol;
            SCIP_Real val;
            int i;

            SCIP_CALL( SCIPdebugGetSol(scip, &debugsol) );

            /* set solution values in the debug solution if it is available */
            if( debugsol != NULL )
            {
               SCIPdebugSolEnable(subscip);

               for( i = 0; i < ncompvars; ++i )
               {
                  SCIP_CALL( SCIPdebugGetSolVal(scip, compvars[i], &val) );
                  SCIP_CALL( SCIPdebugAddSolVal(subscip, subvars[i], val) );
               }
            }
         }
#endif

         /* solve the subproblem and evaluate the result, i.e. apply fixings of variables and remove constraints */
         SCIP_CALL( solveAndEvalSubscip(scip, conshdlrdata, subscip, compvars, subvars, compconss,
               ncompvars, ncompconss, ndelconss, nfixedvars, nchgbds, result, &solved) );

         /* free variable hash map */
         SCIPhashmapFree(&varmap);

         if( solved )
            ++nsolved;

         /* if the component is unbounded or infeasible, this holds for the complete problem as well */
         if( *result == SCIP_UNBOUNDED || *result == SCIP_CUTOFF )
            break;
         /* if there is only one component left, let's solve this in the main SCIP */
         else if( nsolved == ncomponents - 1 )
            break;
      }

      SCIPfreeBufferArray(scip, &subvars);
      SCIPhashmapFree(&consmap);

      SCIP_CALL( SCIPfree(&subscip) );
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &compstartsconss);
   SCIPfreeBufferArray(scip, &compstartsvars);
   SCIPfreeBufferArray(scip, &sortedconss);
   SCIPfreeBufferArray(scip, &sortedvars);

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteComponents)
{  /*lint --e{715}*/
   PROBLEM* problem;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   problem = (PROBLEM*)(*consdata);

   SCIP_CALL( freeProblem(&problem) );

   *consdata = NULL;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxComponents)
{  /*lint --e{715}*/
   assert(result != NULL );

   /* no enforcement is performed, but the callback is needed for all constraint handlers with needscons = FALSE */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockComponents)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

#ifndef NDEBUG
/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolComponents)
{  /*lint --e{715}*/
   assert(nconss == 0);

   return SCIP_OKAY;
}
#endif

#define consEnfolpComponents NULL
#define consEnfopsComponents NULL
#define consCheckComponents NULL

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the components constraint handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrComponents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create components propagator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->subscipdepth = 0;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC, CONSHDLR_ENFOPRIORITY,
         CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpComponents, consEnfopsComponents, consCheckComponents, consLockComponents,
         conshdlrdata) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropComponents,
         CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING));
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolComponents,
         CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING));

   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, conshdlrFreeComponents) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxComponents) );
#ifndef NDEBUG
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolComponents) );
#endif
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyComponents, NULL) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteComponents) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxdepth",
         "maximum depth of a node to run components detection (-1: disable component detection during solving)",
         &conshdlrdata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxintvars",
         "maximum number of integer (or binary) variables to solve a subproblem during presolving (-1: unlimited)",
         &conshdlrdata->maxintvars, TRUE, DEFAULT_MAXINTVARS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/minsize",
         "minimum absolute size (in terms of variables) to solve a component individually during branch-and-bound",
         &conshdlrdata->minsize, TRUE, DEFAULT_MINSIZE, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/minrelsize",
         "minimum relative size (in terms of variables) to solve a component individually during branch-and-bound",
         &conshdlrdata->minrelsize, TRUE, DEFAULT_MINRELSIZE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "constraints/" CONSHDLR_NAME "/nodelimit",
         "maximum number of nodes to be solved in subproblems during presolving",
         &conshdlrdata->nodelimit, FALSE, DEFAULT_NODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/intfactor",
         "the weight of an integer variable compared to binary variables",
         &conshdlrdata->intfactor, FALSE, DEFAULT_INTFACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/feastolfactor",
         "factor to increase the feasibility tolerance of the main SCIP in all sub-SCIPs, default value 1.0",
         &conshdlrdata->feastolfactor, TRUE, DEFAULT_FEASTOLFACTOR, 0.0, 1000000.0, NULL, NULL) );


   return SCIP_OKAY;
}
