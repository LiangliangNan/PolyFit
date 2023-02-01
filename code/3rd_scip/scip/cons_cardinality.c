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

/**@file   cons_cardinality.c
 * @brief  constraint handler for cardinality constraints
 * @author Tobias Fischer
 *
 * This constraint handler handles cardinality constraints of the form
 * \f[
 *   |\mbox{supp}(x)| \leq b
 * \f]
 * with integer right-hand side \f$b\f$. Here, \f$|\mbox{supp}(x)|\f$ denotes the number of nonzero entries of the
 * vector \f$x\f$.
 *
 * Note that cardinality constraints generalize special ordered set of type one (SOS1) constraints in which \f$b = 1\f$.
 *
 * The implementation of this constraint handler is based on@n
 * "On the Structure of Linear Programs with Overlapping Cardinality Constraints"@n
 * T. Fischer and M. E. Pfetsch, Tech. rep., 2016
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_cardinality.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include <string.h>
#include <ctype.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "cardinality"
#define CONSHDLR_DESC          "cardinality constraint handler"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            10 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in
                                         *   (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_FAST

/* branching rules */
#define DEFAULT_BRANCHBALANCED    FALSE /**< whether to use balanced instead of unbalanced branching */
#define DEFAULT_BALANCEDDEPTH        20 /**< maximum depth for using balanced branching (-1: no limit) */
#define DEFAULT_BALANCEDCUTOFF      2.0 /**< determines that balanced branching is only used if the branching cut off value
                                         *   w.r.t. the current LP solution is greater than a given value */

/* event handler properties */
#define EVENTHDLR_NAME         "cardinality"
#define EVENTHDLR_DESC         "bound change event handler for cardinality constraints"


/** constraint data for cardinality constraints */
struct SCIP_ConsData
{
   SCIP_CONS*            cons;               /**< cardinality constraint */
   int                   cardval;            /**< number of variables that the constraint allows to be nonzero */
   int                   nvars;              /**< number of variables in the constraint */
   int                   maxvars;            /**< maximal number of variables (= size of storage) */
   int                   ntreatnonzeros;     /**< number of variables in constraint that are either known to be nonzero
                                              *   (because zero is not in variable domain) or may be treated as nonzero */
   SCIP_EVENTDATA**      eventdatascurrent;  /**< event datas for current bound change events */
   SCIP_VAR**            eventvarscurrent;   /**< event variables for current bound change events */
   int                   neventdatascurrent; /**< number of current bound change events */
   SCIP_EVENTDATA**      eventdatas;         /**< event data array for bound change events */
   SCIP_VAR**            vars;               /**< variables in constraint */
   SCIP_VAR**            indvars;            /**< indicator variables that indicate which variables may be treated as
                                              *   nonzero in cardinality constraint */
   SCIP_Real*            weights;            /**< weights determining the order (ascending), or NULL if not used */
   SCIP_ROW*             rowlb;              /**< row corresponding to lower bounds, or NULL if not yet created */
   SCIP_ROW*             rowub;              /**< row corresponding to upper bounds, or NULL if not yet created */
};

/** cardinality constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_HASHMAP*         varhash;            /**< hash map from implied variable to (binary) indicator variable */
   SCIP_Bool             branchbalanced;     /**< whether to use balanced instead of unbalanced branching */
   int                   balanceddepth;      /**< maximum depth for using balanced branching (-1: no limit) */
   SCIP_Real             balancedcutoff;     /**< determines that balanced branching is only used if the branching cut off
                                              *   value w.r.t. the current LP solution is greater than a given value */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
};

/** event data for bound changes events */
struct SCIP_EventData
{
   SCIP_CONSDATA*        consdata;           /**< cardinality constraint data to process the bound change for */
   SCIP_VAR*             var;                /**< implied variable */
   SCIP_VAR*             indvar;             /**< indicator variable */
   unsigned int          pos:30;             /**< position in constraint */
   unsigned int          varmarked:1;        /**< whether implied variable is marked for propagation */
   unsigned int          indvarmarked:1;     /**< whether indicator variable is marked for propagation */
};

/** catches bound change events for a variable and its indicator variable */
static
SCIP_RETCODE catchVarEventCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for bound change events */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR*             var,                /**< implied variable */
   SCIP_VAR*             indvar,             /**< indicator variable */
   int                   pos,                /**< position in constraint */
   SCIP_EVENTDATA**      eventdata           /**< pointer to store event data for bound change events */
   )
{
   assert(eventhdlr != NULL);
   assert(consdata != NULL);
   assert(var != NULL);
   assert(indvar != NULL);
   assert(pos >= 0);

   /* create event data of indicator variable */
   SCIP_CALL( SCIPallocBlockMemory(scip, eventdata) );

   (*eventdata)->consdata = consdata;
   (*eventdata)->var = var;
   (*eventdata)->indvar = indvar;
   (*eventdata)->varmarked = FALSE;
   (*eventdata)->indvarmarked = FALSE;
   (*eventdata)->pos = (unsigned int)pos;

   /* catch bound change events of each variable */
   SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, *eventdata, NULL) );
   SCIP_CALL( SCIPcatchVarEvent(scip, indvar, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, *eventdata, NULL) );

   return SCIP_OKAY;
}

/** drops bound change events for a variable and its indicator variable */
static
SCIP_RETCODE dropVarEventCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for bound change events */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR*             var,                /**< implied variable */
   SCIP_VAR*             indvar,             /**< indicator variable */
   SCIP_EVENTDATA**      eventdata           /**< pointer of event data for bound change events */
   )
{
   assert(eventhdlr != NULL);
   assert(consdata != NULL);
   assert(var != NULL);
   assert(indvar != NULL);
   assert(eventdata != NULL);

   /* drop bound change events of each variable */
   SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, *eventdata, -1) );
   SCIP_CALL( SCIPdropVarEvent(scip, indvar, SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr, *eventdata, -1) );

   /* free event data of indicator variable */
   SCIPfreeBlockMemory(scip, eventdata);
   *eventdata = NULL;

   return SCIP_OKAY;
}

/** fix variable in given node to 0 or add constraint if variable is multi-aggregated
 *
 *  @todo Try to handle multi-aggregated variables as in \ref fixVariableZero() below.
 */
static
SCIP_RETCODE fixVariableZeroNode(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0 */
   SCIP_NODE*            node,               /**< node */
   SCIP_Bool*            infeasible          /**< pointer to store if fixing is infeasible */
   )
{
   /* if variable cannot be nonzero */
   *infeasible = FALSE;
   if( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* if variable is multi-aggregated */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_CONS* cons;
      SCIP_Real val;

      val = 1.0;

      if( !SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || !SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIPdebugMsg(scip, "creating constraint to force multi-aggregated variable <%s> to 0.\n", SCIPvarGetName(var));

         /* we have to insert a local constraint var = 0 */
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "branch", 1, &var, &val, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
               TRUE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddConsNode(scip, node, cons, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   else
   {
      if ( ! SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) )
      {
         SCIP_CALL( SCIPchgVarLbNode(scip, node, var, 0.0) );
      }
      if ( ! SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIP_CALL( SCIPchgVarUbNode(scip, node, var, 0.0) );
      }
   }

   return SCIP_OKAY;
}

/** try to fix variable to 0
 *
 *  Try to treat fixing by special consideration of multiaggregated variables. For a multi-aggregation
 *  \f[
 *  x = \sum_{i=1}^n \alpha_i x_i + c,
 *  \f]
 *  we can express the fixing \f$x = 0\f$ by fixing all \f$x_i\f$ to 0 if \f$c = 0\f$ and the lower bounds of \f$x_i\f$
 *  are nonnegative if \f$\alpha_i > 0\f$ or the upper bounds are nonpositive if \f$\alpha_i < 0\f$.
 */
static
SCIP_RETCODE fixVariableZero(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable to be fixed to 0*/
   SCIP_Bool*            infeasible,         /**< if fixing is infeasible */
   SCIP_Bool*            tightened           /**< if fixing was performed */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(tightened != NULL);

   *infeasible = FALSE;
   *tightened = FALSE;

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Real aggrconst;

      /* if constant is 0 */
      aggrconst = SCIPvarGetMultaggrConstant(var);
      if( SCIPisZero(scip, aggrconst) )
      {
         SCIP_VAR** aggrvars;
         SCIP_Real* aggrvals;
         SCIP_Bool allnonnegative = TRUE;
         int naggrvars;
         int i;

         SCIP_CALL( SCIPflattenVarAggregationGraph(scip, var) );

         /* check whether all variables are "nonnegative" */
         naggrvars = SCIPvarGetMultaggrNVars(var);
         aggrvars = SCIPvarGetMultaggrVars(var);
         aggrvals = SCIPvarGetMultaggrScalars(var);
         for( i = 0; i < naggrvars; ++i )
         {
            if( (SCIPisPositive(scip, aggrvals[i]) && SCIPisNegative(scip, SCIPvarGetLbLocal(aggrvars[i]))) ||
                 (SCIPisNegative(scip, aggrvals[i]) && SCIPisPositive(scip, SCIPvarGetUbLocal(aggrvars[i]))) )
            {
               allnonnegative = FALSE;
               break;
            }
         }

         if( allnonnegative )
         {
            /* all variables are nonnegative -> fix variables */
            for( i = 0; i < naggrvars; ++i )
            {
               SCIP_Bool fixed;
               SCIP_CALL( SCIPfixVar(scip, aggrvars[i], 0.0, infeasible, &fixed) );
               if( *infeasible )
                  return SCIP_OKAY;
               *tightened = *tightened || fixed;
            }
         }
      }
   }
   else
   {
      SCIP_CALL( SCIPfixVar(scip, var, 0.0, infeasible, tightened) );
   }

   return SCIP_OKAY;
}

/** add lock on variable */
static
SCIP_RETCODE lockVariableCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             indvar              /**< indicator variable */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   /* rounding down == bad if lb < 0, rounding up == bad if ub > 0 */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)),
         SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var))) );
   SCIP_CALL( SCIPlockVarCons(scip, indvar, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/* remove lock on variable */
static
SCIP_RETCODE unlockVariableCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             indvar              /**< indicator variable */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   /* rounding down == bad if lb < 0, rounding up == bad if ub > 0 */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)),
         SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var))) );
   SCIP_CALL( SCIPunlockVarCons(scip, indvar, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** ensures that the vars and weights array can store at least num entries */
static
SCIP_RETCODE consdataEnsurevarsSizeCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   num,                /**< minimum number of entries to store */
   SCIP_Bool             reserveweights      /**< whether the weights array is handled */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->maxvars);

   if( num > consdata->maxvars )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->maxvars, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->indvars, consdata->maxvars, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->maxvars, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventdatascurrent, 4*consdata->maxvars, 4*newsize) );/*lint !e647*/
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventvarscurrent, 4*consdata->maxvars, 4*newsize) );/*lint !e647*/

      if ( reserveweights )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->weights, consdata->maxvars, newsize) );
      }
      consdata->maxvars = newsize;
   }
   assert(num <= consdata->maxvars);

   return SCIP_OKAY;
}

/** handle new variable that was added to the constraint
 *
 *  We perform the following steps:
 *
 *  - catch bound change events of variable.
 *  - update rounding locks of variable.
 *  - don't allow multiaggregation of variable, since this cannot be handled by branching in the current implementation
 *  - update lower and upper bound row, i.e., the linear representations of the cardinality constraints
 */
static
SCIP_RETCODE handleNewVariableCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             indvar,             /**< indicator variable to indicate whether variable may be treated as
                                              *   nonzero in cardinality constraint */
   int                   pos,                /**< position in constraint */
   SCIP_Bool             transformed,        /**< whether original variable was transformed */
   SCIP_EVENTDATA**      eventdata           /**< pointer to store event data for bound change events */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(conshdlrdata != NULL);
   assert(var != NULL);

   /* if we are in transformed problem, catch the variable's events */
   if( transformed )
   {
      /* catch bound change events of variable */
      SCIP_CALL( catchVarEventCardinality(scip, conshdlrdata->eventhdlr, consdata, var, indvar, pos, eventdata) );
      assert(eventdata != NULL );

      /* if the variable is fixed to nonzero */
      assert(consdata->ntreatnonzeros >= 0 );
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(indvar), 1.0) )
         ++consdata->ntreatnonzeros;
   }

   /* branching on multiaggregated variables does not seem to work well, so avoid it */
   SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, indvar) );

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockVariableCardinality(scip, cons, var, indvar) );

   /* add the new coefficient to the upper bound LP row, if necessary */
   if( consdata->rowub != NULL && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var))
       && !SCIPisZero(scip, SCIPvarGetUbGlobal(var)) )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rowub, var, 1.0/SCIPvarGetUbGlobal(var)) );
   }

   /* add the new coefficient to the lower bound LP row, if necessary */
   if( consdata->rowlb != NULL && !SCIPisInfinity(scip, SCIPvarGetLbGlobal(var))
       && !SCIPisZero(scip, SCIPvarGetLbGlobal(var)) )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rowlb, var, 1.0/SCIPvarGetLbGlobal(var)) );
   }

   return SCIP_OKAY;
}

/** adds a variable to a cardinality constraint, at position given by weight - ascending order */
static
SCIP_RETCODE addVarCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_VAR*             indvar,             /**< indicator variable to indicate whether variable may be treated as nonzero
                                              *   in cardinality constraint (or NULL) */
   SCIP_Real             weight              /**< weight to determine position */
   )
{
   SCIP_EVENTDATA* eventdata = NULL;
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;
   int pos;

   assert(var != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL );

   if( consdata->weights == NULL && consdata->maxvars > 0 )
   {
      SCIPerrorMessage("cannot add variable to cardinality constraint <%s> that does not contain weights.\n",
         SCIPconsGetName(cons));
      return SCIP_INVALIDCALL;
   }

   /* check indicator variable */
   if( indvar == NULL )
   {
      if( conshdlrdata->varhash == NULL )
      {
         /* set up hash map */
         SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->varhash, SCIPblkmem(scip), SCIPgetNTotalVars(scip)) );
      }

      /* check whether an indicator variable already exists for implied variable */
      if( SCIPhashmapExists(conshdlrdata->varhash, var) )
      {
         assert((SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, var) != NULL);
         indvar = (SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, var);
         assert(indvar != NULL);
      }
      else
      {
         /* if implied variable is binary, then it is also not necessary to create an indicator variable */
         if( SCIPvarIsBinary(var) )
            indvar = var;
         else
         {
            char varname[SCIP_MAXSTRLEN];
            SCIP_VAR* newvar;

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ind_%s", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, FALSE, FALSE,
                  NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, newvar) );
            indvar = newvar;

            SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
         }
         assert(indvar != NULL);

         /* insert implied variable to hash map */
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varhash, var, (void*) indvar) );/*lint !e571*/
         assert(indvar == (SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, var));
         assert(SCIPhashmapExists(conshdlrdata->varhash, var));
      }
   }

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
      SCIP_CALL( SCIPgetTransformedVar(scip, indvar, &indvar) );
   }
   assert(var != NULL);
   assert(indvar != NULL);
   assert(transformed == SCIPvarIsTransformed(var));
   assert(transformed == SCIPvarIsTransformed(indvar));

   /* ensure that the new information can be storend in the constraint data */
   SCIP_CALL( consdataEnsurevarsSizeCardinality(scip, consdata, consdata->nvars + 1, TRUE) );
   assert(consdata->weights != NULL);
   assert(consdata->maxvars >= consdata->nvars+1);

   /* move other variables, if necessary */
   for( pos = consdata->nvars; pos >= 1; --pos )
   {
      if( consdata->weights[pos-1] > weight )
      {
         consdata->vars[pos] = consdata->vars[pos-1];
         consdata->indvars[pos] = consdata->indvars[pos-1];
         consdata->eventdatas[pos] = consdata->eventdatas[pos-1];
         consdata->weights[pos] = consdata->weights[pos-1];

         if( consdata->eventdatas[pos] != NULL )
         {
            consdata->eventdatas[pos]->pos = (unsigned int)pos;
         }
      }
      else
         break;
   }
   assert(0 <= pos && pos <= consdata->nvars);

   /* handle the new variable */
   SCIP_CALL( handleNewVariableCardinality(scip, cons, consdata, conshdlrdata, var, indvar, pos, transformed, &eventdata) );
   assert(! transformed || eventdata != NULL);

   /* insert variable */
   consdata->vars[pos] = var;
   consdata->indvars[pos] = indvar;
   consdata->eventdatas[pos] = eventdata;
   consdata->weights[pos] = weight;
   ++consdata->nvars;

   return SCIP_OKAY;
}

/** appends a variable to a cardinality constraint */
static
SCIP_RETCODE appendVarCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_VAR*             indvar              /**< indicator variable to indicate whether variable may be treated as nonzero
                                              *   in cardinality constraint */
   )
{
   SCIP_EVENTDATA* eventdata = NULL;
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check indicator variable */
   if( indvar == NULL )
   {
      if( conshdlrdata->varhash == NULL )
      {
         /* set up hash map */
         SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->varhash, SCIPblkmem(scip), SCIPgetNTotalVars(scip)) );
      }

      /* check whether an indicator variable already exists for implied variable */
      if( SCIPhashmapExists(conshdlrdata->varhash, var) )
      {
         assert((SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, var) != NULL);
         indvar = (SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, var);
         assert(indvar != NULL);
      }
      else
      {
         /* if implied variable is binary, then it is also not necessary to create an indicator variable */
         if( SCIPvarIsBinary(var) )
            indvar = var;
         else
         {
            char varname[SCIP_MAXSTRLEN];
            SCIP_VAR* newvar;

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ind_%s", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateVar(scip, &newvar, varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, FALSE, FALSE,
                  NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, newvar) );
            indvar = newvar;

            SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
         }
         assert(indvar != NULL);

         /* insert implied variable to hash map */
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varhash, var, (void*) indvar) );/*lint !e571*/
         assert(indvar == (SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, var));
         assert(SCIPhashmapExists(conshdlrdata->varhash, var));
      }
   }

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
      SCIP_CALL( SCIPgetTransformedVar(scip, indvar, &indvar) );
   }
   assert(var != NULL);
   assert(indvar != NULL);
   assert(transformed == SCIPvarIsTransformed(var));
   assert(transformed == SCIPvarIsTransformed(indvar));

   /* ensure that the new information can be stored in the constraint data */
   SCIP_CALL( consdataEnsurevarsSizeCardinality(scip, consdata, consdata->nvars + 1, FALSE) );

   /* handle the new variable */
   SCIP_CALL( handleNewVariableCardinality(scip, cons, consdata, conshdlrdata, var, indvar, consdata->nvars, transformed,
        &eventdata) );
   assert(!transformed || eventdata != NULL);

   /* insert variable */
   consdata->vars[consdata->nvars] = var;
   consdata->indvars[consdata->nvars] = indvar;
   consdata->eventdatas[consdata->nvars] = eventdata;

   if( consdata->weights != NULL && consdata->nvars > 0 )
      consdata->weights[consdata->nvars] = consdata->weights[consdata->nvars-1] + 1.0;
   ++consdata->nvars;

   assert(consdata->weights != NULL || consdata->nvars > 0);

   return SCIP_OKAY;
}

/** deletes a variable from a cardinality constraint */
static
SCIP_RETCODE deleteVarCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< corresponding event handler */
   int                   pos                 /**< position of variable in array */
   )
{ /*lint --e{679}*/
   int j;

   assert(0 <= pos && pos < consdata->nvars);

   /* remove lock of variable */
   SCIP_CALL( unlockVariableCardinality(scip, cons, consdata->vars[pos], consdata->indvars[pos]) );

   /* drop events on indicator variable and implied variable */
   SCIP_CALL( dropVarEventCardinality(scip, eventhdlr, consdata, consdata->vars[pos], consdata->indvars[pos],
        &consdata->eventdatas[pos]) );

   /* update number of variables that may be treated as nonzero */
   if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->indvars[pos]), 1.0) )
      --(consdata->ntreatnonzeros);

   /* delete variable - need to copy since order is important */
   for( j = pos; j < consdata->nvars-1; ++j )
   {
      consdata->vars[j] = consdata->vars[j+1];
      consdata->indvars[j] = consdata->indvars[j+1];
      consdata->eventdatas[j] = consdata->eventdatas[j+1];
      if( consdata->weights != NULL )
         consdata->weights[j] = consdata->weights[j+1];

      consdata->eventdatas[j]->pos = (unsigned int)j;
   }
   --consdata->nvars;

   return SCIP_OKAY;
}

/** for each indicator variable sets solution to 1.0 if the solution value of the implied variable is nonzero */
static
SCIP_RETCODE polishPrimalSolution(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to be enforced (or NULL) */
   SCIP_SOL*             primsol             /**< primal solution */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** indvars;
   SCIP_VAR** vars;
   int nvars;
   int c;

   /* check each constraint */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons;
      int j;

      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      nvars = consdata->nvars;
      vars = consdata->vars;
      indvars = consdata->indvars;

      for( j = 0; j < nvars; ++j )
      {
         if( SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, vars[j])) )
         {
            SCIP_CALL( SCIPsetSolVal(scip, primsol, indvars[j], 0.0) );
         }
         else
         {
            SCIP_CALL( SCIPsetSolVal(scip, primsol, indvars[j], 1.0) );
         }
      }
   }

   return SCIP_OKAY;
}

/** unmark variables that are marked for propagation */
static
SCIP_RETCODE consdataUnmarkEventdataVars(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP_EVENTDATA** eventdatas;
   int nvars;
   int j;

   eventdatas = consdata->eventdatas;
   nvars = consdata->nvars;
   assert(eventdatas != NULL);

   for( j = 0; j < nvars; ++j )
   {
      SCIP_EVENTDATA* eventdata;

      eventdata = eventdatas[j];
      eventdata->varmarked = FALSE;
      eventdata->indvarmarked = FALSE;
   }

   return SCIP_OKAY;
}

/** perform one presolving round
 *
 *  We perform the following presolving steps:
 *
 *  - If the bounds of some variable force it to be nonzero, we can
 *    fix all other variables to zero and remove the cardinality constraints
 *    that contain it.
 *  - If a variable is fixed to zero, we can remove the variable.
 *  - If a variable appears twice, it can be fixed to 0.
 *  - We substitute appregated variables.
 */
static
SCIP_RETCODE presolRoundCardinality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   SCIP_Bool*            success,            /**< whether we performed a successful reduction */
   int*                  ndelconss,          /**< number of deleted constraints */
   int*                  nupgdconss,         /**< number of upgraded constraints */
   int*                  nfixedvars,         /**< number of fixed variables */
   int*                  nremovedvars        /**< number of variables removed */
   )
{
   SCIP_VAR** indvars;
   SCIP_VAR** vars;
   SCIP_Bool allvarsbinary;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(eventhdlr != NULL);
   assert(cutoff != NULL);
   assert(success != NULL);
   assert(ndelconss != NULL);
   assert(nfixedvars != NULL);
   assert(nremovedvars != NULL);

   *cutoff = FALSE;
   *success = FALSE;

   SCIPdebugMsg(scip, "Presolving cardinality constraint <%s>.\n", SCIPconsGetName(cons) );

   /* reset number of events stored for propagation, since presolving already performs a
    * complete propagation of all variables */
   consdata->neventdatascurrent = 0;
   SCIP_CALL( consdataUnmarkEventdataVars(consdata) );

   j = 0;
   allvarsbinary = TRUE;
   vars = consdata->vars;
   indvars = consdata->indvars;

   /* check for variables fixed to 0 and bounds that fix a variable to be nonzero */
   while ( j < consdata->nvars )
   {
      int l;
      SCIP_VAR* var;
      SCIP_VAR* oldvar;
      SCIP_VAR* indvar;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real indlb;
      SCIP_Real indub;
      SCIP_Real scalar;
      SCIP_Real constant;

      scalar = 1.0;
      constant = 0.0;

      /* check for aggregation: if the constant is zero the variable is zero iff the aggregated
       * variable is 0 */
      var = vars[j];
      indvar = indvars[j];
      oldvar = var;
      SCIP_CALL( SCIPgetProbvarSum(scip, &var, &scalar, &constant) );

      /* if constant is zero and we get a different variable, substitute variable */
      if( SCIPisZero(scip, constant) && !SCIPisZero(scip, scalar) && var != vars[j] )
      {
         SCIPdebugMsg(scip, "substituted variable <%s> by <%s>.\n", SCIPvarGetName(vars[j]), SCIPvarGetName(var));

         /* we reuse the same indicator variable for the new variable */
         SCIP_CALL( dropVarEventCardinality(scip, eventhdlr, consdata, consdata->vars[j], consdata->indvars[j],
              &consdata->eventdatas[j]) );
         SCIP_CALL( catchVarEventCardinality(scip, eventhdlr, consdata, var, consdata->indvars[j], j,
              &consdata->eventdatas[j]) );
         assert(consdata->eventdatas[j] != NULL);

         /* change the rounding locks */
         SCIP_CALL( unlockVariableCardinality(scip, cons, consdata->vars[j], consdata->indvars[j]) );
         SCIP_CALL( lockVariableCardinality(scip, cons, var, consdata->indvars[j]) );

         /* update event data */
         consdata->eventdatas[j]->var = var;

         vars[j] = var;
      }
      assert(var == vars[j]);

      /* check whether the variable appears again later */
      for( l = j+1; l < consdata->nvars; ++l )
      {
         if( var == vars[l] || oldvar == vars[l] )
         {
            SCIPdebugMsg(scip, "variable <%s> appears twice in constraint <%s>.\n", SCIPvarGetName(vars[j]),
                 SCIPconsGetName(cons));
            return SCIP_INVALIDDATA;
         }
      }

      /* get bounds of variable */
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      /* if the variable is fixed to nonzero */
      if( SCIPisFeasPositive(scip, lb) || SCIPisFeasNegative(scip, ub) )
      {
         assert(SCIPvarIsBinary(indvar));

         /* fix (binary) indicator variable to 1.0 (the cardinality constraint will then be modified below) */
         SCIP_CALL( SCIPfixVar(scip, indvar, 1.0, &infeasible, &fixed) );
         if( infeasible )
         {
            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         if( fixed )
         {
            SCIPdebugMsg(scip, "fixed binary variable <%s> to 1.0.\n", SCIPvarGetName(indvar));
            ++(*nfixedvars);
         }
      }

      /* if the variable is fixed to 0 */
      if( SCIPisFeasZero(scip, lb) && SCIPisFeasZero(scip, ub) )
      {
         assert(SCIPvarIsBinary(indvar));

         /* fix (binary) indicator variable to 0.0, if possible (the cardinality constraint will then be modified below)
          * note that an infeasibility implies no cut off */
         SCIP_CALL( SCIPfixVar(scip, indvar, 0.0, &infeasible, &fixed) );
         if( fixed )
         {
            SCIPdebugMsg(scip, "fixed binary variable <%s> to 0.0.\n", SCIPvarGetName(indvar));
            ++(*nfixedvars);
         }
      }

      /* get bounds of indicator variable */
      indlb = SCIPvarGetLbLocal(indvar);
      indub = SCIPvarGetUbLocal(indvar);

      /* if the variable may be treated as nonzero */
      if( SCIPisFeasEQ(scip, indlb, 1.0) )
      {
         assert(indub == 1.0);

         /* modify row and delete variable */
         SCIP_CALL( deleteVarCardinality(scip, cons, consdata, eventhdlr, j) );
         SCIPdebugMsg(scip, "deleting variable <%s> from constraint <%s>, since it may be treated as nonzero.\n",
            SCIPvarGetName(var), SCIPconsGetName(cons));
         --(consdata->cardval);
         ++(*nremovedvars);
      }
      /* if the indicator variable is fixed to 0 */
      else if( SCIPisFeasEQ(scip, indub, 0.0) )
      {
         assert(indlb == 0.0);

         /* fix variable to 0.0 */
         SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );
         if( infeasible )
         {
            *cutoff = TRUE;
            return SCIP_OKAY;
         }
         if( fixed )
         {
            SCIPdebugMsg(scip, "fixed variable <%s> to 0.0.\n", SCIPvarGetName(var));
            ++(*nfixedvars);
         }

         /* delete variable */
         SCIP_CALL( deleteVarCardinality(scip, cons, consdata, eventhdlr, j) );
         SCIPdebugMsg(scip, "deleting variable <%s> from constraint <%s>, since it is fixed to 0.\n", SCIPvarGetName(var),
              SCIPconsGetName(cons));
         ++(*nremovedvars);
      }
      else
      {
         /* check whether all variables are binary */
         if( !SCIPvarIsBinary(var) )
            allvarsbinary = FALSE;

         ++j;
      }
   }

   /* if the cardinality value is smaller than 0, then the problem is infeasible */
   if( consdata->cardval < 0 )
   {
      SCIPdebugMsg(scip, "The problem is infeasible: more variables have bounds that keep them from being 0 than allowed.\n");

      *cutoff = TRUE;
      return SCIP_OKAY;
   }
   /* else if the cardinality value is 0 */
   else if( consdata->cardval == 0 )
   {
      /* fix all variables of the constraint to 0 */
      for( j = 0; j < consdata->nvars; ++j )
      {
         SCIP_CALL( SCIPfixVar(scip, consdata->vars[j], 0.0, &infeasible, &fixed) );
         if( infeasible )
         {
            *cutoff = TRUE;
            return SCIP_OKAY;
         }
         if( fixed )
         {
            SCIPdebugMsg(scip, "fixed variable <%s> to 0.0.\n", SCIPvarGetName(consdata->vars[j]));
            ++(*nfixedvars);
         }
      }
   }

   /* if the cardinality constraint is redundant */
   if( consdata->nvars <= consdata->cardval )
   {
      SCIPdebugMsg(scip, "Deleting cardinality constraint <%s> with <%d> variables and cardinality value <%d>.\n",
         SCIPconsGetName(cons), consdata->nvars, consdata->cardval);

      /* delete constraint */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);
      *success = TRUE;
      return SCIP_OKAY;
   }
   else
   {
      /* if all variables are binary create a knapsack constraint */
      if( allvarsbinary )
      {
         SCIP_CONS* knapsackcons;
         SCIP_Longint* vals;

         SCIP_CALL( SCIPallocBufferArray(scip, &vals, consdata->nvars) );
         for( j = 0; j < consdata->nvars; ++j )
            vals[j] = 1;

         /* create, add, and release the knapsack constraint */
         SCIP_CALL( SCIPcreateConsKnapsack(scip, &knapsackcons, SCIPconsGetName(cons), consdata->nvars, consdata->vars,
              vals, (SCIP_Longint) consdata->cardval, SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
              SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
              SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
              SCIPconsIsStickingAtNode(cons)) );/*lint !e524*/
         SCIP_CALL( SCIPaddCons(scip, knapsackcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &knapsackcons) );

         SCIPfreeBufferArray(scip, &vals);

         SCIPdebugMsg(scip, "Upgrading cardinality constraint <%s> to knapsack constraint.\n", SCIPconsGetName(cons));

         /* remove the cardinality constraint globally */
         assert(!SCIPconsIsModifiable(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );
         ++(*nupgdconss);
         *success = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** propagates a cardinality constraint and its variables
 *
 *  The number 'ntreatnonzeros' that is assigned to the constraint data returns the number of variables that are either
 *  known to be nonzero or can be treated as nonzero. We say that a variable is known to be nonzero, if zero is outside
 *  the domain of this variable. A variable can be treated as nonzero, if its corresponding indicator variable 'indvar' is
 *  fixed to 1.0, e.g., by branching.
 *
 *  We perform the following propagation steps:
 *
 *  - if the number 'ntreatnonzeros' is greater than the cardinality value of the constraint, then the current subproblem
 *    is marked as infeasible.
 *  - if the cardinality constraint is saturated, i.e., the number 'ntreatnonzeros' is equal to the cardinality value of
 *    the constraint, then fix all the other variables of the constraint to zero.
 *  - remove the cardinality constraint locally if all variables are either fixed to zero or can be treated as nonzero.
 *  - if a (binary) indicator variable is fixed to zero, then fix the corresponding implied variable to zero.
 *  - if zero is outside of the domain of an implied variable, then fix the corresponding indicator variable to one.
 */
static
SCIP_RETCODE propCardinality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool*            cutoff,             /**< whether a cutoff happened */
   int*                  nchgdomain          /**< number of domain changes */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(cutoff != NULL);
   assert(nchgdomain != NULL);

   *cutoff = FALSE;

   /* if more variables may be treated as nonzero than allowed */
   if( consdata->ntreatnonzeros > consdata->cardval )
   {
      SCIPdebugMsg(scip, "the node is infeasible, more than %d variables are fixed to be nonzero.\n", consdata->cardval);
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      *cutoff = TRUE;

      return SCIP_OKAY;
   }

   /* if number of nonzeros is saturated */
   if( consdata->ntreatnonzeros == consdata->cardval )
   {
      SCIP_VAR** vars;
      SCIP_VAR** indvars;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;
      SCIP_Bool allvarfixed;
      int nvars;
      int cnt = 0;
      int j;

      nvars = consdata->nvars;
      vars = consdata->vars;
      indvars = consdata->indvars;
      assert(vars != NULL);
      assert(indvars != NULL);

      /* fix free variables to zero */
      allvarfixed = TRUE;
      for( j = 0; j < nvars; ++j )
      {
         /* if variable is implied to be treated as nonzero */
         if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(indvars[j]), 1.0) )
            ++cnt;
         /* else fix variable to zero if not done already */
         else
         {
            SCIP_VAR* var;

            var = vars[j];

            /* fix variable */
            if( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var)) )
            {
               SCIP_CALL( fixVariableZero(scip, var, &infeasible, &tightened) );
               if( infeasible )
               {
                  SCIPdebugMsg(scip, "the node is infeasible, more than %d variables are fixed to be nonzero.\n",
                     consdata->cardval);
                  SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  *cutoff = TRUE;

                  return SCIP_OKAY;
               }

               if( tightened )
               {
                  SCIPdebugMsg(scip, "fixed variable <%s> to 0, since constraint <%s> with cardinality value %d is \
                     saturated.\n", SCIPvarGetName(var), SCIPconsGetName(cons), consdata->cardval);
                  ++(*nchgdomain);
               }
               else
                  allvarfixed = FALSE;
            }
         }
      }
      assert(cnt == consdata->ntreatnonzeros);

      /* reset constraint age counter */
      if( *nchgdomain > 0 )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }

      /* delete constraint locally */
      if( allvarfixed )
      {
         assert(!SCIPconsIsModifiable(cons));
         SCIP_CALL( SCIPdelConsLocal(scip, cons) );

         return SCIP_OKAY;
      }
   }

   /* if relevant bound change events happened */
   if( consdata->neventdatascurrent > 0 )
   {
      SCIP_EVENTDATA** eventdatas;
      SCIP_VAR** eventvars;
      int neventdatas;
      int j;

      neventdatas = consdata->neventdatascurrent;
      eventvars = consdata->eventvarscurrent;
      eventdatas = consdata->eventdatascurrent;
      assert(eventdatas != NULL && eventvars != NULL);

      for( j = 0; j < neventdatas; ++j )
      {
         SCIP_EVENTDATA* eventdata;
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         SCIP_VAR* var;

         eventdata = eventdatas[j];
         var = eventvars[j];
         assert(var != NULL && eventdata != NULL);
         assert(eventdata->var != NULL);
         assert(eventdata->indvar != NULL);
         assert(var == eventdata->var || var == eventdata->indvar);
         assert(SCIPvarIsBinary(eventdata->indvar));

         /* if variable is an indicator variable */
         if( eventdata->indvar == var )
         {
            assert(eventdata->indvarmarked);

            /* if variable is fixed to zero */
            if( SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
            {
               SCIP_VAR* implvar;

               implvar = eventdata->var;

               /* fix implied variable to zero if not done already */
               if( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(implvar)) ||
                     SCIPisFeasPositive(scip, SCIPvarGetUbLocal(implvar)) )
               {
                  SCIP_CALL( fixVariableZero(scip, implvar, &infeasible, &tightened) );

                  if( infeasible )
                  {
                     SCIPdebugMsg(scip, "the node is infeasible, indicator variable %s is fixed to zero although implied "
                        "variable %s is nonzero.\n", SCIPvarGetName(var), SCIPvarGetName(implvar));
                     SCIP_CALL( SCIPresetConsAge(scip, cons) );
                     *cutoff = TRUE;

                     return SCIP_OKAY;
                  }

                  if( tightened )
                  {
                     SCIPdebugMsg(scip, "fixed variable <%s> to 0, since indicator variable %s is 0.\n",
                        SCIPvarGetName(implvar), SCIPvarGetName(var));
                     ++(*nchgdomain);
                  }
               }
            }
            eventdata->indvarmarked = FALSE;
         }
         /* else if variable is an implied variable */
         else
         {
            assert(eventdata->var == var);
            assert(eventdata->varmarked);

            /* if variable is is nonzero */
            if( SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var)) || SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
            {
               SCIP_VAR* indvar;

               indvar = eventdata->indvar;
               assert(SCIPvarIsBinary(indvar));

               /* fix indicator variable to 1.0 if not done already */
               if( !SCIPisFeasEQ(scip, SCIPvarGetLbLocal(indvar), 1.0) )
               {
                  /* if fixing is infeasible */
                  if( SCIPvarGetUbLocal(indvar) != 1.0 )
                  {
                     SCIPdebugMsg(scip, "the node is infeasible, implied variable %s is fixed to nonzero "
                        "although indicator variable %s is 0.\n", SCIPvarGetName(var), SCIPvarGetName(indvar));
                     SCIP_CALL( SCIPresetConsAge(scip, cons) );
                     *cutoff = TRUE;

                     return SCIP_OKAY;
                  }
                  SCIP_CALL( SCIPchgVarLb(scip, indvar, 1.0) );
                  SCIPdebugMsg(scip, "fixed variable <%s> to 1.0, since implied variable %s is nonzero.\n",
                     SCIPvarGetName(indvar), SCIPvarGetName(var));
                  ++(*nchgdomain);
               }
            }
            eventdata->varmarked = FALSE;
         }
      }
   }
   consdata->neventdatascurrent = 0;

   return SCIP_OKAY;
}

/** apply unbalanced branching (see the function \ref enforceCardinality() for further information) */
static
SCIP_RETCODE branchUnbalancedCardinality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_SOL*             sol,                /**< solution to be enforced (or NULL) */
   SCIP_CONS*            branchcons,         /**< cardinality constraint */
   SCIP_VAR**            vars,               /**< variables of constraint */
   SCIP_VAR**            indvars,            /**< indicator variables */
   int                   nvars,              /**< number of variables of constraint */
   int                   cardval,            /**< cardinality value of constraint */
   int                   branchnnonzero,     /**< number of variables that are fixed to be nonzero */
   int                   branchpos           /**< position in array 'vars' */
   )
{
   SCIP_Bool infeasible;
   SCIP_NODE* node1;
   SCIP_NODE* node2;

   /* check whether the variable selected for branching has a nonzero LP solution */
   assert(!SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, vars[branchpos])));
   assert(SCIPisFeasZero(scip, SCIPvarGetLbLocal(indvars[branchpos])));
   assert(SCIPisFeasEQ(scip, SCIPvarGetUbLocal(indvars[branchpos]), 1.0));

   /* create branches */
   SCIPdebugMsg(scip, "apply unbalanced branching on variable <%s> of constraint <%s>.\n",
        SCIPvarGetName(indvars[branchpos]), SCIPconsGetName(branchcons));

   /* create node 1 */

   /* calculate node selection and objective estimate for node 1 */
   SCIP_CALL( SCIPcreateChild(scip, &node1, SCIPcalcNodeselPriority(scip, vars[branchpos], SCIP_BRANCHDIR_DOWNWARDS, 0.0),
         SCIPcalcChildEstimate(scip, vars[branchpos], 0.0) ) );

   /* fix branching variable to zero */
   SCIP_CALL( fixVariableZeroNode(scip, vars[branchpos], node1, &infeasible) );
   assert(! infeasible);

   /* create node 2 */

   /* if the new number of nonzero variables is equal to the number of allowed nonzero variables;
    * i.e. cardinality constraint is saturated */
   assert(branchnnonzero + 1 <= cardval);
   if( branchnnonzero + 1 == cardval )
   {
      SCIP_Real nodeselest;
      SCIP_Real objest;
      int cnt;
      int j;

      /* calculate node selection and objective estimate for node 2 */
      nodeselest = 0.0;
      objest = 0.0;
      cnt = 0;
      for( j = 0; j < nvars; ++j )
      {
         /* we only consider variables in constraint that are not the branching variable and are not fixed to nonzero */
         if( j != branchpos && SCIPvarGetLbLocal(indvars[j]) != 1.0 && !SCIPisFeasPositive(scip, SCIPvarGetLbLocal(vars[j]))
            && !SCIPisFeasNegative(scip, SCIPvarGetUbLocal(vars[j]))
            )
         {
            nodeselest += SCIPcalcNodeselPriority(scip, vars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
            objest += SCIPcalcChildEstimate(scip, vars[j], 0.0);
            ++cnt;
         }
      }
      assert(cnt == nvars - (1 + branchnnonzero));
      assert(cnt > 0);

      /* take the average of the individual estimates */
      objest = objest / (SCIP_Real) cnt;

      /* create node 2 */
      SCIP_CALL( SCIPcreateChild(scip, &node2, nodeselest, objest) );

      /* indicate that branching variable may be treated as nonzero */
      SCIP_CALL( SCIPchgVarLbNode(scip, node2, indvars[branchpos], 1.0) );

      /* fix variables to zero since cardinality constraint is now saturated */
      for( j = 0; j < nvars; ++j )
      {
         /* we only consider variables in constraint that are not the branching variable and are not fixed to nonzero */
         if( j != branchpos && SCIPvarGetLbLocal(indvars[j]) != 1.0
            && !SCIPisFeasPositive(scip, SCIPvarGetLbLocal(vars[j]))
            && !SCIPisFeasNegative(scip, SCIPvarGetUbLocal(vars[j]))
            )
         {
            SCIP_CALL( fixVariableZeroNode(scip, vars[j], node2, &infeasible) );
            assert(!infeasible);
         }
      }
   }
   else
   {
      /* calculate node selection estimate for node 2 */
      SCIP_CALL( SCIPcreateChild(scip, &node2, 0.0, SCIPgetLocalTransEstimate(scip)) );

      /* indicate that branching variable may be treated as nonzero */
      SCIP_CALL( SCIPchgVarLbNode(scip, node2, indvars[branchpos], 1.0) );
   }

   return SCIP_OKAY;
}

/** apply balanced branching (see the function enforceCardinality() for further information) */
static
SCIP_RETCODE branchBalancedCardinality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SOL*             sol,                /**< solution to be enforced (or NULL) */
   SCIP_CONS*            branchcons,         /**< cardinality constraint */
   SCIP_VAR**            vars,               /**< variables of constraint */
   SCIP_VAR**            indvars,            /**< indicator variables */
   int                   nvars,              /**< number of variables of constraint */
   int                   cardval,            /**< cardinality value of constraint */
   int                   branchnnonzero,     /**< number of variables that are fixed to be nonzero */
   int                   branchpos,          /**< position in array 'vars' */
   SCIP_Real             balancedcutoff      /**< cut off value for deciding whether to apply balanced branching */
   )
{
   SCIP_VAR** branchvars;
   SCIP_VAR** branchindvars;
   int nbranchvars;
   SCIP_Real splitval1;
   SCIP_Real splitval2;
   SCIP_Real weight1;
   SCIP_Real weight2;
   SCIP_Real sum1;
   SCIP_Real sum2;
   SCIP_Real w;
   int newcardval;
   int nnonzero;
   int nzero;
   int nbuffer;
   int ind;
   int cnt;
   int j;

   /* check parameters */
   if( SCIPconshdlrGetSepaFreq(conshdlr) != 1 )
   {
      SCIPerrorMessage("balanced branching is only possible if separation frequency of constraint handler is 1.\n");
      return SCIP_PARAMETERWRONGVAL;
   }

   cnt = 0;
   nzero = 0;
   nnonzero = 0;
   nbranchvars = 0;

   weight1 = 0.0;
   weight2 = 0.0;
   sum1 = 0.0;
   sum2 = 0.0;

   /* allocate buffer arrays */
   nbuffer = nvars-branchnnonzero;
   SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, nbuffer) );
   SCIP_CALL( SCIPallocBufferArray(scip, &branchindvars, nbuffer) );

   /* compute weight */
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;

      var = vars[j];

      /* if(binary) indicator variable is not fixed to 1.0 */
      if( SCIPvarGetLbLocal(indvars[j]) != 1.0 && !SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var))
           && !SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var)) )
      {
         /* if implied variable is not already fixed to zero */
         if( !SCIPisFeasZero(scip, SCIPvarGetLbLocal(var)) || !SCIPisFeasZero(scip, SCIPvarGetUbLocal(var)) )
         {
            SCIP_Real val = REALABS(SCIPgetSolVal(scip, sol, var));

            weight1 += val * (SCIP_Real) (j - (nnonzero + nzero));
            weight2 += val;
            branchindvars[nbranchvars] = indvars[j];
            branchvars[nbranchvars++] = var;

            if( !SCIPisFeasZero(scip, val) )
               ++cnt;
         }
         else
            ++nzero;
      }
      else
         ++nnonzero;
   }
   assert(nnonzero == branchnnonzero);
   assert(nbranchvars <= nvars - branchnnonzero);

   assert(cnt >= cardval-nnonzero);
   assert(!SCIPisFeasZero(scip, weight2));
   w = weight1/weight2;  /*lint !e414*/

   ind = (int)SCIPfloor(scip, w);
   assert(0 <= ind && ind < nbranchvars-1);

   /* compute LP sums */
   for( j = 0; j <= ind; ++j )
   {
      SCIP_Real val;

      val = SCIPgetSolVal(scip, sol, branchvars[j]);

      if( SCIPisFeasPositive(scip, val) )
      {
         assert(SCIPisFeasPositive(scip, SCIPvarGetUbLocal(branchvars[j])));
         sum1 += val / SCIPvarGetUbLocal(branchvars[j]);
      }
      else if( SCIPisFeasNegative(scip, val) )
      {
         assert(SCIPisFeasNegative(scip, SCIPvarGetLbLocal(branchvars[j])));
         sum1 += val / SCIPvarGetLbLocal(branchvars[j]);
      }
   }
   for( j = ind+1; j < nbranchvars; ++j )
   {
      SCIP_Real val;

      val = SCIPgetSolVal(scip, sol, branchvars[j]);

      if( SCIPisFeasPositive(scip, val) )
      {
         assert(SCIPisFeasPositive(scip, SCIPvarGetUbLocal(branchvars[j])));
         sum2 += val/SCIPvarGetUbLocal(branchvars[j]);
      }
      else if( SCIPisFeasNegative(scip, val) )
      {
         assert(SCIPisFeasNegative(scip, SCIPvarGetLbLocal(branchvars[j])));
         sum2 += val/SCIPvarGetLbLocal(branchvars[j]);
      }
   }

   /* compute cardinality values of branching constraints */
   newcardval = cardval - nnonzero;
   splitval1 = sum1 + (SCIP_Real)newcardval - sum2 - 1.0;/*lint !e834*/
   splitval1 = SCIPfloor(scip, splitval1/2);
   splitval1 = MAX(splitval1, 0);
   assert((int)splitval1 >= 0);
   assert((int)splitval1 <= MIN(newcardval-1, ind));
   splitval2 = (SCIP_Real)(newcardval-1);
   splitval2 -= splitval1;

   /* the lower or upper LP row of each branching constraint should cut off the current LP solution
    * if this is not the case, then use unbalanced branching */
   if ( !SCIPisFeasLT(scip, (SCIP_Real) splitval1 + balancedcutoff, sum1) ||
         !SCIPisFeasLT(scip, (SCIP_Real) splitval2 + balancedcutoff, sum2) )
   {
      SCIP_CALL( branchUnbalancedCardinality(scip, sol, branchcons, vars, indvars, nvars, cardval,
           branchnnonzero, branchpos) );
   }
   else
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_NODE* node1;
      SCIP_NODE* node2;
      SCIP_CONS* cons1;
      SCIP_CONS* cons2;

      SCIPdebugMsg(scip, "apply balanced branching on constraint <%s>.\n", SCIPconsGetName(branchcons));

      if( SCIPisFeasZero(scip, splitval1) )
      {
         SCIP_Bool infeasible;
         SCIP_Real nodeselest;
         SCIP_Real objest;

         nodeselest = 0.0;
         objest = 0.0;

         /* calculate node selection and objective estimate for node */
         for( j = 0; j <= ind; ++j )
         {
            nodeselest += SCIPcalcNodeselPriority(scip, branchvars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
            objest += SCIPcalcChildEstimate(scip, branchvars[j], 0.0);
         }

         /* take the average of the individual estimates */
         objest = objest/(SCIP_Real)(ind + 1.0);

         /* create node 1 */
         SCIP_CALL( SCIPcreateChild(scip, &node1, nodeselest, objest) );

         for( j = 0; j <= ind; ++j )
         {
            SCIP_CALL( fixVariableZeroNode(scip, branchvars[j], node1, &infeasible) );
            assert(!infeasible);
         }
      }
      else
      {
         /* calculate node selection and objective estimate for node */
         SCIP_CALL( SCIPcreateChild(scip, &node1, 0.0, SCIPgetLocalTransEstimate(scip)) );

         /* create branching constraint for node 1 */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "brleft_#%d", SCIPgetNNodes(scip));
         SCIP_CALL( SCIPcreateConsCardinality(scip, &cons1, name, ind+1, branchvars, (int)splitval1, branchindvars, NULL,
               FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE) );

         /* add constraint to node */
         SCIP_CALL( SCIPaddConsNode(scip, node1, cons1, NULL) );

         /* release constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
      }

      if( SCIPisFeasZero(scip, splitval2) )
      {
         SCIP_Bool infeasible;
         SCIP_Real nodeselest;
         SCIP_Real objest;

         nodeselest = 0.0;
         objest = 0.0;

         /* calculate node selection and objective estimate for node */
         for( j = ind+1; j < nbranchvars; ++j )
         {
            nodeselest += SCIPcalcNodeselPriority(scip, branchvars[j], SCIP_BRANCHDIR_DOWNWARDS, 0.0);
            objest += SCIPcalcChildEstimate(scip, branchvars[j], 0.0);
         }
         assert(nbranchvars - (ind + 1) > 0);

         /* take the average of the individual estimates */
         objest = objest/((SCIP_Real)(nbranchvars - (ind + 1)));/*lint !e414*/

         /* create node 1 */
         SCIP_CALL( SCIPcreateChild(scip, &node2, nodeselest, objest) );

         for( j = ind+1; j < nbranchvars; ++j )
         {
            SCIP_CALL( fixVariableZeroNode(scip, branchvars[j], node2, &infeasible) );
            assert(!infeasible);
         }
      }
      else
      {
         /* calculate node selection and objective estimate for node */
         SCIP_CALL( SCIPcreateChild(scip, &node2, 0.0, SCIPgetLocalTransEstimate(scip)) );

         /* shift the second half of variables */
         cnt = 0;
         for( j = ind+1; j < nbranchvars; ++j )
         {
            branchvars[cnt] = branchvars[j];
            branchindvars[cnt++] = branchindvars[j];
         }
         assert(cnt == nbranchvars - (ind + 1));

         /* create branching constraint for node 2 */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "brright_#%d", SCIPgetNNodes(scip));
         SCIP_CALL( SCIPcreateConsCardinality(scip, &cons2, name, cnt, branchvars, (int)splitval2, branchindvars, NULL,
               FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE) );

         /* add constraint to node */
         SCIP_CALL( SCIPaddConsNode(scip, node2, cons2, NULL) );

         /* release constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &branchindvars);
   SCIPfreeBufferArray(scip, &branchvars);

   return SCIP_OKAY;
}

/** enforcement method
 *
 *  We check whether the current solution is feasible. If not, the cardinality constraints can be enforced by different
 *  branching rules:
 *
 *  - Unbalanced branching: Branch on the neighborhood of a single variable \f$i\f$, i.e., in one branch \f$x_i\f$ is
 *    fixed to zero and in the other we modify cardinality constraints \f$|\mbox{supp}(x)| \leq k\f$ with \f$i\in D\f$ to
 *    \f$|\mbox{supp}(x_{D\setminus i}) \leq k-1\f$
 *
 *  - Balanced branching: First, choose a cardinality constraint \f$|\mbox{supp}(x_D) \leq k\f$ that is violated by the
 *    current LP solution. Then, we compute \f$W = \sum_{j=1}^n |x_i|\f$ and \f$w = \sum_{j=1}^n j\, |x_i|\f$. Next,
 *    search for the index \f$r\f$ that satisfies
 *    \f[
 *        r \leq \frac{w}{W} < r+1.
 *    \f]
 *    Choose a number \f$s\f$ with \f$0\leq s < \min\{k, r\}\f$. The branches are then
 *    \f[
 *        |\mbox{supp}(x_{d_1}, \ldots, x_{d_r})| \leq s \qquad \mbox{and}\qquad
 *        |\mbox{supp}(x_{d_{r+1}}, \ldots, x_{d_{n}})| \leq k-s-1,
 *  \f]
 *  where \f$d_1, \ldots, d_n\f$ are the elements of the set \f$D\f$.
 *
 * The branching constraint is chosen by the largest sum of variable values.
 */
static
SCIP_RETCODE enforceCardinality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SOL*             sol,                /**< solution to be enforced (or NULL) */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< indicator constraints */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* branchcons;
   SCIP_Real maxweight;
   SCIP_VAR** indvars;
   SCIP_VAR** vars;
   int nvars;
   int cardval;

   SCIP_Bool branchbalanced = FALSE;
   SCIP_Bool branchallpos = TRUE;
   SCIP_Bool branchallneg = TRUE;
   int branchnnonzero;
   int branchpos;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(result != NULL);

   maxweight = -SCIP_REAL_MAX;
   branchcons = NULL;
   branchnnonzero = -1;
   branchpos = -1;

   SCIPdebugMsg(scip, "Enforcing cardinality constraints <%s>.\n", SCIPconshdlrGetName(conshdlr) );
   *result = SCIP_FEASIBLE;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* search for a constraint with largest violation; from this constraint, we select the variable with largest LP value */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_Bool cutoff;
      SCIP_Real weight;
      SCIP_Real maxval;
      SCIP_Bool allpos = TRUE;
      SCIP_Bool allneg = TRUE;
      int nnonzero;    /* number of variables that are currently deactivated in constraint */
      int pos;         /* position of variable with largest LP solution value */
      int nchgdomain;
      int cnt;
      int j;

      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      nchgdomain = 0;
      cnt = 0;
      nnonzero = 0;
      pos = -1;
      nvars = consdata->nvars;
      vars = consdata->vars;
      indvars = consdata->indvars;
      cardval = consdata->cardval;

      /* do nothing if there are not enough variables - this is usually eliminated by preprocessing */
      if( nvars < 2 )
         continue;

      /* first perform propagation (it might happen that standard propagation is turned off) */
      SCIP_CALL( propCardinality(scip, cons, consdata, &cutoff, &nchgdomain) );

      SCIPdebugMsg(scip, "propagating <%s> in enforcing (cutoff: %u, domain reductions: %d).\n",
           SCIPconsGetName(cons), cutoff, nchgdomain);
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( nchgdomain > 0 )
      {
         *result = SCIP_REDUCEDDOM;
         return SCIP_OKAY;
      }
      assert(nchgdomain == 0);

      /* check constraint */
      weight = 0.0;
      maxval = -SCIPinfinity(scip);

      for( j = 0; j < nvars; ++j )
      {
         SCIP_VAR* var;

         /* check whether indicator variable is zero, but variable in cardinality constraint is not fixed to zero;
          * if the variable is not multiaggregated this case should already be handled in propagation */
         if( SCIPvarGetUbLocal(indvars[j]) == 0.0 &&
             (
                !SCIPisFeasZero(scip, SCIPvarGetLbLocal(vars[j])) || !SCIPisFeasZero(scip, SCIPvarGetUbLocal(vars[j]))
             )
           )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         assert(SCIPvarGetStatus(indvars[j]) != SCIP_VARSTATUS_MULTAGGR);

         var = vars[j];

         /* variable is not fixed to nonzero */
         if( SCIPvarGetLbLocal(indvars[j]) != 1.0
             && !SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var))
             && !SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var))
           )
         {
            SCIP_Real val;

            val = SCIPgetSolVal(scip, sol, var);
            if( SCIPisFeasPositive(scip, val))
               allneg = FALSE;
            else if( SCIPisFeasNegative(scip, val))
               allpos = FALSE;
            val = REALABS(val);

            if( !SCIPisFeasZero(scip, val) )
            {
               /* determine maximum nonzero-variable solution value */
               if( SCIPisFeasGT(scip, val, maxval) )
               {
                  pos = j;
                  maxval = val;
               }

               weight += val;
               ++cnt;
            }
         }
         else
            ++nnonzero;
      }
      weight -= cardval;
      weight += nnonzero;

      /* if we detected a cut off */
      if( nnonzero > cardval )
      {
         SCIPdebugMsg(scip, "Detected cut off: constraint <%s> has %d many variables that can be treated as nonzero, \
            although only %d many are feasible.\n", SCIPconsGetName(cons), nnonzero, cardval);
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      /* else if domain can be reduced (since node 2 created in branchUnbalancedCardinality() would be infeasible) */
      else if( cnt > 0 && nnonzero + 1 > cardval )
      {
         SCIP_Bool infeasible;
         int v;

         for( v = 0; v < nvars; ++v )
         {
            SCIP_VAR* var;

            var = vars[v];

            /* variable is not fixed to nonzero */
            if( !SCIPisFeasEQ(scip, SCIPvarGetLbLocal(indvars[v]), 1.0)
                && !SCIPisFeasPositive(scip, SCIPvarGetLbLocal(var))
                && !SCIPisFeasNegative(scip, SCIPvarGetUbLocal(var))
              )
            {
               SCIP_CALL( fixVariableZeroNode(scip, var, SCIPgetCurrentNode(scip), &infeasible) );
               assert(!infeasible);
               SCIPdebugMsg(scip, "detected domain reduction in enforcing: fixed variable <%s> to zero.\n", SCIPvarGetName(var));
            }
         }

         *result = SCIP_REDUCEDDOM;
         return SCIP_OKAY;
      }

      /* if constraint is violated */
      if( cnt > cardval - nnonzero && weight > maxweight )
      {
         maxweight = weight;
         branchcons = cons;
         branchnnonzero = nnonzero;
         branchpos = pos;
         branchallneg = allneg;
         branchallpos = allpos;
      }
   }

   /* if all constraints are feasible */
   if( branchcons == NULL )
   {
      SCIP_SOL* primsol;
      SCIP_Bool success;

      /* polish primal solution */
      SCIP_CALL( SCIPcreateSolCopy(scip, &primsol, sol) );
      SCIP_CALL( polishPrimalSolution(scip, conss, nconss, sol, primsol) );
      SCIP_CALL( SCIPtrySol(scip, primsol, FALSE, TRUE, FALSE, TRUE, FALSE, &success) );
      SCIP_CALL( SCIPfreeSol(scip, &primsol) );

      SCIPdebugMsg(scip, "All cardinality constraints are feasible.\n");
      return SCIP_OKAY;
   }
   assert(branchnnonzero >= 0);
   assert(branchpos >= 0);

   /* get data for branching or domain reduction */
   consdata = SCIPconsGetData(branchcons);
   assert(consdata != NULL);
   nvars = consdata->nvars;
   vars = consdata->vars;
   indvars = consdata->indvars;
   cardval = consdata->cardval;

   /* we only use balanced branching if either the lower or the upper bound row of the branching constraint is known
    * to be tight or violated */
   if( conshdlrdata->branchbalanced && !SCIPisFeasNegative(scip, maxweight) && ( branchallneg || branchallpos )
       && (conshdlrdata->balanceddepth == -1 || SCIPgetDepth(scip) <= conshdlrdata->balanceddepth)
     )
   {
         branchbalanced = TRUE;
   }

   /* apply branching rule */
   if( branchbalanced )
   {
      SCIP_CALL( branchBalancedCardinality(scip, conshdlr, sol, branchcons, vars, indvars, nvars, cardval, branchnnonzero, branchpos,
           conshdlrdata->balancedcutoff) );
   }
   else
   {
      SCIP_CALL( branchUnbalancedCardinality(scip, sol, branchcons, vars, indvars, nvars, cardval, branchnnonzero,
           branchpos) );
   }

   SCIP_CALL( SCIPresetConsAge(scip, branchcons) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** Generate row
 *
 *  We generate the row corresponding to the following simple valid inequalities:
 *  \f[
 *         \frac{x_1}{u_1} + \ldots + \frac{x_n}{u_n} \leq k\qquad\mbox{and}\qquad
 *         \frac{x_1}{\ell_1} + \ldots + \frac{x_n}{\ell_1} \leq k,
 *  \f]
 *  where \f$\ell_1, \ldots, \ell_n\f$ and \f$u_1, \ldots, u_n\f$ are the nonzero and finite lower and upper bounds of
 *  the variables \f$x_1, \ldots, x_n\f$ and k is the right hand side of the cardinality constraint. If at least k upper
 *  bounds < 0 or a lower bounds > 0, the constraint itself is redundant, so the cut is not applied (lower bounds > 0
 *  and upper bounds < 0 are usually detected in presolving or propagation). Infinite bounds and zero are skipped. Thus
 *  \f$\ell_1, \ldots, \ell_n\f$ are all negative, which results in the \f$\leq\f$ inequality.
 *
 *  Note that in fact, any mixture of nonzero finite lower and upper bounds would lead to a valid inequality as
 *  above. However, usually either the lower or upper bound is nonzero. Thus, the above inequalities are the most
 *  interesting.
 */
static
SCIP_RETCODE generateRowCardinality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             local,              /**< produce local cut? */
   SCIP_ROW**            rowlb,              /**< output: row for lower bounds (or NULL if not needed) */
   SCIP_ROW**            rowub               /**< output: row for upper bounds (or NULL if not needed) */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real val;
   int nvars;
   int cnt;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->vars != NULL);
   assert(consdata->indvars != NULL);

   nvars = consdata->nvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );

   /* take care of upper bounds */
   if( rowub != NULL )
   {
      int cardval;

      cnt = 0;
      cardval = consdata->cardval;
      for( j = 0; j < nvars; ++j )
      {
         if( local )
            val = SCIPvarGetLbLocal(consdata->vars[j]);
         else
            val = SCIPvarGetUbGlobal(consdata->vars[j]);

         /* if a variable may be treated as nonzero, then update cardinality value */
         if( SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->indvars[j]), 1.0) )
         {
            --cardval;
            continue;
         }

         if( !SCIPisInfinity(scip, val) && !SCIPisZero(scip, val) && !SCIPisNegative(scip, val) )
         {
            assert(consdata->vars[j] != NULL);
            vars[cnt] = consdata->vars[j];
            vals[cnt++] = 1.0/val;
         }
      }
      assert(cardval >= 0);

      /* if cut is meaningful */
      if( cnt > cardval )
      {
         /* create upper bound inequality if at least two of the bounds are finite and nonzero */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cardub#%s", SCIPconsGetName(cons));
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowub, conshdlr, name, -SCIPinfinity(scip), (SCIP_Real)cardval,
              local, TRUE, FALSE) );
         SCIP_CALL( SCIPaddVarsToRow(scip, *rowub, cnt, vars, vals) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowub, NULL) ) );
      }
   }

   /* take care of lower bounds */
   if( rowlb != NULL )
   {
      int cardval;

      cnt = 0;
      cardval = consdata->cardval;
      for( j = 0; j < nvars; ++j )
      {
         if( local )
            val = SCIPvarGetLbLocal(consdata->vars[j]);
         else
            val = SCIPvarGetLbGlobal(consdata->vars[j]);

         /* if a variable may be treated as nonzero, then update cardinality value */
         if( SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->indvars[j]), 1.0) )
         {
            --cardval;
            continue;
         }

         if( !SCIPisInfinity(scip, -val) && !SCIPisZero(scip, val) && !SCIPisPositive(scip, val) )
         {
            assert(consdata->vars[j] != NULL);
            vars[cnt] = consdata->vars[j];
            vals[cnt++] = 1.0/val;
         }
      }
      assert(cardval >= 0);

      /* if cut is meaningful */
      if( cnt > cardval )
      {
         /* create lower bound inequality if at least two of the bounds are finite and nonzero */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cardlb#%s", SCIPconsGetName(cons));
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowlb, conshdlr, name, -SCIPinfinity(scip), (SCIP_Real)cardval,
              local, TRUE, FALSE) );
         SCIP_CALL( SCIPaddVarsToRow(scip, *rowlb, nvars, vars, vals) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *rowlb, NULL) ) );
      }
   }

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** initialize or separate bound inequalities from cardinality constraints
 *  (see the function \ref generateRowCardinality() for an explanation of bound inequalities)
 */
static
SCIP_RETCODE initsepaBoundInequalityFromCardinality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< cardinality constraints */
   int                   nconss,             /**< number of cardinality constraints */
   SCIP_SOL*             sol,                /**< LP solution to be separated (or NULL) */
   SCIP_Bool             solvedinitlp,       /**< TRUE if initial LP relaxation at a node is solved */
   int*                  ngen,               /**< pointer to store number of cuts generated (or NULL) */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff occurred */
   )
{
   int cnt = 0;
   int c;

   assert(scip != NULL);
   assert(conss != NULL);

   *cutoff = FALSE;

   for( c = nconss-1; c >= 0; --c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_ROW* rowub = NULL;
      SCIP_ROW* rowlb = NULL;
      SCIP_Bool release = FALSE;

      assert(conss != NULL);
      assert(conss[c] != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( solvedinitlp )
         SCIPdebugMsg(scip, "Separating inequalities for cardinality constraint <%s>.\n", SCIPconsGetName(conss[c]) );
      else
         SCIPdebugMsg(scip, "Checking for initial rows for cardinality constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* generate rows associated to cardinality constraint; the rows are stored in the constraint data
       * if they are globally valid */
      if( SCIPconsIsLocal(conss[c]) )
      {
         SCIP_CALL( generateRowCardinality(scip, conshdlr, conss[c], TRUE, &rowlb, &rowub) );
         release = TRUE;
      }
      else
      {
         if( consdata->rowub == NULL || consdata->rowlb == NULL )
         {
            SCIP_CALL( generateRowCardinality(scip, conshdlr, conss[c], FALSE,
                 (consdata->rowlb == NULL) ? &consdata->rowlb : NULL,
                 (consdata->rowub == NULL) ? &consdata->rowub : NULL) );/*lint !e826*/
         }
         rowub = consdata->rowub;
         rowlb = consdata->rowlb;
      }

      /* put corresponding rows into LP */
      if( rowub != NULL && !SCIProwIsInLP(rowub) && ( solvedinitlp || SCIPisCutEfficacious(scip, sol, rowub) ) )
      {
         assert(SCIPisInfinity(scip, -SCIProwGetLhs(rowub)));
         assert(SCIPisLE(scip, SCIProwGetRhs(rowub), (SCIP_Real)consdata->cardval));

         SCIP_CALL( SCIPaddRow(scip, rowub, FALSE, cutoff) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rowub, NULL) ) );

         if( solvedinitlp )
         {
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
         }
         ++cnt;
      }

      if( ! (*cutoff) && rowlb != NULL && !SCIProwIsInLP(rowlb)
          && ( solvedinitlp || SCIPisCutEfficacious(scip, sol, rowlb) )
        )
      {
         assert(SCIPisInfinity(scip, -SCIProwGetLhs(rowlb)));
         assert(SCIPisLE(scip, SCIProwGetRhs(rowlb), (SCIP_Real)consdata->cardval));

         SCIP_CALL( SCIPaddRow(scip, rowlb, FALSE, cutoff) );
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, rowlb, NULL) ) );

         if( solvedinitlp )
         {
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
         }
         ++cnt;
      }

      if( release )
      {
         if( rowlb != NULL )
         {
            SCIP_CALL( SCIPreleaseRow(scip, &rowlb) );
         }
         if( rowub != NULL )
         {
            SCIP_CALL( SCIPreleaseRow(scip, &rowub) );
         }
      }

      if( *cutoff )
         break;
   }

   /* store number of generated cuts */
   if( ngen != NULL )
      *ngen = cnt;

   return SCIP_OKAY;
}

/** separates cardinality constraints for arbitrary solutions */
static
SCIP_RETCODE separateCardinality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SOL*             sol,                /**< solution to be separated (or NULL) */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< cardinality constraints */
   SCIP_RESULT*          result              /**< result */
   )
{
   SCIP_Bool cutoff;
   int ngen = 0;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( nconss == 0 )
      return SCIP_OKAY;

   /* only separate cuts if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* separate bound inequalities from cardinality constraints */
   SCIP_CALL( initsepaBoundInequalityFromCardinality(scip, conshdlr, conss, nconss, sol, TRUE, &ngen, &cutoff) );
   if( cutoff )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* evaluate results */
   if( ngen > 0 )
      *result = SCIP_SEPARATED;
   SCIPdebugMsg(scip, "Separated %d bound inequalities.\n", ngen);

   return SCIP_OKAY;
}

/* ---------------------------- constraint handler callback methods ----------------------*/

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyCardinality)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrCardinality(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeCardinality)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free hash map */
   if( conshdlrdata->varhash != NULL )
   {
      SCIPhashmapFree(&conshdlrdata->varhash);
   }

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolCardinality)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check each constraint */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      assert(conss != NULL);
      assert(conss[c] != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      SCIPdebugMsg(scip, "Exiting cardinality constraint <%s>.\n", SCIPconsGetName(conss[c]) );

      /* free rows */
      if( consdata->rowub != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rowub) );
      }
      if( consdata->rowlb != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rowlb) );
      }
   }

   /* free hash map */
   if( conshdlrdata->varhash != NULL )
   {
      SCIPhashmapFree(&conshdlrdata->varhash);
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteCardinality)
{ /*lint --e{737, 647}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIPdebugMsg(scip, "Deleting cardinality constraint <%s>.\n", SCIPconsGetName(cons) );

   /* drop events on transformed variables */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int j;

      /* get constraint handler data */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      for( j = 0; j < (*consdata)->nvars; ++j )
      {
         SCIP_CALL( dropVarEventCardinality(scip, conshdlrdata->eventhdlr, *consdata, (*consdata)->vars[j],
              (*consdata)->indvars[j], &(*consdata)->eventdatas[j]) );
         assert((*consdata)->eventdatas[j] == NULL);
      }
   }

   if( (*consdata)->weights != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->weights, (*consdata)->maxvars);
   }
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventdatas, (*consdata)->maxvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventvarscurrent, 4 * (*consdata)->maxvars);/*lint !e647*/
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->eventdatascurrent, 4 * (*consdata)->maxvars);/*lint !e647*/
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->indvars, (*consdata)->maxvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->maxvars);

   /* free rows */
   if( (*consdata)->rowub != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->rowub) );
   }
   if( (*consdata)->rowlb != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->rowlb) );
   }
   assert((*consdata)->rowub == NULL);
   assert((*consdata)->rowlb == NULL);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransCardinality)
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   char s[SCIP_MAXSTRLEN];
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   SCIPdebugMsg(scip, "Transforming cardinality constraint: <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->nvars > 0);
   assert(sourcedata->nvars <= sourcedata->maxvars);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->cons = NULL;
   consdata->nvars = sourcedata->nvars;
   consdata->maxvars = sourcedata->nvars;
   consdata->cardval = sourcedata->cardval;
   consdata->rowub = NULL;
   consdata->rowlb = NULL;
   consdata->eventdatascurrent = NULL;
   consdata->neventdatascurrent = 0;
   consdata->ntreatnonzeros = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, consdata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->indvars, consdata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventdatas, consdata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventdatascurrent, 4*consdata->nvars) );/*lint !e647*/
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventvarscurrent, 4*consdata->nvars) );/*lint !e647*/

   /* if weights were used */
   if( sourcedata->weights != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->weights, sourcedata->weights, consdata->nvars) );
   }
   else
      consdata->weights = NULL;

   for( j = 0; j < sourcedata->nvars; ++j )
   {
      assert(sourcedata->vars[j] != 0);
      assert(sourcedata->indvars[j] != 0);
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->vars[j], &(consdata->vars[j])) );
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->indvars[j], &(consdata->indvars[j])) );

      /* if variable is fixed to be nonzero */
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->indvars[j]), 1.0) )
         ++(consdata->ntreatnonzeros);
   }

   /* create transformed constraint with the same flags */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   consdata->cons = *targetcons;
   assert(consdata->cons != NULL);

   /* catch bound change events on variable */
   for( j = 0; j < consdata->nvars; ++j )
   {
      SCIP_CALL( catchVarEventCardinality(scip, conshdlrdata->eventhdlr, consdata,
           consdata->vars[j], consdata->indvars[j], j, &consdata->eventdatas[j]) );
      assert(consdata->eventdatas[j] != NULL);
   }

#ifdef SCIP_DEBUG
   if( SCIPisGT(scip, (SCIP_Real)consdata->ntreatnonzeros, consdata->cardval) )
   {
      SCIPdebugMsg(scip, "constraint <%s> has %d variables fixed to be nonzero, allthough the constraint allows \
           only %d nonzero variables\n", SCIPconsGetName(*targetcons), consdata->ntreatnonzeros, consdata->cardval);
   }
#endif

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolCardinality)
{  /*lint --e{715}*/
   /* cppcheck-suppress unassignedVariable */
   int oldnfixedvars;
   /* cppcheck-suppress unassignedVariable */
   int oldndelconss;
   /* cppcheck-suppress unassignedVariable */
   int oldnupgdconss;
   int nremovedvars;
   SCIP_EVENTHDLR* eventhdlr;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Presolving cardinality constraints.\n");

   *result = SCIP_DIDNOTRUN;
   SCIPdebug( oldnfixedvars = *nfixedvars; )
   SCIPdebug( oldndelconss = *ndelconss; )
   SCIPdebug( oldnupgdconss = *nupgdconss; )
   nremovedvars = 0;

   /* only run if success if possible */
   if( nrounds == 0 || nnewfixedvars > 0 || nnewaggrvars > 0 )
   {
      /* get constraint handler data */
      assert(SCIPconshdlrGetData(conshdlr) != NULL);
      eventhdlr = SCIPconshdlrGetData(conshdlr)->eventhdlr;
      assert(eventhdlr != NULL);

      *result = SCIP_DIDNOTFIND;

      /* check each constraint */
      for( c = 0; c < nconss; ++c )
      {
         SCIP_CONSDATA* consdata;
         SCIP_CONS* cons;
         SCIP_Bool cutoff;
         SCIP_Bool success;

         assert(conss != NULL);
         assert(conss[c] != NULL);
         cons = conss[c];
         consdata = SCIPconsGetData(cons);

         assert(consdata != NULL);
         assert(consdata->nvars >= 0);
         assert(consdata->nvars <= consdata->maxvars);
         assert(!SCIPconsIsModifiable(cons));

         /* perform one presolving round */
         SCIP_CALL( presolRoundCardinality(scip, cons, consdata, eventhdlr, &cutoff, &success,
              ndelconss, nupgdconss, nfixedvars, &nremovedvars) );

         if( cutoff )
         {
            SCIPdebugMsg(scip, "presolving detected cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         if( success )
            *result = SCIP_SUCCESS;
      }
   }
   (*nchgcoefs) += nremovedvars;

   SCIPdebugMsg(scip, "presolving fixed %d variables, removed %d variables, deleted %d constraints, \
        and upgraded %d constraints.\n", *nfixedvars - oldnfixedvars, nremovedvars, *ndelconss - oldndelconss,
        *nupgdconss - oldnupgdconss);

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpCardinality)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   /* checking for initial rows for cardinality constraints */
   SCIP_CALL( initsepaBoundInequalityFromCardinality(scip, conshdlr, conss, nconss, NULL, FALSE, NULL, &cutoff) );
   assert(!cutoff);

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpCardinality)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(result != NULL);

   SCIP_CALL( separateCardinality(scip, conshdlr, NULL, nconss, conss, result) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolCardinality)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(result != NULL);

   SCIP_CALL( separateCardinality(scip, conshdlr, sol, nconss, conss, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpCardinality)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   SCIP_CALL( enforceCardinality(scip, conshdlr, NULL, nconss, conss, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxCardinality)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIP_CALL( enforceCardinality(scip, conshdlr, sol, nconss, conss, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsCardinality)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   SCIP_CALL( enforceCardinality(scip, conshdlr, NULL, nconss, conss, result) );

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions
 *
 *  We simply check whether more variables than allowed are nonzero in the given solution.
 */
static
SCIP_DECL_CONSCHECK(consCheckCardinality)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* check each constraint */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      int cardval;
      int j;
      int cnt;

      cnt = 0;
      assert(conss[c] != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      cardval = consdata->cardval;
      SCIPdebugMsg(scip, "Checking cardinality constraint <%s>.\n", SCIPconsGetName(conss[c]));

      /* check all variables */
      for( j = 0; j < consdata->nvars; ++j )
      {
         /* if variable is nonzero */
         if( !SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[j])) )
         {
            ++cnt;

            /* if more variables than allowed are nonzero */
            if( cnt > cardval )
            {
               SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
               *result = SCIP_INFEASIBLE;

               if( printreason )
               {
                  int l;

                  SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
                  SCIPinfoMessage(scip, NULL, ";\nviolation: ");

                  for( l = 0; l < consdata->nvars; ++l )
                  {
                     /* if variable is nonzero */
                     if( !SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vars[l])) )
                     {
                        SCIPinfoMessage(scip, NULL, "<%s> = %.15g ",
                           SCIPvarGetName(consdata->vars[l]), SCIPgetSolVal(scip, sol, consdata->vars[l]));
                     }
                  }
                  SCIPinfoMessage(scip, NULL, "\n");
               }
               if( sol != NULL )
                  SCIPupdateSolConsViolation(scip, sol, 1.0, 1.0);
               return SCIP_OKAY;
            }
         }
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropCardinality)
{  /*lint --e{715}*/
   int nchgdomain = 0;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   assert(SCIPisTransformed(scip));

   /* check each constraint */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      SCIP_Bool cutoff;

      *result = SCIP_DIDNOTFIND;
      assert(conss[c] != NULL);
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      SCIPdebugMsg(scip, "Propagating cardinality constraint <%s>.\n", SCIPconsGetName(cons) );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( propCardinality(scip, cons, consdata, &cutoff, &nchgdomain) );
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   SCIPdebugMsg(scip, "Propagated %d domains.\n", nchgdomain);
   if( nchgdomain > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler
 *
 *  Let lb and ub be the lower and upper bounds of a
 *  variable. Preprocessing usually makes sure that lb <= 0 <= ub.
 *
 *  - If lb < 0 then rounding down may violate the constraint.
 *  - If ub > 0 then rounding up may violated the constraint.
 *  - If lb > 0 or ub < 0 then the rhs of the constraint can be updated and the variable
 *    can be removed from the constraint. Thus, we do not have to deal with it here.
 *  - If lb == 0 then rounding down does not violate the constraint.
 *  - If ub == 0 then rounding up does not violate the constraint.
 */
static
SCIP_DECL_CONSLOCK(consLockCardinality)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   SCIP_VAR** indvars;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMsg(scip, "Locking constraint <%s>.\n", SCIPconsGetName(cons));

   vars = consdata->vars;
   indvars = consdata->indvars;
   nvars = consdata->nvars;
   assert(vars != NULL);

   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      SCIP_VAR* indvar;
      var = vars[j];
      indvar = indvars[j];

      /* if lower bound is negative, rounding down may violate constraint */
      if( SCIPisFeasNegative(scip, SCIPvarGetLbLocal(var)) )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, var, nlockspos, nlocksneg) );
      }

      /* additionally: if upper bound is positive, rounding up may violate constraint */
      if( SCIPisFeasPositive(scip, SCIPvarGetUbLocal(var)) )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, var, nlocksneg, nlockspos) );
      }

      /* add lock on indicator variable; @todo write constraint handler to handle down locks */
      SCIP_CALL( SCIPaddVarLocks(scip, indvar, nlockspos, nlockspos) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintCardinality)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( j = 0; j < consdata->nvars; ++j )
   {
      if( j > 0 )
         SCIPinfoMessage(scip, file, ", ");
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[j], FALSE) );
      if( consdata->weights == NULL )
         SCIPinfoMessage(scip, file, " (%d)", j+1);
      else
         SCIPinfoMessage(scip, file, " (%3.2f)", consdata->weights[j]);
   }
   SCIPinfoMessage(scip, file, " <= %d", consdata->cardval);

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyCardinality)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourceconsdata;
   SCIP_VAR** sourcevars;
   SCIP_VAR** targetvars;
   SCIP_VAR** sourceindvars;
   SCIP_VAR** targetindvars;
   SCIP_Real* sourceweights;
   SCIP_Real* targetweights;
   const char* consname;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);
   assert(SCIPisTransformed(sourcescip));
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0);

   *valid = TRUE;

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIPdebugMsg(scip, "Copying cardinality constraint <%s> ...\n", consname);

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get variables and weights of the source constraint */
   nvars = sourceconsdata->nvars;

   if( nvars == 0 )
      return SCIP_OKAY;

   sourcevars = sourceconsdata->vars;
   assert(sourcevars != NULL);
   sourceindvars = sourceconsdata->indvars;
   assert(sourceindvars != NULL);
   sourceweights = sourceconsdata->weights;
   assert(sourceweights != NULL);

   /* duplicate variable array */
   SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(sourcescip, &targetindvars, nvars) );
   SCIP_CALL( SCIPduplicateBufferArray(sourcescip, &targetweights, sourceweights, nvars) );

   /* get copied variables in target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      assert(sourcevars[v] != NULL);
      assert(sourceindvars[v] != NULL);

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &(targetvars[v]), varmap, consmap, global, valid) );
      if( *valid )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourceindvars[v], &(targetindvars[v]), varmap, consmap, global, valid) );
      }
   }

    /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsCardinality(scip, cons, consname, nvars, targetvars, sourceconsdata->cardval, targetindvars,
            targetweights, initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(sourcescip, &targetweights);
   SCIPfreeBufferArray(sourcescip, &targetindvars);
   SCIPfreeBufferArray(sourcescip, &targetvars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseCardinality)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_Real weight;
   int cardval;
   const char* s;
   char* t;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert(cons != NULL);
   assert(success != NULL);

   *success = TRUE;
   s = str;

   /* create empty cardinality constraint */
   SCIP_CALL( SCIPcreateConsCardinality(scip, cons, name, 0, NULL, 0, NULL, NULL, initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );

   /* loop through string */
   do
   {
      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, s, &var, &t) );
      s = t;

      /* skip until beginning of weight */
      while ( *s != '\0' && *s != '(' )
         ++s;
 
      if ( *s == '\0' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected weight at input: %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      /* skip '(' */
      ++s;

      /* find weight */
      weight = strtod(s, &t);
      if ( t == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error during parsing of the weight: %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      s = t;

      /* skip white space, ',', and ')' */
      while ( *s != '\0' && ( isspace((unsigned char)*s) ||  *s == ',' || *s == ')' ) )
         ++s;

      /* add variable */
      SCIP_CALL( SCIPaddVarCardinality(scip, *cons, var, NULL, weight) );

      /* check if there is a '<=' */
      if ( *s == '<' && *(s+1) == '='  )
      {
         s = s + 2;

         /* skip white space */
         while ( isspace((unsigned char)*s) )
            ++s;

         cardval = (int)strtod(s, &t);
         if ( t == NULL )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error during parsing of the cardinality restriction value: %s\n", s);
            *success = FALSE;
            return SCIP_OKAY;
         }
         s = t;
              
         SCIP_CALL( SCIPchgCardvalCardinality(scip, *cons, cardval));
      }
   }
   while ( *s != '\0' );
  
   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsCardinality)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsCardinality)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * update the number of variables fixed to be nonzero
 * update the bound constraints
 */
static
SCIP_DECL_EVENTEXEC(eventExecCardinality)
{
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSDATA* consdata;
   SCIP_Real oldbound;
   SCIP_Real newbound;
   SCIP_VAR* var;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   consdata = eventdata->consdata;
   assert(consdata != NULL);
   assert(0 <= consdata->ntreatnonzeros && consdata->ntreatnonzeros <= consdata->nvars);
   assert(consdata->eventdatascurrent != NULL);
   assert(consdata->eventvarscurrent != NULL);

   var = SCIPeventGetVar(event);
   assert(var != NULL);
   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);
   eventtype = SCIPeventGetType(event);

#ifdef SCIP_DEBUG
   if( ( eventdata->varmarked && var == eventdata->var) || ( eventdata->indvarmarked && var == eventdata->indvar)  )
   {
      int i;

      for( i = 0; i < consdata->neventdatascurrent; ++i )
      {
         if( var == consdata->eventvarscurrent[i] )
         {
            break;
         }
      }
      assert(i < consdata->neventdatascurrent);
   }
#endif

   /* if variable is an indicator variable */
   if( var == eventdata->indvar )
   {
      assert(SCIPvarIsBinary(var));
      assert(consdata->cons != NULL);

      if( eventtype == SCIP_EVENTTYPE_LBTIGHTENED || eventtype == SCIP_EVENTTYPE_LBRELAXED )
      {
         if( eventtype == SCIP_EVENTTYPE_LBTIGHTENED )
            ++(consdata->ntreatnonzeros);
         else if( eventtype == SCIP_EVENTTYPE_LBRELAXED )
            --(consdata->ntreatnonzeros);

         assert(0 <= consdata->ntreatnonzeros && consdata->ntreatnonzeros <= consdata->nvars);
      }
      else if( eventtype == SCIP_EVENTTYPE_UBTIGHTENED && ! eventdata->indvarmarked )
      {
         assert(oldbound == 1.0 && newbound == 0.0 );

         /* save event data for propagation */
         consdata->eventdatascurrent[consdata->neventdatascurrent] = eventdata;
         consdata->eventvarscurrent[consdata->neventdatascurrent] = var;
         ++consdata->neventdatascurrent;
         eventdata->indvarmarked = TRUE;
         assert(consdata->neventdatascurrent <= 4 * consdata->maxvars);
         assert(var == eventdata->indvar );
      }
   }

   /* if variable is an implied variable,
    * notice that the case consdata->var = consdata->indvar is possible */
   if( var == eventdata->var && ! eventdata->varmarked )
   {
      if( eventtype == SCIP_EVENTTYPE_LBTIGHTENED )
      {
         /* if variable is now fixed to be nonzero */
         if( !SCIPisFeasPositive(scip, oldbound) && SCIPisFeasPositive(scip, newbound) )
         {
            /* save event data for propagation */
            consdata->eventdatascurrent[consdata->neventdatascurrent] = eventdata;
            consdata->eventvarscurrent[consdata->neventdatascurrent] = var;
            ++consdata->neventdatascurrent;
            eventdata->varmarked = TRUE;
            assert(consdata->neventdatascurrent <= 4 * consdata->maxvars );
            assert(var == eventdata->var );
         }
      }
      else if( eventtype == SCIP_EVENTTYPE_UBTIGHTENED )
      {
         /* if variable is now fixed to be nonzero */
         if( !SCIPisFeasNegative(scip, oldbound) && SCIPisFeasNegative(scip, newbound) )
         {
            /* save event data for propagation */
            consdata->eventdatascurrent[consdata->neventdatascurrent] = eventdata;
            consdata->eventvarscurrent[consdata->neventdatascurrent] = var;
            ++consdata->neventdatascurrent;
            eventdata->varmarked = TRUE;
            assert(consdata->neventdatascurrent <= 4 * consdata->maxvars );
            assert(var == eventdata->var);
         }
      }
   }
   assert(0 <= consdata->ntreatnonzeros && consdata->ntreatnonzeros <= consdata->nvars);

   SCIPdebugMsg(scip, "event exec cons <%s>: changed bound of variable <%s> from %f to %f (ntreatnonzeros: %d).\n",
        SCIPconsGetName(consdata->cons), SCIPvarGetName(SCIPeventGetVar(event)),
        oldbound, newbound, consdata->ntreatnonzeros);

   return SCIP_OKAY;
}

/* ---------------- Constraint specific interface methods ---------------- */

/** creates the handler for cardinality constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrCardinality(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->eventhdlr = NULL;
   conshdlrdata->varhash = NULL;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
        eventExecCardinality, NULL) );
   if( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for cardinality constraints not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpCardinality, consEnfopsCardinality, consCheckCardinality, consLockCardinality, conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyCardinality, consCopyCardinality) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteCardinality) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolCardinality) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeCardinality) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsCardinality) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsCardinality) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpCardinality) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseCardinality) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolCardinality, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintCardinality) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropCardinality, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
        CONSHDLR_PROP_TIMING) );
   /*SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropCardinality) ); @todo: implement repropagation */
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpCardinality, consSepasolCardinality, CONSHDLR_SEPAFREQ,
        CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransCardinality) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxCardinality) );

   /* add cardinality constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/branchbalanced",
         "whether to use balanced instead of unbalanced branching",
         &conshdlrdata->branchbalanced, TRUE, DEFAULT_BRANCHBALANCED, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/balanceddepth",
         "maximum depth for using balanced branching (-1: no limit)",
         &conshdlrdata->balanceddepth, TRUE, DEFAULT_BALANCEDDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/balancedcutoff",
         "determines that balanced branching is only used if the branching cut off value \
         w.r.t. the current LP solution is greater than a given value",
         &conshdlrdata->balancedcutoff, TRUE, DEFAULT_BALANCEDCUTOFF, 0.01, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a cardinality constraint
 *
 *  We set the constraint to not be modifable. If the weights are non
 *  NULL, the variables are ordered according to these weights (in
 *  ascending order).
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method \ref SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   int                   cardval,            /**< number of variables allowed to be nonzero */
   SCIP_VAR**            indvars,            /**< indicator variables indicating which variables may be treated as nonzero
                                              *   in cardinality constraint, or NULL if new indicator variables should be
                                              *   introduced automatically */
   SCIP_Real*            weights,            /**< weights determining the variable order, or NULL if variables should be
                                                  ordered in the same way they were added to the constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_Bool modifiable;
   SCIP_Bool transformed;
   int v;

   modifiable = FALSE;

   /* find the cardinality constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("<%s> constraint handler not found\n", CONSHDLR_NAME);
      return SCIP_PLUGINNOTFOUND;
   }

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED;

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->cons = NULL;
   consdata->vars = NULL;
   consdata->indvars = NULL;
   consdata->eventdatas = NULL;
   consdata->nvars = nvars;
   consdata->cardval = cardval;
   consdata->maxvars = nvars;
   consdata->rowub = NULL;
   consdata->rowlb = NULL;
   consdata->eventdatascurrent = NULL;
   consdata->eventvarscurrent = NULL;
   consdata->neventdatascurrent = 0;
   consdata->ntreatnonzeros = transformed ? 0 : -1;
   consdata->weights = NULL;

   if( nvars > 0 )
   {
      /* duplicate memory for implied variables */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->vars, vars, nvars) );

      /* create indicator variables if not present */
      if( indvars != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->indvars, indvars, nvars) );
      }
      else
      {
         if( conshdlrdata->varhash == NULL )
         {
            /* set up hash map */
            SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->varhash, SCIPblkmem(scip), SCIPgetNTotalVars(scip)) );
         }

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->indvars, nvars) );
         for( v = 0; v < nvars; ++v )
         {
            SCIP_VAR* implvar;

            implvar = vars[v];
            assert(implvar != NULL);

            /* check whether an indicator variable already exists for implied variable */
            if( SCIPhashmapExists(conshdlrdata->varhash, implvar) )
            {
               assert((SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, implvar) != NULL);
               consdata->indvars[v] = (SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, implvar);
            }
            else
            {
               /* if implied variable is binary, then it is not necessary to create an indicator variable */
               if( SCIPvarIsBinary(implvar) )
                  consdata->indvars[v] = implvar;
               else
               {
                  char varname[SCIP_MAXSTRLEN];
                  SCIP_VAR* var;

                  (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ind_%s", SCIPvarGetName(vars[v]));
                  SCIP_CALL( SCIPcreateVar(scip, &var, varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, FALSE, FALSE,
                        NULL, NULL, NULL, NULL, NULL) );
                  SCIP_CALL( SCIPaddVar(scip, var) );
                  consdata->indvars[v] = var;
                  SCIP_CALL( SCIPreleaseVar(scip, &var) );
               }

               /* insert implied variable to hash map */
               SCIP_CALL( SCIPhashmapInsert(conshdlrdata->varhash, implvar, (void*) consdata->indvars[v]) );/*lint !e571*/
               assert(consdata->indvars[v] == (SCIP_VAR*) SCIPhashmapGetImage(conshdlrdata->varhash, implvar));
               assert(SCIPhashmapExists(conshdlrdata->varhash, implvar));
            }
         }
      }

      /* allocate block memory */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventdatascurrent, 4*nvars) );/*lint !e647*/
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventvarscurrent, 4*nvars) );/*lint !e647*/
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventdatas, nvars) );

      /* check weights */
      if( weights != NULL )
      {
         int* dummy;

         /* store weights */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->weights, weights, nvars) );

         /* create dummy array to make code compatible with SCIP 3.2.0
          * (the function SCIPsortRealPtrPtr() is not available) */
         SCIP_CALL( SCIPallocBufferArray(scip, &dummy, nvars) );
         for( v = 0; v < nvars; ++v )
            dummy[v] = 0;

         /* sort variables - ascending order */
         SCIPsortRealPtrPtrInt(consdata->weights, (void**)consdata->vars, (void**)consdata->indvars, dummy, nvars);

         SCIPfreeBufferArray(scip, &dummy);
      }
   }
   else
   {
      assert(weights == NULL);
   }

   /* create cardinality constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   consdata->cons = *cons;
   assert(consdata->cons != NULL);

   /* replace original variables by transformed variables in transformed constraint, add locks, and catch events */
   for( v = nvars - 1; v >= 0; --v )
   {
      /* always use transformed variables in transformed constraints */
      if( transformed )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, consdata->vars[v], &(consdata->vars[v])) );
         SCIP_CALL( SCIPgetTransformedVar(scip, consdata->indvars[v], &(consdata->indvars[v])) );
      }
      assert(consdata->vars[v] != NULL);
      assert(consdata->indvars[v] != NULL);
      assert(transformed == SCIPvarIsTransformed(consdata->vars[v]));
      assert(transformed == SCIPvarIsTransformed(consdata->indvars[v]));

      /* handle the new variable */
      SCIP_CALL( handleNewVariableCardinality(scip, *cons, consdata, conshdlrdata, consdata->vars[v],
           consdata->indvars[v], v, transformed, &consdata->eventdatas[v]) );
      assert(! transformed || consdata->eventdatas[v] != NULL);
   }

   return SCIP_OKAY;
}

/** creates and captures a cardinality constraint with all constraint flags set to their default values.
 *
 *  @warning Do NOT set the constraint to be modifiable manually, because this might lead
 *  to wrong results as the variable array will not be resorted
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method \ref SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   int                   cardval,            /**< number of variables allowed to be nonzero */
   SCIP_VAR**            indvars,            /**< indicator variables indicating which variables may be treated as nonzero
                                              *   in cardinality constraint, or NULL if new indicator variables should be
                                              *   introduced automatically */
   SCIP_Real*            weights             /**< weights determining the variable order, or NULL if variables should be
                                              *   ordered in the same way they were added to the constraint */
   )
{
   SCIP_CALL( SCIPcreateConsCardinality(scip, cons, name, nvars, vars, cardval, indvars, weights, TRUE, TRUE, TRUE, TRUE,
        TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** changes cardinality value of cardinality constraint (i.e., right hand side of cardinality constraint) */
SCIP_RETCODE  SCIPchgCardvalCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pointer to hold the created constraint */
   int                   cardval             /**< number of variables allowed to be nonzero */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cardinality constraint.\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMsg(scip, "modify right hand side of cardinality constraint from <%i> to <%i>\n", consdata->cardval, cardval);

   /* create constraint data */
   consdata->cardval = cardval;

   return SCIP_OKAY;
}

/** adds variable to cardinality constraint, the position is determined by the given weight */
SCIP_RETCODE SCIPaddVarCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_VAR*             indvar,             /**< indicator variable indicating whether variable may be treated as nonzero
                                              *   in cardinality constraint (or NULL if this variable should be created
                                              *   automatically) */
   SCIP_Real             weight              /**< weight determining position of variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);
   assert(var != NULL);
   assert(cons != NULL);

   SCIPdebugMsg(scip, "adding variable <%s> to constraint <%s> with weight %g\n", SCIPvarGetName(var),
        SCIPconsGetName(cons), weight);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cardinality constraint.\n");
      return SCIP_INVALIDDATA;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( addVarCardinality(scip, cons, conshdlrdata, var, indvar, weight) );

   return SCIP_OKAY;
}

/** appends variable to cardinality constraint */
SCIP_RETCODE SCIPappendVarCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to add to the constraint */
   SCIP_VAR*             indvar              /**< indicator variable indicating whether variable may be treated as nonzero
                                              *   in cardinality constraint (or NULL if this variable should be created
                                              *   automatically) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);
   assert(var != NULL);
   assert(cons != NULL);

   SCIPdebugMsg(scip, "appending variable <%s> to constraint <%s>\n", SCIPvarGetName(var), SCIPconsGetName(cons));

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cardinality constraint.\n");
      return SCIP_INVALIDDATA;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( appendVarCardinality(scip, cons, conshdlrdata, var, indvar) );

   return SCIP_OKAY;
}

/** gets number of variables in cardinality constraint */
int SCIPgetNVarsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cardinality constraint.\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in cardinality constraint */
SCIP_VAR** SCIPgetVarsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cardinality constraint.\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets cardinality value of cardinality constraint (i.e., right hand side of cardinality constraint) */
int SCIPgetCardvalCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cardinality constraint.\n");
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->cardval;
}

/** gets array of weights in cardinality constraint (or NULL if not existent) */
SCIP_Real* SCIPgetWeightsCardinality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a cardinality constraint.\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->weights;
}
