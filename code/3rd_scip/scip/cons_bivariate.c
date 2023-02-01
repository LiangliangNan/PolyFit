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

/**@file   cons_bivariate.c
 * @brief  constraint handler for bivariate nonlinear constraints \f$\textrm{lhs} \leq f(x,y) + c z \leq \textrm{rhs}\f$
 * @author Martin Ballerstein
 * @author Dennis Michaels
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>

#include "scip/cons_bivariate.h"
#include "scip/cons_linear.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/debug.h"
#include "nlpi/nlpi.h"
#include "nlpi/exprinterpret.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "bivariate"
#define CONSHDLR_DESC          "constraint handler for constraints of the form lhs <= f(x,y) + c*z <= rhs where f(x,y) is a bivariate function"
#define CONSHDLR_SEPAPRIORITY         5 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -55 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -3600000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_FAST
#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define INTERVALINFTY             1E+43 /**< value for infinity in interval operations */
#define NEWTONMAXITER              1000 /**< maximal number of iterations in newton method */
#define INITLPMAXVARVAL          1000.0 /**< maximal absolute value of variable for still generating a linearization cut at that point in initlp */

#define QUADCONSUPGD_PRIORITY      5000 /**< priority of the constraint handler for upgrading of quadratic constraints */
#define NONLINCONSUPGD_PRIORITY   10000 /**< priority of the constraint handler for upgrading of nonlinear constraints */

/* activate the following define to get output on number of bivariate constraints for each convexity-type during INITSOL */
/* #define TYPESTATISTICS */

/*
 * Data structures
 */

/** data structure to cache data used for separation of convex-concave constraints */
struct SepaData_ConvexConcave
{
   SCIP_Bool             linearinx;          /**< whether the function is linear in x */
   SCIP_Bool             lineariny;          /**< whether the function is linear in y */
   SCIP_EXPRTREE*        f_yfixed;           /**< expression tree for f(x,yfixed) */
   SCIP_EXPRTREE*        f_neg_swapped;      /**< expression tree for -f(y,x) */
   SCIP_EXPRTREE*        f_neg_swapped_yfixed;/**< expression tree for -f(y,xfixed) */
   SCIP_EXPRTREE*        vred;               /**< expression tree for vred to underestimate  f(x,y) */
   SCIP_EXPRTREE*        vred_neg_swapped;   /**< expression tree for vred to underestimate -f(y,x) */
};
/** data structure to cache data used for separation of convex-concave constraints */
typedef struct SepaData_ConvexConcave SEPADATA_CONVEXCONCAVE;

/** constraint data for bivariate constraints */
struct SCIP_ConsData
{
   SCIP_EXPRTREE*        f;                  /**< expression tree of bivariate function f(x,y) */
   SCIP_BIVAR_CONVEXITY  convextype;         /**< kind of convexity of f(x,y) */
   SCIP_VAR*             z;                  /**< linear variable */
   SCIP_Real             zcoef;              /**< coefficient of linear variable */
   SCIP_Real             lhs;                /**< left hand side */
   SCIP_Real             rhs;                /**< right hand side */

   SCIP_Real             activity;           /**< activity of bivariate function w.r.t. current solution */
   SCIP_Real             lhsviol;            /**< violation of left hand side in current solution */
   SCIP_Real             rhsviol;            /**< violation of left hand side in current solution */

   unsigned int          mayincreasez:1;     /**< whether z can be increased without harming other constraints */
   unsigned int          maydecreasez:1;     /**< whether z can be decreased without harming other constraints */
   int                   eventfilterpos;     /**< position of z var events in SCIP event filter */

   SCIP_EXPRGRAPHNODE*   exprgraphnode;      /**< node in expression graph corresponding to bivariate function */

   SEPADATA_CONVEXCONCAVE sepaconvexconcave; /**< separation data for convex-concave constraints */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EXPRINT*         exprinterpreter;    /**< expression interpreter (computer gradients and hessians) */

   SCIP_Real             cutmaxrange;        /**< maximal range (maximal coef / minimal coef) of a cut in order to be added to LP */
   SCIP_Bool             linfeasshift;       /**< whether to make solutions in check feasible if possible */
   int                   maxproprounds;      /**< limit on number of propagation rounds for a single constraint within one round of SCIP propagation */
   int                   ninitlprefpoints;   /**< number of reference points in each direction where to compute linear support for envelope in LP initialization */
   SCIP_Bool             enfocutsremovable;  /**< are cuts added during enforcement removable from the LP in the same node? */

   SCIP_EVENTHDLR*       linvareventhdlr;    /**< handler for linear variable bound change events */
   SCIP_EVENTHDLR*       nonlinvareventhdlr; /**< handler for nonlinear variable bound change events */
   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subNLP heuristic */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the TRYSOL heuristic, if available */
   int                   newsoleventfilterpos;/**< filter position of new solution event handler, if catched */

   SCIP_EXPRGRAPH*       exprgraph;          /**< expression graph */
   SCIP_Bool             isremovedfixings;   /**< whether variable fixations have been removed from the expression graph */
   SCIP_Bool             ispropagated;       /**< whether the bounds on the variables in the expression graph have been propagated */
   SCIP*                 scip;               /**< SCIP data structure, needed in expression graph callbacks */

   SCIP_NODE*            lastenfonode;       /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenforounds;        /**< counter on number of enforcement rounds for the current node */
};


/*
 * Local methods
 */

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/** processes bound tightening event */
static
SCIP_DECL_EVENTEXEC(processLinearVarEvent)
{
   SCIP_CONS* cons;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_BOUNDTIGHTENED);

   cons = (SCIP_CONS*) eventdata;
   assert(cons != NULL);

   SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );

   return SCIP_OKAY;
}

/** catches variable bound change events on the linear variable in a bivariate constraint */
static
SCIP_RETCODE catchLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsEnabled(cons));
   assert(SCIPconsIsTransformed(cons));

   assert(SCIPconsGetHdlr(cons) != NULL);
   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->linvareventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->z == NULL )
      return SCIP_OKAY;
   assert(consdata->eventfilterpos == -1);

   eventtype = SCIP_EVENTTYPE_DISABLED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest */
      if( consdata->zcoef > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;
      else
         eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest */
      if( consdata->zcoef > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
      else
         eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;
   }

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->z, eventtype, conshdlrdata->linvareventhdlr, (SCIP_EVENTDATA*)cons, &consdata->eventfilterpos) );

   SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );

   return SCIP_OKAY;
}

/** drops variable bound change events on the linear variable in a bivariate constraint */
static
SCIP_RETCODE dropLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   assert(SCIPconsGetHdlr(cons) != NULL);
   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->linvareventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->z == NULL )
      return SCIP_OKAY;
   assert(consdata->eventfilterpos >= 0);

   eventtype = SCIP_EVENTTYPE_DISABLED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest */
      if( consdata->zcoef > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;
      else
         eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest */
      if( consdata->zcoef > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
      else
         eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;
   }

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->z, eventtype, conshdlrdata->linvareventhdlr, (SCIP_EVENTDATA*)cons, consdata->eventfilterpos) );
   consdata->eventfilterpos = -1;

   return SCIP_OKAY;
}


/** processes bound change events for variables in expression graph */
static
SCIP_DECL_EVENTEXEC(processNonlinearVarEvent)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   conshdlrdata = (SCIP_CONSHDLRDATA*)SCIPeventhdlrGetData(eventhdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   eventtype = SCIPeventGetType(event);
   assert( eventtype & (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED) );

   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      SCIPdebugMsg(scip, "changed %s bound on expression graph variable <%s> from %g to %g\n",
         (eventtype & SCIP_EVENTTYPE_LBCHANGED) ? "lower" : "upper",
         SCIPvarGetName(SCIPeventGetVar(event)), SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));

      if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
         conshdlrdata->ispropagated = FALSE;

      /* update variable bound in expression graph
       * @todo should we add epsilon to variable range?
       */
      if( eventtype & SCIP_EVENTTYPE_LBCHANGED )
         SCIPexprgraphSetVarNodeLb(conshdlrdata->exprgraph, (SCIP_EXPRGRAPHNODE*)eventdata,
            -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -SCIPeventGetNewbound(event)));  /*lint !e666*/
      else
         SCIPexprgraphSetVarNodeUb(conshdlrdata->exprgraph, (SCIP_EXPRGRAPHNODE*)eventdata,
            +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  SCIPeventGetNewbound(event)));  /*lint !e666*/
   }
   else
   {
      assert(eventtype & SCIP_EVENTTYPE_VARFIXED);
      conshdlrdata->isremovedfixings = FALSE;
   }

   return SCIP_OKAY;
}

/** callback method for variable addition in expression graph */
static
SCIP_DECL_EXPRGRAPHVARADDED( exprgraphVarAdded )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTERVAL varbounds;
   SCIP_VAR* var_;

   assert(exprgraph != NULL);
   assert(var != NULL);
   assert(varnode != NULL);

   var_ = (SCIP_VAR*)var;

   conshdlrdata = (SCIP_CONSHDLRDATA*)userdata;
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph == exprgraph);

   /* catch variable bound change events */
   SCIP_CALL( SCIPcatchVarEvent(conshdlrdata->scip, (SCIP_VAR*)var, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, conshdlrdata->nonlinvareventhdlr, (SCIP_EVENTDATA*)varnode, NULL) );
   SCIPdebugMessage("catch boundchange events on new expression graph variable <%s>\n", SCIPvarGetName(var_));

   /* set current bounds in expression graph */
   SCIPintervalSetBounds(&varbounds,
      -infty2infty(SCIPinfinity(conshdlrdata->scip), INTERVALINFTY, -MIN(SCIPvarGetLbLocal(var_), SCIPvarGetUbLocal(var_))),  /*lint !e666*/
      +infty2infty(SCIPinfinity(conshdlrdata->scip), INTERVALINFTY,  MAX(SCIPvarGetLbLocal(var_), SCIPvarGetUbLocal(var_)))   /*lint !e666*/
      );
   SCIPexprgraphSetVarNodeBounds(exprgraph, varnode, varbounds);

   SCIP_CALL( SCIPaddVarLocks(conshdlrdata->scip, var_, 1, 1) );
   SCIPdebugMessage("increased up- and downlocks of variable <%s>\n", SCIPvarGetName(var_));

   conshdlrdata->isremovedfixings &= SCIPvarIsActive(var_);
   conshdlrdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** callback method for variable removal in expression graph */
static
SCIP_DECL_EXPRGRAPHVARREMOVE( exprgraphVarRemove )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR* var_;

   assert(exprgraph != NULL);
   assert(var != NULL);
   assert(varnode != NULL);

   var_ = (SCIP_VAR*)var;

   conshdlrdata = (SCIP_CONSHDLRDATA*)userdata;
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph == exprgraph);

   SCIP_CALL( SCIPdropVarEvent(conshdlrdata->scip, var_, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, conshdlrdata->nonlinvareventhdlr, (SCIP_EVENTDATA*)varnode, -1) );
   SCIPdebugMessage("drop boundchange events on expression graph variable <%s>\n", SCIPvarGetName(var_));

   SCIP_CALL( SCIPaddVarLocks(conshdlrdata->scip, var_, -1, -1) );
   SCIPdebugMessage("decreased up- and downlocks of variable <%s>\n", SCIPvarGetName(var_));

   return SCIP_OKAY;
}

/** locks linear variable in a constraint */
static
SCIP_RETCODE lockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to lock a variable */
   SCIP_VAR*             var,                /**< variable to lock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** unlocks linear variable in a constraint */
static
SCIP_RETCODE unlockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to unlock a variable */
   SCIP_VAR*             var,                /**< variable to unlock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** resolves variable fixations and aggregations in a constraint */
static
SCIP_RETCODE removeFixedVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint where to remove fixed variables */
   SCIP_Bool*            ischanged,          /**< buffer to store whether something was changed in the constraint */
   SCIP_Bool*            isupgraded          /**< buffer to store whether the constraint has been upgraded (and deleted) */
   )
{
#ifndef NDEBUG
   SCIP_CONSHDLRDATA* conshdlrdata;
#endif
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* substexpr[2];
   SCIP_VAR* var;
   SCIP_VAR* vars[2];
   SCIP_Real coef;
   SCIP_Real constant;
   int i;

   assert(conshdlr != NULL);
   assert(scip != NULL);
   assert(cons != NULL);
   assert(ischanged != NULL);
   assert(isupgraded != NULL);

#ifndef NDEBUG
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
#endif

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->f != NULL);

   *ischanged = FALSE;
   *isupgraded = FALSE;

   if( consdata->z != NULL && !SCIPvarIsActive(consdata->z) && SCIPvarGetStatus(consdata->z) != SCIP_VARSTATUS_MULTAGGR )
   {
      /* replace z by active or multaggr. variable */

      /* drop events on z, unlock and release variable */
      SCIP_CALL( dropLinearVarEvents(scip, cons) );
      SCIP_CALL( unlockLinearVariable(scip, cons, consdata->z, consdata->zcoef) );

      /* replace by new variable, or NULL */
      constant = 0.0;
      SCIP_CALL( SCIPgetProbvarSum(scip, &consdata->z, &consdata->zcoef, &constant) );
      if( consdata->zcoef == 0.0 )
         consdata->z = NULL;
      if( constant != 0.0 && !SCIPisInfinity(scip, -consdata->lhs) )
         consdata->lhs -= constant;
      if( constant != 0.0 && !SCIPisInfinity(scip,  consdata->rhs) )
         consdata->rhs -= constant;

      if( consdata->z != NULL )
      {
         /* catch events on new z, lock and capture variable, mark as not to multaggr */
         SCIP_CALL( catchLinearVarEvents(scip, cons) );
         SCIP_CALL( lockLinearVariable(scip, cons, consdata->z, consdata->zcoef) );
         if( SCIPvarIsActive(consdata->z) )
         {
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, consdata->z) );
         }
      }

      *ischanged = TRUE;
   }

   assert(SCIPexprtreeGetNVars(consdata->f) == 2);
   vars[0] = SCIPexprtreeGetVars(consdata->f)[0];
   vars[1] = SCIPexprtreeGetVars(consdata->f)[1];

   if( vars[0] == NULL || vars[1] == NULL )
      return SCIP_INVALIDDATA;

   if( SCIPvarGetStatus(SCIPvarGetProbvar(vars[0])) == SCIP_VARSTATUS_FIXED ||
      SCIPvarGetStatus(SCIPvarGetProbvar(vars[1])) == SCIP_VARSTATUS_FIXED ||
      SCIPvarGetProbvar(vars[0]) == SCIPvarGetProbvar(vars[1]) )
   {
      /* if number of variable reduces, then upgrade to nonlinear constraint
       * except if we are in the exit-presolving stage, where upgrading is not allowed
       * in the latter case, we just do nothing, which may not be most efficient, but should still work
       */
      SCIP_EXPRTREE* tree;
      SCIP_CONS* nlcons;

      if( SCIPgetStage(scip) == SCIP_STAGE_EXITPRESOLVE )
         return SCIP_OKAY;

      SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &tree, consdata->f) );

      for( i = 0; i < 2; ++i )
      {
         substexpr[i] = NULL;

         var = vars[i];
         if( (SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR) )
            continue;

         coef = 1.0;
         constant = 0.0;
         SCIP_CALL( SCIPgetProbvarSum(scip, &var, &coef, &constant) );

         if( coef == 0.0 )
         {
            /* replace var_i by constant in expression tree */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &substexpr[i], SCIP_EXPR_CONST, constant) );
            vars[i] = NULL;
         }
         else if( coef == 1.0 && constant == 0.0 )
         {
            /* do not need to change expression tree, just store new variable in tree */
            substexpr[i] = NULL;
            vars[i] = var;
         }
         else
         {
            /* replace var_i by coef * var_i + constant in expression tree */
            SCIP_EXPR* child;

            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &child, SCIP_EXPR_VARIDX, i) );
            SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &substexpr[i], 1, &child, &coef, constant) );
            vars[i] = var;
         }
      }

      assert(substexpr[0] != NULL || substexpr[1] != NULL);

      SCIP_CALL( SCIPexprtreeSubstituteVars(tree, substexpr) );
      if( substexpr[0] != NULL )
         SCIPexprFreeDeep(SCIPblkmem(scip), &substexpr[0]);
      if( substexpr[1] != NULL )
         SCIPexprFreeDeep(SCIPblkmem(scip), &substexpr[1]);

      /* if variable 0 has been remove or is the same as variable 1, reindex 1 to 0 */
      if( (vars[0] == NULL || vars[0] == vars[1]) && vars[1] != NULL )
      {
         int reindex[2];

         reindex[0] = 0;
         reindex[1] = 0;
         SCIPexprReindexVars(SCIPexprtreeGetRoot(tree), reindex);
         vars[0] = vars[1];
         vars[1] = NULL;
      }

      /* update variables array in tree */
      assert(vars[1] == NULL || vars[0] != NULL);
      SCIP_CALL( SCIPexprtreeSetVars(tree, vars[0] == NULL ? 0 : (vars[1] == NULL ? 1 : 2), vars) );

      SCIP_CALL( SCIPcreateConsNonlinear(scip, &nlcons, SCIPconsGetName(cons),
            consdata->z != NULL ? 1 : 0, consdata->z != NULL ? &consdata->z : NULL, &consdata->zcoef,
            1, &tree, NULL, consdata->lhs, consdata->rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );  /*lint !e826*/
      SCIP_CALL( SCIPaddCons(scip, nlcons) );
      SCIPdebugMsg(scip, "upgraded to"); SCIPdebugPrintCons(scip, nlcons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &nlcons) );

      *isupgraded = TRUE;

      SCIP_CALL( SCIPexprtreeFree(&tree) );

      return SCIP_OKAY;
   }

   for( i = 0; i < 2; ++i )
   {
      substexpr[i] = NULL;

      var = vars[i];
      if( SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
         continue;

      coef = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPgetProbvarSum(scip, &var, &coef, &constant) );
      assert(coef != 0.0); /* fixed vars should have been handled above */

      if( coef == 1.0 && constant == 0.0 )
      {
         /* do not need to change expression tree, just store new variable in tree */
         substexpr[i] = NULL;
         vars[i] = var;
      }
      else
      {
         /* replace var_i by coef * var_i + constant in expression tree */
         SCIP_EXPR* child;

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &child, SCIP_EXPR_VARIDX, i) );
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &substexpr[i], 1, &child, &coef, constant) );
         vars[i] = var;
      }

      /* update variables array in tree for next operation */
      SCIP_CALL( SCIPexprtreeSetVars(consdata->f, 2, vars) );

      /* mark that variables in constraint should not be multiaggregated (bad for bound tightening and branching) */
      if( SCIPvarIsActive(vars[0]) )
      {
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, vars[0]) );
      }
      if( SCIPvarIsActive(vars[1]) )
      {
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, vars[1]) );
      }

      *ischanged = TRUE;
   }

   /* update expression tree, if necessary */
   if( substexpr[0] != NULL || substexpr[1] != NULL )
   {
      SCIP_CALL( SCIPexprtreeSubstituteVars(consdata->f, substexpr) );
      if( substexpr[0] != NULL )
         SCIPexprFreeDeep(SCIPblkmem(scip), &substexpr[0]);
      if( substexpr[1] != NULL )
         SCIPexprFreeDeep(SCIPblkmem(scip), &substexpr[1]);
   }

   return SCIP_OKAY;
}

/** removes fixed variables from expression graph */
static
SCIP_RETCODE removeFixedNonlinearVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int nvars;
   int varssize;
   SCIP_Real constant;
   int i;
   int requsize;
   SCIPdebug( int j );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   if( conshdlrdata->isremovedfixings )
      return SCIP_OKAY;

   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, varssize) );

   i = 0;
   while( i < SCIPexprgraphGetNVars(conshdlrdata->exprgraph) )
   {
      var = (SCIP_VAR*) SCIPexprgraphGetVars(conshdlrdata->exprgraph)[i];
      if( SCIPvarIsActive(var) )
      {
         ++i;
         continue;
      }

      vars[0]  = var;
      coefs[0] = 1.0;
      constant = 0.0;
      nvars = 1;
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, coefs, &nvars, varssize, &constant, &requsize, TRUE) );

      if( requsize > varssize )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars,  requsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, requsize) );
         varssize = requsize;
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, coefs, &nvars, varssize, &constant, &requsize, TRUE) );
         assert(requsize <= varssize);
      }

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "replace fixed variable <%s> by %g", SCIPvarGetName(var), constant);
      for( j = 0; j < nvars; ++j )
      {
         SCIPdebugMsgPrint(scip, " %+g <%s>", coefs[j], SCIPvarGetName(vars[j]));
      }
      SCIPdebugMsgPrint(scip, "\n");
#endif

      SCIP_CALL( SCIPexprgraphReplaceVarByLinearSum(conshdlrdata->exprgraph, var, nvars, coefs, (void**)vars, constant) );

      i = 0;
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &coefs);

   conshdlrdata->isremovedfixings = TRUE;

   return SCIP_OKAY;
}

/** computes violation of a constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol                 /**< solution or NULL if LP solution should be used */
   )
{  /*lint --e{666}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real xyvals[2];
   SCIP_Real zval = 0.0;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real ylb;
   SCIP_Real yub;
   SCIP_Real absviol;
   SCIP_Real relviol;
   SCIP_VAR* x;
   SCIP_VAR* y;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprinterpreter != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPexprtreeGetInterpreterData(consdata->f) == NULL )
   {
      SCIP_CALL( SCIPexprintCompile(conshdlrdata->exprinterpreter, consdata->f) );
   }

   x = SCIPexprtreeGetVars(consdata->f)[0];
   y = SCIPexprtreeGetVars(consdata->f)[1];

   xyvals[0] = SCIPgetSolVal(scip, sol, x);
   xyvals[1] = SCIPgetSolVal(scip, sol, y);
   if( consdata->z != NULL )
      zval = SCIPgetSolVal(scip, sol, consdata->z);

   /* @todo proper handling of variables at infinity
    * for now, just say infeasible and keep fingers crossed
    */
   if( SCIPisInfinity(scip, REALABS(xyvals[0])) )
   {
      consdata->lhsviol = consdata->rhsviol = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   if( SCIPisInfinity(scip, REALABS(xyvals[1])) )
   {
      consdata->lhsviol = consdata->rhsviol = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   /* project point onto box if from LP or very close to bounds to avoid eval error when function is not defined slightly outside bounds */
   xlb = SCIPvarGetLbGlobal(x);
   xub = SCIPvarGetUbGlobal(x);
   ylb = SCIPvarGetLbGlobal(y);
   yub = SCIPvarGetUbGlobal(y);
   /* @todo handle case where variables are outside of bounds as in other constraint handlers, see also #627 */
   if( sol == NULL )
   {
      assert(SCIPisFeasGE(scip, xyvals[0], xlb));
      assert(SCIPisFeasLE(scip, xyvals[0], xub));
      xyvals[0] = MAX(xlb, MIN(xub, xyvals[0]));

      assert(SCIPisFeasGE(scip, xyvals[1], ylb));
      assert(SCIPisFeasLE(scip, xyvals[1], yub));
      xyvals[1] = MAX(ylb, MIN(yub, xyvals[1]));

      if( consdata->z != NULL )
      {
         assert(SCIPisFeasGE(scip, zval, SCIPvarGetLbLocal(consdata->z)));
         assert(SCIPisFeasLE(scip, zval, SCIPvarGetUbLocal(consdata->z)));
         zval = MAX(SCIPvarGetLbLocal(consdata->z), MIN(SCIPvarGetUbLocal(consdata->z), zval));
      }
   }
   else
   {
      if( SCIPisEQ(scip, xyvals[0], xlb) || SCIPisEQ(scip, xyvals[0], xub) )
         xyvals[0] = MAX(xlb, MIN(xub, xyvals[0]));
      if( SCIPisEQ(scip, xyvals[1], ylb) || SCIPisEQ(scip, xyvals[1], yub) )
         xyvals[1] = MAX(ylb, MIN(yub, xyvals[1]));
   }

   /* compute activity of constraint */
   SCIP_CALL( SCIPexprintEval(conshdlrdata->exprinterpreter, consdata->f, xyvals, &consdata->activity) );

   /* point is outside the domain of f */
   if( !SCIPisFinite(consdata->activity) )
   {
       consdata->lhsviol = consdata->rhsviol = SCIPinfinity(scip);
       return SCIP_OKAY;
   }

   if( consdata->z != NULL )
      consdata->activity += consdata->zcoef * zval;

   /* compute violation of constraint sides */
   absviol = 0.0;
   relviol = 0.0;
   if( consdata->activity < consdata->lhs && !SCIPisInfinity(scip, -consdata->lhs) )
   {
      consdata->lhsviol = consdata->lhs - consdata->activity;
      absviol = consdata->lhsviol;
      relviol = SCIPrelDiff(consdata->lhs, consdata->activity);
   }
   else
      consdata->lhsviol = 0.0;

   if( consdata->activity > consdata->rhs && !SCIPisInfinity(scip,  consdata->rhs) )
   {
      consdata->rhsviol = consdata->activity - consdata->rhs;
      absviol = consdata->rhsviol;
      relviol = SCIPrelDiff(consdata->activity, consdata->rhs);
   }
   else
      consdata->rhsviol = 0.0;

   if( sol != NULL )
      SCIPupdateSolConsViolation(scip, sol, absviol, relviol);

   return SCIP_OKAY;
}

/** computes violation of a set of constraints */
static
SCIP_RETCODE computeViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_CONS**           maxviolcon          /**< buffer to store constraint with largest violation, or NULL if solution is feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      viol;
   SCIP_Real      maxviol;
   int            c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(maxviolcon != NULL);

   *maxviolcon = NULL;

   maxviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      SCIP_CALL( computeViolation(scip, conshdlr, conss[c], sol) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      viol = MAX(consdata->lhsviol, consdata->rhsviol);
      if( viol > maxviol && SCIPisGT(scip, viol, SCIPfeastol(scip)) )
      {
         maxviol = viol;
         *maxviolcon = conss[c];
      }
   }

   return SCIP_OKAY;
}

/** setup vred(s;x0,y0,ylb,yub) for a given f(x,y) for computing a convex-concave underestimator
 * vred(s;x0,y0,ylb,yub) = (yub-y0)/(yub-ylb) f((yub-ylb)/(yub-y0)x0 - (y0-ylb)/(yub-y0)*s, ylb) + (y0-ylb)/(yub-ylb) f(s,yub)
 */
static
SCIP_RETCODE initSepaDataCreateVred(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE**       vred,               /**< buffer where to store exprtree for vred */
   SCIP_EXPRTREE*        f                   /**< function f(x,y) for which vred should be setup */
   )
{
   SCIP_EXPR* subst[2];
   SCIP_Real minusone;
   SCIP_EXPR* e1;
   SCIP_EXPR* e2;
   SCIP_EXPR* e3;
   SCIP_EXPR* e4;
   SCIP_EXPR* e5;
   SCIP_EXPR* e6;
   SCIP_EXPR* arg1;
   SCIP_EXPR* arg2;
   SCIP_EXPR* vredexpr;

   assert(scip != NULL);
   assert(vred != NULL);
   assert(f != NULL);
   assert(SCIPexprGetOperator(SCIPexprtreeGetRoot(f)) != SCIP_EXPR_VARIDX);  /* substitute cannot substitute the root node, but f should not be a single variable anyway */

   /* setup vred(s;x0,y0,ylb,yub) for computing a convex-concave underestimator in the case where y is not at one of its bounds
    * vred(s;x0,y0,ylb,yub) = (yub-y0)/(yub-ylb) f((yub-ylb)/(yub-y0)x0 - (y0-ylb)/(yub-y0)*s, ylb) + (y0-ylb)/(yub-ylb) f(s,yub)
    */
   /* create expression for x0(yub-ylb)/(yub-y0) */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_PARAM, 2) ); /* ylb */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e2, SCIP_EXPR_PARAM, 3) ); /* yub */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e3, SCIP_EXPR_MINUS, e2, e1) ); /* yub-ylb */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_PARAM, 0) ); /* x0 */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e3, SCIP_EXPR_MUL, e1, e3) ); /* x0(yub-ylb) */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_PARAM, 1) ); /* y0 */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e2, SCIP_EXPR_PARAM, 3) ); /* yub */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e4, SCIP_EXPR_MINUS, e2, e1) ); /* yub-y0 */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e5, SCIP_EXPR_DIV, e3, e4) ); /* x0(yub-ylb)/(yub-y0) */

   /* create expression for s(y0-ylb)/(yub-y0) */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_PARAM, 1) ); /* y0 */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e2, SCIP_EXPR_PARAM, 2) ); /* ylb */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e3, SCIP_EXPR_MINUS, e1, e2) ); /* y0-ylb */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_VARIDX, 0) ); /* s */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e3, SCIP_EXPR_MUL, e1, e3) ); /* s(y0-ylb) */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_PARAM, 1) ); /* y0 */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e2, SCIP_EXPR_PARAM, 3) ); /* yub */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e4, SCIP_EXPR_MINUS, e2, e1) ); /* yub-y0 */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e6, SCIP_EXPR_DIV, e3, e4) ); /* s(y0-ylb)/(yub-y0) */

   /* create expression for (yub-ylb)/(yub-y0)x0 - (y0-ylb)/(yub-y0)*s */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[0], SCIP_EXPR_MINUS, e5, e6) );

   /* create expression for ylb */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_PARAM, 2) );

   /* create expression for f((yub-ylb)/(yub-y0)x0 - (y0-ylb)/(yub-y0)*s, ylb) */
   SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &arg1, SCIPexprtreeGetRoot(f)) );
   SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), arg1, subst) );
   SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
   SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);

   /* create expression for f(s,yub) */
   SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &arg2, SCIPexprtreeGetRoot(f)) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_PARAM, 3) );
   SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), arg2, subst) );
   SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);

   /* create expression for (yub-y0)/(yub-ylb) */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_PARAM, 1) ); /* y0 */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e2, SCIP_EXPR_PARAM, 3) ); /* yub */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e3, SCIP_EXPR_MINUS, e2, e1) ); /* yub-y0 */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_PARAM, 2) ); /* ylb */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e2, SCIP_EXPR_PARAM, 3) ); /* yub */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e4, SCIP_EXPR_MINUS, e2, e1) ); /* yub-ylb */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e5, SCIP_EXPR_DIV, e3, e4) ); /* (yub-y0)/(yub-ylb) */

   /* create expression for 1 - (yub-y0)/(yub-ylb) */
   minusone = -1.0;
   SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e1, e5) ); /* (yub-y0)/(yub-ylb) */
   SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &e6, 1, &e1, &minusone, 1.0) ); /* 1 - (yub-y0)/(yub-ylb) */

   /* create expression for vred = e5*arg1 + e6*arg2 */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_MUL, e5, arg1) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e2, SCIP_EXPR_MUL, e6, arg2) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &vredexpr, SCIP_EXPR_PLUS, e1, e2) );

   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), vred, vredexpr, 1, 4, NULL) );

   return SCIP_OKAY;
}

/** initializes separation data */
static
SCIP_RETCODE initSepaData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->f != NULL);

   switch( consdata->convextype )
   {
   case SCIP_BIVAR_CONVEX_CONCAVE:
   {
      SCIP_VAR** xy;
      SCIP_Real ref[2];
      SCIP_Bool sparsity[4];

      if( SCIPexprtreeGetInterpreterData(consdata->f) == NULL )
      {
         SCIP_CALL( SCIPexprintCompile(exprinterpreter, consdata->f) );
      }

      xy = SCIPexprtreeGetVars(consdata->f);
      assert(xy != NULL);

      /* check if the function is linear in x or y */
      ref[0] = MIN(MAX(SCIPvarGetLbLocal(xy[0]), 0.0), SCIPvarGetUbLocal(xy[0]));  /*lint !e666*/
      ref[1] = MIN(MAX(SCIPvarGetLbLocal(xy[1]), 0.0), SCIPvarGetUbLocal(xy[1]));  /*lint !e666*/

      SCIP_CALL( SCIPexprintHessianSparsityDense(exprinterpreter, consdata->f, ref, sparsity) );

      consdata->sepaconvexconcave.linearinx = !sparsity[0];
      consdata->sepaconvexconcave.lineariny = !sparsity[3];

      if( !consdata->sepaconvexconcave.linearinx && !SCIPisInfinity(scip,  consdata->rhs) )
      {
         SCIP_EXPR* subst[2];
         SCIP_Real one;

         /* setup f(x,yfixed) for computing a convex-concave underestimator in the case where y is at one of its bounds */
         SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &consdata->sepaconvexconcave.f_yfixed, consdata->f) );

         /* x stays x, nothing to substitute
          * y is substituted by SCIP_EXPR_PARAM
          */
         subst[0] = NULL;
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_PARAM, 0) );

         /* make y a parameter */
         SCIP_CALL( SCIPexprtreeSubstituteVars(consdata->sepaconvexconcave.f_yfixed, subst) );

         /* reset variables array to {x} and parameters array to {y} */
         one = 1.0;
         SCIP_CALL( SCIPexprtreeSetVars(consdata->sepaconvexconcave.f_yfixed, 1, &xy[0]) );
         SCIP_CALL( SCIPexprtreeSetParams(consdata->sepaconvexconcave.f_yfixed, 1, &one) );

         /* free subst[1] */
         SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);

         SCIP_CALL( SCIPexprintCompile(exprinterpreter, consdata->sepaconvexconcave.f_yfixed) );

         /* setup vred(s;x0,y0,ylb,yub) for computing a convex-concave underestimator in the case where y is not at one of its bounds
          * vred(s;x0,y0,ylb,yub) = (yub-y0)/(yub-ylb) f((yub-ylb)/(yub-y0)x0 - (y0-ylb)/(yub-y0)*s, ylb) + (y0-ylb)/(yub-ylb) f(s,yub)
          */
         SCIP_CALL( initSepaDataCreateVred(scip, &consdata->sepaconvexconcave.vred, consdata->f) );
         SCIP_CALL( SCIPexprintCompile(exprinterpreter, consdata->sepaconvexconcave.vred) );
      }
      else
      {
         consdata->sepaconvexconcave.f_yfixed = NULL;
         consdata->sepaconvexconcave.vred = NULL;
      }

      if( !consdata->sepaconvexconcave.lineariny && !SCIPisInfinity(scip, -consdata->lhs) )
      {
         /* if we have a left hand side and are not linear y in, then we may need to call
          * generateConvexConcaveUnderestimator for -f with swapped variables
          */
         SCIP_EXPR* minusf;
         SCIP_EXPR* fcopy;
         SCIP_VAR*  vars[2];
         int        reindex[2];
         SCIP_Real  minusone;
         SCIP_Real  one;
         SCIP_EXPR* subst[2];

         /* create expression for -f */
         minusone = -1.0;
         SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &fcopy, SCIPexprtreeGetRoot(consdata->f)) );
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &minusf, 1, &fcopy, &minusone, 0.0) );

         /* reindex/swap variables */
         reindex[0] = 1;
         reindex[1] = 0;
         SCIPexprReindexVars(minusf, reindex);

         /* create expression tree for -f(y,x) */
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &consdata->sepaconvexconcave.f_neg_swapped, minusf, 2, 0, NULL) );

         vars[0] = xy[1];
         vars[1] = xy[0];
         SCIP_CALL( SCIPexprtreeSetVars(consdata->sepaconvexconcave.f_neg_swapped, 2, vars) );

         SCIP_CALL( SCIPexprintCompile(exprinterpreter, consdata->sepaconvexconcave.f_neg_swapped) );

         /* setup -f(y, xfixed) for computing a convex-concave overestimator in the case where x is at on of it's bounds */
         SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &consdata->sepaconvexconcave.f_neg_swapped_yfixed, consdata->sepaconvexconcave.f_neg_swapped) );

         /* y stays y, nothing to substitute
          * x is substituted by SCIP_EXPR_PARAM
          */
         subst[0] = NULL;
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_PARAM, 0) );

         /* make x a parameter */
         SCIP_CALL( SCIPexprtreeSubstituteVars(consdata->sepaconvexconcave.f_neg_swapped_yfixed, subst) );

         /* reset variables array to {y} and parameters array to {x} */
         one = 1.0;
         SCIP_CALL( SCIPexprtreeSetVars(consdata->sepaconvexconcave.f_neg_swapped_yfixed, 1, &xy[1]) );
         SCIP_CALL( SCIPexprtreeSetParams(consdata->sepaconvexconcave.f_neg_swapped_yfixed, 1, &one) );

         /* free subst[1] */
         SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);

         SCIP_CALL( SCIPexprintCompile(exprinterpreter, consdata->sepaconvexconcave.f_neg_swapped_yfixed) );

         /* setup vred(s;y0,x0,xlb,xub) for computing a convex-concave underestimator in the case where x is not at one of its bounds */
         SCIP_CALL( initSepaDataCreateVred(scip, &consdata->sepaconvexconcave.vred_neg_swapped, consdata->sepaconvexconcave.f_neg_swapped) );
         SCIP_CALL( SCIPexprintCompile(exprinterpreter, consdata->sepaconvexconcave.vred_neg_swapped) );
      }
      else
      {
         consdata->sepaconvexconcave.f_neg_swapped = NULL;
         consdata->sepaconvexconcave.f_neg_swapped_yfixed = NULL;
         consdata->sepaconvexconcave.vred_neg_swapped = NULL;
      }

      break;
   }

   default: ;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** frees separation data */
static
SCIP_RETCODE freeSepaData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->f != NULL);

   switch( consdata->convextype )
   {
   case SCIP_BIVAR_CONVEX_CONCAVE:
   {
      if( consdata->sepaconvexconcave.f_yfixed != NULL )
      {
         SCIP_CALL( SCIPexprtreeFree(&consdata->sepaconvexconcave.f_yfixed) );
      }
      if( consdata->sepaconvexconcave.f_neg_swapped != NULL )
      {
         SCIP_CALL( SCIPexprtreeFree(&consdata->sepaconvexconcave.f_neg_swapped) );
      }
      if( consdata->sepaconvexconcave.f_neg_swapped_yfixed != NULL )
      {
         SCIP_CALL( SCIPexprtreeFree(&consdata->sepaconvexconcave.f_neg_swapped_yfixed) );
      }
      if( consdata->sepaconvexconcave.vred != NULL )
      {
         SCIP_CALL( SCIPexprtreeFree(&consdata->sepaconvexconcave.vred) );
      }
      if( consdata->sepaconvexconcave.vred_neg_swapped != NULL )
      {
         SCIP_CALL( SCIPexprtreeFree(&consdata->sepaconvexconcave.vred_neg_swapped) );
      }
      break;
   }

   default: ;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** perturbs a value w.r.t. bounds */
static
void perturb(
   SCIP_Real*            val,                /**< value to perturb on input; perturbed value on output */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub,                 /**< upper bound */
   SCIP_Real             amount              /**< relative amount of perturbation */
   )
{
   SCIP_Real range;
   SCIP_Real mid;

   assert(val != NULL);

   range = ub - lb;
   mid = 0.5 * (lb + ub);

   if( *val < mid )
      *val += MIN(1.0, amount * range);
   else
      *val -= MIN(1.0, amount * range);
}

/** solves an equation f'(s) = constant for a univariate convex or concave function f with respect to bounds on s
 * if there is no s between the bounds such that f'(s) = constant, then it returns the closest bound (and still claims success)
 */
static
SCIP_RETCODE solveDerivativeEquation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< expression tree for f(s) */
   SCIP_Real             targetvalue,        /**< target value for derivative */
   SCIP_Real             lb,                 /**< lower bound on variable */
   SCIP_Real             ub,                 /**< upper bound on variable */
   SCIP_Real*            val,                /**< buffer to store solution value */
   SCIP_Bool*            success             /**< buffer to indicate whether a solution has been found */
   )
{
   SCIP_Real fval;
   SCIP_Real grad;
   SCIP_Real hess;
   SCIP_Real s;
   SCIP_Real nexts;
   SCIP_Real step;
   int iter;

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f != NULL);
   assert(SCIPexprtreeGetInterpreterData(f) != NULL);
   assert(SCIPexprtreeGetNVars(f) == 1);
   assert(val != NULL);
   assert(success != NULL);

   if( SCIPisEQ(scip, lb, ub) )
   {
      *val = lb;
      *success = TRUE;
      return SCIP_OKAY;
   }

   *success = FALSE;

   iter = 0;

   /* start at 0.0, projected onto interior of interval
    * we don't want to start at a bound, because we would not recognize if hessian is 0.0 then
    */
   s = MIN(MAX(0.0, lb), ub);
   perturb(&s, lb, ub, 0.1);

   while( ++iter < NEWTONMAXITER )
   {
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, &s, TRUE, &fval, &grad) );

      /* SCIPdebugMsg(scip, "s = %.20g [%g,%g] f(s) = %g grad = %g\n", s, lb, ub, fval, grad); */

      if( !SCIPisFinite(grad) )
      {
         /* if f cannot be differentiated at s, perturb s to some other point close by
          * for that, we perturb by 0.1 * 2^{-iter}, if iter <= 65, otherwise by 1e-20
          * if that amount is too small to get a change in s, we increase by a factor of 2
          */
         SCIP_Real amount;
         SCIP_Real sold;

         sold = s;
         amount = iter <= 65 ? 0.1 / (1u<<iter) : 1e-20; /*lint !e790*/
         do
         {
            perturb(&s, lb, ub, amount);
            amount *= 2.0;
         } while( s == sold ); /*lint !e777*/

         SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, &s, TRUE, &fval, &grad) );

         /* SCIPdebugMsg(scip, "s = %.20g [%g,%g] f(s) = %g grad = %g (perturbed by %g)\n", s, lb, ub, fval, grad, iter <= 65 ? 0.1 / (1<<iter) : 1e-20); */

         assert(SCIPisFinite(grad));
      }

      if( SCIPisRelEQ(scip, grad, targetvalue) )
      {
         /* if grad is targetvalue (w.r.t. epsilon), then we are done */
         *val = s;
         *success = TRUE;
         break;
      }

      SCIP_CALL( SCIPexprintHessianDense(exprinterpreter, f, &s, FALSE, &fval, &hess) ); /* coverity ignore ARRAY_VS_SINGLETON warning */

      /* SCIPdebugMsg(scip, "s = %.20g [%g,%g] f(s) = %g hess = %g\n", s, lb, ub, fval, hess); */

      if( !SCIPisFinite(hess) )
      {
          SCIP_Real smod;
          SCIP_Real smodval;

          /* if f cannot be two times differentiated at s, take the Hessian from another point close by */
          smod = s;
          perturb(&smod, lb, ub, 0.01);
          SCIP_CALL( SCIPexprintHessianDense(exprinterpreter, f, &smod, TRUE, &smodval, &hess) );

          assert(SCIPisFinite(hess));
      }

      /* next iterate would be s - (grad - targetvalue) / hess */

      if( SCIPisEQ(scip, s, lb) && (grad - targetvalue) * hess >= 0 )
      {
         /* if we are on the left boundary and would go left (or stay), then stop
          * (multiply instead of divide by hess for the case that hess is zero and since only the sign matters
          */
         *val = lb;
         *success = TRUE;
         break;
      }

      if( SCIPisEQ(scip, s, ub) && (grad - targetvalue) * hess <= 0 )
      {
         /* similar, if we are on the right boundary and would go right (or stay), then stop */
         *val = ub;
         *success = TRUE;
         break;
      }

      if( SCIPisZero(scip, hess) )
      {
         /* hmm, stationary point, don't know how to continue; thus, give up */
         break;
      }

      if( SCIPisZero(scip, (grad - targetvalue) / hess) && SCIPisFeasEQ(scip, grad, targetvalue) )
      {
         /* if grad is targetvalue (w.r.t. feastol) and step length would be almost 0, then we are also done */
         *val = s;
         *success = TRUE;
         break;
      }

      /* @todo we could also implement a damped Newton method if the step is too large */
      step = (grad - targetvalue) / hess;
      assert(step != 0.0);

      nexts = s - step;
      while( s == nexts ) /*lint !e777*/
      {
         /* if steplength is so tiny that there is no change in s, go by 1e-9 into given direction */
         step *= 2.0;
         nexts = s - step;
      }
      assert(nexts != s); /*lint !e777*/
      s = nexts;

      if( s < lb )
         s = lb;
      else if( s > ub )
         s = ub;
   }

   return SCIP_OKAY;
}

/** generates a cut for f(x,y) + c*z <= rhs with f(x,y) being convex or 1-convex with x or y fixed or convex-concave with y fixed
 * f(x0, y0) + <grad, (x,y)-(x0,y0)> + c*z <= rhs, where grad is gradient of f in (x0, y0)
 */
static
SCIP_RETCODE generateLinearizationCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            x0y0,               /**< value of x and y variables where to generate cut */
   SCIP_Bool             newxy,              /**< whether the last evaluation of f(x,y) with the expression interpreter was at (x0, y0) */
   SCIP_ROW**            row                 /**< storage for cut */
   )
{
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_CONSDATA* consdata;
   char           rowname[SCIP_MAXSTRLEN];
   SCIP_Real      fval;
   SCIP_Real      fgrad[2];
   SCIP_Real      rhs;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisInfinity(scip, consdata->rhs));
   assert(newxy || SCIPexprtreeGetInterpreterData(consdata->f) != NULL);

   /* compile expression if evaluated the first time; can only happen if newxy is FALSE */
   if( newxy && SCIPexprtreeGetInterpreterData(consdata->f) == NULL )
   {
      SCIP_CALL( SCIPexprintCompile(exprinterpreter, consdata->f) );
   }

   x = SCIPexprtreeGetVars(consdata->f)[0];
   y = SCIPexprtreeGetVars(consdata->f)[1];

   assert(consdata->convextype == SCIP_BIVAR_ALLCONVEX ||
      (consdata->convextype == SCIP_BIVAR_1CONVEX_INDEFINITE && (SCIPisEQ(scip, SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x)) || SCIPisEQ(scip, SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y)))) ||
      (consdata->convextype == SCIP_BIVAR_CONVEX_CONCAVE && SCIPisEQ(scip, SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y))) );

   /* compute f(x,y) and gradient of f in (x, y) */
   SCIP_CALL( SCIPexprintGrad(exprinterpreter, consdata->f, x0y0, newxy, &fval, fgrad) );

   if( !SCIPisFinite(fval) || !SCIPisFinite(fgrad[0]) || !SCIPisFinite(fgrad[1]) )
   {
      perturb(&x0y0[0], SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), 0.001);
      perturb(&x0y0[1], SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y), 0.001);

      SCIP_CALL( SCIPexprintGrad(exprinterpreter, consdata->f, x0y0, TRUE, &fval, fgrad) );

      if( !SCIPisFinite(fval) || !SCIPisFinite(fgrad[0]) || !SCIPisFinite(fgrad[1]) )
      {
         SCIPdebugMsg(scip, "could not evaluate f at given reference point and perturbed one");
         *row = NULL;
         return SCIP_OKAY;
      }
   }

   rhs = consdata->rhs - fval + fgrad[0] * x0y0[0] + fgrad[1] * x0y0[1];

   /* setup SCIP row */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_linearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, row, SCIPconsGetHdlr(cons), rowname, -SCIPinfinity(scip), rhs, FALSE, FALSE /* modifiable */, TRUE /* removable */) );

   SCIP_CALL( SCIPaddVarsToRow(scip, *row, 2, SCIPexprtreeGetVars(consdata->f), fgrad) );

   if( consdata->z != NULL )
      SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->z, consdata->zcoef) );

   return SCIP_OKAY;
}

/** given a convex (concave, resp.) bivariate function, computes an over- (under-, resp.) estimating hyperplane
 *  does not succeed if some variable is unbounded or both variables are fixed
 */
static
SCIP_RETCODE generateEstimatingHyperplane(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expression interpreter */
   SCIP_EXPRTREE*        f,                  /**< bivariate function to compute under or overestimator for */
   SCIP_Bool             doover,             /**< whether to compute an overestimator (TRUE) or an underestimator (FALSE) */
   SCIP_Real*            x0y0,               /**< reference values for nonlinear variables */
   SCIP_Real*            coefx,              /**< coefficient of x in estimator */
   SCIP_Real*            coefy,              /**< coefficient of y in estimator */
   SCIP_Real*            constant,           /**< constant part of estimator */
   SCIP_Bool*            success             /**< pointer to indicate whether coefficients where successfully computed */
   )
{
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      ylb;
   SCIP_Real      yub;

   SCIP_Real      p1[2];
   SCIP_Real      p2[2];
   SCIP_Real      p3[2];
   SCIP_Real      p4[2];
   SCIP_Real      p1val;
   SCIP_Real      p2val;
   SCIP_Real      p3val;
   SCIP_Real      p4val;

   SCIP_Real      alpha;
   SCIP_Real      beta;
   SCIP_Real      gamma_;
   SCIP_Real      delta;

   SCIP_Bool      tryother;

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f    != NULL);
   assert(x0y0 != NULL);
   assert(coefx != NULL);
   assert(coefy != NULL);
   assert(constant != NULL);
   assert(success != NULL);

   *success = FALSE;

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);
   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   /* reference point should not be outside of bounds */
   assert(SCIPisLE(scip, xlb, x0y0[0]));
   assert(SCIPisGE(scip, xub, x0y0[0]));
   assert(SCIPisLE(scip, ylb, x0y0[1]));
   assert(SCIPisGE(scip, yub, x0y0[1]));

   if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) || SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub) )
   {
      SCIPdebugMsg(scip, "skip estimating hyperplane since <%s> or <%s> is unbounded\n", SCIPvarGetName(x), SCIPvarGetName(y));
      return SCIP_OKAY;
   }

   if( SCIPisEQ(scip, xlb, xub) && SCIPisEQ(scip, ylb, yub) )
   {
      SCIPdebugMsg(scip, "skip estimating hyperplane since both <%s> and <%s> are fixed\n", SCIPvarGetName(x), SCIPvarGetName(y));
      return SCIP_OKAY;
   }

   /* unten links */
   p1[0] = xlb;
   p1[1] = ylb;

   /* unten rechts */
   p2[0] = xub;
   p2[1] = ylb;

   /* oben rechts */
   p3[0] = xub;
   p3[1] = yub;

   /* oben links */
   p4[0] = xlb;
   p4[1] = yub;

   if( SCIPisEQ(scip, xlb, xub) )
   {
      /* secant between p1 and p4: p1val + [(p4val - p1val) / (yub - ylb)] * (y - ylb) */
      assert(!SCIPisEQ(scip, ylb, yub));

      SCIP_CALL( SCIPexprintEval(exprinterpreter, f, p1, &p1val) );
      SCIP_CALL( SCIPexprintEval(exprinterpreter, f, p4, &p4val) );

      if( !SCIPisFinite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !SCIPisFinite(p4val) || SCIPisInfinity(scip, REALABS(p4val)) )
      {
         SCIPdebugMsg(scip, "skip hyperplane since function cannot be evaluated\n");
         return SCIP_OKAY;
      }

      *coefx = 0.0;
      *coefy = (p4val - p1val) / (yub - ylb);
      *constant = p1val - *coefy * ylb;

      *success = TRUE;

      return SCIP_OKAY;
   }

   if( SCIPisEQ(scip, ylb, yub) )
   {
      /* secant between p1 and p2: p1val + [(p2val - p1val) / (xub - xlb)] * (x - xlb) */
      assert(!SCIPisEQ(scip, xlb, xub));

      SCIP_CALL( SCIPexprintEval(exprinterpreter, f, p1, &p1val) );
      SCIP_CALL( SCIPexprintEval(exprinterpreter, f, p2, &p2val) );

      if( !SCIPisFinite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !SCIPisFinite(p2val) || SCIPisInfinity(scip, REALABS(p2val)) )
      {
         SCIPdebugMsg(scip, "skip hyperplane since function cannot be evaluated\n");
         return SCIP_OKAY;
      }

      *coefx = (p2val - p1val) / (xub - xlb);
      *coefy = 0.0;
      *constant = p1val - *coefx * xlb;

      *success = TRUE;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPexprintEval(exprinterpreter, f, p1, &p1val) );
   SCIP_CALL( SCIPexprintEval(exprinterpreter, f, p2, &p2val) );
   SCIP_CALL( SCIPexprintEval(exprinterpreter, f, p3, &p3val) );
   SCIP_CALL( SCIPexprintEval(exprinterpreter, f, p4, &p4val) );

   /* if we want an underestimator, flip f(x,y), i.e., do as if we compute an overestimator for -f(x,y) */
   if( !doover )
   {
      p1val = -p1val;
      p2val = -p2val;
      p3val = -p3val;
      p4val = -p4val;
   }

   SCIPdebugMsg(scip, "p1 = (%g, %g), f(p1) = %g\n", p1[0], p1[1], p1val);
   SCIPdebugMsg(scip, "p2 = (%g, %g), f(p2) = %g\n", p2[0], p2[1], p2val);
   SCIPdebugMsg(scip, "p3 = (%g, %g), f(p3) = %g\n", p3[0], p3[1], p3val);
   SCIPdebugMsg(scip, "p4 = (%g, %g), f(p4) = %g\n", p4[0], p4[1], p4val);

   if( !SCIPisFinite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !SCIPisFinite(p2val) || SCIPisInfinity(scip, REALABS(p2val)) ||
      ! SCIPisFinite(p3val) || SCIPisInfinity(scip, REALABS(p3val)) || !SCIPisFinite(p4val) || SCIPisInfinity(scip, REALABS(p4val)) )
   {
      SCIPdebugMsg(scip, "skip hyperplane since function cannot be evaluated\n");
      return SCIP_OKAY;
   }

   /* compute coefficients alpha, beta, gamma (>0), delta such that
    *   alpha*x + beta*y + gamma*z = delta
    * is satisfied by at least three of the corner points (p1,f(p1)), ..., (p4,f(p4)) and
    * the fourth corner point lies below this hyperplane.
    * Since we assume that f is convex, we then know that all points (x,y,f(x,y)) are below this hyperplane, i.e.,
    *    alpha*x + beta*y - delta <= -gamma * f(x,y),
    * or, equivalently,
    *   -alpha/gamma*x - beta/gamma*y + delta/gamma >= f(x,y).
    */

   tryother = FALSE;
   if( x0y0[1] <= ylb + (yub - ylb)/(xub - xlb) * (x0y0[0] - xlb) )
   {
      SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1], p3val, &alpha,
            &beta, &gamma_, &delta) );

      assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
      assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
      assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));

      /* if hyperplane through p1,p2,p3 does not overestimate f(p4), then it must be the other variant */
      if( SCIPisInfinity(scip, delta) || alpha * p4[0] + beta * p4[1] + gamma_ * p4val > delta )
         tryother = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1], p4val, &alpha,
            &beta, &gamma_, &delta) );

      assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
      assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
      assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

      /* if hyperplane through p1,p3,p4 does not overestimate f(p2), then it must be the other variant */
      if( SCIPisInfinity(scip, delta) || alpha * p2[0] + beta * p2[1] + gamma_ * p2val > delta )
         tryother = TRUE;
   }

   if( tryother )
   {
      if( x0y0[1] <= yub + (ylb - yub)/(xub - xlb) * (x0y0[0] - xlb) )
      {
         SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );

         /* hyperplane should be above (p3,f(p3)) and other points should lie on hyperplane */
         assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
         assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
         assert(SCIPisInfinity(scip, delta) || SCIPisFeasLE(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
         assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));
      }
      else
      {
         SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );

         /* hyperplane should be above (p1,f(p1)) and other points should lie on hyperplane */
         assert(SCIPisInfinity(scip, delta) || SCIPisFeasLE(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
         assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
         assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
         assert(SCIPisInfinity(scip, delta) || SCIPisFeasEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));
      }
   }

   SCIPdebugMsg(scip, "alpha = %g, beta = %g, gamma = %g, delta = %g\n", alpha, beta, gamma_, delta);

   /* check if bad luck: should not happen if xlb != xub and ylb != yub and numerics are fine */
   if( SCIPisInfinity(scip, delta) || SCIPisZero(scip, gamma_) )
      return SCIP_OKAY;
   assert(!SCIPisNegative(scip, gamma_));

   /* flip hyperplane */
   if( !doover )
      gamma_ = -gamma_;

   *coefx    = -alpha / gamma_;
   *coefy    = -beta  / gamma_;
   *constant =  delta / gamma_;

   *success = TRUE;

   return SCIP_OKAY;
}

/** generates a cut for lhs <= f(x,y) + c*z with f(x,y) being convex */
static
SCIP_RETCODE generateOverestimatingHyperplaneCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            x0y0,               /**< reference values for nonlinear variables */
   SCIP_ROW**            row                 /**< storage for cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real coefs[2];
   SCIP_Real constant = SCIP_INVALID;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row != NULL);

   *row = NULL;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( generateEstimatingHyperplane(scip, exprinterpreter, consdata->f, TRUE, x0y0, &coefs[0], &coefs[1], &constant, &success) );

   if( success )
   {
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      assert(SCIPisFinite(coefs[0]));
      assert(SCIPisFinite(coefs[1]));
      assert(SCIPisFinite(constant));

      SCIP_CALL( SCIPcreateRowCons(scip, row, SCIPconsGetHdlr(cons), "bivaroveresthyperplanecut", 0, NULL, NULL, consdata->lhs - constant, SCIPinfinity(scip), TRUE, FALSE, TRUE) );

      SCIP_CALL( SCIPaddVarsToRow(scip, *row, 2, SCIPexprtreeGetVars(consdata->f), coefs) );
      if( consdata->z != NULL )
         SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->z, consdata->zcoef) );
   }
   else
   {
      SCIPdebugMsg(scip, "failed to compute overestimator for all-convex constraint <%s>\n", SCIPconsGetName(cons));
   }

   return SCIP_OKAY;
}

/** generates a linear underestimator for f(x,y)
 * when the generators of the underestimating segment
 * are contained in y=ylb and y=yub.
 * Generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that
 * alpha * x + beta * y - delta <= gamma * f(x,y)
 */
static
SCIP_RETCODE generateUnderestimatorParallelYFacets(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< function f(x,y) */
   SCIP_Real*            xyref,              /**< reference values for x and y */
   SCIP_Real             cutcoeff[4],        /**< cut coefficients alpha, beta, gamma, delta */
   SCIP_Real*            convenvvalue,       /**< function value of the convex envelope */
   SCIP_Bool*            success             /**< buffer to store whether coefficients were successfully computed */
   )
{
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xval;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      yval;
   SCIP_Real      ylb;
   SCIP_Real      yub;

   SCIP_Real      t;
   SCIP_EXPR*     vred;
   SCIP_EXPRTREE* vredtree;
   SCIP_EXPR*     e1;
   SCIP_EXPR*     e2;
   SCIP_EXPR*     tmp;
   SCIP_EXPR*     tmp2;
   SCIP_EXPR*     subst[2];

   SCIP_Real      sval;
   SCIP_Real      slb;
   SCIP_Real      sub;
   SCIP_Real      rval;

   SCIP_Real      frval;
   SCIP_Real      fsval;
   SCIP_Real      x0y0[2];
   SCIP_Real      grad[2];

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f != NULL);
   assert(xyref != NULL);
   assert(success != NULL);

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   xval = xyref[0];
   yval = xyref[1];

   *success = FALSE;

   /* check that variables are not unbounded or fixed and reference point is in interior */
   assert(!SCIPisInfinity(scip, -xlb));
   assert(!SCIPisInfinity(scip,  xub));
   assert(!SCIPisInfinity(scip, -ylb));
   assert(!SCIPisInfinity(scip,  yub));
   assert(!SCIPisEQ(scip,xlb,xub));
   assert(!SCIPisEQ(scip,ylb,yub));
   assert(!SCIPisEQ(scip,xlb,xval));
   assert(!SCIPisEQ(scip,xub,xval));
   assert(!SCIPisEQ(scip,ylb,yval));
   assert(!SCIPisEQ(scip,yub,yval));

   SCIPdebugMsg(scip, "f(%s, %s) = ", SCIPvarGetName(x), SCIPvarGetName(y));
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(f, SCIPgetMessagehdlr(scip), NULL) ) );
   SCIPdebugMsgPrint(scip, "\n");

   t = (yub - yval) / (yub - ylb);

   /* construct v_red(s) := t f(1/t xval + (1-1/t) s, ylb) + (1-t) f(s, yub) */

   /* construct e1 := f(1/t xval + (1-1/t) s, ylb) */
   SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e1, SCIPexprtreeGetRoot(f)) );          /* e1 = f(x,y) */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_VARIDX, 0) );             /* tmp  = s */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp2, SCIP_EXPR_CONST, 1.0 - 1.0 / t) );  /* tmp2 = 1-1/t */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_MUL, tmp, tmp2) );        /* tmp  = (1-1/t)*s */
   if( xval != 0.0 )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp2, SCIP_EXPR_CONST, 1/t*xval) );    /* tmp2 = 1/t*xval */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_PLUS, tmp, tmp2) );    /* tmp = 1/t*xval + (1-1/t)*s */
   }
   subst[0] = tmp;

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_CONST, ylb) );        /* tmp = ylb */

   assert(SCIPexprGetOperator(e1) != SCIP_EXPR_VARIDX);  /* substitute cannot substitute the root node, but f should not be a single variable anyway */
   SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e1, subst) );                      /* e1 = f(1/t*xval + (1-1/t)*s, ylb) */

   SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
   SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);

   /* construct e2 := f(s, yub) */
   SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e2, SCIPexprtreeGetRoot(f)) );          /* e2 = f(x,y) */

   subst[0] = NULL;

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_CONST, yub) );

   assert(SCIPexprGetOperator(e2) != SCIP_EXPR_VARIDX);  /* substitute cannot substitute the root node, but f should not be a single variable anyway */
   SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e2, subst) );                      /* e2 = f(s,yub) */

   SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);

   /* construct vred := t * e1 + (1-t) * e2 */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, t) );               /* tmp = t */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e1, SCIP_EXPR_MUL, e1, tmp) );            /* e1 = t * f(1/t*xval+(1-1/t)*s,ylb) */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, 1.0 - t) );         /* tmp = 1 - t */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &e2, SCIP_EXPR_MUL, e2, tmp) );            /* e2 = (1-t) * f(s, yub) */

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &vred, SCIP_EXPR_PLUS, e1, e2) );
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &vredtree, vred, 1, 0, NULL) );

   SCIP_CALL( SCIPexprintCompile(exprinterpreter, vredtree) );

   /* compute bounds on s */
   slb = (yval - yub) / (ylb - yval) * (xval / t - xub);
   sub = (yval - yub) / (ylb - yval) * (xval / t - xlb);
   if( slb < xlb )
      slb = xlb;
   if( sub > xub )
      sub = xub;

   /* find s in [slb, sub] such that vred'(s) = 0 */
   SCIP_CALL( solveDerivativeEquation(scip, exprinterpreter, vredtree, 0.0, slb, sub, &sval, success) );

   SCIP_CALL( SCIPexprtreeFree(&vredtree) );

   if( *success == FALSE )
   {
      /* something went wrong when computing s */
      return SCIP_OKAY;
   }

   /* compute r from s */
   rval = 1.0 / t * xval + (1.0 - 1.0 / t) * sval;
   assert(SCIPisFeasGE(scip, rval, xlb));
   assert(SCIPisFeasLE(scip, rval, xub));
   rval = MAX(xlb, MIN(rval, xub));

   /* compute f(sval, yub) */
   x0y0[0] = sval;
   x0y0[1] = yub;
   SCIP_CALL( SCIPexprtreeEval(f, x0y0, &fsval) );

   /* compute f(rval, ylb) */
   x0y0[0] = rval;
   x0y0[1] = ylb;
   SCIP_CALL( SCIPexprtreeEval(f, x0y0, &frval) );

   if( !SCIPisEQ(scip, sval, xlb) && !SCIPisEQ(scip, sval, xub) )
   {
      x0y0[0] = sval;
      x0y0[1] = yub;

      /* compute f'(xbar, ybar) */
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fsval, grad) );
   }
   else if( !SCIPisEQ(scip, rval, xlb) && !SCIPisEQ(scip, rval, xub) )
   {
      x0y0[0] = rval;
      x0y0[1] = ylb;

      /* compute f'(xbar, ybar) */
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &frval, grad) );
   }
   else
   {
      /* rare case
       * both points (sval, yub) and (rval, ylb) should yield valid inequality
       * for now, just take the first one, if differentiable, otherwise second one */
      x0y0[0] = sval;
      x0y0[1] = yub;

      /* compute f'(xbar, ybar) */
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fsval, grad) );

      if( !SCIPisFinite(grad[0]) )
      {
         x0y0[0] = rval;
         x0y0[1] = ylb;

         /* compute f'(xbar, ybar) */
         SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &frval, grad) );
      }
   }

   /* compute vred(s) = t * f(rval, ylb) + (1-t) * f(s, yub) */
   /* SCIP_CALL( SCIPexprtreeEval(vredtree, &sval, &vredval) ); */
   *convenvvalue = t * frval + (1.0 - t) * fsval;

   SCIPdebugMsg(scip, "Parallel: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
   SCIPdebugMsg(scip, "Parallel: r=%g in [%g,%g], s=%g in [%g,%g], f(r,ylb)=%g, f(xlb,s)=%g\n",rval,xlb,xub,sval,ylb,yub,frval,fsval);
   SCIPdebugMsg(scip, "(r,ylb)=(%g,%g), (s,yub)=(%g,%g), vredval=%g\n",rval,ylb,sval,yub,*convenvvalue);

   if( !SCIPisFinite(grad[0]) || SCIPisInfinity(scip, REALABS(grad[0])) )
   {
      SCIPdebugMsg(scip, "f not differentiable in (x0,y0) w.r.t. x\n");
      return SCIP_OKAY;
   }

   /* compute cut coefficients */
   cutcoeff[0]   = (yub - ylb) * grad[0];
   cutcoeff[1]   = fsval - frval - (sval - rval) * grad[0];
   cutcoeff[2]   = yub - ylb;
   cutcoeff[3]   = cutcoeff[0] * xval + cutcoeff[1] * yval - cutcoeff[2] * *convenvvalue;

   SCIPdebugMsg(scip, "Parallel: cutcoeff[0]=%g, cutcoeff[1]=%g, cutcoeff[2]=1.0, cutcoeff[3]=%g\n",cutcoeff[0]/cutcoeff[2],cutcoeff[1]/cutcoeff[2],cutcoeff[3]/cutcoeff[2]);

   *success = TRUE;

   return SCIP_OKAY;
}


/** generates a linear underestimator for f(x,y)
 * with f(x,y) being convex in x and convex in y.
 * The segmenent connects orthogonal facets: Either (x=l_x,y=l_y)
 * or (x=u_x,y=u_y).
 * generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that
 * alpha * x + beta * y - delta <= gamma * f(x,y)
 */
static
SCIP_RETCODE generateOrthogonal_lx_ly_Underestimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< function f(x,y) */
   SCIP_Real*            xyref,              /**< reference values for x and y */
   SCIP_Real             cutcoeff[4],        /**< cut coefficients alpha, beta, gamma, delta */
   SCIP_Real*            convenvvalue,       /**< function value of the convex envelope */
   SCIP_Bool*            success             /**< buffer to store whether coefficients were successfully computed */
   )
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_Real xval;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real yval;
   SCIP_Real ylb;
   SCIP_Real yub;

   SCIP_Real x0y0[2];

   SCIP_EXPR* vred;
   SCIP_EXPRTREE* vredtree;
   SCIP_EXPR* e1;
   SCIP_EXPR* e2;
   SCIP_EXPR* tmp;
   SCIP_EXPR* expr;
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;
   SCIP_EXPR* subst[2];

   SCIP_Real tval, tlb, tub;
   SCIP_Real sval;
   SCIP_Real rval;

   SCIP_Real frval,fsval;
   SCIP_Real grad_rval[2];
   SCIP_Real grad_sval[2];

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f != NULL);
   assert(convenvvalue != NULL);
   assert(success != NULL);

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   xval = xyref[0];
   yval = xyref[1];

   /* check that variables are not unbounded or fixed and reference point is in interior */
   assert(!SCIPisInfinity(scip, -xlb));
   assert(!SCIPisInfinity(scip,  xub));
   assert(!SCIPisInfinity(scip, -ylb));
   assert(!SCIPisInfinity(scip,  yub));
   assert(!SCIPisEQ(scip,xlb,xub));
   assert(!SCIPisEQ(scip,ylb,yub));
   assert(!SCIPisEQ(scip,xlb,xval));
   assert(!SCIPisEQ(scip,xub,xval));
   assert(!SCIPisEQ(scip,ylb,yval));
   assert(!SCIPisEQ(scip,yub,yval));

   *success = FALSE;

   SCIPdebugMsg(scip, "f(%s, %s) = ", SCIPvarGetName(x), SCIPvarGetName(y));
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(f, SCIPgetMessagehdlr(scip), NULL) ) );
   SCIPdebugMsgPrint(scip, "\n");
   SCIPdebugMsg(scip, "%s[%g,%g] = %g  %s[%g,%g] = %g\n", SCIPvarGetName(x), xlb, xub, xval, SCIPvarGetName(y), ylb, yub, yval);

   /* check in which triangle the point (xval,yval) lies */
   if( yval <= (ylb-yub) / (xub-xlb) * (xval-xlb) + yub )
   {
      /* (xval,yval) lies in lower left triangle, i.e. region A_1 */
      /* construct v_red(t) := t f( xlb, (yval-(1-t)ylb)/t ) + (1-t)*f( (xval-xlb*t)/(1-t), ylb ) */

      /* construct e1 := f(xlb, ylb + (yval-ylb)/t) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );           /* expr = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, yval-ylb) );     /* tmp  = yval-ylb */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, tmp, expr) );      /* expr = (yval-ylb) / t */
      if( ylb != 0.0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, ylb) );       /* tmp = ylb */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PLUS, expr, tmp) );  /* expr = ylb + (yval-ylb) / t */
      }
      subst[1] = expr;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[0], SCIP_EXPR_CONST, xlb) );      /* subst[0] = xlb */

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e1, SCIPexprtreeGetRoot(f)) );        /* e1 = f(x,y) */
      assert(SCIPexprGetOperator(e1) != SCIP_EXPR_VARIDX);  /* expr substitute vars cannot substitute the root node, but f should not be a single variable anyway */
      SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e1, subst) );                    /* e1 = f(xlb, ylb + (yval-ylb)/t) */

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);


      /* construct e2 := f((xval-xlb*t)/(1-t), ylb) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_VARIDX, 0) );          /* expr1 = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, 1.0) );         /* tmp   = 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_MINUS, tmp, expr1) );  /* expr1 = 1-t */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_VARIDX, 0) );          /* expr2 = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, xlb) );         /* tmp   = xlb */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MUL, expr2, tmp) );    /* expr2 = xlb * t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, xval) );        /* tmp   = xval */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MINUS, tmp, expr2) );  /* expr2 = xval - xlb * t */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_DIV, expr2, expr1) );  /* expr =  (xval-t*xlb)/(1-t) */
      subst[0] = expr;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_CONST, ylb) );      /* subst[0] = ylb */

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e2, SCIPexprtreeGetRoot(f)) );        /* e2 = f(x,y) */
      assert(SCIPexprGetOperator(e2) != SCIP_EXPR_VARIDX);  /* expr substitute vars cannot substitute the root node, but f should not be a single variable anyway */
      SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e2, subst) );                    /* e2 = f((xval-xlb*t)/(1-t), ylb) */

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);


      /* construct vred := t * e1 + (1-t) * e2 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_VARIDX, 0) );          /* expr  = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_MUL, expr, e1) );      /* expr1 = t * e1*/

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );           /* expr  = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, 1.0) );          /* tmp   = 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_MINUS, tmp, expr) );    /* expr  = 1 - t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MUL, expr, e2) );      /* expr2 = (1-t) * e2 */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &vred, SCIP_EXPR_PLUS, expr1, expr2) );
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &vredtree, vred, 1, 0, NULL) );
      SCIP_CALL( SCIPexprintCompile(exprinterpreter, vredtree) );

      /* compute bounds on t */
      tlb = (yval-ylb)/(yub-ylb);
      tub = (xub-xval)/(xub-xlb);

      /* find t in [lambalb, tub] such that vred'(t) = 0 */
      SCIP_CALL( solveDerivativeEquation(scip, exprinterpreter, vredtree, 0.0, tlb, tub, &tval, success) );

      /* computing the cut coefficients */
      if( *success == FALSE )
      {
         /* something went wrong when computing s */
         SCIP_CALL( SCIPexprtreeFree(&vredtree) );
         return SCIP_OKAY;
      }

      /* compute r and s from tval */
      rval = (yval-(1-tval)*ylb)/tval;
      rval = MAX(ylb, MIN(yub, rval));
      sval = (xval-xlb*tval)/(1-tval);
      sval = MAX(xlb, MIN(xub, sval));

      SCIPdebugMsg(scip, "LowerLeft: t[%g,%g] = %g -> r = %g, s = %g\n",tlb,tub,tval,rval,sval);

      /* compute vred(tval) */
      SCIP_CALL( SCIPexprtreeEval(vredtree, &tval, convenvvalue) );

      SCIP_CALL( SCIPexprtreeFree(&vredtree) );

      /* compute f(s, ylb) and f'(s, ylb) */
      x0y0[0] = sval;
      x0y0[1] = ylb;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fsval, grad_sval) );

      /* compute f(xlb, r) and f'(xlb,r) */
      x0y0[0] = xlb;
      x0y0[1] = rval;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &frval, grad_rval) );

      /* generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that
       * alpha * x + beta * y - delta <= gamma * f(x,y)
       * cf. Section 2.5.2 Aux.prob. 2 case (ii)
       */
      if( !SCIPisEQ(scip, sval, xub) )
      {
         /* use the x-axis to determine the second direction */
         if( !SCIPisFinite(grad_sval[0]) || SCIPisInfinity(scip, REALABS(grad_sval[0])) )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (rval-ylb) * grad_sval[0];
         cutcoeff[1] = (sval-xlb) * grad_sval[0] + frval - fsval;
         cutcoeff[2] = rval-ylb;
         cutcoeff[3] = cutcoeff[0]*xlb+cutcoeff[1]*rval-cutcoeff[2]*frval;
      }
      else if( !SCIPisEQ(scip,rval,yub) )
      {
         /* use the y-axis to determine the second direction */
         if( !SCIPisFinite(grad_rval[1]) || SCIPisInfinity(scip, REALABS(grad_rval[1])) )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (rval-ylb)*grad_rval[1]+fsval-frval;
         cutcoeff[1] = (sval-xlb)*grad_rval[1];
         cutcoeff[2] = sval-xlb;
         cutcoeff[3] = cutcoeff[0]*xlb+cutcoeff[1]*rval-cutcoeff[2]*frval;
      }
      else
      {
         /* the point lies on the segment between (xlb,yub) and (xub,ylb) */
         if( !SCIPisFinite(grad_sval[0]) || !SCIPisFinite(grad_rval[0]) || SCIPisInfinity(scip, REALABS(MIN(grad_sval[0],grad_rval[0]))) )
         {
            /* FIXME maybe it is sufficient to have one of them finite, using that one for the MIN below? */
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (rval-ylb)* MIN(grad_sval[0],grad_rval[0]);
         cutcoeff[1] = (sval-xlb)* MIN(grad_sval[0],grad_rval[0])+frval-fsval;
         cutcoeff[2] = (rval-ylb);
         cutcoeff[3] = cutcoeff[0]*xlb+cutcoeff[1]*rval-cutcoeff[2]*frval;
      }

      SCIPdebugMsg(scip, "LowerLeft: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
      SCIPdebugMsg(scip, "LowerLeft: r=%g in [%g,%g], s=%g in [%g,%g], f(s,ylb)=%g, f(xlb,r)=%g\n",rval,xlb,xub,sval,ylb,yub,fsval,frval);
      SCIPdebugMsg(scip, "(s,ylb)=(%g,%g) (xlb,r)=(%g,%g) t=%g, vredval=%g\n",sval,ylb,xlb,rval,tval,*convenvvalue);
      SCIPdebugMsg(scip, "LowerLeft: cutcoeff[0]=%g, cutcoeff[1]=%g,cutcoeff[2]=%g,cutcoeff[3]=%g\n",cutcoeff[0],cutcoeff[1],cutcoeff[2],cutcoeff[3]);
   }
   else
   {
      /* (xval,yval) lies in the upper right triangle, i.e region A_2 */
      /* construct v_red(t) := t f( xub, yub + (yval-yub)/t ) + (1-t)*f((xval-xub*t)/(1-t), yub) */

      /* construct e1 := f(xub, yub+(yval-yub)/t) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );           /* expr     =  t*/
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, yval-yub) );     /* tmp  = yval-yub*/
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, tmp, expr) );      /* expr = (yval-yub) / t */
      if( yub != 0.0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, yub) );       /* tmp = yub */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PLUS, expr, tmp) );  /* expr = yub + (yval-yub)/t */
      }
      subst[1] = expr;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[0], SCIP_EXPR_CONST, xub) );      /* tmp = xub */

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e1, SCIPexprtreeGetRoot(f)) );        /* e1 = f(x,y) */
      assert(SCIPexprGetOperator(e1) != SCIP_EXPR_VARIDX); /* cannot substitute root */
      SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e1, subst) );                    /* e1 = f(xub, yub + (yval-yub)/t) */

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);

      /* construct e2 := f((xval-t*xub)/(1-t), yub) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_VARIDX, 0) );          /* expr1 = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, 1.0) );         /* tmp = 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_MINUS, tmp, expr1) );  /* expr1 = 1-t */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_VARIDX, 0) );          /* expr2 = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, xub) );         /* tmp   = xub */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MUL, expr2, tmp) );    /* expr2 = xub * t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, xval) );        /* tmp   = xval */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MINUS, tmp, expr2) );  /* expr2 = xval - xub * t */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, expr2, expr1) );   /* expr =  (xval-t*xub)/(1-t) */
      subst[0] = expr;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_CONST, yub) );      /* tmp = yub */

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e2, SCIPexprtreeGetRoot(f)) );        /* e2 = f(x,y) */
      assert(SCIPexprGetOperator(e2) != SCIP_EXPR_VARIDX); /* cannot substitute root */
      SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e2, subst) );                    /* e2 =  f((xval-t*xub)/(1-t), yub) */

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);

      /* construct vred := t * e1 + (1-t) * e2 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_VARIDX, 0) );          /* expr  = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_MUL, e1, expr) );      /* expr1 = t * e1*/

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_VARIDX, 0) );          /* expr  = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, 1.0) );         /* tmp   = 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MINUS, tmp, expr) );   /* expr  = 1-t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MUL, e2, expr) );      /* expr2 = (1-t) * e2*/


      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &vred, SCIP_EXPR_PLUS, expr1, expr2) );
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &vredtree, vred, 1, 0, NULL) );
      SCIP_CALL( SCIPexprintCompile(exprinterpreter, vredtree) );

      /* compute bounds on t */
      tlb = (yub-yval)/(yub-ylb);
      tub = (xval-xlb)/(xub-xlb);

      /* find t in [tlb, tub] such that vred'(t) = 0 */
      SCIP_CALL( solveDerivativeEquation(scip, exprinterpreter, vredtree, 0.0, tlb, tub, &tval, success) );

      SCIP_CALL( SCIPexprtreeFree(&vredtree) );

      if( *success == FALSE )
      {
         /* something went wrong when computing s */
         return SCIP_OKAY;
      }

      /* computing the cut coefficients */

      /* compute r and s from tval */
      rval = (yval-(1.0-tval)*yub)/tval;
      assert(SCIPisFeasGE(scip, rval, ylb));
      assert(SCIPisFeasLE(scip, rval, yub));
      rval = MAX(ylb, MIN(yub, rval));

      sval = (xval-xub*tval)/(1.0-tval);
      assert(SCIPisFeasGE(scip, sval, xlb));
      assert(SCIPisFeasLE(scip, sval, xub));
      sval = MAX(xlb, MIN(xub, sval));

      /* compute f(xub,r) and f'(xub,r) */
      x0y0[0] = xub;
      x0y0[1] = rval;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &frval, grad_rval) );

      /* compute f(s,yub) and f'(s,yub) */
      x0y0[0] = sval;
      x0y0[1] = yub;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fsval, grad_sval) );

      /* compute vred(tval) */
      *convenvvalue = tval * frval + (1.0-tval) * fsval;

      /* generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that
       * alpha * x + beta * y - delta <= gamma * f(x,y) */

      if( !SCIPisEQ(scip, sval, xlb) )
      {
         /* use the x-axis to determine the second direction */
         if( !SCIPisFinite(grad_sval[0]) || SCIPisInfinity(scip, REALABS(grad_sval[0])) )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }

         cutcoeff[0] = (yub-rval)*grad_sval[0];
         cutcoeff[1] = (xub-sval)*grad_sval[0]+fsval-frval;
         cutcoeff[2] = yub-rval;
         cutcoeff[3] = cutcoeff[0]*sval+cutcoeff[1]*yub-cutcoeff[2]*fsval;
      }
      else if( !SCIPisEQ(scip,rval,ylb) )
      {
         /* use the y-axis to determine the second direction */
         if( !SCIPisFinite(grad_rval[1]) || SCIPisInfinity(scip, REALABS(grad_rval[1])) )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (yub-rval)*grad_rval[1]+frval-fsval;
         cutcoeff[1] = (xub-sval)*grad_rval[1];
         cutcoeff[2] = xub-sval;
         cutcoeff[3] = cutcoeff[0]*sval+cutcoeff[1]*yub-cutcoeff[2]*fsval;
      }
      else
      {
         /* the point lies on the segment between (xlb,yub) and (xub,ylb)
          * due to numerics, we get into this case here instead in the LowerLeft
          */
         assert(SCIPisFeasLE(scip, yval, (ylb-yub) / (xub-xlb) * (xval-xlb) + yub));
         if( !SCIPisFinite(grad_sval[0]) || !SCIPisFinite(grad_rval[0]) || SCIPisInfinity(scip, REALABS(MIN(grad_sval[0],grad_rval[0]))) )
         {
            /* FIXME maybe it is sufficient to have one of them finite, using that one for the MIN below? */
            *success = FALSE;
            return SCIP_OKAY;
         }

         cutcoeff[0] = (yub-rval)*MIN(grad_sval[0],grad_rval[0]);
         cutcoeff[1] = (xub-sval)*MIN(grad_sval[0],grad_rval[0])+fsval-frval;
         cutcoeff[2] = xub-sval;
         cutcoeff[3] = cutcoeff[0]*sval+cutcoeff[1]*yub-cutcoeff[2]*fsval;
      }

      SCIPdebugMsg(scip, "UpperRight: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
      SCIPdebugMsg(scip, "UpperRight: r=%g in [%g,%g], s=%g in [%g,%g], f(r,yub)=%g, f(xub,s)=%g\n",rval,xlb,xub,sval,ylb,yub,frval,fsval);
      SCIPdebugMsg(scip, "(s,yub)=(%g,%g) (xub,r)=(%g,%g) t=%g, vredval=%g\n",sval,yub,xub,rval,tval,*convenvvalue);
      SCIPdebugMsg(scip, "UpperRight: cutcoeff[0]=%g, cutcoeff[1]=%g, cutcoeff[2]=%g, cutcoeff[3]=%g\n",cutcoeff[0],cutcoeff[1],cutcoeff[2],cutcoeff[3]);
   }

   return SCIP_OKAY;
}

/** generates a linear underestimator for f(x,y) with f(x,y) being convex in x and convex in y
 * generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that
 * alpha * x + beta * y - delta <= gamma * f(x,y)
 */
static
SCIP_RETCODE generateOrthogonal_lx_uy_Underestimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< function f(x,y) */
   SCIP_Real*            xyref,              /**< reference values for x and y */
   SCIP_Real             cutcoeff[4],        /**< cut coefficients alpha, beta, gamma, delta */
   SCIP_Real*            convenvvalue,       /**< function value of the convex envelope */
   SCIP_Bool*            success             /**< buffer to store whether coefficients were successfully computed */
   )
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_Real xval;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real yval;
   SCIP_Real ylb;
   SCIP_Real yub;
   SCIP_Real x0y0[2];

   SCIP_EXPR* vred;
   SCIP_EXPRTREE* vredtree;
   SCIP_EXPR* e1;
   SCIP_EXPR* e2;
   SCIP_EXPR* tmp;
   SCIP_EXPR* expr;
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;
   SCIP_EXPR* subst[2];

   SCIP_Real tval;
   SCIP_Real tlb;
   SCIP_Real tub;
   SCIP_Real sval;
   SCIP_Real rval;

   SCIP_Real frval;
   SCIP_Real fsval;
   SCIP_Real grad_rval[2];
   SCIP_Real grad_sval[2];

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f != NULL);
   assert(convenvvalue != NULL);
   assert(success != NULL);

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   xval = xyref[0];
   yval = xyref[1];

   /* check that variables are not unbounded or fixed and reference point is in interior */
   assert(!SCIPisInfinity(scip, -xlb));
   assert(!SCIPisInfinity(scip,  xub));
   assert(!SCIPisInfinity(scip, -ylb));
   assert(!SCIPisInfinity(scip,  yub));
   assert(!SCIPisEQ(scip,xlb,xub));
   assert(!SCIPisEQ(scip,ylb,yub));
   assert(!SCIPisEQ(scip,xlb,xval));
   assert(!SCIPisEQ(scip,xub,xval));
   assert(!SCIPisEQ(scip,ylb,yval));
   assert(!SCIPisEQ(scip,yub,yval));

   *success = FALSE;

   SCIPdebugMsg(scip, "f(%s, %s) = ", SCIPvarGetName(x), SCIPvarGetName(y));
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(f, SCIPgetMessagehdlr(scip), NULL) ) );
   SCIPdebugMsgPrint(scip, "\n");

   /* check in which triangle the point (xval,yval) lies */
   if( yval <= (yub-ylb)/(xub-xlb)*(xval-xlb)+ylb )
   {
      /* lower right triangle, i.e. region A_2 */
      /* construct v_red(t) := t f( xub+(xval-xub)/t, ylb ) + (1-t)*f( xub, (yval-ylb*t)/(1-t)) */

      /* construct e1:= f(xub+(xval-xub)/t, ylb) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );           /* expr = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, xval-xub) );     /* tmp  = xval-xub */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, tmp, expr) );      /* expr = (xval-xub)/t */
      if( xub != 0.0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, xub) );        /* tmp = xub */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PLUS, expr, tmp) );   /* expr = xub + (xval-xub)/t */
      }
      subst[0] = expr;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_CONST, ylb) );       /* subst[1] = ylb */

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e1, SCIPexprtreeGetRoot(f)) );         /* e1 = f(x,y) */
      assert(SCIPexprGetOperator(e1) != SCIP_EXPR_VARIDX);
      SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e1, subst) );                     /* e1 = f(xub + (xval-xub)/t, ylb) */

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);


      /* construct e2 := f(xub, (yval-t*ylb)/(1-t)) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_VARIDX, 0) );          /* expr1 = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, 1.0) );         /* tmp   = 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_MINUS, tmp, expr1) );  /* expr1 = 1-t */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_VARIDX, 0) );          /* expr2 = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, ylb) );         /* tmp   = ylb */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MUL, expr2, tmp) );    /* expr2 = ylb * t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, yval) );        /* tmp   = yval */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MINUS, tmp, expr2) );  /* expr2 = yval - ylb * t */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, expr2, expr1) );   /* expr =  (yval-t*ylb)/(1-t) */
      subst[1] = expr;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[0], SCIP_EXPR_CONST, xub) );      /* subst[0] = xub */

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e2, SCIPexprtreeGetRoot(f)) );        /* e2 = f(x,y) */
      assert(SCIPexprGetOperator(e2) != SCIP_EXPR_VARIDX);
      SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e2, subst) );                    /* e2 = f(xub, (yval-t*ylb)/(1-t)) */

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);


      /* construct vred := t * e1 + (1-t) * e2 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );           /* expr  = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_MUL, e1, expr) );      /* expr1 = t * e1*/


      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_VARIDX, 0) );          /* expr  = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, 1.0) );         /* tmp   = 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MINUS, tmp, expr) );   /* expr  = 1-t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MUL, e2, expr) );      /* expr2 = (1-t) * e2*/


      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &vred, SCIP_EXPR_PLUS, expr1, expr2) );
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &vredtree, vred, 1, 0, NULL) );
      SCIP_CALL( SCIPexprintCompile(exprinterpreter, vredtree) );


      /* compute bounds on t */
      tlb = (xub-xval)/(xub-xlb);
      tub = (yub-yval)/(yub-ylb);

      /* find t in [tlb, tub] such that vred'(t) = 0 */
      SCIP_CALL( solveDerivativeEquation(scip, exprinterpreter, vredtree, 0.0, tlb, tub, &tval, success) );

      if( *success == FALSE )
      {
         /* something went wrong when computing t */
         SCIP_CALL( SCIPexprtreeFree(&vredtree) );
         return SCIP_OKAY;
      }

      /* computing the cut coefficients */

      /* compute r and s from tval */
      rval = xub+(xval-xub)/tval;
      rval = MAX(xlb, MIN(xub, rval));
      sval = (yval-tval*ylb)/(1-tval);
      sval = MAX(ylb, MIN(yub, sval));

      /* compute vred(tval) */
      SCIP_CALL( SCIPexprtreeEval(vredtree, &tval, convenvvalue) );

      SCIP_CALL( SCIPexprtreeFree(&vredtree) );

      /* compute f(r, ylb) and f'(r, ylb) */
      x0y0[0] = rval;
      x0y0[1] = ylb;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &frval, grad_rval) );

      /* compute f(xub, s) and f'(xub,s) */
      x0y0[0] = xub;
      x0y0[1] = sval;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fsval, grad_sval) );

      /* generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that
       * alpha * x + beta * y - delta <= gamma * f(x,y) */
      if( !(SCIPisEQ(scip,rval,xlb)) )
      {
         /* take the slope along the x-axis and the slope between the points */
         if( !SCIPisFinite(grad_rval[0]) || SCIPisInfinity(scip, REALABS(grad_rval[0])) )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (sval-ylb)*grad_rval[0];
         cutcoeff[1] = (rval-xub)*grad_rval[0]-frval+fsval;
         cutcoeff[2] = sval-ylb;
         cutcoeff[3] = cutcoeff[0]*xub+cutcoeff[1]*sval-cutcoeff[2]*fsval;
      }
      else if( !(SCIPisEQ(scip,sval,yub)) )
      {
         /* take the slope along the y-axis and the slope between the points */
         if( !SCIPisFinite(grad_sval[1]) || SCIPisInfinity(scip, REALABS(grad_sval[1])) )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (ylb-sval)*grad_sval[1]-frval+fsval;
         cutcoeff[1] = (xub-rval)*grad_sval[1];
         cutcoeff[2] = xub-rval;
         cutcoeff[3] = cutcoeff[0]*xub+cutcoeff[1]*sval-cutcoeff[2]*fsval;
      }
      else
      {
         /* the point lies on the segment between (xlb,yub) and (xub,ylb) */
         if( !SCIPisFinite(grad_sval[0]) || !SCIPisFinite(grad_rval[0]) || SCIPisInfinity(scip, REALABS(MIN(grad_sval[0],grad_rval[0]))) )
         {
            /* FIXME maybe it is sufficient to have one of them finite, using that one for the MIN below? */
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (sval-ylb)*MIN(grad_sval[0],grad_rval[0]);
         cutcoeff[1] = (rval-xub)*MIN(grad_sval[0],grad_rval[0])+fsval-frval;
         cutcoeff[2] = sval-ylb;
         cutcoeff[3] = cutcoeff[0]*xub+cutcoeff[1]*sval-cutcoeff[2]*fsval;
      }


      SCIPdebugMsg(scip, "LowerRight: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
      SCIPdebugMsg(scip, "LowerRight: t=%g in [%g,%g], r=%g in [%g,%g], s=%g in [%g,%g]\n",tval,tlb,tub,rval,xlb,xub,sval,ylb,yub);
      SCIPdebugMsg(scip, "LowerRight: (r,ylb)=(%g,%g) (xub,sval)=(%g,%g) vredval=%g\n",rval,ylb,xub,sval,*convenvvalue);
      SCIPdebugMsg(scip, "LowerRight: cutcoeff[0]=%g, cutcoeff[1]=%g,cutcoeff[2]=1.0,cutcoeff[3]=%g\n",cutcoeff[0]/cutcoeff[2],cutcoeff[1]/cutcoeff[2],cutcoeff[3]/cutcoeff[2]);

   }
   else
   {
      /* (xval,yval) lie in the upper left triangle, i.e. region A_1 */
      /* construct v_red(t) := t f( xlb+(xval-xlb)/t, yub ) + (1-t)*f( xlb, (yval-yub*t)/(1-t) )  */

      /* construct e1:= f(xlb+(xval-xlb)/t, yub) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );           /* expr = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, xval-xlb) );     /* tmp  = xval-xlb */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, tmp, expr) );      /* expr = (xval-xlb)/lambda */
      if( xlb != 0.0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp, SCIP_EXPR_CONST, xlb) );        /* tmp = xlb */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_PLUS, expr, tmp) ); /* expr = xlb + (xval-xlb)/t */
      }
      subst[0] = expr;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_CONST, yub) );      /* subst[1] = yub */

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e1, SCIPexprtreeGetRoot(f)) );        /* e1 = f(x,y) */
      assert(SCIPexprGetOperator(e1) != SCIP_EXPR_VARIDX);
      SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e1, subst) );                    /* e1 = f(xlb + (xval-xlb)/t, yub) */

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);


      /* construct e2 := f(xlb, (yval-t*yub)/(1-t) ) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_VARIDX, 0) );          /* expr1 = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, 1.0) );         /* tmp   = 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_MINUS, tmp, expr1) );  /* expr1 = 1-t */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_VARIDX, 0) );          /* expr2 = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, yub) );         /* tmp   = yub */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MUL, expr2, tmp) );    /* expr2 = yub * t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,   SCIP_EXPR_CONST, yval) );        /* tmp   = yval */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MINUS, tmp, expr2) );  /* expr2 = yval - yub * t */

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, expr2, expr1) );   /* expr =  (yval-t*yub)/(1-t) */
      subst[1] = expr;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[0], SCIP_EXPR_CONST, xlb) );      /* subst[0] = xlb */

      SCIP_CALL( SCIPexprCopyDeep(SCIPblkmem(scip), &e2, SCIPexprtreeGetRoot(f)) );        /* e2 = f(x,y) */
      SCIP_CALL( SCIPexprSubstituteVars(SCIPblkmem(scip), e2, subst) );                    /* e2 = f( xlb , (yval-t*yub)/(1-t) ) */

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);


      /* construct vred := t * e1 + (1-t) * e2 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_VARIDX, 0) );          /* expr  = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr1, SCIP_EXPR_MUL, e1, expr) );      /* expr1 = t * e1*/


      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );           /* expr  = t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &tmp,  SCIP_EXPR_CONST, 1.0) );          /* tmp   = 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_MINUS, tmp, expr) );    /* expr  = 1-t */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_MUL, e2, expr) );      /* expr2 = (1-t) * e2*/


      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &vred, SCIP_EXPR_PLUS, expr1, expr2) );
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &vredtree, vred, 1, 0, NULL) );
      SCIP_CALL( SCIPexprintCompile(exprinterpreter, vredtree) );


      /* compute bounds on lambda */
      tlb = (xval-xlb)/(xub-xlb);
      tub = (yval-ylb)/(yub-ylb);

      /* find t in [tlb, tub] such that vred'(t) = 0 */
      SCIP_CALL( solveDerivativeEquation(scip, exprinterpreter, vredtree, 0.0, tlb, tub, &tval, success) );

      if( *success == FALSE )
      {
         /* something went wrong when computing s */
         SCIP_CALL( SCIPexprtreeFree(&vredtree) );
         return SCIP_OKAY;
      }

      /* computing the cut coefficients */

      /* compute r and s from tval */
      rval = xlb+(xval-xlb)/tval;
      rval = MAX(xlb, MIN(xub, rval));
      sval = (yval-tval*yub)/(1-tval);
      sval = MAX(ylb, MIN(yub, sval));

      /* compute vred(tval) */
      SCIP_CALL( SCIPexprtreeEval(vredtree, &tval, convenvvalue) );

      SCIP_CALL( SCIPexprtreeFree(&vredtree) );

      /* compute f(r, yub) and f'(r, yub) */
      x0y0[0] = rval;
      x0y0[1] = yub;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &frval, grad_rval) );

      /* compute f(xlb, s) and f'(xlb, s) */
      x0y0[0] = xlb;
      x0y0[1] = sval;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fsval, grad_sval) );

      /* generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that
       * alpha * x + beta * y - delta <= gamma * f(x,y) */
      if( !SCIPisEQ(scip,rval,xub) )
      {
         /* take the slope along the x-axis and the slope between the points */
         if( !SCIPisFinite(grad_rval[0]) || SCIPisInfinity(scip, REALABS(grad_rval[0])) )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (yub-sval)*grad_rval[0];
         cutcoeff[1] = (xlb-rval)*grad_rval[0]-fsval+frval;
         cutcoeff[2] = yub-sval;
         cutcoeff[3] = cutcoeff[0]*xlb+cutcoeff[1]*sval-cutcoeff[2]*fsval;
      }
      else if( !SCIPisEQ(scip,sval,ylb) )
      {
         /* take the slope along the y-axis and the slope between the points */
         if( !SCIPisFinite(grad_sval[1]) || SCIPisInfinity(scip, REALABS(grad_sval[1])) )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (sval-yub)*grad_sval[1]-fsval+frval;
         cutcoeff[1] = (rval-xlb)*grad_sval[1];
         cutcoeff[2] = rval-xlb;
         cutcoeff[3] = cutcoeff[0]*xlb+cutcoeff[1]*sval-cutcoeff[2]*fsval;
      }
      else
      {
         /* the point lies on the segment between (xlb,yub) and (xub,ylb) */
         if( !SCIPisFinite(grad_sval[0]) || !SCIPisFinite(grad_rval[0]) || SCIPisInfinity(scip, REALABS(MIN(grad_rval[0],grad_sval[0]))) )
         {
            /* FIXME maybe it is sufficient to have one of them finite, using that one for the MIN below? */
            *success = FALSE;
            return SCIP_OKAY;
         }
         cutcoeff[0] = (yub-sval)*MIN(grad_rval[0],grad_sval[0]);
         cutcoeff[1] = (xlb-rval)*MIN(grad_rval[0],grad_sval[0])-fsval+frval;
         cutcoeff[2] = yub-sval;
         cutcoeff[3] = cutcoeff[0]*xlb+cutcoeff[1]*sval-cutcoeff[2]*fsval;
      }

      SCIPdebugMsg(scip, "UpperLeft: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
      SCIPdebugMsg(scip, "UpperLeft: r=%g in [%g,%g], s=%g in [%g,%g], f(r,yub)=%g, f(xlb,s)=%g\n",rval,xlb,xub,sval,ylb,yub,frval,fsval);
      SCIPdebugMsg(scip, "t=%g in [%g,%g], (r,yub)=(%g,%g) (xlb,sval)=(%g,%g) vredval=%g\n",tval,tlb,tub,rval,yub,xlb,sval,*convenvvalue);
      SCIPdebugMsg(scip, "UpperLeft: cutcoeff[0]=%g, cutcoeff[1]=%g,cutcoeff[2]=1.0,cutcoeff[3]=%g\n",cutcoeff[0]/cutcoeff[2],cutcoeff[1]/cutcoeff[2],cutcoeff[3]/cutcoeff[2]);
   }

   return SCIP_OKAY;
}


/** generates a linear underestimator for f(x,y) with f(x,y) being STRICTLY convex in x and concave in y
 *  generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that alpha * x + beta * y - delta <= gamma * f(x,y)
 */
static
SCIP_RETCODE generateConvexConcaveUnderestimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< function f(x,y) */
   SCIP_EXPRTREE*        f_yfixed,           /**< function f(x;y) with x variable and y parameter */
   SCIP_EXPRTREE*        vred,               /**< function vred(s;x0,y0,ylb,yub) */
   SCIP_Real             xyref[2],           /**< reference values for (x,y) */
   SCIP_Real             cutcoeff[4],        /**< cut coefficients alpha, beta, gamma, delta */
   SCIP_Real*            convenvvalue,       /**< function value of the convex envelope */
   SCIP_Bool*            success             /**< buffer to store whether coefficients were successfully computed */
   )
{
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xval;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      yval;
   SCIP_Real      ylb;
   SCIP_Real      yub;

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f != NULL);
   assert(success != NULL);
   assert(xyref != NULL);

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   xval = xyref[0];
   yval = xyref[1];

   /* reference point should not be outside of bounds */
   assert(SCIPisLE(scip, xlb, xval));
   assert(SCIPisGE(scip, xub, xval));
   assert(SCIPisLE(scip, ylb, yval));
   assert(SCIPisGE(scip, yub, yval));

   *success = FALSE;

   if( SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub) )
   {
      SCIPdebugMsg(scip, "skip convex-concave underestimator, since y is unbounded\n");
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "f(%s, %s) = ", SCIPvarGetName(x), SCIPvarGetName(y));
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(f, SCIPgetMessagehdlr(scip), NULL) ) );
   SCIPdebugMsgPrint(scip, "\n");

   if( SCIPisEQ(scip, xlb, xub) )
   {
      /* x is fixed, so function is now concave -> generate secant between (x, ylb) and (x, yub) */
      SCIP_Real xy[2];
      SCIP_Real f_ylb;
      SCIP_Real f_yub;
      SCIP_Real slope;

      if( SCIPisEQ(scip, ylb, yub) )
      {
         SCIPdebugMsg(scip, "skip convex-concave underestimator, since both x and y are fixed\n");
         return SCIP_OKAY;
      }

      /* get f(x, ylb) */
      xy[0] = xlb;
      xy[1] = ylb;
      SCIP_CALL( SCIPexprintEval(exprinterpreter, f, xy, &f_ylb) );

      if( !SCIPisFinite(f_ylb) )
      {
         SCIPdebugMsg(scip, "cannot evaluate function at (xlb, ylb)\n");
         return SCIP_OKAY;
      }

      /* get f(x, yub) */
      xy[1] = yub;
      SCIP_CALL( SCIPexprintEval(exprinterpreter, f, xy, &f_yub) );

      if( !SCIPisFinite(f_yub) )
      {
         SCIPdebugMsg(scip, "cannot evaluate function at (xlb, yub)\n");
         return SCIP_OKAY;
      }

      slope = (f_yub - f_ylb) / (yub - ylb);

      /* secant is f(x,ylb) + slope * (y - ylb) <= f(x,y)*/

      cutcoeff[0]   = 0.0;                    /* coefficient of x == 0 */
      cutcoeff[1]   = slope;                  /* coefficient of y == slope */
      cutcoeff[2]   = 1.0;                    /* coefficient of f(x,y) == 1.0 */
      cutcoeff[3]   = -(f_ylb - slope * ylb); /* constant == -(f(x,ylb) - slope * ylb) */
      *convenvvalue = f_ylb+slope*(yval-ylb);

      *success = TRUE;
      return SCIP_OKAY;
   }

   if( SCIPisEQ(scip, ylb, yub) )
   {
      /* y is fixed, so function is now convex -> linearize in (xval, ylb) */
      SCIP_Real xy[2];
      SCIP_Real grad[2];
      SCIP_Real fval;

      xy[0] = xval;
      xy[1] = ylb;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, xy, TRUE, &fval, grad) );

      if( !SCIPisFinite(fval) || !SCIPisFinite(grad[0]) || SCIPisInfinity(scip, REALABS(grad[0])) )
      {
         perturb(&xval, xlb, xub, 0.001);
         xy[0] = xval;

         SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, xy, TRUE, &fval, grad) );

         if( !SCIPisFinite(fval) || !SCIPisFinite(grad[0]) || SCIPisInfinity(scip, REALABS(grad[0])) )
         {
            SCIPdebugMsg(scip, "cannot evaluate function or derivative in (xval,ylb), also after perturbation\n");
            return SCIP_OKAY;
         }
      }

      /* linearization is f(xval,ylb) + df/dx(xval,ylb) * (x - xval) <= f(x,y) */

      cutcoeff[0]   = grad[0];                  /* coefficient of x == gradient in x */
      cutcoeff[1]   = 0.0;                      /* coefficient of y == 0 */
      cutcoeff[2]   = 1.0;                      /* coefficient of f(x,y) == 1.0 */
      cutcoeff[3]   = -(fval - grad[0] * xval); /* constant == -(f(xval,ylb) - grad * xval) */
      *convenvvalue = fval;

      *success = TRUE;
      return SCIP_OKAY;
   }

   /* compute coefficients of a valid underestimating hyperplane */

   if( SCIPisFeasEQ(scip, xlb, xval) || SCIPisFeasEQ(scip, xub, xval) )
   {
      /* x is at it's lower or upper bound */
      SCIP_Real x0y0[2];
      SCIP_Real gradylb[2];
      SCIP_Real gradyub[2];
      SCIP_Real fvalylb;
      SCIP_Real fvalyub;

      xval = SCIPisFeasEQ(scip, xlb, xval) ? xlb : xub;

      /* compute f'(xval, ylb) and f'(xval, yub) */
      x0y0[0] = xval;
      x0y0[1] = ylb;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fvalylb, gradylb) );

      x0y0[1] = yub;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fvalyub, gradyub) );

      if( !SCIPisFinite(gradylb[0]) || !SCIPisFinite(gradyub[0]) || !SCIPisFinite(fvalylb) || !SCIPisFinite(fvalyub) ||
         SCIPisInfinity(scip, REALABS(gradylb[0])) || SCIPisInfinity(scip, REALABS(gradyub[0])) )
      {
         /* move xval inside domain and continue below, hope this will work better */
         perturb(&xval, xlb, xub, 0.001);
      }
      else
      {
         /* setup cut coefficients */
         if( xval == xlb ) /*lint !e777*/
            cutcoeff[0]   = (yub - ylb) * MIN(gradylb[0], gradyub[0]);/* coefficient of x */
         else
            cutcoeff[0]   = (yub - ylb) * MAX(gradylb[0], gradyub[0]);/* coefficient of x */
         cutcoeff[1]   = fvalyub - fvalylb;                           /* coefficient of y */
         cutcoeff[2]   = yub - ylb;                                   /* coefficient of f(x,y) */
         cutcoeff[3]   = cutcoeff[0] * xval + cutcoeff[1] * ylb - cutcoeff[2] * fvalylb;   /* constant */
         *convenvvalue = fvalylb;

         SCIPdebugMsg(scip, "alpha: %g, beta: %g, gamma: 1.0, delta: %g\n",
            cutcoeff[0]/cutcoeff[2], cutcoeff[1]/cutcoeff[2], cutcoeff[3]/cutcoeff[2]);

         *success = TRUE;
         return SCIP_OKAY;
      }
   }

   if( SCIPisFeasEQ(scip, ylb, yval) )
   {
      /* y is at it's lower bound */
      SCIP_Real x0y0[2];
      SCIP_Real grad[2];
      SCIP_Real xtilde;
      SCIP_Real fval, ftilde;

      /* these two cases should have been handled above */
      assert(!SCIPisEQ(scip, xlb, xval));
      assert(!SCIPisEQ(scip, xub, xval));

      assert(f_yfixed != NULL);

      /* compute f(xval, ylb) and f'(xval, ylb) */
      x0y0[0] = xval;
      x0y0[1] = ylb;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fval, grad) );

      if( !SCIPisFinite(fval) || !SCIPisFinite(grad[0]) || SCIPisInfinity(scip, REALABS(grad[0])) )
      {
         /* move yval inside domain and continue below, hope this will work better */
         perturb(&yval, ylb, yub, 0.001);
      }
      else
      {
         /* setup f(x,yub) */
         SCIPexprtreeSetParamVal(f_yfixed, 0, yub);
         SCIP_CALL( SCIPexprintNewParametrization(exprinterpreter, f_yfixed) );

         SCIPdebugMsg(scip, "f(x,yub) = ");
         SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(f_yfixed, SCIPgetMessagehdlr(scip), NULL) ) );
         SCIPdebugMsgPrint(scip, "\n");

         assert(SCIPexprtreeGetNVars(f_yfixed) == 1);

         /* find xtilde in [xlb, xub] such that f'(xtilde,yub) = f'(xval,ylb) */
         SCIP_CALL( solveDerivativeEquation(scip, exprinterpreter, f_yfixed, grad[0], xlb, xub, &xtilde, success) );

         if( !*success )
         {
            SCIP_Real fxlb;
            SCIP_Real fxub;

            /* if we could not find an xtilde such that f'(xtilde,yub) = f'(xval,ylb), then probably because f'(x,yub) is constant
             * in this case, choose xtilde from {xlb, xub} such that it maximizes f'(xtilde, yub) - grad[0]*xtilde
             */
            SCIP_CALL( SCIPexprintEval(exprinterpreter, f_yfixed, &xlb, &fxlb) ); /* coverity ignore ARRAY_VS_SINGLETON warning */
            SCIP_CALL( SCIPexprintEval(exprinterpreter, f_yfixed, &xub, &fxub) ); /* coverity ignore ARRAY_VS_SINGLETON warning */

            SCIPdebugMsg(scip, "couldn't solve deriv equ, compare f(%g,%g) - %g*%g = %g and f(%g,%g) - %g*%g = %g\n",
               xlb, ylb, grad[0], xlb, fxlb - grad[0] * xlb,
               xub, ylb, grad[0], xub, fxub - grad[0] * xub);

            if( SCIPisFinite(fxlb) && SCIPisFinite(fxub) )
            {
               if( fxlb - grad[0] * xlb > fxub - grad[0] * xub )
                  xtilde = xlb;
               else
                  xtilde = xub;
               *success = TRUE;
            }
            else
            {
               /* move yval inside domain and continue below, hope this will work better */
               perturb(&yval, ylb, yub, 0.001);
            }
         }

         if( *success )
         {
            /* compute f(xtilde, yub) */
            SCIP_CALL( SCIPexprintEval(exprinterpreter, f_yfixed, &xtilde, &ftilde) ); /* coverity ignore ARRAY_VS_SINGLETON warning */

            SCIPdebugMsg(scip, "xtilde = %g, f(%g,%g) = %g\n", xtilde, xtilde, yub, ftilde);

            /* setup cut coefficients */
            cutcoeff[0]   = (yub - ylb) * grad[0];                       /* coefficient of x */
            cutcoeff[1]   = ftilde - fval - grad[0] * (xtilde - xval);   /* coefficient of y */
            cutcoeff[2]   = yub - ylb;                                   /* coefficient of f(x,y) */
            cutcoeff[3]   = cutcoeff[0] * xval + cutcoeff[1] * ylb - cutcoeff[2] * fval;   /* constant */
            *convenvvalue = fval;

            SCIPdebugMsg(scip, "alpha: %g, beta: %g, gamma: %g, delta: %g\n", cutcoeff[0], cutcoeff[1], cutcoeff[2], cutcoeff[3]);

            return SCIP_OKAY;
         }
      }
   }

   if( SCIPisFeasEQ(scip, yval, yub) )
   {
      /* y is at it's upper bound */
      SCIP_Real x0y0[2];
      SCIP_Real grad[2];
      SCIP_Real fval;
      SCIP_Real xtilde;
      SCIP_Real ftilde;

      assert(f_yfixed != NULL);

      /* compute f(xval, yub) and f'(xval, yub) */
      x0y0[0] = xval;
      x0y0[1] = yub;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fval, grad) );

      if( !SCIPisFinite(fval) || !SCIPisFinite(grad[0]) || SCIPisInfinity(scip, REALABS(grad[0])) )
      {
         /* move yval inside domain and continue below, hope this will work better */
         perturb(&yval, ylb, yub, 0.001);
      }
      else
      {
         /* setup f(x,ylb) */
         SCIPexprtreeSetParamVal(f_yfixed, 0, ylb);
         SCIP_CALL( SCIPexprintNewParametrization(exprinterpreter, f_yfixed) );

         assert(SCIPexprtreeGetNVars(f_yfixed) == 1);

         /* find xtilde in [xlb, xub] such that f'(x,ylb) = f'(xval,yub) */
         SCIP_CALL( solveDerivativeEquation(scip, exprinterpreter, f_yfixed, grad[0], xlb, xub, &xtilde, success) );

         if( !*success )
         {
            SCIP_Real fxlb;
            SCIP_Real fxub;

            /* if we could not find an xtilde such that f'(xtilde,ylb) = f'(xval,yub), then probably because f'(x,ylb) is constant
             * in this case, choose xtilde from {xlb, xub} such that it maximizes f'(xtilde, yub) - grad[0]*xtilde
             */
            SCIP_CALL( SCIPexprintEval(exprinterpreter, f_yfixed, &xlb, &fxlb) ); /* coverity ignore ARRAY_VS_SINGLETON warning */
            SCIP_CALL( SCIPexprintEval(exprinterpreter, f_yfixed, &xub, &fxub) ); /* coverity ignore ARRAY_VS_SINGLETON warning */

            SCIPdebugMsg(scip, "couldn't solve deriv equ, compare f(%g,%g) - %g*%g = %g and f(%g,%g) - %g*%g = %g\n",
               xlb, yub, grad[0], xlb, fxlb - grad[0] * xlb,
               xub, yub, grad[0], xub, fxub - grad[0] * xub);

            if( SCIPisFinite(fxlb) && SCIPisFinite(fxub) )
            {
               if( fxlb - grad[0] * xlb < fxub - grad[0] * xub )
                  xtilde = xlb;
               else
                  xtilde = xub;
               *success = TRUE;
            }
            else
            {
               /* move yval inside domain and continue below, hope this will work better */
               perturb(&yval, ylb, yub, 0.001);
            }
         }

         if( *success )
         {
            /* compute f(xtilde, yub) */
            SCIP_CALL( SCIPexprintEval(exprinterpreter, f_yfixed, &xtilde, &ftilde) );  /* coverity ignore ARRAY_VS_SINGLETON warning */

            SCIPdebugMsg(scip, "xtilde = %g, f(%g,%g) = %g\n", xtilde, xtilde, ylb, ftilde);

            /* set up cut coefficients */
            cutcoeff[0]   = (yub - ylb) * grad[0];
            cutcoeff[1]   = grad[0] * (xtilde - xval) - ftilde + fval;
            cutcoeff[2]   = yub - ylb;
            cutcoeff[3]   = cutcoeff[0] * xval + cutcoeff[1] * yub - cutcoeff[2] * fval;
            *convenvvalue = fval;

            SCIPdebugMsg(scip, "alpha: %g, beta: %g, gamma: %g, delta: %g\n", cutcoeff[0], cutcoeff[1], cutcoeff[2], cutcoeff[3]);

            return SCIP_OKAY;
         }
      }
   }

   {
      /* x and y are somewhere between the bounds,
       * -> envelope is generated from f(x,y) in y=ylb and in y=yub
       */
      SCIP_Real paramvals[4];
#ifdef SCIP_DEBUG
      const char* paramnames[4] = {"x0", "y0", "ylb", "yub"};
#endif
      SCIP_Real t;
      SCIP_Real slb;
      SCIP_Real sub;
      SCIP_Real sval;
      SCIP_Real rval;
      SCIP_Real fsval;
      SCIP_Real frval;
      SCIP_Real grad[2];
      SCIP_Real x0y0[2];

      assert(vred != NULL);

      /* check that variables are not unbounded or fixed and reference point is in interior
       * @todo it should also work if x is unbounded, or? */
      /* assert(!SCIPisInfinity(scip, -xlb));
         assert(!SCIPisInfinity(scip,  xub)); */
      assert(!SCIPisInfinity(scip, -ylb));
      assert(!SCIPisInfinity(scip,  yub));

      /* update parameter values in vred */
      paramvals[0] = xval;
      paramvals[1] = yval;
      paramvals[2] = ylb;
      paramvals[3] = yub;
      SCIP_CALL( SCIPexprtreeSetParams(vred, 4, paramvals) );
      SCIP_CALL( SCIPexprintNewParametrization(exprinterpreter, vred) );

      SCIPdebugMsg(scip, "vred(s;x0,y0,ylb,yub) = ");
      SCIPdebug( SCIPexprtreePrint(vred, SCIPgetMessagehdlr(scip), NULL, NULL, paramnames) );
      SCIPdebugMsgPrint(scip, "\n");

      /* compute bounds on s */
      t = (yub - yval) / (yub - ylb);
      if( !SCIPisInfinity(scip, xub) )
         slb = (yval - yub) / (ylb - yval) * (xval / t - xub);
      else
         slb = -SCIPinfinity(scip);
      if( !SCIPisInfinity(scip, xlb) )
         sub = (yval - yub) / (ylb - yval) * (xval / t - xlb);
      else
         sub =  SCIPinfinity(scip);
      if( slb < xlb )
         slb = xlb;
      if( sub > xub )
         sub = xub;

      /* find s in [slb, sub] such that vred'(s) = 0 */
      SCIP_CALL( solveDerivativeEquation(scip, exprinterpreter, vred, 0.0, slb, sub, &sval, success) );
      assert(!*success || !SCIPisInfinity(scip, REALABS(sval)));

      if( *success )
      {
         /* compute r from s */
         rval = xval / t + (1.0 - 1.0 / t) * sval;
         assert(SCIPisFeasGE(scip, rval, xlb));
         assert(SCIPisFeasLE(scip, rval, xub));
         rval = MAX(xlb, MIN(rval, xub));

         /* compute f(sval, yub) */
         x0y0[0] = sval;
         x0y0[1] = yub;
         SCIP_CALL( SCIPexprtreeEval(f, x0y0, &fsval) );

         /* compute f(rval, ylb) */
         x0y0[0] = rval;
         x0y0[1] = ylb;
         SCIP_CALL( SCIPexprtreeEval(f, x0y0, &frval) );

         if( !SCIPisEQ(scip, sval, xlb) && !SCIPisEQ(scip, sval, xub) )
         {
            x0y0[0] = sval;
            x0y0[1] = yub;

            /* compute f'(xbar, ybar) */
            SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fsval, grad) );
         }
         else if( !SCIPisEQ(scip, rval, xlb) && !SCIPisEQ(scip, rval, xub) )
         {
            x0y0[0] = rval;
            x0y0[1] = ylb;

            /* compute f'(xbar, ybar) */
            SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &frval, grad) );
         }
         else
         {
            /* rare case
             * both points (sval, yub) and (rval, ylb) should yield valid inequality
             * for now, just take the first one, if differentiable, otherwise second one
             */
            x0y0[0] = sval;
            x0y0[1] = yub;

            /* compute f'(xbar, ybar) */
            SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fsval, grad) );

            if( !SCIPisFinite(grad[0]) )
            {
               x0y0[0] = rval;
               x0y0[1] = ylb;

               /* compute new f'(xbar, ybar) */
               SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &frval, grad) );
            }
         }

         /* compute vred(s) = t * f(rval, ylb) + (1-t) * f(sval, yub) */
         *convenvvalue = t * frval + (1.0 - t) * fsval;

         SCIPdebugMsg(scip, "Parallel: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
         SCIPdebugMsg(scip, "Parallel: r=%g s=%g in [%g,%g], y in [%g,%g], f(r,ylb)=%g, f(xlb,s)=%g\n",rval,sval,xlb,xub,ylb,yub,frval,fsval);
         SCIPdebugMsg(scip, "(r,ylb)=(%g,%g), (s,yub)=(%g,%g), vredval=%g\n",rval,ylb,sval,yub,*convenvvalue);

         if( !SCIPisFinite(grad[0]) || SCIPisInfinity(scip, REALABS(grad[0])) )
         {
            SCIPdebugMsg(scip, "f not differentiable at (x0,y0) w.r.t. x\n");
            *success = FALSE;
            return SCIP_OKAY;
         }

         /* compute cut coefficients */
         cutcoeff[0]   = (yub - ylb) * grad[0];
         cutcoeff[1]   = fsval - frval - (sval - rval) * grad[0];
         cutcoeff[2]   = yub - ylb;
         cutcoeff[3]   = cutcoeff[0] * xval + cutcoeff[1] * yval - cutcoeff[2] * *convenvvalue;

         SCIPdebugMsg(scip, "Parallel: cutcoeff[0]=%g, cutcoeff[1]=%g,cutcoeff[2]=1.0,cutcoeff[3]=%g\n",cutcoeff[0]/cutcoeff[2],cutcoeff[1]/cutcoeff[2],cutcoeff[3]/cutcoeff[2]);
      }
   }

   return SCIP_OKAY;
}


/** generates a cut for one side of lhs <= f(x,y) + c*z <= rhs with f(x,y) being convex in x and concave in y */
static
SCIP_RETCODE generateConvexConcaveEstimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             xyref[2],           /**< reference values for nonlinear variables */
   SCIP_SIDETYPE         violside,           /**< for which side of constraint to find a cut */
   SCIP_ROW**            row                 /**< storage for cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      cutcoeff[4];
   SCIP_Real      dummy;
   SCIP_Bool      success;
   SCIP_Real      coefs[2];
   char           cutname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->f != NULL);
   assert(consdata->convextype == SCIP_BIVAR_CONVEX_CONCAVE);

   *row = NULL;

   SCIPdebugMsg(scip, "generate %sestimator for convex-concave constraint <%s>\n",
      (violside == SCIP_SIDETYPE_LEFT ? "over" : "under"), SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   if( violside == SCIP_SIDETYPE_LEFT )
   {
      /* need overestimator */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      if( consdata->sepaconvexconcave.lineariny )
      {
         /* f is strictly convex in x and linear in y -> overestimator is polyhedral */
         SCIP_Real constant;

         SCIP_CALL( generateEstimatingHyperplane(scip, exprinterpreter, consdata->f, TRUE, xyref, &coefs[0], &coefs[1], &constant, &success) );

         if( success )
         {
            assert(SCIPisFinite(coefs[0]));
            assert(SCIPisFinite(coefs[1]));
            assert(SCIPisFinite(constant));

            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_overesthyperplanecut_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));
            SCIP_CALL( SCIPcreateRowCons(scip, row, SCIPconsGetHdlr(cons), cutname, 0, NULL, NULL, consdata->lhs - constant, SCIPinfinity(scip), TRUE, FALSE, TRUE) );

            SCIP_CALL( SCIPaddVarsToRow(scip, *row, 2, SCIPexprtreeGetVars(consdata->f), coefs) );
            if( consdata->z != NULL )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->z, consdata->zcoef) );
            }
         }
      }
      else
      {
         SCIP_Real xyref_[2];

         /* f is strictly concave in y -> can compute overestimator by applying generateConvexConcaveUnderstimator on -f(y,x) */
         assert(consdata->sepaconvexconcave.f_neg_swapped != NULL);

         xyref_[0] = xyref[1];
         xyref_[1] = xyref[0];
         SCIP_CALL( generateConvexConcaveUnderestimator(scip, exprinterpreter, consdata->sepaconvexconcave.f_neg_swapped, consdata->sepaconvexconcave.f_neg_swapped_yfixed, consdata->sepaconvexconcave.vred_neg_swapped, xyref_, cutcoeff, &dummy, &success) );

         if( success )
         {
            assert(SCIPisFinite(cutcoeff[0]));
            assert(SCIPisFinite(cutcoeff[1]));
            assert(SCIPisFinite(cutcoeff[2]));
            assert(SCIPisFinite(cutcoeff[3]));
            assert(SCIPisPositive(scip, cutcoeff[2])); /* assert gamma > 0 */

            /* construct row from cut coefficients (alpha, beta, gamma, delta)
             * coefficients are such that alpha * y + beta * x - gamma * (-f(x,y)) <= delta,
             * i.e., gamma * f(x,y) <= delta - alpha * y - beta * x
             * -> lhs <= f(x,y) + c*z <= delta/gamma - alpha/gamma * y - beta/gamma * x + c*z
             */
            coefs[0] = -cutcoeff[1] / cutcoeff[2];
            coefs[1] = -cutcoeff[0] / cutcoeff[2];
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_convexconcaveoverest_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, row, SCIPconsGetHdlr(cons), cutname, consdata->lhs - cutcoeff[3]/cutcoeff[2], SCIPinfinity(scip),
                  TRUE, FALSE /* modifiable */, TRUE /* removable */) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *row, 2, SCIPexprtreeGetVars(consdata->f), coefs) );
            if( consdata->z != NULL )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->z, consdata->zcoef) );
            }
         }
      }
   }
   else
   {
      /* need underestimator */
      assert(violside == SCIP_SIDETYPE_RIGHT);
      assert(!SCIPisInfinity(scip, consdata->rhs));

      if( consdata->sepaconvexconcave.linearinx )
      {
         /* f is linear in x and strictly concave in y -> underestimator is polyhedral  */
         SCIP_Real constant;

         SCIP_CALL( generateEstimatingHyperplane(scip, exprinterpreter, consdata->f, FALSE, xyref, &coefs[0], &coefs[1], &constant, &success) );

         if( success )
         {
            assert(SCIPisFinite(coefs[0]));
            assert(SCIPisFinite(coefs[1]));
            assert(SCIPisFinite(constant));

            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_underesthyperplanecut_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));
            SCIP_CALL( SCIPcreateRowCons(scip, row, SCIPconsGetHdlr(cons), cutname, 0, NULL, NULL, -SCIPinfinity(scip), consdata->rhs - constant, TRUE, FALSE, TRUE) );

            SCIP_CALL( SCIPaddVarsToRow(scip, *row, 2, SCIPexprtreeGetVars(consdata->f), coefs) );
            if( consdata->z != NULL )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->z, consdata->zcoef) );
            }
         }
      }
      else
      {
         /* f is strictly convex in x -> can compute underestimator by applying generateConvexConcaveUnderstimator */
         assert(!consdata->sepaconvexconcave.linearinx); /* generateConvexConcaveUnderestimator assumes that if f is strictly convex in x */

         SCIP_CALL( generateConvexConcaveUnderestimator(scip, exprinterpreter, consdata->f, consdata->sepaconvexconcave.f_yfixed, consdata->sepaconvexconcave.vred, xyref, cutcoeff, &dummy, &success) );

         if( success )
         {
            assert(SCIPisFinite(cutcoeff[0]));
            assert(SCIPisFinite(cutcoeff[1]));
            assert(SCIPisFinite(cutcoeff[2]));
            assert(SCIPisFinite(cutcoeff[3]));
            assert(SCIPisPositive(scip, cutcoeff[2])); /* assert gamma > 0 */

            /* construct row from cut coefficients (alpha, beta, gamma, delta)
             * coefficients are such that alpha * x + beta * y - gamma * f(x,y) <= delta,
             * i.e., alpha/gamma * x + beta/gamma * y - delta/gamma <= f(x,y)
             * -> alpha/gamma * x + beta/gamma * y - delta/gamma + c*z <= f(x,y) + c*z <= rhs
             */

            coefs[0] = cutcoeff[0] / cutcoeff[2];
            coefs[1] = cutcoeff[1] / cutcoeff[2];
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_convexconcaveunderest_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, row, SCIPconsGetHdlr(cons), cutname, -SCIPinfinity(scip), consdata->rhs + cutcoeff[3]/cutcoeff[2],
                  TRUE, FALSE /* modifiable */, TRUE /* removable */) );
            SCIP_CALL( SCIPaddVarsToRow(scip, *row, 2, SCIPexprtreeGetVars(consdata->f), coefs) );
            if( consdata->z != NULL )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->z, consdata->zcoef) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** computes an underestimating hyperplane for functions that are convex in x and y if the point to cut off lies on the boundary */
static
SCIP_RETCODE lifting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< function f(x,y) */
   SCIP_Real             xval,               /**< current x value */
   SCIP_Real             yval,               /**< current y value */
   SCIP_Real             xlb,                /**< lower bound x */
   SCIP_Real             xub,                /**< upper bound x */
   SCIP_Real             ylb,                /**< lower bound y */
   SCIP_Real             yub,                /**< upper bound y */
   int                   min_max,            /**< min=-1  max=1 */
   SCIP_Real             cutcoeff[4],        /**< returns the lifting coefficient*/
   SCIP_Real*            convenvvalue,       /**< value of the convex envelope at (xval,yval) */
   SCIP_Bool*            success             /**< buffer to indicate whether lifting was successful */
   )
{
   int idx; /* indicates which variable is at the boundary */

   SCIP_Real mu;
   SCIP_Real fval;
   SCIP_Real grad[2];

   SCIP_Real x0y0[2];
   SCIP_Real f_lb;
   SCIP_Real f_ub;
   SCIP_Real grad_lb[2];
   SCIP_Real grad_ub[2];

   assert(SCIPisEQ(scip,xlb,xub) || SCIPisEQ(scip,ylb,yub));
   assert(success != NULL);

   *success = FALSE;
   idx = SCIPisEQ(scip, xlb, xub) ? 0 : 1;

   /* determine mu
    * if f is bivariate quadratic then f_x(xlb,yval) is linear in yval
    * thus the minimum is attained at the lower or the upper bound
    */
   x0y0[0] = xlb;
   x0y0[1] = ylb;
   SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &f_lb, grad_lb) );
   if( !SCIPisFinite(grad_lb[idx]) )
      return SCIP_OKAY;

   x0y0[0] = xub;
   x0y0[1] = yub;
   SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &f_ub, grad_ub) );
   if( !SCIPisFinite(grad_ub[idx]) )
      return SCIP_OKAY;

   /* if min_max=-1 choose min( grad_lb[idx], grad_ub[idx] )
    * if min_max= 1 choose max( grad_lb[idx], grad_ub[idx] )
    */
   if( min_max * (grad_lb[idx] - grad_ub[idx]) >= 0 )
      mu = grad_lb[idx];
   else
      mu = grad_ub[idx];

   /* determine coefficients for the hyperplane */
   x0y0[0] = xval;
   x0y0[1] = yval;
   SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, x0y0, TRUE, &fval, grad) );

   if( idx == 0 )
   {
      if( !SCIPisFinite(grad[1]) || SCIPisInfinity(scip, REALABS(grad[1])) )
         return SCIP_OKAY;
      cutcoeff[0] = mu;
      cutcoeff[1] = grad[1];
   }
   else
   {
      assert(idx == 1);
      if( !SCIPisFinite(grad[0]) || SCIPisInfinity(scip, REALABS(grad[0])) )
         return SCIP_OKAY;
      cutcoeff[0] = grad[0];
      cutcoeff[1] = mu;
   }
   cutcoeff[2] = 1;
   cutcoeff[3] = -(fval-cutcoeff[0]*xval-cutcoeff[1]*yval);
   *convenvvalue = fval;
   *success = TRUE;

   return SCIP_OKAY;
}

/** generate a linear underestimator for f(x,y) with f(x,y) being convex in x and convex in y and the point to cut off lies on the boundary
 * generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that alpha * x + beta * y - delta <= gamma * f(x,y)
 */
static
SCIP_RETCODE generate1ConvexIndefiniteUnderestimatorAtBoundary(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< function f(x,y) */
   SCIP_Real             xyref[2],           /**< reference values for x and y */
   SCIP_Real             cutcoeff[4],        /**< cut coefficients alpha, beta, gamma, delta */
   SCIP_Real*            convenvvalue,       /**< function value of the convex envelope */
   SCIP_Bool*            success             /**< buffer to store whether coefficients were successfully computed */
   )
{
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_Real xval;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real yval;
   SCIP_Real ylb;
   SCIP_Real yub;

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f != NULL);
   assert(convenvvalue != NULL);
   assert(success != NULL);

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   *success = FALSE;

   SCIPdebugMsg(scip, "f(%s, %s) = ", SCIPvarGetName(x), SCIPvarGetName(y));
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(f, SCIPgetMessagehdlr(scip), NULL) ) );
   SCIPdebugMsgPrint(scip, "\n");

   xval = xyref[0];
   yval = xyref[1];

   SCIPdebugMsg(scip, "xval=%g in [%g,%g], yval=%g in [%g,%g]\n",xval,xlb,xub,yval,ylb,yub);

   if( SCIPisEQ(scip, ylb, yub) )
   {
      /* y is fixed, so function is now convex -> linearize in (xval, ylb) */
      SCIP_Real xy[2];
      SCIP_Real grad[2];
      SCIP_Real fval;

      xy[0] = xval;
      xy[1] = ylb;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, xy, TRUE, &fval, grad) );
      if( !SCIPisFinite(grad[0]) || SCIPisInfinity(scip, REALABS(grad[0])) )
         return SCIP_OKAY;

      /* linearization is f(xval,ylb) + df/dx(xval,ylb) * (x - xval) <= f(x,y) */

      cutcoeff[0] = grad[0];                  /* coefficient of x == gradient in x */
      cutcoeff[1] = 0.0;                      /* coefficient of y == 0 */
      cutcoeff[2] = 1.0;                      /* coefficient of f(x,y) == 1.0 */
      cutcoeff[3] = -(fval - grad[0] * xval); /* constant == -(f(xval,ylb) - grad * xval) */

      *success = TRUE;
      return SCIP_OKAY;
   }

   if( SCIPisEQ(scip, xlb, xub) )
   {
      /* x is fixed, so function is now convex -> linearize in (xlb, yval) */
      SCIP_Real xy[2];
      SCIP_Real grad[2];
      SCIP_Real fval;

      xy[0] = xlb;
      xy[1] = yval;
      SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, xy, TRUE, &fval, grad) );
      if( !SCIPisFinite(grad[1]) || SCIPisInfinity(scip, REALABS(grad[1])) )
         return SCIP_OKAY;

      /* linearization is f(xlb,yval) + df/dy(xlb,yval) * (y - yval) <= f(x,y) */

      cutcoeff[0] = 0.0;                      /* coefficient of x == 0.0 */
      cutcoeff[1] = grad[1];                  /* coefficient of y == gradient in y */
      cutcoeff[2] = 1.0;                      /* coefficient of f(x,y) == 1.0 */
      cutcoeff[3] = -(fval - grad[1] * yval); /* constant == -(f(xlb,yval) - grad * yval) */

      *success = TRUE;
      return SCIP_OKAY;
   }

   /* check if the points lie on a boundary */
   if( SCIPisFeasEQ(scip, xlb, xval) )
   {
      /* apply a lifting and exploit that the function is convex in x and y
       * Idea: f(xlb,y) + mu (x-xlb) <= f(x,y)
       * determine mu with mu <= min_{x,y} ( f(x,y)-f(xlb,y) ) / (x-xlb)
       * f is convex in x: mu<= min_{y} f_x(xlb,y)
       *
       * mu (x-lb) + f_y(xlb,yval) * y <= f(x,y)
       */
      xval = xlb;

      SCIP_CALL( lifting(scip,exprinterpreter,f,xval,yval,xlb,xlb,ylb,yub,-1,cutcoeff,convenvvalue,success) );

      if( !*success )
         return SCIP_OKAY;

      SCIPdebugMsg(scip, "Boundary x=lb: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
      SCIPdebugMsg(scip, "convenvvalue = %g\n",*convenvvalue);
      SCIPdebugMsg(scip, "cutcoeff[0]=%g, cutcoeff[1]=%g,cutcoeff[2]=%g,cutcoeff[3]=%g\n",
         cutcoeff[0],cutcoeff[1],cutcoeff[2],cutcoeff[3]);

      return SCIP_OKAY;
   }

   if( SCIPisFeasEQ(scip, ylb, yval) )
   {
      yval = ylb;

      SCIP_CALL( lifting(scip,exprinterpreter,f,xval,yval,xlb,xub,ylb,ylb,-1,cutcoeff,convenvvalue,success) );

      if( !*success )
         return SCIP_OKAY;

      SCIPdebugMsg(scip, "Boundary y=lb: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
      SCIPdebugMsg(scip, "convenvvalue = %g\n",*convenvvalue);
      SCIPdebugMsg(scip, "cutcoeff[0]=%g, cutcoeff[1]=%g,cutcoeff[2]=%g,cutcoeff[3]=%g\n",
         cutcoeff[0],cutcoeff[1],cutcoeff[2],cutcoeff[3]);

      return SCIP_OKAY;
   }

   if( SCIPisFeasEQ(scip, xub, xval) )
   {
      /* apply a lifting and exploit that the function is convex in x and y
       * Idea: f(xlb,y) + mu (xub-x) <= f(x,y)
       * determine mu with mu <= min_{x,y} ( f(x,y)-f(xub,y) ) / (xub-x)
       * f is convex in x: -1 * mu >= min_{y} f_x(xub,y)
       *
       * mu (xub-x)    + f_y(xub,yval) * y <= f(x,y)
       * -mu*x -mu*xub + f_y(xub,yval) * y <= f(x,y)
       */
      xval = xub;

      SCIP_CALL( lifting(scip,exprinterpreter,f,xval,yval,xub,xub,ylb,yub,1,cutcoeff,convenvvalue,success) );

      if( !*success )
         return SCIP_OKAY;

      SCIPdebugMsg(scip, "Boundary x=ub: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
      SCIPdebugMsg(scip, "convenvvalue = %g\n",*convenvvalue);
      SCIPdebugMsg(scip, "cutcoeff[0]=%g, cutcoeff[1]=%g,cutcoeff[2]=%g,cutcoeff[3]=%g\n",
         cutcoeff[0],cutcoeff[1],cutcoeff[2],cutcoeff[3]);

      return SCIP_OKAY;
   }

   if( SCIPisFeasEQ(scip, yub, yval) )
   {
      yval = yub;

      SCIP_CALL( lifting(scip,exprinterpreter,f,xval,yval,xlb,xub,yub,yub,1,cutcoeff,convenvvalue,success) );

      if( !*success )
         return SCIP_OKAY;

      SCIPdebugMsg(scip, "Boundary y=ub: Cut of (xval,yval)=(%g,%g)\n",xval,yval);
      SCIPdebugMsg(scip, "convenvvalue = %g\n",*convenvvalue);
      SCIPdebugMsg(scip, "cutcoeff[0]=%g, cutcoeff[1]=%g,cutcoeff[2]=%g,cutcoeff[3]=%g\n",
         cutcoeff[0],cutcoeff[1],cutcoeff[2],cutcoeff[3]);

      return SCIP_OKAY;
   }

   /* (xval,yval) lies in the interior */
   SCIPerrorMessage("Tries to compute underestimator for a point at the boundary. But point is not on the boundary!\n");
   return SCIP_ERROR;
}

/** generates a linear underestimator for f(x,y) with f(x,y) being convex in x and convex in y but indefinite
 *  This is for the case where the cone of the concave directions is (R_+ x R_-) union (R_\- x R_+).
 *  We consider two cases:
 *    a) the underestimating segmenent connects parallel facets
 *    b) the underestimating segmenent connects orthogonal facets where
 *       x=l_x, y=l_y and x=u_x, y=u_y
 *  We ensure that the parallel facets are the horizontal with y=l_y and y=u_y
 *  We compute the objective value of the two problems.
 *  The smaller objective value corresponds to the convex envelope.
 *  The supporting hyperplane is then constructed at the this point.
 */
static
SCIP_RETCODE generate1ConvexIndefiniteUnderestimatorInTheInteriorPatternA(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< function f(x,y) */
   SCIP_Real             xyref[2],           /**< reference values for x and y */
   SCIP_Real             cutcoeff[4],        /**< cut coefficients alpha, beta, gamma, delta */
   SCIP_Real*            convenvvalue,       /**< function value of the convex envelope */
   SCIP_Bool*            success             /**< buffer to store whether coefficients were successfully computed */
   )
{
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      ylb;
   SCIP_Real      yub;
   SCIP_Real      xub_ylb[2];
   SCIP_Real      xlb_yub[2];
   SCIP_Real      grad_xub_ylb[2];
   SCIP_Real      grad_xlb_yub[2];
   SCIP_Real      fval_xub_ylb;
   SCIP_Real      fval_xlb_yub;

   SCIP_Real      all_cutcoeff[2][4];
   SCIP_Real      all_convenvvalue[2];
   SCIP_Bool      all_success[2];

   SCIP_Real      lowest;
   int            lowestidx;
   int            i;

   SCIP_EXPRTREE* fswapped;
   SCIP_VAR*      vars[2];
   SCIP_Bool      swapped;
   SCIP_Real      swap_buffer;
   SCIP_EXPR*     subst[2];

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f != NULL);
   assert(convenvvalue != NULL);
   assert(success != NULL);

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   *success = FALSE;

   xub_ylb[0] = xub;
   xub_ylb[1] = ylb;
   xlb_yub[0] = xlb;
   xlb_yub[1] = yub;

   SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, xub_ylb, TRUE, &fval_xub_ylb, grad_xub_ylb) );
   SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, xlb_yub, TRUE, &fval_xlb_yub, grad_xlb_yub) );

   if( !SCIPisFinite(fval_xub_ylb) || SCIPisInfinity(scip, REALABS(fval_xub_ylb)) || !SCIPisFinite(fval_xlb_yub) || SCIPisInfinity(scip, REALABS(fval_xlb_yub)) )
   {
      SCIPdebugMsg(scip, "skip 1-convex underestimator since function cannot be evaluated\n");
      return SCIP_OKAY;
   }

   if( !SCIPisFinite(grad_xub_ylb[0]) || !SCIPisFinite(grad_xlb_yub[1]) )
   {
      SCIPdebugMsg(scip, "skip 1-convex underestimator since function cannot be differentiated\n");
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "f(%s, %s) = ", SCIPvarGetName(x), SCIPvarGetName(y));
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(f, SCIPgetMessagehdlr(scip), NULL) ) );
   SCIPdebugMsgPrint(scip, "\n");

   SCIPdebugMsg(scip, "xval=%g in [%g,%g], yval=%g in [%g,%g]\n", xyref[0], xlb, xub, xyref[1], ylb, yub);

   /* assure (xub-xlb)*f_x(xub,ylb) - (yub-ylb)*f_y(xlb,yub) >= f(xub,ylb) - f(xlb,yub) */
   /* f_y(xlb,yub)*(ylb-yub)* + f(xlb,yub) >= f_x(xub,ylb)*(xub-xlb) + f(xub,ylb) */
   if( fval_xub_ylb-fval_xlb_yub <= (xub-xlb)*grad_xub_ylb[0]-(yub-ylb)*grad_xlb_yub[1] )
   {
      swapped = 0;
   }
   else
   {
      /* swap the variables */
      swapped = 1;

      vars[0] = SCIPexprtreeGetVars(f)[1];
      vars[1] = SCIPexprtreeGetVars(f)[0];

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[0], SCIP_EXPR_VARIDX, 1) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_VARIDX, 0) );

      SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &fswapped, f) );
      SCIP_CALL( SCIPexprtreeSubstituteVars(fswapped, subst) );
      SCIP_CALL( SCIPexprtreeSetVars(fswapped, 2, vars) );
      SCIP_CALL( SCIPexprintCompile(exprinterpreter, fswapped) );

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);
   }

   if( swapped == 0 )
   {
      /* assume (xval,yval) lie in A1 (lower left triangle) or A2 (upper right triangle) */
      SCIP_CALL( generateOrthogonal_lx_ly_Underestimator(scip, exprinterpreter, f, xyref, all_cutcoeff[0], &all_convenvvalue[0], &all_success[0]) );
      /* assume (xval,yval) lie in A3 */
      SCIP_CALL( generateUnderestimatorParallelYFacets(scip, exprinterpreter, f, xyref, all_cutcoeff[1], &all_convenvvalue[1], &all_success[1]) );
   }
   else
   {
      SCIP_Real xyref_[2];

      assert(swapped == 1);

      xyref_[0] = xyref[1];
      xyref_[1] = xyref[0];

      /* assume (xval,yval) lie in A1 (lower left triangle) or A2 (upper right triangle) */
      SCIP_CALL( generateOrthogonal_lx_ly_Underestimator(scip, exprinterpreter, fswapped, xyref_, all_cutcoeff[0], &all_convenvvalue[0], &all_success[0]) );  /*lint !e644*/
      /* assume (xval,yval) lie in A3 */
      SCIP_CALL( generateUnderestimatorParallelYFacets(scip, exprinterpreter, fswapped, xyref_, all_cutcoeff[1], &all_convenvvalue[1], &all_success[1]) );

      /* swap back */
      swap_buffer        = all_cutcoeff[0][0];
      all_cutcoeff[0][0] = all_cutcoeff[0][1];
      all_cutcoeff[0][1] = swap_buffer;

      swap_buffer        = all_cutcoeff[1][0];
      all_cutcoeff[1][0] = all_cutcoeff[1][1];
      all_cutcoeff[1][1] = swap_buffer;

      SCIP_CALL( SCIPexprtreeFree(&fswapped) );
   }

   /* Select the underestimator with the lowest convex envelope */
   SCIPdebugMsg(scip, "\n");
   SCIPdebugMsg(scip, "Triangulation: convenvvalue=%g\n", all_convenvvalue[0]);
   SCIPdebugMsg(scip, "Parallel    Y: convenvvalue=%g\n", all_convenvvalue[1]);

   lowest = SCIPinfinity(scip);
   lowestidx = -1;

   if( all_success[0] && all_success[1] )
   {
      *success = TRUE;
      for( i = 0; i < 2; ++i )
      {
         assert(SCIPisFinite(all_cutcoeff[i][0]));
         assert(SCIPisFinite(all_cutcoeff[i][1]));
         assert(SCIPisFinite(all_cutcoeff[i][2]));
         assert(SCIPisFinite(all_cutcoeff[i][3]));

         if( all_convenvvalue[i] < lowest )
         {
            /* if all_convenvvalue[0] == all_convenvalue[1], take all_convenvvalue[0] */
            lowest = all_convenvvalue[i];
            lowestidx = i;
         }
      }
      assert(lowestidx >= 0);

      *convenvvalue = all_convenvvalue[lowestidx];
      cutcoeff[0] = all_cutcoeff[lowestidx][0];
      cutcoeff[1] = all_cutcoeff[lowestidx][1];
      cutcoeff[2] = all_cutcoeff[lowestidx][2];
      cutcoeff[3] = all_cutcoeff[lowestidx][3];
      assert(SCIPisPositive(scip, cutcoeff[2])); /* assert gamma > 0 */
   }
   else
   {
      *success = FALSE;
   }

   return SCIP_OKAY;
}


/** generates a linear underestimator for f(x,y) with f(x,y) being convex in x and convex in y but indefinite
 *  This is for the case where the cone of the concave directions is (R_+ x R_+) union (R_- x R_-).
 *  We consider two cases:
 *    a) the underestimating segmenent connects parallel facets
 *    b) the underestimating segmenent connects orthogonal facets where
 *       x=l_x, y=u_y and x=u_x, y=l_y
 *  We ensure that the parallel facets are the horizontal with y=l_y and y=u_y
 *  We compute the objective value of the two problems.
 *  The smaller objective value corresponds to the convex envelope.
 *  The supporting hyperplane is then constructed at the this point.
 *  Generates coefficients cutcoeff = (alpha, beta, gamma, delta), such that alpha * x + beta * y - delta <= gamma * f(x,y)
 */
static
SCIP_RETCODE generate1ConvexIndefiniteUnderestimatorInTheInteriorPatternB(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_EXPRTREE*        f,                  /**< function f(x,y) */
   SCIP_Real             xyref[2],           /**< reference values for x and y */
   SCIP_Real             cutcoeff[4],        /**< cut coefficients alpha, beta, gamma, delta */
   SCIP_Real*            convenvvalue,       /**< function value of the convex envelope */
   SCIP_Bool*            success             /**< buffer to store whether coefficients were successfully computed */
   )
{
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      ylb;
   SCIP_Real      yub;
   SCIP_Real      xlb_ylb[2];
   SCIP_Real      xub_yub[2];
   SCIP_Real      grad_xlb_ylb[2];
   SCIP_Real      grad_xub_yub[2];
   SCIP_Real      fval_xlb_ylb;
   SCIP_Real      fval_xub_yub;

   SCIP_Real      all_cutcoeff[2][4];
   SCIP_Real      all_convenvvalue[2];
   SCIP_Bool      all_success[2];

   SCIP_Real      lowest;
   int            lowestidx;
   int            i;

   SCIP_EXPRTREE* fswapped;
   SCIP_VAR*      vars[2];
   SCIP_Bool      swapped;
   SCIP_Real      swap_buffer;
   SCIP_EXPR*     subst[2];

   assert(scip != NULL);
   assert(exprinterpreter != NULL);
   assert(f != NULL);
   assert(convenvvalue != NULL);
   assert(success != NULL);

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   *success = FALSE;

   SCIPdebugMsg(scip, "f(%s, %s) = ", SCIPvarGetName(x), SCIPvarGetName(y));
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(f, SCIPgetMessagehdlr(scip), NULL) ) );
   SCIPdebugMsgPrint(scip, "\n");

   xlb_ylb[0] = xlb;
   xlb_ylb[1] = ylb;
   xub_yub[0] = xub;
   xub_yub[1] = yub;

   SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, xlb_ylb, TRUE, &fval_xlb_ylb, grad_xlb_ylb) );
   SCIP_CALL( SCIPexprintGrad(exprinterpreter, f, xub_yub, TRUE, &fval_xub_yub, grad_xub_yub) );

   if( !SCIPisFinite(fval_xlb_ylb) || SCIPisInfinity(scip, REALABS(fval_xlb_ylb)) || !SCIPisFinite(fval_xub_yub) || SCIPisInfinity(scip, REALABS(fval_xub_yub)) )
   {
      SCIPdebugMsg(scip, "skip 1-convex underestimator since function cannot be evaluated\n");
      return SCIP_OKAY;
   }

   if( !SCIPisFinite(grad_xlb_ylb[1]) || !SCIPisFinite(grad_xub_yub[0]) )
   {
      SCIPdebugMsg(scip, "skip 1-convex underestimator since function cannot be differentiated\n");
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "xval=%g in [%g,%g], yval=%g in [%g,%g]\n",xyref[0],xlb,xub,xyref[1],ylb,yub);

   /* assure f_y(xlb,ylb)*(yub-ylb)* + f(xlb,ylb) >= f_x(xub,yub)*(xlb-xub) + f(xub,yub) */
   if( SCIPisGE( scip, fval_xlb_ylb+(yub-ylb)*grad_xlb_ylb[1], fval_xub_yub+(xlb-xub)*grad_xub_yub[0] ) )
   {
      swapped = 0;
   }
   else
   {
      /* swap the variables */
      swapped = 1;

      vars[0] = SCIPexprtreeGetVars(f)[1];
      vars[1] = SCIPexprtreeGetVars(f)[0];

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[0], SCIP_EXPR_VARIDX, 1) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &subst[1], SCIP_EXPR_VARIDX, 0) );

      SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &fswapped, f) );
      SCIP_CALL( SCIPexprtreeSubstituteVars(fswapped, subst) );
      SCIP_CALL( SCIPexprtreeSetVars(fswapped, 2, vars) );
      SCIP_CALL( SCIPexprintCompile(exprinterpreter, fswapped) );

      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[0]);
      SCIPexprFreeDeep(SCIPblkmem(scip), &subst[1]);
   }

   if( swapped == 0 )
   {
      /* assume (xval,yval) lie in A1 (lower left triangle) or A2 (upper right triangle) */
      SCIP_CALL( generateOrthogonal_lx_uy_Underestimator(scip, exprinterpreter, f, xyref, all_cutcoeff[0], &all_convenvvalue[0], &all_success[0]) );
      /* assume (xval,yval) lie in A3*/
      SCIP_CALL( generateUnderestimatorParallelYFacets(scip, exprinterpreter, f, xyref, all_cutcoeff[1], &all_convenvvalue[1], &all_success[1]) );
   }
   else
   {
      SCIP_Real xyref_[2];

      assert(swapped == 1);

      xyref_[0] = xyref[1];
      xyref_[1] = xyref[0];
      /* assume (xval,yval) lie in A1 (upper left triangle) or A2 (lower left triangle) */
      SCIP_CALL( generateOrthogonal_lx_uy_Underestimator(scip, exprinterpreter, fswapped, xyref_, all_cutcoeff[0], &all_convenvvalue[0], &all_success[0]) );  /*lint !e644*/
      /* assume (xval,yval) lie in A3 */
      SCIP_CALL( generateUnderestimatorParallelYFacets(scip, exprinterpreter, fswapped, xyref_, all_cutcoeff[1], &all_convenvvalue[1], &all_success[1]) );

      /* swap back */
      swap_buffer        = all_cutcoeff[0][0];
      all_cutcoeff[0][0] = all_cutcoeff[0][1];
      all_cutcoeff[0][1] = swap_buffer;

      swap_buffer        = all_cutcoeff[1][0];
      all_cutcoeff[1][0] = all_cutcoeff[1][1];
      all_cutcoeff[1][1] = swap_buffer;

      SCIP_CALL( SCIPexprtreeFree(&fswapped) );
   }

   /* select the underestimator with the lowest convex envelope */
   SCIPdebugMsg(scip, "\n");
   SCIPdebugMsg(scip, "Triangulation: convenvvalue=%g\n", all_convenvvalue[0]);
   SCIPdebugMsg(scip, "Parallel    Y: convenvvalue=%g\n", all_convenvvalue[1]);

   lowest = SCIPinfinity(scip);
   lowestidx = -1;

   if( all_success[0] && all_success[1] )
   {
      *success = TRUE;
      for( i = 0; i < 2; ++i )
      {
         assert(SCIPisFinite(all_cutcoeff[i][0]));
         assert(SCIPisFinite(all_cutcoeff[i][1]));
         assert(SCIPisFinite(all_cutcoeff[i][2]));
         assert(SCIPisFinite(all_cutcoeff[i][3]));

         /* if all_convenvvalue[0]==all_convenvalue[1], take all_convenvvalue[0] */
         if( all_convenvvalue[i] < lowest )
         {
            lowest = all_convenvvalue[i];
            lowestidx = i;
         }
      }
      assert(lowestidx >= 0);

      *convenvvalue = all_convenvvalue[lowestidx];
      cutcoeff[0] = all_cutcoeff[lowestidx][0];
      cutcoeff[1] = all_cutcoeff[lowestidx][1];
      cutcoeff[2] = all_cutcoeff[lowestidx][2];
      cutcoeff[3] = all_cutcoeff[lowestidx][3];
      assert(SCIPisPositive(scip, cutcoeff[2])); /* assert gamma > 0 */
   }
   else
   {
      *success = FALSE;
   }

   return SCIP_OKAY;
}


/** generates a linear underestimator for f(x,y) with f(x,y) being convex in x and convex in y but indefinite
 * generate coefficients cutcoeff = (alpha, beta, gamma, delta), such that alpha * x + beta * y - delta <= gamma * f(x,y)
 * 1. If the point lies on the boundary we apply the lifting technique.
 * 2. If the point lies in the interior we check the pattern of
 *    the concave directions and compute the corresponding underestimators.
 */
static
SCIP_RETCODE generate1ConvexIndefiniteUnderestimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            xyref,              /**< reference values for x and y */
   SCIP_ROW**            row                 /**< storage for cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* f;
   SCIP_Real      cutcoeff[4];
   SCIP_Bool      success;
   SCIP_Real      rhs;
   SCIP_Real      convenvvalue;

   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      ylb;
   SCIP_Real      yub;
   SCIP_Real      xy_mid[2];
   SCIP_Real      fval_mid;
   SCIP_Real      hess[4];

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(consdata->convextype == SCIP_BIVAR_1CONVEX_INDEFINITE);

   assert(!SCIPisInfinity(scip, consdata->rhs));

   f = consdata->f;

   x = SCIPexprtreeGetVars(f)[0];
   y = SCIPexprtreeGetVars(f)[1];

   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   xy_mid[0] = 0.5 * (xlb+xub);
   xy_mid[1] = 0.5 * (ylb+yub);

   /* assert that the bounds are finite */
   if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) || SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub) )
   {
      SCIPdebugMsg(scip, "skip underestimate for 1-convex indefinite constraint <%s> since <%s> or <%s> is unbounded\n", SCIPconsGetName(cons), SCIPvarGetName(x), SCIPvarGetName(y));
      return SCIP_OKAY;
   }

   success = FALSE;
   cutcoeff[0] = SCIP_INVALID;
   cutcoeff[1] = SCIP_INVALID;
   cutcoeff[2] = SCIP_INVALID;
   cutcoeff[3] = SCIP_INVALID;

   /* (xval,yval) lie on a boundary */
   if( SCIPisFeasEQ(scip,xyref[0],xlb) || SCIPisFeasEQ(scip,xyref[0],xub) || SCIPisFeasEQ(scip,xyref[1],ylb) || SCIPisFeasEQ(scip,xyref[1],yub) )
   {
      SCIP_CALL( generate1ConvexIndefiniteUnderestimatorAtBoundary(scip, exprinterpreter, f, xyref, cutcoeff, &convenvvalue, &success) );

      if( !success )
      {
         /* maybe f is not differentiable on boundary, so move reference point into interior
          * we do this here w.r.t. both coordinates
          */
         perturb(&xyref[0], xlb, xub, 0.001);
         perturb(&xyref[1], ylb, yub, 0.001);
      }
   }

   if( !success )
   {
      /* xyref lies in the interior */
      /* check the pattern of the concave directions */
      SCIP_CALL( SCIPexprintHessianDense(exprinterpreter, f, xy_mid, TRUE, &fval_mid, hess) );
      assert(SCIPisFinite(hess[1]));

      if( hess[1] > 0.0 )
      {
         /* Pattern A: (R>=0 x R<=0) union (R<=0 x R>=0)*/
         SCIPdebugMsg(scip, "Pattern A\n");
         SCIP_CALL( generate1ConvexIndefiniteUnderestimatorInTheInteriorPatternA(scip, exprinterpreter, f, xyref, cutcoeff, &convenvvalue, &success) );
      }
      else
      {
         /* Pattern B: (R>=0 x R>=0) union (R<=0 x R <=0)*/
         SCIPdebugMsg(scip, "Pattern B\n");
         SCIP_CALL( generate1ConvexIndefiniteUnderestimatorInTheInteriorPatternB(scip, exprinterpreter, f, xyref, cutcoeff, &convenvvalue, &success) );
      }
   }

   if( !success )
   {
      /* bad luck */
      *row = NULL;
      return SCIP_OKAY;
   }


   /* construct row from cut coefficients (alpha, beta, gamma, delta)
    * coefficients are such that alpha * x + beta * y - gamma * f(x,y) <= delta,
    * i.e., alpha/gamma * x + beta/gamma * y - delta/gamma <= f(x,y)
    * -> alpha/gamma * x + beta/gamma * y - delta/gamma + c*z <= f(x,y) + c*z <= rhs
    */

   assert(cutcoeff[0] != SCIP_INVALID); /*lint !e777*/
   assert(cutcoeff[1] != SCIP_INVALID); /*lint !e777*/
   assert(cutcoeff[2] != SCIP_INVALID); /*lint !e777*/
   assert(cutcoeff[3] != SCIP_INVALID); /*lint !e777*/
   assert(SCIPisFinite(cutcoeff[0]));
   assert(SCIPisFinite(cutcoeff[1]));
   assert(SCIPisFinite(cutcoeff[2]));
   assert(SCIPisFinite(cutcoeff[3]));
   assert(SCIPisPositive(scip, cutcoeff[2])); /* assert gamma > 0 */

   if( SCIPisInfinity(scip, REALABS(cutcoeff[0]/cutcoeff[2])) ||
      SCIPisInfinity( scip, REALABS(cutcoeff[1]/cutcoeff[2])) ||
      SCIPisInfinity( scip, REALABS(cutcoeff[3]/cutcoeff[2])) )
   {
      *row = NULL;
      return SCIP_OKAY;
   }

   rhs = consdata->rhs + cutcoeff[3]/cutcoeff[2];
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, row, SCIPconsGetHdlr(cons), "1ConvexUnderest", -SCIPinfinity(scip), rhs,
         TRUE, FALSE /* modifiable */, TRUE /* removable */) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, SCIPexprtreeGetVars(consdata->f)[0], cutcoeff[0] / cutcoeff[2]) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, SCIPexprtreeGetVars(consdata->f)[1], cutcoeff[1] / cutcoeff[2]) );
   if( consdata->z != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->z, consdata->zcoef) );
   }

   return SCIP_OKAY;
}

/** generates a cut */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprinterpreter,    /**< expressions interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_SIDETYPE         violside,           /**< for which side of constraint we want to generate a cut */
   SCIP_Real             cutmaxrange,        /**< bound on cut coef range */
   SCIP_ROW**            row                 /**< storage for cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_Real x0y0[2];

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *row = NULL;

   x = SCIPexprtreeGetVars(consdata->f)[0];
   y = SCIPexprtreeGetVars(consdata->f)[1];

   x0y0[0] = SCIPgetSolVal(scip, sol, x);
   x0y0[1] = SCIPgetSolVal(scip, sol, y);

   assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(x), x0y0[0]));
   assert(SCIPisFeasGE(scip, SCIPvarGetUbLocal(x), x0y0[0]));
   assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(y), x0y0[1]));
   assert(SCIPisFeasGE(scip, SCIPvarGetUbLocal(y), x0y0[1]));

   /* project into box */
   x0y0[0] = MIN(MAX(SCIPvarGetLbLocal(x),x0y0[0]),SCIPvarGetUbLocal(x));  /*lint !e666*/
   x0y0[1] = MIN(MAX(SCIPvarGetLbLocal(y),x0y0[1]),SCIPvarGetUbLocal(y));  /*lint !e666*/

   SCIPdebugMsgPrint(scip, "\n");
   SCIPdebugMsg(scip, "generate cut for constraint <%s> with %s hand side violated by %g\n", SCIPconsGetName(cons), violside == SCIP_SIDETYPE_LEFT ? "left" : "right", violside == SCIP_SIDETYPE_LEFT ? consdata->lhsviol : consdata->rhsviol);
   SCIPdebugMsg(scip, "convextype = %d\n",consdata->convextype);
   SCIPdebugMsg(scip, "%s = %g with bounds [%g, %g], %s = %g with bounds [%g, %g]",
      SCIPvarGetName(x), SCIPgetSolVal(scip, sol, x), SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x),
      SCIPvarGetName(y), SCIPgetSolVal(scip, sol, y), SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y));
   if( consdata->z != NULL )
      SCIPdebugMsgPrint(scip, ", %s = %g with bounds [%g, %g]", SCIPvarGetName(consdata->z), SCIPgetSolVal(scip, sol, consdata->z), SCIPvarGetLbLocal(consdata->z), SCIPvarGetUbLocal(consdata->z));
   SCIPdebugMsgPrint(scip, "\n");
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIPdebugMsgPrint(scip, "\n");

   switch( consdata->convextype )
   {
   case SCIP_BIVAR_ALLCONVEX:
   {
      if( violside == SCIP_SIDETYPE_RIGHT )
      {
         /* rhs is violated */
         SCIP_CALL( generateLinearizationCut(scip, exprinterpreter, cons, x0y0, FALSE, row) );
      }
      else
      {
         /* lhs is violated */
         SCIP_CALL( generateOverestimatingHyperplaneCut(scip, exprinterpreter, cons, x0y0, row) );
      }

      break;
   }

   case SCIP_BIVAR_CONVEX_CONCAVE:
   {
      SCIP_CALL( generateConvexConcaveEstimator(scip, exprinterpreter, cons, x0y0, violside, row) );
      break;
   }

   case SCIP_BIVAR_1CONVEX_INDEFINITE:
   {
      if( violside == SCIP_SIDETYPE_RIGHT )
      {
         /* rhs is violated */
         SCIP_CALL( generate1ConvexIndefiniteUnderestimator(scip, exprinterpreter, cons, x0y0, row) );
      }
      else
      {
         /* lhs is violated */
         SCIP_CALL( generateOverestimatingHyperplaneCut(scip, exprinterpreter, cons, x0y0, row) );
      }
      break;
   }
   default:
   {
      SCIPdebugMsg(scip, "cut generation for convexity type not implemented\n");
   }
   }  /*lint !e788*/

   if( *row == NULL )
      return SCIP_OKAY;

   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *row, NULL) ) );

   /* check numerics */
   {
      SCIP_Real mincoef;
      SCIP_Real maxcoef;

      mincoef = SCIPgetRowMinCoef(scip, *row);
      maxcoef = SCIPgetRowMaxCoef(scip, *row);

      while( maxcoef / mincoef > cutmaxrange )
      {
         SCIP_VAR* var;
         SCIP_Real coef;
         SCIP_Real constant;
         int j;

         /* if range of coefficients is bad, find very small coefficients and make them zero */
         SCIPdebugMsg(scip, "cut coefficients for constraint <%s> have very large range: mincoef = %g maxcoef = %g\n", SCIPconsGetName(cons), mincoef, maxcoef);

         /* if minimal coefficient is given by z, then give up (probably the maximal coefficient is the problem) */
         if( mincoef == consdata->zcoef )  /*lint !e777*/
         {
            SCIPdebugMsg(scip, "could not eliminate small coefficient, since it comes from linear part\n");
            break;
         }

         constant = 0.0;
         for( j = 0; j < SCIProwGetNNonz(*row); ++j )
         {
            coef = SCIProwGetVals(*row)[j];
            if( !SCIPisEQ(scip, REALABS(coef), mincoef) )
               continue;

            var = SCIPcolGetVar(SCIProwGetCols(*row)[j]);
            assert(var != NULL);

            /* try to eliminate coefficient with minimal absolute value by weakening cut and try again */
            if( ((coef > 0.0 && violside == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && violside == SCIP_SIDETYPE_LEFT)) && !SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
            {
               SCIPdebugMsg(scip, "eliminate coefficient %g for <%s> = %g [%g, %g]\n", coef, SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

               constant += coef * (SCIProwIsLocal(*row) ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var));
               SCIP_CALL( SCIPaddVarToRow(scip, *row, var, -coef) );
               continue;
            }

            if( ((coef < 0.0 && violside == SCIP_SIDETYPE_RIGHT) || (coef > 0.0 && violside == SCIP_SIDETYPE_LEFT)) && !SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
            {
               SCIPdebugMsg(scip, "eliminate coefficient %g for <%s> = %g [%g, %g]\n", coef, SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

               constant += coef * (SCIProwIsLocal(*row) ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var));
               SCIP_CALL( SCIPaddVarToRow(scip, *row, var, -coef) );
               continue;
            }

            break;
         }

         if( j < SCIProwGetNNonz(*row) )
         {
            SCIPdebugMsg(scip, "could not eliminate small coefficient\n");
            SCIP_CALL( SCIPreleaseRow(scip, row) );
            break;
         }

         if( violside == SCIP_SIDETYPE_LEFT )
         {
            SCIP_CALL( SCIPchgRowLhs(scip, *row, SCIProwGetLhs(*row) - constant) );
         }
         else
         {
            SCIP_CALL( SCIPchgRowRhs(scip, *row, SCIProwGetRhs(*row) - constant) );
         }

         /* update min/max coefficient */
         mincoef = SCIPgetRowMinCoef(scip, *row);
         maxcoef = SCIPgetRowMaxCoef(scip, *row);
      };

      /* avoid numerically very bad cuts */
      if( maxcoef / mincoef > cutmaxrange )
      {
         SCIPdebugMsg(scip, "drop row for constraint <%s> because range of coefficients is too large: mincoef = %g, maxcoef = %g -> range = %g\n",
            SCIPconsGetName(cons), mincoef, maxcoef, maxcoef / mincoef);
      }

      if( *row != NULL &&
         (  (violside == SCIP_SIDETYPE_LEFT  && SCIPisInfinity(scip, -SCIProwGetLhs(*row))) ||
            (violside == SCIP_SIDETYPE_RIGHT && SCIPisInfinity(scip,  SCIProwGetRhs(*row)))) )
      {
         SCIPdebugMsg(scip, "drop row for constraint <%s> because of very large side: %g\n", SCIPconsGetName(cons), violside == SCIP_SIDETYPE_LEFT ? -SCIProwGetLhs(*row) : SCIProwGetRhs(*row));
         SCIP_CALL( SCIPreleaseRow(scip, row) );
      }
   }

   return SCIP_OKAY;
}

/** returns whether one side of a constraint function is convex w.r.t. local bounds
 * i.e., if side == RIGHT, then returns whether constraint function is convex  w.r.t. local bounds
 * and   if side == LEFT,  then returns whether constraint function is concave w.r.t. local bounds
 */
static
SCIP_Bool isConvexLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         side                /**< constraint side to consider */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** xy;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->f != NULL);

   switch( consdata->convextype )
   {
   case SCIP_BIVAR_ALLCONVEX:
      /* always convex w.r.t. right hand side and concave w.r.t. left hand side */
      return side == SCIP_SIDETYPE_RIGHT;

   case SCIP_BIVAR_1CONVEX_INDEFINITE:
   {
      /* always not convex w.r.t. left hand side */
      if( side == SCIP_SIDETYPE_LEFT )
         return FALSE;

      xy = SCIPexprtreeGetVars(consdata->f);
      assert(xy != NULL);

      /* convex w.r.t. right hand side if one of the variables is fixed */
      return SCIPisEQ(scip, SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0])) ||
         SCIPisEQ(scip, SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1]));
   }

   case SCIP_BIVAR_CONVEX_CONCAVE:
   {
      xy = SCIPexprtreeGetVars(consdata->f);
      assert(xy != NULL);

      /* convex w.r.t. right hand side if y is fixed and
       * convex w.r.t. left  hand side if x is fixed */
      return (side == SCIP_SIDETYPE_RIGHT && SCIPisEQ(scip, SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1]))) ||
         (side == SCIP_SIDETYPE_LEFT && SCIPisEQ(scip, SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0])));
   }

   default:
      return FALSE;
   }  /*lint !e788*/
}

#ifdef SCIP_DEBUG
static
void printEstimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         side,               /**< violated side of constraint */
   SCIP_ROW*             row                 /**< row */
   )
{
   SCIP_CONSDATA* consdata;
   const char* varnames[2] = {"x", "y"};
   SCIP_VAR* x;
   SCIP_VAR* y;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   x = SCIPexprtreeGetVars(consdata->f)[0];
   y = SCIPexprtreeGetVars(consdata->f)[1];

   SCIPinfoMessage(scip, NULL, "splot [%g:%g] [%g:%g] ", SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y));
   SCIPexprtreePrint(consdata->f, SCIPgetMessagehdlr(scip), NULL, varnames, NULL);
   SCIPinfoMessage(scip, NULL, "%+g", side == SCIP_SIDETYPE_LEFT ? consdata->lhs : consdata->rhs);

   SCIPinfoMessage(scip, NULL, ", %g", SCIPisInfinity(scip, SCIProwGetRhs(row)) ? -SCIProwGetLhs(row) : -SCIProwGetRhs(row));
   for( i = 0; i < SCIProwGetNNonz(row); ++i )
   {
      SCIP_VAR* var;

      var = SCIPcolGetVar(SCIProwGetCols(row)[i]);
      if( var != x && var != y )
         continue;

      SCIPinfoMessage(scip, NULL, "%+g * %s", SCIProwGetVals(row)[i], var == x ? "x" : "y");
   }

   SCIPinfoMessage(scip, NULL, ", \"< echo '%g %g %g'\" with circles", SCIPgetSolVal(scip, sol, x), SCIPgetSolVal(scip, sol, y), consdata->activity);

   SCIPinfoMessage(scip, NULL, "\n");
}
#endif

/** tries to separate solution or LP solution by a linear cut
 *
 *  assumes that constraint violations have been computed
 */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Real             minefficacy,        /**< minimal efficacy of a cut if it should be added to the LP */
   SCIP_Bool             inenforcement,      /**< whether we are in constraint enforcement */
   SCIP_RESULT*          result,             /**< result of separation */
   SCIP_Real*            bestefficacy        /**< buffer to store best efficacy of a cut that was added to the LP, if found; or NULL if not of interest */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_SIDETYPE      violside;
   SCIP_Real          feasibility;
   SCIP_Real          efficacy;
   int                c;
   SCIP_ROW*          row;

   assert(scip         != NULL);
   assert(conshdlr     != NULL);
   assert(conss        != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(result       != NULL);

   *result = SCIP_FEASIBLE;

   if( bestefficacy != NULL )
      *bestefficacy = 0.0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         /* we are not feasible anymore */
         if( *result == SCIP_FEASIBLE )
            *result = SCIP_DIDNOTFIND;

         violside = SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT;

         /* generate cut */
         SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], sol, violside, conshdlrdata->cutmaxrange, &row) );
         if( row == NULL ) /* failed to generate cut */
            continue;

         if( sol == NULL )
            feasibility = SCIPgetRowLPFeasibility(scip, row);
         else
            feasibility = SCIPgetRowSolFeasibility(scip, row, sol);
         efficacy = -feasibility;

         SCIPdebug( printEstimator(scip, sol, conss[c], violside, row) );

         /* if cut is strong enough or it's weak but we separate on a convex function and accept weak cuts there, add cut to SCIP */
         if( (SCIPisGT(scip, efficacy, minefficacy) ||
              (inenforcement && SCIPisGT(scip, efficacy, SCIPfeastol(scip)) && isConvexLocal(scip, conss[c], violside))) &&
             SCIPisCutApplicable(scip, row) )
         {
            SCIP_Bool infeasible;

            /* cut cuts off solution sufficiently */
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, "cut for constraint <%s> is infeasible -> cutoff.\n", SCIPconsGetName(conss[c]));
               *result = SCIP_CUTOFF;
            }
            else
            {
               SCIPdebugMsg(scip, "added cut with efficacy %g for constraint <%s> violated by %g\n", efficacy, SCIPconsGetName(conss[c]), MAX(consdata->lhsviol, consdata->rhsviol));
               *result = SCIP_SEPARATED;
            }
            if( bestefficacy != NULL && efficacy > *bestefficacy )
               *bestefficacy = efficacy;

            /* mark row as not removable from LP for current node, if in enforcement */
            if( inenforcement && !conshdlrdata->enfocutsremovable )
               SCIPmarkRowNotRemovableLocal(scip, row);
         }
         else
         {
            SCIPdebugMsg(scip, "abandon cut since efficacy %g is too small or not applicable\n", efficacy);
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }

      if( *result == SCIP_CUTOFF )
         break;

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */
      if( c >= nusefulconss && *result == SCIP_FEASIBLE )
         break;
   }

   return SCIP_OKAY;
}

/** processes the event that a new primal solution has been found adds linearizations of all-convex constraints to the cutpool */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS**    conss;
   int            nconss;
   SCIP_CONSDATA* consdata;
   int            c;
   SCIP_SOL*      sol;
   SCIP_ROW*      row;
   SCIP_Real      x0y0[2];

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) != 0);

   conshdlr = (SCIP_CONSHDLR*)eventdata;

   nconss = SCIPconshdlrGetNConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   /* we are only interested in solution coming from some heuristic other than trysol, but not from the tree
    * the reason for ignoring trysol solutions is that they may come from an NLP solve in sepalp, where we already added linearizations,
    * or are from the tree, but postprocessed via proposeFeasibleSolution
    */
   if( SCIPsolGetHeur(sol) == NULL || SCIPsolGetHeur(sol) == conshdlrdata->trysolheur )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   SCIPdebugMsg(scip, "catched new sol event %" SCIP_EVENTTYPE_FORMAT " from heur <%s>; have %d conss\n", SCIPeventGetType(event), SCIPheurGetName(SCIPsolGetHeur(sol)), nconss);

   row = NULL;

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsLocal(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->convextype == SCIP_BIVAR_ALLCONVEX && !SCIPisInfinity(scip, consdata->rhs) )
      {
         SCIP_CALL( SCIPgetSolVals(scip, sol, 2, SCIPexprtreeGetVars(consdata->f), x0y0) );
         SCIP_CALL( generateLinearizationCut(scip, conshdlrdata->exprinterpreter, conss[c], x0y0, TRUE, &row) );
      }
      else
         continue;

      if( row == NULL )
         continue;

      assert(!SCIProwIsLocal(row));

      SCIP_CALL( SCIPaddPoolCut(scip, row) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** registers unfixed variables in nonlinear terms of violated constraints as external branching candidates
 * We score the variables by their gap between the convex envelope and the bivariate function in the current (x,y).
 * This value is given by the constraint violation, since we assume that cuts have been generated which support
 * the convex envelope in the LP.
 */
static
SCIP_RETCODE registerBranchingVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR**     xy;
   int            c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      SCIPdebugMsg(scip, "cons <%s> violation: %g %g\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);

      xy = SCIPexprtreeGetVars(consdata->f);
      assert(xy != NULL);

      /* @todo prefer binary before continuous, prefer unbounded before bounded */

      switch( consdata->convextype )
      {
      case SCIP_BIVAR_CONVEX_CONCAVE:
      {
         /* need to branch on the variable in which function is concave (or linear) */
         if( !SCIPisFeasZero(scip, consdata->lhsviol) )
         {
            /* regarding left hand side, we are concave in x and convex in y, so branch on x, if not fixed */
            if( !SCIPisEQ(scip, SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0])) )
            {
               SCIPdebugMsg(scip, "register variable x = <%s>[%g,%g] in convex-concave <%s> with violation %g %g\n", SCIPvarGetName(xy[0]), SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0]), SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
               SCIP_CALL( SCIPaddExternBranchCand(scip, xy[0], consdata->lhsviol, SCIP_INVALID) );
               ++*nnotify;
            }
         }
         if( !SCIPisFeasZero(scip, consdata->rhsviol) )
         {
            /* regarding right hand side, we are convex in x and concave in y, so branch on y, if not fixed */
            if( !SCIPisEQ(scip, SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1])) )
            {
               SCIPdebugMsg(scip, "register variable y = <%s>[%g,%g] in convex-concave <%s> with violation %g %g\n", SCIPvarGetName(xy[1]), SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1]), SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
               SCIP_CALL( SCIPaddExternBranchCand(scip, xy[1], consdata->lhsviol, SCIP_INVALID) );
               ++*nnotify;
            }
         }
         break;
      }

      case SCIP_BIVAR_1CONVEX_INDEFINITE:
      {
         if( !SCIPisFeasZero(scip, consdata->rhsviol) )
            if( SCIPisEQ(scip, SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0])) || SCIPisEQ(scip, SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1])) )
               break;

         /* register both variables, if not fixed */
         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0])) )
         {
            SCIPdebugMsg(scip, "register variable x = <%s>[%g,%g] in 1-convex <%s> with violation %g %g\n", SCIPvarGetName(xy[0]), SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0]), SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
            SCIP_CALL( SCIPaddExternBranchCand(scip, xy[0], consdata->lhsviol, SCIP_INVALID) );
            ++*nnotify;
         }

         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1])) )
         {
            SCIPdebugMsg(scip, "register variable y = <%s>[%g,%g] in 1-convex <%s> with violation %g %g\n", SCIPvarGetName(xy[1]), SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1]), SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
            SCIP_CALL( SCIPaddExternBranchCand(scip, xy[1], consdata->lhsviol, SCIP_INVALID) );
            ++*nnotify;
         }

         break;
      }

      case SCIP_BIVAR_ALLCONVEX:
      {
         if( SCIPisFeasZero(scip, consdata->lhsviol) )
            continue;
      }  /*lint -fallthrough*/

      default:
      {
         /* register both variables, if not fixed */
         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0])) )
         {
            SCIPdebugMsg(scip, "register variable x = <%s>[%g,%g] in allconvex <%s> with violation %g %g\n", SCIPvarGetName(xy[0]), SCIPvarGetLbLocal(xy[0]), SCIPvarGetUbLocal(xy[0]), SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
            SCIP_CALL( SCIPaddExternBranchCand(scip, xy[0], consdata->lhsviol, SCIP_INVALID) );
            ++*nnotify;
         }

         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1])) )
         {
            SCIPdebugMsg(scip, "register variable y = <%s>[%g,%g] in allconvex <%s> with violation %g %g\n", SCIPvarGetName(xy[1]), SCIPvarGetLbLocal(xy[1]), SCIPvarGetUbLocal(xy[1]), SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);
            SCIP_CALL( SCIPaddExternBranchCand(scip, xy[1], consdata->lhsviol, SCIP_INVALID) );
            ++*nnotify;
         }
      }
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** registers a nonlinear variable from a violated constraint as branching candidate that has a large absolute value in the relaxation */
static
SCIP_RETCODE registerLargeRelaxValueVariableForBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_VAR**            brvar               /**< buffer to store branching variable */
   )
{
   SCIP_CONSDATA*      consdata;
   SCIP_VAR*           var;
   SCIP_Real           val;
   SCIP_Real           brvarval;
   int                 i;
   int                 c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   *brvar = NULL;
   brvarval = -1.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->f != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( i = 0; i < 2; ++i )
      {
         var = SCIPexprtreeGetVars(consdata->f)[i];
         /* do not propose fixed variables */
         if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            continue;
         val = SCIPgetSolVal(scip, sol, var);
         if( REALABS(val) > brvarval )
         {
            brvarval = REALABS(val);
            *brvar = var;
         }
      }
   }

   if( *brvar != NULL )
   {
      SCIP_CALL( SCIPaddExternBranchCand(scip, *brvar, brvarval, SCIP_INVALID) );
   }

   return SCIP_OKAY;
}

/** enforces violated bivariate constraints where both nonlinear variables can be assumed to be fixed
 * apply a bound change to the remaining linear variable, or recognizing infeasibility
 */
static
SCIP_RETCODE enforceViolatedFixedNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            reduceddom,         /**< whether a domain has been reduced */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_CONSDATA*      consdata;
   SCIP_INTERVAL       nonlinact;
   SCIP_Real           lhs;
   SCIP_Real           rhs;
   int                 c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);
   assert(reduceddom != NULL);
   assert(infeasible != NULL);

   *reduceddom = FALSE;
   *infeasible = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      /* get activity for f(x,y) */
      SCIP_CALL( SCIPevalExprtreeLocalBounds(scip, consdata->f, SCIPinfinity(scip), &nonlinact) );
      assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), nonlinact));

      /* if all variables are fixed (at least up to epsilson), then the activity of the nonlinear part should be bounded */
      assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(nonlinact)));
      assert(!SCIPisInfinity(scip,  SCIPintervalGetSup(nonlinact)));

      if( !SCIPisInfinity(scip, -consdata->lhs) )
         lhs = consdata->lhs - SCIPintervalGetSup(nonlinact);
      else
         lhs = -SCIPinfinity(scip);

      if( !SCIPisInfinity(scip,  consdata->rhs) )
         rhs = consdata->rhs - SCIPintervalGetInf(nonlinact);
      else
         rhs = SCIPinfinity(scip);

      if( consdata->z != NULL )
      {
         SCIP_Bool tightened;
         SCIP_Real coef;

         coef = consdata->zcoef;
         assert(!SCIPisZero(scip, coef));

         SCIPdebugMsg(scip, "Linear constraint with one variable: %g <= %g <%s> <= %g\n", lhs, coef, SCIPvarGetName(consdata->z), rhs);

         /* possibly correct lhs/rhs */
         if( coef >= 0.0 )
         {
            if( !SCIPisInfinity(scip, -lhs) )
               lhs /= coef;
            if( !SCIPisInfinity(scip,  rhs) )
               rhs /= coef;
         }
         else
         {
            SCIP_Real h;
            h = rhs;
            if( !SCIPisInfinity(scip, -lhs) )
               rhs = lhs/coef;
            else
               rhs = SCIPinfinity(scip);

            if( !SCIPisInfinity(scip,  h) )
               lhs = h/coef;
            else
               lhs = -SCIPinfinity(scip);
         }
         SCIPdebugMsg(scip, "Linear constraint is a bound: %g <= <%s> <= %g\n", lhs, SCIPvarGetName(consdata->z), rhs);

         if( !SCIPisInfinity(scip, -lhs) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, consdata->z, lhs, TRUE, infeasible, &tightened) );
            if( *infeasible )
            {
               SCIPdebugMsg(scip, "Lower bound leads to infeasibility.\n");
               return SCIP_OKAY;
            }
            if( tightened )
            {
               SCIPdebugMsg(scip, "Lower bound changed.\n");
               *reduceddom = TRUE;
               return SCIP_OKAY;
            }
         }

         if( !SCIPisInfinity(scip, rhs) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, consdata->z, rhs, TRUE, infeasible, &tightened) );
            if( *infeasible )
            {
               SCIPdebugMsg(scip, "Upper bound leads to infeasibility.\n");
               return SCIP_OKAY;
            }
            if( tightened )
            {
               SCIPdebugMsg(scip, "Upper bound changed.\n");
               *reduceddom = TRUE;
               return SCIP_OKAY;
            }
         }
      }
      else
      {
         /* no variable, thus check feasibility of lhs <= 0.0 <= rhs */
         *infeasible = SCIPisFeasGT(scip, lhs, 0.0) || SCIPisFeasLT(scip, rhs, 0.0);
      }
   }

   return SCIP_OKAY;
}

/** tightens bounds on a variable to given interval */
static
SCIP_RETCODE propagateBoundsTightenVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which bounds to tighten */
   SCIP_INTERVAL         bounds,             /**< new bounds */
   SCIP_CONS*            cons,               /**< constraint that is propagated */
   SCIP_RESULT*          result,             /**< pointer where to update the result of the propagation call */
   int*                  nchgbds             /**< buffer where to add the the number of changed bounds */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;
   SCIP_Real bnd;

   assert(scip    != NULL);
   assert(var     != NULL);
   assert(result  != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   if( SCIPintervalIsPositiveInfinity(SCIPinfinity(scip), bounds) ||
      SCIPintervalIsNegativeInfinity(SCIPinfinity(scip), bounds) ||
      SCIPintervalIsEmpty(SCIPinfinity(scip), bounds) )
   {
      /* domain outside [-infty, +infty] or empty -> declare node infeasible */
      SCIPdebugMsg(scip, "found <%s> infeasible due to domain propagation for variable <%s>\n", cons != NULL ? SCIPconsGetName(cons) : "???", SCIPvarGetName(var));  /*lint !e585*/
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -SCIPintervalGetInf(bounds)) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      bnd = SCIPadjustedVarLb(scip, var, SCIPintervalGetInf(bounds));
      SCIP_CALL( SCIPtightenVarLb(scip, var, bnd, FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMsg(scip, "found <%s> infeasible due to domain propagation for variable <%s>\n", cons != NULL ? SCIPconsGetName(cons) : "???", SCIPvarGetName(var));  /*lint !e585*/
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( tightened )
      {
         SCIPdebugMsg(scip, "tightened lower bound of variable <%s> in constraint <%s> to %g\n", SCIPvarGetName(var), cons != NULL ? SCIPconsGetName(cons) : "???", SCIPvarGetLbLocal(var));  /*lint !e585*/
         ++*nchgbds;
         *result = SCIP_REDUCEDDOM;
      }
   }

   if( !SCIPisInfinity(scip, SCIPintervalGetSup(bounds)) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      bnd = SCIPadjustedVarLb(scip, var, SCIPintervalGetSup(bounds));
      SCIP_CALL( SCIPtightenVarUb(scip, var, bnd, FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMsg(scip, "found <%s> infeasible due to domain propagation for variable <%s>\n", cons != NULL ? SCIPconsGetName(cons) : "???", SCIPvarGetName(var));  /*lint !e585*/
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( tightened )
      {
         SCIPdebugMsg(scip, "tightened upper bound of variable <%s> in constraint <%s> to %g\n", SCIPvarGetName(var), cons != NULL ? SCIPconsGetName(cons) : "???", SCIPvarGetUbLocal(var));  /*lint !e585*/
         ++*nchgbds;
         *result = SCIP_REDUCEDDOM;
      }
   }

   return SCIP_OKAY;
}

/** tightens bounds of z in a single bivariate constraint
 * checks for redundancy and infeasibility
 */
static
SCIP_RETCODE propagateBoundsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation call */
   int*                  nchgbds,            /**< buffer where to add the the number of changed bounds */
   SCIP_Bool*            redundant           /**< buffer where to store whether constraint has been found to be redundant */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL  consbounds;    /* left and right side of constraint */
   SCIP_INTERVAL  ftermactivity; /* activity of f(x,y) */
   SCIP_INTERVAL  ztermactivity; /* activity of c*z */
   SCIP_INTERVAL  consactivity;  /* activity of f(x,y) + c*z */
   SCIP_INTERVAL  tmp;
   SCIP_Bool      cutoff;

   assert(scip    != NULL);
   assert(cons    != NULL);
   assert(result  != NULL);
   assert(nchgbds != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->exprgraphnode != NULL);

   *result = SCIP_DIDNOTRUN;
   *redundant = FALSE;

   /* extend interval by epsilon to avoid cutoff in forward propagation if constraint is only almost feasible */
   SCIPintervalSetBounds(&consbounds,
      -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -consdata->lhs+SCIPepsilon(scip)),    /*lint !e666*/
      +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  consdata->rhs+SCIPepsilon(scip)) );  /*lint !e666*/

   /* get activity for f(x,y) */
   ftermactivity = SCIPexprgraphGetNodeBounds(consdata->exprgraphnode);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), ftermactivity) );

   /* get activity for c*z */
   if( consdata->z != NULL )
   {
      SCIPintervalSetBounds(&ztermactivity,
         -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -MIN(SCIPvarGetLbLocal(consdata->z), SCIPvarGetUbLocal(consdata->z))),   /*lint !e666*/
         +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  MAX(SCIPvarGetLbLocal(consdata->z), SCIPvarGetUbLocal(consdata->z))));  /*lint !e666*/
      SCIPintervalMulScalar(INTERVALINFTY, &ztermactivity, ztermactivity, consdata->zcoef);
   }
   else
   {
      SCIPintervalSet(&ztermactivity, 0.0);
   }

   /* get activity for f(x,y)+c*z */
   SCIPintervalAdd(INTERVALINFTY, &consactivity, ftermactivity, ztermactivity);

   /* check redundancy */
   if( SCIPintervalIsSubsetEQ(INTERVALINFTY, consactivity, consbounds) )
   {
      SCIPdebugMsg(scip, "found constraint <%s> to be redundant: sides: [%g, %g], activity: [%g, %g]\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity));
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   /* check infeasibility */
   if( SCIPintervalAreDisjoint(consbounds, consactivity) )
   {
      SCIPdebugMsg(scip, "found constraint <%s> to be infeasible; sides: [%g, %g], activity: [%g, %g], infeas: %g\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity),
         MAX(consdata->lhs - SCIPintervalGetSup(consactivity), SCIPintervalGetInf(consactivity) - consdata->rhs));  /*lint !e666*/
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* try to tighten bounds on z */
   if( consdata->z != NULL )
   {
      *result = SCIP_DIDNOTFIND;

      /* compute ([lhs, rhs] - f([xlb,xub], [ylb,yub])) / zcoef */
      SCIPintervalSub(INTERVALINFTY, &tmp, consbounds, ftermactivity);
      SCIPintervalDivScalar(INTERVALINFTY, &tmp, tmp, consdata->zcoef);

      SCIP_CALL( propagateBoundsTightenVar(scip, consdata->z, tmp, cons, result, nchgbds) );

      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;

      if( *result == SCIP_SUCCESS )
      {
         SCIPintervalSetBounds(&ztermactivity,
            -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -MIN(SCIPvarGetLbLocal(consdata->z), SCIPvarGetUbLocal(consdata->z))),   /*lint !e666*/
            +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  MAX(SCIPvarGetLbLocal(consdata->z), SCIPvarGetUbLocal(consdata->z))));  /*lint !e666*/
         SCIPintervalMulScalar(INTERVALINFTY, &ztermactivity, ztermactivity, consdata->zcoef);
      }
   }

   /* set bounds for exprgraphnode = [lhs,rhs] - c*z */
   SCIPintervalSub(INTERVALINFTY, &tmp, consbounds, ztermactivity);
   SCIPexprgraphTightenNodeBounds(conshdlrdata->exprgraph, consdata->exprgraphnode, tmp, 0.05, INTERVALINFTY, &cutoff);
   if( cutoff )
   {
      SCIPdebugMsg(scip, "found constraint <%s> infeasible%s\n", SCIPconsGetName(cons), SCIPinProbing(scip) ? " in probing" : "");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** calls domain propagation for a set of constraints */
static
SCIP_RETCODE propagateBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation calls */
   int*                  nchgbds,            /**< buffer where to add the the number of changed bounds */
   int*                  ndelconss           /**< buffer where to increase if a constraint was deleted (locally) due to redundancy */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT propresult;
   SCIP_Bool   redundant;
   SCIP_Bool   domainerror;
   int         roundnr;
   SCIP_Bool   success;
   int         nvars;
   SCIP_VAR**  vars;
   SCIP_EXPRGRAPHNODE** varnodes;
   SCIP_Bool   cutoff;
   int         c;
   int         i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   *result = SCIP_DIDNOTRUN;

   if( nconss == 0 )
      return SCIP_OKAY;

   if( conshdlrdata->ispropagated )
   {
      /* check whether there was also no tightening in the bounds of the linear variables
       * @todo put this in processLinearVarEvent
       */
      for( c = 0; c < nconss; ++c )
      {
         assert(conss[c] != NULL);  /*lint !e613*/

         if( SCIPconsIsMarkedPropagate(conss[c]) )  /*lint !e613*/
            break;
      }
      if( c == nconss )
         return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   roundnr = 0;
   do
   {
      success = FALSE;

      SCIPdebugMsg(scip, "starting domain propagation round %d for %d constraints\n", roundnr, nconss);

      conshdlrdata->ispropagated = TRUE;

      /* propagate variable bounds through expression graph
       * roundnr == 0 clears remainings from a previous backward propagation
       * @todo could give FALSE if no linear variable in the constraints had been relaxed since last time
       */
      SCIP_CALL( SCIPexprgraphPropagateVarBounds(conshdlrdata->exprgraph, INTERVALINFTY, roundnr == 0, &domainerror) );

      if( domainerror )
      {
         SCIPdebugMsg(scip, "current bounds out of domain for some expression, do cutoff\n");
         *result = SCIP_CUTOFF;
         break;
      }

      /* check for redundancy and infeasibility of constraints
       * tighten bounds on linear variables
       * setup bounds for expression graph nodes */
      for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
      {
         assert(conss != NULL);
         if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) )
            continue;

         SCIP_CALL( propagateBoundsCons(scip, conshdlr, conss[c], &propresult, nchgbds, &redundant) );
         if( propresult != SCIP_DIDNOTFIND && propresult != SCIP_DIDNOTRUN )
         {
            *result = propresult;
            success = TRUE;
         }
         if( redundant )
         {
            SCIPdebugMsg(scip, "delete redundant constraint <%s> locally\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
            ++*ndelconss;
         }

         SCIP_CALL( SCIPunmarkConsPropagate(scip, conss[c]) );
      }
      if( *result == SCIP_CUTOFF )
         break;

      /* propagate backward through expression graph */
      SCIPdebugMsg(scip, "start backward propagation in expression graph\n");

      /* compute bound tightenings for nonlinear variables */
      SCIPexprgraphPropagateNodeBounds(conshdlrdata->exprgraph, INTERVALINFTY, 0.05, &cutoff);

      if( cutoff )
      {
         SCIPdebugMsg(scip, "backward propagation found problem infeasible%s\n", SCIPinProbing(scip) ? " in probing" : "");
         *result = SCIP_CUTOFF;
         break;
      }

      /* put back new bounds into SCIP variables */
      nvars = SCIPexprgraphGetNVars(conshdlrdata->exprgraph);
      vars  = (SCIP_VAR**)SCIPexprgraphGetVars(conshdlrdata->exprgraph);
      varnodes = SCIPexprgraphGetVarNodes(conshdlrdata->exprgraph);
      propresult = SCIP_DIDNOTFIND;
      for( i = 0; i < nvars && propresult != SCIP_CUTOFF; ++i )
      {
         SCIP_CALL( propagateBoundsTightenVar(scip, vars[i], SCIPexprgraphGetNodeBounds(varnodes[i]), NULL, &propresult, nchgbds) );
      }
      if( propresult != SCIP_DIDNOTFIND )
      {
         *result = propresult;
         success = TRUE;
      }
   }
   while( success && *result != SCIP_CUTOFF && ++roundnr < conshdlrdata->maxproprounds );

   return SCIP_OKAY;
}


/** Given a solution where every bivariate constraint is either feasible or can be made feasible by
 * moving the linear variable, construct the corresponding feasible solution and pass it to the trysol heuristic.
 * The method assumes that this is always possible and that not all constraints are feasible already.
 */
static
SCIP_RETCODE proposeFeasibleSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to process */
   SCIP_Bool*            success             /**< buffer to store whether we succeeded to construct a solution that satisfies all provided constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_SOL* newsol;
   SCIP_VAR* var;
   int c;
   SCIP_Real viol;
   SCIP_Real delta;
   SCIP_Real gap;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(success != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->trysolheur != NULL);

   *success = FALSE;

   /* don't propose new solutions if not in presolve or solving */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
      return SCIP_OKAY;

   if( sol != NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &newsol, sol) );
   }
   else
   {
      SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
   }
   SCIP_CALL( SCIPunlinkSol(scip, newsol) );

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* recompute violation of solution in case solution has changed
       * get absolution violation and sign */
      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conshdlr, conss[c], newsol) );  /*lint !e613*/
         viol = consdata->lhs - consdata->activity;
      }
      else if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conshdlr, conss[c], newsol) );  /*lint !e613*/
         viol = consdata->rhs - consdata->activity;
      }
      else
         continue; /* constraint is satisfied */

      assert(viol != 0.0);
      if( consdata->mayincreasez &&
         ((viol > 0.0 && consdata->zcoef > 0.0) || (viol < 0.0 && consdata->zcoef < 0.0)) )
      {
         /* have variable where increasing makes the constraint less violated */
         var = consdata->z;
         /* compute how much we would like to increase var */
         delta = viol / consdata->zcoef;
         assert(delta > 0.0);
         /* if var has an upper bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            gap = SCIPvarGetUbGlobal(var) - SCIPgetSolVal(scip, newsol, var);
            delta = MIN(MAX(0.0, gap), delta);
         }
         if( SCIPisPositive(scip, delta) )
         {
            /* if variable is integral, round delta up so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPceil(scip, delta);

            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            SCIPdebugMsg(scip, "increase <%s> by %g to %g\n", SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->zcoef * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      assert(viol != 0.0);
      if( consdata->maydecreasez &&
         ((viol > 0.0 && consdata->zcoef < 0.0) || (viol < 0.0 && consdata->zcoef > 0.0)) )
      {
         /* have variable where decreasing makes constraint less violated */
         var = consdata->z;
         /* compute how much we would like to decrease var */
         delta = viol / consdata->zcoef;
         assert(delta < 0.0);
         /* if var has a lower bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
         {
            gap = SCIPgetSolVal(scip, newsol, var) - SCIPvarGetLbGlobal(var);
            delta = MAX(MIN(0.0, gap), delta);
         }
         if( SCIPisNegative(scip, delta) )
         {
            /* if variable is integral, round delta down so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPfloor(scip, delta);
            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            SCIPdebugMsg(scip, "increase <%s> by %g to %g\n", SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->zcoef * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      /* still here... so maybe we could not make constraint feasible due to variable bounds
       * check if we are feasible w.r.t. (relative) feasibility tolerance */
      SCIP_CALL( computeViolation(scip, conshdlr, conss[c], newsol) );  /*lint !e613*/
      /* if still violated, we give up */
      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         break;

      /* if objective value is not better than current upper bound, we give up */
      if( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) && !SCIPisSumLT(scip, SCIPgetSolTransObj(scip, newsol), SCIPgetUpperbound(scip)) )
         break;
   }

   /* if we have a solution that should satisfy all nonlinear constraints and has a better objective than the current upper bound,
    * then pass it to the trysol heuristic */
   if( c == nconss )
   {
      SCIPdebugMsg(scip, "pass solution with objective value %g to trysol heuristic\n", SCIPgetSolTransObj(scip, newsol));

      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, newsol) );
      *success = TRUE;
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/** creates bivariate constraint from quadratic constraint data of the form
 * lhs <= xsqrcoef * x^2 + xlincoef * x + ysqrcoef * y^2 + ylincoef * y + bilincoef * x*y + zcoef * z <= rhs
 */
static
SCIP_RETCODE createConsFromQuadTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            srccons,            /**< source constraint to take attributes from */
   SCIP_CONS**           cons,               /**< pointer to store new constraint */
   const char*           name,               /**< name of new constraint */
   SCIP_VAR*             x,                  /**< first nonlinear variable */
   SCIP_VAR*             y,                  /**< second nonlinear variable */
   SCIP_VAR*             z,                  /**< linear variable, can be NULL */
   SCIP_Real             coefxx,             /**< coefficient of x^2 */
   SCIP_Real             coefx,              /**< coefficient of x */
   SCIP_Real             coefyy,             /**< coefficient of y^2 */
   SCIP_Real             coefy,              /**< coefficient of y */
   SCIP_Real             coefxy,             /**< coefficient of x*y */
   SCIP_Real             coefz,              /**< coefficient of z */
   SCIP_Real             lhs,                /**< left-hand-side */
   SCIP_Real             rhs                 /**< right-hand-side */
   )
{
   SCIP_Real mult;
   SCIP_VAR* xy[2];
   SCIP_BIVAR_CONVEXITY convextype;
   SCIP_EXPR* e;
   SCIP_EXPRTREE* exprtree;

   SCIP_EXPR* children[2];
   SCIP_Real  lincoefs[2];
   SCIP_QUADELEM quadelems[3];
   int nquadelems;

   assert(scip != NULL);
   assert(srccons != NULL);
   assert(cons != NULL);
   assert(name != NULL);

   assert(x != NULL);
   assert(y != NULL);
   assert(SCIPisLE(scip, lhs, rhs));

   if( coefxx >= 0 && coefyy >= 0 && 4 * coefxx * coefyy >= coefxy * coefxy )
   {
      /* quadratic term is convex in both variables (jointly) */
      mult =  1.0;
      convextype = SCIP_BIVAR_ALLCONVEX;
   }
   else if( coefxx <= 0 && coefyy <= 0 && 4 * coefxx * coefyy >= coefxy * coefxy )
   {
      /* quadratic term is concave in both variables (jointly) */
      mult = -1.0;
      convextype = SCIP_BIVAR_ALLCONVEX;
   }
   else if( coefxx > 0 && coefyy > 0 )
   {
      /* indefinite but 1-convex */
      assert(4 * coefxx * coefyy < coefxy * coefxy); /* assert indefiniteness */
      mult =  1.0;
      convextype = SCIP_BIVAR_1CONVEX_INDEFINITE;
   }
   else if( coefxx < 0 && coefyy < 0 )
   {
      /* indefinite but 1-convex */
      assert(4 * coefxx * coefyy < coefxy * coefxy); /* assert indefiniteness */
      mult = -1.0;
      convextype = SCIP_BIVAR_1CONVEX_INDEFINITE;
   }
   else
   {
      /* convex in one variable and concave in other variable */
      assert(coefxx * coefyy <= 0);
      convextype = SCIP_BIVAR_CONVEX_CONCAVE;
      if( coefxx != 0.0 )
      {
         /* if coefxx < 0 (and thus coefyy >= 0) f(x,y) is concave in x and convex in y
          * but we need convex in x and concave in y, thus we multiply by -1
          */
         if( coefxx < 0.0 )
            mult =  -1.0;
         else
            mult =  1.0;
      }
      else if( coefyy != 0.0 )
      {
         /* coefxx == 0.0 */
         /* if coefyy < 0 (and coefxx == 0) f(x,y) is concave in y and convex in x
          * otherwise we convert to convex in y and concave in x by multiplying by -1
          */
         if( coefyy < 0.0 )
            mult =  1.0;
         else
            mult =  -1.0;
      }
      else
      {
         /* coefxx == 0.0 && coefyy == 0.0 && coefxy != 0.0 */
         assert(coefxy != 0.0);
         mult = 1.0;
      }
   }

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[0], SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[1], SCIP_EXPR_VARIDX, 1) );

   lincoefs[0] = coefx * mult;
   lincoefs[1] = coefy * mult;

   nquadelems = 0;
   if( coefxx != 0.0 )
   {
      quadelems[nquadelems].idx1 = 0;
      quadelems[nquadelems].idx2 = 0;
      quadelems[nquadelems].coef = coefxx * mult;
      ++nquadelems;
   }
   if( coefyy != 0.0 )
   {
      quadelems[nquadelems].idx1 = 1;
      quadelems[nquadelems].idx2 = 1;
      quadelems[nquadelems].coef = coefyy * mult;
      ++nquadelems;
   }
   if( coefxy != 0.0 )
   {
      quadelems[nquadelems].idx1 = 0;
      quadelems[nquadelems].idx2 = 1;
      quadelems[nquadelems].coef = coefxy * mult;
      ++nquadelems;
   }

   SCIP_CALL( SCIPexprCreateQuadratic(SCIPblkmem(scip), &e, 2, children, 0.0, (coefx != 0.0 || coefy != 0.0) ? lincoefs : NULL, nquadelems, quadelems) );  /*lint !e826*/
   assert(e != NULL);

   xy[0] = x;
   xy[1] = y;

   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, e, 2, 0, NULL) );
   SCIP_CALL( SCIPexprtreeSetVars(exprtree, 2, xy) );

   if( mult == -1.0 )
   {
      SCIP_Real tmp;
      tmp = lhs;
      lhs = -rhs;
      rhs = -tmp;
      coefz = -coefz;
   }
   else
   {
      assert(mult == 1.0);
   }

   SCIPdebugMsg(scip, "upgrading constraint <%s> to bivariate constraint <%s> with convexity type %d\n", SCIPconsGetName(srccons), name, convextype);

   SCIP_CALL( SCIPcreateConsBivariate(scip, cons, name,
         exprtree, convextype, z, coefz, lhs, rhs,
         SCIPconsIsInitial(srccons), SCIPconsIsSeparated(srccons), SCIPconsIsEnforced(srccons),
         SCIPconsIsChecked(srccons), SCIPconsIsPropagated(srccons), SCIPconsIsLocal(srccons),
         SCIPconsIsModifiable(srccons), SCIPconsIsDynamic(srccons), SCIPconsIsRemovable(srccons),
         SCIPconsIsStickingAtNode(srccons)) );
   SCIPdebugPrintCons(scip, *cons, NULL);

   SCIP_CALL( SCIPexprtreeFree(&exprtree) );

   return SCIP_OKAY;
}

/** creates expression tree for monomial of the form coef * x^p * y^q with x >= 0 and y >= 0 and checks its convexity type */
static
SCIP_RETCODE createExprtreeFromMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             coef,               /**< monomial coefficient */
   SCIP_Real             p,                  /**< exponent of x */
   SCIP_Real             q,                  /**< exponent of y */
   SCIP_EXPRTREE**       exprtree,           /**< buffer to store pointer to expression tree */
   SCIP_Real*            mult,               /**< buffer to store multiplicator for generated expression tree */
   SCIP_BIVAR_CONVEXITY* convextype          /**< buffer to store convexity type of expression tree */
   )
{
   SCIP_Bool swapvars;
   SCIP_EXPR* children[2];
   int childidxs[2];
   SCIP_Real exponents[2];
   SCIP_VAR* vars[2];
   SCIP_EXPR* e;
   SCIP_EXPRDATA_MONOMIAL* monomial;

   assert(scip != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(!SCIPisZero(scip, coef));
   assert(!SCIPisZero(scip, p));
   assert(!SCIPisZero(scip, q));
   assert(exprtree != NULL);
   assert(mult != NULL);
   assert(convextype != NULL);

   /* determine convexity type, and whether to negate monomial or swap variables */
   *mult = coef < 0.0 ? -1.0 : 1.0;  /* for the check, assume that monomial has positive coefficient */
   swapvars = FALSE;
   *convextype = SCIP_BIVAR_UNKNOWN;
   if( (p + q >= 1.0 && ((p > 1.0 && q < 0.0) || (p < 0.0 && q > 1.0))) ||
      (p < 0.0 && q < 0.0) )
   {
      *convextype = SCIP_BIVAR_ALLCONVEX;
   }
   else if( (p > 1.0 && q > 1.0) || (p + q < 1.0 && ((p > 1.0 && q < 0.0) || (p < 0.0 && q > 1.0))) )
   {
      *convextype = SCIP_BIVAR_1CONVEX_INDEFINITE;
   }
   else if( (p < 0.0 || p > 1.0) && q > 0.0 && q < 1.0 )
   {
      *convextype = SCIP_BIVAR_CONVEX_CONCAVE;
   }
   else if( (p < 0.0 || p > 1.0) && q == 1.0 )
   {
      *mult *= -1.0;
      swapvars = TRUE;
      *convextype = SCIP_BIVAR_CONVEX_CONCAVE;
   }
   else if( (q < 0.0 || q > 1.0) && p > 0.0 && p <= 1.0 )
   {
      swapvars = TRUE;
      *convextype = SCIP_BIVAR_CONVEX_CONCAVE;
   }
   else if( p > 0.0 && p < 1.0 && q > 0.0 && q < 1.0 && p + q > 1.0 )
   {
      *mult *= -1.0;
      *convextype = SCIP_BIVAR_1CONVEX_INDEFINITE;
   }
   else if( p == 1.0 && q > 0.0 && q < 1.0 )
   {
      *convextype = SCIP_BIVAR_CONVEX_CONCAVE;
   }
   else if( q == 1.0 && p > 0.0 && p < 1.0 )
   {
      swapvars = TRUE;
      *convextype = SCIP_BIVAR_CONVEX_CONCAVE;
   }
   else if( p == 1.0 && q == 1.0 )
   {
      *convextype = SCIP_BIVAR_CONVEX_CONCAVE;
   }
   else if( p > 0.0 && p < 1.0 && q > 0.0 && q < 1.0 && p + q <= 1.0 )
   {
      *mult *= -1.0;
      *convextype = SCIP_BIVAR_ALLCONVEX;
   }
   assert(*convextype != SCIP_BIVAR_UNKNOWN); /* there should be no case where this can still happen */

   /* setup expression tree */
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[0], SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[1], SCIP_EXPR_VARIDX, 1) );
   childidxs[0] = 0;
   childidxs[1] = 1;
   if( !swapvars )
   {
      exponents[0] = p;
      exponents[1] = q;
      vars[0] = x;
      vars[1] = y;
   }
   else
   {
      exponents[0] = q;
      exponents[1] = p;
      vars[0] = y;
      vars[1] = x;
   }
   SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomial, *mult*coef, 2, childidxs, exponents) );

   SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), &e, 2, children, 1, &monomial, 0.0, FALSE) );
   assert( e != NULL );

   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), exprtree, e, 2, 0, NULL) );
   SCIP_CALL( SCIPexprtreeSetVars(*exprtree, 2, vars) );

   return SCIP_OKAY;
}

/** creates bivariate constraint from monomial of the form coef * x^p * y^q with x >= 0 and y >= 0
 * lhs <= coef * x^p * y^q + zcoef * z <= rhs
 */
static
SCIP_RETCODE createConsFromMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            srccons,            /**< source constraint to take attributes from, or NULL */
   SCIP_CONS**           cons,               /**< pointer to store new constraint */
   const char*           name,               /**< name of new constraint */
   SCIP_VAR*             x,                  /**< first nonlinear variable */
   SCIP_VAR*             y,                  /**< second nonlinear variable */
   SCIP_VAR*             z,                  /**< linear variable, can be NULL */
   SCIP_Real             coef,               /**< monomial coefficient */
   SCIP_Real             p,                  /**< exponent of x */
   SCIP_Real             q,                  /**< exponent of y */
   SCIP_Real             zcoef,              /**< coefficient of z */
   SCIP_Real             lhs,                /**< left-hand-side */
   SCIP_Real             rhs                 /**< right-hand-side */
   )
{
   SCIP_Real mult;
   SCIP_BIVAR_CONVEXITY convextype;
   SCIP_EXPRTREE* exprtree;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(name != NULL);

   assert(x != NULL);
   assert(y != NULL);
   assert(!SCIPisZero(scip, coef));
   assert(!SCIPisZero(scip, p));
   assert(!SCIPisZero(scip, q));
   assert(SCIPisLE(scip, lhs, rhs));

   SCIP_CALL( createExprtreeFromMonomial(scip, x, y, coef, p, q, &exprtree, &mult, &convextype) );

   if( mult == -1.0 )
   {
      SCIP_Real tmp;
      tmp = lhs;
      lhs = -rhs;
      rhs = -tmp;
      zcoef = -zcoef;
   }
   else
   {
      assert(mult == 1.0);
   }

   SCIPdebugMsg(scip, "upgrading monomial %g<%s>^%g<%s>^%g from constraint <%s> to bivariate constraint with convexity type %d\n",  /*lint !e585*/
      coef, SCIPvarGetName(x), p, SCIPvarGetName(y), q, srccons != NULL ? SCIPconsGetName(srccons) : "n/a", convextype);  /*lint !e585*/

   if( srccons != NULL )
   {
      SCIP_CALL( SCIPcreateConsBivariate(scip, cons, name,
            exprtree, convextype, z, zcoef, lhs, rhs,
            SCIPconsIsInitial(srccons), SCIPconsIsSeparated(srccons), SCIPconsIsEnforced(srccons),
            SCIPconsIsChecked(srccons), SCIPconsIsPropagated(srccons), SCIPconsIsLocal(srccons),
            SCIPconsIsModifiable(srccons), SCIPconsIsDynamic(srccons), SCIPconsIsRemovable(srccons),
            SCIPconsIsStickingAtNode(srccons)) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsBivariate(scip, cons, name,
            exprtree, convextype, z, zcoef, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   }
   SCIPdebugPrintCons(scip, *cons, NULL);

   SCIP_CALL( SCIPexprtreeFree(&exprtree) );

   return SCIP_OKAY;
}

/** helper function to enforce constraints */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_CONS*         maxviolcons;
   SCIP_Real          maxviol;
   SCIP_RESULT        propresult;
   SCIP_RESULT        separateresult;
   int                dummy;
   int                nnotify;
   SCIP_Real          sepaefficacy;
   SCIP_Real          minefficacy;
   SCIP_Real          leastpossibleefficacy;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, sol, &maxviolcons) );
   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   /* if we are above the 100'th enforcement round for this node, something is strange (maybe the relaxation does not
    * think that the cuts we add are violated, or we do ECP on a high-dimensional convex function) in this case, check
    * if some limit is hit or SCIP should stop for some other reason and terminate enforcement by creating a dummy node
    * (in optimized more, returning SCIP_INFEASIBLE in *result would be sufficient, but in debug mode this would give an
    * assert in scip.c) the reason to wait for 100 rounds is to avoid calls to SCIPisStopped in normal runs, which may
    * be expensive we only increment nenforounds until 101 to avoid an overflow
    */
   if( conshdlrdata->lastenfonode == SCIPgetCurrentNode(scip) )
   {
      if( conshdlrdata->nenforounds > 100 )
      {
         if( SCIPisStopped(scip) )
         {
            SCIP_NODE* child;

            SCIP_CALL( SCIPcreateChild(scip, &child, 1.0, SCIPnodeGetEstimate(SCIPgetCurrentNode(scip))) );
            *result = SCIP_BRANCHED;

            return SCIP_OKAY;
         }
      }
      else
         ++conshdlrdata->nenforounds;
   }
   else
   {
      conshdlrdata->lastenfonode = SCIPgetCurrentNode(scip);
      conshdlrdata->nenforounds = 0;
   }

   consdata = SCIPconsGetData(maxviolcons);
   assert(consdata != NULL);
   maxviol = consdata->lhsviol + consdata->rhsviol;
   assert(SCIPisGT(scip, maxviol, SCIPfeastol(scip)));

   SCIPdebugMsg(scip, "enforcement with max violation %g in cons <%s> for %s solution\n", maxviol, SCIPconsGetName(maxviolcons),
         sol == NULL ? "LP" : "relaxation");

   /* run domain propagation */
   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, &dummy, &dummy) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we would like a cut that is efficient enough that it is not redundant in the LP (>lpfeastol)
    * however, if the maximal violation is very small, also the best cut efficacy cannot be large
    * thus, in the latter case, we are also happy if the efficacy is at least, say, 75% of the maximal violation
    * but in any case we need an efficacy that is at least lpfeastol
    */
   minefficacy = MIN(0.75*maxviol, 2.0 * SCIPlpfeastol(scip));  /*lint !e666*/
   minefficacy = MAX(minefficacy, SCIPlpfeastol(scip));  /*lint !e666*/
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, minefficacy, TRUE, &separateresult,
         &sepaefficacy) );
   if( separateresult == SCIP_SEPARATED || separateresult == SCIP_CUTOFF )
   {
      SCIPdebugMessage("separation succeeded (bestefficacy = %g, minefficacy = %g, cutoff = %d)\n", sepaefficacy,
         minefficacy, separateresult == SCIP_CUTOFF);
      *result = separateresult;
      return SCIP_OKAY;
   }

   /* we are not feasible, the whole node is not infeasible, and we cannot find a good cut
    * -> collect variables for branching
    */

   SCIPdebugMsg(scip, "separation failed (bestefficacy = %g < %g = minefficacy ); max viol: %g\n", sepaefficacy,
      minefficacy, maxviol);

   /* find branching candidates */
   SCIP_CALL( registerBranchingVariables(scip, conss, nconss, &nnotify) );

   leastpossibleefficacy = SCIPlpfeastol(scip);
   if( nnotify == 0 && !solinfeasible && minefficacy > leastpossibleefficacy )
   {
      /* fallback 1: we also have no branching candidates, so try to find a weak cut */
      SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, leastpossibleefficacy, TRUE,
            &separateresult, &sepaefficacy) );
      if( separateresult == SCIP_SEPARATED || separateresult == SCIP_CUTOFF )
      {
         *result = separateresult;
         return SCIP_OKAY;
      }
   }

   if( nnotify == 0 && !solinfeasible )
   {
      /* fallback 2: separation probably failed because of numerical difficulties with a convex constraint;
       * if noone declared solution infeasible yet and we had not even found a weak cut, try to resolve by branching
       */
      SCIP_VAR* brvar = NULL;
      SCIP_CALL( registerLargeRelaxValueVariableForBranching(scip, conss, nconss, sol, &brvar) );
      if( brvar == NULL )
      {
         /* fallback 3: all nonlinear variables in all violated constraints seem to be fixed -> treat as linear
          * constraint in one variable
          */
         SCIP_Bool reduceddom;
         SCIP_Bool infeasible;

         SCIP_CALL( enforceViolatedFixedNonlinear(scip, conss, nconss, &reduceddom, &infeasible) );
         /* if the linear constraints are actually feasible, then adding them and returning SCIP_CONSADDED confuses SCIP
          * when it enforces the new constraints again and nothing resolves the infeasiblity that we declare here thus,
          * we only add them if considered violated, and otherwise claim the solution is feasible (but print a warning)
          */
         if ( infeasible )
            *result = SCIP_CUTOFF;
         else if ( reduceddom )
            *result = SCIP_REDUCEDDOM;
         else
         {
            *result = SCIP_FEASIBLE;
            SCIPwarningMessage(scip, "could not enforce feasibility by separating or branching; declaring solution with viol %g as feasible\n", maxviol);
         }
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMsg(scip, "Could not find any usual branching variable candidate. Proposed variable <%s> with LP value %g for branching.\n",
            SCIPvarGetName(brvar), SCIPgetSolVal(scip, sol, brvar));
         nnotify = 1;
      }
   }

   assert(*result == SCIP_INFEASIBLE && (solinfeasible || nnotify > 0));
   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBivariate)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   /* assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0); */

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBivariate(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip     != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprinterpreter != NULL);
   assert(conshdlrdata->exprgraph != NULL);
   assert(SCIPexprgraphGetNVars(conshdlrdata->exprgraph) == 0);

   /* free expression graph */
   SCIP_CALL( SCIPexprgraphFree(&conshdlrdata->exprgraph) );

   if( conshdlrdata->exprinterpreter != NULL )
   {
      SCIP_CALL( SCIPexprintFree(&conshdlrdata->exprinterpreter) );
   }

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreBivariate)
{  /*lint --e{715}*/
   SCIP_CONSDATA*     consdata;
   int                c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* reset may{in,de}creasez to FALSE in case some values are still set from a previous solve round */
      consdata->mayincreasez = FALSE;
      consdata->maydecreasez = FALSE;

      /* mark the constraint to be propagated */
      SCIP_CALL( SCIPmarkConsPropagate(scip, conss[c]) );  /*lint !e613*/
   }

   return SCIP_OKAY;
}

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                c;
   SCIP_Bool          changed;
   SCIP_Bool          upgraded;
#ifndef NDEBUG
   SCIP_CONSDATA*     consdata;
#endif

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdata->isremovedfixings )
   {
      SCIP_CALL( removeFixedNonlinearVariables(scip, conshdlr) );
      assert(conshdlrdata->isremovedfixings);
      /* @todo call expression graph simplifier? */
   }

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);  /* for flexelint */
      assert(conss[c] != NULL);

      /* make sure variable fixations have been resolved */
      SCIP_CALL( removeFixedVariables(scip, conshdlr, conss[c], &changed, &upgraded) );
      assert(!upgraded);

#ifndef NDEBUG
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      assert(consdata->f != NULL);
      assert(SCIPexprtreeGetNVars(consdata->f) == 2);
      assert(consdata->z == NULL || SCIPvarIsActive(consdata->z) || SCIPvarGetStatus(consdata->z) == SCIP_VARSTATUS_MULTAGGR);
#endif

      /* tell SCIP that we have something nonlinear */
      if( SCIPconsIsAdded(conss[c]) )
         SCIPenableNLP(scip);
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c;
#ifdef TYPESTATISTICS
   int                nconvextypeslhs[(int)SCIP_BIVAR_UNKNOWN+1];
   int                nconvextypesrhs[(int)SCIP_BIVAR_UNKNOWN+1];
#endif

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

#ifdef TYPESTATISTICS
   BMSclearMemoryArray(nconvextypeslhs, (int)SCIP_BIVAR_UNKNOWN+1);
   BMSclearMemoryArray(nconvextypesrhs, (int)SCIP_BIVAR_UNKNOWN+1);
#endif

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* check if linear variable can be rounded up or down without harming other constraints */
      if( consdata->z != NULL )
      {
         int poslock;
         int neglock;

         if( consdata->zcoef > 0.0 )
         {
            poslock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
            neglock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
         }
         else
         {
            poslock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
            neglock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
         }

         if( SCIPvarGetNLocksDown(consdata->z) - neglock == 0 )
         {
            /* for c*z + f(x,y) \in [lhs, rhs], we can decrease z without harming other constraints */
            consdata->maydecreasez = TRUE;
            SCIPdebugMsg(scip, "may decrease <%s> to become feasible\n", SCIPvarGetName(consdata->z));
         }

         if( SCIPvarGetNLocksDown(consdata->z) - poslock == 0 )
         {
            /* for c*x + f(x,y) \in [lhs, rhs], we can increase x without harming other constraints */
            consdata->mayincreasez = TRUE;
            SCIPdebugMsg(scip, "may increase <%s> to become feasible\n", SCIPvarGetName(consdata->z));
         }
      }

      /* add nlrow respresentation to NLP, if NLP had been constructed */
      if( SCIPisNLPConstructed(scip) && SCIPconsIsEnabled(conss[c]) )  /*lint !e613*/
      {
         SCIP_NLROW* nlrow;

         SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[c]), 0.0,
               consdata->z != NULL ? 1 : 0, consdata->z != NULL ? &consdata->z : NULL, &consdata->zcoef,
               0, NULL, 0, NULL,
               consdata->f, consdata->lhs, consdata->rhs,
               consdata->convextype == SCIP_BIVAR_ALLCONVEX ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_UNKNOWN) );  /*lint !e826 !e613*/

         SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
         SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
      }

      /* initialize data for cut generation */
      SCIP_CALL( initSepaData(scip, conshdlrdata->exprinterpreter, conss[c]) );  /*lint !e613*/

#ifdef TYPESTATISTICS
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         ++nconvextypeslhs[consdata->convextype];
      if( !SCIPisInfinity(scip,  consdata->rhs) )
         ++nconvextypesrhs[consdata->convextype];
#endif
   }

   conshdlrdata->newsoleventfilterpos = -1;
   if( nconss != 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );

#ifdef TYPESTATISTICS
      for( c = 0; c <= (int)SCIP_BIVAR_UNKNOWN; ++c )
      {
         const char* typename;
         switch( c )
         {
         case SCIP_BIVAR_ALLCONVEX:
            typename = "allconvex";
            break;
         case SCIP_BIVAR_1CONVEX_INDEFINITE:
            typename = "1-convex";
            break;
         case SCIP_BIVAR_CONVEX_CONCAVE:
            typename = "convex-concave";
            break;
         case SCIP_BIVAR_UNKNOWN:
         default:
            typename = "unknown";
            break;
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%4d left and %4d right bivariate constraints of type [%s]\n", nconvextypeslhs[c], nconvextypesrhs[c], typename);
      }
#endif
   }

   /* reset counter */
   conshdlrdata->lastenfonode = NULL;
   conshdlrdata->nenforounds = 0;

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->newsoleventfilterpos >= 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, conshdlrdata->newsoleventfilterpos) );
      conshdlrdata->newsoleventfilterpos = -1;
   }

   for( c = 0; c < nconss; ++c )
   {
      /* free data for cut generation */
      assert(conss[c] != NULL);  /*lint !e613*/

      SCIP_CALL( freeSepaData(scip, conss[c]) );  /*lint !e613*/
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteBivariate)
{  /*lint --e{715}*/
#ifndef NDEBUG
   SCIP_CONSHDLRDATA* conshdlrdata;
#endif

   assert(scip != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);

#ifndef NDEBUG
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
#endif

   /* expression should have been removed from expression graph when constraint was deactivated */
   assert((*consdata)->exprgraphnode == NULL);

   if( (*consdata)->f != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&(*consdata)->f) );
   }

   SCIPfreeBlockMemory(scip, consdata);
   *consdata = NULL;

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransBivariate)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   SCIP_VAR* targetvars[2];

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( SCIPduplicateBlockMemory(scip, &targetdata, sourcedata) );
   assert(targetdata->eventfilterpos == -1);

   assert(sourcedata->f != NULL);
   SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &targetdata->f, sourcedata->f) );
   SCIP_CALL( SCIPgetTransformedVars(scip, 2, SCIPexprtreeGetVars(sourcedata->f), targetvars) );
   SCIP_CALL( SCIPexprtreeSetVars(targetdata->f, 2, targetvars) );

   if( sourcedata->z != NULL )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->z, &targetdata->z) );
   }

   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_ROW*          row1;
   SCIP_ROW*          row2;
   SCIP_Real          xy[2];
   int                c;
   int                i;
   int                ix;
   int                iy;
   int                nref;
   SCIP_Real          lb[2];
   SCIP_Real          ub[2];
   SCIP_Bool          unbounded[2];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *infeasible = FALSE;

   nref = conshdlrdata->ninitlprefpoints;

   if( nref == 0 )
   {
      SCIPdebugMsg(scip, "skip LP initialization since ninitlprefpoints is 0\n");
      return SCIP_OKAY;
   }

   row1 = NULL;
   row2 = NULL;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);
      assert(consdata->f != NULL);

      if( SCIPexprtreeGetInterpreterData(consdata->f) == NULL )
      {
         SCIP_CALL( SCIPexprintCompile(conshdlrdata->exprinterpreter, consdata->f) );
      }

      /* create a bounded rectangle in which we take reference points for initial cut generation
       * For a missing bound, we either reflect the other bound at 0.0 if finite and on the right side,
       * or double the other bound if on the same side but not 0.0, or set it to +/-1000.0.
       */
      for( i = 0; i < 2; ++i )
      {
         lb[i] = SCIPvarGetLbGlobal(SCIPexprtreeGetVars(consdata->f)[i]);
         ub[i] = SCIPvarGetUbGlobal(SCIPexprtreeGetVars(consdata->f)[i]);

         unbounded[i] = FALSE;
         if( SCIPisInfinity(scip, -lb[i]) )
         {
            unbounded[i] = TRUE;
            ub[i] = MIN(INITLPMAXVARVAL, ub[i]);
            if( SCIPisPositive(scip, ub[i]) )
               lb[i] = -ub[i];
            else if( SCIPisZero(scip, ub[i]) )
               lb[i] = -INITLPMAXVARVAL;
            else
               lb[i] = 2.0 * ub[i];
         }
         else if( SCIPisInfinity(scip, ub[i]) )
         {
            unbounded[i] = TRUE;
            assert(!SCIPisInfinity(scip, -lb[i]));
            lb[i] = MAX(-INITLPMAXVARVAL, lb[i]);
            if( SCIPisNegative(scip, lb[i]) )
               ub[i] = -lb[i];
            else if( SCIPisZero(scip, lb[i]) )
               ub[i] = INITLPMAXVARVAL;
            else
               ub[i] = 2.0 * lb[i];
         }
      }

      for( ix = 0; ix < nref; ++ix )
      {
         if( nref > 1 )
            xy[0] = lb[0] + ix * (ub[0] - lb[0]) / (nref - 1.0);
         else
            xy[0] = (lb[0] + ub[0]) / 2.0;

         for( iy = 0; iy < nref; ++iy )
         {
            if( nref > 1 )
               xy[1] = lb[1] + iy * (ub[1] - lb[1]) / (nref - 1.0);
            else
               xy[1] = (lb[1] + ub[1]) / 2.0;

            SCIPdebugMsg(scip, "cons <%s>: generate cuts for <%s> = %g [%g,%g], <%s> = %g [%g,%g]\n",
               SCIPconsGetName(conss[c]),  /*lint !e613*/
               SCIPvarGetName(SCIPexprtreeGetVars(consdata->f)[0]), xy[0],
               SCIPvarGetLbGlobal(SCIPexprtreeGetVars(consdata->f)[0]), SCIPvarGetUbGlobal(SCIPexprtreeGetVars(consdata->f)[0]),
               SCIPvarGetName(SCIPexprtreeGetVars(consdata->f)[1]), xy[1],
               SCIPvarGetLbGlobal(SCIPexprtreeGetVars(consdata->f)[1]), SCIPvarGetUbGlobal(SCIPexprtreeGetVars(consdata->f)[1])
               );

            /* try to generate one cut for each side */
            switch( consdata->convextype )
            {
            case SCIP_BIVAR_ALLCONVEX:
            {
               if( !SCIPisInfinity(scip, -consdata->lhs) && !unbounded[0] && !unbounded[1] && (ix == 0 || ix == nref-1) && (iy == 0 || iy == nref-1) )
               {
                  /* lhs is finite and both variables are bounded, so can do overest. hyperplane
                   * do this only for corner points, since we can get at most two cuts out of it
                   * @todo generate only two cuts instead of four
                   */
                  SCIP_CALL( generateOverestimatingHyperplaneCut(scip, conshdlrdata->exprinterpreter, conss[c], xy, &row1) );  /*lint !e613*/
               }
               if( !SCIPisInfinity(scip,  consdata->rhs) )
               {
                  /* rhs is finite */
                  SCIP_CALL( generateLinearizationCut(scip, conshdlrdata->exprinterpreter, conss[c], xy, TRUE, &row2) );  /*lint !e613*/
               }
               break;
            }

            case SCIP_BIVAR_CONVEX_CONCAVE:
            {
               if( !SCIPisInfinity(scip, -consdata->lhs) && !unbounded[0])
               {
                  /* lhs is finite and x is bounded */
                  SCIP_CALL( generateConvexConcaveEstimator(scip, conshdlrdata->exprinterpreter, conss[c], xy, SCIP_SIDETYPE_LEFT, &row1) );  /*lint !e613*/
               }
               if( !SCIPisInfinity(scip,  consdata->rhs) && !unbounded[1])
               {
                  /* rhs is finite and y is bounded */
                  SCIP_CALL( generateConvexConcaveEstimator(scip, conshdlrdata->exprinterpreter, conss[c], xy, SCIP_SIDETYPE_RIGHT, &row2) );  /*lint !e613*/
               }
               break;
            }

            case SCIP_BIVAR_1CONVEX_INDEFINITE:
            {
               if( !SCIPisInfinity(scip, -consdata->lhs) && !unbounded[0] && !unbounded[1] && (ix == 0 || ix == nref-1) && (iy == 0 || iy == nref-1) )
               {
                  /* lhs is finite and both variables are bounded
                   * do this only for corner points, since we can get at most two cuts out of it
                   * @todo generate only two cuts instead of four
                   */
                  SCIP_CALL( generateOverestimatingHyperplaneCut(scip, conshdlrdata->exprinterpreter, conss[c], xy, &row1) );  /*lint !e613*/
               }
               if( !SCIPisInfinity(scip,  consdata->rhs) && !unbounded[0] && !unbounded[1] )
               { /* rhs is finite and both variables are bounded */
                  SCIP_CALL( generate1ConvexIndefiniteUnderestimator(scip, conshdlrdata->exprinterpreter, conss[c], xy, &row2) );  /*lint !e613*/
               }
               break;
            }

            default:
            {
               SCIPwarningMessage(scip, "initlp for convexity type %d not implemented\n", consdata->convextype);
            }
            }  /*lint !e788*/

            /* check numerics */
            if( row1 != NULL )
            {
               if( SCIPgetRowMaxCoef(scip, row1) / SCIPgetRowMinCoef(scip, row1) > conshdlrdata->cutmaxrange )
               {
                  SCIPdebugMsg(scip, "drop row1 for constraint <%s> because range of coefficients is too large: mincoef = %g, maxcoef = %g -> range = %g\n",
                     SCIPconsGetName(conss[c]), SCIPgetRowMinCoef(scip, row1), SCIPgetRowMaxCoef(scip, row1), SCIPgetRowMaxCoef(scip, row1) / SCIPgetRowMinCoef(scip, row1));  /*lint !e613*/
               }
               else if( SCIPisInfinity(scip, -SCIProwGetLhs(row1)) )
               {
                  /* row1 should be a cut with finite lhs, but infinite rhs */
                  assert(SCIPisInfinity(scip, SCIProwGetRhs(row1)));
                  SCIPdebugMsg(scip, "drop row1 for constraint <%s> because of very large lhs: %g\n", SCIPconsGetName(conss[c]), SCIProwGetLhs(row1));  /*lint !e613*/
               }
               /* add row to LP */
               else
               {
                  SCIP_CALL( SCIPaddRow(scip, row1, FALSE /* forcecut */, infeasible) );
                  SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row1, NULL) ) );
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row1) );
            }

            if( row2 != NULL )
            {
               if( SCIPgetRowMaxCoef(scip, row2) / SCIPgetRowMinCoef(scip, row2) > conshdlrdata->cutmaxrange )
               {
                  SCIPdebugMsg(scip, "drop row2 for constraint <%s> because range of coefficients is too large: mincoef = %g, maxcoef = %g -> range = %g\n",
                     SCIPconsGetName(conss[c]), SCIPgetRowMinCoef(scip, row2), SCIPgetRowMaxCoef(scip, row2), SCIPgetRowMaxCoef(scip, row2) / SCIPgetRowMinCoef(scip, row2));  /*lint !e613*/
               }
               else if( SCIPisInfinity(scip, SCIProwGetRhs(row2)) )
               {
                  /* row2 should be a cut with finite rhs, but infinite lhs */
                  assert(SCIPisInfinity(scip, SCIProwGetRhs(row2)));
                  SCIPdebugMsg(scip, "drop row2 for constraint <%s> because of very large rhs: %g\n", SCIPconsGetName(conss[c]), SCIProwGetLhs(row2));  /*lint !e613*/
               }
               /* add row to LP */
               else if( !(*infeasible) )
               {
                  SCIP_CALL( SCIPaddRow(scip, row2, FALSE /* forcecut */, infeasible) );
                  SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row2, NULL) ) );
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row2) );
            }

            if( *infeasible )
               return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpBivariate)
{  /*lint --e{715}*/
   SCIP_CONS*         maxviolcon;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, NULL, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   /* @todo add separation of convex (only?) constraints in nlp relaxation solution */

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, SCIPgetSepaMinEfficacy(scip), FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolBivariate)
{  /*lint --e{715}*/
   SCIP_CONS*         maxviolcon;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(sol      != NULL);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, sol, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, SCIPgetSepaMinEfficacy(scip), FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBivariate)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, solinfeasible, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxBivariate)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, solinfeasible, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBivariate)
{  /*lint --e{715}*/
   SCIP_CONS*         maxviolcons;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        propresult;
   SCIP_VAR*          var;
   int                nnotify;
   int                dummy;
   int                c;
   int                i;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, NULL, &maxviolcons) );
   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   SCIPdebugMsg(scip, "enfops with max violation in cons <%s>\n", SCIPconsGetName(maxviolcons));

   /* run domain propagation */
   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, &dummy, &dummy) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we are not feasible and we cannot proof that the whole node is infeasible
    * -> collect all variables in violated constraints for branching
    */

   nnotify = 0;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->f != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      /* if nonlinear variables are fixed, z should be propagated such that the constraint becomes feasible,
       * so there should be no branching on z necessary
       */
      if( consdata->z != NULL && !SCIPisRelEQ(scip, SCIPvarGetLbLocal(consdata->z), SCIPvarGetUbLocal(consdata->z)) )
      {
         SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->z, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
         ++nnotify;
      }

      for( i = 0; i < 2; ++i )
      {
         var = SCIPexprtreeGetVars(consdata->f)[i];
         if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++nnotify;
         }
      }
   }

   if( nnotify == 0 )
   {
      SCIPdebugMsg(scip, "All variables in violated constraints fixed (up to epsilon). Cannot find branching candidate. Forcing solution of LP.\n");
      *result = SCIP_SOLVELP;
   }

   assert(*result == SCIP_SOLVELP || (*result == SCIP_INFEASIBLE && nnotify > 0));
   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   int                c;
   SCIP_Bool          maypropfeasible; /* whether we may be able to propose a feasible solution */

   assert(scip   != NULL);
   assert(conss  != NULL || nconss == 0);
   assert(sol    != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   maxviol = 0.0;
   maypropfeasible = conshdlrdata->linfeasshift && (conshdlrdata->trysolheur != NULL);
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( computeViolation(scip, conshdlr, conss[c], sol) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         *result = SCIP_INFEASIBLE;
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g (scaled: %.15g)\n", consdata->lhs - consdata->activity, consdata->lhsviol);
            }
            if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g (scaled: %.15g)\n", consdata->activity - consdata->rhs, consdata->rhsviol);
            }
         }

         if( (conshdlrdata->subnlpheur == NULL || sol == NULL) && !maypropfeasible && !completely )
            return SCIP_OKAY;

         if( consdata->lhsviol > maxviol || consdata->rhsviol > maxviol )
            maxviol = consdata->lhsviol + consdata->rhsviol;

         /* do not try to shift linear variables if activity is at infinity (leads to setting variable to infinity in solution, which is not allowed) */
         if( maypropfeasible && SCIPisInfinity(scip, REALABS(consdata->activity)) )
            maypropfeasible = FALSE;

         if( maypropfeasible )
         {
            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               /* check if the linear variable may help to get the left hand side satisfied
                * if not, then we cannot get feasible */
               if( !(consdata->mayincreasez && consdata->zcoef > 0.0) && !(consdata->maydecreasez && consdata->zcoef < 0.0) )
                  maypropfeasible = FALSE;
            }
            else
            {
               assert(SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)));
               /* check if the linear variable may help to get the right hand side satisfied
                * if not, then we cannot get feasible */
               if( !(consdata->mayincreasez && consdata->zcoef < 0.0) && !(consdata->maydecreasez && consdata->zcoef > 0.0) )
                  maypropfeasible = FALSE;
            }
         }
      }
   }

   if( *result == SCIP_INFEASIBLE && maypropfeasible )
   {
      SCIP_Bool success;

      SCIP_CALL( proposeFeasibleSolution(scip, conshdlr, conss, nconss, sol, &success) );

      /* do not pass solution to NLP heuristic if we made it feasible this way */
      if( success )
         return SCIP_OKAY;
   }

   if( *result == SCIP_INFEASIBLE && conshdlrdata->subnlpheur != NULL && sol != NULL && !SCIPisInfinity(scip, maxviol) )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropBivariate)
{  /*lint --e{715}*/
   int dummy;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, result, &dummy, &dummy) );

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolBivariate)
{  /*lint --e{715}*/
#ifndef NDEBUG
   SCIP_CONSDATA*     consdata;
#endif
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT        propresult;
   SCIP_Bool          havechange;
   SCIP_Bool          upgraded;
   int                c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   if( !conshdlrdata->isremovedfixings )
   {
      SCIP_CALL( removeFixedNonlinearVariables(scip, conshdlr) );
      assert(conshdlrdata->isremovedfixings);
   }
   /* @todo call expression graph simplifier, if not done yet? */

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);

#ifndef NDEBUG
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
#endif

      SCIPdebugMsg(scip, "process constraint <%s>\n", SCIPconsGetName(conss[c]));
      SCIPdebugPrintCons(scip, conss[c], NULL);

      havechange = FALSE;

      SCIP_CALL( removeFixedVariables(scip, conshdlr, conss[c], &havechange, &upgraded) );
      if( upgraded )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         ++*nupgdconss;
         continue;
      }
      if( havechange )
      {
         SCIPdebugMsg(scip, "removed fixed variables -> ");
         SCIPdebugPrintCons(scip, conss[c], NULL);
      }
   }

   /* run domain propagation */
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, nchgbds, ndelconss) );
   switch( propresult )
   {
   case SCIP_REDUCEDDOM:
      *result = SCIP_SUCCESS;
      break;
   case SCIP_CUTOFF:
      SCIPdebugMsg(scip, "propagation says problem is infeasible in presolve\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   default:
      assert(propresult == SCIP_DIDNOTFIND || propresult == SCIP_DIDNOTRUN);
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBivariate)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->z != NULL )
   {
      if( consdata->zcoef > 0 )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->z, nlockspos, nlocksneg) );
         }
         if( !SCIPisInfinity(scip,  consdata->rhs) )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->z, nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->z, nlocksneg, nlockspos) );
         }
         if( !SCIPisInfinity(scip,  consdata->rhs) )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->z, nlockspos, nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool exprtreeisnew;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->exprgraphnode == NULL);

   SCIPdebugMsg(scip, "activate %scons <%s>\n", SCIPconsIsTransformed(cons) ? "transformed " : "", SCIPconsGetName(cons));

   /* add exprtree to expression graph */
   SCIP_CALL( SCIPexprgraphAddExprtreeSum(conshdlrdata->exprgraph, 1, &consdata->f, NULL, &consdata->exprgraphnode, &exprtreeisnew) );
   assert(consdata->exprgraphnode != NULL);

   /* mark that variables in constraint should not be multiaggregated (bad for bound tightening and branching) */
   if( SCIPvarIsActive(SCIPexprtreeGetVars(consdata->f)[0]) )
   {
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, SCIPexprtreeGetVars(consdata->f)[0]) );
   }
   if( SCIPvarIsActive(SCIPexprtreeGetVars(consdata->f)[1]) )
   {
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, SCIPexprtreeGetVars(consdata->f)[1]) );
   }
   if( consdata->z != NULL && SCIPvarIsActive(consdata->z) )
   {
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, consdata->z) );
   }

   return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->exprgraphnode != NULL);

   SCIPdebugMsg(scip, "deactivate %scons <%s>\n", SCIPconsIsTransformed(cons) ? "transformed " : "", SCIPconsGetName(cons));

   SCIP_CALL( SCIPexprgraphReleaseNode(conshdlrdata->exprgraph, &consdata->exprgraphnode) );

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));
   assert(SCIPconsIsActive(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->exprgraphnode != NULL);

   SCIPdebugMsg(scip, "enable %scons <%s>\n", SCIPconsIsTransformed(cons) ? "transformed " : "", SCIPconsGetName(cons));

   /* enable node of expression in expression graph */
   SCIPexprgraphEnableNode(conshdlrdata->exprgraph, consdata->exprgraphnode);

   /* enable event catching for linear variables */
   SCIP_CALL( catchLinearVarEvents(scip, cons) );

   return SCIP_OKAY;
}

/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableBivariate)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->exprgraphnode != NULL);

   SCIPdebugMsg(scip, "disable %scons <%s>\n", SCIPconsIsTransformed(cons) ? "transformed " : "", SCIPconsGetName(cons));

   /* disable node of expression in expression graph */
   SCIPexprgraphDisableNode(conshdlrdata->exprgraph, consdata->exprgraphnode);

   SCIP_CALL( dropLinearVarEvents(scip, cons) );

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintBivariate)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   SCIP_CALL( SCIPexprtreePrintWithNames(consdata->f, SCIPgetMessagehdlr(scip), file) );

   if( consdata->z != NULL )
   {
      SCIPinfoMessage(scip, file, "%+.15g", consdata->zcoef);
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->z, TRUE) );
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   }
   else
   {
      SCIPinfoMessage(scip, file, " [free]");
   }

   /* print convexity type, if known */
   switch( consdata->convextype )
   {
   case SCIP_BIVAR_ALLCONVEX:
      SCIPinfoMessage(scip, file, " [allconvex]");
      break;
   case SCIP_BIVAR_1CONVEX_INDEFINITE:
      SCIPinfoMessage(scip, file, " [1-convex]");
      break;
   case SCIP_BIVAR_CONVEX_CONCAVE:
      SCIPinfoMessage(scip, file, " [convex-concave]");
      break;
   default: ;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyBivariate)
{  /*lint --e{715}*/
   SCIP_CONSDATA*    consdata;
   SCIP_EXPRTREE*    f;
   SCIP_VAR*         xy[2];
   SCIP_VAR*         z;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourceconshdlr != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);

   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);
   assert(consdata->f != NULL);

   *valid = TRUE;

   if( consdata->z != NULL )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->z, &z, varmap, consmap, global, valid) );
      assert(!*valid || z != NULL);
   }
   else
      z = NULL;

   if( *valid )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, SCIPexprtreeGetVars(consdata->f)[0], &xy[0], varmap, consmap, global, valid) );
      assert(!*valid || xy[0] != NULL);
   }

   if( *valid )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, SCIPexprtreeGetVars(consdata->f)[1], &xy[1], varmap, consmap, global, valid) );
      assert(!*valid || xy[1] != NULL);
   }

   if( *valid )
   {
      SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &f, consdata->f) );
      SCIP_CALL( SCIPexprtreeSetVars(f, 2, xy) );
   }
   else
      f = NULL;

   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsBivariate(scip, cons, name ? name : SCIPconsGetName(sourcecons),
            f, consdata->convextype, z, consdata->zcoef, consdata->lhs, consdata->rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   if( f != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&f) );
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsBivariate)
{  /*lint --e{715}*/

   if( varssize < 3 )
      (*success) = FALSE;
   else
   {
      SCIP_CONSDATA* consdata;

      assert(cons != NULL);
      assert(vars != NULL);

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      vars[0] = SCIPexprtreeGetVars(consdata->f)[0];
      vars[1] = SCIPexprtreeGetVars(consdata->f)[1];
      vars[2] = consdata->z;
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsBivariate)
{  /*lint --e{715}*/

   (*nvars) = 3;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * Quadratic constraint upgrading
 */

/** tries to upgrade a quadratic constraint into a bivariate constraint */
static
SCIP_DECL_QUADCONSUPGD(quadconsUpgdBivariate)
{  /*lint --e{715}*/
   SCIP_QUADVARTERM* quadvarterms;
   SCIP_BILINTERM* bilinterms;
   int nquadvarterms;
   int nbilinterms;

   SCIP_VAR* x;
   SCIP_VAR* y;

   SCIP_Real coefxx;
   SCIP_Real coefxy;
   SCIP_Real coefyy;
   SCIP_Real coefx;
   SCIP_Real coefy;

   SCIP_Real zcoef;
   SCIP_VAR* z;

   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   /* not interested in univariate case */
   if( nbinquad + nintquad + ncontquad < 2 )
      return SCIP_OKAY;

   if( SCIPgetNBilinTermsQuadratic(scip, cons) == 0 )
      return SCIP_OKAY;

   quadvarterms  = SCIPgetQuadVarTermsQuadratic(scip, cons);
   nquadvarterms = SCIPgetNQuadVarTermsQuadratic(scip, cons);
   bilinterms  = SCIPgetBilinTermsQuadratic(scip, cons);
   nbilinterms = SCIPgetNBilinTermsQuadratic(scip, cons);

   if( nquadvarterms == 2 && SCIPgetNLinearVarsQuadratic(scip, cons) <= 1 )
   {
      x = quadvarterms[0].var;
      y = quadvarterms[1].var;

      coefxx = quadvarterms[0].sqrcoef;
      coefyy = quadvarterms[1].sqrcoef;

      /* only one bilinear term -> not interesting for us */
      if( coefxx == 0.0 && coefyy == 0.0 )
         return SCIP_OKAY;

      /* two square terms without bilinear term -> also not interesting for us */
      if( nbilinterms == 0 )
         return SCIP_OKAY;
      assert(nbilinterms == 1);

      assert(bilinterms[0].var1 == x || bilinterms[0].var1 == y);
      assert(bilinterms[0].var2 == x || bilinterms[0].var2 == y);

      coefxy = bilinterms[0].coef;

      coefx = quadvarterms[0].lincoef;
      coefy = quadvarterms[1].lincoef;

      if( SCIPgetNLinearVarsQuadratic(scip, cons) )
      {
         assert(SCIPgetNLinearVarsQuadratic(scip, cons) == 1);
         zcoef = SCIPgetCoefsLinearVarsQuadratic(scip, cons)[0];
         z = SCIPgetLinearVarsQuadratic(scip, cons)[0];
      }
      else
      {
         z = NULL;
         zcoef = 0.0;
      }

      if( upgdconsssize < 1 )
      {
         *nupgdconss = -1;
         return SCIP_OKAY;
      }

      SCIP_CALL( createConsFromQuadTerm(scip, cons, &upgdconss[0], SCIPconsGetName(cons),
            x, y, z, coefxx, coefx, coefyy, coefy, coefxy, zcoef, SCIPgetLhsQuadratic(scip, cons), SCIPgetRhsQuadratic(scip, cons)) );
      *nupgdconss = 1;
   }
   else
   {
      SCIP_CONS* quadcons;
      SCIP_Bool upgdlhs;
      SCIP_Bool upgdrhs;
      SCIP_Bool keeporig;
      SCIP_Bool* marked;
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR* auxvar;
      int xpos;
      int ypos;
      int pos;
      int i;

      /* needs to check curvature, which might be expensive */
      if( (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 && nquadvarterms > 10 )
         return SCIP_OKAY;
      if( (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 && nquadvarterms > 50 )
         return SCIP_OKAY;

      /* check if we find at least one bilinear term for which we would create a bivariate constraint
       * thus, we search for a variable that has a square term and is involved in at least one bivariate term */
      for( i = 0; i < nquadvarterms; ++i )
         if( quadvarterms[i].sqrcoef != 0.0 && quadvarterms[i].nadjbilin > 0 )
            break;

      /* if nothing found, then don't try upgrade and return */
      if( i == nquadvarterms )
         return SCIP_OKAY;

      /* check which constraint side we want to upgrade and whether to keep some
       * we want to upgrade those that are nonconvex */
      SCIP_CALL( SCIPcheckCurvatureQuadratic(scip, cons) );
      upgdlhs = FALSE;
      upgdrhs = FALSE;
      keeporig = FALSE;
      if( !SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) )
      {
         if( SCIPisConcaveQuadratic(scip, cons) )
            keeporig = TRUE;
         else
            upgdlhs = TRUE;
      }
      if( !SCIPisInfinity(scip,  SCIPgetRhsQuadratic(scip, cons)) )
      {
         if( SCIPisConvexQuadratic(scip, cons) )
            keeporig = TRUE;
         else
            upgdrhs = TRUE;
      }

      /* if nothing to upgrade, then return */
      if( !upgdlhs && !upgdrhs )
         return SCIP_OKAY;

      /* require enough space here already, so we do not create and add aux vars that we cannot get rid of easily later */
      if( upgdconsssize < nbilinterms + 1 + (keeporig ? 1 : 0) )
      {
         *nupgdconss = -(nbilinterms + 1 + (keeporig ? 1 : 0));
         return SCIP_OKAY;
      }

      /* initial remaining quadratic constraint: take linear part and constraint sides from original constraint */
      SCIP_CALL( SCIPcreateConsQuadratic(scip, &quadcons, SCIPconsGetName(cons),
            SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
            0, NULL, NULL, NULL,
            upgdlhs ? SCIPgetLhsQuadratic(scip, cons) : -SCIPinfinity(scip),
            upgdrhs ? SCIPgetRhsQuadratic(scip, cons) :  SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );

      /* remember for each quadratic variable whether its linear and square part has been moved into a bivariate constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &marked, nquadvarterms) );
      BMSclearMemoryArray(marked, SCIPgetNQuadVarTermsQuadratic(scip,cons));

      /* @todo what is a good partition of a number of quadratic terms into bivariate terms? */

      /* check for each bilinear term, whether we want to create a bivariate constraint for it and associated square terms */
      for( i = 0; i < nbilinterms; ++i )
      {
         assert(bilinterms[i].coef != 0.0);

         x = bilinterms[i].var1;
         y = bilinterms[i].var2;

         SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, cons, x, &xpos) );
         assert(xpos >= 0);
         assert(xpos < nquadvarterms);
         assert(quadvarterms[xpos].var == x);

         SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, cons, y, &ypos) );
         assert(ypos >= 0);
         assert(ypos < nquadvarterms);
         assert(quadvarterms[ypos].var == y);

         coefxx = marked[xpos] ? 0.0 : quadvarterms[xpos].sqrcoef;
         coefyy = marked[ypos] ? 0.0 : quadvarterms[ypos].sqrcoef;

         /* if there are no square terms, then do not upgrade bilinear term to bivariate constraint
          * thus, add bivariate term to quadcons and continue
          */
         if( coefxx == 0.0 && coefyy == 0.0 )
         {
            /* check if x and y already are in quadcons and add if not there yet */
            SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, quadcons, x, &pos) );
            if( pos == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarQuadratic(scip, quadcons, x, 0.0, 0.0) );
            }
            SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, quadcons, y, &pos) );
            if( pos == -1 )
            {
               SCIP_CALL( SCIPaddQuadVarQuadratic(scip, quadcons, y, 0.0, 0.0) );
            }

            SCIP_CALL( SCIPaddBilinTermQuadratic(scip, quadcons, x, y, bilinterms[i].coef) );

            continue;
         }

         coefx = marked[xpos] ? 0.0 : quadvarterms[xpos].lincoef;
         coefy = marked[ypos] ? 0.0 : quadvarterms[ypos].lincoef;
         coefxy = bilinterms[i].coef;

         /* create new auxiliary variable for bilinear quad. term in x and y */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_auxvar%d", SCIPconsGetName(cons), *nupgdconss);
         SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, auxvar) );

         /* add 1*auxvar to quadcons */
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, quadcons, auxvar, 1.0) );

         /* create new bivariate constraint */
         assert(*nupgdconss < upgdconsssize);
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_auxcons%d", SCIPconsGetName(cons), *nupgdconss);
         SCIP_CALL( createConsFromQuadTerm(scip, cons, &upgdconss[*nupgdconss], name,
               x, y, auxvar, coefxx, coefx, coefyy, coefy, coefxy, -1.0,
               SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) ? -SCIPinfinity(scip) : 0.0,
               SCIPisInfinity(scip,  SCIPgetRhsQuadratic(scip, cons)) ?  SCIPinfinity(scip) : 0.0) );
         /* need to enforce new constraints, as relation auxvar = f(x,y) is not redundant, even if original constraint is */
         SCIP_CALL( SCIPsetConsEnforced(scip, upgdconss[*nupgdconss], TRUE) );
         SCIP_CALL( SCIPsetConsChecked(scip, upgdconss[*nupgdconss], TRUE) );
         ++*nupgdconss;

         /* compute value of auxvar in debug solution */
#ifdef WITH_DEBUG_SOLUTION
         if( SCIPdebugIsMainscip(scip) )
         {
            SCIP_Real xval;
            SCIP_Real yval;
            SCIP_CALL( SCIPdebugGetSolVal(scip, x, &xval) );
            SCIP_CALL( SCIPdebugGetSolVal(scip, y, &yval) );
            SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, coefxx * xval * xval + coefyy * yval * yval + coefxy * xval * yval + coefx * xval + coefy * yval) );
         }
#endif

         SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

         marked[xpos] = TRUE;
         marked[ypos] = TRUE;
      }

      if( *nupgdconss == 0 )
      {
         /* if no constraints created, then forget also quadcons and do no upgrade */
         SCIP_CALL( SCIPreleaseCons(scip, &quadcons) );
      }
      else
      {
         /* complete quadcons: check for unmarked quadvarterms and add their linear and square coefficients to quadcons */
         for( i = 0; i < nquadvarterms; ++i )
         {
            if( marked[i] )
               continue;

            x = quadvarterms[i].var;

            /* check if variable is already in quadcons
             * if the variable appears in a bilinear term, then this term should have been added to quadcons above, so the variable is there
             */
            pos = -1;
            if( quadvarterms[i].nadjbilin > 0 )
            {
               SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, quadcons, x, &pos) );
            }

            /* create new quad var or add existing quad var */
            if( quadvarterms[i].sqrcoef != 0.0 )
            {
               if( pos == -1 )
               {
                  SCIP_CALL( SCIPaddQuadVarQuadratic(scip, quadcons, x, quadvarterms[i].lincoef, quadvarterms[i].sqrcoef) );
               }
               else
               {
                  SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quadcons, x, quadvarterms[i].sqrcoef) );
                  SCIP_CALL( SCIPaddQuadVarLinearCoefQuadratic(scip, quadcons, x, quadvarterms[i].lincoef) );
               }
            }
            else if( quadvarterms[i].lincoef != 0.0 )
            {
               /* if no square term and no quadratic variable term, then add to linear part */
               SCIP_CALL( SCIPaddLinearVarQuadratic(scip, quadcons, x, quadvarterms[i].lincoef) );
            }
         }

         /* add quadcons to set of upgrade constraints */
         assert(*nupgdconss < upgdconsssize);
         upgdconss[*nupgdconss] = quadcons;
         ++*nupgdconss;

         SCIPdebugPrintCons(scip, quadcons, NULL);

         if( keeporig )
         {
            assert(*nupgdconss < upgdconsssize);
            /* copy of original quadratic constraint with one of the sides relaxed */
            SCIP_CALL( SCIPcreateConsQuadratic2(scip, &upgdconss[*nupgdconss], SCIPconsGetName(cons),
                  SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
                  SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
                  SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
                  upgdlhs ? -SCIPinfinity(scip) : SCIPgetLhsQuadratic(scip, cons),
                  upgdrhs ?  SCIPinfinity(scip) : SCIPgetRhsQuadratic(scip, cons),
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
            ++*nupgdconss;
         }
      }

      SCIPfreeBufferArray(scip, &marked);
   }

   return SCIP_OKAY;
}


/*
 * Nonlinear constraint upgrading
 */

/** tries to reformulate a expression graph node that is a monomial in two variables */
static
SCIP_DECL_EXPRGRAPHNODEREFORM(exprgraphnodeReformBivariate)
{
   SCIP_EXPRDATA_MONOMIAL* monomial;
   SCIP_CONS* cons;
   SCIP_VAR* auxvar;
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_Real expx;
   SCIP_Real expy;

   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(naddcons != NULL);
   assert(reformnode != NULL);

   *reformnode = NULL;

   /* could also upgrade bivariate quadratic, but if we don't then node will appear in cons_quadratic later, from which we also upgrade...
    * @todo could also upgrade x/y from EXPR_DIV */
   if( SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_POLYNOMIAL )
      return SCIP_OKAY;

   /* sums of monomials are split up by reformulation, so wait that this happened */
   if( SCIPexprgraphGetNodePolynomialNMonomials(node) != 1 )
      return SCIP_OKAY;

   /* we are only interested in monomials that are not convex or concave, since cons_nonlinear can handle these the same was as we do */
   if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
      return SCIP_OKAY;

   monomial = SCIPexprgraphGetNodePolynomialMonomials(node)[0];
   assert(monomial != NULL);

   /* @todo we could also do some more complex reformulation for n-variate monomials, something better than what reformMonomial in cons_nonlinear is doing */
   if( SCIPexprGetMonomialNFactors(monomial) != 2 )
      return SCIP_OKAY;
   assert(SCIPexprgraphGetNodeNChildren(node) == 2);

   expx = SCIPexprGetMonomialExponents(monomial)[0];
   expy = SCIPexprGetMonomialExponents(monomial)[1];

   /* no interest in upgrading x*y -> let cons_quadratic do this */
   if( SCIPisEQ(scip, expx, 1.0) && SCIPisEQ(scip, expy, 1.0) )
      return SCIP_OKAY;

   /* so far only support variables as arguments @todo could allow more here, e.g., f(x)^pg(y)^q */
   if( SCIPexprgraphGetNodeOperator(SCIPexprgraphGetNodeChildren(node)[0]) != SCIP_EXPR_VARIDX ||
      (SCIPexprgraphGetNodeOperator(SCIPexprgraphGetNodeChildren(node)[1]) != SCIP_EXPR_VARIDX) )
      return SCIP_OKAY;

   x = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[0]);
   y = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[1]);
   assert(x != y);

   /* so far only allow positive x and y @todo could also allow x<0 or y<0 */
   if( SCIPisNegative(scip, SCIPvarGetLbGlobal(x)) || SCIPisNegative(scip, SCIPvarGetLbGlobal(y)) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "reformulate bivariate monomial in node %p\n", (void*)node);

   /* create auxiliary variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%dbv", *naddcons);
   SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, SCIPexprgraphGetNodeBounds(node).inf, SCIPexprgraphGetNodeBounds(node).sup,
         0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar) );

   /* create bivariate constraint */
   SCIP_CALL( createConsFromMonomial(scip, NULL, &cons, name, x, y, auxvar,
         SCIPexprGetMonomialCoef(monomial), expx, expy, -1.0, -SCIPexprgraphGetNodePolynomialConstant(node), -SCIPexprgraphGetNodePolynomialConstant(node)) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   ++*naddcons;

   /* add auxvar to exprgraph and return it in reformnode */
   SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, reformnode) );

   /* set value of auxvar and reformnode in debug solution */
#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(node));
      SCIPexprgraphSetVarNodeValue(*reformnode, SCIPexprgraphGetNodeVal(node));
   }
#endif

   SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for bivariate constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBivariate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create bivariate constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpBivariate, consEnfopsBivariate, consCheckBivariate, consLockBivariate,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveBivariate) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyBivariate, consCopyBivariate) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveBivariate) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteBivariate) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableBivariate) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableBivariate) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitBivariate) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreBivariate) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolBivariate) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeBivariate) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsBivariate) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsBivariate) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitBivariate) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreBivariate) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolBivariate) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpBivariate) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolBivariate, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintBivariate) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropBivariate, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpBivariate, consSepasolBivariate, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransBivariate) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxBivariate) );

   /* include the quadratic constraint upgrade in the quadratic constraint handler */
   SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, quadconsUpgdBivariate, QUADCONSUPGD_PRIORITY, FALSE, CONSHDLR_NAME) );

   /* include the quadratic constraint upgrade in the quadratic constraint handler */
   SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, NULL, exprgraphnodeReformBivariate, NONLINCONSUPGD_PRIORITY, FALSE, CONSHDLR_NAME) );

   /* add bivariate constraint handler parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/cutmaxrange",
         "maximal coef range of a cut (maximal coefficient divided by minimal coefficient) in order to be added to LP relaxation",
         &conshdlrdata->cutmaxrange, TRUE, 1e+7, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/linfeasshift",
         "whether to try to make solutions in check function feasible by shifting a linear variable (esp. useful if constraint was actually objective function)",
         &conshdlrdata->linfeasshift, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a single constraint within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 1, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/ninitlprefpoints",
         "number of reference points in each direction where to compute linear support for envelope in LP initialization",
         &conshdlrdata->ninitlprefpoints, FALSE, 3, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/enfocutsremovable",
         "are cuts added during enforcement removable from the LP in the same node?",
         &conshdlrdata->enfocutsremovable, TRUE, FALSE, NULL, NULL) );

   conshdlrdata->linvareventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(conshdlrdata->linvareventhdlr), CONSHDLR_NAME"_boundchange", "signals a bound tightening in a linear variable to a bivariate constraint",
         processLinearVarEvent, NULL) );
   assert(conshdlrdata->linvareventhdlr != NULL);

   conshdlrdata->nonlinvareventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(conshdlrdata->nonlinvareventhdlr), CONSHDLR_NAME"_boundchange2", "signals a bound change in a nonlinear variable to the bivariate constraint handler",
         processNonlinearVarEvent, (SCIP_EVENTHDLRDATA*)conshdlrdata) );
   assert(conshdlrdata->nonlinvareventhdlr != NULL);

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, CONSHDLR_NAME"_newsolution", "handles the event that a new primal solution has been found",
         processNewSolutionEvent, NULL) );

   /* create expression interpreter */
   SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &conshdlrdata->exprinterpreter) );

   /* create expression graph */
   SCIP_CALL( SCIPexprgraphCreate(SCIPblkmem(scip), &conshdlrdata->exprgraph, -1, -1,
         exprgraphVarAdded, exprgraphVarRemove, NULL, (void*)conshdlrdata) );
   conshdlrdata->isremovedfixings = TRUE;
   conshdlrdata->ispropagated = TRUE;

   conshdlrdata->scip = scip;

   return SCIP_OKAY;
}

/** creates and captures a bivariate constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPRTREE*        f,                  /**< expression tree specifying bivariate function f(x,y) */
   SCIP_BIVAR_CONVEXITY  convextype,         /**< kind of convexity of f(x,y) */
   SCIP_VAR*             z,                  /**< linear variable in constraint */
   SCIP_Real             zcoef,              /**< coefficient of linear variable */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
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
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(f != NULL);
   assert(!SCIPisInfinity(scip, REALABS(zcoef)));
   assert(modifiable == FALSE); /* we do not support column generation */

   /* find the bivariate constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("bivariate constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   BMSclearMemory(consdata);

   SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &consdata->f, f) );
   consdata->convextype = convextype;
   consdata->z          = z;
   consdata->zcoef      = zcoef;
   consdata->lhs        = lhs;
   consdata->rhs        = rhs;

   assert(SCIPexprtreeGetNVars(consdata->f) == 2);
   assert(SCIPexprtreeGetVars(consdata->f) != NULL);
   assert(SCIPexprtreeGetVars(consdata->f)[0] != NULL);
   assert(SCIPexprtreeGetVars(consdata->f)[1] != NULL);

   /* mark that variable events are not catched so far */
   consdata->eventfilterpos = -1;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures an absolute power constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsBivariate(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsBivariate() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPRTREE*        f,                  /**< expression tree specifying bivariate function f(x,y) */
   SCIP_BIVAR_CONVEXITY  convextype,         /**< kind of convexity of f(x,y) */
   SCIP_VAR*             z,                  /**< linear variable in constraint */
   SCIP_Real             zcoef,              /**< coefficient of linear variable */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsBivariate(scip, cons, name, f, convextype, z, zcoef, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** gets the linear variable of a bivariate constraint, or NULL if no such variable */
SCIP_VAR* SCIPgetLinearVarBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->z;
}

/** gets the coefficients of the linear variable of a bivariate constraint */
SCIP_Real SCIPgetLinearCoefBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->zcoef;
}

/** gets the expression tree of a bivariate constraint */
SCIP_EXPRTREE* SCIPgetExprtreeBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->f;
}

/** gets the left hand side of a bivariate constraint */
SCIP_Real SCIPgetLhsBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lhs;
}

/** gets the right hand side of a bivariate constraint */
SCIP_Real SCIPgetRhsBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhs;
}
