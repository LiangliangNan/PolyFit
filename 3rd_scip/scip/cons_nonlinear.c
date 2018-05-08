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

/**@file   cons_nonlinear.c
 * @brief  constraint handler for nonlinear constraints \f$\textrm{lhs} \leq \sum_{i=1}^n a_ix_i + \sum_{j=1}^m c_jf_j(x) \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 * @author Ingmar Vierhaus (consparse)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_nonlinear.h"
#define SCIP_PRIVATE_ROWPREP
#include "scip/cons_quadratic.h"  /* for SCIP_ROWPREP */
#include "scip/cons_linear.h"
#include "scip/heur_trysol.h"
#include "scip/heur_subnlp.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi_ipopt.h"
#include "scip/debug.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "nonlinear"
#define CONSHDLR_DESC          "constraint handler for nonlinear constraints"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -60 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000010 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_ALWAYS /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */

#define INTERVALINFTY             1E+43 /**< value for infinity in interval operations */
#define BOUNDTIGHTENING_MINSTRENGTH 0.05/**< minimal required bound tightening strength in expression graph domain tightening for propagating bound change */
#define INITLPMAXVARVAL          1000.0 /**< maximal absolute value of variable for still generating a linearization cut at that point in initlp */

/*
 * Data structures
 */

/** event data for linear variable bound change events */
struct LinVarEventData
{
   SCIP_CONSHDLRDATA*    conshdlrdata;       /**< the constraint handler data */
   SCIP_CONS*            cons;               /**< the constraint */
   int                   varidx;             /**< the index of the linear variable which bound change is catched */
   int                   filterpos;          /**< position of eventdata in SCIP's event filter */
};
typedef struct LinVarEventData LINVAREVENTDATA;

/** constraint data for nonlinear constraints */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of constraint */
   SCIP_Real             rhs;                /**< right hand side of constraint */

   int                   nlinvars;           /**< number of linear variables */
   int                   linvarssize;        /**< length of linear variable arrays */
   SCIP_VAR**            linvars;            /**< linear variables */
   SCIP_Real*            lincoefs;           /**< coefficients of linear variables */
   LINVAREVENTDATA**     lineventdata;       /**< eventdata for bound change of linear variable */

   int                   nexprtrees;         /**< number of expression trees */
   SCIP_Real*            nonlincoefs;        /**< coefficients of expression trees */
   SCIP_EXPRTREE**       exprtrees;          /**< nonlinear part of constraint */
   SCIP_EXPRCURV*        curvatures;         /**< curvature of each expression tree (taking nonlincoefs into account) */
   SCIP_EXPRGRAPHNODE*   exprgraphnode;      /**< node in expression graph corresponding to expression tree of this constraint */
   SCIP_EXPRCURV         curvature;          /**< curvature of complete nonlinear part, if checked */

   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */

   unsigned int          linvarssorted:1;    /**< are the linear variables already sorted? */
   unsigned int          linvarsmerged:1;    /**< are equal linear variables already merged? */

   unsigned int          iscurvchecked:1;    /**< is nonlinear function checked on convexity or concavity ? */
   unsigned int          isremovedfixingslin:1; /**< did we removed fixed/aggr/multiaggr variables in linear part? */
   unsigned int          ispresolved:1;      /**< did we checked for possibilities of upgrading or implicit integer variables? */
   unsigned int          forcebackprop:1;    /**< should we force to run the backward propagation on our subgraph in the exprgraph? */

   SCIP_Real             minlinactivity;     /**< sum of minimal activities of all linear terms with finite minimal activity */
   SCIP_Real             maxlinactivity;     /**< sum of maximal activities of all linear terms with finite maximal activity */
   int                   minlinactivityinf;  /**< number of linear terms with infinite minimal activity */
   int                   maxlinactivityinf;  /**< number of linear terms with infinity maximal activity */
   SCIP_Real             activity;           /**< activity of constraint function w.r.t. current solution */
   SCIP_Real             lhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */

   int                   linvar_maydecrease; /**< index of a variable in linvars that may be decreased without making any other constraint infeasible, or -1 if none */
   int                   linvar_mayincrease; /**< index of a variable in linvars that may be increased without making any other constraint infeasible, or -1 if none */

   SCIP_Real             lincoefsmin;        /**< maximal absolute value of coefficients in linear part, only available in solving stage */
   SCIP_Real             lincoefsmax;        /**< minimal absolute value of coefficients in linear part, only available in solving stage */
   unsigned int          ncuts;              /**< number of cuts created for this constraint so far */
};

/** nonlinear constraint update method */
struct SCIP_NlConsUpgrade
{
   SCIP_DECL_NONLINCONSUPGD((*nlconsupgd));  /**< method to call for upgrading nonlinear constraint */
   SCIP_DECL_EXPRGRAPHNODEREFORM((*nodereform));/**< method to call for reformulating an expression graph node */
   int                   priority;           /**< priority of upgrading method */
   SCIP_Bool             active;             /**< is upgrading enabled */
};
typedef struct SCIP_NlConsUpgrade SCIP_NLCONSUPGRADE;

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EXPRINT*         exprinterpreter;    /**< expression interpreter to compute gradients */

   SCIP_Real             cutmaxrange;        /**< maximal range (maximal coef / minimal coef) of a cut in order to be added to LP */
   SCIP_Bool             linfeasshift;       /**< whether to make solutions in check feasible if possible */
   SCIP_Bool             checkconvexexpensive;/**< whether to apply expensive curvature checking methods */
   SCIP_Bool             assumeconvex;       /**< whether functions in inequalities should be assumed to be convex */
   int                   maxproprounds;      /**< limit on number of propagation rounds for a single constraint within one round of SCIP propagation */
   SCIP_Bool             reformulate;        /**< whether to reformulate expression graph */
   int                   maxexpansionexponent;/**< maximal exponent where still expanding non-monomial polynomials in expression simplification */
   SCIP_Real             sepanlpmincont;     /**< minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation */
   SCIP_Bool             enfocutsremovable;  /**< are cuts added during enforcement removable from the LP in the same node? */

   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subNLP heuristic, if available */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the TRYSOL heuristic, if available */
   SCIP_EVENTHDLR*       linvareventhdlr;    /**< our handler for linear variable bound change events */
   SCIP_EVENTHDLR*       nonlinvareventhdlr; /**< our handler for nonlinear variable bound change events */
   int                   newsoleventfilterpos;/**< filter position of new solution event handler, if catched */

   SCIP_NLCONSUPGRADE**  nlconsupgrades;     /**< nonlinear constraint upgrade methods for specializing nonlinear constraints */
   int                   nlconsupgradessize; /**< size of nlconsupgrade array */
   int                   nnlconsupgrades;    /**< number of nonlinear constraint upgrade methods */

   SCIP_EXPRGRAPH*       exprgraph;          /**< expression graph */
   SCIP*                 scip;               /**< SCIP pointer for use in expression graph callbacks */
   unsigned int          isremovedfixings:1; /**< have fixed variables been removed in the expression graph? */
   unsigned int          ispropagated:1;     /**< have current bounds of linear variables in constraints and variables in expression graph been propagated? */
   unsigned int          isreformulated:1;   /**< has expression graph been reformulated? */
   unsigned int          sepanlp:1;          /**< has a linearization in the NLP relaxation been added? */
   int                   naddedreformconss;  /**< number of constraints added via reformulation */
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

/* catches variable bound change events on a linear variable in a nonlinear constraint */
static
SCIP_RETCODE catchLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   linvarpos           /**< position of variable in linear variables array */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   LINVAREVENTDATA* eventdata;
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

   assert(linvarpos >= 0);
   assert(linvarpos < consdata->nlinvars);

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventdata) );

   eventtype = SCIP_EVENTTYPE_VARFIXED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
   }

   eventdata->conshdlrdata = conshdlrdata;
   eventdata->cons = cons;
   eventdata->varidx = linvarpos;
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->linvars[linvarpos], eventtype, conshdlrdata->linvareventhdlr, (SCIP_EVENTDATA*)eventdata, &eventdata->filterpos) );

   /* ensure lineventdata array is existing */
   if( consdata->lineventdata == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize) );
   }

   consdata->lineventdata[linvarpos] = eventdata;

   /* since bound changes were not catched before, a possibly stored linear activity may have become outdated, so set to invalid */
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;

   /* mark constraint for propagation */
   SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );

   return SCIP_OKAY;
}

/* drops variable bound change events on a linear variable in a nonlinear constraint */
static
SCIP_RETCODE dropLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   linvarpos           /**< position of variable in linear variables array */
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

   assert(linvarpos >= 0);
   assert(linvarpos < consdata->nlinvars);
   assert(consdata->lineventdata != NULL);
   assert(consdata->lineventdata[linvarpos] != NULL);
   assert(consdata->lineventdata[linvarpos]->cons == cons);
   assert(consdata->lineventdata[linvarpos]->varidx == linvarpos);
   assert(consdata->lineventdata[linvarpos]->filterpos >= 0);

   eventtype = SCIP_EVENTTYPE_VARFIXED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
   }

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->linvars[linvarpos], eventtype, conshdlrdata->linvareventhdlr, (SCIP_EVENTDATA*)consdata->lineventdata[linvarpos], consdata->lineventdata[linvarpos]->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->lineventdata[linvarpos]);  /*lint !e866*/

   return SCIP_OKAY;
}

/** locks a linear variable in a constraint */
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

/** unlocks a linear variable in a constraint */
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

/** computes the minimal and maximal activity for the linear part in a constraint data
 *  only sums up terms that contribute finite values
 *  gives the number of terms that contribute infinite values
 *  only computes those activities where the corresponding side of the constraint is finite
 */
static
void consdataUpdateLinearActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{  /*lint --e{666}*/
   SCIP_ROUNDMODE prevroundmode;
   int       i;
   SCIP_Real bnd;

   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->minlinactivity != SCIP_INVALID && consdata->maxlinactivity != SCIP_INVALID )  /*lint !e777*/
   {
      /* activities should be uptodate */
      assert(consdata->minlinactivityinf >= 0);
      assert(consdata->maxlinactivityinf >= 0);
      return;
   }

   consdata->minlinactivityinf = 0;
   consdata->maxlinactivityinf = 0;

   /* if lhs is -infinite, then we do not compute a maximal activity, so we set it to  infinity
    * if rhs is  infinite, then we do not compute a minimal activity, so we set it to -infinity
    */
   consdata->minlinactivity = SCIPisInfinity(scip,  consdata->rhs) ? -INTERVALINFTY : 0.0;
   consdata->maxlinactivity = SCIPisInfinity(scip, -consdata->lhs) ?  INTERVALINFTY : 0.0;

   if( consdata->nlinvars == 0 )
      return;

   /* if the activities computed here should be still uptodate after boundchanges,
    * variable events need to be catched */
   assert(consdata->lineventdata != NULL);

   prevroundmode = SCIPintervalGetRoundingMode();

   if( !SCIPisInfinity(scip,  consdata->rhs) )
   {
      /* compute minimal activity only if there is a finite right hand side */
      SCIPintervalSetRoundingModeDownwards();

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         assert(SCIPvarGetLbLocal(consdata->linvars[i]) <= SCIPvarGetUbLocal(consdata->linvars[i]));
         assert(consdata->lineventdata[i] != NULL);
         if( consdata->lincoefs[i] >= 0.0 )
         {
            bnd = SCIPvarGetLbLocal(consdata->linvars[i]);
            if( SCIPisInfinity(scip, -bnd) )
            {
               ++consdata->minlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, bnd)); /* do not like variables that are fixed at +infinity */
         }
         else
         {
            bnd = SCIPvarGetUbLocal(consdata->linvars[i]);
            if( SCIPisInfinity(scip,  bnd) )
            {
               ++consdata->minlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, -bnd)); /* do not like variables that are fixed at -infinity */
         }
         consdata->minlinactivity += consdata->lincoefs[i] * bnd;
      }
   }

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* compute maximal activity only if there is a finite left hand side */
      SCIPintervalSetRoundingModeUpwards();

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         assert(consdata->lineventdata[i] != NULL);
         assert(SCIPvarGetLbLocal(consdata->linvars[i]) <= SCIPvarGetUbLocal(consdata->linvars[i]));
         if( consdata->lincoefs[i] >= 0.0 )
         {
            bnd = SCIPvarGetUbLocal(consdata->linvars[i]);
            if( SCIPisInfinity(scip,  bnd) )
            {
               ++consdata->maxlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, -bnd)); /* do not like variables that are fixed at -infinity */
         }
         else
         {
            bnd = SCIPvarGetLbLocal(consdata->linvars[i]);
            if( SCIPisInfinity(scip, -bnd) )
            {
               ++consdata->maxlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, bnd)); /* do not like variables that are fixed at +infinity */
         }
         consdata->maxlinactivity += consdata->lincoefs[i] * bnd;
      }
   }
   assert(consdata->minlinactivity <= consdata->maxlinactivity || consdata->minlinactivityinf > 0 || consdata->maxlinactivityinf > 0);

   SCIPintervalSetRoundingMode(prevroundmode);
}

/** update the linear activities after a change in the lower bound of a variable */
static
void consdataUpdateLinearActivityLbChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             oldbnd,             /**< previous lower bound of variable */
   SCIP_Real             newbnd              /**< new lower bound of variable */
   )
{
   SCIP_ROUNDMODE prevroundmode;

   assert(scip != NULL);
   assert(consdata != NULL);
   /* we can't deal with lower bounds at infinity */
   assert(!SCIPisInfinity(scip, oldbnd));
   assert(!SCIPisInfinity(scip, newbnd));

   /* assume lhs <= a*x + y <= rhs, then the following boundchanges can be deduced:
    * a > 0:  y <= rhs - a*lb(x),  y >= lhs - a*ub(x)
    * a < 0:  y <= rhs - a*ub(x),  y >= lhs - a*lb(x)
    */

   if( coef > 0.0 )
   {
      /* we should only be called if rhs is finite */
      assert(!SCIPisInfinity(scip, consdata->rhs));

      /* we have no min activities computed so far, so cannot update */
      if( consdata->minlinactivity == SCIP_INVALID )  /*lint !e777*/
         return;

      assert(consdata->minlinactivity > -INTERVALINFTY);

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      /* update min activity */
      if( SCIPisInfinity(scip, -oldbnd) )
      {
         --consdata->minlinactivityinf;
         assert(consdata->minlinactivityinf >= 0);
      }
      else
      {
         consdata->minlinactivity += SCIPintervalNegateReal(coef) * oldbnd;
      }

      if( SCIPisInfinity(scip, -newbnd) )
      {
         ++consdata->minlinactivityinf;
      }
      else
      {
         consdata->minlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
   else
   {
      /* we should only be called if lhs is finite */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      /* we have no max activities computed so far, so cannot update */
      if( consdata->maxlinactivity == SCIP_INVALID )  /*lint !e777*/
         return;

      assert(consdata->maxlinactivity < INTERVALINFTY);

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();

      /* update max activity */
      if( SCIPisInfinity(scip, -oldbnd) )
      {
         --consdata->maxlinactivityinf;
         assert(consdata->maxlinactivityinf >= 0);
      }
      else
      {
         consdata->maxlinactivity += SCIPintervalNegateReal(coef) * oldbnd;
      }

      if( SCIPisInfinity(scip, -newbnd) )
      {
         ++consdata->maxlinactivityinf;
      }
      else
      {
         consdata->maxlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
}

/** update the linear activities after a change in the upper bound of a variable */
static
void consdataUpdateLinearActivityUbChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             oldbnd,             /**< previous lower bound of variable */
   SCIP_Real             newbnd              /**< new lower bound of variable */
   )
{
   SCIP_ROUNDMODE prevroundmode;

   assert(scip != NULL);
   assert(consdata != NULL);
   /* we can't deal with upper bounds at -infinity */
   assert(!SCIPisInfinity(scip, -oldbnd));
   assert(!SCIPisInfinity(scip, -newbnd));

   /* assume lhs <= a*x + y <= rhs, then the following boundchanges can be deduced:
    * a > 0:  y <= rhs - a*lb(x),  y >= lhs - a*ub(x)
    * a < 0:  y <= rhs - a*ub(x),  y >= lhs - a*lb(x)
    */
   if( coef > 0.0 )
   {
      /* we should only be called if lhs is finite */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      /* we have no max activities computed so far, so cannot update */
      if( consdata->maxlinactivity == SCIP_INVALID )  /*lint !e777*/
         return;

      assert(consdata->maxlinactivity < INTERVALINFTY);

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();

      /* update max activity */
      if( SCIPisInfinity(scip, oldbnd) )
      {
         --consdata->maxlinactivityinf;
         assert(consdata->maxlinactivityinf >= 0);
      }
      else
      {
         consdata->maxlinactivity += SCIPintervalNegateReal(coef) * oldbnd;
      }

      if( SCIPisInfinity(scip, newbnd) )
      {
         ++consdata->maxlinactivityinf;
      }
      else
      {
         consdata->maxlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
   else
   {
      /* we should only be called if rhs is finite */
      assert(!SCIPisInfinity(scip, consdata->rhs));

      /* we have no min activities computed so far, so cannot update */
      if( consdata->minlinactivity == SCIP_INVALID )  /*lint !e777*/
         return;

      assert(consdata->minlinactivity > -INTERVALINFTY);

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      /* update min activity */
      if( SCIPisInfinity(scip, oldbnd) )
      {
         --consdata->minlinactivityinf;
         assert(consdata->minlinactivityinf >= 0);
      }
      else
      {
         consdata->minlinactivity += SCIPintervalNegateReal(coef) * oldbnd;
      }

      if( SCIPisInfinity(scip, newbnd) )
      {
         ++consdata->minlinactivityinf;
      }
      else
      {
         consdata->minlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processLinearVarEvent)
{
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   int varidx;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   cons = ((LINVAREVENTDATA*)eventdata)->cons;
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   varidx = ((LINVAREVENTDATA*)eventdata)->varidx;
   assert(varidx >= 0);
   assert(varidx < consdata->nlinvars);

   eventtype = SCIPeventGetType(event);

   if( eventtype & SCIP_EVENTTYPE_VARFIXED )
   {
      consdata->isremovedfixingslin = FALSE;
   }

   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      /* update activity bounds for linear terms */
      if( eventtype & SCIP_EVENTTYPE_LBCHANGED )
         consdataUpdateLinearActivityLbChange(scip, consdata, consdata->lincoefs[varidx], SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));
      else
         consdataUpdateLinearActivityUbChange(scip, consdata, consdata->lincoefs[varidx], SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));

      if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
      {
         assert(((LINVAREVENTDATA*)eventdata)->conshdlrdata != NULL);
         ((LINVAREVENTDATA*)eventdata)->conshdlrdata->ispropagated = FALSE;

         /* mark constraint for propagation */
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
      }
   }

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
      SCIP_Real newbd;

      SCIPdebugMsg(scip, "changed %s bound on expression graph variable <%s> from %g to %g\n",
         (eventtype & SCIP_EVENTTYPE_LBCHANGED) ? "lower" : "upper",
         SCIPvarGetName(SCIPeventGetVar(event)), SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));

      if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
         conshdlrdata->ispropagated = FALSE;
      /* @todo a global bound tightening may yield in convex/concave curvatures, maybe the iscurvcheck flag should be reset? */

      /* update variable bound in expression graph, with epsilon added */
      if( eventtype & SCIP_EVENTTYPE_LBCHANGED )
      {
         newbd = -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -SCIPeventGetNewbound(event));  /*lint !e666*/
         /* if newbd in [0,eps], then relax to 0.0, otherwise relax by -epsilon
          * this is to ensure that an original positive variable does not get negative by this, which may have an adverse effect on convexity recoginition, for example */
         if( newbd >= 0.0 && newbd <= SCIPepsilon(scip) )
            newbd = 0.0;
         else
            newbd -= SCIPepsilon(scip);
         SCIPexprgraphSetVarNodeLb(conshdlrdata->exprgraph, (SCIP_EXPRGRAPHNODE*)eventdata, newbd);
      }
      else
      {
         newbd = +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  SCIPeventGetNewbound(event));  /*lint !e666*/
         /* if newbd in [-eps,0], then relax to 0.0, otherwise relax by +epsilon */
         if( newbd >= -SCIPepsilon(scip) && newbd <= 0.0 )
            newbd = 0.0;
         else
            newbd += SCIPepsilon(scip);
         SCIPexprgraphSetVarNodeUb(conshdlrdata->exprgraph, (SCIP_EXPRGRAPHNODE*)eventdata, newbd);
      }
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

   SCIP_CALL( SCIPcaptureVar(conshdlrdata->scip, var_) );
   SCIPdebugMessage("capture variable <%s>\n", SCIPvarGetName(var_));

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

   SCIPdebugMessage("release variable <%s>\n", SCIPvarGetName(var_));
   SCIP_CALL( SCIPreleaseVar(conshdlrdata->scip, &var_) );

   return SCIP_OKAY;
}

/* adds expression trees to constraint */
static
SCIP_RETCODE consdataAddExprtrees(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< nonlinear constraint data */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            coefs,              /**< coefficients of expression trees, or NULL if all 1.0 */
   SCIP_Bool             copytrees           /**< whether trees should be copied or ownership should be assumed */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->exprtrees != NULL || consdata->nexprtrees == 0);

   if( nexprtrees == 0 )
      return SCIP_OKAY;

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispresolved  = FALSE;
   consdata->curvature = SCIP_EXPRCURV_UNKNOWN;
   consdata->iscurvchecked = FALSE;

   if( consdata->nexprtrees == 0 )
   {
      assert(consdata->exprtrees   == NULL);
      assert(consdata->nonlincoefs == NULL);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->exprtrees,   nexprtrees) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->nonlincoefs, nexprtrees) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->curvatures,  nexprtrees) );
   }
   else
   {
      assert(consdata->exprtrees   != NULL);
      assert(consdata->nonlincoefs != NULL);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->exprtrees,   consdata->nexprtrees, consdata->nexprtrees + nexprtrees) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nonlincoefs, consdata->nexprtrees, consdata->nexprtrees + nexprtrees) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->curvatures,  consdata->nexprtrees, consdata->nexprtrees + nexprtrees) );
   }

   for( i = 0; i < nexprtrees; ++i )
   {
      assert(exprtrees[i] != NULL);
      /* the expression tree need to have SCIP_VAR*'s stored */
      assert(SCIPexprtreeGetNVars(exprtrees[i]) == 0 || SCIPexprtreeGetVars(exprtrees[i]) != NULL);

      if( copytrees )
      {
         SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &consdata->exprtrees[consdata->nexprtrees + i], exprtrees[i]) );
      }
      else
      {
         consdata->exprtrees[consdata->nexprtrees + i] = exprtrees[i];
      }

      consdata->nonlincoefs[consdata->nexprtrees + i] = (coefs != NULL ? coefs[i] : 1.0);
      consdata->curvatures[consdata->nexprtrees + i] = SCIP_EXPRCURV_UNKNOWN;
   }
   consdata->nexprtrees += nexprtrees;

   return SCIP_OKAY;
}

/* sets expression trees of constraints, clears previously ones */
static
SCIP_RETCODE consdataSetExprtrees(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< nonlinear constraint data */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            coefs,              /**< coefficients of expression trees, or NULL if all 1.0 */
   SCIP_Bool             copytrees           /**< whether trees should be copied or ownership should be assumed */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->exprtrees != NULL || consdata->nexprtrees == 0);

   /* clear existing expression trees */
   if( consdata->nexprtrees > 0 )
   {
      for( i = 0; i < consdata->nexprtrees; ++i )
      {
         assert(consdata->exprtrees[i] != NULL);
         SCIP_CALL( SCIPexprtreeFree(&consdata->exprtrees[i]) );
      }

      /* invalidate activity information */
      consdata->activity = SCIP_INVALID;

      /* invalidate nonlinear row */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }

      consdata->ispresolved  = FALSE;
      consdata->curvature = SCIP_EXPRCURV_LINEAR;
      consdata->iscurvchecked = TRUE;

      SCIPfreeBlockMemoryArray(scip, &consdata->exprtrees,   consdata->nexprtrees);
      SCIPfreeBlockMemoryArray(scip, &consdata->nonlincoefs, consdata->nexprtrees);
      SCIPfreeBlockMemoryArray(scip, &consdata->curvatures,  consdata->nexprtrees);
      consdata->nexprtrees = 0;
   }

   SCIP_CALL( consdataAddExprtrees(scip, consdata, nexprtrees, exprtrees, coefs, copytrees) );

   return SCIP_OKAY;
}

/** ensures, that linear vars and coefs arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureLinearVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< nonlinear constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nlinvars <= consdata->linvarssize);

   if( num > consdata->linvarssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->linvars,  consdata->linvarssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->lincoefs, consdata->linvarssize, newsize) );
      if( consdata->lineventdata != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize, newsize) );
      }
      consdata->linvarssize = newsize;
   }
   assert(num <= consdata->linvarssize);

   return SCIP_OKAY;
}

/** creates constraint data structure */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< a buffer to store pointer to new constraint data */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            lincoefs,           /**< array of coefficients of linear variables */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            nonlincoefs,        /**< coefficients of expression trees */
   SCIP_Bool             capturevars         /**< whether we should capture variables */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);

   assert(nlinvars == 0 || linvars  != NULL);
   assert(nlinvars == 0 || lincoefs != NULL);
   assert(nexprtrees == 0 || exprtrees != NULL);
   assert(nexprtrees == 0 || nonlincoefs != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   BMSclearMemory(*consdata);

   (*consdata)->minlinactivity = SCIP_INVALID;
   (*consdata)->maxlinactivity = SCIP_INVALID;
   (*consdata)->minlinactivityinf = -1;
   (*consdata)->maxlinactivityinf = -1;

   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;

   if( nlinvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->linvars,  linvars,  nlinvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->lincoefs, lincoefs, nlinvars) );
      (*consdata)->nlinvars = nlinvars;
      (*consdata)->linvarssize = nlinvars;

      if( capturevars )
         for( i = 0; i < nlinvars; ++i )
         {
            SCIP_CALL( SCIPcaptureVar(scip, linvars[i]) );
         }
   }
   else
   {
      (*consdata)->linvarssorted = TRUE;
      (*consdata)->linvarsmerged = TRUE;
   }

   SCIP_CALL( consdataSetExprtrees(scip, *consdata, nexprtrees, exprtrees, nonlincoefs, TRUE) );

   (*consdata)->linvar_maydecrease = -1;
   (*consdata)->linvar_mayincrease = -1;

   (*consdata)->activity = SCIP_INVALID;
   (*consdata)->lhsviol  = SCIPisInfinity(scip, -lhs) ? 0.0 : SCIP_INVALID;
   (*consdata)->rhsviol  = SCIPisInfinity(scip,  rhs) ? 0.0 : SCIP_INVALID;

   return SCIP_OKAY;
}

/** creates empty constraint data structure */
static
SCIP_RETCODE consdataCreateEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< a buffer to store pointer to new constraint data */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   BMSclearMemory(*consdata);

   (*consdata)->lhs = -SCIPinfinity(scip);
   (*consdata)->rhs =  SCIPinfinity(scip);

   (*consdata)->linvarssorted    = TRUE;
   (*consdata)->linvarsmerged    = TRUE;

   (*consdata)->isremovedfixingslin = TRUE;

   (*consdata)->linvar_maydecrease = -1;
   (*consdata)->linvar_mayincrease = -1;

   (*consdata)->minlinactivity = SCIP_INVALID;
   (*consdata)->maxlinactivity = SCIP_INVALID;
   (*consdata)->minlinactivityinf = -1;
   (*consdata)->maxlinactivityinf = -1;

   (*consdata)->curvature = SCIP_EXPRCURV_LINEAR;
   (*consdata)->iscurvchecked = TRUE;

   (*consdata)->ncuts = 0;

   return SCIP_OKAY;
}

/** frees constraint data structure */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data to free */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release linear variables and free linear part */
   if( (*consdata)->linvarssize > 0 )
   {
      for( i = 0; i < (*consdata)->nlinvars; ++i )
      {
         assert((*consdata)->lineventdata == NULL || (*consdata)->lineventdata[i] == NULL); /* variable events should have been dropped earlier */
         SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->linvars[i]) );
      }
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->linvars,  (*consdata)->linvarssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->lincoefs, (*consdata)->linvarssize);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lineventdata, (*consdata)->linvarssize);
   }
   assert((*consdata)->linvars == NULL);
   assert((*consdata)->lincoefs == NULL);
   assert((*consdata)->lineventdata == NULL);

   if( (*consdata)->nexprtrees > 0 )
   {
      assert((*consdata)->exprtrees   != NULL);
      assert((*consdata)->nonlincoefs != NULL);
      assert((*consdata)->curvatures  != NULL);
      for( i = 0; i < (*consdata)->nexprtrees; ++i )
      {
         assert((*consdata)->exprtrees[i] != NULL);
         SCIP_CALL( SCIPexprtreeFree(&(*consdata)->exprtrees[i]) );
         assert((*consdata)->exprtrees[i] == NULL);
      }
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->exprtrees,   (*consdata)->nexprtrees);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->nonlincoefs, (*consdata)->nexprtrees);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->curvatures,  (*consdata)->nexprtrees);
   }
   assert((*consdata)->exprtrees   == NULL);
   assert((*consdata)->nonlincoefs == NULL);
   assert((*consdata)->curvatures  == NULL);

   /* free nonlinear row representation */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeBlockMemory(scip, consdata);
   *consdata = NULL;

   return SCIP_OKAY;
}

/** sorts linear part of constraint data */
static
void consdataSortLinearVars(
   SCIP_CONSDATA*        consdata            /**< nonlinear constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->linvarssorted )
      return;

   if( consdata->nlinvars <= 1 )
   {
      consdata->linvarssorted = TRUE;
      return;
   }

   if( consdata->lineventdata == NULL )
   {
      SCIPsortPtrReal((void**)consdata->linvars, consdata->lincoefs, SCIPvarComp, consdata->nlinvars);
   }
   else
   {
      int i;

      SCIPsortPtrPtrReal((void**)consdata->linvars, (void**)consdata->lineventdata, consdata->lincoefs, SCIPvarComp, consdata->nlinvars);

      /* update variable indices in event data */
      for( i = 0; i < consdata->nlinvars; ++i )
         if( consdata->lineventdata[i] != NULL )
            consdata->lineventdata[i]->varidx = i;
   }

   consdata->linvarssorted = TRUE;
}

/* this function is currently not needed, but also to nice to be deleted, so it is only deactivated */
#ifdef SCIP_DISABLED_CODE
/** returns the position of variable in the linear coefficients array of a constraint, or -1 if not found */
static
int consdataFindLinearVar(
   SCIP_CONSDATA*        consdata,           /**< nonlinear constraint data */
   SCIP_VAR*             var                 /**< variable to search for */
   )
{
   int pos;

   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->nlinvars == 0 )
      return -1;

   consdataSortLinearVars(consdata);

   if( !SCIPsortedvecFindPtr((void**)consdata->linvars, SCIPvarComp, (void*)var, consdata->nlinvars, &pos) )
      pos = -1;

   return pos;
}
#endif

/** moves a linear variable from one position to another */
static
void consdataMoveLinearVar(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   oldpos,             /**< position of variable that shall be moved */
   int                   newpos              /**< new position of variable */
   )
{
   assert(consdata != NULL);
   assert(oldpos >= 0);
   assert(oldpos < consdata->nlinvars);
   assert(newpos >= 0);
   assert(newpos < consdata->linvarssize);

   if( newpos == oldpos )
      return;

   consdata->linvars [newpos] = consdata->linvars [oldpos];
   consdata->lincoefs[newpos] = consdata->lincoefs[oldpos];

   if( consdata->lineventdata != NULL )
   {
      assert(newpos >= consdata->nlinvars || consdata->lineventdata[newpos] == NULL);

      consdata->lineventdata[newpos] = consdata->lineventdata[oldpos];
      consdata->lineventdata[newpos]->varidx = newpos;

      consdata->lineventdata[oldpos] = NULL;
   }

   consdata->linvarssorted = FALSE;
}

/** adds linear coefficient in nonlinear constraint */
static
SCIP_RETCODE addLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             coef                /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);

   /* ignore coefficient if it is nearly zero */
   if( SCIPisZero(scip, coef) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars+1) );
   consdata->linvars [consdata->nlinvars] = var;
   consdata->lincoefs[consdata->nlinvars] = coef;

   ++consdata->nlinvars;

   /* catch variable events */
   if( SCIPconsIsEnabled(cons) )
   {
      /* catch bound change events of variable */
      SCIP_CALL( catchLinearVarEvents(scip, cons, consdata->nlinvars-1) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockLinearVariable(scip, cons, var, coef) );

   /* capture new variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   consdata->ispresolved = FALSE;
   consdata->isremovedfixingslin = consdata->isremovedfixingslin && SCIPvarIsActive(var);
   if( consdata->nlinvars == 1 )
      consdata->linvarssorted = TRUE;
   else
      consdata->linvarssorted = consdata->linvarssorted &&
         (SCIPvarCompare(consdata->linvars[consdata->nlinvars-2], consdata->linvars[consdata->nlinvars-1]) == -1);
   /* always set too FALSE since the new linear variable should be checked if already existing as quad var term */
   consdata->linvarsmerged = FALSE;

   return SCIP_OKAY;
}

/** deletes linear coefficient at given position from nonlinear constraint data */
static
SCIP_RETCODE delLinearCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nlinvars);

   var  = consdata->linvars[pos];
   coef = consdata->lincoefs[pos];
   assert(var != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockLinearVariable(scip, cons, var, coef) );

   /* if constraint is enabled, drop the events on the variable */
   if( SCIPconsIsEnabled(cons) )
   {
      /* drop bound change events of variable */
      SCIP_CALL( dropLinearVarEvents(scip, cons, pos) );
   }

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->linvars[pos]) );

   /* move the last variable to the free slot */
   consdataMoveLinearVar(consdata, consdata->nlinvars-1, pos);

   --consdata->nlinvars;

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispresolved  = FALSE;

   return SCIP_OKAY;
}

/** changes linear coefficient value at given position of nonlinear constraint */
static
SCIP_RETCODE chgLinearCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   int                   pos,                /**< position of linear coefficient to change */
   SCIP_Real             newcoef             /**< new value of linear coefficient */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisZero(scip, newcoef));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos);
   assert(pos < consdata->nlinvars);
   assert(!SCIPisZero(scip, newcoef));

   var = consdata->linvars[pos];
   coef = consdata->lincoefs[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* if necessary, remove the rounding locks and event catching of the variable */
   if( newcoef * coef < 0.0 )
   {
      if( SCIPconsIsLocked(cons) )
      {
         assert(SCIPconsIsTransformed(cons));

         /* remove rounding locks for variable with old coefficient */
         SCIP_CALL( unlockLinearVariable(scip, cons, var, coef) );
      }

      if( SCIPconsIsEnabled(cons) )
      {
         /* drop bound change events of variable */
         SCIP_CALL( dropLinearVarEvents(scip, cons, pos) );
      }
   }

   /* change the coefficient */
   consdata->lincoefs[pos] = newcoef;

   /* if necessary, install the rounding locks and event catching of the variable again */
   if( newcoef * coef < 0.0 )
   {
      if( SCIPconsIsLocked(cons) )
      {
         /* install rounding locks for variable with new coefficient */
         SCIP_CALL( lockLinearVariable(scip, cons, var, newcoef) );
      }

      if( SCIPconsIsEnabled(cons) )
      {
         /* catch bound change events of variable */
         SCIP_CALL( catchLinearVarEvents(scip, cons, pos) );
      }
   }

   consdata->ispresolved = FALSE;

   return SCIP_OKAY;
}


/* merges entries with same linear variable into one entry and cleans up entries with coefficient 0.0 */
static
SCIP_RETCODE mergeAndCleanLinearVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< nonlinear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real newcoef;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( consdata->linvarsmerged )
      return SCIP_OKAY;

   if( consdata->nlinvars == 0 )
   {
      consdata->linvarsmerged = TRUE;
      return SCIP_OKAY;
   }

   i = 0;
   while( i < consdata->nlinvars )
   {
      /* make sure linear variables are sorted (do this in every round, since we may move variables around) */
      consdataSortLinearVars(consdata);

      /* sum up coefficients that correspond to variable i */
      newcoef = consdata->lincoefs[i];
      for( j = i+1; j < consdata->nlinvars && consdata->linvars[i] == consdata->linvars[j]; ++j )
         newcoef += consdata->lincoefs[j];
      /* delete the additional variables in backward order */
      for( j = j-1; j > i; --j )
      {
         SCIP_CALL( delLinearCoefPos(scip, cons, j) );
      }

      /* delete also entry at position i, if it became zero (or was zero before) */
      if( SCIPisZero(scip, newcoef) )
      {
         SCIP_CALL( delLinearCoefPos(scip, cons, i) );
      }
      else
      {
         SCIP_CALL( chgLinearCoefPos(scip, cons, i, newcoef) );
         ++i;
      }
   }

   consdata->linvarsmerged = TRUE;

   return SCIP_OKAY;
}

/** removes fixes (or aggregated) linear variables from a nonlinear constraint */
static
SCIP_RETCODE removeFixedLinearVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< nonlinearconstraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real coef;
   SCIP_Real offset;
   SCIP_VAR* var;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( !consdata->isremovedfixingslin )
   {
      i = 0;
      while( i < consdata->nlinvars )
      {
         var = consdata->linvars[i];

         if( SCIPvarIsActive(var) )
         {
            ++i;
            continue;
         }

         coef = consdata->lincoefs[i];
         offset = 0.0;

         SCIP_CALL( SCIPgetProbvarSum(scip, &var, &coef, &offset) );

         SCIPdebugMsg(scip, "  linear term %g*<%s> is replaced by %g * <%s> + %g\n", consdata->lincoefs[i], SCIPvarGetName(consdata->linvars[i]), coef, SCIPvarGetName(var), offset);

         /* delete previous variable (this will move another variable to position i) */
         SCIP_CALL( delLinearCoefPos(scip, cons, i) );

         /* put constant part into bounds */
         if( offset != 0.0 )
         {
            if( !SCIPisInfinity(scip, -consdata->lhs) )
            {
               consdata->lhs -= offset;
               assert(!SCIPisInfinity(scip, REALABS(consdata->lhs)));
            }
            if( !SCIPisInfinity(scip,  consdata->rhs) )
            {
               consdata->rhs -= offset;
               assert(!SCIPisInfinity(scip, REALABS(consdata->rhs)));
            }
         }

         /* nothing left to do if variable had been fixed */
         if( coef == 0.0 )
            continue;

         /* if GetProbvar gave a linear variable, just add it
          * if it's a multilinear variable, add it's disaggregated variables */
         if( SCIPvarIsActive(var) )
         {
            SCIP_CALL( addLinearCoef(scip, cons, var, coef) );
         }
         else
         {
            int        naggrs;
            SCIP_VAR** aggrvars;
            SCIP_Real* aggrscalars;
            SCIP_Real  aggrconstant;

            assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

            naggrs = SCIPvarGetMultaggrNVars(var);
            aggrvars = SCIPvarGetMultaggrVars(var);
            aggrscalars = SCIPvarGetMultaggrScalars(var);
            aggrconstant = SCIPvarGetMultaggrConstant(var);

            SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars + naggrs) );

            for( j = 0; j < naggrs; ++j )
            {
               SCIP_CALL( addLinearCoef(scip, cons, aggrvars[j], coef * aggrscalars[j]) );
            }

            if( aggrconstant != 0.0 )
            {
               if( !SCIPisInfinity(scip, -consdata->lhs) )
               {
                  consdata->lhs -= coef * aggrconstant;
                  assert(!SCIPisInfinity(scip, REALABS(consdata->lhs)));
               }
               if( !SCIPisInfinity(scip,  consdata->rhs) )
               {
                  consdata->rhs -= coef * aggrconstant;
                  assert(!SCIPisInfinity(scip, REALABS(consdata->rhs)));
               }
            }
         }
      }

      SCIP_CALL( mergeAndCleanLinearVars(scip, cons) );

      consdata->isremovedfixingslin = TRUE;
   }

   SCIPdebugMsg(scip, "removed fixations of linear variables from <%s>\n  -> ", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

#ifndef NDEBUG
   for( i = 0; i < consdata->nlinvars; ++i )
      assert(SCIPvarIsActive(consdata->linvars[i]));
#endif

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
      var = (SCIP_VAR*)SCIPexprgraphGetVars(conshdlrdata->exprgraph)[i];
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

/** moves constant and linear part from expression graph node into constraint sides and linear part
 * frees expression graph node if expression is constant or linear */
static
SCIP_RETCODE splitOffLinearPart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_Bool*            infeasible          /**< pointer to store whether the problem is infeasible or not */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   SCIP_Real constant;
   int linvarssize;
   int nlinvars;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *infeasible = FALSE;

   if( consdata->exprgraphnode == NULL )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   if( SCIPexprgraphGetNodeOperator(consdata->exprgraphnode) == SCIP_EXPR_SUM && SCIPexprgraphGetNodeNChildren(consdata->exprgraphnode) == 1 )
   {
      /* it can happen that simplifies leaves a node that is an identity w.r.t only child, if that node is the root of an expression, used by a constraint
       * replace exprgraphnode by its child then
       */
      SCIP_EXPRGRAPHNODE* child;

      child = SCIPexprgraphGetNodeChildren(consdata->exprgraphnode)[0];
      assert(child != NULL);

      SCIPexprgraphCaptureNode(child);
      SCIP_CALL( SCIPexprgraphReleaseNode(conshdlrdata->exprgraph, &consdata->exprgraphnode) );
      consdata->exprgraphnode = child;
   }

   /* number of children of expression graph node is a good upper estimate on number of linear variables */
   linvarssize = MAX(SCIPexprgraphGetNodeNChildren(consdata->exprgraphnode), 1);  /*lint !e666*/
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars,  linvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, linvarssize) );

   /* get linear and constant part from expression graph node
    * releases expression graph node if not uses otherwise */
   SCIP_CALL( SCIPexprgraphNodeSplitOffLinear(conshdlrdata->exprgraph, &consdata->exprgraphnode, linvarssize, &nlinvars, (void**)linvars, lincoefs, &constant) );

   if( SCIPisInfinity(scip, constant) )
   {
      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         /* setting constraint lhs to -infinity; this may change linear variable locks and events */
         for( i = 0; i < consdata->nlinvars; ++i )
         {
            if( SCIPconsIsLocked(cons) )
            {
               SCIP_CALL( unlockLinearVariable(scip, cons, consdata->linvars[i], consdata->lincoefs[i]) );
            }
            if( SCIPconsIsEnabled(cons) )
            {
               SCIP_CALL( dropLinearVarEvents(scip, cons, i) );
            }
         }

         consdata->lhs = -SCIPinfinity(scip);

         for( i = 0; i < consdata->nlinvars; ++i )
         {
            if( SCIPconsIsEnabled(cons) )
            {
               SCIP_CALL( catchLinearVarEvents(scip, cons, i) );
            }
            if( SCIPconsIsLocked(cons) )
            {
               SCIP_CALL( lockLinearVariable(scip, cons, consdata->linvars[i], consdata->lincoefs[i]) );
            }
         }
      }

      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         *infeasible = TRUE;
         goto TERMINATE;
      }
   }
   else if( SCIPisInfinity(scip, -constant) )
   {
      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         /* setting constraint rhs to infinity; this may change linear variable locks and events */
         for( i = 0; i < consdata->nlinvars; ++i )
         {
            if( SCIPconsIsLocked(cons) )
            {
               SCIP_CALL( unlockLinearVariable(scip, cons, consdata->linvars[i], consdata->lincoefs[i]) );
            }
            if( SCIPconsIsEnabled(cons) )
            {
               SCIP_CALL( dropLinearVarEvents(scip, cons, i) );
            }
         }

         consdata->rhs = SCIPinfinity(scip);

         for( i = 0; i < consdata->nlinvars; ++i )
         {
            if( SCIPconsIsEnabled(cons) )
            {
               SCIP_CALL( catchLinearVarEvents(scip, cons, i) );
            }
            if( SCIPconsIsLocked(cons) )
            {
               SCIP_CALL( lockLinearVariable(scip, cons, consdata->linvars[i], consdata->lincoefs[i]) );
            }
         }
      }
      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         *infeasible = TRUE;
         goto TERMINATE;
      }
   }
   else if( constant != 0.0 )
   {
      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         consdata->lhs -= constant;
         assert(!SCIPisInfinity(scip, REALABS(consdata->lhs)));
      }
      if( !SCIPisInfinity(scip,  consdata->rhs) )
      {
         consdata->rhs -= constant;
         assert(!SCIPisInfinity(scip, REALABS(consdata->rhs)));
      }
   }

TERMINATE:
   for( i = 0; i < nlinvars; ++i )
   {
      SCIP_CALL( addLinearCoef(scip, cons, linvars[i], lincoefs[i]) );
   }

   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &lincoefs);

   /* @todo linear variables that are also children of exprgraphnode could be moved into the expression graph for certain nonlinear operators (quadratic, polynomial), since that may allow better bound tightening */

   return SCIP_OKAY;
}

/** create a nonlinear row representation of the constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< nonlinear constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   if( consdata->nexprtrees == 0 )
   {
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            consdata->nlinvars, consdata->linvars, consdata->lincoefs,
            0, NULL, 0, NULL,
            NULL, consdata->lhs, consdata->rhs,
            consdata->curvature) );
   }
   else if( consdata->nexprtrees == 1 && consdata->nonlincoefs[0] == 1.0 )
   {
      assert(consdata->exprtrees[0] != NULL);
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            consdata->nlinvars, consdata->linvars, consdata->lincoefs,
            0, NULL, 0, NULL,
            consdata->exprtrees[0], consdata->lhs, consdata->rhs,
            consdata->curvature) );
   }
   else
   {
      /* since expression trees may share variable, we cannot easily sum them up,
       * but we can request a single expression tree from the expression graph
       */
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_EXPRTREE* exprtree;

      assert(consdata->exprgraphnode != NULL); /* since nexprtrees > 0 */
      conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
      assert(conshdlrdata != NULL);

      SCIP_CALL( SCIPexprgraphGetTree(conshdlrdata->exprgraph, consdata->exprgraphnode, &exprtree) );
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            consdata->nlinvars, consdata->linvars, consdata->lincoefs,
            0, NULL, 0, NULL,
            exprtree, consdata->lhs, consdata->rhs,
            consdata->curvature) );
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** tries to automatically convert a nonlinear constraint (or a part of it) into a more specific and more specialized constraint */
static
SCIP_RETCODE presolveUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_Bool*            upgraded,           /**< buffer to store whether constraint was upgraded */
   int*                  nupgdconss,         /**< buffer to increase if constraint was upgraded */
   int*                  naddconss           /**< buffer to increase with number of additional constraints created during upgrade */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS** upgdconss;
   int upgdconsssize;
   int nupgdconss_;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsModifiable(cons));
   assert(upgraded   != NULL);
   assert(nupgdconss != NULL);
   assert(naddconss  != NULL);

   *upgraded = FALSE;

   nupgdconss_ = 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if there are no upgrade methods, we can stop */
   if( conshdlrdata->nnlconsupgrades == 0 )
      return SCIP_OKAY;

   /* set debug solution in expression graph and evaluate nodes, so reformulation methods can compute debug solution values for new auxiliary variables */
#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      SCIP_Real* varvals;

      SCIP_CALL( SCIPallocBufferArray(scip, &varvals, SCIPexprgraphGetNVars(conshdlrdata->exprgraph)) );

      for( i = 0; i < SCIPexprgraphGetNVars(conshdlrdata->exprgraph); ++i )
         SCIP_CALL( SCIPdebugGetSolVal(scip, (SCIP_VAR*)SCIPexprgraphGetVars(conshdlrdata->exprgraph)[i], &varvals[i]) );

      SCIP_CALL( SCIPexprgraphEval(conshdlrdata->exprgraph, varvals) );

      SCIPfreeBufferArray(scip, &varvals);
   }
#endif

   upgdconsssize = 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &upgdconss, upgdconsssize) );

   /* call the upgrading methods */

   SCIPdebugMsg(scip, "upgrading nonlinear constraint <%s> (up to %d upgrade methods):\n",
      SCIPconsGetName(cons), conshdlrdata->nnlconsupgrades);
   SCIPdebugPrintCons(scip, cons, NULL);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nnlconsupgrades; ++i )
   {
      if( !conshdlrdata->nlconsupgrades[i]->active )
         continue;
      if( conshdlrdata->nlconsupgrades[i]->nlconsupgd == NULL )
         continue;

      SCIP_CALL( conshdlrdata->nlconsupgrades[i]->nlconsupgd(scip, cons, &nupgdconss_, upgdconss, upgdconsssize) );

      while( nupgdconss_ < 0 )
      {
         /* upgrade function requires more memory: resize upgdconss and call again */
         assert(-nupgdconss_ > upgdconsssize);
         upgdconsssize = -nupgdconss_;
         SCIP_CALL( SCIPreallocBufferArray(scip, &upgdconss, -nupgdconss_) );

         SCIP_CALL( conshdlrdata->nlconsupgrades[i]->nlconsupgd(scip, cons, &nupgdconss_, upgdconss, upgdconsssize) );

         assert(nupgdconss_ != 0);
      }

      if( nupgdconss_ > 0 )
      {
         /* got upgrade */
         int j;

         SCIPdebugMsg(scip, " -> upgraded to %d constraints:\n", nupgdconss_);

         /* add the upgraded constraints to the problem and forget them */
         for( j = 0; j < nupgdconss_; ++j )
         {
            SCIPdebugMsgPrint(scip, "\t");
            SCIPdebugPrintCons(scip, upgdconss[j], NULL);

            SCIP_CALL( SCIPaddCons(scip, upgdconss[j]) );      /*lint !e613*/
            SCIP_CALL( SCIPreleaseCons(scip, &upgdconss[j]) ); /*lint !e613*/
         }

         /* count the first upgrade constraint as constraint upgrade and the remaining ones as added constraints */
         *nupgdconss += 1;
         *naddconss += nupgdconss_ - 1;
         *upgraded = TRUE;

         /* delete upgraded constraint */
         SCIPdebugMsg(scip, "delete constraint <%s> after upgrade\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );

         break;
      }
   }

   SCIPfreeBufferArray(scip, &upgdconss);

   return SCIP_OKAY;
}

/** checks a nonlinear constraint for convexity and/or concavity */
static
SCIP_RETCODE checkCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_Bool             expensivechecks,    /**< whether also expensive checks should be executed */
   SCIP_Bool             assumeconvex        /**< whether to assume convexity in inequalities */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL* varbounds;
   int varboundssize;
   SCIP_VAR* var;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->iscurvchecked )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Checking curvature of constraint <%s>\n", SCIPconsGetName(cons));

   consdata->curvature = SCIP_EXPRCURV_LINEAR;
   consdata->iscurvchecked = TRUE;

   varbounds = NULL;
   varboundssize = 0;

   for( i = 0; i < consdata->nexprtrees; ++i )
   {
      assert(consdata->exprtrees[i] != NULL);
      assert(SCIPexprtreeGetNVars(consdata->exprtrees[i]) > 0 );

      if( assumeconvex )
      {
         /* for constraints a*f(x) <= rhs, we assume that it is convex */
         if( SCIPisInfinity(scip, -consdata->lhs) )
            consdata->curvatures[i] = SCIP_EXPRCURV_CONVEX;

         /* for constraints lhs <= a*f(x), we assume that it is concave */
         if( SCIPisInfinity(scip,  consdata->rhs) )
            consdata->curvatures[i] = SCIP_EXPRCURV_CONCAVE;
      }
      else
      {
         if( varboundssize == 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &varbounds, SCIPexprtreeGetNVars(consdata->exprtrees[i])) );
            varboundssize = SCIPexprtreeGetNVars(consdata->exprtrees[i]);
         }
         else if( varboundssize < SCIPexprtreeGetNVars(consdata->exprtrees[i]) )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &varbounds, SCIPexprtreeGetNVars(consdata->exprtrees[i])) );
            varboundssize = SCIPexprtreeGetNVars(consdata->exprtrees[i]);
         }
         assert(varbounds != NULL);

         for( j = 0; j < SCIPexprtreeGetNVars(consdata->exprtrees[i]); ++j )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[i])[j];
            SCIPintervalSetBounds(&varbounds[j],
               -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -MIN(SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))),    /*lint !e666*/
               +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  MAX(SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))) );  /*lint !e666*/
         }

         SCIP_CALL( SCIPexprtreeCheckCurvature(consdata->exprtrees[i], INTERVALINFTY, varbounds, &consdata->curvatures[i], NULL) );
         consdata->curvatures[i] = SCIPexprcurvMultiply(consdata->nonlincoefs[i], consdata->curvatures[i]);

         if( consdata->curvatures[i] == SCIP_EXPRCURV_UNKNOWN && SCIPconshdlrGetData(SCIPconsGetHdlr(cons))->isreformulated && SCIPexprGetOperator(SCIPexprtreeGetRoot(consdata->exprtrees[i])) != SCIP_EXPR_USER )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "indefinite expression tree in constraint <%s>\n", SCIPconsGetName(cons));
            SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(consdata->exprtrees[i], SCIPgetMessagehdlr(scip), NULL) ) );
            SCIPdebugMsgPrint(scip, "\n");
         }
      }

      /* @todo implement some more expensive checks */

      consdata->curvature = SCIPexprcurvAdd(consdata->curvature, consdata->curvatures[i]);

      SCIPdebugMsg(scip, "-> tree %d with coef %g is %s -> nonlinear part is %s\n", i, consdata->nonlincoefs[i], SCIPexprcurvGetName(consdata->curvatures[i]), SCIPexprcurvGetName(consdata->curvature));
   }

   SCIPfreeBufferArrayNull(scip, &varbounds);

   return SCIP_OKAY;
}  /*lint !e715*/

/* replaces a node by another node in expression graph
 * moves all parents of node to replacement
 * replaces all exprgraphnode's in constraints that are node by replacement
 * node may be freed, if captured only by given constraints
 */
static
SCIP_RETCODE reformReplaceNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  node,               /**< pointer to node to be replaced in expression graph */
   SCIP_EXPRGRAPHNODE*   replacement,        /**< node which takes node's place */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int c;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(*node != NULL);
   assert(replacement != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( SCIPexprgraphMoveNodeParents(exprgraph, node, replacement) );

   /* node was not captured by any constraint */
   if( *node == NULL )
      return SCIP_OKAY;

   /* if node still exists, then because it is captured by some constraint (it should not have parents anymore)
    * thus, look into the given constraints and replace their exprgraphnode by replacement
    * @todo may be expensive when this is done more often, esp. when we know that node will not be freed due to an added auxiliary constraint
    */
   assert(*node == NULL || SCIPexprgraphGetNodeNParents(*node) == 0);
   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( consdata->exprgraphnode == *node )
      {
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, &consdata->exprgraphnode) );
         consdata->exprgraphnode = replacement;
         SCIPexprgraphCaptureNode(replacement);

         /* since we change the node, also the constraint changes, so ensure that it is presolved again */
         consdata->ispresolved = FALSE;
      }
   }
   *node = NULL;

   return SCIP_OKAY;
}

/** creates a new auxiliary variable and a new auxiliary nonlinear constraint connecting the var and a given node
 * node is replaced by new new auxiliary variables node in all parents of node in expression graph and in all constraints that use node
 */
static
SCIP_RETCODE reformNode2Var(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_CONS**           conss,              /**< constraints where to update exprgraphnode */
   int                   nconss,             /**< number of constraints */
   int*                  naddcons,           /**< counter to increase when constraint is added */
   SCIP_Bool             donotmultaggr       /**< whether to mark auxiliary variable as not to multiaggregate */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* auxvar;
   SCIP_CONS* auxcons;
   SCIP_EXPRGRAPHNODE* auxvarnode;
   SCIP_INTERVAL bounds;
   SCIP_Real minusone;
   SCIP_Bool cutoff;

   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(naddcons != NULL);
   assert(SCIPexprgraphGetNodeDepth(node) >= 1); /* do not turn vars or consts into new vars */

   bounds = SCIPexprgraphGetNodeBounds(node);
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);

   SCIPdebugMsg(scip, "add auxiliary variable and constraint %s for node %p(%d,%d)\n", name, (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node));

   SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, SCIPintervalGetInf(bounds), SCIPintervalGetSup(bounds), 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar) );
   SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, &auxvarnode) );
#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      /* store debug sol value of node as value for auxvar in debug solution and as value for auxvarnode */
      SCIPexprgraphSetVarNodeValue(auxvarnode, SCIPexprgraphGetNodeVal(node));
      SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(node)) );
   }
#endif

   if( donotmultaggr )
   {
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, auxvar) );
   }

   /* set also bounds of auxvarnode to bounds, so it is available for new parent nodes (currently node->parents)
    * when updating their curvature information; avoid having to run domain propagation through exprgraph
    */
   SCIPexprgraphTightenNodeBounds(exprgraph, auxvarnode, bounds, BOUNDTIGHTENING_MINSTRENGTH, INTERVALINFTY, &cutoff);
   assert(!cutoff); /* we tightenend bounds from [-inf,+inf] to bounds, this should not be infeasible */

   /* add new constraint auxvar == node */
   minusone = -1.0;
   SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 1, &auxvar, &minusone, node, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
         FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, auxcons) );

   /* move parents of node in expression graph to auxvarnode
    * replace node by auxvarnode in constraints that use node */
   SCIP_CALL( reformReplaceNode(exprgraph, &node, auxvarnode, conss, nconss) );

   SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

   ++*naddcons;

   return SCIP_OKAY;
}

/** ensures that all children of a node have at least a given curvature by adding auxiliary variables */
static
SCIP_RETCODE reformEnsureChildrenMinCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_EXPRCURV         mincurv,            /**< minimal desired curvature */
   SCIP_CONS**           conss,              /**< constraints to check whether they use one of the replaced nodes */
   int                   nconss,             /**< number of constraints to check */
   int*                  naddcons            /**< counter to increase when constraint is added */
   )
{
   SCIP_EXPRGRAPHNODE* child;
   SCIP_Bool needupdate;

   int i;
   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(naddcons != NULL);
   assert(SCIPexprgraphGetNodeDepth(node) >= 1); /* do not turn vars or consts into new vars */
   assert(mincurv != SCIP_EXPRCURV_UNKNOWN); /* this is trivial and makes no sense */

   needupdate = FALSE; /* whether we need to update curvature of node */

   for( i = 0; i < SCIPexprgraphGetNodeNChildren(node); ++i )
   {
      child = SCIPexprgraphGetNodeChildren(node)[i];
      assert(child != NULL);

      if( (SCIPexprgraphGetNodeCurvature(child) & mincurv) != mincurv )
      {
         SCIPdebugMsg(scip, "add auxiliary variable for child %p(%d,%d) with curvature %s\n",
            (void*)child, SCIPexprgraphGetNodeDepth(child), SCIPexprgraphGetNodePosition(child), SCIPexprcurvGetName(SCIPexprgraphGetNodeCurvature(child)) );

         SCIP_CALL( reformNode2Var(scip, exprgraph, child, conss, nconss, naddcons, FALSE) );
         needupdate = TRUE;

         /* i'th child of node should now be a variable */
         assert(SCIPexprgraphGetNodeChildren(node)[i] != child);
         assert(SCIPexprgraphGetNodeOperator(SCIPexprgraphGetNodeChildren(node)[i]) == SCIP_EXPR_VARIDX);
      }

      assert(SCIPexprgraphGetNodeCurvature(SCIPexprgraphGetNodeChildren(node)[i]) & mincurv);
   }

   if( needupdate )
   {
      SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
      assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(node)));
   }

   return SCIP_OKAY;
}

/** reformulates a monomial by adding auxiliary variables and constraints for bilinear terms */
static
SCIP_RETCODE reformMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   int                   nfactors,           /**< number of factors */
   SCIP_EXPRGRAPHNODE**  factors,            /**< factors */
   SCIP_Real*            exponents,          /**< exponents, or NULL if all 1.0 */
   SCIP_EXPRGRAPHNODE**  resultnode,         /**< buffer to store node which represents the reformulated monomial */
   SCIP_Bool             createauxcons,      /**< whether to create auxiliary var/cons */
   int                   mindepth,           /**< minimal depth of new nodes in expression graph, or -1 */
   int*                  naddcons            /**< buffer to increase by number of added cons */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* auxvar;
   SCIP_CONS* auxcons;
   SCIP_Real minusone;

   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(nfactors > 0);
   assert(factors != NULL);
   assert(resultnode != NULL);
   assert(naddcons != NULL);

   /* factors are just one node */
   if( nfactors == 1 && (exponents == NULL || exponents[0] == 1.0) )
   {
      *resultnode = factors[0];
      return SCIP_OKAY;
   }

   /* only one factor, but with exponent < 0.0 and factor has mixed sign, e.g., x^(-3)
    * reformulate as auxvar * factor^(-exponent) = 1 and return the node for auxvar in resultnode
    */
   if( nfactors == 1 && exponents[0] < 0.0 && SCIPexprgraphGetNodeBounds(factors[0]).inf < 0.0 && SCIPexprgraphGetNodeBounds(factors[0]).sup > 0.0 )  /*lint !e613*/
   {
      SCIP_EXPRGRAPHNODE* auxnode;
      SCIP_EXPRGRAPHNODE* reformfactors[2];
      SCIP_Real reformexp[2];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);
      SCIPdebugMsg(scip, "add auxiliary variable and constraint %s\n", name);

      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );
      SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, resultnode) );

#ifdef WITH_DEBUG_SOLUTION
      /* store debug sol value of node as value for auxvar in debug solution and as value for resultnode */
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real debugval;
         debugval = pow(SCIPexprgraphGetNodeVal(factors[0]), exponents[0]);
         SCIPexprgraphSetVarNodeValue(*resultnode, debugval);
         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
      }
#endif

      /* increase naddcons before next call to reformMonomial, to avoid duplicate variable names, which is bad for debugging */
      ++*naddcons;

      /* add reformulation for resultnode(=auxvar) * factor^(-exponent) = 1.0
       * if exponent != -1.0, then factor^(-exponent) should be moved into extra variable
       * finally one should get an EXPR_MUL node */
      reformfactors[0] = *resultnode;
      reformfactors[1] = factors[0];
      reformexp[0] = 1.0;
      reformexp[1] = -exponents[0];  /*lint !e613*/
      SCIP_CALL( reformMonomial(scip, exprgraph, 2, reformfactors, reformexp, &auxnode, FALSE, mindepth, naddcons) );

      SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 0, NULL, NULL, auxnode, 1.0, 1.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, auxcons) );

      SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

      return SCIP_OKAY;
   }

   /* only one factor, but with exponent != 1.0 */
   if( nfactors == 1 )
   {
      /* create some power expression node, if not existing already */
      SCIP_EXPRGRAPHNODE* expnode;
      SCIP_EXPRGRAPHNODE* parent;
      int p;

      assert(exponents != NULL);

      /* check if there is already a node for factors[0]^exponents[0] */
      expnode = NULL;
      for( p = 0; p < SCIPexprgraphGetNodeNParents(factors[0]); ++p)
      {
         parent = SCIPexprgraphGetNodeParents(factors[0])[p];
         if( SCIPisIntegral(scip, exponents[0]) &&
            SCIPexprgraphGetNodeOperator(parent) == SCIP_EXPR_INTPOWER &&
            SCIPexprgraphGetNodeIntPowerExponent(parent) == (int)SCIPround(scip, exponents[0]) )
         {
            expnode = parent;
            break;
         }
         if( SCIPexprgraphGetNodeOperator(parent) == SCIP_EXPR_REALPOWER &&
            SCIPisEQ(scip, SCIPexprgraphGetNodeRealPowerExponent(parent), exponents[0]) )
         {
            expnode = parent;
         }
      }
      if( expnode == NULL )
      {
         if( SCIPisIntegral(scip, exponents[0]) )
            SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &expnode, SCIP_EXPR_INTPOWER, (int)SCIPround(scip, exponents[0])) );
         else
            SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &expnode, SCIP_EXPR_REALPOWER, exponents[0]) );

         SCIP_CALL( SCIPexprgraphAddNode(exprgraph, expnode, mindepth, 1, &factors[0]) );
         SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(expnode, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
         assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(expnode)));
      }

      if( createauxcons )
      {
         /* @todo if there was already a node for factors[0]^exponents[0], then there may have been also been already an auxiliary variable and constraint (-> ex7_3_4) */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);
         SCIPdebugMsg(scip, "add auxiliary variable and constraint %s\n", name);

         SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, auxvar) );
         SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, resultnode) );

#ifdef WITH_DEBUG_SOLUTION
         if( SCIPdebugIsMainscip(scip) )
         {
            SCIPexprgraphSetVarNodeValue(*resultnode, SCIPexprgraphGetNodeVal(expnode));
            SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(expnode)) );
         }
#endif

         /* add new constraint resultnode(=auxvar) = expnode */
         minusone = -1.0;
         SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 1, &auxvar, &minusone, expnode, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
               FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, auxcons) );

         SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

         ++*naddcons;
      }
      else
      {
         *resultnode = expnode;
      }

      return SCIP_OKAY;
   }

   if( nfactors == 2 && exponents != NULL && exponents[0] != 1.0 && exponents[0] == exponents[1] )  /*lint !e777*/
   {
      /* factor0^exponent * factor1^exponent with exponent != 1.0, reform as (factor0*factor1)^exponent */
      SCIP_EXPRGRAPHNODE* productnode;

      /* create node for factor0*factor1 */
      SCIP_CALL( reformMonomial(scip, exprgraph, 2, factors, NULL, &productnode, TRUE, mindepth, naddcons) );

      /* create node for productnode^exponents[0] by just calling this method again */
      SCIP_CALL( reformMonomial(scip, exprgraph, 1, &productnode, &exponents[0], resultnode, createauxcons, mindepth, naddcons) );

      return SCIP_OKAY;
   }

   if( nfactors == 2 && exponents != NULL && exponents[0] == -exponents[1] )  /*lint !e777*/
   {
      /* factor0^exponent * factor1^(-exponent), reform as (factor0/factor1)^exponent or (factor1/factor0)^(-exponent) */
      SCIP_EXPRGRAPHNODE* auxvarnode;
      SCIP_EXPRGRAPHNODE* auxconsnode;
      SCIP_EXPRGRAPHNODE* leftright[2];
      SCIP_Real absexp;

      /* create variable and constraint for factor0 = auxvar * factor1 (if exponent > 0) or factor1 = auxvar * factor0 (if exponent < 0) */

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);
      SCIPdebugMsg(scip, "add auxiliary variable and constraint %s\n", name);

      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );
      SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, &auxvarnode) );

#ifdef WITH_DEBUG_SOLUTION
      /* store debug sol value of node as value for auxvar in debug solution and as value for resultnode */
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real debugval;
         if( exponents[0] > 0.0 )
            debugval = SCIPexprgraphGetNodeVal(factors[0]) / SCIPexprgraphGetNodeVal(factors[1]);
         else
            debugval = SCIPexprgraphGetNodeVal(factors[1]) / SCIPexprgraphGetNodeVal(factors[0]);
         SCIPexprgraphSetVarNodeValue(auxvarnode, debugval);
         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
      }
#endif

      /* add new constraint resultnode(= auxvar) * factor1 - factor0 == 0 (exponent > 0) or auxvar * factor0 - factor1 == 0 (exponent < 0) */
      leftright[0] = auxvarnode;
      leftright[1] = exponents[0] > 0.0 ? factors[1] : factors[0];

      SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &auxconsnode, SCIP_EXPR_MUL, NULL) );
      SCIP_CALL( SCIPexprgraphAddNode(exprgraph, auxconsnode, -1, 2, leftright) );

      leftright[0] = auxconsnode;
      leftright[1] = exponents[0] > 0.0 ? factors[0] : factors[1];

      SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &auxconsnode, SCIP_EXPR_MINUS, NULL) );
      SCIP_CALL( SCIPexprgraphAddNode(exprgraph, auxconsnode, -1, 2, leftright) );

      SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 0, NULL, NULL, auxconsnode, 0.0, 0.0,
         TRUE, TRUE, TRUE, TRUE, TRUE,
         FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, auxcons) );

      SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

      ++*naddcons;

      /* create node for auxvarnode^abs(exponents[0]) by just calling this method again */
      absexp = fabs(exponents[0]);
      SCIP_CALL( reformMonomial(scip, exprgraph, 1, &auxvarnode, &absexp, resultnode, createauxcons, mindepth, naddcons) );

      return SCIP_OKAY;
   }

   /* @todo if nfactors > 2, assemble groups of factors with same exponent and replace these by a single variable first */

   {
      /* at least two factors */
      /* create auxvar for left half (recursively) and auxvar for right half (recursively) and maybe new auxvar for product */
      /* @todo it may be enough to replace single factors in a monomial to get it convex or concave, see Westerlund et.al. */

      SCIP_EXPRGRAPHNODE* productnode;
      SCIP_EXPRGRAPHNODE* leftright[2]; /* {left, right} */
      SCIP_EXPRGRAPHNODE* parent;
      int half;
      int p;

      half = nfactors / 2;
      assert(half > 0);
      assert(half < nfactors);

      SCIP_CALL( reformMonomial(scip, exprgraph, half, factors, exponents, &leftright[0], TRUE, mindepth, naddcons) );
      SCIP_CALL( reformMonomial(scip, exprgraph, nfactors-half, &factors[half], exponents != NULL ? &exponents[half] : NULL, &leftright[1], TRUE, mindepth, naddcons) );  /*lint !e826*/

      /* check if there is already a node for left * right */
      productnode = NULL;
      for( p = 0; p < SCIPexprgraphGetNodeNParents(leftright[0]); ++p)
      {
         parent = SCIPexprgraphGetNodeParents(factors[0])[p];
         if( SCIPexprgraphGetNodeOperator(parent) != SCIP_EXPR_MUL )
            continue;

         assert(SCIPexprgraphGetNodeNChildren(parent) == 2);
         if( (SCIPexprgraphGetNodeChildren(parent)[0] == leftright[0] && SCIPexprgraphGetNodeChildren(parent)[1] == leftright[1]) ||
            ( SCIPexprgraphGetNodeChildren(parent)[0] == leftright[1] && SCIPexprgraphGetNodeChildren(parent)[1] == leftright[0]) )
         {
            productnode = parent;
            break;
         }
      }
      if( productnode == NULL )
      {
         /* create node for left * right */
         SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &productnode, SCIP_EXPR_MUL, NULL) );
         SCIP_CALL( SCIPexprgraphAddNode(exprgraph, productnode, mindepth, 2, leftright) );
         SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(productnode, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
         assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(productnode)));
      }

      if( createauxcons )
      {
         /* @todo if there was already a node for factors[0]^exponents[0], then there may have been also been already an auxiliary variable and constraint (-> ex7_3_4) */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);
         SCIPdebugMsg(scip, "add auxiliary variable and constraint %s\n", name);

         SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, auxvar) );
         SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, resultnode) );

#ifdef WITH_DEBUG_SOLUTION
         if( SCIPdebugIsMainscip(scip) )
         {
            SCIPexprgraphSetVarNodeValue(*resultnode, SCIPexprgraphGetNodeVal(productnode));
            SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(productnode)) );
         }
#endif

         /* add new constraint resultnode(= auxvar) == left * right */
         minusone = -1.0;
         SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 1, &auxvar, &minusone, productnode, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
               FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, auxcons) );

         SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

         ++*naddcons;
      }
      else
      {
         *resultnode = productnode;
      }
   }

   return SCIP_OKAY;
}

/** reformulates expression graph into a form so that for each node under- and overestimators could be computed
 * similar to factorable reformulation in other global solvers, but sometimes does not split up complex operands (like quadratic)
 */
static
SCIP_RETCODE reformulate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  naddcons            /**< buffer to increase by the number of added constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPRGRAPH* exprgraph;
   SCIP_EXPRGRAPHNODE* node;
   SCIP_EXPRGRAPHNODE** children;
   SCIP_EXPRGRAPHNODE* reformnode;
   SCIP_Bool havenonlinparent;
   SCIP_Bool domainerror;
   int nchildren;
   int c;
   int d;
   int i;
   int u;
#ifndef NDEBUG
   int j;
#endif

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(naddcons != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
   assert(!SCIPinProbing(scip));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->isreformulated )
   {
      SCIPdebugMsg(scip, "skip reformulation, already done\n");
      return SCIP_OKAY;
   }

   exprgraph = conshdlrdata->exprgraph;

   /* make sure current variable bounds are variable nodes of exprgraph */
   SCIP_CALL( SCIPexprgraphPropagateVarBounds(exprgraph, INTERVALINFTY, FALSE, &domainerror) );
   assert(!domainerror); /* should have been found by domain propagation */

   /* set debug solution in expression graph and evaluate nodes, so we can compute debug solution values for auxiliary variables */
#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      SCIP_Real* varvals;

      SCIP_CALL( SCIPallocBufferArray(scip, &varvals, SCIPexprgraphGetNVars(exprgraph)) );

      for( i = 0; i < SCIPexprgraphGetNVars(exprgraph); ++i )
         SCIP_CALL( SCIPdebugGetSolVal(scip, (SCIP_VAR*)SCIPexprgraphGetVars(exprgraph)[i], &varvals[i]) );

      SCIP_CALL( SCIPexprgraphEval(exprgraph, varvals) );

      SCIPfreeBufferArray(scip, &varvals);
   }
#endif

   for( d = 1; d < SCIPexprgraphGetDepth(exprgraph); ++d )
   {
      i = 0;
      while( i < SCIPexprgraphGetNNodes(exprgraph)[d] )
      {
         node = SCIPexprgraphGetNodes(exprgraph)[d][i];
         assert(node != NULL);

         /* skip disabled nodes, they should be removed soon */
         if( !SCIPexprgraphIsNodeEnabled(node) )
         {
            ++i;
            continue;
         }

         /* make sure bounds and curvature of node are uptodate */
         SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
         assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(node)));

         /* try external reformulation methods */
         for( u = 0; u < conshdlrdata->nnlconsupgrades; ++u )
         {
            if( conshdlrdata->nlconsupgrades[u]->nodereform != NULL && conshdlrdata->nlconsupgrades[u]->active )
            {
               SCIP_CALL( conshdlrdata->nlconsupgrades[u]->nodereform(scip, exprgraph, node, naddcons, &reformnode) );
               if( reformnode == NULL )
                  continue;

               SCIPdebugMsg(scip, "external nodereform reformulated node %p(%d,%d), replace by %p\n",
                  (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node), (void*)reformnode);

               SCIP_CALL( reformReplaceNode(exprgraph, &node, reformnode, conss, nconss) );
               SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(reformnode, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
               assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(reformnode)));

               break;
            }
         }
         /* if node has been reformulated, continue with next node without increasing i */
         if( u < conshdlrdata->nnlconsupgrades )
            continue;

         /* leave nodes that are known to be convex/concave/linear as they are */
         if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
         {
            SCIPdebugMsg(scip, "skip reformulating node %p(%d,%d) = ", (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node));
            SCIPdebug( SCIPexprgraphPrintNode(node, SCIPgetMessagehdlr(scip), NULL) );
            SCIPdebugMsgPrint(scip, ", curv = %s\n", SCIPexprcurvGetName(SCIPexprgraphGetNodeCurvature(node)));
            ++i;
            continue;
         }

         /* get flag whether node has a nonlinear parent
          * we want to know whether the current node will be at the top of the tree after the next simplification run
          * due to the tricky reformulation of polynomials below, this may not be the case yet
          */
         havenonlinparent = SCIPexprgraphHasNodeNonlinearAncestor(node);

         /* take action */
         assert(SCIPexprgraphGetNodeCurvature(node) == SCIP_EXPRCURV_UNKNOWN);
         SCIPdebugMsg(scip, "think about reformulating %s node %p(%d,%d) = ", SCIPexpropGetName(SCIPexprgraphGetNodeOperator(node)), (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node));
         SCIPdebug( SCIPexprgraphPrintNode(node, SCIPgetMessagehdlr(scip), NULL) );
         SCIPdebugMsgPrint(scip, "\n");

         children  = SCIPexprgraphGetNodeChildren(node);
         nchildren = SCIPexprgraphGetNodeNChildren(node);
         assert(children != NULL || nchildren == 0);

#ifndef NDEBUG
         /* at this place, all children nodes should have a known curvature, except if they only appear only linearly in constraints
          * the latter for cases where constraints with unknown curvature are upgraded to other constraint handler that can handle these (quadratic, signpower,...)
          */
         for( j = 0; j < nchildren; ++j )
         {
            assert(children[j] != NULL);  /*lint !e613*/
            if( havenonlinparent ||
               (  SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_PLUS  &&
                  SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_MINUS &&
                  SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_SUM   &&
                  SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_LINEAR) )
               assert(SCIPexprgraphGetNodeCurvature(children[j]) != SCIP_EXPRCURV_UNKNOWN || SCIPexprgraphGetNodeOperator(children[j]) == SCIP_EXPR_USER);  /*lint !e613*/
         }
#endif


         switch( SCIPexprgraphGetNodeOperator(node) )
         {
         case SCIP_EXPR_VARIDX:
         case SCIP_EXPR_PARAM:
         case SCIP_EXPR_CONST:
            SCIPerrorMessage("node with operator %d cannot have unknown curvature\n", SCIPexprgraphGetNodeOperator(node));
            SCIPABORT();
            break;  /*lint !e527*/

            /* linear operands */
         case SCIP_EXPR_PLUS:
         case SCIP_EXPR_MINUS:
         case SCIP_EXPR_SUM:
         case SCIP_EXPR_LINEAR:
            /* children have conflicting curvature, we can handle such sums in cons_nonlinear
             * thus, turn node into variable, if it has nonlinear parents */
            if( havenonlinparent )
            {
               SCIP_CALL( reformNode2Var(scip, exprgraph, node, conss, nconss, naddcons, FALSE) );
               assert(node != NULL);
               assert(SCIPexprgraphGetNodeNParents(node) == 0); /* node should now be at top of graph */
            }
            ++i;
            break;

         /* quadratic operands */
         case SCIP_EXPR_MUL:
         case SCIP_EXPR_QUADRATIC:
         {
            SCIP_EXPRGRAPHNODE* child;
            SCIP_Bool needupdate;

            /* ensure all children are linear, so next simplifier run makes sure all children will be variables (by distributing the product)
             * however, that will not work for user-expressions, so we should also ensure that they are none (@todo as they are linear, they could actually be replaced by a regular linear expression)
             */
            SCIPdebugMessage("ensure children are linear\n");
            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_LINEAR, conss, nconss, naddcons) );

            needupdate = FALSE;  /* whether we need to update curvature of node */
            for( c = 0; c < SCIPexprgraphGetNodeNChildren(node); ++c )
            {
               child = SCIPexprgraphGetNodeChildren(node)[c];
               assert(child != NULL);

               if( SCIPexprgraphGetNodeCurvature(child) != SCIP_EXPRCURV_LINEAR || SCIPexprgraphGetNodeOperator(child) == SCIP_EXPR_USER )
               {
                  SCIPdebugMessage("add auxiliary variable for child %p(%d,%d) with curvature %s operator %s\n",
                     (void*)child, SCIPexprgraphGetNodeDepth(child), SCIPexprgraphGetNodePosition(child), SCIPexprcurvGetName(SCIPexprgraphGetNodeCurvature(child)), SCIPexpropGetName(SCIPexprgraphGetNodeOperator(child)) );

                  SCIP_CALL( reformNode2Var(scip, exprgraph, child, conss, nconss, naddcons, FALSE) );
                  needupdate = TRUE;

                  /* c'th child of node should now be a variable */
                  assert(SCIPexprgraphGetNodeChildren(node)[c] != child);
                  assert(SCIPexprgraphGetNodeOperator(SCIPexprgraphGetNodeChildren(node)[c]) == SCIP_EXPR_VARIDX);
               }
            }
            if( needupdate )
            {
               SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
               assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(node)));
            }

            if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
            {
               /* if curvature is now known then we are done */
               ++i;
               break;
            }

            /* if we have nonlinear parents or a sibling, then add auxiliary variable for this node, so an upgrade to cons_quadratic should take place
             * we assume that siblings are non-linear and non-quadratic, which should be the case if simplifier was run, and also if this node was created during reformulating a polynomial
             * @todo we could also add auxvars for the sibling nodes, e.g., if there is only one
             * @todo if sibling nodes are quadratic (or even linear) due to reformulation, then we do not need to reform here... (-> nvs16)
             *       maybe this step should not be done here at all if havenonlinparent is FALSE? e.g., move into upgrade from quadratic?
             */
            if( havenonlinparent || SCIPexprgraphHasNodeSibling(node) )
            {
               SCIP_CALL( reformNode2Var(scip, exprgraph, node, conss, nconss, naddcons, FALSE) );
               assert(node != NULL);
               assert(SCIPexprgraphGetNodeNParents(node) == 0); /* node should now be at top of graph, so it can be upgraded by cons_quadratic */
               break;
            }

            ++i;
            break;
         }

         case SCIP_EXPR_DIV:
         {
            /* reformulate as bilinear term
             * note that in the reformulation, a zero in the denominator forces the nominator to be zero too, but the auxiliary can be arbitrary
             */
            SCIP_EXPRGRAPHNODE* auxvarnode;
            SCIP_EXPRGRAPHNODE* auxnode;
            SCIP_EXPRGRAPHNODE* auxchildren[3];
            SCIP_Real           lincoefs[3];
            SCIP_QUADELEM       quadelem;
            SCIP_VAR*           auxvar;
            SCIP_CONS*          auxcons;
            char                name[SCIP_MAXSTRLEN];
            SCIP_INTERVAL       bounds;

            assert(children != NULL);
            assert(nchildren == 2);

            bounds = SCIPexprgraphGetNodeBounds(node);
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);

            SCIPdebugMsg(scip, "add auxiliary variable %s for division in node %p(%d,%d)\n", name, (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node));

            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, SCIPintervalGetInf(bounds), SCIPintervalGetSup(bounds), 0.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );
            SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, &auxvarnode) );

#ifdef WITH_DEBUG_SOLUTION
            if( SCIPdebugIsMainscip(scip) )
            {
               SCIP_Real debugval;
               debugval = SCIPexprgraphGetNodeVal(children[0]) / SCIPexprgraphGetNodeVal(children[1]);
               SCIPexprgraphSetVarNodeValue(auxvarnode, debugval);
               SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
            }
#endif

            /* add new constraint auxvar * child[1] - child[0] == 0 */
            auxchildren[0] = children[0];  /*lint !e613*/
            auxchildren[1] = children[1];  /*lint !e613*/
            auxchildren[2] = auxvarnode;

            lincoefs[0] = -1.0;
            lincoefs[1] =  0.0;
            lincoefs[2] =  0.0;

            quadelem.idx1 = 1;
            quadelem.idx2 = 2;
            quadelem.coef = 1.0;

            SCIP_CALL( SCIPexprgraphCreateNodeQuadratic(SCIPblkmem(scip), &auxnode, 3, lincoefs, 1, &quadelem, 0.0) );
            SCIP_CALL( SCIPexprgraphAddNode(exprgraph, auxnode, -1, 3, auxchildren) );

            SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 0, NULL, NULL, auxnode, 0.0, 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );

            /* replace node by auxvarnode in graph and constraints that use it */
            SCIP_CALL( reformReplaceNode(exprgraph, &node, auxvarnode, conss, nconss) );

            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

            ++*naddcons;

            /* do not increase i, since node was removed and not necessarily replaced here */
            break;
         }

         case SCIP_EXPR_MIN:
         {
            /* make sure that both children are concave, because min of concave functions is concave */
            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_CONCAVE, conss, nconss, naddcons) );
            assert(SCIPexprgraphGetNodeCurvature(node) & SCIP_EXPRCURV_CONCAVE);
            ++i;
            break;
         }

         case SCIP_EXPR_MAX:
         {
            /* make sure that both children are convex, because max of convex functions is convex */
            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_CONVEX, conss, nconss, naddcons) );
            assert(SCIPexprgraphGetNodeCurvature(node) & SCIP_EXPRCURV_CONVEX);
            ++i;
            break;

         }

         case SCIP_EXPR_INTPOWER:
         {
            assert(nchildren == 1);

            /* for an intpower with mixed sign in the base and negative exponent, we reformulate similar as for EXPR_DIV */
            if( SCIPexprgraphGetNodeIntPowerExponent(node) < 0 && SCIPintervalGetInf(SCIPexprgraphGetNodeBounds(children[0])) < 0.0 && SCIPintervalGetSup(SCIPexprgraphGetNodeBounds(children[0])) > 0.0 )  /*lint !e613*/
            {
               SCIP_EXPRGRAPHNODE* auxvarnode;
               SCIP_Real exponent;

               /* if we have something like x^(-3) with mixed sign for x, then add auxvar and reform as auxvar*x^3 = 1 via reformMonomial */
               exponent = (SCIP_Real)SCIPexprgraphGetNodeIntPowerExponent(node);
               SCIP_CALL( reformMonomial(scip, exprgraph, 1, children, &exponent, &auxvarnode, TRUE, SCIPexprgraphGetNodeDepth(node), naddcons) );
               /* replace node by auxvarnode */
               SCIP_CALL( reformReplaceNode(exprgraph, &node, auxvarnode, conss, nconss) );
               break;
            }

            /* otherwise, we continue as for other univariate operands */
         }   /*lint -fallthrough*/

         /* univariate operands where the child does not have bounds and curvature from which we can deduce curvature of this node,
          * but where we can do more if the child is linear
          * thus, turn child into auxiliary variable
          */
         case SCIP_EXPR_SQUARE:
         case SCIP_EXPR_SQRT:
         case SCIP_EXPR_EXP:
         case SCIP_EXPR_LOG:
         case SCIP_EXPR_ABS:
         case SCIP_EXPR_REALPOWER:
         case SCIP_EXPR_SIGNPOWER:
         {
            assert(nchildren == 1);

            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_LINEAR, conss, nconss, naddcons) );

            if( SCIPexprgraphGetNodeCurvature(node) == SCIP_EXPRCURV_UNKNOWN )
            {
               /* the only case where f(x) for a linear term x is indefinite here is if f is intpower or signpower and x has mixed sign */
               assert(SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_INTPOWER || SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_SIGNPOWER);
               assert(SCIPintervalGetInf(SCIPexprgraphGetNodeBounds(children[0])) < 0.0);  /*lint !e613*/
               assert(SCIPintervalGetSup(SCIPexprgraphGetNodeBounds(children[0])) > 0.0);  /*lint !e613*/
            }

            /* update curvature of node */
            SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
            assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(node)));

            if( SCIPexprgraphGetNodeCurvature(node) == SCIP_EXPRCURV_UNKNOWN )
            {
               /* if intpower and signpower with positive exponent and a mixed sign in the child bounds still does not give a curvature,
                * we can do more if we make this node the root of a nonlinear constraints expression node, so it can be upgraded by cons_signpower
                * of course, this is only required if the node is still intermediate
                *
                * an intpower with negative exponent should have been handled above
                * for signpower, we assume the exponent is > 1.0
                */
               assert(SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_INTPOWER || SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_SIGNPOWER);
               assert(SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_INTPOWER  || SCIPexprgraphGetNodeIntPowerExponent(node) > 1);
               assert(SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_SIGNPOWER || SCIPexprgraphGetNodeSignPowerExponent(node) > 1.0);
               if( havenonlinparent )
               {
                  SCIP_CALL( reformNode2Var(scip, exprgraph, node, conss, nconss, naddcons, FALSE) );
                  assert(node != NULL); /* it should be used by some auxiliary constraint now */
                  assert(SCIPexprgraphGetNodeNParents(node) == 0); /* node should now be at top of graph (and used by new auxiliary constraint) */
               }
            }
            ++i;

            break;
         }

         case SCIP_EXPR_SIN:
         case SCIP_EXPR_COS:
         case SCIP_EXPR_TAN:
         case SCIP_EXPR_SIGN:
            /* case SCIP_EXPR_ERF   : */
            /* case SCIP_EXPR_ERFI  : */
         {
            SCIPerrorMessage("no support for trigonometric or sign operands yet\n");
            return SCIP_ERROR;
         }

         case SCIP_EXPR_PRODUCT:
         {
            /* ensure all children are linear */
            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_LINEAR, conss, nconss, naddcons) );
            if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
            {
               ++i;
               break;
            }

            /* if curvature is still unknown (quite likely), then turn into a cascade of bilinear terms
             * if node has parents, then ensure that it has a known curvature, otherwise we are also fine with a node that is a product of two (aux)variables */
            SCIP_CALL( reformMonomial(scip, exprgraph, nchildren, children, NULL, &reformnode, havenonlinparent, SCIPexprgraphGetNodeDepth(node), naddcons) );

            /* replace node by reformnode in graph and in all constraints that use it */
            SCIP_CALL( reformReplaceNode(exprgraph, &node, reformnode, conss, nconss) );

            /* do not increase i, since node was removed and not necessarily replaced here */
            break;
         }

         case SCIP_EXPR_POLYNOMIAL:
         {
            /* if polynomial has several monomials, replace by a sum of nodes each having a single monomial and one that has all linear and quadratic monomials
             * if polynomial has only a single monomial, then reformulate that one
             */
            SCIP_EXPRDATA_MONOMIAL** monomials;
            SCIP_EXPRDATA_MONOMIAL* monomial;
            int nmonomials;
            SCIP_Real* exponents;
            SCIP_Real coef;
            int* childidxs;
            int nfactors;
            int f;
            SCIP_INTERVAL childbounds;
            SCIP_EXPRCURV childcurv;
            SCIP_Bool modified;

            monomials  = SCIPexprgraphGetNodePolynomialMonomials(node);
            nmonomials = SCIPexprgraphGetNodePolynomialNMonomials(node);
            assert(nmonomials >= 1); /* constant polynomials should have been simplified away */

            if( nmonomials > 1 )
            {
               SCIP_EXPRGRAPHNODE* sumnode;
               SCIP_Real constant;
               int nquadelems;
               SCIP_QUADELEM* quadelems;
               SCIP_Real* lincoefs;
               int nmonomialnodes;
               SCIP_EXPRGRAPHNODE** childrennew;
               SCIP_EXPRGRAPHNODE** monomialnodes;
               int m;

               /* @todo if a monomial is a factor of another monomial, then we could (and should?) replace it there by the node we create for it here -> ex7_2_1
                * @todo factorizing the polynomial could be beneficial
                */

               /* constant part of polynomials, to add to first monomialnode, if any, or quadratic or linear part */
               constant = SCIPexprgraphGetNodePolynomialConstant(node);

               /* coefficients from linear monomials */
               lincoefs = NULL;

               /* quadratic elements */
               nquadelems = 0;
               quadelems = NULL;

               /* expression graph nodes representing single higher-degree monomials, and single node with linear and/or quadratic monomials */
               nmonomialnodes = 0;
               SCIP_CALL( SCIPallocBufferArray(scip, &monomialnodes, nmonomials) );

               /* children of new monomial nodes that are setup */
               childrennew = NULL;

               for( m = 0; m < nmonomials; ++m )
               {
                  monomial = monomials[m];
                  assert(monomial != NULL);

                  coef = SCIPexprGetMonomialCoef(monomial);
                  exponents = SCIPexprGetMonomialExponents(monomial);
                  childidxs = SCIPexprGetMonomialChildIndices(monomial);
                  nfactors = SCIPexprGetMonomialNFactors(monomial);
                  assert(nfactors >= 1); /* constant monomials should have been simplified away */
                  assert(coef != 0.0);  /* zero-monomials should have been simplified away */

                  if( nfactors == 1 && exponents[0] == 1.0 )
                  {
                     /* linear monomial */
                     if( lincoefs == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nchildren) );
                        BMSclearMemoryArray(lincoefs, nchildren);
                     }
                     assert(0 <= childidxs[0] && childidxs[0] < nchildren);
                     assert(lincoefs[childidxs[0]] == 0.0); /* monomials should have been merged */
                     lincoefs[childidxs[0]] = coef;
                  }
                  else if( nfactors == 1 && exponents[0] == 2.0 )
                  {
                     /* square monomial */
                     if( quadelems == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nmonomials) );
                     }
                     quadelems[nquadelems].idx1 = childidxs[0];
                     quadelems[nquadelems].idx2 = childidxs[0];
                     quadelems[nquadelems].coef = coef;
                     ++nquadelems;
                  }
                  else if( nfactors == 2 && exponents[0] == 1.0 && exponents[1] == 1.0 )
                  {
                     /* bilinear monomial */
                     if( quadelems == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nmonomials) );
                     }
                     if( childidxs[0] < childidxs[1] )
                     {
                        quadelems[nquadelems].idx1 = childidxs[0];
                        quadelems[nquadelems].idx2 = childidxs[1];
                     }
                     else
                     {
                        quadelems[nquadelems].idx1 = childidxs[1];
                        quadelems[nquadelems].idx2 = childidxs[0];
                     }
                     quadelems[nquadelems].coef = coef;
                     ++nquadelems;
                  }
                  else
                  {
                     /* general monomial -> pass into separate expression graph node */
                     SCIP_EXPRDATA_MONOMIAL* monomialnew;

                     /* create new node for this monomial, children will be those associated with factors */
                     SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomialnew, coef, nfactors, NULL, exponents) );
                     SCIP_CALL( SCIPexprgraphCreateNodePolynomial(SCIPblkmem(scip), &monomialnodes[nmonomialnodes], 1, &monomialnew, constant, FALSE) );
                     constant = 0.0;

                     if( childrennew == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &childrennew, nchildren) );
                     }
                     assert(nfactors <= nchildren);
                     for( f = 0; f < nfactors; ++f )
                        childrennew[f] = children[childidxs[f]];  /*lint !e613*/

                     /* add new node to same depth as this node, so we will reformulate it during this run
                      * no need to refresh bounds/curvature here, since that will be done when we reach this node next */
                     SCIP_CALL( SCIPexprgraphAddNode(exprgraph, monomialnodes[nmonomialnodes], SCIPexprgraphGetNodeDepth(node), nfactors, childrennew) );

                     ++nmonomialnodes;
                  }
               }
               /* should have had at least one linear, quadratic, or general monomial */
               assert(lincoefs != NULL || nquadelems > 0 || nmonomialnodes > 0);

               if( nquadelems > 0 )
               {
                  /* create and add additional node for quadratic and linear part, simplifier should take care of removing unused children later */
                  SCIP_CALL( SCIPexprgraphCreateNodeQuadratic(SCIPblkmem(scip), &monomialnodes[nmonomialnodes], nchildren, lincoefs, nquadelems, quadelems, constant) );
                  constant = 0.0;
                  SCIP_CALL( SCIPexprgraphAddNode(exprgraph, monomialnodes[nmonomialnodes], SCIPexprgraphGetNodeDepth(node), nchildren, children) );
                  ++nmonomialnodes;
               }
               else if( lincoefs != NULL )
               {
                  /* create additional node for linear part, simplifier should take care of removing unused children later */
                  SCIP_CALL( SCIPexprgraphCreateNodeLinear(SCIPblkmem(scip), &monomialnodes[nmonomialnodes], nchildren, lincoefs, constant) );
                  constant = 0.0;
                  SCIP_CALL( SCIPexprgraphAddNode(exprgraph, monomialnodes[nmonomialnodes], SCIPexprgraphGetNodeDepth(node), nchildren, children) );
                  ++nmonomialnodes;
               }
               assert(constant == 0.0); /* the constant should have been used somewhere */

               SCIPfreeBufferArrayNull(scip, &lincoefs);
               SCIPfreeBufferArrayNull(scip, &quadelems);
               SCIPfreeBufferArrayNull(scip, &childrennew);

               assert(nmonomialnodes > 0);
               if( nmonomialnodes > 1 )
               {
                  /* add node for sum of monomials to expression graph */
                  SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &sumnode, nmonomialnodes == 2 ? SCIP_EXPR_PLUS : SCIP_EXPR_SUM) );
                  SCIP_CALL( SCIPexprgraphAddNode(exprgraph, sumnode, -1, nmonomialnodes, monomialnodes) );
               }
               else
               {
                  /* if only one monomial, then because polynomial was linear or quadratic... */
                  assert(SCIPexprgraphGetNodeOperator(monomialnodes[0]) == SCIP_EXPR_LINEAR || SCIPexprgraphGetNodeOperator(monomialnodes[0]) == SCIP_EXPR_QUADRATIC);
                  sumnode = monomialnodes[0];
               }
               SCIPfreeBufferArray(scip, &monomialnodes);

               /* replace node by sumnode, and we are done */
               SCIP_CALL( reformReplaceNode(exprgraph, &node, sumnode, conss, nconss) );

               SCIPdebugMsg(scip, "splitup polynomial into sum of %d nodes\n", nmonomialnodes);

               break;
            }

            /* reformulate a monomial such that it becomes convex or concave, if necessary */

            monomial = monomials[0];
            assert(monomial != NULL);

            coef = SCIPexprGetMonomialCoef(monomial);
            exponents = SCIPexprGetMonomialExponents(monomial);
            childidxs = SCIPexprGetMonomialChildIndices(monomial);
            nfactors = SCIPexprGetMonomialNFactors(monomial);
            assert(nfactors >= 1); /* constant monomials should have been simplified away */
            assert(coef != 0.0);  /* zero-monomials should have been simplified away */
            assert(children != NULL);

            /* check if we make monomial convex or concave by making a child linear */
            modified = FALSE;
            if( nfactors == 1 )
            {
               /* ensure that the child of an univariate monomial is linear if its current (bounds,curvature) yields an unknown curvature for the monomial
                * and with linear child it had a known curvature (rules out x^a, a negative, x not linear) */
               childcurv = SCIPexprgraphGetNodeCurvature(children[childidxs[0]]);  /*lint !e613*/
               childbounds = SCIPexprgraphGetNodeBounds(children[childidxs[0]]);  /*lint !e613*/
               assert(SCIPexprcurvPower(childbounds, childcurv, exponents[0]) == SCIP_EXPRCURV_UNKNOWN); /* this is exactly the curvature of the node, which is unknown */

               /* if monomial were convex or concave if child were linear, then make child linear */
               if( SCIPexprcurvPower(childbounds, SCIP_EXPRCURV_LINEAR, exponents[0]) != SCIP_EXPRCURV_UNKNOWN )
               {
                  assert(childcurv != SCIP_EXPRCURV_LINEAR);
                  SCIPdebugMsg(scip, "reform child %d (univar. monomial) with curv %s into var\n", childidxs[0], SCIPexprcurvGetName(childcurv));
                  SCIP_CALL( reformNode2Var(scip, exprgraph, children[childidxs[0]], conss, nconss, naddcons, FALSE) );  /*lint !e613*/
                  modified = TRUE;
               }
            }
            else
            {
               /* check if the conditions on the exponents allow for a convex or concave monomial, assuming that the children are linear
                * if one of these conditions is fulfilled but a children curvature does not fit, then make these children linear
                */
               int nnegative;
               int npositive;
               SCIP_Real sum;
               SCIP_Bool expcurvpos;
               SCIP_Bool expcurvneg;
               SCIP_EXPRCURV desiredcurv;

               nnegative = 0; /* number of negative exponents */
               npositive = 0; /* number of positive exponents */
               sum = 0.0;     /* sum of exponents */
               expcurvpos = TRUE; /* whether exp_j * f_j''(x) >= 0 for all factors (assuming f_j >= 0) */
               expcurvneg = TRUE; /* whether exp_j * f_j''(x) <= 0 for all factors (assuming f_j >= 0) */

               /* ensure that none of the children have unknown curvature */
               for( c = 0; c < SCIPexprgraphGetNodeNChildren(node); ++c )
               {
                  childcurv = SCIPexprgraphGetNodeCurvature(children[c]);  /*lint !e613*/
                  if( childcurv == SCIP_EXPRCURV_UNKNOWN )
                  {
                     SCIPdebugMessage("reform child %d with unknown curvature into var\n", c);
                     SCIP_CALL( reformNode2Var(scip, exprgraph, children[c], conss, nconss, naddcons, FALSE) );  /*lint !e613*/
                     modified = TRUE;
                  }
               }
               if( modified )
               {
                  /* refresh curvature information in node, since we changed children */
                  SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
                  assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(node)));

                  modified = FALSE;
               }

               for( f = 0; f < nfactors; ++f )
               {
                  childcurv = SCIPexprgraphGetNodeCurvature(children[childidxs[f]]);  /*lint !e613*/
                  assert(childcurv != SCIP_EXPRCURV_UNKNOWN);
                  childbounds = SCIPexprgraphGetNodeBounds(children[childidxs[f]]);  /*lint !e613*/
                  if( childbounds.inf < 0.0 && childbounds.sup > 0.0 )
                     break;

                  if( exponents[f] < 0.0 )
                     ++nnegative;
                  else
                     ++npositive;
                  sum += exponents[f];

                  /* negate curvature if factor is negative */
                  if( childbounds.inf < 0.0 )
                     childcurv = SCIPexprcurvNegate(childcurv);

                  /* check if exp_j * checkcurv is convex (>= 0) and/or concave */
                  childcurv = SCIPexprcurvMultiply(exponents[f], childcurv);
                  if( !(childcurv & SCIP_EXPRCURV_CONVEX) )
                     expcurvpos = FALSE;
                  if( !(childcurv & SCIP_EXPRCURV_CONCAVE) )
                     expcurvneg = FALSE;
               }

               /* if some child can be both positive and negative, then nothing we can do here to get the monomial convex or concave
                * otherwise (i.e., f == nfactors), look further */
               desiredcurv = SCIP_EXPRCURV_UNKNOWN;
               if( f == nfactors )
               {
                  /* if all factors are linear, then a product f_j^exp_j with f_j >= 0 is convex if
                   * - all exponents are negative, or
                   * - all except one exponent j* are negative and exp_j* >= 1 - sum_{j!=j*}exp_j, but the latter is equivalent to sum_j exp_j >= 1
                   * further, the product is concave if
                   * - all exponents are positive and the sum of exponents is <= 1.0
                   *
                   * if factors are nonlinear, then we require additionally, that for convexity
                   * - each factor is convex if exp_j >= 0, or concave if exp_j <= 0, i.e., exp_j*f_j'' >= 0
                   * and for concavity, we require that
                   * - all factors are concave, i.e., exp_j*f_j'' <= 0
                   */

                  if( nnegative == nfactors || (nnegative == nfactors-1 && SCIPisGE(scip, sum, 1.0)) )
                  {
                     /* if exponents are such that we can be convex, but children curvature does not fit, make some children linear */
                     SCIPdebugMsg(scip, "%d-variate monomial is convex (modulo sign), child curv fits = %u\n", nfactors, expcurvpos);
                     /* since current node curvature is set to unknown, there must be such a child, since otherwise the node curvature had to be convex */
                     assert(!expcurvpos);
                     desiredcurv = SCIP_EXPRCURV_CONVEX;
                  }
                  else if( npositive == nfactors && SCIPisLE(scip, sum, 1.0) )
                  {
                     /* if exponents are such that we can be concave, but children curvature does not fit, make some children linear */
                     SCIPdebugMsg(scip, "%d-variate monomial is concave (modulo sign), child curv fits = %u\n", nfactors, expcurvneg);
                     /* since current node curvature is set to unknown, there must be such a child, since otherwise the node curvature had to be concave */
                     assert(!expcurvneg);
                     desiredcurv = SCIP_EXPRCURV_CONCAVE;
                  }
                  else
                  {
                     /* exponents are such that monomial is neither convex nor concave even if children were linear
                      * thus, reformulate monomial below
                      */
                  }
               }

               if( desiredcurv != SCIP_EXPRCURV_UNKNOWN )
               {
                  for( f = 0; f < nfactors; ++f )
                  {
                     childcurv = SCIPexprgraphGetNodeCurvature(children[childidxs[f]]);  /*lint !e613*/
                     assert(childcurv != SCIP_EXPRCURV_UNKNOWN);
                     childbounds = SCIPexprgraphGetNodeBounds(children[childidxs[f]]);  /*lint !e613*/
                     assert(childbounds.inf >= 0.0 || childbounds.sup <= 0.0);

                     /* negate curvature if factor is negative */
                     if( childbounds.inf < 0.0 )
                        childcurv = SCIPexprcurvNegate(childcurv);

                     /* check if exp_j * checkcurv is convex (>= 0) and/or concave */
                     childcurv = SCIPexprcurvMultiply(SCIPexprGetMonomialExponents(monomial)[f], childcurv);
                     if( (desiredcurv == SCIP_EXPRCURV_CONVEX  && !(childcurv & SCIP_EXPRCURV_CONVEX )) ||
                        (desiredcurv == SCIP_EXPRCURV_CONCAVE && !(childcurv & SCIP_EXPRCURV_CONCAVE)) )
                     {
                        SCIPdebugMsg(scip, "reform child %d (factor %d) with curv %s into var\n",
                           childidxs[f], f, SCIPexprcurvGetName(SCIPexprgraphGetNodeCurvature(children[childidxs[f]])));  /*lint !e613*/
                        SCIP_CALL( reformNode2Var(scip, exprgraph, children[childidxs[f]], conss, nconss, naddcons, FALSE) );  /*lint !e613*/
                        modified = TRUE;
                     }
                  }
               }
            }

            if( modified )
            {
               /* refresh curvature information in node, since we changed children, it should be convex or concave now */
               SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
               assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(node)));
               assert(SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN);

               /* we are done and can proceed with the next node */
               ++i;
               break;
            }

            /* monomial can only have unknown curvature here, if it has several factors
             * or is of form x^a with x both negative and positive and a an odd or negative integer (-> INTPOWER expression)
             */
            assert(nfactors > 1 ||
               (SCIPexprgraphGetNodeBounds(children[childidxs[0]]).inf < 0.0 && SCIPexprgraphGetNodeBounds(children[childidxs[0]]).sup > 0.0 &&
                  SCIPisIntegral(scip, exponents[0]) && (exponents[0] < 0.0 || ((int)SCIPround(scip, exponents[0]) % 2 != 0)))
               );  /*lint !e613*/

            /* bilinear monomials should not come up here, since simplifier should have turned them into quadratic expression nodes */
            assert(!(nfactors == 2 && exponents[0] == 1.0 && exponents[1] == 1.0));

            /* reform monomial if it is a product, or we need it to be on the top of the graph, or if it of the form x^a with a < 0.0 (and thus x having mixed sign, see assert above)
             * thus, in the case x^a with a an odd positive integer we assume that cons_signpower will do something */
            if( nfactors > 1 || havenonlinparent || exponents[0] < 0.0 )
            {
               SCIP_EXPRGRAPHNODE* auxnode;
               SCIP_EXPRGRAPHNODE** factors;

               if( nfactors > 1 )
               {
                  SCIP_CALL( SCIPallocBufferArray(scip, &factors, nfactors) );
                  for( f = 0; f < nfactors; ++f )
                     factors[f] = children[childidxs[f]];  /*lint !e613*/
               }
               else
                  factors = &children[childidxs[0]];  /*lint !e613*/

               SCIPdebugMsg(scip, "reform monomial node, create auxvar = %u\n", havenonlinparent);
               /* get new auxnode for monomial
                * if node has parents and monomial is of indefinite form x^a, then also create auxvar for it, since otherwise we create a auxnode with unknown curvature
                * note, that the case x^a with positive and odd a will still give an indefinite node (without parents), where we assume that signpower will pick it up at some point
                */
               SCIP_CALL( reformMonomial(scip, exprgraph, nfactors, factors, exponents, &auxnode, havenonlinparent, SCIPexprgraphGetNodeDepth(node), naddcons) );

               if( nfactors > 1 )
               {
                  SCIPfreeBufferArray(scip, &factors);
               }

               /* create node for monomialcoef * auxnode + monomialconstant, if not identical to auxnode */
               if( SCIPexprgraphGetNodePolynomialConstant(node) != 0.0 || coef != 1.0 )
               {
                  SCIP_EXPRGRAPHNODE* replnode;

                  SCIP_CALL( SCIPexprgraphCreateNodeLinear(SCIPblkmem(scip), &replnode, 1, &coef, SCIPexprgraphGetNodePolynomialConstant(node)) );
                  SCIP_CALL( SCIPexprgraphAddNode(exprgraph, replnode, -1, 1, &auxnode) );
                  auxnode = replnode;
               }

               /* replace node by auxnode and refresh its curvature */
               SCIP_CALL( reformReplaceNode(exprgraph, &node, auxnode, conss, nconss) );
               SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(auxnode, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
               assert(!SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(auxnode)));

               break;
            }
            else
            {
               SCIPdebugMsg(scip, "no reformulation of monomial node, assume signpower will take care of it\n");
            }

            ++i;
            break;
         }

         case SCIP_EXPR_USER:
         {
            /* ensure all children are linear */
            SCIP_CALL( reformEnsureChildrenMinCurvature( scip, exprgraph, node, SCIP_EXPRCURV_LINEAR, conss, nconss, naddcons ) );

            /* unknown curvature can be handled by user estimator callback or interval gradient */
            /*
            if( SCIPexprgraphGetNodeCurvature( node ) == SCIP_EXPRCURV_UNKNOWN )
            {
               SCIPerrorMessage("user expression with unknown curvature not supported\n");
               return SCIP_ERROR;
            }
            */

            ++i;
            break;
         }

         case SCIP_EXPR_LAST:
            SCIPABORT();
            break;
         }
      }
   }

   /* for constraints with concave f(g(x)) with linear g:R^n -> R, n>1, reformulate to get a univariate concave function, since this is easier to underestimate
    * @todo this does not work yet for sums of functions other than polynomials
    */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_EXPRGRAPHNODE* multivarnode;
      SCIP_EXPRCURV curv;

      assert(conss[c] != NULL);  /*lint !e613*/

      /* skip constraints that are to be deleted */
      if( SCIPconsIsDeleted(conss[c]) )  /*lint !e613*/
         continue;

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( consdata->exprgraphnode == NULL )
         continue;

      /* after reformulation, force a round of backpropagation in expression graph for all constraints,
       * since new variables (nlreform*) may now be used in existing constraints and we want domain restrictions
       * of operators propagated for these variables
       */
      consdata->forcebackprop = TRUE;

      if( SCIPexprgraphGetNodeOperator(consdata->exprgraphnode) == SCIP_EXPR_POLYNOMIAL )
      {
         SCIP_EXPRDATA_MONOMIAL* monomial;
         int m;
         int f;

         for( m = 0; m < SCIPexprgraphGetNodePolynomialNMonomials(consdata->exprgraphnode); ++m )
         {
            SCIP_CALL( SCIPexprgraphGetNodePolynomialMonomialCurvature(consdata->exprgraphnode, m, INTERVALINFTY, &curv) );

            monomial = SCIPexprgraphGetNodePolynomialMonomials(consdata->exprgraphnode)[m];
            assert(monomial != NULL);

            /* if nothing concave, then continue */
            if( (SCIPisInfinity(scip,  consdata->rhs) || curv != SCIP_EXPRCURV_CONCAVE) &&
               ( SCIPisInfinity(scip, -consdata->lhs) || curv != SCIP_EXPRCURV_CONVEX) )
               continue;

            for( f = 0; f < SCIPexprGetMonomialNFactors(monomial); ++f )
            {
               multivarnode = SCIPexprgraphGetNodeChildren(consdata->exprgraphnode)[SCIPexprGetMonomialChildIndices(monomial)[f]];

               /* search for a descendant of node that has > 1 children
                * after simplifier run, there should be no constant expressions left
                */
               while( SCIPexprgraphGetNodeNChildren(multivarnode) == 1 )
                  multivarnode = SCIPexprgraphGetNodeChildren(multivarnode)[0];

               /* if node expression is obviously univariate, then continue */
               if( SCIPexprgraphGetNodeNChildren(multivarnode) == 0 )
               {
                  assert(SCIPexprgraphGetNodeOperator(multivarnode) == SCIP_EXPR_CONST || SCIPexprgraphGetNodeOperator(multivarnode) == SCIP_EXPR_VARIDX);
                  continue;
               }

               /* if multivarnode is a linear expression, then replace this by an auxiliary variable/node
                * mark auxiliary variable as not to multiaggregate, so SCIP cannot undo what we just did
                */
               if( SCIPexprgraphGetNodeCurvature(multivarnode) == SCIP_EXPRCURV_LINEAR )
               {
                  SCIPdebugMsg(scip, "replace linear multivariate node %p(%d,%d) in expression of cons <%s> by auxvar\n",
                     (void*)multivarnode, SCIPexprgraphGetNodeDepth(multivarnode), SCIPexprgraphGetNodePosition(multivarnode), SCIPconsGetName(conss[c]));  /*lint !e613*/
                  SCIPdebugPrintCons(scip, conss[c], NULL);  /*lint !e613*/
                  SCIP_CALL( reformNode2Var(scip, exprgraph, multivarnode, conss, nconss, naddcons, TRUE) );
               }
            }
         }
      }
      else
      {
         curv = SCIPexprgraphGetNodeCurvature(consdata->exprgraphnode);

         /* if nothing concave, then continue */
         if( (SCIPisInfinity(scip,  consdata->rhs) || curv != SCIP_EXPRCURV_CONCAVE) &&
            ( SCIPisInfinity(scip, -consdata->lhs) || curv != SCIP_EXPRCURV_CONVEX) )
            continue;

         /* search for a descendant of node that has > 1 children
          * after simplifier run, there should be no constant expressions left
          */
         multivarnode = consdata->exprgraphnode;
         while( SCIPexprgraphGetNodeNChildren(multivarnode) == 1 )
            multivarnode = SCIPexprgraphGetNodeChildren(multivarnode)[0];

         /* if node expression is obviously univariate, then continue */
         if( SCIPexprgraphGetNodeNChildren(multivarnode) == 0 )
         {
            assert(SCIPexprgraphGetNodeOperator(multivarnode) == SCIP_EXPR_CONST || SCIPexprgraphGetNodeOperator(multivarnode) == SCIP_EXPR_VARIDX);
            continue;
         }

         /* if node itself is multivariate, then continue */
         if( multivarnode == consdata->exprgraphnode )
            continue;

         /* if multivarnode is a linear expression, then replace this by an auxiliary variable/node
          * mark auxiliary variable as not to multiaggregate, so SCIP cannot undo what we just did
          */
         if( SCIPexprgraphGetNodeCurvature(multivarnode) == SCIP_EXPRCURV_LINEAR )
         {
            SCIPdebugMsg(scip, "replace linear multivariate node %p(%d,%d) in expression of cons <%s> by auxvar\n",
               (void*)multivarnode, SCIPexprgraphGetNodeDepth(multivarnode), SCIPexprgraphGetNodePosition(multivarnode), SCIPconsGetName(conss[c]));  /*lint !e613*/
            SCIPdebugPrintCons(scip, conss[c], NULL);  /*lint !e613*/
            SCIP_CALL( reformNode2Var(scip, exprgraph, multivarnode, conss, nconss, naddcons, TRUE) );
         }
      }
   }

   conshdlrdata->isreformulated = TRUE;

   return SCIP_OKAY;
}

/** computes activity and violation of a constraint
 *
 * During presolving and if the constraint is active, it is assumes that SCIPexprgraphEval has been called for sol before.
 *
 * If a solution is found to violate the variable bounds, then violation calculation is stopped and solviolbounds is set to TRUE.
 */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool*            solviolbounds       /**< buffer to indicate whether solution is found to violate variable bounds by more than feastol */
   )
{  /*lint --e{666}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real varval;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(solviolbounds != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprinterpreter != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->activity = 0.0;
   consdata->lhsviol = 0.0;
   consdata->rhsviol = 0.0;
   varval = 0.0;
   *solviolbounds = FALSE;

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_Real activity;

      var = consdata->linvars[i];
      varval = SCIPgetSolVal(scip, sol, var);

      /* project onto local box, in case the LP solution is slightly outside the bounds (which is not our job to enforce) */
      if( sol == NULL )
      {
         /* with non-initial columns, this might fail because variables can shortly be a column variable before entering the LP and have value 0.0 in this case */
         if( (!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) && !SCIPisFeasGE(scip, varval, SCIPvarGetLbLocal(var))) ||
             (!SCIPisInfinity(scip,  SCIPvarGetUbLocal(var)) && !SCIPisFeasLE(scip, varval, SCIPvarGetUbLocal(var))) )
         {
            *solviolbounds = TRUE;
            return SCIP_OKAY;
         }
         varval = MAX(SCIPvarGetLbLocal(var), MIN(SCIPvarGetUbLocal(var), varval));
      }
      activity = consdata->lincoefs[i] * varval;

      /* the contribution of a variable with |varval| = +inf is +inf when activity > 0.0, -inf when activity < 0.0, and
       * 0.0 otherwise
       */
      if( SCIPisInfinity(scip, REALABS(varval)) )
      {
         if( activity > 0.0 && !SCIPisInfinity(scip, consdata->rhs) )
         {
            consdata->activity = SCIPinfinity(scip);
            consdata->rhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }

         if( activity < 0.0 && !SCIPisInfinity(scip, -consdata->lhs) )
         {
            consdata->activity = -SCIPinfinity(scip);
            consdata->lhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }
      }

      consdata->activity += activity;
   }

   for( i = 0; i < consdata->nexprtrees; ++i )
   {
      SCIP_Real activity;
      SCIP_Real val;
      int nvars;

      /* compile expression tree, if not done before */
      if( SCIPexprtreeGetInterpreterData(consdata->exprtrees[i]) == NULL )
      {
         SCIP_CALL( SCIPexprintCompile(conshdlrdata->exprinterpreter, consdata->exprtrees[i]) );
      }

      nvars = SCIPexprtreeGetNVars(consdata->exprtrees[i]);

      if( nvars == 1 )
      {
         /* in the not so unusual case that an expression has only one variable, we do not need to extra allocate memory */
         var = SCIPexprtreeGetVars(consdata->exprtrees[i])[0];
         varval = SCIPgetSolVal(scip, sol, var);

         /* project onto local box, in case the LP solution is slightly outside the bounds (and then cannot be evaluated) */
         if( sol == NULL )
         {
            /* with non-initial columns, this might fail because variables can shortly be a column variable before entering the LP and have value 0.0 in this case */
            if( (!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) && !SCIPisFeasGE(scip, varval, SCIPvarGetLbLocal(var))) ||
                (!SCIPisInfinity(scip,  SCIPvarGetUbLocal(var)) && !SCIPisFeasLE(scip, varval, SCIPvarGetUbLocal(var))) )
            {
               *solviolbounds = TRUE;
               return SCIP_OKAY;
            }
            varval = MAX(SCIPvarGetLbLocal(var), MIN(SCIPvarGetUbLocal(var), varval));
         }

         SCIP_CALL( SCIPexprintEval(conshdlrdata->exprinterpreter, consdata->exprtrees[i], &varval, &val) ); /* coverity ignore ARRAY_VS_SINGLETON warning */
      }
      else
      {
         SCIP_Real* x;
         int j;

         SCIP_CALL( SCIPallocBufferArray(scip, &x, nvars) );

         for( j = 0; j < nvars; ++j )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[i])[j];
            varval = SCIPgetSolVal(scip, sol, var);

            /* project onto local box, in case the LP solution is slightly outside the bounds (and then cannot be evaluated) */
            if( sol == NULL )
            {
               /* with non-initial columns, this might fail because variables can shortly be a column variable before entering the LP and have value 0.0 in this case */
               if( (!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) && !SCIPisFeasGE(scip, varval, SCIPvarGetLbLocal(var))) ||
                   (!SCIPisInfinity(scip,  SCIPvarGetUbLocal(var)) && !SCIPisFeasLE(scip, varval, SCIPvarGetUbLocal(var))) )
               {
                  *solviolbounds = TRUE;
                  SCIPfreeBufferArray(scip, &x);
                  return SCIP_OKAY;
               }
               varval = MAX(SCIPvarGetLbLocal(var), MIN(SCIPvarGetUbLocal(var), varval));
            }

            x[j] = varval;
         }

         SCIP_CALL( SCIPexprintEval(conshdlrdata->exprinterpreter, consdata->exprtrees[i], x, &val) );

         SCIPfreeBufferArray(scip, &x);
      }

      /* set the activity to infinity if a function evaluation was not valid (e.g., sqrt(-1) ) */
      if( !SCIPisFinite(val) )
      {
         consdata->activity = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }

      /* the contribution of an expression with |val| = +inf is +inf when its coefficient is > 0.0, -inf when its coefficient is < 0.0, and
       * 0.0 otherwise
       */
      activity = consdata->nonlincoefs[i] * val;
      if( SCIPisInfinity(scip, REALABS(val)) )
      {
         if( activity > 0.0 && !SCIPisInfinity(scip, consdata->rhs) )
         {
            consdata->activity = SCIPinfinity(scip);
            consdata->rhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }

         if( activity < 0.0  && !SCIPisInfinity(scip, -consdata->lhs) )
         {
            consdata->activity = -SCIPinfinity(scip);
            consdata->lhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }
      }

      consdata->activity += activity;
   }

   if( consdata->nexprtrees == 0 && consdata->exprgraphnode != NULL )
   {
      SCIP_Real val;

      assert(SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) <= SCIP_STAGE_EXITPRESOLVE);

      val = SCIPexprgraphGetNodeVal(consdata->exprgraphnode);
      assert(val != SCIP_INVALID);  /*lint !e777*/

      /* set the activity to infinity if a function evaluation was not valid (e.g., sqrt(-1) ) */
      if( !SCIPisFinite(val) )
      {
         consdata->activity = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }

      if( SCIPisInfinity(scip, val) && !SCIPisInfinity(scip, consdata->rhs) )
      {
         consdata->activity = SCIPinfinity(scip);
         consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      else if( SCIPisInfinity(scip, -val) && !SCIPisInfinity(scip, -consdata->lhs) )
      {
         consdata->activity = -SCIPinfinity(scip);
         consdata->lhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }

      consdata->activity += val;
   }

   if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisGT(scip, consdata->lhs - consdata->activity, SCIPfeastol(scip)) )
      consdata->lhsviol = consdata->lhs - consdata->activity;
   else
      consdata->lhsviol = 0.0;

   if( !SCIPisInfinity(scip,  consdata->rhs) && SCIPisGT(scip, consdata->activity - consdata->rhs, SCIPfeastol(scip)) )
      consdata->rhsviol = consdata->activity - consdata->rhs;
   else
      consdata->rhsviol = 0.0;

   /* update absolute and relative violation of the solution */
   if( sol != NULL )
   {
      SCIP_Real absviol;
      SCIP_Real relviol;
      SCIP_Real lhsrelviol;
      SCIP_Real rhsrelviol;

      absviol = MAX(consdata->lhsviol, consdata->rhsviol);

      lhsrelviol = SCIPrelDiff(consdata->lhs, consdata->activity);
      rhsrelviol = SCIPrelDiff(consdata->activity, consdata->rhs);
      relviol = MAX(lhsrelviol, rhsrelviol);

      SCIPupdateSolConsViolation(scip, sol, absviol, relviol);
   }

   return SCIP_OKAY;
}

/** computes violation of a set of constraints
 *
 * If the solution is found to violate bounds of some variable in some constraint, then violation computation is stopped and solviolbounds is set to TRUE.
 */
static
SCIP_RETCODE computeViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool*            solviolbounds,      /**< buffer to indicate whether solution violates bounds of some variable by more than feastol */
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
   assert(solviolbounds != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) <= SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Real* varvals;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->exprgraph != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &varvals, SCIPexprgraphGetNVars(conshdlrdata->exprgraph)) );
      SCIP_CALL( SCIPgetSolVals(scip, sol, SCIPexprgraphGetNVars(conshdlrdata->exprgraph), (SCIP_VAR**)SCIPexprgraphGetVars(conshdlrdata->exprgraph), varvals) );

      SCIP_CALL( SCIPexprgraphEval(conshdlrdata->exprgraph, varvals) );

      SCIPfreeBufferArray(scip, &varvals);
   }

   *maxviolcon = NULL;

   maxviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      SCIP_CALL( computeViolation(scip, conshdlr, conss[c], sol, solviolbounds) );

      /* stop if solution violates bounds */
      if( *solviolbounds )
         break;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      viol = MAX(consdata->lhsviol, consdata->rhsviol);
      if( viol > maxviol && SCIPisGT(scip, viol, SCIPfeastol(scip)) )
      {
         maxviol = viol;
         *maxviolcon = conss[c];
      }

      /* SCIPdebugMsg(scip, "constraint <%s> violated by (%g, %g), activity = %g\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->activity); */
   }

   return SCIP_OKAY;
}

/** adds linearization of a constraints expression tree in reference point to a row */
static
SCIP_RETCODE addLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a linearization should be added */
   SCIP_Real*            x,                  /**< value of expression tree variables where to generate cut */
   SCIP_Bool             newx,               /**< whether the last evaluation of the expression with the expression interpreter was not at x */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to add linearization */
   SCIP_Bool*            success             /**< buffer to store whether a linearization was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real treecoef;
   SCIP_Real val;
   SCIP_Real* grad;
   SCIP_Real constant = 0.0;
   SCIP_Bool perturbedx;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x    != NULL);
   assert(rowprep  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);
   assert(newx || SCIPexprtreeGetInterpreterData(exprtree) != NULL);

   treecoef = consdata->nonlincoefs[exprtreeidx];

   *success = FALSE;

   /* compile expression if evaluated the first time; can only happen if newx is FALSE */
   if( newx && SCIPexprtreeGetInterpreterData(exprtree) == NULL )
   {
      SCIP_CALL( SCIPexprintCompile(exprint, exprtree) );
   }

   nvars = SCIPexprtreeGetNVars(exprtree);
   SCIP_CALL( SCIPallocBufferArray(scip, &grad, nvars) );

   perturbedx = FALSE;
   do
   {
      /* get value and gradient */
      SCIP_CALL( SCIPexprintGrad(exprint, exprtree, x, newx, &val, grad) );
      if( SCIPisFinite(val) && !SCIPisInfinity(scip, REALABS(val)) )
      {
         val *= treecoef;
         /* check gradient entries and compute constant f(refx) - grad * refx */
         constant = val;
         for( i = 0; i < nvars; ++i )
         {
            if( !SCIPisFinite(grad[i]) || SCIPisInfinity(scip, grad[i]) || SCIPisInfinity(scip, -grad[i]) )
               break;

            grad[i] *= treecoef;
            constant -= grad[i] * x[i];

            /* try to perturb x if the constant is too large */
            if( SCIPisInfinity(scip, REALABS(constant)) )
               break;

            /* coefficients smaller than epsilon are rounded to 0.0 when added to row, this can be wrong if variable value is very large (bad numerics)
             * in this case, set gradient to 0.0 here, but modify constant so that cut is still valid (if possible)
             * i.e., estimate grad[i]*x >= grad[i] * bound(x) or grad[i]*x <= grad[i] * bound(x), depending on whether we compute an underestimator (convex) or an overestimator (concave)
             * if required bound of x is not finite, then give up
             */
            if( grad[i] != 0.0 && SCIPisZero(scip, grad[i]) )
            {
               SCIP_VAR* var;
               SCIP_Real xbnd;

               var = SCIPexprtreeGetVars(exprtree)[i];
               if( consdata->curvatures[exprtreeidx] & SCIP_EXPRCURV_CONVEX )
               {
                  xbnd = grad[i] > 0.0 ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var);
               }
               else
               {
                  assert(consdata->curvatures[exprtreeidx] & SCIP_EXPRCURV_CONCAVE);
                  xbnd = grad[i] > 0.0 ? SCIPvarGetUbGlobal(var) : SCIPvarGetLbGlobal(var);
               }
               if( !SCIPisInfinity(scip, REALABS(xbnd)) )
               {
                  SCIPdebugMsg(scip, "var <%s> [%g,%g] has tiny gradient %g, replace coefficient by constant %g\n",
                     SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), grad[i], grad[i] * xbnd);
                  constant += grad[i] * xbnd;
                  grad[i] = 0.0;
               }
               else
               {
                  *success = FALSE;
                  SCIPdebugMsg(scip, "skipping linearization, var <%s> [%g,%g] has tiny gradient %g but no finite bound in this direction\n",
                     SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), grad[i]);
                  SCIPfreeBufferArray(scip, &grad);
                  return SCIP_OKAY;
               }
            }
         }

         if( i == nvars )
            break;
      }

      SCIPdebugMsg(scip, "got nonfinite value in evaluation or gradient of <%s>: ", SCIPconsGetName(cons));
      if( !perturbedx )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         SCIPdebugMsgPrint(scip, "perturbing reference point and trying again\n");
         for( i = 0; i < nvars; ++i )
         {
            lb = SCIPvarGetLbGlobal(SCIPexprtreeGetVars(exprtree)[i]);
            ub = SCIPvarGetUbGlobal(SCIPexprtreeGetVars(exprtree)[i]);
            if( SCIPisEQ(scip, x[i], lb) )
               x[i] += MIN(0.9*(ub-lb), i*SCIPfeastol(scip));  /*lint !e666*/
            else if( SCIPisEQ(scip, x[i], ub) )
               x[i] -= MIN(0.9*(ub-lb), i*SCIPfeastol(scip));  /*lint !e666*/
            else
               x[i] += MIN3(0.9*(ub-x[i]), 0.9*(x[i]-lb), i*SCIPfeastol(scip)) * (i%2 != 0 ? -1.0 : 1.0);  /*lint !e666*/
         }
         newx = TRUE;
         perturbedx = TRUE;
      }
      else
      {
         SCIPdebugMsgPrint(scip, "skipping linearization\n");
         SCIPfreeBufferArray(scip, &grad);
         return SCIP_OKAY;
      }
   }
   while( TRUE );  /*lint !e506*/

   /* add linearization to SCIP row */
   SCIPaddRowprepConstant(rowprep, constant);
   SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, nvars, SCIPexprtreeGetVars(exprtree), grad) );

   *success = TRUE;

   SCIPfreeBufferArray(scip, &grad);

   SCIPdebugMsg(scip, "added linearization for tree %d of constraint <%s>\n", exprtreeidx, SCIPconsGetName(cons));
   SCIPdebug( SCIPprintRowprep(scip, rowprep, NULL) );

   return SCIP_OKAY;
}

/** adds secant of a constraints univariate expression tree in reference point to a row */
static
SCIP_RETCODE addConcaveEstimatorUnivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to add estimator */
   SCIP_Bool*            success             /**< buffer to store whether a secant was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real      treecoef;
   SCIP_VAR*      var;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      vallb;
   SCIP_Real      valub;
   SCIP_Real      slope;
   SCIP_Real      constant;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(rowprep  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);
   assert(SCIPexprtreeGetNVars(exprtree) == 1);

   treecoef = consdata->nonlincoefs[exprtreeidx];

   *success = FALSE;

   var = SCIPexprtreeGetVars(exprtree)[0];
   xlb = SCIPvarGetLbLocal(var);
   xub = SCIPvarGetUbLocal(var);

   /* if variable is unbounded, then cannot really compute secant */
   if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
   {
      SCIPdebugMsg(scip, "skip secant for tree %d of constraint <%s> since variable is unbounded\n", exprtreeidx, SCIPconsGetName(cons));
      return SCIP_OKAY;
   }
   assert(SCIPisLE(scip, xlb, xub));

   SCIP_CALL( SCIPexprtreeEval(exprtree, &xlb, &vallb) ); /* coverity ignore ARRAY_VS_SINGLETON warning */
   if( !SCIPisFinite(vallb) || SCIPisInfinity(scip, REALABS(vallb)) )
   {
      SCIPdebugMsg(scip, "skip secant for tree %d of constraint <%s> since function cannot be evaluated in lower bound\n", exprtreeidx, SCIPconsGetName(cons));
      return SCIP_OKAY;
   }
   vallb *= treecoef;

   SCIP_CALL( SCIPexprtreeEval(exprtree, &xub, &valub) ); /* coverity ignore ARRAY_VS_SINGLETON warning */
   if( !SCIPisFinite(valub) || SCIPisInfinity(scip, REALABS(valub)) )
   {
      SCIPdebugMsg(scip, "skip secant for tree %d of constraint <%s> since function cannot be evaluated in upper bound\n", exprtreeidx, SCIPconsGetName(cons));
      return SCIP_OKAY;
   }
   valub *= treecoef;

   if( SCIPisEQ(scip, xlb, xub) )
   {
      slope = 0.0;
      /* choose most conservative value for the cut */
      if( rowprep->sidetype == SCIP_SIDETYPE_LEFT )
         constant = MAX(vallb, valub);
      else
         constant = MIN(vallb, valub);
   }
   else
   {
      slope = (valub - vallb) / (xub - xlb);
      constant = vallb - slope * xlb;
   }

   /* add secant to SCIP row */
   SCIPaddRowprepConstant(rowprep, constant);
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, slope) );

   *success = TRUE;

   SCIPdebugMsg(scip, "added secant for tree %d of constraint <%s>, slope = %g\n", exprtreeidx, SCIPconsGetName(cons), slope);
   SCIPdebug( SCIPprintRowprep(scip, rowprep, NULL) );

   return SCIP_OKAY;
}

/** adds estimator of a constraints bivariate expression tree to a row
 * a reference point is given to decide which hyperplane to choose
 */
static
SCIP_RETCODE addConcaveEstimatorBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_Real*            ref,                /**< reference values of expression tree variables where to generate cut */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to add estimator */
   SCIP_Bool*            success             /**< buffer to store whether a secant was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real      treecoef;
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      ylb;
   SCIP_Real      yub;

   SCIP_Real      coefx;
   SCIP_Real      coefy;
   SCIP_Real      constant;

   SCIP_Real      p1[2];
   SCIP_Real      p2[2];
   SCIP_Real      p3[2];
   SCIP_Real      p4[2];
   SCIP_Real      p1val, p2val, p3val, p4val;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);
   assert(SCIPexprtreeGetNVars(exprtree) == 2);

   treecoef = consdata->nonlincoefs[exprtreeidx];

   *success = FALSE;

   x = SCIPexprtreeGetVars(exprtree)[0];
   y = SCIPexprtreeGetVars(exprtree)[1];
   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);
   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) || SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub) )
   {
      SCIPdebugMsg(scip, "skip bivariate secant since <%s> or <%s> is unbounded\n", SCIPvarGetName(x), SCIPvarGetName(y));
      return SCIP_OKAY;
   }

   /* reference point should not be outside of bounds */
   assert(SCIPisFeasLE(scip, xlb, ref[0]));
   assert(SCIPisFeasGE(scip, xub, ref[0]));
   ref[0] = MIN(xub, MAX(xlb, ref[0]));
   assert(SCIPisFeasLE(scip, ylb, ref[1]));
   assert(SCIPisFeasGE(scip, yub, ref[1]));
   ref[1] = MIN(yub, MAX(ylb, ref[1]));

   /* lower left */
   p1[0] = xlb;
   p1[1] = ylb;

   /* lower right */
   p2[0] = xub;
   p2[1] = ylb;

   /* upper right */
   p3[0] = xub;
   p3[1] = yub;

   /* upper left */
   p4[0] = xlb;
   p4[1] = yub;

   if( SCIPisEQ(scip, xlb, xub) && SCIPisEQ(scip, ylb, yub) )
   {
      SCIP_CALL( SCIPexprtreeEval(exprtree, p1, &p1val) );

      if( !SCIPisFinite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) )
      {
         SCIPdebugMsg(scip, "skip secant for tree %d of constraint <%s> since function cannot be evaluated\n", exprtreeidx, SCIPconsGetName(cons));
         return SCIP_OKAY;
      }

      p1val *= treecoef;

      coefx = 0.0;
      coefy = 0.0;
      constant = p1val;
   }
   else if( SCIPisEQ(scip, xlb, xub) )
   {
      /* secant between p1 and p4: p1val + [(p4val - p1val) / (yub - ylb)] * (y - ylb) */
      assert(!SCIPisEQ(scip, ylb, yub));

      SCIP_CALL( SCIPexprtreeEval(exprtree, p1, &p1val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p4, &p4val) );
      if( !SCIPisFinite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !SCIPisFinite(p4val) || SCIPisInfinity(scip, REALABS(p4val)) )
      {
         SCIPdebugMsg(scip, "skip secant for tree %d of constraint <%s> since function cannot be evaluated\n", exprtreeidx, SCIPconsGetName(cons));
         return SCIP_OKAY;
      }
      p1val *= treecoef;
      p4val *= treecoef;

      coefx = 0.0;
      coefy = (p4val - p1val) / (yub - ylb);
      constant = p1val - coefy * ylb;
   }
   else if( SCIPisEQ(scip, ylb, yub) )
   {
      /* secant between p1 and p2: p1val + [(p2val - p1val) / (xub - xlb)] * (x - xlb) */
      assert(!SCIPisEQ(scip, xlb, xub));

      SCIP_CALL( SCIPexprtreeEval(exprtree, p1, &p1val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p2, &p2val) );
      if( !SCIPisFinite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !SCIPisFinite(p2val) || SCIPisInfinity(scip, REALABS(p2val)) )
      {
         SCIPdebugMsg(scip, "skip secant for tree %d of constraint <%s> since function cannot be evaluated\n", exprtreeidx, SCIPconsGetName(cons));
         return SCIP_OKAY;
      }

      p1val *= treecoef;
      p2val *= treecoef;

      coefx = (p2val - p1val) / (xub - xlb);
      coefy = 0.0;
      constant = p1val - coefx * xlb;
   }
   else
   {
      SCIP_Real alpha, beta, gamma_, delta;
      SCIP_Bool tryother;
      SCIP_Bool doover;

      /* if function is convex, then we want an overestimator, otherwise we want an underestimator */
      assert(consdata->curvatures[exprtreeidx] == SCIP_EXPRCURV_CONVEX || consdata->curvatures[exprtreeidx] == SCIP_EXPRCURV_CONCAVE);
      doover = (consdata->curvatures[exprtreeidx] & SCIP_EXPRCURV_CONVEX);  /*lint !e641*/

      SCIP_CALL( SCIPexprtreeEval(exprtree, p1, &p1val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p2, &p2val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p3, &p3val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p4, &p4val) );
      if( !SCIPisFinite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !SCIPisFinite(p2val) || SCIPisInfinity(scip, REALABS(p2val)) ||
         ! SCIPisFinite(p3val) || SCIPisInfinity(scip, REALABS(p3val)) || !SCIPisFinite(p4val) || SCIPisInfinity(scip, REALABS(p4val)) )
      {
         SCIPdebugMsg(scip, "skip secant for tree %d of constraint <%s> since function cannot be evaluated\n", exprtreeidx, SCIPconsGetName(cons));
         return SCIP_OKAY;
      }
      p1val *= treecoef;
      p2val *= treecoef;
      p3val *= treecoef;
      p4val *= treecoef;

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

      /* Compute coefficients alpha, beta, gamma (>0), delta such that
       *   alpha*x + beta*y + gamma*z = delta
       * is satisfied by at least three of the corner points (p1,f(p1)), ..., (p4,f(p4)) and
       * the fourth corner point lies below this hyperplane.
       * Since we assume that f is convex, we then know that all points (x,y,f(x,y)) are below this hyperplane, i.e.,
       *    alpha*x + beta*y - delta <= -gamma * f(x,y),
       * or, equivalently,
       *   -alpha/gamma*x - beta/gamma*y + delta/gamma >= f(x,y).
       */

      tryother = FALSE;
      if( ref[1] <= ylb + (yub - ylb)/(xub - xlb) * (ref[0] - xlb) )
      {
         SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1], p3val,
               &alpha, &beta, &gamma_, &delta) );

         assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
         assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
         assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));

         /* if hyperplane through p1,p2,p3 does not overestimate f(p4), then it must be the other variant */
         if( alpha * p4[0] + beta * p4[1] + gamma_ * p4val > delta )
            tryother = TRUE;
         else if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
            (      !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta /gamma_)) )
         {
            /* if numerically bad, take alternative hyperplane */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1],
                  p4val, &alpha, &beta, &gamma_, &delta) );

            assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
            assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
            assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

            /* if hyperplane through p1,p3,p4 does not overestimate f(p2), then it must be the other variant */
            if( alpha * p2[0] + beta * p2[1] + gamma_ * p2val > delta )
               tryother = TRUE;
         }
      }
      else
      {
         SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );

         assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
         assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
         assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

         /* if hyperplane through p1,p3,p4 does not overestimate f(p2), then it must be the other variant */
         if( alpha * p2[0] + beta * p2[1] + gamma_ * p2val > delta )
            tryother = TRUE;
         else if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
            (      !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta /gamma_)) )
         {
            /* if numerically bad, take alternative */
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1],
                  p3val, &alpha, &beta, &gamma_, &delta) );

            assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
            assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
            assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));

            /* if hyperplane through p1,p2,p3 does not overestimate f(p4), then it must be the other variant */
            if( alpha * p4[0] + beta * p4[1] + gamma_ * p4val > delta )
               tryother = TRUE;
         }
      }

      if( tryother )
      {
         if( ref[1] <= yub + (ylb - yub)/(xub - xlb) * (ref[0] - xlb) )
         {
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1],
                  p4val, &alpha, &beta, &gamma_, &delta) );

            /* hyperplane should be above (p3,f(p3)) and other points should lie on hyperplane */
            assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
            assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
            assert(SCIPisRelLE(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
            assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

            if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
               ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta /gamma_)) )
            {
               /* if numerically bad, take alternative */
               SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1],
                     p4val, &alpha, &beta, &gamma_, &delta) );

               /* hyperplane should be above (p1,f(p1)) and other points should lie on hyperplane */
               assert(SCIPisRelLE(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
               assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
               assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
               assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));
            }
         }
         else
         {
            SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1],
                  p4val, &alpha, &beta, &gamma_, &delta) );

            /* hyperplane should be above (p1,f(p1)) and other points should lie on hyperplane */
            assert(SCIPisRelLE(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
            assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
            assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
            assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

            if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
               ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta /gamma_)) )
            {
               /* if numerically bad, take alternative */
               SCIP_CALL( SCIPcomputeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1],
                     p4val, &alpha, &beta, &gamma_, &delta) );

               /* hyperplane should be above (p3,f(p3)) and other points should lie on hyperplane */
               assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
               assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
               assert(SCIPisRelLE(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
               assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));
            }
         }
      }

      SCIPdebugMsg(scip, "alpha = %g, beta = %g, gamma = %g, delta = %g\n", alpha, beta, gamma_, delta);

      /* check if bad luck: should not happen if xlb != xub and ylb != yub and numerics are fine */
      if( SCIPisZero(scip, gamma_) )
         return SCIP_OKAY;
      assert(!SCIPisNegative(scip, gamma_));

      /* flip hyperplane */
      if( !doover )
         gamma_ = -gamma_;

      coefx    = -alpha / gamma_;
      coefy    = -beta  / gamma_;
      constant =  delta / gamma_;

      /* if we loose coefficients because division by gamma makes them < SCIPepsilon(scip), then better not generate a cut here */
      if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, coefx)) ||
         ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, coefy)) )
      {
         SCIPdebugMsg(scip, "skip bivar secant for <%s> tree %d due to bad numerics\n", SCIPconsGetName(cons), exprtreeidx);
         return SCIP_OKAY;
      }
   }

   /* add hyperplane coefs to SCIP row */
   SCIPaddRowprepConstant(rowprep, constant);
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, x, coefx) );
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, y, coefy) );

   *success = TRUE;

   SCIPdebugMsg(scip, "added bivariate secant for tree %d of constraint <%s>\n", exprtreeidx, SCIPconsGetName(cons));
   SCIPdebug( SCIPprintRowprep(scip, rowprep, NULL) );

   return SCIP_OKAY;
}

/** internal method using an auxiliary LPI, see addConcaveEstimatorMultivariate() */
static
SCIP_RETCODE _addConcaveEstimatorMultivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LPI*             lpi,                /**< auxiliary LPI */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_Real*            ref,                /**< reference values of expression tree variables where to generate cut */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to add estimator */
   SCIP_VAR**            vars,               /**< variables of the constraint */
   SCIP_EXPRTREE*        exprtree,           /**< expression tree of constraint */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             doupper,            /**< should an upper estimator be computed */
   SCIP_Bool*            success             /**< buffer to store whether a secant was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* val;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* corner;
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   int* beg;
   int* ind;
   SCIP_Real lpobj;
   int ncols;
   int nrows;
   int nnonz;
   SCIP_Real funcval;
   SCIP_Real treecoef;

   int i;
   int j;
   int idx;

   SCIP_RETCODE lpret;

   assert(lpi != NULL);
   assert(nvars <= 10);

   consdata = SCIPconsGetData(cons);
   treecoef = consdata->nonlincoefs[exprtreeidx];

   /* columns are cut coefficients plus constant */
   ncols = nvars + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   corner = lb; /* will not use lb and corner simultaneously, so can share memory */

   /* one row for each corner of domain, i.e., 2^nvars many */
   nrows = (int)(1u << nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nrows) );

   /* the coefficients matrix will have at most ncols * nrows many nonzeros */
   nnonz = nrows * ncols;
   SCIP_CALL( SCIPallocBufferArray(scip, &beg, nrows+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind, nnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnonz) );

   /* setup LP data */
   idx = 0;
   for( i = 0; i < nrows; ++i )
   {
      /* assemble corner point */
      SCIPdebugMsg(scip, "f(");
      for( j = 0; j < nvars; ++j )
      {
         /* if j'th bit of row index i is set, then take upper bound on var j, otherwise lower bound var j
          * we check this by shifting i for j positions to the right and checking whether the j'th bit is set */
         if( ((unsigned int)i >> j) & 0x1 )
            corner[j] = SCIPvarGetUbLocal(vars[j]);
         else
            corner[j] = SCIPvarGetLbLocal(vars[j]);
         SCIPdebugMsgPrint(scip, "%g, ", corner[j]);
         assert(!SCIPisInfinity(scip, REALABS(corner[j])));
      }

      /* evaluate function in current corner */
      SCIP_CALL( SCIPexprtreeEval(exprtree, corner, &funcval) );
      SCIPdebugMsgPrint(scip, ") = %g\n", funcval);

      if( !SCIPisFinite(funcval) || SCIPisInfinity(scip, REALABS(funcval)) )
      {
         SCIPdebugMsg(scip, "cannot compute underestimator for concave because constaint <%s> cannot be evaluated\n", SCIPconsGetName(cons));
         goto TERMINATE;
      }

      funcval *= treecoef;

      if( !doupper )
      {
         lhs[i] = -SCIPlpiInfinity(lpi);
         rhs[i] = funcval;
      }
      else
      {
         lhs[i] = funcval;
         rhs[i] = SCIPlpiInfinity(lpi);
      }

      /* add nonzeros of corner to matrix */
      beg[i] = idx;
      for( j = 0; j < nvars; ++j )
      {
         if( corner[j] != 0.0 )
         {
            ind[idx] = j;
            val[idx] = corner[j];
            ++idx;
         }
      }

      /* coefficient for constant is 1.0 */
      val[idx] = 1.0;
      ind[idx] = nvars;
      ++idx;
   }
   nnonz = idx;
   beg[nrows] = nnonz;

   for( j = 0; j < ncols; ++j )
   {
      lb[j] = -SCIPlpiInfinity(lpi);
      ub[j] =  SCIPlpiInfinity(lpi);
   }

   /* objective coefficients are reference points, and an additional 1.0 for the constant */
   BMScopyMemoryArray(obj, ref, nvars);
   obj[nvars] = 1.0;

   /* get function value in reference point, so we can use this as a cutoff */
   SCIP_CALL( SCIPexprtreeEval(exprtree, ref, &funcval) );
   funcval *= treecoef;

   SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddRows(lpi, nrows, lhs, rhs, NULL, nnonz, beg, ind, val) );

   /* make use of this convenient features, since for us nrows >> ncols */
   /*SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_ROWREPSWITCH, 5.0) ); */
   /* get accurate coefficients */
   SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_FEASTOL, SCIPfeastol(scip)/100.0) );
   SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, funcval) );
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, 10 * nvars) );
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_SCALING, 1) );
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FROMSCRATCH, 1) );

   /* SCIPdebug( SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1) ) ); */

   lpret = SCIPlpiSolveDual(lpi);
   if( lpret != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "solving auxiliary LP for underestimator of concave function returned %d\n", lpret);
      goto TERMINATE;
   }

   if( !SCIPlpiIsPrimalFeasible(lpi) )
   {
      SCIPdebugMsg(scip, "failed to find feasible solution for auxiliary LP for underestimator of concave function, iterlimexc = %u, cutoff = %u, unbounded = %u\n", SCIPlpiIsIterlimExc(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsPrimalUnbounded(lpi));
      goto TERMINATE;
   }
   /* should be either solved to optimality, or the objective or iteration limit be hit */
   assert(SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsIterlimExc(lpi));

   /* setup row coefficient, reuse obj array to store LP sol values */
   SCIP_CALL( SCIPlpiGetSol(lpi, &lpobj, obj, NULL, NULL, NULL) );

   /* check that computed hyperplane is on right side of function in refpoint
    * if numerics is very bad (e.g., st_e32), then even this can happen */
   if( (!doupper && SCIPisFeasGT(scip, lpobj, funcval)) || (doupper && SCIPisFeasGT(scip, funcval, lpobj)) )
   {
      SCIPwarningMessage(scip, "computed cut does not underestimate concave function in refpoint\n");
      goto TERMINATE;
   }
   assert( doupper || SCIPisFeasLE(scip, lpobj, funcval) );
   assert(!doupper || SCIPisFeasLE(scip, funcval, lpobj) );

   /* add estimator to rowprep */
   SCIPaddRowprepConstant(rowprep, obj[nvars]);
   SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, nvars, vars, obj) );

   *success = TRUE;

TERMINATE:
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &beg);
   SCIPfreeBufferArray(scip, &ind);
   SCIPfreeBufferArray(scip, &val);

   return SCIP_OKAY;
}



/** adds estimator of a constraints multivariate expression tree to a row
 * Given concave function f(x) and reference point ref.
 * Let (v_i: i=1,...,n) be corner points of current domain of x.
 * Find (coef,constant) such that <coef,v_i> + constant <= f(v_i) (cut validity) and
 * such that <coef, ref> + constant is maximized (cut efficacy).
 * Then <coef, x> + constant <= f(x) for all x in current domain.
 *
 * Similar to compute an overestimator for a convex function f(x).
 * Find (coef,constant) such that <coef,v_i> + constant >= f(v_i) and
 * such that <coef, ref> + constant is minimized.
 * Then <coef, x> + constant >= f(x) for all x in current domain.
 */
static
SCIP_RETCODE addConcaveEstimatorMultivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_Real*            ref,                /**< reference values of expression tree variables where to generate cut */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to add estimator */
   SCIP_Bool*            success             /**< buffer to store whether a secant was succefully added to the row */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_LPI* lpi;
   int nvars;
   int j;
   SCIP_Bool doupper;

   SCIP_RETCODE retcode;

   static SCIP_Bool warned_highdim_concave = FALSE;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);

   nvars = SCIPexprtreeGetNVars(exprtree);
   assert(nvars >= 2);

   *success = FALSE;

   /* size of LP is exponential in number of variables of tree, so do only for small trees */
   if( nvars > 10 )
   {
      if( !warned_highdim_concave )
      {
         SCIPwarningMessage(scip, "concave function in constraint <%s> too high-dimensional to compute underestimator\n", SCIPconsGetName(cons));
         warned_highdim_concave = TRUE;
      }
      return SCIP_OKAY;
   }

   vars = SCIPexprtreeGetVars(exprtree);

   /* check whether bounds are finite
    * make sure reference point is strictly within bounds
    * otherwise we can easily get an unbounded LP below, e.g., with instances like ex6_2_* from GlobalLib
    */
   for( j = 0; j < nvars; ++j )
   {
      if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(vars[j])) || SCIPisInfinity(scip, SCIPvarGetUbLocal(vars[j])) )
      {
         SCIPdebugMsg(scip, "cannot compute underestimator for concave because variable <%s> is unbounded\n", SCIPvarGetName(vars[j]));
         return SCIP_OKAY;
      }
      assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(vars[j]), ref[j]));
      assert(SCIPisFeasGE(scip, SCIPvarGetUbLocal(vars[j]), ref[j]));
      ref[j] = MIN(SCIPvarGetUbLocal(vars[j]), MAX(SCIPvarGetLbLocal(vars[j]), ref[j]));  /*lint !e666*/
   }

   /* create empty auxiliary LP and decide its objective sense */
   assert(consdata->curvatures[exprtreeidx] == SCIP_EXPRCURV_CONVEX || consdata->curvatures[exprtreeidx] == SCIP_EXPRCURV_CONCAVE);
   doupper = (consdata->curvatures[exprtreeidx] & SCIP_EXPRCURV_CONVEX);  /*lint !e641*/
   SCIP_CALL( SCIPlpiCreate(&lpi, SCIPgetMessagehdlr(scip), "concaveunderest", doupper ? SCIP_OBJSEN_MINIMIZE : SCIP_OBJSEN_MAXIMIZE) );
   if( lpi == NULL )
   {
      SCIPerrorMessage("failed to create auxiliary LP\n");
      return SCIP_ERROR;
   }

   /* capture the retcode, free the LPI afterwards */
   retcode = _addConcaveEstimatorMultivariate(scip, lpi, cons, exprtreeidx, ref, rowprep, vars, exprtree, nvars, doupper, success);

   assert(lpi != NULL);
   SCIP_CALL( SCIPlpiFree(&lpi) );

   return retcode;
}

/** Computes the linear coeffs and the constant in a linear expression
 * both scaled by a given scalar value.
 * The coeffs of the variables will be stored in the given array at
 * their variable index.
 * The constant of the given linear expression will be added to the given
 * buffer.
 */
static
SCIP_RETCODE getCoeffsAndConstantFromLinearExpr(
   SCIP_EXPR*            expr,               /**< the linear expression */
   SCIP_Real             scalar,             /**< the scalar value, i.e. the coeff of the given expression */
   SCIP_Real*            varcoeffs,          /**< buffer array to store the computed coefficients */
   SCIP_Real*            constant            /**< buffer to hold the constant value of the given expression */
   )
{
   switch( SCIPexprGetOperator( expr ) )
   {
   case SCIP_EXPR_VARIDX: /* set coeff for this variable to current scalar */
   {
      /* TODO: can a linear expression contain the same variable twice?
       * if yes varcoeffs need to be initialized to zero before calling this function
       * and coeff must not be overridden but summed up instead. */
      varcoeffs[SCIPexprGetOpIndex( expr )] = scalar;
      return SCIP_OKAY;
   }

   case SCIP_EXPR_CONST:
   {
      /* constant value increases */
      *constant += scalar * SCIPexprGetOpReal( expr );
      return SCIP_OKAY;
   }

   case SCIP_EXPR_MUL: /* need to find the constant part of the muliplication and then recurse  */
   {
      SCIP_EXPR** children;
      children = SCIPexprGetChildren( expr );

      /* first part is constant */
      if( SCIPexprGetOperator( children[0] ) == SCIP_EXPR_CONST )
      {
         SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[1], scalar * SCIPexprGetOpReal( children[0] ), varcoeffs, constant ) );
         return SCIP_OKAY;
      }

      /* second part is constant */
      if( SCIPexprGetOperator( children[1] ) == SCIP_EXPR_CONST )
      {
         SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[0], scalar * SCIPexprGetOpReal( children[1] ), varcoeffs, constant ) );
         return SCIP_OKAY;
      }

      /* nonlinear -> break out to error case  */
      break;
   }

   case SCIP_EXPR_PLUS: /* just recurse */
   {
      SCIP_EXPR** children;
      children = SCIPexprGetChildren( expr );
      SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[0], scalar, varcoeffs, constant ) );
      SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[1], scalar, varcoeffs, constant ) );
      return SCIP_OKAY;
   }

   case SCIP_EXPR_MINUS: /* recursion on second child is called with negated scalar */
   {
      SCIP_EXPR** children;
      children = SCIPexprGetChildren( expr );
      SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[0], scalar, varcoeffs, constant ) );
      SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[1], -scalar, varcoeffs, constant ) );
      return SCIP_OKAY;
   }

   case SCIP_EXPR_SUM: /* just recurse */
   {
      SCIP_EXPR** children;
      int nchildren;
      int c;

      children = SCIPexprGetChildren(expr);
      nchildren = SCIPexprGetNChildren(expr);

      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[c], scalar, varcoeffs, constant ) );
      }

      return SCIP_OKAY;
   }

   case SCIP_EXPR_LINEAR: /* add scaled constant and recurse on children with their coeff multiplied into scalar */
   {
      SCIP_Real* childCoeffs;
      SCIP_EXPR** children;
      int i;

      *constant += scalar * SCIPexprGetLinearConstant( expr );

      children = SCIPexprGetChildren( expr );
      childCoeffs = SCIPexprGetLinearCoefs( expr );

      for( i = 0; i < SCIPexprGetNChildren( expr ); ++i )
      {
         SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[i], scalar * childCoeffs[i], varcoeffs, constant ) );
      }

      return SCIP_OKAY;
   }

   default:
      break;
   } /*lint !e788*/

   SCIPerrorMessage( "Cannot extract linear coefficients from expressions with operator %d %s\n", SCIPexprGetOperator( expr ), SCIPexpropGetName(SCIPexprGetOperator( expr )));
   SCIPABORT();
   return SCIP_ERROR; /*lint !e527*/
}

/** adds estimator from user callback of a constraints user expression tree to a row
 */
static
SCIP_RETCODE addUserEstimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree an estimator should be added */
   SCIP_Real*            x,                  /**< value of expression tree variables where to generate cut */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to add estimator */
   SCIP_Bool*            success             /**< buffer to store whether an estimator was succefully added to the rowprep */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR** children;
   SCIP_VAR** vars;
   SCIP_Real* params;
   SCIP_INTERVAL* varbounds;

   SCIP_INTERVAL* childbounds;
   SCIP_Real* childvals;
   SCIP_Real* childcoeffs;

   SCIP_Real constant;
   SCIP_Real treecoef;
   int nvars;
   int nchildren;
   int i;

   consdata = SCIPconsGetData( cons );
   assert( consdata != NULL );
   assert( exprtreeidx >= 0 );
   assert( exprtreeidx < consdata->nexprtrees );
   assert( consdata->exprtrees != NULL );
   assert( rowprep != NULL );
   assert( success != NULL );

   exprtree = consdata->exprtrees[exprtreeidx];
   assert( exprtree != NULL );
   assert( SCIPexprGetOperator(SCIPexprtreeGetRoot(exprtree)) == SCIP_EXPR_USER );

   /* if user did not implement estimator callback, then we cannot do anything */
   if( !SCIPexprHasUserEstimator(SCIPexprtreeGetRoot(exprtree)) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   params = SCIPexprtreeGetParamVals( exprtree );
   nvars = SCIPexprtreeGetNVars( exprtree );
   vars = SCIPexprtreeGetVars( exprtree );
   nchildren = SCIPexprGetNChildren( SCIPexprtreeGetRoot( exprtree ) );
   children = SCIPexprGetChildren( SCIPexprtreeGetRoot( exprtree ) );

   /* Get bounds of variables */
   SCIP_CALL( SCIPallocBufferArray( scip, &varbounds, nchildren ) );

   for( i = 0; i < nvars; ++i )
   {
      double lb = SCIPvarGetLbLocal( vars[i] );
      double ub = SCIPvarGetUbLocal( vars[i] );
      SCIPintervalSetBounds( &varbounds[i],
                             -infty2infty( SCIPinfinity( scip ), INTERVALINFTY, -MIN( lb, ub ) ),
                             +infty2infty( SCIPinfinity( scip ), INTERVALINFTY,  MAX( lb, ub ) ) );
   }

   /* Compute bounds and solution value for the user expressions children */
   SCIP_CALL( SCIPallocBufferArray( scip, &childcoeffs, nchildren ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &childbounds, nchildren ) );
   SCIP_CALL( SCIPallocBufferArray( scip, &childvals, nchildren ) );

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( SCIPexprEval( children[i], x, params, &childvals[i] ) );
      SCIP_CALL( SCIPexprEvalInt( children[i], INTERVALINFTY, varbounds, params, &childbounds[i] ) );
   }

   /* varbounds not needed any longer */
   SCIPfreeBufferArray( scip, &varbounds );

   /* call estimator for user expressions to compute coeffs and constant for the user expressions children */
   SCIP_CALL( SCIPexprEstimateUser( SCIPexprtreeGetRoot( exprtree ), INTERVALINFTY, childvals, childbounds, overestimate, childcoeffs, &constant, success ) );

   if( *success )
   {
      SCIP_Real* varcoeffs;
      SCIP_CALL( SCIPallocBufferArray( scip, &varcoeffs, nvars ) );

      treecoef = consdata->nonlincoefs[exprtreeidx];
      constant *= treecoef;

      for( i = 0; i < nchildren; ++i )
      {
         SCIP_CALL( getCoeffsAndConstantFromLinearExpr( children[i], childcoeffs[i]*treecoef, varcoeffs, &constant ) );
      }

      SCIPaddRowprepConstant(rowprep, constant);
      SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, nvars, vars, varcoeffs) );

      SCIPfreeBufferArray( scip, &varcoeffs );
   }

   SCIPfreeBufferArray( scip, &childcoeffs );
   SCIPfreeBufferArray( scip, &childbounds );
   SCIPfreeBufferArray( scip, &childvals );

   return SCIP_OKAY;
}

/** adds estimator from interval gradient of a constraints univariate expression tree to a row
 * a reference point is used to decide in which corner to generate the cut
 */
static
SCIP_RETCODE addIntervalGradientEstimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_Real*            x,                  /**< value of expression tree variables where to generate cut */
   SCIP_Bool             newx,               /**< whether the last evaluation of the expression with the expression interpreter was not at x */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_ROWPREP*         rowprep,            /**< rowprep where to add estimator */
   SCIP_Bool*            success             /**< buffer to store whether an estimator was succefully added to the rowprep */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real treecoef;
   SCIP_Real* coefs;
   SCIP_Real constant;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_INTERVAL* box;
   SCIP_INTERVAL* intgrad;
   SCIP_INTERVAL intval;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x    != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);
   assert(newx || SCIPexprtreeGetInterpreterData(exprtree) != NULL);

   *success = FALSE;

   /* skip interval gradient if expression interpreter cannot compute interval gradients */
   if( !(SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_INTGRADIENT) )
      return SCIP_OKAY;

   nvars = SCIPexprtreeGetNVars(exprtree);
   vars = SCIPexprtreeGetVars(exprtree);

   box = NULL;
   intgrad = NULL;
   coefs = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &box, nvars) );

   /* move reference point to bounds, setup box */
   for( i = 0; i < nvars; ++i )
   {
      lb  = SCIPvarGetLbLocal(vars[i]);
      ub  = SCIPvarGetUbLocal(vars[i]);
      if( SCIPisInfinity(scip, -lb) )
      {
         if( SCIPisInfinity(scip, ub) )
         {
            SCIPdebugMsg(scip, "skip interval gradient estimator for constraint <%s> because variable <%s> is still unbounded.\n", SCIPconsGetName(cons), SCIPvarGetName(vars[i]));
            goto INTGRADESTIMATOR_CLEANUP;
         }
         x[i] = ub;
      }
      else
      {
         if( SCIPisInfinity(scip, ub) )
            x[i] = lb;
         else
            x[i] = (2.0*x[i] < lb+ub) ? lb : ub;
      }
      SCIPintervalSetBounds(&box[i],
         -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -MIN(lb, ub)),
         +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  MAX(lb, ub)));
   }

   /* compile expression if evaluated the first time; can only happen if newx is FALSE */
   if( newx && SCIPexprtreeGetInterpreterData(exprtree) == NULL )
   {
      SCIP_CALL( SCIPexprintCompile(exprint, exprtree) );
   }

   /* evaluate in reference point */
   SCIP_CALL( SCIPexprintEval(exprint, exprtree, x, &val) );
   if( !SCIPisFinite(val) )
   {
      SCIPdebugMsg(scip, "Got nonfinite function value from evaluation of constraint %s tree %d. skipping interval gradient estimator.\n", SCIPconsGetName(cons), exprtreeidx);
      goto INTGRADESTIMATOR_CLEANUP;
   }

   treecoef = consdata->nonlincoefs[exprtreeidx];
   val *= treecoef;
   constant = val;

   /* compute interval gradient */
   SCIP_CALL( SCIPallocBufferArray(scip, &intgrad, nvars) );
   SCIP_CALL( SCIPexprintGradInt(exprint, exprtree, INTERVALINFTY, box, TRUE, &intval, intgrad) );
   SCIPintervalMulScalar(INTERVALINFTY, &intval, intval, treecoef);

   /* printf("nvars %d side %d xref = %g x = [%g,%g] intval = [%g,%g] intgrad = [%g,%g]\n", nvars, side, x[0],
      box[0].inf, box[0].sup, intval.inf, intval.sup, intgrad[0].inf, intgrad[0].sup); */

   /* compute coefficients and constant */
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      val = x[i];
      lb  = SCIPintervalGetInf(box[i]);
      ub  = SCIPintervalGetSup(box[i]);

      SCIPintervalMulScalar(INTERVALINFTY, &intgrad[i], intgrad[i], treecoef);

      if( SCIPisEQ(scip, lb, ub) )
         coefs[i] = 0.0;
      else if( (overestimate && val == ub) ||  /*lint !e777*/
         (!overestimate && val == lb) )   /*lint !e777*/
         coefs[i] = SCIPintervalGetInf(intgrad[i]);
      else
         coefs[i] = SCIPintervalGetSup(intgrad[i]);

      if( SCIPisZero(scip, coefs[i]) )
         continue;

      if( SCIPisInfinity(scip, -coefs[i]) || SCIPisInfinity(scip, coefs[i]) )
      {
         SCIPdebugMsg(scip, "skip intgrad estimator because of infinite interval bound\n");
         goto INTGRADESTIMATOR_CLEANUP;
      }

      constant -= coefs[i] * val;
   }

   /* add interval gradient estimator to row */
   SCIPaddRowprepConstant(rowprep, constant);
   SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, nvars, vars, coefs) );

 INTGRADESTIMATOR_CLEANUP:
   SCIPfreeBufferArrayNull(scip, &box);
   SCIPfreeBufferArrayNull(scip, &intgrad);
   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** generates a cut based on linearization (if convex), secant (if concave), or intervalgradient (if indefinite)
 */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real**           ref,                /**< reference point for each exprtree, or NULL if sol should be used */
   SCIP_SOL*             sol,                /**< reference solution where cut should be generated, or NULL if LP solution should be used */
   SCIP_Bool             newsol,             /**< whether the last evaluation of the expression with the expression interpreter was not at sol */
   SCIP_SIDETYPE         side,               /**< for which side a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real             minviol,            /**< minimal absolute violation we try to achieve */
   SCIP_Real             maxrange,           /**< maximal range allowed */
   SCIP_Bool             expensivecurvchecks,/**< whether also expensive checks should be executed */
   SCIP_Bool             assumeconvex        /**< whether to assume convexity in inequalities */
   )
{
   SCIP_ROWPREP* rowprep;
   SCIP_CONSDATA* consdata;
   SCIP_Bool success;
   SCIP_Real* x;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   SCIPdebugMsg(scip, "constructing cut for %s hand side of constraint <%s>\n", side == SCIP_SIDETYPE_LEFT ? "left" : "right", SCIPconsGetName(cons));

   SCIP_CALL( checkCurvature(scip, cons, expensivecurvchecks, assumeconvex) );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nexprtrees == 0 )
   {
      char rowname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%u", SCIPconsGetName(cons), ++(consdata->ncuts));

      /* if we are actually linear, add the constraint as row to the LP */
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, row, SCIPconsGetHdlr(cons), rowname, consdata->lhs, consdata->rhs, SCIPconsIsLocal(cons), FALSE , TRUE) );
      SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, side,
      !(side == SCIP_SIDETYPE_LEFT  && (consdata->curvature & SCIP_EXPRCURV_CONCAVE)) &&
      !(side == SCIP_SIDETYPE_RIGHT && (consdata->curvature & SCIP_EXPRCURV_CONVEX ))) );
   SCIPaddRowprepSide(rowprep, side == SCIP_SIDETYPE_LEFT  ? consdata->lhs : consdata->rhs);
   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%s_%u", SCIPconsGetName(cons), ++(consdata->ncuts));

   if( ref == NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &x, SCIPexprtreeGetNVars(consdata->exprtrees[0])) );
   }

   success = TRUE;
   for( i = 0; i < consdata->nexprtrees; ++i )
   {
      if( ref == NULL )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &x, SCIPexprtreeGetNVars(consdata->exprtrees[i])) );  /*lint !e644*/
         SCIP_CALL( SCIPgetSolVals(scip, sol, SCIPexprtreeGetNVars(consdata->exprtrees[i]), SCIPexprtreeGetVars(consdata->exprtrees[i]), x) );
      }
      else
      {
         x = ref[i];
      }

      if( (side == SCIP_SIDETYPE_LEFT && (consdata->curvatures[i] & SCIP_EXPRCURV_CONCAVE)) ||
         (side == SCIP_SIDETYPE_RIGHT && (consdata->curvatures[i] & SCIP_EXPRCURV_CONVEX )) )
      {
         SCIP_CALL( addLinearization(scip, exprint, cons, i, x, newsol, rowprep, &success) );
      }
      else if( (side == SCIP_SIDETYPE_LEFT  && (consdata->curvatures[i] & SCIP_EXPRCURV_CONVEX)) ||
         (      side == SCIP_SIDETYPE_RIGHT && (consdata->curvatures[i] & SCIP_EXPRCURV_CONCAVE)) )
      {
         switch( SCIPexprtreeGetNVars(consdata->exprtrees[i]) )
         {
         case 1:
            SCIP_CALL( addConcaveEstimatorUnivariate(scip, cons, i, rowprep, &success) );
            break;

         case 2:
            SCIP_CALL( addConcaveEstimatorBivariate(scip, cons, i, x, rowprep, &success) );
            break;

         default:
            SCIP_CALL( addConcaveEstimatorMultivariate(scip, cons, i, x, rowprep, &success) );
            break;
         }
         if( !success )
         {
            SCIPdebugMsg(scip, "failed to generate polyhedral estimator for %d-dim concave function in exprtree %d, fall back to intervalgradient cut\n", SCIPexprtreeGetNVars(consdata->exprtrees[i]), i);
            SCIP_CALL( addIntervalGradientEstimator(scip, exprint, cons, i, x, newsol, side == SCIP_SIDETYPE_LEFT, rowprep, &success) );
         }
      }
      else if( SCIPexprGetOperator( SCIPexprtreeGetRoot( consdata->exprtrees[i] ) ) == SCIP_EXPR_USER )
      {
         SCIP_CALL( addUserEstimator( scip, cons, i, x, side == SCIP_SIDETYPE_LEFT, rowprep, &success ) );
         if( !success ) /* the user estimation may not be implemented -> try interval estimator */
         {
            SCIP_CALL( addIntervalGradientEstimator(scip, exprint, cons, i, x, newsol, side == SCIP_SIDETYPE_LEFT, rowprep, &success) );
         }
      }
      else
      {
         SCIP_CALL( addIntervalGradientEstimator(scip, exprint, cons, i, x, newsol, side == SCIP_SIDETYPE_LEFT, rowprep, &success) );
      }

      if( !success )
         break;
   }

   if( ref == NULL )
   {
      SCIPfreeBufferArray(scip, &x);
   }

   if( success )
   {
      SCIP_Real coefrange;

      /* add coefficients for linear variables */
      SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );

      /* merge terms in same variable */
      SCIPmergeRowprepTerms(scip, rowprep);

      /* cleanup row */
      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, maxrange, minviol, &coefrange, NULL) );

      /* check that coefficient range is ok */
      success = coefrange <= maxrange;

      /* check that side is finite */ /*lint --e{514} */
      success &= !SCIPisInfinity(scip, REALABS(rowprep->side));

      /* check whether maximal coef is finite, if any */
      success &= (rowprep->nvars == 0) || !SCIPisInfinity(scip, REALABS(rowprep->coefs[0]));
   }

   if( success )
   {
      SCIP_CALL( SCIPgetRowprepRowCons(scip, row, rowprep, SCIPconsGetHdlr(cons)) );
   }
   else
      *row = NULL;

   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** tries to separate solution or LP solution by a linear cut
 *
 *  assumes that constraint violations have been computed
 */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Bool             newsol,             /**< have the constraints just been evaluated at this point? */
   SCIP_Real             minefficacy,        /**< minimal efficacy of a cut if it should be added to the LP */
   SCIP_Bool             inenforcement,      /**< whether we are in constraint enforcement */
   SCIP_RESULT*          result,             /**< result of separation */
   SCIP_Real*            bestefficacy        /**< buffer to store best efficacy of a cut that was added to the LP, if found; or NULL if not of interest */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          efficacy;
   SCIP_Real          feasibility;
   SCIP_SIDETYPE      violside;
   int                c;
   SCIP_ROW*          row;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( bestefficacy != NULL )
      *bestefficacy = 0.0;

   for( c = 0, violside = SCIP_SIDETYPE_LEFT; c < nconss; c = (violside == SCIP_SIDETYPE_RIGHT ? c+1 : c), violside = (violside == SCIP_SIDETYPE_LEFT ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT) )
   {
      assert(conss != NULL);

      /* skip constraints that are not enabled */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* if side violside of cons c not violated, then continue to next side or next cons */
      if( !SCIPisGT(scip, violside == SCIP_SIDETYPE_LEFT ? consdata->lhsviol : consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      /* we are not feasible anymore */
      if( *result == SCIP_FEASIBLE )
         *result = SCIP_DIDNOTFIND;

      /* generate cut
       * if function is defined at sol (activity<infinity) and constraint is violated, then expression interpreter should have evaluated at sol to get gradient before
       */
      SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], NULL, sol, newsol || SCIPisInfinity(scip, consdata->activity), violside, &row, minefficacy, conshdlrdata->cutmaxrange, conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );

      if( row == NULL ) /* failed to generate cut */
         continue;

      if( sol == NULL )
         feasibility = SCIPgetRowLPFeasibility(scip, row);
      else
         feasibility = SCIPgetRowSolFeasibility(scip, row, sol);
      efficacy = -feasibility;

      if( SCIPisGT(scip, efficacy, minefficacy) && SCIPisCutApplicable(scip, row) )
      {
         SCIP_Bool infeasible;

         /* cut cuts off solution */
         SCIP_CALL( SCIPaddRow(scip, row, FALSE /* forcecut */, &infeasible) );
         if ( infeasible )
            *result = SCIP_CUTOFF;
         else
            *result = SCIP_SEPARATED;

         SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );

         SCIPdebugMsg(scip, "add cut with efficacy %g for constraint <%s> violated by %g\n", efficacy, SCIPconsGetName(conss[c]), MAX(consdata->lhsviol, consdata->rhsviol));
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

         if( bestefficacy != NULL && efficacy > *bestefficacy )
            *bestefficacy = efficacy;

         /* mark row as not removable from LP for current node, if in enforcement */
         if( inenforcement && !conshdlrdata->enfocutsremovable )
            SCIPmarkRowNotRemovableLocal(scip, row);
      }
      else
      {
         SCIPdebugMsg(scip, "drop cut since efficacy %g is too small (< %g)\n", efficacy, minefficacy);
      }

      SCIP_CALL( SCIPreleaseRow (scip, &row) );

      if ( *result == SCIP_CUTOFF )
         break;

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */
      if( c >= nusefulconss && *result == SCIP_SEPARATED )
         break;
   }

   return SCIP_OKAY;
}

/** adds linearizations cuts for convex constraints w.r.t. a given reference point to cutpool and sepastore
 * if separatedlpsol is not NULL, then a cut that separates the LP solution is added to the sepastore and is forced to enter the LP
 * if separatedlpsol is not NULL, but cut does not separate the LP solution, then it is added to the cutpool only
 * if separatedlpsol is NULL, then cut is added to cutpool only
 */
static
SCIP_RETCODE addLinearizationCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             ref,                /**< reference point where to linearize, or NULL for LP solution */
   SCIP_Bool*            separatedlpsol,     /**< buffer to store whether a cut that separates the current LP solution was found and added to LP, or NULL if adding to cutpool only */
   SCIP_Real             minefficacy         /**< minimal efficacy of a cut when checking for separation of LP solution */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool addedtolp;
   SCIP_ROW* row;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( separatedlpsol != NULL )
      *separatedlpsol = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      if( SCIPconsIsLocal(conss[c]) )  /*lint !e613*/
         continue;

      SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* if we cannot linearize, then skip constraint */
      if( (!(consdata->curvature & SCIP_EXPRCURV_CONVEX)  || SCIPisInfinity(scip,  consdata->rhs)) &&
         ( !(consdata->curvature & SCIP_EXPRCURV_CONCAVE) || SCIPisInfinity(scip, -consdata->lhs)) )
         continue;

      SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], NULL, ref, TRUE,
            (consdata->curvature & SCIP_EXPRCURV_CONVEX) ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT,
            &row, minefficacy, conshdlrdata->cutmaxrange, FALSE, FALSE) );  /*lint !e613*/

      if( row == NULL )
         continue;

      addedtolp = FALSE;

      /* if caller wants, then check if cut separates LP solution and add to sepastore if so */
      if( separatedlpsol != NULL )
      {
         if( -SCIPgetRowLPFeasibility(scip, row) >= minefficacy )
         {
            SCIP_Bool infeasible;

            *separatedlpsol = TRUE;
            addedtolp = TRUE;
            SCIP_CALL( SCIPaddRow(scip, row, TRUE, &infeasible) );
            assert( ! infeasible );
         }
      }

      if( !SCIProwIsLocal(row) && !addedtolp )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** processes the event that a new primal solution has been found */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS**    conss;
   int            nconss;
   SCIP_SOL*      sol;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) != 0);

   conshdlr = (SCIP_CONSHDLR*)eventdata;

   nconss = SCIPconshdlrGetNConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* we are only interested in solution coming from some heuristic other than trysol, but not from the tree
    * the reason for ignoring trysol solutions is that they may come from an NLP solve in sepalp, where we already added linearizations,
    * or are from the tree, but postprocessed via proposeFeasibleSolution
    */
   if( SCIPsolGetHeur(sol) == NULL || SCIPsolGetHeur(sol) == conshdlrdata->trysolheur )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   SCIPdebugMsg(scip, "catched new sol event %" SCIP_EVENTTYPE_FORMAT " from heur <%s>; have %d conss\n", SCIPeventGetType(event), SCIPheurGetName(SCIPsolGetHeur(sol)), nconss);

   SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, sol, NULL, 0.0) );

   return SCIP_OKAY;
}

/** registers unfixed variables in nonlinear terms of violated constraints as external branching candidates */
static
SCIP_RETCODE registerBranchingVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int c;
   int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->nexprtrees == 0 )
         continue;

      /* do not branch on violation of convex constraint */
      if( (!SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || (consdata->curvature & SCIP_EXPRCURV_CONCAVE)) &&
         ( !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) || (consdata->curvature & SCIP_EXPRCURV_CONVEX )) )
         continue;
      SCIPdebugMsg(scip, "cons <%s> violation: %g %g  curvature: %s\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, SCIPexprcurvGetName(consdata->curvature));

      for( i = 0; i < consdata->nexprtrees; ++i )
      {
         /* skip convex summands */
         if( (!SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || (consdata->curvatures[i] & SCIP_EXPRCURV_CONCAVE)) &&
            ( !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) || (consdata->curvatures[i] & SCIP_EXPRCURV_CONVEX )) )
            continue;

         for( j = 0; j < SCIPexprtreeGetNVars(consdata->exprtrees[i]); ++j )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[i])[j];
            assert(var != NULL);

            if( SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            {
               SCIPdebugMsg(scip, "ignore fixed variable <%s>[%g, %g], width %g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var));
               continue;
            }

            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++*nnotify;
         }
      }
   }

   SCIPdebugMsg(scip, "registered %d branching candidates\n", *nnotify);

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
   int i;
   int j;
   int c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   *brvar = NULL;
   brvarval = -1.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]); ++i )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[j])[i];
            /* do not propose fixed variables */
            if( SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
               continue;
            val = SCIPgetSolVal(scip, sol, var);
            if( REALABS(val) > brvarval )
            {
               brvarval = ABS(val);
               *brvar = var;
            }
         }
      }
   }

   if( *brvar != NULL )
   {
      SCIP_CALL( SCIPaddExternBranchCand(scip, *brvar, brvarval, SCIP_INVALID) );
   }

   return SCIP_OKAY;
}

/** replaces violated nonlinear constraints where all nonlinear variables are almost fixed by linear constraints
 * only adds constraint if it is violated in current solution
 * first tries to fix almost fixed variables
 */
static
SCIP_RETCODE replaceViolatedByLinearConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            addedcons,          /**< buffer to store whether a linear constraint was added */
   SCIP_Bool*            reduceddom,         /**< whether a domain has been reduced */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_CONS*          cons;
   SCIP_CONSDATA*      consdata;
   SCIP_Real           lhs;
   SCIP_Real           rhs;
   SCIP_Real           lb;
   SCIP_Real           ub;
   SCIP_RESULT         checkresult;
   SCIP_VAR*           var;
   SCIP_Bool           tightened;
   int c;
   int i;
   int v;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);
   assert(addedcons != NULL);
   assert(reduceddom != NULL);
   assert(infeasible != NULL);

   *addedcons = FALSE;
   *reduceddom = FALSE;
   *infeasible = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      lhs = consdata->lhs;
      rhs = consdata->rhs;

      for( i = 0; i < consdata->nexprtrees; ++i )
      {
         SCIP_INTERVAL nonlinactivity;

         /* check whether there are almost fixed nonlinear variables that can be fixed */
         for( v = 0; v < SCIPexprtreeGetNVars(consdata->exprtrees[i]); ++v )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[i])[v];

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);
            assert(SCIPisRelEQ(scip, lb, ub)); /* variable should be almost fixed */

            assert(!SCIPisInfinity(scip, -lb));
            assert(!SCIPisInfinity(scip,  ub));

            if( !SCIPisEQ(scip, lb, ub) )
            {
               /* try to fix variable */
               SCIP_CALL( SCIPtightenVarLb(scip, var, (lb+ub)/2.0, TRUE, infeasible, &tightened) );
               if( *infeasible )
               {
                  SCIPdebugMsg(scip, "Fixing almost fixed variable <%s> lead to infeasibility.\n", SCIPvarGetName(var));
                  return SCIP_OKAY;
               }
               if( tightened )
               {
                  SCIPdebugMsg(scip, "Tightened lower bound of almost fixed variable <%s>.\n", SCIPvarGetName(var));
                  *reduceddom = TRUE;
               }

               SCIP_CALL( SCIPtightenVarUb(scip, var, (lb+ub)/2.0, TRUE, infeasible, &tightened) );
               if( *infeasible )
               {
                  SCIPdebugMsg(scip, "Fixing almost fixed variable <%s> lead to infeasibility.\n", SCIPvarGetName(var));
                  return SCIP_OKAY;
               }
               if( tightened )
               {
                  SCIPdebugMsg(scip, "Tightened upper bound of almost fixed variable <%s>.\n", SCIPvarGetName(var));
                  *reduceddom = TRUE;
               }
            }
         }

         SCIP_CALL( SCIPevalExprtreeLocalBounds(scip, consdata->exprtrees[i], INTERVALINFTY, &nonlinactivity) );
         SCIPintervalMulScalar(INTERVALINFTY, &nonlinactivity, nonlinactivity, consdata->nonlincoefs[i]);

         if( !SCIPisInfinity(scip, -lhs) )
         {
            if( SCIPintervalGetSup(nonlinactivity) >= INTERVALINFTY )
               lhs = -SCIPinfinity(scip);
            else if( SCIPintervalGetSup(nonlinactivity) <= -INTERVALINFTY )
            {
               /* lhs <= [...,-infinity] + ...  will never be feasible */
               *infeasible = TRUE;
               return SCIP_OKAY;
            }
            else
               lhs -= SCIPintervalGetSup(nonlinactivity);
         }

         if( !SCIPisInfinity(scip,  rhs) )
         {
            if( SCIPintervalGetInf(nonlinactivity) <= -INTERVALINFTY )
               rhs = SCIPinfinity(scip);
            else if( SCIPintervalGetInf(nonlinactivity) >= INTERVALINFTY )
            {
               /* [infinity,...] + ... <= rhs will never be feasible */
               *infeasible = TRUE;
               return SCIP_OKAY;
            }
            else
               rhs -= SCIPintervalGetInf(nonlinactivity);
         }
      }

      /* if some nonlinear variable was fixed now, then restart node (next enfo round) */
      if( *reduceddom )
         return SCIP_OKAY;

      /* check if we have a bound change */
      if ( consdata->nlinvars == 0 )
      {
         assert(SCIPisFeasLE(scip, lhs, rhs));
      }
      else if ( consdata->nlinvars == 1 )
      {
         SCIP_Real coef;

         coef = *consdata->lincoefs;
         SCIPdebugMsg(scip, "Linear constraint with one variable: %g <= %g <%s> <= %g\n", lhs, coef, SCIPvarGetName(*consdata->linvars), rhs);

         /* possibly correct lhs/rhs */
         assert( ! SCIPisZero(scip, coef) );
         if ( coef >= 0.0 )
         {
            if ( ! SCIPisInfinity(scip, -lhs) )
               lhs /= coef;
            if ( ! SCIPisInfinity(scip, rhs) )
               rhs /= coef;
         }
         else
         {
            SCIP_Real h;
            h = rhs;
            if ( ! SCIPisInfinity(scip, -lhs) )
               rhs = lhs/coef;
            else
               rhs = SCIPinfinity(scip);

            if ( ! SCIPisInfinity(scip, h) )
               lhs = h/coef;
            else
               lhs = -SCIPinfinity(scip);
         }
         SCIPdebugMsg(scip, "Linear constraint is a bound: %g <= <%s> <= %g\n", lhs, SCIPvarGetName(*consdata->linvars), rhs);

         /* cut off the node if SCIP needs to tight the lb/ub to +/-inf */
         if( SCIPisInfinity(scip, lhs) || SCIPisInfinity(scip, -rhs) )
         {
            *infeasible = TRUE;
            assert(consdata->linvars[0] != NULL);
            SCIPwarningMessage(scip, "Activity of nonlinear part is beyond SCIP's value for infinity. To enforce "
               "the constraint %s SCIP needs to tight bounds of %s to a value beyond +/- infinity. Please check if "
               "finite bounds can be added.\n", SCIPconsGetName(conss[c]), SCIPvarGetName(consdata->linvars[0]));
            return SCIP_OKAY;
         }

         if ( ! SCIPisInfinity(scip, -lhs) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, *consdata->linvars, lhs, TRUE, infeasible, &tightened) );
            if ( *infeasible )
            {
               SCIPdebugMsg(scip, "Lower bound leads to infeasibility.\n");
               return SCIP_OKAY;
            }
            if ( tightened )
            {
               SCIPdebugMsg(scip, "Lower bound changed.\n");
               *reduceddom = TRUE;
               return SCIP_OKAY;
            }
         }

         if ( ! SCIPisInfinity(scip, rhs) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, *consdata->linvars, rhs, TRUE, infeasible, &tightened) );
            if ( *infeasible )
            {
               SCIPdebugMsg(scip, "Upper bound leads to infeasibility.\n");
               return SCIP_OKAY;
            }
            if ( tightened )
            {
               SCIPdebugMsg(scip, "Upper bound changed.\n");
               *reduceddom = TRUE;
               return SCIP_OKAY;
            }
         }
      }
      else
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, SCIPconsGetName(conss[c]),
               consdata->nlinvars, consdata->linvars, consdata->lincoefs, lhs, rhs,
               SCIPconsIsInitial(conss[c]), SCIPconsIsSeparated(conss[c]), SCIPconsIsEnforced(conss[c]),
               SCIPconsIsChecked(conss[c]), SCIPconsIsPropagated(conss[c]),  TRUE,
               SCIPconsIsModifiable(conss[c]), SCIPconsIsDynamic(conss[c]), SCIPconsIsRemovable(conss[c]),
               SCIPconsIsStickingAtNode(conss[c])) );

         SCIPdebugMsg(scip, "replace violated nonlinear constraint <%s> by linear constraint after all nonlinear vars have been fixed\n", SCIPconsGetName(conss[c]) );
         SCIPdebugPrintCons(scip, conss[c], NULL);
         SCIPdebugPrintCons(scip, cons, NULL);

         SCIP_CALL( SCIPcheckCons(scip, cons, NULL, FALSE, FALSE, FALSE, &checkresult) );

         if( checkresult != SCIP_INFEASIBLE && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIPdebugMsg(scip, "linear constraint is feasible, thus do not add\n");
         }
         else
         {
            SCIP_CALL( SCIPaddConsLocal(scip, cons, NULL) );
            *addedcons = TRUE;
         }
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
      SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
   }

   return SCIP_OKAY;
}

/* tightens a lower bound on a variable and checks the result */
static
SCIP_RETCODE propagateBoundsTightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where we currently propagate, or NULL if tightening is from expression graph */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             bnd,                /**< new lower bound for variable */
   SCIP_RESULT*          result,             /**< result to update if there was a tightening or cutoff */
   int*                  nchgbds             /**< counter to increase if a bound was tightened */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(bnd > -INTERVALINFTY);
   assert(var != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   if( SCIPisInfinity(scip, bnd) )
   { /* domain will be outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      return SCIP_OKAY;
   }

   /* new lower bound is very low (between -INTERVALINFTY and -SCIPinfinity()) */
   if( SCIPisInfinity(scip, -bnd) )
      return SCIP_OKAY;

   bnd = SCIPadjustedVarLb(scip, var, bnd);
   SCIP_CALL( SCIPtightenVarLb(scip, var, bnd, FALSE, &infeas, &tightened) );
   if( infeas )
   {
      SCIPdebugMsg(scip, "%sfound constraint <%s> infeasible due to tightened lower bound %g for variable <%s>\n", SCIPinProbing(scip) ? "in probing " : "", cons != NULL ? SCIPconsGetName(cons) : "??", bnd, SCIPvarGetName(var));  /*lint !e585*/
      *result = SCIP_CUTOFF;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      return SCIP_OKAY;
   }
   if( tightened )
   {
      SCIPdebugMsg(scip, "%stightened lower bound of variable <%s> in constraint <%s> to %.20g\n", SCIPinProbing(scip) ? "in probing " : "", SCIPvarGetName(var), cons != NULL ? SCIPconsGetName(cons) : "??", bnd);  /*lint !e585*/
      ++*nchgbds;
      *result = SCIP_REDUCEDDOM;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/* tightens an upper bound on a variable and checks the result */
static
SCIP_RETCODE propagateBoundsTightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where we currently propagate, or NULL if tightening is from expression graph */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             bnd,                /**< new upper bound for variable */
   SCIP_RESULT*          result,             /**< result to update if there was a tightening or cutoff */
   int*                  nchgbds             /**< counter to increase if a bound was tightened */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(bnd < INTERVALINFTY);
   assert(var != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   if( SCIPisInfinity(scip, -bnd) )
   { /* domain will be outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      return SCIP_OKAY;
   }

   /* new upper bound is very high (between SCIPinfinity() and INTERVALINFTY) */
   if( SCIPisInfinity(scip, bnd) )
      return SCIP_OKAY;

   bnd = SCIPadjustedVarUb(scip, var, bnd);
   SCIP_CALL( SCIPtightenVarUb(scip, var, bnd, FALSE, &infeas, &tightened) );
   if( infeas )
   {
      SCIPdebugMsg(scip, "%sfound constraint <%s> infeasible due to tightened upper bound %g for variable <%s>\n", SCIPinProbing(scip) ? "in probing " : "", cons != NULL ? SCIPconsGetName(cons) : "??", bnd, SCIPvarGetName(var));  /*lint !e585*/
      *result = SCIP_CUTOFF;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      return SCIP_OKAY;
   }
   if( tightened )
   {
      SCIPdebugMsg(scip, "%stightened upper bound of variable <%s> in constraint <%s> to %g\n", SCIPinProbing(scip) ? "in probing " : "", SCIPvarGetName(var), cons != NULL ? SCIPconsGetName(cons) : "??", bnd);  /*lint !e585*/
      ++*nchgbds;
      *result = SCIP_REDUCEDDOM;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** tightens bounds of linear variables in a single nonlinear constraint */
static
SCIP_RETCODE propagateBoundsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation call */
   int*                  nchgbds,            /**< buffer where to add the the number of changed bounds */
   SCIP_Bool*            redundant           /**< buffer where to store whether constraint has been found to be redundant */
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA*     consdata;
   SCIP_INTERVAL      consbounds;    /* lower and upper bounds of constraint */
   SCIP_INTERVAL      consactivity;  /* activity of linear plus nonlinear part */
   SCIP_VAR*          var;
   SCIP_INTERVAL      rhs;           /* right hand side of nonlinear equation */
   SCIP_ROUNDMODE     roundmode;
   SCIP_Real          bnd;
   int                i;
   SCIP_INTERVAL      nonlinactivity;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTRUN;
   *redundant = FALSE;

   if( !SCIPconsIsMarkedPropagate(cons) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "start linear vars domain propagation for constraint <%s>\n", SCIPconsGetName(cons));

   /* unmark constraint for propagation */
   SCIP_CALL( SCIPunmarkConsPropagate(scip, cons) );

   /* make sure we have activity of linear term */
   consdataUpdateLinearActivity(scip, consdata);
   assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777*/
   assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777*/
   assert(consdata->minlinactivityinf >= 0);
   assert(consdata->maxlinactivityinf >= 0);
   assert(consdata->exprgraphnode != NULL || consdata->nexprtrees == 0);

   /* get activity of nonlinear part, should have been updated in propagateBounds */
   if( consdata->exprgraphnode != NULL )
   {
      nonlinactivity = SCIPexprgraphGetNodeBounds(consdata->exprgraphnode);
   }
   else
   {
      SCIPintervalSet(&nonlinactivity, 0.0);
   }
   assert(!SCIPintervalIsEmpty(INTERVALINFTY, nonlinactivity) );

   /* get activity of constraint function */
   SCIPintervalSetBounds(&consactivity, consdata->minlinactivityinf > 0 ? -INTERVALINFTY : consdata->minlinactivity, consdata->maxlinactivityinf > 0 ? INTERVALINFTY : consdata->maxlinactivity);
   SCIPintervalAdd(INTERVALINFTY, &consactivity, consactivity, nonlinactivity);

   /* check infeasibility */
   if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisGT(scip, consdata->lhs-SCIPfeastol(scip), SCIPintervalGetSup(consactivity))) ||
       (!SCIPisInfinity(scip,  consdata->rhs) && SCIPisLT(scip, consdata->rhs+SCIPfeastol(scip), SCIPintervalGetInf(consactivity))) )
   {
      SCIPdebugMsg(scip, "found constraint <%s> to be infeasible; sides: [%g, %g], activity: [%g, %g], infeas: %.20g\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity),
         MAX(consdata->lhs - SCIPintervalGetSup(consactivity), SCIPintervalGetInf(consactivity) - consdata->rhs));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIPintervalSetBounds(&consbounds,
      -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -consdata->lhs + SCIPepsilon(scip)),
      +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  consdata->rhs + SCIPepsilon(scip)));

   /* check redundancy */
   if( SCIPintervalIsSubsetEQ(INTERVALINFTY, consactivity, consbounds) )
   {
      SCIPdebugMsg(scip, "found constraint <%s> to be redundant: sides: [%g, %g], activity: [%g, %g]\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity));
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   /* propagate linear part in rhs = consbounds - nonlinactivity (use the one from consdata, since that includes infinities) */
   SCIPintervalSub(INTERVALINFTY, &rhs, consbounds, nonlinactivity);
   if( !SCIPintervalIsEntire(INTERVALINFTY, rhs) )
   {
      SCIP_Real coef;

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         coef = consdata->lincoefs[i];
         var  = consdata->linvars[i];

         /* skip fixed variables
          * @todo is that a good or a bad idea?
          *   we can't expect much more tightening, but may detect infeasiblity, but shouldn't the check on the constraints activity detect that?
          */
         if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            continue;

         if( coef > 0.0 )
         {
            if( SCIPintervalGetSup(rhs) < INTERVALINFTY )
            {
               assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777*/
               /* try to tighten the upper bound on var x */
               if( consdata->minlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
                  /* tighten upper bound on x to (rhs.sup - (minlinactivity - coef * xlb)) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = SCIPintervalGetSup(rhs);
                  bnd -= consdata->minlinactivity;
                  bnd += coef * SCIPvarGetLbLocal(var);
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->minlinactivityinf == 1 && SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  /* x was the variable that made the minimal linear activity equal -infinity, so
                   * we tighten upper bound on x to just (rhs.sup - minlinactivity) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = SCIPintervalGetSup(rhs);
                  bnd -= consdata->minlinactivity;
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the minimal activity is -infinity and x is not solely responsible for this */
            }

            if( SCIPintervalGetInf(rhs) > -INTERVALINFTY )
            {
               assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777*/
               /* try to tighten the lower bound on var x */
               if( consdata->maxlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  /* tighten lower bound on x to (rhs.inf - (maxlinactivity - coef * xub)) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = SCIPintervalGetInf(rhs);
                  bnd -= consdata->maxlinactivity;
                  bnd += coef * SCIPvarGetUbLocal(var);
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->maxlinactivityinf == 1 && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal infinity, so
                   * we tighten upper bound on x to just (rhs.inf - maxlinactivity) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = SCIPintervalGetInf(rhs);
                  bnd -= consdata->maxlinactivity;
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the maximal activity is +infinity and x is not solely responsible for this */
            }
         }
         else
         {
            assert(coef < 0.0 );
            if( SCIPintervalGetInf(rhs) > -INTERVALINFTY )
            {
               assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777*/
               /* try to tighten the upper bound on var x */
               if( consdata->maxlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetLbLocal(var)));
                  /* compute upper bound on x to (maxlinactivity - coef * xlb) - rhs.inf / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = consdata->maxlinactivity;
                  bnd += (-coef) * SCIPvarGetLbLocal(var);
                  bnd -= SCIPintervalGetInf(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->maxlinactivityinf == 1 && SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal infinity, so
                   * we tighten upper bound on x to just (maxlinactivity - rhs.inf) / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = consdata->maxlinactivity;
                  bnd -= SCIPintervalGetInf(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the maximal activity is infinity and x is not solely responsible for this */
            }

            if( SCIPintervalGetSup(rhs) < INTERVALINFTY )
            {
               assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777*/
               /* try to tighten the lower bound on var x */
               if( consdata->minlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  /* compute lower bound on x to (minlinactivity - coef * xub) - rhs.sup / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = consdata->minlinactivity;
                  bnd += (-coef) * SCIPvarGetUbLocal(var);
                  bnd -= SCIPintervalGetSup(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->minlinactivityinf == 1 && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal -infinity, so
                   * we tighten lower bound on x to just (minlinactivity - rhs.sup) / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = consdata->minlinactivity;
                  bnd -= SCIPintervalGetSup(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the minimal activity is -infinity and x is not solely responsible for this */
            }
         }
      }
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** propagate constraints sides minus linear activity into nonlinear variables */
static
SCIP_RETCODE propagateConstraintSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation calls */
   int*                  nchgbds             /**< buffer where to add the number of changed bounds */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int         nvars;
   SCIP_VAR**  vars;
   SCIP_EXPRGRAPHNODE** varnodes;
   SCIP_INTERVAL bounds;
   SCIP_Bool   cutoff;
   SCIP_ROUNDMODE roundmode;
   int         c;
   int         i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   SCIPdebugMsg(scip, "start backward propagation in expression graph\n");

#ifdef SCIP_OUTPUT
   {
      FILE* file;
      file = fopen("exprgraph_propconss1.dot", "w");
      if( file != NULL )
      {
         SCIP_CALL( SCIPexprgraphPrintDot(conshdlrdata->exprgraph, SCIPgetMessagehdlr(scip), file, NULL) );
         fclose(file);
      }
   }
#endif

   /* put constraint sides less linear activity into expression graph nodes
    * also add a [-feastol,feastol] range around constraint sides to cope with numerics */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->exprgraphnode == NULL )
         continue;

      /* skip (just) deleted or disabled constraints */
      if( SCIPconsIsDeleted(conss[c]) || !SCIPconsIsEnabled(conss[c]) )
         continue;

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      if( !SCIPisInfinity(scip, -consdata->lhs) && consdata->maxlinactivityinf == 0 )
         bounds.inf = consdata->lhs - consdata->maxlinactivity - SCIPfeastol(scip);
      else
         bounds.inf = -INTERVALINFTY;

      if( !SCIPisInfinity(scip,  consdata->rhs) && consdata->minlinactivityinf == 0 )
         bounds.sup = SCIPintervalNegateReal(consdata->minlinactivity - consdata->rhs - SCIPfeastol(scip));
      else
         bounds.sup =  INTERVALINFTY;

      SCIPintervalSetRoundingMode(roundmode);

      /* if we want the expression graph to propagate the bounds in any case, we set minstrength to a negative value */
      SCIPexprgraphTightenNodeBounds(conshdlrdata->exprgraph, consdata->exprgraphnode, bounds,
         consdata->forcebackprop ? -1.0 : BOUNDTIGHTENING_MINSTRENGTH, INTERVALINFTY, &cutoff);
      consdata->forcebackprop = FALSE; /* do this only once */

      if( cutoff )
      {
         SCIPdebugMsg(scip, "found constraint <%s> infeasible%s\n", SCIPconsGetName(conss[c]), SCIPinProbing(scip) ? " in probing" : "");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   /* compute bound tightenings for nonlinear variables */
   SCIPexprgraphPropagateNodeBounds(conshdlrdata->exprgraph, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, &cutoff);

#ifdef SCIP_OUTPUT
   {
      FILE* file;
      file = fopen("exprgraph_propconss2.dot", "w");
      if( file != NULL )
      {
         SCIP_CALL( SCIPexprgraphPrintDot(conshdlrdata->exprgraph, SCIPgetMessagehdlr(scip), file, NULL) );
         fclose(file);
      }
   }
#endif

   if( cutoff )
   {
      SCIPdebugMsg(scip, "backward propagation found problem infeasible%s\n", SCIPinProbing(scip) ? " in probing" : "");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* put tighter bounds into variables */
   nvars = SCIPexprgraphGetNVars(conshdlrdata->exprgraph);
   vars  = (SCIP_VAR**)SCIPexprgraphGetVars(conshdlrdata->exprgraph);
   varnodes = SCIPexprgraphGetVarNodes(conshdlrdata->exprgraph);

   /* put back new bounds into SCIP variables */
   for( i = 0; i < nvars && *result != SCIP_CUTOFF; ++i )
   {
      if( !SCIPisInfinity(scip, -SCIPintervalGetInf(SCIPexprgraphGetNodeBounds(varnodes[i]))) )
      {
         SCIP_CALL( propagateBoundsTightenVarLb(scip, NULL, vars[i], SCIPintervalGetInf(SCIPexprgraphGetNodeBounds(varnodes[i])), result, nchgbds) );
      }
      if( *result != SCIP_CUTOFF && !SCIPisInfinity(scip,  SCIPintervalGetSup(SCIPexprgraphGetNodeBounds(varnodes[i]))) )
      {
         SCIP_CALL( propagateBoundsTightenVarUb(scip, NULL, vars[i], SCIPintervalGetSup(SCIPexprgraphGetNodeBounds(varnodes[i])), result, nchgbds) );
      }
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
   SCIP_Bool             needclear,          /**< whether we may need to clear remainings from a previous backward propagation */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation calls */
   int*                  nchgbds,            /**< buffer where to add the the number of changed bounds */
   int*                  ndelconss           /**< buffer where to increase if a constraint was deleted (locally) due to redundancy */
   )
{
#ifndef NDEBUG
   SCIP_CONSDATA* consdata;
#endif
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT propresult;
   SCIP_Bool   domainerror;
   SCIP_Bool   redundant;
   int         roundnr;
   SCIP_Bool   success;
   int         c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   if( nconss == 0 || conshdlrdata->ispropagated )
   {
      *result = SCIP_DIDNOTRUN;
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
      SCIP_CALL( SCIPexprgraphPropagateVarBounds(conshdlrdata->exprgraph, INTERVALINFTY, (roundnr == 0) && needclear, &domainerror) );

#ifdef SCIP_OUTPUT
      {
         FILE* file;
         file = fopen("exprgraph_propvars.dot", "w");
         if( file != NULL )
         {
            SCIP_CALL( SCIPexprgraphPrintDot(conshdlrdata->exprgraph, SCIPgetMessagehdlr(scip), file, NULL) );
            fclose(file);
         }
      }
#endif

      if( domainerror )
      {
         SCIPdebugMsg(scip, "current bounds out of domain for some expression, do cutoff\n");
         *result = SCIP_CUTOFF;
         break;
      }

      /* check for redundancy and infeasibility of constraints, tighten bounds on linear variables */
      for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
      {
         assert(conss != NULL);
         if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) )
            continue;
         assert(SCIPconsIsActive(conss[c]));

#ifndef NDEBUG
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         assert(consdata->exprgraphnode == NULL || !SCIPintervalIsEmpty(INTERVALINFTY, SCIPexprgraphGetNodeBounds(consdata->exprgraphnode)));
#endif

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
      }

      /* propagate backward through expression graph */
      if( *result != SCIP_CUTOFF )
      {
         propresult = SCIP_DIDNOTFIND;
         SCIP_CALL( propagateConstraintSides(scip, conshdlr, conss, nconss, &propresult, nchgbds) );

         if( propresult != SCIP_DIDNOTFIND )
         {
            *result = propresult;
            success = TRUE;
         }
      }
   }
   while( success && *result != SCIP_CUTOFF && ++roundnr < conshdlrdata->maxproprounds );

   return SCIP_OKAY;
}

/* checks for a linear variable that can be increase or decreased without harming feasibility */
static
void consdataFindUnlockedLinearVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int i;
   int poslock;
   int neglock;

   consdata->linvar_maydecrease = -1;
   consdata->linvar_mayincrease = -1;

   /* check for a linear variable that can be increase or decreased without harming feasibility
    * setup lincoefsmin, lincoefsmax */
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      /* compute locks of i'th linear variable */
      assert(consdata->lincoefs[i] != 0.0);
      if( consdata->lincoefs[i] > 0.0 )
      {
         poslock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
         neglock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
      }
      else
      {
         poslock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
         neglock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
      }

      if( SCIPvarGetNLocksDown(consdata->linvars[i]) - neglock == 0 )
      {
         /* for a*x + q(y) \in [lhs, rhs], we can decrease x without harming other constraints */
         /* if we have already one candidate, then take the one where the loss in the objective function is less */
         if( (consdata->linvar_maydecrease < 0) ||
            (SCIPvarGetObj(consdata->linvars[consdata->linvar_maydecrease]) / consdata->lincoefs[consdata->linvar_maydecrease] > SCIPvarGetObj(consdata->linvars[i]) / consdata->lincoefs[i]) )
            consdata->linvar_maydecrease = i;
      }

      if( SCIPvarGetNLocksDown(consdata->linvars[i]) - poslock == 0 )
      {
         /* for a*x + q(y) \in [lhs, rhs], we can increase x without harm */
         /* if we have already one candidate, then take the one where the loss in the objective function is less */
         if( (consdata->linvar_mayincrease < 0) ||
            (SCIPvarGetObj(consdata->linvars[consdata->linvar_mayincrease]) / consdata->lincoefs[consdata->linvar_mayincrease] > SCIPvarGetObj(consdata->linvars[i]) / consdata->lincoefs[i]) )
            consdata->linvar_mayincrease = i;
      }
   }

#ifdef SCIP_DEBUG
   if( consdata->linvar_mayincrease >= 0 )
   {
      SCIPdebugMsg(scip, "may increase <%s> to become feasible\n", SCIPvarGetName(consdata->linvars[consdata->linvar_mayincrease]));
   }
   if( consdata->linvar_maydecrease >= 0 )
   {
      SCIPdebugMsg(scip, "may decrease <%s> to become feasible\n", SCIPvarGetName(consdata->linvars[consdata->linvar_maydecrease]));
   }
#endif
}

/** Given a solution where every nonlinear constraint is either feasible or can be made feasible by
 * moving a linear variable, construct the corresponding feasible solution and pass it to the trysol heuristic.
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
   SCIP_Bool solviolbounds;

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
         SCIP_CALL( computeViolation(scip, conshdlr, conss[c], newsol, &solviolbounds) );  /*lint !e613*/
         assert(!solviolbounds);
         viol = consdata->lhs - consdata->activity;
      }
      else if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conshdlr, conss[c], newsol, &solviolbounds) );  /*lint !e613*/
         assert(!solviolbounds);
         viol = consdata->rhs - consdata->activity;
      }
      else
         continue; /* constraint is satisfied */

      assert(viol != 0.0);
      if( consdata->linvar_mayincrease >= 0 &&
         ((  viol > 0.0 && consdata->lincoefs[consdata->linvar_mayincrease] > 0.0) ||
            (viol < 0.0 && consdata->lincoefs[consdata->linvar_mayincrease] < 0.0)) )
      {
         /* have variable where increasing makes the constraint less violated */
         var = consdata->linvars[consdata->linvar_mayincrease];
         /* compute how much we would like to increase var */
         delta = viol / consdata->lincoefs[consdata->linvar_mayincrease];
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
            viol -= consdata->lincoefs[consdata->linvar_mayincrease] * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      assert(viol != 0.0);
      if( consdata->linvar_maydecrease >= 0 &&
         ((  viol > 0.0 && consdata->lincoefs[consdata->linvar_maydecrease] < 0.0) ||
            (viol < 0.0 && consdata->lincoefs[consdata->linvar_maydecrease] > 0.0)) )
      {
         /* have variable where decreasing makes constraint less violated */
         var = consdata->linvars[consdata->linvar_maydecrease];
         /* compute how much we would like to decrease var */
         delta = viol / consdata->lincoefs[consdata->linvar_maydecrease];
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
            viol -= consdata->lincoefs[consdata->linvar_maydecrease] * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      /* still here... so probably we could not make constraint feasible due to variable bounds, thus give up */
      break;
   }

   /* if we have a solution that should satisfy all nonlinear constraints and has a better objective than the current upper bound,
    * then pass it to the trysol heuristic */
   if( c == nconss && (SCIPisInfinity(scip, SCIPgetUpperbound(scip)) || SCIPisSumLT(scip, SCIPgetSolTransObj(scip, newsol), SCIPgetUpperbound(scip))) )
   {
      SCIPdebugMsg(scip, "pass solution with objective value %g to trysol heuristic\n", SCIPgetSolTransObj(scip, newsol));

      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, newsol) );
      *success = TRUE;
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

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
   SCIP_Bool          solviolbounds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, sol, &solviolbounds, &maxviolcons) );

   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   if( solviolbounds )
   {
      /* if LP solution violates variable bounds, then this should be because a row was added that
       * introduced this variable newly to the LP, in which case it gets value 0.0; the row should
       * have been added to resolve an infeasibility, so solinfeasible should be TRUE
       * see also issue #627
       */
      assert(solinfeasible);
      /* however, if solinfeasible is actually not TRUE, then better cut off the node to avoid that SCIP
       * stops because infeasible cannot be resolved */ /*lint --e{774} */
      if( !solinfeasible )
         *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   consdata = SCIPconsGetData(maxviolcons);
   assert(consdata != NULL);

   maxviol = consdata->lhsviol + consdata->rhsviol;
   assert(SCIPisGT(scip, maxviol, SCIPfeastol(scip)));

   SCIPdebugMsg(scip, "enforcement with max violation %g in cons <%s> for %s solution\n", maxviol, SCIPconsGetName(maxviolcons),
         sol == NULL ? "LP" : "relaxation");

   /* we propagate and separate constraints only if they are active and enforcing by branching only does not seem much effective */
   assert(SCIPconsIsActive(maxviolcons));

   /* if we are above the 100'th enforcement round for this node, something is strange
    * (maybe the LP does not think that the cuts we add are violated, or we do ECP on a high-dimensional convex function)
    * in this case, check if some limit is hit or SCIP should stop for some other reason and terminate enforcement by creating a dummy node
    * (in optimized more, returning SCIP_INFEASIBLE in *result would be sufficient, but in debug mode this would give an assert in scip.c)
    * the reason to wait for 100 rounds is to avoid calls to SCIPisStopped in normal runs, which may be expensive
    * we only increment nenfolprounds until 101 to avoid an overflow
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

   /* run domain propagation */
   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, TRUE, &propresult, &dummy, &dummy) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we would like a cut that is efficient enough that it is not redundant in the LP (>lpfeastol)
    * however, we also don't want very weak cuts, so try to reach at least feastol (=lpfeastol by default, though)
    */
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, TRUE /* because computeviolation projects point onto box */, SCIPfeastol(scip), TRUE, &separateresult, &sepaefficacy) );
   if( separateresult == SCIP_CUTOFF )
   {
      SCIPdebugMsg(scip, "separation found cutoff\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }
   if( separateresult == SCIP_SEPARATED )
   {
      SCIPdebugMsg(scip, "separation succeeded (bestefficacy = %g, minefficacy = %g)\n", sepaefficacy, SCIPfeastol(scip));
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }

   /* we are not feasible, the whole node is not infeasible, and we cannot find a good cut
    * -> collect variables for branching
    */
   SCIPdebugMsg(scip, "separation failed (bestefficacy = %g < %g = minefficacy ); max viol: %g\n", sepaefficacy, SCIPfeastol(scip), maxviol);

   /* find branching candidates */
   SCIP_CALL( registerBranchingVariables(scip, conshdlr, conss, nconss, &nnotify) );

   if( nnotify == 0 && !solinfeasible && SCIPfeastol(scip) > SCIPlpfeastol(scip) )
   {
      /* fallback 1: we also have no branching candidates, so try to find a weak cut */
      SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, FALSE, SCIPlpfeastol(scip), TRUE, &separateresult, &sepaefficacy) );
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
         /* fallback 3: all nonlinear variables in all violated constraints seem to be fixed -> replace by linear constraints */
         SCIP_Bool addedcons;
         SCIP_Bool reduceddom;
         SCIP_Bool infeasible;

         SCIPdebugMsg(scip, "All nonlinear variables seem to be fixed. Replace remaining violated nonlinear constraints by linear constraints.\n");
         SCIP_CALL( replaceViolatedByLinearConstraints(scip, conss, nconss, &addedcons, &reduceddom, &infeasible) );
         /* if the linear constraints are actually feasible, then adding them and returning SCIP_CONSADDED confuses SCIP
          * when it enforces the new constraints again and nothing resolves the infeasiblity that we declare here thus,
          * we only add them if considered violated, and otherwise claim the solution is feasible (but print a
          * warning) */
         if ( infeasible )
            *result = SCIP_CUTOFF;
         else if ( addedcons )
            *result = SCIP_CONSADDED;
         else if ( reduceddom )
            *result = SCIP_REDUCEDDOM;
         else
         {
            *result = SCIP_FEASIBLE;
            SCIPwarningMessage(scip, "could not enforce feasibility by separating or branching; declaring solution with viol %g as feasible\n", maxviol);
            assert(!SCIPisInfinity(scip, maxviol));
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
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyNonlinear)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   /* assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0); */

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprinterpreter != NULL);
   assert(conshdlrdata->exprgraph != NULL);
   assert(SCIPexprgraphGetNVars(conshdlrdata->exprgraph) == 0);

   /* free expression graph */
   SCIP_CALL( SCIPexprgraphFree(&conshdlrdata->exprgraph) );

   /* free upgrade functions */
   for( i = 0; i < conshdlrdata->nnlconsupgrades; ++i )
   {
      assert(conshdlrdata->nlconsupgrades[i] != NULL);
      SCIPfreeBlockMemory(scip, &conshdlrdata->nlconsupgrades[i]);  /*lint !e866*/
   }
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->nlconsupgrades, conshdlrdata->nlconsupgradessize);

   /* free expressions interpreter */
   SCIP_CALL( SCIPexprintFree(&conshdlrdata->exprinterpreter) );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");

   /* reset counter, since we have a new problem */
   conshdlrdata->naddedreformconss = 0;

#ifdef SCIP_OUTPUT
   {
      FILE* file;
      file = fopen("exprgraph_init.dot", "w");
      if( file != NULL )
      {
         SCIP_CALL( SCIPexprgraphPrintDot(conshdlrdata->exprgraph, SCIPgetMessagehdlr(scip), file, NULL) );
         fclose(file);
      }
   }
#endif

   return SCIP_OKAY;
}  /*lint !e715*/

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   for( c = 0; c < nconss; ++c )
   {
      /* skip not yet active constraints */
      if( !SCIPconsIsActive(conss[c]) )  /*lint !e613*/
         continue;

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* forget expression trees */
      assert(consdata->nexprtrees == 0 || consdata->exprgraphnode != NULL);
      SCIP_CALL( consdataSetExprtrees(scip, consdata, 0, NULL, NULL, FALSE) );

      /* mark constraint for propagation */
      SCIP_CALL( SCIPmarkConsPropagate(scip, conss[c]) );  /*lint !e613*/
   }

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          havegraphchange;
   SCIP_Bool          havechange;
   SCIP_Bool          domainerror;
#ifndef NDEBUG
   int i;
   int j;
#endif
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   havegraphchange = FALSE;

   if( !conshdlrdata->isremovedfixings )
   {
      SCIP_CALL( removeFixedNonlinearVariables(scip, conshdlr) );
      assert(conshdlrdata->isremovedfixings);

      havegraphchange = TRUE;
   }

   /* if undefined expressions in exprgraph (very unlikely), we will hopefully recognize this during domain propagation later (if it involved an active constraint) */
   SCIP_CALL( SCIPexprgraphSimplify(conshdlrdata->exprgraph, SCIPgetMessagehdlr(scip), SCIPepsilon(scip), conshdlrdata->maxexpansionexponent, &havechange, &domainerror) );
   SCIPdebugMsg(scip, "expression graph simplifier found %schange, domain error = %u\n", havechange ? "" : "no ", domainerror);
   havegraphchange |= havechange;

   /* some of the methods below will not work if there was a domain error (#1148, point 3) */
   if( domainerror )
      return SCIP_OKAY;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);

      /* skip inactive constraints */
      if( !SCIPconsIsActive(conss[c]) )
         continue;
      assert(SCIPconsIsAdded(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->isremovedfixingslin )
      {
         SCIP_CALL( removeFixedLinearVariables(scip, conss[c]) );
      }

      if( !consdata->ispresolved || havegraphchange )
      {
         SCIP_Bool infeasible;
         SCIP_CALL( splitOffLinearPart(scip, conshdlr, conss[c], &infeasible) );

         /* the infeasibility should have been detected during presolve */
         assert(!infeasible);
      }

      SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );

      assert(consdata->isremovedfixingslin);
      assert(consdata->linvarsmerged);
#ifndef NDEBUG
      for( i = 0; i < consdata->nlinvars; ++i )
         assert(SCIPvarIsActive(consdata->linvars[i]));
#endif

      if( consdata->exprgraphnode != NULL )
      {
         /* get expression trees from expression graph */
         SCIP_EXPRTREE** exprtrees;
         SCIP_Real* coefs;
         int nexprtrees;
         int exprtreessize;

         exprtreessize = SCIPexprgraphGetSumTreesNSummands(consdata->exprgraphnode);

         SCIP_CALL( SCIPallocBufferArray(scip, &exprtrees, exprtreessize) );
         SCIP_CALL( SCIPallocBufferArray(scip, &coefs,     exprtreessize) );

         SCIP_CALL( SCIPexprgraphGetSumTrees(conshdlrdata->exprgraph, consdata->exprgraphnode,
               exprtreessize, &nexprtrees, exprtrees, coefs) );
         assert(nexprtrees > 0);

         SCIP_CALL( consdataSetExprtrees(scip, consdata, nexprtrees, exprtrees, coefs, FALSE) );

         SCIPfreeBufferArray(scip, &exprtrees);
         SCIPfreeBufferArray(scip, &coefs);

         assert(consdata->nexprtrees > 0 );
#ifndef NDEBUG
         for( j = 0; j < consdata->nexprtrees; ++j )
            for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]); ++i )
               assert(SCIPvarIsActive(SCIPexprtreeGetVars(consdata->exprtrees[j])[i]));
#endif

         /* tell SCIP that we have something nonlinear */
         SCIPenableNLP(scip);
      }
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* check for a linear variable that can be increase or decreased without harming feasibility */
      consdataFindUnlockedLinearVar(scip, consdata);

      /* setup lincoefsmin, lincoefsmax */
      consdata->lincoefsmin = SCIPinfinity(scip);
      consdata->lincoefsmax = 0.0;
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         consdata->lincoefsmin = MIN(consdata->lincoefsmin, REALABS(consdata->lincoefs[i]));  /*lint !e666*/
         consdata->lincoefsmax = MAX(consdata->lincoefsmax, REALABS(consdata->lincoefs[i]));  /*lint !e666*/
      }

      /* add nlrow respresentation to NLP, if NLP had been constructed */
      if( SCIPisNLPConstructed(scip) && SCIPconsIsEnabled(conss[c]) )
      {
         if( consdata->nlrow == NULL )
         {
            /* compute curvature for the nonlinear constraint if not done yet */
            SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );

            SCIP_CALL( createNlRow(scip, conss[c]) );
            assert(consdata->nlrow != NULL);
         }
         SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
      }
   }

   conshdlrdata->newsoleventfilterpos = -1;
   if( nconss != 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );
   }

   /* reset flags and counters */
   conshdlrdata->sepanlp = FALSE;
   conshdlrdata->lastenfonode = NULL;
   conshdlrdata->nenforounds = 0;

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
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
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* free nonlinear row representation */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }
   }

   return SCIP_OKAY;
}  /*lint !e715*/


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteNonlinear)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsActive(cons));
   assert(consdata != NULL);
   assert(SCIPconsGetData(cons) == *consdata);

   SCIPdebugMsg(scip, "consDelete for cons <%s>\n", SCIPconsGetName(cons));

   /* expression should have been removed from expression graph when constraint was deactivated */
   assert((*consdata)->exprgraphnode == NULL);

   SCIP_CALL( consdataFree(scip, consdata) );

   assert(*consdata == NULL);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransNonlinear)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int            i;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( consdataCreate(scip, &targetdata,
         sourcedata->lhs, sourcedata->rhs,
         sourcedata->nlinvars, sourcedata->linvars, sourcedata->lincoefs,
         sourcedata->nexprtrees, sourcedata->exprtrees, sourcedata->nonlincoefs,
         FALSE) );

   /* copy information on curvature, if known in original constraint */
   if( sourcedata->iscurvchecked && sourcedata->nexprtrees > 0 )
   {
      BMScopyMemoryArray(targetdata->curvatures, sourcedata->curvatures, sourcedata->nexprtrees);
      targetdata->curvature = sourcedata->curvature;
      targetdata->iscurvchecked = TRUE;
   }

   for( i = 0; i < targetdata->nlinvars; ++i )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->linvars[i], &targetdata->linvars[i]) );
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->linvars[i]) );
   }

   for( i = 0; i < targetdata->nexprtrees; ++i )
   {
      SCIP_CALL( SCIPgetExprtreeTransformedVars(scip, targetdata->exprtrees[i]) );
   }

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   SCIPdebugMsg(scip, "created transformed nonlinear constraint ");
   SCIPdebugPrintCons(scip, *targetcons, NULL);

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_ROW*          row;
   int                c;
   SCIP_Real**        x;
   int                nvars;
   int                i;
   int                j;
   SCIP_VAR*          var;
   SCIP_Bool          haveunboundedvar;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *infeasible = FALSE;

   for( c = 0; c < nconss && !(*infeasible); ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      row = NULL;

      if( consdata->nexprtrees == 0 )
      {
         assert(consdata->exprgraphnode == NULL);
         /* if we are actually linear, add the constraint as row to the LP */
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(conss[c]), SCIPconsGetName(conss[c]), consdata->lhs, consdata->rhs,
               SCIPconsIsLocal(conss[c]), FALSE , TRUE) );  /*lint !e613*/
         SCIP_CALL( SCIPaddVarsToRow(scip, row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
         SCIP_CALL( SCIPreleaseRow (scip, &row) );
         continue;
      }

      /* setup reference points for each exprtree */
      SCIP_CALL( SCIPallocBufferArray(scip, &x, consdata->nexprtrees) );
      haveunboundedvar = FALSE;
      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         nvars = SCIPexprtreeGetNVars(consdata->exprtrees[j]);

         SCIP_CALL( SCIPallocBufferArray(scip, &x[j], nvars) );  /*lint !e866*/
         for( i = 0; i < nvars; ++i )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[j])[i];
            assert(var != NULL);
            /* use midpoint as reference value, if both bounds are finite
             * otherwise use 0.0, projected on bounds
             */
            if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
            {
               if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
               {
                  x[j][i] = 0.0;
                  haveunboundedvar = TRUE;
               }
               else
                  x[j][i] = MIN(0.0, SCIPvarGetUbGlobal(var));  /*lint !e666*/
            }
            else
            {
               if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
                  x[j][i] = MAX(0.0, SCIPvarGetLbGlobal(var));  /*lint !e666*/
               else
               {
                  x[j][i] = (SCIPvarGetLbGlobal(var) + SCIPvarGetUbGlobal(var)) / 2.0;
                  /* shift refpoint into [-INITLPMAXVARVAL, INITLPMAXVARVAL], if bounds allow */
                  if( x[j][i] < -INITLPMAXVARVAL && SCIPvarGetUbGlobal(var) >= -INITLPMAXVARVAL )
                     x[j][i] = -INITLPMAXVARVAL;
                  else if( x[j][i] > INITLPMAXVARVAL && SCIPvarGetLbGlobal(var) <= INITLPMAXVARVAL )
                     x[j][i] =  INITLPMAXVARVAL;
               }
            }
         }
      }

      /* for inequalities that are convex or that have bounded variables, try to generate a cut */
      if( !SCIPisInfinity(scip,  consdata->rhs) && ((consdata->curvature & SCIP_EXPRCURV_CONVEX)  || !haveunboundedvar) )
      {
         SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], x, NULL, TRUE, SCIP_SIDETYPE_RIGHT, &row,
               -SCIPinfinity(scip), conshdlrdata->cutmaxrange, FALSE, FALSE) );  /*lint !e613*/

         if( row != NULL )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE /* forcecut */, infeasible) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }
      }

      if( !(*infeasible) && !SCIPisInfinity(scip, -consdata->lhs) &&
         ((consdata->curvature & SCIP_EXPRCURV_CONCAVE) || !haveunboundedvar) )
      {
         SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], x, NULL, TRUE, SCIP_SIDETYPE_LEFT, &row,
            -SCIPinfinity(scip), conshdlrdata->cutmaxrange, FALSE, FALSE) );  /*lint !e613*/

         if( row != NULL )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE /* forcecut */, infeasible) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }
      }

      /* @todo could add more linearizations for convex or multivariate concave inequ. */

      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         SCIPfreeBufferArray(scip, &x[j]);
      }
      SCIPfreeBufferArray(scip, &x);
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          solviolbounds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, NULL, &solviolbounds, &maxviolcon) );

   /* it can happen here that the solution violates some bound - we then just don't separate, see also discussion in issue #627 */
   if( solviolbounds )
      return SCIP_OKAY;

   /* nothing violated -> nothing to separate */
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   /* at root, check if we want to solve the NLP relaxation and use its solutions as reference point
    * if there is something convex, then linearizing in the solution of the NLP relaxation can be very useful
    */
   if( SCIPgetDepth(scip) == 0 && !conshdlrdata->sepanlp &&
      (SCIPgetNContVars(scip) >= conshdlrdata->sepanlpmincont * SCIPgetNVars(scip) || (SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY && conshdlrdata->sepanlpmincont <= 1.0)) &&
      SCIPisNLPConstructed(scip) && SCIPgetNNlpis(scip) > 0 )
   {
      SCIP_CONSDATA* consdata;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Bool solvednlp;   /* whether we invoked an NLP solve here */
      int c;

      solstat = SCIPgetNLPSolstat(scip);
      solvednlp = FALSE;
      if( solstat == SCIP_NLPSOLSTAT_UNKNOWN )
      {
         /* NLP is not solved yet, so we might want to do this
          * but first check whether there is a violated constraint side which corresponds to a convex function
          * @todo put this check into initsol and update via consenable/consdisable
          */
         for( c = 0; c < nconss; ++c )
         {
            assert(conss[c] != NULL);  /*lint !e613*/

            consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
            assert(consdata != NULL);

            /* skip feasible constraints */
            if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
               continue;

            /* make sure curvature has been checked */
            SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );  /*lint !e613*/

            if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && (consdata->curvature & SCIP_EXPRCURV_CONVEX )) ||
               ( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && (consdata->curvature & SCIP_EXPRCURV_CONCAVE)) )
               break;
         }

         if( c < nconss )
         {
            /* try to solve NLP and update solstat */

            /* ensure linear conss are in NLP */
            if( conshdlrdata->subnlpheur != NULL )
            {
               SCIP_CALL( SCIPaddLinearConsToNlpHeurSubNlp(scip, conshdlrdata->subnlpheur, TRUE, TRUE) );
            }

            /* set LP solution as starting values, if available */
            if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );
            }

            /* SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, 1) ); */
            SCIP_CALL( SCIPsolveNLP(scip) );

            solstat = SCIPgetNLPSolstat(scip);
            SCIPdebugMsg(scip, "solved NLP relax, solution status: %d\n", solstat);

            solvednlp = TRUE;
         }
      }

      conshdlrdata->sepanlp = TRUE;

      if( solstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE )
      {
         SCIPdebugMsg(scip, "NLP relaxation is globally infeasible, thus can cutoff node\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         /* if we have feasible NLP solution, generate linearization cuts there */
         SCIP_Bool lpsolseparated;
         SCIP_SOL* nlpsol;

         SCIP_CALL( SCIPcreateNLPSol(scip, &nlpsol, NULL) );
         assert(nlpsol != NULL);

         /* if we solved the NLP and solution is integral, then pass it to trysol heuristic */
         if( solvednlp && conshdlrdata->trysolheur != NULL )
         {
            int nfracvars;

            nfracvars = 0;
            if( SCIPgetNBinVars(scip) > 0 || SCIPgetNIntVars(scip) > 0 )
            {
               SCIP_CALL( SCIPgetNLPFracVars(scip, NULL, NULL, NULL, &nfracvars, NULL) );
            }

            if( nfracvars == 0 )
            {
               SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, nlpsol) );
            }
         }

         SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, nlpsol, &lpsolseparated, SCIPgetSepaMinEfficacy(scip)) );

         SCIP_CALL( SCIPfreeSol(scip, &nlpsol) );

         /* if a cut that separated the LP solution was added, then return, otherwise continue with usual separation in LP solution */
         if( lpsolseparated )
         {
            SCIPdebugMsg(scip, "linearization cuts separate LP solution\n");

            *result = SCIP_SEPARATED;

            return SCIP_OKAY;
         }
      }
   }
   /* if we do not want to try solving the NLP, or have no NLP, or have no NLP solver, or solving the NLP failed,
    * or separating with NLP solution as reference point failed, then try (again) with LP solution as reference point
    */

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, TRUE, SCIPgetSepaMinEfficacy(scip), FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolNonlinear)
{
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          solviolbounds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(sol != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, sol, &solviolbounds, &maxviolcon) );

   /* odd, if this happens for non-LP solutions, but luckily we can just give up here */
   if( solviolbounds )
      return SCIP_OKAY;

   /* nothing violated -> nothing to separate */
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   /* computeViolations already evaluated all constraints, so can pass newsol = FALSE here
    * in contrast to Sepalp, a sol != NULL is not projected onto the box in computeViolation
    */
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, FALSE, SCIPgetSepaMinEfficacy(scip), FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpNonlinear)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, solinfeasible, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxNonlinear)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, solinfeasible, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsNonlinear)
{
   SCIP_CONS*         maxviolcons;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        propresult;
   SCIP_VAR*          var;
   int                dummy;
   int                nnotify;
   int                c;
   int                i;
   int                j;
   SCIP_Bool          solviolbounds;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, NULL, &solviolbounds, &maxviolcons) );

   /* we enforce a pseudo-solution, which should be within (read: at) bounds by definition */
   assert(!solviolbounds);

   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   SCIPdebugMsg(scip, "enfops with max violation in cons <%s>\n", SCIPconsGetName(maxviolcons));

   /* we propagate constraints only if they are active and enforcing by branching only does not seem much effective */
   assert(SCIPconsIsActive(maxviolcons));

   /* run domain propagation */
   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, TRUE, &propresult, &dummy, &dummy) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* We are not feasible and we cannot prove that the whole node is infeasible -> collect all variables in violated
    * constraints for branching. */
   nnotify = 0;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         var = consdata->linvars[i];
         if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++nnotify;
         }
      }

      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]); ++i )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[j])[i];
            if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            {
               SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
               ++nnotify;
            }
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
}  /*lint !e715*/


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   int                c;
   SCIP_Bool          maypropfeasible; /* whether we may be able to propose a feasible solution */
   SCIP_Bool          solviolbounds;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   /* during presolve, we do not have exprtrees in the constraints, but we can get values from the expression graph, if we have evaluated it */
   if( SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) <= SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_Real* varvals;

      assert(conshdlrdata->exprgraph != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &varvals, SCIPexprgraphGetNVars(conshdlrdata->exprgraph)) );
      SCIP_CALL( SCIPgetSolVals(scip, sol, SCIPexprgraphGetNVars(conshdlrdata->exprgraph), (SCIP_VAR**)SCIPexprgraphGetVars(conshdlrdata->exprgraph), varvals) );

      SCIP_CALL( SCIPexprgraphEval(conshdlrdata->exprgraph, varvals) );

      SCIPfreeBufferArray(scip, &varvals);
   }

   /* @todo adapt proposeFeasibleSolution to function also during presolving */
   maxviol = 0.0;
   maypropfeasible = conshdlrdata->linfeasshift && (conshdlrdata->trysolheur != NULL) &&
      SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED &&
      (SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) > SCIP_STAGE_EXITPRESOLVE) &&
      SCIPgetStage(scip) <= SCIP_STAGE_SOLVING;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( computeViolation(scip, conshdlr, conss[c], sol, &solviolbounds) );
      assert(!solviolbounds);  /* see also issue #627 */

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         *result = SCIP_INFEASIBLE;
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");
            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
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
            maxviol = MAX(consdata->lhsviol, consdata->rhsviol);

         /* do not try to shift linear variables if activity is at infinity (leads to setting variable to infinity in solution, which is not allowed) */
         if( maypropfeasible && SCIPisInfinity(scip, REALABS(consdata->activity)) )
            maypropfeasible = FALSE;

         if( maypropfeasible )
         {
            /* update information on linear variables that may be in- or decreased */
            if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
               consdataFindUnlockedLinearVar(scip, consdata);

            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               /* check if there is a variable which may help to get the left hand side satisfied
                * if there is no such var, then we cannot get feasible */
               if( !(consdata->linvar_mayincrease >= 0 && consdata->lincoefs[consdata->linvar_mayincrease] > 0.0) &&
                  ! (consdata->linvar_maydecrease >= 0 && consdata->lincoefs[consdata->linvar_maydecrease] < 0.0) )
                  maypropfeasible = FALSE;
            }
            else
            {
               assert(SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)));
               /* check if there is a variable which may help to get the right hand side satisfied
                * if there is no such var, then we cannot get feasible */
               if( !(consdata->linvar_mayincrease >= 0 && consdata->lincoefs[consdata->linvar_mayincrease] < 0.0) &&
                  ! (consdata->linvar_maydecrease >= 0 && consdata->lincoefs[consdata->linvar_maydecrease] > 0.0) )
                  maypropfeasible = FALSE;
            }
         }
      }
      else
      {
         /* SCIPdebugMsg(scip, "constraint <%s> is feasible (%g, %g) in check, activity = %g, sides = [%g, %g]\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->activity, consdata->lhs, consdata->rhs); */
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
}  /*lint !e715*/


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropNonlinear)
{
   int dummy;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nmarkedconss, TRUE, result, &dummy, &dummy) );

   return SCIP_OKAY;
}  /*lint !e715*/

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        propresult;
   SCIP_Bool          havechange;
   SCIP_Bool          domainerror;
   SCIP_Bool          havegraphchange;
   SCIP_Bool          tryupgrades;
   int                c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   havegraphchange = FALSE;

   if( !conshdlrdata->isremovedfixings )
   {
      SCIP_CALL( removeFixedNonlinearVariables(scip, conshdlr) );
      assert(conshdlrdata->isremovedfixings);

      havegraphchange = TRUE;
   }

   SCIP_CALL( SCIPexprgraphSimplify(conshdlrdata->exprgraph, SCIPgetMessagehdlr(scip), SCIPepsilon(scip), conshdlrdata->maxexpansionexponent, &havechange, &domainerror) );
   SCIPdebugMsg(scip, "expression graph simplifier found %schange, domain error = %u\n", havechange ? "" : "no ", domainerror);

   /* if simplifier found some undefined expression, then declare problem as infeasible
    * usually, this should be discovered during domain propagation already, but since that is using interval arithmetics,
    *   it may overestimate in a way that actually undefined expressions still get a value assigned (e.g., 0^(-1) = [-inf,inf])
    */
   if( domainerror )
      *result = SCIP_CUTOFF;

   havegraphchange |= havechange;

   /* if graph has changed, then we will try upgrades, otherwise we only do for changing or not-yet-presolved constraints */
   tryupgrades = havegraphchange;

   /* remove fix vars, do some algebraic manipulation, etc; this loop need to finish, even if a cutoff is found because data
    * might be unconsistent otherwise (i.e. some asserts might pop later, e.g. exitpresol, etc)
    */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      SCIPdebugMsg(scip, "process constraint <%s>\n", SCIPconsGetName(conss[c]));
      SCIPdebugPrintCons(scip, conss[c], NULL);

      havechange = FALSE;

      if( !consdata->isremovedfixingslin )
      {
         SCIP_CALL( removeFixedLinearVariables(scip, conss[c]) );
         assert(consdata->isremovedfixingslin);
         havechange = TRUE;
      }

      /* the reductions below require the constraint nonlinear function to be in the expression graph, which is only the
       * case for active constraints
       */
      if( !SCIPconsIsActive(conss[c]) )
         continue;

      if( !consdata->ispresolved || havegraphchange )
      {
         SCIP_Bool infeasible;

         SCIP_CALL( splitOffLinearPart(scip, conshdlr, conss[c], &infeasible) );

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            continue;
         }
      }

      if( consdata->nlinvars == 0 && consdata->exprgraphnode == NULL )
      {
         /* all variables fixed or removed, constraint function is 0.0 now */
         if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasPositive(scip, consdata->lhs)) ||
            ( !SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasNegative(scip, consdata->rhs)) )
         {
            /* left hand side positive or right hand side negative */
            SCIPdebugMsg(scip, "constraint <%s> is constant and infeasible\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            *result = SCIP_CUTOFF;
         }
         else
         {
            /* left and right hand side are consistent */
            SCIPdebugMsg(scip, "constraint <%s> is constant and feasible, deleting\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            ++*ndelconss;

            if( *result != SCIP_CUTOFF )
               *result = SCIP_SUCCESS;
            continue;
         }
      }

      /* remember that we want to call upgrade methods for the current constraint */
      if( havechange )
         consdata->ispresolved = FALSE;

      /* if a constraint is not finished presolving yet, then we will try upgrade methods */
      if( !consdata->ispresolved )
         tryupgrades = TRUE;
   }

   /* if a cutoff was found, return; data is consistent at this point */
   if( *result == SCIP_CUTOFF )
      return SCIP_OKAY;

   if( tryupgrades )
   {
      /* upgrade methods may look at expression graph bounds, which are not present in the first presolving round yet and may be invalid in later rounds (e.g., due to probing) */
      SCIP_CALL( SCIPexprgraphPropagateVarBounds(conshdlrdata->exprgraph, INTERVALINFTY, TRUE, &domainerror) );

      if( domainerror )
      {
         SCIPdebugMsg(scip, "propagating variable bounds through expression graph found that some expressions cannot be evaluated w.r.t. current bounds, thus cutoff\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      for( c = 0; c < nconss; ++c )
      {
         consdata = SCIPconsGetData(conss[c]);  /*lint !e794*/
         assert(consdata != NULL);

         /* call upgrade methods if constraint was not presolved, has been changed, or the expression graph has changed */
         if( !consdata->ispresolved || havegraphchange )
         {
            SCIP_Bool upgraded;

            SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss) );  /*lint !e794*/
            if( upgraded )
            {
               *result = SCIP_SUCCESS;
               continue;
            }
         }

         consdata->ispresolved = TRUE;
      }
   }

   /* run domain propagation (if updated bounds in graph above, then can skip cleanup) */
   if( (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 )
   {
      SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, !tryupgrades, &propresult, nchgbds, ndelconss) );
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
   }

   if( conshdlrdata->reformulate && !conshdlrdata->assumeconvex )
   {
      /* if other presolvers did not find enough changes for another presolving round,
       * then try the reformulations (replacing products with binaries, disaggregation, setting default variable bounds)
       * otherwise, we wait with these
       */
      if( SCIPisPresolveFinished(scip) || (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
      {
         int naddconssbefore;

         SCIPdebugMsg(scip, "reformulating expression graph\n");

         naddconssbefore = conshdlrdata->naddedreformconss;
         SCIP_CALL( reformulate(scip, conshdlr, conss, nconss, &conshdlrdata->naddedreformconss) );

         if( conshdlrdata->naddedreformconss > naddconssbefore )
         {
            *result = SCIP_SUCCESS;
            *naddconss += conshdlrdata->naddedreformconss - naddconssbefore;

            /* if expression graph changed, ensure that we apply all presolving techniques (esp. upgrades) in next round again */
            for( c = 0; c < nconss; ++c )
            {
               assert(conss[c] != NULL);  /*lint !e794*/

               consdata = SCIPconsGetData(conss[c]);  /*lint !e794*/
               assert(consdata != NULL);

               consdata->ispresolved = FALSE;
            }
         }
      }
   }

   return SCIP_OKAY;
}  /*lint !e715*/


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockNonlinear)
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool      havelhs;
   SCIP_Bool      haverhs;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);

   /* variable locking for nonlinear part is done w.r.t. variables in the expression graph
    * since only active constraints have their nonlinear part in the expression graph, we can lock only active constraints
    */
   assert(SCIPconsIsActive(cons) || SCIPconsIsDeleted(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   havelhs = !SCIPisInfinity(scip, -consdata->lhs);
   haverhs = !SCIPisInfinity(scip,  consdata->rhs);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( consdata->lincoefs[i] > 0 )
      {
         if( havelhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
         if( haverhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( havelhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
         if( haverhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveNonlinear)
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

   SCIPdebugMsg(scip, "activate cons <%s>\n", SCIPconsGetName(cons));

   if( consdata->nexprtrees > 0 )
   {
      SCIP_Bool exprtreeisnew;

      assert(consdata->exprgraphnode == NULL);

      /* add exprtrees to expression graph */
      SCIP_CALL( SCIPexprgraphAddExprtreeSum(conshdlrdata->exprgraph, consdata->nexprtrees, consdata->exprtrees, consdata->nonlincoefs, &consdata->exprgraphnode, &exprtreeisnew) );
      assert(consdata->exprgraphnode != NULL);
      /* @todo do something with exprtreeisnew? */

      /* if during presolving, then forget expression trees */
      if( SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) < SCIP_STAGE_EXITPRESOLVE )
      {
         SCIP_CALL( consdataSetExprtrees(scip, consdata, 0, NULL, NULL, FALSE) );
      }

      /* remember that we should run reformulation again */
      conshdlrdata->isreformulated = FALSE;

      /* remember that we should force backward propagation on our subgraph propagating the next time,
       * so possible domain restrictions are propagated into variable bounds
       */
      consdata->forcebackprop = TRUE;
   }
   else if( consdata->exprgraphnode != NULL )
   {
      /* if constraint already comes with node in expression graph, then also remember that we should run reformulation again */
      conshdlrdata->isreformulated = FALSE;

      /* remember that we should force backward propagation on our subgraph propagating the next time,
       * so possible domain restrictions are propagated into variable bounds
       */
      consdata->forcebackprop = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveNonlinear)
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
   assert(consdata->exprgraphnode != NULL || consdata->nexprtrees == 0);

   SCIPdebugMsg(scip, "deactivate cons <%s>\n", SCIPconsGetName(cons));

   if( consdata->exprgraphnode != NULL )
   {
      if( consdata->nexprtrees == 0 )
      {
         /* during presolving, the exprtrees in the constraint are removed, so put them back before releasing the exprgraphnode */
         SCIP_EXPRTREE* exprtree;

         /* if only presolve is run and problem is found infeasible there, then constraints may not be deactivated there, but in a later call to freeTransform */
         /* @todo if infeasible in presolve, will constraints be deactivated still in presolving stage, or in exitpre? */
         assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) <= SCIP_STAGE_EXITPRESOLVE || SCIPgetStage(scip) == SCIP_STAGE_FREETRANS);

         SCIP_CALL( SCIPexprgraphGetTree(conshdlrdata->exprgraph, consdata->exprgraphnode, &exprtree) );
         SCIP_CALL( consdataSetExprtrees(scip, consdata, 1, &exprtree, NULL, FALSE) );
      }

      SCIP_CALL( SCIPexprgraphReleaseNode(conshdlrdata->exprgraph, &consdata->exprgraphnode) );
   }

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int i;

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

   SCIPdebugMsg(scip, "enable cons <%s>\n", SCIPconsGetName(cons));

   if( consdata->exprgraphnode != NULL )
   {
      /* enable node of expression in expression graph */
      SCIPexprgraphEnableNode(conshdlrdata->exprgraph, consdata->exprgraphnode);
   }

   /* enable event catching for linear variables */
   consdata->isremovedfixingslin = TRUE;
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( catchLinearVarEvents(scip, cons, i) );

      consdata->isremovedfixingslin = consdata->isremovedfixingslin && SCIPvarIsActive(consdata->linvars[i]);
   }

   return SCIP_OKAY;
}

/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lineventdata != NULL || consdata->nlinvars == 0);

   SCIPdebugMsg(scip, "disable cons <%s>\n", SCIPconsGetName(cons));

   /* disable node of expression in expression graph */
   if( consdata->exprgraphnode != NULL )
   {
      SCIPexprgraphDisableNode(conshdlrdata->exprgraph, consdata->exprgraphnode);
   }

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( dropLinearVarEvents(scip, cons, i) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintNonlinear)
{
   SCIP_CONSDATA* consdata;
   int            j;

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
   if( consdata->nlinvars == 0 && consdata->nexprtrees == 0 && consdata->exprgraphnode == 0 )
   {
      SCIPinfoMessage(scip, file, "0 ");
   }
   else
   {
      if( consdata->nexprtrees > 0 )
      {
         for( j = 0; j < consdata->nexprtrees; ++j )
         {
            if( j > 0 || consdata->nonlincoefs[j] != 1.0 )
               SCIPinfoMessage(scip, file, " %+.20g ", consdata->nonlincoefs[j]);
            SCIP_CALL( SCIPexprtreePrintWithNames(consdata->exprtrees[j], SCIPgetMessagehdlr(scip), file) );
         }
      }
      else if( consdata->exprgraphnode != NULL )
      {
         SCIP_CONSHDLRDATA* conshdlrdata;
         SCIP_EXPRTREE* tree;

         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);
         SCIP_CALL( SCIPexprgraphGetTree(conshdlrdata->exprgraph, consdata->exprgraphnode, &tree) );

         SCIP_CALL( SCIPexprtreePrintWithNames(tree, SCIPgetMessagehdlr(scip), file) );

         SCIP_CALL( SCIPexprtreeFree(&tree) );
      }

      for( j = 0; j < consdata->nlinvars; ++j )
      {
         SCIPinfoMessage(scip, file, "%+.15g<%s>[%c] ", consdata->lincoefs[j], SCIPvarGetName(consdata->linvars[j]),
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_BINARY ? 'B' :
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
      }
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

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyNonlinear)
{
   SCIP_CONSDATA*    consdata;
   SCIP_CONSDATA*    targetconsdata;
   SCIP_VAR**        linvars;
   SCIP_Real*        nonlincoefs;
   SCIP_EXPRTREE**   exprtrees;
   int               nexprtrees;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourceconshdlr != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);

   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);

   linvars = NULL;
   exprtrees = NULL;

   *valid = TRUE;

   if( consdata->nlinvars != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &linvars, consdata->nlinvars) );
      for( i = 0; i < consdata->nlinvars && *valid; ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->linvars[i], &linvars[i], varmap, consmap, global, valid) );
         assert(!*valid || linvars[i] != NULL);
      }
   }

   nexprtrees = 0;
   nonlincoefs = NULL;

   if( *valid && consdata->nexprtrees > 0 )
   {
      SCIP_VAR** nonlinvars;

      nonlincoefs = consdata->nonlincoefs;
      nexprtrees = consdata->nexprtrees;

      SCIP_CALL( SCIPallocBufferArray(sourcescip, &exprtrees, nexprtrees) );
      BMSclearMemoryArray(exprtrees, nexprtrees);
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &nonlinvars, SCIPexprtreeGetNVars(consdata->exprtrees[0])) );

      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         SCIP_CALL( SCIPreallocBufferArray(sourcescip, &nonlinvars, SCIPexprtreeGetNVars(consdata->exprtrees[j])) );
         for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]) && *valid; ++i )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, SCIPexprtreeGetVars(consdata->exprtrees[j])[i], &nonlinvars[i], varmap, consmap, global, valid) );
            assert(!*valid || nonlinvars[i] != NULL);
         }

         if( *valid )
         {
            SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &exprtrees[j], consdata->exprtrees[j]) );
            SCIP_CALL( SCIPexprtreeSetVars(exprtrees[j], SCIPexprtreeGetNVars(consdata->exprtrees[j]), nonlinvars) );
         }
         else
            break;
      }

      SCIPfreeBufferArray(sourcescip, &nonlinvars);
   }

   if( *valid && consdata->nexprtrees == 0 && consdata->exprgraphnode != NULL )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_VAR** nonlinvars;

      conshdlrdata = SCIPconshdlrGetData(sourceconshdlr);

      nexprtrees = 1;
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &exprtrees, 1) );

      SCIP_CALL( SCIPexprgraphGetTree(conshdlrdata->exprgraph, consdata->exprgraphnode, &exprtrees[0]) );

      nonlinvars = SCIPexprtreeGetVars(exprtrees[0]);
      for( i = 0; i < SCIPexprtreeGetNVars(exprtrees[0]); ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, nonlinvars[i], &nonlinvars[i], varmap, consmap, global, valid) );
         assert(!*valid || nonlinvars[i] != NULL);
      }
   }

   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, name ? name : SCIPconsGetName(sourcecons),
            consdata->nlinvars, linvars, consdata->lincoefs,
            nexprtrees, exprtrees, nonlincoefs,
            consdata->lhs, consdata->rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

      /* copy information on curvature */
      targetconsdata = SCIPconsGetData(*cons);
      targetconsdata->curvature     = consdata->curvature;
      targetconsdata->iscurvchecked = consdata->iscurvchecked && global; /* if the copy is local, then curvature may change (get stronger) */
   }

   SCIPfreeBufferArrayNull(sourcescip, &linvars);
   if( exprtrees != NULL )
   {
      for( j = 0; j < nexprtrees; ++j )
      {
         if( exprtrees[j] != NULL )
         {
            SCIP_CALL( SCIPexprtreeFree(&exprtrees[j]) );
         }
      }
      SCIPfreeBufferArray(sourcescip, &exprtrees);
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int cnt;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *success = TRUE;

   if( varssize < consdata->nlinvars )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   BMScopyMemoryArray(vars, consdata->linvars, consdata->nlinvars);
   cnt = consdata->nlinvars;

   if( consdata->exprgraphnode != NULL )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int* varsusage;
      int i;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &varsusage, SCIPexprgraphGetNVars(conshdlrdata->exprgraph)) );

      SCIPexprgraphGetSubtreeVarsUsage(conshdlrdata->exprgraph, consdata->exprgraphnode, varsusage);

      for( i = 0; i < SCIPexprgraphGetNVars(conshdlrdata->exprgraph); ++i )
      {
         if( varsusage[i] == 0 )
            continue;

         if( cnt >= varssize )
         {
            *success = FALSE;
            break;
         }

         vars[cnt] = (SCIP_VAR*)(SCIPexprgraphGetVars(conshdlrdata->exprgraph)[i]);
         ++cnt;
      }

      SCIPfreeBufferArray(scip, &varsusage);
   }
   else
   {
      SCIP_VAR** exprvars;
      int nexprvars;
      int e;

      for( e = 0; e < consdata->nexprtrees; ++e )
      {
         exprvars  = SCIPexprtreeGetVars(consdata->exprtrees[e]);
         nexprvars = SCIPexprtreeGetNVars(consdata->exprtrees[e]);
         assert(exprvars != NULL || nexprvars == 0);

         if( cnt + nexprvars > varssize )
         {
            *success = FALSE;
            break;
         }

         BMScopyMemoryArray(&vars[cnt], exprvars, nexprvars);  /*lint !e866*/
         cnt += nexprvars;
      }
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *nvars = consdata->nlinvars;

   if( consdata->exprgraphnode != NULL )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      int* varsusage;
      int i;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &varsusage, SCIPexprgraphGetNVars(conshdlrdata->exprgraph)) );

      SCIPexprgraphGetSubtreeVarsUsage(conshdlrdata->exprgraph, consdata->exprgraphnode, varsusage);

      for( i = 0; i < SCIPexprgraphGetNVars(conshdlrdata->exprgraph); ++i )
         if( varsusage[i] > 0 )
            ++*nvars;

      SCIPfreeBufferArray(scip, &varsusage);
   }
   else
   {
      int e;

      for( e = 0; e < consdata->nexprtrees; ++e )
         *nvars += SCIPexprtreeGetNVars(consdata->exprtrees[e]);
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseNonlinear)
{  /*lint --e{715}*/
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* expr;
   SCIP_VAR** exprvars;
   SCIP_RETCODE retcode;
   int        nvars;
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   const char* endptr;
   char*       nonconstendptr;
   const char* exprstart;
   const char* exprlastchar;
   int* varnames;
   int* curvarname;
   int varnameslength;
   int i;

   SCIPdebugMsg(scip, "cons_nonlinear::consparse parsing %s\n",str);

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* return if string empty */
   if( !*str )
      return SCIP_OKAY;

   endptr = str;

   expr = NULL;
   nvars = 0;

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   /* parse constraint to get lhs, rhs, and expression in between (from cons_linear.c::consparse, but parsing whole string first, then getting expression) */

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, &nonconstendptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_READERROR;
      }
      endptr = nonconstendptr;

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the first coefficient, but not a left-hand-side */
         lhs = -SCIPinfinity(scip);
      }
      else
      {
         /* it was indeed a left-hand-side, so continue parsing after it */
         str = endptr + 2;

         /* ignore whitespace */
         while( isspace((unsigned char)*str) )
            ++str;
      }
   }

   /* Move endptr forward until we find end of expression */
   while( !(strncmp(endptr, "[free]", 6) == 0)    &&
          !(endptr[0] == '<' && endptr[1] == '=') &&
          !(endptr[0] == '=' && endptr[1] == '=') &&
          !(endptr[0] == '>' && endptr[1] == '=') &&
          !(endptr[0] == '\0') )
      ++endptr;

   exprstart = str;
   exprlastchar = endptr - 1;

   *success = FALSE;
   str = endptr;

   /* check for left or right hand side */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for free constraint */
   if( strncmp(str, "[free]", 6) == 0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIPerrorMessage("cannot have left hand side and [free] status \n");
         return SCIP_OKAY;
      }
      (*success) = TRUE;
   }
   else
   {
      switch( *str )
      {
         case '<':
            *success = SCIPstrToRealValue(str+2, &rhs, &nonconstendptr);
            break;
         case '=':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have == on rhs if there was a <= on lhs\n");
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &rhs, &nonconstendptr);
               lhs = rhs;
            }
            break;
         case '>':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have => on rhs if there was a <= on lhs\n");
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &lhs, &nonconstendptr);
               break;
            }
         case '\0':
            *success = TRUE;
            break;
         default:
            SCIPerrorMessage("unexpected character %c\n", *str);
            return SCIP_OKAY;
      }
   }

   /* alloc some space for variable names incl. indices; shouldn't be longer than expression string, and we even give it sizeof(int) times this length (plus 5) */
   varnameslength = (int) (exprlastchar - exprstart) + 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &varnames, varnameslength) );

   /* parse expression */
   retcode = SCIPexprParse(SCIPblkmem(scip), SCIPgetMessagehdlr(scip), &expr, exprstart, exprlastchar, &nvars, varnames, varnameslength);

   if( retcode != SCIP_OKAY )
   {
      SCIPfreeBufferArray(scip, &varnames);
      return retcode;
   }

   /* get SCIP variables corresponding to variable names stored in varnames buffer */
   SCIP_CALL( SCIPallocBufferArray(scip, &exprvars, nvars) );

   assert( retcode == SCIP_OKAY );
   curvarname = varnames;
   for( i = 0; i < nvars; ++i )
   {
      assert(*curvarname == i);
      ++curvarname;

      exprvars[i] = SCIPfindVar(scip, (char*)curvarname);
      if( exprvars[i] == NULL )
      {
         SCIPerrorMessage("Unknown SCIP variable <%s> encountered in expression.\n", (char*)curvarname);
         retcode = SCIP_READERROR;
         goto TERMINATE;
      }

      curvarname += (strlen((char*)curvarname) + 1)/sizeof(int) + 1;
   }

   /* create expression tree */
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, nvars, 0, NULL) );
   SCIP_CALL( SCIPexprtreeSetVars(exprtree, nvars, exprvars) );

   /* create constraint */
   SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, name,
      0, NULL, NULL,
      1, &exprtree, NULL,
      lhs, rhs,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   SCIPdebugMsg(scip, "created nonlinear constraint:\n");
   SCIPdebugPrintCons(scip, *cons, NULL);

   SCIP_CALL( SCIPexprtreeFree(&exprtree) );

 TERMINATE:
   SCIPfreeBufferArray(scip, &exprvars);
   SCIPfreeBufferArray(scip, &varnames);

   return retcode;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for nonlinear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create nonlinear constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpNonlinear, consEnfopsNonlinear, consCheckNonlinear, consLockNonlinear,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveNonlinear) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyNonlinear, consCopyNonlinear) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveNonlinear) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteNonlinear) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableNonlinear) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableNonlinear) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitNonlinear) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreNonlinear) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolNonlinear) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeNonlinear) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsNonlinear) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsNonlinear) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitNonlinear) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreNonlinear) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolNonlinear) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpNonlinear) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolNonlinear, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintNonlinear) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropNonlinear, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpNonlinear, consSepasolNonlinear, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransNonlinear) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseNonlinear) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxNonlinear) );

   /* add nonlinear constraint handler parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/cutmaxrange",
         "maximal coef range of a cut (maximal coefficient divided by minimal coefficient) in order to be added to LP relaxation",
         &conshdlrdata->cutmaxrange, FALSE, 1e+7, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/linfeasshift",
         "whether to try to make solutions in check function feasible by shifting a linear variable (esp. useful if constraint was actually objective function)",
         &conshdlrdata->linfeasshift, FALSE, TRUE, NULL, NULL) );

#if 0 /* don't have any expensive checks yet, so we disable this parameter for now */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkconvexexpensive",
         "whether to apply expensive curvature checking methods",
         &conshdlrdata->checkconvexexpensive, FALSE, TRUE, NULL, NULL) );
#endif

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/assumeconvex",
         "whether to assume that nonlinear functions in inequalities (<=) are convex (disables reformulation)",
         &conshdlrdata->assumeconvex, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a single constraint within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 1, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformulate",
         "whether to reformulate expression graph",
         &conshdlrdata->reformulate, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxexpansionexponent",
         "maximal exponent where still expanding non-monomial polynomials in expression simplification",
         &conshdlrdata->maxexpansionexponent, TRUE, 2, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/sepanlpmincont",
         "minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation",
         &conshdlrdata->sepanlpmincont, FALSE, 1.0, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/enfocutsremovable",
         "are cuts added during enforcement removable from the LP in the same node?",
         &conshdlrdata->enfocutsremovable, TRUE, FALSE, NULL, NULL) );

   conshdlrdata->linvareventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(conshdlrdata->linvareventhdlr), CONSHDLR_NAME"_boundchange", "signals a bound change to a nonlinear constraint",
         processLinearVarEvent, NULL) );
   assert(conshdlrdata->linvareventhdlr != NULL);

   conshdlrdata->nonlinvareventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(conshdlrdata->nonlinvareventhdlr), CONSHDLR_NAME"_boundchange2", "signals a bound change to a nonlinear constraint handler",
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

/** includes a nonlinear constraint upgrade method into the nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNonlinconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_NONLINCONSUPGD((*nonlinconsupgd)),/**< method to call for upgrading nonlinear constraint, or NULL */
   SCIP_DECL_EXPRGRAPHNODEREFORM((*nodereform)),/**< method to call for reformulating expression graph node, or NULL */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*        conshdlr;
   SCIP_CONSHDLRDATA*    conshdlrdata;
   SCIP_NLCONSUPGRADE*   nlconsupgrade;
   char                  paramname[SCIP_MAXSTRLEN];
   char                  paramdesc[SCIP_MAXSTRLEN];
   int                   i;

   assert(conshdlrname != NULL );

   /* ignore empty upgrade functions */
   if( nonlinconsupgd == NULL && nodereform == NULL )
      return SCIP_OKAY;

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether upgrade method exists already */
   for( i = conshdlrdata->nnlconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->nlconsupgrades[i]->nlconsupgd == nonlinconsupgd && conshdlrdata->nlconsupgrades[i]->nodereform == nodereform)
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage(scip, "Try to add already known upgrade method pair (%p,%p) for constraint handler <%s>.\n", nonlinconsupgd, nodereform, conshdlrname); /*lint !e611*/
#endif
         return SCIP_OKAY;
      }
   }

   /* create a nonlinear constraint upgrade data object */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlconsupgrade) );
   nlconsupgrade->nlconsupgd = nonlinconsupgd;
   nlconsupgrade->nodereform = nodereform;
   nlconsupgrade->priority   = priority;
   nlconsupgrade->active     = active;

   /* insert nonlinear constraint upgrade method into constraint handler data */
   assert(conshdlrdata->nnlconsupgrades <= conshdlrdata->nlconsupgradessize);
   if( conshdlrdata->nnlconsupgrades+1 > conshdlrdata->nlconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->nnlconsupgrades+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->nlconsupgrades, conshdlrdata->nnlconsupgrades, newsize) );
      conshdlrdata->nlconsupgradessize = newsize;
   }
   assert(conshdlrdata->nnlconsupgrades+1 <= conshdlrdata->nlconsupgradessize);

   for( i = conshdlrdata->nnlconsupgrades; i > 0 && conshdlrdata->nlconsupgrades[i-1]->priority < nlconsupgrade->priority; --i )
      conshdlrdata->nlconsupgrades[i] = conshdlrdata->nlconsupgrades[i-1];
   assert(0 <= i && i <= conshdlrdata->nnlconsupgrades);
   conshdlrdata->nlconsupgrades[i] = nlconsupgrade;
   conshdlrdata->nnlconsupgrades++;

   /* adds parameter to turn on and off the upgrading step */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/" CONSHDLR_NAME "/upgrade/%s", conshdlrname);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable nonlinear upgrading for constraint handler <%s>", conshdlrname);
   SCIP_CALL( SCIPaddBoolParam(scip,
         paramname, paramdesc,
         &nlconsupgrade->active, FALSE, active, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 *  this variant takes expression trees as input
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   int                   nexprtrees,         /**< number of expression trees for nonlinear part of constraint */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees for nonlinear part of constraint */
   SCIP_Real*            nonlincoefs,        /**< coefficients for expression trees for nonlinear part, or NULL if all 1.0 */
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
   int i;

   assert(linvars  != NULL || nlinvars == 0);
   assert(lincoefs != NULL || nlinvars == 0);
   assert(exprtrees   != NULL || nexprtrees == 0);
   assert(modifiable == FALSE); /* we do not support column generation */

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreateEmpty(scip, &consdata) );

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   /* add linear variables */
   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, nlinvars) );
   for( i = 0; i < nlinvars; ++i )
   {
      if( SCIPisZero(scip, lincoefs[i]) )  /*lint !e613*/
         continue;

      SCIP_CALL( addLinearCoef(scip, *cons, linvars[i], lincoefs[i]) );  /*lint !e613*/
   }

   /* set expression trees */
   SCIP_CALL( consdataSetExprtrees(scip, consdata, nexprtrees, exprtrees, nonlincoefs, TRUE) );

   SCIPdebugMsg(scip, "created nonlinear constraint ");
   SCIPdebugPrintCons(scip, *cons, NULL);

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsNonlinear(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  this variant takes expression trees as input
 *
 *  @see SCIPcreateConsNonlinear() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   int                   nexprtrees,         /**< number of expression trees for nonlinear part of constraint */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees for nonlinear part of constraint */
   SCIP_Real*            nonlincoefs,        /**< coefficients for expression trees for nonlinear part, or NULL if all 1.0 */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, name, nlinvars, linvars, lincoefs, nexprtrees, exprtrees,
         nonlincoefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 * this variant takes a node of the expression graph as input and can only be used during presolving
 * it is assumed that the nonlinear constraint will be added to the transformed problem short after creation
 * the given exprgraphnode is captured in this method
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsNonlinear2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   SCIP_EXPRGRAPHNODE*   exprgraphnode,      /**< expression graph node associated to nonlinear expression */
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
   int i;

   assert(modifiable == FALSE); /* we do not support column generation */
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreateEmpty(scip, &consdata) );

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   /* add linear variables */
   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, nlinvars) );
   for( i = 0; i < nlinvars; ++i )
   {
      if( SCIPisZero(scip, lincoefs[i]) )
         continue;

      SCIP_CALL( addLinearCoef(scip, *cons, linvars[i], lincoefs[i]) );
   }

   /* set expression graph node */
   if( exprgraphnode != NULL )
   {
      consdata->exprgraphnode = exprgraphnode;
      consdata->curvature = SCIP_EXPRCURV_UNKNOWN;
      consdata->iscurvchecked = FALSE;
      consdata->activity = SCIP_INVALID;
      SCIPexprgraphCaptureNode(exprgraphnode);
   }

   SCIPdebugMsg(scip, "created nonlinear constraint ");
   SCIPdebugPrintCons(scip, *cons, NULL);

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsNonlinear(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  this variant takes a node of the expression graph as input and can only be used during presolving
 *  it is assumed that the nonlinear constraint will be added to the transformed problem short after creation
 *  the given exprgraphnode is captured in this method
 *
 *  @see SCIPcreateConsNonlinear() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicNonlinear2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   SCIP_EXPRGRAPHNODE*   exprgraphnode,      /**< expression graph node associated to nonlinear expression */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsNonlinear2(scip, cons, name, nlinvars, linvars, lincoefs, exprgraphnode, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** adds a linear variable with coefficient to a nonlinear constraint */
SCIP_RETCODE SCIPaddLinearVarNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient of variable */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   SCIP_CALL( addLinearCoef(scip, cons, var, coef) );

   return SCIP_OKAY;
}

/** sets the expression trees in a nonlinear constraint
 * constraint must not be active yet
 */
SCIP_RETCODE SCIPsetExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< new expression trees, or NULL if nexprtrees is 0 */
   SCIP_Real*            coefs               /**< coefficients of expression trees, or NULL if all 1.0 */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsActive(cons));
   assert(SCIPconsGetData(cons) != NULL);
   assert(exprtrees != NULL || nexprtrees == 0);

   SCIP_CALL( consdataSetExprtrees(scip, SCIPconsGetData(cons), nexprtrees, exprtrees, coefs, TRUE) );

   return SCIP_OKAY;
}

/** adds expression trees to a nonlinear constraint
 * constraint must not be active yet
 */
SCIP_RETCODE SCIPaddExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< new expression trees, or NULL if nexprtrees is 0 */
   SCIP_Real*            coefs               /**< coefficients of expression trees, or NULL if all 1.0 */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsActive(cons));
   assert(SCIPconsGetData(cons) != NULL);
   assert(exprtrees != NULL || nexprtrees == 0);

   SCIP_CALL( consdataAddExprtrees(scip, SCIPconsGetData(cons), nexprtrees, exprtrees, coefs, TRUE) );

   return SCIP_OKAY;
}

/** gets the nonlinear constraint as a nonlinear row representation */
SCIP_RETCODE SCIPgetNlRowNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(nlrow != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow == NULL )
   {
      SCIP_CALL( createNlRow(scip, cons) );
   }
   assert(consdata->nlrow != NULL);
   *nlrow = consdata->nlrow;

   return SCIP_OKAY;
}

/** gets the number of variables in the linear term of a nonlinear constraint */
int SCIPgetNLinearVarsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nlinvars;
}

/** gets the variables in the linear part of a nonlinear constraint */
SCIP_VAR** SCIPgetLinearVarsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->linvars;
}

/** gets the coefficients in the linear part of a nonlinear constraint */
SCIP_Real* SCIPgetLinearCoefsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lincoefs;
}

/** gets the number of expression trees of a nonlinear constraint */
int SCIPgetNExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) > SCIP_STAGE_EXITPRESOLVE);

   return SCIPconsGetData(cons)->nexprtrees;
}

/** gets the expression trees of a nonlinear constraint */
SCIP_EXPRTREE** SCIPgetExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) > SCIP_STAGE_EXITPRESOLVE);

   return SCIPconsGetData(cons)->exprtrees;
}

/** gets the coefficients of the expression trees of a nonlinear constraint */
SCIP_Real* SCIPgetExprtreeCoefsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) > SCIP_STAGE_EXITPRESOLVE);

   return SCIPconsGetData(cons)->nonlincoefs;
}

/** gets the expression graph node of a nonlinear constraint */
SCIP_EXPRGRAPHNODE* SCIPgetExprgraphNodeNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->exprgraphnode;
}

/** gets the left hand side of a nonlinear constraint */
SCIP_Real SCIPgetLhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lhs;
}

/** gets the right hand side of a nonlinear constraint */
SCIP_Real SCIPgetRhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhs;
}

/** check the function of a nonlinear constraint for convexity/concavity, if not done yet */
SCIP_RETCODE SCIPcheckCurvatureNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( checkCurvature(scip, cons, conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );

   return SCIP_OKAY;
}

/** gets the curvature of the nonlinear function of a nonlinear constraint */
SCIP_RETCODE SCIPgetCurvatureNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             checkcurv,          /**< whether to check constraint curvature, if not checked before */
   SCIP_EXPRCURV*        curvature           /**< pointer to store curvature of constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(curvature != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( checkcurv && !consdata->iscurvchecked )
   {
      SCIP_CALL( checkCurvature(scip, cons, conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );
   }

   *curvature = consdata->curvature;

   return SCIP_OKAY;
}

/** gets the curvature of the expression trees (multiplied by their coefficient) of a nonlinear constraint */
SCIP_RETCODE SCIPgetExprtreeCurvaturesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             checkcurv,          /**< whether to check constraint curvature, if not checked before */
   SCIP_EXPRCURV**       curvatures          /**< buffer to store curvatures of exprtrees */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(curvatures != NULL);
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) > SCIP_STAGE_EXITPRESOLVE);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   assert(SCIPconsGetData(cons) != NULL);

   if( checkcurv && !consdata->iscurvchecked )
   {
      SCIP_CALL( checkCurvature(scip, cons, conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );
   }

   *curvatures = consdata->curvatures;

   return SCIP_OKAY;
}

/** computes the violation of a nonlinear constraint by a solution */
SCIP_RETCODE SCIPgetViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution which violation to calculate, or NULL for LP solution */
   SCIP_Real*            violation           /**< pointer to store violation of constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_Bool      solviolbounds;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violation != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) <= SCIP_STAGE_EXITPRESOLVE && SCIPconsIsActive(cons) )
   {
      /* @todo make available */
      SCIPwarningMessage(scip, "SCIPgetViolationNonlinear is not available for active constraints during presolve.\n");
      *violation = SCIP_INVALID;
      return SCIP_OKAY;
   }

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   SCIP_CALL( computeViolation(scip, conshdlr, cons, sol, &solviolbounds) );

   if( solviolbounds )
   {
      SCIPerrorMessage("Solution passed to SCIPgetViolationNonlinear() does not satisfy variable bounds.\n");
      return SCIP_ERROR;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violation = MAX(consdata->lhsviol, consdata->rhsviol);

   return SCIP_OKAY;
}

/** get index of a linear variable of a nonlinear constraint that may be decreased without making any other constraint infeasible, or -1 if none */
int SCIPgetLinvarMayDecreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->linvar_mayincrease == -1 && consdata->linvar_maydecrease == -1 )
      consdataFindUnlockedLinearVar(scip, consdata);

   return consdata->linvar_maydecrease;
}

/** get index of a linear variable of a nonlinear constraint that may be increased without making any other constraint infeasible, or -1 if none */
int SCIPgetLinvarMayIncreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->linvar_mayincrease == -1 && consdata->linvar_maydecrease == -1 )
      consdataFindUnlockedLinearVar(scip, consdata);

   return consdata->linvar_mayincrease;
}

/** gets expression graph of nonlinear constraint handler */
SCIP_EXPRGRAPH* SCIPgetExprgraphNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   return conshdlrdata->exprgraph;
}

/** given three points, constructs coefficient of equation for hyperplane generated by these three points
 * Three points a, b, and c are given.
 * Computes coefficients alpha, beta, gamma, and delta, such that a, b, and c, satisfy
 * alpha * x1 + beta * x2 + gamma * x3 = delta and gamma >= 0.0.
 */
SCIP_RETCODE SCIPcomputeHyperplaneThreePoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             a1,                 /**< first coordinate of a */
   SCIP_Real             a2,                 /**< second coordinate of a */
   SCIP_Real             a3,                 /**< third coordinate of a */
   SCIP_Real             b1,                 /**< first coordinate of b */
   SCIP_Real             b2,                 /**< second coordinate of b */
   SCIP_Real             b3,                 /**< third coordinate of b */
   SCIP_Real             c1,                 /**< first coordinate of c */
   SCIP_Real             c2,                 /**< second coordinate of c */
   SCIP_Real             c3,                 /**< third coordinate of c */
   SCIP_Real*            alpha,              /**< coefficient of first coordinate */
   SCIP_Real*            beta,               /**< coefficient of second coordinate */
   SCIP_Real*            gamma_,             /**< coefficient of third coordinate */
   SCIP_Real*            delta               /**< constant right-hand side */
   )
{
   assert(scip != NULL);
   assert(alpha != NULL);
   assert(beta  != NULL);
   assert(gamma_ != NULL);
   assert(delta != NULL);

   *alpha  = -b3*c2 + a3*(-b2+c2) + a2*(b3-c3) + b2*c3;
   *beta   = -(-b3*c1 + a3*(-b1+c1) + a1*(b3-c3) + b1*c3);
   *gamma_ = -a2*b1 + a1*b2 + a2*c1 - b2*c1 - a1*c2 + b1*c2;
   *delta  = -a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3;

   /* check if hyperplane contains all three points (necessary because of numerical troubles) */
   if( !SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 + *gamma_ * a3, *delta) ||
      !SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 + *gamma_ * b3, *delta) ||
      !SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 + *gamma_ * c3, *delta) )
   {
      SCIP_Real m[9];
      SCIP_Real rhs[3];
      SCIP_Real x[3];
      SCIP_Bool success;

      /* initialize matrix column-wise */
      m[0] = a1;
      m[1] = b1;
      m[2] = c1;
      m[3] = a2;
      m[4] = b2;
      m[5] = c2;
      m[6] = a3;
      m[7] = b3;
      m[8] = c3;

      rhs[0] = 1.0;
      rhs[1] = 1.0;
      rhs[2] = 1.0;

      SCIPdebugMsg(scip, "numerical troubles - try to solve the linear system via an LU factorization\n");

      /* solve the linear problem */
      SCIP_CALL( SCIPsolveLinearProb(3, m, rhs, x, &success) );

      *delta  = rhs[0];
      *alpha  = x[0];
      *beta   = x[1];
      *gamma_ = x[2];

      /* set all coefficients to zero if one of the points is not contained in the hyperplane; this ensures that we do
       * not add a cut to SCIP and that all assertions are trivially fulfilled
       */
      if( !success || !SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 + *gamma_ * a3, *delta) ||
         !SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 + *gamma_ * b3, *delta) ||
         !SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 + *gamma_ * c3, *delta) )
      {
         SCIPdebugMsg(scip, "could not resolve numerical difficulties\n");
         *delta  = 0.0;
         *alpha  = 0.0;
         *beta   = 0.0;
         *gamma_ = 0.0;
      }
   }

   if( *gamma_ < 0.0 )
   {
      *alpha  = -*alpha;
      *beta   = -*beta;
      *gamma_  = -*gamma_;
      *delta  = -*delta;
   }

   return SCIP_OKAY;
}
