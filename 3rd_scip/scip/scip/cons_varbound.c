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

/**@file   cons_varbound.c
 * @brief  Constraint handler for variable bound constraints \f$lhs \le x + c y \le rhs\f$.
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Michael Winkler
 * @author Gerald Gamrath
 * @author Stefan Heinz
 *
 * This constraint handler handles a special type of linear constraints, namely variable bound constraints.
 * A variable bound constraint has the form
 * \f[
 *   lhs \leq x + c y \leq rhs
 * \f]
 * with coefficient \f$c \in Q\f$, \f$lhs\in Q \cup \{-\infty\}\f$, \f$rhs\in Q \cup \{\infty\}\f$,
 * and decision variables \f$x\f$ (non-binary) and \f$y\f$ (binary or integer).
 *
 * @note Although x must be non-binary when the constraint is created, it can happen that x is upgraded to a binary
 *       variable, e.g. due to aggregations or bound changes in presolving.
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "scip/cons_varbound.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"


/**@name Constraint handler properties
 *
 * @{
 */

/* constraint handler properties */
#define CONSHDLR_NAME          "varbound"
#define CONSHDLR_DESC          "variable bounds  lhs <= x + c*y <= rhs, x non-binary, y non-continuous"
#define CONSHDLR_SEPAPRIORITY   +900000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -500000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -500000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING            (SCIP_PRESOLTIMING_FAST | SCIP_PRESOLTIMING_MEDIUM)
#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "varbound"
#define EVENTHDLR_DESC         "bound change event handler for variable bound constraints"

#define LINCONSUPGD_PRIORITY     +50000 /**< priority of the constraint handler for upgrading of linear constraints */

/**@} */

/**@name Default parameter values
 *
 * @{
 */

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define DEFAULT_MAXLPCOEF         1e+09 /**< maximum coefficient in varbound constraint to be added as a row into LP */
#define DEFAULT_USEBDWIDENING      TRUE /**< should bound widening be used to initialize conflict analysis? */


#define MAXSCALEDCOEF            1000LL /**< maximal coefficient value after scaling */

/**@} */

/** variable bound constraint data */
struct SCIP_ConsData
{
   SCIP_Real             vbdcoef;            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs;                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs;                /**< right hand side of variable bound inequality */
   SCIP_VAR*             var;                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar;             /**< binary, integer or implicit integer bounding variable y */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   unsigned int          presolved:1;        /**< is the variable bound constraint already presolved? */
   unsigned int          varboundsadded:1;   /**< are the globally valid variable bounds added? */
   unsigned int          changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int          tightened:1;        /**< were the vbdcoef and all sides already tightened? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Real             maxlpcoef;          /**< maximum coefficient in varbound constraint to be added as a row into LP */
   SCIP_Bool             usebdwidening;      /**< should bound widening be used to in conflict analysis? */
};

/** Propagation rules */
enum Proprule
{
   PROPRULE_1,                          /**< left hand side and bounds on y -> lower bound on x */
   PROPRULE_2,                          /**< left hand side and upper bound on x -> bound on y */
   PROPRULE_3,                          /**< right hand side and bounds on y -> upper bound on x */
   PROPRULE_4                           /**< right hand side and lower bound on x -> bound on y */
};
typedef enum Proprule PROPRULE;


/**@name Local methods
 *
 */

/** compares two varbound constraints   cons1: \f$ lhs1 \le x1 + c1 y1 \le rhs1 \f$   and   cons2: \f$ lhs2 \le x2 + c2 y2 \le rhs2 \f$
 *  w.r.t. the indices of the contained variables
 *
 *  returns -1 if:
 *  - the index of x1 is smaller than the index of x2 or
 *  - x1 = x2 and the index of y1 is smaller than the index of y2 or
 *  - x1 = x2 and y1 = y2 and cons2 was recently changed, but cons1 not
 *
 *  returns 0 if x1 = x2, y1 = y2, and the changed status of both constraints is the same
 *
 *  and returns +1 otherwise
 */
static
SCIP_DECL_SORTPTRCOMP(consVarboundComp)
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   consdata1 = SCIPconsGetData((SCIP_CONS*) elem1);
   consdata2 = SCIPconsGetData((SCIP_CONS*) elem2);

   assert(consdata1 != NULL);
   assert(consdata2 != NULL);

   /* comparison is done over 3 ordered criteria:
    *  (i) variable index of variable 1
    *  (ii) variable index of variable 2.
    *  (iii) changed status
    */
   if( SCIPvarGetIndex(consdata1->var) < SCIPvarGetIndex(consdata2->var)
      || (SCIPvarGetIndex(consdata1->var) == SCIPvarGetIndex(consdata2->var)
         && SCIPvarGetIndex(consdata1->vbdvar) < SCIPvarGetIndex(consdata2->vbdvar))
      || (SCIPvarGetIndex(consdata1->var) == SCIPvarGetIndex(consdata2->var)
         && SCIPvarGetIndex(consdata1->vbdvar) == SCIPvarGetIndex(consdata2->vbdvar)
         && !consdata1->changed && consdata2->changed) )
      return -1;
   else if( SCIPvarGetIndex(consdata1->var) == SCIPvarGetIndex(consdata2->var)
      && SCIPvarGetIndex(consdata1->vbdvar) == SCIPvarGetIndex(consdata2->vbdvar)
      && (consdata1->changed == consdata2->changed) )
      return 0;
   else
      return +1;
}

/** creates constraint handler data for varbound constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, conshdlrdata) );

   /* set event handler for bound change events */
   (*conshdlrdata)->eventhdlr = eventhdlr;

   return SCIP_OKAY;
}

/** frees constraint handler data for varbound constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, conshdlrdata);
}

/** catches events for variables
 *
 *  @todo if lhs or rhs is infinite, catch only changes of the bound that could lead to propagation
 */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);
   assert(eventhdlr != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->var, SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, (SCIP_EVENTDATA*)cons, NULL) );
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vbdvar, SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, (SCIP_EVENTDATA*)cons, NULL) );

   return SCIP_OKAY;
}

/** drops events for variables */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   SCIP_CONSDATA* consdata;
   assert(cons != NULL);
   assert(eventhdlr != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->var, SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, (SCIP_EVENTDATA*)cons, -1) );
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vbdvar, SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_VARFIXED, eventhdlr, (SCIP_EVENTDATA*)cons, -1) );

   return SCIP_OKAY;
}

/** creates a variable bound constraint data object */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the variable bound constraint data */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs,                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs                 /**< right hand side of variable bound inequality */
   )
{
   assert(consdata != NULL);
   assert(SCIPvarGetType(vbdvar) != SCIP_VARTYPE_CONTINUOUS);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -rhs) )
      rhs = -SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, lhs) )
      lhs = SCIPinfinity(scip);

   if( SCIPisGT(scip, lhs, rhs) )
   {
      SCIPerrorMessage("left hand side of varbound constraint greater than right hand side\n");
      SCIPerrorMessage(" -> lhs=%g, rhs=%g\n", lhs, rhs);
      return SCIP_INVALIDDATA;
   }

   if( SCIPisZero(scip, vbdcoef) )
   {
      SCIPerrorMessage("varbound coefficient must be different to zero.\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPisInfinity(scip, vbdcoef) )
      vbdcoef = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -vbdcoef) )
      vbdcoef = -SCIPinfinity(scip);

   (*consdata)->var = var;
   (*consdata)->vbdvar = vbdvar;
   (*consdata)->vbdcoef = vbdcoef;
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;
   (*consdata)->row = NULL;
   (*consdata)->presolved = FALSE;
   (*consdata)->varboundsadded = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->tightened = FALSE;

   /* if we are in the transformed problem, get transformed variables, add variable bound information, and catch events */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->var, &(*consdata)->var) );
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vbdvar, &(*consdata)->vbdvar) );

#ifndef NDEBUG
      assert(SCIPvarGetStatus(SCIPvarGetProbvar((*consdata)->var)) != SCIP_VARSTATUS_MULTAGGR);
      assert(SCIPvarGetStatus(SCIPvarGetProbvar((*consdata)->vbdvar)) != SCIP_VARSTATUS_MULTAGGR);
#endif
   }

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->var) );
   SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vbdvar) );

   return SCIP_OKAY;
}

/** frees a variable bound constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the variable bound constraint */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* release variables */
   SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->var) );
   SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vbdvar) );

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** creates LP row corresponding to variable bound constraint */
static
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< variable bound constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->row, SCIPconsGetHdlr(cons), SCIPconsGetName(cons), consdata->lhs, consdata->rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->var, 1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->vbdvar, consdata->vbdcoef) );

   return SCIP_OKAY;
}

/** adds linear relaxation of variable bound constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility was detected */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* find the variable bound constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("variable bound constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   assert(SCIPvarGetType(consdata->vbdvar) != SCIP_VARTYPE_CONTINUOUS);

   /* check whether the coefficient is too large to put the row into the LP */
   if( SCIPisGT(scip, REALABS(consdata->vbdcoef), conshdlrdata->maxlpcoef) )
      return SCIP_OKAY;

   if( consdata->row == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMsg(scip, "adding relaxation of variable bound constraint <%s>: ", SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, consdata->row, NULL)) );
      SCIP_CALL( SCIPaddRow(scip, consdata->row, FALSE, infeasible) );
   }

   return SCIP_OKAY;
}

/** returns whether the given solution is feasible for the given variable bound constraint */
static
SCIP_Bool checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows         /**< Do constraints represented by rows in the current LP have to be checked? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real solval;
   SCIP_Real absviol;
   SCIP_Real relviol;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMsg(scip, "checking variable bound constraint <%s> for feasibility of solution %p (lprows=%u)\n",
      SCIPconsGetName(cons), (void*)sol, checklprows);

   solval = SCIPgetSolVal(scip, sol, consdata->var);

   if( SCIPisFeasZero(scip, SCIPgetSolVal(scip, sol, consdata->vbdvar)) && (!SCIPisFeasLE(scip, solval, consdata->rhs) || !SCIPisFeasGE(scip, solval, consdata->lhs)) )
      return FALSE;


   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      SCIP_Real sum;
      SCIP_Real lhsrelviol;
      SCIP_Real rhsrelviol;

      sum = solval + consdata->vbdcoef * SCIPgetSolVal(scip, sol, consdata->vbdvar);

      /* calculate constraint violation and update it in solution */
      absviol = MAX(consdata->lhs - sum, sum - consdata->rhs);
      lhsrelviol = SCIPrelDiff(consdata->lhs, sum);
      rhsrelviol = SCIPrelDiff(sum, consdata->rhs);
      relviol = MAX(lhsrelviol, rhsrelviol);
      if( sol != NULL )
         SCIPupdateSolLPConsViolation(scip, sol, absviol, relviol);

      return (SCIPisInfinity(scip, -consdata->lhs) || SCIPisFeasGE(scip, sum, consdata->lhs))
         && (SCIPisInfinity(scip, consdata->rhs) || SCIPisFeasLE(scip, sum, consdata->rhs));
   }
   else
      return TRUE;
}


/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) left hand side and bounds on y -> lower bound on x
 *   (2) left hand side and upper bound on x -> bound on y
 *   (3) right hand side and bounds on y -> upper bound on x
 *   (4) right hand side and lower bound on x -> bound on y
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   PROPRULE              proprule,           /**< propagation rule that deduced the bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Real             inferbd,            /**< inference bound which needs to be explained */
   SCIP_Bool             usebdwidening       /**< should bound widening be used to in conflict analysis? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* vbdvar;
   SCIP_VAR* var;
   SCIP_Real vbdcoef;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, consdata->vbdcoef));

   var = consdata->var;
   assert(var != NULL);

   vbdvar = consdata->vbdvar;
   assert(vbdvar != NULL);

   vbdcoef = consdata->vbdcoef;
   assert(!SCIPisZero(scip, vbdcoef));

   switch( proprule )
   {
   case PROPRULE_1:
      /* lhs <= x + c*y: left hand side and bounds on y -> lower bound on x */
      assert(infervar == var);
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      if( usebdwidening )
      {
         SCIP_Real relaxedbd;

         /* For integer variables, we can reduce the inferbound by 1 - z * eps, because this will be adjusted
          * to the bound we need; however, we need to choose z large enough to prevent numerical troubles due to
          * too small and too large vbdcoef values.
          * If inferbound has a large value, adding z*eps might be lost due to fixed precision floating point
          * arithmetics, so we explicitly check this here.
          */
         if( SCIPvarIsIntegral(var) && inferbd < SCIPgetHugeValue(scip) * SCIPfeastol(scip)
               && REALABS(consdata->lhs) < SCIPgetHugeValue(scip) * SCIPfeastol(scip) )
            relaxedbd = (consdata->lhs - (inferbd - 1.0 + 2.0* SCIPfeastol(scip))) / vbdcoef;
         else
            relaxedbd = (consdata->lhs - inferbd) / vbdcoef;

         /* check the computed relaxed lower/upper bound is a proper reason for the inference bound which has to be explained */
         assert(SCIPisEQ(scip, inferbd, SCIPadjustedVarLb(scip, var, consdata->lhs - relaxedbd * vbdcoef)));

         if( vbdcoef > 0.0 )
         {
            /* decrease the computed relaxed upper bound by an epsilon; this ensures that we get the actual
             * inference bound due to the integrality condition of the variable bound variable
             */
            relaxedbd -= SCIPfeastol(scip);
            SCIP_CALL( SCIPaddConflictRelaxedUb(scip, vbdvar, bdchgidx, relaxedbd) );
         }
         else
         {
            /* increase the computed relaxed lower bound by an epsilon; this ensures that we get the actual inference
             * bound due to the integrality condition of the variable bound variable
             */
            relaxedbd += SCIPfeastol(scip);
            SCIP_CALL( SCIPaddConflictRelaxedLb(scip, vbdvar, bdchgidx, relaxedbd) );
         }
      }
      else
      {
         if( vbdcoef > 0.0 )
         {
            SCIP_CALL( SCIPaddConflictUb(scip, vbdvar, bdchgidx) );
         }
         else
         {
            SCIP_CALL( SCIPaddConflictLb(scip, vbdvar, bdchgidx) );
         }
      }

      break;

   case PROPRULE_2:
      /* lhs <= x + c*y: left hand side and upper bound on x -> bound on y */
      assert(infervar == vbdvar);
      assert(SCIPvarGetType(vbdvar) != SCIP_VARTYPE_CONTINUOUS);
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      if( usebdwidening )
      {
         SCIP_Real relaxedub;

         /* compute the relaxed upper bound of the variable which would be sufficient to reach one less (greater) than the
          * inference bound
          */
         if( vbdcoef > 0.0 )
         {
            /* For integer variables, we can reduce the inferbound by 1-z*eps, because this will be adjusted
             * to the bound we need; however, we need to choose z large enough to prevent numerical troubles due to
             * too small and too large vbdcoef values.
             * If inferbound has a large value, adding z*eps might be lost due to fixed precision floating point
             * arithmetics, so we explicitly check this here.
             */
            if( SCIPvarIsIntegral(var) && inferbd < SCIPgetHugeValue(scip) * SCIPfeastol(scip)
                  && REALABS(consdata->rhs) < SCIPgetHugeValue(scip) * SCIPfeastol(scip) )
               relaxedub = consdata->lhs - (inferbd - 1.0 + 2.0 * SCIPfeastol(scip)) * vbdcoef;
            else
               relaxedub = consdata->lhs - inferbd * vbdcoef;

            /* check the computed relaxed upper bound is a proper reason for the inference bound which has to be explained */
            assert(SCIPisEQ(scip, inferbd, SCIPadjustedVarLb(scip, vbdvar, (consdata->lhs - relaxedub) / vbdcoef)));
         }
         else
         {
            /* For integer variables, we can reduce the inferbound by 1-z*eps, because this will be adjusted
             * to the bound we need; however, we need to choose z large enough to prevent numerical troubles due to
             * too small and too large vbdcoef values.
             * If inferbound has a large value, adding z*eps might be lost due to fixed precision floating point
             * arithmetics, so we explicitly check this here.
             */
            if( SCIPvarIsIntegral(var) && inferbd < SCIPgetHugeValue(scip) * SCIPfeastol(scip)
                  && REALABS(consdata->lhs) < SCIPgetHugeValue(scip) * SCIPfeastol(scip) )
               relaxedub = consdata->lhs - (inferbd + 1.0 - 2.0 * SCIPfeastol(scip)) * vbdcoef;
            else
               relaxedub = consdata->lhs - inferbd * vbdcoef;

            /* check the computed relaxed upper bound is a proper reason for the inference bound which has to be explained */
            assert(SCIPisEQ(scip, inferbd, SCIPadjustedVarUb(scip, vbdvar, (consdata->lhs - relaxedub) / vbdcoef)));
         }

         /* decrease the computed relaxed upper bound by an epsilon; this ensures that we get the actual inference bound due
          * to the integrality condition of the variable bound variable
          */
         relaxedub -= SCIPfeastol(scip);

         SCIP_CALL( SCIPaddConflictRelaxedUb(scip, var, bdchgidx, relaxedub) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
      }

      break;

   case PROPRULE_3:
      /* x + c*y <= rhs: right hand side and bounds on y -> upper bound on x */
      assert(infervar == var);
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      assert(!SCIPisInfinity(scip, consdata->rhs));

      if( usebdwidening )
      {
         SCIP_Real relaxedbd;

         /* For integer variables, we can reduce the inferbound by 1-z*eps, because this will be adjusted
          * to the bound we need; however, we need to choose z large enough to prevent numerical troubles due to
          * too small and too large vbdcoef values.
          * If inferbound has a large value, adding z*eps might be lost due to fixed precision floating point
          * arithmetics, so we explicitly check this here.
          */
         if( SCIPvarIsIntegral(var) && inferbd < SCIPgetHugeValue(scip) * SCIPfeastol(scip)
               && REALABS(consdata->rhs) < SCIPgetHugeValue(scip) * SCIPfeastol(scip) )
            relaxedbd = (consdata->rhs - (inferbd + 1.0 - 2.0 * SCIPfeastol(scip))) / vbdcoef;
         else
            relaxedbd = (consdata->rhs - inferbd) / vbdcoef;

         /* check the computed relaxed lower/upper bound is a proper reason for the inference bound which has to be explained */
         assert(SCIPisEQ(scip, inferbd, SCIPadjustedVarUb(scip, var, consdata->rhs - relaxedbd * vbdcoef)));

         if( vbdcoef > 0.0 )
         {
            /* increase the computed relaxed lower bound by an epsilon; this ensures that we get the actual inference bound due
             * to the integrality condition of the variable bound variable
             */
            relaxedbd += SCIPfeastol(scip);
            SCIP_CALL( SCIPaddConflictRelaxedLb(scip, vbdvar, bdchgidx, relaxedbd) );
         }
         else
         {
            /* decrease the computed relaxed upper bound by an epsilon; this ensures that we get the actual inference bound due
             * to the integrality condition of the variable bound variable
             */
            relaxedbd -= SCIPfeastol(scip);
            SCIP_CALL( SCIPaddConflictRelaxedUb(scip, vbdvar, bdchgidx, relaxedbd) );
         }
      }
      else
      {
         if( vbdcoef > 0.0 )
         {
            SCIP_CALL( SCIPaddConflictLb(scip, vbdvar, bdchgidx) );
         }
         else
         {
            SCIP_CALL( SCIPaddConflictUb(scip, vbdvar, bdchgidx) );
         }
      }

      break;

   case PROPRULE_4:
      /* x + c*y <= rhs: right hand side and lower bound on x -> bound on y */
      assert(infervar == vbdvar);
      assert(SCIPvarGetType(vbdvar) != SCIP_VARTYPE_CONTINUOUS);
      assert(!SCIPisInfinity(scip, consdata->rhs));

      if( usebdwidening )
      {
         SCIP_Real relaxedlb;

         /* compute the relaxed lower bound of the variable which would be sufficient to reach one greater (less) than the
          * inference bound
          */
         if( vbdcoef > 0.0 )
         {
            /* For integer variables, we can reduce the inferbound by 1-z*eps, because this will be adjusted
             * to the bound we need; however, we need to choose z large enough to prevent numerical troubles due to
             * too small and too large vbdcoef values.
             * If inferbound has a large value, adding z*eps might be lost due to fixed precision floating point
             * arithmetics, so we explicitly check this here.
             */
            if( SCIPvarIsIntegral(var) && inferbd < SCIPgetHugeValue(scip) * SCIPfeastol(scip)
                  && REALABS(consdata->rhs) < SCIPgetHugeValue(scip) * SCIPfeastol(scip) )
               relaxedlb = consdata->rhs - (inferbd + 1.0 - 2.0 * SCIPfeastol(scip)) * vbdcoef;
            else
               relaxedlb = consdata->rhs - inferbd * vbdcoef;

            /* check the computed relaxed lower bound is a proper reason for the inference bound which has to be explained */
            assert(SCIPisEQ(scip, inferbd, SCIPadjustedVarUb(scip, vbdvar, (consdata->rhs - relaxedlb) / vbdcoef)));
         }
         else
         {
            /* For integer variables, we can reduce the inferbound by 1-z*eps, because this will be adjusted
             * to the bound we need; however, we need to choose z large enough to prevent numerical troubles due to
             * too small and too large vbdcoef values.
             * If inferbound has a large value, adding z*eps might be lost due to fixed precision floating point
             * arithmetics, so we explicitly check this here.
             */
            if( SCIPvarIsIntegral(var) && inferbd < SCIPgetHugeValue(scip) * SCIPfeastol(scip)
                  && REALABS(consdata->lhs) < SCIPgetHugeValue(scip) * SCIPfeastol(scip) )
               relaxedlb = consdata->rhs - (inferbd - 1.0 + 2.0 * SCIPfeastol(scip)) * vbdcoef;
            else
               relaxedlb = consdata->rhs - inferbd * vbdcoef;

            /* check the computed relaxed lower bound is a proper reason for the inference bound which has to be explained */
            assert(SCIPisEQ(scip, inferbd, SCIPadjustedVarLb(scip, vbdvar, (consdata->rhs - relaxedlb) / vbdcoef)));
         }

         /* increase the computed relaxed lower bound by an epsilon; this ensures that we get the actual inference bound due
          * to the integrality condition of the variable bound variable
          */
         relaxedlb += SCIPfeastol(scip);

         SCIP_CALL( SCIPaddConflictRelaxedLb(scip, var, bdchgidx, relaxedlb) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
      }

      break;

   default:
      SCIPerrorMessage("invalid inference information %d in variable bound constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** analyze infeasibility */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   SCIP_Real             inferbd,            /**< bound which led to infeasibility */
   PROPRULE              proprule,           /**< propagation rule that deduced the bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_Bool             usebdwidening       /**< should bound widening be used to in conflict analysis? */
   )
{
   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   /* add the bound which got violated */
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      if( usebdwidening )
      {
         SCIP_Real relaxedub;

         /* adjust lower bound */
         inferbd = SCIPadjustedVarLb(scip, infervar, inferbd);

         /* compute a relaxed upper bound which would be sufficient to be still infeasible */
         if( SCIPvarIsIntegral(infervar) )
            relaxedub = inferbd - 1.0;
         else
         {
            SCIP_CONSDATA* consdata;
            SCIP_Real abscoef;

            consdata = SCIPconsGetData(cons);
            assert(consdata != NULL);

            /* vbdvar can never be of non-integral type */
            assert(infervar == consdata->var);

            abscoef = REALABS(consdata->vbdcoef);

            /* due to resolving a the propagation and dividing by the vbdcoef we need to make sure the the relaxed bound
             * is big enough, therefore we multiply here with the vbdcoef
             *
             * @note it does not matter if we deceed the current local upper bound, because SCIPaddConflictRelaxedUb()
             *       is correcting the bound afterwards
             */
            relaxedub = inferbd - 2*SCIPfeastol(scip) * MAX(1, abscoef);
         }

         /* try to relax inference variable upper bound such that the infeasibility is still given */
         SCIP_CALL( SCIPaddConflictRelaxedUb(scip, infervar, NULL, relaxedub) );

         /* collect the upper bound which is reported to the conflict analysis */
         inferbd = SCIPgetConflictVarUb(scip, infervar);

         /* adjust inference bound with respect to the upper bound reported to the conflict analysis */
         if( SCIPvarIsIntegral(infervar) )
            inferbd = inferbd + 1.0;
         else
         {
            SCIP_CONSDATA* consdata;
            SCIP_Real abscoef;

            consdata = SCIPconsGetData(cons);
            assert(consdata != NULL);

            /* vbdvar can never be of non-integral type */
            assert(infervar == consdata->var);

            abscoef = REALABS(consdata->vbdcoef);

            /* due to resolving a the propagation and dividing by the vbdcoef we need to make sure the the relaxed bound
             * is big enough, therefore we multiply here with the vbdcoef
             */
            inferbd = inferbd + 2*SCIPfeastol(scip) * MAX(1, abscoef);
         }
      }
      else
      {
         SCIP_CALL( SCIPaddConflictUb(scip, infervar, NULL) );
      }
   }
   else
   {
      if( usebdwidening )
      {
         SCIP_Real relaxedlb;

         assert(boundtype == SCIP_BOUNDTYPE_UPPER);

         /* adjust upper bound */
         inferbd = SCIPadjustedVarUb(scip, infervar, inferbd);

         /* compute a relaxed lower bound which would be sufficient to be still infeasible */
         if( SCIPvarIsIntegral(infervar) )
            relaxedlb = inferbd + 1.0;
         else
         {
            SCIP_CONSDATA* consdata;
            SCIP_Real abscoef;

            consdata = SCIPconsGetData(cons);
            assert(consdata != NULL);

            /* vbdvar can never be of non-integral type */
            assert(infervar == consdata->var);

            abscoef = REALABS(consdata->vbdcoef);

            /* due to resolving a the propagation and dividing by the vbdcoef we need to make sure the the relaxed bound
             * is big enough, therefore we multiply here with the vbdcoef
             *
             * @note it does not matter if we exceed the current local lower bound, because SCIPaddConflictRelaxedLb()
             *       is correcting the bound afterwards
             */
            relaxedlb = inferbd + 2*SCIPfeastol(scip) * MAX(1, abscoef);
         }

         /* try to relax inference variable upper bound such that the infeasibility is still given */
         SCIP_CALL( SCIPaddConflictRelaxedLb(scip, infervar, NULL, relaxedlb) );

         /* collect the lower bound which is reported to the conflict analysis */
         inferbd = SCIPgetConflictVarLb(scip, infervar);

         /* adjust inference bound with respect to the lower bound reported to the conflict analysis */
         if( SCIPvarIsIntegral(infervar) )
            inferbd = inferbd - 1.0;
         else
         {
            SCIP_CONSDATA* consdata;
            SCIP_Real abscoef;

            consdata = SCIPconsGetData(cons);
            assert(consdata != NULL);

            /* vbdvar can never be of non-integral type */
            assert(infervar == consdata->var);

            abscoef = REALABS(consdata->vbdcoef);

            /* due to resolving a the propagation and dividing by the vbdcoef we need to make sure the the relaxed bound
             * is big enough, therefore we multiply here with the vbdcoef
             */
            inferbd = inferbd - 2*SCIPfeastol(scip) * MAX(1, abscoef);
         }
      }
      else
      {
         SCIP_CALL( SCIPaddConflictLb(scip, infervar, NULL) );
      }
   }

   /* add the reason for the violated of the bound */
   SCIP_CALL( resolvePropagation(scip, cons, infervar, proprule, boundtype, NULL, inferbd, usebdwidening) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** separates the given variable bound constraint */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_Bool             usebdwidening,      /**< should bound widening be used to in conflict analysis? */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* vbdvar;
   SCIP_VAR* var;
   SCIP_Real vbdcoef;
   SCIP_Real feasibility;

   assert(cons != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* find the variable bound constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("variable bound constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMsg(scip, "separating variable bound constraint <%s>\n", SCIPconsGetName(cons));

   var = consdata->var;
   vbdvar = consdata->vbdvar;
   vbdcoef = consdata->vbdcoef;
   assert(SCIPvarGetType(vbdvar) != SCIP_VARTYPE_CONTINUOUS);

   if( SCIPvarGetLbLocal(vbdvar) + 0.5 > SCIPvarGetUbLocal(vbdvar) )
   {
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vbdvar), SCIPvarGetUbLocal(vbdvar)));

      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         SCIP_Real newlb;
         SCIP_Bool cutoff;
         SCIP_Bool tightened;

         newlb = consdata->lhs - vbdcoef * SCIPvarGetLbLocal(vbdvar);

         SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, (int)PROPRULE_1, TRUE,
               &cutoff, &tightened) );

         if( cutoff )
         {
            assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(var)));

            /* analyze infeasibility */
            SCIP_CALL( analyzeConflict(scip, cons, var, newlb, PROPRULE_1, SCIP_BOUNDTYPE_LOWER, usebdwidening) );
            *result = SCIP_CUTOFF;

            return SCIP_OKAY;
         }
         else if( tightened )
         {
            *result = SCIP_REDUCEDDOM;
         }
      }

      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         SCIP_Real newub;
         SCIP_Bool cutoff;
         SCIP_Bool tightened;

         newub = consdata->rhs - vbdcoef * SCIPvarGetLbLocal(vbdvar);

         SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, (int)PROPRULE_3, TRUE,
               &cutoff, &tightened) );

         if( cutoff )
         {
            assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(var)));

            /* analyze infeasibility */
            SCIP_CALL( analyzeConflict(scip, cons, var, newub, PROPRULE_3, SCIP_BOUNDTYPE_UPPER, usebdwidening) );
            *result = SCIP_CUTOFF;

            return SCIP_OKAY;
         }
         else if( tightened )
         {
            *result = SCIP_REDUCEDDOM;
         }
      }
   }

   /* if we already changed a bound or the coefficient is too large to put the row into the LP, stop here */
   if( *result == SCIP_REDUCEDDOM )
      return SCIP_OKAY;

   /* check constraint for feasibility and create row if constraint is violated */
   if( !checkCons(scip, cons, sol, (sol != NULL)) )
   {
      /* create LP relaxation if not yet existing */
      if( consdata->row == NULL )
      {
         SCIP_CALL( createRelaxation(scip, cons) );
      }
      assert(consdata->row != NULL);

      /* check non-LP rows for feasibility and add them as cut, if violated */
      if( !SCIProwIsInLP(consdata->row) )
      {
         feasibility = SCIPgetRowSolFeasibility(scip, consdata->row, sol);
         if( SCIPisFeasNegative(scip, feasibility) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPaddRow(scip, consdata->row, FALSE, &infeasible) );
            if ( infeasible )
               *result = SCIP_CUTOFF;
            else
               *result = SCIP_SEPARATED;
         }
      }
   }

   return SCIP_OKAY;
}

/** sets left hand side of varbound constraint */
static
SCIP_RETCODE chgLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, lhs));

   /* adjust value to not be smaller than -inf */
   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->var != NULL && consdata->vbdvar != NULL);
   assert(!SCIPisInfinity(scip, consdata->lhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->lhs, lhs) )
      return SCIP_OKAY;

   assert(consdata->row == NULL);

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */
   if( SCIPisEQ(scip, lhs, consdata->rhs) )
      consdata->rhs = lhs;

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      /* the left hand side switched from -infinity to a non-infinite value -> install rounding locks */
      if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, -lhs) )
      {
	 SCIP_CALL( SCIPlockVarCons(scip, consdata->var, cons, TRUE, FALSE) );

	 if( SCIPisPositive(scip, consdata->vbdcoef) )
	 {
	    SCIP_CALL( SCIPlockVarCons(scip, consdata->vbdvar, cons, TRUE, FALSE) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPlockVarCons(scip, consdata->vbdvar, cons, FALSE, TRUE) );
	 }
      }
      /* the left hand side switched from a non-infinite value to -infinity -> remove rounding locks */
      else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, -lhs) )
      {
	 SCIP_CALL( SCIPunlockVarCons(scip, consdata->var, cons, TRUE, FALSE) );
	 if( SCIPisPositive(scip, consdata->vbdcoef) )
	 {
	    SCIP_CALL( SCIPunlockVarCons(scip, consdata->vbdvar, cons, TRUE, FALSE) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPunlockVarCons(scip, consdata->vbdvar, cons, FALSE, TRUE) );
	 }
      }
   }

   /* if left hand side got tighter, we want to do additional presolving on this constraint */
   if( SCIPisLT(scip, consdata->lhs, lhs) )
   {
      consdata->varboundsadded = FALSE;
      consdata->tightened = FALSE;

      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   consdata->presolved = FALSE;
   consdata->lhs = lhs;
   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** sets right hand side of varbound constraint */
static
SCIP_RETCODE chgRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, -rhs));

   /* adjust value to not be larger than inf */
   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->var != NULL && consdata->vbdvar != NULL);
   assert(!SCIPisInfinity(scip, -consdata->rhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->rhs, rhs) )
      return SCIP_OKAY;

   assert(consdata->row == NULL);

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */
   if( SCIPisEQ(scip, rhs, consdata->lhs) )
      consdata->lhs = rhs;

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
      if( SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, rhs) )
      {
	 SCIP_CALL( SCIPlockVarCons(scip, consdata->var, cons, FALSE, TRUE) );

	 if( SCIPisPositive(scip, consdata->vbdcoef) )
	 {
	    SCIP_CALL( SCIPlockVarCons(scip, consdata->vbdvar, cons, FALSE, TRUE) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPlockVarCons(scip, consdata->vbdvar, cons, TRUE, FALSE) );
	 }
      }
      /* the right hand side switched from a non-infinite value to infinity -> remove rounding locks */
      else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, rhs) )
      {
	 SCIP_CALL( SCIPunlockVarCons(scip, consdata->var, cons, FALSE, TRUE) );
	 if( SCIPisPositive(scip, consdata->vbdcoef) )
	 {
	    SCIP_CALL( SCIPunlockVarCons(scip, consdata->vbdvar, cons, FALSE, TRUE) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPunlockVarCons(scip, consdata->vbdvar, cons, TRUE, FALSE) );
	 }
      }
   }

   /* if right hand side got tighter, we want to do additional presolving on this constraint */
   if( SCIPisGT(scip, consdata->rhs, rhs) )
   {
      consdata->varboundsadded = FALSE;
      consdata->tightened = FALSE;

      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   consdata->presolved = FALSE;
   consdata->rhs = rhs;
   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** propagation method for variable bound constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_Bool             usebdwidening,      /**< should bound widening be used to in conflict analysis? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  nchgsides,          /**< pointer to count number of side changes */
   int*                  ndelconss           /**< pointer to count number of deleted constraints, or NULL */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real ylb;
   SCIP_Real yub;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Bool tightened;
   SCIP_Bool tightenedround;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMsg(scip, "propagating variable bound constraint <%s>: %.15g <= <%s>[%.9g, %.9g] + %.15g<%s>[%.9g, %.9g] <= %.15g\n",
      SCIPconsGetName(cons), consdata->lhs, SCIPvarGetName(consdata->var), SCIPvarGetLbLocal(consdata->var),
      SCIPvarGetUbLocal(consdata->var), consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar),
      SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar), consdata->rhs);

   *cutoff = FALSE;

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   /* get current bounds of variables */
   xlb = SCIPvarGetLbLocal(consdata->var);
   xub = SCIPvarGetUbLocal(consdata->var);
   ylb = SCIPvarGetLbLocal(consdata->vbdvar);
   yub = SCIPvarGetUbLocal(consdata->vbdvar);

   /* it can happen that constraint is of form lhs <= x <= rhs */
   if( SCIPisZero(scip, consdata->vbdcoef) && SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
   {
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      SCIP_CALL( SCIPfixVar(scip, consdata->var, consdata->lhs, &infeasible, &fixed) );

      if( infeasible )
      {
         SCIPdebugMsg(scip, "> constraint <%s> is infeasible.\n", SCIPconsGetName(cons));
         *cutoff = TRUE;
         return SCIP_OKAY;
      }
   }

   /* tighten bounds of variables as long as possible */
   do
   {
      tightenedround = FALSE;

      /* propagate left hand side inequality: lhs <= x + c*y */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         assert(!(*cutoff));

         /* propagate bounds on x:
          *  (1) left hand side and bounds on y -> lower bound on x
          */
         if( SCIPvarGetStatus(consdata->var) != SCIP_VARSTATUS_MULTAGGR ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               if( !SCIPisInfinity(scip, yub) )
                  newlb = SCIPadjustedVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * yub);
               else
                  newlb = -SCIPinfinity(scip);
            }
            else
            {
               if( !SCIPisInfinity(scip, -ylb) )
                  newlb = SCIPadjustedVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * ylb);
               else
                  newlb = -SCIPinfinity(scip);
            }

            SCIP_CALL( SCIPinferVarLbCons(scip, consdata->var, newlb, cons, (int)PROPRULE_1, yub < ylb + 0.5, cutoff, &tightened) );

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "cutoff while tightening <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->var), xlb, xub, newlb, xub);
               assert( SCIPisInfinity(scip, newlb) || SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->var)) );

               SCIP_CALL( SCIPresetConsAge(scip, cons) );

               /* analyze infeasibility */
               SCIP_CALL( analyzeConflict(scip, cons, consdata->var, newlb, PROPRULE_1, SCIP_BOUNDTYPE_LOWER, usebdwidening) );
               break;
            }

            if( tightened )
            {
               SCIPdebugMsg(scip, " -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->var), xlb, xub, newlb, xub);
               tightenedround = TRUE;
               (*nchgbds)++;
               SCIP_CALL( SCIPresetConsAge(scip, cons) );
            }
            xlb = SCIPvarGetLbLocal(consdata->var);
         }

         assert(!*cutoff);

         /* propagate bounds on y:
          *  (2) left hand side and upper bound on x -> bound on y
          */
         if( SCIPvarGetStatus(consdata->vbdvar) != SCIP_VARSTATUS_MULTAGGR && !SCIPisInfinity(scip, xub) ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               newlb = SCIPadjustedVarLb(scip, consdata->vbdvar, (consdata->lhs - xub)/consdata->vbdcoef);
               if( newlb > ylb + 0.5 )
               {
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vbdvar, newlb, cons, (int)PROPRULE_2, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     SCIPdebugMsg(scip, "cutoff while tightening <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
                     assert( SCIPisInfinity(scip, newlb) || SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->vbdvar)) );

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, newlb, PROPRULE_2, SCIP_BOUNDTYPE_LOWER, usebdwidening) );
                     break;
                  }

                  if( tightened )
                  {
                     SCIPdebugMsg(scip, " -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  ylb = SCIPvarGetLbLocal(consdata->vbdvar);
               }
            }
            else
            {
               newub = SCIPadjustedVarUb(scip, consdata->vbdvar, (consdata->lhs - xub)/consdata->vbdcoef);
               if( newub < yub - 0.5 )
               {
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->vbdvar, newub, cons, (int)PROPRULE_2, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     SCIPdebugMsg(scip, "cutoff while tightening <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
                     assert( SCIPisInfinity(scip, -newub) || SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->vbdvar)) );

                     SCIP_CALL( SCIPresetConsAge(scip, cons) );

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, newub, PROPRULE_2, SCIP_BOUNDTYPE_UPPER, usebdwidening) );
                     break;
                  }

                  if( tightened )
                  {
                     SCIPdebugMsg(scip, " -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
                     tightenedround = TRUE;
                     (*nchgbds)++;
                     SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  }
                  yub = SCIPvarGetUbLocal(consdata->vbdvar);
               }
            }
         }
      }

      assert(!*cutoff);

      /* propagate right hand side inequality: x + c*y <= rhs */
      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         /* propagate bounds on x:
          *  (3) right hand side and bounds on y -> upper bound on x
          */
         if( SCIPvarGetStatus(consdata->var) != SCIP_VARSTATUS_MULTAGGR ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               if( !SCIPisInfinity(scip, -ylb) )
                  newub = SCIPadjustedVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * ylb);
               else
                  newub = SCIPinfinity(scip);
            }
            else
            {
               if( !SCIPisInfinity(scip, yub) )
                  newub = SCIPadjustedVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * yub);
               else
                  newub = SCIPinfinity(scip);
            }

            SCIP_CALL( SCIPinferVarUbCons(scip, consdata->var, newub, cons, (int)PROPRULE_3, yub < ylb + 0.5, cutoff, &tightened) );

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "cutoff while tightening <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->var), xlb, xub, xlb, newub);
               assert( SCIPisInfinity(scip, -newub) || SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->var)) );

               SCIP_CALL( SCIPresetConsAge(scip, cons) );

               /* analyze infeasibility */
               SCIP_CALL( analyzeConflict(scip, cons, consdata->var, newub, PROPRULE_3, SCIP_BOUNDTYPE_UPPER, usebdwidening) );
               break;
            }

            if( tightened )
            {
               SCIPdebugMsg(scip, " -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->var), xlb, xub, xlb, newub);
               tightenedround = TRUE;
               (*nchgbds)++;
               SCIP_CALL( SCIPresetConsAge(scip, cons) );
            }
            xub = SCIPvarGetUbLocal(consdata->var);
         }

         assert(!*cutoff);

         /* propagate bounds on y:
          *  (4) right hand side and lower bound on x -> bound on y
          */
         if( SCIPvarGetStatus(consdata->vbdvar) != SCIP_VARSTATUS_MULTAGGR && !SCIPisInfinity(scip, -xlb) ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               newub = SCIPadjustedVarUb(scip, consdata->vbdvar, (consdata->rhs - xlb)/consdata->vbdcoef);
               if( newub < yub - 0.5 )
               {
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->vbdvar, newub, cons, (int)PROPRULE_4, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     SCIPdebugMsg(scip, "cutoff while tightening <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
                     assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->vbdvar)));

                     SCIP_CALL( SCIPresetConsAge(scip, cons) );

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, newub, PROPRULE_4, SCIP_BOUNDTYPE_UPPER, usebdwidening) );
                     break;
                  }

                  if( tightened )
                  {
                     SCIPdebugMsg(scip, " -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
                     tightenedround = TRUE;
                     (*nchgbds)++;
                     SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  }
                  yub = SCIPvarGetUbLocal(consdata->vbdvar);
               }
            }
            else
            {
               newlb = SCIPadjustedVarLb(scip, consdata->vbdvar, (consdata->rhs - xlb)/consdata->vbdcoef);
               if( newlb > ylb + 0.5 )
               {
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vbdvar, newlb, cons, (int)PROPRULE_4, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     SCIPdebugMsg(scip, "cutoff while tightening <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
                     assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->vbdvar)));

                     SCIP_CALL( SCIPresetConsAge(scip, cons) );

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, newlb, PROPRULE_4, SCIP_BOUNDTYPE_LOWER, usebdwidening) );
                     break;
                  }

                  if( tightened )
                  {
                     SCIPdebugMsg(scip, " -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n", SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
                     tightenedround = TRUE;
                     (*nchgbds)++;
                     SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  }
                  ylb = SCIPvarGetLbLocal(consdata->vbdvar);
               }
            }
         }
      }
      assert(!(*cutoff));
   }
   while( tightenedround );

   /* check for redundant sides */
   if( !(*cutoff) && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && !SCIPinProbing(scip) )
   {
      /* check left hand side for redundancy */
      if( !SCIPisInfinity(scip, -consdata->lhs) &&
         ((consdata->vbdcoef > 0.0 && SCIPisFeasGE(scip, xlb + consdata->vbdcoef * ylb, consdata->lhs))
            || (consdata->vbdcoef < 0.0 && SCIPisFeasGE(scip, xlb + consdata->vbdcoef * yub, consdata->lhs))) )
      {
         SCIPdebugMsg(scip, "left hand side of variable bound constraint <%s> is redundant\n", SCIPconsGetName(cons));

         SCIP_CALL( chgLhs(scip, cons, -SCIPinfinity(scip)) );
         ++(*nchgsides);
      }

      /* check right hand side for redundancy */
      if( !SCIPisInfinity(scip, consdata->rhs) &&
         ((consdata->vbdcoef > 0.0 && SCIPisFeasLE(scip, xub + consdata->vbdcoef * yub, consdata->rhs))
            || (consdata->vbdcoef < 0.0 && SCIPisFeasLE(scip, xub + consdata->vbdcoef * ylb, consdata->rhs))) )
      {
         SCIPdebugMsg(scip, "right hand side of variable bound constraint <%s> is redundant\n", SCIPconsGetName(cons));

         SCIP_CALL( chgRhs(scip, cons, SCIPinfinity(scip)) );
         ++(*nchgsides);
      }
   }
   /* check varbound constraint for redundancy */
   if( !(*cutoff) && (SCIPisInfinity(scip, -consdata->lhs)
         || (consdata->vbdcoef > 0.0 && SCIPisFeasGE(scip, xlb + consdata->vbdcoef * ylb, consdata->lhs))
         || (consdata->vbdcoef < 0.0 && SCIPisFeasGE(scip, xlb + consdata->vbdcoef * yub, consdata->lhs)))
      && (SCIPisInfinity(scip, consdata->rhs)
         || (consdata->vbdcoef > 0.0 && SCIPisFeasLE(scip, xub + consdata->vbdcoef * yub, consdata->rhs))
         || (consdata->vbdcoef < 0.0 && SCIPisFeasLE(scip, xub + consdata->vbdcoef * ylb, consdata->rhs))) )
   {
      SCIPdebugMsg(scip, "variable bound constraint <%s> is redundant: <%s>[%.15g,%.15g], <%s>[%.15g,%.15g]\n",
         SCIPconsGetName(cons),
         SCIPvarGetName(consdata->var), SCIPvarGetLbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var),
         SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar));
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

      /* this did not seem to help but should be tested again, there might also still be a bug in there */
#ifdef SCIP_DISABLED_CODE
      /* local duality fixing of variables in the constraint */
      if( !SCIPisNegative(scip, SCIPvarGetObj(consdata->vbdvar)) && SCIPvarGetNLocksDown(consdata->vbdvar) == 1
         && !SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->vbdvar))
         && SCIPisFeasLT(scip, SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar))
         && ((consdata->vbdcoef > 0.0 && !SCIPisInfinity(scip, -consdata->lhs))
            || (consdata->vbdcoef < 0.0 && !SCIPisInfinity(scip, consdata->rhs))) )
      {
         SCIPdebugMsg(scip, " --> fixing <%s>[%.15g,%.15g] to %.15g\n", SCIPvarGetName(consdata->vbdvar),
            SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar), SCIPvarGetLbLocal(consdata->vbdvar));
         SCIP_CALL( SCIPchgVarUb(scip, consdata->vbdvar, SCIPvarGetLbLocal(consdata->vbdvar)) );
      }
      else if( !SCIPisPositive(scip, SCIPvarGetObj(consdata->vbdvar)) && SCIPvarGetNLocksUp(consdata->vbdvar) == 1
         && !SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->vbdvar))
         && SCIPisFeasLT(scip, SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar))
         && ((consdata->vbdcoef < 0.0 && !SCIPisInfinity(scip, -consdata->lhs))
            || (consdata->vbdcoef > 0.0 && !SCIPisInfinity(scip, consdata->rhs))) )
      {
         SCIPdebugMsg(scip, " --> fixing <%s>[%.15g,%.15g] to %.15g\n", SCIPvarGetName(consdata->vbdvar),
            SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar));
         SCIP_CALL( SCIPchgVarLb(scip, consdata->vbdvar, SCIPvarGetUbLocal(consdata->vbdvar)) );
      }
      if( !SCIPisNegative(scip, SCIPvarGetObj(consdata->var)) && SCIPvarGetNLocksDown(consdata->var) == 1
         && !SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->var))
         && SCIPisFeasLT(scip, SCIPvarGetLbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var))
         && !SCIPisInfinity(scip, -consdata->lhs) )
      {
         SCIPdebugMsg(scip, " --> fixing <%s>[%.15g,%.15g] to %.15g\n", SCIPvarGetName(consdata->var),
            SCIPvarGetLbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var), SCIPvarGetLbLocal(consdata->var));
         SCIP_CALL( SCIPchgVarUb(scip, consdata->var, SCIPvarGetLbLocal(consdata->var)) );
      }
      else if( !SCIPisPositive(scip, SCIPvarGetObj(consdata->var)) && SCIPvarGetNLocksUp(consdata->var) == 1
         && !SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->var))
         && SCIPisFeasLT(scip, SCIPvarGetLbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var))
         && !SCIPisInfinity(scip, consdata->rhs) )
      {
         SCIPdebugMsg(scip, " --> fixing <%s>[%.15g,%.15g] to %.15g\n", SCIPvarGetName(consdata->var),
            SCIPvarGetLbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var));
         SCIP_CALL( SCIPchgVarLb(scip, consdata->var, SCIPvarGetUbLocal(consdata->var)) );
      }
#endif
      if( ndelconss != NULL )
         (*ndelconss)++;
   }

   SCIP_CALL( SCIPunmarkConsPropagate(scip, cons) );

   return SCIP_OKAY;
}

/* check whether one constraints side is redundant to another constraints side by calculating extreme values for
 * variables
 */
static
void checkRedundancySide(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             coef0,              /**< coefficient c0 of bounding variable y for constraint 0 */
   SCIP_Real             coef1,              /**< coefficient c1 of bounding variable y for constraint 1 */
   SCIP_Real             side0,              /**< one side of variable bound inequality for constraint 0 */
   SCIP_Real             side1,              /**< one side of variable bound inequality for constraint 1 */
   SCIP_Bool*            sideequal,          /**< pointer to store if both constraints have the same redundancy on the
                                              *   given side */
   SCIP_Bool*            cons0sidered,       /**< pointer to store if side of constraint 0 is redundant */
   SCIP_Bool*            cons1sidered,       /**< pointer to store if side of constraint 1 is redundant */
   SCIP_Bool             islhs               /**< do we check the left or the right hand side */
   )
{
   SCIP_Real lbvar;
   SCIP_Real ubvar;
   SCIP_Real lbvbdvar;
   SCIP_Real ubvbdvar;
   SCIP_Real boundxlb1;
   SCIP_Real boundxlb2;
   SCIP_Real boundylb1;
   SCIP_Real boundylb2;
   SCIP_Real boundxub1;
   SCIP_Real boundxub2;
   SCIP_Real boundyub1;
   SCIP_Real boundyub2;
   SCIP_Real boundvaluex1;
   SCIP_Real boundvaluex2;
   SCIP_Real boundvaluey1;
   SCIP_Real boundvaluey2;
   SCIP_Real valuex1;
   SCIP_Real valuex2;
   SCIP_Real valuey1;
   SCIP_Real valuey2;
   SCIP_Bool* redundant0;
   SCIP_Bool* redundant1;
   SCIP_Real eps = SCIPepsilon(scip);

   assert(scip != NULL);
   assert(var != NULL);
   assert(vbdvar != NULL);
   assert(sideequal != NULL);
   assert(cons0sidered != NULL);
   assert(cons1sidered != NULL);

   *cons0sidered = SCIPisInfinity(scip, REALABS(side0));
   *cons1sidered = SCIPisInfinity(scip, REALABS(side1));
   *sideequal = FALSE;

   if( islhs )
   {
      redundant0 = cons1sidered;
      redundant1 = cons0sidered;
   }
   else
   {
      redundant0 = cons0sidered;
      redundant1 = cons1sidered;
   }

   lbvar = SCIPvarGetLbGlobal(var);
   ubvar = SCIPvarGetUbGlobal(var);
   lbvbdvar = SCIPvarGetLbGlobal(vbdvar);
   ubvbdvar = SCIPvarGetUbGlobal(vbdvar);

   /* if both constraint have this side */
   if( !*redundant0 && !*redundant1 )
   {
      /* calculate extreme values, which are reached by setting the other variable to their lower/upper bound */
      boundxlb1 = side0 - lbvbdvar*coef0;
      boundxlb2 = side1 - lbvbdvar*coef1;
      boundylb1 = (side0 - lbvar)/coef0;
      boundylb2 = (side1 - lbvar)/coef1;

      boundxub1 = side0 - ubvbdvar*coef0;
      boundxub2 = side1 - ubvbdvar*coef1;
      boundyub1 = (side0 - ubvar)/coef0;
      boundyub2 = (side1 - ubvar)/coef1;

      if( islhs )
      {
	 boundvaluex1 = MAX(boundxlb1, boundxlb2);
	 boundvaluex2 = MAX(boundxub1, boundxub2);
      }
      else
      {
	 boundvaluex1 = MIN(boundxlb1, boundxlb2);
	 boundvaluex2 = MIN(boundxub1, boundxub2);
      }

      /* calculate important values for variables */
      if( SCIPisPositive(scip, coef0) )
      {
         valuex1 = MIN(boundvaluex1, ubvar);
         valuex1 = MAX(valuex1, lbvar);
         valuex2 = MAX(boundvaluex2, lbvar);
         valuex2 = MIN(valuex2, ubvar);

         /* if variable is of integral type make values integral too */
         if( SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS )
         {
            if( !SCIPisFeasIntegral(scip, valuex1) )
               valuex1 = SCIPfeasFloor(scip, valuex1);
            if( !SCIPisFeasIntegral(scip, valuex2) )
               valuex2 = SCIPfeasCeil(scip, valuex2);
         }
      }
      else
      {
         valuex1 = MAX(boundvaluex1, lbvar);
         valuex1 = MIN(valuex1, ubvar);
         valuex2 = MIN(boundvaluex2, ubvar);
         valuex2 = MAX(valuex2, lbvar);

         /* if variable is of integral type make values integral too */
         if( SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS )
         {
            if( !SCIPisFeasIntegral(scip, valuex1) )
               valuex1 = SCIPfeasCeil(scip, valuex1);
            if( !SCIPisFeasIntegral(scip, valuex2) )
               valuex2 = SCIPfeasFloor(scip, valuex2);
         }
      }

      /* calculate resulting values of variable y by setting x to valuex1 */
      valuey1 = (side0 - valuex1)/coef0;
      valuey2 = (side1 - valuex1)/coef1;

      /* determine redundancy of one constraints side */
      if( valuey1 - valuey2 <= eps )
         *sideequal = TRUE;
      else if( SCIPisPositive(scip, coef0) )
      {
         if( valuey1 < valuey2 )
            *redundant1 = TRUE;
         else
            *redundant0 = TRUE;
      }
      else
      {
         if( valuey1 < valuey2 )
            *redundant0 = TRUE;
         else
            *redundant1 = TRUE;
      }

      /* calculate resulting values of variable y by setting x to valuex2 */
      valuey1 = (side0 - valuex2)/coef0;
      valuey2 = (side1 - valuex2)/coef1;

      /* determine redundancy of one constraints side by checking for the first valuex2 */
      if( SCIPisPositive(scip, coef0) )
      {
         /* if both constraints are weaker than the other on one value, we have no redundancy */
         if( (*redundant1 && valuey1 > valuey2) || (*redundant0 && valuey1 < valuey2) )
         {
            *sideequal = FALSE;
            *redundant0 = FALSE;
            *redundant1 = FALSE;
            return;
         }
         else if( *sideequal )
         {
            if( valuey1 + eps < valuey2 )
            {
               *sideequal = FALSE;
               *redundant1 = TRUE;
            }
            else if( valuey1 + eps > valuey2 )
            {
               *sideequal = FALSE;
               *redundant0 = TRUE;
            }
         }
      }
      else
      {
         /* if both constraints are weaker than the other one on one value, we have no redundancy */
         if( (*redundant1 && valuey1 < valuey2) || (*redundant0 && valuey1 > valuey2) )
         {
            *sideequal = FALSE;
            *redundant0 = FALSE;
            *redundant1 = FALSE;
            return;
         }
         else if( *sideequal )
         {
            if( valuey1 + eps < valuey2 )
            {
               *sideequal = FALSE;
               *redundant0 = TRUE;
            }
            else if( valuey1 + eps > valuey2 )
            {
               *sideequal = FALSE;
               *redundant1 = TRUE;
            }
         }
      }
      assert(*sideequal || *redundant0 || *redundant1);

      /* calculate feasibility domain values for variable y concerning these both constraints */
      if( SCIPisPositive(scip, coef0) )
      {
	 if( islhs )
	 {
	    boundvaluey1 = MAX(boundylb1, boundylb2);
	    boundvaluey2 = MAX(boundyub1, boundyub2);
	 }
	 else
	 {
	    boundvaluey1 = MIN(boundylb1, boundylb2);
	    boundvaluey2 = MIN(boundyub1, boundyub2);
	 }

         valuey1 = MIN(boundvaluey1, ubvbdvar);
         valuey1 = MAX(valuey1, lbvbdvar);
         valuey2 = MAX(boundvaluey2, lbvbdvar);
         valuey2 = MIN(valuey2, ubvbdvar);

         if( !SCIPisFeasIntegral(scip, valuey1) )
            valuey1 = SCIPfeasFloor(scip, valuey1);
         if( !SCIPisFeasIntegral(scip, valuey2) )
            valuey2 = SCIPfeasCeil(scip, valuey2);
      }
      else
      {
	 if( islhs )
	 {
	    boundvaluey1 = MIN(boundylb1, boundylb2);
	    boundvaluey2 = MIN(boundyub1, boundyub2);
	 }
	 else
	 {
	    boundvaluey1 = MAX(boundylb1, boundylb2);
	    boundvaluey2 = MAX(boundyub1, boundyub2);
	 }

         valuey1 = MAX(boundvaluey1, lbvbdvar);
         valuey1 = MIN(valuey1, ubvbdvar);
         valuey2 = MIN(boundvaluey2, ubvbdvar);
         valuey2 = MAX(valuey2, lbvbdvar);

         /* if variable is of integral type make values integral too */
         if( !SCIPisFeasIntegral(scip, valuey1) )
            valuey1 = SCIPfeasCeil(scip, valuey1);
         if( !SCIPisFeasIntegral(scip, valuey2) )
            valuey2 = SCIPfeasFloor(scip, valuey2);
      }

      /* calculate resulting values of variable x by setting y to valuey1 */
      valuex1 = side0 - valuey1*coef0;
      valuex2 = side1 - valuey1*coef1;

      /* determine redundancy of one constraints side by checking for the first valuey1 */
      if( (*redundant1 && valuex1 > valuex2) || (*redundant0 && valuex1 < valuex2) )
      {
         *sideequal = FALSE;
         *redundant0 = FALSE;
         *redundant1 = FALSE;
         return;
      }
      if( *sideequal )
      {
         if( valuex1 + eps < valuex2 )
         {
            *sideequal = FALSE;
            *redundant1 = TRUE;
         }
         else if( valuex1 + eps > valuex2 )
         {
            *sideequal = FALSE;
            *redundant0 = TRUE;
         }
      }

      /* calculate resulting values of variable x by setting y to valuey2 */
      valuex1 = side0 - valuey2*coef0;
      valuex2 = side1 - valuey2*coef1;

      /* determine redundancy of one constraints side by checking for the first valuey1 */
      if( (*redundant1 && valuex1 > valuex2) || (*redundant0 && valuex1 < valuex2) )
      {
         *sideequal = FALSE;
         *redundant0 = FALSE;
         *redundant1 = FALSE;
         return;
      }
      if( *sideequal )
      {
         if( valuex1 + eps < valuex2 )
         {
            *sideequal = FALSE;
            *redundant1 = TRUE;
         }
         else if( valuex1 + eps > valuex2 )
         {
            *sideequal = FALSE;
            *redundant0 = TRUE;
         }
      }
      assert(*redundant0 || *redundant1 || *sideequal);
   }
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint
 *
 *  we will order all constraint to have constraints with same variables next to each other to speed up presolving
 *
 *  consider two constraints like lhs1 <= x + b1*y <= rhs1 and lhs2 <= x + b2*y <= rhs2
 *  we are doing the following presolving steps:
 *
 *  if( b1 == b2 )
 *      newlhs = MAX(lhs1, lhs2)
 *      newrhs = MIN(rhs1, rhs2)
 *      updateSides
 *      delete one constraint
 *  else if( ((b1 > 0) == (b2 > 0)) && (lhs1 != -inf && lhs2 != -inf) || (rhs1 != inf && rhs2 != inf) )
 *
 *       (i.e. both constraint have either a valid lhs or a valid rhs and infinity is on the same side and the
 *             coeffcients have the same size )
 *
 *      if( y is binary )
 *          if( lhs1 != -inf )
 *              newlhs = MAX(lhs1, lhs2)
 *              newb = newlhs - MAX(lhs1 - b1, lhs2 - b2)
 *          else
 *              newrhs = MIN(lhs1, lhs2)
 *              newb = newrhs - MIN(rhs1 - b1, rhs2 - b2)
 *          updateSidesAndCoef
 *          delete one constraint
 *      else
 *          we calculate possible values for both variables and check which constraint is tighter
 *  else
 *      nothing possible
 *
 *  We also try to tighten bounds in the case of two constraints lhs1 <= x + b1*y <= rhs1 and lhs2 <= y + b2*x <= rhs2.
 *  Eliminiating one variable and inserting into the second yields the following bounds:
 *  If b2 > 0:
 *     (1 - b1 * b2) * y >= lhs2 - b2 * rhs1
 *     (1 - b1 * b2) * y <= rhs2 - b2 * lhs1
 *  If b2 < 0:
 *     (1 - b1 * b2) * y >= lhs2 - b2 * lhs1
 *     (1 - b1 * b2) * y <= rhs2 - b2 * rhs1
 *  The case of x is similar.
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
   )
{
   SCIP_CONS** sortedconss;
   int c;
   int s;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

   /* create our temporary working array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortedconss, conss, nconss) );

   /* sort all constraints, so that all constraints with same variables stand next to each other */
   SCIPsortPtr((void**)sortedconss, consVarboundComp, nconss);

   /* check all constraints for redundancy */
   for( c = nconss - 1; c > 0 && !(*cutoff); --c )
   {
      SCIP_CONS* cons0;
      SCIP_CONSDATA* consdata0;

      cons0 = sortedconss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      consdata0 = SCIPconsGetData(cons0);
      assert(consdata0 != NULL);
      assert(consdata0->var != NULL);
      assert(consdata0->vbdvar != NULL);

      /* do not check for already redundant constraints */
      assert(!SCIPisZero(scip, consdata0->vbdcoef));
      assert(!SCIPisInfinity(scip, -consdata0->lhs) || !SCIPisInfinity(scip, consdata0->rhs));

      if( !consdata0->changed )
         continue;

      consdata0->changed = FALSE;

      for( s = c - 1; s >= 0; --s )
      {
         SCIP_CONS* cons1;
         SCIP_CONSDATA* consdata1;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real coef;
         SCIP_Bool deletecons1;

         cons1 = sortedconss[s];

         if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
            continue;

         consdata1 = SCIPconsGetData(cons1);
         assert(consdata1 != NULL);
         assert(consdata1->var != NULL);
         assert(consdata1->vbdvar != NULL);

         /* do not check for already redundant constraints */
         assert(!SCIPisZero(scip, consdata1->vbdcoef));
         assert(!SCIPisInfinity(scip, -consdata1->lhs) || !SCIPisInfinity(scip, consdata1->rhs));

         lhs = consdata0->lhs;
         rhs = consdata0->rhs;
         coef = consdata0->vbdcoef;

         /* check for propagation in the case: lhs1 <= x + b1*y <= rhs1 and lhs2 <= y + b2*x <= rhs2. */
         if ( consdata0->var == consdata1->vbdvar && consdata0->vbdvar == consdata1->var &&
            !SCIPisFeasZero(scip, 1.0 - coef * consdata1->vbdcoef) )
         {
            SCIP_Bool tightened = FALSE;
            SCIP_Real bnd = SCIP_UNKNOWN;
            SCIP_Real scalar;
            SCIP_Real newbnd;

            scalar = (1.0 - coef * consdata1->vbdcoef);

            assert( ! SCIPisInfinity(scip, REALABS(scalar)) );
            assert( ! SCIPisZero(scip, consdata0->vbdcoef) );
            assert( ! SCIPisZero(scip, consdata1->vbdcoef) );

            /* lower bounds for consdata0->var */
            if ( ! SCIPisInfinity(scip, -lhs) )
            {
               if ( SCIPisPositive(scip, coef) )
               {
                  if ( ! SCIPisInfinity(scip, consdata1->rhs) )
                     bnd = (lhs - coef * consdata1->rhs)/scalar;
               }
               else
               {
                  assert( SCIPisNegative(scip, coef) );
                  if ( ! SCIPisInfinity(scip, consdata1->lhs) )
                     bnd = (lhs - coef * consdata1->lhs)/scalar;
               }

               if ( bnd != SCIP_UNKNOWN ) /*lint !e777*/
               {
                  if ( SCIPisFeasPositive(scip, scalar) )
                  {
                     newbnd = SCIPadjustedVarLb(scip, consdata0->var, bnd);
                     SCIP_CALL( SCIPtightenVarLb(scip, consdata0->var, newbnd, FALSE, cutoff, &tightened) );
                     if ( tightened )
                     {
                        SCIPdebugMsg(scip, "<%s>, <%s> -> tightened lower bound: <%s> >= %.15g\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1),
                           SCIPvarGetName(consdata0->var), SCIPvarGetLbGlobal(consdata0->var));
                        (*nchgbds)++;
                     }
                  }
                  else if ( SCIPisFeasNegative(scip, scalar) )
                  {
                     newbnd = SCIPadjustedVarUb(scip, consdata0->var, bnd);
                     SCIP_CALL( SCIPtightenVarUb(scip, consdata0->var, newbnd, FALSE, cutoff, &tightened) );
                     if ( tightened )
                     {
                        SCIPdebugMsg(scip, "<%s>, <%s> -> tightened upper bound: <%s> >= %.15g\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1),
                           SCIPvarGetName(consdata0->var), SCIPvarGetUbGlobal(consdata0->var));
                        (*nchgbds)++;
                     }
                  }
               }
            }

            /* upper bound for consdata0>var */
            if ( ! SCIPisInfinity(scip, rhs) )
            {
               bnd = SCIP_UNKNOWN;
               if ( SCIPisPositive(scip, coef) )
               {
                  if ( ! SCIPisInfinity(scip, consdata1->lhs) )
                     bnd = (rhs - coef * consdata1->lhs)/scalar;
               }
               else
               {
                  assert( SCIPisNegative(scip, coef) );
                  if ( ! SCIPisInfinity(scip, consdata1->rhs) )
                     bnd = (rhs - coef * consdata1->rhs)/scalar;
               }

               if ( bnd != SCIP_UNKNOWN ) /*lint !e777*/
               {
                  if ( SCIPisFeasPositive(scip, scalar) )
                  {
                     newbnd = SCIPadjustedVarUb(scip, consdata0->var, bnd);
                     SCIP_CALL( SCIPtightenVarUb(scip, consdata0->var, newbnd, FALSE, cutoff, &tightened) );
                     if ( tightened )
                     {
                        SCIPdebugMsg(scip, "<%s>, <%s> -> tightened upper bound: <%s> >= %.15g\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1),
                           SCIPvarGetName(consdata0->var), SCIPvarGetUbGlobal(consdata0->var));
                        (*nchgbds)++;
                     }
                  }
                  else if ( SCIPisFeasNegative(scip, scalar) )
                  {
                     newbnd = SCIPadjustedVarLb(scip, consdata0->var, bnd);
                     SCIP_CALL( SCIPtightenVarLb(scip, consdata0->var, newbnd, FALSE, cutoff, &tightened) );
                     if ( tightened )
                     {
                        SCIPdebugMsg(scip, "<%s>, <%s> -> tightened lower bound: <%s> >= %.15g\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1),
                           SCIPvarGetName(consdata0->var), SCIPvarGetLbGlobal(consdata0->var));
                        (*nchgbds)++;
                     }
                  }
               }
            }


            /* lower bounds for consdata1->var */
            if ( ! SCIPisInfinity(scip, -consdata1->lhs) )
            {
               bnd = SCIP_UNKNOWN;
               if ( SCIPisPositive(scip, consdata1->vbdcoef) )
               {
                  if ( ! SCIPisInfinity(scip, rhs) )
                     bnd = (consdata1->lhs - consdata1->vbdcoef * rhs)/scalar;
               }
               else
               {
                  assert( SCIPisNegative(scip, consdata1->vbdcoef) );
                  if ( ! SCIPisInfinity(scip, lhs) )
                     bnd = (consdata1->lhs - consdata1->vbdcoef * lhs)/scalar;
               }

               if ( bnd != SCIP_UNKNOWN ) /*lint !e777*/
               {
                  if ( SCIPisFeasPositive(scip, scalar) )
                  {
                     newbnd = SCIPadjustedVarLb(scip, consdata1->var, bnd);
                     SCIP_CALL( SCIPtightenVarLb(scip, consdata1->var, newbnd, FALSE, cutoff, &tightened) );
                     if ( tightened )
                     {
                        SCIPdebugMsg(scip, "<%s>, <%s> -> tightened lower bound: <%s> >= %.15g\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1),
                           SCIPvarGetName(consdata1->var), SCIPvarGetLbGlobal(consdata1->var));
                        (*nchgbds)++;
                     }
                  }
                  else if ( SCIPisFeasNegative(scip, scalar) )
                  {
                     newbnd = SCIPadjustedVarUb(scip, consdata1->var, bnd);
                     SCIP_CALL( SCIPtightenVarUb(scip, consdata1->var, newbnd, FALSE, cutoff, &tightened) );
                     if ( tightened )
                     {
                        SCIPdebugMsg(scip, "<%s>, <%s> -> tightened upper bound: <%s> >= %.15g\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1),
                           SCIPvarGetName(consdata1->var), SCIPvarGetUbGlobal(consdata1->var));
                        (*nchgbds)++;
                     }
                  }
               }
            }

            /* upper bound for consdata1->var */
            if ( ! SCIPisInfinity(scip, consdata1->rhs) )
            {
               bnd = SCIP_UNKNOWN;
               if ( SCIPisPositive(scip, consdata1->vbdcoef) )
               {
                  if ( ! SCIPisInfinity(scip, lhs) )
                     bnd = (consdata1->rhs - consdata1->vbdcoef * lhs)/scalar;
               }
               else
               {
                  assert( SCIPisNegative(scip, consdata1->vbdcoef) );
                  if ( ! SCIPisInfinity(scip, rhs) )
                     bnd = (consdata1->rhs - consdata1->vbdcoef * rhs)/scalar;
               }

               if ( bnd != SCIP_UNKNOWN ) /*lint !e777*/
               {
                  if ( SCIPisFeasPositive(scip, scalar) )
                  {
                     newbnd = SCIPadjustedVarUb(scip, consdata1->var, bnd);
                     SCIP_CALL( SCIPtightenVarUb(scip, consdata1->var, newbnd, FALSE, cutoff, &tightened) );
                     if ( tightened )
                     {
                        SCIPdebugMsg(scip, "<%s>, <%s> -> tightened upper bound: <%s> >= %.15g\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1),
                           SCIPvarGetName(consdata1->var), SCIPvarGetUbGlobal(consdata1->var));
                        (*nchgbds)++;
                     }
                  }
                  else if ( SCIPisFeasNegative(scip, scalar) )
                  {
                     newbnd = SCIPadjustedVarLb(scip, consdata1->var, bnd);
                     SCIP_CALL( SCIPtightenVarLb(scip, consdata1->var, newbnd, FALSE, cutoff, &tightened) );
                     if ( tightened )
                     {
                        SCIPdebugMsg(scip, "<%s>, <%s> -> tightened lower bound: <%s> >= %.15g\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1),
                           SCIPvarGetName(consdata1->var), SCIPvarGetLbGlobal(consdata1->var));
                        (*nchgbds)++;
                     }
                  }
               }
            }
         }

         /* check for equal variables */
         if( consdata0->var != consdata1->var || consdata0->vbdvar != consdata1->vbdvar )
            break;

         /* mark constraint1 for deletion if possible */
         deletecons1 = TRUE;

         /* the coefficients of both constraints are equal */
         if( SCIPisEQ(scip, coef, consdata1->vbdcoef) )
         {
            lhs = MAX(consdata1->lhs, lhs);
            rhs = MIN(consdata1->rhs, rhs);
         }
         /* now only one side and in both constraints the same side should be infinity and the vbdvar should be binary
          * then we neither do not need to have the same side nor the same coefficient
          */
         else if( SCIPvarIsBinary(consdata0->vbdvar)
            && (SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs))
            && (SCIPisInfinity(scip, -consdata1->lhs) || SCIPisInfinity(scip, consdata1->rhs))
            && (SCIPisInfinity(scip, -lhs) == SCIPisInfinity(scip, -consdata1->lhs)) )
         {
            /* lhs <= x + b*y <= +inf */
            if( !SCIPisInfinity(scip, -lhs) )
            {
               lhs = MAX(consdata1->lhs, lhs);
               coef = lhs - MAX(consdata1->lhs - consdata1->vbdcoef, consdata0->lhs - coef);
            }
            /* -inf <= x + b*y <= rhs */
            else
            {
               rhs = MIN(consdata1->rhs, rhs);
               coef = rhs - MIN(consdata1->rhs - consdata1->vbdcoef, consdata0->rhs - coef);
            }

            SCIP_CALL( SCIPmarkConsPropagate(scip, cons0) );
         }
         else if( SCIPisPositive(scip, coef) == SCIPisPositive(scip, consdata1->vbdcoef)
            && ((!SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, -consdata1->lhs))
               || (!SCIPisInfinity(scip, rhs) && !SCIPisInfinity(scip, consdata1->rhs))) )
         {
            SCIP_Bool cons0lhsred;
            SCIP_Bool cons0rhsred;
            SCIP_Bool cons1lhsred;
            SCIP_Bool cons1rhsred;
            SCIP_Bool lhsequal;
            SCIP_Bool rhsequal;

            assert(!SCIPisInfinity(scip, lhs));
            assert(!SCIPisInfinity(scip, consdata1->lhs));
            assert(!SCIPisInfinity(scip, -rhs));
            assert(!SCIPisInfinity(scip, -consdata1->rhs));

            /* check if a left hand side of one constraints is redundant */
            checkRedundancySide(scip, consdata0->var, consdata0->vbdvar, coef, consdata1->vbdcoef, lhs, consdata1->lhs, &lhsequal, &cons0lhsred, &cons1lhsred, TRUE);

            /* check if a right hand side of one constraints is redundant */
            checkRedundancySide(scip, consdata0->var, consdata0->vbdvar, coef, consdata1->vbdcoef, rhs, consdata1->rhs, &rhsequal, &cons0rhsred, &cons1rhsred, FALSE);

            /* if cons0 is redundant, update cons1 and delete cons0 */
            if( (lhsequal || cons0lhsred) && (rhsequal || cons0rhsred) )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( SCIPupdateConsFlags(scip, cons1, cons0) );

               SCIPdebugMsg(scip, "constraint: ");
               SCIPdebugPrintCons(scip, cons0, NULL);
               SCIPdebugMsg(scip, "is redundant to constraint: ");
               SCIPdebugPrintCons(scip, cons1, NULL);

               SCIP_CALL( SCIPdelCons(scip, cons0) );
               ++(*ndelconss);

               /* get next cons0 */
               break;
            }
            /* if cons1 is redundant, update cons0 and delete cons1 */
            else if( cons1lhsred && cons1rhsred )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );

               SCIPdebugMsg(scip, "constraint: ");
               SCIPdebugPrintCons(scip, cons1, NULL);
               SCIPdebugMsg(scip, "is redundant to constraint: ");
               SCIPdebugPrintCons(scip, cons0, NULL);

               SCIP_CALL( SCIPdelCons(scip, cons1) );
               ++(*ndelconss);

               /* get next cons1 */
               continue;
            }
            /* if left hand side of cons0 is redundant set it to -infinity */
            else if( (lhsequal || cons0lhsred) && !SCIPisInfinity(scip, -lhs) )
            {
               lhs = -SCIPinfinity(scip);

	       /* if right hand side of cons1 is redundant too, set it to infinity */
	       if( cons1rhsred && !SCIPisInfinity(scip, consdata1->rhs) )
	       {
		  SCIP_CALL( chgRhs(scip, cons1, SCIPinfinity(scip)) );
		  ++(*nchgsides);

		  SCIPdebugMsg(scip, "deleted rhs of constraint: ");
		  SCIPdebugPrintCons(scip, cons1, NULL);
		  SCIPdebugMsg(scip, "due to constraint: ");
		  SCIPdebugPrintCons(scip, cons0, NULL);
	       }

               /* later on we cannot not want to delete cons1 */
               deletecons1 = FALSE;
            }
            /* if right hand side of cons0 is redundant set it to infinity */
            else if( (rhsequal || cons0rhsred) && !SCIPisInfinity(scip, rhs) )
            {
               rhs = SCIPinfinity(scip);

	       /* if left hand side of cons1 is redundant too, set it to -infinity */
	       if( cons1lhsred && !SCIPisInfinity(scip, -consdata1->lhs) )
	       {
		  SCIP_CALL( chgLhs(scip, cons1, -SCIPinfinity(scip)) );
		  ++(*nchgsides);

		  SCIPdebugMsg(scip, "deleted lhs of constraint: ");
		  SCIPdebugPrintCons(scip, cons1, NULL);
		  SCIPdebugMsg(scip, "due to constraint: ");
		  SCIPdebugPrintCons(scip, cons0, NULL);
	       }

               /* later on we cannot not want to delete cons1 */
               deletecons1 = FALSE;
            }
            /* if left hand side of cons1 is redundant set it to -infinity */
            else if( cons1lhsred && !SCIPisInfinity(scip, -consdata1->lhs) )
            {
	       SCIP_CALL( chgLhs(scip, cons1, -SCIPinfinity(scip)) );
	       ++(*nchgsides);

               SCIPdebugMsg(scip, "deleted lhs of constraint: ");
               SCIPdebugPrintCons(scip, cons1, NULL);
               SCIPdebugMsg(scip, "due to constraint: ");
               SCIPdebugPrintCons(scip, cons0, NULL);

               continue;
            }
            /* if right hand side of cons1 is redundant set it to infinity */
            else if( cons1rhsred && !SCIPisInfinity(scip, consdata1->rhs) )
            {
	       SCIP_CALL( chgRhs(scip, cons1, SCIPinfinity(scip)) );
	       ++(*nchgsides);

               SCIPdebugMsg(scip, "deleted rhs of constraint: ");
               SCIPdebugPrintCons(scip, cons1, NULL);
               SCIPdebugMsg(scip, "due to constraint: ");
               SCIPdebugPrintCons(scip, cons0, NULL);

               continue;
            }
            else /* nothing was redundant */
               continue;
         }
         else
         {
            /* there is no redundancy in both constraints with same variables */
            continue;
         }

         if( SCIPisFeasLT(scip, rhs, lhs) )
         {
            SCIPdebugMsg(scip, "constraint <%s> and <%s> lead to infeasibility due to their sides\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* ensure that lhs <= rhs holds without tolerances as we only allow such rows to enter the LP */
         if( lhs > rhs )
         {
            rhs = (lhs + rhs)/2;
            lhs = rhs;
         }

         /* we decide to let constraint cons0 stay, so update data structure consdata0 */

         /* update coefficient of cons0 */

         /* special case if new coefficient becomes zero, both constraints are redundant but we may tighten the bounds */
         if( SCIPisZero(scip, coef) )
         {
            SCIP_Bool infeasible;
            SCIP_Bool tightened;

            SCIPdebugMsg(scip, "constraint: ");
            SCIPdebugPrintCons(scip, cons1, NULL);
            SCIPdebugMsg(scip, "and constraint: ");
            SCIPdebugPrintCons(scip, cons0, NULL);
            SCIPdebugMsg(scip, "are both redundant and lead to bounding of <%s> in [%g, %g]\n", SCIPvarGetName(consdata0->var), lhs, rhs);

            /* delete cons1 */
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            ++(*ndelconss);

            /* update upper bound if possible
             *
             * @note we need to force the bound change since we are deleting the constraint afterwards
             */
            SCIP_CALL( SCIPtightenVarUb(scip, consdata0->var, rhs, TRUE, &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened )
               ++(*nchgbds);

            /* update lower bound if possible
             *
             * @note we need to force the bound change since we are deleting the constraint afterwards
             */
            SCIP_CALL( SCIPtightenVarLb(scip, consdata0->var, lhs, TRUE, &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened )
               ++(*nchgbds);

            /* delete cons0 */
            SCIP_CALL( SCIPdelCons(scip, cons0) );
            ++(*ndelconss);

            /* get next cons0 */
            break;
         }

         SCIPdebugMsg(scip, "constraint: ");
         SCIPdebugPrintCons(scip, cons1, NULL);
         SCIPdebugMsg(scip, "and constraint: ");
         SCIPdebugPrintCons(scip, cons0, NULL);

         /* if sign of coefficient switches, update the rounding locks of the variable */
         if( SCIPconsIsLocked(cons0) && consdata0->vbdcoef * coef < 0.0 )
         {
            assert(SCIPconsIsTransformed(cons0));

            /* remove rounding locks for variable with old coefficient and install rounding locks for variable with new
             * coefficient
             */
            if( SCIPisPositive(scip, consdata0->vbdcoef) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, consdata0->vbdvar, cons0,
                     !SCIPisInfinity(scip, -consdata0->lhs), !SCIPisInfinity(scip, consdata0->rhs)) );
               SCIP_CALL( SCIPlockVarCons(scip, consdata0->vbdvar, cons0,
                     !SCIPisInfinity(scip, consdata0->rhs), !SCIPisInfinity(scip, -consdata0->lhs)) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, consdata0->vbdvar, cons0,
                     !SCIPisInfinity(scip, consdata0->rhs), !SCIPisInfinity(scip, -consdata0->lhs)) );
               SCIP_CALL( SCIPlockVarCons(scip, consdata0->vbdvar, cons0,
                     !SCIPisInfinity(scip, -consdata0->lhs), !SCIPisInfinity(scip, consdata0->rhs)) );
            }
         }

         /* now change the coefficient */
         if( !SCIPisEQ(scip, consdata0->vbdcoef, coef) )
         {
            ++(*nchgcoefs);

            /* mark to add new varbound information */
            consdata0->varboundsadded = FALSE;
	    consdata0->tightened = FALSE;
	    consdata0->presolved = FALSE;
	    consdata0->changed = FALSE;

	    consdata0->vbdcoef = coef;

            SCIP_CALL( SCIPmarkConsPropagate(scip, cons0) );
         }

         /* update lhs and rhs of cons0 */
         if( !SCIPisEQ(scip, consdata0->lhs, lhs) )
         {
	    SCIP_CALL( chgLhs(scip, cons0, lhs) );
	    ++(*nchgsides);
	 }
         if( !SCIPisEQ(scip, consdata0->rhs, rhs) )
         {
	    SCIP_CALL( chgRhs(scip, cons0, rhs) );
	    ++(*nchgsides);
	 }

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );

         SCIPdebugMsg(scip, "lead to new constraint: ");
         SCIPdebugPrintCons(scip, cons0, NULL);

	 /* if cons1 is still marked for deletion, delete it */
         if( deletecons1 )
         {
	    /* delete cons1 */
	    SCIP_CALL( SCIPdelCons(scip, cons1) );
	    ++(*ndelconss);
	 }

         assert(SCIPconsIsActive(cons0));
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &sortedconss);

   return SCIP_OKAY;
}

/** for all varbound constraints with two integer variables make the coefficients integral */
static
SCIP_RETCODE prettifyConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
   )
{
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

   /* if we cannot find any constraint for prettifying, stop */
   if( SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip) < 1 )
      return SCIP_OKAY;

   for( c = nconss - 1; c >= 0; --c )
   {
      assert(conss != NULL);

      if( SCIPconsIsDeleted(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* check for integer variables and one coefficient with an absolute value smaller than 1 */
      /* @note: we allow that the variable type of the bounded variable can be smaller than the variable type of the
       *        bounding variable
       */
      if( (SCIPvarGetType(consdata->var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(consdata->var) == SCIP_VARTYPE_INTEGER
	    || SCIPvarGetType(consdata->var) == SCIP_VARTYPE_IMPLINT)
	 && (SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_IMPLINT)
	 && SCIPisLT(scip, REALABS(consdata->vbdcoef), 1.0) )
      {
         SCIP_Real epsilon;
         SCIP_Longint nominator;
         SCIP_Longint denominator;
         SCIP_Longint maxmult;
         SCIP_Bool success;

         epsilon = SCIPepsilon(scip) * 0.9;  /* slightly decrease epsilon to be safe in rational conversion below */
         maxmult = (SCIP_Longint)(SCIPfeastol(scip)/epsilon + SCIPfeastol(scip));
         maxmult = MIN(maxmult, MAXSCALEDCOEF);

         success = SCIPrealToRational(consdata->vbdcoef, -epsilon, epsilon , maxmult, &nominator, &denominator);

         if( success )
         {
            /* it is possible that the dominator is a multiple of the nominator */
            if( SCIPisIntegral(scip, (SCIP_Real) denominator / (SCIP_Real) nominator) )
            {
               denominator /= nominator;
               nominator = 1;
            }

            success = success && (denominator <= maxmult);

            /* scale the constraint denominator/nominator */
            if( success && ABS(denominator) > 1 && nominator == 1)
            {
               SCIP_VAR* swapvar;

               /* print constraint before scaling */
               SCIPdebugPrintCons(scip, conss[c], NULL);

               assert(SCIPisEQ(scip, consdata->vbdcoef * denominator, 1.0));

               /* need to switch sides if coefficient is smaller then 0 */
               if( consdata->vbdcoef < 0 )
               {
                  assert(denominator < 0);

                  /* compute new sides */

                  /* only right hand side exists */
                  if( SCIPisInfinity(scip, -consdata->lhs) )
                  {
                     consdata->lhs = consdata->rhs * denominator;
                     assert(!SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->lhs));

                     consdata->rhs = SCIPinfinity(scip);
                  }
                  /* only left hand side exists */
                  else if( SCIPisInfinity(scip, consdata->rhs) )
                  {
                     consdata->rhs = consdata->lhs * denominator;
                     assert(!SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->rhs));

                     consdata->lhs = -SCIPinfinity(scip);
                  }
                  /* both sides exist */
                  else
                  {
                     SCIP_Real tmp;

                     tmp = consdata->lhs;
                     consdata->lhs = consdata->rhs * denominator;
                     consdata->rhs = tmp * denominator;
		     consdata->tightened = FALSE;

                     assert(!SCIPisInfinity(scip, consdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs));
                     assert(SCIPisGE(scip, consdata->rhs, consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs));
                  }
                  *nchgsides += 2;
               }
               /* coefficient > 0 */
               else
               {
                  assert(denominator > 0);

                  /* compute new left hand side */
                  if( !SCIPisInfinity(scip, -consdata->lhs) )
                  {
                     consdata->lhs *= denominator;
                     assert(!SCIPisInfinity(scip, consdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs));
                     ++(*nchgsides);
                  }

                  /* compute new right hand side */
                  if( !SCIPisInfinity(scip, consdata->rhs) )
                  {
                     consdata->rhs *= denominator;
                     assert(!SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->rhs));
                     ++(*nchgsides);
                  }

                  assert(SCIPisGE(scip, consdata->rhs, consdata->lhs));
               }

               /* swap both variables */
               swapvar = consdata->var;
               consdata->var = consdata->vbdvar;
               consdata->vbdvar = swapvar;

               /* swap coefficient */
               consdata->vbdcoef = (SCIP_Real)denominator;
               ++(*nchgcoefs);

               /* mark to add new varbound information */
               consdata->varboundsadded = FALSE;
	       consdata->tightened = FALSE;

               /* print constraint after scaling */
               SCIPdebugMsg(scip, "transformed into:");
               SCIPdebugPrintCons(scip, conss[c], NULL);
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** replaces fixed and aggregated variables in variable bound constraint by active problem variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_Bool*            cutoff,             /**< pointer to store whether an infeasibility was detected */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  naddconss           /**< pointer to count number of added constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real varscalar;
   SCIP_Real varconstant;
   SCIP_VAR* vbdvar;
   SCIP_Real vbdvarscalar;
   SCIP_Real vbdvarconstant;
   SCIP_Bool varschanged;
   SCIP_Bool redundant;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);

   *cutoff = FALSE;
   redundant = FALSE;

   /* the variable bound constraint is: lhs <= x + c*y <= rhs */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get active problem variables of x and y */
   var = consdata->var;
   varscalar = 1.0;
   varconstant = 0.0;
   SCIP_CALL( SCIPgetProbvarSum(scip, &var, &varscalar, &varconstant) );
   vbdvar = consdata->vbdvar;
   vbdvarscalar = 1.0;
   vbdvarconstant = 0.0;
   SCIP_CALL( SCIPgetProbvarSum(scip, &vbdvar, &vbdvarscalar, &vbdvarconstant) );
   varschanged = (var != consdata->var || vbdvar != consdata->vbdvar);

   /* if the variables are equal, the variable bound constraint reduces to standard bounds on the single variable */
   if( var == vbdvar && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      SCIPdebugMsg(scip, "variable bound constraint <%s> has equal variable and vbd variable <%s>\n",
         SCIPconsGetName(cons), SCIPvarGetName(var));

      /*      lhs <= a1*z + b1 + c(a2*z + b2) <= rhs
       * <=>  lhs <= (a1 + c*a2)z + (b1 + c*b2) <= rhs
       */
      scalar = varscalar + consdata->vbdcoef * vbdvarscalar;
      constant = varconstant + consdata->vbdcoef * vbdvarconstant;
      if( SCIPisZero(scip, scalar) )
      {
         /* no variable is left: the constraint is redundant or infeasible */
         if( SCIPisFeasLT(scip, constant, consdata->lhs) || SCIPisFeasGT(scip, constant, consdata->rhs) )
            *cutoff = TRUE;
      }
      else if( scalar > 0.0 )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;

            SCIP_CALL( SCIPtightenVarLb(scip, var, (consdata->lhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMsg(scip, " -> tightened lower bound: <%s> >= %.15g\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
               (*nchgbds)++;
            }
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;

            SCIP_CALL( SCIPtightenVarUb(scip, var, (consdata->rhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMsg(scip, " -> tightened upper bound: <%s> <= %.15g\n", SCIPvarGetName(var), SCIPvarGetUbGlobal(var));
               (*nchgbds)++;
            }
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;

            SCIP_CALL( SCIPtightenVarUb(scip, var, (consdata->lhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMsg(scip, " -> tightened upper bound: <%s> <= %.15g\n", SCIPvarGetName(var), SCIPvarGetUbGlobal(var));
               (*nchgbds)++;
            }
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;

            SCIP_CALL( SCIPtightenVarLb(scip, var, (consdata->rhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMsg(scip, " -> tightened lower bound: <%s> >= %.15g\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
               (*nchgbds)++;
            }
         }
      }
      redundant = TRUE;
   }
   else
   {
      /* if the variables should be replaced, drop the events and catch the events on the new variables afterwards */
      if( varschanged )
      {
         SCIP_CALL( dropEvents(scip, cons, eventhdlr) );
      }

      /* apply aggregation on x */
      if( SCIPisZero(scip, varscalar) )
      {
         SCIPdebugMsg(scip, "variable bound constraint <%s>: variable <%s> is fixed to %.15g\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->var), varconstant);

         /* cannot change bounds on multi-aggregated variables */
         if( SCIPvarGetStatus(vbdvar) != SCIP_VARSTATUS_MULTAGGR )
         {
            /* x is fixed to varconstant: update bounds of y and delete the variable bound constraint */
            if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
            {
               if( consdata->vbdcoef > 0.0 )
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarLb(scip, consdata->vbdvar, (consdata->lhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMsg(scip, " -> tightened lower bound: <%s> >= %.15g\n", SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
               else
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vbdvar, (consdata->lhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMsg(scip, " -> tightened upper bound: <%s> <= %.15g\n", SCIPvarGetName(consdata->vbdvar), SCIPvarGetUbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
            {
               if( consdata->vbdcoef > 0.0 )
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vbdvar, (consdata->rhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMsg(scip, " -> tightened upper bound: <%s> <= %.15g\n", SCIPvarGetName(consdata->vbdvar), SCIPvarGetUbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
               else
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarLb(scip, consdata->vbdvar, (consdata->rhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMsg(scip, " -> tightened lower bound: <%s> >= %.15g\n", SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
            }
            redundant = TRUE;
         }
      }
      else if( var != consdata->var )
      {
         /* replace aggregated variable x in the constraint by its aggregation */
         if( varscalar > 0.0 )
         {
            /* lhs := (lhs - varconstant) / varscalar
             * rhs := (rhs - varconstant) / varscalar
             * c   := c / varscalar
             */
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs = (consdata->lhs - varconstant)/varscalar;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               consdata->rhs = (consdata->rhs - varconstant)/varscalar;
            consdata->vbdcoef /= varscalar;

            /* try to avoid numerical troubles */
            if( SCIPisIntegral(scip, consdata->vbdcoef) )
               consdata->vbdcoef = SCIPround(scip, consdata->vbdcoef);

            consdata->tightened = FALSE;
         }
         else
         {
            SCIP_Real lhs;

            assert(varscalar != 0.0);

            /* lhs := (rhs - varconstant) / varscalar
             * rhs := (lhs - varconstant) / varscalar
             * c   := c / varscalar
             */
            lhs = consdata->lhs;
            consdata->lhs = -consdata->rhs;
            consdata->rhs = -lhs;
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs = (consdata->lhs + varconstant)/(-varscalar);
            if( !SCIPisInfinity(scip, consdata->rhs) )
               consdata->rhs = (consdata->rhs + varconstant)/(-varscalar);
            consdata->vbdcoef /= varscalar;

            /* try to avoid numerical troubles */
            if( SCIPisIntegral(scip, consdata->vbdcoef) )
               consdata->vbdcoef = SCIPround(scip, consdata->vbdcoef);

            consdata->tightened = FALSE;
         }
         /* release old variable */
         SCIP_CALL( SCIPreleaseVar(scip, &(consdata->var)) );
         consdata->var = var;
         /* capture new variable */
         SCIP_CALL( SCIPcaptureVar(scip, consdata->var) );
      }

      /* apply aggregation on y */
      if( SCIPisZero(scip, vbdvarscalar) )
      {
         SCIPdebugMsg(scip, "variable bound constraint <%s>: vbd variable <%s> is fixed to %.15g\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vbdvar), vbdvarconstant);

         /* cannot change bounds on multi-aggregated variables */
         if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
         {
            /* y is fixed to vbdvarconstant: update bounds of x and delete the variable bound constraint */
            if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
            {
               SCIP_Bool tightened;

               SCIP_CALL( SCIPtightenVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * vbdvarconstant,
                     TRUE, cutoff, &tightened) );
               if( tightened )
               {
                  SCIPdebugMsg(scip, " -> tightened lower bound: <%s> >= %.15g\n", SCIPvarGetName(consdata->var), SCIPvarGetLbGlobal(consdata->var));
                  (*nchgbds)++;
               }
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
            {
               SCIP_Bool tightened;

               SCIP_CALL( SCIPtightenVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * vbdvarconstant,
                     TRUE, cutoff, &tightened) );
               if( tightened )
               {
                  SCIPdebugMsg(scip, " -> tightened upper bound: <%s> <= %.15g\n", SCIPvarGetName(consdata->var), SCIPvarGetUbGlobal(consdata->var));
                  (*nchgbds)++;
               }
            }
            redundant = TRUE;
         }
      }
      else if( vbdvar != consdata->vbdvar )
      {
         /* replace aggregated variable y in the constraint by its aggregation:
          * lhs := lhs - c * vbdvarconstant
          * rhs := rhs - c * vbdvarconstant
          * c   := c * vbdvarscalar
          */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhs -= consdata->vbdcoef * vbdvarconstant;
         if( !SCIPisInfinity(scip, consdata->rhs) )
            consdata->rhs -= consdata->vbdcoef * vbdvarconstant;
         consdata->vbdcoef *= vbdvarscalar;

         consdata->tightened = FALSE;

         /* release old variable */
         SCIP_CALL( SCIPreleaseVar(scip, &(consdata->vbdvar)) );
         consdata->vbdvar = vbdvar;
         /* capture new variable */
         SCIP_CALL( SCIPcaptureVar(scip, consdata->vbdvar) );
      }

      /* catch the events again on the new variables */
      if( varschanged )
      {
         SCIP_CALL( catchEvents(scip, cons, eventhdlr) );
      }
   }

   /* mark constraint changed, if a variable was exchanged */
   if( varschanged )
   {
      consdata->changed = TRUE;
   }

   /* active multi aggregations are now resolved by creating a new linear constraint */
   if( !(*cutoff) && !redundant && (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(vbdvar) == SCIP_VARSTATUS_MULTAGGR) )
   {
      SCIP_CONS* newcons;
      SCIP_Real lhs;
      SCIP_Real rhs;

      lhs = consdata->lhs;
      rhs = consdata->rhs;

      /* create upgraded linear constraint */
      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, SCIPconsGetName(cons), 0, NULL, NULL, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

      /* if var was fixed, then the case that vbdvar was multi-aggregated, was not yet resolved */
      if( var != consdata->var )
      {
	 assert(SCIPvarGetStatus(vbdvar) == SCIP_VARSTATUS_MULTAGGR);
	 assert(SCIPisZero(scip, varscalar)); /* this means that var was fixed */

	 /* add offset that results of the fixed variable */
	 if( SCIPisZero(scip, varconstant) != 0 )
	 {
	    if( !SCIPisInfinity(scip, rhs) )
	    {
	       SCIP_CALL( SCIPchgRhsLinear(scip, newcons, rhs - varconstant) );
	    }
	    if( !SCIPisInfinity(scip, -lhs) )
	    {
	       SCIP_CALL( SCIPchgLhsLinear(scip, newcons, lhs - varconstant) );
	    }
	 }
      }
      else
      {
	 assert(var == consdata->var);

	 SCIP_CALL( SCIPaddCoefLinear(scip, newcons, consdata->var, 1.0) );
      }

      /* if vbdvar was fixed, then the case that var was multi-aggregated, was not yet resolved */
      if( vbdvar != consdata->vbdvar )
      {
	 assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
	 assert(SCIPisZero(scip, vbdvarscalar)); /* this means that var was fixed */

	 /* add offset that results of the fixed variable */
	 if( SCIPisZero(scip, vbdvarconstant) != 0 )
	 {
	    if( !SCIPisInfinity(scip, rhs) )
	    {
	       SCIP_CALL( SCIPchgRhsLinear(scip, newcons, rhs - vbdvarconstant) );
	    }
	    if( !SCIPisInfinity(scip, -lhs) )
	    {
	       SCIP_CALL( SCIPchgLhsLinear(scip, newcons, lhs - vbdvarconstant) );
	    }
	 }
      }
      else
      {
	 assert(vbdvar == consdata->vbdvar);

	 SCIP_CALL( SCIPaddCoefLinear(scip, newcons, consdata->vbdvar, consdata->vbdcoef) );
      }

      SCIP_CALL( SCIPaddCons(scip, newcons) );

      SCIPdebugMsg(scip, "resolved multi aggregation in varbound constraint <%s> by creating a new linear constraint\n", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, newcons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

      redundant = TRUE;
      ++(*naddconss);
   }

   /* delete a redundant constraint */
   if( !(*cutoff) && redundant )
   {
      SCIPdebugMsg(scip, " -> variable bound constraint <%s> is redundant\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss)++;
   }

   return SCIP_OKAY;
}

/** tightens variable bound coefficient by inspecting the global bounds of the involved variables; note: this is also
 *  performed by the linear constraint handler - only necessary if the user directly creates variable bound constraints
 */
static
SCIP_RETCODE tightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   int*                  nchgsides,          /**< pointer to count the number of left and right hand sides */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count number of bound changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real oldcoef;
   int oldnchgcoefs;
   int oldnchgsides;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* tightening already done */
   if( consdata->tightened )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "tightening coefficients on variable bound constraint <%s>\n", SCIPconsGetName(cons));

   consdata->tightened = TRUE;

   /* if values and variable are integral the sides should it be too */
   if( SCIPvarGetType(consdata->var) <= SCIP_VARTYPE_IMPLINT
      && SCIPvarGetType(consdata->vbdvar) <= SCIP_VARTYPE_IMPLINT
      && SCIPisIntegral(scip, consdata->vbdcoef) )
   {
      if( !SCIPisIntegral(scip, consdata->lhs) )
      {
         consdata->lhs = SCIPfeasCeil(scip, consdata->lhs);
         ++(*nchgsides);
         consdata->changed = TRUE;
      }
      if( !SCIPisIntegral(scip, consdata->rhs) )
      {
         consdata->rhs = SCIPfeasFloor(scip, consdata->rhs);
         ++(*nchgsides);
         consdata->changed = TRUE;
      }
   }

   /* coefficient tightening only works for binary bound variable */
   if( !SCIPvarIsBinary(consdata->vbdvar) )
      return SCIP_OKAY;

   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;
   oldcoef = consdata->vbdcoef;

   /* coefficients tightening when all variables are integer */
   /* we consider the following varbound constraint: lhs <= x + b*y <= rhs (sides are possibly infinity)
    * y should always be binary and x of integral type and b not integral, we also need at least one side with infinity
    * or not integral value.
    *
    * 1. if( (lhs is integral and not -infinity) and ((rhs is infinity) or (b - floor(b) <= rhs - floor(rhs))) ):
    *
    *        lhs <= x + b*y <= rhs =>   lhs <= x + floor(b)*y <= floor(rhs)
    *
    * 2. if( (rhs is integral and not infinity) and ((lhs is -infinity) or (b - floor(b) >= lhs - floor(lhs))) ):
    *
    *        lhs <= x + b*y <= rhs   =>  ceil(lhs) <= x + ceil(b)*y <= rhs
    *
    * 3. if( ((lhs is -infinity) or (b - floor(b) >= lhs - floor(lhs)))
    *       and ((rhs is infinity) or (b - floor(b) > rhs - floor(rhs))) ):
    *
    *        lhs <= x + b*y <= rhs  =>   ceil(lhs) <= x + ceil(b)*y <= floor(rhs)
    *
    * 4. if( ((lhs is -infinity) or (b - floor(b) < lhs - floor(lhs)))
    *       and ((rhs is infinity) or (b - floor(b) <= rhs - floor(rhs))) ):
    *
    *        lhs <= x + b*y <= rhs  =>   ceil(lhs) <= x + floor(b)*y <= floor(rhs)
    *
    * 5. if( (lhs is not integral) or (rhs is not integral) )
    *
    *       if (lhs is not -infinity)
    *          if (b - floor(b) < lhs - floor(lhs)):
    *
    *             lhs <= x + b*y  =>   ceil(lhs) <= x + b*y
    *
    *          else if (b - floor(b) > lhs - floor(lhs)):
    *
    *             lhs <= x + b*y  =>   floor(lhs) + b - floor(b) <= x + b*y
    *
    *       if (rhs is not infinity)
    *          if (b - floor(b) < rhs - floor(rhs)):
    *
    *             x + b*y <= rhs  =>   x + b*y <= floor(rhs) + b - floor(b)
    *
    *          else if (b - floor(b) > rhs - floor(rhs)):
    *
    *             x + b*y <= rhs  =>   x + b*y <= floor(rhs)
    */
   if( (SCIPvarGetType(consdata->var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(consdata->var) == SCIP_VARTYPE_IMPLINT || SCIPvarGetType(consdata->var) == SCIP_VARTYPE_BINARY)
      && !SCIPisIntegral(scip, consdata->vbdcoef)
      && (!SCIPisIntegral(scip, consdata->lhs) || SCIPisInfinity(scip, -consdata->lhs)
	 || !SCIPisIntegral(scip, consdata->rhs) || SCIPisInfinity(scip, consdata->rhs)) )
   {
      /* infinity should be an integral value */
      assert(!SCIPisInfinity(scip, -consdata->lhs) || SCIPisIntegral(scip, consdata->lhs));
      assert(!SCIPisInfinity(scip, consdata->rhs) || SCIPisIntegral(scip, consdata->rhs));

      /* should not be a redundant constraint */
      assert(!SCIPisInfinity(scip, consdata->rhs) || !SCIPisInfinity(scip, -consdata->lhs));

      /* case 1 */
      if( SCIPisIntegral(scip, consdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs) &&
         (SCIPisInfinity(scip, consdata->rhs) || SCIPisFeasLE(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfeasFloor(scip, consdata->rhs))) )
      {
         consdata->vbdcoef = SCIPfeasFloor(scip, consdata->vbdcoef);
         ++(*nchgcoefs);

         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            consdata->rhs = SCIPfeasFloor(scip, consdata->rhs);
            ++(*nchgsides);
         }
      }
      /* case 2 */
      else if( SCIPisIntegral(scip, consdata->rhs) && !SCIPisInfinity(scip, consdata->rhs) &&
         (SCIPisInfinity(scip, -consdata->lhs) || SCIPisFeasGE(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfeasFloor(scip, consdata->lhs))) )

      {
         consdata->vbdcoef = SCIPfeasCeil(scip, consdata->vbdcoef);
         ++(*nchgcoefs);

         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            if( !SCIPisIntegral(scip, consdata->lhs) )
               ++(*nchgsides);

            consdata->lhs = SCIPfeasCeil(scip, consdata->lhs);
         }
      }
      /* case 3 */
      else if( (SCIPisInfinity(scip, -consdata->lhs) || SCIPisFeasGE(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfeasFloor(scip, consdata->lhs))) && (SCIPisInfinity(scip, consdata->rhs) || SCIPisFeasGT(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfeasFloor(scip, consdata->rhs))) )
      {
         consdata->vbdcoef = SCIPfeasCeil(scip, consdata->vbdcoef);
         ++(*nchgcoefs);

         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            if( !SCIPisIntegral(scip, consdata->lhs) )
               ++(*nchgsides);

            consdata->lhs = SCIPfeasCeil(scip, consdata->lhs);
         }
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            if( !SCIPisIntegral(scip, consdata->rhs) )
               ++(*nchgsides);

            consdata->rhs = SCIPfeasFloor(scip, consdata->rhs);
         }
      }
      /* case 4 */
      else if( (SCIPisInfinity(scip, -consdata->lhs) || SCIPisFeasLT(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfeasFloor(scip, consdata->lhs))) && (SCIPisInfinity(scip, consdata->rhs) || SCIPisFeasLE(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfeasFloor(scip, consdata->rhs))) )
      {
         consdata->vbdcoef = SCIPfeasFloor(scip, consdata->vbdcoef);
         ++(*nchgcoefs);

         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            if( !SCIPisIntegral(scip, consdata->lhs) )
               ++(*nchgsides);

            consdata->lhs = SCIPfeasCeil(scip, consdata->lhs);
         }
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            if( !SCIPisIntegral(scip, consdata->rhs) )
               ++(*nchgsides);

            consdata->rhs = SCIPfeasFloor(scip, consdata->rhs);
         }
      }
      /* case 5 */
      if( !SCIPisFeasIntegral(scip, consdata->lhs) || !SCIPisFeasIntegral(scip, consdata->rhs) )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            if( SCIPisFeasLT(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfeasFloor(scip, consdata->lhs)) )
            {
               consdata->lhs = SCIPfeasCeil(scip, consdata->lhs);
               ++(*nchgsides);
            }
            else if( SCIPisFeasGT(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfeasFloor(scip, consdata->lhs)) )
            {
               consdata->lhs = SCIPfeasFloor(scip, consdata->lhs) + (consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef));
               ++(*nchgsides);
            }
         }
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            if( SCIPisFeasLT(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfeasFloor(scip, consdata->rhs)) )
            {
               consdata->rhs = SCIPfeasFloor(scip, consdata->rhs) + (consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef));
               ++(*nchgsides);
            }
            else if( SCIPisFeasGT(scip, consdata->vbdcoef - SCIPfeasFloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfeasFloor(scip, consdata->rhs)) )
            {
               consdata->rhs = SCIPfeasFloor(scip, consdata->rhs);
               ++(*nchgsides);
            }
         }
      }
   }

   /* check if due to tightening the constraint got redundant */
   if( SCIPisZero(scip, consdata->vbdcoef) )
   {
      /* we have to make sure that the induced bound(s) is (are) actually applied;
       * if the relative change is too small, this may have been skipped in propagation
       */
      if( SCIPisLT(scip, SCIPvarGetLbGlobal(consdata->var), consdata->lhs) )
      {
         SCIP_Bool tightened;

         SCIP_CALL( SCIPtightenVarLbGlobal(scip, consdata->var, consdata->lhs, TRUE, cutoff, &tightened) );

         if( tightened )
         {
            SCIPdebugMsg(scip, " -> tighten domain of <%s> to [%.15g,%.15g]\n", SCIPvarGetName(consdata->var),
               SCIPvarGetLbGlobal(consdata->var), SCIPvarGetUbGlobal(consdata->var));
            (*nchgbds)++;
         }
      }
      if( SCIPisGT(scip, SCIPvarGetUbGlobal(consdata->var), consdata->rhs) )
      {
         SCIP_Bool tightened;

         SCIP_CALL( SCIPtightenVarUbGlobal(scip, consdata->var, consdata->rhs, TRUE, cutoff, &tightened) );

         if( tightened )
         {
            SCIPdebugMsg(scip, " -> tighten domain of <%s> to [%.15g,%.15g]\n", SCIPvarGetName(consdata->var),
               SCIPvarGetLbGlobal(consdata->var), SCIPvarGetUbGlobal(consdata->var));
            (*nchgbds)++;
         }
      }

      SCIPdebugMsg(scip, " -> variable bound constraint <%s> is redundant\n", SCIPconsGetName(cons));

      /* in order to correctly update the rounding locks, we need the coefficient to have the same sign as before the
       * coefficient tightening
       */
      consdata->vbdcoef = oldcoef;

      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);

      return SCIP_OKAY;
   }

   /* get bounds of variable x */
   xlb = SCIPvarGetLbGlobal(consdata->var);
   xub = SCIPvarGetUbGlobal(consdata->var);

   /* it can happen that var is not of varstatus SCIP_VARSTATUS_FIXED but the bounds are equal, in this case we need to
    * stop
    */
   if( SCIPisEQ(scip, xlb, xub) )
      return SCIP_OKAY;

   /* modification of coefficient checking for slack in constraints */
   if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* lhs <= x + c*y <= rhs  =>  lhs - c*y <= x <= rhs - c*y */
      if( consdata->vbdcoef > 0.0 && SCIPisFeasGT(scip, xlb, consdata->lhs - consdata->vbdcoef) && SCIPisFeasLT(scip, xub, consdata->rhs) )
      {
         SCIP_Real newcoef;
         SCIP_Real newrhs;
         SCIP_Real oldrhs;

         oldrhs = consdata->rhs;

         /* constraint has positive slack for both non-restricting cases y = 0, or y = 1, respectively
          * -> modify coefficients such that constraint is tight in at least one of the non-restricting cases
          * -> c' = MAX(c - rhs + xub, lhs - xlb), rhs' = rhs - c + c'
          */
         newcoef = MAX(consdata->vbdcoef - consdata->rhs + xub, consdata->lhs - xlb);
         newrhs = consdata->rhs - consdata->vbdcoef + newcoef;

         SCIPdebugMsg(scip, "tighten varbound %.15g <= <%s>[%.15g,%.15g] %+.15g<%s> <= %.15g to %.15g <= <%s> %+.15g<%s> <= %.15g\n",
            consdata->lhs, SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
            consdata->lhs, SCIPvarGetName(consdata->var), newcoef, SCIPvarGetName(consdata->vbdvar), newrhs);

         /* we cannot allow that the coefficient changes the sign because of the rounding locks */
         assert(consdata->vbdcoef * newcoef > 0);

         consdata->vbdcoef = newcoef;
         consdata->rhs = newrhs;
         (*nchgcoefs)++;
         (*nchgsides)++;

         /* some of the cases 1. to 5. might be applicable after changing the rhs to an integral value; one example is
          * the varbound constraint 0.225 <= x - 1.225 y <= 0.775 for which none of the above cases apply but after
          * tightening the lhs to 0.0 it is possible to reduce the rhs by applying the 1. reduction
          */
         if( !SCIPisFeasIntegral(scip, oldrhs) && SCIPisFeasIntegral(scip, newrhs) )
         {
            consdata->tightened = FALSE;
            SCIP_CALL( tightenCoefs(scip, cons, nchgcoefs, nchgsides, ndelconss, cutoff, nchgbds) );
            assert(consdata->tightened);
         }
         else
            consdata->tightened = (SCIPisIntegral(scip, consdata->vbdcoef) && SCIPisIntegral(scip, consdata->rhs));
      }
      else if( consdata->vbdcoef < 0.0 && SCIPisFeasGT(scip, xlb, consdata->lhs) && SCIPisFeasLT(scip, xub, consdata->rhs - consdata->vbdcoef) )
      {
         SCIP_Real newcoef;
         SCIP_Real newlhs;
         SCIP_Real oldlhs;

         oldlhs = consdata->lhs;

         /* constraint has positive slack for both non-restricting cases y = 0, or y = 1, respectively
          * -> modify coefficients such that constraint is tight in at least one of the non-restricting cases
          * -> c' = MIN(c - lhs + xlb, rhs - xub), lhs' = lhs - c + c'
          */
         newcoef = MIN(consdata->vbdcoef - consdata->lhs + xlb, consdata->rhs - xub);
         newlhs = consdata->lhs - consdata->vbdcoef + newcoef;

         SCIPdebugMsg(scip, "tighten varbound %.15g <= <%s>[%.15g,%.15g] %+.15g<%s> <= %.15g to %.15g <= <%s> %+.15g<%s> <= %.15g\n",
            consdata->lhs, SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
            newlhs, SCIPvarGetName(consdata->var), newcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs);

         /* we cannot allow that the coefficient changes the sign because of the rounding locks */
         assert(consdata->vbdcoef * newcoef > 0);

         consdata->vbdcoef = newcoef;
         consdata->lhs = newlhs;
         (*nchgcoefs)++;
         (*nchgsides)++;

         /* some of the cases 1. to 5. might be applicable after changing the rhs to an integral value; one example is
          * the varbound constraint 0.225 <= x - 1.225 y <= 0.775 for which none of the above cases apply but after
          * tightening the lhs to 0.0 it is possible to reduce the rhs by applying the 1. reduction
          */
         if( !SCIPisFeasIntegral(scip, oldlhs) && SCIPisFeasIntegral(scip, newlhs) )
         {
            consdata->tightened = FALSE;
            SCIP_CALL( tightenCoefs(scip, cons, nchgcoefs, nchgsides, ndelconss, cutoff, nchgbds) );
            assert(consdata->tightened);
         }
         else
            consdata->tightened = (SCIPisIntegral(scip, consdata->vbdcoef) && SCIPisIntegral(scip, consdata->lhs));
      }
   }
   else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs) )
   {
      /* lhs <= x + c*y  =>  x >= lhs - c*y */
      if( consdata->vbdcoef > 0.0 && SCIPisFeasGT(scip, xlb, consdata->lhs - consdata->vbdcoef) )
      {
         /* constraint has positive slack for the non-restricting case y = 1
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 1 and equivalent in the restricting case y = 0
          * -> c' = lhs - xlb
          */
         SCIPdebugMsg(scip, "tighten binary VLB <%s>[%.15g,%.15g] %+.15g<%s> >= %.15g to <%s> %+.15g<%s> >= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs,
            SCIPvarGetName(consdata->var), consdata->lhs - xlb, SCIPvarGetName(consdata->vbdvar), consdata->lhs);

         /* we cannot allow that the coefficient changes the sign because of the rounding locks */
         assert(consdata->vbdcoef * (consdata->lhs - xlb) > 0);

         consdata->vbdcoef = consdata->lhs - xlb;
         (*nchgcoefs)++;
      }
      else if( consdata->vbdcoef < 0.0 && SCIPisFeasGT(scip, xlb, consdata->lhs) )
      {
         /* constraint has positive slack for the non-restricting case y = 0
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 0 and equivalent in the restricting case y = 1
          * -> c' = c - lhs + xlb, lhs' = xlb
          */
         SCIPdebugMsg(scip, "tighten binary VLB <%s>[%.15g,%.15g] %+.15g<%s> >= %.15g to <%s> %+.15g<%s> >= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs,
            SCIPvarGetName(consdata->var), consdata->vbdcoef - consdata->lhs + xlb, SCIPvarGetName(consdata->vbdvar), xlb);

         /* we cannot allow that the coefficient changes the sign because of the rounding locks */
         assert(consdata->vbdcoef * (consdata->vbdcoef - consdata->lhs + xlb) > 0);

         consdata->vbdcoef = consdata->vbdcoef - consdata->lhs + xlb;
         consdata->lhs = xlb;
         (*nchgcoefs)++;
         (*nchgsides)++;
      }
   }
   else if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* x + c*y <= rhs  =>  x <= rhs - c*y */
      if( consdata->vbdcoef < 0.0 && SCIPisFeasLT(scip, xub, consdata->rhs - consdata->vbdcoef) )
      {
         /* constraint has positive slack for the non-restricting case y = 1
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 1 and equivalent in the restricting case y = 0
          * -> c' = rhs - xub
          */
         SCIPdebugMsg(scip, "tighten binary VUB <%s>[%.15g,%.15g] %+.15g<%s> <= %.15g to <%s> %+.15g<%s> <= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
            SCIPvarGetName(consdata->var), consdata->rhs - xub, SCIPvarGetName(consdata->vbdvar), consdata->rhs);

         /* we cannot allow that the coefficient changes the sign because of the rounding locks */
         assert(consdata->vbdcoef * (consdata->rhs - xub) > 0);

         consdata->vbdcoef = consdata->rhs - xub;
         (*nchgcoefs)++;
      }
      else if( consdata->vbdcoef > 0.0 && SCIPisFeasLT(scip, xub, consdata->rhs) )
      {
         /* constraint has positive slack for the non-restricting case y = 0
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 0 and equivalent in the restricting case y = 1
          * -> c' = c - rhs + xub, rhs' = xub
          */
         SCIPdebugMsg(scip, "tighten binary VUB <%s>[%.15g,%.15g] %+.15g<%s> <= %.15g to <%s> %+.15g<%s> <= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
            SCIPvarGetName(consdata->var), consdata->vbdcoef - consdata->rhs + xub, SCIPvarGetName(consdata->vbdvar), xub);

         /* we cannot allow that the coefficient changes the sign because of the rounding locks */
         assert(consdata->vbdcoef * (consdata->vbdcoef - consdata->rhs + xub) > 0);

         consdata->vbdcoef = consdata->vbdcoef - consdata->rhs + xub;
         consdata->rhs = xub;
         (*nchgcoefs)++;
         (*nchgsides)++;
      }
   }

   /* if something a coefficient or side of the varbound constraint was changed, ensure that the variable lower or
    * upper bounds of the variables are informed
    */
   if( *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
   {
      consdata->varboundsadded = FALSE;
      consdata->changed = TRUE;
   }

   return SCIP_OKAY;
}

/** check if we can upgrade to a set-packing constraint */
static
SCIP_RETCODE upgradeConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  naggrvars,          /**< pointer to count the number of aggregated variables */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   int*                  nchgsides,          /**< pointer to count the number of left and right hand sides */
   int*                  ndelconss,          /**< pointer to count the number of deleted constraints */
   int*                  naddconss           /**< pointer to count the number of added constraints */
   )
{
   SCIP_VAR* vars[2];
   SCIP_CONS* newcons;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(conss != NULL || nconss == 0);
   assert(cutoff != NULL);
   assert(naggrvars != NULL);
   assert(nchgbds != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);

   /* if we cannot find any constraint for upgrading, stop */
   if( SCIPgetNBinVars(scip) + SCIPgetNImplVars(scip) <= 1 )
      return SCIP_OKAY;

   if( nconss == 0 )
      return SCIP_OKAY;

   assert(conss != NULL);

   for( c = nconss - 1; c >= 0; --c )
   {
      cons = conss[c];
      assert(cons != NULL);

      if( !SCIPconsIsActive(cons) )
	 continue;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      assert(SCIPisLE(scip, consdata->lhs, consdata->rhs));

      if( !consdata->presolved )
      {
         /* incorporate fixings and aggregations in constraint */
         SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr, cutoff, nchgbds, ndelconss, naddconss) );

         if( *cutoff )
            return SCIP_OKAY;
         if( !SCIPconsIsActive(cons) )
            continue;
      }

      if( SCIPconsIsMarkedPropagate(cons) )
      {
         /* propagate constraint */
         SCIP_CALL( propagateCons(scip, cons, conshdlrdata->usebdwidening, cutoff, nchgbds, nchgsides, ndelconss) );

         if( *cutoff )
            return SCIP_OKAY;
         if( !SCIPconsIsActive(cons) )
            continue;
      }

      if( !consdata->tightened )
      {
         /* tighten variable bound coefficient */
         SCIP_CALL( tightenCoefs(scip, cons, nchgcoefs, nchgsides, ndelconss, cutoff, nchgbds) );

         if( *cutoff )
            return SCIP_OKAY;
         if( !SCIPconsIsActive(cons) )
            continue;

         assert(SCIPisLE(scip, consdata->lhs, consdata->rhs));
      }

      /* check if both variables are of binary type */
      if( SCIPvarIsBinary(consdata->vbdvar) && SCIPvarIsBinary(consdata->var) )
      {
	 /* coefficient and sides should be tightened and we assume that the constraint is not redundant */
	 assert(SCIPisEQ(scip, REALABS(consdata->vbdcoef), 1.0));
	 assert(SCIPisZero(scip, consdata->rhs) || SCIPisEQ(scip, consdata->rhs, 1.0) || SCIPisInfinity(scip, consdata->rhs));
	 assert(SCIPisZero(scip, consdata->lhs) || SCIPisEQ(scip, consdata->lhs, 1.0) || SCIPisInfinity(scip, -consdata->lhs));
	 assert(!SCIPisInfinity(scip, consdata->rhs) || !SCIPisInfinity(scip, -consdata->lhs));

	 /* the case x + y <= 1 or x + y >= 1 */
	 if( consdata->vbdcoef > 0.0 )
	 {
	    if( SCIPisEQ(scip, consdata->rhs, 1.0) )
	    {
	       /* check for aggregations like x + y == 1 */
	       if( SCIPisEQ(scip, consdata->lhs, 1.0) )
	       {
		  SCIP_Bool infeasible;
		  SCIP_Bool redundant;
		  SCIP_Bool aggregated;

		  SCIPdebugMsg(scip, "varbound constraint <%s>: aggregate <%s> + <%s> == 1\n",
		     SCIPconsGetName(cons), SCIPvarGetName(consdata->var), SCIPvarGetName(consdata->vbdvar));

		  /* aggregate both variables */
		  SCIP_CALL( SCIPaggregateVars(scip, consdata->var, consdata->vbdvar, 1.0, 1.0, 1.0, &infeasible, &redundant, &aggregated) );
		  assert(!infeasible);
		  ++(*naggrvars);

		  SCIP_CALL( SCIPdelCons(scip, cons) );
		  ++(*ndelconss);

		  continue;
	       }
	       assert(consdata->lhs < 0.5);

	       vars[0] = consdata->var;
	       vars[1] = consdata->vbdvar;
	    }
	    else
	    {
	       assert(SCIPisEQ(scip, consdata->lhs, 1.0));

	       SCIP_CALL( SCIPgetNegatedVar(scip, consdata->var, &vars[0]) );
	       SCIP_CALL( SCIPgetNegatedVar(scip, consdata->vbdvar, &vars[1]) );
	    }
	 }
	 /* the case x - y <= 0 or x - y >= 0 */
	 else
	 {
	    /* the case x - y <= 0 */
	    if( SCIPisZero(scip, consdata->rhs) )
	    {
	       /* check for aggregations like x - y == 0 */
	       if( SCIPisZero(scip, consdata->lhs) )
	       {
		  SCIP_Bool infeasible;
		  SCIP_Bool redundant;
		  SCIP_Bool aggregated;

		  SCIPdebugMsg(scip, "varbound constraint <%s>: aggregate <%s> - <%s> == 0\n",
		     SCIPconsGetName(cons), SCIPvarGetName(consdata->var), SCIPvarGetName(consdata->vbdvar));

		  /* aggregate both variables */
		  SCIP_CALL( SCIPaggregateVars(scip, consdata->var, consdata->vbdvar, 1.0, -1.0, 0.0, &infeasible, &redundant, &aggregated) );
		  assert(!infeasible);
		  ++(*naggrvars);

		  SCIP_CALL( SCIPdelCons(scip, cons) );
		  ++(*ndelconss);

		  continue;
	       }
	       assert(consdata->lhs < -0.5);

	       vars[0] = consdata->var;
	       SCIP_CALL( SCIPgetNegatedVar(scip, consdata->vbdvar, &vars[1]) );
	    }
	    /* the case x - y >= 0 */
	    else
	    {
	       assert(SCIPisZero(scip, consdata->lhs));

	       SCIP_CALL( SCIPgetNegatedVar(scip, consdata->var, &vars[0]) );
	       vars[1] = consdata->vbdvar;
	    }
	 }

	 SCIP_CALL( SCIPcreateConsSetpack(scip, &newcons, SCIPconsGetName(cons), 2, vars,
	       SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
	       SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
	       SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
	       SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

	 SCIP_CALL( SCIPaddCons(scip, newcons) );
	 SCIPdebugMsg(scip, "upgraded varbound constraint <%s> to a set-packing constraint\n", SCIPconsGetName(cons));
	 SCIPdebugPrintCons(scip, newcons, NULL);

	 SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
	 ++(*naddconss);

	 SCIP_CALL( SCIPdelCons(scip, cons) );
	 ++(*ndelconss);
      }
   }

   return SCIP_OKAY;
}

/**@} */


/**@name Linear constraint upgrading
 *
 */

/** tries to upgrade a linear constraint into a variable bound constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdVarbound)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to a variable bound constraint  lhs <= x + a*y <= rhs
    * - there are exactly two variables
    * - one of the variables is non-binary (called the bounded variable x)
    * - one of the variables is non-continuous (called the bounding variable y)
    */
   upgrade = (nvars == 2) && (nposbin + nnegbin <= 1) && (nposcont + nnegcont <= 1);

   if( upgrade )
   {
      SCIP_VAR* var;
      SCIP_VAR* vbdvar;
      SCIP_Real vbdcoef;
      SCIP_Real vbdlhs;
      SCIP_Real vbdrhs;
      int vbdind;

      /* decide which variable we want to use as bounding variable y */
      if( SCIPvarGetType(vars[0]) < SCIPvarGetType(vars[1]) )
         vbdind = 0;
      else if( SCIPvarGetType(vars[0]) > SCIPvarGetType(vars[1]) )
         vbdind = 1;
      else if( SCIPisIntegral(scip, vals[0]) && !SCIPisIntegral(scip, vals[1]) )
         vbdind = 0;
      else if( !SCIPisIntegral(scip, vals[0]) && SCIPisIntegral(scip, vals[1]) )
         vbdind = 1;
      else if( REALABS(REALABS(vals[0]) - 1.0) < REALABS(REALABS(vals[1]) - 1.0) )
         vbdind = 1;
      else
         vbdind = 0;

      /* do not upgrade when it is numerical unstable */
      if( SCIPisZero(scip, vals[vbdind]/vals[1-vbdind]) )
         return SCIP_OKAY;

      SCIPdebugMsg(scip, "upgrading constraint <%s> to variable bound constraint\n", SCIPconsGetName(cons));

      var = vars[1-vbdind];
      vbdvar = vars[vbdind];

      assert(!SCIPisZero(scip, vals[1-vbdind]));
      vbdcoef = vals[vbdind]/vals[1-vbdind];

      if( vals[1-vbdind] > 0.0 )
      {
         vbdlhs = SCIPisInfinity(scip, -lhs) ? -SCIPinfinity(scip) : lhs/vals[1-vbdind];
         vbdrhs = SCIPisInfinity(scip, rhs) ? SCIPinfinity(scip) : rhs/vals[1-vbdind];
      }
      else
      {
         vbdlhs = SCIPisInfinity(scip, rhs) ? -SCIPinfinity(scip) : rhs/vals[1-vbdind];
         vbdrhs = SCIPisInfinity(scip, -lhs) ? SCIPinfinity(scip) : lhs/vals[1-vbdind];
      }

      /* create the bin variable bound constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsVarbound(scip, upgdcons, SCIPconsGetName(cons), var, vbdvar, vbdcoef, vbdlhs, vbdrhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyVarbound)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* drop events */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( dropEvents(scip, cons, conshdlrdata->eventhdlr) );
   }

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata,
         sourcedata->var, sourcedata->vbdvar, sourcedata->vbdcoef, sourcedata->lhs, sourcedata->rhs) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch events for variables */
   SCIP_CALL( catchEvents(scip, *targetcons, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpVarbound)
{  /*lint --e{715}*/
   int i;

   *infeasible = FALSE;

   for( i = 0; i < nconss && !(*infeasible); i++ )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i], infeasible) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], conshdlrdata->usebdwidening, NULL, result) );
   }

   /* separate remaining constraints */
   for( i = nusefulconss; i < nconss && *result == SCIP_DIDNOTFIND; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], conshdlrdata->usebdwidening, NULL, result) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], conshdlrdata->usebdwidening, sol, result) );
   }

   /* separate remaining constraints */
   for( i = nusefulconss; i < nconss && *result == SCIP_DIDNOTFIND; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], conshdlrdata->usebdwidening, sol, result) );
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, FALSE) )
      {
         assert((*result) == SCIP_INFEASIBLE || (*result) == SCIP_FEASIBLE);
         (*result) = SCIP_INFEASIBLE;

         SCIP_CALL( SCIPresetConsAge(scip, conss[i]) );

         SCIP_CALL( separateCons(scip, conss[i], conshdlrdata->usebdwidening, NULL, result) );
         assert((*result) != SCIP_FEASIBLE);

         if( (*result) != SCIP_INFEASIBLE )
            break;
      }
      else
      {
         /* increase age of constraint */
         SCIP_CALL( SCIPincConsAge(scip, conss[i]) );
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], sol, FALSE) )
      {
         assert((*result) == SCIP_INFEASIBLE || (*result) == SCIP_FEASIBLE);
         (*result) = SCIP_INFEASIBLE;

         SCIP_CALL( SCIPresetConsAge(scip, conss[i]) );

         SCIP_CALL( separateCons(scip, conss[i], conshdlrdata->usebdwidening, sol, result) );
         assert((*result) != SCIP_FEASIBLE);

         if( (*result) != SCIP_INFEASIBLE )
            break;
      }
      else
      {
         /* increase age of constraint */
         SCIP_CALL( SCIPincConsAge(scip, conss[i]) );
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, TRUE) )
      {
         SCIP_CALL( SCIPresetConsAge(scip, conss[i]) );

         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
      else
      {
         /* increase age of constraint */
         SCIP_CALL( SCIPincConsAge(scip, conss[i]) );
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckVarbound)
{  /*lint --e{715}*/
   int i;

   *result = SCIP_FEASIBLE;

   for( i = 0; i < nconss && (*result == SCIP_FEASIBLE || completely); i++ )
   {
      if( !checkCons(scip, conss[i], sol, checklprows) )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            SCIP_CONSDATA* consdata;
            SCIP_Real sum;

            consdata = SCIPconsGetData(conss[i]);
            assert( consdata != NULL );

            sum = SCIPgetSolVal(scip, sol, consdata->var) + consdata->vbdcoef * SCIPgetSolVal(scip, sol, consdata->vbdvar);

            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");

            if( !SCIPisFeasGE(scip, sum, consdata->lhs) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhs - sum);
            }
            if( !SCIPisFeasLE(scip, sum, consdata->rhs) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", sum - consdata->rhs);
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int nchgbds;
   int nchgsides;
   int i;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   nchgbds = 0;

   SCIPdebugMsg(scip, "propagating %d variable bound constraints\n", nmarkedconss);

   /* process constraints marked for propagation */
   for( i = 0; i < nmarkedconss && !cutoff; i++ )
   {
      SCIP_CALL( propagateCons(scip, conss[i], conshdlrdata->usebdwidening, &cutoff, &nchgbds, &nchgsides, NULL) );
   }

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   int oldnchgbds;
   int oldndelconss;
   int oldnaddconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int oldnaggrvars;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnaddconss = *naddconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;
   oldnaggrvars = *naggrvars;

   for( i = 0; i < nconss; i++ )
   {
      cons = conss[i];
      assert(cons != NULL);

      assert(!SCIPconsIsModifiable(cons));

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( i % 1000 == 0 && SCIPisStopped(scip) )
         break;

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->presolved = FALSE;

      if( consdata->presolved )
         continue;
      consdata->presolved = TRUE;

      /* incorporate fixings and aggregations in constraint */
      SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr, &cutoff, nchgbds, ndelconss, naddconss) );

      if( cutoff )
         break;
      if( !SCIPconsIsActive(cons) )
         continue;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata->usebdwidening, &cutoff, nchgbds, nchgsides, ndelconss) );

      if( cutoff )
         break;
      if( !SCIPconsIsActive(cons) )
         continue;

      /* tighten variable bound coefficient */
      SCIP_CALL( tightenCoefs(scip, cons, nchgcoefs, nchgsides, ndelconss, &cutoff, nchgbds) );
      if( cutoff )
         break;
      if( !SCIPconsIsActive(cons) )
         continue;

      /* informs once variable x about a globally valid variable lower or upper bound */
      if( !consdata->varboundsadded )
      {
         SCIP_Bool infeasible;
         int nlocalchgbds;
         int localoldnchgbds;

         localoldnchgbds = *nchgbds;

         /* if lhs is finite, we have a variable lower bound: lhs <= x + c*y  =>  x >= -c*y + lhs */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIPdebugMsg(scip, "adding variable lower bound <%s> >= %g<%s> + %g (and potentially also <%s> %s %g<%s> + %g)\n",
               SCIPvarGetName(consdata->var), -consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs,
               SCIPvarGetName(consdata->vbdvar), (consdata->vbdcoef > 0 ? ">=" : "<="), 1.0/-consdata->vbdcoef,
               SCIPvarGetName(consdata->var), consdata->lhs/consdata->vbdcoef);

            SCIP_CALL( SCIPaddVarVlb(scip, consdata->var, consdata->vbdvar, -consdata->vbdcoef, consdata->lhs,
                  &infeasible, &nlocalchgbds) );
            assert(!infeasible);

            *nchgbds += nlocalchgbds;
         }

         /* if rhs is finite, we have a variable upper bound: x + c*y <= rhs  =>  x <= -c*y + rhs */
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIPdebugMsg(scip, "adding variable upper bound <%s> <= %g<%s> + %g (and potentially also <%s> %s %g<%s> + %g)\n",
               SCIPvarGetName(consdata->var), -consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
               SCIPvarGetName(consdata->vbdvar), (consdata->vbdcoef > 0 ? "<=" : ">="), 1.0/-consdata->vbdcoef,
               SCIPvarGetName(consdata->var), consdata->rhs/consdata->vbdcoef);

            SCIP_CALL( SCIPaddVarVub(scip, consdata->var, consdata->vbdvar, -consdata->vbdcoef, consdata->rhs,
                  &infeasible, &nlocalchgbds) );
            assert(!infeasible);

            *nchgbds += nlocalchgbds;
         }
         consdata->varboundsadded = TRUE;

         if( *nchgbds > localoldnchgbds )
         {
            /* tighten variable bound coefficient */
            SCIP_CALL( tightenCoefs(scip, cons, nchgcoefs, nchgsides, ndelconss, &cutoff, nchgbds) );
            if( cutoff )
               break;
         }
      }
   }

   if( !cutoff )
   {
      /* for varbound constraint with two integer variables make coefficients integral */
      SCIP_CALL( prettifyConss(scip, conss, nconss, nchgcoefs, nchgsides) );

      /* check if we can upgrade to a set-packing constraint */
      SCIP_CALL( upgradeConss(scip, conshdlrdata, conss, nconss, &cutoff, naggrvars, nchgbds, nchgcoefs, nchgsides, ndelconss, naddconss) );

      if( !cutoff && conshdlrdata->presolpairwise && (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 )
      {
	 /* preprocess pairs of variable bound constraints */
	 SCIP_CALL( preprocessConstraintPairs(scip, conss, nconss, &cutoff, nchgbds, ndelconss, nchgcoefs, nchgsides) );
      }
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nchgbds > oldnchgbds || *ndelconss > oldndelconss || *naddconss > oldnaddconss
      || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides || *naggrvars > oldnaggrvars )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, boundtype, bdchgidx, relaxedbd, conshdlrdata->usebdwidening) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->var, nlockspos, nlocksneg) );
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlockspos, nlocksneg) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlocksneg, nlockspos) );
      }
   }

   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->var, nlocksneg, nlockspos) );
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlocksneg, nlockspos) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlockspos, nlocksneg) );
      }
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   SCIPinfoMessage(scip, file, "<%s>[%c] %+.15g<%s>[%c]", SCIPvarGetName(consdata->var), 
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_BINARY ? SCIP_VARTYPE_BINARY_CHAR :
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_INTEGER ? SCIP_VARTYPE_INTEGER_CHAR :
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_IMPLINT ? SCIP_VARTYPE_IMPLINT_CHAR : SCIP_VARTYPE_CONTINUOUS_CHAR,
      consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar),
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_BINARY ? SCIP_VARTYPE_BINARY_CHAR :
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_INTEGER ? SCIP_VARTYPE_INTEGER_CHAR :
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_IMPLINT ? SCIP_VARTYPE_IMPLINT_CHAR : SCIP_VARTYPE_CONTINUOUS_CHAR);

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyVarbound)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   const char* consname;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, 2) );

   vars[0] = SCIPgetVarVarbound(sourcescip, sourcecons);
   vars[1] = SCIPgetVbdvarVarbound(sourcescip, sourcecons);

   coefs[0] = 1.0;
   coefs[1] = SCIPgetVbdcoefVarbound(sourcescip, sourcecons);

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   /* copy the varbound using the linear constraint copy method */
   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, 2, vars, coefs,
         SCIPgetLhsVarbound(sourcescip, sourcecons), SCIPgetRhsVarbound(sourcescip, sourcecons), varmap, consmap, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );

   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseVarbound)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real lhs;
   SCIP_Real rhs;
   char* endstr;
   int requiredsize;
   int nvars;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   (*success) = FALSE;

   /* return of string empty */
   if( !*str )
      return SCIP_OKAY;

   /* ignore whitespace */
   while( isspace(*str) )
      ++str;

   if( isdigit(str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit(str[1])) )
   {
      if( !SCIPstrToRealValue(str, &lhs, &endstr) )
      {
         SCIPerrorMessage("error parsing left hand side\n");
         return SCIP_OKAY;
      }

      /* ignore whitespace */
      while( isspace(*endstr) )
         ++endstr;

      if( endstr[0] != '<' || endstr[1] != '=' )
      {
         SCIPerrorMessage("missing \"<=\" after left hand side(, found %c%c)\n", endstr[0], endstr[1]);
         return SCIP_OKAY;
      }

      SCIPdebugMsg(scip, "found left hand side <%g>\n", lhs);

      /* it was indeed a left-hand-side, so continue parsing after it */
      str = endstr + 2;
   }

   /* pares x + c*y as linear sum */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, 2) );

   /* parse linear sum to get variables and coefficients */
   SCIP_CALL( SCIPparseVarsLinearsum(scip, str, vars, coefs, &nvars, 2, &requiredsize, &endstr, success) );

   if( requiredsize == 2 && *success )
   {
      SCIP_Bool foundvalue;
      SCIP_Real value;

      assert(nvars == 2);
      assert(SCIPisEQ(scip, coefs[0], 1.0));

      SCIPdebugMsg(scip, "found linear sum <%s> + %g <%s>\n", SCIPvarGetName(vars[0]), coefs[1], SCIPvarGetName(vars[1]));

      /* ignore whitespace */
      while( isspace(*endstr) )
         ++endstr;

      str = endstr;

      foundvalue = SCIPstrToRealValue(str+2, &value, &endstr);

      if( foundvalue )
      {
         /* search for end of linear sum: either '<=', '>=', '==', or '[free]' */
         switch( *str )
         {
         case '<':
            assert(str[1] == '=');
            rhs = value;
            break;
         case '=':
            assert(str[1] == '=');
            assert(SCIPisInfinity(scip, -lhs));
            lhs = value;
            rhs = value;
            break;
         case '>':
            assert(str[1] == '=');
            assert(SCIPisInfinity(scip, -lhs));
            lhs = value;
            break;
         default:
            SCIPerrorMessage("missing relation symbol after linear sum\n");
            *success = FALSE;
         }
      }
      else if( strncmp(str, "[free]", 6) != 0 )
         *success = FALSE;
   }

   if( *success )
   {
      SCIP_CALL( SCIPcreateConsVarbound(scip, cons, name, vars[0], vars[1], coefs[1], lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsVarbound)
{  /*lint --e{715}*/

   if( varssize < 2 )
      (*success) = FALSE;
   else
   {
      SCIP_CONSDATA* consdata;
      assert(cons != NULL);
      assert(vars != NULL);

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      vars[0] = consdata->var;
      vars[1] = consdata->vbdvar;
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsVarbound)
{  /*lint --e{715}*/
   (*nvars) = 2;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * Event Handler
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVarbound)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;

   assert(event != NULL);
   cons = (SCIP_CONS*)eventdata;
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_VARFIXED )
   {
      consdata->presolved = FALSE;
   }
   else
   {
      assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_BOUNDTIGHTENED) != 0);

      consdata->presolved = FALSE;
      consdata->tightened = FALSE;

      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** creates the handler for variable bound constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrVarbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_CONSHDLR* conshdlr;

   /* include event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecVarbound, NULL) );

   /* create variable bound constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpVarbound, consEnfopsVarbound, consCheckVarbound, consLockVarbound,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyVarbound, consCopyVarbound) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteVarbound) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolVarbound) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeVarbound) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsVarbound) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsVarbound) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpVarbound) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseVarbound) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolVarbound, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintVarbound) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropVarbound, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropVarbound) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpVarbound, consSepasolVarbound, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransVarbound) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxVarbound) );

   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint to varbound constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdVarbound, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }

   /* add varbound constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/maxlpcoef",
         "maximum coefficient in varbound constraint to be added as a row into LP",
         &conshdlrdata->maxlpcoef, TRUE, DEFAULT_MAXLPCOEF, 0.0, 1e+20, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/usebdwidening", "should bound widening be used in conflict analysis?",
         &conshdlrdata->usebdwidening, FALSE, DEFAULT_USEBDWIDENING, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a variable bound constraint: lhs <= x + c*y <= rhs
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs,                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs,                /**< right hand side of variable bound inequality */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   /* find the variable bound constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("variable bound constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, var, vbdvar, vbdcoef, lhs, rhs) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   if( SCIPisTransformed(scip) )
   {
      /* catch events for variables */
      SCIP_CALL( catchEvents(scip, *cons, conshdlrdata->eventhdlr) );
   }

   return SCIP_OKAY;
}

/** creates and captures a variable bound constraint: lhs <= x + c*y <= rhs
 *  with all constraint flags set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs,                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs                 /**< right hand side of variable bound inequality */
   )
{
   SCIP_CALL( SCIPcreateConsVarbound(scip, cons, name, var, vbdvar,vbdcoef, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** gets left hand side of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetLhsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetRhsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** gets bounded variable x of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_VAR* SCIPgetVarVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->var;
}

/** gets bounding variable y of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_VAR* SCIPgetVbdvarVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vbdvar;
}

/** gets bound coefficient c of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetVbdcoefVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vbdcoef;
}

/** gets the dual solution of the variable bound constraint in the current LP */
SCIP_Real SCIPgetDualsolVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual Farkas value of the variable bound constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given variable bound constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}

/**@} */
