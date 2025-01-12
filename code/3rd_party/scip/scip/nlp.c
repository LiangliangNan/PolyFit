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

/**@file   nlp.c
 * @ingroup OTHER_CFILES
 * @brief  NLP management methods
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 *
 *  In NLP management, we have to distinguish between the current NLP and the NLPI problem
 *  stored in the NLP solver. All NLP methods affect the current NLP only.
 *  Before solving the current NLP with the NLP solver, the NLP solvers data
 *  has to be updated to the current NLP with a call to SCIPnlpFlush().
 *
 *  @todo handle linear rows from LP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/nlpi.h"
#include "scip/pub_expr.h"
#include "scip/expr.h"
#include "scip/expr_varidx.h"
#include "scip/clock.h"
#include "scip/event.h"
#include "scip/nlp.h"
#include "scip/primal.h"
#include "scip/pub_event.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_nlp.h"
#include "scip/pub_var.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/struct_nlp.h"
/* to get nlp, set, ... in event handling and mapvar2varidx */
#include "scip/struct_scip.h"
/* to get value of parameter "nlp/solver" and nlpis array and to get access to set->lp for releasing a variable */
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/var.h"
#include <string.h>

/* defines */

#define EVENTHDLR_NAME   "nlpEventHdlr"      /**< name of NLP event handler that catches variable events */
#define EVENTHDLR_DESC   "handles all events necessary for maintaining NLP data"  /**< description of NLP event handler */
#define ADDNAMESTONLPI   0                   /**< whether to give variable and row names to NLPI */

/*lint -e440*/
/*lint -e441*/
/*lint -e777*/

#ifdef __cplusplus
extern "C" {
#endif

/* avoid inclusion of scip.h */ /*lint -e{2701}*/
BMS_BLKMEM* SCIPblkmem(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

/*
 * forward declarations
 */

/** NLP event handler execution method */
static
SCIP_DECL_EVENTEXEC( eventExecNlp );

/** announces, that a row of the NLP was modified
 *
 * adjusts status of current solution;
 * calling method has to ensure that change is passed on to the NLPI!
 */
static
SCIP_RETCODE nlpRowChanged(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row which was changed */
   );

/*
 * private NLP nonlinear row methods
 */

/** announces, that the given linear coefficient in the constraint matrix changed */
static
SCIP_RETCODE nlrowLinearCoefChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var,                /**< variable which coefficient changed */
   SCIP_Real             coef,               /**< new coefficient of variable, 0.0 if deleted */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);
   assert(var   != NULL);

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      /* update NLPI problem, if row is in NLPI already */
      if( nlrow->nlpiindex >= 0 )
      {
         int idx;

         /* get index of variable in NLPI */
         assert(SCIPhashmapExists(nlp->varhash, var));
         idx = SCIPhashmapGetImageInt(nlp->varhash, var);
         assert(idx >= 0 && idx < nlp->nvars);

         idx = nlp->varmap_nlp2nlpi[idx];
         assert(idx >= 0 && idx < nlp->nvars_solver);

         /* change coefficient in NLPI problem */
         SCIP_CALL( SCIPnlpiChgLinearCoefs(set, nlp->solver, nlp->problem, nlrow->nlpiindex, 1, &idx, &coef) );
      }
   }

   return SCIP_OKAY;
}

/** create varidx expression for var expression
 *
 * called when expr is duplicated for addition to NLPI
 */
static
SCIP_DECL_EXPR_MAPEXPR(mapvar2varidx)
{
   SCIP_NLP* nlp;
   int nlpidx;

   assert(sourcescip != NULL);
   assert(sourcescip == targetscip);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);
   assert(*targetexpr == NULL);
   assert(mapexprdata != NULL);

   nlp = (SCIP_NLP*)mapexprdata;

   /* do not provide map if not variable */
   if( !SCIPexprIsVar(sourcescip->set, sourceexpr) )
      return SCIP_OKAY;

   assert(SCIPvarIsActive(SCIPgetVarExprVar(sourceexpr)));  /* because we simplified exprs */

   assert(SCIPhashmapExists(nlp->varhash, SCIPgetVarExprVar(sourceexpr)));
   nlpidx = SCIPhashmapGetImageInt(nlp->varhash, SCIPgetVarExprVar(sourceexpr));
   assert(nlpidx < nlp->nvars);

   assert(nlp->varmap_nlp2nlpi[nlpidx] >= 0);
   assert(nlp->varmap_nlp2nlpi[nlpidx] < nlp->nvars_solver);
   SCIP_CALL( SCIPcreateExprVaridx(targetscip, targetexpr, nlp->varmap_nlp2nlpi[nlpidx], ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/** announces, that an expression changed */
static
SCIP_RETCODE nlrowExprChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      if( nlrow->nlpiindex >= 0 )
      {
         /* change expression tree in NLPI problem */
         SCIP_EXPR* nlpiexpr;

         SCIP_CALL( SCIPexprCopy(set, stat, blkmem, set, stat, blkmem, nlrow->expr, &nlpiexpr, mapvar2varidx, (void*)nlp, NULL, NULL) );
         SCIP_CALL( SCIPnlpiChgExpr(set, nlp->solver, nlp->problem, nlrow->nlpiindex, nlpiexpr) );
         SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &nlpiexpr) );
      }
   }

   return SCIP_OKAY;
}

/** notifies nonlinear row, that its sides were changed */
static
SCIP_RETCODE nlrowSideChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      if( nlrow->nlpiindex >= 0 )
      {
         SCIP_Real lhs;
         SCIP_Real rhs;

         /* change sides in NLPI problem */
         lhs = nlrow->lhs;
         rhs = nlrow->rhs;
         if( !SCIPsetIsInfinity(set, -lhs) )
            lhs -= nlrow->constant;
         if( !SCIPsetIsInfinity(set,  rhs) )
            rhs -= nlrow->constant;

         SCIP_CALL( SCIPnlpiChgConsSides(set, nlp->solver, nlp->problem, 1, &nlrow->nlpiindex, &lhs, &rhs) );
      }
   }

   return SCIP_OKAY;
}

/** notifies nonlinear row, that its constant was changed */
static
SCIP_RETCODE nlrowConstantChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      if( nlrow->nlpiindex >= 0 )
      {
         SCIP_Real lhs;
         SCIP_Real rhs;

         lhs = nlrow->lhs;
         rhs = nlrow->rhs;
         if( !SCIPsetIsInfinity(set, -lhs) )
            lhs -= nlrow->constant;
         if( !SCIPsetIsInfinity(set,  rhs) )
            rhs -= nlrow->constant;

         /* change sides in NLPI problem */
         SCIP_CALL( SCIPnlpiChgConsSides(set, nlp->solver, nlp->problem, 1, &nlrow->nlpiindex, &lhs, &rhs) );
      }
   }

   return SCIP_OKAY;
}

/** sorts linear part of row entries such that lower variable indices precede higher ones */
static
void nlrowSortLinear(
   SCIP_NLROW*           nlrow               /**< nonlinear row to be sorted */
   )
{
   assert(nlrow != NULL);

   /* check, if row is already sorted in the LP part, or if the sorting should be delayed */
   if( nlrow->linvarssorted )
      return;

   /* sort linear coefficients */
   SCIPsortPtrReal((void**)nlrow->linvars, nlrow->lincoefs, SCIPvarComp, nlrow->nlinvars);

   nlrow->linvarssorted = TRUE;
}

/** searches linear variable in nonlinear row, returns position in linvars vector or -1 if not found */
static
int nlrowSearchLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be searched in */
   SCIP_VAR*             var                 /**< variable to be searched for */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);

   if( nlrow->nlinvars == 0 )
      return -1;

   nlrowSortLinear(nlrow);
   if( !SCIPsortedvecFindPtr((void**)nlrow->linvars, SCIPvarComp, (void*)var, nlrow->nlinvars, &pos) )
      return -1;

   return pos;
}

/** moves a coefficient in a nonlinear row to a different place, and updates all corresponding data structures */
static
void nlrowMoveLinearCoef(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   oldpos,             /**< old position of coefficient */
   int                   newpos              /**< new position of coefficient */
   )
{
   assert(nlrow != NULL);
   assert(0 <= oldpos && oldpos < nlrow->nlinvars);
   assert(0 <= newpos && newpos < nlrow->nlinvars);
   assert(nlrow->linvars[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   nlrow->linvars[newpos]  = nlrow->linvars[oldpos];
   nlrow->lincoefs[newpos] = nlrow->lincoefs[oldpos];

   /* update sorted flags */
   nlrow->linvarssorted = FALSE;
}

/** adds a previously non existing linear coefficient to a nonlinear row */
static
SCIP_RETCODE nlrowAddLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value of coefficient */
   )
{
   int pos;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(var    != NULL);
   assert(coef   != 0.0);

   /* assert that only active variables are added once the row is in the NLP */
   assert(nlrow->nlpindex == -1 || SCIPvarIsActive(var) );

   SCIP_CALL( SCIPnlrowEnsureLinearSize(nlrow, blkmem, set, nlrow->nlinvars+1) );
   assert(nlrow->linvars  != NULL);
   assert(nlrow->lincoefs != NULL);

   pos = nlrow->nlinvars;
   nlrow->nlinvars++;

   /* insert the variable */
   nlrow->linvars [pos] = var;
   nlrow->lincoefs[pos] = coef;

   SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, var, coef, nlp) );

   /* update sorted flag */
   if( pos > 0 && SCIPvarCompare(nlrow->linvars[pos-1], nlrow->linvars[pos]) > 0 )
      nlrow->linvarssorted = FALSE;

   SCIPsetDebugMsg(set, "added linear coefficient %g * <%s> at position %d to nonlinear row <%s>\n",
      coef, SCIPvarGetName(var), pos, nlrow->name);

   return SCIP_OKAY;
}

#ifdef SCIP_DISABLED_CODE
/** adds a linear coefficient to a nonlinear row
 * if the variable exists in the linear part of the row already, the coefficients are added
 * otherwise the variable is added to the row */
static
SCIP_RETCODE nlrowAddToLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef,               /**< value of coefficient */
   SCIP_Bool             removefixed         /**< whether to disaggregate var before adding */
   )
{
   int pos;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(var    != NULL);

   if( removefixed && !SCIPvarIsActive(var) )
   {
      SCIP_Real constant;

      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &coef, &constant) );
      if( constant != 0.0 )
      {
         nlrow->constant += constant;
         SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
      }

      if( SCIPsetIsZero(set, coef) )
         return SCIP_OKAY;

      if( !SCIPvarIsActive(var) )
      {
         int j;

         /* if var is still not active, then it is multi-aggregated */
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

         if( SCIPvarGetMultaggrConstant(var) != 0.0 )
         {
            nlrow->constant += coef * SCIPvarGetMultaggrConstant(var);
            SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
         }

         for( j = 0; j < SCIPvarGetMultaggrNVars(var); ++j )
         {
            SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, SCIPvarGetMultaggrVars(var)[j], SCIPvarGetMultaggrScalars(var)[j] * coef, TRUE) );
         }

         return SCIP_OKAY;
      }
   }
   else if( SCIPsetIsZero(set, coef) )
      return SCIP_OKAY;

   assert(!removefixed || SCIPvarIsActive(var));

   pos = nlrowSearchLinearCoef(nlrow, var);

   if( pos == -1 )
   {
      /* add as new coefficient */
      SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, var, coef) );
   }
   else
   {
      assert(pos >= 0);
      assert(pos <  nlrow->nlinvars);
      assert(nlrow->linvars[pos] == var);

      /* add to previously existing coefficient */
      nlrow->lincoefs[pos] += coef;
   }

   return SCIP_OKAY;
}
#endif

/** deletes coefficient at given position from row */
static
SCIP_RETCODE nlrowDelLinearCoefPos(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos                 /**< position in row vector to delete */
   )
{
   SCIP_VAR* var;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars[pos] != NULL);

   var = nlrow->linvars[pos];

   /* move last coefficient to position of empty slot (should set sorted flag to FALSE, if not last variable was deleted) */
   nlrowMoveLinearCoef(nlrow, nlrow->nlinvars-1, pos);
   nlrow->nlinvars--;
   assert(pos == nlrow->nlinvars || nlrow->linvarssorted == FALSE);

   SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, var, 0.0, nlp) );

   return SCIP_OKAY;
}

/** changes a coefficient at given position of a nonlinear row */
static
SCIP_RETCODE nlrowChgLinearCoefPos(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos,                /**< position in row vector to change */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   assert(nlrow != NULL);
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars != NULL);
   assert(nlrow->linvars[pos] != NULL);

   if( SCIPsetIsZero(set, coef) )
   {
      /* delete existing coefficient */
      SCIP_CALL( nlrowDelLinearCoefPos(nlrow, set, stat, nlp, pos) );
   }
   else if( !SCIPsetIsEQ(set, nlrow->lincoefs[pos], coef) )
   {
      /* change existing coefficient */
      nlrow->lincoefs[pos] = coef;
      SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, nlrow->linvars[pos], coef, nlp) );
   }

   return SCIP_OKAY;
}

/** calculates minimal and maximal activity of row w.r.t. the variable's bounds */
static
SCIP_RETCODE nlrowCalcActivityBounds(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   SCIP_Real inf;
   SCIP_INTERVAL activity;
   SCIP_INTERVAL bounds;
   int i;

   assert(nlrow != NULL);
   assert(set   != NULL);
   assert(stat  != NULL);

   inf = SCIPsetInfinity(set);

   /* calculate activity bounds */
   SCIPintervalSet(&activity, nlrow->constant);
   for( i = 0; i < nlrow->nlinvars && !SCIPintervalIsEntire(inf, activity); ++i )
   {
      SCIPintervalSetBounds(&bounds, SCIPvarGetLbLocal(nlrow->linvars[i]), SCIPvarGetUbLocal(nlrow->linvars[i]));
      SCIPintervalMulScalar(inf, &bounds, bounds, nlrow->lincoefs[i]);
      SCIPintervalAdd(inf, &activity, activity, bounds);
   }

   if( nlrow->expr != NULL && !SCIPintervalIsEntire(inf, activity) )
   {
      SCIP_CALL( SCIPexprEvalActivity(set, stat, blkmem, nlrow->expr) );
      SCIPintervalAdd(inf, &activity, activity, SCIPexprGetActivity(nlrow->expr));
   }

   nlrow->minactivity = SCIPintervalGetInf(activity);
   nlrow->maxactivity = SCIPintervalGetSup(activity);

   nlrow->validactivitybdsdomchg = stat->domchgcount;

   return SCIP_OKAY;
}

/** makes sure that there is no fixed variable at position pos of the linear part of a nonlinear row
 *
 * a fixed variable is replaced with the corresponding constant or disaggregated term
 */
static
SCIP_RETCODE nlrowRemoveFixedLinearCoefPos(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos                 /**< position of variable in linear variables array */
   )
{
   SCIP_Real oldconstant;
   SCIP_VAR* var;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(pos >= 0);
   assert(pos <  nlrow->nlinvars);

   var = nlrow->linvars[pos];

   if( SCIPvarIsActive(var) )
      return SCIP_OKAY;

   oldconstant = nlrow->constant;

   /* replace fixed, aggregated, or negated variable */
   SCIP_CALL( SCIPvarGetProbvarSum( &nlrow->linvars[pos], set, &nlrow->lincoefs[pos], &nlrow->constant) );

   /* if var had been fixed, entry should be removed from row */
   if( nlrow->lincoefs[pos] == 0.0 )
   {
      nlrowMoveLinearCoef(nlrow, nlrow->nlinvars-1, pos);
      nlrow->nlinvars--;

      if( pos < nlrow->nlinvars )
      {
         SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, pos) );
      }

      return SCIP_OKAY;
   }
   nlrow->linvarssorted = FALSE;

   /* notify nlrow that coefficient of var is now 0.0 in row */
   SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, var, 0.0, nlp) );

   /* notify nlrow that constant of row has changed */
   if( oldconstant != nlrow->constant )
      SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );

   if( SCIPvarIsActive(nlrow->linvars[pos]) )
   {
      /* if var was aggregated or negated, notify nlrow about new coefficient */
      SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, nlrow->linvars[pos], nlrow->lincoefs[pos], nlp) );
   }
   else
   {
      SCIP_Real coef;
      int i;

      /* if not removed or active, the new variable should be multi-aggregated */
      assert(SCIPvarGetStatus(nlrow->linvars[pos]) == SCIP_VARSTATUS_MULTAGGR);

      var  = nlrow->linvars[pos];
      coef = nlrow->lincoefs[pos];

      /* remove the variable from the row */
      SCIP_CALL( nlrowDelLinearCoefPos(nlrow, set, stat, nlp, pos) );

      /* add multi-aggregated term to row */
      if( SCIPvarGetMultaggrConstant(var) != 0.0 )
      {
         nlrow->constant += coef * SCIPvarGetMultaggrConstant(var);
         SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
      }
      SCIP_CALL( SCIPnlrowEnsureLinearSize(nlrow, blkmem, set, nlrow->nlinvars + SCIPvarGetMultaggrNVars(var)) );
      for( i = 0; i < SCIPvarGetMultaggrNVars(var); ++i )
      {
         if( SCIPsetIsZero(set, coef * SCIPvarGetMultaggrScalars(var)[i]) )
            continue;
         SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, SCIPvarGetMultaggrVars(var)[i], coef * SCIPvarGetMultaggrScalars(var)[i]) );
         assert(SCIPvarGetMultaggrVars(var)[i] == nlrow->linvars[nlrow->nlinvars-1]);
         if( !SCIPvarIsActive(SCIPvarGetMultaggrVars(var)[i]) )
         {
            /* if newly added variable is fixed, replace it now */
            SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, nlrow->nlinvars-1) );
         }
      }

      /* due to nlrowDelLinearCoefPos, an inactive variable may have moved to position pos
       * if that is the case, call ourself recursively
       */
      if( pos < nlrow->nlinvars && !SCIPvarIsActive(nlrow->linvars[pos]) )
      {
         SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, pos) );
      }
   }

   return SCIP_OKAY;
}

/** removes fixed variables from the linear part of a nonlinear row */
static
SCIP_RETCODE nlrowRemoveFixedLinearCoefs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   int i;
   int oldlen;

   assert(nlrow != NULL);
   assert(nlrow->linvars != NULL || nlrow->nlinvars == 0);

   oldlen = nlrow->nlinvars;
   for( i = 0; i < MIN(oldlen, nlrow->nlinvars); ++i )
   {
      assert(nlrow->linvars[i] != NULL);
      SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, i) );
   }

   return SCIP_OKAY;
}

/** removes fixed variables from expression of a nonlinear row */
static
SCIP_RETCODE nlrowSimplifyExpr(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   SCIP_EXPR* simplified;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   if( nlrow->expr == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPexprSimplify(set, stat, blkmem, nlrow->expr, &simplified, &changed, &infeasible, NULL, NULL) );
   assert(!infeasible);

   if( !changed )
   {
      SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &simplified) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &nlrow->expr) );
   nlrow->expr = simplified;

   if( SCIPexprIsValue(set, nlrow->expr) )
   {
      /* if expression tree is constant, remove it */
      SCIP_CALL( SCIPnlrowChgConstant(nlrow, set, stat, nlp, nlrow->constant + SCIPgetValueExprValue(nlrow->expr)) );

      SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &nlrow->expr) );
   }

   SCIP_CALL( nlrowExprChanged(nlrow, blkmem, set, stat, nlp) );

   return SCIP_OKAY;
}

/** removes fixed variable from nonlinear row */
static
SCIP_RETCODE nlrowRemoveFixedVar(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var                 /**< variable that had been fixed */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);
   assert(!SCIPvarIsActive(var));

   /* search for variable in linear part and remove if existing */
   pos = nlrowSearchLinearCoef(nlrow, var);
   if( pos >= 0 )
   {
      SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, pos) );
   }

   /* search for variable in nonlinear part and remove all fixed variables in expression if existing
    * TODO only call simplify if var appears in expr, but currently we don't store the vars in a separate array
    */
   if( nlrow->expr != NULL )
   {
      SCIP_CALL( nlrowSimplifyExpr(nlrow, blkmem, set, stat, nlp) );
   }

   return SCIP_OKAY;
}

/*
 * public NLP nonlinear row methods
 */
/**@addtogroup PublicNLRowMethods
 *
 * @{
 */

/** create a new nonlinear row
 *
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreate(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             constant,           /**< constant */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< linear variables, or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< linear coefficients, or NULL if nlinvars == 0 */
   SCIP_EXPR*            expr,               /**< expression, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_EXPRCURV         curvature           /**< curvature of the nonlinear row */
   )
{
#ifndef NDEBUG
   int i;
#endif

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(name   != NULL);
   assert(!SCIPsetIsInfinity(set, ABS(constant)));
   assert(nlinvars   == 0 || linvars   != NULL);
   assert(nlinvars   == 0 || lincoefs  != NULL);
   assert(SCIPsetIsRelLE(set, lhs, rhs));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, nlrow) );

   /* constant part */
   assert(!SCIPsetIsInfinity(set, REALABS(constant)));
   (*nlrow)->constant = constant;

#ifndef NDEBUG
   for( i = 0; i < nlinvars; ++i )
   {
      assert(linvars[i] != NULL);
      assert(!SCIPsetIsInfinity(set, REALABS(lincoefs[i])));
   }
#endif

   /* linear part */
   (*nlrow)->nlinvars = nlinvars;
   (*nlrow)->linvarssize = nlinvars;
   if( nlinvars > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->linvars,  linvars,  nlinvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->lincoefs, lincoefs, nlinvars) );
      (*nlrow)->linvarssorted = FALSE;
   }
   else
   {
      (*nlrow)->linvars  = NULL;
      (*nlrow)->lincoefs = NULL;
      (*nlrow)->linvarssorted = TRUE;
   }

   /* nonlinear part */
   if( expr != NULL )
   {
      /* TODO preserve common subexpressions, or at least use only one varexpr per var */
      SCIP_CALL( SCIPexprCopy(set, stat, blkmem, set, stat, blkmem, expr, &(*nlrow)->expr, NULL, NULL, NULL, NULL) );
   }
   else
   {
      (*nlrow)->expr = NULL;
   }

   /* left and right hand sides, asserted above that lhs is relatively less equal than rhs */
   (*nlrow)->lhs = MIN(lhs, rhs);
   (*nlrow)->rhs = MAX(lhs, rhs);

   /* miscellaneous */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->name, name, strlen(name)+1) );
   (*nlrow)->activity = SCIP_INVALID;
   (*nlrow)->validactivitynlp = FALSE;
   (*nlrow)->pseudoactivity = SCIP_INVALID;
   (*nlrow)->validpsactivitydomchg = FALSE;
   (*nlrow)->minactivity = SCIP_INVALID;
   (*nlrow)->maxactivity = SCIP_INVALID;
   (*nlrow)->validactivitybdsdomchg = FALSE;
   (*nlrow)->nlpindex = -1;
   (*nlrow)->nlpiindex = -1;
   (*nlrow)->nuses = 0;
   (*nlrow)->dualsol = 0.0;
   (*nlrow)->curvature = curvature;

   /* capture the nonlinear row */
   SCIPnlrowCapture(*nlrow);

   return SCIP_OKAY;
}

/** create a nonlinear row that is a copy of a given row
 *
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreateCopy(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           sourcenlrow         /**< nonlinear row to copy */
   )
{
   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(sourcenlrow != NULL);

   SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, stat, sourcenlrow->name,
         sourcenlrow->constant,
         sourcenlrow->nlinvars, sourcenlrow->linvars, sourcenlrow->lincoefs,
         sourcenlrow->expr,
         sourcenlrow->lhs, sourcenlrow->rhs, sourcenlrow->curvature) );

   (*nlrow)->linvarssorted          = sourcenlrow->linvarssorted;
   (*nlrow)->activity               = sourcenlrow->activity;
   (*nlrow)->validactivitynlp       = sourcenlrow->validactivitynlp;
   (*nlrow)->pseudoactivity         = sourcenlrow->pseudoactivity;
   (*nlrow)->validpsactivitydomchg  = sourcenlrow->validpsactivitydomchg;
   (*nlrow)->minactivity            = sourcenlrow->minactivity;
   (*nlrow)->maxactivity            = sourcenlrow->maxactivity;
   (*nlrow)->validactivitybdsdomchg = sourcenlrow->validactivitybdsdomchg;

   return SCIP_OKAY;
}

/** create a new nonlinear row from a linear row
 *
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreateFromRow(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_ROW*             row                 /**< the linear row to copy */
   )
{
   int rownz;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(row    != NULL);

   rownz = SCIProwGetNNonz(row);

   if( rownz > 1 )
   {
      SCIP_VAR** rowvars;
      int i;

      SCIP_CALL( SCIPsetAllocBufferArray(set, &rowvars, rownz) );

      for( i = 0; i < rownz; ++i )
      {
         rowvars[i] = SCIPcolGetVar(SCIProwGetCols(row)[i]);
         assert(rowvars[i] != NULL);
      }

      SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, stat, SCIProwGetName(row),
            SCIProwGetConstant(row),
            rownz, rowvars, SCIProwGetVals(row), NULL,
            SCIProwGetLhs(row), SCIProwGetRhs(row),
            SCIP_EXPRCURV_LINEAR) );

      SCIPsetFreeBufferArray(set, &rowvars);
   }
   else if( rownz == 1 )
   {
      SCIP_VAR* rowvar;

      rowvar = SCIPcolGetVar(SCIProwGetCols(row)[0]);

      SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, stat, SCIProwGetName(row),
            SCIProwGetConstant(row),
            1, &rowvar, SCIProwGetVals(row), NULL,
            SCIProwGetLhs(row), SCIProwGetRhs(row),
            SCIP_EXPRCURV_LINEAR) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, stat, SCIProwGetName(row),
            SCIProwGetConstant(row),
            0, NULL, NULL, NULL,
            SCIProwGetLhs(row), SCIProwGetRhs(row),
            SCIP_EXPRCURV_LINEAR) );
   }

   return SCIP_OKAY;   
}

/** output nonlinear row to file stream */
SCIP_RETCODE SCIPnlrowPrint(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   assert(nlrow != NULL);

   /* print row name */
   if( nlrow->name != NULL && nlrow->name[0] != '\0' )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, "%s: ", nlrow->name);
   }

   /* print left hand side */
   SCIPmessageFPrintInfo(messagehdlr, file, "%.15g <= ", nlrow->lhs);

   /* print constant */
   SCIPmessageFPrintInfo(messagehdlr, file, "%.15g ", nlrow->constant);

   /* print linear coefficients */
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);
      assert(SCIPvarGetName(nlrow->linvars[i]) != NULL);
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.15g<%s> ", nlrow->lincoefs[i], SCIPvarGetName(nlrow->linvars[i]));
   }

   /* print nonlinear part */
   if( nlrow->expr != NULL )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, " + ");
      SCIP_CALL( SCIPexprPrint(set, stat, blkmem, messagehdlr, file, nlrow->expr) );
   }

   /* print right hand side */
   SCIPmessageFPrintInfo(messagehdlr, file, " <= %.15g\n", nlrow->rhs);

   return SCIP_OKAY;
}

/** increases usage counter of nonlinear row */
void SCIPnlrowCapture(
   SCIP_NLROW*           nlrow               /**< nonlinear row to capture */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->nuses >= 0);

   SCIPdebugMessage("capture nonlinear row <%s> with nuses=%d\n", nlrow->name, nlrow->nuses);
   nlrow->nuses++;
}

/** decreases usage counter of nonlinear row */
SCIP_RETCODE SCIPnlrowRelease(
   SCIP_NLROW**          nlrow,              /**< nonlinear row to free */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(blkmem != NULL);
   assert(nlrow  != NULL);
   assert(*nlrow != NULL);
   assert((*nlrow)->nuses >= 1);

   SCIPsetDebugMsg(set, "release nonlinear row <%s> with nuses=%d\n", (*nlrow)->name, (*nlrow)->nuses);
   (*nlrow)->nuses--;
   if( (*nlrow)->nuses > 0 )
   {
      *nlrow = NULL;
      return SCIP_OKAY;
   }

   /* free row */

   assert((*nlrow)->nuses == 0);
   assert((*nlrow)->nlpindex == -1);
   assert((*nlrow)->nlpiindex == -1);

   /* linear part */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->linvars,   (*nlrow)->linvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->lincoefs,  (*nlrow)->linvarssize);

   /* nonlinear part */
   if( (*nlrow)->expr != NULL )
   {
      SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &(*nlrow)->expr) );
   }

   /* miscellaneous */
   BMSfreeBlockMemoryArray(blkmem, &(*nlrow)->name, strlen((*nlrow)->name)+1);

   BMSfreeBlockMemory(blkmem, nlrow);

   return SCIP_OKAY;
}

/** ensures, that linear coefficient array of nonlinear row can store at least num entries */
SCIP_RETCODE SCIPnlrowEnsureLinearSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->nlinvars <= nlrow->linvarssize);

   if( num > nlrow->linvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->linvars,  nlrow->linvarssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->lincoefs, nlrow->linvarssize, newsize) );
      nlrow->linvarssize = newsize;
   }
   assert(num <= nlrow->linvarssize);

   return SCIP_OKAY;
}

/** adds a previously non existing linear coefficient to a nonlinear row */
SCIP_RETCODE SCIPnlrowAddLinearCoef(
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   /* if row is in NLP already, make sure that only active variables are added */
   if( nlrow->nlpindex >= 0 )
   {
      SCIP_Real constant;

      /* get corresponding active or multi-aggregated variable */
      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &val, &constant) );

      /* add constant */
      SCIP_CALL( SCIPnlrowChgConstant(nlrow, set, stat, nlp, nlrow->constant + constant) );

      if( val == 0.0 )
         /* var has been fixed */
         return SCIP_OKAY;

      if( !SCIPvarIsActive(var) )
      {
         /* var should be multi-aggregated, so call this function recursively */
         int i;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
         for( i = 0; i < SCIPvarGetMultaggrNVars(var); ++i )
         {
            SCIP_CALL( SCIPnlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, SCIPvarGetMultaggrVars(var)[i], SCIPvarGetMultaggrScalars(var)[i] * val) );
         }
         return SCIP_OKAY;
      }

      /* var is active, so can go on like normal */
   }

   SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, var, val) );

   return SCIP_OKAY;
}

/** deletes linear coefficient from nonlinear row */
SCIP_RETCODE SCIPnlrowDelLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);

   /* if the row is in the NLP already, we can only have active variables, so var should also be active; in non-debug mode, one gets an error below */
   assert(nlrow->nlpindex == -1 || SCIPvarIsActive(var) );

   /* search the position of the variable in the row's variable vector */
   pos = nlrowSearchLinearCoef(nlrow, var);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for variable <%s> doesn't exist in nonlinear row <%s>\n", SCIPvarGetName(var), nlrow->name);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars[pos] == var);

   /* delete the variable from the row's variable vector */
   SCIP_CALL( nlrowDelLinearCoefPos(nlrow, set, stat, nlp, pos) );

   return SCIP_OKAY;
}

/** changes or adds a linear coefficient to a nonlinear row */
SCIP_RETCODE SCIPnlrowChgLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(nlp != NULL);
   assert(var != NULL);

   /* search the position of the variable in the row's linvars vector */
   pos = nlrowSearchLinearCoef(nlrow, var);

   /* check, if column already exists in the row's linear variables vector */
   if( pos == -1 )
   {
      if( !SCIPsetIsZero(set, coef) )
      {
         /* add previously not existing coefficient */
         SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, var, coef) );
      }
   }
   else
   {
      /* change the coefficient in the row */
      SCIP_CALL( nlrowChgLinearCoefPos(nlrow, set, stat, nlp, pos, coef) );
   }

   return SCIP_OKAY;
}

/** replaces or deletes an expression in a nonlinear row */
SCIP_RETCODE SCIPnlrowChgExpr(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_EXPR*            expr                /**< new expression */
   )
{
   assert(nlrow  != NULL);
   assert(blkmem != NULL);

   /* free previous expression tree */
   if( nlrow->expr != NULL )
   {
      SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &nlrow->expr) );
      assert(nlrow->expr == NULL);
   }

   /* adds new expression tree */
   if( expr != NULL )
   {
      /* TODO preserve common subexpressions, or at least use only one varexpr per var */
      SCIP_CALL( SCIPexprCopy(set, stat, blkmem, set, stat, blkmem, expr, &nlrow->expr, NULL, NULL, NULL, NULL) );

      /* if row is already in NLP, ensure that expr has only active variables */
      if( nlrow->nlpindex >= 0 )
      {
         SCIP_EXPR* simplified;
         SCIP_Bool changed;
         SCIP_Bool infeasible;

         SCIP_CALL( SCIPexprSimplify(set, stat, blkmem, nlrow->expr, &simplified, &changed, &infeasible, NULL, NULL) );
         assert(!infeasible);

         SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &nlrow->expr) );
         nlrow->expr = simplified;
      }
   }

   /* notify row about the change */
   SCIP_CALL( nlrowExprChanged(nlrow, blkmem, set, stat, nlp) );

   return SCIP_OKAY;
}

/** changes constant of nonlinear row */
SCIP_RETCODE SCIPnlrowChgConstant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             constant            /**< new constant */
   )
{
   assert(nlrow != NULL);

   if( !SCIPsetIsEQ(set, nlrow->constant, constant) )
   {
      nlrow->constant = constant;
      SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
   }

   return SCIP_OKAY;
}

/** changes left hand side of nonlinear row */
SCIP_RETCODE SCIPnlrowChgLhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(nlrow != NULL);

   if( !SCIPsetIsEQ(set, nlrow->lhs, lhs) )
   {
      nlrow->lhs = lhs;
      SCIP_CALL( nlrowSideChanged(nlrow, set, stat, nlp) );
   }

   return SCIP_OKAY;
}

/** changes right hand side of nonlinear row */
SCIP_RETCODE SCIPnlrowChgRhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   assert(nlrow != NULL);

   if( !SCIPsetIsEQ(set, nlrow->rhs, rhs) )
   {
      nlrow->rhs = rhs;
      SCIP_CALL( nlrowSideChanged(nlrow, set, stat, nlp) );
   }

   return SCIP_OKAY;
}

/** removes (or substitutes) all fixed, negated, aggregated, multi-aggregated variables from the linear and nonlinear part of a nonlinear row and simplifies its expression */
SCIP_RETCODE SCIPnlrowSimplify(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   SCIP_CALL( nlrowRemoveFixedLinearCoefs(nlrow, blkmem, set, stat, nlp) );
   SCIP_CALL( nlrowSimplifyExpr(nlrow, blkmem, set, stat, nlp) );

   return SCIP_OKAY;
}

/** recalculates the current activity of a nonlinear row in the current NLP solution */
SCIP_RETCODE SCIPnlrowRecalcNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   int i;

   assert(nlrow  != NULL);
   assert(stat   != NULL);
   assert(nlp    != NULL);

   if( nlp->solstat > SCIP_NLPSOLSTAT_LOCINFEASIBLE )
   {
      SCIPerrorMessage("do not have NLP solution for computing NLP activity\n");
      return SCIP_ERROR;
   }

   nlrow->activity = nlrow->constant;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);
      assert(SCIPvarGetNLPSol(nlrow->linvars[i]) < SCIP_INVALID);

      nlrow->activity += nlrow->lincoefs[i] * SCIPvarGetNLPSol(nlrow->linvars[i]);
   }

   if( nlrow->expr != NULL )
   {
      SCIP_SOL* sol;

      SCIP_CALL( SCIPsolCreateNLPSol(&sol, blkmem, set, stat, primal, tree, nlp, NULL) );

      SCIP_CALL( SCIPexprEval(set, stat, blkmem, nlrow->expr, sol, 0L) );
      if( SCIPexprGetEvalValue(nlrow->expr) == SCIP_INVALID )
         nlrow->activity = SCIP_INVALID;
      else
         nlrow->activity += SCIPexprGetEvalValue(nlrow->expr);

      SCIP_CALL( SCIPsolFree(&sol, blkmem, primal) );
   }

   nlrow->validactivitynlp = stat->nnlps;

   return SCIP_OKAY;
}

/** gives the activity of a nonlinear row in the current NLP solution */
SCIP_RETCODE SCIPnlrowGetNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   assert(nlrow  != NULL);
   assert(stat   != NULL);
   assert(activity != NULL);

   assert(nlrow->validactivitynlp <= stat->nnlps);

   if( nlrow->validactivitynlp != stat->nnlps )
   {
      SCIP_CALL( SCIPnlrowRecalcNLPActivity(nlrow, blkmem, set, stat, primal, tree, nlp) );
   }
   assert(nlrow->validactivitynlp == stat->nnlps);
   assert(nlrow->activity < SCIP_INVALID);

   *activity = nlrow->activity;

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row in the current NLP solution: negative value means infeasibility */
SCIP_RETCODE SCIPnlrowGetNLPFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_Real activity;

   assert(nlrow != NULL);
   assert(feasibility != NULL);

   SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, blkmem, set, stat, primal, tree, nlp, &activity) );
   *feasibility = MIN(nlrow->rhs - activity, activity - nlrow->lhs);

   return SCIP_OKAY;
}

/** calculates the current pseudo activity of a nonlinear row */
SCIP_RETCODE SCIPnlrowRecalcPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< SCIP problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< SCIP LP */
   )
{
   SCIP_Real val1;
   int i;

   assert(nlrow  != NULL);
   assert(stat   != NULL);

   nlrow->pseudoactivity = nlrow->constant;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);

      val1 = SCIPvarGetBestBoundLocal(nlrow->linvars[i]);
      nlrow->pseudoactivity += nlrow->lincoefs[i] * val1;
   }

   if( nlrow->expr != NULL )
   {
      SCIP_SOL* sol;

      SCIP_CALL( SCIPsolCreatePseudoSol(&sol, blkmem, set, stat, prob, primal, tree, lp, NULL) );

      SCIP_CALL( SCIPexprEval(set, stat, blkmem, nlrow->expr, sol, 0L) );
      if( SCIPexprGetEvalValue(nlrow->expr) == SCIP_INVALID )
         nlrow->pseudoactivity = SCIP_INVALID;
      else
         nlrow->pseudoactivity += SCIPexprGetEvalValue(nlrow->expr);

      SCIP_CALL( SCIPsolFree(&sol, blkmem, primal) );
   }

   nlrow->validpsactivitydomchg = stat->domchgcount;

   return SCIP_OKAY;
}

/** returns the pseudo activity of a nonlinear row in the current pseudo solution */
SCIP_RETCODE SCIPnlrowGetPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< SCIP problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< SCIP LP */
   SCIP_Real*            pseudoactivity      /**< buffer to store pseudo activity value */
   )
{
   assert(nlrow != NULL);
   assert(stat  != NULL);
   assert(pseudoactivity != NULL);
   assert(nlrow->validpsactivitydomchg <= stat->domchgcount);

   /* check, if pseudo activity has to be calculated */
   if( nlrow->validpsactivitydomchg != stat->domchgcount )
   {
      SCIP_CALL( SCIPnlrowRecalcPseudoActivity(nlrow, blkmem, set, stat, prob, primal, tree, lp) );
   }
   assert(nlrow->validpsactivitydomchg == stat->domchgcount);
   assert(nlrow->pseudoactivity < SCIP_INVALID);

   *pseudoactivity = nlrow->pseudoactivity;

   return SCIP_OKAY;
}

/** returns the pseudo feasibility of a nonlinear row in the current pseudo solution: negative value means infeasibility */
SCIP_RETCODE SCIPnlrowGetPseudoFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< SCIP problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< SCIP LP */
   SCIP_Real*            pseudofeasibility   /**< buffer to store pseudo feasibility value */
   )
{
   SCIP_Real pseudoactivity;

   assert(nlrow != NULL);
   assert(stat  != NULL);
   assert(pseudofeasibility != NULL);

   SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, blkmem, set, stat, prob, primal, tree, lp, &pseudoactivity) );
   *pseudofeasibility = MIN(nlrow->rhs - pseudoactivity, pseudoactivity - nlrow->lhs);

   return SCIP_OKAY;
}

/** returns the activity of a nonlinear row for a given solution */
SCIP_RETCODE SCIPnlrowGetSolActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   SCIP_Real inf;
   SCIP_Real val1;
   int i;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(activity != NULL);

   *activity = nlrow->constant;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);

      val1 = SCIPsolGetVal(sol, set, stat, nlrow->linvars[i]);
      if( val1 == SCIP_UNKNOWN )
      {
         *activity = SCIP_INVALID;
         return SCIP_OKAY;
      }
      *activity += nlrow->lincoefs[i] * val1;
   }

   if( nlrow->expr != NULL )
   {
      SCIP_CALL( SCIPexprEval(set, stat, blkmem, nlrow->expr, sol, 0L) );
      if( SCIPexprGetEvalValue(nlrow->expr) == SCIP_INVALID )
         *activity = SCIP_INVALID;
      else
         *activity += SCIPexprGetEvalValue(nlrow->expr);
   }

   inf = SCIPsetInfinity(set);
   *activity = MAX(*activity, -inf);
   *activity = MIN(*activity, +inf);

   return SCIP_OKAY;
}

/** returns the feasibility of a nonlinear row for the given solution */
SCIP_RETCODE SCIPnlrowGetSolFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_Real activity;

   assert(nlrow != NULL);
   assert(feasibility != NULL);

   SCIP_CALL( SCIPnlrowGetSolActivity(nlrow, blkmem, set, stat, sol, &activity) );

   *feasibility = MIN(nlrow->rhs - activity, activity - nlrow->lhs);

   return SCIP_OKAY;
}

/** returns the minimal activity of a nonlinear row w.r.t. the variables' bounds */
SCIP_RETCODE SCIPnlrowGetActivityBounds(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Real*            minactivity,        /**< buffer to store minimal activity, or NULL */
   SCIP_Real*            maxactivity         /**< buffer to store maximal activity, or NULL */
   )
{
   assert(nlrow != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(nlrow->validactivitybdsdomchg <= stat->domchgcount);

   /* check, if activity bounds has to be calculated */
   if( nlrow->validactivitybdsdomchg != stat->domchgcount )
   {
      SCIP_CALL( nlrowCalcActivityBounds(nlrow, blkmem, set, stat) );
   }
   assert(nlrow->validactivitybdsdomchg == stat->domchgcount);
   assert(nlrow->minactivity < SCIP_INVALID);
   assert(nlrow->maxactivity < SCIP_INVALID);

   if( minactivity != NULL )
      *minactivity = nlrow->minactivity;
   if( maxactivity != NULL )
      *maxactivity = nlrow->maxactivity;

   return SCIP_OKAY;
}

/** returns whether the nonlinear row is redundant w.r.t. the variables' bounds */
SCIP_RETCODE SCIPnlrowIsRedundant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Bool*            isredundant         /**< buffer to store whether row is redundant */
   )
{
   SCIP_Real minactivity;
   SCIP_Real maxactivity;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(isredundant != NULL);

   SCIP_CALL( SCIPnlrowGetActivityBounds(nlrow, blkmem, set, stat, &minactivity, &maxactivity) );

   *isredundant = TRUE;
   if( (!SCIPsetIsInfinity(set, -nlrow->lhs) && SCIPsetIsFeasLT(set, minactivity, nlrow->lhs)) ||
      ( !SCIPsetIsInfinity(set,  nlrow->rhs) && SCIPsetIsFeasGT(set, maxactivity, nlrow->rhs)) )
      *isredundant = FALSE;

   return SCIP_OKAY;
}

#ifdef NDEBUG
/* Undo the defines from pub_nlhdlr.h, which exist if NDEBUG is defined. */
#undef SCIPnlrowGetConstant
#undef SCIPnlrowGetNLinearVars
#undef SCIPnlrowGetLinearVars
#undef SCIPnlrowGetLinearCoefs
#undef SCIPnlrowGetExpr
#undef SCIPnlrowGetLhs
#undef SCIPnlrowGetRhs
#undef SCIPnlrowGetCurvature
#undef SCIPnlrowSetCurvature
#undef SCIPnlrowGetName
#undef SCIPnlrowGetNLPPos
#undef SCIPnlrowIsInNLP
#undef SCIPnlrowGetDualsol
#endif

/** gets constant */
SCIP_Real SCIPnlrowGetConstant(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->constant;
}

/** gets number of variables of linear part */
int SCIPnlrowGetNLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlinvars;
}

/** gets array with variables of linear part */
SCIP_VAR** SCIPnlrowGetLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->linvars;
}

/** gets array with coefficients in linear part */
SCIP_Real* SCIPnlrowGetLinearCoefs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->lincoefs;
}

/** gets expression */
SCIP_EXPR* SCIPnlrowGetExpr(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->expr;
}

/** returns the left hand side of a nonlinear row */
SCIP_Real SCIPnlrowGetLhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->lhs;
}

/** returns the right hand side of a nonlinear row */
SCIP_Real SCIPnlrowGetRhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->rhs;
}

/** returns the curvature of a nonlinear row */
SCIP_EXPRCURV SCIPnlrowGetCurvature(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);
   return nlrow->curvature;
}

/** sets the curvature of a nonlinear row */
void SCIPnlrowSetCurvature(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_EXPRCURV         curvature           /**< curvature of NLP row */
   )
{
   assert(nlrow != NULL);
   nlrow->curvature = curvature;
}

/** returns the name of a nonlinear row */
const char* SCIPnlrowGetName(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->name;
}

/** gets position of a nonlinear row in current NLP, or -1 if not in NLP */
int SCIPnlrowGetNLPPos(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlpindex;
}

/** returns TRUE iff row is member of current NLP */
SCIP_Bool SCIPnlrowIsInNLP(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlpindex != -1;
}

/** gets the dual NLP solution of a nlrow
 *
 * for a ranged constraint, the dual value is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_Real SCIPnlrowGetDualsol(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlpiindex >= 0 ? nlrow->dualsol : 0.0;
}

/** @} */

/*
 * local NLP methods
 */

/** announces, that a row of the NLP was modified
 * adjusts status of current solution
 * calling method has to ensure that change is passed to the NLPI!
 */ /*lint -e{715}*/
static
SCIP_RETCODE nlpRowChanged(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row which was changed */
   )
{  /*lint --e{715}*/
   assert(nlp != NULL);
   assert(nlrow != NULL);
   assert(!nlp->indiving);
   assert(nlrow->nlpindex >= 0);

   /* nlrow is a row in the NLP, so changes effect feasibility */
   /* if we have a feasible NLP solution and it satisfies the modified row, then it is still feasible
    * if the NLP was globally or locally infeasible or unbounded, then this may not be the case anymore
    */
   if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      /* TODO bring this back? then everything that may call nlpRowChanged will need to pass on blkmem, primal, tree as well
      SCIP_Real feasibility;
      SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, blkmem, set, stat, primal, tree, nlp, &feasibility) );
      if( !SCIPsetIsFeasNegative(set, feasibility) )
         nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
      else */
      nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
   }
   else
   {
      nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
   }

   return SCIP_OKAY;
}

/** adds a set of nonlinear rows to the NLP and captures them */
static
SCIP_RETCODE nlpAddNlRows(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   nnlrows,            /**< number of nonlinear rows to add */
   SCIP_NLROW**          nlrows              /**< nonlinear rows to add */
   )
{
#ifndef NDEBUG
   int i;
#endif
   int j;
   SCIP_NLROW* nlrow;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlrows != NULL || nnlrows == 0);
   assert(!nlp->indiving);

   SCIP_CALL( SCIPnlpEnsureNlRowsSize(nlp, blkmem, set, nlp->nnlrows + nnlrows) );

   for( j = 0; j < nnlrows; ++j )
   {
      nlrow = nlrows[j];  /*lint !e613*/

      /* assert that row is not in NLP (or even NLPI) yet */
      assert(nlrow->nlpindex == -1);
      assert(nlrow->nlpiindex == -1);

      /* make sure there are only active variables in row */
      SCIP_CALL( SCIPnlrowSimplify(nlrow, blkmem, set, stat, nlp) );

#ifndef NDEBUG
      /* assert that variables of row are in NLP */
      for( i = 0; i < nlrow->nlinvars; ++i )
         assert(SCIPhashmapExists(nlp->varhash, nlrow->linvars[i]));

      if( nlrow->expr != NULL )
      {
         SCIP_EXPRITER* it;
         SCIP_EXPR* expr;

         SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
         SCIP_CALL( SCIPexpriterInit(it, nlrow->expr, SCIP_EXPRITER_DFS, TRUE) );
         for( expr = nlrow->expr; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
            assert(!SCIPexprIsVar(set, expr) || SCIPhashmapExists(nlp->varhash, SCIPgetVarExprVar(expr)));
         SCIPexpriterFree(&it);
      }
#endif

      /* add row to NLP and capture it */
      nlp->nlrows[nlp->nnlrows + j] = nlrow;
      nlrow->nlpindex = nlp->nnlrows + j;

      SCIPnlrowCapture(nlrow);

      /* if we have a feasible NLP solution and it satisfies the new solution, then it is still feasible
       * if the NLP was globally or locally infeasible, then it stays that way
       * if the NLP was unbounded, then this may not be the case anymore
       */
      if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         /* TODO bring this back? then everything that may call nlpAddNlRows will need to pass on primal, tree as well
         SCIP_Real feasibility;
         SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, blkmem, set, stat, primal, tree, nlp, &feasibility) );
         if( !SCIPsetIsFeasNegative(set, feasibility) )
            nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
         */
         nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
      }
      else if( nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
      {
         nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
      }
   }

   nlp->nnlrows += nnlrows;
   nlp->nunflushednlrowadd += nnlrows;

   return SCIP_OKAY;
}

/** moves a nonlinear row to a different place, and updates all corresponding data structures */
static
void nlpMoveNlrow(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   int                   oldpos,             /**< old position of nonlinear row */
   int                   newpos              /**< new position of nonlinear row */
   )
{
   assert(nlp != NULL);
   assert(0 <= oldpos && oldpos < nlp->nnlrows);
   assert(0 <= newpos && newpos < nlp->nnlrows);
   assert(nlp->nlrows[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   nlp->nlrows[newpos] = nlp->nlrows[oldpos];
   nlp->nlrows[newpos]->nlpindex = newpos;

   /* update nlpi to nlp row index mapping */
   if( nlp->nlrows[newpos]->nlpiindex >= 0 )
   {
      assert(nlp->nlrowmap_nlpi2nlp != NULL);
      assert(nlp->nlrows[newpos]->nlpiindex < nlp->sizenlrows_solver);
      nlp->nlrowmap_nlpi2nlp[nlp->nlrows[newpos]->nlpiindex] = newpos;
   }
}

/** deletes nonlinear row with given position from NLP */
static
SCIP_RETCODE nlpDelNlRowPos(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   pos                 /**< position of nonlinear row that is to be removed */
   )
{
   SCIP_NLROW* nlrow;

   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(pos >= 0);
   assert(pos < nlp->nnlrows);
   assert(!nlp->indiving);
   assert(nlp->nlrows != NULL);

   nlrow = nlp->nlrows[pos];
   assert(nlrow != NULL);
   assert(nlrow->nlpindex == pos);

   /* if row is in NLPI, then mark that it has to be removed in the next flush
    * if row was not in NLPI yet, then we have one unflushed nlrow addition less */
   if( nlrow->nlpiindex >= 0 )
   {
      assert(nlrow->nlpiindex < nlp->nnlrows_solver);
      nlp->nlrowmap_nlpi2nlp[nlrow->nlpiindex] = -1;
      nlrow->nlpiindex = -1;
      ++nlp->nunflushednlrowdel;
   }
   else
   {
      assert(nlrow->nlpiindex == -1);
      --nlp->nunflushednlrowadd;
   }

   /* move NLP row from the end to pos and mark nlrow to be not in NLP anymore */
   nlpMoveNlrow(nlp, nlp->nnlrows-1, pos);
   nlrow->nlpindex = -1;

   /* forget about restriction */
   SCIP_CALL( SCIPnlrowRelease(&nlrow, blkmem, set, stat) );
   --nlp->nnlrows;

   if( nlp->solstat < SCIP_NLPSOLSTAT_LOCOPT )
      nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
   else if( nlp->solstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE )
      nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;

   return SCIP_OKAY; /*lint !e438*/
}

/** updates bounds on a variable in the NLPI problem */
static
SCIP_RETCODE nlpUpdateVarBounds(
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable which bounds have changed */
   SCIP_Bool             tightened           /**< whether the bound change was a bound tightening */
   )
{
   int pos;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));

   /* original variable bounds are ignored during diving
    * (all variable bounds are reset to their current value in exitDiving) */
   if( nlp->indiving )
      return SCIP_OKAY;

   /* get position of variable in NLP */
   pos = SCIPhashmapGetImageInt(nlp->varhash, var);

   /* if variable not in NLPI yet, nothing to do */
   if( nlp->varmap_nlp2nlpi[pos] == -1 )
      return SCIP_OKAY;

   /* update bounds in NLPI problem */
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   pos = nlp->varmap_nlp2nlpi[pos];
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   SCIP_CALL( SCIPnlpiChgVarBounds(set, nlp->solver, nlp->problem, 1, &pos, &lb, &ub) );

   /* if we have a feasible NLP solution and it satisfies the new bounds, then it is still feasible
    * if the NLP was globally or locally infeasible and we tightened a bound, then it stays that way
    * if the NLP was unbounded and we tightened a bound, then this may not be the case anymore
    */
   if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      if( !tightened ||
         ((SCIPsetIsInfinity(set, -lb) || SCIPsetIsFeasLE(set, lb, SCIPvarGetNLPSol(var))) &&
          (SCIPsetIsInfinity(set,  ub) || SCIPsetIsFeasGE(set, ub, SCIPvarGetNLPSol(var)))) )
         nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
      else
         nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
   }
   else if( !tightened || nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
   {
      nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
   }

   return SCIP_OKAY;
}

/** updates coefficient of a variable in the objective */
static
SCIP_RETCODE nlpUpdateObjCoef(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_VAR*             var                 /**< variable which bounds have changed */
   )
{
   int pos;
   int objidx;
   SCIP_Real coef;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));

   /* if the objective in the NLPI is not up to date, then we do not need to do something here */
   if( !nlp->objflushed )
      return SCIP_OKAY;

   /* original objective is ignored during diving
    * we just need to remember that at end of diving we have to flush the objective */
   if( nlp->indiving )
   {
      nlp->objflushed = FALSE;
      return SCIP_OKAY;
   }

   /* get position of variable in NLP and objective coefficient */
   pos  = SCIPhashmapGetImageInt(nlp->varhash, var);
   assert(nlp->varmap_nlp2nlpi[pos] == -1 || nlp->solver != NULL);

   /* actually we only need to remember flushing the objective if we also have an NLPI */
   if( nlp->solver == NULL )
      return SCIP_OKAY;

   coef = SCIPvarGetObj(var);

   /* if variable not in NLPI yet, then we only need to remember to update the objective after variable additions were flushed */
   if( nlp->varmap_nlp2nlpi[pos] == -1 && coef != 0.0 )
   {
      nlp->objflushed = FALSE;

      return SCIP_OKAY;
   }

   /* if we are here, then the objective in the NLPI is up to date,
    * we keep it this way by changing the coefficient of var in the NLPI problem objective */
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   pos = nlp->varmap_nlp2nlpi[pos];
   objidx = -1;
   SCIP_CALL( SCIPnlpiChgLinearCoefs(set, nlp->solver, nlp->problem, objidx, 1, &pos, &coef) );

   /* if we had a solution and it was locally (or globally) optimal, then now we can only be sure that it is still feasible */
   if( nlp->solstat < SCIP_NLPSOLSTAT_FEASIBLE )
      nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;

   return SCIP_OKAY;
}

/** adds new variables to the NLP */
static
SCIP_RETCODE nlpAddVars(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables to add */
   SCIP_VAR**            vars                /**< variable to add to NLP */
   )
{
   int i;
   SCIP_VAR* var;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(vars   != NULL || nvars == 0);
   assert(!nlp->indiving || nvars == 0);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPnlpEnsureVarsSize(nlp, blkmem, set, nlp->nvars + nvars) );
   assert(nlp->sizevars >= nlp->nvars + nvars);

   for( i = 0; i < nvars; ++i )
   {
      var = vars[i];  /*lint !e613*/

      assert(SCIPvarIsTransformed(var));
      assert(SCIPvarIsActive(var));
      assert(!SCIPhashmapExists(nlp->varhash, var));

      SCIPvarCapture(var);

      nlp->vars[nlp->nvars+i]            = var;
      nlp->varmap_nlp2nlpi[nlp->nvars+i] = -1;
      SCIP_CALL( SCIPhashmapInsertInt(nlp->varhash, var, nlp->nvars+i) );

      nlp->varlbdualvals[nlp->nvars+i]   = 0.0;
      nlp->varubdualvals[nlp->nvars+i]   = 0.0;

      /* update objective, if necessary (new variables have coefficient 0.0 anyway) */
      if( SCIPvarGetObj(var) != 0.0 )
      {
         SCIP_CALL( nlpUpdateObjCoef(set, nlp, var) );
      }

      /* let's keep the previous initial guess and set it for the new variable to the best bound
       * (since there can be no row that uses this variable yet, this seems a good guess) */
      if( nlp->haveinitguess )
      {
         assert(nlp->initialguess != NULL);

         nlp->initialguess[nlp->nvars+i] = SCIPvarGetBestBoundLocal(var);
      }

      /* if we have a feasible NLP solution, then it remains feasible
       * but we have to update the objective function
       */
      if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         SCIP_CALL( SCIPvarSetNLPSol(var, set, SCIPvarGetBestBoundLocal(var)) );
         nlp->primalsolobjval += SCIPvarGetObj(var) * SCIPvarGetBestBoundLocal(var);
         nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
      }

      /* catch events on variable */
      SCIP_CALL( SCIPvarCatchEvent(var, blkmem, set, \
            SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_OBJCHANGED, \
            nlp->eventhdlr, (SCIP_EVENTDATA*)nlp, NULL) ); /* @todo should store event filter position in nlp? */
   }

   nlp->nvars += nvars;
   nlp->nunflushedvaradd += nvars;

   return SCIP_OKAY;
}

/** moves a variable to a different place, and updates all corresponding data structures */
static
SCIP_RETCODE nlpMoveVar(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   int                   oldpos,             /**< old position of variable */
   int                   newpos              /**< new position of variable */
   )
{
   int nlpipos;

   assert(nlp != NULL);
   assert(0 <= oldpos && oldpos < nlp->nvars);
   assert(0 <= newpos && newpos < nlp->nvars);
   assert(nlp->vars[oldpos] != NULL);

   if( oldpos == newpos )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapSetImageInt(nlp->varhash, nlp->vars[oldpos], newpos) );
   nlp->vars[newpos]            = nlp->vars[oldpos];
   nlp->varmap_nlp2nlpi[newpos] = nlp->varmap_nlp2nlpi[oldpos];
   nlp->varlbdualvals[newpos]   = nlp->varlbdualvals[oldpos];
   nlp->varubdualvals[newpos]   = nlp->varubdualvals[oldpos];
   if( nlp->initialguess != NULL )
      nlp->initialguess[newpos] = nlp->initialguess[oldpos];

   nlpipos = nlp->varmap_nlp2nlpi[newpos];
   if( nlpipos > 0 )
      nlp->varmap_nlpi2nlp[nlpipos] = newpos;

   return SCIP_OKAY;
}

/** deletes variable with given position from NLP */
static
SCIP_RETCODE nlpDelVarPos(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< SCIP LP, needed if a column-variable is freed */
   int                   pos                 /**< position of nonlinear row that is to be removed */
   )
{
   SCIP_VAR* var;
#ifndef NDEBUG
   int i;
#endif
   int nlpipos;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(pos >= 0);
   assert(pos < nlp->nvars);
   assert(!nlp->indiving);

   var = nlp->vars[pos];
   assert(var != NULL);

#ifndef NDEBUG
   /* assert that variable is not used by any nonlinear row */
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      int j;
      SCIP_NLROW* nlrow;

      nlrow = nlp->nlrows[i];
      assert(nlrow != NULL);

      /* use nlrowSearchLinearCoef only if already sorted, since otherwise we may change the solving process slightly */
      if( nlrow->linvarssorted )
         assert( nlrowSearchLinearCoef(nlrow, var) == -1 );
      else
         for( j = 0; j < nlrow->nlinvars; ++j )
            assert( nlrow->linvars[j] != var );

      if( nlrow->expr != NULL )
      {
         SCIP_EXPRITER* it;
         SCIP_EXPR* expr;
         SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
         SCIP_CALL( SCIPexpriterInit(it, nlrow->expr, SCIP_EXPRITER_DFS, TRUE) );
         for( expr = nlrow->expr; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
            assert(!SCIPexprIsVar(set, expr) || SCIPgetVarExprVar(expr) != var);
         SCIPexpriterFree(&it);
      }
   }
#endif

   /* if we had a feasible solution, then adjust objective function value
    * if NLP was unbounded before, then maybe it is not anymore */
   if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      nlp->primalsolobjval -= SCIPvarGetObj(var) * SCIPvarGetNLPSol(var);
   else if( nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
      nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;

   /* if variable is in NLPI problem, mark that we have to remember to delete it there
    * if it was not in the NLPI yet, then we have one unflushed var addition less now */
   nlpipos = nlp->varmap_nlp2nlpi[pos];
   if( nlpipos >= 0 )
   {
      assert(nlpipos < nlp->nvars_solver);

      nlp->varmap_nlpi2nlp[nlpipos] = -1;
      ++nlp->nunflushedvardel;
   }
   else
      --nlp->nunflushedvaradd;

   /* drop events on variable */
   SCIP_CALL( SCIPvarDropEvent(var, blkmem, set, \
         SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_OBJCHANGED, \
         nlp->eventhdlr, (SCIP_EVENTDATA*)nlp, -1) );

   /* move variable from end to pos */
   SCIP_CALL( nlpMoveVar(nlp, nlp->nvars-1, pos) );

   /* forget about variable */
   SCIP_CALL( SCIPhashmapRemove(nlp->varhash, var) );
   SCIP_CALL( SCIPvarRelease(&var, blkmem, set, eventqueue, lp) );
   --nlp->nvars;

   return SCIP_OKAY;
}

/** notifies NLP that a variable was fixed, so it is removed from objective, all rows, and the NLP variables */
static
SCIP_RETCODE nlpRemoveFixedVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< SCIP LP, needed to release variable */
   SCIP_VAR*             var                 /**< variable that has been fixed */
   )
{
   int i;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(!SCIPvarIsActive(var));
   assert(!nlp->indiving);
   assert(SCIPhashmapExists(nlp->varhash, var));

   /* remove var from all rows */
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      SCIP_CALL( nlrowRemoveFixedVar(nlp->nlrows[i], blkmem, set, stat, nlp, var) );
   }

   /* remove variable from NLP */
   SCIP_CALL( SCIPnlpDelVar(nlp, blkmem, set, stat, eventqueue, lp, var) );

   return SCIP_OKAY;
}

/** creates arrays with NLPI variable indices of linear variables in a nonlinear row */
static
SCIP_RETCODE nlpSetupNlpiIndices(
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   int**                 linidxs             /**< buffer to store pointer to NLPI indices of linear variables */
   )
{
   int i;
   SCIP_VAR* var;

   assert(nlp    != NULL);
   assert(set    != NULL);
   assert(nlrow  != NULL);
   assert(linidxs   != NULL);

   /* get indices of variables in linear part of row */
   if( nlrow->nlinvars > 0 )
   {
      assert(nlrow->linvars  != NULL);
      assert(nlrow->lincoefs != NULL);

      SCIP_CALL( SCIPsetAllocBufferArray(set, linidxs, nlrow->nlinvars) );

      for( i = 0; i < nlrow->nlinvars; ++i )
      {
         var = nlrow->linvars[i];
         assert(var != NULL);
         assert(SCIPvarIsActive(var)); /* at this point, there should be only active variables in the row */

         assert(SCIPhashmapExists(nlp->varhash, var));
         (*linidxs)[i] = nlp->varmap_nlp2nlpi[SCIPhashmapGetImageInt(nlp->varhash, var)];
         assert((*linidxs)[i] >= 0);
      }
   }
   else
      *linidxs = NULL;

   return SCIP_OKAY;
}

/** ensures, that NLPI variables array of NLP can store at least num entries */
static
SCIP_RETCODE nlpEnsureVarsSolverSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nvars_solver <= nlp->sizevars_solver);

   if( num > nlp->sizevars_solver )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varmap_nlpi2nlp, nlp->sizevars_solver, newsize) );

      nlp->sizevars_solver = newsize;
   }
   assert(num <= nlp->sizevars_solver);

   return SCIP_OKAY;
}

/** ensures, that NLPI nonlinear rows array of NLP can store at least num entries */
static
SCIP_RETCODE nlpEnsureNlRowsSolverSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nnlrows_solver <= nlp->sizenlrows_solver);

   if( num > nlp->sizenlrows_solver )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->nlrowmap_nlpi2nlp, nlp->sizenlrows_solver, newsize) );

      nlp->sizenlrows_solver = newsize;
   }
   assert(num <= nlp->sizenlrows_solver);

   return SCIP_OKAY;
}

/** deletes rows from the NLPI problem that have been marked as to remove */
static
SCIP_RETCODE nlpFlushNlRowDeletions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int         j;
   int         c;      /* counts the number of rows to delete */
   int*        rowset; /* marks which rows to delete and stores new indices */
   SCIP_NLROW* nlrow;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushednlrowdel >= 0);
   assert(!nlp->indiving);

   if( nlp->nunflushednlrowdel == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending removals of nonlinear rows */
      for( j = 0; j < nlp->nnlrows_solver; ++j )
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
#endif
      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* create marker which rows have to be deleted */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rowset, nlp->nnlrows_solver) );
   c = 0;
   for( j = 0; j < nlp->nnlrows_solver; ++j )
   {
      if( nlp->nlrowmap_nlpi2nlp[j] == -1 )
      {
         rowset[j] = 1;
         ++c;
      }
      else
         rowset[j] = 0;
   }
   assert(c == nlp->nunflushednlrowdel);

   /* remove rows from NLPI problem */
   SCIP_CALL( SCIPnlpiDelConsSet(set, nlp->solver, nlp->problem, rowset, nlp->nnlrows_solver) );

   /* update NLPI row indices */
   for( j = 0; j < nlp->nnlrows_solver; ++j )
   {
      assert(rowset[j] <= j); /* we assume that the NLP solver did not move a row behind its previous position!! */
      if( rowset[j] < 0 )
      {
         /* assert that row was marked as deleted */
         assert(nlp->nlrowmap_nlpi2nlp[j] == -1);
      }
      else if( rowset[j] < j )
      {
         /* nlrow at position j moved (forward) to position rowset[j] */
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
         assert(nlp->nlrowmap_nlpi2nlp[j] < nlp->nnlrows);

         nlrow = nlp->nlrows[nlp->nlrowmap_nlpi2nlp[j]];
         assert(nlrow->nlpiindex == j);

         /* there should be no row at the new position already */
         assert(nlp->nlrowmap_nlpi2nlp[rowset[j]] == -1);

         nlrow->nlpiindex = rowset[j];
         nlp->nlrowmap_nlpi2nlp[rowset[j]] = nlrow->nlpindex;
      }
      else
      {
         /* row j stays at position j */
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
         assert(nlp->nlrowmap_nlpi2nlp[j] < nlp->nnlrows);
         assert(nlp->nlrows[nlp->nlrowmap_nlpi2nlp[j]]->nlpiindex == j);
      }
   }
   nlp->nnlrows_solver -= c;
   nlp->nunflushednlrowdel = 0;

   /* cleanup */
   SCIPsetFreeBufferArray(set, &rowset);

   return SCIP_OKAY;
}

/** deletes variables from the NLPI problem that have been marked as to remove
 *
 * assumes that there are no pending row deletions (nlpFlushNlRowDeletions() should be called first)
 */
static
SCIP_RETCODE nlpFlushVarDeletions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int  i;
   int  c;      /* counter on number of variables to remove in solver */
   int* colset; /* marks which variables to delete and stores new indices */

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushedvardel >= 0);
   assert(nlp->nunflushednlrowdel == 0);
   assert(!nlp->indiving);

   if( nlp->nunflushedvardel == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending removals of variables */
      for( i = 0; i < nlp->nvars_solver; ++i )
         assert(nlp->varmap_nlpi2nlp[i] >= 0);
#endif
      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* create marker which variables have to be deleted */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &colset, nlp->nvars_solver) );
   c = 0;
   for( i = 0; i < nlp->nvars_solver; ++i )
   {
      if( nlp->varmap_nlpi2nlp[i] == -1 )
      {
         colset[i] = 1;
         ++c;
      }
      else
         colset[i] = 0;
   }
   assert(c == nlp->nunflushedvardel);

   /* delete variables from NLPI problem */
   SCIP_CALL( SCIPnlpiDelVarSet(set, nlp->solver, nlp->problem, colset, nlp->nvars_solver) );

   /* update NLPI variable indices */
   for( i = 0; i < nlp->nvars_solver; ++i )
   {
      assert(colset[i] <= i); /* we assume that the NLP solver did not move a variable behind its previous position!! */
      if( colset[i] < 0 )
      {
         /* assert that variable was marked as deleted */
         assert(nlp->varmap_nlpi2nlp[i] == -1);
      }
      else if( colset[i] < i)
      {
         /* variable at position i moved (forward) to position colset[i] */
         int varpos;

         varpos = nlp->varmap_nlpi2nlp[i]; /* position of variable i in NLP */
         assert(varpos >= 0);
         assert(varpos < nlp->nvars);
         assert(nlp->varmap_nlp2nlpi[varpos] == i);

         /* there should be no variable at the new position already */
         assert(nlp->varmap_nlpi2nlp[colset[i]] == -1);

         nlp->varmap_nlp2nlpi[varpos] = colset[i];
         nlp->varmap_nlpi2nlp[colset[i]] = varpos;
      }
      else
      {
         /* variable i stays at position i */
         assert(nlp->varmap_nlpi2nlp[i] >= 0);
         assert(nlp->varmap_nlpi2nlp[i] < nlp->nvars);
         assert(nlp->varmap_nlp2nlpi[nlp->varmap_nlpi2nlp[i]] == i);
      }
   }

   nlp->nvars_solver -= c;
   nlp->nunflushedvardel = 0;

   /* cleanup */
   SCIPsetFreeBufferArray(set, &colset);

   return SCIP_OKAY;
}

/** adds nonlinear rows to NLPI problem that have been added to NLP before
 *
 * assumes that there are no pending variable additions or deletions (nlpFlushVarDeletions() and nlpFlushVarAdditions() should be called first)
 */
static
SCIP_RETCODE nlpFlushNlRowAdditions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   int c, i;
   SCIP_NLROW* nlrow;
   SCIP_Real*  lhss;
   SCIP_Real*  rhss;
   int*        nlinvars;
   int**       linidxs;
   SCIP_Real** lincoefs;
   SCIP_EXPR** exprs;
   const char** names;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushednlrowadd >= 0);
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->nunflushedvardel == 0);
   assert(!nlp->indiving);

   if( nlp->nunflushednlrowadd == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending additions of variables */
      for( i = 0; i < nlp->nnlrows; ++i )
         assert(nlp->nlrows[i]->nlpiindex >= 0);
#endif
      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( nlpEnsureNlRowsSolverSize(nlp, blkmem, set, nlp->nnlrows_solver + nlp->nunflushednlrowadd) );

   SCIP_CALL( SCIPsetAllocBufferArray(set, &lhss,     nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rhss,     nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &nlinvars, nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &linidxs,  nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lincoefs, nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &exprs,    nlp->nunflushednlrowadd) );
#if ADDNAMESTONLPI
   SCIP_CALL( SCIPsetAllocBufferArray(set, &names,    nlp->nunflushednlrowadd) );
#else
   names = NULL;
#endif

   c = 0;
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      nlrow = nlp->nlrows[i];
      assert(nlrow != NULL);

      /* skip nonlinear rows already in NLPI problem */
      if( nlrow->nlpiindex >= 0 )
         continue;
      assert(c < nlp->nunflushednlrowadd);

      /* get indices in NLPI */
      SCIP_CALL( nlpSetupNlpiIndices(nlp, set, nlrow, &linidxs[c]) );
      assert(linidxs[c] != NULL || nlrow->nlinvars == 0);

      nlp->nlrowmap_nlpi2nlp[nlp->nnlrows_solver+c] = i;
      nlrow->nlpiindex = nlp->nnlrows_solver+c;

      lhss[c] = nlrow->lhs;
      rhss[c] = nlrow->rhs;
      if( nlrow->constant != 0.0 )
      {
         if( !SCIPsetIsInfinity(set, -nlrow->lhs) )
            lhss[c] -= nlrow->constant;
         if( !SCIPsetIsInfinity(set,  nlrow->rhs) )
            rhss[c] -= nlrow->constant;
      }
      if( rhss[c] < lhss[c] )
      {
         assert(SCIPsetIsEQ(set, lhss[c], rhss[c]));
         rhss[c] = lhss[c];
      }

      nlinvars[c] = nlrow->nlinvars;
      lincoefs[c] = nlrow->lincoefs;

      if( nlrow->expr != NULL )
      {
         /* create copy of expr that uses varidx expressions corresponding to variables indices in NLPI */
         SCIP_CALL( SCIPexprCopy(set, stat, blkmem, set, stat, blkmem, nlrow->expr, &exprs[c], mapvar2varidx, (void*)nlp, NULL, NULL) );
      }
      else
         exprs[c] = NULL;

#if ADDNAMESTONLPI
      names[c]      = nlrow->name;
#endif

      ++c;

#ifdef NDEBUG
      /* have c vars to add already, there can be no more */
      if( c == nlp->nunflushednlrowadd )
         break;
#endif
   }
   assert(c == nlp->nunflushednlrowadd);

   nlp->nnlrows_solver += c;

   SCIP_CALL( SCIPnlpiAddConstraints(set, nlp->solver, nlp->problem, c, lhss, rhss,
         nlinvars, linidxs, lincoefs,
         exprs,
         names) );

   for( c = nlp->nunflushednlrowadd - 1; c >= 0 ; --c )
   {
      if( linidxs[c] != NULL )
         SCIPsetFreeBufferArray(set, &linidxs[c]);
      if( exprs[c] != NULL )
      {
         SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &exprs[c]) );
      }
   }

#if ADDNAMESTONLPI
   SCIPsetFreeBufferArray(set, &names);
#endif
   SCIPsetFreeBufferArray(set, &exprs);
   SCIPsetFreeBufferArray(set, &lincoefs);
   SCIPsetFreeBufferArray(set, &linidxs);
   SCIPsetFreeBufferArray(set, &nlinvars);
   SCIPsetFreeBufferArray(set, &rhss);
   SCIPsetFreeBufferArray(set, &lhss);

   nlp->nunflushednlrowadd = 0;

   return SCIP_OKAY;
}


/** adds variables to NLPI problem that have been added to NLP before
 *
 * may set nlp->objflushed to FALSE if a variable with nonzero objective coefficient is added to the NLPI problem
 */
static
SCIP_RETCODE nlpFlushVarAdditions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i, c;
   SCIP_Real*  lbs;
   SCIP_Real*  ubs;
   const char** names;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushedvaradd >= 0);
   assert(!nlp->indiving);

   if( nlp->nunflushedvaradd == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending additions of variables */
      for( i = 0; i < nlp->nvars; ++i )
         assert(nlp->varmap_nlp2nlpi[i] >= 0);
#endif
      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( nlpEnsureVarsSolverSize(nlp, blkmem, set, nlp->nvars_solver + nlp->nunflushedvaradd) );

   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbs,   nlp->nunflushedvaradd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubs,   nlp->nunflushedvaradd) );
#if ADDNAMESTONLPI
   SCIP_CALL( SCIPsetAllocBufferArray(set, &names, nlp->nunflushedvaradd) );
#else
   names = NULL;
#endif

   c = 0;
   for( i = 0; i < nlp->nvars; ++i )
   {
      /* skip variables already in NLPI problem */
      if( nlp->varmap_nlp2nlpi[i] >= 0 )
         continue;
      assert(c < nlp->nunflushedvaradd);

      nlp->varmap_nlpi2nlp[nlp->nvars_solver+c] = i;
      nlp->varmap_nlp2nlpi[i] = nlp->nvars_solver+c;
      lbs[c]   = SCIPvarGetLbLocal(nlp->vars[i]);
      ubs[c]   = SCIPvarGetUbLocal(nlp->vars[i]);
#if ADDNAMESTONLPI
      names[c] = SCIPvarGetName(nlp->vars[i]);
#endif
      ++c;

      /* if the new variable has a nonzero objective coefficient, then the objective need to be updated */
      if( !SCIPsetIsZero(set, SCIPvarGetObj(nlp->vars[i])) )
         nlp->objflushed = FALSE;

#ifdef NDEBUG
      /* have c vars to add already, there can be no more */
      if( c == nlp->nunflushedvaradd )
         break;
#endif
   }
   assert(c == nlp->nunflushedvaradd);

   nlp->nvars_solver += c;

   SCIP_CALL( SCIPnlpiAddVars(set, nlp->solver, nlp->problem, c, lbs, ubs, names) );

#if ADDNAMESTONLPI
   SCIPsetFreeBufferArray(set, &names);
#endif
   SCIPsetFreeBufferArray(set, &ubs);
   SCIPsetFreeBufferArray(set, &lbs);

   nlp->nunflushedvaradd = 0;

   return SCIP_OKAY;
}

/** updates the objective in the NLPI problem, if necessary
 *
 * assumes that there are no unflushed variable additions or deletions (nlpFlushVarDeletions() and nlpFlushVarAdditions() should be called first)
 */
static
SCIP_RETCODE nlpFlushObjective(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int*       linindices;
   SCIP_Real* lincoefs;
   SCIP_Real  coef;
   int        i;
   int        nz;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->nunflushedvardel == 0);
   assert(!nlp->indiving);

   if( nlp->objflushed )
      return SCIP_OKAY;

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* assemble coefficients */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &linindices, nlp->nvars_solver) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lincoefs,   nlp->nvars_solver) );

   nz = 0;
   for( i = 0; i < nlp->nvars_solver; ++i )
   {
      assert(nlp->varmap_nlpi2nlp[i] >= 0); /* there should be no variable deletions pending */

      coef = SCIPvarGetObj(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
      if( SCIPsetIsZero(set, coef) )
         continue;

      linindices[nz] = i;
      lincoefs[nz]   = coef;
      ++nz;
   }

   SCIP_CALL( SCIPnlpiSetObjective(set, nlp->solver, nlp->problem,
         nz, linindices, lincoefs,
         NULL,
         0.0) );

   SCIPsetFreeBufferArray(set, &lincoefs);
   SCIPsetFreeBufferArray(set, &linindices);

   nlp->objflushed = TRUE;

   return SCIP_OKAY;
}

/** solves the NLP (or diving NLP), assuming it has been flushed already */
static
SCIP_RETCODE nlpSolve(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLPPARAM*        nlpparam            /**< NLP solve parameters */
   )
{
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   if( nlp->solver == NULL )
   {
      SCIPmessagePrintWarning(messagehdlr, "Attempted to solve NLP, but no solver available.\n");

      nlp->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlp->termstat = SCIP_NLPTERMSTAT_OTHER;

      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* set initial guess, if available and warmstart hasn't been enabled
    * when using the NLP, passing a dual solution with the initguess is not available at the moment (TODO),
    * so a warmstart has to start from the last solution stored in the NLPI
    */
   if( nlp->haveinitguess && !nlpparam->warmstart )
   {
      /* @todo should we not set it if we had set it already? (initguessflushed...) */
      SCIP_Real* initialguess_solver;
      int nlpidx;

      assert(nlp->initialguess != NULL);

      SCIP_CALL( SCIPsetAllocBufferArray(set, &initialguess_solver, nlp->nvars_solver) );

      for( i = 0; i < nlp->nvars_solver; ++i )
      {
         nlpidx = nlp->varmap_nlpi2nlp[i];
         assert(nlpidx >= 0);
         assert(nlpidx < nlp->nvars);

         initialguess_solver[i] = nlp->initialguess[nlpidx];
      }
      SCIP_CALL( SCIPnlpiSetInitialGuess(set, nlp->solver, nlp->problem, initialguess_solver, NULL, NULL, NULL) );

      SCIPsetFreeBufferArray(set, &initialguess_solver);
   }

   /* let NLP solver do his work */
   SCIPclockStart(stat->nlpsoltime, set);

   SCIP_CALL( SCIPnlpiSolve(set, stat, nlp->solver, nlp->problem, nlpparam) );

   SCIPclockStop(stat->nlpsoltime, set);
   ++stat->nnlps;

   nlp->termstat = SCIPnlpiGetTermstat(set, nlp->solver, nlp->problem);
   nlp->solstat  = SCIPnlpiGetSolstat(set, nlp->solver, nlp->problem);
   switch( nlp->solstat )
   {
   case SCIP_NLPSOLSTAT_GLOBOPT:
   case SCIP_NLPSOLSTAT_LOCOPT:
   case SCIP_NLPSOLSTAT_FEASIBLE:
   case SCIP_NLPSOLSTAT_LOCINFEASIBLE:
   {
      SCIP_Real* primalvals;
      SCIP_Real* nlrowdualvals;
      SCIP_Real* varlbdualvals;
      SCIP_Real* varubdualvals;

      primalvals    = NULL;
      nlrowdualvals = NULL;
      varlbdualvals = NULL;
      varubdualvals = NULL;

      /* get NLP solution */
      SCIP_CALL( SCIPnlpiGetSolution(set, nlp->solver, nlp->problem, &primalvals, &nlrowdualvals, &varlbdualvals, &varubdualvals, NULL) );
      assert(primalvals != NULL || nlp->nvars == 0);
      assert((varlbdualvals != NULL) == (varubdualvals != NULL)); /* if there are duals for one bound, then there should also be duals for the other bound */

      /* store solution primal values in variable and evaluate objective function */
      if( nlp->indiving && nlp->divingobj != NULL )
      {
         for( i = 0; i < nlp->nvars; ++i )
         {
            SCIP_CALL( SCIPvarSetNLPSol(nlp->vars[i], set, primalvals[nlp->varmap_nlp2nlpi[i]]) );  /*lint !e613 */
         }

         /* evaluate modified diving objective */
         SCIP_CALL( SCIPnlrowGetNLPActivity(nlp->divingobj, blkmem, set, stat, primal, tree, nlp, &nlp->primalsolobjval) );
      }
      else
      {
         /* evaluate SCIP objective function */
         nlp->primalsolobjval = 0.0;
         for( i = 0; i < nlp->nvars; ++i )
         {
            SCIP_Real solval = primalvals[nlp->varmap_nlp2nlpi[i]];  /*lint !e613 */

            /* do a quick assert that variable bounds are satisfied, if feasibility is claimed */
            assert(SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(nlp->vars[i])) ||
               SCIPsetIsFeasGE(set, solval, SCIPvarGetLbLocal(nlp->vars[i])) || nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE);
            assert(SCIPsetIsInfinity(set, SCIPvarGetUbLocal(nlp->vars[i])) ||
               SCIPsetIsFeasLE(set, solval, SCIPvarGetUbLocal(nlp->vars[i])) || nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE);

            SCIP_CALL( SCIPvarSetNLPSol(nlp->vars[i], set, solval) );  /*lint !e613 */
            nlp->primalsolobjval += SCIPvarGetObj(nlp->vars[i]) * solval;  /*lint !e613 */
         }
      }

      /* store solution dual values in nlrows and variables */
      for( i = 0; i < nlp->nnlrows; ++i )
      {
         assert(nlp->nlrows[i]->nlpiindex >= 0); /* NLP was flushed before solve, so all nlrows should be in there */

         nlp->nlrows[i]->dualsol = nlrowdualvals != NULL ? nlrowdualvals[nlp->nlrows[i]->nlpiindex] : 0.0;

         /* SCIPsetDebugMsg(set, "dual of nlrow <%s> = %g\n", nlp->nlrows[i]->name, nlp->nlrows[i]->dualsol); */
      }
      assert(nlp->varlbdualvals != NULL || nlp->nvars == 0);
      assert(nlp->varubdualvals != NULL || nlp->nvars == 0);
      if( varlbdualvals != NULL )
      {
         for( i = 0; i < nlp->nvars; ++i )
         {
            assert(nlp->varmap_nlp2nlpi[i] >= 0); /* NLP was flushed before solve, so all vars should be in there */

            nlp->varlbdualvals[i] = varlbdualvals[nlp->varmap_nlp2nlpi[i]];
            nlp->varubdualvals[i] = varubdualvals[nlp->varmap_nlp2nlpi[i]];

            /* SCIPsetDebugMsg(set, "duals of var <%s> = %g %g\n", SCIPvarGetName(nlp->vars[i]), nlp->varlbdualvals[i], nlp->varubdualvals[i]); */
         }
      }
      else if( nlp->nvars > 0 )
      {
         BMSclearMemoryArray(nlp->varlbdualvals, nlp->nvars);
         BMSclearMemoryArray(nlp->varubdualvals, nlp->nvars);
      }

      break;
   }
   default:
      nlp->primalsolobjval = SCIP_INVALID;
      break;
   } /*lint !e788*/

   return SCIP_OKAY;
}

/** assembles list of fractional variables in last NLP solution */
static
SCIP_RETCODE nlpCalcFracVars(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(nlp->validfracvars <= stat->nnlps);
   assert(SCIPnlpHasSolution(nlp));

   SCIPsetDebugMsg(set, "calculating NLP fractional variables: validfracvars=%" SCIP_LONGINT_FORMAT ", nnlps=%" SCIP_LONGINT_FORMAT "\n", nlp->validfracvars, stat->nnlps);

   if( nlp->solstat > SCIP_NLPSOLSTAT_LOCINFEASIBLE )
   {
      nlp->nfracvars     = 0;
      nlp->npriofracvars = 0;
      nlp->validfracvars = stat->nnlps;

      SCIPsetDebugMsg(set, "NLP globally infeasible, unbounded, or worse -> no solution values -> no fractional variables\n");
      return SCIP_OKAY;
   }

   /* check, if the current NLP fractional variables array is invalid */
   if( nlp->validfracvars < stat->nnlps )
   {
      SCIP_VAR* var;
      SCIP_Real primsol;
      SCIP_Real frac;
      int branchpriority;
      int insertpos;
      int maxpriority;
      int i;

      SCIPsetDebugMsg(set, " -> recalculating NLP fractional variables\n");

      if( nlp->fracvarssize == 0 )
      {
         assert(nlp->fracvars     == NULL);
         assert(nlp->fracvarssol  == NULL);
         assert(nlp->fracvarsfrac == NULL);
         nlp->fracvarssize = 5;
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlp->fracvars,     nlp->fracvarssize) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlp->fracvarssol,  nlp->fracvarssize) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlp->fracvarsfrac, nlp->fracvarssize) );
      }

      maxpriority = INT_MIN;
      nlp->nfracvars = 0;
      nlp->npriofracvars = 0;
      for( i = 0; i < nlp->nvars; ++i )
      {
         var = nlp->vars[i];
         assert(var != NULL);

         primsol = SCIPvarGetNLPSol(var);
         assert(primsol < SCIP_INVALID);

         /* consider only binary and integer variables */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY && SCIPvarGetType(var) != SCIP_VARTYPE_INTEGER )
            continue;

         /* ignore fixed variables (due to numerics, it is possible, that the NLP solution of a fixed integer variable
          * (with large fixed value) is fractional in terms of absolute feasibility measure)
          */
         if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
            continue;

         /* check, if the LP solution value is fractional */
         frac = SCIPsetFeasFrac(set, primsol);

         /* The fractionality should not be smaller than -feastol, however, if the primsol is large enough
          * and close to an integer, fixed precision floating point arithmetic might give us values slightly
          * smaller than -feastol. Originally, the "frac >= -feastol"-check was within SCIPsetIsFeasFracIntegral(),
          * however, we relaxed it to "frac >= -2*feastol" and have the stricter check here for small-enough primsols.
          */
         assert(SCIPsetIsGE(set, frac, -SCIPsetFeastol(set)) || (primsol > 1e14 * SCIPsetFeastol(set)));

         if( SCIPsetIsFeasFracIntegral(set, frac) )
            continue;

         /* ensure enough space in fracvars arrays */
         if( nlp->fracvarssize <= nlp->nfracvars )
         {
            int newsize;

            newsize = SCIPsetCalcMemGrowSize(set, nlp->nfracvars + 1);
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->fracvars,     nlp->fracvarssize, newsize) );
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->fracvarssol,  nlp->fracvarssize, newsize) );
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->fracvarsfrac, nlp->fracvarssize, newsize) );
            nlp->fracvarssize = newsize;
         }
         assert(nlp->nfracvars < nlp->fracvarssize);
         assert(nlp->fracvars     != NULL);
         assert(nlp->fracvarssol  != NULL);
         assert(nlp->fracvarsfrac != NULL);

         /* insert candidate in candidate list */
         branchpriority = SCIPvarGetBranchPriority(var);
         insertpos = nlp->nfracvars;
         nlp->nfracvars++;
         if( branchpriority > maxpriority )
         {
            /* candidate has higher priority than the current maximum:
             * move it to the front and declare it to be the single best candidate
             */
            if( insertpos != 0 )
            {
               nlp->fracvars[insertpos]     = nlp->fracvars[0];
               nlp->fracvarssol[insertpos]  = nlp->fracvarssol[0];
               nlp->fracvarsfrac[insertpos] = nlp->fracvarsfrac[0];
               insertpos = 0;
            }
            nlp->npriofracvars = 1;
            maxpriority = branchpriority;
         }
         else if( branchpriority == maxpriority )
         {
            /* candidate has equal priority as the current maximum:
             * move away the first non-maximal priority candidate, move the current candidate to the correct
             * slot (binaries first) and increase the number of maximal priority candidates
             */
            if( insertpos != nlp->npriofracvars )
            {
               nlp->fracvars[insertpos]     = nlp->fracvars[nlp->npriofracvars];
               nlp->fracvarssol[insertpos]  = nlp->fracvarssol[nlp->npriofracvars];
               nlp->fracvarsfrac[insertpos] = nlp->fracvarsfrac[nlp->npriofracvars];
               insertpos = nlp->npriofracvars;
            }
            ++nlp->npriofracvars;
         }
         nlp->fracvars[insertpos]     = var;
         nlp->fracvarssol[insertpos]  = primsol;
         nlp->fracvarsfrac[insertpos] = frac;

         SCIPsetDebugMsg(set, " -> candidate %d: var=<%s>, sol=%g, frac=%g, prio=%d (max: %d) -> pos %d\n",
            nlp->nfracvars, SCIPvarGetName(var), primsol, frac, branchpriority, maxpriority, insertpos);
      }

      nlp->validfracvars = stat->nnlps;
   }
   assert(0 <= nlp->npriofracvars);
   assert(nlp->npriofracvars <= nlp->nfracvars);

   SCIPsetDebugMsg(set, " -> %d fractional variables (%d of maximal priority)\n", nlp->nfracvars, nlp->npriofracvars);

   return SCIP_OKAY;
}

/** event handling for variable events */
static
SCIP_DECL_EVENTEXEC(eventExecNlp)
{
   SCIP_EVENTTYPE etype;
   SCIP_VAR*      var;

   assert(scip      != NULL);
   assert(eventhdlr != NULL);
   assert(event     != NULL);
   assert(eventdata != NULL);

   assert((SCIP_NLP*)eventdata == scip->nlp);

   etype = SCIPeventGetType(event);
   var   = SCIPeventGetVar(event);

   if( SCIP_EVENTTYPE_VARADDED & etype )
   {
      SCIPdebugMessage("-> handling varadd event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( SCIPnlpAddVar(scip->nlp, SCIPblkmem(scip), scip->set, var) );
   }
   else if( SCIP_EVENTTYPE_VARDELETED & etype )
   {
      SCIPdebugMessage("-> handling vardel event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( SCIPnlpDelVar(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, scip->eventqueue, scip->lp, var) );
   }
   else if( SCIP_EVENTTYPE_VARFIXED & etype )
   {
      /* variable was fixed, aggregated, or multi-aggregated */
      /* TODO is this ever happening? that is, can we have changes in a variable status during solve? */
      SCIPdebugMessage("-> handling variable fixation event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( nlpRemoveFixedVar(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, scip->eventqueue, scip->lp, var) );
   }
   else if( SCIP_EVENTTYPE_BOUNDCHANGED & etype )
   {
      SCIPdebugMessage("-> handling bound changed event %" SCIP_EVENTTYPE_FORMAT ", variable <%s>\n", etype, SCIPvarGetName(var) );
      SCIP_CALL( nlpUpdateVarBounds(scip->nlp, scip->set, var, (SCIP_Bool)(SCIP_EVENTTYPE_BOUNDTIGHTENED & etype)) );
   }
   else if( SCIP_EVENTTYPE_OBJCHANGED & etype )
   {
      SCIPdebugMessage("-> handling objchg event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( nlpUpdateObjCoef(scip->set, scip->nlp, var) );
   }
   else
   {
      SCIPerrorMessage("unexpected event %" SCIP_EVENTTYPE_FORMAT " on variable <%s>\n", etype, SCIPvarGetName(var) );
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/*
 * public NLP methods
 */

/** includes event handler that is used by NLP */
SCIP_RETCODE SCIPnlpInclude(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(set->stage == SCIP_STAGE_INIT);

   /* check whether event handler is already present */
   if( SCIPsetFindEventhdlr(set, EVENTHDLR_NAME) != NULL )
   {
      SCIPerrorMessage("event handler <" EVENTHDLR_NAME "> already included.\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPeventhdlrCreate(&eventhdlr, set, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventExecNlp, NULL) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(set, eventhdlr) );

   return SCIP_OKAY;
} /*lint !e715*/

/** construct a new empty NLP */
SCIP_RETCODE SCIPnlpCreate(
   SCIP_NLP**            nlp,                /**< NLP handler, call by reference */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< problem name */
   int                   nvars_estimate      /**< an estimate on the number of variables that may be added to the NLP later */
   )
{
   assert(nlp  != NULL);
   assert(blkmem != NULL);
   assert(set  != NULL);
   assert(stat != NULL);
   assert(name != NULL);

   SCIP_ALLOC( BMSallocMemory(nlp) );

   /* select NLP solver (if any available) and setup problem */
   if( set->nnlpis > 0 )
   {
      assert(set->nlp_solver != NULL);
      if( set->nlp_solver[0] == '\0' )
      { /* take solver with highest priority */
         assert(set->nlpis != NULL);

         /* sort the NLPIs if necessary */
         if( !set->nlpissorted )
            SCIPsetSortNlpis(set);

         (*nlp)->solver = set->nlpis[0];
      }
      else
      { /* find user specified NLP solver */
         (*nlp)->solver = SCIPsetFindNlpi(set, set->nlp_solver);
         if( (*nlp)->solver == NULL )
         {
            SCIPerrorMessage("Selected NLP solver <%s> not available.\n", set->nlp_solver);
            return SCIP_PLUGINNOTFOUND;
         }
      }
      assert((*nlp)->solver != NULL);
      SCIP_CALL( SCIPnlpiCreateProblem(set, (*nlp)->solver, &(*nlp)->problem, name) );
   }
   else
   {
      /* maybe someone wanna use the NLP just to collect nonlinearities, but is not necessarily interesting on solving
       * so we allow this and just continue */
      (*nlp)->solver = NULL;
      (*nlp)->problem = NULL;
   }

   /* status */
   (*nlp)->nunflushedvaradd   = 0;
   (*nlp)->nunflushedvardel   = 0;
   (*nlp)->nunflushednlrowadd = 0;
   (*nlp)->nunflushednlrowdel = 0;
   (*nlp)->indiving   = FALSE;

   /* variables in problem and NLPI problem */
   (*nlp)->nvars = 0;
   (*nlp)->sizevars = 0;
   (*nlp)->vars = NULL;
   SCIP_CALL( SCIPhashmapCreate(&(*nlp)->varhash, blkmem, nvars_estimate) );

   (*nlp)->nvars_solver = 0;
   (*nlp)->sizevars_solver = 0;
   (*nlp)->varmap_nlp2nlpi = NULL;
   (*nlp)->varmap_nlpi2nlp = NULL;

   /* nonlinear rows in problem and NLPI problem */
   (*nlp)->nnlrows = 0;
   (*nlp)->sizenlrows = 0;
   (*nlp)->nlrows = NULL;

   (*nlp)->nnlrows_solver = 0;
   (*nlp)->sizenlrows_solver = 0;
   (*nlp)->nlrowmap_nlpi2nlp = NULL;

   /* objective function */
   (*nlp)->objflushed = TRUE;
   (*nlp)->divingobj = NULL;

   /* initial guess */
   (*nlp)->haveinitguess = FALSE;
   (*nlp)->initialguess = NULL;

   /* solution of NLP */
   (*nlp)->primalsolobjval = SCIP_INVALID;
   (*nlp)->solstat         = SCIP_NLPSOLSTAT_UNKNOWN;
   (*nlp)->termstat        = SCIP_NLPTERMSTAT_OTHER;
   (*nlp)->varlbdualvals   = NULL;
   (*nlp)->varubdualvals   = NULL;

   /* event handling: catch variable addition and deletion events */
   (*nlp)->eventhdlr = SCIPsetFindEventhdlr(set, EVENTHDLR_NAME);
   if( (*nlp)->eventhdlr == NULL )
   {
      SCIPerrorMessage("NLP eventhandler <" EVENTHDLR_NAME "> not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   SCIP_CALL( SCIPeventfilterAdd(set->scip->eventfilter, blkmem, set,
         SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_VARDELETED,
         (*nlp)->eventhdlr, (SCIP_EVENTDATA*)(*nlp), &(*nlp)->globalfilterpos) );

   /* fractional variables in last NLP solution */
   (*nlp)->fracvars     = NULL;
   (*nlp)->fracvarssol  = NULL;
   (*nlp)->fracvarsfrac = NULL;
   (*nlp)->nfracvars     = 0;
   (*nlp)->npriofracvars = 0;
   (*nlp)->fracvarssize  = 0;
   (*nlp)->validfracvars = -1;

   /* miscellaneous */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlp)->name, name, strlen(name)+1) );

   return SCIP_OKAY;
}

/** frees NLP data object */
SCIP_RETCODE SCIPnlpFree(
   SCIP_NLP**            nlp,                /**< pointer to NLP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   )
{
   assert(nlp    != NULL);
   assert(*nlp   != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   /* drop fractional variables */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->fracvars,     (*nlp)->fracvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->fracvarssol,  (*nlp)->fracvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->fracvarsfrac, (*nlp)->fracvarssize);

   /* drop global events (variable addition and deletion) */
   SCIP_CALL( SCIPeventfilterDel(set->scip->eventfilter, blkmem, set,
         SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_VARDELETED,
         (*nlp)->eventhdlr, (SCIP_EVENTDATA*)(*nlp), (*nlp)->globalfilterpos) );

   SCIP_CALL( SCIPnlpReset(*nlp, blkmem, set, stat, eventqueue, lp) );
   assert((*nlp)->nnlrows == 0);
   assert((*nlp)->nnlrows_solver == 0);
   assert((*nlp)->nvars == 0);
   assert((*nlp)->nvars_solver == 0);
   assert((*nlp)->initialguess == NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*nlp)->name, strlen((*nlp)->name)+1);

   /* free nonlinear rows arrays */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->nlrowmap_nlpi2nlp, (*nlp)->sizenlrows_solver);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->nlrows, (*nlp)->sizenlrows);

   /* free variables arrays */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varmap_nlp2nlpi, (*nlp)->sizevars);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varmap_nlpi2nlp, (*nlp)->sizevars_solver);
   SCIPhashmapFree(&(*nlp)->varhash);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->vars, (*nlp)->sizevars);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varlbdualvals, (*nlp)->sizevars);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varubdualvals, (*nlp)->sizevars);

   /* free NLPI problem */
   if( (*nlp)->problem != NULL )
   {
      SCIP_CALL( SCIPnlpiFreeProblem(set, (*nlp)->solver, &(*nlp)->problem) );
   }

   /* free NLP data structure */
   BMSfreeMemory(nlp);

   return SCIP_OKAY;
}

/** resets the NLP to the empty NLP by removing all variables and rows from NLP,
 *  releasing all rows, and flushing the changes to the NLP solver
 */
SCIP_RETCODE SCIPnlpReset(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   )
{
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   if( nlp->indiving )
   {
      SCIP_CALL( SCIPnlpEndDive(nlp, blkmem, set, stat) );
   }

   nlp->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   nlp->termstat = SCIP_NLPTERMSTAT_OTHER;

   BMSfreeBlockMemoryArrayNull(blkmem, &nlp->initialguess, nlp->sizevars);
   nlp->haveinitguess = FALSE;

   for(i = nlp->nnlrows - 1; i >= 0; --i)
   {
      SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, stat, i) );
   }

   for(i = nlp->nvars - 1; i >= 0; --i)
   {
      SCIP_CALL( nlpDelVarPos(nlp, blkmem, set, stat, eventqueue, lp, i) );
   }

   SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set, stat) );

   return SCIP_OKAY;
}

/** currently a dummy function that always returns TRUE */
SCIP_Bool SCIPnlpHasCurrentNodeNLP(
   SCIP_NLP*             nlp                 /**< NLP data */
   )
{
   assert(nlp != NULL);
   return TRUE;
} /*lint !e715*/

/** ensures, that variables array of NLP can store at least num entries */
SCIP_RETCODE SCIPnlpEnsureVarsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nvars <= nlp->sizevars);

   if( num > nlp->sizevars )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->vars,            nlp->sizevars, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varmap_nlp2nlpi, nlp->sizevars, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varlbdualvals,   nlp->sizevars, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varubdualvals,   nlp->sizevars, newsize) );
      if( nlp->initialguess != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->initialguess, nlp->sizevars, newsize) );
      }

      nlp->sizevars = newsize;
   }
   assert(num <= nlp->sizevars);

   return SCIP_OKAY;
}

/** adds a variable to the NLP and captures the variable */
SCIP_RETCODE SCIPnlpAddVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPhashmapExists(nlp->varhash, var));

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add variable during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddVars(nlp, blkmem, set, 1, &var) );

   return SCIP_OKAY;
}

/** adds a set of variables to the NLP and captures the variables */
SCIP_RETCODE SCIPnlpAddVars(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables to add */
   SCIP_VAR**            vars                /**< variables to add */
   )
{
   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(vars != NULL || nvars == 0);

   if( nlp->indiving && nvars > 0)
   {
      SCIPerrorMessage("cannot add variables during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddVars(nlp, blkmem, set, nvars, vars) );

   return SCIP_OKAY;
}

/** deletes a variable from the NLP and releases the variable */
SCIP_RETCODE SCIPnlpDelVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< SCIP LP, needed to release variable */
   SCIP_VAR*             var                 /**< variable */
   )
{
   int varpos;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(var    != NULL);

   if( !SCIPhashmapExists(nlp->varhash, var) )
   {
      SCIPerrorMessage("variable <%s> not found in NLP, cannot delete\n", SCIPvarGetName(var));
      return SCIP_ERROR;
   }

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot delete variable during NLP diving\n");
      return SCIP_ERROR;
   }

   varpos = SCIPhashmapGetImageInt(nlp->varhash, var);

   SCIP_CALL( nlpDelVarPos(nlp, blkmem, set, stat, eventqueue, lp, varpos) );

   return SCIP_OKAY;
}

/** ensures, that nonlinear rows array of NLP can store at least num entries */
SCIP_RETCODE SCIPnlpEnsureNlRowsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nnlrows <= nlp->sizenlrows);

   if( num > nlp->sizenlrows )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->nlrows, nlp->sizenlrows, newsize) );

      nlp->sizenlrows = newsize;
   }
   assert(num <= nlp->sizenlrows);

   return SCIP_OKAY;
}

/** adds a nonlinear row to the NLP and captures it
 *
 * all variables of the row need to be present in the NLP
 */
SCIP_RETCODE SCIPnlpAddNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   assert(nlp   != NULL);
   assert(nlrow != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add row during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddNlRows(nlp, blkmem, set, stat, 1, &nlrow) );

   return SCIP_OKAY;
}

/** adds nonlinear rows to the NLP and captures them
 *
 * all variables of the row need to be present in the NLP
 */
SCIP_RETCODE SCIPnlpAddNlRows(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   nnlrows,            /**< number of rows to add */
   SCIP_NLROW**          nlrows              /**< rows to add */
   )
{
   assert(nlp    != NULL);
   assert(nlrows != NULL || nnlrows == 0);

   if( nnlrows == 0 )
      return SCIP_OKAY;

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add rows during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddNlRows(nlp, blkmem, set, stat, nnlrows, nlrows) );

   return SCIP_OKAY;
}

/** deletes a nonlinear row from the NLP
 *
 * does nothing if nonlinear row is not in NLP
 */
SCIP_RETCODE SCIPnlpDelNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlrow  != NULL);

   /* if row not in NLP, nothing to do */
   if( nlrow->nlpindex == -1 )
      return SCIP_OKAY;

   assert(nlrow->nlpindex >= 0);
   assert(nlrow->nlpindex < nlp->nnlrows);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot delete row during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, stat, nlrow->nlpindex) );

   return SCIP_OKAY;
}

/** applies all cached changes to the NLP solver */
SCIP_RETCODE SCIPnlpFlush(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot flush NLP during NLP diving\n");
      return SCIP_ERROR;
   }

   /* flush removals of nonlinear rows and variables */
   SCIP_CALL( nlpFlushNlRowDeletions(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushVarDeletions(nlp, blkmem, set) );
   assert(nlp->nunflushednlrowdel == 0);
   assert(nlp->nunflushedvardel   == 0);

   /* flush addition of variables, objective, and addition of rows */
   SCIP_CALL( nlpFlushVarAdditions(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushObjective(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushNlRowAdditions(nlp, blkmem, set, stat) );
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->objflushed == TRUE);
   assert(nlp->nunflushednlrowadd == 0);

   assert(nlp->nvars   == nlp->nvars_solver);
   assert(nlp->nnlrows == nlp->nnlrows_solver);

   return SCIP_OKAY;
}

/** solves the NLP or diving NLP */
SCIP_RETCODE SCIPnlpSolve(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NLPPARAM*        nlpparam            /**< NLP solve parameters */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   if( !nlp->indiving )
   {
      SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set, stat) );
   }

   SCIP_CALL( nlpSolve(nlp, blkmem, set, messagehdlr, stat, primal, tree, nlpparam) );

   return SCIP_OKAY;
}

/** gets objective value of current NLP */
SCIP_Real SCIPnlpGetObjval(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->primalsolobjval;
}

/** gives current pseudo objective value */
SCIP_RETCODE SCIPnlpGetPseudoObjval(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< SCIP problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< SCIP LP */
   SCIP_Real*            pseudoobjval        /**< buffer to store pseudo objective value */
   )
{
   assert(nlp != NULL);
   assert(pseudoobjval != NULL);

   if( nlp->divingobj != NULL )
   {
      assert(nlp->indiving);
      SCIP_CALL( SCIPnlrowGetPseudoActivity(nlp->divingobj, blkmem, set, stat, prob, primal, tree, lp, pseudoobjval) );
   }
   else
   {
      int i;

      *pseudoobjval = 0.0;
      for( i = 0; i < nlp->nvars; ++i )
         *pseudoobjval += SCIPvarGetObj(nlp->vars[i]) * SCIPvarGetBestBoundLocal(nlp->vars[i]);
   }

   return SCIP_OKAY;
}

/** gets fractional variables of last NLP solution along with solution values and fractionalities
 */
SCIP_RETCODE SCIPnlpGetFracVars(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR***           fracvars,           /**< pointer to store the array of NLP fractional variables, or NULL */
   SCIP_Real**           fracvarssol,        /**< pointer to store the array of NLP fractional variables solution values, or NULL */
   SCIP_Real**           fracvarsfrac,       /**< pointer to store the array of NLP fractional variables fractionalities, or NULL */
   int*                  nfracvars,          /**< pointer to store the number of NLP fractional variables , or NULL */
   int*                  npriofracvars       /**< pointer to store the number of NLP fractional variables with maximal branching priority, or NULL */
   )
{
   assert(nlp != NULL);

   SCIP_CALL( nlpCalcFracVars(nlp, blkmem, set, stat) );
   assert(nlp->fracvars     != NULL);
   assert(nlp->fracvarssol  != NULL);
   assert(nlp->fracvarsfrac != NULL);

   if( fracvars != NULL )
      *fracvars = nlp->fracvars;
   if( fracvarssol != NULL )
      *fracvarssol = nlp->fracvarssol;
   if( fracvarsfrac != NULL )
      *fracvarsfrac = nlp->fracvarsfrac;
   if( nfracvars != NULL )
      *nfracvars = nlp->nfracvars;
   if( npriofracvars != NULL )
      *npriofracvars = nlp->npriofracvars;

   return SCIP_OKAY;
}

/** removes all redundant nonlinear rows */
SCIP_RETCODE SCIPnlpRemoveRedundantNlRows(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_NLPSOLSTAT solstatus;
   SCIP_Bool isredundant;
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   if( nlp->nnlrows == 0 )
      return SCIP_OKAY;

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot remove redundant rows during NLP diving\n");
      return SCIP_ERROR;
   }

   /* removing redundant rows should not change the solution status, so we reset it at the end */
   solstatus = nlp->solstat;

   for( i = 0; i < nlp->nnlrows; ++i )
   {
      SCIP_CALL( SCIPnlrowIsRedundant(nlp->nlrows[i], blkmem, set, stat, &isredundant) );
      if( isredundant )
      {
         SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, stat, i) );
      }
   }

   nlp->solstat = solstatus;

   return SCIP_OKAY;
}

/** set initial guess (approximate primal solution) for next solve
 *
 *  array initguess must be NULL or have length at least SCIPnlpGetNVars()
 */
SCIP_RETCODE SCIPnlpSetInitialGuess(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_Real*            initguess           /**< new initial guess, or NULL to clear previous one */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   /* if user wants to let NLP solver choose start point, then invalidate current initial guess both in NLP and in NLPI */
   if( initguess == NULL )
   {
      nlp->haveinitguess = FALSE;
      SCIP_CALL( SCIPnlpiSetInitialGuess(set, nlp->solver, nlp->problem, NULL, NULL, NULL, NULL) );
      return SCIP_OKAY;
   }

   if( nlp->initialguess != NULL )
   {
      BMScopyMemoryArray(nlp->initialguess, initguess, nlp->nvars);
   }
   else
   {
      assert( nlp->sizevars >= nlp->nvars );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlp->initialguess, nlp->sizevars) );
      BMScopyMemoryArray(nlp->initialguess, initguess, nlp->nvars);
   }
   nlp->haveinitguess = TRUE;

   return SCIP_OKAY;
}

/** writes NLP to a file */
SCIP_RETCODE SCIPnlpWrite(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           fname               /**< file name */
   )
{
   SCIP_RETCODE retcode = SCIP_OKAY;
   FILE* file;
   int i;

   assert(nlp != NULL);

   if( fname != NULL )
   {
      file = fopen(fname, "w");
      if( file == NULL )
      {
         SCIPerrorMessage("could not open file <%s> for writing\n", fname);
         return SCIP_FILECREATEERROR;
      }
   }
   else
      file = stdout;

   SCIPmessageFPrintInfo(messagehdlr, file, "STATISTICS\n");
   SCIPmessageFPrintInfo(messagehdlr, file, "  NLP name: %s\n", nlp->name);
   SCIPmessageFPrintInfo(messagehdlr, file, "  Variables: %d\n", nlp->nvars);
   SCIPmessageFPrintInfo(messagehdlr, file, "  Rows: %d\n", nlp->nnlrows);

   SCIPmessageFPrintInfo(messagehdlr, file, "VARIABLES\n");
   for( i = 0; i < nlp->nvars; ++i )
   {
      SCIP_CALL( SCIPvarPrint(nlp->vars[i], set, messagehdlr, file) );
   }

   SCIPmessageFPrintInfo(messagehdlr, file, "NONLINEAR ROWS\n");
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, "  ");
      SCIP_CALL_TERMINATE( retcode, SCIPnlrowPrint(nlp->nlrows[i], blkmem, set, stat, messagehdlr, file), TERMINATE );
   }

 TERMINATE:
   if( fname != NULL )
   {
      fclose(file);
   }

   return retcode;
}

/** gets array with variables of the NLP */
SCIP_VAR** SCIPnlpGetVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->vars;
}

/** gets current number of variables in NLP */
int SCIPnlpGetNVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nvars;
}

/** computes for each variables the number of NLP rows in which the variable appears in a nonlinear var */
SCIP_RETCODE SCIPnlpGetVarsNonlinearity(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   int*                  nlcount             /**< an array of length at least SCIPnlpGetNVars() to store nonlinearity counts of variables */
   )
{
   SCIP_NLROW* nlrow;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   int varidx;
   int c;

   assert(nlp != NULL);
   assert(nlcount != NULL || nlp->nvars == 0);

   BMSclearMemoryArray(nlcount, nlp->nvars);

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );

   for( c = 0; c < nlp->nnlrows; ++c )
   {
      nlrow = nlp->nlrows[c];
      assert(nlrow != NULL);

      if( nlrow->expr == NULL )
         continue;

      SCIP_CALL( SCIPexpriterInit(it, nlrow->expr, SCIP_EXPRITER_DFS, FALSE) );
      for( expr = nlrow->expr; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         if( !SCIPexprIsVar(set, expr) )
            continue;

         assert(SCIPhashmapExists(nlp->varhash, SCIPgetVarExprVar(expr)));

         varidx = SCIPhashmapGetImageInt(nlp->varhash, SCIPgetVarExprVar(expr));
         assert(varidx < nlp->nvars);
         assert(nlcount != NULL);
         ++nlcount[varidx];
      }
   }

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}


/** indicates whether there exists a row that contains a continuous variable in a nonlinear term
 *
 * @note The method may have to touch every row and nonlinear term to compute its result.
 */
SCIP_RETCODE SCIPnlpHasContinuousNonlinearity(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Bool*            result              /**< buffer to store whether continuous variable present in an expression of any row */
   )
{
   SCIP_NLROW* nlrow;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   int c;

   assert(nlp != NULL);

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   *result = FALSE;
   for( c = 0; c < nlp->nnlrows && !*result; ++c )
   {
      nlrow = nlp->nlrows[c];
      assert(nlrow != NULL);

      if( nlrow->expr == NULL )
         continue;

      for( expr = SCIPexpriterRestartDFS(it, nlrow->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         if( SCIPexprIsVar(set, expr) && SCIPvarGetType(SCIPgetVarExprVar(expr)) == SCIP_VARTYPE_CONTINUOUS )
         {
            *result = TRUE;
            break;
         }
      }
   }

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** gives dual solution values associated with lower bounds of NLP variables */
SCIP_Real* SCIPnlpGetVarsLbDualsol(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->varlbdualvals;
}

/** gives dual solution values associated with upper bounds of NLP variables */
SCIP_Real* SCIPnlpGetVarsUbDualsol(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->varubdualvals;
}

/** gets array with nonlinear rows of the NLP */
SCIP_NLROW** SCIPnlpGetNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nlrows;
}

/** gets current number of nonlinear rows in NLP */
int SCIPnlpGetNNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nnlrows;
}

/** gets the NLP solver interface */
SCIP_NLPI* SCIPnlpGetNLPI(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solver;
}

/** gets the NLP problem in the solver interface */
SCIP_NLPIPROBLEM* SCIPnlpGetNLPIProblem(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->problem;
}

/** indicates whether NLP is currently in diving mode */
SCIP_Bool SCIPnlpIsDiving(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->indiving;
}

/** gets solution status of current NLP */
SCIP_NLPSOLSTAT SCIPnlpGetSolstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solstat;
}

/** gets termination status of last NLP solve */
SCIP_NLPTERMSTAT SCIPnlpGetTermstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->termstat;
}

/** gives statistics (number of iterations, solving time, ...) of last NLP solve */
SCIP_RETCODE SCIPnlpGetStatistics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   )
{
   assert(nlp != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(statistics != NULL);

   SCIP_CALL( SCIPnlpiGetStatistics(set, nlp->solver, nlp->problem, statistics) );

   return SCIP_OKAY;
}

/** indicates whether a solution for the current NLP is available
 *
 * The solution may be optimal, feasible, or infeasible.
 * Thus, returns whether the NLP solution status is at most \ref SCIP_NLPSOLSTAT_LOCINFEASIBLE.
 */
SCIP_Bool SCIPnlpHasSolution(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solstat <= SCIP_NLPSOLSTAT_LOCINFEASIBLE;
}

/*
 * NLP diving methods
 */

/** signals start of diving */
SCIP_RETCODE SCIPnlpStartDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(nlp != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("NLP is already in diving mode\n");
      return SCIP_ERROR;
   }

   if( nlp->solver == NULL )
   {
      /* In diving mode we do not cache changes but put them directly in the NLPI problem, which does not exist if there is no solver.
       * So we forbid diving of no solver is available. */
      SCIPerrorMessage("Cannot start diving if no NLP solver is available\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set, stat) );

   nlp->indiving = TRUE;

   return SCIP_OKAY;
}

/** resets the bound and objective changes made during diving and disables diving mode */
SCIP_RETCODE SCIPnlpEndDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   int i;
   int* varidx;
   SCIP_Real* varlb;
   SCIP_Real* varub;

   assert(nlp != NULL);
   assert(set != NULL);
   assert(nlp->nvars == nlp->nvars_solver);

   if( !nlp->indiving )
   {
      SCIPerrorMessage("NLP not in diving mode, cannot end dive\n");
      return SCIP_ERROR;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* reset variable bounds in NLPI problem to their current values */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varidx, nlp->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varlb,  nlp->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varub,  nlp->nvars) );
   for( i = 0; i < nlp->nvars; ++i )
   {
      varidx[i] = i;
      varlb[i] = SCIPvarGetLbLocal(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
      varub[i] = SCIPvarGetUbLocal(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
   }

   SCIP_CALL( SCIPnlpiChgVarBounds(set, nlp->solver, nlp->problem, nlp->nvars, varidx, varlb, varub) );

   SCIPsetFreeBufferArray(set, &varidx);
   SCIPsetFreeBufferArray(set, &varlb);
   SCIPsetFreeBufferArray(set, &varub);

   /* clear diving objective, if one was used (i.e., if SCIPnlpChgVarObjDive had been called)
    * the objective in the NLPI will be reset in the next flush */
   if( nlp->divingobj != NULL )
   {
      SCIP_CALL( SCIPnlrowRelease(&nlp->divingobj, blkmem, set, stat) );
      assert(nlp->divingobj == NULL);
      assert(nlp->objflushed == FALSE);
   }

   /* we do not have a valid solution anymore */
   nlp->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   nlp->termstat = SCIP_NLPTERMSTAT_OTHER;
   nlp->primalsolobjval = SCIP_INVALID;

   nlp->indiving = FALSE;

   return SCIP_OKAY;
}

/** changes coefficient of variable in diving NLP */
SCIP_RETCODE SCIPnlpChgVarObjDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             coef                /**< new linear coefficient of variable in objective */
   )
{
   int pos;
   int objidx;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));
   assert(nlp->indiving);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* get position of variable in NLPI problem */
   pos = SCIPhashmapGetImageInt(nlp->varhash, var);
   pos = nlp->varmap_nlp2nlpi[pos];
   assert(pos >= 0);

   /* set coefficient in NLPI problem objective */
   objidx = -1;
   SCIP_CALL( SCIPnlpiChgLinearCoefs(set, nlp->solver, nlp->problem, objidx, 1, &pos, &coef) );

   /* create an nlrow that holds the diving objective, if not done yet */
   if( nlp->divingobj == NULL )
   {
      SCIP_Real* coefs;
      int        i;

      SCIP_CALL( SCIPsetAllocBufferArray(set, &coefs, nlp->nvars) );
      for( i = 0; i < nlp->nvars; ++i )
         coefs[i] = SCIPvarGetObj(nlp->vars[i]);

      SCIP_CALL( SCIPnlrowCreate(&nlp->divingobj, blkmem, set, stat, "divingobj",
            0.0, nlp->nvars, nlp->vars, coefs, NULL,
            -SCIPsetInfinity(set), SCIPsetInfinity(set),
            SCIP_EXPRCURV_LINEAR) );

      SCIPsetFreeBufferArray(set, &coefs);
   }
   assert(nlp->divingobj != NULL);

   /* modify coefficient in diving objective */
   SCIP_CALL( SCIPnlrowChgLinearCoef(nlp->divingobj, blkmem, set, stat, nlp, var, coef) );

   /* remember that we have to store objective after diving ended */
   nlp->objflushed = FALSE;

   return SCIP_OKAY;
}

/** changes bounds of variable in diving NLP */
SCIP_RETCODE SCIPnlpChgVarBoundsDive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             lb,                 /**< new lower bound of variable */
   SCIP_Real             ub                  /**< new upper bound of variable */
   )
{
   int pos;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));
   assert(nlp->indiving);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* get position of variable in NLPI problem */
   pos = SCIPhashmapGetImageInt(nlp->varhash, var);
   pos = nlp->varmap_nlp2nlpi[pos];
   assert(pos >= 0);

   /* set new bounds in NLPI */
   SCIP_CALL( SCIPnlpiChgVarBounds(set, nlp->solver, nlp->problem, 1, &pos, &lb, &ub) );

   return SCIP_OKAY;
}

/** changes bounds of a set of variables in diving NLP */
SCIP_RETCODE SCIPnlpChgVarsBoundsDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables which bounds to change */
   SCIP_VAR**            vars,               /**< variables which bounds to change */
   SCIP_Real*            lbs,                /**< new lower bounds of variables */
   SCIP_Real*            ubs                 /**< new upper bounds of variables */
   )
{
   int i;
   int* poss;

   assert(nlp  != NULL);
   assert(vars != NULL || nvars == 0);
   assert(nlp->indiving);
   assert(lbs  != NULL || nvars == 0);
   assert(ubs  != NULL || nvars == 0);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPsetAllocBufferArray(set, &poss, nvars) );

   for( i = 0; i < nvars; ++i )
   {
      assert(SCIPhashmapExists(nlp->varhash, vars[i]));  /*lint !e613*/

      /* get position of variable in NLPI problem */
      poss[i] = SCIPhashmapGetImageInt(nlp->varhash, vars[i]);   /*lint !e613*/
      poss[i] = nlp->varmap_nlp2nlpi[poss[i]];
      assert(poss[i] >= 0);
   }

   /* set new bounds in NLPI */
   SCIP_CALL( SCIPnlpiChgVarBounds(set, nlp->solver, nlp->problem, nvars, poss, lbs, ubs) );

   SCIPsetFreeBufferArray(set, &poss);

   return SCIP_OKAY;
}

/** returns whether the objective function has been changed during diving */
SCIP_Bool SCIPnlpIsDivingObjChanged(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   return nlp->divingobj != NULL;
}
