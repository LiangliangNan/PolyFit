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

/**@file   primal.c
 * @brief  methods for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/visual.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/disp.h"
#include "scip/pub_message.h"


/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that sols array can store at least num entries */
static
SCIP_RETCODE ensureSolsSize(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(primal->nsols <= primal->solssize);

   if( num > primal->solssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&primal->sols, newsize) );
      primal->solssize = newsize;
   }
   assert(num <= primal->solssize);

   return SCIP_OKAY;
}

/** ensures, that partialsols array can store at least num entries */
static
SCIP_RETCODE ensurePartialsolsSize(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(primal->npartialsols <= primal->partialsolssize);

   if( num > primal->partialsolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      newsize = MIN(newsize, set->limit_maxorigsol);

      SCIP_ALLOC( BMSreallocMemoryArray(&primal->partialsols, newsize) );
      primal->partialsolssize = newsize;
   }
   assert(num <= primal->partialsolssize);

   return SCIP_OKAY;
}

/** ensures, that existingsols array can store at least num entries */
static
SCIP_RETCODE ensureExistingsolsSize(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(primal->nexistingsols <= primal->existingsolssize);

   if( num > primal->existingsolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&primal->existingsols, newsize) );
      primal->existingsolssize = newsize;
   }
   assert(num <= primal->existingsolssize);

   return SCIP_OKAY;
}

/** creates primal data */
SCIP_RETCODE SCIPprimalCreate(
   SCIP_PRIMAL**         primal              /**< pointer to primal data */
   )
{
   assert(primal != NULL);

   SCIP_ALLOC( BMSallocMemory(primal) );
   (*primal)->sols = NULL;
   (*primal)->partialsols = NULL;
   (*primal)->existingsols = NULL;
   (*primal)->currentsol = NULL;
   (*primal)->primalray = NULL;
   (*primal)->solssize = 0;
   (*primal)->partialsolssize = 0;
   (*primal)->nsols = 0;
   (*primal)->npartialsols = 0;
   (*primal)->existingsolssize = 0;
   (*primal)->nexistingsols = 0;
   (*primal)->nsolsfound = 0;
   (*primal)->nlimsolsfound = 0;
   (*primal)->nbestsolsfound = 0;
   (*primal)->nlimbestsolsfound = 0;
   (*primal)->upperbound = SCIP_INVALID;
   (*primal)->cutoffbound = SCIP_INVALID;
   (*primal)->updateviolations = TRUE;

   return SCIP_OKAY;
}

/** frees primal data */
SCIP_RETCODE SCIPprimalFree(
   SCIP_PRIMAL**         primal,             /**< pointer to primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int s;

   assert(primal != NULL);
   assert(*primal != NULL);

   /* free temporary solution for storing current solution */
   if( (*primal)->currentsol != NULL )
   {
      SCIP_CALL( SCIPsolFree(&(*primal)->currentsol, blkmem, *primal) );
   }

   /* free solution for storing primal ray */
   if( (*primal)->primalray != NULL )
   {
      SCIP_CALL( SCIPsolFree(&(*primal)->primalray, blkmem, *primal) );
   }

   /* free feasible primal CIP solutions */
   for( s = 0; s < (*primal)->nsols; ++s )
   {
      SCIP_CALL( SCIPsolFree(&(*primal)->sols[s], blkmem, *primal) );
   }
   /* free partial CIP solutions */
   for( s = 0; s < (*primal)->npartialsols; ++s )
   {
      SCIP_CALL( SCIPsolFree(&(*primal)->partialsols[s], blkmem, *primal) );
   }
   assert((*primal)->nexistingsols == 0);

   BMSfreeMemoryArrayNull(&(*primal)->sols);
   BMSfreeMemoryArrayNull(&(*primal)->partialsols);
   BMSfreeMemoryArrayNull(&(*primal)->existingsols);
   BMSfreeMemory(primal);

   return SCIP_OKAY;
}

/** clears primal data */
SCIP_RETCODE SCIPprimalClear(
   SCIP_PRIMAL**         primal,             /**< pointer to primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int s;

   assert(primal != NULL);
   assert(*primal != NULL);

   /* free temporary solution for storing current solution */
   if( (*primal)->currentsol != NULL )
   {
      SCIP_CALL( SCIPsolFree(&(*primal)->currentsol, blkmem, *primal) );
   }

   /* free solution for storing primal ray */
   if( (*primal)->primalray != NULL )
   {
      SCIP_CALL( SCIPsolFree(&(*primal)->primalray, blkmem, *primal) );
   }

   /* free feasible primal CIP solutions */
   for( s = 0; s < (*primal)->nsols; ++s )
   {
      SCIP_CALL( SCIPsolFree(&(*primal)->sols[s], blkmem, *primal) );
   }

   (*primal)->currentsol = NULL;
   (*primal)->primalray = NULL;
   (*primal)->nsols = 0;
   (*primal)->nsolsfound = 0;
   (*primal)->nlimsolsfound = 0;
   (*primal)->nbestsolsfound = 0;
   (*primal)->nlimbestsolsfound = 0;
   (*primal)->upperbound = SCIP_INVALID;
   (*primal)->cutoffbound = SCIP_INVALID;
   (*primal)->updateviolations = TRUE;

   return SCIP_OKAY;
}

/** sets the cutoff bound in primal data and in LP solver */
static
SCIP_RETCODE primalSetCutoffbound(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   )
{
   assert(primal != NULL);
   assert(cutoffbound <= SCIPsetInfinity(set));
   assert(primal->upperbound == SCIP_INVALID || SCIPsetIsLE(set, cutoffbound, primal->upperbound)); /*lint !e777*/
   assert(!SCIPtreeInRepropagation(tree));

   SCIPsetDebugMsg(set, "changing cutoff bound from %g to %g\n", primal->cutoffbound, cutoffbound);

   primal->cutoffbound = MIN(cutoffbound, primal->upperbound); /* get rid of numerical issues */

   /* set cut off value in LP solver */
   SCIP_CALL( SCIPlpSetCutoffbound(lp, set, prob, primal->cutoffbound) );

   /* cut off leaves of the tree */
   SCIP_CALL( SCIPtreeCutoff(tree, reopt, blkmem, set, stat, eventqueue, lp, primal->cutoffbound) );

   return SCIP_OKAY;
}

/** sets the cutoff bound in primal data and in LP solver */
SCIP_RETCODE SCIPprimalSetCutoffbound(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound,        /**< new cutoff bound */
   SCIP_Bool             useforobjlimit      /**< should the cutoff bound be used to update the objective limit, if
                                              *   better? */
   )
{
   assert(primal != NULL);
   assert(cutoffbound <= SCIPsetInfinity(set));
   assert(cutoffbound <= primal->upperbound);
   assert(transprob != NULL);
   assert(origprob != NULL);

   if( cutoffbound < primal->cutoffbound )
   {
      if( useforobjlimit )
      {
         SCIP_Real objval;

         objval = SCIPprobExternObjval(transprob, origprob, set, cutoffbound);

         if( objval < SCIPprobGetObjlim(origprob, set) )
         {
            SCIPsetDebugMsg(set, "changing cutoff bound from %g to %g changes objective limit from %g to %g\n",
               primal->cutoffbound, cutoffbound, SCIPprobGetObjlim(origprob, set), objval);
            SCIPprobSetObjlim(origprob, objval);
            SCIPprobSetObjlim(transprob, objval);
         }
      }

      /* update cutoff bound */
      SCIP_CALL( primalSetCutoffbound(primal, blkmem, set, stat, transprob, eventqueue, tree, reopt, lp, cutoffbound) );
   }
   else if( cutoffbound > primal->cutoffbound )
   {
      SCIPerrorMessage("invalid increase in cutoff bound\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** sets upper bound in primal data and in LP solver */
static
SCIP_RETCODE primalSetUpperbound(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             upperbound          /**< new upper bound */
   )
{
   SCIP_Real cutoffbound;

   assert(primal != NULL);
   assert(stat != NULL);
   assert(upperbound <= SCIPsetInfinity(set));
   assert(upperbound <= primal->upperbound || stat->nnodes == 0);

   SCIPsetDebugMsg(set, "changing upper bound from %g to %g\n", primal->upperbound, upperbound);

   primal->upperbound = upperbound;

   /* if objective value is always integral, the cutoff bound can be reduced to nearly the previous integer number */
   if( SCIPprobIsObjIntegral(prob) && !SCIPsetIsInfinity(set, upperbound) )
   {
      SCIP_Real delta;

      delta = SCIPsetCutoffbounddelta(set);

      cutoffbound = SCIPsetFeasCeil(set, upperbound) - (1.0 - delta);
      cutoffbound = MIN(cutoffbound, upperbound); /* SCIPsetFeasCeil() can increase bound by almost 1.0 due to numerics
                                                   * and very large upperbound value */
   }
   else
      cutoffbound = upperbound;

   /* update cutoff bound */
   if( cutoffbound < primal->cutoffbound )
   {
      SCIP_CALL( primalSetCutoffbound(primal, blkmem, set, stat, prob, eventqueue, tree, reopt, lp, cutoffbound) );
   }

   /* update upper bound in visualization output */
   if( SCIPtreeGetCurrentDepth(tree) >= 0 )
   {
      SCIPvisualUpperbound(stat->visual, set, stat, primal->upperbound);
   }

   return SCIP_OKAY;
}

/** sets upper bound in primal data and in LP solver */
SCIP_RETCODE SCIPprimalSetUpperbound(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             upperbound          /**< new upper bound */
   )
{
   assert(primal != NULL);
   assert(upperbound <= SCIPsetInfinity(set));

   if( upperbound < primal->upperbound )
   {
      /* update primal bound */
      SCIP_CALL( primalSetUpperbound(primal, blkmem, set, stat, eventqueue, prob, tree, reopt, lp, upperbound) );
   }
   else if( upperbound > primal->upperbound )
   {
      SCIPerrorMessage("invalid increase in upper bound\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** updates upper bound and cutoff bound in primal data after a tightening of the problem's objective limit */
SCIP_RETCODE SCIPprimalUpdateObjlimit(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Real objlimit;
   SCIP_Real inf;

   assert(primal != NULL);

   /* get internal objective limit */
   objlimit = SCIPprobInternObjval(transprob, origprob, set, SCIPprobGetObjlim(origprob, set));
   inf = SCIPsetInfinity(set);
   objlimit = MIN(objlimit, inf);

   /* update the cutoff bound */
   if( objlimit < primal->cutoffbound )
   {
      SCIP_CALL( primalSetCutoffbound(primal, blkmem, set, stat, transprob, eventqueue, tree, reopt, lp, objlimit) );
   }

   /* set new upper bound (and decrease cutoff bound, if objective value is always integral) */
   if( objlimit < primal->upperbound )
   {
      SCIP_CALL( primalSetUpperbound(primal, blkmem, set, stat, eventqueue, transprob, tree, reopt, lp, objlimit) );
   }

   return SCIP_OKAY;
}

/** recalculates upper bound and cutoff bound in primal data after a change of the problem's objective offset */
SCIP_RETCODE SCIPprimalUpdateObjoffset(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_SOL* sol;
   SCIP_Real upperbound;
   SCIP_Real objval;
   SCIP_Real inf;
   int i;
   int j;

   assert(primal != NULL);
   assert(SCIPsetGetStage(set) <= SCIP_STAGE_PRESOLVED);

   /* recalculate internal objective limit */
   upperbound = SCIPprobInternObjval(transprob, origprob, set, SCIPprobGetObjlim(origprob, set));
   inf = SCIPsetInfinity(set);
   upperbound = MIN(upperbound, inf);

   /* resort current primal solutions */
   for( i = 1; i < primal->nsols; ++i )
   {
      sol = primal->sols[i];
      objval = SCIPsolGetObj(sol, set, transprob, origprob);
      for( j = i; j > 0 && objval < SCIPsolGetObj(primal->sols[j-1], set, transprob, origprob); --j )
         primal->sols[j] = primal->sols[j-1];
      primal->sols[j] = sol;
   }

   /* compare objective limit to currently best solution */
   if( primal->nsols > 0 )
   {
      SCIP_Real obj;

      assert(SCIPsolIsOriginal(primal->sols[0]));
      obj = SCIPsolGetObj(primal->sols[0], set, transprob, origprob);

      upperbound = MIN(upperbound, obj);
   }

   /* invalidate old upper bound */
   SCIP_CALL( primalSetUpperbound(primal, blkmem, set, stat, eventqueue, transprob, tree, reopt, lp, SCIPsetInfinity(set)) );

   /* reset the cutoff bound
    *
    * @note we might need to relax the bound since in presolving the objective correction of an
    *       aggregation is still in progress
    */
   SCIP_CALL( primalSetCutoffbound(primal, blkmem, set, stat, transprob, eventqueue, tree, reopt, lp, upperbound) );

   /* set new upper bound (and decrease cutoff bound, if objective value is always integral) */
   SCIP_CALL( primalSetUpperbound(primal, blkmem, set, stat, eventqueue, transprob, tree, reopt, lp, upperbound) );

   return SCIP_OKAY;
}

/** adds additional objective offset in original space to all existing solution (in original space) */
void SCIPprimalAddOrigObjoffset(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             addval              /**< additional objective offset in original space */
   )
{
   int i;

   assert(primal != NULL);
   assert(set != NULL);
   assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM);

#ifndef NDEBUG
   assert(primal->nsols == 0 || SCIPsolGetOrigin(primal->sols[0]) == SCIP_SOLORIGIN_ORIGINAL);

   /* check current order of primal solutions */
   for( i = 1; i < primal->nsols; ++i )
   {
      assert(SCIPsolGetOrigin(primal->sols[i]) == SCIP_SOLORIGIN_ORIGINAL);
      assert(SCIPsetIsLE(set, SCIPsolGetOrigObj(primal->sols[i-1]), SCIPsolGetOrigObj(primal->sols[i])));
   }
#endif

   /* check current order of primal solutions */
   for( i = 0; i < primal->nexistingsols; ++i )
   {
      assert(primal->existingsols[i] != NULL);
      SCIPsolOrigAddObjval(primal->existingsols[i], addval);
   }
}

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 */
SCIP_Bool SCIPprimalUpperboundIsSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob            /**< original problem data */
   )
{
   assert(primal != NULL);

   return (primal->nsols > 0 && SCIPsetIsEQ(set, primal->upperbound, SCIPsolGetObj(primal->sols[0], set, transprob, origprob)));
}

/** returns the primal ray thats proves unboundedness */
SCIP_SOL* SCIPprimalGetRay(
   SCIP_PRIMAL*          primal              /**< primal data */
   )
{
   assert(primal != NULL);

   return primal->primalray;
}

/** update the primal ray thats proves unboundedness */
SCIP_RETCODE SCIPprimalUpdateRay(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_SOL*             primalray,          /**< the new primal ray */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(primal != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(primalray != NULL);
   assert(blkmem != NULL);

   /* clear previously stored primal ray, if any */
   if( primal->primalray != NULL )
   {
      SCIP_CALL( SCIPsolFree(&primal->primalray, blkmem, primal) );
   }

   assert(primal->primalray == NULL);

   SCIP_CALL( SCIPsolCopy(&primal->primalray, blkmem, set, stat, primal, primalray) );

   return SCIP_OKAY;
}

/** adds primal solution to solution storage at given position */
static
SCIP_RETCODE primalAddSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL**            solptr,             /**< pointer to primal CIP solution */
   int                   insertpos,          /**< position in solution storage to add solution to */
   SCIP_Bool             replace             /**< should the solution at insertpos be replaced by the new solution? */
   )
{
   SCIP_SOL* sol;
   /* cppcheck-suppress unassignedVariable */
   SCIP_EVENT event;
   SCIP_Real obj;
   int pos;

   assert(primal != NULL);
   assert(set != NULL);
   assert(solptr != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(0 <= insertpos && insertpos < set->limit_maxsol);
   assert(tree == NULL || !SCIPtreeInRepropagation(tree));

   sol = *solptr;
   assert(sol != NULL);
   obj = SCIPsolGetObj(sol, set, transprob, origprob);

   SCIPsetDebugMsg(set, "insert primal solution %p with obj %g at position %d (replace=%u):\n",
      (void*)sol, obj, insertpos, replace);

   SCIPdebug( SCIP_CALL( SCIPsolPrint(sol, set, messagehdlr, stat, transprob, NULL, NULL, FALSE, FALSE) ) );

#if 0 /* this is not a valid debug check, but can be used to track down numerical troubles */
#ifndef NDEBUG
   /* check solution again completely
    * it fail for different reasons:
    * - in the LP solver, the feasibility tolerance is a relative measure against the row's norm
    * - in SCIP, the feasibility tolerance is a relative measure against the row's rhs/lhs
    * - the rhs/lhs of a row might drastically change during presolving when variables are fixed or (multi-)aggregated
    */
   if( !SCIPsolIsOriginal(sol) )
   {
      SCIP_Bool feasible;

      SCIP_CALL( SCIPsolCheck(sol, set, messagehdlr, blkmem, stat, transprob, TRUE, TRUE, TRUE, TRUE, &feasible) );

      if( !feasible )
      {
         SCIPerrorMessage("infeasible solution accepted:\n");
         SCIP_CALL( SCIPsolPrint(sol, set, messagehdlr, stat, origprob, transprob, NULL, FALSE, FALSE) );
      }
      assert(feasible);
   }
#endif
#endif

   /* completely fill the solution's own value array to unlink it from the LP or pseudo solution */
   SCIP_CALL( SCIPsolUnlink(sol, set, transprob) );

   /* allocate memory for solution storage */
   SCIP_CALL( ensureSolsSize(primal, set, set->limit_maxsol) );

   /* if set->limit_maxsol was decreased in the meantime, free all solutions exceeding the limit */
   for( pos = set->limit_maxsol; pos < primal->nsols; ++pos )
   {
      SCIP_CALL( SCIPsolFree(&primal->sols[pos], blkmem, primal) );
   }
   primal->nsols = MIN(primal->nsols, set->limit_maxsol);

   /* if the solution should replace an existing one, free this solution, otherwise,
    * free the last solution if the solution storage is full;
    */
   if( replace )
   {
      SCIP_CALL( SCIPsolTransform(primal->sols[insertpos], solptr, blkmem, set, primal) );
      sol = primal->sols[insertpos];
   }
   else
   {
      if( primal->nsols == set->limit_maxsol )
      {
         SCIP_CALL( SCIPsolFree(&primal->sols[set->limit_maxsol - 1], blkmem, primal) );
      }
      else
      {
         primal->nsols = primal->nsols + 1;
         assert(primal->nsols <= set->limit_maxsol);
      }

      /* move all solutions with worse objective value than the new solution */
      for( pos = primal->nsols-1; pos > insertpos; --pos )
         primal->sols[pos] = primal->sols[pos-1];

      /* insert solution at correct position */
      assert(0 <= insertpos && insertpos < primal->nsols);
      primal->sols[insertpos] = sol;
      primal->nsolsfound++;

      /* check if solution is better than objective limit */
      if( SCIPsetIsFeasLE(set, obj, SCIPprobInternObjval(transprob, origprob, set, SCIPprobGetObjlim(origprob, set))) )
         primal->nlimsolsfound++;
   }

   /* if its the first primal solution, store the relevant statistics */
   if( primal->nsolsfound == 1 )
   {
      SCIP_Real primalsolval;

      stat->nnodesbeforefirst = SCIPsolGetNodenum(sol);
      stat->nrunsbeforefirst = SCIPsolGetRunnum(sol);
      stat->firstprimalheur = SCIPsolGetHeur(sol);
      stat->firstprimaltime = SCIPsolGetTime(sol);
      stat->firstprimaldepth = SCIPsolGetDepth(sol);

      primalsolval = obj;
      stat->firstprimalbound = SCIPprobExternObjval(transprob, origprob, set, primalsolval);

      SCIPsetDebugMsg(set, "First Solution stored in problem specific statistics.\n");
      SCIPsetDebugMsg(set, "-> %" SCIP_LONGINT_FORMAT " nodes, %d runs, %.2g time, %d depth, %.15g objective\n", stat->nnodesbeforefirst, stat->nrunsbeforefirst,
         stat->firstprimaltime, stat->firstprimaldepth, stat->firstprimalbound);
   }

   SCIPsetDebugMsg(set, " -> stored at position %d of %d solutions, found %" SCIP_LONGINT_FORMAT " solutions\n",
      insertpos, primal->nsols, primal->nsolsfound);

   /* update the solution value sums in variables */
   if( !SCIPsolIsOriginal(sol) )
   {
      SCIPsolUpdateVarsum(sol, set, stat, transprob,
         (SCIP_Real)(primal->nsols - insertpos)/(SCIP_Real)(2.0*primal->nsols - 1.0));
   }

   /* change color of node in visualization output */
   SCIPvisualFoundSolution(stat->visual, set, stat, SCIPtreeGetCurrentNode(tree), insertpos == 0 ? TRUE : FALSE, sol);

   /* check, if the global upper bound has to be updated */
   if( obj < primal->cutoffbound && insertpos == 0 )
   {
      /* update the upper bound */
      SCIP_CALL( SCIPprimalSetUpperbound(primal, blkmem, set, stat, eventqueue, transprob, tree, reopt, lp, obj) );

      /* issue BESTSOLFOUND event */
      SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_BESTSOLFOUND) );
      primal->nbestsolsfound++;
      stat->bestsolnode = stat->nnodes;
   }
   else
   {
      /* issue POORSOLFOUND event */
      SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_POORSOLFOUND) );
   }
   SCIP_CALL( SCIPeventChgSol(&event, sol) );
   SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

   /* display node information line */
   if( insertpos == 0 && !replace && set->stage >= SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( SCIPdispPrintLine(set, messagehdlr, stat, NULL, TRUE, TRUE) );
   }

   /* if an original solution was added during solving, try to transfer it to the transformed space */
   if( SCIPsolIsOriginal(sol) && SCIPsetGetStage(set) == SCIP_STAGE_SOLVING && set->misc_transorigsols )
   {
      SCIP_Bool added;

      SCIP_CALL( SCIPprimalTransformSol(primal, sol, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt,
            lp, eventqueue, eventfilter, NULL, NULL, 0, &added) );

      SCIPsetDebugMsg(set, "original solution %p was successfully transferred to the transformed problem space\n",
         (void*)sol);

   }

   return SCIP_OKAY;
}

/** adds primal solution to solution storage at given position */
static
SCIP_RETCODE primalAddOrigSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< original problem data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   insertpos           /**< position in solution storage to add solution to */
   )
{
   int pos;

   assert(primal != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(sol != NULL);
   assert(0 <= insertpos && insertpos < set->limit_maxorigsol);
   assert(!set->reopt_enable);

   SCIPsetDebugMsg(set, "insert primal solution candidate %p with obj %g at position %d:\n", (void*)sol, SCIPsolGetOrigObj(sol), insertpos);

   /* allocate memory for solution storage */
   SCIP_CALL( ensureSolsSize(primal, set, set->limit_maxorigsol) );

   /* if the solution storage is full, free the last solution(s)
    * more than one solution may be freed, if set->limit_maxorigsol was decreased in the meantime
    */
   for( pos = set->limit_maxorigsol-1; pos < primal->nsols; ++pos )
   {
      SCIP_CALL( SCIPsolFree(&primal->sols[pos], blkmem, primal) );
   }

   /* insert solution at correct position */
   primal->nsols = MIN(primal->nsols+1, set->limit_maxorigsol);
   for( pos = primal->nsols-1; pos > insertpos; --pos )
      primal->sols[pos] = primal->sols[pos-1];

   assert(0 <= insertpos && insertpos < primal->nsols);
   primal->sols[insertpos] = sol;
   primal->nsolsfound++;

   /* check if solution is better than objective limit */
   if( SCIPsetIsFeasLE(set, SCIPsolGetOrigObj(sol), SCIPprobGetObjlim(prob, set)) )
      primal->nlimsolsfound++;

   SCIPsetDebugMsg(set, " -> stored at position %d of %d solutions, found %" SCIP_LONGINT_FORMAT " solutions\n",
      insertpos, primal->nsols, primal->nsolsfound);

   return SCIP_OKAY;
}

/** adds primal solution to solution storage */
static
SCIP_RETCODE primalAddOrigPartialSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< original problem data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{  /*lint --e{715}*/
   assert(primal != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(sol != NULL);

   if( primal->npartialsols >= set->limit_maxorigsol )
   {
      SCIPerrorMessage("Cannot add partial solution to storage: limit reached.\n");
      return SCIP_INVALIDCALL;
   }

   SCIPsetDebugMsg(set, "insert partial solution candidate %p:\n", (void*)sol);

   /* allocate memory for solution storage */
   SCIP_CALL( ensurePartialsolsSize(primal, set, primal->npartialsols+1) );

   primal->partialsols[primal->npartialsols] = sol;
   ++primal->npartialsols;

   return SCIP_OKAY;
}

/** uses binary search to find position in solution storage */
static
int primalSearchSolPos(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_SOL*             sol                 /**< primal solution to search position for */
   )
{
   SCIP_SOL** sols;
   SCIP_Real obj;
   SCIP_Real middleobj;
   int left;
   int right;
   int middle;

   assert(primal != NULL);

   obj = SCIPsolGetObj(sol, set, transprob, origprob);
   sols = primal->sols;

   left = -1;
   right = primal->nsols;
   while( left < right-1 )
   {
      middle = (left+right)/2;
      assert(left < middle && middle < right);
      assert(0 <= middle && middle < primal->nsols);

      middleobj = SCIPsolGetObj(sols[middle], set, transprob, origprob);

      if( obj < middleobj )
         right = middle;
      else
         left = middle;
   }
   assert(left == right-1);

   /* prefer solutions that live in the transformed space */
   if( !SCIPsolIsOriginal(sol) )
   {
      while( right > 0 && SCIPsolIsOriginal(sols[right-1])
         && SCIPsetIsEQ(set, SCIPsolGetObj(sols[right-1], set, transprob, origprob), obj) )
         --right;
   }

   return right;
}

/** uses binary search to find position in solution storage */
static
int primalSearchOrigSolPos(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sol                 /**< primal solution to search position for */
   )
{
   SCIP_Real obj;
   SCIP_Real middleobj;
   int left;
   int right;
   int middle;

   assert(primal != NULL);

   obj = SCIPsolGetOrigObj(sol);

   left = -1;
   right = primal->nsols;
   while( left < right-1 )
   {
      middle = (left+right)/2;
      assert(left < middle && middle < right);
      assert(0 <= middle && middle < primal->nsols);
      middleobj = SCIPsolGetOrigObj(primal->sols[middle]);
      if( obj < middleobj )
         right = middle;
      else
         left = middle;
   }
   assert(left == right-1);

   return right;
}

/** returns whether the given primal solution is already existent in the solution storage */
static
SCIP_Bool primalExistsSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_SOL*             sol,                /**< primal solution to search position for */
   int*                  insertpos,          /**< pointer to insertion position returned by primalSearchSolPos(); the
                                              *   position might be changed if an existing solution should be replaced */
   SCIP_Bool*            replace             /**< pointer to store whether the solution at insertpos should be replaced */
   )
{
   SCIP_Real obj;
   int i;

   assert(primal != NULL);
   assert(insertpos != NULL);
   assert(replace != NULL);
   assert(0 <= (*insertpos) && (*insertpos) <= primal->nsols);

   obj = SCIPsolGetObj(sol, set, transprob, origprob);

   assert(primal->sols != NULL || primal->nsols == 0);
   assert(primal->sols != NULL || (*insertpos) == 0);

   /* search in the better solutions */
   for( i = (*insertpos)-1; i >= 0; --i )
   {
      SCIP_Real solobj;

      solobj = SCIPsolGetObj(primal->sols[i], set, transprob, origprob);

      /* due to transferring the objective value of transformed solutions to the original space, small numerical errors might occur
       * which can lead to SCIPsetIsLE() failing in case of high absolute numbers
       */
      assert(SCIPsetIsLE(set, solobj, obj) || (REALABS(obj) > 1e+13 * SCIPsetEpsilon(set) && SCIPsetIsFeasLE(set, solobj, obj)));

      if( SCIPsetIsLT(set, solobj, obj) )
         break;

      if( SCIPsolsAreEqual(sol, primal->sols[i], set, stat, origprob, transprob) )
      {
         if( SCIPsolIsOriginal(primal->sols[i]) && !SCIPsolIsOriginal(sol) )
         {
            (*insertpos) = i;
            (*replace) = TRUE;
         }
         return TRUE;
      }
   }

   /* search in the worse solutions */
   for( i = (*insertpos); i < primal->nsols; ++i )
   {
      SCIP_Real solobj;

      solobj = SCIPsolGetObj(primal->sols[i], set, transprob, origprob);

      /* due to transferring the objective value of transformed solutions to the original space, small numerical errors might occur
       * which can lead to SCIPsetIsLE() failing in case of high absolute numbers
       */
      assert( SCIPsetIsGE(set, solobj, obj) || (REALABS(obj) > 1e+13 * SCIPsetEpsilon(set) && SCIPsetIsFeasGE(set, solobj, obj)));

      if( SCIPsetIsGT(set, solobj, obj) )
         break;

      if( SCIPsolsAreEqual(sol, primal->sols[i], set, stat, origprob, transprob) )
      {
         if( SCIPsolIsOriginal(primal->sols[i]) && !SCIPsolIsOriginal(sol) )
         {
            (*insertpos) = i;
            (*replace) = TRUE;
         }
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the given primal solution is already existent in the original solution candidate storage */
static
SCIP_Bool primalExistsOrigSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< original problem */
   SCIP_SOL*             sol,                /**< primal solution to search position for */
   int                   insertpos           /**< insertion position returned by primalSearchOrigSolPos() */
   )
{
   SCIP_Real obj;
   int i;

   assert(primal != NULL);
   assert(0 <= insertpos && insertpos <= primal->nsols);

   obj = SCIPsolGetOrigObj(sol);

   /* search in the better solutions */
   for( i = insertpos-1; i >= 0; --i )
   {
      SCIP_Real solobj;

      solobj = SCIPsolGetOrigObj(primal->sols[i]);
      assert( SCIPsetIsLE(set, solobj, obj) );

      if( SCIPsetIsLT(set, solobj, obj) )
         break;

      if( SCIPsolsAreEqual(sol, primal->sols[i], set, stat, prob, NULL) )
         return TRUE;
   }

   /* search in the worse solutions */
   for( i = insertpos; i < primal->nsols; ++i )
   {
      SCIP_Real solobj;

      solobj = SCIPsolGetOrigObj(primal->sols[i]);
      assert( SCIPsetIsGE(set, solobj, obj) );

      if( SCIPsetIsGT(set, solobj, obj) )
         break;

      if( SCIPsolsAreEqual(sol, primal->sols[i], set, stat, prob, NULL) )
         return TRUE;
   }

   return FALSE;
}

/** check if we are willing to check the solution for feasibility */
static
SCIP_Bool solOfInterest(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int*                  insertpos,          /**< pointer to store the insert position of that solution */
   SCIP_Bool*            replace             /**< pointer to store whether the solution at insertpos should be replaced
                                              *   (e.g., because it lives in the original space) */
   )
{
   SCIP_Real obj;

   obj = SCIPsolGetObj(sol, set, transprob, origprob);

   /* check if we are willing to check worse solutions; a solution is better if the objective is smaller than the
    * current cutoff bound; solutions with infinite objective value are never accepted
    */
   if( (!set->misc_improvingsols || obj < primal->cutoffbound) && !SCIPsetIsInfinity(set, obj) )
   {
      /* find insert position for the solution */
      (*insertpos) = primalSearchSolPos(primal, set, transprob, origprob, sol);
      (*replace) = FALSE;

      /* the solution should be added, if the insertpos is smaller than the maximum number of solutions to be stored
       * and it does not already exist or it does exist, but the existing solution should be replaced by the new one
       */
      if( (*insertpos) < set->limit_maxsol &&
         (!primalExistsSol(primal, set, stat, origprob, transprob, sol, insertpos, replace) || (*replace)) )
         return TRUE;
   }

   return FALSE;
}

/** check if we are willing to store the solution candidate for later checking */
static
SCIP_Bool origsolOfInterest(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int*                  insertpos           /**< pointer to store the insert position of that solution */
   )
{
   assert(SCIPsolIsOriginal(sol));

   /* find insert position for the solution */
   (*insertpos) = primalSearchOrigSolPos(primal, sol);

   if( !set->reopt_enable && (*insertpos) < set->limit_maxorigsol && !primalExistsOrigSol(primal, set, stat, origprob, sol, *insertpos) )
      return TRUE;

   return FALSE;
}

/** adds primal solution to solution storage by copying it */
SCIP_RETCODE SCIPprimalAddSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_Bool replace;
   int insertpos;

   assert(primal != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(messagehdlr != NULL);
   assert(stat != NULL);
   assert(origprob != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(eventqueue != NULL);
   assert(eventfilter != NULL);
   assert(sol != NULL);
   assert(stored != NULL);

   insertpos = -1;

   assert(!SCIPsolIsPartial(sol));

   if( solOfInterest(primal, set, stat, origprob, transprob, sol, &insertpos, &replace) )
   {
      SCIP_SOL* solcopy;
#ifdef SCIP_MORE_DEBUG
      int i;
#endif

      assert(insertpos >= 0 && insertpos < set->limit_maxsol);

      /* create a copy of the solution */
      SCIP_CALL( SCIPsolCopy(&solcopy, blkmem, set, stat, primal, sol) );

      /* insert copied solution into solution storage */
      SCIP_CALL( primalAddSol(primal, blkmem, set, messagehdlr, stat, origprob, transprob,
            tree, reopt, lp, eventqueue, eventfilter, &solcopy, insertpos, replace) );
#ifdef SCIP_MORE_DEBUG
      for( i = 0; i < primal->nsols - 1; ++i )
      {
         assert(SCIPsetIsLE(set, SCIPsolGetObj(primal->sols[i], set, transprob, origprob), SCIPsolGetObj(primal->sols[i+1], set, transprob, origprob)));
      }
#endif
      *stored = TRUE;
   }
   else
      *stored = FALSE;

   return SCIP_OKAY;
}

/** adds primal solution to solution storage, frees the solution afterwards */
SCIP_RETCODE SCIPprimalAddSolFree(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_Bool replace;
   int insertpos;

   assert(primal != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(stored != NULL);

   insertpos = -1;

   if( solOfInterest(primal, set, stat, origprob, transprob, *sol, &insertpos, &replace) )
   {
      assert(insertpos >= 0 && insertpos < set->limit_maxsol);

      /* insert solution into solution storage */
      SCIP_CALL( primalAddSol(primal, blkmem, set, messagehdlr, stat, origprob, transprob,
            tree, reopt, lp, eventqueue, eventfilter, sol, insertpos, replace) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;

      *stored = TRUE;
   }
   else
   {
      /* the solution is too bad -> free it immediately */
      SCIP_CALL( SCIPsolFree(sol, blkmem, primal) );

      *stored = FALSE;
   }
   assert(*sol == NULL);

   return SCIP_OKAY;
}

/** adds primal solution to solution candidate storage of original problem space */
SCIP_RETCODE SCIPprimalAddOrigSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< original problem data */
   SCIP_SOL*             sol,                /**< primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   int insertpos;

   assert(primal != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(sol != NULL);
   assert(SCIPsolIsOriginal(sol));
   assert(stored != NULL);

   insertpos = -1;

   if( SCIPsolIsPartial(sol) )
   {
      SCIP_SOL* solcopy;

      /* create a copy of the solution */
      SCIP_CALL( SCIPsolCopy(&solcopy, blkmem, set, stat, primal, sol) );

      SCIP_CALL( primalAddOrigPartialSol(primal, blkmem, set, prob, solcopy) );

      *stored = TRUE;
   }
   else if( origsolOfInterest(primal, set, stat, prob, sol, &insertpos) )
   {
      SCIP_SOL* solcopy;

      assert(insertpos >= 0 && insertpos < set->limit_maxorigsol);
      assert(!set->reopt_enable);

      /* create a copy of the solution */
      SCIP_CALL( SCIPsolCopy(&solcopy, blkmem, set, stat, primal, sol) );

      /* insert solution into solution storage */
      SCIP_CALL( primalAddOrigSol(primal, blkmem, set, prob, solcopy, insertpos) );

      *stored = TRUE;
   }
   else
      *stored = FALSE;

   return SCIP_OKAY;
}

/** adds primal solution to solution candidate storage of original problem space, frees the solution afterwards */
SCIP_RETCODE SCIPprimalAddOrigSolFree(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< original problem data */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   int insertpos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(SCIPsolIsOriginal(*sol));
   assert(stored != NULL);

   insertpos = -1;

   if( SCIPsolIsPartial(*sol) )
   {
      /* insert solution into solution storage */
      SCIP_CALL( primalAddOrigPartialSol(primal, blkmem, set, prob, *sol) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;

      *stored = TRUE;
   }
   else if( origsolOfInterest(primal, set, stat, prob, *sol, &insertpos) )
   {
      assert(insertpos >= 0 && insertpos < set->limit_maxorigsol);
      assert(!set->reopt_enable);

      /* insert solution into solution storage */
      SCIP_CALL( primalAddOrigSol(primal, blkmem, set, prob, *sol, insertpos) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;

      *stored = TRUE;
   }
   else
   {
      /* the solution is too bad -> free it immediately */
      SCIP_CALL( SCIPsolFree(sol, blkmem, primal) );

      *stored = FALSE;
   }
   assert(*sol == NULL);

   return SCIP_OKAY;
}

/** links temporary solution of primal data to current solution */
static
SCIP_RETCODE primalLinkCurrentSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            prob,               /**< transformed problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(primal != NULL);

   if( primal->currentsol == NULL )
   {
      SCIP_CALL( SCIPsolCreateCurrentSol(&primal->currentsol, blkmem, set, stat, prob, primal, tree, lp, heur) );
   }
   else
   {
      SCIP_CALL( SCIPsolLinkCurrentSol(primal->currentsol, set, stat, prob, tree, lp) );
      SCIPsolSetHeur(primal->currentsol, heur);
   }

   return SCIP_OKAY;
}

/** adds current LP/pseudo solution to solution storage */
SCIP_RETCODE SCIPprimalAddCurrentSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution (or NULL if it's from the tree) */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   assert(primal != NULL);

   /* link temporary solution to current solution */
   SCIP_CALL( primalLinkCurrentSol(primal, blkmem, set, stat, transprob, tree, lp, heur) );

   /* add solution to solution storage */
   SCIP_CALL( SCIPprimalAddSol(primal, blkmem, set, messagehdlr, stat, origprob, transprob,
         tree, reopt, lp, eventqueue, eventfilter, primal->currentsol, stored) );

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage by copying it */
SCIP_RETCODE SCIPprimalTrySol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   SCIP_Bool feasible;
   SCIP_Bool replace;
   int insertpos;

   assert(primal != NULL);
   assert(set != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(sol != NULL);
   assert(stored != NULL);

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || set->misc_exactsolve;

   insertpos = -1;

   if( solOfInterest(primal, set, stat, origprob, transprob, sol, &insertpos, &replace) )
   {
      /* check solution for feasibility */
      SCIP_CALL( SCIPsolCheck(sol, set, messagehdlr, blkmem, stat, transprob, printreason, completely, checkbounds,
            checkintegrality, checklprows, &feasible) );
   }
   else
      feasible = FALSE;

   if( feasible )
   {
      SCIP_SOL* solcopy;

      assert(insertpos >= 0 && insertpos < set->limit_maxsol);

      /* create a copy of the solution */
      SCIP_CALL( SCIPsolCopy(&solcopy, blkmem, set, stat, primal, sol) );

      /* insert copied solution into solution storage */
      SCIP_CALL( primalAddSol(primal, blkmem, set, messagehdlr, stat, origprob, transprob,
            tree, reopt, lp, eventqueue, eventfilter, &solcopy, insertpos, replace) );

      *stored = TRUE;
   }
   else
      *stored = FALSE;

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
SCIP_RETCODE SCIPprimalTrySolFree(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool             printreason,        /**< Should all the reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   )
{
   SCIP_Bool feasible;
   SCIP_Bool replace;
   int insertpos;

   assert(primal != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(tree != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(stored != NULL);

   *stored = FALSE;

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || set->misc_exactsolve;

   insertpos = -1;

   if( solOfInterest(primal, set, stat, origprob, transprob, *sol, &insertpos, &replace) )
   {
      /* check solution for feasibility */
      SCIP_CALL( SCIPsolCheck(*sol, set, messagehdlr, blkmem, stat, transprob, printreason, completely, checkbounds,
            checkintegrality, checklprows, &feasible) );
   }
   else
      feasible = FALSE;

   if( feasible )
   {
      assert(insertpos >= 0 && insertpos < set->limit_maxsol);

      /* insert solution into solution storage */
      SCIP_CALL( primalAddSol(primal, blkmem, set, messagehdlr, stat, origprob, transprob,
            tree, reopt, lp, eventqueue, eventfilter, sol, insertpos, replace) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;
      *stored = TRUE;
   }
   else
   {
      /* the solution is too bad or infeasible -> free it immediately */
      SCIP_CALL( SCIPsolFree(sol, blkmem, primal) );
      *stored = FALSE;
   }
   assert(*sol == NULL);

   return SCIP_OKAY;
}

/** checks current LP/pseudo solution; if feasible, adds it to storage */
SCIP_RETCODE SCIPprimalTryCurrentSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution (or NULL if it's from the tree) */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   assert(primal != NULL);

   /* link temporary solution to current solution */
   SCIP_CALL( primalLinkCurrentSol(primal, blkmem, set, stat, transprob, tree, lp, heur) );

   /* add solution to solution storage */
   SCIP_CALL( SCIPprimalTrySol(primal, blkmem, set, messagehdlr, stat, origprob, transprob,
         tree, reopt, lp, eventqueue, eventfilter, primal->currentsol,
         printreason, completely, FALSE, checkintegrality, checklprows, stored) );

   return SCIP_OKAY;
}

/** inserts solution into the global array of all existing primal solutions */
SCIP_RETCODE SCIPprimalSolCreated(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(primal != NULL);
   assert(sol != NULL);
   assert(SCIPsolGetPrimalIndex(sol) == -1);

   /* allocate memory for solution storage */
   SCIP_CALL( ensureExistingsolsSize(primal, set, primal->nexistingsols+1) );

   /* append solution */
   SCIPsolSetPrimalIndex(sol, primal->nexistingsols);
   primal->existingsols[primal->nexistingsols] = sol;
   primal->nexistingsols++;

   return SCIP_OKAY;
}

/** removes solution from the global array of all existing primal solutions */
void SCIPprimalSolFreed(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   int idx;

   assert(primal != NULL);
   assert(sol != NULL);

#ifndef NDEBUG
   for( idx = 0; idx < primal->nexistingsols; ++idx )
   {
      assert(idx == SCIPsolGetPrimalIndex(primal->existingsols[idx]));
   }
#endif

   /* remove solution */
   idx = SCIPsolGetPrimalIndex(sol);
   assert(0 <= idx && idx < primal->nexistingsols);
   assert(sol == primal->existingsols[idx]);
   if( idx < primal->nexistingsols-1 )
   {
      primal->existingsols[idx] = primal->existingsols[primal->nexistingsols-1];
      SCIPsolSetPrimalIndex(primal->existingsols[idx], idx);
   }
   primal->nexistingsols--;
}

/** updates all existing primal solutions after a change in a variable's objective value */
void SCIPprimalUpdateVarObj(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             oldobj,             /**< old objective value */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   int i;

   assert(primal != NULL);

   for( i = 0; i < primal->nexistingsols; ++i )
   {
      if( !SCIPsolIsOriginal(primal->existingsols[i]) )
         SCIPsolUpdateVarObj(primal->existingsols[i], var, oldobj, newobj);
   }
}

/** retransforms all existing solutions to original problem space
 *
 * @note as a side effect, the objective value of the solutions can change (numerical errors)
 * so we update the objective cutoff value and upper bound accordingly
 */
SCIP_RETCODE SCIPprimalRetransformSolutions(
   SCIP_PRIMAL*          primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   SCIP_Bool hasinfval;
   int i;

   assert(primal != NULL);

   for( i = 0; i < primal->nsols; ++i )
   {
      if( SCIPsolGetOrigin(primal->sols[i]) == SCIP_SOLORIGIN_ZERO )
      {
         SCIP_CALL( SCIPsolRetransform(primal->sols[i], set, stat, origprob, transprob, &hasinfval) );
      }
   }

   /* check if the global upper bound has to be updated
    * @todo we do not inform anybody about this change; if this leads to some
    * problem, a possible solution is to issue a BESTSOLFOUND event
    */
   if( primal->nsols > 0 )
   {
      SCIP_Real obj;

      obj = SCIPsolGetObj(primal->sols[0], set, transprob, origprob);
      if( obj < primal->cutoffbound )
      {
         /* update the upper bound */
         SCIP_CALL( SCIPprimalSetUpperbound(primal, blkmem, set, stat, eventqueue, transprob, tree, reopt, lp, obj) );
      }
   }

   return SCIP_OKAY;
}

/** tries to transform original solution to the transformed problem space */
SCIP_RETCODE SCIPprimalTransformSol(
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_SOL*             sol,                /**< primal solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_Real*            solvals,            /**< array for internal use to store solution values, or NULL;
                                              *   if the method is called multiple times in a row, an array with size >=
                                              *   number of active variables should be given for performance reasons */
   SCIP_Bool*            solvalset,          /**< array for internal use to store which solution values were set, or NULL;
                                              *   if the method is called multiple times in a row, an array with size >=
                                              *   number of active variables should be given for performance reasons */
   int                   solvalssize,        /**< size of solvals and solvalset arrays, should be >= number of active
                                              *   variables */
   SCIP_Bool*            added               /**< pointer to store whether the solution was added */
   )
{
   SCIP_VAR** origvars;
   SCIP_VAR** transvars;
   SCIP_VAR* var;
   SCIP_Real* localsolvals;
   SCIP_Bool* localsolvalset;
   SCIP_Real solval;
   SCIP_Real scalar;
   SCIP_Real constant;
   SCIP_Bool localarrays;
   SCIP_Bool feasible;
   int norigvars;
   int ntransvars;
   int nvarsset;
   int v;

   assert(origprob != NULL);
   assert(transprob != NULL);
   assert(SCIPsolIsOriginal(sol));
   assert(solvalssize == 0 || solvals != NULL);
   assert(solvalssize == 0 || solvalset != NULL);

   origvars = origprob->vars;
   norigvars = origprob->nvars;
   transvars = transprob->vars;
   ntransvars = transprob->nvars;
   assert(solvalssize == 0 || solvalssize >= ntransvars);

   SCIPsetDebugMsg(set, "try to transfer original solution %p with objective %g into the transformed problem space\n",
      (void*)sol, SCIPsolGetOrigObj(sol));

   /* if no solvals and solvalset arrays are given, allocate local ones, otherwise use the given ones */
   localarrays = (solvalssize == 0);
   if( localarrays )
   {
      SCIP_CALL( SCIPsetAllocBufferArray(set, &localsolvals, ntransvars) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &localsolvalset, ntransvars) );
   }
   else
   {
      localsolvals = solvals;
      localsolvalset = solvalset;
   }

   BMSclearMemoryArray(localsolvalset, ntransvars);
   feasible = TRUE;
   (*added) = FALSE;
   nvarsset = 0;

   /* for each original variable, get the corresponding active, fixed or multi-aggregated variable;
    * if it resolves to an active variable, we set its solution value or check whether an already stored solution value
    * is consistent; if it resolves to a fixed variable, we check that the fixing matches the original solution value;
    * multi-aggregated variables are skipped, because their value is defined by setting solution values for the active
    * variables, anyway
    */
   for( v = 0; v < norigvars && feasible; ++v )
   {
      var = origvars[v];

      solval = SCIPsolGetVal(sol, set, stat, var);

      /* get corresponding active, fixed, or multi-aggregated variable */
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &scalar, &constant) );
      assert(SCIPvarIsActive(var) || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

      /* check whether the fixing corresponds to the solution value of the original variable */
      if( scalar == 0.0 )
      {
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED ||
            (SCIPsetIsInfinity(set, constant) || SCIPsetIsInfinity(set, -constant)));

         if( !SCIPsetIsEQ(set, solval, constant) )
         {
            SCIPsetDebugMsg(set, "original variable <%s> (solval=%g) resolves to fixed variable <%s> (original solval=%g)\n",
               SCIPvarGetName(origvars[v]), solval, SCIPvarGetName(var), constant);
            feasible = FALSE;
         }
      }
      else if( SCIPvarIsActive(var) )
      {
         /* if we already assigned a solution value to the transformed variable, check that it corresponds to the
          * value obtained from the currently regarded original variable
          */
         if( localsolvalset[SCIPvarGetProbindex(var)] )
         {
            if( !SCIPsetIsEQ(set, solval, scalar * localsolvals[SCIPvarGetProbindex(var)] + constant) )
            {
               SCIPsetDebugMsg(set, "original variable <%s> (solval=%g) resolves to active variable <%s> with assigned solval %g (original solval=%g)\n",
                  SCIPvarGetName(origvars[v]), solval, SCIPvarGetName(var), localsolvals[SCIPvarGetProbindex(var)],
                  scalar * localsolvals[SCIPvarGetProbindex(var)] + constant);
               feasible = FALSE;
            }
         }
         /* assign solution value to the transformed variable */
         else
         {
            assert(scalar != 0.0);

            localsolvals[SCIPvarGetProbindex(var)] = (solval - constant) / scalar;
            localsolvalset[SCIPvarGetProbindex(var)] = TRUE;
            ++nvarsset;
         }
      }
#ifndef NDEBUG
      /* we do not have to handle multi-aggregated variables here, since by assigning values to all active variabes,
       * we implicitly assign values to the multi-aggregated variables, too
       */
      else
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
#endif
   }

   /* if the solution values of fixed and active variables lead to no contradiction, construct solution and try it */
   if( feasible )
   {
      SCIP_SOL* transsol;

      SCIP_CALL( SCIPsolCreate(&transsol, blkmem, set, stat, primal, tree, SCIPsolGetHeur(sol)) );

      /* set solution values for variables to which we assigned a value */
      for( v = 0; v < ntransvars; ++v )
      {
         if( localsolvalset[v] )
         {
            SCIP_CALL( SCIPsolSetVal(transsol, set, stat, tree, transvars[v], localsolvals[v]) );
         }
      }

      SCIP_CALL( SCIPprimalTrySolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob,
            tree, reopt, lp, eventqueue, eventfilter, &transsol, FALSE, FALSE, TRUE, TRUE, TRUE, added) );

      SCIPsetDebugMsg(set, "solution transferred, %d/%d active variables set (stored=%u)\n", nvarsset, ntransvars, *added);
   }
   else
      (*added) = FALSE;

   /* free local arrays, if needed */
   if( localarrays )
   {
      SCIPsetFreeBufferArray(set, &localsolvalset);
      SCIPsetFreeBufferArray(set, &localsolvals);
   }

   return SCIP_OKAY;
}


/** is the updating of violations enabled for this problem? */
SCIP_Bool SCIPprimalUpdateViolations(
   SCIP_PRIMAL*          primal              /**< problem data */
   )
{
   assert(primal != NULL);

   return primal->updateviolations;
}

/** set whether the updating of violations is turned on */
void SCIPprimalSetUpdateViolations(
   SCIP_PRIMAL*          primal,             /**< problem data */
   SCIP_Bool             updateviolations    /**< marks whether the updating of violations is turned on */
   )
{
   assert(primal != NULL);

   primal->updateviolations = updateviolations;
}
