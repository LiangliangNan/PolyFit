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

/**@file   solve.c
 * @brief  main solving loop and node processing
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/visual.h"
#include "scip/interrupt.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/mem.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/pricestore.h"
#include "scip/sepastore.h"
#include "scip/cutpool.h"
#include "scip/solve.h"
#include "scip/scip.h"
#include "scip/branch.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/disp.h"
#include "scip/heur.h"
#include "scip/nodesel.h"
#include "scip/pricer.h"
#include "scip/relax.h"
#include "scip/sepa.h"
#include "scip/prop.h"
#include "scip/pub_misc.h"
#include "scip/debug.h"
#include "scip/concurrent.h"
#include "scip/syncstore.h"


#define MAXNLPERRORS  10                /**< maximal number of LP error loops in a single node */
#define MAXNCLOCKSKIPS 64               /**< maximum number of SCIPsolveIsStopped() calls without checking the clock */
#define NINITCALLS 1000L                /**< minimum number of calls to SCIPsolveIsStopped() prior to dynamic clock skips */
#define SAFETYFACTOR 1e-2               /**< the probability that SCIP skips the clock call after the time limit has already been reached */

/** returns whether the solving process will be / was stopped before proving optimality;
 *  if the solving process was stopped, stores the reason as status in stat
 */
SCIP_Bool SCIPsolveIsStopped(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_Bool             checknodelimits     /**< should the node limits be involved in the check? */
   )
{
   assert(set != NULL);
   assert(stat != NULL);

   /* increase the number of calls to this method */
   SCIPstatIncrement(stat, set, nisstoppedcalls);

   /* in case lowerbound >= upperbound, we do not want to terminate with SCIP_STATUS_GAPLIMIT but with the ordinary
    * SCIP_STATUS_OPTIMAL/INFEASIBLE/...
    */
   if( set->stage >= SCIP_STAGE_SOLVING && SCIPsetIsLE(set, SCIPgetUpperbound(set->scip), SCIPgetLowerbound(set->scip)) )
      return TRUE;

   /* if some limit has been changed since the last call, we reset the status */
   if( set->limitchanged )
   {
      stat->status = SCIP_STATUS_UNKNOWN;
      set->limitchanged = FALSE;
   }

   if( SCIPinterrupted() || stat->userinterrupt )
   {
      stat->status = SCIP_STATUS_USERINTERRUPT;
      stat->userinterrupt = FALSE;

      /* only reset the interrupted counter if this is the main SCIP catching CTRL-C */
      if( set->misc_catchctrlc )
      {
         SCIPresetInterrupted();
      }
   }
   /* only measure the clock if a time limit is set */
   else if( set->istimelimitfinite )
   {
      /* check if we have already called this function sufficiently often for a valid estimation of its average call interval */
      if( stat->nclockskipsleft <= 0 || stat->nisstoppedcalls < NINITCALLS )
      {
         SCIP_Real currtime = SCIPclockGetTime(stat->solvingtime);

         /* use the measured time to update the average time interval between two calls to this method */
         if( set->time_rareclockcheck && stat->nisstoppedcalls >= NINITCALLS )
         {
            SCIP_Real avgisstoppedfreq;
            int nclockskips = MAXNCLOCKSKIPS;

            avgisstoppedfreq = currtime / stat->nisstoppedcalls;

            /* if we are approaching the time limit, reset the number of clock skips to 0 */
            if( (SAFETYFACTOR * (set->limit_time - currtime) / (avgisstoppedfreq + 1e-6)) < nclockskips )
               nclockskips = 0;

            stat->nclockskipsleft = nclockskips;
         }
         else
            stat->nclockskipsleft = 0;

         /* set the status if the time limit was hit */
         if( currtime >= set->limit_time )
         {
            stat->status = SCIP_STATUS_TIMELIMIT;
            return TRUE;
         }
      }
      else if( SCIPclockGetLastTime(stat->solvingtime) >= set->limit_time )
      {
         /* use information if clock has been updated more recently */
         stat->status = SCIP_STATUS_TIMELIMIT;
         return TRUE;
      }
      else
         --stat->nclockskipsleft;
   }
   if( SCIPgetConcurrentMemTotal(set->scip) >= set->limit_memory*1048576.0 - stat->externmemestim * (1.0 + SCIPgetNConcurrentSolvers(set->scip)) )
      stat->status = SCIP_STATUS_MEMLIMIT;
   else if( SCIPgetNLimSolsFound(set->scip) > 0
      && (SCIPsetIsLT(set, SCIPgetGap(set->scip), set->limit_gap)
         || SCIPsetIsLT(set, SCIPgetUpperbound(set->scip) - SCIPgetLowerbound(set->scip), set->limit_absgap)) )
      stat->status = SCIP_STATUS_GAPLIMIT;
   else if( set->limit_solutions >= 0 && set->stage >= SCIP_STAGE_PRESOLVED
      && SCIPgetNLimSolsFound(set->scip) >= set->limit_solutions )
      stat->status = SCIP_STATUS_SOLLIMIT;
   else if( set->limit_bestsol >= 0 && set->stage >= SCIP_STAGE_PRESOLVED
      && SCIPgetNBestSolsFound(set->scip) >= set->limit_bestsol )
      stat->status = SCIP_STATUS_BESTSOLLIMIT;
   else if( checknodelimits && set->limit_nodes >= 0 && stat->nnodes >= set->limit_nodes )
      stat->status = SCIP_STATUS_NODELIMIT;
   else if( checknodelimits && set->limit_totalnodes >= 0 && stat->ntotalnodes >= set->limit_totalnodes )
      stat->status = SCIP_STATUS_TOTALNODELIMIT;
   else if( checknodelimits && set->limit_stallnodes >= 0 && stat->nnodes >= stat->bestsolnode + set->limit_stallnodes )
      stat->status = SCIP_STATUS_STALLNODELIMIT;

   /* If stat->status was initialized to SCIP_STATUS_NODELIMIT or SCIP_STATUS_STALLNODELIMIT due to a previous call to SCIPsolveIsStopped(,,TRUE),
    * in the case of checknodelimits == FALSE, we do not want to report here that the solve will be stopped due to a nodelimit.
    */
   if( !checknodelimits )
      return SCIPsyncstoreSolveIsStopped(SCIPgetSyncstore(set->scip)) || (stat->status != SCIP_STATUS_UNKNOWN && stat->status != SCIP_STATUS_NODELIMIT && stat->status != SCIP_STATUS_TOTALNODELIMIT && stat->status != SCIP_STATUS_STALLNODELIMIT);
   else
      return SCIPsyncstoreSolveIsStopped(SCIPgetSyncstore(set->scip)) || (stat->status != SCIP_STATUS_UNKNOWN);
}

/** calls primal heuristics */
SCIP_RETCODE SCIPprimalHeuristics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree, or NULL if called during presolving */
   SCIP_LP*              lp,                 /**< LP data, or NULL if called during presolving or propagation */
   SCIP_NODE*            nextnode,           /**< next node that will be processed, or NULL if no more nodes left
                                              *   (only needed when calling after node heuristics) */
   SCIP_HEURTIMING       heurtiming,         /**< current point in the node solving process */
   SCIP_Bool             nodeinfeasible,     /**< was the current node already detected to be infeasible? */
   SCIP_Bool*            foundsol,           /**< pointer to store whether a solution has been found */
   SCIP_Bool*            unbounded           /**< pointer to store whether an unbounded ray was found in the LP */
   )
{  /*lint --e{715}*/

   SCIP_RESULT result;
   SCIP_Longint oldnbestsolsfound;
   SCIP_Real lowerbound;
   int ndelayedheurs;
   int depth;
   int lpstateforkdepth;
   int h;
#ifndef NDEBUG
   SCIP_Bool inprobing;
   SCIP_Bool indiving;
#endif

   assert(set != NULL);
   assert(primal != NULL);
   assert(tree != NULL || heurtiming == SCIP_HEURTIMING_BEFOREPRESOL || heurtiming == SCIP_HEURTIMING_DURINGPRESOLLOOP);
   assert(lp != NULL || heurtiming == SCIP_HEURTIMING_BEFOREPRESOL || heurtiming == SCIP_HEURTIMING_DURINGPRESOLLOOP
      || heurtiming == SCIP_HEURTIMING_AFTERPROPLOOP);
   assert(heurtiming == SCIP_HEURTIMING_BEFORENODE || heurtiming == SCIP_HEURTIMING_DURINGLPLOOP
      || heurtiming == SCIP_HEURTIMING_AFTERLPLOOP || heurtiming == SCIP_HEURTIMING_AFTERNODE
      || heurtiming == SCIP_HEURTIMING_DURINGPRICINGLOOP || heurtiming == SCIP_HEURTIMING_BEFOREPRESOL
      || heurtiming == SCIP_HEURTIMING_DURINGPRESOLLOOP || heurtiming == SCIP_HEURTIMING_AFTERPROPLOOP
      || heurtiming == (SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE));
   assert(heurtiming != SCIP_HEURTIMING_AFTERNODE || (nextnode == NULL) == (SCIPtreeGetNNodes(tree) == 0));
   assert(foundsol != NULL);

   *foundsol = FALSE;

   /* nothing to do, if no heuristics are available, or if the branch-and-bound process is finished */
   if( set->nheurs == 0 || (heurtiming == SCIP_HEURTIMING_AFTERNODE && nextnode == NULL) )
      return SCIP_OKAY;

   /* do not continue if we reached a time limit */
   if( SCIPsolveIsStopped(set, stat, FALSE) )
      return SCIP_OKAY;

   /* sort heuristics by priority, but move the delayed heuristics to the front */
   SCIPsetSortHeurs(set);

   /* specialize the AFTERNODE timing flag */
   if( (heurtiming & SCIP_HEURTIMING_AFTERNODE) == SCIP_HEURTIMING_AFTERNODE )
   {
      SCIP_Bool plunging;
      SCIP_Bool pseudonode;

      /* clear the AFTERNODE flags and replace them by the right ones */
      heurtiming &= ~SCIP_HEURTIMING_AFTERNODE;

      /* we are in plunging mode iff the next node is a sibling or a child, and no leaf */
      assert(nextnode == NULL
         || SCIPnodeGetType(nextnode) == SCIP_NODETYPE_SIBLING
         || SCIPnodeGetType(nextnode) == SCIP_NODETYPE_CHILD
         || SCIPnodeGetType(nextnode) == SCIP_NODETYPE_LEAF);
      plunging = (nextnode != NULL && SCIPnodeGetType(nextnode) != SCIP_NODETYPE_LEAF);
      pseudonode = !SCIPtreeHasFocusNodeLP(tree);
      if( plunging && SCIPtreeGetCurrentDepth(tree) > 0 ) /* call plunging heuristics also at root node */
      {
         if( !pseudonode )
            heurtiming |= SCIP_HEURTIMING_AFTERLPNODE;
         else
            heurtiming |= SCIP_HEURTIMING_AFTERPSEUDONODE;
      }
      else
      {
         if( !pseudonode )
            heurtiming |= SCIP_HEURTIMING_AFTERLPPLUNGE | SCIP_HEURTIMING_AFTERLPNODE;
         else
            heurtiming |= SCIP_HEURTIMING_AFTERPSEUDOPLUNGE | SCIP_HEURTIMING_AFTERPSEUDONODE;
      }
   }

   /* initialize the tree related data, if we are not in presolving */
   if( heurtiming == SCIP_HEURTIMING_BEFOREPRESOL || heurtiming == SCIP_HEURTIMING_DURINGPRESOLLOOP )
   {
      depth = -1;
      lpstateforkdepth = -1;

      SCIPsetDebugMsg(set, "calling primal heuristics %s presolving\n",
         heurtiming == SCIP_HEURTIMING_BEFOREPRESOL ? "before" : "during");
   }
   else
   {
      assert(tree != NULL); /* for lint */
      depth = SCIPtreeGetFocusDepth(tree);
      lpstateforkdepth = (tree->focuslpstatefork != NULL ? SCIPnodeGetDepth(tree->focuslpstatefork) : -1);

      SCIPsetDebugMsg(set, "calling primal heuristics in depth %d (timing: %u)\n", depth, heurtiming);
   }

   /* call heuristics */
   ndelayedheurs = 0;
   oldnbestsolsfound = primal->nbestsolsfound;

#ifndef NDEBUG
   /* remember old probing and diving status */
   inprobing = tree != NULL && SCIPtreeProbing(tree);
   indiving = lp != NULL && SCIPlpDiving(lp);

   /* heuristics should currently not be called in diving mode */
   assert(!indiving);
#endif

   /* collect lower bound of current node */
   if( tree !=  NULL )
   {
      assert(SCIPtreeGetFocusNode(tree) != NULL);
      lowerbound = SCIPnodeGetLowerbound(SCIPtreeGetFocusNode(tree));
   }
   else if( lp != NULL )
      lowerbound = SCIPlpGetPseudoObjval(lp, set, prob);
   else
      lowerbound = -SCIPsetInfinity(set);

   for( h = 0; h < set->nheurs; ++h )
   {
#ifndef NDEBUG
      size_t nusedbuffer = BMSgetNUsedBufferMemory(SCIPbuffer(set->scip));
#endif
      /* it might happen that a diving heuristic renders the previously solved node LP invalid
       * such that additional calls to LP heuristics will fail; better abort the loop in this case
       */
      if( lp != NULL && lp->resolvelperror)
         break;

#ifdef SCIP_DEBUG
      {
         SCIP_Bool delayed;
         if( SCIPheurShouldBeExecuted(set->heurs[h], depth, lpstateforkdepth, heurtiming, &delayed) )
         {
            SCIPsetDebugMsg(set, " -> executing heuristic <%s> with priority %d\n",
               SCIPheurGetName(set->heurs[h]), SCIPheurGetPriority(set->heurs[h]));
         }
      }
#endif

      SCIP_CALL( SCIPheurExec(set->heurs[h], set, primal, depth, lpstateforkdepth, heurtiming, nodeinfeasible,
            &ndelayedheurs, &result) );

#ifndef NDEBUG
      if( BMSgetNUsedBufferMemory(SCIPbuffer(set->scip)) > nusedbuffer )
      {
         SCIPerrorMessage("Buffer not completely freed after executing heuristic <%s>\n", SCIPheurGetName(set->heurs[h]));
         SCIPABORT();
      }
#endif

      /* if the new solution cuts off the current node due to a new primal solution (via the cutoff bound) interrupt
       * calling the remaining heuristics
       */
      if( (result == SCIP_FOUNDSOL && lowerbound > primal->cutoffbound) || SCIPsolveIsStopped(set, stat, FALSE) )
         break;

      /* check if the problem is proven to be unbounded, currently this happens only in reoptimization */
      if( result == SCIP_UNBOUNDED )
      {
         *unbounded = TRUE;
         break;
      }

      /* make sure that heuristic did not change probing or diving status */
      assert(tree == NULL || inprobing == SCIPtreeProbing(tree));
      assert(lp == NULL || indiving == SCIPlpDiving(lp));
   }
   assert(0 <= ndelayedheurs && ndelayedheurs <= set->nheurs);

   *foundsol = (primal->nbestsolsfound > oldnbestsolsfound);

   return SCIP_OKAY;
}

/** applies one round of propagation */
static
SCIP_RETCODE propagationRound(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   depth,              /**< depth level to use for propagator frequency checks */
   SCIP_Bool             fullpropagation,    /**< should all constraints be propagated (or only new ones)? */
   SCIP_Bool             onlydelayed,        /**< should only delayed propagators be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a propagator was delayed */
   SCIP_Bool*            propagain,          /**< pointer to store whether propagation should be applied again */
   SCIP_PROPTIMING       timingmask,         /**< timing mask to decide which propagators are executed */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            postpone            /**< pointer to store whether the node should be postponed */
   )
{  /*lint --e{715}*/
   SCIP_RESULT result;
   SCIP_Bool abortoncutoff;
   int i;

   assert(set != NULL);
   assert(delayed != NULL);
   assert(propagain != NULL);
   assert(cutoff != NULL);
   assert(postpone != NULL);

   *delayed = FALSE;
   *propagain = FALSE;

   /* sort propagators */
   SCIPsetSortProps(set);

   /* check if we want to abort on a cutoff; if we are not in the solving stage (e.g., in presolving), we want to abort
    * anyway
    */
   abortoncutoff = set->prop_abortoncutoff || (set->stage != SCIP_STAGE_SOLVING);

   /* call additional propagators with nonnegative priority */
   for( i = 0; i < set->nprops && !(*postpone) && (!(*cutoff) || !abortoncutoff); ++i )
   {
#ifndef NDEBUG
      size_t nusedbuffer = BMSgetNUsedBufferMemory(SCIPbuffer(set->scip));
#endif
      /* timing needs to fit */
      if( (SCIPpropGetTimingmask(set->props[i]) & timingmask) == 0 )
         continue;

      if( SCIPpropGetPriority(set->props[i]) < 0 )
         continue;

      if( onlydelayed && !SCIPpropWasDelayed(set->props[i]) )
         continue;

      SCIPsetDebugMsg(set, "calling propagator <%s>\n", SCIPpropGetName(set->props[i]));

      SCIP_CALL( SCIPpropExec(set->props[i], set, stat, depth, onlydelayed, tree->sbprobing, timingmask, &result) );

#ifndef NDEBUG
      if( BMSgetNUsedBufferMemory(SCIPbuffer(set->scip)) > nusedbuffer )
      {
         SCIPerrorMessage("Buffer not completely freed after executing propagator <%s>\n", SCIPpropGetName(set->props[i]));
         SCIPABORT();
      }
#endif

      *delayed = *delayed || (result == SCIP_DELAYED);
      *propagain = *propagain || (result == SCIP_REDUCEDDOM);

      /* beside the result pointer of the propagator we have to check if an internal cutoff was detected; this can
       * happen when a global bound change was applied which is globally valid and leads locally (for the current node
       * and others) to an infeasible problem;
       */
      *cutoff = *cutoff || (result == SCIP_CUTOFF) || (tree->cutoffdepth <= SCIPtreeGetCurrentDepth(tree));
      *postpone = (result == SCIP_DELAYNODE) && !(*cutoff);

      if( result == SCIP_CUTOFF )
      {
         SCIPsetDebugMsg(set, " -> propagator <%s> detected cutoff\n", SCIPpropGetName(set->props[i]));
      }

      /* if we work off the delayed propagators, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_REDUCEDDOM )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* propagate constraints */
   for( i = 0; i < set->nconshdlrs && !(*postpone) && (!(*cutoff) || !abortoncutoff); ++i )
   {
      /* timing needs to fit */
      if( (SCIPconshdlrGetPropTiming(set->conshdlrs[i]) & timingmask) == 0 )
         continue;

      if( onlydelayed && !SCIPconshdlrWasPropagationDelayed(set->conshdlrs[i]) )
         continue;

      SCIPsetDebugMsg(set, "calling propagation method of constraint handler <%s>\n", SCIPconshdlrGetName(set->conshdlrs[i]));

      SCIP_CALL( SCIPconshdlrPropagate(set->conshdlrs[i], blkmem, set, stat, depth, fullpropagation, onlydelayed,
            tree->sbprobing, timingmask, &result) );
      *delayed = *delayed || (result == SCIP_DELAYED);
      *propagain = *propagain || (result == SCIP_REDUCEDDOM);

      /* beside the result pointer of the propagator we have to check if an internal cutoff was detected; this can
       * happen when a global bound change was applied which is globally valid and leads locally (for the current node
       * and others) to an infeasible problem;
       */
      *cutoff = *cutoff || (result == SCIP_CUTOFF) || (tree->cutoffdepth <= SCIPtreeGetCurrentDepth(tree));
      *postpone = (result == SCIP_DELAYNODE) && !(*cutoff);

      if( result == SCIP_CUTOFF )
      {
         SCIPsetDebugMsg(set, " -> constraint handler <%s> detected cutoff in propagation\n",
            SCIPconshdlrGetName(set->conshdlrs[i]));
      }

      /* if we work off the delayed propagators, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_REDUCEDDOM )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* call additional propagators with negative priority */
   for( i = 0; i < set->nprops && !(*postpone) && (!(*cutoff) || !abortoncutoff); ++i )
   {
      /* timing needs to fit */
      if( (SCIPpropGetTimingmask(set->props[i]) & timingmask) == 0 )
         continue;

      if( SCIPpropGetPriority(set->props[i]) >= 0 )
         continue;

      if( onlydelayed && !SCIPpropWasDelayed(set->props[i]) )
         continue;

      SCIPsetDebugMsg(set, "calling propagator <%s>\n", SCIPpropGetName(set->props[i]));

      SCIP_CALL( SCIPpropExec(set->props[i], set, stat, depth, onlydelayed, tree->sbprobing, timingmask, &result) );
      *delayed = *delayed || (result == SCIP_DELAYED);
      *propagain = *propagain || (result == SCIP_REDUCEDDOM);

      /* beside the result pointer of the propagator we have to check if an internal cutoff was detected; this can
       * happen when a global bound change was applied which is globally valid and leads locally (for the current node
       * and others) to an infeasible problem;
       */
      *cutoff = *cutoff || (result == SCIP_CUTOFF) || (tree->cutoffdepth <= SCIPtreeGetCurrentDepth(tree));
      *postpone = (result == SCIP_DELAYNODE) && !(*cutoff);

      if( result == SCIP_CUTOFF )
      {
         SCIPsetDebugMsg(set, " -> propagator <%s> detected cutoff\n", SCIPpropGetName(set->props[i]));
      }

      /* if we work off the delayed propagators, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_REDUCEDDOM )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** applies domain propagation on current node */
static
SCIP_RETCODE propagateDomains(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   depth,              /**< depth level to use for propagator frequency checks */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Bool             fullpropagation,    /**< should all constraints be propagated (or only new ones)? */
   SCIP_PROPTIMING       timingmask,         /**< timing mask to decide which propagators are executed */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            postpone            /**< pointer to store whether the node should be postponed */
   )
{
   SCIP_NODE* node;
   SCIP_Bool delayed;
   SCIP_Bool propagain;
   int propround;

   assert(set != NULL);
   assert(tree != NULL);
   assert(depth >= 0);
   assert(cutoff != NULL);

   node = SCIPtreeGetCurrentNode(tree);
   assert(node != NULL && SCIPnodeIsActive(node));
   assert(SCIPnodeGetType(node) == SCIP_NODETYPE_FOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_REFOCUSNODE
      || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

   /* adjust maximal number of propagation rounds */
   if( maxproprounds == 0 )
      maxproprounds = (depth == 0 ? set->prop_maxroundsroot : set->prop_maxrounds);
   if( maxproprounds == -1 )
      maxproprounds = INT_MAX;

   SCIPsetDebugMsg(set, "domain propagation of node %p in depth %d (using depth %d, maxrounds %d, proptiming %u)\n",
      (void*)node, SCIPnodeGetDepth(node), depth, maxproprounds, timingmask);

   /* propagate as long new bound changes were found and the maximal number of propagation rounds is not exceeded */
   *cutoff = FALSE;
   *postpone = FALSE;
   propround = 0;
   propagain = TRUE;
   while( propagain && !(*cutoff) && !(*postpone) && propround < maxproprounds && !SCIPsolveIsStopped(set, stat, FALSE) )
   {
      propround++;

      /* perform the propagation round by calling the propagators and constraint handlers */
      SCIP_CALL( propagationRound(blkmem, set, stat, primal, tree, depth, fullpropagation, FALSE, &delayed, &propagain, timingmask, cutoff, postpone) );

      /* if the propagation will be terminated, call the delayed propagators */
      while( delayed && (!propagain || propround >= maxproprounds) && !(*cutoff) )
      {
         /* call the delayed propagators and constraint handlers */
         SCIP_CALL( propagationRound(blkmem, set, stat, primal, tree, depth, fullpropagation, TRUE, &delayed, &propagain, timingmask, cutoff, postpone) );
      }

      /* if a reduction was found, we want to do another full propagation round (even if the propagator only claimed
       * to have done a domain reduction without applying a domain change)
       */
      fullpropagation = TRUE;
   }

   /* mark the node to be completely propagated in the current repropagation subtree level */
   SCIPnodeMarkPropagated(node, tree);

   if( *cutoff )
   {
      SCIPsetDebugMsg(set, " --> domain propagation of node %p finished: cutoff!\n", (void*)node);
   }

   return SCIP_OKAY;
}

/** applies domain propagation on current node and flushes the conflict store afterwards */
SCIP_RETCODE SCIPpropagateDomains(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   int                   depth,              /**< depth level to use for propagator frequency checks */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_PROPTIMING       timingmask,         /**< timing mask to decide which propagators are executed */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_Bool postpone;

   /* apply domain propagation */
   SCIP_CALL( propagateDomains(blkmem, set, stat, primal, tree, depth, maxproprounds, TRUE, timingmask, cutoff, &postpone) );

   /* flush the conflict set storage */
   SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable) );

   return SCIP_OKAY;
}

/** returns whether the given variable with the old LP solution value should lead to an update of the pseudo cost entry */
static
SCIP_Bool isPseudocostUpdateValid(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             oldlpsolval,        /**< solution value of variable in old LP */
   SCIP_Bool             updateintegers,     /**< whether to update pseudo costs for integer variables */
   SCIP_Bool             updatecontinuous    /**< whether to update pseudo costs for continuous variables */
   )
{
   SCIP_Real newlpsolval;

   assert(var != NULL);

   if( !updatecontinuous && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      return FALSE;

   if( !updateintegers && SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      return FALSE;

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && set->branch_lpgainnorm != 'l' )
   {
      /* if the variable is fixed at +/- infinity or it has an unbounded domain, then the domain-based update strategies will not work */
      if( SCIPsetIsInfinity(set, REALABS(SCIPvarGetLbLocal(var))) || SCIPsetIsInfinity(set, REALABS(SCIPvarGetUbLocal(var))) )
         return FALSE;

      /* @todo if set->branch_lpgainnorm == 's', then we would need to know then domain before branching
       * since this is difficult to get, we don't check for unboundedness here and let the pscost update fail later
       * however, this makes the weights used to spread a pseudo cost update over all domain changes inaccurate
       */

      return TRUE;
   }

   /* if the old LP solution value is unknown, the pseudo cost update cannot be performed */
   if( oldlpsolval >= SCIP_INVALID )
      return FALSE;

   /* the bound change on the given variable was responsible for the gain in the dual bound, if the variable's
    * old solution value is outside the current bounds, and the new solution value is equal to the bound
    * closest to the old solution value
    */

   /* find out, which of the current bounds is violated by the old LP solution value */
   if( SCIPsetIsLT(set, oldlpsolval, SCIPvarGetLbLocal(var)) )
   {
      newlpsolval = SCIPvarGetLPSol(var);
      return SCIPsetIsEQ(set, newlpsolval, SCIPvarGetLbLocal(var));
   }
   else if( SCIPsetIsGT(set, oldlpsolval, SCIPvarGetUbLocal(var)) )
   {
      newlpsolval = SCIPvarGetLPSol(var);
      return SCIPsetIsEQ(set, newlpsolval, SCIPvarGetUbLocal(var));
   }
   else
      return FALSE;
}

/** pseudo cost flag stored in the variables to mark them for the pseudo cost update */
enum PseudocostFlag
{
   PSEUDOCOST_NONE     = 0,             /**< variable's bounds were not changed */
   PSEUDOCOST_IGNORE   = 1,             /**< bound changes on variable should be ignored for pseudo cost updates */
   PSEUDOCOST_UPDATE   = 2              /**< pseudo cost value of variable should be updated */
};
typedef enum PseudocostFlag PSEUDOCOSTFLAG;

/** updates the variable's pseudo cost values after the node's initial LP was solved */
static
SCIP_RETCODE updatePseudocost(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             updateintegers,     /**< whether to update pseudo costs for integer variables */
   SCIP_Bool             updatecontinuous    /**< whether to update pseudo costs for continuous variables */
   )
{
   SCIP_NODE* focusnode;
   int actdepth;

   assert(lp != NULL);
   assert(tree != NULL);
   assert(tree->path != NULL);

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(SCIPnodeIsActive(focusnode));
   assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   actdepth = SCIPnodeGetDepth(focusnode);
   assert(tree->path[actdepth] == focusnode);

   if( (updateintegers || updatecontinuous) && lp->solved && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && tree->focuslpstatefork != NULL )
   {
      SCIP_BOUNDCHG** updates;
      SCIP_NODE* node;
      SCIP_VAR* var;
      SCIP_Real weight;
      SCIP_Real lpgain;
      int nupdates;
      int nvalidupdates;
      int d;
      int i;

      assert(SCIPnodeIsActive(tree->focuslpstatefork));
      assert(tree->path[tree->focuslpstatefork->depth] == tree->focuslpstatefork);

      /* get a buffer for the collected bound changes; start with a size twice as large as the number of nodes between
       * current node and LP fork
       */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &updates, (int)(2*(actdepth - tree->focuslpstatefork->depth))) );
      nupdates = 0;
      nvalidupdates = 0;

      /* search the nodes from LP fork down to current node for bound changes in between; move in this direction,
       * because the bound changes closer to the LP fork are more likely to have a valid LP solution information
       * attached; collect the bound changes for pseudo cost value updates and mark the corresponding variables such
       * that they are not updated twice in case of more than one bound change on the same variable
       */
      for( d = tree->focuslpstatefork->depth+1; d <= actdepth; ++d )
      {
         node = tree->path[d];

         if( node->domchg != NULL )
         {
            SCIP_BOUNDCHG* boundchgs;
            int nboundchgs;

            boundchgs = node->domchg->domchgbound.boundchgs;
            nboundchgs = node->domchg->domchgbound.nboundchgs;
            for( i = 0; i < nboundchgs; ++i )
            {
               var = boundchgs[i].var;
               assert(var != NULL);

               /* we even collect redundant bound changes, since they were not redundant in the LP branching decision
                * and therefore should be regarded in the pseudocost updates
                *
                * however, if the variable is continuous and we normalize the pseudo costs by the domain reduction,
                * then getting the variable bound before the branching is not possible by looking at the variables branching information (since redundant branchings are not applied)
                * thus, in this case we ignore the boundchange
                */
               if( (SCIP_BOUNDCHGTYPE)boundchgs[i].boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING &&
                   (PSEUDOCOSTFLAG)var->pseudocostflag == PSEUDOCOST_NONE
                 )
               {
                  /* remember the bound change and mark the variable */
                  SCIP_CALL( SCIPsetReallocBufferArray(set, &updates, nupdates+1) );
                  updates[nupdates] = &boundchgs[i];
                  nupdates++;

                  /* check, if the bound change would lead to a valid pseudo cost update
                   * and see comment above (however, ...) */
                  if( isPseudocostUpdateValid(var, set, boundchgs[i].data.branchingdata.lpsolval, updateintegers, updatecontinuous) &&
                      (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS || !boundchgs[i].redundant || set->branch_lpgainnorm != 'd')
                    )
                  {
                     var->pseudocostflag = PSEUDOCOST_UPDATE; /*lint !e641*/
                     nvalidupdates++;
                  }
                  else
                     var->pseudocostflag = PSEUDOCOST_IGNORE; /*lint !e641*/
               }
            }
         }
      }

      /* update the pseudo cost values and reset the variables' flags; assume, that the responsibility for the dual gain
       * is equally spread on all bound changes that lead to valid pseudo cost updates
       */
      assert(SCIPnodeGetType(tree->focuslpstatefork) == SCIP_NODETYPE_FORK);
      weight = (nvalidupdates > 0 ? 1.0 / (SCIP_Real)nvalidupdates : 1.0);
      lpgain = (SCIPlpGetObjval(lp, set, prob) - tree->focuslpstatefork->data.fork->lpobjval) * weight;
      lpgain = MAX(lpgain, 0.0);

      for( i = 0; i < nupdates; ++i )
      {
         assert((SCIP_BOUNDCHGTYPE)updates[i]->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING);

         var = updates[i]->var;
         assert(var != NULL);
         assert((PSEUDOCOSTFLAG)var->pseudocostflag != PSEUDOCOST_NONE);

         if( (PSEUDOCOSTFLAG)var->pseudocostflag == PSEUDOCOST_UPDATE )
         {
            if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS || set->branch_lpgainnorm == 'l' )
            {
               SCIPsetDebugMsg(set, "updating pseudocosts of <%s>: sol: %g -> %g, LP: %e -> %e => solvaldelta = %g, gain=%g, weight: %g\n",
                  SCIPvarGetName(var), updates[i]->data.branchingdata.lpsolval, SCIPvarGetLPSol(var),
                  tree->focuslpstatefork->data.fork->lpobjval, SCIPlpGetObjval(lp, set, prob),
                  SCIPvarGetLPSol(var) - updates[i]->data.branchingdata.lpsolval, lpgain, weight);
               SCIP_CALL( SCIPvarUpdatePseudocost(var, set, stat,
                  SCIPvarGetLPSol(var) - updates[i]->data.branchingdata.lpsolval, lpgain, weight) );
            }
            else
            {
               /* set->branch_lpgainnorm == 'd':
                * For continuous variables, we want to pseudocosts to be the average of the gain in the LP value
                * if the domain is reduced from x% of its original width to y% of its original (e.g., global) width, i.e.,
                * to be the average of LPgain / (oldwidth/origwidth - newwidth/origwidth) = LPgain * origwidth / (oldwidth - newwidth).
                * Then an expected improvement in the LP value by a reduction of the domain width
                * from x% to y% of its original width can be computed by pseudocost * (oldwidth - newwidth) / origwidth.
                * Since the original width cancels out, we can also define the pseudocosts as average of LPgain / (oldwidth - newwidth)
                * and compute the expected improvement as pseudocost * (oldwidth - newwidth).
                *
                * Let var have bounds [a,c] before the branching and assume we branched on some value b.
                * b is given by updates[i]->newbound.
                *
                * If updates[i]->boundtype = upper, then node corresponds to the child [a,b].
                * Thus, we have oldwidth = c-a, newwidth = b-a, and oldwidth - newwidth = c-b.
                * To get c (the previous upper bound), we look into the var->ubchginfos array.
                *
                * If updates[i]->boundtype = lower, then node corresponds to the child [b,c].
                * Thus, we have oldwidth = c-a, newwidth = c-b, and oldwidth - newwidth = b-a.
                * To get c (the previous lower bound), we look into the var->lbchginfos array.
                */
               SCIP_BDCHGINFO* bdchginfo;
               SCIP_Real oldbound;
               SCIP_Real delta;
               int j;
               int nbdchginfos;

               assert(set->branch_lpgainnorm == 'd' || set->branch_lpgainnorm == 's');

               oldbound = SCIP_INVALID;

               if( set->branch_lpgainnorm == 'd' )
               {
                  assert(!updates[i]->redundant);

                  if( (SCIP_BOUNDTYPE)updates[i]->boundtype == SCIP_BOUNDTYPE_UPPER )
                  {
                     nbdchginfos = SCIPvarGetNBdchgInfosUb(var);

                     /* walk backwards through bound change information array to find the bound change corresponding to branching in updates[i]
                      * usually it will be the first one we look at */
                     for( j = nbdchginfos-1; j >= 0; --j )
                     {
                        bdchginfo = SCIPvarGetBdchgInfoUb(var, j);

                        if( bdchginfo->oldbound > updates[i]->newbound )
                        {
                           /* first boundchange which upper bound is above the upper bound set by the branching in updates[i]
                            * if bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING, then this should be exactly the bound change that we are looking for
                            * if bdchginfo->boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING, then this should be because the branching domain change has not been applied to the variable due to redundancy
                            * in this case, i.e., if there was another boundchange coming from somewhere else, I am not sure whether oldbound is an accurate value to compute the old domain size, so we skip the pseudocosts update
                            */
                           if( (SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
                           {
                              assert(bdchginfo->newbound == updates[i]->newbound); /*lint !e777*/
                              oldbound = bdchginfo->oldbound;
                           }
                           else
                              assert(updates[i]->redundant);

                           break;
                        }
                     }
                     /* if the bound change was redundant (e.g., due to a change in the global bound), then it was not applied, so there exists no corresponding bound change info
                      * if it is not redundant, then we should have found at least one corresponding boundchange */
                     assert(j >= 0 || updates[i]->redundant);
                     if( oldbound != SCIP_INVALID ) /*lint !e777*/
                     {
                        assert(!SCIPsetIsInfinity(set, -oldbound)); /* branching on a variable fixed to -infinity does not make sense */
                        assert(!SCIPsetIsInfinity(set, updates[i]->newbound)); /* branching to infinity does not make sense */

                        /* if the old upper bound is at infinity or the new upper bound is at -infinity, then we say the delta (c-b) is infinity */
                        if( SCIPsetIsInfinity(set, oldbound) || SCIPsetIsInfinity(set, -updates[i]->newbound) )
                           delta = SCIP_INVALID;
                        else
                           delta = updates[i]->newbound - oldbound;
                     }
                     else
                        delta = SCIP_INVALID;

                  }
                  else
                  {
                     assert((SCIP_BOUNDTYPE)updates[i]->boundtype == SCIP_BOUNDTYPE_LOWER);
                     nbdchginfos = SCIPvarGetNBdchgInfosLb(var);

                     /* walk backwards through bound change information array to find the bound change corresponding to branching in updates[i]
                      * usually it will be the first one we look at */
                     for( j = nbdchginfos-1; j >= 0; --j )
                     {
                        bdchginfo = SCIPvarGetBdchgInfoLb(var, j);

                        if( bdchginfo->oldbound < updates[i]->newbound )
                        {
                           /* first boundchange which lower bound is below the lower bound set by the branching in updates[i]
                            * if bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING, then this should be exactly the bound change that we are looking for
                            * if bdchginfo->boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING, then this should be because the branching domain change has not been applied to the variable due to redundancy
                            * in this case, i.e., if there was another boundchange coming from somewhere else, I am not sure whether oldbound is an accurate value to compute the old domain size, so we skip the pseudocosts update
                            */
                           if( (SCIP_BOUNDCHGTYPE)bdchginfo->boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING )
                           {
                              assert(bdchginfo->newbound == updates[i]->newbound); /*lint !e777*/
                              oldbound = bdchginfo->oldbound;
                           }
                           else
                              assert(updates[i]->redundant);

                           break;
                        }
                     }
                     /* if the bound change was redundant (e.g., due to a change in the global bound), then it was not applied, so there exists no corresponding bound change info
                      * if it is not redundant, then we should have found at least one corresponding boundchange */
                     assert(j >= 0 || updates[i]->redundant);
                     if( oldbound != SCIP_INVALID ) /*lint !e777*/
                     {
                        assert(!SCIPsetIsInfinity(set, oldbound)); /* branching on a variable fixed to +infinity does not make sense */
                        assert(!SCIPsetIsInfinity(set, -updates[i]->newbound)); /* branching to infinity does not make sense */

                        /* if the old lower bound is at -infinity or the new lower bound is at +infinity, then we say the delta (b-a) is infinity */
                        if( SCIPsetIsInfinity(set, -oldbound) || SCIPsetIsInfinity(set, updates[i]->newbound) )
                           delta = SCIP_INVALID;
                        else
                           delta = updates[i]->newbound - oldbound;
                     }
                     else
                        delta = SCIP_INVALID;
                  }
               }
               else
               {
                  /* set->branch_lpgainnorm == 's':
                   * Here, we divide the LPgain by the reduction in the sibling node.
                   *
                   * If updates[i]->boundtype = upper, then node corresponds to the child [a,b].
                   * Thus, we have oldwidth = c-a, newwidth = c-b, and oldwidth - newwidth = b-a.
                   * Conveniently, we just use the current lower bound for a (it may have been tightened, though).
                   *
                   * If updates[i]->boundtype = lower, then node corresponds to the child [b,a].
                   * Thus, we have oldwidth = c-a, newwidth = b-a, and oldwidth - newwidth = c-b.
                   * Conveniently, we just use the current upper bound for c (it may have been tightened, though).
                   */
                  if( (SCIP_BOUNDTYPE)updates[i]->boundtype == SCIP_BOUNDTYPE_UPPER )
                  {
                     assert(!SCIPsetIsInfinity(set, updates[i]->newbound)); /* branching on a variable fixed to +infinity does not make sense */
                     assert(!SCIPsetIsInfinity(set, SCIPvarGetLbLocal(var))); /* branching to infinity does not make sense */
                     if( SCIPsetIsInfinity(set, -updates[i]->newbound) || SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)) )
                        delta = SCIP_INVALID;
                     else
                        delta = updates[i]->newbound - SCIPvarGetLbLocal(var);
                  }
                  else
                  {
                     assert((SCIP_BOUNDTYPE)updates[i]->boundtype == SCIP_BOUNDTYPE_LOWER);
                     assert(!SCIPsetIsInfinity(set, -updates[i]->newbound)); /* branching on a variable fixed to -infinity does not make sense */
                     assert(!SCIPsetIsInfinity(set, -SCIPvarGetUbLocal(var))); /* branching to -infinity does not make sense */
                     if( SCIPsetIsInfinity(set, updates[i]->newbound) || SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var)) )
                        delta = SCIP_INVALID;
                     else
                        delta = -(SCIPvarGetUbLocal(var) - updates[i]->newbound);
                  }
               }

               if( delta != SCIP_INVALID ) /*lint !e777*/
               {
                  SCIPsetDebugMsg(set, "updating pseudocosts of <%s> with strategy %c: domain: [%g,%g] -> [%g,%g], LP: %e -> %e => "
                     "delta = %g, gain=%g, weight: %g\n",
                     SCIPvarGetName(var), set->branch_lpgainnorm,
                     (SCIP_BOUNDTYPE)updates[i]->boundtype == SCIP_BOUNDTYPE_UPPER ? SCIPvarGetLbLocal(var) : oldbound,
                     (SCIP_BOUNDTYPE)updates[i]->boundtype == SCIP_BOUNDTYPE_UPPER ? oldbound : SCIPvarGetUbLocal(var),
                     (SCIP_BOUNDTYPE)updates[i]->boundtype == SCIP_BOUNDTYPE_UPPER ? SCIPvarGetLbLocal(var) : updates[i]->newbound,
                     (SCIP_BOUNDTYPE)updates[i]->boundtype == SCIP_BOUNDTYPE_UPPER ? updates[i]->newbound : SCIPvarGetUbLocal(var),
                     tree->focuslpstatefork->lowerbound, SCIPlpGetObjval(lp, set, prob),
                     delta, lpgain, weight);

                  SCIP_CALL( SCIPvarUpdatePseudocost(var, set, stat, delta, lpgain, weight) );
               }
            }
         }
         var->pseudocostflag = PSEUDOCOST_NONE; /*lint !e641*/
      }

      /* free the buffer for the collected bound changes */
      SCIPsetFreeBufferArray(set, &updates);
   }

   return SCIP_OKAY;
}

/** updates the estimated value of a primal feasible solution for the focus node after the LP was solved */
static
SCIP_RETCODE updateEstimate(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   SCIP_NODE* focusnode;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   SCIP_Real estimate;
   int nlpcands;
   int i;

   /* estimate is only available if LP was solved to optimality */
   if( !SCIPtreeHasFocusNodeLP(tree) || SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL || !SCIPlpIsRelax(lp) )
      return SCIP_OKAY;

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(focusnode != NULL);

   /* get the fractional variables */
   SCIP_CALL( SCIPbranchcandGetLPCands(branchcand, set, stat, lp, &lpcands, NULL, &lpcandsfrac, &nlpcands, NULL, NULL) );

   /* calculate the estimate: lowerbound + sum(min{f_j * pscdown_j, (1-f_j) * pscup_j}) */
   estimate = SCIPnodeGetLowerbound(focusnode);

   /* an infinite lower bound implies an infinite estimate */
   if( SCIPsetIsInfinity(set, estimate) )
   {
      SCIPnodeSetEstimate(focusnode, set, estimate);
      return SCIP_OKAY;
   }

   for( i = 0; i < nlpcands; ++i )
   {
      SCIP_Real pscdown;
      SCIP_Real pscup;

      pscdown = SCIPvarGetPseudocost(lpcands[i], stat, 0.0-lpcandsfrac[i]);
      pscup = SCIPvarGetPseudocost(lpcands[i], stat, 1.0-lpcandsfrac[i]);
      estimate += MIN(pscdown, pscup);
   }
   SCIPnodeSetEstimate(focusnode, set, estimate);

   return SCIP_OKAY;
}

/** puts all constraints with initial flag TRUE into the LP */
SCIP_RETCODE SCIPinitConssLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cutpool */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             root,               /**< is this the initial root LP? */
   SCIP_Bool             firstsubtreeinit,   /**< is this the first call in the current subtree after jumping through the tree? */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   int h;

   assert(set != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* inform separation storage, that LP is now filled with initial data */
   SCIPsepastoreStartInitialLP(sepastore);

   /* add LP relaxations of all initial constraints to LP */
   SCIPsetDebugMsg(set, "init LP: initial rows\n");
   for( h = 0; h < set->nconshdlrs && !(*cutoff); ++h )
   {
      SCIP_CALL( SCIPconshdlrInitLP(set->conshdlrs[h], blkmem, set, stat, tree, firstsubtreeinit, cutoff) );
   }

   if( set->reopt_enable && set->reopt_usecuts && firstsubtreeinit && !(*cutoff) )
   {
      /* add stored cuts from last reoptimization run */
      SCIP_CALL( SCIPreoptApplyCuts(reopt, tree->focusnode, sepastore, cutpool, blkmem, set, stat, eventqueue,
            eventfilter, lp, root) );
   }

   if( !(*cutoff) )
   {
      /* apply cuts */
      SCIP_CALL( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
            eventqueue, eventfilter, cliquetable, root, SCIP_EFFICIACYCHOICE_LP, cutoff) );
   }
   else
   {
      /* the current node will be cut off; we clear the sepastore */
      SCIP_CALL( SCIPsepastoreClearCuts(sepastore, blkmem, set, eventqueue, eventfilter, lp) );
   }


   /* inform separation storage, that initial LP setup is now finished */
   SCIPsepastoreEndInitialLP(sepastore);

  return SCIP_OKAY;
}

/** constructs the initial LP of the current node */
static
SCIP_RETCODE initLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             root,               /**< is this the initial root LP? */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_VAR* var;
   int v;

   assert(set != NULL);
   assert(transprob != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* at the root node, we have to add the initial variables as columns */
   if( root )
   {
      assert(SCIPlpGetNCols(lp) == 0);
      assert(SCIPlpGetNRows(lp) == 0);
      assert(lp->nremovablecols == 0);
      assert(lp->nremovablerows == 0);

      /* inform pricing storage, that LP is now filled with initial data */
      SCIPpricestoreStartInitialLP(pricestore);

      /* add all initial variables to LP */
      SCIPsetDebugMsg(set, "init LP: initial columns\n");
      for( v = 0; v < transprob->nvars && !(*cutoff); ++v )
      {
         var = transprob->vars[v];
         assert(SCIPvarGetProbindex(var) >= 0);

         if( SCIPvarIsInitial(var) )
         {
            SCIP_CALL( SCIPpricestoreAddVar(pricestore, blkmem, set, eventqueue, lp, var, 0.0, TRUE) );
         }

         /* check for empty domains (necessary if no presolving was performed) */
         if( SCIPsetIsGT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            *cutoff = TRUE;
      }
      assert(lp->nremovablecols == 0);
      SCIP_CALL( SCIPpricestoreApplyVars(pricestore, blkmem, set, stat, eventqueue, transprob, tree, lp) );

      /* inform pricing storage, that initial LP setup is now finished */
      SCIPpricestoreEndInitialLP(pricestore);
   }

   if( *cutoff )
      return SCIP_OKAY;

   /* put all initial constraints into the LP */
   /* @todo check whether we jumped through the tree */
   SCIP_CALL( SCIPinitConssLP(blkmem, set, sepastore, cutpool, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
         eventfilter, cliquetable, root, TRUE, cutoff) );

   return SCIP_OKAY;
}

/** constructs the LP of the current node, but does not load the LP state and warmstart information  */
SCIP_RETCODE SCIPconstructCurrentLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cutpool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             newinitconss,       /**< do we have to add new initial constraints? */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_Bool initroot = FALSE;

   assert(tree != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   if( !SCIPtreeIsFocusNodeLPConstructed(tree) )
   {
      /* inform separation storage, that LP is now filled with initial data */
      SCIPsepastoreStartInitialLP(sepastore);

      if( tree->correctlpdepth >= 0 )
      {
         int i;

         for( i = tree->pathnlprows[tree->correctlpdepth]; i < lp->nrows; ++i )
         {
            /* keep all active global cuts that where applied in the previous node in the lp */
            if( !lp->rows[i]->local && lp->rows[i]->age == 0 )
            {
               SCIP_CALL( SCIPsepastoreAddCut(sepastore, blkmem, set, stat, eventqueue, eventfilter, lp, lp->rows[i], TRUE, (SCIPtreeGetCurrentDepth(tree) == 0), cutoff) );
            }
         }
      }

      if( !(*cutoff) )
      {
         /* load the LP into the solver and load the LP state */
         SCIPsetDebugMsg(set, "loading LP\n");
         SCIP_CALL( SCIPtreeLoadLP(tree, blkmem, set, eventqueue, eventfilter, lp, &initroot) );
         assert(initroot || SCIPnodeGetDepth(SCIPtreeGetFocusNode(tree)) > 0);
         assert(SCIPtreeIsFocusNodeLPConstructed(tree));

         /* apply cuts */
         SCIP_CALL( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
                                           eventqueue, eventfilter, cliquetable, (SCIPtreeGetCurrentDepth(tree) == 0), SCIP_EFFICIACYCHOICE_LP, cutoff) );
      }
      else
      {
         /* the current node will be cut off; we clear the sepastore */
         SCIP_CALL( SCIPsepastoreClearCuts(sepastore, blkmem, set, eventqueue, eventfilter, lp) );
      }

      /* inform separation storage, that initial LP setup is now finished */
      SCIPsepastoreEndInitialLP(sepastore);

      if( !(*cutoff) )
      {
         /* setup initial LP relaxation of node */
         SCIP_CALL( initLP(blkmem, set, stat, transprob, origprob, tree, reopt, lp, pricestore, sepastore, cutpool, branchcand,
                           eventqueue, eventfilter, cliquetable, initroot, cutoff) );
      }
   }
   else if( newinitconss )
   {
      SCIP_CALL( SCIPinitConssLP(blkmem, set, sepastore, cutpool, stat, transprob,
            origprob, tree, reopt, lp, branchcand, eventqueue, eventfilter, cliquetable, FALSE, FALSE,
            cutoff) );
   }

   return SCIP_OKAY;
}

/** updates the primal ray stored in primal data
 * clears previously stored primal ray, if existing and there was no LP error
 * stores current primal ray, if LP is unbounded and there has been no error
 */
static
SCIP_RETCODE updatePrimalRay(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             lperror             /**< has there been an LP error? */
   )
{
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(primal != NULL);
   assert(tree != NULL);
   assert(lp != NULL);

   if( lperror )
      return SCIP_OKAY;

   /* clear previously stored primal ray, if any */
   if( primal->primalray != NULL )
   {
      SCIP_CALL( SCIPsolFree(&primal->primalray, blkmem, primal) );
   }

   /* store unbounded ray, if LP is unbounded */
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIP_VAR** vars;
      SCIP_Real* ray;
      int nvars;
      int i;

      SCIPsetDebugMsg(set, "LP is unbounded, store primal ray\n");

      vars = prob->vars;
      nvars = prob->nvars;

      /* get buffer memory for storing the ray and load the ray values into it */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &ray, nvars) );
      BMSclearMemoryArray(ray, nvars);
      SCIP_CALL( SCIPlpGetPrimalRay(lp, set, ray) );

      /* create solution to store the primal ray in */
      assert(primal->primalray == NULL);
      SCIP_CALL( SCIPsolCreate(&primal->primalray, blkmem, set, stat, primal, tree, NULL) );

      /* set values of all active variable in the solution that represents the primal ray */
      for( i = 0; i < nvars; i++ )
      {
         SCIP_CALL( SCIPsolSetVal(primal->primalray, set, stat, tree, vars[i], ray[i]) );
      }

      SCIPdebug( SCIP_CALL( SCIPprintRay(set->scip, primal->primalray, NULL, FALSE) ) );

      /* free memory for buffering the ray values */
      SCIPsetFreeBufferArray(set, &ray);
   }

   return SCIP_OKAY;
}

/** load and solve the initial LP of a node */
static
SCIP_RETCODE solveNodeInitialLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cutpool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             newinitconss,       /**< do we have to add new initial constraints? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved error in LP solving occured */
   )
{
   /* initializing variables for compiler warnings, which are not correct */
   SCIP_Real starttime = 0.0;
   SCIP_Longint nlpiterations = 0;
   SCIP_NODE* focusnode;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);
   assert(lperror != NULL);
   assert(SCIPtreeGetFocusNode(tree) != NULL);
   assert(SCIPnodeGetType(SCIPtreeGetFocusNode(tree)) == SCIP_NODETYPE_FOCUSNODE);

   *cutoff = FALSE;
   *lperror = FALSE;

   /* load the LP into the solver */
   SCIP_CALL( SCIPconstructCurrentLP(blkmem, set, stat, transprob, origprob, tree, reopt, lp, pricestore, sepastore, cutpool,
         branchcand, eventqueue, eventfilter, cliquetable, newinitconss, cutoff) );

   if( *cutoff )
      return SCIP_OKAY;

   /* load the LP state */
   SCIP_CALL( SCIPtreeLoadLPState(tree, blkmem, set, stat, eventqueue, lp) );

   focusnode = SCIPtreeGetFocusNode(tree);

   /* store current LP iteration count and solving time if we are at the root node */
   if( focusnode->depth == 0 )
   {
      nlpiterations = stat->nlpiterations;
      starttime = SCIPclockGetTime(stat->solvingtime);
   }

   /* solve initial LP */
   SCIPsetDebugMsg(set, "node: solve initial LP\n");
   SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob,
         SCIPnodeGetDepth(SCIPtreeGetFocusNode(tree)) == 0 ? set->lp_rootiterlim : set->lp_iterlim, TRUE, TRUE, FALSE, lperror) );
   assert(lp->flushed);
   assert(lp->solved || *lperror);

   /* save time for very first LP in root node */
   if ( stat->nnodelps == 0 && focusnode->depth == 0 )
   {
      stat->firstlptime = SCIPclockGetTime(stat->solvingtime) - starttime;
   }

   /* remove previous primal ray, store new one if LP is unbounded */
   SCIP_CALL( updatePrimalRay(blkmem, set, stat, transprob, primal, tree, lp, *lperror) );

   if( !(*lperror) )
   {
      /* cppcheck-suppress unassignedVariable */
      SCIP_EVENT event;

      if( SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_ITERLIMIT && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_TIMELIMIT )
      {
         /* issue FIRSTLPSOLVED event */
         SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_FIRSTLPSOLVED) );
         SCIP_CALL( SCIPeventChgNode(&event, SCIPtreeGetFocusNode(tree)) );
         SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
      }

      /* update pseudo cost values for integer variables (always) and for continuous variables (if not delayed) */
      SCIP_CALL( updatePseudocost(set, stat, transprob, tree, lp, TRUE, !set->branch_delaypscost) );

      /* update lower bound of current node w.r.t. initial lp */
      assert(!(*cutoff));
      if( (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY
	    || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT)
	 && SCIPprobAllColsInLP(transprob, set, lp) && SCIPlpIsRelax(lp) )
      {
	 SCIP_CALL( SCIPnodeUpdateLowerboundLP(focusnode, set, stat, tree, transprob, origprob, lp) );

         /* if this is the first LP solved at the root, store its iteration count and solution value */
         if( stat->nnodelps == 0 && focusnode->depth == 0 )
         {
            SCIP_Real lowerbound;

            assert(stat->nrootfirstlpiterations == 0);
            stat->nrootfirstlpiterations = stat->nlpiterations - nlpiterations;

            if( set->misc_exactsolve )
            {
               SCIP_CALL( SCIPlpGetProvedLowerbound(lp, set, &lowerbound) );
            }
            else
               lowerbound = SCIPlpGetObjval(lp, set, transprob);

            stat->firstlpdualbound = SCIPprobExternObjval(transprob, origprob, set, lowerbound);
         }
      }
   }

   return SCIP_OKAY;
}

/** makes sure the LP is flushed and solved */
static
SCIP_RETCODE separationRoundResolveLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved error in LP solving occured */
   SCIP_Bool*            mustsepa,           /**< pointer to store TRUE if additional separation rounds should be performed */
   SCIP_Bool*            mustprice           /**< pointer to store TRUE if additional pricing rounds should be performed */
   )
{
   assert(lp != NULL);
   assert(lperror != NULL);
   assert(mustsepa != NULL);
   assert(mustprice != NULL);

   /* if bound changes were applied in the separation round, we have to resolve the LP */
   if( !lp->flushed )
   {
      /* solve LP (with dual simplex) */
      SCIPsetDebugMsg(set, "separation: resolve LP\n");
      SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, prob, set->lp_iterlim, FALSE, TRUE, FALSE, lperror) );
      assert(lp->flushed);
      assert(lp->solved || *lperror);
      *mustsepa = TRUE;
      *mustprice = TRUE;

      /* remove previous primal ray, store new one if LP is unbounded */
      SCIP_CALL( updatePrimalRay(blkmem, set, stat, prob, primal, tree, lp, *lperror) );
   }

   return SCIP_OKAY;
}

/** applies one round of LP separation */
static
SCIP_RETCODE separationRoundLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   actdepth,           /**< current depth in the tree */
   SCIP_Real             bounddist,          /**< current relative distance of local dual bound to global dual bound */
   SCIP_Bool             allowlocal,         /**< should the separators be asked to separate local cuts */
   SCIP_Bool             onlydelayed,        /**< should only delayed separators be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            enoughcuts,         /**< pointer to store whether enough cuts have been found this round */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved error in LP solving occured */
   SCIP_Bool*            mustsepa,           /**< pointer to store TRUE if additional separation rounds should be performed */
   SCIP_Bool*            mustprice           /**< pointer to store TRUE if additional pricing rounds should be performed */
   )
{
   SCIP_RESULT result;
   int i;
   SCIP_Bool consadded;
   SCIP_Bool root;

   assert(set != NULL);
   assert(lp != NULL);
   assert(set->conshdlrs_sepa != NULL);
   assert(delayed != NULL);
   assert(enoughcuts != NULL);
   assert(cutoff != NULL);
   assert(lperror != NULL);

   root = (actdepth == 0);
   *delayed = FALSE;
   *enoughcuts = (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root));
   *lperror = FALSE;
   consadded = FALSE;

   SCIPsetDebugMsg(set, "calling separators on LP solution in depth %d (onlydelayed: %u)\n", actdepth, onlydelayed);

   /* sort separators by priority */
   SCIPsetSortSepas(set);

   /* call LP separators with nonnegative priority */
   for( i = 0; i < set->nsepas && !(*cutoff) && !(*lperror) && !(*enoughcuts) && lp->flushed && lp->solved
           && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
        ++i )
   {
#ifndef NDEBUG
      size_t nusedbuffer = BMSgetNUsedBufferMemory(SCIPbuffer(set->scip));
#endif

      if( SCIPsepaGetPriority(set->sepas[i]) < 0 )
         continue;

      if( onlydelayed && !SCIPsepaWasLPDelayed(set->sepas[i]) )
         continue;

      SCIPsetDebugMsg(set, " -> executing separator <%s> with priority %d\n",
         SCIPsepaGetName(set->sepas[i]), SCIPsepaGetPriority(set->sepas[i]));
      SCIP_CALL( SCIPsepaExecLP(set->sepas[i], set, stat, sepastore, actdepth, bounddist, allowlocal, onlydelayed, &result) );
#ifndef NDEBUG
      if( BMSgetNUsedBufferMemory(SCIPbuffer(set->scip)) > nusedbuffer )
      {
         SCIPerrorMessage("Buffer not completely freed after executing separator <%s>\n", SCIPsepaGetName(set->sepas[i]));
         SCIPABORT();
      }
#endif
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
      *delayed = *delayed || (result == SCIP_DELAYED);

      if( !(*cutoff) )
      {
         /* make sure the LP is solved (after adding bound changes, LP has to be flushed and resolved) */
         SCIP_CALL( separationRoundResolveLP(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, prob, primal, tree, lp, lperror, mustsepa, mustprice) );
      }
      else
      {
         SCIPsetDebugMsg(set, " -> separator <%s> detected cutoff\n", SCIPsepaGetName(set->sepas[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED || result == SCIP_NEWROUND) )
      {
         SCIPsetDebugMsg(set, " -> delayed separator <%s> found a cut\n", SCIPsepaGetName(set->sepas[i]));
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* try separating constraints of the constraint handlers */
   for( i = 0; i < set->nconshdlrs && !(*cutoff) && !(*lperror) && !(*enoughcuts) && lp->flushed && lp->solved
           && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
        ++i )
   {
      if( onlydelayed && !SCIPconshdlrWasLPSeparationDelayed(set->conshdlrs_sepa[i]) )
         continue;

      SCIPsetDebugMsg(set, " -> executing separation of constraint handler <%s> with priority %d\n",
         SCIPconshdlrGetName(set->conshdlrs_sepa[i]), SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i]));
      SCIP_CALL( SCIPconshdlrSeparateLP(set->conshdlrs_sepa[i], blkmem, set, stat, sepastore, actdepth, onlydelayed,
            &result) );

      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
      *delayed = *delayed || (result == SCIP_DELAYED);

      if( !(*cutoff) )
      {
         /* make sure the LP is solved (after adding bound changes, LP has to be flushed and resolved) */
         SCIP_CALL( separationRoundResolveLP(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, prob, primal, tree, lp, lperror, mustsepa, mustprice) );
      }
      else
      {
         SCIPsetDebugMsg(set, " -> constraint handler <%s> detected cutoff in separation\n", SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED || result == SCIP_NEWROUND) )
      {
         SCIPsetDebugMsg(set, " -> delayed constraint handler <%s> found a cut\n",
            SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* call LP separators with negative priority */
   for( i = 0; i < set->nsepas && !(*cutoff) && !(*lperror) && !(*enoughcuts) && lp->flushed && lp->solved
           && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
        ++i )
   {
      if( SCIPsepaGetPriority(set->sepas[i]) >= 0 )
         continue;

      if( onlydelayed && !SCIPsepaWasLPDelayed(set->sepas[i]) )
         continue;

      SCIPsetDebugMsg(set, " -> executing separator <%s> with priority %d\n",
         SCIPsepaGetName(set->sepas[i]), SCIPsepaGetPriority(set->sepas[i]));
      SCIP_CALL( SCIPsepaExecLP(set->sepas[i], set, stat, sepastore, actdepth, bounddist, allowlocal, onlydelayed, &result) );

      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
      *delayed = *delayed || (result == SCIP_DELAYED);

      if( !(*cutoff) )
      {
         /* make sure the LP is solved (after adding bound changes, LP has to be flushed and resolved) */
         SCIP_CALL( separationRoundResolveLP(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, prob, primal, tree, lp, lperror, mustsepa, mustprice) );
      }
      else
      {
         SCIPsetDebugMsg(set, " -> separator <%s> detected cutoff\n", SCIPsepaGetName(set->sepas[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED || result == SCIP_NEWROUND) )
      {
         SCIPsetDebugMsg(set, " -> delayed separator <%s> found a cut\n", SCIPsepaGetName(set->sepas[i]));
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* process the constraints that were added during this separation round */
   while( consadded )
   {
      assert(!onlydelayed);
      consadded = FALSE;

      for( i = 0; i < set->nconshdlrs && !(*cutoff) && !(*lperror) && !(*enoughcuts) && lp->flushed && lp->solved
              && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
           ++i )
      {
         SCIPsetDebugMsg(set, " -> executing separation of constraint handler <%s> with priority %d\n",
            SCIPconshdlrGetName(set->conshdlrs_sepa[i]), SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i]));
         SCIP_CALL( SCIPconshdlrSeparateLP(set->conshdlrs_sepa[i], blkmem, set, stat, sepastore, actdepth, onlydelayed,
            &result) );

         *cutoff = *cutoff || (result == SCIP_CUTOFF);
         consadded = consadded || (result == SCIP_CONSADDED);
         *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
         *delayed = *delayed || (result == SCIP_DELAYED);

         if( !(*cutoff) )
         {
            /* make sure the LP is solved (after adding bound changes, LP has to be flushed and resolved) */
            SCIP_CALL( separationRoundResolveLP(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, prob, primal, tree, lp, lperror, mustsepa, mustprice) );
         }
         else
         {
            SCIPsetDebugMsg(set, " -> constraint handler <%s> detected cutoff in separation\n", SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
         }
      }
   }

   SCIPsetDebugMsg(set, " -> separation round finished: delayed=%u, enoughcuts=%u, lpflushed=%u, cutoff=%u\n",
      *delayed, *enoughcuts, lp->flushed, *cutoff);

   return SCIP_OKAY;
}

/** applies one round of separation on the given primal solution */
static
SCIP_RETCODE separationRoundSol(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   int                   actdepth,           /**< current depth in the tree */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             onlydelayed,        /**< should only delayed separators be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            enoughcuts,         /**< pointer to store whether enough cuts have been found this round */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_RESULT result;
   int i;
   SCIP_Bool consadded;
   SCIP_Bool root;

   assert(set != NULL);
   assert(set->conshdlrs_sepa != NULL);
   assert(delayed != NULL);
   assert(enoughcuts != NULL);
   assert(cutoff != NULL);

   *delayed = FALSE;
   *enoughcuts = FALSE;
   consadded = FALSE;
   root = (actdepth == 0);

   SCIPsetDebugMsg(set, "calling separators on primal solution in depth %d (onlydelayed: %u)\n", actdepth, onlydelayed);

   /* sort separators by priority */
   SCIPsetSortSepas(set);

   /* call separators with nonnegative priority */
   for( i = 0; i < set->nsepas && !(*cutoff) && !(*enoughcuts) && !SCIPsolveIsStopped(set, stat, FALSE); ++i )
   {
      if( SCIPsepaGetPriority(set->sepas[i]) < 0 )
         continue;

      if( onlydelayed && !SCIPsepaWasSolDelayed(set->sepas[i]) )
         continue;

      SCIP_CALL( SCIPsepaExecSol(set->sepas[i], set, stat, sepastore, sol, actdepth, allowlocal, onlydelayed, &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPsetDebugMsg(set, " -> separator <%s> detected cutoff\n", SCIPsepaGetName(set->sepas[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED || result == SCIP_NEWROUND) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* try separating constraints of the constraint handlers */
   for( i = 0; i < set->nconshdlrs && !(*cutoff) && !(*enoughcuts) && !SCIPsolveIsStopped(set, stat, FALSE); ++i )
   {
      if( onlydelayed && !SCIPconshdlrWasSolSeparationDelayed(set->conshdlrs_sepa[i]) )
         continue;

      SCIP_CALL( SCIPconshdlrSeparateSol(set->conshdlrs_sepa[i], blkmem, set, stat, sepastore, sol, actdepth, onlydelayed,
            &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPsetDebugMsg(set, " -> constraint handler <%s> detected cutoff in separation\n",
            SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED || result == SCIP_NEWROUND) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* call separators with negative priority */
   for( i = 0; i < set->nsepas && !(*cutoff) && !(*enoughcuts) && !SCIPsolveIsStopped(set, stat, FALSE); ++i )
   {
      if( SCIPsepaGetPriority(set->sepas[i]) >= 0 )
         continue;

      if( onlydelayed && !SCIPsepaWasSolDelayed(set->sepas[i]) )
         continue;

      SCIP_CALL( SCIPsepaExecSol(set->sepas[i], set, stat, sepastore, sol, actdepth, allowlocal, onlydelayed, &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      consadded = consadded || (result == SCIP_CONSADDED);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
      *delayed = *delayed || (result == SCIP_DELAYED);
      if( *cutoff )
      {
         SCIPsetDebugMsg(set, " -> separator <%s> detected cutoff\n", SCIPsepaGetName(set->sepas[i]));
      }

      /* if we work off the delayed separators, we stop immediately if a cut was found */
      if( onlydelayed && (result == SCIP_CONSADDED || result == SCIP_REDUCEDDOM || result == SCIP_SEPARATED || result == SCIP_NEWROUND) )
      {
         *delayed = TRUE;
         return SCIP_OKAY;
      }
   }

   /* process the constraints that were added during this separation round */
   while( consadded )
   {
      assert(!onlydelayed);
      consadded = FALSE;

      for( i = 0; i < set->nconshdlrs && !(*cutoff) && !(*enoughcuts) && !SCIPsolveIsStopped(set, stat, FALSE); ++i )
      {
         SCIP_CALL( SCIPconshdlrSeparateSol(set->conshdlrs_sepa[i], blkmem, set, stat, sepastore, sol, actdepth, onlydelayed, &result) );
         *cutoff = *cutoff || (result == SCIP_CUTOFF);
         consadded = consadded || (result == SCIP_CONSADDED);
         *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
         *delayed = *delayed || (result == SCIP_DELAYED);
         if( *cutoff )
         {
            SCIPsetDebugMsg(set, " -> constraint handler <%s> detected cutoff in separation\n",
               SCIPconshdlrGetName(set->conshdlrs_sepa[i]));
         }
      }
   }

   SCIPsetDebugMsg(set, " -> separation round finished: delayed=%u, enoughcuts=%u, cutoff=%u\n",
      *delayed, *enoughcuts, *cutoff);

   return SCIP_OKAY;
}

/** applies one round of separation on the given primal solution or on the LP solution */
SCIP_RETCODE SCIPseparationRound(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   int                   actdepth,           /**< current depth in the tree */
   SCIP_Bool             allowlocal,         /**< should the separator be asked to separate local cuts */
   SCIP_Bool             onlydelayed,        /**< should only delayed separators be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_Bool enoughcuts;

   assert(delayed != NULL);
   assert(cutoff != NULL);

   *delayed = FALSE;
   *cutoff = FALSE;
   enoughcuts = FALSE;

   if( sol == NULL )
   {
      SCIP_Bool lperror;
      SCIP_Bool mustsepa;
      SCIP_Bool mustprice;

      /* apply a separation round on the LP solution */
      lperror = FALSE;
      mustsepa = FALSE;
      mustprice = FALSE;
      SCIP_CALL( separationRoundLP(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, prob, primal, tree, lp, sepastore, \
                                   actdepth, 0.0, allowlocal, onlydelayed, delayed, &enoughcuts, cutoff, \
                                   &lperror, &mustsepa, &mustprice) );
   }
   else
   {
      /* apply a separation round on the given primal solution */
      SCIP_CALL( separationRoundSol(blkmem, set, stat, sepastore, sol, actdepth, allowlocal, onlydelayed, delayed, &enoughcuts, cutoff) );
   }

   return SCIP_OKAY;
}

/** solves the current LP completely with pricing in new variables */
SCIP_RETCODE SCIPpriceLoop(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cutpool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             pretendroot,        /**< should the pricers be called as if we are at the root node? */
   SCIP_Bool             displayinfo,        /**< should info lines be displayed after each pricing round? */
   int                   maxpricerounds,     /**< maximal number of pricing rounds (-1: no limit);
                                              *   a finite limit means that the LP might not be solved to optimality! */
   int*                  npricedcolvars,     /**< pointer to store number of column variables after problem vars were priced */
   SCIP_Bool*            mustsepa,           /**< pointer to store TRUE if a separation round should follow */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved error in LP solving occured */
   SCIP_Bool*            aborted             /**< pointer to store whether the pricing was aborted and the lower bound must
                                              *   not be used */
   )
{
   SCIP_NODE* currentnode;
   int npricerounds;
   SCIP_Bool mustprice;
   SCIP_Bool cutoff;
   SCIP_Bool unbounded;

   assert(transprob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(npricedcolvars != NULL);
   assert(mustsepa != NULL);
   assert(lperror != NULL);
   assert(aborted != NULL);

   currentnode = SCIPtreeGetCurrentNode(tree);
   assert(currentnode == SCIPtreeGetFocusNode(tree) || SCIPtreeProbing(tree));
   *npricedcolvars = transprob->ncolvars;
   *lperror = FALSE;
   *aborted = FALSE;

   /* if the LP is unbounded, we don't need to price */
   mustprice = (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL
      || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE
      || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT);

   /* if all the variables are already in the LP, we don't need to price */
   mustprice = mustprice && !SCIPprobAllColsInLP(transprob, set, lp);

   /* check if infinite number of pricing rounds should be used */
   if( maxpricerounds == -1 )
      maxpricerounds = INT_MAX;

   /* pricing (has to be done completely to get a valid lower bound) */
   npricerounds = 0;
   while( !(*lperror) && mustprice && npricerounds < maxpricerounds )
   {
      SCIP_Bool enoughvars;
      SCIP_RESULT result;
      SCIP_Real lb;
      SCIP_Bool foundsol;
      SCIP_Bool stopearly;
      SCIP_Bool stoppricing;
      int p;

      assert(lp->flushed);
      assert(lp->solved);
      assert(SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY);

      /* check if pricing loop should be aborted */
      if( SCIPsolveIsStopped(set, stat, FALSE) )
      {
         /* do not print the warning message if we stopped because the problem is solved */
         if( !SCIPsetIsLE(set, SCIPgetUpperbound(set->scip), SCIPgetLowerbound(set->scip)) )
            SCIPmessagePrintWarning(messagehdlr, "pricing has been interrupted -- LP of current node is invalid\n");

         *aborted = TRUE;
         break;
      }

      /* call primal heuristics which are callable during pricing */
      SCIP_CALL( SCIPprimalHeuristics(set, stat, transprob, primal, tree, lp, NULL, SCIP_HEURTIMING_DURINGPRICINGLOOP,
            FALSE, &foundsol, &unbounded) );

      /* price problem variables */
      SCIPsetDebugMsg(set, "problem variable pricing\n");
      assert(SCIPpricestoreGetNVars(pricestore) == 0);
      assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
      SCIP_CALL( SCIPpricestoreAddProbVars(pricestore, blkmem, set, stat, transprob, tree, lp, branchcand, eventqueue) );
      *npricedcolvars = transprob->ncolvars;

      /* call external pricers to create additional problem variables */
      SCIPsetDebugMsg(set, "external variable pricing\n");

      /* sort pricer algorithms by priority */
      SCIPsetSortPricers(set);

      /* call external pricer algorithms, that are active for the current problem */
      enoughvars = (SCIPsetGetPriceMaxvars(set, pretendroot) == INT_MAX ?
            FALSE : SCIPpricestoreGetNVars(pricestore) >= SCIPsetGetPriceMaxvars(set, pretendroot) + 1);
      stoppricing = FALSE;
      for( p = 0; p < set->nactivepricers && !enoughvars; ++p )
      {
         SCIP_CALL( SCIPpricerExec(set->pricers[p], set, transprob, lp, pricestore, &lb, &stopearly, &result) );
         assert(result == SCIP_DIDNOTRUN || result == SCIP_SUCCESS);
         SCIPsetDebugMsg(set, "pricing: pricer %s returned result = %s, lowerbound = %f\n",
            SCIPpricerGetName(set->pricers[p]), (result == SCIP_DIDNOTRUN ? "didnotrun" : "success"), lb);
         enoughvars = enoughvars || (SCIPsetGetPriceMaxvars(set, pretendroot) == INT_MAX ?
               FALSE : SCIPpricestoreGetNVars(pricestore) >= SCIPsetGetPriceMaxvars(set, pretendroot) + 1);
         *aborted = ( (*aborted) || (result == SCIP_DIDNOTRUN) );

         /* set stoppricing to TRUE, of the first pricer wants to stop pricing */
         if( p == 0 && stopearly )
            stoppricing = TRUE;

         /* stoppricing only remains TRUE, if all other pricers want to stop pricing as well */
         if( stoppricing && !stopearly )
            stoppricing = FALSE;

         /* update lower bound w.r.t. the lower bound given by the pricer */
         SCIPnodeUpdateLowerbound(currentnode, stat, set, tree, transprob, origprob, lb);
         SCIPsetDebugMsg(set, " -> new lower bound given by pricer %s: %g\n", SCIPpricerGetName(set->pricers[p]), lb);
      }

      /* apply the priced variables to the LP */
      SCIP_CALL( SCIPpricestoreApplyVars(pricestore, blkmem, set, stat, eventqueue, transprob, tree, lp) );
      assert(SCIPpricestoreGetNVars(pricestore) == 0);
      assert(!lp->flushed || lp->solved);
      mustprice = !lp->flushed || (transprob->ncolvars != *npricedcolvars);
      *mustsepa = *mustsepa || !lp->flushed;

      /* after adding columns, the LP should be primal feasible such that the primal simplex is applicable;
       * if LP was infeasible, we have to use dual simplex
       */
      SCIPsetDebugMsg(set, "pricing: solve LP\n");
      SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob, -1LL, FALSE, TRUE, FALSE, lperror) );
      assert(lp->flushed);
      assert(lp->solved || *lperror);

      /* reset bounds temporarily set by pricer to their original values */
      SCIPsetDebugMsg(set, "pricing: reset bounds\n");
      SCIP_CALL( SCIPpricestoreResetBounds(pricestore, blkmem, set, stat, lp, branchcand, eventqueue) );
      assert(SCIPpricestoreGetNVars(pricestore) == 0);
      assert(SCIPpricestoreGetNBoundResets(pricestore) == 0);
      assert(!lp->flushed || lp->solved || *lperror);

      /* put all initial constraints into the LP */
      SCIP_CALL( SCIPinitConssLP(blkmem, set, sepastore, cutpool, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
            eventfilter, cliquetable, FALSE, FALSE, &cutoff) );
      assert(cutoff == FALSE);

      mustprice = mustprice || !lp->flushed || (transprob->ncolvars != *npricedcolvars);
      *mustsepa = *mustsepa || !lp->flushed;

      /* if all pricers wanted to stop pricing, do not do another pricing round (LP value is no valid dual bound in this case) */
      if( stoppricing )
      {
         SCIPsetDebugMsg(set, "pricing: stop pricing and perform early branching\n");
         mustprice = FALSE;
         *aborted = TRUE;
      }

      /* solve LP again after resetting bounds and adding new initial constraints (with dual simplex) */
      SCIPsetDebugMsg(set, "pricing: solve LP after resetting bounds and adding new initial constraints\n");
      SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob, -1LL, FALSE, FALSE, FALSE, lperror) );
      assert(lp->flushed);
      assert(lp->solved || *lperror);

      /* remove previous primal ray, store new one if LP is unbounded */
      SCIP_CALL( updatePrimalRay(blkmem, set, stat, transprob, primal, tree, lp, *lperror) );

      /* increase pricing round counter */
      stat->npricerounds++;
      npricerounds++;

      /* display node information line */
      if( displayinfo && mustprice )
      {
         if( (SCIP_VERBLEVEL)set->disp_verblevel >= SCIP_VERBLEVEL_FULL
             || ((SCIP_VERBLEVEL)set->disp_verblevel >= SCIP_VERBLEVEL_HIGH && npricerounds % 100 == 1) )
         {
            SCIP_CALL( SCIPdispPrintLine(set, messagehdlr, stat, NULL, TRUE, TRUE) );
         }
      }

      /* if the LP is unbounded, we can stop pricing */
      mustprice = mustprice &&
         (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE
          || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT );

      /* if the lower bound is already higher than the cutoff bound, we can stop pricing */
      mustprice = mustprice && SCIPsetIsLT(set, SCIPnodeGetLowerbound(currentnode), primal->cutoffbound);
   }
   assert(lp->flushed);
   assert(lp->solved || *lperror);

   *aborted = ( (*aborted) || (*lperror) || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_NOTSOLVED
      || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_ERROR || npricerounds == maxpricerounds );

   /* set information, whether the current lp is a valid relaxation of the current problem */
   SCIPlpSetIsRelax(lp, !(*aborted));

   return SCIP_OKAY;
}

/** separates cuts of the cut pool */
static
SCIP_RETCODE cutpoolSeparate(
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_Bool             cutpoolisdelayed,   /**< is the cutpool delayed (count cuts found)? */
   SCIP_Bool             root,               /**< are we at the root node? */
   int                   actdepth,           /**< the depth of the focus node */
   SCIP_Bool*            enoughcuts,         /**< pointer to store if enough cuts were found in current separation round */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   if( (set->sepa_poolfreq == 0 && actdepth == 0)
      || (set->sepa_poolfreq > 0 && actdepth % set->sepa_poolfreq == 0) )
   {
      SCIP_RESULT result;

      SCIP_CALL( SCIPcutpoolSeparate(cutpool, blkmem, set, stat, eventqueue, eventfilter, lp, sepastore, NULL, cutpoolisdelayed, root, &result) );
      *cutoff = *cutoff || (result == SCIP_CUTOFF);
      *enoughcuts = *enoughcuts || (SCIPsepastoreGetNCuts(sepastore) >= 2 * (SCIP_Longint)SCIPsetGetSepaMaxcuts(set, root)) || (result == SCIP_NEWROUND);
   }

   return SCIP_OKAY;
}

/** solve the current LP of a node with a price-and-cut loop */
static
SCIP_RETCODE priceAndCutLoop(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_MEM*             mem,                /**< block memory pools */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_CUTPOOL*         delayedcutpool,     /**< global delayed cut pool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             fullseparation,     /**< are we in the first prop-and-cut-and-price loop? */
   SCIP_Bool*            propagateagain,     /**< pointer to store whether we want to propagate again */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            unbounded,          /**< pointer to store whether an unbounded ray was found in the LP */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved error in LP solving occured */
   SCIP_Bool*            pricingaborted      /**< pointer to store whether the pricing was aborted and the lower bound must
                                              *   not be used */
   )
{
   SCIP_NODE* focusnode;
   /* cppcheck-suppress unassignedVariable */
   SCIP_EVENT event;
   SCIP_LPSOLSTAT stalllpsolstat;
   SCIP_Real loclowerbound;
   SCIP_Real glblowerbound;
   SCIP_Real bounddist;
   SCIP_Real stalllpobjval;
   SCIP_Bool separate;
   SCIP_Bool mustprice;
   SCIP_Bool mustsepa;
   SCIP_Bool delayedsepa;
   SCIP_Bool root;
   SCIP_Bool allowlocal;
   int maxseparounds;
   int maxincseparounds;
   int nsepastallrounds;
   int maxnsepastallrounds;
   int stallnfracs;
   int actdepth;
   int npricedcolvars;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(pricestore != NULL);
   assert(sepastore != NULL);
   assert(cutpool != NULL);
   assert(delayedcutpool != NULL);
   assert(primal != NULL);
   assert(cutoff != NULL);
   assert(unbounded != NULL);
   assert(lperror != NULL);

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(focusnode != NULL);
   assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   actdepth = SCIPnodeGetDepth(focusnode);
   root = (actdepth == 0);

   /* check, if we want to separate at this node */
   loclowerbound = SCIPnodeGetLowerbound(focusnode);
   glblowerbound = SCIPtreeGetLowerbound(tree, set);
   assert(primal->cutoffbound > glblowerbound);
   bounddist = (loclowerbound - glblowerbound)/(primal->cutoffbound - glblowerbound);
   allowlocal = SCIPsetIsLE(set, bounddist, set->sepa_maxlocalbounddist);
   separate = (set->sepa_maxruns == -1 || stat->nruns <= set->sepa_maxruns);

   /* get maximal number of separation rounds */
   maxseparounds = (root ? set->sepa_maxroundsroot : set->sepa_maxrounds);
   if( maxseparounds == -1 )
      maxseparounds = INT_MAX;
   if( stat->nruns > 1 && root && set->sepa_maxroundsrootsubrun >= 0 )
      maxseparounds = MIN(maxseparounds, set->sepa_maxroundsrootsubrun);
   if( !fullseparation && set->sepa_maxaddrounds >= 0 )
      maxseparounds = MIN(maxseparounds, stat->nseparounds + set->sepa_maxaddrounds);
   maxnsepastallrounds = root ? set->sepa_maxstallroundsroot : set->sepa_maxstallrounds;
   if( maxnsepastallrounds == -1 )
      maxnsepastallrounds = INT_MAX;
   maxincseparounds = set->sepa_maxincrounds;
   if( maxincseparounds == -1 )
      maxincseparounds = INT_MAX;

   /* solve initial LP of price-and-cut loop */
   /* @todo check if LP is always already solved, because of calling solveNodeInitialLP() in solveNodeLP()? */
   SCIPsetDebugMsg(set, "node: solve LP with price and cut\n");
   SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem,  stat, eventqueue, eventfilter, transprob,
         set->lp_iterlim, FALSE, TRUE, FALSE, lperror) );
   assert(lp->flushed);
   assert(lp->solved || *lperror);

   /* remove previous primal ray, store new one if LP is unbounded */
   SCIP_CALL( updatePrimalRay(blkmem, set, stat, transprob, primal, tree, lp, *lperror) );

   /* price-and-cut loop */
   npricedcolvars = transprob->ncolvars;
   mustprice = TRUE;
   mustsepa = separate;
   delayedsepa = FALSE;
   *cutoff = FALSE;
   *unbounded = (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
   nsepastallrounds = 0;
   stalllpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   stalllpobjval = SCIP_REAL_MIN;
   stallnfracs = INT_MAX;
   lp->installing = FALSE;
   while( !(*cutoff) && !(*unbounded) && !(*lperror) && (mustprice || mustsepa || delayedsepa) )
   {
      SCIPsetDebugMsg(set, "-------- node solving loop --------\n");
      assert(lp->flushed);
      assert(lp->solved);

      /* solve the LP with pricing in new variables */
      while( mustprice && !(*lperror) )
      {
         SCIP_CALL( SCIPpriceLoop(blkmem, set, messagehdlr, stat, transprob, origprob, primal, tree, reopt, lp,
               pricestore, sepastore, cutpool, branchcand, eventqueue, eventfilter, cliquetable, root, root, -1, &npricedcolvars,
               &mustsepa, lperror, pricingaborted) );

         mustprice = FALSE;

         assert(lp->flushed);
         assert(lp->solved || *lperror);

         /* update lower bound w.r.t. the LP solution */
         if( !(*lperror) && !(*pricingaborted) && SCIPlpIsRelax(lp) )
         {
            SCIP_CALL( SCIPnodeUpdateLowerboundLP(focusnode, set, stat, tree, transprob, origprob, lp) );
            SCIPsetDebugMsg(set, " -> new lower bound: %g (LP status: %d, LP obj: %g)\n",
               SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set, transprob));

            /* update node estimate */
            SCIP_CALL( updateEstimate(set, stat, tree, lp, branchcand) );

            if( root && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
               SCIPprobUpdateBestRootSol(transprob, set, stat, lp);
         }
         else
         {
            SCIPsetDebugMsg(set, " -> error solving LP or pricing aborted. keeping old bound: %g\n", SCIPnodeGetLowerbound(focusnode));
         }

         /* display node information line for root node */
         if( root && (SCIP_VERBLEVEL)set->disp_verblevel >= SCIP_VERBLEVEL_HIGH )
         {
            SCIP_CALL( SCIPdispPrintLine(set, messagehdlr, stat, NULL, TRUE, TRUE) );
         }

         if( !(*lperror) )
         {
            /* call propagators that are applicable during LP solving loop only if the node is not cut off */
            if( SCIPsetIsLT(set, SCIPnodeGetLowerbound(focusnode), primal->cutoffbound) )
            {
               SCIP_Longint oldnboundchgs;
               SCIP_Longint oldninitconssadded;
               SCIP_Bool postpone;

               oldnboundchgs = stat->nboundchgs;
               oldninitconssadded = stat->ninitconssadded;

               SCIPsetDebugMsg(set, " -> LP solved: call propagators that are applicable during LP solving loop\n");

               SCIP_CALL( propagateDomains(blkmem, set, stat, primal, tree, SCIPtreeGetCurrentDepth(tree), 0, FALSE,
                     SCIP_PROPTIMING_DURINGLPLOOP, cutoff, &postpone) );
               assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);
               assert(!postpone);

               if( stat->ninitconssadded != oldninitconssadded )
               {
                  SCIPsetDebugMsg(set, "new initial constraints added during propagation: old=%" SCIP_LONGINT_FORMAT ", new=%" SCIP_LONGINT_FORMAT "\n", oldninitconssadded, stat->ninitconssadded);

                  SCIP_CALL( SCIPinitConssLP(blkmem, set, sepastore, cutpool, stat, transprob, origprob, tree, reopt, lp,
                        branchcand, eventqueue, eventfilter, cliquetable, FALSE, FALSE, cutoff) );
               }

               if( !(*cutoff) && !(*unbounded) )
               {
                  /* if we found something, solve LP again */
                  if( !lp->flushed )
                  {
                     SCIPsetDebugMsg(set, "    -> found reduction: resolve LP\n");

                     /* in the root node, remove redundant rows permanently from the LP */
                     if( root )
                     {
                        SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );
                        SCIP_CALL( SCIPlpRemoveRedundantRows(lp, blkmem, set, stat, eventqueue, eventfilter) );
                     }

                     /* resolve LP */
                     SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob,
                           set->lp_iterlim, FALSE, TRUE, FALSE, lperror) );
                     assert(lp->flushed);
                     assert(lp->solved || *lperror);

                     /* remove previous primal ray, store new one if LP is unbounded */
                     SCIP_CALL( updatePrimalRay(blkmem, set, stat, transprob, primal, tree, lp, *lperror) );

                     mustprice = TRUE;
                     *propagateagain = TRUE;
                  }
                  /* propagation might have changed the best bound of loose variables, thereby changing the loose objective
                   * value which is added to the LP value; because of the loose status, the LP might not be reoptimized,
                   * but the lower bound of the node needs to be updated
                   */
                  else if( stat->nboundchgs > oldnboundchgs )
                  {
                     *propagateagain = TRUE;

                     if( lp->solved && SCIPprobAllColsInLP(transprob, set, lp) && SCIPlpIsRelax(lp) )
                     {
                        assert(lp->flushed);
                        assert(lp->solved);

                        SCIP_CALL( SCIPnodeUpdateLowerboundLP(focusnode, set, stat, tree, transprob, origprob, lp) );
                        SCIPsetDebugMsg(set, " -> new lower bound: %g (LP status: %d, LP obj: %g)\n",
                           SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set, transprob));

                        /* update node estimate */
                        SCIP_CALL( updateEstimate(set, stat, tree, lp, branchcand) );

                        if( root && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
                           SCIPprobUpdateBestRootSol(transprob, set, stat, lp);
                     }
                  }
               }
            }
         }

         /* call primal heuristics that are applicable during node LP solving loop */
         if( !*cutoff && !(*unbounded) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIP_Bool foundsol;

            SCIP_CALL( SCIPprimalHeuristics(set, stat, transprob, primal, tree, lp, NULL, SCIP_HEURTIMING_DURINGLPLOOP,
                  FALSE, &foundsol, unbounded) );
            assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

            *lperror = *lperror || lp->resolvelperror;
         }
      }
      assert(lp->flushed || *cutoff || *unbounded);
      assert(lp->solved || *lperror || *cutoff || *unbounded);

      /* check, if we exceeded the separation round limit */
      mustsepa = mustsepa
         && stat->nseparounds < maxseparounds
         && stat->nincseparounds < maxincseparounds
         && nsepastallrounds < maxnsepastallrounds
         && !(*cutoff);

      /* if separators were delayed, we want to apply a final separation round with the delayed separators */
      delayedsepa = delayedsepa && !mustsepa && !(*cutoff); /* if regular separation applies, we ignore delayed separators */
      mustsepa = mustsepa || delayedsepa;

      if( mustsepa )
      {
         /* if the LP is infeasible, unbounded, exceeded the objective limit or a global performance limit was reached,
          * we don't need to separate cuts
          * (the global limits are only checked at the root node in order to not query system time too often)
          */
         if( !separate || (*cutoff) || (*unbounded)
             || (SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY)
             || SCIPsetIsGE(set, SCIPnodeGetLowerbound(focusnode), primal->cutoffbound)
             || (root && SCIPsolveIsStopped(set, stat, FALSE)) )
         {
            mustsepa = FALSE;
            delayedsepa = FALSE;
         }
         else
            assert(!(*lperror));
      }

      /* separation (needs not to be done completely, because we just want to increase the lower bound) */
      if( mustsepa )
      {
         SCIP_Longint olddomchgcount;
         SCIP_Longint oldninitconssadded;
         int oldnlprows;
         SCIP_Bool enoughcuts;

         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

         olddomchgcount = stat->domchgcount;
         oldninitconssadded = stat->ninitconssadded;

         mustsepa = FALSE;
         enoughcuts = (SCIPsetGetSepaMaxcuts(set, root) == 0);
         oldnlprows = lp->nrows;

         /* global cut pool separation */
         if( !enoughcuts && !delayedsepa )
         {
            SCIP_CALL( cutpoolSeparate(cutpool, blkmem, set, stat, eventqueue, eventfilter, lp, sepastore, FALSE, root,
                  actdepth, &enoughcuts, cutoff) );

            if( *cutoff )
            {
               SCIPsetDebugMsg(set, " -> global cut pool detected cutoff\n");
            }
         }
         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

         /* constraint separation */
         SCIPsetDebugMsg(set, "constraint separation\n");

         /* separate constraints and LP */
         if( !(*cutoff) && !(*lperror) && !enoughcuts && lp->solved )
         {
            /* apply a separation round */
            SCIP_CALL( separationRoundLP(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, transprob, primal, tree,
                  lp, sepastore, actdepth, bounddist, allowlocal, delayedsepa,
                  &delayedsepa, &enoughcuts, cutoff, lperror, &mustsepa, &mustprice) );
            assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

            /* if we are close to the stall round limit, also call the delayed separators */
            if( !(*cutoff) && !(*lperror) && !enoughcuts && lp->solved
               && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY)
               && nsepastallrounds >= maxnsepastallrounds-1 && delayedsepa )
            {
               SCIP_CALL( separationRoundLP(blkmem, set, messagehdlr, stat, eventqueue, eventfilter, transprob, primal,
                     tree, lp, sepastore, actdepth, bounddist, allowlocal, delayedsepa,
                     &delayedsepa, &enoughcuts, cutoff, lperror, &mustsepa, &mustprice) );
               assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);
            }
         }

         /* call global cut pool separation again since separators may add cuts to the pool instead of the sepastore */
         if( !(*cutoff) && !(*lperror) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && !enoughcuts && lp->solved )
         {
            SCIP_CALL( cutpoolSeparate(cutpool, blkmem, set, stat, eventqueue, eventfilter, lp, sepastore, FALSE, root,
                  actdepth, &enoughcuts, cutoff) );

            if( *cutoff )
            {
               SCIPsetDebugMsg(set, " -> global cut pool detected cutoff\n");
            }
         }

         /* delayed global cut pool separation */
         if( !(*cutoff) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && SCIPsepastoreGetNCuts(sepastore) == 0 )
         {
            assert( !(*lperror) );

            SCIP_CALL( cutpoolSeparate(delayedcutpool, blkmem, set, stat, eventqueue, eventfilter, lp, sepastore, TRUE,
                  root, actdepth, &enoughcuts, cutoff) );

            if( *cutoff )
            {
               SCIPsetDebugMsg(set, " -> delayed global cut pool detected cutoff\n");
            }
            assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL);
            assert(lp->flushed);
            assert(lp->solved);
         }

         assert(*cutoff || *lperror || SCIPlpIsSolved(lp));
         assert(!SCIPlpIsSolved(lp)
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_ITERLIMIT
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_TIMELIMIT);

         if( *cutoff || *lperror
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT
            || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_ITERLIMIT  || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_TIMELIMIT )
         {
            /* the found cuts are of no use, because the node is infeasible anyway (or we have an error in the LP) */
            SCIP_CALL( SCIPsepastoreClearCuts(sepastore, blkmem, set, eventqueue, eventfilter, lp) );
         }
         else
         {
            /* apply found cuts */
            SCIP_CALL( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp,
                  branchcand, eventqueue, eventfilter, cliquetable, root, SCIP_EFFICIACYCHOICE_LP, cutoff) );

            if( !(*cutoff) )
            {
               mustprice = mustprice || !lp->flushed || (transprob->ncolvars != npricedcolvars);
               mustsepa = mustsepa || !lp->flushed;

               /* if a new bound change (e.g. a cut with only one column) was found, propagate domains again */
               if( stat->domchgcount != olddomchgcount )
               {
                  SCIPsetDebugMsg(set, " -> separation changed bound: propagate again\n");

                  *propagateagain = TRUE;

                  /* in the root node, remove redundant rows permanently from the LP */
                  if( root )
                  {
                     SCIP_CALL( SCIPlpFlush(lp, blkmem, set, eventqueue) );
                     SCIP_CALL( SCIPlpRemoveRedundantRows(lp, blkmem, set, stat, eventqueue, eventfilter) );
                  }
               }

               if( stat->ninitconssadded != oldninitconssadded )
               {
                  SCIPsetDebugMsg(set, "new initial constraints added during propagation: old=%" SCIP_LONGINT_FORMAT ", new=%" SCIP_LONGINT_FORMAT "\n",
                        oldninitconssadded, stat->ninitconssadded);

                  SCIP_CALL( SCIPinitConssLP(blkmem, set, sepastore, cutpool, stat, transprob, origprob, tree, reopt, lp,
                        branchcand, eventqueue, eventfilter, cliquetable, FALSE, FALSE, cutoff) );
               }

               if( !(*cutoff) )
               {
                  SCIP_Real lpobjval;

                  /* solve LP (with dual simplex) */
                  SCIPsetDebugMsg(set, "separation: solve LP\n");
                  SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob,
                        set->lp_iterlim, FALSE, TRUE, FALSE, lperror) );
                  assert(lp->flushed);
                  assert(lp->solved || *lperror);

                  /* remove previous primal ray, store new one if LP is unbounded */
                  SCIP_CALL( updatePrimalRay(blkmem, set, stat, transprob, primal, tree, lp, *lperror) );

                  if( !(*lperror) )
                  {
                     SCIP_Bool stalling;

                     /* propagation might have changed the best bound of loose variables, thereby changing the loose objective value
                      * which is added to the LP value; because of the loose status, the LP might not be reoptimized, but the lower
                      * bound of the node needs to be updated
                      */
                     if( stat->domchgcount != olddomchgcount && (!mustprice || mustsepa) && !(*cutoff)
                        && SCIPprobAllColsInLP(transprob, set, lp) && SCIPlpIsRelax(lp) )
                     {
                        SCIP_CALL( SCIPnodeUpdateLowerboundLP(focusnode, set, stat, tree, transprob, origprob, lp) );
                        SCIPsetDebugMsg(set, " -> new lower bound: %g (LP status: %d, LP obj: %g)\n",
                           SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set, transprob));

                        /* update node estimate */
                        SCIP_CALL( updateEstimate(set, stat, tree, lp, branchcand) );

                        if( root && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
                           SCIPprobUpdateBestRootSol(transprob, set, stat, lp);
                     }

                     /* check if we are stalling
                      * If we have an LP solution, then we are stalling if
                      *   we had an LP solution before and
                      *   the LP value did not improve and
                      *   the number of fractional variables did not decrease.
                      * If we do not have an LP solution, then we are stalling if the solution status of the LP did not change.
                      */
                     if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
                     {
                        SCIP_Real objreldiff;
                        int nfracs;

                        SCIP_CALL( SCIPbranchcandGetLPCands(branchcand, set, stat, lp, NULL, NULL, NULL, &nfracs, NULL,
                              NULL) );
                        lpobjval = SCIPlpGetObjval(lp, set, transprob);

                        objreldiff = SCIPrelDiff(lpobjval, stalllpobjval);
                        SCIPsetDebugMsg(set, " -> LP bound moved from %g to %g (reldiff: %g)\n",
                           stalllpobjval, lpobjval, objreldiff);

                        stalling = (stalllpsolstat == SCIP_LPSOLSTAT_OPTIMAL &&
                            objreldiff <= 1e-04 &&
                            nfracs >= (0.9 - 0.1 * nsepastallrounds) * stallnfracs);

                        stalllpobjval = lpobjval;
                        stallnfracs = nfracs;
                     }
                     else
                     {
                        stalling = (stalllpsolstat == SCIPlpGetSolstat(lp));
                     }

                     if( !stalling )
                     {
                        nsepastallrounds = 0;
                        lp->installing = FALSE;
                     }
                     else
                     {
                        nsepastallrounds++;
                     }
                     stalllpsolstat = SCIPlpGetSolstat(lp);

                     /* tell LP that we are (close to) stalling */
                     if( nsepastallrounds >= maxnsepastallrounds-2 )
                        lp->installing = TRUE;
                     SCIPsetDebugMsg(set, " -> nsepastallrounds=%d/%d\n", nsepastallrounds, maxnsepastallrounds);
                  }
               }
            }
         }
         assert(*cutoff || *lperror || (lp->flushed && lp->solved)); /* cutoff: LP may be unsolved due to bound changes */

         SCIPsetDebugMsg(set, "separation round %d/%d finished (%d/%d stall rounds): mustprice=%u, mustsepa=%u, delayedsepa=%u, propagateagain=%u\n",
            stat->nseparounds, maxseparounds, nsepastallrounds, maxnsepastallrounds, mustprice, mustsepa, delayedsepa, *propagateagain);

         /* increase separation round counter */
         stat->nseparounds++;

         /* if size of lp relaxation increased also count this round separately */
         if( lp->nrows > oldnlprows )
            ++stat->nincseparounds;
         else if( lp->nrows < oldnlprows )
            stat->nincseparounds = 0;
      }
   }

   if( root && nsepastallrounds >= maxnsepastallrounds )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
         "Truncate separation round because of stalling (%d stall rounds).\n", maxnsepastallrounds);
   }

   if( !*lperror )
   {
      /* update pseudo cost values for continuous variables, if it should be delayed */
      SCIP_CALL( updatePseudocost(set, stat, transprob, tree, lp, FALSE, set->branch_delaypscost) );
   }

   /* update lower bound w.r.t. the LP solution */
   if( *cutoff )
   {
      SCIPnodeUpdateLowerbound(focusnode, stat, set, tree, transprob, origprob, SCIPsetInfinity(set));
   }
   else if( !(*lperror) )
   {
      assert(lp->flushed);
      assert(lp->solved);

      if( SCIPlpIsRelax(lp) )
      {
         SCIP_CALL( SCIPnodeUpdateLowerboundLP(focusnode, set, stat, tree, transprob, origprob, lp) );
      }

      /* update node estimate */
      SCIP_CALL( updateEstimate(set, stat, tree, lp, branchcand) );

      /* issue LPSOLVED event */
      if( SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_ITERLIMIT && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_TIMELIMIT )
      {
         SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_LPSOLVED) );
         SCIP_CALL( SCIPeventChgNode(&event, focusnode) );
         SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
      }

      /* if the LP is a relaxation and we are not solving exactly, then we may analyze an infeasible or bound exceeding
       * LP (not necessary in the root node) and cut off the current node
       */
      if( !set->misc_exactsolve && !root && SCIPlpIsRelax(lp) && SCIPprobAllColsInLP(transprob, set, lp)
         && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT) )
      {
         SCIP_CALL( SCIPconflictAnalyzeLP(conflict, conflictstore, blkmem, set, stat, transprob, origprob, tree, reopt,
               lp, branchcand, eventqueue, cliquetable, NULL) );
         *cutoff = TRUE;
      }
   }
   /* check for unboundedness */
   if( !(*lperror) )
   {
      *unbounded = (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
      /* assert(!(*unbounded) || root); */ /* unboundedness can only happen in the root node; no, of course it can also happens in the tree if a branching did not help to resolve unboundedness */
   }

   lp->installing = FALSE;

   SCIPsetDebugMsg(set, " -> final lower bound: %g (LP status: %d, LP obj: %g)\n",
      SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp),
      (*cutoff || *unbounded) ? SCIPsetInfinity(set) : *lperror ? -SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set, transprob));

   return SCIP_OKAY;
}

/** updates the current lower bound with the pseudo objective value, cuts off node by bounding, and applies conflict
 *  analysis if the pseudo objective lead to the cutoff
 */
static
SCIP_RETCODE applyBounding(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if the node can be cut off */
   )
{
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(primal != NULL);
   assert(cutoff != NULL);

   if( !(*cutoff) )
   {
      SCIP_NODE* focusnode;
      SCIP_Real pseudoobjval;

      /* get current focus node */
      focusnode = SCIPtreeGetFocusNode(tree);

      /* update lower bound w.r.t. the pseudo solution */
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set, transprob);
      SCIPnodeUpdateLowerbound(focusnode, stat, set, tree, transprob, origprob, pseudoobjval);
      SCIPsetDebugMsg(set, " -> lower bound: %g [%g] (pseudoobj: %g [%g]), cutoff bound: %g [%g]\n",
         SCIPnodeGetLowerbound(focusnode), SCIPprobExternObjval(transprob, origprob, set, SCIPnodeGetLowerbound(focusnode)) + SCIPgetOrigObjoffset(set->scip),
         pseudoobjval, SCIPprobExternObjval(transprob, origprob, set, pseudoobjval) + SCIPgetOrigObjoffset(set->scip),
         primal->cutoffbound, SCIPprobExternObjval(transprob, origprob, set, primal->cutoffbound) + SCIPgetOrigObjoffset(set->scip));

      /* check for infeasible node by bounding */
      if( (set->misc_exactsolve && SCIPnodeGetLowerbound(focusnode) >= primal->cutoffbound)
         || (!set->misc_exactsolve && SCIPsetIsGE(set, SCIPnodeGetLowerbound(focusnode), primal->cutoffbound)) )
      {
         SCIPsetDebugMsg(set, "node is cut off by bounding (lower=%g, upper=%g)\n",
            SCIPnodeGetLowerbound(focusnode), primal->cutoffbound);
         SCIPnodeUpdateLowerbound(focusnode, stat, set, tree, transprob, origprob, SCIPsetInfinity(set));
         *cutoff = TRUE;

         /* call pseudo conflict analysis, if the node is cut off due to the pseudo objective value */
         if( pseudoobjval >= primal->cutoffbound && !SCIPsetIsInfinity(set, primal->cutoffbound) && !SCIPsetIsInfinity(set, -pseudoobjval) )
         {
            SCIP_CALL( SCIPconflictAnalyzePseudo(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable, NULL) );
         }
      }
   }

   return SCIP_OKAY;
}

/** marks all relaxators to be unsolved */
static
void markRelaxsUnsolved(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   )
{
   int r;

   assert(set != NULL);
   assert(relaxation != NULL);

   SCIPrelaxationSetSolValid(relaxation, FALSE, FALSE);

   for( r = 0; r < set->nrelaxs; ++r )
      SCIPrelaxMarkUnsolved(set->relaxs[r]);
}

/** solves the current node's LP in a price-and-cut loop */
static
SCIP_RETCODE solveNodeLP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_MEM*             mem,                /**< block memory pools */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_RELAXATION*      relaxation,         /**< relaxators */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_CUTPOOL*         delayedcutpool,     /**< global delayed cut pool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             initiallpsolved,    /**< was the initial LP already solved? */
   SCIP_Bool             fullseparation,     /**< are we in the first prop-and-cut-and-price loop? */
   SCIP_Bool             newinitconss,       /**< do we have to add new initial constraints? */
   SCIP_Bool*            propagateagain,     /**< pointer to store whether we want to propagate again */
   SCIP_Bool*            solverelaxagain,    /**< pointer to store TRUE, if the external relaxators should be called again */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            unbounded,          /**< pointer to store TRUE, if an unbounded ray was found in the LP */
   SCIP_Bool*            lperror,            /**< pointer to store TRUE, if an unresolved error in LP solving occured */
   SCIP_Bool*            pricingaborted      /**< pointer to store TRUE, if the pricing was aborted and the lower bound
                                              *   must not be used */
   )
{
   SCIP_Longint nlpiterations;
   SCIP_Longint nlps;

   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeHasFocusNodeLP(tree));
   assert(cutoff != NULL);
   assert(unbounded != NULL);
   assert(lperror != NULL);
   assert(*cutoff == FALSE);
   assert(*unbounded == FALSE);
   assert(*lperror == FALSE);

   nlps = stat->nlps;
   nlpiterations = stat->nlpiterations;

   if( !initiallpsolved )
   {
      /* load and solve the initial LP of the node */
      SCIP_CALL( solveNodeInitialLP(blkmem, set, messagehdlr, stat, transprob, origprob, primal, tree, reopt, lp,
            pricestore, sepastore, cutpool, branchcand, eventfilter, eventqueue, cliquetable, newinitconss, cutoff, lperror) );

      assert(*cutoff || *lperror || (lp->flushed && lp->solved));
      SCIPsetDebugMsg(set, "price-and-cut-loop: initial LP status: %d, LP obj: %g\n",
         SCIPlpGetSolstat(lp),
         *cutoff ? SCIPsetInfinity(set) : *lperror ? -SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set, transprob));

      /* update initial LP iteration counter */
      stat->ninitlps += stat->nlps - nlps;
      stat->ninitlpiterations += stat->nlpiterations - nlpiterations;

      /* In the root node, we try if the initial LP solution is feasible to avoid expensive setup of data structures in
       * separators; in case the root LP is aborted, e.g., by hitting the time limit, we do not check the LP solution
       * since the corresponding data structures have not been updated.
       */
      if( SCIPtreeGetCurrentDepth(tree) == 0 && !(*cutoff) && !(*lperror)
         && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY)
         && !SCIPsolveIsStopped(set, stat, FALSE) )
      {
         SCIP_Bool checklprows;
         SCIP_Bool stored;
         SCIP_SOL* sol;
         SCIP_Longint oldnbestsolsfound = primal->nbestsolsfound;

         SCIP_CALL( SCIPsolCreateLPSol(&sol, blkmem, set, stat, transprob, primal, tree, lp, NULL) );

         if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
            checklprows = FALSE;
         else
            checklprows = TRUE;

#ifndef NDEBUG
         /* in the debug mode we want to explicitly check if the solution is feasible if it was stored */
         SCIP_CALL( SCIPprimalTrySol(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
               eventqueue, eventfilter, sol, FALSE, FALSE, TRUE, TRUE, checklprows, &stored) );

         if( stored )
         {
            SCIP_Bool feasible;

            SCIP_CALL( SCIPsolCheck(sol, set, messagehdlr, blkmem, stat, transprob, FALSE, FALSE, TRUE, TRUE,
                  checklprows, &feasible) );
            assert(feasible);
         }

         SCIP_CALL( SCIPsolFree(&sol, blkmem, primal) );
#else
         SCIP_CALL( SCIPprimalTrySolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
               eventqueue, eventfilter, &sol, FALSE, FALSE, TRUE, TRUE, checklprows, &stored) );
#endif
         if( stored )
         {
            stat->nlpsolsfound++;

            if( primal->nbestsolsfound != oldnbestsolsfound )
            {
               stat->nlpbestsolsfound++;
               SCIPstoreSolutionGap(set->scip);
            }

            if( set->reopt_enable )
            {
               assert(reopt != NULL);
               SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, tree->focusnode, SCIP_EVENTTYPE_NODEFEASIBLE, lp,
                     SCIPlpGetSolstat(lp), tree->root == tree->focusnode, TRUE, tree->focusnode->lowerbound,
                     tree->effectiverootdepth) );
            }
         }

         if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
            *unbounded = TRUE;
      }
   }
   else
   {
      SCIP_CALL( SCIPinitConssLP(blkmem, set, sepastore, cutpool, stat, transprob,
            origprob, tree, reopt, lp, branchcand, eventqueue, eventfilter, cliquetable, FALSE, FALSE,
            cutoff) );
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   /* check for infeasible node by bounding */
   SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue, conflict, cliquetable, cutoff) );
#ifdef SCIP_DEBUG
   if( *cutoff )
   {
      if( SCIPtreeGetCurrentDepth(tree) == 0 )
      {
         SCIPsetDebugMsg(set, "solution cuts off root node, stop solution process\n");
      }
      else
      {
         SCIPsetDebugMsg(set, "solution cuts off node\n");
      }
   }
#endif

   if( !(*cutoff) && !(*lperror) )
   {
      SCIP_Longint oldninitconssadded;
      SCIP_Longint oldnboundchgs;
      SCIP_Longint olddomchgcount;
      int oldnpricedvars;
      int oldncutsapplied;

      oldnpricedvars = transprob->ncolvars;
      oldninitconssadded = stat->ninitconssadded;
      oldncutsapplied = SCIPsepastoreGetNCutsApplied(sepastore);
      oldnboundchgs = stat->nboundchgs;
      olddomchgcount = stat->domchgcount;

      /* solve the LP with price-and-cut*/
      SCIP_CALL( priceAndCutLoop(blkmem, set, messagehdlr, stat, mem, transprob, origprob,  primal, tree, reopt, lp,
            pricestore, sepastore, cutpool, delayedcutpool, branchcand, conflict, conflictstore, eventfilter,
            eventqueue, cliquetable, fullseparation, propagateagain, cutoff, unbounded, lperror, pricingaborted) );

      /* check if the problem was changed and the relaxation needs to be resolved */
      if( (transprob->ncolvars != oldnpricedvars) || (stat->ninitconssadded != oldninitconssadded) ||
            (SCIPsepastoreGetNCutsApplied(sepastore) != oldncutsapplied) || (stat->nboundchgs != oldnboundchgs) ||
            (stat->domchgcount != olddomchgcount) )
      {
         *solverelaxagain = TRUE;
         markRelaxsUnsolved(set, relaxation);
      }

   }
   assert(*cutoff || *lperror || (lp->flushed && lp->solved));

   /* if there is no LP error, then *unbounded should be TRUE, iff the LP solution status is unboundedray */
   assert(*lperror || ((SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY) == *unbounded));

   /* If pricing was aborted while solving the LP of the node and the node cannot be cut off due to the lower bound computed by the pricer,
   *  the solving of the LP might be stopped due to the objective limit, but the node may not be cut off, since the LP objective
   *  is not a feasible lower bound for the solutions in the current subtree.
   *  In this case, the LP has to be solved to optimality by temporarily removing the cutoff bound.
   */
   if( (*pricingaborted) && (SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_ITERLIMIT)
      && !(*cutoff) )
   {
      SCIP_Real tmpcutoff;

      /* temporarily disable cutoffbound, which also disables the objective limit */
      tmpcutoff = lp->cutoffbound;
      lp->cutoffbound = SCIPlpiInfinity(SCIPlpGetLPI(lp));

      lp->solved = FALSE;
      SCIP_CALL( SCIPlpSolveAndEval(lp, set, messagehdlr, blkmem, stat, eventqueue, eventfilter, transprob, -1LL, FALSE, FALSE, FALSE, lperror) );

      /* reinstall old cutoff bound */
      lp->cutoffbound = tmpcutoff;

      SCIPsetDebugMsg(set, "re-optimized LP without cutoff bound: LP status: %d, LP obj: %g\n",
         SCIPlpGetSolstat(lp), *lperror ? -SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set, transprob));

      /* lp solstat should not be objlimit, since the cutoff bound was removed temporarily */
      assert(SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OBJLIMIT);
      /* lp solstat should not be unboundedray, since the lp was dual feasible */
      assert(SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY);
      /* there should be no primal ray, since the lp was dual feasible */
      assert(primal->primalray == NULL);
      if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE )
      {
         *cutoff = TRUE;
      }
   }

   assert(!(*pricingaborted) || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL /* cppcheck-suppress assertWithSideEffect */
      || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_NOTSOLVED || SCIPsolveIsStopped(set, stat, FALSE) || (*cutoff));

   assert(*cutoff || *lperror || (lp->flushed && lp->solved));

   /* update node's LP iteration counter */
   stat->nnodelps += stat->nlps - nlps;
   stat->nnodelpiterations += stat->nlpiterations - nlpiterations;

   /* update number of root node LPs and iterations if the root node was processed */
   if( SCIPnodeGetDepth(tree->focusnode) == 0 )
   {
      stat->nrootlps += stat->nlps - nlps;
      stat->nrootlpiterations += stat->nlpiterations - nlpiterations;
   }

   return SCIP_OKAY;
}

/** calls relaxators */
static
SCIP_RETCODE solveNodeRelax(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_RELAXATION*      relaxation,         /**< relaxators */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             beforelp,           /**< should the relaxators with non-negative or negative priority be called? */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   SCIP_Bool*            solvelpagain,       /**< pointer to store TRUE, if the node's LP has to be solved again */
   SCIP_Bool*            solverelaxagain,    /**< pointer to store TRUE, if the external relaxators should be called
                                              *   again */
   SCIP_Bool*            relaxcalled         /**< pointer to store TRUE, if at least one relaxator was called (unmodified
                                              *   otherwise) */
   )
{
   SCIP_RESULT result;
   SCIP_Real lowerbound;
   int r;

   assert(set != NULL);
   assert(relaxation != NULL);
   assert(cutoff != NULL);
   assert(solvelpagain != NULL);
   assert(propagateagain != NULL);
   assert(solverelaxagain != NULL);
   assert(relaxcalled != NULL);
   assert(!(*cutoff));

   /* sort by priority */
   SCIPsetSortRelaxs(set);

   for( r = 0; r < set->nrelaxs && !(*cutoff); ++r )
   {
      if( beforelp != (SCIPrelaxGetPriority(set->relaxs[r]) >= 0) )
         continue;

      *relaxcalled = TRUE;

      lowerbound = -SCIPsetInfinity(set);

      SCIP_CALL( SCIPrelaxExec(set->relaxs[r], set, stat, depth, &lowerbound, &result) );

      switch( result )
      {
      case SCIP_CUTOFF:
         *cutoff = TRUE;
         SCIPsetDebugMsg(set, " -> relaxator <%s> detected cutoff\n", SCIPrelaxGetName(set->relaxs[r]));
         /* @todo does it make sense to proceed if the node is proven to be infeasible? */
         return SCIP_OKAY;

      case SCIP_CONSADDED:
         *solvelpagain = TRUE;   /* the separation for new constraints should be called */
         *propagateagain = TRUE; /* the propagation for new constraints should be called */
         break;

      case SCIP_REDUCEDDOM:
         *solvelpagain = TRUE;
         *propagateagain = TRUE;
         break;

      case SCIP_SEPARATED:
         *solvelpagain = TRUE;
         break;

      case SCIP_SUSPENDED:
         *solverelaxagain = TRUE;
         break;

      case SCIP_SUCCESS:
      case SCIP_DIDNOTRUN:
         break;

      default:
         SCIPerrorMessage("invalid result code <%d> of relaxator <%s>\n", result, SCIPrelaxGetName(set->relaxs[r]));
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/

      if( result != SCIP_CUTOFF && result != SCIP_DIDNOTRUN && result != SCIP_SUSPENDED )
      {
         SCIP_NODE* focusnode;

         focusnode = SCIPtreeGetFocusNode(tree);
         assert(focusnode != NULL);
         assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);

         /* update lower bound w.r.t. the lower bound given by the relaxator */
         SCIPnodeUpdateLowerbound(focusnode, stat, set, tree, transprob, origprob, lowerbound);
         SCIPsetDebugMsg(set, " -> new lower bound given by relaxator %s: %g\n",
            SCIPrelaxGetName(set->relaxs[r]), lowerbound);
      }
   }

   return SCIP_OKAY;
}

/** enforces constraints by branching, separation, or domain reduction */
static
SCIP_RETCODE enforceConstraints(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_Bool*            branched,           /**< pointer to store whether a branching was created */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the LP/pseudo solution is infeasible */
   SCIP_Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   SCIP_Bool*            solvelpagain,       /**< pointer to store TRUE, if the node's LP has to be solved again */
   SCIP_Bool*            solverelaxagain,    /**< pointer to store TRUE, if the external relaxators should be called again */
   SCIP_Bool             forced              /**< should enforcement of pseudo solution be forced? */
   )
{
   SCIP_RESULT result;
   SCIP_SOL* relaxsol = NULL;
   SCIP_Real pseudoobjval;
   SCIP_Bool resolved;
   SCIP_Bool objinfeasible;
   SCIP_Bool enforcerelaxsol;
   int h;

   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(SCIPtreeGetFocusNode(tree) != NULL);
   assert(branched != NULL);
   assert(cutoff != NULL);
   assert(infeasible != NULL);
   assert(propagateagain != NULL);
   assert(solvelpagain != NULL);
   assert(solverelaxagain != NULL);
   assert(!(*cutoff));
   assert(!(*propagateagain));
   assert(!(*solvelpagain));
   assert(!(*solverelaxagain));

   *branched = FALSE;
   /**@todo avoid checking the same pseudosolution twice */


   /* enforce (best) relaxation solution if the LP has a worse objective value */
   enforcerelaxsol = SCIPrelaxationIsSolValid(relaxation) && SCIPrelaxationIsLpIncludedForSol(relaxation) && (!SCIPtreeHasFocusNodeLP(tree)
         || SCIPsetIsGT(set, SCIPrelaxationGetSolObj(relaxation), SCIPlpGetObjval(lp, set, prob)));

   /* check if all constraint handlers implement the enforelax callback, otherwise enforce the LP solution */
   for( h = 0; h < set->nconshdlrs && enforcerelaxsol; ++h )
   {
      if( set->conshdlrs_enfo[h]->consenforelax == NULL && ((! set->conshdlrs_enfo[h]->needscons) ||
            (set->conshdlrs_enfo[h]->nconss > 0)) )
      {
         SCIP_VERBLEVEL verblevel;

         enforcerelaxsol = FALSE;

         verblevel = SCIP_VERBLEVEL_FULL;

         if( !stat->disableenforelaxmsg && set->disp_verblevel == SCIP_VERBLEVEL_HIGH )
         {
            verblevel = SCIP_VERBLEVEL_HIGH;

            /* remember that the disable relaxation enforcement message was posted and only post it again if the
             * verblevel is SCIP_VERBLEVEL_FULL
             */
            stat->disableenforelaxmsg = TRUE;
         }
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, verblevel, "Disable enforcement of relaxation solutions"
               " since constraint handler %s does not implement enforelax-callback\n",
               SCIPconshdlrGetName(set->conshdlrs_enfo[h]));
      }
   }

   /* enforce constraints by branching, applying additional cutting planes (if LP is being processed),
    * introducing new constraints, or tighten the domains
    */
#ifndef SCIP_NDEBUG
   if( enforcerelaxsol )
   {
      SCIPsetDebugMsg(set, "enforcing constraints on relaxation solution\n");
   }
   else
   {
      SCIPsetDebugMsg(set, "enforcing constraints on %s solution\n", SCIPtreeHasFocusNodeLP(tree) ? "LP" : "pseudo");
   }
#endif

   /* check, if the solution is infeasible anyway due to it's objective value */
   if( SCIPtreeHasFocusNodeLP(tree) || enforcerelaxsol )
      objinfeasible = FALSE;
   else
   {
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set, prob);
      objinfeasible = SCIPsetIsFeasLT(set, pseudoobjval, SCIPnodeGetLowerbound(SCIPtreeGetFocusNode(tree)));
   }

   /* during constraint enforcement, generated cuts should enter the LP in any case; otherwise, a constraint handler
    * would fail to enforce its constraints if it relies on the modification of the LP relaxation
    */
   SCIPsepastoreStartForceCuts(sepastore);

   /* enforce constraints until a handler resolved an infeasibility with cutting off the node, branching,
    * reducing a domain, or separating a cut
    * if a constraint handler introduced new constraints to enforce his constraints, the newly added constraints
    * have to be enforced themselves
    */
   resolved = FALSE;

   /* in case the relaxation solution should be enforced, we need to create the corresponding solution for the enforelax callbacks */
   if( enforcerelaxsol )
   {
      SCIP_CALL( SCIPsolCreateRelaxSol(&relaxsol, blkmem, set, stat, primal, tree, relaxation, NULL) );
   }

   for( h = 0; h < set->nconshdlrs && !resolved; ++h )
   {
      assert(SCIPsepastoreGetNCuts(sepastore) == 0); /* otherwise, the LP should have been resolved first */

      /* enforce LP, pseudo, or relaxation solution */
      if( enforcerelaxsol )
      {
         SCIPsetDebugMsg(set, "enforce relaxation solution with value %g\n", SCIPrelaxationGetSolObj(relaxation));

         SCIP_CALL( SCIPconshdlrEnforceRelaxSol(set->conshdlrs_enfo[h], blkmem, set, stat, tree, sepastore,
               relaxsol, *infeasible, &result) );
      }
      else if( SCIPtreeHasFocusNodeLP(tree) )
      {

         SCIPsetDebugMsg(set, "enforce LP solution with value %g\n", SCIPlpGetObjval(lp, set, prob));

         assert(lp->flushed);
         assert(lp->solved);
         assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);
         SCIP_CALL( SCIPconshdlrEnforceLPSol(set->conshdlrs_enfo[h], blkmem, set, stat, tree, sepastore, *infeasible,
               &result) );
      }
      else
      {
         SCIP_CALL( SCIPconshdlrEnforcePseudoSol(set->conshdlrs_enfo[h], blkmem, set, stat, tree, branchcand, *infeasible,
               objinfeasible, forced, &result) );
         if( SCIPsepastoreGetNCuts(sepastore) != 0 )
         {
            SCIPerrorMessage("pseudo enforcing method of constraint handler <%s> separated cuts\n",
               SCIPconshdlrGetName(set->conshdlrs_enfo[h]));
            return SCIP_INVALIDRESULT;
         }
      }
      SCIPsetDebugMsg(set, "enforcing of <%s> returned result %d\n", SCIPconshdlrGetName(set->conshdlrs_enfo[h]), result);

      switch( result )
      {
      case SCIP_CUTOFF:
         assert(tree->nchildren == 0);
         *cutoff = TRUE;
         *infeasible = TRUE;
         resolved = TRUE;
         SCIPsetDebugMsg(set, " -> constraint handler <%s> detected cutoff in enforcement\n",
            SCIPconshdlrGetName(set->conshdlrs_enfo[h]));
         break;

      case SCIP_CONSADDED:
         assert(tree->nchildren == 0);
         *infeasible = TRUE;
         *propagateagain = TRUE; /* the propagation for new constraints should be called */
         *solvelpagain = TRUE;   /* the separation for new constraints should be called */
         *solverelaxagain = TRUE;
         markRelaxsUnsolved(set, relaxation);
         resolved = TRUE;
         break;

      case SCIP_REDUCEDDOM:
         assert(tree->nchildren == 0);
         *infeasible = TRUE;
         *propagateagain = TRUE;
         *solvelpagain = TRUE;
         *solverelaxagain = TRUE;
         markRelaxsUnsolved(set, relaxation);
         resolved = TRUE;
         break;

      case SCIP_SEPARATED:
         assert(tree->nchildren == 0);
         assert(SCIPsepastoreGetNCuts(sepastore) > 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE;
         *solverelaxagain = TRUE;
         markRelaxsUnsolved(set, relaxation);
         resolved = TRUE;
         break;

      case SCIP_BRANCHED:
         assert(tree->nchildren >= 1);
         assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         *branched = TRUE;
         resolved = TRUE;

         /* increase the number of internal nodes */
         stat->ninternalnodes++;
         stat->ntotalinternalnodes++;
         break;

      case SCIP_SOLVELP:
         assert(!SCIPtreeHasFocusNodeLP(tree));
         assert(tree->nchildren == 0);
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         *solvelpagain = TRUE;
         resolved = TRUE;
         SCIPtreeSetFocusNodeLP(tree, TRUE); /* the node's LP must be solved */
         break;

      case SCIP_INFEASIBLE:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         *infeasible = TRUE;
         break;

      case SCIP_FEASIBLE:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         break;

      case SCIP_DIDNOTRUN:
         assert(tree->nchildren == 0);
         assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
         assert(SCIPsepastoreGetNCuts(sepastore) == 0);
         assert(objinfeasible);
         *infeasible = TRUE;
         break;

      default:
         SCIPerrorMessage("invalid result code <%d> from enforcing method of constraint handler <%s>\n",
            result, SCIPconshdlrGetName(set->conshdlrs_enfo[h]));
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/

      /* the enforcement method may add a primal solution, after which the LP status could be set to
       * objective limit reached
       */
      if( SCIPtreeHasFocusNodeLP(tree) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT )
      {
         *cutoff = TRUE;
         *infeasible = TRUE;
         resolved = TRUE;
         SCIPsetDebugMsg(set, " -> LP exceeded objective limit\n");

         /* If we used the probing mode during branching, it might happen that we added a constraint or global bound
          * and returned SCIP_CONSADDED or SCIP_REDUCEDDOM, but when reoptimizing the LP after ending the probing mode,
          * this leads to hitting the objective limit. In this case, we do not need to propagate or solve the LP again.
          */
         *propagateagain = FALSE;
         *solvelpagain = FALSE;
      }

      assert(!(*branched) || (resolved && !(*cutoff) && *infeasible && !(*propagateagain) && !(*solvelpagain)));
      assert(!(*cutoff) || (resolved && !(*branched) && *infeasible && !(*propagateagain) && !(*solvelpagain)));
      assert(*infeasible || (!resolved && !(*branched) && !(*cutoff) && !(*propagateagain) && !(*solvelpagain)));
      assert(!(*propagateagain) || (resolved && !(*branched) && !(*cutoff) && *infeasible));
      assert(!(*solvelpagain) || (resolved && !(*branched) && !(*cutoff) && *infeasible));
   }
   assert(!objinfeasible || *infeasible);
   assert(resolved == (*branched || *cutoff || *propagateagain || *solvelpagain));
   assert(*cutoff || *solvelpagain || SCIPsepastoreGetNCuts(sepastore) == 0);

   /* in case the relaxation solution was enforced, free the created solution */
   if( enforcerelaxsol )
   {
      SCIP_CALL( SCIPsolFree(&relaxsol, blkmem, primal) );
   }

   /* deactivate the cut forcing of the constraint enforcement */
   SCIPsepastoreEndForceCuts(sepastore);

   SCIPsetDebugMsg(set, " -> enforcing result: branched=%u, cutoff=%u, infeasible=%u, propagateagain=%u, solvelpagain=%u, resolved=%u\n",
      *branched, *cutoff, *infeasible, *propagateagain, *solvelpagain, resolved);

   return SCIP_OKAY;
}

/** applies the cuts stored in the separation store, or clears the store if the node can be cut off */
static
SCIP_RETCODE applyCuts(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reotimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_RELAXATION*      relaxation,         /**< relaxators */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             root,               /**< is this the initial root LP? */
   SCIP_EFFICIACYCHOICE  efficiacychoice,    /**< type of solution to base efficiacy computation on */
   SCIP_Bool*            cutoff,             /**< pointer to whether the node can be cut off */
   SCIP_Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   SCIP_Bool*            solvelpagain,       /**< pointer to store TRUE, if the node's LP has to be solved again */
   SCIP_Bool*            solverelaxagain     /**< pointer to store TRUE, if the node's relaxation has to be solved again */
   )
{
   assert(stat != NULL);
   assert(cutoff != NULL);
   assert(propagateagain != NULL);
   assert(solvelpagain != NULL);

   if( *cutoff )
   {
      /* the found cuts are of no use, because the node is infeasible anyway (or we have an error in the LP) */
      SCIP_CALL( SCIPsepastoreClearCuts(sepastore, blkmem, set, eventqueue, eventfilter, lp) );
   }
   else if( SCIPsepastoreGetNCuts(sepastore) > 0 )
   {
      SCIP_Longint olddomchgcount;
      int oldncutsapplied;

      olddomchgcount = stat->domchgcount;
      oldncutsapplied = SCIPsepastoreGetNCutsApplied(sepastore);
      SCIP_CALL( SCIPsepastoreApplyCuts(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
            eventqueue, eventfilter, cliquetable, root, efficiacychoice, cutoff) );
      *propagateagain = *propagateagain || (stat->domchgcount != olddomchgcount);
      *solvelpagain = TRUE;
      if( (stat->domchgcount != olddomchgcount) || (SCIPsepastoreGetNCutsApplied(sepastore) != oldncutsapplied) )
      {
         *solverelaxagain = TRUE;
         markRelaxsUnsolved(set, relaxation);
      }
   }

   return SCIP_OKAY;
}

/** updates the cutoff, propagateagain, and solverelaxagain status of the current solving loop */
static
void updateLoopStatus(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   depth,              /**< depth of current node */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   SCIP_Bool*            propagateagain,     /**< pointer to store TRUE, if domain propagation should be applied again */
   SCIP_Bool*            solverelaxagain     /**< pointer to store TRUE, if at least one relaxator should be called again */
   )
{
   SCIP_NODE* focusnode;
   int r;

   assert(set != NULL);
   assert(stat != NULL);
   assert(cutoff != NULL);
   assert(propagateagain != NULL);
   assert(solverelaxagain != NULL);

   /* check, if the path was cutoff */
   *cutoff = *cutoff || (tree->cutoffdepth <= depth);

   /* check if branching was already performed */
   if( tree->nchildren == 0 )
   {
      /* check, if the focus node should be repropagated */
      focusnode = SCIPtreeGetFocusNode(tree);
      *propagateagain = *propagateagain || SCIPnodeIsPropagatedAgain(focusnode);

      /* check, if one of the external relaxations should be solved again */
      for( r = 0; r < set->nrelaxs && !(*solverelaxagain); ++r )
         *solverelaxagain = *solverelaxagain || ( !SCIPrelaxIsSolved(set->relaxs[r], stat) );
   }
   else
   {
      /* if branching was performed, avoid another node loop iteration */
      *propagateagain = FALSE;
      *solverelaxagain = FALSE;
   }
}

/** propagate domains and solve relaxation and lp */
static
SCIP_RETCODE propAndSolve(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_MEM*             mem,                /**< block memory pools */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_CUTPOOL*         delayedcutpool,     /**< global delayed cut pool */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_NODE*            focusnode,          /**< focused node */
   int                   actdepth,           /**< depth in the b&b tree */
   SCIP_Bool             propagate,          /**< should we propagate */
   SCIP_Bool             solvelp,            /**< should we solve the lp */
   SCIP_Bool             solverelax,         /**< should we solve the relaxation */
   SCIP_Bool             forcedlpsolve,      /**< is there a need for a solve lp */
   SCIP_Bool             initiallpsolved,    /**< was the initial LP already solved? */
   SCIP_Bool             fullseparation,     /**< are we in the first prop-and-cut-and-price loop? */
   SCIP_Longint*         afterlpproplps,     /**< pointer to store the last LP count for which AFTERLP propagation was performed */
   SCIP_HEURTIMING*      heurtiming,         /**< timing for running heuristics after propagation call */
   int*                  nlperrors,          /**< pointer to store the number of lp errors */
   SCIP_Bool*            fullpropagation,    /**< pointer to store whether we want to do a fullpropagation next time */
   SCIP_Bool*            propagateagain,     /**< pointer to store whether we want to propagate again */
   SCIP_Bool*            lpsolved,           /**< pointer to store whether the lp was solved */
   SCIP_Bool*            relaxcalled,        /**< pointer to store whether a relaxator was called; initialized with last loop's result */
   SCIP_Bool*            solvelpagain,       /**< pointer to store whether we want to solve the lp again */
   SCIP_Bool*            solverelaxagain,    /**< pointer to store whether we want to solve the relaxation again */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            postpone,           /**< pointer to store whether the node should be postponed */
   SCIP_Bool*            unbounded,          /**< pointer to store whether the focus node is unbounded */
   SCIP_Bool*            stopped,            /**< pointer to store whether solving was interrupted */
   SCIP_Bool*            lperror,            /**< pointer to store TRUE, if an unresolved error in LP solving occured */
   SCIP_Bool*            pricingaborted,     /**< pointer to store TRUE, if the pricing was aborted and the lower bound must not be used */
   SCIP_Bool*            forcedenforcement   /**< pointer to store whether the enforcement of pseudo solution should be forced */
   )
{
   SCIP_Bool newinitconss;

   assert(set != NULL);
   assert(stat != NULL);
   assert(origprob != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(primal != NULL);
   assert(pricestore != NULL);
   assert(sepastore != NULL);
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);
   assert(branchcand != NULL);
   assert(cutpool != NULL);
   assert(delayedcutpool != NULL);
   assert(conflict != NULL);
   assert(SCIPconflictGetNConflicts(conflict) == 0);
   assert(eventfilter != NULL);
   assert(eventqueue != NULL);
   assert(focusnode != NULL);
   assert(heurtiming != NULL);
   assert(nlperrors != NULL);
   assert(fullpropagation != NULL);
   assert(propagateagain != NULL);
   assert(afterlpproplps != NULL);
   assert(lpsolved != NULL);
   assert(solvelpagain != NULL);
   assert(solverelaxagain != NULL);
   assert(cutoff != NULL);
   assert(postpone != NULL);
   assert(unbounded != NULL);
   assert(lperror != NULL);
   assert(pricingaborted != NULL);
   assert(forcedenforcement != NULL);

   newinitconss = FALSE;

   if( !(*cutoff) && !(*postpone) )
   {
      SCIP_Longint oldninitconssadded;
      SCIP_Longint oldnboundchgs;
      SCIP_Bool lpwasflushed;

      lpwasflushed = lp->flushed;
      oldnboundchgs = stat->nboundchgs;
      oldninitconssadded = stat->ninitconssadded;

      /* call after LP propagators */
      if( ((*afterlpproplps) < stat->nnodelps && (*lpsolved)) || (*relaxcalled) )
      {
         SCIP_CALL( propagateDomains(blkmem, set, stat, primal, tree, SCIPtreeGetCurrentDepth(tree), 0, *fullpropagation,
               SCIP_PROPTIMING_AFTERLPLOOP, cutoff, postpone) );
         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

         /* check, if the path was cutoff */
         *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);
         *afterlpproplps = stat->nnodelps;
         propagate = propagate || (stat->nboundchgs > oldnboundchgs);
      }

      /* call before LP propagators */
      if( propagate && !(*cutoff) )
      {
         SCIP_CALL( propagateDomains(blkmem, set, stat, primal, tree, SCIPtreeGetCurrentDepth(tree), 0, *fullpropagation,
               SCIP_PROPTIMING_BEFORELP, cutoff, postpone) );
         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);
      }

      newinitconss = (stat->ninitconssadded != oldninitconssadded);
      *fullpropagation = FALSE;

      /* check, if the path was cutoff */
      *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

      /* if the LP was flushed and is now no longer flushed, a bound change occurred, and the LP has to be resolved;
       * we also have to solve the LP if new intial constraints were added which need to be added to the LP
       */
      solvelp = solvelp || (lpwasflushed && (!lp->flushed || newinitconss));
      solverelax = solverelax || newinitconss;

      /* the number of bound changes was increased by the propagation call, thus the relaxation should be solved again */
      if( stat->nboundchgs > oldnboundchgs )
      {
         /* propagation might have changed the best bound of loose variables, thereby changing the loose objective value
          * which is added to the LP value; because of the loose status, the LP might not be reoptimized, but the lower
          * bound of the node needs to be updated
          */
         if( !solvelp && lp->flushed && lp->solved && SCIPprobAllColsInLP(transprob, set, lp) && SCIPlpIsRelax(lp) )
         {
            SCIP_CALL( SCIPnodeUpdateLowerboundLP(focusnode, set, stat, tree, transprob, origprob, lp) );
            SCIPsetDebugMsg(set, " -> new lower bound: %g (LP status: %d, LP obj: %g)\n",
               SCIPnodeGetLowerbound(focusnode), SCIPlpGetSolstat(lp), SCIPlpGetObjval(lp, set, transprob));

            if( SCIPtreeHasFocusNodeLP(tree) )
            {
               /* update node estimate */
               SCIP_CALL( updateEstimate(set, stat, tree, lp, branchcand) );

               if( actdepth == 0 && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL )
                  SCIPprobUpdateBestRootSol(transprob, set, stat, lp);
            }
         }

         solverelax = TRUE;
         markRelaxsUnsolved(set, relaxation);
      }

      /* update lower bound with the pseudo objective value, and cut off node by bounding */
      SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue,
            conflict, cliquetable, cutoff) );
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   if( *postpone )
      return SCIP_OKAY;

   /* call primal heuristics that are applicable after propagation loop before lp solve;
    * the first time we go here, we call the before node heuristics instead
    */
   if( !(*cutoff) && !SCIPtreeProbing(tree) )
   {
      /* if the heuristics find a new incumbent solution, propagate again */
      SCIP_CALL( SCIPprimalHeuristics(set, stat, transprob, primal, tree, lp, NULL, *heurtiming,
            FALSE, propagateagain, unbounded) );
      assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

      *heurtiming = SCIP_HEURTIMING_AFTERPROPLOOP;

      /* check if primal heuristics found a solution and we therefore reached a solution limit */
      if( SCIPsolveIsStopped(set, stat, FALSE) )
      {
         /* cppcheck-suppress unassignedVariable */
         SCIP_NODE* node;

         /* we reached a solution limit and do not want to continue the processing of the current node, but in order to
          * allow restarting the optimization process later, we need to create a "branching" with only one child node that
          * is a copy of the focusnode
          */
         SCIPtreeSetFocusNodeLP(tree, FALSE);
         SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, 1.0, focusnode->estimate) );
         assert(tree->nchildren >= 1);
         *stopped = TRUE;
         return SCIP_OKAY;
      }

      /* if diving produced an LP error, switch back to non-LP node */
      if( lp->resolvelperror )
      {
         SCIPtreeSetFocusNodeLP(tree, FALSE);
         lp->resolvelperror = FALSE;
      }

      if( *propagateagain )
      {
         *solvelpagain = solvelp;
         *solverelaxagain = solverelax;

         return SCIP_OKAY;
      }
   }

   /* solve external relaxations with non-negative priority */
   *relaxcalled = FALSE;
   if( solverelax && !(*cutoff) )
   {

      /* clear the storage of external branching candidates */
      SCIPbranchcandClearExternCands(branchcand);

      SCIP_CALL( solveNodeRelax(set, stat, tree, relaxation, transprob, origprob, actdepth, TRUE,
            cutoff, propagateagain, solvelpagain, solverelaxagain, relaxcalled) );
      assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

      /* check, if the path was cutoff */
      *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

      /* apply found cuts */
      SCIP_CALL( applyCuts(blkmem, set, stat, transprob, origprob, tree, reopt, lp, relaxation, sepastore, branchcand,
            eventqueue, eventfilter, cliquetable, (actdepth == 0), SCIP_EFFICIACYCHOICE_RELAX, cutoff, propagateagain,
            solvelpagain, solverelaxagain) );

      /* update lower bound with the pseudo objective value, and cut off node by bounding */
      SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue, conflict, cliquetable, cutoff) );
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   /* check, if we want to solve the LP at this node */
   if( solvelp && !(*cutoff) && SCIPtreeHasFocusNodeLP(tree) )
   {
      *lperror = FALSE;
      *unbounded = FALSE;

      /* solve the node's LP */
      SCIP_CALL( solveNodeLP(blkmem, set, messagehdlr, stat, mem, origprob, transprob, primal, tree, reopt, lp, relaxation, pricestore,
            sepastore, cutpool, delayedcutpool, branchcand, conflict, conflictstore, eventfilter, eventqueue, cliquetable,
            initiallpsolved, fullseparation, newinitconss, propagateagain, solverelaxagain, cutoff, unbounded, lperror, pricingaborted) );

      *lpsolved = TRUE;
      *solvelpagain = FALSE;
      SCIPsetDebugMsg(set, " -> LP status: %d, LP obj: %g, iter: %" SCIP_LONGINT_FORMAT ", count: %" SCIP_LONGINT_FORMAT "\n",
         SCIPlpGetSolstat(lp),
         *cutoff ? SCIPsetInfinity(set) : (*lperror ? -SCIPsetInfinity(set) : SCIPlpGetObjval(lp, set, transprob)),
         stat->nlpiterations, stat->lpcount);

      /* check, if the path was cutoff */
      *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

      /* if an error occured during LP solving, switch to pseudo solution */
      if( *lperror )
      {
         if( forcedlpsolve )
         {
            SCIPerrorMessage("(node %" SCIP_LONGINT_FORMAT ") unresolved numerical troubles in LP %" SCIP_LONGINT_FORMAT " cannot be dealt with\n",
               stat->nnodes, stat->nlps);
            return SCIP_LPERROR;
         }
         SCIPtreeSetFocusNodeLP(tree, FALSE);
         ++(*nlperrors);
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, actdepth == 0 ? SCIP_VERBLEVEL_HIGH : SCIP_VERBLEVEL_FULL,
            "(node %" SCIP_LONGINT_FORMAT ") unresolved numerical troubles in LP %" SCIP_LONGINT_FORMAT " -- using pseudo solution instead (loop %d)\n",
            stat->nnodes, stat->nlps, *nlperrors);
      }

      if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_TIMELIMIT || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_ITERLIMIT )
      {
         SCIPtreeSetFocusNodeLP(tree, FALSE);
         *forcedenforcement = TRUE;

         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, actdepth == 0 ? SCIP_VERBLEVEL_HIGH : SCIP_VERBLEVEL_FULL,
            "(node %" SCIP_LONGINT_FORMAT ") LP solver hit %s limit in LP %" SCIP_LONGINT_FORMAT " -- using pseudo solution instead\n",
            stat->nnodes, SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_TIMELIMIT ? "time" : "iteration", stat->nlps);
      }

      if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
      {
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %" SCIP_LONGINT_FORMAT ") LP relaxation is unbounded (LP %" SCIP_LONGINT_FORMAT ")\n", stat->nnodes, stat->nlps);
      }

      /* if we solve exactly, the LP claims to be infeasible but the infeasibility could not be proved,
       * we have to forget about the LP and use the pseudo solution instead
       */
      if( !(*cutoff) && !(*lperror) && (set->misc_exactsolve || *pricingaborted) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE
         && SCIPnodeGetLowerbound(focusnode) < primal->cutoffbound )
      {
         if( SCIPbranchcandGetNPseudoCands(branchcand) == 0 && transprob->ncontvars > 0 )
         {
            SCIPerrorMessage("(node %" SCIP_LONGINT_FORMAT ") could not prove infeasibility of LP %" SCIP_LONGINT_FORMAT " (exactsolve=%u, pricingaborted=%u), all variables are fixed, %d continuous vars\n",
               stat->nnodes, stat->nlps, set->misc_exactsolve, *pricingaborted, transprob->ncontvars);
            SCIPerrorMessage("(node %" SCIP_LONGINT_FORMAT ")  -> have to call PerPlex() (feature not yet implemented)\n", stat->nnodes);
            /**@todo call PerPlex */
            return SCIP_LPERROR;
         }
         else
         {
            SCIPtreeSetFocusNodeLP(tree, FALSE);
            *forcedenforcement = TRUE;

            SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "(node %" SCIP_LONGINT_FORMAT ") could not prove infeasibility of LP %" SCIP_LONGINT_FORMAT " (exactsolve=%u, pricingaborted=%u) -- using pseudo solution (%d unfixed vars) instead\n",
               stat->nnodes, stat->nlps, set->misc_exactsolve, *pricingaborted, SCIPbranchcandGetNPseudoCands(branchcand));
         }
      }

      /* update lower bound with the pseudo objective value, and cut off node by bounding */
      SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue, conflict,
            cliquetable, cutoff) );
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);
   assert(*cutoff || !SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));

   /* reset solverelaxagain if no relaxations were solved up to this point (the LP-changes are already included in
    * relaxators called after the LP)
    */
   *solverelaxagain = *solverelaxagain && *relaxcalled;

   /* solve external relaxations with negative priority */
   if( solverelax && !(*cutoff) )
   {
      SCIP_CALL( solveNodeRelax(set, stat, tree, relaxation, transprob, origprob, actdepth, FALSE, cutoff,
            propagateagain, solvelpagain, solverelaxagain, relaxcalled) );
      assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

      /* check, if the path was cutoff */
      *cutoff = *cutoff || (tree->cutoffdepth <= actdepth);

      /* apply found cuts */
      SCIP_CALL( applyCuts(blkmem, set, stat, transprob, origprob, tree, reopt, lp, relaxation, sepastore, branchcand,
            eventqueue, eventfilter, cliquetable, (actdepth == 0), SCIP_EFFICIACYCHOICE_RELAX, cutoff, propagateagain,
            solvelpagain, solverelaxagain) );

      /* update lower bound with the pseudo objective value, and cut off node by bounding */
      SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue, conflict,
            cliquetable, cutoff) );
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);

   return SCIP_OKAY;
}

/** check if a restart can be performed */
#ifndef NDEBUG
static
SCIP_Bool restartAllowed(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(set != NULL);
   assert(stat != NULL);

   return (set->nactivepricers == 0 && !set->reopt_enable
         && (set->presol_maxrestarts == -1 || stat->nruns <= set->presol_maxrestarts)
         && (set->limit_restarts == -1 || stat->nruns <= set->limit_restarts));
}
#else
#define restartAllowed(set,stat)             ((set)->nactivepricers == 0 && !set->reopt_enable && ((set)->presol_maxrestarts == -1 || (stat)->nruns <= (set)->presol_maxrestarts) \
                                                && (set->limit_restarts == -1 || stat->nruns <= set->limit_restarts))
#endif

/** solves the focus node */
static
SCIP_RETCODE solveNode(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_MEM*             mem,                /**< block memory pools */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_CUTPOOL*         delayedcutpool,     /**< global delayed cut pool */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   SCIP_Bool*            postpone,           /**< pointer to store whether the node should be postponed */
   SCIP_Bool*            unbounded,          /**< pointer to store whether the focus node is unbounded */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the focus node's solution is infeasible */
   SCIP_Bool*            restart,            /**< should solving process be started again with presolving? */
   SCIP_Bool*            afternodeheur,      /**< pointer to store whether AFTERNODE heuristics were already called */
   SCIP_Bool*            stopped             /**< pointer to store whether solving was interrupted */
   )
{
   SCIP_NODE* focusnode;
   SCIP_Longint lastdomchgcount;
   SCIP_Longint afterlpproplps;
   SCIP_Real restartfac;
   SCIP_Longint lastlpcount;
   SCIP_HEURTIMING heurtiming;
   int actdepth;
   int nlperrors;
   int nloops;
   SCIP_Bool foundsol = FALSE;
   SCIP_Bool focusnodehaslp;
   SCIP_Bool lpsolved;
   SCIP_Bool initiallpsolved;
   SCIP_Bool fullseparation;
   SCIP_Bool solverelaxagain;
   SCIP_Bool solvelpagain;
   SCIP_Bool propagateagain;
   SCIP_Bool fullpropagation;
   SCIP_Bool branched;
   SCIP_Bool forcedlpsolve;
   SCIP_Bool wasforcedlpsolve;
   SCIP_Bool pricingaborted;

   assert(set != NULL);
   assert(stat != NULL);
   assert(origprob != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);
   assert(primal != NULL);
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);
   assert(SCIPconflictGetNConflicts(conflict) == 0);
   assert(cutoff != NULL);
   assert(postpone != NULL);
   assert(unbounded != NULL);
   assert(infeasible != NULL);
   assert(restart != NULL);
   assert(afternodeheur != NULL);

   *cutoff = FALSE;
   *postpone = FALSE;
   *unbounded = FALSE;
   *infeasible = FALSE;
   *restart = FALSE;
   *afternodeheur = FALSE;
   *stopped = FALSE;
   pricingaborted = FALSE;

   focusnode = SCIPtreeGetFocusNode(tree);
   assert(focusnode != NULL);
   assert(SCIPnodeGetType(focusnode) == SCIP_NODETYPE_FOCUSNODE);
   actdepth = SCIPnodeGetDepth(focusnode);

   /* invalidate relaxation solution */
   SCIPrelaxationSetSolValid(relaxation, FALSE, FALSE);

   /* clear the storage of external branching candidates */
   SCIPbranchcandClearExternCands(branchcand);

   SCIPsetDebugMsg(set, "Processing node %" SCIP_LONGINT_FORMAT " in depth %d, %d siblings\n",
      stat->nnodes, actdepth, tree->nsiblings);
   SCIPsetDebugMsg(set, "current pseudosolution: obj=%g\n", SCIPlpGetPseudoObjval(lp, set, transprob));
   /*debug(SCIPprobPrintPseudoSol(transprob, set));*/

   /* check, if we want to solve the LP at the selected node:
    * - solve the LP, if the lp solve depth and frequency demand solving
    * - solve the root LP, if the LP solve frequency is set to 0
    * - solve the root LP, if there are continuous variables present
    * - don't solve the node if its cut off by the pseudo objective value anyway
    */
   focusnodehaslp = (set->lp_solvedepth == -1 || actdepth <= set->lp_solvedepth);
   focusnodehaslp = focusnodehaslp && (set->lp_solvefreq >= 1 && actdepth % set->lp_solvefreq == 0);
   focusnodehaslp = focusnodehaslp || (actdepth == 0 && set->lp_solvefreq == 0);
   focusnodehaslp = focusnodehaslp && SCIPsetIsLT(set, SCIPlpGetPseudoObjval(lp, set, transprob), primal->cutoffbound);
   focusnodehaslp = set->reopt_enable ? focusnodehaslp && SCIPreoptGetSolveLP(reopt, set, focusnode) : focusnodehaslp;
   SCIPtreeSetFocusNodeLP(tree, focusnodehaslp);

   /* external node solving loop:
    *  - propagate domains
    *  - solve SCIP_LP
    *  - enforce constraints
    * if a constraint handler adds constraints to enforce its own constraints, both, propagation and LP solving
    * is applied again (if applicable on current node); however, if the new constraints don't have the enforce flag set,
    * it is possible, that the current infeasible solution is not cut off; in this case, we have to declare the solution
    * infeasible and perform a branching
    */
   lastdomchgcount = stat->domchgcount;
   lastlpcount = stat->lpcount;
   initiallpsolved = FALSE;
   fullseparation = TRUE;
   heurtiming = SCIP_HEURTIMING_BEFORENODE;
   nlperrors = 0;
   stat->npricerounds = 0;
   stat->nseparounds = 0;
   stat->nincseparounds = 0;
   solverelaxagain = TRUE;
   solvelpagain = TRUE;
   propagateagain = TRUE;
   fullpropagation = TRUE;
   forcedlpsolve = FALSE;
   nloops = 0;

   while( !(*cutoff) && !(*postpone) && (solverelaxagain || solvelpagain || propagateagain) && nlperrors < MAXNLPERRORS && !(*restart) )
   {
      SCIP_Bool lperror;
      SCIP_Bool solverelax;
      SCIP_Bool solvelp;
      SCIP_Bool propagate;
      SCIP_Bool forcedenforcement;
      SCIP_Bool relaxcalled;

      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      *unbounded = FALSE;
      *infeasible = FALSE;
      foundsol = FALSE;

      nloops++;
      lperror = FALSE;
      lpsolved = FALSE;
      relaxcalled = FALSE;
      forcedenforcement = FALSE;
      afterlpproplps = -1L;

      while( !lperror && !(*cutoff) && (propagateagain || solvelpagain || solverelaxagain
            || (afterlpproplps < stat->nnodelps && lpsolved) || relaxcalled) )
      {
         solverelax = solverelaxagain;
         solverelaxagain = FALSE;
         solvelp = solvelpagain;
         solvelpagain = FALSE;
         propagate = propagateagain;
         propagateagain = FALSE;

         /* update lower bound with the pseudo objective value, and cut off node by bounding */
         SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue,
               conflict, cliquetable, cutoff) );

         /* propagate domains before lp solving and solve relaxation and lp */
         SCIPsetDebugMsg(set, " -> node solving loop: call propagators that are applicable before%s LP is solved\n",
            lpsolved ? " and after" : "");
         SCIP_CALL( propAndSolve(blkmem, set, messagehdlr, stat, mem, origprob, transprob, primal, tree, reopt, lp,
               relaxation, pricestore, sepastore, branchcand, cutpool, delayedcutpool, conflict, conflictstore, eventfilter,
               eventqueue, cliquetable, focusnode, actdepth, propagate, solvelp, solverelax, forcedlpsolve, initiallpsolved,
               fullseparation, &afterlpproplps, &heurtiming, &nlperrors, &fullpropagation, &propagateagain, &lpsolved, &relaxcalled,
               &solvelpagain, &solverelaxagain, cutoff, postpone, unbounded, stopped, &lperror, &pricingaborted, &forcedenforcement) );
         initiallpsolved |= lpsolved;

         /* time or solution limit was hit and we already created a dummy child node to terminate fast */
         if( *stopped )
            return SCIP_OKAY;

      }
      fullseparation = FALSE;

      /* update the cutoff, propagateagain, and solverelaxagain status of current solving loop */
      updateLoopStatus(set, stat, tree, actdepth, cutoff, &propagateagain, &solverelaxagain);

      /* call primal heuristics that should be applied after the LP relaxation of the node was solved;
       * if this is the first loop of the root node, call also AFTERNODE heuristics already here, since they might help
       * to improve the primal bound, thereby producing additional reduced cost strengthenings and strong branching
       * bound fixings which also might lead to a restart
       */
      if( !(*postpone) && (!(*cutoff) || SCIPtreeGetNNodes(tree) > 0) )
      {
         if( actdepth == 0 && !(*afternodeheur) )
         {
            SCIP_CALL( SCIPprimalHeuristics(set, stat, transprob, primal, tree, lp, NULL,
                  SCIP_HEURTIMING_AFTERLPLOOP | SCIP_HEURTIMING_AFTERNODE, *cutoff, &foundsol, unbounded) );
            *afternodeheur = TRUE; /* the AFTERNODE heuristics should not be called again after the node */
         }
         else if( lpsolved || SCIPrelaxationIsSolValid(relaxation) )
         {
            SCIP_CALL( SCIPprimalHeuristics(set, stat, transprob, primal, tree, lp, NULL, SCIP_HEURTIMING_AFTERLPLOOP,
                  *cutoff, &foundsol, unbounded) );
         }
         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

         /* heuristics might have found a solution or set the cutoff bound such that the current node is cut off */
         SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue, conflict, cliquetable, cutoff) );
      }

      /* check if heuristics leave us with an invalid LP */
      if( lp->resolvelperror )
      {
         if( forcedlpsolve )
         {
            SCIPerrorMessage("(node %" SCIP_LONGINT_FORMAT ") unresolved numerical troubles in LP %" SCIP_LONGINT_FORMAT " cannot be dealt with\n",
               stat->nnodes, stat->nlps);
            return SCIP_LPERROR;
         }
         SCIPtreeSetFocusNodeLP(tree, FALSE);
         lp->resolvelperror = FALSE;
         nlperrors++;
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "(node %" SCIP_LONGINT_FORMAT ") unresolved numerical troubles in LP %" SCIP_LONGINT_FORMAT " -- using pseudo solution instead (loop %d)\n",
            stat->nnodes, stat->nlps, nlperrors);
      }

      if( pricingaborted && !(*cutoff) && SCIPlpGetSolstat(lp) != SCIP_LPSOLSTAT_OPTIMAL )
      {

         SCIPtreeSetFocusNodeLP(tree, FALSE);

         /* if we just ran into the time limit this is not really a numerical trouble;
          * however, if this is not the case, we print messages about numerical troubles in the current LP
          */
         if( !SCIPsolveIsStopped(set, stat, FALSE) )
         {
            if( forcedlpsolve )
            {
               SCIPerrorMessage("(node %" SCIP_LONGINT_FORMAT ") unresolved numerical troubles in LP %" SCIP_LONGINT_FORMAT " cannot be dealt with\n",
                  stat->nnodes, stat->nlps);
               return SCIP_LPERROR;
            }
            nlperrors++;
            SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "(node %" SCIP_LONGINT_FORMAT ") unresolved numerical troubles in LP %" SCIP_LONGINT_FORMAT " -- using pseudo solution instead (loop %d)\n",
               stat->nnodes, stat->nlps, nlperrors);
         }
      }

      /* if an improved solution was found, propagate and solve the relaxations again */
      if( foundsol )
      {
         propagateagain = TRUE;
         solvelpagain = TRUE;
         solverelaxagain = TRUE;
         markRelaxsUnsolved(set, relaxation);
      }

      /* enforce constraints */
      branched = FALSE;
      if( !(*postpone) && !(*cutoff) && !solverelaxagain && !solvelpagain && !propagateagain )
      {
         /* if the solution changed since the last enforcement, we have to completely reenforce it; otherwise, we
          * only have to enforce the additional constraints added in the last enforcement, but keep the infeasible
          * flag TRUE in order to not declare the infeasible solution feasible due to disregarding the already
          * enforced constraints
          */
         if( lastdomchgcount != stat->domchgcount || lastlpcount != stat->lpcount )
         {
            lastdomchgcount = stat->domchgcount;
            lastlpcount = stat->lpcount;
            *infeasible = FALSE;
         }

         /* call constraint enforcement */
         SCIP_CALL( enforceConstraints(blkmem, set, messagehdlr, stat, transprob, primal, tree, lp, relaxation, sepastore,
               branchcand, &branched, cutoff, infeasible, &propagateagain, &solvelpagain, &solverelaxagain,
               forcedenforcement) );
         assert(branched == (tree->nchildren > 0));
         assert(!branched || (!(*cutoff) && *infeasible && !propagateagain && !solvelpagain));
         assert(!(*cutoff) || (!branched && *infeasible && !propagateagain && !solvelpagain));
         assert(*infeasible || (!branched && !(*cutoff) && !propagateagain && !solvelpagain));
         assert(!propagateagain || (!branched && !(*cutoff) && *infeasible));
         assert(!solvelpagain || (!branched && !(*cutoff) && *infeasible));

         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

         /* apply found cuts */
         SCIP_CALL( applyCuts(blkmem, set, stat, transprob, origprob, tree, reopt, lp, relaxation, sepastore, branchcand,
               eventqueue, eventfilter, cliquetable, (actdepth == 0), SCIP_EFFICIACYCHOICE_LP, cutoff, &propagateagain,
               &solvelpagain, &solverelaxagain) );

         /* update lower bound with the pseudo objective value, and cut off node by bounding */
         SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue, conflict, cliquetable, cutoff) );

         /* update the cutoff, propagateagain, and solverelaxagain status of current solving loop */
         updateLoopStatus(set, stat, tree, actdepth, cutoff, &propagateagain, &solverelaxagain);
      }
      assert(SCIPsepastoreGetNCuts(sepastore) == 0);

      /* The enforcement detected no infeasibility, so, no branching was performed,
       * but the pricing was aborted and the current feasible solution does not have to be the
       * best solution in the current subtree --> we have to do a pseudo branching,
       * so we set infeasible TRUE and add the current solution to the solution pool
       */
      if( pricingaborted && !(*infeasible) && !(*cutoff) && !(*postpone) )
      {
         SCIP_Longint oldnbestsolsfound = primal->nbestsolsfound;
         SCIP_SOL* sol;
         SCIP_Bool stored;

         /* in case the relaxation was enforced add this solution, otherwise decide between LP and pseudosol */
         if( SCIPrelaxationIsSolValid(relaxation) && SCIPrelaxationIsLpIncludedForSol(relaxation) && (!SCIPtreeHasFocusNodeLP(tree)
               || SCIPsetIsGT(set, SCIPrelaxationGetSolObj(relaxation), SCIPlpGetObjval(lp, set, transprob))) )
         {
            SCIP_SOL* relaxsol;

            SCIP_CALL( SCIPsolCreateRelaxSol(&relaxsol, blkmem, set, stat, primal, tree, relaxation, NULL) );

            SCIP_CALL( SCIPprimalTrySol(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
                  eventqueue, eventfilter, relaxsol, FALSE, FALSE, TRUE, TRUE, TRUE,
                  &stored) );

            SCIP_CALL( SCIPsolFree(&relaxsol, blkmem, primal) );

            if( stored )
            {
               stat->nrelaxsolsfound++;

               if( primal->nbestsolsfound != oldnbestsolsfound )
               {
                  stat->nrelaxbestsolsfound++;
                  SCIPstoreSolutionGap(set->scip);
               }
            }
         }
         else
         {
            SCIP_CALL( SCIPsolCreateCurrentSol(&sol, blkmem, set, stat, transprob, primal, tree, lp, NULL) );
            SCIP_CALL( SCIPprimalTrySolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
                  eventqueue, eventfilter, &sol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );

            if( stored )
            {
               stat->nlpsolsfound++;

               if( primal->nbestsolsfound != oldnbestsolsfound )
               {
                  stat->nlpbestsolsfound++;
                  SCIPstoreSolutionGap(set->scip);
               }
            }
         }

         *infeasible = TRUE;
      }

      /* if the node is infeasible, but no constraint handler could resolve the infeasibility
       * -> branch on LP, external candidates, or the pseudo solution
       * -> e.g. select non-fixed binary or integer variable x with value x', create three
       *    sons: x <= x'-1, x = x', and x >= x'+1.
       *    In the left and right branch, the current solution is cut off. In the middle
       *    branch, the constraints can hopefully reduce domains of other variables to cut
       *    off the current solution.
       * In LP branching, we cannot allow adding constraints, because this does not necessary change the LP and can
       * therefore lead to an infinite loop.
       */
      wasforcedlpsolve = forcedlpsolve;
      forcedlpsolve = FALSE;
      if( (*infeasible) && !(*cutoff) && !(*postpone)
         && (!(*unbounded) || SCIPbranchcandGetNExternCands(branchcand) > 0 || SCIPbranchcandGetNPseudoCands(branchcand) > 0)
         && !solverelaxagain && !solvelpagain && !propagateagain && !branched )
      {
         SCIP_RESULT result = SCIP_DIDNOTRUN;
         int nlpcands = 0;

         if( SCIPtreeHasFocusNodeLP(tree) )
         {
            SCIP_CALL( SCIPbranchcandGetLPCands(branchcand, set, stat, lp, NULL, NULL, NULL, &nlpcands, NULL, NULL) );
         }

         if ( nlpcands > 0 || SCIPbranchcandGetNExternCands(branchcand) > 0 )
         {
            /* If there are LP candidates and their maximal priority is at least the maximal priority of the external
             * candidates, then branch on the LP candidates. Note that due to implicit integer variables,
             * SCIPbranchcandGetLPMaxPrio(branchcand) might be finite and SCIPbranchcandGetNPrioLPCands(branchcand) > 0,
             * but nlpcands == 0. */
            if ( SCIPbranchcandGetLPMaxPrio(branchcand) >= SCIPbranchcandGetExternMaxPrio(branchcand) && nlpcands > 0 )
            {
               assert( SCIPbranchcandGetNPrioLPCands(branchcand) > 0 );
               assert( nlpcands > 0 );

               /* branch on LP solution */
               SCIPsetDebugMsg(set, "infeasibility in depth %d was not resolved: branch on LP solution with %d fractionals\n",
                  SCIPnodeGetDepth(focusnode), nlpcands);
               SCIP_CALL( SCIPbranchExecLP(blkmem, set, stat, transprob, origprob, tree, reopt, lp, sepastore, branchcand,
                     eventqueue, primal->cutoffbound, FALSE, &result) );
               assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);
               assert(result != SCIP_DIDNOTRUN && result != SCIP_DIDNOTFIND);
            }
            else
            {
               assert( SCIPbranchcandGetNPrioExternCands(branchcand) > 0 );
               assert( SCIPbranchcandGetNExternCands(branchcand) > 0 );

               /* branch on external candidates */
               SCIPsetDebugMsg(set, "infeasibility in depth %d was not resolved: branch on %d external branching candidates.\n",
                  SCIPnodeGetDepth(focusnode), SCIPbranchcandGetNExternCands(branchcand));
               SCIP_CALL( SCIPbranchExecExtern(blkmem, set, stat, transprob, origprob, tree, reopt, lp, sepastore, branchcand,
                     eventqueue, primal->cutoffbound, TRUE, &result) );
               assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);
            }
         }

         if( result == SCIP_DIDNOTRUN || result == SCIP_DIDNOTFIND )
         {
            /* branch on pseudo solution */
            SCIPsetDebugMsg(set, "infeasibility in depth %d was not resolved: branch on pseudo solution with %d unfixed integers\n",
               SCIPnodeGetDepth(focusnode), SCIPbranchcandGetNPseudoCands(branchcand));
            SCIP_CALL( SCIPbranchExecPseudo(blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
                  primal->cutoffbound, TRUE, &result) );
            assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);
         }

         /* SCIP cannot guarantee convergence if it is necessary to branch on unbounded variables */
         if( result == SCIP_BRANCHED )
         {
            SCIP_VAR* var = stat->lastbranchvar;

            if( var != NULL && !stat->branchedunbdvar && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
               && (SCIPsetIsInfinity(set, -SCIPvarGetLbLocal(var)) || SCIPsetIsInfinity(set, SCIPvarGetUbLocal(var))) )
            {
               SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
                  "Starting spatial branch-and-bound on unbounded variable <%s> ([%g,%g]) - cannot guarantee finite termination.\n",
                  SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
                  stat->branchedunbdvar = TRUE;
            }
         }

         switch( result )
         {
         case SCIP_CUTOFF:
            assert(tree->nchildren == 0);
            *cutoff = TRUE;
            SCIPsetDebugMsg(set, " -> branching rule detected cutoff\n");
            break;
         case SCIP_CONSADDED:
            assert(tree->nchildren == 0);
            if( nlpcands > 0 )
            {
               SCIPerrorMessage("LP branching rule added constraint, which was not allowed this time\n");
               return SCIP_INVALIDRESULT;
            }
            propagateagain = TRUE;
            solvelpagain = TRUE;
            solverelaxagain = TRUE;
            markRelaxsUnsolved(set, relaxation);
            break;
         case SCIP_REDUCEDDOM:
            assert(tree->nchildren == 0);
            propagateagain = TRUE;
            solvelpagain = TRUE;
            solverelaxagain = TRUE;
            markRelaxsUnsolved(set, relaxation);
            break;
         case SCIP_SEPARATED:
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) > 0);
            solvelpagain = TRUE;
            solverelaxagain = TRUE;
            markRelaxsUnsolved(set, relaxation);
            break;
         case SCIP_BRANCHED:
            assert(tree->nchildren >= 1);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            branched = TRUE;

            /* increase the number of internal nodes */
            stat->ninternalnodes++;
            stat->ntotalinternalnodes++;
            break;
         case SCIP_DIDNOTFIND: /*lint -fallthrough*/
         case SCIP_DIDNOTRUN:
            /* all integer variables in the infeasible solution are fixed,
             * - if no continuous variables exist and all variables are known, the infeasible pseudo solution is completely
             *   fixed, and the node can be cut off
             * - if at least one continuous variable exists or we do not know all variables due to external pricers, we
             *   cannot resolve the infeasibility by branching -> solve LP (and maybe price in additional variables)
             */
            assert(tree->nchildren == 0);
            assert(SCIPsepastoreGetNCuts(sepastore) == 0);
            assert(SCIPbranchcandGetNPseudoCands(branchcand) == 0);

            if( transprob->ncontvars == 0 && set->nactivepricers == 0 )
            {
               *cutoff = TRUE;
               SCIPsetDebugMsg(set, " -> cutoff because all variables are fixed in current node\n");
            }
            else
            {
               /* feasible LP solutions with all integers fixed must be feasible
                * if also no external branching candidates were available
                */
               assert(!SCIPtreeHasFocusNodeLP(tree) || pricingaborted);

               if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_TIMELIMIT || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_ITERLIMIT || SCIPsolveIsStopped(set, stat, FALSE) )
               {
                  SCIP_NODE* node;

                  /* as we hit the time or iteration limit or another interrupt (e.g., gap limit), we do not want to solve the LP again.
                   * in order to terminate correctly, we create a "branching" with only one child node
                   * that is a copy of the focusnode
                   */
                  SCIP_CALL( SCIPnodeCreateChild(&node, blkmem, set, stat, tree, 1.0, focusnode->estimate) );
                  assert(tree->nchildren >= 1);
                  assert(SCIPsepastoreGetNCuts(sepastore) == 0);
                  branched = TRUE;
               }
               else
               {
                  SCIP_VERBLEVEL verblevel;

                  if( pricingaborted )
                  {
                     SCIPerrorMessage("pricing was aborted, but no branching could be created!\n");
                     return SCIP_INVALIDRESULT;
                  }

                  if( wasforcedlpsolve )
                  {
                     assert(SCIPtreeHasFocusNodeLP(tree));
                     SCIPerrorMessage("LP was solved, all integers fixed, some constraint still infeasible, but no branching could be created!\n");
                     return SCIP_INVALIDRESULT;
                  }

                  verblevel = SCIP_VERBLEVEL_FULL;

                  if( !tree->forcinglpmessage && set->disp_verblevel == SCIP_VERBLEVEL_HIGH )
                  {
                     verblevel = SCIP_VERBLEVEL_HIGH;

                     /* remember that the forcing LP solving message was posted and do only post it again if the
                      * verblevel is SCIP_VERBLEVEL_FULL
                      */
                     tree->forcinglpmessage = TRUE;
                  }

                  SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, verblevel,
                     "(node %" SCIP_LONGINT_FORMAT ") forcing the solution of an LP (last LP %" SCIP_LONGINT_FORMAT ")...\n", stat->nnodes, stat->nlps);

                  /* solve the LP in the next loop */
                  SCIPtreeSetFocusNodeLP(tree, TRUE);
                  solvelpagain = TRUE;
                  forcedlpsolve = TRUE; /* this LP must be solved without error - otherwise we have to abort */
               }
            }
            break;
         default:
            SCIPerrorMessage("invalid result code <%d> from SCIPbranchLP(), SCIPbranchExt() or SCIPbranchPseudo()\n", result);
            return SCIP_INVALIDRESULT;
         }  /*lint !e788*/
         assert(*cutoff || solvelpagain || propagateagain || branched); /* something must have been done */
         assert(!(*cutoff) || (!solvelpagain && !propagateagain && !branched));
         assert(!solvelpagain || (!(*cutoff) && !branched));
         assert(!propagateagain || (!(*cutoff) && !branched));
         assert(!branched || (!solvelpagain && !propagateagain));
         assert(branched == (tree->nchildren > 0));

         /* apply found cuts */
         SCIP_CALL( applyCuts(blkmem, set, stat, transprob, origprob, tree, reopt, lp, relaxation, sepastore, branchcand,
               eventqueue, eventfilter, cliquetable, (actdepth == 0), SCIP_EFFICIACYCHOICE_LP, cutoff, &propagateagain,
               &solvelpagain, &solverelaxagain) );

         /* update lower bound with the pseudo objective value, and cut off node by bounding */
         SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue, conflict, cliquetable, cutoff) );

         /* update the cutoff, propagateagain, and solverelaxagain status of current solving loop */
         updateLoopStatus(set, stat, tree, actdepth, cutoff, &propagateagain, &solverelaxagain);
      }

      /* check for immediate restart */
      *restart = *restart || (actdepth == 0 && restartAllowed(set, stat) && (stat->userrestart
	    || (stat->nrootintfixingsrun > set->presol_immrestartfac * (transprob->nvars - transprob->ncontvars)
	       && (stat->nruns == 1 || transprob->nvars <= (1.0-set->presol_restartminred) * stat->prevrunnvars))) );

      SCIPsetDebugMsg(set, "node solving iteration %d finished: cutoff=%u, postpone=%u, propagateagain=%u, solverelaxagain=%u, solvelpagain=%u, nlperrors=%d, restart=%u\n",
         nloops, *cutoff, *postpone, propagateagain, solverelaxagain, solvelpagain, nlperrors, *restart);
   }
   assert(SCIPsepastoreGetNCuts(sepastore) == 0);
   assert(*cutoff || SCIPconflictGetNConflicts(conflict) == 0);

   /* flush the conflict set storage */
   SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue, cliquetable) );

   /* check for too many LP errors */
   if( nlperrors >= MAXNLPERRORS )
   {
      SCIPerrorMessage("(node %" SCIP_LONGINT_FORMAT ") unresolved numerical troubles in LP %" SCIP_LONGINT_FORMAT " -- aborting\n", stat->nnodes, stat->nlps);
      return SCIP_LPERROR;
   }

   /* check for final restart */
   restartfac = set->presol_subrestartfac;
   if( actdepth == 0 )
      restartfac = MIN(restartfac, set->presol_restartfac);
   *restart = *restart || (restartAllowed(set, stat) && (stat->userrestart
	 || (stat->nrootintfixingsrun > restartfac * (transprob->nvars - transprob->ncontvars)
	    && (stat->nruns == 1 || transprob->nvars <= (1.0-set->presol_restartminred) * stat->prevrunnvars))) );

   /* remember the last root LP solution */
   if( actdepth == 0 && !(*cutoff) && !(*unbounded) && !(*postpone) )
   {
      /* the root pseudo objective value and pseudo objective value should be equal in the root node */
      assert(SCIPsetIsFeasEQ(set, SCIPlpGetGlobalPseudoObjval(lp, set, transprob), SCIPlpGetPseudoObjval(lp, set, transprob)));

      SCIPprobStoreRootSol(transprob, set, stat, lp, SCIPtreeHasFocusNodeLP(tree));
   }

   /* check for cutoff */
   if( *cutoff )
   {
      SCIPsetDebugMsg(set, "node is cut off\n");

      SCIPnodeUpdateLowerbound(focusnode, stat, set, tree, transprob, origprob, SCIPsetInfinity(set));
      *infeasible = TRUE;
      SCIP_CALL( SCIPdebugRemoveNode(blkmem, set, focusnode) ); /*lint !e506 !e774*/
   }
   else if( !(*unbounded) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL && lp->looseobjvalinf == 0 )
   {
      /* update the regression statistic nlpbranchcands and LP objective value  */
      int nlpbranchcands;
      SCIP_Real lpobjval;

      /* get number of LP candidate variables */
      SCIP_CALL( SCIPbranchcandGetLPCands(branchcand, set, stat, lp, NULL, NULL, NULL, &nlpbranchcands, NULL, NULL) );

      /* get LP objective value */
      lpobjval = SCIPlpGetObjval(lp, set, transprob);
      assert(lpobjval != SCIP_INVALID); /*lint !e777*/

      /* add the observation to the regression */
      SCIPregressionAddObservation(stat->regressioncandsobjval, (SCIP_Real)nlpbranchcands, lpobjval);
   }

   return SCIP_OKAY;
}

/** if feasible, adds current solution to the solution storage */
static
SCIP_RETCODE addCurrentSolution(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_Bool             checksol            /**< should the solution be checked? */
   )
{
   SCIP_Longint oldnbestsolsfound = primal->nbestsolsfound;
   SCIP_SOL* sol;
   SCIP_Bool foundsol;

   /* found a feasible solution */
   if( SCIPrelaxationIsSolValid(relaxation) && SCIPrelaxationIsLpIncludedForSol(relaxation) && (!SCIPtreeHasFocusNodeLP(tree)
         || SCIPsetIsGT(set, SCIPrelaxationGetSolObj(relaxation), SCIPlpGetObjval(lp, set, transprob))) )
   {
      /* start clock for relaxation solutions */
      SCIPclockStart(stat->relaxsoltime, set);

      SCIP_CALL( SCIPsolCreateRelaxSol(&sol, blkmem, set, stat, primal, tree, relaxation, NULL) );

      SCIPsetDebugMsg(set, "found relaxation solution with objective %f\n", SCIPsolGetObj(sol, set, transprob, origprob));

      if( checksol || set->misc_exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         SCIP_CALL( SCIPprimalTrySolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
               eventqueue, eventfilter, &sol, FALSE, FALSE, TRUE, TRUE, TRUE, &foundsol) );
      }
      else
      {
         SCIP_CALL( SCIPprimalAddSolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
               eventqueue, eventfilter, &sol, &foundsol) );
      }

      if( foundsol )
      {
         stat->nrelaxsolsfound++;

         if( primal->nbestsolsfound != oldnbestsolsfound )
         {
            stat->nrelaxbestsolsfound++;
            SCIPstoreSolutionGap(set->scip);
         }
      }

      /* stop clock for relaxation solutions */
      SCIPclockStop(stat->relaxsoltime, set);

   }
   else if( SCIPtreeHasFocusNodeLP(tree) )
   {
      assert(lp->primalfeasible);

      /* start clock for LP solutions */
      SCIPclockStart(stat->lpsoltime, set);

      /* add solution to storage */
      SCIP_CALL( SCIPsolCreateLPSol(&sol, blkmem, set, stat, transprob, primal, tree, lp, NULL) );

      SCIPsetDebugMsg(set, "found lp solution with objective %f\n", SCIPsolGetObj(sol, set, transprob, origprob));

      if( checksol || set->misc_exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         SCIP_CALL( SCIPprimalTrySolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
               eventqueue, eventfilter, &sol, FALSE, FALSE, TRUE, TRUE, TRUE, &foundsol) );
      }
      else
      {
         SCIP_CALL( SCIPprimalAddSolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
               eventqueue, eventfilter, &sol, &foundsol) );
      }

      if( foundsol )
      {
         stat->nlpsolsfound++;

         if( primal->nbestsolsfound != oldnbestsolsfound )
         {
            stat->nlpbestsolsfound++;
            SCIPstoreSolutionGap(set->scip);
         }
      }

      /* stop clock for LP solutions */
      SCIPclockStop(stat->lpsoltime, set);
   }
   else
   {
      /* start clock for pseudo solutions */
      SCIPclockStart(stat->pseudosoltime, set);

      /* add solution to storage */
      SCIP_CALL( SCIPsolCreatePseudoSol(&sol, blkmem, set, stat, transprob, primal, tree, lp, NULL) );

      SCIPsetDebugMsg(set, "found pseudo solution with objective %f\n", SCIPsolGetObj(sol, set, transprob, origprob));

      if( checksol || set->misc_exactsolve )
      {
         /* if we want to solve exactly, we have to check the solution exactly again */
         SCIP_CALL( SCIPprimalTrySolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
               eventqueue, eventfilter, &sol, FALSE, FALSE, TRUE, TRUE, TRUE, &foundsol) );
      }
      else
      {
         SCIP_CALL( SCIPprimalAddSolFree(primal, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt, lp,
               eventqueue, eventfilter, &sol, &foundsol) );
      }

      /* stop clock for pseudo solutions */
      SCIPclockStop(stat->pseudosoltime, set);

      if( foundsol )
      {
         stat->npssolsfound++;

         if( primal->nbestsolsfound != oldnbestsolsfound )
         {
            stat->npsbestsolsfound++;
            SCIPstoreSolutionGap(set->scip);
         }
      }
   }

   return SCIP_OKAY;
}

/** main solving loop */
SCIP_RETCODE SCIPsolveCIP(
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_MEM*             mem,                /**< block memory pools */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cut pool */
   SCIP_CUTPOOL*         delayedcutpool,     /**< global delayed cut pool */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            restart             /**< should solving process be started again with presolving? */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODE* focusnode;
   SCIP_NODE* nextnode;
   SCIP_EVENT event;
   SCIP_Real restartfac;
   SCIP_Real restartconfnum;
   int nnodes;
   int depth;
   SCIP_Bool cutoff;
   SCIP_Bool postpone;
   SCIP_Bool unbounded;
   SCIP_Bool infeasible;
   SCIP_Bool foundsol;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(pricestore != NULL);
   assert(sepastore != NULL);
   assert(branchcand != NULL);
   assert(cutpool != NULL);
   assert(delayedcutpool != NULL);
   assert(primal != NULL);
   assert(eventfilter != NULL);
   assert(eventqueue != NULL);
   assert(restart != NULL);

   /* check for immediate restart (if problem solving marked to be restarted was aborted) */
   restartfac = set->presol_subrestartfac;
   if( SCIPtreeGetCurrentDepth(tree) == 0 )
      restartfac = MIN(restartfac, set->presol_restartfac);
   *restart = restartAllowed(set, stat) && (stat->userrestart
      || (stat->nrootintfixingsrun > restartfac * (transprob->nvars - transprob->ncontvars)
	 && (stat->nruns == 1 || transprob->nvars <= (1.0-set->presol_restartminred) * stat->prevrunnvars)) );

   /* calculate the number of successful conflict analysis calls that should trigger a restart */
   if( set->conf_restartnum > 0 )
   {
      int i;

      restartconfnum = (SCIP_Real)set->conf_restartnum;
      for( i = 0; i < stat->nconfrestarts; ++i )
         restartconfnum *= set->conf_restartfac;
   }
   else
      restartconfnum = SCIP_REAL_MAX;
   assert(restartconfnum >= 0.0);

   /* switch status to UNKNOWN */
   stat->status = SCIP_STATUS_UNKNOWN;

   focusnode = NULL;
   nextnode = NULL;
   unbounded = FALSE;
   postpone = FALSE;

   while( !SCIPsolveIsStopped(set, stat, TRUE) && !(*restart) )
   {
      SCIP_Longint nsuccessconflicts;
      SCIP_Bool afternodeheur;
      SCIP_Bool stopped;
      SCIP_Bool branched;

      assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

      foundsol = FALSE;
      infeasible = FALSE;

      do
      {
         /* update the memory saving flag, switch algorithms respectively */
         SCIPstatUpdateMemsaveMode(stat, set, messagehdlr, mem);

         /* get the current node selector */
         nodesel = SCIPsetGetNodesel(set, stat);

         /* inform tree about the current node selector */
         SCIP_CALL( SCIPtreeSetNodesel(tree, set, messagehdlr, stat, nodesel) );

         /* the next node was usually already selected in the previous solving loop before the primal heuristics were
          * called, because they need to know, if the next node will be a child/sibling (plunging) or not;
          * if the heuristics found a new best solution that cut off some of the nodes, the node selector must be called
          * again, because the selected next node may be invalid due to cut off
          */
         if( nextnode == NULL )
         {
            /* select next node to process */
            SCIP_CALL( SCIPnodeselSelect(nodesel, set, &nextnode) );
         }
         focusnode = nextnode;
         nextnode = NULL;
         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

         /* start node activation timer */
         SCIPclockStart(stat->nodeactivationtime, set);

         /* focus selected node */
         SCIP_CALL( SCIPnodeFocus(&focusnode, blkmem, set, messagehdlr, stat, transprob, origprob, primal, tree, reopt,
               lp, branchcand, conflict, conflictstore, eventfilter, eventqueue, cliquetable, &cutoff, FALSE, FALSE) );
         if( cutoff )
            stat->ndelayedcutoffs++;

         /* stop node activation timer */
         SCIPclockStop(stat->nodeactivationtime, set);

         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);
      }
      while( cutoff ); /* select new node, if the current one was located in a cut off subtree */

      assert(SCIPtreeGetCurrentNode(tree) == focusnode);
      assert(SCIPtreeGetFocusNode(tree) == focusnode);

      /* if no more node was selected, we finished optimization */
      if( focusnode == NULL )
      {
         assert(SCIPtreeGetNNodes(tree) == 0);
         break;
      }

      /* update maxdepth and node count statistics */
      depth = SCIPnodeGetDepth(focusnode);
      stat->maxdepth = MAX(stat->maxdepth, depth);
      stat->maxtotaldepth = MAX(stat->maxtotaldepth, depth);
      stat->nnodes++;
      stat->ntotalnodes++;

      /* update reference bound statistic, if available */
      if( SCIPsetIsGE(set, SCIPnodeGetLowerbound(focusnode), stat->referencebound) )
         stat->nnodesaboverefbound++;

      /* issue NODEFOCUSED event */
      SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFOCUSED) );
      SCIP_CALL( SCIPeventChgNode(&event, focusnode) );
      SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

      /* solve focus node */
      SCIP_CALL( solveNode(blkmem, set, messagehdlr, stat, mem, origprob, transprob, primal, tree, reopt, lp, relaxation,
            pricestore, sepastore, branchcand, cutpool, delayedcutpool, conflict, conflictstore, eventfilter, eventqueue,
            cliquetable, &cutoff, &postpone, &unbounded, &infeasible, restart, &afternodeheur, &stopped) );
      assert(!cutoff || infeasible);
      assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);
      assert(SCIPtreeGetCurrentNode(tree) == focusnode);
      assert(SCIPtreeGetFocusNode(tree) == focusnode);

      branched = (tree->nchildren > 0);

      if( stopped )
         break;

      /* check for restart */
      if( !(*restart) && !postpone )
      {
         /* change color of node in visualization */
         SCIPvisualSolvedNode(stat->visual, set, stat, focusnode);

         /* check, if the current solution is feasible */
         if( !infeasible )
         {
            SCIP_Bool feasible;

            assert(!SCIPtreeHasFocusNodeLP(tree) || (lp->flushed && lp->solved));
            assert(!cutoff);

            /* in the unbounded case, we check the solution w.r.t. the original problem, because we do not want to rely
             * on the LP feasibility and integrality is not checked for unbounded solutions, anyway
             */
            if( unbounded )
            {
               SCIP_SOL* sol;

               if( SCIPrelaxationIsSolValid(relaxation) && SCIPrelaxationIsLpIncludedForSol(relaxation) && (!SCIPtreeHasFocusNodeLP(tree)
                     || SCIPsetIsGT(set, SCIPrelaxationGetSolObj(relaxation), SCIPlpGetObjval(lp, set, transprob))) )
               {
                  SCIP_CALL( SCIPsolCreateRelaxSol(&sol, blkmem, set, stat, primal, tree, relaxation, NULL) );
               }
               else if( SCIPtreeHasFocusNodeLP(tree) )
               {
                  SCIP_CALL( SCIPsolCreateLPSol(&sol, blkmem, set, stat, transprob, primal, tree, lp, NULL) );
               }
               else
               {
                  SCIP_CALL( SCIPsolCreatePseudoSol(&sol, blkmem, set, stat, transprob, primal, tree, lp, NULL) );
               }
               SCIP_CALL( SCIPcheckSolOrig(set->scip, sol, &feasible, FALSE, FALSE) );

               SCIP_CALL( SCIPsolFree(&sol, blkmem, primal) );
            }
            else
               feasible = TRUE;

            /* node solution is feasible: add it to the solution store */
            if( feasible )
            {
               SCIP_CALL( addCurrentSolution(blkmem, set, messagehdlr, stat, origprob, transprob, primal, relaxation, tree, reopt,
                     lp, eventqueue, eventfilter, FALSE) );

               /* update the cutoff pointer if the new solution made the cutoff bound equal to the lower bound */
               SCIP_CALL( applyBounding(blkmem, set, stat, transprob, origprob, primal, tree, reopt, lp, branchcand, eventqueue, conflict, cliquetable, &cutoff) );


               /* increment number of feasible leaf nodes */
               stat->nfeasleaves++;

               /* issue NODEFEASIBLE event */
               SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEFEASIBLE) );
               SCIP_CALL( SCIPeventChgNode(&event, focusnode) );
               SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

               if( set->reopt_enable )
               {
                  assert(reopt != NULL);
                  SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, focusnode, SCIP_EVENTTYPE_NODEFEASIBLE, lp,
                        SCIPlpGetSolstat(lp), tree->root == focusnode, tree->focusnode == focusnode,
                        focusnode->lowerbound, tree->effectiverootdepth) );
               }
            }
         }
         else if( !unbounded )
         {
            /* node solution is not feasible */
            if( !branched )
            {
               assert(tree->nchildren == 0);

               /* change color of node in visualization output */
               SCIPvisualCutoffNode(stat->visual, set, stat, focusnode, TRUE);

               /* issue NODEINFEASIBLE event */
               SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEINFEASIBLE) );

               /* we only increase the number of objective leaf nodes if we hit the LP objective limit; we might have also
                * hit the objective limit at a node that is actually infeasible, or a dual reduction led to an infeasibility prior
                * to LP solving such that the node will be marked as infeasible */
               if( SCIPtreeHasCurrentNodeLP(tree) && SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OBJLIMIT )
                  stat->nobjleaves++;
               else
                  stat->ninfeasleaves++;

               if( set->reopt_enable )
               {
                  assert(reopt != NULL);
                  SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, focusnode, SCIP_EVENTTYPE_NODEINFEASIBLE, lp,
                        SCIPlpGetSolstat(lp), tree->root == focusnode, tree->focusnode == focusnode,
                        focusnode->lowerbound, tree->effectiverootdepth) );
               }

               /* increase the cutoff counter of the branching variable */
               if( stat->lastbranchvar != NULL )
               {
                  SCIP_CALL( SCIPvarIncCutoffSum(stat->lastbranchvar, blkmem, set, stat, stat->lastbranchdir, stat->lastbranchvalue, 1.0) );
               }
               /**@todo if last branching variable is unknown, retrieve it from the nodes' boundchg arrays */
            }
            else
            {
               assert(tree->nchildren > 0);

               /* issue NODEBRANCHED event */
               SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_NODEBRANCHED) );

               if( set->reopt_enable )
               {
                  assert(reopt != NULL);
                  SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, focusnode, SCIP_EVENTTYPE_NODEBRANCHED, lp,
                        SCIPlpGetSolstat(lp), tree->root == focusnode, tree->focusnode == focusnode,
                        focusnode->lowerbound, tree->effectiverootdepth) );
               }
            }
            SCIP_CALL( SCIPeventChgNode(&event, focusnode) );
            SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );
         }
         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

         /* if no branching was created, the node was not cut off, but its lower bound is still smaller than
          * the cutoff bound, we have to branch on a non-fixed variable;
          * this can happen, if we want to solve exactly, the current solution was declared feasible by the
          * constraint enforcement, but in exact solution checking it was found out to be infeasible;
          * in this case, no branching would have been generated by the enforcement of constraints, but we
          * have to further investigate the current sub tree;
          * note that we must noch check tree->nchildren > 0 here to determine whether we branched, we rather
          * check it directly after solveNode() and store the result, because an event handler might impose a
          * new cutoff bound (as is the case in ParaSCIP)
          */
         if( !cutoff && !unbounded && !branched && SCIPnodeGetLowerbound(focusnode) < primal->cutoffbound )
         {
            SCIP_RESULT result;

            assert(set->misc_exactsolve);

            do
            {
               result = SCIP_DIDNOTRUN;
               if( SCIPbranchcandGetNPseudoCands(branchcand) == 0 )
               {
                  if( transprob->ncontvars > 0 )
                  {
                     /**@todo call PerPlex */
                     SCIPerrorMessage("cannot branch on all-fixed LP -- have to call PerPlex instead\n");
                  }
               }
               else
               {
                  SCIP_CALL( SCIPbranchExecPseudo(blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
                        eventqueue, primal->cutoffbound, FALSE, &result) );
                  assert(result != SCIP_DIDNOTRUN && result != SCIP_DIDNOTFIND);
               }
            }
            while( result == SCIP_REDUCEDDOM );
         }
         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

         /* select node to process in next solving loop; the primal heuristics need to know whether a child/sibling
          * (plunging) will be selected as next node or not
          */
         SCIP_CALL( SCIPnodeselSelect(nodesel, set, &nextnode) );
         assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

         /* call primal heuristics that should be applied after the node was solved */
         nnodes = SCIPtreeGetNNodes(tree);
         stopped = SCIPsolveIsStopped(set, stat, FALSE);
         if( !afternodeheur && (!cutoff || nnodes > 0) && !stopped )
         {
            SCIP_CALL( SCIPprimalHeuristics(set, stat, transprob, primal, tree, lp, nextnode, SCIP_HEURTIMING_AFTERNODE,
                  cutoff, &foundsol, &unbounded) );
            assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

            stopped = SCIPsolveIsStopped(set, stat, FALSE);
         }

         /* if the heuristics found a new best solution that cut off some of the nodes, the node selector must be called
          * again, because the selected next node may be invalid due to cut off
          */
         assert(!tree->cutoffdelayed);

         if( nnodes != SCIPtreeGetNNodes(tree) || stopped )
            nextnode = NULL;
      }
      else if( !infeasible && !postpone )
      {
         /* The current solution was not proven to be infeasible, but due to the restart, this does not mean that it is
          * feasible, we might just have skipped the check. Thus, we try to add it to the solution store, but check it
          * again.
          */
         SCIP_CALL( addCurrentSolution(blkmem, set, messagehdlr, stat, origprob, transprob, primal, relaxation, tree, reopt, lp,
               eventqueue, eventfilter, TRUE) );

         if( set->reopt_enable )
         {
            assert(reopt != NULL);
            SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, focusnode, SCIP_EVENTTYPE_NODEFEASIBLE, lp,
                  SCIPlpGetSolstat(lp), tree->root == focusnode, tree->focusnode == focusnode, focusnode->lowerbound,
                  tree->effectiverootdepth) );
         }
      }
      /* we want to put the focusnode back into the leaf queue if it was postponed */
      else if( postpone )
      {
         SCIP_NODE* newfocusnode = NULL;

         /* @todo should we re-install the old focus node in order to somehow set the forks more clever? */
         SCIP_CALL( SCIPnodeFocus(&newfocusnode, blkmem, set, messagehdlr, stat, transprob, origprob, primal, tree, reopt, lp,
               branchcand, conflict, conflictstore, eventfilter, eventqueue, cliquetable, &cutoff, TRUE, FALSE) );
      }
      /* compute number of successfully applied conflicts */
      nsuccessconflicts = SCIPconflictGetNPropSuccess(conflict) + SCIPconflictGetNInfeasibleLPSuccess(conflict)
         + SCIPconflictGetNBoundexceedingLPSuccess(conflict) + SCIPconflictGetNStrongbranchSuccess(conflict)
         + SCIPconflictGetNPseudoSuccess(conflict);

      /* trigger restart due to conflicts and the restart parameters allow another restart */
      if( nsuccessconflicts >= restartconfnum && restartAllowed(set, stat) )
      {
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "(run %d, node %" SCIP_LONGINT_FORMAT ") restarting after %" SCIP_LONGINT_FORMAT " successful conflict analysis calls\n",
            stat->nruns, stat->nnodes, nsuccessconflicts);
         *restart = TRUE;

         stat->nconfrestarts++;
      }

      /* restart if the userrestart was set to true, we have still some nodes left and the restart parameters allow
       * another restart
       */
      *restart = *restart || (stat->userrestart && SCIPtreeGetNNodes(tree) > 0 && restartAllowed(set, stat));
      if( restartAllowed(set, stat) && set->limit_autorestartnodes == stat->nnodes && stat->ntotalnodes - stat->nruns + 1 == set->limit_autorestartnodes )
      {
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
               "(run %d, node %" SCIP_LONGINT_FORMAT ") restarting: triggering parameter controlled restart)\n",
               stat->nruns, stat->nnodes);
         *restart = TRUE;
      }
      /* if restart limit was exceeded, change the status; if status is different from unknown, ie some other limit was
       * hit, leave it unchanged
       */
      if( *restart && stat->status == SCIP_STATUS_UNKNOWN && set->limit_restarts >= 0 && stat->nruns > set->limit_restarts )
      {
         *restart = FALSE;
         stat->status = SCIP_STATUS_RESTARTLIMIT;
      }

      /* display node information line */
      SCIP_CALL( SCIPdispPrintLine(set, messagehdlr, stat, NULL, (SCIPnodeGetDepth(focusnode) == 0) && infeasible && !foundsol, TRUE) );

      SCIPsetDebugMsg(set, "Processing of node %" SCIP_LONGINT_FORMAT " in depth %d finished. %d siblings, %d children, %d leaves left\n",
         stat->nnodes, SCIPnodeGetDepth(focusnode), tree->nsiblings, tree->nchildren, SCIPtreeGetNLeaves(tree));
      SCIPsetDebugMsg(set, "**********************************************************************\n");
   }

   /* update the primal-dual integral if node or time limits were hit or an interruption signal was called */
   if( SCIPsolveIsStopped(set, stat, TRUE) )
   {
      SCIPstatUpdatePrimalDualIntegral(stat, set, transprob, origprob, SCIPsetInfinity(set), -SCIPsetInfinity(set));
   }

   assert(BMSgetNUsedBufferMemory(mem->buffer) == 0);

   SCIPsetDebugMsg(set, "Problem solving finished with status %u (restart=%u, userrestart=%u)\n", stat->status, *restart, stat->userrestart);

   /* cuts off nodes with lower bound is not better than given cutoff bound, manually; this necessary to ensure that
    * SCIP terminates with a proper solve stage
    */
   SCIP_CALL( SCIPtreeCutoff(tree, reopt, blkmem, set, stat, eventqueue, lp, primal->cutoffbound) );

   /* if the current node is the only remaining node, and if its lower bound exceeds the upper bound, we have
    * to delete it manually in order to get to the SOLVED stage instead of thinking, that only the gap limit
    * was reached (this may happen, if the current node is the one defining the global lower bound and a
    * feasible solution with the same value was found at this node)
    */
   if( tree->focusnode != NULL && SCIPtreeGetNNodes(tree) == 0
      && SCIPsetIsGE(set, tree->focusnode->lowerbound, primal->cutoffbound) )
   {
      if( set->reopt_enable )
      {
         assert(reopt != NULL);
         SCIP_CALL( SCIPreoptCheckCutoff(reopt, set, blkmem, tree->focusnode, SCIP_EVENTTYPE_NODEINFEASIBLE, lp,
               SCIPlpGetSolstat(lp), tree->root == focusnode, tree->focusnode == focusnode, tree->focusnode->lowerbound,
               tree->effectiverootdepth) );
      }

      focusnode = NULL;
      SCIP_CALL( SCIPnodeFocus(&focusnode, blkmem, set, messagehdlr, stat, transprob, origprob, primal, tree, reopt, lp,
            branchcand, conflict, conflictstore, eventfilter, eventqueue, cliquetable, &cutoff, FALSE, FALSE) );
   }

   /* check whether we finished solving */
   if( SCIPtreeGetNNodes(tree) == 0 && SCIPtreeGetCurrentNode(tree) == NULL )
   {
      /* no restart necessary */
      *restart = FALSE;

      /* set the solution status */
      if( unbounded )
      {
         if( primal->nsols > 0 )
         {
            /* switch status to UNBOUNDED */
            stat->status = SCIP_STATUS_UNBOUNDED;
         }
         else
         {
            /* switch status to INFORUNB */
            stat->status = SCIP_STATUS_INFORUNBD;
         }
      }
      else if( primal->nlimsolsfound == 0 )
      {
         assert(primal->nsols == 0 || SCIPsetIsGT(set, SCIPsolGetObj(primal->sols[0], set, transprob, origprob),
               SCIPprobInternObjval(transprob, origprob, set, SCIPprobGetObjlim(transprob, set))));

         /* switch status to INFEASIBLE */
         stat->status = SCIP_STATUS_INFEASIBLE;
      }
      else
      {
         /* switch status to OPTIMAL */
         stat->status = SCIP_STATUS_OPTIMAL;
      }
   }

   return SCIP_OKAY;
}
