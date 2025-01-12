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

/**@file   sepastore.c
 * @ingroup OTHER_CFILES
 * @brief  methods for storing separated cuts
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/sepastore.h"
#include "scip/event.h"
#include "scip/sepa.h"
#include "scip/cons.h"
#include "scip/debug.h"
#include "scip/scip.h"
#include "scip/cuts.h"
#include "scip/cutsel.h"
#include "scip/struct_event.h"
#include "scip/struct_sepastore.h"
#include "scip/misc.h"



/*
 * dynamic memory arrays
 */

/** resizes cuts and score arrays to be able to store at least num entries */
static
SCIP_RETCODE sepastoreEnsureCutsMem(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(sepastore != NULL);
   assert(set != NULL);

   if( num > sepastore->cutssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&sepastore->cuts, newsize) );
      sepastore->cutssize = newsize;
   }
   assert(num <= sepastore->cutssize);

   return SCIP_OKAY;
}

/** creates separation storage */
SCIP_RETCODE SCIPsepastoreCreate(
   SCIP_SEPASTORE**      sepastore,          /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepastore != NULL);

   SCIP_ALLOC( BMSallocMemory(sepastore) );

   (*sepastore)->cuts = NULL;
   (*sepastore)->cutssize = 0;
   (*sepastore)->ncuts = 0;
   (*sepastore)->nforcedcuts = 0;
   (*sepastore)->ncutsadded = 0;
   (*sepastore)->ncutsaddedviapool = 0;
   (*sepastore)->ncutsaddeddirect = 0;
   (*sepastore)->ncutsfoundround = 0;
   (*sepastore)->ncutsapplied = 0;
   (*sepastore)->initiallp = FALSE;
   (*sepastore)->forcecuts = FALSE;

   SCIP_CALL( SCIPrandomCreate(&(*sepastore)->randnumgen, blkmem, (unsigned int)SCIPsetInitializeRandomSeed(set, 0x5EED)) );

   return SCIP_OKAY;
}

/** frees separation storage */
SCIP_RETCODE SCIPsepastoreFree(
   SCIP_SEPASTORE**      sepastore,          /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(sepastore != NULL);
   assert(*sepastore != NULL);
   assert((*sepastore)->ncuts == 0);

   SCIPrandomFree(&(*sepastore)->randnumgen, blkmem);
   BMSfreeMemoryArrayNull(&(*sepastore)->cuts);
   BMSfreeMemory(sepastore);

   return SCIP_OKAY;
}

/** informs separation storage that the setup of the initial LP starts now */
void SCIPsepastoreStartInitialLP(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(!sepastore->initiallp);
   assert(sepastore->ncuts == 0);

   sepastore->initiallp = TRUE;
}

/** informs separation storage that the setup of the initial LP is now finished */
void SCIPsepastoreEndInitialLP(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(sepastore->initiallp);
   assert(sepastore->ncuts == 0);

   sepastore->initiallp = FALSE;
}

/** informs separation storage that the following cuts should be used in any case */
void SCIPsepastoreStartForceCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(!sepastore->forcecuts);

   sepastore->forcecuts = TRUE;
}

/** informs separation storage that the following cuts should no longer be used in any case */
void SCIPsepastoreEndForceCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);
   assert(sepastore->forcecuts);

   sepastore->forcecuts = FALSE;
}

/** checks cut for redundancy due to activity bounds */
static
SCIP_Bool sepastoreIsCutRedundant(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_ROW*             cut                 /**< separated cut */
   )
{
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(sepastore != NULL);
   assert(cut != NULL);

   /* modifiable cuts cannot be declared redundant, since we don't know all coefficients */
   if( SCIProwIsModifiable(cut) )
      return FALSE;

   /* check for activity redundancy */
   lhs = SCIProwGetLhs(cut);
   rhs = SCIProwGetRhs(cut);
   minactivity = SCIProwGetMinActivity(cut, set, stat);
   maxactivity = SCIProwGetMaxActivity(cut, set, stat);

   if( (SCIPsetIsInfinity(set, -lhs) || SCIPsetIsLE(set, lhs, minactivity)) &&
       (SCIPsetIsInfinity(set, rhs) || SCIPsetIsLE(set, maxactivity, rhs)) )
   {
      SCIPsetDebugMsg(set, "ignoring activity redundant cut <%s> (sides=[%g,%g], act=[%g,%g])\n",
         SCIProwGetName(cut), lhs, rhs, minactivity, maxactivity);
      /*SCIPdebug(SCIProwPrint(cut, set->scip->messagehdlr, NULL));*/
      return TRUE;
   }

   return FALSE;
}

/** checks cut for redundancy or infeasibility due to activity bounds */
static
SCIP_Bool sepastoreIsCutRedundantOrInfeasible(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool*            infeasible          /**< pointer to store whether the cut has been detected to be infeasible */
   )
{
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(sepastore != NULL);
   assert(cut != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* modifiable cuts cannot be declared redundant or infeasible, since we don't know all coefficients */
   if( SCIProwIsModifiable(cut) )
      return FALSE;

   /* check for activity redundancy */
   lhs = SCIProwGetLhs(cut);
   rhs = SCIProwGetRhs(cut);
   minactivity = SCIProwGetMinActivity(cut, set, stat);
   maxactivity = SCIProwGetMaxActivity(cut, set, stat);

   if( (SCIPsetIsInfinity(set, -lhs) || SCIPsetIsLE(set, lhs, minactivity)) &&
       (SCIPsetIsInfinity(set, rhs) || SCIPsetIsLE(set, maxactivity, rhs)) )
   {
      SCIPsetDebugMsg(set, "ignoring activity redundant cut <%s> (sides=[%g,%g], act=[%g,%g])\n",
         SCIProwGetName(cut), lhs, rhs, minactivity, maxactivity);
      /*SCIPdebug(SCIProwPrint(cut, set->scip->messagehdlr, NULL));*/
      return TRUE;
   }

   if( (!SCIPsetIsInfinity(set,  rhs) && SCIPsetIsFeasPositive(set, minactivity - rhs)) ||
       (!SCIPsetIsInfinity(set, -lhs) && SCIPsetIsFeasNegative(set, maxactivity - lhs)) )
   {
      SCIPsetDebugMsg(set, "cut <%s> is infeasible (sides=[%g,%g], act=[%g,%g])\n",
         SCIProwGetName(cut), lhs, rhs, minactivity, maxactivity);
      /*SCIPdebug(SCIProwPrint(cut, set->scip->messagehdlr, NULL));*/
      *infeasible = TRUE;
      return TRUE;
   }

   return FALSE;
}

/** checks whether a cut with only one variable can be applied as boundchange
 *
 *  This is the case if the bound change would prove infeasibility (w.r.t feastol), or if the new bound is at least
 *  epsilon better than the old bound.  In the latter case, also the opposite bound has to be taken into account.
 */
static
SCIP_Bool sepastoreIsBdchgApplicable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             cut                 /**< cut with a single variable */
   )
{
   SCIP_COL** cols;
   SCIP_Real* vals;
   SCIP_VAR* var;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool local;
   SCIP_Real oldlb;
   SCIP_Real oldub;

   assert(set != NULL);
   assert(cut != NULL);
   assert(!SCIProwIsModifiable(cut));
   assert(SCIProwGetNNonz(cut) == 1);

   /* get the single variable and its coefficient of the cut */
   cols = SCIProwGetCols(cut);
   assert(cols != NULL);

   var = SCIPcolGetVar(cols[0]);
   vals = SCIProwGetVals(cut);
   assert(vals != NULL);
   assert(!SCIPsetIsZero(set, vals[0]));

   /* if the coefficient is nearly zero, we better ignore this cut for numerical reasons */
   if( SCIPsetIsFeasZero(set, vals[0]) )
      return FALSE;

   local = SCIProwIsLocal(cut);

   oldlb = local ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
   oldub = local ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);

   /* get the left hand side of the cut and convert it to a bound */
   lhs = SCIProwGetLhs(cut);
   if( !SCIPsetIsInfinity(set, -lhs) )
   {
      lhs -= SCIProwGetConstant(cut);
      if( vals[0] > 0.0 )
      {
         /* coefficient is positive -> lhs corresponds to lower bound */
         SCIP_Real newlb;

         newlb = lhs/vals[0];
         SCIPvarAdjustLb(var, set, &newlb);

         /* bound changes that improve the bound sufficiently are applicable */
         if( SCIPsetIsFeasGT(set, newlb, oldub) || SCIPsetIsGT(set, MIN(newlb, oldub), oldlb) )
            return TRUE;
      }
      else
      {
         /* coefficient is negative -> lhs corresponds to upper bound */
         SCIP_Real newub;

         newub = lhs/vals[0];
         SCIPvarAdjustUb(var, set, &newub);

         /* bound changes that improve the bound sufficiently are applicable */
         if( SCIPsetIsFeasLT(set, newub, oldlb) || SCIPsetIsLT(set, MAX(newub, oldlb), oldub) )
            return TRUE;
      }
   }

   /* get the right hand side of the cut and convert it to a bound */
   rhs = SCIProwGetRhs(cut);
   if( !SCIPsetIsInfinity(set, rhs) )
   {
      rhs -= SCIProwGetConstant(cut);
      if( vals[0] > 0.0 )
      {
         /* coefficient is positive -> rhs corresponds to upper bound */
         SCIP_Real newub;

         newub = rhs/vals[0];
         SCIPvarAdjustUb(var, set, &newub);

         /* bound changes that improve the bound sufficiently are applicable */
         if( SCIPsetIsFeasLT(set, newub, oldlb) || SCIPsetIsLT(set, MAX(newub, oldlb), oldub) )
            return TRUE;
      }
      else
      {
         /* coefficient is negative -> rhs corresponds to lower bound */
         SCIP_Real newlb;

         newlb = rhs/vals[0];
         SCIPvarAdjustLb(var, set, &newlb);

         /* bound changes that improve the bound sufficiently are applicable */
         if( SCIPsetIsFeasGT(set, newlb, oldub) || SCIPsetIsGT(set, MIN(newlb, oldub), oldlb) )
            return TRUE;
      }
   }

   return FALSE;
}

/** removes a non-forced cut from the separation storage */
static
SCIP_RETCODE sepastoreDelCut(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   int                   pos                 /**< position of cut to delete */
   )
{
   assert(sepastore != NULL);
   assert(sepastore->cuts != NULL);
   assert(sepastore->nforcedcuts <= pos && pos < sepastore->ncuts);

   /* check, if the row deletions from separation storage events are tracked if so, issue ROWDELETEDSEPA event */
   if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDSEPA) != 0 )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowDeletedSepa(&event, blkmem, sepastore->cuts[pos]) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
   }

   /* release the row */
   if( !sepastore->initiallp )
   {
      sepastore->ncutsadded--;
      if( sepastore->cuts[pos]->fromcutpool )
         sepastore->ncutsaddedviapool--;
      else
         sepastore->ncutsaddeddirect--;

      if( (SCIP_ROWORIGINTYPE) sepastore->cuts[pos]->origintype == SCIP_ROWORIGINTYPE_SEPA )
      {
         SCIP_SEPA* sepa;
         sepa = SCIProwGetOriginSepa(sepastore->cuts[pos]);
         SCIPsepaDecNCutsAdded(sepa, sepastore->cuts[pos]->fromcutpool);
      }
   }
   SCIP_CALL( SCIProwRelease(&sepastore->cuts[pos], blkmem, set, lp) );

   /* move last cut to the empty position */
   sepastore->cuts[pos] = sepastore->cuts[sepastore->ncuts-1];
   sepastore->ncuts--;

   return SCIP_OKAY;
}

/** adds cut to separation storage and captures it */
SCIP_RETCODE SCIPsepastoreAddCut(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_Bool*            infeasible          /**< pointer to store whether the cut is infeasible */
   )
{
   SCIP_Bool redundant;
   int pos;

   assert(sepastore != NULL);
   assert(sepastore->nforcedcuts <= sepastore->ncuts);
   assert(set != NULL);
   assert(cut != NULL);
   assert(!SCIPsetIsInfinity(set, -SCIProwGetLhs(cut)) || !SCIPsetIsInfinity(set, SCIProwGetRhs(cut)));
   assert(eventqueue != NULL);
   assert(eventfilter != NULL);

   /* debug: check cut for feasibility */
   SCIP_CALL( SCIPdebugCheckRow(set, cut) ); /*lint !e506 !e774*/

   /* the cut will be forced to enter the LP if the dual must be collected and the initial LP is being constructed */
   forcecut = forcecut || (set->lp_alwaysgetduals && sepastore->initiallp);

   /* in the root node, every local cut is a global cut, and global cuts are nicer in many ways ... */
   if( root && SCIProwIsLocal(cut) )
   {
      SCIPsetDebugMsg(set, "change local flag of cut <%s> to FALSE due to addition in root node\n", SCIProwGetName(cut));

      SCIP_CALL( SCIProwChgLocal(cut, FALSE) );

      assert(!SCIProwIsLocal(cut));
   }

   /* check cut for redundancy or infeasibility */
   redundant = sepastoreIsCutRedundantOrInfeasible(sepastore, set, stat, cut, infeasible);
   /* Note that we add infeasible rows in any case, since we cannot be sure whether the return values are handled
    * correctly. In this way, the LP becomes infeasible. */

   /* in each separation round, make sure that at least one (even redundant) cut enters the LP to avoid cycling */
   if( !forcecut && sepastore->ncuts > 0 && redundant )
      return SCIP_OKAY;

   /* if only one cut is currently present in sepastore, it could be redundant; in this case, it can now be removed
    * again, because now a non redundant cut enters the sepastore */
   if( sepastore->ncuts == 1 && sepastoreIsCutRedundant(sepastore, set, stat, sepastore->cuts[0]) )
   {
      /* check, if the row deletions from separation storage events are tracked if so, issue ROWDELETEDSEPA event */
      if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDSEPA) != 0 )
      {
         SCIP_EVENT* event;

         SCIP_CALL( SCIPeventCreateRowDeletedSepa(&event, blkmem, sepastore->cuts[0]) );
         SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
      }
      /* update statistics of total number of found cuts */
      if( !sepastore->initiallp )
      {
         sepastore->ncutsadded--;
         if( sepastore->cuts[0]->fromcutpool )
            sepastore->ncutsaddedviapool--;
         else
            sepastore->ncutsaddeddirect--;

         if( (SCIP_ROWORIGINTYPE) sepastore->cuts[0]->origintype == SCIP_ROWORIGINTYPE_SEPA )
         {
            SCIP_SEPA* sepa;
            sepa = SCIProwGetOriginSepa(sepastore->cuts[0]);
            SCIPsepaDecNCutsAdded(sepa, sepastore->cuts[0]->fromcutpool);
         }
      }

      SCIP_CALL( SCIProwRelease(&sepastore->cuts[0], blkmem, set, lp) );
      sepastore->ncuts = 0;
      sepastore->nforcedcuts = 0;
   }

   /* a cut is forced to enter the LP if
    *  - we construct the initial LP, or
    *  - it has infinite score factor, or
    *  - it is a bound change that can be applied
    * if it is a non-forced cut and no cuts should be added, abort
    */
   forcecut = forcecut || sepastore->initiallp || sepastore->forcecuts || (!SCIProwIsModifiable(cut) && SCIProwGetNNonz(cut) == 1 && sepastoreIsBdchgApplicable(set, cut));
   if( !forcecut && SCIPsetGetSepaMaxcuts(set, root) == 0 )
      return SCIP_OKAY;

   /* get enough memory to store the cut */
   SCIP_CALL( sepastoreEnsureCutsMem(sepastore, set, sepastore->ncuts+1) );
   assert(sepastore->ncuts < sepastore->cutssize);

   SCIPsetDebugMsg(set, "adding cut <%s> to separation storage of size %d (forcecut=%u, len=%d)\n",
      SCIProwGetName(cut), sepastore->ncuts, forcecut, SCIProwGetNNonz(cut));
   /*SCIP_CALL( SCIPprintRow(set->scip, cut, NULL) );*/

   /* capture the cut */
   SCIProwCapture(cut);

   /* add cut to arrays */
   if( forcecut )
   {
      /* make room at the beginning of the array for forced cut */
      pos = sepastore->nforcedcuts;
      sepastore->cuts[sepastore->ncuts] = sepastore->cuts[pos];
      sepastore->nforcedcuts++;
   }
   else
      pos = sepastore->ncuts;

   sepastore->cuts[pos] = cut;

   /* update statistics of total number of found cuts */
   if( !sepastore->initiallp )
   {
      sepastore->ncutsadded++;
      sepastore->ncutsfoundround++;
      if( cut->fromcutpool )
         sepastore->ncutsaddedviapool++;
      else
         sepastore->ncutsaddeddirect++;

      if( (SCIP_ROWORIGINTYPE) cut->origintype == SCIP_ROWORIGINTYPE_SEPA )
      {
         SCIP_SEPA* sepa;

         sepa = SCIProwGetOriginSepa(cut);
         SCIPsepaIncNCutsAdded(sepa, cut->fromcutpool);
      }
   }
   sepastore->ncuts++;

   /* check, if the row addition to separation storage events are tracked if so, issue ROWADDEDSEPA event */
   if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWADDEDSEPA) != 0 )
   {
      SCIP_EVENT* event;

      SCIP_CALL( SCIPeventCreateRowAddedSepa(&event, blkmem, cut) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
   }

   /* If the duals need to be collected, then the infeasible flag is set to FALSE. This ensures that the LP is solved */
   if( set->lp_alwaysgetduals && sepastore->initiallp )
      (*infeasible) = FALSE;

   return SCIP_OKAY;
}

/** applies a lower bound change */
static
SCIP_RETCODE sepastoreApplyLb(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             bound,              /**< new lower bound of variable */
   SCIP_Bool             local,              /**< is it a local bound change? (otherwise global) */
   SCIP_Bool*            applied,            /**< pointer to store whether the domain change was applied */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if an infeasibility has been detected */
   )
{
   assert(sepastore != NULL);
   assert(cutoff != NULL);
   assert(applied != NULL);

   /* adjust bound to the one that would be applied, so the SCIPsetIsGT check below is more reliable */
   SCIPvarAdjustLb(var, set, &bound);

   if( local )
   {
      /* apply the local bound change or detect a cutoff */
      if( SCIPsetIsGT(set, bound, SCIPvarGetLbLocal(var)) )
      {
         SCIPsetDebugMsg(set, " -> applying bound change: <%s>: [%.15g,%.15g] -> [%.15g,%.15g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), bound, SCIPvarGetUbLocal(var));

         /* changing the lower bound to a value >= SCIPinfinity should result in a cutoff,
          * since "infinite" values in solutions are reserved for another meaning
          */
         if( !SCIPsetIsInfinity(set, bound) && SCIPsetIsFeasLE(set, bound, SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(tree), blkmem, set, stat, transprob, origprob, tree,
                  reopt, lp, branchcand, eventqueue, cliquetable, var, bound, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         else
            *cutoff = TRUE;

         *applied = TRUE;
      }
      else
      {
         SCIPsetDebugMsg(set, " -> ignoring bound change: <%s>: [%.15g,%.15g] -> [%.15g,%.15g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), bound, SCIPvarGetUbLocal(var));
      }
   }
   else
   {
      /* apply the global bound change or detect a global cutoff which means we can cutoff the root node */
      if( SCIPsetIsGT(set, bound, SCIPvarGetLbGlobal(var)) )
      {
         SCIPsetDebugMsg(set, " -> applying global bound change: <%s>: [%.15g,%.15g] -> [%.15g,%.15g]\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), bound, SCIPvarGetUbGlobal(var));

         /* changing the lower bound to a value >= SCIPinfinity should result in a cutoff,
          * since "infinite" values in solutions are reserved for another meaning
          */
         if( !SCIPsetIsInfinity(set, bound) && SCIPsetIsFeasLE(set, bound, SCIPvarGetUbGlobal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob, tree, reopt,
                  lp, branchcand, eventqueue, cliquetable, var, bound, SCIP_BOUNDTYPE_LOWER, FALSE) );
         }
         else
         {
            /* we are done with solving since a global bound change is infeasible */
            SCIP_CALL( SCIPnodeCutoff(SCIPtreeGetRootNode(tree), set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
            *cutoff = TRUE;
         }

         *applied = TRUE;
      }
      else
      {
         SCIPsetDebugMsg(set, " -> ignoring global bound change: <%s>: [%.15g,%.15g] -> [%.15g,%.15g]\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), bound, SCIPvarGetUbGlobal(var));
      }
   }

   return SCIP_OKAY;
}

/** applies an upper bound change */
static
SCIP_RETCODE sepastoreApplyUb(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             bound,              /**< new upper bound of variable */
   SCIP_Bool             local,              /**< is it a local bound change? (otherwise global) */
   SCIP_Bool*            applied,            /**< pointer to store whether the domain change was applied */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if an infeasibility has been detected */
   )
{
   assert(sepastore != NULL);
   assert(cutoff != NULL);
   assert(applied != NULL);

   /* adjust bound to the one that would be applied, so the SCIPsetIsGT check below is more reliable */
   SCIPvarAdjustUb(var, set, &bound);

   if( local )
   {
      /* apply the local bound change or detect a cutoff */
      if( SCIPsetIsLT(set, bound, SCIPvarGetUbLocal(var)) )
      {
         SCIPsetDebugMsg(set, " -> applying bound change: <%s>: [%.15g,%.15g] -> [%.15g,%.15g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var), bound);

         /* changing the upper bound to a value <= -SCIPinfinity should result in a cutoff,
          * since "infinite" values in solutions are reserved for another meaning
          */
         if( !SCIPsetIsInfinity(set, -bound) && SCIPsetIsFeasGE(set, bound, SCIPvarGetLbLocal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(tree), blkmem, set, stat, transprob, origprob, tree,
                  reopt, lp, branchcand, eventqueue, cliquetable, var, bound, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else
            *cutoff = TRUE;

         *applied = TRUE;
      }
      else
      {
         SCIPsetDebugMsg(set, " -> ignoring bound change: <%s>: [%.15g,%.15g] -> [%.15g,%.15g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLbLocal(var), bound);
      }
   }
   else
   {
      /* apply the global bound change or detect a global cutoff which means we can cutoff the root node */
      if( SCIPsetIsLT(set, bound, SCIPvarGetUbGlobal(var)) )
      {
         SCIPsetDebugMsg(set, " -> applying global bound change: <%s>: [%.15g,%.15g] -> [%.15g,%.15g]\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), SCIPvarGetLbGlobal(var), bound);

         /* changing the upper bound to a value <= -SCIPinfinity should result in a cutoff,
          * since "infinite" values in solutions are reserved for another meaning
          */
         if( !SCIPsetIsInfinity(set, -bound) && SCIPsetIsFeasGE(set, bound, SCIPvarGetLbGlobal(var)) )
         {
            SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetRootNode(tree), blkmem, set, stat, transprob, origprob, tree, reopt,
                  lp, branchcand, eventqueue, cliquetable, var, bound, SCIP_BOUNDTYPE_UPPER, FALSE) );
         }
         else
         {
            /* we are done with solving since a global bound change is infeasible */
            SCIP_CALL( SCIPnodeCutoff(SCIPtreeGetRootNode(tree), set, stat, tree, transprob, origprob, reopt, lp, blkmem) );
            *cutoff = TRUE;
         }

         *applied = TRUE;
      }
      else
      {
         SCIPsetDebugMsg(set, " -> ignoring global bound change: <%s>: [%.15g,%.15g] -> [%.15g,%.15g]\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), SCIPvarGetLbGlobal(var), bound);
      }
   }

   return SCIP_OKAY;
}

/** applies a cut that is a bound change directly as bound change instead of adding it as row to the LP */
static
SCIP_RETCODE sepastoreApplyBdchg(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_ROW*             cut,                /**< cut with a single variable */
   SCIP_Bool*            applied,            /**< pointer to store whether the domain change was applied */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   SCIP_COL** cols;
   SCIP_Real* vals;
   SCIP_VAR* var;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool local;

   assert(sepastore != NULL);
   assert(!SCIProwIsModifiable(cut));
   assert(SCIProwGetNNonz(cut) == 1);
   assert(cutoff != NULL);
   assert(applied != NULL);

   *applied = FALSE;
   *cutoff = FALSE;

   /* get the single variable and its coefficient of the cut */
   cols = SCIProwGetCols(cut);
   assert(cols != NULL);

   var = SCIPcolGetVar(cols[0]);
   vals = SCIProwGetVals(cut);
   assert(vals != NULL);
   assert(!SCIPsetIsZero(set, vals[0]));

   /* if the coefficient is nearly zero, we better ignore this cut for numerical reasons */
   if( SCIPsetIsFeasZero(set, vals[0]) )
      return SCIP_OKAY;

   local = SCIProwIsLocal(cut);

   /* get the left hand side of the cut and convert it to a bound */
   lhs = SCIProwGetLhs(cut);
   if( !SCIPsetIsInfinity(set, -lhs) )
   {
      lhs -= SCIProwGetConstant(cut);
      if( vals[0] > 0.0 )
      {
         /* coefficient is positive -> lhs corresponds to lower bound */
         SCIP_CALL( sepastoreApplyLb(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               cliquetable, var, lhs/vals[0], local, applied, cutoff) );
      }
      else
      {
         /* coefficient is negative -> lhs corresponds to upper bound */
         SCIP_CALL( sepastoreApplyUb(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               cliquetable, var, lhs/vals[0], local, applied, cutoff) );
      }
   }

   /* get the right hand side of the cut and convert it to a bound */
   rhs = SCIProwGetRhs(cut);
   if( !SCIPsetIsInfinity(set, rhs) )
   {
      rhs -= SCIProwGetConstant(cut);
      if( vals[0] > 0.0 )
      {
         /* coefficient is positive -> rhs corresponds to upper bound */
         SCIP_CALL( sepastoreApplyUb(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               cliquetable, var, rhs/vals[0], local, applied, cutoff) );
      }
      else
      {
         /* coefficient is negative -> rhs corresponds to lower bound */
         SCIP_CALL( sepastoreApplyLb(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand, eventqueue,
               cliquetable, var, rhs/vals[0], local, applied, cutoff) );
      }
   }

   /* count the bound change as applied cut */
   if( *applied && !sepastore->initiallp )
      sepastore->ncutsapplied++;

   return SCIP_OKAY;
}

/** applies the given cut to the LP and updates the orthogonalities and scores of remaining cuts */
static
SCIP_RETCODE sepastoreApplyCut(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_ROW*             cut,                /**< cut to apply to the LP */
   int                   depth,              /**< depth of current node */
   int*                  ncutsapplied        /**< pointer to count the number of applied cuts */
   )
{
   assert(sepastore != NULL);
   assert(ncutsapplied != NULL);

   /* a row could have been added twice to the separation store; add it only once! */
   if( !SCIProwIsInLP(cut) )
   {
      /* add cut to the LP and capture it */
      SCIP_CALL( SCIPlpAddRow(lp, blkmem, set, eventqueue, eventfilter, cut, depth) );

      /* update statistics -> only if we are not in the initial lp (cuts are only counted if added during run) */
      if( !sepastore->initiallp )
      {
         sepastore->ncutsapplied++;

         /* increase count of applied cuts for origins of row */
         /* TODO: adjust conshdlr statistics to mirror cut statistics */
         switch ( (SCIP_ROWORIGINTYPE) cut->origintype )
         {
         case SCIP_ROWORIGINTYPE_CONSHDLR:
            assert( cut->origin != NULL );
            SCIPconshdlrIncNAppliedCuts((SCIP_CONSHDLR*) cut->origin);
            break;
         case SCIP_ROWORIGINTYPE_CONS:
            assert( cut->origin != NULL );
            SCIPconshdlrIncNAppliedCuts(SCIPconsGetHdlr((SCIP_CONS*)cut->origin));
            break;
         case SCIP_ROWORIGINTYPE_SEPA:
            assert( cut->origin != NULL );
            SCIPsepaIncNCutsApplied((SCIP_SEPA*)cut->origin, cut->fromcutpool);
            break;
         case SCIP_ROWORIGINTYPE_UNSPEC:
         case SCIP_ROWORIGINTYPE_REOPT:
            /* do nothing - cannot update statistics */
            break;
         default:
            SCIPerrorMessage("unknown type of row origin.\n");
            return SCIP_INVALIDDATA;
         }
      }

      (*ncutsapplied)++;
   }

   return SCIP_OKAY;
}

/** adds cuts to the LP and clears separation storage */
SCIP_RETCODE SCIPsepastoreApplyCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_EFFICIACYCHOICE  efficiacychoice,    /**< type of solution to base efficiacy computation on */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   SCIP_NODE* node;
   int maxsepacuts;
   int ncutsapplied;
   int nselectedcuts;
   int depth;
   int i;

   assert(sepastore != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(lp != NULL);
   assert(cutoff != NULL);

   SCIP_UNUSED(efficiacychoice);

   SCIPsetDebugMsg(set, "selecting from %d cuts\n", sepastore->ncuts);

   /* get maximal number of cuts to add to the LP */
   maxsepacuts = SCIPsetGetSepaMaxcuts(set, root);

   /* if all cuts are forced cuts, no selection is required */
   if( sepastore->nforcedcuts >= MIN(sepastore->ncuts, maxsepacuts) )
   {
      nselectedcuts = sepastore->nforcedcuts;
   }
   else
   {
      /* call cut selection algorithms */
      nselectedcuts = 0;
      SCIP_CALL( SCIPcutselsSelect(set, sepastore->cuts, sepastore->ncuts, sepastore->nforcedcuts, root, sepastore->initiallp, maxsepacuts,
            &nselectedcuts) );
      assert(nselectedcuts + sepastore->nforcedcuts <= maxsepacuts);

      /* note that cut selector statistics are updated here also when in probing mode; this may lead to an offset with
       * separator/constraint handler statistics */
      /* @todo do not update cutselector statistics if SCIPtreeProbing(scip->tree) */
      nselectedcuts += sepastore->nforcedcuts;
   }

   /*
    * apply all selected cuts
    */
   ncutsapplied = 0;
   *cutoff = FALSE;

   node = SCIPtreeGetCurrentNode(tree);
   assert(node != NULL);

   /* get depth of current node */
   depth = SCIPnodeGetDepth(node);

   for( i = 0; i < nselectedcuts && !(*cutoff); i++ )
   {
      SCIP_ROW* cut;

      cut = sepastore->cuts[i];

      if( i < sepastore->nforcedcuts || SCIPsetIsFeasPositive(set, SCIProwGetLPEfficacy(cut, set, stat, lp)) )
      {
         SCIP_Bool applied = FALSE;

         /* if the cut is a bound change (i.e. a row with only one variable), add it as bound change instead of LP row */
         if( !SCIProwIsModifiable(cut) && SCIProwGetNNonz(cut) == 1 )
         {
            SCIPsetDebugMsg(set, " -> applying forced cut <%s> as boundchange\n", SCIProwGetName(cut));
            SCIP_CALL( sepastoreApplyBdchg(sepastore, blkmem, set, stat, transprob, origprob, tree, reopt, lp, branchcand,
                  eventqueue, cliquetable, cut, &applied, cutoff) );
            assert(applied || !sepastoreIsBdchgApplicable(set, cut));
         }

         if( !applied )
         {
            /* add cut to the LP and update orthogonalities */
            SCIPsetDebugMsg(set, " -> applying%s cut <%s>\n", (i < sepastore->nforcedcuts) ? " forced" : "", SCIProwGetName(cut));
            /*SCIPdebug( SCIProwPrint(cut, set->scip->messagehdlr, NULL));*/
            SCIP_CALL( sepastoreApplyCut(sepastore, blkmem, set, eventqueue, eventfilter, lp, cut, depth, &ncutsapplied) );
         }
      }
   }

   /* clear the separation storage and reset statistics for separation round */
   SCIP_CALL( SCIPsepastoreClearCuts(sepastore, blkmem, set, eventqueue, eventfilter, lp) );

   return SCIP_OKAY;
}

/** clears the separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPsepastoreClearCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp                  /**< LP data */
   )
{
   int c;

   assert(sepastore != NULL);

   SCIPsetDebugMsg(set, "clearing %d cuts\n", sepastore->ncuts);

   /* release cuts */
   for( c = 0; c < sepastore->ncuts; ++c )
   {
      /* check, if the row deletions from separation storage events are tracked if so, issue ROWDELETEDSEPA event */
      if( eventfilter->len > 0 && (eventfilter->eventmask & SCIP_EVENTTYPE_ROWDELETEDSEPA) != 0 )
      {
         SCIP_EVENT* event;

         SCIP_CALL( SCIPeventCreateRowDeletedSepa(&event, blkmem, sepastore->cuts[c]) );
         SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
      }

      SCIP_CALL( SCIProwRelease(&sepastore->cuts[c], blkmem, set, lp) );
   }

   /* reset counters */
   sepastore->ncuts = 0;
   sepastore->nforcedcuts = 0;
   sepastore->ncutsfoundround = 0;

   /* if we have just finished the initial LP construction, free the (potentially large) cuts array */
   if( sepastore->initiallp )
   {
      BMSfreeMemoryArrayNull(&sepastore->cuts);
      sepastore->cutssize = 0;
   }

   return SCIP_OKAY;
}

/** removes cuts that are inefficacious w.r.t. the current LP solution from separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPsepastoreRemoveInefficaciousCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_EFFICIACYCHOICE  efficiacychoice     /**< type of solution to base efficiacy computation on */
   )
{
   int cnt = 0;
   int c;

   assert( sepastore != NULL );

   /* check non-forced cuts only */
   c = sepastore->nforcedcuts;
   while( c < sepastore->ncuts )
   {
      SCIP_Real cutefficacy;

      /* calculate cut's efficacy */
      switch ( efficiacychoice )
      {
         case SCIP_EFFICIACYCHOICE_LP:
            cutefficacy = SCIProwGetLPEfficacy(sepastore->cuts[c], set, stat, lp);
            break;
         case SCIP_EFFICIACYCHOICE_RELAX:
            cutefficacy = SCIProwGetRelaxEfficacy(sepastore->cuts[c], set, stat);
            break;
         case SCIP_EFFICIACYCHOICE_NLP:
            cutefficacy = SCIProwGetNLPEfficacy(sepastore->cuts[c], set, stat);
            break;
         default:
            SCIPerrorMessage("Invalid efficiacy choice.\n");
            return SCIP_INVALIDCALL;
      }

      if( !SCIPsetIsEfficacious(set, root, cutefficacy) )
      {
         SCIP_CALL( sepastoreDelCut(sepastore, blkmem, set, eventqueue, eventfilter, lp, c) );
         ++cnt;
      }
      else
         ++c;
   }
   SCIPsetDebugMsg(set, "removed %d non-efficacious cuts\n", cnt);

   return SCIP_OKAY;
}

/** indicates whether a cut is applicable
 *
 *  A cut is applicable if it is modifiable, not a bound change, or a bound change that changes bounds by at least epsilon.
 */
SCIP_Bool SCIPsepastoreIsCutApplicable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             cut                 /**< cut to check */
   )
{
   return SCIProwIsModifiable(cut) || SCIProwGetNNonz(cut) != 1 || sepastoreIsBdchgApplicable(set, cut);
}

/** get cuts in the separation storage */
SCIP_ROW** SCIPsepastoreGetCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->cuts;
}

/** get number of cuts in the separation storage */
int SCIPsepastoreGetNCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncuts;
}

/** gets the total number of cutting planes added to the separation storage;
 *  this is equal to the sum of added cuts directly and via the pool. */
int SCIPsepastoreGetNCutsAdded(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsadded;
}

/** gets the number of cutting planes added to the separation storage from the cut pool */
int SCIPsepastoreGetNCutsAddedViaPool(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsaddedviapool;
}

/** gets the number of cutting planes added to the separation storage directly */
int SCIPsepastoreGetNCutsAddedDirect(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsaddeddirect;
}

/** get number of cuts found so far in current separation round */
int SCIPsepastoreGetNCutsFoundRound(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsfoundround;
}

/** gets the total number of cutting planes applied to the LP */
int SCIPsepastoreGetNCutsApplied(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   )
{
   assert(sepastore != NULL);

   return sepastore->ncutsapplied;
}
