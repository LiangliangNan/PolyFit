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

/**@file   stat.c
 * @brief  methods for problem statistics
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Gerald Gamrath
 * @author Marc Pfetsch
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/prob.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/visual.h"
#include "scip/mem.h"
#include "scip/var.h"
#include "scip/history.h"
#include "scip/concsolver.h"



/** creates problem statistics data */
SCIP_RETCODE SCIPstatCreate(
   SCIP_STAT**           stat,               /**< pointer to problem statistics data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL */
   SCIP_PROB*            origprob,           /**< original problem, or NULL */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert(stat != NULL);
   assert(set != NULL);

   SCIP_ALLOC( BMSallocMemory(stat) );

   SCIP_CALL( SCIPclockCreate(&(*stat)->solvingtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->solvingtimeoverall, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->presolvingtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->presolvingtimeoverall, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->primallptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->duallptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->lexduallptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->barrierlptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->divinglptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->strongbranchtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->conflictlptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->lpsoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->relaxsoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->pseudosoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->sbsoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->nodeactivationtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->nlpsoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->copyclock, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->strongpropclock, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->reoptupdatetime, SCIP_CLOCKTYPE_DEFAULT) );

   /* turn statistic timing on or off, depending on the user parameter */
   SCIPstatEnableOrDisableStatClocks(*stat, set->time_statistictiming);

   SCIP_CALL( SCIPhistoryCreate(&(*stat)->glbhistory, blkmem) );
   SCIP_CALL( SCIPhistoryCreate(&(*stat)->glbhistorycrun, blkmem) );
   SCIP_CALL( SCIPvisualCreate(&(*stat)->visual, messagehdlr) );

   SCIP_CALL( SCIPregressionCreate(&(*stat)->regressioncandsobjval) );

   (*stat)->status = SCIP_STATUS_UNKNOWN;
   (*stat)->marked_nvaridx = 0;
   (*stat)->marked_ncolidx = 0;
   (*stat)->marked_nrowidx = 0;
   (*stat)->userinterrupt = FALSE;
   (*stat)->userrestart = FALSE;
   (*stat)->inrestart = FALSE;
   (*stat)->collectvarhistory = TRUE;
   (*stat)->performpresol = FALSE;
   (*stat)->branchedunbdvar = FALSE;
   (*stat)->disableenforelaxmsg = FALSE;
   (*stat)->subscipdepth = 0;
   (*stat)->detertimecnt = 0.0;
   (*stat)->nreoptruns = 0;

   SCIPstatReset(*stat, set, transprob, origprob);

   return SCIP_OKAY;
}

/** frees problem statistics data */
SCIP_RETCODE SCIPstatFree(
   SCIP_STAT**           stat,               /**< pointer to problem statistics data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(stat != NULL);
   assert(*stat != NULL);

   SCIPclockFree(&(*stat)->solvingtime);
   SCIPclockFree(&(*stat)->solvingtimeoverall);
   SCIPclockFree(&(*stat)->presolvingtime);
   SCIPclockFree(&(*stat)->presolvingtimeoverall);
   SCIPclockFree(&(*stat)->primallptime);
   SCIPclockFree(&(*stat)->duallptime);
   SCIPclockFree(&(*stat)->lexduallptime);
   SCIPclockFree(&(*stat)->barrierlptime);
   SCIPclockFree(&(*stat)->divinglptime);
   SCIPclockFree(&(*stat)->strongbranchtime);
   SCIPclockFree(&(*stat)->conflictlptime);
   SCIPclockFree(&(*stat)->lpsoltime);
   SCIPclockFree(&(*stat)->relaxsoltime);
   SCIPclockFree(&(*stat)->pseudosoltime);
   SCIPclockFree(&(*stat)->sbsoltime);
   SCIPclockFree(&(*stat)->nodeactivationtime);
   SCIPclockFree(&(*stat)->nlpsoltime);
   SCIPclockFree(&(*stat)->copyclock);
   SCIPclockFree(&(*stat)->strongpropclock);
   SCIPclockFree(&(*stat)->reoptupdatetime);

   SCIPhistoryFree(&(*stat)->glbhistory, blkmem);
   SCIPhistoryFree(&(*stat)->glbhistorycrun, blkmem);
   SCIPvisualFree(&(*stat)->visual);

   SCIPregressionFree(&(*stat)->regressioncandsobjval);

   BMSfreeMemory(stat);

   return SCIP_OKAY;
}

/** diables the collection of any statistic for a variable */
void SCIPstatDisableVarHistory(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->collectvarhistory = FALSE;
}

/** enables the collection of statistics for a variable */
void SCIPstatEnableVarHistory(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->collectvarhistory = TRUE;
}

/** marks statistics to be able to reset them when solving process is freed */
void SCIPstatMark(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->marked_nvaridx = stat->nvaridx;
   stat->marked_ncolidx = stat->ncolidx;
   stat->marked_nrowidx = stat->nrowidx;
}

/** reset statistics to the data before solving started */
void SCIPstatReset(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL */
   SCIP_PROB*            origprob            /**< original problem, or NULL */
   )
{
   assert(stat != NULL);
   assert(stat->marked_nvaridx >= 0);
   assert(stat->marked_ncolidx >= 0);
   assert(stat->marked_nrowidx >= 0);

   SCIPclockReset(stat->solvingtime);
   SCIPclockReset(stat->presolvingtime);
   SCIPclockReset(stat->primallptime);
   SCIPclockReset(stat->duallptime);
   SCIPclockReset(stat->lexduallptime);
   SCIPclockReset(stat->barrierlptime);
   SCIPclockReset(stat->divinglptime);
   SCIPclockReset(stat->strongbranchtime);
   SCIPclockReset(stat->conflictlptime);
   SCIPclockReset(stat->lpsoltime);
   SCIPclockReset(stat->relaxsoltime);
   SCIPclockReset(stat->pseudosoltime);
   SCIPclockReset(stat->sbsoltime);
   SCIPclockReset(stat->nodeactivationtime);
   SCIPclockReset(stat->nlpsoltime);
   SCIPclockReset(stat->copyclock);
   SCIPclockReset(stat->strongpropclock);

   SCIPhistoryReset(stat->glbhistory);

   stat->lastsblpsolstats[0] = stat->lastsblpsolstats[1] = SCIP_LPSOLSTAT_NOTSOLVED;

   stat->vsidsweight = 1.0;
   stat->nlpiterations = 0;
   stat->nrootlpiterations = 0;
   stat->nrootfirstlpiterations = 0;
   stat->nprimallpiterations = 0;
   stat->nduallpiterations = 0;
   stat->nlexduallpiterations = 0;
   stat->nbarrierlpiterations = 0;
   stat->nprimalresolvelpiterations = 0;
   stat->ndualresolvelpiterations = 0;
   stat->nlexdualresolvelpiterations = 0;
   stat->nnodelpiterations = 0;
   stat->ninitlpiterations = 0;
   stat->ndivinglpiterations = 0;
   stat->nsbdivinglpiterations = 0;
   stat->nsblpiterations = 0;
   stat->nsbtimesiterlimhit = 0L;
   stat->nrootsblpiterations = 0;
   stat->nconflictlpiterations = 0;
   stat->ntotalnodes = 0;
   stat->ntotalinternalnodes = 0;
   stat->ntotalnodesmerged = 0;
   stat->ncreatednodes = 0;
   stat->nlpsolsfound = 0;
   stat->nrelaxsolsfound = 0;
   stat->npssolsfound = 0;
   stat->nsbsolsfound = 0;
   stat->nlpbestsolsfound = 0;
   stat->nrelaxbestsolsfound = 0;
   stat->npsbestsolsfound = 0;
   stat->nsbbestsolsfound = 0;
   stat->nexternalsolsfound = 0;
   stat->domchgcount = 0;
   stat->nboundchgs = 0;
   stat->nholechgs = 0;
   stat->nprobboundchgs = 0;
   stat->nprobholechgs = 0;
   stat->nsbdowndomchgs = 0;
   stat->nsbupdomchgs = 0;
   stat->nruns = 0;
   stat->nconfrestarts = 0;
   stat->nrootboundchgs = 0;
   stat->nrootintfixings = 0;
   stat->prevrunnvars = 0;
   stat->nvaridx = stat->marked_nvaridx;
   stat->ncolidx = stat->marked_ncolidx;
   stat->nrowidx = stat->marked_nrowidx;
   stat->nnz = 0;
   stat->avgnnz = 0;
   stat->lpcount = 0;
   stat->relaxcount = 0;
   stat->nlps = 0;
   stat->nrootlps = 0;
   stat->nprimallps = 0;
   stat->nprimalzeroitlps = 0;
   stat->nduallps = 0;
   stat->ndualzeroitlps = 0;
   stat->nlexduallps = 0;
   stat->nbarrierlps = 0;
   stat->nbarrierzeroitlps = 0;
   stat->nprimalresolvelps = 0;
   stat->ndualresolvelps = 0;
   stat->nlexdualresolvelps = 0;
   stat->nnodelps = 0;
   stat->nisstoppedcalls = 0;
   stat->ninitlps = 0;
   stat->ndivinglps = 0;
   stat->nsbdivinglps = 0;
   stat->nnumtroublelpmsgs = 0;
   stat->nstrongbranchs = 0;
   stat->nrootstrongbranchs = 0;
   stat->nconflictlps = 0;
   stat->nnlps = 0;
   stat->maxtotaldepth = -1;
   stat->nactiveconss = 0;
   stat->nenabledconss = 0;
   stat->solindex = 0;
   stat->memsavemode = FALSE;
   stat->nnodesbeforefirst = -1;
   stat->ninitconssadded = 0;
   stat->nactiveconssadded = 0;
   stat->externmemestim = 0;
   stat->nincseparounds = 0;
   stat->nrunsbeforefirst = -1;
   stat->firstprimalheur = NULL;
   stat->firstprimaltime = SCIP_DEFAULT_INFINITY;
   stat->firstprimalbound = SCIP_DEFAULT_INFINITY;
   stat->firstsolgap = SCIP_DEFAULT_INFINITY;
   stat->lastsolgap = SCIP_DEFAULT_INFINITY;
   stat->primalzeroittime = 0.0;
   stat->dualzeroittime = 0.0;
   stat->barrierzeroittime = 0.0;
   stat->maxcopytime = SCIP_REAL_MIN;
   stat->mincopytime = SCIP_REAL_MAX;
   stat->firstlptime = 0.0;
   stat->firstlpdualbound = SCIP_UNKNOWN;
   stat->ncopies = 0;
   stat->nclockskipsleft = 0;
   stat->marked_nvaridx = -1;
   stat->marked_ncolidx = -1;
   stat->marked_nrowidx = -1;
   stat->branchedunbdvar = FALSE;
   stat->bestefficacy = 0.0;
   stat->minefficacyfac = 0.5;
   stat->ncutpoolfails = 0;

   stat->ndivesetlpiterations = 0;
   stat->ndivesetcalls = 0;
   stat->ndivesetlps = 0;
   stat->totaldivesetdepth = 0;

   SCIPstatResetImplications(stat);
   SCIPstatResetPresolving(stat, set, transprob, origprob);
   SCIPstatResetPrimalDualIntegral(stat, set, FALSE);
}

/** reset implication counter */
void SCIPstatResetImplications(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->nimplications = 0;
}

/** reset presolving and current run specific statistics */
void SCIPstatResetPresolving(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL if not yet existing */
   SCIP_PROB*            origprob            /**< original problem, or NULL */
   )
{
   assert(stat != NULL);

   stat->npresolrounds = 0;
   stat->npresolroundsfast = 0;
   stat->npresolroundsmed = 0;
   stat->npresolroundsext = 0;
   stat->npresolfixedvars = 0;
   stat->npresolaggrvars = 0;
   stat->npresolchgvartypes = 0;
   stat->npresolchgbds = 0;
   stat->npresoladdholes = 0;
   stat->npresoldelconss = 0;
   stat->npresoladdconss = 0;
   stat->npresolupgdconss = 0;
   stat->npresolchgcoefs = 0;
   stat->npresolchgsides = 0;

   SCIPstatResetCurrentRun(stat, set, transprob, origprob, FALSE);
}

/** reset primal-dual integral */
void SCIPstatResetPrimalDualIntegral(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             partialreset        /**< should time and integral value be kept? (in combination with no statistical
                                              *   reset, integrals are added for each problem to be solved) */
   )
{
   assert(stat != NULL);

   stat->previousgap = 100.0;
   stat->lastprimalbound = SCIP_UNKNOWN;
   stat->lastdualbound = SCIP_UNKNOWN;
   stat->lastlowerbound = -SCIPsetInfinity(set);
   stat->lastupperbound = SCIPsetInfinity(set);

   /* partial resets keep the integral value and previous evaluation time */
   if( !partialreset )
   {
      stat->previntegralevaltime = 0.0;
      stat->primaldualintegral = 0.0;
   }
}

/** update the primal-dual integral statistic. method accepts + and - SCIPsetInfinity() as values for
 *  upper and lower bound, respectively
 */
void SCIPstatUpdatePrimalDualIntegral(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_Real             upperbound,         /**< current upper bound in transformed problem, or infinity */
   SCIP_Real             lowerbound          /**< current lower bound in transformed space, or -infinity */
   )
{
   SCIP_Real currentgap;
   SCIP_Real solvingtime;
   SCIP_Real primalbound;
   SCIP_Real dualbound;

   assert(stat != NULL);
   assert(set != NULL);

   solvingtime = SCIPclockGetTime(stat->solvingtime);
   assert(solvingtime >= stat->previntegralevaltime);

   if( !SCIPsetIsInfinity(set, upperbound) ) /*lint !e777*/
   {
      /* get value in original space for gap calculation */
      primalbound = SCIPprobExternObjval(transprob, origprob, set, upperbound);

      if( SCIPsetIsZero(set, primalbound) )
         primalbound = 0.0;
   }
   else
   {
      /* no new upper bound: use stored values from last update */
      upperbound = stat->lastupperbound;
      primalbound = stat->lastprimalbound;
      assert(SCIPsetIsZero(set, primalbound) == (primalbound == 0.0)); /*lint !e777*/
   }

   if( !SCIPsetIsInfinity(set, -lowerbound) ) /*lint !e777*/
   {
      /* get value in original space for gap calculation */
      dualbound = SCIPprobExternObjval(transprob, origprob, set, lowerbound);

      if( SCIPsetIsZero(set, dualbound) )
         dualbound = 0.0;
   }
   else
   {
      /* no new lower bound: use stored values from last update */
      lowerbound = stat->lastlowerbound;
      dualbound = stat->lastdualbound;
      assert(SCIPsetIsZero(set, dualbound) == (dualbound == 0.0)); /*lint !e777*/
   }

   /* computation of the gap, special cases are handled first */
   if( primalbound == SCIP_UNKNOWN || dualbound == SCIP_UNKNOWN ) /*lint !e777*/
      currentgap = 100.0;
   /* the gap is 0.0 if bounds coincide */
   else if( SCIPsetIsGE(set, lowerbound, upperbound) || SCIPsetIsEQ(set, primalbound, dualbound) )
      currentgap = 0.0;
   /* the gap is 100.0 if bounds have different signs */
   else if( primalbound * dualbound <= 0.0 ) /*lint !e777*/
      currentgap = 100.0;
   else if( !SCIPsetIsInfinity(set, REALABS(primalbound)) && !SCIPsetIsInfinity(set, REALABS(dualbound)) )
   {
      SCIP_Real absprim = REALABS(primalbound);
      SCIP_Real absdual = REALABS(dualbound);

      /* The gap in the definition of the primal-dual integral differs from the default SCIP gap function.
       * Here, the MAX(primalbound, dualbound) is taken for gap quotient in order to ensure a gap <= 100.
       */
      currentgap = 100.0 * REALABS(primalbound - dualbound) / MAX(absprim, absdual);
      assert(SCIPsetIsLE(set, currentgap, 100.0));
   }
   else
      currentgap = 100.0;

   /* if primal and dual bound have opposite signs, the gap always evaluates to 100.0% */
   assert(currentgap == 0.0 || currentgap == 100.0 || SCIPsetIsGE(set, primalbound * dualbound, 0.0));

   /* update the integral based on previous information */
   stat->primaldualintegral += (solvingtime - stat->previntegralevaltime) * stat->previousgap;

   /* update all relevant information for next evaluation */
   stat->previousgap = currentgap;
   stat->previntegralevaltime = solvingtime;
   stat->lastprimalbound = primalbound;
   stat->lastdualbound = dualbound;
   stat->lastlowerbound = lowerbound;
   stat->lastupperbound = upperbound;
}

/** update and return the primal-dual integral statistic */
SCIP_Real SCIPstatGetPrimalDualIntegral(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob            /**< original problem */
   )
{
   assert(stat != NULL);
   assert(set != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);

   /* update the primal dual integral first */
   SCIPstatUpdatePrimalDualIntegral(stat, set, transprob, origprob, SCIPsetInfinity(set), -SCIPsetInfinity(set));

   return stat->primaldualintegral;
}

/** reset current branch and bound run specific statistics */
void SCIPstatResetCurrentRun(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL */
   SCIP_PROB*            origprob,           /**< original problem, or NULL */
   SCIP_Bool             solved              /**< is problem already solved? */
   )
{
   assert(stat != NULL);

   stat->nnodes = 0;
   stat->ninternalnodes = 0;
   stat->ncreatednodesrun = 0;
   stat->nactivatednodes = 0;
   stat->ndeactivatednodes = 0;
   stat->nbacktracks = 0;
   stat->ndelayedcutoffs = 0;
   stat->nreprops = 0;
   stat->nrepropboundchgs = 0;
   stat->nrepropcutoffs = 0;
   stat->lastdivenode = 0;
   stat->lastconflictnode = 0;
   stat->bestsolnode = 0;
   stat->rootlowerbound = SCIP_REAL_MIN;
   stat->lastbranchvalue = SCIP_UNKNOWN;
   stat->rootlpbestestimate = SCIP_INVALID;
   stat->lastbranchvar = NULL;
   stat->lastbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
   stat->nrootboundchgsrun = 0;
   stat->nrootintfixingsrun = 0;
   stat->npricerounds = 0;
   stat->nseparounds = 0;
   stat->maxdepth = -1;
   stat->plungedepth = 0;
   stat->nobjleaves = 0;
   stat->ninfeasleaves = 0;
   stat->nfeasleaves = 0;
   stat->branchedunbdvar = FALSE;
   stat->nnumtroublelpmsgs = 0;

   stat->nearlybacktracks = 0;
   stat->nnodesaboverefbound = 0;

   assert(transprob == NULL || origprob != NULL);
   /* calculate the reference bound in transformed space from the reference value */
   if( transprob != NULL && !SCIPsetIsInfinity(set, SCIPsetGetReferencevalue(set)) )
      stat->referencebound = SCIPprobInternObjval(transprob, origprob, set, SCIPsetGetReferencevalue(set));
   else
      stat->referencebound = SCIPsetInfinity(set);


   if( !solved )
      stat->status = SCIP_STATUS_UNKNOWN;

   SCIPhistoryReset(stat->glbhistorycrun);

   SCIPregressionReset(stat->regressioncandsobjval);

   SCIPstatResetDisplay(stat);
}

/** resets display statistics, such that a new header line is displayed before the next display line */
void SCIPstatResetDisplay(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->lastdispnode = 0;
   stat->ndisplines = 0;
}

/** increases LP count, such that all lazy updates depending on the LP are enforced again */
void SCIPstatEnforceLPUpdates(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->lpcount++;
}

/** depending on the current memory usage, switches mode flag to standard or memory saving mode */
void SCIPstatUpdateMemsaveMode(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_MEM*             mem                 /**< block memory pools */
   )
{
   assert(stat != NULL);
   assert(set != NULL);

   if( SCIPsetIsLT(set, set->mem_savefac, 1.0) )
   {
      SCIP_Longint memused;

      memused = SCIPmemGetUsed(mem);
      if( !stat->memsavemode && memused >= set->mem_savefac * set->limit_memory * 1024.0 * 1024.0 )
      {
         /* switch to memory saving mode */
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "(node %" SCIP_LONGINT_FORMAT ") switching to memory saving mode (mem: %.1fM/%.1fM)\n",
            stat->nnodes, (SCIP_Real)memused/(1024.0*1024.0), set->limit_memory);
         stat->memsavemode = TRUE;
         set->nodesel = NULL;
      }
      else if( stat->memsavemode && memused < 0.5 * set->mem_savefac * set->limit_memory * 1024.0 * 1024.0 )
      {
         /* switch to standard mode */
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "(node %" SCIP_LONGINT_FORMAT ") switching to standard mode (mem: %.1fM/%.1fM)\n",
            stat->nnodes, (SCIP_Real)memused/(1024.0*1024.0), set->limit_memory);
         stat->memsavemode = FALSE;
         set->nodesel = NULL;
      }
   }
   else
      stat->memsavemode = FALSE;
}

/** returns the estimated number of bytes used by extern software, e.g., the LP solver */
SCIP_Longint SCIPstatGetMemExternEstim(
   SCIP_STAT*            stat                /**< dynamic SCIP statistics */
   )
{
   return stat->externmemestim;
}

/** enables or disables all statistic clocks of \p stat concerning LP execution time, strong branching time, etc.
 *
 *  @note: The (pre-)solving time clocks which are relevant for the output during (pre-)solving
 *         are not affected by this method
 *
 *  @see: For completely disabling all timing of SCIP, consider setting the parameter timing/enabled to FALSE
 */
void SCIPstatEnableOrDisableStatClocks(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   SCIP_Bool             enable              /**< should the LP clocks be enabled? */
   )
{
   assert(stat != NULL);

   SCIPclockEnableOrDisable(stat->primallptime, enable);
   SCIPclockEnableOrDisable(stat->duallptime, enable);
   SCIPclockEnableOrDisable(stat->lexduallptime, enable);
   SCIPclockEnableOrDisable(stat->barrierlptime, enable);
   SCIPclockEnableOrDisable(stat->divinglptime, enable);
   SCIPclockEnableOrDisable(stat->strongbranchtime, enable);
   SCIPclockEnableOrDisable(stat->conflictlptime, enable);
   SCIPclockEnableOrDisable(stat->lpsoltime, enable);
   SCIPclockEnableOrDisable(stat->relaxsoltime, enable);
   SCIPclockEnableOrDisable(stat->pseudosoltime, enable);
   SCIPclockEnableOrDisable(stat->sbsoltime, enable);
   SCIPclockEnableOrDisable(stat->nodeactivationtime, enable);
   SCIPclockEnableOrDisable(stat->nlpsoltime, enable);
   SCIPclockEnableOrDisable(stat->copyclock, enable);
   SCIPclockEnableOrDisable(stat->strongpropclock, enable);
}

/** recompute root LP best-estimate from scratch */
void SCIPstatComputeRootLPBestEstimate(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             rootlpobjval,       /**< root LP objective value */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars               /**< number of variables */
   )
{
   int v;
   stat->rootlpbestestimate = rootlpobjval;

   /* compute best-estimate contribution for every variable */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real rootlpsol;
      SCIP_Real varminpseudoscore;

      /* stop at the first continuous variable */
      if( !SCIPvarIsIntegral(vars[v]) )
         break;

      rootlpsol = SCIPvarGetRootSol(vars[v]);
      varminpseudoscore = SCIPvarGetMinPseudocostScore(vars[v], stat, set, rootlpsol);
      assert(varminpseudoscore >= 0);
      stat->rootlpbestestimate += varminpseudoscore;

      SCIPstatDebugMsg(stat, "Root LP Estimate initialization: <%s> + %15.9f\n", SCIPvarGetName(vars[v]), varminpseudoscore);
   }
}

/** update root LP best-estimate with changed variable pseudo-costs */
SCIP_RETCODE SCIPstatUpdateVarRootLPBestEstimate(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable with changed pseudo costs */
   SCIP_Real             oldrootpscostscore  /**< old minimum pseudo cost score of variable */
   )
{
   SCIP_Real rootlpsol;
   SCIP_Real varminpseudoscore;

   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE );

   /* entire root LP best-estimate must be computed from scratch first */
   if( stat->rootlpbestestimate == SCIP_INVALID ) /*lint !e777*/
      return SCIP_OKAY;

   rootlpsol = SCIPvarGetRootSol(var);

   /* LP root estimate only works for variables with fractional LP root solution */
   if( SCIPsetIsFeasIntegral(set, rootlpsol) )
      return SCIP_OKAY;

   /* subtract old pseudo cost contribution and add new contribution afterwards */
   stat->rootlpbestestimate -= oldrootpscostscore;

   varminpseudoscore = SCIPvarGetMinPseudocostScore(var, stat, set, rootlpsol);
   assert(varminpseudoscore >= 0.0);
   stat->rootlpbestestimate += varminpseudoscore;

   SCIPstatDebugMsg(stat, "Root LP estimate update: <%s> - %15.9f + %15.9f\n", SCIPvarGetName(var), oldrootpscostscore, varminpseudoscore);

   return SCIP_OKAY;
}

/** prints a debug message */
void SCIPstatPrintDebugMessage(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert( sourcefile != NULL );
   assert( stat != NULL );

   if ( stat->subscipdepth > 0 )
      printf("%d: [%s:%d] debug: ", stat->subscipdepth, sourcefile, sourceline);
   else
      printf("[%s:%d] debug: ", sourcefile, sourceline);

   va_start(ap, formatstr); /*lint !e838*/
   printf(formatstr, ap);
   va_end(ap);
}

/** prints a debug message without precode */
EXTERN
void SCIPstatDebugMessagePrint(
   SCIP_STAT*            stat,               /**< SCIP statistics */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{  /*lint --e{715}*/
   va_list ap;

   va_start(ap, formatstr); /*lint !e838*/
   printf(formatstr, ap);
   va_end(ap);
}
